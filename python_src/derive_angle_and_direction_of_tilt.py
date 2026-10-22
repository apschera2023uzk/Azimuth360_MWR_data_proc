#! /usr/bin/env python3

##############################################################################
# 1st Necessary Modules
##############################################################################

import xarray as xr
import glob
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import sys
sys.path.append('/home/aki/pyrtlib')
from pyrtlib.climatology import AtmosphericProfiles as atmp
from pyrtlib.tb_spectrum import TbCloudRTE
from pyrtlib.utils import ppmv2gkg, mr2rh
import shutil
import subprocess
import matplotlib.ticker as ticker
import matplotlib
from scipy.optimize import minimize_scalar
import warnings

##############################################################################
# 2nd Params:
##############################################################################

matplotlib.use("Agg")
first_col = np.arange(0, 360, 5)
second_col = (first_col + 180) % 360
azi_pairs = np.column_stack([first_col, second_col])
n_sigmas_rfi = 1.3

##############################################################################
# 3rd Argparse
##############################################################################

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="This script derives the instrument tilt and direction from a NetCDF file of TBs of 360° Azimuth scans at scan frequency."
    )
    parser.add_argument(
        "--infile", "-i",
        type=str,
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_JOYCE_Joyhat_202408_obs.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_sinthern_may26_obs.nc"),
        default=os.path.expanduser("~/PhD_data/scans/MWR_scans_vettweiss_may26_obs.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_airport_may26_obs.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_RAO_Foghat_202105_08_obs.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_aachen_may26_obs.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_juelich_may26_obs.nc"),
        help="Input NetCDF file which is already resampled to scan frequency."
    )
    # ./derive_angle_and_direction_of_tilt.py --infile /home/aki/PhD_data/scans/MWR_scans_RAO_Foghat_202105_08.nc
    parser.add_argument(
        "--rttov", "-rt",
        type=str,
        default=os.path.expanduser("~/RTTOV-gb"),
        help="Runscript for RTTOV-gb."
    )
    return parser.parse_args()

##############################################################################
# 4th Functions
##############################################################################

def clausius_clapeyron_liq(temp_celsius):
    # Sättigungsdampfdruck für eine Temperatur in °C
    # es returned in Pa
    # https://en.wikipedia.org/wiki/Latent_heat - enthalpiewerte
    L = 2.5e6
    esl = 610.78 * np.exp(L / 462 * (1/273.15 - 1/(273.15+temp_celsius)))
    return esl

###############################################################################

def rh2ppmv(RH=70, abs_T=273.15+15, p=101325):
    es = clausius_clapeyron_liq(abs_T-273.15)
    e = es * RH / 100
    ppmv = 1000000*e / p
    return ppmv

###############################################################################

def get_TB_opposite_differences(ds, azi_pairs=azi_pairs, i_elev=0):
    tb_da = ds["tb"].isel(elevation=i_elev).mean(dim="time")
    pair_labels = []
    dtbs = np.full((len(ds["azimuth"].values), 14), np.nan)
    
    for i, azi in enumerate(ds["azimuth"].values):
        if i>35:
            anti_index = int(i- (len(ds["azimuth"])/2))
        else:
            anti_index = int(i+ (len(ds["azimuth"])/2))
            
        # Calculate TB difference for every pair of opposite angles:
        dtbs[i, :] =  tb_da.values[i,:]- tb_da.values[anti_index,:]        
        pair_labels.append(str(ds["azimuth"].values[i])+" - "+str(ds["azimuth"].values[anti_index]))
    
    return ds["azimuth"].values, dtbs, pair_labels

###############################################################################

def write1profile2str(t_array, ppmv_array,length_value,\
        p_array, liquid_array, height_in_km=0., deg_lat=50.,\
        zenith_angle=0., clear_sky_bool=True):
    string = ""
    
    if clear_sky_bool:
        liquid_array = np.array([0.]*len(liquid_array))
    
    for value in p_array:
        string+=f"{value:8.4f}\n"
    for value in t_array:
        string+=f"{value:6.3f}\n"
    for value in ppmv_array:
        string+=f"{value:9.4f}\n"
    for value in liquid_array:
        string+=f"{value:12.6E}\n"
    string+=f"{t_array[-1]:10.4f}{p_array[-1]:10.2f}\n"
    string+=f"{height_in_km:6.3f}{deg_lat:6.1f}\n"
    string+=f"{zenith_angle:6.4f}\n"
        
    return string
    
###############################################################################

def get_rttov_outputs(rttovgb_outfile=\
        "/home/aki/RTTOV-gb/rttov_test/test_example_k.1/output_example_k.dat.gfortran-openmp",\
                   n_azis=0,n_levels=0):

    print("Reading in RTTOV-gb output from: ", rttovgb_outfile)
    
    tbs = np.full((n_azis,14), np.nan)
    
    switch = False
    tb_string = ""
    switch_count = 0

    file = open(rttovgb_outfile, "r")
    
    # Read in TBs: 
    for i, line in enumerate(file.readlines()):
            
        if "Profile      " in line:
            prof_idx=int(line.split(" ")[-1])-1
        
        if switch and switch_count<2:
            switch_count+= 1
            tb_string+= line
        elif "CALCULATED BRIGHTNESS TEMPERATURES (K):" in line:
            switch = True
        elif switch:
            switch = False
            liste = tb_string.split(" ")
            tbs_rs = [float(s.strip("\n")) for s in liste if s.strip() != ""]
            tb_string = ""
            switch_count = 0
            tbs[prof_idx,:] = np.array(tbs_rs)
            
    file.close()
    print("Finished reading RTTOV-gb output from: ", rttovgb_outfile)
    return tbs #, trans, trans_by_lev, jacs_by_lev

###############################################################################

def calc_TBs4_72_cases(args, angs=np.array([30]*len(azi_pairs)), model="RTTOV-gb"):
    # Determine prams:
    z, p, d, t, md = atmp.gl_atm(atm=1) # midlatitude summer!
    gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
    rh = mr2rh(p, t, gkg)[0] / 100      

    if model=="RTTOV-gb":
        z_in = z*1000
        ppmv = np.full((len(rh)), np.nan)
        for i, rh1 in enumerate(rh):
            ppmv[i] = rh2ppmv(RH=rh1*100, abs_T=t[i], p=p[i]*100)

        # Processing for several angles:
        profiles = ""
        for ang in angs:
            profile1 = write1profile2str(t[::-1], ppmv[::-1],len(ppmv), p[::-1],\
                    np.array([0.]*len(ppmv)), height_in_km=z_in[1], deg_lat=50.,\
                zenith_angle=(90-ang), clear_sky_bool=True)
            profiles+=profile1

        # After loops - save results:
        dir_name = os.path.dirname(args.infile)
        outfile = os.path.expanduser("~/prof_plev.dat")
        out = open(outfile, "w")
        out.write(profiles)
        out.close()

        # Modify runscript copy prof_plev.dat and run RTTOV-gb:
        nlevels=len(ppmv)
        shutil.copy(outfile, args.rttov+"/rttov_test/test_example_k.1/")
        with open(args.rttov+"/rttov_test/run_apschera.sh", "r") as f:
            lines = f.readlines()
            # Zeile 30 (Index 29) ersetzen
        lines[28] = f"NPROF="+str(len(angs))+"\n"    
        lines[29] = f"NLEVELS={nlevels}\n"
        with open(args.rttov+"/rttov_test/run_apschera.sh", "w") as f:
            f.writelines(lines) 
        subprocess.run(["bash", args.rttov+"/rttov_test/run_apschera.sh",\
                "ARCH=gfortran-openmp"], cwd=args.rttov+"/rttov_test/") 

        # Read outputs:
        tbs = get_rttov_outputs(rttovgb_outfile=\
            args.rttov+"/rttov_test/test_example_k.1/output_example_k.dat.gfortran-openmp",\
            n_azis=len(angs),n_levels=nlevels)

    return tbs
            
###############################################################################

def calc_TB_diffs_360deg4tilt(azis,pair,args, ang0=30, tilt=0.5):

    # Calculate Tbs for all 72 Azimuths:
    angs = []
    for i, azi in enumerate(azis):
        grad = azi-pair[0]
        angs.append(ang0+tilt*np.cos(np.deg2rad(grad)))
    tbs = calc_TBs4_72_cases(args, angs=np.array(angs), model="RTTOV-gb")
    
    # Calculate Tb differences for all angle pairs:
    dtbs_mod = np.full((len(azis), 14), np.nan)
    for i, azi in enumerate(azis):
        if i>35:
            anti_index = int(i- (len(azis)/2))
        else:
            anti_index = int(i+ (len(azis)/2))    
        dtbs_mod[i, :] =  tbs[i,:]- tbs[anti_index,:]    
        
    return dtbs_mod

###############################################################################


def calc_estimated_tilt(azis, dtbs, pair_labels, args, elevation=30,
                        azi_pairs=azi_pairs, delta=0.1, chan_weights=None,
                        k_clip=3.0, max_drop_frac=0.2, n_iter=5,
                        tilt_max=3.0, refine=True, return_keep=False, verbose=True):
    """
    Estimate instrument tilt from an azimuth scan.

    1) Linearised fit: for every pair the TB response to a small tilt (delta) is computed
       once; the best tilt per pair follows analytically (least-squares projection).
    2) Robust outlier removal: azimuths whose RMS residual exceeds median + k_clip*MAD
       are excluded, refit, repeat until stable (max. max_drop_frac of azimuths dropped).
    3) Optional nonlinear refinement of the tilt with the full model (bounded Brent)
       for the best pair and its neighbours.

    Returns: tilt [deg], angle_pair (element of azi_pairs), modelled dTBs (n_azi, n_chan)
             (+ keep mask per azimuth if return_keep=True)
    """
    dtbs = np.array(dtbs, dtype=float)          # copy: caller's array is not modified
    n_azi, n_chan = dtbs.shape
    n_pairs = len(azi_pairs)
    w = np.ones(n_chan) if chan_weights is None else np.asarray(chan_weights, float)
    azis_arr = np.asarray(azis)

    def label(i):
        if pair_labels is not None and len(pair_labels) == n_pairs:
            return pair_labels[i]
        return azi_pairs[i]

    def model(pair, tilt):
        return np.asarray(calc_TB_diffs_360deg4tilt(azis, pair, args,
                                                    ang0=elevation, tilt=tilt), float)

    # ------------------------------------------------------------------
    # 1) Linear responses: dTB(pair, tau) ~ m0 + tau * R[pair]
    # ------------------------------------------------------------------
    m0 = model(azi_pairs[0], 0.0)    # untilted model (if tilt=0 fails, use e.g. 1e-4)
    R = np.stack([(model(p, delta) - m0) / delta for p in azi_pairs])   # (n_pairs, n_azi, n_chan)

    def linear_fit(keep):
        y = np.where(keep[:, None], (dtbs - m0) * w, np.nan)
        taus = np.zeros(n_pairs)
        cost = np.full(n_pairs, np.inf)
        for i in range(n_pairs):
            x = R[i] * w
            ok = np.isfinite(x) & np.isfinite(y)
            if not ok.any() or np.sum(x[ok] ** 2) == 0:
                continue
            tau = np.clip(np.sum(x[ok] * y[ok]) / np.sum(x[ok] ** 2), 0.0, tilt_max)
            taus[i] = tau
            cost[i] = np.mean((y[ok] - tau * x[ok]) ** 2)
        return taus, cost

    # ------------------------------------------------------------------
    # 2) Iterative robust clipping of outlier azimuths
    # ------------------------------------------------------------------
    keep = np.isfinite(dtbs).any(axis=1)
    n_min = int(np.ceil((1 - max_drop_frac) * keep.sum()))

    for it in range(n_iter):
        taus, cost = linear_fit(keep)
        i_best = int(np.argmin(cost))
        fit = m0 + taus[i_best] * R[i_best]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            r_azi = np.sqrt(np.nanmean(((dtbs - fit) * w) ** 2, axis=1))
        valid = np.isfinite(r_azi)
        med = np.nanmedian(r_azi)
        mad = 1.4826 * np.nanmedian(np.abs(r_azi - med))

        new_keep = valid & (r_azi <= med + k_clip * mad) if mad > 0 else valid.copy()
        if new_keep.sum() < n_min:                       # cap number of dropped azimuths
            new_keep = np.zeros(n_azi, bool)
            new_keep[np.argsort(np.where(valid, r_azi, np.inf))[:n_min]] = True

        if verbose:
            print(f"*** iter {it}: pair={label(i_best)}, tilt(lin)={taus[i_best]:.3f}°, "
                  f"dropped azimuths: {azis_arr[~new_keep & valid]}")
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep

    taus, cost = linear_fit(keep)                        # final linear fit on final mask
    i_best = int(np.argmin(cost))
    tilt, i_pair = taus[i_best], i_best

    # ------------------------------------------------------------------
    # 3) Nonlinear refinement with the full model (best pair + neighbours)
    # ------------------------------------------------------------------
    if refine:
        mask = keep[:, None] & np.isfinite(dtbs)

        def cost_full(pair, tau):
            r = (dtbs - model(pair, tau)) * w
            return np.nanmean(r[mask] ** 2)

        best_cost = np.inf
        for i in range(max(0, i_best - 1), min(n_pairs, i_best + 2)):
            lo, hi = max(0.0, taus[i] - 0.3), min(tilt_max, taus[i] + 0.3)
            res = minimize_scalar(lambda t: cost_full(azi_pairs[i], t),
                                  bounds=(lo, hi), method="bounded",
                                  options={"xatol": 1e-3})
            if verbose:
                print(f"    refine pair={label(i)}: tilt={res.x:.3f}°, cost={res.fun:.4g}")
            if res.fun < best_cost:
                best_cost, tilt, i_pair = res.fun, res.x, i

    angle_pair = azi_pairs[i_pair]
    dtbs_mod = model(angle_pair, tilt)

    if verbose:
        print(f"*** Final: tilt={tilt:.3f}°, pair={label(i_pair)}, "
              f"used {keep.sum()}/{n_azi} azimuths")

    if return_keep:
        return tilt, angle_pair, dtbs_mod, keep
    return tilt, angle_pair, dtbs_mod

##############################################################################

FLAG_VARS = ["comb_flag", "obs_flag", "rfi_flag"]

def clear_dataset(ds0, flag_vars=FLAG_VARS, verbose=True):
    """Return a copy of ds0 where all numeric data variables are set to NaN
    wherever any flag in `flag_vars` is active (flag != 0).
    A flag is applied to a variable only if all flag dims are contained in the variable's dims."""
    ds = ds0.copy()

    # boolean "bad" masks; fillna(0) treats missing flag values as clean
    bad_masks = {}
    for f in flag_vars:
        if f not in ds0:
            print(f"  [clear_dataset] '{f}' not in dataset, skipped")
            continue
        bad_masks[f] = ds0[f].fillna(0) != 0       # robust for int, float and bool flags

    used = set()
    for var in ds.data_vars:
        if var in flag_vars:
            continue
        da = ds[var]
        if not np.issubdtype(da.dtype, np.number):  # skip bool, strings, datetimes
            continue

        mask = None
        for f, bad in bad_masks.items():
            if set(bad.dims) <= set(da.dims):
                mask = bad if mask is None else (mask | bad)
                used.add(f)
        if mask is None:
            continue

        ds[var] = da.where(~mask).transpose(*da.dims)
        if verbose:
            frac = float(mask.broadcast_like(da).mean())
            print(f"  [clear_dataset] {var:20s} {100 * frac:5.1f} % set to NaN")

    for f in set(bad_masks) - used:
        print(f"  [clear_dataset] WARNING: '{f}' {bad_masks[f].dims} matched no variable "
              f"- check dimension names!")
    return ds

##############################################################################

def plot_tilt_comparison(azis, dtbs, dtbs_mod, tilt, angle_pair,\
        save_path=os.path.expanduser("~/PhD_plots/tilts/tilt_plot_"),tag="",\
        elevation=90):

    grad  = azis
    theta = np.deg2rad(grad)
    short_tag = tag.split(".")[0]

    # ── Style ────────────────────────────────────────────────────────────────
    plt.rcParams.update({
        "font.family":      "serif",
        "font.serif":       ["Times New Roman", "DejaVu Serif"],
        "font.size":        13,
        "axes.linewidth":   1.2,
        "axes.spines.top":  False,
        "axes.spines.right":False,
        "xtick.direction":  "in",
        "ytick.direction":  "in",
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "legend.frameon":   True,
        "legend.framealpha":0.9,
        "legend.edgecolor": "0.7",
    })

    fig, ax = plt.subplots(figsize=(16, 9))

    # ── Data ─────────────────────────────────────────────────────────────────
    y_meas = np.nanmean(dtbs,     axis=1)   # (azimuth,)
    y_mod  = np.nanmean(dtbs_mod, axis=1)

    ax.plot(grad, y_meas,
            color="#1f4e79", lw=1.8, marker="o", ms=3.5,
            label="Observed TB difference")
    ax.plot(grad, y_mod,
            color="#c0392b", lw=2.0, ls="--",
            label=f"Modelled tilt signal ({tilt:.2f}°)")

    ax.axhline(0, color="0.4", lw=0.8, ls=":")

    # ── Axes ─────────────────────────────────────────────────────────────────
    ax.set_xlim(grad.min(), grad.max())
    ax.set_xlabel("Azimuth angle [°]", fontsize=14)
    ax.set_ylabel(r"Mean $\Delta T_\mathrm{B}$ [K]", fontsize=14)

    ax.set_xticks(grad[::5])
    ax.set_xticklabels(pair_labels[::5], rotation=45, ha='right', fontsize=10)

    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

    # ── Annotations ──────────────────────────────────────────────────────────
    ax.set_title(
        f"Azimuthal brightness temperature differences and estimated instrument tilt ({short_tag} / elevation: {elevation}°)\n"
        f"Tilt magnitude: {tilt:.2f}°   |   "
        f"Identified tilt angles: {angle_pair[0]:.0f}° (high) to {angle_pair[1]:.0f}° (low) | Tilt axes: {(angle_pair[0]+90)%360}° - {(angle_pair[1]+90)%360}°",
        fontsize=13, pad=12, loc="left"
    )

    ax.legend(loc="upper right", fontsize=12)

    # ── Tilt-axis marker ─────────────────────────────────────────────────────
    for ang in angle_pair:
        ax.axvline(ang, color="#c0392b", lw=0.8, ls=":", alpha=0.6)

    fig.tight_layout()
    fig.savefig(save_path+short_tag+"_"+str(elevation)+".png",dpi=300,\
        bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {save_path}"+short_tag+"_"+str(elevation)+".png")

    ########
    # ══════════════════════════════════════════════════════════════════════════
    # 2nd plot: observed differences corrected by the modelled tilt signal
    # ══════════════════════════════════════════════════════════════════════════
    y_corr = np.nanmean(dtbs - dtbs_mod, axis=1)            # (azimuth,)

    rms_before = np.sqrt(np.nanmean(y_meas ** 2))
    rms_after  = np.sqrt(np.nanmean(y_corr ** 2))

    fig, ax = plt.subplots(figsize=(16, 9))

    ax.plot(grad, y_meas,
            color="#1f4e79", lw=1.2, marker="o", ms=2.5, alpha=0.35,
            label=f"Observed TB difference (RMS {rms_before:.2f} K)")
    ax.plot(grad, y_corr,
            color="#117a65", lw=1.8, marker="o", ms=3.5,
            label=f"Tilt-corrected TB difference (RMS {rms_after:.2f} K)")

    ax.axhline(0, color="0.4", lw=0.8, ls=":")

    # ── Axes ─────────────────────────────────────────────────────────────────
    ax.set_xlim(grad.min(), grad.max())
    ax.set_xlabel("Azimuth angle [°]", fontsize=14)
    ax.set_ylabel(r"Mean $\Delta T_\mathrm{B}$ [K]", fontsize=14)

    ax.set_xticks(grad[::5])
    ax.set_xticklabels(pair_labels[::5], rotation=45, ha='right', fontsize=10)

    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

    # ── Annotations ──────────────────────────────────────────────────────────
    ax.set_title(
        f"Azimuthal brightness temperature differences after tilt correction ({short_tag} / elevation: {elevation}°)\n"
        f"Removed tilt signal: {tilt:.2f}°   |   "
        f"RMS reduced from {rms_before:.2f} K to {rms_after:.2f} K "
        f"({100 * (1 - rms_after / rms_before):.0f} %)",
        fontsize=13, pad=12, loc="left"
    )

    ax.legend(loc="upper right", fontsize=12)

    for ang in angle_pair:
        ax.axvline(ang, color="#c0392b", lw=0.8, ls=":", alpha=0.6)

    fig.tight_layout()
    fig.savefig(save_path+short_tag+"_"+str(elevation)+"_corr.png", dpi=300,\
        bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {save_path}"+short_tag+"_"+str(elevation)+"_corr.png")

##############################################################################
# 5th Main code:
##############################################################################

if __name__=="__main__":    
    args = parse_arguments()

    ds0 = xr.open_dataset(args.infile)

    tag=(args.infile).split("/")[-1]
    short_tag = tag.split(".")[0]
    dtbs_mod_all_elevs = np.full((len(ds0["elevation"]), len(ds0["azimuth"]),14), np.nan)

    # Apply Cloud, RFI and obstacle flags:
    ds = clear_dataset(ds0)
    
    ###
    # Tilt derivation:    
    # i_elev=0 # 30° at Joyhat viti / 15° at Vettweiss vitII
    for i_elev, ele in enumerate(ds["elevation"].values):
        azis,dtbs, pair_labels = get_TB_opposite_differences(ds, i_elev=i_elev)
        tilt, angle_pair, dtbs_mod = calc_estimated_tilt(azis,dtbs, pair_labels,\
            args, elevation=ds["elevation"].values[i_elev])

        # Create a comparison plot: => def function
        plot_tilt_comparison(azis, dtbs, dtbs_mod, tilt, angle_pair,\
            tag=tag, elevation=ds["elevation"].values[i_elev])

        # Datafile: => def function
        with open(os.path.expanduser("~/PhD_plots/tilts/tilt_detection_"+\
                short_tag+"_"+str(ele)+".txt"), "w") as f:
            f.writelines("Tilt: "+str(tilt)+"°\n") 
            f.writelines("Angle_pair: "+str(angle_pair)+"°\n")

        # Add this elevation to correction mask:
        dtbs_mod_all_elevs[i_elev,:,:] = -dtbs_mod
        # Later the correction should use tilt derived from all elevations together and not by themselves...
        #######

    ####
    # ADDING CORRECTION:
    ds["tilt_corr"] = xr.DataArray(
        dtbs_mod_all_elevs,
        dims=["elevation", "azimuth", "N_Channels"],
        coords={
            "elevation":  ds["elevation"],
            "azimuth":    ds["azimuth"],
            "N_Channels":    ds["N_Channels"],
        },
        attrs={"long_name": "tilt correction", "units": "K",
               "comment": "Add this mask from TBs to receive tilt corrected values"}
    )
    '''
    tb_tilt_corr = np.full((len(ds0["time"]),len(ds0["elevation"]),\
        len(ds0["azimuth"]),14), np.nan)
    for i in range(len(ds0["time"])):
        tb_tilt_corr[i,:,:,:] = tb[i,:,:,:]+dtbs_mod_all_elevs
    '''
    tb_tilt_corr = ds["tb"].values + dtbs_mod_all_elevs            # (time, ele, azi, 14) + (ele, azi, 14)
    ds["tb_tilt_corr"] = xr.DataArray(
        tb_tilt_corr,
        dims=ds["tb"].dims, 
        coords=ds["tb"].coords,
        attrs={"long_name": "tilt correction", "units": "K",
               "comment": "Add this mask from TBs to receive tilt corrected values"}
    )
    # Remember: That TBs in l1 are not on this grid - you have to apply the correction to
    # tb (time x freq) and elevation (time) later!    

    site_file = args.infile
    ds.to_netcdf(site_file.replace("_obs.nc", "_tiltcorr.nc"))
    ds.close()
    del ds

    ##############
    # 0th add corrected dataset with correction mask NetCDF...
    # Tilt timeseries!!!
    # 3. get another profile from a radiosonde instead of pure clim! => or mean profile of these TBs...
    # Understand current algorithm to find tilt :D























