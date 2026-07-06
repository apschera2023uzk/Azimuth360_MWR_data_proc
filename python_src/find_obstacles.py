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
import gc

##############################################################################
# 2nd Params:
##############################################################################

matplotlib.use("Agg")
first_col = np.arange(0, 360, 5)
second_col = (first_col + 180) % 360
azi_pairs = np.column_stack([first_col, second_col])
n_sigmas_rfi = 1.3
lowest_pos_sigma = 2.5 # add 2. nearly all wrong detection in tophat vanished; at 2.5 all!
min_std_fac = 1.
mod_std_fac = 1.5
sigsig_fac = 3.
lowest_pos_sigsig = 20.



vitII = [os.path.expanduser("~/PhD_data/scans/MWR_scans_sinthern_may26.nc"),\
    os.path.expanduser("~/PhD_data/scans/MWR_scans_vettweiss_may26.nc"),\
    os.path.expanduser("~/PhD_data/scans/MWR_scans_airport_may26.nc"),\
    os.path.expanduser("~/PhD_data/scans/MWR_scans_aachen_may26.nc"),\
    os.path.expanduser("~/PhD_data/scans/MWR_scans_juelich_may26.nc"),\
    # os.path.expanduser("~/PhD_data/scans/MWR_scans_JOYCE_Tophat_202510_12.nc"),\
    os.path.expanduser("~/PhD_data/scans/MWR_scans_RAO_Foghat_202105_08.nc"),\
    os.path.expanduser("~/PhD_data/scans/MWR_scans_foghat_may26.nc")]
'''
vitII = [os.path.expanduser("~/PhD_data/scans/MWR_scans_airport_may26.nc")]
# vitII = [os.path.expanduser("~/PhD_data/scans/MWR_scans_JOYCE_Tophat_202510_12.nc")]
# vitII = [os.path.expanduser("~/PhD_data/scans/MWR_scans_aachen_may26.nc")]
# vitII = [os.path.expanduser("~/PhD_data/scans/MWR_scans_foghat_may26.nc")]
'''
##############################################################################
# 3rd Argparse
##############################################################################

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="This script finds obstacles at different Azimuths and elevations."
    )
    '''
    parser.add_argument(
        "--infile", "-i",
        type=str,
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_JOYCE_Tophat_202510_12.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_sinthern_may26.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_vettweiss_may26.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_airport_may26.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_RAO_Foghat_202105_08.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_aachen_may26.nc"),
        default=os.path.expanduser("~/PhD_data/scans/MWR_scans_juelich_may26.nc"),
        help="Input NetCDF file which is already resampled to scan frequency."
    )
    '''
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

def clear_dataset(ds_in):
    
    n_time = ds_in.sizes["time"]
    clear_mask = np.ones(n_time, dtype=bool)  # alle True = behalten

    if "liquid_cloud_flag" in ds_in:
        cf = ds_in["liquid_cloud_flag"].values
        clear_mask &= (cf == 0)
    else:
        print("Info: liquid_cloud_flag nicht im Dataset — kein Cloud-Filter.")

    if "rainfall_rate" in ds_in:
        rr = ds_in["rainfall_rate"].values
        fill = 9.96921e+36
        rr_clean = np.where(np.abs(rr) >= fill, 0.0, rr)  # FillValue als 0 behandeln
        clear_mask &= (rr_clean == 0)
    else:
        print("Info: rainfall_rate nicht im Dataset — kein Regen-Filter.")

    n_clear = np.sum(clear_mask)
    print(f"Zeitschritte gesamt: {n_time} | nach Filter: {n_clear} "
          f"({100*n_clear/n_time:.1f}% behalten)")

    ds_out = ds_in.isel(time=clear_mask)
    return ds_out

##############################################################################

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
                   n_azis=1,n_levels=0):

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


def get_TBref4elev(args, elevation=30, model="RTTOV-gb"):
    # Determine prams:
    z, p, d, t, md = atmp.gl_atm(atm=1) # midlatitude summer!
    gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
    rh = mr2rh(p, t, gkg)[0] / 100      

    if model=="RTTOV-gb":
        z_in = z*1000
        ppmv = np.full((len(rh)), np.nan)
        for i, rh1 in enumerate(rh):
            ppmv[i] = rh2ppmv(RH=rh1*100, abs_T=t[i], p=p[i]*100)
        profile1 = write1profile2str(t[::-1], ppmv[::-1],len(ppmv), p[::-1],\
                    np.array([0.]*len(ppmv)), height_in_km=z_in[1], deg_lat=50.,\
                zenith_angle=(90-elevation), clear_sky_bool=True)

        # After loops - save results:
        outfile = os.path.expanduser("~/prof_plev.dat")
        out = open(outfile, "w")
        out.write(profile1)
        out.close()

        # Modify runscript copy prof_plev.dat and run RTTOV-gb:
        nlevels=len(ppmv)
        shutil.copy(outfile, args.rttov+"/rttov_test/test_example_k.1/")
        with open(args.rttov+"/rttov_test/run_apschera.sh", "r") as f:
            lines = f.readlines()
            # Zeile 30 (Index 29) ersetzen
        lines[28] = f"NPROF="+str(1)+"\n"    
        lines[29] = f"NLEVELS={nlevels}\n"
        with open(args.rttov+"/rttov_test/run_apschera.sh", "w") as f:
            f.writelines(lines) 
        subprocess.run(["bash", args.rttov+"/rttov_test/run_apschera.sh",\
                "ARCH=gfortran-openmp"], cwd=args.rttov+"/rttov_test/") 

        # Read outputs:
        tbs = get_rttov_outputs(rttovgb_outfile=\
            args.rttov+"/rttov_test/test_example_k.1/output_example_k.dat.gfortran-openmp",\
            n_levels=nlevels)

    return tbs

##############################################################################

def determine_obstacle_probability(ds, i_elev,
        tbs_mod=np.array([54.83,53.18,46.26,33.72,29.82,25.41,23.36,
                          121.3,166.14,263.94,287.9,291.96,292.31,292.48])):

    # Einmal laden — danach nur numpy:
    tb = ds["tb"].isel(elevation=i_elev).values  # (time, azimuth, N_Channels)

    # Referenz: Minimum über Azimuthe, gemittelt über Zeit:
    ref_min = np.nanmin(np.nanmean(tb, axis=0), axis=0)  # (N_Channels,)

    # Alles vektorisiert — keine Schleife mehr:
    tb_mean = np.nanmean(tb, axis=0)   # (azimuth, N_Channels)
    tb_std  = np.nanstd(tb,  axis=0)   # (azimuth, N_Channels)

    overmin_azi = tb_mean - ref_min[np.newaxis, :]          # (azimuth, N_Channels)
    overmod_azi = tb_mean - tbs_mod[np.newaxis, :]          # (azimuth, N_Channels)
    azi_stds    = tb_std                                     # (azimuth, N_Channels)
    p90_min = np.nanpercentile(overmin_azi, 90)
    p90_mod = np.nanpercentile(overmod_azi, 90)
    p90_std = np.nanpercentile(azi_stds, 90)
    p10_min = np.nanpercentile(overmin_azi, 10)
    p10_mod = np.nanpercentile(overmod_azi, 10)
    p10_std = np.nanpercentile(azi_stds, 10)


    return overmin_azi, overmod_azi, azi_stds, p90_min, p90_mod, p90_std, p10_std

##############################################################################

def plot_obstacle_overview(overmin_ele_azi, overmod_ele_azi, std_ele_azi,
                           obstacle_flag, azimuth, elevations,
                           outpath="/home/aki/PhD_plots/", tag=""):
    """
    4-panel colorplot: x=azimuth, y=elevation, color=mean over K-band (Ch 1-7)
    One single plot covering all elevations and azimuths.
    """
    # Mean over K-band channels (index 0-6):
    overmin_2d  = np.nanmean(overmin_ele_azi[:, :, 0:7],  axis=2)  # (elev, azi)
    overmod_2d  = np.nanmean(overmod_ele_azi[:, :, 0:7],  axis=2)
    std_2d      = np.nanmean(std_ele_azi[:,    :, 0:7],   axis=2)
    flag_2d     = np.nanmax( obstacle_flag[:,  :, 0:7],   axis=2)  # flag if ANY channel flagged

    fig, axes = plt.subplots(2, 2, figsize=(18, 10))
    fig.suptitle(f"Obstacle detection — all elevations — K-band mean (Ch 1–7)  [{tag}]",
                 fontsize=14, fontweight="bold")

    plots = [
        (axes[0, 0], overmin_2d,        "Deviation from min TB mean [K]",   "RdYlGn_r", 0, 100),
        (axes[0, 1], overmod_2d,        "Deviation from RTTOV model [K]",   "bwr", -80, 80),
        (axes[1, 0], std_2d,            "Temporal std of TB [K]",            "viridis",   0,    40),
        (axes[1, 1], flag_2d.astype(float), "Obstacle flag (1=obstacle)",    "Reds",      0,    1),
    ]

    for ax, data, title, cmap, vmin, vmax in plots:
        im = ax.pcolormesh(azimuth, elevations, data,
                           cmap=cmap, vmin=vmin, vmax=vmax,
                           shading="auto")
        plt.colorbar(im, ax=ax)
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("Azimuth [°]")
        ax.set_ylabel("Elevation [°]")
        ax.set_xticks(np.arange(0, 361, 30))
        ax.set_yticks(elevations)

    plt.tight_layout()
    plt.savefig(outpath+f"obstacle_plot_{tag}.png", dpi=200, bbox_inches="tight")
    plt.close()

##############################################################################
# 5th Main code:
##############################################################################

if __name__=="__main__":    
    args = parse_arguments()

    # To process other inputs replace site file and loop with args.infile...

    for site_file in vitII:

        ds = xr.open_dataset(site_file, mask_and_scale=False)
        overmin_ele_azi = np.full((len(ds["elevation"]), len(ds["azimuth"]),14), np.nan)
        overmod_ele_azi = np.full((len(ds["elevation"]), len(ds["azimuth"]),14), np.nan)
        std_ele_azi = np.full((len(ds["elevation"]), len(ds["azimuth"]),14), np.nan)
        obstacle_flag = np.full((len(ds["elevation"]), len(ds["azimuth"]),14), np.nan) 

        # Calculate model TBs by elevation:
        tbs_mod_cache = {}
        for ele in ds["elevation"].values:
            tbs_mod_cache[ele] = get_TBref4elev(args, elevation=ele, model="RTTOV-gb")

        ################
        # LOOP BY ELEVS STARTS:
        for i_elev, ele in enumerate(ds["elevation"].values):
            print("*****")
            print(i_elev, ele)

            tbs_mod = tbs_mod_cache[ele]
            overmin_azi, overmod_azi, azi_stds, p90_min, p90_mod, p90_std, p10_std =\
                determine_obstacle_probability(\
                    ds, i_elev, tbs_mod=tbs_mod)

            # Save results:
            overmin_ele_azi[i_elev, :,:] = overmin_azi
            overmod_ele_azi[i_elev, :,:] = overmod_azi
            std_ele_azi[i_elev, :,:] = azi_stds

            # Calculate thresholds:
            sigma_of_mins = np.nanstd(np.nanmean(overmin_azi[:,1:7], axis=1))
            mean_of_mins = np.nanmean(np.nanmean(overmin_azi[:,1:7], axis=1))
            threshold_min = max(mean_of_mins+sigma_of_mins*min_std_fac, mean_of_mins+lowest_pos_sigma)
            sigma_of_mods = np.nanstd(np.nanmean(overmod_azi[:,1:7], axis=1))
            mean_of_mods = np.nanmean(np.nanmean(overmod_azi[:,1:7], axis=1))
            threshold_mod = mean_of_mods+sigma_of_mods*mod_std_fac
            sigma_of_sigmas = np.nanstd(np.nanmean(azi_stds[:,1:7], axis=1))
            mean_of_sigmas = np.nanmean(np.nanmean(azi_stds[:,1:7], axis=1))
            threshold_sig = min(mean_of_sigmas-sigma_of_sigmas*sigsig_fac,\
                mean_of_mins-lowest_pos_sigsig )
            threshold_sig_high = max(mean_of_sigmas+sigma_of_sigmas*sigsig_fac,\
                mean_of_mins+lowest_pos_sigsig )

            # New thresholds:
            thrs_min = max(p90_min, mean_of_mins+lowest_pos_sigma)
            thrs_mod = max(p90_mod, mean_of_mods+lowest_pos_sigma)
            thrs_std_hgh = max(p90_std, mean_of_sigmas+lowest_pos_sigsig)
            thrs_std_low = min(p10_std, mean_of_sigmas-lowest_pos_sigsig/4)

            #################################
            # 1st Write a flag into an array
            '''
            obstacle_flag[i_elev, :, :] = np.where(
                (azi_stds    < thrs_std_low),
                1, 0
            )
            '''

            obstacle_flag[i_elev, :, :] = np.where(
                (overmin_azi > threshold_min) &
                (overmod_azi > threshold_mod), # |
                # (azi_stds    < thrs_std_low) | 
                # (azi_stds > thrs_std_hgh),
                1, 0
            )
            
            
            # Last time I was busy here...
            # For airport data somehow all three flags work fine by themselves, but not together...
            # Why???
            #################
        # LOOP BY ELEV ENDS
        ###################

        ####
        # 2nd create colorplots of stds; overshots and flag
        plot_obstacle_overview(
            overmin_ele_azi, overmod_ele_azi, std_ele_azi,
            obstacle_flag,
            azimuth=ds["azimuth"].values,
            elevations=ds["elevation"].values,
            tag=os.path.basename(site_file).split(".")[0]
        )       

        ###
        # 3rd Write obstacles flag into input file:
        ds["obs_flag"] = xr.DataArray(
            obstacle_flag,
            dims=["elevation", "azimuth", "N_Channels"],
            coords={
                "elevation":  ds["elevation"],
                "azimuth":    ds["azimuth"],
                "N_Channels": ds["N_Channels"],
            },
            attrs={"long_name": "Obstacle flag", "units": "1",
                   "comment": "1=obstacle detected, 0=clean"}
        )

        ds.to_netcdf(site_file.replace(".nc", "_obs.nc"))
        ds.close()
        del ds
        del overmin_ele_azi, overmod_ele_azi, std_ele_azi, obstacle_flag
        gc.collect()

        #############
        # Optional:
        # Checke aussortierenden Algorithmus für Liquid und cloud....
        # Kommt bei Tophat shcon aus resample script nicht mit...
        # Zunächst prozentuale Wahrscheinlichkeit...



















