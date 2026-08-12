#! /usr/bin/env python3

##############################################################################
# 1st Necessary Modules
##############################################################################

import xarray as xr
import glob
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib
import pandas as pd
from datetime import datetime

##############################################################################
# 2nd Params:
##############################################################################

matplotlib.use("Agg")
max_elev_azi_diff = 0.05 #°
azimuths = np.arange(0.,355.1,5.) # Interpolate between these!
elevations = np.array([90, 70, 50,45,40,35, 30, 25, 20, 19.2, 15, 14.4, 11.4, 10, 8.4,  6.6,\
                               5.4, 5, 4.8, 4.2])
FILL_VALUE = -9999

# Alle Höhen-Spalten, die als "gültige Wolkenhöhe" in Frage kommen:
HEIGHT_COLS = ["z_first", "z_second", "z_third", "z_vis",
              "z_low", "z_mid", "z_high"]

##############################################################################
# 3rd Argparse
##############################################################################

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Add pattern of MWR l1 files with scan TBs. And the program returns a NetCDF of all scans in scan frequency."
    )
    parser.add_argument(
        "--in_pattern", "-i",
        type=str,
        # default=os.path.expanduser("~/PhD_data/tophat_joyce_2025/2025/*/sups_joy_mwr00_l1_tb_p00_*.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/joyhat_raw_jun_jul_aug_sep/MWR_1C01_*.nc"),
        # default=os.path.expanduser("~/PhD_data/FESSTVaL_14GB/foghat/l1/*/*/fval_uzk_mwr00_l1_tb*.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/vitII_site_eval/aachen_may26/*/MWR_1C01_aachen_*.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/vitII_site_eval/sinthern_may26/*/MWR_1C01_sinthern_*.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/vitII_site_eval/vettweiss_may26/*/MWR_1C01_vettweiss_*.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/vitII_site_eval/airport_may26/*/MWR_1C01_airport_*.nc"),
        default=os.path.expanduser("~/PhD_data/scans/vitII_site_eval/juelich_may26/sups_joy_mwr00_l1_tb_p00_*.nc"),   
        # default=os.path.expanduser("~/PhD_data/scans/vitII_site_eval/foghat_may26/MWR_1C01_*.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/vitII_site_eval/mecken_jun26/MWR_1C01_*.nc"),
        help="Pattern of MWR output files with TBs of scans."
    )
    parser.add_argument(
        "--outfile", "-o",
        type=str,
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_JOYCE_Tophat_202510_12.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_JOYCE_Joyhat_202406_09.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_sinthern_may26.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_airport_may26.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_vettweiss_may26.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_aachen_may26.nc"),
        default=os.path.expanduser("~/PhD_data/scans/MWR_scans_juelich_may26.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_foghat_may26.nc"),
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_mechat_jun26.nc"),
        help="NetCDF Output file path."
    )
    return parser.parse_args()

##############################################################################
# 4th Functions
##############################################################################

def extract_lowest_cloud_base(folder, pattern="*.txt"):
    """
    Reads all ceilometer files matching `pattern` in `folder`, extracts
    timestamp + lowest valid cloud base height (ignoring FILL_VALUE=-9999),
    and returns one combined, sorted time series.
    """
    files = sorted(glob.glob(os.path.join(folder, pattern)))
    if not files:
        raise ValueError(f"No files found matching {pattern} in {folder}")

    all_series = []

    for f in files:
        df = pd.read_csv(
            f,
            sep="\t",
            comment="#",
            header=None,
            names=["t", "ta", "tb", "N_first", "N_second", "N_third", "N_vis",
                  "N_low", "N_mid", "N_high",
                  "z_first", "z_second", "z_third", "z_vis",
                  "z_low", "z_mid", "z_high"],
        )

        if df.empty:
            continue

        # Zeitstempel parsen:
        times = pd.to_datetime(df["t"], format="%d.%m.%Y %H:%M:%S")

        # Alle Höhen-Spalten: FillValue -> NaN, dann Minimum über die Zeile:
        heights = df[HEIGHT_COLS].replace(FILL_VALUE, np.nan)
        lowest_cbh = heights.min(axis=1, skipna=True)

        s = pd.Series(lowest_cbh.values, index=times)
        all_series.append(s)

    combined = pd.concat(all_series).sort_index()
    combined.name = "lowest_cloud_base_m"
    return combined

##############################################################################

def interpolate_azimuths(ds, ele_var="ele", tb_var="tb"):    
    
    for i, elev in enumerate(ds["elevation"].values):
        ds[tb_var][:,i,:,:] = ds[tb_var].isel(elevation=i)\
            .interpolate_na(dim="azimuth", method="linear")

    for i, elev in enumerate(ds["elevation"].values):
        ds["irt"][:,i,:] = ds["irt"].isel(elevation=i)\
            .interpolate_na(dim="azimuth", method="linear")
   
    return ds

###############################################################################

def determine_scan_slices(ds_in, ele_var="ele", azi_var="azi",\
                          max_elev_azi_diff=max_elev_azi_diff, max_ele=85):
    # Output: list of lists, because each scan is one list of indices!
    time_indices_list_list = []
    BL_time_indices_list_list = []
    scan_switch = False
    ele_old = 90
    
    for i, timestep in enumerate(ds_in["time"].values):
        ele_curr = ds_in[ele_var].values[i]
        azi_curr = ds_in[azi_var].values[i]

        # Values are clearly invalid:
        if ele_curr<0 or ele_curr>180:
            continue

        # Scan is going on:
        if scan_switch and ele_curr<max_ele:
            scan_list.append(i)
            
        # Scan starts:
        elif ele_curr<max_ele:
            scan_switch = True
            scan_list = []
            scan_list.append(i)

        # Scan ends:
        elif max_ele<=ele_curr and scan_switch:
            scan_switch = False
            if len(scan_list)>10: # distinguish BL scans with 10 elevs from Azis at 360°
                #############
                # print("len(scan_list): ", len(scan_list))
                #############
                time_indices_list_list.append(scan_list)
            else:
                BL_time_indices_list_list.append(scan_list)
                
        else:
            continue
            
        ele_old = ds_in[ele_var].values[i]
        azi_old = ds_in[azi_var].values[i]
    
    return time_indices_list_list, BL_time_indices_list_list

###############################################################################

def running_std_10min(ds,  timeslice,channel_idx=6):
    """
    Berechnet die rollende 10-Minuten-Standardabweichung für einen Kanal
    des tb-DataArrays, und extrahiert davon nur die gewünschten Zeitpunkte
    (timeslice = Array von Indizes).
    """
    # tb-Zeitreihe für den gewünschten Kanal als pandas Series
    # (Zeitindex bleibt identisch zur Originalauflösung!)
    times = ds["time"].values
    tb_channel = ds["tb"].values[:,channel_idx] 

    s = pd.Series(tb_channel, index=pd.DatetimeIndex(times))

    # Rollendes 10-Minuten-Fenster, zentriert:
    std_10min = s.rolling("10min", center=True, min_periods=1).std()
    # std_10min.values[timeslice] within intervall
    
    std10_before = std_10min.values[timeslice[0]-20:timeslice[0]-10]
    std10_after = std_10min.values[timeslice[0]+10:timeslice[0]+20]

    # Nur die gewünschten Zeitpunkte (timeslice-Indizes) zurückgeben:
    return std10_before, std10_after

###############################################################################
##############################################################################
# Registry describing how to interpret known flag-like variables.
# "mode012"  -> priority-based combination: any 1 (cloudy) dominates,
#               else any 2 (undefined) dominates, else 0 (clear)
# "bitwise"  -> bitwise OR across all valid values in the timeslice
#               (works for plain 0/1 flags too, since OR of 0/1 == logical OR)
# "mean"     -> not a flag, but a continuous quantity (e.g. rain rate) that
#               gets averaged over the timeslice instead of combined logically
##############################################################################

FLAG_REGISTRY = {
    "liquid_cloud_flag": {"kind": "mode012", "fill": np.nan},
    "flag":              {"kind": "bitwise", "fill": 0},
    "quality_flag":      {"kind": "bitwise", "fill": np.nan},
    "rain_flag":         {"kind": "bitwise", "fill": np.nan},
}

RATE_REGISTRY = {
    "rainfall_rate": {"fill": np.nan},
}


def _extract_fill(ds_old, var, fallback):
    """Prefer the variable's own _FillValue attribute, else use fallback."""
    if var in ds_old and "_FillValue" in ds_old[var].attrs:
        return ds_old[var].attrs["_FillValue"]
    return fallback


def _combine_flag_for_timeslice(vals, kind, fill):
    """Combine all valid flag values within one timeslice into a single value,
    following the rule: a single flagged sample sets the whole scan's flag."""
    valid = vals[~np.isnan(vals) & (vals != fill)] if np.issubdtype(vals.dtype, np.floating) \
        else vals[vals != fill]

    if len(valid) == 0:
        return np.nan

    valid = valid.astype(int)

    if kind == "mode012":
        if 1 in valid:
            return 1
        elif 2 in valid:
            return 2
        else:
            return 0

    elif kind == "bitwise":
        combined = 0
        for v in valid:
            combined |= v
        return combined

    else:
        raise ValueError(f"Unknown flag kind: {kind}")


def _detect_available_flags(ds_old):
    """
    Scan ds_old for known flag/rate variables (from the registries), plus any
    additional '*_flag' variables not explicitly registered (auto-detected,
    treated generically as bitwise flags).
    Returns dict: {var_name: {"kind": ..., "fill": ...}} for flags
           dict: {var_name: {"fill": ...}} for rate-like variables
    """
    found_flags = {}
    for var, meta in FLAG_REGISTRY.items():
        if var in ds_old.data_vars and "time" in ds_old[var].dims:
            found_flags[var] = {
                "kind": meta["kind"],
                "fill": _extract_fill(ds_old, var, meta["fill"]),
            }

    # Auto-detect any further "*_flag" variables not already registered:
    for var in ds_old.data_vars:
        if var not in found_flags and "flag" in var.lower() and "time" in ds_old[var].dims:
            found_flags[var] = {
                "kind": "bitwise",
                "fill": _extract_fill(ds_old, var, np.nan),
            }

    found_rates = {}
    for var, meta in RATE_REGISTRY.items():
        if var in ds_old.data_vars and "time" in ds_old[var].dims:
            found_rates[var] = {"fill": _extract_fill(ds_old, var, meta["fill"])}

    return found_flags, found_rates

##############################################################################

def std10_before_after_zenith(ds_old, timeslice, ele_var="ele",
                              tb_var="tb", channel_idx=6,
                              elevations=elevations, window="10min"):
    """
    Computes the std of TB (channel_idx, default=31GHz/ch7) over the 10min
    window directly before and directly after the scan, using ONLY samples
    with elevation angle == 90° (zenith). Returns (std_before, std_after),
    each a single float (NaN if no zenith samples found in that window).
    """
    scan_start = ds_old["time"].values[timeslice[0]]
    scan_end   = ds_old["time"].values[timeslice[-1]]

    win = np.timedelta64(10, "m")
    before_mask_time = (ds_old["time"].values >= scan_start - win) & \
                      (ds_old["time"].values <  scan_start)
    after_mask_time  = (ds_old["time"].values >  scan_end) & \
                      (ds_old["time"].values <= scan_end + win)

    zenith_mask = np.abs(ds_old[ele_var].values - 90.0) < max_elev_azi_diff

    before_mask = before_mask_time & zenith_mask
    after_mask  = after_mask_time  & zenith_mask

    tb_vals = ds_old[tb_var].values[:, channel_idx]

    std_before = np.nanstd(tb_vals[before_mask]) if before_mask.any() else np.nan
    std_after  = np.nanstd(tb_vals[after_mask])  if after_mask.any()  else np.nan

    return std_before, std_after

##############################################################################

def determine_data_in_time_for_scanset(ds_old, time_indices_list_list,
        tb_var="tb", ele_var="ele", azi_var="azi",
        elevations=elevations, azimuths=azimuths, cloudbases=None):

    n_scans = len(time_indices_list_list)

    # ── TB array (unverändert) ──────────────────────────────────────────────
    tbs = np.full((n_scans, len(elevations), len(azimuths), 14), np.nan)
    irt = np.full((n_scans, len(elevations), len(azimuths)), np.nan)

    # ── Auto-detect which flags/rates actually exist in this dataset ─────────
    available_flags, available_rates = _detect_available_flags(ds_old)

    # Output containers — one array per detected variable, NaN-filled if
    # nothing at all was found for that name:
    flag_outputs = {name: np.full(n_scans, np.nan) for name in available_flags}
    rate_outputs = {name: np.full(n_scans, np.nan) for name in available_rates}

    time_array = []
    std10_before_list = []
    std10_after_list  = []
    ir_cloudflag = []
    cloudbase_list = [] 

    # All scans within MWR file:
    for i, timeslice in enumerate(time_indices_list_list):

        # ── Generic flag aggregation ──────────────────────────────────────────
        for var, meta in available_flags.items():
            vals = ds_old[var].values[timeslice]
            flag_outputs[var][i] = _combine_flag_for_timeslice(
                vals, meta["kind"], meta["fill"])

        # ── Generic rate aggregation (mean, not logical combination) ───────────
        for var, meta in available_rates.items():
            vals = ds_old[var].values[timeslice]
            fill = meta["fill"]
            valid = vals[~np.isnan(vals) & (vals != fill)]
            if len(valid) > 0:
                rate_outputs[var][i] = np.nanmean(valid)

        ###############
        # 1st: 31 GHz std, 10min before and after scan (zenith-only samples):
        std10_before, std10_after = std10_before_after_zenith(
            ds_old, timeslice, ele_var=ele_var, tb_var=tb_var, elevations=elevations)
        std10_before_list.append(std10_before)
        std10_after_list.append(std10_after)

        #########
        # 2nd cloud flag based on T_ir values (if they are there):
        ir_threshold_K = 273.15 - 30  # -30°C in Kelvin
        if "tb_irp" in ds_old:
            ir_var = "tb_irp"
        elif "irt" in ds_old:
            ir_var = "irt"
        else: 
            print("No IR TB var found! Processing expects IR TBs and will crash!")

        ####
        # Always:
        # Average of time over timeslice:
        times = ds_old["time"].values[timeslice]
        times_ns = times.astype("int64")
        mean_ns  = np.nanmean(times_ns)
        mean_time = mean_ns.astype("datetime64[ns]")
        time_array.append(mean_time)

        ###############
        # 3rd: Ceilometer cloud base — mean over the timeslice's time window:
        if cloudbases is not None:
            t_start = times.min()
            t_end   = times.max()
            mask = (cloudbases.index >= t_start) & (cloudbases.index <= t_end)
            vals_in_window = cloudbases.values[mask]
            if len(vals_in_window) > 0 and not np.all(np.isnan(vals_in_window)):
                cloudbase_list.append(np.nanmean(vals_in_window))
            else:
                cloudbase_list.append(np.nan)
        else:
            cloudbase_list.append(np.nan)

        # Check one scan of MWR in:
        for j in timeslice:
            k = np.nanargmin(np.abs(elevations - ds_old[ele_var].values[j]))
            m = np.nanargmin(np.abs(azimuths - ds_old[azi_var].values[j]))
            tbs[i, k, m, :] = ds_old[tb_var].values[j, :]
            irt[i, k, m] = ds_old[ir_var].values[j, 0]

    return (time_array, tbs, flag_outputs, rate_outputs,
            np.array(std10_before_list), np.array(std10_after_list),\
            irt, np.array(cloudbase_list))

###############################################################################

def determine_ds_vars4elev_azi_TB(ds_in):
    
    ele_candidates = ["elevation_angle", "ele", "elevation"]
    azi_candidates = ["azimuth_angle", "azi", "azimuth"]
    tb_candidates = ["tb"]
    
    for ele in ele_candidates:
        if ele in ds_in.data_vars:
            ele_var=ele
    
    for azi in azi_candidates:
        if azi in ds_in.data_vars:
            azi_var=azi
            
    for tb in tb_candidates:
        if tb in ds_in.data_vars:
            tb_var=tb
    
    return ele_var, azi_var, tb_var

###############################################################################

def create_scan_ds(time_array, tbs, flag_outputs, rate_outputs,
                   std10_before, std10_after, irt,cloudbases,
                   elevations=elevations, azimuths=azimuths,
                   ele_var="ele", tb_var="tb", azi_var="azi"):
    # tbs -> DataArray, damit wir bequem über dims prüfen können
    tbs_da = xr.DataArray(
        tbs,
        dims=("time", "elevation", "azimuth", "N_Channels"),
        coords={
            "time": time_array,
            "elevation": elevations,
            "azimuth": azimuths,
            "N_Channels": np.arange(14) + 1,
        },
    )
    has_data_time    = ~np.isnan(tbs_da).all(dim=("elevation", "azimuth", "N_Channels"))
    has_data_elev    = ~np.isnan(tbs_da).all(dim=("time", "azimuth", "N_Channels"))
    tbs_clean = tbs_da.sel(
        time=tbs_da.time[has_data_time],
        elevation=tbs_da.elevation[has_data_elev],
        azimuth=tbs_da.azimuth,
    )
    ds_out = xr.Dataset(
        data_vars={
            "tb": (("time", "elevation", "azimuth", "N_Channels"), tbs_clean.data),
        },
        coords={
            "time": tbs_clean.coords["time"],
            "elevation": tbs_clean.coords["elevation"],
            "azimuth": tbs_clean.coords["azimuth"],
            "N_Channels": tbs_clean.coords["N_Channels"],
        },
    )

    ####
    # Generic flag variables:
    for var_name, flag_vals in flag_outputs.items():
        flag_masked = np.asarray(flag_vals)[has_data_time.values]
        if np.all(np.isnan(flag_masked)):
            continue
        ds_out[var_name] = xr.DataArray(
            flag_masked, dims=("time",),
            coords={"time": tbs_clean.coords["time"]},
            attrs={"units": "1", "long_name": f"{var_name} (auto-aggregated per scan)",
                  "_FillValue": np.nan,
                  "comment": "A single flagged sample within the scan sets the whole scan's flag value."}
        )

    ####
    # Generic rate variables:
    for var_name, rate_vals in rate_outputs.items():
        rate_masked = np.asarray(rate_vals)[has_data_time.values]
        if np.all(np.isnan(rate_masked)):
            continue
        ds_out[var_name] = xr.DataArray(
            rate_masked, dims=("time",),
            coords={"time": tbs_clean.coords["time"]},
            attrs={"units": "unknown", "long_name": f"{var_name} (mean-aggregated per scan)",
                  "_FillValue": np.nan}
        )

    ###
    # std10_before / std10_after — replaces the old std10_cloud_flag:
    std_before_masked = np.asarray(std10_before)[has_data_time.values]
    if not np.all(np.isnan(std_before_masked.astype(float))):
        ds_out["std10_before"] = xr.DataArray(
            std_before_masked, dims=("time",),
            coords={"time": tbs_clean.coords["time"]},
            attrs={
                "units": "K",
                "long_name": "Std of 31 GHz TB (zenith-only) over 10min before scan",
                "comment": "Computed only from samples with elevation angle == 90°.",
            }
        )

    std_after_masked = np.asarray(std10_after)[has_data_time.values]
    if not np.all(np.isnan(std_after_masked.astype(float))):
        ds_out["std10_after"] = xr.DataArray(
            std_after_masked, dims=("time",),
            coords={"time": tbs_clean.coords["time"]},
            attrs={
                "units": "K",
                "long_name": "Std of 31 GHz TB (zenith-only) over 10min after scan",
                "comment": "Computed only from samples with elevation angle == 90°.",
            }
        )

    cloudbases_masked = np.asarray(cloudbases)[has_data_time.values]
    if not np.all(np.isnan(cloudbases_masked.astype(float))):
        ds_out["cloudbases"] = xr.DataArray(
            cloudbases_masked, dims=("time",),
            coords={"time": tbs_clean.coords["time"]},
            attrs={
                "units": "m",
                "long_name": "Cloudbases from Ceilometer",
                "comment": "Nan for clear sky, otherwise cbh",
            }
        )


    ###
    # Add IR brightness temperature (irt), same spatial structure as tb:
    irt_da = xr.DataArray(
        irt,
        dims=("time", "elevation", "azimuth"),
        coords={
            "time": time_array,
            "elevation": elevations,      # volle 20 Elevationen hier
            "azimuth": azimuths,
        },
    )

    # WICHTIG: mit der GLEICHEN Maske filtern wie tbs_clean:
    irt_clean = irt_da.sel(
        time=irt_da.time[has_data_time],
        elevation=irt_da.elevation[has_data_elev],   # ← reduziert auf 4
        azimuth=irt_da.azimuth,
    )

    ds_out["irt"] = xr.DataArray(
        irt_clean.data,                               # jetzt shape (n_time, 4, n_azi)
        dims=("time", "elevation", "azimuth"),
        coords={
            "time": tbs_clean.coords["time"],
            "elevation": tbs_clean.coords["elevation"],  # 4 Elevationen
            "azimuth": tbs_clean.coords["azimuth"],
        },
        attrs={"units": "K", "long_name": "Infrared brightness temperature"}
    )

    ds_out = interpolate_azimuths(ds_out, ele_var=ele_var, tb_var=tb_var)

    return ds_out

###############################################################################

def determine_data_in_time_for_BLset(ds_old, BL_time_indices_list_list,\
                                tb_var="tb", ele_var="ele", azi_var="azi",\
                                elevations=elevations, azimuths=azimuths):
    
    # Create empty TB array with dims: (time,elev,azi,ch)
    tbs = np.full((len(BL_time_indices_list_list),len(elevations),\
                   14), np.nan)
    flags = np.full(len(BL_time_indices_list_list), -2147483647, dtype=int)
    rainfall = np.full(len(BL_time_indices_list_list), np.nan, dtype=float)
    time_array = []
    
    # All scans within MWR file:
    for i, timeslice in enumerate(BL_time_indices_list_list):

        ####
        # Only for Jülich / Vital I:
        if ele_var == "elevation_angle": 

            # Flag: nimm den häufigsten Wert im Timeslice (Mode)
            flag_vals = ds_old["liquid_cloud_flag"].values[timeslice]
            flag_vals_valid = flag_vals[flag_vals != -2147483647]

            if len(flag_vals_valid) > 0:
                if 1 in flag_vals_valid:
                    flags[i] = 1   # cloudy dominiert
                elif 2 in flag_vals_valid:
                    flags[i] = 2   # undefined
                else:
                    flags[i] = 0   # nur dann clear

            # Rainrate:
            rain_vals = ds_old["rainfall_rate"].values[timeslice]
            fill = 9.96921e+36
            rain_vals_valid = rain_vals[rain_vals != fill]
            if len(rain_vals_valid) > 0:
                rainfall[i] = np.nanmean(rain_vals_valid)

        ###
        # Only for FESSTVaL Foghat / RAO:
        else:
            # Bitwise OR: wenn irgendein Zeitschritt ein Bit hat, bleibt es gesetzt
            flag_vals = ds_old["flag"].values[timeslice]
            fill = 0  # _FillValue ist 0s

            ###########################
            # Alternate Nan-Filter:
            # flag_vals_valid = flag_vals[flag_vals != fill].astype(int)
            flag_vals_valid = flag_vals[~np.isnan(flag_vals) & (flag_vals != fill)].astype(int)
            ########################

            if len(flag_vals_valid) > 0:
                combined = 0
                for fv in flag_vals_valid:
                    combined |= int(fv)
                flags[i] = combined
            else:
                flags[i] = 0  # kein gültiger Wert → kein Flag

        ####
        # Always:
        # Average of time over timeslice:
        times = ds_old["time"].values[timeslice]
        times_ns = times.astype("int64")  
        mean_ns  = np.nanmean(times_ns)
        mean_time = mean_ns.astype("datetime64[ns]")
        time_array.append(mean_time)
        
        # Check one scan of MWR in:
        for j in timeslice:
            k = np.nanargmin(np.abs(elevations-ds_old[ele_var].values[j]))
            # m = np.nanargmin(np.abs(azimuths-ds_old[azi_var].values[j]))
            tbs[i, k, :] = ds_old[tb_var].values[j, :]

    # plt.figure()
    # plt.pcolormesh(time_array, elevations, tbs[:,:,0].T)
    # plt.savefig("test_BL.png")

    return time_array, tbs, flags, rainfall

###############################################################################

def create_BL_ds(ds_old, time_array, tbs, flags, rainfall,
                 elevations=elevations,
                 ele_var="ele",
                 tb_var="tb",
                 out_file=None):
    """
    Creates an xarray Dataset for Boundary Layer scans (varying elevation,
    fixed azimuth ~0°/North) and optionally saves it to NetCDF.

    Parameters
    ----------
    ds_old      : xr.Dataset   — original input dataset (for metadata/coords)
    time_array  : list         — mean datetime64 per scan
    tbs         : np.ndarray   — shape (n_scans, n_elevations, 14)
    flags       : np.ndarray   — shape (n_scans,) int cloud/quality flag
    rainfall    : np.ndarray   — shape (n_scans,) float rainfall rate
    elevations  : np.ndarray   — elevation angles used as coordinate
    ele_var     : str          — name of elevation variable in ds_old
    tb_var      : str          — name of TB variable in ds_old
    out_file    : str or None  — path to output NetCDF (None = don't save)

    Returns
    -------
    ds_bl : xr.Dataset
    """

    time_array = np.array(time_array, dtype="datetime64[ns]")
    n_time, n_ele, n_chan = tbs.shape

    # ── Frequency coordinate from original dataset ────────────────────────────
    if "frequency" in ds_old.coords:
        freqs = ds_old["frequency"].values
    elif "n_freq" in ds_old.dims:
        freqs = np.arange(n_chan, dtype=float)
    else:
        freqs = np.arange(n_chan, dtype=float)

    # ── Latitude / Longitude / Altitude scalars from original ─────────────────
    def scalar(var):
        if var in ds_old:
            v = ds_old[var].values.flat[0]
            return float(v) if not np.isnan(float(v)) else np.nan
        return np.nan

    lat = scalar("latitude")
    lon = scalar("longitude")
    alt = scalar("altitude")

    # ── Build Dataset ─────────────────────────────────────────────────────────
    ds_bl = xr.Dataset(
        data_vars={
            # TBs: (time, elevation, frequency)
            "tb": xr.DataArray(
                tbs.astype(np.float32),
                dims=["time", "elevation", "frequency"],
                attrs={
                    "units"    : "K",
                    "long_name": "Microwave brightness temperature",
                    "standard_name": "brightness_temperature",
                }
            ),
            # Cloud / quality flag: (time,)
            "liquid_cloud_flag": xr.DataArray(
                flags.astype(np.int32),
                dims=["time"],
                attrs={
                    "long_name"   : "Liquid cloud flag",
                    "flag_values" : "0 1 2",
                    "flag_meanings": "clear cloudy undefined",
                    "_FillValue"  : -2147483647,
                }
            ),
            # Rainfall rate: (time,)
            "rainfall_rate": xr.DataArray(
                rainfall.astype(np.float32),
                dims=["time"],
                attrs={
                    "units"    : "mm/h",
                    "long_name": "Rainfall rate",
                }
            ),
            # Azimuth: constant 0° (North) for all BL scans
            "azimuth_angle": xr.DataArray(
                np.zeros(n_time, dtype=np.float32),
                dims=["time"],
                attrs={
                    "units"    : "degrees",
                    "long_name": "Azimuth angle (fixed North for BL scans)",
                }
            ),
            # Static scalars:
            "latitude": xr.DataArray(
                np.full(n_time, lat, dtype=np.float32), dims=["time"],
                attrs={"units": "degrees_north", "long_name": "Latitude"}
            ),
            "longitude": xr.DataArray(
                np.full(n_time, lon, dtype=np.float32), dims=["time"],
                attrs={"units": "degrees_east", "long_name": "Longitude"}
            ),
            "altitude": xr.DataArray(
                np.full(n_time, alt, dtype=np.float32), dims=["time"],
                attrs={"units": "m", "long_name": "Altitude above sea level"}
            ),
        },
        coords={
            "time"     : ("time",      time_array),
            "elevation": ("elevation", elevations.astype(np.float32),
                          {"units": "degrees", "long_name": "Elevation angle"}),
            "frequency": ("frequency", freqs.astype(np.float32),
                          {"units": "GHz", "long_name": "Frequency"}),
        },
        attrs={
            "title"      : "MWR Boundary Layer elevation scans",
            "institution": ds_old.attrs.get("institution", ""),
            "source"     : ds_old.attrs.get("source", ""),
            "scan_type"  : "BL_elevation_scan",
            "azimuth_fixed_deg": "0.0 (North)",
            "history"    : f"Created by preproc_resample2scan_freq.py",
            "Conventions": "CF-1.8",
        }
    )

    # ── Save to NetCDF ────────────────────────────────────────────────────────
    if out_file is not None:
        os.makedirs(os.path.dirname(os.path.abspath(out_file)), exist_ok=True)
        ds_bl.to_netcdf(out_file)
        print(f"  Saved BL dataset: {out_file}")

    return ds_bl

###############################################################################

def resample_mwr_ds_on_scan_freq(ds_old, cloudbases=None):
    
    # 0th Determine azi, ele and tb vars:
    ele_var, azi_var, tb_var = determine_ds_vars4elev_azi_TB(ds_old)
    
    # 1st determine timeslices with const elev<90 and variable azi.
    time_indices_list_list, BL_time_indices_list_list =\
        determine_scan_slices(ds_old, ele_var=ele_var, azi_var=azi_var)
    
    #########################
    # 2nd calc mean timestamp and mean measurements for these timeslices
    time_array, tbs, flag_outputs, rate_outputs, std10_before, std10_after,\
            irt, cloudbases =\
            determine_data_in_time_for_scanset(ds_old,\
            time_indices_list_list, tb_var=tb_var, ele_var=ele_var,\
            azi_var=azi_var, cloudbases=cloudbases) 
    time_array_bl, tbs_bl, flags_bl, rainfall_bl =\
            determine_data_in_time_for_BLset(ds_old, BL_time_indices_list_list,\
            tb_var=tb_var, ele_var=ele_var, azi_var=azi_var)

    # 3rd Create new dataset of scans:
    ds_new = create_scan_ds(time_array, tbs, flag_outputs, rate_outputs,
                            std10_before, std10_after, irt,cloudbases,
                            elevations=elevations, azimuths=azimuths,
                            ele_var=ele_var, tb_var=tb_var, azi_var=azi_var)
    ds_bl = create_BL_ds(ds_old,time_array_bl, tbs_bl, flags_bl, rainfall_bl,
                 elevations=elevations,
                 ele_var=ele_var,
                 tb_var=tb_var,
                 out_file=None)

    return ds_new, ds_bl

###############################################################################

def drop_allnan_time(ds, tb_var="tb"):
    """Remove timesteps where the TB variable is entirely NaN."""
    dims_to_check = [d for d in ds[tb_var].dims if d != "time"]
    has_data = ~ds[tb_var].isnull().all(dim=dims_to_check)
    n_before = ds.sizes["time"]
    ds = ds.isel(time=has_data.values)
    n_after = ds.sizes["time"]
    print(f"  Dropped {n_before - n_after} all-NaN timesteps "
          f"({n_before} → {n_after})")
    return ds

##############################################################################
# 5th Main code:
##############################################################################

if __name__ == "__main__":
    args = parse_arguments()
    files = sorted(glob.glob(args.in_pattern))
    n = len(files)
    tmp_files    = []
    tmp_files_bl = []

    cloudbases = extract_lowest_cloud_base(\
            os.path.dirname(args.in_pattern), "ceilo/*_cloudcov.txt")

    print("*********************")
    print(cloudbases)

    for i, file in enumerate(files):
        print(f"Read file {i} of {n}")
        ds = xr.open_dataset(file)
        ds_resamp, ds_bl = resample_mwr_ds_on_scan_freq(ds, cloudbases=cloudbases)

        # ── Azimuth scan dataset ──────────────────────────────────────────────
        if ds_resamp.sizes["time"] == 0:
            print(f"  → no valid azimuth scans in file {i} — skipping")
        else:
            ds_resamp = ds_resamp.assign_coords(
                time=ds_resamp["time"].astype("datetime64[ns]"))
            tmp_path = args.outfile + f".tmp_{i:04d}.nc"
            ds_resamp.to_netcdf(tmp_path)
            tmp_files.append(tmp_path)
            ds_resamp.close()

        # ── BL scan dataset ───────────────────────────────────────────────────
        if ds_bl is None or ds_bl.sizes["time"] == 0:
            print(f"  → no valid BL scans in file {i} — skipping")
        else:
            ds_bl = ds_bl.assign_coords(
                time=ds_bl["time"].astype("datetime64[ns]"))
            tmp_path_bl = args.outfile + f".tmp_bl_{i:04d}.nc"
            ds_bl.to_netcdf(tmp_path_bl)
            tmp_files_bl.append(tmp_path_bl)
            ds_bl.close()

        ds.close()

    # ── Concatenate azimuth scans ─────────────────────────────────────────────
    if tmp_files:
        print("Concatenating azimuth scans...")
        ds_final = xr.open_mfdataset(tmp_files, combine="nested", concat_dim="time")
        ds_final = drop_allnan_time(ds_final, tb_var="tb")   # ← neu
        ds_final.to_netcdf(args.outfile)
        ds_final.close()
        for tmp in tmp_files:
            os.remove(tmp)
        print("Done:", args.outfile)
    else:
        print("No azimuth scan data to write.")

    # ── Concatenate BL scans ──────────────────────────────────────────────────
    if tmp_files_bl:
        outfile_bl = args.outfile.replace(".nc", "_BL.nc")
        print("Concatenating BL scans...")
        ds_final_bl = xr.open_mfdataset(tmp_files_bl, combine="nested", concat_dim="time")
        ds_final_bl = drop_allnan_time(ds_final_bl, tb_var="tb")   # ← neu
        ds_final_bl.to_netcdf(outfile_bl)
        ds_final_bl.close()
        for tmp in tmp_files_bl:
            os.remove(tmp)
        print("Done:", outfile_bl)
    else:
        print("No BL scan data to write.")

###################################################

