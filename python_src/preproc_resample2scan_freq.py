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

##############################################################################
# 2nd Params:
##############################################################################

matplotlib.use("Agg")
max_elev_azi_diff = 0.05 #°
azimuths = np.arange(0.,355.1,5.) # Interpolate between these!
elevations = np.array([90, 70, 50,45,40,35, 30, 25, 20, 19.2, 15, 14.4, 11.4, 10, 8.4,  6.6,\
                               5.4, 5, 4.8, 4.2])

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
        # default=os.path.expanduser("~/PhD_data/scans/vitII_site_eval/juelich_may26/sups_joy_mwr00_l1_tb_p00_*.nc"),   
        default=os.path.expanduser("~/PhD_data/scans/vitII_site_eval/foghat_may26/MWR_1C01_*.nc"),
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
        # default=os.path.expanduser("~/PhD_data/scans/MWR_scans_juelich_may26.nc"),
        default=os.path.expanduser("~/PhD_data/scans/MWR_scans_foghat_may26.nc"),
        help="NetCDF Output file path."
    )
    return parser.parse_args()

##############################################################################
# 4th Functions
##############################################################################

def interpolate_azimuths(ds, ele_var="ele", tb_var="tb"):    
    
    for i, elev in enumerate(ds["elevation"].values):
        ds[tb_var][:,i,:,:] = ds[tb_var].isel(elevation=i)\
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

def determine_data_in_time_for_scanset(ds_old, time_indices_list_list,\
                                tb_var="tb", ele_var="ele", azi_var="azi",\
                                elevations=elevations, azimuths=azimuths):
    
    # Create empty TB array with dims: (time,elev,azi,ch)
    tbs = np.full((len(time_indices_list_list),len(elevations),\
                   len(azimuths), 14), np.nan)
    flags = np.full(len(time_indices_list_list), -2147483647, dtype=int)
    rainfall = np.full(len(time_indices_list_list), np.nan, dtype=float)
    time_array = []
    
    # All scans within MWR file:
    for i, timeslice in enumerate(time_indices_list_list):

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
            m = np.nanargmin(np.abs(azimuths-ds_old[azi_var].values[j]))
            tbs[i, k, m,:] = ds_old[tb_var].values[j, :]

    return time_array, tbs, flags, rainfall

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

def create_scan_ds(time_array, tbs, flags, rainfall, elevations=elevations,\
                    azimuths=azimuths,\
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

    # Maske: wo ist irgendwo ein echter Wert?
    # über N_Channels mitteln (oder any), dann über time prüfen
    has_data_time    = ~np.isnan(tbs_da).all(dim=("elevation", "azimuth", "N_Channels"))
    has_data_elev    = ~np.isnan(tbs_da).all(dim=("time", "azimuth", "N_Channels"))
    # has_data_azimuth = ~np.isnan(tbs_da).all(dim=("time", "elevation", "N_Channels"))

    # nur die Slices behalten, die irgendwo Daten haben
    tbs_clean = tbs_da.sel(
        time=tbs_da.time[has_data_time],
        elevation=tbs_da.elevation[has_data_elev],
        azimuth=tbs_da.azimuth,
    )

    # Neues Dataset aus dem bereinigten DataArray
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
    # Only for Jülich / Vital I:
    if ele_var == "elevation_angle": 
        # Add liquid cloud flag:
        flag_da = xr.DataArray(
            flags[has_data_time.values],   # gleiche Zeitmaske anwenden
            dims=("time",),
            coords={"time": tbs_clean.coords["time"]},
            attrs={
                "units": "1",
                "long_name": "Liquid cloud flag",
                "_FillValue": -2147483647,
                "comment": "Flag meaning: no liquid cloud (0), liquid cloud present (1), undefined (2)"
            }
        )
        ds_out["liquid_cloud_flag"] = flag_da

        # Add rainrate:
        rain_da = xr.DataArray(
            rainfall[has_data_time.values],
            dims=("time",),
            coords={"time": tbs_clean.coords["time"]},
            attrs={
                "units": "m s-1",
                "long_name": "Rainfall rate",
                "standard_name": "rainfall_rate",
                "_FillValue": 9.96921e+36,
            }
        )
        ds_out["rainfall_rate"] = rain_da

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

def resample_mwr_ds_on_scan_freq(ds_old):
    
    # 0th Determine azi, ele and tb vars:
    ele_var, azi_var, tb_var = determine_ds_vars4elev_azi_TB(ds_old)
    
    # 1st determine timeslices with const elev<90 and variable azi.
    time_indices_list_list, BL_time_indices_list_list =\
        determine_scan_slices(ds_old, ele_var=ele_var, azi_var=azi_var)
    
    # 2nd calc mean timestamp and mean measurements for these timeslices
    time_array, tbs, flags, rainfall = determine_data_in_time_for_scanset(ds_old,\
            time_indices_list_list, tb_var=tb_var, ele_var=ele_var,\
            azi_var=azi_var)
    time_array_bl, tbs_bl, flags_bl, rainfall_bl =\
            determine_data_in_time_for_BLset(ds_old, BL_time_indices_list_list,\
            tb_var=tb_var, ele_var=ele_var, azi_var=azi_var)

    # 3rd Create new dataset of scans:
    ds_new = create_scan_ds(time_array, tbs, flags, rainfall, elevations=elevations,\
            azimuths=azimuths, ele_var=ele_var, tb_var=tb_var, azi_var=azi_var)
    ds_bl = create_BL_ds(ds_old,time_array_bl, tbs_bl, flags_bl, rainfall_bl,
                 elevations=elevations,
                 ele_var=ele_var,
                 tb_var=tb_var,
                 out_file=None)

    return ds_new, ds_bl

##############################################################################
# 5th Main code:
##############################################################################
'''
# Claude variant against RAM shortage:
if __name__ == "__main__":
    args = parse_arguments()
    files = sorted(glob.glob(args.in_pattern))
    n = len(files)

    tmp_files = []
    for i, file in enumerate(files):
        print(f"Read file {i} of {n}")
        ds = xr.open_dataset(file)
        ds_resamp, ds_bl = resample_mwr_ds_on_scan_freq(ds)

        ###########
        # Break here when working on timeslice detection:
        # break
        ###########

        ###
        # Leave out empty datasets:
        if ds_resamp.sizes["time"] == 0:
            print(f"  → no valid scans in file {i} — skipping")
            ds.close()
            ds_resamp.close()
            continue
        ###

        ds_resamp = ds_resamp.assign_coords(
            time=ds_resamp["time"].astype("datetime64[ns]"))

        # Direkt als temp-Datei schreiben statt in RAM halten:
        tmp_path = args.outfile + f".tmp_{i:04d}.nc"
        ds_resamp.to_netcdf(tmp_path)
        tmp_files.append(tmp_path)
        ds.close()
        ds_resamp.close()

    # Am Ende lazy einlesen und zusammenfügen:
    print("Concatenating...")
    # ds_final = xr.open_mfdataset(tmp_files, combine="by_coords")
    ds_final = xr.open_mfdataset(tmp_files, combine="nested", concat_dim="time")
    ds_final.to_netcdf(args.outfile)

    # Temp-Dateien aufräumen:
    for tmp in tmp_files:
        os.remove(tmp)

    print("Done:", args.outfile)
'''
if __name__ == "__main__":
    args = parse_arguments()
    files = sorted(glob.glob(args.in_pattern))
    n = len(files)
    tmp_files    = []
    tmp_files_bl = []

    for i, file in enumerate(files):
        print(f"Read file {i} of {n}")
        ds = xr.open_dataset(file)
        ds_resamp, ds_bl = resample_mwr_ds_on_scan_freq(ds)

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
        ds_final_bl.to_netcdf(outfile_bl)
        ds_final_bl.close()
        for tmp in tmp_files_bl:
            os.remove(tmp)
        print("Done:", outfile_bl)
    else:
        print("No BL scan data to write.")


###################################################

