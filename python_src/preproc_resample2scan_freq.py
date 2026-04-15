#! /usr/bin/env python3

##############################################################################
# 1st Necessary Modules
##############################################################################

import xarray as xr
import glob
import numpy as np
import argparse
import os

##############################################################################
# 2nd Params:
##############################################################################

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
        default=os.path.expanduser("~/PhD_data/scans/joyhat_raw_jun_jul_aug_sep/MWR_1C01_*.nc"),
        # default=os.path.expanduser("~/PhD_data/FESSTVaL_14GB/foghat/l1/*/*/fval_uzk_mwr00_l1_tb*.nc"),
        help="Pattern of MWR output files with TBs of scans."
    )
    parser.add_argument(
        "--outfile", "-o",
        type=str,
        default=os.path.expanduser("~/PhD_data/scans/MWR_scans_JOYCE_Joyhat_202408.nc"),
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
                          max_elev_azi_diff=max_elev_azi_diff):
    # Output: list of lists, because each scan is one list of indices!
    time_indices_list_list = []
    scan_switch = False
    ele_old = 90
    
    for i, timestep in enumerate(ds_in["time"].values):
        ele_curr = ds_in[ele_var].values[i]

        if scan_switch and ele_curr<90:
            if abs(ele_curr-ele_old)<max_elev_azi_diff:
                scan_list.append(i)
            else:
                scan_switch = False
                if len(scan_list)>3:
                    time_indices_list_list.append(scan_list)
        elif ele_curr<90:
            scan_switch = True
            scan_list = []
            scan_list.append(i)
        elif abs(ele_curr-ele_old)>max_elev_azi_diff and scan_switch:
            scan_switch = False
            if len(scan_list)>3:
                time_indices_list_list.append(scan_list)
        else:
            continue
            
        ele_old = ds_in[ele_var].values[i]
        azi_old = ds_in[azi_var].values[i]
    
    return time_indices_list_list

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
            flag_vals_valid = flag_vals[flag_vals != fill].astype(int)
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

def resample_mwr_ds_on_scan_freq(ds_old):
    
    # 0th Determine azi, ele and tb vars:
    ele_var, azi_var, tb_var = determine_ds_vars4elev_azi_TB(ds_old)
    
    # 1st determine timeslices with const elev<90 and variable azi.
    time_indices_list_list = determine_scan_slices(ds_old,\
                        ele_var=ele_var, azi_var=azi_var)
    
    # 2nd calc mean timestamp and mean measurements for these timeslices
    time_array, tbs, flags, rainfall = determine_data_in_time_for_scanset(ds_old,\
            time_indices_list_list, tb_var=tb_var, ele_var=ele_var,\
            azi_var=azi_var)

    # 3rd Create new dataset of scans:
    ds_new = create_scan_ds(time_array, tbs, flags, rainfall, elevations=elevations,\
            azimuths=azimuths, ele_var=ele_var, tb_var=tb_var, azi_var=azi_var)

    return ds_new

##############################################################################
# 5th Main code:
##############################################################################

if __name__=="__main__":    
    args = parse_arguments()
    files = glob.glob(args.in_pattern)
    n = len(files)
    ds_list = []

    for i, file in enumerate(files):
        print("Read file ", i, " of ", n)
        ds = xr.open_dataset(file)
        ds_resamp = resample_mwr_ds_on_scan_freq(ds)
        ds_list.append(ds_resamp)

    # Fix 1 — alle time-Koordinaten auf ns casten vor dem concat:
    ds_list = [ds.assign_coords(time=ds["time"].astype("datetime64[ns]"))
               for ds in ds_list]
    ds_joy = xr.concat(ds_list, dim="time")

    print(ds_joy)
    print(ds_joy["time"])

    ds_joy.to_netcdf(args.outfile)#, format="NETCDF4_CLASSIC")


###################################################

