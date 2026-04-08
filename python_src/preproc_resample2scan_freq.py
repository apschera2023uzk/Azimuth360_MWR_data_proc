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
    time_array = []
    
    # All scans within MWR file:
    for i, timeslice in enumerate(time_indices_list_list):
        
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
    
    return time_array, tbs

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

def create_scan_ds(time_array, tbs, elevations=elevations, azimuths=azimuths,\
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
    time_array, tbs = determine_data_in_time_for_scanset(ds_old,\
            time_indices_list_list, tb_var=tb_var, ele_var=ele_var,\
            azi_var=azi_var)

    # 3rd Create new dataset of scans:
    ds_new = create_scan_ds(time_array, tbs, elevations=elevations,\
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
'''
netcdf MWR_1C01_XXX_20240805 {
dimensions:
	time = 57743 ;
	frequency = 14 ;
	receiver_nb = 2 ;
	ir_wavelength = 1 ;
	bnds = 2 ;
	t_amb_nb = 2 ;
variables:
	int time(time) ;
		time:_FillValue = -2147483647 ;
		time:units = "seconds since 1970-01-01 00:00:00.000" ;
		time:long_name = "Time (UTC) of the measurement" ;
		time:comment = "Time indication of samples is at end of integration-time" ;
	int time_bnds(time, bnds) ;
		time_bnds:_FillValue = -2147483647 ;
		time_bnds:units = "seconds since 1970-01-01 00:00:00.000" ;
		time_bnds:long_name = "Start and end time (UTC) of the measurements" ;
	float latitude(time) ;
		latitude:_FillValue = 9.96921e+36f ;
		latitude:units = "degrees_north" ;
		latitude:long_name = "Latitude of measurement station" ;
		latitude:standard_name = "latitude" ;
	float longitude(time) ;
		longitude:_FillValue = 9.96921e+36f ;
		longitude:units = "degrees_east" ;
		longitude:long_name = "Longitude of measurement station" ;
		longitude:standard_name = "longitude" ;
	float altitude(time) ;
		altitude:_FillValue = 9.96921e+36f ;
		altitude:units = "m" ;
		altitude:long_name = "Altitude above mean sea level of measurement station" ;
		altitude:standard_name = "altitude" ;
	float frequency(frequency) ;
		frequency:_FillValue = 9.96921e+36f ;
		frequency:units = "GHz" ;
		frequency:long_name = "Nominal centre frequency of microwave channels" ;
		frequency:standard_name = "radiation_frequency" ;
		frequency:comment = "1) For double-sideband receivers, frequency corresponds to the\n",
			"local oscillator frequency whereas the radio frequency of the upper/lower\n",
			"sideband is frequency+/-sideband_IF_separation. 2) In case of known\n",
			"offset between the real and the nominal frequency of some channels,\n",
			"frequency+freq_shift gives more accurate values." ;
	int receiver_nb(receiver_nb) ;
		receiver_nb:_FillValue = -2147483647 ;
		receiver_nb:units = "1" ;
		receiver_nb:long_name = "Microwave receiver number" ;
	int receiver(frequency) ;
		receiver:_FillValue = -2147483647 ;
		receiver:units = "1" ;
		receiver:long_name = "Corresponding microwave receiver for each channel" ;
	float bandwidth(frequency) ;
		bandwidth:_FillValue = 9.96921e+36f ;
		bandwidth:units = "GHz" ;
		bandwidth:long_name = "Bandwidth of microwave channels" ;
	int n_sidebands(receiver_nb) ;
		n_sidebands:_FillValue = -2147483647 ;
		n_sidebands:units = "1" ;
		n_sidebands:long_name = "Number of sidebands" ;
		n_sidebands:comment = "0: direct-detection receivers, 1: single-sideband,\n",
			"2: double-sideband. The frequency separation of sidebands\n",
			"is indicated in sideband_IF_separation." ;
	float sideband_IF_separation(frequency) ;
		sideband_IF_separation:_FillValue = 9.96921e+36f ;
		sideband_IF_separation:units = "GHz" ;
		sideband_IF_separation:long_name = "Sideband IF separation" ;
		sideband_IF_separation:comment = "For double sideband channels, this is the positive and negative\n",
			"IF range distance of the two band passes around the centre frequency\n",
			"(which is the LO frqeuency)" ;
	float freq_shift(frequency) ;
		freq_shift:_FillValue = 9.96921e+36f ;
		freq_shift:units = "GHz" ;
		freq_shift:long_name = "Frequency shift of the microwave channels" ;
		freq_shift:comment = "For more accurate frequency values use frequency + freq_shift." ;
	float tb(time, frequency) ;
		tb:_FillValue = 9.96921e+36f ;
		tb:units = "K" ;
		tb:long_name = "Microwave brightness temperature" ;
		tb:standard_name = "brightness_temperature" ;
	float azimuth_angle(time) ;
		azimuth_angle:_FillValue = 9.96921e+36f ;
		azimuth_angle:units = "degree" ;
		azimuth_angle:long_name = "Azimuth angle" ;
		azimuth_angle:standard_name = "sensor_azimuth_angle" ;
		azimuth_angle:comment = "0=North, 90=East, 180=South, 270=West" ;
	float elevation_angle(time) ;
		elevation_angle:_FillValue = 9.96921e+36f ;
		elevation_angle:units = "degree" ;
		elevation_angle:long_name = "Sensor elevation angle" ;
		elevation_angle:comment = "0=horizon, 90=zenith" ;
	int quality_flag(time, frequency) ;
		quality_flag:_FillValue = -2147483647 ;
		quality_flag:units = "1" ;
		quality_flag:long_name = "Quality flag" ;
		quality_flag:definition = "\n",
			"Bit 1: missing_tb\n",
			"Bit 2: tb_below_threshold\n",
			"Bit 3: tb_above_threshold\n",
			"Bit 4: spectral_consistency_above_threshold\n",
			"Bit 5: receiver_sanity_failed\n",
			"Bit 6: rain_detected\n",
			"Bit 7: sun_moon_in_beam\n",
			"Bit 8: tb_offset_above_threshold" ;
		quality_flag:comment = "0 indicates data with good quality according to applied tests.\n",
			"The list of (not) applied tests is encoded in quality_flag_status" ;
	int quality_flag_status(time, frequency) ;
		quality_flag_status:_FillValue = -2147483647 ;
		quality_flag_status:units = "1" ;
		quality_flag_status:long_name = "Quality flag status" ;
		quality_flag_status:definition = "\n",
			"Bit 1: missing_tb_not_checked\n",
			"Bit 2: tb_lower_threshold_not_checked\n",
			"Bit 3: tb_upper_threshold_not_checked\n",
			"Bit 4: spectral_consistency_not_checked\n",
			"Bit 5: receiver_sanity_not_checked\n",
			"Bit 6: rain_not_checked\n",
			"Bit 7: sun_moon_in_beam_not_checked\n",
			"Bit 8: tb_offset_not_checked" ;
		quality_flag_status:comment = "Checks not executed in determination of quality_flag.\n",
			"0 indicates quality check has been applied." ;
	int liquid_cloud_flag(time) ;
		liquid_cloud_flag:_FillValue = -2147483647 ;
		liquid_cloud_flag:units = "1" ;
		liquid_cloud_flag:long_name = "Liquid cloud flag" ;
		liquid_cloud_flag:comment = "Flag meaning: no liquid cloud (0), liquid cloud present (1),\n",
			"undefined (2)" ;
	int liquid_cloud_flag_status(time) ;
		liquid_cloud_flag_status:_FillValue = -2147483647 ;
		liquid_cloud_flag_status:units = "1" ;
		liquid_cloud_flag_status:long_name = "Liquid cloud flag status" ;
		liquid_cloud_flag_status:comment = "Flag meaning: using mwr and ir (0), using mwr only (1), other (2)" ;
	int pointing_flag(time) ;
		pointing_flag:_FillValue = -2147483647 ;
		pointing_flag:units = "1" ;
		pointing_flag:long_name = "Pointing flag" ;
		pointing_flag:comment = "Flag indicating a single pointing (staring = 0)\n",
			"or multiple pointing (scanning = 1) observation sequence" ;
	float t_amb(time, t_amb_nb) ;
		t_amb:_FillValue = 9.96921e+36f ;
		t_amb:units = "K" ;
		t_amb:long_name = "Ambient target temperature" ;
	float t_rec(time, receiver_nb) ;
		t_rec:_FillValue = 9.96921e+36f ;
		t_rec:units = "K" ;
		t_rec:long_name = "Receiver physical temperature" ;
	float t_sta(time, receiver_nb) ;
		t_sta:_FillValue = 9.96921e+36f ;
		t_sta:units = "K" ;
		t_sta:long_name = "Receiver temperature stability" ;
	float tb_spectrum(time, frequency) ;
		tb_spectrum:_FillValue = 9.96921e+36f ;
		tb_spectrum:units = "K" ;
		tb_spectrum:long_name = "Retrieved brightness temperature spectrum" ;
	float ir_wavelength(ir_wavelength) ;
		ir_wavelength:_FillValue = 9.96921e+36f ;
		ir_wavelength:units = "m" ;
		ir_wavelength:long_name = "Wavelength of infrared channels" ;
		ir_wavelength:standard_name = "sensor_band_central_radiation_wavelength" ;
	float ir_bandwidth ;
		ir_bandwidth:_FillValue = 9.96921e+36f ;
		ir_bandwidth:units = "m" ;
		ir_bandwidth:long_name = "Bandwidth of infrared channels" ;
		ir_bandwidth:comment = "Channel centre frequency." ;
	float ir_beamwidth ;
		ir_beamwidth:_FillValue = 9.96921e+36f ;
		ir_beamwidth:units = "degree" ;
		ir_beamwidth:long_name = "Beam width of the infrared radiometer" ;
	float irt(time, ir_wavelength) ;
		irt:_FillValue = 9.96921e+36f ;
		irt:units = "K" ;
		irt:long_name = "Infrared brightness temperatures" ;
	float ir_azimuth_angle(time) ;
		ir_azimuth_angle:_FillValue = 9.96921e+36f ;
		ir_azimuth_angle:units = "degree" ;
		ir_azimuth_angle:long_name = "Infrared sensor azimuth angle" ;
		ir_azimuth_angle:standard_name = "sensor_azimuth_angle" ;
		ir_azimuth_angle:comment = "0=North, 90=East, 180=South, 270=West" ;
	float ir_elevation_angle(time) ;
		ir_elevation_angle:_FillValue = 9.96921e+36f ;
		ir_elevation_angle:units = "degree" ;
		ir_elevation_angle:long_name = "Infrared sensor elevation angle" ;
		ir_elevation_angle:comment = "0=horizon, 90=zenith" ;
	float air_temperature(time) ;
		air_temperature:_FillValue = 9.96921e+36f ;
		air_temperature:units = "K" ;
		air_temperature:long_name = "Air temperature" ;
		air_temperature:standard_name = "air_temperature" ;
	float relative_humidity(time) ;
		relative_humidity:_FillValue = 9.96921e+36f ;
		relative_humidity:units = "1" ;
		relative_humidity:long_name = "Relative humidity" ;
		relative_humidity:standard_name = "relative_humidity" ;
	float air_pressure(time) ;
		air_pressure:_FillValue = 9.96921e+36f ;
		air_pressure:units = "Pa" ;
		air_pressure:long_name = "Air pressure" ;
		air_pressure:standard_name = "air_pressure" ;
	float rainfall_rate(time) ;
		rainfall_rate:_FillValue = 9.96921e+36f ;
		rainfall_rate:units = "m s-1" ;
		rainfall_rate:long_name = "Rainfall rate" ;
		rainfall_rate:standard_name = "rainfall_rate" ;
	float wind_direction(time) ;
		wind_direction:_FillValue = 9.96921e+36f ;
		wind_direction:units = "degree" ;
		wind_direction:long_name = "Wind direction" ;
		wind_direction:standard_name = "wind_from_direction" ;
	float wind_speed(time) ;
		wind_speed:_FillValue = 9.96921e+36f ;
		wind_speed:units = "m s-1" ;
		wind_speed:long_name = "Wind speed" ;
		wind_speed:standard_name = "wind_speed" ;
	int met_quality_flag(time) ;
		met_quality_flag:_FillValue = -2147483647 ;
		met_quality_flag:units = "1" ;
		met_quality_flag:long_name = "Meteorological data quality flag" ;
		met_quality_flag:definition = "\n",
			"Bit 1: low_quality_air_temperature\n",
			"Bit 2: low_quality_relative_humidity\n",
			"Bit 3: low_quality_air_pressure\n",
			"Bit 4: low_quality_rainfall_rate\n",
			"Bit 5: low_quality_wind_direction\n",
			"Bit 6: low_quality_wind_speed" ;
		met_quality_flag:comment = "0=ok, 1=problem. Note: should also be set to 1\n",
			"if corresponding sensor not available" ;
'''
# Quality flags from here????
