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
import matplotlib

##############################################################################
# 2nd Params:
##############################################################################

matplotlib.use("Agg")
first_col = np.arange(0, 360, 5)
second_col = (first_col + 180) % 360
azi_pairs = np.column_stack([first_col, second_col])

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
        default=os.path.expanduser("~/PhD_data/scans/MWR_scans_JOYCE_Joyhat_202408.nc"),
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
    '''
    trans = np.full((batch_size, 14,10,2), np.nan)
    # time, channel, elevation, crop
    trans_by_lev=np.full((batch_size,n_levels, 14, 10,2 ), np.nan)
    # time, level, channel, elevation, crop
    jacs_by_lev=np.full((batch_size, n_levels, 14, 10,2, 4), np.nan)
    # time, level, channel, elevation, crop, variable (last one removed in ds)
    '''
    
    switch = False
    tb_string = ""
    switch_count = 0
    '''
    switch_t = False
    sw_trans_by_lev=False
    sw_jacs_by_lev=False
    switcht_count = 0
    sw_c_trans_by_lev = 0
    sw_c_jacs_by_lev = 0
    trans_string = ""
    string_jc_by_lev=""
    string_tr_by_lev = ""
    '''

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
            # print("prof_idx: ",prof_idx)
            # print("TBs rs: ", np.array(tbs_rs))

            
    '''
        # Read in Tot Transmittances:     
        if switch_t and switcht_count<2:
            switcht_count+= 1
            trans_string+= line
        elif "CALCULATED SURFACE TO SPACE TRANSMITTANCE:" in line:
            switch_t = True
        elif switch_t:
            switch_t = False
            liste = trans_string.split(" ")
            tot_trans_by_chan = [float(s.strip("\n")) for s in liste if s.strip() != ""]
            trans_string = ""
            switcht_count = 0
            trans[rs_time_idx, :, ele_idx, crop_idx]=np.array(tot_trans_by_chan)
          
        # Read in Lev Transmittances:   
        if sw_trans_by_lev and sw_c_trans_by_lev<n_levels:
            sw_c_trans_by_lev+= 1
            string_tr_by_lev+= line
        elif "Level to surface transmittances for channels" in line:
            sw_trans_by_lev = True
        elif sw_trans_by_lev:
            sw_trans_by_lev = False
            liste = string_tr_by_lev.split("\n")[1:]  #.split(" ")
            for j, line in enumerate(liste):
                list_of_numbers = line.split(" ")
                trans_by_lev1 = [float(s) for s in list_of_numbers if s.strip() != "" and s.strip() != "**"]
                if len(trans_by_lev1)<4:
                    break
                elif 4<=len(trans_by_lev1)<5:
                    trans_by_lev[rs_time_idx, j,10:, ele_idx, crop_idx]=\
                        np.array(trans_by_lev1)
                elif 5<=len(trans_by_lev1)<6:
                    trans_by_lev[rs_time_idx, j,10:, ele_idx, crop_idx] =\
                        np.array(trans_by_lev1[1:])                  
                elif j<99:
                    trans_by_lev[rs_time_idx, j,:10, ele_idx, crop_idx] =\
                        np.array(trans_by_lev1[1:])
                else:
                    trans_by_lev[rs_time_idx, j,:10, ele_idx, crop_idx] =\
                        np.array(trans_by_lev1)            
            string_tr_by_lev = ""
            sw_c_trans_by_lev = 0              
        
        # Read in Jacobians somehow:       
        if "Channel        " in line:
            sw_jacs_by_lev = True
            ch_idx = int(line.split("Channel")[-1])-1
        if sw_jacs_by_lev and sw_c_jacs_by_lev<n_levels+3:
            sw_c_jacs_by_lev+= 1
            string_jc_by_lev+= line
        elif sw_jacs_by_lev:
            liste = string_jc_by_lev.split("\n")[3:n_levels+3]
            for j, line in enumerate(liste):
                werte = line.split()
                jacs_by_lev[rs_time_idx, j, ch_idx, ele_idx, crop_idx, :]=\
                    werte[1:]
            string_jc_by_lev=""
            sw_jacs_by_lev = False
            sw_c_jacs_by_lev= 0
    '''
            
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

def calc_estimated_tilt(azis,dtbs, pair_labels, args, elevation=30,\
        mask=None, azi_pairs=azi_pairs):

    ###
    # 1st round: 
    tilts = np.arange(0.01, 2, 0.25) # 20 => ca 8.
    all_dtbs_mod = np.full((len(azi_pairs), len(tilts), len(dtbs), 14), np.nan)
    sqr_diffs = np.full((len(azi_pairs), len(tilts)), np.nan)
    for i, pair in enumerate(azi_pairs):
        for j, i_tilt in enumerate(tilts):
            dtbs_mod = calc_TB_diffs_360deg4tilt(azis,pair,args,\
                    ang0=elevation, tilt=i_tilt)
            all_dtbs_mod[i,j,:,:] = dtbs_mod
            sqr_diffs[i,j] = np.nanmean((dtbs[:,mask]-dtbs_mod[:,mask])**2)

    ###
    #2nd round:
    flat_index = np.nanargmin(sqr_diffs)
    i_pair, i_tilt = np.unravel_index(flat_index, sqr_diffs.shape)
    tilt       = tilts[i_tilt]
    angle_pair = azi_pairs[i_pair]
    print("First tilt: ", tilt)
    print("First angle: ", angle_pair )
    start = np.max([tilt-0.25, 0])
    star_index_pair = np.max([0, i_pair-1])
    new_pairs = azi_pairs[star_index_pair : i_pair+1]
    tilts = np.arange(start, tilt+0.25, 0.01)
    print("tilts: ", tilts)
    print("New pairs: ", new_pairs)
    all_dtbs_mod = np.full((len(new_pairs), len(tilts), len(dtbs), 14), np.nan)
    sqr_diffs = np.full((len(new_pairs), len(tilts)), np.nan)
    for i, pair in enumerate(new_pairs):
        for j, i_tilt in enumerate(tilts):
            dtbs_mod = calc_TB_diffs_360deg4tilt(azis,pair,args,\
                    ang0=elevation, tilt=i_tilt)
            all_dtbs_mod[i,j,:,:] = dtbs_mod
            sqr_diffs[i,j] = np.nanmean((dtbs[:,mask]-dtbs_mod[:,mask])**2)

    # Get minimum difference between tilted model and measured data:
    flat_index = np.nanargmin(sqr_diffs)
    i_pair, i_tilt = np.unravel_index(flat_index, sqr_diffs.shape)
    tilt       = tilts[i_tilt]
    angle_pair = new_pairs[i_pair]

    '''
    ############################################
    tilts = []
    sqr_diffs = []
    i_tilt = 0.3
    d_tilt = 0.001
    all_dtbs_mod = np.full((len(azi_pairs), len(dtbs), 14), np.nan)
    for i, pair in enumerate(azi_pairs):
        for j in range(50):
            dtbs_mod = calc_TB_diffs_360deg4tilt(azis,pair,args,\
                    ang0=elevation, tilt=i_tilt+d_tilt)
            sqr_diff2 = np.nanmean((dtbs[:,mask]-dtbs_mod[:,mask])**2)

            dtbs_mod = calc_TB_diffs_360deg4tilt(azis,pair,args,\
                    ang0=elevation, tilt=i_tilt)
            sqr_diff1 = np.nanmean((dtbs[:,mask]-dtbs_mod[:,mask])**2)
            
            if (sqr_diff1 / ((sqr_diff2-sqr_diff1)/d_tilt))<=d_tilt:
                print("Found zero for pair - Iterations:", j)
                break
            else:
                print("Iteration: ", j)
                print("sqr_diff: ", sqr_diff1)
                print("tilt: ", i_tilt)
                print("Gradient: ", ((sqr_diff2-sqr_diff1)/d_tilt))
                i_tilt = i_tilt - sqr_diff1 / ((sqr_diff2-sqr_diff1)/d_tilt)
            ##########################
            # gradient = (sqr_diff2 - sqr_diff1) / d_tilt
            # i_tilt = i_tilt - learning_rate * gradient   
            ###########################
        all_dtbs_mod[i,:,:] = dtbs_mod
        sqr_diffs.append(sqr_diffs)
        tilts.append(i_tilt)

    # Get minimum difference between tilted model and measured data:
    i_pair = np.nanargmin(sqr_diffs)
    tilt       = tilts[i_pair]
    angle_pair = azi_pairs[i_pair]
    ###########################################
        
    
    # Get minimum difference between tilted model and measured data:
    flat_index = np.nanargmin(sqr_diffs)
    i_pair, i_tilt = np.unravel_index(flat_index, sqr_diffs.shape)
    tilt       = tilts[i_tilt]
    angle_pair = azi_pairs[i_pair]
    '''
    return tilt, angle_pair, all_dtbs_mod[i_pair, i_tilt ,:,:]

##############################################################################

def get_channel_mask(dtbs):
    maxima = np.nanmax(dtbs, axis=0)
    minima = np.nanmin(dtbs, axis=0)
    idx_high = np.argsort(maxima)[-2:]   # die 2 größten
    idx_low  = np.argsort(minima)[:2]    # die 2 kleinsten
    exclude = np.union1d(idx_high, idx_low)
    mask = np.ones(dtbs.shape[1], dtype=bool)
    mask[exclude] = False
    return mask

##############################################################################
# 5th Main code:
##############################################################################

if __name__=="__main__":    
    args = parse_arguments()
    ds = xr.open_dataset(args.infile)

    i_elev=0 # 30 ° Elevation!
    azis,dtbs, pair_labels = get_TB_opposite_differences(ds, i_elev=i_elev)
    mask = get_channel_mask(dtbs)
    tilt, angle_pair, dtbs_mod = calc_estimated_tilt(azis,dtbs, pair_labels,\
        args, elevation=ds["elevation"].values[i_elev],mask=mask)

    # Create a comparison plot:
    grad = azis
    theta = np.deg2rad(grad)  # Grad → Radiant
    plt.figure(figsize=(15,10))
    plt.title(f"Detected 360° Azimuth instrument tilt: Angle {tilt:.2f}°"
          f" Between high and low angle: {angle_pair}°")
    plt.plot(theta[::1], np.nanmean(dtbs[:, mask], axis=1), label="Instrument")
    plt.plot(theta[::1], np.nanmean(dtbs_mod[:,mask],axis=1), label="Model")
    plt.xticks(theta[::5], pair_labels[::5], rotation=45, ha='right')
    plt.legend()
    plt.savefig(os.path.expanduser("~/tilt_plot.png"))

    # Datafile:
    with open(os.path.expanduser("~/tilt_detection.txt"), "w") as f:
        f.writelines("Tilt: "+str(tilt)+"°\n") 
        f.writelines("Angle_pair: "+str(angle_pair)+"°\n")
        f.writelines("TB differences: "+str(dtbs_mod))

    #################

    # Tilt timeseries!!!

    # CLEANER DATASETS
    # 1. Only use clear sky TBs!!!
    # Exclude RFI Interferences!!! (Tobis Paper???


    # CORRECTION!!!
    # For correction i probably have to save the dtbs_mod to subtract them from
    # a measurement...
    # 0th add corrected dataset with correction mask NetCDF...


    ##############
    # 2. What can I do to exclude general Water Vapor features of the landscape...?
    # 3. get another profile from a radiosonde instead of pure clim!
    # 4. Add newtons method to determine tilt value instead of monte carlo!























