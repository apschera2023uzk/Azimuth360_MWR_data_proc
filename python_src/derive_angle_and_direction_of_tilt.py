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

##############################################################################
# 2nd Params:
##############################################################################

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
    '''
    parser.add_argument(
        "--outfile", "-o",
        type=str,
        default=os.path.expanduser("~/PhD_data/scans/MWR_scans_JOYCE_Joyhat_202408.nc"),
        help="NetCDF Output file path."
    )
    '''
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
    string+=f"{zenith_angle:6.1f}\n"
        
    return string

##############################################################################

def write_combined_input_prof_file(t_array, ppmv_array,length_value, p_array,height_in_km=0., deg_lat=50.,\
                                   filename="prof_plev.dat", zenith_angle=0.):
    with open(filename, "w") as file:
        # print("pressure levels: ", length_value)
        for value in p_array:
            file.write(f"{value:8.4f}\n")  # eingerückt, 4 Nachkommastellen
        for value in t_array:
            file.write(f"{value:6.3f}\n")  # eingerückt, 4 Nachkommastellen
        for value in ppmv_array:
            file.write(f"{value:9.4f}\n")  # eingerückt, 4 Nachkommastellen
        for value in range(len(ppmv_array)):
            file.write(f"{0.:12.6E}\n")  # eingerückt, 4 Nachkommastellen
        file.write(f"{t_array[-11]:10.4f}{p_array[-1]:10.2f}\n")
        file.write(f"{height_in_km:6.1f}{deg_lat:6.1f}\n")
        file.write(f"{zenith_angle:6.1f}\n")
        
        return 0

##############################################################################

def get_radisonde_tbs_of_file(rttovgb_outfile):
    # Diese Funktion läuft einwandfrei.
    switch = False
    switch_count = 0
    tb_string = ""
    file = open(rttovgb_outfile, "r")
    for i, line in enumerate(file.readlines()):
        if switch and switch_count<2:
            switch_count+= 1
            tb_string+= line
        elif "CALCULATED BRIGHTNESS TEMPERATURES (K):" in line:
            switch = True
    liste = tb_string.split(" ")
    tbs_rs = [float(s.strip("\n")) for s in liste if s.strip() != ""]
    
    file.close()
    return np.array(tbs_rs)
    
##############################################################################

def calc_TBs4_one_case(ang=90, model="RTTOV-gb"):
    # Determine prams:
    z, p, d, t, md = atmp.gl_atm(atm=1) # midlatitude summer!
    gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
    rh = mr2rh(p, t, gkg)[0] / 100      

    if model=="RTTOV-gb":
        z_in = z*1000
        ppmv = np.full((len(rh)), np.nan)
        for i, rh1 in enumerate(rh):
            ppmv[i] = rh2ppmv(RH=rh1*100, abs_T=t[i], p=p[i]*100)

            
        write_combined_input_prof_file(t[::-1], ppmv[::-1],len(ppmv), p[::-1],\
                    height_in_km=z_in[1], deg_lat=50.,\
                    filename="prof_plev.dat", zenith_angle=(90-ang))


        nlevels=len(ppmv)
        shutil.copy("prof_plev.dat", "/home/aki/RTTOV-gb/rttov_test/test_example_k.1/")
        with open("/home/aki/RTTOV-gb/rttov_test/run_apschera.sh", "r") as f:
            lines = f.readlines()
            # Zeile 30 (Index 29) ersetzen
        lines[28] = f"NPROF=1\n"    
        lines[29] = f"NLEVELS={nlevels}\n"
        with open("/home/aki/RTTOV-gb/rttov_test/run_apschera.sh", "w") as f:
            f.writelines(lines) 
        subprocess.run(["bash", "/home/aki/RTTOV-gb/rttov_test/run_apschera.sh", "ARCH=gfortran-openmp"],\
                           cwd="/home/aki//RTTOV-gb/rttov_test/") 
        
        tbs = get_radisonde_tbs_of_file("/home/aki/RTTOV-gb/rttov_test/test_example_k.1/output_example_k.dat.gfortran")
    
    return tbs

###############################################################################

def calc_TBs4_72_cases(angs=np.array([30]*len(azi_pairs)), model="RTTOV-gb"):
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
                    height_in_km=z_in[1], deg_lat=50.,\
                zenith_angle=(90-ang), clear_sky_bool=True)
            profiles+=profile1

        # After loops - save results:
        dir_name = os.path.dirname(args.input)
        outfile = "prof_plev.dat"
        out = open(outfile, "w")
        out.write(profiles)
        out.close()

        ########################
        # Here I was just working....
        # Trying to create an Input file and use it that has 72 entries instead of 1...
        # I think the file exists now... It just needs to be processed in RTTOV-gb
        # And the results have to be returned in the expected way in the function above...


            
###############################################################################

def calc_TB_diffs_360deg4tilt(azis,pair, ang0=30, tilt=0.5):
    
    ########################################
    # Calculate TBs for all azimuth angles:
    '''
    tbs = np.full((len(azis), 14), np.nan)
    for i, azi in enumerate(azis):
        grad = azi-pair[0]
        ang = ang0+tilt*np.cos(np.deg2rad(grad))
        tbs[i,:] = calc_TBs4_one_case(ang=ang)
    '''
    #######################################
    # Alternate:
    angs = []
    for i, azi in enumerate(azis):
        grad = azi-pair[0]
        angs.append(ang0+tilt*np.cos(np.deg2rad(grad)))
    #######################################
    
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

def calc_estimated_tilt(azis,dtbs, pair_labels, elevation=30, azi_pairs=azi_pairs):

    tilts = np.arange(0, 2, 0.1) # 200
    sqr_diffs = np.full((len(azi_pairs), len(tilts)), np.nan)
    for i, pair in enumerate(azi_pairs):
        for j, i_tilt in enumerate(tilts):
            dtbs_mod = calc_TB_diffs_360deg4tilt(azis,pair, ang0=elevation, tilt=i_tilt)
            ######################
            sqr_diffs[i,j] = np.nanmean((dtbs-dtbs_mod)**2) # enthalten NaNs???
            ########################
                        
    # Get minimum difference between tilted model and measured data:
    flat_index = np.nanargmin(sqr_diffs)
    i_pair, i_tilt = np.unravel_index(flat_index, sqr_diffs.shape)
    tilt       = tilts[i_tilt]
    angle_pair = azi_pairs[i_pair]
    
    return tilt, angle_pair

##############################################################################
# 5th Main code:
##############################################################################

if __name__=="__main__":    
    args = parse_arguments()
    ds = xr.open_dataset(args.infile)

    i_elev=0 # 30 ° Elevation!
    azis,dtbs, pair_labels = get_TB_opposite_differences(ds, i_elev=i_elev)
    tilt, angle_pair = calc_estimated_tilt(azis,dtbs, pair_labels,\
                                    elevation=ds["elevation"].values[i_elev])

    #################
    # Possible improvements:
    # 1st Handle data from several angles (the full 360 deg azimuth scan)
    # together in one prof_plev.dat file (less reading and writing)
    # 2nd Create a plot at the end of the run of the modelled TB sin wave 
    # and data TB sin wave.
    # 2.5th Write an output file with all data! (Input to next step: Correction!)
    # 3rd: Maybe exclude channels which strongly deviate from all others
    # 4th What can I do to exclude general Water Vapor features of the landscape...?























