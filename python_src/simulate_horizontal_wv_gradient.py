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
fs = 20
plt.rc('font', size=fs) 
plt.style.use('seaborn-poster')
matplotlib.use("Agg")
first_col = np.arange(0, 360, 5)
second_col = (first_col + 180) % 360
azi_pairs = np.column_stack([first_col, second_col])
azimuths = np.arange(0.,355.1,5.)
elevations = np.array([90., 30, 19.2, 14.4, 11.4, 8.4,  6.6,  5.4, 4.8,  4.2])

##############################################################################
# 3rd Argparse
##############################################################################

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="This script derives the instrument tilt and direction from a NetCDF file of TBs of 360° Azimuth scans at scan frequency."
    )
    # ./derive_angle_and_direction_of_tilt.py --infile /home/aki/PhD_data/scans/MWR_scans_RAO_Foghat_202105_08.nc
    parser.add_argument(
        "--rttov", "-rt",
        type=str,
        default=os.path.expanduser("~/RTTOV-gb"),
        help="Runscript for RTTOV-gb."
    )
    parser.add_argument(
        "--plots", "-p",
        type=str,
        default=os.path.expanduser("~/PhD_plots/2026"),
        help="Where to save plots."
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

##############################################################################

def ensure_folder_exists(base_path, folder_name):
    # Join the base path with the folder name
    folder_path = os.path.join(base_path, folder_name)
    # Create the directory if it does not exist
    os.makedirs(folder_path, exist_ok=True)

    return os.path.abspath(folder_path)

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

def determine_x_value(azi, elevation, z):
    r = z / np.tan(np.deg2rad(elevation))     
    x = r * np.sin(np.deg2rad(azi))
    return x

###############################################################################

def wv_factor4x(x, maxfac=1.15, minfac=0.85, distance=40000):
    # linear increase and constant values on the borders!!!
    m = (maxfac-minfac)/distance
    b = minfac + distance*m/2
    if x>distance/2:
        factor = maxfac
    elif x<-distance/2:
        factor = minfac
    else:
        factor = m * x + b
    return factor

###############################################################################

def calc_TBs4_72_cases(args, azis=azimuths, elevation=10, model="RTTOV-gb",\
        maxfac=1.15, minfac=0.85):

    # Determine params:
    z, p, d, t, md = atmp.gl_atm(atm=1) # midlatitude summer!
    gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
    rh = mr2rh(p, t, gkg)[0] / 100      

    if model=="RTTOV-gb":
        z_in = z*1000
        ppmv = np.full((len(rh), len(azis)), np.nan)
        for i, (rh1, z1) in enumerate(zip(rh,z_in)):
            for j, azi in enumerate(azis): 
                # scale ppmv here depending on angle and distance!!!
                x = determine_x_value(azi, elevation, z1)
                factor = wv_factor4x(x, maxfac=maxfac, minfac=minfac)
                ppmv[i,j] = rh2ppmv(RH=rh1*100, abs_T=t[i], p=p[i]*100)*factor

                '''
                print("***************")
                print("azi: ", azi)
                print("z => x", z1, " => ", x)
                print("factor: ", factor)
                print("ppmv: ", ppmv[i,j])
                '''
        # Processing for 72 Azimuths
        profiles = ""
        for j, azi in enumerate(azis):
            profile1 = write1profile2str(t[::-1], ppmv[::-1,j],len(z_in), p[::-1],\
                    np.array([0.]*len(ppmv)), height_in_km=z_in[1], deg_lat=50.,\
                zenith_angle=(90-elevation), clear_sky_bool=True)
            profiles+=profile1

        # After loops - save results:
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
        lines[28] = f"NPROF="+str(len(azis))+"\n"    
        lines[29] = f"NLEVELS={nlevels}\n"
        with open(args.rttov+"/rttov_test/run_apschera.sh", "w") as f:
            f.writelines(lines) 
        subprocess.run(["bash", args.rttov+"/rttov_test/run_apschera.sh",\
                "ARCH=gfortran-openmp"], cwd=args.rttov+"/rttov_test/") 

        # Read outputs:
        tbs = get_rttov_outputs(rttovgb_outfile=\
            args.rttov+"/rttov_test/test_example_k.1/output_example_k.dat.gfortran-openmp",\
            n_azis=len(azis),n_levels=nlevels)

    return tbs

##############################################################################

def plot_wv_function(args, maxfac=1.15, minfac=0.85):
    base_path = args.plots
    outpath = ensure_folder_exists(base_path,"sim_hor_grad")

    xs = np.linspace(-40000,40000,100)
    ys = []
    for x in xs:
        ys.append( wv_factor4x(x, maxfac=maxfac, minfac=minfac))

    plt.figure()
    plt.title(f"{(maxfac-1)*100:.0f} %)")
    plt.plot(xs,ys)
    plt.tight_layout()
    plt.savefig(outpath+f"/x_factor_{(maxfac-1)*100:.0f}_percent.png",
                dpi=150)

    return 0

##############################################################################

def plot_wv_colorplot(args, lvl=1, resol=100, maxfac=1.15, minfac=0.85):
    base_path = args.plots
    outpath = ensure_folder_exists(base_path,"sim_hor_grad")

    # Determine params:
    z, p, d, t, md = atmp.gl_atm(atm=1) # midlatitude summer!
    gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)      
    xs = np.linspace(-40000,40000,resol)
    ys = np.linspace(-40000,40000,resol)
    facs = []
    for x in xs:
        facs.append(wv_factor4x(x, maxfac=maxfac, minfac=minfac))

    c = np.array(facs) * float(gkg[lvl])
    z_title = z[lvl]*1000
    c_array = np.tile(c,(resol,1))
    vmin = np.nanmin(c_array)
    vmax = np.nanmax(c_array)
    ticks = np.linspace(vmin, vmax, 6)  # 6 gleichmäßige Ticks

    ###
    # Lineplot - somehow gkg is not just one array....:    
    plt.figure
    plt.plot(xs, c)
    plt.title(f"Water Vapor gradient in g/kg\n WV perturbation: {(maxfac-1)*100:.0f} %")
    plt.axvline(0, color="black")
    plt.ylim(6, 11)
    plt.ylabel("WV mixing ratio [g/kg]")
    plt.xlabel("x [m]")
    plt.tight_layout()
    plt.savefig(outpath+f"/lineplot_wv_{(maxfac-1)*100:.0f}_percent.png",
                dpi=150)

    # Colorplot:
    plt.figure(figsize=(15,10))
    pm =  plt.pcolormesh(xs, ys, c_array, vmin=vmin, vmax=vmax)
    cb = plt.colorbar(pm, label="WV mixing ratio [g/kg]")
    cb.set_ticks(ticks)
    cb.set_ticklabels([f"{t:.2f}" for t in ticks])
    plt.axhline(0, color="black")
    plt.axvline(0, color="black")
    plt.title(f"Horizontal water vapor gradient (z = {z_title} m)\n WV perturbation: {(maxfac-1)*100:.0f} %")
    plt.scatter([0],[0], marker="X", color="red", linewidth=3,\
        label="MWR position")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.tight_layout()
    plt.savefig(outpath+f"/colorplot_wv_{(maxfac-1)*100}_percent.png",
                dpi=150)

    # DIAL-Colorplot:
    plt.figure(figsize=(15,10))
    dial_array = c_array
    dial_array[:,:] = c_array[50,50]
    pm =  plt.pcolormesh(xs, ys, dial_array, vmin=vmin, vmax=vmax)
    cb = plt.colorbar(pm, label="WV mixing ratio [g/kg]")
    cb.set_ticks(ticks)
    cb.set_ticklabels([f"{t:.2f}" for t in ticks])
    plt.axhline(0, color="black")
    plt.axvline(0, color="black")
    plt.title(f"Horizontal water vapor gradient (seen by DIAL / (z = {z_title} m)\n WV perturbation: {(maxfac-1)*100:.0f} %")
    plt.scatter([0],[0], marker="X", color="red", linewidth=3,\
        label="MWR position")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.tight_layout()
    plt.savefig(outpath+f"/colorplot_wv_DIAL_{(maxfac-1)*100:.0f}.png",
                dpi=150)
    return 0

##############################################################################

def create_360deg_TBplot(tbs, azis=azimuths, elevation="10",\
        maxfac=0.85, minfac=1.15):
    base_path = args.plots
    outpath = ensure_folder_exists(base_path,"sim_hor_grad")

    theta = np.deg2rad(azis)  # Grad → Radiant
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(12, 12))
    for i in range(7):
        ax.plot(theta[::1], tbs[::-1, i], markersize=6, label=f"ch {i}")
    ax.grid(True)
    ax.set_theta_zero_location('N')
    ax.set_title(f'Azimuth Scan (K-band channels {(maxfac-1)*100}_percent pert)\nBrightness temperatures for {elevation}° elevation')
    ax.set_xticks(theta[::5])           
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath+f"/Azimuth_TBs_{elevation}_{(maxfac-1)*100}_percent.png",
                dpi=150)

    return 0

###############################################################################

def get_TB_opposite_differences(tbs, azis=azimuths, azi_pairs=azi_pairs):
    pair_labels = []
    dtbs = np.full((len(azis), 14), np.nan)
    
    for i, azi in enumerate(azis):
        if i>35:
            anti_index = int(i- (len(azis)/2))
        else:
            anti_index = int(i+ (len(azis)/2))
            
        # Calculate TB difference for every pair of opposite angles:
        dtbs[i, :] =  tbs[i,:]- tbs[anti_index,:]        
        pair_labels.append(str(azimuths[i])+" - "+str(azimuths[anti_index]))
    
    return azis, dtbs, pair_labels

##############################################################################

def create_TBdiff_plot(azis, dtbs, pair_labels, elevation=10, maxfac=1.5):
    base_path = args.plots
    outpath = ensure_folder_exists(base_path, "sim_hor_grad")

    theta = np.deg2rad(azis)

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(12, 12))

    for i in range(7):
        ax.plot(theta, dtbs[:, i], marker='o', markersize=6, label=f"Ch {i+1}")

    ax.grid(True)
    ax.set_theta_zero_location('N')
    ax.set_title(f'Mean TB difference between opposite pointing azimuth angles\n'
                 f'Elevation: {elevation}°, WV gradient: {(maxfac-1)*100:.0f}%\n(K-Band)')
    ax.set_xticks(theta[::5])
    ax.set_xticklabels(pair_labels[::5])
    ax.legend(loc="lower right", fontsize=10)

    plt.tight_layout()
    plt.savefig(outpath + f"/Azimuth_TBdiffs_{elevation}_{(maxfac-1)*100:.0f}_percent.png",
                dpi=150)
    plt.close()
    return 0

##############################################################################
# 5th Main code:
##############################################################################

if __name__=="__main__":    
    args = parse_arguments()


    minfacs = np.array([1-0.05,1-0.1,1-0.15,1-0.2, 1-0.25])
    maxfacs = np.array([1+0.05,1+0.1,1+0.15,1+0.2, 1+0.25])
    for minfac, maxfac in zip(minfacs, maxfacs):
        for elevation in elevations:
            tbs = calc_TBs4_72_cases(args, elevation=elevation, maxfac=maxfac, minfac=minfac)
            create_360deg_TBplot(tbs, elevation=elevation, maxfac=maxfac,\
                minfac=minfac)

            azis, dtbs, pair_labels = get_TB_opposite_differences(tbs,\
                azis=azimuths, azi_pairs=azi_pairs)

            create_TBdiff_plot(azis, dtbs, pair_labels, elevation=elevation,\
                maxfac=maxfac)

        # Plot some stuff!!!
        plot_wv_function(args, maxfac=maxfac, minfac=minfac)
        plot_wv_colorplot(args, maxfac=maxfac, minfac=minfac)




















