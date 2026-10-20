#! /bin/bash

# Tophat Winter:
./python_src/preproc_resample2scan_freq.py -i "/home/aki/PhD_data/tophat_joyce_2025/2025/*/sups_joy_mwr00_l1_tb_p00_*.nc" -o "/home/aki/PhD_data/scans/MWR_scans_JOYCE_Tophat_202510_12.nc" ;
# Joyhat Vital I
./python_src/preproc_resample2scan_freq.py -i "/home/aki/PhD_data/scans/joyhat_raw_jun_jul_aug_sep/MWR_1C01_*.nc" -o "/home/aki/PhD_data/scans/MWR_scans_JOYCE_Joyhat_202406_09.nc" ;

###
# Vital II
# 1st Aachen
./python_src/preproc_resample2scan_freq.py -i "/home/aki/PhD_data/scans/vitII_site_eval/aachen_may26/MWR_1C01_aachen_*.nc" -o "/home/aki/PhD_data/scans/MWR_scans_aachen_may26.nc" ;
# 2nd Sinthern:
./python_src/preproc_resample2scan_freq.py -i "/home/aki/PhD_data/scans/vitII_site_eval/sinthern_may26/MWR_1C01_sinthern_*.nc" -o "/home/aki/PhD_data/scans/MWR_scans_sinthern_may26.nc" ;
# 3rd: Vettweiss:
./python_src/preproc_resample2scan_freq.py -i "/home/aki/PhD_data/scans/vitII_site_eval/vettweiss_may26/MWR_1C01_vettweiss_*.nc" -o "/home/aki/PhD_data/scans/MWR_scans_vettweiss_may26.nc" ;
# 4th: Joyce:
./python_src/preproc_resample2scan_freq.py -i "/home/aki/PhD_data/scans/vitII_site_eval/juelich_may26/sups_joy_mwr00_l1_tb_p00_*.nc" -o "/home/aki/PhD_data/scans/MWR_scans_juelich_may26.nc" ;
# 5th: Airport:
./python_src/preproc_resample2scan_freq.py -i "/home/aki/PhD_data/scans/vitII_site_eval/airport_may26/*/MWR_1C01_airport_*.nc" -o "/home/aki/PhD_data/scans/MWR_scans_airport_may26.nc" ;
# 6th Cologne / Biozentrum:
./python_src/preproc_resample2scan_freq.py -i "/home/aki/PhD_data/scans/vitII_site_eval/foghat_may26/MWR_1C01_*.nc" -o "/home/aki/PhD_data/scans/MWR_scans_foghat_may26.nc" ;
# 7th: Meckenheim:
./python_src/preproc_resample2scan_freq.py -i "/home/aki/PhD_data/scans/vitII_site_eval/mecken_jun26/MWR_1C01_*.nc" -o "/home/aki/PhD_data/scans/MWR_scans_mechat_jun26.nc" ;

########################
./python_src/find_obstacles.py
