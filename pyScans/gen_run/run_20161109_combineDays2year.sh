#!/bin/sh

echo 
echo '========================='
echo 'BEGINNING OF SCRIPT'
echo ''

dir_pyScans=/home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.3/pyScans/

#sub_dir=LB_L2_20160919_Wallis_LargeLissajous_samplerate_SCANSPEC_1440min_38.0degs_93.0min_57.0degs_0.1rpm_365day_nside256_10Hz/
#sub_dir=LB_L2_20160919_Wallis_SmallLissajous_samplerate_SCANSPEC_1440min_56.0degs_93.0min_39.0degs_0.1rpm_365day_nside256_10Hz/
#sub_dir=LB_L2_20160919_Wallis_MediumLissajous2_samplerate_SCANSPEC_1440min_50.0degs_93.0min_50.0degs_0.1rpm_365day_nside256_10Hz/
#sub_dir=LB_L2_20160919_LB2_samplerate_SCANSPEC_1440min_65.0degs_93.0min_30.0degs_0.1rpm_365day_nside256_10Hz/
#sub_dir=LB_L2_20160919_Wallis_MediumLissajous_samplerate_SCANSPEC_1440min_50.0degs_93.0min_45.0degs_0.1rpm_365day_nside256_10Hz/

#sub_dir=LB_L2_20161108_example_samplerate_SCANSPEC_1440min_65.0degs_93.0min_30.0degs_0.3rpm_365day_nside256_10Hz
#sub_dir=LB_L2_20161108_example2_samplerate_SCANSPEC_1440min_65.0degs_93.0min_30.0degs_0.1rpm_365day_nside256_20Hz
#sub_dir=LB_L2_20161108_Wallis_LargeLissajous_30arcmin_samplerate_SCANSPEC_1440min_38.0degs_900.0min_57.0degs_0.361rpm_365day_nside256_15Hz
#sub_dir=LB_L2_20161108_WMAP_30arcmin_samplerate_SCANSPEC_1440min_22.5degs_60.0min_70.0degs_0.465rpm_365day_nside256_21Hz
#sub_dir=LB_L2_20161108_Wallis_Medium_Lissajous_30arcmin_samplerate_SCANSPEC_1440min_50.0degs_2400.0min_45.0degs_0.423rpm_365day_nside256_15Hz
#sub_dir=LB_L2_20161108_Wallis_Small_Lissajous_30arcmin_samplerate_SCANSPEC_1440min_56.0degs_7800.0min_39.0degs_0.48rpm_365day_nside256_15Hz
#sub_dir=LB_L2_20161108_Wallis_Medium_Lissajous2_30arcmin_samplerate_SCANSPEC_1440min_50.0degs_2400.0min_50.0degs_0.423rpm_365day_nside256_15Hz
#sub_dir=LB_L2_20161108_EPIC_30arcmin_samplerate_SCANSPEC_1440min_45.0degs_180.0min_50.0degs_1.0rpm_365day_nside256_37Hz
#sub_dir=LB_L2_20161108_Planck_30arcmin_samplerate_SCANSPEC_1440min_7.5degs_259200.0min_85.0degs_1.0rpm_365day_nside256_48Hz
##sub_dir=LB_L2_20161108_Wallis_SmallLissajous_samplerate_SCANSPEC_1440min_65.0degs_1440.0min_30.0degs_0.1rpm_365day_nside256_10Hz

filename_db=$dir_pyScans/dataout/$sub_dir/ptg/LBPTG_latlon_eclip

#filename_out=$dir_pyScans/dataout/$sub_dir/ptg/sample_boresight_1day
#python $dir_pyScans/src/main_combined_day2all.py $filename_db 'select * from LBSimPtg where id < 2' $filename_out

filename_out=$dir_pyScans/dataout/$sub_dir/ptg/sample_boresight_365days
python $dir_pyScans/src/main_combined_day2all.py $filename_db 'select * from LBSimPtg where id < 366' $filename_out &

ls ${filename_out}*

echo 
echo 'END OF SCRIPT'
echo '========================='
echo ''

exit