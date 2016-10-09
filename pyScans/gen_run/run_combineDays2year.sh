#!/bin/sh

dir_pyScans=/home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.3/pyScans/

#sub_dir=LB_L2_20160919_Wallis_LargeLissajous_samplerate_SCANSPEC_1440min_38.0degs_93.0min_57.0degs_0.1rpm_365day_nside256_10Hz/
#sub_dir=LB_L2_20160919_Wallis_SmallLissajous_samplerate_SCANSPEC_1440min_56.0degs_93.0min_39.0degs_0.1rpm_365day_nside256_10Hz/
#sub_dir=LB_L2_20160919_Wallis_MediumLissajous2_samplerate_SCANSPEC_1440min_50.0degs_93.0min_50.0degs_0.1rpm_365day_nside256_10Hz/
#sub_dir=LB_L2_20160919_LB2_samplerate_SCANSPEC_1440min_65.0degs_93.0min_30.0degs_0.1rpm_365day_nside256_10Hz/
sub_dir=LB_L2_20160919_Wallis_MediumLissajous_samplerate_SCANSPEC_1440min_50.0degs_93.0min_45.0degs_0.1rpm_365day_nside256_10Hz/

filename_db=$dir_pyScans/dataout/$sub_dir/ptg/LBPTG_latlon_eclip

#filename_out=$dir_pyScans/dataout/$sub_dir/ptg/sample_boresight_1day
#python $dir_pyScans/src/main_combined_day2all.py $filename_db 'select * from LBSimPtg where id < 2' $filename_out

filename_out=$dir_pyScans/dataout/$sub_dir/ptg/sample_boresight_365days
python $dir_pyScans/src/main_combined_day2all.py $filename_db 'select * from LBSimPtg where id < 366' $filename_out &
