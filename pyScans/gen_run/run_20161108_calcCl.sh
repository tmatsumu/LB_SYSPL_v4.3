#!/bin/sh

dir=/home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.3/pyScans/
dir_src=$dir/src/

dir_fits=/home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.3/pyScans/dataout/
sub_dir=LB_L2_20160919_LB2_samplerate_SCANSPEC_1440min_65.0degs_93.0min_30.0degs_0.1rpm_365day_nside256_10Hz
fits_filename=LB_L2_20160919_LB2_samplerate_SCANSPEC_1440min_65.0degs_93.0min_30.0degs_0.1rpm_365day_nside256_10Hz_nhits

python $dir_src/main_CalAnafast_from_fits.py $dir_fits/$sub_dir/ps/ $fits_filename 256