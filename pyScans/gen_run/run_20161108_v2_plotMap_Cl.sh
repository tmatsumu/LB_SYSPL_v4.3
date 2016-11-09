#!/bin/sh

dir_in=/home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.3/Simulator/data_proj/SimedMaps/RunLog/gainval_v001/coadd_map/coadd_map1/
dir_src=/home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.3/pyScans/src/
dir_data=$dir_in/
filename=ClHout_anafast
#filename=LB_L2_20160919_Wallis_MediumLissajous_samplerate_SCANSPEC_1440min_50.0degs_93.0min_45.0degs_0.1rpm_365day_nside256_10Hz
#filename=LB_L2_20160919_Wallis_LargeLissajous_samplerate_SCANSPEC_1440min_38.0degs_93.0min_57.0degs_0.1rpm_365day_nside256_10Hz
#filename=LB_L2_20160919_Wallis_MediumLissajous2_samplerate_SCANSPEC_1440min_50.0degs_93.0min_50.0degs_0.1rpm_365day_nside256_10Hz
#filename=LB_L2_20160919_Wallis_SmallLissajous_samplerate_SCANSPEC_1440min_56.0degs_93.0min_39.0degs_0.1rpm_365day_nside256_10Hz
#filename=LB_L2_20161031_Wallis_SmallLissajous_samplerate_SCANSPEC_1440min_45.0degs_93.0min_50.0degs_0.1rpm_365day_nside256_10Hz

max=1000
min=0

python $dir_src/main_plotMollviewMap_from_fits.py $dir_data/ mapH $max $min
python $dir_src/main_plotAnafastCls_1rawfits.py $dir_data/ $filename 0

display  $dir_data/$filename/ps/${filename}_nhits.png &
display  $dir_data/${filename}.png &

