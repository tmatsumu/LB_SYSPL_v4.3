#!/bin/sh

dir_data=/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v3.2/Simulator/data_proj/SimedMaps/RunLog/
#xml_filename=
dir_out=/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v3.2/Simulator/debug/tmp/
beam_FWHM_arcmin=32.
option=plot_anafast
nside=256

dir_in=test

fname_out=20140717_140GHz_TQU
dir_in0=${dir_data}/${fname_out}/coadd_map/coadd_map1/

fname_out=20140717_140GHz_TQU_sigma1degs
dir_in1=${dir_data}/${fname_out}/coadd_map/coadd_map1/

fname_out=20140717_140GHz_TQU_sigma2degs
dir_in2=${dir_data}/${fname_out}/coadd_map/coadd_map1/

fname_out=20140717_140GHz_TQU_sigma5degs
dir_in5=${dir_data}/${fname_out}/coadd_map/coadd_map1/

fname_out=20140717_140GHz_TQU_sigma7degs
dir_in7=${dir_data}/${fname_out}/coadd_map/coadd_map1/

fname_out=20140717_140GHz_TQU_sigma10degs
dir_in10=${dir_data}/${fname_out}/coadd_map/coadd_map1/

fname_out=20140717_140GHz_TQU_sigma15degs
dir_in15=${dir_data}/${fname_out}/coadd_map/coadd_map1/

fname_out=20140717_140GHz_TQU_combined
python 20140723_main_plotAnafastCls_sum.py $dir_in $dir_in0 $dir_in1 $dir_in2 $dir_in5 $dir_in7 $dir_in10 $dir_in15 $dir_out $beam_FWHM_arcmin $nside $option $fname_out
