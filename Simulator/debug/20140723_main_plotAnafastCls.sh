#!/bin/sh

dir_data=/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v3.2/Simulator/data_proj/SimedMaps/RunLog/
#xml_filename=
dir_out=/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v3.2/Simulator/debug/tmp/
beam_FWHM_arcmin=32.
option=plot_anafast
nside=256

fname_out=20140717_140GHz_TQU_sigma1degs
dir_in=${dir_data}/${fname_out}/coadd_map/coadd_map1/
mean=0.10348617
python 20140723_main_plotAnafastCls.py $dir_in $dir_out $beam_FWHM_arcmin $nside $option $mean $fname_out

fname_out=20140717_140GHz_TQU_sigma2degs
dir_in=${dir_data}/${fname_out}/coadd_map/coadd_map1/
mean=-0.14944824
python 20140723_main_plotAnafastCls.py $dir_in $dir_out $beam_FWHM_arcmin $nside $option $mean $fname_out

fname_out=20140717_140GHz_TQU_sigma5degs
dir_in=${dir_data}/${fname_out}/coadd_map/coadd_map1/
mean=0.93833654
python 20140723_main_plotAnafastCls.py $dir_in $dir_out $beam_FWHM_arcmin $nside $option $mean $fname_out

fname_out=20140717_140GHz_TQU_sigma7degs
dir_in=${dir_data}/${fname_out}/coadd_map/coadd_map1/
mean=-0.14599284
python 20140723_main_plotAnafastCls.py $dir_in $dir_out $beam_FWHM_arcmin $nside $option $mean $fname_out

fname_out=20140717_140GHz_TQU_sigma10degs
dir_in=${dir_data}/${fname_out}/coadd_map/coadd_map1/
mean=-1.23017429
python 20140723_main_plotAnafastCls.py $dir_in $dir_out $beam_FWHM_arcmin $nside $option $mean $fname_out

fname_out=20140717_140GHz_TQU_sigma15degs
dir_in=${dir_data}/${fname_out}/coadd_map/coadd_map1/
mean=-2.2545263
python 20140723_main_plotAnafastCls.py $dir_in $dir_out $beam_FWHM_arcmin $nside $option $mean $fname_out

