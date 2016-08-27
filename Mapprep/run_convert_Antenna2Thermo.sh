#!/bin/sh

nu_GHz=60
dir_in=/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM/processed_maps/${nu_GHz}GHz/
bsub -q e python convert_Antenna2Thermo.py $dir_in $nu_GHz

nu_GHz=78
dir_in=/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM/processed_maps/${nu_GHz}GHz/
bsub -q e python convert_Antenna2Thermo.py $dir_in $nu_GHz

nu_GHz=100
dir_in=/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM/processed_maps/${nu_GHz}GHz/
bsub -q e python convert_Antenna2Thermo.py $dir_in $nu_GHz

nu_GHz=140
dir_in=/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM/processed_maps/${nu_GHz}GHz/
bsub -q e python convert_Antenna2Thermo.py $dir_in $nu_GHz

nu_GHz=195
dir_in=/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM/processed_maps/${nu_GHz}GHz/
bsub -q e python convert_Antenna2Thermo.py $dir_in $nu_GHz

nu_GHz=280
dir_in=/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM/processed_maps/${nu_GHz}GHz/
bsub -q e python convert_Antenna2Thermo.py $dir_in $nu_GHz 