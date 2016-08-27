#!/bin/sh

echo ''
echo '##############################'
echo ''
echo 'Initialization in pyScan, T.Matsumura 2014-10-2'
echo ''
echo ' Execute below then you know what to do '
echo '> init.sh help'
echo ''
echo 'This script creats the symbolic link for pyScans'
echo ' and symbolic link to the pointing data'
echo '##############################'
echo ''

dir_LBSYSver=/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v4.2_release/
dir_data=/group/cmb/litebird/simdata/Scans/dataout/

if [[ $1 = "help" ]];then
    echo '==========================='
    echo ' ###### INSTRUCTION #######'
    echo 'change the following directory names in init.sh'
    echo '   dir_LBSYSver '
    echo '   dir_data (leave this as it is if you want to write the output pointing to the shared pointing directory)'
    echo ''
    echo '   Once you change the directory names, do below'
    echo '> ./init.sh exec'
    echo ''
fi


if [[ $1 = "exec" ]];then
    echo ''
    echo 'add the sympolic link of codes/matsumulib.py to pyScans/src'
    ln -s ${dir}/codes/matsumulib.py  ${dir}/pyScans/src/matsumulib.py
    echo 'add the sympolic link of codes/matsumulib.py to pyScans/src'

    ln -s ${dir_data} ${dir}/pyScans/dataout
    echo ''
    echo 'DONE'
fi