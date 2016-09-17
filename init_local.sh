#!/bin/bash

echo ''
echo '=============================================='
echo '  initialization, T. Matsumura 2016-9-2 '
echo '> init_local.sh help'
echo '=============================================='
echo ''

# change this to your version name
LB_SYSPL_version=LB_SYSPL_v4.2

# change the following directory names
#   - dir_mm should be your working directory
#   - dir_proj should be your working data directory, 
#        the easiest change is to simply change the account name from tmatsumu to your account name
#   - dir_ptg should be left as it is. The LiteBIRD pointing is read from the shared simulated data directory
#   - dir_mapin should be left as it is with the same reason as dir_ptg.

dir_default=/Users/tomotake_matsumura/Documents/Projects/LiteBIRD/20160827_LBSYSPLv4.2/
dir_simulator=$dir_default/${LB_SYSPL_version}
dir_proj=$dir_defualt/../data/dataout/simulator/
dir_ptg=$dir_defualt/../data/dataout/ptg/
dir_mapin=$dir_defualt/../data/datain/


if [[ $1 = "help" ]];then
    echo ''
    echo '======================='
    echo ' ***  INSTRUCTION ***'
    echo ''
    echo '> emacs init_local.sh or vi init_local.sh'
    echo '  change the following directory names in the init.sh '
    echo '    LB_SYSPL_version'
    echo '    dir_simulator'
    echo '    (dir_proj) see Note1'
    echo '    (dir_mapin) see Note2'
    echo ''
    echo '     Note1: change dir_ptg only if you want to create a symbolic link to your own pointing dir. If you leave the path as it is, it will create the symbolic link to the shared pointing dir.'
    echo '     Note2: change dir_mapin only if you want to create a symbolic link to your own map dir. If you leave the path as it is, it will create the symbolic link to the shared map dir.'
    echo ' '
    echo '  then do below'
    echo '> init.sh mkdir'
    echo ''
    echo '  currently the set directories in this init.sh are'
    echo '    dir_simulator: '+$dir_simulator
    echo '    dir_proj: '+$dir_proj
    echo '    dir_ptg: '+$dir_ptg
    echo '    dir_mapin: '+$dir_mapin
    echo ''
    echo '  Check to see if all the dir paths are changed.'
    echo ''
    echo '  Once all the action items are done , please do '
    echo '> cd Simulator'
    echo '> ./init.sh'
    echo ''
fi 

if [[ $1 = "mkdir" ]];then
    rm data_proj
    mkdir -p $dir_proj
    mkdir -p $dir_ptg
    mkdir -p $dir_mapin
    ln -s $dir_ptg data_ptg
    ln -s $dir_mapin data_mapin
    ln -s $dir_proj data_proj

fi