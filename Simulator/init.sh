#!/bin/sh

dir_simulator=/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v4.2/Simulator/
dir_proj=/group/cmb/litebird/usr/tmatsumu/LB_SYSPL_v4.2/
dir_db=${dir_proj}DB

echo ''
echo '##############################'
echo 'TYPE '
echo '$ ./init.sh help'
echo ''
echo 'currently set directories are'
echo '  dir_simulator='$dir_simulator
echo '  dir_proj='$dir_proj
echo ''
echo 'change the first two lines of this script for your own use'
echo '  dir_simulator should be a path to LB_SYSPL_v4.2_release/Simulator/ ' 
echo '  dir_proj should be a path to /group/cmb/litebird/usr/your_name/LB_SYSPLv4.2_release/ '
echo '##############################'
echo ''

if [[ $1 = "help" ]];then
    echo ''
    echo '========================================='
    echo '  *** INSTRUCTION *** '
    echo '             by T. Matsumura 2012-10-02'
    echo '   revised by TM on 2013-9-24 specifically for kekcc usage'    
    echo '   revised by TM on 2013-12-20 specifically for kekcc usage'    
    echo '   revised by TM on 2013-12-28 specifically for kekcc usage'    
    echo '   revised by TM on 2014-10-2 specifically for kekcc usage' 
    echo ''
    echo 'The following instruction is generally ONE TIME ONLY action right after checking out LB_SYSPL_v*.* from the svn/git repository'
    echo ''
    echo '1. change the directory name at the first two lines of this file'
    echo ''
    echo '2. make a necessary directories by '
    echo '  > ./init.sh mkdir'
    echo ''
    echo '3. copy the shared database to your own working database directory'
    echo '  > ./init.sh cp_db'
    echo ''
    echo '4. make global_par.py file in gen_run by doing '
    echo '  > ./init.sh gen_global_par'
    echo ''
    echo '5. in src/ directory'
    echo '   > python setup.py build'
    echo '   > cp build/lib.linux-x86_64-2.7/_lib_mapmaker.so ./'
    echo '   if _lib_mapmaker.so file already exits, do not overwrite.'
    echo ''
fi

if [[ $1 = "mkdir" ]];then
    rm data_proj
    rm src_fpdb/fp_db
#    rm src_obsdb/obs_db
    rm src_sqdb/lib_mapmaker.py

    mkdir -p ${dir_proj}/SimedMaps/RunLog
    mkdir -p ${dir_proj}/DB
    mkdir -p ${dir_proj}/DB/fp_db
    mkdir -p ${dir_proj}/DB/gain_db
    mkdir -p ${dir_proj}/DB/obs_db
    mkdir -p ${dir_proj}/mapin
    mkdir -p ${dir_proj}/clin

    ln -s ${dir_proj}/DB/fp_db ${dir_simulator}/src_fpdb/fp_db 
#    ln -s ${dir_proj}DB/obs_db ${dir_simulator}/src_obsdb/obs_db 
    ln -s ${dir_proj} data_proj
    ln -s ${dir_simulator}/src/lib_mapmaker.py ${dir_simulator}/src_sqdb/lib_mapmaker.py

    filename=${dir_simulator}/gen_run/global_par.py
    echo 'import sys' > $filename
    echo ''
    echo 'dir_simulator="'${dir_simulator}'"' >> $filename
    echo 'dir_proj="'${dir_proj}'"' >> $filename
    echo 'dir_db="'${dir_proj}/DB'"' >> $filename
    echo 'if sys.argv[1] == "dir_simulator": print dir_simulator' >> $filename
    echo 'if sys.argv[1] == "dir_proj": print dir_proj' >> $filename
    echo 'if sys.argv[1] == "dir_db": print dir_db' >> $filename
    echo ""
    echo "DONE mkdir"
    echo ""
fi

if [[ $1 = "cp_db" ]];then
    dir_share=/home/cmb/tmatsumu/group_lb/simdata
#    cp $dir_share/DB/fp_db/LB_HFW_example_v2.db ${dir_proj}DB/fp_db
#    cp $dir_share/DB/gain_db/relgain_flat_bias_0.db ${dir_proj}DB/gain_db
    cp $dir_share/DB/fp_db/* ${dir_proj}/DB/fp_db
    cp $dir_share/DB/gain_db/* ${dir_proj}/DB/gain_db
    cp $dir_share/DB/obs_db/LB_L2_20131212_samplerate_SCANSPEC_1440min_65.0degs_90.0min_30.0degs_0.1rpm_365day_nside256_10Hz.db ${dir_proj}/DB/obs_db
    cp /group/cmb/litebird/simdata/clin/reference/standard_lensedtotCls.txt ${dir_proj}/clin
    echo ""
    echo "DONE cp_db"
    echo ""
fi

if [[ $1 = "gen_global_par" ]];then
    echo 'import sys' > $dir_simulator/gen_run/global_par.py
    echo 'dir_simulator="'$dir_simulator'"' >> $dir_simulator/gen_run/global_par.py
    echo 'dir_proj="'$dir_proj'"' >> $dir_simulator/gen_run/global_par.py
    echo 'dir_db="'$dir_db'"' >> $dir_simulator/gen_run/global_par.py
    echo 'if sys.argv[1] == "dir_simulator": print dir_simulator' >> $dir_simulator/gen_run/global_par.py
    echo 'if sys.argv[1] == "dir_proj": print dir_proj' >> $dir_simulator/gen_run/global_par.py
    echo 'if sys.argv[1] == "dir_db": print dir_db' >> $dir_simulator/gen_run/global_par.py
    echo ""
    echo "DONE gen_global_par"
    echo ""
fi