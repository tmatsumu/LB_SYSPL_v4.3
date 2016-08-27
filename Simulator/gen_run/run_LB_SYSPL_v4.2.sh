#!/bin/sh

date=$2

dir_ntp=`python global_par.py dir_simulator`
dir_run=${dir_ntp}/gen_run
dir_src=${dir_ntp}/src
dir_xml=${dir_ntp}/xml
dir_db=`python global_par.py dir_db`

sqlite_command_obsdb=$3
sqlite_command_fpdb=$4

echo $i
xml_file=xml_par_${date}.xml

if [[ $1 = "init" ]];then
    echo 'init'
    echo python $dir_run/gen_mm_xml.py par_${date} $dir_xml/${xml_file}
    python $dir_run/gen_mm_xml.py par_${date} $dir_xml/${xml_file} &
    echo 'DONE'
fi

if [[ $1 = "mm" ]];then
    echo '[MAPMAKER]'
    python $dir_src/main_mapmaker_dist.py $dir_xml/${xml_file}  "${sqlite_command_obsdb}" "${sqlite_command_fpdb}" NULL &
    python $dir_src/clean_log.py $dir_xml/${xml_file} mm
    echo 'DONE'
fi

if [[ $1 == "extTODmm" ]];then
    echo '[external TOD MAPMAKER] ', $1
    echo '   TOD dir name:', $5
    python $dir_src/main_mapmaker_dist.py $dir_xml/${xml_file}  "${sqlite_command_obsdb}" "${sqlite_command_fpdb}" extTODmm $5 &
    echo 'DONE'
fi

if [[ $1 = "coadd" ]];then
    echo '[COADD]'
    python $dir_src/main_coaddDB.py $dir_xml/${xml_file} "${sqlite_command_obsdb}" 
    python $dir_src/main_coaddmaps.py $dir_xml/${xml_file} "${sqlite_command_obsdb}" &
    python $dir_src/clean_log.py $dir_xml/${xml_file} coadd
    echo 'DONE'
fi

if [[ $1 = "extTODmm_coadd" ]];then
    echo '[COADD]'
    python $dir_src/main_coaddDB.py $dir_xml/${xml_file} "${sqlite_command_obsdb}" 
    python $dir_src/main_coaddmaps_ext.py $dir_xml/${xml_file} "${sqlite_command_obsdb}" $5 &
    python $dir_src/clean_log.py $dir_xml/${xml_file} coadd
    echo 'DONE'
fi

if [[ $1 = "clean_npy" ]];then
    echo '[CLEAN_NPY]'
    python $dir_src/clean_log.py $dir_xml/${xml_file} clean_npy
    echo 'DONE'
fi

if [[ $1 = "clean_npz" ]];then
    echo '[CLEAN_NPZ]'
    python $dir_src/clean_log.py $dir_xml/${xml_file} clean_npz
    echo 'DONE'
fi

if [[ $1 = "clean_tod" ]];then
    echo '[CLEAN_TOD]'
    python $dir_src/clean_log.py $dir_xml/${xml_file} clean_tod
    echo 'DONE'
fi

if [[ $1 = "anafast" ]];then
    echo '[ANAFAST]'
    python $dir_src/main_anafast.py $dir_xml/${xml_file} "${sqlite_command_obsdb}" anafast $5 &
    echo 'DONE'
fi

if [[ $1 = "plot_anafast" ]];then
    echo '[PLOT_ANAFAST]'
    python $dir_src/main_plotAnafastCls.py $dir_xml/${xml_file} "${sqlite_command_obsdb}" plot_anafast $5 $6 &
    echo 'DONE'
fi

if [[ $1 = "extTODmm_anafast" ]];then
    echo '[ANAFAST]'
    python $dir_src/main_anafast.py $dir_xml/${xml_file} "${sqlite_command_obsdb}" extTODmm_anafast $5 &
    echo 'DONE'
fi

if [[ $1 = "extTODmm_plot_anafast" ]];then
    echo '[PLOT_ANAFAST]'
    python $dir_src/main_plotAnafastCls.py $dir_xml/${xml_file} "${sqlite_command_obsdb}" extTODmm_plot_anafast $5 $6 &
    echo 'DONE'
fi

if [[ $1 = "gen_mappng" ]];then
    echo '[GEN_MAPPNG]'
    python $dir_src/main_gen_mappng.py $dir_xml/${xml_file} "${sqlite_command_obsdb}" gen_mappng $5 $6 &
    echo 'DONE'
fi

if [[ $1 = "cp_png" ]];then
    echo '[CP_PNG]'
    python $dir_src/main_gen_mappng.py $dir_xml/${xml_file} "${sqlite_command_obsdb}" cp_png &
    echo 'DONE'
fi

if [[ $1 = "cal_map_stats" ]];then    
    python $dir_src/cal_map_stats.py $dir_xml/${xml_file} "${sqlite_command_obsdb}" &
    echo 'DONE'
fi

if [[ $1 = "xpure_window" ]];then
    echo '[XPURE_WINDOW]'
    python $dir_src/main_xpure.py $dir_xml/${xml_file} xpure_window "${sqlite_command_obsdb}" &
    python $dir_src/clean_log.py $dir_xml/${xml_file} xpure_window
    echo 'DONE'
fi

if [[ $1 = "plot_xpure_window" ]];then    
    echo '[PLOT_XPURE_WINDOW]'
    python $dir_src/main_xpure.py $dir_xml/${xml_file} plot_xpure_window "${sqlite_command_obsdb}" &
    echo 'DONE'
fi    


if [[ $1 = "xpure_cl" ]];then
    echo '[XPURE_Cl]'
    python $dir_src/main_xpure.py $dir_xml/${xml_file} xpure_cl "${sqlite_command_obsdb}" &
    echo 'DONE'
fi

if [[ $1 = "plot_xpure_cl" ]];then
    echo '[PLOT_XPURE_Cl]'
    python $dir_src/main_xpure.py $dir_xml/${xml_file} plot_xpure_cl "${sqlite_command_obsdb}" &
    python $dir_src/clean_log.py $dir_xml/${xml_file}  xpure_cl
    echo 'DONE'
fi

if [[ $1 = "mm_db" ]];then
    python $dir_run/make_piperundb.py par_${date} $dir_db/sim_run.db
    echo 'DONE'
fi

if [[ $1 = "help" ]];then
    echo ' init'
    echo ' mm'
    echo ' coadd'
    echo ''
    echo ' extTODmm'
    echo ' extTODmm_coadd your_toddir_name'
    echo ''
    echo ' anafast'
    echo ' plot_anafast'
    echo ''
    echo ' extTODmm_anafast'
    echo ' extTODmm_plot_anafast'
    echo ''
    echo ' gen_mappng'
    echo ' cal_map_stats'
#    echo ' xpure_window'
#    echo ' plot_xpure_window'
#    echo ' xpure_cl'
#    echo ' plot_xpure_cl'
    echo ' mm_db'
fi
