#!/bin/sh

#dir=/home/cmb/tmatsumu/develop/LiteBIRD/projects/20130924_LBBasicMM/Mapmake/src_fpdb
dir=/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v3.2/Simulator/src_fpdb

#filename=${dir}/fp_db/LB_HFW_example
filename=${dir}/LB_HFW_example_sigma15degs
skyview=True
rm ${filename}.db
python simedFPDB4LB_v4.py $filename $skyview

