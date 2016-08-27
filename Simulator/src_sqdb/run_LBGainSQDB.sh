#!/bin/sh

dir_mm=/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v3.0/Simulator/src/
dir_db=/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v3.0/Simulator/data_proj/DB/

#-----
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0.db bias 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p05.db bias 0.05 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p01.db bias 0.01 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p005.db bias 0.005 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p001.db bias 0.001 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p0005.db bias 0.0005 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p0001.db bias 0.0001 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p00005.db bias 0.00005 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p00001.db bias 0.00001 0. 

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomr_0p05.db random_r 0.05 0. 0. 0. 0.
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomr_0p01.db random_r 0.01 0. 0. 0. 0.
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomr_0p005.db random_r 0.005 0. 0. 0. 0.
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomr_0p001.db random_r 0.001 0. 0. 0. 0.
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomr_0p0005.db random_r 0.0005 0. 0. 0. 0.

#-----
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p05_0p0017_1p5.db 1of_c 0.05 0.0017 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p01_0p0017_1p5.db 1of_c 0.01 0.0017 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p005_0p0017_1p5.db 1of_c 0.005 0.0017 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p001_0p0017_1p5.db 1of_c 0.001 0.0017 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p0005_0p0017_1p5.db 1of_c 0.0005 0.0017 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p0001_0p0017_1p5.db 1of_c 0.0001 0.0017 1.5 0. 0. 0. 

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p05_0p0017_1p5.db 1of_r 0.05 0.0017 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p01_0p0017_1p5.db 1of_r 0.01 0.0017 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p005_0p0017_1p5.db 1of_r 0.005 0.0017 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p001_0p0017_1p5.db 1of_r 0.001 0.0017 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p0005_0p0017_1p5.db 1of_r 0.0005 0.0017 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p0001_0p0017_1p5.db 1of_r 0.0001 0.0017 1.5 0. 0. 0. 

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p05_0p00018_1p5.db 1of_c 0.05 0.00018 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p01_0p00018_1p5.db 1of_c 0.01 0.00018 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p005_0p00018_1p5.db 1of_c 0.005 0.00018 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p001_0p00018_1p5.db 1of_c 0.001 0.00018 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p0005_0p00018_1p5.db 1of_c 0.0005 0.00018 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p0001_0p00018_1p5.db 1of_c 0.0001 0.00018 1.5 0. 0. 0. 

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p05_0p00018_1p5.db 1of_r 0.05 0.00018 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p01_0p00018_1p5.db 1of_r 0.01 0.00018 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p005_0p00018_1p5.db 1of_r 0.005 0.00018 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p001_0p00018_1p5.db 1of_r 0.001 0.00018 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p0005_0p00018_1p5.db 1of_r 0.0005 0.00018 1.5 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p0001_0p00018_1p5.db 1of_r 0.0001 0.00018 1.5 0. 0. 0. 

#-----
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p05_0p0017_1.db 1of_c 0.05 0.0017 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p01_0p0017_1.db 1of_c 0.01 0.0017 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p005_0p0017_1.db 1of_c 0.005 0.0017 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p001_0p0017_1.db 1of_c 0.001 0.0017 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p0005_0p0017_1.db 1of_c 0.0005 0.0017 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p0001_0p0017_1.db 1of_c 0.0001 0.0017 1. 0. 0. 0. 

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p05_0p0017_1.db 1of_r 0.05 0.0017 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p01_0p0017_1.db 1of_r 0.01 0.0017 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p005_0p0017_1.db 1of_r 0.005 0.0017 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p001_0p0017_1.db 1of_r 0.001 0.0017 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p0005_0p0017_1.db 1of_r 0.0005 0.0017 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p0001_0p0017_1.db 1of_r 0.0001 0.0017 1. 0. 0. 0. 

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p05_0p00018_1.db 1of_c 0.05 0.00018 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p01_0p00018_1.db 1of_c 0.01 0.00018 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p005_0p00018_1.db 1of_c 0.005 0.00018 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p001_0p00018_1.db 1of_c 0.001 0.00018 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p0005_0p00018_1.db 1of_c 0.0005 0.00018 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p0001_0p00018_1.db 1of_c 0.0001 0.00018 1. 0. 0. 0. 

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p05_0p00018_1.db 1of_r 0.05 0.00018 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p01_0p00018_1.db 1of_r 0.01 0.00018 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p005_0p00018_1.db 1of_r 0.005 0.00018 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p001_0p00018_1.db 1of_r 0.001 0.00018 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p0005_0p00018_1.db 1of_r 0.0005 0.00018 1. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p0001_0p00018_1.db 1of_r 0.0001 0.00018 1. 0. 0. 0. 

#-----
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p05_0p0017_2.db 1of_c 0.05 0.0017 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p01_0p0017_2.db 1of_c 0.01 0.0017 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p005_0p0017_2.db 1of_c 0.005 0.0017 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p001_0p0017_2.db 1of_c 0.001 0.0017 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p0005_0p0017_2.db 1of_c 0.0005 0.0017 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p0001_0p0017_2.db 1of_c 0.0001 0.0017 2. 0. 0. 0. 

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p05_0p0017_2.db 1of_r 0.05 0.0017 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p01_0p0017_2.db 1of_r 0.01 0.0017 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p005_0p0017_2.db 1of_r 0.005 0.0017 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p001_0p0017_2.db 1of_r 0.001 0.0017 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p0005_0p0017_2.db 1of_r 0.0005 0.0017 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p0001_0p0017_2.db 1of_r 0.0001 0.0017 2. 0. 0. 0. 

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p05_0p00018_2.db 1of_c 0.05 0.00018 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p01_0p00018_2.db 1of_c 0.01 0.00018 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p005_0p00018_2.db 1of_c 0.005 0.00018 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p001_0p00018_2.db 1of_c 0.001 0.00018 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p0005_0p00018_2.db 1of_c 0.0005 0.00018 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p0001_0p00018_2.db 1of_c 0.0001 0.00018 2. 0. 0. 0. 

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p05_0p00018_2.db 1of_r 0.05 0.00018 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p01_0p00018_2.db 1of_r 0.01 0.00018 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p005_0p00018_2.db 1of_r 0.005 0.00018 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p001_0p00018_2.db 1of_r 0.001 0.00018 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p0005_0p00018_2.db 1of_r 0.0005 0.00018 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p0001_0p00018_2.db 1of_r 0.0001 0.00018 2. 0. 0. 0. 
exit

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0.db bias 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p05.db bias 0.05 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p01.db bias 0.01 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p005.db bias 0.005 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p001.db bias 0.001 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p0005.db bias 0.0005 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p0001.db bias 0.0001 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p00005.db bias 0.00005 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_bias_0p00001.db bias 0.00001 0. 

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomr_0p05.db random_r 0.05 0. 0. 0. 0.
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomr_0p01.db random_r 0.01 0. 0. 0. 0.
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomr_0p005.db random_r 0.005 0. 0. 0. 0.
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomr_0p001.db random_r 0.001 0. 0. 0. 0.
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomr_0p0005.db random_r 0.0005 0. 0. 0. 0.

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomc_0p05.db random_c 0.05 0. 0. 0. 0.
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomc_0p01.db random_c 0.01 0. 0. 0. 0.
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomc_0p005.db random_c 0.005 0. 0. 0. 0.
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomc_0p001.db random_c 0.001 0. 0. 0. 0.
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_randomc_0p0005.db random_c 0.0005 0. 0. 0. 0.

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p05_0p0002_2.db 1of_c 0.05 0.0002 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p01_0p0002_2.db 1of_c 0.01 0.0002 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p005_0p0002_2.db 1of_c 0.005 0.0002 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p001_0p0002_2.db 1of_c 0.001 0.0002 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofc_0p0005_0p0002_2.db 1of_c 0.0005 0.0002 2. 0. 0. 0. 

python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p05_0p0002_2.db 1of_r 0.05 0.0002 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p01_0p0002_2.db 1of_r 0.01 0.0002 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p005_0p0002_2.db 1of_r 0.005 0.0002 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p001_0p0002_2.db 1of_r 0.001 0.0002 2. 0. 0. 0. 
python $dir_mm/gen_LBGainSQDB.py $dir_db/fp_db/LB_HFW_example.db $dir_db/gain_db/relgain_flat_1ofr_0p0005_0p0002_2.db 1of_r 0.0005 0.0002 2. 0. 0. 0. 

