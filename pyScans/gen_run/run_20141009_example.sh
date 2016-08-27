#!/bin/sh

dir_bin=/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v3.2_release/pyScans

# LiteBIRD sidelobe source 

# (option of shell input)
# (title)
# (sampleing [Hz])
# (sec in a day)
# (precession angle [degs])
# (precession period [min])
# (spin axis angle [degs])
# (spin axis period [rpm])
# (sim run period [days])
# (nside for projection)

date=20141009_example

python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 10 86400. 65. 93. 30. 0.1 365 256
