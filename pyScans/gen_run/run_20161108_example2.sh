#!/bin/sh

dir_bin=/home/cmb/tmatsumu/LB_SYSPL/LB_SYSPL_v4.3/pyScans/

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

date=20161108_example2

python ${dir_bin}/src/run_scan_todgen_c.mod.py y LB_L2_${date}_samplerate 20 86400. 65. 93. 30. 0.1 365 256
