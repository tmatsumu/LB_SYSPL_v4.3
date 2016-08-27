#!/bin/sh

waittime=1s

# 
run_name=20160818_140GHz_TQU_7
sleep $waittime | ./run_LB_SYSPL_v4.1.sh $1 $run_name 'select * from LBSimPtg where id > 180 and id < 366' 'select * from detector where detid < 370' $2 $3
#sleep $waittime | ./run_LB_SYSPL_v4.1.sh $1 $run_name 'select * from LBSimPtg where id < 2' 'select * from detector where detid < 10' $2 $3

exit
