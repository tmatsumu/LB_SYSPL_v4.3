import numpy as np
import lib_PB1SQLDB as lib_sq

dir_db = '/scratch/scratchdirs/tmatsumu/sim/PB1_NTP/DB/'
filename = dir_db + 'beamprm_polangle_gradation_x1.db'
out = lib_sq.read_FPDB()
out.filename = filename
BeamParams = out.read_BeamParams()

num = len(BeamParams['boloid'])

flag = np.zeros(num,dtype='int')
for i in range(0,num):
    if ((BeamParams['boloname'][i][0] == 'D') or (BeamParams['boloname'][i][-1:] == 'd')): 
        flag[i] = 1
        print BeamParams['boloname'][i], flag[i]
    else:
        flag[i] = 0
#        print BeamParams['boloname'][i], flag[i]

lib_sq.make_boloflag("test.db",BeamParams['boloid'],flag)
