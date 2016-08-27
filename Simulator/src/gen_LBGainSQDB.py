import numpy as np
import sqlite3 as sq
import sys
import os
#import lib_LBSQLDB as lib_sq
import lib_mapmaker as lmm

'''
in v3.0, modify gen_LBGainSQDB
   - gain_type is now
     bias: bias
     random_c: random in time and coherent over detectors
     random_r: random in time and random over detector
     1of_c: 1of in time and coherent over detectors (one seed for all the detectors)
     1of_r: 1of in time and random over detector (different seed for all the detectors)
'''


filename_fpdb = sys.argv[1]
db_name = sys.argv[2]
gain_type = sys.argv[3]
val1 = float(sys.argv[4])
val2 = float(sys.argv[5])

if (('1of_r' in gain_type) or ('1of_c' in gain_type)):
    val3 = float(sys.argv[6])
    val4 = float(sys.argv[7])
    val5 = float(sys.argv[8])
    val6 = float(sys.argv[9])

'''
python gen_LBGainSQDBpy filename_fpdb gain_type va1 va2
   for gain_type = bias,
       val1=fractional gain bias for t
       val2=fractional gain bias for b
   for gain_type = random_r or random_c,
       val1=fractional gain bias for t
       val2=fractional gain bias for b
   for gain_type = 1of_r, 1of_c
       val1=fractional gain error for t
       val2=knee for t
       val3=alpha for t
       val4=fractional gain error for b
       val5=knee for b
       val6=alpha for b
'''

readDB = lmm.readDB()
readDB.filename = filename_fpdb
readDB.sq_command = 'select * from detector;'
fpdb_simin = readDB.read_LBFPDB()
num_det = len(fpdb_simin['detid'])

params1 = np.zeros(num_det)
params2 = np.zeros(num_det)
params3 = np.zeros(num_det)

print ''
print '[gen_LBGainSQDB.py] Reading ', filename_fpdb
print '[gen_LBGainSQDB.py] output dbname ', db_name
print '[gen_LBGainSQDB.py] gain_type ', gain_type

if 'bias' in gain_type:
    ind_t = np.where(np.array(fpdb_simin['pair'][:]) == 't')
    params1[ind_t[0]] = val1
    ind_b = np.where(np.array(fpdb_simin['pair'][:]) == 'b')
    params1[ind_b[0]] = val2

if (('random_c' in gain_type) or ('random_c' in gain_type)):
    ind_t = np.where(np.array(fpdb_simin['pair'][:]) == 't')
    params1[ind_t[0]] = val1
    ind_b = np.where(np.array(fpdb_simin['pair'][:]) == 'b')
    params1[ind_b[0]] = val2

if (('1of_r' in gain_type) or ('1of_c' in gain_type)):
    ind_t = np.where(np.array(fpdb_simin['pair'][:]) == 't')
    params1[ind_t[0]] = val1
    params2[ind_t[0]] = val2
    params3[ind_t[0]] = val3
    ind_b = np.where(np.array(fpdb_simin['pair'][:]) == 'b')
    params1[ind_b[0]] = val4
    params2[ind_b[0]] = val5
    params3[ind_b[0]] = val6

if os.path.exists(db_name):
    print 'DB already exits'
else:
    conn = sq.connect(db_name)
    c = conn.cursor()
    c.execute('create table GainParams (detid integer, gain_type text, params1 real, params2 real, params3 real )')
    for i in range(0,num_det):
        list_entries = ( int(fpdb_simin['detid'][i]),
                         str(gain_type),
                         float(params1[i]),
                         float(params2[i]),
                         float(params3[i]))
        c.execute('insert into GainParams values (?,?,?,?,?)',list_entries)
    conn.commit()
    c.close()
