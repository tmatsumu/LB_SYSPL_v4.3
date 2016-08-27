import numpy as np
import pylab as py
import sqlite3 as sq
import lib_PB1SQLDB as libsq
import sys

#--------------------------------------------------------------------------------------
#--- define the basic parameters
dir = '/scratch/scratchdirs/tmatsumu/sim/PB1_NTP/DB/'

boloid = libsq.read_FPDB()
boloid.filename = dir+'pb1_boloid.db'
boloid_arr = boloid.read_boloid()
num_boloid = len(boloid_arr['boloid'])

gen = libsq.lib_GainDB()
gen.sqlite_command = 'select * from pb1_gain_el_nods where run_id==692;'
gain = gen.read_GainDB(dir+'pb1_gain_el_nods.db','el_nods')
num_gain = len(gain['run_id'])

for i in range(num_gain):
    ind = np.where(gain['boloid'][i]==np.array(boloid_arr['boloid']))
    if boloid_arr['pair'][ind[0]] == 't':
        gain['g_begin_el_nods'][i] = 1.
        gain['g_end_el_nods'][i] = 1.

    if boloid_arr['pair'][ind[0]] == 'b':
        gain['g_begin_el_nods'][i] = 1.
        gain['g_end_el_nods'][i] = 1.

gen.create_simGainDB(dir+'pb1_gain_sim_t1_b1.db',gain,'el_nods','sim')

sys.exit()
