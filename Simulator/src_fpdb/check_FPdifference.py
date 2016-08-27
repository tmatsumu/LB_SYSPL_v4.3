import numpy as np
import pylab as py
import sqlite3 as sq
import lib_mapmaker as lmm
import matsumulib as mylib
pi = np.pi
radeg = (180./pi)

figname = 'sigma7degs'

readDB = lmm.readDB()
sq_command = 'select * from detector where detid < 370'
readDB.sq_command = sq_command
readDB.filename = 'tmp/LB_HFW_example_original.db'
fpdb_simin = readDB.read_LBFPDB()

readDB.filename = 'tmp/LB_HFW_example_'+figname+'.db'
fpdb_orig = readDB.read_LBFPDB()

#print np.array(fpdb_simin['polang'])*radeg
#print np.array(fpdb_orig['polang'])*radeg

diff = np.array(fpdb_simin['polang'])*radeg - np.array(fpdb_orig['polang'])*radeg

parout = mylib.plot_hist(diff, 30, fit=True, init_auto=True, xtitle='degrees')
py.savefig('/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v3.2/Simulator/src_fpdb/hist_'+figname+'.png')

