import numpy as np
import pylab as py
import sqlite3 as sq
import lib_mapmaker as lmm
import matsumulib as mylib
import sys

#+++++++++++++++++++++++++++++++++++++++++++++++++++++

pi = np.pi
radeg = (180./pi)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++

filename = sys.argv[1]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++

readDB = lmm.readDB()
sq_command = 'select * from detector where detid < 370'
readDB.sq_command = sq_command
readDB.filename = filename
fpdb_simin = readDB.read_LBFPDB()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++
# self.BeamParams = {'detid':detid,'pixel':pixel,'xpos':xpos,'ypos':ypos,'polang':polang,'wafer':wafer,'detname':detname,'pair':pair}
num = len(fpdb_simin['polang'])

print 'Original values in the DB'
for i in range(num):
    print np.array(fpdb_simin['detid'][i]), np.array(fpdb_simin['pixel'][i]), np.array(fpdb_simin['xpos'][i]), np.array(fpdb_simin['ypos'][i]), np.array(fpdb_simin['polang'][i]), np.array(fpdb_simin['wafer'][i]), np.array(fpdb_simin['detname'][i]), np.array(fpdb_simin['pair'][i])

print ''
print 'Angles are convereted to degree'
for i in range(num):
    print np.array(fpdb_simin['detid'][i]), np.array(fpdb_simin['pixel'][i]), np.array(fpdb_simin['xpos'][i])*radeg, np.array(fpdb_simin['ypos'][i])*radeg, np.array(fpdb_simin['polang'][i])*radeg, np.array(fpdb_simin['wafer'][i]), np.array(fpdb_simin['detname'][i]), np.array(fpdb_simin['pair'][i])


