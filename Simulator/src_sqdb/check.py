import numpy as np
import sqlite3 as sq
import sys
import lib_PB1SQLDB as libsq

class read_DB():
    def __init__(self):
        self.filename = 'tmp.db'
        
    def read_BeamParams(self):
        conn = sq.connect(self.filename)
        c = conn.cursor()
        c.execute('select * from BeamParams')
        boloid=[]; boloname=[]; xpos=[]; ypos=[]; polang=[]; poleff=[]
        sigma_x=[]; sigma_y=[]; amp=[]; beam_tilt=[]
        for ar in c:
            boloid.append(int(ar[0]))
            boloname.append(str(ar[1]))
            xpos.append(float(ar[2]))
            ypos.append(float(ar[3]))
            polang.append(float(ar[4]))
            poleff.append(float(ar[5]))
            sigma_x.append(float(ar[6]))
            sigma_y.append(float(ar[7]))
            amp.append(float(ar[8]))
            beam_tilt.append(float(ar[9]))
        c.close()
        self.BeamParams = {'boloid':boloid,'boloname':boloname,'xpos':xpos,'ypos':ypos,
                           'polang':polang,'poleff':poleff,'sigma_x':sigma_x,'sigma_y':sigma_y,
                           'amp':amp,'beam_tilt':beam_tilt}
        return self.BeamParams



read_db = libsq.read_DB()
read_db.filename = '/scratch/scratchdirs/tmatsumu/sim/PB1_NTP/DB/beamprm_20120530_031419_hwp112.5.db'
beam1 = read_db.read_BeamParams_selective([1,1,0,0,0,0,0,0,0,0])

read_db = libsq.read_DB()
read_db.filename = '/scratch/scratchdirs/tmatsumu/sim/PB1_NTP/DB/pb1_fpdb_ver0.db'
beam2 = read_db.read_BeamParams_selective([1,1,0,0,0,0,0,0,0,0])

num = len(beam2['boloid'])

for i in range(num):
    ind = np.where(beam2['boloid'][i] == np.array(beam1['boloid']))
#    print ind[0]
    if ( beam1['boloname'][ind[0]] != beam2['boloname'][i] ):
        print beam2['boloid'][i], beam1['boloname'][ind[0]], beam2['boloname'][i]
