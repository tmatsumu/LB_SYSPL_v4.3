import numpy as np
import pylab as py
import sqlite3 as sq
import sys
#import lib_PB1SQLDB as libsq
import matsumulib as mylib

fname_db1 = sys.argv[1]
fname_db2 = sys.argv[2]
fname_db3 = sys.argv[3]
fname_png = sys.argv[4]

def read_BeamParams(filename):
    conn = sq.connect(filename)
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
    BeamParams = {'boloid':boloid,'boloname':boloname,'xpos':xpos,'ypos':ypos,
                  'polang':polang,'poleff':poleff,'sigma_x':sigma_x,'sigma_y':sigma_y,
                  'amp':amp,'beam_tilt':beam_tilt}
    return BeamParams

def read_boloid(filename):
    print filename
    conn = sq.connect(filename)
    c = conn.cursor()
    c.execute('select * from pb1_boloid')
    boloid=[]; boloname=[]; pair=[]; pixelname=[]; wafer=[]; pixel=[]; board=[]; squid=[];
    for ar in c:
        boloid.append(int(ar[0]))
        boloname.append(str(ar[1]))
        pair.append(str(ar[2]))
        pixelname.append(str(ar[3]))
        wafer.append(str(ar[4]))
        pixel.append(int(ar[5]))
        board.append(int(ar[6]))
        squid.append(str(ar[7]))
    c.close()
    boloid_dict = {'boloid':boloid,'boloname':boloname,'pair':pair,'pixelname':pixelname,'wafer':wafer,'pixel':pixel,'board':board,'squid':squid}
    return boloid_dict


boloid = read_boloid(fname_db1)
beamparams_0 = read_BeamParams(fname_db2)
beamparams = read_BeamParams(fname_db3)

#ind = np.where(np.array(boloid['wafer'])=='9.2')
#ind2 = np.where(beamparams['boloid'] == boloid['boloid'][ind[0]])

#diff_ang =  np.array(beamparams['polang'])[ind[0]] -  np.array(beamparams_0['polang'])[ind[0]]
diff_ang =  np.array(beamparams['polang']) -  np.array(beamparams_0['polang'])
#print diff_ang, len(diff_ang)
mylib.plot_hist(diff_ang,30,init_auto=True,fit=True)
py.savefig(fname_png)
