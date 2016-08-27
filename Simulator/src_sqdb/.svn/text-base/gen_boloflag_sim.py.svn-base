import numpy as np
import pylab as py
import healpy as h
import sqlite3 as sq
import fileinput
import sys

'''
read the simmed FPDB.db
read the experimentally determined FPDB and find the flagged bolo
generate the flag to the simed FPDB.db and generate the new DB
'''


def read_FPDB_fits(filename):
    fpdb = h.mrdfits(filename,hdu=1)
    return fpdb


class readDB():
    def __init__(self):
        self.xml = {}
        self.sq_command = 'select * from pb1_observation'
        self.filename = 'tmp'

    def read_pb1_observation(self):
        self.filename = self.xml['db_ces']
        conn = sq.connect(self.filename)
        c = conn.cursor()
        c.execute(self.sq_command)
        id=[];run_id=[];run_subid=[];dir_ptg=[];first_mjd=[];last_mjd=[]
        for ar in c:
            id.append(int(ar[0]))
            run_id.append(int(ar[1]))
            run_subid.append(int(ar[2]))
            dir_ptg.append(str(ar[3]))
            first_mjd.append(float(ar[4]))
            last_mjd.append(float(ar[5]))
        c.close()
        self.CESdb = {'id':id,'run_id':run_id,'run_subid':run_subid,'dir_ptg':dir_ptg,
                           'first_mjd':first_mjd,'last_mjd':last_mjd}
        return self.CESdb

    def read_BeamParams(self):
#        self.filename = self.xml['file_fpdb_mmin']
        conn = sq.connect(self.filename)
        c = conn.cursor()
#        c.execute('select * from BeamParams')
        c.execute(self.sq_command)
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


    def read_boloid(self):
        self.filename = self.xml['file_boloid']
        conn = sq.connect(self.filename)
        c = conn.cursor()
#        c.execute('select * from pb1_boloid')
        c.execute(self.sq_command)
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
        self.boloid_dict = {'boloid':boloid,'boloname':boloname,'pair':pair,'pixelname':pixelname,'wafer':wafer,
                            'pixel':pixel,'board':board,'squid':squid}
        return self.boloid_dict

    def read_boloid_selected(self):
        self.filename = self.xml['file_flag_pixel']
        print self.filename
        conn = sq.connect(self.filename)
        c = conn.cursor()
        c.execute(self.sq_command)
        boloid=[]; flag=[]
        for ar in c:
            boloid.append(int(ar[0]))
            flag.append(int(ar[1]))
        c.close()
        self.boloid_dict = {'boloid':boloid,'flag':flag}
        return self.boloid_dict

    def display_all(self,db_dict):
        keys = db_dict.keys()
        num = len(db_dict[keys[0]])
        print keys
        for i in range(num):
            tmp = []
            for j in keys: tmp.append(db_dict[j][i])
#            print tmp
        print keys

def gen_boloid_selected(filename,boloid,flag):
    nb = len(boloid)
    conn = sq.connect(filename)
    c = conn.cursor()
    c.execute('create table bolo_flag (boloid integer, flag integer)')
    for i in range(0,nb):
        list_entries = ( int(boloid[i]), int(flag[i]) )
        c.execute('insert into bolo_flag values (?,?)',list_entries)
    conn.commit()
    c.close()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dir = '/project/projectdirs/polar/data/ChileData/calibration/'
#filename = dir+'pm-chile-8.fits'
filename = dir+'beamprm-18obs-hwmap0507-filtp7-rf7-pmChile5.fits'
#filename = dir+'polangle.fits'
fpdb_ab = read_FPDB_fits(filename)

nb = len(fpdb_ab[0])
print nb
#for i in range(nb):
#    print fpdb_ab[1][i], fpdb_ab[2][i], fpdb_ab[3][i], fpdb_ab[4][i]
    
is_nan = np.isnan(fpdb_ab[1])
ind_nan = np.where(is_nan == True)
flag_orig = np.zeros(nb,dtype='int')
flag_orig[ind_nan[0]] = 1
#++++++++++++++++++++++++++++++++++++++++++

dir = '/project/projectdirs/polar/user/chinoney/data/level1/run0707.130/'
filename = 'bolo_names.npy'
bolo_names = np.load(dir+filename)
#print bolo_names
#++++++++++++++++++++++++++++++++++++++++++

dir = '/scratch/scratchdirs/tmatsumu/sim/PB1_NTP/DB/'
filename = dir+'pb1_fpdb_ver0.db'
read_sqDB = readDB()
read_sqDB.sq_command = 'select * from BeamParams;'
read_sqDB.filename = filename
fpdb = read_sqDB.read_BeamParams()
nb_out = len(fpdb['boloid'])

out_boloid = []
out_flag = np.zeros(nb_out,dtype='int')

for i in range(nb_out):
    ind = np.where(fpdb['boloname'][i] == bolo_names)
#    flag_new[i] = flag_orig[ind[0]]
    out_boloid.append(fpdb['boloid'][i])
#    out_flag[i] = flag_orig[ind[0]]
    out_flag[i] = 0 # comment out above and add this line on 2012-9-14
#    print out_boloid[i], out_flag[i]
#++++++++++++++++++++++++++++++++++++++++++

filename_out = '/scratch/scratchdirs/tmatsumu/sim/PB1_NTP/DB/boloflag_sim.db'
gen_boloid_selected(filename_out,out_boloid,out_flag)

