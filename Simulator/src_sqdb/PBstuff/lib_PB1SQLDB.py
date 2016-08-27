#!/usr/local/bin/env python

import numpy as np
import sqlite3 as sq
import fileinput
import sys
import os

'''
   lib_PB1SQLDB.py

   This library defines classes and functions to manage the database related I/O in python.

   Written by T. Matsumura
   2012-5-23: initial version is placed on NERSC

   Example)
     case1:
       $ import lib_PB1SQLDB as lib_pb1sqldb
       $ lib_pb1sqldb.io_example_boloid(database_name)
     case2:
       $ import lib_PB1SQLDB as lib_pb1sqldb
       $ lib_pb1sqldb.io_example_fpdb(fpdb_name)
'''

#run = int(sys.argv[1])
#filedbout = sys.argv[2]

class read_fpdbfile():
    def __init__(self):
        self.filename = 'tmp.txt'
        
    def read_FPDB(self):
        ch = []
        az = []
        el = []
        amp = []
        sig_x = []
        sig_y = []
        theta_tilt = []
        
        i = 0
        for line in fileinput.input(self.filename):
            ar = line.split()
            if ((len(ar)>1) & (i>3)):
                ch.append(int(ar[0]))
                az.append(float(ar[1]))
                el.append(float(ar[2]))
                amp.append(float(ar[3]))
                sig_x.append(float(ar[4]))
                sig_y.append(float(ar[5]))
                theta_tilt.append(float(ar[6]))
            i += 1
        print "[READ FPDB]: End reading "+self.filename
        FPDBlist = {'ch': np.array(ch), 'xpos': np.array(az), 'ypos': np.array(el), 'amp': np.array(amp),
                    'sig_x': np.array(sig_x), 'sig_y': np.array(sig_y), 'theta_tilt': np.array(theta_tilt)}
        num = len(FPDBlist['ch'])
        return FPDBlist, num

    def read_Beamparam(self):
        ch = []
        az = []
        el = []
        polang = []
        poleff = []
        boloid = []
        
        i = 0
        for line in fileinput.input(self.filename):
            ar = line.split()
            if (len(ar)>1):
                ch.append(int(ar[0]))
                az.append(float(ar[1]))
                el.append(float(ar[2]))
                polang.append(float(ar[3]))
                poleff.append(float(ar[4]))
                boloid.append(str(ar[5]))
            i += 1
        print "[READ FPDB]: End reading "+self.filename
        FPDBlist = {'id': np.array(ch), 'xpos': np.array(az), 'ypox': np.array(el), 'polang': np.array(polang),
                    'poleff': np.array(poleff), 'boloid': np.array(boloid)}
        num = len(FPDBlist['id'])
        return FPDBlist, num

    def read_simedFPDB(self):
        boloname = []
        pixelname = []
        wafer = []
        pair = []
        pixel = []
        az = []
        el = []
        polang = []
        poleff = []
        sigma_x = []
        sigma_y = []
        beam_theta = []
        
        i = 0
        for line in fileinput.input(self.filename):
            ar = line.split()
            if ((len(ar)>1) and (i>0)):
                boloname.append(ar[0])
                pixelname.append(ar[1])
                wafer.append(ar[2])
                pair.append(ar[3])
                pixel.append(ar[4])
                az.append(ar[5])
                el.append(ar[6])
                polang.append(ar[7])
                poleff.append(ar[8])
                sigma_x.append(ar[9])
                sigma_y.append(ar[10])
                beam_theta.append(ar[11])
            i += 1
        print "[READ FPDB]: End reading "+self.filename
        FPDBlist = {'boloname': np.array(boloname), 'pixelname': np.array(pixelname), 
                    'wafer': np.array(wafer), 'pair': np.array(pair),
                    'pixel': np.array(pixel), 'xpos': np.array(az), 'ypos': np.array(el), 
                    'polang': np.array(polang), 'poleff': np.array(poleff),
                    'sigma_x': np.array(sigma_x), 'sigma_y': np.array(sigma_y), 'amp':np.ones(len(boloname)),  
                    'beam_tilt': np.array(beam_theta)}
        num = len(FPDBlist['boloname'])
        return FPDBlist, num


def make_boloflag(db_name,boloid,flag):
    num = len(boloid)
    if os.path.exists(db_name):
        print 'DB already exits'
        return -1
    else:
        conn = sq.connect(db_name)
        c = conn.cursor()
        c.execute('create table boloselect_flag (boloid integer, flag integer)')
        for i in range(0,num): 
            list_entries = ( int(boloid[i]), int(flag[i]) ) 
            c.execute('insert into boloselect_flag values (?,?)',list_entries)
        conn.commit()
        c.close()        

def read_boloflag(filenanme):
    conn = sq.connect(filename)
    c = conn.cursor()
    c.execute('select * from boloselect_flag')
    boloid=[]; flag=[]                                                                                           
    for ar in c:
        boloid.append(int(ar[0]))
        flag.append(int(ar[1]))
    c.close()
    boloflag = {'boloid':boloid,'flag':flag}
    return boloflag

class construct_fpdb():
    def __init__(self):
        self.db_name = 'tmp.db'
        self.beam_params = {}
        self.num = 1
       
    def make_FPDB(self):
        if os.path.exists(self.db_name):
            print 'DB already exits'
            return -1
        else:
            conn = sq.connect(self.db_name)
            c = conn.cursor()
            c.execute('create table BeamParams (boloid integer, boloname text, xpos real, ypos real, polang real, poleff real, sigma_x real, sigma_y real, amp real, beam_tilt real)')
            for i in range(0,self.num): 
                list_entries = ( int(self.beam_params['boloid'][i]), 
                                 str(self.beam_params['boloname'][i]), 
                                 float(self.beam_params['xpos'][i]),
                                 float(self.beam_params['ypos'][i]), 
                                 float(self.beam_params['polang'][i]), 
                                 float(self.beam_params['poleff'][i]), 
                                 float(self.beam_params['sigma_x'][i]),
                                 float(self.beam_params['sigma_y'][i]),
                                 float(self.beam_params['amp'][i]),
                                 float(self.beam_params['beam_tilt'][i]))
                c.execute('insert into BeamParams values (?,?,?,?,?,?,?,?,?,?)',list_entries)
            conn.commit()
            c.close()

    def make_FPDB_Syspolang(self,amp,type):
#        print self.db_name, amp, type
        if os.path.exists(self.db_name):
            print 'DB already exits'
            return -1
        else:
            conn = sq.connect(self.db_name)
            c = conn.cursor()
            c.execute('create table BeamParams (boloid integer, boloname text, xpos real, ypos real, polang real, poleff real, sigma_x real, sigma_y real, amp real, beam_tilt real)')
            amp_tmp = np.zeros(self.num)
            if 'random' in type: 
                amp_tmp = np.float_(self.beam_params['polang'])+amp*np.random.normal(0.,1.,self.num)
            if 'shift' in type: 
                amp_tmp = np.float_(self.beam_params['polang'])+amp

            for i in range(0,self.num): 
                list_entries = ( int(self.beam_params['boloid'][i]), 
                                 str(self.beam_params['boloname'][i]), 
                                 float(self.beam_params['xpos'][i]),
                                 float(self.beam_params['ypos'][i]), 
                                 float(amp_tmp[i]),
                                 float(self.beam_params['poleff'][i]), 
                                 float(self.beam_params['sigma_x'][i]),
                                 float(self.beam_params['sigma_y'][i]),
                                 float(self.beam_params['amp'][i]),
                                 float(self.beam_params['beam_tilt'][i]))
                c.execute('insert into BeamParams values (?,?,?,?,?,?,?,?,?,?)',list_entries)
            conn.commit()
            c.close()
#            print self.db_name, amp, type

    def make_FPDB_Sysptg(self,amp_x,amp_y,type):
        if os.path.exists(self.db_name):
            print 'DB already exits'
            return -1
        else:
            conn = sq.connect(self.db_name)
            c = conn.cursor()
            c.execute('create table BeamParams (boloid integer, boloname text, xpos real, ypos real, polang real, poleff real, sigma_x real, sigma_y real, amp real, beam_tilt real)')
            if 'diffptg' in type:
                amp = amp_x*np.random.normal(0.,1.,self.num)
                amp_x_tmp = np.float_(self.beam_params['xpos'])*(1.+amp)
                amp_y_tmp = np.float_(self.beam_params['ypos'])*(1.+amp)

            if 'xy_shift' in type:
                pass
            if 'boresight_shift' in type:
                amp_x_tmp = np.float_(self.beam_params['xpos'])+amp_x
                amp_y_tmp = np.float_(self.beam_params['ypos'])+amp_x

            for i in range(0,self.num): 
                list_entries = ( int(self.beam_params['boloid'][i]), 
                                 str(self.beam_params['boloname'][i]), 
                                 amp_x_tmp[i],
                                 amp_y_tmp[i],
                                 float(self.beam_params['polang'][i]), 
                                 float(self.beam_params['poleff'][i]), 
                                 float(self.beam_params['sigma_x'][i]),
                                 float(self.beam_params['sigma_y'][i]),
                                 float(self.beam_params['amp'][i]),
                                 float(self.beam_params['beam_tilt'][i]))
                c.execute('insert into BeamParams values (?,?,?,?,?,?,?,?,?,?)',list_entries)
            conn.commit()
            c.close()

    def make_FPDB_Syspolang_part(self,amp,type,wafer_name):
        num_all = len(self.beam_params['wafer'])
        ind = np.where(wafer_name == np.array(self.beam_params['wafer']))
        ind = ind[0]
        num = len(ind)
#        print num, ind
#        print self.beam_params['wafer']
        print np.array(self.beam_params['wafer'])[ind], num
#        print self.db_name, amp, type
        if os.path.exists(self.db_name):
            print 'DB already exits'
            return -1
        else:
            print self.db_name
            conn = sq.connect(self.db_name)
            c = conn.cursor()
            c.execute('create table BeamParams (boloid integer, boloname text, xpos real, ypos real, polang real, poleff real, sigma_x real, sigma_y real, amp real, beam_tilt real)')
            amp_tmp = np.zeros(num)
            if 'random' in type: 
                rand_tmp = amp*np.random.normal(0.,1.,num)
                np.save('rand_tmp',rand_tmp)
                amp_tmp = np.float_(self.beam_params['polang'])
                amp_tmp[ind] = np.float_(self.beam_params['polang'][ind])+rand_tmp
            if 'shift' in type: 
                amp_tmp = np.float_(self.beam_params['polang'])+amp

            for i in range(0,num): 
                list_entries = ( int(self.beam_params['boloid'][ind[i]]), 
                                 str(self.beam_params['boloname'][ind[i]]), 
                                 float(self.beam_params['xpos'][ind[i]]),
                                 float(self.beam_params['ypos'][ind[i]]), 
                                 float(amp_tmp[ind[i]]),
                                 float(self.beam_params['poleff'][ind[i]]), 
                                 float(self.beam_params['sigma_x'][ind[i]]),
                                 float(self.beam_params['sigma_y'][ind[i]]),
                                 float(self.beam_params['amp'][ind[i]]),
                                 float(self.beam_params['beam_tilt'][ind[i]]))
                c.execute('insert into BeamParams values (?,?,?,?,?,?,?,?,?,?)',list_entries)
            conn.commit()
            c.close()

    def make_FPDB_Sysptg_part(self,amp_x,amp_y,type,wafer_name):
        num_all = len(self.beam_params['wafer'])
        ind = np.where(wafer_name == np.array(self.beam_params['wafer']))
        ind = ind[0]
        num = len(ind)
        print np.array(self.beam_params['wafer'])[ind], num
#        print self.db_name, amp, type
        if os.path.exists(self.db_name):
            print 'DB already exits'
            return -1
        else:
            print self.db_name
            conn = sq.connect(self.db_name)
            c = conn.cursor()
            c.execute('create table BeamParams (boloid integer, boloname text, xpos real, ypos real, polang real, poleff real, sigma_x real, sigma_y real, amp real, beam_tilt real)')

            amp_x_tmp = np.float_(self.beam_params['xpos'])
            amp_y_tmp = np.float_(self.beam_params['ypos'])

            if 'diffptg' in type:
                amp = amp_x*np.random.normal(0.,1.,num)
                amp_x_tmp[ind] = np.float_(self.beam_params['xpos'][ind])*(1.+amp)
                amp_y_tmp[ind] = np.float_(self.beam_params['ypos'][ind])*(1.+amp)

            if 'xy_shift' in type:
                pass
            if 'boresight_shift' in type:
                amp_x_tmp[ind] = np.float_(self.beam_params['xpos'][ind])+amp_x
                amp_y_tmp[ind] = np.float_(self.beam_params['ypos'][ind])+amp_x

            for i in range(0,num): 
                list_entries = ( int(self.beam_params['boloid'][ind[i]]), 
                                 str(self.beam_params['boloname'][ind[i]]), 
                                 amp_x_tmp[ind[i]],
                                 amp_y_tmp[ind[i]],
                                 float(self.beam_params['polang'][ind[i]]), 
                                 float(self.beam_params['poleff'][ind[i]]), 
                                 float(self.beam_params['sigma_x'][ind[i]]),
                                 float(self.beam_params['sigma_y'][ind[i]]),
                                 float(self.beam_params['amp'][ind[i]]),
                                 float(self.beam_params['beam_tilt'][ind[i]]))
                c.execute('insert into BeamParams values (?,?,?,?,?,?,?,?,?,?)',list_entries)
            conn.commit()
            c.close()
#            print self.db_name, amp, type

class read_FPDB():
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


    def read_BeamParams_selective(self,select):
        conn = sq.connect(self.filename)
        c = conn.cursor()
        c.execute('select * from BeamParams')
        boloid=[]; boloname=[]; xpos=[]; ypos=[]; polang=[]; poleff=[]
        sigma_x=[]; sigma_y=[]; amp=[]; beam_tilt=[]
        for ar in c:
            if select[0]==1: boloid.append(int(ar[0]))
            if select[1]==1: boloname.append(str(ar[1]))
            if select[2]==1: 
                if ar[2]==None: xpos.append(None)
                else: xpos.append(float(ar[2]))
            if select[3]==1: 
                if ar[3]==None: ypos.append(None)
                else: ypos.append(float(ar[3]))
            if select[4]==1: 
                if ar[4]==None: polang.append(None)
                else: polang.append(float(ar[4]))
            if select[5]==1: 
                if ar[5]==None: poleff.append(None)
                else: poleff.append(float(ar[5]))
            if select[6]==1: 
                if ar[6]==None: sigma_x.append(None)
                else: sigma_x.append(float(ar[6]))
            if select[7]==1: 
                if ar[7]==None: sigma_y.append(None)
                else: sigma_y.append(float(ar[7]))
            if select[8]==1: 
                if ar[8]==None: amp.append(None)
                else: amp.append(float(ar[8]))
            if select[9]==1: 
                if ar[9]==None: beam_tilt.append(None)
                else: beam_tilt.append(float(ar[9]))
        c.close()
        self.BeamParams = {'boloid':boloid,'boloname':boloname,'xpos':xpos,'ypos':ypos,
                           'polang':polang,'poleff':poleff,'sigma_x':sigma_x,'sigma_y':sigma_y,
                           'amp':amp,'beam_tilt':beam_tilt,'select':select}
        return self.BeamParams

    def read_boloid(self):
        conn = sq.connect(self.filename)
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
        self.boloid_dict = {'boloid':boloid,'boloname':boloname,'pair':pair,'pixelname':pixelname,'wafer':wafer,'pixel':pixel,'board':board,'squid':squid}
        return self.boloid_dict

    def display_all(self,db_dict):
        keys = db_dict.keys()
        num = len(db_dict[keys[0]])
        print keys
        for i in range(num):
            tmp = []
            for j in keys: tmp.append(db_dict[j][i])
            print tmp
        print keys

class construct_muellerdb():
    def __init__(self):
        self.db_name = 'tmp.db'
        self.params = {}
        self.fpdb = {}
        self.num = 1
        
    def make_MuellerDB(self):
        if os.path.exists(self.db_name):
            print 'DB already exits'
            return -1
        else:
            conn = sq.connect(self.db_name)
            c = conn.cursor()
            c.execute('create table MuellerParams (boloid integer, boloname text, T1 real, T2 real, trans real, rho real, cosdel real)')            
            for i in range(0,self.num): 
                list_entries = ( int(self.fpdb['boloid'][i]),
                                 str(self.fpdb['boloname'][i]),
                                 float(self.params['T1'][i]), 
                                 float(self.params['T2'][i]), 
                                 float(self.params['trans'][i]),
                                 float(self.params['rho'][i]), 
                                 float(self.params['cosdel'][i]) )
                c.execute('insert into MuellerParams values (?,?,?,?,?,?,?)',list_entries)
            conn.commit()
            c.close()

def gen_MuellerDB_tomo(filename):
    out = read_fpdbfile()
    out.filename = 'simfullfpdb_skview.txt'
    fpdb, nb = out.read_simedFPDB() 

    read = read_FPDB()
    read.filename = 'pb1_boloid.db'
    db_ref = read.read_boloid()
    boloname = []
    boloid = []
    wafer = []
    for i in range(nb):
        ind = np.where( fpdb['boloname'][i] == np.array(db_ref['boloname']) )
        if len(ind[[0][0]]) !=0 :
            boloid.append( db_ref['boloid'][ind[0][0]] )
            wafer.append( db_ref['wafer'][ind[0][0]])
            boloname.append( db_ref['boloname'][ind[0][0]])

    fpdb['boloid'] = boloid
    fpdb['wafer'] = wafer
    fpdb['boloname'] = boloname 
    
    params = {'T1': np.zeros(nb), \
                  'T2': np.zeros(nb), \
                  'trans': np.zeros(nb), \
                  'rho': np.zeros(nb), \
                  'cosdel': np.zeros(nb) }
    
    for i in range(nb):
        params['T1'][i] = 1.
        params['T2'][i] = 1.
        params['trans'][i] = 1.
        params['rho'][i] = 0.
        params['cosdel'][i] = -1.

    mueller_class = construct_muellerdb()
    mueller_class.db_name = filename
    mueller_class.params = params
    mueller_class.fpdb = fpdb
    mueller_class.num = nb
    mueller_class.make_MuellerDB()
    return 0


class lib_GainDB():
    def _init_(self):
        self.sqlite_command = '.schema'

    def _help_(self, text):
        if text == 'read_GainDB':
            print 'select * from pb1_gain_el_nods;'
        if text == 'gain_type':
            print 'el_nods'
            print 'planet'
            print 'planet_w_stimulator'
            print 'planet_w_stimulator_on_2013-03-23'

    def read_GainDB(self,filename,gain_type):
        conn = sq.connect(filename)
        c = conn.cursor()
#        c.execute('select * from pb1_gain_el_nods;')
        c.execute(self.sqlite_command)
        run_id=[]; run_subid=[]; boloid=[]; 
        g_begin=[]; g_end=[]; 
        available=[]; pixel_available=[] 
        for ar in c:
            run_id.append(int(ar[0]))
            run_subid.append(int(ar[1]))
            boloid.append(int(ar[2]))
            g_begin.append(float(ar[3]))
            g_end.append(float(ar[4]))
            available.append(int(ar[5]))
            pixel_available.append(int(ar[6]))
        c.close()
        gain = {'run_id':run_id,'run_subid':run_subid,'boloid':boloid, \
                    'g_begin_'+gain_type:g_begin, 'g_end_'+gain_type:g_end, \
                    'available_'+gain_type: available, 'pixel_available_'+gain_type:pixel_available}
        return gain
        
    def create_simGainDB(self,filename,gain_in,gain_type_in,gain_type_out):
        num = len(gain_in['run_id'])
        if os.path.exists(filename):
            print 'DB already exits'
            return -1
        else:
            conn = sq.connect(filename)
            c = conn.cursor()
            c.execute('create table pb1_gain_'+gain_type_out+' (run_id integer, run_subid integer, boloid integer, g_begin_'+gain_type_out+' real, g_end_'+gain_type_out+' real, available_'+gain_type_out+' integer, pixel_available_'+gain_type_out+' integer)')
            for i in range(0,num): 
                list_entries = ( int(gain_in['run_id'][i]), 
                                 int(gain_in['run_subid'][i]), 
                                 int(gain_in['boloid'][i]),
                                 float(gain_in['g_begin_'+gain_type_in][i]), 
                                 float(gain_in['g_end_'+gain_type_in][i]),
                                 int(gain_in['available_'+gain_type_in][i]), 
                                 int(gain_in['pixel_available_'+gain_type_in][i]) )
                c.execute('insert into pb1_gain_'+gain_type_out+' values (?,?,?,?,?,?,?)',list_entries)
            conn.commit()
            c.close()
        

def io_example_boloid(filename):
    read = read_FPDB()
    read.filename = filename
    db_new = read.read_boloid()
    read.display_all(db_new)

def io_example_fpdb(filename):
    read = read_FPDB()
    read.filename = filename
    db_new = read.read_BeamParams()
    read.display_all(db_new)

#if __name__=="__main__":
##############################################################################################################
#    gen_FPDB_tomo()
#    filename = '/scratch/scratchdirs/tmatsumu/sim/PB1_NTP/DB/pb1_muellermatrix_test.db'
#    gen_MuellerDB_tomo(filename)
#    filename = '/raid/users/tmatsumu/work/PBI/PBI_Chile/Database/database/pb1_boloid.db'
#    io_example_boloid(filename)
#
#    filename = '/raid/users/tmatsumu/work/PBI/PBI_Chile/Database/database/pb1_fpdb_ver0.db'
#    io_example_fpdb(filename)
