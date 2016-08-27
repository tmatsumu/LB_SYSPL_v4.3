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

pi = np.pi
#run = int(sys.argv[1])
#filedbout = sys.argv[2]

class read_file():
    def __init__(self):
        self.filename = 'tmp.txt'
        
    def read_FPDBtxt(self):
        ch = []
        pix = []
        xpos_deg = []
        ypos_deg = []
        polang_deg = []
        wafer = []
        boloid = []
        
        i = 0
        for line in fileinput.input(self.filename):
            ar = line.split()
            if ((len(ar)>1) & (i>3)):
                ch.append(int(ar[0]))
                pix.append(int(ar[1]))
                xpos_deg.append(float(ar[2]))
                ypos_deg.append(float(ar[3]))
                polang_deg.append(float(ar[4]))
                wafer.append(ar[5])
                boloid.append(ar[6])
            i += 1
        print "[READ FPDB]: End reading "+self.filename
        FPDBlist = {'ch': np.array(ch), 'pix': np.array(pix), 'xpos': np.array(xpos_deg), 'ypos': np.array(ypos_deg),
                    'polang': np.array(polang_deg), 'wafer': np.array(wafer), 'boloid': np.array(boloid)}
        num = len(FPDBlist['ch'])
        return FPDBlist, num

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

    def make_FPDB_LB(self):
        if os.path.exists(self.db_name):
            print 'DB already exits'
            return -1
        else:
            conn = sq.connect(self.db_name)
            c = conn.cursor()
            c.execute('create table BeamParams (boloid integer, boloname text, xpos real, ypos real, polang real, poleff real, sigma_x real, sigma_y real, amp real, beam_tilt real)')
            for i in range(0,self.num): 
                list_entries = ( int(self.beam_params['ch'][i]), 
                                 str(self.beam_params['boloid'][i]), 
                                 float(self.beam_params['xpos'][i])/180.*pi,
                                 float(self.beam_params['ypos'][i])/180.*pi, 
                                 float(self.beam_params['polang'][i])/180.*pi, 
                                 float(1.), 
                                 float(30./60./180.*pi),
                                 float(30./60./180.*pi),
                                 float(1.),
                                 float(0.))
                c.execute('insert into BeamParams values (?,?,?,?,?,?,?,?,?,?)',list_entries)
            conn.commit()
            c.close()

class lib_GainDB():
    def _init_(self):
        self.sqlite_command = '.schema'

    def _help_(self, text):
        if text == 'read_GainDB':
            print 'select * from GainParams;'
        if text == 'gain_type':
            print 'el_nods'
            print 'planet'
            print 'planet_w_stimulator'
            print 'planet_w_stimulator_on_2013-03-23'

    def read_GainDB(self,filename):
        conn = sq.connect(filename)
        c = conn.cursor()
#        c.execute('select * from pb1_gain_el_nods;')
        print filename
        c.execute(self.sqlite_command)
#        run_id=[]; run_subid=[]; boloid=[]; 
#        g_begin=[]; g_end=[]; 
#        available=[]; pixel_available=[] 
#        for ar in c:
#            run_id.append(int(ar[0]))
#            run_subid.append(int(ar[1]))
#            boloid.append(int(ar[2]))
#            g_begin.append(float(ar[3]))
#            g_end.append(float(ar[4]))
#            available.append(int(ar[5]))
#            pixel_available.append(int(ar[6]))
#        c.close()
#        gain = {'run_id':run_id,'run_subid':run_subid,'boloid':boloid, \
#                    'g_begin_'+gain_type:g_begin, 'g_end_'+gain_type:g_end, \
#                    'available_'+gain_type: available, 'pixel_available_'+gain_type:pixel_available}
#        return gain
        detid=[];  run_type=[]; params1=[]; params2=[]; params3=[]
        for ar in c:
            detid.append(int(ar[0]))
            run_type.append(ar[1])
            params1.append(float(ar[2]))
            params2.append(float(ar[3]))
            params3.append(float(ar[4]))
        c.close()
        gain = {'detid':detid,'gain_type':run_type,'params1':params1,'params2':params2, 'params3':params3,\
                    'help': 'detid, gain_type, params1, params2, params3' }
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


def io_example_boloid(filename):
    read = read_DB()
    read.filename = filename
    db_new = read.read_boloid()
    read.display_all(db_new)

def io_example_fpdb(filename):
    read = read_DB()
    read.filename = filename
    db_new = read.read_BeamParams()
    read.display_all(db_new)


def gen_LB_FPDB(filename):
    rf = read_file()
    rf.filename = filename
    fpdblist, nb = rf.read_FPDBtxt()

    cf = construct_fpdb()
    cf.db_name = 'tmp.db'
    cf.beam_params = fpdblist
    cf.num = nb
    cf.make_FPDB_LB()


def gen_MuellerDB_tomo(filename):
    out = read_file()
    out.filename = 'simfullfpdb_skview.txt'
    fpdb, nb = out.read_simedFPDB() 

    read = read_DB()
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
        params['T2'][i] = 0.9
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

#if __name__=="__main__":
##############################################################################################################
#filename = '/global/homes/t/tmatsumu/develop/PBI/repo_test/PB1_NTP_develop/src_sqdb/LBFPver1_simfullfpdb_skview.txt'
#gen_LB_FPDB(filename)
#    gen_FPDB_tomo()
#    filename = '/scratch/scratchdirs/tmatsumu/sim/PB1_NTP/DB/pb1_muellermatrix_test.db'
#    gen_MuellerDB_tomo(filename)
#    filename = '/raid/users/tmatsumu/work/PBI/PBI_Chile/Database/database/pb1_boloid.db'
#    io_example_boloid(filename)
#
#    filename = '/raid/users/tmatsumu/work/PBI/PBI_Chile/Database/database/pb1_fpdb_ver0.db'
#    io_example_fpdb(filename)
