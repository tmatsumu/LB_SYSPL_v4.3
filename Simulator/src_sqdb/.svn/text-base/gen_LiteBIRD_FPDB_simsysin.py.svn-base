import numpy as np
import sqlite3 as sq
import fileinput
import sys
import os

run = int(sys.argv[1])
filedbout = sys.argv[2]


def gen_FPDB_tomo():
    out = read_file()
#    out.filename = '/raid/users/tmatsumu/work/PBI/PBI_Chile/Database/data/simfullfpdb_skview.txt'
    out.filename = 'LBFPver1_simfullfpdb_skview.txt'
    fpdb, num = out.read_simedFPDB()

    read = read_DB()
#    read.filename = '/raid/users/tmatsumu/work/PBI/PBI_Chile/Database/database/pb1_boloid.db'
    read.filename = 'pb1_boloid.db'
    db_ref = read.read_boloid()

    boloname = []
    boloid = []
    wafer = []
    for i in range(num):
        ind = np.where( fpdb['boloname'][i] == np.array(db_ref['boloname']) )
#        ind = np.where( (fpdb['boloname'][i] == np.array(db_ref['boloname'])) & (np.array(db_ref['wafer']) == '9.4') )
        if len(ind[[0][0]]) !=0: 
            boloid.append( db_ref['boloid'][ind[0][0]] )
            wafer.append( db_ref['wafer'][ind[0][0]])
            boloname.append( db_ref['boloname'][ind[0][0]])
#            print boloid[i], wafer[i]
    fpdb['boloid'] = boloid
    fpdb['wafer'] = wafer
    fpdb['boloname'] = boloname
#    sys.exit()

    option_createDB = True
    if option_createDB:
        gen = construct_fpdb()
#        gen.db_name = '/raid/users/tmatsumu/work/PBI/PBI_Chile/Database/database/tmp1.db'
        gen.db_name = filedbout
        gen.beam_params = fpdb
        gen.num = num
        beamFWHM = 30./60. # [deg]
#        gen.make_FPDB()
        if run==1: gen.make_FPDB_Syspolang(1.,'shift')
        if run==2: gen.make_FPDB_Syspolang(0.5,'shift')
        if run==3: gen.make_FPDB_Syspolang(0.25,'shift')

        if run==4: gen.make_FPDB_Syspolang(10.,'random')
        if run==5: gen.make_FPDB_Syspolang(5.,'random')
        if run==6: gen.make_FPDB_Syspolang(4.,'random')
        if run==7: gen.make_FPDB_Syspolang(3.,'random')
        if run==8: gen.make_FPDB_Syspolang(2.,'random')
        if run==9: gen.make_FPDB_Syspolang(1.,'random')
        if run==10: gen.make_FPDB_Syspolang(.5,'random')

        if run==11: gen.make_FPDB_Sysptg(beamFWHM*2.0,beamFWHM*2.0,'diffptg')
        if run==12: gen.make_FPDB_Sysptg(beamFWHM*1.0,beamFWHM*1.0,'diffptg')
        if run==13: gen.make_FPDB_Sysptg(beamFWHM*0.5,beamFWHM*0.5,'diffptg')
        if run==14: gen.make_FPDB_Sysptg(beamFWHM*0.3,beamFWHM*0.3,'diffptg')
        if run==15: gen.make_FPDB_Sysptg(beamFWHM*0.1,beamFWHM*0.1,'diffptg')
        if run==16: gen.make_FPDB_Sysptg(beamFWHM*0.05,beamFWHM*0.05,'diffptg')
        if run==17: gen.make_FPDB_Sysptg(beamFWHM*0.01,beamFWHM*0.01,'diffptg')

        if run==18: gen.make_FPDB_Sysptg(beamFWHM*2.0,beamFWHM*2.0,'pixel_random')
        if run==19: gen.make_FPDB_Sysptg(beamFWHM*1.0,beamFWHM*1.0,'pixel_random')
        if run==20: gen.make_FPDB_Sysptg(beamFWHM*0.5,beamFWHM*0.5,'pixel_random')
        if run==21: gen.make_FPDB_Sysptg(beamFWHM*0.3,beamFWHM*0.3,'pixel_random')
        if run==22: gen.make_FPDB_Sysptg(beamFWHM*0.1,beamFWHM*0.1,'pixel_random')
        if run==23: gen.make_FPDB_Sysptg(beamFWHM*0.05,beamFWHM*0.05,'pixel_random')
        if run==24: gen.make_FPDB_Sysptg(beamFWHM*0.01,beamFWHM*0.01,'pixel_random')

        if run==25: gen.make_FPDB_Sysptg(beamFWHM*2.0,beamFWHM*2.0,'boresight_shift')
        if run==26: gen.make_FPDB_Sysptg(beamFWHM*1.0,beamFWHM*1.0,'boresight_shift')
        if run==27: gen.make_FPDB_Sysptg(beamFWHM*0.5,beamFWHM*0.5,'boresight_shift')
        if run==28: gen.make_FPDB_Sysptg(beamFWHM*0.3,beamFWHM*0.3,'boresight_shift')
        if run==29: gen.make_FPDB_Sysptg(beamFWHM*0.1,beamFWHM*0.1,'boresight_shift')
        if run==30: gen.make_FPDB_Sysptg(beamFWHM*0.05,beamFWHM*0.05,'boresight_shift')
        if run==31: gen.make_FPDB_Sysptg(beamFWHM*0.01,beamFWHM*0.01,'boresight_shift')

        if run==32: gen.make_FPDB_Syspolang_part(5.,'random','10.2')
        if run==33: gen.make_FPDB_Sysptg_part(beamFWHM*0.05,beamFWHM*0.05,'diffptg','10.2')
        if run==34: gen.make_FPDB_Sysptg_part(beamFWHM*0.05,beamFWHM*0.05,'boresight_shift','10.2')
        if run==35: gen.make_FPDB_Sysptg_part(beamFWHM*0.0,beamFWHM*0.0,'boresight_shift','10.2')

        if run==36: gen.make_FPDB_Syspolang(8.,'random')
        if run==37: gen.make_FPDB_Syspolang(9.,'random')

    
    read = read_DB()
#    read.filename = '/raid/users/tmatsumu/work/PBI/PBI_Chile/Database/database/tmp1.db'
    read.filename = filedbout
    db_new = read.read_BeamParams()
    read.display_all(db_new)

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

##############################################################################################################
    gen_FPDB_tomo()
#    filename = '/scratch/scratchdirs/tmatsumu/sim/PB1_NTP/DB/pb1_muellermatrix_test.db'
#    gen_MuellerDB_tomo(filename)
#    filename = '/raid/users/tmatsumu/work/PBI/PBI_Chile/Database/database/pb1_boloid.db'
#    io_example_boloid(filename)
#
#    filename = '/raid/users/tmatsumu/work/PBI/PBI_Chile/Database/database/pb1_fpdb_ver0.db'
#    io_example_fpdb(filename)
