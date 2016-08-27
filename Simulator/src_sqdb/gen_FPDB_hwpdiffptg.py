import numpy as np
import pylab as py
import sqlite3 as sq
import lib_PB1SQLDB as libsq
import sys


#--------------------------------------------------------------------------------------
#--- define the basic parameters
dir = '/scratch/scratchdirs/tmatsumu/sim/PB1_NTP/DB/'
i_th = int(sys.argv[1])
factor=int(sys.argv[2])

#--------------------------------------------------------------------------------------
#--- read the simulated pb1_fpdb_ver0.db
read_db = libsq.read_DB()
read_db.filename=dir+'pb1_fpdb_ver0.db'
beam_params = read_db.read_BeamParams_selective([1,1,1,1,1,1,1,1,1,1])
num_bolo = len(beam_params['boloid'])
#print beam_params
print num_bolo

#--------------------------------------------------------------------------------------
#--- read the nishino diff hwp file
dbfile_nishino = ['beamprm_20120530_031419_hwp112.5.db', \
                      'beamprm_20120630_011427_hwp0.0.db', \
                      'beamprm_20120629_012206_hwp157.5.db', \
                      'beamprm_20120623_014144_hwp22.5.db', \
                      'beamprm_20120729_234559_hwp135.0.db', \
                      'beamprm_20120726_235949_hwp67.5.db', \
                      'beamprm_20120726_001053_hwp45.0.db', \
                      'beamprm_20120827_192720_hwp90.0.db' ]
i_dbfile_nishino = dbfile_nishino[i_th]
read_db.filename=dir+i_dbfile_nishino
beam_params_diffhwp = read_db.read_BeamParams_selective([1,1,1,1,0,0,0,0,0,0])
num_nishino = len(beam_params_diffhwp['boloid'])
#print beam_params_diffhwp
print num_nishino

#--------------------------------------------------------------------------------------
#--- shuffle the two DB
newflag = np.ones(num_nishino,dtype='int')
flag = []
xpos_ave = []; ypos_ave = []
xpos_t = []; ypos_t = []
xpos_b = []; ypos_b = []
xpos = []; ypos = []

newboloid = []; newpolang=[]; newboloname=[]

for i in range(num_bolo):
    pixelname = beam_params['boloname'][i]
    ind_t = np.where(np.array(beam_params_diffhwp['boloname']) == pixelname[:-1]+'t')
    ind_b = np.where(np.array(beam_params_diffhwp['boloname']) == pixelname[:-1]+'b')
    ind = np.where(np.array(beam_params_diffhwp['boloname']) == pixelname)

    if len(ind[0])>1: 
        print 'there are more than two identical boloname', len(ind[0])
        break

    if ( (beam_params_diffhwp['xpos'][ind_t[0]] == None) 
        or (beam_params_diffhwp['ypos'][ind_t[0]] == None) 
        or (beam_params_diffhwp['xpos'][ind_b[0]] == None) 
        or (beam_params_diffhwp['ypos'][ind_b[0]] == None) ): 
        xpos_ave.append(0.)
        ypos_ave.append(0.)
        xpos.append(0.)
        ypos.append(0.)
        newflag[beam_params_diffhwp['boloid'][ind_t[0]]] = 1
        newflag[beam_params_diffhwp['boloid'][ind_b[0]]] = 1
        newboloid.append(int(beam_params_diffhwp['boloid'][ind[0]]))
        newpolang.append(float(beam_params['polang'][i]))
        newboloname.append(pixelname)
        flag.append(1)
    else: 
        xpos_ave.append((float(beam_params_diffhwp['xpos'][ind_t[0]])+float(beam_params_diffhwp['xpos'][ind_b[0]]))*0.5)
        ypos_ave.append((float(beam_params_diffhwp['ypos'][ind_t[0]])+float(beam_params_diffhwp['ypos'][ind_b[0]]))*0.5)
#        xpos.append(float(beam_params_diffhwp['xpos'][ind[0]]))
#        ypos.append(float(beam_params_diffhwp['ypos'][ind[0]]))
        xpos.append( (float(beam_params_diffhwp['xpos'][ind[0]]))*factor  
                     + float(xpos_ave[i])*(1.-factor) )
        ypos.append( (float(beam_params_diffhwp['ypos'][ind[0]]))*factor
                     + float(ypos_ave[i])*(1.-factor) )
        newflag[beam_params_diffhwp['boloid'][ind[0]]] = 0
        newboloid.append(int(beam_params_diffhwp['boloid'][ind[0]]))
        newpolang.append(float(beam_params['polang'][i]))
        newboloname.append(pixelname)
        flag.append(0)

#################################################################################################3
#out_beamparams = {'boloid': np.array(beam_params['boloid']), 
out_beamparams = {'boloid': np.array(newboloid), 
#                  'boloname': np.array(beam_params['boloname']),
                  'boloname': np.array(newboloname),
                  'xpos': np.array(xpos_ave), 
                  'ypos': np.array(ypos_ave),
                  'polang': np.array(newpolang), 
                  'poleff': np.array(beam_params['poleff']), 
                  'sigma_x': np.array(beam_params['sigma_x']), 
                  'sigma_y': np.array(beam_params['sigma_y']), 
                  'amp': np.array(beam_params['amp']) ,
                  'beam_tilt': np.array(beam_params['beam_tilt']) }

con_fpdb = libsq.construct_fpdb()  
con_fpdb.db_name = dir+i_dbfile_nishino+'_simin_ave_fact'+str(factor)+'.db'
con_fpdb.beam_params = out_beamparams
con_fpdb.num = num_bolo
con_fpdb.make_FPDB()

#################################################################################################3
#out_beamparams = {'boloid': np.array(beam_params['boloid']), 
out_beamparams = {'boloid': np.array(newboloid), 
#                  'boloname': np.array(beam_params['boloname']),
                  'boloname': np.array(newboloname),
                  'xpos': np.array(xpos),
                  'ypos': np.array(ypos),
#                  'polang': np.array(beam_params['polang']), 
                  'polang': np.array(newpolang), 
                  'poleff': np.array(beam_params['poleff']), 
                  'sigma_x': np.array(beam_params['sigma_x']), 
                  'sigma_y': np.array(beam_params['sigma_y']), 
                  'amp': np.array(beam_params['amp']) ,
                  'beam_tilt': np.array(beam_params['beam_tilt']) }

con_fpdb = libsq.construct_fpdb()  
con_fpdb.db_name = dir+i_dbfile_nishino+'_simin_diffptg_fact'+str(factor)+'.db'
con_fpdb.beam_params = out_beamparams
con_fpdb.num = num_bolo
con_fpdb.make_FPDB()

#################################################################################################3
libsq.make_boloflag(dir+i_dbfile_nishino+'_simin_flag_fact'+str(factor)+'.db', np.array(beam_params_diffhwp['boloid']), newflag)

#################################################################################################3
tmp1 = np.where(np.array(beam_params['boloname']) == '10.1_16b')
tmp2 = np.where(np.array(beam_params_diffhwp['boloname']) == '10.1_16b')
print tmp1, tmp2
print beam_params['boloid'][tmp1[0]], beam_params_diffhwp['boloid'][tmp2[0]], beam_params['boloname'][tmp1[0]], beam_params_diffhwp['boloname'][tmp2[0]], xpos[tmp1[0]], flag[tmp1[0]], newflag[tmp2[0]]  
