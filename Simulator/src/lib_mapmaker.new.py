import numpy as np
import healpy as h
import pylab as py
import math as m
from numpy import fft
from numpy import random
import Gen_Mapmaking_inputs as gmi
import ctypes
import fileinput
import sys
import os
import re
import time
import glob
import ReadMapMakeXml as rxml
import string
#import AnalysisBackend.misc.util as util
#from multiprocessing import Pool
#from multiprocessing import Process
import copy
import _lib_mapmaker as lib_c
#from guppy import hpy; hpy=hpy()
import pickle
import sqlite3 as sq
from scipy.weave import inline

pi = np.pi
radeg = (180./pi)

#######################################################################################################
#######################################################################################################
def info(title, runtime_init):
    if runtime_init != 0: 
        print ''
        print "________________________________________"
        print title
        print '     module name:', __name__
        print '  parent process:', os.getppid()
        print '      process id:', os.getpid()
        print '        run time:', time.time()-runtime_init
        print "________________________________________"
        print ''
#######################################################################################################
#######################################################################################################
# write the text file
#   input: filename, array
#   output: output file
def write_txtf(fname,array):
    nb = len(array)
    f = open(fname, "w")
    for i in range(0,nb):
        f.write('%f \n' % (array[i]))
    print '[WRITE TEXT FILE]: ', fname
    f.close()
#######################################################################################################
#######################################################################################################
# write the text file
#   input: filename, array
#   output: output file
def write_txt2f(fname,array1,array2):
    nb = len(array1)
    f = open(fname, "w")
    for i in range(0,nb):
        f.write('%f %f \n' % (array1[i], array2[i]))
    print '[WRITE TEXT FILE]: ', fname
    f.close()                          
#######################################################################################################
#######################################################################################################
# write the text file
#   input: filename, array
#   output: output file
def write_txt2i(fname,array1,array2):
    nb = len(array1)
    f = open(fname, "w")
    for i in range(0,nb):
        f.write('%d %d \n' % (array1[i], array2[i]))
    print '[WRITE TEXT FILE]: ', fname
    f.close()                          
#######################################################################################################
#######################################################################################################
# write the text file
#   input: filename, array
#   output: output file
def write_txt3f(fname,array1,array2,array3):
    nb = len(array1)
    f = open(fname, "w")
    for i in range(0,nb):
        f.write('%f %f %f \n' % (array1[i], array2[i],array3[i]))
    print '[WRITE TEXT FILE]: ', fname
    f.close()                          
#######################################################################################################
#######################################################################################################
# read the DB information
#   input: xml file
#   output: FP data, pixel list, CES list
class DataSelection():
    def __init__(self):
        self.xml = xml
    def read_CESlist(self):
        pass
    def read_Pixellist_sim(self):
        pass
    def read_Pixellist(self):
        pass
    def read_flag_Relgainlist(self):
        pass
    def read_flag_HWPangles(self):
        pass

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

    def read_boloid_flag(self):
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
#        flag_in = np.zeros(len(flag),dtype='int')  ## ADDED TO DISABLE THE FLAG OPTION: ### REMOVE IT ###
#        self.boloid_dict = {'boloid':boloid,'flag':flag_in.tolist()} # added 2012-9-10 above and this lines
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

class shuffle_pixlist():
    def __init__(self):
        self.boloid = {}
        self.boloid_flag = {}

    def gen_pixlist(self): 
        tmp_flag = np.array(self.boloid_flag['flag'])
        ind = np.where(tmp_flag == 0) # exclude if the flag is non-zero.
        num = len(ind[0])
        print '[shuffle_pixlist] the # of bolos', num
        tmp_boloid_selected=[]; tmp_pair_selected=[]; tmp_pixel_selected=[]; tmp_id_selected=[]
        for i in range(0,num):
            tmp_boloid_selected.append(self.boloid['boloid'][ind[0][i]])
            tmp_pair_selected.append(self.boloid['pair'][ind[0][i]])
            tmp_pixel_selected.append(self.boloid['pixel'][ind[0][i]])
            tmp_id_selected.append(int(i))
        # so far, the unflagged bolos are selected.

        tmp_t=[];tmp_b=[];done_pix=[-1];ind_tmp=[];ind_tmp.append(0);idx_pix=[]
        for i in range(0,num):
            if 't' == tmp_pair_selected[i]:
                ind = np.where(tmp_pixel_selected[i] == np.array(tmp_pixel_selected))
#                tmp_t = tmp_boloid_selected[ind[0][0]]
#                tmp_b = tmp_boloid_selected[ind[0][1]]
                tmp_t = tmp_id_selected[ind[0][0]]
                tmp_b = tmp_id_selected[ind[0][1]]
                tmp_p = tmp_pixel_selected[ind[0][0]]
                if 't' in tmp_pair_selected[ind[0][0]]:
                    idx_pix.append([tmp_t,tmp_b,tmp_p])
                if 't' in tmp_pair_selected[ind[0][1]]:
                    idx_pix.append([tmp_b,tmp_t,tmp_p])
        print '[lib_mapmaker.py] the number of bolometers to process', num
        return idx_pix

    def gen_selected_fpdb(self,beampar):
        tmp_flag = np.array(self.boloid_flag['flag'])
        tmp_boloid = np.array(self.boloid_flag['boloid'])
        ind = np.where(tmp_flag == 0)
        ind = np.array(ind[0])
        num = len(ind)
#        print num
        boloid_selected = {}
        keys = ['boloid','boloname','xpos','ypos','polang','poleff','sigma_x','sigma_y','amp','beam_tilt']
        dtypes = ['int','str','float','float','float','float','float','float','float','float']
        num_keys = len(keys)
        for i in range(0,num):
            idx = np.where(beampar['boloid'] == tmp_boloid[ind[i]])
#            print i, idx, len(idx)
            idx = idx[0]
            if i==0:
                for j in range(0,num_keys):
                    boloid_selected[keys[j]] = np.zeros(num,dtype=dtypes[j])
#            print '.', idx, idx[0]
            for j in range(0,num_keys):
                boloid_selected[keys[j]][i] = beampar[keys[j]][idx]     ## added 2012-9-10 idx -> idx[0]
        return boloid_selected
#######################################################################################################
#######################################################################################################
# generate the scanset dir
#
#
def Gen_scansetdirs(xml_filename):
    print ""
    print "[Gen_scansetdirs]: START", xml_filename
    # xml_filename = sys.argv[1]
    xml_input = rxml.Get_Mapmake_Inputs(xml_filename)

    fname_ptg = xml_input["file_input_ptg"]
        
    date_i = xml_input["date_i"]
    date_f = xml_input["date_f"]
    
    fileNames = glob.glob(fname_ptg+"/observation_*")

    nb = len(fileNames)
#    for i in range(0,nb): print fileNames[i]
    
    nb_tmp = len(fname_ptg)
    fileNames_date = np.zeros(nb,int)
    for i in range(0,nb): fileNames_date[i] = int(fileNames[i][(nb_tmp+13):]) 

    ind = np.where((fileNames_date >= date_i) & (fileNames_date <= date_f))

    nb_obs = len(ind[0])
    dir_obs = []
    dir_obs_part = []
    for i in range(0,nb_obs):
        dir_obs.append(fname_ptg+"/observation_"+str(fileNames_date[ind[0][i]]))
        dir_obs_part.append("observation_"+str(fileNames_date[ind[0][i]]))
    
    file_input_pointing = []
    file_input_pointing_part = []
    for i in range(0,nb_obs):
        fileNames = glob.glob(dir_obs[i]+'/scan*')
        nb = len(fileNames)
#        print fileNames
#        print nb, nb_obs
        for j in range(0,nb): 
            file_input_pointing.append(fileNames[j])
            file_input_pointing_part.append(fileNames[j][(nb_tmp+13):])
            
#    print ""
#    print "  scan setdirs..."
    nb_tmp = nb_obs*nb
#    for i in range(0,nb_tmp): print "   "+file_input_pointing[i]

    xml_input["dirs_pointing"] = file_input_pointing
    xml_input["dirs_pointing_part"] = file_input_pointing_part

#    sys.exit()
    print "[Gen_scansetdirs]: END"
    return xml_input
#######################################################################################################
#######################################################################################################
# complex array
#  input: two arrays, real and imaginary
#  output: complex array
def complex_arr(arr1,arr2):
    nb1 = np.size(arr1)
    nb2 = np.size(arr2)
    
    if nb1 == nb2:
        out_arr = np.zeros(nb1,complex)
        out_arr.real=arr1
        out_arr.imag=arr2
        
    if nb1!=nb2 and nb1==1:
        out_arr = np.zeros(nb2,complex)
        i = arange(1,nb2+1)/arange(1,nb2+1)
        out_arr.real=i*arr1
        out_arr.imag=arr2

    if nb1!=nb2 and nb2==1:
        out_arr = np.zeros(nb1,complex)
        i = arange(1,nb1+1)/arange(1,nb1+1)
        out_arr.real=arr1
        out_arr.imag=i*arr2
        
    if nb1==1 and nb2==1:
        out_arr = np.zeros(nb1,complex)
        out_arr.real=arr1
        out_arr.imag=arr2
            
    return np.array(out_arr)
#######################################################################################################
#######################################################################################################
# read pointing directory list
#  input: filename
#  output: directory names
def read_ptgdirlist(filename):
    print ""
    print "[READ FPDB]: START reading "+filename

    dirs = []
    for line in fileinput.input(filename):
#        ar = line.split()
        dirs.append(int(ar[0]))
    print "[READ FPDB]: End reading "
    return dirs
#######################################################################################################
#######################################################################################################
# read the focal plane data base
#  input: filename
#  output: ch, az, el, amp, sig_x, sig_y, theta_tilt
def read_MuellerMatrix(filename, option_debug=False):
    if option_debug:
        MuellerMatrix = {'T1':float(1.), 
                         'T2':float(0.), 
                         'trans':float(1.), 
                         'rho':float(0.), 
                         'cosdel':float(-1.)}
        return MuellerMatrix

    if not option_debug:
        print "[READ FPDB]: START reading "+filename
        conn = sq.connect(filename)
        c = conn.cursor()
        c.execute('select * from MuellerParams;')
        boloid = []; boloname = []; T1 = []; T2 = []; trans = []; rho = []; cosdel =[]
        for ar in c:
            boloid.append(ar[0])
            boloname.append(ar[1])
            T1.append(float(ar[2]))
            T2.append(float(ar[3]))
            trans.append(float(ar[4]))
            rho.append(float(ar[5]))
            cosdel.append(float(ar[6]))
        c.close()

        MuellerMatrix = { 'boloid':boloid, \
                              'boloname':boloname, \
                              'T1':np.array(T1), \
                              'T2':np.array(T2), \
                              'trans':np.array(trans), \
                              'rho':np.array(rho), \
                              'cosdel':np.array(cosdel)}
        print "[READ MuellerMatrix]: End reading "
        return MuellerMatrix
#######################################################################################################
#######################################################################################################
# read the focal plane data base
#  input: filename
#  output: ch, az, el, amp, sig_x, sig_y, theta_tilt
def read_FPDB1(filename):
    print "[READ FPDB]: START reading "+filename

    ch = []
    az = []
    el = []
    amp = []
    sig_x = []
    sig_y = []
    theta_tilt = []

    i = 0
    for line in fileinput.input(filename):
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
    print "[READ FPDB]: End reading "
    FPDBlist = {'ch': np.array(ch), 'az': np.array(az), 'el': np.array(el), 'amp': np.array(amp),
               'sig_x': np.array(sig_x), 'sig_y': np.array(sig_y), 'theta_tilt': np.array(theta_tilt)}
    return FPDBlist
######################################################################################################
#######################################################################################################
# read the focal plane data base
#  input: filename
#  output: ch, az, el, amp, sig_x, sig_y, amp, poleff, polang, wafer, pix, torb, flag
def read_FPDB2(filename):
    print ""
    print "[READ FPDB]: START reading "+filename

    ch = []
    az = []
    el = []
    sig_x = []
    sig_y = []
    amp = []
    polang = []
    poleff = []
    wafer = []
    pix = []
    torb = []
    flag = []
    i = 0
    for line in fileinput.input(filename):
        ar = line.split()
        if ((len(ar)>1) and (i>0)):
            ch.append(int(ar[0]))
            az.append(float(ar[1]))
            el.append(float(ar[2]))
            sig_x.append(float(ar[3]))
            sig_y.append(float(ar[4]))
            amp.append(float(ar[5]))
            polang.append(float(ar[6]))
            poleff.append(float(ar[7]))
            wafer.append(str(ar[8]))
            pix.append(int(ar[9]))
            torb.append(str(ar[10]))
            flag.append(int(ar[11]))
        i += 1
    print "[READ FPDB]: End reading "
    FPDBlist = {'ch':np.array(ch),'az':np.array(az),'el': np.array(el),'sig_x':np.array(sig_x),'sig_y':np.array(sig_y),
                'amp': np.array(amp),'poleff':np.array(poleff),'polang':np.array(polang),
                'wafer':np.array(wafer),'pix':np.array(pix),'torb': np.array(torb),'flag':np.array(flag)}
    return FPDBlist
######################################################################################################
#######################################################################################################
# read FPDBlist
#  input: flag
#  output: 
def shuffle_pixlist_old(FPDBlist,flag_id):
    ch = FPDBlist["ch"]
    flag = FPDBlist["flag"]
    pix = FPDBlist["pix"]
    torb = FPDBlist["torb"]
    wafer = FPDBlist["wafer"]
#    wafer_id = ["8.2.0","8.2.1"]

    pix_id = []
    torb_id = []
    idx_id = []
#    for i in wafer_id:
#        ind = np.where((flag==flag_id) & (wafer==i))
    ind = np.where(flag==flag_id)
    ind = ind[0]

# commented below
#    for j in range(min(pix[ind]),max(pix[ind])):
#        idx = np.where((j==pix[ind]))
#        idx = idx[0]
#        if (len(idx)==2):
#            if "t" in torb[ind[idx[0]]]:
#                tmp1 = pix[ind[idx]]
#                tmp2 = torb[ind[idx]]
#                tmp3 = ch[ind[idx]]
#                pix_id.append(tmp1)
#                torb_id.append(tmp2)
#                idx_id.append(tmp3)
#            if "b" in torb[ind[idx[0]]]:
#                tmp1 = pix[ind[idx[::-1]]] 
#                tmp2 = torb[ind[idx[::-1]]] 
#                tmp3 = ch[ind[idx[::-1]]]
#                pix_id.append(tmp1)
#                torb_id.append(tmp2)
#                idx_id.append(tmp3)
    wafer_sel = wafer[ind]
    wafer_arr = []
    wafer_arr.append(wafer_sel[0])
    j = 0 # added
    for i in range(0,len(ind)):
        if np.array(wafer_arr[j]) != str(wafer_sel[i]):
#            print i, j, np.array(wafer_arr[j]), str(wafer_sel[i])
            wafer_arr.append(wafer_sel[i])
            j += 1
    num_wafer = len(wafer_arr)
#    print num_wafer, wafer_arr
#    for i in range(0,len(ind)):
#        if len(wafer_arr)==1:
#            if np.array(wafer_arr[j]) != str(wafer_sel[i]):
#                wafer_arr.append(wafer_sel[i])
#        else:
#            ind_arr = np.where((wafer_arr !=  wafer_sel[i]))
#            if len(ind_arr)!=0: wafer_arr.append(wafer_sel[i])
#    num_wafer = len(wafer_arr)

#    sys.exit()

    if num_wafer == 1:
        for j in range(min(pix[ind]),max(pix[ind])):
            idx = np.where((j==pix[ind]))
            idx = idx[0]
            if (len(idx)==2):
                if "t" in torb[ind[idx[0]]]:
                    tmp1 = pix[ind[idx]]
                    tmp2 = torb[ind[idx]]
                    tmp3 = ch[ind[idx]]
                    pix_id.append(tmp1)
                    torb_id.append(tmp2)
                    idx_id.append(tmp3)
                if "b" in torb[ind[idx[0]]]:
                    tmp1 = pix[ind[idx[::-1]]] 
                    tmp2 = torb[ind[idx[::-1]]] 
                    tmp3 = ch[ind[idx[::-1]]]
                    pix_id.append(tmp1)
                    torb_id.append(tmp2)
                    idx_id.append(tmp3)
    else:
        for i in wafer_arr:
            for j in range(min(pix[ind]),max(pix[ind])):
                idx = np.where((j==pix[ind]) & (i==wafer[ind]))
#                print idx, j, pix[ind], i, wafer[ind]
                idx = idx[0]
                if "t" in torb[ind[idx[0]]]:
                    tmp1 = pix[ind[idx]]
                    tmp2 = torb[ind[idx]]
                    tmp3 = ch[ind[idx]]
                    pix_id.append(tmp1)
                    torb_id.append(tmp2)
                    idx_id.append(tmp3)
                if "b" in torb[ind[idx[0]]]:
                    tmp1 = pix[ind[idx[::-1]]] 
                    tmp2 = torb[ind[idx[::-1]]] 
                    tmp3 = ch[ind[idx[::-1]]]
                    pix_id.append(tmp1)
                    torb_id.append(tmp2)
                    idx_id.append(tmp3)
#                print pix[idx[0]], wafer[idx[0]]
#    sys.exit()
#    print idx_id, len(idx_id)
#    sys.exit()
    out = np.array(idx_id) #-1
    return out
######################################################################################################
#######################################################################################################
# read the pointing file
#  input: boresight fits filename
#  output: {ra, dec, pa, hwp}
def read_ptg(filename):
    print ""
    print "[READ PTG]: BEGIN reading "+filename 
    ptg = h.mrdfits(filename,hdu=1)
    ptg_package = {'pa':ptg[0], 'flag':ptg[1], 'dec':ptg[2], 'ra':ptg[3], 'hwp':ptg[4]}
#    print 'HWP angle [degs]', ptg[4]/pi*180.
    print "[READ PTG]: END reading "
    return ptg_package
######################################################################################################
#######################################################################################################
# read the level2 pointing file
#  input: boresight fits filename
#  output: {ra, dec, pa, hwp}
def read_level2_ptg(filename):
    print ""
    print "[READ PTG]: BEGIN reading "+filename 
    ptg = np.load(filename)
    ptg_package = {'pa':ptg[2], 'dec':ptg[1], 'ra':ptg[0], 'hwp':ptg[3]}    
    print 'HWP angle [degs]', ptg[3]/pi*180.
    print "[READ PTG]: END reading "
    return ptg_package
######################################################################################################
#######################################################################################################
# read pointing index
#  input: filename
#  output: 
def read_ptg_idx(filename):
    print ""
    print "[READ Pointing Flag text]: Start reading "+filename
    idx = []
    idx_s = []
    idx_e = []
    direc = []
    flag = []

    i = 0
    for line in fileinput.input(filename):
        ar = line.split()
        if ((len(ar)>1) & (i>0)):
            idx.append(int(ar[0]))
            idx_s.append(int(ar[1]))
            idx_e.append(int(ar[2]))
            direc.append(int(ar[3]))
            flag.append(int(ar[4]))
        i += 1
    print "[READ Pointing Flag text]: End reading "
    ptg_flag_package = {'idx': np.array(idx), 's_idx': np.array(idx_s), 'e_idx': np.array(idx_e), 'direc': np.array(direc), 'flag_hscan': np.array(flag)}
    return ptg_flag_package
######################################################################################################
#######################################################################################################
# read pointing index
#  input: filename
#  output: 
def read_level2_ptg_idx(filename):
    print ""
    print "[READ Pointing Flag text]: Start reading "+filename
    idx = []
    idx_s = []
    idx_e = []
    direc = []
    flag = []
    vel_mean = []
    frac = []

    conn = sq.connect(filename)
    c = conn.cursor()
    c.execute('select * from ptgflag;') 
    for ar in c:
        idx.append(int(ar[0]))
        idx_s.append(int(ar[1]))
        idx_e.append(int(ar[2]))
        direc.append(int(ar[3]))
        flag.append(int(ar[4]))
#        vel_mean.append(float(ar[5]))
#        frac.append(float(ar[6]))
    c.close()

    print "[READ Pointing Flag text]: End reading "
    ptg_flag_package = {'idx': np.array(idx), 's_idx': np.array(idx_s), 'e_idx': np.array(idx_e), 
                        'direc': np.array(direc), 'flag_hscan': np.array(flag)}
#                        'vel_mean':vel_mean, 'frac':frac}
    return ptg_flag_package
######################################################################################################
#######################################################################################################
# flag the portion of the scan going outside of the specified range
def flag_hscan(ptg, idx, xml_input):
    nbData = len(idx['idx'])
    phi_c = float(xml_input['patch_ra'])/180.*pi
    theta_c = pi/2.-float(xml_input['patch_dec'])/180.*pi
    width = float(xml_input['patch_width'])/180.*pi

    for i in range(nbData):
        idx_s = idx['s_idx'][i]
        idx_e = idx['e_idx'][i]

        phi = ptg['ra'][idx_s:idx_e]
        theta = pi/2. - ptg['dec'][idx_s:idx_e]

        if ( ((phi_c < pi) & (phi_c > width))
             or ((phi_c>=pi) & (abs(2.*pi-phi_c)>width) )):
            pass
        if ( (phi_c < pi) & (phi_c < width) ):
            ind = np.where(phi>pi)
            phi[ind[0]] = phi[ind[0]]-2.*pi
        if ( (phi_c > pi) & (abs(2.*pi-phi_c) < width) ):
            ind = np.where(phi<pi)
            phi[ind[0]] = phi[ind[0]]+2.*pi

        ind_part = np.where(((theta>theta_c-width/2.) & (theta<theta_c+width/2.))
                            & ((phi>phi_c-width/2.) & (phi<phi_c+width/2.)))
        if len(ind_part[0]) != len(theta):
            idx['flag_hscan'][i]=1
            print '[lib_mapmaker.py/flag_hscan] flag the half-scan', idx['flag_hscan'][i]
#        ind_ra = np.where((ptg['ra']>ra-w/2.) & (ptg['ra']<ra+w/2.))
#        ind_dec = np.where((ptg['dec']>dec-w/2.) & (ptg['dec']<dec+w/2.))
#        if len(ind_ra[0])!=0: 
#            idx['flag_hscan'][i]=1
#        if len(ind_dec[0])!=0: 
#            idx['flag_hscan'][i]=1
    return idx
######################################################################################################
#######################################################################################################
# read simulation noise 
#  input: filename
#  output: 
def read_simnoise(filename,option_debug=False):
    if option_debug:
        num = 1511
        zero = np.zeros(num)
        one = np.ones(num)
        simnoise_package = {'idx': zero, 'net': zero, 'fknee': one, 'power': zero, 'seed': zero}
        return simnoise_package

    print ""
    print "[READ simulation noise]: Start reading "+filename
    idx = []
    net = []
    fknee = []
    pow = []
    seed = []

    i = 0
    for line in fileinput.input(filename):
        ar = line.split()
        if ((len(ar)>1) & (i>0)):
            idx.append(int(ar[0]))
            net.append(float(ar[1]))
            fknee.append(float(ar[2]))
            pow.append(float(ar[3]))
            seed.append(int(ar[4]))
        i += 1
    print "[READ Pointing Flag text]: End reading "
    simnoise_package = {'idx': np.array(idx), 'net': np.array(net), 'fknee': np.array(fknee), 'power': np.array(pow), 'seed': np.array(seed)}
    return simnoise_package
######################################################################################################
#######################################################################################################
# read simulation noise 
#  input: filename
#  output: 
#def read_relgain(filename,option_debug=False):
#    if option_debug:
#        relgain_package = {'net': 0., 
#                           'fknee': 1., 
#                           'power': 1., 
#                           'randshift': 1.} # 0: random, 1: shift
#        return relgain_package
#        
#    print ""
#    print "[READ simulation noise]: Start reading "+filename
#    net = []
#    fknee = []
#    pow = []
#    randshift = []
#    i = 0
#    for line in fileinput.input(filename):
#        ar = line.split()
#        if ((len(ar)>1) & (i>0)):
#            net.append(float(ar[0]))
#            fknee.append(float(ar[1]))
#            pow.append(float(ar[2]))
#            randshift.append(int(ar[3]))
#        i += 1
#    print "[READ Pointing Flag text]: End reading "
#    relgain_package = {'net': net, 'fknee': fknee, 'power': pow, 'randshift':randshift}
#    return relgain_package
######################################################################################################
#######################################################################################################
# read simulation relative gain
#  input: filename
#  output: 
def read_relgain(filename):
    conn = sq.connect(filename)
    c = conn.cursor()
    c.execute('select * from pb1_relgain;') 
    for ar in c:
        boloid.append(int(ar[0]))
        boloname.append(ar[1])
        relgain.append(float(ar[2]))
    c.close()

    print "[READ real_relgaindb]: End reading "
    relgain_package = {'boloid': np.array(boloid), 
                       'boloname': np.array(boloname), 
                       'relgain': np.array(relgain)}
    return relgain_package
#######################################################################################################
#######################################################################################################
# generate the pixel pointing from boresight pointing
#  input: pixel number, pixel list, pointing package, fpdb_list
#  output: 
def boreptg2pixptg_c(i_pix,pix_list,pointing_package,fpdb_list,time0,option_silent):
    if option_silent==False: print ""
    if option_silent==False: print "=================================="
    if option_silent==False: info("[BOREPTG2PIXPTG] START, pixel#="+str(i_pix), time0)

    b_ra = pointing_package["ra"]
    b_dec = pointing_package["dec"]
    b_dk = pointing_package["pa"]
    b_hwp = pointing_package["hwp"]
    
#    az = fpdb_list['az']*pi/180.
#    el = fpdb_list['el']*pi/180.
    az = fpdb_list['xpos']*pi/180.
    el = fpdb_list['ypos']*pi/180.
#    theta_tilt = fpdb_list['theta_tilt']*pi/180.
    polang = fpdb_list['polang']*pi/180.

#    print fpdb_list
#    print az
#    print el

    code = '''
    for (i=0;i<2;i++):
        ii = int(pix_list[i_pix][i]) 

        r_psb = sqrt(az[ii]*az[ii]+el[ii]*el[ii])
        th_psb = atan(el[ii]/az[ii])
        chi_psb = -(th_psb-polang[ii])
        /* Calculate dec first  */
        theta0 = pi/2. - b_dec
        theta = cos(r_psb)*cos(theta0) + sin(r_psb)*sin(theta0)*cos(b_dk+th_psb-pi/2.)
/*        ind = np.where(theta > 1.)
        theta[ind[0]] = 1.
        ind = np.where(theta < -1.)
        theta[ind[0]] = -1. */
        theta = acos(theta)
        dec_psb = pi/2.-theta
    
        /* Now do RA */
        phi = ( cos(r_psb) - cos(theta)*cos(theta0) ) / ( sin(theta)*sin(theta0) )
/*        ind = np.where(phi > 1.)
        phi[ind[0]] = 1.
        ind = np.where(phi < -1.)
        phi[ind[0]] = -1. */
        phi = acos(phi) 
        ind = np.where(np.sin(th_psb-b_dk-pi/2.) < 0.0)
        phi[ind[0]] = -1.*phi[ind[0]]
        ra_psb = phi+b_ra

/*        ind = np.where(ra_psb < 0.0) */
/*        ra_psb[ind[0]] = ra_psb[ind[0]] + 2.*pi */
        while (ra_psb < 0.0):
            ra_psb = ra_psb + 2.*pi
        ra_psb = ra_psb % (2.*pi)

        if (r_psb == 0.0):
            psi_ra = chi_psb + th_psb - b_dk - pi/2.;
            while (psi_ra < 2.*pi):
                psi_ra = psi_ra + 2.*pi
/*            ind = np.where(psi_ra < 2.*pi) */
/*            psi_ra[ind[0]] = psi_ra[ind[0]] + 2.*pi  */
            psi_ra = psi_ra % pi


        beta = ( cos(theta0) - cos(theta)*cos(r_psb) ) / ( sin(theta)*sin(r_psb) )
        if (beta > 1.0): beta = 1.0
        if (beta < -1.0): beta = -1.0

/*        ind = np.where(beta > 1.0)  */
/*        beta[ind[0]] = 1.  */
/*        ind = np.where(beta < -1.0)  */
/*        beta[ind[0]] = -1.  */
        beta = acos(beta)
        if (sin(th_psb-b_dk-pi/2.) < 0.0): beta = -1*beta
/*        ind = np.where(np.sin(th_psb-b_dk-pi/2.) < 0.0)  */
/*        beta[ind[0]] = -1.*beta[ind[0]]   */
        psi_ra = ( pi - beta + chi_psb) % pi

        if (i==0): top_ptg = {'ra': np.array(ra_psb), 'dec':np.array(dec_psb), 'psi':np.array(psi_ra)}
        if (i==1): bot_ptg = {'ra': np.array(ra_psb), 'dec':np.array(dec_psb), 'psi':np.array(psi_ra)}
'''

#    pointing = {"top_ptg": top_ptg, "bot_ptg": bot_ptg}
#    return pointing
    if option_silent==False: print "[BOREPTG2PIXPTG] END, pixel#=", str(i_pix)
    if option_silent==False: print "=================================="
    return top_ptg, bot_ptg
#######################################################################################################
#######################################################################################################
# generate the pixel pointing from boresight pointing
#  input: pixel number, pixel list, pointing package, fpdb_list
#  output: 
def boreptg2pixptg(i_pix,pix_list,pointing_package,fpdb_list,time0,option_silent):
    if option_silent==False: print ""
    if option_silent==False: print "=================================="
    if option_silent==False: info("[BOREPTG2PIXPTG] START, pixel#="+str(i_pix), time0)

    b_ra = pointing_package["ra"]
    b_dec = pointing_package["dec"]
    b_dk = pointing_package["pa"]
    b_hwp = pointing_package["hwp"]
    
#    az = fpdb_list['az']*pi/180.
#    el = fpdb_list['el']*pi/180.
    az = fpdb_list['xpos']*pi/180.
    el = fpdb_list['ypos']*pi/180.
#    theta_tilt = fpdb_list['theta_tilt']*pi/180.
    polang = fpdb_list['polang']*pi/180.

#    print fpdb_list
#    print az
#    print el

    for i in range(0,2):
        ii = int(pix_list[i_pix][i]) #fpdb_list['boloid'][i]]

        r_psb = np.sqrt(az[ii]**2.+el[ii]**2.)
        th_psb = np.arctan(el[ii]/az[ii])
        chi_psb = -(th_psb-polang[ii])
#        print type(az), type(el), len(az), len(el)
#        print i, i_pix, ii, az[ii], el[ii], polang[ii], r_psb, th_psb, chi_psb
        # chi_psb = polang[ii]

        # Calculate dec first  
        theta0 = pi/2. - b_dec
#        theta = np.cos(r_psb)*np.cos(theta0) + np.sin(r_psb)*np.sin(theta0)*np.cos(th_psb-b_dk-pi/2.)
        theta = np.cos(r_psb)*np.cos(theta0) + np.sin(r_psb)*np.sin(theta0)*np.cos(b_dk+th_psb-pi/2.)
        ind = np.where(theta > 1.)
        theta[ind[0]] = 1.
        ind = np.where(theta < -1.)
        theta[ind[0]] = -1.
        theta = np.arccos(theta)
        dec_psb = pi/2.-theta
    
        # Now do RA */
        phi = ( np.cos(r_psb) - np.cos(theta)*np.cos(theta0) ) / ( np.sin(theta)*np.sin(theta0) )
        ind = np.where(phi > 1.)
        phi[ind[0]] = 1.
        ind = np.where(phi < -1.)
        phi[ind[0]] = -1.
        phi = np.arccos(phi) # watch out for acos domain
        ind = np.where(np.sin(th_psb-b_dk-pi/2.) < 0.0)
        phi[ind[0]] = -1.*phi[ind[0]]
        ra_psb = phi+b_ra
        ind = np.where(ra_psb < 0.0)
        ra_psb[ind[0]] = ra_psb[ind[0]] + 2.*pi
#        while (ra_psb < 0.0):
#            ra_psb = ra_psb + 2.*pi
        ra_psb = ra_psb % (2.*pi)

        if (r_psb == 0.0):
            psi_ra = chi_psb + th_psb - b_dk - pi/2.;
#            while (psi_ra < 2.*pi):
#                psi_ra = psi_ra + 2.*pi
            ind = np.where(psi_ra < 2.*pi)
            psi_ra[ind[0]] = psi_ra[ind[0]] + 2.*pi                
            psi_ra = psi_ra % pi
#            print "one before return"
#            return

        beta = ( np.cos(theta0) - np.cos(theta)*np.cos(r_psb) ) / ( np.sin(theta)*np.sin(r_psb) )
#        if (beta > 1.0): beta = 1.0
#        if (beta < -1.0): beta = -1.0
        ind = np.where(beta > 1.0)
        beta[ind[0]] = 1.
        ind = np.where(beta < -1.0)
        beta[ind[0]] = -1.
        beta = np.arccos(beta)
#        if (sin(th_psb-b_dk-pi/2.) < 0.0): beta = -1*beta
        ind = np.where(np.sin(th_psb-b_dk-pi/2.) < 0.0)
        beta[ind[0]] = -1.*beta[ind[0]]
        psi_ra = ( pi - beta + chi_psb) % pi

        if (i==0): top_ptg = {'ra':np.array(ra_psb), 'dec':np.array(dec_psb), 'psi':np.array(psi_ra)}
        if (i==1): bot_ptg = {'ra':np.array(ra_psb), 'dec':np.array(dec_psb), 'psi':np.array(psi_ra)}

#    pointing = {"top_ptg": top_ptg, "bot_ptg": bot_ptg}
#    return pointing
    if option_silent==False: print "[BOREPTG2PIXPTG] END, pixel#=", str(i_pix)
    if option_silent==False: print "=================================="
    return top_ptg, bot_ptg
#######################################################################################################
#######################################################################################################
# generate the pixel pointing from boresight pointing
#  input: pixel number, pixel list, pointing package, fpdb_list
#  output: 
def boreptg2pixptg_inC(i_pix,pix_list,pointing_package,fpdb_list,time0,option_silent):
    if option_silent==False: print ""
    if option_silent==False: print "=================================="
    if option_silent==False: info("[BOREPTG2PIXPTG] START, pixel#="+str(i_pix), time0)

    b_ra = pointing_package["ra"]
    b_dec = pointing_package["dec"]
    b_dk = pointing_package["pa"]
    b_hwp = pointing_package["hwp"]
    
    az = fpdb_list['xpos']*pi/180.
    el = fpdb_list['ypos']*pi/180.
    polang = fpdb_list['polang']*pi/180.

    for i in range(0,2):
        ii = int(pix_list[i_pix][i]) 

        r_psb = np.sqrt(az[ii]**2.+el[ii]**2.)
        th_psb = np.arctan(el[ii]/az[ii])
        chi_psb = -(th_psb-polang[ii])
        # Calculate dec first  
        theta0 = pi/2. - b_dec
#        theta = np.cos(r_psb)*np.cos(theta0) + np.sin(r_psb)*np.sin(theta0)*np.cos(th_psb-b_dk-pi/2.)
        theta = np.cos(r_psb)*np.cos(theta0) + np.sin(r_psb)*np.sin(theta0)*np.cos(b_dk+th_psb-pi/2.)
        ind = np.where(theta > 1.)
        theta[ind[0]] = 1.
        ind = np.where(theta < -1.)
        theta[ind[0]] = -1.
        theta = np.arccos(theta)
        dec_psb = pi/2.-theta
    
        # Now do RA */
        phi = ( np.cos(r_psb) - np.cos(theta)*np.cos(theta0) ) / ( np.sin(theta)*np.sin(theta0) )
        ind = np.where(phi > 1.)
        phi[ind[0]] = 1.
        ind = np.where(phi < -1.)
        phi[ind[0]] = -1.
        phi = np.arccos(phi) # watch out for acos domain
        ind = np.where(np.sin(th_psb-b_dk-pi/2.) < 0.0)
        phi[ind[0]] = -1.*phi[ind[0]]
        ra_psb = phi+b_ra
        ind = np.where(ra_psb < 0.0)
        ra_psb[ind[0]] = ra_psb[ind[0]] + 2.*pi
        ra_psb = ra_psb % (2.*pi)

        if (r_psb == 0.0):
            psi_ra = chi_psb + th_psb - b_dk - pi/2.;
            ind = np.where(psi_ra < 2.*pi)
            psi_ra[ind[0]] = psi_ra[ind[0]] + 2.*pi                
            psi_ra = psi_ra % pi

        beta = ( np.cos(theta0) - np.cos(theta)*np.cos(r_psb) ) / ( np.sin(theta)*np.sin(r_psb) )
        ind = np.where(beta > 1.0)
        beta[ind[0]] = 1.
        ind = np.where(beta < -1.0)
        beta[ind[0]] = -1.
        beta = np.arccos(beta)
        ind = np.where(np.sin(th_psb-b_dk-pi/2.) < 0.0)
        beta[ind[0]] = -1.*beta[ind[0]]
        psi_ra = ( pi - beta + chi_psb) % pi

        if (i==0): top_ptg = {'ra': np.array(ra_psb), 'dec':np.array(dec_psb), 'psi':np.array(psi_ra)}
        if (i==1): bot_ptg = {'ra': np.array(ra_psb), 'dec':np.array(dec_psb), 'psi':np.array(psi_ra)}

    if option_silent==False: print "[BOREPTG2PIXPTG] END, pixel#=", str(i_pix)
    if option_silent==False: print "=================================="
    return top_ptg, bot_ptg
#######################################################################################################
#######################################################################################################
def euler_astrolib(ai,bi,select, FK4=False, radian=False):
    '''
    ; NAME:
    ;     EULER
    ; PURPOSE:
    ;     Transform between Galactic, celestial, and ecliptic coordinates.
    ; EXPLANATION:
    ;     Use the procedure ASTRO to use this routine interactively
    ;
    ; CALLING SEQUENCE:
    ;      EULER, AI, BI, AO, BO, [ SELECT, /FK4 ]
    ;
    ; INPUTS:
    ;       AI - Input Longitude in DEGREES, scalar or vector.  If only two
    ;               parameters are supplied, then  AI and BI will be modified to
    ;               contain the output longitude and latitude.
    ;       BI - Input Latitude in DEGREES
    ;
    ; OPTIONAL INPUT:
    ;       SELECT - Integer (1-6) specifying type of coordinate transformation.
    ;
    ;      SELECT   From          To        |   SELECT      From            To
    ;       1     RA-Dec (2000)  Galactic   |     4       Ecliptic      RA-Dec
    ;       2     Galactic       RA-DEC     |     5       Ecliptic      Galactic
    ;       3     RA-Dec         Ecliptic   |     6       Galactic      Ecliptic
    ;
    ;      If omitted, program will prompt for the value of SELECT
    ;      Celestial coordinates (RA, Dec) should be given in equinox J2000
    ;      unless the /FK4 keyword is set.
    ; OUTPUTS:
    ;       AO - Output Longitude in DEGREES
    ;       BO - Output Latitude in DEGREES
    ;
    ; INPUT KEYWORD:
    ;       /FK4 - If this keyword is set and non-zero, then input and output
    ;             celestial and ecliptic coordinates should be given in equinox
    ;             B1950.
    ;
    ; NOTES:
    ;       EULER was changed in December 1998 to use J2000 coordinates as the
    ;       default, ** and may be incompatible with earlier versions***.
    ; REVISION HISTORY:
    ;       Written W. Landsman,  February 1987
    ;       Adapted from Fortran by Daryl Yentis NRL
    ;       Converted to IDL V5.0   W. Landsman   September 1997
    ;       Made J2000 the default, added /FK4 keyword  W. Landsman December 1998
    ;-
    '''
    
    twopi   =   2.0*pi
    fourpi  =   4.0*pi
    deg_to_rad = 180.0/pi
    if radian: deg_to_rad = 1.
    
#;   J2000 coordinate conversions are based on the following constants
#;  eps = 23.4392911111d              Obliquity of the ecliptic
#;  alphaG = 192.85948d               Right Ascension of Galactic North Pole
#;  deltaG = 27.12825d                Declination of Galactic North Pole
#;  lomega = 32.93192d                Galactic longitude of celestial equator
#;  alphaE = 180.02322d              Ecliptic longitude of Galactic North Pole
#;  deltaE = 29.811438523d            Ecliptic latitude of Galactic North Pole
#;  Eomega  = 6.3839743d              Galactic longitude of ecliptic equator

    if FK4:
        equinox = '(B1950)'
        psi   = np.array([ 0.57595865315, 4.9261918136, 0.00000000000, 0.0000000000, 0.11129056012, 4.7005372834])
        stheta = np.array([ 0.88781538514,-0.88781538514, 0.39788119938,-0.39788119938, 0.86766174755,-0.86766174755])
        ctheta = np.array([ 0.46019978478, 0.46019978478, 0.91743694670, 0.91743694670, 0.49715499774, 0.49715499774])
        phi  = np.array([ 4.9261918136,  0.57595865315, 0.0000000000, 0.00000000000, 4.7005372834, 0.11129056012])
        
    equinox = '(J2000)'
    psi   =  np.array([ 0.57477043300, 4.9368292465, 0.00000000000, 0.0000000000, 0.11142137093, 4.71279419371])
    stheta = np.array([ 0.88998808748,-0.88998808748, 0.39777715593,-0.39777715593, 0.86766622025,-0.86766622025])
    ctheta = np.array([ 0.45598377618, 0.45598377618, 0.91748206207, 0.91748206207, 0.49714719172, 0.49714719172])
    phi  = np.array([ 4.9368292465,  0.57477043300, 0.0000000000, 0.00000000000, 4.71279419371, 0.11142137093])

    i  = select - 1                         # IDL offset
    a  = ai/deg_to_rad - phi[i]
    b = bi/deg_to_rad
    sb = np.sin(b)
    cb = np.cos(b)
    cbsa = cb * np.sin(a)
    b  = -stheta[i] * cbsa + ctheta[i] * sb
    bo = np.arcsin(b)*deg_to_rad
    
    a =  np.arctan2( ctheta[i] * cbsa + stheta[i] * sb, cb * np.cos(a) )
    ao = ( (a+psi[i]+fourpi) % twopi) * deg_to_rad

    return ao, bo
#######################################################################################################
#######################################################################################################
def Read_bolos():
    pass
#######################################################################################################
#######################################################################################################
#    top_tod, bot_tod = SignalGen(i_pix, packedPTG, SimInputs, SysInputs)
# GENERATE THE SIGNAL
def SignalGen(i_pix, pix_list, top_ptg, bot_ptg, hwpang, SimInputs, MuellerMatrix, time0, option_TQU, option_silent):
    if option_silent==False: print ""
    if option_silent==False: print "=================================="
    if option_silent==False: info("[SIGNALGEN] START, pixel#="+str(i_pix), time0)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ind = SimInputs['ind_band_in']
    idx_min_in = SimInputs['ind_min_in']
    idx_max_in = SimInputs['ind_max_in']

    nside_in = SimInputs['nside']
    mapT = SimInputs['T']
    if option_TQU == "TQU":
        mapQ = SimInputs['Q']
        mapU = SimInputs['U']

    ii = int(pix_list[i_pix][0])
    T1 = float(MuellerMatrix['T1'][ii])
    T2 = float(MuellerMatrix['T2'][ii])
    trans = float(MuellerMatrix['trans'][ii])
    rho = float(MuellerMatrix['rho'][ii])
    cosdel = float(MuellerMatrix['cosdel'][ii])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
# TOP
    top_ra = top_ptg['ra']
    top_dec = pi/2. -  top_ptg['dec']
    top_psi = top_ptg['psi']

    top_ipix = h.ang2pix(nside_in, top_dec, top_ra)
    top_ipix = np.int_(top_ipix-idx_min_in)
    if option_TQU=="T":
        M11 = 0.5*( (T1+T2)*trans + rho*(T1-T2)*np.cos(2.*hwpang-2.*top_psi) )
        top_tod = M11*mapT[ind[top_ipix]]
    if option_TQU=="TQU":
        M11 = 0.5*( (T1+T2)*trans + rho*(T1-T2)*np.cos(2.*hwpang-2.*top_psi) )
        M12 = 0.5*( (T1+T2)*rho*np.cos(2.*hwpang) \
                        + 0.5*(T1-T2)*((trans+cosdel)*np.cos(2.*top_psi)+(trans-cosdel)*np.cos(4.*hwpang-2.*top_psi)))
        M13 = 0.5*(-(T1+T2)*rho*np.sin(2.*hwpang) \
                        + 0.5*(T1-T2)*((trans+cosdel)*np.cos(2.*top_psi)+(trans-cosdel)*np.sin(4.*hwpang-2.*top_psi)))
        top_tod = M11*mapT[ind[top_ipix]]+M12*mapQ[ind[top_ipix]]+M13*mapU[ind[top_ipix]]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BOTTOM
    ii = int(pix_list[i_pix][1]) 
    T1 = float(MuellerMatrix['T1'][ii])
    T2 = float(MuellerMatrix['T2'][ii])
    trans = float(MuellerMatrix['trans'][ii])
    rho = float(MuellerMatrix['rho'][ii])
    cosdel = float(MuellerMatrix['cosdel'][ii])
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    bot_ra = bot_ptg['ra']
    bot_dec = pi/2. -  bot_ptg['dec']
    bot_psi = bot_ptg['psi']
    
    bot_ipix = h.ang2pix(nside_in, bot_dec, bot_ra)
    bot_ipix = np.int_(bot_ipix-idx_min_in)
    if option_TQU=="T":
        M11 = 0.5*( (T1+T2)*trans + rho*(T1-T2)*np.cos(2.*hwpang-2.*bot_psi) )
        bot_tod = M11*mapT[ind[bot_ipix]]
    if option_TQU=="TQU":
        M11 = 0.5*( (T1+T2)*trans + rho*(T1-T2)*np.cos(2.*hwpang-2.*bot_psi) )
        M12 = 0.5*( (T1+T2)*rho*np.cos(2.*hwpang) \
                        + 0.5*(T1-T2)*((trans+cosdel)*np.cos(2.*bot_psi)+(trans-cosdel)*np.cos(4.*hwpang-2.*bot_psi)))
        M13 = 0.5*(-(T1+T2)*rho*np.sin(2.*hwpang) \
                        + 0.5*(T1-T2)*((trans+cosdel)*np.cos(2.*bot_psi)+(trans-cosdel)*np.sin(4.*hwpang-2.*bot_psi)))
        bot_tod = M11*mapT[ind[bot_ipix]]+M12*mapQ[ind[bot_ipix]]+M13*mapU[ind[bot_ipix]]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    if option_silent==False: print "[SIGNALGEN]: END, ", str(i_pix)
    if option_silent==False: print "=================================="
    return np.array(top_tod), np.array(bot_tod)
#######################################################################################################
#######################################################################################################
#def NoiseGen_auto(nbData, NoiseInput, seed):
def NoiseGen_auto(nbData, net, knee, power, seed, time0, option_silent):
    if option_silent==False: info("[NoiseGen_auto]", time0)
    return 0.,0.
#    print ""
#    print "[NoiseGen_auto]: START"
    
#    a0 = NoiseInput['white_noise']
#    a0 = NoiseInput['net']
#    a1 = NoiseInput['knee']
#    a2 = NoiseInput['noise_spectrum']
#    a2 = NoiseInput['power']

#    print '>>>', net, knee, power, seed

#    nbf = int(nbData/2.)+1
#    freq = np.arange(nbf)
#    psd = net*(1.+(knee/freq)**power)
#    idx = range(1,nbf-1)
#    real_psd = np.concatenate(( np.array([psd[1]]), psd[idx], np.array([psd[nbf-1]]), psd[idx[::-1]] ))
#    imag_psd = np.concatenate(( np.array([0.]),   psd[idx], np.array([0.]),        -psd[idx[::-1]] ))
#
#    for i in range(0,2):
#        random.seed(seed+2*i)
#        f_rand = random.randn(2*(nbData+1))
#  #        print "len", len(f_rand), nbData, type(f_rand)
#  #        print len(real_psd), len(imag_psd), len(f_rand[0:nbData]), len(f_rand[nbData+1:(nbData+1)*2-1])
#  #        print type(real_psd), type(f_rand[0:nbData])
#  #        print len(real_psd), len(f_rand[0:nbData]), len(imag_psd), len(f_rand[nbData+1:(nbData+1)*2-1])
#        real_psdout = real_psd*f_rand[0:nbData]
#        imag_psdout = imag_psd*f_rand[nbData+1:(nbData+1)*2-1]
#        psd_complex = complex_arr(real_psdout,imag_psdout)
#        if (i==0): top_noise = fft.fft(psd_complex)
#        if (i==1): bot_noise = fft.fft(psd_complex)
#
#  #    print "[NoiseGen_auto]: END"
#    return np.array(top_noise.real), np.array(bot_noise.real)
#######################################################################################################
#######################################################################################################
#def SignalFilt(top_tod, bot_tod, top_gain, bot_gain, poly_n):
#    sum_tod = 0.5*(top_tod*top_gain + bot_tod*top_gain)
#    dif_tod = 0.5*(top_tod*bot_gain - bot_tod*bot_gain)
def RelgainCorr(top_tod, bot_tod, RelgainNoiseInput, time0):
    if option_silent==False: info("[RelgainCorr]", time0)
    frac = RelgainNoiseInput['net']
    knee = RelgainNoiseInput['fknee']
    power = RelgainNoiseInput['power']
    randshift = RelgainNoiseInput['randshift']

    nbData = len(top_tod)

    if randshift[0] == 0:
        seed = int(time.time()*1e6)
        random.seed(seed)
        relgain = frac[0]*random.randn(nbData)
        top_tod = top_tod*(1.+relgain)
        
        seed = int(2.*time.time()*1e6)
        random.seed(seed)
        relgain = frac[0]*random.randn(nbData)
        bot_tod = bot_tod*(1.+relgain)
        
    if randshift[0] == 1:
        bot_tod = bot_tod*(1.+frac[0])

    return top_tod, bot_tod
#######################################################################################################
#######################################################################################################
#def SignalFilt(top_tod, bot_tod, top_gain, bot_gain, poly_n):
#    sum_tod = 0.5*(top_tod*top_gain + bot_tod*top_gain)
#    dif_tod = 0.5*(top_tod*bot_gain - bot_tod*bot_gain)
def SignalFilt(top_tod, bot_tod, poly_n, time0, option_silent):
    
    nb = len(top_tod)
    sum_tod = 0.5*(top_tod + bot_tod)
    dif_tod = 0.5*(top_tod - bot_tod)
    if (poly_n == -1):
        sum_fit = 0.
        dif_fit = 0.
    if (poly_n == 0):
        sum_fit = np.mean(sum_tod)
        dif_fit = np.mean(dif_tod)
    sum_fit = np.zeros(nb)
    dif_fit = np.zeros(nb)
    if (poly_n > 0):
        idx = np.arange(0,nb)
        sum_fitpar = np.polyfit(idx, sum_tod, poly_n)
        dif_fitpar = np.polyfit(idx, dif_tod, poly_n)
        for i in range(poly_n,-1,-1):
            sum_fit += sum_fit + sum_fitpar[poly_n-i]*idx**i
            dif_fit += dif_fit + dif_fitpar[poly_n-i]*idx**i
    sum_tod_filt = sum_tod - sum_fit
    dif_tod_filt = dif_tod - dif_fit

    if (poly_n != -1):
        sum_std = np.std(sum_tod_filt)
        dif_std = np.std(dif_tod_filt)

    if (poly_n == -1):
        sum_std = 1.
        dif_std = 1.
    
#    if sum_std == 0: sum_std = 1.e20
#    if dif_std == 0: dif_std = 1.e20

    packedTOD = {'sum_tod':sum_tod_filt, 'dif_tod':dif_tod_filt, 'sum_std':sum_std, 'dif_std':dif_std}
    return packedTOD
#######################################################################################################
#######################################################################################################
def ConcatenateTOD_init():
    sum_tod = [];     dif_tod = [];
    sum_std = [];     dif_std = [];
    packedTOD = {'sum_tod':sum_tod, 'dif_tod':dif_tod, 'sum_std':sum_std, 'dif_std':dif_std}
    return packedTOD
#######################################################################################################
#######################################################################################################
def ConcatenateTOD(packedTOD_out, packedTOD, time0, option_silent):
    if option_silent==False: info("[ConcatenateTOD]", time0)
    sum_tod_out = packedTOD_out["sum_tod"]
    dif_tod_out = packedTOD_out["dif_tod"]
    sum_tod = packedTOD["sum_tod"]
    dif_tod = packedTOD["dif_tod"]
    nb = len(sum_tod)

    sum_tod_out_tmp = np.concatenate((sum_tod_out,sum_tod))
    dif_tod_out_tmp = np.concatenate((dif_tod_out,dif_tod))

    sum_std_out = packedTOD_out["sum_std"]
    dif_std_out = packedTOD_out["dif_std"]
    sum_std = packedTOD["sum_std"]
    dif_std = packedTOD["dif_std"]
#    print sum_std_out, sum_std
    sum_std_out_tmp = np.concatenate((np.array(sum_std_out),np.ones(nb)*sum_std))
    dif_std_out_tmp = np.concatenate((np.array(dif_std_out),np.ones(nb)*dif_std))

    packedTOD_ = {'sum_tod':sum_tod_out_tmp, 'dif_tod':dif_tod_out_tmp, 'sum_std':sum_std_out_tmp, 'dif_std':dif_std_out_tmp}
    return packedTOD_
#######################################################################################################
#######################################################################################################
def ConcatenatePTG_init():
    top_ra = []; top_dec = []; top_pa = [];
    bot_ra = []; bot_dec = []; bot_pa = [];
    hwpang = []
    packedPTG = {"top_ra":top_ra, "top_dec":top_dec, "top_pa":top_pa, "bot_ra":bot_ra, "bot_dec":bot_dec, "bot_pa":bot_pa, "hwp":hwpang}
    return packedPTG
#######################################################################################################
#######################################################################################################
def ConcatenatePTG(packedPTG_out, top_ptg, bot_ptg, hwpang, i_tod, time0, option_silent):
    if option_silent==False: info("[ConcatenatePTG]", time0)
    top_ra = top_ptg["ra"]
    top_dec = top_ptg["dec"]
    top_pa = top_ptg["psi"]
    bot_ra = bot_ptg["ra"]
    bot_dec = bot_ptg["dec"]
    bot_pa = bot_ptg["psi"]
     
    top_ra_ = packedPTG_out["top_ra"]
    top_dec_ = packedPTG_out["top_dec"]
    top_pa_ = packedPTG_out["top_pa"]
    bot_ra_ = packedPTG_out["bot_ra"]
    bot_dec_ = packedPTG_out["bot_dec"]
    bot_pa_ = packedPTG_out["bot_pa"]
    hwpang_ = packedPTG_out["hwp"]

    top_ra_out = np.concatenate((top_ra_, top_ra[i_tod]))
    top_dec_out = np.concatenate((top_dec_, top_dec[i_tod]))
    top_pa_out = np.concatenate((top_pa_, top_pa[i_tod]))
    bot_ra_out = np.concatenate((bot_ra_, bot_ra[i_tod]))
    bot_dec_out = np.concatenate((bot_dec_, bot_dec[i_tod]))
    bot_pa_out = np.concatenate((bot_pa_, bot_pa[i_tod]))
    hwpang_out = np.concatenate((hwpang_,hwpang))
    
    packedPTG_out_ = {"top_ra":top_ra_out, "top_dec":top_dec_out, "top_pa":top_pa_out,
                      "bot_ra":bot_ra_out, "bot_dec":bot_dec_out, "bot_pa":bot_pa_out, "hwp":hwpang_out}
    return packedPTG_out_
#######################################################################################################
#######################################################################################################
def Read_SimMap_part(filename, theta_c, phi_c, width, nside_out, option_TQU, option_silent):
    if option_silent=='Y': option_silent=True
    if option_silent=='N': option_silent=False
    if option_silent==False: print "[READ SimMap_part]: START reading "+filename
    map = h.read_map(filename,field=0)
    npix = len(map)
    nside = h.npix2nside(npix)
    ipix = np.arange(npix,dtype='int')
    theta, phi = h.pix2ang(nside,ipix)
    if ( ((phi_c < pi) & (phi_c > width))         
         or ((phi_c>=pi) & (abs(2.*pi-phi_c)>width) )):
        pass
    if ( (phi_c < pi) & (phi_c < width) ):
        ind = np.where(phi>pi)
        phi[ind[0]] = phi[ind[0]]-2.*pi
    if ( (phi_c > pi) & (abs(2.*pi-phi_c) < width) ):
        ind = np.where(phi<pi)
        phi[ind[0]] = phi[ind[0]]+2.*pi
    ind_part_in = np.where(((theta>theta_c-width/2.) & (theta<theta_c+width/2.))
                   & ((phi>phi_c-width/2.) & (phi<phi_c+width/2.)))

    theta=0; phi=0; ipix=0
    mapout_T = map[ind_part_in[0]];    map = 0

    if option_TQU=='TQU':
        map = h.read_map(filename,field=1)
        mapout_Q = map[ind_part_in[0]];    map = 0
        map = h.read_map(filename,field=2)
        mapout_U = map[ind_part_in[0]];    map = 0

    ind_all = np.ones(npix,dtype='int')*(-1)
    ind_all[ind_part_in[0]] = np.arange(len(ind_part_in[0]),dtype='int')
    ind_band_in = ind_all[min(ind_part_in[0]):max(ind_part_in[0])+1]
    print '[Read_SimMap_part] part map, nside, ind_all', len(mapout_T), nside, len(ind_all)

    idx_min_out = min(ind_part_in[0])
    if nside != nside_out:
        print '[Read_SimMap_part]', nside, nside_out
        npix_out = h.nside2npix(nside_out)
        ipix = np.arange(npix_out,dtype='int')
        theta, phi = h.pix2ang(nside_out,ipix)
        if ( ((phi_c < pi) & (phi_c > width))         
             or ((phi_c>=pi) & (abs(2.*pi-phi_c)>width) )):
            pass
        if ( (phi_c < pi) & (phi_c < width) ):
            ind = np.where(phi>pi)
            phi[ind[0]] = phi[ind[0]]-2.*pi
        if ( (phi_c > pi) & (abs(2.*pi-phi_c) < width) ):
            ind = np.where(phi<pi)
            phi[ind[0]] = phi[ind[0]]+2.*pi
        ind_part_out = np.where(((theta>theta_c-width/2.) & (theta<theta_c+width/2.))
                                & ((phi>phi_c-width/2.) & (phi<phi_c+width/2.)))
        theta=0; phi=0; ipix=0
        idx_min_out = min(ind_part_out[0])

        ind_all_out = np.ones(npix_out,dtype='int')*(-1)
        ind_all_out[ind_part_out[0]] = np.arange(len(ind_part_out[0]),dtype='int')
        ind_band_out = ind_all_out[min(ind_part_out[0]):max(ind_part_out[0])+1]
    else:
        ind_band_out = ind_band_in
        ind_part_out = ind_part_in

    if option_TQU=='TQU':
        mapout = {'T':mapout_T, 'Q': mapout_Q, 'U': mapout_U, 'nside': nside, \
                      'ind_band_in': np.int_(ind_band_in), 'ind_max_in':max(ind_part_in[0]), 'ind_min_in':min(ind_part_in[0]), \
                      'ind_band_out':ind_band_out, 'ind_min_out':idx_min_out}
    if option_TQU=='T':
        mapout = {'T':mapout_T, 'nside': nside, \
                      'ind_band_in': np.int_(ind_band_in), 'ind_max_in':max(ind_part_in[0]), 'ind_min_in':min(ind_part_in[0]), \
                      'ind_band_out':ind_band_out, 'ind_min_out':idx_min_out}

    if option_silent==False: print "[READ SimMap_part]: END reading "+filename
    return mapout, ind_part_out, ind_band_out, idx_min_out
#######################################################################################################
#######################################################################################################
def Read_SimMap_part_old(filename, theta_c, phi_c, width, option_silent):
    if option_silent=='Y': option_silent=True
    if option_silent=='N': option_silent=False
    if option_silent==False: print "[READ SimMap_part]: START reading "+filename 
    map = h.read_map(filename,field=0)
    npix = len(map)
    nside = h.npix2nside(npix)
    ipix = np.arange(npix,dtype='int')
    theta, phi = h.pix2ang(nside,ipix)

#    print theta_c, phi_c, width, pi, 2.*pi-phi_c
#    if (phi_c > pi): print phi_c*radeg

    if ( ((phi_c < pi) & (phi_c > width)) 
         or ((phi_c>=pi) & (abs(2.*pi-phi_c)>width) )):
        pass
    if ( (phi_c < pi) & (phi_c < width) ):
        ind = np.where(phi>pi)
        phi[ind[0]] = phi[ind[0]]-2.*pi
    if ( (phi_c > pi) & (abs(2.*pi-phi_c) < width) ):
        ind = np.where(phi<pi)
        phi[ind[0]] = phi[ind[0]]+2.*pi

    ind_part = np.where(((theta>theta_c-width/2.) & (theta<theta_c+width/2.)) 
                   & ((phi>phi_c-width/2.) & (phi<phi_c+width/2.)))
#    print 'ind_part', len(ind_part[0]) 
    theta=0; phi=0; ipix=0

    mapout_T = map[ind_part[0]];    map = 0
    map = h.read_map(filename,field=1)
    mapout_Q = map[ind_part[0]];    map = 0
    map = h.read_map(filename,field=2)
    mapout_U = map[ind_part[0]];    map = 0

    ind_all = np.ones(npix,dtype='int')*(-1)
    ind_all[ind_part[0]] = np.arange(len(ind_part[0]),dtype='int')

#    print len(mapout_T), len(mapout_Q), len(mapout_U)
#    print len(ind_all), len(ind)
#    sys.exit()
    print '[Read_SimMap_part] part map, nside, ind_all', len(mapout_T), nside, len(ind_all)
    mapout = {'T':mapout_T, 'Q': mapout_Q, 'U': mapout_U, 'nside': nside, 'ind_all': np.int_(ind_all)}
    if option_silent==False: print "[READ SimMap_part]: END reading "+filename 
    return mapout, ind_part[0]
#######################################################################################################
#######################################################################################################
def maps_init_part(nside,npix,ind,ind_min):
    print '[maps_init_part] npix, ind', npix, len(ind)
    In = np.zeros(npix)
    Id = np.zeros(npix)
    AA = np.zeros(npix)
    BB = np.zeros(npix)
    AB = np.zeros(npix)
    Ad = np.zeros(npix)
    Bd = np.zeros(npix)
    H = np.zeros(npix)
#    print '[maps_init_part] creating map_pack' 
    map_pack = {"In":In, "Id":Id, "AA":AA, "BB":BB, "AB":AB, "Ad": Ad, "Bd":Bd, "H":H, "ind_band_out":ind, "nside":nside, "ind_min":ind_min}
#    print '[maps_init_part] just before return' 
    return map_pack
#######################################################################################################
#######################################################################################################
def maps_zeros_part(packedMap):
#    print '[maps_free_part] start'
    npix = len(packedMap['In'])
#    print '[maps_free_part] start1'
    In = np.zeros(npix)
    Id = np.zeros(npix)
    AA = np.zeros(npix)
    BB = np.zeros(npix)
    AB = np.zeros(npix)
    Ad = np.zeros(npix)
    Bd = np.zeros(npix)
    H = np.zeros(npix)
#######################################################################################################
#######################################################################################################
def Cal_ConditionNum_QUpix(AA,BB,AB):
    m = np.zeros((2,2))
    m[0][0] = AA
    m[0][1] = AB
    m[1][0] = m[0][1]
    m[1][1] = BB

    if ( (np.isnan(AA)) or (np.isnan(BB)) or (np.isnan(AB)) ):
        print '[Cal_ConditionNum_QUpix] AA, BB, AB', np.isnan(AA), np.isnan(BB), np.isnan(AB)
#    print m
    w, v = np.linalg.eig(m)
    condition_number = np.max(w)/np.min(w)
    return condition_number
#######################################################################################################
#######################################################################################################
def MapMake_part(i_pix, nside, packedTOD, packedPTG, packedMap, idx_min_out, time0, option_silent, option_pixelmapio):
#def MapMake_part(i_pix, nside, packedTOD, packedPTG, ind_all, time0, option_silent):
    if option_silent==False: print "=================================="
    if option_silent==False: info("[MAPMAKE] START..., pixel#=" + str(i_pix), time0)
#    info("[MAPMAKE] START..., pixel#=" + str(i_pix), time0)

    pi = np.pi
    
    top_ra = packedPTG["top_ra"]
    top_dec = packedPTG["top_dec"]
    top_psi = packedPTG["top_pa"]

    bot_ra = packedPTG["bot_ra"]
    bot_dec = packedPTG["bot_dec"]
    bot_psi = packedPTG["bot_pa"]

    hwpang = packedPTG["hwp"]
    
    sum_dat = packedTOD["sum_tod"]
    dif_dat = packedTOD["dif_tod"]

    sum_std = packedTOD["sum_std"]
    dif_std = packedTOD["dif_std"]

    p_weight = 1./sum_std**2.
    n_weight = 1./dif_std**2.

    if option_silent==False: info("[MAPMAKE] mid1, pixel#=" + str(i_pix), time0)
#    info("[MAPMAKE] mid1, pixel#=" + str(i_pix), time0)

    top_gamma = 1. # gamma[2*i_pix]
    bot_gamma = 1. # gamma[2*i_pix+1]

    ra = 0.5*(top_ra + bot_ra)
    dec = 0.5*(top_dec + bot_dec)
    alpha = 0.5*(top_gamma*np.cos(4.*hwpang-2.*top_psi) - bot_gamma*np.cos(4.*hwpang-2.*bot_psi))
    beta = 0.5*(top_gamma*np.sin(4.*hwpang-2.*top_psi) - bot_gamma*np.sin(4.*hwpang-2.*bot_psi))

    i_hpx = h.ang2pix(nside, pi/2.-dec, ra)
    nb = len(sum_dat)

    if option_silent==False: info("[MAPMAKE] mid1 before forloop, pixel#=" + str(i_pix) + ', nb'+str(nb), time0)
    nb_map = len(packedMap['Id'])
    ind_band_out = packedMap['ind_band_out']
#    MapProjection_forLoopInC_part_inline(i_hpx,sum_dat,dif_dat,alpha,beta,p_weight,n_weight,packedMap) 
    out =  lib_c.MapProjection_forLoopInC_part(int(nb_map),
                                               np.float_(i_hpx-idx_min_out),
                                               np.float_(sum_dat),
                                               np.float_(dif_dat),
                                               np.float_(alpha),
                                               np.float_(beta),
                                               np.float_(p_weight),
                                               np.float_(n_weight),
                                               np.float_(ind_band_out)) 
    
    if option_pixelmapio==True:
        packedMap["In"] = out[0]
        packedMap["Id"] = out[1]
        packedMap["AA"] = out[2]
        packedMap["BB"] = out[3]
        packedMap["AB"] = out[4]
        packedMap["Ad"] = out[5]
        packedMap["Bd"] = out[6]
        packedMap["H"] = out[7]

    if option_pixelmapio==False:
        packedMap["In"] += out[0]
        packedMap["Id"] += out[1]
        packedMap["AA"] += out[2]
        packedMap["BB"] += out[3]
        packedMap["AB"] += out[4]
        packedMap["Ad"] += out[5]
        packedMap["Bd"] += out[6]
        packedMap["H"] += out[7]

    i_hpx=0;
    top_ra = 0; top_dec = 0; top_psi = 0;
    bot_ra = 0; bot_dec = 0; bot_psi = 0;
    hwpang = 0;
    sum_dat=0; dif_dat=0
    alpha=0;   beta=0; ra=0; dec=0
    p_weight=0;  n_weight=0
    
    if option_silent==False: info("[MAPMAKE] one line before END, pixel#=" + str(i_pix), time0)
    if option_silent==False: info("[MAPMAKE] END, pixel#=" + str(i_pix), time0)
#    info("[MAPMAKE] END, pixel#=" + str(i_pix), time0)
    if option_silent==False: print "===================================================================================="
    return packedMap
#######################################################################################################
#######################################################################################################
#def main_mapmaker(SimInput,fpdb_list,NoiseInput,ptg_package,ptg_idx,pix_list,SysInputs,outdir_ptg,outdir):
def main_mapmaker(inputs):
    run_type = inputs["run_type"]
    runtime_init = inputs["runtime_init"]
    option_silent = inputs["silent"]

    if option_silent=='Y': option_silent=True
    if option_silent=='N': option_silent=False
    option_pixelmapio=inputs["pixelmapio"]
    if option_pixelmapio=='Y': option_pixelmapio=True
    if option_pixelmapio=='N': option_pixelmapio=False

#    print option_silent, inputs["silent"]
#    if option_silent==False: print "it is not silent"
    if option_silent==False: print ""
    if option_silent==False: print "===================================================================================="
    if option_silent==False: info("[MAIN_MAPMAKER] START...", runtime_init)
    run_id = inputs["run_id"]

    if run_type=='sim_realptg' or run_type=='sim':
        option_TQU = inputs["TQU"]
        SimInputs = inputs["SimInput"]
        idx_min_in = SimInputs['ind_min_in']
        idx_min_out = SimInputs['ind_min_out']
        fpdb_list_simgen = inputs["fpdb_list_simgen"]
        NoiseInput = inputs["NoiseInput"]
        RelgainNoiseInput = inputs["NoiseInput"]

    fpdb_list_mmin = inputs["fpdb_list_mmin"]
    MuellerMatrix = inputs["muellermatrix"]
    Relgain = inputs["relgain"]
    pix_list = inputs["pix_list"]
    ptg_package = inputs["ptg_package"]
    ptg_idx = inputs["ptg_idx"]
    out_dir = inputs["out_dir"]
    out_dir_ptg = inputs["out_dir_ptg"]
    poly_n = inputs["poly"]

    nsideout = inputs["nside"]
    npix = h.nside2npix(nsideout)
#    packedMap = maps_init(nsideout)
#    packedMap = maps_init_part(nsideout,len(SimInputs["T"]),SimInputs["ind_all"])
    top_ptg = {}
    bot_ptg = {}
    n_pix = len(pix_list)

    net = NoiseInput['net']
    knee = NoiseInput['fknee']
    power = NoiseInput['power']
    seed = NoiseInput['seed']

    hwpang = ptg_package['hwp']
    packedMap_arr = {}
    
    out_dir_arr = []
    pix_arr = []

    packedMap = maps_init_part(nsideout,len(SimInputs["T"]),SimInputs["ind_band_out"],SimInputs["ind_min_out"])
#    print 'end maps_init_part'
    if option_silent==False: info("[MAIN_MAPMAKER] before for loop...", runtime_init)
#    info("[MAIN_MAPMAKER] before for loop...", runtime_init)
    for i_pix in range(0,n_pix):
#        packedMap = maps_init_part(nsideout,len(SimInputs["T"]),SimInputs["ind_all"])
        if option_silent==False: print "[MAIN_MAPMAKER]: DETECTOR #", pix_list[i_pix]
#        info("[MAIN_MAPMAKER] START FOR LOOP i_pix="+str(i_pix)+"...", runtime_init)
        top_ptg, bot_ptg = boreptg2pixptg(i_pix,pix_list,ptg_package,fpdb_list_simgen,0,option_silent)
#        info("[MAIN_MAPMAKER] after 1st boreptg2pixptg...", runtime_init)
       
        s_idx = ptg_idx["s_idx"]; e_idx = ptg_idx["e_idx"];
        direc = ptg_idx["direc"]; scan_flag = ptg_idx["flag_hscan"];  
        
        hscan_idx = np.where(scan_flag == 0)

#        info("[MAIN_MAPMAKER] before signalgen...", runtime_init)
        top_tod, bot_tod = SignalGen(i_pix, pix_list, top_ptg, bot_ptg, hwpang, SimInputs, MuellerMatrix,0, option_TQU, option_silent)
#        info("[MAIN_MAPMAKER] after signalgen...", runtime_init)
        packedTODout = ConcatenateTOD_init()
        packedPTGout = ConcatenatePTG_init()
        top_noise = []
        bot_noise = []

#        if i_pix==30: sys.exit()
#        np.savez('/global/homes/t/tmatsumu/develop/PBI/repo_test/PB1_NTP_develop/gen_run/tmp_'+str(i_pix),\
#                     top_tod,bot_tod,top_ptg['ra'], bot_ptg['ra'], top_ptg['dec'], bot_ptg['dec'])
#        if i_pix==30: sys.exit()
#        info("[MAIN_MAPMAKER] before 2nd boreptg2pixptg...", runtime_init)
        top_ptg, bot_ptg = boreptg2pixptg(i_pix,pix_list,ptg_package,fpdb_list_mmin,0,option_silent)
#        info("[MAIN_MAPMAKER] after 2nd boreptg2pixptg and just before i_hscan for loop...", runtime_init)

        for i_hscan in hscan_idx[0]:
            seed = NoiseInput['seed']
            nbData = e_idx[i_hscan] - s_idx[i_hscan] + 1
            if (nbData%2 == 1): nbData = nbData-1 
            i_tod = s_idx[i_hscan] + np.array(range(0,nbData))
 
#            top_noise, bot_noise = NoiseGen_auto(nbData,net[i_pix],knee[i_pix],power[i_pix],seed[i_pix]*i_tmp*i_hscan*i_pix,0,option_silent)
            top_dat = top_tod[i_tod] #+ top_noise
            bot_dat = bot_tod[i_tod] #+ bot_noise

#            if (i_hscan > 0): break
#            top_gain, bot_gain = NoiseGen_auto(nbData,gain_net[i_pix],gain_knee[i_pix],gain_power[i_pix],gain_seed[i_pix]*i_tmp*i_hscan*i_pix+1e5)
#            packedTOD_ = SignalFilt(top_dat, bot_dat, top_gain, bot_gain, poly_n)
#            top_dat, bot_dat = RelgainCorr(top_dat, bot_dat, Relgain, time)
            packedTOD_ = SignalFilt(top_dat, bot_dat, poly_n, 0,option_silent)
            packedTODout = ConcatenateTOD(packedTODout,packedTOD_,0,option_silent)
            packedPTGout = ConcatenatePTG(packedPTGout,top_ptg,bot_ptg,hwpang[i_tod],i_tod,0,option_silent)
#        del(top_dat); del(bot_dat); del(packedTOD_)
#        info("[MAIN_MAPMAKER] end of the i_hscan for loop...", runtime_init)

        if option_silent==False: info('[MAIN_MAPMAKER: pre-MapMake func] runtime', runtime_init) #,option_silent)
        MapMake_part(i_pix, nsideout, packedTODout, packedPTGout, packedMap, idx_min_out, runtime_init, option_silent, option_pixelmapio)
        if option_silent==False: info('[MAIN_MAPMAKER: post-MapMake func] runtime', runtime_init) #,option_silent)
#        del(packedTODout); del(packedPTGout)

#        ind = np.where(np.abs(packedMap["Id"]) > np.max(packedMap["Id"])*1.e-3)
        tmpH = np.copy(packedMap["H"])
        ind_select = np.where( tmpH > 0 )
        tmpH = 0
        if option_silent==False: print "[MAIN_MAPMAKER: OUTPUT DIRECTORY]: ", out_dir, str(i_pix)
#        info('[MAIN_MAPMAKER: pre-save npy obj] runtime', runtime_init)
        if option_pixelmapio==True:
            np.savez(out_dir+ '/map_'+str(i_pix),
                     ind_select[0],
                     packedMap["In"][ind_select[0]],
                     packedMap["Id"][ind_select[0]],
                     packedMap["AA"][ind_select[0]],
                     packedMap["BB"][ind_select[0]],
                     packedMap["AB"][ind_select[0]],
                     packedMap["Ad"][ind_select[0]],
                     packedMap["Bd"][ind_select[0]],
                     packedMap["H"][ind_select[0]] )
            if i_pix != n_pix-1: maps_zeros_part(packedMap)
        if option_silent==False: info('[MAIN_MAPMAKER: post-save npy obj] runtime', runtime_init)

    if option_pixelmapio==False:
        tmpH = np.copy(packedMap["H"])
        ind_select = np.where( tmpH > 0 )
        tmpH = 0
        np.savez(out_dir+ '/map_all',
                 ind_select[0],
                 packedMap["In"][ind_select[0]],
                 packedMap["Id"][ind_select[0]],
                 packedMap["AA"][ind_select[0]],
                 packedMap["BB"][ind_select[0]],
                 packedMap["AB"][ind_select[0]],
                 packedMap["Ad"][ind_select[0]],
                 packedMap["Bd"][ind_select[0]],
                 packedMap["H"][ind_select[0]] )

#    info('[MAIN_MAPMAKER: post-save npy obj] runtime', runtime_init)
    if option_silent==False: print "[main_mapmaker.py]: End of the for loop"
#    print "[main_mapmaker.py]: End of the for loop, ", i_pix, '/', n_pix

    if option_silent==False: print "[MAIN_MAPMAKER]: generating the mapfilelist.db"
#    print "[MAIN_MAPMAKER]: generating the mapfilelist.db"
    db_name = out_dir+'/mapfilelist.db'
    os.system("rm -f "+db_name)
    conn = sq.connect(db_name)
    c = conn.cursor()
    c.execute('create table mapfilelist (run_id integer, outdir text, pix integer)')
#    print "[MAIN_MAPMAKER]: generating the mapfilelist.db, 2"
    if option_pixelmapio==True: 
        for i_pix in range(0,n_pix):
            list_entries = (run_id, out_dir, pix_list[i_pix][2])
            c.execute('insert into mapfilelist values (?,?,?)', list_entries)
    if option_pixelmapio==False: 
        list_entries = (run_id, out_dir, -1)
        c.execute('insert into mapfilelist values (?,?,?)', list_entries)
#    print "[MAIN_MAPMAKER]: generating the mapfilelist.db, 3"
    conn.commit()
    c.close()  
#    print "[MAIN_MAPMAKER]: generating the mapfilelist.db, 4"

    if option_silent==False: info("[MAIN_MAPMAKER] END for a scanset...", runtime_init)
#    info("[MAIN_MAPMAKER] END for a scanset...", runtime_init)
#    print "===================================================================================="

#######################################################################################################
#######################################################################################################

if __name__ == '__main__': _main()
