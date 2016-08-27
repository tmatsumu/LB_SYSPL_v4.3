import numpy as np
import healpy as h
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
import lib_mapmaker as lmm
import time
import AnalysisBackend.misc.util as util

pi = np.pi


#######################################################################################################
#######################################################################################################
print ""
print ""
print ""
print "#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#"
print "####  POLARBEAR generate pointing in galactic coordinate                        ####"
print "####   written by T. Matsumura (2011-5-17)                                      ####"
print "####     use as, >python gen_ptg_cel2gal.py xml_input                           ####"
print "####     revision: 2011-6-1, switch from .fits to .npy                          ####"
print "####     revision: 2011-11-13, add the flag field                               ####"
print "#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#"
print ""
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
runtime_init = time.time()

xml_filename = sys.argv[1]
xml_input = lmm.Gen_scansetdirs(xml_filename)

dir_simulator = xml_input["dir_simulator"] 
runlog = xml_input["xml_file"]

ptg_dirs = xml_input["dirs_pointing"]
ptg_dirs_part = xml_input["dirs_pointing_part"]
nb_scanset = len(ptg_dirs)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#lmm.info('########### Parent Process ID ############', runtime_init)
print nb_scanset
for i_scanset in range(0,nb_scanset):
    print "i_scanset", i_scanset+1, '/', nb_scanset
    ptg_package = lmm.read_ptg(ptg_dirs[i_scanset]+'/pointing.fits')
    glon,glat = lmm.euler_astrolib(ptg_package['ra'],ptg_package['dec'],1,radian=True)

    fn = ptg_dirs[i_scanset]+'/pointing_g.fits'
    print 'writing the fits to ', fn
    MJD_init = 0.
    print 'MJD', MJD_init, 'it used to be', 55414.61274250192
    fsample = 100.
    headerdict = {'rate': str(fsample), 'startmjd': str(MJD_init), 'Coordinate': 'G'}
    data = {}
    data['ra'] = glon
    data['dec'] = glat
    data['pa'] = ptg_package['pa']
    data['hwp'] = ptg_package['hwp']
    data['flag'] = ptg_package['flag'] # added 2011-11-13
    print fn
    util.savetofits(fn, data, headerdict)
