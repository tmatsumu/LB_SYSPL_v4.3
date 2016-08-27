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
from multiprocessing import Pool
from multiprocessing import Process
import copy
import lib_mapmaker as lmm
import time
import pickle
import gen_qsub_PsuedoCl as gqp

pi = np.pi

xml_filename = sys.argv[1]
xml_input = lmm.Gen_scansetdirs(xml_filename)
walltime= sys.argv[2]

dir_simulator = xml_input["dir_simulator"]
runlog = xml_input["xml_file"]
nside = xml_input["nside"]

#os.popen("mkdir "+dir_simulator+"/RunJob/"+runlog)

f_qsub = open(dir_simulator+"/RunJob/"+runlog+"/qsub_combine_"+runlog+".sh", "w")
f_qsub.write('%s\n' % ('#PBS -S /bin/bash'))
f_qsub.write('%s\n' % ('#PBS -l nodes=1:ppn=1'))
f_qsub.write('%s\n' % ('#PBS -l walltime='+walltime))
f_qsub.write('%s\n' % ('#PBS -N '+runlog+'_combine'))
f_qsub.write('%s\n' % ('#PBS -A mp107'))
f_qsub.write('%s\n' % ('#PBS -q regular'))
f_qsub.write('%s\n' % ('#PBS -o '+runlog+'_combine.log'))
f_qsub.write('%s\n' % ('#PBS -V'))
f_qsub.write('%s\n' % ('#PBS -l gres=project:1'))
f_qsub.write('%s\n' % ('cd $PBS_O_WORKDIR'))
# main_mapcombiner.py
f_qsub.write('%s\n' % ('module load cmb'))
f_qsub.write('%s\n' % ('python '+dir_simulator+'/src/main_mapcombiner_part.py '+xml_filename+' ' \
                       +xml_input["dir_simedmap"]+'/mapT '\
                       +xml_input["dir_simedmap"]+'/mapQ '\
                       +xml_input["dir_simedmap"]+'/mapU '\
                       +xml_input["dir_simedmap"]+'/mask '\
                       +xml_input["dir_simedmap"]+'/Tweight '\
                       +xml_input["dir_simedmap"]+'/conditionN'))

f_qsub.write('%s\n' % ('python '+dir_simulator+'/viewer/viewmap.py '+xml_input["dir_simedmap"]+'/ '+'200e-6'+' 20e-6'))

# fits_sep2oneTQU.py
f_qsub.write('%s\n' % ('python '+dir_simulator+'/src/fits_sep2oneTQU.py '  \
                                               +xml_input["dir_simedmap"]+'/mapT.fits ' \
                                               +xml_input["dir_simedmap"]+'/mapQ.fits ' \
                                               +xml_input["dir_simedmap"]+'/mapU.fits ' \
                                               +xml_input["dir_simedmap"]+'/mapTQU.fits'))
#for i_scanset in range(0,nb_scanset):
dir_xpure = "/project/projectdirs/polar/user/tmatsumu/sim/Pseudo_Cl/xpure/data"
f_qsub.write('%s\n' % ("mkdir "+dir_xpure+"/"+runlog))
f_qsub.write('%s\n' % ("cp "+xml_input["dir_simedmap"]+"/*.fits "+dir_xpure+"/"+runlog))
f_qsub.close()

filename_xpurewindow = dir_simulator+"/RunJob/"+runlog+"/PseudoCl_xpurewindow_"+runlog+".sh"
fname_nhits = dir_xpure+'/'+runlog+'/Tweight_CNexcluded.fits'
fname_binmask = dir_xpure+'/'+runlog+'/mask_CNexcluded.fits'
gqp.gen_qsub_xpure_window(filename_xpurewindow, runlog, fname_nhits, fname_binmask, int(nside), dir_xpure)


f_sh = open(dir_simulator+"/RunJob/"+runlog+"/PseudoCl_xpure_"+runlog+".sh", "w")
# genPBS_xpure_in.sh
f_sh.write('%s\n' % ('#!/bin/sh'))
f_sh.write('%s\n' % ('/global/homes/t/tmatsumu/develop/PBI/Pseudo_Cl/xpure/RunJob/Xpure/genPBS_xpure_in.sh '+'Simed '+runlog+' '+str(nside)+' '+str(int(nside)*2)))
f_sh.close()
os.popen("chmod 744 "+dir_simulator+"/RunJob/"+runlog+"/PseudoCl_xpure_"+runlog+".sh")


sys.exit()
