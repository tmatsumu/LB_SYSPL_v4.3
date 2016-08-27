import numpy as np
import pylab as py
import healpy as h
import matsumulib as mylib
import sys
import os
import glob

import importlib
param_module=sys.argv[1]
g = importlib.import_module(param_module)
  
pi = np.pi
radeg = (180./pi)

title = g.title
sample_rate = g.sample_rate
total_time = g.total_time
theta_antisun = g.theta_antisun
freq_antisun = g.freq_antisun
theta_boresight = g.theta_boresight
freq_boresight = g.freq_boresight
ydays = g.ydays
nside = g.nside
today_julian = g.today_julian

filename = title+'_scan_'+str(int(total_time/60.))+'min'+'_'+str(int(ydays))+'day_nside'+str(int(nside))+'_'+str(int(sample_rate))+'Hz'

path_fits = glob.glob('dataout/'+filename+'/*.fits')
nb_fits = len(path_fits)

path_txt = glob.glob('dataout/'+filename+'/*.txt')
nb_txt = len(path_txt)

i = 0
time, dec, ra, angle, beta = mylib.read_txt5f(path_txt[i])

py.subplot(311)
py.plot(time-time[0],dec*180./pi,'.',markersize=0.5)
py.xlabel('time [s]')
py.ylabel('DEC [degrees]')
py.xlim([0,max(time-time[0])])

py.subplot(312)
py.plot(time-time[0],ra*180./pi,'.',markersize=0.5)
py.xlabel('time [s]')
py.ylabel('RA [degrees]')
py.xlim([0,max(time-time[0])])

py.subplot(313)
py.plot(time-time[0],angle*180./pi,'.',markersize=0.5)
py.plot(time-time[0],beta*180./pi,'.',markersize=0.5)
py.xlabel('time [s]')
py.ylabel('angle [degrees]')
py.xlim([0,max(time-time[0])])
py.savefig('dataout/'+filename+'/tmp1')

py.subplot(111)
py.plot(ra*180./pi,dec*180./pi,'.',markersize=0.5)
py.xlabel('RA [degrees]')
py.ylabel('DEC [degrees]')
py.savefig('dataout/'+filename+'/tmp2')

map = h.read_map('dataout/'+filename+'/'+filename+'_nhits.fits')
h.mollview(map, coord='E') #, rot=[290,70])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig('dataout/'+filename+'/tmp3')

