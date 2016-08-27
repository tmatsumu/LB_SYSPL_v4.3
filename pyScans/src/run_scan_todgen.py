import numpy as np
import pylab as py
import healpy as h
import matsumulib as mylib
import lib_LBScan as liblb
import time
import os
import sys

import global_params as gl
sys.path.append(gl.path_inputparams)

import importlib
param_module=sys.argv[1]
g = importlib.import_module(param_module)
  
pi = np.pi
radeg = (180./pi)

runtime_i = time.time()

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
option_gen_ptgtxt = g.option_gen_ptgtxt

filename = title+'_scan_'+str(int(total_time/60.))+'min'+'_'+str(int(ydays))+'day_nside'+str(int(nside))+'_'+str(int(sample_rate))+'Hz'
os.system('mkdir -p dataout/'+filename)
dir_out = 'dataout/'+filename
 
today_julian_i = today_julian
  
npix = h.nside2npix(nside)
nhits = np.zeros(npix)
cos_r1 = np.zeros(npix)
sin_r1 = np.zeros(npix)
cos_r2 = np.zeros(npix)
sin_r2 = np.zeros(npix)
cos_r4 = np.zeros(npix)
sin_r4 = np.zeros(npix)
  
for i in range(0,ydays):
    print 'current date: ', i, '/', int(ydays)
    liblb.gen_scan(theta_antisun, theta_boresight, freq_antisun, freq_boresight, total_time, today_julian, sample_rate, dir_out, filename, title, nside, runtime_i, option_gen_ptgtxt)
    today_julian += 1.
    
    nhits += h.read_map('dataout/'+filename+'/nhits_tmp.fits')
    cos_r1 += h.read_map('dataout/'+filename+'/cos_r1_tmp.fits')
    sin_r1 += h.read_map('dataout/'+filename+'/sin_r1_tmp.fits')
    cos_r2 += h.read_map('dataout/'+filename+'/cos_r2_tmp.fits')
    sin_r2 += h.read_map('dataout/'+filename+'/sin_r2_tmp.fits')
    cos_r4 += h.read_map('dataout/'+filename+'/cos_r4_tmp.fits')
    sin_r4 += h.read_map('dataout/'+filename+'/sin_r4_tmp.fits')

    if (i%20 == 0): h.write_map('dataout/'+filename+'/nhits_'+str(ydays)+'_tmp.fits',nhits)
    if (i%20 == 0): h.write_map('dataout/'+filename+'/cos_r1_'+str(ydays)+'_tmp.fits',cos_r1)
    if (i%20 == 0): h.write_map('dataout/'+filename+'/sin_r1_'+str(ydays)+'_tmp.fits',sin_r1)
    if (i%20 == 0): h.write_map('dataout/'+filename+'/cos_r2_'+str(ydays)+'_tmp.fits',cos_r2)
    if (i%20 == 0): h.write_map('dataout/'+filename+'/sin_r2_'+str(ydays)+'_tmp.fits',sin_r2)
    if (i%20 == 0): h.write_map('dataout/'+filename+'/cos_r4_'+str(ydays)+'_tmp.fits',cos_r4)
    if (i%20 == 0): h.write_map('dataout/'+filename+'/sin_r4_'+str(ydays)+'_tmp.fits',sin_r4)

r1 = np.sqrt((cos_r1/np.float_(nhits))**2. + (sin_r1/np.float_(nhits))**2.)
r2 = np.sqrt((cos_r2/np.float_(nhits))**2. + (sin_r2/np.float_(nhits))**2.)
r4 = np.sqrt((cos_r4/np.float_(nhits))**2. + (sin_r4/np.float_(nhits))**2.)

filename_out = dir_out+'/'+filename+'_nhits'
h.write_map(filename_out+'.fits', nhits)
h.mollview(nhits, title='Nobs, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')


filename_out = dir_out+'/'+filename+'_cos_r1'
h.write_map(filename_out+'.fits', cos_r1)
h.mollview(cos_r1, title='cos_r1, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/'+filename+'_sin_r1'
h.write_map(filename_out+'.fits', sin_r1)
h.mollview(sin_r1, title='sin_r1, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/'+filename+'_r1'
h.write_map(filename_out+'.fits', r1)
h.mollview(r1, title='r1, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/'+filename+'_cos_r2'
h.write_map(filename_out+'.fits', cos_r2)
h.mollview(cos_r2, title='cos_r2, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/'+filename+'_sin_r2'
h.write_map(filename_out+'.fits', sin_r2)
h.mollview(sin_r2, title='sin_r2, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/'+filename+'_r2'
h.write_map(filename_out+'.fits', r2)
h.mollview(r2, title='r2, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')


filename_out = dir_out+'/'+filename+'_cos_r4'
h.write_map(filename_out+'.fits', cos_r4)
h.mollview(cos_r4, title='cos_r4, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/'+filename+'_sin_r4'
h.write_map(filename_out+'.fits', sin_r4)
h.mollview(sin_r4, title='sin_r4, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/'+filename+'_r4'
h.write_map(filename_out+'.fits', r4)
h.mollview(r4, title='r4, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

print '[run_scan_todgen.py] END of the program'
print ''
