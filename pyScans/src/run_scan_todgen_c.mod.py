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
sys.path.append(gl.path_src)
 
pi = np.pi
radeg = (180./pi)

runtime_i = time.time()

option_sh = sys.argv[1]
if option_sh!='y':
    import importlib
    param_module=sys.argv[2]
    g = importlib.import_module(param_module)

    title = g.title
    sample_rate = g.sample_rate
    total_time = g.total_time
    theta_antisun = g.theta_antisun
    freq_antisun = 1./(g.freq_antisun*60.)
    theta_boresight = g.theta_boresight
    freq_boresight = g.freq_boresight
    ydays = g.ydays
    nside = g.nside
    today_julian = g.today_julian
    option_gen_ptg = g.option_gen_ptg
    dir_out = gl.path_out

if option_sh=='y':
    title = sys.argv[2]
    sample_rate = np.float(sys.argv[3])
    total_time = np.float(sys.argv[4])
    theta_antisun = np.float(sys.argv[5])/radeg
    freq_antisun = 1./(np.float(sys.argv[6])*60.)
    theta_boresight = np.float(sys.argv[7])/radeg
    freq_boresight = np.float(sys.argv[8])/60.
    ydays = np.int(sys.argv[9])
    nside = np.int(sys.argv[10])
    today_julian = gl.today_julian
    option_gen_ptg = gl.option_gen_ptg
    dir_out = gl.path_out

# run summary
print '            title:', title, type(title)
print '      smaple_rate:', sample_rate, type(sample_rate)
print '       total_time:', total_time, type(total_time)
print '    theta_antisun:', theta_antisun*radeg, type(theta_antisun)
print '     freq_antisun:', freq_antisun, type(freq_antisun)
print '  theta_boresight:', theta_boresight*radeg, type(theta_boresight)
print '   freq_boresight:', freq_boresight, type(freq_boresight)
print '            ydays:', ydays, type(ydays)
print '            nside:', nside, type(nside)
print '     today_julian:', today_julian, type(today_julian)
print '   option_gen_ptg:', option_gen_ptg, type(option_gen_ptg)
print '          dir_out:', dir_out, type(dir_out)
print '' 

filename = title+'_SCANSPEC_'+str(int(total_time/60.))+'min_'+str(theta_antisun*radeg)+'degs_'+str(1./freq_antisun/60.)+'min_'+str(theta_boresight*radeg)+'degs_'+str(freq_boresight*60.)+'rpm_'+str(int(ydays))+'day_nside'+str(int(nside))+'_'+str(int(sample_rate))+'Hz'
os.system('mkdir -p '+dir_out+filename)
os.system('mkdir -p '+dir_out+filename+'/fits')
os.system('mkdir -p '+dir_out+filename+'/png')
os.system('mkdir -p '+dir_out+filename+'/ps')
os.system('mkdir -p '+dir_out+filename+'/processed_data')
os.system('mkdir -p '+dir_out+filename+'/ptg')
dir_out = dir_out+filename

f = open(dir_out+"/ptg/runlog.txt", "w")
f.write('%s %s \n' % ('title, ', title))
f.write('%s %f \n' % ('smaple_rate [Hz], ', sample_rate))
f.write('%s %f \n' % ('total_time [sec] per day, ', total_time))
f.write('%s %f \n' % ('theta_antisun [degs], ', theta_antisun*radeg))
f.write('%s %f \n' % ('freq antisun [Hz], ', freq_antisun))
f.write('%s %f \n' % ('theta_boresight [degs], ', theta_boresight*radeg))
f.write('%s %f \n' % ('freq_boresight [Hz], ', freq_boresight))
f.write('%s %d \n' % ('ydays [days], ', ydays))
f.write('%s %d \n' % ('nside, ', nside))
f.write('%s %f \n' % ('today_julian (starting date), ', today_julian))
f.write('%s %s \n' % ('option_gen_ptg, ', option_gen_ptg))
f.write('%s %s \n' % ('dir_out, ', dir_out))
f.close() 

np.savez(dir_out+"/ptg/runlog",
         title=title,
         sample_rate=sample_rate,
         total_tme=total_time,
         theta_antisun=theta_antisun,
         freq_antisun=freq_antisun,
         theta_boresight=theta_boresight,
         freq_boresight=freq_boresight,
         ydays=ydays,
         nside=nside,
         today_julian=today_julian,
         option_gen_ptg=option_gen_ptg,
         dir_out=dir_out,
         help='keywords are \n title \n sample_rate \n total_time \n theta_antisun \n freq_antisun \n theta_boresight \n freq_boresight \n ydays \n nside \n today_julian \n option_gen_ptg \n dir_out')
 
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
    print 'current date: ', i+1, '/', int(ydays)
    mapout = liblb.gen_scan_c_mod(theta_antisun, theta_boresight, freq_antisun, freq_boresight, \
                                      total_time, today_julian, sample_rate, dir_out, filename, \
                                      title, nside, runtime_i, option_gen_ptg)

    today_julian += 1.

    nhits += mapout[0] 
    cos_r1 += mapout[1]
    sin_r1 += mapout[2]
    cos_r2 += mapout[3]
    sin_r2 += mapout[4]
    cos_r4 += mapout[5]
    sin_r4 += mapout[6]

    if ((i)%14 == 0): 
        h.write_map(dir_out+'/fits/nhits_'+str(i+1)+'_tmp.fits',nhits)
        h.write_map(dir_out+'/fits/cos_r1_'+str(i+1)+'_tmp.fits',cos_r1)
        h.write_map(dir_out+'/fits/sin_r1_'+str(i+1)+'_tmp.fits',sin_r1)
        h.write_map(dir_out+'/fits/cos_r2_'+str(i+1)+'_tmp.fits',cos_r2)
        h.write_map(dir_out+'/fits/sin_r2_'+str(i+1)+'_tmp.fits',sin_r2)
        h.write_map(dir_out+'/fits/cos_r4_'+str(i+1)+'_tmp.fits',cos_r4)
        h.write_map(dir_out+'/fits/sin_r4_'+str(i+1)+'_tmp.fits',sin_r4)

    if (i == ydays-1): 
        h.write_map(dir_out+'/fits/nhits_'+str(i+1)+'_tmp.fits',nhits)
        h.write_map(dir_out+'/fits/cos_r1_'+str(i+1)+'_tmp.fits',cos_r1)
        h.write_map(dir_out+'/fits/sin_r1_'+str(i+1)+'_tmp.fits',sin_r1)
        h.write_map(dir_out+'/fits/cos_r2_'+str(i+1)+'_tmp.fits',cos_r2)
        h.write_map(dir_out+'/fits/sin_r2_'+str(i+1)+'_tmp.fits',sin_r2)
        h.write_map(dir_out+'/fits/cos_r4_'+str(i+1)+'_tmp.fits',cos_r4)
        h.write_map(dir_out+'/fits/sin_r4_'+str(i+1)+'_tmp.fits',sin_r4)
        
r1 = np.sqrt((cos_r1/np.float_(nhits))**2. + (sin_r1/np.float_(nhits))**2.)
r2 = np.sqrt((cos_r2/np.float_(nhits))**2. + (sin_r2/np.float_(nhits))**2.)
r4 = np.sqrt((cos_r4/np.float_(nhits))**2. + (sin_r4/np.float_(nhits))**2.)

filename_out = dir_out+'/ps/'+filename+'_nhits'
h.write_map(filename_out+'.fits', nhits)
h.mollview(nhits, title='Nobs, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')


filename_out = dir_out+'/ps/'+filename+'_cos_r1'
h.write_map(filename_out+'.fits', cos_r1)
h.mollview(cos_r1, title='cos_r1, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/ps/'+filename+'_sin_r1'
h.write_map(filename_out+'.fits', sin_r1)
h.mollview(sin_r1, title='sin_r1, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/ps/'+filename+'_r1'
h.write_map(filename_out+'.fits', r1)
h.mollview(r1, title='r1, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/ps/'+filename+'_cos_r2'
h.write_map(filename_out+'.fits', cos_r2)
h.mollview(cos_r2, title='cos_r2, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/ps/'+filename+'_sin_r2'
h.write_map(filename_out+'.fits', sin_r2)
h.mollview(sin_r2, title='sin_r2, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/ps/'+filename+'_r2'
h.write_map(filename_out+'.fits', r2)
h.mollview(r2, title='r2, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')


filename_out = dir_out+'/ps/'+filename+'_cos_r4'
h.write_map(filename_out+'.fits', cos_r4)
h.mollview(cos_r4, title='cos_r4, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/ps/'+filename+'_sin_r4'
h.write_map(filename_out+'.fits', sin_r4)
h.mollview(sin_r4, title='sin_r4, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

filename_out = dir_out+'/ps/'+filename+'_r4'
h.write_map(filename_out+'.fits', r4)
h.mollview(r4, title='r4, '+filename_out, rot=[0.,0.])
h.graticule(dpar=10,dmer=10,coord='E')
py.savefig(filename_out+'.ps')

print '[run_scan_todgen_c.py] END of the program'
print ''
