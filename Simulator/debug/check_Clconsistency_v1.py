import numpy as np
import pylab as py
import healpy as h
import os
import sys

pi = np.pi 
radeg = (180./pi)

'''
   check if the cl from the map and the original cl are consistent to each other.
   2014-2-3, T. Matsumura
'''


#+++++++++
# READ MAP

dir_map = '/group/cmb/litebird/simdata/Maps/CMB/tmp/'
fname_map = 'cmb_map_band40.fits'
map = h.read_map(dir_map+'/'+fname_map,field=0)
npix = len(map)
nside = h.npix2nside(npix)
print 'nside = ', nside

#+++++++++
# COMPUTE THE ORIGINAL Cl

fname_cl = 'tmp/Cl_tmp1.fits'
fname_anafast = 'tmp/tmp_cl1'
#os.system('rm -rf '+fname_cl)
f_ana = open(fname_anafast+'.param', 'w')
f_ana.write('infile = %s\n' % (dir_map+'/'+fname_map))
f_ana.write('outfile = %s\n' % (fname_cl))
f_ana.write('simul_type = 2 \n')
f_ana.write('nlmax = %d \n' % (3*nside-1))
f_ana.close()

os.system('bsub -q e -o '+fname_anafast+'.o'+' anafast -d '+fname_anafast+'.param')
print fname_cl
cl = h.read_cl(fname_cl)
ell = np.arange(len(cl[0]))
const = ell*(ell+1.)/2./pi

py.plot(ell,const*cl[0],'-b')
py.plot(ell,const*cl[1],'--b')
py.plot(ell,const*cl[2],'-.b')
#py.plot(ell,const*cl[3],'.')

#+++++++++
# READ THE ORIGINAL MAP AND SMOOTH

fwhm_arcmin=30.
fwhm_in=fwhm_arcmin/60./radeg
map_Tin = h.read_map(dir_map+'/'+fname_map,field=0)
map_Qin = h.read_map(dir_map+'/'+fname_map,field=1)
map_Uin = h.read_map(dir_map+'/'+fname_map,field=2)
map_in = np.array([map_Tin,map_Qin,map_Uin])
map_out = h.smoothing(map_in, fwhm=fwhm_in)
fname_map_smooth = dir_map+'/cmb_map_band40_b'+str(fwhm_arcmin)+'arcmin.fits'
h.write_map(fname_map_smooth,map_out)

#+++++++++
# COMPUTE THE SMOOTHED Cl

fname_cl = 'tmp/Cl_tmp2.fits'
#os.system('rm -rf '+fname_cl)
fname_anafast = 'tmp/tmp_cl2'
f_ana = open(fname_anafast+'.param', 'w')
f_ana.write('infile = %s\n' % (fname_map_smooth))
f_ana.write('outfile = %s\n' % (fname_cl))
f_ana.write('simul_type = 2 \n')
f_ana.write('nlmax = %d \n' % (3*nside-1))
f_ana.close()

os.system('bsub -q e -o '+fname_anafast+'.o'+' anafast -d '+fname_anafast+'.param')
print fname_cl
cl = h.read_cl(fname_cl)
ell = np.arange(len(cl[0]))

py.plot(ell,const*cl[0],'-m')
py.plot(ell,const*cl[1],'--m')
py.plot(ell,const*cl[2],'-.m')
#py.plot(ell,const*cl[3],'.')

#+++++++++
# COMPUTE THE BEAM AND DIVITED THE SMOOTHED Cl BY THE BEAM

beam = h.gauss_beam(fwhm_in,lmax=nside)
pixwin = h.pixwin(nside)
beam = beam*pixwin[0:len(beam)]
print len(cl[0]), len(beam), len(pixwin)

py.plot(ell,const*cl[0]/beam**2,'.c')
py.plot(ell,const*cl[1]/beam**2,'.c')
py.plot(ell,const*cl[2]/beam**2,'.c')

#+++++++++
# UD_GRADE THE MAPS TO NSIDE_OUT

nside_out = 512
map_out_ud = h.ud_grade(map_out,nside_out=nside_out)
fname_map_ud = dir_map+'/cmb_map_band40_b'+str(fwhm_arcmin)+'arcmin_nside'+str(nside_out)+'.fits'
h.write_map(fname_map_ud,map_out_ud)

#+++++++++
# COMPUTE THE SMOOTHED Cl

fname_cl = 'tmp/Cl_tmp3.fits'
#os.system('rm -rf '+fname_cl)
fname_anafast = 'tmp/tmp_cl3'
print fname_map_ud
f_ana = open(fname_anafast+'.param', 'w')
f_ana.write('infile = %s\n' % (fname_map_ud))
f_ana.write('outfile = %s\n' % (fname_cl))
f_ana.write('simul_type = 2 \n')
f_ana.write('nlmax = %d \n' % (3*nside_out-1))
f_ana.close()

os.system('bsub -q e -o '+fname_anafast+'.o'+' anafast -d '+fname_anafast+'.param')
print fname_cl
cl = h.read_cl(fname_cl)
ell = np.arange(len(cl[0]))
const = ell*(ell+1.)/2./pi

beam = h.gauss_beam(fwhm_in,lmax=2*nside_out)
pixwin = h.pixwin(nside_out)
beam = beam*pixwin[0:len(beam)]
print len(cl[0]), len(beam), len(pixwin)

py.plot(ell,const*cl[0],'-r')
py.plot(ell,const*cl[1],'--r')
py.plot(ell,const*cl[2],'-.r')
#py.plot(ell,const*cl[3],'.')

py.plot(ell,const*cl[0]/beam**2,'.g')
py.plot(ell,const*cl[1]/beam**2,'.g')
py.plot(ell,const*cl[2]/beam**2,'.g')
#py.plot(ell,const*cl[3],'.')


py.loglog()
py.savefig('tmp/tmp.png')
os.system('display tmp/tmp.png &')
