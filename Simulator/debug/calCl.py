import numpy as np
import pylab as py
import healpy as h
import sys
import os

dir_here='/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v3.2/Simulator/debug/'

pi = np.pi
radeg = (180./pi)

filename = sys.argv[1]
beam_FWHM_arcmin = float(sys.argv[2])
name_out = sys.argv[3]

map = h.read_map(filename,field=0)
npix = len(map)
nside = h.npix2nside(npix)

print ''
print '==========='
print 'filename: ', filename
print 'beam_FHWM_arcmin: ', beam_FWHM_arcmin
print 'nside: ', nside

fname_anafast = dir_here+'tmp/'+name_out
os.system('rm -rf '+fname_anafast+'.params')
os.system('rm -rf '+fname_anafast+'.o')
os.system('rm -rf '+fname_anafast+'_Cl.fits')
os.system('rm -rf '+fname_anafast+'_Cl.png')

f_ana = open(fname_anafast+'.params', 'w')
f_ana.write('infile = %s\n' % (filename))
f_ana.write('outfile = %s\n' % (fname_anafast+'_Cl.fits'))
f_ana.write('simul_type = 2 \n')
f_ana.write('nlmax = %d \n' % (nside*2))
f_ana.close()

#os.system('bsub -q e -o '+fname_anafast+'.o'+' anafast -d '+fname_anafast+'.params')
os.system('anafast -d '+fname_anafast+'.params')
os.system('sleep 1m ')

Clin = h.read_cl(fname_anafast+'_Cl.fits')

pixwin = h.pixwin(nside)
B_l_tmp = h.gauss_beam(beam_FWHM_arcmin/60./radeg,lmax=nside*5)
B_l = B_l_tmp[0:len(Clin[0])]*pixwin[0:len(Clin[0])]

num = len(Clin[0])
ell = np.arange(num)
print 'lmax: ', num
print ''
print '==========='

py.subplot(121)
py.plot(ell, ell*(ell+1.)/2./np.pi*Clin[0]/B_l**2, label='TT')
py.plot(ell, ell*(ell+1.)/2./np.pi*Clin[1]/B_l**2, label='EE')
py.plot(ell, ell*(ell+1.)/2./np.pi*Clin[2]/B_l**2, label='BB')
py.loglog()
py.xlim([2,num*2])
py.ylim([min(ell[2:]*(ell[2:]+1.)/2./np.pi*Clin[2][2:]/B_l[2:]**2),max(ell[2:]*(ell[2:]+1.)/2./np.pi*Clin[0][2:]/B_l[2:]**2)])
#py.xlim([1,num*2])
#py.ylim([min(ell[1:]*(ell[1:]+1.)/2./np.pi*Clin[2][1:]/B_l[1:]**2),max(ell[1:]*(ell[1:]+1.)/2./np.pi*Clin[0][1:]/B_l[1:]**2)])
py.legend(loc='best',prop={'size':8})
py.xlabel('$l$')
py.ylabel('$l(l+1)/2\pi C_l$ [$\mu K^2$]')
py.title(name_out)
py.savefig(fname_anafast+'_Cl.png')
os.system('display '+fname_anafast+'_Cl.png'+' &')
