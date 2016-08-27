import numpy as np
import pylab as py
import healpy as h
import os
import sys

pi = np.pi
radeg = (180./pi)

num_file = int(sys.argv[1])
filename_in = sys.argv[2]
filename0 = sys.argv[3]
filename1 = sys.argv[4]
filename2 = sys.argv[5]
filename3 = sys.argv[6]
filename4 = sys.argv[7]
if num_file == 5:
    filename_out = sys.argv[8]
    beam_FWHM_arcmin=float(sys.argv[9])
if num_file == 6:
    filename5 = sys.argv[8]
    filename_out = sys.argv[9]
    beam_FWHM_arcmin=float(sys.argv[10])

clin = h.read_cl(filename_in)
cl0 = h.read_cl(filename0)
cl1 = h.read_cl(filename1)
cl2 = h.read_cl(filename2)
cl3 = h.read_cl(filename3)
cl4 = h.read_cl(filename4)
if num_file == 6:
    cl5 = h.read_cl(filename5)

#+++++++++++++++++++++++++++++++++

nside = 512
pixwin = h.pixwin(nside)
B_l = h.gauss_beam(beam_FWHM_arcmin/60./radeg,lmax=nside)
B_l = B_l*pixwin[0:nside+1]

#+++++++++++++++++++++++++++++++++

num_ell = len(clin[0])
ell = np.arange(num_ell)
const = ell*(ell+1.)/2./pi

#+++++++++++++++++++++++++++++++++
py.figure(figsize=(12,8))
py.subplot(221)
py.plot(ell, const*clin[0]/B_l**2,'r-')
py.plot(ell, const*clin[1]/B_l**2,'g-')
py.plot(ell, const*clin[2]/B_l**2,'b-')

py.plot(ell, const*cl0[1]/B_l**2,'y.')
py.plot(ell, const*cl0[2]/B_l**2,'y.')

py.plot(ell, const*cl1[1]/B_l**2,'m.')
py.plot(ell, const*cl1[2]/B_l**2,'m.')

py.plot(ell, const*cl2[1]/B_l**2,'c.')
py.plot(ell, const*cl2[2]/B_l**2,'c.')

py.plot(ell, const*cl3[1]/B_l**2,'b.')
py.plot(ell, const*cl3[2]/B_l**2,'b.')

py.plot(ell, const*cl4[1]/B_l**2,'g.')
py.plot(ell, const*cl4[2]/B_l**2,'g.')

if num_file == 6:
    py.plot(ell, const*cl5[1]/B_l**2,'k.')
    py.plot(ell, const*cl5[2]/B_l**2,'k.')

py.loglog()
py.xlim([2,1000])
#py.xlabel('$l$')
py.ylabel('$l(l+1)C_l/2\pi$ [$\mu K^2$]')
py.title('EE,BB')
#+++++++++++++++++++++++++++++++++
py.subplot(223)
py.plot(ell, const*clin[3]/B_l**2,'b-')
py.plot(ell, const*cl0[3]/B_l**2,'y.') 
py.plot(ell, const*cl1[3]/B_l**2,'m.')
py.plot(ell, const*cl2[3]/B_l**2,'c.')
py.plot(ell, const*cl3[3]/B_l**2,'b.')
py.plot(ell, const*cl4[3]/B_l**2,'g.')
if num_file == 6:
    py.plot(ell, const*cl5[3]/B_l**2,'k.')
py.title('TE')
py.semilogx()
py.xlim([2,1000])
py.xlabel('$l$')
py.ylabel('$l(l+1)C_l/2\pi$ [$\mu K^2$]')

#+++++++++++++++++++++++++++++++++
py.subplot(224)
py.plot(ell, const*clin[4]/B_l**2,'b-')

py.plot(ell, const*cl0[4]/B_l**2,'y.')
py.plot(ell, const*cl1[4]/B_l**2,'m.')
py.plot(ell, const*cl2[4]/B_l**2,'c.')
py.plot(ell, const*cl3[4]/B_l**2,'b.')
py.plot(ell, const*cl4[4]/B_l**2,'g.')
if num_file == 6:
    py.plot(ell, const*cl5[4]/B_l**2,'k.')

py.title('TB')
py.semilogx()
py.xlim([2,1000])
py.xlabel('$l$')
#+++++++++++++++++++++++++++++++++

py.subplot(222)
py.plot(ell, const*clin[5]/B_l**2,'b-')
py.plot(ell, const*cl0[5]/B_l**2,'y.')
py.plot(ell, const*cl1[5]/B_l**2,'m.')
py.plot(ell, const*cl2[5]/B_l**2,'c.')
py.plot(ell, const*cl3[5]/B_l**2,'b.')
py.plot(ell, const*cl4[5]/B_l**2,'g.')
if num_file == 6:
    py.plot(ell, const*cl5[5]/B_l**2,'k.')

py.semilogx()
py.xlim([2,1000])
#py.xlabel('$l$')
py.title('EB')

py.savefig(filename_out+'.png')
#os.system('display '+filename_out+'.png')

#+++++++++++++++++++++++++++++++++
py.figure(figsize=(12,8))
py.subplot(121)
py.plot(ell, const*clin[0]/B_l**2,'r-')
py.plot(ell, const*clin[1]/B_l**2,'g-')
py.plot(ell, const*clin[2]/B_l**2,'b-')

py.plot(ell, const*cl0[1]/B_l**2,'y.')
py.plot(ell, const*cl0[2]/B_l**2,'y.')

py.plot(ell, const*cl1[1]/B_l**2,'m.')
py.plot(ell, const*cl1[2]/B_l**2,'m.')

py.plot(ell, const*cl2[1]/B_l**2,'c.')
py.plot(ell, const*cl2[2]/B_l**2,'c.')

py.plot(ell, const*cl3[1]/B_l**2,'b.')
py.plot(ell, const*cl3[2]/B_l**2,'b.')

py.plot(ell, const*cl4[1]/B_l**2,'g.')
py.plot(ell, const*cl4[2]/B_l**2,'g.')

if num_file == 6:
    py.plot(ell, const*cl5[1]/B_l**2,'k.')
    py.plot(ell, const*cl5[2]/B_l**2,'k.')

py.loglog()
py.xlim([2,1000])
#py.xlabel('$l$')
py.ylabel('$l(l+1)C_l/2\pi$ [$\mu K^2$]')
py.title('EE,BB')
py.savefig(filename_out+'_EEBB.png')
