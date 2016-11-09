import numpy as np
import pylab as py
import healpy as h
import os
import sys

pi = np.pi
radeg = (180./pi)

def gen_GaussianBl(ell,beam_FWHM_arcmin):
    sigma = beam_FWHM_arcmin/60./radeg/np.sqrt(8.*np.sqrt(2.))
    B_l = np.exp(-ell*(ell+1.)*sigma**2)
    return B_l

const = 1.

beam_FWHM_arcmin = float(sys.argv[3])

dir_in = sys.argv[1]
fname_Clin = sys.argv[2]
Clin = h.read_cl(dir_in+'/'+fname_Clin+'.fits')

num = len(Clin)
ell = np.arange(num)

py.subplot(121)

py.plot(ell, ell*(ell+1.)/2./np.pi*Clin*const, '-r')

py.loglog()
py.xlim([2,num*2])
py.ylim([1e-5,1e5])
py.legend(loc='best',prop={'size':8})
py.xlabel('$l$')
py.ylabel('$l(l+1)/2\pi C_l$ [$\mu K^2$]')

fname_png = 'Cl_anafast.png'
py.savefig(dir_in+'/'+fname_Clin+'.png')
