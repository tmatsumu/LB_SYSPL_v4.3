import numpy as np
import pylab as py
import healpy as h
import ReadMapMakeXml as rxml
import os
import sys

pi = np.pi
radeg = (180./pi)

def read_cltxt(filename):
    import fileinput
    ell = []
    TT = []
    TE = []
    EE = []
    BB = []
    filelines = fileinput.input(filename)
    i=0
    for line in filelines:
        if i>0:
            ar = line.split()
            ell.append(int(ar[0]))
            TT.append(float(ar[1]))
            TE.append(float(ar[2]))
            EE.append(float(ar[3]))
            BB.append(float(ar[4]))
        i+=1
    return np.array(ell), np.array(TT), np.array(TE), np.array(EE), np.array(BB)

def gen_GaussianBl(ell,beam_FWHM_arcmin):
    sigma = beam_FWHM_arcmin/60./radeg/np.sqrt(8.*np.sqrt(2.))
    B_l = np.exp(-ell*(ell+1.)*sigma**2)
    return B_l

const = 1.

beam_FWHM_arcmin = float(sys.argv[2])

fname_Clin = sys.argv[1]
Clin = h.read_cl(fname_Clin)

num = len(Clin)
ell = np.arange(num)

py.subplot(121)

py.plot(ell, ell*(ell+1.)/2./np.pi*Clin*const, '-r')

py.loglog()
py.xlim([2,num*2])
py.ylim([1e-7,1e3])
py.legend(loc='best',prop={'size':8})
py.xlabel('$l$')
py.ylabel('$l(l+1)/2\pi C_l$ [$\mu K^2$]')

fname_png = 'Cl_anafast.png'
py.savefig(fname_png)
