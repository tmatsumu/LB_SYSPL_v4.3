import numpy as np
import pylab as py
import healpy as h
import matsumulib as mylib
import os
import sys

pi = np.pi

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
        if i>=0:
            ar = line.split()
            ell.append(int(ar[0]))
            TT.append(float(ar[1]))
            TE.append(float(ar[2]))
            EE.append(float(ar[3]))
            BB.append(float(ar[4]))
        i+=1
    return np.array(ell), np.array(TT), np.array(TE), np.array(EE), np.array(BB)

dir = '/group/cmb/litebird/simdata/Maps/CMB/CMB/cosmo/standard/'
filename = 'standard_lensedtotCls.txt'
ell,TT,EE,BB,TE = read_cltxt(dir+'/'+filename)
py.plot(ell,TT)
py.plot(ell,TE)
py.plot(ell,EE)
py.plot(ell,BB)
py.loglog()
py.savefig('tmp.png')
os.system('display tmp.png &')

print ell

prefact = ell*(ell+1.)/(2.*pi)
Clin = [TT/prefact,EE/prefact,np.zeros(len(TT)),TE/prefact]
#h.write_cl(dir+'/standard_lensedtotCls_noBB.fits',Clin)
h.mwrfits(dir+'/standard_lensedtotCls_noBB.fits', Clin, hdu=1, colnames=['TT','EE','BB','TE'])

sys.exit()


