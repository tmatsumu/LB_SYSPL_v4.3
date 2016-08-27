import numpy as np
import pylab as py
import matsumulib as mylib
import os

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
py.show()
