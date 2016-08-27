import numpy as np
import pylab as py
import lib_gain as lib_g
import os
import matsumulib as mylib

nbData = 20000
fsample = 150.
net = 50.
fknee = 0.1
power = 2.
seed = 1
model = 1

tod = lib_g.NoiseGen_auto(nbData,fsample,net,fknee,power,seed,model)
time = np.arange(nbData)/float(nbData)/fsample

freq, psd = mylib.calPSD(tod,fsample,4)
if model == 1:
    ps_model = net**2 * fsample*  (1.+ (fknee/freq)**power)
if model == 2:
    ps_model = net**2 * fsample*  (fknee/freq)**power

py.subplot(221)
py.plot(time, tod)

py.subplot(222)
par = mylib.plot_hist(tod,30,init_auto=True,fit=True)
py.title('$\sigma$='+str(par[2]))

py.subplot(212)
py.plot(freq, psd/np.sqrt(2.))
py.plot(freq, np.sqrt(ps_model))

py.ylim([1e0,1e4])
py.loglog()
py.savefig('tmp.png')
os.system('display tmp.png')

