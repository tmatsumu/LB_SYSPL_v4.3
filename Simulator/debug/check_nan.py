import numpy as np
import pylab as py
import healpy as h
import os

pi = np.pi
radeg = (180./pi)

dir='/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v3.1/Simulator/data_proj/SimedMaps/RunLog/20140206_run40/'

#out=np.load(dir+'day0/tod/tod0.npz')
out=np.load(dir+'day26/tod/tod82.npz')

top_tod = out['top_tod']
bot_tod = out['bot_tod']
top_ptg = out['top_ptg_out'].item()
bot_ptg = out['bot_ptg_out'].item()

num = len(top_ptg['psi'])
num_step = 6
print num, num/num_step
print ''
py.figure(figsize=(15,8))
for i in range(0,num_step):
    py.subplot(num_step,1,i+1)
    print i, i*num/num_step, (i+1)*num/num_step-1, num/num_step
    py.plot(np.arange(num/num_step)+i*num/num_step,top_ptg['psi'][i*num/num_step:(i+1)*num/num_step]*radeg)
    py.plot(np.arange(num/num_step)+i*num/num_step,bot_ptg['psi'][i*num/num_step:(i+1)*num/num_step]*radeg)
    py.xlim([i*num/num_step-5000,(i+1)*num/num_step+5000])
    py.grid()
py.savefig('tmp/tmp1.png')
os.system('display tmp/tmp1.png &')

print ''
dif_psi = top_ptg['psi']-bot_ptg['psi']

py.figure(figsize=(15,8))
for i in range(0,num_step):
    py.subplot(num_step,1,i+1)
    print i, i*num/num_step, (i+1)*num/num_step-1
    py.plot(np.arange(num/num_step)+i*num/num_step,dif_psi[i*num/num_step:(i+1)*num/num_step]*radeg)
    py.xlim([i*num/num_step-5000,(i+1)*num/num_step+5000])
    py.grid()
py.savefig('tmp/tmp2.png')
os.system('display tmp/tmp2.png &')


#+++++++
#num_i = 863000
#num_d =   1000
num_i = 335000
num_d =   1000
py.subplot(2,1,1)
py.plot(top_ptg['psi'][num_i:num_i+num_d]*radeg)
py.plot(bot_ptg['psi'][num_i:num_i+num_d]*radeg)
py.plot(top_ptg['psi'][num_i:num_i+num_d]*radeg,'.')
py.plot(bot_ptg['psi'][num_i:num_i+num_d]*radeg,'.')
py.xlim([-100,num_d+100])
py.grid()

py.subplot(2,1,2)
py.plot(dif_psi[num_i:num_i+num_d]*radeg)
py.plot(dif_psi[num_i:num_i+num_d]*radeg,'.')
py.xlim([-100,num_d+100])
py.grid()
py.savefig('tmp/tmp3.png')
os.system('display tmp/tmp3.png &')

