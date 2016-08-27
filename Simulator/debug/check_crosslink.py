import numpy as np
import pylab as py
import healpy as h
import glob
import sys
import os
import matsumulib as mylib

pi = np.pi
radeg = (180./pi)

dir_in = './'
fileNames = glob.glob(dir_in+"/tmp_*.npz")

ind = []
alpha = []; beta = []
dif_dat = []; sum_dat = []
ra = []; dec = []
top_psi = []; bot_psi = []
i_pix = []

for i_file in fileNames:
    out = np.load(i_file)

    ind = np.where( ((out['dec'] > 25.825/radeg) & (out['dec'] < 26./radeg)) \
                        & ((out['ra'] > 35.825/radeg) & (out['ra'] < 36./radeg)) )
    print i_file, len(ind[0])
    if len(ind[0])==0: continue
    ra = np.hstack((ra,out['ra'][ind]))
    dec = np.hstack((dec,out['dec'][ind]))
    alpha = np.hstack((alpha,out['alpha'][ind]))
    beta = np.hstack((beta,out['beta'][ind]))
    dif_dat = np.hstack((dif_dat,out['dif_dat'][ind]))
    sum_dat = np.hstack((sum_dat,out['sum_dat'][ind]))
    top_psi = np.hstack((top_psi,out['top_psi'][ind]))
    bot_psi = np.hstack((bot_psi,out['bot_psi'][ind]))

num = len(alpha)

top_tod = 10.*np.ones(num)
bot_tod = 10.*np.ones(num)*0.0001
sum_dat = 0.5*(top_tod+bot_tod)
dis_dat = 0.5*(top_tod-bot_tod)

top_psi = np.random.uniform(0.,pi,num)
bot_psi = top_psi + 90./radeg

alpha = 0.5*(np.cos(-2.*top_psi) - np.cos(-2.*bot_psi))
beta = 0.5*(np.sin(-2.*top_psi) - np.sin(-2.*bot_psi))

print ''
print 'num', num
print '<A>, <A^2>, <A>^2', np.sum(alpha), np.sum(alpha**2), np.sum(alpha)**2
print '<B>, <B^2>, <B>^2', np.sum(beta), np.sum(beta**2), np.sum(beta)**2
print '<AB>, <A><B>', np.sum(alpha*beta), np.sum(alpha)*np.sum(beta)
det = np.sum(alpha*alpha)*np.sum(beta*beta)+2.*np.sum(alpha)*np.sum(beta)*np.sum(alpha*beta) \
    - np.sum(alpha*alpha)*np.sum(beta)**2 - np.sum(alpha)**2*np.sum(beta*beta)-np.sum(alpha*beta)**2
dalpha = dif_dat*alpha
dbeta = dif_dat*beta
print '<dA>, <dB>', np.sum(dalpha), np.sum(dbeta)
print '<d>', np.sum(dif_dat)
print 'mean_difftod', np.mean(dif_dat)
print 'det', det

print ''
print  (np.sum(alpha*alpha)*np.sum(beta*beta)-np.sum(alpha*beta)**2)/det, (-np.sum(alpha)*np.sum(beta*beta)+np.sum(beta)*np.sum(beta*alpha))/det, (np.sum(alpha)*np.sum(beta*alpha)-np.sum(beta)*np.sum(alpha*alpha))/det
print (-np.sum(alpha)*np.sum(beta*beta)+np.sum(beta)*np.sum(beta*alpha))/det, (np.sum(beta*beta)-np.sum(beta)**2)/det, (-np.sum(beta*alpha)+np.sum(alpha)*np.sum(beta))/det
print (np.sum(alpha*beta)*np.sum(alpha) - np.sum(beta)*np.sum(alpha*alpha) )/det, (-np.sum(beta*alpha)+np.sum(alpha)*np.sum(beta))/det, (np.sum(alpha*alpha)-np.sum(alpha)**2)/det
print ''
T = 2.*np.sum(sum_dat)/np.float(len(sum_dat))
S = 2./det*( (np.sum(alpha*alpha)*np.sum(beta*beta)-np.sum(alpha*beta)**2)*np.sum(dif_dat) \
                 +(-np.sum(alpha)*np.sum(beta*beta)+np.sum(beta)*np.sum(beta*alpha))*np.sum(dalpha) \
                 +(np.sum(alpha)*np.sum(beta*alpha)-np.sum(beta)*np.sum(alpha*alpha))*np.sum(dbeta))
Q = 2./det*( (-np.sum(alpha)*np.sum(beta*beta)+np.sum(beta)*np.sum(beta*alpha))*np.sum(dif_dat) \
                 +(np.sum(beta*beta)-np.sum(beta)**2)*np.sum(dalpha) \
                 +(-np.sum(beta*alpha)+np.sum(alpha)*np.sum(beta))*np.sum(dbeta))
U = 2./det*( (np.sum(alpha*beta)*np.sum(alpha) - np.sum(beta)*np.sum(alpha*alpha) )*np.sum(dif_dat) \
                 +(-np.sum(beta*alpha)+np.sum(alpha)*np.sum(beta))*np.sum(dalpha) \
                 +(np.sum(alpha*alpha)-np.sum(alpha)**2)*np.sum(dbeta))
print 'T', T
print 'S', S
print 'Q', Q
print 'U', U
print 'P/T', np.sqrt(Q**2+U**2)/T
print ''

py.figure(figsize=(15,8))

py.subplot(241)
py.plot(top_psi*radeg,'.')
py.plot(bot_psi*radeg,'.')
py.title('top, bot psi')

py.subplot(242)
py.plot(top_psi*radeg-bot_psi*radeg,'.')
py.title('dif_psi')

py.subplot(243)
py.plot(alpha,'b.')
py.plot(beta,'r.')
py.title('alpha, beta')

py.subplot(244)
mylib.plot_hist(alpha,30,fit=True,init_auto=True)
mylib.plot_hist(beta,30,fit=True,init_auto=True)

py.subplot(245)
py.plot(sum_dat,'.')
py.title('sum_dat')

py.subplot(246)
py.plot(dif_dat,'.')
py.title('dif_dat')

py.subplot(247)
py.plot(ra*radeg,dec*radeg,'.')
py.title('ra,dec')

py.savefig('tmp.png')
os.system('display tmp.png &')
