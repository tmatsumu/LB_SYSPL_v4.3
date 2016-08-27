import numpy as np
import pylab as py
import healpy as h
import glob
import sys
import os
import matsumulib as mylib

pi = np.pi
radeg = (180./pi)

num = 100
top_psi = np.random.uniform(0.,pi,num)
bot_psi = top_psi + 90./radeg

Tin = 10.
Qin = 0.
Uin = 0.
del_gain = 0.0001
#top_tod = 10.*np.ones(num)
#bot_tod = 10.*np.ones(num)*(1.-del_gain)
top_tod = Tin*np.ones(num) + Qin*np.cos(-2.*top_psi) + Uin*np.sin(-2.*top_psi)
bot_tod = Tin*np.ones(num)*(1.-del_gain) + Qin*np.cos(-2.*bot_psi) + Uin*np.sin(-2.*bot_psi)
sum_dat = 0.5*(top_tod+bot_tod)
dif_dat = 0.5*(top_tod-bot_tod)

alpha = 0.5*(np.cos(-2.*top_psi) - np.cos(-2.*bot_psi))
beta = 0.5*(np.sin(-2.*top_psi) - np.sin(-2.*bot_psi))

print ''
print 'num,', num
print '(sum A, sum A^2, sum A^2)', np.sum(alpha), np.sum(alpha**2), np.sum(alpha)**2
print '(sum B, sum B^2, sum B^2)', np.sum(beta), np.sum(beta**2), np.sum(beta)**2
print '(sum AB, sum A sum B)', np.sum(alpha*beta), np.sum(alpha)*np.sum(beta)
det = num*np.sum(alpha*alpha)*np.sum(beta*beta)+2.*np.sum(alpha)*np.sum(beta)*np.sum(alpha*beta) \
    - np.sum(alpha*alpha)*np.sum(beta)**2 - np.sum(alpha)**2*np.sum(beta*beta)-num*np.sum(alpha*beta)**2
dalpha = dif_dat*alpha
dbeta = dif_dat*beta
print '(sum dA, sum dB)', np.sum(dalpha), np.sum(dbeta)
print 'sum d', np.sum(dif_dat)
print 'mean_difftod', np.mean(dif_dat)
print 'det', det
print ''
print 'inverse matrix'
print  num*(np.sum(alpha*alpha)*np.sum(beta*beta)-np.sum(alpha*beta)**2)/det
print (-np.sum(alpha)*np.sum(beta*beta)+np.sum(beta)*np.sum(beta*alpha))/det, (num*np.sum(beta*beta)-np.sum(beta)**2)/det
print (np.sum(alpha*beta)*np.sum(alpha) - np.sum(beta)*np.sum(alpha*alpha) )/det, num*(-np.sum(beta*alpha)+np.sum(alpha)*np.sum(beta))/det, (num*np.sum(alpha*alpha)-np.sum(alpha)**2)/det
print ''
T = np.sum(sum_dat)/np.float(len(sum_dat))
S = 1./det*( (np.sum(alpha*alpha)*np.sum(beta*beta)-np.sum(alpha*beta)**2)*np.sum(dif_dat) \
                 +(-np.sum(alpha)*np.sum(beta*beta)+np.sum(beta)*np.sum(beta*alpha))*np.sum(dalpha) \
                 +(np.sum(alpha)*np.sum(beta*alpha)-np.sum(beta)*np.sum(alpha*alpha))*np.sum(dbeta))
Q = 1./det*( (-np.sum(alpha)*np.sum(beta*beta)+np.sum(beta)*np.sum(beta*alpha))*np.sum(dif_dat) \
                 +(num*np.sum(beta*beta)-np.sum(beta)**2)*np.sum(dalpha) \
                 +(-num*np.sum(beta*alpha)+np.sum(alpha)*np.sum(beta))*np.sum(dbeta))
U = 1./det*( (np.sum(alpha*beta)*np.sum(alpha) - np.sum(beta)*np.sum(alpha*alpha) )*np.sum(dif_dat) \
                 +(-num*np.sum(beta*alpha)+np.sum(alpha)*np.sum(beta))*np.sum(dalpha) \
                 +(num*np.sum(alpha*alpha)-np.sum(alpha)**2)*np.sum(dbeta))
print 'del_gain', del_gain
print 'T', T, Tin, (T-Tin)/Tin
print 'S', S
print 'Q', Q, Qin, (Q-Qin)/Qin
print 'U', U, Uin, (U-Uin)/Uin
print 'P/T', np.sqrt(Q**2+U**2)/T, np.sqrt(Qin**2+Uin**2)/Tin, (np.sqrt(Q**2+U**2)/T - np.sqrt(Qin**2+Uin**2)/Tin)/(np.sqrt(Qin**2+Uin**2)/Tin)
print ''

py.figure(figsize=(10,8))

py.subplot(231)
py.plot(top_psi*radeg,'.')
py.plot(bot_psi*radeg,'.')
py.title('top, bot psi')

py.subplot(232)
py.plot(top_psi*radeg-bot_psi*radeg,'.')
py.title('dif_psi')

py.subplot(233)
py.plot(alpha,'b.')
py.plot(beta,'r.')
py.title('alpha, beta')

py.subplot(234)
mylib.plot_hist(alpha,30,fit=True,init_auto=True)
mylib.plot_hist(beta,30,fit=True,init_auto=True)
py.title('hist of sin(2alpha) and cos(2alpha)')

py.subplot(235)
py.plot(sum_dat,'.')
py.title('sum_dat')

py.subplot(236)
py.plot(dif_dat,'.')
py.title('dif_dat')

#py.subplot(247)
#py.plot(ra*radeg,dec*radeg,'.')
#py.title('ra,dec')

py.savefig('tmp.png')
os.system('display tmp.png &')
