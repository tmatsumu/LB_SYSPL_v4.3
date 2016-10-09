import numpy as np
import pylab as py
import healpy as h
from scipy.optimize import curve_fit
import os
import sys


def func(x, a, b):
    return a * x**b

pi = np.pi
radeg = (180./pi)

const = 1.

dir = '/group/cmb/litebird/usr/tmatsumu/LB_SYSPL_v4.3/SimedMaps/RunLog/'

#sim_option='bias'
#sim_option='random_r'
sim_option='1of_r'

subdir1 = '/gainval_v001/coadd_map/coadd_map1/'

if sim_option == 'bias':
    subdir2 = '/gainval_v002/coadd_map/coadd_map1/'
    subdir3 = '/gainval_v003/coadd_map/coadd_map1/'
    subdir4 = '/gainval_v004/coadd_map/coadd_map1/'
    subdir5 = '/gainval_v005/coadd_map/coadd_map1/'
    subdir6 = '/gainval_v006/coadd_map/coadd_map1/'
if sim_option == 'random_r':
    subdir2 = '/gainval_v012/coadd_map/coadd_map1/'
    subdir3 = '/gainval_v013/coadd_map/coadd_map1/'
    subdir4 = '/gainval_v014/coadd_map/coadd_map1/'
    subdir5 = '/gainval_v015/coadd_map/coadd_map1/'
    subdir6 = '/gainval_v016/coadd_map/coadd_map1/'
if sim_option == '1of_r':
    subdir2 = '/gainval_v022/coadd_map/coadd_map1/'
    subdir3 = '/gainval_v023/coadd_map/coadd_map1/'
    subdir4 = '/gainval_v024/coadd_map/coadd_map1/'
    subdir5 = '/gainval_v025/coadd_map/coadd_map1/'
    subdir6 = '/gainval_v026/coadd_map/coadd_map1/'

Clin =  h.read_cl(dir+subdir1+'Clin_anafast.fits')
Cl1 =  h.read_cl(dir+subdir1+'Clout_anafast.fits')
Cl2 =  h.read_cl(dir+subdir2+'Clout_anafast.fits')
Cl3 =  h.read_cl(dir+subdir3+'Clout_anafast.fits')
Cl4 =  h.read_cl(dir+subdir4+'Clout_anafast.fits')
Cl5 =  h.read_cl(dir+subdir5+'Clout_anafast.fits')
Cl6 =  h.read_cl(dir+subdir6+'Clout_anafast.fits')

num = len(Clin[0])
ell = np.arange(num)

py.figure(0,figsize=(15,8))
py.subplot(121)
py.plot(ell,ell*(ell+1.)/(2.*pi)*Clin[0],'-g')
py.plot(ell,ell*(ell+1.)/(2.*pi)*Clin[1],'-b')
py.plot(ell,ell*(ell+1.)/(2.*pi)*Clin[2],'-r')

py.plot(ell,ell*(ell+1.)/(2.*pi)*Cl6[2],'.b', label='$g_b=5\cdot10^{-2}$')
py.plot(ell,ell*(ell+1.)/(2.*pi)*Cl5[2],'.g', label='$g_b=1\cdot10^{-2}$')
py.plot(ell,ell*(ell+1.)/(2.*pi)*Cl4[2],'.k', label='$g_b=5\cdot10^{-3}$')
py.plot(ell,ell*(ell+1.)/(2.*pi)*Cl3[2],'.y', label='$g_b=1\cdot10^{-3}$')
py.plot(ell,ell*(ell+1.)/(2.*pi)*Cl2[2],'.m', label='$g_b=5\cdot10^{-4}$')

py.plot(ell,ell*(ell+1.)/(2.*pi)*Clin[2]/100.,'--r', label='$10^{-2}\cdot C_{\ell}^{Lens}$')

py.xlim([2,1000])
py.ylim([1e-8,1e4])
py.legend(loc='best',prop={'size':15})
py.xlabel('$\ell$', fontsize=17)
py.ylabel('$\ell(\ell+1)/2\pi C_{\ell}$ [$\mu$K$^2$]', fontsize=17)
py.xticks( color = 'k', size = 17)
py.yticks( color = 'k', size = 17)
py.loglog()

gain_arr = np.array([5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 0])
Cl_ave6 = ell*(ell+1.)/(2.*pi)*Cl6[2]
Cl_ave5 = ell*(ell+1.)/(2.*pi)*Cl5[2]
Cl_ave4 = ell*(ell+1.)/(2.*pi)*Cl4[2]
Cl_ave3 = ell*(ell+1.)/(2.*pi)*Cl3[2]
Cl_ave2 = ell*(ell+1.)/(2.*pi)*Cl2[2]
Cl_ave = ell*(ell+1.)/(2.*pi)*Clin[2]
Cl_ave_arr = np.array( [np.mean(Cl_ave6[0:10]), np.mean(Cl_ave5[0:10]), np.mean(Cl_ave4[0:10]), np.mean(Cl_ave3[0:10]), np.mean(Cl_ave2[0:10]), np.mean(Cl_ave[0:10]) ])
Cl_std_arr = np.array( [np.std(Cl_ave6[0:10]), np.std(Cl_ave5[0:10]), np.std(Cl_ave4[0:10]), np.std(Cl_ave3[0:10]), np.std(Cl_ave2[0:10]), np.std(Cl_ave[0:10]) ])

py.subplot(122)

x_tmp = np.hstack((np.array([0]), gain_arr[0:5]))
y_tmp = np.hstack((np.array([0]), Cl_ave_arr[0:5] ))
par = np.polyfit(x_tmp,y_tmp,1)
par, pcov = curve_fit(func, x_tmp, y_tmp)
print par
py.errorbar(gain_arr[0:5], Cl_ave_arr[0:5], Cl_std_arr[0:5], fmt='ob')
py.plot(x_tmp,func(x_tmp,par[0],par[1]),'--b')
py.plot(gain_arr[0:5], Cl_ave_arr[5]*np.ones(len(gain_arr[0:5])), 'r--', label='mean of $C_{\ell}^{Lens}$')
py.plot(gain_arr[0:5], Cl_ave_arr[5]*np.ones(len(gain_arr[0:5]))/100., 'g--', label='mean of $10^{-2}\cdot C_{\ell}^{Lens}$')
py.xlabel('Gain, g', fontsize=17)
py.ylabel('mean of $\ell(\ell+1)/2\pi C_{l<10}$', fontsize=17)
py.title('$%e \\times g^{%1.3f}$' % (par[0], par[1]))
py.xticks( color = 'k', size = 17)
py.yticks( color = 'k', size = 17)
py.loglog()

#py.show()
py.savefig(sim_option+'.png')
