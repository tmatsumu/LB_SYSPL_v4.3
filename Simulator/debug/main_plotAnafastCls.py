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

xml_filename = sys.argv[1]
xml_input = rxml.Get_Mapmake_Inputs(xml_filename)
if sys.argv[3] == 'plot_anafast':
    name_extTOD = 'coadd_map1'
if sys.argv[3] == 'extTODmm_plot_anafast':
    if sys.argv[5] == '': 
        print '[main_anafast.py] WARNING! add the directory name for extTODmm'
        sys.exit()
    name_extTOD = 'coadd_map1_'+sys.argv[5]
beam_FWHM_arcmin = float(sys.argv[4])

fname_Clin = xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/Clin_anafast.fits'
Clin = h.read_cl(fname_Clin)
fname_Clout = xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/Clout_anafast.fits'
Clout = h.read_cl(fname_Clout)
fname_ClHout = xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/ClHout_anafast.fits'
ClHout = h.read_cl(fname_ClHout)
fname_ClSout = xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/ClSout_anafast.fits'
ClSout = h.read_cl(fname_ClSout)

num = len(Clin[0])
pixwin = h.pixwin(xml_input["nside"])
B_l_tmp = h.gauss_beam(beam_FWHM_arcmin/60./radeg,lmax=xml_input["nside"]*4)
B_l = B_l_tmp[0:num]*pixwin[0:num]
ell = np.arange(num)

py.subplot(121)

py.plot(ell, ell*(ell+1.)/2./np.pi*Clin[0]*const/B_l**2, '-r', label='TTin')
py.plot(ell, ell*(ell+1.)/2./np.pi*Clin[1]*const/B_l**2, '-g', label='EEin')
py.plot(ell, ell*(ell+1.)/2./np.pi*Clin[2]*const/B_l**2, '-b', label='BBin')
py.plot(ell, ell*(ell+1.)/2./np.pi*Clout[0]*const/B_l**2, '.m', label='TT')
py.plot(ell, ell*(ell+1.)/2./np.pi*Clout[1]*const/B_l**2, '.c', label='EE')
py.plot(ell, ell*(ell+1.)/2./np.pi*Clout[2]*const/B_l**2, '.r', label='BB')
py.plot(ell, ell*(ell+1.)/2./np.pi*ClSout*const/B_l**2, '.y', label='S')
py.loglog()
py.xlim([2,num*2])
py.legend(loc='best',prop={'size':8})
py.xlabel('$l$')
py.ylabel('$l(l+1)/2\pi C_l$ [$\mu K^2$]')

py.subplot(122)
py.plot(ell, np.abs(ell*(ell+1.)/2./np.pi*Clin[3]*const/B_l**2), '-r', label='TEin')
py.plot(ell, np.abs(ell*(ell+1.)/2./np.pi*Clin[4]*const/B_l**2), '-g', label='TBin')
py.plot(ell, np.abs(ell*(ell+1.)/2./np.pi*Clin[5]*const/B_l**2), '-b', label='EBin')
py.plot(ell, np.abs(ell*(ell+1.)/2./np.pi*Clout[3]*const/B_l**2), '.m', label='TE')
py.plot(ell, np.abs(ell*(ell+1.)/2./np.pi*Clout[4]*const/B_l**2), '.c', label='TB')
py.plot(ell, np.abs(ell*(ell+1.)/2./np.pi*Clout[5]*const/B_l**2), '.r', label='EB')
py.loglog()
py.xlim([2,num*2])
py.legend(loc='best',prop={'size':8})
py.xlabel('$l$')
#py.ylabel('$l(l+1)/2\pi C_l$ [$\mu K^2$]')

fname_png = xml_input["dir_simedmap"]+'/png/'+name_extTOD+'/Cl_anafast.png'
py.savefig(fname_png)
