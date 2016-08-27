import numpy as np
import fileinput
import healpy as h
import sys
import pylab as py
import lib_mapmaker as lmm
from ReadMapMakeXml import *
import time as time
import os

def read_fitslist(filename):
    dirs = []
    idx = []
    for line in fileinput.input(filename):
        ar = line.split()
        if (len(ar)>1):
            dirs.append(ar[0])
            idx.append(ar[1])
    print "[READ fitslist]: End reading "+filename
    return dirs, np.array(idx)


##################################################################################
time0 = time.time()
print ""
print ""
print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print "++++++++ START: main_mapcombiner.py +++++++++++++++++++++++++++++++++++++++"
print time.time()-time0
xml_file = sys.argv[1]
xml_input = lmm.Gen_scansetdirs(xml_file)
outdir = xml_input["dir_simedmap"]

selected = sys.argv[2]
#filename = sys.argv[1]
dirs, idx = read_fitslist(outdir+"/fitsfilelist_"+selected+".txt")

#fnameT = outdir+'/'+sys.argv[2]
#fnameQ = outdir+'/'+sys.argv[3]
#fnameU = outdir+'/'+sys.argv[4]
fnameT = sys.argv[3]+"_"+selected
fnameQ = sys.argv[4]+"_"+selected
fnameU = sys.argv[5]+"_"+selected
fnameM = sys.argv[6]+"_"+selected
fnameW = sys.argv[7]+"_"+selected
fnameCN = sys.argv[8]+"_"+selected

nside = xml_input["nside"]
npix = h.nside2npix(nside)
#tmp = np.zeros(npix)

mapTn_sum=np.zeros(npix); mapTd_sum=np.zeros(npix);
mapAA_sum=np.zeros(npix); mapBB_sum=np.zeros(npix); mapAB_sum=np.zeros(npix);
mapAd_sum=np.zeros(npix); mapBd_sum=np.zeros(npix)

mapT=np.zeros(npix); mapQ=np.zeros(npix); mapU=np.zeros(npix)

nb = len(idx)
ii = 0
#for i in range(0,1):
for i in range(0,nb):
#    print ""
    print "[MAIN_MAPCOMBINER]: ", dirs[i], ii, '/', nb
    ii += 1

#    print ''
#    print 'read ipix', npix, npix*4.*0./1.e6, 'Mbyte'
    ind = np.load(dirs[i]+'/ipix_'+idx[i]+'.npy')

#    print 'read Tn', npix, npix*4./1.e6, 'Mbyte'
#    mapTn = h.read_map(dirs[i]+'/mapTn_'+idx[i]+'.fits')
    mapTn = np.load(dirs[i]+'/mapTn_'+idx[i]+'.npy')
#    ind = np.where(np.isnan(mapTn))
#    mapTn[ind[0]] = 0.

#    print 'read Td', npix*2, npix*4.*2./1.e6, 'Mbyte'
    mapTd = np.load(dirs[i]+'/mapTd_'+idx[i]+'.npy')
#    mapTd = h.read_map(dirs[i]+'/mapTd_'+idx[i]+'.fits')
#    ind = np.where(np.isnan(mapTd))
#    mapTd[ind[0]] = 0.

#    print 'read AA', npix*3, npix*4.*3./1.e6, 'Mbyte'
    mapAA = np.load(dirs[i]+'/mapAA_'+idx[i]+'.npy')
#    mapAA = h.read_map(dirs[i]+'/mapAA_'+idx[i]+'.fits')
#    ind = np.where(np.isnan(mapAA))
#    mapAA[ind[0]] = 0.

#    print 'read BB', npix*4, npix*4.*4./1.e6, 'Mbyte'
    mapBB = np.load(dirs[i]+'/mapBB_'+idx[i]+'.npy')
#    mapBB = h.read_map(dirs[i]+'/mapBB_'+idx[i]+'.fits')
#    ind = np.where(np.isnan(mapBB))
#    mapBB[ind[0]] = 0.

#    print 'read AB', npix*5, npix*4.*5./1.e6, 'Mbyte'
    mapAB = np.load(dirs[i]+'/mapAB_'+idx[i]+'.npy')
#    mapAB = h.read_map(dirs[i]+'/mapAB_'+idx[i]+'.fits')
#    ind = np.where(np.isnan(mapAB))
#    mapAB[ind[0]] = 0.

#    print 'read Ad', npix*6, npix*4.*6./1.e6, 'Mbyte'
    mapAd = np.load(dirs[i]+'/mapAd_'+idx[i]+'.npy')
#    mapAd = h.read_map(dirs[i]+'/mapAd_'+idx[i]+'.fits')
#    ind = np.where(np.isnan(mapAd))
#    mapAd[ind[0]] = 0.

#    print 'read Bd', npix*7, npix*4.*7./1.e6, 'Mbyte'
    mapBd = np.load(dirs[i]+'/mapBd_'+idx[i]+'.npy')
#    mapBd = h.read_map(dirs[i]+'/mapBd_'+idx[i]+'.fits')
#    ind = np.where(np.isnan(mapBd))
#    mapBd[ind[0]] = 0.
#    print len(mapTn_sum), len(mapTn)
    
    mapTn_sum[ind] += mapTn;    mapTd_sum[ind] += mapTd;
    mapAA_sum[ind] += mapAA;    mapBB_sum[ind] += mapBB;    mapAB_sum[ind] += mapAB;
    mapAd_sum[ind] += mapAd;    mapBd_sum[ind] += mapBd;
    del(mapTn); del(mapTd); del(mapAA); del(mapBB); del(mapAB); del(mapAd); del(mapBd);
    if ((i % 10) == 0): print i,'/', nb, time.time()-time0
#        np.save(outdir+'/mapTn_tmp', mapTn_sum)
#        np.save(outdir+'/mapTd_tmp', mapTd_sum)
#        np.save(outdir+'/mapAA_tmp', mapAA_sum)
#        np.save(outdir+'/mapAB_tmp', mapAB_sum)
#        np.save(outdir+'/mapBB_tmp', mapBB_sum)
#        np.save(outdir+'/mapAd_tmp', mapAd_sum)
#        np.save(outdir+'/mapBd_tmp', mapBd_sum)
#        np.save(outdir+'/log_tmp', [i, idx[i]])

#h.mollview(mapTn_sum)
#h.mollview(mapTd_sum)

mapT = 2.*mapTn_sum/mapTd_sum
mapQ = 2./(mapAA_sum*mapBB_sum-mapAB_sum**2) * (mapBB_sum*mapAd_sum - mapAB_sum*mapBd_sum)
mapU = 2./(mapAA_sum*mapBB_sum-mapAB_sum**2) * (-mapAB_sum*mapAd_sum + mapAA_sum*mapBd_sum)

h.write_map(fnameT+'.nan.fits',mapT)
h.write_map(fnameQ+'.nan.fits',mapQ)
h.write_map(fnameU+'.nan.fits',mapU)
h.write_map(fnameW+'.fits',mapTd_sum)

#############################################################

indT = np.isnan(mapT)
indQ = np.isnan(mapQ)
indU = np.isnan(mapU)

indT2 = np.where(indT == True)
indQ2 = np.where(indQ == True)
indU2 = np.where(indU == True)

mapT[indT2] = 0.
mapQ[indQ2] = 0.
mapU[indU2] = 0.

h.write_map(fnameT+'.fits',mapT)
h.write_map(fnameQ+'.fits',mapQ)
h.write_map(fnameU+'.fits',mapU)

ind = np.where(mapTd_sum > 0)
cn = np.zeros(npix)
for i in ind[0]:
    cn[i] = lmm.Cal_ConditionNum_QUpix(mapAA_sum[i],mapBB_sum[i],mapAB_sum[i])
h.write_map(fnameCN+'.fits',cn)

mask = np.zeros(npix,int)
mask[ind[0]] = 1
h.write_map(fnameM+'.fits',mask)

ind = np.where((mapTd_sum > 0) & (np.abs(cn) < 1.e3))
mask = np.zeros(npix,int)
mask[ind[0]] = 1
h.write_map(fnameM+'_CNexcluded.fits',mask)

ind = np.where(np.abs(cn) > 1.e3)
mapTd_sum[ind[0]] = 0
h.write_map(fnameW+'_CNexcluded.fits',mapTd_sum)
#############################################################

os.popen('mkdir -p '+outdir+'/npy')
os.popen('mv '+outdir+'/*.npy '+outdir+'/npy/')

h.mollview(fnameT+'.nan.fits', title='T')
h.graticule(dpar=10,dmer=10,coord='C')
py.savefig(outdir+'/mapT.png', dpi=DPI)

h.mollview(fnameQ+'.nan.fits', title='Q')
h.graticule(dpar=10,dmer=10,coord='C')
py.savefig(outdir+'/mapQ.png', dpi=DPI)

h.mollview(fnameU+'.nan.fits', title='U')
h.graticule(dpar=10,dmer=10,coord='C')
py.savefig(outdir+'/mapU.png', dpi=DPI)

h.cartview(fnameT+'.nan.fits',coord='C',rot=[-15,-35],lonra=[-20,20],latra=[-20,20],title='T')
h.graticule(dpar=10,dmer=10,coord='C')
py.savefig(outdir+'/map_cart_T.png', dpi=DPI)

h.cartview(fnameQ+'.nan.fits',coord='C',rot=[-15,-35],lonra=[-20,20],latra=[-20,20],title='Q')
h.graticule(dpar=10,dmer=10,coord='C')
py.savefig(outdir+'/map_cart_Q.png', dpi=DPI)

h.cartview(fnameU+'.nan.fits',coord='C',rot=[-15,-35],lonra=[-20,20],latra=[-20,20],title='U')
h.graticule(dpar=10,dmer=10,coord='C')
py.savefig(outdir+'/map_cart_U.png', dpi=DPI)

#py.show()
