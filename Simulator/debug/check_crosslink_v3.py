import numpy as np
import pylab as py
import healpy as h
import sys 

'''
   readin the map from 
       Simulator/data_proj/SimedMaps/RunLog/
   and compute the crosslink and other crosslink related files
   2016-8-18, wrtten by T. Matsumura

'''

dir_in = sys.argv[1]

filename_H = 'mapH'
filename_arr = ['mapAA_CNexcluded', 'mapBB_CNexcluded', 'mapAB_CNexcluded']

mapH = h.read_map(dir_in+'/'+filename_H+'.fits')
num = len(filename_arr)
for i in range(num):
    mapin = h.read_map(dir_in+'/'+filename_arr[i]+'.fits')
    h.mollview(mapin/np.float_(mapH), title=filename_arr[i],max=max(mapin/np.float_(mapH)),min=min(mapin/np.float_(mapH)))
    py.savefig(filename_arr[i]+'.png')


mapA = h.read_map(dir_in+'/'+'mapA_CNexcluded.fits')
mapB = h.read_map(dir_in+'/'+'mapB_CNexcluded.fits')
crosslink = (mapA**2 + mapB**2)/np.float_(mapH)**2

h.mollview(crosslink, title='crosslink',max=1,min=0)
py.savefig('crosslink.png')

h.mollview(mapH, title='mapH',max=max(mapH),min=0)
py.savefig('mapH.png')
