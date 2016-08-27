import numpy as np
import pylab as py
import healpy as h
import sys 

'''

'''

filename1 = '/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v4.1_release/Simulator/data_proj/SimedMaps/RunLog/20160818_140GHz_TQU_1/coadd_map/coadd_map1/mapAA_CNexcluded'
filename2 = '/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v4.1_release/Simulator/data_proj/SimedMaps/RunLog/20160818_140GHz_TQU_2/coadd_map/coadd_map1/mapAA_CNexcluded'
mapin1 = h.read_map(filename1+'.fits')
mapin2 = h.read_map(filename2+'.fits')
h.mollview(mapin1/mapin2, max=2000, min=1)
py.savefig('AAratio.png')

#h.mollview(mapin1/mapin2,norm='log')
#py.savefig('AAratio_log.png')
