import numpy as np
import pylab as py
import healpy as h
import os
import sys
import lib_mapprep as lib_m

pi = np.pi 
radeg = (180./pi)

'''
   generate the CMB map seen by sidelobe
   this includes not only CMB but dipole
   2014-2-4, T. Matsumura
'''

#+++++++++
# READ MAP

filename_CMB = sys.argv[1]
filename_out = sys.argv[2]
nside_out = int(sys.argv[3])
nu_GHz = float(sys.argv[4])
fwhm_in = float(sys.argv[5])
db = float(sys.argv[6])
supression = lib_m.db_convert(db)
option_dipole = sys.argv[7]
option_unit = sys.argv[8]
option_FG = sys.argv[9]

print ''
print '+++++++++++++++++++++++++++++++'
print filename_CMB
print filename_out
print nside_out
print fwhm_in/60., 'degs'
print db, 'db'
print supression
print option_dipole
print option_unit
print ''

mapT_CMB = h.read_map(filename_CMB+'.fits',field=0)
mapQ_CMB = h.read_map(filename_CMB+'.fits',field=1)
mapU_CMB = h.read_map(filename_CMB+'.fits',field=2)

npix = len(mapT_CMB)
nside = h.npix2nside(npix)

if ((option_unit=='RJ') or (option_unit=='Antenna')):
    A2T = 1.
if option_unit=='thermo':
    A2T = lib_m.Antenna2Thermo(nu_GHz)
    
mapT_CMB = mapT_CMB * A2T
mapQ_CMB = mapQ_CMB * A2T
mapU_CMB = mapU_CMB * A2T

if option_dipole=='y':
    mapDipole = lib_m.gen_dipole(nside,'uK')
    if option_unit=='thermo': pass
    if ((option_unit=='RJ') or (option_unit=='Antenna')):
        mapDipole = mapDipole/A2T

if option_dipole=='n':
    mapDipole = np.zeros(npix)

if option_FG=='n':
    map_in = np.array([ (mapT_CMB+mapDipole)*supression, \
                            (mapQ_CMB)*supression, \
                            (mapU_CMB)*supression])
    mapT_CMB=mapQ_CMB=mapU_CMB=mapDipole=0
    
if option_FG=='y':
    filename_synch = sys.argv[10]
    filename_dust = sys.argv[11]
    mapT_synch = h.read_map(filename_synch+'.fits',field=0)
    mapQ_synch = h.read_map(filename_synch+'.fits',field=1)
    mapU_synch = h.read_map(filename_synch+'.fits',field=2)
    
    mapT_dust = h.read_map(filename_dust+'.fits',field=0)
    mapQ_dust = h.read_map(filename_dust+'.fits',field=1)
    mapU_dust = h.read_map(filename_dust+'.fits',field=2)

    mapT_synch=mapT_synch*A2T
    mapQ_synch=mapQ_synch*A2T
    mapU_synch=mapU_synch*A2T

    mapT_dust=mapT_dust*A2T
    mapQ_dust=mapQ_dust*A2T
    mapU_dust=mapU_dust*A2T

    map_in = np.array([ (mapT_CMB+mapT_synch+mapT_dust+mapDipole)*supression, \
                            (mapQ_CMB+mapQ_synch+mapQ_dust)*supression, \
                            (mapU_CMB+mapU_synch+mapU_dust)*supression])
    mapT_CMB=mapT_dust=mapT_synch=mapQ_CMB=mapQ_dust=mapQ_synch=mapU_CMB=mapU_dust=mapU_synch=0
    

print 'nside = ', nside

fwhm_in = fwhm_in/radeg/60.
map_out = h.smoothing(map_in, fwhm=fwhm_in)
map_out_ud = h.ud_grade(map_out,nside_out=nside_out)
h.write_map(filename_out+'_'+option_unit+'.fits',map_out_ud)

h.mollview(map_out_ud[0],title='T_'+option_unit)
py.savefig(filename_out+'_'+option_unit+'_T.png')

h.mollview(map_out_ud[1],title='Q_'+option_unit)
py.savefig(filename_out+'_'+option_unit+'_Q.png')

h.mollview(map_out_ud[2],title='U_'+option_unit)
py.savefig(filename_out+'_'+option_unit+'_U.png')

os.system('')
