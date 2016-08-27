import numpy as np
import healpy as h
import lib_mapprep as lib_m
import sys

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/000uK/lensed/r00001/k10/0051-0100/PSM_OUTPUT/skyinbands/PSM_IDEAL/band10/cmb_map_band10.fits'
map60_o = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/000uK/lensed/r00001/k10/0051-0100/PSM_OUTPUT/skyinbands/PSM_IDEAL/band28/cmb_map_band28.fits'
map78_o = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/000uK/lensed/r00001/k10/0051-0100/PSM_OUTPUT/skyinbands/PSM_IDEAL/band50/cmb_map_band50.fits'
map100_o = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/000uK/lensed/r00001/k10/0101-0150/PSM_OUTPUT/skyinbands/PSM_IDEAL/band40/cmb_map_band40.fits'
map140_o = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/000uK/lensed/r00001/k10/0151-0200/PSM_OUTPUT/skyinbands/PSM_IDEAL/band45/cmb_map_band45.fits'
map195_o = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/000uK/lensed/r00001/k10/0251-0300/PSM_OUTPUT/skyinbands/PSM_IDEAL/band30/cmb_map_band30.fits'
map280_o = h.read_map(filename, field=0)

#++++++++

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/processed_maps/60GHz/thermo/fits/cmb_map_band10_b0arcmin_ud1024_thermo.fits'
map60 = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/processed_maps/78GHz/thermo/fits/cmb_map_band28_b0arcmin_ud1024_thermo.fits'
map78 = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/processed_maps/100GHz/thermo/fits/cmb_map_band50_b0arcmin_ud1024_thermo.fits'
map100 = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/processed_maps/140GHz/thermo/fits/cmb_map_band40_b0arcmin_ud1024_thermo.fits'
map140 = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/processed_maps/195GHz/thermo/fits/cmb_map_band45_b0arcmin_ud1024_thermo.fits'
map195 = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/processed_maps/280GHz/thermo/fits/cmb_map_band30_b0arcmin_ud1024_thermo.fits'
map280 = h.read_map(filename, field=0)

#++++++++

A2T_60 = lib_m.Antenna2Thermo(60)
A2T_78 = lib_m.Antenna2Thermo(78)
A2T_100 = lib_m.Antenna2Thermo(100)
A2T_140 = lib_m.Antenna2Thermo(140)
A2T_195 = lib_m.Antenna2Thermo(195)
A2T_280 = lib_m.Antenna2Thermo(280)

print 'Antenna to thermodynamics'
print (A2T_60)
print (A2T_78)
print (A2T_100)
print (A2T_140)
print (A2T_195)
print (A2T_280)
print ''
print 'Antenna to thermodynamics, ratio w.r.t. 60GHz'
print (A2T_60/A2T_60)
print (A2T_60/A2T_78)
print (A2T_60/A2T_100)
print (A2T_60/A2T_140)
print (A2T_60/A2T_195)
print (A2T_60/A2T_280)
print ''
print 'RJ unit, map(*)/map(60GHz)'
print np.mean(map60_o/map60_o)
print np.mean(map78_o/map60_o)
print np.mean(map100_o/map60_o)
print np.mean(map140_o/map60_o)
print np.mean(map195_o/map60_o)
print np.mean(map280_o/map60_o)
print ''
print 'RMS amplitude of the map in RJ'
print np.std(map60_o)
print np.std(map78_o)
print np.std(map100_o)
print np.std(map140_o)
print np.std(map195_o)
print np.std(map280_o)
print ''
print 'RMS amplitude of the map in thermodynamics'
print np.std(map60)
print np.std(map78)
print np.std(map100)
print np.std(map140)
print np.std(map195)
print np.std(map280)
print ''
print 'thermodynamic unit, mean map(*)/map(60GHz)'
print np.mean(map60/map60)
print np.mean(map78/map60)
print np.mean(map100/map60)
print np.mean(map140/map60)
print np.mean(map195/map60)
print np.mean(map280/map60)

#++++++++
sys.exit()
filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM/processed_maps/60GHz/thermo/fits/cmb_map_band10_b75arcmin_ud512_thermo.fits'
map60_s = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM/processed_maps/78GHz/thermo/fits/cmb_map_band28_b0arcmin_ud1024_thermo.fits'
map78_s = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM/processed_maps/100GHz/thermo/fits/cmb_map_band50_b0arcmin_ud1024_thermo.fits'
map100_s = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM/processed_maps/140GHz/thermo/fits/cmb_map_band40_b0arcmin_ud1024_thermo.fits'
map140_s = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM/processed_maps/195GHz/thermo/fits/cmb_map_band45_b0arcmin_ud1024_thermo.fits'
map195_s = h.read_map(filename, field=0)

filename='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM/processed_maps/280GHz/thermo/fits/cmb_map_band30_b0arcmin_ud1024_thermo.fits'
map280_s = h.read_map(filename, field=0)

#++++++++


