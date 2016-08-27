import numpy as np
import pylab as py
import healpy as h
import lib_maputil as lib_mu
import lib_dbutil as lib_du
import ReadMapMakeXml as rxml
import sys
import os

xml_filename = sys.argv[1]
xml_input = rxml.Get_Mapmake_Inputs(xml_filename)  
sqlite_command = sys.argv[2]
Tmax = float(sys.argv[3])
Pmax = float(sys.argv[4])

out_db = lib_du.read_coaddDB()
out_db.xml_input = xml_input
out_db.sqlite_command = sqlite_command
out_dir, dir_coadd = out_db.read_coaddDB()

lon = float(xml_input['patch_ra'])
lat = float(xml_input['patch_dec'])
width = float(xml_input['patch_width'])

out = lib_mu.map_util()
out.out_dir = out_dir
out.dir_coadd = dir_coadd
out.def_filename()
out.make_figures(lon,lat,width,Tmax,Pmax,option_part=True)
