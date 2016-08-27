import numpy as np
import lib_PB1SQLDB as lib
import sys

pi = np.pi

dir = sys.argv[1]
filename = sys.argv[2]
filename_out = sys.argv[3]
amp_arcsec = float(sys.argv[4])

readDB = lib.read_DB()
readDB.filename = dir+filename 
beams = readDB.read_BeamParams()

nb = len(beams['boloid'])

xpos_out = []
for i in range(0,nb):
    xpos = beams['xpos'][i]
    if beams['boloname'][i][-1::] == 't':
        xpos = xpos+0.5*amp_arcsec/3600. /180.*pi
    if beams['boloname'][i][-1::] == 'b':
        xpos = xpos-0.5*amp_arcsec/3600. /180.*pi
    beams['xpos'][i] = xpos

makeDB=lib.construct_fpdb()
makeDB.db_name = dir+filename_out
makeDB.beam_params = beams
makeDB.num = nb
makeDB.make_FPDB()
