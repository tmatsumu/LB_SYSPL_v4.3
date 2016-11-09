import numpy as np
import pylab as py
import healpy as h
import sys
import os

dir_in = sys.argv[1]
fits_filename = sys.argv[2]
nside = int(sys.argv[3])
#+++++++++++++++++++++++++++++++++++++++++++++++++

#os.system('rm -f '+dir_in+'/Clout_anafast.fits')
f_ana = open(fits_filename+'.param', 'w')
f_ana.write('infile = %s\n' % (dir_in + '/' + fits_filename+'.fits'))
f_ana.write('outfile = %s\n' % (dir_in+'/Cl_'+fits_filename+'.fits'))
f_ana.write('simul_type = 1 \n')
f_ana.write('nlmax = %d \n' % (nside*2))
#f_ana.write('maskfile = %s \n' % (xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/mapH.fits'))
f_ana.close()

os.system('bsub -q s -o '+fits_filename+'.o'+' anafast -d '+fits_filename+'.param')

#+++++++++++++++++++++++++++++++++++++++++++++++++

