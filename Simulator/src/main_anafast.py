import numpy as np
import pylab as py
import healpy as h
import ReadMapMakeXml as rxml
import sys
import os

xml_filename = sys.argv[1]
xml_input = rxml.Get_Mapmake_Inputs(xml_filename)
if sys.argv[3] == 'anafast':
    name_extTOD = 'coadd_map1'
if sys.argv[3] == 'extTODmm_anafast':
    if sys.argv[4] == '': 
        print '[main_anafast.py] WARNING! add the directory name for extTODmm'
        sys.exit()
    name_extTOD = 'coadd_map1_'+sys.argv[4]

#+++++++++++++++++++++++++++++++++++++++++++++++++

os.system('rm -f '+xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/Clout_anafast.fits')
fname_anafast = xml_input['dir_simulator']+'/RunJob/'+xml_input['runID']+'/anafast_'+xml_input['runID']+'_out'
f_ana = open(fname_anafast+'.param', 'w')
f_ana.write('infile = %s\n' % (xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/mapTQU.fits'))
f_ana.write('outfile = %s\n' % (xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/Clout_anafast.fits'))
f_ana.write('simul_type = 2 \n')
f_ana.write('nlmax = %d \n' % (xml_input['nside']*2))
#f_ana.write('maskfile = %s \n' % (xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/mapH.fits'))
f_ana.close()

os.system('bsub -q e -o '+fname_anafast+'.o'+' anafast -d '+fname_anafast+'.param')

#+++++++++++++++++++++++++++++++++++++++++++++++++

os.system('rm -f '+xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/ClHout_anafast.fits')
fname_anafast = xml_input['dir_simulator']+'/RunJob/'+xml_input['runID']+'/anafast_'+xml_input['runID']+'_Hout'
f_ana = open(fname_anafast+'.param', 'w')
f_ana.write('infile = %s\n' % (xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/mapH.fits'))
f_ana.write('outfile = %s\n' % (xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/ClHout_anafast.fits'))
f_ana.write('simul_type = 1 \n')
f_ana.write('nlmax = %d \n' % (xml_input['nside']*2))
f_ana.close()

os.system('bsub -q e -o '+fname_anafast+'.o'+' anafast -d '+fname_anafast+'.param')

#+++++++++++++++++++++++++++++++++++++++++++++++++

os.system('rm -f '+xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/ClSout_anafast.fits')
fname_anafast = xml_input['dir_simulator']+'/RunJob/'+xml_input['runID']+'/anafast_'+xml_input['runID']+'_Sout'
f_ana = open(fname_anafast+'.param', 'w')
f_ana.write('infile = %s\n' % (xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/mapS.fits'))
f_ana.write('outfile = %s\n' % (xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/ClSout_anafast.fits'))
f_ana.write('simul_type = 1 \n')
f_ana.write('nlmax = %d \n' % (xml_input['nside']*2))
f_ana.close()

os.system('bsub -q e -o '+fname_anafast+'.o'+' anafast -d '+fname_anafast+'.param')

#+++++++++++++++++++++++++++++++++++++++++++++++++

os.system('rm -f '+xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/Clin_anafast.fits')
fname_anafast = xml_input['dir_simulator']+'/RunJob/'+xml_input['runID']+'/anafast_'+xml_input['runID']+'_in'
f_ana = open(fname_anafast+'.param', 'w')
f_ana.write('infile = %s\n' % (xml_input["file_input_maps"]+'.fits'))
f_ana.write('outfile = %s\n' % (xml_input["dir_simedmap"]+'/coadd_map/'+name_extTOD+'/Clin_anafast.fits'))
f_ana.write('simul_type = 2 \n')
f_ana.write('nlmax = %d \n' % (xml_input['nside']*2))
f_ana.close()

os.system('bsub -q e -o '+fname_anafast+'.o'+' anafast -d '+fname_anafast+'.param')

