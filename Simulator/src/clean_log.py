import os
import sys
import ReadMapMakeXml as rxml

option = sys.argv[2]

xml_filename = sys.argv[1]
print "[clean_log.py]: xml_filename=", xml_filename
xml_input = rxml.Get_Mapmake_Inputs(xml_filename)

if 'mm' in option:
    os.system('mkdir '+xml_input['dir_simedmap']+'/log')
    os.system('mkdir '+xml_input['dir_simedmap']+'/par')
    os.system('mv PB_mm_*.e* '+xml_input['dir_simedmap']+'/log')
    
if 'coadd' in option:
    os.system('mv PB_coadd_*.e* '+xml_input['dir_simedmap']+'/log')
    os.system('rm -rf '+xml_input['dir_simedmap']+'/coadd_map/*.npy')

if 'clean_npy' in option:
    print 'rm -f '+xml_input['dir_simedmap']+'/input_npy/*.npy'
    os.system('rm -f '+xml_input['dir_simedmap']+'/input_npy/*.npy')    

if 'clean_npz' in option:
    print '[clean_npz]'
    print 'rm -f '+xml_input['dir_simedmap']+'/day*/*.npz'
    os.system('rm -f '+xml_input['dir_simedmap']+'/day*/*.npz')

if 'clean_tod' in option:
    print '[clean_tod]'
    print 'rm -rf '+xml_input['dir_simedmap']+'/day*/tod'
    os.system('rm -rf '+xml_input['dir_simedmap']+'/day*/tod')

if 'xpure_window' in option:
    os.system('mv '+xml_input['dir_simulator']+'/gen_run/xpure_window_'+xml_input['runID']+'.e* '+xml_input['dir_simedmap']+'/log')
    os.system('mv '+xml_input['dir_simulator']+'/gen_run/output_*_'+xml_input['runID']+'.log '+xml_input['dir_simedmap']+'/log')
    os.system('mv '+xml_input['dir_simulator']+'/gen_run/param_*_'+xml_input['runID']+'.par '+xml_input['dir_simedmap']+'/par')

if 'xpure_cl' in option:
    pass
#    os.system('mv '+xml_input['dir_simulator']+'/gen_run/param_*_'+xml_input['runID']+'')

if 'clean' in option:
    pass
