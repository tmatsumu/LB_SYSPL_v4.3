import numpy as np
import lib_qsub as lib_q
import ReadMapMakeXml as rxml
import lib_qsub as lib_q
import lib_dbutil as lib_db
import os
import sys

xml_filename = sys.argv[1]
xml_input = rxml.Get_Mapmake_Inputs(xml_filename)
sqlite_command = sys.argv[2] 
exe_type = sys.argv[3]

if exe_type == 'gen_mappng':
    Tmax = sys.argv[4]
    Pmax = sys.argv[5]

    out = lib_q.gen_mappng()
    out.runID = xml_input['runID']
    out.dir_simulator = xml_input['dir_simulator']
    out.sqlite_command = sqlite_command
    out.xml_filename = xml_filename
    out.Tmax = Tmax
    out.Pmax = Pmax
    filename_SH = out.gen_SH()
    filename_PBS = out.gen_PBS()
    os.system('qsub '+filename_PBS)

if exe_type == 'cp_png':
    out = lib_db.read_coaddDB()
    out.xml_input = xml_input
    out.sqlite_command = sqlite_command
    outdirs, dir_coadds = out.read_coaddDB() 
    os.system('mkdir -p ~/www/tmp/'+xml_input['runID']+'/'+dir_coadds)
    os.system('cp '+outdirs+'/png/'+dir_coadds+'/*.png '+'~/www/tmp/'+xml_input['runID']+'/'+dir_coadds)
    print ''
    print 'ls -hlrt ~/www/tmp/'+xml_input['runID']+'/'+dir_coadds
    os.system('ls -hlrt ~/www/tmp/'+xml_input['runID']+'/'+dir_coadds)
    os.system('chmod 777 '+'~/www/tmp/'+xml_input['runID'] )
    os.system('chmod 777 '+'~/www/tmp/'+xml_input['runID']+'/'+dir_coadds)
    os.system('chmod 777 '+'~/www/tmp/'+xml_input['runID']+'/'+dir_coadds+'/*' )
    print ''
