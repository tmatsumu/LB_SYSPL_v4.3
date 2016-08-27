import numpy as np
import fileinput
import sys
import pylab as py
import lib_mapmaker as lmm
import lib_qsub as lqsub
import ReadMapMakeXml as rxml
import glob
import time as time
import sqlite3 as sq
import os

xml_filename = sys.argv[1]
select_ces_in = sys.argv[2]

def gen_mapfilelist():
    print ""
    print "  [main_coaddDB.py:gen_mapfilelist]: START", xml_filename
    xml_input = rxml.Get_Mapmake_Inputs(xml_filename) 
    dir_out = xml_input["dir_simedmap"]
    
    # os.system('mkdir '+dir_out+'/coadd_map')
#    dir_ces = glob.glob(dir_out+'/ces*')
    dir_ces = glob.glob(dir_out+'/day*')
    nb_ces = len(dir_ces)
    
    run_id = []; path = []; pix = []
    for i in range(nb_ces):
        print '  [main_coaddDB.py:gen_mapfilelist] ', dir_ces[i]+'/mapfilelist.db'
        conn = sq.connect(dir_ces[i]+'/mapfilelist.db')
        c = conn.cursor()
        c.execute('select * from mapfilelist')
        for ar in c:
            run_id.append(ar[0])
            path.append(ar[1])
            pix.append(ar[2])
        c.close()

    readDB = lmm.readDB()
    readDB.xml=xml_input
    readDB.sq_command = select_ces_in
#    readDB.sq_command = xml_input['sqlite_command_ces']
#    CES = readDB.read_pb1_observation() 
    CES = readDB.read_lb_observation() 

    os.system('rm -f '+dir_out+'/mapfilelist_all.db')
#    nb = len(pix)
    nb = len(CES['id'])
    conn = sq.connect(dir_out+'/mapfilelist_all.db')
    c = conn.cursor()
#    c.execute('create table mapfilelist (id integer, run_id integer, run_subid integer, dir_ptg text, outdir text, pix integer)')
    c.execute('create table mapfilelist (id integer, juliantime real, ourdir text)')
    for i in range(0,nb):
#        list_entries = (run_id[i], path[i], pix[i])
#        list_entries = (CES['id'][i], run_id[i], CES['run_subid'][i], CES['dir_ptg'][i], path[i], pix[i])
        list_entries = (CES['id'][i], CES['juliantime'][i], path[i])
        c.execute('insert into mapfilelist values (?,?,?)', list_entries)
    conn.commit()
    c.close()
    print ""

def gen_coaddmap_db():
    print ""
    print "  [main_coaddDB.py:gen_coaddmap_db]: START", xml_filename
    xml_input = rxml.Get_Mapmake_Inputs(xml_filename)
    print xml_input["dir_simedmap"]+'/coadd_map.db'

    print "  [main_coaddDB.py:gen_coaddmap_db] Read database"
    if os.path.exists(xml_input["dir_simedmap"]+'/coadd_map.db'):
        conn = sq.connect(xml_input["dir_simedmap"]+'/coadd_map.db')
        c = conn.cursor()
        c.execute('select * from coadd_map_db')
        id=[]; sys_run_name=[]; outdir=[]; dir_coadd=[]; select_ces=[]
        for ar in c:
            id.append(int(ar[0]))
            sys_run_name.append(str(ar[1]))
            outdir.append(str(ar[2]))
            dir_coadd.append(str(ar[3]))
            select_ces.append(str(ar[4]))
        c.close()
        ind = np.where( (np.array(select_ces) == select_ces_in) )

        if len(ind[0]) == 1:
            print "  [main_coaddDB.py:gen_coaddmap_db] There exists the same entry in the DB."
        if len(ind[0]) == 0:
            conn = sq.connect(xml_input["dir_simedmap"]+'/coadd_map.db')
            c = conn.cursor()
            list_entries = (None, 
                            xml_input["runID"], 
                            xml_input["dir_simedmap"],
                            'coadd_map'+str(max(id)+1), 
                            select_ces_in)
            c.execute('insert into coadd_map_db values (?,?,?,?,?)', list_entries)
            conn.commit()
            c.close()
            os.system("mkdir "+xml_input["dir_simedmap"]+"/coadd_map/coadd_map"+str(max(id)+1))

    if not os.path.exists(xml_input["dir_simedmap"]+'/coadd_map.db'):
        print ''
        print '  [main_coaddDB.py] NO '+xml_input["dir_simedmap"]+'/coadd_map.db, creating....'
        conn = sq.connect(xml_input["dir_simedmap"]+'/coadd_map.db')
        c = conn.cursor()
        c.execute('create table coadd_map_db (id integer primary key, sys_run_name text, dir_out text, dir_coadd text, select_ces text)')
        list_entries = (None, 
                        xml_input["runID"], 
                        xml_input["dir_simedmap"],
                        'coadd_map1',
                        select_ces_in)
        c.execute('insert into coadd_map_db values (?,?,?,?,?)', list_entries)
        conn.commit()
        c.close()
        os.system("mkdir "+xml_input["dir_simedmap"]+"/coadd_map/coadd_map1")
    print ""

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print ''
print ''
print '[main_coaddDB.py] begin gen_mapfilelist()'
gen_mapfilelist()
print '[main_coaddDB.py] end gen_mapfilelist()'
print ''

print ''
print '[main_coaddDB.py] begin gen_coaddmap_db()'
gen_coaddmap_db()
print '[main_coaddDB.py] end gen_coaddmap_db()'
print ''
print ''

