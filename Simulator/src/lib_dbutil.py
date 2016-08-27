import numpy as np
import pylab as py
import healpy as h
import sqlite3 as sq
import os

class read_coaddDB():
    def _init_(self):
        self.xml_input = 1
        self.sqlite_command = '.schema'

    def read_coaddDB(self):
        print "[main_coaddmaps.py:read_coaddDB] Read coadd database"
        if os.path.exists(self.xml_input["dir_simedmap"]+'/coadd_map.db'):
            conn = sq.connect(self.xml_input["dir_simedmap"]+'/coadd_map.db')
            c = conn.cursor()
#            c.execute(self.sqlite_command)
            c.execute('select * from coadd_map_db')
            id=[]; sys_run_name=[]; outdir=[]; dir_coadd=[]; select_ces=[]
            for ar in c:
                id.append(ar[0])
                sys_run_name.append(ar[1])
                outdir.append(ar[2])
                dir_coadd.append(ar[3])
                select_ces.append(ar[4])
            c.close()
            ind = np.where( (np.array(select_ces) == self.sqlite_command) )
            self.outdirs = outdir[ind[0]]
            self.dir_coadds = dir_coadd[ind[0]]
            return self.outdirs, self.dir_coadds
