import numpy as np
import sqlite3 as sq
import ReadMapMakeXml as rxml
import glob
import sys
import os

print '[main_coaddmaps.py] beginning of the file'

xml_filename = sys.argv[1]

xml_input = rxml.Get_Mapmake_Inputs(xml_filename)

nside = str(xml_input['nside'])
out_dir = xml_input["dir_simedmap"] 
runID = xml_input['runID']
#n_core = int(xml_input["core"])
#n_cpu = int(xml_input["cpu"])
dir_simulator = xml_input["dir_simulator"]
machine = xml_input["machine"]
mode = xml_input["debug"]

sqlite_command = sys.argv[2]
name_extTOD = sys.argv[3]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class generate_coadd_PBS_SH():
    def __init__(self):
        pass

    def read_DBmapfilelist(self):
        os.system('mkdir '+out_dir+'/coadd_map')
        conn = sq.connect(out_dir+'/mapfilelist_all.db')
        c = conn.cursor()
        c.execute('select * from mapfilelist;')
        self.runNum_arr = []; self.path_arr = []; self.pix_arr = []
        for ar in c:
            self.runNum_arr.append(ar[0])
            self.path_arr.append(ar[1])
            self.pix_arr.append(ar[2])
        c.close()
        self.nb = len(self.pix_arr)

        f_sh = open(out_dir+'/coadd_map/README', 'w')
        f_sh.write('select the data to coadd as  \n')
        f_sh.write('%s\n' % (out_dir+'/mapfilelist_all.db'))
        f_sh.write('%s\n' % (sqlite_command))
        f_sh.close()
        return self.runNum_arr, self.path_arr, self.pix_arr

    def read_coaddDB(self):
        print "  [main_coaddmaps.py:read_coaddDB] Read coadd database"
        if os.path.exists(xml_input["dir_simedmap"]+'/coadd_map.db'):
            conn = sq.connect(xml_input["dir_simedmap"]+'/coadd_map.db')
            c = conn.cursor()
            c.execute('select * from coadd_map_db')
            id=[]; sys_run_name=[]; outdir=[]; dir_coadd=[]; select_ces=[]
            for ar in c:
                id.append(ar[0])
                sys_run_name.append(ar[1])
                outdir.append(ar[2])
                dir_coadd.append(ar[3])
                select_ces.append(ar[4])
            c.close()
#            print select_ces
#            print np.array(select_ces)
#            print sqlite_command
#            ind = np.where( ((select_ces) == sqlite_command) )
#            print ind
#            print ind[0]
#            self.outdirs = outdir[ind[0]]
#            self.dir_coadds = dir_coadd[ind[0]]
            self.outdirs = outdir
            self.dir_coadds = dir_coadd

    def generate_coadd_PBS_SH(self):
        nb_qsub = 0
        out_file = []
        for i in range(0,nb_qsub+1):
            print dir_simulator+'/RunJob/'+runID 
            dir_gen = dir_simulator+'/RunJob/'+runID
            f_sh = open(dir_gen+'/BSUB_'+runID+'_'+self.dir_coadds[0]+'_'+name_extTOD+'.sh', 'w')
            if machine=="kekcc":
                f_sh.write('%s\n' % ('#!/bin/bash'))
                f_sh.write('%s\n' % ('bsub -q e -o '+dir_gen+'/SH_'+runID+'_'+self.dir_coadds[0]+'_'+name_extTOD+'.o '+dir_gen+'/SH_'+runID+'_'+self.dir_coadds[0]+'_'+name_extTOD+'.sh'))
                f_sh.write('%s\n' % ('wait'))
#                f_sh.write('%s\n' % ('bsub -q e -o '+dir_gen+'/SH_'+runID+'_'+self.dir_coadds[0]+'_TQU.o '+dir_gen+'/SH_'+runID+'_'+self.dir_coadds[0]+'_TQU.sh'))
#                f_sh.write('%s\n' % ('rm -f '+self.outdirs+'/coadd_map/'+self.dir_coadds+'/*.npy'))

#            if machine=='hopper':
#                f_sh.write('%s\n' % ('#PBS -S /bin/bash'))
#                f_sh.write('%s\n' % ('#PBS -l mppwidth='+str(n_core_tot)))
#                f_sh.write('%s\n' % ('#PBS -l mppwidth=1'))
#                f_sh.write('%s\n' % ('#PBS -N PB_'+self.dir_coadds+'_'+runID+'_'+str(i)))
#                f_sh.write('%s\n' % ('#PBS -A mp107'))
#            if machine=='carver':
#                f_sh.write('%s\n' % ('#PBS -S /bin/bash'))
#                f_sh.write('%s\n' % ('#PBS -l nodes='+str(n_cpu)+':ppn='+str(n_core)))
#                f_sh.write('%s\n' % ('#PBS -N PB_'+self.dir_coadds+'_'+runID+'_'+str(i)))
#                f_sh.write('%s\n' % ('#PBS -A mp107'))
#            if mode == 'debug': 
#                f_sh.write('%s\n' % ('#PBS -q debug'))
#                f_sh.write('%s\n' % ('#PBS -l walltime=01:00:00')) 
#            else:
#                f_sh.write('%s\n' % ('#PBS -q regular'))
#                f_sh.write('%s\n' % ('#PBS -l walltime=01:00:00')) 
#
#            f_sh.write('%s\n' % ('#PBS -o '+self.dir_coadds+'_'+runID+'.log'))
#            f_sh.write('%s\n' % ('#PBS -V'))
#            f_sh.write('%s\n' % ('#PBS -j eo'))
#            f_sh.write('%s\n' % ('cd $PBS_O_WORKDIR'))
    
            out_file.append(dir_gen+'/BSUB_'+runID+'_'+self.dir_coadds[0]+'_'+name_extTOD+'.sh')

#            if machine=='hopper':
#                f_sh.write('%s\n' % ('aprun -n 1 -d 1 -a xt '+dir_gen+'/SH_'+runID+'_'+self.dir_coadds+'.sh'))
#                f_sh.write('%s\n' % ('aprun -n 1 -d 1 -a xt '+dir_gen+'/SH_'+runID+'_'+self.dir_coadds+'_TQU.sh'))
#                f_sh.write('%s\n' % ('rm -f '+self.outdirs+'/coadd_map/'+self.dir_coadds+'/*.npy'))
#                f_sh.write('aprun -n '+str(n_cpu)+' -d '+str(n_core)+' -a xt '+dir_gen+'/SH_'+runID+'_coadd_'+str(i)+'.sh')
#                f_sh.write('%s\n' % ('python '+dir_simulator+'/src_viewer/viewmap.py '+out_dir+'/coadd_map/'+' '+'200e-6'+' 20e-6'))
#                f_sh.write('%s\n' % ('aprun -n 1 -d 1 -a xt python '+dir_simulator+'/src/fits_sep2oneTQU.py '  \
#                                         +self.outdirs+'/coadd_map/'+self.dir_coadds+'/mapT.fits ' \
#                                         +self.outdirs+'/coadd_map/'+self.dir_coadds+'/mapQ.fits ' \
#                                         +self.outdirs+'/coadd_map/'+self.dir_coadds+'/mapU.fits ' \
#                                         +self.outdirs+'/coadd_map/'+self.dir_coadds+'/mapTQU.fits'))
#            if machine=='carver':
#                f_sh.write(dir_gen+'/SH_'+runID+'_'+self.dir_coadds+'.sh')
#            f_sh.close()

            f_qsub = open(dir_gen+'/SH_'+runID+'_'+self.dir_coadds[0]+'_'+name_extTOD+'.sh', 'w')
            f_qsub.write('%s\n' % ('#!/bin/sh'))
#            f_qsub.write('%s\n' % ('python '+dir_simulator+'/src/run_coaddmaps.py '+out_dir+' '+nside+' "'+sqlite_command+'"& '))
            f_qsub.write('%s\n' % ('python '+dir_simulator+'/src/run_coaddmaps.py '+self.outdirs[0]+' '+self.dir_coadds[0]+'_'+name_extTOD+' '+nside+' "'+sqlite_command+'" '+xml_filename+' & '))
            f_qsub.write('%s\n' % ('wait'))
            f_qsub.close()

            f_TQU = open(dir_gen+'/SH_'+runID+'_'+self.dir_coadds[0]+'_'+name_extTOD+'_TQU.sh', 'w')
            f_TQU.write('%s\n' % ('#!/bin/sh'))
            f_TQU.write('%s\n' % ('python '+dir_simulator+'/src/fits_sep2oneTQU.py ' \
                                      +self.outdirs[0]+'/coadd_map/'+self.dir_coadds[0]+'_'+name_extTOD+'/mapT.fits ' \
                                      +self.outdirs[0]+'/coadd_map/'+self.dir_coadds[0]+'_'+name_extTOD+'/mapQ.fits ' \
                                      +self.outdirs[0]+'/coadd_map/'+self.dir_coadds[0]+'_'+name_extTOD+'/mapU.fits ' \
                                      +self.outdirs[0]+'/coadd_map/'+self.dir_coadds[0]+'_'+name_extTOD+'/mapTQU.fits')) 
            f_TQU.write('%s\n' % ('rm -f '+self.outdirs[0]+'/coadd_map/'+self.dir_coadds[0]+'_'+name_extTOD+'/*.npy'))  
            f_TQU.write('%s\n' % ('wait'))
            f_TQU.close()
            
            os.system('chmod 744 '+dir_gen+'/SH_'+runID+'_'+self.dir_coadds[0]+'_'+name_extTOD+'.sh')
            os.system('chmod 744 '+dir_gen+'/SH_'+runID+'_'+self.dir_coadds[0]+'_'+name_extTOD+'_TQU.sh')
            
        return out_file

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print ''
print ''
print '[main_coaddmaps.py] begin generate_coadd_PBS_SH()'
out = generate_coadd_PBS_SH()
out.read_DBmapfilelist()
out.read_coaddDB()
out_file = out.generate_coadd_PBS_SH()
print '[main_coaddmaps.py] end generate_coadd_PBS_SH()'
print ''

nb_qsub = len(out_file)
for i in range(nb_qsub):
#    os.system('qsub '+out_file[i])
    print ''
    print '[main_coaddmaps.py] bsub '+out_file[i]
    os.system('chmod 744 '+out_file[i])
    os.system(out_file[i])

print '[main_coaddmaps.py] end of the file'    

    
