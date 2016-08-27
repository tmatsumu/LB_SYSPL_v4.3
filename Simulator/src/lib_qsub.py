import numpy as np
import os
import sys

def log_round(n):
    x = np.log(n)/np.log(2.)
    if (x - int(x)) == 0: pass
    if (x - int(x)) > 0: x = int(x)+1
    return 2**int(x)

class gen_qsub():
    def __def__(self):
#        self.i_scanset = 0
#        self.i_scanset_arr = np.array([1,2])
        self.file_inputnpy = ['tmp','tmp']
        self.dir_simulator = 'tmp'
        self.runID = 'tmp'
        self.mode = ''
        self.dir_out = 'tmp'
        self.machine = 'hopper'
        self.debug = ''

#    def gen_qsub_mm(self,n_core,n_cpu):
    def gen_qsub_mm(self):
        nb = len(self.file_inputnpy)

        out_file = []

        for i in range(0,nb):
            dir_gen = self.dir_simulator+"/RunJob/"+self.runID
            print '[lib_qsub.py:gen_qsub:'+self.machine+']', dir_gen
            if not os.path.exists(dir_gen):
                os.system('mkdir -p '+dir_gen)
                print '[lib_qsub.py:gen_qsub_mm]', 'the directory does not exist', dir_gen
                print '[lib_qsub.py:gen_qsub_mm]', 'creating...'

            f_sh = open(dir_gen+"/PBS_mm_"+self.runID+"_"+str(i)+".sh", "w")
            print '[lib_qsub.py:gen_qsub_mm] Generating PBS file', i, '/',nb, dir_gen+"/PBS_mm_"+self.runID+"_"+str(i)+".sh"

        if self.machine=='shell':
            f_shell = open(dir_gen+"/SH_"+self.runID+"_all.sh", "w")
            for i in range(0,nb):
                f_shell.write(''+self.dir_simulator+'/RunJob/'+self.runID+'/SH_'+self.runID+'_'+str(i)+'.sh > '+self.dir_simulator+'/RunJob/'+self.runID+'/log.txt &')
            f_shell.close()

        if self.machine=='kekcc':
            out_file.append(dir_gen+"/SH_"+self.runID+"_all.sh")
            f_shell = open(dir_gen+"/SH_"+self.runID+"_all.sh", "w")
            f_shell.write('#!/bin/sh \n')
            f_shell.write('')
            for i in range(0,nb):
                f_shell.write('bsub -q s -o '+self.dir_simulator+'/RunJob/'+self.runID+'/SH_'+self.runID+'_'+str(i)+'.o '+self.dir_simulator+'/RunJob/'+self.runID+'/SH_'+self.runID+'_'+str(i)+'.sh \n')
            f_shell.close()
            os.system('chmod 744 '+self.dir_simulator+'/RunJob/'+self.runID+'/SH_'+self.runID+'_'+str(i)+'.sh')
            os.system('chmod 744 '+self.dir_simulator+'/RunJob/'+self.runID+'/SH_'+self.runID+'_all.sh')

        for i in range(0,nb):
            f_qsub = open(dir_gen+"/SH_"+self.runID+"_"+str(i)+".sh", "w")
            print ''
            print '[lib_qsub.py:gen_qsub_mm] Generating SH file', i, '/',nb, dir_gen+"/SH_"+self.runID+"_"+str(i)+".sh"
            f_qsub.write('%s\n' % ('#!/bin/sh'))
            f_qsub.write('%s\n' % ('python '+self.dir_simulator+'/src/run_mapmaker.py '+self.file_inputnpy[i]+'.npy &'))
            f_qsub.write('%s\n' % ('wait'))
            f_qsub.close()
            os.system('chmod 744 '+self.dir_simulator+'/RunJob/'+self.runID+'/SH_'+self.runID+'_'+str(i)+'.sh')
        return out_file

class gen_mappng():
    def _init_(self):
        self.runID = 'test'
        self.dir_simulator = './'
        self.sqlite_command = '.schema'
        self.xml_filename = 'test'
        self.Tmax = Tmax
        self.Pmax = Pmax

    def gen_SH(self):
        dir_gen = self.dir_simulator+"/RunJob/"+self.runID
        dir_src = self.dir_simulator+"/src/"
        filename_SH = dir_gen+"/SH_gen_mappng_"+self.runID+".sh"
        f_sh = open(filename_SH, "w")
        f_sh.write('%s\n' % ('#!/bin/sh'))
        f_sh.write('%s\n' % ('python '+dir_src+'/gen_mappng.py '+self.xml_filename+' "'+self.sqlite_command+'" '+self.Tmax+' '+self.Pmax+' &'))
        f_sh.close()
        os.system('chmod 744 '+filename_SH)
        return filename_SH

    def gen_PBS(self):
        dir_gen = self.dir_simulator+"/RunJob/"+self.runID
        filename_PBS = dir_gen+"/PBS_gen_mappng_"+self.runID+".sh"
        f_qsub = open(filename_PBS, "w")
        f_qsub.write('%s\n' % ('#PBS -S /bin/bash'))
        f_qsub.write('%s\n' % ('#PBS -l mppwidth=1'))
        f_qsub.write('%s\n' % ('#PBS -N mappng_'+self.runID+'_PB'))
        f_qsub.write('%s\n' % ('#PBS -A mp107'))
        f_qsub.write('%s\n' % ('#PBS -o mappng_'+self.runID+'.log'))
        f_qsub.write('%s\n' % ('#PBS -q regular'))
        f_qsub.write('%s\n' % ('#PBS -l walltime=00:10:00')) 
        f_qsub.write('%s\n' % ('#PBS -V'))
        f_qsub.write('%s\n' % ('#PBS -j eo'))
        f_qsub.write('%s\n' % ('cd $PBS_O_WORKDIR'))
#        f_qsub.write('%s\n' % ('aprun -n 1 -N 1 -d 1 -a xt '+filename_PBS))
        f_qsub.write('%s\n' % ('aprun -n 1 -d 1 -a xt '+filename_PBS))
        f_qsub.close()
        os.system('chmod 744 '+filename_PBS)
        return filename_PBS

