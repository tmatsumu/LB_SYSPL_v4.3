import numpy as np
import sqlite3 as sq
import sys
import os

def write_txt(fname,array,option_init=False):
    if option_init: f = open(fname, "w")
    f = open(fname, "a")
    f.write('%s \n' % (array))
    f.close()

def write_xml(fname,p):
    write_txt(fname,'<?xml version="1.0" encoding="UTF-8"?>',option_init=True)
    write_txt(fname,'<Mapmaking>')
    write_txt(fname,'  <log>')
    write_txt(fname,'    <runID>'+p.runID+'</runID>')
    write_txt(fname,'    <dir_simulator>'+p.dir_simulator+'</dir_simulator>')
    write_txt(fname,'    <debug>'+p.debug+'</debug>')
    write_txt(fname,'    <machine>'+p.machine+'</machine>')
    write_txt(fname,'    <silent>'+p.silent+'</silent>')
#    write_txt(fname,'    <core>'+str(p.core)+'</core>')
#    write_txt(fname,'    <cpu>'+str(p.cpu)+'</cpu>')
    write_txt(fname,'  </log>')
    write_txt(fname,' ')
#    write_txt(fname,'  <date>')
#    write_txt(fname,'    <date_i>'+str(p.date_i)+'</date_i>')
#    write_txt(fname,'    <date_f>'+str(p.date_f)+'</date_f>')
#    write_txt(fname,'  </date>')
#    write_txt(fname,' ')
    write_txt(fname,'  <basicpar>')
#    write_txt(fname,'    <fsample>'+str(p.fsample)+'</fsample>')
    write_txt(fname,'    <nside>'+str(p.nside)+'</nside>')
    write_txt(fname,'    <run_type>'+p.run_type+'</run_type>')
    write_txt(fname,'    <coord>'+p.coord+'</coord>')
    write_txt(fname,'    <pixelmapio>'+p.pixelmapio+'</pixelmapio>')
    write_txt(fname,'    <gen_tod>'+p.gen_tod+'</gen_tod>')
    write_txt(fname,'    <TQU>'+str(p.TQU)+'</TQU>')
    write_txt(fname,'  </basicpar>')
    write_txt(fname,' ')
    write_txt(fname,'  <simulations>')
    if (p.TQU=='T') | (p.TQU=='TQU'):
        write_txt(fname,'    <file_input_maps>'+p.file_input_maps+'</file_input_maps>')
    if p.TQU=='TQU2':
        write_txt(fname,'    <file_input_maps>'+p.file_input_maps+'</file_input_maps>')
        write_txt(fname,'    <file_input_maps2>'+p.file_input_maps2+'</file_input_maps2>')
    if p.TQU=='T_bandpassmismatch':
        write_txt(fname,'    <file_input_maps>'+p.file_input_maps+'</file_input_maps>')
        write_txt(fname,'    <file_input_map_dust>'+p.file_input_map_dust+'</file_input_map_dust>')
        write_txt(fname,'    <file_input_map_synch>'+p.file_input_map_synch+'</file_input_map_synch>')
        write_txt(fname,'    <file_input_bandpassmismatch>'+p.file_input_bandpassmismatch+'</file_input_bandpassmismatch>')
#    write_txt(fname,'    <file_input_ptg>'+p.file_input_ptg+'</file_input_ptg>')
    write_txt(fname,'    <file_input_noise>'+p.file_input_noise+'</file_input_noise>')
    write_txt(fname,'    <file_input_fpdb>'+p.file_input_fpdb+'</file_input_fpdb>')
    write_txt(fname,'    <file_input_muellermatrix>'+p.file_input_muellermatrix+'</file_input_muellermatrix>')
#    write_txt(fname,'    <patch_ra>'+str(p.patch_ra)+'</patch_ra>')
#    write_txt(fname,'    <patch_dec>'+str(p.patch_dec)+'</patch_dec>')
#    write_txt(fname,'    <patch_width>'+str(p.patch_width)+'</patch_width>')
    write_txt(fname,'  </simulations>')
    write_txt(fname,' ')
    write_txt(fname,'  <database>')
    write_txt(fname,'    <file_fpdb_mmin>'+p.file_fpdb_mmin+'</file_fpdb_mmin>')
    write_txt(fname,'    <file_muellermatrix>'+p.file_muellermatrix+'</file_muellermatrix>')
#    write_txt(fname,'    <file_relgain>'+p.file_relgain+'</file_relgain>')
#    write_txt(fname,'    <file_boloid>'+p.file_boloid+'</file_boloid>')
#    write_txt(fname,'  </database>')
#    write_txt(fname,' ')
#    write_txt(fname,'  <data_db>')
 #   write_txt(fname,'    <file_flag_pixel>'+p.file_flag_pixel+'</file_flag_pixel>')
    write_txt(fname,'    <db_ces>'+p.db_ces+'</db_ces>')
    write_txt(fname,'    <db_ces2>'+p.db_ces2+'</db_ces2>')
#    write_txt(fname,'    <sqlite_command_ces>'+p.sqlite_command_ces+'</sqlite_command_ces>')
    write_txt(fname,'    <db_gain>'+p.db_gain+'</db_gain>')
    write_txt(fname,'    <gain_type>'+p.gain_type+'</gain_type>')
    write_txt(fname,'    <gain_corr>'+p.gain_type+'</gain_corr>')
#    write_txt(fname,'    <db_HWPangles>'+p.db_ces+'</db_HWPangles>')
#    write_txt(fname,'    <RorL>'+p.db_ces+'</RorL>')
    write_txt(fname,'  </database>')
    write_txt(fname,' ')
    write_txt(fname,' <filtering_choice>')
    write_txt(fname,'    <filter_choice>'+str(p.filter_choice)+'</filter_choice>')
    write_txt(fname,'    <poly>'+str(p.poly)+'</poly>')
#    write_txt(fname,'    <file_noisefft>'+str(p.file_noisefft)+'</file_noisefft>')
    write_txt(fname,' </filtering_choice>')
    write_txt(fname,' ')
    write_txt(fname,'  <output_file>')
    write_txt(fname,'    <dir_simedmap>'+p.dir_simedmap+'</dir_simedmap>')
    write_txt(fname,'    <dir_combinedmap>'+p.dir_combinedmap+'</dir_combinedmap>')
    write_txt(fname,'    <Tn_map>map_Tn.fits</Tn_map>')
    write_txt(fname,'    <Td_map>map_Td.fits</Td_map>')
    write_txt(fname,'    <AA_map>map_AA.fits</AA_map>')
    write_txt(fname,'    <BB_map>map_BB.fits</BB_map>')
    write_txt(fname,'    <AB_map>map_AB.fits</AB_map>')
    write_txt(fname,'    <Ad_map>map_Ad.fits</Ad_map>')
    write_txt(fname,'    <Bd_map>map_Bd.fits</Bd_map>')
    write_txt(fname,'  </output_file>')
    write_txt(fname,' ')
#    write_txt(fname,'  <xpure>')
#    write_txt(fname,'    <dir_xpure>'+p.dir_xpure+'</dir_xpure>')
#    write_txt(fname,'    <radius>'+str(p.radius)+'</radius>')
#    write_txt(fname,'  </xpure>')
    write_txt(fname,' ')
    write_txt(fname,'</Mapmaking>')


class select_boloid():
    def __init__(self):
        self.file_boloid = 'tmp'
        self.file_flag_db = 'nosuchfile'

    def read_boloid(self,sq_command):
        conn = sq.connect(self.file_boloid)
        c = conn.cursor()
        c.execute(sq_command)
        boloid=[]; boloname=[]; pair=[]; pixelname=[]; wafer=[]; pixel=[]; board=[]; squid=[];
        for ar in c:
            boloid.append(int(ar[0]))
            boloname.append(str(ar[1]))
            pair.append(str(ar[2]))
            pixelname.append(str(ar[3]))
            wafer.append(str(ar[4]))
            pixel.append(int(ar[5]))
            board.append(int(ar[6]))
            squid.append(str(ar[7]))
        c.close()
        self.boloid_dict = {'boloid':boloid,'boloname':boloname,'pair':pair,'pixelname':pixelname,'wafer':wafer,
                            'pixel':pixel,'board':board,'squid':squid}
        return self.boloid_dict

    def display_all(self,db_dict):
        keys = db_dict.keys()
        num = len(db_dict[keys[0]])
        print keys
        for i in range(num):
            tmp = []
            for j in keys: tmp.append(db_dict[j][i])
            print tmp
        print keys

    def gen_boloflag(self,sq_command_select):
        sq_command = "select * from pb1_boloid"
        boloid = self.read_boloid(sq_command)
#        self.display_all(boloid)
        self.num_all = len(boloid['boloid'])

        boloid_selected = self.read_boloid(sq_command_select)
#        self.display_all(boloid_selected)
        num = len(boloid_selected['boloid'])
        print 'selected bolo #', num

        tmp = np.array(boloid['boloname'])
        tmp_selected = np.array(boloid_selected['boloname'])
        
        self.flag = np.ones(self.num_all,dtype='int')
        for i in range(0,num):
            ind = np.where(tmp_selected[i] == tmp)
            self.flag[ind[0]] = 0
        return self.flag

    def makeDB_flag(self):
        if os.path.exists(self.file_flag_db):
            print 'DB already exits'
            os.system('rm '+self.file_flag_db)
#            return -1
        else:
            conn = sq.connect(self.file_flag_db)
            c = conn.cursor()
            c.execute('create table boloselect_flag (boloid integer, flag integer)')
            for i in range(0,self.num_all):
                list_entries = ( i, int(self.flag[i]) )
                c.execute('insert into boloselect_flag values (?,?)',list_entries)
            conn.commit()
            c.close()

