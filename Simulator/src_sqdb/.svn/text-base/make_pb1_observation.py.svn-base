import numpy as np
import sqlite3 as sq
import fileinput
import sys
import os

def read_table(filename):
    import fileinput
    arr1=[]; arr2=[]; arr3=[]; arr4=[];
    arr5=[]; arr6=[]; arr7=[]; arr8=[];
    arr9=[]; arr10=[]; arr11=[]; arr12=[];
    arr13=[]; arr14=[]; arr15=[]; arr16=[];
    filelines = fileinput.input(filename)
    i=0
    for line in filelines:
        if i>=0:
            ar = line.split()
            arr1.append(int(ar[0])); arr2.append(float(ar[1])); arr3.append(float(ar[2])); arr4.append(float(ar[3]))
            arr5.append(float(ar[4])); arr6.append(int(ar[5])); arr7.append(int(ar[6])); arr8.append(int(ar[7]))
            arr9.append(float(ar[8])); arr10.append(float(ar[9])); arr11.append(float(ar[10])); arr12.append(float(ar[11]))
            arr13.append(float(ar[12])); arr14.append(float(ar[13])); arr15.append(float(ar[14])); arr16.append(float(ar[15]))
        i+=1
    return np.array(arr1),np.array(arr2),np.array(arr3),np.array(arr4),np.array(arr5),np.array(arr6),np.array(arr7),np.array(arr8),np.array(arr9),np.array(arr10),np.array(arr11),np.array(arr12),np.array(arr13),np.array(arr14),np.array(arr15),np.array(arr16)

class construct_table():
    def __init__(self):
        self.db_name = 'tmp.db'
        self.dir_ptg = {}
        self.first_mjd = {} 
        self.last_mjd = {} 
        self.num = 1
        
    def make_CESdb(self):
        if os.path.exists(self.db_name):
            print 'DB already exits'
            return -1
        else:
            conn = sq.connect(self.db_name)
            c = conn.cursor()
            c.execute('create table pb1_observation (id integer, run_id integer, run_subid integer, dir_ptg text, first_mjd real, last_mjd real)')
            for i in range(0,self.num): 
                list_entries = ( int(i),
                                 int(i/5),
                                 int(i%5),
                                 self.dir_ptg[i],
                                 self.first_mjd[i],
                                 self.last_mjd[i])
                c.execute('insert into pb1_observation values (?,?,?,?,?,?)',list_entries)
            conn.commit()
            c.close()

    def make_Flag_pixel(self):
        num = 1511
        if os.path.exists(self.db_name):
            print 'DB already exits'
            return -1
        else:
            conn = sq.connect(self.db_name)
            c = conn.cursor()
            c.execute('create table Flag_pixel (flag integer)')
            for i in range(0,self.num): 
                list_entries = ( int(0))
                c.execute('insert into CESdb values (?)',list_entries)
            conn.commit()
            c.close()

class read_DB():
    def __init__(self):
        self.filename = 'tmp.db'
        self.sq_command = 'select * from pb1_observation'

    def read_pb1_observation_fake(self):
        conn = sq.connect(self.filename)
        c = conn.cursor()
        c.execute(self.sq_command)
        id=[];run_id=[];run_subid=[];dir_ptg=[];first_mjd=[];last_mjd=[]
        for ar in c:
            id.append(int(ar[0]))
            run_id.append(int(ar[1]))
            run_subid.append(int(ar[2]))
            dir_ptg.append(str(ar[3]))
            first_mjd.append(float(ar[4]))
            last_mjd.append(float(ar[5]))
        c.close()
        self.CESdb = {'id':id,'run_id':run_id,'run_subid':run_subid,'dir_ptg':dir_ptg,
                           'first_mjd':first_mjd,'last_mjd':last_mjd}
        return self.CESdb

    def display_all(self,db_dict):
        keys = db_dict.keys()
        num = len(db_dict[keys[0]])
        print keys
        for i in range(num):
            tmp = []
            for j in keys: tmp.append(db_dict[j][i])
            print tmp
        print keys

def io_example_pb1_observation_fake(filename):
    read = read_DB()
    read.filename = filename
    db_new = read.read_pb1_observation_fake()
    read.display_all(db_new)


####################################################################################################################################
#dir_ptgdata = '/home/tmatsumu/data_sim/PBI/ScanStrategy/RunLog/1333707630.59_id20120406_candidate_100Hz_lst23_7days_noHWP_el40'
#dir_ptgdata = '/project/projectdirs/polar/user/tmatsumu/sim/ScanStrategy/RunLog/1340002291.46_id20120618_candidate_100Hz_lst23_7days_noHWP_el40'
#dir_ptgdata = '/project/projectdirs/polar/pipeline/pb1_ntp/ExampleData/1340002291.46_id20120618_candidate_100Hz_lst23_7days_noHWP_el40'
#dir_ptgdata = '/project/projectdirs/polar/user/tmatsumu/sim/ScanStrategy/RunLog/1343282418.39_id20120726_candidate_100Hz_lst23_7days_noHWP_el40'
#dir_ptgdata = '/project/projectdirs/polar/user/tmatsumu/sim/ScanStrategy/RunLog/1357279950.23_id20130104_10Hz_lst23_7days_HWP11p25_el40'
dir_ptgdata = '/project/projectdirs/polar/user/tmatsumu/sim/ScanStrategy/RunLog/1357354562.46_id20130104_50Hz_lst23_7days_HWP11p25_el40'
filename_table = dir_ptgdata+'/table.txt'
a1,a2,period,a4,a5,a6,a7,a8,a9,a10,a11,first_mjd,a13,a14,a15,a16=read_table(filename_table)

dir_date = ['observation_20120701','observation_20120702','observation_20120703','observation_20120704','observation_20120705','observation_20120706','observation_20120707']
#dir_date = ['observation_20120401','observation_20120402','observation_20120403','observation_20120404','observation_20120405','observation_20120406','observation_20120407','observation_20120408','observation_20120409','observation_20120410','observation_20120411','observation_20120412','observation_20120413']
dir_ces = ['scan0','scan1','scan2','scan3','scan4']

num_date = len(dir_date)
num_ces = len(dir_ces)

dir_out = []
for i in range(0,num_date):
    for j in range(0,num_ces):
        dir_out.append(dir_ptgdata+'/'+dir_date[i]+'/'+dir_ces[j])

print dir_out
db_filename = './pb1_observation_test15x15.db'
make_db = construct_table()
make_db.db_name = db_filename
make_db.dir_ptg = dir_out
make_db.first_mjd = first_mjd
make_db.last_mjd = first_mjd+period/3600./24.
make_db.num = len(dir_out)
print len(dir_out)
make_db.make_CESdb()

io_example_pb1_observation_fake(db_filename)
