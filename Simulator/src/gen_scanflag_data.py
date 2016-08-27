from numpy import *
from healpy import *
from pylab import *
from ReadMapMakeXml import *
import lib_mapmaker as lmm
import glob

print ""
print ""
print "++++++++++++++++++++++++++++++++++++++++++++++++++++"
print "++++++++++++++++++++++++++++++++++++++++++++++++++++"
print "+++++++ GENERATE SCAN FLAG +++++++++++++++++++++++++"
print ''


xml_filename = sys.argv[1]
print "[Gen_scansetdirs]: START", xml_filename
xml_input = lmm.Gen_scansetdirs(xml_filename)
ptg_dirs = xml_input["dirs_pointing"]
ptg_dirs_part = xml_input["dirs_pointing_part"]
fname_ptg = xml_input["dir_InputPtg"]
date_i = xml_input["date_i"]
date_f = xml_input["date_f"]

fileNames = glob.glob(fname_ptg+"/observation_*")

nb = len(fileNames)
for i in range(0,nb): print fileNames[i]

nb_tmp = len(fname_ptg)
fileNames_date = np.zeros(nb,int)
for i in range(0,nb): fileNames_date[i] = int(fileNames[i][(nb_tmp+13):])
ind = np.where((fileNames_date >= date_i) & (fileNames_date <= date_f))


nb_obs = len(ind[0])
dir_obs = []
dir_obs_part = []
for i in range(0,nb_obs):
    dir_obs.append(fname_ptg+"/observation_"+str(fileNames_date[ind[0][i]]))
    dir_obs_part.append("observation_"+str(fileNames_date[ind[0][i]]))
    
file_input_pointing = []
file_input_pointing_part = []
for ii in range(0,nb_obs):
    fileNames_tmp = glob.glob(dir_obs[ii]+'/scan*')
    print fileNames_tmp
    nb = len(fileNames_tmp)
    #    print fileNames
    for j in range(0,nb):
        file_input_pointing.append(fileNames_tmp[j])
        file_input_pointing_part.append(fileNames_tmp[j][(nb_tmp+13):])

            
            

nb_scanset = len(ptg_dirs)
# +++++++++++++++++++++++++++++++++++++++++++++
# input files

nb = len(file_input_pointing)
for i in range(0,nb):
    filename = str(file_input_pointing[i])
    fileout = str(file_input_pointing[i])

    # +++++++++++++++++++++++++++++++++++++++++++++
    # read pointing.fits
    ptg = mrdfits(file_input_pointing[i]+'/pointing.fits', hdu=1)
    dec = ptg[1]
    
    # +++++++++++++++++++++++++++++++++++++++++++++
    # compute the # of samples from dec array 
    nb = len(dec)
    
    # +++++++++++++++++++++++++++++++++++++++++++++
    # generate an idex array [0,1,...,nb-1] [1,2,...,nb]
    idx1 = arange(0,nb-1) 
    idx2 = idx1 + 1
    
    # +++++++++++++++++++++++++++++++++++++++++++++
    # compute the transient
    #  and compute the index of the half scan
    del_dec = dec[idx2] - dec[idx1]
    deriv_max = max(del_dec)
    del_dec_norm = del_dec/deriv_max
    ind_p = where(del_dec_norm > 0.95)
    ind_n = where(del_dec_norm < -0.95)
    ind_p = ind_p[0]
    ind_n = ind_n[0]
    
    # +++++++++++++++++++++++++++++++++++++++++++++
    # indices of right and left scans
    start_right = ind_p+1
    end_right = ind_n
    
    start_left = ind_n+1
    end_left = ind_p
    
    nb_p = len(ind_p)
    nb_n = len(ind_n)
    nb = nb_p + nb_n # the total number of samples in two half scans (back and forth)
    
    idx_r = concatenate((start_right,end_right))
    idx_l = concatenate((start_left,end_left))
    
    if (idx_l[0] < idx_r[0]):
        idx_s = concatenate((start_left,start_right))
        idx_e = concatenate((end_left,end_right))
        nb = len(idx_s)
        sign = arange(0,nb)
        direc = (-1.)**(sign+1)
    if (idx_l[0] > idx_r[0]):
        idx_s = concatenate((start_right,start_left))
        idx_e = concatenate((end_right,end_left))
        nb = len(idx_s)
        sign = arange(0,nb)
        direc = (-1.)**sign

    ind_s = idx_s.argsort()
    ind_e = idx_e.argsort()

    print ''
    print 'Generate ', fileout+'/pointing_flag.txt'
    f = open(fileout+'/pointing_flag.txt', "w")
    f.write('%s\n' % ('scan idx, idx_start, idx_end, direction, flag')) 
    
    if (idx_s[ind_s[0]] > idx_e[ind_e[0]]):
        if ((nb % 2) == 1): nb = nb - 1
        #    print 'scan idx, idx_start, idx_end, direction, flag' 
        for i in range(0,nb-1):
            #        print i, idx_s[ind_s[i]], idx_e[ind_e[i+1]], int(direc[i+1]), 0
            f.write('%d   %d   %d   %d   %d\n' % (i, idx_s[ind_s[i]], idx_e[ind_e[i+1]], int(direc[i+1]), 0))
    if (idx_s[ind_s[0]] < idx_e[ind_e[0]]):
        #    print 'scan idx, idx_start, idx_end, direction, flag' 
        for i in range(0,nb):
            #        print i, idx_s[ind_s[i]], idx_e[ind_e[i]], int(direc[i]), 0
            f.write('%d   %d   %d   %d   %d\n' % (i, idx_s[ind_s[i]], idx_e[ind_e[i]], int(direc[i]), 0))
    f.close()

plot((dec-mean(dec))/max(dec-mean(dec)),'.')
plot(ind_n,(dec[ind_n]-mean(dec))/max(dec-mean(dec)),'.')
plot(ind_p,(dec[ind_p]-mean(dec))/max(dec-mean(dec)),'.')
plot(ind_n+1,(dec[ind_n+1]-mean(dec))/max(dec-mean(dec)),'*')
plot(ind_p+1,(dec[ind_p+1]-mean(dec))/max(dec-mean(dec)),'*')
show()
