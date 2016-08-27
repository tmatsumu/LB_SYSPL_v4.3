import numpy as np
import healpy as h
import pylab as py
import fileinput
import sys

pi = np.pi

######################################################################################################
#######################################################################################################
# read the focal plane data base
#  input: filename
#  output: ch, az, el, amp, sig_x, sig_y, amp, poleff, polang, wafer, pix, torb, flag
def read_FPDB2(filename):
    print ""
    print "[READ FPDB]: START reading "+filename
    
    ch = []
    az = []
    el = []
    sig_x = []
    sig_y = []
    amp = []
    polang = []
    poleff = []
    wafer = []
    pix = []
    torb = []
    flag = []
    i = 0
    for line in fileinput.input(filename):
        ar = line.split()
        if ((len(ar)>1) and (i>0)):
            ch.append(int(ar[0]))
            az.append(float(ar[1]))
            el.append(float(ar[2]))
            sig_x.append(float(ar[3]))
            sig_y.append(float(ar[4]))
            amp.append(float(ar[5]))
            polang.append(float(ar[6]))
            poleff.append(float(ar[7]))
            wafer.append(str(ar[8]))
            pix.append(int(ar[9]))
            torb.append(str(ar[10]))
            flag.append(int(ar[11]))
        i += 1
    print "[READ FPDB]: End reading "
    FPDBlist = {'ch':np.array(ch),'az':np.array(az),'el': np.array(el),
                'sig_x':np.array(sig_x),'sig_y':np.array(sig_y),
                'amp': np.array(amp),'poleff':np.array(poleff),'polang':np.array(polang),
                'wafer':np.array(wafer),'pix':np.array(pix),
                'torb': np.array(torb),'flag':np.array(flag)}
    return FPDBlist  

def gen_polvec(az,el,polang,amp):
    polang = polang/180.*pi
    az_out1 = az+amp*np.cos(polang)
    az_out2 = az-amp*np.cos(polang)
    el_out1 = el+amp*np.sin(polang)
    el_out2 = el-amp*np.sin(polang)
    return np.array([az_out1,az_out2]), np.array([el_out1,el_out2])

######################################################################################################################
######################################################################################################################

F = py.gcf()
DPI = F.get_dpi()
DefaultSize = F.get_size_inches()
F.set_size_inches( (DefaultSize[0], DefaultSize[1]) )
py.ioff()

###########################################################
#dir = '/global/homes/t/tmatsumu/develop/PBI/MapMaker/NaiveMapmaker/FPDB/'

dir = sys.argv[1]
fpdb_fname1 = sys.argv[2]
fpdb_list1 = read_FPDB2(dir+fpdb_fname1)

amp = float(sys.argv[3])
###########################################################
az1 = fpdb_list1['az']
el1 = fpdb_list1['el']
polang1 = fpdb_list1['polang']
flag1 = fpdb_list1['flag']
torb1 = fpdb_list1['torb']
ind_t1 = np.where((flag1 == 0) & ('t' == torb1))
ind_b1 = np.where((flag1 == 0) & ('b' == torb1))
ind1 = np.where(flag1 == 0)
ind1 = ind1[0]; ind_t1 = ind_t1[0]; ind_b1 = ind_b1[0]
nb1 = len(ind1)

###########################################################

print nb1

###########################################################
polvec_x1, polvec_y1 = gen_polvec(az1,el1,polang1,amp)
print len(polvec_x1[0,:])
print len(polvec_x1[:,0])

###########################################################
py.plot(az1[ind1],el1[ind1],'.')
for i in ind_t1:
    py.plot(polvec_x1[:,i],polvec_y1[:,i],'b')
for i in ind_b1:
    py.plot(polvec_x1[:,i],polvec_y1[:,i],'r')
py.title(fpdb_fname1)
py.xlabel('az [degs]')
py.ylabel('el [degs]')
py.savefig(dir+fpdb_fname1+'.png', dpi=DPI)

py.show()

