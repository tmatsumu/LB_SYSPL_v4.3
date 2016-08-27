import numpy as np
import fileinput
import healpy as h
import sys
import pylab as py
import lib_mapmaker as lmm
from ReadMapMakeXml import *
import time as time
import os

def read_fitslist(filename):
    dirs = []
    idx = []
    for line in fileinput.input(filename):
        ar = line.split()
        if (len(ar)>1):
            dirs.append(ar[0])
            idx.append(ar[1])
    print "[READ fitslist]: End reading "+filename
    return dirs, np.array(idx)


outdir = sys.argv[1]

dirs, idx = read_fitslist(outdir+"/fitsfilelist.txt")
nb = len(idx)

f = open(outdir+"/fitsfilelist_7days.txt", "w")
for i in range(0,nb):
#    if ('20111101' in dirs[i]):
#    if (('20111101' in dirs[i]) | ('20111102' in dirs[i])):
#    if (('20111101' in dirs[i]) | ('20111102' in dirs[i]) | ('20111103' in dirs[i])):
#    if (('20111101' in dirs[i]) | ('20111102' in dirs[i]) | ('20111103' in dirs[i]) | ('20111104' in dirs[i])):
#    if (('20111101' in dirs[i]) | ('20111102' in dirs[i]) | ('20111103' in dirs[i]) | ('20111104' in dirs[i]) | ('20111105' in dirs[i])):
#    if (('20111101' in dirs[i]) | ('20111102' in dirs[i]) | ('20111103' in dirs[i]) | ('20111104' in dirs[i]) | ('20111105' in dirs[i]) | ('20111106' in dirs[i])):
    if (('20111101' in dirs[i]) | ('20111102' in dirs[i]) | ('20111103' in dirs[i]) | ('20111104' in dirs[i]) | ('20111105' in dirs[i]) | ('20111106' in dirs[i]) | ('20111107' in dirs[i])):
        f.write('%s %d \n' % (dirs[i], int(idx[i])))
f.close()
        
