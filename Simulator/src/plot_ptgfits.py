import pylab as py
import numpy as ny
import lib_mapmaker as lib
import sys


'''
python plot_ptgfits.py runID date/scan0 

'''

dir = '/global/scratch/sd/tmatsumu/ScanStrategy/RunLog/'
runID = sys.argv[1]
date = sys.argv[2]
num1 = int(sys.argv[3])
num2 = int(sys.argv[4])

ptg = lib.read_ptg(dir+runID+'/'+date+'/pointing.fits')

py.subplot(511)
py.plot(ptg['flag'][num1:num2])
py.ylim([-.2,1.2])
py.subplot(512)
py.plot(ptg['ra'][num1:num2])
py.subplot(513)
py.plot(ptg['dec'][num1:num2])
py.subplot(514)
py.plot(ptg['pa'][num1:num2])
py.subplot(515)
py.plot(ptg['hwp'][num1:num2])
py.show()
