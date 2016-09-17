import numpy as np
import fileinput
import sys
#import pylab as py
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as py

import healpy as h
import lib_mapmaker as lmm
import sqlite3 as sq
import time as time
import os
import ReadMapMakeXml as rxml

out_dir = sys.argv[1]
dir_coadd = sys.argv[2]
nside = int(sys.argv[3]) 
sqlite_command = sys.argv[4]
xml_filename = sys.argv[5]
xml_input = rxml.Get_Mapmake_Inputs(xml_filename)

class coadd():
    def __init__(self):
        self.fnameH = out_dir+'/coadd_map/'+dir_coadd+'/mapH'
        self.fnameT = out_dir+'/coadd_map/'+dir_coadd+'/mapT' 
        self.fnameQ = out_dir+'/coadd_map/'+dir_coadd+'/mapQ' 
        self.fnameU = out_dir+'/coadd_map/'+dir_coadd+'/mapU' 
        self.fnameS = out_dir+'/coadd_map/'+dir_coadd+'/mapS' 
        self.fnameW = out_dir+'/coadd_map/'+dir_coadd+'/mapW' 
        self.fnameM = out_dir+'/coadd_map/'+dir_coadd+'/mapM'
        self.fnameCN = out_dir+'/coadd_map/'+dir_coadd+'/mapCN'
        self.fnameA = out_dir+'/coadd_map/'+dir_coadd+'/mapA'
        self.fnameB = out_dir+'/coadd_map/'+dir_coadd+'/mapB'
        self.fnameAA = out_dir+'/coadd_map/'+dir_coadd+'/mapAA'
        self.fnameBB = out_dir+'/coadd_map/'+dir_coadd+'/mapBB'
        self.fnameAB = out_dir+'/coadd_map/'+dir_coadd+'/mapAB'
        self.fnameAd = out_dir+'/coadd_map/'+dir_coadd+'/mapAd'
        self.fnameBd = out_dir+'/coadd_map/'+dir_coadd+'/mapBd'
        self.fnamed = out_dir+'/coadd_map/'+dir_coadd+'/mapd'
        self.fnameTQU = out_dir+'/coadd_map/'+dir_coadd+'/mapTQU'
        self.nside = nside

##################################################################################
    def coadd_maps(self,path,idx):
        time0 = time.time()
        print ""
        print ""
        print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        print "++++++++ START: main_mapcombiner.py +++++++++++++++++++++++++++++++++++++++"
        print time.time()-time0

#        tmp = np.load(out_dir+'/ind_template.npz')        
#        keys = tmp.files
#        ind_part_out = tmp[keys[1]]
#        ind_band_out = tmp[keys[0]]
#        ind_min = tmp[keys[3]]
#        nside_out = tmp[keys[2]]
#        ind_part_tmp = np.where(ind_band_out != -1)
#        npix = len(ind_part_tmp[0])
        npix = h.nside2npix(nside)
       
        mapTn_s=np.zeros(npix); mapTd_s=np.zeros(npix);
        mapA_s=np.zeros(npix); mapB_s=np.zeros(npix);
        mapAA_s=np.zeros(npix); mapBB_s=np.zeros(npix); mapAB_s=np.zeros(npix);
        mapAd_s=np.zeros(npix); mapBd_s=np.zeros(npix); mapd_s=np.zeros(npix);
        mapH_s=np.zeros(npix)
        nb = len(idx)
           
        if xml_input['pixelmapio']=='Y': option_mapio_each=True
        if xml_input['pixelmapio']=='N': option_mapio_each=False

        print option_mapio_each
        print nb
        for i in range(0,nb):
            if option_mapio_each==True:
                print "[MAIN_MAPCOMBINER]: ", i, path[i], idx[i]-1, '/', nb
                print path[i]+'/map_'+str(idx[i]-1)+'.npz'
                out = np.load(path[i]+'/map_'+str(idx[i]-1)+'.npz')
            if option_mapio_each==False:
                print "[MAIN_MAPCOMBINER]: ", i, path[i], i, '/', nb
                print path[i]+'/map_all.npz'
                out = np.load(path[i]+'/map_all.npz')

            ind = out['ind']
            mapTn_s[ind] += out['In'] #mapTn
            mapTd_s[ind] += out['Id'] #mapTd
            mapA_s[ind] += out['A'] #mapA
            mapB_s[ind] += out['B'] #mapB
            mapAA_s[ind] += out['AA'] #mapAA
            mapBB_s[ind] += out['BB'] #mapBB
            mapAB_s[ind] += out['AB'] #mapAB
            mapAd_s[ind] += out['Ad'] #mapAd
            mapBd_s[ind] += out['Bd'] #mapBd
            mapd_s[ind] += out['d'] #mapBd
            mapH_s[ind] += out['H'] #mapBd

            out.close()
            del(ind); del(out)

        mapT=np.zeros(npix); 
        mapQ=np.zeros(npix); mapU=np.zeros(npix); mapS=np.zeros(npix)
        mapT = 2.*mapTn_s/mapTd_s
        print 'mapT ,n,d, nd', mapTn_s, mapTd_s, mapT

        S_params = False
        if S_params==True:
            det = mapH_s*mapA_s*mapB_s*mapAB_s + mapAA_s*mapBB_s + mapAB_s*mapA_s*mapB_s \
                - mapH_s*mapAB_s**2 - mapB_s**2*mapAA_s - mapA_s**2*mapBB_s
            mapS = 2./det*(  (mapAA_s*mapBB_s-mapAB_s**2)   *mapd_s - (mapA_s*mapBB_s-mapAB_s*mapB_s)*mapAd_s + (mapA_s*mapAB_s-mapB_s*mapAA_s)*mapBd_s )
            mapQ = 2./det*(  (mapAB_s*mapB_s-mapA_s*mapBB_s)*mapd_s - (mapB_s**2-mapBB_s*mapH_s)     *mapAd_s + (mapB_s*mapA_s-mapAB_s*mapH_s) *mapBd_s )
            mapU = 2./det*( -(mapAA_s*mapB_s-mapAB_s*mapA_s)*mapd_s + (mapA_s*mapB_s-mapAB_s*mapH_s) *mapAd_s - (mapA_s**2-mapAA_s*mapH_s)     *mapBd_s )
        if S_params==False:
            mapQ = 2./(mapAA_s*mapBB_s-mapAB_s**2) * (mapBB_s*mapAd_s - mapAB_s*mapBd_s)
            mapU = 2./(mapAA_s*mapBB_s-mapAB_s**2) * (-mapAB_s*mapAd_s + mapAA_s*mapBd_s)

        np.save(self.fnameT+'_nan',mapT)
        np.save(self.fnameQ+'_nan',mapQ)
        np.save(self.fnameU+'_nan',mapU)
        np.save(self.fnameS+'_nan',mapS)
        np.save(self.fnameW+'',mapTd_s)
        
#############################################################
        indT = np.isnan(mapT)
        indQ = np.isnan(mapQ)
        indU = np.isnan(mapU)
        indH = np.isnan(mapH_s)
        indQ_ = np.isinf(mapQ)
        indU_ = np.isinf(mapU)

        indT2 = np.where(indT == True)
        indQ2 = np.where(indQ == True)
        indU2 = np.where(indU == True)
        indH2 = np.where(indH == True)
        
        mapT[indT2] = 0.
        mapQ[indQ2] = 0.
        mapU[indU2] = 0.
        mapH_s[indH2] = 0.
        
        tmp_A_nan = np.isnan(mapA_s)
        tmp_B_nan = np.isnan(mapB_s)

        tmp_AA_nan = np.isnan(mapAA_s)
        tmp_BB_nan = np.isnan(mapBB_s)
        tmp_AB_nan = np.isnan(mapAB_s)
        tmp_Ad_nan = np.isnan(mapAd_s)
        tmp_Bd_nan = np.isnan(mapBd_s)
        
        tmp_A_inf = np.isinf(mapA_s)
        tmp_B_inf = np.isinf(mapB_s)
        tmp_AA_inf = np.isinf(mapAA_s)
        tmp_BB_inf = np.isinf(mapBB_s)
        tmp_AB_inf = np.isinf(mapAB_s)
        tmp_Ad_inf = np.isinf(mapAd_s)
        tmp_Bd_inf = np.isinf(mapBd_s)

        np.save(out_dir+'/coadd_map/'+dir_coadd+'/mapAsum',mapA_s)
        np.save(out_dir+'/coadd_map/'+dir_coadd+'/mapBsum',mapB_s)
        np.save(out_dir+'/coadd_map/'+dir_coadd+'/mapAAsum',mapAA_s)
        np.save(out_dir+'/coadd_map/'+dir_coadd+'/mapBBsum',mapBB_s)
        np.save(out_dir+'/coadd_map/'+dir_coadd+'/mapABsum',mapAB_s)
        np.save(out_dir+'/coadd_map/'+dir_coadd+'/mapAdsum',mapAd_s)
        np.save(out_dir+'/coadd_map/'+dir_coadd+'/mapBdsum',mapBd_s)
        np.save(out_dir+'/coadd_map/'+dir_coadd+'/mapdsum',mapd_s)
#        self.ind = np.where(self.M==1 & ((infQ==False) & (infU==False))) 
        #        ind = np.where(mapH_s > 0)
        ind = np.where( (mapH_s > 0) &  
                        ((tmp_A_nan==False) & (tmp_B_nan==False) & 
                         (tmp_AA_nan==False) & (tmp_BB_nan==False) & (tmp_AB_nan==False) &  
                         (tmp_Ad_nan==False) & (tmp_Bd_nan==False) &
                         (tmp_A_inf==False) & (tmp_B_inf==False) & 
                         (tmp_AA_inf==False) & (tmp_BB_inf==False) & (tmp_AB_inf==False) &  
                         (tmp_Ad_inf==False) & (tmp_Bd_inf==False) & 
                         (indQ_==False) & (indU_==False) ) )
        
        ind_inf = np.where( (tmp_A_nan==True) | (tmp_B_nan==True) |
                            (tmp_AA_nan==True) | (tmp_BB_nan==True) | (tmp_AB_nan==True) |
                            (tmp_Ad_nan==True) | (tmp_Bd_nan==True) |
                            (tmp_A_inf==True) | (tmp_B_inf==True) |
                            (tmp_AA_inf==True) | (tmp_BB_inf==True) | (tmp_AB_inf==True) |
                            (tmp_Ad_inf==True) | (tmp_Bd_inf==True) |
                            (indQ_==True) | (indU_==True) )
        
        mapQ[ind_inf[0]] = 0.
        mapU[ind_inf[0]] = 0.
        np.save(self.fnameT,mapT)
        np.save(self.fnameQ,mapQ)
        np.save(self.fnameU,mapU)
        np.save(self.fnameS,mapS)
        np.save(self.fnameH,mapH_s)

        h.write_map(self.fnameT+'.fits',mapT)
        h.write_map(self.fnameQ+'.fits',mapQ)
        h.write_map(self.fnameU+'.fits',mapU)
        h.write_map(self.fnameS+'.fits',mapS)
        h.write_map(self.fnameH+'.fits',mapH_s)

        mapTQU = [mapT,mapQ,mapU]
        h.write_map(self.fnameTQU+'.fits',mapTQU)
        del(mapTQU)

        h.mollview(mapT, title='T',max=np.abs(min(mapT)),min=-np.abs(min(mapT)))
        py.savefig(out_dir+'/png/'+dir_coadd+'/mapT.png')
        h.mollview(mapQ, title='Q',max=np.abs(min(mapT))*0.05, min=-np.abs(min(mapT))*0.05)
        py.savefig(out_dir+'/png/'+dir_coadd+'/mapQ.png')
        h.mollview(mapU, title='U',max=np.abs(min(mapT))*0.05, min=-np.abs(min(mapT))*0.05)
        py.savefig(out_dir+'/png/'+dir_coadd+'/mapU.png')
        h.mollview(mapS, title='S',max=np.abs(min(mapT))*0.05, min=-np.abs(min(mapT))*0.05)
        py.savefig(out_dir+'/png/'+dir_coadd+'/mapS.png')
        h.mollview(mapH_s, title='H',max=max(mapH_s),min=min(mapH_s))
        py.savefig(out_dir+'/png/'+dir_coadd+'/mapH.png')

        # this is the case for Q and U maps are zero.
        ind_zero = np.where(mapQ == 0)
        print len(mapQ), len(ind_zero[0])
        if (len(ind_zero[0]) == len(mapQ)):

            ind = np.where( (mapH_s>0) )
            mask = np.zeros(len(mapQ),dtype='int')
            mask[ind[0]] = 1

            np.save(self.fnameM,mask)
            h.write_map(self.fnameM+'.fits',mask)
            h.mollview(mask, title='M',max=max(mask),min=min(mask))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mask.png')

            np.save(self.fnameM+'_CNexcluded',mask)
            h.write_map(self.fnameM+'_CNexcluded.fits',mask)
            h.mollview(mask, title='CNexcluded',max=max(mask),min=min(mask))
            py.savefig(out_dir+'/png/'+dir_coadd+'/CNexcluded.png')

            np.save(self.fnameW+'_CNexcluded',mapTd_s)
            h.write_map(self.fnameW+'_CNexcluded.fits',mapTd_s)
            h.mollview(mapTd_s, title='mapTd_s',max=max(mapTd_s),min=min(mapTd_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapTd_s.png')

            np.save(self.fnameA+'_CNexcluded',mapA_s)
            h.write_map(self.fnameA+'_CNexcluded.fits',mapA_s)
            h.mollview(mapA_s, title='mapA_s',max=max(mapA_s),min=min(mapA_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapA_s.png')
            h.mollview(mapA_s/np.float_(mapH_s), title='mapA_s/H',max=max(mapA_s/np.float_(mapH_s)),min=min(mapA_s/np.float_(mapH_s)))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapA_s_normH.png')

            np.save(self.fnameB+'_CNexcluded',mapB_s)
            h.write_map(self.fnameB+'_CNexcluded.fits',mapB_s)
            h.mollview(mapB_s, title='mapB_s',max=max(mapB_s),min=min(mapB_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapB_s.png')
            h.mollview(mapB_s/np.float_(mapH_s), title='mapB_s/H',max=max(mapB_s/np.float_(mapH_s)),min=min(mapB_s/np.float_(mapH_s)))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapB_s_normH.png')

            np.save(self.fnameAA+'_CNexcluded',mapAA_s)
            h.write_map(self.fnameAA+'_CNexcluded.fits',mapAA_s)
            h.mollview(mapAA_s, title='mapAA_s',max=max(mapAA_s),min=min(mapAA_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapAA_s.png')
            h.mollview(mapAA_s/np.float_(mapH_s), title='mapAA_s/H',max=max(mapAA_s/np.float_(mapH_s)),min=min(mapAA_s/np.float_(mapH_s)))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapAA_s_normH.png')

            np.save(self.fnameBB+'_CNexcluded',mapBB_s)
            h.write_map(self.fnameBB+'_CNexcluded.fits',mapBB_s)
            h.mollview(mapBB_s, title='mapBB_s',max=max(mapBB_s),min=min(mapBB_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapBB_s.png')
            h.mollview(mapBB_s/np.float_(mapH_s), title='mapBB_s/H',max=max(mapBB_s/np.float_(mapH_s)),min=min(mapBB_s/np.float_(mapH_s)))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapBB_s_normH.png')

            np.save(self.fnameAB+'_CNexcluded',mapAB_s)
            h.write_map(self.fnameAB+'_CNexcluded.fits',mapAB_s)
            h.mollview(mapAB_s, title='mapAB_s',max=max(mapAB_s),min=min(mapAB_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapAB_s.png')
            h.mollview(mapAB_s/np.float_(mapH_s), title='mapAB_s/H',max=max(mapAB_s/np.float_(mapH_s)),min=min(mapAB_s/np.float_(mapH_s)))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapAB_s_normH.png')

            np.save(self.fnameAd+'_CNexcluded',mapAd_s)
            h.write_map(self.fnameAd+'_CNexcluded.fits',mapAd_s)
            h.mollview(mapAd_s, title='mapAd_s',max=max(mapAd_s),min=min(mapAd_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapAd_s.png')

            np.save(self.fnameBd+'_CNexcluded',mapBd_s)
            h.write_map(self.fnameBd+'_CNexcluded.fits',mapBd_s)
            h.mollview(mapBd_s, title='mapBd_s',max=max(mapBd_s),min=min(mapBd_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapBd_s.png')

            np.save(self.fnamed+'_CNexcluded',mapd_s)
            h.write_map(self.fnamed+'_CNexcluded.fits',mapd_s)
            h.mollview(mapd_s, title='mapd_s',max=max(mapd_s),min=min(mapd_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapd_s.png')

        cn = np.ones(npix)*1e10
        for i in ind[0]:
            cn[i] = lmm.Cal_ConditionNum_QUpix(mapAA_s[i],mapBB_s[i],mapAB_s[i])
        np.save(self.fnameCN,cn)
        h.write_map(self.fnameCN+'_CNexcluded.fits',cn)
        h.mollview(cn, title='CN',max=max(cn),min=min(cn))
        py.savefig(out_dir+'/png/'+dir_coadd+'/mapCN.png')

        if (len(ind_zero[0]) != len(mapQ)):
            print "POLARIZATION BASED MASK"
            mask = np.zeros(npix,int)
            mask[ind[0]] = 1
            np.save(self.fnameM,mask)
            ind = np.where((mapH_s > 0) & (np.abs(cn) < 1.e3) & 
                           ((tmp_A_nan==False) & (tmp_B_nan==False) &
                            (tmp_AA_nan==False) & (tmp_BB_nan==False) & (tmp_AB_nan==False) &
                            (tmp_Ad_nan==False) & (tmp_Bd_nan==False) & 
                            (tmp_A_inf==False) & (tmp_B_inf==False) &
                            (tmp_AA_inf==False) & (tmp_BB_inf==False) & (tmp_AB_inf==False) &
                            (tmp_Ad_inf==False) & (tmp_Bd_inf==False) &
                            (indQ_ == False) & (indU_ == False ) ))
            mask = np.zeros(npix,int)
            mask[ind[0]] = 1
            np.save(self.fnameM+'_CNexcluded_pol',mask)
            h.write_map(self.fnameM+'_CNexcluded_pol.fits',mask)
            h.mollview(mask, title='mask (pol)',max=max(mask),min=min(mask))
            py.savefig(out_dir+'/png/'+dir_coadd+'/maskpol.png')

            np.save(self.fnameW+'_CNexcluded',mapTd_s)
            h.write_map(self.fnameW+'_CNexcluded.fits',mapTd_s)
            h.mollview(mask, title='mapTd_s',max=max(mapTd_s),min=min(mapTd_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapTd_s.png')

            np.save(self.fnameA+'_CNexcluded',mapA_s)
            h.write_map(self.fnameA+'_CNexcluded.fits',mapA_s)
            h.mollview(mapA_s, title='mapA_s',max=max(mapA_s),min=min(mapA_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapA_s.png')
            h.mollview(mapA_s/np.float_(mapH_s), title='mapA_s/H',max=max(mapA_s/np.float_(mapH_s)),min=min(mapA_s/np.float_(mapH_s)))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapA_s_normH.png')

            np.save(self.fnameB+'_CNexcluded',mapB_s)
            h.write_map(self.fnameB+'_CNexcluded.fits',mapB_s)
            h.mollview(mapB_s, title='mapB_s',max=max(mapB_s),min=min(mapB_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapB_s.png')
            h.mollview(mapB_s/np.float_(mapH_s), title='mapB_s/H',max=max(mapB_s/np.float_(mapH_s)),min=min(mapB_s/np.float_(mapH_s)))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapB_s_normH.png')

            np.save(self.fnameAA+'_CNexcluded',mapAA_s)
            h.write_map(self.fnameAA+'_CNexcluded.fits',mapAA_s)
            h.mollview(mapAA_s, title='mapAA_s',max=max(mapAA_s),min=min(mapAA_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapAA_s.png')
            h.mollview(mapAA_s/np.float_(mapH_s), title='mapAA_s/H',max=max(mapAA_s/np.float_(mapH_s)),min=min(mapAA_s/np.float_(mapH_s)))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapAA_s_normH.png')

            np.save(self.fnameBB+'_CNexcluded',mapBB_s)
            h.write_map(self.fnameBB+'_CNexcluded.fits',mapBB_s)
            h.mollview(mapBB_s, title='mapBB_s',max=max(mapBB_s),min=min(mapBB_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapBB_s.png')
            h.mollview(mapBB_s/np.float_(mapH_s), title='mapBB_s/H',max=max(mapBB_s/np.float_(mapH_s)),min=min(mapBB_s/np.float_(mapH_s)))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapBB_s_normH.png')

            np.save(self.fnameAB+'_CNexcluded',mapAB_s)
            h.write_map(self.fnameAB+'_CNexcluded.fits',mapAB_s)
            h.mollview(mapAB_s, title='mapAB_s',max=max(mapAB_s),min=min(mapAB_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapAB_s.png')
            h.mollview(mapAB_s/np.float_(mapH_s), title='mapAB_s/H',max=max(mapAB_s/np.float_(mapH_s)),min=min(mapAB_s/np.float_(mapH_s)))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapAB_s_normH.png')

            np.save(self.fnameAd+'_CNexcluded',mapAd_s)
            h.write_map(self.fnameAd+'_CNexcluded.fits',mapAd_s)
            h.mollview(mapAd_s, title='mapAd_s',max=max(mapAd_s),min=min(mapAd_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapAd_s.png')

            np.save(self.fnameBd+'_CNexcluded',mapBd_s)
            h.write_map(self.fnameBd+'_CNexcluded.fits',mapBd_s)
            h.mollview(mapBd_s, title='mapBd_s',max=max(mapBd_s),min=min(mapBd_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapBd_s.png')

            np.save(self.fnamed+'_CNexcluded',mapd_s)
            h.write_map(self.fnamed+'_CNexcluded.fits',mapd_s)
            h.mollview(mapd_s, title='mapd_s',max=max(mapd_s),min=min(mapd_s))
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapd_s.png')

        del(mapA_s); del(mapB_s)
        del(mapAA_s); del(mapBB_s); del(mapAB_s); del(mapAd_s); del(mapBd_s); del(mapH_s)
        del(mapd_s)
        del(tmp_A_nan); del(tmp_B_nan)
        del(tmp_AA_nan);       del(tmp_BB_nan);       del(tmp_AB_nan);
        del(tmp_Ad_nan);       del(tmp_Bd_nan)
        del(tmp_A_inf);        del(tmp_B_inf)
        del(tmp_AA_inf);       del(tmp_BB_inf);       del(tmp_AB_inf);
        del(tmp_Ad_inf);       del(tmp_Bd_inf)
        del(mask);
        del(mapTd_s); del(mapTn_s)
        del(mapT);del(mapQ);del(mapU); del(mapS)

#############################################################
    def coadd_maps_part2full(self):
#        tmp = np.load(out_dir+'/ind_template.npz')
#        keys = tmp.files
#        ind_part_out = tmp[keys[1]]
#        ind_band_out = tmp[keys[0]]
#        ind_min = tmp[keys[3]]
#        nside_out = tmp[keys[2]]
#
#        print len(ind_part_out)
#        print len(ind_band_out)
#        print ind_min
#        print nside_out
#
        npix = h.nside2npix(nside)
#        npix_band = len(ind_band_out)
#
#        ind_part = np.where(ind_band_out != -1)

        fname = [self.fnameT, self.fnameQ, self.fnameU, self.fnameH, self.fnameM+'_CNexcluded']
        nb = len(fname)
        for i in range(0,nb):
            map_part = np.load(fname[i]+'_part.npy')
            map = np.zeros(npix)
#            map[ind_band_out[ind_part[0]]+ind_min] = map_part
            map[ind_part_out] = map_part
            h.write_map(fname[i]+'.fits',map)
            map=0

        fname = self.fnameH
        map_part = np.load(fname+'_part.npy')
        map = np.zeros(npix) 
        map[ind_part_out] = map_part 
        h.write_map(fname+'_inverseNoise.fits',map)
        map=0

#############################################################
    def make_figures(self,lon,lat,width, option_part=False):
#        os.popen('mkdir -p '+out_dir+'/png/'+dir_coadd)
#        os.popen('mv '+out_dir+'/*.npy '+out_dir+'/npy/')

        if option_part:
            map = h.read_map(self.fnameH+'.fits')
            h.mollview(map, title='H',max=max(map),min=0)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapH.png')
            
            h.cartview(map,coord='C',rot=[lon,lat],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],\
                           title='H', max=max(map),min=0)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapH_cart.png')


            map = h.read_map(self.fnameM+'_CNexcluded.fits')
            h.mollview(map, title='M CNexcluded',max=max(map),min=0)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapM_CNexcluded.png')
            
            h.cartview(map,coord='C',rot=[lon,lat],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],\
                           title='M CNexcluded', max=max(map),min=0)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapM_cart_CNexcluded.png')

            ind = np.where(map == 1)
            ind_non = np.where(map != 1)
            

            map = h.read_map(self.fnameT+'.fits')
            h.mollview(map, title='T', max=300e-6, min=-300e-6)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapT.png')
            
            h.cartview(map,coord='C',rot=[lon,lat],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],\
                           title='T', max=300e-6, min=-300e-6)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapT_cart.png')


            map = h.read_map(self.fnameQ+'.fits')
            h.mollview(map, title='Q',max=10e-6,min=-10e-6)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapQ.png')

            h.cartview(map,coord='C',rot=[lon,lat],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],\
                           title='Q', max=10e-6, min=-10e-6)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapQ_cart.png')


            map = h.read_map(self.fnameU+'.fits')
            h.mollview(map, title='U',max=10e-6,min=-10e-6)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapU.png')

            h.cartview(map,coord='C',rot=[lon,lat],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],\
                           title='U', max=10e-6,min=-10e-6)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapU_cart.png')


        if not option_part:
            map = h.read_map(self.fnameT+'.nan.fits')
            h.mollview(map, title='T')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapT.png')
            
            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='T')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/map_cart_T.png')
            
            
            map = h.read_map(self.fnameQ+'.nan.fits')
            h.mollview(map, title='Q')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapQ.png')
            
            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='Q')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/map_cart_Q.png')
            
            
            map = h.read_map(self.fnameU+'.nan.fits')
            h.mollview(map, title='U')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapU.png')
            
            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='U')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/map_cart_U.png')
            
            map = h.read_map(self.fnameH+'.fits')
            h.mollview(map, title='H')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapH.png')
            
            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='H')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/map_cart_H.png')
            
            
            map = h.read_map(self.fnameCN+'.fits')
            h.mollview(map, title='CN')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapCN.png')
            
            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='CN')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/map_cart_CN.png')
            
            
            map = h.read_map(self.fnameM+'.fits')
            h.mollview(map, title='M')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapM.png')
            
            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='M')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/map_cart_M.png')
            
            
            map = h.read_map(self.fnameM+'_CNexcluded.fits')
            h.mollview(map, title='M_CNexcluded')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapM_CNexcluded.png')
            
            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='M_CNexcluded')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/map_cart_M_CNexcluded.png')
            

            map = h.read_map(self.fnameW+'.fits')
            h.mollview(map, title='W')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapW.png')
            
            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='W')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/map_cart_W.png')
            
            
            map = h.read_map(self.fnameW+'_CNexcluded.fits')
            h.mollview(map, title='W_CNexcluded')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/mapW_CNexcluded.png')
            
            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='W_CNexcluded')
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(out_dir+'/png/'+dir_coadd+'/map_cart_W_CNexcluded.png')
#############################################################
    def make_figures_part(self):
        tmp = np.load(out_dir+'/ind_template.npz')
        ind_all = tmp[0]
        ind_all = 0
        ind_part = tmp[1]
        theta, phi = h.ang2pix(nside,ind_part)

        npix = len(ind_all)
        npix_part = len(ind_part)

        fname = [self.fnameT, self.fnameQ, self.fnameU, self.fnameH]
        nb = len(fname)
        for i in range(0,nb):
            map_part = np.load(fname[i]+'_part')
            map = np.zeros(npix)
            map[ind_part] = map_part
            h.write_map(fname[i]+'.fits')
            map=0

##########################################################################################################################
##########################################################################################################################
def read_DBmapfilelist():
    conn = sq.connect(out_dir+'/mapfilelist_all.db')
    c = conn.cursor()
#    c.execute('select * from mapfilelist')
    print '[read_DBmapfilelist]', sqlite_command
#    c.execute(sqlite_command)
    c.execute('select * from mapfilelist;')
#    runNum_arr = []; path_arr = []; pix_arr = []
#    id=[]; run_id=[]; run_subid=[]; dir_ptg=[]; outdir=[]; pix=[] 
    id=[]; juliantime=[]; path=[]
    for ar in c:
#        runNum_arr.append(ar[0])
#        path_arr.append(ar[1])
#        pix_arr.append(ar[2])
        id.append(ar[0])
        juliantime.append(ar[1])
        path.append(ar[2])
#        run_id.append(ar[1])
#        run_subid.append(ar[2])
#        dir_ptg.append(ar[3])
#        outdir.append(ar[4])
#        pix.append(ar[5])
    c.close()
#    nb = len(pix_arr)
#    return runNum_arr, path_arr, pix_arr 
#    return id, run_id, run_subid, dir_ptg, outdir, pix
    return id, juliantime, path 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print out_dir+'/coadd_map/'+dir_coadd
os.system('mkdir -p '+out_dir+'/coadd_map/'+dir_coadd)
os.system('mkdir -p '+out_dir+'/png/'+dir_coadd)
#id, run_id, run_subid, dir_ptg, outdir, pix = read_DBmapfilelist()
id, juliantime, path = read_DBmapfilelist()

coadd = coadd()
coadd.coadd_maps(path,id)

#coadd.coadd_maps_part2full()
#lat = float(xml_input['patch_dec'])
#lon = float(xml_input['patch_ra'])
#width = float(xml_input['patch_width'])
#coadd.make_figures(lon,lat, width, option_part=True)


#  LocalWords:  mapTd
