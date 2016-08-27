import numpy as np
import fileinput
import sys
import pylab as py
import healpy as h
import lib_mapmaker as lmm
import sqlite3 as sq
import time as time
import os
import ReadMapMakeXml as rxml

#out_dir = sys.argv[1]
#dir_coadd = sys.argv[2]
#nside = int(sys.argv[3]) 
#sqlite_command = sys.argv[4]
#xml_filename = sys.argv[5]
#xml_input = rxml.Get_Mapmake_Inputs(xml_filename)

class map_util():
    def __init__(self):
        self.out_dir = './'
        self.dir_coadd = './'
        self.nside = 16
        self.sqlite_command = '.schema'
        self.xml_filename = 'test.xml'

    def def_filename(self):
        self.fnameH = self.out_dir+'/coadd_map/'+self.dir_coadd+'/mapH'
        self.fnameT = self.out_dir+'/coadd_map/'+self.dir_coadd+'/mapT' 
        self.fnameQ = self.out_dir+'/coadd_map/'+self.dir_coadd+'/mapQ' 
        self.fnameU = self.out_dir+'/coadd_map/'+self.dir_coadd+'/mapU' 
        self.fnameW = self.out_dir+'/coadd_map/'+self.dir_coadd+'/mapW' 
        self.fnameM = self.out_dir+'/coadd_map/'+self.dir_coadd+'/mapM'
        self.fnameCN = self.out_dir+'/coadd_map/'+self.dir_coadd+'/mapCN'

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

        tmp = np.load(self.out_dir+'/ind_template.npz')        
        keys = tmp.files
        ind_part_out = tmp[keys[1]]
        ind_band_out = tmp[keys[0]]
        ind_min = tmp[keys[3]]
#        nside_out = tmp[keys[2]]
        ind_part_tmp = np.where(ind_band_out != -1)
        npix = len(ind_part_tmp[0])
        
        mapTn_sum=np.zeros(npix); mapTd_sum=np.zeros(npix);
        mapAA_sum=np.zeros(npix); mapBB_sum=np.zeros(npix); mapAB_sum=np.zeros(npix);
        mapAd_sum=np.zeros(npix); mapBd_sum=np.zeros(npix)
        mapH_sum=np.zeros(npix)
        nb = len(idx)
           
        if xml_input['pixelmapio']=='Y': option_mapio_each=True
        if xml_input['pixelmapio']=='N': option_mapio_each=False

        for i in range(0,nb):
            if option_mapio_each==True:
                print "[MAIN_MAPCOMBINER]: ", i, path[i], idx[i]-1, '/', nb
                print path[i]+'/map_'+str(idx[i]-1)+'.npz'
                out = np.load(path[i]+'/map_'+str(idx[i]-1)+'.npz')
            if option_mapio_each==False:
                print "[MAIN_MAPCOMBINER]: ", i, path[i], i, '/', nb
                print path[i]+'/map_all.npz'
                out = np.load(path[i]+'/map_all.npz')
            keys = out.files
            ind = out[keys[1]]
            mapTn_sum[ind] += out[keys[0]] #mapTn
            mapTd_sum[ind] += out[keys[3]] #mapTd
            mapAA_sum[ind] += out[keys[2]] #mapAA
            mapBB_sum[ind] += out[keys[5]] #mapBB
            mapAB_sum[ind] += out[keys[4]] #mapAB
            mapAd_sum[ind] += out[keys[7]] #mapAd
            mapBd_sum[ind] += out[keys[6]] #mapBd
            mapH_sum[ind] += out[keys[8]] #mapBd

            out.close()
            del(ind); del(out); del(keys); 

        mapT=np.zeros(npix); 
        mapQ=np.zeros(npix); mapU=np.zeros(npix)
        mapT = 2.*mapTn_sum/mapTd_sum
#        print mapAA_sum,mapBB_sum, mapAB_sum**2, mapBB_sum,mapAd_sum, mapAB_sum,mapBd_sum
#        print mapAA_sum*mapBB_sum, mapAB_sum**2, mapBB_sum*mapAd_sum, mapAB_sum*mapBd_sum
        mapQ = 2./(mapAA_sum*mapBB_sum-mapAB_sum**2) * (mapBB_sum*mapAd_sum - mapAB_sum*mapBd_sum)
        mapU = 2./(mapAA_sum*mapBB_sum-mapAB_sum**2) * (-mapAB_sum*mapAd_sum + mapAA_sum*mapBd_sum)
        
        np.save(self.fnameT+'_nan_part',mapT)
        np.save(self.fnameQ+'_nan_part',mapQ)
        np.save(self.fnameU+'_nan_part',mapU)
        np.save(self.fnameW+'_part',mapTd_sum)
        
#############################################################
        indT = np.isnan(mapT)
        indQ = np.isnan(mapQ)
        indU = np.isnan(mapU)
        indH = np.isnan(mapH_sum)
        indQ_ = np.isinf(mapQ)
        indU_ = np.isinf(mapU)

        indT2 = np.where(indT == True)
        indQ2 = np.where(indQ == True)
        indU2 = np.where(indU == True)
        indH2 = np.where(indH == True)
        
        mapT[indT2] = 0.
        mapQ[indQ2] = 0.
        mapU[indU2] = 0.
        mapH_sum[indH2] = 0.
        
        tmp_AA_nan = np.isnan(mapAA_sum)
        tmp_BB_nan = np.isnan(mapBB_sum)
        tmp_AB_nan = np.isnan(mapAB_sum)
        tmp_AA_inf = np.isinf(mapAA_sum)
        tmp_BB_inf = np.isinf(mapBB_sum)
        tmp_AB_inf = np.isinf(mapAB_sum)

        np.save(self.out_dir+'/coadd_map/'+self.dir_coadd+'/mapAAsum_part',mapAA_sum)
        np.save(self.out_dir+'/coadd_map/'+self.dir_coadd+'/mapBBsum_part',mapBB_sum)
        np.save(self.out_dir+'/coadd_map/'+self.dir_coadd+'/mapABsum_part',mapAB_sum)
#        self.ind = np.where(self.M==1 & ((infQ==False) & (infU==False))) 
        #        ind = np.where(mapH_sum > 0)
        ind = np.where((mapH_sum > 0) &  
                       ((tmp_AA_nan==False) & (tmp_BB_nan==False) & (tmp_AB_nan==False) &  
                        (tmp_AA_inf==False) & (tmp_BB_inf==False) & (tmp_AB_inf==False) &  
                        (indQ_==False) & (indU_==False) ))

        ind_inf = np.where( (tmp_AA_nan==True) | (tmp_BB_nan==True) | (tmp_AB_nan==True) |
                            (tmp_AA_inf==True) | (tmp_BB_inf==True) | (tmp_AB_inf==True) |
                            (indQ_==True) | (indU_==True) )
        
        mapQ[ind_inf[0]] = 0.
        mapU[ind_inf[0]] = 0.
        np.save(self.fnameT+'_part',mapT)
        np.save(self.fnameQ+'_part',mapQ)
        np.save(self.fnameU+'_part',mapU)
        np.save(self.fnameH+'_part',mapH_sum)

        # the condition above doesn't work well when input Q and U maps are zeros.
        # this is the case for Q and U maps are zero.
        ind_zero = np.where(mapQ == 0)
        print len(mapQ), len(ind_zero[0])
        if (len(ind_zero[0]) == len(mapQ)):
            ind = np.where( (mapH_sum>0) )
            mask = np.zeros(len(mapQ),dtype='int')
            mask[ind[0]] = 1
            np.save(self.fnameM+'_part',mask)
            np.save(self.fnameM+'_CNexcluded_part',mask)
#            ind = np.where(np.abs(cn) < 1.e3)
#            mapTd_sum[ind[0]] = 0
            np.save(self.fnameW+'_CNexcluded_part',mapTd_sum)

        cn = np.ones(npix)*1e10
        for i in ind[0]:
            cn[i] = lmm.Cal_ConditionNum_QUpix(mapAA_sum[i],mapBB_sum[i],mapAB_sum[i])
        np.save(self.fnameCN+'_part',cn)

        if (len(ind_zero[0]) != len(mapQ)):
            print "POLARIZATION BASED MASK"
            mask = np.zeros(npix,int)
            mask[ind[0]] = 1
            np.save(self.fnameM+'_part',mask)
            ind = np.where((mapH_sum > 0) & (np.abs(cn) < 1.e3) & 
                           ((tmp_AA_nan==False) & (tmp_BB_nan==False) & (tmp_AB_nan==False) &
                            (tmp_AA_inf==False) & (tmp_BB_inf==False) & (tmp_AB_inf==False) &
                            (indQ_ == False) & (indU_ == False ) ))
            mask = np.zeros(npix,int)
            mask[ind[0]] = 1
            np.save(self.fnameM+'_CNexcluded_part',mask)
#        ind = np.where(np.abs(cn) < 1.e3)
#        mapTd_sum[ind[0]] = 0
            np.save(self.fnameW+'_CNexcluded_part',mapTd_sum)

        del(mapAA_sum); del(mapBB_sum); del(mapAB_sum); del(mapAd_sum); del(mapBd_sum); del(mapH_sum)
        del(tmp_AA_nan);       del(tmp_BB_nan);       del(tmp_AB_nan);
        del(tmp_AA_inf);       del(tmp_BB_inf);       del(tmp_AB_inf);
        del(mask);
        del(mapTd_sum); del(mapTn_sum)
        del(mapT);del(mapQ);del(mapU);

#############################################################
    def coadd_maps_part2full(self):
        tmp = np.load(self.out_dir+'/ind_template.npz')
        keys = tmp.files
        ind_part_out = tmp[keys[1]]
        ind_band_out = tmp[keys[0]]
        ind_min = tmp[keys[3]]
        nside_out = tmp[keys[2]]

        print len(ind_part_out)
        print len(ind_band_out)
        print ind_min
        print nside_out

        npix = h.nside2npix(nside_out)
        npix_band = len(ind_band_out)

        ind_part = np.where(ind_band_out != -1)

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
    def make_figures(self,lon,lat,width,Tmax,Pmax,option_part=False):
#        os.popen('mkdir -p '+self.out_dir+'/png/'+self.dir_coadd)
#        os.popen('mv '+self.out_dir+'/*.npy '+self.out_dir+'/npy/')

        if option_part:
            map = h.read_map(self.fnameH+'.fits')
            h.mollview(map, title='H',max=max(map),min=0)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapH.png')
            
            h.cartview(map,coord='C',rot=[lon,lat],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],\
                           title='H', max=max(map),min=0)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapH_cart.png')


            map = h.read_map(self.fnameM+'_CNexcluded.fits')
            h.mollview(map, title='M CNexcluded',max=max(map),min=0)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapM_CNexcluded.png')
            
            h.cartview(map,coord='C',rot=[lon,lat],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],\
                           title='M CNexcluded', max=max(map),min=0)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapM_cart_CNexcluded.png')

            ind = np.where(map == 1)
            ind_non = np.where(map != 1)
            

            map = h.read_map(self.fnameT+'.fits')
            h.mollview(map, title='T', max=Tmax, min=-Tmax)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapT.png')
            
            h.cartview(map,coord='C',rot=[lon,lat],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],\
                           title='T', max=Tmax, min=-Tmax)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapT_cart.png')


            map = h.read_map(self.fnameQ+'.fits')
            h.mollview(map, title='Q',max=Pmax, min=-Pmax)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapQ.png')

            h.cartview(map,coord='C',rot=[lon,lat],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],\
                           title='Q', max=Pmax, min=-Pmax)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapQ_cart.png')


            map = h.read_map(self.fnameU+'.fits')
            h.mollview(map, title='U',max=Pmax,min=-Pmax)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapU.png')

            h.cartview(map,coord='C',rot=[lon,lat],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],\
                           title='U', max=Pmax,min=-Pmax)
            h.graticule(dpar=10,dmer=10,coord='C')
            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapU_cart.png')


#        if not option_part:
#            map = h.read_map(self.fnameT+'.nan.fits')
#            h.mollview(map, title='T')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapT.png')
#            
#            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='T')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/map_cart_T.png')
#            
#            
#            map = h.read_map(self.fnameQ+'.nan.fits')
#            h.mollview(map, title='Q')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapQ.png')
#            
#            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='Q')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/map_cart_Q.png')
#            
#            
#            map = h.read_map(self.fnameU+'.nan.fits')
#            h.mollview(map, title='U')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapU.png')
#            
#            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='U')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/map_cart_U.png')
#            
#            map = h.read_map(self.fnameH+'.fits')
#            h.mollview(map, title='H')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapH.png')
#            
#            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='H')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/map_cart_H.png')
#            
#            
#            map = h.read_map(self.fnameCN+'.fits')
#            h.mollview(map, title='CN')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapCN.png')
#            
#            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='CN')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/map_cart_CN.png')
#            
#            
#            map = h.read_map(self.fnameM+'.fits')
#            h.mollview(map, title='M')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapM.png')
#            
#            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='M')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/map_cart_M.png')
#            
#            
#            map = h.read_map(self.fnameM+'_CNexcluded.fits')
#            h.mollview(map, title='M_CNexcluded')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapM_CNexcluded.png')
#            
#            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='M_CNexcluded')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/map_cart_M_CNexcluded.png')
#
#
#            map = h.read_map(self.fnameW+'.fits')
#            h.mollview(map, title='W')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapW.png')
#            
#            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='W')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/map_cart_W.png')
#            
#            
#            map = h.read_map(self.fnameW+'_CNexcluded.fits')
#            h.mollview(map, title='W_CNexcluded')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/mapW_CNexcluded.png')
#            
#            h.cartview(map,coord='C',rot=[ra,lon],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='W_CNexcluded')
#            h.graticule(dpar=10,dmer=10,coord='C')
#            py.savefig(self.out_dir+'/png/'+self.dir_coadd+'/map_cart_W_CNexcluded.png')
#############################################################
    def make_figures_part(self):
        tmp = np.load(self.out_dir+'/ind_template.npz')
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
    conn = sq.connect(self.out_dir+'/mapfilelist_all.db')
    c = conn.cursor()
#    c.execute('select * from mapfilelist')
    print '[read_DBmapfilelist]', sqlite_command
    c.execute(sqlite_command)
#    runNum_arr = []; path_arr = []; pix_arr = []
    id=[]; run_id=[]; run_subid=[]; dir_ptg=[]; outdir=[]; pix=[] 
    for ar in c:
#        runNum_arr.append(ar[0])
#        path_arr.append(ar[1])
#        pix_arr.append(ar[2])
        id.append(ar[0])
        run_id.append(ar[1])
        run_subid.append(ar[2])
        dir_ptg.append(ar[3])
        outdir.append(ar[4])
        pix.append(ar[5])
    c.close()
#    nb = len(pix_arr)
#    return runNum_arr, path_arr, pix_arr 
    return id, run_id, run_subid, dir_ptg, outdir, pix

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
