import numpy as np
import pylab as py
import healpy as h
import sqlite3 as sq
from scipy.optimize import fmin
import ReadMapMakeXml as rxml
import sys
import os

xml_filename = sys.argv[1]
xml_input = rxml.Get_Mapmake_Inputs(xml_filename)
dir_out = xml_input["dir_simedmap"]
print '[cal_map_stats]', dir_out
sqlite_command = sys.argv[2]

def read_coaddDB(dir_simedmap):
    if os.path.exists(dir_simedmap+'/coadd_map.db'):
        conn = sq.connect(dir_simedmap+'/coadd_map.db')
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
    else:
        print '[read_coaddDB] no '+dir_simedmap+'/coadd_map.db'
    return sys_run_name, outdir, dir_coadd, select_ces

                                                      
class cal_map_stats():
    def __init__(self):
        self.outdir = './tmp'
        self.dir_coadd = './tmp'

    def read_coaddmaps(self):
        self.T = h.read_map(self.outdir+'/coadd_map/'+self.dir_coadd+'/mapT.fits')
        self.Q = h.read_map(self.outdir+'/coadd_map/'+self.dir_coadd+'/mapQ.fits')
        self.U = h.read_map(self.outdir+'/coadd_map/'+self.dir_coadd+'/mapU.fits')
        self.M = h.read_map(self.outdir+'/coadd_map/'+self.dir_coadd+'/mapM_CNexcluded.fits')

    def cal_stats(self):
        self.nbpix = len(self.M)
        infQ = np.isinf(self.Q)
        infU = np.isinf(self.U)

        self.ind = np.where(self.M==1 & ((infQ==False) & (infU==False)))
        meanT = np.mean(self.T[self.ind[0]])
        meanQ = np.mean(self.Q[self.ind[0]])
        meanU = np.mean(self.U[self.ind[0]])

        for i in range(0,5):
            if i ==0:
                self.stdT = np.std(self.T[self.ind[0]])
                self.stdQ = np.std(self.Q[self.ind[0]])
                self.stdU = np.std(self.U[self.ind[0]])
            else:
                self.stdT = np.std(self.T[self.indT[0]])
                self.stdQ = np.std(self.Q[self.indQ[0]])
                self.stdU = np.std(self.U[self.indU[0]])
                
            self.indT = np.where((self.M==1) & ((infQ==False) & (infU==False)) & (np.abs(self.T)<20.*self.stdT))
            self.indQ = np.where((self.M==1) & ((infQ==False) & (infU==False)) & (np.abs(self.Q)<20.*self.stdQ))
            self.indU = np.where((self.M==1) & ((infQ==False) & (infU==False)) & (np.abs(self.U)<20.*self.stdU))
            
        map_stats = [meanT,meanQ,meanU,self.stdT,self.stdQ,self.stdU]

#        self.indT_ol = np.where(self.T[self.ind[0]] > 3.*self.stdT)
#        self.indQ_ol = np.where(self.Q[self.ind[0]] > 3.*self.stdQ)
#        self.indU_ol = np.where(self.U[ind[0]] > 3.*self.stdU)

        return map_stats

    def gen_plots(self,lon,lat,width):
        mapT_ol = np.zeros(self.nbpix,dtype='int')
#        if self.indT[0] != 0: mapT_ol[self.indT[0]] = 0
        h.cartview(self.T/self.stdT,max=5,min=-5,rot=[lon,lat],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='T S/N')
        h.graticule(dpar=10,dmer=10,coord='C')
        py.savefig(self.outdir+'/png/'+self.dir_coadd+'/mapT_sn.png')

        mapQ_ol = np.zeros(self.nbpix,dtype='int')
#        if self.indQ_ol[0] != 0: mapQ_ol[self.indQ_ol[0]] = 1
        h.cartview(self.Q/self.stdQ,max=5,min=-5,rot=[lon,lat],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='Q S/N')
        h.graticule(dpar=10,dmer=10,coord='C')
        py.savefig(self.outdir+'/png/'+self.dir_coadd+'/mapQ_sn.png')

        mapU_ol = np.zeros(self.nbpix,dtype='int')
#        if self.indU_ol[0] != 0: mapU_ol[self.indU_ol[0]] = 1
        h.cartview(self.U/self.stdU,max=5,min=-5,rot=[lon,lat],lonra=[-width/2.,width/2.],latra=[-width/2.,width/2.],title='U S/N')
        h.graticule(dpar=10,dmer=10,coord='C')
        py.savefig(self.outdir+'/png/'+self.dir_coadd+'/mapU_sn.png')

        py.figure(10)
        par = plot_hist(self.T[self.indT[0]], 40, init_auto=True, xtitle='T', \
                            fname=self.outdir+'/png/'+self.dir_coadd+'/histT.png')
        py.figure(11)
        par = plot_hist(self.T[self.indT[0]], 40, init_auto=True, xtitle='T', \
                            fname=self.outdir+'/png/'+self.dir_coadd+'/histT_log.png',option_ylog=True)
        py.figure(12)
        par = plot_hist(self.Q[self.indQ[0]], 40, init_auto=True, xtitle='Q', \
                            fname=self.outdir+'/png/'+self.dir_coadd+'/histQ.png')
        py.figure(13)
        par = plot_hist(self.Q[self.indQ[0]], 40, init_auto=True, xtitle='Q', \
                            fname=self.outdir+'/png/'+self.dir_coadd+'/histQ_log.png',option_ylog=True)
        py.figure(14)
        par = plot_hist(self.U[self.indU[0]], 40, init_auto=True, xtitle='U', \
                            fname=self.outdir+'/png/'+self.dir_coadd+'/histU.png')
        py.figure(15)
        par = plot_hist(self.U[self.indU[0]], 40, init_auto=True, xtitle='U', \
                            fname=self.outdir+'/png/'+self.dir_coadd+'/histU_log.png',option_ylog=True)

    def clear_maps(self):
        self.T = 0
        self.Q = 0
        self.U = 0
        self.M = 0

        
def plot_hist(x,nbin,par=-1,init_auto=False,xtitle=-1, \
                  no_plot=False,normed=False,option_ylog=False,fname=''):
    """
    plot_hist.py: plot histogram and fit with a 2D gaussian
     inputs
         x: input
         nbin: number of bin
     options
         par: initial guess of parmaeters (amp,mu,sigma)
         fit: True/False
         init_auto: True/False (auto initial guess)
         xtitle: xtitle
     output:
         fit parameters
    """
    # the histogram of the data
    non, bins, patches = py.hist(x, nbin, histtype='step', normed=normed)#, normed=1, facecolor='green', alpha=0.75)
    
    bincenters = 0.5*(bins[1:]+bins[:-1])
    
    func_gauss = lambda p, xin: p[0]*np.exp(-(xin-p[1])**2/(2.*p[2]**2))
    chi_nosigma = lambda p, xin, d: ((func_gauss(p,xin)-d)**2).sum()
    
    print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    print '+++  Fit the histogram with Gaussian +++'
    if init_auto: par0 = [np.max(non),np.median(x),np.std(x)]
    if init_auto == False: par0 = par
    print 'initial guess:', par0
    x = np.arange(min(bincenters),max(bincenters),(max(bincenters)-min(bincenters))/500.)
    par = fmin(chi_nosigma, par0, args=(bincenters,non), maxiter=10000, maxfun=10000, xtol=0.01)
    py.plot(x,func_gauss(par,x),'r', linewidth=1)
    if not option_ylog:
        py.text(-400e-6,par[0]*0.9, 'A'+str(par[0]))
        py.text(-400e-6,par[0]*0.8, '$\mu$'+str(par[1]))
        py.text(-400e-6,par[0]*0.7, '$\sigma$'+str(par[2]))
    if option_ylog: 
        py.semilogy()
        py.text(-400e-6,par[0]*0.8, 'A'+str(par[0]))
        py.text(-400e-6,par[0]*0.4, '$\mu$'+str(par[1]))
        py.text(-400e-6,par[0]*0.05, '$\sigma$'+str(par[2]))
    py.xlabel(xtitle)
    if xtitle=='T': py.xlim([-500e-6,500e-6])
    if xtitle!='T': py.xlim([-30e-6,30e-6])
    py.savefig(fname)
    print 'fitted parameters:', par
    print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    return np.array(par)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sys_run_name, outdir, dir_coadd, select_ces = read_coaddDB(dir_out)
nb_coadd = len(dir_coadd)
ind = np.where( ( np.array(select_ces) == sqlite_command) )

if len(ind[0])==1:
    lon = float(xml_input["patch_ra"])
    lat = float(xml_input["patch_dec"])
    width = float(xml_input["patch_width"]) 
    cal_coaddmap_stats = cal_map_stats()
    cal_coaddmap_stats.outdir = outdir[ind[0]]
    cal_coaddmap_stats.dir_coadd = dir_coadd[ind[0]]
    cal_coaddmap_stats.read_coaddmaps()
    cal_coaddmap_stats.cal_stats()
    cal_coaddmap_stats.gen_plots(lon,lat,width)

if len(ind[0])==0:
    print "[cal_map_stats.py] no DB entry in coadd_map.db"
