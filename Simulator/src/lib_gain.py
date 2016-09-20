import numpy as np
import sqlite3 as sq
import sys
import os
import lib_LBSQLDB as lib_sq
import numpy.fft as fft

pi = np.pi

class gen_gain4mm():
    def _init_(self):
        self.sqlite_command = '.schema'
        self.filename = './tmp'
        self.gain_type = 'sim'
        self.gain_corr = 'sim'
        self.boloid = np.zeros(1)
        self.pix_list = pix_list
        self.num_det = 1

    def prep_relgain4mm_sim1(self):
        gain_in = np.arange(self.num_det,dtype='int')
        gain = np.ones(self.num_det)
        gain_mean = float(1.)
        gain_dict = {'boloid':gain_in, 'gain':gain, 'gain_mean':gain_mean}
        return gain_dict

#    def prep_relgain4mm(self):
#        Gain_cl = lib_sq.lib_GainDB()
#        Gain_cl.sqlite_command = self.sqlite_command
#        gain_in = Gain_cl.read_GainDB(self.filename,self.gain_type)
#        
#        num_gain = max(self.boloid)+1
#
#        gain_mean_tmp=0.; num_gain_eff = 0
#        gain = np.zeros(num_gain)
#
#        num_pix=len(self.pix_list[:][:])
#        pix_t=[]; pix_b=[]
#        for i in range(num_pix):
#            pix_t.append(self.pix_list[i][0])
#            pix_b.append(self.pix_list[i][1])
#
#        num_id = len(gain_in['boloid'])
#        for i in range(num_id):
#            ind_t = np.where( np.array(pix_t)==gain_in['boloid'][i])
#            if len(ind_t[0])!=0:
#                gain[pix_t[ind_t[0]]] = 0.5*(gain_in['g_begin_'+self.gain_type][i] + gain_in['g_end_'+self.gain_type][i])
#                num_gain_eff += 1
#                gain_mean_tmp += gain[pix_t[ind_t[0]]]
#            ind_b = np.where( np.array(pix_b)==gain_in['boloid'][i])
#            if len(ind_b[0])!=0:
#                gain[pix_b[ind_b[0]]] = 0.5*(gain_in['g_begin_'+self.gain_type][i] + gain_in['g_end_'+self.gain_type][i])
#                num_gain_eff += 1
#                gain_mean_tmp += gain[pix_b[ind_b[0]]]
#
#        gain_mean = float(gain_mean_tmp/float(num_gain_eff))
#        gain_dict = {'boloid':gain_in['boloid'], 'gain':gain, 'gain_mean':gain_mean}
#        return gain_dict

    def prep_relgain4mm(self,num_subscan,fsample):
        Gain_cl = lib_sq.lib_GainDB()
        Gain_cl.sqlite_command = self.sqlite_command
        gain_in = Gain_cl.read_GainDB(self.filename)
        
        num_det = len(gain_in['detid'])
        gain_out = np.zeros((num_det,num_subscan))

        seed_tmp = np.random.uniform(0,2e16,1)
        seed = int(seed_tmp[0])
        for i in range(num_det):
#            if 'flat_bias' in gain_in['gain_type'][i]:
#                gain_out[i][:] = 1.+ np.ones(num_subscan)*gain_in['params1'][i]
#            if 'flat_random' in gain_in['gain_type'][i]:
#                gain_out[i][:] = 1.+ np.random.normal(0.,1.,num_subscan)*gain_in['params1'][i]
#            if 'flat_1of' in gain_in['gain_type'][i]:
#                tmp = NoiseGen_auto(num_subscan,fsample, \
#                                        np.array(gain_in['params1'][i]), \
#                                        np.array(gain_in['params2'][i]), \
#                                        np.array(gain_in['params3'][i]), \
#                                        int(np.random.uniform(1e-10,1e10,1)))
#                gain_out[i][:] = 1. + tmp

            if 'bias' in self.gain_type: #gain_in['gain_type'][i]:
                gain_out[i][:] = 1.+ np.ones(num_subscan)*gain_in['params1'][i]
            if 'random_r' in self.gain_type: #gain_in['gain_type'][i]:
                gain_out[i][:] = 1.+ np.random.normal(0.,1.,num_subscan)*gain_in['params1'][i]
            if 'random_c' in self.gain_type: #gain_in['gain_type'][i]:
                np.random.seed(seed)
                gain_out[i][:] = 1.+ np.random.normal(0.,1.,num_subscan)*gain_in['params1'][i]
            if '1of_r' in self.gain_type: #gain_in['gain_type'][i]:
                tmp = NoiseGen_auto(num_subscan,fsample, \
                                        np.array(gain_in['params1'][i]), \
                                        np.array(gain_in['params2'][i]), \
                                        np.array(gain_in['params3'][i]), \
                                        int(np.random.uniform(1e-10,1e10,1)), \
                                        2)
                gain_out[i][:] = 1. + tmp
            if '1of_c' in self.gain_type: #gain_in['gain_type'][i]:
                np.random.seed(seed)
                tmp = NoiseGen_auto(num_subscan,fsample, \
                                        np.array(gain_in['params1'][i]), \
                                        np.array(gain_in['params2'][i]), \
                                        np.array(gain_in['params3'][i]), \
                                        int(np.random.uniform(1e-10,1e10,1)), \
                                        2)
                gain_out[i][:] = 1. + tmp

        return gain_out

def NoiseGen_auto(nbData,fsample,net,fknee,power,seed,model):
    '''
       generate the noise tod based on the 1of model
       nbData: number of data
       fsample: sampling rate in Hz
       net: noise level in unit of Vrtsec, uKrtsec
       knee: 1of knee in Hz
       power: 1of index
       seed: seed for noise generation
       
       the noise model 1 is sigma = net^2*fsample * [1+(fknee/f)^alpha]
       the noise model 2 is sigma = net^2*fsample * (fknee/f)^alpha

       return TOD with sigma of net * sqrt(fsample)
    '''
    nbf = int((nbData)/2.)+1
    freq = np.arange(nbf)/np.double(nbf-1)*(fsample)/2.
    if model == 1:
        psd = np.sqrt(net**2*fsample * (1.+(fknee/freq)**power))
    if model == 2:
        psd = np.sqrt(net**2*fsample * (fknee/freq)**power)
    # if nbData = 10000                                                                                                  
    #        0, 1     2,          4999,     5000,    5001,          , 9999                                               
    #        0, 1,    2,    ..., N/2-1,      N/2,   N/2+1,          , N-1                                                
    #        0, 1,    2,    ..., nbf-2,    nbf-1,     nbf,          , 2*(nbf-1)                                          
    #        0, (        idx         ),                                                                                  
    #        0, 1/NT, 2/NT, ..., (N/2-1)/NT, 1/(2T), -(N/2-1)/NT, ..., -1/NT                                             
    idx = range(1,nbf-1) #                                                                                  
    real_psd = np.concatenate(( psd[idx], np.array([psd[nbf-1]]), psd[idx[::-1]] ))
    real_psd = np.hstack([0,real_psd])
    
    np.random.seed(seed)
    f_rand = np.random.uniform(low=0.0, high=2.*pi, size=nbData)
    real_psdout = real_psd*np.cos(f_rand)
    imag_psdout = real_psd*np.sin(f_rand)
    psd_complex = complex_arr(real_psdout,imag_psdout)
    top_noise = fft.ifft(psd_complex)
    return np.array(top_noise.real)*np.sqrt(2.*nbData)


def complex_arr(arr1,arr2):
    nb1 = np.size(arr1)
    nb2 = np.size(arr2)

    if nb1 == nb2:
        out_arr = np.zeros(nb1,complex)
        out_arr.real=arr1
        out_arr.imag=arr2

    if nb1!=nb2 and nb1==1:
        out_arr = np.zeros(nb2,complex)
        i = np.arange(1,nb2+1)/np.arange(1,nb2+1)
        out_arr.real=i*arr1
        out_arr.imag=arr2

    if nb1!=nb2 and nb2==1:
        out_arr = np.zeros(nb1,complex)
        i = np.arange(1,nb1+1)/np.arange(1,nb1+1)
        out_arr.real=arr1
        out_arr.imag=i*arr2

    if nb1==1 and nb2==1:
        out_arr = np.eros(nb1,complex)
        out_arr.real=arr1
        out_arr.imag=arr2

    return np.array(out_arr)
