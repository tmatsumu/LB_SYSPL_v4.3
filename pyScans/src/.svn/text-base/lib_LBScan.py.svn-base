import numpy as np
import pylab as py
import healpy as h
import ephem as ephem
import time
import sys
import matsumulib as mylib
import lib_LBScan_c as liblbscanc
import sqlite3 as sq
import os 

# define basic parameters
pi = np.pi
radeg = (180./pi)

def write_ptgdb(db_name,juliantime,filename):
    if os.path.exists(db_name):
        print 'DB already exists'
        conn = sq.connect(db_name)
        c = conn.cursor()
        list_entries = (float(juliantime), filename)
        c.execute('insert into LBSimPtg(juliantime,filename) values (?,?)', list_entries)
        conn.commit()
        c.close()
    else:
        conn = sq.connect(db_name)
        c = conn.cursor()
        c.execute('create table LBSimPtg (id integer primary key autoincrement, juliantime real, filename text)')
        list_entries = (float(juliantime), filename)
        c.execute('insert into LBSimPtg(juliantime,filename) values (?,?)', list_entries)
        conn.commit()
        c.close()

def write_ptgidxdb(db_name,juliantime,samplerate,freq_boresight):
    num = int(24*3600*freq_boresight)
#    print (24*3600*freq_boresight), float(num), (24*3600*freq_boresight)-float(num)
#    if (24*3600*freq_boresight)-float(num)>0.5: num+=1 
#    print num, 24*3600*freq_boresight
#    sys.exit()
    sample_period = int(samplerate/freq_boresight)
    for i in range(0,num):
        idx_s = i*sample_period
        idx_e = (i+1)*sample_period-1
        if i==num-1:
            idx_e = int(samplerate*24.*3600.)-1
        julian_s = juliantime + idx_s/samplerate/24./3600.
        julian_e = juliantime + idx_e/samplerate/24./3600.
        if os.path.exists(db_name):
            print ''
            print 'DB already exists!'
            print 'Please delete the file or check the file name.'
            print 'Currently overwriting ', db_name
            print ''
            conn = sq.connect(db_name)
            c = conn.cursor()
            list_entries = (float(julian_s), float(julian_e), idx_s, idx_e)
            c.execute('insert into LBSimPtgIdx(julian_s,julian_e,idx_s,idx_e) values (?,?,?,?)', list_entries)
            conn.commit()
            c.close()
        else:
            conn = sq.connect(db_name)
            c = conn.cursor()
            c.execute('create table LBSimPtgIdx (id integer primary key autoincrement, julian_s real, julian_e real, idx_s integer, idx_e integer)')
            list_entries = (float(julian_s), float(julian_e), idx_s, idx_e)
            c.execute('insert into LBSimPtgIdx(julian_s,julian_e,idx_s,idx_e) values (?,?,?,?)', list_entries)
            conn.commit()
            c.close()
    
def rot_xy(phi):
    out = np.matrix([[np.cos(phi),np.sin(phi),0.], [-np.sin(phi),np.cos(phi),0.], [0.,0.,1.]])
    return out

def rot_xz(theta):
    out = np.matrix([[np.cos(theta),0.,np.sin(theta)], [0., 1.,0.], [-np.sin(theta),0.,np.cos(theta)]])
    return out

def sun_position_quick(DJD):
    freq = 1./((365.25)*24.*60.*60.)
    phi = mylib.wraparound_2pi(2.*pi*freq*DJD*(8.64*1.e4))
    return phi, np.zeros(len(phi))

def sun_position(DJD):
    nb = len(DJD)
    ra_arr = []; dec_arr = []
    for i in range(nb):
        time = ephem.Date(DJD[i])
        sun = ephem.Sun(time)
        ra_arr.append(float(repr(sun.ra)))
        dec_arr.append(float(repr(sun.dec)))
    return np.array(ra_arr), np.array(dec_arr)

def moon_position(DJD):
    nb = len(DJD)
    ra_arr = []; dec_arr = []
    for i in range(nb):
        time = ephem.Date(DJD[i])
        moon = ephem.Moon(time)
        ra_arr.append(float(repr(moon.ra)))
        dec_arr.append(float(repr(moon.dec)))
    return np.array(ra_arr), np.array(dec_arr)

def convert_Gregorian2Julian(year, month, day, hour, minute, sec):
    a = int((14-month)/12)
    y = year + 4800 - a
    m = month + 12*a - 3
    JDN = day+int((153*m+2)/5)+365*y+int(y/4)-int(y/100)+int(y/400)-32045
    JD = JDN + (hour-12)/24. + minute/1440. + sec/86400.
    return np.array(JD)

def convert_Julian2Dublin(JD):
    DJD = JD - 2415020.
    return np.array(DJD)

def ang2pos_3D(theta,phi):
    x = np.sin(theta)*np.cos(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(theta)
    return x, y, z

def ang2_scandirec_3D(x,y,z):
    xx = np.diff(x)
    yy = np.diff(y)
    zz = np.diff(z)
    return np.array([xx, yy, zz])

def ang2_deriv_theta_3D(theta,phi):
    nb = len(theta)
    x =  np.cos(theta)*np.cos(phi)
    y =  np.cos(theta)*np.sin(phi)
    z = -np.sin(theta)
    return np.array([x[0:nb-1],y[0:nb-1],z[0:nb-1]])

def ang2_deriv_phi_3D(theta,phi):
    nb = len(theta)
    x = -np.sin(phi)
    y =  np.cos(phi)
    z =  np.zeros(nb)
    return np.array([x[0:nb-1],y[0:nb-1],z[0:nb-1]])

def Cal_scanangle(point,theta,phi):
    alpha = np.arccos(point*theta/np.sqrt(point)/np.sqrt(theta))
    beta  = np.arccos(point*phi/np.sqrt(point)/np.sqrt(phi))
    nb = len(alpha)
    for i in range(nb):
        if ((alpha[i] >= 0./180.*pi) and (alpha[i] < 90./180.*pi) 
            and (beta[i] >= 0./180.*pi) and (beta[i] < 90./180.*pi)):
            alpha[i] = alpha[i]
            beta[i] = pi/2. - beta[i]

        if ((alpha[i] >= 90./180.*pi) and (alpha[i] < 180./180.*pi) 
            and (beta[i] >= 0./180.*pi) and (beta[i] < 90./180.*pi)):
            alpha[i] = alpha[i]
            beta[i] = pi/2. + beta[i]
        
        if ((alpha[i] >= 0./180.*pi) and (alpha[i] < 90./180.*pi) 
            and (beta[i] >= 90./180.*pi) and (beta[i] < 180./180.*pi)):
            alpha[i] = 2.*pi-alpha[i]
            beta[i] = pi/2. + beta[i]
     
        if ((alpha[i] >= 90./180.*pi) and (alpha[i] < 180./180.*pi) 
            and (beta[i] >= 90./180.*pi) and (beta[i] < 180./180.*pi)):
            alpha[i] = 2.*pi-alpha[i]
            beta[i] = 3.*pi/2. - beta[i]

        cos_clfom1[ipix[i]] += np.cos(alpha[i])
        sin_clfom1[ipix[i]] += np.sin(alpha[i])
        cos_clfom2[ipix[i]] += np.cos(2.*alpha[i])
        sin_clfom2[ipix[i]] += np.sin(2.*alpha[i])
        cos_clfom4[ipix[i]] += np.cos(4.*alpha[i])
        sin_clfom4[ipix[i]] += np.sin(4.*alpha[i])

    return cos_clfom1, sin_clfom1, cos_clfom2, sin_clfom2, cos_clfom4, sin_clfom4, alpha
            

#; ----------------------------------------------------------------------------
def gen_scan( theta_antisun, theta_boresight, freq_antisun, freq_boresight, \
                  total_time, today_julian, sample_rate, dir_out, filename, title, nside, runtime_i, option_gen_ptg):
    
    # define time related variables
    n_time = long(total_time*sample_rate)
    print '[lib_LBscan.py] The number of samples in a day', n_time
    sec2date = (1.e-4/8.64)
    time_julian = np.arange(n_time)/(n_time-1.)*total_time*sec2date + today_julian # time in julian [day]
    
    # cal sun position
    DJD = mylib.convert_Julian2Dublin(time_julian)
    ra_sp, dec_sp = sun_position(DJD)

    # convert ra/dec to ecliptic coordinate
    lon_spe, lat_spe = mylib.EULER( ra_sp*radeg, dec_sp*radeg, 3, 'J2000')
    ra_sp = 0
    dec_sp = 0
   
    time_i = time_julian /sec2date # time in [sec]
    
    # convert the sun vector to anti-sun vector
    lat_aspe = -lat_spe/radeg
    lon_aspe = lon_spe/radeg + pi
    
    lat_spe = 0
    lon_spe = 0
    
#    idxwrap = np.where( lon_aspe > 2.*pi )
#    if len(idxwrap[0]) == 1:
#        if idxwrap[0] == -1: lon_aspe = lon_aspe
#        if idxwrap[0] != -1: lon_aspe[idxwrap[0]] = lon_aspe[idxwrap[0]] - 2.*pi
#    else: lon_aspe[idxwrap[0]] = lon_aspe[idxwrap[0]] - 2.*pi
    
    # from ecliptic lat, lon convention to theta, phi convention
    theta_asp = pi/2.-lat_aspe
    phi_asp = lon_aspe

    x_asp,y_asp,z_asp = ang2pos_3D(theta_asp,phi_asp)

    phi_antisun = mylib.wraparound_2pi(2.*pi*freq_antisun*time_i)
    phi_boresight = mylib.wraparound_2pi(2.*pi*freq_boresight*time_i)
    
    theta_out = np.zeros(n_time)
    phi_out = np.zeros(n_time)

    print '[lib_LBscan.py,gen_scan] current run time: before matrix', (time.time()-runtime_i)/60., 'min'
    
    for i in range(0,n_time):
        if y_asp[i] >= 0: rel_phi = np.arccos(x_asp[i])
        if y_asp[i] < 0: rel_phi = -np.arccos(x_asp[i]) + 2.*pi

        p_asp = [[x_asp[i]],[y_asp[i]],[z_asp[i]]]        

        p_out = rot_xy(-rel_phi) \
            * rot_xz(-pi/2.) \
            * rot_xy(-phi_antisun[i]) \
            * rot_xz(-theta_antisun) \
            * rot_xy(-phi_boresight[i]) \
            * rot_xz(-theta_boresight) \
            * rot_xz(pi/2.) \
            * rot_xy(rel_phi) \
            * p_asp              
            #            * rot_xz(theta_antisun) \
            #            * rot_xy(phi_antisun[i]) \
            #            * rot_xy(-phi_antisun[i]) \
            #            * rot_xz(-theta_antisun) \

        theta_out[i] = np.arctan2(np.sqrt(p_out[0]**2+p_out[1]**2),p_out[2]) 
        phi_out[i] = np.arctan2(p_out[1],p_out[0])+pi

    print 'current run time: after matrix', (time.time()-runtime_i)/60., 'min'
    print 'after the first for loop'

#    py.subplot(311)
#    py.plot(phi_out/pi*180.,theta_out/pi*180.,'.')
#    py.subplot(312)
#    py.plot(theta_out/pi*180.)
#    py.subplot(313)
#    py.plot(phi_out/pi*180.)
#    py.show()
#    sys.exit()

    p_asp = 0
    x_asp = 0
    y_asp = 0
    z_asp = 0
    lon_aspe = 0
    lat_aspe = 0
  
    ipix = h.ang2pix( nside, theta_out, phi_out)
    
    nbPix_scan = len(ipix)-1
    nbPix = h.nside2npix(nside)
        
    dpvec = np.zeros((3,nbPix_scan))
    dtvec = np.zeros((3,nbPix_scan))
    dphivec = np.zeros((3,nbPix_scan))
    
    xp,yp,zp = ang2pos_3D(theta_out,phi_out)
    dpvec = ang2_scandirec_3D(xp,yp,zp)
    dtvec = ang2_deriv_theta_3D(theta_out,phi_out)
    dphivec = ang2_deriv_phi_3D(theta_out,phi_out)
    
    alpha = np.zeros(nbPix_scan)
    alpha = np.arccos( (dpvec[0,:]*dtvec[0,:]+dpvec[1,:]*dtvec[1,:]+dpvec[2,:]*dtvec[2,:])
                       /np.sqrt(dpvec[0,:]*dpvec[0,:]+dpvec[1,:]*dpvec[1,:]+dpvec[2,:]*dpvec[2,:])
                       /np.sqrt(dtvec[0,:]*dtvec[0,:]+dtvec[1,:]*dtvec[1,:]+dtvec[2,:]*dtvec[2,:]) )
    
    beta = np.zeros(nbPix_scan)
    beta = np.arccos( (dpvec[0,:]*dphivec[0,:]+dpvec[1,:]*dphivec[1,:]+dpvec[2,:]*dphivec[2,:])
                      /np.sqrt(dpvec[0,:]*dpvec[0,:]+dpvec[1,:]*dpvec[1,:]+dpvec[2,:]*dpvec[2,:])
                      /np.sqrt(dphivec[0,:]*dphivec[0,:]+dphivec[1,:]*dphivec[1,:]+dphivec[2,:]*dphivec[2,:]) )
    
    print 'current run time: before summing up ', (time.time()-runtime_i)/60., 'min'
    
    Nhits = np.zeros(nbPix)
    clfom = np.zeros(nbPix)
    cos_r1 = np.zeros(nbPix)
    sin_r1 = np.zeros(nbPix)
    cos_r2 = np.zeros(nbPix)
    sin_r2 = np.zeros(nbPix)
    cos_r4 = np.zeros(nbPix)
    sin_r4 = np.zeros(nbPix)

    for i in range(0,nbPix_scan):
        Nhits[ipix[i]] += 1.
#        if ((alpha[i] >= 0./180.*pi) and (alpha[i] < 90./180.*pi) 
#            and (beta[i] >= 0./180.*pi) and (beta[i] < 90./180.*pi)):
#            alpha[i] = alpha[i]
#
#        if ((alpha[i] >= 90./180.*pi) and (alpha[i] < 180./180.*pi) 
#            and (beta[i] >= 0./180.*pi) and (beta[i] < 90./180.*pi)):
#            alpha[i] = alpha[i]
#
#        if ((alpha[i] >= 0./180.*pi) and (alpha[i] < 90./180.*pi) 
#            and (beta[i] >= 90./180.*pi) and (beta[i] < 180./180.*pi)):
#            alpha[i] = - alpha[i]
#     
#        if ((alpha[i] >= 90./180.*pi) and (alpha[i] < 180./180.*pi) 
#            and (beta[i] >= 90./180.*pi) and (beta[i] < 180./180.*pi)):
#            alpha[i] = - alpha[i]

#        print ''
#        print alpha[i]/pi*180., beta[i]/pi*180., alpha[i]/pi*180.+beta[i]/pi*180., alpha[i]/pi*180.-beta[i]/pi*180.
#        if (alpha[i] >= 180./180.*pi): alpha[i] = pi - alpha[i]
        if ((beta[i] >= 90./180.*pi) and (beta[i] < 180./180.*pi)):
            alpha[i] = - alpha[i]
        
#        print alpha[i]/pi*180., beta[i]/pi*180., alpha[i]/pi*180.+beta[i]/pi*180., alpha[i]/pi*180.-beta[i]/pi*180.
#        if i == 10000: return
        cos_r1[ipix[i]] += np.cos(alpha[i])
        sin_r1[ipix[i]] += np.sin(alpha[i])
        cos_r2[ipix[i]] += np.cos(2.*alpha[i])
        sin_r2[ipix[i]] += np.sin(2.*alpha[i])
        cos_r4[ipix[i]] += np.cos(4.*alpha[i])
        sin_r4[ipix[i]] += np.sin(4.*alpha[i])

    ############################################################################3
    if option_gen_ptg == True:
        fname=dir_out+'/ptg/'+title+'_'+str(int(today_julian))+'_latlon_eclip.txt'
        tmpstr = time_julian[0]-2400000.5
        time_arr = np.arange(len(alpha))
        print len(theta_out), len(phi_out), len(alpha), len(time_arr)
        mylib.write_txt5(fname,time_arr/float(sample_rate)+tmpstr, pi/2.-theta_out, phi_out, alpha, beta)
    ############################################################################3
     
    runtime_m = time.time()
    print 'current run time: after summing up', (runtime_m-runtime_i)/60., 'min'
    h.write_map( 'dataout/'+filename+'/fits/nhits_tmp.fits', Nhits)
    h.write_map( 'dataout/'+filename+'/fits/cos_r1_tmp.fits', cos_r1)
    h.write_map( 'dataout/'+filename+'/fits/sin_r1_tmp.fits', sin_r1)
    h.write_map( 'dataout/'+filename+'/fits/cos_r2_tmp.fits', cos_r2)
    h.write_map( 'dataout/'+filename+'/fits/sin_r2_tmp.fits', sin_r2)
    h.write_map( 'dataout/'+filename+'/fits/cos_r4_tmp.fits', cos_r4)
    h.write_map( 'dataout/'+filename+'/fits/sin_r4_tmp.fits', sin_r4)
    
    runtime_f = time.time()
    print 'final run time: ', (runtime_f-runtime_i)/60., 'min'
    print '[lib_LBscan.py/gen_scan] End of gen_scan()' 
    print ''


#; ----------------------------------------------------------------------------
def gen_scan_c( theta_antisun, theta_boresight, freq_antisun, freq_boresight, \
                    total_time, today_julian, sample_rate, dir_out, filename, \
                    title, nside, runtime_i, option_gen_ptg):
    
    # define time related variables
    n_time = long(total_time*sample_rate)
    print '[lib_LBscan.py] The number of samples in a day', n_time
    sec2date = (1.e-4/8.64)
    time_julian = np.arange(n_time)/(n_time-1.)*total_time*sec2date + today_julian # time in julian [day]
    
    # cal sun position
    DJD = mylib.convert_Julian2Dublin(time_julian)
    ra_sp, dec_sp = sun_position(DJD)

    # convert ra/dec to ecliptic coordinate
    lon_spe, lat_spe = mylib.EULER( ra_sp*radeg, dec_sp*radeg, 3, 'J2000')
    ra_sp = 0
    dec_sp = 0
   
    time_i = time_julian /sec2date # time in [sec]
    
    # convert the sun vector to anti-sun vector
    lat_aspe = -lat_spe/radeg
    lon_aspe = lon_spe/radeg + pi
    
    lat_spe = 0
    lon_spe = 0
    
#    idxwrap = np.where( lon_aspe > 2.*pi )
#    if len(idxwrap[0]) == 1:
#        if idxwrap[0] == -1: lon_aspe = lon_aspe
#        if idxwrap[0] != -1: lon_aspe[idxwrap[0]] = lon_aspe[idxwrap[0]] - 2.*pi
#    else: lon_aspe[idxwrap[0]] = lon_aspe[idxwrap[0]] - 2.*pi
    
    # from ecliptic lat, lon convention to theta, phi convention
    theta_asp = pi/2.-lat_aspe
    phi_asp = lon_aspe

    x_asp,y_asp,z_asp = ang2pos_3D(theta_asp,phi_asp)
#    for i in range(0,len(x_asp)):
#        print '>', theta_asp[i]/pi*180., phi_asp[i]/pi*180., x_asp[i], y_asp[i], z_asp[i]

    phi_antisun = mylib.wraparound_2pi(2.*pi*freq_antisun*time_i)
    phi_boresight = mylib.wraparound_2pi(2.*pi*freq_boresight*time_i)
    
#    for i in range(0,len(x_asp)):
#        print '>', theta_asp[i]/pi*180., phi_asp[i]/pi*180., x_asp[i], y_asp[i], z_asp[i], phi_antisun[i]/pi*180., phi_boresight[i]/pi*180.

 #   py.subplot(211)
 #   py.plot(phi_antisun/pi*180.)
 #   py.subplot(212)
 #   py.plot(phi_boresight/pi*180.)
 #   py.show()

#    theta_out = np.zeros(n_time)
#    phi_out = np.zeros(n_time)

    print '[lib_LBscan.py,gen_scan] current run time: before matrix', (time.time()-runtime_i)/60., 'min'
    
#    print theta_asp, phi_asp
#    print theta_antisun/pi*180., phi_antisun
#    print theta_boresight/pi*180., phi_boresight

#    py.subplot(211)
#    py.plot(phi_antisun)
#    py.subplot(212)
#    py.plot(phi_boresight)
#    py.show()
#    sys.exit()

    p_out = liblbscanc.LB_rotmatrix_multi(theta_asp,
                                          phi_asp,
                                          theta_antisun,
                                          phi_antisun,
                                          theta_boresight,
                                          phi_boresight)
    #    for i in range(0,n_time):
#        if y_asp[i] >= 0: rel_phi = np.arccos(x_asp[i])
#        if y_asp[i] < 0: rel_phi = -np.arccos(x_asp[i]) + 2.*pi
#        
#        p_asp = [[x_asp[i]],[y_asp[i]],[z_asp[i]]]
#        p_out = rot_xy(-rel_phi) * rot_xz(-pi/2.) \
#            * rot_xy(-phi_antisun[i]) * rot_xz(-theta_antisun) \
#            * rot_xy(-phi_boresight[i]) * rot_xz(-theta_boresight) \
#            * rot_xz(theta_antisun) * rot_xy(phi_antisun[i]) \
#            * rot_xy(-phi_antisun[i]) * rot_xz(-theta_antisun) \
#            * rot_xz(pi/2.) * rot_xy(rel_phi) \
#            * p_asp              

    theta_out = np.arctan2(np.sqrt(p_out[0,:]**2+p_out[1,:]**2),p_out[2,:]) 
    phi_out = np.arctan2(p_out[1,:],p_out[0,:])+pi

#    py.subplot(311)
#    py.plot(phi_out/pi*180.,theta_out/pi*180.,'.')
#    py.subplot(312)
#    py.plot(theta_out/pi*180.)
#    py.subplot(313)
#    py.plot(phi_out/pi*180.)
#    py.show()
#    sys.exit()

    print 'current run time: after matrix', (time.time()-runtime_i)/60., 'min'
    print 'after the first for loop'
    
#    p_asp = 0
#    x_asp = 0
#    y_asp = 0
#    z_asp = 0
    lon_aspe = 0
    lat_aspe = 0
  
    ipix = h.ang2pix( nside, theta_out, phi_out)
    
    nbPix_scan = len(ipix)-1
    nbPix = h.nside2npix(nside)
        
    dpvec = np.zeros((3,nbPix_scan))
    dtvec = np.zeros((3,nbPix_scan))
    dphivec = np.zeros((3,nbPix_scan))
    
    xp,yp,zp = ang2pos_3D(theta_out,phi_out)
    dpvec = ang2_scandirec_3D(xp,yp,zp)
    dtvec = ang2_deriv_theta_3D(theta_out,phi_out)
    dphivec = ang2_deriv_phi_3D(theta_out,phi_out)
    
    alpha = np.zeros(nbPix_scan)
    alpha = np.arccos( (dpvec[0,:]*dtvec[0,:]+dpvec[1,:]*dtvec[1,:]+dpvec[2,:]*dtvec[2,:])
                       /np.sqrt(dpvec[0,:]*dpvec[0,:]+dpvec[1,:]*dpvec[1,:]+dpvec[2,:]*dpvec[2,:])
                       /np.sqrt(dtvec[0,:]*dtvec[0,:]+dtvec[1,:]*dtvec[1,:]+dtvec[2,:]*dtvec[2,:]) )
    
    beta = np.zeros(nbPix_scan)
    beta = np.arccos( (dpvec[0,:]*dphivec[0,:]+dpvec[1,:]*dphivec[1,:]+dpvec[2,:]*dphivec[2,:])
                      /np.sqrt(dpvec[0,:]*dpvec[0,:]+dpvec[1,:]*dpvec[1,:]+dpvec[2,:]*dpvec[2,:])
                      /np.sqrt(dphivec[0,:]*dphivec[0,:]+dphivec[1,:]*dphivec[1,:]+dphivec[2,:]*dphivec[2,:]) )
    
    print 'current run time: before summing up ', (time.time()-runtime_i)/60., 'min'
    
    Nhits = np.zeros(nbPix)
    clfom = np.zeros(nbPix)
    cos_r1 = np.zeros(nbPix)
    sin_r1 = np.zeros(nbPix)
    cos_r2 = np.zeros(nbPix)
    sin_r2 = np.zeros(nbPix)
    cos_r4 = np.zeros(nbPix)
    sin_r4 = np.zeros(nbPix)

    for i in range(0,nbPix_scan):
        Nhits[ipix[i]] += 1.
#        if ((alpha[i] >= 0./180.*pi) and (alpha[i] < 90./180.*pi) 
#            and (beta[i] >= 0./180.*pi) and (beta[i] < 90./180.*pi)):
#            alpha[i] = alpha[i]
#
#        if ((alpha[i] >= 90./180.*pi) and (alpha[i] < 180./180.*pi) 
#            and (beta[i] >= 0./180.*pi) and (beta[i] < 90./180.*pi)):
#            alpha[i] = alpha[i]
#
#        if ((alpha[i] >= 0./180.*pi) and (alpha[i] < 90./180.*pi) 
#            and (beta[i] >= 90./180.*pi) and (beta[i] < 180./180.*pi)):
#            alpha[i] = - alpha[i]
#     
#        if ((alpha[i] >= 90./180.*pi) and (alpha[i] < 180./180.*pi) 
#            and (beta[i] >= 90./180.*pi) and (beta[i] < 180./180.*pi)):
#            alpha[i] = - alpha[i]

#        print ''
#        print alpha[i]/pi*180., beta[i]/pi*180., alpha[i]/pi*180.+beta[i]/pi*180., alpha[i]/pi*180.-beta[i]/pi*180.
#        if (alpha[i] >= 180./180.*pi): alpha[i] = pi - alpha[i]
        if ((beta[i] >= 90./180.*pi) and (beta[i] < 180./180.*pi)):
            alpha[i] = - alpha[i]
        
#        print alpha[i]/pi*180., beta[i]/pi*180., alpha[i]/pi*180.+beta[i]/pi*180., alpha[i]/pi*180.-beta[i]/pi*180.
#        if i == 10000: return
        cos_r1[ipix[i]] += np.cos(alpha[i])
        sin_r1[ipix[i]] += np.sin(alpha[i])
        cos_r2[ipix[i]] += np.cos(2.*alpha[i])
        sin_r2[ipix[i]] += np.sin(2.*alpha[i])
        cos_r4[ipix[i]] += np.cos(4.*alpha[i])
        sin_r4[ipix[i]] += np.sin(4.*alpha[i])

    ############################################################################3
    if option_gen_ptg == True:
        fname=dir_out+'/ptg/'+title+'_'+str(int(today_julian))+'_latlon_eclip.txt'
        tmpstr = time_julian[0]-2400000.5
        time_arr = np.arange(len(alpha))
        print len(theta_out), len(phi_out), len(alpha), len(time_arr)
        mylib.write_txt5(fname,time_arr/float(sample_rate)+tmpstr, pi/2.-theta_out, phi_out, alpha, beta)
    ############################################################################3
     
    runtime_m = time.time()
    print 'current run time: after summing up', (runtime_m-runtime_i)/60., 'min'
    h.write_map( 'dataout/'+filename+'/fits/nhits_tmp.fits', Nhits)
    h.write_map( 'dataout/'+filename+'/fits/cos_r1_tmp.fits', cos_r1)
    h.write_map( 'dataout/'+filename+'/fits/sin_r1_tmp.fits', sin_r1)
    h.write_map( 'dataout/'+filename+'/fits/cos_r2_tmp.fits', cos_r2)
    h.write_map( 'dataout/'+filename+'/fits/sin_r2_tmp.fits', sin_r2)
    h.write_map( 'dataout/'+filename+'/fits/cos_r4_tmp.fits', cos_r4)
    h.write_map( 'dataout/'+filename+'/fits/sin_r4_tmp.fits', sin_r4)
    
    runtime_f = time.time()
    print 'final run time: ', (runtime_f-runtime_i)/60., 'min'
    print '[lib_LBscan.py/gen_scan_c] End of gen_scan_c()' 
    print ''



#; ----------------------------------------------------------------------------
def gen_scan_c_mod( theta_antisun, theta_boresight, freq_antisun, freq_boresight, \
                        total_time, today_julian, sample_rate, dir_out, filename, \
                        title, nside, runtime_i, option_gen_ptg):
    
    print '[lib_LBscan.py,gen_scan_c_mod] initial', (time.time()-runtime_i)/60., 'min'    
    
    # define time related variables
    n_time = long(total_time*sample_rate)
    print '[lib_LBscan.py,gen_scan_c_mod] The number of samples in a day', n_time
    sec2date = (1.e-4/8.64)
    time_julian = np.arange(n_time)/(n_time-1.)*total_time*sec2date + today_julian # time in julian [day]
    
    # cal sun position
    print '[lib_LBscan.py,gen_scan_c_mod] before sun_position_quick', (time.time()-runtime_i)/60., 'min'    
    DJD = mylib.convert_Julian2Dublin(time_julian)
    phi_asp, theta_asp = sun_position_quick(DJD)
    DJD = 0;
    
    time_i = time_julian /sec2date # time in [sec]
    
    # from ecliptic lat, lon convention to theta, phi convention
    theta_asp = pi/2.-theta_asp

    phi_antisun = mylib.wraparound_2pi(2.*pi*freq_antisun*time_i)
    phi_boresight = mylib.wraparound_2pi(2.*pi*freq_boresight*time_i)
    
    print '[lib_LBscan.py,gen_scan_c_mod] before LB_rotmatrix_multi', (time.time()-runtime_i)/60., 'min'    
    p_out = liblbscanc.LB_rotmatrix_multi(theta_asp, phi_asp,
                                          theta_antisun, phi_antisun,
                                          theta_boresight, phi_boresight)

    theta_out = np.arctan2(np.sqrt(p_out[0,:]**2+p_out[1,:]**2),p_out[2,:]) 
    phi_out = np.arctan2(p_out[1,:],p_out[0,:])+pi
    theta_asp=0; phi_asp=0;
    p_out=0; time_i=0;
    phi_antisun=0; phi_boresight=0;

    ipix = h.ang2pix( nside, theta_out, phi_out)
    
    nbPix_scan = len(ipix)-1
    nbPix = h.nside2npix(nside)
        
    dpvec = np.zeros((3,nbPix_scan))
    dtvec = np.zeros((3,nbPix_scan))
    dphivec = np.zeros((3,nbPix_scan))
    
    xp,yp,zp = ang2pos_3D(theta_out,phi_out)
    dpvec = ang2_scandirec_3D(xp,yp,zp)
    dtvec = ang2_deriv_theta_3D(theta_out,phi_out)
    dphivec = ang2_deriv_phi_3D(theta_out,phi_out)
    
    alpha = np.zeros(nbPix_scan)
    alpha = np.arccos( (dpvec[0,:]*dtvec[0,:]+dpvec[1,:]*dtvec[1,:]+dpvec[2,:]*dtvec[2,:])
                       /np.sqrt(dpvec[0,:]*dpvec[0,:]+dpvec[1,:]*dpvec[1,:]+dpvec[2,:]*dpvec[2,:])
                       /np.sqrt(dtvec[0,:]*dtvec[0,:]+dtvec[1,:]*dtvec[1,:]+dtvec[2,:]*dtvec[2,:]) )
    
    beta = np.zeros(nbPix_scan)
    beta = np.arccos( (dpvec[0,:]*dphivec[0,:]+dpvec[1,:]*dphivec[1,:]+dpvec[2,:]*dphivec[2,:])
                      /np.sqrt(dpvec[0,:]*dpvec[0,:]+dpvec[1,:]*dpvec[1,:]+dpvec[2,:]*dpvec[2,:])
                      /np.sqrt(dphivec[0,:]*dphivec[0,:]+dphivec[1,:]*dphivec[1,:]+dphivec[2,:]*dphivec[2,:]) )
    
    print '[lib_LBscan.py,gen_scan_c_mod] before summing up', (time.time()-runtime_i)/60., 'min'    
    
    ind = np.where((beta >= 90./180.*pi) & (beta < 180./180.*pi))
    alpha[ind[0]] = - alpha[ind[0]]

    ############################################################################3
    if option_gen_ptg == True:
        fname=dir_out+'/ptg/LBPTG_latlon_eclip_'+str(int(today_julian))
        tmpstr = time_julian[0]-2400000.5
        time_arr = np.arange(len(alpha)+1)
        print len(theta_out), len(phi_out), len(alpha), len(time_arr)
#        outdict = {'julian_t0':time_julian[0], 'samplerate': sample_rate, 'lat': pi/2.-theta_out, 'lon': phi_out, 'pa': alpha, 'help': 'This dictionary in npy contains with the keyword of \n julian_t0 (julian time of the first sample), \n samplerate (samplerate in [Hz]), \n lon (ecliptic coordinate in [rad]) \n lat (ecliptic coordinate in [rad]), \n pa (alpha angle [rad] wrt theta hat angle) \n This pointing file is generated by the gen_scan_c_mod in lib_LBScan.py written by T.Matsumura.'}
        np.savez(fname, \
                     julian_t0=time_julian[0], \
                     sample_rate= sample_rate, \
                     lat=pi/2.-theta_out, \
                     lon=phi_out, \
                     pa=alpha, \
                     help='This dictionary in npy contains with the keyword of \n julian_t0 (julian time of the first sample), \n samplerate (sample_rate in [Hz]), \n lon (ecliptic coordinate in [rad]) \n lat (ecliptic coordinate in [rad]), \n pa (alpha angle [rad] wrt theta hat angle) \n This pointing file is generated by the gen_scan_c_mod in lib_LBScan.py written by T.Matsumura.')
        db_name = dir_out+'/ptg/LBPTG_latlon_eclip'
        fname_tmp='LBPTG_latlon_eclip_'+str(int(today_julian))
        write_ptgdb(db_name+'.db',time_julian[0],dir_out+'/ptg/'+fname_tmp)
        write_ptgidxdb(db_name+'_'+str(int(today_julian))+'.db',today_julian,sample_rate,freq_boresight)
        del(time_arr)
    del(time_julian)
    ############################################################################3

    beta=0; dphivec=0; dpvec=0; dtvec=0; theta_out=0; phi_out=0
    mapout = liblbscanc.Maps_summingup(nbPix, np.float_(ipix), alpha)
    print '[lib_LBscan.py,gen_scan_c_mod] before writing text file ['+str(option_gen_ptg)+']', (time.time()-runtime_i)/60., 'min'   
    ipix=0; alpha=0
    runtime_f = time.time()
    print '[lib_LBscan.py,gen_scan_c_mod] END of gen_scan_c_mod ', (time.time()-runtime_i)/60., 'min'    
    print ''
    return mapout



#; ----------------------------------------------------------------------------
def gen_scan_c_mod_moon( theta_antisun, theta_boresight, freq_antisun, freq_boresight, \
                         total_time, today_julian, sample_rate, dir_out, filename, \
                         title, nside, runtime_i, option_gen_ptg):
    
    print '[lib_LBscan.py,gen_scan_c_mod] initial', (time.time()-runtime_i)/60., 'min'    
    
    # define time related variables
    n_time = long(total_time*sample_rate)
    print '[lib_LBscan.py,gen_scan_c_mod] The number of samples in a day', n_time
    sec2date = (1.e-4/8.64)
    time_julian = np.arange(n_time)/(n_time-1.)*total_time*sec2date + today_julian # time in julian [day]

    # cal sun position
    print '[lib_LBscan.py,gen_scan_c_mod] before sun_position_quick', (time.time()-runtime_i)/60., 'min'    
    DJD = mylib.convert_Julian2Dublin(time_julian)

    lon_spe, lat_spe = sun_position_quick(DJD)
    ra_mp, dec_mp = moon_position(DJD)
    print today_julian, time_julian, DJD

    DJD = 0;

    # convert ra/dec to ecliptic coordinate
#    lon_spe, lat_spe = mylib.EULER( ra_sp*radeg, dec_sp*radeg, 3, 'J2000')
    lon_mpe, lat_mpe = mylib.EULER( ra_mp*radeg, dec_mp*radeg, 3, 'J2000')
    ra_mp = 0; dec_mp = 0
   
    time_i = time_julian /sec2date # time in [sec]
    
    # convert the sun vector to anti-sun vector
#    lat_aspe = -lat_spe/radeg
#    lon_aspe = lon_spe/radeg + pi
    lat_aspe = -lat_spe
    lon_aspe = lon_spe + pi
    
    lat_spe = 0
    lon_spe = 0
    
    # from ecliptic lat, lon convention to theta, phi convention
    theta_asp = pi/2.-lat_aspe
    phi_asp = lon_aspe

    x_asp,y_asp,z_asp = ang2pos_3D(theta_asp,phi_asp)

    phi_antisun = mylib.wraparound_2pi(2.*pi*freq_antisun*time_i)
    phi_boresight = mylib.wraparound_2pi(2.*pi*freq_boresight*time_i)

    print '[lib_LBscan.py,gen_scan_c_mod] before LB_rotmatrix_multi', (time.time()-runtime_i)/60., 'min'    
    p_out = liblbscanc.LB_rotmatrix_multi(theta_asp, phi_asp,
                                          theta_antisun, phi_antisun,
                                          theta_boresight, phi_boresight)

    theta_out = np.arctan2(np.sqrt(p_out[0,:]**2+p_out[1,:]**2),p_out[2,:]) 
    phi_out = np.arctan2(p_out[1,:],p_out[0,:])+pi
    theta_asp=0; phi_asp=0;
    p_out=0; time_i=0;
    phi_antisun=0; phi_boresight=0;

    ipix = h.ang2pix( nside, theta_out, phi_out)
    
    nbPix_scan = len(ipix)-1
    nbPix = h.nside2npix(nside)
        
    dpvec = np.zeros((3,nbPix_scan))
    dtvec = np.zeros((3,nbPix_scan))
    dphivec = np.zeros((3,nbPix_scan))
    
    xp,yp,zp = ang2pos_3D(theta_out,phi_out)
    dpvec = ang2_scandirec_3D(xp,yp,zp)
    dtvec = ang2_deriv_theta_3D(theta_out,phi_out)
    dphivec = ang2_deriv_phi_3D(theta_out,phi_out)
    
    alpha = np.zeros(nbPix_scan)
    alpha = np.arccos( (dpvec[0,:]*dtvec[0,:]+dpvec[1,:]*dtvec[1,:]+dpvec[2,:]*dtvec[2,:])
                       /np.sqrt(dpvec[0,:]*dpvec[0,:]+dpvec[1,:]*dpvec[1,:]+dpvec[2,:]*dpvec[2,:])
                       /np.sqrt(dtvec[0,:]*dtvec[0,:]+dtvec[1,:]*dtvec[1,:]+dtvec[2,:]*dtvec[2,:]) )
    
    beta = np.zeros(nbPix_scan)
    beta = np.arccos( (dpvec[0,:]*dphivec[0,:]+dpvec[1,:]*dphivec[1,:]+dpvec[2,:]*dphivec[2,:])
                      /np.sqrt(dpvec[0,:]*dpvec[0,:]+dpvec[1,:]*dpvec[1,:]+dpvec[2,:]*dpvec[2,:])
                      /np.sqrt(dphivec[0,:]*dphivec[0,:]+dphivec[1,:]*dphivec[1,:]+dphivec[2,:]*dphivec[2,:]) )
    
    print '[lib_LBscan.py,gen_scan_c_mod] before summing up', (time.time()-runtime_i)/60., 'min'    
    
    ind = np.where((beta >= 90./180.*pi) & (beta < 180./180.*pi))
    alpha[ind[0]] = - alpha[ind[0]]

    ############################################################################3
    if option_gen_ptg == True:
        fname=dir_out+'/'+title+'_'+str(int(today_julian))+'_latlon_eclip'
        tmpstr = time_julian[0]-2400000.5
#        time_arr = np.arange(len(alpha)+1)
#        print len(theta_out), len(phi_out), len(alpha), len(time_arr)
#        np.savez(fname,time_julian[0],time_arr/float(sample_rate)+tmpstr*86400., pi/2.-theta_out, phi_out, alpha, beta)
        np.savez(fname,time_julian[0],float(sample_rate), pi/2.-theta_out, phi_out, alpha, beta)
        np.savez(fname+'_moon', lat_mpe/radeg, lon_mpe/radeg)
#        del(time_arr)
    del(time_julian)
    ############################################################################3

    beta=0; dphivec=0; dpvec=0; dtvec=0; theta_out=0; phi_out=0
    mapout = liblbscanc.Maps_summingup(nbPix, np.float_(ipix), alpha)
    print '[lib_LBscan.py,gen_scan_c_mod] before writing text file ['+str(option_gen_ptg)+']', (time.time()-runtime_i)/60., 'min'   
    ipix=0; alpha=0
    runtime_f = time.time()
    print '[lib_LBscan.py,gen_scan_c_mod] END of gen_scan_c_mod ', (time.time()-runtime_i)/60., 'min'    
    print ''
    return mapout


#; ----------------------------------------------------------------------------
def gen_scan_c_mod_moon_ellip( theta_antisun_long, theta_antisun_short, theta_boresight, freq_antisun, freq_boresight, \
                                   total_time, today_julian, sample_rate, dir_out, filename, \
                                   title, nside, runtime_i, option_gen_ptg):
    
    print '[lib_LBscan.py,gen_scan_c_mod] initial', (time.time()-runtime_i)/60., 'min'    
    
    # define time related variables
    n_time = long(total_time*sample_rate)
    print '[lib_LBscan.py,gen_scan_c_mod] The number of samples in a day', n_time
    sec2date = (1.e-4/8.64)
    time_julian = np.arange(n_time)/(n_time-1.)*total_time*sec2date + today_julian # time in julian [day]

    # cal sun position
    print '[lib_LBscan.py,gen_scan_c_mod] before sun_position_quick', (time.time()-runtime_i)/60., 'min'    
    DJD = mylib.convert_Julian2Dublin(time_julian)

    lon_spe, lat_spe = sun_position_quick(DJD)
    ra_mp, dec_mp = moon_position(DJD)
    print today_julian, time_julian, DJD

    DJD = 0;

    # convert ra/dec to ecliptic coordinate
#    lon_spe, lat_spe = mylib.EULER( ra_sp*radeg, dec_sp*radeg, 3, 'J2000')
    lon_mpe, lat_mpe = mylib.EULER( ra_mp*radeg, dec_mp*radeg, 3, 'J2000')
    ra_mp = 0; dec_mp = 0
   
    time_i = time_julian /sec2date # time in [sec]
    
    # convert the sun vector to anti-sun vector
#    lat_aspe = -lat_spe/radeg
#    lon_aspe = lon_spe/radeg + pi
    lat_aspe = -lat_spe
    lon_aspe = lon_spe + pi
    
    lat_spe = 0
    lon_spe = 0
    
    # from ecliptic lat, lon convention to theta, phi convention
    theta_asp = pi/2.-lat_aspe
    phi_asp = lon_aspe

    x_asp,y_asp,z_asp = ang2pos_3D(theta_asp,phi_asp)

    phi_antisun = mylib.wraparound_2pi(2.*pi*freq_antisun*time_i)
    phi_boresight = mylib.wraparound_2pi(2.*pi*freq_boresight*time_i)
#    theta_antisun_arr = 0.5*(theta_antisun_long - theta_antisun_short) * np.cos(mylib.wraparound_2pi(2.*pi*freq_antisun*time_i)) \
#        + 0.5*(theta_antisun_long+theta_antisun_short)

    theta_antisun_arr = 0.5*(70./180.*pi - 5./180.*pi) * np.cos(mylib.wraparound_2pi(2.*pi*freq_antisun*time_i)) \
        + 0.5*(70./180.*pi + 5./180.*pi)
    theta_boresight = 0.

    print '[lib_LBscan.py,gen_scan_c_mod] before LB_rotmatrix_multi', (time.time()-runtime_i)/60., 'min'    
    p_out = liblbscanc.LB_rotmatrix_multi_ellip(theta_asp, phi_asp,
                                                theta_antisun_arr, phi_antisun,
                                                theta_boresight, phi_boresight)

    theta_out = np.arctan2(np.sqrt(p_out[0,:]**2+p_out[1,:]**2),p_out[2,:]) 
    phi_out = np.arctan2(p_out[1,:],p_out[0,:])+pi
    theta_asp=0; phi_asp=0;
    p_out=0; time_i=0;
    phi_antisun=0; phi_boresight=0;

    ipix = h.ang2pix( nside, theta_out, phi_out)
    
    nbPix_scan = len(ipix)-1
    nbPix = h.nside2npix(nside)
        
    dpvec = np.zeros((3,nbPix_scan))
    dtvec = np.zeros((3,nbPix_scan))
    dphivec = np.zeros((3,nbPix_scan))
    
    xp,yp,zp = ang2pos_3D(theta_out,phi_out)
    dpvec = ang2_scandirec_3D(xp,yp,zp)
    dtvec = ang2_deriv_theta_3D(theta_out,phi_out)
    dphivec = ang2_deriv_phi_3D(theta_out,phi_out)
    
    alpha = np.zeros(nbPix_scan)
    alpha = np.arccos( (dpvec[0,:]*dtvec[0,:]+dpvec[1,:]*dtvec[1,:]+dpvec[2,:]*dtvec[2,:])
                       /np.sqrt(dpvec[0,:]*dpvec[0,:]+dpvec[1,:]*dpvec[1,:]+dpvec[2,:]*dpvec[2,:])
                       /np.sqrt(dtvec[0,:]*dtvec[0,:]+dtvec[1,:]*dtvec[1,:]+dtvec[2,:]*dtvec[2,:]) )
    
    beta = np.zeros(nbPix_scan)
    beta = np.arccos( (dpvec[0,:]*dphivec[0,:]+dpvec[1,:]*dphivec[1,:]+dpvec[2,:]*dphivec[2,:])
                      /np.sqrt(dpvec[0,:]*dpvec[0,:]+dpvec[1,:]*dpvec[1,:]+dpvec[2,:]*dpvec[2,:])
                      /np.sqrt(dphivec[0,:]*dphivec[0,:]+dphivec[1,:]*dphivec[1,:]+dphivec[2,:]*dphivec[2,:]) )
    
    print '[lib_LBscan.py,gen_scan_c_mod] before summing up', (time.time()-runtime_i)/60., 'min'    
    
    ind = np.where((beta >= 90./180.*pi) & (beta < 180./180.*pi))
    alpha[ind[0]] = - alpha[ind[0]]

    ############################################################################3
    if option_gen_ptg == True:
        fname=dir_out+'/'+title+'_'+str(int(today_julian))+'_latlon_eclip'
        tmpstr = time_julian[0]-2400000.5
#        time_arr = np.arange(len(alpha)+1)
#        print len(theta_out), len(phi_out), len(alpha), len(time_arr)
#        np.savez(fname,time_julian[0],time_arr/float(sample_rate)+tmpstr*86400., pi/2.-theta_out, phi_out, alpha, beta)
        np.savez(fname,time_julian[0],float(sample_rate), pi/2.-theta_out, phi_out, alpha, beta)
        np.savez(fname+'_moon', lat_mpe/radeg, lon_mpe/radeg)
#        del(time_arr)
    del(time_julian)
    ############################################################################3

    beta=0; dphivec=0; dpvec=0; dtvec=0; theta_out=0; phi_out=0
    mapout = liblbscanc.Maps_summingup(nbPix, np.float_(ipix), alpha)
    print '[lib_LBscan.py,gen_scan_c_mod] before writing text file ['+str(option_gen_ptg)+']', (time.time()-runtime_i)/60., 'min'   
    ipix=0; alpha=0
    runtime_f = time.time()
    print '[lib_LBscan.py,gen_scan_c_mod] END of gen_scan_c_mod ', (time.time()-runtime_i)/60., 'min'    
    print ''
    return mapout



