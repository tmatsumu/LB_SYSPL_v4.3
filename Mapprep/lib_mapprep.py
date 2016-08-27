import numpy as np
import pylab as py
import healpy as h

pi = np.pi
radeg = (180./pi)

def Antenna2Thermo(nu_GHz):
    Tcmb = 2.725
    h = 6.62606957e-34
    k_b = 1.3806488e-23
    nu = nu_GHz*1.e9

    x = (h*nu)/k_b/Tcmb
    return (np.exp(x)-1.)**2 / (x**2*np.exp(x))

def gen_dipole(nside,unit):
    # from WMAP 1st yaer Bennett et al
    l = 263.85 # +/-0.1 degs
    b = 48.25 # +/-0.04 degs
    amp_mK = 3.346 # +/-0.017 mK
    if unit == 'uK': amp_unit = amp_mK * 1.e3; print unit
    if unit == 'mK': amp_unit = amp_mK; print unit
    if unit == 'K': amp_unit = amp_mK * 1.e-3; print unit

    npix = h.nside2npix(nside)
    ipix = range(npix)
    theta, phi = h.pix2ang(nside,ipix)
    dipole = np.cos(theta)

    xyz = h.ang2vec(theta,phi)
    theta0 = pi/2.-b/radeg;  phi0 = l/radeg
    x = np.cos(phi0)*np.cos(theta0)*xyz[:,0] + np.sin(phi0)*np.cos(theta0)*xyz[:,1] - np.sin(theta0)*xyz[:,2]
    y = - np.sin(phi0)*xyz[:,0]+np.cos(phi0)*xyz[:,1]
    z = np.sin(theta0)*np.cos(phi0)*xyz[:,0] + np.sin(phi0)*np.sin(theta0)*xyz[:,1] + np.cos(theta0)*xyz[:,2]

    ipix = h.vec2pix(nside,x,y,z)

    dipole_out = amp_unit*dipole[ipix]
    return dipole_out

def db_convert(db):
    return 10.**(db/10.)

def gen_LBGaussBeam_FWHMarcmin(nu_GHz):
    if nu_GHz==60: return 75.
    if nu_GHz==78: return 60.
    if nu_GHz==100: return 45.
    if nu_GHz==140: return 32.
    if nu_GHz==195: return 23.
    if nu_GHz==280: return 16.

