import numpy as np
import pylab as py

pi = np.pi
radeg = (180./pi)

def boresight2pixel(theta_b,phi_b,dk_b,r_fp,theta_fp):
    phi = np.atan(-np.sin(theta_fp+dk_b)*np.sin(r_fp),np.cos(theta_b)*np.cos(theta_fp+dk_b)*np.sin(r_fp)+np.sin(theta_b)*np.cos(r_fp))
    theta = -np.sin(theta_b)*np.cos(theta_fp+dk_b)*np.sin(r_fp)+np.cos(theta_b)*np.cos(r_fp)
    return theta, phi

r_fp = 10./radeg
theta_fp = 30./radeg

theta_b = 30./radeg
phi_b = 50./radeg
dk_b = 0./radeg


