import numpy as np
import pylab as py
import healpy as h
from scipy.optimize import fmin

def chisq_nonoise(data,model):
    chi = np.sum((data-model)**2)
    return chi

def chisq_noise(data,model,sigma):
    chi = np.sum((data-model)**2/sigma**2)
    return chi

def chisq_dipole_gain_fit(par,tod,model):
    chisq = tod - par*model
    return chisq

def dipole_gain_anafit(top_tod, bot_tod, top_dipole, bot_dipole):
    
    top_par = np.sum(top_tod*top_dipole)/np.sum(top_dipole**2)
    bot_par = np.sum(bot_tod*bot_dipole)/np.sum(bot_dipole**2)

#    np.savez('tmp.npz',top_tod=top_tod,top_par=top_par,top_dipole=top_dipole, \
#                 bot_tod=bot_tod,bot_par=bot_par,bot_dipole=bot_dipole )
#    sys.exit()

    return top_tod/top_par-top_dipole, bot_tod/bot_par-bot_dipole
