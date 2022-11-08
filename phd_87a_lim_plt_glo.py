# Dennis Alp 2016-01-27
# Find models that are consistent with oberved limits.
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/Dropbox/bin/phd_87a_uplim_plt_lim.py

import numpy as np
import os
import pdb
from glob import glob
from datetime import date

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import shift
from pyraf import iraf
from astropy.io import fits

from phd_87a_red import cor_red

################################################################
#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

################################################################
# Parameters
# Logistics
WRK_DIR = "/Users/silver/Dropbox/phd/projects/87a/uplim/glo"
os.chdir(WRK_DIR) #Move to designated directory

distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e16
Lsun = 3.826e33 # erg s-1
cc = 299792458 # m s-1
kpc = 3.086e21 # cm

# Redenning
RV = 3.1
EBV = 0.19

################################################################
# Help functions
def lam2nu(lam):
    return cc/lam

################################################################
# Data
nu = []
fn = []
co = []

########
# callingham16, entire remnant
tmp = 1e9*np.array([0.076, 0.084, 0.092, 0.099, 0.107, 0.115, 0.123, 0.130, 0.143, 0.150, 0.158, 0.166, 0.174, 0.181, 0.189, 0.197, 0.204, 0.212, 0.219, 0.227, 1.375, 1.375, 2.351, 2.351, 4.788, 4.788, 8.642, 8.642])
nu.append(tmp)
tmp = 1e-23*np.array([5.1, 4.9, 4.7, 4.6, 4.5, 4.2, 4.0, 3.9, 3.6, 3.4, 3.3, 3.1, 3.0, 2.9, 2.8, 2.7, 2.5, 2.5, 2.4, 2.3, 0.58, 0.58, 0.43, 0.42, 0.28, 0.30, 0.18, 0.17])
fn.append(tmp)
co.append(0*np.ones(len(tmp)))

# zanardo14, estimate
tmp = 1e9*np.linspace(102, 672, 1000)
nu.append(tmp)
tmp = 3e-26*np.ones(1000)
fn.append(tmp)
co.append(1*np.ones(len(tmp)))

# potter09, estimate 0.3 mJy, 3sigma 0.9 mJy
tmp = [36.2e9]
nu.append(tmp)
tmp = [0.3e-26]
fn.append(tmp)
co.append(2*np.ones(len(tmp)))

# ng08, 3sigma limit
tmp = [9e9]
nu.append(tmp)
tmp = [0.3e-26]
fn.append(tmp)
co.append(3*np.ones(len(tmp)))

# lakicevic12, 2sigma limit
tmp = [94e9]
nu.append(tmp)
tmp = [1e-26]
fn.append(tmp)
co.append(4*np.ones(len(tmp)))

# zanardo13, loose upper limit
tmp = [44e9]
nu.append(tmp)
tmp = [2.2e-26]
fn.append(tmp)
co.append(5*np.ones(len(tmp)))

# matsuura15, entire remnant, dust except 500 microns, which is 3sigma limit
tmp = lam2nu(1e-6*np.array([70, 100, 160, 250, 350, 500]))
nu.append(tmp)
tmp = 1e-26*np.array([45.4, 82.4, 153., 110.7, 69.3, 60.])
fn.append(tmp)
co.append(6*np.ones(len(tmp)))

# arendt16, entire remnant
tmp = lam2nu(1e-6*np.array([3.6, 4.5, 5.8, 8., 24.]))
nu.append(tmp)
tmp = 1e-26*np.array([1.52, 2.17, 4.08, 13.61, 75.7])
fn.append(tmp)
co.append(7*np.ones(len(tmp)))

# frank16, 90% limit
tmp = 1e17*np.linspace(4.836, 24.18, 1000)
nu.append(tmp)
tmp = 3.1e33/(4*np.pi*(51.2*kpc)**2)
tmp = tmp/(2*(np.sqrt(2.418e18)-np.sqrt(4.836e17)))
tmp = tmp*np.array(nu[-1])**-0.5
fn.append(tmp)
co.append(8*np.ones(len(tmp)))
#print np.sum(np.array(fn[-1])*(nu[-1][1]-nu[-1][0]))*4*np.pi*(51.2*kpc)**2

# grebenev12, entire remnant
tmp = 1e18*np.linspace(4.836, 14.51, 1000)
nu.append(tmp)
tmp = 3e35/(4*np.pi*(51.2*kpc)**2)
tmp = tmp/(-10*(1.451e19**-0.1-4.836e18**-0.1))
tmp = tmp*np.array(nu[-1])**-1.1
fn.append(tmp)
co.append(9*np.ones(len(tmp)))

########
# New limits
# ALMA
tmp = [213e9]
nu.append(tmp)
tmp = [4.1e-26]
fn.append(tmp)
co.append(10*np.ones(len(tmp)))

# SINFONI
tmp = cc*1e10/np.array([23275., 21300., 18262.5, 17512.5, 15475.])
nu.append(tmp)
tmp = 1e-19*np.array([1.2, 1.8, 0.83, 1.4, 2.5])*np.array([23275., 21300., 18262.5, 17512.5, 15475.])/tmp
fn.append(tmp)
co.append(11*np.ones(len(tmp)))

# STIS
tmp = cc/(1e-10*np.linspace(5300, 10000, 1000))
nu.append(tmp)
tmp = 7.80191789426e-15*np.linspace(5300, 10000, 1000)**1.05/(cc*1e10)
fn.append(tmp)
co.append(11*np.ones(len(tmp)))

# WFC3/UVIS
tmp = cc*1e10/np.array([6255.5, 4330.5, 8074., 5334., 3359.5, 2382.5])
nu.append(tmp)
tmp = 1e-18*np.array([5.5, 7.1, 2.8, 6.7, 9.5, 40.])*np.array([6255.5, 4330.5, 8074., 5334., 3359.5, 2382.5])/tmp
fn.append(tmp)
co.append(12*np.ones(len(tmp)))

########
# finalize
nu = np.concatenate(nu)
fn = np.concatenate(fn)
co = np.concatenate(co)

################################################################
# Initiate plot
fig = plt.figure(figsize=(5, 3.75))
ax1 = fig.add_subplot(111)

ax1.loglog(nu, nu*fn, '.')
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("$\\nu F_\\nu$ (erg s$^{-1}$ cm$^{-2}$)")

# This is wrong!
#ax2 = ax1.twiny()
#def tick_function(X):
#    V = 1e10*cc/X
#    return ["%.3f" % z for z in V]
#
##new_tick_locations = np.linspace(0,1,ax1.get_xticks().size)
##ax2.set_xscale(ax1.get_xscale())
##ax2.set_xlim(new_tick_locations)
##ax2.set_xticks(tick_function(new_tick_locations))
##ax2.set_xticklabels(ax1.get_xticklabels())
#ax2.set_xscale(ax1.get_xscale())
#ax2.set_xlim(ax1.get_xlim())
#ax2.set_xticks(ax1.get_xticks())
#ax2.set_xticklabels(1e10*3e8/ax1.get_xticks())
#import matplotlib.ticker as ticker
#ax2.xaxis.set_major_formatter(ticker.LogFormatterMathtext())
#ax2.set_xlabel("Wavelength (\AA{})")

fig.savefig('/Users/silver/Dropbox/phd/projects/87a/uplim/art/figs/glo_lim.pdf',bbox_inches='tight', pad_inches=0.03)
plt.show()
pdb.set_trace()
