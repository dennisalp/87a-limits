# Dennis Alp 2018-04-17
# Gets the beam-corrected limits. This does not require the limits to
# be from the same point.
# time python /Users/silver/box/bin/phd_87a_lim_get_lim_alm.py

import numpy as np
import os
import pdb
from glob import glob
from datetime import date

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import shift
from astropy.io import fits

#########
# Parameters
# Logistics
WRK_DIR = '/Users/silver/dat/alm/87a/lim/'
BOX = [577,597,582,602]
limits = 'limits_alma_new.npy'

os.chdir(WRK_DIR) #Move to designated directory
files = sorted(glob('alma*.fits')) #Put ALMA last
limits = np.load(limits)

# Source location 0-indexed
XCOO = 592.158345778
YCOO = 586.775275464
VKICK = 0.8e6
distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e16 # m SI!!!

#########
# Help functions
# Just reasonable plot of sky images, good for debugging
def sky_plt(image):
    plt.imshow(image,interpolation='nearest', cmap='afmhot', origin='lower')
    plt.show()

# Days for scaling
SNDATE = date(1987, 2, 23)
def get_days(fname): # 1994-09-24_drz_al2.fits
    yr = int(fname[-23:-19])
    month = int(fname[-18:-16].lstrip("0"))
    day = int(fname[-15:-13].lstrip("0"))
    return (date(yr, month, day)-SNDATE).days*24*60*60

def inside_kick(yy, xx, key):
    telapsed = get_days(key)
    offset = np.rad2deg(np.arctan(VKICK*telapsed/(distance*1e3*pc)))*60*60*1000/25
    
    dely = (BOX[0]+yy/2)-YCOO
    delx = (BOX[2]+xx/2)-XCOO
    return np.sqrt(delx**2+dely**2) < offset

max_val = [0, 0, 0]
for ii, key in enumerate(files): # Colors
    for yy in range(0, limits.shape[0]):
        for xx in range(0, limits.shape[1]):
            if inside_kick(yy, xx, key):
                tmp = limits[yy, xx, ii]
                if tmp > max_val[0]:
                    max_val = [tmp, xx, yy]
    print(key, max_val)
    max_val = [0, 0, 0]

pdb.set_trace()
