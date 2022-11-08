# Dennis Alp 2016-10-20
# Align images using IRAF/PyRAF
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/$USER/Dropbox/bin/sci_iraf_87a_fold.py

import numpy as np
import os
import pdb

from glob import glob
from pyraf import iraf
from astropy.io import fits
from scipy.misc import imresize
from scipy.interpolate import griddata
from datetime import date

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#########
# Parameters
DAT_DIR = "/Users/silver/Dropbox/phd/data/hst/87a/all/"
WRK_DIR = "/Users/silver/Dropbox/phd/projects/87a/uplim/fold/"

os.chdir(WRK_DIR) #Move to designated directory
files = glob(DAT_DIR+'/*.fits') #Find all files.
NN = len(files)
X0 = 593.119765996 # In pixels, so indexed from 1
Y0 = 587.715824525
XMAX = 1150
YMAX = 1115
SNDATE = date(1987, 2, 23)
REF_DATE = float((date(2016, 6, 8)-SNDATE).days)
INTERPOLATION = 'nearest'

XX = np.arange(0, XMAX)
YY = np.arange(0, YMAX)
XG, YG = np.meshgrid(XX, YY)
XG = XG.reshape(XMAX*YMAX)[:,np.newaxis]
YG = YG.reshape(XMAX*YMAX)[:,np.newaxis]
folded = np.zeros((YMAX, XMAX)) # Flipped dimension since it is an image

def get_days(fname): # 1994-09-24_drz_al2.fits
    yr = int(fname[-23:-19])
    month = int(fname[-18:-16].lstrip("0"))
    day = int(fname[-15:-13].lstrip("0"))
    return (date(yr, month, day)-SNDATE).days

for ff in files:#[::15]:
    days = get_days(ff)
    print ff, "DAY", days
    
    ff = fits.open(ff)[0].data
    mag = days/REF_DATE
#    plt.imshow(ff, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(ff)/1000000., vmax=np.amax(ff)))
#    plt.show()
    flux = np.sum(ff[510:540,790:820])
    print flux
    
    ff = griddata(np.hstack((XG, YG)), ff.reshape(XMAX*YMAX), np.hstack(((XG-X0)*mag+X0, (YG-Y0)*mag+Y0)), method='linear', fill_value=0)
    ff = ff.reshape((YMAX, XMAX)) # Flipped dimension since it is an image
#    plt.imshow(ff, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(ff)/1000000., vmax=np.amax(ff)))
#    plt.show()

#    normalization = np.sum(ff)/flux
#    print normalization
    folded = folded+ff/flux#/normalization
#    plt.imshow(folded, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(ff)/1000000., vmax=np.amax(ff)))
#    plt.show()

hdu = fits.PrimaryHDU(folded)
hdulist = fits.HDUList([hdu])
hdulist.writeto('folded_hst_normed.fits')
