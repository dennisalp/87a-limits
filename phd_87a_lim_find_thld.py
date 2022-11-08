# Dennis Alp 2017-03-28
# Find the threshold/sigma for ALMA images
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/Dropbox/bin/phd_87a_uplim_find_thld.py

import numpy as np
import os
import pdb
import time
from glob import glob
from datetime import date

import matplotlib.pyplot as plt
import scipy.stats as scs
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import rotate
from scipy.ndimage.interpolation import zoom
from pyraf import iraf
from astropy.io import fits

#########
# Parameters
# Logistics
WRK_DIR = "/Users/silver/Dropbox/phd/projects/87a/uplim/tmp"
os.chdir(WRK_DIR) #Move to designated directory
files = sorted(glob('c23*.fits')) #Find all files.
FSTARS = 'stars.dat'
DECONV = 'DECONV.fits'

sig = np.linspace(2e-5, 40e-5, 60)
det = np.zeros(sig.shape[0])

# Utility for making fits image out of a numpy array
def mk_fits(image, output):
    hdu = fits.PrimaryHDU(image)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(output, clobber=True)
    hdulist.close()

for key in files:
    for idx, ss in enumerate(sig):
        scale = fits.getval(key,'CDELT2')
        minor = fits.getval(key, 'BMIN')/scale
        major = fits.getval(key, 'BMAJ')/scale
        fwhm = (major+minor)/2

        angle = fits.getval(key,'BPA')

        img = fits.open(key)[0].data[0,0]
        img = rotate(img, angle, axes=(1, 0), reshape=False, order=1, mode='constant', cval=0.0, prefilter=False)
        img = zoom(img, (1, major/minor), order=1, mode='constant', cval=0.0, prefilter=False)
        mk_fits(img, DECONV)
        
        iraf.daofind(DECONV, output=FSTARS, fwhmpsf=fwhm, sigma=ss, threshold=3, datamin=-0.1, datamax=9999, verify=False, Stdout=1)
        stscoo=iraf.pdump(FSTARS,'XCENTER,YCENTER','yes',Stdout=1)
        det[idx] = len(stscoo)

    dim = fits.open(key)[0].data[0,0].shape
    expected = dim[0]*dim[1]/fwhm**2*scs.norm(0,1).sf(3)*2
    print key, fwhm, expected
    print griddata(det, sig, expected)
    plt.semilogy(sig, det+1)
    plt.savefig(key.replace('fits','pdf'), bbox_inches='tight', pad_inches=0.01)
    plt.show()
