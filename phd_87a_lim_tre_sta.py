# Dennis Alp 2018-04-09
# Gemini/T-ReCS, stack exposures 
# time python /Users/silver/box/bin/phd_87a_lim_tre_sta.py

import numpy as np
import os
import pdb
import time
from glob import glob
from datetime import date

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy.ndimage.interpolation import rotate
from scipy.ndimage.interpolation import zoom
from astropy.io import fits
from astropy.coordinates import SkyCoord

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

#########
# Help functions
# Just reasonable plot of sky images, good for debugging
def sky_plt(image):
    plt.imshow(image,interpolation='nearest', cmap='afmhot', origin='lower')
    plt.show()

# Utility for making fits image out of a numpy array
def mk_fits(image, output):
    hdu = fits.PrimaryHDU(image)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(output, clobber=True)
    hdulist.close()

def bin_plt(xx, yy):
    nx = np.max(xx)
    ny = np.max(yy)
    img, xedges, yedges = np.histogram2d(xx, yy, bins=bin_cou, range=src_loc)
    plt.imshow(img, interpolation='nearest', origin='low')
    plt.show()

def get_cen(xx, yy):
    xi = (xx > src_loc[0][0]) & (xx < src_loc[0][1])
    yi = (yy > src_loc[1][0]) & (yy < src_loc[1][1])
    return np.mean(xx[xi]), np.mean(yy[yi])

def coords2pix(ff, ra, de):
    raref = fits.getval(ff,'TCRVL11', 1)
    radel = fits.getval(ff,'TCDLT11', 1)
    rapix = fits.getval(ff,'TCRPX11', 1)
    deref = fits.getval(ff,'TCRVL12', 1)
    dedel = fits.getval(ff,'TCDLT12', 1)
    depix = fits.getval(ff,'TCRPX12', 1)
    x0 = (ra-raref)*np.cos(np.deg2rad(de))/radel+rapix-1
    y0 = (de-deref)/dedel+depix-1
    return x0, y0

# Fit ellipses for all observations
def help_ell(dummy, x0, y0, aa, bb, alpha, mag, sig, tilt, phi):
#     aa, bb, y0, x0, alpha, mag, sig = AA, BB, Y0, X0, ALPHA, 1., 1.
    xx = np.arange(0, nx)
    yy = np.arange(0, ny)[:,None]
    xx = xx-x0
    yy = yy-y0
    uu = xx*np.cos(alpha)-yy*np.sin(alpha)
    vv = xx*np.sin(alpha)+yy*np.cos(alpha)
    rr = (uu/aa)**2 + (vv/bb)**2
#    print x0, y0, aa, bb, alpha, mag, sig, tilt, phi
    global last_ell
    last_ell = abs(mag)*np.exp(-np.abs((np.sqrt(rr)-1)/sig)**2)*(1+tilt*np.cos(np.arctan2(vv,uu)+phi))
#    pdb.set_trace()
    return last_ell.ravel()

#########
# Parameters
# Logistics
path = '/Users/silver/dat/gem/87a/tre/red/'
out_dir = '/Users/silver/dat/gem/87a/tre/pro/'
out = ['001_2006-12-22_n__.fits']
groups = [
    ['rS20061223S0053.fits', 'rS20061222S0061.fits', 'rS20061222S0060.fits', 'rS20061222S0059.fits', 'rS20061222S0058.fits', 'rS20061222S0057.fits']
    ]
coords = [
    [[181, 104], [181, 106], [180, 104], [177, 104], [176, 103], [174, 103]]
    ]

#For flux calibration one can directly compare T-ReCS raw difference
#images that have the same time on source per saveset. If the images
#differ in the FRMTIME, FRMCOADD, or CHPCOADD, values then one would
#need to scale to a common saveset exposure time before comparing
#images.
standard = [['HD45669', 'HD42540', 'HD42540', 'HD42540', 'HD42540', 'HD42540']]
flux = {'HD42540': 5.745000/1.944E6,
        'HD45669': 6.694000/2.316E6}
    
# Allocations
fov = 50 # even
img = np.zeros((fov, fov))



################################################################
# Main
for ii, group in enumerate(groups):
    for jj, ff in enumerate(group):
        dat = fits.open(path + ff)
        print('FRMTIME =', dat[0].header['FRMTIME'], '\tFRMCOADD = ', dat[0].header['FRMCOADD'], '\tCHPCOADD = ', dat[0].header['CHPCOADD'])
        dat = dat[1].data

        xmin = int(coords[ii][jj][0] - fov/2)
        xmax = int(coords[ii][jj][0] + fov/2)
        ymin = int(coords[ii][jj][1] - fov/2)
        ymax = int(coords[ii][jj][1] + fov/2)

#        x0 = np.average(np.linspace(0, fov-1), weights=np.sum(dat[ymin:ymax, xmin:xmax], axis=0)) + xmin
#        y0 = np.average(np.linspace(0, fov-1), weights=np.sum(dat[ymin:ymax, xmin:xmax], axis=1)) + ymin
        ny, nx = dat.shape
        guess = np.array([coords[ii][jj][0], coords[ii][jj][1], 8, 6, 0.1, 40, 0.4, 0., 0.])
        pars, covar = curve_fit(help_ell, np.arange(0, dat.size), dat.ravel()-np.mean(dat), guess, maxfev=20000)
        x0 = pars[0]
        y0 = pars[1]
        print(x0, y0)
        
        xmin = int(np.round(x0 - fov/2))
        xmax = int(np.round(x0 + fov/2))
        ymin = int(np.round(y0 - fov/2))
        ymax = int(np.round(y0 + fov/2))
        i2jy = flux[standard[ii][jj]] # Intensity to Jansky
        img += dat[ymin:ymax, xmin:xmax]*i2jy
        

    img /= (jj+1)
    mk_fits(img, out_dir + out[ii])

# Show residual of stacked fit
ny, nx = img.shape
guess = np.array([fov/2, fov/2, 8, 6, 0.1, 1e-4, 0.4, 0., 0.])
pars, covar = curve_fit(help_ell, np.arange(0, img.size), img.ravel()-np.mean(img), guess, maxfev=20000)
res = img-last_ell
sky_plt(res)
pdb.set_trace()
    
