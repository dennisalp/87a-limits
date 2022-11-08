# Dennis Alp 2018-01-22
# Just make a pretty plot of an ALMA image.
# time python /Users/silver/box/bin/phd_87a_uplim_plt_alm.py

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

def coords2pix(image, ra, de):
    raref = fits.getval(image,'CRVAL1')
    radel = fits.getval(image,'CDELT1')
    rapix = fits.getval(image,'CRPIX1')
    deref = fits.getval(image,'CRVAL2')
    dedel = fits.getval(image,'CDELT2')
    depix = fits.getval(image,'CRPIX2')
    xx = (ra-raref)*np.cos(np.deg2rad(YCOO))/radel+rapix-1
    yy = (de-deref)/dedel+depix-1
    return xx, yy

# Days for scaling
SNDATE = date(1987, 2, 23)
def get_days(yr, month, day):
    return (date(yr, month, day)-SNDATE).days*24*60*60

def mk_beam(img_path):
    scale = fits.getval(img_path,'CDELT2')
    minor = fits.getval(img_path,'BMIN')/scale
    major = fits.getval(img_path,'BMAJ')/scale
    angle = fits.getval(img_path,'BPA')
    return major, minor, angle

########
# Plots
def plt_blu(img):
    plt.figure()
    plt.imshow(img, origin='lower', interpolation='nearest', cmap='gray')
    plt.contour(np.log(blu)+2.8, levels=[0.5, 0.8, 1.2, 2, 3], cmap='viridis')
    plt.plot(-left, -down, '.k')

def plt_ibl(img):
    plt.figure()
    plt.imshow(np.log(blu)+2.8, origin='lower', interpolation='nearest', cmap='afmhot', vmin=0.5, vmax=2)
    plt.contour(img, levels=[0.005, 0.008, 0.011, 0.014, 0.017, 0.020, 0.023], cmap='viridis')
    plt.plot(-left, -down, '.k')
    
def pro_plt(img, img_path):
    fig = plt.figure(figsize=(5, 3.75))
    ax = fig.gca()
    
    #Maps
    plt.imshow(1000*img, origin='lower', interpolation='nearest', cmap='viridis', vmin=-0.11, vmax=0.16)
    cb = plt.colorbar()
    cb.set_label('Intensity (mJy beam$^{-1}$)')
    bl2 = blu.copy()
    bl2[:130, :30] = 1e-10
    bl2[:75, :75] = 1e-10
    bl2[:50, :100] = 1e-10
    plt.contour(np.log(bl2)+2.8, levels=[0.6, 1.6, 3], colors='w', linewidths=0.75)
    
    # Compass
    ax.arrow(11, 11, 0, 50, head_width=4, head_length=4, fc='w', ec='w', linewidth=1)
    ax.arrow(11, 11, 50, 0, head_width=4, head_length=4, fc='w', ec='w', linewidth=1)
    ax.annotate("W", xy=(70, 7), color='w')
    ax.annotate("N", xy=(5, 72), color='w')

    # Scale
    plt.plot([360, 376], [300, 300], 'w', lw=2)
    plt.annotate("$0.''1$", xy=(353, 308), color='w')
    
    # Beam
    major, minor, angle = mk_beam(img_path)
    ax.add_patch(patches.Ellipse((370, 30), minor, major, angle=angle, facecolor='w', linewidth=0, zorder=3))
#    plt.plot(-left, -down, '.w')

    # Search region
    day = get_days(2014, 9, 2)
    rad = 8e7*day/(51.2*3.086e21)
    rad = rad/(2*np.pi)*360*60*60*1000
    ellipse = patches.Ellipse(xy=(-left, -down), width=2*rad/6, height=2*rad/6, edgecolor='w', fc='None', lw=2, ls=':')
    ax.add_artist(ellipse)

    ax.set_axis_off()
    fig.savefig('/Users/silver/box/phd/pro/87a/lim/art/figs/alm.pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)
    
#########
# Parameters
# Logistics
a21_path = '/Users/silver/dat/alm/87a/larsson/alm21_2014-09-02_sum_una.fits'
a54_path = '/Users/silver/dat/alm/87a/larsson/alm54_2014-09-02_sum_una.fits'
a65_path = '/Users/silver/dat/alm/87a/larsson/alm65_2014-09-02_sum_una.fits'
a47_path = '/Users/silver/dat/alm/87a/lim/alma_2014-09-02_247_247.fits'
a33_path = '/Users/silver/dat/alm/87a/lim/alma_2014-09-02_233_233.fits'
a13_path = '/Users/silver/dat/alm/87a/lim/alma_2014-09-02_212_213.fits'
key_path = '/Users/silver/dat/alm/87a/lim/alma_2014-09-02_247_247.fits'
blu_path = '/Users/silver/dat/hst/87a/all/w3_b_2014-06-15_drz_al2.fits'
red_path = '/Users/silver/dat/hst/87a/all/w3_r_2014-06-15_drz_al2.fits'

cc = SkyCoord('05h35m27.9875s', '-69d16m11.107s', frame='icrs')
CRAB = SkyCoord('05h34m31.97232s', '+22d00m52.0690s', frame='fk5')
CASA = SkyCoord('23h23m27.943s', '58d48m42.51s', frame='fk5')
CRAX = CRAB.icrs.ra.degree
CRAY = CRAB.icrs.dec.degree
XCOO = cc.ra.degree
YCOO = cc.dec.degree

#XCOO = cc.fk5.ra.degree
#YCOO = cc.fk5.dec.degree
VKICK = 0.8e6
distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e18
mag = 25/6.    

down = -170
up = 170
left = -205
right = 200

# Get the coordinates, Hubble and ALMA
xh = np.int(np.round(592.118767092*mag))
yh = np.int(np.round(586.714283999*mag))
xa, ya = coords2pix(key_path, XCOO, YCOO)
xa = np.int(np.round(xa))
ya = np.int(np.round(ya))


#########
# Load all data
a21 = fits.open(a21_path)[0].data[ya+down:ya+up, xa+left:xa+right]
a54 = fits.open(a54_path)[0].data[ya+down:ya+up, xa+left:xa+right]
a65 = fits.open(a65_path)[0].data[ya+down:ya+up, xa+left:xa+right]
a47 = fits.open(a47_path)[0].data[0,0, ya+down:ya+up, xa+left:xa+right]
a33 = fits.open(a33_path)[0].data[0,0, ya+down:ya+up, xa+left:xa+right]
a13 = fits.open(a13_path)[0].data[0,0, ya+down:ya+up, xa+left:xa+right]
blu = fits.open(blu_path)[0].data
blu = zoom(blu, mag)[yh+down:yh+up, xh+left:xh+right]
red = fits.open(red_path)[0].data
red = zoom(red, mag)[yh+down:yh+up, xh+left:xh+right]
    
#plt_blu(a21)
#plt_blu(a54)
#plt_blu(a65)
##plt_blu(a47)
##plt_blu(a33)
##plt_blu(a13)
#
#plt_ibl(a21)
#plt_ibl(a54)
#plt_ibl(a65)
##plt_ibl(a47)
##plt_ibl(a33)
##plt_ibl(a13)

pro_plt(a13, a13_path)

plt.show()
pdb.set_trace()
