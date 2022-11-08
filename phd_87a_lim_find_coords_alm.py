# Dennis Alp 2016-11-08
# Fit stuff to SN 1987A remnant and determine center, ALMA
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python -i /Users/silver/box/bin/phd_87a_uplim_find_coords_alm.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
import os
import time
import pdb

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units
from scipy.optimize import curve_fit
from scipy.interpolate import griddata
from scipy.stats import sem
from glob import glob

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

#########
# Parameters
DAT_DIR = "/Users/silver/box/phd/data/alm/87a/larsson/"
WRK_DIR = "/Users/silver/box/phd/projects/87a/uplim/"
RINGEST = "/Users/silver/box/phd/data/alm/87a/larsson/co21_sum.fits"
THICKNESS = 62.5 # Thickness of Gaussian radial profile of ellipse
NX = 600
NY = 600
INTERPOLATION = 'nearest'

os.chdir(WRK_DIR) #Move to designated directory
files = glob(DAT_DIR+'/*sum.fits') #Find all files.
NN = len(files)
coords = []
boxes = np.array([140, 475, 120, 540])

#########
# Ellipse parameters http://math.stackexchange.com/questions/426150/what-is-the-general-\
# equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate
X0 = 324
Y0 = 300
AA = 117
BB = 83
ALPHA = 0.15
subnx = boxes[3]-boxes[2]
subny = boxes[1]-boxes[0]
 
def gen_mask(AA, BB, X0, Y0, ALPHA):
    xx = np.arange(0,NX)
    yy = np.arange(0,NY)[:,None]
    xx = xx-X0
    yy = yy-Y0
    uu = xx*np.cos(ALPHA)-yy*np.sin(ALPHA)
    vv = xx*np.sin(ALPHA)+yy*np.cos(ALPHA)
    return ((uu/AA)**2 + (vv/BB)**2 > 0.5) & ((uu/AA)**2 + (vv/BB)**2 < 2.7) # True for points inside
mask = gen_mask(AA, BB, X0, Y0, ALPHA)
 
#########
# Fit ellipses for all observations
def help_ell(dummy, x0, y0, aa, bb, alpha, mag, sig, tilt, phi):
#     aa, bb, y0, x0, alpha, mag, sig = AA, BB, Y0, X0, ALPHA, 1., 1.
    xx = np.arange(0,subnx)
    yy = np.arange(0,subny)[:,None]
    xx = xx-x0
    yy = yy-y0
    uu = xx*np.cos(alpha)-yy*np.sin(alpha)
    vv = xx*np.sin(alpha)+yy*np.cos(alpha)
    rr = (uu/aa)**2 + (vv/bb)**2
#     print x0, y0, aa, bb, alpha, mag, sig, tilt, phi
    global last_ell
    last_ell = abs(mag)*np.exp(-np.abs(np.sqrt(rr)-1)**2/abs(sig))*(1+tilt*np.cos(np.arctan2(vv,uu)+phi))
    return np.reshape(last_ell, subnx*subny)
 
for img_path in files:
    img = fits.open(img_path)[0].data
    img = np.where(mask, img, 0.)[boxes[0]:boxes[1],boxes[2]:boxes[3]]
    volume = np.sum(np.abs(img))
 
# Fit
    guess = np.array([X0-boxes[2], Y0-boxes[0], AA, BB, ALPHA, 1., 0.05, 0.25, 0.])
    pars, covar = curve_fit(help_ell, np.arange(0, subnx*subny), np.reshape(img, subnx*subny), guess)
    pars[0] = pars[0] + boxes[2]
    pars[1] = pars[1] + boxes[0]
    res = np.sum(np.abs(img-last_ell)) / volume
    print("{0:17.13f} {1:17.13f} {2:17.13f}".format(pars[0], pars[1], res))
    coords.append([pars[0], pars[1]])
    print img_path, pars[2]*6, pars[3]*6
 
# Plots
    fig = plt.figure()
    plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(img)/100, vmax=np.amax(img)))
    fig.savefig('coord_fits/alm/' + img_path[-13:-5] + '_ell_img.pdf',bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)
    
    fig = plt.figure()
    plt.imshow(img-last_ell, cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
    fig.savefig('coord_fits/alm/' + img_path[-13:-5] + '_ell_res.pdf',bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)

    #########
# Gaussian parameters
X0 = 322
Y0 = 296
AA = 0.02 # Amplitude
SIG = 42.
gbox = np.array([229, 365, 260, 390])
subnx = gbox[3]-gbox[2]
subny = gbox[1]-gbox[0]
 
def gen_mask(X0, Y0, SIG):
    xx = np.arange(0,NX)
    yy = np.arange(0,NY)[:,None]
    xx = xx-X0
    yy = yy-Y0
    return np.sqrt(xx**2+yy**2) < 1.85*SIG
 
#########
# Fit Gaussians for all observations
def help_gauss(dummy, x0, y0, aa, sig):
    xx = np.arange(0,subnx)
    yy = np.arange(0,subny)[:,None]
    xx = xx-x0
    yy = yy-y0
    rr = xx**2+yy**2
    global last_gau
    last_gau = aa*np.exp(-rr/sig**2)
    return np.reshape(last_gau, subnx*subny)
 
for img_path in files:
# Gauss
    img = fits.open(img_path)[0].data
    mask = gen_mask(X0, Y0, SIG)
    img = np.where(mask, img, 0.)[gbox[0]:gbox[1],gbox[2]:gbox[3]]
    volume = np.sum(np.abs(img))
 
    guess = np.array([X0-gbox[2], Y0-gbox[0], AA, SIG])
    pars, covar = curve_fit(help_gauss, np.arange(0, subnx*subny), np.reshape(img, subnx*subny), guess)
    pars[0] = pars[0] + gbox[2]
    pars[1] = pars[1] + gbox[0]
    res = np.sum(np.abs(img-last_gau)) / volume
    print("{0:17.13f} {1:17.13f} {2:17.13f} Gauss".format(pars[0], pars[1], res))
    coords.append([pars[0], pars[1]])
 
# Plots
    fig = plt.figure()
    plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
    fig.savefig('coord_fits/alm/' + img_path[-13:-5] + '_gau_img.pdf',bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)
 
    fig = plt.figure()
    plt.imshow(img-last_gau, cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
    fig.savefig('coord_fits/alm/' + img_path[-13:-5] + '_gau_res.pdf',bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)
    
#########
# Plot and wrap-up
img = fits.open(RINGEST)[0].data[boxes[0]:boxes[1],boxes[2]:boxes[3]]
coords.append([330.20284-1, 294.53336-1])
coords = np.array(coords)

fig = plt.figure()
col = ['or','og','ob','oc','om','oy','ow']

for i in range(0,coords.shape[0]):
    plt.plot(coords[i,0]-boxes[2], coords[i,1]-boxes[0],col[i])

plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(img)/100, vmax=np.amax(img)))
fig.savefig('coord_fits/alm/' + RINGEST[-13:-5] + '_coords.pdf',bbox_inches='tight', pad_inches=0.01)
plt.show()
plt.close(fig)

ra = np.mean(coords[0:3,0])
dec= np.mean(coords[0:3,1])

dec= -6.9269742E+01+(dec-300)*1.666666666667E-06
ra = 8.38667500E+01-(ra-300) *1.666666666667E-06/np.cos(np.deg2rad(dec))
coo=SkyCoord(ra, dec, frame='icrs', unit='deg')
print coo.to_string('hmsdms')
