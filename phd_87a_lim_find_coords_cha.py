# Dennis Alp 2017-07-10
# Fit stuff to SN 1987A remnant and determine center, Chandra
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python -i /Users/silver/box/bin/phd_87a_uplim_find_coords_cha.py

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
from scipy.ndimage.interpolation import zoom
from scipy.stats import sem
from scipy.misc import factorial
import scipy.stats as sts
from glob import glob

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

#########
# Parameters
WRK_DIR = "/Users/silver/box/phd/pro/87a/lim/cha/16756/repro"
img_path = "acisf16756_repro_img.fits"
psf_path = "sim.psf"
THICKNESS = 15 # Thickness of Gaussian radial profile of ellipse
NX = 640
NY = 640
INTERPOLATION = 'nearest'

os.chdir(WRK_DIR) #Move to designated directory

#########
# Ellipse parameters http://math.stackexchange.com/questions/426150/what-is-the-general-\
# equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate
X0 = 320
Y0 = 320
AA = 15
BB = 14
ALPHA = 0.15
smajor = 9.154668411754782
sminor = 6.659880901513138
 
#########
# Fit ellipses for all observations
def gen_mask():
    xx = np.arange(0, NX)
    yy = np.arange(0, NY)[:,None]
    xx = xx-3.18322892e+02 
    yy = yy-3.21661990e+02 
    alpha = 1.30411412e-01
    uu = xx*np.cos(alpha)-yy*np.sin(alpha)
    vv = xx*np.sin(alpha)+yy*np.cos(alpha)
    aa = smajor
    bb = sminor
    rr = (uu/aa)**2 + (vv/bb)**2
    return rr > 1

def help_ell(dummy, x0, y0, aa, bb, alpha, mag, sig, tilt, phi):
#     aa, bb, y0, x0, alpha, mag, sig = AA, BB, Y0, X0, ALPHA, 1., 1.
    xx = np.arange(0,NX)
    yy = np.arange(0,NY)[:,None]
    xx = xx-x0
    yy = yy-y0
    uu = xx*np.cos(alpha)-yy*np.sin(alpha)
    vv = xx*np.sin(alpha)+yy*np.cos(alpha)
    rr = (uu/aa)**2 + (vv/bb)**2
#    print x0, y0, aa, bb, alpha, mag, sig, tilt, phi
    global last_ell, fold_ell
    last_ell = np.exp(-np.abs(np.sqrt(rr)-1)**2/abs(sig))*(1+tilt*np.cos(np.arctan2(vv,uu)+phi))

    # Now forward fold it with the psf
    last_ell = np.where(last_ell > 0, last_ell, 0)
    tstfld = np.zeros(last_ell.shape)
    for xi in range(300,340):
        for yi in range(300,340):
            tstfld[yi-nn:yi+nn+1, xi-nn: xi+nn+1] += last_ell[yi, xi]*psf
    bkgimg = img.copy()
    bkgimg[210:430, 210:430] = 0
    bkgflx = np.sum(bkgimg)/(640.**2-220.**2)
    tstfld = tstfld/np.sum(tstfld)*np.sum(img-bkgflx)+bkgflx
    fold_ell = tstfld
    return np.reshape(np.where(mask, tstfld, 0), NX*NY)



# Utility for making fits image out of a numpy array
def mk_fits(image, output):
    hdu = fits.PrimaryHDU(image)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(output, clobber=True)
    hdulist.close()
    
def poisson(k, lamb):
    """poisson pdf, parameter lamb is the fit parameter"""
    return (lamb**k/factorial(k)) * np.exp(-lamb)


def negLogLikelihood(model, data):
    lnl = -np.sum(np.log(poisson(data, model)))
    return lnl

def get_err():
    mag = 10
    err = pars.copy()
    help_ell(0, err[0],err[1],err[2],err[3],err[4],err[5],err[6],err[7],err[8])
    sav_ell = fold_ell.copy()
    nerr = 100
    x1 = np.zeros(nerr)
    y1 = np.zeros(nerr)
    fs = np.zeros(nerr)
    for ii in range(0, nerr):
        print ii
        poifld = sts.poisson.rvs(fold_ell)
        guess = np.array([err[0],err[1],err[2],err[3],err[4],err[5],err[6],err[7],err[8]])
        temp, covar = curve_fit(help_ell, np.arange(0, NX*NY), np.reshape(poifld, NX*NY), guess)
        x1[ii] = temp[0]
        y1[ii] = temp[1]

        zoomed = zoom(poifld, mag, order=1)
        xx = np.arange(0,mag*NX)
        yy = np.arange(0,mag*NY)[:,None]
        x0=mag*pars[0]
        y0=mag*pars[1]
        xx = xx-x0
        yy = yy-y0
        uu = xx*np.cos(temp[4])-yy*np.sin(temp[4])
        vv = xx*np.sin(temp[4])+yy*np.cos(temp[4])

        aa = mag*smajor
        bb = mag*sminor
        rr = (uu/aa)**2 + (vv/bb)**2
        fs[ii] = np.sum(np.where(rr <= 1, zoomed, 0))

    aa = 2032.7452595559723
    bb = 1478.7910083532158
    rr = (uu/aa)**2 + (vv/bb)**2
    outer = np.sum(np.where(rr <= 1, zoomed, 0))
    fs = fs/(outer-fs)

    np.save('fit_fs', fs)
    np.save('fit_xx', x1)
    np.save('fit_yy', y1)
    
    pdb.set_trace()

#np.median(fs)
#0.072624048507333047
#(Pdb) np.percentile(fs, 50-68.2689492137086/2)
#0.067401584039103749
#(Pdb) np.percentile(fs, 50+68.2689492137086/2)
#0.083936234873874749
#(Pdb) np.std(np.sqrt(x1**2+y1**2))*50
#33.636866958753515

def radial_profile(temp):
    x0 = temp.shape[0]/2.-0.5
    y0 = temp.shape[1]/2.-0.5
    yy, xx = np.indices((temp.shape))
    rr = np.sqrt(2)*np.sqrt((xx - x0)**2 + (yy - y0)**2)
    rr = np.round(rr).astype(np.int)

    tbin = np.bincount(rr.ravel(), temp.ravel())
    nr = np.bincount(rr.ravel())
    radialprofile = tbin / nr
    return radialprofile


################################################################
# Do stuff
img = fits.open(img_path)[0].data
psf = fits.open(psf_path)[0].data
volume = np.sum(np.abs(img))
nn = psf.shape[0]/2

# Fit
mask = gen_mask()
guess = np.array([X0, Y0, AA, BB, ALPHA, 20., 0.004, 0.25, 0.])
pars, covar = curve_fit(help_ell, np.arange(0, NX*NY), np.reshape(np.where(mask, img, 0), NX*NY), guess)
res = np.sum(np.abs(img-fold_ell)) / volume
print pars, np.rad2deg(pars[4]), res

#get_err()

# Plots
fig = plt.figure()
plt.imshow(img-last_ell, cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
plt.plot(pars[0], pars[1], 'wo')

fig = plt.figure()
plt.imshow(last_ell, cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
plt.plot(pars[0], pars[1], 'wo')

fig = plt.figure()
plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
plt.plot(pars[0], pars[1], 'wo')



################################################################
# Compute the spread light factor of the ER into the ejecta
help_ell(0, pars[0],pars[1],pars[2],pars[3],pars[4],pars[5],pars[6],pars[7],pars[8])

# Try different hot spots
best_lnl = 1e9
#for ii in range(57572,57573):
for ii in range(4659, 4660): # 12347.828991 -1 4659 0.0700929070407
#for ii in range(9416, 9417):
    np.random.seed(ii)
    tstfld = np.zeros(last_ell.shape)

    # Fiddle with the unfolded model
    # 9416 12348.1347178
#    lim = np.random.rand()
#    test_ell = np.where(np.random.rand(NY, NX) > lim, last_ell, 0)

    test_ell = sts.norm.rvs(last_ell, last_ell/3)
    test_ell = np.where(test_ell < 0, 0, test_ell)
    lim = -1
    
    # 2D convolve, forward fold
    for xi in range(300,340):
        for yi in range(300,340):
            tstfld[yi-nn:yi+nn+1, xi-nn: xi+nn+1] += test_ell[yi, xi]*psf
    bkgimg = img.copy()
    bkgimg[210:430, 210:430] = 0
    bkgflx = np.sum(bkgimg)/(640.**2-220.**2)
    tstfld = tstfld/np.sum(tstfld)*np.sum(img-bkgflx)+bkgflx
    
    lnl = negLogLikelihood(tstfld, img)
    print ii, lnl
    if lnl < best_lnl:
        best_lnl = lnl
        best_lim = lim
        best_see = ii
        best_ell = test_ell
        folded = tstfld

# Best realization
poifld = sts.poisson.rvs(folded)
mk_fits(folded, 'folded.fits')

# Goodness
ng = 1000
goodness = np.zeros(ng)
for ii in range(0, ng):
    poifld = sts.poisson.rvs(folded)
    goodness[ii] = negLogLikelihood(folded, poifld)
    print ii, goodness[ii]
print 'Goodness', np.sum(goodness < best_lnl)

# Prepare psf
mag = 10
alpha = pars[4]
zoopsf = zoom(psf, mag, order=1)
xxp = np.arange(0,mag*psf.shape[1])
yyp = np.arange(0,mag*psf.shape[0])[:,None]
x0p = psf.shape[1]*mag/2
y0p = psf.shape[0]*mag/2
xxp = xxp-x0p
yyp = yyp-y0p
uup = xxp*np.cos(alpha)-yyp*np.sin(alpha)
vvp = xxp*np.sin(alpha)+yyp*np.cos(alpha)

# FWHM
rad_prof = radial_profile(zoopsf)
rmax = 200
fwhmfit = np.sqrt(2)*griddata(rad_prof[1:rmax],np.arange(1,rmax),np.amax(rad_prof[1:])/2,method='linear')
#(Pdb) fwhmfit*5
#690.09775615028445

# Compute the factor of spread light
zoomed = zoom(folded, mag, order=1)
xx = np.arange(0,mag*NX)
yy = np.arange(0,mag*NY)[:,None]
x0=mag*pars[0]
y0=mag*pars[1]
xx = xx-x0
yy = yy-y0
uu = xx*np.cos(alpha)-yy*np.sin(alpha)
vv = xx*np.sin(alpha)+yy*np.cos(alpha)

best_snr = 0
# Get optimal SNR for inner region
for ii in np.linspace(0.95, 1.05, 25):
    aa = ii*mag*smajor
    bb = ii*mag*sminor
    rr = (uup/aa)**2 + (vvp/bb)**2
    psfcts = np.sum(np.where(rr <= 1, zoopsf, 0))
    
    rr = (uu/aa)**2 + (vv/bb)**2
    inner = np.sum(np.where(rr <= 1, zoomed, 0))
    test_snr = psfcts/np.sqrt(inner)
    
    if test_snr > best_snr:
        best_snr = test_snr
        best_rad = ii
        best_cor = np.sum(zoopsf)/psfcts
        print best_snr, best_rad, best_cor, 1/best_cor, 'Best'
    else:
        print test_snr, ii, np.sum(zoopsf)/psfcts

aa = 2032.7452595559723
bb = 1478.7910083532158
rr = (uu/aa)**2 + (vv/bb)**2
outer = np.sum(np.where(rr <= 1, zoomed, 0))
        
aa = best_rad*mag*smajor
bb = best_rad*mag*sminor
rr = (uu/aa)**2 + (vv/bb)**2
inner = np.sum(np.where(rr <= 1, zoomed, 0))

spread_light = inner/(outer-inner)
print best_lnl, best_lim, best_see, spread_light

# Plots
fig = plt.figure()
plt.imshow(folded[290:350, 290:350], cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
fig = plt.figure()
plt.imshow(poifld[290:350, 290:350], cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
fig = plt.figure()
plt.imshow(best_ell[290:350, 290:350], cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
fig = plt.figure()
plt.imshow(last_ell[290:350, 290:350], cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
fig = plt.figure()
plt.imshow(img[290:350, 290:350], cmap='afmhot', origin='lower', interpolation=INTERPOLATION)

fig = plt.figure()
plt.imshow(np.where(mask,folded,0)[290:350, 290:350], cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
fig = plt.figure()
plt.imshow(np.where(mask,img,0)[290:350, 290:350], cmap='afmhot', origin='lower', interpolation=INTERPOLATION)

# MC of fit
fig = plt.figure()
plt.hist(goodness, 100)
plt.axvline(x=best_lnl, color='k', linewidth=2)

########
# Sky map
xlo = 290
ylo = 295
fig = plt.figure(figsize=(2.4375, 2.4375))
ax = plt.Axes(fig, [0., 0., 1., 1.])
# Center
#ax.plot(pars[0]-xlo, pars[1]-ylo, 'wo', mew=0, ms=4)
# Scale
ax.plot([45, 55], [48, 48], 'w', lw=2)
ax.annotate("$0.''5$", xy=(46.9, 49.6), color="w")

ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(img[ylo:350, xlo:350], cmap='afmhot', origin='lower', interpolation=INTERPOLATION, aspect='equal')
ellipse = patches.Ellipse(xy=(pars[0]-xlo, pars[1]-ylo), width=2*smajor, height=2*sminor, angle=-np.rad2deg(alpha), edgecolor='w', fc='None', lw=2)
ax.add_artist(ellipse)
ax.arrow(2.2, 2.2, 0, 8, head_width=1, head_length=1, fc='w', ec='w', linewidth=1.3)
ax.arrow(2.2, 2.2, 8, 0, head_width=1, head_length=1, fc='w', ec='w', linewidth=1.3)
ax.annotate("W", xy=(12.8,1), color="w")
ax.annotate("N", xy=(1,12.8), color="w")
extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
fig.savefig('/Users/silver/box/phd/pro/87a/lim/art/figs/cha_sky.pdf', bbox_inches=extent, pad_inches=0.03)

########
# Mod map
xlo = 290
ylo = 295
fig = plt.figure(figsize=(2.4375, 2.4375))
ax = plt.Axes(fig, [0., 0., 1., 1.])
# Center
#ax.plot(pars[0]-xlo, pars[1]-ylo, 'wo', mew=0, ms=4)
# Scale
ax.plot([45, 55], [48, 48], 'w', lw=2)
ax.annotate("$0.''5$", xy=(46.9, 49.6), color="w")

ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(folded[ylo:350, xlo:350], cmap='afmhot', origin='lower', interpolation=INTERPOLATION, aspect='equal', vmax=np.max(img)/1.6)
ellipse = patches.Ellipse(xy=(pars[0]-xlo, pars[1]-ylo), width=2*smajor, height=2*sminor, angle=-np.rad2deg(alpha), edgecolor='w', fc='None', lw=2)
ax.add_artist(ellipse)
ax.arrow(2.2, 2.2, 0, 8, head_width=1, head_length=1, fc='w', ec='w', linewidth=1.3)
ax.arrow(2.2, 2.2, 8, 0, head_width=1, head_length=1, fc='w', ec='w', linewidth=1.3)
ax.annotate("W", xy=(12.8,1), color="w")
ax.annotate("N", xy=(1,12.8), color="w")
extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
fig.savefig('/Users/silver/box/phd/pro/87a/lim/art/figs/cha_mod.pdf', bbox_inches=extent, pad_inches=0.03)
plt.show()


#[  3.18303954e+02   3.21822521e+02   1.76931736e+01   1.26778055e+01
#   1.29024229e-01   2.00000000e+01   3.41000983e-03   2.70100376e-01
#  -5.80954422e-01] 7.39254376436 0.593753642882
#12353.4160056 0.0472917451798 118380 0.0705455772733
