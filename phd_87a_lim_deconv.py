# Dennis Alp 2017-06-02
# Fold an image through a response and try to deconvolve it
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/Dropbox/bin/phd_87a_uplim_deconv.py

import numpy as np
import os
import pdb
import time
from glob import glob
from datetime import date

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.ndimage.interpolation import zoom
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import shift
from pyraf import iraf
from astropy.io import fits



################################################################
# Parameters
# Logistics
WRK_DIR = "/Users/silver/Dropbox/phd/projects/87a/uplim/deconv"
os.chdir(WRK_DIR) #Move to designated directory
rm = glob('*.fits')
for ii in rm:
   os.remove(ii)

# Parameters
NX = 200
NY = 200
GRANULARITY = 0.99
RSP_SIG = 20.
TOT_CTS = 1e9
ZOOM = 1/20.
PIX_WTH = 5 # Linear size of a CCD pixel
PIX_FRC = 5 # Size of the final pixel size

# IRAF parameters
FAKE = "FAKE.fits" # Name of tampered image, handle with care!
COPY = "COPY.fits" # Make a new copy of original .fits but with different data type since daofind has issues
PSF = "PSF.fits"
FSTARS = 'stars.dat'



################################################################
# Help functions
# Just reasonable plot of sky images, good for debugging
def sky_plt(image):
    plt.figure()
    plt.imshow(image,interpolation='nearest', cmap='afmhot', origin='lower')

# Utility for making fits image out of a numpy array
def mk_fits(image, output):
    hdu = fits.PrimaryHDU(image)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(output, clobber=True)
    hdulist.close()
    
# For FWHM determination from psf see
def radial_profile(data):
    x0 = subnx/2.-0.5
    y0 = subny/2.-0.5
    yy, xx = np.indices((data.shape))
    rr = np.sqrt(2)*np.sqrt((xx - x0)**2 + (yy - y0)**2)
    rr = np.round(rr).astype(np.int)

    tbin = np.bincount(rr.ravel(), data.ravel())
    nr = np.bincount(rr.ravel())
    radialprofile = tbin / nr
    return radialprofile

# Make the inner ring of SN 1987A
def mk_ring(aa, bb, alpha, sig, tilt, phi):
    xx = np.arange(0, NX)
    yy = np.arange(0, NY)[:,None]
    xx = xx-NX/2
    yy = yy-NY/2
    uu = xx*np.cos(alpha)-yy*np.sin(alpha)
    vv = xx*np.sin(alpha)+yy*np.cos(alpha)
    rr = (uu/aa)**2 + (vv/bb)**2

    res = np.exp(-np.abs(np.sqrt(rr)-1)**2/abs(sig))*(1+tilt*np.cos(np.arctan2(vv,uu)+phi))
    res /= np.amax(res)
    screen = np.random.random((NY, NX))
#    res = np.where(screen < GRANULARITY, 0, res)
    res = np.where(np.exp(-np.abs(np.sqrt(rr)-1)**2/abs(sig)) > 0.99975, res, 0)
    return res

# Fold an image through a response
def fold_map(map_true):
    def mk_rsp(x0, y0):
        xx = np.arange(0, NX)
        yy = np.arange(0, NY)[:,None]
        xx = xx-x0
        yy = yy-y0
        rr = xx**2+yy**2
        return np.exp(-rr/RSP_SIG**2)
    map_rsp = np.zeros(map_true.shape)
    for xx in range(0,NX):
        print xx
        for yy in range(0,NY):
            rsp = mk_rsp(xx, yy)
            map_rsp[yy, xx] = np.sum(rsp*map_true)
    return map_rsp

# Transform the folded image into as it would be observed
def obs_map(map_rsp):
    map_obs = zoom(map_rsp, ZOOM, order=0, mode='constant', cval=0.0, prefilter=False)
    map_obs = TOT_CTS/np.sum(map_obs)*map_obs
    map_obs = np.random.poisson(map_obs)
    map_obs = zoom(map_obs, 1/ZOOM, order=0, mode='constant', cval=0.0, prefilter=False)
    return map_obs.astype('f')

# Transform the folded image into as it would be observed
def drz_map(map_rsp):
    map_drz = np.zeros((NY/PIX_FRC, NX/PIX_FRC))
    for ix, xx in enumerate(range(PIX_FRC/2, NX, PIX_FRC)):
        for iy, yy in enumerate(range(PIX_FRC/2, NY, PIX_FRC)):
            xlo = np.max((xx-PIX_WTH/2, 0))
            xhi = xx+PIX_WTH/2
            ylo = np.max((yy-PIX_WTH/2, 0))
            yhi = yy+PIX_WTH/2
            map_drz[iy, ix] = np.mean(map_rsp[ylo:yhi, xlo:xhi])
#            pdb.set_trace()
                                                  
    map_drz = TOT_CTS/np.sum(map_drz)*map_drz
    map_drz = np.random.poisson(map_drz)
    return map_drz.astype('f')

# Deconvolve the image using Richardson-Lucy
def deconv_map(sky_map, output):
    def mk_psf():
        nn = np.int(6*RSP_SIG*ZOOM+1)
        xx = np.arange(0, nn)
        yy = np.arange(0, nn)[:,None]
        xx = xx-nn/2
        yy = yy-nn/2
        rr = xx**2+yy**2
        psf = np.exp(-rr/(RSP_SIG*ZOOM)**2)
        mk_fits(psf, PSF)

    def mk_psf():
        xx = np.arange(0, NX)
        yy = np.arange(0, NY)[:,None]
        xx = xx-NX/2
        yy = yy-NY/2
        rr = xx**2+yy**2
        psf = np.exp(-rr/RSP_SIG**2)
        mk_fits(psf, PSF)

    def mk_psf():
        nn = np.int(6*RSP_SIG/PIX_FRC+1)
        xx = np.arange(0, nn)
        yy = np.arange(0, nn)[:,None]
        xx = xx-nn/2
        yy = yy-nn/2
        rr = xx**2+yy**2
        psf = np.exp(-rr/(RSP_SIG/PIX_FRC)**2)
        mk_fits(psf, PSF)

    def mk_psf():
        xx = np.arange(0, NX)
        yy = np.arange(0, NY)[:,None]
        xx = xx-NX/2
        yy = yy-NY/2
        rr = xx**2+yy**2
        psf = np.exp(-rr/(RSP_SIG)**2)
        
        psf_drz = np.zeros((NY/PIX_FRC, NX/PIX_FRC))
        for ix, xx in enumerate(range(PIX_FRC/2, NX, PIX_FRC)):
            for iy, yy in enumerate(range(PIX_FRC/2, NY, PIX_FRC)):
                xlo = np.max((xx-PIX_WTH/2, 0))
                xhi = xx+PIX_WTH/2
                ylo = np.max((yy-PIX_WTH/2, 0))
                yhi = yy+PIX_WTH/2
                psf_drz[iy, ix] = np.mean(psf[ylo:yhi, xlo:xhi])
        mk_fits(psf_drz, PSF)
        
    mk_psf()
    iraf.lucy(input=sky_map, psf=PSF, output=output, adu=1, noise=0, niter=24, accel_method='none', nsave=8, limchisq=0.0001)
    return fits.open('deconv_copy.fits')[0].data


    
################################################################
iraf.stsdas()
iraf.analysis()
iraf.restore()

#3.12575493e+01   2.28780350e+01   8.79246740e-02  2.27276768e-02   2.41630537e-01  -1.35412431e+00
map_true = mk_ring(3.1e+01, 2.3e+01, 8.8e-02, 2.3e-02, 2.4e-01, -1.4e+00)
#map_true[NY/2, NX/2] = 2.
map_rsp = fold_map(map_true)
#map_obs = obs_map(map_rsp)
map_drz = drz_map(map_rsp)
mk_fits(map_drz, COPY)
map_deconv = deconv_map(COPY, 'deconv_copy.fits')

sky_plt(map_true)
sky_plt(map_rsp)
sky_plt(map_drz)
sky_plt(map_deconv)
plt.show()
pdb.set_trace()

########
psf = fits.open(img_path+".see.2.fits")[0].data
subny = psf.shape[0]
subnx = psf.shape[1]    
rad_prof = radial_profile(psf)
rmax = np.int(np.ceil(fwhmpsf[img_path])+3)
fwhmfit = np.sqrt(2)*griddata(rad_prof[0:rmax],np.arange(0,rmax),np.amax(rad_prof)/2,method='linear')
print "Check match iraf, fit:",fwhmpsf[key],fwhmfit
sigmaobs = sigma[key]

subny = int(subny/2)
subnx = int(subnx/2)
psf = shift(psf, SHIFTS, order=1, mode='constant', cval=0.0, prefilter=False)
psf = psf[subny-APERTURE:subny+APERTURE+1,subnx-APERTURE:subnx+APERTURE+1]

for ii in range(BOX[0],BOX[1]):
    for jj in range(BOX[2],BOX[3]):
        # Do the binary search for a limit at a point
        unfluxed = bin_search(COPY, psf, np.floor(ii), np.floor(jj))
        limits[ii-BOX[0],jj-BOX[2],img_nr] = unfluxed
        print "Limit:", key, jj, ii, limits[ii-BOX[0],jj-BOX[2],img_nr]

img_nr += 1
print ""
