# Dennis Alp 2016-12-07
# Find the limits in the ejecta of SN 1987A.
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/box/bin/phd_87a_lim_find_lim.py

import numpy as np
import os
import pdb
import time
from glob import glob
from datetime import date

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import rotate
from scipy.ndimage.interpolation import zoom
from pyraf import iraf
from astropy.io import fits

#########
# Parameters
# Logistics
WRK_DIR = "/Users/silver/dat/alm/87a/lim"
os.chdir(WRK_DIR) #Move to designated directory
files = sorted(glob('alma*.fits')) #Find all files.
BIN_SEARCH_STOP = 0.01
RTRSH = np.sqrt(2) # Threshold radius for detection when inserting artificial sources
RTRSH = 6
SHRINK = 1 # Float, or 1 if undesired
OUTPUT = "limits_alma"

NY = 40
NX = 40
PSFY = 50
PSFX = 50

XCOO = 83.86661458333334
YCOO = -69.2635813888889
YCOO = -69.26975194444444
VKICK = 0.4e6
distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e16

# IRAF parameters
FAKE = "FAKE.fits" # Name of tampered image, handle with care!
COPY = "COPY.fits" # Make a new copy of original .fits but with different data type since daofind has issues
FSTARS = 'stars.dat'

sigma = {
# From find_thld
"alma_2014-09-02_212_213.fits": 4.22545340531e-05,
"alma_2014-09-02_212_247.fits": 3.16911325604e-05,
"alma_2014-09-02_233_233.fits": 5.92417818956e-05,
"alma_2014-09-02_247_247.fits": 4.43683329131e-05}

## From find_thld
#sigma = {
#"alma_2014-09-02_212_247.fits": 3.09466631569e-05,
#"alma_2014-09-02_233_233.fits": 7.83525388361e-05,
#"alma_2014-09-02_247_247.fits": 5.04554094034e-05}

# Heuristic
#"alma_2014-09-02_212_213.fits": 5e-5 }
# 3 samples in original image: 1.7262078e-05 1.5970476e-05 1.6458962e-05
# Based on 3 samples (fwhm~1): 1.0509892964625505e-05 9.194606625325014e-06 1.0962775941903835e-05
#"alma_2014-09-02_212_213.fits": 1.02e-5 }
# Based on 3 samples (fwhm~2): 7.6486692149898068e-06 9.4231245812237393e-06 7.1015945260976804e-06
#"alma_2014-09-02_212_213.fits": 8.06e-6 }
#"alma_2014-09-02_212_213.fits": 1.6e-5 }
    
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

# Makes beam image of data from .fits header
def mk_beam(image):
    scale = fits.getval(image,'CDELT2')
    minor = fits.getval(image,'BMIN')/scale
    major = fits.getval(image,'BMAJ')/scale
    angle = np.deg2rad(fits.getval(image,'BPA')+90)

    xx = np.arange(0,PSFX)
    yy = np.arange(0,PSFY)[np.newaxis].T
    x0 = PSFX/2.-0.5
    y0 = PSFY/2.-0.5
    fwhm2sig = 1/(2*np.sqrt(2*(np.log(2))))
    sx = major*fwhm2sig
    sy = minor*fwhm2sig

    uu = (xx-x0)*np.cos(angle)+(yy-y0)*np.sin(angle)
    vv = (x0-xx)*np.sin(angle)+(yy-y0)*np.cos(angle)
    psf = np.exp(-0.5*(uu**2/sx**2+vv**2/sy**2))
    fwhm = (major+minor)/2
    return psf, fwhm

# convert Jy beam^-1 to Jy
def beam_corr(image, flux):
    scale = fits.getval(image,'CDELT2')
    minor = fits.getval(image,'BMIN')/scale
    major = fits.getval(image,'BMAJ')/scale
    fwhm2sig = 1/(2*np.sqrt(2*(np.log(2))))
    sx = major*fwhm2sig
    sy = minor*fwhm2sig
    return flux/(2*np.pi*sx*sy)

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

# transform image such that beams are circular, suited for daofind
def pre_trans(img):
    minor = fits.getval(key,'BMIN')
    major = fits.getval(key,'BMAJ')
    angle = fits.getval(key,'BPA')
    
    img = rotate(img, angle, axes=(1, 0), reshape=False, order=1, mode='constant', cval=0.0, prefilter=False)
    img = zoom(img, (1, major/minor), order=1, mode='constant', cval=0.0, prefilter=False)
    img = zoom(img, 1./SHRINK, order=1, mode='constant', cval=0.0, prefilter=False)
#    sky_plt(img)
#    pdb.set_trace()
    return img

# transform daofind coords back into obs frame
def post_trans(testx, testy):
    testx = float(testx)
    testy = float(testy)
    testx -= 1
    testy -= 1

    testx *= SHRINK
    testy *= SHRINK
    
    minor = fits.getval(key,'BMIN')
    major = fits.getval(key,'BMAJ')
    angle = -np.deg2rad(fits.getval(key,'BPA'))
    testx = testx * minor/major

    nimg = fits.open(key)[0].data[0,0].shape
    nx = nimg[0]/2.
    ny = nimg[1]/2.

    uu = (testx-nx)*np.cos(angle)+(testy-ny)*np.sin(angle)
    vv = (nx-testx)*np.sin(angle)+(testy-ny)*np.cos(angle)
    return uu+nx, vv+ny

# adding but removing counts so that to make it just fit, i.e. max(ejecta,artificial)
def hide_star(imgcp, psfcp, flux, xstart, xstop, ystart, ystop):
    norm = np.sum(psfcp)
    psfcp = flux/norm * psfcp
    test_psf = np.zeros(imgcp.shape)
    test_psf[ystart:ystop, xstart:xstop] = psfcp
    if np.all(test_psf[(ystart+ystop)/2, (xstart+xstop)/2] < imgcp[(ystart+ystop)/2, (xstart+xstop)/2]):
        imgcp = pre_trans(imgcp)
        mk_fits(imgcp, FAKE)
        return False
    
    test_img = np.where(test_psf > imgcp, test_psf, imgcp)
    test_img = pre_trans(test_img)
    mk_fits(test_img, FAKE)
    return True

# tries to find hidden star using daofind
def try_find(xx, yy):
    # Find psf
    """
    BUG: In rare circumstances daofind may abort with a "pixel file truncation
        error" when trying to read back the convolved images it has just 
    	written. This only occurs on certain sized images and is due to
    	the interaction of the read, write and boundary extension in image
    	i/o. For example daofind works fine on a 640 by 1024 image but fails on
    	one that is 641 by 1025 pixels.
    
    STATUS:	The problem is fixed in 2.9. The solution was to add an imflush
    	statement to flush the image i/o buffers after the write was
    	complete and before initiating the read. There is no workaround.
    	Contact the IRAF group for a fix.
    
    https://github.com/joequant/iraf/blob/master/local/bugs.log
    """
    try:
        iraf.daofind(FAKE,output=FSTARS,fwhmpsf=fwhm/SHRINK,sigma=sigmaobs,threshold=3,datamin=-0.1,datamax=9999,verify=False, Stdout=1)
    except Exception:
        print key, 'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG'
        return True

    stscoo=iraf.pdump(FSTARS,'XCENTER,YCENTER','yes',Stdout=1)
    
    # Check if found
#    print 'TCE', len(stscoo)
    for xy in stscoo:
        testx, testy = xy.split()
        testx, testy = post_trans(testx, testy)
        delx = testx-xx
        dely = testy-yy
        if np.sqrt(delx**2+dely**2) < fwhm/2.:
            return True

    return False
    
# this uses binary search to find a limit by hiding and daofinding
def bin_search(img_path, psf, yy, xx, key):
    xstart = xx-PSFX/2
    xstop  = xx+PSFX/2
    ystart = yy-PSFY/2
    ystop  = yy+PSFY/2
    img = fits.open(img_path)[0].data[0,0]
    
    upper = np.sum(img[ystart:ystop, xstart:xstop])/8 # 8 is heuristic
    
    up = 2*upper
    lo = 1e-12 # Suppress divide by 0 warning
    oflux = -1
    flux = 1e-12
    
    # Catch detections of stuff in ejecta
    hide_star(img.copy(), psf.copy(), flux, xstart, xstop, ystart, ystop)
    if try_find(xx, yy):
        print "WARNING: DETECTION"
        while abs(flux/oflux-1) > BIN_SEARCH_STOP:
            if hide_star(img.copy(), psf.copy(), flux, xstart, xstop, ystart, ystop):
                up = flux
            else:
                lo = flux
            oflux = flux
            flux = (up + lo)/2.
    else:
        while abs(flux/oflux-1) > BIN_SEARCH_STOP:
#            print flux, oflux, flux/oflux
#            sky_plt(fits.open(FAKE)[0].data)
            hide_star(img.copy(), psf.copy(), flux, xstart, xstop, ystart, ystop)
            
            if try_find(xx, yy):
                up = flux
            else:
                lo = flux
            oflux = flux
            flux = (up + lo)/2.

    return beam_corr(key, flux)


#########
# Definitions
limits = np.zeros((NY, NX, len(files)))
img_nr = 0
iraf.stsdas()
iraf.hst_calib()

#########
# Insert psf and try to find it
for img_path in files[:]:
    key = img_path[-28:]
    print key
    # To cast data type, this is SUPER IMPORTANT. Daofind throws fits otherwise
#    mk_fits(fits.open(img_path)[0].data.astype('>f8'), COPY)
    mk_fits(fits.open(img_path)[0].data.astype(np.float64), COPY)

    # Prepare the psf
    # Find the FWHM by fitting to the see result from IRAF
    psf, fwhm = mk_beam(key)
    sigmaobs = sigma[key]

    iraf.daofind(COPY,output=FSTARS,fwhmpsf=fwhm/SHRINK,sigma=sigmaobs,threshold=3,datamin=-0.1,datamax=9999,verify=False, Stdout=1)
        
    for ii in range(0, NY):
        for jj in range(0, NX): # TODO FROM HERE
            # Do the binary search for a limit at a point
            yy = YCOO-(NY/2.-ii)*0.0125/3600
            xx = XCOO+(NX/2.-jj)*0.0125/3600/np.cos(np.deg2rad(YCOO))
            xx, yy = coords2pix(key, xx, yy)
            xx = np.int(np.round(xx))
            yy = np.int(np.round(yy))
            
            unfluxed = bin_search(COPY, psf, yy, xx, key)
            limits[ii, jj, img_nr] = unfluxed
            print "Limit:", key, jj, ii, limits[ii, jj, img_nr]

    img_nr += 1
    np.save(OUTPUT + '_' + key + "new.npy", limits)
    print ""

np.save(OUTPUT + "_new.npy", limits)
pdb.set_trace()
# Asymmetric beam
#Limit: alma_2014-09-02_212_213.fits 0 0 0.00447974503114
#Limit: alma_2014-09-02_212_213.fits 1 0 0.00428657071673
#Limit: alma_2014-09-02_212_213.fits 0 1 0.00472904913374
#Limit: alma_2014-09-02_212_213.fits 1 1 0.00473177796754

# Transform to circular beam
#Limit: alma_2014-09-02_212_213.fits 0 0 0.00427765127039
#Limit: alma_2014-09-02_212_213.fits 1 0 0.00428657071673
#Limit: alma_2014-09-02_212_213.fits 0 1 0.00431482585204
#Limit: alma_2014-09-02_212_213.fits 1 1 0.00432327915028
