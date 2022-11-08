# Dennis Alp 2016-12-01
# Find the count (poisson) noise in drizzled HST images.
# This is specifically for the ejecta in SN 1987A. 
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/Dropbox/bin/phd_87a_uplim_find_sig.py

import numpy as np
import os
import pdb

from glob import glob
from pyraf import iraf
from astropy.io import fits
from datetime import date

#########
# Parameters
# Pixel sizes for HST
nir = 0.1282500028610229
uvis = 0.03962000086903572
wfpc2 = 0.046
PIXSCALE = 0.025

# Some data for SN 1987A
distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e16
cc = 299792458
kick = 1600e3
SNDATE = date(1987, 2, 23)

def get_days(fname): # 1994-09-24_drz_al2.fits
    fname = fname.replace("0306-00-00","2006-12-06") # Set the folded to latest time
    fname = fname.replace("0709-00-00","2009-04-29")
    fname = fname.replace("0916-00-00","2016-06-08")
    
    yr = int(fname[-23:-19])
    month = int(fname[-18:-16].lstrip("0"))
    day = int(fname[-15:-13].lstrip("0"))
    print "Date from filename:", yr, month, day
    return (date(yr, month, day)-SNDATE).days*24*60*60

# Logistics
WRK_DIR = "/Users/silver/Dropbox/phd/projects/87a/uplim/sig"
os.chdir(WRK_DIR) #Move to designated directory
files = glob('*.fits') #Find all files.

# Image data 0-index
XCOO=592.158345778
YCOO=586.775275464
NX = 1150
NY = 1115

# Sky sigma for all observations
sky_sig = {
"w31102011-01-05_drz_al2.fits": 0.016   ,
"w31602011-01-05_drz_al2.fits": 0.024   ,
"na_h_2010-10-26_map_al2.fits": 2.2     ,
"na_k_2012-12-14_map_al2.fits": 2.2     ,
"w32252009-12-13_drz_al2.fits": 0.0047  ,
"w33362009-12-13_drz_al2.fits": 0.0045  ,
"w35022016-06-08_drz_al2.fits": 0.0033  ,
"w35552009-12-13_drz_al2.fits": 0.012   ,
"w36452016-06-08_drz_al2.fits": 0.0022  ,
"w38142009-12-13_drz_al2.fits": 0.01    ,
"ac_b_2003-01-05_drz_al2.fits": 0.01    ,
"ac_b_2003-08-12_drz_al2.fits": 0.01375 ,
"ac_b_2003-11-28_drz_al2.fits": 0.011   ,
"ac_b_2004-12-15_drz_al2.fits": 0.01    ,
"ac_b_2005-04-02_drz_al2.fits": 0.01    ,
"ac_b_2006-12-06_drz_al2.fits": 0.01    ,
"w2_b_2007-05-12_drz_al2.fits": 0.00075 ,
"w2_b_2008-02-19_drz_al2.fits": 0.00067 ,
"w2_b_2009-04-29_drz_al2.fits": 0.00058 ,
"w3_b_2009-12-12_drz_al2.fits": 0.005   ,
"w3_b_2011-01-05_drz_al2.fits": 0.0065  ,
"w3_b_2013-02-06_drz_al2.fits": 0.00375 ,
"w3_b_2014-06-15_drz_al2.fits": 0.003375,
"w3_b_2015-05-24_drz_al2.fits": 0.0045  ,
"w3_b_2016-06-08_drz_al2.fits": 0.00975 ,
"ac_b_0306-00-00_drz_al2.fits": 0.004   ,
"w2_b_0709-00-00_drz_al2.fits": 0.000233,
"w2_b_1994-09-24_drz_al2.fits": 0.0006  ,
"w3_b_0916-00-00_drz_al2.fits": 0.002   ,
"ac_r_2003-01-05_drz_al2.fits": 0.015   ,
"ac_r_2003-08-12_drz_al2.fits": 0.02125 ,
"ac_r_2003-11-28_drz_al2.fits": 0.015   ,
"ac_r_2005-09-26_drz_al2.fits": 0.004   ,
"ac_r_2006-04-15_drz_al2.fits": 0.0145  ,
"ac_r_2006-04-29_drz_al2.fits": 0.03    ,
"ac_r_2006-12-06_drz_al2.fits": 0.016   ,
"w2_r_2007-05-12_drz_al2.fits": 0.00062 ,
"w2_r_2008-02-19_drz_al2.fits": 0.0007  ,
"w2_r_2009-04-29_drz_al2.fits": 0.0007  ,
"w3_r_2009-12-12_drz_al2.fits": 0.005   ,
"w3_r_2011-01-05_drz_al2.fits": 0.00625 ,
"w3_r_2013-02-06_drz_al2.fits": 0.006875,
"w3_r_2014-06-15_drz_al2.fits": 0.0054  ,
"w3_r_2015-05-24_drz_al2.fits": 0.0055  ,
"w3_r_2016-06-08_drz_al2.fits": 0.013   ,
"ac_r_0306-00-00_drz_al2.fits": 0.0033  ,
"w2_r_0709-00-00_drz_al2.fits": 0.00033 ,
"w2_r_1994-09-24_drz_al2.fits": 0.000875,
"w3_r_0916-00-00_drz_al2.fits": 0.002   ,
"w322r2009-12-13_drz_al2.fits": 0.0057  ,
"w322x2009-12-13_drz_al2.fits": 0.0057  ,
"w333r2009-12-13_drz_al2.fits": 0.0055  ,
"w333x2009-12-13_drz_al2.fits": 0.0055  ,
"w355x2009-12-13_drz_al2.fits": 0.012   ,
"w355r2009-12-13_drz_al2.fits": 0.012   ,
"w36572016-06-08_drz_al2.fits": 0.0029  ,
"w381x2009-12-13_drz_al2.fits": 0.012   ,
"w381r2009-12-13_drz_al2.fits": 0.012   ,
"w3_bx2015-05-24_drz_al2.fits": 0.005   ,
"w3_br2015-05-24_drz_al2.fits": 0.005   ,
"w3_rr2015-05-24_drz_al2.fits": 0.0064  ,
"w3_rx2015-05-24_drz_al2.fits": 0.0064  }


# Misc
INTERPOLATION = 'nearest'
save_sig = []

#########
# Help functions
def mask_img(img, offset):
    y,x = np.ogrid[-YCOO:NY-YCOO, -XCOO:NX-XCOO]
    mask = x*x + y*y <= offset**2
    return img[mask]

#########
# Loop files and find noise level
for img_path in files:
    print "\n",img_path
    img = fits.open(img_path)[0].data
    exposure = fits.getval(img_path,'EXPTIME')

    if 'ac_' in img_path:
        imgscale = np.sqrt(0.028*0.025)
        gain = 2.
    elif 'na_' in img_path:
        imgscale = fits.getval(img_path,'HIERARCH ESO INS PIXSCALE')
        gain = 12.1
        exposure = exposure*18 if 'na_h' in img_path else exposure*21
    elif 'w2_' in img_path:
        imgscale = wfpc2
        gain = 7.
    elif 'w31' in img_path: # IR
        imgscale = fits.getval(img_path,'D001ISCL')
        gain = 2.5
    else:
        imgscale = uvis # UVIS
        gain = 1.5
    print "Img scale and gain:", imgscale, gain

    try:
        print "Date from .fits:   ", fits.getval(img_path,'DATE-OBS')
    except KeyError:
        print "WARNING: DATE-OBS key not found."
    
    telapsed = get_days(img_path)
    offset = np.rad2deg(np.arctan(kick*telapsed/(distance*1e3*pc)))*60*60*1000/25
    print "Offset from center in pixels:", offset

    region = mask_img(img, offset)
    ovrsamp = imgscale**2/PIXSCALE**2
    if ovrsamp < 1:
        ovrsamp = 1
    region = np.amax(region)*exposure*ovrsamp/gain
    ej_sig = np.sqrt(region)
    ej_sig = ej_sig/exposure/ovrsamp*gain
    sigma = np.sqrt(ej_sig**2+sky_sig[img_path[-28:]]**2)
    print "Most electrons:", region, region/exposure/ovrsamp
    print "Ejecta sigma:", ej_sig*exposure*ovrsamp, ej_sig
    print "Exposure, oversampling:", exposure, ovrsamp
    print "Ejecta tot:", sigma, ej_sig, sky_sig[img_path[-28:]], ej_sig/sky_sig[img_path[-28:]]
    save_sig.append(sigma)

#    plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(img)/1000, vmax=np.amax(img)))
#    plt.imshow(img)
#    plt.show()

print ""
counter=0
for ff in files:
    print ff,":",save_sig[counter],","
    counter+=1
