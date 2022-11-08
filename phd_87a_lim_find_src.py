# Dennis Alp 2016-12-07
# Find the limits in the ejecta of SN 1987A.
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/Dropbox/bin/phd_87a_uplim_find_src.py

import numpy as np
import os
import pdb

from glob import glob
from pyraf import iraf
from astropy.io import fits
from datetime import date
from scipy.interpolate import griddata

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#########
# Parameters
# Logistics
WRK_DIR = "/Users/silver/Dropbox/phd/projects/87a/uplim/src"
os.chdir(WRK_DIR) #Move to designated directory
files = glob('*al2.fits') #Find all files.
# Easier to copy all .fits in folders, therefore only glob psfs and
# strip first thing in loop.

REFERENCE = "w3_r_2015-05-24_drz_al2.fits"
THRESHOLD = 3

# Source location 1-indexed
XCOO = 593.158345778
YCOO = 587.775275464
VKICK = 1.6e6
distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e16

# Cut/subsets
# [554:619,549:634] python 0ind
# [573:599,580:606] python 0ind
BOX = [573,599,580,606]
DAOSUB = "[550:634,555:619]"
DAOSUB = "[581:606,574:599]"

# IRAF parameters
fstars = 'stars.dat'
fwhmpsf = {
"w31102011-01-05_drz_al2.fits": 6.5 ,
"w31602011-01-05_drz_al2.fits": 7.2 ,
"na_h_2010-10-26_map_al2.fits": 4   ,
"na_k_2012-12-14_map_al2.fits": 3.4 ,
"w32252009-12-13_drz_al2.fits": 3.7 ,
"w33362009-12-13_drz_al2.fits": 3.7 ,
"w35022016-06-08_drz_al2.fits": 3.6 ,
"w35552009-12-13_drz_al2.fits": 3.5 ,
"w36452016-06-08_drz_al2.fits": 3.3 ,
"w38142009-12-13_drz_al2.fits": 3.7 ,
"ac_b_2003-01-05_drz_al2.fits": 2.4 ,
"ac_b_2003-08-12_drz_al2.fits": 2.4 ,
"ac_b_2003-11-28_drz_al2.fits": 2.4 ,
"ac_b_2004-12-15_drz_al2.fits": 2.4 ,
"ac_b_2005-04-02_drz_al2.fits": 2.4 ,
"ac_b_2006-12-06_drz_al2.fits": 2.5 ,
"w2_b_2007-05-12_drz_al2.fits": 3.5 ,
"w2_b_2008-02-19_drz_al2.fits": 3.5 ,
"w2_b_2009-04-29_drz_al2.fits": 3.5 ,
"w3_b_2009-12-12_drz_al2.fits": 3.25,
"w3_b_2011-01-05_drz_al2.fits": 3.6 ,
"w3_b_2013-02-06_drz_al2.fits": 3.5 ,
"w3_b_2014-06-15_drz_al2.fits": 3.5 ,
"w3_b_2015-05-24_drz_al2.fits": 3.6 ,
"w3_b_2016-06-08_drz_al2.fits": 3.4 ,
"ac_b_0306-00-00_drz_al2.fits": 2.4 ,
"w2_b_0709-00-00_drz_al2.fits": 3.4 ,
"w2_b_1994-09-24_drz_al2.fits": 3.4 ,
"w3_b_0916-00-00_drz_al2.fits": 3.5 ,
"ac_r_2003-01-05_drz_al2.fits": 2.7 ,
"ac_r_2003-08-12_drz_al2.fits": 2.7 ,
"ac_r_2003-11-28_drz_al2.fits": 2.7 ,
"ac_r_2005-09-26_drz_al2.fits": 2.7 ,
"ac_r_2006-04-15_drz_al2.fits": 2.6 ,
"ac_r_2006-04-29_drz_al2.fits": 2.3 ,
"ac_r_2006-12-06_drz_al2.fits": 2.7 ,
"w2_r_2007-05-12_drz_al2.fits": 3.3 ,
"w2_r_2008-02-19_drz_al2.fits": 3.3 ,
"w2_r_2009-04-29_drz_al2.fits": 3.3 ,
"w3_r_2009-12-12_drz_al2.fits": 3.3 ,
"w3_r_2011-01-05_drz_al2.fits": 3.6 ,
"w3_r_2013-02-06_drz_al2.fits": 3.6 ,
"w3_r_2014-06-15_drz_al2.fits": 3.3 ,
"w3_r_2015-05-24_drz_al2.fits": 3.3 ,
"w3_r_2016-06-08_drz_al2.fits": 3.3 ,
"ac_r_0306-00-00_drz_al2.fits": 2.6 ,
"w2_r_0709-00-00_drz_al2.fits": 3.4 ,
"w2_r_1994-09-24_drz_al2.fits": 3.3 ,
"w3_r_0916-00-00_drz_al2.fits": 3.3 ,
"w322r2009-12-13_drz_al2.fits": 3.4 ,
"w322x2009-12-13_drz_al2.fits": 3.3 ,
"w333r2009-12-13_drz_al2.fits": 3.4 ,
"w333x2009-12-13_drz_al2.fits": 3.4 ,
"w355r2009-12-13_drz_al2.fits": 3.3 ,
"w355x2009-12-13_drz_al2.fits": 3.2 ,
"w36572016-06-08_drz_al2.fits": 3.4 ,
"w381r2009-12-13_drz_al2.fits": 3.4 ,
"w381x2009-12-13_drz_al2.fits": 3.4 ,
"w3_br2015-05-24_drz_al2.fits": 3.6 ,
"w3_bx2015-05-24_drz_al2.fits": 3.5 ,
"w3_rr2015-05-24_drz_al2.fits": 3.3 ,
"w3_rx2015-05-24_drz_al2.fits": 3.2 }

# These considers both sky background and ejecta Poisson
# Computed by find_sig
sigma = {
"ac_b_0306-00-00_drz_al2.fits": 0.00701834489262 ,
"ac_b_2003-01-05_drz_al2.fits": 0.0176363972296 ,
"ac_b_2003-08-12_drz_al2.fits": 0.0227780756507 ,
"ac_b_2003-11-28_drz_al2.fits": 0.0170122807654 ,
"ac_b_2004-12-15_drz_al2.fits": 0.0162289850833 ,
"ac_b_2005-04-02_drz_al2.fits": 0.0186378333763 ,
"ac_b_2006-12-06_drz_al2.fits": 0.0163254725161 ,
"ac_r_0306-00-00_drz_al2.fits": 0.00701116831843 ,
"ac_r_2003-01-05_drz_al2.fits": 0.0303011700681 ,
"ac_r_2003-08-12_drz_al2.fits": 0.0407432131762 ,
"ac_r_2003-11-28_drz_al2.fits": 0.0311758177859 ,
"ac_r_2005-09-26_drz_al2.fits": 0.00842018783725 ,
"ac_r_2006-04-15_drz_al2.fits": 0.0286161094516 ,
"ac_r_2006-04-29_drz_al2.fits": 0.0436876260986 ,
"ac_r_2006-12-06_drz_al2.fits": 0.0296670558917 ,
"na_h_2010-10-26_map_al2.fits": 2.27982611406 ,
"na_k_2012-12-14_map_al2.fits": 2.2356313251 ,
"w2_b_0709-00-00_drz_al2.fits": 0.000553888968355 ,
"w2_b_1994-09-24_drz_al2.fits": 0.00428215404537 ,
"w2_b_2007-05-12_drz_al2.fits": 0.00108902480225 ,
"w2_b_2008-02-19_drz_al2.fits": 0.00111122441586 ,
"w2_b_2009-04-29_drz_al2.fits": 0.00112299110325 ,
"w2_r_0709-00-00_drz_al2.fits": 0.00176941685158 ,
"w2_r_1994-09-24_drz_al2.fits": 0.0127718854892 ,
"w2_r_2007-05-12_drz_al2.fits": 0.00262361709162 ,
"w2_r_2008-02-19_drz_al2.fits": 0.003379339168 ,
"w2_r_2009-04-29_drz_al2.fits": 0.00359374133242 ,
"w31102011-01-05_drz_al2.fits": 0.0235934291968 ,
"w31602011-01-05_drz_al2.fits": 0.0250481916918 ,
"w32252009-12-13_drz_al2.fits": 0.00731436276443 ,
"w322r2009-12-13_drz_al2.fits": 0.00819395071846 ,
"w322x2009-12-13_drz_al2.fits": 0.00814496244355 ,
"w33362009-12-13_drz_al2.fits": 0.00873584764837 ,
"w333r2009-12-13_drz_al2.fits": 0.00931990955467 ,
"w333x2009-12-13_drz_al2.fits": 0.00930611798337 ,
"w35022016-06-08_drz_al2.fits": 0.00397044276575 ,
"w35552009-12-13_drz_al2.fits": 0.0290518148953 ,
"w355r2009-12-13_drz_al2.fits": 0.0287894390672 ,
"w355x2009-12-13_drz_al2.fits": 0.0288194041376 ,
"w36452016-06-08_drz_al2.fits": 0.00402281544854 ,
"w36572016-06-08_drz_al2.fits": 0.00832793530841 ,
"w38142009-12-13_drz_al2.fits": 0.0239838675949 ,
"w381r2009-12-13_drz_al2.fits": 0.0249604464494 ,
"w381x2009-12-13_drz_al2.fits": 0.0249613407583 ,
"w3_b_0916-00-00_drz_al2.fits": 0.00433296341128 ,
"w3_b_2009-12-12_drz_al2.fits": 0.0125557885671 ,
"w3_b_2011-01-05_drz_al2.fits": 0.0108696330795 ,
"w3_b_2013-02-06_drz_al2.fits": 0.00970687833124 ,
"w3_b_2014-06-15_drz_al2.fits": 0.00924723033467 ,
"w3_b_2015-05-24_drz_al2.fits": 0.0094900377064 ,
"w3_b_2016-06-08_drz_al2.fits": 0.0154812646506 ,
"w3_br2015-05-24_drz_al2.fits": 0.00974715488903 ,
"w3_bx2015-05-24_drz_al2.fits": 0.00975035655871 ,
"w3_r_0916-00-00_drz_al2.fits": 0.00691686632371 ,
"w3_r_2009-12-12_drz_al2.fits": 0.0122142571276 ,
"w3_r_2011-01-05_drz_al2.fits": 0.0204315629735 ,
"w3_r_2013-02-06_drz_al2.fits": 0.018637654414 ,
"w3_r_2014-06-15_drz_al2.fits": 0.0179594227444 ,
"w3_r_2015-05-24_drz_al2.fits": 0.0177253756152 ,
"w3_r_2016-06-08_drz_al2.fits": 0.0273143813317 ,
"w3_rr2015-05-24_drz_al2.fits": 0.0180942252429 ,
"w3_rx2015-05-24_drz_al2.fits": 0.0181124092083 }

# Definitions
detections = []
unscdet = [] # Unscaled detections
labels = []
    
#########
# Help functions
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

# Days for scaling
SNDATE = date(1987, 2, 23)
def get_days(fname): # 1994-09-24_drz_al2.fits
#    fname = fname.replace("0306-00-00","2006-12-06") # Set the folded to latest time
#    fname = fname.replace("0709-00-00","2009-04-29")
#    fname = fname.replace("0916-00-00","2016-06-08")
    fname = fname.replace("0306-00-00","2004-12-06") # Set the folded to latest time
    fname = fname.replace("0709-00-00","2008-04-29")
    fname = fname.replace("0916-00-00","2013-01-01")

    
    yr = int(fname[-23:-19])
    month = int(fname[-18:-16].lstrip("0"))
    day = int(fname[-15:-13].lstrip("0"))
    print "Date from filename:", yr, month, day
    return (date(yr, month, day)-SNDATE).days*24*60*60
REFT = get_days(REFERENCE)

# Insert psf and try to find it
for img in files[0:]:
    print img
    # Prepare the psf
    # Find the FWHM by fitting to the see result from IRAF
    try:
        psf = fits.open(img+".see.2.fits")[0].data
        subnx = psf.shape[0]
        subny = psf.shape[1]
        rad_prof = radial_profile(psf)
        rmax = np.int(np.ceil(fwhmpsf[img])+3)
        fwhmfit = np.sqrt(2)*griddata(rad_prof[0:rmax],np.arange(0,rmax),np.amax(rad_prof)/2,method='linear')
        print "Check match iraf, fit:",fwhmpsf[img],fwhmfit
    except Exception:
        print "PSF not found, using ", fwhmpsf[img]
        fwhmfit = fwhmpsf[img]

    # Get kick distance
    telapsed = get_days(img)
    offset = np.rad2deg(np.arctan(VKICK*telapsed/(distance*1e3*pc)))*60*60*1000/25
    print "Offset from center in pixels:", offset
    
    # Find psf
    iraf.daofind(img + DAOSUB,output=fstars,fwhmpsf=fwhmfit,sigma=sigma[img],\
                     threshold=THRESHOLD,datamin=-0.1,datamax=9999,verify=False)
    stscoo=iraf.pdump(fstars,'XCENTER,YCENTER','yes',Stdout=1)

    # Save coordinates of detections
    sources = 0
    for xy in stscoo:
        testx, testy = xy.split()
        testx = float(testx)
        testy = float(testy)
        rkick = np.sqrt((testx+BOX[2]-XCOO)**2+(testy+BOX[0]-YCOO)**2)
#        print rkick
        if rkick < offset:
            sources += 1
            scaledx = (BOX[2]+testx-XCOO)*REFT/telapsed+XCOO-BOX[2]
            scaledy = (BOX[0]+testy-YCOO)*REFT/telapsed+YCOO-BOX[0]
            detections.append([scaledx, scaledy])
            unscdet.append([testx, testy])
            print testx, testy
            labels.append(img[-26:-13])
    print sources,"\n"

# Plot detections
sky_img = fits.open(REFERENCE)[0].data
sky_img = sky_img[BOX[0]:BOX[1],BOX[2]:BOX[3]]
detections = np.array(detections).astype(np.float)
unscdet = np.array(unscdet).astype(np.float)

offset = np.rad2deg(np.arctan(VKICK*REFT/(distance*1e3*pc)))*60*60*1000/25

plt.imshow(sky_img,interpolation='nearest', cmap='afmhot', origin='lower')
theta = np.arange(0,2*np.pi,0.001)
plt.plot(detections[:,0]-1,detections[:,1]-1,'og')
plt.plot(unscdet[:,0]-1,unscdet[:,1]-1,'.b')
plt.plot(XCOO-1-BOX[2],YCOO-1-BOX[0],'wo')
plt.plot(XCOO-1-BOX[2]+offset*np.cos(theta),YCOO-1-BOX[0]+offset*np.sin(theta),'k')

lab_count = 0
for i in unscdet:
    plt.annotate(labels[lab_count], xy=i-1) #, textcoords='data')
    lab_count += 1
plt.show()
