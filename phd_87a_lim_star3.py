# Dennis Alp 2018-01-23
# Just get a light curve of star 3.
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/box/bin/phd_87a_uplim_star3.py

import numpy as np
import os
import pdb
from glob import glob
from datetime import date

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.integrate import simps
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import shift
from pyraf import iraf
from astropy.io import fits

from phd_87a_red import red_mod, cor_red

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

#########
# Parameters
# Logistics
WRK_DIR = "/Users/silver/dat/hst/87a/all/"
BOX = [577,597,582,602]

os.chdir(WRK_DIR) #Move to designated directory
files = sorted(glob('*.fits'))

bandpar_arg = {
"w3_b_2009-12-12_drz_al2.fits": 'wfc3,uvis2,f438w,aper#0.2'       ,
"w3_b_2011-01-05_drz_al2.fits": 'wfc3,uvis2,f438w,aper#0.2'       ,
"w3_b_2013-02-06_drz_al2.fits": 'wfc3,uvis2,f438w,aper#0.2'       ,
"w3_b_2014-06-15_drz_al2.fits": 'wfc3,uvis2,f438w,aper#0.2'       ,
"w3_b_2015-05-24_drz_al2.fits": 'wfc3,uvis2,f438w,aper#0.2'       ,
"w3_b_2016-06-08_drz_al2.fits": 'wfc3,uvis2,f438w,aper#0.2'       ,
"w3_r_2009-12-12_drz_al2.fits": 'wfc3,uvis2,f625w,aper#0.2'       ,
"w3_r_2011-01-05_drz_al2.fits": 'wfc3,uvis2,f625w,aper#0.2'       ,
"w3_r_2013-02-06_drz_al2.fits": 'wfc3,uvis2,f625w,aper#0.2'       ,
"w3_r_2014-06-15_drz_al2.fits": 'wfc3,uvis2,f625w,aper#0.2'       ,
"w3_r_2015-05-24_drz_al2.fits": 'wfc3,uvis2,f625w,aper#0.2'       ,
"w3_r_2016-06-08_drz_al2.fits": 'wfc3,uvis2,f625w,aper#0.2'       ,
"w2_b_2009-04-29_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_1995-03-05_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_1996-02-06_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_1994-09-24_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_1999-01-07_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_1997-07-10_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_1998-02-06_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_1999-04-21_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_2000-02-02_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_2000-06-16_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_2000-11-13_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_2001-03-23_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_2007-05-12_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_2008-02-19_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_2001-12-07_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_r_2009-04-29_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_1994-09-24_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_1995-03-05_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_1996-02-06_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_2001-12-07_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_2007-05-12_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_2008-02-19_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_1997-07-10_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_1998-02-06_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_1999-01-07_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_1999-04-21_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_2000-02-02_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_2000-06-16_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_2000-11-13_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_2001-03-23_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"acs_b_2003-01-05_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#54075',
"acs_b_2003-08-12_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#54075',
"acs_b_2003-11-28_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#54075',
"acs_b_2004-12-15_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#54075',
"acs_b_2005-04-02_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#54075',
"acs_b_2006-04-15_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#54075',
"acs_b_2006-12-06_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#54075',
"acs_r_2003-11-28_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#54075',
"acs_r_2005-09-26_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#54075',
"acs_r_2003-01-05_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#54075',
"acs_r_2003-08-12_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#54075',
"acs_r_2005-09-28_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#54075',
"acs_r_2006-12-06_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#54075',
"acs_r_2006-04-15_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#54075',
"acs_r_2006-04-29_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#54075',
}

#########
# Help functions
def sky_plt(image):
    plt.figure()
    plt.imshow(image,interpolation='nearest', cmap='afmhot', origin='lower')
    
def get_bandpar(key):
    filter_arg = bandpar_arg[key]
    dat_band = iraf.bandpar(filter_arg, Stdout=1)
    uresp = float( dat_band[1].split()[1])
    avgwl = float( dat_band[7].split()[1])
    rectw = float(dat_band[10].split()[1])
    return uresp, avgwl, rectw

# Days for scaling
SNDATE = date(1987, 2, 23)
def get_days(fname): # 1994-09-24_drz_al2.fits
    yr = int(fname[-23:-19])
    month = int(fname[-18:-16].lstrip("0"))
    day = int(fname[-15:-13].lstrip("0"))
    return (date(yr, month, day)-SNDATE).days

tt = np.zeros(len(files))
fluxes = np.zeros(len(files)) # this is not specific flux, so essentially uresp folded over throughput
uresp = np.zeros(len(files))
avgwl = np.zeros(len(files))
rectw = np.zeros(len(files))
red = np.array([False if '_b_' in ff else True for ff in files])

iraf.stsdas()
iraf.hst_calib()

for Nfilt, key in enumerate(files): # Colors
    uresp[Nfilt], avgwl[Nfilt], rectw[Nfilt] = get_bandpar(key)
    dat = fits.open(key)[0].data[540:570, 510:555]
    dat[7:24, 15:32] = 0
    fluxes[Nfilt] = uresp[Nfilt]*np.sum(dat)
    tt[Nfilt] = get_days(key)
    low = np.percentile(dat, 1)
    dat = np.where(dat < low, low, dat)
    sky_plt(np.log10(dat))

plt.figure()
plt.plot(tt[~red], fluxes[~red], '.b')
plt.plot(tt[red], fluxes[red], '.r')
plt.ylim([0, 3e-16])
plt.show()
pdb.set_trace()
    
