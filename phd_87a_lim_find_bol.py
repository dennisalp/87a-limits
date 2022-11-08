# Dennis Alp 2016-01-27
# Find point with highest optical bolometric luminosity given the upper limits.
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/Dropbox/bin/phd_87a_uplim_find_bol.py

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
WRK_DIR = "/Users/silver/box/phd/projects/87a/lim/lim"
BB_FILE = "bb.dat"
BOX = [577,597,582,602]
SET_KICK = 25 # Set kick region to a fix time in units of years
NACO_INF = 32 # Infinite aperture in pixels
NACO_H_AVGWL = 16600
NACO_H_RECTW = 3300
NACO_K_AVGWL = 21800
NACO_K_RECTW = 3500
# http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
# From Cohen et al. (2003)
REF_MAG_H = 1.133E-10 # Note that this is erg not W and A-1 not um-1
REF_MAG_K = 4.283E-11
STAR2_MASS_H = 15.186 # 2MASS magnitudes
STAR2_MASS_K = 15.233
STAR3_MASS_H = 15.446
STAR3_MASS_K = 15.174

# Reddening and model
RV = 3.1
EBV = 0.19
T_bb = 1e6
#T_bb = 5800

os.chdir(WRK_DIR) #Move to designated directory
files = sorted(glob('*al2.fits')) #Find all files. (except ALMA)
files += sorted(glob('alma*.fits')) #Put ALMA last
limits = sorted(glob('*400.0.npy')) + ['limits_alma.npy'] #Find all files.

# Source location 0-indexed
XCOO = 592.158345778
YCOO = 586.775275464
VKICK = 0.4e6
distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e16

# Please, only activate observations from the same epoch!
# Has to do with kick region
active = [] 
##active.append("ac_b_0306-00-00_drz_al2.fits")
##active.append("ac_r_0306-00-00_drz_al2.fits")
#
#active.append("ac_b_2003-11-28_drz_al2.fits")
##active.append("ac_r_2003-11-28_drz_al2.fits")
##active.append("ac_r_2005-09-26_drz_al2.fits")
##active.append("ac_b_2006-12-06_drz_al2.fits")
##active.append("ac_r_2006-12-06_drz_al2.fits")
#
active.append("na_h_2010-10-26_map_al2.fits")
active.append("na_k_2012-12-14_map_al2.fits")
#
##active.append("w2_b_0709-00-00_drz_al2.fits")
##active.append("w2_r_0709-00-00_drz_al2.fits")
#
##active.append("w2_b_1994-09-24_drz_al2.fits")
##active.append("w2_r_1994-09-24_drz_al2.fits")
#
##active.append("w32252009-12-13_drz_al2.fits")
active.append("w322r2009-12-13_drz_al2.fits")
##active.append("w322x2009-12-13_drz_al2.fits")
##active.append("w33362009-12-13_drz_al2.fits")
active.append("w333r2009-12-13_drz_al2.fits")
##active.append("w333x2009-12-13_drz_al2.fits")
##active.append("w35552009-12-13_drz_al2.fits")
active.append("w355r2009-12-13_drz_al2.fits")
##active.append("w355x2009-12-13_drz_al2.fits")
##active.append("w38142009-12-13_drz_al2.fits")
active.append("w381r2009-12-13_drz_al2.fits")
##active.append("w381x2009-12-13_drz_al2.fits")
active.append("w3_b_2009-12-12_drz_al2.fits")
active.append("w3_r_2009-12-12_drz_al2.fits")
#
active.append("w31102011-01-05_drz_al2.fits")
active.append("w31602011-01-05_drz_al2.fits")
#
##active.append("w3_b_0916-00-00_drz_al2.fits")
##active.append("w3_r_0916-00-00_drz_al2.fits")
#
##active.append("w3_b_2015-05-24_drz_al2.fits")
active.append("w3_br2015-05-24_drz_al2.fits")
##active.append("w3_bx2015-05-24_drz_al2.fits")
##active.append("w3_r_2015-05-24_drz_al2.fits")
active.append("w3_rr2015-05-24_drz_al2.fits")
##active.append("w3_rx2015-05-24_drz_al2.fits")
#
#active.append("w35022016-06-08_drz_al2.fits")
#active.append("w36452016-06-08_drz_al2.fits")
#active.append("w36572016-06-08_drz_al2.fits")
#active.append("alma_2014-09-02_212_213.fits")

bandpar_arg = {
"ac_b_0306-00-00_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#53200',
"ac_b_2003-11-28_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#52971',
"ac_b_2006-12-06_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#54075',
"ac_r_0306-00-00_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#53639', # This fold is dominated by r_2005-09-26
"ac_r_2003-11-28_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#52971',
"ac_r_2005-09-26_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#53639',
"ac_r_2006-12-06_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#54075',
"alma_2014-09-02_212_213.fits": None                              ,
"na_h_2010-10-26_map_al2.fits": None                              ,
"na_k_2012-12-14_map_al2.fits": None                              ,
"w2_b_0709-00-00_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_b_1994-09-24_drz_al2.fits": 'wfpc2,1,A2D7,f439w'              ,
"w2_r_0709-00-00_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w2_r_1994-09-24_drz_al2.fits": 'wfpc2,1,A2D7,f675w'              ,
"w31102011-01-05_drz_al2.fits": 'wfc3,ir,f110w,aper#0.2'          ,
"w31602011-01-05_drz_al2.fits": 'wfc3,ir,f160w,aper#0.2'          ,
"w32252009-12-13_drz_al2.fits": 'wfc3,uvis1,f225w,aper#0.2'       ,
"w322r2009-12-13_drz_al2.fits": 'wfc3,uvis1,f225w,aper#0.2'       ,
"w322x2009-12-13_drz_al2.fits": 'wfc3,uvis1,f225w,aper#0.2'       ,
"w33362009-12-13_drz_al2.fits": 'wfc3,uvis1,f336w,aper#0.2'       ,
"w333r2009-12-13_drz_al2.fits": 'wfc3,uvis1,f336w,aper#0.2'       ,
"w333x2009-12-13_drz_al2.fits": 'wfc3,uvis1,f336w,aper#0.2'       ,
"w35022016-06-08_drz_al2.fits": 'wfc3,uvis2,f502n,aper#0.2'       ,
"w35552009-12-13_drz_al2.fits": 'wfc3,uvis1,f555w,aper#0.2'       ,
"w355r2009-12-13_drz_al2.fits": 'wfc3,uvis1,f555w,aper#0.2'       ,
"w355x2009-12-13_drz_al2.fits": 'wfc3,uvis1,f555w,aper#0.2'       ,
"w36452016-06-08_drz_al2.fits": 'wfc3,uvis2,f645n,aper#0.2'       ,
"w36572016-06-08_drz_al2.fits": 'wfc3,uvis2,f657n,aper#0.2'       ,
"w38142009-12-13_drz_al2.fits": 'wfc3,uvis1,f814w,aper#0.2'       ,
"w381r2009-12-13_drz_al2.fits": 'wfc3,uvis1,f814w,aper#0.2'       ,
"w381x2009-12-13_drz_al2.fits": 'wfc3,uvis1,f814w,aper#0.2'       ,
"w3_b_0916-00-00_drz_al2.fits": 'wfc3,uvis2,f438w,aper#0.2'       ,
"w3_b_2009-12-12_drz_al2.fits": 'wfc3,uvis1,f438w,aper#0.2'       ,
"w3_b_2015-05-24_drz_al2.fits": 'wfc3,uvis2,f438w,aper#0.2'       ,
"w3_br2015-05-24_drz_al2.fits": 'wfc3,uvis2,f438w,aper#0.2'       ,
"w3_bx2015-05-24_drz_al2.fits": 'wfc3,uvis2,f438w,aper#0.2'       ,
"w3_r_0916-00-00_drz_al2.fits": 'wfc3,uvis2,f625w,aper#0.2'       ,
"w3_r_2009-12-12_drz_al2.fits": 'wfc3,uvis1,f625w,aper#0.2'       ,
"w3_r_2015-05-24_drz_al2.fits": 'wfc3,uvis2,f625w,aper#0.2'       ,
"w3_rr2015-05-24_drz_al2.fits": 'wfc3,uvis2,f625w,aper#0.2'       ,
"w3_rx2015-05-24_drz_al2.fits": 'wfc3,uvis2,f625w,aper#0.2'       }

#########
# Help functions
# Just reasonable plot of sky images, good for debugging
def sky_plt(image):
    plt.imshow(image,interpolation='nearest', cmap='afmhot', origin='lower')
    plt.show()

def get_bb(TT):
    bb_out = open(BB_FILE,"w")
    kB = 1.380658 # Boltzmann constant
    hh = 6.6260755e-11 # Planck constant
    cc = 2.99792458e18 # Speed of light
    la = np.logspace(1.5,4.5,10000) # Wavelength Angstrom
    
    II = 1/la**5*1/(np.exp(hh*cc/(la*kB*TT))-1)
    etot = 1./15.*(np.pi*kB*TT/(cc*hh))**4
    II /= etot

    red_mod(la, II, RV, EBV)
    np.savetxt(bb_out, np.c_[la, II])
    bb_out.close()
    return BB_FILE

# Days for scaling
SNDATE = date(1987, 2, 23)
def get_days(fname): # 1994-09-24_drz_al2.fits
    if SET_KICK > 0:
        return SET_KICK*365*24*60*60
    fname = fname.replace("0306-00-00","2004-12-06") # Set the folded to latest time
    fname = fname.replace("0709-00-00","2008-04-29")
    fname = fname.replace("0916-00-00","2013-01-01")
    
    yr = int(fname[-23:-19])
    month = int(fname[-18:-16].lstrip("0"))
    day = int(fname[-15:-13].lstrip("0"))
    return (date(yr, month, day)-SNDATE).days*24*60*60

def inside_kick(Nbox, yy, xx, key):
    telapsed = get_days(key)
    offset = np.rad2deg(np.arctan(VKICK*telapsed/(distance*1e3*pc)))*60*60*1000/25
    
    xshift = float(limits[Nbox][13:16])
    yshift = float(limits[Nbox][ 8:11])
    
    dely = (BOX[0]+yy+yshift)-YCOO
    delx = (BOX[2]+xx+xshift)-XCOO
    return np.sqrt(delx**2+dely**2) < offset

def fold_nacoh_derp(): # I overcomplicated things here obsoloete probably
    nacoh = np.loadtxt('filters/H_CONICA.dat')
    mod = np.loadtxt('bb.dat')
    
    lam  = nacoh[:,0]
    tran = nacoh[:,1]
    mlam = mod[:,0]
    flux = mod[:,1]
    
    inside = mlam > NACO_H_AVGWL-NACO_H_RECTW/2
    inside = inside & (mlam < NACO_H_AVGWL+NACO_H_RECTW/2)
    incident = simps(flux[inside], mlam[inside])
    spec = griddata(mlam, flux, lam*1e4)
    folded = spec*tran
    adu = incident/uresp[Nfilt]

def fold_nacoh(ursp):
    mod = np.loadtxt('bb.dat')
    mlam = mod[:,0]
    flux = mod[:,1]
    
    inside = mlam > NACO_H_AVGWL-NACO_H_RECTW/2
    inside = inside & (mlam < NACO_H_AVGWL+NACO_H_RECTW/2)
    incident = simps(flux[inside], mlam[inside])
    return incident/(ursp*NACO_H_RECTW)

def fold_nacok(ursp):
    mod = np.loadtxt('bb.dat')
    mlam = mod[:,0]
    flux = mod[:,1]
    
    inside = mlam > NACO_K_AVGWL-NACO_K_RECTW/2
    inside = inside & (mlam < NACO_K_AVGWL+NACO_K_RECTW/2)
    incident = simps(flux[inside], mlam[inside])
    return incident/(ursp*NACO_K_RECTW)
    
def get_nah_rsp():
    img = fits.open('na_h_2010-10-26_map_al2.fits')[0].data
    star2x = 667.96
    star2y = 672.88
    star3x = 533.03
    star3y = 555.66
    
    grid = np.indices(img.shape)
    dely = grid[0]-star2y
    delx = grid[1]-star2x
    inside = dely**2+delx**2 < NACO_INF**2
    adu = np.sum(img[inside])
    resp2 = STAR2_FLUX_H/adu

    dely = grid[0]-star3y
    delx = grid[1]-star3x
    inside = dely**2+delx**2 < NACO_INF**2
    adu = np.sum(img[inside])
    resp3 = STAR3_FLUX_H/adu
    print resp2/resp3, "H"
    return resp2

def get_nak_rsp():
    img = fits.open('na_k_2012-12-14_map_al2.fits')[0].data
    star2x = 667.96 # Same coords verified ok
    star2y = 672.88
    star3x = 533.03
    star3y = 555.66
    
    grid = np.indices(img.shape)
    dely = grid[0]-star2y
    delx = grid[1]-star2x
    inside = dely**2+delx**2 < NACO_INF**2
    adu = np.sum(img[inside])
    resp2 = STAR2_FLUX_K/adu

    dely = grid[0]-star3y
    delx = grid[1]-star3x
    inside = dely**2+delx**2 < NACO_INF**2
    adu = np.sum(img[inside])
    resp3 = STAR3_FLUX_K/adu
    print resp2/resp3, "K"
    return resp2
    
def get_bandpar(key):
    filter_arg = bandpar_arg[key]
    if filter_arg == None:
        if 'na_h' in key:
            uresp = get_nah_rsp()
            avgwl = NACO_H_AVGWL
            rectw = NACO_H_RECTW
        elif 'na_k' in key:
            uresp = get_nak_rsp()
            avgwl = NACO_K_AVGWL
            rectw = NACO_K_RECTW
        elif '212_213' in key:
            uresp = 1e-26 # UNITS ARE WRONG HERE ALMA IS PER HERTZ, OPTICAL IS PER ANGSTROM
            avgwl = 1e8
            rectw = 1e4
        return uresp, avgwl, rectw

    dat_band = iraf.bandpar(filter_arg, Stdout=1)
    uresp = float( dat_band[1].split()[1])
    avgwl = float( dat_band[7].split()[1])
    rectw = float(dat_band[10].split()[1])

    return uresp, avgwl, rectw

def merge_boxes():
    alm = np.load(limits[4])

    shift0 = np.load(limits[0])
    atmp = alm[0::2, 0::2]
    shift0 = np.concatenate((shift0, atmp), axis=2)

    shift1 = np.load(limits[1])
    atmp = alm[0::2, 1::2]
    shift1 = np.concatenate((shift1, atmp), axis=2)

    shift2 = np.load(limits[2])
    atmp = alm[1::2, 0::2]
    shift2 = np.concatenate((shift2, atmp), axis=2)

    shift3 = np.load(limits[3])
    atmp = alm[1::2, 1::2]
    shift3 = np.concatenate((shift3, atmp), axis=2)

    return [shift0, shift1, shift2, shift3]
    
lim_boxes = merge_boxes()

#########
models = []
models.append(get_bb(T_bb))

# Prepare units
STAR2_ESO_K = STAR2_MASS_K + 0.044
STAR3_ESO_K = STAR3_MASS_K + 0.044
STAR2_ESO_H = (STAR2_MASS_H-STAR2_MASS_K-0.035)/1.163+STAR2_MASS_K
STAR3_ESO_H = (STAR3_MASS_H-STAR3_MASS_K-0.035)/1.163+STAR3_MASS_K
STAR2_FLUX_H = REF_MAG_H * 10**(-0.4*STAR2_ESO_H)
STAR2_FLUX_K = REF_MAG_K * 10**(-0.4*STAR2_ESO_K)
STAR3_FLUX_H = REF_MAG_H * 10**(-0.4*STAR3_ESO_H)
STAR3_FLUX_K = REF_MAG_K * 10**(-0.4*STAR3_ESO_K)

# Definitions
img_nr = 0
phy_lim = np.zeros((len(models),5)) # physical limits, will contain limiting position and filter for each model
fluxes = np.zeros(len(files)) # this is not specific flux, so essentially uresp folded over throughput
uresp = np.zeros(len(files))
avgwl = np.zeros(len(files))
rectw = np.zeros(len(files))
iraf.stsdas()
iraf.hst_calib()

#########
# Insert psf and try to find it
for Nmod, mod in enumerate(models):
    for Nfilt, key in enumerate(files): # Colors
        uresp[Nfilt], avgwl[Nfilt], rectw[Nfilt] = get_bandpar(key)
        if 'na_h' in key:
            fluxes[Nfilt] = fold_nacoh(uresp[Nfilt])
        elif 'na_k' in key:
            fluxes[Nfilt] = fold_nacok(uresp[Nfilt])
        elif 'alma' in key:
            fluxes[Nfilt] = 1e-26*1e9 # 1 mJy in cgs over 1 GHz
        else:
            phot = iraf.calcphot(bandpar_arg[key],mod,"counts",Stdout=1)
            fluxes[Nfilt] = float(phot[6].split()[1])
        
    for Nbox, lim_box in enumerate(lim_boxes): # These are the four fractional pixel shifts
        for yy in range(0,lim_box.shape[0]):
            for xx in range(0, lim_box.shape[1]):
                faintest = [np.inf] # Most limiting filter at this point
                for Nfilt, key in enumerate(files): # Colors
                    if not key in active:
                        continue # Skip unwanted filters, not the most effective way #YOLO
                    if inside_kick(Nbox, yy, xx, key):
                        temp = lim_box[yy,xx,Nfilt]/fluxes[Nfilt]
                        print key, temp
                        # Maximum allowed flux in any filter
                        if temp < faintest[0]:
                            faintest = [temp, Nfilt, Nbox, BOX[0]+yy+1, BOX[2]+xx+1] # 1-index

                # Upper limit is the least restrictive point
                if faintest[0] > phy_lim[Nmod,0] and not faintest[0] is np.inf:
                    phy_lim[Nmod,:] = np.array(faintest)

yy = phy_lim[0,3]
xx = phy_lim[0,4]
print "Limits:"
print "Flux:",   phy_lim[0,0], "erg s-1 cm-2"
print "Filter:", files[int(phy_lim[0,1])]
print "Shift:",  limits[int(phy_lim[0,2])]
print "x:", xx
print "y:", yy
print "1st order L:", (distance*1e3*pc*100)**2*np.pi*4*phy_lim[0,0], "erg s-1"

########
# Plot
dereddening = []
fig = plt.figure(figsize=(5, 3.75))
yy = int(yy-1-BOX[0])
xx = int(xx-1-BOX[2])
lim_box = lim_boxes[int(phy_lim[0,2])]
for Nfilt, key in enumerate(files):
    if not key in active:
        continue
    lo = avgwl[Nfilt]-rectw[Nfilt]/2.
    up = avgwl[Nfilt]+rectw[Nfilt]/2.
    correction = cor_red(avgwl[Nfilt], RV, EBV)
    dereddening.append(correction)
    plt_lim = uresp[Nfilt]*correction*np.array([lim_box[yy,xx,Nfilt],lim_box[yy,xx,Nfilt]])
    plt.loglog([lo,up],plt_lim,'k')
    plt_lim = uresp[Nfilt]*np.array([lim_box[yy,xx,Nfilt],lim_box[yy,xx,Nfilt]])
    plt.loglog([lo,up],plt_lim,'r')

dereddening = np.array(dereddening)
print 'Deredden', dereddening
pdb.set_trace()
mod_spec = np.loadtxt(models[0])
plt.loglog(mod_spec[:,0], phy_lim[0,0]*mod_spec[:,1])
plt.ylabel('Flux density (erg\ s$^{-1}$\ cm$^{-2}$\ \AA$^{-1}$)')
fig.savefig('lim.pdf',bbox_inches='tight', pad_inches=0.01)
plt.show()
