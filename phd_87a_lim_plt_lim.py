# Dennis Alp 2016-01-27
# Find models that are consistent with oberved limits.
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/box/bin/phd_87a_uplim_plt_lim.py

import numpy as np
import os
import pdb
from glob import glob
from datetime import date

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import shift
from pyraf import iraf
from astropy.io import fits

from phd_87a_red import cor_red

#########
# Parameters
# Logistics
WRK_DIR = '/Users/silver/box/phd/pro/87a/lim/lim'
BOX = [577,597,582,602]
EPOCH = "not_loaded"
# 1 is 2009
# 2 is 2015
# 3 is 2006
# 4 is 2005 red only
# 5 HST IR
# 6 NACO

os.chdir(WRK_DIR) #Move to designated directory
files = sorted(glob('*al2.fits')) #Find all files. (except ALMA)
files += sorted(glob('alma*.fits')) #Put ALMA last
limits = sorted(glob('*800.0.npy')) + ['limits_alma.npy'] #Find all files.

# Source location 0-indexed
XCOO = 592.158345778
YCOO = 586.775275464
VKICK = 0.8e6
distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e16 # m SI!!!
sb = 5.67051e-5 # Stefan-Boltzmann
Lsun = 3.826e33 # erg s-1
Rsun = 6.96e10 # cm
hh = 6.626e-27
cc = 2.99792458e10
kB = 1.380658e-16

# Redenning
RV = 3.1
EBV = 0.19

# Parameters for calibration of VLT NACO
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

STAR2_ESO_K = STAR2_MASS_K + 0.044
STAR3_ESO_K = STAR3_MASS_K + 0.044
STAR2_ESO_H = (STAR2_MASS_H-STAR2_MASS_K-0.035)/1.163+STAR2_MASS_K
STAR3_ESO_H = (STAR3_MASS_H-STAR3_MASS_K-0.035)/1.163+STAR3_MASS_K
STAR2_FLUX_H = REF_MAG_H * 10**(-0.4*STAR2_ESO_H)
STAR2_FLUX_K = REF_MAG_K * 10**(-0.4*STAR2_ESO_K)
STAR3_FLUX_H = REF_MAG_H * 10**(-0.4*STAR3_ESO_H)
STAR3_FLUX_K = REF_MAG_K * 10**(-0.4*STAR3_ESO_K)

# Please, only activate observations from the same epoch!
# Has to do with kick region
active = [] 
#active.append("ac_b_0306-00-00_drz_al2.fits")
#active.append("ac_r_0306-00-00_drz_al2.fits")

#active.append("ac_b_2003-11-28_drz_al2.fits")
#active.append("ac_r_2003-11-28_drz_al2.fits")
#active.append("ac_r_2005-09-26_drz_al2.fits")
#active.append("ac_b_2006-12-06_drz_al2.fits")
#active.append("ac_r_2006-12-06_drz_al2.fits")

#active.append("na_h_2010-10-26_map_al2.fits")
#active.append("na_k_2012-12-14_map_al2.fits")

#active.append("w2_b_0709-00-00_drz_al2.fits")
#active.append("w2_r_0709-00-00_drz_al2.fits")

#active.append("w2_b_1994-09-24_drz_al2.fits")
#active.append("w2_r_1994-09-24_drz_al2.fits")

##active.append("w32252009-12-13_drz_al2.fits")
#active.append("w322r2009-12-13_drz_al2.fits")
##active.append("w322x2009-12-13_drz_al2.fits")
##active.append("w33362009-12-13_drz_al2.fits")
#active.append("w333r2009-12-13_drz_al2.fits")
##active.append("w333x2009-12-13_drz_al2.fits")
##active.append("w35552009-12-13_drz_al2.fits")
#active.append("w355r2009-12-13_drz_al2.fits")
##active.append("w355x2009-12-13_drz_al2.fits")
##active.append("w38142009-12-13_drz_al2.fits")
#active.append("w381r2009-12-13_drz_al2.fits")
##active.append("w381x2009-12-13_drz_al2.fits")
#active.append("w3_b_2009-12-12_drz_al2.fits")
#active.append("w3_r_2009-12-12_drz_al2.fits")

#active.append("w31102011-01-05_drz_al2.fits")
#active.append("w31602011-01-05_drz_al2.fits")

#active.append("w3_b_0916-00-00_drz_al2.fits")
#active.append("w3_r_0916-00-00_drz_al2.fits")

##active.append("w3_b_2015-05-24_drz_al2.fits")
#active.append("w3_br2015-05-24_drz_al2.fits")
##active.append("w3_bx2015-05-24_drz_al2.fits")
##active.append("w3_r_2015-05-24_drz_al2.fits")
#active.append("w3_rr2015-05-24_drz_al2.fits")
##active.append("w3_rx2015-05-24_drz_al2.fits")

#active.append("w35022016-06-08_drz_al2.fits")
#active.append("w36452016-06-08_drz_al2.fits")
#active.append("w36572016-06-08_drz_al2.fits")

active.append("alma_2014-09-02_212_213.fits")
#active.append("alma_2014-09-02_212_247.fits")
active.append("alma_2014-09-02_233_233.fits")
active.append("alma_2014-09-02_247_247.fits")

bandpar_arg = {
"ac_b_0306-00-00_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#53200',
"ac_b_2003-11-28_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#52971',
"ac_b_2006-12-06_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#54075',
"ac_r_0306-00-00_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#53639', # This fold is dominated by r_2005-09-26
"ac_r_2003-11-28_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#52971',
"ac_r_2005-09-26_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#53639',
"ac_r_2006-12-06_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#54075',
"alma_2014-09-02_212_213.fits": None                              ,
"alma_2014-09-02_233_233.fits": None                              ,
"alma_2014-09-02_247_247.fits": None                              ,
"alma_2014-09-02_212_247.fits": None                              ,
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

# Days for scaling
SNDATE = date(1987, 2, 23)
def get_days(fname): # 1994-09-24_drz_al2.fits
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

def get_bandpar(key):
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
        elif 'alma' in key:
            uresp = 1.
            avgwl = 1.
            if '212_213' in key:
                rectw = 1/4.22545340531e-05
            elif '212_247' in key:
                rectw = 1/3.16911325604e-05
            elif '233_233' in key:
                rectw = 1/5.92417818956e-05
            elif '247_247' in key:
                rectw = 1/4.43683329131e-05
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
# Definitions
img_nr = 0
phy_lim = np.zeros((len(files),4)) # physical limits, will contain limiting position and filter for each model
bol_lim = np.zeros(4)
uresp = np.zeros(len(files))
avgwl = np.zeros(len(files))
rectw = np.zeros(len(files))
iraf.stsdas()
iraf.hst_calib()
REFERENCE = active[0]

#########
for Nfilt, key in enumerate(files): # Colors
    uresp[Nfilt], avgwl[Nfilt], rectw[Nfilt] = get_bandpar(key)

for Nbox, lim_box in enumerate(lim_boxes): # These are the four fractional pixel shifts
    lim_box = lim_box*cor_red(avgwl, RV, EBV)
    for yy in range(0,lim_box.shape[0]):
        for xx in range(0, lim_box.shape[1]):
            bol_tmp = 0
            for Nfilt, key in enumerate(files): # Colors
                if not key in active:
                    continue # Skip unwanted filters, not the most effective way #YOLO
                if inside_kick(Nbox, yy, xx, key):
                    temp = lim_box[yy, xx, Nfilt]*uresp[Nfilt]
                    bol_tmp += temp*rectw[Nfilt]
                    # Maximum allowed flux in any filter
                    if temp > phy_lim[Nfilt, 0]:
                        phy_lim[Nfilt,:] = np.array([temp, yy, xx, Nbox])

            if bol_tmp > bol_lim[0]:
                bol_lim = np.array([bol_tmp, yy, xx, Nbox])

act_ind =  np.where(phy_lim[:,0] > 0)
lim_box = lim_boxes[int(bol_lim[3])]*cor_red(avgwl, RV, EBV)
print 'Free points', phy_lim[act_ind,0]
print 'Max point', lim_box[int(bol_lim[1]), int(bol_lim[2]), act_ind]*uresp[act_ind]
print 'Ratio', phy_lim[act_ind,0]/lim_box[int(bol_lim[1]), int(bol_lim[2]), act_ind]/uresp[act_ind]
print 'Bolometric lim', bol_lim
print 'Bolometric of free points', np.sum(phy_lim[act_ind,0]*rectw[act_ind])
print 'Bolometric ratio', np.sum(phy_lim[act_ind,0]*rectw[act_ind])/bol_lim[0]
print 'Coords of free points', phy_lim[act_ind,1:3]
print "1st order L:", (distance*1e3*pc*100)**2*np.pi*4*bol_lim[0], "erg s-1"
print "1st order L:", (distance*1e3*pc*100)**2*np.pi*4*bol_lim[0]/Lsun, "L_Sun"
#pdb.set_trace()
# Save for merged plot
save = avgwl[act_ind]
save = np.vstack((save, rectw[act_ind]))
save = np.vstack((save, lim_box[int(bol_lim[1]), int(bol_lim[2]), act_ind]*uresp[act_ind]))
np.save("phy_lim_" + EPOCH + ".npy", save.T)



################################################################
# Plot
def get_shifts(Nbox):
    xx = limits[Nbox][13:16]
    yy = limits[Nbox][8:11]
    return float(yy), float(xx)

########
# Series of plots of all positions
theta = np.arange(0,2*np.pi,0.001)
for Nfilt in range(0, phy_lim.shape[0]):
    key = files[Nfilt]
    if key in active:
        Nbox = int(phy_lim[Nfilt, 3])
        yy, xx = get_shifts(Nbox)
        
        sky_img = fits.open(key)[0].data
        sky_img = fits.open(key)[0].data[0,0] # ALMA
        BOX = [200,400,200,400] # ALMA
        sky_img = sky_img[BOX[0]:BOX[1],BOX[2]:BOX[3]]
        plt.imshow(sky_img,interpolation='nearest', cmap='afmhot', origin='lower')

        telapsed = get_days(key)
        offset = np.rad2deg(np.arctan(VKICK*telapsed/(distance*1e3*pc)))*60*60*1000/25
        plt.plot(XCOO-BOX[2]+offset*np.cos(theta),YCOO-BOX[0]+offset*np.sin(theta),'k')
        plt.plot(XCOO-BOX[2],YCOO-BOX[0],'wo')
        
        plt.plot(phy_lim[Nfilt, 2]+xx, phy_lim[Nfilt,1]+yy,'ok')
        plt.savefig("plt/" + key[:-5] + ".pdf",bbox_inches='tight', pad_inches=0.01)
#        plt.show()
        plt.close()

        plt.imshow(lim_boxes[Nbox][:,:,Nfilt],interpolation='nearest', cmap='afmhot', origin='lower')

        plt.plot(XCOO-BOX[2]+offset*np.cos(theta),YCOO-BOX[0]+offset*np.sin(theta),'k')
        plt.plot(XCOO-BOX[2],YCOO-BOX[0],'wo')
        
        plt.plot(phy_lim[Nfilt, 2]+xx, phy_lim[Nfilt,1]+yy,'ok')
        plt.savefig("plt/" + key[:-8] + "lim.pdf",bbox_inches='tight', pad_inches=0.01)
#        plt.show()
        plt.close()

########
# Single plot for all positions        
for Nfilt in range(0, phy_lim.shape[0]):
    key = files[Nfilt]
    if key in active:
        Nbox = int(phy_lim[Nfilt, 3])
        yy, xx = get_shifts(Nbox)
        plt.plot(phy_lim[Nfilt, 2] + xx, phy_lim[Nfilt,1] + yy,'ok')
        plt.annotate(files[Nfilt][:9], xy=phy_lim[Nfilt,-2:0:-1]+0.5*np.random.uniform(0,1))
        
sky_img = fits.open(key)[0].data
sky_img = fits.open(key)[0].data[0,0] # ALMA
sky_img = sky_img[BOX[0]:BOX[1],BOX[2]:BOX[3]]
plt.imshow(sky_img,interpolation='nearest', cmap='afmhot', origin='lower')

telapsed = get_days(REFERENCE)
offset = np.rad2deg(np.arctan(VKICK*telapsed/(distance*1e3*pc)))*60*60*1000/25
plt.plot(XCOO-BOX[2]+offset*np.cos(theta),YCOO-BOX[0]+offset*np.sin(theta),'k')
plt.plot(XCOO-BOX[2],YCOO-BOX[0],'wo')
plt.savefig("plt/0all.pdf",bbox_inches='tight', pad_inches=0.01)
plt.close()

########
# 400
# [[  8.93462541e-18,  4.12659688e-18,  3.43879761e-18,  2.31490265e-18,  4.65800603e-18,  4.98632241e-18]] 2.54403225278e-14
# [[  8.73213613e-18,  4.12659688e-18,  3.43879761e-18,  2.31490265e-18,  4.65800603e-18,  4.98632241e-18]]
# [  2.53462622e-14,  9.00000000e+00,  1.10000000e+01,  3.00000000e+00]

#Free points [[  8.73213613e-18   3.63281605e-18   2.78510049e-18   1.86208513e-18
#    3.85958065e-18   3.91173408e-18]]
#Max point [[  8.63108211e-18   3.28626317e-18   2.64981604e-18   1.85295819e-18
#    3.85958065e-18   3.91173408e-18]]
#Ratio [[ 1.01170815  1.105455    1.05105428  1.00492561  1.          1.        ]]
#Bolometric lim [  2.08458267e-14   1.00000000e+01   1.00000000e+01   3.00000000e+00]
#Bolometric of free points 2.12961431687e-14
#Bolometric ratio 1.02160223394

########
# 1600
# [[  1.16648484e-17,  4.93478707e-18,  4.29162343e-18,  2.50121543e-18,  5.83323388e-18,  5.60758175e-18]] 3.03753840639e-14
# [[  1.01708400e-17,  4.80565965e-18,  4.18568515e-18,  2.41943526e-18,  5.64131190e-18,  5.44864396e-18]]
# [  2.89705265e-14,  1.30000000e+01,  1.10000000e+01,  3.00000000e+00]

#Free points [[  1.16648484e-17   4.82320793e-18   4.29162343e-18   2.35531331e-18
#    5.37549957e-18   5.17344639e-18]]
#Max point [[  1.00671266e-17   4.15364658e-18   3.93600119e-18   1.98569563e-18
#    4.10027057e-18   5.17344639e-18]]
#Ratio [[ 1.15870683  1.16119844  1.09035115  1.18614014  1.31101094  1.        ]]
#Bolometric lim [  2.61670120e-14   1.10000000e+01   1.50000000e+01   3.00000000e+00]
#Bolometric of free points 2.91712393719e-14
#Bolometric ratio 1.11480972093

########
# 800
# [[  1.06323619e-17,  4.82320793e-18,  3.82540306e-18,  2.31490265e-18,  4.93614394e-18,  5.29225003e-18]] 2.78076673059e-14
# [[  9.78479882e-18,  4.42684378e-18,  3.82540306e-18,  2.20002549e-18,  4.93614394e-18,  5.29225003e-18]]
# [  2.70283603e-14,  1.10000000e+01,  1.20000000e+01,  3.00000000e+00]

#Free points [[  9.78479882e-18   4.82320793e-18   3.82540306e-18   2.01514940e-18
#    4.54123635e-18   4.21471766e-18]]
#Max point [[  9.78479882e-18   3.93103728e-18   3.82540306e-18   2.01514940e-18
#    4.05622263e-18   4.21062229e-18]]
#Ratio [[ 1.          1.22695553  1.          1.          1.11957276  1.00097263]]
#Bolometric lim [  2.43577294e-14   1.10000000e+01   1.20000000e+01   3.00000000e+00]
#Bolometric of free points 2.51196656791e-14
#Bolometric ratio 1.03128108776

################################################################
#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
# Merge several epochs and plot the image limits
files = sorted(glob('phy_lim_??.npy')) #Find all files.
dereddening = []

lims = np.zeros((0,3))
for ff in files:
    if ff == 'phy_lim_06.npy':
#        print 'Adding 10% to NACO values.'
        tmp = np.load(ff)
        tmp[:, 2] = 1.*tmp[:,2]
        lims = np.vstack((lims, tmp))
    else:
        lims = np.vstack((lims, np.load(ff)))

# Initiate plot
fig = plt.figure(figsize=(5, 3.75))
col = [(0.5, 0., 0.5), (0.5, 0., 1.), (0., 1., 0.), (0.5, 0., 0.), (0., 0., 1.), (1., 0., 0.), (0., 0., 1.), (1., 0., 0.), 'k','k','k','k','k','k','k','k','k','k','k']
out = []
for idx, lim in enumerate(lims):
    if idx in [8, 9, 10]:
        continue
    lo = lim[0]-lim[1]/2.
    up = lim[0]+lim[1]/2.
    plt_lim = lim[2]*np.ones(2)
#    plt.semilogy([lo,up], plt_lim, 'r')

    correction = cor_red(lim[0], RV, EBV)
    dereddening.append(correction)
    plt_lim *= correction
    plt.loglog([lo,up], plt_lim, lw=3, color = col[idx])
    out.append([lo, plt_lim[0]])
    out.append([up, plt_lim[1]])

# Fudge spectrum for bolometric limit
plt.loglog([1000, 1000, 3000, 3000], [1e-99, 4e-17, 4e-17, 1e-99], color='#bcbd22', ls='-.', lw=3, zorder=-100)

# Companion
TT = 5780
AA = 4*np.pi*Rsun**2
lam = np.linspace(1000e-8, 10e-4, 1000)
companion2 = AA*sb*TT**4
companion2 = companion2/(4*np.pi*(distance*1e3*pc*100)**2)
companion = 2*hh*cc**2/lam**5*1/(np.exp(hh*cc/(lam*kB*TT))-1)
companion = AA*np.pi*companion
companion = companion/(4*np.pi*(distance*1e3*pc*100)**2)
plt.loglog(lam/1e-8, companion*1e-8, color='#ff7f0e', ls='--', lw=3)

TT = 6300
AA = 4*np.pi*(1.4*Rsun)**2
lam = np.linspace(1000e-8, 10e-4, 1000)
companion2 = AA*sb*TT**4
companion2 = companion2/(4*np.pi*(distance*1e3*pc*100)**2)
companion = 2*hh*cc**2/lam**5*1/(np.exp(hh*cc/(lam*kB*TT))-1)
companion = AA*np.pi*companion
companion = companion/(4*np.pi*(distance*1e3*pc*100)**2)
#plt.loglog(lam/1e-8, companion*1e-8, color='#bcbd22', ls='-.', lw=3)

TT = 8000
AA = 4*np.pi*(0.6*Rsun)**2
lam = np.linspace(1000e-8, 10e-4, 1000)
companion2 = AA*sb*TT**4
companion2 = companion2/(4*np.pi*(distance*1e3*pc*100)**2)
companion = 2*hh*cc**2/lam**5*1/(np.exp(hh*cc/(lam*kB*TT))-1)
companion = AA*np.pi*companion
companion = companion/(4*np.pi*(distance*1e3*pc*100)**2)
#plt.loglog(lam/1e-8, companion*1e-8, color='#e377c2', ls=':', lw=3)

## Netron star
#GG = 6.67259e-8 # cm3 g-1 s-2
#Msun = 1.989e33 # g
#RR = 1e6
#MM = 1.4
#RS = 2*GG*MM*Msun/cc**2
#gr = np.sqrt(1-RS/RR)
#
#TT = 4e6
#AA = 4*np.pi*RR**2
#lam = np.linspace(1000e-8, 10e-4, 1000)
#companion2 = AA*sb*TT**4*gr**2
#companion2 = companion2/(4*np.pi*(distance*1e3*pc*100)**2)
#companion = 2*hh*cc**2/lam**5*1/(np.exp(hh*cc/(lam*kB*TT))-1)
#companion = AA*np.pi*companion*gr**2
#companion = companion/(4*np.pi*(distance*1e3*pc*100)**2)
#plt.loglog(lam/1e-8, companion*1e-8, color='#e377c2', ls=':', lw=3)

# STIS
sts_par = np.array([7.4822863429878992e-15, -0.9527313929920449])
sts_int = np.arange(5300,10000)
sts_val = sts_par[0]*sts_int**sts_par[1]
plt.plot(sts_int, sts_val, color=(0.6,0.6,0.6), lw=3)
out.append([sts_int[0], sts_val[0]])
out.append([sts_int[-1], sts_val[-1]])
# SINFONI
#    Final bin lims flx [  2.43663528e-19   1.39144179e-19   1.82898534e-19   1.20305289e-19]
snf_lim = np.array([2.43663528e-19, 1.39144179e-19, 1.82898534e-19, 1.20305289e-19])
snf_int = np.array([[15150, 15800],
       [16900, 18125],
       [19875, 22725],
       [22725, 23825]])

for snf_ii in range(0, snf_int.shape[0]):
    plt.plot(snf_int[snf_ii], snf_lim[snf_ii]*np.ones(2), color=(0.6,0.6,0.6), lw=3)
    out.append([snf_int[snf_ii][0], snf_lim[snf_ii]])
    out.append([snf_int[snf_ii][1], snf_lim[snf_ii]])

np.savetxt('/Users/silver/box/phd/pro/87a/pre/lim/lim.txt', out)

# Cosmetics
plt.xlim([9e2, 3e4])
plt.ylim([1e-19, 5e-17])
plt.gca().set_xticks([1e3, 2e3, 5e3, 1e4, 2e4])
import matplotlib
plt.gca().get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
dereddening = np.array(dereddening)
print 'Deredden', dereddening
plt.ylabel('Flux density (erg\ s$^{-1}$\ cm$^{-2}$\ \AA$^{-1}$)')
plt.xlabel('Wavelength (\AA)')
fig.savefig('/Users/silver/box/phd/pro/87a/lim/art/figs/img_lim.pdf',bbox_inches='tight', pad_inches=0.03)
plt.show()

# For table values
tab_der = cor_red(lims[:,0], RV, EBV)
tab_lo = lims[:,0]-lims[:,1]/2.
tab_up = lims[:,0]+lims[:,1]/2.
tab_flx = lims[:,2]*tab_der
tab_lum = (distance*1e3*pc*100)**2*np.pi*4*tab_flx*(tab_up-tab_lo)
tab_sun = tab_lum/Lsun
tab_ord = [0, 1, 4, 2, 5, 3, 11, 12, 6, 7, 13, 14]
print 'FOR TABLE VALUES ON LAST LINES IN SCRIPT'

pdb.set_trace()
