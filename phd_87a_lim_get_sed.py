# Dennis Alp 2017-02-13
# Get the spectral energy distribution of several images
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/Dropbox/bin/phd_87a_uplim_get_sed.py

import numpy as np
import os
import pdb
from glob import glob
from datetime import date

import matplotlib.pyplot as plt
from pyraf import iraf
from astropy.io import fits

#########
# Parameters
# Logistics
WRK_DIR = "/Users/silver/Dropbox/phd/projects/87a/uplim/sed"
BOX = [577,597,582,602]

os.chdir(WRK_DIR) #Move to designated directory
files = sorted(glob('*al2.fits')) #Find all files.

# Source location 0-indexed
BLOCK = 4
XCOO = 592.158345778
YCOO = 586.775275464
VKICK = 1.6e6
distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e16

bandpar_arg = {
"ac_b_0306-00-00_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#53200',
"ac_b_2003-11-28_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#52971',
"ac_b_2006-12-06_drz_al2.fits": 'acs,hrc,f435w,aper#0.2,mjd#54075',
"ac_r_0306-00-00_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#53639', # This fold is dominated by r_2005-09-26
"ac_r_2003-11-28_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#52971',
"ac_r_2005-09-26_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#53639',
"ac_r_2006-12-06_drz_al2.fits": 'acs,hrc,f625w,aper#0.2,mjd#54075',
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

sigma = {
"ac_b_0306-00-00_drz_al2.fits": 0.00701834489262 ,
"ac_b_2003-11-28_drz_al2.fits": 0.0170122807654 ,
"ac_b_2006-12-06_drz_al2.fits": 0.0163254725161 ,
"ac_r_0306-00-00_drz_al2.fits": 0.00701116831843 ,
"ac_r_2003-11-28_drz_al2.fits": 0.0311758177859 ,
"ac_r_2005-09-26_drz_al2.fits": 0.00842018783725 ,
"ac_r_2006-12-06_drz_al2.fits": 0.0296670558917 ,
"na_h_2010-10-26_map_al2.fits": 2.27982611406 ,
"na_k_2012-12-14_map_al2.fits": 2.2356313251 ,
"w2_b_0709-00-00_drz_al2.fits": 0.000553888968355 ,
"w2_b_1994-09-24_drz_al2.fits": 0.00428215404537 ,
"w2_r_0709-00-00_drz_al2.fits": 0.00176941685158 ,
"w2_r_1994-09-24_drz_al2.fits": 0.0127718854892 ,
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
"w3_b_2015-05-24_drz_al2.fits": 0.0094900377064 ,
"w3_br2015-05-24_drz_al2.fits": 0.00974715488903 ,
"w3_bx2015-05-24_drz_al2.fits": 0.00975035655871 ,
"w3_r_0916-00-00_drz_al2.fits": 0.00691686632371 ,
"w3_r_2009-12-12_drz_al2.fits": 0.0122142571276 ,
"w3_r_2015-05-24_drz_al2.fits": 0.0177253756152 ,
"w3_rr2015-05-24_drz_al2.fits": 0.0180942252429 ,
"w3_rx2015-05-24_drz_al2.fits": 0.0181124092083 }

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

def inside_kick(yy, xx, key):
    telapsed = get_days(key)
    offset = np.rad2deg(np.arctan(VKICK*telapsed/(distance*1e3*pc)))*60*60*1000/25
    
    dely = (BOX[0]+yy)-YCOO
    delx = (BOX[2]+xx)-XCOO
    return np.sqrt(delx**2+dely**2) < offset

def get_bandpar(key):
    filter_arg = bandpar_arg[key]
    if filter_arg == None:
        uresp = 1e-20
        if 'na_h' in key:
            avgwl = 16600
            rectw = 3300
        elif 'na_k' in key:
            avgwl = 21800
            rectw = 3500
        return uresp, avgwl, rectw

    dat_band = iraf.bandpar(filter_arg, Stdout=1)
    uresp = float( dat_band[1].split()[1])
    avgwl = float( dat_band[7].split()[1])
    rectw = float(dat_band[10].split()[1])

    return uresp, avgwl, rectw

#########
# Definitions
uresp = np.zeros(len(files))
avgwl = np.zeros(len(files))
rectw = np.zeros(len(files))
seds  = np.zeros((BOX[1]-BOX[0], BOX[3]-BOX[2], len(files)))
chi2map= np.zeros((BOX[1]-BOX[0], BOX[3]-BOX[2]))
iraf.stsdas()
iraf.hst_calib()

#########
for Nfilt, key in enumerate(files): # Colors
    uresp[Nfilt], avgwl[Nfilt], rectw[Nfilt] = get_bandpar(key)
    seds[:,:,Nfilt] = fits.open(key)[0].data[BOX[0]:BOX[1],BOX[2]:BOX[3]]

norm = np.sum(seds, axis=2)
#for ii in range(0, len(files)):
#    seds[:,:,ii] /= seds[:,:,5]
#    seds[:,:,ii] /= norm
    
order = avgwl.argsort()
griy, grix = np.indices(chi2map.shape)
inside = inside_kick(griy, grix, key)
ejecta = np.mean(uresp*seds[inside], axis=(0))
ej_std = np.std(uresp*seds[inside], axis=(0))
ebar = []
for key in files:
    ebar.append(sigma[key])
ebar = np.array(ebar)

for yy in range(0,BOX[1]-BOX[0]):
    for xx in range(0,BOX[3]-BOX[2]):
        for Nfilt, key in enumerate(files):
            if inside_kick(yy, xx, key):
                lo = avgwl[Nfilt]-rectw[Nfilt]/2.
                up = avgwl[Nfilt]+rectw[Nfilt]/2.
                val = seds[yy,xx,Nfilt]
                plt_lim = uresp[Nfilt]*np.array([val,val])
#                plt.plot(np.array([lo,up])+xx*5000,plt_lim+4e-19*yy,'k')
        if inside_kick(yy, xx, key):
#            plt.errorbar(avgwl[order]+xx*5000, uresp[order]*seds[yy,xx,order]+4e-19*yy, yerr=uresp[order]*ebar[order], color='k')
            plt.plot(avgwl[order]+xx*6250, (uresp[order]*seds[yy,xx,order]-ejecta[order])/ej_std[order]+5*yy, color='k')
            
            chi2map[yy,xx] = np.sum(((uresp*seds[yy,xx,:]-ejecta)/ej_std)**2)

plt.show()

plt.plot(XCOO-BOX[2], YCOO-BOX[0], 'go')
sky_plt(chi2map)

########
# Grouped
plt.gca().set_xticks(np.arange(0,BOX[3]-BOX[2],BLOCK)-0.5)
plt.gca().set_yticks(np.arange(0,BOX[1]-BOX[0],BLOCK)-0.5)
plt.grid(linewidth=3, color="g", linestyle="-")
sky_plt(fits.open("w3_r_2009-12-12_drz_al2.fits")[0].data[BOX[0]:BOX[1],BOX[2]:BOX[3]])

ejecta = np.mean(uresp*seds, axis=(0,1))

for yy in range(0,BOX[1]-BOX[0], BLOCK):
    for xx in range(0,BOX[3]-BOX[2], BLOCK):
        val = np.mean(seds[yy:yy+BLOCK , xx:xx+BLOCK,:], axis=(0,1))
        plt.plot(avgwl[order]+xx*1250,  uresp[order]*val[order]+1e-19*yy, color='k')
        plt.plot(avgwl[order]+xx*1250,            ejecta[order]+1e-19*yy, color='g', linestyle="--")
plt.show()

ejecta /= np.sum(ejecta)
for yy in range(0,BOX[1]-BOX[0], BLOCK):
    for xx in range(0,BOX[3]-BOX[2], BLOCK):
        val = np.mean(seds[yy:yy+BLOCK , xx:xx+BLOCK,:], axis=(0,1))
        val = uresp[order]*val[order]
        val /= np.sum(val)
        plt.plot(avgwl[order]+xx*1250, (val-ejecta[order])+0.03*yy, color='k')
        plt.plot(avgwl[order]+xx*1250,np.zeros(ejecta.shape[0])+0.03*yy, color='g', linestyle="--")
plt.show()

plt.gca().set_xticks(np.arange(0,BOX[3]-BOX[2],BLOCK)-0.5)
plt.gca().set_yticks(np.arange(0,BOX[1]-BOX[0],BLOCK)-0.5)
plt.grid(linewidth=3, color="g", linestyle="-")
sky_plt(fits.open("w322r2009-12-13_drz_al2.fits")[0].data[BOX[0]:BOX[1],BOX[2]:BOX[3]])
