# Dennis Alp 2017-03-24
# Align STIS images and extract useful data
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/$USER/box/bin/phd_87a_uplim_ali.py

#Power law parameters: [  2.33970706e-16  -6.03139744e-01]
#chi2: 75.8293728242 	dof: 57 red: 1.33033987411
#chi2: 74.4410518432 	dof: 57 red: 1.30598336567
#Observed flux:	5.07842629326e-15
#Spectral flux density:	1.08051623261e-18
#Luminosity:	1.59286888927e+33
#Spectral luminosity density:	3.38908274313e+29
#Solar luminosities:	0.416327467139

#Power law parameters: [  2.08459079e-15  -8.87125361e-01]
#chi2: 76.2700985884 	dof: 57 red: 1.33807190506
#chi2: 75.141600385 	dof: 57 red: 1.31827369096
#Observed flux:	3.6119491375e-15
#Spectral flux density:	7.68499816489e-19
#Luminosity:	1.13290241475e+33
#Spectral luminosity density:	2.41043066968e+29
#Solar luminosities:	0.296106224451

#Power law parameters: [  9.05565294e-18  -1.28948884e-01]
#chi2: 75.4055708189 	dof: 57 red: 1.32290475121
#chi2: 72.0512268283 	dof: 57 red: 1.26405661102
#Observed flux:	1.34663464996e-14
#Spectral flux density:	2.86518010629e-18
#Luminosity:	4.22377389228e+33
#Spectral luminosity density:	8.98675296229e+29
#Solar luminosities:	1.1039659938

import os
import pdb
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from astropy.io import fits

from phd_87a_red import cor_red

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

#########
# Help functions
# Just reasonable plot of sky images, good for debugging
def sky_plt(image):
    plt.imshow(image,interpolation='nearest', cmap='afmhot', origin='lower')
    plt.show()

def mk_fits(image, error, output):
    hdulist = fits.HDUList()
    hdulist.append(fits.ImageHDU(image))
    hdulist.append(fits.ImageHDU(error))
    hdulist.writeto(output, clobber=True)
    hdulist.close()

def gauss(yy, aa, y0, sig):
    return aa*np.exp(-(yy-y0)**2/(2*sig**2))

def get_ringy(idx):
    slit = int(np.round(specX[idx]))-2 # lower edge of pixel
    cut = np.mean(ref_img[:, slit:slit+4], axis=1) # slice along y in image
    yy = np.arange(0, cut.shape[0])
    cut = np.where(np.abs(yy-ringY[0]) < 10, cut, 0)
    par, cov = curve_fit(gauss, yy, cut, p0=[10, ringY[0], 3])
    bot = par[1]

#    plt.plot(yy, cut, 'b')
#    plt.plot(yy, gauss(yy, *par), 'r')
#    plt.show()

    cut = np.mean(ref_img[:, slit:slit+4], axis=1) # slice along y in image
    yy = np.arange(0, cut.shape[0])
    cut = np.where(np.abs(yy-ringY[1]) < 6, cut, 0)
    par, cov = curve_fit(gauss, yy, cut, p0=[10, ringY[1], 3])
    top = par[1]

    return (bot, top)

   
#########
# HST STIS Parameters
WRK_DIR = "/Users/silver/dat/hst/87a/stis/"
REFERENCE = "/Users/silver/box/phd/pro/87a/lim/data/w3_rr2015-05-24_drz_al2.fits"
REFERENCE = "/Users/silver/box/phd/pro/87a/lim/data/w3_r_2014-06-15_drz_al2.fits"
REF_STAR = (66.93, 288.88) # 0-indexed, imexamined the reference image
dec = -69-16/60.-18.36/3600 # Declination at 1 (reference for the observation)
RA_OFFSET = np.array([-2.4435, -2.4625, -2.4815, -2.5005, -2.5195]) # in seconds
RA_OFFSET = RA_OFFSET*np.cos(np.deg2rad(dec))*15 # to arcsec
DEL_LAM = 4.882 # Difference in A in x axis, from .fits
REF_LAM = 600 # Reference x pixel (for unit). 0 indexed
REF_VAL = 7751 # Reference value in A

# Reddening
RV = 3.1
EBV = 0.19
distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e16
Lsun = 3.826e33 # erg s-1

# Source location 0-indexed
XCOO = 592.158345778
YCOO = 586.775275464
ringY = (563, 610) # y coordinate of the ring

# y coordinate of some lines in STIS, 0-index! using imexam
lineYup = [[240.37, 240.44, 240.12, 240.14, 240.09, 240.26, 240.01, 240.08],
           [240.34, 240.24, 240.15, 240.16, 240.16, 240.25, 240.13, 240.10],
           [240.62, 240.62, 240.47, 240.51, 240.54, 240.54, 240.37, 240.37, 240.41, 240.41, 240.41, 240.44],
           [240.60, 240.72, 240.57, 240.67, 240.65, 240.68, 240.62, 240.61, 240.52],
           [240.23, 240.26, 240.12, 240.16, 240.13, 240.17, 240.13, 240.16]]

lineYlo = [[217.94, 219.01, 219.01, 217.91, 217.54, 217.79],
           [216.91, 216.89, 216.92, 217.11, 216.51, 216.57, 216.99],
           [217.50, 217.50, 217.49, 217.54, 217.48, 217.53, 217.30, 217.36, 217.38, 217.41],
           [217.66, 217.68, 217.72, 217.71, 217.71, 217.56, 217.51, 217.54, 217.52, 217.56, 217.60, 217.50],
           [217.60, 217.52, 217.55, 217.60, 217.42, 217.56, 217.46, 217.54, 217.53]]
    

################################################################
# Align images
ref_img = fits.open(REFERENCE)[0].data
os.chdir(WRK_DIR)
files = sorted(glob("*al2.fits")) #Find all files.
nn = len(files)
fit_ringy = []

specX = REF_STAR[0]-RA_OFFSET/0.025 # x positions of spectra wrt SN 1987A imalign
specY = np.zeros(nn) # y...
slit = np.diff(RA_OFFSET)[0]/0.025 # Exact slit, not needed, always +-1 pixel
yy = np.zeros(len(files))
nn = fits.open(files[0])[0].data.shape[1]

bkg_sum = np.zeros((len(files), 7, nn, 2))
bkg_err = np.zeros((len(files), 7, nn, 2))

for idx, key in enumerate(files):
# Following lines for inverting Y
#    spc = fits.open(key)[1].data
#    err = fits.open(key)[2].data
#    unit_conv = 0.10087169*fits.getval(key, 'CD2_2', ext=1)*3600
#    mk_fits(unit_conv*spc[::-1,:], unit_conv*err[::-1,:], key.replace('una', 'al2'))

# Concerning ejecta, i.e. the spectra we want
    fit_ringy.append(get_ringy(idx))
    yup = np.mean(np.array(lineYup[idx]))
    ylo = np.mean(np.array(lineYlo[idx]))

    delpix = fit_ringy[-1][1]-fit_ringy[-1][0]
    check = 25*delpix/(yup-ylo)
    print "Check pix size:", check
    
    frac_ring = (YCOO-fit_ringy[-1][0])/delpix
    yy[idx] = (yup-ylo)*frac_ring+ylo
    print 'Fraction of ring:', frac_ring, yy[idx]

# Concerning the ring, i.e. for background subtraction
    bkg = fits.open(key)[0].data
    ber = fits.open(key)[1].data
    yup = int(np.round(yup))
    ylo = int(np.round(ylo))
    bkg_sum[idx, :, :, 0] = bkg[yup-3:yup+4,:]
    bkg_err[idx, :, :, 0] = ber[yup-3:yup+4,:]**2
    bkg_sum[idx, :, :, 1] = bkg[ylo-3:ylo+4,:]
    bkg_err[idx, :, :, 1] = ber[ylo-3:ylo+4,:]**2

# Finalize background
bkg = np.mean(bkg_sum, axis=(0,1,3))
ber = np.sqrt(np.mean(bkg_err, axis=(0,1,3)))

# Finalize y coordinate
print 'Average:', np.mean(yy)
print 'SEM:', np.std(yy)/np.sqrt(yy.shape[0]-1)
yy = int(np.round(np.mean(yy)))



################################################################
# Prepare spectra
DEL_LAM = 4.882 # Difference in Angstrom in x axis, from .fits
REF_LAM = 600 # Reference x pixel (for unit). 0 indexed
REF_VAL = 7751 # Reference value in Angstrom

cc = 299792.458
vv = 287.
lam = np.arange(0, bkg.shape[0])-REF_LAM
lam = REF_VAL+lam*DEL_LAM
# Dereden
correction = cor_red(lam, RV, EBV)
lam = lam*(1-vv/cc) # Do this after assuming ISM is not comoving with ejecta
bkg *= correction
ber *= correction

SPC_EXT = 15 # Spectral extraction height

model = np.loadtxt('/Users/silver/box/phd/pro/87a/lim/misc/SN1987Amodel_8years.txt')
mod_lam = model[:,0]
mod_flx = model[:,1]*np.exp(-28/85.0)/np.exp(-8/85.0)
bkg_wgh = np.array([0.65, 0.55, 0.65, 0.55, 0.65])
avg_wgh = np.mean(bkg_wgh)#3.8

all_spc = np.zeros((len(files), nn))
all_err = np.zeros((len(files), nn))

for idx, key in enumerate(files):
    spc = fits.open(key)[0].data
    err = fits.open(key)[1].data

    # Ejecta spectra
    lo = yy-SPC_EXT/2
    hi = yy+SPC_EXT/2+1
    spc = np.sum(spc[lo:hi,:], axis=0)
    err = np.sqrt(np.sum(err[lo:hi,:]**2, axis=0))#/err[lo:hi,:].shape[0]

    # Deredden
    all_spc[idx] = correction*spc
    all_err[idx] = correction*err
    
    # Plots
#    sub = spc - bkg_wgh[idx]*bkg/np.amax(bkg)*np.amax(spc)
#    plt.plot(lam, sub)
#    plt.plot(lam, bkg_wgh[idx]*bkg/np.amax(bkg)*np.amax(spc))
#    plt.plot(mod_lam, mod_flx)
#    plt.show()
#
#avg_spc = np.sum(all_spc, axis=0)
#avg_err = np.sqrt(np.sum(all_err**2, axis=0))
#avg_sub = avg_spc - avg_wgh*bkg/np.amax(bkg)*np.amax(avg_spc)
#
#plt.errorbar(lam, avg_sub, yerr=avg_err)
#plt.plot(lam, avg_wgh*bkg/np.amax(bkg)*np.amax(avg_spc))
#plt.plot(mod_lam, mod_flx)
#plt.show()

################################################################
# Signal data
def draw_regions():
    for reg in regions:
        plt.gca().axvspan(reg[0], reg[1], alpha=0.2, color='k')
    
SPC_EXT = 4 # Spectral extraction height
regions = [[5450, 5550], [6025, 6100], [6800, 7100], [7500, 7700]] # Spectral regions used for fitting
regions = [[6025, 6100], [6850, 6950], [7550, 7650]] # Spectral regions used for fitting
#regions = [[6025, 6100], [7550, 7650]] # Spectral regions used for fitting
model = np.loadtxt('/Users/silver/box/phd/pro/87a/lim/misc/SN1987Amodel_8years.txt')
mod_lam = model[:,0]
mod_flx = model[:,1]*np.exp(-28/85.0)/np.exp(-8/85.0)
avg_wgh = 0.24

all_spc = np.zeros((len(files[1:3]), nn))
all_err = np.zeros((len(files[1:3]), nn))

for idx, key in enumerate(files[1:3]):
    spc = fits.open(key)[0].data
    err = fits.open(key)[1].data

    # Ejecta spectra
    lo = yy-SPC_EXT/2
    hi = yy+SPC_EXT/2+1
    spc = np.sum(spc[lo:hi,:], axis=0)
    err = np.sqrt(np.sum(err[lo:hi,:]**2, axis=0))

# Sky subtraction
    sky = fits.open(key)[0].data
    ser = fits.open(key)[1].data

## Verify sky
#    sky_slice1 = np.convolve(np.mean(sky[280:290,:], axis=0), np.ones(30), 'same')
#    sky_slice2 = np.convolve(np.mean(sky[110:160,:], axis=0), np.ones(30), 'same')
#    sky_slice3 = np.convolve(np.mean(sky[302:317,:], axis=0), np.ones(30), 'same')
#    plt.plot(lam, sky_slice1)
#    plt.plot(lam, sky_slice2)
#    plt.plot(lam, sky_slice3)
#    plt.show()
#
## More rigorous verification of sky
#    sky_rows = np.r_[110:160, 280:290, 302:317]
#    sky_all  = np.zeros((sky_rows.size, sky.shape[1]))
#    ser_all  = np.zeros((sky_rows.size, sky.shape[1]))
#    for row_i, row in enumerate(sky_rows):
#        sky_all[row_i, :] = sky[row,:]
#        ser_all[row_i, :] = ser[row,:]
#
#    sky_all = np.std (sky_all, axis=0)
#    ser_all = np.mean(ser_all, axis=0)
#    plt.plot(sky_all/ser_all)
#    plt.show()
#    pdb.set_trace()
       
    sky = np.sum(sky[280:290,:], axis=0)+np.sum(sky[110:160,:], axis=0)+np.sum(sky[302:317,:], axis=0)
    ser = np.sqrt(np.sum(np.vstack((ser[280:290,:], ser[110:160,:], ser[302:317,:]))**2 , axis=0))
    sky = sky/75.*(SPC_EXT+1-(avg_wgh/len(files[1:3])))
    ser = ser/75.*(SPC_EXT+1-(avg_wgh/len(files[1:3])))

#    plt.errorbar(lam, sky, yerr=ser)
#    plt.plot(lam, spc)
#    plt.plot(lam, spc-sky)
#    plt.show()

    err = np.sqrt(ser**2+err**2)
    spc = spc-sky
#    plt.plot(ser)
#    plt.plot(err)
#    plt.plot(np.sqrt(ser**2+err**2))
#    plt.show()

# Deredden
    all_spc[idx] = correction*spc
    all_err[idx] = correction*err

avg_spc = np.sum(all_spc, axis=0)
avg_err = np.sqrt(np.sum(all_err**2, axis=0))
wgh_bkg = avg_wgh*bkg
avg_sub = avg_spc - wgh_bkg

########
# Fit to the spectrum
def fit_pow():
    def help_func(xx, CC, alpha):
        return CC*xx**alpha
        
    line_free = [] # Index of line free regions
    for reg in regions:
        line_free.append((lam > reg[0]) & (lam < reg[1]))
    line_free = np.logical_or.reduce(np.array(line_free))

    fitx = lam[line_free]
    fity = avg_sub[line_free]
    fite = avg_err[line_free]
    guess = (1e-17, 0)
    params, covar = curve_fit(help_func, fitx, fity, guess, sigma=fite)

    chi2 = np.sum((fity-params[0]*fitx**params[1])**2/fite**2)
    print "Power law parameters:", params
    print 'chi2:', chi2, '\tdof:', fitx.size, 'red:', chi2/fitx.size

    # Push to 3 sigma upper limit
    steppar_chi2 = chi2
    global steppar_CC
    steppar_CC = params[0]
    print "Test increase:"
    while steppar_chi2 < chi2+7.74049755:
        steppar_CC *= 1.000001
        steppar_chi2 = np.sum((fity-steppar_CC*fitx**params[1])**2/fite**2)
    print steppar_CC, params[0], steppar_CC/params[0], steppar_chi2-chi2
        
    polycoef = np.polyfit(fitx, fity, 3)
    chi2 = np.sum((fity-np.polyval(polycoef, fitx))**2/fite**2)
    print 'chi2:', chi2, '\tdof:', fitx.size, 'red:', chi2/fitx.size    
    return params
    
params = fit_pow()

########
# Print data
lam_min = 5300.
lam_max = 10000.
obs_flx =steppar_CC/(params[1]+1)*(lam_max**(params[1]+1)-lam_min**(params[1]+1))
obs_flx_den = obs_flx/(lam_max-lam_min)
tot_lum = (distance*1e3*pc*100)**2*np.pi*4*obs_flx
tot_lum_den = tot_lum/(lam_max-lam_min)
print 'Observed flux:\t', obs_flx
print 'Spectral flux density:\t', obs_flx_den
print 'Luminosity:\t', tot_lum
print 'Spectral luminosity density:\t', tot_lum_den
print 'Solar luminosities:\t', tot_lum/Lsun
print 'Adjusted params:\t', steppar_CC*1e4**params[1]

########
# Plain plot
show = (mod_lam > 4800) & (mod_lam < 10500)
plt.errorbar(lam, avg_sub, yerr=avg_err)
plt.plot(lam, wgh_bkg, '#1f77b4')
plt.plot(mod_lam[show], 0.1*mod_flx[show], '#ff7f0e')
plt.plot(lam, params[0]*lam**params[1], 'k--')
plt.plot(lam, steppar_CC*lam**params[1], 'k')
draw_regions()
plt.show()

########
# Smooth the spectrum
def plot_smoothed():
    scale = 1e19
    fig = plt.figure(figsize=(5, 3.75))
    kernel = np.array([0.25, 0.5, 0.75, 1., 0.75, 0.5, 0.25])
    kernel = np.ones(11)
    kernel /= np.sum(kernel)
    smt_sub = np.convolve(avg_sub, kernel, 'same')

    plt.plot(mod_lam[show], scale*np.amax(smt_sub)*mod_flx[show]/np.amax(mod_flx[show]), '#ff7f0e', alpha=1)
    plt.plot(lam, scale*smt_sub, '#1f77b4', alpha=1)
    plt.plot(lam, scale*params[0]*lam**params[1], 'k--', lw=3)
    plt.plot(lam, scale*steppar_CC*lam**params[1], 'k', lw=3)
    plt.xlabel('Wavelength (\AA)')
    plt.ylabel('Flux density ($10^{-19}$~erg\ s$^{-1}$\ cm$^{-2}$\ \AA$^{-1}$)')
    plt.xlim([lam_min, lam_max])
    plt.ylim([-10, 150])
    draw_regions()
    fig.savefig('/Users/silver/box/phd/pro/87a/lim/art/figs/spc_lim.pdf', bbox_inches='tight', pad_inches=0.03)
    plt.show()

plot_smoothed()

################################################################
# Compare models
model2 = np.loadtxt('/Users/silver/box/phd/pro/87a/lim/misc/SN1987A_20y_NIR.dat.txt')
mod2_lam = model2[:,0]
mod2_flx = model2[:,1]*np.exp(-28/85.0)/np.exp(-20/85.0)
#plt.plot(mod_lam, mod_flx)
#plt.plot(mod2_lam, mod2_flx)
#plt.show()

pdb.set_trace()
