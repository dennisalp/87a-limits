# Dennis Alp 2017-04-13
# Fit stuff to SN 1987A remnant and determine center, VLT SINFONI
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python -i /Users/silver/box/bin/phd_87a_uplim_find_coords_snf.py

from datetime import date

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
from scipy.stats import sem
from glob import glob

from phd_87a_red import cor_red

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rcParams.update({'figure.max_open_warning': 0})

#########
# Parameters
DAT_DIR = "/Users/silver/dat/vlt/87a/sin/"
WRK_DIR = "/Users/silver/dat/vlt/87a/sin/"
RINGEST = "/Users/silver/dat/vlt/87a/sin/kcube.fits"
INTERPOLATION = 'nearest'

VKICK = 0.8e6
distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e16
Lsun = 3.826e33 # erg s-1

# Reddening
RV = 3.1
EBV = 0.19

# Systematic velocity
cc = 299792.458
vv = 287.

# ADD 10% BECAUSE OF FLUX CALIBRATION UNCERTAINTY?
# SPYROMILIO SAYS NO
# "Almost always, the case of too god a chi-squared fit is that the experimenter,
# in a "fit" of conservativism, has overestimated his or her errors."
flx_cal_err = 1.

# For line removal
ejecta_v = 2500.
# Wavelengths from Kjaer+2007 and Fransson+2016 (last few H2 lines)
kjaer07_lam = np.array([14982.5, 15098.6, 15146.1, 15203.7, 15270.5, 15349.4, 15451.2, 15570.7, 15714.2, 15894.7, 16009.8, 16123.2, 16423.6, 16443.8, 16451.8, 16653.4, 16785.1, 16822.7, 17017.1, 17128.9, 17376.9, 17466.1, 20085.5, 20169.6, 20477.9, 20600.3, 21141.2, 21346.2, 21633.1, 21675.4, 21905.2, 22065.9, 22202.6, 22259.3, 22458.0, 22552.5, 23103.0, 23238.1, 15000.0, 16200.0, 17300.0, 17500.0, 19530.0, 20340.0, 20710.0, 21220.0, 22230.0, 22480.0, 23510.0, 24000.0])
# Fluxes from Kjaer+2007 and Fransson+2016 (last few H2 lines)
kjaer07_flx = np.sqrt(np.array([2, 5, 3, 3, 4, 19, 5, 5, 7, 10, 11, 12, 12, 98, 98, 6, 12, 15, 15, 3, 21, 6, 5, 8, 8, 108, 8, 3, 4, 58, 3, 2, 2, 4, 5, 4, 3, 3, 5, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]))*13.6 # Make max value correspond to 5000 km s-1 in units of pixels

os.chdir(WRK_DIR) #Move to designated directory
files = sorted(glob(DAT_DIR+'/sf_*.fits')) #Find all files.
reg_paths = sorted(glob(DAT_DIR+'?_bad_pix.txt'))
NN = len(files)
coords = []
all_lam = []
all_spc = []
all_pix = []
all_fre = []
all_sky = []
bkg_wgh = [1/300., 1/750.] # Background weights
bkg_wgh = [3.75, 1.5] # Background weights
boxes = np.array([
    [27, 56, 10, 50],
    [31, 62, 16, 56]])

# Sky background boxes
#Bin 0
#Flux 1.41806226196e-16 Specific flux density: 2.18163424917e-19
#chi2: 70.9248581754 . dof: 61 . reduced: 1.16270259304
#Bin 1
#Flux 1.90538602131e-16 Specific flux density: 1.43802718589e-19
#chi2: 359.15033458 . dof: 310 . reduced: 1.15854946639
#Bin 2
#Flux 1.51942369433e-17 Specific flux density: 5.52517707029e-20
#chi2: 93.6585329869 . dof: 89 . reduced: 1.05234306727
#Bin 3
#Flux 5.17982628151e-16 Specific flux density: 1.81748290579e-19
#chi2: 425.9652543 . dof: 385 . reduced: 1.10640325792
#Bin 4
#Flux 1.0791303431e-16 Specific flux density: 9.8102758464e-20
#chi2: 337.64083544 . dof: 336 . reduced: 1.00488343881
#
#Band 0
#Linear flux: 5.53653645137e-16 . Specific flux density: 1.70354967734e-19
#chi2: 525.408043163 . dof: 460 . reduced: 1.14219139818
#Quadratic flux: 5.50262934626e-16 . Specific flux density: 1.69311672193e-19
#chi2: 519.114270741 . dof: 460 . reduced: 1.12850928422
#Power law flux: 5.46977309473e-16 . Specific flux density: 1.68300710607e-19
#chi2: 538.593307467 . dof: 460 . reduced: 1.17085501623
#Band 1
#Linear flux: 6.5991860549e-16 . Specific flux density: 1.6706800139e-19
#chi2: 805.474127918 . dof: 721 . reduced: 1.11716245204
#Quadratic flux: 6.22950847915e-16 . Specific flux density: 1.57709075422e-19
#chi2: 788.585503713 . dof: 721 . reduced: 1.09373856271
#Power law flux: 6.7581968872e-16 . Specific flux density: 1.71093592081e-19
#chi2: 809.54456406 . dof: 721 . reduced: 1.12280799453
skies = np.array([[
    [ 9, 13, 17, 43],
    [54, 58, 54, 63]],
     [[59, 63, 59, 63],
    [19, 25, 22, 62]]])
sky_label = 'sky1'

#Bin 0
#Flux 1.3064742216e-16 Specific flux density: 2.00996034092e-19
#chi2: 74.1423387413 . dof: 61 . reduced: 1.21544817609
#Bin 1
#Flux 1.75482170422e-16 Specific flux density: 1.32439373903e-19
#chi2: 362.049845138 . dof: 310 . reduced: 1.16790272625
#Bin 2
#Flux 1.43946879568e-17 Specific flux density: 5.23443198429e-20
#chi2: 93.909686606 . dof: 89 . reduced: 1.05516501804
#Bin 3
#Flux 4.79066049551e-16 Specific flux density: 1.6809335072e-19
#chi2: 438.691853109 . dof: 385 . reduced: 1.13945935872
#Bin 4
#Flux 9.87631039202e-17 Specific flux density: 8.97846399275e-20
#chi2: 337.887783417 . dof: 336 . reduced: 1.00561840303
#
#Band 0
#Linear flux: 5.10121510662e-16 . Specific flux density: 1.56960464819e-19
#chi2: 532.844864913 . dof: 460 . reduced: 1.15835840198
#Quadratic flux: 5.07445331295e-16 . Specific flux density: 1.56137025014e-19
#chi2: 529.832646643 . dof: 460 . reduced: 1.1518101014
#Power law flux: 5.11546801022e-16 . Specific flux density: 1.57399015699e-19
#chi2: 541.149911484 . dof: 460 . reduced: 1.17641285105
#Band 1
#Linear flux: 6.10802692419e-16 . Specific flux density: 1.54633593017e-19
#chi2: 804.673333729 . dof: 721 . reduced: 1.11605178048
#Quadratic flux: 5.79163023111e-16 . Specific flux density: 1.46623550155e-19
#chi2: 794.370742619 . dof: 721 . reduced: 1.10176247243
#Power law flux: 6.25698770355e-16 . Specific flux density: 1.58404751989e-19
#chi2: 808.705927985 . dof: 721 . reduced: 1.1216448377
#skies = np.array([[
#    [ 14, 18, 17, 43],
#    [59, 63, 54, 63]],
#     [[54, 58, 59, 63],
#    [19, 25, 22, 32]]])
#sky_label = 'sky2'

#Bin 0
#Flux 1.43574238729e-16 Specific flux density: 2.20883444198e-19
#chi2: 70.7459253282 . dof: 61 . reduced: 1.15976926768
#Bin 1
#Flux 1.89699317632e-16 Specific flux density: 1.43169296326e-19
#chi2: 355.751809459 . dof: 310 . reduced: 1.14758648213
#Bin 2
#Flux 1.36753919736e-17 Specific flux density: 4.9728698086e-20
#chi2: 94.2859939458 . dof: 89 . reduced: 1.0593931904
#Bin 3
#Flux 4.98387770969e-16 Specific flux density: 1.74872902094e-19
#chi2: 414.884720698 . dof: 385 . reduced: 1.07762265116
#Bin 4
#Flux 9.34723928474e-17 Specific flux density: 8.49749025885e-20
#chi2: 342.061335826 . dof: 336 . reduced: 1.01803968996
#
#Band 0
#Linear flux: 5.53812246373e-16 . Specific flux density: 1.70403768115e-19
#chi2: 533.382940678 . dof: 460 . reduced: 1.15952813191
#Quadratic flux: 5.50639459726e-16 . Specific flux density: 1.6942752607e-19
#chi2: 524.408910012 . dof: 460 . reduced: 1.14001936959
#Power law flux: 5.4595360455e-16 . Specific flux density: 1.67985724477e-19
#chi2: 548.621467356 . dof: 460 . reduced: 1.19265536382
#Band 1
#Linear flux: 6.27716431807e-16 . Specific flux density: 1.58915552356e-19
#chi2: 802.868052957 . dof: 721 . reduced: 1.11354792366
#Quadratic flux: 5.75834664501e-16 . Specific flux density: 1.45780927722e-19
#chi2: 779.460762464 . dof: 721 . reduced: 1.0810828883
#Power law flux: 6.44316947675e-16 . Specific flux density: 1.63118214601e-19
#chi2: 805.446454352 . dof: 721 . reduced: 1.11712406984
#skies = np.array([[
#    [ 9, 13, 43, 58],
#    [57, 62, 5, 10]],
#     [[54, 58, 9, 14],
#    [19, 25, 59, 63]]])
#sky_label = 'sky3'

#########
# Ellipse parameters http://math.stackexchange.com/questions/426150/what-is-the-general-\
# equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate
X0 = np.array([27.7018, 34.6286])-boxes[:,2]
Y0 = np.array([41.6443, 46.3902])-boxes[:,0]
AA = 14
BB = 10
ALPHA = 0.15
subnx = boxes[:,3]-boxes[:,2]
subny = boxes[:,1]-boxes[:,0]
 
def gen_mask(AA, BB, X0, Y0, ALPHA):
    xx = np.arange(0,NX)
    yy = np.arange(0,NY)[:,None]
    xx = xx-X0
    yy = yy-Y0
    uu = xx*np.cos(ALPHA)-yy*np.sin(ALPHA)
    vv = xx*np.sin(ALPHA)+yy*np.cos(ALPHA)
    return ((uu/AA)**2 + (vv/BB)**2 > 0.5) & ((uu/AA)**2 + (vv/BB)**2 < 2.7) # True for points inside

# Just reasonable plot of sky images, good for debugging
def sky_plt(image):
    plt.figure()
    plt.imshow(image,interpolation='nearest', cmap='afmhot', origin='lower')
#    plt.show()

# Days for scaling
SNDATE = date(1987, 2, 23)
def get_days(fname): # 1994-09-24_drz_al2.fits
    yr = int(fname[-23:-19])
    month = int(fname[-18:-16].lstrip("0"))
    day = int(fname[-15:-13].lstrip("0"))
    return (date(yr, month, day)-SNDATE).days*24*60*60

def inside_kick(yy, xx, key, YCOO, XCOO):
    telapsed = get_days(key)
    offset = np.rad2deg(np.arctan(VKICK*telapsed/(distance*1e3*pc)))*60*60*1000/50
    dely = yy-YCOO
    delx = xx-XCOO
    rr = np.sqrt(delx**2+dely**2)
    if rr > offset+1/np.sqrt(2): # This selects fractional pixels inside kick
        return 0.
    elif rr < offset-1/np.sqrt(2):
        return 1.
    else:
        xx += np.random.rand(1000000)-0.5
        yy += np.random.rand(1000000)-0.5
        dely = yy-YCOO
        delx = xx-XCOO
        rr = np.sqrt(delx**2+dely**2)
        rr = rr < offset
        return np.sum(rr)/1000000.

def fix_unit(img, img_path):
    if 'sf_h' in img_path:
        return img*1e-17/600.
    if 'sf_k' in img_path:
        eff = np.loadtxt('efficiency.txt')[:,1]
        return img*eff/600.
    
#########
# Fit ellipses for all observations
def help_ell(dummy, x0, y0, aa, bb, alpha, mag, sig, tilt, phi):
#     aa, bb, y0, x0, alpha, mag, sig = AA, BB, Y0, X0, ALPHA, 1., 1.
    xx = np.arange(0,subnx[ii])
    yy = np.arange(0,subny[ii])[:,None]
    xx = xx-x0
    yy = yy-y0
    uu = xx*np.cos(alpha)-yy*np.sin(alpha)
    vv = xx*np.sin(alpha)+yy*np.cos(alpha)
    rr = (uu/aa)**2 + (vv/bb)**2
#    print x0, y0, aa, bb, alpha, mag, sig, tilt, phi
    global last_ell
    last_ell = abs(mag)*np.exp(-np.abs((np.sqrt(rr)-1)/sig)**2)*(1+tilt*np.cos(np.arctan2(vv,uu)+phi))
    return np.reshape(last_ell, subnx[ii]*subny[ii])

def draw_regions(regions):
    tot_out = 0
    for ii in range(0,regions.shape[0]):
        plt.gca().axvspan(regions[ii,0]-0.5, regions[ii,1]-0.5, alpha=0.2, color='k')
        tot_out += regions[ii,1]-regions[ii,0]

    print 'Total outside of regions:', tot_out, 'Inside', regions[ii,1]-tot_out, 'ratio in', 1-tot_out/regions[ii,1]

def get_free(regions, spc):
    idx = np.arange(0,spc.size)
    line_free = [] # Index of line free regions
    for ii in range(0,regions.shape[0]):
        line_free.append((idx >= regions[ii,0]) & (idx < regions[ii,1]))
    idx = np.logical_and.reduce(-np.array(line_free))
    return idx

def help_func(xx, AA, CC, alpha):
    chi2 = AA*spc_ring[line_free]+CC*xx**alpha
    chi2 = (chi2-fity)/fite
    chi2 = np.sum(chi2**2)
    print chi2, AA, CC, alpha
    return AA*spc_ring[line_free]+CC*xx**alpha

def show_func(xx, AA, CC, alpha):
    plt.plot(fitx, fity)
    plt.plot(fitx, AA*spc_ring[line_free]+CC*xx**alpha)
    plt.plot(fitx, AA*spc_ring[line_free])
    plt.plot(fitx, CC*xx**alpha)
#    plt.show()
        
################################################################
# Main part
for ii, img_path in enumerate(files):
# Load image box
    img = fits.open(img_path)[0].data
    img = np.where(np.isnan(img), 0, img)
    img = np.rollaxis(img, 2)
    img = np.rollaxis(img, 2) # Puts spectral info last

# Get wavelengths
    REF_LAM = fits.getval(img_path, 'CRPIX3')-1
    REF_VAL = fits.getval(img_path, 'CRVAL3')*1e4
    DEL_LAM = fits.getval(img_path, 'CDELT3')*1e4
    lam = np.arange(0, img.shape[2])-REF_LAM
    lam = REF_VAL+lam*DEL_LAM

# Deredden
    correction = cor_red(lam, RV, EBV)
    lam = lam*(1-vv/cc) # Do this after assuming ISM is not comoving with ejecta
    img *= correction

# Preserve sky 
    img = fix_unit(img, img_path)
    sky_plt(np.sum(img, axis=2))
    for sky_box in range(0, skies.shape[1]):
        plt.gca().add_patch(patches.Rectangle((skies[ii,sky_box,2]-0.5, skies[ii,sky_box,0]-0.5), skies[ii,sky_box,3]-skies[ii,sky_box,2], skies[ii,sky_box,1]-skies[ii,sky_box,0], fill=False, color='g'))
    
    sky_area = (skies[ii,0,1]-skies[ii,0,0])*(skies[ii,0,3]-skies[ii,0,2])
    sky_area += (skies[ii,1,1]-skies[ii,1,0])*(skies[ii,1,3]-skies[ii,1,2])
    all_sky.append(np.sum(img[skies[ii,0,0]:skies[ii,0,1],skies[ii,0,2]:skies[ii,0,3]], axis=(0, 1)))
    all_sky[ii] += np.sum(img[skies[ii,1,0]:skies[ii,1,1],skies[ii,1,2]:skies[ii,1,3]], axis=(0, 1))
    all_sky[ii] /= sky_area # Because this must be scaled by extraction region
    np.save(sky_label + '_band' + str(ii), np.array(all_sky[ii]))
    img = img[boxes[ii,0]:boxes[ii,1],boxes[ii,2]:boxes[ii,3]]

# Define some stuff
    NX = img.shape[1]
    NY = img.shape[0]
    mask = gen_mask(AA, BB, X0[ii], Y0[ii], ALPHA)

    avg = np.sum(img, axis=2)
    msk = np.where(mask, avg, 0.)
    volume = np.sum(np.abs(msk))

# Merge selection filters
#    for idx, line in enumerate(kjaer07_lam):
#        print (line-REF_VAL)/DEL_LAM+REF_LAM, '\t', (line-REF_VAL)/DEL_LAM+REF_LAM-kjaer07_flx[idx], '\t', (line-REF_VAL)/DEL_LAM+REF_LAM+kjaer07_flx[idx]
        
# Get the ring spectrum for fitting
    mask3d = np.repeat(mask[:,:,np.newaxis], img.shape[2], axis=2)
    rin = np.where(mask3d, img, 0.)
    rin = np.mean(rin, axis=(0,1))
    
# Fit to get coordinates
    guess = np.array([X0[ii], Y0[ii], AA, BB, ALPHA, 1e-16, 3, 0.25, 0.])
    pars, covar = curve_fit(help_ell, np.arange(0, subnx[ii]*subny[ii]), np.reshape(msk, subnx[ii]*subny[ii]), guess, maxfev=20000)
    perr = np.sqrt(np.diag(covar))
    print 'Positional accuracy:', np.sqrt(perr[0]**2+perr[1]**2)*50, 'mas' 
    pars[0] = pars[0] + boxes[ii,2]
    pars[1] = pars[1] + boxes[ii,0]
    res = np.sum(np.abs(msk-last_ell)) / volume
    print("{0:17.13f} {1:17.13f} {2:17.13f}".format(pars[0], pars[1], res))
    coords.append([pars[0], pars[1]])
    
# Filter wavelengths/lines
    spc = np.sum(img, axis=(0,1))
    std = np.std(img, axis=(0,1))
    regions = np.loadtxt(reg_paths[ii])

    plt.figure()
    draw_regions(regions)
    plt.plot(spc)
    plt.plot(1000*std)
#    plt.show() # Reveals all bad pixels
    line_free = get_free(regions, spc)
    np.save('line_free_band' + str(ii), line_free)
    pix_free = np.arange(0, spc.size)[line_free]

    plt.figure()
    draw_regions(regions)
    plt.plot(pix_free , spc[line_free], '.')
#    plt.show()
    cont_map = np.sum(img[:,:,line_free],axis=2)
    
# Get spectrum
    all_fre.append(line_free)
    all_pix.append(pix_free)
    all_lam.append(lam)
    all_spc.append(np.zeros(spc.size))
    ycent = coords[ii][1]-boxes[ii,0]
    xcent = coords[ii][0]-boxes[ii,2]
    for yy in range(0, img.shape[0]):
        for xx in range(0, img.shape[1]):
            wgh = inside_kick(yy, xx, img_path, ycent, xcent)
            all_spc[ii] = all_spc[ii] + wgh*img[yy,xx]

# Subtract spread light
    telapsed = get_days(img_path)
    offset = np.rad2deg(np.arctan(VKICK*telapsed/(distance*1e3*pc)))*60*60*1000/50

    plt.figure()
    plt.plot(all_spc[ii])
    all_spc[ii] -= bkg_wgh[ii]*rin # This removes spread light
    all_spc[ii] -= (np.pi*offset**2-bkg_wgh[ii])*all_sky[ii] # The parenthesis is the search area but with bkg_wgh removed because that is subtracted in the line above 
#    plt.plot(all_spc[ii])
    plt.plot(bkg_wgh[ii]*rin)
#    plt.plot(np.pi*offset**2*all_sky[ii])
    plt.show()
    
# Plots
    plt.figure()
    plt.imshow(msk, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(msk)/100, vmax=np.amax(msk)))
#    plt.show()

    plt.figure()
    plt.imshow(msk-last_ell, cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
#    plt.show()

    plt.figure()
    plt.plot(np.reshape(last_ell, subnx[ii]*subny[ii]))
    plt.plot(np.reshape(msk, subnx[ii]*subny[ii]))
#    plt.show()

#########
# Plot and wrap-up
    plt.figure()
    plt.plot(coords[ii][0]-boxes[ii,2], coords[ii][1]-boxes[ii,0], 'ok')
#    plt.imshow(cont_map, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.max((np.amax(cont_map)/30,1.*np.amin(cont_map))), vmax=np.amax(cont_map)))
    plt.imshow(cont_map, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amin(cont_map), vmax=np.amax(cont_map)))
#    plt.show()

#########
## Correlation map
## This probably wont work, no signal will always be no signal
## Get the ring spectrum for fitting
#    mask3d = np.repeat(mask[:,:,np.newaxis], img.shape[2], axis=2)
#    spc_ring = np.where(mask3d, img, 0.)
#    spc_ring = np.sum(spc_ring, axis=(0,1))
#
#    for yy in range(0, img.shape[0]):
#        for xx in range(0, img.shape[1]):
#            fitx = lam[line_free]
#            fity = img[yy, xx, line_free]
#            fite = std[line_free]
#            guess = (0.0013, 1e8, -2)
#            params, covar = curve_fit(help_func, fitx, fity, guess, sigma=fite, maxfev=10000)
#            print params
#            show_func(fitx, params[0], params[1], params[2])
#            pdb.set_trace()

################################################################
# Load models
model = np.loadtxt('/Users/silver/box/phd/pro/87a/lim/misc/SN1987Amodel_8years.txt')
mod_lam = model[:,0]
mod_flx = model[:,1]*np.exp(-28/85.0)/np.exp(-8/85.0)
model2 = np.loadtxt('/Users/silver/box/phd/pro/87a/lim/misc/SN1987A_20y_NIR.dat.txt')
mod2_lam = model2[:,0]
mod2_flx = model2[:,1]*np.exp(-28/85.0)/np.exp(-20/85.0)

################################################################
# Plot final and get limits
def help_func(xx, CC, alpha):
    return CC*xx**alpha

plt.figure()
for ii in range(0, len(all_spc)):
    # Fit power law
    fitx = all_lam[ii][all_fre[ii]]
    fity = all_spc[ii][all_fre[ii]]
    fite = np.zeros(fitx.size)
    for idx, yval in enumerate(fity):
        run_y = fity[np.max((idx-4,0)):idx+5]
        run_avg = np.mean(run_y)
        run_n = run_y.size
        fite[idx] = np.sqrt(np.sum((run_y-run_avg)**2)/run_n)
        
    guess = (1e-10, -2.)
    params, covar = curve_fit(help_func, fitx, fity, guess, maxfev=20000, sigma=fite)
#    params = np.polyfit(fitx, fity, 16)
    print "Power law parameters:", params
    
    # Plot
    plt.errorbar(fitx, fity, fmt='bo', yerr=fite, alpha=0.2)

    plt.plot(fitx, params[0]*fitx**params[1], 'k')
    chi2 = np.sum((params[0]*fitx**params[1]-fity)**2/fite**2)

#    plt.plot(fitx, np.polyval(params, fitx), 'k')
#    chi2 = np.sum((np.polyval(params, fitx)-fity)**2/fite**2)
    
    print 'chi2:', chi2, '\tdof:', fitx.size, 'red:', chi2/fitx.size
    
all_lam = np.concatenate(all_lam)
all_spc = np.concatenate(all_spc)
all_pix = np.concatenate(all_pix)
all_fre = np.concatenate(all_fre)

#fitx = all_lam[all_fre]
#fity = all_spc[all_fre]
#params, covar = curve_fit(help_func, fitx, fity, guess, maxfev=20000)
#print "Power law parameters:", params
#plt.plot(fitx, params[0]*fitx**params[1], 'k')

plt.figure()
# Just shows all Kjaer+2007 lines (and H2 from Fransson+2016)
for idx, line in enumerate(kjaer07_lam):
    plt.gca().axvspan((1-ejecta_v/cc)*line, (1+ejecta_v/cc)*line, alpha=0.5, color='k')

# Plot with smooth
kernel = np.array([0.25, 0.5, 0.75, 1., 0.75, 0.5, 0.25])
kernel = np.ones(8)
kernel /= np.sum(kernel)
smt_spc = np.convolve(all_spc, kernel, 'same')
smt_lam = np.convolve(all_lam, kernel, 'same')
plt.plot(all_lam, smt_spc)
plt.plot(mod_lam, mod_flx/30)
plt.plot(mod2_lam, mod2_flx/10)

plt.figure()
kernel = np.array([0.25, 0.5, 0.75, 1., 0.75, 0.5, 0.25])
kernel = np.ones(20)
kernel /= np.sum(kernel)
smt_spc = np.convolve(all_spc[all_fre], kernel, 'same')
plt.plot(all_lam[all_fre], smt_spc, '.')
#plt.plot(all_lam[all_fre], all_spc[all_fre], '.')



################################################################
# Some experimental stuff on background properties
lfh = np.load('line_free_band0.npy')
skyh = np.zeros((3, lfh.size))
skyh[0, :] = np.pi*offset**2*np.load('sky1_band0.npy')
skyh[1, :] = np.pi*offset**2*np.load('sky2_band0.npy')
skyh[2, :] = np.pi*offset**2*np.load('sky3_band0.npy')

lfk = np.load('line_free_band1.npy')
skyk = np.zeros((3, lfk.size))
skyk[0, :] = np.pi*offset**2*np.load('sky1_band1.npy')
skyk[1, :] = np.pi*offset**2*np.load('sky2_band1.npy')
skyk[2, :] = np.pi*offset**2*np.load('sky3_band1.npy')

sigh = np.std(skyh[:, lfh], axis=0)
sigh = np.convolve(sigh, np.ones(13)/13., 'same')
sigk = np.std(skyk[:, lfk], axis=0)
sigk = np.convolve(sigk, np.ones(13)/13., 'same')
sig_bkg = np.concatenate((sigh, sigk))

# Mean in each slice
diff = np.diff(all_lam[all_fre])
cuts = np.where((diff[0:-1] > 1.1*DEL_LAM)==True)[0]+1
cuts = np.insert(cuts,0,0) # Inserts 0 at index 0
cuts = np.append(cuts,all_lam[all_fre].size)

cur_lam = np.zeros(cuts.size-1)
cur_spc = np.zeros(cuts.size-1)
cur_err = np.zeros(cuts.size-1)
cur_xer = np.zeros(cuts.size-1)
all_err = np.zeros(all_lam.size)

for i in range(0,len(cuts)-1):
    tmp_lam = all_lam[all_fre][cuts[i]:cuts[i+1]]
    tmp_spc = all_spc[all_fre][cuts[i]:cuts[i+1]]
    cur_err[i] = sem(tmp_spc)
    cur_xer[i] = (np.amax(tmp_lam)-np.amin(tmp_lam))/2.
    cur_lam[i] = np.mean(tmp_lam)
    cur_spc[i] = np.mean(tmp_spc)

    # Also define the error for each point
    cut_idx = np.where(all_fre)[0][cuts[i]:cuts[i+1]]
    all_err[cut_idx] = np.std(tmp_spc)

plt.figure()
#plt.plot(cur_lam, cur_spc)
plt.errorbar(cur_lam, cur_spc, fmt='bo', alpha=0.4, yerr=cur_err, xerr=cur_xer)
#pdb.set_trace()

# INTERVALS ARE DEFINED HERE
intervals = np.array([[15150, 15800], [16900, 18125], [18125, 18400], [19875,22725], [22725, 23825]])
intervals = np.array([[15150, 15800], [16900, 18125], [19875,22725], [22725, 23825]])
def draw_intervals():
    for ii in range(0,intervals.shape[0]):
        plt.gca().axvspan(intervals[ii,0], intervals[ii,1], alpha=0.2, color='k')

plt.figure()
draw_intervals()
plt.errorbar(cur_lam, cur_spc, fmt='bo', alpha=0.4, yerr=cur_err, xerr=cur_xer)

########
# Make binned limits
print ''
plt.close('all')

scale = 1e19
fig = plt.figure(figsize=(5, 3.75))
# Plot the lines
#plt.errorbar(cur_lam, cur_spc, fmt='ok', alpha=0.6, yerr=cur_err, xerr=cur_xer, markersize=0, elinewidth=2, capsize=0)
tmp_spc = np.convolve(all_spc, np.ones(21)/21., 'same')
tmp_in = all_lam < 18500
plt.plot(all_lam[tmp_in], tmp_spc[tmp_in]*scale*flx_cal_err, '#1f77b4', alpha=1)
tmp_in = all_lam > 19500
plt.plot(all_lam[tmp_in], tmp_spc[tmp_in]*scale*flx_cal_err, '#1f77b4', alpha=1)
model2 = np.loadtxt('/Users/silver/box/phd/pro/87a/lim/misc/SN1987A_20y_NIR.dat.txt')
mod2_lam = model2[:,0]
mod2_flx = model2[:,1]*0.046
#plt.plot(mod2_lam, mod2_flx*scale*flx_cal_err, '#ff7f0e', alpha=1, zorder=-1)

bin_lam = np.zeros(intervals.shape[0])
bin_spc = np.zeros(intervals.shape[0])
bin_err = np.zeros(intervals.shape[0])
fin_lim = np.zeros(intervals.shape[0])
for ii in range(0,intervals.shape[0]):
    print 'Bin', ii
    lo = intervals[ii,0]
    up = intervals[ii,1]

    # Bins
    segment = (all_lam[all_fre] > lo) & (all_lam[all_fre] < up)
    tmp_lam = all_lam[all_fre][segment]
    tmp_spc = all_spc[all_fre][segment]
    tmp_err = all_err[all_fre][segment]
    bin_lam[ii] = np.mean(tmp_lam)
    bin_spc[ii] = np.mean(tmp_spc)
    bin_err[ii] = sem(tmp_spc)

    chi2 = np.sum((bin_spc[ii]-tmp_spc)**2/tmp_err**2)
    print 'Flux', bin_spc[ii]*(up-lo), 'Specific flux density:', bin_spc[ii]
    print 'chi2:', chi2, '. dof:', tmp_spc.size, '. reduced:', chi2/tmp_spc.size

    # Push to 3 sigma upper limit
    steppar_chi2 = chi2
    steppar_lim = bin_spc[ii]
    while steppar_chi2 < chi2+7.74049755:
        steppar_lim *= 1.000001
        steppar_chi2 = np.sum((tmp_spc-steppar_lim)**2/tmp_err**2)
    print 'Flux', steppar_lim*(up-lo), 'Specific flux density:', steppar_lim, steppar_chi2-chi2

    steppar_lim *= flx_cal_err
    fin_lim[ii] = steppar_lim
    tmp_yy = bin_spc[ii]*flx_cal_err*np.ones(2)
    steppar_lim *= np.ones(2)
    
    # Plot
    plt.plot([lo,up], tmp_yy*scale, 'k--', lw=3)
    plt.plot([lo,up], steppar_lim*scale, 'k', lw=3)
    tmp_xx = bin_lam[ii]*np.ones(2)
    tmp_yy = [tmp_yy-bin_err[ii], tmp_yy+bin_err[ii]]
#    plt.plot(tmp_xx, tmp_yy, 'k', lw=3)

# Limit H and K Band
in_hh = all_lam[all_fre] < 19000
in_kk = np.logical_not(in_hh)
band_idx =  [in_hh, in_kk]
int_lst = [[0,1],[2,3]]

print ''
for ii, idx in enumerate(band_idx):
    print 'Band', ii
    tmp_lam = all_lam[all_fre][idx]
    tmp_spc = all_spc[all_fre][idx]
    tmp_err = all_err[all_fre][idx]
    lo = intervals[int_lst[ii][0],0]
    hi = intervals[int_lst[ii][1],1]
    val_lam = np.linspace(lo, hi, 256)

    # Linear fit
    lin_coe = np.polyfit(tmp_lam, tmp_spc, 1)
    lin_val = np.polyval(lin_coe, val_lam)
    lin_int = np.polyint(lin_coe, 1)
    lin_int = np.polyval(lin_int, hi)-np.polyval(lin_int, lo)
#    plt.plot(val_lam, lin_val, 'r')
    
    chi2 = np.sum((np.polyval(lin_coe, tmp_lam)-tmp_spc)**2/tmp_err**2)
    print 'Linear flux:', lin_int, '. Specific flux density:', lin_int/(hi-lo)
    print 'chi2:', chi2, '. dof:', tmp_spc.size, '. reduced:', chi2/tmp_spc.size
    
    # Quadratic fit
    qad_coe = np.polyfit(tmp_lam, tmp_spc, 2)
    qad_val = np.polyval(qad_coe, val_lam)
    qad_int = np.polyint(qad_coe, 1)
    qad_int = np.polyval(qad_int, hi)-np.polyval(qad_int, lo)
#    plt.plot(val_lam, qad_val, 'g')

    chi2 = np.sum((np.polyval(qad_coe, tmp_lam)-tmp_spc)**2/tmp_err**2)
    print 'Quadratic flux:', qad_int, '. Specific flux density:', qad_int/(hi-lo)
    print 'chi2:', chi2, '. dof:', tmp_spc.size, '. reduced:', chi2/tmp_spc.size

    # Power law fit
    guess = (1e-10, -2.)
    params, covar = curve_fit(help_func, tmp_lam, tmp_spc, guess, maxfev=20000, sigma=tmp_err)
    pow_val = params[0]*val_lam**params[1]
    pow_int = params[0]/(params[1]+1)*(hi**(params[1]+1)-lo**(params[1]+1))
    pow_int_den = pow_int/(hi-lo)
#    tot_lum = (distance*1e3*pc*100)**2*np.pi*4*obs_flx
#    tot_lum_den = tot_lum/(lam_max-lam_min)
#    plt.plot(val_lam, pow_val, 'b')
    
    chi2 = np.sum((params[0]*tmp_lam**params[1]-tmp_spc)**2/tmp_err**2)
    print 'Power law flux:', pow_int, '. Specific flux density:', pow_int_den
    print 'chi2:', chi2, '. dof:', tmp_spc.size, '. reduced:', chi2/tmp_spc.size

# Plot cosmetics
plt.gca().set_xticks([15000, 17000, 19000, 21000, 23000])
plt.gca().set_xticks([15000, 18000, 21000, 24000])
plt.ylim([-0.3, 5])
plt.xlim([15000, 24000])
plt.xlabel('Wavelength (\AA)')
plt.ylabel('Flux density ($10^{-19}$~erg\ s$^{-1}$\ cm$^{-2}$\ \AA$^{-1}$)')

# Add companion
Rsun = 6.96e10 # cm
hh = 6.626e-27
cc = 2.99792458e10
kB = 1.380658e-16
sb = 5.67051e-5 # Stefan-Boltzmann
TT = 5800
AA = 4*np.pi*Rsun**2
lam = np.linspace(1000e-8, 10e-4, 1000)

companion2 = AA*sb*TT**4
companion2 = companion2/(4*np.pi*(distance*1e3*pc*100)**2)
companion = 2*hh*cc**2/lam**5*1/(np.exp(hh*cc/(lam*kB*TT))-1)
companion = AA*np.pi*companion
companion = companion/(4*np.pi*(distance*1e3*pc*100)**2)
#plt.plot(lam/1e-8, companion*1e-8*1e19, color='#2ca02c', ls='--', zorder=-1000)

fig.savefig('/Users/silver/box/phd/pro/87a/lim/art/figs/nir_lim.pdf', bbox_inches='tight', pad_inches=0.03)

plt.figure()
plt.plot(all_lam, np.convolve(all_spc, np.ones(21)/21., 'same'), 'b', alpha=1)
plt.plot(all_lam, np.convolve(all_spc, np.ones(1), 'same'), 'b', alpha=0.2)
plt.errorbar(cur_lam, cur_spc, fmt='ok', alpha=1, yerr=cur_err, xerr=cur_xer, markersize=0, elinewidth=4, capsize=0)

# INTERVALS ARE HARD CODED HERE ALSO
#h_bol_flx  = (intervals[0,1]-15000)*fin_lim[0]
#h_bol_flx += (intervals[1,0]-intervals[0,1])*np.mean(fin_lim[:2]) # Interpolation in the gap
#h_bol_flx += (intervals[1,1]-intervals[1,0])*fin_lim[1]
#h_bol_flx += (18500-intervals[2,0])*fin_lim[2]
#
#k_bol_flx  = (intervals[3,1]-19500)*fin_lim[3]
#k_bol_flx += (24000-intervals[4,0])*fin_lim[4]
h_bol_flx  = (intervals[0,1]-15000)*fin_lim[0]
h_bol_flx += (intervals[1,0]-intervals[0,1])*np.mean(fin_lim[:2]) # Interpolation in the gap
h_bol_flx += (18500-intervals[1,0])*fin_lim[1]

k_bol_flx  = (intervals[2,1]-19500)*fin_lim[2]
k_bol_flx += (24000-intervals[3,0])*fin_lim[3]

print 'H flux:', h_bol_flx, h_bol_flx/(18500-15000), (distance*1e3*pc*100)**2*np.pi*4*h_bol_flx, (distance*1e3*pc*100)**2*np.pi*4*h_bol_flx/Lsun
print 'K flux:', k_bol_flx, k_bol_flx/(24000-19500), (distance*1e3*pc*100)**2*np.pi*4*k_bol_flx, (distance*1e3*pc*100)**2*np.pi*4*k_bol_flx/Lsun
print 'Final bin lims flx', fin_lim
print 'Final bin lims lum', (distance*1e3*pc*100)**2*np.pi*4*fin_lim*(intervals[:,1]-intervals[:,0])
print 'Final bin lims sun', (distance*1e3*pc*100)**2*np.pi*4*fin_lim*(intervals[:,1]-intervals[:,0])/Lsun
print 'S/B in H band:', np.median(all_spc[all_fre][in_hh]/skyh[0, lfh])
print 'S/B in K band:', np.median(all_spc[all_fre][in_kk]/skyk[0, lfk])
# How to stop
pdb.set_trace()
plt.show()
