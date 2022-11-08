# Dennis Alp 2016-10-13
# Fit stuff to SN 1987A remnant and determine center
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python -i /Users/silver/box/bin/phd_87a_uplim_find_coords.py

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

#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

#########
# Parameters
DAT_DIR = "/Users/silver/dat/hst/87a/all/"
WRK_DIR = "/Users/silver/box/phd/pro/87a/lim/"
RINGEST = "/Users/silver/dat/hst/87a/all/acs_r_2006-12-06_drz_al2.fits"
CHANDRAEST = "/Users/silver/dat/hst/87a/all/w3_r_2015-05-24_drz_al2.fits"
THICKNESS = 15 # Thickness of Gaussian radial profile of ellipse
NX = 1150
NY = 1115
INTERPOLATION = 'nearest'
GAIAX1 = 83.86777124038386
GAIAY1 = 69.26996980653968
GAIAX2 = 83.86512140949159
GAIAY2 = 69.26915041477726
# The above two are the saturated in HST, the below two are fine
GAIAX1 = 83.86846175388055
GAIAY1 = 69.27182546903957
GAIAX2 = 83.86242771177496
GAIAY2 = 69.27017546580065
DR2_GAIAX1 = 83.86846278787391
DR2_GAIAY1 = 69.27182569443815
DR2_GAIAX2 = 83.86242840949296
DR2_GAIAY2 = 69.27017569105270
Xstar1 = 499
Ystar1 = 289
Xstar2 = 806
Ystar2 = 527

os.chdir(WRK_DIR) #Move to designated directory
files = glob(DAT_DIR+'/*.fits') #Find all files.
#pdb.set_trace()
#files = files[4:6]
#files = ["/Users/silver/box/phd/data/hst/87a/all/w2_b_1999-01-07_drz_al2.fits","/Users/silver/box/phd/data/hst/87a/all/w2_b_1999-01-07_drz_al2.fits"]
NN = len(files)
alphs =[]
coords = []
save_stars = []
prod_pars = []
Nspots = 26
hotspot_theta = np.zeros(Nspots)
#                [534,574,509,549],   # Star 1
#                [649,689,649,689],   # Star 2
boxes = np.array([[554,619,549,634], # Ellipse
                [564,604,569,609],   # Gauss. Following three must have same area!
                [269,309,479,519],   # Star 1
                [507,547,786,826],   # Star 2
                [605,612,577,583],   # 1
                [605,609,573,576],   # 2
                [602,606,568,573],   # 3
                [596,602,561,567],   # 4
                [590,595,559,564],   # 5
                [587,591,557,563],   # 6
                [582,588,558,563],   # 7
                [577,582,561,567],   # 8
                [573,577,564,569],   # 9
                [569,574,568,573],   # 10
                [566,572,573,580],   # 11
                [564,568,580,586],   # 12
                [561,566,587,592],   # 13
                [561,566,593,599],   # 14
                [561,566,599,603],   # 15
                [562,568,606,609],   # 16, 17 and 18 omitted due to possible confusion with foreground star!
                [573,579,615,624],   # 19
                [578,582,620,625],   # 20
                [583,589,619,627],   # 21
                [592,596,618,624],   # 22
                [596,600,616,622],   # 23
                [601,607,606,613],   # 24
                [606,610,603,608],   # 25
                [607,612,598,603],   # 26
                [605,615,590,598],   # 27
                [607,613,585,590]])  # 28         

#########
# Ellipse parameters http://math.stackexchange.com/questions/426150/what-is-the-general-\
# equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate
X0 = 592.118767092
Y0 = 586.714283999
AA = 28
BB = 20
ALPHA = 0.15
subnx = boxes[0,3]-boxes[0,2]
subny = boxes[0,1]-boxes[0,0]
 
def gen_mask(AA, BB, X0, Y0, ALPHA):
    xx = np.arange(0,NX)
    yy = np.arange(0,NY)[:,None]
    xx = xx-X0
    yy = yy-Y0
    uu = xx*np.cos(ALPHA)-yy*np.sin(ALPHA)
    vv = xx*np.sin(ALPHA)+yy*np.cos(ALPHA)
    return ((uu/AA)**2 + (vv/BB)**2 > 0.5) & ((uu/AA)**2 + (vv/BB)**2 < 2.7) # True for points inside
mask = gen_mask(AA, BB, X0, Y0, ALPHA)
 
#########
# Fit ellipses for all observations
def help_ell(dummy, x0, y0, aa, bb, alpha, mag, sig, tilt, phi):
#     aa, bb, y0, x0, alpha, mag, sig = AA, BB, Y0, X0, ALPHA, 1., 1.
    xx = np.arange(0,subnx)
    yy = np.arange(0,subny)[:,None]
    xx = xx-x0
    yy = yy-y0
    uu = xx*np.cos(alpha)-yy*np.sin(alpha)
    vv = xx*np.sin(alpha)+yy*np.cos(alpha)
    rr = (uu/aa)**2 + (vv/bb)**2
#     print x0, y0, aa, bb, alpha, mag, sig, tilt, phi
    global last_ell
    last_ell = abs(mag)*np.exp(-np.abs(np.sqrt(rr)-1)**2/abs(sig))*(1+tilt*np.cos(np.arctan2(vv,uu)+phi))
    return np.reshape(last_ell, subnx*subny)
 
for img_path in files:
    img = fits.open(img_path)[0].data
    img = np.where(mask, img, 0.)[boxes[0,0]:boxes[0,1],boxes[0,2]:boxes[0,3]]
    volume = np.sum(np.abs(img))
 
# Fit
    guess = np.array([X0-boxes[0,2], Y0-boxes[0,0], AA, BB, ALPHA, 1., 0.05, 0.25, 0.])
    pars, covar = curve_fit(help_ell, np.arange(0, subnx*subny), np.reshape(img, subnx*subny), guess)
    alphs.append(-pars[4])
    pars[0] = pars[0] + boxes[0,2]
    pars[1] = pars[1] + boxes[0,0]
    res = np.sum(np.abs(img-last_ell)) / volume
    print("{0:17.13f} {1:17.13f} {2:17.13f}".format(pars[0], pars[1], res))
    coords.append([pars[0], pars[1]])
 
# Plots
    fig = plt.figure()
    plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(img)/100, vmax=np.amax(img)))
    fig.savefig('coord_fits/hst/' + img_path[-28:-5] + '_ell_img.pdf',bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)
    
    fig = plt.figure()
    plt.imshow(img-last_ell, cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
    fig.savefig('coord_fits/hst/' + img_path[-28:-5] + '_ell_res.pdf',bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)
    
#########
# Gaussian parameters
X0 = 592.118767092
Y0 = 586.714283999
AA = 1. # Amplitude
SIG = 10.
subnx = boxes[1,3]-boxes[1,2]
subny = boxes[1,1]-boxes[1,0]
 
def gen_mask(X0, Y0, SIG):
    xx = np.arange(0,NX)
    yy = np.arange(0,NY)[:,None]
    xx = xx-X0
    yy = yy-Y0
    return np.sqrt(xx**2+yy**2) < 1.85*SIG
 
#########
# Fit Gaussians for all observations
def help_gauss(dummy, x0, y0, aa, sig):
    xx = np.arange(0,subnx)
    yy = np.arange(0,subny)[:,None]
    xx = xx-x0
    yy = yy-y0
    rr = xx**2+yy**2
    global last_gau
    last_gau = aa*np.exp(-rr/sig**2)
    return np.reshape(last_gau, subnx*subny)
 
for img_path in files:
# Star 1
    img = fits.open(img_path)[0].data
    mask = gen_mask(Xstar1, Ystar1, SIG)
    img = np.where(mask, img, 0.)[boxes[2,0]:boxes[2,1],boxes[2,2]:boxes[2,3]]
    volume = np.sum(np.abs(img))
 
    guess = np.array([Xstar1-boxes[2,2], Ystar1-boxes[2,0], AA, SIG])
    pars, covar = curve_fit(help_gauss, np.arange(0, subnx*subny), np.reshape(img, subnx*subny), guess)
    pars[0] = pars[0] + boxes[2,2]
    pars[1] = pars[1] + boxes[2,0]
    if pars[0] < 496:
        print "WARNING, FIDDLED WITH STAR 1"
        pars[0] = 498.00040341
        pars[1] = 288.16928931
        
    res = np.sum(np.abs(img-last_gau)) / volume
    print("{0:17.13f} {1:17.13f} {2:17.13f} Star 1".format(pars[0], pars[1], res))
    temp = pars[0:2]
 
# Star 2
    img = fits.open(img_path)[0].data
    mask = gen_mask(Xstar2, Ystar2, SIG)
    img = np.where(mask, img, 0.)[boxes[3,0]:boxes[3,1],boxes[3,2]:boxes[3,3]]
    volume = np.sum(np.abs(img))
    
    guess = np.array([Xstar2-boxes[3,2], Ystar2-boxes[3,0], AA, SIG])
    pars, covar = curve_fit(help_gauss, np.arange(0, subnx*subny), np.reshape(img, subnx*subny), guess)
    pars[0] = pars[0] + boxes[3,2]
    pars[1] = pars[1] + boxes[3,0]
    res = np.sum(np.abs(img-last_gau)) / volume
    print("{0:17.13f} {1:17.13f} {2:17.13f} Star 2".format(pars[0], pars[1], res))
    stars = [temp[0], temp[1], pars[0], pars[1]]
    save_stars.append(stars)
 
# Gauss
    img = fits.open(img_path)[0].data
    mask = gen_mask(X0, Y0, SIG)
    img = np.where(mask, img, 0.)[boxes[1,0]:boxes[1,1],boxes[1,2]:boxes[1,3]]
    volume = np.sum(np.abs(img))
 
    guess = np.array([X0-boxes[1,2], Y0-boxes[1,0], AA, SIG])
    pars, covar = curve_fit(help_gauss, np.arange(0, subnx*subny), np.reshape(img, subnx*subny), guess)
    pars[0] = pars[0] + boxes[1,2]
    pars[1] = pars[1] + boxes[1,0]
    res = np.sum(np.abs(img-last_gau)) / volume
    print("{0:17.13f} {1:17.13f} {2:17.13f} Gauss".format(pars[0], pars[1], res))
    coords.append([pars[0], pars[1]])
 
# Plots
    fig = plt.figure()
    plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
    fig.savefig('coord_fits/hst/' + img_path[-28:-5] + '_gau_img.pdf',bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)
 
    fig = plt.figure()
    plt.imshow(img-last_gau, cmap='afmhot', origin='lower', interpolation=INTERPOLATION)
    fig.savefig('coord_fits/hst/' + img_path[-28:-5] + '_gau_res.pdf',bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)
 
# Additional points
    def coord2pix(coord1, coord2):
        ra1 = np.deg2rad(GAIAX1)
        ra2 = np.deg2rad(GAIAX2)
        dec1= np.deg2rad(90+GAIAY1)
        dec2= np.deg2rad(90+GAIAY2)
        xx = (coord1-GAIAX2)/(GAIAX1-GAIAX2)*(stars[0]-stars[2])+stars[2]
        yy = (coord2-GAIAY2)/(GAIAY1-GAIAY2)*(stars[1]-stars[3])+stars[3]
        coords.append([xx,yy])
        
    coord2pix(83.866533, 69.269747)
    coord2pix(83.866333, 69.269750)
    coord2pix(83.866642, 69.269744)
        
#     delr = np.rad2deg(np.arccos(np.cos(dec1)*np.cos(dec2)+np.sin(dec1)*np.sin(dec2)*np.cos(ra1-ra2)))
#     delXb= stars[2]-stars[0]
#     delYb= stars[3]-stars[1]
#     delrb= np.sqrt(delXb**2+delYb**2)
#     print("Check: {0:17} {1:17} {2:17}".format(delr*3600, delrb, delr/delrb*3600))
 
#########
# Sparse ellipse
X0 = 592.118767092
Y0 = 586.714283999
AA = 28
BB = 20
ALPHA = -0.15
SLICES = 32
ring = np.zeros((SLICES,2))
subnx = boxes[0,3]-boxes[0,2]
subny = boxes[0,1]-boxes[0,0]
 
def gen_mask(AA, BB, x0, y0, ALPHA, pie):
    xx = np.arange(0,subnx)
    yy = np.arange(0,subny)[:,None]
    xx = xx-x0
    yy = yy-y0
    uu = xx*np.cos(-ALPHA)-yy*np.sin(-ALPHA)
    vv = xx*np.sin(-ALPHA)+yy*np.cos(-ALPHA)
    theta = np.arctan2(yy, xx)
    return ((uu/AA)**2 + (vv/BB)**2 > 0.5) & ((uu/AA)**2 + (vv/BB)**2 < 2.7) & (theta > pie[0]) & (theta < pie[1])
 
#########
# Fit sparse ellipses for all observations
def help_ell(pts, x0, y0, aa, bb, alpha):
#    aa, bb, y0, x0, alpha = AA, BB, Y0, X0, ALPHA
    pts = pts.reshape((SLICES,2))
    inX = pts[:,0]-x0
    inY = pts[:,1]-y0
    inTheta = np.arctan2(inY, inX)
    
    gridTheta = np.arange(-np.pi, np.pi, 0.01)
    xx = aa*np.cos(gridTheta)
    yy = bb*np.sin(gridTheta)
    uu = xx*np.cos(alpha)-yy*np.sin(alpha)
    vv = xx*np.sin(alpha)+yy*np.cos(alpha)
    rotTheta = np.arctan2(vv, uu)
    uu = uu+x0
    vv = vv+y0
    uu = griddata(rotTheta, uu, inTheta)
    vv = griddata(rotTheta, vv, inTheta)
#    print x0, y0, aa, bb, alpha
    global last_sel
    last_sel = (uu, vv)
    return np.hstack((uu[:,np.newaxis], vv[:,np.newaxis])).reshape(2*SLICES)
 
for img_path in files:
    bup_img = fits.open(img_path)[0].data[boxes[0,0]:boxes[0,1],boxes[0,2]:boxes[0,3]]
    for pslice in range(0, SLICES):
        img = bup_img
        mask = gen_mask(AA, BB, X0-boxes[0,2], Y0-boxes[0,0], ALPHA, [pslice*2*np.pi/SLICES-np.pi,(pslice+1)*2*np.pi/SLICES-np.pi])
        img = np.where(mask, img, 0.)
#        plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(img)/100, vmax=np.amax(img)))
#        plt.show()
        ytemp, xtemp = np.unravel_index(img.argmax(), img.shape)
        ring[pslice,:] = np.array([xtemp+boxes[0,2], ytemp+boxes[0,0]])
 
    guess = np.array([X0, Y0, AA, BB, ALPHA])
    pars, covar = curve_fit(help_ell, ring.reshape(2*SLICES), ring.reshape(2*SLICES), guess)
    alphs.append(pars[4]) ####
    print("{0:17.13f} {1:17.13f}".format(pars[0], pars[1]))
    coords.append([pars[0], pars[1]])
 
    fig = plt.figure()
    plt.plot(ring[:,0]-boxes[0,2], ring[:,1]-boxes[0,0])
    plt.plot(last_sel[0]-boxes[0,2], last_sel[1]-boxes[0,0], 'r')
    img = fits.open(img_path)[0].data[boxes[0,0]:boxes[0,1],boxes[0,2]:boxes[0,3]]
    plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(img)/100, vmax=np.amax(img)))
    fig.savefig('coord_fits/hst/' + img_path[-28:-5] + '_sel_res.pdf',bbox_inches='tight', pad_inches=0.01)
#    plt.show()
    plt.close(fig)

#########
# Hotspots
AA = 2. # Amplitude
SIG = 5.
X0 = 592.118767092
Y0 = 586.714283999

def help_gauss(dummy, x0, y0, aa, sig):
    xx = np.arange(0,subnx)
    yy = np.arange(0,subny)[:,None]
    xx = xx-x0
    yy = yy-y0
    rr = xx**2+yy**2
    global last_gau
    last_gau = aa*np.exp(-rr/sig**2)
    return np.reshape(last_gau, subnx*subny)

#########
# Fit the boxes
fig = plt.figure()
ax = fig.gca()
for box in range(4,boxes.shape[0]):
    subnx = boxes[box,3]-boxes[box,2]
    subny = boxes[box,1]-boxes[box,0]
    boxX0 = (boxes[box,2]+boxes[box,3])/2.
    boxY0 = (boxes[box,0]+boxes[box,1])/2.
    
    img = fits.open(RINGEST)[0].data
    img = img[boxes[box,0]:boxes[box,1],boxes[box,2]:boxes[box,3]]
    volume = np.sum(np.abs(img))

    guess = np.array([boxX0-boxes[box,2], boxY0-boxes[box,0], AA, SIG])
    pars, covar = curve_fit(help_gauss, np.arange(0, subnx*subny), np.reshape(img, subnx*subny), guess)
    pars[0] = pars[0] + boxes[box,2]
    pars[1] = pars[1] + boxes[box,0]
    prod_pars.append(pars[0:2])
    res = np.sum(np.abs(img-last_gau)) / volume
    print("{0:17.13f} {1:17.13f} {2:17.13f} Box".format(pars[0], pars[1], res))
    hotspot_theta[box-4] = np.arctan2(pars[1]-Y0, pars[0]-X0)

# Plots
    ax.add_patch(patches.Rectangle(
        (boxes[box,2]-boxes[0,2], boxes[box,0]-boxes[0,0]),
        boxes[box,3]-boxes[box,2],
        boxes[box,1]-boxes[box,0],
        fill=False))      # remove background
        
    plt.plot(boxes[box,2]-boxes[0,2],boxes[box,0]-boxes[0,0],'b^')
    plt.plot(boxes[box,3]-boxes[0,2],boxes[box,1]-boxes[0,0],'rv')
    plt.plot(pars[0]-boxes[0,2], pars[1]-boxes[0,0], 'go')

img = fits.open(RINGEST)[0].data[boxes[0,0]:boxes[0,1],boxes[0,2]:boxes[0,3]]
plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(img)/100, vmax=np.amax(img)))
fig.savefig('coord_fits/hst/' + RINGEST[-28:-5] + '_fransson15b.pdf',bbox_inches='tight', pad_inches=0.01)
#plt.show()
plt.close(fig)

#########
# Fit radial, 1D Gaussians for all observations
AA = 28
BB = 20
ALPHA = -0.15
subnx = boxes[0,3]-boxes[0,2]
subny = boxes[0,1]-boxes[0,0]
ring = np.zeros((Nspots,2))

def help_1DGauss(pts, rshift, aa, sig):
#    print rshift
    global last_1Dgau
    last_1Dgau = aa*np.exp(-(pts-rshift)**2/sig**2) #-1 # -1 because heuristic rescaling of +0.1 before log to shift reference
    return last_1Dgau

def help_ell(pts, x0, y0, aa, bb, alpha):
#    aa, bb, y0, x0, alpha = AA, BB, Y0, X0, ALPHA
    pts = pts.reshape((Nspots,2))
    inX = pts[:,0]-x0
    inY = pts[:,1]-y0
    inTheta = np.arctan2(inY, inX)
    
    gridTheta = np.arange(-np.pi, np.pi, 0.01)
    xx = aa*np.cos(gridTheta)
    yy = bb*np.sin(gridTheta)
    uu = xx*np.cos(alpha)-yy*np.sin(alpha)
    vv = xx*np.sin(alpha)+yy*np.cos(alpha)
    rotTheta = np.arctan2(vv, uu)
    uu = uu+x0
    vv = vv+y0
    uu = griddata(rotTheta, uu, inTheta)
    vv = griddata(rotTheta, vv, inTheta)
#     print x0, y0, aa, bb, alpha
    global last_spt
    last_spt = (uu, vv)
    return np.hstack((uu[:,np.newaxis], vv[:,np.newaxis])).reshape(2*Nspots)
    
for img_path in files:
# Gauss
    img = fits.open(img_path)[0].data
    img = img[boxes[0,0]:boxes[0,1],boxes[0,2]:boxes[0,3]]

    img_pts = np.hstack((np.tile(np.arange(0,subnx),subny)[:,np.newaxis],np.arange(0,subny).repeat(subnx)[:,np.newaxis]))

    for spt in range(0,len(hotspot_theta)):
        pts = np.arange(15,40,0.01)
        lineX = X0 + pts*np.cos(hotspot_theta[spt]) - boxes[0,2]
        lineY = Y0 + pts*np.sin(hotspot_theta[spt]) - boxes[0,0]
        line = griddata(img_pts, img.reshape(subnx*subny), np.hstack((lineX[:,np.newaxis],lineY[:,np.newaxis])), method='linear', fill_value=-1)

        guess = np.array([25, AA, SIG])
        param_bounds=([15,-np.inf,-np.inf],[40,np.inf,np.inf])
#        plt.plot(line)
#        plt.show()
        try:
            pars, covar = curve_fit(help_1DGauss, pts, line, guess, bounds=param_bounds)
        except:
            print "WARNING 1D GAUSS FIT FAILED"
            pars = guess

        xtemp = X0 + pars[0]*np.cos(hotspot_theta[spt])# - boxes[0,2]
        ytemp = Y0 + pars[0]*np.sin(hotspot_theta[spt])# - boxes[0,0]
        ring[spt,:] = np.array([xtemp, ytemp])

    guess = np.array([X0, Y0, AA, BB, ALPHA])
    try:
        pars, covar = curve_fit(help_ell, ring.reshape(2*Nspots), ring.reshape(2*Nspots), guess)
        alphs.append(pars[4]) ####
        if img_path == RINGEST:
            prod_spt = pars
        if img_path == CHANDRAEST:
            chandra_spt = pars
        print img_path, pars[2]*25, pars[3]*25
    except:
        print "WARNING SPARSE ELLIPSE TO HOTSPOTS FAILED"
        pars = guess

    print("{0:17.13f} {1:17.13f}".format(pars[0], pars[1]))
    coords.append([pars[0], pars[1]])

    fig = plt.figure()
    plt.plot(ring[:,0]-boxes[0,2], ring[:,1]-boxes[0,0])
    plt.plot(last_spt[0]-boxes[0,2], last_spt[1]-boxes[0,0], 'r')
    img = fits.open(img_path)[0].data[boxes[0,0]:boxes[0,1],boxes[0,2]:boxes[0,3]]
    plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(img)/100, vmax=np.amax(img)))
    fig.savefig('coord_fits/hst/' + img_path[-28:-5] + '_spt_res.pdf',bbox_inches='tight', pad_inches=0.01)
#    plt.show()
    plt.close(fig)
    
#########
# Plot and wrap-up
img = fits.open(RINGEST)[0].data[boxes[0,0]:boxes[0,1],boxes[0,2]:boxes[0,3]]
coords = np.array(coords)
stars = np.array(save_stars)

fig = plt.figure()
col2 = ['.b','ob','^b','>b','vb','<b','.r','.g','.y','.m','or','og','oy','om','^r','^g','^y','^m',\
           '>r','>g','>y','>m','vr','vg','vy','vm','<r','<g','<y','<m',\
           '.c','oc','^c','>c','vc','<c','.k','ok','^k','>k','vk','<k']
col = ['ob','or','og','oy','om','oc','ow']

# Select what to keep
selected = np.zeros((33,6))
keep = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,27,28,29,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56]
counter = 0
for i in range(0,coords.shape[0]):
    if i >= NN and  i < 5*NN: # Set plot marker and figure out obs ID
        tempCol = col[np.mod(i-NN,4)+1]
        fileID = (i-NN)/4
    else:
        tempCol = col[i/NN]
        fileID = np.mod(i,NN)

    # Further selection
    if any(substring in files[fileID] for substring in ["1994","1995","1996","1997"]):
        tempCol = tempCol.replace("o","s")
    elif any(substring in files[fileID] for substring in ["1998","1999","2000","2001"]):
        tempCol = tempCol.replace("o","^")

    if i >= NN and  i < 5*NN:
        inverting_condition = 1 # :)
    else:
        if fileID in keep:
            print(counter, fileID) # Save the stuff to keep
            selected[np.mod(counter,33),2*(counter/33)] = coords[i,0]
            selected[np.mod(counter,33),2*(counter/33)+1] = coords[i,1]
            counter += 1
        
    plt.plot(coords[i,0]-boxes[0,2], coords[i,1]-boxes[0,0],tempCol)

print("Elliptic annulus x: mu = {0}, std = {1}, sem = {2}".format(np.mean(selected[:,0]),np.std(selected[:,0]),25*sem(selected[:,0])))
print("Elliptic annulus y: mu = {0}, std = {1}, sem = {2}".format(np.mean(selected[:,1]),np.std(selected[:,1]),25*sem(selected[:,1])))
print(  "Sparse ellipse x: mu = {0}, std = {1}, sem = {2}".format(np.mean(selected[:,2]),np.std(selected[:,2]),25*sem(selected[:,2])))
print(  "Sparse ellipse y: mu = {0}, std = {1}, sem = {2}".format(np.mean(selected[:,3]),np.std(selected[:,3]),25*sem(selected[:,3])))
print(        "Hotspots x: mu = {0}, std = {1}, sem = {2}".format(np.mean(selected[:,4]),np.std(selected[:,4]),25*sem(selected[:,4])))
print(        "Hotspots y: mu = {0}, std = {1}, sem = {2}".format(np.mean(selected[:,5]),np.std(selected[:,5]),25*sem(selected[:,5])))

pdb.set_trace()

# Convert to world coordinates
s0 = (np.sum(stars[:15,0])+np.sum(stars[27:30,0])+np.sum(stars[42:,0]))/33.
s1 = (np.sum(stars[:15,1])+np.sum(stars[27:30,1])+np.sum(stars[42:,1]))/33.
s2 = (np.sum(stars[:15,2])+np.sum(stars[27:30,2])+np.sum(stars[42:,2]))/33.
s3 = (np.sum(stars[:15,3])+np.sum(stars[27:30,3])+np.sum(stars[42:,3]))/33.

ra = (np.mean(selected[:,4])-s2)*(GAIAX1-GAIAX2)/(s0-s2)+GAIAX2
dec = (np.mean(selected[:,5])-s3)*(-GAIAY1+GAIAY2)/(s1-s3)-GAIAY2
spt_center = SkyCoord(ra, dec, frame='icrs', unit='deg')
print("Hotspots: {0} {1} {2}".format(spt_center.to_string('hmsdms'),sem(selected[:,4])*0.025, sem(selected[:,5])*0.025))

ra = (np.mean(selected[:,0])-s2)*(GAIAX1-GAIAX2)/(s0-s2)+GAIAX2
dec = (np.mean(selected[:,1])-s3)*(-GAIAY1+GAIAY2)/(s1-s3)-GAIAY2
ell_center = SkyCoord(ra, dec, frame='icrs', unit='deg')
print("Elliptical annulus: {0} {1} {2}".format(ell_center.to_string('hmsdms'),sem(selected[:,0])*0.025, sem(selected[:,1])*0.025))

#########
# Plots
# coords
plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(img)/100, vmax=np.amax(img)))
fig.savefig('coord_fits/hst/' + RINGEST[-28:-5] + '_coords.pdf',bbox_inches='tight', pad_inches=0.01)
plt.show()
plt.close(fig)

## all coords
#for i in range(0,coords.shape[0]):
#    if i >= NN and  i < 5*NN: # Set plot marker and figure out obs ID
#        tempCol = col[np.mod(i-NN,4)+1]
#        fileID = (i-NN)/4
#    else:
#        tempCol = col[i/NN]
#        fileID = np.mod(i,NN)
#
#    # Further selection
#    if any(substring in files[fileID] for substring in ["1994","1995","1996","1997"]):
#        tempCol = tempCol.replace("o","s")
#    elif any(substring in files[fileID] for substring in ["1998","1999","2000","2001"]):
#        tempCol = tempCol.replace("o","^")
#        
#    plt.plot(coords[i,0]-boxes[0,2], coords[i,1]-boxes[0,0],tempCol)
#
#plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(img)/100, vmax=np.amax(img)))
#fig.savefig('coord_fits/hst/' + RINGEST[-28:-5] + '_allcoo.pdf',bbox_inches='tight', pad_inches=0.01)
#plt.show()
#plt.close(fig)

# zoom2
fig = plt.figure()
ax = fig.gca()
for i in range(0,coords.shape[0]):
    if (np.mod(i,NN) in keep) and (i/NN == 0 or i/NN == 6):
        tempCol = 'ob' if i/NN == 0 else 'sc'
        fileID = np.mod(i,NN)
        if any(substring in files[fileID] for substring in ["1994","1995","1996","1997","1998","1999","2000","2001"]):
            if '_r_' in files[fileID]:
                plt.plot(25*(coords[i,0]-np.mean(selected[:,4])), 25*(coords[i,1]-np.mean(selected[:,5])),tempCol.replace('ob','pr').replace('sc','*m'), markersize=10)
        else:
            tempCol = tempCol.replace('ob','dr').replace('sc','vm') if '_r_' in files[fileID] else tempCol
            plt.plot(25*(coords[i,0]-np.mean(selected[:,4])), 25*(coords[i,1]-np.mean(selected[:,5])),tempCol, markersize=10)

ellipse = patches.Ellipse(xy=(25*(np.mean(selected[:,0])-np.mean(selected[:,4])), 25*(np.mean(selected[:,1])-np.mean(selected[:,5]))), width=2*25*sem(selected[:,0]), height=2*25*sem(selected[:,1]), edgecolor='g', fc='None', lw=6)
ax.add_artist(ellipse)
ellipse.set_zorder(5)
ellipse = patches.Ellipse(xy=(0, 0), width=2*25*sem(selected[:,4]), height=2*25*sem(selected[:,5]), edgecolor='y', fc='None', lw=6)
ax.add_artist(ellipse)
ellipse.set_zorder(5)

ellipse = patches.Ellipse(xy=(25*(np.mean(selected[:,0])-np.mean(selected[:,4])), 25*(np.mean(selected[:,1])-np.mean(selected[:,5]))), width=2*25*np.std(selected[:,0]), height=2*25*np.std(selected[:,1]), edgecolor='g', fc='None', lw=3)
ax.add_artist(ellipse)
ellipse.set_zorder(5)
ellipse = patches.Ellipse(xy=(0, 0), width=2*25*np.std(selected[:,4]), height=2*25*np.std(selected[:,5]), edgecolor='y', fc='None', lw=3)
ax.add_artist(ellipse)
ellipse.set_zorder(5)

fig.savefig('coord_fits/hst/' + RINGEST[-28:-5] + '_zoom2.pdf',bbox_inches='tight', pad_inches=0.01)
plt.show()
plt.close(fig)

# ell_sel
fig = plt.figure()
plt.plot(last_sel[0]-boxes[0,2], last_sel[1]-boxes[0,0], 'g')
plt.imshow(last_ell, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(img)/100, vmax=np.amax(img)))
fig.savefig('coord_fits/hst/' + RINGEST[-28:-5] + '_ell_sel.pdf',bbox_inches='tight', pad_inches=0.01)
#plt.show()
plt.close(fig)

#########
# Products
# gau_spt
from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams.update({'font.size': 22})

fig = plt.figure(figsize=(10, 7.5))
ax = fig.gca()
for box in range(4,boxes.shape[0]):
    ax.add_patch(patches.Rectangle(
        (boxes[box,2]-boxes[0,2]-0.5, boxes[box,0]-boxes[0,0]-0.5),
        boxes[box,3]-boxes[box,2],
        boxes[box,1]-boxes[box,0],
        fill=False, zorder=2))      # remove background
        
#    plt.plot(boxes[box,2]-boxes[0,2],boxes[box,0]-boxes[0,0],'b^')
#    plt.plot(boxes[box,3]-boxes[0,2],boxes[box,1]-boxes[0,0],'rv')
    plt.plot(prod_pars[box-4][0]-boxes[0,2]-0.5, prod_pars[box-4][1]-boxes[0,0]-0.5, 'o', color='#2ca02c', ms=8, mew=0, zorder=4)

# Image with scale
img = fits.open(RINGEST)[0].data[boxes[0,0]:boxes[0,1],boxes[0,2]:boxes[0,3]]
plt.imshow(img, cmap='afmhot', origin='lower', interpolation=INTERPOLATION, norm=LogNorm(vmin=np.amax(img)/300, vmax=np.amax(img)), zorder=0)
plt.plot([78, 82], [58, 58], 'w', lw=2)
plt.annotate("$0.''1$", xy=(77.5, 59.75), color="w")

def pix2icrf(xpix, ypix): # to ICRF
    ra = (xpix-s2)*(GAIAX1-GAIAX2)/(s0-s2)+GAIAX2
    dec = (ypix-s3)*(-GAIAY1+GAIAY2)/(s1-s3)-GAIAY2
    spt_center = SkyCoord(ra, dec, frame='icrs', unit='deg')
    return spt_center.to_string('hmsdms').encode('ascii','ignore')

def coord2pix(coord1, coord2):
    xx = (coord1-GAIAX2)/(GAIAX1-GAIAX2)*(s0-s2)+s2
    yy = (coord2-GAIAY2)/(GAIAY1-GAIAY2)*(s1-s3)+s3
    return xx, yy

# Stis slits
sts_dec = -69-16/60.-18.36/3600 # Declination at 1 (reference for the observation)
sts_off = np.array([-2.4435, -2.4625, -2.4815, -2.5005, -2.5195]) # in seconds
sts_off = np.array([-2.4435, -2.4625, -2.4815, -2.5005, -2.5195, -2.5385]) # in seconds
sts_off = sts_off*np.cos(np.deg2rad(sts_dec))*15 # to arcsec
sts_off = sts_off+0.05
sts_ra  = 5/24.*360+15*35/60.+15*30.517/3600.
for sts_ii in sts_off:
#    sts_xx, sts_yy = coord2pix(sts_ra+sts_ii/3600/np.cos(np.deg2rad(sts_dec)), sts_dec)
    sts_xx = 66.93-sts_ii/0.025
    plt.axvline(x=sts_xx-boxes[0,2]-0.5, linewidth=2, color = '#7f7f7f', zorder=1)

sts_box_xx = 66.93-sts_off/0.025
sts_pix = 0.05071/0.025
plt.plot(sts_box_xx[1:4]-boxes[0,2]-0.5, np.mean(selected[:,5])-2.5*sts_pix-3.04605225/25*np.ones(3)-boxes[0,0]-0.5, color = '#7f7f7f', zorder=1, lw=2)
plt.plot(sts_box_xx[1:4]-boxes[0,2]-0.5, np.mean(selected[:,5])+2.5*sts_pix-3.04605225/25*np.ones(3)-boxes[0,0]-0.5, color = '#7f7f7f', zorder=1, lw=2)

# xlabels
xcent = np.mean(selected[:,4])
xlabel = []
xlabel.append(pix2icrf(xcent-30,0).split()[0][6:-1]+"$^\mathrm{s}$")
xlabel.append(pix2icrf(xcent,   0).split()[0][6:-1]+"$^\mathrm{s}$")
xlabel.append(pix2icrf(xcent+30,0).split()[0][6:-1]+"$^\mathrm{s}$")
xcent -= boxes[0,2]
xcent -= 0.5
plt.xticks([xcent-30, xcent, xcent+30], xlabel, rotation=45)
ax.xaxis.get_major_formatter().set_offset_string("$5^\mathrm{h}\,35^\mathrm{m}$")

# ylabels
ycent = np.mean(selected[:,5])
ylabel = []
ylabel.append(pix2icrf(0,ycent-30).split()[1][7:-1]+"$''$")
ylabel.append(pix2icrf(0,ycent   ).split()[1][7:-1]+"$''$")
ylabel.append(pix2icrf(0,ycent+30).split()[1][7:-1]+"$''$")
ycent -= boxes[0,0]
ycent -= 0.5
plt.yticks([ycent-30, ycent, ycent+30], ylabel, rotation=45)
ax.yaxis.get_major_formatter().set_offset_string("$-69^\circ\,16'$")
ax.axis([0,boxes[0,3]-boxes[0,2]-1,0,boxes[0,1]-boxes[0,0]-1])

# Kick and ellipse
gridTheta = np.arange(-np.pi, np.pi+0.01, 0.001)
xx = prod_spt[2]*np.cos(gridTheta)
yy = prod_spt[3]*np.sin(gridTheta)
uu = xx*np.cos(prod_spt[4])-yy*np.sin(prod_spt[4])
vv = xx*np.sin(prod_spt[4])+yy*np.cos(prod_spt[4])
plt.plot(uu+prod_spt[0]-boxes[0,2]-0.5,vv+prod_spt[1]-boxes[0,0]-0.5,'#1f77b4', lw=2, zorder=3)
plt.plot(xcent+3.6*np.cos(gridTheta), ycent+3.6*np.sin(gridTheta), 'w', lw=2, zorder=4)

# Reynolds stuff
reynolds_optical_ra = 5/24.*360+15*35/60.+15*27.968/3600.
reynolds_radio_ra = 5/24.*360+15*35/60.+15*27.994/3600.
reynolds_optical_dec = 69+16/60.+11.09/3600.
reynolds_radio_dec = 69+16/60.+11.08/3600.

# Correct for proper motion, van_der_marel16
# 1991--2015, 24 years from reynolds95 to gaia16
pm_ra = 24*1.905/25
pm_dec = -24*0.275*2/25 # Last factor of two for rotation of LMC

reynolds_optical_xx, reynolds_optical_yy = coord2pix(reynolds_optical_ra, reynolds_optical_dec)
reynolds_radio_xx, reynolds_radio_yy = coord2pix(reynolds_radio_ra, reynolds_radio_dec)
reynolds_optical_xx -= pm_ra
reynolds_optical_yy -= pm_dec
reynolds_radio_xx -= pm_ra
reynolds_radio_yy -= pm_dec
plt.plot(reynolds_optical_xx-boxes[0,2]-0.5, reynolds_optical_yy-boxes[0,0]-0.5, 'd', color='#bcbd22', ms=9.5, zorder=5)
plt.plot(reynolds_radio_xx-boxes[0,2]-0.5, reynolds_radio_yy-boxes[0,0]-0.5, 's', color='#17becf', ms=7, zorder=6)

# Compass
ax.arrow(2, 2, 0, 8, head_width=1, head_length=1, fc='w', ec='w', linewidth=2)
ax.arrow(2, 2, 8, 0, head_width=1, head_length=1, fc='w', ec='w', linewidth=2)
ax.add_patch(patches.Rectangle((xcent-5/25.,ycent-6/25.),14/25., 12/25.,facecolor='w', linewidth=0, zorder=3))
ax.annotate("W", xy=(12,1), color="w")
ax.annotate("N", xy=(1,12), color="w")

# save
fig.savefig('art/figs/prod_gau_spt.pdf',bbox_inches='tight', pad_inches=0.1)
plt.show()
plt.close(fig)
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

#########
# center
fig = plt.figure(figsize=(10, 7.5))
ax = fig.gca()
for i in range(0,coords.shape[0]):
    if (np.mod(i,NN) in keep) and i/NN == 6:
        tempCol = 'sb'
        fileID = np.mod(i,NN)
        if any(substring in files[fileID] for substring in ["1994","1995","1996","1997","1998","1999","2000","2001"]):
            if '_r_' in files[fileID]:
                plt.plot(25*(coords[i,0]-np.mean(selected[:,4])), 25*(coords[i,1]-np.mean(selected[:,5])),tempCol.replace('sc','*r'), markersize=10)
        else:
            tempCol = tempCol.replace('sb','or') if '_r_' in files[fileID] else tempCol
            plt.plot(25*(coords[i,0]-np.mean(selected[:,4])), 25*(coords[i,1]-np.mean(selected[:,5])),tempCol, ms=10, mew=0)


# sanduleak
sanduX, sanduY = coord2pix(83.86653333333334, 69.26974722222222)
#plt.plot(25*(sanduX-np.mean(selected[:,4])), 25*(sanduY-np.mean(selected[:,5])),"*y", markersize=10)

#ellipse = patches.Ellipse(xy=(25*(np.mean(selected[:,0])-np.mean(selected[:,4])), 25*(np.mean(selected[:,1])-np.mean(selected[:,5]))), width=2*25*np.std(selected[:,0]), height=2*25*np.std(selected[:,1]), edgecolor='#2ca02c', fc='None', lw=3, linestyle='dashed')
#ax.add_artist(ellipse)
#ellipse.set_zorder(5)
#ellipse = patches.Ellipse(xy=(0, 0), width=2*25*np.std(selected[:,4]), height=2*25*np.std(selected[:,5]), edgecolor='k', fc='None', lw=3)
plt.plot(25*(np.mean(selected[:,0])-np.mean(selected[:,4])), 25*(np.mean(selected[:,1])-np.mean(selected[:,5])), c='#2ca02c', marker='p', ms=20, zorder=5)
plt.plot(0, 0, 'xk', ms=20, zorder=5, mew=3)
ellipse = patches.Ellipse(xy=(0, 0), width=2*6, height=2*4, edgecolor='k', fc='None', lw=3)
ax.add_artist(ellipse)
ellipse.set_zorder(5)

# xlabels
shift = 4/25.
xcent = np.mean(selected[:,4])
xlabel = []
ax.set_xlim([-5,9])
xlabel.append(pix2icrf(xcent-shift,0).split()[0][6:-1]+"$^\mathrm{s}$")
xlabel.append(pix2icrf(xcent,      0).split()[0][6:-1]+"$^\mathrm{s}$")
xlabel.append(pix2icrf(xcent+shift,0).split()[0][6:-1]+"$^\mathrm{s}$")
xlabel.append(pix2icrf(xcent+2*shift,0).split()[0][6:-1]+"$^\mathrm{s}$")
#xlabel.append(pix2icrf(xcent+3*shift,0).split()[0][6:-1]+"$^\mathrm{s}$")
xcent -= boxes[0,2]
plt.xticks([-shift*25, 0, shift*25, shift*50], xlabel, rotation=45)
ax.xaxis.get_major_formatter().set_offset_string("$5^\mathrm{h}\,35^\mathrm{m}$")

# ylabels
shift = 5/25.
ycent = np.mean(selected[:,5])
ax.set_ylim([-6,6])
ylabel = []
ylabel.append(pix2icrf(0,ycent-shift).split()[1][7:-1]+"$''$")
ylabel.append(pix2icrf(0,ycent   ).split()[1][7:-1]+"$''$")
ylabel.append(pix2icrf(0,ycent+shift).split()[1][7:-1]+"$''$")
ycent -= boxes[0,0]
plt.yticks([-shift*25, 0, shift*25], ylabel, rotation=45)
ax.yaxis.get_major_formatter().set_offset_string("$-69^\circ\,16'$")
#ax.axis([0,boxes[0,3]-boxes[0,2]-1,0,boxes[0,1]-boxes[0,0]-1])

# Compass
prevx = boxes[0,3]-boxes[0,2]-1
currx = 14
prevy = boxes[0,1]-boxes[0,0]-1
curry = 12
tranx = 1./prevx*currx
trany = 1./prevy*curry

# Scale and arrows
plt.plot([5, 7], [4, 4], 'k', lw=3)
plt.annotate("2~mas", xy=(5.3, 4.3), color='k')
ax.arrow(2*tranx-5, 2*trany-6, 0, 8*trany, head_width=1*tranx, head_length=1*trany, fc='k', ec='k', linewidth=2)
ax.arrow(2*tranx-5, 2*trany-6, 8*tranx, 0, head_width=1*trany, head_length=1*tranx, fc='k', ec='k', linewidth=2)
ax.annotate("W", xy=(12*tranx-5,1*trany-6), color="k")
ax.annotate("N", xy=(1*tranx-5,12*trany-6), color="k")

ax.set_aspect('equal')
fig.savefig('art/figs/prod_center.pdf',bbox_inches='tight', pad_inches=0.01)
plt.show()
plt.close(fig)

pdb.set_trace()
