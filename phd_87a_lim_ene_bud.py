# Dennis Alp 2017-11-21
# Find models that are consistent with oberved limits and energy budget
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/box/bin/phd_87a_uplim_ene_bud.py

import numpy as np
import os
import pdb
from glob import glob
from datetime import date

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from scipy.integrate import simps
from scipy.optimize import curve_fit
from scipy.ndimage.interpolation import shift
from astropy.io import fits

from sympy.stats import P, E, variance, Die, Normal
from sympy import Eq, simplify

from phd_87a_red import cor_red

################################################################
#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

################################################################
# Parameters
# Logistics
#WRK_DIR = "/Users/silver/box/phd/pro/87a/lim/glo"
#os.chdir(WRK_DIR) #Move to designated directory
NN = 8

########
# cgs constants
sb = 5.67051e-5
mp = 1.6726231e-24
distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e18 # cm
kpc = 1e3*pc # cm
DD = distance*kpc
Lsun = 3.826e33 # erg s-1
Msun = 1.989e33 # g
uu = 1.660539040e-24 # g
sT = 6.6524e-25 # cm2
mti = 43.9596901*uu # g
cc = 29979245800. # cm s-1
GG = 6.67259e-8 # cm3 g-1 s-2
hh = 6.6260755e-27  # erg s
kB = 1.380658e-16 # erg K-1
kev2erg = 1.60218e-9 # erg keV-1
kev2hz = 2.417989262503753e+17 # Hz keV-1

kth_blue = (25/255., 84/255., 166/255.)
kth_teal = (36/255.,160/255., 216/255.)

# Redenning
RV = 3.1
EBV = 0.19

################################################################
# Help functions
def tic_hlp(xx):
    return ["$10^{%i}$" % np.log10(tmp) for tmp in xx]

def lam2nu(lam):
    return cc/lam

def nu2lam(nu):
    return cc/nu
    
def Bnu(nu, TT):
    return 2*hh*nu**3/cc**2*1/(np.exp(hh*nu/(kB*TT))-1)

def get_lim(mod):
    grd_mod = griddata(mod_nu, mod, nu, method='linear', fill_value=0)
    return np.amin(nu*fn/grd_mod)*mod

def emod2nufnu(ee, emod):
    ee = ee*kev2hz
    nufnu = emod*kev2erg/kev2hz*ee
    return ee, nufnu

def from_xspec_mod(bins, flx):
    nu = bins[:-1]*kev2hz
    width = np.diff(bins)*kev2hz
    nufnu = flx*bins[:-1]*kev2erg*nu/width
    return nu, nufnu

def mod2nufnu(ee, mod):
    ee = ee*kev2hz
    nufnu = mod*kev2erg/kev2hz*ee**2
    return ee, nufnu

def get_crab_pulsar():
    dat = np.loadtxt('/Users/silver/box/sci/lib/b/buhler14_crab_pulsar.txt')
    nu = dat[:,0]
    nfn = dat[:,1]
    return nu, nfn

def get_crab_nebula():
    dat = np.loadtxt('/Users/silver/box/sci/lib/b/buhler14_crab_nebula.txt')
    nu = dat[:,0]
    nfn = dat[:,1]
    return nu, nfn

def get_cco_cas(): # Cas A, posselt13
    dat = np.loadtxt('/Users/silver/box/sci/lib/p/posselt13.txt')
    ee = dat[:,0]
    emod = dat[:,1]
    flx = np.loadtxt('/Users/silver/box/sci/lib/p/posselt13.flx')
    ene = np.loadtxt('/Users/silver/box/sci/lib/p/posselt13.ene')
    xx, yy = from_xspec_mod(ene, flx)
    nu, nufnu = emod2nufnu(ee, emod)
    return nu, nufnu, xx, yy

# Variance of products
def var_pro(EX, SX, EY, SY):
    return

# Days for scaling
SNDATE = date(1987, 2, 23)
def get_days(fname): # 1994-09-24_drz_al2.fits
    yr = int(fname[0:4])
    month = int(fname[5:7].lstrip("0"))
    day = int(fname[8:10].lstrip("0"))
    return (date(yr, month, day)-SNDATE).days*24*60*60

def get_ene_bud_lim():
    def get_ir_lum():
        # Date, microns, mJy, sd
        dat = np.array([
            [get_days('2011-08-15'), 350,   51,    10], #lakicevic12b, APEX
#            [get_days('2012-08-07'), 870,    5,     1], #indebetouw14, ALMA
#            [get_days('2012-09-20'), 440,   50,    10], #indebetouw14, ALMA
            [get_days('2010-07-01'), 100,  98.3,  8.5], #matsuura15, Herschel
            [get_days('2010-07-01'), 160, 169.6, 11.1], #matsuura15, Herschel
            [get_days('2010-07-01'), 250, 123.3, 13.4], #matsuura15, Herschel
            [get_days('2010-07-01'), 350,  53.8, 18.1], #matsuura15, Herschel
            [get_days('2012-01-13'), 100,  82.4,  4.5], #matsuura15, Herschel
            [get_days('2012-01-13'), 160, 153.0,  9.0], #matsuura15, Herschel
            [get_days('2012-01-13'), 250, 110.7, 25.2], #matsuura15, Herschel
            [get_days('2012-01-13'), 350,  69.3, 22.8]]) #matsuura15, Herschel

        scal_nu = cc/(dat[:,1]*1e-4)
        scal_ff = dat[:,2]*np.exp(-(dat[:,0]-tt)/tau)
        scal_sd = dat[:,3]*np.exp(-(dat[:,0]-tt)/tau)
        
        def planck(nu, CC, TT):
            return nu**3*CC/(np.exp(hh*nu/(kB*TT))-1)

        guess = [1e-33, 22]
        params, covar = curve_fit(planck, scal_nu, scal_ff, guess, sigma=scal_sd)
        nu = np.logspace(10, 13, 1000)
        plt.errorbar(scal_nu, scal_ff, scal_sd, fmt='o')
        plt.loglog(nu, planck(nu, params[0], params[1]))
        plt.show()
        ir_lum = simps(planck(nu, params[0], params[1]), nu)/1e26*4*np.pi*DD**2/Lsun
        
        return mu, sd

    def get_fd():
        tau = np.array([0, 0.5, 0.79, 1, 1.11, 1.5])
        mti = np.array([0.68, 1.07, 1.32, 1.55, 1.68, 2.12])
        coe = np.polyfit(mti, tau, 2)
#        plt.plot(mti, tau)
#        xx = np.linspace(-3, 4, 100)
#        plt.plot(xx, np.polyval(coe, xx))
#        plt.show()
        return coe
    
    # Get initial Ti input
    scale = (51.2/50)**2
    ti_half_life = 58.9 # yr
    Epositron = 596 #keV
    tau = ti_half_life/np.log(2)*365.25*24*3600
    tt = 9090.*24*3600
    ti_ini = scale*1.5e-4*Msun/mti
    ti_rate = ti_ini/tau*np.exp(-tt/tau)
    ti_lum = ti_rate*Epositron*kev2erg
    ti_std = 0.1822/1.5*ti_lum # norm.cdf(1.8, 1.5, 0.1822)-norm.cdf(1.2, 1.5, 0.1822)
    
    scale = np.exp(-575/(85.0*365.26))
    l1 = scale*310
    s1 = scale*27
    l2 = 270
    s2 = 19
    mu_dwek15 = Lsun*s1**2*s2**2/(s1**2+s2**2)*(l1/s1**2+l2/s2**2)
    std_dwek15 = Lsun*np.sqrt(s1**2*s2**2/(s1**2+s2**2))
    scale = (51.2/50)**2
    mu_dwek15 *= scale
    std_dwek15 *= scale
    get_fd()

# Numerical verification in case on does not trust SymPy.
    num_ti_in = np.random.normal(ti_lum, ti_std, 10000000)
    num_heat = np.random.uniform(0.55, 0.85, 10000000)
    num_dust = np.random.uniform(0.5, 0.7, 10000000)
    num_ir = np.random.normal(mu_dwek15, std_dwek15, 10000000)
    
    num_co = (num_ir - (num_heat + (1-num_heat)*num_dust)*num_ti_in)/num_dust
    num_co =  num_ir - (num_heat + (1-num_heat)*num_dust)*num_ti_in
    upp_lim_co = np.percentile(num_co, 99.7300203936740)
    med_lim_co = np.percentile(num_co, 50)
    low_lim_co = np.percentile(num_co, 15.865525393145702)
    hig_lim_co = np.percentile(num_co, 84.1344746068543)
    print('This is the 3-sigma upper limit on a CO:', upp_lim_co, upp_lim_co/Lsun)
    print('Median:', med_lim_co, '+', hig_lim_co-med_lim_co, '-', med_lim_co-low_lim_co)
    print('Median:', med_lim_co/Lsun, '+', (hig_lim_co-med_lim_co)/Lsun, '-', (med_lim_co-low_lim_co)/Lsun)
    print(np.mean(num_co)/Lsun, np.std(num_co)/Lsun, upp_lim_co/Lsun)

    return np.mean(num_co), np.std(num_co), upp_lim_co
    
def get_ene_bud_lim_df():
    uv_lim = 4e-17
    uv_frac = 0.3
    return uv_lim*2000*4*np.pi*DD**2/uv_frac

################################################################
# Data
def get_all_lim():
    nu = []
    fn = []
    co = []
    
    ########
    # callingham16, entire remnant
    tmp = 1e9*np.array([0.076, 0.084, 0.092, 0.099, 0.107, 0.115, 0.123, 0.130, 0.143, 0.150, 0.158, 0.166, 0.174, 0.181, 0.189, 0.197, 0.204, 0.212, 0.219, 0.227, 1.375, 1.375, 2.351, 2.351, 4.788, 4.788, 8.642, 8.642])
    nu.append(tmp)
    tmp = 1e-23*np.array([5.1, 4.9, 4.7, 4.6, 4.5, 4.2, 4.0, 3.9, 3.6, 3.4, 3.3, 3.1, 3.0, 2.9, 2.8, 2.7, 2.5, 2.5, 2.4, 2.3, 0.58, 0.58, 0.43, 0.42, 0.28, 0.30, 0.18, 0.17])
    fn.append(tmp)
    co.append(0*np.ones(len(tmp)))
    
    # zanardo14, estimate, not understandable, subtracted, modeled, something bogus, MEM image
    tmp = np.logspace(np.log10(102e9), np.log10(672e9), NN)
    nu.append(tmp)
    tmp = 3e-26*np.ones(NN)
    fn.append(tmp)
    co.append(1*np.ones(len(tmp)))
    
    # potter09, estimate 0.3 mJy, 3sigma 0.9 mJy, this subtracts a model torus and fits a central point source, unclear if MEM image
    tmp = [36.2e9]
    nu.append(tmp)
    tmp = [0.3e-26]
    fn.append(tmp)
    co.append(2*np.ones(len(tmp)))
    
    # ng11, 3sigma limit, VLBI image
    tmp = [1.7e9]
    nu.append(tmp)
    tmp = [0.3e-26]
    fn.append(tmp)
    co.append(19*np.ones(len(tmp)))
    
    # ng08, 3sigma limit, this subtracts a model torus and fits a central point source, MEM image
    tmp = [9e9]
    nu.append(tmp)
    tmp = [0.3e-26]
    fn.append(tmp)
    co.append(3*np.ones(len(tmp)))
    
    # lakicevic12, 2sigma limit, simply taking the flux in the central region of a MEM image
    tmp = [94e9]
    nu.append(tmp)
    tmp = [1e-26]
    fn.append(tmp)
    co.append(4*np.ones(len(tmp)))
    
    # zanardo13, simply taking the flux in the central region of a MEM image, 2.2 mJy if dust corrected
    tmp = [44e9]
    nu.append(tmp)
    tmp = [2.2e-26]
    fn.append(tmp)
    co.append(5*np.ones(len(tmp)))
    
    # matsuura15, entire remnant, dust except 500 microns, which is 3sigma limit
    tmp = lam2nu(1e-4*np.array([70, 100, 160, 250, 350, 500]))
    nu.append(tmp)
    tmp = 1e-26*np.array([45.4, 82.4, 153., 110.7, 69.3, 60.])
    fn.append(tmp)
    co.append(6*np.ones(len(tmp)))

    # bouchet04, gemini resolved, maybe dust in ejecta
    tmp = [lam2nu(10.335*1e-4)]
    nu.append(tmp)
    tmp = [1e-26*np.array(0.32)]
    fn.append(tmp)
    co.append(17*np.ones(len(tmp)))

    # bouchet04, gemini resolved, 3sigma upper limit, narrower filter than bouchet04
    tmp = [lam2nu(11.655*1e-4)]
    nu.append(tmp)
    tmp = [1e-26*np.array(0.34)]
    fn.append(tmp)
    co.append(18*np.ones(len(tmp)))
        
    # arendt16, entire remnant
    tmp = lam2nu(1e-4*np.array([3.6, 4.5, 5.8, 8., 24., 3.6, 4.5, 5.8, 8., 24.]))
    nu.append(tmp)
    tmp = 1e-26*np.array([1.52, 2.17, 4.08, 13.61, 75.7, 0.99, 1.2, 1.62, 4.69, 26.3])
    fn.append(tmp)
    co.append(7*np.ones(len(tmp)))
    
    # grebenev12, entire remnant
    tmp = np.logspace(np.log10(4.836e18), np.log10(14.51e18), NN)
    nu.append(tmp)
    tmp = 3e35/(4*np.pi*(51.2*kpc)**2)
    tmp = tmp/(-10*(1.451e19**-0.1-4.836e18**-0.1))
    tmp = tmp*np.array(nu[-1])**-1.1
    fn.append(tmp)
    co.append(9*np.ones(len(tmp)))
    #print np.trapz(fn[-1], nu[-1])*4*np.pi*(51.2*kpc)**2    
    
    # ackermann16, entire remnant, 95% limit
    tmp = np.logspace(np.log10(2.418e23), np.log10(2.418e24), NN)
    nu.append(tmp)
    tmp = 1.6021766e-6*7.8e-7
    tmp = tmp/np.log(2.418e24/2.418e23)
    tmp = tmp*np.array(nu[-1])**-1
    fn.append(tmp)
    co.append(15*np.ones(len(tmp)))
    #print np.trapz(fn[-1], nu[-1])/1.6021766e-6
    
    # hess15, entire remnant, 99%
    tmp = np.logspace(np.log10(2.418e26), np.log10(2.418e27), NN)
    nu.append(tmp)
    tmp = 2.2e34/(4*np.pi*(51.2*kpc)**2)
    tmp = tmp*-0.8/(2.418e27**-0.8-2.418e26**-0.8)
    tmp = tmp*np.array(nu[-1])**-1.8
    fn.append(tmp)
    co.append(16*np.ones(len(tmp)))
    #print np.trapz(fn[-1], nu[-1])*4*np.pi*(51.2*kpc)**2

    ########
    # New limits
    # ALMA
    tmp = [213e9, 233e9, 247e9]
    nu.append(tmp)
    tmp = [3.28467956e-26, 4.78012075e-26, 2.29446108e-26] # These are from phd_87a_uplim_plt_lim.py
    tmp = [4.49094e-26, 3.16767e-26, 2.46548e-26] # These are from phd_87a_uplim_plt_lim.py
    tmp = [0.10878717422692565e-26, 0.19854317637027171e-26, 0.11552087612340893e-26] # These are from phd_87a_lim_get_lim_alm.py
    fn.append(tmp)
    co.append(10*np.ones(len(tmp)))
    
    # SINFONI
    tmp = cc*1e8/np.array([23275., 21300., 17512.5, 15475.])
    nu.append(tmp)
#    Final bin lims flx [  2.43663528e-19   1.39144179e-19   1.82898534e-19   1.20305289e-19]
    tmp = 1e-19*np.array([1.2030, 1.8289, 1.3914, 2.4366])*np.array([23275., 21300., 17512.5, 15475.])/tmp
    fn.append(tmp)
    co.append(11*np.ones(len(tmp)))
    
    # STIS
    #tmp = cc/(1e-8*np.linspace(5300, 10000, NN))
    #nu.append(tmp)
    #tmp = 7.80191789426e-15*np.linspace(5300, 10000, NN)**1.05/(cc*1e8)
    #fn.append(tmp)
    #co.append(12*np.ones(len(tmp)))
    tmp = np.logspace(np.log10(cc/5300e-8), np.log10(cc/10000e-8), NN)
    nu.append(tmp)
    tmp = np.logspace(np.log10(7.4823e-15*5300**(2-0.9527)/(1e8*cc)), np.log10(7.4823e-15*10000**(2-0.9527)/(1e8*cc)), NN)
    fn.append(tmp)
    co.append(12*np.ones(len(tmp)))
    
    
    # WFC3/UVIS
    tmp = cc*1e8/np.array([6255.5, 4330.5, 8074., 5334., 3359.5, 2382.5])
    nu.append(tmp)
    tmp = 1e-18*np.array([5.5, 7.1, 2.8, 6.7, 9.5, 40.])*np.array([6255.5, 4330.5, 8074., 5334., 3359.5, 2382.5])/tmp
    fn.append(tmp)
    co.append(13*np.ones(len(tmp)))
    
    # Chandra ACIS, 3.32580E-04 for Gamma = 1.63, 7.45166E-04 for Gamma = 2.108
#    tmp = np.logspace(np.log10(4.836e17), np.log10(24.18e17), NN)
#    nu.append(tmp)
#    tmp = 10**-13.0817
#    tmp = tmp/(2*(np.sqrt(2.418e18)-np.sqrt(4.836e17)))
#    tmp = tmp*np.array(nu[-1])**-0.5
#    fn.append(tmp)
#    co.append(8*np.ones(len(tmp)))

    tmp = np.logspace(np.log10(4.836e17), np.log10(24.18e17), NN)
    nu.append(tmp)
    tmp = np.logspace(np.log10(2), 1, NN)
    tmp = 3.45466e-04*tmp**-0.63*kev2erg
    fn.append(tmp/kev2hz)
    co.append(8*np.ones(len(tmp)))
    
    tmp = np.logspace(np.log10(4.836e17), np.log10(24.18e17), NN)
    nu.append(tmp)
    tmp = np.logspace(np.log10(2), 1, NN)
    tmp = 7.74668e-04*tmp**-1.108*kev2erg
    fn.append(tmp/kev2hz)
    co.append(14*np.ones(len(tmp)))

    ########
    # finalize
    nu = np.concatenate(nu)
    fn = np.concatenate(fn)
    co = np.concatenate(co)
    return nu, fn, co



################################################################
# Overplotting
crp_nu, crp_nFn = get_crab_pulsar()
crp_nFn =  crp_nFn * (2.2/distance)**2

crn_nu, crn_nFn = get_crab_nebula()
crn_nFn =  crn_nFn * (2.2/distance)**2

cas_nu, cas_nFn, cas_xx, cas_yyy = get_cco_cas()
cas_nFn =  cas_nFn * (3.4/distance)**2
cas_yyy =  cas_yyy * (3.4/distance)**2

nu, fn, co = get_all_lim()
ene_bud_lim = get_ene_bud_lim()
ene_bud_lim_df = get_ene_bud_lim_df()
print('A solar luminosity is', Lsun)
print('The energy budget estimate of the CO contribution is', ene_bud_lim[0]/Lsun, 'with a 1 std uncertainty of', ene_bud_lim[1]/Lsun)

########
# Blackbody
RR = 1.2e6
MM = 1.4
RS = 2*GG*MM*Msun/cc**2
gr = np.sqrt(1-RS/RR)
TT = (ene_bud_lim[2]/(gr**2*4*np.pi*RR**2*sb))**0.25
TT_df = (ene_bud_lim_df/(gr**2*4*np.pi*RR**2*sb))**0.25

########
# Accretion
acc_coe = 0.1*cc**2*Msun/(365.25*24*3600)
print('Accretion limited to', ene_bud_lim[2]/acc_coe, 'M_Sun yr-1.')
edd_coe = 4*np.pi*GG*Msun*mp*cc/sT

    
smooth = 5
cra_nu = np.logspace(8, 25.9, 300)
crp_nFn2 = 10**griddata(np.log10(crp_nu), np.log10(crp_nFn), np.log10(cra_nu), method='linear')
crp_nFn2 = np.convolve(crp_nFn2, np.ones(smooth)/smooth, 'same')
crn_nFn2 = 10**griddata(np.log10(crn_nu), np.log10(crn_nFn), np.log10(cra_nu), method='linear')
crn_nFn2 = np.convolve(crn_nFn2, np.ones(smooth)/smooth, 'same')
crn_nFn2 += crp_nFn2
lum_crn = np.trapz(crn_nFn2/cra_nu, cra_nu)*4*np.pi*(distance*kpc)**2
lum_crp = np.trapz(crp_nFn2/cra_nu, cra_nu)*4*np.pi*(distance*kpc)**2



################################################################
# Initiate plot
fig = plt.figure(figsize=(10, 3.75))
ax1 = fig.add_subplot(111)

# Plot limits
tmp_ind = (co==10) | (co==11) | (co==12) | (co==13) | (co==14) | (co==8)
ax1.loglog(nu[tmp_ind], nu[tmp_ind]*fn[tmp_ind], 'v', mew=0)
tmp_ind = (co==0) | (co==6) | (co==7) | (co==9) | (co==16) | (co==15) | (co==17) | (co==18) | (co==19)
ax1.loglog(nu[tmp_ind], nu[tmp_ind]*fn[tmp_ind], 'v', mew=0)
tmp_ind = (co==2) | (co==4) | (co==5) | (co==3)
ax1.loglog(nu[tmp_ind], nu[tmp_ind]*fn[tmp_ind], 'v', alpha=0.6, mew=0, zorder=-2)

# Reconstruct the crabs
smooth = 5
cra_nu2 = np.logspace(7, 29, 500)
crp_nFn4 = 10**griddata(np.log10(crp_nu), np.log10(crp_nFn), np.log10(cra_nu2), method='linear')
crn_nFn4 = 10**griddata(np.log10(crn_nu), np.log10(crn_nFn), np.log10(cra_nu2), method='linear')
# Cut some trash
pind = -np.isnan(crp_nFn4)
nind = -np.isnan(crn_nFn4)
crp_nFn4 = crp_nFn4[pind]
crp_nu2 = cra_nu2[pind]
crn_nFn4[pind] += crp_nFn4
crn_nFn4 = crn_nFn4[nind]
crn_nu2 = cra_nu2[nind]

crp_nFn4 = np.convolve(crp_nFn4, np.ones(smooth)/smooth, 'same')
crn_nFn4 = np.convolve(crn_nFn4, np.ones(smooth)/smooth, 'same')
crp_nFn4 = crp_nFn4[int(smooth/2)+1:-int(smooth/2)-1]
crp_nu2 = crp_nu2[int(smooth/2)+1:-int(smooth/2)-1]
crn_nFn4 = crn_nFn4[int(smooth/2)+1:-int(smooth/2)-1]
crn_nu2 = crn_nu2[int(smooth/2)+1:-int(smooth/2)-1]

crp_cut1 = (crp_nu2 < 2e9)
crp_cut2 = (crp_nu2 > 1e14)

# Plot overplots
#ax1.loglog(crp_nu, crp_nFn, 'o', zorder=-190, mew=0, ms=3)
#ax1.loglog(crn_nu, crn_nFn, 'o', zorder=-190, mew=0, ms=3)

ax1.loglog(crp_nu2[crp_cut1], crp_nFn4[crp_cut1], ls='--', zorder=-200)
ax1.loglog(crp_nu2[crp_cut2], crp_nFn4[crp_cut2], ls='--', color='#d62728', zorder=-200)
ax1.loglog(crn_nu2, crn_nFn4, zorder=-200, ls='-.')
ax1.loglog(cas_xx, cas_yyy, zorder=-200)

# Labels and limits
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("$\\nu F_\\nu$ (Hz~erg s$^{-1}$ cm$^{-2}$~Hz$^{-1}$)")
ax1.set_xlim([  1e7,  1e29])
ax1.set_ylim([1e-19, 1e-10])

# Upper x-axis
up_x = np.logspace(-21, 21, 22)
ax2 = ax1.twiny()
ax2.set_xscale(ax1.get_xscale())
ax2.set_xticks(cc/(up_x*1e-8))
ax2.set_xticklabels(tic_hlp(up_x))
ax2.set_xlim(ax1.get_xlim())
ax2.minorticks_off()
ax2.set_xlabel("Wavelength (\AA{})")

# Right y-axis
ri_y = np.logspace(10, 40, 31)
ax3 = ax1.twinx()
ax3.set_yscale(ax1.get_yscale())
ax3.set_yticks(ri_y/(4*np.pi*DD**2))
ax3.set_yticklabels(tic_hlp(ri_y))
ax3.set_ylim(ax1.get_ylim())
ax3.minorticks_off()
ax3.set_ylabel("$\\nu L_\\nu$ (Hz~erg s$^{-1}$~Hz$^{-1}$)")

# Annotations
ax1.annotate('Crab Pulsar', xy=(1e22, 3e-14), color='#d62728')
ax1.annotate('Crab Nebula', xy=(1e26, 1e-12), color='#9464bd')
ax1.annotate('Cas A CCO', xy=(3e15, 1e-17), color='#8c564b')
fig.savefig('/Users/silver/box/phd/pro/87a/lim/art/figs/ene_bud.pdf',bbox_inches='tight', pad_inches=0.03)


################################################################
# Model limits
#fig = plt.figure(figsize=(10, 3.75))
fig = plt.figure(figsize=(5, 2.4))
ax1 = fig.add_subplot(111)

# Plot limits
tmp_ind = (co==10) | (co==11) | (co==12) | (co==13) | (co==14) | (co==8) 
ax1.loglog(nu[tmp_ind], nu[tmp_ind]*fn[tmp_ind], 's', mew=0)
#tmp_ind = (co==0) | (co==6) | (co==7) | (co==9) | (co==16) | (co==15)
#ax1.loglog(nu[tmp_ind], nu[tmp_ind]*fn[tmp_ind], 's', mew=0)
#tmp_ind = (co==2) | (co==3) | (co==4) | (co==5)
#ax1.loglog(nu[tmp_ind], nu[tmp_ind]*fn[tmp_ind], 's', alpha=0.7, mew=0)

########
# Find limits on Crab for the given limits
#ax1.loglog(cra_nu, crp_nFn2, zorder=-200)
#ax1.loglog(cra_nu, crn_nFn2, zorder=-200)
#ax1.loglog(cas_xx, cas_yyy, zorder=-200)

########
# These are the observed values at the position we have limits
crp_nFn3 = griddata(cra_nu, crp_nFn2, nu, method='linear')
crn_nFn3 = griddata(cra_nu, crn_nFn2, nu, method='linear')
crp_frac = nu*fn/crp_nFn3
crn_frac = nu*fn/crn_nFn3
crp_frac = np.where(np.isnan(crp_frac), 1e99, crp_frac)
crn_frac = np.where(np.isnan(crn_frac), 1e99, crn_frac)
crp_lim_idx = np.argmin(crp_frac)
crn_lim_idx = np.argmin(crn_frac)
crp_lim_fra = crp_frac[crp_lim_idx]
crn_lim_fra = crn_frac[crn_lim_idx]
lum_crp_lim = lum_crp*crp_lim_fra
lum_crn_lim = lum_crn*crn_lim_fra
wav_crp_lim = nu2lam(nu[crp_lim_idx])*1e4
wav_crn_lim = nu2lam(nu[crn_lim_idx])*1e4
print('The limit for the Crab pulsar is L =', lum_crp_lim/Lsun, 'Lsun,', lum_crp_lim, 'at a wavelength of', wav_crp_lim, 'microns, which is', crp_lim_fra, 'of the real Crab pulsar.')
print('The limit for the Crab nebula is L =', lum_crn_lim/Lsun, 'Lsun,', lum_crn_lim, 'at a wavelength of', wav_crn_lim, 'microns, which is', crn_lim_fra, 'of the real Crab nebula.')

crp_ind = (co==8)
xcrp_nFn3 = griddata(cra_nu, crp_nFn2, nu[crp_ind], method='linear')
xcrp_frac = nu[crp_ind]*fn[crp_ind]/xcrp_nFn3
xcrp_frac = np.where(np.isnan(xcrp_frac), 1e99, xcrp_frac)
xcrp_lim_idx = np.argmin(xcrp_frac)
xcrp_lim_fra = xcrp_frac[xcrp_lim_idx]
xlum_crp_lim = lum_crp*xcrp_lim_fra
xkev_crp_lim = nu[crp_ind][xcrp_lim_idx]/kev2hz
print('The limit for the Crab pulsar is L =', xlum_crp_lim/Lsun, 'Lsun,', xlum_crp_lim, 'at an energy of', xkev_crp_lim, 'keV, which is', xcrp_lim_fra, 'of the real Crab pulsar.')

crn_ind = (co==10)
xcrn_nFn3 = griddata(cra_nu, crn_nFn2, nu[crn_ind], method='linear')
xcrn_frac = nu[crn_ind]*fn[crn_ind]/xcrn_nFn3
xcrn_frac = np.where(np.isnan(xcrn_frac), 1e99, xcrn_frac)
xcrn_lim_idx = np.argmin(xcrn_frac)
xcrn_lim_fra = xcrn_frac[xcrn_lim_idx]
xlum_crn_lim = lum_crn*xcrn_lim_fra
print('The limit for the Crab nebula is L =', xlum_crn_lim/Lsun, 'Lsun,', xlum_crn_lim, 'at a frequency of', nu[crn_ind][xcrn_lim_idx]/1e9, 'GHz, which is', xcrn_lim_fra, 'of the real Crab nebula.')

################################################################
# Revised version that accounts for dust using a covering factor, reducing the number of cases from 4 to 2
uvoir = (co==11) | (co==12) | (co==13)
fdust = 0.6
cor_nfn = nu*fn
cor_nfn = cor_nfn[uvoir]/(1-fdust)

dcrn_nFn3 = griddata(cra_nu, crn_nFn2, nu[uvoir], method='linear')
dcrn_frac = cor_nfn/dcrn_nFn3
dcrn_frac = np.where(np.isnan(dcrn_frac), 1e99, dcrn_frac)
dcrn_lim_idx = np.argmin(dcrn_frac)
dcrn_lim_fra = dcrn_frac[dcrn_lim_idx]
dlum_crn_lim = lum_crn*dcrn_lim_fra
dwav_crn_lim = nu2lam(nu[uvoir][dcrn_lim_idx])*1e4
print('The limit for the Crab nebula is L =', dlum_crn_lim/Lsun, 'Lsun,', dlum_crn_lim, 'at a wavelength of', dwav_crn_lim, 'microns, which is', dcrn_lim_fra, 'of the real Crab nebula.')


#har_ind = (co==8)
#
#dcrp_frac = cor_nfn/crp_nFn3
#dlum_crn_lim = lum_crn*dcrn_lim_fra
#dcrp_frac = np.where(np.isnan(dcrp_frac), 1e99, dcrp_frac)
#dcrp_lim_idx = np.argmin(dcrp_frac)
#dcrp_lim_fra = dcrp_frac[dcrp_lim_idx]
#dlum_crp_lim = lum_crp*dcrp_lim_fra
#dwav_crp_lim = nu[dcrp_lim_idx]/kev2hz
#
#har_ind = (co==8)
#dcrn_nFn3 = griddata(cra_nu, crn_nFn2, nu[har_ind], method='linear')
#dcrn_frac = nu[har_ind]*fn[har_ind]/dcrn_nFn3
#dcrn_frac = cor_nfn/dcrn_nFn3
#dcrn_frac = np.where(np.isnan(dcrn_frac), 1e99, dcrn_frac)
#dcrn_lim_idx = np.argmin(dcrn_frac)
#dcrn_lim_fra = dcrn_frac[dcrn_lim_idx]
#dwav_crn_lim = nu2lam(nu[dcrn_lim_idx])*1e4

# Plot the limited Crabs
ax1.loglog(cra_nu, crp_lim_fra*crp_nFn2, zorder=-200)
ax1.loglog(cra_nu, crn_lim_fra*crn_nFn2, zorder=-200)
ax1.loglog(cra_nu, xcrp_lim_fra*crp_nFn2, zorder=-200)
ax1.loglog(cra_nu, xcrn_lim_fra*crn_nFn2, zorder=-200)
ax1.loglog(cra_nu, dcrn_lim_fra*crn_nFn2, zorder=-200)

# Labels and limits
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("$\\nu F_\\nu$ (Hz~erg s$^{-1}$ cm$^{-2}$~Hz$^{-1}$)")
ax1.set_xlim([  1e7,  1e29])
ax1.set_ylim([1e-22, 1e-10])
ax1.set_xlim([  1e8,  1e26])
ax1.set_ylim([1e-17, 1e-10])

# Upper x-axis
up_x = np.logspace(-21, 21, 22)
ax2 = ax1.twiny()
ax2.set_xscale(ax1.get_xscale())
ax2.set_xticks(cc/(up_x*1e-8))
ax2.set_xticklabels(tic_hlp(up_x))
ax2.set_xlim(ax1.get_xlim())
ax2.minorticks_off()
ax2.set_xlabel("Wavelength (\AA{})")

# Right y-axis
ri_y = np.logspace(10, 40, 31)
ax3 = ax1.twinx()
ax3.set_yscale(ax1.get_yscale())
ax3.set_yticks(ri_y/(4*np.pi*DD**2))
ax3.set_yticklabels(tic_hlp(ri_y))
ax3.set_ylim(ax1.get_ylim())
ax3.minorticks_off()
ax3.set_ylabel("$\\nu L_\\nu$ (Hz~erg s$^{-1}$~Hz$^{-1}$)")

fig.savefig('/Users/silver/box/phd/pro/87a/lim/art/figs/ene_bud_lim.pdf',bbox_inches='tight', pad_inches=0.03)



################################################################
# BP-plane
manchester05 = np.loadtxt('/Users/silver/box/sci/lib/m/manchester05_p0_bsurf.dat')
P0 = manchester05[:,0]
Bsurf = manchester05[:,1]

RR = 1e6
PP, BB = np.meshgrid(np.logspace(-3, 1.5, 100), np.logspace(7, 17, 100))
Edot = 2**5*BB**2*RR**6*np.pi**4/(3*cc**3*PP**4)
xx = np.logspace(-3, 1.5, 100)
Blim = np.sqrt(ene_bud_lim[2]*3*cc**3/(2**5*RR**6*np.pi**4))*xx**2
fig = plt.figure(figsize=(5, 3.75))
#fig = plt.figure(figsize=(2.4, 2.4))
ax = plt.gca()
cf = ax.pcolor(np.log10(PP), np.log10(BB), np.log10(Edot), rasterized=True, alpha=0.)

# Plot the regions limited by different luminosities
levels = ene_bud_lim[0]+ene_bud_lim[1]*np.arange(1,4)
levels = ene_bud_lim[2]*np.arange(1,2)
#plt.contour(np.log10(PP), np.log10(BB), np.log10(Edot), np.log10(levels), linewidths=2, colors='#1f77b4', zorder=100)
plt.contour(np.log10(PP), np.log10(BB), np.log10(Edot), np.log10([xlum_crp_lim]), linewidths=2, colors='#ffffff', linestyles='dashed', zorder=100)
plt.contour(np.log10(PP), np.log10(BB), np.log10(Edot), np.log10([lum_crn_lim]), linewidths=2, colors='#ffffff', linestyles='solid', zorder=100)
#plt.plot(np.log10(xx), np.log10(Blim), 'w')
#plt.contour(np.log10(PP), np.log10(BB), np.log10(Edot), np.log10([lum_crp_lim]), linewidths=2, colors='#d62728') # red
#plt.contour(np.log10(PP), np.log10(BB), np.log10(Edot), np.log10([lum_crn_lim]), linewidths=2, colors='#ff7f0e') # orange
#plt.contour(np.log10(PP), np.log10(BB), np.log10(Edot), np.log10([xlum_crp_lim]), linewidths=2, colors='#e377c2') # pink
#plt.contour(np.log10(PP), np.log10(BB), np.log10(Edot), np.log10([xlum_crn_lim]), linewidths=2, colors='#8c564b') # brown

# Plot the full pulsar catalogue along with objects of special interes
plt.plot(np.log10(P0), np.log10(Bsurf), '.', ms=1, color='#000000')
plt.plot(np.log10(0.0333924123), np.log10(3.79e+12), 'd', color='#ffffff', ms=6)

# Extrapolate in time
# Characteristic time
tt = 30*365.26*24*3600
tp =  1.4*Msun*3*cc**3/(40*np.pi**2*RR**4)*PP**2/BB**2
# Luminosity
LL = 2**5*np.pi**4*RR**6/(3*cc**3)*(1+tt/tp)**-2*BB**2/PP**4

# This is just a horrible way of plotting lines
cf = ax.pcolor(np.log10(PP), np.log10(BB), np.log10(LL), rasterized=True, alpha=1)
cb = fig.colorbar(cf)
#plt.contour(np.log10(PP), np.log10(BB), np.log10(LL), np.log10(levels), linewidths=2, colors='#1f77b4', zorder=100, alpha=0.5)
plt.contour(np.log10(PP), np.log10(BB), np.log10(LL), np.log10([xlum_crp_lim]), linewidths=2, colors='#ffffff', linestyles='dashed', zorder=99, alpha=0.7)
plt.contour(np.log10(PP), np.log10(BB), np.log10(LL), np.log10([lum_crn_lim]), linewidths=2, colors='#ffffff', linestyles='solid', zorder=99, alpha=0.7)

plt.gca().annotate('Excluded', xy=(-2.7, 14))
plt.gca().annotate('Allowed', xy=(0, 9), color='#ffffff')
plt.xlabel('log$(P)$ (s)')
plt.ylabel('log$(B)$ (G)')
cb.set_label('log$(L)$ (erg s$^{-1}$)')
fig.savefig('/Users/silver/box/phd/pro/87a/lim/art/figs/bp.pdf',bbox_inches='tight', pad_inches=0.03)

print('Here comes the C_BP values:')
print(np.sqrt(lum_crp_lim*3*cc**3/(2**5*RR**6*np.pi**4)))
print(np.sqrt(lum_crn_lim*3*cc**3/(2**5*RR**6*np.pi**4)))
print(np.sqrt(xlum_crp_lim*3*cc**3/(2**5*RR**6*np.pi**4)))
print(np.sqrt(xlum_crn_lim*3*cc**3/(2**5*RR**6*np.pi**4)))
print(np.sqrt(dlum_crn_lim*3*cc**3/(2**5*RR**6*np.pi**4)))

print('Dust-free bolometric limit C_BP:', np.sqrt(ene_bud_lim_df/(2**5*RR**6*np.pi**4/(3*cc**3))))
print('Dust-obscured bolometric limit C_BP:', np.sqrt(ene_bud_lim[2]/(2**5*RR**6*np.pi**4/(3*cc**3))))

pdb.set_trace()
plt.show()
