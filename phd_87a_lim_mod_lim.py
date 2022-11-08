# Dennis Alp 2016-01-27
# Find models that are consistent with oberved limits.
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/Dropbox/bin/phd_87a_lim_plt_lim.py

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

################################################################
#For LaTeX style font in plots
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)

################################################################
# Parameters
# Logistics
WRK_DIR = "/Users/silver/box/phd/pro/87a/lim/glo"
#os.chdir(WRK_DIR) #Move to designated directory
NN = 8

distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e18 # cm
kpc = 1e3*pc # cm
DD = distance*kpc
Lsun = 3.826e33 # erg s-1
Msun = 1.989e33 # g
cc = 29979245800 # cm s-1
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

def get_crab():
    # radio: lorimer95
    # nir: tziamtzis09
    # uvo: sollerman00
    # xray: kuiper01
    nu = np.array([408e6, 606e6, 925e6, 1408e6, 10**14.147, 10**14.266, 10**14.575, 10**14.5, 10**14.8, 10**15.1, 10**15.4, 10**16.883411, 10**17.383411, 10**17.883411, 10**18.383411, 10**18.883411, 10**19.383411])
    ff = np.array([646e-26, 210e-26, 45.1e-26, 14.4e-26, 10**-25.55, 10**-25.55, 10**-25.611, 10**-25.53, 10**-25.5, 10**-25.42, 10**-25.4, 4e-10/10**16.883411, 7e-10/10**17.383411, 1.1e-9/10**17.883411, 1.3e-9/10**18.383411, 1.3e-9/10**18.883411, 1.3e-9/10**19.383411])
    return nu, nu*ff

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

def get_geminga(): # Geminga, danilenko11
    dat = np.loadtxt('/Users/silver/box/sci/lib/d/danilenko11_geminga.txt')
    nu = dat[:,0]
    fnu = dat[:,1]*1e-29
    return nu, nu*fnu

def get_vela(): # Vela, danilenko11
    dat = np.loadtxt('/Users/silver/box/sci/lib/d/danilenko11_vela.txt')
    nu = dat[:,0]
    fnu = dat[:,1]*1e-29
    return nu, nu*fnu

def get_cco_g15(): # SNR G15.9+0.2, CXOU J1818, klochkov16, 2 MK
    dat = np.loadtxt('/Users/silver/box/sci/lib/k/klochkov16.txt')
    ee = dat[:,0]
    emod = dat[:,1]
    flx = np.loadtxt('/Users/silver/box/sci/lib/k/klochkov16.flx')
    ene = np.loadtxt('/Users/silver/box/sci/lib/k/klochkov16.ene')
    xx, yy = from_xspec_mod(ene, flx)
    nu, nufnu = emod2nufnu(ee, emod)
#    return np.r_[nu, tm], np.r_[nufnu, tmp]
    return nu, nufnu, xx, yy

def get_cco_hes(): # HESS J1731-347, klochkov15, exceptionally hot for the estimated age of ~30 kyr, 2.24 MK
    dat = np.loadtxt('/Users/silver/box/sci/lib/k/klochkov15.txt')
    ee = dat[:,0]
    emod = dat[:,1]
    flx = np.loadtxt('/Users/silver/box/sci/lib/k/klochkov15.flx')
    ene = np.loadtxt('/Users/silver/box/sci/lib/k/klochkov15.ene')
    xx, yy = from_xspec_mod(ene, flx)
    nu, nufnu = emod2nufnu(ee, emod)
    return nu, nufnu, xx, yy

def get_cco_cas(): # Cas A, posselt13
    dat = np.loadtxt('/Users/silver/box/sci/lib/p/posselt13.txt')
    ee = dat[:,0]
    emod = dat[:,1]
    flx = np.loadtxt('/Users/silver/box/sci/lib/p/posselt13.flx')
    ene = np.loadtxt('/Users/silver/box/sci/lib/p/posselt13.ene')
    xx, yy = from_xspec_mod(ene, flx)
    nu, nufnu = emod2nufnu(ee, emod)
    return nu, nufnu, xx, yy

def get_cco_kes(): # PSR J1852+0040 in Kes 79, bogdanov14
    dat = np.loadtxt('/Users/silver/box/sci/lib/b/bogdanov14.txt')
    ee = dat[:,0]
    emod = dat[:,1]
    return emod2nufnu(ee, emod)

def get_lmc_pulsar(): # PSR J0537-6910 in LMC, hess12
    dat = np.loadtxt('/Users/silver/box/sci/lib/h/hess12.txt')
    nu = dat[:,0]*1e9*kev2hz
    nufnu = dat[:,1]
    return nu, nufnu

def get_g3106_pulsar(): # J1400-6325 in G310.6-1.6, martin14
    dat = np.loadtxt('/Users/silver/box/sci/lib/m/martin14_g3106-16.txt')
    nu = dat[:,0]
    nufnu = dat[:,1]
    return nu, nufnu

def get_n157b_pulsar(): # N157B in the LMC, martin14
    dat = np.loadtxt('/Users/silver/box/sci/lib/m/martin14_n157b.txt')
    nu = dat[:,0]
    nufnu = dat[:,1]
    return nu, nufnu

def get_g292_pulsar(): # J1124-5916 in G292.0+0.18, martin14
    dat = np.loadtxt('/Users/silver/box/sci/lib/m/martin14_g2920+18.txt')
    nu = dat[:,0]
    nufnu = dat[:,1]
    return nu, nufnu

def get_g769_pulsar(): # J2022+3842 in G76.9+1.0, martin14
    dat = np.loadtxt('/Users/silver/box/sci/lib/m/martin14_g769+10.txt')
    nu = dat[:,0]
    nufnu = dat[:,1]
    return nu, nufnu

def get_b1509_pulsar(): # PSR B1509-58, abdo10c
    dat = np.loadtxt('/Users/silver/box/sci/lib/a/abdo10c.txt')
    ee = dat[:,0]/1000.*kev2hz
    eemod = dat[:,1]
    return ee, eemod

def get_xsp():
    flx = np.loadtxt('/Users/silver/box/phd/pro/87a/lim/cha/16756/repro/spectra/bbo.flx')
    ene = np.loadtxt('/Users/silver/box/phd/pro/87a/lim/cha/16756/repro/spectra/bbo.ene')
    xx, yy = from_xspec_mod(ene, flx)
    return xx, yy

################################################################
# Data
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

## zanardo14, estimate, not understandable, subtracted, modeled, something bogus, MEM image
#tmp = 1e9*np.linspace(102, 672, NN)
#nu.append(tmp)
#tmp = 3e-26*np.ones(NN)
#fn.append(tmp)
#co.append(1*np.ones(len(tmp)))

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
co.append(3*np.ones(len(tmp)))

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
tmp = [1.4e-26]
fn.append(tmp)
co.append(5*np.ones(len(tmp)))

# matsuura15, entire remnant, dust except 500 microns, which is 3sigma limit
tmp = lam2nu(1e-4*np.array([70, 100, 160, 250, 350, 500]))
nu.append(tmp)
tmp = 1e-26*np.array([45.4, 82.4, 153., 110.7, 69.3, 60.])
fn.append(tmp)
co.append(6*np.ones(len(tmp)))

# arendt16, entire remnant
tmp = lam2nu(1e-4*np.array([3.6, 4.5, 5.8, 8., 24., 3.6, 4.5, 5.8, 8., 24.]))
nu.append(tmp)
tmp = 1e-26*np.array([1.52, 2.17, 4.08, 13.61, 75.7, 0.99, 1.2, 1.62, 4.69, 26.3])
fn.append(tmp)
co.append(7*np.ones(len(tmp)))

# ng09, 90% limit
tmp = np.logspace(np.log10(4.836e17), np.log10(24.18e17), NN)
nu.append(tmp)
tmp = 3.3e34/(4*np.pi*(51.2*kpc)**2)
tmp = tmp/(2*(np.sqrt(2.418e18)-np.sqrt(4.836e17)))
tmp = tmp*np.array(nu[-1])**-0.5
fn.append(tmp)
co.append(14*np.ones(len(tmp)))

# frank16, 90% limit
tmp = np.logspace(np.log10(4.836e17), np.log10(24.18e17), NN)
nu.append(tmp)
tmp = 1.5e34/(4*np.pi*(51.2*kpc)**2)
#tmp = 3.1e33/(4*np.pi*(51.2*kpc)**2)
tmp = tmp/(2*(np.sqrt(2.418e18)-np.sqrt(4.836e17)))
tmp = tmp*np.array(nu[-1])**-0.5
fn.append(tmp)
co.append(8*np.ones(len(tmp)))

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
tmp = [3.28467956e-26, 4.78012075e-26, 2.29446108e-26] # These are from phd_87a_lim_plt_lim.py
tmp = [4.49094e-26, 3.16767e-26, 2.46548e-26] # These are from phd_87a_lim_plt_lim.py
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
tmp = np.logspace(np.log10(7.80191789426e-15*5300**1.05/(cc*1e8)), np.log10(7.80191789426e-15*10000**1.05/(cc*1e8)), NN)
fn.append(tmp)
co.append(12*np.ones(len(tmp)))


# WFC3/UVIS
tmp = cc*1e8/np.array([6255.5, 4330.5, 8074., 5334., 3359.5, 2382.5])
nu.append(tmp)
tmp = 1e-18*np.array([5.5, 7.1, 2.8, 6.7, 9.5, 40.])*np.array([6255.5, 4330.5, 8074., 5334., 3359.5, 2382.5])/tmp
fn.append(tmp)
co.append(13*np.ones(len(tmp)))

# mine, for comparison with frank16, model: con_pha_veq_vps_bkn_con_pha_cfl_pow, 1e-13.0817
tmp = np.logspace(np.log10(4.836e17), np.log10(24.18e17), NN)
nu.append(tmp)
tmp = 10**-13.0817
tmp = tmp/(2*(np.sqrt(2.418e18)-np.sqrt(4.836e17)))
tmp = tmp*np.array(nu[-1])**-0.5
fn.append(tmp)
co.append(17*np.ones(len(tmp)))

########
# finalize
nu = np.concatenate(nu)
fn = np.concatenate(fn)
co = np.concatenate(co)



################################################################
# Overplotting
crp_nu, crp_nFn = get_crab_pulsar()
crp_nFn =  crp_nFn * (2.2/distance)**2

crn_nu, crn_nFn = get_crab_nebula()
crn_nFn =  crn_nFn * (2.2/distance)**2

g15_nu, g15_nFn, g15_xx, g15_yyy = get_cco_g15()
g15_nFn =  g15_nFn * (10/distance)**2
g15_yyy =  g15_yyy * (10/distance)**2

hes_nu, hes_nFn, hes_xx, hes_yyy = get_cco_hes()
hes_nFn =  hes_nFn * (3.2/distance)**2
hes_yyy =  hes_yyy * (3.2/distance)**2

cas_nu, cas_nFn, cas_xx, cas_yyy = get_cco_cas()
cas_nFn =  cas_nFn * (3.4/distance)**2
cas_yyy =  cas_yyy * (3.4/distance)**2

kes_nu, kes_nFn = get_cco_kes()
kes_nFn =  kes_nFn * (7.1/distance)**2

gem_nu, gem_nFn = get_geminga()
gem_nFn =  gem_nFn * (0.25/distance)**2

vel_nu, vel_nFn = get_vela()
vel_nFn =  vel_nFn * (0.287/distance)**2

lmc_nu, lmc_nFn = get_lmc_pulsar()
lmc_nFn =  lmc_nFn * (48.1/distance)**2

g31_nu, g31_nFn = get_g3106_pulsar()
g31_nFn =  g31_nFn * (7/distance)**2

n15_nu, n15_nFn = get_n157b_pulsar()
n15_nFn =  n15_nFn * (48/distance)**2

g29_nu, g29_nFn = get_g292_pulsar()
g29_nFn =  g29_nFn * (6/distance)**2

g76_nu, g76_nFn = get_g769_pulsar()
g76_nFn =  g76_nFn * (10/distance)**2

b15_nu, b15_nFn = get_b1509_pulsar()
b15_nFn =  b15_nFn * (5.2/distance)**2

xsp_xx, xsp_yyy = get_xsp()


################################################################
# Modelling part
RR = 1.3055230790582144*1e6
TT = 3459030.3939
alpha_nu = np.log10(1e6*np.array([408, 606, 925, 1408]))
alpha_ss = np.log10(1e-26*np.array([646, 210, 45.1, 14.1]))
alpha_cc = np.polyfit(alpha_nu, alpha_ss, 1)
alpha = alpha_cc[0] # S(nu) = nu^alpha, lorimer95
Gamma = 1.63 # N(E) = E^-Gamma, willingale01
#RR = 5.5e15
#TT = 140
mod_nu = np.logspace(7, 20, 1024)

# Blackbody
pro_are = 4*np.pi*RR**2 # Surface area
sol_ang = np.pi # Solid angle
mod_bbF  = mod_nu*pro_are*sol_ang*Bnu(mod_nu, TT)/(4*np.pi*DD**2)
mod_bbF2 = mod_nu*pro_are*sol_ang*Bnu(mod_nu, TT/2)/(4*np.pi*DD**2)

# Radio, limit
mod_nuS = mod_nu**(alpha+1)
mod_nuS = np.where((mod_nu < 1.606e9) & (mod_nu > 0.408e9), mod_nuS, 0) # < 1.606 GHz
mod_nuS = get_lim(mod_nuS)
# Crab
mod_nuS = 10**alpha_cc[1]*mod_nu**(alpha+1)
mod_nuS = np.where((mod_nu < 1.606e9) & (mod_nu > 0.408e9), mod_nuS, 0) #  < 1.606 GHz
mod_nuS = mod_nuS * (2.2/distance)**2 # 2.2 kpc fom manchester05

# X-ray
mod_eeF = mod_nu**(-Gamma+2)
mod_eeF = np.where(mod_nu > 2.418e16, mod_eeF, 0) # > 0.1 keV
mod_eeF = get_lim(mod_eeF)

## GR blackbody, just an exercise lol
#MM = 1.5*Msun
#RS = 2*GG*MM/cc**2
#gr = np.sqrt(1-RS/RR)
#T2 = TT/gr # Sort of gravitational well escape
#R2 = RR*gr # Lenth contraction
#pro_ar2 = np.pi*R2**2
#mod_bb2 = mod_nu*gr*pro_ar2*sol_ang*Bnu(mod_nu, T2)/(4*np.pi*DD**2)
#mod_nu2 = mod_nu*gr # Photons escaping a well
#mod_bb2 = mod_bb2*gr # Time dilation
#plt.loglog(mod_nu,mod_bbF)
#plt.loglog(mod_nu2,mod_bb2)
#plt.show()


################################################################
# Initiate plot
fig = plt.figure(figsize=(10, 3.75))
ax1 = fig.add_subplot(111)

# Plot limits
tmp_ind = (co==0) | (co==6) | (co==7) | (co==9) | (co==16) | (co==15)
ax1.loglog(nu[tmp_ind], nu[tmp_ind]*fn[tmp_ind], 'bs', mew=0)
tmp_ind = (co==2) | (co==3) | (co==4) | (co==5) | (co==8)
ax1.loglog(nu[tmp_ind], nu[tmp_ind]*fn[tmp_ind], 'gs', alpha=0.3, mew=0)
tmp_ind = (co==10) | (co==11) | (co==12) | (co==13) | (co==14) | (co==17)
ax1.loglog(nu[tmp_ind], nu[tmp_ind]*fn[tmp_ind], 'rs', mew=0)

# Plot models
#ax1.loglog(mod_nu, mod_bbF, 'k', alpha=0.3)
#ax1.loglog(mod_nu, mod_bbF2, 'k')
#ax1.loglog(mod_nu, mod_nuS, 'm')
#ax1.loglog(mod_nu, mod_eeF, 'y')

## Plot overplots
ax1.loglog(crp_nu, crp_nFn, 'co', zorder=-190, mew=0, ms=3)
ax1.loglog(crn_nu, crn_nFn, 'mo', zorder=-190, mew=0, ms=3)
#ax1.loglog(lmc_nu, lmc_nFn, 'yo', zorder=-100, mew=0, ms=3)
#ax1.loglog(b15_nu, b15_nFn, 'ko', zorder=-100, mew=0, ms=3)

#ax1.loglog(g15_nu, g15_nFn, 'co', zorder=-100, mew=0, ms=3)
#ax1.loglog(g15_xx, g15_yyy,  'b', zorder=-200)
#ax1.loglog(hes_nu, hes_nFn, 'mo', zorder=-100, mew=0, ms=3)
#ax1.loglog(hes_xx, hes_yyy,  'r', zorder=-200)
ax1.loglog(cas_nu, cas_nFn, 'yo', zorder=-100, mew=0, ms=3)
ax1.loglog(cas_xx, cas_yyy,  'g', zorder=-200)
#ax1.loglog(kes_nu, kes_nFn,  'k', zorder=-200)

ax1.loglog(gem_nu, gem_nFn, 'ro', zorder=-100, mew=0, ms=3)
ax1.loglog(vel_nu, vel_nFn, 'ko', zorder=-100, mew=0, ms=3)

#ax1.loglog(g31_nu, g31_nFn, 'co', zorder=-100, mew=0, ms=3)
#ax1.loglog(n15_nu, n15_nFn, 'mo', zorder=-100, mew=0, ms=3)
#ax1.loglog(g29_nu, g29_nFn, 'yo', zorder=-100, mew=0, ms=3)
#ax1.loglog(g76_nu, g76_nFn, 'ko', zorder=-100, mew=0, ms=3)

#ax1.loglog(xsp_xx, xsp_yyy,  'k', zorder=-200)

# Labels and limits
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("$\\nu F_\\nu$ (erg s$^{-1}$ cm$^{-2}$)")
ax1.set_xlim([  1e7,  1e29])
ax1.set_ylim([1e-22, 1e-10])

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
ax3.set_ylabel("$\\nu L_\\nu$ (erg s$^{-1}$)")


fig.savefig('/Users/silver/box/phd/pro/87a/lim/art/figs/mod_lim.pdf',bbox_inches='tight', pad_inches=0.03)

def plt_pos():
    fig = plt.figure(figsize=(12, 3))
    ax1 = fig.add_subplot(111)

    # Plot limits
    tmp_ind = (co==10) | (co==11) | (co==12) | (co==13) | (co==17)
    ax1.loglog(nu[tmp_ind], nu[tmp_ind]*fn[tmp_ind], 'v', color=kth_blue, mew=0, ms=8)
        
    ## Plot overplots
    ax1.loglog(crn_nu, crn_nFn, 'o', mec=kth_teal, mfc='w', zorder=-190, mew=3, ms=8)
    cas_pos = griddata(cas_xx, cas_yyy, cas_nu, method='linear')
    ax1.loglog(cas_nu[::10], cas_pos[::10], 'o', mec=kth_teal, mfc='w', zorder=-100, mew=3, ms=8)
#    ax1.loglog(cas_xx, cas_yyy,  'g', zorder=-200)
    
    # Labels and limits
#    ax1.set_xlabel("Frequency (Hz)")
#    ax1.set_ylabel("$\\nu F_\\nu$ (erg s$^{-1}$ cm$^{-2}$)")
    ax1.set_xlim([  1e7,  1e21])
    ax1.set_ylim([1e-16, 1e-10])
    
    # Upper x-axis
    up_x = np.logspace(-21, 21, 22)
    ax2 = ax1.twiny()
    ax2.set_xscale(ax1.get_xscale())
    ax2.set_xticks(cc/(up_x*1e-8))
    ax2.set_xticklabels(tic_hlp(up_x))
    ax2.set_xlim(ax1.get_xlim())
    ax2.minorticks_off()
#    ax2.set_xlabel("Wavelength (\AA{})")
    
    # Right y-axis
    ri_y = np.logspace(10, 40, 31)
    ax3 = ax1.twinx()
    ax3.set_yscale(ax1.get_yscale())
    ax3.set_yticks(ri_y/(4*np.pi*DD**2))
    ax3.set_yticklabels(tic_hlp(ri_y))
    ax3.set_ylim(ax1.get_ylim())
    ax3.minorticks_off()
#    ax3.set_ylabel("$\\nu L_\\nu$ (erg s$^{-1}$)")
    fig.savefig('/Users/silver/box/sci/doc/posters/2017/astronomdagarna/pos_lim.pdf',bbox_inches='tight', pad_inches=0.03)
    
plt_pos()

plt.show()
