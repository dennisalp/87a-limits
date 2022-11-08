# Dennis Alp 2016-12-07
# Find the limits in the ejecta of SN 1987A.
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/silver/box/bin/phd_87a_uplim_find_lim.py

import numpy as np
import os
import pdb
import time
from glob import glob
from datetime import date

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from scipy.ndimage.interpolation import shift
from pyraf import iraf
from astropy.io import fits

#########
# Parameters
# Logistics
WRK_DIR = "/Users/silver/lim"
os.chdir(WRK_DIR) #Move to designated directory
files = sorted(glob('*al2.fits')) #Find all files.

BIN_SEARCH_STOP = 0.01
RTRSH = 0.5 # Threshold radius for detection when inserting artificial sources
APERTURE = 8 # i.e. 0.2 arcsec
OUTPUT = "limits"
SHIFTS = np.array([0., 0.5]) # In image, [y, x]

# Source location 1-indexed
XCOO = 593.158345778
YCOO = 587.775275464
VKICK = 1.6e6
distance = 51.2 # kpc Panagia et al. (1991)
pc = 3.08567758149137e16

# Cut/subsets
# [554:619,549:634] python 0ind
# [573:599,580:606] python 0ind
BOX = [577,597,582,602]
NY = BOX[1]-BOX[0]
NX = BOX[3]-BOX[2]

# IRAF parameters
FAKE = "FAKE.fits" # Name of tampered image, handle with care!
COPY = "COPY.fits" # Make a new copy of original .fits but with different data type since daofind has issues
FSTARS = 'stars.dat'
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
"ac_b_2003-11-28_drz_al2.fits": 2.4 ,
"ac_b_2006-12-06_drz_al2.fits": 2.5 ,
"w3_b_2009-12-12_drz_al2.fits": 3.25,
"w3_b_2015-05-24_drz_al2.fits": 3.6 ,
"ac_b_0306-00-00_drz_al2.fits": 2.4 ,
"w2_b_0709-00-00_drz_al2.fits": 3.4 ,
"w2_b_1994-09-24_drz_al2.fits": 3.4 ,
"w3_b_0916-00-00_drz_al2.fits": 3.5 ,
"ac_r_2003-11-28_drz_al2.fits": 2.7 ,
"ac_r_2005-09-26_drz_al2.fits": 2.7 ,
"ac_r_2006-12-06_drz_al2.fits": 2.7 ,
"w3_r_2009-12-12_drz_al2.fits": 3.3 ,
"w3_r_2015-05-24_drz_al2.fits": 3.3 ,
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

# This is the correct gain for 1600 km s-1
sigma = {
"ac_b_0306-00-00_drz_al2.fits": 0.00908373987207 ,
"ac_b_2003-01-05_drz_al2.fits": 0.0228491797332 ,
"ac_b_2003-08-12_drz_al2.fits": 0.0291310652174 ,
"ac_b_2003-11-28_drz_al2.fits": 0.0213970884394 ,
"ac_b_2004-12-15_drz_al2.fits": 0.020658168207 ,
"ac_b_2005-04-02_drz_al2.fits": 0.0243872439181 ,
"ac_b_2006-12-06_drz_al2.fits": 0.0208096637587 ,
"ac_r_0306-00-00_drz_al2.fits": 0.00935002472611 ,
"ac_r_2003-01-05_drz_al2.fits": 0.040141273211 ,
"ac_r_2003-08-12_drz_al2.fits": 0.0535579717674 ,
"ac_r_2003-11-28_drz_al2.fits": 0.041459175453 ,
"ac_r_2005-09-26_drz_al2.fits": 0.0112160209713 ,
"ac_r_2006-04-15_drz_al2.fits": 0.0377824488393 ,
"ac_r_2006-04-29_drz_al2.fits": 0.0540112705669 ,
"ac_r_2006-12-06_drz_al2.fits": 0.0387848992594 ,
"na_h_2010-10-26_map_al2.fits": 2.25397237037 ,
"na_k_2012-12-14_map_al2.fits": 2.2206002385 ,
"w2_b_0709-00-00_drz_al2.fits": 0.0013497469855 ,
"w2_b_1994-09-24_drz_al2.fits": 0.0112337839964 ,
"w2_b_2007-05-12_drz_al2.fits": 0.00221964527333 ,
"w2_b_2008-02-19_drz_al2.fits": 0.00243933144874 ,
"w2_b_2009-04-29_drz_al2.fits": 0.00260947564194 ,
"w2_r_0709-00-00_drz_al2.fits": 0.00461112263583 ,
"w2_r_1994-09-24_drz_al2.fits": 0.0337231917625 ,
"w2_r_2007-05-12_drz_al2.fits": 0.00677326852443 ,
"w2_r_2008-02-19_drz_al2.fits": 0.00877493774831 ,
"w2_r_2009-04-29_drz_al2.fits": 0.00935226375538 ,
"w31102011-01-05_drz_al2.fits": 0.0317431056004 ,
"w31602011-01-05_drz_al2.fits": 0.0265429796289 ,
"w32252009-12-13_drz_al2.fits": 0.00831894548452 ,
"w322r2009-12-13_drz_al2.fits": 0.00919055180959 ,
"w322x2009-12-13_drz_al2.fits": 0.00912499971563 ,
"w33362009-12-13_drz_al2.fits": 0.0102150649143 ,
"w333r2009-12-13_drz_al2.fits": 0.0107315456091 ,
"w333x2009-12-13_drz_al2.fits": 0.0107135777348 ,
"w35022016-06-08_drz_al2.fits": 0.00426633609016 ,
"w35552009-12-13_drz_al2.fits": 0.0345544776124 ,
"w355r2009-12-13_drz_al2.fits": 0.0342234963542 ,
"w355x2009-12-13_drz_al2.fits": 0.0342613059043 ,
"w36452016-06-08_drz_al2.fits": 0.00467488675793 ,
"w36572016-06-08_drz_al2.fits": 0.0099913342328 ,
"w38142009-12-13_drz_al2.fits": 0.0285103289567 ,
"w381r2009-12-13_drz_al2.fits": 0.0293689603226 ,
"w381x2009-12-13_drz_al2.fits": 0.02937010042 ,
"w3_b_0916-00-00_drz_al2.fits": 0.00511486636045 ,
"w3_b_2009-12-12_drz_al2.fits": 0.0149656854107 ,
"w3_b_2011-01-05_drz_al2.fits": 0.012493933925 ,
"w3_b_2013-02-06_drz_al2.fits": 0.0115889594186 ,
"w3_b_2014-06-15_drz_al2.fits": 0.0110712054806 ,
"w3_b_2015-05-24_drz_al2.fits": 0.0111788292546 ,
"w3_b_2016-06-08_drz_al2.fits": 0.0176627597723 ,
"w3_br2015-05-24_drz_al2.fits": 0.0114022165672 ,
"w3_bx2015-05-24_drz_al2.fits": 0.0114063219108 ,
"w3_r_0916-00-00_drz_al2.fits": 0.00835251815982 ,
"w3_r_2009-12-12_drz_al2.fits": 0.0145355466277 ,
"w3_r_2011-01-05_drz_al2.fits": 0.0246301014677 ,
"w3_r_2013-02-06_drz_al2.fits": 0.0223027000739 ,
"w3_r_2014-06-15_drz_al2.fits": 0.0216617473434 ,
"w3_r_2015-05-24_drz_al2.fits": 0.0213578653205 ,
"w3_r_2016-06-08_drz_al2.fits": 0.0321654028623 ,
"w3_rr2015-05-24_drz_al2.fits": 0.0216938120373 ,
"w3_rx2015-05-24_drz_al2.fits": 0.0217165616752 }

# This is the correct gain for 400 km s-1
#sigma = {
#"ac_b_0306-00-00_drz_al2.fits": 0.00804003697591 ,
#"ac_b_2003-01-05_drz_al2.fits": 0.0203991667965 ,
#"ac_b_2003-08-12_drz_al2.fits": 0.0254123941445 ,
#"ac_b_2003-11-28_drz_al2.fits": 0.0192891841515 ,
#"ac_b_2004-12-15_drz_al2.fits": 0.0187255617157 ,
#"ac_b_2005-04-02_drz_al2.fits": 0.0209396982412 ,
#"ac_b_2006-12-06_drz_al2.fits": 0.0185964299217 ,
#"ac_r_0306-00-00_drz_al2.fits": 0.00844848526691 ,
#"ac_r_2003-01-05_drz_al2.fits": 0.0381147976978 ,
#"ac_r_2003-08-12_drz_al2.fits": 0.050070286133 ,
#"ac_r_2003-11-28_drz_al2.fits": 0.037321446243 ,
#"ac_r_2005-09-26_drz_al2.fits": 0.0101263626838 ,
#"ac_r_2006-04-15_drz_al2.fits": 0.0338693534419 ,
#"ac_r_2006-04-29_drz_al2.fits": 0.0505959122732 ,
#"ac_r_2006-12-06_drz_al2.fits": 0.0350272800492 ,
#"na_h_2010-10-26_map_al2.fits": 2.24928576173 ,
#"na_k_2012-12-14_map_al2.fits": 2.21908655059 ,
#"w2_b_0709-00-00_drz_al2.fits": 0.00120069159413 ,
#"w2_b_1994-09-24_drz_al2.fits": 0.0107149436258 ,
#"w2_b_2007-05-12_drz_al2.fits": 0.00204408091332 ,
#"w2_b_2008-02-19_drz_al2.fits": 0.00218124791075 ,
#"w2_b_2009-04-29_drz_al2.fits": 0.00230463141915 ,
#"w2_r_0709-00-00_drz_al2.fits": 0.00413801146008 ,
#"w2_r_1994-09-24_drz_al2.fits": 0.0329843435858 ,
#"w2_r_2007-05-12_drz_al2.fits": 0.00600026201078 ,
#"w2_r_2008-02-19_drz_al2.fits": 0.00808990877227 ,
#"w2_r_2009-04-29_drz_al2.fits": 0.0080919337969 ,
#"w31102011-01-05_drz_al2.fits": 0.0302127706271 ,
#"w31602011-01-05_drz_al2.fits": 0.0264220377503 ,
#"w32252009-12-13_drz_al2.fits": 0.00787481502908 ,
#"w322r2009-12-13_drz_al2.fits": 0.00880833767334 ,
#"w322x2009-12-13_drz_al2.fits": 0.00897974499237 ,
#"w33362009-12-13_drz_al2.fits": 0.00915132626068 ,
#"w333r2009-12-13_drz_al2.fits": 0.00970162126774 ,
#"w333x2009-12-13_drz_al2.fits": 0.0097801420786 ,
#"w35022016-06-08_drz_al2.fits": 0.00390669163071 ,
#"w35552009-12-13_drz_al2.fits": 0.0301530573838 ,
#"w355r2009-12-13_drz_al2.fits": 0.0296174753484 ,
#"w355x2009-12-13_drz_al2.fits": 0.0296530847633 ,
#"w36452016-06-08_drz_al2.fits": 0.00371134971152 ,
#"w36572016-06-08_drz_al2.fits": 0.00770612441247 ,
#"w38142009-12-13_drz_al2.fits": 0.0254550332237 ,
#"w381r2009-12-13_drz_al2.fits": 0.0261615164765 ,
#"w381x2009-12-13_drz_al2.fits": 0.0262393487314 ,
#"w3_b_0916-00-00_drz_al2.fits": 0.00437791470828 ,
#"w3_b_2009-12-12_drz_al2.fits": 0.0129506066244 ,
#"w3_b_2011-01-05_drz_al2.fits": 0.0108739941548 ,
#"w3_b_2013-02-06_drz_al2.fits": 0.00988693968902 ,
#"w3_b_2014-06-15_drz_al2.fits": 0.00900820321504 ,
#"w3_b_2015-05-24_drz_al2.fits": 0.00951089921685 ,
#"w3_b_2016-06-08_drz_al2.fits": 0.0156495176159 ,
#"w3_br2015-05-24_drz_al2.fits": 0.00980976550171 ,
#"w3_bx2015-05-24_drz_al2.fits": 0.00985531843612 ,
#"w3_r_0916-00-00_drz_al2.fits": 0.00728052046935 ,
#"w3_r_2009-12-12_drz_al2.fits": 0.0128006002307 ,
#"w3_r_2011-01-05_drz_al2.fits": 0.020720224599 ,
#"w3_r_2013-02-06_drz_al2.fits": 0.0185325257105 ,
#"w3_r_2014-06-15_drz_al2.fits": 0.0172285327037 ,
#"w3_r_2015-05-24_drz_al2.fits": 0.0176451364889 ,
#"w3_r_2016-06-08_drz_al2.fits": 0.0263625246168 ,
#"w3_rr2015-05-24_drz_al2.fits": 0.0179960318673 ,
#"w3_rx2015-05-24_drz_al2.fits": 0.0181149162572 }

# This is the correct gain for 800 km s-1
#sigma = {
#"ac_b_0306-00-00_drz_al2.fits": 0.00835699112049 ,
#"ac_b_2003-01-05_drz_al2.fits": 0.0207189150881 ,
#"ac_b_2003-08-12_drz_al2.fits": 0.0266817330091 ,
#"ac_b_2003-11-28_drz_al2.fits": 0.0197441767194 ,
#"ac_b_2004-12-15_drz_al2.fits": 0.0193651305348 ,
#"ac_b_2005-04-02_drz_al2.fits": 0.0221000007033 ,
#"ac_b_2006-12-06_drz_al2.fits": 0.0191832471222 ,
#"ac_r_0306-00-00_drz_al2.fits": 0.00887877443251 ,
#"ac_r_2003-01-05_drz_al2.fits": 0.0389824735051 ,
#"ac_r_2003-08-12_drz_al2.fits": 0.0511026921974 ,
#"ac_r_2003-11-28_drz_al2.fits": 0.0387432804522 ,
#"ac_r_2005-09-26_drz_al2.fits": 0.0107191210819 ,
#"ac_r_2006-04-15_drz_al2.fits": 0.0353253492111 ,
#"ac_r_2006-04-29_drz_al2.fits": 0.0517702201209 ,
#"ac_r_2006-12-06_drz_al2.fits": 0.0366317055757 ,
#"na_h_2010-10-26_map_al2.fits": 2.25223377716 ,
#"na_k_2012-12-14_map_al2.fits": 2.21986271601 ,
#"w2_b_0709-00-00_drz_al2.fits": 0.00124416803283 ,
#"w2_b_1994-09-24_drz_al2.fits": 0.0112337839964 ,
#"w2_b_2007-05-12_drz_al2.fits": 0.00204408091332 ,
#"w2_b_2008-02-19_drz_al2.fits": 0.00220996927684 ,
#"w2_b_2009-04-29_drz_al2.fits": 0.00250512713014 ,
#"w2_r_0709-00-00_drz_al2.fits": 0.00445188311827 ,
#"w2_r_1994-09-24_drz_al2.fits": 0.0329843435858 ,
#"w2_r_2007-05-12_drz_al2.fits": 0.00640765899687 ,
#"w2_r_2008-02-19_drz_al2.fits": 0.00830794369396 ,
#"w2_r_2009-04-29_drz_al2.fits": 0.00887472555493 ,
#"w31102011-01-05_drz_al2.fits": 0.031333406738 ,
#"w31602011-01-05_drz_al2.fits": 0.0264767650475 ,
#"w32252009-12-13_drz_al2.fits": 0.00805602755352 ,
#"w322r2009-12-13_drz_al2.fits": 0.0088973874362 ,
#"w322x2009-12-13_drz_al2.fits": 0.00897974499237 ,
#"w33362009-12-13_drz_al2.fits": 0.009769511511 ,
#"w333r2009-12-13_drz_al2.fits": 0.0103171249684 ,
#"w333x2009-12-13_drz_al2.fits": 0.0103761467072 ,
#"w35022016-06-08_drz_al2.fits": 0.0040084058677 ,
#"w35552009-12-13_drz_al2.fits": 0.0318457844766 ,
#"w355r2009-12-13_drz_al2.fits": 0.0316919984452 ,
#"w355x2009-12-13_drz_al2.fits": 0.0319366019181 ,
#"w36452016-06-08_drz_al2.fits": 0.00397529530419 ,
#"w36572016-06-08_drz_al2.fits": 0.00886888732454 ,
#"w38142009-12-13_drz_al2.fits": 0.0276292715035 ,
#"w381r2009-12-13_drz_al2.fits": 0.02880011569 ,
#"w381x2009-12-13_drz_al2.fits": 0.029021599921 ,
#"w3_b_0916-00-00_drz_al2.fits": 0.00471698323766 ,
#"w3_b_2009-12-12_drz_al2.fits": 0.0140333690138 ,
#"w3_b_2011-01-05_drz_al2.fits": 0.0114474455416 ,
#"w3_b_2013-02-06_drz_al2.fits": 0.0103525449996 ,
#"w3_b_2014-06-15_drz_al2.fits": 0.00972024891876 ,
#"w3_b_2015-05-24_drz_al2.fits": 0.0104725814586 ,
#"w3_b_2016-06-08_drz_al2.fits": 0.0160167110346 ,
#"w3_br2015-05-24_drz_al2.fits": 0.010692775895 ,
#"w3_bx2015-05-24_drz_al2.fits": 0.0106979365494 ,
#"w3_r_0916-00-00_drz_al2.fits": 0.00792838848224 ,
#"w3_r_2009-12-12_drz_al2.fits": 0.0139054396777 ,
#"w3_r_2011-01-05_drz_al2.fits": 0.0230962473663 ,
#"w3_r_2013-02-06_drz_al2.fits": 0.0205453258392 ,
#"w3_r_2014-06-15_drz_al2.fits": 0.0192353367381 ,
#"w3_r_2015-05-24_drz_al2.fits": 0.0192548271699 ,
#"w3_r_2016-06-08_drz_al2.fits": 0.0286244415553 ,
#"w3_rr2015-05-24_drz_al2.fits": 0.019573106933 ,
#"w3_rx2015-05-24_drz_al2.fits": 0.0196537293536 }
    
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

## The range of wavelengths spanned by the filters, [central, width], in Angstroms
#filter_width = {
#"ac_b_2003-11-28_drz_al2.fits": [ 4297, 1038], # F435W Johnson B http://www.stsci.edu/hst/acs/documents/handbooks/current/c05_imaging2.html
#"ac_r_2003-11-28_drz_al2.fits": [ 6318, 1442], # F625W SDSS r
#"na_h_2010-10-26_map_al2.fits": [16600, 3300], # NaCo H http://www.eso.org/sci/facilities/paranal/instruments/naco/inst/filters.html
#"na_k_2012-12-14_map_al2.fits": [21800, 3500], # NaCo Ks Equivalent https://www.eso.org/sci/facilities/paranal/instruments/naco/doc/VLT-MAN-ESO-14200-2761_v82.pdf
#"w2_b_1994-09-24_drz_al2.fits": [ 4311,  473], # F439W http://www.stsci.edu/hst/wfpc2/documents/wfpc2_filters_archive.html
#"w2_r_1994-09-24_drz_al2.fits": [ 6717,  867], # F675W
#"w31102011-01-05_drz_al2.fits": [11534, 4430], # F110W http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c07_ir06.html
#"w31602011-01-05_drz_al2.fits": [15369, 2683], # F160W
#"w32252009-12-13_drz_al2.fits": [ 2359,  467], # F225W http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c06_uvis06.html
#"w33362009-12-13_drz_al2.fits": [ 3355,  511], # F336W
#"w35022016-06-08_drz_al2.fits": [ 5010,   65], # F502N
#"w35552009-12-13_drz_al2.fits": [ 5308, 1562], # F555W
#"w36452016-06-08_drz_al2.fits": [ 6454,   84], # F645N
#"w38142009-12-13_drz_al2.fits": [ 8024, 1536], # F814W
#"w3_b_2015-05-24_drz_al2.fits": [ 4325,  618], # F438W
#"w3_b_2016-06-08_drz_al2.fits": [ 4325,  618], # F438W
#"w3_r_2015-05-24_drz_al2.fits": [ 6242, 1463], # F625W
#"w3_r_2016-06-08_drz_al2.fits": [ 6242, 1463], # F625W
#}
#
## Convert the instrumental electron per second to physical units
#unit_conv = {
#"ac_b_2003-11-28_drz_al2.fits": 5.332e-19, # .fits probably populated by old values
#"ac_r_2003-11-28_drz_al2.fits": 1.948e-19, # from the online app, from Bohlin (2016), synphot has values without interpolation
#"na_h_2010-10-26_map_al2.fits": 1.e-20, # dummies
#"na_k_2012-12-14_map_al2.fits": 1.e-20, #
#"w2_b_1994-09-24_drz_al2.fits": 2.945e-17, # matches .fits
#"w2_r_1994-09-24_drz_al2.fits": 2.899e-18, # matches .fits
#"w31102011-01-05_drz_al2.fits": 1.5274e-20, # matches .fits
#"w31602011-01-05_drz_al2.fits": 1.9276e-20, # matches .fits
#"w32252009-12-13_drz_al2.fits": 4.3542E-18, # UVIS1 matches nothing
#"w33362009-12-13_drz_al2.fits": 1.4218E-18, # UVIS1 matches nothing
#"w35022016-06-08_drz_al2.fits": 5.5958E-18, # UVIS2 matches PHTFLAM2
#"w35552009-12-13_drz_al2.fits": 1.9812E-19, # UVIS1 matches nothing
#"w36452016-06-08_drz_al2.fits": 3.6232E-18, # UVIS2 matches PHTFLAM2
#"w38142009-12-13_drz_al2.fits": 1.6305E-19, # UVIS1 matches nothing
#"w3_b_2015-05-24_drz_al2.fits": 7.3404E-19, # UVIS2 matches PHTFLAM2
#"w3_b_2016-06-08_drz_al2.fits": 7.3404E-19, # UVIS2 matches PHTFLAM2
#"w3_r_2015-05-24_drz_al2.fits": 1.8740E-19, # UVIS2 matches PHTFLAM2
#"w3_r_2016-06-08_drz_al2.fits": 1.8740E-19, # UVIS2 matches PHTFLAM2
#}
#
## Encircled energy since photflam is given for apertures of 10 pixels (uvis pixels?, unclear)
#aperture = {
#"ac_b_2003-11-28_drz_al2.fits": 1., # https://acszeropoints.stsci.edu/
#"ac_r_2003-11-28_drz_al2.fits": 1., # http://www.stsci.edu/hst/acs/analysis/zeropoints
#"na_h_2010-10-26_map_al2.fits": 1., #
#"na_k_2012-12-14_map_al2.fits": 1., #
#"w2_b_1994-09-24_drz_al2.fits": 1., # http://www.stsci.edu/hst/wfpc2/analysis/wfpc2_photflam.html
#"w2_r_1994-09-24_drz_al2.fits": 1., #
#"w31102011-01-05_drz_al2.fits": 1., # http://www.stsci.edu/hst/wfc3/ir_phot_zpt
#"w31602011-01-05_drz_al2.fits": 1., #
#"w32252009-12-13_drz_al2.fits": 1., # UVIS1 http://www.stsci.edu/hst/wfc3/analysis/uvis_zpts
#"w33362009-12-13_drz_al2.fits": 1., # UVIS1
#"w35022016-06-08_drz_al2.fits": 1., # UVIS2
#"w35552009-12-13_drz_al2.fits": 1., # UVIS1
#"w36452016-06-08_drz_al2.fits": 1., # UVIS2
#"w38142009-12-13_drz_al2.fits": 1., # UVIS1
#"w3_b_2015-05-24_drz_al2.fits": 1., # UVIS2
#"w3_b_2016-06-08_drz_al2.fits": 1., # UVIS2
#"w3_r_2015-05-24_drz_al2.fits": 1., # UVIS2
#"w3_r_2016-06-08_drz_al2.fits": 1., # UVIS2
#}
    
#########
# Help functions
# Just reasonable plot of sky images, good for debugging
def sky_plt(image):
    plt.imshow(image,interpolation='nearest', cmap='afmhot', origin='lower')
    plt.show()

# Utility for making fits image out of a numpy array
def mk_fits(image, output):
    hdu = fits.PrimaryHDU(image)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(output, clobber=True)
    hdulist.close()
    
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

# adding but removing counts so that to make it just fit, i.e. max(ejecta,artificial)
def hide_star(imgcp, psfcp, flux):
    norm = np.sum(psfcp)
    psfcp = flux/norm * psfcp
    test_psf = np.zeros(imgcp.shape)
    xstart = int(NX/2 - APERTURE)
    ystart = int(NY/2 - APERTURE)
    
    test_psf[ystart:ystart+2*APERTURE+1,xstart:xstart+2*APERTURE+1] = psfcp
    if np.all(test_psf < imgcp):
        mk_fits(imgcp, FAKE)
        return False
    
    test_img = np.where(test_psf > imgcp, test_psf, imgcp)
    mk_fits(test_img, FAKE)
    return True

# tries to find hidden star using daofind
def try_find():
    # Find psf
    """
    BUG: In rare circumstances daofind may abort with a "pixel file truncation
        error" when trying to read back the convolved images it has just 
    	written. This only occurs on certain sized images and is due to
    	the interaction of the read, write and boundary extension in image
    	i/o. For example daofind works fine on a 640 by 1024 image but fails on
    	one that is 641 by 1025 pixels.
    
    STATUS:	The problem is fixed in 2.9. The solution was to add an imflush
    	statement to flush the image i/o buffers after the write was
    	complete and before initiating the read. There is no workaround.
    	Contact the IRAF group for a fix.
    
    https://github.com/joequant/iraf/blob/master/local/bugs.log
    """
    iraf.daofind(FAKE,output=FSTARS,fwhmpsf=fwhmfit,sigma=sigmaobs,threshold=3,datamin=-0.1,datamax=9999,verify=False, Stdout=1)

    
    stscoo=iraf.pdump(FSTARS,'XCENTER,YCENTER','yes',Stdout=1)
    
    # Check if found
    for xy in stscoo:
        testx, testy = xy.split()
        testx = float(testx)
        testy = float(testy)
        delx = (testx-1)-NX/2-SHIFTS[1]
        dely = (testy-1)-NY/2-SHIFTS[0]
        if np.sqrt(delx**2+dely**2) < fwhmfit/2.:
            return True

    return False
    
# this uses binary search to find a limit by hiding and daofinding
def bin_search(img_path, psf, yy, xx):
    xstart = int(xx-NX/2)
    xstop  = int(xx+NX/2+1)
    ystart = int(yy-NY/2)
    ystop  = int(yy+NY/2+1)
    img = fits.open(img_path)[0].data[ystart:ystop,xstart:xstop]
    
    xstart = int(NX/2 - APERTURE)
    ystart = int(NY/2 - APERTURE)
    upper = np.sum(img[ystart:ystart+2*APERTURE+1,xstart:xstart+2*APERTURE+1])/2 # 2 is heuristic
    
    up = 2*upper
    lo = 1e-12 # Suppress divide by 0 warning
    oflux = -1
    flux = 1e-12
    
    # Catch detections of stuff in ejecta
    hide_star(img.copy(), psf.copy(), flux)
    if try_find():
        while abs(flux/oflux-1) > BIN_SEARCH_STOP:
            print "WARNING: DETECTION"
            if hide_star(img.copy(), psf.copy(), flux):
                up = flux
            else:
                lo = flux
            oflux = flux
            flux = (up + lo)/2.
    else:
        while abs(flux/oflux-1) > BIN_SEARCH_STOP:
#            print flux, oflux, flux/oflux
            hide_star(img.copy(), psf.copy(), flux)
            
            if try_find():
                up = flux
            else:
                lo = flux
            oflux = flux
            flux = (up + lo)/2.
                
    return flux


#########
# Definitions
limits = np.zeros((NY,NX,len(files)))
img_nr = 0
iraf.stsdas()
iraf.hst_calib()

#########
# Insert psf and try to find it
for img_path in files[0:]:
    key = img_path[-28:]
    print key

    # To cast data type, this is SUPER IMPORTANT. Daofind throws fits otherwise
    mk_fits(fits.open(img_path)[0].data.astype('>f8'), COPY)
    
    # Prepare the psf
    # Find the FWHM by fitting to the see result from IRAF
    psf = fits.open(img_path+".see.2.fits")[0].data
    subny = psf.shape[0]
    subnx = psf.shape[1]    
    rad_prof = radial_profile(psf)
    rmax = np.int(np.ceil(fwhmpsf[img_path])+3)
    fwhmfit = np.sqrt(2)*griddata(rad_prof[0:rmax],np.arange(0,rmax),np.amax(rad_prof)/2,method='linear')
    print "Check match iraf, fit:",fwhmpsf[key],fwhmfit
    sigmaobs = sigma[key]

    subny = int(subny/2)
    subnx = int(subnx/2)
    psf = shift(psf, SHIFTS, order=1, mode='constant', cval=0.0, prefilter=False)
    psf = psf[subny-APERTURE:subny+APERTURE+1,subnx-APERTURE:subnx+APERTURE+1]
    
    for ii in range(BOX[0],BOX[1]):
        for jj in range(BOX[2],BOX[3]):
            # Do the binary search for a limit at a point
            unfluxed = bin_search(COPY, psf, np.floor(ii), np.floor(jj))
            limits[ii-BOX[0],jj-BOX[2],img_nr] = unfluxed
            print "Limit:", key, jj, ii, limits[ii-BOX[0],jj-BOX[2],img_nr]

    img_nr += 1
    print ""

np.save(OUTPUT + "_y" + str(SHIFTS[0]) + "_x" + str(SHIFTS[1]) + "_v" + str(VKICK/1e3) + ".npy", limits)
