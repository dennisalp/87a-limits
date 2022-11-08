# Dennis Alp 2016-10-20
# Align images using IRAF/PyRAF
# On my setup I need to work in the iraf27 environment: source activate iraf27
# Then astropy becomes accesible from both python and python2.7 (both are 2.7.12 anyway)
# time python /Users/$USER/Dropbox/bin/phd_87a_uplim_ali.py

import numpy as np
import os
from shutil import copyfile
import pdb

from glob import glob
from pyraf import iraf
from astropy.io import fits

#########
# Magnify
def magnify():
    iraf.magnify(INPUT, OUTPUT, MAG, MAG, interpolation = "linear", fluxconserve = "yes")

#########
# Pad
def pad():
    for ofile in out_files:
        ff = fits.open(ofile, mode='update')
        padded = np.zeros((NY,NX))
        xx = ff[0].data.shape[1]
        yy = ff[0].data.shape[0]
        xstart = (NX-xx)/2
        ystart = (NY-yy)/2
        padded[ystart:ystart+yy, xstart:xstart+xx] = ff[0].data
        ff[0].data = padded
        ff.close()

#########
# Move
def move(xstart, ystart):
    for ofile in out_files:
        ff = fits.open(ofile, mode='update')
        padded = np.zeros((NY,NX))
        xx = ff[0].data.shape[1]
        yy = ff[0].data.shape[0]
        padded[ystart:ystart+yy, xstart:xstart+xx] = ff[0].data
        ff[0].data = padded
        ff.close()

#########
# Align
def align():
    iraf.imalign(REFERENCE + "," + INPUT, REFERENCE, STARS, output=(REFERENCE[-29:] + "," + OUTPUT), shifts=SHIFTS, bigbox=45,  boxsize=25)

#########
# Geomap+tran
def geo():
    iraf.geomap(MAPS, "geomap.db", 1, 1150, 1, 1115, calctype = "real", interactive = "no")
#    out_files = [item[-10:-5]+"_map.fits" for item in files]
#    OUTPUT = ",".join(out_files)
    print os.listdir(".")
    iraf.geotran(INPUT, OUTPUT.replace("sft","map"), "geomap.db", MAPS)


#########
# Global parameters
NX = 4096
NY = 4096
    
##########
## VLT Parameters
#DAT_DIR = "/Users/silver/Dropbox/phd/data/vlt/87a/1red/"
#WRK_DIR = "/Users/silver/Dropbox/phd/data/vlt/87a/4ali/"
#WRK_DIR = "/Users/silver/Dropbox/phd/data/87a/"
#REFERENCE = "/Users/silver/Dropbox/phd/data/hst/87a/all/acs_r_2006-04-29_drz_al2.fits"
#STARS = "/Users/silver/Dropbox/phd/data/vlt/87a/4ali/close.dat"
#SHIFTS = "/Users/silver/Dropbox/phd/data/vlt/87a/4ali/shifts.dat"
#MAPS = "/Users/silver/Dropbox/phd/data/vlt/87a/4ali/2010h.map,/Users/silver/Dropbox/phd/data/vlt/87a/4ali/2012k.map"
#
#MAG = 13.221/25 # From NACO documentation, CONICA S13: http://www.eso.org/sci/facilities/paranal/instruments/naco/doc/VLT-MAN-ESO-14200-2761_v99.pdf
#MAG = 13.27/25 # From .fits-file, this is better compared to HST
##MAG = 8.9/25 # Testing
#
#os.chdir(WRK_DIR) #Move to designated directory
#map(os.remove, glob("*.fits"))
#map(os.remove, glob("*.db"))
#files = glob(DAT_DIR+"/*.fits") #Find all files.
##maps = glob("*.map") #Find all files.
#INPUT = ",".join(files)
#out_files = [item[-10:] for item in files]
#OUTPUT = ",".join(out_files)
#
#magnify()
#pad()
#INPUT = OUTPUT
#out_files = ["nac_" + item[-6] + "_" + item[-10:-6]+"-00-00_sft_al2.fits" for item in files]
#OUTPUT = ",".join(out_files)
#align()
#geo()

#########
# HST UVIS Parameters
DAT_DIR = "/Users/silver/Dropbox/phd/data/hst/87a/red/"
WRK_DIR = "/Users/silver/Dropbox/phd/data/87a/"
REFERENCE = "/Users/silver/Dropbox/phd/data/hst/87a/all/acs_r_2006-04-29_drz_al2.fits"
STARS =   "/Users/silver/Dropbox/phd/data/hst/87a/red/uv.dat"
SHIFTS =  "/Users/silver/Dropbox/phd/data/hst/87a/red/shifts.dat"

os.chdir(WRK_DIR) #Move to designated directory
map(os.remove, glob("*.fits")) # This removes files in a folder
files = glob(DAT_DIR+"/*.fits") #Find all files.
INPUT = ",".join(files)

#pdb.set_trace()
out_files = [item[-28:] for item in files]
OUTPUT = ",".join(out_files)

MAG = 0.502548188828556*0.0396/0.025 # Because I got the following files in a different scale
for i in range(0,len(files)):
    if any(substring in files[i] for substring in ["225","336","555","814"]):
        iraf.magnify(files[i], out_files[i], MAG, MAG, interpolation = "linear", fluxconserve = "yes")
    else:
        copyfile(files[i], out_files[i])

pad()
INPUT = OUTPUT#.replace(".fits",".fits[1]")
OUTPUT = OUTPUT.replace("una","al2")
align()

##########
## HST IR Parameters
#DAT_DIR = "/Users/silver/Dropbox/phd/data/hst/87a/drz/prd/"
#WRK_DIR = "/Users/silver/Dropbox/phd/data/87a/"
#REFERENCE = "/Users/silver/Dropbox/phd/data/hst/87a/all/acs_r_2006-04-29_drz_al2.fits"
#STARS =   "/Users/silver/Dropbox/phd/data/hst/87a/drz/prd/close.dat"
#SHIFTS =  "/Users/silver/Dropbox/phd/data/hst/87a/drz/prd/shifts.dat"
#
#os.chdir(WRK_DIR) #Move to designated directory
#map(os.remove, glob("*.fits"))
#files = glob(DAT_DIR+"/*.fits") #Find all files.
#INPUT = ",".join(files)
#
##pdb.set_trace()
#out_files = [item[-28:] for item in files]
#OUTPUT = ",".join(out_files)
#
#for i in range(0,len(files)):
#    copyfile(files[i], out_files[i])
#
#pad()
#INPUT = OUTPUT#.replace(".fits",".fits[1]")
#OUTPUT = OUTPUT.replace("una","al2")
#align()

##########
## ALMA Parameters
#NX = 1150
#NY = 1115
#DAT_DIR = "/Users/silver/Dropbox/phd/data/alm/87a/all/"
#WRK_DIR = "/Users/silver/Dropbox/phd/data/87a/"
#
#os.chdir(WRK_DIR) #Move to designated directory
#map(os.remove, glob("*.fits"))
#files = glob(DAT_DIR+"/*una.fits") #Find all files.
#INPUT = ",".join(files)
#
##pdb.set_trace()
#out_files = [item[-29:].replace("una","al2") for item in files]
#OUTPUT = ",".join(out_files)
#
#MAG = 0.006/0.025
#for i in range(0,len(files)):
#    iraf.magnify(files[i], out_files[i], MAG, MAG, interpolation = "linear", fluxconserve = "yes")
#
#xstart = 592.158345778-322.47022880232083*MAG
#ystart = 586.775275464-295.49159774774677*MAG
#move(xstart, ystart)
