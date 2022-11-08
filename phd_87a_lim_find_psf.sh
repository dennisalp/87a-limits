#!/bin/bash
# Dennis Alp 2016-11-17
# Find psf of a star field

# Display the image with display in ds9, then daoedit
# to get interactive FWHM-plot using the i key.

function find_psf {
    printf "!echo ECHOES ARE COMMENTS
iraf.Verbose=2
daophot

!echo INITIAL IDENTIFICATION OF STARS .coo
daofind $img $img.coo.1
$fwhmpsf
$skystd
$threshold
$mindat
$maxdat

!echo INITIAL PHOTOMETRY .mag
phot $img $img.coo.1 $img.mag.1


$ffwhmpsf
$ffwhmpsf
$fwhmpsf
$skystd
$mindat
$maxdat

!echo SELECT STARS FOR PSF BUILD .pst
pstselect $img $img.mag.1 $img.pst.1
$maxnpsf
$psfrad
$fitrad
$mindat
$maxdat

!echo REMOVE SOME UNWANTED STARS, THESE HAVE TO BE MANUALLY SELECTED
!echo THESE ALSO DEPEND ON DAOFIND PARAMS OF COURSE
" > pyraf.cmd
    
    for i in ${reject[@]}
    do
	printf "!sed '/^${i} /d' $img.pst.1 > $img.tmp.1; mv $img.tmp.1 $img.pst.1\n" >> pyraf.cmd
    done

    printf "fields $img.pst.1 2,3 > $img.pts.1

!echo COMPUTE THE PSF .psf .pst.2 .psg
!echo SECOND .pst IS STARS THAT ACTUALLY WERE USED
psf $img $img.mag.1 $img.pst.1 $img.psf.1 $img.pst.2 $img.psg.1 $inter
$analytical
$varorder
$psfrad
$fitrad
$mindat
$maxdat
fields $img.pst.2 2,3 > $img.pts.2

!echo GROUP THE STARS FOR NSTAR .grp
group $img $img.mag.1 $img.psf.1 $img.grp.1 
$psfrad
$fitrad
$critovl
$mindat
$maxdat

!echo NSTAR FITTING .nst .nrj
nstar $img $img.psg.1 $img.psf.1 $img.nst.1 $img.nrj.1
yes
no
yes
$psfrad
$fitrad
$maxgroup
$mindat
$maxdat

!echo SUBSTAR FOR SUBTRACTION VERIFICATION .sub
substar $img $img.nst.1 \"\" $img.psf.1 $img.sub.1
$psfrad
$mindat
$maxdat

!echo SUBSTAR FOR ITERATION .sub.2
substar $img $img.nst.1 $img.pst.2 $img.psf.1 $img.sub.2
$psfrad
$mindat
$maxdat

!echo COMPUTE FIRST PSF ITERATION ON SUBTRACTED IMAGE
psf $img.sub.2 $img.psg.1 $img.pst.2 $img.psf.2 $img.pst.3 $img.psg.2 inter-
$analytical
$varorder
$psfrad
$fitrad
$mindat
$maxdat
fields $img.pst.3 2,3 > $img.pts.3

!echo FIRST ITERATION OF ALLSTAR .als .sub.3 .alj
allstar $img $img.mag.1 $img.psf.2 $img.als.1 $img.arj.1 $img.sub.3
yes
yes
no
$psfrad
$fitrad
$maxgroup
$mindat
$maxdat

!echo JUST TO COMPUTE THE CHI2 FOR PSF STARS ONLY, NOT THE FASTEST METHOD .nst .nrj
nstar $img $img.psg.2 $img.psf.2 $img.nst.2 $img.nrj.2
yes
no
yes
$psfrad
$fitrad
$maxgroup
$mindat
$maxdat

!echo FIRST ITERATION TO PICK UP MISSED STARS, CURRENTLY NOT FOLLOWING UP
daofind $img.sub.3 $img.coo.2
$fwhmpsf
$skystd
$threshold
$mindat
$maxdat

seepsf $img.psf.1 $img.see.1 xpsf=593.158345778 ypsf=587.775275464
seepsf $img.psf.2 $img.see.2 xpsf=593.158345778 ypsf=587.775275464
display $img 1 z1=0 z2=$maxdat ztrans=log
display $img.sub.1.fits 2 z1=0 z2=$maxdat ztrans=log 
display $img.see.1.fits 3 z1=0 z2=$maxdat ztrans=log 
display $img.see.2.fits 4 z1=0 z2=$maxdat ztrans=log 
pdump $img.pst.2 xcenter,ycenter yes | tvmark 1 STDIN col=205 point=1
pdump $img.nst.2 chi yes > chi2.dat
pdump $img.nst.2 chi yes | graph STDIN point+

.exit" >> pyraf.cmd

    pyraf < pyraf.cmd

    printf "import numpy as np
dat = np.loadtxt(\"chi2.dat\")
print \"Average chi2 value:\", np.mean(dat)
" > python.cmd
    python python.cmd
}

#########
# HST 110 2011
function hst110a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/110a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=6.5
    ffwhmpsf=26
    psfrad=27
    fitrad=6.5
    threshold=3
    skystd=0.016
    mindat=-0.1
    maxdat=200
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(237 478 15 107 57 120 416 325 5 487 418)
    find_psf
}

#########
# HST 160 2011
function hst160a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/160a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=7.2
    ffwhmpsf=28.8
    psfrad=29.8
    fitrad=7.2
    threshold=3
    skystd=0.024
    mindat=-0.1
    maxdat=200
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(275 130 9 54 37 62 239 2 184 283 238)
    find_psf
}

#########
# VLT H 2010
function vlt163a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/163a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=4
    ffwhmpsf=16
    psfrad=17
    fitrad=4
    threshold=3
    skystd=2.2
    mindat=-0.1
    maxdat=8000
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(273 214 228 351 250 243 161 175)
    find_psf
}

#########
# VLT K 2012
function vlt219a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/219a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.4
    ffwhmpsf=13.6
    psfrad=14.6
    fitrad=3.4
    threshold=3
    skystd=2.2
    mindat=-0.1
    maxdat=8000
    maxnpsf=16
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(128 139 126)
    find_psf
}

#########
# HST 225 2009
function hst225a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/225a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.7
    ffwhmpsf=14.8
    psfrad=15.8
    fitrad=3.7
    threshold=3
    skystd=0.0047
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(62 18 83 78 6 180 16)
    find_psf
}

#########
# HST 336 2009
function hst336a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/336a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.7
    ffwhmpsf=14.8
    psfrad=15.8
    fitrad=3.7
    threshold=3
    skystd=0.0045
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(172 84 487 219 211 161 31)
    find_psf
}

#########
# HST 502 2009
function hst502a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/502a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.6
    ffwhmpsf=14.4
    psfrad=15.4
    fitrad=3.6
    threshold=3
    skystd=0.0033
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(213 122 81 138 43 6 21 288 298)
    find_psf
}

#########
# HST 555 2009
function hst555a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/555a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.5
    ffwhmpsf=14
    psfrad=15
    fitrad=3.5
    threshold=3
    skystd=0.012
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(231 130 280 17 480 267 223 90 304 422 319)
    find_psf
}

#########
# HST 645 2016
function hst645a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/645a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.3
    ffwhmpsf=13.2
    psfrad=14.2
    fitrad=3.3
    threshold=3
    skystd=0.0022
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(257 160 216 109 185 179 82 246 310)
    find_psf
}

########
# HST 814 2009
function hst814a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/814a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.7
    ffwhmpsf=14.8
    psfrad=15.8
    fitrad=3.7
    threshold=3
    skystd=0.01
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(629 482 690 395 191 424 382 329 569)
    find_psf
}

#########
# HST Blue 2003-01-05
function hstb03a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b03a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.4
    ffwhmpsf=9.6
    psfrad=10.6
    fitrad=2.4
    threshold=3
    skystd=0.01
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(999999999)
    find_psf
}

#########
# HST Blue 2003-08-12
function hstb03b {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b03b
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.4
    ffwhmpsf=9.6
    psfrad=10.6
    fitrad=2.4
    threshold=3
    skystd=0.01375
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(999999999999)
    find_psf
}

#########
# HST Blue 2003-11-28
function hstb03c {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b03c
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.4
    ffwhmpsf=9.6
    psfrad=10.6
    fitrad=2.4
    threshold=3
    skystd=0.011
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(364 463 505)
    find_psf
}

#########
# HST Blue 2004-12-15
function hstb04a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b04a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.4
    ffwhmpsf=9.6
    psfrad=10.6
    fitrad=2.4
    threshold=3
    skystd=0.01
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(999999999999)
    find_psf
}

#########
# HST Blue 2005-04-02
function hstb05a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b05a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.4
    ffwhmpsf=9.6
    psfrad=10.6
    fitrad=2.4
    threshold=3
    skystd=0.01
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(9999999)
    find_psf
}

#########
# HST Blue 2006
function hstb06a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b06a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.5
    ffwhmpsf=10
    psfrad=11
    fitrad=2.5
    threshold=3
    skystd=0.01
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(399 527 569 560 516 839)
    find_psf
}

#########
# HST Blue 2007
function hstb07a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b07a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.5
    ffwhmpsf=14
    psfrad=15
    fitrad=3.5
    threshold=3
    skystd=0.00075
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(999999999)
    find_psf
}

#########
# HST Blue 2008
function hstb08a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b08a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.5
    ffwhmpsf=14
    psfrad=15
    fitrad=3.5
    threshold=3
    skystd=0.00067
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(999999999)
    find_psf
}

#########
# HST Blue 2009-04-29
function hstb09a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b09a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.5
    ffwhmpsf=14
    psfrad=15
    fitrad=3.5
    threshold=3
    skystd=0.00058
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(999999999)
    find_psf
}

#########
# HST Blue 2009-12-12
function hstb09b {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b09b
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.25
    ffwhmpsf=13
    psfrad=14
    fitrad=3.25
    threshold=3
    skystd=0.005
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(480 148 350 336 611 51)
    find_psf
}

#########
# HST Blue 2011
function hstb11a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b11a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.6
    ffwhmpsf=14.4
    psfrad=15.4
    fitrad=3.6
    threshold=3
    skystd=0.0065
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(893 484 347 884 53 526 599)
    reject=(99999999)
    find_psf
}

#########
# HST Blue 2013
function hstb13a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b13a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.5
    ffwhmpsf=14
    psfrad=15
    fitrad=3.5
    threshold=3
    skystd=0.00375
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(893 484 347 884 53 526 599)
    reject=(99999999)
    find_psf
}

#########
# HST Blue 2014
function hstb14a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b14a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.5
    ffwhmpsf=14
    psfrad=15
    fitrad=3.5
    threshold=3
    skystd=0.003375
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(893 484 347 884 53 526 599)
    reject=(99999999)
    find_psf
}

#########
# HST Blue 2015
function hstb15a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b15a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.6
    ffwhmpsf=14.4
    psfrad=15.4
    fitrad=3.6
    threshold=3
    skystd=0.0045
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(904 494 347 54)
    find_psf
}

#########
# HST Blue 2016
function hstb16a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b16a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.4
    ffwhmpsf=13.6
    psfrad=14.6
    fitrad=3.4
    threshold=3
    skystd=0.00975
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(108 202 376 232 225 265 351 10 119 40 160 195 99 136 247 17)
    reject=(99999999)
    find_psf
}

#########
# HST Blue 0306
function hstb36a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b36a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.4
    ffwhmpsf=9.6
    psfrad=10.6
    fitrad=2.4
    threshold=3
    skystd=0.004
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(974 506 660 727 714 20 120 79)
    find_psf
}

#########
# HST Blue 0709
function hstb79a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b79a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.4
    ffwhmpsf=13.6
    psfrad=14.6
    fitrad=3.4
    threshold=3
    skystd=0.000233
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(404 215 150 243 235 209 351 279 195)
    find_psf
}

#########
# HST Blue 1994
function hstb94a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b94a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.4
    ffwhmpsf=13.6
    psfrad=14.6
    fitrad=3.4
    threshold=3
    skystd=0.0006
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(915 798 611 993 1126)
    find_psf
}

#########
# HST Blue 0916
function hstb96a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b96a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.5
    ffwhmpsf=14.0
    psfrad=15.0
    fitrad=3.5
    threshold=3
    skystd=0.002
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(792 221 466 330 917 517 504)
    find_psf
}

#########
# HST Red 2003-01-05
function hstr03a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r03a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.7
    ffwhmpsf=10.8
    psfrad=11.8
    fitrad=2.7
    threshold=3
    skystd=0.015
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(999999)
    find_psf
}

#########
# HST Red 2003-08-12
function hstr03b {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r03b
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.7
    ffwhmpsf=10.8
    psfrad=11.8
    fitrad=2.7
    threshold=3
    skystd=0.02125
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(999999)
    find_psf
}

#########
# HST Red 2003-11-28
function hstr03c {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r03c
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.7
    ffwhmpsf=10.8
    psfrad=11.8
    fitrad=2.7
    threshold=3
    skystd=0.015
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(847 1270 775 606)
    find_psf
}

#########
# HST Red 2005-09-26
function hstr05a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r05a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.7
    ffwhmpsf=10.8
    psfrad=11.8
    fitrad=2.7
    threshold=3
    skystd=0.004
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(102 1216 1118 795 1157 1952)
    find_psf
}

#########
# HST Red 2006-04-15
function hstr06a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r06a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.6
    ffwhmpsf=10.4
    psfrad=11.4
    fitrad=2.6
    threshold=3
    skystd=0.0145
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(9999999)
    find_psf
}

#########
# HST Red 2006-04-29
function hstr06b {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r06b
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.3
    ffwhmpsf=9.2
    psfrad=10.2
    fitrad=2.3
    threshold=3
    skystd=0.03
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(9999999)
    find_psf
}

#########
# HST Red 2006-12-06
function hstr06c {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r06c
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.7
    ffwhmpsf=10.8
    psfrad=11.8
    fitrad=2.7
    threshold=3
    skystd=0.016
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(584 648 574 68 412 559 611)
    find_psf
}

#########
# HST Red 2007
function hstr07a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r07a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.3
    ffwhmpsf=13.2
    psfrad=14.2
    fitrad=3.3
    threshold=3
    skystd=0.00062
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(9999999999)
    find_psf
}

#########
# HST Red 2008
function hstr08a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r08a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.3
    ffwhmpsf=13.2
    psfrad=14.2
    fitrad=3.3
    threshold=3
    skystd=0.0007
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(9999999999)
    find_psf
}

#########
# HST Red 2009-04-29
function hstr09a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r09a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.3
    ffwhmpsf=13.2
    psfrad=14.2
    fitrad=3.3
    threshold=3
    skystd=0.0007
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(9999999)
    find_psf
}

#########
# HST Red 2009-12-12
function hstr09b {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r09b
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.3
    ffwhmpsf=13.2
    psfrad=14.2
    fitrad=3.3
    threshold=3
    skystd=0.005
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(643 401 502 241 435)
    find_psf
}

########
# HST Red 2011
function hstr11a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r11a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.6
    ffwhmpsf=14.4
    psfrad=15.4
    fitrad=3.6
    threshold=3
    skystd=0.00625
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(9999999)
    find_psf
}

########
# HST Red 2013
function hstr13a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r13a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.6
    ffwhmpsf=14.4
    psfrad=15.4
    fitrad=3.6
    threshold=3
    skystd=0.006875
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(9999999)
    find_psf
}

########
# HST Red 2014
function hstr14a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r14a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.3
    ffwhmpsf=13.2
    psfrad=14.2
    fitrad=3.3
    threshold=3
    skystd=0.0054
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(9999999)
    find_psf
}

########
# HST Red 2015
function hstr15a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r15a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.3
    ffwhmpsf=13.2
    psfrad=14.2
    fitrad=3.3
    threshold=3
    skystd=0.0055
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(753 681 1216 823 533)
    find_psf
}

########
# HST Red 2016
function hstr16a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r16a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.3
    ffwhmpsf=13.2
    psfrad=14.2
    fitrad=3.3
    threshold=3
    skystd=0.013
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(355 384 329 410 190 399 468 462 553 558 203 241 16 55 168 437)
    reject=(9999999)
    find_psf
}

########
# HST Red 0306
function hstr36a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r36a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=2.6
    ffwhmpsf=10.4
    psfrad=11.4
    fitrad=2.6
    threshold=3
    skystd=0.0033
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(91 1072 537 1278 1173 846)
    find_psf
}

########
# HST Red 0709
function hstr79a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r79a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.4
    ffwhmpsf=13.6
    psfrad=14.6
    fitrad=3.4
    threshold=3
    skystd=0.00033
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(855 713 518 235 572 621 965)
    find_psf
}

#########
# HST Red 1994
function hstr94a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r94a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.3
    ffwhmpsf=13.2
    psfrad=14.2
    fitrad=3.3
    threshold=3
    skystd=0.000875
    mindat=-0.1
    maxdat=35
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(1341 1073 1535 1674 1166)
    find_psf
}

########
# HST Red 0916
function hstr96a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r96a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.3
    ffwhmpsf=13.2
    psfrad=14.2
    fitrad=3.3
    threshold=3
    skystd=0.002
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(403 928 1195 1016 851 642 1594)
    find_psf
}

################################################################
# HST 225 2009 Reference pixfrac of 0.6
function hst22ra {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/22ra
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.4
    ffwhmpsf=13.6
    psfrad=14.6
    fitrad=3.4
    threshold=3
    skystd=0.0057
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(60 18 83 94 184 6 207 16)
    find_psf
}

#########
# HST 225 2009 Experimental pixfrac of 0.3
function hst22xa {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/22xa
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.3
    ffwhmpsf=13.2
    psfrad=14.2
    fitrad=3.3
    threshold=3
    skystd=0.0057
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(156 74 19 106 119 6 229 257)
    reject=(156 74 106 119 6 229 257)
#    reject=(9999999999)
    find_psf
}

#########
# HST 336 2009 Reference pixfrac of 0.6
function hst33ra {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/33ra
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.4
    ffwhmpsf=13.6
    psfrad=14.6
    fitrad=3.4
    threshold=3
    skystd=0.0055
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(144 380 64 182 173 136 24)
    find_psf
}

#########
# HST 336 2009 Experimental pixfrac of 0.3
function hst33xa {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/33xa
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.4
    ffwhmpsf=13.6
    psfrad=14.6
    fitrad=3.4
    threshold=3
    skystd=0.0055
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(297 164 73 444 211 154 25)
    find_psf
}

#########
# HST 555 2009 Reference pixfrac of 0.6
function hst55ra {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/55ra
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.3
    ffwhmpsf=13.2
    psfrad=14.2
    fitrad=3.3
    threshold=3
    skystd=0.012
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(271 139 320 309 261 344 96 151 535)
    find_psf
}

#########
# HST 555 2009 Experimental pixfrac of 0.3
function hst55xa {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/55xa
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.2
    ffwhmpsf=12.8
    psfrad=13.8
    fitrad=3.2
    threshold=3
    skystd=0.012
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(298 155 354 108 538 530 167 17 46)
    find_psf
}

#########
# HST 657 2016
function hst657a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/657a
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.4
    ffwhmpsf=13.6
    psfrad=14.6
    fitrad=3.4
    threshold=3
    skystd=0.0029
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(360 203 244 295 287 12 107 435 386 11)
    find_psf
}

#########
# HST 814 2009 Reference pixfrac of 0.6
function hst81ra {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/81ra
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.4
    ffwhmpsf=13.6
    psfrad=14.6
    fitrad=3.4
    threshold=3
    skystd=0.012
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(520 396 322 349 141 94 256 462 365)
    find_psf
}

#########
# HST 814 2009 Experimental pixfrac of 0.3
function hst81xa {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/81xa
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.4
    ffwhmpsf=13.6
    psfrad=14.6
    fitrad=3.4
    threshold=3
    skystd=0.012
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(569 427 351 107 157 338 281 378 504)
    find_psf
}

#########
# HST red 2015 Reference pixfrac of 0.7
function hstr15r {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r15r
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.3
    ffwhmpsf=13.2
    psfrad=14.2
    fitrad=3.3
    threshold=3
    skystd=0.0064
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(739 1207 807 526)
    find_psf
}

#########
# HST red 2015 Experimental pixfrac of 0.3
function hstr15x {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r15x
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.2
    ffwhmpsf=12.8
    psfrad=13.8
    fitrad=3.2
    threshold=3
    skystd=0.0064
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(851 885 1411 929 590 912)
    find_psf
}

#########
# HST blu 2015 Reference pixfrac of 0.7
function hstb15r {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b15r
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.6
    ffwhmpsf=14.4
    psfrad=15.4
    fitrad=3.6
    threshold=3
    skystd=0.005
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(839 468 333 55)
    find_psf
}

#########
# HST blu 2015 Experimental pixfrac of 0.5
function hstb15x {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b15x
    rm *.???.?*
    img=$(ls *.fits)
    fwhmpsf=3.5
    ffwhmpsf=14
    psfrad=15
    fitrad=3.5
    threshold=3
    skystd=0.005
    mindat=-0.1
    maxdat=80
    maxnpsf=24
    critovl=1
    maxgroup=60
    analytical=gauss
    varorder=1
    reject=(723 397 290 48)
    find_psf
}

################################################################
source activate iraf27
inter=inter-

#hst110a
#hst160a
#vlt163a
#vlt219a
#hst225a
#hst336a
#hst502a
#hst555a
#hst645a
#hst814a
#
#hstb03c
#hstb06a
#hstb09b
#hstb15a
#
#hstb36a
#hstb79a
#hstb94a
#hstb96a
#
#hstr03c
#hstr05a
#hstr06c
#hstr09b
#hstr15a
#
#hstr36a
#hstr79a
#hstr94a
#hstr96a
#
#hst22ra
#hst22xa
#hst33ra
#hst33xa
#hst55ra
#hst55xa
#hst657a
#hst81ra
#hst81xa
#hstr15r
#hstr15x
#hstb15r
hstb15x
