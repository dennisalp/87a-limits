#!/bin/bash
# Dennis Alp 2016-11-30
# Find the daofind threshold parameter by plotting number
# of stars as a function of threshold.
# This script is closely related to phd_87a_uplim_find_psf.sh

# Display the image with display in ds9, then daoedit
# to get interactive FWHM-plot using the i key.
# imexamine for min, max, and sigma (skystd). 

function find_find {
    printf "!echo ECHOES ARE COMMENTS
iraf.Verbose=2
daophot
" > pyraf.cmd

    num=$(awk 'BEGIN{for(i=1;i<=5;i+=0.125)print i}')
    for threshold in $num
    do
	printf "daofind $img sts_$threshold.coo
$fwhmpsf
$skystd
$threshold
$mindat
$maxdat
!printf \"$threshold \$(awk 'END {print \$NF}' sts_$threshold.coo)"'\\n'"\" >> nstars.dat\n\n" >> pyraf.cmd
    done
    
    printf ".exit" >> pyraf.cmd
    pyraf < pyraf.cmd

    printf "import numpy as np
import matplotlib.pyplot as plt
dat = np.loadtxt(\"nstars.dat\")
plt.plot(dat[:,0], dat[:,1],'ok')
plt.xlabel('Threshold [sigma]')
plt.ylabel('Number of detections')
plt.savefig('nstars.pdf',bbox_inches='tight', pad_inches=0.01)
" > python.cmd
    python python.cmd
}

#########
# Help function, clean up the folder and prepare
function prep {
    printf "" > nstars.dat
    rm *.coo
    img=$(ls *.fits)
    }

#########
# HST 110 2011
function hst110a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/110a
    prep
    fwhmpsf=6.5
    skystd=0.016
    mindat=-0.1
    maxdat=200
    find_find
}

#########
# HST 160 2011
function hst160a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/160a
    prep
    fwhmpsf=7.2
    skystd=0.024
    mindat=-0.1
    maxdat=200
    find_find
}

#########
# VLT H 2010
function vlt163a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/163a
    prep
    fwhmpsf=4
    skystd=2.2
    mindat=-0.1
    maxdat=8000
    find_find
}

#########
# VLT K 2012
function vlt219a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/219a
    prep
    fwhmpsf=3.4
    skystd=2.2
    mindat=-0.1
    maxdat=8000
    find_find
}

#########
# HST 225 2009
function hst225a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/225a
    prep
    fwhmpsf=3.7
    skystd=0.0047
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST 336 2009
function hst336a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/336a
    prep
    fwhmpsf=3.7
    skystd=0.0045
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST 502 2016
function hst502a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/502a
    prep
    fwhmpsf=3.6
    skystd=0.0033
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST 555 2009
function hst555a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/555a
    prep
    fwhmpsf=3.5
    skystd=0.012
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST 645 2016
function hst645a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/645a
    prep
    fwhmpsf=3.3
    skystd=0.0022
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST 814 2009
function hst814a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/814a
    prep
    fwhmpsf=3.7
    skystd=0.01
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2003
function hstb03a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b03a
    prep
    fwhmpsf=2.4
    skystd=0.01
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2003
function hstb03b {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b03b
    prep
    fwhmpsf=2.4
    skystd=0.01375
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2003
function hstb03c {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b03c
    prep
    fwhmpsf=2.4
    skystd=0.011
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2004
function hstb04a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b04a
    prep
    fwhmpsf=2.4
    skystd=0.01
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2005
function hstb05a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b05a
    prep
    fwhmpsf=2.4
    skystd=0.01
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST Blue 2006
function hstb06a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b06a
    prep
    fwhmpsf=2.5
    skystd=0.01
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2007
function hstb07a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b07a
    prep
    fwhmpsf=3.5
    skystd=0.00075
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2008
function hstb08a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b08a
    prep
    fwhmpsf=3.5
    skystd=0.00067
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST Blue 2009
function hstb09a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b09a
    prep
    fwhmpsf=3.5
    skystd=0.00058
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST Blue 2009
function hstb09b {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b09b
    prep
    fwhmpsf=3.25
    skystd=0.005
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2011
function hstb11a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b11a
    prep
    fwhmpsf=3.6
    skystd=0.0065
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2013
function hstb13a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b13a
    prep
    fwhmpsf=3.5
    skystd=0.00375
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2014
function hstb14a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b14a
    prep
    fwhmpsf=3.5
    skystd=0.003375
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2015
function hstb15a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b15a
    prep
    fwhmpsf=3.6
    skystd=0.0045
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2016
function hstb16a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b16a
    prep
    fwhmpsf=3.4
    skystd=0.00975
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST Blue 0306
function hstb36a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b36a
    prep
    fwhmpsf=2.4
    skystd=0.004
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST Blue 0709
function hstb79a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b79a
    prep
    fwhmpsf=3.4
    skystd=0.000233
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 1994
function hstb94a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b94a
    prep
    fwhmpsf=3.4
    skystd=0.0006
    mindat=-0.1
    maxdat=35
    find_find
}

#########
# HST Blue 0916
function hstb96a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b96a
    prep
    fwhmpsf=3.5
    skystd=0.002
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2003
function hstr03a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r03a
    prep
    fwhmpsf=2.7
    skystd=0.015
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2003
function hstr03b {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r03b
    prep
    fwhmpsf=2.7
    skystd=0.02125
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2003
function hstr03c {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r03c
    prep
    fwhmpsf=2.7
    skystd=0.015
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2005
function hstr05a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r05a
    prep
    fwhmpsf=2.7
    skystd=0.004
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2006
function hstr06a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r06a
    prep
    fwhmpsf=2.6
    skystd=0.0145
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2006
function hstr06b {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r06b
    prep
    fwhmpsf=2.3
    skystd=0.03
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2006
function hstr06c {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r06c
    prep
    fwhmpsf=2.7
    skystd=0.016
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2007
function hstr07a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r07a
    prep
    fwhmpsf=3.3
    skystd=0.00062
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2008
function hstr08a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r08a
    prep
    fwhmpsf=3.3
    skystd=0.0007
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2009
function hstr09a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r09a
    prep
    fwhmpsf=3.3
    skystd=0.0007
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2009
function hstr09b {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r09b
    prep
    fwhmpsf=3.3
    skystd=0.005
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2011
function hstr11a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r11a
    prep
    fwhmpsf=3.6
    skystd=0.00625
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2013
function hstr13a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r13a
    prep
    fwhmpsf=3.6
    skystd=0.006875
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2014
function hstr14a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r14a
    prep
    fwhmpsf=3.3
    skystd=0.0054
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2015
function hstr15a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r15a
    prep
    fwhmpsf=3.3
    skystd=0.0055
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2016
function hstr16a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r16a
    prep
    fwhmpsf=3.3
    skystd=0.013
    mindat=-0.1
    maxdat=80
    find_find
}

########
# HST Red 0306
function hstr36a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r36a
    prep
    fwhmpsf=2.6
    skystd=0.0033
    mindat=-0.1
    maxdat=80
    find_find
}

########
# HST Red 0709
function hstr79a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r79a
    prep
    fwhmpsf=3.4
    skystd=0.00033
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 1994
function hstr94a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r94a
    prep
    fwhmpsf=3.3
    skystd=0.000875
    mindat=-0.1
    maxdat=35
    find_find
}

########
# HST Red 0916
function hstr96a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r96a
    prep
    fwhmpsf=3.3
    skystd=0.002
    mindat=-0.1
    maxdat=80
    find_find
}



################################################################
#########
# HST 225 2009 Reference pixfrac of 0.6
function hst22ra {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/22ra
    prep
    fwhmpsf=3.4
    skystd=0.0057
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST 225 2009 Experimental pixfrac of 0.3
function hst22xa {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/22xa
    prep
    fwhmpsf=3.3
    skystd=0.0057
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST 336 2009 Reference pixfrac of 0.6
function hst33ra {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/33ra
    prep
    fwhmpsf=3.4
    skystd=0.0055
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST 336 2009 Experimental pixfrac of 0.3
function hst33xa {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/33xa
    prep
    fwhmpsf=3.4
    skystd=0.0055
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST 555 2009 Reference pixfrac of 0.6
function hst55ra {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/55ra
    prep
    fwhmpsf=3.3
    skystd=0.012
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST 555 2009 Experimental pixfrac of 0.3
function hst55xa {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/55xa
    prep
    fwhmpsf=3.2
    skystd=0.012
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST 657 2016
function hst657a {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/657a
    prep
    fwhmpsf=3.4
    skystd=0.0029
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST 814 2009 Reference pixfrac of 0.6
function hst81ra {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/81ra
    prep
    fwhmpsf=3.4
    skystd=0.012
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST 814 2009 Experimental pixfrac of 0.3
function hst81xa {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/81xa
    prep
    fwhmpsf=3.4
    skystd=0.012
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2015 Reference pixfrac of 0.7
function hstr15r {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r15r
    prep
    fwhmpsf=3.3
    skystd=0.0064
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST red 2015 Experimental pixfrac of 0.3
function hstr15x {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/r15x
    prep
    fwhmpsf=3.2
    skystd=0.0064
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2015 Reference pixfrac of 0.7
function hstb15r {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b15r
    prep
    fwhmpsf=3.6
    skystd=0.005
    mindat=-0.1
    maxdat=80
    find_find
}

#########
# HST blu 2015 Experimental pixfrac of 0.5
function hstb15x {
    cd /Users/silver/Dropbox/phd/projects/87a/uplim/psf/b15x
    prep
    fwhmpsf=3.5
    skystd=0.005
    mindat=-0.1
    maxdat=80
    find_find
}



################################################################
source activate iraf27

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
#hstb03a
#hstb03b
#hstb03c
#hstb04a
#hstb05a
#hstb06a
#hstb07a
#hstb08a
#hstb09a
#hstb09b
#hstb11a
#hstb13a
#hstb14a
#hstb15a
#hstb16a
#
#hstb36a
#hstb79a
#hstb94a
#hstb96a
#
#hstr03a
#hstr03b
#hstr03c
#hstr05a
#hstr06a
#hstr06b
#hstr06c
#hstr07a
#hstr08a
#hstr09a
#hstr09b
#hstr11a
#hstr13a
#hstr14a
#hstr15a
#hstr16a
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
hstr15r
hstr15x
#hstb15r
#hstb15x
