# Reddening based on model by cardelli89

import pdb

import numpy as np

def get_av(rv, ebv):
    return rv*ebv

def ir(lam, spec, rv, av):
    xx = 1e4/lam
    aa =  0.574*xx**1.61
    bb = -0.527*xx**1.61
    
    lo = 10000/1.1
    up = 10000/0.3
    in_range = (lam > lo) & (lam <= up)
#    pdb.set_trace()
#    if in_range.shape[0] == 0:
#        return
    
    al = av*(aa+bb/rv)
    spec[in_range] = 10**(-0.4*al[in_range])*spec[in_range]
 
def onir(lam, spec, rv, av):
    xx = 1e4/lam
    yy = (xx-1.82)
    aa = np.array([0.32999, -0.77530, 0.01979, 0.72085, -0.02427, -0.50447, 0.17699, 1])
    aa = np.polyval(aa, yy)
    bb = np.array([-2.09002, 5.30260, -0.62251, -5.38434, 1.07233, 2.28305, 1.41338, 0])
    bb = np.polyval(bb, yy)

    lo = 10000/3.3
    up = 10000/1.1
    in_range = (lam > lo) & (lam <= up)
#    if in_range.shape[0] == 0:
#        return
    
    al = av*(aa+bb/rv)
    spec[in_range] = 10**(-0.4*al[in_range])*spec[in_range]

def uv(lam, spec, rv, av):
    def Fa(xx):
        interval = (xx <= 8) & (xx >= 5.9)
        return interval*(-0.04473*(xx-5.9)**2-0.009779*(xx-5.9)**3)
    def Fb(xx):
        interval = (xx <= 8) & (xx >= 5.9)
        return interval*(0.2130*(xx-5.9)**2+0.1207*(xx-5.9)**3)
    
    xx = 1e4/lam
    aa = 1.752-0.316*xx-0.104/((xx-4.67)**2+0.341)+Fa(xx)
    bb = -3.090+1.825*xx+1.206/((xx-4.62)**2+0.263)+Fb(xx)

    lo = 10000/8.
    up = 10000/3.3
    in_range = (lam > lo) & (lam <= up)
#    if in_range.shape[0] == 0:
#        return
    
    al = av*(aa+bb/rv)
    spec[in_range] = 10**(-0.4*al[in_range])*spec[in_range]
    
#    import matplotlib.pyplot as plt
#    plt.plot(xx[in_range],al[in_range]/av)
#    plt.show()
#    pdb.set_trace()


def red_mod(lam, spec, rv, ebv):
    av = get_av(rv, ebv)
    ir(lam, spec, rv, av)
    onir(lam, spec, rv, av)
    uv(lam, spec, rv, av)

def cor_red(lam, rv, ebv):
    if lam.size == 1:
        lam = np.ones(2)*lam
        spec = np.ones(2)
        vector = False
    else:
        spec = np.ones(lam.size)
        vector = True
        
    av = get_av(rv, ebv)
    ir(lam, spec, rv, av)
    onir(lam, spec, rv, av)
    uv(lam, spec, rv, av)

    return 1/spec if vector else 1/spec[0]
