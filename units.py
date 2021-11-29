import numpy as np
from numpy import pi

def units(m1,m2,a11,a22,a12,n01):
    # natural constants
    [hbar,a0,Da] = natural_consts()

    if m1 == m2:
        m = m1
        delta_a = a12 + np.sqrt(a11*a22)
        
        # length and time scaling
        xi = np.sqrt((3/(8*pi))*((np.sqrt(a11) + np.sqrt(a22))/(np.abs(delta_a)*np.sqrt(a11)*n01)))
        tau = (3*m/(8*pi))*((np.sqrt(a11) + np.sqrt(a22))/(np.abs(delta_a)*np.sqrt(a11)*n01))
        
    elif m1 != m2:
        # interaction strengths
        g11 = 4*pi*hbar**2*a11/m1
        g12 = 2*pi*hbar**2*a12*(1/m1 + 1/m2)
        g22 = 4*pi*hbar**2*a22/m2
        deltag = g12 + np.sqrt(g11*g22)
        
        # length and time scaling
        xi = hbar*np.sqrt(1.5*(np.sqrt(g22)/m1 + np.sqrt(g11)/m2)/(np.abs(deltag)*np.sqrt(g11)*n01))
        tau = hbar*1.5*(np.sqrt(g11) + np.sqrt(g22))/(np.abs(deltag)*np.sqrt(g11)*n01)
        
    return xi,tau

def natural_consts():
    # natural constants
    hbar = 1.0545718*1e-34
    a0 = 5.29*1e-11 # Bohr radius
    Da = 1.66053906660*1e-27
    return hbar,a0,Da