from numpy import pi, sqrt, absolute

def units(m1, m2, a11, a22, a12, n01):
    """
    time and length scales
    """
    # natural constants
    [hbar,a0,Da] = natural_consts()

    if m1 == m2:
        m = m1
        delta_a = a12 + sqrt(a11*a22)
        
        # length and time scaling
        xi = sqrt((3/(8*pi))*((sqrt(a11) + sqrt(a22))/(absolute(delta_a)*sqrt(a11)*n01)))
        tau = (3*m/(8*pi*hbar))*((sqrt(a11) + sqrt(a22))/(absolute(delta_a)*sqrt(a11)*n01))
        
    elif m1 != m2:
        # interaction strengths
        g11 = 4*pi*hbar**2*a11/m1
        g12 = 2*pi*hbar**2*a12*(1/m1 + 1/m2)
        g22 = 4*pi*hbar**2*a22/m2
        deltag = g12 + sqrt(g11*g22)
        
        # length and time scaling
        xi = hbar*sqrt(1.5*(sqrt(g22)/m1 + sqrt(g11)/m2)/(absolute(deltag)*sqrt(g11)*n01))
        tau = hbar*1.5*(sqrt(g11) + sqrt(g22))/(absolute(deltag)*sqrt(g11)*n01)
        
    return xi, tau

def natural_consts():
    """
    The natural constants used throughout defining the parameters in these equations:
    -> hbar - Planck's constant on 2\pi
    -> a0 - Bohr radius used in defining the scattering length
    -> Da - Dalton used in defining the unified atomic mass
    """
    # natural constants
    hbar = 1.054571817*1e-34
    a0 = 5.29177210903*1e-11
    Da = 1.66053906660*1e-27
    return hbar, a0, Da
