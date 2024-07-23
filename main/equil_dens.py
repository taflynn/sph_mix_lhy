import scipy.integrate as integrate
from scipy import optimize
from numpy import pi, sqrt, tan, cos, linspace, empty, arcsin
from main.units import units, natural_consts

from scipy.special import ellipeinc,ellipe,ellipkinc,ellipk
#from scipy.interpolate import CubicHermiteSpline

def eq_dens_lck(m1,m2,a11,a22,a12):
    if m1 == m2:
        # equilibrium densities (equal masses)
        n01 = (25*pi/1024)*(a12 + sqrt(a11*a22))**2/(a11**1.5*a22*(sqrt(a11) + sqrt(a22))**5)
        n02 = (25*pi/1024)*(a12 + sqrt(a11*a22))**2/(a22**1.5*a11*(sqrt(a11) + sqrt(a22))**5)

    elif m1 != m2:
        z = m2/m1

        g11 = 4*pi*a11/m1
        g22 = 4*pi*a22/m2
        g12 = 2*pi*a12*(1/m1 + 1/m2)

        # define interaction strengths
        deltag = g12 + sqrt(g11*g22)
        x = sqrt(g22/g11) 

        # calculate integral f using the t defined integral
        f = fzx_analytic(z, x)
        
        # equilibrium densities (unequal masses)
        n01 = (25*pi/1024)*(deltag**2/(f**2*a11**3*g11*g22))
        n02 = n01*sqrt(g11/g22)

    return n01, n02

def fzx_analytic(z, x):
    f3d_1 = (-2 - 7*x*z + 2*z**2 + x**2*z**2)*(sqrt(x+z)/(2*sqrt(z)*(z**2 - 1))) 
    f3d_2 = (-2 - 7*x*z + 3*z**2 + 3*x**2*z**2 - 7*x*z**3 - 2*x**2*z**4)*((ellipeinc(arcsin(1/z),-x*z)-ellipe(-x*z))/(2*(z**2 - 1)**1.5))
    f3d_3 = (2 + 8*x*z - 3*z**2 + 6*x**2*z**2 - 2*x*z**3 + x**2*z**4)*((ellipkinc(arcsin(1/z),-x*z)-ellipk(-x*z))/(2*(z**2 - 1)**1.5))
    return f3d_1 + f3d_2 + f3d_3

def dfdx_analytic(z, x):
    dfdx_l1 = (-7*z + 2*x*z**2)*sqrt(x + z)/(2*sqrt(z)*(z**2 - 1))
    dfdx_l2 = (-2 - 7*x*z + 2*z**2 + x**2*z**2)/(4*sqrt(z)*(z**2 - 1)*sqrt(x + z))
    
    dfdx_l3 = (-7*z + 6*x*z**2 - 7*z**3 - 4*x*z**4)*(ellipeinc(arcsin(1/z),-x*z) - ellipe(-x*z))
    dfdx_l4 = (-2 -7*x*z + 3*z**2 + 3*x**2*z**2 - 7*x*z**3 - 2*x**2*z**4)
    dfdx_l5 = ellipeinc(arcsin(1/z),-x*z) - ellipkinc(arcsin(1/z),-x*z) - ellipe(-x*z) + ellipk(-x*z)
    
    dfdx_l6 = (8*z + 12*x*z**2 - 2*z**3 + 2*x*z**4)*(ellipkinc(arcsin(1/z),-x*z) - ellipk(-x*z))
    dfdx_l7 = 2 + 8*x*z - 3*z**2 + 6*x**2*z**2 - 2*x*z**3 + x**2*z**4
    
    dfdx_l8 = ellipeinc(arcsin(1/z),-x*z) - (1+z*x)*ellipkinc(arcsin(1/z),-x*z) - ellipe(-x*z) + (1+z*x)*ellipk(-x*z)
    dfdx_l9 = sqrt(z**2 - 1)/(sqrt(z**2 + z*x))
    # combine lines
    poly_terms = dfdx_l1 + dfdx_l2
    first_ellip_terms = dfdx_l3 + (dfdx_l4*dfdx_l5)/(2*x)
    second_ellip_terms = dfdx_l6 + dfdx_l7*(dfdx_l8/x + dfdx_l9)/(2*(1 + z*x))
    
    total_terms = poly_terms + (first_ellip_terms + second_ellip_terms)/(2*(z**2 - 1)**1.5)
    
    return total_terms

def setup_f_interp(z):
    x = linspace(0.000001, 1.0, 100)
    fz_analytic_dat = empty(len(x))
    dfdx_analytic_dat = empty(len(x))

    for i in range(0, len(x)):
        fz_analytic_dat[i] = fzx_analytic(z, x[i])
        dfdx_analytic_dat[i] = dfdx_analytic(z, x[i])

    f_interp = CubicHermiteSpline(x, fz_analytic_dat, dfdx_analytic_dat)
    df_interp = f_interp.derivative()

    return f_interp, df_interp
