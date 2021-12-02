import numpy as np
import scipy.integrate as integrate
from numpy import pi
from units import units,natural_consts

def eq_dens_lck(m1,m2,a11,a22,a12):
    if m1 == m2:
        # equilibrium densities (equal masses)
        n01 = (25*pi/1024)*(a12+np.sqrt(a11*a22))**2/(a11**1.5*a22*(np.sqrt(a11) + np.sqrt(a22))**5)
        n02 = (25*pi/1024)*(a12+np.sqrt(a11*a22))**2/(a22**1.5*a11*(np.sqrt(a11) + np.sqrt(a22))**5)

    elif m1 != m2:
        [hbar,a0,Da] = natural_consts()
        # define interaction strengths
        g11 = 4*pi*hbar**2*a11/m1
        g12 = 2*pi*hbar**2*a12*(1/m1 + 1/m2)
        g22 = 4*pi*hbar**2*a22/m2
        deltag = g12 + np.sqrt(g11*g22)
        
        # calculate integral f using the t defined integral
        f = (15/32)*integrate.quadrature(lambda t,z,x: Ft(t,z,x),0,pi/2,args=(m2/m1,np.sqrt(g22/g11)),rtol=1e-06,maxiter=5000)[0]
        
        # equilibrium densities (unequal masses)
        n01 = (25*pi/1024)*(deltag**2/(f**2*a11**3*g11*g22))
        n02 = n01*np.sqrt(g11/g22)
    return n01,n02

def Ft(t,z,x):
    """
    To calculate the equilibrium densities of the density-locked mixture with a mass imbalance requires the computing of the
    integral within the function 'f' in the Bose-Bose Lee-Huang-Yang (LHY) term. To solve this integral we use Gaussian
    quadrature and so this function simply defines the integrand of this integral using the tranposition k -> tan(t).
    """
    pet_eq = np.cos(t)**(-2)*(np.tan(t)**2*np.sqrt(0.5*(np.tan(t)**2*(1+x/z) + 0.25*np.tan(t)**4*(1+z**(-2))) + np.sqrt(0.25*((np.tan(t)**2 + 0.25*np.tan(t)**4)-((x/z)*np.tan(t)**2 + 0.25*z**(-2)*np.tan(t)**4))**2 + (x/z)*np.tan(t)**4)) + np.tan(t)**2*np.sqrt(0.5*(np.tan(t)**2*(1+x/z) + 0.25*np.tan(t)**4*(1+z**(-2))) - np.sqrt(0.25*((np.tan(t)**2 + 0.25*np.tan(t)**4)-((x/z)*np.tan(t)**2 + 0.25*z**(-2)*np.tan(t)**4))**2 + (x/z)*np.tan(t)**4))-((1+z)/(2*z))*np.tan(t)**4 - (1+x)*np.tan(t)**2 + (1/(1+z))*((1+x*z)**2 + z*(1+x)**2))
    return pet_eq

def eq_dens_unlck(m1,m2,a11,a22,a12):
    """
    This function calculates the equilibrium densities for each component of the density-unlocked mixture by 
    solving the coupled equations defined in 'fun' function. To solve these coupled equations the 'optimize.root' 
    function is used. This requires:
    (1) The coupled equations in 'fun'
    (2) An initial guess which here are given as the associated density-locked equilibrium densities
    (3) The function also requires an associated Jacobian for efficiency and accuracy of the'hypr' algorithm
    
    In the future this function will also calculate the equilibrium densities for the unequal masses density-unlocked
    mixture.
    """
    if m1 == m2:
        n01_lck = (25*pi/1024)*(a12+np.sqrt(a11*a22))**2/(a11**1.5*a22*(np.sqrt(a11) + np.sqrt(a22))**5)
        n02_lck = (25*pi/1024)*(a12+np.sqrt(a11*a22))**2/(a22**1.5*a11*(np.sqrt(a11) + np.sqrt(a22))**5)
        sol = optimize.root(fun, [n01_lck,n02_lck],args=(a11,a22,a12), jac = jac, method='hybr')
    elif m1!= m2:
        print('Unequal masses, density-unlocked is not ready yet')
    return n01,n02

def fun(x,a11,a22,a12):
    """
    To find the two-component equilibrium densities we must solve the coupled equations specified in this function.
    """
    return [(((a11*x[0] + a12*x[1])-((0.5*a11*x[0]**2 + a12*x[0]*x[1] + 0.5*a22*x[1]**2)/(x[0] + x[1]))) \
            +((64/(15*np.sqrt(pi)))*(a11*x[0] + a22*x[1])**1.5*(2.5*a11 - (a11*x[0] + a22*x[1])/(x[0] + x[1])))), \
            (((a22*x[1] + a12*x[0])-((0.5*a11*x[0]**2 + a12*x[0]*x[1] + 0.5*a22*x[1]**2)/(x[0] + x[1]))) \
            +((64/(15*np.sqrt(pi)))*(a11*x[0] + a22*x[1])**1.5*(2.5*a22 - (a11*x[0] + a22*x[1])/(x[0] + x[1]))))]

def jac(x,a11,a22,a12):
    """
    To find the two-component equilibrium densities we must solve the coupled equations in the 'fun' functions. This function
    defines the associated Jacobian to these functions in order for the numerical solver to have an improved method for solving
    the coupled equations.
    """
    return [[(x[0] + x[1])**(-2)*(0.5*a11*x[0]**2 + a11*x[0]*x[1] + (a11 -a12 + 0.5*a22)*x[1]**2 + (64/(15*np.sqrt(pi)))*(a11*x[0] + a22*x[1])**0.5*(1.5*a11*(x[0] + x[1])*(1.5*a11*x[0] + (2.5*a11 - a22)*x[1]) + (a11*x[0] + a22*x[1])*(a22*x[1] - a11*x[1]))),\
            (x[0] + x[1])**(-2)*(0.5*a11*x[0]**2 + (2*a12 - a22)*x[0]*x[1] - 0.5*a22*x[1]**2 + (64/(15*np.sqrt(pi)))*(a11*x[0] + a22*x[1])**0.5*(1.5*a22*(x[0] + x[1])*(1.5*a11*x[0] + (2.5*a11 - a22)*x[1]) + (a11*x[0] + a22*x[1])*(a11*x[0] - a22*x[0])))], \
            [(x[0] + x[1])**(-2)*(0.5*a22*x[1]**2 + (2*a12 - a11)*x[0]*x[1] - 0.5*a11*x[0]**2 + (64/(15*np.sqrt(pi)))*(a11*x[0] + a22*x[1])**0.5*(1.5*a11*(x[0] + x[1])*((2.5*a22 - a11)*x[0] + 1.5*a22*x[1]) + (a11*x[0] + a22*x[1])*(a22*x[1] - a11*x[1]))), \
            (x[0] + x[1])**(-2)*(0.5*a22*x[1]**2 + a22*x[0]*x[1] + (a22 -a12 + 0.5*a11)*x[0]**2 + (64/(15*np.sqrt(pi)))*(a11*x[0] + a22*x[1])**0.5*(1.5*a22*(x[0] + x[1])*((2.5*a22 - a11)*x[0] + 1.5*a22*x[1]) + (a11*x[0] + a22*x[1])*(a22*x[0] - a11*x[0])))]]