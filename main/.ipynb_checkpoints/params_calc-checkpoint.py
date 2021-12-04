import numpy as np
from numpy import pi
from main.equil_dens import eq_dens_lck
from main.units import units,natural_consts

def params_dens_lck(m1,m2,a11,a22,a12,N):
    """
    This function calculates the necessary dimensionless parameter for the density-locked mixture. 
    
    The outputs are hence:
    (1) N_lck - the dimensionless effective atom number appearing in the Gross-Pitaevskii (GP) equation
    (2) xi,tau - the length and timescales that define the system
    (3) n01,n02 - the equilibrium densities for each respective component
    
    The only inputs thar are neeeded are
    -> m1,m2 - the two component masses
    -> a11,a22 - the two intraspecies scattering lengths
    -> a12 - the interspecies scattering length
    -> N1,N2 - the two component atom numbers
    """
    [hbar,a0,Da] = natural_consts()
    a11 = a11*a0; a22 = a22*a0; a12 = a12*a0
    m1 = m1*Da; m2 = m2*Da
    # equilibrium densities equal masses
    n01,n02 = eq_dens_lck(m1,m2,a11,a22,a12)
    
    # length and time scales
    xi,tau = units(m1,m2,a11,a22,a12,n01)
    
    # density-locked effective atom number
    if m1 == m2:
        N_lck = (N/(n01*xi**3))*(np.sqrt(a22)/(np.sqrt(a11) + np.sqrt(a22)))
    elif m1 != m2:
        # interaction strengths
        g11 = 4*pi*hbar**2*a11/m1
        g12 = 2*pi*hbar**2*a12*(1/m1 + 1/m2)
        g22 = 4*pi*hbar**2*a22/m2
        N_lck = (N/(n01*xi**3))*(np.sqrt(g22)/(np.sqrt(g11) + np.sqrt(g22)))
    
    print('Calculating defining parameters and scales for density locked mixture:')
    print('Effective atom number of density-locked mixture, N_lck = ',N_lck)
    print('Lengthscale, xi = ',xi)
    print('Timescale, tau = ',tau)
    print('Equilibrium density of component 1, n01 = ',n01)
    print('Equilibrium density of component 2, n02 = ',n02)
    print(67*'-')
    return N_lck,xi,tau,n01,n02

def params_dens_ulck(m1,m2,a11,a22,a12,N1,N2):
    """
    This function calculates the necessary dimensionless parameters for the density-unlocked mixture (though this is 
    currently only for the equal masses mixture). The outputs are hence:
    (1) alpha,beta,eta - the three dimensionless parameters appearing in the Gross-Pitaevskii (GP) equation
    (2) xi,tau - the length and timescales that define the system
    (3) n01,n02 - the equilibrium densities for each component
    (4) rho1,rho2 - the density scalings for the two component wavefunctions
    (5) N1,N2 - the rescaled atom numbers for each component
    
    The only inputs thar are neeeded are
    -> m1,m2 - the two component masses
    -> a11,a22 - the two intraspecies scattering lengths
    -> a12 - the interspecies scattering length
    -> N1,N2 - the two component atom numbers
    """
    
    # natural constants
    [hbar,a0,Da] = natural_consts()
    
    if m1 == m2:
        a11 = a11*a0; a22 = a22*a0; a12 = a12*a0
        m1 = m1*Da; m2 = m2*Da
        # equilibrium densities equal masses solving coupled equations
        #n01,n02 = eq_dens_unlck_eqm(m1,m2,a11,a22,a12)
        
        # equilibrium densities equal masses with density-lock
        n01,n02 = eq_dens_lck(m1,m2,a11,a22,a12)
        
        delta_a = a12 + np.sqrt(a11*a22)
        
        alpha = (32/3)*np.sqrt((2/(3*pi))*(np.abs(delta_a)*a11**2.5*n01)/(np.sqrt(a11) + np.sqrt(a22)))
        beta = np.sqrt(a22/a11)
        eta = a12/np.sqrt(a11*a22)
        
        rho1 = (2/3)*((np.abs(delta_a)*n01)/(np.sqrt(a11)*(np.sqrt(a11) + np.sqrt(a22))))
        rho2 = (2/3)*((np.abs(delta_a)*n01)/(np.sqrt(a22)*(np.sqrt(a11) + np.sqrt(a22))))
        
        xi,tau = units(m1,m2,a11,a22,a12,n01)
        
        N1 = N1/(rho1*xi**3)
        N2 = N2/(rho2*xi**3)
        
    elif m1 != m2:
        # equilibrium densities unequal masses
        print('Density-unlocked and unequal masses is not yet available')
    
    print('Calculating defining parameters and scales for density-unlocked mixture:')
    print('Dimensionless parameters of mixture:')
    print('alpha = ',alpha, ', beta = ',beta,', eta = ',eta)
    print('Lengthscale, xi = ',xi)
    print('Timescale, tau = ',tau)
    print('Equilibrium density of component 1, n01 = ',n01)
    print('Equilibrium density of component 2, n02 = ',n02)
    print(67*'-')
    
    return alpha,beta,eta,xi,tau,n01,n02,rho1,rho2,N1,N2