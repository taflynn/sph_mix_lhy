from numpy import pi
import numpy as np

def grid_setup(dr,Nr,DT_COEF):
    """
    This function defines the radial space array, r, by taking as inputs:
    (1) dr - the spatial discretisation
    (2) Nr - the number of grid points
    and also defines a discrete time-step, dt, by allowing the user to input the coefficient in the 
    formula:
    
    dt = DT_COEF*dr**2
    
    This formula is commonly used to choose a time-step, though if DT_COEF >= 0.5 then the function
    will print a statement as this is commonly too big a time-step.
    """
    # box size
    Lr = Nr*dr
    
    # position array with 2 ghost points
    r = np.arange(-1/2,(Nr + 3/2),1)*dr
    
    if DT_COEF >= 0.5:
        print('Time step too big! Should be less than 0.5')
    
    # time step
    dt = DT_COEF*dr**2
    
    print('Initialising system...')
    print('Box size = ',Lr,' with dr = ',dr,' and dt = ',dt)
    print(67*'-')
    
    return Lr,r,dt

def init_wavefun_dens_lck(r,dr,GAUSS_SIGMA,INIT_TYPE):
    """
    This initialisation function defines the initial condition for imaginary time of the density-locked mixture. The
    function inputs are:
    (1) r - the spatial array
    (2) dr - tehe spatial discretisation
    (3) GAUSS_SIGMA - this is a numeric value which is the width parameter for a Gaussian
    (4) INIT_TYPE - this is a string input which allows the user to choose the function form of the initial condition
    
    INIT_TYPE can be specified as:
    -> 'GAUSS' - this gives a Gaussian form
    -> 'S_GAUSS' - this gives a super-Gaussian form (i.e. a Gaussian with a higher-order exponent)
    
    The super-Gaussian form can give a significant benefit in convergence time for the flat-topped droplet form. The 
    function only outputs phi, i.e., the initial wavefunction profile.
    """
    
    print('Density-locked mixture using initial condition:')
    
    if INIT_TYPE == 'GAUSS':
        # Gaussian initial condition
        phi = np.exp(-(r)**2/(2*(GAUSS_SIGMA)**2))
        print('Gaussian')
    elif INIT_TYPE == 'S_GAUSS':
        # Gaussian initial condition
        phi = np.exp(-((r)**2/(2*(GAUSS_SIGMA)**2))**3)
        print('Super-Gaussian')
    # normalise initial condition
    Norm = 4*pi*np.trapz(r**2*np.abs(phi)**2)*dr
    phi = phi/np.sqrt(Norm)
    print(67*'-')
    return phi

def init_wavefun_dens_ulck(r,dr,GAUSS_SIGMA,INIT_TYPE1,INIT_TYPE2):
    """
    This initialisation function defines the initial condition for imaginary time of the density-locked mixture. The
    function inputs are:
    (1) r - the spatial array
    (2) dr - tehe spatial discretisation
    (3) GAUSS_SIGMA - this is a numeric value which is the width parameter for a Gaussian
    (4) INIT_TYPE - this is a string input which allows the user to choose the function form of the initial condition
    
    INIT_TYPE can be specified as:
    -> 'GAUSS' - this gives a Gaussian form
    -> 'S_GAUSS' - this gives a super-Gaussian form (i.e. a Gaussian with a higher-order exponent)
    
    The super-Gaussian form can give a significant benefit in convergence time for the flat-topped droplet form. The 
    function only outputs phi, i.e., the initial wavefunction profile.
    """
    
    print('Density-unlocked mixture using initial condition:') 

    # choosing initial condition (1st-component)
    if INIT_TYPE1 == 'GAUSS':
        # Gaussian initial condition (1st-component)
        psi1 = np.exp(-(r)**2/(2*(GAUSS_SIGMA)**2)) 
        print('Gaussian in 1st component')
    elif INIT_TYPE1 == 'S_GAUSS':
        # Super-Gaussian initial condition
        psi1 = np.exp(-((r)**2/(2*(GAUSS_SIGMA)**2))**3)
        print('Super-Gaussian in 1st component')
    elif INIT_TYPE1 == 'NON_ZERO_TAIL':
        # Super-Gaussian initial condition w/ non-zero tail (1st-component)
        psi1 = np.exp(-((r)**2/(2*(GAUSS_SIGMA)**2))**3) + 0.001*np.ones(len(r))
        print('Super-Gaussian w/ non-zero density tail in 1st component')
    
    # normalise initial condition (1st-component)
    Norm = 4*pi*np.trapz(r**2*abs(psi1)**2)*dr
    psi1 = psi1/np.sqrt(Norm)

    # choosing initial condition (2nd-component)
    if INIT_TYPE2 == 'GAUSS':
        # Gaussian initial condition (2nd-component)
        psi2 = np.exp(-(r)**2/(2*(GAUSS_SIGMA)**2)) 
        print('Gaussian in 2nd component')
    elif INIT_TYPE2 == 'S_GAUSS':
        # Super-Gaussian initial condition (2nd-component)
        psi2 = np.exp(-((r)**2/(2*(GAUSS_SIGMA)**2))**3)
        print('Super-Gaussian in 2nd component')
    elif INIT_TYPE2 == 'NON_ZERO_TAIL':
        # Super-Gaussian initial condition w/ non-zero tail (2nd-component)
        psi2 = np.exp(-((r)**2/(2*(GAUSS_SIGMA)**2))**3) + 0.0001*np.ones(len(r))
        print('Super-Gaussian w/ non-zero density tail in 1st component')

    # normalise initial condition (2nd-component)
    Norm = 4*pi*np.trapz(r**2*abs(psi2)**2)*dr
    psi2 = psi2/np.sqrt(Norm) 
    print(67*'-')
    return psi1,psi2

def potential_dens_lck(r,OMEGA,tau):
    """
    In the density-locked regime we can only impose one trapping potential on the combined densities. 
    Here there is only one option and that is a spherically symmetric harmonic trap in which the 
    trapping can be varied by OMEGA. If no trap is desired, simply set OMEGA = 0.
    """
  
    OMEGA = 2*pi*OMEGA*tau
    
    V = 0.5*OMEGA**2*r**2
    
    if OMEGA == 0:
        print('No trapping potential on density-locked mixture')
    elif OMEGA != 0:
        print('Harmonic trapping potential with frequency, ',OMEGA,', on density-locked mixture')
    print(67*'-')
    return V

def potential_dens_ulck(r,OMEGA1,OMEGA2,tau):
    """
    In the density-unlocked regime we can imposed separate trapping potential on each component. 
    Here there is only one option and that is spherically symmetric harmonic traps in which the 
    trapping can be varied by OMEGA1 and OMEGA2. If no trap is desired, simply set OMEGA1 = OMEGA2 = 0.
    """
    OMEGA1 = 2*pi*OMEGA1*tau
    OMEGA2 = 2*pi*OMEGA2*tau
    
    V1 = 0.5*OMEGA1**2*r**2
    V2 = 0.5*OMEGA2**2*r**2
    
    if OMEGA1 == 0 and OMEGA2 == 0:
        print('No trapping potential on density-unlocked mixture')
    elif OMEGA1 == 0 and OMEGA2 != 0 :
        print('No trapping potential on component 1')
        print('and harmonic trapping potential with frequency, ',OMEGA2,', on component 2')
    elif OMEGA1 != 0 and OMEGA2 == 0 :
        print('Harmonic trapping potential with frequency, ',OMEGA1,', on component 1')
        print('and no trapping potential on component 2')
    elif OMEGA1 != 0 and OMEGA2 != 0:
        print('Harmonic trapping potential with frequency, ',OMEGA1,', on component 1')
        print('and harmonic trapping potential with frequency, ',OMEGA2,', on component 2')
    print(67*'-')
    return V1,V2