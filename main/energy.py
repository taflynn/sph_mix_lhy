import numpy as np
from numpy import pi
def energy_eqm_dens_lck(phi,V,r,dr,N_lck):
    """
    This function calculates the separate energy components of the density-locked mixture (w/ LHY), i.e.,
    (1) E_ke - Kinetic energy
    (2) E_pot - Potential energy due to trapping potential
    (3) E_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    (4) E_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    
    This requires the inputs:
    -> phi - the wavefunction
    -> V - the trapping potential
    -> r - the spatial array
    -> dr- the spatial discretisation
    -> N_lck - the effective atom number used to normalise wavefunction
    """
    dphi_dr = np.gradient(phi,dr)
    E_ke  = (1/N_lck)*4*pi*np.trapz(r**2*0.5*np.abs(dphi_dr)**2)*dr 
    E_pot = (1/N_lck)*4*pi*np.trapz(r**2*V*np.abs(phi)**2)*dr
    E_int = (1/N_lck)*4*pi*np.trapz(r**2*(-1.5*N_lck*np.abs(phi)**4))*dr
    E_lhy = (1/N_lck)*4*pi*np.trapz(r**2*(N_lck**1.5*np.abs(phi)**5))*dr
                          
    return E_ke,E_pot,E_int,E_lhy

def energy_eqm_dens_ulck(psi1,psi2,r,V1,V2,dr,alpha,beta,eta,N1,N2):
    """
    This function calculates the separate energy components of the density-unlocked mixture (w/ LHY) for
    both components of the mixture. Thus for the i-th component we have the outputs:
    (1) Ei_ke - Kinetic energy 
    (2) Ei_pot - Potential energy due to trapping potential
    (3) Ei_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    (4) Ei_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    
    This requires the inputs:
    -> psi1,psi2 - the wavefunctions
    -> V1,V2 - the trapping potentials
    -> r - the spatial array
    -> dr- the spatial discretisation
    -> alpha,beta,eta - the dimensionless parameters defining the Gross-Pitaevskii equations
    -> N1,N2 - the rescaled atom numbers used to normalise wavefunctions
    """
    dpsi1_dr = np.gradient(psi1,dr)
    dpsi2_dr = np.gradient(psi2,dr)
    E_coef = 4*alpha**(-1)/(3*pi**2)
    
    E_ke = E_coef*4*pi*np.trapz(r**2*(0.5*np.abs(dpsi1_dr)**2 + 0.5*beta**(-1)*np.abs(dpsi2_dr)**2))*dr
    E_pot = 0.0
    E_int = E_coef*4*pi*np.trapz(r**2*(0.5*np.abs(psi1)**4 + 0.5*np.abs(psi2)**4 + eta*np.abs(psi1)**2*np.abs(psi2)**2))*dr
    E_lhy = E_coef*4*pi*np.trapz(r**2*(0.4*alpha*(np.abs(psi1)**2 + beta*np.abs(psi2)**2)**2.5))*dr
    
    return E_ke,E_pot,E_int,E_lhy
