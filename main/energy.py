from numpy import pi, gradient, trapz, absolute

def energy_eqm_dens_lck(phi, V, r, dr, N_lck):
    """
    This function calculates the separate energy components of the density-locked mixture (w/ LHY), i.e.,
    * E_ke - Kinetic energy
    * E_pot - Potential energy due to trapping potential
    * E_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    * E_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    """
    dphi_dr = gradient(phi, dr)
    E_ke  = (1/N_lck)*4*pi*trapz(r**2*N_lck*0.5*absolute(dphi_dr)**2)*dr 
    E_pot = (1/N_lck)*4*pi*trapz(r**2*N_lck*V*absolute(phi)**2)*dr
    E_int = (1/N_lck)*4*pi*trapz(r**2*N_lck**2*(-1.5*absolute(phi)**4))*dr
    E_lhy = (1/N_lck)*4*pi*trapz(r**2*N_lck**2.5*absolute(phi)**5)*dr
                          
    return E_ke, E_pot, E_int, E_lhy

def energy_eqm_dens_ulck(psi1, psi2, r, V1, V2, dr, alpha, beta, eta, N1, N2):
    """
    This function calculates the separate energy components of the density-unlocked mixture (w/ LHY) for
    both components of the mixture.
    * E_ke - Kinetic energy 
    * E_pot - Potential energy due to trapping potential
    * E_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    * E_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    """
    dpsi1_dr = gradient(psi1, dr)
    dpsi2_dr = gradient(psi2, dr)
    E_coef = 4*alpha**(-1)/(3*pi**2)
    
    E_ke = E_coef*4*pi*trapz(r**2*(0.5*N1*absolute(dpsi1_dr)**2 + 0.5*N2*beta**(-1)*absolute(dpsi2_dr)**2))*dr
    E_pot = E_coef*4*pi*trapz(r**2*(N1*V1*absolute(psi1)**2 + N2*beta**(-1)*V2*absolute(psi2)**2))*dr
    E_int = E_coef*4*pi*trapz(r**2*(0.5*N1**2*absolute(psi1)**4 + 0.5*N2**2*absolute(psi2)**4 + N1*N2*eta*absolute(psi1)**2*absolute(psi2)**2))*dr
    E_lhy = E_coef*4*pi*trapz(r**2*(0.4*alpha*(N1*absolute(psi1)**2 + N2*beta*absolute(psi2)**2)**2.5))*dr
    
    return E_ke, E_pot, E_int, E_lhy

def energy_uneqm_dens_ulck(psi1, psi2, r, V1, V2, dr, gam, z, alpha, beta, eta, N1, N2):
    """
    This function calculates the separate energy components of the density-unlocked mixture (w/ LHY) for
    both components of the mixture. Thus for the i-th component we have the outputs:
    * E_ke - Kinetic energy 
    * E_pot - Potential energy due to trapping potential
    * E_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    * E_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    """
    dpsi1_dr = gradient(psi1, dr)
    dpsi2_dr = gradient(psi2, dr)
    E_coef = 4*alpha**(-1)/(3*pi**2)
    
    E_ke = E_coef*4*pi*trapz(r**2*(0.5*gam*N1*absolute(dpsi1_dr)**2 + 0.5*(gam/z)*N2*beta**(-1)*absolute(dpsi2_dr)**2))*dr
    E_pot = E_coef*4*pi*trapz(r**2*(N1*V1*absolute(psi1)**2 + N2*beta**(-1)*V2*absolute(psi2)**2))*dr
    E_int = E_coef*4*pi*trapz(r**2*(0.5*N1**2*absolute(psi1)**4 + 0.5*N2**2*absolute(psi2)**4 + N1*N2*eta*absolute(psi1)**2*absolute(psi2)**2))*dr
    E_lhy = E_coef*4*pi*trapz(r**2*(alpha*(N1*absolute(psi1)**2 + N2*z**0.6*beta*absolute(psi2)**2)**2.5))*dr
    
    return E_ke, E_pot, E_int, E_lhy
