from numpy import gradient, trapz, absolute, conj

def mu_dens_lck(phi, r, V, dr, N_lck):
    """
    This function calculates mu via an integration which is defined as the integral of the \phi^* H \phi. 
    """
    dphi_dr = gradient(phi, dr)
    mu = trapz(r**2*(0.5*absolute(dphi_dr)**2
                     + V*absolute(phi)**2
                     - 3*N_lck*absolute(phi)**4
                     + 2.5*N_lck**1.5*absolute(phi)**5)) \
         / trapz(r**2*absolute(phi)**2)
    return mu

def mu_eqm_dens_ulck(psi1, psi2, r, V1, V2, dr, alpha, beta, eta, N1, N2):
    """
    This function calculates both mu1 and mu2 via an integration which is defined as the integral of the 
    \psi1^* H1 \psi1 and \psi2^* H2 \psi2. 
    """
    dpsi1_dr = gradient(psi1, dr)
    mu1 = trapz(r**2*(0.5*absolute(dpsi1_dr)**2 
                      + V1*absolute(psi1)**2
                      + N1*absolute(psi1)**4
                      + N2*eta*absolute(psi1)**2*absolute(psi2)**2
                      + alpha*(N1*absolute(psi1)**2 
                      + N2*beta*absolute(psi2)**2)**1.5*absolute(psi1)**2)) \
           /trapz(r**2*absolute(psi1)**2)
    
    dpsi2_dr = gradient(psi2, dr)
    mu2 = trapz(r**2*(0.5*absolute(dpsi2_dr)**2
                      + V2*absolute(psi2)**2
                      + N2*beta*absolute(psi2)**4
                      + N1*eta*beta*absolute(psi1)**2*absolute(psi2)**2
                      + alpha*beta**2*(N1*absolute(psi1)**2 
                      + N2*beta*absolute(psi2)**2)**1.5*absolute(psi2)**2)) \
           /trapz(r**2*absolute(psi2)**2)
    return mu1, mu2

def mu_uneqm_dens_ulck(psi1, psi2, r, V1, V2, dr, gam, z, alpha, beta, eta, N1, N2):
    """
    This function calculates both mu1 and mu2 via an integration which is defined as the integral of the 
    \psi1^* H1 \psi1 and \psi2^* H2 \psi2. 
    """
    dpsi1_dr = gradient(psi1, dr)
    mu1 = trapz(r**2*(0.5*gam*absolute(dpsi1_dr)**2
                      + V1*absolute(psi1)**2
                      + N1*absolute(psi1)**4
                      + N2*eta*absolute(psi1)**2*absolute(psi2)**2
                      + 2.5*alpha*(N1*absolute(psi1)**2 
                      + N2*beta*z**0.6*absolute(psi2)**2)**1.5*absolute(psi1)**2)) \
          /trapz(r**2*absolute(psi1)**2)
    
    dpsi2_dr = gradient(psi2, dr)
    mu2 = trapz(r**2*(0.5*(gam/z)*absolute(dpsi2_dr)**2
                      + V2*absolute(psi2)**2
                      + N2*beta*absolute(psi2)**4
                      + N1*eta*beta*absolute(psi1)**2*absolute(psi2)**2
                      + 2.5*alpha*beta**2.0*z**0.6*(N1*absolute(psi1)**2 
                      + N2*beta*z**0.6*absolute(psi2)**2)**1.5*absolute(psi2)**2)) \
          /trapz(r**2*absolute(psi2)**2)
    return mu1, mu2

def mu_an_uneqm_dens_ulck(psi1, psi2, r, V1, V2, dr,
                         gam, z, alpha, beta, eta,
                         lhy_ham1, lhy_ham2):

    dpsi1_dr = gradient(psi1, dr)
    mu1 = trapz(r**2*(0.5*gam*absolute(dpsi1_dr)**2 \
                      + V1*absolute(psi1)**2 \
                      + absolute(psi1)**4 \
                      + eta*absolute(psi1)**2*absolute(psi2)**2 \
                      + lhy_ham1*conj(psi1))) \
          /trapz(r**2*absolute(psi1)**2)
    
    dpsi2_dr = gradient(psi2,dr)
    mu2 = trapz(r**2*(0.5*(gam/z)*absolute(dpsi2_dr)**2 \
                      + V2*absolute(psi2)**2 \
                      + beta*absolute(psi2)**4 \
                      + eta*beta*absolute(psi1)**2*absolute(psi2)**2 \
                      + lhy_ham2*conj(psi2))) \
          /trapz(r**2*absolute(psi2)**2)
    return mu1, mu2
