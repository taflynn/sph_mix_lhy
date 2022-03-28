import numpy as np

def mu_dens_lck(phi,r,V,dr,N_lck):
    """
    The chemical potential, mu, of the density-locked mixture is included in the Gross-Pitaevskii (GP) Hamiltonian.
    This function calculates mu via an integration which is defined as the integral of the \phi^* H \phi. 
    
    This requires the inputs:
    -> phi - the wavefunction
    -> V - the trapping potential
    -> r - the spatial array
    -> dr- the spatial discretisation
    -> N_lck - the effective atom number used to normalise wavefunction
    
    The only output is the single chemical potential, mu.
    """
    dphi_dr = np.gradient(phi,dr)
    mu = np.trapz(r**2*(0.5*np.abs(dphi_dr)**2 \
                        + V*np.abs(phi)**2 \
                        - 3*N_lck*np.abs(phi)**4 \
                        + 2.5*N_lck**1.5*np.abs(phi)**5)) \
                        /np.trapz(r**2*np.abs(phi)**2)
    return mu

def mu_eqm_dens_ulck(psi1,psi2,r,V1,V2,dr,alpha,beta,eta,N1,N2):
    """
    The chemical potentials, mu1 and mu2, of the density-unlocked mixture are included in the Gross-Pitaevskii (GP) 
    Hamiltonians. This function calculates both mu1 and mu2 via an integration which is defined as the integral of the 
    \psi1^* H1 \psi1 and \psi2^* H2 \psi2. 
    
    This requires the inputs:
    -> psi1,psi2 - the wavefunctions
    -> V1,V2 - the trapping potentials
    -> r - the spatial array
    -> dr- the spatial discretisation
    -> alpha,beta,eta - the dimensionless parameters defining the Gross-Pitaevskii equations
    -> N1,N2 - the rescaled atom numbers used to normalise wavefunctions
    
    The outputs are the chemical potentials, mu1 and mu2.
    """
    dpsi1_dr = np.gradient(psi1,dr)
    mu1 = np.trapz(r**2*(0.5*np.abs(dpsi1_dr)**2 \
                         + V1*np.abs(psi1)**2 \
                         + N1*np.abs(psi1)**4 \
                         + N2*eta*np.abs(psi1)**2*np.abs(psi2)**2 \
                         + alpha*(N1*np.abs(psi1)**2 
                         + N2*beta*np.abs(psi2)**2)**1.5*np.abs(psi1)**2)) \
                         /np.trapz(r**2*np.abs(psi1)**2)
    
    dpsi2_dr = np.gradient(psi2,dr)
    mu2 = np.trapz(r**2*(0.5*np.abs(dpsi2_dr)**2 \
                         + V2*np.abs(psi2)**2 \
                         + N2*beta*np.abs(psi2)**4 \
                         + N1*eta*beta*np.abs(psi1)**2*np.abs(psi2)**2 \
                         + alpha*beta**2*(N1*np.abs(psi1)**2 
                         + N2*beta*np.abs(psi2)**2)**1.5*np.abs(psi2)**2)) \
                         /np.trapz(r**2*np.abs(psi2)**2)
    return mu1,mu2

def mu_uneqm_dens_ulck(psi1,psi2,r,V1,V2,dr,gam1,gam2,alpha,beta,eta,N1,N2,z):
    """
    The chemical potentials, mu1 and mu2, of the density-unlocked mixture are included in the Gross-Pitaevskii (GP) 
    Hamiltonians. This function calculates both mu1 and mu2 via an integration which is defined as the integral of the 
    \psi1^* H1 \psi1 and \psi2^* H2 \psi2. 
    
    This requires the inputs:
    -> psi1,psi2 - the wavefunctions
    -> V1,V2 - the trapping potentials
    -> r - the spatial array
    -> dr- the spatial discretisation
    -> alpha,beta,eta - the dimensionless parameters defining the Gross-Pitaevskii equations
    -> N1,N2 - the rescaled atom numbers used to normalise wavefunctions
    
    The outputs are the chemical potentials, mu1 and mu2.
    """
    dpsi1_dr = np.gradient(psi1,dr)
    mu1 = np.trapz(r**2*(0.5*gam1*np.abs(dpsi1_dr)**2 \
                         + V1*np.abs(psi1)**2 \
                         + N1*np.abs(psi1)**4 \
                         + N2*eta*np.abs(psi1)**2*np.abs(psi2)**2 \
                         + 1.25*alpha*(1 + z**0.6*((N2*np.abs(psi2)**2)/(N1*np.abs(psi1)**2)))**1.5 \
                         *(5*N1**1.5*np.abs(psi1)**3 + 3*N1**0.5*N2*z**0.6*np.abs(psi1)*np.abs(psi2)**2)*np.abs(psi1)**2)) \
                         /np.trapz(r**2*np.abs(psi1)**2)
    
    dpsi2_dr = np.gradient(psi2,dr)
    mu2 = np.trapz(r**2*(0.5*gam2*np.abs(dpsi2_dr)**2 \
                         + V2*np.abs(psi2)**2 \
                         + N2*beta*np.abs(psi2)**4 \
                         + N1*eta*beta*np.abs(psi1)**2*np.abs(psi2)**2 \
                         + 2.5*N1**1.5*alpha*beta**(1 + z**0.5*((N2*np.abs(psi2)**2)/(N1*np.abs(psi1)**2)))**1.5*np.abs(psi1)**3*np.abs(psi2)**2)) \
                         /np.trapz(r**2*np.abs(psi2)**2)
    return mu1,mu2
