import numpy as np
from scipy.sparse import eye

def ham_dens_lck(phi,V,r,dr,N_lck,mu,IM_REAL):
    """
    The Gross-Pitaevskii (GP) Hamiltonian for the density-locked mixture is defined here. This is to be called within 
    either imaginary time (IM_REAL = 0) or real time(IM_REAL = 1). The Hamiltonian is split into multiple contributions
    from:
    (1) H_ke - Kinetic energy
    (2) H_pot - Potential energy due to trapping potential
    (3) H_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    (4) H_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    (5) H_mu - Chemical potential term (calculated in a separate function)
    
       This requires the inputs:
    -> phi - the wavefunction
    -> V - the trapping potential
    -> r - the spatial array
    -> dr- the spatial discretisation
    -> N_lck - the effective atom number used to normalise wavefunction
    -> mu - the chemical potential, calculated in a separate function, but needed for the H_mu contribution here
    -> IM_REAL - this a 0 or 1 switch used to indicate imaginary or real time as this dictates whether the arrays
                 are designated as real (in imaginary time) or complex (in real time)
    """
    if IM_REAL == 0:
        # INITIALISE ARRAYS
        H_ke = np.zeros(phi.size)
        H_lhy = np.zeros(phi.size)
        H_int = np.zeros(phi.size)
        H_trap = np.zeros(phi.size)
        H_mu = np.zeros(phi.size)
        KE = np.zeros(phi.size)
    elif IM_REAL == 1:
        H_ke = np.zeros(phi.size).astype(complex)
        H_lhy = np.zeros(phi.size).astype(complex)
        H_int = np.zeros(phi.size).astype(complex)
        H_trap = np.zeros(phi.size).astype(complex)
        H_mu = np.zeros(phi.size).astype(complex)
        KE = np.zeros(phi.size).astype(complex)

    # first order derivative in the form of a sparse matrix (centrally defined)
    Dr = (1/(2*dr))*(-1*eye(phi.size-2,phi.size,k=0,dtype=float) + eye(phi.size-2,phi.size,k=2,dtype=float))
    
    # second order derivative in the form of a sparse matrix (centrally defined 3-point formula)
    Dr2 =  (1/dr**2)*(eye(phi.size-2,phi.size,k=0,dtype=float) -2*eye(phi.size-2,phi.size,k=1,dtype=float) + eye(phi.size-2,phi.size,k=2,dtype=float))

    KE[1:-1] = (2/r[1:-1])*(Dr @ phi) + Dr2 @ phi # Kinetic Energy derivatives
    # HAMILTONIAN TERMS
    H_ke[1:-1] = -0.5*KE[1:-1] # KE term
    H_lhy[1:-1] = 2.5*(N_lck**1.5)*np.abs(phi[1:-1])**3*phi[1:-1] # LHY term
    H_int[1:-1] = -3*N_lck*np.abs(phi[1:-1])**2*phi[1:-1] # s-wave term
    H_trap[1:-1] = V[1:-1]*phi[1:-1] # potential term
    H_mu[1:-1] = -mu*phi[1:-1]
    
    return H_ke,H_trap,H_int,H_lhy,H_mu

def eqm_dens_ulck_ham1(psi1,psi2,V1,r,dr,N1,N2,alpha,beta,eta,mu,im_real):
    """
    The Gross-Pitaevskii (GP) Hamiltonian for the 1st-component of the density-unlocked mixture is defined here. 
    This is to be called within  either imaginary time (IM_REAL = 0) or real time(IM_REAL = 1). The Hamiltonian 
    is split into multiple contributions
    from:
    (1) H_ke - Kinetic energy
    (2) H_pot - Potential energy due to trapping potential
    (3) H_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    (4) H_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    (5) H_mu - Chemical potential term (calculated in a separate function)
    
       This requires the inputs:
    -> psi1,psi2 - the wavefunctions
    -> V1 - the trapping potential of component 1
    -> r - the spatial array
    -> dr- the spatial discretisation
    -> N1,N2 - the atom numbers used to normalise the wavefunctions
    -> alpha,beta,eta - the dimensionless parameters used to define the density-unlocked, equal masses Hamiltonian
    -> mu - the chemical potential, calculated in a separate function, but needed for the H_mu contribution here
    -> IM_REAL - this a 0 or 1 switch used to indicate imaginary or real time as this dictates whether the arrays
                 are designated as real (in imaginary time) or complex (in real time)
    """
    if im_real == 0:
        # INITIALISE ARRAYS
        H_ke = np.zeros(psi1.size)
        H_lhy = np.zeros(psi1.size)
        H_int = np.zeros(psi1.size)
        H_trap = np.zeros(psi1.size)
        H_mu = np.zeros(psi1.size)
        KE = np.zeros(psi1.size)
    elif im_real == 1:
        H_ke = np.zeros(psi1.size).astype(complex)
        H_lhy = np.zeros(psi1.size).astype(complex)
        H_int = np.zeros(psi1.size).astype(complex)
        H_trap = np.zeros(psi1.size).astype(complex)
        H_mu = np.zeros(psi1.size).astype(complex)
        KE = np.zeros(psi1.size).astype(complex)

    # DIFFERENTIAL OPERATORS
    Dr = (1/(2*dr))*(-1*eye(psi1.size-2,psi1.size,k=0,dtype=float) + eye(psi1.size-2,psi1.size,k=2,dtype=float))
    Dr2 =  (1/dr**2)*(eye(psi1.size-2,psi1.size,k=0,dtype=float) -2*eye(psi1.size-2,psi1.size,k=1,dtype=float) + eye(psi1.size-2,psi1.size,k=2,dtype=float))

    KE[1:-1] = (2/r[1:-1])*(Dr @ psi1) + Dr2 @ psi1
    H_ke[1:-1] = -0.5*KE[1:-1] # KE term  
    
    H_lhy[1:-1] = alpha*(N1*np.abs(psi1[1:-1])**2 + beta*N2*np.abs(psi2[1:-1])**2)**1.5*psi1[1:-1]

    H_int[1:-1] = (N1*np.abs(psi1[1:-1])**2 + eta*N2*np.abs(psi2[1:-1])**2)*psi1[1:-1]

    H_trap[1:-1] = V1[1:-1]*psi1[1:-1]
    
    H_mu[1:-1] = -mu*psi1[1:-1]

    return H_ke,H_trap,H_int,H_lhy,H_mu
    
def eqm_dens_ulck_ham2(psi1,psi2,V2,r,dr,N1,N2,alpha,beta,eta,mu,im_real):
    """
    The Gross-Pitaevskii (GP) Hamiltonian for the 2nd-component of the density-unlocked mixture is defined here. 
    This is to be called within  either imaginary time (IM_REAL = 0) or real time(IM_REAL = 1). The Hamiltonian 
    is split into multiple contributions
    from:
    (1) H_ke - Kinetic energy
    (2) H_pot - Potential energy due to trapping potential
    (3) H_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    (4) H_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    (5) H_mu - Chemical potential term (calculated in a separate function)
    
       This requires the inputs:
    -> psi1,psi2 - the wavefunctions
    -> V2 - the trapping potential of component 2
    -> r - the spatial array
    -> dr- the spatial discretisation
    -> N1,N2 - the atom numbers used to normalise the wavefunctions
    -> alpha,beta,eta - the dimensionless parameters used to define the density-unlocked, equal masses Hamiltonian
    -> mu - the chemical potential, calculated in a separate function, but needed for the H_mu contribution here
    -> IM_REAL - this a 0 or 1 switch used to indicate imaginary or real time as this dictates whether the arrays
                 are designated as real (in imaginary time) or complex (in real time)
    """
    if im_real == 0:
        # INITIALISE ARRAYS
        H_ke = np.zeros(psi2.size)
        H_lhy = np.zeros(psi2.size)
        H_int = np.zeros(psi2.size)
        H_trap = np.zeros(psi2.size)
        H_mu = np.zeros(psi2.size)
        KE = np.zeros(psi2.size)
    elif im_real == 1:
        # INITIALISE ARRAYS
        H_ke = np.zeros(psi2.size).astype(complex)
        H_lhy = np.zeros(psi2.size).astype(complex)
        H_int = np.zeros(psi2.size).astype(complex)
        H_trap = np.zeros(psi2.size).astype(complex)
        H_mu = np.zeros(psi2.size).astype(complex)
        KE = np.zeros(psi2.size).astype(complex)

    # DIFFERENTIAL OPERATORS
    Dr = (1/(2*dr))*(-1*eye(psi2.size-2,psi2.size,k=0,dtype=float) + eye(psi2.size-2,psi2.size,k=2,dtype=float))
    Dr2 =  (1/dr**2)*(eye(psi2.size-2,psi2.size,k=0,dtype=float) -2*eye(psi2.size-2,psi2.size,k=1,dtype=float) + eye(psi2.size-2,psi2.size,k=2,dtype=float))

    KE[1:-1] = (2/r[1:-1])*(Dr @ psi2) + Dr2 @ psi2
    H_ke[1:-1] = -0.5*KE[1:-1] # KE term  
    
    H_lhy[1:-1] = alpha*beta**2*(N1*np.abs(psi1[1:-1])**2 + beta*N2*np.abs(psi2[1:-1])**2)**1.5*psi2[1:-1]

    H_int[1:-1] = (beta*N2*np.abs(psi2[1:-1])**2 + eta*beta*N1*np.abs(psi1[1:-1])**2)*psi2[1:-1]

    H_trap[1:-1] = V2[1:-1]*psi2[1:-1]

    H_mu[1:-1] = -mu*psi2[1:-1]
    
    return H_ke,H_trap,H_int,H_lhy,H_mu

def ham1_uneqm_dens_ulck(psi1,psi2,V1,r,dr,N1,N2,gam1,gam2,alpha,beta,eta,mu,im_real,z):
    """
    The Gross-Pitaevskii (GP) Hamiltonian for the 1st-component of the density-unlocked mixture is defined here. 
    This is to be called within  either imaginary time (IM_REAL = 0) or real time(IM_REAL = 1). The Hamiltonian 
    is split into multiple contributions
    from:
    (1) H_ke - Kinetic energy
    (2) H_pot - Potential energy due to trapping potential
    (3) H_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    (4) H_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    (5) H_mu - Chemical potential term (calculated in a separate function)
    
       This requires the inputs:
    -> psi1,psi2 - the wavefunctions
    -> V1 - the trapping potential of component 1
    -> r - the spatial array
    -> dr- the spatial discretisation
    -> N1,N2 - the atom numbers used to normalise the wavefunctions
    -> alpha,beta,eta - the dimensionless parameters used to define the density-unlocked, equal masses Hamiltonian
    -> mu - the chemical potential, calculated in a separate function, but needed for the H_mu contribution here
    -> IM_REAL - this a 0 or 1 switch used to indicate imaginary or real time as this dictates whether the arrays
                 are designated as real (in imaginary time) or complex (in real time)
    """
    if im_real == 0:
        # INITIALISE ARRAYS
        H_ke = np.zeros(psi1.size)
        H_lhy = np.zeros(psi1.size)
        H_int = np.zeros(psi1.size)
        H_trap = np.zeros(psi1.size)
        H_mu = np.zeros(psi1.size)
        KE = np.zeros(psi1.size)
    elif im_real == 1:
        H_ke = np.zeros(psi1.size).astype(complex)
        H_lhy = np.zeros(psi1.size).astype(complex)
        H_int = np.zeros(psi1.size).astype(complex)
        H_trap = np.zeros(psi1.size).astype(complex)
        H_mu = np.zeros(psi1.size).astype(complex)
        KE = np.zeros(psi1.size).astype(complex)

    # DIFFERENTIAL OPERATORS
    Dr = (1/(2*dr))*(-1*eye(psi1.size-2,psi1.size,k=0,dtype=float) + eye(psi1.size-2,psi1.size,k=2,dtype=float))
    Dr2 =  (1/dr**2)*(eye(psi1.size-2,psi1.size,k=0,dtype=float) -2*eye(psi1.size-2,psi1.size,k=1,dtype=float) + eye(psi1.size-2,psi1.size,k=2,dtype=float))

    KE[1:-1] = (2/r[1:-1])*(Dr @ psi1) + Dr2 @ psi1
    H_ke[1:-1] = -0.5*gam1*KE[1:-1] # KE term  

    H_lhy[1:-1] = N1**1.5*2.5*alpha*(1 + z**0.6*((N2/N1)*np.abs(psi2[1:-1])**2/np.abs(psi2[1:-1])**2))**1.5*np.abs(psi1[1:-1])**3*psi1[1:-1]

    H_int[1:-1] = (N1*np.abs(psi1[1:-1])**2 + eta*N2*np.abs(psi2[1:-1])**2)*psi1[1:-1]

    H_trap[1:-1] = V1[1:-1]*psi1[1:-1]
    
    H_mu[1:-1] = -mu*psi1[1:-1]

    return H_ke,H_trap,H_int,H_lhy,H_mu

def ham2_uneqm_dens_ulck(psi1,psi2,V2,r,dr,N1,N2,gam1,gam2,alpha,beta,eta,mu,im_real,z):
    """
    The Gross-Pitaevskii (GP) Hamiltonian for the 2nd-component of the density-unlocked mixture is defined here. 
    This is to be called within  either imaginary time (IM_REAL = 0) or real time(IM_REAL = 1). The Hamiltonian 
    is split into multiple contributions
    from:
    (1) H_ke - Kinetic energy
    (2) H_pot - Potential energy due to trapping potential
    (3) H_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    (4) H_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    (5) H_mu - Chemical potential term (calculated in a separate function)
    
       This requires the inputs:
    -> psi1,psi2 - the wavefunctions
    -> V2 - the trapping potential of component 2
    -> r - the spatial array
    -> dr- the spatial discretisation
    -> N1,N2 - the atom numbers used to normalise the wavefunctions
    -> alpha,beta,eta - the dimensionless parameters used to define the density-unlocked, equal masses Hamiltonian
    -> mu - the chemical potential, calculated in a separate function, but needed for the H_mu contribution here
    -> IM_REAL - this a 0 or 1 switch used to indicate imaginary or real time as this dictates whether the arrays
                 are designated as real (in imaginary time) or complex (in real time)
    """
    if im_real == 0:
        # INITIALISE ARRAYS
        H_ke = np.zeros(psi2.size)
        H_lhy = np.zeros(psi2.size)
        H_int = np.zeros(psi2.size)
        H_trap = np.zeros(psi2.size)
        H_mu = np.zeros(psi2.size)
        KE = np.zeros(psi2.size)
    elif im_real == 1:
        # INITIALISE ARRAYS
        H_ke = np.zeros(psi2.size).astype(complex)
        H_lhy = np.zeros(psi2.size).astype(complex)
        H_int = np.zeros(psi2.size).astype(complex)
        H_trap = np.zeros(psi2.size).astype(complex)
        H_mu = np.zeros(psi2.size).astype(complex)
        KE = np.zeros(psi2.size).astype(complex)

    # DIFFERENTIAL OPERATORS
    Dr = (1/(2*dr))*(-1*eye(psi2.size-2,psi2.size,k=0,dtype=float) + eye(psi2.size-2,psi2.size,k=2,dtype=float))
    Dr2 =  (1/dr**2)*(eye(psi2.size-2,psi2.size,k=0,dtype=float) -2*eye(psi2.size-2,psi2.size,k=1,dtype=float) + eye(psi2.size-2,psi2.size,k=2,dtype=float))

    KE[1:-1] = (2/r[1:-1])*(Dr @ psi2) + Dr2 @ psi2
    H_ke[1:-1] = -0.5*gam2*KE[1:-1] # KE term  
    
    H_lhy[1:-1] = 2.5*N1**1.5*alpha*beta*z**0.6*(1 + z**0.6*(N2*np.abs(psi2[1:-1])**2)/(N1*np.abs(psi1[1:-1])**2))**1.5*np.abs(psi1[1:-1])**3*psi2[1:-1]

    H_int[1:-1] = (beta*N2*np.abs(psi2[1:-1])**2 + eta*beta*N1*np.abs(psi1[1:-1])**2)*psi2[1:-1]

    H_trap[1:-1] = V2[1:-1]*psi2[1:-1]

    H_mu[1:-1] = -mu*psi2[1:-1]
    
    return H_ke,H_trap,H_int,H_lhy,H_mu
