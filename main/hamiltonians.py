from numpy import zeros, absolute, vectorize
# from scipy.sparse import eye

def ham_dens_lck(phi, V, r, Dr, Dr2, dr, N_lck, mu, IM_REAL):
    """
    The Gross-Pitaevskii (GP) Hamiltonian for the density-locked mixture is defined here. This is to be called within 
    either imaginary time (IM_REAL = 0) or real time(IM_REAL = 1). The Hamiltonian is split into multiple contributions
    from:
    * H_ke - Kinetic energy
    * H_pot - Potential energy due to trapping potential
    * H_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    * H_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    * H_mu - Chemical potential term (calculated in a separate function)
    """
    if IM_REAL == 0:
        # INITIALISE ARRAYS
        H_ke = zeros(phi.size)
        H_lhy = zeros(phi.size)
        H_int = zeros(phi.size)
        H_trap = zeros(phi.size)
        H_mu = zeros(phi.size)
        KE = zeros(phi.size)
    elif IM_REAL == 1:
        H_ke = zeros(phi.size).astype(complex)
        H_lhy = zeros(phi.size).astype(complex)
        H_int = zeros(phi.size).astype(complex)
        H_trap = zeros(phi.size).astype(complex)
        H_mu = zeros(phi.size).astype(complex)
        KE = zeros(phi.size).astype(complex)

    KE[1:-1] = (2/r[1:-1])*(Dr @ phi) + (Dr2 @ phi) # Kinetic Energy derivatives
    # HAMILTONIAN TERMS
    H_ke[1:-1] = -0.5*KE[1:-1] # KE term
    H_lhy[1:-1] = 2.5*(N_lck**1.5)*absolute(phi[1:-1])**3*phi[1:-1] # LHY term
    H_int[1:-1] = -3*N_lck*absolute(phi[1:-1])**2*phi[1:-1] # s-wave term
    H_trap[1:-1] = V[1:-1]*phi[1:-1] # potential term
    H_mu[1:-1] = -mu*phi[1:-1]
    
    return H_ke, H_trap, H_int, H_lhy, H_mu

def eqm_dens_ulck_ham1(psi1, psi2, V1, r, Dr, Dr2, dr, N1, N2, alpha, beta, eta, mu, im_real):
    """
    The Gross-Pitaevskii (GP) Hamiltonian for the 1st-component of the density-unlocked mixture is defined here. 
    This is to be called within  either imaginary time (IM_REAL = 0) or real time(IM_REAL = 1). The Hamiltonian 
    is split into multiple contributions
    from:
    * H_ke - Kinetic energy
    * H_pot - Potential energy due to trapping potential
    * H_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    * H_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    * H_mu - Chemical potential term (calculated in a separate function)
    """
    if im_real == 0:
        # INITIALISE ARRAYS
        H_ke = zeros(psi1.size)
        H_lhy = zeros(psi1.size)
        H_int = zeros(psi1.size)
        H_trap = zeros(psi1.size)
        H_mu = zeros(psi1.size)
        KE = zeros(psi1.size)
    elif im_real == 1:
        H_ke = zeros(psi1.size).astype(complex)
        H_lhy = zeros(psi1.size).astype(complex)
        H_int = zeros(psi1.size).astype(complex)
        H_trap = zeros(psi1.size).astype(complex)
        H_mu = zeros(psi1.size).astype(complex)
        KE = zeros(psi1.size).astype(complex)

    KE[1:-1] = (2/r[1:-1])*(Dr @ psi1) + (Dr2 @ psi1)
    H_ke[1:-1] = -0.5*KE[1:-1] # KE term  
    
    H_lhy[1:-1] = alpha*(N1*absolute(psi1[1:-1])**2 + beta*N2*absolute(psi2[1:-1])**2)**1.5*psi1[1:-1]

    H_int[1:-1] = (N1*absolute(psi1[1:-1])**2 + eta*N2*absolute(psi2[1:-1])**2)*psi1[1:-1]

    H_trap[1:-1] = V1[1:-1]*psi1[1:-1]
    
    H_mu[1:-1] = -mu*psi1[1:-1]

    return H_ke, H_trap, H_int, H_lhy, H_mu

def eqm_dens_ulck_ham2(psi1, psi2, V2, r, Dr, Dr2, dr, N1, N2, alpha, beta, eta, mu, im_real):
    """
    The Gross-Pitaevskii (GP) Hamiltonian for the 2nd-component of the density-unlocked mixture is defined here. 
    This is to be called within  either imaginary time (IM_REAL = 0) or real time(IM_REAL = 1). The Hamiltonian 
    is split into multiple contributions
    from:
    * H_ke - Kinetic energy
    * H_pot - Potential energy due to trapping potential
    * H_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    * H_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    * H_mu - Chemical potential term (calculated in a separate function)
    """
    if im_real == 0:
        # INITIALISE ARRAYS
        H_ke = zeros(psi2.size)
        H_lhy = zeros(psi2.size)
        H_int = zeros(psi2.size)
        H_trap = zeros(psi2.size)
        H_mu = zeros(psi2.size)
        KE = zeros(psi2.size)
    elif im_real == 1:
        # INITIALISE ARRAYS
        H_ke = zeros(psi2.size).astype(complex)
        H_lhy = zeros(psi2.size).astype(complex)
        H_int = zeros(psi2.size).astype(complex)
        H_trap = zeros(psi2.size).astype(complex)
        H_mu = zeros(psi2.size).astype(complex)
        KE = zeros(psi2.size).astype(complex)

    KE[1:-1] = (2/r[1:-1])*(Dr @ psi2) + (Dr2 @ psi2)
    H_ke[1:-1] = -0.5*KE[1:-1] # KE term  
    
    H_lhy[1:-1] = alpha*beta**2*(N1*absolute(psi1[1:-1])**2 + beta*N2*absolute(psi2[1:-1])**2)**1.5*psi2[1:-1]

    H_int[1:-1] = (beta*N2*absolute(psi2[1:-1])**2 + eta*beta*N1*absolute(psi1[1:-1])**2)*psi2[1:-1]

    H_trap[1:-1] = V2[1:-1]*psi2[1:-1]

    H_mu[1:-1] = -mu*psi2[1:-1]
    
    return H_ke, H_trap, H_int, H_lhy, H_mu

def ham1_uneqm_dens_ulck(psi1, psi2, V1, r, Dr, Dr2, dr, N1, N2, gam, z, alpha, beta, eta, mu, im_real):
    """
    The Gross-Pitaevskii (GP) Hamiltonian for the 1st-component of the density-unlocked mixture is defined here. 
    This is to be called within  either imaginary time (IM_REAL = 0) or real time(IM_REAL = 1). The Hamiltonian 
    is split into multiple contributions
    from:
    * H_ke - Kinetic energy
    * H_pot - Potential energy due to trapping potential
    * H_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    * H_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    * H_mu - Chemical potential term (calculated in a separate function)
    """
    if im_real == 0:
        # INITIALISE ARRAYS
        H_ke = zeros(psi1.size)
        H_lhy = zeros(psi1.size)
        H_int = zeros(psi1.size)
        H_trap = zeros(psi1.size)
        H_mu = zeros(psi1.size)
        H_3bl = zeros(psi1.size)
        KE = zeros(psi1.size)
    elif im_real == 1:
        H_ke = zeros(psi1.size).astype(complex)
        H_lhy = zeros(psi1.size).astype(complex)
        H_int = zeros(psi1.size).astype(complex)
        H_trap = zeros(psi1.size).astype(complex)
        H_mu = zeros(psi1.size).astype(complex)
        H_3bl = zeros(psi1.size).astype(complex)
        KE = zeros(psi1.size).astype(complex)

    KE[1:-1] = (2/r[1:-1])*(Dr @ psi1) + (Dr2 @ psi1)
    H_ke[1:-1] = -0.5*gam*KE[1:-1] # KE term  

    H_lhy[1:-1] = alpha*(N1*absolute(psi1[1:-1])**2 + N2*z**0.6*beta*absolute(psi2[1:-1])**2)**1.5*psi1[1:-1]

    H_int[1:-1] = (N1*absolute(psi1[1:-1])**2 + eta*N2*absolute(psi2[1:-1])**2)*psi1[1:-1]

    H_trap[1:-1] = V1[1:-1]*psi1[1:-1]
    
    H_mu[1:-1] = -mu*psi1[1:-1]

    return H_ke, H_trap, H_int, H_lhy, H_mu

def ham2_uneqm_dens_ulck(psi1, psi2, V2, r, Dr, Dr2, dr, N1, N2, gam, z, alpha, beta, eta, mu, im_real):
    """
    The Gross-Pitaevskii (GP) Hamiltonian for the 2nd-component of the density-unlocked mixture is defined here. 
    This is to be called within  either imaginary time (IM_REAL = 0) or real time(IM_REAL = 1). The Hamiltonian 
    is split into multiple contributions
    from:
    * H_ke - Kinetic energy
    * H_pot - Potential energy due to trapping potential
    * H_int - Interaction energy due to effective interaction term (resulting from two-body interactions)
    * H_lhy - Lee-Huang-Yang energy due to quantum fluctuations
    * H_mu - Chemical potential term (calculated in a separate function)
    """
    if im_real == 0:
        # INITIALISE ARRAYS
        H_ke = zeros(psi2.size)
        H_lhy = zeros(psi2.size)
        H_int = zeros(psi2.size)
        H_trap = zeros(psi2.size)
        H_mu = zeros(psi2.size)
        KE = zeros(psi2.size)
    elif im_real == 1:
        # INITIALISE ARRAYS
        H_ke = zeros(psi2.size).astype(complex)
        H_lhy = zeros(psi2.size).astype(complex)
        H_int = zeros(psi2.size).astype(complex)
        H_trap = zeros(psi2.size).astype(complex)
        H_mu = zeros(psi2.size).astype(complex)
        KE = zeros(psi2.size).astype(complex)

    KE[1:-1] = (2/r[1:-1])*(Dr @ psi2) + (Dr2 @ psi2)
    H_ke[1:-1] = -0.5*(gam/z)*KE[1:-1] # KE term  
    
    H_lhy[1:-1] = alpha*beta**2.0*z**0.6*(N1*absolute(psi1[1:-1])**2 + N2*z**0.6*beta*absolute(psi2[1:-1])**2)**1.5*psi2[1:-1]

    H_int[1:-1] = (beta*N2*absolute(psi2[1:-1])**2 + eta*beta*N1*absolute(psi1[1:-1])**2)*psi2[1:-1]

    H_trap[1:-1] = V2[1:-1]*psi2[1:-1]

    H_mu[1:-1] = -mu*psi2[1:-1]
    
    return H_ke, H_trap, H_int, H_lhy, H_mu

def ham1_an_uneqm_dens_ulck(psi1, psi2, V1, r, Dr, Dr2, dr,
                           gam, z, alpha, beta, eta, 
                           mu, im_real,
                           lhy_ham1):
    
    if im_real == 0:
        # INITIALISE ARRAYS
        H_ke = zeros(psi1.size)
        H_lhy = zeros(psi1.size)
        H_int = zeros(psi1.size)
        H_trap = zeros(psi1.size)
        H_mu = zeros(psi1.size)
        KE = zeros(psi1.size)
    elif im_real == 1:
        H_ke = zeros(psi1.size).astype(complex)
        H_lhy = zeros(psi1.size).astype(complex)
        H_int = zeros(psi1.size).astype(complex)
        H_trap = zeros(psi1.size).astype(complex)
        H_mu = zeros(psi1.size).astype(complex)
        KE = zeros(psi1.size).astype(complex)
        
    KE[1:-1] = (2/r[1:-1])*(Dr @ psi1) + (Dr2 @ psi1)
    H_ke[1:-1] = -0.5*gam*KE[1:-1] # KE term  

    H_lhy[1:-1] = lhy_ham1[1:-1]

    H_int[1:-1] = (absolute(psi1[1:-1])**2 
                   + eta*absolute(psi2[1:-1])**2)*psi1[1:-1]

    H_trap[1:-1] = V1[1:-1]*psi1[1:-1]
    
    H_mu[1:-1] = -mu*psi1[1:-1]
    
    return H_ke, H_trap, H_int, H_lhy, H_mu

def ham2_an_uneqm_dens_ulck(psi1, psi2, V2, r, Dr, Dr2, dr,
                           gam, z, alpha, beta, eta, 
                           mu, im_real,
                           lhy_ham2):
    
    if im_real == 0:
        # INITIALISE ARRAYS
        H_ke = zeros(psi1.size)
        H_lhy = zeros(psi1.size)
        H_int = zeros(psi1.size)
        H_trap = zeros(psi1.size)
        H_mu = zeros(psi1.size)
        KE = zeros(psi1.size)
    elif im_real == 1:
        H_ke = zeros(psi1.size).astype(complex)
        H_lhy = zeros(psi1.size).astype(complex)
        H_int = zeros(psi1.size).astype(complex)
        H_trap = zeros(psi1.size).astype(complex)
        H_mu = zeros(psi1.size).astype(complex)
        KE = zeros(psi1.size).astype(complex)
        
    KE[1:-1] = (2/r[1:-1])*(Dr @ psi2) + Dr2 @ psi2
    H_ke[1:-1] = -0.5*(gam/z)*KE[1:-1] # KE term 
    
    H_lhy[1:-1] = lhy_ham2[1:-1]

    H_int[1:-1] = (beta*absolute(psi2[1:-1])**2 
                   + eta*beta*absolute(psi1[1:-1])**2)*psi2[1:-1]

    H_trap[1:-1] = V2[1:-1]*psi2[1:-1]

    H_mu[1:-1] = -mu*psi2[1:-1]

    return H_ke, H_trap, H_int, H_lhy, H_mu

def lhy_category(psi1, psi2, beta, tol, f_interp, df_interp):
    if absolute(psi1)**2 <= tol:
        if absolute(psi2)**2 <= tol:
            x = 0
            f_category = 'std'
            
        elif absolute(psi2)**2 > tol:
            x = (1/beta)*absolute(psi1)**2/absolute(psi2)**2
            f_category = 'inv'
            
    elif absolute(psi1)**2 > tol and absolute(psi2)**2 > tol and beta*(absolute(psi2)**2/absolute(psi1)**2) > 1:
        x = (1/beta)*absolute(psi1)**2/absolute(psi2)**2
        f_category = 'inv'
        
    else:
        x = beta*absolute(psi2)**2/absolute(psi1)**2
        f_category = 'std'
        
    f = f_interp(x)
    dfdx = df_interp(x)

    return x, f, dfdx, f_category


def lhy_ham1_uneqm(psi1, psi2, alpha, beta, f, dfdx):
    return alpha*(f*absolute(psi1)**3 
                  - 0.4*beta*dfdx*absolute(psi1)*absolute(psi2)**2)*psi1

def lhy_ham2_uneqm(psi1, psi2, alpha, beta, f, dfdx):
    return 0.4*alpha*beta**2*dfdx*absolute(psi1)**3*psi2

def inv_lhy_ham1_uneqm(psi1, psi2, alpha, beta, f, dfdx):
    return 0.4*alpha*beta**1.5*dfdx*absolute(psi2)**3*psi1

def inv_lhy_ham2_uneqm(psi1, psi2, alpha, beta, f, dfdx):
    return alpha*beta**2.5*(beta*f*absolute(psi2)**3 
                            - 0.4*dfdx*absolute(psi1)**2*absolute(psi2))*psi2

def ham_checker(f_category, psi1, psi2, alpha, beta, f, dfdx):   
    if f_category == 'std':
        combined_lhy_ham1 = lhy_ham1_uneqm(psi1, psi2, alpha, beta, f, dfdx)
        combined_lhy_ham2 = lhy_ham2_uneqm(psi1, psi2, alpha, beta, f, dfdx)
        
    elif f_category == 'inv':
        combined_lhy_ham1 = inv_lhy_ham1_uneqm(psi1, psi2, alpha, beta, f, dfdx)
        combined_lhy_ham2 = inv_lhy_ham2_uneqm(psi1, psi2, alpha, beta, f, dfdx)
        
    return combined_lhy_ham1, combined_lhy_ham2


def ham_lhy_constructor(psi1, psi2, alpha, beta, f_interp, df_interp):

    vect_ham_checker = vectorize(ham_checker)
    vect_lhy_category = vectorize(lhy_category)

    [x, f, dfdx, f_category] = vect_lhy_category(psi1, psi2, beta, 1e-10,
                                                 f_interp, df_interp)
    
    [combined_lhy_ham1, combined_lhy_ham2] = vect_ham_checker(f_category, psi1, psi2, 
                                                              alpha, beta, f, dfdx)
    
    return combined_lhy_ham1, combined_lhy_ham2
