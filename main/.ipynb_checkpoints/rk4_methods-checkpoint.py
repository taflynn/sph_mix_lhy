from main.hamiltonians import ham_dens_lck
from main.hamiltonians import eqm_dens_ulck_ham1,eqm_dens_ulck_ham2
from main.boundary import bc
from main.chem_pot import mu_dens_lck,mu_eqm_dens_ulck
from main.energy import energy_eqm_dens_ulck,energy_eqm_dens_lck
from numpy import pi

import numpy as np

import matplotlib.pyplot as plt

def rk4_dens_lck(r,phi,V,N_lck,dr,dt,T_STEPS,T_SAVE,IM_REAL,BC_TYPE):
    """
    The rk4* functions are the main body of this package. They contain the Runge-Kutta 4th-order time-stepping methods
    for solving the Gross-Pitaevskii (GP) equations. 
    
    The outputs of this function are:
    (1) phi - the final wavefunction after imaginary or real time propagation
    (2) mu - the associated chemical potential of the density-locked mixture
    (3) t_array - a 1D array of the times saved at every T_SAVE time-steps
    (4) spacetime - an 2D array of the wavefunction saved at every T_SAVE time-steps
    (5) E_array - a 1D array of the energy of the mixture saved at every T_SAVE time-steps 
    
    This function solves the density-locked GP equation, and takes the following inputs:
    -> r - the spatial array
    -> phi - the initial guess wavefunction
    -> V - the trapping potential
    -> N_lck - the effective atom number used to normalise wavefunction
    -> dr- the spatial discretisation
    -> dt - the time discretisation
    -> T_STEPS - the number of time-steps (either real or imaginary time-steps)
    -> T_SAVE - the number of time-steps between saving snapshots of the data
    -> IM_REAL - a switch that informs the function of whether imaginary (IM_REAL = 0) or real (IM_REAL) time is required
    -> BC_TYPE - a switch that informs the function of which boundary conditions to implement (see 'boundary.py' for
                 further details
    """
    if IM_REAL == 0:
        # initalise Runge-Kutta arrays
        k1 = np.zeros(phi.size)
        k2 = np.zeros(phi.size)
        k3 = np.zeros(phi.size)
        k4 = np.zeros(phi.size)
        
        # initialise array to save wavefunction snapshots
        spacetime = np.zeros((r.size,(T_STEPS//T_SAVE)))
        
    elif IM_REAL == 1:
        # for real time, specify a complex time step
        dt = 1j*dt
        phi = phi.astype(complex)

        # initialise Runge-Kutta arrays
        k1 = np.zeros(phi.size).astype(complex)
        k2 = np.zeros(phi.size).astype(complex)
        k3 = np.zeros(phi.size).astype(complex)
        k4 = np.zeros(phi.size).astype(complex)
        
        # initialise array to save wavefunction snapshots
        spacetime = np.zeros((r.size,(T_STEPS//T_SAVE))).astype(complex)
        
    # initialise data saving arrays
    t_array = np.zeros((T_STEPS//T_SAVE))
    E_array = np.zeros((T_STEPS//T_SAVE))
    
    # initialise variables used as counters or for convergence calculations
    t = 0.0
    previous_mode = 1.0
    mu = 0.0
    mu_prev = 1.0
    E_prev = 1.0
    
    for l in range(0,T_STEPS):
        # k1 calculation
        [H_ke,H_trap,H_int,H_lhy,H_mu] = ham_dens_lck(phi,V,r,dr,N_lck,mu,IM_REAL)

        k1[1:-1] = -dt*(H_ke[1:-1] + H_trap[1:-1] + H_lhy[1:-1] + H_int[1:-1] + H_mu[1:-1])

        # apply boundary conditions
        k1 = bc(k1,BC_TYPE)

        # k2 calculation
        [H_ke,H_trap,H_int,H_lhy,H_mu] = ham_dens_lck(phi + k1/2.0,V,r,dr,N_lck,mu,IM_REAL)
        
        k2[1:-1] = -dt*(H_ke[1:-1] + H_trap[1:-1] + H_lhy[1:-1] + H_int[1:-1] + H_mu[1:-1])

        # apply boundary conditions
        k2 = bc(k2,BC_TYPE)

        # k3 calculation
        [H_ke,H_trap,H_int,H_lhy,H_mu] = ham_dens_lck(phi + k2/2.0,V,r,dr,N_lck,mu,IM_REAL)

        k3[1:-1] = -dt*(H_ke[1:-1] + H_trap[1:-1] + H_lhy[1:-1] + H_int[1:-1] + H_mu[1:-1])

        # apply boundary conditions
        k3 = bc(k3,BC_TYPE)

        #k4 calculation
        [H_ke,H_trap,H_int,H_lhy,H_mu] = ham_dens_lck(phi + k3,V,r,dr,N_lck,mu,IM_REAL)
        
        k4[1:-1] = -dt*(H_ke[1:-1] + H_trap[1:-1] + H_lhy[1:-1] + H_int[1:-1] + H_mu[1:-1])

        # apply boundary conditions
        k4 = bc(k4,BC_TYPE)

        # final Runge-Kutta step    
        phi[1:-1] = phi[1:-1] + (1.0/6.0)*(k1[1:-1] + 2*k2[1:-1] + 2*k3[1:-1] + k4[1:-1])
        
        # apply boundary conditions
        phi = bc(phi,BC_TYPE)
        
        if (IM_REAL == 0): 
            # normalise wavefunction in imaginary time
            Norm = 4*pi*np.trapz(r**2*abs(phi)**2)*dr
            phi = phi/np.sqrt(Norm)
            
            # calculate chemical potential
            mu = mu_dens_lck(phi,r,V,dr,N_lck)
        
            # tolerance between successive chemical potentials
            mu_tol = np.abs(mu - mu_prev)/np.abs(mu_prev)
            
            # iterate chemical potential
            mu_prev = mu
        
        if (IM_REAL == 0 and l % T_SAVE == 0): 
            # convergence of max density
            current_mode = N_lck*np.abs(phi[1])**2
            tol_mode = np.abs((current_mode - previous_mode)/previous_mode)
            previous_mode = current_mode
            
            # convergence of energy
            E_ke,E_pot,E_int,E_lhy = energy_eqm_dens_lck(phi,V,r,dr,N_lck)
            E_total = E_ke + E_pot + E_int + E_lhy
            E_tol = np.abs(E_total - E_prev)/np.abs(E_prev)
            E_prev = E_total
            
            # printed output
            print('-'*21,'Max Densty Convergence','-'*21)
            print('l = ',l,'(Percentage of imaginary time done = ',100*(l/T_STEPS),'%)')
            print('Tolerance between successive max densities = ',tol_mode)
            print('-'*17,'Chemical Potential Convergence','-'*17)
            print('Chemical potential = ',mu)
            print('Tolerance of chemical potential = ',mu_tol)
            print('-'*23,'Energy Convergence','-'*23)            
            print('Total energy = ',E_total)
            print('Tolerance between successive total energies = ',E_tol)
            print('-'*66)
            
        elif (IM_REAL == 1):
            print('l = ',l,'(Percentage of real time done = ',100*(l/T_STEPS),'%)')
            
        # save data
        if (l % T_SAVE == 0):
            # save energy
            E_array[l//T_SAVE] = E_total
            
            # save wavefunction
            spacetime[:,l//T_SAVE] = np.sqrt(N_lck)*phi
            
            # save current time
            t_array[l//T_SAVE] = t
            
            # plotting density
            plt.plot(r,N_lck*np.abs(phi)**2)
            plt.xlim((0,np.max(r)))
            plt.ylim(0,1.2*N_lck*np.abs(phi[1])**2)
            plt.xlabel(r'$r$')
            plt.ylabel(r'$n_0(r)$')
            plt.close()    
        
        # iterate time
        t = t + dt
    
    if IM_REAL == 0:
        print('Imaginary time finished')
    if IM_REAL == 1:
        print('Real time finished')
    return phi,mu,t_array,spacetime,E_array

def rk4_dens_ulck(r,psi1,psi2,V1,V2,alpha,beta,eta,N1,N2,dr,dt,T_STEPS,T_SAVE,IM_REAL,BC_TYPE):
    """
    The rk4* functions are the main body of this package. They contain the Runge-Kutta 4th-order time-stepping methods
    for solving the Gross-Pitaevskii (GP) equations. 
    
    The outputs of this function are:
    (1) psi1,psi2 - the final wavefunctions after imaginary or real time propagation
    (2) mu1,mu2 - the associated chemical potentials of the density-unlocked mixture
    (3) t_array - a 1D array of the times saved at every T_SAVE time-steps
    (4) spacetime1,spacetime2 - 2D arrays of the wavefunctions saved at every T_SAVE time-steps
    (5) E1_array,E2_array - 1D arrays of the energy of the mixture saved at every T_SAVE time-steps 
    
    This function solves the density-locked GP equation, and takes the following inputs:
    -> r - the spatial array
    -> psi1,psi2 - the initial guess wavefunctions
    -> V1,V2 - the trapping potentials
    -> alpha,beta,eta - the dimensionless parameters defining the GP equations
    -> N1,N2 - the rescaled atom numbers used to normalise wavefunctions
    -> dr- the spatial discretisation
    -> dt - the time discretisation
    -> T_STEPS - the number of time-steps (either real or imaginary time-steps)
    -> T_SAVE - the number of time-steps between saving snapshots of the data
    -> IM_REAL - a switch that informs the function of whether imaginary (IM_REAL = 0) or real (IM_REAL) time is required
    -> BC_TYPE - a switch that informs the function of which boundary conditions to implement (see 'boundary.py' for
                 further details
    """
    if IM_REAL == 0:
        # initalise Runge-Kutta arrays
        k1_1 = np.zeros(psi1.size)
        k2_1 = np.zeros(psi1.size)
        k3_1 = np.zeros(psi1.size)
        k4_1 = np.zeros(psi1.size)
        k1_2 = np.zeros(psi2.size)
        k2_2 = np.zeros(psi2.size)
        k3_2 = np.zeros(psi2.size)
        k4_2 = np.zeros(psi2.size)
        
        # initialise array to save wavefunction snapshots
        spacetime1 = np.zeros((r.size,(T_STEPS//T_SAVE)))
        spacetime2 = np.zeros((r.size,(T_STEPS//T_SAVE)))
        
        t = 0.0
        
    elif IM_REAL == 1:
        # for real time, specify a complex time step
        dt = 1j*dt
        
        print('dt =',dt)
        
        psi1 = psi1.astype(complex)
        psi2 = psi2.astype(complex)
        
        k1_1 = np.zeros(psi1.size).astype(complex)
        k2_1 = np.zeros(psi1.size).astype(complex)
        k3_1 = np.zeros(psi1.size).astype(complex)
        k4_1 = np.zeros(psi1.size).astype(complex)
        k1_2 = np.zeros(psi2.size).astype(complex)
        k2_2 = np.zeros(psi2.size).astype(complex)
        k3_2 = np.zeros(psi2.size).astype(complex)
        k4_2 = np.zeros(psi2.size).astype(complex)
        
        # initialise array to save wavefunction snapshots
        spacetime1 = np.zeros((r.size,(T_STEPS//T_SAVE))).astype(complex)
        spacetime2 = np.zeros((r.size,(T_STEPS//T_SAVE))).astype(complex)
        
        t = 0.0 + 0.0j
        
    E1_array = np.zeros((T_STEPS//T_SAVE))  
    E2_array = np.zeros((T_STEPS//T_SAVE))    
    t_array = np.zeros((T_STEPS//T_SAVE))
    
    # initialise variables used as counters or for convergence calculations    
    previous_mode = 1.0
    mu1 = 0.0
    mu2 = 0.0
    mu1_prev = 1.0
    mu2_prev = 1.0
    mu1_tol = 1.0
    mu2_tol = 1.0
    E_prev = 1.0
    
    for l in range(0,T_STEPS):
        # k1 CALCULATION FOR PSI1
        [H_ke1,H_trap1,H_int1,H_lhy1,H_mu1] = eqm_dens_ulck_ham1(psi1,psi2,V1,r,dr,N1,N2,alpha,beta,eta,mu1,IM_REAL)

        k1_1[1:-1] = -dt*(H_ke1[1:-1] + H_trap1[1:-1] + H_lhy1[1:-1] + H_int1[1:-1] + H_mu1[1:-1])

        # BOUNDARY CONDITION
        k1_1 = bc(k1_1,BC_TYPE)

        # k1 CALCULATION FOR PSI2
        [H_ke2,H_trap2,H_int2,H_lhy2,H_mu2] = eqm_dens_ulck_ham2(psi1,psi2,V2,r,dr,N1,N2,alpha,beta,eta,mu2,IM_REAL)

        k1_2[1:-1] = -dt*(H_ke2[1:-1] + H_trap2[1:-1] + H_lhy2[1:-1] + H_int2[1:-1] + H_mu2[1:-1])

        # BOUNDARY CONDITION
        k1_2 = bc(k1_2,BC_TYPE)

        # k2 CALCULATION FOR PSI1
        [H_ke1,H_trap1,H_int1,H_lhy1,H_mu1] = eqm_dens_ulck_ham1(psi1+k1_1/2.0,psi2+k1_2/2.0,V1,r,dr,N1,N2,alpha,beta,eta,mu1,IM_REAL)

        k2_1[1:-1] = -dt*(H_ke1[1:-1] + H_trap1[1:-1] + H_lhy1[1:-1] + H_int1[1:-1] + H_mu1[1:-1])

        # BOUNDARY CONDITION
        k2_1 = bc(k2_1,BC_TYPE)

        # k2 CALCULATION FOR PSI2
        [H_ke2,H_trap2,H_int2,H_lhy2,H_mu2] = eqm_dens_ulck_ham2(psi1+k1_1/2.0,psi2+k1_2/2.0,V2,r,dr,N1,N2,alpha,beta,eta,mu2,IM_REAL)

        k2_2[1:-1] = -dt*(H_ke2[1:-1] + H_trap2[1:-1] + H_lhy2[1:-1] + H_int2[1:-1] + H_mu2[1:-1])

        # BOUNDARY CONDITION
        k2_2 = bc(k2_2,BC_TYPE)

        # k3 CALCULATION FOR PSI1
        [H_ke1,H_trap1,H_int1,H_lhy1,H_mu1] = eqm_dens_ulck_ham1(psi1+k2_1/2.0,psi2+k2_2/2.0,V1,r,dr,N1,N2,alpha,beta,eta,mu1,IM_REAL)

        k3_1[1:-1] = -dt*(H_ke1[1:-1] + H_trap1[1:-1] + H_lhy1[1:-1] + H_int1[1:-1] + H_mu1[1:-1])

        # BOUNDARY CONDITION
        k3_1 = bc(k3_1,BC_TYPE)

        # k3 CALCULATION FOR PSI2
        [H_ke2,H_trap2,H_int2,H_lhy2,H_mu2] = eqm_dens_ulck_ham2(psi1+k2_1/2.0,psi2+k2_2/2.0,V2,r,dr,N1,N2,alpha,beta,eta,mu2,IM_REAL)

        k3_2[1:-1] = -dt*(H_ke2[1:-1] + H_trap2[1:-1] + H_lhy2[1:-1] + H_int2[1:-1] + H_mu2[1:-1])

        # BOUNDARY CONDITION
        k3_2 = bc(k3_2,BC_TYPE)

        # k4 CALCULATION FOR PSI1
        [H_ke1,H_trap1,H_int1,H_lhy1,H_mu2] = eqm_dens_ulck_ham1(psi1+k3_1,psi2+k3_2,V1,r,dr,N1,N2,alpha,beta,eta,mu1,IM_REAL)

        k4_1[1:-1] = -dt*(H_ke1[1:-1] + H_trap1[1:-1] + H_lhy1[1:-1] + H_int1[1:-1] + H_mu1[1:-1])

        # BOUNDARY CONDITION
        k4_1 = bc(k4_1,BC_TYPE)

        # k4 CALCULATION FOR PSI2
        [H_ke2,H_trap2,H_int2,H_lhy2,H_mu2] = eqm_dens_ulck_ham2(psi1+k3_1,psi2+k3_2,V2,r,dr,N1,N2,alpha,beta,eta,mu2,IM_REAL)

        k4_2[1:-1] = -dt*(H_ke2[1:-1] + H_trap2[1:-1] + H_lhy2[1:-1] + H_int2[1:-1] + H_mu2[1:-1])

        # BOUNDARY CONDITION
        k4_2 = bc(k4_2,BC_TYPE)

        # FINAL RUNGE-KUTTA STEP FOR PSI1   
        psi1[1:-1] = psi1[1:-1] + (1./6)*(k1_1[1:-1] + 2*k2_1[1:-1] + 2*k3_1[1:-1] + k4_1[1:-1])

         # BOUNDARY CONDITION
        psi1 = bc(psi1,BC_TYPE)

        # FINAL RUNGE-KUTTA STEP FOR PSI2
        psi2[1:-1] = psi2[1:-1] + (1./6)*(k1_2[1:-1] + 2*k2_2[1:-1] + 2*k3_2[1:-1] + k4_2[1:-1])

        # BOUNDARY CONDITION
        psi2 = bc(psi2,BC_TYPE)
        
        if (IM_REAL == 0):
            # normalise wavefunction
            Norm1 = 4*pi*np.trapz(r**2*abs(psi1)**2)*dr
            psi1 = psi1/np.sqrt(Norm1) 
            Norm2 = 4*pi*np.trapz(r**2*abs(psi2)**2)*dr
            psi2 = psi2/np.sqrt(Norm2) 
            
            # convergence of chemical potentials
            mu1, mu2 = mu_eqm_dens_ulck(psi1,psi2,r,V1,V2,dr,alpha,beta,eta,N1,N2)
            mu1_tol = np.abs(mu1 - mu1_prev)/np.abs(mu1_prev)
            mu2_tol = np.abs(mu2 - mu2_prev)/np.abs(mu2_prev)
            mu1_prev = mu1
            mu2_prev = mu2
         
        if (IM_REAL == 0 and l % T_SAVE == 0):    
            # convgerence of max densities
            current_mode = N1*np.abs(psi1[1])**2 + N2*np.abs(psi2[1])**2
            tol_mode = np.abs((current_mode - previous_mode)/previous_mode)
            previous_mode = current_mode
            
            # convergence of energies
            E1_ke,E1_pot,E1_int,E1_lhy,E2_ke,E2_pot,E2_int,E2_lhy = energy_eqm_dens_ulck(psi1,psi2,r,V1,V2,dr,alpha,beta,eta,N1,N2)
            E1_total = E1_ke + E1_pot + E1_int + E1_lhy
            E2_total = E2_ke + E2_pot + E2_int + E2_lhy
            E_total = E1_total + E2_total
            E_tol = np.abs(E_total - E_prev)/np.abs(E_prev)
            E_prev = E_total
            
            # print convergence outputs
            print('l = ',l,'(Percentage of imaginary time done = ',100*(l/T_STEPS),'%)')
            print('-'*21,'Max Density Convergence','-'*21)
            print('Max density of component 1 = ',N1*np.abs(psi1[1])**2)
            print('Max density of component 2 = ',N2*np.abs(psi2[1])**2)
            print('Density at max radius of component 1 = ',N1*np.abs(psi1[-1])**2)
            print('Density at max radius of component 2 = ',N2*np.abs(psi2[-1])**2)
            print('Tolerance between successive max densities = ',tol_mode)
            print('-'*17,'Chemical Potential Convergence','-'*17)
            print('Chemical potential of component 1 = ',mu1)
            print('Chemical potential of component 2 = ',mu2)
            print('Tolerance of chemical potential 1 = ',mu1_tol)
            print('Tolerance of chemical potential 2 = ',mu2_tol)
            print('-'*23,'Energy Convergence','-'*23)            
            print('Total energy of component 1 = ',E1_total)
            print('Total energy of component 2 = ',E2_total)
            print('Total combined energies = ',E_total)
            print('Tolerance between successive total energies = ',E_tol)
            print('-'*66)
            
        elif (IM_REAL == 1 and l % T_SAVE == 0):
            print('l = ',l,'(Percentage of real time done = ',100*(l/T_STEPS),'%)')
        
        # save data
        if (l % T_SAVE == 0):
            if IM_REAL == 0:
                # save energies
                E1_array[l//T_SAVE] = E1_total
                E2_array[l//T_SAVE] = E2_total
                # save current time
                t_array[l//T_SAVE] = t.real
                
            elif IM_REAL == 1:
                # save energies
                E1_array[l//T_SAVE] = 0.0
                E2_array[l//T_SAVE] = 0.0
                # save current time
                t_array[l//T_SAVE] = t.imag

            # save wavefunctions
            spacetime1[:,l//T_SAVE] = np.sqrt(N1)*psi1
            spacetime2[:,l//T_SAVE] = np.sqrt(N2)*psi2

            # plotting densities
            plt.plot(r,N1*np.abs(psi1)**2,r,N2*np.abs(psi2)**2,r,N1*np.abs(psi1)**2 + N2*np.abs(psi2)**2)
            plt.xlim((0,np.max(r)))
            plt.ylim(0,1.2*(N1*np.abs(psi1[1])**2 + N2*np.abs(psi2[1])**2))
            plt.xlabel(r'$r$')
            plt.ylabel(r'$n_0(r)$')
            plt.legend((r'$|\psi_1|^2$',r'$|\psi_2|^2$',r'$|\psi_1|^2 + |\psi_2|^2$'))
            plt.close()    

        t = t + dt
        
    if IM_REAL == 0:
        print('Imaginary time finished')
    if IM_REAL == 1:
        print('Real time finished')
        
    return psi1,psi2,mu1,mu2,t_array,spacetime1,spacetime2,E1_array,E2_array
