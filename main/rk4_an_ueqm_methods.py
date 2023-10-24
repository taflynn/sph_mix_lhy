from main.hamiltonians import ham1_an_uneqm_dens_ulck, ham2_an_uneqm_dens_ulck, ham_lhy_constructor

from main.boundary import bc

from main.chem_pot import mu_an_uneqm_dens_ulck

# from main.energy import energy_eqm_dens_ulck, energy_eqm_dens_lck, energy_uneqm_dens_ulck

from numpy import pi, trapz, absolute, zeros, sqrt, savetxt

def rk4_an_uneqm_dens_ulck(r, psi1, psi2, V1, V2, 
                           gam, z, alpha, beta, eta, N1, N2,
                           Dr, Dr2, dr, dt,
                           T_STEPS, T_SAVE, 
                           IM_REAL, BC_TYPE, PATH,
                           f_interp, df_interp):
    if IM_REAL == 0:
        # initalise Runge-Kutta arrays
        k1_1 = zeros(psi1.size)
        k2_1 = zeros(psi1.size)
        k3_1 = zeros(psi1.size)
        k4_1 = zeros(psi1.size)
        k1_2 = zeros(psi2.size)
        k2_2 = zeros(psi2.size)
        k3_2 = zeros(psi2.size)
        k4_2 = zeros(psi2.size)
        
        E_array = zeros(100)   
        t_array = zeros(100)

        t = 0.0
        mu1 = -1.0
        mu2 = -1.0

    elif IM_REAL == 1:
        # for real time, specify a complex time step
        dt = 1j*dt
        
        print(f'dt = {dt}')
        
        psi1 = psi1.astype(complex)
        psi2 = psi2.astype(complex)
        
        k1_1 = zeros(psi1.size).astype(complex)
        k2_1 = zeros(psi1.size).astype(complex)
        k3_1 = zeros(psi1.size).astype(complex)
        k4_1 = zeros(psi1.size).astype(complex)
        k1_2 = zeros(psi2.size).astype(complex)
        k2_2 = zeros(psi2.size).astype(complex)
        k3_2 = zeros(psi2.size).astype(complex)
        k4_2 = zeros(psi2.size).astype(complex)
        
        E_array = zeros((T_STEPS//T_SAVE))   
        t_array = zeros((T_STEPS//T_SAVE))
        
        t = 0.0 + 0.0j

        mu1 = 0.0
        mu2 = 0.0
    
    # initialise variables used as counters or for convergence calculations    
    previous_mode = 1.0
    mu1_prev = 1.0
    mu2_prev = 1.0
    mu1_tol = 1.0
    mu2_tol = 1.0
    E_prev = 1.0
    frame = 1
    
    for l in range(0,T_STEPS):

        
        # k1 CALCULATION FOR PSI1
        [lhy_ham1, lhy_ham2] = ham_lhy_constructor(psi1, psi2, alpha, beta, f_interp, df_interp)
        [H_ke1, H_trap1, H_int1, H_lhy1, H_mu1] = ham1_an_uneqm_dens_ulck(psi1, psi2, V1, r, Dr, Dr2, dr,
                                                                          gam, z, alpha, beta, eta, 
                                                                          mu1, IM_REAL,
                                                                          lhy_ham1)
    
        k1_1[1:-1] = -dt*(H_ke1[1:-1] + H_trap1[1:-1] + H_lhy1[1:-1] + H_int1[1:-1] + H_mu1[1:-1])

        k1_1 = bc(k1_1, BC_TYPE)

        # k1 CALCULATION FOR PSI2
        [H_ke2, H_trap2, H_int2, H_lhy2, H_mu2] = ham2_an_uneqm_dens_ulck(psi1, psi2, V2, r, Dr, Dr2, dr,
                                                                          gam, z, alpha, beta, eta, 
                                                                          mu2, IM_REAL,
                                                                          lhy_ham2)
    
        k1_2[1:-1] = -dt*(H_ke2[1:-1] + H_trap2[1:-1] + H_lhy2[1:-1] + H_int2[1:-1] + H_mu2[1:-1])

        k1_2 = bc(k1_2, BC_TYPE)
    
        # k2 CALCULATION FOR PSI1
        [lhy_ham1, lhy_ham2] = ham_lhy_constructor(psi1+k1_1/2.0, psi2+k1_2/2.0, alpha, beta, f_interp, df_interp)
        [H_ke1, H_trap1, H_int1, H_lhy1, H_mu1] = ham1_an_uneqm_dens_ulck(psi1+k1_1/2.0, psi2+k1_2/2.0, V1, r, Dr, Dr2, dr,
                                                                          gam, z, alpha, beta, eta, 
                                                                          mu1, IM_REAL,
                                                                          lhy_ham1)
    
        k2_1[1:-1] = -dt*(H_ke1[1:-1] + H_trap1[1:-1] + H_lhy1[1:-1] + H_int1[1:-1] + H_mu1[1:-1])

        k2_1 = bc(k2_1, BC_TYPE)

        # k2 CALCULATION FOR PSI2
        [H_ke2, H_trap2, H_int2, H_lhy2, H_mu2] = ham2_an_uneqm_dens_ulck(psi1+k1_1/2.0, psi2+k1_2/2.0, V2, r, Dr, Dr2, dr,
                                                                          gam, z, alpha, beta, eta, 
                                                                          mu2, IM_REAL,
                                                                          lhy_ham2)
    
        k2_2[1:-1] = -dt*(H_ke2[1:-1] + H_trap2[1:-1] + H_lhy2[1:-1] + H_int2[1:-1] + H_mu2[1:-1])

        k2_2 = bc(k2_2, BC_TYPE)

        # k3 CALCULATION FOR PSI1
        [lhy_ham1, lhy_ham2] = ham_lhy_constructor(psi1+k2_1/2.0, psi2+k2_2/2.0, alpha, beta, f_interp, df_interp)
        [H_ke1, H_trap1, H_int1, H_lhy1, H_mu1] = ham1_an_uneqm_dens_ulck(psi1+k2_1/2.0, psi2+k2_2/2.0, V1, r, Dr, Dr2, dr,
                                                                          gam, z, alpha, beta, eta, 
                                                                          mu1, IM_REAL,
                                                                          lhy_ham1)
    
        k3_1[1:-1] = -dt*(H_ke1[1:-1] + H_trap1[1:-1] + H_lhy1[1:-1] + H_int1[1:-1] + H_mu1[1:-1])

        # BOUNDARY CONDITION
        k3_1 = bc(k3_1, BC_TYPE)

        # k3 CALCULATION FOR PSI2
        [H_ke2, H_trap2, H_int2, H_lhy2, H_mu2] = ham2_an_uneqm_dens_ulck(psi1+k2_1/2.0, psi2+k2_2/2.0, V2, r, Dr, Dr2, dr,
                                                                          gam, z, alpha, beta, eta, 
                                                                          mu2, IM_REAL,
                                                                          lhy_ham2)
    
        k3_2[1:-1] = -dt*(H_ke2[1:-1] + H_trap2[1:-1] + H_lhy2[1:-1] + H_int2[1:-1] + H_mu2[1:-1])

        # BOUNDARY CONDITION
        k3_2 = bc(k3_2, BC_TYPE)

        # k4 CALCULATION FOR PSI1
        [lhy_ham1, lhy_ham2] = ham_lhy_constructor(psi1+k3_1, psi2+k3_2, alpha, beta, f_interp, df_interp)
        [H_ke1, H_trap1, H_int1, H_lhy1, H_mu1] = ham1_an_uneqm_dens_ulck(psi1+k3_1, psi2+k3_2, V1, r, Dr, Dr2, dr,
                                                                          gam, z, alpha, beta, eta, 
                                                                          mu1, IM_REAL,
                                                                          lhy_ham1)
    
        k4_1[1:-1] = -dt*(H_ke1[1:-1] + H_trap1[1:-1] + H_lhy1[1:-1] + H_int1[1:-1] + H_mu1[1:-1])

        # BOUNDARY CONDITION
        k4_1 = bc(k4_1, BC_TYPE)

        # k4 CALCULATION FOR PSI2
        [H_ke2, H_trap2, H_int2, H_lhy2, H_mu2] = ham2_an_uneqm_dens_ulck(psi1+k3_1, psi2+k3_2, V2, r, Dr, Dr2, dr,
                                                                          gam, z, alpha, beta, eta, 
                                                                          mu2, IM_REAL,
                                                                          lhy_ham2)
    
        k4_2[1:-1] = -dt*(H_ke2[1:-1] + H_trap2[1:-1] + H_lhy2[1:-1] + H_int2[1:-1] + H_mu2[1:-1])

        # BOUNDARY CONDITION
        k4_2 = bc(k4_2, BC_TYPE)

        # FINAL RUNGE-KUTTA STEP FOR PSI1   
        psi1[1:-1] = psi1[1:-1] + (1./6)*(k1_1[1:-1] + 2*k2_1[1:-1] + 2*k3_1[1:-1] + k4_1[1:-1])

        # BOUNDARY CONDITION
        psi1 = bc(psi1, BC_TYPE)

        # FINAL RUNGE-KUTTA STEP FOR PSI2
        psi2[1:-1] = psi2[1:-1] + (1./6)*(k1_2[1:-1] + 2*k2_2[1:-1] + 2*k3_2[1:-1] + k4_2[1:-1])

        # BOUNDARY CONDITION
        psi2 = bc(psi2, BC_TYPE)

        if (IM_REAL == 0):
            # normalise wavefunction
            Norm1 = 4*pi*trapz(r**2*abs(psi1)**2)*dr
            psi1 = psi1*sqrt(N1)/sqrt(Norm1) 
            Norm2 = 4*pi*trapz(r**2*abs(psi2)**2)*dr
            psi2 = psi2*sqrt(N2)/sqrt(Norm2) 
            
            # convergence of chemical potentials
            [lhy_ham1, lhy_ham2] = ham_lhy_constructor(psi1, psi2, alpha, beta, f_interp, df_interp)
            [mu1, mu2] = mu_an_uneqm_dens_ulck(psi1, psi2, r, V1, V2, dr,
                                               gam, z, alpha, beta, eta,
                                               lhy_ham1, lhy_ham2)

            mu1_tol = absolute(mu1 - mu1_prev)/absolute(mu1_prev)
            mu2_tol = absolute(mu2 - mu2_prev)/absolute(mu2_prev)
            mu1_prev = mu1
            mu2_prev = mu2
         
        if (IM_REAL == 0 and l % (T_STEPS//100) == 0):    
            # convgerence of max densities
            current_mode = absolute(psi1[1])**2 + absolute(psi2[1])**2
            tol_mode = absolute((current_mode - previous_mode)/previous_mode)
            previous_mode = current_mode
            
            # convergence of energies
            # E_ke, E_pot, E_int, E_lhy = energy_uneqm_dens_ulck(psi1, psi2, r, V1, V2, dr, gam, z, alpha, beta, eta, N1, N2)
            # E_total = E_ke + E_pot + E_int + E_lhy
            # E_tol = absolute(E_total - E_prev)/absolute(E_prev)
            # E_prev = E_total
            
            # save energies
            # E_array[l//(T_STEPS//100)] = E_total
            # save current time
            t_array[l//(T_STEPS//100)] = t.real
            
            # print convergence outputs
            print('l = ', l, ' (Percentage of imaginary time done = ',100*(l/T_STEPS),'%)')
            print('-'*21,'Max Density Convergence','-'*21)
            print('Max density of component 1 = ', absolute(psi1[1])**2)
            print('Max density of component 2 = ', absolute(psi2[1])**2)
            print('Density at max radius of component 1 = ', absolute(psi1[-1])**2)
            print('Density at max radius of component 2 = ', absolute(psi2[-1])**2)
            print(f'Tolerance between successive max densities = {tol_mode}')
            print('-'*17,'Chemical Potential Convergence','-'*17)
            print(f'Chemical potential of component 1 = {mu1}')
            print(f'Chemical potential of component 2 = {mu2}')
            print(f'Tolerance of chemical potential 1 = {mu1_tol}')
            print(f'Tolerance of chemical potential 2 = {mu2_tol}')
            # print('-'*23,'Energy Convergence','-'*23)       
            # print(f'Total combined energies = {E_total}')
            # print(f'Tolerance between successive total energies = {E_tol}')
            print('-'*66)
            
        elif (IM_REAL == 1 and l % T_SAVE == 0):
            # save energies
            E_array[l//T_SAVE] = 0.0
            # save current time
            t_array[l//T_SAVE] = t.imag

            savetxt(PATH + 'psi1_re_t' + str(frame) + '.txt', psi1, delimiter=',', fmt='%18.16f')
            savetxt(PATH + 'psi2_re_t' + str(frame) + '.txt', psi2, delimiter=',', fmt='%18.16f')

            frame = frame + 1
            
            print('l = ', l,' (Percentage of real time done = ',100*(l/T_STEPS),'%)')
            print('-'*21,'Max Densities','-'*21)
            print('Max density of component 1 = ', absolute(psi1[1])**2)
            print('Max density of component 2 = ', absolute(psi2[1])**2)
            print('Density at max radius of component 1 = ', absolute(psi1[-1])**2)
            print('Density at max radius of component 2 = ', absolute(psi2[-1])**2)

        t += dt
        
    if IM_REAL == 0:
        print('Imaginary time finished')
    if IM_REAL == 1:
        print('Real time finished')
        
    return psi1, psi2, mu1, mu2, t_array, E_array
