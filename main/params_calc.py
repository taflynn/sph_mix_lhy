from numpy import pi, sqrt, absolute
from main.equil_dens import eq_dens_lck,eq_dens_ulck
from main.units import units,natural_consts

def params_dens_lck(m1, m2, a11, a22, a12, N):
    """
    This function calculates the necessary dimensionless parameter for the density-locked mixture. 
    
    The outputs are hence:
    (1) N_lck - the dimensionless effective atom number appearing in the Gross-Pitaevskii (GP) equation
    (2) xi, tau - the length and timescales that define the system
    (3) n01, n02 - the equilibrium densities for each respective component
    
    The only inputs thar are neeeded are
    -> m1,m2 - the two component masses
    -> a11,a22 - the two intraspecies scattering lengths
    -> a12 - the interspecies scattering length
    -> N1,N2 - the two component atom numbers
    """
    [hbar, a0, Da] = natural_consts()
    a11 = a11*a0; a22 = a22*a0; a12 = a12*a0
    m1 = m1*Da; m2 = m2*Da
    # equilibrium densities equal masses
    n01, n02 = eq_dens_lck(m1, m2, a11, a22, a12)
    
    # length and time scales
    xi, tau = units(m1, m2, a11, a22, a12, n01)
    
    # density-locked effective atom number
    if m1 == m2:
    	
        delta_a = a12 + sqrt(a11*a22)
        dim_pot = (m1*tau*xi**2)/hbar
    
        N_lck = (N/(n01*xi**3))*(sqrt(a22)/(sqrt(a11) + sqrt(a22)))
        print('Balanced experimental N1 = ', N*(sqrt(a22)/(sqrt(a11) + sqrt(a22))))
        print('Balanced experimental N2 = ', N*(sqrt(a11)/(sqrt(a11) + sqrt(a22))))
    
    elif m1 != m2:
        # interaction strengths
        g11 = 4*pi*hbar**2*a11/m1
        g12 = 2*pi*hbar**2*a12*(1/m1 + 1/m2)
        g22 = 4*pi*hbar**2*a22/m2
        
        dim_pot = (m1*m2/(m1 + m2))*(tau*xi**2)/hbar
        
        N_lck = (N/(n01*xi**3))*(sqrt(g22)/(sqrt(g11) + sqrt(g22)))
        print('Balanced experimental N1 = ', N*(sqrt(g22)/(sqrt(g11) + sqrt(g22))))
        print('Balanced experimental N2 = ', N*(sqrt(g11)/(sqrt(g11) + sqrt(g22))))
 
    print('Calculating defining parameters and scales for density locked mixture:')
    print(f'Effective atom number of density-locked mixture, N_lck = {N_lck}')
    print(f'Lengthscale, xi = {xi}')
    print(f'Timescale, tau = {tau}')
    print(f'Equilibrium density of component 1, n01 = {n01}')
    print(f'Equilibrium density of component 2, n02 = {n02}')
    print(67*'-')
    return N_lck, xi, tau, n01, n02, dim_pot

def params_dens_ulck(m1, m2, a11, a22, a12, N1, N2, BALANCE):
    """
    This function calculates the necessary dimensionless parameters for the density-unlocked mixture (though this is 
    currently only for the equal masses mixture). The outputs are hence:
    * alpha, beta, eta - the three dimensionless parameters appearing in the Gross-Pitaevskii (GP) equation
    * xi, tau - the length and timescales that define the system
    * n01, n02 - the equilibrium densities for each component
    * rho1, rho2 - the density scalings for the two component wavefunctions
    * N1, N2 - the rescaled atom numbers for each component
    
    The only inputs that are neeeded are
    -> m1, m2 - the two component masses
    -> a11, a22 - the two intraspecies scattering lengths
    -> a12 - the interspecies scattering length
    -> N1, N2 - the two component atom numbers
    """
    
    # natural constants
    [hbar, a0, Da] = natural_consts()
    
    if m1 == m2:
        a11 = a11*a0; a22 = a22*a0; a12 = a12*a0
        m1 = m1*Da; m2 = m2*Da
        # equilibrium densities equal masses solving coupled equations
        n01, n02 = eq_dens_ulck(m1, m2, a11, a22, a12)
        
        # equilibrium densities equal masses with density-lock
        #n01,n02 = eq_dens_lck(m1,m2,a11,a22,a12)
        
        delta_a = a12 + sqrt(a11*a22)
        
        alpha = (32/3)*sqrt((2/(3*pi))*(absolute(delta_a)*a11**2.5*n01)/(sqrt(a11) + sqrt(a22)))
        beta = sqrt(a22/a11)
        eta = a12/sqrt(a11*a22)
        
        rho1 = (2/3)*((absolute(delta_a)*n01)/(sqrt(a11)*(sqrt(a11) + sqrt(a22))))
        rho2 = (2/3)*((absolute(delta_a)*n01)/(sqrt(a22)*(sqrt(a11) + sqrt(a22))))
        
        xi, tau = units(m1, m2, a11, a22, a12, n01)
        
        if BALANCE == 1:
            print('Balancing N2 to the value of N1')
            N2 = N1*sqrt(a11/a22)
        elif BALANCE == 2:
            print('Balancing N1 to the value of N2')
            N1 = N2*sqrt(a22/a11)

        dim_pot = (m1*tau*xi**2)/hbar

        N1 = N1/(rho1*xi**3)
        N2 = N2/(rho2*xi**3)
        print('Calculating defining parameters and scales for density-unlocked mixture (equal masses):')
        print('Dimensionless parameters of mixture:')
        print(f'alpha = {alpha}, beta = {beta}, eta = {eta}')
        print(f'Lengthscale, xi = {xi}')
        print(f'Timescale, tau = {tau}')
        print(f'Equilibrium density of component 1, n01 = {n01}')
        print(f'Equilibrium density of component 2, n02 = {n02}')
        print(67*'-')
        return alpha, beta, eta, xi, tau, n01, n02, rho1, rho2, N1, N2, dim_pot
        
    elif m1 != m2:
        # equilibrium densities unequal masses
        z = m2/m1
        a11 = a11*a0; a22 = a22*a0; a12 = a12*a0
        m1 = m1*Da; m2 = m2*Da    
        
        # equilibrium densities equal masses with density-lock
        n01,n02 = eq_dens_lck(m1,m2,a11,a22,a12)
        
        # defining effective interaction strenghts
        g11 = 4*pi*hbar**2*a11/m1
        g22 = 4*pi*hbar**2*a22/m2
        g12 = 2*pi*hbar**2*a12*(1/m1 + 1/m2)
        
        delta_g = g12 + sqrt(g11*g22)
        
        gam = (sqrt(g11) + sqrt(g22))/(m1*(sqrt(g11)/m2 + sqrt(g22)/m1))
        alpha = (8/(15*pi**2))*sqrt((2/3)*(m1/hbar**2)**3*(absolute(delta_g)*g11**2.5*n01)/(sqrt(g11) + sqrt(g22)))
        beta = sqrt(g22/g11)
        eta = g12/sqrt(g11*g22)
        
        rho1 = (2/3)*((absolute(delta_g)*n01)/(sqrt(g11)*(sqrt(g11) + sqrt(g22))))
        rho2 = (2/3)*((absolute(delta_g)*n01)/(sqrt(g22)*(sqrt(g11) + sqrt(g22))))

        xi, tau = units(m1, m2, a11, a22, a12, n01)
        
        if BALANCE == 1:
            print('Balancing N2 to the value of N1')
            N2 = N1*sqrt(g11/g22)
            print(f'Experimental N2 = {N2}')
        elif BALANCE == 2:
            print('Balancing N1 to the value of N2')
            N1 = N2*sqrt(g22/g11)
            print(f'Experimental N1 = {N1}')

        dim_pot1 = (m1*xi**2*tau)/hbar#(sqrt(g22) + (m1/m2)*sqrt(g11))/(sqrt(g11) + sqrt(g22))
        dim_pot2 = (m2*xi**2*tau)/hbar#((m2/m1)*sqrt(g22) + sqrt(g11))/(sqrt(g11) + sqrt(g22))

        N1 = N1/(rho1*xi**3)
        N2 = N2/(rho2*xi**3)
        print(f'Calculating defining parameters and scales for density-unlocked mixture (mass ratio z = {z}):')
        print('Dimensionless parameters of mixture:')
        print(f'gamma = {gam}')
        print(f'alpha = {alpha}, beta = {beta}, eta = {eta}')
        print(f'Lengthscale, xi = {xi}')
        print(f'Timescale, tau = {tau}')
        print(f'Equilibrium density of component 1, n01 = {n01}')
        print(f'Equilibrium density of component 2, n02 = {n02}')
        print(67*'-')
        return gam, alpha, beta, eta, xi, tau, n01, n02, rho1, rho2, N1, N2, z, dim_pot1, dim_pot2
