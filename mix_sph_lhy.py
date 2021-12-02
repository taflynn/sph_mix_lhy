import json

from init import grid_setup,potential_dens_lck,potential_dens_ulck,init_wavefun_dens_ulck,init_wavefun_dens_lck
from params_calc import params_dens_lck, params_dens_ulck
from rk4_methods import rk4_dens_lck, rk4_dens_ulck

def time(json_input):
    
    f = open(json_input,"r")

    setup = json.loads(f.read())
    
    # locking?
    DENS_LCK = setup['DENS_LCK']
    
    # experimental parameters
    m1 = setup['m1']
    m2 = setup['m2']
    a11 = setup['a11']
    a22 = setup['a22']
    a12 = setup['a12']
    if DENS_LCK == 1:
        N = setup['N']
    elif DENS_LCK == 0:
        N1 = setup['N1']
        N2 = setup['N2']
    
    # grid
    dr = setup['dr']
    Nr = setup['Nr']
    
    # boundary conditions
    BC_TYPE = setup['BC_TYPE']
    
    # time-stepping
    DT_COEF = setup['DT_COEF']
    IM_T_STEPS = setup['IM_T_STEPS']
    RE_T_STEPS = setup['RE_T_STEPS']
    T_SAVE = setup['T_SAVE']
    
    # configure grid
    [Lr,r,dt] = grid_setup(dr,Nr,DT_COEF)
    
    if DENS_LCK == 1:
        # setup initial wavefunctions 
        phi_0 = init_wavefun_dens_lck(r,dr,setup['GAUSS_SIGMA'],setup['INIT_TYPE'])
        
        # setup trapping potential
        V = potential_dens_lck(r,setup['OMEGA'])
        
        # theoretical parameters
        [N_lck,xi,tau,n01,n02] = params_dens_lck(m1,m2,a11,a22,a12,N)

        if IM_T_STEPS > 0:
            # imaginary time
            phi_im,mu_im,t_array_im,spacetime_im,E_array_im = rk4_dens_lck(r,phi_0,V,N_lck,dr,dt,IM_T_STEPS,T_SAVE,0,BC_TYPE)
        
        if RE_T_STEPS > 0:
            if BREATH == 1:
                phi_im = phi_im*np.exp(1.0j*lamb*r**2)
            # real time
            phi_re,mu_re,t_array_re,spacetime_re,E_array_re = rk4_dens_lck(r,phi_im,V,N_lck,dr,dt,RE_T_STEPS,T_SAVE,1,BC_TYPE)
        
        # writing data into dictionary
        if IM_T_STEPS > 0 and RE_T_STEPS > 0:
            mix_data = {
            'r' : r,
            't_array_im' : t_array_im,
            't_array_re' : t_array_re,
            'phi_im' : phi_im,
            'phi_re' : phi_re,
            'mu_im' : mu_im,
            'E_array_im' : E_array_im,
            'spacetime_im' : spacetime_im,
            'spacetime_re' : spacetime_re
            }
        elif IM_T_STEPS > 0 and RE_T_STEPS == 0:
            mix_data = {
            'r' : r,
            't_array_im' : t_array_im,
            'phi_im' : phi_im,
            'mu_im' : mu_im,
            'E_array_im' : E_array_im,
            'spacetime_im' : spacetime_im,
            }
        elif IM_T_STEPS == 0 and RE_STEPS > 0 :
            mix_data = {
            'r' : r,
            't_array_re' : t_array_re,
            'phi_re' : phi_re,
            'spacetime_re' : spacetime_re               
            } 

    elif DENS_LCK == 0:
        if m1 == m2:
            # setup initial wavefunctions 
            [psi1_0,psi2_0] = init_wavefun_dens_ulck(r,dr,setup['GAUSS_SIGMA'],setup['INIT_TYPE1'],setup['INIT_TYPE2'])

            # setup trapping potentials
            [V1,V2] = potential_dens_ulck(r,setup['OMEGA1'],setup['OMEGA2'])

            # theoretical parameters
            alpha,beta,eta,xi,tau,n01,n02,rho1,rho2,N1_rescale,N2_rescale \
            = params_dens_ulck(m1,m2,a11,a22,a12,N1,N2)
            
            if IM_T_STEPS > 0:
                # imaginary time
                psi1_im,psi2_im,mu1_im,mu2_im,t_array_im,spacetime1_im,spacetime2_im,E1_array_im,E2_array_im \
            = rk4_dens_ulck(r,psi1_0,psi2_0,V1,V2,alpha,beta,eta,N1_rescale,N2_rescale,dr,dt,IM_T_STEPS,T_SAVE,0,BC_TYPE)
            
            if RE_T_STEPS > 0:
                if setup['BREATH1'] == 1:
                    psi1_im = psi1_im*np.exp(1.0j*lamb*r**2)
                if setup['BREATH2'] == 1:
                    psi2_im = psi2_im*np.exp(1.0j*lamb*r**2)
                    
                # real time
            psi1_re,psi2_re,mu1_re,mu2_re,t_array_re,spacetime1_re,spacetime2_re,E1_array_re,E2_array_re \
            = rk4_dens_ulck(r,psi1_im,psi2_im,V1,V2,alpha,beta,eta,N1_rescale,N2_rescale,dr,dt,RE_T_STEPS,T_SAVE,1,BC_TYPE)
        
        elif m1 != m2:
            print('Density-unlocked, unequal masses is still being worked on!')
    
        # writing data into dictionary
        if IM_T_STEPS > 0 and RE_T_STEPS > 0:
            mix_data = {
            'r' : r,
            't_array_im' : t_array_im,
            't_array_re' : t_array_re,
            'psi1_im' : psi1_im,
            'psi1_re' : psi1_re,
            'psi2_im' : psi2_im,
            'psi2_re' : psi2_re,
            'mu1_im' : mu1_im,
            'mu2_im' : mu2_im,
            'E1_array_im' : E1_array_im,
            'E2_array_im' : E2_array_im,
            'spacetime1_im' : spacetime1_im,
            'spacetime2_im' : spacetime2_im,
            'spacetime1_re' : spacetime1_re,
            'spacetime2_re' : spacetime2_re
            }
        elif IM_T_STEPS > 0 and RE_T_STEPS == 0:
            mix_data = {
            'r' : r,
            't_array_im' : t_array_im,
            'psi1_im' : psi1_im,
            'psi2_im' : psi2_im,
            'mu1_im' : mu1_im,
            'mu2_im' : mu2_im,
            'E1_array_im' : E1_array_im,
            'E2_array_im' : E2_array_im,
            'spacetime1_im' : spacetime1_im,
            'spacetime2_im' : spacetime2_im
            }
        elif IM_T_STEPS == 0 and RE_STEPS > 0 :
            mix_data = {
            'r' : r,
            't_array_re' : t_array_re,
            'psi1_re' : psi1_re,
            'psi2_re' : psi2_re,
            'E1_array_im' : E1_array_im,
            'E2_array_im' : E2_array_im,
            'spacetime1_re' : spacetime1_re,
            'spacetime2_re' : spacetime2_re
            }
        
    return mix_data
"""
    # parse arguments
    p = argparse.ArgumentParser()
    # REQUIRED:
    # masses
    p.add_argument('--m1', dest ='m1', required = True, nargs=1)
    p.add_argument('--m2', dest ='m2', required = True, nargs=1)
    # scattering lengths
    p.add_argument('--a11', dest ='a11', required = True, nargs=1)
    p.add_argument('--a22', dest ='a22', required = True, nargs=1)
    p.add_argument('--a12', dest ='a12', required = True, nargs=1)
    # grid spacing
    p.add_argument('--dr', dest ='dr', required = True, nargs=1)
    p.add_argument('--Nr', dest ='Nr', required = True, nargs=1)
    # density-locked?
    p.add_argument('--DENS_LCK', dest ='DENS_LCK', required = True, nargs=1)
    # time-stepping
    p.add_argument('--DT_COEF', dest ='DT_COEF', required = True, nargs=1)
    p.add_argument('--IM_T_STEPS', dest ='IM_T_STEPS', required = True, nargs=1)
    p.add_argument('--RE_T_STEPS', dest ='RE_T_STEPS', required = True, nargs=1)
    p.add_argument('--T_SAVE', dest ='T_SAVE', required = True, nargs=1)
    # atom numbers
    p.add_argument('N1', dest ='N1', required = True, nargs=1)
    p.add_argument('N2', dest ='N2', required = True, nargs=1)
    
    # NOT REQUIRED
    # boundary conditions
    p.add_argument('BC_TYPE', dest ='BC_TYPE', required = False, default=0, nargs=1)
    
    # initial condition
    p.add_argument('GAUSS_SIGMA', dest ='GAUSS_SIGMA', required = False, default=10.0, nargs=1)
    if DENS_LCK == 1:
        p.add_argument('INIT_TYPE', dest ='INIT_TYPE', required = False, default='S_GAUSS', nargs=1)
        # trapping potential
        p.add_argument('OMEGA', dest ='OMEGA', required = False, default=0.0, nargs=1)
    elif DENS_LCK == 0:
        p.add_argument('INIT_TYPE1', dest ='INIT_TYPE1', required = False, default='NON_ZERO_TAIL', nargs=1)
        p.add_argument('INIT_TYPE2', dest ='INIT_TYPE2', required = False, default='S_GAUSS', nargs=1)
        # trapping potential
        p.add_argument('OMEGA1', dest ='OMEGA1', required = False, default=0.0, nargs=1)
        p.add_argument('OMEGA2', dest ='OMEGA2', required = False, default=0.0, nargs=1)
    
    # breathing mode perturbation
    p.add_argument('lamb', dest ='lamb', required = False, default=0.0, nargs=1)

    args = p.parse_args()
"""