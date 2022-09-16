import json

import numpy as np

from main.init import grid_setup,potential_dens_lck,potential_dens_ulck,init_wavefun_dens_ulck,init_wavefun_dens_lck
from main.params_calc import params_dens_lck, params_dens_ulck
from main.rk4_methods import rk4_dens_lck, rk4_eqm_dens_ulck, rk4_uneqm_dens_ulck
from main.boundary import absorb_bc_dens_lck,absorb_bc_dens_ulck

def time(json_input,path):
    
    f = open(json_input,"r")

    setup = json.loads(f.read())
    
    f.close()
    
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
        
        # theoretical parameters
        [N_lck,xi,tau,n01,n02] = params_dens_lck(m1,m2,a11,a22,a12,N)
       
        # setup trapping potential
        V = potential_dens_lck(r,setup['OMEGA'],tau)

        if IM_T_STEPS > 0:
            # imaginary time
            phi_im,mu_im,t_array_im,spacetime_im,E_array_im = rk4_dens_lck(r,phi_0,V,N_lck,dr,dt,IM_T_STEPS,T_SAVE,0,BC_TYPE)
        
            # save theoretical parameters into dictionary
            theory_params = {
            'N_lck':N_lck,
            'xi':xi,
            'tau':tau,
            'n01':n01,
            'n02':n02,
            'mu':mu_im
            }

        if RE_T_STEPS > 0:
            epsilon = setup['epsilon']
            if setup['BREATH'] == 1:
                phi_im = phi_im*np.exp(1.0j*epsilon*r**2)
            if setup['ABSORB_BC'] == 1:
                V = absorb_bc_dens_lck(r,setup['ABS_HEIGHT'],setup['ABS_SLOPE'],setup['ABS_POS']) 
            # real time
            phi_re,mu_re,t_array_re,spacetime_re,E_array_re = rk4_dens_lck(r,phi_im,V,N_lck,dr,dt,RE_T_STEPS,T_SAVE,1,BC_TYPE)
       
        # writing data into dictionary
        if IM_T_STEPS > 0 and RE_T_STEPS > 0:
            mix_data = {
            'r' : r,
            't_array_im' : t_array_im,
            't_array_re' : t_array_re,
            'phi_im' : N_lck**0.5*phi_im,
            'phi_re' : N_lck**0.5*phi_re,
            'mu_im' : mu_im,
            'E_array_im' : E_array_im,
            'spacetime_im' : spacetime_im,
            'spacetime_re' : spacetime_re
            }
        elif IM_T_STEPS > 0 and RE_T_STEPS == 0:
            mix_data = {
            'r' : r,
            't_array_im' : t_array_im,
            'phi_im' : N_lck**0.5*phi_im,
            'mu_im' : mu_im,
            'E_array_im' : E_array_im,
            'spacetime_im' : spacetime_im,
            }
        elif IM_T_STEPS == 0 and RE_STEPS > 0 :
            mix_data = {
            'r' : r,
            't_array_re' : t_array_re,
            'phi_re' : N_lck**0.5*phi_re,
            'spacetime_re' : spacetime_re               
            } 

    elif DENS_LCK == 0:
        # setup initial wavefunctions 
        [psi1_0,psi2_0] = init_wavefun_dens_ulck(r,dr,setup['GAUSS_SIGMA1'],setup['GAUSS_SIGMA2'],setup['INIT_TYPE1'],setup['INIT_TYPE2'])
        if m1 == m2:

            # theoretical parameters
            alpha,beta,eta,xi,tau,n01,n02,rho1,rho2,N1_rescale,N2_rescale \
            = params_dens_ulck(m1,m2,a11,a22,a12,N1,N2,setup['BALANCE'])
            
            # setup trapping potentials
            [V1,V2] = potential_dens_ulck(r,setup['OMEGA1'],setup['OMEGA2'],tau)
            
            if IM_T_STEPS > 0:
                # imaginary time
                psi1_gs,psi2_gs,mu1_im,mu2_im,t_array_im,spacetime1_im,spacetime2_im,E_array_im \
            = rk4_eqm_dens_ulck(r,psi1_0,psi2_0,V1,V2,alpha,beta,eta,N1_rescale,N2_rescale,dr,dt,IM_T_STEPS,T_SAVE,0,BC_TYPE,path)
            psi1_im = psi1_gs
            psi2_im = psi2_gs
            
            # save theoretical parameters into dictionary
            theory_params = {
            'alpha':alpha,
            'beta':beta,
            'eta':eta,
            'xi':xi,
            'tau':tau,
            'n01':n01,
            'n02':n02,
            'rho1':rho1,
            'rho2':rho2,
            'N1':N1_rescale,
            'N2':N2_rescale,
            'mu1':mu1_im,
            'mu2':mu2_im
            }

            if RE_T_STEPS > 0:
                epsilon = setup['epsilon']
                if setup['BREATH1'] == 1:
                    psi1_im = psi1_gs*np.exp(1.0j*epsilon*r**2)
                if setup['BREATH2'] == 1:
                    psi2_im = psi2_gs*np.exp(1.0j*epsilon*r**2)
                if setup['ABSORB_BC'] == 1:
                    V1,V2 = absorb_bc_dens_ulck(r,setup['ABS_HEIGHT'],setup['ABS_SLOPE'],setup['ABS_POS'],setup['ABS_COMP']) 
                # real time
            psi1_re,psi2_re,mu1_re,mu2_re,t_array_re,spacetime1_re,spacetime2_re,E1_array_re \
            = rk4_eqm_dens_ulck(r,psi1_im,psi2_im,V1,V2,alpha,beta,eta,N1_rescale,N2_rescale,dr,dt,RE_T_STEPS,T_SAVE,1,BC_TYPE,path)
        
        elif m1 != m2:
            # theoretical parameters
            gam1,gam2,alpha,beta,eta,xi,tau,n01,n02,rho1,rho2,N1_rescale,N2_rescale,z \
            = params_dens_ulck(m1,m2,a11,a22,a12,N1,N2,setup['BALANCE'])
            
            # save theoretical parameters into dictionary
            theory_params = {
            'gamma1':gam1,
            'gamma2':gam2,
            'alpha':alpha,
            'beta':beta,
            'eta':eta,
            'xi':xi,
            'tau':tau,
            'n01':n01,
            'n02':n02,
            'rho1':rho1,
            'rho2':rho2,
            'N1':N1_rescale,
            'N2':N2_rescale,
            'z':z
            }
            
            # setup trapping potentials
            [V1,V2] = potential_dens_ulck(r,setup['OMEGA1'],setup['OMEGA2'],tau)
            
            if IM_T_STEPS > 0:
                # imaginary time
                psi1_gs,psi2_gs,mu1_im,mu2_im,t_array_im,spacetime1_im,spacetime2_im,E_array_im \
            = rk4_uneqm_dens_ulck(r,psi1_0,psi2_0,V1,V2,gam1,gam2,alpha,beta,eta,N1_rescale,N2_rescale,dr,dt,IM_T_STEPS,T_SAVE,0,BC_TYPE,z)
            psi1_im = psi1_gs
            psi2_im = psi2_gs
            
            # save theoretical parameters into dictionary
            theory_params = {
            'gamma1':gam1,
            'gamma2':gam2,
            'alpha':alpha,
            'beta':beta,
            'eta':eta,
            'xi':xi,
            'tau':tau,
            'n01':n01,
            'n02':n02,
            'rho1':rho1,
            'rho2':rho2,
            'N1':N1_rescale,
            'N2':N2_rescale,
            'mu1':mu1_im,
            'mu2':mu2_im,
            'z':z
            }

            if RE_T_STEPS > 0:
                epsilon = setup['epsilon']
                if setup['BREATH1'] == 1:
                    psi1_im = psi1_gs*np.exp(1.0j*epsilon*r**2)
                if setup['BREATH2'] == 1:
                    psi2_im = psi2_gs*np.exp(1.0j*epsilon*r**2)
                if setup['ABSORB_BC'] == 1:
                    V1,V2 = absorb_bc_dens_ulck(r,setup['ABS_HEIGHT'],setup['ABS_SLOPE'],setup['ABS_POS'],setup['ABS_COMP']) 
                # real time
            psi1_re,psi2_re,mu1_re,mu2_re,t_array_re,spacetime1_re,spacetime2_re,E1_array_re \
            = rk4_uneqm_dens_ulck(r,psi1_im,psi2_im,V1,V2,gam1,gam2,alpha,beta,eta,N1_rescale,N2_rescale,dr,dt,RE_T_STEPS,T_SAVE,1,BC_TYPE,z)

        # writing data into dictionary
        if IM_T_STEPS > 0 and RE_T_STEPS > 0:
            mix_data = {
            'r' : r,
            't_array_im' : t_array_im,
            't_array_re' : t_array_re,
            'psi1_im' : N1_rescale**0.5*psi1_gs,
            'psi1_re' : N1_rescale**0.5*psi1_re,
            'psi2_im' : N2_rescale**0.5*psi2_gs,
            'psi2_re' : N2_rescale**0.5*psi2_re,
            'mu1_im' : mu1_im,
            'mu2_im' : mu2_im,
            'E_array_im' : E_array_im,
            'spacetime1_im' : spacetime1_im,
            'spacetime2_im' : spacetime2_im,
            'spacetime1_re' : spacetime1_re,
            'spacetime2_re' : spacetime2_re
            }
        elif IM_T_STEPS > 0 and RE_T_STEPS == 0:
            mix_data = {
            'r' : r,
            't_array_im' : t_array_im,
            'psi1_im' : N1_rescale**0.5*psi1_gs,
            'psi2_im' : N2_rescale**0.5*psi2_gs,
            'mu1_im' : mu1_im,
            'mu2_im' : mu2_im,
            'E_array_im' : E_array_im,
            'spacetime1_im' : spacetime1_im,
            'spacetime2_im' : spacetime2_im
            }
        elif IM_T_STEPS == 0 and RE_STEPS > 0 :
            mix_data = {
            'r' : r,
            't_array_re' : t_array_re,
            'psi1_re' : N1_rescale**0.5*psi1_re,
            'psi2_re' : N2_rescale**0.5*psi2_re,
            'E_array_re' : E_array_re,
            'spacetime1_re' : spacetime1_re,
            'spacetime2_re' : spacetime2_re
            }
        
    return mix_data,theory_params
