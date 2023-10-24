import json

import numpy as np

from main.init import grid_setup, derivatives, potential_dens_ulck, init_wavefun_dens_ulck
from main.params_calc import params_dens_ulck
from main.equil_dens import setup_f_interp
from main.rk4_an_ueqm_methods import rk4_an_uneqm_dens_ulck
from main.boundary import absorb_bc_dens_ulck

def time(json_input,PATH):
    
    f = open(json_input,"r")

    setup = json.loads(f.read())
    
    f.close()
    
    # experimental parameters
    m1 = setup['m1']
    m2 = setup['m2']
    a11 = setup['a11']
    a22 = setup['a22']
    a12 = setup['a12']

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
    [Lr, r, dt] = grid_setup(dr, Nr, DT_COEF)
    [Dr, Dr2] = derivatives(dr, r)
    
    # setup initial wavefunctions 
    [psi1_0, psi2_0] = init_wavefun_dens_ulck(r, dr, 
                                              setup['GAUSS_SIGMA1'], setup['GAUSS_SIGMA2'],
                                              setup['INIT_TYPE1'], setup['INIT_TYPE2'])
    # theoretical parameters
    gam, alpha, beta, eta, xi, tau, n01, n02, rho1, rho2, N1_rescale, N2_rescale, z, dim_pot1, dim_pot2 \
    = params_dens_ulck(m1, m2, a11, a22, a12, N1, N2, setup['BALANCE'])
            
    # save theoretical parameters into dictionary
    theory_params = {
    'gamma':gam,
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
    [V1, V2] = potential_dens_ulck(r, setup['OMEGA1'], setup['OMEGA2'], tau, dim_pot1, dim_pot2)

    np.savetxt(PATH + 'potential1.txt', V1, delimiter=',', fmt='%18.16f')
    np.savetxt(PATH + 'potential2.txt', V2, delimiter=',', fmt='%18.16f')
            
    [f_interp, df_interp] = setup_f_interp(z)
    if IM_T_STEPS > 0:
        # imaginary time
        psi1_gs, psi2_gs, mu1_im, mu2_im, t_array_im, E_array_im \
        = rk4_an_uneqm_dens_ulck(r, psi1_0, psi2_0, V1, V2,
                                 gam, z, alpha, beta, eta, 
                                 N1_rescale, N2_rescale,
                                 Dr, Dr2, dr, dt,
                                 IM_T_STEPS, T_SAVE,
                                 0, BC_TYPE, PATH,
                                 f_interp, df_interp)

    psi1_im = psi1_gs
    psi2_im = psi2_gs
            
    # save theoretical parameters into dictionary
    theory_params = {
    'gamma':gam,
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
            V1, V2 = absorb_bc_dens_ulck(r, setup['ABS_HEIGHT'], setup['ABS_SLOPE'], setup['ABS_POS'], setup['ABS_COMP']) 

    
        psi1_re, psi2_re, mu1_re, mu2_re, t_array_re, E1_array_re \
        = rk4_an_uneqm_dens_ulck(r, psi1_im, psi2_im, V1, V2,
                                 gam, z, alpha, beta, eta, 
                                 N1_rescale, N2_rescale,
                                 Dr, Dr2, dr, dt,
                                 RE_T_STEPS, T_SAVE,
                                 1, BC_TYPE, PATH,
                                 f_interp, df_interp)

    #writing data into dictionary
    if IM_T_STEPS > 0 and RE_T_STEPS > 0:
        mix_data = {
        'r' : r,
        't_array_im' : t_array_im,
        't_array_re' : t_array_re,
        'psi1_im' : psi1_gs,
        'psi1_re' : psi1_re,
        'psi2_im' : psi2_gs,
        'psi2_re' : psi2_re,
        'mu1_im' : mu1_im,
        'mu2_im' : mu2_im,
        'E_array_im' : E_array_im
        }
    elif IM_T_STEPS > 0 and RE_T_STEPS == 0:
        mix_data = {
        'r' : r,
        't_array_im' : t_array_im,
        'psi1_im' : psi1_gs,
        'psi2_im' : psi2_gs,
        'mu1_im' : mu1_im,
        'mu2_im' : mu2_im,
        'E_array_im' : E_array_im
        }
    elif IM_T_STEPS == 0 and RE_STEPS > 0 :
        mix_data = {
        'r' : r,
        't_array_re' : t_array_re,
        'psi1_re' : psi1_re,
        'psi2_re' : psi2_re,
        'E_array_re' : E_array_re
        }
        
    return mix_data, theory_params
