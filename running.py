from mix_sph_lhy import time

import matplotlib.pyplot as plt

import numpy as np

import json

import os

fname = 'config_dens_ulck.json'
dirarg = 'imbalance'

path = os.path.join('./data',dirarg+'/')
if not os.path.isdir(path):
    os.mkdir(path)

f = open(fname,"r")
setup = json.loads(f.read())
f.close()
mix_data = time(fname)

# extracting data from dictionary
r = mix_data['r']
if setup['DENS_LCK'] == 1:
    if setup['IM_T_STEPS'] > 0 and setup['RE_T_STEPS'] > 0:
        t_array_im = mix_data['t_array_im']
        t_array_re = mix_data['t_array_re']
        phi_im = mix_data['phi_im']
        phi_re = mix_data['phi_re']
        mu_im = mix_data['mu_im']
        E_array_im = mix_data['E_array_im']
        spacetime_im = mix_data['spacetime_im']
        spacetime_re = mix_data['spacetime_re']
        
        # output files
        im_dens = np.column_stack((r,np.abs(phi_im)**2))
        np.savetxt(path + 'imag_fin_dens.csv',im_dens,delimiter=',',fmt='%18.16f')
        re_dens = np.column_stack((r,np.abs(phi_re)**2))
        np.savetxt(path + 'real_fin_dens.csv',re_dens,delimiter=',',fmt='%18.16f')
        E_im = np.column_stack((t_array_im,E_array_im))
        np.savetxt(path + 'energy_imag.csv',E_im,delimiter=',',fmt='%18.16f')
        
        [T_IM,R] = np.meshgrid(t_array_im,r)
        
        plt.pcolormesh(T_IM,R,np.abs(spacetime_im)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'imag_tot_dens.png',dpi='figure')
        plt.close()    


        [R,T_RE] = np.meshgrid(r,t_array_re)
        
        plt.pcolormesh(T_RE,R,np.abs(spacetime_re)**2,shading='gourand')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'real_tot_dens.png',dpi='figure')
        plt.close()    

    elif setup['IM_T_STEPS'] > 0 and setup['RE_T_STEPS'] == 0:
        t_array_im = mix_data['t_array_im']
        phi_im = mix_data['phi_im']
        mu_im = mix_data['mu_im']
        E_array_im = mix_data['E_array_im']
        spacetime_im = mix_data['spacetime_im']
        
        # output files
        im_dens = np.column_stack((r,np.abs(phi_im)**2))
        np.savetxt(path + 'imag_fin_dens.csv',im_dens,delimiter=',',fmt='%18.16f')
        E_im = np.column_stack((t_array_im,E_array_im))
        np.savetxt(path + 'energy_imag.csv',E_im,delimiter=',',fmt='%18.16f')
        plt.close()    

        [T_IM,R] = np.meshgrid(t_array_im,r)

        plt.pcolormesh(T_IM,R,np.abs(spacetime_im)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'imag_tot_dens.png',dpi='figure')
        plt.close()    

    elif setup['IM_T_STEPS'] == 0 and setup['RE_T_STEPS'] > 0:
        t_array_re = mix_data['t_array_re']
        phi_re = mix_data['phi_re']
        spacetime_re = mix_data['spacetime_re']
        
        # output files
        re_dens = np.column_stack((r,np.abs(phi_re)**2))
        np.savetxt(path + 'real_fin_dens.csv',re_dens,delimiter=',',fmt='%18.16f')

        [T_RE,R] = np.meshgrid(t_array_re,r)

        plt.pcolormesh(T_RE,R,np.abs(spacetime_re)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'real_tot_dens.png',dpi='figure')
        plt.close()    

elif setup['DENS_LCK'] == 0:
    if setup['IM_T_STEPS'] > 0 and setup['RE_T_STEPS'] > 0:    
        t_array_im = mix_data['t_array_im']
        t_array_re = mix_data['t_array_re']
        psi1_im = mix_data['psi1_im']
        psi1_re = mix_data['psi1_re']
        psi2_im = mix_data['psi2_im']
        psi2_re = mix_data['psi2_re']
        mu1_im = mix_data['mu1_im']
        mu2_im = mix_data['mu2_im']
        E1_array_im = mix_data['E1_array_im']
        E2_array_im = mix_data['E2_array_im']
        spacetime1_im = mix_data['spacetime1_im']
        spacetime2_im = mix_data['spacetime2_im']
        spacetime1_re = mix_data['spacetime1_re']
        spacetime2_re = mix_data['spacetime2_re']
        
        # output files
        im_dens1 = np.column_stack((r,np.abs(psi1_im)**2))
        np.savetxt(path + 'imag_fin_dens1.csv',im_dens1,delimiter=',',fmt='%18.16f')
        re_dens1 = np.column_stack((r,np.abs(psi1_re)**2))
        np.savetxt(path + 'real_fin_dens1.csv',re_dens1,delimiter=',',fmt='%18.16f')
        
        im_dens2 = np.column_stack((r,np.abs(psi2_im)**2))
        np.savetxt(path + 'imag_fin_dens2.csv',im_dens2,delimiter=',',fmt='%18.16f')
        re_dens2 = np.column_stack((r,np.abs(psi2_re)**2))
        np.savetxt(path + 'real_fin_dens2.csv',re_dens2,delimiter=',',fmt='%18.16f')
        
        im_dens_tot = np.column_stack((r,np.abs(psi1_im)**2) + np.abs(psi2_im)**2)
        np.savetxt(path + 'imag_fin_dens_tot.csv',im_dens_tot,delimiter=',',fmt='%18.16f')
        re_dens_tot = np.column_stack((r,np.abs(psi1_re)**2) + np.abs(psi2_re)**2)
        np.savetxt(path + 'real_fin_dens_tot.csv',re_dens_tot,delimiter=',',fmt='%18.16f')
        
        E_tot_im = np.column_stack((t_array_im,E1_array_im + E2_array_im))
        np.savetxt(path + 'tot_energy_imag.csv',E_im,delimiter=',',fmt='%18.16f')
        
        [T_IM,R] = np.meshgrid(t_array_im,r)

        plt.pcolormesh(T_IM,R,np.abs(spacetime1_im)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'imag_tot_dens1.png',dpi='figure')
        plt.close()    
        
        plt.pcolormesh(T_IM,R,np.abs(spacetime2_im)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'imag_tot_dens2.png',dpi='figure')
        plt.close()    
        
        plt.pcolormesh(T_IM,R,np.abs(spacetime1_im)**2 + np.abs(spacetime2_im)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'imag_tot_dens_tot.png',dpi='figure')
        plt.close()    
        
        [T_RE,R] = np.meshgrid(t_array_re,r)
        
        plt.pcolormesh(T_RE,R,np.abs(spacetime1_re)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'real_tot_dens1.png',dpi='figure')
        plt.close()    
        
        plt.pcolormesh(T_RE,R,np.abs(spacetime2_re)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'real_tot_dens2.png',dpi='figure')
        plt.close()    
        
        plt.pcolormesh(T_RE,R,np.abs(spacetime_1_re)**2 + np.abs(spacetime2_re)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'real_tot_dens_tot.png',dpi='figure')
        plt.close()    
       
    elif setup['IM_T_STEPS'] > 0 and setup['RE_T_STEPS'] == 0:    
        t_array_im = mix_data['t_array_im']
        psi1_im = mix_data['psi1_im']
        psi2_im = mix_data['psi2_im']
        mu1_im = mix_data['mu1_im']
        mu2_im = mix_data['mu2_im']
        E1_array_im = mix_data['E1_array_im']
        E2_array_im = mix_data['E2_array_im']
        spacetime1_im = mix_data['spacetime1_im']
        spacetime2_im = mix_data['spacetime2_im']
        
        # output files
        im_dens1 = np.column_stack((r,np.abs(psi1_im)**2))
        np.savetxt(path + 'imag_fin_dens1.csv',im_dens1,delimiter=',',fmt='%18.16f')
        
        im_dens2 = np.column_stack((r,np.abs(psi2_im)**2))
        np.savetxt(path + 'imag_fin_dens2.csv',im_dens2,delimiter=',',fmt='%18.16f')
        
        im_dens_tot = np.column_stack((r,np.abs(psi1_im)**2) + np.abs(psi2_im)**2)
        np.savetxt(path + 'imag_fin_dens_tot.csv',im_dens_tot,delimiter=',',fmt='%18.16f')

        E_tot_im = np.column_stack((t_array_im,E1_array_im + E2_array_im))
        np.savetxt(path + 'tot_energy_imag.csv',E_tot_im,delimiter=',',fmt='%18.16f')
        
        [T_IM,R] = np.meshgrid(t_array_im,r)
        plt.close()
        
        plt.pcolormesh(T_IM,R,np.abs(spacetime1_im)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'imag_tot_dens1.png',dpi='figure')
        plt.close()    
        
        plt.pcolormesh(T_IM,R,np.abs(spacetime2_im)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'imag_tot_dens2.png',dpi='figure')
        plt.close()    
        
        plt.pcolormesh(T_IM,R,np.abs(spacetime1_im)**2 + np.abs(spacetime2_im)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'imag_tot_dens_tot.png',dpi='figure')
        plt.close()    
        
    elif setup['IM_T_STEPS'] == 0 and setup['RE_T_STEPS'] > 0:    
        t_array_re = mix_data['t_array_re']
        psi1_re = mix_data['psi1_re']
        psi2_re = mix_data['psi2_re']
        spacetime1_re = mix_data['spacetime1_re']
        spacetime2_re = mix_data['spacetime2_re']

        # output files
        re_dens1 = np.column_stack((r,np.abs(psi1_re)**2))
        np.savetxt(path + 'real_fin_dens1.csv',re_dens,delimiter=',',fmt='%18.16f')
        
        re_dens2 = np.column_stack((r,np.abs(psi2_re)**2))
        np.savetxt(path + 'real_fin_dens1.csv',re_dens,delimiter=',',fmt='%18.16f')

        re_dens_tot = np.column_stack((r,np.abs(psi1_re)**2) + np.abs(psi2_re)**2)
        np.savetxt(path + 'real_fin_dens_tot.csv',re_dens_tot,delimiter=',',fmt='%18.16f')

        [T_RE,R] = np.meshgrid(t_array_re,r)
        
        plt.pcolormesh(T_RE,R,np.abs(spacetime1_re)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'real_tot_dens1.png',dpi='figure')
        plt.close()    
        
        plt.pcolormesh(T_RE,R,np.abs(spacetime2_re)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'real_tot_dens2.png',dpi='figure')
        plt.close()    
        
        plt.pcolormesh(T_RE,R,np.abs(spacetime1_re)**2 + np.abs(spacetime2_re)**2,shading='gouraud')
        plt.xlabel('t')
        plt.ylabel('r')
        plt.savefig(path + 'real_tot_dens2.png',dpi='figure')
        plt.close()    
