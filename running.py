from mix_sph_lhy import time

import matplotlib.pyplot as plt
import numpy as np
import json
import os
import sys
import shutil

# select data file to load in
fname = 'config_dens_lck.json'

# give name to directory storing simulation data
dirarg = 'pert_lck_absorb'

# create directory for data storage
path = os.path.join('./data',dirarg+'/')
if not os.path.isdir(path):
    os.mkdir(path)

# file to save print outputs to
sys.stdout=open(path+'sim.out',"w")

# class to reflush print commands 
class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def writelines(self, datas):
       self.stream.writelines(datas)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)

sys.stdout = Unbuffered(sys.stdout)

shutil.copy2(fname,path)

# load in data for simulation
f = open(fname,"r")
setup = json.loads(f.read())
f.close()

# run simulation
mix_data = time(fname)

# extracting data from dictionary
r = mix_data['r']
# save spacetime in a format that it is easy to reload into Python
with open(path + 'r_array.npy', 'wb') as f:
    np.save(f,r)


# data-saving of density-lock mixture
if setup['DENS_LCK'] == 1:
    # both imaginary and real time simulation
    if setup['IM_T_STEPS'] > 0 and setup['RE_T_STEPS'] > 0:
        # load in data from dictionary
        t_array_im = mix_data['t_array_im']
        t_array_re = mix_data['t_array_re']
        phi_im = mix_data['phi_im']
        phi_re = mix_data['phi_re']
        mu_im = mix_data['mu_im']
        E_array_im = mix_data['E_array_im']
        spacetime_im = mix_data['spacetime_im']
        spacetime_re = mix_data['spacetime_re']
        
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'imag_t_array.npy', 'wb') as f:
            np.save(f,t_array_im)
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'real_t_array.npy', 'wb') as f:
            np.save(f,t_array_re)

        # imaginary time final density
        plt.plot(r,np.abs(phi_im)**2)
        plt.xlabel(r'$r$')
        plt.ylabel(r'$|\phi|^2$')
        plt.xlim((0,setup['Nr']*setup['dr']))
        plt.ylim((0,1.2*np.max(np.abs(phi_im)**2)))
        plt.savefig(path + 'imag_fin_dens.png',dpi='figure')
        plt.close()    
        # save final imaginary time density in csv file
        im_dens = np.column_stack((r,np.abs(phi_im)**2))
        np.savetxt(path + 'imag_fin_dens.csv',im_dens,delimiter=',',fmt='%18.16f')
        # real time final density
        plt.plot(r,np.abs(phi_im)**2)
        plt.xlabel(r'$r$')
        plt.ylabel(r'$|\phi|^2$')
        plt.xlim((0,setup['Nr']*setup['dr']))
        plt.ylim((0,1.2*np.max(np.abs(phi_re)**2)))
        plt.savefig(path + 'real_fin_dens.png',dpi='figure')
        plt.close()    
        # save final real time density in csv file
        re_dens = np.column_stack((r,np.abs(phi_re)**2))
        np.savetxt(path + 'real_fin_dens.csv',re_dens,delimiter=',',fmt='%18.16f')
        # total energy in imaginary time 
        E_im = np.column_stack((t_array_im,E_array_im))
        np.savetxt(path + 'energy_imag.csv',E_im,delimiter=',',fmt='%18.16f')
        
        # meshgrid imaginary time and spatial arrays
        [T_IM,R] = np.meshgrid(t_array_im,r)
        # save total image of snapshots from imaginary time 
        plt.pcolormesh(T_IM,R,np.abs(spacetime_im)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\phi|^2$')
        plt.savefig(path + 'imag_spacetime_dens.png',dpi='figure')
        plt.close()    
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'imag_spacetime_wav.npy', 'wb') as f:
            np.save(f, spacetime_im)
        # meshgrid real time and spatial arrays
        [T_RE,R] = np.meshgrid(t_array_re,r)
        # save total image of snapshots from real time
        plt.pcolormesh(T_RE,R,np.abs(spacetime_re)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\phi|^2$')
        plt.savefig(path + 'real_spacetime_dens.png',dpi='figure')
        plt.close()    
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'real_spacetime_wav.npy', 'wb') as f:
            np.save(f, spacetime_re)

    # only imaginary time simulation
    elif setup['IM_T_STEPS'] > 0 and setup['RE_T_STEPS'] == 0:
        # load in data from dictionary
        t_array_im = mix_data['t_array_im']
        phi_im = mix_data['phi_im']
        mu_im = mix_data['mu_im']
        E_array_im = mix_data['E_array_im']
        spacetime_im = mix_data['spacetime_im']
        
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'imag_t_array.npy', 'wb') as f:
            np.save(f,t_array_im)
        
        # imaginary time final density
        plt.plot(r,np.abs(phi_im)**2)
        plt.xlabel(r'$r$')
        plt.ylabel(r'$|\phi|^2$')
        plt.xlim((0,setup['Nr']*setup['dr']))
        plt.ylim((0,1.2*np.max(np.abs(phi_im)**2)))
        plt.savefig(path + 'imag_fin_dens.png',dpi='figure')
        plt.close()    
        # save final imaginary time density in csv file
        im_dens = np.column_stack((r,np.abs(phi_im)**2))
        np.savetxt(path + 'imag_fin_dens.csv',im_dens,delimiter=',',fmt='%18.16f')
        # total energy in imaginary time 
        E_im = np.column_stack((t_array_im,E_array_im))
        np.savetxt(path + 'energy_imag.csv',E_im,delimiter=',',fmt='%18.16f')

        # meshgrid imaginary time and spatial arrays
        [T_IM,R] = np.meshgrid(t_array_im,r)
        # save total image of snapshots from imaginary time 
        plt.pcolormesh(T_IM,R,np.abs(spacetime_im)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\phi|^2$')
        plt.savefig(path + 'imag_spacetime_dens.png',dpi='figure')
        plt.close()    
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'imag_spacetime_wav.npy', 'wb') as f:
            np.save(f, spacetime_im)
    
    # only real time simulation
    elif setup['IM_T_STEPS'] == 0 and setup['RE_T_STEPS'] > 0:
        # load in data from dictionary
        t_array_re = mix_data['t_array_re']
        phi_re = mix_data['phi_re']
        spacetime_re = mix_data['spacetime_re']
        
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'real_t_array.npy', 'wb') as f:
            np.save(f,t_array_re)
        
        # real time final density
        plt.plot(r,np.abs(phi_re)**2)
        plt.xlabel(r'$r$')
        plt.ylabel(r'$|\phi|^2$')
        plt.xlim((0,setup['Nr']*setup['dr']))
        plt.ylim((0,1.2*np.max(np.abs(phi_re)**2)))
        plt.savefig(path + 'real_fin_dens.png',dpi='figure')
        plt.close()   
        # save final real time density in csv file
        re_dens = np.column_stack((r,np.abs(phi_re)**2))
        np.savetxt(path + 'real_fin_dens.csv',re_dens,delimiter=',',fmt='%18.16f')

        # meshgrid real time and spatial arrays
        [T_RE,R] = np.meshgrid(t_array_re,r)
        # save total image of snapshots from real time 
        plt.pcolormesh(T_RE,R,np.abs(spacetime_re)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)	
        cbar = plt.colorbar()
        cbar.set_label(r'$|\phi|^2$')
        plt.savefig(path + 'real_spacetime_dens.png',dpi='figure')
        plt.close()    
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'real_spacetime_wav.npy', 'wb') as f:
            np.save(f, spacetime_re)

# data saving of density-unlocked mixture
elif setup['DENS_LCK'] == 0:
    # both imaginary and real time simulation
    if setup['IM_T_STEPS'] > 0 and setup['RE_T_STEPS'] > 0: 
        # load in data from dictionary
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
        
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'imag_t_array.npy', 'wb') as f:
            np.save(f,t_array_im)
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'real_t_array.npy', 'wb') as f:
            np.save(f,t_array_re)
        
        # imaginary time final density of component 1
        im_dens1 = np.column_stack((r,np.abs(psi1_im)**2))
        np.savetxt(path + 'imag_fin_dens1.csv',im_dens1,delimiter=',',fmt='%18.16f')
        # imaginary time final density of component 2
        im_dens2 = np.column_stack((r,np.abs(psi2_im)**2))
        np.savetxt(path + 'imag_fin_dens2.csv',im_dens2,delimiter=',',fmt='%18.16f')
        # imaginary time final total density 
        im_dens_tot = np.column_stack((r,np.abs(psi1_im)**2) + np.abs(psi2_im)**2)
        np.savetxt(path + 'imag_fin_dens_tot.csv',im_dens_tot,delimiter=',',fmt='%18.16f')
        # imaginary time final density
        plt.plot(r,np.abs(psi1_im)**2,r,np.abs(psi2_im)**2,r,np.abs(psi1_im)**2 + np.abs(psi2_im)**2)
        plt.xlabel(r'$r$')
        plt.ylabel(r'$|\psi|^2$')
        plt.xlim((0,setup['Nr']*setup['dr']))
        plt.ylim((0,1.2*np.max(np.abs(psi1_im)**2 + np.abs(psi2_im)**2)))
        plt.legend((r'$|\psi_1|^2$',r'$|\psi_2|^2$',r'$|\psi_{1+2}|^2$'))
        plt.savefig(path + 'imag_fin_dens.png',dpi='figure')
        plt.close()   
        # real time final density of component 1
        re_dens1 = np.column_stack((r,np.abs(psi1_re)**2))
        np.savetxt(path + 'real_fin_dens1.csv',re_dens1,delimiter=',',fmt='%18.16f')
        # real time final density of component 2
        re_dens2 = np.column_stack((r,np.abs(psi2_re)**2))
        np.savetxt(path + 'real_fin_dens2.csv',re_dens2,delimiter=',',fmt='%18.16f')
        # real time final total density
        re_dens_tot = np.column_stack((r,np.abs(psi1_re)**2) + np.abs(psi2_re)**2)
        np.savetxt(path + 'real_fin_dens_tot.csv',re_dens_tot,delimiter=',',fmt='%18.16f')
        # real time final density
        plt.plot(r,np.abs(psi1_re)**2,r,np.abs(psi2_re)**2,r,np.abs(psi1_re)**2 + np.abs(psi2_re)**2)
        plt.xlabel(r'$r$')
        plt.ylabel(r'$|\psi|^2$')
        plt.xlim((0,setup['Nr']*setup['dr']))
        plt.ylim((0,1.2*np.max(np.abs(psi1_re)**2 + np.abs(psi2_re)**2)))
        plt.legend((r'$|\psi_1|^2$',r'$|\psi_2|^2$',r'$|\psi_{1+2}|^2$'))
        plt.savefig(path + 'real_fin_dens.png',dpi='figure')
        plt.close()
        # total energy in imaginary time
        E_tot_im = np.column_stack((t_array_im,E1_array_im + E2_array_im))
        np.savetxt(path + 'tot_energy_imag.csv',E_tot_im,delimiter=',',fmt='%18.16f')
        
        # meshgrid imaginary time and spatial arrays
        [T_IM,R] = np.meshgrid(t_array_im,r)
        # save image of component 1 density snapshots from imaginary time 
        plt.pcolormesh(T_IM,R,np.abs(spacetime1_im)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)		
        cbar = plt.colorbar()
        cbar.set_label(r'$|\psi_1|^2$')
        plt.savefig(path + 'imag_spacetime_dens1.png',dpi='figure')
        plt.close()    
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'imag_spacetime_wav1.npy', 'wb') as f:
            np.save(f, spacetime1_im)
        # save image of component 2 density snapshots from imaginary time 
        plt.pcolormesh(T_IM,R,np.abs(spacetime2_im)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\psi_2|^2$')
        plt.savefig(path + 'imag_spacetime_dens2.png',dpi='figure')
        cbar = plt.colorbar()
        plt.close()    
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'imag_spacetime_wav2.npy', 'wb') as f:
            np.save(f, spacetime2_im)
        # save image of total density snapshots from imaginary time 
        plt.pcolormesh(T_IM,R,np.abs(spacetime1_im)**2 + np.abs(spacetime2_im)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\psi_{1+2}|^2$')
        plt.savefig(path + 'imag_spacetime_dens_tot.png',dpi='figure')
        plt.close()    
        
        # meshgrid real time and spatial arrays
        [T_RE,R] = np.meshgrid(t_array_re,r)
        # save image of component 1 density snapshots from real time 
        plt.pcolormesh(T_RE,R,np.abs(spacetime1_re)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\psi_1|^2$')
        plt.savefig(path + 'real_spacetime_dens1.png',dpi='figure')
        plt.close()    
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'real_spacetime_wav1.npy', 'wb') as f:
            np.save(f, spacetime1_re)
        # save image of component 2 density snapshots from real time 
        plt.pcolormesh(T_RE,R,np.abs(spacetime2_re)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\psi_2|^2$')
        plt.savefig(path + 'real_spacetime_dens2.png',dpi='figure')
        plt.close()    
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'real_spacetime_wav2.npy', 'wb') as f:
            np.save(f, spacetime2_re)
        # save image of total density snapshots from real time 
        plt.pcolormesh(T_RE,R,np.abs(spacetime1_re)**2 + np.abs(spacetime2_re)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\psi_{1+2}|^2$')
        plt.savefig(path + 'real_spacetime_dens_tot.png',dpi='figure')
        plt.close()    
    
    # only imaginary time simulation
    elif setup['IM_T_STEPS'] > 0 and setup['RE_T_STEPS'] == 0: 
        # load in data from dictionary
        t_array_im = mix_data['t_array_im']
        psi1_im = mix_data['psi1_im']
        psi2_im = mix_data['psi2_im']
        mu1_im = mix_data['mu1_im']
        mu2_im = mix_data['mu2_im']
        E1_array_im = mix_data['E1_array_im']
        E2_array_im = mix_data['E2_array_im']
        spacetime1_im = mix_data['spacetime1_im']
        spacetime2_im = mix_data['spacetime2_im']
       
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'imag_t_array.npy', 'wb') as f:
            np.save(f,t_array_im)
        
        # imaginary time final density of component 1
        im_dens1 = np.column_stack((r,np.abs(psi1_im)**2))
        np.savetxt(path + 'imag_fin_dens1.csv',im_dens1,delimiter=',',fmt='%18.16f')
        # imaginary time final density of component 2
        im_dens2 = np.column_stack((r,np.abs(psi2_im)**2))
        np.savetxt(path + 'imag_fin_dens2.csv',im_dens2,delimiter=',',fmt='%18.16f')
        # imaginary time final total density
        im_dens_tot = np.column_stack((r,np.abs(psi1_im)**2 + np.abs(psi2_im)**2))
        np.savetxt(path + 'imag_fin_dens_tot.csv',im_dens_tot,delimiter=',',fmt='%18.16f')
        # imaginary time final density
        plt.plot(r,np.abs(psi1_im)**2,r,np.abs(psi2_im)**2,r,np.abs(psi1_im)**2 + np.abs(psi2_im)**2)
        plt.xlabel(r'$r$')
        plt.ylabel(r'$|\psi|^2$')
        plt.xlim((0,setup['Nr']*setup['dr']))
        plt.ylim((0,1.2*np.max(np.abs(psi1_im)**2 + np.abs(psi2_im)**2)))
        plt.legend((r'$|\psi_1|^2$',r'$|\psi_2|^2$',r'$|\psi_{1+2}|^2$'))
        plt.savefig(path + 'imag_fin_dens.png',dpi='figure')
        plt.close()   
        # total energy in imaginary time
        E_tot_im = np.column_stack((t_array_im,E1_array_im + E2_array_im))
        np.savetxt(path + 'tot_energy_imag.csv',E_tot_im,delimiter=',',fmt='%18.16f') 
        plt.close()

        # meshgrid imaginary time and spatial arrays
        [T_IM,R] = np.meshgrid(t_array_im,r)
        # save image of component 1 density snapshots from imaginary time 
        plt.pcolormesh(T_IM,R,np.abs(spacetime1_im)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\psi_1|^2$')
        plt.savefig(path + 'imag_spacetime_dens1.png',dpi='figure')
        plt.close()    
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'imag_spacetime_wav1.npy', 'wb') as f:
            np.save(f, spacetime1_im)
        # save image of component 2 density snapshots from imaginary time 
        plt.pcolormesh(T_IM,R,np.abs(spacetime2_im)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\psi_2|^2$')
        plt.savefig(path + 'imag_spacetime_dens2.png',dpi='figure')
        plt.close()    
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'imag_spacetime_wav2.npy', 'wb') as f:
            np.save(f, spacetime2_im)
        # save image of total density snapshots from imaginary time 
        plt.pcolormesh(T_IM,R,np.abs(spacetime1_im)**2 + np.abs(spacetime2_im)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\psi_{1+2}|^2$')
        plt.savefig(path + 'imag_spacetime_dens_tot.png',dpi='figure')
        plt.close()    
        
    # only real time simulation
    elif setup['IM_T_STEPS'] == 0 and setup['RE_T_STEPS'] > 0:    
        # load in data from dictionary
        t_array_re = mix_data['t_array_re']
        psi1_re = mix_data['psi1_re']
        psi2_re = mix_data['psi2_re']
        spacetime1_re = mix_data['spacetime1_re']
        spacetime2_re = mix_data['spacetime2_re']

        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'real_t_array.npy', 'wb') as f:
            np.save(f,t_array_re)
        
        # real time final density of component 1
        re_dens1 = np.column_stack((r,np.abs(psi1_re)**2))
        np.savetxt(path + 'real_fin_dens1.csv',re_dens1,delimiter=',',fmt='%18.16f')
        # real time final density of component 2
        re_dens2 = np.column_stack((r,np.abs(psi2_re)**2))
        np.savetxt(path + 'real_fin_dens2.csv',re_dens2,delimiter=',',fmt='%18.16f')
        # real time final total density
        re_dens_tot = np.column_stack((r,np.abs(psi1_re)**2) + np.abs(psi2_re)**2)
        np.savetxt(path + 'real_fin_dens_tot.csv',re_dens_tot,delimiter=',',fmt='%18.16f')
        # real time final density
        plt.plot(r,np.abs(psi1_re)**2,r,np.abs(psi2_re)**2,r,np.abs(psi1_re)**2 + np.abs(psi2_re)**2)
        plt.xlabel(r'$r$')
        plt.ylabel(r'$|\psi|^2$')
        plt.xlim((0,setup['Nr']*setup['dr']))
        plt.ylim((0,1.2*np.max(np.abs(psi1_re)**2 + np.abs(psi2_re)**2)))
        plt.legend((r'$|\psi_1|^2$',r'$|\psi_2|^2$',r'$|\psi_{1+2}|^2$'))
        plt.savefig(path + 'real_fin_dens.png',dpi='figure')
        plt.close()

        # meshgrid real time and spatial arrays
        [T_RE,R] = np.meshgrid(t_array_re,r)
        # save image of component 1 density snapshots from real time 
        plt.pcolormesh(T_RE,R,np.abs(spacetime1_re)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\psi_1|^2$')
        plt.savefig(path + 'real_spacetime_dens1.png',dpi='figure')
        plt.close()    
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'real_spacetime_wav1.npy', 'wb') as f:
            np.save(f, spacetime1_re)
        # save image of component 2 density snapshots from real time 
        plt.pcolormesh(T_RE,R,np.abs(spacetime2_re)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\psi_2|^2$')
        plt.savefig(path + 'real_spacetime_dens2.png',dpi='figure')
        plt.close()    
        # save spacetime in a format that it is easy to reload into Python
        with open(path + 'real_spacetime_wav2.npy', 'wb') as f:
            np.save(f, spacetime2_re)
        # save image of total density snapshots from real time 
        plt.pcolormesh(T_RE,R,np.abs(spacetime1_re)**2 + np.abs(spacetime2_re)**2,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.ylabel(r'$r$')
        plt.clim(0.0,None)
        cbar = plt.colorbar()
        cbar.set_label(r'$|\psi_{1+2}|^2$')
        plt.savefig(path + 'real_spacetime_dens_tot.png',dpi='figure')
        plt.close()    

sys.stdout.close()
