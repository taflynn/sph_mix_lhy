from uneqm_sph_lhy import time

import matplotlib.pyplot as plt
from numpy import save, absolute, amax, column_stack, savetxt, maximum
import json
import os
import sys
import shutil
import argparse

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

def main(dirarg,fname):
    # create directory for data storage
    path = os.path.join('./data',dirarg+'/')
    if not os.path.isdir(path):
        os.mkdir(path)

    # file to save print outputs to
    sys.stdout=open(path+'sim.out',"w")
    sys.stdout = Unbuffered(sys.stdout)
    shutil.copy2(fname,path)

    # load in data for simulation
    f = open(fname,"r")
    setup = json.loads(f.read())
    f.close()

    # run simulation
    [mix_data,theory_params] = time(fname,path)

    # save theoretical parameters dictionary to a json file
    with open(path + 'theory_params.json', 'w') as f:
        json.dump(theory_params, f, indent=0)

    # extracting data from dictionary
    r = mix_data['r']
    # save spacetime in a format that it is easy to reload into Python
    with open(path + 'r_array.npy', 'wb') as f:
        save(f,r)
    
    # data saving of density-unlocked mixture
    if setup['DENS_LCK'] == 0:
        # imaginary time simulation
        if setup['IM_T_STEPS'] > 0: 
            # load in data from dictionary
            t_array_im = mix_data['t_array_im']
            psi1_im = mix_data['psi1_im']
            psi2_im = mix_data['psi2_im']
            mu1_im = mix_data['mu1_im']
            mu2_im = mix_data['mu2_im']
            E_array_im = mix_data['E_array_im']

            # save spacetime in a format that it is easy to reload into Python
            with open(path + 'imag_t_array.npy', 'wb') as f:
                save(f,t_array_im)

            # imaginary time final density of component 1
            im_dens1 = column_stack((r,absolute(psi1_im)**2))
            savetxt(path + 'imag_fin_dens1.csv',im_dens1,delimiter=',',fmt='%18.16f')
            # imaginary time final density of component 2
            im_dens2 = column_stack((r,absolute(psi2_im)**2))
            savetxt(path + 'imag_fin_dens2.csv',im_dens2,delimiter=',',fmt='%18.16f')
            # imaginary time final density
            plt.plot(r,absolute(psi1_im)**2,r,absolute(psi2_im)**2)
            plt.xlabel(r'$r$')
            plt.ylabel(r'$|\psi|^2$')
            plt.xlim((0,setup['Nr']*setup['dr']))
            plt.ylim((0,1.2*maximum(absolute(psi1_im[1])**2,absolute(psi2_im[1])**2)))
            plt.legend((r'$|\psi_1|^2$',r'$|\psi_2|^2$'))
            plt.savefig(path + 'imag_fin_dens.png',dpi='figure')
            plt.close()   
            # total energy in imaginary time
            E_tot_im = column_stack((t_array_im,E_array_im))
            savetxt(path + 'tot_energy_imag.csv',E_tot_im,delimiter=',',fmt='%18.16f') 

        # real time simulation
        elif setup['RE_T_STEPS'] > 0:    
            # load in data from dictionary
            t_array_re = mix_data['t_array_re']
            psi1_re = mix_data['psi1_re']
            psi2_re = mix_data['psi2_re']

            # save spacetime in a format that it is easy to reload into Python
            with open(path + 'real_t_array.npy', 'wb') as f:
                save(f,t_array_re)

            # real time final density of component 1
            re_dens1 = column_stack((r,absolute(psi1_re)**2))
            savetxt(path + 'real_fin_dens1.csv',re_dens1,delimiter=',',fmt='%18.16f')
            # real time final density of component 2
            re_dens2 = column_stack((r,absolute(psi2_re)**2))
            savetxt(path + 'real_fin_dens2.csv',re_dens2,delimiter=',',fmt='%18.16f')
            # real time final density
            plt.plot(r,absolute(psi1_re)**2,r,absolute(psi2_re)**2)
            plt.xlabel(r'$r$')
            plt.ylabel(r'$|\psi|^2$')
            plt.xlim((0,setup['Nr']*setup['dr']))
            plt.ylim((0,1.2*maximum(absolute(psi1_re[1])**2,absolute(psi2_re[1])**2)))
            plt.legend((r'$|\psi_1|^2$',r'$|\psi_2|^2$'))
            plt.savefig(path + 'real_fin_dens.png',dpi='figure')
            plt.close()

    sys.stdout.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Running mixture (w/ LHY) simulation and outputting data')
    parser.add_argument('--write_path','-wp',
            dest = 'WRITE_PATH',
            type = str,
            required = True,
            nargs = 1)
    parser.add_argument('--read_path','-rp',
            dest = 'READ_PATH',
            type = str,
            required = True,
            nargs = 1)
    args = parser.parse_args()
    main(args.WRITE_PATH[0],args.READ_PATH[0])
   
