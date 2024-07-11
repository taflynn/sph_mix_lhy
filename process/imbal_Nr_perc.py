## Import packages
import numpy as np
import json
import os
import sys
import shutil
import argparse
from scipy.optimize import curve_fit

sys.path.insert(1, '../')

from main.params_calc import params_dens_lck
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

def run_ulck_process(dirarg, num_sims, imbal_size, which_maj):
    if imbal_size == 'Nr': 
        sim_style = 'box_size'
    elif imbal_size == 'perc':
        sim_style = 'imbal'
    else:
        print('Error: Need to specify style of simulation, i.e., imbalance or droplet size')
        sim_style = 'N/A'

    # generating empty arrays for saving fitted parameters and confidence intervals
    ## difference in densities
    del_dens_out_array = np.empty((num_sims, 1))
    ## end densities
    end_dens1_out_array = np.empty((num_sims, 1))
    end_dens2_out_array = np.empty((num_sims, 1))
    near_end_dens1_out_array = np.empty((num_sims, 1))
    near_end_dens2_out_array = np.empty((num_sims, 1))
    ## energies
    energy_out_array = np.empty((num_sims, 1))
    ## chemical potentials
    mu1_out_array = np.empty((num_sims, 1))
    mu2_out_array = np.empty((num_sims, 1))
    ## independent variables
    Nr_or_perc_array = np.empty((num_sims, 1))
    Ns_array = np.empty((num_sims, 1))
    
    # file to save print outputs to
    sys.stdout=open('../data/' + dirarg + 'saved/' +'process.out', "w")
    sys.stdout = Unbuffered(sys.stdout)

    for i in range(0, num_sims):
        # load in simulation parameters
        fname = 'config_dens_ulck' + str(i + 1) + '.json'
        print(i)

        # read in data
        f = open('../data/' + dirarg + str(i + 1)  + '/' + fname, "r")
        setup = json.loads(f.read())
        f.close()
    
        dr = setup['dr']

        # read in theory parameters
        f = open('../data/' + dirarg + str(i + 1) + '/theory_params.json', "r")
        theory = json.loads(f.read())
        f.close()

        if imbal_size == 'perc':
            if which_maj == 1:
                Nr_or_perc_array[i] = ((setup['N1'] - setup['N2'])/setup['N2'])*100.0
                Ns_array[i] = (theory['N1'] - theory['N2'])*10**(-3)
                print('1st-component is majority')
            elif which_maj == 2:
                Nr_or_perc_array[i] = ((setup['N1'] - setup['N2'])/setup['N2'])*100.0
                Ns_array[i] = (theory['N2'] - theory['N1'])*10**(-3)
                print('2nd-component is majority')
        elif imbal_size == 'Nr':
            Nr_or_perc_array[i] = 1/(dr*setup['Nr'])
        
        # reading in ground states
        r_dens1 = np.loadtxt('../data/' + dirarg + str(i + 1) +  '/imag_fin_dens1.csv', delimiter=",")
        r_array = r_dens1[:, 0]
        dens1 = r_dens1[:, 1]
        r_dens2 = np.loadtxt('../data/' + dirarg + str(i + 1) + '/imag_fin_dens2.csv', delimiter=",")
        dens2 = r_dens2[:, 1]
        t_eng = np.loadtxt('../data/' + dirarg + str(i + 1) + '/tot_energy_imag.csv', delimiter=",")

        ## difference in central densities
        del_dens_out_array[i] = dens1[1] - dens2[1]        
        ## end densities
        end_dens1_out_array[i] = dens1[-2] 
        end_dens2_out_array[i] = dens2[-2]    
        ## densities at 0.75Lr
        near_edge = np.where(r_array>=(0.75*dr*setup['Nr']))[0][0]
        near_end_dens1_out_array[i] = dens1[near_edge]    
        near_end_dens2_out_array[i] = dens2[near_edge]    
        ## energies
        energy_out_array[i] = t_eng[-1, 1]
        ## chemical potentials
        mu1_out_array[i] = theory['mu1']
        mu2_out_array[i] = theory['mu2']
 
    # concatenating data into 2D array of omega and gamma for component 1 
    dens_data = np.column_stack((Ns_array, del_dens_out_array))
    dens1_end_data = np.column_stack((Ns_array, end_dens1_out_array))
    dens2_end_data = np.column_stack((Ns_array, end_dens2_out_array))
    eng_data = np.column_stack((Ns_array, energy_out_array))
    mu1_data = np.column_stack((Ns_array, mu1_out_array))
    mu2_data = np.column_stack((Ns_array, mu2_out_array))
    
    # new outputs, as function of N1
    del_dens_v_N = np.column_stack((Ns_array, del_dens_out_array))
    near_end1_v_N = np.column_stack((Ns_array, near_end_dens1_out_array))

    # saving arrays of omega's and confidence intervals
    np.savetxt('../data/' + dirarg + 'saved/' + 'del_dens_' + sim_style + '.csv', dens_data, delimiter=',', fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'end_dens1_' + sim_style + '.csv', dens1_end_data, delimiter=',', fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'end_dens2_' + sim_style + '.csv', dens2_end_data, delimiter=',', fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'energy_' + sim_style + '.csv', eng_data, delimiter=',', fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'mu1_' + sim_style + '.csv', mu1_data, delimiter=',', fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'mu2_' + sim_style + '.csv', mu2_data, delimiter=',', fmt='%18.16f')

    np.savetxt('../data/' + dirarg + 'saved/' + 'N_v_del_dens_' + sim_style + '.csv', del_dens_v_N, delimiter=',', fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'N_v_near_end1_' + sim_style + '.csv', near_end1_v_N, delimiter=',', fmt='%18.16f')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Process data of density-unlocked mixture simulation')
    parser.add_argument('--read_path','-rp',
            dest = 'READ_PATH',
            type = str,
            required = True,
            nargs = 1)
    parser.add_argument('--num_sims','-ns',
            dest = 'NUM_SIMS',
            type = int,
            required = True,
            nargs = 1)
    parser.add_argument('--sim_type','-st',
            dest = 'IMBAL_SIZE',
            type = str,
            required = True,
            nargs = 1)
    parser.add_argument('--which_maj','-wm',
            dest = 'WHICH_COMP_IS_MAJ',
            type = int,
            required = True,
            nargs = 1)
    args = parser.parse_args()
    run_ulck_process(args.READ_PATH[0], args.NUM_SIMS[0], args.IMBAL_SIZE[0], args.WHICH_COMP_IS_MAJ[0])
