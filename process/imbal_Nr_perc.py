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

def run_ulck_process(dirarg,num_sims,imbal_size):
    if imbal_size == 'Nr': 
        sim_style = 'box_size'
    elif imbal_size == 'perc':
        sim_style = 'imbal'
    else:
        print('Error: Need to specify style of simulation, i.e., imbalance or droplet size')
        sim_style = 'N/A'

    # generating empty arrays for saving fitted parameters and confidence intervals
    del_dens_out_array = np.empty((num_sims,1))
    energy_out_array = np.empty((num_sims,1))
    Nr_or_perc_array = np.empty((num_sims,1))
    
    # file to save print outputs to
    sys.stdout=open('../data/' + dirarg + 'saved/' +'process.out',"w")
    sys.stdout = Unbuffered(sys.stdout)

    for i in range(0,num_sims):
        # load in simulation parameters
        fname = 'config_dens_ulck' + str(i + 1) + '.json'
        print(i)

        # read in data
        f = open('../data/' + dirarg + str(i + 1)  + '/' + fname,"r")
        setup = json.loads(f.read())
        f.close()
      
        dr = setup['dr']

        # read in theory parameters
        f = open('../data/' + dirarg + str(i + 1) + '/theory_params.json',"r")
        theory = json.loads(f.read())
        f.close()

        if imbal_size == 'perc':
            Nr_or_perc_array[i] = ((setup['N1']-setup['N2'])/setup['N2'])*100.0
        elif imbal_size == 'Nr':
            Nr_or_perc_array[i] = 1/(dr*setup['Nr'])
        
        # reading in ground states
        r_dens1 = np.loadtxt('../data/' + dirarg + str(i + 1) +  '/imag_fin_dens1.csv', delimiter=",")
        dens1 = r_dens1[:,1]
        r_dens2 = np.loadtxt('../data/' + dirarg + str(i + 1) + '/imag_fin_dens2.csv', delimiter=",")
        dens2 = r_dens2[:,1]
        t_eng = np.loadtxt('../data/' + dirarg + str(i + 1) + '/tot_energy_imag.csv', delimiter=",")

        del_dens_out_array[i] = dens1[1] - dens2[1]        
        energy_out_array[i] = t_eng[-1,1]
 
    # concatenating data into 2D array of omega and gamma for component 1 
    dens_data = np.column_stack((Nr_or_perc_array,del_dens_out_array))
    eng_data = np.column_stack((Nr_or_perc_array,energy_out_array))
    
    # saving arrays of omega's and confidence intervals
    np.savetxt('../data/' + dirarg + 'saved/' + 'del_dens_' + sim_style + '.csv',dens_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'energy_' + sim_style + '.csv',eng_data,delimiter=',',fmt='%18.16f')

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
    args = parser.parse_args()
    run_ulck_process(args.READ_PATH[0],args.NUM_SIMS[0],args.IMBAL_SIZE[0])
