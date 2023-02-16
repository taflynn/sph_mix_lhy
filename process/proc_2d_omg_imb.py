## Import packages
import numpy as np
import json
import os
import sys
import shutil
import argparse

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

def run_ulck_process(dirarg, omg_num, imb_num):

    # generating empty arrays saving extracted data
    centre_dens1_out_array = np.empty((imb_num, omg_num))
    centre_dens2_out_array = np.empty((imb_num, omg_num))

    end_dens1_out_array = np.empty((imb_num, omg_num))
    end_dens2_out_array = np.empty((imb_num, omg_num))

    energy_out_array = np.empty((imb_num, omg_num))

    mu1_out_array = np.empty((imb_num, omg_num))
    mu2_out_array = np.empty((imb_num, omg_num))

    omg_array = np.empty((omg_num))
    imb_array = np.empty((imb_num))
    
    # file to save print outputs to
    sys.stdout=open('../data/' + dirarg + 'saved/' +'process.out',"w")
    sys.stdout = Unbuffered(sys.stdout)

    job = 1    

    # iterate across 2D array
    for j in range(0, omg_num):
        for i in range(0, imb_num):
            # load in simulation parameters
            fname = 'config_dens_ulck' + str(job) + '.json'
            print(i)

            # read in data
            f = open('../data/' + dirarg + str(job)  + '/' + fname, "r")
            setup = json.loads(f.read())
            f.close()
      
            dr = setup['dr']

            # read in theory parameters
            f = open('../data/' + dirarg + str(job) + '/theory_params.json', "r")
            theory = json.loads(f.read())
            f.close()
        
            # reading in ground states
            r_dens1 = np.loadtxt('../data/' + dirarg + str(job) +  '/imag_fin_dens1.csv', delimiter=",")
            dens1 = r_dens1[:,1]
            r_dens2 = np.loadtxt('../data/' + dirarg + str(job) + '/imag_fin_dens2.csv', delimiter=",")
            dens2 = r_dens2[:,1]
            t_eng = np.loadtxt('../data/' + dirarg + str(job) + '/tot_energy_imag.csv', delimiter=",")

            centre_dens1_out_array[i, j] = dens1[1]
            centre_dens2_out_array[i, j] = dens2[1]        
            end_dens1_out_array[i, j] = dens1[-2]    
            end_dens2_out_array[i, j] = dens2[-2]    
            energy_out_array[i, j] = t_eng[-1,1]
            mu1_out_array[i, j] = theory['mu1']
            mu2_out_array[i, j] = theory['mu2']
           
            if j == 0:
              imb_array[i] = ((setup['N1'] - setup['N2'])/setup['N2'])*100.0
            job += 1
        omg_array[j] = setup['OMEGA1']
 
    # saving arrays of omega's and confidence intervals
    np.savetxt('../data/' + dirarg + 'saved/' + 'centre_dens1.csv', centre_dens1_out_array, delimiter=',', fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'centre_dens2.csv', centre_dens2_out_array, delimiter=',', fmt='%18.16f')

    np.savetxt('../data/' + dirarg + 'saved/' + 'end_dens1.csv', end_dens1_out_array, delimiter=',', fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'end_dens2.csv', end_dens2_out_array, delimiter=',', fmt='%18.16f')

    np.savetxt('../data/' + dirarg + 'saved/' + 'energy.csv', energy_out_array, delimiter=',', fmt='%18.16f')

    np.savetxt('../data/' + dirarg + 'saved/' + 'mu1.csv', mu1_out_array, delimiter=',', fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'mu2.csv', mu2_out_array, delimiter=',', fmt='%18.16f')

    np.savetxt('../data/' + dirarg + 'saved/' + 'omegas.csv', omg_array, delimiter=',', fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'imbalances.csv', imb_array, delimiter=',', fmt='%18.16f')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Process data across 2D parameter space of trapping freq. and imbalance')
    parser.add_argument('--read_path', '-rp',
            dest = 'READ_PATH',
            type = str,
            required = True,
            nargs = 1)
    parser.add_argument('--omg_num', '-on',
            dest = 'NUM_OF_OMEGAS',
            type = int,
            required = True,
            nargs = 1)
    parser.add_argument('--imb_num', '-in',
            dest = 'NUM_OF_IMBALANCES',
            type = int,
            required = True,
            nargs = 1)
    args = parser.parse_args()
    run_ulck_process(args.READ_PATH[0], args.NUM_OF_OMEGAS[0], args.NUM_OF_IMBALANCES[0])
