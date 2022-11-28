## Import packages
import numpy as np
from numpy import pi
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

def run_ulck_process(dirarg,num_sims):
    # generating empty arrays for saving fitted parameters and confidence intervals
    acs_array = np.empty((num_sims,1))
    mu_cs_array = np.empty((num_sims,1))
    mu_yb_array = np.empty((num_sims,1))
    width_cs_array = np.empty((num_sims,1))
    width_yb_array = np.empty((num_sims,1))
    
    # file to save print outputs to
    sys.stdout=open('../data/' + dirarg + 'output/' +'process.out',"w")
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
        Nr = setup['Nr']

        # read in theory parameters
        f = open('../data/' + dirarg + str(i + 1) + '/theory_params.json',"r")
        theory = json.loads(f.read())
        f.close()

        acs_array[i] = setup['a11']

        mu_cs_array[i] = theory['mu1']
        mu_yb_array[i] = theory['mu2']
        
        # reading in ground states
        r_dens1 = np.loadtxt('../data/' + dirarg + str(i + 1) +  '/imag_fin_dens1.csv', delimiter=",")
        dens1 = r_dens1[:,1]
        r = r_dens1[:,0]
        r_dens2 = np.loadtxt('../data/' + dirarg + str(i + 1) + '/imag_fin_dens2.csv', delimiter=",")
        dens2 = r_dens2[:,1]

        if r[dens1<(0.01*np.max(dens1))].size == 0:
            width_cs_array[i] = Nr*dr
        elif r[dens1<(0.01*np.max(dens1))].size > 0:
            width_cs_array[i] = r[dens1<(0.01*np.max(dens1))][0]
        if r[dens2<(0.01*np.max(dens2))].size == 0:
            width_yb_array[i] = Nr*dr
        elif r[dens2<(0.01*np.max(dens2))].size > 0:
            width_yb_array[i] = r[dens2<(0.01*np.max(dens2))][0]
 
    # concatenating arrays of mu_cs, mu_yb, width_cs, width_yb
    mu_cs_data = np.column_stack((acs_array,mu_cs_array))
    mu_yb_data = np.column_stack((acs_array,mu_yb_array))
    width_cs_data = np.column_stack((acs_array,width_cs_array))
    width_yb_data = np.column_stack((acs_array,width_yb_array))
    
    # saving arrays of mu_cs, mu_yb, width_cs, width_yb
    np.savetxt('../data/' + dirarg + 'output/' + 'cs_mu.csv',mu_cs_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'output/' + 'yb_mu.csv',mu_yb_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'output/' + 'cs_width.csv',width_cs_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'output/' + 'yb_width.csv',width_yb_data,delimiter=',',fmt='%18.16f')

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
    args = parser.parse_args()
    run_ulck_process(args.READ_PATH[0],args.NUM_SIMS[0])
