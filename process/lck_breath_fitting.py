## Import packages
import numpy as np
import matplotlib.pyplot as plt
import json
import os
import sys
import shutil
import matplotlib.animation as animation
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
# Define a general damped sine function
    def damp_sin_func(x,a,b,c,d,f):
        return a*np.exp(-f*x)*np.sin(b*x + c) + d
    # function to fit the 
    def curve_fitting(t_array,centre_dens):
        # Initial guess parameters for the curve_fit function
        A = np.max(centre_dens) - centre_dens[0]  # initial guess of amplitude
        B = 0.3 # very rough guess for frequency
        C = 0.0 # default input value for phase shift
        D =  centre_dens[0] # initial guess for the y shift of the sine curve
        F = 0.01
    
        # Extracting the fitted parameter values from the curve_fit function
        popt,pcov = curve_fit(damp_sin_func,t_array[0:-1],centre_dens[0:-1],p0=[A,B,C,D,F])
        # Return only the second parameter as this is the breathing mode frequency
        return popt,pcov
    # extract values of density arrays and time arrays after the first max of the oscillating density
    def turning_points(array):
        ''' turning_points(array) -> min_indices, max_indices
        Finds the turning points within an 1D array and returns the indices of the minimum and 
        maximum turning points in two separate lists.
        '''
        idx_max, idx_min = [], []
        if (len(array) < 3): 
            return idx_min, idx_max

        NEUTRAL, RISING, FALLING = range(3)
        def get_state(a, b):
            if a < b: return RISING
            if a > b: return FALLING
            return NEUTRAL

        ps = get_state(array[0], array[1])
        begin = 1
        for i in range(2, len(array)):
            s = get_state(array[i - 1], array[i])
            if s != NEUTRAL:
                if ps != NEUTRAL and ps != s:
                    if s == FALLING: 
                        idx_max.append((begin + i - 1) // 2)
                    else:
                        idx_min.append((begin + i - 1) // 2)
                begin = i
                ps = s
        return idx_min, idx_max
    
    # generating empty arrays for saving fitted parameters and confidence intervals
    omega0_array = np.empty((num_sims,1))
    lowconf_omega0 = np.empty((num_sims,1))
    highconf_omega0 = np.empty((num_sims,1))
    
    gamma_array = np.empty((num_sims,1))
    lowconf_gamma = np.empty((num_sims,1))
    highconf_gamma = np.empty((num_sims,1))
    
    mu_array = np.empty((num_sims,1))
    
    size_array = np.empty((num_sims,1))
    
    # file to save print outputs to
    sys.stdout=open('../data/' + dirarg + 'saved/' +'process.out',"w")
    sys.stdout = Unbuffered(sys.stdout)

    for i in range(0,num_sims):
        # load in simulation parameters
        fname = 'config_dens_lck' + str(i + 1) + '.json'
        print(i)

        # read in data
        f = open('../data/' + dirarg + str(i + 1)  + '/' + fname,"r")
        setup = json.loads(f.read())
        f.close()

        # grid spacing
        dr = setup['dr']

        # read in theory parameters
        f = open('../data/' + dirarg + str(i + 1) + '/theory_params.json',"r")
        theory = json.loads(f.read())
        f.close()
        
        N_lck = theory['N_lck']
        size_array[i] = (N_lck - 18.65)**0.25

        # main data read in from NumPy files
        r = np.load('../data/' + dirarg + str(i+1) + '/r_array.npy')
        t_real = np.load('../data/' + dirarg + str(i+1) + '/real_t_array.npy')
        phi = np.load('../data/' + dirarg + str(i+1) + '/real_spacetime_wav.npy')
        
        # extract central density values
        centre_n = (np.abs(phi[1,:])**2)
        
        # cutoff the first transient of the oscillations 
        #[idx_min,idx_max] = turning_points(centre_n) 
        #if centre_n[idx_max[0]] < centre_n[0]:
        #    extract = idx_max[1]
        #else:
        #    extract = idx_max[0]
        if i == 0:
            extract = 0
        elif i > 0:
            extract = (np.where(centre_n>centre_n[0]))[0][0]
        cut_n = centre_n[extract:]
        cut_t = t_real[extract:] - t_real[extract]

        # fit central density oscillations to damped sine curve
        [fitted_params,cov] = curve_fitting(cut_t,cut_n)
        
        # extract breathing mode frequency and decay rate of oscillation
        omega0_array[i] = np.abs(fitted_params[1])
        lowconf_omega0[i] = np.abs(fitted_params[1]) - 2*np.sqrt(cov[1,1])
        highconf_omega0[i] = np.abs(fitted_params[1]) + 2*np.sqrt(cov[1,1])
        
        gamma_array[i] = fitted_params[-1]
        lowconf_gamma[i] = fitted_params[-1] - 2*np.sqrt(cov[-1,-1])
        highconf_gamma[i] = fitted_params[-1] + 2*np.sqrt(cov[-1,-1])
        
        mu_array[i] = -1*theory['mu']

    # concatenating data into 2D array of omega and gamma 
    omega0_data = np.column_stack((size_array,omega0_array))
    lowconf_omega0_data = np.column_stack((size_array,lowconf_omega0))
    highconf_omega0_data = np.column_stack((size_array,highconf_omega0))

    gamma_data = np.column_stack((size_array,gamma_array))
    lowconf_gamma_data = np.column_stack((size_array,lowconf_gamma))
    highconf_gamma_data = np.column_stack((size_array,highconf_gamma))

    mu_data = np.column_stack((size_array,mu_array))
 
    # saving arrays of omega's and confidence intervals
    np.savetxt('../data/' + dirarg + 'saved/' + 'omega0_size.csv',omega0_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'omega0_low_size.csv',lowconf_omega0_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'omega0_high_size.csv',highconf_omega0_data,delimiter=',',fmt='%18.16f')

    np.savetxt('../data/' + dirarg + 'saved/' + 'gamma_size.csv',gamma_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'gamma_low_size.csv',lowconf_gamma_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'gamma_high_size.csv',highconf_gamma_data,delimiter=',',fmt='%18.16f')
    
    np.savetxt('../data/' + dirarg + 'saved/' + 'mu_size.csv',mu_data,delimiter=',',fmt='%18.16f')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Process data of density-locked mixture simulation')
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
