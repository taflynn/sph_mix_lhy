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

def run_ulck_process(dirarg,num_sims,imbal_size):
    if imbal_size == 'IMBAL': 
        sim_style = 'imbal'
    elif imbal_size == 'SIZE':
        sim_style = 'size'
    else:
        print('Error: Need to specify style of simulation, i.e., imbalance or droplet size')
        sim_style = 'N/A'

    # Define a general damped sine function
    def damp_sin_func(x,a,b,c,d,f):
        return a*np.exp(-f*x)*np.sin(b*x + c) + d
    # function to fit the 
    def curve_fitting(t_array,centre_dens):
        # Initial guess parameters for the curve_fit function
        A = np.max(centre_dens) # initial guess of amplitude
        B = 0.5 # very rough guess for frequency
        C = np.pi/2 # default input value for phase shift
        D = (np.max(centre_dens)-np.min(centre_dens))/2.0 # initial guess for the y shift of the sine curve
        F = 0.1
    
        # Extracting the fitted parameter values from the curve_fit function
        popt, cov = curve_fit(damp_sin_func,t_array[0:-1],centre_dens[0:-1],p0=[A,B,C,D,F])
        return popt,cov
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
    omega01_array = np.empty((num_sims,1))
    lowconf_omega01 = np.empty((num_sims,1))
    highconf_omega01 = np.empty((num_sims,1))
    
    omega02_array = np.empty((num_sims,1))
    lowconf_omega02 = np.empty((num_sims,1))
    highconf_omega02 = np.empty((num_sims,1))
    
    gamma1_array = np.empty((num_sims,1))
    lowconf_gamma1 = np.empty((num_sims,1))
    highconf_gamma1 = np.empty((num_sims,1))
    
    gamma2_array = np.empty((num_sims,1))
    lowconf_gamma2 = np.empty((num_sims,1))
    highconf_gamma2 = np.empty((num_sims,1))
    
    imb_size_array = np.empty((num_sims,1))
    
    # file to save print outputs to
    sys.stdout=open('../data/' + dirarg + 'saved/' +'process.out',"w")
    sys.stdout = Unbuffered(sys.stdout)

    for i in range(1,num_sims):
        # load in simulation parameters
        fname = 'config_dens_ulck' + str(i + 1) + '.json'
        print(i)

        # read in data
        f = open('../data/' + dirarg + str(i + 1)  + '/' + fname,"r")
        setup = json.loads(f.read())
        f.close()

        # grid spacing
        dr = setup['dr']

        if imbal_size == 'IMBAL':
            imb_size_array[i] = ((setup['N1']-setup['N2'])/setup['N2'])*100.0
        elif imbal_size == 'SIZE':
            N1 = setup['N1']
            N2 = setup['N2']
            perc_diff = np.abs(N1-N2)/N2
            N1 = N1 - (N2*perc_diff)
            N_tot = N1 + N2
            N_lck,xi,tau,n01,n02 = params_dens_lck(setup['m1'],setup['m2'],setup['a11'],setup['a22'],setup['a12'],N_tot) 
            imb_size_array[i] = (N_lck - 18.65)**0.25

        # main data read in from NumPy files
        r = np.load('../data/' + dirarg + str(i+1) + '/r_array.npy')
        t_real = np.load('../data/' + dirarg + str(i+1) + '/real_t_array.npy')
        psi1 = np.load('../data/' + dirarg + str(i+1) + '/real_spacetime_wav1.npy')
        psi2 = np.load('../data/' + dirarg + str(i+1) + '/real_spacetime_wav2.npy')
        
        # extract central density values
        centre_n1 = (np.abs(psi1[1,:])**2)
        centre_n2 = (np.abs(psi2[1,:])**2)
        
        # cutoff the first transient of the oscillations 
        [idx_min,idx_max] = turning_points(centre_n1) 
        if centre_n1[idx_max[0]] < centre_n1[0]:
            extract = idx_max[1]
        else:
            extract = idx_max[0]
        cut_n1 = centre_n1[extract:]
        cut_n2 = centre_n2[extract:]
        cut_t = t_real[extract:]

        # fit central density oscillations to damped sine curve
        [fitted_params1,cov1] = curve_fitting(cut_t,cut_n1)
        [fitted_params2,cov2] = curve_fitting(cut_t,cut_n2)
        
        # extract breathing mode frequency and decay rate of oscillation for component 1
        omega01_array[i] = fitted_params1[1]
        lowconf_omega01[i] = fitted_params1[1] - 2*np.sqrt(cov1[1,1])
        highconf_omega01[i] = fitted_params1[1] + 2*np.sqrt(cov1[1,1])
        
        gamma1_array[i] = fitted_params1[-1]
        lowconf_gamma1[i] = fitted_params1[-1] - 2*np.sqrt(cov1[-1,-1])
        highconf_gamma1[i] = fitted_params1[-1] + 2*np.sqrt(cov1[-1,-1])
        
        # extract breathing mode frequency and decay rate of oscillation for component 2
        omega02_array[i] = fitted_params2[1]
        lowconf_omega02[i] = fitted_params2[1] - 2*np.sqrt(cov2[1,1])
        highconf_omega02[i] = fitted_params2[1] + 2*np.sqrt(cov2[1,1])
        
        gamma2_array[i] = fitted_params2[-1]
        lowconf_gamma2[i] = fitted_params2[-1] - 2*np.sqrt(cov2[-1,-1])
        highconf_gamma2[i] = fitted_params2[-1] + 2*np.sqrt(cov2[-1,-1])
    
    # concatenating data into 2D array of omega and gamma for component 1 
    omega01_data = np.column_stack((imb_size_array,omega01_array))
    lowconf_omega01_data = np.column_stack((imb_size_array,lowconf_omega01))
    highconf_omega01_data = np.column_stack((imb_size_array,highconf_omega01))

    gamma1_data = np.column_stack((imb_size_array,gamma1_array))
    lowconf_gamma1_data = np.column_stack((imb_size_array,lowconf_gamma1))
    highconf_gamma1_data = np.column_stack((imb_size_array,highconf_gamma1))
    
    # concatenating data into 2D array of omega and gamma for component 2
    omega02_data = np.column_stack((imb_size_array,omega02_array))
    lowconf_omega02_data = np.column_stack((imb_size_array,lowconf_omega02))
    highconf_omega02_data = np.column_stack((imb_size_array,highconf_omega02))
    
    gamma2_data = np.column_stack((imb_size_array,gamma2_array))
    lowconf_gamma2_data = np.column_stack((imb_size_array,lowconf_gamma2))
    highconf_gamma2_data = np.column_stack((imb_size_array,highconf_gamma2))
    
    # saving arrays of omega's and confidence intervals
    np.savetxt('../data/' + dirarg + 'saved/' + 'omega01_' + sim_style + '.csv',omega01_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'omega01_low_' + sim_style + '.csv',lowconf_omega01_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'omega01_high_' + sim_style + '.csv',highconf_omega01_data,delimiter=',',fmt='%18.16f')

    
    np.savetxt('../data/' + dirarg + 'saved/' + 'omega02_' + sim_style + '.csv',omega02_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'omega02_low_' + sim_style + '.csv',lowconf_omega02_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'omega02_high_' + sim_style + '.csv',highconf_omega02_data,delimiter=',',fmt='%18.16f')

    
    np.savetxt('../data/' + dirarg + 'saved/' + 'gamma1_' + sim_style + '.csv',gamma1_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'gamma1_low_' + sim_style + '.csv',lowconf_gamma1_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'gamma1_high_' + sim_style + '.csv',highconf_gamma1_data,delimiter=',',fmt='%18.16f')

    
    np.savetxt('../data/' + dirarg + 'saved/' + 'gamma2_' + sim_style + '.csv',gamma2_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'gamma2_low_' + sim_style + '.csv',lowconf_gamma2_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'gamma2_high_' + sim_style + '.csv',highconf_gamma2_data,delimiter=',',fmt='%18.16f')

    
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
