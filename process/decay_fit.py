## Import packages
from numpy import savetxt, empty, exp, sin, load, loadtxt, absolute, column_stack, sqrt, pi, mean
import json
import os
import sys
import shutil
import argparse
from scipy.optimize import curve_fit

sys.path.insert(1, '../')

from main.params_calc import params_dens_ulck
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

    # Define a general damped sine function
    def damp_sin_func(x, a, b, c, d, f):
        return a*exp(-f*x)*sin(b*x + c) + d

    def curve_fitting(t_array, centre_dens, A, B, C, D, F):
        # Extracting the fitted parameter values from the curve_fit function
        popt, popv = curve_fit(damp_sin_func, t_array[0:-1], centre_dens[0:-1],
                               p0=[A, B, C, D, F])
        return popt, popv

    # extract values of density arrays and time arrays after the first max of the oscillating density
    def turning_points(array):
        ''' turning_points(array) -> min_indices, max_indices
        Finds the turning points within an 1D array and returns the indices of the
        minimum and  maximum turning points in two separate lists.
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
    omega01_array = empty((num_sims, 1))
    lowconf_omega01 = empty((num_sims, 1))
    highconf_omega01 = empty((num_sims, 1))
    
    omega02_array = empty((num_sims, 1))
    lowconf_omega02 = empty((num_sims, 1))
    highconf_omega02 = empty((num_sims, 1))
    
    gamma1_array = empty((num_sims, 1))
    lowconf_gamma1 = empty((num_sims, 1))
    highconf_gamma1 = empty((num_sims, 1))
    
    gamma2_array = empty((num_sims, 1))
    lowconf_gamma2 = empty((num_sims, 1))
    highconf_gamma2 = empty((num_sims, 1))
    
    mu1_array = empty((num_sims, 1))
    mu2_array = empty((num_sims, 1))
    
    imb_size_array = empty((num_sims, 1))
    
    # file to save print outputs to
    sys.stdout=open('../data/' + dirarg + 'saved/' +'process.out', "w")
    sys.stdout = Unbuffered(sys.stdout)

    for i in range(0, num_sims):
        # load in simulation parameters
        fname = 'config_dens_ulck' + str(i+1) + '.json'
        print(i)

        # read in data
        f = open('../data/' + dirarg + str(i+1)  + '/' + fname, "r")
        setup = json.loads(f.read())
        f.close()

        # grid spacing
        dr = setup['dr']
        
        # read in theory parameters
        f = open('../data/' + dirarg + str(i+1) + '/theory_params.json', "r")
        theory = json.loads(f.read())
        f.close()
        [alpha, beta, eta, xi, tau, _, _, _, _, N1, N2, dim_pot] = params_dens_ulck(setup['m1'], setup['m2'], 
                                                                                    setup['a11'], setup['a22'], setup['a12'],
                                                                                    setup['N1'], setup['N2'], setup['BALANCE'])

        theory_omg1 = 2*pi*setup['OMEGA1']*sqrt(dim_pot)
        theory_omg2 = 2*pi*setup['OMEGA2']*sqrt(dim_pot)
        
        print('Theoretical omega 1 = ', theory_omg1)
        print('Theoretical omega 2 = ', theory_omg2)
        
        t_cutoff = 0.5 # could add this is an input variable if necessary!
        cut_trap_period = t_cutoff*0.5*(2*pi/theory_omg1)

        # read in data
        t_centre_n1 = loadtxt('../data/' + dirarg + str(i+1) + '/n01_centred.csv', delimiter=",")
        t_centre_n2 = loadtxt('../data/' + dirarg + str(i+1) + '/n02_centred.csv', delimiter=",")
        
        if imbal_size == 'IMBAL':
            imb_size_array[i] = theory['N1']-theory['N2']
        elif imbal_size == 'SIZE':
            N1 = setup['N1']
            N2 = setup['N2']
            perc_diff = np.abs(N1-N2)/N2
            N1 = N1 - (N2*perc_diff)
            N_tot = N1 + N2
            N_lck, xi, tau, n01, n02 = params_dens_lck(setup['m1'], setup['m2'],
                                                       setup['a11'],setup['a22'],
                                                       setup['a12'],N_tot) 
            imb_size_array[i] = (N_lck - 18.65)**0.25

        #extract = (np.where(centre_n1>centre_n1[0]))[0][0]
        t_real = t_centre_n1[:, 0]
        cut_n1 = t_centre_n1[t_real < cut_trap_period, 1]
        cut_n2 = t_centre_n2[t_real < cut_trap_period, 1]
        cut_t = t_centre_n1[t_real < cut_trap_period, 0]

        # fit central density oscillations to damped sine curve
        [fitted_params1, cov1] = curve_fitting(cut_t, cut_n1,
                                               absolute(cut_n1.max() - cut_n1.min()),
                                               0.3, pi, mean(cut_n1), 0.03)
        [fitted_params2, cov2] = curve_fitting(cut_t, cut_n2,
                                               absolute(cut_n2.max() - cut_n2.min()),
                                               0.3, pi, mean(cut_n1), 0.03)
        
        # extract breathing mode frequency and decay rate of oscillation for component 1
        omega01_array[i] = absolute(fitted_params1[1])
        lowconf_omega01[i] = absolute(fitted_params1[1]) - 2*sqrt(cov1[1, 1])
        highconf_omega01[i] = absolute(fitted_params1[1]) + 2*sqrt(cov1[1, 1])
        
        gamma1_array[i] = fitted_params1[-1]
        lowconf_gamma1[i] = fitted_params1[-1] - 2*sqrt(cov1[-1, -1])
        highconf_gamma1[i] = fitted_params1[-1] + 2*sqrt(cov1[-1, -1])
        
        # extract breathing mode frequency and decay rate of oscillation for component 2
        omega02_array[i] = absolute(fitted_params2[1])
        lowconf_omega02[i] = absolute(fitted_params2[1]) - 2*sqrt(cov2[1, 1])
        highconf_omega02[i] = absolute(fitted_params2[1]) + 2*sqrt(cov2[1, 1])
        
        gamma2_array[i] = fitted_params2[-1]
        lowconf_gamma2[i] = fitted_params2[-1] - 2*sqrt(cov2[-1, -1])
        highconf_gamma2[i] = fitted_params2[-1] + 2*sqrt(cov2[-1, -1])
    
        mu1_array[i] = -1*theory['mu1']
        mu2_array[i] = -1*theory['mu2']
    
    # concatenating data into 2D array of omega and gamma for component 1 
    omega01_data = column_stack((imb_size_array, omega01_array))
    lowconf_omega01_data = column_stack((imb_size_array, lowconf_omega01))
    highconf_omega01_data = column_stack((imb_size_array, highconf_omega01))

    gamma1_data = column_stack((imb_size_array, gamma1_array))
    lowconf_gamma1_data = column_stack((imb_size_array, lowconf_gamma1))
    highconf_gamma1_data = column_stack((imb_size_array, highconf_gamma1))
    
    # concatenating data into 2D array of omega and gamma for component 2
    omega02_data = column_stack((imb_size_array, omega02_array))
    lowconf_omega02_data = column_stack((imb_size_array, lowconf_omega02))
    highconf_omega02_data = column_stack((imb_size_array, highconf_omega02))
    
    gamma2_data = column_stack((imb_size_array, gamma2_array))
    lowconf_gamma2_data = column_stack((imb_size_array, lowconf_gamma2))
    highconf_gamma2_data = column_stack((imb_size_array, highconf_gamma2))
    
    mu1_data = column_stack((imb_size_array, mu1_array))
    mu2_data = column_stack((imb_size_array, mu2_array))

    # saving arrays of omega's and confidence intervals
    ## omega 1's
    savetxt('../data/' + dirarg + 'saved/' + 'omega01.csv',
            omega01_data, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved/' + 'omega01_low.csv',
            lowconf_omega01_data, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved/' + 'omega01_high.csv',
            highconf_omega01_data, delimiter=',', fmt='%18.16f')
    
    ## omega 2's
    savetxt('../data/' + dirarg + 'saved/' + 'omega02.csv',
            omega02_data,delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved/' + 'omega02_low.csv',
            lowconf_omega02_data, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved/' + 'omega02_high.csv',
            highconf_omega02_data, delimiter=',', fmt='%18.16f')
    
    ## gamma 1's
    savetxt('../data/' + dirarg + 'saved/' + 'gamma1.csv',
            gamma1_data, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved/' + 'gamma1_low.csv',
            lowconf_gamma1_data, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved/' + 'gamma1_high.csv',
            highconf_gamma1_data, delimiter=',', fmt='%18.16f')
    
    ## gamma 2's
    savetxt('../data/' + dirarg + 'saved/' + 'gamma2.csv',
            gamma2_data, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved/' + 'gamma2_low.csv',
            lowconf_gamma2_data, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved/' + 'gamma2_high.csv',
            highconf_gamma2_data, delimiter=',', fmt='%18.16f')

    ## mu's
    savetxt('../data/' + dirarg + 'saved/' + 'mu1.csv', 
            mu1_data, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved/' + 'mu2.csv', 
            mu2_data, delimiter=',', fmt='%18.16f')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Fitting decay rate of breathing modes')
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
    run_ulck_process(args.READ_PATH[0], args.NUM_SIMS[0], args.IMBAL_SIZE[0])
