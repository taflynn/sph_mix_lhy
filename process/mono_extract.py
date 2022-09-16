## Import packages
import numpy as np
from numpy import fft
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

def mu_omg_jobs(dirarg,job_no):
    dirarg = dirarg + str(job_no)

    fname = 'config_dens_ulck' + str(job_no) + '.json'
    
    # load in data for simulation
    f = open('../data/' + dirarg + '/' + fname,"r")
    setup = json.loads(f.read())
    f.close()

    dr = setup['dr']
    Nr = setup['Nr']
    dt = setup['DT_COEF']*dr**2

    r = np.arange(-1/2,(Nr + 3/2),1)*dr

    frames = setup['RE_T_STEPS']//setup['T_SAVE']

    n01 = np.empty(frames)
    n02 = np.empty(frames)
    times = np.empty(frames)
    
    # load in and save relevant data from each time slice 
    for i in range(0,frames):
        psi1 = np.loadtxt('../data/' + dirarg + '/psi1_re_t' + str(i+1) + '.txt',dtype=complex)
        psi2 = np.loadtxt('../data/' + dirarg + '/psi2_re_t' + str(i+1) + '.txt',dtype=complex)
        n01[i] = np.abs(psi1[1])**2
        n02[i] = np.abs(psi2[1])**2
        times[i] = dt*i*setup['T_SAVE']
    
    # use central density values to extract breathing mode frequencies
    extract = (np.where(n01>n01[0]))[0][0]
    upper = -1

    filter_n1 = n01[extract:upper]
    filter_n2 = n02[extract:upper]
    filter_t = times[extract:upper]# - t_real[extract]

    signal1 = fft.fft(filter_n1 - np.mean(filter_n1))
    signal2 = fft.fft(filter_n2 - np.mean(filter_n2))
    
    freqs = fft.fftshift(fft.fftfreq(len(signal1),d=1/filter_t.size)*(2*np.pi/(filter_t[-1])))

    power1 = np.abs(fft.fftshift(signal1))**2
    power2 = np.abs(fft.fftshift(signal2))**2 
    omega1 = np.abs(np.max(freqs[power1 == np.max(power1)]))
    omega2 = np.abs(np.max(freqs[power2 == np.max(power2)]))
    if omega2 < 0.15:
        omega2 = np.abs(np.max(freqs[power2 == np.max(power2)]))
        freqs_filt = freqs[freqs>0.1]
        power2_filt = power2[freqs>0.1]
        omega2 = np.abs(np.max(freqs_filt[power2_filt == np.max(power2_filt)]))

    perc_imb = ((setup['N1']-setup['N2'])/setup['N2'])*100.0

    theory_fname = 'theory_params.json'
    
    # load in data for simulation
    theory_vals = open('../data/' + dirarg + '/' + theory_fname,"r")
    theories = json.loads(theory_vals.read())
    theory_vals.close()
    mu1 = theories['mu1']
    mu2 = theories['mu2']

    return perc_imb,mu1,mu2,omega1,omega2

def run_ulck_process(dirarg,num_sims):
    # generating empty arrays for saving fitted parameters and confidence intervals
    mu1_array = np.empty((num_sims,1))
    mu2_array = np.empty((num_sims,1))
    omg1_array = np.empty((num_sims,1))
    omg2_array = np.empty((num_sims,1))
    perc_array = np.empty((num_sims,1))
    
    # file to save print outputs to
    sys.stdout=open('../data/' + dirarg + 'saved/' +'process.out',"w")
    sys.stdout = Unbuffered(sys.stdout)

    for i in range(1,num_sims+1):
        
        [perc_imb,mu1,mu2,omg1,omg2] = mu_omg_jobs(dirarg,i)
        perc_array[i-1] = perc_imb
        mu1_array[i-1] = -mu1
        mu2_array[i-1] = -mu2
        omg1_array[i-1] = omg1
        omg2_array[i-1] = omg2
        print('Completed job ' + str(i) + 'of' + str(num_sims))
 
    # concatenating data into 2D array of omega and gamma for component 1 
    mu1_data = np.column_stack((perc_array,mu1_array))
    mu2_data = np.column_stack((perc_array,mu2_array))
    omg1_data = np.column_stack((perc_array,omg1_array))
    omg2_data = np.column_stack((perc_array,omg2_array))
    
    # saving arrays of omega's and confidence intervals
    np.savetxt('../data/' + dirarg + 'saved/' + 'chempot1.csv',mu1_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'chempot2.csv',mu2_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'omega1.csv',omg1_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + 'saved/' + 'omega2.csv',omg2_data,delimiter=',',fmt='%18.16f')

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
