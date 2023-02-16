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

def centre_dens(dirarg,job_no):
    if job_no > 0:
        fname = 'config_dens_ulck' + str(job_no) + '.json'
    if job_no == 0:
        fname = 'config_dens_ulck.json'
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
        print(i)
        times[i] = dt*i*setup['T_SAVE']
    
    # use central density values to extract breathing mode frequencies
    extract = 0#(np.where(n01>n01[0]))[0][0]
    upper = -1

    filter_n1 = n01[extract:upper] - np.mean(n01[extract:upper])
    filter_n2 = n02[extract:upper] - np.mean(n02[extract:upper])
    filter_t = times[extract:upper]# - t_real[extract]
    
    psi1_data = np.column_stack((filter_t,filter_n1))
    psi2_data = np.column_stack((filter_t,filter_n2))
    # saving arrays of omega's and confidence intervals
    np.savetxt('../data/' + dirarg + '/n01.csv',psi1_data,delimiter=',',fmt='%18.16f')
    np.savetxt('../data/' + dirarg + '/n02.csv',psi2_data,delimiter=',',fmt='%18.16f')

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Process data of density-unlocked mixture simulation')
    parser.add_argument('--read_path','-rp',
            dest = 'READ_PATH',
            type = str,
            required = True,
            nargs = 1)
    parser.add_argument('--job_no','-jn',
            dest = 'JOB_NO',
            type = int,
            required = True,
            nargs = 1)
    args = parser.parse_args()
    centre_dens(args.READ_PATH[0],args.JOB_NO[0])
