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
        fname = 'config_dens_lck' + str(job_no) + '.json'
    if job_no == 0:
        fname = 'config_dens_lck.json'
    # load in data for simulation
    f = open('../data/' + dirarg + '/' + fname,"r")
    setup = json.loads(f.read())
    f.close()

    dr = setup['dr']
    Nr = setup['Nr']
    dt = setup['DT_COEF']*dr**2

    r = np.arange(-1/2,(Nr + 3/2),1)*dr

    frames = setup['RE_T_STEPS']//setup['T_SAVE']

    n0 = np.empty(frames)
    times = np.empty(frames)
    
    # load in and save relevant data from each time slice 
    for i in range(0,frames):
        phi = np.loadtxt('../data/' + dirarg + '/phi_re_t' + str(i+1) + '.txt',dtype=complex)
        n0[i] = np.abs(phi[1])**2
        times[i] = dt*i*setup['T_SAVE']
    
    # use central density values to extract breathing mode frequencies
    extract = (np.where(n0>n0[0]))[0][0]
    upper = -1

    filter_n = n0[extract:upper] - np.mean(n0[extract:upper])
    filter_t = times[extract:upper]# - t_real[extract]
    
    phi_data = np.column_stack((filter_t,filter_n))
    # saving arrays of omega's and confidence intervals
    np.savetxt('../data/' + dirarg + '/n0.csv',phi_data,delimiter=',',fmt='%18.16f')

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
