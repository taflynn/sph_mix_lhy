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

    dens1 = np.empty((Nr+2,frames))
    dens2 = np.empty((Nr+2,frames))
    times = np.empty(frames)
    
    # load in and save relevant data from each time slice 
    for i in range(0,frames):
        psi1 = np.loadtxt('../data/' + dirarg + '/psi1_re_t' + str(i+1) + '.txt',dtype=complex)
        psi2 = np.loadtxt('../data/' + dirarg + '/psi2_re_t' + str(i+1) + '.txt',dtype=complex)
        dens1[:,i] = np.abs(psi1)**2
        dens2[:,i] = np.abs(psi2)**2
        times[i] = dt*i*setup['T_SAVE']
    
    # save density arrays
    with open('../data/' + dirarg + '/' + 'real_spacetime_wav1.npy', 'wb') as f:
        np.save(f,dens1)
    with open('../data/' + dirarg + '/' + 'real_spacetime_wav2.npy', 'wb') as f:
        np.save(f,dens2)

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
