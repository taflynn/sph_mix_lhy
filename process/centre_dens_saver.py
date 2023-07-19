## Import packages
from numpy import empty, absolute, fft, loadtxt, mean, column_stack, arange, savetxt
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

def centre_dens(dirarg, job_no, job_frames):
    if job_no > 0:
        fname = 'config_dens_ulck' + str(job_no) + '.json'
        dirarg = dirarg + str(job_no)
    elif job_no == 0:
        fname = 'config_dens_ulck.json'

    # load in data for simulation
    f = open('../data/' + dirarg + '/' + fname,"r")
    setup = json.loads(f.read())
    f.close()

    dr = setup['dr']
    Nr = setup['Nr']
    dt = setup['DT_COEF']*dr**2

    r = arange(-1/2, (Nr + 3/2), 1)*dr

    if job_frames == 0:
        # if input the number of frames as 0, then use config.json
        frames = setup['RE_T_STEPS']//setup['T_SAVE']
    else:
        # otherwise use the specified number of frames
        frames = job_frames

    n01 = empty(frames)
    n02 = empty(frames)
    del_n = empty(frames)
    times = empty(frames)
    
    # load in and save relevant data from each time slice 
    for i in range(0, frames):
        psi1 = loadtxt('../data/' + dirarg + '/psi1_re_t' + str(i+1) + '.txt', dtype=complex)
        psi2 = loadtxt('../data/' + dirarg + '/psi2_re_t' + str(i+1) + '.txt', dtype=complex)
        n01[i] = absolute(psi1[1])**2
        n02[i] = absolute(psi2[1])**2
        del_n[i] = absolute(psi1[1])**2 - absolute(psi2[1])**2
        times[i] = dt*i*setup['T_SAVE']
    
    # use central density values to extract breathing mode frequencies
    extract = 0 #(np.where(n01>n01[0]))[0][0]
    upper = -1

    # the raw central density data
    n1_data = column_stack((times, n01))
    n2_data = column_stack((times, n02))
    del_n_data = column_stack((times, del_n))    

    # shifted central density data (by the means)
    filter_n1 = n01[extract:upper] - mean(n01[extract:upper])
    filter_n2 = n02[extract:upper] - mean(n02[extract:upper])
    filter_t = times[extract:upper] # - t_real[extract]

    centred_n1_data = column_stack((filter_t, filter_n1))
    centred_n2_data = column_stack((filter_t, filter_n2))

    # saving arrays of central densities
    savetxt('../data/' + dirarg + '/n01_uncentred.csv', n1_data, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + '/n02_uncentred.csv', n2_data, delimiter=',', fmt='%18.16f')

    savetxt('../data/' + dirarg + '/del_n0s.csv', del_n_data, delimiter=',', fmt='%18.16f')

    savetxt('../data/' + dirarg + '/n01_centred.csv', centred_n1_data, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + '/n02_centred.csv', centred_n2_data, delimiter=',', fmt='%18.16f')

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
    parser.add_argument('--frames_no','-fn',
            dest = 'FRAMES',
            type = int,
            required = True,
            nargs = 1)
    args = parser.parse_args()
    centre_dens(args.READ_PATH[0], args.JOB_NO[0], args.FRAMES[0])
