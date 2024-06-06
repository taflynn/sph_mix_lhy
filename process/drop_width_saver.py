## Import packages
from numpy import empty, absolute, fft, loadtxt, mean, column_stack, arange, savetxt, trapz, pi, zeros
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

def main(dirarg, no_jobs):

    omg1s = zeros(no_jobs)

    n1_widths = zeros(no_jobs)
    n2_widths = zeros(no_jobs)

    for i in range(1, no_jobs + 1):
        dirarg_job = dirarg + str(i)
        print("job directory is: ", dirarg_job)
        n1_widths[i-1], n2_widths[i-1], omg1s[i-1] = drop_width(dirarg_job, i)
        print("Finished job", i)
    
    n1_data = column_stack((omg1s, n1_widths))
    n2_data = column_stack((omg1s, n2_widths))

    # saving arrays of central densities
    savetxt('../data/' + dirarg + 'saved/n1_widths_vs_omg1s.csv', n1_data, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved/n2_widths_vs_omg1s.csv', n2_data, delimiter=',', fmt='%18.16f')

    return 0

def drop_width(dirarg, job_no):
    if job_no > 0:
        fname = 'config_dens_ulck' + str(job_no) + '.json'
    elif job_no == 0:
        fname = 'config_dens_ulck.json'

    # load in data for simulation
    f = open('../data/' + dirarg + '/' + fname,"r")
    setup = json.loads(f.read())
    f.close()

    dr = setup['dr']
    Nr = setup['Nr']
    dt = setup['DT_COEF']*dr**2
    omgr1 = setup['OMEGA1']

    r = arange(-1/2, (Nr + 3/2), 1)*dr
        
    dens1 = loadtxt('../data/' + dirarg + '/imag_fin_dens1.csv', delimiter=',', dtype=float)
    dens2 = loadtxt('../data/' + dirarg + '/imag_fin_dens2.csv', delimiter=',', dtype=float)
    n1_width = 4*pi*trapz(r**4*dens1[:, 1])*dr
    n2_width = 4*pi*trapz(r**4*dens2[:, 1])*dr
    
    return n1_width, n2_width, omgr1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Compute droplet widths across trapping potentials')
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
    main(args.READ_PATH[0], args.JOB_NO[0])
