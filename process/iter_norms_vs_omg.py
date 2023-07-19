## Import packages
from numpy import load, loadtxt, savetxt, column_stack, absolute, zeros, empty, trapz, pi, sqrt
import json
import os
import sys
import shutil
import argparse
from norms_vs_t import load_in_cal_norms

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

def omg_norms_iter(dirarg, no_of_sims, rad_cutoff):
    
    omega_1_dat = empty(no_of_sims, dtype=float)
    omega_2_dat = empty(no_of_sims, dtype=float)
    harm_len_1_dat = empty(no_of_sims, dtype=float)
    harm_len_2_dat = empty(no_of_sims, dtype=float)
    norms_1_lim = empty(no_of_sims, dtype=float)
    norms_2_lim = empty(no_of_sims, dtype=float)
    deln_r0_lim = empty(no_of_sims, dtype=float)
    
    # file to save print outputs to
    sys.stdout=open('../data/' + dirarg + 'saved_radcut' + str(rad_cutoff)  + '/processing.out', "w")
    sys.stdout = Unbuffered(sys.stdout) 
   
    for i in range(0, no_of_sims):
        dirarg_current = dirarg + str(i+1)
        print('current simulation: ', dirarg_current)
        [psi1_cat, psi2_cat, norms_1, norms_2, deln_r0, times, theory_omg1, theory_omg2] = load_in_cal_norms(dirarg_current, i+1, rad_cutoff)
        norms_1_lim[i] = norms_1[-1]
        norms_2_lim[i] = norms_2[-1]
        deln_r0_lim[i] = deln_r0[-1]
        omega_1_dat[i] = theory_omg1
        omega_2_dat[i] = theory_omg2
        harm_len_1_dat[i] = 1.0/sqrt(theory_omg1)
        harm_len_2_dat[i] = 1.0/sqrt(theory_omg2)
        print('finished simulation: ', dirarg_current)
        print((i+1)/no_of_sims*100, '% completed')
    
    print('all omega 1 = ', omega_1_dat)
    print('all omega 2 = ', omega_2_dat)

    # new outputs, as function of N1
    norms1_vs_omg = column_stack((omega_1_dat, norms_1_lim))
    norms2_vs_omg = column_stack((omega_1_dat, norms_2_lim))
    deln_r0_vs_omg = column_stack((omega_1_dat, deln_r0_lim))

    norms1_vs_len = column_stack((harm_len_1_dat, norms_1_lim))
    norms2_vs_len = column_stack((harm_len_1_dat, norms_2_lim))
    deln_r0_vs_len = column_stack((harm_len_1_dat, deln_r0_lim))

    print('writing results')
    # saving arrays of omega's and confidence intervals
    savetxt('../data/' + dirarg + 'saved_radcut' + str(rad_cutoff) + '/norms1_vs_omg.csv', 
            norms1_vs_omg, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved_radcut' + str(rad_cutoff) + '/norms2_vs_omg.csv',
            norms2_vs_omg, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved_radcut' + str(rad_cutoff) + '/deln_r0_vs_omg.csv',
            deln_r0_vs_omg, delimiter=',', fmt='%18.16f')

    savetxt('../data/' + dirarg + 'saved_radcut' + str(rad_cutoff) + '/norms1_vs_len.csv',
            norms1_vs_len, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved_radcut' + str(rad_cutoff) + '/norms2_vs_len.csv',
            norms2_vs_len, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + 'saved_radcut' + str(rad_cutoff) + '/deln_r0_vs_len.csv',
            deln_r0_vs_len, delimiter=',', fmt='%18.16f')

    print('complete!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Extract norms and central density difference in real time')
    parser.add_argument('--read_path','-rp',
            dest = 'READ_PATH',
            type = str,
            required = True,
            nargs = 1)
    parser.add_argument('--num_of_sims','-ns',
            dest = 'SIM_NUM',
            type = int,
            required = True,
            nargs = 1)
    parser.add_argument('--rad_cut','-rc',
            dest = 'RAD_CUT',
            type = float,
            required = True,
            nargs = 1)
    args = parser.parse_args()
    omg_norms_iter(args.READ_PATH[0], args.SIM_NUM[0], args.RAD_CUT[0])
