## Import packages
from numpy import load, loadtxt, savetxt, column_stack, absolute, zeros, trapz, pi, sqrt
import json
import os
import sys
import shutil
import argparse

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

def load_in_cal_norms(dirarg, sim_num, rad_cutoff):
    # load in config file for simulation
    if sim_num == 0:
        fname = 'config_dens_ulck.json'
    else:
        fname = 'config_dens_ulck' + str(sim_num) + '.json'    
    f = open('../data/' + dirarg + '/' + fname,"r")
    setup = json.loads(f.read())
    f.close()

    no_of_frames = setup['RE_T_STEPS']//setup['T_SAVE']
    dr = setup['dr']
    r = load('../data/' + dirarg + '/r_array.npy')

    dt = setup['DT_COEF']*dr**2.0

    [alpha, beta, eta, xi, tau, _, _, _, _, N1, N2, dim_pot] = params_dens_ulck(setup['m1'], setup['m2'], 
                                                                                setup['a11'], setup['a22'], setup['a12'],
                                                                                setup['N1'], setup['N2'], setup['BALANCE'])

    theory_omg1 = 2*pi*setup['OMEGA1']*sqrt(dim_pot)
    theory_omg2 = 2*pi*setup['OMEGA2']*sqrt(dim_pot)
    
    print('Theoretical omega 1 = ', theory_omg1)
    print('Theoretical omega 2 = ', theory_omg2)

    psi1_saved = zeros((setup['Nr']+2, no_of_frames), dtype=complex)
    psi2_saved = zeros((setup['Nr']+2, no_of_frames), dtype=complex)

    norms_1 = zeros(no_of_frames, dtype=float)
    norms_2 = zeros(no_of_frames, dtype=float)
    deln_r0 = zeros(no_of_frames, dtype=float)
    times = zeros(no_of_frames, dtype=float)

    for i in range(0, no_of_frames):
        psi1_saved[:, i] = loadtxt('../data/' + dirarg + '/psi1_re_t' + str(i+1) + '.txt',
                                  dtype=complex)
        psi2_saved[:, i] = loadtxt('../data/' + dirarg + '/psi2_re_t' + str(i+1) + '.txt',
                                  dtype=complex)
        R_1 = r[(absolute(psi1_saved[:, i])**2.0 
                 > rad_cutoff*(absolute(psi1_saved[:, i])**2.0).max())][-1]
        R_2 = r[(absolute(psi2_saved[:, i])**2.0 
                 > rad_cutoff*(absolute(psi2_saved[:, i])**2.0).max())][-1]
        norms_1[i] = 4*pi*trapz(r[r<R_1]**2.0*absolute(psi1_saved[r<R_1, i])**2.0)*dr
        norms_2[i] = 4*pi*trapz(r[r<R_2]**2.0*absolute(psi2_saved[r<R_2, i])**2.0)*dr    
        deln_r0[i] = absolute(psi1_saved[1, i])**2.0 - absolute(psi2_saved[1, i])**2.0
        times[i] = dt*i*setup['T_SAVE']

    return psi1_saved, psi2_saved, norms_1, norms_2, deln_r0, times, theory_omg1, theory_omg2

def run_norms_process(dirarg, sim_num, rad_cutoff):

    [psi1_cat, psi2_cat, norms_1, norms_2, deln_r0, times, theory_omg1, theory_omg2] = load_in_cal_norms(dirarg, sim_num, rad_cutoff)
    
    # new outputs, as function of N1
    norms1_vs_t = column_stack((times, norms_1))
    norms2_vs_t = column_stack((times, norms_2))
    deln_r0_vs_t = column_stack((times, deln_r0))

    # saving arrays of omega's and confidence intervals
    savetxt('../data/' + dirarg + '/norms1_vs_t.csv',
            norms1_vs_t, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + '/norms2_vs_t.csv',
            norms2_vs_t, delimiter=',', fmt='%18.16f')
    savetxt('../data/' + dirarg + '/deln_r0_vs_t.csv',
            deln_r0_vs_t, delimiter=',', fmt='%18.16f')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Extract norms and central density difference in real time')
    parser.add_argument('--read_path','-rp',
            dest = 'READ_PATH',
            type = str,
            required = True,
            nargs = 1)
    parser.add_argument('--sim_num','-sn',
            dest = 'SIM_NUM',
            type = int,
            required = True,
            nargs = 1)
    parser.add_argument('--rad_cut','-rc',
            dest = 'RAD_CUTOFF',
            type = float,
            required = True,
            nargs = 1)
    args = parser.parse_args()
    run_norms_process(args.READ_PATH[0], args.SIM_NUM[0], args.RAD_CUTOFF[0])
