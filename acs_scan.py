import json

import numpy as np

import argparse

from numpy import pi

import shutil

import os

def acs_scan(acs_start,acs_step,acs_end):
   
    acs_array = np.arange(acs_start,acs_end+1,acs_step)

    # open and load the experimental parameters file
    f = open('./config_dens_ulck.json','r')
    exp = json.loads(f.read())
    f.close()
    
    # simulation counter for specifying job array
    dir_counter = 1

    for i in acs_array:
        # overwrite a11 with local a11
        exp['a11'] = i

        # save this overwritten experimental file to a temporary file
        with open('./config_dens_ulck'+str(dir_counter)+'.json','w') as f:
            json.dump(exp,f,indent=4)

        # iterate the job for the job array
        dir_counter+=1

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'experimenal -> theory parameters')
    parser.add_argument('--acs_start','-start',
            dest = 'START',
            type = float,
            required = True,
            nargs = 1)
    parser.add_argument('--acs_step','-step',
            dest = 'STEP',
            type = float,
            required = True,
            nargs = 1)
    parser.add_argument('--acs_end','-end',
            dest = 'END',
            type = float,
            required = True,
            nargs = 1)
    args = parser.parse_args()
    acs_scan(args.START[0],args.STEP[0],args.END[0])
