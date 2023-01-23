## Import packages
import numpy as np
import os
import sys
import shutil
import argparse

from centre_dens_saver import centre_dens

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

    for i in range(1, no_jobs + 1):
        dirarg_job = dirarg + str(i)
        print("job directory is: ", dirarg_job)
        centre_dens(dirarg_job, i)
        print("Finished job", i)
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Process data of density-unlocked mixture simulation')
    parser.add_argument('--read_path','-rp',
            dest = 'READ_PATH',
            type = str,
            required = True,
            nargs = 1)
    parser.add_argument('--no_of_jobs','-nj',
            dest = 'JOB_NO',
            type = int,
            required = True,
            nargs = 1)
    args = parser.parse_args()
    main(args.READ_PATH[0], args.JOB_NO[0])
