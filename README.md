# SPHERICAL BOSE-BOSE MIXTURE W/ LHY #

##Intro:
* Overview: These codes are to solve the Bose-Bose mixture, with the Lee-Huang-Yang (LHY) correction, in a spherically symmetric geometry (i.e. the density of the mixture is a function of radius only).

* Time-stepping scheme: The code uses a 4th-order Runge-Kutta Scheme.

* Boundary conditions: The code has options for boundary conditions in the 'boundary.py' file though the only current options are: 1. Neumann boundary conditions at both ends; 2. Neumann boundary condition at centre of the box w/ density set to zero at the far end of the box.

## Code Structure:
### Main files and directories
* The main body of the code is within the 'main' directory.
* json_init is a directory which contains two files: default.config_dens_lck.json and default.config_dens_ulck.json. It is best to copy these into the sph_mix_lhy directory and rename them config_dens_lck.json and config_dens_ulck.json. These are the config files for each simulation and will be described below
* mix_sph_lhy.py is the file which reads in the simulation data and runs the codes in main
* running.py is the file which runs mix_sph_lhy.py and saves most of the data and files that each simulation produces
### HPC workflow files
These are files that are designed to be ran on HPC clusters (using the SLURM queueing system):
* drop_run.sh is the shell scipt to run a single simulation via SLURM
* drop_array.sh is the shell script to run a job array of simulations via SLURM
* generator_imbal.sh is the shell script to generate config files for a job array of simulations for various atom number imbalances
* generator_size.sh is the shell script to generate config files for a job array of simulations for various droplet sizes
### Processing files
The directory './process' contains files for processing output of simulations, particularly for real time breathing mode simulations:
* imbal_breath_fitting.py is a file which reads in multiple simulations (from a job array) and fits the oscillation of the central density of a droplet and outputs values of the breathing mode frequency and breathing mode decay rate
* breath_fit.sh is the shell script to run the imbal_breath_fitting.py file on the SLURM queueing system
* process_ulck.py is a file which generates movies and other output for an individual simulation

## Describing a mixture:
A Bose-Bose mixture is here defined by:
* a11 - intraspecies scattering length for component 1
* a12 - interspecies scattering length
* a22 - intraspecies scattering length for component 2
* m1 - atomic mass of component 1
* m2 - atomic mass of component 2
We also require an atom number for the system and this is where the first separation occurs. 

As mention the mixture can be density-locked (DENS_LCK = TRUE) or density-unlocked (DENS_LCK = FALSE). When in the density-locked regime, you need only specify the total number of atoms, N. However, when in the density-unlocked regime, you must specify the atom number for both components (i.e. N1 and N2).

One other important way to specify your system is the trapping potential. In the spherically symmetric case we are somewhat limited so the main choices are: 1. No potential; 2. Harmonic potential. The only subtlty is that in the density-locked case you can only specify harmonic trap frequency (OMEGA) but in the density-unlocked case you can impose different trapping potentials on each component (via OMEGA1 and OMEGA2).

## Configuration for a simulation
