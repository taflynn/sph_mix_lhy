# SPHERICAL BOSE-BOSE MIXTURE W/ LHY #

## Intro:
* Overview: These codes are to solve the Bose-Bose mixture, with the Lee-Huang-Yang (LHY) correction, in a spherically symmetric geometry (i.e. the density of the mixture is a function of radius only).

* Time-stepping scheme: The code uses a 4th-order Runge-Kutta Scheme.

* Boundary conditions: The code has options for boundary conditions in the `boundary.py` file though the only current options are: 1. Neumann boundary conditions at both ends; 2. Neumann boundary condition at centre of the box with density set to zero at the far end of the box.

## Code Structure:
### Main files and directories
* The main body of the code is within the 'main' directory.
* `mix_sph_lhy/json_init` is a directory which contains two files: `default.config_dens_lck.json` and `default.config_dens_ulck.json`. It is best to copy these into the `sph_mix_lhy` directory and rename them `config_dens_lck.json` and `config_dens_ulck.json`. These are the configuration files for each simulation and will be described below
* `mix_sph_lhy.py` is the file which reads in the simulation data and runs the codes in main
* `running.py` is the file which runs `mix_sph_lhy.py` and saves most of the data and files that each simulation produces
### HPC workflow files
These are files that are designed to be ran on HPC clusters (using the SLURM queueing system):
* `drop_run.sh` is the shell scipt to run a single simulation via SLURM
* `drop_array.sh` is the shell script to run a job array of simulations via SLURM
* `generator_imbal.sh` is the shell script to generate config files for a job array of simulations for various atom number imbalances
* `generator_size.sh` is the shell script to generate config files for a job array of simulations for various droplet sizes
### Processing files
The directory `mix_sph_lhy/process` contains files for processing output of simulations, particularly for real time breathing mode simulations:
* `imbal_breath_fitting.py` is a file which reads in multiple simulations (from a job array) and fits the oscillation of the central density of a droplet and outputs values of the breathing mode frequency and breathing mode decay rate
* `breath_fit.sh` is the shell script to run the `imbal_breath_fitting.py` file on the SLURM queueing system
* `process_ulck.py` is a file which generates movies and other output for an individual simulation

## Describing a mixture:
A Bose-Bose mixture is here defined by:
* `a11` - intraspecies scattering length for component 1
* `a12` - interspecies scattering length
* `a22` - intraspecies scattering length for component 2
* `m1` - atomic mass of component 1
* `m2` - atomic mass of component 2
We also require an atom number for the system and this is where the first separation occurs. 

As mention the mixture can be density-locked (`DENS_LCK = 1`) or density-unlocked (`DENS_LCK = 0`). When in the density-locked regime, you need only specify the total number of atoms, `N`. However, when in the density-unlocked regime, you must specify the atom number for both components (i.e. `N1` and `N2`).

One other important way to specify your system is the trapping potential. In the spherically symmetric case we are somewhat limited so the main choices are: 1. No potential; 2. Harmonic potential. The only subtlty is that in the density-locked case you can only specify harmonic trap frequency (`OMEGA`) but in the density-unlocked case you can impose different trapping potentials on each component (via `OMEGA1` and `OMEGA2`).

## Configuration for a simulation:
There are two possibly simulation config files: 1) `config_dens_lck.json` - for a density-locked mixture; 2) `config_dens_ulck.json` - for density-unlocked mixture. As mentioned before, copy the default versions of these files from the 'sph_mix_lhy/json_init' directory and rename them to remove the `default.` prefix. These config files give you some flexibility in what to look at with these simulations
### Density-locked:
* The above experimental parameters for defining the mixture
* dr - the grid spacing
* Nr - the number of grid points (hence this also changes the box size, Lr = Nr\*dr)
* BC_TYPE - choice of boundary conditions
    * 0 - Neumann boundary conditions at both r = 0 and r = Lr
    * 1 - Neumann boundary conditions at r = 0 and Dirichlet at r = Lr
* DENS_LCK - is by default set to 1 in the config_dens_lck.json file
    * 0 - Unlocked mixture
    * 1 - Locked mixture
* DT_COEF - the size of the time-step set by dt = DT_COEFx\*dr^2)
* IM_T_STEPS - the number of imaginary time steps
* RE_T_STEPS - the number of real time steps
* T_SAVE - the number of steps between data saving
* GAUSS_SIGMA - the width of the initial condition
* OMEGA - the trapping potential strength
* INIT_TYPE - the form of the initial wavefunction guess
    * 'GAUSS' - Gaussian wavefunction profile
    * 'S_GAUSS' - Super-Gaussian wavefunction profile
* BREATH - sets whether to trigger a breathing mode for real time
    * 0 - No triggering of breathing mode
    * 1 - Trigger breathing mode using phase imprint exp(i\lambda r^2)
* ABSORB_BC - sets whether to turn on absorbing boundary conditions
    * 0 - Do not apply absorbing boundary conditions
    * 1 - Apply absorbing boundary conditions
* ABS_HEIGHT - sets height of tanh(r) function for absorbing boundary conditions
* ABS_SLOPE - sets the slope of tanh(r) function for absorbing boundary conditions
* ABS_POS - sets the position of tanh(r) function for absorbing boundary conditions
### Density-unlocked:
As above but with some slight differences, such as:
* OMEGA1,OMEGA2 - two trapping potentials, one for each component 
* INIT_TYPE1,INIT_TYPE2 - two initial guess functions for the wavefunction (for an imbalanced mixture it is often favourable to set the minority component to 'S_GAUSS' and the majority component to 'NON_ZERO_TAIL'
    * 'GAUSS' - Gaussian wavefunction profile
    * 'S_GAUSS' - Super-Gaussian wavefunction profile
    * 'NON_ZERO_TAIL' - Super-Gaussian + a small constant density mimicing the unbound density of an imbalanced mixture
* BREATH1,BREATH2 - two switches for triggering a breathing mode in real time
    * 0 - No triggering of breathing mode in the ith - component
    * 1 - Trigger breathing mode using phase imprint exp(i\lambda r^2) in the ith - component
* ABS_COMP - this flag allows for choosing with component to apply absorbing boundary conditions to
    * 'BOTH' - apply absorbing boundary conditions to both components
    * 'FIRST' - apply absorbing boundary conditions to component 1
    * 'SECOND' - apply absorbing boundary conditions to component 2
* T_DEPEN_POT - (soon to be removed) this switch allows for varying the absorbing boundary conditions in real time
    * 0 - Do not vary absorbing boundary conditions in real time
    * 1 - Turn off absorbing boundary conditions at some point in real time (this is vague as this feature is due to be removed)

## Simulating mixture on a single process
In order to run this a single process on a single process, the following command is needed:
`python3 running.py -wp <path-to-data-dir> -rp <sim_config_file>`
where `<sim_config_file>` is either `config_dens_lck.json` or `config_dens_ulck.json` 
