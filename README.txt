--- SPHERICAL BOSE-BOSE MIXTURE W/ LHY ---

(1) Intro:
-> Overview: These codes are to solve the Bose-Bose mixture, with the Lee-Huang-Yang (LHY) correction, in a spherically symmetric geometry (i.e. the density of the mixture is a function of radius only).

-> Time-stepping scheme: The code uses a 4th-order Runge-Kutta Scheme.

-> Boundary conditions: The code has options for boundary conditions in the 'boundary.py' file though the only current options are: 1. Neumann boundary conditions at both ends; 2. Neumann boundary condition at centre of the box w/ density set to zero at the far end of the box.

(2) Overview:

(3) Describing a mixture:
A Bose-Bose mixture is here defined by:
* a11 - intraspecies scattering length for component 1
* a12 - interspecies scattering length
* a22 - intraspecies scattering length for component 2
* m1 - atomic mass of component 1
* m2 - atomic mass of component 2
We also require an atom number for the system and this is where the first separation occurs. 

As mention the mixture can be density-locked (DENS_LCK = TRUE) or density-unlocked (DENS_LCK = FALSE). When in the density-locked regime, you need only specify the total number of atoms, N. However, when in the density-unlocked regime, you must specify the atom number for both components (i.e. N1 and N2).

One other important way to specify your system is the trapping potential. In the spherically symmetric case we are somewhat limited so the main choices are: 1. No potential; 2. Harmonic potential. The only subtlty is that in the density-locked case you can only specify harmonic trap frequency (OMEGA) but in the density-unlocked case you can impose different trapping potentials on each component (via OMEGA1 and OMEGA2).