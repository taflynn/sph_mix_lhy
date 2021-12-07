import numpy as np

def bc(x,BC_TYPE):
    """
    This function applies boundary conditions to the radial box. The function takes two inputs:
    (1) x - The array to apply the boundary conditions too (e.g. the wavefunction)
    (2) BC_TYPE - The type of boundary conditions
    
    There are only two options for boundary conditions:
    (0) BC_TYPE = 0: This corresponds to Neumann boundary conditions at both ends of the box, i.e.,
        we set the first derivative of the array equal to zero at both ends of the box
    (1) BC_TYPE = 1: This corresponds to a Neumann boundary condition at the inner end of the box
        and fixing the array to zero at the far end of the box
    """
    # Neumann boundary conditions at both ends
    if BC_TYPE == 0:
        x[0] = x[1]
        x[-1] = x[-2]
    # Neumann boundary condition at centre with zero at far end
    if BC_TYPE == 1:
        x[0] = x[1]
        x[-1] = 0
    return x

def absorb_bc_dens_lck(r,ABS_COEF):
    """
    This function calculates the imaginary potential needed to apply absorbing boundary conditions to the density-locked mixture.
    There function takes 2 inputs:
    (1) r - the spatial array
    (2) ABS_COEF - the coefficient of the imaginary potential

    The function then just outputs one potential, V.
    """
    V = 1.0j*ABS_COEF*(np.tanh(r - 0.9*np.max(r))+1)
    return V

def absorb_bc_dens_ulck(r,ABS_COEF,ABS_COMP):
    """
    This function calculates the imaginary potential needed to apply absorbing boundary conditions to the density-unlocked mixture.
    There function takes 3 inputs:
    (1) r - the spatial array
    (2) ABS_COEF - the coefficient of the imaginary potential
    (3) ABS_COMP - the component the boundary conditions should be administered too e.g.
        'BOTH' - apply conditions to both components
        'FIRST' - apply condition to 1st-component
        'SECOND' - apply condition to 2nd-component

    The function then just outputs two potentials, V1 and V2.
    """
    if ABS_COMP == 'BOTH':
        V1 = 1.0j*ABS_COEF*(np.tanh(r - 0.75*np.max(r))+1)
        V2 = 1.0j*ABS_COEF*(np.tanh(r - 0.75*np.max(r))+1)
    if ABS_COMP == 'FIRST':
        V1 = 1.0j*ABS_COEF*(np.tanh(r - 0.75*np.max(r))+1)
        V2 = np.zeros(len(r))
    if ABS_COMP == 'SECOND':
        V1 = np.zeros(len(r))
        V2 = 1.0j*ABS_COEF*(np.tanh(r - 0.75*np.max(r))+1)
    return V1,V2
