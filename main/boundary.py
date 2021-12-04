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
