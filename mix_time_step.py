def imag_time(m1,m2,a11,a22,a12,N,DENS_LCK,dr,Nr,dt):
    if DENS_LCK == 1:
        N_lck = dens_lck_params(m1,m2,a11,a22,a12)

        psi,tol,spacetime = rk4_dens_lck(r,psi,N_lck)

    elseif DENS_LCK == 0:
        if m1 == m2:
            [alpha,beta,eta] = dens_unlck_params(m1,m2,a11,a22,a12)
        else if m1 != m2:
            print('Density-unlocked, unequal masses is still being worked on!')
    return r,psi,psi1,psi2,tol,spacetime
