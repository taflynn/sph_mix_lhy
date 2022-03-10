## Import packages
import numpy as np
import matplotlib.pyplot as plt
import json
import matplotlib.animation as animation
import argparse


def run_ulck_process(dirarg):
    fname = 'config_dens_ulck.json'
    ## Read in data
    f = open('../data/' + dirarg + '/' + fname,"r")
    setup = json.loads(f.read())
    f.close()

    # grid spacing
    dr = setup['dr']

    # main data read in from NumPy files
    r = np.load('../data/' + dirarg + '/r_array.npy')
    t_real = np.load('../data/' + dirarg + '/real_t_array.npy')
    psi1 = np.load('../data/' + dirarg + '/real_spacetime_wav1.npy')
    psi2 = np.load('../data/' + dirarg + '/real_spacetime_wav2.npy')

    ## Plotting average droplet width as a function of time
    width1 = np.zeros((len(t_real)))
    width2 = np.zeros((len(t_real)))
    for i in range(0,len(t_real)):
        width1[i] = 4*np.pi*np.trapz(r**4*np.abs(psi1[:,i])**2,dx=dr)
        width2[i] = 4*np.pi*np.trapz(r**4*np.abs(psi2[:,i])**2,dx=dr)    
    plt.figure(figsize=(8,6))
    plt.plot(t_real,width1,t_real,width2,'--')
    plt.xlim((np.min(t_real),np.max(t_real)))
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\langle r^2 \rangle$')
    plt.legend((r'$\psi_1$',r'$\psi_2$'))
    plt.savefig('../data/' + dirarg + "/width.png",dpi=400)
    plt.close()

    ## Plotting central droplet density as function of time
    plt.plot(t_real,np.abs(psi1[1,:])**2)
    plt.plot(t_real,np.abs(psi2[1,:])**2)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$n_i^{(0)}$')
    plt.xlim((np.min(t_real),np.max(t_real)))
    plt.legend((r'$n_1^{(0)}$',r'$n_2^{(0)}$'))
    plt.savefig('../data/' + dirarg + "/centre_dens.png",dpi=400)
    plt.close()

    ## Animate the two components in time
    fig, ax = plt.subplots()
    line1, = ax.plot(r,np.abs(psi1[:,0])**2)
    line2, = ax.plot(r,np.abs(psi2[:,0])**2)

    def animate(i):
        line1.set_ydata(np.abs(psi1[:,i])**2)
        line2.set_ydata(np.abs(psi2[:,i])**2)
        plt.title('t = '+str(t_real[i]))
        plt.xlim((0,16))
        plt.ylim((0,None))
        plt.xlabel(r'$r$')
        plt.ylabel(r'$|\psi_{i}|^2$')
        plt.legend((r'$|\psi_1|^2$',r'$|\psi_2|^2$'))
        return line1,line2

    ani = animation.FuncAnimation(
        fig, animate, interval=10, blit=True, save_count=len(t_real))

    writer = animation.FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    ani.save('../data/' + dirarg + "/movie.avi", writer=writer,dpi=400)

    ## Animate the deviations from ground state (1st component)
    fig, ax = plt.subplots()
    line, = ax.plot(r,np.abs(psi1[:,0])**2 - np.abs(psi1[:,0])**2,linewidth=1)

    def animate(i):
        line.set_ydata(np.abs(np.abs(psi1[:,i])**2 - np.abs(psi1[:,0])**2))
        plt.yscale("log")
        plt.title('t = '+str(t_real[i]))
        plt.xlim((0,np.max(r)))
        plt.ylim((1e-16,1e-2))
        plt.xlabel(r'$r$')
        plt.ylabel(r'$|n_{1} - n_1^{0}|$')
        return line

    ani = animation.FuncAnimation(
        fig, animate, interval=10, blit=False, save_count=len(t_real))

    writer = animation.FFMpegWriter(fps=15, metadata=dict(artist='Me'),bitrate=1800)
    ani.save('../data/' + dirarg + "/dens1_dev_movie.avi", writer=writer,dpi=400)

    ## Animate the deviations from ground state (2nd component)
    fig, ax = plt.subplots()
    line, = ax.plot(r,np.abs(psi2[:,0])**2 - np.abs(psi2[:,0])**2)

    def animate(i):
        line.set_ydata(np.abs(np.abs(psi2[:,i])**2 - np.abs(psi2[:,0])**2))
        plt.yscale("log")
        plt.title('t = '+str(t_real[i]))
        plt.xlim((0,np.max(r)))
        plt.ylim((1e-16,1e-2))
        plt.xlabel(r'$r$')
        plt.ylabel(r'$|n_{2} - n_2^{0}|$')
        return line

    ani = animation.FuncAnimation(
        fig, animate, interval=10, blit=False, save_count=len(t_real))

    writer = animation.FFMpegWriter(fps=15, metadata=dict(artist='Me'),bitrate=1800)
    ani.save('../data/' + dirarg + "/dens2_dev_movie.avi", writer=writer,dpi=400)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Process data of density-unlocked mixture simulation')
    parser.add_argument('--read_path','-rp',
            dest = 'READ_PATH',
            type = str,
            required = True,
            nargs = 1)
    args = parser.parse_args()
    run_ulck_process(args.READ_PATH[0])
