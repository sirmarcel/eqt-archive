import ase 
import ase.io, ase.neighborlist
import numpy as np
import scipy as sp
import argparse
import matplotlib.pyplot as plt
from tqdm import tqdm
from ipi import read_output, read_trajectory

def autocorrelate(x, dt=1, xbar=None, normalize=True):
    """ Computes the autocorrelation function of a trajectory. 
    It can be given the exact average as a parameter"""
    if xbar is None:
        xbar = x.mean()
    acf = np.correlate(x - xbar, x-xbar, mode='same')
    if normalize:
        acf /= ((x-xbar)**2).sum()
    
    # fourier transform, with windowing
    ftacf = np.fft.hfft( (acf*np.blackman(len(acf)) )[len(x)//2:]) [:len(x)//2+1]
    print(len(acf), len(ftacf))
    
    omega = (
        np.arange(len(ftacf))
        / len(ftacf)
        * (np.pi / dt)
    )
    
    return np.arange(len(x)//2+1)*dt, acf[len(x)//2:], omega, ftacf

def pots_acf(filename, dt=1, progress=(lambda x: x)):
    """
    Computes the time correlation function and its FT, for both
    the potential and the difference between baseline and averaged 
    potentials.
    
    param filename: string, filename of the delta trajectory 
    param dt: float, time step of the delta trajectory
    """
    
    data, desc = read_output(filename)
    acf_pot = autocorrelate(data["pot_component_raw(0)"], dt=dt, normalize=False)   
    acf_delta = autocorrelate(data["pot_component_raw(1)"]-data["pot_component_raw(0)"], dt=dt, normalize=False)
    return dict(potential=acf_pot, delta=acf_delta)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Computes the pair correlation function and the orientational correlation
function, resolved into a longitudinal and transverse components.""")
    
    parser.add_argument("filename", type=str, help="The name of the file to process")
    parser.add_argument("--dt", type=float, default=1.0, help="The trajectory time step")    
    
    args = parser.parse_args()
    acfs = pots_acf(args.filename, dt=args.dt, progress=tqdm)

    # prints out the grid file
    np.savetxt("acf.dat", np.hstack([
                 acfs["potential"][0][:,np.newaxis], acfs["potential"][1][:,np.newaxis], acfs["delta"][1][:,np.newaxis] 
                 ]),
                 header = "   t     c_VV     c_DD")
                 
    np.savetxt("ftacf.dat", np.hstack([
                 acfs["potential"][2][:,np.newaxis], acfs["potential"][3][:,np.newaxis], acfs["delta"][3][:,np.newaxis] 
                 ]),
                 header = "   omega     chat_VV     chat_DD")
    
    # makes a figure for quick inspection
    fig, axs = plt.subplots(2,1,figsize=(5,3), constrained_layout=True)
    axs[0].plot(acfs["potential"][0], acfs["potential"][1], 'k')
    axs[0].plot(acfs["delta"][0], acfs["delta"][1], 'r')    
    axs[0].set_xlabel(r"$t$ / fs")
    axs[0].set_ylabel(r"$c(t)$")
    axs[1].semilogy(acfs["potential"][2], acfs["potential"][3], 'k')
    axs[1].semilogy(acfs["delta"][3], acfs["delta"][3], 'r')              
    axs[1].set_xlabel(r"$\omega$ / fs^{-1}")
    axs[1].set_ylabel(r"$\hat{c}(\omega)$")
    
    fig.savefig("acf.png")

