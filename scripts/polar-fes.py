import ase 
import ase.io
import numpy as np
import scipy as sp
import argparse
import matplotlib.pyplot as plt
from tqdm import tqdm

def polar_kde(tp, gt, gp, bw=0.1):
    """ Computes a KDE for polar coordinates. 
    
    param tp: np.array((i,2)) array of (theta, phi) values for the samples
    param gt: np.array(ntheta) grid of cos(theta) values
    param gp: np.array(nphi) grid of phi values
    param bw: float bandwidth
    
    returns the grid of properly-symmetrized Gaussian values, accumulated over the samples
    """
    
    # computes the theta grid. this should have hard walls at [-1,1] which we implement with folded Gaussians
    gx = np.exp(-0.5*((tp[:,0][:,np.newaxis]-gt[np.newaxis,:]) / bw)**2)
    gx += np.exp(-0.5*((2-tp[:,0][:,np.newaxis]-gt[np.newaxis,:]) / bw)**2)
    gx += np.exp(-0.5*((-2-tp[:,0][:,np.newaxis]-gt[np.newaxis,:]) / bw)**2)
    gx /= np.sqrt(2*np.pi)*bw
    
    # computes the phi grid. this is periodic in [0, 2*pi] which we implement with one periodic image
    gy = np.exp(-0.5*((tp[:,1][:,np.newaxis]-gp[np.newaxis,:]) / bw)**2)
    gy += np.exp(-0.5*((2*np.pi+tp[:,1][:,np.newaxis]-gp[np.newaxis,:]) / bw)**2)
    gy += np.exp(-0.5*((-2*np.pi+tp[:,1][:,np.newaxis]-gp[np.newaxis,:]) / bw)**2)
    gy /= np.sqrt(2*np.pi)*bw
    
    return np.einsum("it,ip->tp", gx, gy)

def orientational_distribution(filename, kde_size=100, bandwidth=0.1, block_size=200, progress=(lambda x: x)):
    """
    Computes the distribution of the absolute orientation of water molecules in an ASE-readable file.
    The distribution is computed over the cosine of the polar angle and the azimuthal angle, and 
    divided by 4pi so that in a uniformly-distributed system it'll have a value of 1 everywhere.
    
    param filename: string, filename
    param kde_size: int, the grid will be (kde_size+1, 2*kde_size+1) along the (cos(theta), phi) directions
    param bandwidth: the bandwidth for the density estimation
    """
    
    grid = np.zeros(shape=(kde_size+1, 2*kde_size+1))
    grid_costheta = np.linspace(-1, 1, kde_size+1)
    grid_phi = np.linspace(np.pi*0.5/kde_size, 2*np.pi*(1-0.5/kde_size), 2*kde_size+1)    
    frame = ase.io.read(filename, "0")
    thetaphi = np.zeros(shape=(len(frame)//3, 2))
    nsamples = 0
    blocks = []
    for i, frame in progress(enumerate(ase.io.iread(filename))):
    
        # computes in one go the orientation vector for every water molecule (mid-point between H minus O)
        u = 0.5*(frame.positions[1::3]+frame.positions[2::3])-frame.positions[0::3]
        u /= np.sqrt((u**2).sum(axis=1))[:, np.newaxis]
        
        # converts in cos theta, phi
        thetaphi[:,0] = u[:,2]
        thetaphi[:,1] = np.arctan2(u[:,1], u[:,0])
        
        # computes the KDE for this frame and accumulates
        grid += polar_kde(thetaphi, grid_costheta, grid_phi, bandwidth)
        nsamples += len(thetaphi)
        
        # store blocks to compute uncertainty on the distribution
        if (i+1)%block_size == 0:
            blocks.append(grid*4*np.pi/nsamples)
            grid *= 0
            nsamples = 0
            
    if len(blocks)==0:  # if there's just one block, store it for further processing
        blocks.append(grid*4*np.pi/nsamples)
        
    blocks = np.asarray(blocks)
    print(blocks.shape)
    
    return (grid_phi, grid_costheta, blocks.mean(axis=0), blocks.std(axis=0)/np.sqrt(blocks.shape[0]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Computes the orientational distributions in a file corresponding to a MD trajectory of water. Lots of hardcoded assumptions: molecules should be 'intact' (i.e. not split over periodic boundaries), file should be in xyz or ase readable format, and atoms should be listed as OHH groups where each group is one water molecule.")
    
    parser.add_argument("filename", type=str, help="The name of the file to process")
    parser.add_argument("--bw", type=float, default=0.1,
                        help="Bandwidth for the density estimation (default: 0.1)")
    parser.add_argument("--ngrid", type=int, default=100,
                        help="Number of grid points for the kde (default: 100)")
    parser.add_argument("--block", type=int, default=100,
                        help="Blocking analysis size (default: 100)")                        
    parser.add_argument("--prefix", type=str, default="oc",
                        help="Prefix for the output files (default: 'oc')")

    args = parser.parse_args()
    grids = orientational_distribution(args.filename, args.ngrid, 
                                args.bw, args.block, progress=tqdm)

    # prints out the grid file
    np.savetxt(args.prefix+".dat", grids[2])    
    
    # makes a figure for quick inspection
    fig, axs = plt.subplots(2,1,figsize=(5,6), sharex=True, constrained_layout=True)
    cp0 = axs[0].contourf(grids[0], grids[1], grids[2], levels=20, cmap="viridis")
    cp1 = axs[1].contourf(grids[0], grids[1], grids[3], levels=20, cmap="viridis")
    axs[1].set_xlabel(r"$\phi$")
    axs[1].set_ylabel(r"$\cos\theta$")
    axs[0].set_ylabel(r"$\cos\theta$")
    
    fig.colorbar(cp0, ax=axs[0], orientation='vertical', label=r"$P(\theta,\phi)$")
    fig.colorbar(cp1, ax=axs[1], orientation='vertical', label=r"$\Delta P(\theta,\phi)$")   
    fig.savefig(args.prefix+".png")

