import ase 
import ase.io, ase.neighborlist
import numpy as np
import scipy as sp
import argparse
import matplotlib.pyplot as plt
from tqdm import tqdm
   
def correlations(filename, block_size=200, dt=1, progress=(lambda x: x)):
    """
    Computes dynamical correlation functions - MSD and orientational correlations.
    
    param filename: string, filename
    param block_size: int, maximum time lag considered
    
    """
    
    # init grids
    msd = np.zeros(shape=(block_size))
    dipc = np.zeros(shape=(block_size))
    
    # geometric scaling factor for the g(r)
    frame = ase.io.read(filename, "0")
    
    nsamples = 0
    # computes multiple blocks so we can estimate statistical errors
    block_x = []
    block_u = []
    block_msd = []
    block_dipc = []
        
    nframes = 0
    
    for iframe, frame in progress(enumerate(ase.io.iread(filename))):
    
        # computes in one go the orientation vector for every water molecule (mid-point between H minus O)
        u = 0.5*(frame.positions[1::3]+frame.positions[2::3])-frame.positions[0::3]
        u /= np.sqrt((u**2).sum(axis=1))[:, np.newaxis]
        
        # we only need oxygens to make correlations
        del frame[np.where(frame.symbols=="H")[0]]

        x = frame.positions
        
        block_u.append(u)
        block_x.append(x)
                             
        # store blocks to compute uncertainty on the distribution
        if (iframe+1)%block_size == 0:
            dx = np.asarray(block_x) - block_x[0][np.newaxis,:,:]
            msd = (dx**2).sum(axis=-1).mean(axis=-1)
            du = np.asarray(block_u) * block_u[0][np.newaxis,:,:]
            dipc = du.sum(axis=-1).mean(axis=-1)
            
            block_msd.append(msd)
            block_dipc.append(dipc)
            block_u = []
            block_x = []
            
    if len(block_msd)==0:  # if there's just one block, store it for further processing
        dx = np.asarray(block_x) - block_x[0][np.newaxis,:,:]
        msd = (dx**2).sum(axis=-1).mean(axis=-1)
        du = np.asarray(block_u) * block_u[0][np.newaxis,:,:]
        dipc = du.sum(axis=-1).mean(axis=-1)
        
        block_msd.append(msd)
        block_dipc.append(dipc)
        
    
    block_msd = np.asarray(block_msd)
    block_dipc = np.asarray(block_dipc)

    return (
            dt*np.arange(block_msd.shape[1]),
            block_msd.mean(axis=0), block_msd.std(axis=0)/np.sqrt(block_msd.shape[0]),
            block_dipc.mean(axis=0), block_dipc.std(axis=0)/np.sqrt(block_dipc.shape[0]),
            block_dipc.mean(axis=0), block_dipc.std(axis=0)/np.sqrt(block_dipc.shape[0])
            )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Computes the pair correlation function and the orientational correlation
function, resolved into a longitudinal and transverse components.""")
    
    parser.add_argument("filename", type=str, help="The name of the file to process")
    parser.add_argument("--dt", type=float, default=10,
                        help="Frame stride (default: 10)")
    parser.add_argument("--block", type=int, default=100,
                        help="Blocking analysis size (default: 100)")                        
    parser.add_argument("--prefix", type=str, default="cft",
                        help="Prefix for the output files (default: 'cft')")

    args = parser.parse_args()
    grids = correlations(args.filename, block_size=args.block, dt=args.dt, progress=tqdm)

    # prints out the grid file
    np.savetxt(args.prefix+".dat", np.vstack(grids).T, 
            header = "  t     msd(t)    dmsd(t)     duu(t)    dcuu(t)")
    
    # makes a figure for quick inspection
    fig, axs = plt.subplots(2,1,figsize=(5,6), sharex=True, constrained_layout=True)
    axs[0].plot(grids[0], grids[1],'k')
    axs[1].plot(grids[0], grids[3], 'k')
    axs[1].set_xlabel(r"$t$ / fs")
    axs[1].set_ylabel(r"$c_{uu}(r)$")
    axs[0].set_ylabel(r"$MSD / \AA^2$")
    
    fig.savefig(args.prefix+".png")

