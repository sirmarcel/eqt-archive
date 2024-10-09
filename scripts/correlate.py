import ase 
import ase.io, ase.neighborlist
import numpy as np
import scipy as sp
import argparse
import matplotlib.pyplot as plt
from tqdm import tqdm

def correlations(filename, kde_size=100, r_max=8, bandwidth=0.1, block_size=200, progress=(lambda x: x)):
    """
    Computes the distribution of the absolute orientation of water molecules in an ASE-readable file.
    The distribution is computed over the cosine of the polar angle and the azimuthal angle, and 
    divided by 4pi so that in a uniformly-distributed system it'll have a value of 1 everywhere.
    
    param filename: string, filename
    param kde_size: int, the grid will be (kde_size+1, 2*kde_size+1) along the (cos(theta), phi) directions
    param bandwidth: the bandwidth for the density estimation
    """
    
    # init grids
    grid_r = np.linspace(0, r_max, kde_size+1)
    grid_gr = np.zeros(shape=(kde_size+1))
    grid_dl = np.zeros(shape=(kde_size+1))
    grid_dt = np.zeros(shape=(kde_size+1))    
    
    # geometric scaling factor for the g(r)
    ir2 = 4 * np.pi * grid_r**2 * np.sqrt(2*np.pi*bandwidth**2)
    ir2[0] = 1 # avoids nan at the origin
    ir2 = 1/ir2 
    frame = ase.io.read(filename, "0")
    
    # nb - this builds a half neigbor list so some atoms might have no neighbors
    nl = ase.neighborlist.NeighborList([r_max/2]*(len(frame)//3), self_interaction=False)
    
    nsamples = 0
    # computes multiple blocks so we can estimate statistical errors
    blocks_gr = []
    blocks_dl = []
    blocks_dt = []    
    nframes = 0
    for iframe, frame in progress(enumerate(ase.io.iread(filename))):
    
        # computes in one go the orientation vector for every water molecule (mid-point between H minus O)
        u = 0.5*(frame.positions[1::3]+frame.positions[2::3])-frame.positions[0::3]
        u /= np.sqrt((u**2).sum(axis=1))[:, np.newaxis]
        
        # we only need oxygens to make correlations
        del frame[np.where(frame.symbols=="H")[0]]
        nl.update(frame)
        
        nrho = (0.5*len(frame)*(len(frame)-1)/frame.cell.volume)
                
        for i in range(len(frame)):
            indices, offsets = nl.get_neighbors(i)
            if len(indices) == 0:
                continue
            dr = frame.positions[indices] + offsets @ frame.get_cell() - frame.positions[i]
            r = np.sqrt((dr**2).sum(axis=1))
            
            gr = np.exp(-0.5*((r[:,np.newaxis]-grid_r[np.newaxis,:]) / bandwidth)**2)
            grid_gr += gr.sum(axis=0)
            
            # longitudinal and transverse dipole correlations
            rhat = dr/r[:,np.newaxis]
            uijl = (rhat@u[i])*(rhat*u[indices]).sum(axis=1)
            uijt = u[indices]@u[i] - uijl  # transverse component is full product minus the part projected on rhat
            
            dl = gr*uijl[:, np.newaxis]
            dt = gr*uijt[:, np.newaxis]            
            
            grid_dl += dl.sum(axis=0)            
            grid_dt += dt.sum(axis=0)
            
            nsamples += len(np.where(r<r_max)[0])
        
        # store blocks to compute uncertainty on the distribution
        if (iframe+1)%block_size == 0:
            blocks_gr.append(grid_gr*ir2/nrho/block_size)
            blocks_dl.append(grid_dl*ir2/nrho/block_size)
            blocks_dt.append(grid_dt*ir2/nrho/block_size)            
            grid_gr *= 0
            grid_dl *= 0            
            grid_dt *= 0                        
            nsamples = 0
        nframes = nframes + 1
            
    if len(blocks_gr)==0:  # if there's just one block, store it for further processing
        blocks_gr.append(grid_gr*ir2/nrho/nframes)
        blocks_dl.append(grid_dl*ir2/nrho/nframes)
        blocks_dt.append(grid_dt*ir2/nrho/nframes)
    
    blocks_gr = np.asarray(blocks_gr)
    blocks_dl = np.asarray(blocks_dl)
    blocks_dt = np.asarray(blocks_dt)            
    
    return (grid_r,
            blocks_gr.mean(axis=0), blocks_gr.std(axis=0)/np.sqrt(blocks_gr.shape[0]),
            blocks_dl.mean(axis=0), blocks_dl.std(axis=0)/np.sqrt(blocks_dl.shape[0]),
            blocks_dt.mean(axis=0), blocks_dt.std(axis=0)/np.sqrt(blocks_dt.shape[0])
            )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Computes the pair correlation function and the orientational correlation
function, resolved into a longitudinal and transverse components.""")
    
    parser.add_argument("filename", type=str, help="The name of the file to process")
    parser.add_argument("--rmax", type=float, default=8,
                        help="Maximum distance considered (default: 8)")
    parser.add_argument("--bw", type=float, default=0.1,
                        help="Bandwidth for the density estimation (default: 0.1)")
    parser.add_argument("--ngrid", type=int, default=200,
                        help="Number of grid points for the kde (default: 200)")
    parser.add_argument("--block", type=int, default=100,
                        help="Blocking analysis size (default: 100)")                        
    parser.add_argument("--prefix", type=str, default="gr",
                        help="Prefix for the output files (default: 'gr')")

    args = parser.parse_args()
    grids = correlations(args.filename, args.ngrid, args.rmax, 
                                args.bw, args.block, progress=tqdm)

    # prints out the grid file
    np.savetxt(args.prefix+".dat", np.vstack(grids).T, 
            header = "  r     g(r)    dg(r)     cLuu(r)    dcLuu(r)      cTuu(r)    dcTuu(r)")
    
    # makes a figure for quick inspection
    fig, axs = plt.subplots(2,1,figsize=(5,6), sharex=True, constrained_layout=True)
    axs[0].plot(grids[0], grids[1],'k')
    axs[1].plot(grids[0], grids[3]/grids[1], 'b', label="longitudinal")
    axs[1].plot(grids[0], grids[5]/grids[1], 'r', label="transverse")    
    axs[1].set_xlabel(r"$r$ / Å")
    axs[1].set_ylabel(r"$c_{uu}(r)$")
    axs[1].legend()
    axs[1].set_ylim(-0.1,0.5)
    axs[0].set_ylabel(r"$g(r)$")
    
    fig.savefig(args.prefix+".png")

