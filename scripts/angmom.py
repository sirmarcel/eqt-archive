import ase 
import ase.io
import numpy as np
import scipy as sp
import argparse
import matplotlib.pyplot as plt
from ipi import read_output, read_trajectory
from tqdm import tqdm


def conserved_stuff(prefix, progress=(lambda x: x)):
    """
    Computes (or extracts from the trajectory) conserved quantities for a water molecule in vacuum.
    Everything is hardcoded, because this is the way
    
    param predix: string, file prefix
    """

    # trajectory outputs
    out, info = read_output(f"{prefix}.out")

    # angular momentum (assumes linear momentum is zero!)
    t_pos = ase.io.read(f"{prefix}.pos_0.ase", ":", format="xyz")
    t_mom = read_trajectory(f"{prefix}.mom_0.xyz")
    t_for = read_trajectory(f"{prefix}.for_0.xyz")    
   
    angmom = np.zeros((len(t_pos), 3))
    torque = np.zeros((len(t_pos), 3))
    m = np.array([15.9994, 1.00794, 1.00794])
    for i, (fp, fm, ff) in enumerate(zip(t_pos, t_mom, t_for)):
       pos = fp.positions / 0.529177 # to a.u.
       # subtracts center of mass
       com = (pos[0]*15.9994 + pos[1]*1.00794 + pos[2]*1.00794)/m.sum()
       pos -= com
       mom = fm.positions  # i-PI stores momenta and forces in the position field of the xy, and in a.u.
       ffor = ff.arrays['forces']*0.019446904 # back to a.u. 
       angmom[i] = np.sum([np.cross(pos[i], mom[i]) for i in range(3)], axis=0)
       torque[i] = np.sum([np.cross(pos[i], ffor[i]) for i in range(3)], axis=0)
    
    # time in a.u. ...
    return dict(conserved=out["conserved"], angmom=angmom, torque=torque, time=out["time"])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Computes conserved energy and angular momentum for a water molecule. Everything is hardcoded!")
    
    parser.add_argument("prefix", type=str, help="The prefix of the i-PI run")

    args = parser.parse_args()
    stuff = conserved_stuff(args.prefix, progress=tqdm)
    # prints out the grid file, converting to "friendly" units - ps, eV, eV*ps
    np.savetxt(args.prefix+".dat", np.hstack([stuff["time"][:, np.newaxis]*2.4188843e-05, stuff["conserved"][:, np.newaxis]*27.211386, stuff["angmom"]*0.00065821196, stuff["torque"]*27.211386 ]) , 
               header = "t/ps   H/eV     Lx/(eV*ps)     Ly/(eV*ps)    Lz/(eV*ps)   taux/eV   tauy/eV   tauz /eV" )

    exit(0)
    
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

