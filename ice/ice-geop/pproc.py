#!/usr/bin/python

import numpy as np
from ipi import read_output, read_trajectory

fmse = 0.0
xmse = 0.0
nat = 0
emol_pet = []
emol_aug = []
for i in range(2,11):
    prefix=f"ice-geop-{i:02d}"
    print(f"\n------------------------------\n Analyzing {prefix}")
    data_pet, desc = read_output(f"{prefix}.out")
    data_aug, desc = read_output(f"{prefix}_aug.out")    
    force_pet = read_trajectory(f"{prefix}.for_0.xyz")
    pos_pet = read_trajectory(f"{prefix}.pos_0.extxyz")
    force_aug = read_trajectory(f"{prefix}_aug.for_0.xyz")
    pos_aug = read_trajectory(f"{prefix}_aug.pos_0.extxyz")
    print("|F(init)|" , np.sqrt((force_aug[0].arrays["forces"]**2).sum()/len(pos_pet[0])))
    print("Delta E (init): ",  (data_pet["potential"][0] - data_aug["potential"][0])/ len(pos_pet[0]) )
    print("Delta Delta E: ",  (data_pet["potential"][-1] - data_pet["potential"][0] - (data_aug["potential"][-1] -  data_aug["potential"][0]))/ len(pos_pet[0]) )
    print("DeltaCOM", (pos_aug[-1].positions-pos_pet[-1].positions).sum(axis=0))
    print("Delta x", np.sqrt(np.sum((pos_aug[-2].positions-pos_aug[-1].positions)**2)/len(pos_pet[0])))
    emol_pet.append(data_pet["potential"][-1]*3/len(pos_pet[0]))
    emol_aug.append(data_aug["potential"][-1]*3/len(pos_pet[0]))
    fmse += np.sum((force_aug[0].arrays["forces"]-force_pet[0].arrays["forces"])**2)
    xmse += np.sum((pos_aug[-1].positions-pos_pet[-1].positions)**2)
    nat += len(pos_pet[0])

emol_pet = np.array(emol_pet)*27.211386
emol_aug = np.array(emol_aug)*27.211386
print("E per mol. ", emol_aug)
print("     std", np.std(emol_aug))
print(f"RMSE on relative energies: {np.sqrt( np.sum( ((emol_pet - emol_pet.mean()) - (emol_aug-emol_aug.mean()) )**2  ) )  } ")
print(f"RMSE force error on initial structures: {np.sqrt(fmse*3/nat)} eV/Å")
print(f"RMSE position on final structures: {np.sqrt(xmse*3/nat)} Å")
