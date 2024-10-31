# `model/`

This subfolder contains the production model used in the manuscript. (Internally named `pet_3`, earlier versions were used for testing.)

It contains the checkpoint output of the PET code, which consists of:

- `all_species.npy`: All species used during training
- `hypers_used.yaml`: All hyper-parameters for training (note that the model was trained in two stages, this contains the hypers for the second stage)
- `model_state_dict`: `torch` state dict for the model, renamed from `best_val_mae_forces_model_state_dict`, i.e., this contains the weights for the model with the best MAE (mean absolute error) for the forces
- `self_contributions.npy`: Per-element energy contributions (not used here)
- `summary.txt`: Summary metrics
- `note.txt`: Internal note about the provenance of the model (initialls trained with slow LR decay, then "polished" with rapid LR decay)

The input scripts in the `md_templates` folder expect this folder.

## Dataset

For training, we used the dataset from https://github.com/BingqingCheng/neural-network-potential-for-water-revPBE0-D3/tree/master/training-set, the splits are given in `permutation.npy` in this folder.

The splits used in this work can be obtained as:

```python
import numpy as np
import ase.io

structures = ase.io.read('dataset_1593.xyz', index = ':')

for struc in structures:
    struc.positions = struc.positions * 0.529177210903
    struc.cell = struc.cell * 0.529177210903

    struc.info['energy'] = struc.info['TotEnergy'] * 27.211386245988
    struc.info.pop('TotEnergy', None)

    struc.arrays['forces'] = struc.arrays['force'] * 51.42208619083232
    struc.arrays.pop('force', None)

permutation = np.load("permutation.npy")

structures = [structures[index] for index in permutation]
N_TEST = 190
N_VAL = 100
N_TRAIN = len(structures) - N_TEST - N_VAL

ase.io.write(f'water_train_{N_TRAIN}.xyz', structures[0:N_TRAIN])
ase.io.write(f'water_val.xyz', structures[N_TRAIN : N_TRAIN + N_VAL])
ase.io.write(f'water_test.xyz', structures[N_TRAIN + N_VAL:])

```