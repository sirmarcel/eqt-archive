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