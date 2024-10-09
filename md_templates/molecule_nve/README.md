This runs NVE MD for a single water molecule to check angular momentum conservation for an "unconstrained" model.

`input-mol.xml` runs plain MD, `input-mol-delta.xml` additionally expects an "augmented", i.e., rotationally averaged model, in the background to check the delta between the two.

MD is run using i-pi (https://github.com/i-pi/i-pi/) using the PET (https://github.com/spozdn/pet) driver. The `run.sh` script expects the model to be stored in `../model` with the state dict named `model_state_dict`. This can be adapted as needed. Same for `run-delta.sh`

The `angmom.py` script is designed to check angular momentum.  `timecf.py` checks the correlations between augmented and unaugmented energies.