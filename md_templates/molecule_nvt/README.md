This runs NVT MD for a single water molecule to check angular momentum conservation for an "unconstrained" model. We use multiple time stepping (MTS) to execute the augmented model every 10 timesteps.

MD is run using i-pi (https://github.com/i-pi/i-pi/) using the PET (https://github.com/spozdn/pet) driver. The `run.sh` script expects the model to be stored in `../model` with the state dict named `model_state_dict`. This can be adapted as needed.

The `polar-fes.py` script generates the orientational distribution.

Note: For practical reasons, production calculations were run in parallel (one long run of about 5ns, and shorter runs of ~0.5ns each), making sure to differ initial geometries and initialisation seed. Resulting trajectories were then stitched together (accounting for post-processing block size) to obtain the data reported in the paper. We checked correlation times to ensure that this yields reasonable results.
