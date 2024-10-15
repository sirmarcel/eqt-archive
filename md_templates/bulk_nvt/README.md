This runs MTS NVT MD for bulk water.

- `init-512_h2o.xyz`: Starting geometry (512 water molecules)
- `input-mts.xml`: i-PI input file
- `run.sh`: Example bash script on how to execute a run

MD is run using i-PI (https://github.com/i-pi/i-pi/) using the PET (https://github.com/spozdn/pet) driver. The `run.sh` script expects the model to be stored in `../model` with the state dict named `model_state_dict`. This can be adapted as needed.

The `polar-fes.py` script generates the orientational distribution, `tcorrelate.py` and `correlate.py` the dynamic and static correlations (respectively). 

Note: Similar to the molecule NVT case, production runs were not run sequentially, but in parallel, differing initial geometry and seeds. A total of 5.05ns for the unaugmented and augmented models were run each (so 10ns of total MD for the experiment). Runs were combined taking into account the analysis block size.