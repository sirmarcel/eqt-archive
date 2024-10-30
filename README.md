# Data Archive for "Probing the effects of broken symmetries in machine learning"

This collect the input files, scripts, and data for the [publication](https://doi.org/10.1088/2632-2153/ad86a0):

```
Title: Probing the effects of broken symmetries in machine learning
Authors: Marcel F. Langer, Sergey N. Pozdnyakov, and Michele Ceriotti
in Mach. Learn.: Sci. Technol. 5 04LT01
DOI: 10.1088/2632-2153/ad86a0
Preprint: arxiv:2406.17747 (2024)
```

The data in this archive is stored in the Materials Cloud Archive under [DOI:10.24435/materialscloud:kz-3b](https://doi.org/10.24435/materialscloud:kz-3b), and mirrored at [github.com/sirmarcel/eqt-archive](https://github.com/sirmarcel/eqt-archive).

***

The following items are contained within:

- `equivator-figs.ipynb`: Notebook to generate all figures in the publication (also contains some additional analysis used in the manuscript)
- `figures/`: Figures from the publication
- `data/`: The raw data used to generate the figures
- `md_templates/`: Input files for (gas, liquid) MD runs
- `scripts/`: Various Python scripts to process MD runs into the data used for publication
- `ice/`: Experiment for water ice
- `model/`: Contains the PET model used for this work

Additional `README.md` files are supplied to give more detailed explanations.

The main tools used in this work are the PET model, which is implemented in this repository: https://github.com/spozdn/pet, and the I-Pi code, obtained from https://github.com/i-pi/i-pi/ (`@99fcf1`). Experiments used a custom branch of PET, which can be found under https://github.com/lab-cosmo/pet/tree/neighbors_convert_cpp (`@e2631cc`) -- the functionality will be merged into main PET eventually.


***

To download the repository:

- On GitHub, you simply can run `git clone`, provided that [`git-lfs`](https://git-lfs.com) is enabled.
- On Materials Cloud Archive, the entire contents listed above are provided as a `.zip` file.

The files in this repository are licensed under the [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode) license.