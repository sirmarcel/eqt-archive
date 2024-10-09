# Data Archive for "Probing the effects of broken symmetries in machine learning"

This collect the input files, scripts, and data for the publication:

```
Title: Probing the effects of broken symmetries in machine learning
Authors: Marcel F. Langer, Sergey N. Pozdnyakov, and Michele Ceriotti
Publication: TBD
Preprint: arxiv:2406.17747 (2024)
```

The data in this archive is stored in the Materials Cloud under TBD, and mirrored at [github.com/sirmarcel/eqt-archive](https://github.com/sirmarcel/eqt-archive).

***

The following items are contained within:

- `equivator-figs.ipynb`: Notebook to generate all figures in the publication (also contains some additional analysis used in the manuscript)
- `figures/`: Figures from the publication
- `data/`: The raw data used to generate the figures
- `md_templates/`: Input files for (gas, liquid) MD runs
- `scripts/`: Various Python scripts to process MD runs into the data used for publication
- `ice/`: Experiment for water ice
- `models/`: Contains the checkpoint (called `pet_3` internally) used for this work

Additional `README.md` files are supplied to give more detailed explanations.

The main tools used in this work are the PET model, which is implemented in this repository: https://github.com/spozdn/pet, and the I-Pi code, obtained from https://github.com/i-pi/i-pi/ . Experiments used a custom branch of PET, which can be found under https://github.com/lab-cosmo/pet/tree/neighbors_convert_cpp -- the functionality will be merged into main PET eventually.
