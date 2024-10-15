# `data/`

This contains all the raw data for the figures, including the SI. In particular:

```
# H2O molecule
# ... NVE (Figure 1)
data/nve_pet3.dat
data/nve_pet3-aug2i.dat
data/nve_pet3-delta_powerspec.dat
data/mm-*.dat

# ... NVT (Figure 1)
data/pfes_pet3.dat
data/pfes_pet3-2i.dat

# bulk water (NVT)
# ... orientational "free energy" (Figure 2, SI Figure 3)
data/bulk_pfes_pet3.dat
data/bulk_pfes_pet3-2i.dat
data/bulk_pfes_pet3rnd.dat

# ... static correlations (Figure 2, SI Figure 3)
data/bulk_corr_pet3.dat
data/bulk_corr_pet3-2i.dat
data/bulk_corr_pet3rnd.dat

# ... dynamical correlations (Figure 3, SI Figure 4)
data/bulk_tcorr_pet3.dat
data/bulk_tcorr_pet3-2i.dat
data/bulk_tcorr_pet3rnd.dat

# ice (no figure)
data/base-nrg.dat
data/pet_3_2_w_100_lr_1e-5-nrg.dat

# SI: error convergence (Figure 1)
data/test_errors/*
data/nvt_errors/*
data/mol_errors/*

# SI: MD stability (Figure 2)
data/nvt.out
data/mts-md.out
data/nvt-rndaug.out
```