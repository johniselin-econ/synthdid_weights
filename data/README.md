This package contains the following data sets:

1. California proposition 99 (`california_prop99.csv`)
This data is from Abadie, Diamond, and Hainmueller (2010). The raw data is in MATLAB format from https://web.stanford.edu/~jhain/synthpage.html and is preprocessed here to a long panel, and saved as a `;` delimited CSV to work with R's `data()` function.

2. CPS (`CPS.csv`) and Penn World Table (`PENN.csv`) panels used in the placebo studies of Arkhangelsky et al. (2021).

Paper-specific analysis data (the ACA Medicaid expansion application, Borgschulte & Vogler 2020 replication) lives in `paper/data/`, not here, so it is not bundled into the package. Those files are produced by the scripts in `scripts/` (see `scripts/county_dataset_creation.R` etc.).
