# BITS paper code
[![DOI](https://zenodo.org/badge/910080608.svg)](https://zenodo.org/badge/latestdoi/910080608)

This repository contains the `R` code used for obtaining the empirical results
presented in the manuscript that introduces BITS (Boosting Interaction Tree Stumps).
The `BITS` software package can be found in [another repository](https://github.com/michlau/BITS).

The code in this repository is structured as follows:

```
.
├── common
│   └── prep.R [Common code for data generation and splitting and method evaluation]
├── real-data-app
│   ├── dorothea_data.R [Preparing Dorothea data]
│   ├── dorothea_results.R [Processing Dorothea results]
│   ├── dorothea_run.R [Running method evaluation on Dorothea data]
│   ├── salia_results.R [Processing SALIA results]
│   └── salia_run.R [Running method evaluation on SALIA data]
└── simulation-study
    ├── discarded_terms.R [Evaluating how many terms can be discarded by BITS]
    ├── sim_results.R [Processing simulation study results]
    ├── sim_run.R [Running method evaluation on simulated data]
    ├── terms.R [Evaluating term identification performance]
    └── time_boosting.R [Measuring time needed for model fitting and prediction]
```


