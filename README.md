This repository contains data and code for Goodman et al. (2022), **Shifting fish distributions impact predation intensity in a sub-Arctic ecosystem**, in review in the journal *Ecography*.

**Authors:** Maurice C. Goodman , Gemma Carroll, Stephanie Brodie, Arnaud Grüss, James T. Thorson, Stan  Kotwicki, Kirstin Holsman, Becca Selden, Elliott L. Hazen1, & Giulio A. De Leo

The relevant R scripts (contained in :file_folder: `src` ), in the order in which they should be sourced, are: 

- :file_folder: `predation_models` - this folder contains [VAST](https://github.com/James-Thorson-NOAA/VAST) models developed by Arnaud Grüss to jointly model predator biomass catch rate and stomach contents data (as prey-mass-per-predator-biomass), in order to derive estimates of predation across space and time. [See Grüss et al. (2020)](https://onlinelibrary.wiley.com/doi/full/10.1111/faf.12457). These scripts contain additional code to couch the model code in this project, sample from model predictive distributions, and plot and save model output. These scripts should be run first, in any order.
- `VAST-delta-models` - this script runs Bernoulli-Gamma delta models to estimate spatiotemporal biomass catch rate distributions for juvenile pollock and the predators examined in this study.
- `overlap_predation.R` - this script uses output from the :file_folder: `predation_models` and `VAST_delta-models.R` scripts to total and relative predation estimates against annual overlap metrics for each predator, and to produce plots used in the manuscript.
- `arrowtooth_example_plot.R` - this script produces figure 1.

The scripts `VAST_functions.R` and `range_overlap_functions.R`, which contain functions used and sourced by other scripts, do not need to be sourced themselves.
