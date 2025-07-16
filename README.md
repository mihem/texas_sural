# Texas Sural Nerve Single-Cell RNA-seq Analysis

## Overview

This project analyzes snRNAseq data from sural nerves.

## Analysis Pipeline

The analysis is implemented using the [targets](https://docs.ropensci.org/targets/) R package for reproducible workflows. The pipeline consists of four main modules:


### Environment Management via `rix`
- Use `generate_env.R` to set up R environment.
```
Rscript generate_env.R
```
- Nix environment specified in `default.nix`
- Build the nix environment with `nix-build` and `nix-shell`

### Targets pipeline


```r
# Load targets
library(targets)

# Run entire pipeline
tar_make()

# Check pipeline status
tar_visnetwork()
```


### Key Outputs
- stored in `results/` directory