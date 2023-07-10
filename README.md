# AMM7_surge

[![DOI](https://zenodo.org/badge/293564924.svg)](https://zenodo.org/badge/latestdoi/293564924)

A 7km resolution Atlantic Margin Model 2D surge configuration.

The configuration is based on v4.0.4 of the [NEMO](https://www.nemo-ocean.eu) modelling framework.

The configuration recipe has be written for the [ARCHER2](https://www.archer2.ac.uk) UK National Supercomputing Service. 

The recipe describes how to configure and run a tides-only example.

---

## Repository File Hierarchy

### DOCS

A recipe on how to build and run AMM7_surge model.

### EXP_tideonly

An experiment directory for a (FES2014) tide only demonstration simulation.

### INPUTS

Store for external forcing files (e.g. tides, meteorology) and domain configuration file. Also store for boundary condition setup namelist file. *The domain configuration file can be downloaded from elsewhere. The tidal boundaries can be generated from this recipe or downloaded elsewhere*.

### ARCH

Store for architecture build files.

### MY_SRC

Store for FORTRAN modification to NEMO checkout from NEMO repository.


---

## Setting up AMM7 surge model

To run the AMM7 surge model follow the [AMM7_surge recipe](docs/AMM7_SURGE_build_and_run.rst). This recipe forms a template of how to obtain, compile and run the code.

As noted additional files are required to run the code. To run a full surge model meteorological forcing is required. For simplicity this demonstration simulation is configured to run without meteorological forcing (otherwise requiring sea level pressure and 10m winds). Tidal boundary conditions are also required - these can be generated following this [recipe](docs/generate_tidal_boundaries.rst) from the [docs](docs) folder, or downloaded elsewhere. Finally a domain configuration file is required - this can be generated following NEMO guidelines or downloaded elsewhere.
