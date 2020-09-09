# AMM7_surge

A 7km resolution Atlantic Margin Model 2D surge configuration using the [NEMO](https://www.nemo-ocean.eu) modelling framework.  

This configuration recipe has be written with the [ARCHER](https://www.archer.ac.uk) HPC INTEL environment in mind.

---

## Setting up AMM7 surge model

See [AMM7_surge recipe](docs/AMM7_SURGE_build_and_run.rst)

This recipe forms a template of how to obtain and compile the code. Additional files are required to run the code.


---

## File Hierarchy


### EXP_tideonly

FES2014 tide only simulation

### INPUTS

Store for external forcing files and domain configuration file.

### ARCH

Store for architecture build files

### MY_SRC

Store for FORTRAN modification to NEMO checkout from NEMO repository.

### DOCS

Recipe on how to build and run AMM7_surge model.
