# AMM7_surge

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10605585.svg)](https://doi.org/10.5281/zenodo.10605585)

A 7km resolution Atlantic Margin Model 2D surge configuration.

The configuration is based on v4.0.4 of the [NEMO](https://www.nemo-ocean.eu) modelling framework.

The configuration recipe has be written for the [ARCHER2](https://www.archer2.ac.uk) UK National Supercomputing Service. 

The recipe describes how to configure and run a tides-only and a tide+meteorolgy example.

---

## Repository File Hierarchy

### DOCS

A recipe on how to build and run AMM7_surge model.

### INPUTS

A store for external forcing files (e.g. tides, meteorology) and domain configuration file. Also store for boundary condition setup namelist file. *The domain configuration file can be downloaded from elsewhere. The tidal boundaries can be generated from this recipe or downloaded elsewhere*.

### ARCH

Store for architecture build files.

### NEMO_4.0.4_surge

Parent directory for source code and experiment set up files

### NEMO_4.0.4_surge/cfgs/AMM7_SURGE/MY_SRC

Store for FORTRAN modification to NEMO checkout from NEMO repository.

### NEMO_4.0.4_surge/cfgs/AMM7_SURGE/EXP_NOWIND_DEMO

An experiment directory for a tide only (FES2014) demonstration simulation.

### NEMO_4.0.4_surge/cfgs/AMM7_SURGE/EXP_ERA5_DEMO

An experiment directory for a tide (FES2014) + surface wind and sea level pressure forced (ERA5) demonstration simulation.

---

## Setting up AMM7 surge model

To run the AMM7 surge model follow the [AMM7_surge recipe](docs/AMM7_SURGE_build_and_run.rst). This recipe forms a template of how to obtain, compile and run the code.

In the demonstration examples are given how to run on ARCHER2 and also how to run in a Singularity container (which has been independently tested on a macbook and a linux server).

As noted additional files are required to run the code. These include: a domain file, boundary tides, sea level pressure and 10m winds. The are expected to be copied into the INPUTS folder. Example files are downloadable from ...


---


## Configuration Input Files

|  **Input** | **Download Location** |
|-------------- | -------------- |
| **Domain_cfg.nc** | https://gws-access.jasmin.ac.uk/public/jmmp/AMM7_surge/domain_cfg.nc |
| **Open ocean boundary coordinates.bdy.nc** | http://gws-access.jasmin.ac.uk/public/jmmp/AMM7_surge/coordinates.bdy.nc |
| **Bottom friction 2D scaling bfr_coef.nc** | http://gws-access.jasmin.ac.uk/public/jmmp/AMM7_surge/bfr_coef.nc |

Example extraction into INPUTS directory:

```
wget https://gws-access.jasmin.ac.uk/public/jmmp/AMM7_surge/domain_cfg.nc     -O INPUTS/amm7_surge_domain_cfg.nc
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7_surge/coordinates.bdy.nc -O INPUTS/coordinates.bdy.nc
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7_surge/bfr_coef.nc        -O INPUTS/bfr_coef.nc
```
---

## Sample Forcing Files

| **Forcing** | **Download Location** |
|-------------- | ------------------|
| **Surface boundary** | http://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/SBC/ |
| **Tide** | https://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/TIDE/FES/ |

Example extraction into INPUTS directory:

```
cd INPUTS/fluxes
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/SBC/ERA5_U10_y2017.nc .
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/SBC/ERA5_V10_y2017.nc .
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/SBC/ERA5_MSL_y2017.nc .
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/SBC/ERA5_LSM.nc       .
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/SBC/weights_era5_bicubic.nc .

cd INPUTS/bdydta
wget https://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/TIDE/FES/AMM7_surge_bdytide_rotT_M2_grid_U.nc .
wget https://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/TIDE/FES/AMM7_surge_bdytide_rotT_M2_grid_V.nc .
wget https://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/TIDE/FES/AMM7_surge_bdytide_rotT_M2_grid_T.nc .
wget https://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/TIDE/FES/AMM7_surge_bdytide_rotT_S2_grid_U.nc .
wget https://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/TIDE/FES/AMM7_surge_bdytide_rotT_S2_grid_V.nc .
wget https://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/TIDE/FES/AMM7_surge_bdytide_rotT_S2_grid_T.nc .
wget https://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/TIDE/FES/AMM7_surge_bdytide_rotT_K2_grid_U.nc .
wget https://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/TIDE/FES/AMM7_surge_bdytide_rotT_K2_grid_V.nc .
wget https://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/TIDE/FES/AMM7_surge_bdytide_rotT_K2_grid_T.nc .
```
