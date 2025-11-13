# AMM7_surge

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10605585.svg)](https://doi.org/10.5281/zenodo.10605585)

A 7km resolution Atlantic Margin Model 2D surge configuration, based on v4.0.4 of the [NEMO](https://www.nemo-ocean.eu) modelling framework.

The configuration has example demonstrations written for the [ARCHER2](https://www.archer2.ac.uk) UK National Supercomputing Service, or to be run in a Singularity container.

This release includes a tutorial describing how to configure and run a tides-only and a tide+meteorolgy example.

This collaboration was facilitated by the [National Partnership for Ocean Prediction](https://oceanprediction.org/) and the [Joint Marine Modelling Programme](https://www.metoffice.gov.uk/research/approach/collaboration/joint-marine-modelling-programme).

---

## Repository File Hierarchy

### DOCS

A recipe on how to build and run AMM7_surge model on the NOC HPC cluster Anemone.

### INPUTS

A store for external forcing files (e.g. tides, meteorology) and domain configuration file. Also store for boundary condition setup namelist file. *The domain configuration file can be downloaded from elsewhere. The tidal boundaries can be generated from this recipe or downloaded elsewhere*.

### ARCH

Store for architecture build files.

### NEMO_4.0.4_surge

Parent directory for source code and experiment set up files

### NEMO_4.0.4_surge/cfgs/AMM7_SURGE/MY_SRC

Store for FORTRAN modification to NEMO checkout from NEMO repository.

### NEMO_4.0.4_surge/cfgs/AMM7_SURGE/EXP_ERA5_DEMO

An experiment directory for a tide (FES2014) + surface wind and sea level pressure forced (ERA5) demonstration simulation.

---

### Setting up AMM7 surge model

To run the AMM7 surge model on Anemone follow the recipe below.

```
ssh anemone

module purge
module load NEMO/prg-env
 
# Navigate to a suitable directory
# cd my_work_dir

export CONFIG=AMM7_SURGE
export WORK=$PWD
export WDIR=$WORK/$CONFIG
export INPUTS=$WDIR/INPUTS
export CDIR=$WDIR/NEMO_4.0.4_surge/cfgs # compile location is down here
export EXP=$CDIR/$CONFIG/EXP_ERA5_DEMO

# clone repo
git clone -b feature/v4.0.4 https://github.com/JMMP-Group/AMM7_surge.git $CONFIG

# cp arch file for compiling nemo
cd $CDIR/..
sed "s|XXX_XIOS_XXX|$WDIR/xios|g"  ../ARCH/ANEMONE/arch-anemone-ifort-impi.fcm > arch/NOC/arch-anemone-ifort-impi.fcm
cd $WDIR

# extract xios and compile
svn co https://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5@1964 xios
cd xios
cp ../ARCH/ANEMONE/XIOS/* arch/
./make_xios --full --prod --arch anemone-ifort-impi --netcdf_lib netcdf4_par --job 10

# compile nemo
cd $CDIR
echo "AMM7_SURGE OCE" >> ref_cfgs.txt
cd $CDIR/..
./makenemo -m anemone-ifort-impi -r AMM7_SURGE -j 16

# tidy up a messy exp directory
cd $EXP
rm xios_server.exe
rm nemo
rm bdydta
ln -s ../../../../xios/bin/xios_server.exe
ln -s ../EXP00/nemo

# grab inputs
cd ../
mkdir INPUTS
mkdir INPUTS/fluxes
mkdir INPUTS/bdydta

wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7_surge/domain_cfg.nc     -O INPUTS/amm7_surge_domain_cfg.nc
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7_surge/coordinates.bdy.nc -O INPUTS/coordinates.bdy.nc
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7_surge/bfr_coef.nc        -O INPUTS/bfr_coef.nc

cd INPUTS/fluxes
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/SBC/ERA5_U10_y2017.nc .
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/SBC/ERA5_V10_y2017.nc .
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/SBC/ERA5_MSL_y2017.nc .
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/SBC/ERA5_LSM.nc       .
wget http://gws-access.jasmin.ac.uk/public/jmmp/AMM7/inputs/SBC/weights_era5_bicubic.nc .

cd ../../
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

cd $EXP
ln -s ../INPUTS/
ln -s INPUTS/bdydta/
ln -s INPUTS/fluxes
ln -s ../INPUTS/amm7_surge_domain_cfg.nc
ln -s ../INPUTS/bfr_coef.nc
ln -s ../INPUTS/coordinates.bdy.nc

# submit a run
sbatch runscript_anemone.slurm
```
