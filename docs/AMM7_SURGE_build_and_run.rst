
.. _AMM7_SURGE_build_and_run-label:

************************************************
Setting up and running the NEMO AMM7 surge model
************************************************

This recipe describes how to build NEMO and XIOS appropriate for a surge model.

Example domain file (i.e. the grid and bathymetry), tidal and meterology forcing are shared here. Guidance for generating these files can be found elsewhere (see references).
To run NEMO you need to build NEMO and XIOS with the appropriate libraries. In my experience building the netcdf libraries is problematic.

Here we present a solution that runs on ARCHER2.ac.uk.
We also present a solution that uses containers so that the libraries can be more easily controlled. But you have to install Singularity, or equivalent. 

A basic "operating system for NEMO" container can be downloaded from  https://github.com/NOC-MSM/CoNES/releases/download/0.0.2/nemo_baseOS.sif, or better still you can follow instructions on the CoNES repository to build you own version. You will place this in your INPUTS directory. You can then open a command line shell in the container and it is like running on a new machine with all the libraries and compilers sorted out.

The source code in the directory `NEMO_4.0.4_surge` was made available from::

  https://code.metoffice.gov.uk/svn/nemo/NEMO/branches/UKMO/NEMO_4.0.4_surge/

In this processes the source code was stripped back to make a light weight configuration. In previous and future versions the surge configuration will be shared as edits to version controlled official releases. This surge configuration falls between the gaps as NEMO transitioned from svn to git and so a pragmatic solution is given.

Binaries files for the example forcings and domain file (best kept separate from git repositories) are available for download - links are on the README.
Download these and copy them into the `INPUTS` directory for the git cloned surge repository.

The remaining setup instructions are machine dependent. Two examples are given: 1) ARCHER2 HPC; 2) Using a Singularity container.

1) Building and running on ARCHER2 HPC
======================================

As an example the following process has been used on ARCHER2 using XIOS2 and the CRAY compiler (valid at 22nd Jan 2024)

Login to ARCHER2 ::

  ssh -l $USER login.archer2.ac.uk

Define some paths ::

  export CONFIG=AMM7_SURGE
  export WORK=/work/n01/n01
  export WDIR=/work/n01/n01/$USER/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export CDIR=$WDIR/NEMO_4.0.4_surge/cfgs # compile location is down here
  export EXP=$CDIR/$CONFIG/EXP_ERA5_DEMO

Clone the repository ::

  cd $WORK/$USER
  git clone -b feature/v4.0.4 https://github.com/JMMP-Group/AMM7_surge.git $CONFIG


Get compiler option files using a shared XIOS2 install::

  cd $CDIR/..
  cp /work/n01/shared/nemo/ARCH/*4.2.fcm arch/NOC/.

Load some modules::

  module swap craype-network-ofi craype-network-ucx
  module swap cray-mpich cray-mpich-ucx
  module load cray-hdf5-parallel/1.12.2.1
  module load cray-netcdf-hdf5parallel/4.9.0.1


Compile NEMO::

  cd $CDIR
  echo "AMM7_SURGE OCE" >> ref_cfgs.txt
  cd $CDIR/..
  ./makenemo -m X86_ARCHER2-Cray_4.2 -r AMM7_SURGE -j 16


Link executables to experiment directory (first remove any old existing placeholder links, as appropriate)::

  cd $EXP
  ln -s /work/n01/shared/nemo/XIOS2_Cray/bin/xios_server.exe $EXP/xios_server.exe
  ln -s $CDIR/$CONFIG/BLD/bin/nemo.exe $EXP/nemo

(N.B. sometimes the executable is expected to be called `opa` or `nemo.exe`)


Populate the INPUTS directory according to the REAME. Then make a link between binaries and where they are expected to be found (first remove any old existing placeholder links, as appropriate)::

    cd $EXP
    ln -s /work/n01/n01/shared/CO_AMM7/TIDE/FES               $EXP/bdydta
    ln -s /work/n01/n01/shared/CO_AMM15/INPUTS/forcing/era5   $EXP/fluxes
    ln -s $INPUTS/coordinates.bdy.nc $EXP/coordinates.bdy.nc
    ln -s $INPUTS/bfr_coef.nc        $EXP/bfr_coef.nc
    ln -s $INPUTS/amm7_surge_domain_cfg.nc $EXP/amm7_surge_domain_cfg.nc


Finally we are ready to submit a run script job from the experiment directory.
Edit the runscript.

Submit::

  cd $EXP
  sbatch runscript.slurm

Sea surface height is output every 15 mins::

  AMMSRG_met_15mi_20170101_20170206_grid_V.nc
  AMMSRG_met_15mi_20170101_20170206_grid_U.nc
  AMMSRG_met_15mi_20170101_20170206_grid_T.nc


2) Building and running in a Singularity Container
==================================================

There are many reasons why ARCHER2 is not suitable. Here is the workflow to get the code running in a Singularity container, on the assumption that you have already got Singularity installed.

Set up some paths::

  export CONFIG=AMM7_SURGE
  export WORK=/work/$USER/TEST
  export WDIR=$WORK/$CONFIG
  export GIT_DIR=$WORK/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export CDIR=$WDIR/NEMO_4.0.4_surge/cfgs # compile location is down here
  export XIOS_DIR=$WORK/XIOS2
  export EXP=$CDIR/$CONFIG/EXP_NOWIND_DEMO

This workflow includes the building of XIOS. The idea is to use a container with a controlled operating system and prebuilt libraries so that you can be confident that the NEMO and XIOS programs will compile::

  cd $WORK
  wget https://github.com/NOC-MSM/CoNES/releases/download/0.0.2/nemo_baseOS.sif  # 297Mb
  chmod u+x nemo_baseOS.sif
  singularity shell ./nemo_baseOS.sif



Set up some library paths that have been preprepared::

  PATH=$PATH:/opt/mpi/install/bin:/opt/hdf5/install/bin
  LD_LIBRARY_PATH=/opt/hdf5/install/lib:$LD_LIBRARY_PATH



Clone the configuration repository (and select the appropriate branch)::

  git clone -b feature/v4.0.4 https://github.com/JMMP-Group/AMM7_surge.git $CONFIG




Clone the XIOS repository, and copy in the arch files::

  cd $WORK
  svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS2/trunk XIOS2
  cd $XIOS_DIR
  cp $GIT_DIR/ARCH/SINGULARITY/xios/* arch/.


Compile::

  ./make_xios --full --debug --arch singularity --netcdf_lib netcdf4_par -j 8

NB ``./make_xios --full --prod --arch singularity --netcdf_lib netcdf4_par -j 8`` does not work...

This builds the ``$XIOS_DIR/bin/xios_server.exe`` executable and libraries, which need to be linked into the NEMO builds.

Edit the NEMO arch files to point to new XIOS builds::

  sed -i "s?XXX_XIOS_DIR_XXX?$XIOS_DIR?g" $GIT_DIR/ARCH/SINGULARITY/nemo/arch-singularity.fcm 



Copy arch files for NEMO build into place::
  
  cp $GIT_DIR/ARCH/SINGULARITY/nemo/*.fcm $CDIR/../arch/.



Compile NEMO, as before::

  cd $CDIR
  echo "AMM7_SURGE OCE" >> ref_cfgs.txt
  cd $CDIR/..
  ./makenemo -m singularity -r AMM7_SURGE -j 16


Link executables to experiment directory (first remove any old existing placeholder links, as appropriate)::

  ln -s $XIOS_DIR/bin/xios_server.exe $EXP/xios_server.exe
  ln -s $CDIR/$CONFIG/BLD/bin/nemo.exe $EXP/nemo

(N.B. sometimes the executable is expected to be called `opa` or `nemo.exe`)


Make a link between binaries and where they are expected to be found (first remove any old existing placeholder links, as appropriate)::

    ln -s $INPUTS/bdydta                   $EXP/bdydta
    ln -s $INPUTS/fluxes                   $EXP/fluxes   # Not needed for no-wind example
    ln -s $INPUTS/coordinates.bdy.nc       $EXP/coordinates.bdy.nc
    ln -s $INPUTS/bfr_coef.nc              $EXP/bfr_coef.nc
    ln -s $INPUTS/amm7_surge_domain_cfg.nc $EXP/amm7_surge_domain_cfg.nc


Run the configuration::

  mpirun -n 1 ./nemo : -n 1 ./xios_server.exe








************************************************
Generate tidal boundary conditions
************************************************

The tidal boundary conditions were generated from the FES2014 tidal model with a tool called ``PyBDY`` <https://github.com/NOC-MSM/pyBDY>
The boundary data are stored in ``$INPUTS``. Data are provided for this configuration. Notes for generating tidal files for other configurations can be found e.g. https://github.com/JMMP-Group/SEVERN-SWOT/wiki/04.-Make-tidal-boundary-conditions


************************************************
Generate surface forcing
************************************************

The surge model requires 10m wind velocity and atmospheric pressure. As a demonstration some example data is provided that has been processed from the ERA5 dataset. Data were processed using the tool ``pySBC`` <https://github.com/NOC-MSM/pySBC>

************************************************
2D bottom friction coefficients
************************************************

This surge configuration uses 2D spatially varying bottom friction coefficients. These were generated following Warder and Piggot (2022) 
and converted into depth vary drag coefficents. The map is implemented as a scaling on the existing drag coefficient.

************************************************
Useful references
************************************************

* SEAsia wiki notes (https://zenodo.org/record/6483231)
* Julian Mak's NEMO notes: https://nemo-related.readthedocs.io/en/latest/
* Collated guidance in Polton et al (2023). Reproducible and relocatable regional ocean modelling: fundamentals and practices. DOI: https://doi.org/10.5194/gmd-16-1481-2023
* Polton, J. A., Wise, A., O'Neill, C. K., & O'Dea, E. (2020). AMM7-surge: A 7km resolution Atlantic Margin Model surge configuration using NEMOv3.6 (v0.0.9). Zenodo. https://doi.org/10.5281/zenodo.4022310
* Old notes on generating boundary conditions from previous AMM7_surge release: https://github.com/JMMP-Group/AMM7_surge/blob/v0.0.9/docs/generate_tidal_boundaries.rst
* Warder and Piggott (2022). Optimal experiment design for a bottom friction parameter estimation problem. GEM - International Journal on Geomathematics  DOI: https://doi.org/10.1007/s13137-022-00196-4

