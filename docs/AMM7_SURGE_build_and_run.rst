
.. _AMM7_SURGE_build_and_run:

*****************************************
Setting up the NEMO AMM7 surge model
*****************************************

This recipe principly describes how to build NEMO and XIOS approprite for a
 surge model. The method for generating the tidal boundary conditions is
  outlined and exhaustive details on how to generate the grid and bathymetry
  are omitted since these require data from external sources.

1) Get NEMO codebase
====================

Login to ARCHER ::

  ssh -l $USER login.archer.ac.uk

Define some paths ::

  export CONFIG=AMM7_SURGE
  export WORK=/work/n01/n01
  export WDIR=/work/n01/n01/$USER/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES
  export CDIR=$WDIR/dev_r8814_surge_modelling_Nemo4/CONFIG
  export TDIR=$WDIR/dev_r8814_surge_modelling_Nemo4/TOOLS
  export EXP=$CDIR/$CONFIG/EXP_tideonly


Clone the repository ::

  cd $WORK/$USER
  git clone https://github.com/JMMP-Group/AMM7_surge.git $CONFIG

Get the code::

  cd $WDIR
  svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/UKMO/dev_r8814_surge_modelling_Nemo4/NEMOGCM dev_r8814_surge_modelling_Nemo4

Make a link between where the inputs files are and where the model expects them ::

    ln -s $INPUTS $EXP/bdydta

Put files from git repo into ``MY_SRC``::

  rsync -vt MY_SRC/* $CDIR/$CONFIG/MY_SRC

Add files to the experiment directory. Demonstrate with tide-only::

  rsync -vt EXP_tideonly/* $CDIR/$CONFIG/EXP_tideonly


NB Have added a couple of extra lines into the field_def files. This is a glitch in the surge code,
because it doesn't expect you to not care about the winds::

  vi $EXP/field_def_nemo-opa.xml
  line 338
  <field id="wspd"         long_name="wind speed module"        standard_name="wind_speed"     unit="m/s" />                                                          unit="m/s"                            />
  <field id="uwnd"         long_name="u component of wind"       unit="m/s"         />
  <field id="vwnd"         long_name="v component of wind"       unit="m/s"        />


Load modules ::

  module load cray-netcdf-hdf5parallel/4.4.1.1
  module load cray-hdf5-parallel/1.10.0.1
  module swap PrgEnv-cray PrgEnv-intel/5.2.82

2) Build XIOS2 @ r1242
======================

Note when NEMO (nemo.exe / opa) is compiled it is done with reference to a particular version of
XIOS. So on NEMO run time the version of XIOS that built xios_server.exe must be compatible with the
version of XIOS that built nemo.exe / opa.

Modules::

  module load cray-netcdf-hdf5parallel/4.4.1.1
  module load cray-hdf5-parallel/1.10.0.1
  module swap PrgEnv-cray PrgEnv-intel/5.2.82

Download XIOS2 and prep::

  cd $WORK/$USER
  svn co -r1242 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk xios-2.0_r1242
  cd xios-2.0_r1242
  cp $WORK/$USER/$CONFIG/ARCH/arch-XC30_ARCHER* arch/.

Implement make command::

  ./make_xios --full --prod --arch XC30_ARCHER --netcdf_lib netcdf4_par --job 8

Link the xios-2.0_r1242 to a generic XIOS directory name::

  ln -s  $WORK/$USER/xios-2.0_r1242  $WORK/$USER/XIOS

Link xios executable to the EXP directory::

  ln -s  $WORK/$USER/xios-2.0_r1242/bin/xios_server.exe $EXP/xios_server.exe



3) Build NEMO
==============

Already got NEMO branch ::

    #cd $WDIR
    #svn co http://forge.ipsl.jussieu.fr/nemo/svn/branches/UKMO/dev_r8814_surge_modelling_Nemo4/NEMOGCM dev_r8814_surge_modelling_Nemo4



Copy files required to build ``nemo.exe``. Or get it from git repo. Or get it here.
Set the compile flags (will use the FES tide) ::

  vi $CDIR/$CONFIG/cpp_AMM7_SURGE.fcm
  bld::tool::fppkeys  key_nosignedzero key_diainstant key_mpp_mpi key_iomput  \
                      key_diaharm_fast key_FES14_tides

Put the HPC compiler file (from the git repo) in the correct place ::

  rsync -vt $WDIR/ARCH/arch-XC_ARCHER_INTEL.fcm $CDIR/../ARCH/.


Make NEMO ::

  cd $CDIR
  ./makenemo -n $CONFIG -m XC_ARCHER_INTEL -j 10


Copy executable to experiment directory ::

  ln -s $CDIR/$CONFIG/BLD/bin/nemo.exe $EXP/opa



4) Generate a domain configuration file
========================================

Copy a domain file that holds all the coordinates and domain discretisation.
This files is called ``domain_cfg.nc``. The generation of this file is not
described here. Obtain the file E.g. ::

  cd /projects/jcomp/fred/SURGE/AMM7_INPUTS
  scp amm7_surge_domain_cfg.nc jelt@login.archer.ac.uk:$INPUTS/domain_cfg.nc
  ln -s $INPUTS/domain_cfg.nc $EXP/.


5) Generate tidal boundary conditions
======================================

The tidal boundary conditions were generated from the FES2014 tidal model with
 a tool called PyNEMO. At this time the version of PyNEMO did not support outputting only 2D tidal forcing,
so some of the error checking for 3D boundary conditions is not needed but has
to be satisfied. This is how it was done. A new version of PyNEMO now exists.
The boundary data are stored in ``$INPUTS``

See :ref:`generate_tidal_boundaries`.


6) Run NEMO
===========

Submit a run script job from the experiment directory ::

  cd $EXP
  mkdir Restart_files
  qsub runscript



/work/n01/n01/jelt/AMM7_SURGE/dev_r8814_surge_modelling_Nemo4/CONFIG/AMM7_SURGE/EXP_tideonly
