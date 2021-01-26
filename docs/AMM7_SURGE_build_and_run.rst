
.. _AMM7_SURGE_build_and_run-label:

************************************************
Setting up and running the NEMO AMM7 surge model
************************************************

This recipe describes how to build NEMO and XIOS appropriate for a surge model.
For simplicity, in this example meteorological forcing is switched off.
The method for generating the tidal boundary conditions is given. The method for
generating the domain file (i.e. the grid and bathymetry) is not given. However
the tidal boundary conditions and domain configuration file can be downloaded elsewhere.


1) Get NEMO codebase
====================

Login to ARCHER ::

  ssh -l $USER login.archer2.ac.uk

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

Add files to the experiment directory. This demonstration is tide-only::

  rsync -vt EXP_tideonly/* $CDIR/$CONFIG/EXP_tideonly


NB Have added a couple of extra lines into the field_def files. This is a glitch in the surge code,
because it doesn't expect you to not care about the winds::

  vi $EXP/field_def_nemo-opa.xml
  line 338
  <field id="wspd"         long_name="wind speed module"        standard_name="wind_speed"     unit="m/s" />                                                          unit="m/s"                            />
  <field id="uwnd"         long_name="u component of wind"       unit="m/s"         />
  <field id="vwnd"         long_name="v component of wind"       unit="m/s"        />


Load some modules::

  module unload craype-network-ofi
  module unload cray-mpich
  module load craype-network-ucx
  module load cray-mpich-ucx
  module load libfabric
  module load cray-hdf5-parallel
  module load cray-netcdf-hdf5parallel
  module load gcc

2) Build XIOS2.5 @ r2022
========================

Note when NEMO (nemo.exe / opa) is compiled it is done with reference to a particular version of
XIOS. So on NEMO run time the version of XIOS that built xios_server.exe must be compatible with the
version of XIOS that built nemo.exe / opa.


Download XIOS2.5 and prep::

  cd $WORK/$USER
  svn co -r2022 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5/  xios-2.5_r2022
  cd xios-2.5_r2022

Make a mod (line 894). Though you might need to run the ``make_xios`` command
once first to unpack the tar files::

  vi tools/FCM/lib/Fcm/Config.pm
  FC_MODSEARCH => '-J',              # FC flag, specify "module" path

Copy architecture files from git repo::

  cp $WORK/$USER/$CONFIG/ARCH/XIOS/arch-X86_ARCHER2-Cray* arch/.

Implement make command::

  ./make_xios --prod --arch X86_ARCHER2-Cray --netcdf_lib netcdf4_par --job 16 --full

Link the xios-2.5_r2022 to a generic XIOS directory name::

  ln -s  $WORK/$USER/xios-2.5_r2022  $WORK/$USER/XIOS

Link xios executable to the EXP directory::

  ln -s  $WORK/$USER/xios-2.5_r2022/bin/xios_server.exe $EXP/xios_server.exe



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

Put the HPC compiler file (from the git repo) in the correct place (this
currently uses xios2.5 from acc) ::

  rsync -vt $WDIR/ARCH/arch-X86_ARCHER2-Cray.fcm $CDIR/../ARCH/.


Make a mod (line 894). Though you might need to run the ``make_xios`` command
once first to unpack the tar files::

  vi $WDIR/dev_r8814_surge_modelling_Nemo4/EXTERNAL/fcm/lib/Fcm/Config.pm
  FC_MODSEARCH => '-J',              # FC flag, specify "module" path

Make NEMO ::

  cd $CDIR
  ./makenemo -n $CONFIG  -m X86_ARCHER2-Cray -j 16

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

The tidal boundary conditions were generated from the FES2014 tidal model with a tool called PyNEMO.
At this time the version of PyNEMO did not support outputting only 2D tidal forcing,
so some of the error checking for 3D boundary conditions is not needed but has
to be satisfied. This is how it was done. A new version of PyNEMO now exists.
The boundary data are stored in ``$INPUTS``

See **generate tidal boundaries** page.

6) Summary of external requirements
===================================

To successfully run NEMO will expect a ``coordinates.bdy.nc`` file in `$INPUTS`
(generated by PyNEMO) it will also expect boundary files of the type::

  AMM7_surge_bdytide_rotT_*.nc
  amm7_bdytide_*.nc

E.g. ``AMM7_surge_bdytide_rotT_M2_grid_V.nc`` and ``amm7_bdytide_M2_grid_T.nc``

There must also be a ``domain_cfg.nc`` domain file in ``$EXP``.


7) Run NEMO
===========

Finally we are ready to submit a run script job from the experiment directory.

Make the runscript (to be downloaded from repo but not settled on processor
split yet). For example, to run with 4 xios servers (a maximum of 2 per node),
each with sole occupancy of a 16-core NUMA region and 96 ocean cores, spaced
with an idle core in between each, use::

  cd $EXP
  /work/n01/shared/acc/mkslurm -S 4 -s 16 -m 2 -C 96 -c 2 > runscript.slurm

(rename executable in script from "nemo" to "opa")

Submit::

  cd $EXP
  mkdir Restart_files
  sbatch runscript.slurm

Sea surface height is output every 15 mins.
