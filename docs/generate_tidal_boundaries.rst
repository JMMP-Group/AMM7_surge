
.. _generate_tidal_boundaries-label:

Generate tidal boundary conditions
==================================

The boundary conditions were generated with a tool called PyNEMO. At this time
the version of PyNEMO did not support outputting only 2D tidal forcing,
so some of the error checking for 3D boundary conditions is not needed but has
to be satisfied. This is how it was done. A new version of PyNEMO now exists.



Install PyNEMO
==============

Method for installing on a centos box. NB At this time PyNEMO was briefly called ``nrct``::

  export CONFIG=AMM7_SURGE
  export WORK=/work
  export WDIR=$WORK/$USER/NEMO/$CONFIG
  export INPUTS=$WDIR/INPUTS
  export START_FILES=$WDIR/START_FILES

  cd $WORK/$USER
  mkdir $WDIR
  module load anaconda/2.1.0  # Want python2
  conda create --name nrct_env scipy=0.16.0 numpy=1.9.2 matplotlib=1.4.3 basemap=1.0.7 netcdf4=1.1.9 libgfortran=1.0.0
  source activate nrct_env
  conda install -c https://conda.anaconda.org/conda-forge seawater=3.3.4 # Note had to add https path
  conda install -c https://conda.anaconda.org/srikanthnagella thredds_crawler
  conda install -c https://conda.anaconda.org/srikanthnagella pyjnius

Find java object by doing a ``which java`` and then following the trail
``find  /usr/lib/jvm/jre-1.7.0-openjdk.x86_64/ -name libjvm.so -print``
::

  cd $WORK/$USER
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH
  unset SSH_ASKPASS # Didn't need this on ARCHER...
  git clone https://jpolton@bitbucket.org/jdha/nrct.git nrct  # Give jpolton@bitbucket passwd
  cd $WORK/$USER/nrct/Python

For generation of FES tides need to checkout the appropriate branch ::

    git fetch
    git checkout Generalise-tide-input

You have to manually set the TPXO or FES data source in ``Python/pynemo/tide/nemo_bdy_tide3.py``
Check this. Then build PyNEMO ::

  python setup.py build
  export PYTHONPATH=/login/$USER/.conda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH
  python setup.py install --prefix ~/.conda/envs/nrct_env
  cd $INPUTS

If everything else is ready (it isn't) you could then just do: ``pynemo -s namelist.bdy``


Create inputs files for PyNEMO
==============================

PyNEMO is controlled by a ``namelist.bdy`` file and which
has two input files: ``inputs_src.ncml`` which points to the data source and
``inputs_dst.ncml`` which remaps some variable names in the destination files.

Copy across some parent mesh files and a mask file (even though they are not
used. This is because this old version of PyNEMO didn't anticipate tide-only usage)::

  cp ../../SEAsia/INPUTS/mesh_?gr_src.nc $INPUTS/.
  cp ../../SEAsia/INPUTS/mask_src.nc $INPUTS/.
  cp ../../SEAsia/INPUTS/inputs_dst.ncml $INPUTS/.
  cp ../../SEAsia/INPUTS/cut_inputs_src.ncml $INPUTS/.


If I don't make a boundary mask then it doesn't work... This can also be done with
the PyNEMO GUI. The mask variable takes values (-1 mask, 1 wet, 0 land). Get a
template from ``domain_cfg.nc`` and then modify as desired around the boundary ::

  module load nco/gcc/4.4.2.ncwa
  rm -f bdy_mask.nc tmp[12].nc
  ncks -v top_level domain_cfg.nc tmp1.nc
  ncrename -h -v top_level,mask tmp1.nc tmp2.nc
  ncwa -a t tmp2.nc bdy_mask.nc
  rm -f tmp[12].nc

Then in ipython::

  import netCDF4, numpy
  dset = netCDF4.Dataset('bdy_mask.nc','a')
  dset.variables['mask'][0,:]  = -1     # Southern boundary
  dset.variables['mask'][-1,:] = -1    # Northern boundary
  dset.variables['mask'][:,-1] = -1    # Eastern boundary
  dset.variables['mask'][:,0] = -1        # Western boundary
  dset.close()



Make a bathymetry file from envelope bathymetry variable ``hbatt`` (I think
this is OK to do..)::

  module load nco/gcc/4.4.2.ncwa
  rm -f hbatt.nc tmp1.nc tmp2.nc
  ncks -v hbatt, nav_lat, nav_lon domain_cfg.nc tmp1.nc
  ncrename -h -v hbatt,Bathymetry tmp1.nc tmp2.nc
  ncwa -a t tmp2.nc hbatt.nc


FES2014 tidal data is used as the tidal data source. This is clumsily set in
``nemo_bdy_tide3.py`` before pynemo is built, though the attached
``namelist.bdy`` has redundant references to TPXO. The FES2014 files will need
to be obtained from https://datastore.cls.fr/catalogues/fes2014-tide-model/


Run PyNEMO
==========

Generate the boundary conditions with PyNEMO
::

  module load anaconda/2.1.0  # Want python2
  source activate nrct_env
  cd $INPUTS
  export LD_LIBRARY_PATH=/usr/lib/jvm/jre-1.7.0-openjdk.x86_64/lib/amd64/server:$LD_LIBRARY_PATH
  export PYTHONPATH=/login/$USER/.conda/envs/nrct_env/lib/python2.7/site-packages/:$PYTHONPATH

  pynemo -s namelist.bdy


This creates::

  coordinates.bdy.nc
  AMM7_surge_bdytide_rotT_NU2_grid_T.nc
  AMM7_surge_bdytide_rotT_O1_grid_T.nc
  AMM7_surge_bdytide_rotT_P1_grid_T.nc
  AMM7_surge_bdytide_rotT_Q1_grid_T.nc
  AMM7_surge_bdytide_rotT_MTM_grid_T.nc
  AMM7_surge_bdytide_rotT_MU2_grid_T.nc
  AMM7_surge_bdytide_rotT_N2_grid_T.nc
  AMM7_surge_bdytide_rotT_N4_grid_T.nc
  AMM7_surge_bdytide_rotT_R2_grid_T.nc
  AMM7_surge_bdytide_rotT_S1_grid_T.nc
  AMM7_surge_bdytide_rotT_2N2_grid_T.nc
  AMM7_surge_bdytide_rotT_J1_grid_T.nc
  AMM7_surge_bdytide_rotT_EPS2_grid_T.nc
  AMM7_surge_bdytide_rotT_K2_grid_T.nc
  AMM7_surge_bdytide_rotT_K1_grid_T.nc
  AMM7_surge_bdytide_rotT_LA2_grid_T.nc
  AMM7_surge_bdytide_rotT_L2_grid_T.nc
  AMM7_surge_bdytide_rotT_M3_grid_T.nc
  AMM7_surge_bdytide_rotT_M2_grid_T.nc
  AMM7_surge_bdytide_rotT_M6_grid_T.nc
  AMM7_surge_bdytide_rotT_M4_grid_T.nc
  AMM7_surge_bdytide_rotT_MF_grid_T.nc
  AMM7_surge_bdytide_rotT_M8_grid_T.nc
  AMM7_surge_bdytide_rotT_MM_grid_T.nc
  AMM7_surge_bdytide_rotT_MKS2_grid_T.nc
  AMM7_surge_bdytide_rotT_MS4_grid_T.nc
  AMM7_surge_bdytide_rotT_MN4_grid_T.nc
  AMM7_surge_bdytide_rotT_MSQM_grid_T.nc
  AMM7_surge_bdytide_rotT_MSF_grid_T.nc
  AMM7_surge_bdytide_rotT_S4_grid_T.nc
  AMM7_surge_bdytide_rotT_S2_grid_T.nc
  AMM7_surge_bdytide_rotT_T2_grid_T.nc
  AMM7_surge_bdytide_rotT_SSA_grid_T.nc
  AMM7_surge_bdytide_rotT_SA_grid_T.nc
  AMM7_surge_bdytide_rotT_NU2_grid_U.nc
  AMM7_surge_bdytide_rotT_O1_grid_U.nc
  AMM7_surge_bdytide_rotT_P1_grid_U.nc
  AMM7_surge_bdytide_rotT_Q1_grid_U.nc
  AMM7_surge_bdytide_rotT_MTM_grid_U.nc
  AMM7_surge_bdytide_rotT_MU2_grid_U.nc
  AMM7_surge_bdytide_rotT_N2_grid_U.nc
  AMM7_surge_bdytide_rotT_N4_grid_U.nc
  AMM7_surge_bdytide_rotT_R2_grid_U.nc
  AMM7_surge_bdytide_rotT_S1_grid_U.nc
  AMM7_surge_bdytide_rotT_2N2_grid_U.nc
  AMM7_surge_bdytide_rotT_J1_grid_U.nc
  AMM7_surge_bdytide_rotT_EPS2_grid_U.nc
  AMM7_surge_bdytide_rotT_K2_grid_U.nc
  AMM7_surge_bdytide_rotT_K1_grid_U.nc
  AMM7_surge_bdytide_rotT_LA2_grid_U.nc
  AMM7_surge_bdytide_rotT_L2_grid_U.nc
  AMM7_surge_bdytide_rotT_M3_grid_U.nc
  AMM7_surge_bdytide_rotT_M2_grid_U.nc
  AMM7_surge_bdytide_rotT_M6_grid_U.nc
  AMM7_surge_bdytide_rotT_M4_grid_U.nc
  AMM7_surge_bdytide_rotT_MF_grid_U.nc
  AMM7_surge_bdytide_rotT_M8_grid_U.nc
  AMM7_surge_bdytide_rotT_MM_grid_U.nc
  AMM7_surge_bdytide_rotT_MKS2_grid_U.nc
  AMM7_surge_bdytide_rotT_MS4_grid_U.nc
  AMM7_surge_bdytide_rotT_MN4_grid_U.nc
  AMM7_surge_bdytide_rotT_MSQM_grid_U.nc
  AMM7_surge_bdytide_rotT_MSF_grid_U.nc
  AMM7_surge_bdytide_rotT_S4_grid_U.nc
  AMM7_surge_bdytide_rotT_S2_grid_U.nc
  AMM7_surge_bdytide_rotT_T2_grid_U.nc
  AMM7_surge_bdytide_rotT_SSA_grid_U.nc
  AMM7_surge_bdytide_rotT_SA_grid_U.nc
  AMM7_surge_bdytide_rotT_NU2_grid_V.nc
  AMM7_surge_bdytide_rotT_O1_grid_V.nc
  AMM7_surge_bdytide_rotT_P1_grid_V.nc
  AMM7_surge_bdytide_rotT_Q1_grid_V.nc
  AMM7_surge_bdytide_rotT_MTM_grid_V.nc
  AMM7_surge_bdytide_rotT_MU2_grid_V.nc
  AMM7_surge_bdytide_rotT_N2_grid_V.nc
  AMM7_surge_bdytide_rotT_N4_grid_V.nc
  AMM7_surge_bdytide_rotT_R2_grid_V.nc
  AMM7_surge_bdytide_rotT_S1_grid_V.nc
  AMM7_surge_bdytide_rotT_2N2_grid_V.nc
  AMM7_surge_bdytide_rotT_J1_grid_V.nc
  AMM7_surge_bdytide_rotT_EPS2_grid_V.nc
  AMM7_surge_bdytide_rotT_K2_grid_V.nc
  AMM7_surge_bdytide_rotT_K1_grid_V.nc
  AMM7_surge_bdytide_rotT_LA2_grid_V.nc
  AMM7_surge_bdytide_rotT_L2_grid_V.nc
  AMM7_surge_bdytide_rotT_M3_grid_V.nc
  AMM7_surge_bdytide_rotT_M2_grid_V.nc
  AMM7_surge_bdytide_rotT_M6_grid_V.nc
  AMM7_surge_bdytide_rotT_M4_grid_V.nc
  AMM7_surge_bdytide_rotT_MF_grid_V.nc
  AMM7_surge_bdytide_rotT_M8_grid_V.nc
  AMM7_surge_bdytide_rotT_MM_grid_V.nc
  AMM7_surge_bdytide_rotT_MKS2_grid_V.nc
  AMM7_surge_bdytide_rotT_MS4_grid_V.nc
  AMM7_surge_bdytide_rotT_MN4_grid_V.nc
  AMM7_surge_bdytide_rotT_MSQM_grid_V.nc
  AMM7_surge_bdytide_rotT_MSF_grid_V.nc
  AMM7_surge_bdytide_rotT_S4_grid_V.nc
  AMM7_surge_bdytide_rotT_S2_grid_V.nc
  AMM7_surge_bdytide_rotT_T2_grid_V.nc
  AMM7_surge_bdytide_rotT_SSA_grid_V.nc
  AMM7_surge_bdytide_rotT_SA_grid_V.nc


Copy these files into ``$INPUTS`` on ARCHER.
