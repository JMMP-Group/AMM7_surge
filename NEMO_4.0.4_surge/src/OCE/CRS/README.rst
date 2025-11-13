**********************************
On line biogeochemistry coarsening
**********************************

.. todo::



.. contents::
   :local:

Presentation
============

A capacity of coarsening physics to force a BGC model coupled to NEMO has been developed.
This capacity allow to run 'online' a BGC model coupled to OCE-SI3 with a lower resolution,
to reduce the CPU cost of the BGC model, while preserving the effective resolution of the dynamics.

A presentation is available [attachment:crs_wiki_1.1.pdf​ here], where the methodology is presented.

What is available and working for now in this version
=====================================================

[To be completed]

Description of the successful validation tests
==============================================

[To be completed]

What is not working yet with on line coarsening of biogeochemistry
==================================================================

[To be completed]

''should include precise explanation on MPI decomposition problems too''

How to set up and use on line biogeochemistry
=============================================

Extract the on line biogeochemistry branch
------------------------------------------

To get the appropriate source code with the on line coarsening of biogeochemistry feature:

.. code-block:: console

	$ svn co https://forge.ipsl.jussieu.fr/nemo/browser/NEMO/branches/2018/dev_r5003_MERCATOR6_CRS


How to activate coarsening?
---------------------------

To activate the coarsening, ``key_crs`` should be added to list of CPP keys.
This key will only activate the coarsening of dynamics.

Some parameters are available in the namelist_cfg:

.. code-block:: fortran

	               !   passive tracer coarsened online simulations
	!-----------------------------------------------------------------------
	   nn_factx    = 3         !  Reduction factor of x-direction
	   nn_facty    = 3         !  Reduction factor of y-direction
	   nn_msh_crs  = 0         !  create (=1) a mesh file or not (=0)
	   nn_crs_kz   = 3         ! 0, volume-weighted MEAN of KZ
	                           ! 1, MAX of KZ
	                           ! 2, MIN of KZ
	                           ! 3, 10^(MEAN(LOG(KZ))
	                           ! 4, MEDIANE of KZ
	   ln_crs_wn   = .false.   ! wn coarsened (T) or computed using horizontal divergence ( F )
	                           !                           !
	   ln_crs_top = .true.     !coarsening online for the bio
	/

- Only ``nn_factx = 3`` is available and the coarsening only works for grids with a T-pivot point for
  the north-fold lateral boundary condition (ORCA025, ORCA12, ORCA36, ...).
- ``nn_msh_crs = 1`` will activate the generation of the coarsened grid meshmask.
- ``nn_crs_kz`` is the operator to coarsen the vertical mixing coefficient.
- ``ln_crs_wn``

  - when ``key_vvl`` is activated, this logical has no effect;
    the coarsened vertical velocities are computed using horizontal divergence.
  - when ``key_vvl`` is not activated,

    - coarsened vertical velocities are computed using horizontal divergence (``ln_crs_wn = .false.``)
    - or coarsened vertical velocities are computed with an average operator (``ln_crs_wn = .true.``)
- ``ln_crs_top = .true.``: should be activated to run BCG model in coarsened space;
  so only works when ``key_top`` is in the cpp list and eventually ``key_pisces`` or ``key_my_trc``.

Choice of operator to coarsene KZ
---------------------------------

A sensiblity test has been done with an Age tracer to compare the different operators.
The 3 and 4 options seems to provide the best results.

Some results can be found [xxx here]

Example of xml files to output coarsened variables with XIOS
------------------------------------------------------------

In the [attachment:iodef.xml iodef.xml]  file, a "nemo" context is defined and
some variable defined in [attachment:file_def.xml file_def.xml] are writted on the ocean-dynamic grid.
To write variables on the coarsened grid, and in particular the passive tracers,
a "nemo_crs" context should be defined in [attachment:iodef.xml iodef.xml] and
the associated variable are listed in [attachment:file_crs_def.xml file_crs_def.xml ].

Passive tracers tracers initial conditions
------------------------------------------

When initial conditions are provided in NetCDF files, the field might be:

- on the coarsened grid
- or they can be on another grid and
  interpolated `on-the-fly <http://forge.ipsl.jussieu.fr/nemo/wiki/Users/SetupNewConfiguration/Weight-creator>`_.
  Example of namelist for PISCES :

	.. code-block:: fortran

		!-----------------------------------------------------------------------
		&namtrc_dta      !    Initialisation from data input file
		!-----------------------------------------------------------------------
		!
		   sn_trcdta(1)  = 'DIC_REG1'        ,        -12        ,  'DIC'     ,    .false.   , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
		   sn_trcdta(2)  = 'ALK_REG1'        ,        -12        ,  'ALK'     ,    .false.   , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
		   sn_trcdta(3)  = 'O2_REG1'         ,        -1         ,  'O2'      ,    .true.    , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
		   sn_trcdta(5)  = 'PO4_REG1'        ,        -1         ,  'PO4'     ,    .true.    , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
		   sn_trcdta(7)  = 'Si_REG1'         ,        -1         ,  'Si'      ,    .true.    , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
	   	sn_trcdta(10) = 'DOC_REG1'        ,        -12        ,  'DOC'     ,    .false.   , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
		   sn_trcdta(14) = 'Fe_REG1'         ,        -12        ,  'Fe'      ,    .false.   , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
		   sn_trcdta(23) = 'NO3_REG1'        ,        -1         ,  'NO3'     ,    .true.    , .true. , 'yearly'  , 'reshape_REG1toeORCA075_bilin.nc'       , ''   , ''
	   	rn_trfac(1)   =   1.0e-06  !  multiplicative factor
		   rn_trfac(2)   =   1.0e-06  !  -      -      -     -
		   rn_trfac(3)   =  44.6e-06  !  -      -      -     -
	   	rn_trfac(5)   = 122.0e-06  !  -      -      -     -
		   rn_trfac(7)   =   1.0e-06  !  -      -      -     -
		   rn_trfac(10)  =   1.0e-06  !  -      -      -     -
	   	rn_trfac(14)  =   1.0e-06  !  -      -      -     -
		   rn_trfac(23)  =   7.6e-06  !  -      -      -     -

	   	cn_dir        =  './'      !  root directory for the location of the data files

PISCES forcing files
--------------------

They might be on the coarsened grid.

Perspectives
============

For the future, a few options are on the table to implement coarsening for biogeochemistry in 4.0 and
future releases.
Those will be discussed in Autumn 2018
