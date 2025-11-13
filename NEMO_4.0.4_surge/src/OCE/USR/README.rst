******************************
Setting up a new configuration
******************************

.. todo::



.. contents::
   :local:

Starting from an existing configuration
=======================================

There are three options to build a new configuration from an existing one.

Option 1: Duplicate an existing configuration
---------------------------------------------

The NEMO so-called Reference Configurations cover a number of major features for NEMO setup
(global, regional, 1D, using embedded zoom with AGRIF...)

One can create a new configuration by duplicating one of the reference configurations
(``ORCA2_ICE_PISCES`` in the following example)

.. code-block:: console

   $ ./makenemo –n 'ORCA2_ICE_PISCES_MINE' -r 'ORCA2_ICE_PISCES' -m 'my_arch'

Option 2: Duplicate with differences
------------------------------------

Create and compile a new configuration based on a reference configuration
(``ORCA2_ICE_PISCES`` in the following example) but with different pre-processor options.
For this either add ``add_key`` or ``del_key`` keys as required; e.g.

.. code-block:: console

   $ ./makenemo –n 'ORCA2_ICE_PISCES_MINE' -r 'ORCA2_ICE_PISCES' -m 'my_arch' del_key 'key_iomput' add_key 'key_diahth'

Option 3: Use the SIREN tools to subset an existing model
---------------------------------------------------------

Define a regional configuration which is a {sub,super}-set of an existing configuration.

This last option employs the SIREN software tools that are included in the standard distribution.
The software is written in Fortran 95 and available in the :file:`./tools/SIREN` directory.
SIREN allows you to create your own regional configuration embedded in a wider one.

SIREN is a set of programs to create all the input files you need to
run a NEMO regional configuration.

:Demo:     Set of GLORYS files (GLObal ReanalYSis on the ORCA025 grid),
           as well as examples of namelists are available `here`_.
:Doc:      :forge:`chrome/site/doc/SIREN/html/index.html`
:Support:  Any questions or comments regarding the use of SIREN should be posted in
  :forge:`the corresponding forum <discussion/forum/2>`.

.. _here: https://prodn.idris.fr/thredds/catalog/ipsl_public/rron463/catalog.html

Creating a completely new configuration
=======================================

From NEMO version 4.0 there are two ways to build configurations from scratch.
The appropriate method to use depends largely on the target configuration.
Method 1 is for more complex/realistic global or regional configurations and
method 2 is intended for simpler, idealised configurations whose
domains and characteristics can be described in simple geometries and formulae.

Option 1: Create and use a domain configuration file
----------------------------------------------------

This method is used by each of the reference configurations,
so that downloading their input files linked to their description can help.
Although starting from scratch,
it is advisable to create the directory structure to house your new configuration by
duplicating the closest reference configuration to your target application.
For example, if your application requires both ocean ice and passive tracers,
then use the ``ORCA2_ICE_PISCES`` as template,
and execute following command to build your ``MY_NEW_CONFIG`` configuration:

.. code-block:: sh

   $ ./makenemo –n 'MY_NEW_CONFIG' -r 'ORCA2_ICE_PISCES' -m 'my_arch'

where ``MY_NEW_CONFIG`` can be substituted with
a suitably descriptive name for your new configuration.

The purpose of this step is simply to create and populate the appropriate :file:`WORK`,
:file:`MY_SRC` and :file:`EXP00` subdirectories for your new configuration.
Other choices for the base reference configuration might be

:GYRE:  If your target application is ocean-only
:AMM12: If your target application is regional with open boundaries

All the domain information for your new configuration will be contained within
a netcdf file called :file:`domain_cfg.nc` which you will need to create and
place in the :file:`./cfgs/MY_NEW_CONFIG/EXP00` sub-directory.
Firstly though, ensure that your configuration is set to use such a file by checking that

.. code-block:: fortran

   ln_read_cfg = .true.

in :file:`./cfgs/MY_NEW_CONFIG/EXP00/namelist_cfg`

Create the :file:`domain_cfg.nc` file which must contain the following fields

.. code-block:: c

   /* configuration name, configuration resolution                 */
   int    ORCA, ORCA_index
   /* global domain sizes                                          */
   int    jpiglo, jpjglo, jpkglo
   /* lateral global domain b.c.                                   */
   int    jperio
   /* flags for z-coord, z-coord with partial steps and s-coord    */
   int    ln_zco, ln_zps, ln_sco
   /* flag  for ice shelf cavities                                 */
   int    ln_isfcav
   /* geographic position                                          */
   double glamt, glamu, glamv, glamf
   /* geographic position                                          */
   double gphit, gphiu, gphiv, gphif
   /* Coriolis parameter (if not on the sphere)                    */
   double iff, ff_f, ff_t
   /* horizontal scale factors                                     */
   double e1t, e1u, e1v, e1f
   /* horizontal scale factors                                     */
   double e2t, e2u, e2v, e2f
   /* U and V surfaces (if grid size reduction in some straits)    */
   double ie1e2u_v, e1e2u, e1e2v
   /* reference vertical scale factors at T and W points           */
   double e3t_1d, e3w_1d
   /* vertical scale factors 3D coordinate at T,U,V,F and W points */
   double e3t_0, e3u_0, e3v_0, e3f_0, e3w_0
   /* vertical scale factors 3D coordinate at UW and VW points     */
   double e3uw_0, e3vw_0
   /* last wet T-points, 1st wet T-points (for ice shelf cavities) */
   int    bottom_level, top_level

There are two options for creating a :file:`domain_cfg.nc` file:

- Users can use tools of their own choice to build a :file:`domain_cfg.nc` with all mandatory fields.
- Users can adapt and apply the supplied tool available in :file:`./tools/DOMAINcfg`.
  This tool is based on code extracted from NEMO version 3.6 and will allow similar choices for
  the horizontal and vertical grids that were available internally to that version.
  See :ref:`tools <DOMAINcfg>` for details.

Option 2: Adapt the usr_def configuration module of NEMO for you own purposes
-----------------------------------------------------------------------------

This method is intended for configuring easily simple/idealised configurations which
are often used as demonstrators or for process evaluation and comparison.
This method can be used whenever the domain geometry has a simple mathematical description and
the ocean initial state and boundary forcing is described analytically.
As a start, consider the case of starting a completely new ocean-only test case based on
the ``LOCK_EXCHANGE`` example.

.. note::

   We probably need an even more basic example than this with only one namelist and
   minimal changes to the usrdef modules

Firstly, construct the directory structure, starting in the :file:`cfgs` directory:

.. code-block:: console

   $ ./makenemo -n 'MY_NEW_TEST' -t 'LOCK_EXCHANGE' -m 'my_arch'

where the ``-t`` option has been used to locate the new configuration in
the :file:`tests` subdirectory
(it is recommended practice to keep full configurations and idealised cases clearly distinguishable).
This command will create (amongst others) the following files and directories::

   ./tests/MY_NEW_TEST:
   BLD  EXP00  MY_SRC WORK  cpp_MY_NEW_TEST.fcm

   ./tests/MY_NEW_TEST/EXP00:
   context_nemo.xml  domain_def_nemo.xml  field_def_nemo-oce.xml  file_def_nemo-oce.xml  iodef.xml
   namelist_cfg      namelist_ref

   ./tests/MY_NEW_TEST/MY_SRC:
   usrdef_hgr.F90  usrdef_nam.F90  usrdef_zgr.F90  usrdef_istate.F90  usrdef_sbc.F90  zdfini.F90

The key to setting up an idealised configuration lies in
adapting a small set of short Fortran 90 modules which
should be dropped into the :file:`MY_SRC` directory.
Here the ``LOCK_EXCHANGE`` example is using 5 such routines but the full set that is available in
the :file:`src/OCE/USR` directory is::

   ./src/OCE/USR:
   usrdef_closea.F90  usrdef_fmask.F90  usrdef_hgr.F90  usrdef_istate.F90
   usrdef_nam.F90     usrdef_sbc.F90    usrdef_zgr.F90

Before discussing these in more detail it is worth noting the various namelist controls that
engage the different user-defined aspects.
These controls are set using two new logical switches or are implied by the settings of existing ones.
For example, the mandatory requirement for an idealised configuration is to provide routines which
define the horizontal and vertical domains.
Templates for these are provided in the :file:`usrdef_hgr.F90` and :file:`usrdef_zgr.F90` modules.
The application of these modules is activated whenever:

.. code-block:: fortran

   ln_read_cfg = .false.

in any configuration's :file:`namelist_cfg` file.
This setting also activates the reading of an optional ``&nam_usrdef`` namelist which can be used to
supply configuration specific settings.
These need to be declared and read in the :file:`usrdef_nam.F90` module.

Another explicit control is available in the ``&namsbc`` namelist which
activates the use of analytical forcing.
With

.. code-block:: fortran

   ln_usr = .true.

Other usrdef modules are activated by less explicit means.
For example, code in :file:`usrdef_istate.F90` is used to
define initial temperature and salinity fields if

.. code-block:: fortran

   ln_tsd_init   = .false.

in the ``&namtsd`` namelist.
The remaining modules, namely :file:`usrdef_closea.F90` :file:`usrdef_fmask.F90` are specific to
ORCA configurations and set local variations of some specific fields for
the various resolutions of the global models.
They do not need to be considered here in the context of idealised cases but
it is worth noting that all configuration specific code has now been isolated in the usrdef modules.
In the case of these last two modules, they are activated only if an ORCA configuration is detected.
Currently,
this requires a specific integer variable named ``ORCA`` to be set in a :file:`domain_cfg.nc` file.

.. note::

   This would be less confusing if the ``cn_cfg`` string is read directly as
   a character attribue from the :file:`domain_cfg.nc`.

So, in most cases, the set up of idealised model configurations can be completed by
copying the template routines from :file:`./src/OCE/USR` into
your new :file:`./cfgs/MY_NEW_TEST/MY_SRC` directory and
editing the appropriate modules as needed.
The default set are those used for the GYRE reference configuration.
The contents of :file:`MY_SRC` directories from other idealised configurations may provide
more convenient templates if they share common characteristics with your target application.

Whatever the starting point,
it should not require too many changes or additional lines of code to produce routines in
:file:`./src/OCE/USR` that define analytically the domain,
the initial state and the surface boundary conditions for your new configuration.

To summarize, the base set of modules is:

:usrdef_hgr.F90:    Define horizontal grid
:usrdef_zgr.F90:    Define vertical grid
:usrdef_sbc.F90:    Provides at each time-step the surface boundary condition,
   i.e. the momentum, heat and freshwater fluxes
:usrdef_istate.F90: Defines initialization of the dynamics and tracers
:usrdef_nam.F90:    Configuration-specific namelist processing to
   set any associated run-time parameters

with two specialised ORCA modules
(not related to idealised configurations but used to isolate configuration specific code that
is used in ORCA2 reference configurations and established global configurations using
the ORCA tripolar grid):

:usrdef_fmask.F90:  only used in ORCA configurations for
   alteration of f-point land/ocean mask in some straits
:usrdef_closea.F90: only used in ORCA configurations for
   specific treatments associated with closed seas

From version 4.0, the NEMO release includes a :file:`tests` subdirectory containing available and
up to date :doc:`test cases <tests>` build by the community.
These will not be fully supported as are NEMO reference configurations,
but should provide a source of raw material.
