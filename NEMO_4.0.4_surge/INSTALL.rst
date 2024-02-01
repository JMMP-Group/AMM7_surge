*******************
Build the framework
*******************

.. todo::



.. contents::
   :local:

Prerequisites
=============

| The NEMO source code is written in *Fortran 95* and
  some of its prerequisite tools and libraries are already included in the download.
| It contains the AGRIF_ preprocessing program ``conv``; the FCM_ build system and
  the IOIPSL_ library for parts of the output.

System dependencies
-------------------

In the first place the other requirements should be provided natively by your system or
can be installed from the official repositories of your Unix-like distribution:

- *Perl* interpreter
- *Fortran* compiler (``ifort``, ``gfortran``, ``pgfortran``, ...),
- *Message Passing Interface (MPI)* implementation (e.g. |OpenMPI|_ or |MPICH|_).
- |NetCDF|_ library with its underlying |HDF|_

**NEMO, by default, takes advantage of some MPI features introduced into the MPI-3 standard.**

.. hint::

   The MPI implementation is not strictly essential
   since it is possible to compile and run NEMO on a single processor.
   However most realistic configurations will require the parallel capabilities of NEMO and
   these use the MPI standard.

.. note::

   On older systems, that do not support MPI-3 features,
   the ``key_mpi2`` preprocessor key should be used at compile time.
   This will limit MPI features to those defined within the MPI-2 standard
   (but will lose some performance benefits).

.. |OpenMPI| replace:: *OpenMPI*
.. _OpenMPI: https://www.open-mpi.org
.. |MPICH|   replace:: *MPICH*
.. _MPICH:   https://www.mpich.org
.. |NetCDF|  replace:: *Network Common Data Form (NetCDF)*
.. _NetCDF:  https://www.unidata.ucar.edu
.. |HDF|     replace:: *Hierarchical Data Form (HDF)*
.. _HDF:     https://www.hdfgroup.org

Specifics for NetCDF and HDF
----------------------------

NetCDF and HDF versions from official repositories may have not been compiled with MPI support.
However access to all the options available with the XIOS IO-server will require
the parallelism of these libraries.

| **To satisfy these requirements, it is common to have to compile from source
  in this order HDF (C library) then NetCDF (C and Fortran libraries)**
| It is also necessary to compile these libraries with the same version of the MPI implementation that
  both NEMO and XIOS (see below) have been compiled and linked with.

.. hint::

   | It is difficult to define the options for the compilation as
     they differ from one architecture to another according to
     the hardware used and the software installed.
   | The following is provided without any warranty

   .. code-block:: console

      $ ./configure [--{enable-fortran,disable-shared,enable-parallel}] ...

   It is recommended to build the tests ``--enable-parallel-tests`` and run them with ``make check``

Particular versions of these libraries may have their own restrictions.
State the following requirements for netCDF-4 support:

.. caution::

   | When building NetCDF-C library versions older than 4.4.1, use only HDF5 1.8.x versions.
   | Combining older NetCDF-C versions with newer HDF5 1.10 versions will create superblock 3 files
     that are not readable by lots of older software.

Extract and install XIOS
========================

With the sole exception of running NEMO in mono-processor mode
(in which case output options are limited to those supported by the ``IOIPSL`` library),
diagnostic outputs from NEMO are handled by the third party ``XIOS`` library.
It can be used in two different modes:

:*attached*:  Every NEMO process also acts as a XIOS server
:*detached*:  Every NEMO process runs as a XIOS client.
  Output is collected and collated by external, stand-alone XIOS server processors.

Instructions on how to install XIOS can be found on its :xios:`wiki<>`.

.. hint::

   It is recommended to use XIOS 2.5 release.
   This version should be more stable (in terms of future code changes) than the XIOS trunk.
   It is also the one used by the NEMO system team when testing all developments and new releases.

   This particular version has its own branch and can be checked out with:

   .. code:: console

      $ svn co https://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5

Download and install the NEMO code
==================================

Checkout the NEMO sources
-------------------------

.. code:: console

   $ svn co https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-4.0-HEAD

Description of 1\ :sup:`st` level tree structure
------------------------------------------------

+---------------+----------------------------------------+
| :file:`arch`  | Compilation settings                   |
+---------------+----------------------------------------+
| :file:`cfgs`  | :doc:`Reference configurations <cfgs>` |
+---------------+----------------------------------------+
| :file:`doc`   | :doc:`Documentation <doc>`             |
+---------------+----------------------------------------+
| :file:`ext`   | Dependencies included                  |
|               | (``AGRIF``, ``FCM`` & ``IOIPSL``)      |
+---------------+----------------------------------------+
| :file:`mk`    | Compilation scripts                    |
+---------------+----------------------------------------+
| :file:`src`   | :doc:`Modelling routines <src>`        |
+---------------+----------------------------------------+
| :file:`tests` | :doc:`Test cases <tests>`              |
|               | (unsupported)                          |
+---------------+----------------------------------------+
| :file:`tools` | :doc:`Utilities <tools>`               |
|               | to {pre,post}process data              |
+---------------+----------------------------------------+

Setup your architecture configuration file
------------------------------------------

All compiler options in NEMO are controlled using files in :file:`./arch/arch-'my_arch'.fcm` where
``my_arch`` is the name of the computing architecture
(generally following the pattern ``HPCC-compiler`` or ``OS-compiler``).
It is recommended to copy and rename an configuration file from an architecture similar to your owns.
You will need to set appropriate values for all of the variables in the file.
In particular the FCM variables:
``%NCDF_HOME``; ``%HDF5_HOME`` and ``%XIOS_HOME`` should be set to
the installation directories used for XIOS installation

.. code-block:: sh

   %NCDF_HOME    /usr/local/path/to/netcdf
   %HDF5_HOME    /usr/local/path/to/hdf5
   %XIOS_HOME    /home/$( whoami )/path/to/xios-2.5
   %OASIS_HOME   /home/$( whoami )/path/to/oasis

Create and compile a new configuration
======================================

The main script to {re}compile and create executable is called :file:`makenemo` located at
the root of the working copy.
It is used to identify the routines you need from the source code, to build the makefile and run it.
As an example, compile a :file:`MY_GYRE` configuration from GYRE with 'my_arch':

.. code-block:: sh

   ./makenemo –m 'my_arch' –r GYRE -n 'MY_GYRE'

Then at the end of the configuration compilation,
:file:`MY_GYRE` directory will have the following structure.

+------------+----------------------------------------------------------------------------+
| Directory  | Purpose                                                                    |
+============+============================================================================+
| ``BLD``    | BuiLD folder: target executable, headers, libs, preprocessed routines, ... |
+------------+----------------------------------------------------------------------------+
| ``EXP00``  | Run   folder: link to executable, namelists, ``*.xml`` and IOs             |
+------------+----------------------------------------------------------------------------+
| ``EXPREF`` | Files under version control only for :doc:`official configurations <cfgs>` |
+------------+----------------------------------------------------------------------------+
| ``MY_SRC`` | New routines or modified copies of NEMO sources                            |
+------------+----------------------------------------------------------------------------+
| ``WORK``   | Links to all raw routines from :file:`./src` considered                    |
+------------+----------------------------------------------------------------------------+

After successful execution of :file:`makenemo` command,
the executable called `nemo` is available in the :file:`EXP00` directory

More :file:`makenemo` options
-----------------------------

``makenemo`` has several other options that can control which source files are selected and
the operation of the build process itself.

.. literalinclude:: ../../../makenemo
   :language: text
   :lines: 119-143
   :caption: Output of ``makenemo -h``

These options can be useful for maintaining several code versions with only minor differences but
they should be used sparingly.
Note however the ``-j`` option which should be used more routinely to speed up the build process.
For example:

.. code-block:: sh

        ./makenemo –m 'my_arch' –r GYRE -n 'MY_GYRE' -j 8

will compile up to 8 processes simultaneously.

Default behaviour
-----------------

At the first use,
you need the ``-m`` option to specify the architecture configuration file
(compiler and its options, routines and libraries to include),
then for next compilation, it is assumed you will be using the same compiler.
If the ``-n`` option is not specified the last compiled configuration will be used.

Tools used during the process
-----------------------------

* :file:`functions.sh`: bash functions used by ``makenemo``, for instance to create the WORK directory
* :file:`cfg.txt`     : text list of configurations and source directories
* :file:`bld.cfg`     : FCM rules for compilation

Examples
--------

.. literalinclude:: ../../../makenemo
   :language: text
   :lines: 146-153

Running the model
=================

Once :file:`makenemo` has run successfully,
the ``nemo`` executable is available in :file:`./cfgs/MY_CONFIG/EXP00`.
For the reference configurations, the :file:`EXP00` folder also contains the initial input files
(namelists, ``*.xml`` files for the IOs, ...).
If the configuration needs other input files, they have to be placed here.

.. code-block:: sh

   cd 'MY_CONFIG'/EXP00
   mpirun -n $NPROCS ./nemo   # $NPROCS is the number of processes
                              # mpirun is your MPI wrapper

Viewing and changing list of active CPP keys
============================================

For a given configuration (here called ``MY_CONFIG``),
the list of active CPP keys can be found in :file:`./cfgs/'MYCONFIG'/cpp_MY_CONFIG.fcm`

This text file can be edited by hand or with :file:`makenemo` to change the list of active CPP keys.
Once changed, one needs to recompile ``nemo`` in order for this change to be taken in account.
Note that most NEMO configurations will need to specify the following CPP keys:
``key_iomput`` for IOs and ``key_mpp_mpi`` for parallelism.
