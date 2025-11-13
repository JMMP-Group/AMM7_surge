MODULE obs_profiles
   !!=====================================================================
   !!                       ***  MODULE  obs_profiles  ***
   !! Observation diagnostics: Storage space for profile observations
   !!                          arrays and additional flags etc.
   !!=====================================================================
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

   
   !! * Modules used 
   USE obs_profiles_def ! Definition of profile data types and tools

   IMPLICIT NONE
   
   SAVE

   !! * Routine accessibility
   PRIVATE

   PUBLIC nprofsets, nprofvars, nprofextr, profdata, prodatqc
   PUBLIC nvelosets, nvelovars, nveloextr, velodata, veldatqc
   
   !! * Shared Module variables
   INTEGER :: nprofsets                    ! Total number of profile data sets
   INTEGER :: nprofvars                    ! Total number of variables for profiles
   INTEGER :: nprofextr                    ! Extra fields for each variable
   TYPE(obs_prof), POINTER ::  profdata(:) ! Initial profile data
   TYPE(obs_prof), POINTER ::  prodatqc(:) ! Profile data after quality control

   INTEGER :: nvelosets                     ! Total number of velocity profile data sets
   INTEGER :: nvelovars                     ! Total number of variables for profiles
   INTEGER :: nveloextr                     ! Extra fields for each variable
   TYPE(obs_prof), POINTER ::  velodata(:)  ! Initial velocity profile data
   TYPE(obs_prof), POINTER ::  veldatqc(:)  ! Velocity profile data after quality control
END MODULE obs_profiles
