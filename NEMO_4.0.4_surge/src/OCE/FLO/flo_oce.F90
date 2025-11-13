MODULE flo_oce
   !!======================================================================
   !!                     ***  MODULE flo_oce  ***
   !! lagrangian floats :   define in memory all floats parameters and variables
   !!======================================================================
   !! History :   OPA  ! 1999-10  (CLIPPER projet)
   !!   NEMO      1.0  ! 2002-11  (G. Madec, A. Bozec)  F90: Free form and module
   !!----------------------------------------------------------------------
   USE par_oce         ! ocean parameters
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PUBLIC

   PUBLIC   flo_oce_alloc   ! Routine called in floats.F90

   !! float parameters
   !! ----------------
   LOGICAL, PUBLIC ::   ln_floats   !: Activate floats or not
   INTEGER, PUBLIC ::   jpnfl       !: total number of floats during the run
   INTEGER, PUBLIC ::   jpnnewflo   !: number of floats added in a new run
   INTEGER, PUBLIC ::   jpnrstflo   !: number of floats for the restart

   !! float variables
   !! ---------------
   INTEGER , PUBLIC, ALLOCATABLE, DIMENSION(:) ::   nisobfl   !: =0 for a isobar float , =1 for a float following the w velocity
   INTEGER , PUBLIC, ALLOCATABLE, DIMENSION(:) ::   ngrpfl    !: number to identify searcher group
   INTEGER , PUBLIC, ALLOCATABLE, DIMENSION(:) ::   nfloat    !: number to identify searcher group

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   flxx , flyy , flzz    !: long, lat, depth of float (decimal degree, m >0)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   tpifl, tpjfl, tpkfl   !: (i,j,k) indices of float position

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wb   !: vertical velocity at previous time step (m s-1).
   
   !                                   !! * namelist namflo : langrangian floats *
   LOGICAL, PUBLIC  ::   ln_rstflo      !: T/F float restart 
   LOGICAL, PUBLIC  ::   ln_argo        !: T/F argo type floats
   LOGICAL, PUBLIC  ::   ln_flork4      !: T/F 4th order Runge-Kutta
   LOGICAL, PUBLIC  ::   ln_ariane      !: handle ariane input/output convention
   LOGICAL, PUBLIC  ::   ln_flo_ascii   !: write in ascii (T) or in Netcdf (F)

   INTEGER, PUBLIC  ::   nn_writefl     !: frequency of float output file 
   INTEGER, PUBLIC  ::   nn_stockfl     !: frequency of float restart file

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION flo_oce_alloc()
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION flo_oce_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( wb(jpi,jpj,jpk) , nfloat(jpnfl) , nisobfl(jpnfl) , ngrpfl(jpnfl) , &
         &      flxx(jpnfl)     , flyy(jpnfl)   , flzz(jpnfl)    ,                 & 
         &      tpifl(jpnfl)    , tpjfl(jpnfl)  , tpkfl(jpnfl)   , STAT=flo_oce_alloc )
      !
      CALL mpp_sum ( 'flo_oce', flo_oce_alloc )
      IF( flo_oce_alloc /= 0 )   CALL ctl_stop( 'STOP', 'flo_oce_alloc: failed to allocate arrays' )
   END FUNCTION flo_oce_alloc

   !!======================================================================
END MODULE flo_oce
