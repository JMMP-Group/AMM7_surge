MODULE obs_inter_h2d
   !!======================================================================
   !!                       ***  MODULE obs_inter_h2d   ***
   !! Observation diagnostics: Perform the horizontal interpolation
   !!                          from model grid to observation location
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   obs_int_h2d     : Horizontal interpolation to the observation point
   !!   obs_int_h2d_ds1 : Distance-weighted interpolation                 (n2dint=0)
   !!   obs_int_h2d_ds2 : Distance-weighted interpolation (small angle)   (n2dint=1)
   !!   obs_int_h2d_bil : Bilinear interpolation (geographical grid)      (n2dint=2)
   !!   obs_int_h2d_bir : Bilinear remapping interpolation (general grid) (n2dint=3)
   !!   obs_int_h2d_pol : Polynomial interpolation                        (n2dint=4)
   !!   bil_wgt         : Compute weights for bilinear remapping
   !!   lu_invmat       : Invert a matrix using LU decomposition
   !!   lu_decomp       : LU decomposition
   !!   lu_backsb       : LU decomposition - back substitution
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind, ONLY : &  ! Precision variables
      & wp
   USE phycst,   ONLY : &  ! Physical constants
      & rad,  &
      & rpi
   USE in_out_manager
   USE obs_const, ONLY : &
      & obfillflt		! Fillvalue
   USE obs_utils           ! Utility functions
      USE lib_mpp,ONLY : &
      & ctl_warn, ctl_stop

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE obs_int_h2d_ds1, & ! Distance-weighted interpolation               
      &    obs_int_h2d_ds2, & ! Distance-weighted interpolation (small angle) 
      &    obs_int_h2d_bil, & ! Bilinear interpolation (geographical grid)    
      &    obs_int_h2d_bir, & ! Bilinear remapping interpolation (general grid)
      &    obs_int_h2d_pol, & ! Polynomial interpolation                       
      &    lu_invmat,       & ! Invert a matrix using LU decomposition
      &    lu_decomp,       & ! LU decomposition
      &    lu_backsb,       & ! LU decomposition - back substitution
      &    bil_wgt            ! Compute weights for bilinear remapping
   PUBLIC obs_int_h2d,      & ! Horizontal interpolation to the observation point
      &   obs_int_h2d_init    ! Set up weights and vertical mask

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS
 
#include "obsinter_h2d.h90"

END MODULE obs_inter_h2d
