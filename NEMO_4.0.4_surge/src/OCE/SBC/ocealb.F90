MODULE ocealb
   !!======================================================================
   !!                       ***  MODULE  ocealb  ***
   !! Ocean forcing:  bulk thermohaline forcing of the ocean
   !!=====================================================================
   !! History :
   !!   NEMO     4.0  ! 2017-07  (C. Rousset) Split ocean and ice albedos
   !!----------------------------------------------------------------------
   !!   oce_alb    : albedo for ocean (clear and overcast skies)
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   oce_alb   ! routine called by sbccpl
  
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE oce_alb( palb_os , palb_cs )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE oce_alb  ***
      !! 
      !! ** Purpose :   Computation of the albedo of the ocean
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   palb_os   !  albedo of ocean under overcast sky
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   palb_cs   !  albedo of ocean under clear sky
      !!
      REAL(wp) ::   zcoef 
      REAL(wp) ::   rmue = 0.40    !  cosine of local solar altitude
      !!----------------------------------------------------------------------
      !
      zcoef = 0.05 / ( 1.1 * rmue**1.4 + 0.15 )   ! Parameterization of Briegled and Ramanathan, 1982
      palb_cs(:,:) = zcoef 
      palb_os(:,:) = 0.06                         ! Parameterization of Kondratyev, 1969 and Payne, 1972
      !
   END SUBROUTINE oce_alb

   !!======================================================================
END MODULE ocealb
