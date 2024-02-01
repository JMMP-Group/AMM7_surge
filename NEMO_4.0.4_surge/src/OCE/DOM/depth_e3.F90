MODULE depth_e3
   !!======================================================================
   !!                       ***  MODULE  depth_e3  ***
   !!
   !! zgr : vertical coordinate system 
   !!======================================================================
   !! History :  4.0  ! 2016-11  (S. Flavoni, G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   depth_to_e3   : use the depth of t- and w-points to calculate e3t & e3w
   !!                   (generic interface for 1D and 3D fields)
   !!   e3_to_depth   : use e3t & e3w to calculate the depth of t- and w-points
   !!                   (generic interface for 1D and 3D fields)
   !!---------------------------------------------------------------------
   USE oce               ! ocean variables
   USE dom_oce           ! ocean domain
   !
   USE in_out_manager    ! I/O manager
   USE lbclnk            ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp           ! distributed memory computing library
   USE timing            ! Timing

   IMPLICIT NONE
   PRIVATE
  
   INTERFACE depth_to_e3
      MODULE PROCEDURE depth_to_e3_1d, depth_to_e3_3d
   END INTERFACE

   INTERFACE e3_to_depth
      MODULE PROCEDURE e3_to_depth_1d, e3_to_depth_3d
   END INTERFACE

   PUBLIC   depth_to_e3        ! called by usrdef_zgr
   PUBLIC   e3_to_depth        ! called by domzgr.F90
      
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS             

   SUBROUTINE depth_to_e3_1d( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE depth_to_e3_1d  ***
      !!
      !! ** Purpose :   compute e3t & e3w scale factors from t- & w-depths of model levels
      !!
      !! ** Method  :   The scale factors are given by the discrete derivative 
      !!              of the depth:
      !!                               e3w(jk) = dk[ dept_1d ] 
      !!                               e3t(jk) = dk[ depw_1d ]
      !!              with, at top and bottom :
      !!                      e3w( 1 ) = 2 * ( dept( 1 ) - depw( 1 ) )
      !!                      e3t(jpk) = 2 * ( dept(jpk) - depw(jpk) )   
      !!
      !! ** Action  : - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(in   ) ::   pdept_1d, pdepw_1d   ! depths          [m]
      REAL(wp), DIMENSION(:), INTENT(  out) ::   pe3t_1d , pe3w_1d    ! e3.=dk[depth]   [m]
      !
      INTEGER  ::   jk           ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      ! use pdep. at w- and t-points to compute e3. (e3. = dk[depth])
      !
      pe3w_1d( 1 ) = 2._wp * ( pdept_1d(1) - pdepw_1d(1) ) 
      DO jk = 1, jpkm1
         pe3w_1d(jk+1) = pdept_1d(jk+1) - pdept_1d(jk) 
         pe3t_1d(jk  ) = pdepw_1d(jk+1) - pdepw_1d(jk) 
      END DO
      pe3t_1d(jpk) = 2._wp * ( pdept_1d(jpk) - pdepw_1d(jpk) )
      !
   END SUBROUTINE depth_to_e3_1d
   
      
   SUBROUTINE depth_to_e3_3d( pdept_3d, pdepw_3d, pe3t_3d, pe3w_3d )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE depth_to_e3_3d  ***
      !!
      !! ** Purpose :   compute e3t & e3w scale factors from t- & w-depths of model levels
      !!
      !! ** Method  :   The scale factors are given by the discrete derivative 
      !!              of the depth:
      !!                               e3w(jk) = dk[ dept_1d ] 
      !!                               e3t(jk) = dk[ depw_1d ]
      !!              with, at top and bottom :
      !!                      e3w( 1 ) = 2 * ( dept( 1 ) - depw( 1 ) )
      !!                      e3t(jpk) = 2 * ( dept(jpk) - depw(jpk) )   
      !!
      !! ** Action  : - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdept_3d, pdepw_3d   ! depth           [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t_3d , pe3w_3d    ! e3.=dk[depth]   [m]
      !
      INTEGER  ::   jk           ! dummy loop indices
      !!----------------------------------------------------------------------      
      pe3w_3d(:,:, 1 ) = 2._wp * ( pdept_3d(:,:,1) - pdepw_3d(:,:,1) ) 
      DO jk = 1, jpkm1
         pe3w_3d(:,:,jk+1) = pdept_3d(:,:,jk+1) - pdept_3d(:,:,jk) 
         pe3t_3d(:,:,jk  ) = pdepw_3d(:,:,jk+1) - pdepw_3d(:,:,jk) 
      END DO
      pe3t_3d(:,:,jpk) = 2._wp * ( pdept_3d(:,:,jpk) - pdepw_3d(:,:,jpk) )   
      !
   END SUBROUTINE depth_to_e3_3d


   SUBROUTINE e3_to_depth_1d( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE e3_to_depth_1d  ***
      !!
      !! ** Purpose :   compute t- & w-depths of model levels from e3t & e3w scale factors
      !!
      !! ** Method  :   The t- & w-depth are given by the summation of e3w & e3t, resp. 
      !!
      !! ** Action  : - pe3t_1d, pe3w_1d : scale factor of t- and w-point (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(in   ) ::   pe3t_1d , pe3w_1d    ! vert. scale factors   [m]
      REAL(wp), DIMENSION(:), INTENT(  out) ::   pdept_1d, pdepw_1d   ! depth = SUM( e3 )     [m]
      !
      INTEGER  ::   jk           ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      pdepw_1d(1) = 0.0_wp
      pdept_1d(1) = 0.5_wp * pe3w_1d(1)
      DO jk = 2, jpk
         pdepw_1d(jk) = pdepw_1d(jk-1) + pe3t_1d(jk-1) 
         pdept_1d(jk) = pdept_1d(jk-1) + pe3w_1d(jk  ) 
      END DO
      !
   END SUBROUTINE e3_to_depth_1d
   
      
   SUBROUTINE e3_to_depth_3d( pe3t_3d, pe3w_3d, pdept_3d, pdepw_3d )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE e3_to_depth_3d  ***
      !!
      !! ** Purpose :   compute t- & w-depths of model levels from e3t & e3w scale factors
      !!
      !! ** Method  :   The t- & w-depth are given by the summation of e3w & e3t, resp. 
      !!
      !! ** Action  : - pe3t_1d, pe3w_1d : scale factor of t- and w-point (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pe3t_3d , pe3w_3d    ! vert. scale factors   [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept_3d, pdepw_3d   ! depth = SUM( e3 )     [m]
      !
      INTEGER  ::   jk           ! dummy loop indices
      !!----------------------------------------------------------------------      
      !
      pdepw_3d(:,:,1) = 0.0_wp
      pdept_3d(:,:,1) = 0.5_wp * pe3w_3d(:,:,1)
      DO jk = 2, jpk
         pdepw_3d(:,:,jk) = pdepw_3d(:,:,jk-1) + pe3t_3d(:,:,jk-1) 
         pdept_3d(:,:,jk) = pdept_3d(:,:,jk-1) + pe3w_3d(:,:,jk  ) 
      END DO
      !
   END SUBROUTINE e3_to_depth_3d

   !!======================================================================
END MODULE depth_e3
