MODULE crslbclnk
   !!======================================================================
   !!                       ***  MODULE  crslbclnk  ***
   !!               A temporary solution for lbclnk for coarsened grid.
   !! Ocean        : lateral boundary conditions for grid coarsening
   !!=====================================================================
   !! History :   ! 2012-06  (J. Simeon, G. Madec, C. Ethe, C. Calone)     Original code
   !!----------------------------------------------------------------------
   USE par_kind, ONLY: wp
   USE dom_oce
   USE crs
   !
   USE lbclnk
   USE in_out_manager
   
   INTERFACE crs_lbc_lnk
      MODULE PROCEDURE crs_lbc_lnk_3d, crs_lbc_lnk_2d
   END INTERFACE
   
   PUBLIC crs_lbc_lnk
   
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE crs_lbc_lnk_3d( pt3d1, cd_type1, psgn, kfillmode, pfillval  )
      !!---------------------------------------------------------------------
      !!                  ***  SUBROUTINE crs_lbc_lnk  ***
      !!
      !! ** Purpose :   set lateral boundary conditions for coarsened grid
      !!
      !! ** Method  :   Swap domain indices from full to coarse domain
      !!                before arguments are passed directly to lbc_lnk.
      !!                Upon exiting, switch back to full domain indices.
      !!----------------------------------------------------------------------
      CHARACTER(len=1)                        , INTENT(in   ) ::   cd_type1 ! grid type
      REAL(wp)                                , INTENT(in   ) ::   psgn     ! control of the sign
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(inout) ::   pt3d1    ! 3D array on which the lbc is applied
      INTEGER                     , OPTIONAL  , INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = cst)
      REAL(wp)                    , OPTIONAL  , INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      !
      LOGICAL  ::   ll_grid_crs
      !!----------------------------------------------------------------------
      !
      ll_grid_crs = ( jpi == jpi_crs )
      !
      IF( .NOT.ll_grid_crs )   CALL dom_grid_crs   ! Save the parent grid information  & Switch to coarse grid domain
      !
      CALL lbc_lnk( 'crslbclnk', pt3d1, cd_type1, psgn, kfillmode, pfillval )
      !
      IF( .NOT.ll_grid_crs )   CALL dom_grid_glo   ! Return to parent grid domain
      !
   END SUBROUTINE crs_lbc_lnk_3d
   
   
   SUBROUTINE crs_lbc_lnk_2d(pt2d, cd_type, psgn, kfillmode, pfillval )
      !!---------------------------------------------------------------------
      !!                  ***  SUBROUTINE crs_lbc_lnk  ***
      !!
      !! ** Purpose :   set lateral boundary conditions for coarsened grid
      !!
      !! ** Method  :   Swap domain indices from full to coarse domain
      !!                before arguments are passed directly to lbc_lnk.
      !!                Upon exiting, switch back to full domain indices.
      !!----------------------------------------------------------------------
      CHARACTER(len=1)                    , INTENT(in   ) ::   cd_type  ! grid type
      REAL(wp)                            , INTENT(in   ) ::   psgn     ! control of the sign
      REAL(wp), DIMENSION(jpi_crs,jpj_crs), INTENT(inout) ::   pt2d     ! 2D array on which the lbc is applied
      INTEGER                 , OPTIONAL  , INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = constant)
      REAL(wp)                , OPTIONAL  , INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      !      
      LOGICAL  ::   ll_grid_crs
      !!----------------------------------------------------------------------
      !
      ll_grid_crs = ( jpi == jpi_crs )
      !
      IF( .NOT.ll_grid_crs )   CALL dom_grid_crs   ! Save the parent grid information  & Switch to coarse grid domain
      !
      CALL lbc_lnk( 'crslbclnk', pt2d, cd_type, psgn, kfillmode, pfillval )
      !
      IF( .NOT.ll_grid_crs )   CALL dom_grid_glo   ! Return to parent grid domain
      !
   END SUBROUTINE crs_lbc_lnk_2d

   !!======================================================================
END MODULE crslbclnk
