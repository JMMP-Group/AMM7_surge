MODULE diaprod
! Requires key_iom_put
# if defined key_iomput
   !!======================================================================
   !!                     ***  MODULE  diaprod  ***
   !! Ocean diagnostics :  write ocean product diagnostics
   !!=====================================================================
   !! History :  3.4  ! 2012  (D. Storkey)  Original code
   !!            4.0  ! 2019  (D. Storkey)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_prod      : calculate and write out product diagnostics
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE domvvl          ! for thickness weighted diagnostics if key_vvl
   USE eosbn2          ! equation of state  (eos call)
   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager
   USE iom
   USE ioipsl
   USE lib_mpp         ! MPP library
   USE timing          ! preformance summary

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_prod                 ! routines called by step.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.4 , NEMO Consortium (2012)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_prod( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_prod  ***
      !!                   
      !! ** Purpose :   Write out product diagnostics (uT, vS etc.)
      !!
      !! ** Method  :  use iom_put
      !!               Product diagnostics are not thickness-weighted in this routine.
      !!               They should be thickness-weighted using XIOS if key_vvl is set. 
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER                      ::   ji, jj, jk              ! dummy loop indices
      REAL(wp)                     ::   zztmp, zztmpx, zztmpy   ! 
      !!
      REAL(wp), POINTER, DIMENSION(:,:)   :: z2d      ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) :: z3d      ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zrhop    ! potential density
      !!----------------------------------------------------------------------
      ! 
      IF( ln_timing )   CALL timing_start('dia_prod')
      ! 
      ALLOCATE( z2d(jpi,jpj), z3d(jpi,jpj,jpk), zrhop(jpi,jpj,jpk) )
      !

      IF( iom_use("urhop") .OR. iom_use("vrhop") .OR. iom_use("wrhop") &
#if ! defined key_diaar5
     &  .OR. iom_use("rhop") &
#endif
     & ) THEN
         CALL eos( tsn, z3d, zrhop )                 ! now in situ and potential density
         zrhop(:,:,:) = zrhop(:,:,:)-1000.e0         ! reference potential density to 1000 to avoid precision issues in rhop2 calculation
         zrhop(:,:,jpk) = 0._wp
#if ! defined key_diaar5
         CALL iom_put( 'rhop', zrhop )
#else
         ! If key_diaar5 set then there is already an iom_put call to output rhop.
         ! Really should be a standard diagnostics option?
#endif
      ENDIF

      IF( iom_use("ut") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z3d(ji,jj,jk) = un(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_tem) + tsn(ji+1,jj,jk,jp_tem) )
               END DO
            END DO
         END DO
         CALL iom_put( "ut", z3d )                  ! product of temperature and zonal velocity at U points
      ENDIF

      IF( iom_use("vt") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z3d(ji,jj,jk) = vn(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_tem) + tsn(ji,jj+1,jk,jp_tem) )
               END DO
            END DO
         END DO
         CALL iom_put( "vt", z3d )                  ! product of temperature and meridional velocity at V points
      ENDIF

      IF( iom_use("wt") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               z3d(ji,jj,1) = wn(ji,jj,1) * tsn(ji,jj,1,jp_tem)
            END DO
         END DO
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z3d(ji,jj,jk) = wn(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk-1,jp_tem) + tsn(ji,jj,jk,jp_tem) )
               END DO
            END DO
         END DO
         CALL iom_put( "wt", z3d )                  ! product of temperature and vertical velocity at W points
      ENDIF

      IF( iom_use("us") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z3d(ji,jj,jk) = un(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_sal) + tsn(ji+1,jj,jk,jp_sal) )
               END DO
            END DO
         END DO
         CALL iom_put( "us", z3d )                  ! product of salinity and zonal velocity at U points
      ENDIF

      IF( iom_use("vs") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z3d(ji,jj,jk) = vn(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_sal) + tsn(ji,jj+1,jk,jp_sal) )
               END DO
            END DO
         END DO
         CALL iom_put( "vs", z3d )                  ! product of salinity and meridional velocity at V points
      ENDIF

      IF( iom_use("ws") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               z3d(ji,jj,1) = wn(ji,jj,1) * tsn(ji,jj,1,jp_sal)
            END DO
         END DO
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z3d(ji,jj,jk) = wn(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk-1,jp_sal) + tsn(ji,jj,jk,jp_sal) )
               END DO
            END DO
         END DO
         CALL iom_put( "ws", z3d )                  ! product of salinity and vertical velocity at W points
      ENDIF

      IF( iom_use("uv") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z3d(ji,jj,jk) = 0.25 * ( un(ji-1,jj,jk) + un(ji,jj,jk) ) * ( vn(ji,jj-1,jk) + vn(ji,jj,jk) ) 
               END DO
            END DO
         END DO
         CALL iom_put( "uv", z3d )                  ! product of zonal velocity and meridional velocity at T points
      ENDIF

      IF( iom_use("uw") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               z3d(ji,jj,1) = 0.5 * ( wn(ji,jj,1) + wn(ji+1,jj,1) ) * un(ji,jj,1) 
            END DO
         END DO
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z3d(ji,jj,jk) = 0.25 * ( wn(ji,jj,jk) + wn(ji+1,jj,jk) ) * ( un(ji,jj,jk-1) + un(ji,jj,jk) ) 
               END DO
            END DO
         END DO
         CALL iom_put( "uw", z3d )                  ! product of zonal velocity and vertical velocity at UW points
      ENDIF

      IF( iom_use("vw") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               z3d(ji,jj,1) = 0.5 * ( wn(ji,jj,1) + wn(ji,jj+1,1) ) * vn(ji,jj,1) 
            END DO
         END DO
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z3d(ji,jj,jk) = 0.25 * ( wn(ji,jj,jk) + wn(ji,jj+1,jk) ) * ( vn(ji,jj,jk-1) + vn(ji,jj,jk) ) 
               END DO
            END DO
         END DO
         CALL iom_put( "vw", z3d )                  ! product of meriodional velocity and vertical velocity at VW points
      ENDIF

      IF( iom_use("urhop") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z3d(ji,jj,jk) = un(ji,jj,jk) * 0.5 * ( zrhop(ji,jj,jk) + zrhop(ji+1,jj,jk) )
               END DO
            END DO
         END DO
         CALL iom_put( "urhop", z3d )                  ! product of density and zonal velocity at U points
      ENDIF

      IF( iom_use("vrhop") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z3d(ji,jj,jk) = vn(ji,jj,jk) * 0.5 * ( zrhop(ji,jj,jk) + zrhop(ji,jj+1,jk) )
               END DO
            END DO
         END DO
         CALL iom_put( "vrhop", z3d )                  ! product of density and meridional velocity at V points
      ENDIF

      IF( iom_use("wrhop") ) THEN
         z3d(:,:,:) = 0.e0 
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               z3d(ji,jj,1) = wn(ji,jj,1) * zrhop(ji,jj,1)
            END DO
         END DO
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z3d(ji,jj,jk) = wn(ji,jj,jk) * 0.5 * ( zrhop(ji,jj,jk-1) + zrhop(ji,jj,jk) )
               END DO
            END DO
         END DO
         CALL iom_put( "wrhop", z3d )                  ! product of density and vertical velocity at W points
      ENDIF

      !
      DEALLOCATE( z2d, z3d, zrhop )
      !
      IF( ln_timing )   CALL timing_stop('dia_prod')
      !
   END SUBROUTINE dia_prod
#else
   !!----------------------------------------------------------------------
   !!   Default option :                                         NO diaprod
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER :: lk_diaprod = .FALSE.   ! coupled flag
CONTAINS
   SUBROUTINE dia_prod( kt )   ! Empty routine
      INTEGER ::   kt
      WRITE(*,*) 'dia_prod: You should not have seen this print! error?', kt
   END SUBROUTINE dia_prod
#endif
   !!======================================================================
END MODULE diaprod
