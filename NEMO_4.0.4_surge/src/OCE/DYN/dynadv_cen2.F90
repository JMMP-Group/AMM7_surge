MODULE dynadv_cen2
   !!======================================================================
   !!                       ***  MODULE  dynadv  ***
   !! Ocean dynamics: Update the momentum trend with the flux form advection
   !!                 using a 2nd order centred scheme
   !!======================================================================
   !! History :  2.0  ! 2006-08  (G. Madec, S. Theetten)  Original code
   !!            3.2  ! 2009-07  (R. Benshila)  Suppression of rigid-lid option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_adv_cen2  : flux form momentum advection (ln_dynadv_cen2=T) using a 2nd order centred scheme  
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_adv_cen2   ! routine called by step.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_adv_cen2( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_cen2  ***
      !!
      !! ** Purpose :   Compute the now momentum advection trend in flux form
      !!              and the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time) 
      !!
      !! ** Action  :   (ua,va) updated with the now vorticity term trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zfu_t, zfu_f, zfu_uw, zfu
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zfv_t, zfv_f, zfv_vw, zfv, zfw
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_adv_cen2 : 2nd order flux form momentum advection'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      IF( l_trddyn ) THEN           ! trends: store the input trends
         zfu_uw(:,:,:) = ua(:,:,:)
         zfv_vw(:,:,:) = va(:,:,:)
      ENDIF
      !
      !                             !==  Horizontal advection  ==!
      !
      DO jk = 1, jpkm1                    ! horizontal transport
         zfu(:,:,jk) = 0.25_wp * e2u(:,:) * e3u_n(:,:,jk) * un(:,:,jk)
         zfv(:,:,jk) = 0.25_wp * e1v(:,:) * e3v_n(:,:,jk) * vn(:,:,jk)
         DO jj = 1, jpjm1                 ! horizontal momentum fluxes (at T- and F-point)
            DO ji = 1, fs_jpim1   ! vector opt.
               zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj,jk) ) * ( un(ji,jj,jk) + un(ji+1,jj  ,jk) )
               zfv_f(ji  ,jj  ,jk) = ( zfv(ji,jj,jk) + zfv(ji+1,jj,jk) ) * ( un(ji,jj,jk) + un(ji  ,jj+1,jk) )
               zfu_f(ji  ,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji,jj+1,jk) ) * ( vn(ji,jj,jk) + vn(ji+1,jj  ,jk) )
               zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji,jj+1,jk) ) * ( vn(ji,jj,jk) + vn(ji  ,jj+1,jk) )
            END DO
         END DO
         DO jj = 2, jpjm1                 ! divergence of horizontal momentum fluxes
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - (  zfu_t(ji+1,jj,jk) - zfu_t(ji,jj  ,jk)    &
                  &                           + zfv_f(ji  ,jj,jk) - zfv_f(ji,jj-1,jk)  ) * r1_e1e2u(ji,jj) / e3u_n(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) - (  zfu_f(ji,jj  ,jk) - zfu_f(ji-1,jj,jk)    &
                  &                           + zfv_t(ji,jj+1,jk) - zfv_t(ji  ,jj,jk)  ) * r1_e1e2v(ji,jj) / e3v_n(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      IF( l_trddyn ) THEN           ! trends: send trend to trddyn for diagnostic
         zfu_uw(:,:,:) = ua(:,:,:) - zfu_uw(:,:,:)
         zfv_vw(:,:,:) = va(:,:,:) - zfv_vw(:,:,:)
         CALL trd_dyn( zfu_uw, zfv_vw, jpdyn_keg, kt )
         zfu_t(:,:,:) = ua(:,:,:)
         zfv_t(:,:,:) = va(:,:,:)
      ENDIF
      !
      !                             !==  Vertical advection  ==!
      !
      DO jj = 2, jpjm1                    ! surface/bottom advective fluxes set to zero
         DO ji = fs_2, fs_jpim1
            zfu_uw(ji,jj,jpk) = 0._wp   ;   zfv_vw(ji,jj,jpk) = 0._wp
            zfu_uw(ji,jj, 1 ) = 0._wp   ;   zfv_vw(ji,jj, 1 ) = 0._wp
         END DO
      END DO
      IF( ln_linssh ) THEN                ! linear free surface: advection through the surface
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zfu_uw(ji,jj,1) = 0.5_wp * ( e1e2t(ji,jj) * wn(ji,jj,1) + e1e2t(ji+1,jj) * wn(ji+1,jj,1) ) * un(ji,jj,1)
               zfv_vw(ji,jj,1) = 0.5_wp * ( e1e2t(ji,jj) * wn(ji,jj,1) + e1e2t(ji,jj+1) * wn(ji,jj+1,1) ) * vn(ji,jj,1)
            END DO
         END DO
      ENDIF
      DO jk = 2, jpkm1                    ! interior advective fluxes
         DO jj = 2, jpj                       ! 1/4 * Vertical transport
            DO ji = 2, jpi
               zfw(ji,jj,jk) = 0.25_wp * e1e2t(ji,jj) * wn(ji,jj,jk)
            END DO
         END DO
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zfu_uw(ji,jj,jk) = ( zfw(ji,jj,jk) + zfw(ji+1,jj  ,jk) ) * ( un(ji,jj,jk) + un(ji,jj,jk-1) )
               zfv_vw(ji,jj,jk) = ( zfw(ji,jj,jk) + zfw(ji  ,jj+1,jk) ) * ( vn(ji,jj,jk) + vn(ji,jj,jk-1) )
            END DO
         END DO
      END DO
      DO jk = 1, jpkm1                    ! divergence of vertical momentum flux divergence
         DO jj = 2, jpjm1 
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - ( zfu_uw(ji,jj,jk) - zfu_uw(ji,jj,jk+1) ) * r1_e1e2u(ji,jj) / e3u_n(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) - ( zfv_vw(ji,jj,jk) - zfv_vw(ji,jj,jk+1) ) * r1_e1e2v(ji,jj) / e3v_n(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      IF( l_trddyn ) THEN                 ! trends: send trend to trddyn for diagnostic
         zfu_t(:,:,:) = ua(:,:,:) - zfu_t(:,:,:)
         zfv_t(:,:,:) = va(:,:,:) - zfv_t(:,:,:)
         CALL trd_dyn( zfu_t, zfv_t, jpdyn_zad, kt )
      ENDIF
      !                                   ! Control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' cen2 adv - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=           ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
   END SUBROUTINE dyn_adv_cen2

   !!==============================================================================
END MODULE dynadv_cen2
