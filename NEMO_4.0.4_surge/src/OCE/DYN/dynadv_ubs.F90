MODULE dynadv_ubs
   !!======================================================================
   !!                       ***  MODULE  dynadv_ubs  ***
   !! Ocean dynamics: Update the momentum trend with the flux form advection
   !!                 trend using a 3rd order upstream biased scheme
   !!======================================================================
   !! History :  2.0  ! 2006-08  (R. Benshila, L. Debreu)  Original code
   !!            3.2  ! 2009-07  (R. Benshila)  Suppression of rigid-lid option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_adv_ubs   : flux form momentum advection using    (ln_dynadv=T)
   !!                   an 3rd order Upstream Biased Scheme or Quick scheme
   !!                   combined with 2nd or 4th order finite differences 
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   REAL(wp), PARAMETER :: gamma1 = 1._wp/3._wp  ! =1/4 quick      ; =1/3  3rd order UBS
   REAL(wp), PARAMETER :: gamma2 = 1._wp/32._wp ! =0   2nd order  ; =1/32 4th order centred

   PUBLIC   dyn_adv_ubs   ! routine called by step.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_adv_ubs( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_ubs  ***
      !!
      !! ** Purpose :   Compute the now momentum advection trend in flux form
      !!              and the general trend of the momentum equation.
      !!
      !! ** Method  :   The scheme is the one implemeted in ROMS. It depends 
      !!      on two parameter gamma1 and gamma2. The former control the 
      !!      upstream baised part of the scheme and the later the centred 
      !!      part:     gamma1 = 0    pure centered  (no diffusive part)
      !!                       = 1/4  Quick scheme
      !!                       = 1/3  3rd order Upstream biased scheme
      !!                gamma2 = 0    2nd order finite differencing 
      !!                       = 1/32 4th order finite differencing
      !!      For stability reasons, the first term of the fluxes which cor-
      !!      responds to a second order centered scheme is evaluated using  
      !!      the now velocity (centered in time) while the second term which  
      !!      is the diffusive part of the scheme, is evaluated using the 
      !!      before velocity (forward in time). 
      !!      Default value (hard coded in the begining of the module) are 
      !!      gamma1=1/3 and gamma2=1/32.
      !!
      !! ** Action : - (ua,va) updated with the 3D advective momentum trends
      !!
      !! Reference : Shchepetkin & McWilliams, 2005, Ocean Modelling. 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zui, zvj, zfuj, zfvi, zl_u, zl_v   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk)   ::   zfu_t, zfu_f, zfu_uw, zfu
      REAL(wp), DIMENSION(jpi,jpj,jpk)   ::   zfv_t, zfv_f, zfv_vw, zfv, zfw
      REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   zlu_uu, zlu_uv
      REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   zlv_vv, zlv_vu
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_adv_ubs : UBS flux form momentum advection'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      zfu_t(:,:,:) = 0._wp
      zfv_t(:,:,:) = 0._wp
      zfu_f(:,:,:) = 0._wp
      zfv_f(:,:,:) = 0._wp
      !
      zlu_uu(:,:,:,:) = 0._wp
      zlv_vv(:,:,:,:) = 0._wp 
      zlu_uv(:,:,:,:) = 0._wp 
      zlv_vu(:,:,:,:) = 0._wp 
      !
      IF( l_trddyn ) THEN           ! trends: store the input trends
         zfu_uw(:,:,:) = ua(:,:,:)
         zfv_vw(:,:,:) = va(:,:,:)
      ENDIF
      !                                      ! =========================== !
      DO jk = 1, jpkm1                       !  Laplacian of the velocity  !
         !                                   ! =========================== !
         !                                         ! horizontal volume fluxes
         zfu(:,:,jk) = e2u(:,:) * e3u_n(:,:,jk) * un(:,:,jk)
         zfv(:,:,jk) = e1v(:,:) * e3v_n(:,:,jk) * vn(:,:,jk)
         !            
         DO jj = 2, jpjm1                          ! laplacian
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zlu_uu(ji,jj,jk,1) = ( ub (ji+1,jj  ,jk) - 2.*ub (ji,jj,jk) + ub (ji-1,jj  ,jk) ) * umask(ji,jj,jk)
               zlv_vv(ji,jj,jk,1) = ( vb (ji  ,jj+1,jk) - 2.*vb (ji,jj,jk) + vb (ji  ,jj-1,jk) ) * vmask(ji,jj,jk)
               zlu_uv(ji,jj,jk,1) = ( ub (ji  ,jj+1,jk) - ub (ji  ,jj  ,jk) ) * fmask(ji  ,jj  ,jk)   &
                  &               - ( ub (ji  ,jj  ,jk) - ub (ji  ,jj-1,jk) ) * fmask(ji  ,jj-1,jk)
               zlv_vu(ji,jj,jk,1) = ( vb (ji+1,jj  ,jk) - vb (ji  ,jj  ,jk) ) * fmask(ji  ,jj  ,jk)   &
                  &               - ( vb (ji  ,jj  ,jk) - vb (ji-1,jj  ,jk) ) * fmask(ji-1,jj  ,jk)
               !
               zlu_uu(ji,jj,jk,2) = ( zfu(ji+1,jj  ,jk) - 2.*zfu(ji,jj,jk) + zfu(ji-1,jj  ,jk) ) * umask(ji,jj,jk)
               zlv_vv(ji,jj,jk,2) = ( zfv(ji  ,jj+1,jk) - 2.*zfv(ji,jj,jk) + zfv(ji  ,jj-1,jk) ) * vmask(ji,jj,jk)
               zlu_uv(ji,jj,jk,2) = ( zfu(ji  ,jj+1,jk) - zfu(ji  ,jj  ,jk) ) * fmask(ji  ,jj  ,jk)   &
                  &               - ( zfu(ji  ,jj  ,jk) - zfu(ji  ,jj-1,jk) ) * fmask(ji  ,jj-1,jk)
               zlv_vu(ji,jj,jk,2) = ( zfv(ji+1,jj  ,jk) - zfv(ji  ,jj  ,jk) ) * fmask(ji  ,jj  ,jk)   &
                  &               - ( zfv(ji  ,jj  ,jk) - zfv(ji-1,jj  ,jk) ) * fmask(ji-1,jj  ,jk)
            END DO
         END DO
      END DO
      CALL lbc_lnk_multi( 'dynadv_ubs', zlu_uu(:,:,:,1), 'U', 1. , zlu_uv(:,:,:,1), 'U', 1.,  &
                      &   zlu_uu(:,:,:,2), 'U', 1. , zlu_uv(:,:,:,2), 'U', 1.,  & 
                      &   zlv_vv(:,:,:,1), 'V', 1. , zlv_vu(:,:,:,1), 'V', 1.,  &
                      &   zlv_vv(:,:,:,2), 'V', 1. , zlv_vu(:,:,:,2), 'V', 1.   )
      !
      !                                      ! ====================== !
      !                                      !  Horizontal advection  !
      DO jk = 1, jpkm1                       ! ====================== !
         !                                         ! horizontal volume fluxes
         zfu(:,:,jk) = 0.25_wp * e2u(:,:) * e3u_n(:,:,jk) * un(:,:,jk)
         zfv(:,:,jk) = 0.25_wp * e1v(:,:) * e3v_n(:,:,jk) * vn(:,:,jk)
         !
         DO jj = 1, jpjm1                          ! horizontal momentum fluxes at T- and F-point
            DO ji = 1, fs_jpim1   ! vector opt.
               zui = ( un(ji,jj,jk) + un(ji+1,jj  ,jk) )
               zvj = ( vn(ji,jj,jk) + vn(ji  ,jj+1,jk) )
               !
               IF( zui > 0 ) THEN   ;   zl_u = zlu_uu(ji  ,jj,jk,1)
               ELSE                 ;   zl_u = zlu_uu(ji+1,jj,jk,1)
               ENDIF
               IF( zvj > 0 ) THEN   ;   zl_v = zlv_vv(ji,jj  ,jk,1)
               ELSE                 ;   zl_v = zlv_vv(ji,jj+1,jk,1)
               ENDIF
               !
               zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj  ,jk)                               &
                  &                    - gamma2 * ( zlu_uu(ji,jj,jk,2) + zlu_uu(ji+1,jj  ,jk,2) )  )   &
                  &                * ( zui - gamma1 * zl_u)
               zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji  ,jj+1,jk)                               &
                  &                    - gamma2 * ( zlv_vv(ji,jj,jk,2) + zlv_vv(ji  ,jj+1,jk,2) )  )   &
                  &                * ( zvj - gamma1 * zl_v)
               !
               zfuj = ( zfu(ji,jj,jk) + zfu(ji  ,jj+1,jk) )
               zfvi = ( zfv(ji,jj,jk) + zfv(ji+1,jj  ,jk) )
               IF( zfuj > 0 ) THEN   ;    zl_v = zlv_vu( ji  ,jj  ,jk,1)
               ELSE                  ;    zl_v = zlv_vu( ji+1,jj,jk,1)
               ENDIF
               IF( zfvi > 0 ) THEN   ;    zl_u = zlu_uv( ji,jj  ,jk,1)
               ELSE                  ;    zl_u = zlu_uv( ji,jj+1,jk,1)
               ENDIF
               !
               zfv_f(ji  ,jj  ,jk) = ( zfvi - gamma2 * ( zlv_vu(ji,jj,jk,2) + zlv_vu(ji+1,jj  ,jk,2) )  )   &
                  &                * ( un(ji,jj,jk) + un(ji  ,jj+1,jk) - gamma1 * zl_u )
               zfu_f(ji  ,jj  ,jk) = ( zfuj - gamma2 * ( zlu_uv(ji,jj,jk,2) + zlu_uv(ji  ,jj+1,jk,2) )  )   &
                  &                * ( vn(ji,jj,jk) + vn(ji+1,jj  ,jk) - gamma1 * zl_v )
            END DO
         END DO
         DO jj = 2, jpjm1                          ! divergence of horizontal momentum fluxes
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - (  zfu_t(ji+1,jj,jk) - zfu_t(ji,jj  ,jk)    &
                  &                           + zfv_f(ji  ,jj,jk) - zfv_f(ji,jj-1,jk)  ) * r1_e1e2u(ji,jj) / e3u_n(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) - (  zfu_f(ji,jj  ,jk) - zfu_f(ji-1,jj,jk)    &
                  &                           + zfv_t(ji,jj+1,jk) - zfv_t(ji  ,jj,jk)  ) * r1_e1e2v(ji,jj) / e3v_n(ji,jj,jk)
            END DO
         END DO
      END DO
      IF( l_trddyn ) THEN                          ! trends: send trends to trddyn for diagnostic
         zfu_uw(:,:,:) = ua(:,:,:) - zfu_uw(:,:,:)
         zfv_vw(:,:,:) = va(:,:,:) - zfv_vw(:,:,:)
         CALL trd_dyn( zfu_uw, zfv_vw, jpdyn_keg, kt )
         zfu_t(:,:,:) = ua(:,:,:)
         zfv_t(:,:,:) = va(:,:,:)
      ENDIF
      !                                      ! ==================== !
      !                                      !  Vertical advection  !
      !                                      ! ==================== !
      DO jj = 2, jpjm1                             ! surface/bottom advective fluxes set to zero                  
         DO ji = fs_2, fs_jpim1
            zfu_uw(ji,jj,jpk) = 0._wp
            zfv_vw(ji,jj,jpk) = 0._wp
            zfu_uw(ji,jj, 1 ) = 0._wp
            zfv_vw(ji,jj, 1 ) = 0._wp
         END DO
      END DO
      IF( ln_linssh ) THEN                         ! constant volume : advection through the surface
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zfu_uw(ji,jj,1) = 0.5_wp * ( e1e2t(ji,jj) * wn(ji,jj,1) + e1e2t(ji+1,jj) * wn(ji+1,jj,1) ) * un(ji,jj,1)
               zfv_vw(ji,jj,1) = 0.5_wp * ( e1e2t(ji,jj) * wn(ji,jj,1) + e1e2t(ji,jj+1) * wn(ji,jj+1,1) ) * vn(ji,jj,1)
            END DO
         END DO
      ENDIF
      DO jk = 2, jpkm1                          ! interior fluxes
         DO jj = 2, jpj
            DO ji = 2, jpi
               zfw(ji,jj,jk) = 0.25_wp * e1e2t(ji,jj) * wn(ji,jj,jk)
            END DO
         END DO
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zfu_uw(ji,jj,jk) = ( zfw(ji,jj,jk)+ zfw(ji+1,jj,jk) ) * ( un(ji,jj,jk) + un(ji,jj,jk-1) )
               zfv_vw(ji,jj,jk) = ( zfw(ji,jj,jk)+ zfw(ji,jj+1,jk) ) * ( vn(ji,jj,jk) + vn(ji,jj,jk-1) )
            END DO
         END DO
      END DO
      DO jk = 1, jpkm1                          ! divergence of vertical momentum flux divergence
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) =  ua(ji,jj,jk) - ( zfu_uw(ji,jj,jk) - zfu_uw(ji,jj,jk+1) ) * r1_e1e2u(ji,jj) / e3u_n(ji,jj,jk)
               va(ji,jj,jk) =  va(ji,jj,jk) - ( zfv_vw(ji,jj,jk) - zfv_vw(ji,jj,jk+1) ) * r1_e1e2v(ji,jj) / e3v_n(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      IF( l_trddyn ) THEN                       ! save the vertical advection trend for diagnostic
         zfu_t(:,:,:) = ua(:,:,:) - zfu_t(:,:,:)
         zfv_t(:,:,:) = va(:,:,:) - zfv_t(:,:,:)
         CALL trd_dyn( zfu_t, zfv_t, jpdyn_zad, kt )
      ENDIF
      !                                         ! Control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' ubs2 adv - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=           ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
   END SUBROUTINE dyn_adv_ubs

   !!==============================================================================
END MODULE dynadv_ubs
