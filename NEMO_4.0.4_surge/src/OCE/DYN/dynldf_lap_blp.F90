MODULE dynldf_lap_blp
   !!======================================================================
   !!                   ***  MODULE  dynldf_lap_blp  ***
   !! Ocean dynamics:  lateral viscosity trend (laplacian and bilaplacian)
   !!======================================================================
   !! History : 3.7  ! 2014-01  (G. Madec, S. Masson)  Original code, re-entrant laplacian
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_ldf_lap   : update the momentum trend with the lateral viscosity using an iso-level   laplacian operator
   !!   dyn_ldf_blp   : update the momentum trend with the lateral viscosity using an iso-level bilaplacian operator
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE ldfdyn         ! lateral diffusion: eddy viscosity coef.
   USE ldfslp         ! iso-neutral slopes 
   USE zdf_oce        ! ocean vertical physics
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_ldf_lap  ! called by dynldf.F90
   PUBLIC dyn_ldf_blp  ! called by dynldf.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_ldf_lap( kt, pub, pvb, pua, pva, kpass )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_lap  ***
      !!                       
      !! ** Purpose :   Compute the before horizontal momentum diffusive 
      !!      trend and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The Laplacian operator apply on horizontal velocity is 
      !!      writen as :   grad_h( ahmt div_h(U )) - curl_h( ahmf curl_z(U) ) 
      !!
      !! ** Action : - pua, pva increased by the harmonic operator applied on pub, pvb.
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   kpass      ! =1/2 first or second passage
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pub, pvb   ! before velocity  [m/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pua, pva   ! velocity trend   [m/s2]
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zsign        ! local scalars
      REAL(wp) ::   zua, zva     ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   zcur, zdiv
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_ldf : iso-level harmonic (laplacian) operator, pass=', kpass
         WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !
      IF( kpass == 1 ) THEN   ;   zsign =  1._wp      ! bilaplacian operator require a minus sign
      ELSE                    ;   zsign = -1._wp      !  (eddy viscosity coef. >0)
      ENDIF
      !
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO jj = 2, jpj
            DO ji = fs_2, jpi   ! vector opt.
               !                                      ! ahm * e3 * curl  (computed from 1 to jpim1/jpjm1)
               zcur(ji-1,jj-1) = ahmf(ji-1,jj-1,jk) * e3f_n(ji-1,jj-1,jk) * r1_e1e2f(ji-1,jj-1)       &   ! ahmf already * by fmask
                  &     * (  e2v(ji  ,jj-1) * pvb(ji  ,jj-1,jk) - e2v(ji-1,jj-1) * pvb(ji-1,jj-1,jk)  &
                  &        - e1u(ji-1,jj  ) * pub(ji-1,jj  ,jk) + e1u(ji-1,jj-1) * pub(ji-1,jj-1,jk)  )
               !                                      ! ahm * div        (computed from 2 to jpi/jpj)
               zdiv(ji,jj)     = ahmt(ji,jj,jk) * r1_e1e2t(ji,jj) / e3t_b(ji,jj,jk)                   &   ! ahmt already * by tmask
                  &     * (  e2u(ji,jj)*e3u_b(ji,jj,jk) * pub(ji,jj,jk) - e2u(ji-1,jj)*e3u_b(ji-1,jj,jk) * pub(ji-1,jj,jk)  &
                  &        + e1v(ji,jj)*e3v_b(ji,jj,jk) * pvb(ji,jj,jk) - e1v(ji,jj-1)*e3v_b(ji,jj-1,jk) * pvb(ji,jj-1,jk)  )
            END DO  
         END DO  
         !
         DO jj = 2, jpjm1                             ! - curl( curl) + grad( div )
            DO ji = fs_2, fs_jpim1   ! vector opt.
               pua(ji,jj,jk) = pua(ji,jj,jk) + zsign * umask(ji,jj,jk) * (       &    ! * by umask is mandatory for dyn_ldf_blp use
                  &              - ( zcur(ji  ,jj) - zcur(ji,jj-1) ) * r1_e2u(ji,jj) / e3u_n(ji,jj,jk)   &
                  &              + ( zdiv(ji+1,jj) - zdiv(ji,jj  ) ) * r1_e1u(ji,jj)                     )
                  !
               pva(ji,jj,jk) = pva(ji,jj,jk) + zsign * vmask(ji,jj,jk) * (       &    ! * by vmask is mandatory for dyn_ldf_blp use
                  &                ( zcur(ji,jj  ) - zcur(ji-1,jj) ) * r1_e1v(ji,jj) / e3v_n(ji,jj,jk)   &
                  &              + ( zdiv(ji,jj+1) - zdiv(ji  ,jj) ) * r1_e2v(ji,jj)                     )
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      !
   END SUBROUTINE dyn_ldf_lap


   SUBROUTINE dyn_ldf_blp( kt, pub, pvb, pua, pva )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dyn_ldf_blp  ***
      !!                    
      !! ** Purpose :   Compute the before lateral momentum viscous trend 
      !!              and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The lateral viscous trends is provided by a bilaplacian
      !!      operator applied to before field (forward in time).
      !!      It is computed by two successive calls to dyn_ldf_lap routine
      !!
      !! ** Action :   pta   updated with the before rotated bilaplacian diffusion
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt         ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pub, pvb   ! before velocity fields
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pua, pva   ! momentum trend
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zulap, zvlap   ! laplacian at u- and v-point
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_ldf_blp : bilaplacian operator momentum '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      zulap(:,:,:) = 0._wp
      zvlap(:,:,:) = 0._wp
      !
      CALL dyn_ldf_lap( kt, pub, pvb, zulap, zvlap, 1 )   ! rotated laplacian applied to ptb (output in zlap)
      !
      CALL lbc_lnk_multi( 'dynldf_lap_blp', zulap, 'U', -1., zvlap, 'V', -1. )             ! Lateral boundary conditions
      !
      CALL dyn_ldf_lap( kt, zulap, zvlap, pua, pva, 2 )   ! rotated laplacian applied to zlap (output in pta)
      !
   END SUBROUTINE dyn_ldf_blp

   !!======================================================================
END MODULE dynldf_lap_blp
