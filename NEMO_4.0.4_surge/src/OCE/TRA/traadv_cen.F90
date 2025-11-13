MODULE traadv_cen
   !!======================================================================
   !!                     ***  MODULE  traadv_cen  ***
   !! Ocean  tracers:   advective trend (2nd/4th order centered)
   !!======================================================================
   !! History :  3.7  ! 2014-05  (G. Madec)  original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv_cen   : update the tracer trend with the advection trends using a centered or scheme (2nd or 4th order)
   !!                   NB: on the vertical it is actually a 4th order COMPACT scheme which is used
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE eosbn2         ! equation of state
   USE traadv_fct     ! acces to routine interp_4th_cpt 
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers 
   USE diaptr         ! poleward transport diagnostics
   USE diaar5         ! AR5 diagnostics
   !
   USE in_out_manager ! I/O manager
   USE iom            ! IOM library
   USE trc_oce        ! share passive tracers/Ocean variables
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_cen   ! called by traadv.F90
   
   REAL(wp) ::   r1_6 = 1._wp / 6._wp   ! =1/6

   LOGICAL ::   l_trd   ! flag to compute trends
   LOGICAL ::   l_ptr   ! flag to compute poleward transport
   LOGICAL ::   l_hst   ! flag to compute heat/salt transport

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv_cen( kt, kit000, cdtype, pun, pvn, pwn,     &
      &                                        ptn, pta, kjpt, kn_cen_h, kn_cen_v ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_cen  ***
      !!                 
      !! ** Purpose :   Compute the now trend due to the advection of tracers
      !!      and add it to the general trend of passive tracer equations.
      !!
      !! ** Method  :   The advection is evaluated by a 2nd or 4th order scheme
      !!               using now fields (leap-frog scheme). 
      !!       kn_cen_h = 2  ==>> 2nd order centered scheme on the horizontal
      !!                = 4  ==>> 4th order    -        -       -      -
      !!       kn_cen_v = 2  ==>> 2nd order centered scheme on the vertical
      !!                = 4  ==>> 4th order COMPACT  scheme     -      -
      !!
      !! ** Action : - update pta  with the now advective tracer trends
      !!             - send trends to trdtra module for further diagnostcs (l_trdtra=T)
      !!             - htr_adv, str_adv : poleward advective heat and salt transport (ln_diaptr=T)
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt              ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt            ! number of tracers
      INTEGER                              , INTENT(in   ) ::   kn_cen_h        ! =2/4 (2nd or 4th order scheme)
      INTEGER                              , INTENT(in   ) ::   kn_cen_v        ! =2/4 (2nd or 4th order scheme)
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun, pvn, pwn   ! 3 ocean velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptn             ! now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta             ! tracer trend 
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      INTEGER  ::   ierr             ! local integer
      REAL(wp) ::   zC2t_u, zC4t_u   ! local scalars
      REAL(wp) ::   zC2t_v, zC4t_v   !   -      -
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zwx, zwy, zwz, ztu, ztv, ztw
      !!----------------------------------------------------------------------
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_cen : centered advection scheme on ', cdtype, ' order h/v =', kn_cen_h,'/', kn_cen_v
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~ '
      ENDIF
      !                          ! set local switches
      l_trd = .FALSE.
      l_hst = .FALSE.
      l_ptr = .FALSE.
      IF( ( cdtype == 'TRA' .AND. l_trdtra ) .OR. ( cdtype == 'TRC' .AND. l_trdtrc ) )        l_trd = .TRUE.
      IF(   cdtype == 'TRA' .AND. ln_diaptr )                                                 l_ptr = .TRUE. 
      IF(   cdtype == 'TRA' .AND. ( iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. &
         &                          iom_use("uadv_salttr") .OR. iom_use("vadv_salttr")  ) )   l_hst = .TRUE.
      !
      !                    
      zwz(:,:, 1 ) = 0._wp       ! surface & bottom vertical flux set to zero for all tracers
      zwz(:,:,jpk) = 0._wp
      !
      DO jn = 1, kjpt            !==  loop over the tracers  ==!
         !
         SELECT CASE( kn_cen_h )       !--  Horizontal fluxes  --!
         !
         CASE(  2  )                         !* 2nd order centered
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.
                     zwx(ji,jj,jk) = 0.5_wp * pun(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji+1,jj  ,jk,jn) )
                     zwy(ji,jj,jk) = 0.5_wp * pvn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji  ,jj+1,jk,jn) )
                  END DO
               END DO
            END DO
            !
         CASE(  4  )                         !* 4th order centered
            ztu(:,:,jpk) = 0._wp                   ! Bottom value : flux set to zero
            ztv(:,:,jpk) = 0._wp
            DO jk = 1, jpkm1                       ! masked gradient
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     ztu(ji,jj,jk) = ( ptn(ji+1,jj  ,jk,jn) - ptn(ji,jj,jk,jn) ) * umask(ji,jj,jk)
                     ztv(ji,jj,jk) = ( ptn(ji  ,jj+1,jk,jn) - ptn(ji,jj,jk,jn) ) * vmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            CALL lbc_lnk_multi( 'traadv_cen', ztu, 'U', -1. , ztv, 'V', -1. )   ! Lateral boundary cond.
            !
            DO jk = 1, jpkm1                       ! Horizontal advective fluxes
               DO jj = 2, jpjm1
                  DO ji = 2, fs_jpim1   ! vector opt.
                     zC2t_u = ptn(ji,jj,jk,jn) + ptn(ji+1,jj  ,jk,jn)   ! C2 interpolation of T at u- & v-points (x2)
                     zC2t_v = ptn(ji,jj,jk,jn) + ptn(ji  ,jj+1,jk,jn)
                     !                                                  ! C4 interpolation of T at u- & v-points (x2)
                     zC4t_u =  zC2t_u + r1_6 * ( ztu(ji-1,jj,jk) - ztu(ji+1,jj,jk) )
                     zC4t_v =  zC2t_v + r1_6 * ( ztv(ji,jj-1,jk) - ztv(ji,jj+1,jk) )
                     !                                                  ! C4 fluxes
                     zwx(ji,jj,jk) =  0.5_wp * pun(ji,jj,jk) * zC4t_u
                     zwy(ji,jj,jk) =  0.5_wp * pvn(ji,jj,jk) * zC4t_v
                  END DO
               END DO
            END DO         
            CALL lbc_lnk_multi( 'traadv_cen', zwx, 'U', -1. , zwy, 'V', -1. )
            !
         CASE DEFAULT
            CALL ctl_stop( 'traadv_cen: wrong value for nn_cen' )
         END SELECT
         !
         SELECT CASE( kn_cen_v )       !--  Vertical fluxes  --!   (interior)
         !
         CASE(  2  )                         !* 2nd order centered
            DO jk = 2, jpk
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     zwz(ji,jj,jk) = 0.5 * pwn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji,jj,jk-1,jn) ) * wmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
         CASE(  4  )                         !* 4th order compact
            CALL interp_4th_cpt( ptn(:,:,:,jn) , ztw )      ! ztw = interpolated value of T at w-point
            DO jk = 2, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     zwz(ji,jj,jk) = pwn(ji,jj,jk) * ztw(ji,jj,jk) * wmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
         END SELECT
         !
         IF( ln_linssh ) THEN                !* top value   (linear free surf. only as zwz is multiplied by wmask)
            IF( ln_isfcav ) THEN                  ! ice-shelf cavities (top of the ocean)
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zwz(ji,jj, mikt(ji,jj) ) = pwn(ji,jj,mikt(ji,jj)) * ptn(ji,jj,mikt(ji,jj),jn) 
                  END DO
               END DO   
            ELSE                                   ! no ice-shelf cavities (only ocean surface)
               zwz(:,:,1) = pwn(:,:,1) * ptn(:,:,1,jn)
            ENDIF
         ENDIF
         !               
         DO jk = 1, jpkm1              !--  Divergence of advective fluxes  --!
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn)    &
                     &             - (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )    &
                     &                + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  )    &
                     &                + zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1)  ) * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
               END DO
            END DO
         END DO
         !                             ! trend diagnostics
         IF( l_trd ) THEN
            CALL trd_tra( kt, cdtype, jn, jptra_xad, zwx, pun, ptn(:,:,:,jn) )
            CALL trd_tra( kt, cdtype, jn, jptra_yad, zwy, pvn, ptn(:,:,:,jn) )
            CALL trd_tra( kt, cdtype, jn, jptra_zad, zwz, pwn, ptn(:,:,:,jn) )
         END IF
         !                                 ! "Poleward" heat and salt transports 
         IF( l_ptr )   CALL dia_ptr_hst( jn, 'adv', zwy(:,:,:) )
         !                                 !  heat and salt transport
         IF( l_hst )   CALL dia_ar5_hst( jn, 'adv', zwx(:,:,:), zwy(:,:,:) )
         !
      END DO
      !
   END SUBROUTINE tra_adv_cen
   
   !!======================================================================
END MODULE traadv_cen
