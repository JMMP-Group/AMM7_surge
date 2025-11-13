MODULE traadv_ubs
   !!==============================================================================
   !!                       ***  MODULE  traadv_ubs  ***
   !! Ocean active tracers:  horizontal & vertical advective trend
   !!==============================================================================
   !! History :  1.0  !  2006-08  (L. Debreu, R. Benshila)  Original code
   !!            3.3  !  2010-05  (C. Ethe, G. Madec)  merge TRC-TRA + switch from velocity to transport
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv_ubs : update the tracer trend with the horizontal
   !!                 advection trends using a third order biaised scheme  
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE dom_oce        ! ocean space and time domain
   USE trc_oce        ! share passive tracers/Ocean variables
   USE trd_oce        ! trends: ocean variables
   USE traadv_fct      ! acces to routine interp_4th_cpt 
   USE trdtra         ! trends manager: tracers 
   USE diaptr         ! poleward transport diagnostics
   USE diaar5         ! AR5 diagnostics
   !
   USE iom            ! I/O library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! massively parallel library
   USE lbclnk         ! ocean lateral boundary condition (or mpp link)
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_ubs   ! routine called by traadv module

   LOGICAL :: l_trd   ! flag to compute trends
   LOGICAL :: l_ptr   ! flag to compute poleward transport
   LOGICAL :: l_hst   ! flag to compute heat transport


   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv_ubs( kt, kit000, cdtype, p2dt, pun, pvn, pwn,          &
      &                                                ptb, ptn, pta, kjpt, kn_ubs_v )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_ubs  ***
      !!                 
      !! ** Purpose :   Compute the now trend due to the advection of tracers
      !!      and add it to the general trend of passive tracer equations.
      !!
      !! ** Method  :   The 3rd order Upstream Biased Scheme (UBS) is based on an
      !!      upstream-biased parabolic interpolation (Shchepetkin and McWilliams 2005)
      !!      It is only used in the horizontal direction.
      !!      For example the i-component of the advective fluxes are given by :
      !!                !  e2u e3u un ( mi(Tn) - zltu(i  ) )   if un(i) >= 0
      !!          ztu = !  or 
      !!                !  e2u e3u un ( mi(Tn) - zltu(i+1) )   if un(i) < 0
      !!      where zltu is the second derivative of the before temperature field:
      !!          zltu = 1/e3t di[ e2u e3u / e1u di[Tb] ]
      !!        This results in a dissipatively dominant (i.e. hyper-diffusive) 
      !!      truncation error. The overall performance of the advection scheme 
      !!      is similar to that reported in (Farrow and Stevens, 1995). 
      !!        For stability reasons, the first term of the fluxes which corresponds
      !!      to a second order centered scheme is evaluated using the now velocity 
      !!      (centered in time) while the second term which is the diffusive part 
      !!      of the scheme, is evaluated using the before velocity (forward in time). 
      !!      Note that UBS is not positive. Do not use it on passive tracers.
      !!                On the vertical, the advection is evaluated using a FCT scheme,
      !!      as the UBS have been found to be too diffusive. 
      !!                kn_ubs_v argument controles whether the FCT is based on 
      !!      a 2nd order centrered scheme (kn_ubs_v=2) or on a 4th order compact 
      !!      scheme (kn_ubs_v=4).
      !!
      !! ** Action : - update pta  with the now advective tracer trends
      !!             - send trends to trdtra module for further diagnostcs (l_trdtra=T)
      !!             - htr_adv, str_adv : poleward advective heat and salt transport (ln_diaptr=T)
      !!
      !! Reference : Shchepetkin, A. F., J. C. McWilliams, 2005, Ocean Modelling, 9, 347-404. 
      !!             Farrow, D.E., Stevens, D.P., 1995, J. Phys. Ocean. 25, 1731Ð1741. 
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt              ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt            ! number of tracers
      INTEGER                              , INTENT(in   ) ::   kn_ubs_v        ! number of tracers
      REAL(wp)                             , INTENT(in   ) ::   p2dt            ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun, pvn, pwn   ! 3 ocean transport components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb, ptn        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta             ! tracer trend 
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::   ztra, zbtr, zcoef                       ! local scalars
      REAL(wp) ::   zfp_ui, zfm_ui, zcenut, ztak, zfp_wk, zfm_wk   !   -      -
      REAL(wp) ::   zfp_vj, zfm_vj, zcenvt, zeeu, zeev, z_hdivn    !   -      -
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ztu, ztv, zltu, zltv, zti, ztw   ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_ubs :  horizontal UBS advection scheme on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      l_trd = .FALSE.
      l_hst = .FALSE.
      l_ptr = .FALSE.
      IF( ( cdtype == 'TRA' .AND. l_trdtra ) .OR. ( cdtype == 'TRC' .AND. l_trdtrc ) )      l_trd = .TRUE.
      IF(   cdtype == 'TRA' .AND. ln_diaptr )                                               l_ptr = .TRUE. 
      IF(   cdtype == 'TRA' .AND. ( iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. &
         &                          iom_use("uadv_salttr") .OR. iom_use("vadv_salttr")  ) ) l_hst = .TRUE.
      !
      ztw (:,:, 1 ) = 0._wp      ! surface & bottom value : set to zero for all tracers
      zltu(:,:,jpk) = 0._wp   ;   zltv(:,:,jpk) = 0._wp
      ztw (:,:,jpk) = 0._wp   ;   zti (:,:,jpk) = 0._wp
      !
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         !                                              
         DO jk = 1, jpkm1        !==  horizontal laplacian of before tracer ==!
            DO jj = 1, jpjm1              ! First derivative (masked gradient)
               DO ji = 1, fs_jpim1   ! vector opt.
                  zeeu = e2_e1u(ji,jj) * e3u_n(ji,jj,jk) * umask(ji,jj,jk)
                  zeev = e1_e2v(ji,jj) * e3v_n(ji,jj,jk) * vmask(ji,jj,jk)
                  ztu(ji,jj,jk) = zeeu * ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) )
                  ztv(ji,jj,jk) = zeev * ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) )
               END DO
            END DO
            DO jj = 2, jpjm1              ! Second derivative (divergence)
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zcoef = 1._wp / ( 6._wp * e3t_n(ji,jj,jk) )
                  zltu(ji,jj,jk) = (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk)  ) * zcoef
                  zltv(ji,jj,jk) = (  ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  ) * zcoef
               END DO
            END DO
            !                                    
         END DO         
         CALL lbc_lnk( 'traadv_ubs', zltu, 'T', 1. )   ;    CALL lbc_lnk( 'traadv_ubs', zltv, 'T', 1. )   ! Lateral boundary cond. (unchanged sgn)
         !    
         DO jk = 1, jpkm1        !==  Horizontal advective fluxes  ==!     (UBS)
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zfp_ui = pun(ji,jj,jk) + ABS( pun(ji,jj,jk) )      ! upstream transport (x2)
                  zfm_ui = pun(ji,jj,jk) - ABS( pun(ji,jj,jk) )
                  zfp_vj = pvn(ji,jj,jk) + ABS( pvn(ji,jj,jk) )
                  zfm_vj = pvn(ji,jj,jk) - ABS( pvn(ji,jj,jk) )
                  !                                                  ! 2nd order centered advective fluxes (x2)
                  zcenut = pun(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji+1,jj  ,jk,jn) )
                  zcenvt = pvn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji  ,jj+1,jk,jn) )
                  !                                                  ! UBS advective fluxes
                  ztu(ji,jj,jk) = 0.5 * ( zcenut - zfp_ui * zltu(ji,jj,jk) - zfm_ui * zltu(ji+1,jj,jk) )
                  ztv(ji,jj,jk) = 0.5 * ( zcenvt - zfp_vj * zltv(ji,jj,jk) - zfm_vj * zltv(ji,jj+1,jk) )
               END DO
            END DO
         END DO         
         !
         zltu(:,:,:) = pta(:,:,:,jn)      ! store the initial trends before its update
         !
         DO jk = 1, jpkm1        !==  add the horizontal advective trend  ==!
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn)                        &
                     &             - (  ztu(ji,jj,jk) - ztu(ji-1,jj  ,jk)    &
                     &                + ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk)  ) * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
               END DO
            END DO
            !                                             
         END DO
         !
         zltu(:,:,:) = pta(:,:,:,jn) - zltu(:,:,:)    ! Horizontal advective trend used in vertical 2nd order FCT case
         !                                            ! and/or in trend diagnostic (l_trd=T) 
         !                
         IF( l_trd ) THEN                  ! trend diagnostics
             CALL trd_tra( kt, cdtype, jn, jptra_xad, ztu, pun, ptn(:,:,:,jn) )
             CALL trd_tra( kt, cdtype, jn, jptra_yad, ztv, pvn, ptn(:,:,:,jn) )
         END IF
         !     
         !                                ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( l_ptr )  CALL dia_ptr_hst( jn, 'adv', ztv(:,:,:) )
         !                                !  heati/salt transport
         IF( l_hst )  CALL dia_ar5_hst( jn, 'adv', ztu(:,:,:), ztv(:,:,:) )
         !
         !
         !                       !== vertical advective trend  ==!
         !
         SELECT CASE( kn_ubs_v )       ! select the vertical advection scheme
         !
         CASE(  2  )                   ! 2nd order FCT 
            !         
            IF( l_trd )   zltv(:,:,:) = pta(:,:,:,jn)          ! store pta if trend diag.
            !
            !                          !*  upstream advection with initial mass fluxes & intermediate update  ==!
            DO jk = 2, jpkm1                 ! Interior value (w-masked)
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zfp_wk = pwn(ji,jj,jk) + ABS( pwn(ji,jj,jk) )
                     zfm_wk = pwn(ji,jj,jk) - ABS( pwn(ji,jj,jk) )
                     ztw(ji,jj,jk) = 0.5_wp * (  zfp_wk * ptb(ji,jj,jk,jn) + zfm_wk * ptb(ji,jj,jk-1,jn)  ) * wmask(ji,jj,jk)
                  END DO
               END DO
            END DO 
            IF( ln_linssh ) THEN             ! top ocean value (only in linear free surface as ztw has been w-masked)
               IF( ln_isfcav ) THEN                ! top of the ice-shelf cavities and at the ocean surface
                  DO jj = 1, jpj
                     DO ji = 1, jpi
                        ztw(ji,jj, mikt(ji,jj) ) = pwn(ji,jj,mikt(ji,jj)) * ptb(ji,jj,mikt(ji,jj),jn)   ! linear free surface 
                     END DO
                  END DO   
               ELSE                                ! no cavities: only at the ocean surface
                  ztw(:,:,1) = pwn(:,:,1) * ptb(:,:,1,jn)
               ENDIF
            ENDIF
            !
            DO jk = 1, jpkm1           !* trend and after field with monotonic scheme
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     ztak = - ( ztw(ji,jj,jk) - ztw(ji,jj,jk+1) ) * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
                     pta(ji,jj,jk,jn) =   pta(ji,jj,jk,jn) +  ztak 
                     zti(ji,jj,jk)    = ( ptb(ji,jj,jk,jn) + p2dt * ( ztak + zltu(ji,jj,jk) ) ) * tmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            CALL lbc_lnk( 'traadv_ubs', zti, 'T', 1. )      ! Lateral boundary conditions on zti, zsi   (unchanged sign)
            !
            !                          !*  anti-diffusive flux : high order minus low order
            DO jk = 2, jpkm1        ! Interior value  (w-masked)
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     ztw(ji,jj,jk) = (   0.5_wp * pwn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji,jj,jk-1,jn) )   &
                        &              - ztw(ji,jj,jk)   ) * wmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            !                                            ! top ocean value: high order == upstream  ==>>  zwz=0
            IF( ln_linssh )   ztw(:,:, 1 ) = 0._wp       ! only ocean surface as interior zwz values have been w-masked
            !
            CALL nonosc_z( ptb(:,:,:,jn), ztw, zti, p2dt )      !  monotonicity algorithm
            !
         CASE(  4  )                               ! 4th order COMPACT
            CALL interp_4th_cpt( ptn(:,:,:,jn) , ztw )         ! 4th order compact interpolation of T at w-point
            DO jk = 2, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     ztw(ji,jj,jk) = pwn(ji,jj,jk) * ztw(ji,jj,jk) * wmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            IF( ln_linssh )   ztw(:,:, 1 ) = pwn(:,:,1) * ptn(:,:,1,jn)     !!gm ISF & 4th COMPACT doesn't work
            !
         END SELECT
         !
         DO jk = 1, jpkm1        !  final trend with corrected fluxes
            DO jj = 2, jpjm1 
               DO ji = fs_2, fs_jpim1   ! vector opt.   
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) - ( ztw(ji,jj,jk) - ztw(ji,jj,jk+1) ) * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
               END DO
            END DO
         END DO
         !
         IF( l_trd )  THEN       ! vertical advective trend diagnostics
            DO jk = 1, jpkm1                       ! (compute -w.dk[ptn]= -dk[w.ptn] + ptn.dk[w])
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     zltv(ji,jj,jk) = pta(ji,jj,jk,jn) - zltv(ji,jj,jk)                          &
                        &           + ptn(ji,jj,jk,jn) * (  pwn(ji,jj,jk) - pwn(ji,jj,jk+1)  )   &
                        &                              * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
                  END DO
               END DO
            END DO
            CALL trd_tra( kt, cdtype, jn, jptra_zad, zltv )
         ENDIF
         !
      END DO
      !
   END SUBROUTINE tra_adv_ubs


   SUBROUTINE nonosc_z( pbef, pcc, paft, p2dt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nonosc_z  ***
      !!     
      !! **  Purpose :   compute monotonic tracer fluxes from the upstream 
      !!       scheme and the before field by a nonoscillatory algorithm 
      !!
      !! **  Method  :   ... ???
      !!       warning : pbef and paft must be masked, but the boundaries
      !!       conditions on the fluxes are not necessary zalezak (1979)
      !!       drange (1995) multi-dimensional forward-in-time and upstream-
      !!       in-space based differencing for fluid
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )                          ::   p2dt   ! tracer time-step
      REAL(wp),                DIMENSION (jpi,jpj,jpk) ::   pbef   ! before field
      REAL(wp), INTENT(inout), DIMENSION (jpi,jpj,jpk) ::   paft   ! after field
      REAL(wp), INTENT(inout), DIMENSION (jpi,jpj,jpk) ::   pcc    ! monotonic flux in the k direction
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikm1         ! local integer
      REAL(wp) ::   zpos, zneg, zbt, za, zb, zc, zbig, zrtrn   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zbetup, zbetdo     ! 3D workspace
      !!----------------------------------------------------------------------
      !
      zbig  = 1.e+40_wp
      zrtrn = 1.e-15_wp
      zbetup(:,:,:) = 0._wp   ;   zbetdo(:,:,:) = 0._wp
      !
      ! Search local extrema
      ! --------------------
      !                    ! large negative value (-zbig) inside land
      pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:) - zbig * ( 1.e0 - tmask(:,:,:) )
      paft(:,:,:) = paft(:,:,:) * tmask(:,:,:) - zbig * ( 1.e0 - tmask(:,:,:) )
      !
      DO jk = 1, jpkm1     ! search maximum in neighbourhood
         ikm1 = MAX(jk-1,1)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zbetup(ji,jj,jk) = MAX(  pbef(ji  ,jj  ,jk  ), paft(ji  ,jj  ,jk  ),   &
                  &                     pbef(ji  ,jj  ,ikm1), pbef(ji  ,jj  ,jk+1),   &
                  &                     paft(ji  ,jj  ,ikm1), paft(ji  ,jj  ,jk+1)  )
            END DO
         END DO
      END DO
      !                    ! large positive value (+zbig) inside land
      pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:) + zbig * ( 1.e0 - tmask(:,:,:) )
      paft(:,:,:) = paft(:,:,:) * tmask(:,:,:) + zbig * ( 1.e0 - tmask(:,:,:) )
      !
      DO jk = 1, jpkm1     ! search minimum in neighbourhood
         ikm1 = MAX(jk-1,1)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zbetdo(ji,jj,jk) = MIN(  pbef(ji  ,jj  ,jk  ), paft(ji  ,jj  ,jk  ),   &
                  &                     pbef(ji  ,jj  ,ikm1), pbef(ji  ,jj  ,jk+1),   &
                  &                     paft(ji  ,jj  ,ikm1), paft(ji  ,jj  ,jk+1)  )
            END DO
         END DO
      END DO
      !                    ! restore masked values to zero
      pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:)
      paft(:,:,:) = paft(:,:,:) * tmask(:,:,:)
      !
      ! Positive and negative part of fluxes and beta terms
      ! ---------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! positive & negative part of the flux
               zpos = MAX( 0., pcc(ji  ,jj  ,jk+1) ) - MIN( 0., pcc(ji  ,jj  ,jk  ) )
               zneg = MAX( 0., pcc(ji  ,jj  ,jk  ) ) - MIN( 0., pcc(ji  ,jj  ,jk+1) )
               ! up & down beta terms
               zbt = e1e2t(ji,jj) * e3t_n(ji,jj,jk) / p2dt
               zbetup(ji,jj,jk) = ( zbetup(ji,jj,jk) - paft(ji,jj,jk) ) / (zpos+zrtrn) * zbt
               zbetdo(ji,jj,jk) = ( paft(ji,jj,jk) - zbetdo(ji,jj,jk) ) / (zneg+zrtrn) * zbt
            END DO
         END DO
      END DO
      !
      ! monotonic flux in the k direction, i.e. pcc
      ! -------------------------------------------
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               za = MIN( 1., zbetdo(ji,jj,jk), zbetup(ji,jj,jk-1) )
               zb = MIN( 1., zbetup(ji,jj,jk), zbetdo(ji,jj,jk-1) )
               zc = 0.5 * ( 1.e0 + SIGN( 1.e0, pcc(ji,jj,jk) ) )
               pcc(ji,jj,jk) = pcc(ji,jj,jk) * ( zc * za + ( 1.e0 - zc) * zb )
            END DO
         END DO
      END DO
      !
   END SUBROUTINE nonosc_z

   !!======================================================================
END MODULE traadv_ubs
