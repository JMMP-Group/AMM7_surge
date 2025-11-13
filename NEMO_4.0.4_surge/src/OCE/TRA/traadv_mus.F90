MODULE traadv_mus
   !!======================================================================
   !!                       ***  MODULE  traadv_mus  ***
   !! Ocean  tracers:  horizontal & vertical advective trend
   !!======================================================================
   !! History :       !  2000-06  (A.Estublier)  for passive tracers
   !!                 !  2001-08  (E.Durand, G.Madec)  adapted for T & S
   !!   NEMO     1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!            3.2  !  2010-05  (C. Ethe, G. Madec)  merge TRC-TRA + switch from velocity to transport
   !!            3.4  !  2012-06  (P. Oddo, M. Vichi) include the upstream where needed
   !!            3.7  !  2015-09  (G. Madec) add the ice-shelf cavities boundary condition
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv_mus   : update the tracer trend with the horizontal
   !!                   and vertical advection trends using MUSCL scheme
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE trc_oce        ! share passive tracers/Ocean variables
   USE dom_oce        ! ocean space and time domain
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! tracers trends manager
   USE sbcrnf         ! river runoffs
   USE diaptr         ! poleward transport diagnostics
   USE diaar5         ! AR5 diagnostics

   !
   USE iom            ! XIOS library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing
   USE lbclnk         ! ocean lateral boundary condition (or mpp link) 
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_mus   ! routine called by traadv.F90
   
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   upsmsk   !: mixed upstream/centered scheme near some straits
   !                                                           !  and in closed seas (orca 2 and 1 configurations)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xind     !: mixed upstream/centered index
   
   LOGICAL  ::   l_trd   ! flag to compute trends
   LOGICAL  ::   l_ptr   ! flag to compute poleward transport
   LOGICAL  ::   l_hst   ! flag to compute heat/salt transport

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv_mus( kt, kit000, cdtype, p2dt, pun, pvn, pwn,             &
      &                                              ptb, pta, kjpt, ld_msc_ups )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE tra_adv_mus  ***
      !!
      !! ** Purpose :   Compute the now trend due to total advection of tracers
      !!              using a MUSCL scheme (Monotone Upstream-centered Scheme for
      !!              Conservation Laws) and add it to the general tracer trend.
      !!
      !! ** Method  : MUSCL scheme plus centered scheme at ocean boundaries
      !!              ld_msc_ups=T : 
      !!
      !! ** Action : - update pta  with the now advective tracer trends
      !!             - send trends to trdtra module for further diagnostcs (l_trdtra=T)
      !!             - htr_adv, str_adv : poleward advective heat and salt transport (ln_diaptr=T)
      !!
      !! References : Estubier, A., and M. Levy, Notes Techn. Pole de Modelisation
      !!              IPSL, Sept. 2000 (http://www.lodyc.jussieu.fr/opa)
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt              ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt            ! number of tracers
      LOGICAL                              , INTENT(in   ) ::   ld_msc_ups      ! use upstream scheme within muscl
      REAL(wp)                             , INTENT(in   ) ::   p2dt            ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun, pvn, pwn   ! 3 ocean velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb             ! before tracer field
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta             ! tracer trend 
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      INTEGER  ::   ierr             ! local integer
      REAL(wp) ::   zu, z0u, zzwx, zw , zalpha   ! local scalars
      REAL(wp) ::   zv, z0v, zzwy, z0w           !   -      -
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zwx, zslpx   ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zwy, zslpy   ! -      - 
      !!----------------------------------------------------------------------
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv : MUSCL advection scheme on ', cdtype
         IF(lwp) WRITE(numout,*) '        : mixed up-stream           ', ld_msc_ups
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         IF(lwp) WRITE(numout,*)
         !
         ! Upstream / MUSCL scheme indicator
         !
         ALLOCATE( xind(jpi,jpj,jpk), STAT=ierr )
         xind(:,:,:) = 1._wp              ! set equal to 1 where up-stream is not needed
         !
         IF( ld_msc_ups ) THEN            ! define the upstream indicator (if asked)
            ALLOCATE( upsmsk(jpi,jpj), STAT=ierr )
            upsmsk(:,:) = 0._wp                             ! not upstream by default
            !
            DO jk = 1, jpkm1
               xind(:,:,jk) = 1._wp                              &                 ! =>1 where up-stream is not needed
                  &         - MAX ( rnfmsk(:,:) * rnfmsk_z(jk),  &                 ! =>0 near runoff mouths (& closed sea outflows)
                  &                 upsmsk(:,:)                ) * tmask(:,:,jk)   ! =>0 in some user defined area
            END DO
         ENDIF 
         !
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
      DO jn = 1, kjpt            !==  loop over the tracers  ==!
         !
         !                          !* Horizontal advective fluxes
         !
         !                                !-- first guess of the slopes
         zwx(:,:,jpk) = 0._wp                   ! bottom values
         zwy(:,:,jpk) = 0._wp  
         DO jk = 1, jpkm1                       ! interior values
            DO jj = 1, jpjm1      
               DO ji = 1, fs_jpim1   ! vector opt.
                  zwx(ji,jj,jk) = umask(ji,jj,jk) * ( ptb(ji+1,jj,jk,jn) - ptb(ji,jj,jk,jn) )
                  zwy(ji,jj,jk) = vmask(ji,jj,jk) * ( ptb(ji,jj+1,jk,jn) - ptb(ji,jj,jk,jn) )
               END DO
           END DO
         END DO
         ! lateral boundary conditions   (changed sign)
         CALL lbc_lnk_multi( 'traadv_mus', zwx, 'U', -1. , zwy, 'V', -1. )
         !                                !-- Slopes of tracer
         zslpx(:,:,jpk) = 0._wp                 ! bottom values
         zslpy(:,:,jpk) = 0._wp
         DO jk = 1, jpkm1                       ! interior values
            DO jj = 2, jpj
               DO ji = fs_2, jpi   ! vector opt.
                  zslpx(ji,jj,jk) =                    ( zwx(ji,jj,jk) + zwx(ji-1,jj  ,jk) )   &
                     &            * ( 0.25 + SIGN( 0.25, zwx(ji,jj,jk) * zwx(ji-1,jj  ,jk) ) )
                  zslpy(ji,jj,jk) =                    ( zwy(ji,jj,jk) + zwy(ji  ,jj-1,jk) )   &
                     &            * ( 0.25 + SIGN( 0.25, zwy(ji,jj,jk) * zwy(ji  ,jj-1,jk) ) )
               END DO
            END DO
         END DO
         !
         DO jk = 1, jpkm1                 !-- Slopes limitation
            DO jj = 2, jpj
               DO ji = fs_2, jpi   ! vector opt.
                  zslpx(ji,jj,jk) = SIGN( 1., zslpx(ji,jj,jk) ) * MIN(    ABS( zslpx(ji  ,jj,jk) ),   &
                     &                                                 2.*ABS( zwx  (ji-1,jj,jk) ),   &
                     &                                                 2.*ABS( zwx  (ji  ,jj,jk) ) )
                  zslpy(ji,jj,jk) = SIGN( 1., zslpy(ji,jj,jk) ) * MIN(    ABS( zslpy(ji,jj  ,jk) ),   &
                     &                                                 2.*ABS( zwy  (ji,jj-1,jk) ),   &
                     &                                                 2.*ABS( zwy  (ji,jj  ,jk) ) )
               END DO
           END DO
         END DO
         !
         DO jk = 1, jpkm1                 !-- MUSCL horizontal advective fluxes
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! MUSCL fluxes
                  z0u = SIGN( 0.5, pun(ji,jj,jk) )
                  zalpha = 0.5 - z0u
                  zu  = z0u - 0.5 * pun(ji,jj,jk) * p2dt * r1_e1e2u(ji,jj) / e3u_n(ji,jj,jk)
                  zzwx = ptb(ji+1,jj,jk,jn) + xind(ji,jj,jk) * zu * zslpx(ji+1,jj,jk)
                  zzwy = ptb(ji  ,jj,jk,jn) + xind(ji,jj,jk) * zu * zslpx(ji  ,jj,jk)
                  zwx(ji,jj,jk) = pun(ji,jj,jk) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
                  !
                  z0v = SIGN( 0.5, pvn(ji,jj,jk) )
                  zalpha = 0.5 - z0v
                  zv  = z0v - 0.5 * pvn(ji,jj,jk) * p2dt * r1_e1e2v(ji,jj) / e3v_n(ji,jj,jk)
                  zzwx = ptb(ji,jj+1,jk,jn) + xind(ji,jj,jk) * zv * zslpy(ji,jj+1,jk)
                  zzwy = ptb(ji,jj  ,jk,jn) + xind(ji,jj,jk) * zv * zslpy(ji,jj  ,jk)
                  zwy(ji,jj,jk) = pvn(ji,jj,jk) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
               END DO
            END DO
         END DO
         CALL lbc_lnk_multi( 'traadv_mus', zwx, 'U', -1. , zwy, 'V', -1. )   ! lateral boundary conditions   (changed sign)
         !
         DO jk = 1, jpkm1                 !-- Tracer advective trend
            DO jj = 2, jpjm1      
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) - ( zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )       &
                  &                                     + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  ) )     &
                  &                                   * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
               END DO
           END DO
         END DO        
         !                                ! trend diagnostics
         IF( l_trd )  THEN
            CALL trd_tra( kt, cdtype, jn, jptra_xad, zwx, pun, ptb(:,:,:,jn) )
            CALL trd_tra( kt, cdtype, jn, jptra_yad, zwy, pvn, ptb(:,:,:,jn) )
         END IF
         !                                 ! "Poleward" heat and salt transports 
         IF( l_ptr )  CALL dia_ptr_hst( jn, 'adv', zwy(:,:,:) )
         !                                 !  heat transport
         IF( l_hst )  CALL dia_ar5_hst( jn, 'adv', zwx(:,:,:), zwy(:,:,:) )
         !
         !                          !* Vertical advective fluxes
         !
         !                                !-- first guess of the slopes
         zwx(:,:, 1 ) = 0._wp                   ! surface & bottom boundary conditions
         zwx(:,:,jpk) = 0._wp
         DO jk = 2, jpkm1                       ! interior values
            zwx(:,:,jk) = tmask(:,:,jk) * ( ptb(:,:,jk-1,jn) - ptb(:,:,jk,jn) )
         END DO
         !                                !-- Slopes of tracer
         zslpx(:,:,1) = 0._wp                   ! surface values
         DO jk = 2, jpkm1                       ! interior value
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zslpx(ji,jj,jk) =                     ( zwx(ji,jj,jk) + zwx(ji,jj,jk+1) )  &
                     &            * (  0.25 + SIGN( 0.25, zwx(ji,jj,jk) * zwx(ji,jj,jk+1) )  )
               END DO
            END DO
         END DO
         DO jk = 2, jpkm1                 !-- Slopes limitation
            DO jj = 1, jpj                      ! interior values
               DO ji = 1, jpi
                  zslpx(ji,jj,jk) = SIGN( 1., zslpx(ji,jj,jk) ) * MIN(    ABS( zslpx(ji,jj,jk  ) ),   &
                     &                                                 2.*ABS( zwx  (ji,jj,jk+1) ),   &
                     &                                                 2.*ABS( zwx  (ji,jj,jk  ) )  )
               END DO
            END DO
         END DO
         DO jk = 1, jpk-2                 !-- vertical advective flux
            DO jj = 2, jpjm1      
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  z0w = SIGN( 0.5, pwn(ji,jj,jk+1) )
                  zalpha = 0.5 + z0w
                  zw  = z0w - 0.5 * pwn(ji,jj,jk+1) * p2dt * r1_e1e2t(ji,jj) / e3w_n(ji,jj,jk+1)
                  zzwx = ptb(ji,jj,jk+1,jn) + xind(ji,jj,jk) * zw * zslpx(ji,jj,jk+1)
                  zzwy = ptb(ji,jj,jk  ,jn) + xind(ji,jj,jk) * zw * zslpx(ji,jj,jk  )
                  zwx(ji,jj,jk+1) = pwn(ji,jj,jk+1) * ( zalpha * zzwx + (1.-zalpha) * zzwy ) * wmask(ji,jj,jk)
               END DO 
            END DO
         END DO
         IF( ln_linssh ) THEN                   ! top values, linear free surface only
            IF( ln_isfcav ) THEN                      ! ice-shelf cavities (top of the ocean)
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zwx(ji,jj, mikt(ji,jj) ) = pwn(ji,jj,mikt(ji,jj)) * ptb(ji,jj,mikt(ji,jj),jn)
                  END DO
               END DO   
            ELSE                                      ! no cavities: only at the ocean surface
               zwx(:,:,1) = pwn(:,:,1) * ptb(:,:,1,jn)
            ENDIF
         ENDIF
         !
         DO jk = 1, jpkm1                 !-- vertical advective trend
            DO jj = 2, jpjm1      
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  pta(ji,jj,jk,jn) =  pta(ji,jj,jk,jn) - ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) ) * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
               END DO
            END DO
         END DO
         !                                ! send trends for diagnostic
         IF( l_trd )  CALL trd_tra( kt, cdtype, jn, jptra_zad, zwx, pwn, ptb(:,:,:,jn) )
         !
      END DO                     ! end of tracer loop
      !
   END SUBROUTINE tra_adv_mus

   !!======================================================================
END MODULE traadv_mus
