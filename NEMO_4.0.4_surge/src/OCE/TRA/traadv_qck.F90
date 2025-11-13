MODULE traadv_qck
   !!==============================================================================
   !!                       ***  MODULE  traadv_qck  ***
   !! Ocean tracers:  horizontal & vertical advective trend
   !!==============================================================================
   !! History :  3.0  !  2008-07  (G. Reffray)  Original code
   !!            3.3  !  2010-05  (C.Ethe, G. Madec)  merge TRC-TRA + switch from velocity to transport
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv_qck    : update the tracer trend with the horizontal advection
   !!                    trends using a 3rd order finite difference scheme
   !!   tra_adv_qck_i  : apply QUICK scheme in i-direction
   !!   tra_adv_qck_j  : apply QUICK scheme in j-direction
   !!   tra_adv_cen2_k : 2nd centered scheme for the vertical advection
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE trc_oce         ! share passive tracers/Ocean variables
   USE trd_oce         ! trends: ocean variables
   USE trdtra          ! trends manager: tracers 
   USE diaptr          ! poleward transport diagnostics
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing
   USE lbclnk          ! ocean lateral boundary condition (or mpp link)
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_qck   ! routine called by step.F90

   REAL(wp) :: r1_6 = 1./ 6.   ! 1/6 ratio

   LOGICAL  ::   l_trd   ! flag to compute trends
   LOGICAL  ::   l_ptr   ! flag to compute poleward transport


   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv_qck ( kt, kit000, cdtype, p2dt, pun, pvn, pwn,      &
      &                                       ptb, ptn, pta, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_qck  ***
      !!
      !! ** Purpose :   Compute the now trend due to the advection of tracers
      !!      and add it to the general trend of passive tracer equations.
      !!
      !! ** Method :   The advection is evaluated by a third order scheme
      !!             For a positive velocity u :              u(i)>0
      !!                                          |--FU--|--FC--|--FD--|------|
      !!                                             i-1    i      i+1   i+2
      !!
      !!             For a negative velocity u :              u(i)<0
      !!                                          |------|--FD--|--FC--|--FU--|
      !!                                             i-1    i      i+1   i+2
      !!             where  FU is the second upwind point
      !!                    FD is the first douwning point
      !!                    FC is the central point (or the first upwind point)
      !!
      !!      Flux(i) = u(i) * { 0.5(FC+FD)  -0.5C(i)(FD-FC)  -((1-C(i))/6)(FU+FD-2FC) }
      !!                with C(i)=|u(i)|dx(i)/dt (=Courant number)
      !!
      !!         dt = 2*rdtra and the scalar values are tb and sb
      !!
      !!       On the vertical, the simple centered scheme used ptn
      !!
      !!               The fluxes are bounded by the ULTIMATE limiter to
      !!             guarantee the monotonicity of the solution and to
      !!            prevent the appearance of spurious numerical oscillations
      !!
      !! ** Action : - update pta  with the now advective tracer trends
      !!             - send trends to trdtra module for further diagnostcs (l_trdtra=T)
      !!             - htr_adv, str_adv : poleward advective heat and salt transport (ln_diaptr=T)
      !!
      !! ** Reference : Leonard (1979, 1991)
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt              ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt            ! number of tracers
      REAL(wp)                             , INTENT(in   ) ::   p2dt            ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun, pvn, pwn   ! 3 ocean velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb, ptn        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta             ! tracer trend 
      !!----------------------------------------------------------------------
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_qck : 3rd order quickest advection scheme on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
         IF(lwp) WRITE(numout,*)
      ENDIF
      !
      l_trd = .FALSE.
      l_ptr = .FALSE.
      IF( ( cdtype == 'TRA' .AND. l_trdtra ) .OR. ( cdtype == 'TRC' .AND. l_trdtrc ) )      l_trd = .TRUE.
      IF(   cdtype == 'TRA' .AND. ln_diaptr )                                               l_ptr = .TRUE. 
      !
      !
      !        ! horizontal fluxes are computed with the QUICKEST + ULTIMATE scheme
      CALL tra_adv_qck_i( kt, cdtype, p2dt, pun, ptb, ptn, pta, kjpt ) 
      CALL tra_adv_qck_j( kt, cdtype, p2dt, pvn, ptb, ptn, pta, kjpt ) 

      !        ! vertical fluxes are computed with the 2nd order centered scheme
      CALL tra_adv_cen2_k( kt, cdtype, pwn,         ptn, pta, kjpt )
      !
   END SUBROUTINE tra_adv_qck


   SUBROUTINE tra_adv_qck_i( kt, cdtype, p2dt, pun,                  &
      &                                        ptb, ptn, pta, kjpt   )
      !!----------------------------------------------------------------------
      !!
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      REAL(wp)                             , INTENT(in   ) ::   p2dt       ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun        ! i-velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb, ptn   ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta        ! tracer trend 
      !!
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::   ztra, zbtr, zdir, zdx, zmsk   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zwx, zfu, zfc, zfd
      !----------------------------------------------------------------------
      !
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         zfu(:,:,:) = 0._wp     ;   zfc(:,:,:) = 0._wp 
         zfd(:,:,:) = 0._wp     ;   zwx(:,:,:) = 0._wp   
         !
!!gm why not using a SHIFT instruction...
         DO jk = 1, jpkm1     !--- Computation of the ustream and downstream value of the tracer and the mask
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zfc(ji,jj,jk) = ptb(ji-1,jj,jk,jn)        ! Upstream   in the x-direction for the tracer
                  zfd(ji,jj,jk) = ptb(ji+1,jj,jk,jn)        ! Downstream in the x-direction for the tracer
               END DO
            END DO
         END DO
         CALL lbc_lnk_multi( 'traadv_qck', zfc(:,:,:), 'T', 1. , zfd(:,:,:), 'T', 1. )   ! Lateral boundary conditions 
         
         !
         ! Horizontal advective fluxes
         ! ---------------------------
         DO jk = 1, jpkm1                             
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.         
                  zdir = 0.5 + SIGN( 0.5, pun(ji,jj,jk) )   ! if pun > 0 : zdir = 1 otherwise zdir = 0 
                  zfu(ji,jj,jk) = zdir * zfc(ji,jj,jk ) + ( 1. - zdir ) * zfd(ji+1,jj,jk)  ! FU in the x-direction for T 
               END DO
            END DO
         END DO
         !
         DO jk = 1, jpkm1  
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.   
                  zdir = 0.5 + SIGN( 0.5, pun(ji,jj,jk) )   ! if pun > 0 : zdir = 1 otherwise zdir = 0 
                  zdx = ( zdir * e1t(ji,jj) + ( 1. - zdir ) * e1t(ji+1,jj) ) * e2u(ji,jj) * e3u_n(ji,jj,jk)
                  zwx(ji,jj,jk)  = ABS( pun(ji,jj,jk) ) * p2dt / zdx    ! (0<zc_cfl<1 : Courant number on x-direction)
                  zfc(ji,jj,jk)  = zdir * ptb(ji  ,jj,jk,jn) + ( 1. - zdir ) * ptb(ji+1,jj,jk,jn)  ! FC in the x-direction for T
                  zfd(ji,jj,jk)  = zdir * ptb(ji+1,jj,jk,jn) + ( 1. - zdir ) * ptb(ji  ,jj,jk,jn)  ! FD in the x-direction for T
               END DO
            END DO
         END DO 
         !--- Lateral boundary conditions 
         CALL lbc_lnk_multi( 'traadv_qck', zfu(:,:,:), 'T', 1. , zfd(:,:,:), 'T', 1., zfc(:,:,:), 'T', 1.,  zwx(:,:,:), 'T', 1. )

         !--- QUICKEST scheme
         CALL quickest( zfu, zfd, zfc, zwx )
         !
         ! Mask at the T-points in the x-direction (mask=0 or mask=1)
         DO jk = 1, jpkm1  
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.               
                  zfu(ji,jj,jk) = tmask(ji-1,jj,jk) + tmask(ji,jj,jk) + tmask(ji+1,jj,jk) - 2.
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'traadv_qck', zfu(:,:,:), 'T', 1. )      ! Lateral boundary conditions 

         !
         ! Tracer flux on the x-direction
         DO jk = 1, jpkm1  
            !
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.               
                  zdir = 0.5 + SIGN( 0.5, pun(ji,jj,jk) )   ! if pun > 0 : zdir = 1 otherwise zdir = 0 
                  !--- If the second ustream point is a land point
                  !--- the flux is computed by the 1st order UPWIND scheme
                  zmsk = zdir * zfu(ji,jj,jk) + ( 1. - zdir ) * zfu(ji+1,jj,jk)
                  zwx(ji,jj,jk) = zmsk * zwx(ji,jj,jk) + ( 1. - zmsk ) * zfc(ji,jj,jk)
                  zwx(ji,jj,jk) = zwx(ji,jj,jk) * pun(ji,jj,jk)
               END DO
            END DO
         END DO
         !
         CALL lbc_lnk( 'traadv_qck', zwx(:,:,:), 'T', 1. ) ! Lateral boundary conditions
         !
         ! Computation of the trend
         DO jk = 1, jpkm1  
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.  
                  zbtr = r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
                  ! horizontal advective trends
                  ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji-1,jj,jk) )
                  !--- add it to the general tracer trends
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + ztra
               END DO
            END DO
         END DO
         !                                 ! trend diagnostics
         IF( l_trd )   CALL trd_tra( kt, cdtype, jn, jptra_xad, zwx, pun, ptn(:,:,:,jn) )
         !
      END DO
      !
   END SUBROUTINE tra_adv_qck_i


   SUBROUTINE tra_adv_qck_j( kt, cdtype, p2dt, pvn,                &
      &                                        ptb, ptn, pta, kjpt )
      !!----------------------------------------------------------------------
      !!
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      REAL(wp)                             , INTENT(in   ) ::   p2dt       ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pvn        ! j-velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb, ptn   ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta        ! tracer trend 
      !!
      INTEGER  :: ji, jj, jk, jn                ! dummy loop indices
      REAL(wp) :: ztra, zbtr, zdir, zdx, zmsk   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zwy, zfu, zfc, zfd   ! 3D workspace
      !----------------------------------------------------------------------
      !
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         zfu(:,:,:) = 0.0     ;   zfc(:,:,:) = 0.0  
         zfd(:,:,:) = 0.0     ;   zwy(:,:,:) = 0.0     
         !                                                  
         DO jk = 1, jpkm1                                
            !                                             
            !--- Computation of the ustream and downstream value of the tracer and the mask
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! Upstream in the x-direction for the tracer
                  zfc(ji,jj,jk) = ptb(ji,jj-1,jk,jn)
                  ! Downstream in the x-direction for the tracer
                  zfd(ji,jj,jk) = ptb(ji,jj+1,jk,jn)
               END DO
            END DO
         END DO
         CALL lbc_lnk_multi( 'traadv_qck', zfc(:,:,:), 'T', 1. , zfd(:,:,:), 'T', 1. )   ! Lateral boundary conditions 

         
         !
         ! Horizontal advective fluxes
         ! ---------------------------
         !
         DO jk = 1, jpkm1                             
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.         
                  zdir = 0.5 + SIGN( 0.5, pvn(ji,jj,jk) )   ! if pun > 0 : zdir = 1 otherwise zdir = 0 
                  zfu(ji,jj,jk) = zdir * zfc(ji,jj,jk ) + ( 1. - zdir ) * zfd(ji,jj+1,jk)  ! FU in the x-direction for T 
               END DO
            END DO
         END DO
         !
         DO jk = 1, jpkm1  
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.   
                  zdir = 0.5 + SIGN( 0.5, pvn(ji,jj,jk) )   ! if pun > 0 : zdir = 1 otherwise zdir = 0 
                  zdx = ( zdir * e2t(ji,jj) + ( 1. - zdir ) * e2t(ji,jj+1) ) * e1v(ji,jj) * e3v_n(ji,jj,jk)
                  zwy(ji,jj,jk)  = ABS( pvn(ji,jj,jk) ) * p2dt / zdx    ! (0<zc_cfl<1 : Courant number on x-direction)
                  zfc(ji,jj,jk)  = zdir * ptb(ji,jj  ,jk,jn) + ( 1. - zdir ) * ptb(ji,jj+1,jk,jn)  ! FC in the x-direction for T
                  zfd(ji,jj,jk)  = zdir * ptb(ji,jj+1,jk,jn) + ( 1. - zdir ) * ptb(ji,jj  ,jk,jn)  ! FD in the x-direction for T
               END DO
            END DO
         END DO

         !--- Lateral boundary conditions 
         CALL lbc_lnk_multi( 'traadv_qck', zfu(:,:,:), 'T', 1. , zfd(:,:,:), 'T', 1., zfc(:,:,:), 'T', 1., zwy(:,:,:), 'T', 1. )

         !--- QUICKEST scheme
         CALL quickest( zfu, zfd, zfc, zwy )
         !
         ! Mask at the T-points in the x-direction (mask=0 or mask=1)
         DO jk = 1, jpkm1  
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.               
                  zfu(ji,jj,jk) = tmask(ji,jj-1,jk) + tmask(ji,jj,jk) + tmask(ji,jj+1,jk) - 2.
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'traadv_qck', zfu(:,:,:), 'T', 1. )    !--- Lateral boundary conditions 
         !
         ! Tracer flux on the x-direction
         DO jk = 1, jpkm1  
            !
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.               
                  zdir = 0.5 + SIGN( 0.5, pvn(ji,jj,jk) )   ! if pun > 0 : zdir = 1 otherwise zdir = 0 
                  !--- If the second ustream point is a land point
                  !--- the flux is computed by the 1st order UPWIND scheme
                  zmsk = zdir * zfu(ji,jj,jk) + ( 1. - zdir ) * zfu(ji,jj+1,jk)
                  zwy(ji,jj,jk) = zmsk * zwy(ji,jj,jk) + ( 1. - zmsk ) * zfc(ji,jj,jk)
                  zwy(ji,jj,jk) = zwy(ji,jj,jk) * pvn(ji,jj,jk)
               END DO
            END DO
         END DO
         !
         CALL lbc_lnk( 'traadv_qck', zwy(:,:,:), 'T', 1. ) ! Lateral boundary conditions
         !
         ! Computation of the trend
         DO jk = 1, jpkm1  
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.  
                  zbtr = r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
                  ! horizontal advective trends
                  ztra = - zbtr * ( zwy(ji,jj,jk) - zwy(ji,jj-1,jk) )
                  !--- add it to the general tracer trends
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + ztra
               END DO
            END DO
         END DO
         !                                 ! trend diagnostics
         IF( l_trd )   CALL trd_tra( kt, cdtype, jn, jptra_yad, zwy, pvn, ptn(:,:,:,jn) )
         !                                 ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( l_ptr )   CALL dia_ptr_hst( jn, 'adv', zwy(:,:,:) )
         !
      END DO
      !
   END SUBROUTINE tra_adv_qck_j


   SUBROUTINE tra_adv_cen2_k( kt, cdtype, pwn,           &
     &                                    ptn, pta, kjpt )
      !!----------------------------------------------------------------------
      !!
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt       ! ocean time-step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt     ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pwn      ! vertical velocity 
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptn      ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta      ! tracer trend 
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zwz   ! 3D workspace
      !!----------------------------------------------------------------------
      !
      zwz(:,:, 1 ) = 0._wp       ! surface & bottom values set to zero for all tracers
      zwz(:,:,jpk) = 0._wp
      !
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         !
         DO jk = 2, jpkm1                    !* Interior point   (w-masked 2nd order centered flux)
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zwz(ji,jj,jk) = 0.5 * pwn(ji,jj,jk) * ( ptn(ji,jj,jk-1,jn) + ptn(ji,jj,jk,jn) ) * wmask(ji,jj,jk)
               END DO
            END DO
         END DO
         IF( ln_linssh ) THEN                !* top value   (only in linear free surf. as zwz is multiplied by wmask)
            IF( ln_isfcav ) THEN                  ! ice-shelf cavities (top of the ocean)
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zwz(ji,jj, mikt(ji,jj) ) = pwn(ji,jj,mikt(ji,jj)) * ptn(ji,jj,mikt(ji,jj),jn)   ! linear free surface 
                  END DO
               END DO   
            ELSE                                   ! no ocean cavities (only ocean surface)
               zwz(:,:,1) = pwn(:,:,1) * ptn(:,:,1,jn)
            ENDIF
         ENDIF
         !
         DO jk = 1, jpkm1          !==  Tracer flux divergence added to the general trend  ==!
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) - ( zwz(ji,jj,jk) - zwz(ji,jj,jk+1) )   &
                     &                                * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
               END DO
            END DO
         END DO
         !                                 ! Send trends for diagnostic
         IF( l_trd )  CALL trd_tra( kt, cdtype, jn, jptra_zad, zwz, pwn, ptn(:,:,:,jn) )
         !
      END DO
      !
   END SUBROUTINE tra_adv_cen2_k


   SUBROUTINE quickest( pfu, pfd, pfc, puc )
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose :  Computation of advective flux with Quickest scheme
      !!
      !! ** Method :   
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pfu   ! second upwind point
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pfd   ! first douwning point
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pfc   ! the central point (or the first upwind point)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   puc   ! input as Courant number ; output as flux
      !!
      INTEGER  ::  ji, jj, jk               ! dummy loop indices 
      REAL(wp) ::  zcoef1, zcoef2, zcoef3   ! local scalars          
      REAL(wp) ::  zc, zcurv, zfho          !   -      -
      !----------------------------------------------------------------------
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zc     = puc(ji,jj,jk)                         ! Courant number
               zcurv  = pfd(ji,jj,jk) + pfu(ji,jj,jk) - 2. * pfc(ji,jj,jk)
               zcoef1 = 0.5 *      ( pfc(ji,jj,jk) + pfd(ji,jj,jk) )
               zcoef2 = 0.5 * zc * ( pfd(ji,jj,jk) - pfc(ji,jj,jk) )
               zcoef3 = ( 1. - ( zc * zc ) ) * r1_6 * zcurv
               zfho   = zcoef1 - zcoef2 - zcoef3              !  phi_f QUICKEST 
               !
               zcoef1 = pfd(ji,jj,jk) - pfu(ji,jj,jk)
               zcoef2 = ABS( zcoef1 )
               zcoef3 = ABS( zcurv )
               IF( zcoef3 >= zcoef2 ) THEN
                  zfho = pfc(ji,jj,jk) 
               ELSE
                  zcoef3 = pfu(ji,jj,jk) + ( ( pfc(ji,jj,jk) - pfu(ji,jj,jk) ) / MAX( zc, 1.e-9 ) )    ! phi_REF
                  IF( zcoef1 >= 0. ) THEN
                     zfho = MAX( pfc(ji,jj,jk), zfho ) 
                     zfho = MIN( zfho, MIN( zcoef3, pfd(ji,jj,jk) ) ) 
                  ELSE
                     zfho = MIN( pfc(ji,jj,jk), zfho ) 
                     zfho = MAX( zfho, MAX( zcoef3, pfd(ji,jj,jk) ) ) 
                  ENDIF
               ENDIF
               puc(ji,jj,jk) = zfho
            END DO
         END DO
      END DO
      !
   END SUBROUTINE quickest

   !!======================================================================
END MODULE traadv_qck
