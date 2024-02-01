MODULE traldf_triad
   !!======================================================================
   !!                   ***  MODULE  traldf_triad  ***
   !! Ocean  tracers:  horizontal component of the lateral tracer mixing trend
   !!======================================================================
   !! History :  3.3  ! 2010-10  (G. Nurser, C. Harris, G. Madec)  Griffies operator (original code)
   !!            3.7  ! 2013-12  (F. Lemarie, G. Madec)  triad operator (Griffies) + Method of Stabilizing Correction
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_ldf_triad : update the tracer trend with the iso-neutral laplacian triad-operator
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE trc_oce        ! share passive tracers/Ocean variables
   USE zdf_oce        ! ocean vertical physics
   USE ldftra         ! lateral physics: eddy diffusivity
   USE ldfslp         ! lateral physics: iso-neutral slopes
   USE traldf_iso     ! lateral diffusion (Madec operator)         (tra_ldf_iso routine)
   USE diaptr         ! poleward transport diagnostics
   USE diaar5         ! AR5 diagnostics
   USE zpshde         ! partial step: hor. derivative     (zps_hde routine)
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf_triad   ! routine called by traldf.F90

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, SAVE ::   zdkt3d   !: vertical tracer gradient at 2 levels

   LOGICAL  ::   l_ptr   ! flag to compute poleward transport
   LOGICAL  ::   l_hst   ! flag to compute heat transport


   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

  SUBROUTINE tra_ldf_triad( kt, kit000, cdtype, pahu, pahv, pgu , pgv ,   &
      &                                                     pgui, pgvi,   &
      &                                         ptb , ptbb, pta , kjpt, kpass )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_triad  ***
      !!
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive
      !!      trend for a laplacian tensor (ezxcept the dz[ dz[.] ] term) and
      !!      add it to the general trend of tracer equation.
      !!
      !! ** Method  :   The horizontal component of the lateral diffusive trends
      !!      is provided by a 2nd order operator rotated along neural or geopo-
      !!      tential surfaces to which an eddy induced advection can be added
      !!      It is computed using before fields (forward in time) and isopyc-
      !!      nal or geopotential slopes computed in routine ldfslp.
      !!
      !!      see documentation for the desciption
      !!
      !! ** Action :   pta   updated with the before rotated diffusion
      !!               ah_wslp2 ....
      !!               akz   stabilizing vertical diffusivity coefficient (used in trazdf_imp)
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000     ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      INTEGER                              , INTENT(in   ) ::   kpass      ! =1/2 first or second passage
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pahu, pahv ! eddy diffusivity at u- and v-points  [m2/s]
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(in   ) ::   pgu , pgv  ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,    kjpt), INTENT(in   ) ::   pgui, pgvi ! tracer gradient at top   levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb        ! tracer (kpass=1) or laplacian of tracer (kpass=2)
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptbb       ! tracer (only used in kpass=2)
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta        ! tracer trend
      !
      INTEGER  ::  ji, jj, jk, jn   ! dummy loop indices
      INTEGER  ::  ip,jp,kp         ! dummy loop indices
      INTEGER  ::  ierr            ! local integer
      REAL(wp) ::  zmsku, zabe1, zcof1, zcoef3          ! local scalars
      REAL(wp) ::  zmskv, zabe2, zcof2, zcoef4          !   -      -
      REAL(wp) ::  zcoef0, ze3w_2, zsign, z2dt, z1_2dt  !   -      -
      !
      REAL(wp) ::   zslope_skew, zslope_iso, zslope2, zbu, zbv
      REAL(wp) ::   ze1ur, ze2vr, ze3wr, zdxt, zdyt, zdzt
      REAL(wp) ::   zah, zah_slp, zaei_slp
      REAL(wp), DIMENSION(jpi,jpj    ) ::   z2d                                              ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zdit, zdjt, zftu, zftv, ztfw, zpsi_uw, zpsi_vw   ! 3D     -
      !!----------------------------------------------------------------------
      !
      IF( .NOT.ALLOCATED(zdkt3d) )  THEN
         ALLOCATE( zdkt3d(jpi,jpj,0:1) , STAT=ierr )
         CALL mpp_sum ( 'traldf_triad', ierr )
         IF( ierr > 0 )   CALL ctl_stop('STOP', 'tra_ldf_triad: unable to allocate arrays')
      ENDIF
     !
      IF( kpass == 1 .AND. kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_triad : rotated laplacian diffusion operator on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !   
      l_hst = .FALSE.
      l_ptr = .FALSE.
      IF( cdtype == 'TRA' .AND. ln_diaptr )                                                 l_ptr = .TRUE. 
      IF( cdtype == 'TRA' .AND. ( iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. &
         &                        iom_use("uadv_salttr") .OR. iom_use("vadv_salttr")  ) )   l_hst = .TRUE.
      !
      !                                                        ! set time step size (Euler/Leapfrog)
      IF( neuler == 0 .AND. kt == kit000 ) THEN   ;   z2dt =     rdt      ! at nit000   (Euler)
      ELSE                                        ;   z2dt = 2.* rdt      !             (Leapfrog)
      ENDIF
      z1_2dt = 1._wp / z2dt
      !
      IF( kpass == 1 ) THEN   ;   zsign =  1._wp      ! bilaplacian operator require a minus sign (eddy diffusivity >0)
      ELSE                    ;   zsign = -1._wp
      ENDIF
      !    
      !!----------------------------------------------------------------------
      !!   0 - calculate  ah_wslp2, akz, and optionally zpsi_uw, zpsi_vw
      !!----------------------------------------------------------------------
      !
      IF( kpass == 1 ) THEN         !==  first pass only  and whatever the tracer is  ==!
         !
         akz     (:,:,:) = 0._wp      
         ah_wslp2(:,:,:) = 0._wp
         IF( ln_ldfeiv_dia ) THEN
            zpsi_uw(:,:,:) = 0._wp
            zpsi_vw(:,:,:) = 0._wp
         ENDIF
         !
         DO ip = 0, 1                            ! i-k triads
            DO kp = 0, 1
               DO jk = 1, jpkm1
                  DO jj = 1, jpjm1
                     DO ji = 1, fs_jpim1
                        ze3wr = 1._wp / e3w_n(ji+ip,jj,jk+kp)
                        zbu   = e1e2u(ji,jj) * e3u_n(ji,jj,jk)
                        zah   = 0.25_wp * pahu(ji,jj,jk)
                        zslope_skew = triadi_g(ji+ip,jj,jk,1-ip,kp)
                        ! Subtract s-coordinate slope at t-points to give slope rel to s-surfaces (do this by *adding* gradient of depth)
                        zslope2 = zslope_skew + ( gdept_n(ji+1,jj,jk) - gdept_n(ji,jj,jk) ) * r1_e1u(ji,jj) * umask(ji,jj,jk+kp)
                        zslope2 = zslope2 *zslope2
                        ah_wslp2(ji+ip,jj,jk+kp) = ah_wslp2(ji+ip,jj,jk+kp) + zah * zbu * ze3wr * r1_e1e2t(ji+ip,jj) * zslope2
                        akz     (ji+ip,jj,jk+kp) = akz     (ji+ip,jj,jk+kp) + zah * r1_e1u(ji,jj)       &
                           &                                                      * r1_e1u(ji,jj) * umask(ji,jj,jk+kp)
                           !
                       IF( ln_ldfeiv_dia )   zpsi_uw(ji,jj,jk+kp) = zpsi_uw(ji,jj,jk+kp)   &
                           &                                       + 0.25_wp * aeiu(ji,jj,jk) * e2u(ji,jj) * zslope_skew
                     END DO
                  END DO
               END DO
            END DO
         END DO
         !
         DO jp = 0, 1                            ! j-k triads 
            DO kp = 0, 1
               DO jk = 1, jpkm1
                  DO jj = 1, jpjm1
                     DO ji = 1, fs_jpim1
                        ze3wr = 1.0_wp / e3w_n(ji,jj+jp,jk+kp)
                        zbv   = e1e2v(ji,jj) * e3v_n(ji,jj,jk)
                        zah   = 0.25_wp * pahv(ji,jj,jk)
                        zslope_skew = triadj_g(ji,jj+jp,jk,1-jp,kp)
                        ! Subtract s-coordinate slope at t-points to give slope rel to s surfaces
                        !    (do this by *adding* gradient of depth)
                        zslope2 = zslope_skew + ( gdept_n(ji,jj+1,jk) - gdept_n(ji,jj,jk) ) * r1_e2v(ji,jj) * vmask(ji,jj,jk+kp)
                        zslope2 = zslope2 * zslope2
                        ah_wslp2(ji,jj+jp,jk+kp) = ah_wslp2(ji,jj+jp,jk+kp) + zah * zbv * ze3wr * r1_e1e2t(ji,jj+jp) * zslope2
                        akz     (ji,jj+jp,jk+kp) = akz     (ji,jj+jp,jk+kp) + zah * r1_e2v(ji,jj)     &
                           &                                                      * r1_e2v(ji,jj) * vmask(ji,jj,jk+kp)
                        !
                        IF( ln_ldfeiv_dia )   zpsi_vw(ji,jj,jk+kp) = zpsi_vw(ji,jj,jk+kp)   &
                           &                                       + 0.25 * aeiv(ji,jj,jk) * e1v(ji,jj) * zslope_skew
                     END DO
                  END DO
               END DO
            END DO
         END DO
         !
         IF( ln_traldf_msc ) THEN                ! stabilizing vertical diffusivity coefficient
            !
            IF( ln_traldf_blp ) THEN                ! bilaplacian operator
               DO jk = 2, jpkm1
                  DO jj = 1, jpjm1
                     DO ji = 1, fs_jpim1
                        akz(ji,jj,jk) = 16._wp * ah_wslp2(ji,jj,jk)   &
                           &          * (  akz(ji,jj,jk) + ah_wslp2(ji,jj,jk) / ( e3w_n(ji,jj,jk) * e3w_n(ji,jj,jk) )  )
                     END DO
                  END DO
               END DO
            ELSEIF( ln_traldf_lap ) THEN              ! laplacian operator
               DO jk = 2, jpkm1
                  DO jj = 1, jpjm1
                     DO ji = 1, fs_jpim1
                        ze3w_2 = e3w_n(ji,jj,jk) * e3w_n(ji,jj,jk)
                        zcoef0 = z2dt * (  akz(ji,jj,jk) + ah_wslp2(ji,jj,jk) / ze3w_2  )
                        akz(ji,jj,jk) = MAX( zcoef0 - 0.5_wp , 0._wp ) * ze3w_2 * z1_2dt
                     END DO
                  END DO
               END DO
           ENDIF
           !
         ELSE                                    ! 33 flux set to zero with akz=ah_wslp2 ==>> computed in full implicit
            akz(:,:,:) = ah_wslp2(:,:,:)      
         ENDIF
         !
         IF( ln_ldfeiv_dia .AND. cdtype == 'TRA' )   CALL ldf_eiv_dia( zpsi_uw, zpsi_vw )
         !
      ENDIF                                  !==  end 1st pass only  ==!
      !
      !                                                           ! ===========
      DO jn = 1, kjpt                                             ! tracer loop
         !                                                        ! ===========
         ! Zero fluxes for each tracer
!!gm  this should probably be done outside the jn loop
         ztfw(:,:,:) = 0._wp
         zftu(:,:,:) = 0._wp
         zftv(:,:,:) = 0._wp
         !
         DO jk = 1, jpkm1        !==  before lateral T & S gradients at T-level jk  ==!
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zdit(ji,jj,jk) = ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) ) * umask(ji,jj,jk)
                  zdjt(ji,jj,jk) = ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
         IF( ln_zps .AND. l_grad_zps ) THEN    ! partial steps: correction at top/bottom ocean level
            DO jj = 1, jpjm1                       ! bottom level
               DO ji = 1, fs_jpim1   ! vector opt.
                  zdit(ji,jj,mbku(ji,jj)) = pgu(ji,jj,jn)
                  zdjt(ji,jj,mbkv(ji,jj)) = pgv(ji,jj,jn)
               END DO
            END DO
            IF( ln_isfcav ) THEN                   ! top level (ocean cavities only)
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.
                     IF( miku(ji,jj)  > 1 )   zdit(ji,jj,miku(ji,jj) ) = pgui(ji,jj,jn) 
                     IF( mikv(ji,jj)  > 1 )   zdjt(ji,jj,mikv(ji,jj) ) = pgvi(ji,jj,jn) 
                  END DO
               END DO
            ENDIF
         ENDIF
         !
         !!----------------------------------------------------------------------
         !!   II - horizontal trend  (full)
         !!----------------------------------------------------------------------
         !
         DO jk = 1, jpkm1
            !                    !==  Vertical tracer gradient at level jk and jk+1
            zdkt3d(:,:,1) = ( ptb(:,:,jk,jn) - ptb(:,:,jk+1,jn) ) * tmask(:,:,jk+1)
            !
            !                    ! surface boundary condition: zdkt3d(jk=0)=zdkt3d(jk=1)
            IF( jk == 1 ) THEN   ;   zdkt3d(:,:,0) = zdkt3d(:,:,1)
            ELSE                 ;   zdkt3d(:,:,0) = ( ptb(:,:,jk-1,jn) - ptb(:,:,jk,jn) ) * tmask(:,:,jk)
            ENDIF
            !
            zaei_slp = 0._wp
            !
            IF( ln_botmix_triad ) THEN
               DO ip = 0, 1              !==  Horizontal & vertical fluxes
                  DO kp = 0, 1
                     DO jj = 1, jpjm1
                        DO ji = 1, fs_jpim1
                           ze1ur = r1_e1u(ji,jj)
                           zdxt  = zdit(ji,jj,jk) * ze1ur
                           ze3wr = 1._wp / e3w_n(ji+ip,jj,jk+kp)
                           zdzt  = zdkt3d(ji+ip,jj,kp) * ze3wr
                           zslope_skew = triadi_g(ji+ip,jj,jk,1-ip,kp)
                           zslope_iso  = triadi  (ji+ip,jj,jk,1-ip,kp)
                           !
                           zbu = 0.25_wp * e1e2u(ji,jj) * e3u_n(ji,jj,jk)
                           ! ln_botmix_triad is .T. don't mask zah for bottom half cells    !!gm ?????   ahu is masked....
                           zah = pahu(ji,jj,jk)
                           zah_slp  = zah * zslope_iso
                           IF( ln_ldfeiv )   zaei_slp = aeiu(ji,jj,jk) * zslope_skew
                           zftu(ji   ,jj,jk   ) = zftu(ji   ,jj,jk   ) - ( zah * zdxt + (zah_slp - zaei_slp) * zdzt ) * zbu * ze1ur
                           ztfw(ji+ip,jj,jk+kp) = ztfw(ji+ip,jj,jk+kp) - ( zah_slp + zaei_slp) * zdxt                 * zbu * ze3wr
                        END DO
                     END DO
                  END DO
               END DO
               !
               DO jp = 0, 1
                  DO kp = 0, 1
                     DO jj = 1, jpjm1
                        DO ji = 1, fs_jpim1
                           ze2vr = r1_e2v(ji,jj)
                           zdyt  = zdjt(ji,jj,jk) * ze2vr
                           ze3wr = 1._wp / e3w_n(ji,jj+jp,jk+kp)
                           zdzt  = zdkt3d(ji,jj+jp,kp) * ze3wr
                           zslope_skew = triadj_g(ji,jj+jp,jk,1-jp,kp)
                           zslope_iso  = triadj(ji,jj+jp,jk,1-jp,kp)
                           zbv = 0.25_wp * e1e2v(ji,jj) * e3v_n(ji,jj,jk)
                           ! ln_botmix_triad is .T. don't mask zah for bottom half cells    !!gm ?????  ahv is masked...
                           zah = pahv(ji,jj,jk)
                           zah_slp = zah * zslope_iso
                           IF( ln_ldfeiv )   zaei_slp = aeiv(ji,jj,jk) * zslope_skew
                           zftv(ji,jj   ,jk   ) = zftv(ji,jj   ,jk   ) - ( zah * zdyt + (zah_slp - zaei_slp) * zdzt ) * zbv * ze2vr
                           ztfw(ji,jj+jp,jk+kp) = ztfw(ji,jj+jp,jk+kp) - ( zah_slp + zaei_slp ) * zdyt                * zbv * ze3wr
                        END DO
                     END DO
                  END DO
               END DO
               !
            ELSE
               !
               DO ip = 0, 1               !==  Horizontal & vertical fluxes
                  DO kp = 0, 1
                     DO jj = 1, jpjm1
                        DO ji = 1, fs_jpim1
                           ze1ur = r1_e1u(ji,jj)
                           zdxt  = zdit(ji,jj,jk) * ze1ur
                           ze3wr = 1._wp / e3w_n(ji+ip,jj,jk+kp)
                           zdzt  = zdkt3d(ji+ip,jj,kp) * ze3wr
                           zslope_skew = triadi_g(ji+ip,jj,jk,1-ip,kp)
                           zslope_iso  = triadi(ji+ip,jj,jk,1-ip,kp)
                           !
                           zbu = 0.25_wp * e1e2u(ji,jj) * e3u_n(ji,jj,jk)
                           ! ln_botmix_triad is .F. mask zah for bottom half cells
                           zah = pahu(ji,jj,jk) * umask(ji,jj,jk+kp)         ! pahu(ji+ip,jj,jk)   ===>>  ????
                           zah_slp  = zah * zslope_iso
                           IF( ln_ldfeiv )   zaei_slp = aeiu(ji,jj,jk) * zslope_skew        ! aeit(ji+ip,jj,jk)*zslope_skew
                           zftu(ji   ,jj,jk   ) = zftu(ji   ,jj,jk   ) - ( zah * zdxt + (zah_slp - zaei_slp) * zdzt ) * zbu * ze1ur
                           ztfw(ji+ip,jj,jk+kp) = ztfw(ji+ip,jj,jk+kp) - (zah_slp + zaei_slp) * zdxt * zbu * ze3wr
                        END DO
                     END DO
                  END DO
               END DO
               !
               DO jp = 0, 1
                  DO kp = 0, 1
                     DO jj = 1, jpjm1
                        DO ji = 1, fs_jpim1
                           ze2vr = r1_e2v(ji,jj)
                           zdyt  = zdjt(ji,jj,jk) * ze2vr
                           ze3wr = 1._wp / e3w_n(ji,jj+jp,jk+kp)
                           zdzt  = zdkt3d(ji,jj+jp,kp) * ze3wr
                           zslope_skew = triadj_g(ji,jj+jp,jk,1-jp,kp)
                           zslope_iso  = triadj(ji,jj+jp,jk,1-jp,kp)
                           zbv = 0.25_wp * e1e2v(ji,jj) * e3v_n(ji,jj,jk)
                           ! ln_botmix_triad is .F. mask zah for bottom half cells
                           zah = pahv(ji,jj,jk) * vmask(ji,jj,jk+kp)         ! pahv(ji,jj+jp,jk)  ????
                           zah_slp = zah * zslope_iso
                           IF( ln_ldfeiv )   zaei_slp = aeiv(ji,jj,jk) * zslope_skew        ! aeit(ji,jj+jp,jk)*zslope_skew
                           zftv(ji,jj,jk) = zftv(ji,jj,jk) - ( zah * zdyt + (zah_slp - zaei_slp) * zdzt ) * zbv * ze2vr
                           ztfw(ji,jj+jp,jk+kp) = ztfw(ji,jj+jp,jk+kp) - (zah_slp + zaei_slp) * zdyt * zbv * ze3wr
                        END DO
                     END DO
                  END DO
               END DO
            ENDIF
            !                             !==  horizontal divergence and add to the general trend  ==!
            DO jj = 2 , jpjm1
               DO ji = fs_2, fs_jpim1  ! vector opt.
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + zsign * (  zftu(ji-1,jj,jk) - zftu(ji,jj,jk)       &
                     &                                           + zftv(ji,jj-1,jk) - zftv(ji,jj,jk)   )   &
                     &                                        / (  e1e2t(ji,jj) * e3t_n(ji,jj,jk)  )
               END DO
            END DO
            !
         END DO
         !
         !                                !==  add the vertical 33 flux  ==!
         IF( ln_traldf_lap ) THEN               ! laplacian case: eddy coef = ah_wslp2 - akz
            DO jk = 2, jpkm1       
               DO jj = 1, jpjm1
                  DO ji = fs_2, fs_jpim1
                     ztfw(ji,jj,jk) = ztfw(ji,jj,jk) - e1e2t(ji,jj) / e3w_n(ji,jj,jk) * tmask(ji,jj,jk)   &
                        &                            * ( ah_wslp2(ji,jj,jk) - akz(ji,jj,jk) )             &
                        &                            * ( ptb(ji,jj,jk-1,jn) - ptb(ji,jj,jk,jn) )
                  END DO
               END DO
            END DO
         ELSE                                   ! bilaplacian 
            SELECT CASE( kpass )
            CASE(  1  )                            ! 1st pass : eddy coef = ah_wslp2
               DO jk = 2, jpkm1 
                  DO jj = 1, jpjm1
                     DO ji = fs_2, fs_jpim1
                        ztfw(ji,jj,jk) = ztfw(ji,jj,jk) - e1e2t(ji,jj) / e3w_n(ji,jj,jk) * tmask(ji,jj,jk)             &
                           &                            * ah_wslp2(ji,jj,jk) * ( ptb(ji,jj,jk-1,jn) - ptb(ji,jj,jk,jn) )
                     END DO
                  END DO
               END DO 
            CASE(  2  )                            ! 2nd pass : eddy flux = ah_wslp2 and akz applied on ptb  and ptbb gradients, resp.
               DO jk = 2, jpkm1 
                  DO jj = 1, jpjm1
                     DO ji = fs_2, fs_jpim1
                        ztfw(ji,jj,jk) = ztfw(ji,jj,jk) - e1e2t(ji,jj) / e3w_n(ji,jj,jk) * tmask(ji,jj,jk)                      &
                           &                            * (  ah_wslp2(ji,jj,jk) * ( ptb (ji,jj,jk-1,jn) - ptb (ji,jj,jk,jn) )   &
                           &                               + akz     (ji,jj,jk) * ( ptbb(ji,jj,jk-1,jn) - ptbb(ji,jj,jk,jn) )   )
                     END DO
                  END DO
               END DO
            END SELECT 
         ENDIF
         !
         DO jk = 1, jpkm1                 !==  Divergence of vertical fluxes added to pta  ==!
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1  ! vector opt.
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + zsign * (  ztfw(ji,jj,jk+1) - ztfw(ji,jj,jk)  )   &
                     &                                        / ( e1e2t(ji,jj) * e3t_n(ji,jj,jk) )
               END DO
            END DO
         END DO
         !
         IF( ( kpass == 1 .AND. ln_traldf_lap ) .OR.  &     !==  first pass only (  laplacian)  ==!
             ( kpass == 2 .AND. ln_traldf_blp ) ) THEN      !==  2nd   pass      (bilaplacian)  ==!
            !
            !                          ! "Poleward" diffusive heat or salt transports (T-S case only)
            IF( l_ptr )  CALL dia_ptr_hst( jn, 'ldf', zftv(:,:,:)  )
            !                          ! Diffusive heat transports
            IF( l_hst )  CALL dia_ar5_hst( jn, 'ldf', zftu(:,:,:), zftv(:,:,:) )
            !
         ENDIF                                                    !== end pass selection  ==!
         !
         !                                                        ! ===============
      END DO                                                      ! end tracer loop
      !                                                           ! ===============
   END SUBROUTINE tra_ldf_triad

   !!==============================================================================
END MODULE traldf_triad
