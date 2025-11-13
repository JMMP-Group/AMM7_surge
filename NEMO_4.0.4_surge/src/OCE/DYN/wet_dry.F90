MODULE wet_dry

   !! includes updates to namelist namwad for diagnostic outputs of ROMS wetting and drying

   !!==============================================================================
   !!                       ***  MODULE  wet_dry  ***
   !! Wetting and drying includes initialisation routine and routines to
   !! compute and apply flux limiters and preserve water depth positivity
   !! only effects if wetting/drying is on (ln_wd_il == .true. or ln_wd_dl==.true. )
   !!==============================================================================
   !! History :  3.6  ! 2014-09  ((H.Liu)  Original code
   !!                 ! will add the runoff and periodic BC case later
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   wad_init      : initialisation of wetting and drying
   !!   wad_lmt       : horizontal flux limiter and limited velocity when wetting and drying happens
   !!   wad_lmt_bt    : same as wad_lmt for the barotropic stepping (dynspg_ts)
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce  , ONLY: ln_rnf   ! surface boundary condition: ocean
   USE sbcrnf         ! river runoff 
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! MPP library
   USE timing         ! timing of the main modules

   IMPLICIT NONE
   PRIVATE

   !!----------------------------------------------------------------------
   !! critical depths,filters, limiters,and masks for  Wetting and Drying
   !! ---------------------------------------------------------------------

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   wdmask   !: u- and v- limiter 
   !                                                           !  (can include negative depths)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   wdramp, wdrampu, wdrampv !: for hpg limiting

   LOGICAL,  PUBLIC  ::   ln_wd_il    !: Wetting/drying il activation switch (T:on,F:off)
   LOGICAL,  PUBLIC  ::   ln_wd_dl    !: Wetting/drying dl activation switch (T:on,F:off)
   REAL(wp), PUBLIC  ::   rn_wdmin0   !: depth at which wetting/drying starts
   REAL(wp), PUBLIC  ::   rn_wdmin1   !: minimum water depth on dried cells
   REAL(wp), PUBLIC  ::   r_rn_wdmin1 !: 1/minimum water depth on dried cells 
   REAL(wp), PUBLIC  ::   rn_wdmin2   !: tolerance of minimum water depth on dried cells
   REAL(wp), PUBLIC  ::   rn_wd_sbcdep   !: Depth at which to taper sbc fluxes
   REAL(wp), PUBLIC  ::   rn_wd_sbcfra   !: Fraction of SBC at taper depth
   REAL(wp), PUBLIC  ::   rn_wdld     !: land elevation below which wetting/drying will be considered
   INTEGER , PUBLIC  ::   nn_wdit     !: maximum number of iteration for W/D limiter
   LOGICAL,  PUBLIC  ::   ln_wd_dl_bc !: DL scheme: True implies 3D velocities are set to the barotropic values at points 
                                      !: where the flow is from wet points on less than half the barotropic sub-steps  
   LOGICAL,  PUBLIC  ::  ln_wd_dl_rmp !: use a ramp for the dl flux limiter between 2 rn_wdmin1 and rn_wdmin1 (rather than a cut-off at rn_wdmin1)      
   REAL(wp), PUBLIC  ::   ssh_ref     !: height of z=0 with respect to the geoid; 

   LOGICAL,  PUBLIC  ::   ll_wd       !: Wetting/drying activation switch if either ln_wd_il or ln_wd_dl

   PUBLIC   wad_init                  ! initialisation routine called by step.F90
   PUBLIC   wad_lmt                   ! routine called by sshwzv.F90
   PUBLIC   wad_lmt_bt                ! routine called by dynspg_ts.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE wad_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE wad_init  ***
      !!                    
      !! ** Purpose :   read wetting and drying namelist and print the variables.
      !!
      !! ** input   : - namwad namelist
      !!----------------------------------------------------------------------
      INTEGER  ::   ios, ierr   ! Local integer
      !!
      NAMELIST/namwad/ ln_wd_il, ln_wd_dl   , rn_wdmin0, rn_wdmin1, rn_wdmin2, rn_wdld,   &
         &             nn_wdit , ln_wd_dl_bc, ln_wd_dl_rmp, rn_wd_sbcdep,rn_wd_sbcfra
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namwad in reference namelist : Parameters for Wetting/Drying
      READ  ( numnam_ref, namwad, IOSTAT = ios, ERR = 905)
905   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namwad in reference namelist' ) 
      REWIND( numnam_cfg )              ! Namelist namwad in configuration namelist : Parameters for Wetting/Drying
      READ  ( numnam_cfg, namwad, IOSTAT = ios, ERR = 906)
906   IF( ios >  0 )   CALL ctl_nam ( ios , 'namwad in configuration namelist' )
      IF(lwm) WRITE ( numond, namwad )
      !
      IF( rn_wd_sbcfra>=1 )   CALL ctl_stop( 'STOP', 'rn_wd_sbcfra >=1 : must be < 1' )
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'wad_init : Wetting and drying initialization through namelist read'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '   Namelist namwad'
         WRITE(numout,*) '      Logical for Iter Lim wd option   ln_wd_il     = ', ln_wd_il
         WRITE(numout,*) '      Logical for Dir. Lim wd option   ln_wd_dl     = ', ln_wd_dl
         WRITE(numout,*) '      Depth at which wet/drying starts rn_wdmin0    = ', rn_wdmin0
         WRITE(numout,*) '      Minimum wet depth on dried cells rn_wdmin1    = ', rn_wdmin1
         WRITE(numout,*) '      Tolerance of min wet depth       rn_wdmin2    = ', rn_wdmin2
         WRITE(numout,*) '      land elevation threshold         rn_wdld      = ', rn_wdld
         WRITE(numout,*) '      Max iteration for W/D limiter    nn_wdit      = ', nn_wdit
         WRITE(numout,*) '      T => baroclinic u,v=0 at dry pts: ln_wd_dl_bc = ', ln_wd_dl_bc     
         WRITE(numout,*) '      use a ramp for rwd limiter:  ln_wd_dl_rwd_rmp = ', ln_wd_dl_rmp
         WRITE(numout,*) '      cut off depth sbc for wd   rn_wd_sbcdep       = ', rn_wd_sbcdep
         WRITE(numout,*) '      fraction to start sbc wgt rn_wd_sbcfra        = ', rn_wd_sbcfra
      ENDIF
      IF( .NOT. ln_read_cfg ) THEN
         IF(lwp) WRITE(numout,*) '      No configuration file so seting ssh_ref to zero  '
         ssh_ref=0._wp
      ENDIF

      r_rn_wdmin1 = 1 / rn_wdmin1
      ll_wd = .FALSE.
      IF( ln_wd_il .OR. ln_wd_dl ) THEN
         ll_wd = .TRUE.
         ALLOCATE( wdmask(jpi,jpj),   STAT=ierr )
         ALLOCATE( wdramp(jpi,jpj), wdrampu(jpi,jpj), wdrampv(jpi,jpj), STAT=ierr ) 
         IF( ierr /= 0 ) CALL ctl_stop('STOP', 'wad_init : Array allocation error')
      ENDIF
      !
   END SUBROUTINE wad_init


   SUBROUTINE wad_lmt( sshb1, sshemp, z2dt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE wad_lmt  ***
      !!                    
      !! ** Purpose :   generate flux limiters for wetting/drying
      !!
      !! ** Method  : - Prevent negative depth occurring (Not ready for Agrif) 
      !!
      !! ** Action  : - calculate flux limiter and W/D flag
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   sshb1        !!gm DOCTOR names: should start with p !
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   sshemp
      REAL(wp)                , INTENT(in   ) ::   z2dt
      !
      INTEGER  ::   ji, jj, jk, jk1     ! dummy loop indices
      INTEGER  ::   jflag               ! local scalar
      REAL(wp) ::   zcoef, zdep1, zdep2 ! local scalars
      REAL(wp) ::   zzflxp, zzflxn      ! local scalars
      REAL(wp) ::   zdepwd              ! local scalar, always wet cell depth
      REAL(wp) ::   ztmp                ! local scalars
      REAL(wp),  DIMENSION(jpi,jpj) ::   zwdlmtu, zwdlmtv   ! W/D flux limiters
      REAL(wp),  DIMENSION(jpi,jpj) ::   zflxp  ,  zflxn    ! local 2D workspace
      REAL(wp),  DIMENSION(jpi,jpj) ::   zflxu  ,  zflxv    ! local 2D workspace
      REAL(wp),  DIMENSION(jpi,jpj) ::   zflxu1 , zflxv1    ! local 2D workspace
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('wad_lmt')      !
      !
      DO jk = 1, jpkm1
         un(:,:,jk) = un(:,:,jk) * zwdlmtu(:,:) 
         vn(:,:,jk) = vn(:,:,jk) * zwdlmtv(:,:) 
      END DO
      jflag  = 0
      zdepwd = 50._wp      ! maximum depth on which that W/D could possibly happen
      !
      zflxp(:,:)   = 0._wp
      zflxn(:,:)   = 0._wp
      zflxu(:,:)   = 0._wp
      zflxv(:,:)   = 0._wp
      !
      zwdlmtu(:,:) = 1._wp
      zwdlmtv(:,:) = 1._wp
      !
      DO jk = 1, jpkm1     ! Horizontal Flux in u and v direction
         zflxu(:,:) = zflxu(:,:) + e3u_n(:,:,jk) * un(:,:,jk) * umask(:,:,jk)
         zflxv(:,:) = zflxv(:,:) + e3v_n(:,:,jk) * vn(:,:,jk) * vmask(:,:,jk)
      END DO
      zflxu(:,:) = zflxu(:,:) * e2u(:,:)
      zflxv(:,:) = zflxv(:,:) * e1v(:,:)
      !
      wdmask(:,:) = 1._wp
      DO jj = 2, jpj
         DO ji = 2, jpi 
            !
            IF( tmask(ji,jj,1)        < 0.5_wp )   CYCLE    ! we don't care about land cells
            IF( ht_0(ji,jj) - ssh_ref > zdepwd )   CYCLE    ! and cells which are unlikely to dry
            !
            zflxp(ji,jj) = MAX( zflxu(ji,jj) , 0._wp ) - MIN( zflxu(ji-1,jj  ) , 0._wp )   &
               &         + MAX( zflxv(ji,jj) , 0._wp ) - MIN( zflxv(ji,  jj-1) , 0._wp ) 
            zflxn(ji,jj) = MIN( zflxu(ji,jj) , 0._wp ) - MAX( zflxu(ji-1,jj  ) , 0._wp )   &
               &         + MIN( zflxv(ji,jj) , 0._wp ) - MAX( zflxv(ji,  jj-1) , 0._wp ) 
            !
            zdep2 = ht_0(ji,jj) + sshb1(ji,jj) - rn_wdmin1
            IF( zdep2 <= 0._wp ) THEN     ! add more safty, but not necessary
               sshb1(ji,jj) = rn_wdmin1 - ht_0(ji,jj)
               IF(zflxu(ji,  jj) > 0._wp) zwdlmtu(ji  ,jj) = 0._wp
               IF(zflxu(ji-1,jj) < 0._wp) zwdlmtu(ji-1,jj) = 0._wp
               IF(zflxv(ji,  jj) > 0._wp) zwdlmtv(ji  ,jj) = 0._wp
               IF(zflxv(ji,jj-1) < 0._wp) zwdlmtv(ji,jj-1) = 0._wp 
               wdmask(ji,jj) = 0._wp
            END IF
         END DO
      END DO
      !
      !           ! HPG limiter from jholt
      wdramp(:,:) = min((ht_0(:,:) + sshb1(:,:) - rn_wdmin1)/(rn_wdmin0 - rn_wdmin1),1.0_wp)
      !jth assume don't need a lbc_lnk here
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
            wdrampu(ji,jj) = MIN( wdramp(ji,jj) , wdramp(ji+1,jj) )
            wdrampv(ji,jj) = MIN( wdramp(ji,jj) , wdramp(ji,jj+1) )
         END DO
      END DO
      !           ! end HPG limiter
      !
      !
      DO jk1 = 1, nn_wdit + 1      !==  start limiter iterations  ==!
         !
         zflxu1(:,:) = zflxu(:,:) * zwdlmtu(:,:)
         zflxv1(:,:) = zflxv(:,:) * zwdlmtv(:,:)
         jflag = 0     ! flag indicating if any further iterations are needed
         !
         DO jj = 2, jpj
            DO ji = 2, jpi 
               IF( tmask(ji, jj, 1) < 0.5_wp )   CYCLE 
               IF( ht_0(ji,jj)      > zdepwd )   CYCLE
               !
               ztmp = e1e2t(ji,jj)
               !
               zzflxp = MAX( zflxu1(ji,jj) , 0._wp ) - MIN( zflxu1(ji-1,jj  ) , 0._wp)   &
                  &   + MAX( zflxv1(ji,jj) , 0._wp ) - MIN( zflxv1(ji,  jj-1) , 0._wp) 
               zzflxn = MIN( zflxu1(ji,jj) , 0._wp ) - MAX( zflxu1(ji-1,jj  ) , 0._wp)   &
                  &   + MIN( zflxv1(ji,jj) , 0._wp ) - MAX( zflxv1(ji,  jj-1) , 0._wp) 
               !
               zdep1 = (zzflxp + zzflxn) * z2dt / ztmp
               zdep2 = ht_0(ji,jj) + sshb1(ji,jj) - rn_wdmin1 - z2dt * sshemp(ji,jj)
               !
               IF( zdep1 > zdep2 ) THEN
                  wdmask(ji, jj) = 0._wp
                  zcoef = ( ( zdep2 - rn_wdmin2 ) * ztmp - zzflxn * z2dt ) / ( zflxp(ji,jj) * z2dt )
                  !zcoef = ( ( zdep2 - rn_wdmin2 ) * ztmp - zzflxn * z2dt ) / ( zzflxp * z2dt )
                  ! flag if the limiter has been used but stop flagging if the only
                  ! changes have zeroed the coefficient since further iterations will
                  ! not change anything
                  IF( zcoef > 0._wp ) THEN   ;   jflag = 1 
                  ELSE                       ;   zcoef = 0._wp
                  ENDIF
                  IF( jk1 > nn_wdit )   zcoef = 0._wp
                  IF( zflxu1(ji  ,jj  ) > 0._wp )   zwdlmtu(ji  ,jj  ) = zcoef
                  IF( zflxu1(ji-1,jj  ) < 0._wp )   zwdlmtu(ji-1,jj  ) = zcoef
                  IF( zflxv1(ji  ,jj  ) > 0._wp )   zwdlmtv(ji  ,jj  ) = zcoef
                  IF( zflxv1(ji  ,jj-1) < 0._wp )   zwdlmtv(ji  ,jj-1) = zcoef
               ENDIF
            END DO
         END DO
         CALL lbc_lnk_multi( 'wet_dry', zwdlmtu, 'U', 1., zwdlmtv, 'V', 1. )
         !
         CALL mpp_max('wet_dry', jflag)   !max over the global domain
         !
         IF( jflag == 0 )   EXIT
         !
      END DO  ! jk1 loop
      !
      DO jk = 1, jpkm1
         un(:,:,jk) = un(:,:,jk) * zwdlmtu(:,:) 
         vn(:,:,jk) = vn(:,:,jk) * zwdlmtv(:,:) 
      END DO
      un_b(:,:) = un_b(:,:) * zwdlmtu(:, :)
      vn_b(:,:) = vn_b(:,:) * zwdlmtv(:, :)
      !
!!gm TO BE SUPPRESSED ?  these lbc_lnk are useless since zwdlmtu and zwdlmtv are defined everywhere !
      CALL lbc_lnk_multi( 'wet_dry', un  , 'U', -1., vn  , 'V', -1. )
      CALL lbc_lnk_multi( 'wet_dry', un_b, 'U', -1., vn_b, 'V', -1. )
!!gm
      !
      IF(jflag == 1 .AND. lwp)   WRITE(numout,*) 'Need more iterations in wad_lmt!!!'
      !
      !IF( ln_rnf      )   CALL sbc_rnf_div( hdivn )          ! runoffs (update hdivn field)
      !
      IF( ln_timing )   CALL timing_stop('wad_lmt')      !
      !
   END SUBROUTINE wad_lmt


   SUBROUTINE wad_lmt_bt( zflxu, zflxv, sshn_e, zssh_frc, rdtbt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE wad_lmt  ***
      !!                    
      !! ** Purpose :   limiting flux in the barotropic stepping (dynspg_ts)
      !!
      !! ** Method  : - Prevent negative depth occurring (Not ready for Agrif) 
      !!
      !! ** Action  : - calculate flux limiter and W/D flag
      !!----------------------------------------------------------------------
      REAL(wp)                , INTENT(in   ) ::   rdtbt    ! ocean time-step index
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   zflxu,  zflxv, sshn_e, zssh_frc  
      !
      INTEGER  ::   ji, jj, jk, jk1     ! dummy loop indices
      INTEGER  ::   jflag               ! local integer
      REAL(wp) ::   z2dt
      REAL(wp) ::   zcoef, zdep1, zdep2 ! local scalars
      REAL(wp) ::   zzflxp, zzflxn      ! local scalars
      REAL(wp) ::   zdepwd              ! local scalar, always wet cell depth
      REAL(wp) ::   ztmp                ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   zwdlmtu, zwdlmtv         !: W/D flux limiters
      REAL(wp), DIMENSION(jpi,jpj) ::   zflxp,  zflxn            ! local 2D workspace
      REAL(wp), DIMENSION(jpi,jpj) ::   zflxu1, zflxv1           ! local 2D workspace
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('wad_lmt_bt')      !
      !
      jflag  = 0
      zdepwd = 50._wp   ! maximum depth that ocean cells can have W/D processes
      !
      z2dt = rdtbt   
      !
      zflxp(:,:)   = 0._wp
      zflxn(:,:)   = 0._wp
      zwdlmtu(:,:) = 1._wp
      zwdlmtv(:,:) = 1._wp
      !
      DO jj = 2, jpj      ! Horizontal Flux in u and v direction
         DO ji = 2, jpi 
            !
            IF( tmask(ji, jj, 1 ) < 0.5_wp) CYCLE   ! we don't care about land cells
            IF( ht_0(ji,jj) > zdepwd )      CYCLE   ! and cells which are unlikely to dry
            !
            zflxp(ji,jj) = MAX( zflxu(ji,jj) , 0._wp ) - MIN( zflxu(ji-1,jj  ) , 0._wp )   &
               &         + MAX( zflxv(ji,jj) , 0._wp ) - MIN( zflxv(ji,  jj-1) , 0._wp ) 
            zflxn(ji,jj) = MIN( zflxu(ji,jj) , 0._wp ) - MAX( zflxu(ji-1,jj  ) , 0._wp )   &
               &         + MIN( zflxv(ji,jj) , 0._wp ) - MAX( zflxv(ji,  jj-1) , 0._wp ) 
            !
            zdep2 = ht_0(ji,jj) + sshn_e(ji,jj) - rn_wdmin1
            IF( zdep2 <= 0._wp ) THEN  !add more safety, but not necessary
              sshn_e(ji,jj) = rn_wdmin1 - ht_0(ji,jj)
              IF( zflxu(ji  ,jj  ) > 0._wp)   zwdlmtu(ji  ,jj  ) = 0._wp
              IF( zflxu(ji-1,jj  ) < 0._wp)   zwdlmtu(ji-1,jj  ) = 0._wp
              IF( zflxv(ji  ,jj  ) > 0._wp)   zwdlmtv(ji  ,jj  ) = 0._wp
              IF( zflxv(ji  ,jj-1) < 0._wp)   zwdlmtv(ji  ,jj-1) = 0._wp 
            ENDIF
         END DO
      END DO
      !
      DO jk1 = 1, nn_wdit + 1      !! start limiter iterations 
         !
         zflxu1(:,:) = zflxu(:,:) * zwdlmtu(:,:)
         zflxv1(:,:) = zflxv(:,:) * zwdlmtv(:,:)
         jflag = 0     ! flag indicating if any further iterations are needed
         !
         DO jj = 2, jpj
            DO ji = 2, jpi 
               !
               IF( tmask(ji, jj, 1 ) < 0.5_wp )   CYCLE 
               IF( ht_0(ji,jj)       > zdepwd )   CYCLE
               !
               ztmp = e1e2t(ji,jj)
               !
               zzflxp = max(zflxu1(ji,jj), 0._wp) - min(zflxu1(ji-1,jj),   0._wp)   &
                  &   + max(zflxv1(ji,jj), 0._wp) - min(zflxv1(ji,  jj-1), 0._wp) 
               zzflxn = min(zflxu1(ji,jj), 0._wp) - max(zflxu1(ji-1,jj),   0._wp)   &
                  &   + min(zflxv1(ji,jj), 0._wp) - max(zflxv1(ji,  jj-1), 0._wp) 
          
               zdep1 = (zzflxp + zzflxn) * z2dt / ztmp
               zdep2 = ht_0(ji,jj) + sshn_e(ji,jj) - rn_wdmin1 - z2dt * zssh_frc(ji,jj)
          
               IF(zdep1 > zdep2) THEN
                 zcoef = ( ( zdep2 - rn_wdmin2 ) * ztmp - zzflxn * z2dt ) / ( zflxp(ji,jj) * z2dt )
                 !zcoef = ( ( zdep2 - rn_wdmin2 ) * ztmp - zzflxn * z2dt ) / ( zzflxp * z2dt )
                 ! flag if the limiter has been used but stop flagging if the only
                 ! changes have zeroed the coefficient since further iterations will
                 ! not change anything
                 IF( zcoef > 0._wp ) THEN
                    jflag = 1 
                 ELSE
                    zcoef = 0._wp
                 ENDIF
                 IF(jk1 > nn_wdit) zcoef = 0._wp
                 IF(zflxu1(ji,  jj) > 0._wp) zwdlmtu(ji  ,jj) = zcoef
                 IF(zflxu1(ji-1,jj) < 0._wp) zwdlmtu(ji-1,jj) = zcoef
                 IF(zflxv1(ji,  jj) > 0._wp) zwdlmtv(ji  ,jj) = zcoef
                 IF(zflxv1(ji,jj-1) < 0._wp) zwdlmtv(ji,jj-1) = zcoef
               END IF
            END DO ! ji loop
         END DO  ! jj loop
         !
         CALL lbc_lnk_multi( 'wet_dry', zwdlmtu, 'U', 1., zwdlmtv, 'V', 1. )
         !
         CALL mpp_max('wet_dry', jflag)   !max over the global domain
         !
         IF(jflag == 0)   EXIT
         !
      END DO  ! jk1 loop
      !
      zflxu(:,:) = zflxu(:,:) * zwdlmtu(:, :) 
      zflxv(:,:) = zflxv(:,:) * zwdlmtv(:, :) 
      !
!!gm THIS lbc_lnk is useless since it is already done at the end of the jk1-loop
      CALL lbc_lnk_multi( 'wet_dry', zflxu, 'U', -1., zflxv, 'V', -1. )
!!gm end
      !
      IF( jflag == 1 .AND. lwp )   WRITE(numout,*) 'Need more iterations in wad_lmt_bt!!!'
      !
      !IF( ln_rnf      )   CALL sbc_rnf_div( hdivn )          ! runoffs (update hdivn field)
      !
      IF( ln_timing )   CALL timing_stop('wad_lmt_bt')      !
      !
   END SUBROUTINE wad_lmt_bt

   !!==============================================================================
END MODULE wet_dry
