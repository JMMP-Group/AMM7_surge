MODULE dynspg_ts

   !! Includes ROMS wd scheme with diagnostic outputs ; un and ua updates are commented out ! 

   !!======================================================================
   !!                   ***  MODULE  dynspg_ts  ***
   !! Ocean dynamics:  surface pressure gradient trend, split-explicit scheme
   !!======================================================================
   !! History :   1.0  ! 2004-12  (L. Bessieres, G. Madec)  Original code
   !!              -   ! 2005-11  (V. Garnier, G. Madec)  optimization
   !!              -   ! 2006-08  (S. Masson)  distributed restart using iom
   !!             2.0  ! 2007-07  (D. Storkey) calls to BDY routines
   !!              -   ! 2008-01  (R. Benshila)  change averaging method
   !!             3.2  ! 2009-07  (R. Benshila, G. Madec) Complete revisit associated to vvl reactivation
   !!             3.3  ! 2010-09  (D. Storkey, E. O'Dea) update for BDY for Shelf configurations
   !!             3.3  ! 2011-03  (R. Benshila, R. Hordoir, P. Oddo) update calculation of ub_b
   !!             3.5  ! 2013-07  (J. Chanut) Switch to Forward-backward time stepping
   !!             3.6  ! 2013-11  (A. Coward) Update for z-tilde compatibility
   !!             3.7  ! 2015-11  (J. Chanut) free surface simplification
   !!              -   ! 2016-12  (G. Madec, E. Clementi) update for Stoke-Drift divergence
   !!             4.0  ! 2017-05  (G. Madec)  drag coef. defined at t-point (zdfdrg.F90)
   !!---------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_spg_ts     : compute surface pressure gradient trend using a time-splitting scheme 
   !!   dyn_spg_ts_init: initialisation of the time-splitting scheme
   !!   ts_wgt         : set time-splitting weights for temporal averaging (or not)
   !!   ts_rst         : read/write time-splitting fields in restart file
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean
   USE zdf_oce         ! vertical physics: variables
   USE zdfdrg          ! vertical physics: top/bottom drag coef.
   USE sbcisf          ! ice shelf variable (fwfisf)
   USE sbcapr          ! surface boundary condition: atmospheric pressure
   USE dynadv    , ONLY: ln_dynadv_vec
   USE dynvor          ! vortivity scheme indicators
   USE phycst          ! physical constants
   USE dynvor          ! vorticity term
   USE wet_dry         ! wetting/drying flux limter
   USE bdy_oce         ! open boundary
   USE bdyvol          ! open boundary volume conservation
   USE bdytides        ! open boundary condition data
   USE bdydyn2d        ! open boundary conditions on barotropic variables
   USE sbctide         ! tides
   USE updtide         ! tide potential
   USE sbcwave         ! surface wave
#if defined key_agrif
   USE agrif_oce_interp ! agrif
   USE agrif_oce
#endif
#if defined key_asminc   
   USE asminc          ! Assimilation increment
#endif
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE iom             ! IOM library
   USE restart         ! only for lrst_oce

   USE iom   ! to remove

   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_spg_ts        ! called by dyn_spg 
   PUBLIC dyn_spg_ts_init   !    -    - dyn_spg_init

   !! Time filtered arrays at baroclinic time step:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   un_adv , vn_adv   !: Advection vel. at "now" barocl. step
   !
   INTEGER, SAVE :: icycle      ! Number of barotropic sub-steps for each internal step nn_baro <= 2.5 nn_baro
   REAL(wp),SAVE :: rdtbt       ! Barotropic time step
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)   ::   wgtbtp1, wgtbtp2   ! 1st & 2nd weights used in time filtering of barotropic fields
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zwz                ! ff_f/h at F points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ftnw, ftne         ! triad of coriolis parameter
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ftsw, ftse         ! (only used with een vorticity scheme)

   REAL(wp) ::   r1_12 = 1._wp / 12._wp   ! local ratios
   REAL(wp) ::   r1_8  = 0.125_wp         !
   REAL(wp) ::   r1_4  = 0.25_wp          !
   REAL(wp) ::   r1_2  = 0.5_wp           !

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dyn_spg_ts_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_ts_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(3)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( wgtbtp1(3*nn_baro), wgtbtp2(3*nn_baro), zwz(jpi,jpj), STAT=ierr(1) )
      IF( ln_dynvor_een .OR. ln_dynvor_eeT )   &
         &     ALLOCATE( ftnw(jpi,jpj) , ftne(jpi,jpj) , ftsw(jpi,jpj) , ftse(jpi,jpj), STAT=ierr(2)   )
         !
      ALLOCATE( un_adv(jpi,jpj), vn_adv(jpi,jpj)                    , STAT=ierr(3) )
      !
      dyn_spg_ts_alloc = MAXVAL( ierr(:) )
      !
      CALL mpp_sum( 'dynspg_ts', dyn_spg_ts_alloc )
      IF( dyn_spg_ts_alloc /= 0 )   CALL ctl_stop( 'STOP', 'dyn_spg_ts_alloc: failed to allocate arrays' )
      !
   END FUNCTION dyn_spg_ts_alloc


   SUBROUTINE dyn_spg_ts( kt )
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose : - Compute the now trend due to the explicit time stepping
      !!              of the quasi-linear barotropic system, and add it to the
      !!              general momentum trend. 
      !!
      !! ** Method  : - split-explicit schem (time splitting) :
      !!      Barotropic variables are advanced from internal time steps
      !!      "n"   to "n+1" if ln_bt_fw=T
      !!      or from 
      !!      "n-1" to "n+1" if ln_bt_fw=F
      !!      thanks to a generalized forward-backward time stepping (see ref. below).
      !!
      !! ** Action :
      !!      -Update the filtered free surface at step "n+1"      : ssha
      !!      -Update filtered barotropic velocities at step "n+1" : ua_b, va_b
      !!      -Compute barotropic advective fluxes at step "n"     : un_adv, vn_adv
      !!      These are used to advect tracers and are compliant with discrete
      !!      continuity equation taken at the baroclinic time steps. This 
      !!      ensures tracers conservation.
      !!      - (ua, va) momentum trend updated with barotropic component.
      !!
      !! References : Shchepetkin and McWilliams, Ocean Modelling, 2005. 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk, jn        ! dummy loop indices
      LOGICAL  ::   ll_fw_start           ! =T : forward integration 
      LOGICAL  ::   ll_init               ! =T : special startup of 2d equations
      INTEGER  ::   noffset               ! local integers  : time offset for bdy update
      REAL(wp) ::   r1_2dt_b, z1_hu, z1_hv          ! local scalars
      REAL(wp) ::   za0, za1, za2, za3              !   -      -
      REAL(wp) ::   zztmp, zldg               !   -      -
      REAL(wp) ::   zhu_bck, zhv_bck, zhdiv         !   -      -
      REAL(wp) ::   zun_save, zvn_save              !   -      -
      REAL(wp), DIMENSION(jpi,jpj) :: zu_trd, zu_frc, zu_spg, zssh_frc
      REAL(wp), DIMENSION(jpi,jpj) :: zv_trd, zv_frc, zv_spg
      REAL(wp), DIMENSION(jpi,jpj) :: zsshu_a, zhup2_e, zhtp2_e
      REAL(wp), DIMENSION(jpi,jpj) :: zsshv_a, zhvp2_e, zsshp2_e
      REAL(wp), DIMENSION(jpi,jpj) :: zCdU_u, zCdU_v   ! top/bottom stress at u- & v-points
      REAL(wp), DIMENSION(jpi,jpj) :: zhU, zhV         ! fluxes
      !
      REAL(wp) ::   zwdramp                     ! local scalar - only used if ln_wd_dl = .True. 

      INTEGER  :: iwdg, jwdg, kwdg   ! short-hand values for the indices of the output point

      REAL(wp) ::   zepsilon, zgamma            !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zcpx, zcpy   ! Wetting/Dying gravity filter coef.
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ztwdmask, zuwdmask, zvwdmask ! ROMS wetting and drying masks at t,u,v points
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zuwdav2, zvwdav2    ! averages over the sub-steps of zuwdmask and zvwdmask
      !!----------------------------------------------------------------------
      !
      IF( ln_wd_il ) ALLOCATE( zcpx(jpi,jpj), zcpy(jpi,jpj) )
      !                                         !* Allocate temporary arrays
      IF( ln_wd_dl ) ALLOCATE( ztwdmask(jpi,jpj), zuwdmask(jpi,jpj), zvwdmask(jpi,jpj), zuwdav2(jpi,jpj), zvwdav2(jpi,jpj))
      !
      zwdramp = r_rn_wdmin1               ! simplest ramp 
!     zwdramp = 1._wp / (rn_wdmin2 - rn_wdmin1) ! more general ramp
      !                                         ! inverse of baroclinic time step 
      IF( kt == nit000 .AND. neuler == 0 ) THEN   ;   r1_2dt_b = 1._wp / (         rdt )
      ELSE                                        ;   r1_2dt_b = 1._wp / ( 2._wp * rdt )
      ENDIF
      !
      ll_init     = ln_bt_av                    ! if no time averaging, then no specific restart 
      ll_fw_start = .FALSE.
      !                                         ! time offset in steps for bdy data update
      IF( .NOT.ln_bt_fw ) THEN   ;   noffset = - nn_baro
      ELSE                       ;   noffset =   0 
      ENDIF
      !
      IF( kt == nit000 ) THEN                   !* initialisation
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_ts : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~   free surface with time splitting'
         IF(lwp) WRITE(numout,*)
         !
         IF( neuler == 0 )   ll_init=.TRUE.
         !
         IF( ln_bt_fw .OR. neuler == 0 ) THEN
            ll_fw_start =.TRUE.
            noffset     = 0
         ELSE
            ll_fw_start =.FALSE.
         ENDIF
         !                    ! Set averaging weights and cycle length:
         CALL ts_wgt( ln_bt_av, ll_fw_start, icycle, wgtbtp1, wgtbtp2 )
         !
      ENDIF
      !
      ! If forward start at previous time step, and centered integration, 
      ! then update averaging weights:
      IF (.NOT.ln_bt_fw .AND.( neuler==0 .AND. kt==nit000+1 ) ) THEN
         ll_fw_start=.FALSE.
         CALL ts_wgt( ln_bt_av, ll_fw_start, icycle, wgtbtp1, wgtbtp2 )
      ENDIF
      !
                          
      ! -----------------------------------------------------------------------------
      !  Phase 1 : Coupling between general trend and barotropic estimates (1st step)
      ! -----------------------------------------------------------------------------
      !      
      !
      !                                   !=  zu_frc =  1/H e3*d/dt(Ua)  =!  (Vertical mean of Ua, the 3D trends)
      !                                   !  ---------------------------  !
      zu_frc(:,:) = SUM( e3u_n(:,:,:) * ua(:,:,:) * umask(:,:,:) , DIM=3 ) * r1_hu_n(:,:)
      zv_frc(:,:) = SUM( e3v_n(:,:,:) * va(:,:,:) * vmask(:,:,:) , DIM=3 ) * r1_hv_n(:,:)
      !
      !
      !                                   !=  Ua => baroclinic trend  =!   (remove its vertical mean)
      DO jk = 1, jpkm1                    !  ------------------------  !
         ua(:,:,jk) = ( ua(:,:,jk) - zu_frc(:,:) ) * umask(:,:,jk)
         va(:,:,jk) = ( va(:,:,jk) - zv_frc(:,:) ) * vmask(:,:,jk)
      END DO
      
!!gm  Question here when removing the Vertically integrated trends, we remove the vertically integrated NL trends on momentum....
!!gm  Is it correct to do so ?   I think so...
      
      !                                   !=  remove 2D Coriolis and pressure gradient trends  =!
      !                                   !  -------------------------------------------------  !
      !
      IF( kt == nit000 .OR. .NOT. ln_linssh )   CALL dyn_cor_2D_init   ! Set zwz, the barotropic Coriolis force coefficient
      !       ! recompute zwz = f/depth  at every time step for (.NOT.ln_linssh) as the water colomn height changes
      !
      !                                         !* 2D Coriolis trends
      zhU(:,:) = un_b(:,:) * hu_n(:,:) * e2u(:,:)        ! now fluxes 
      zhV(:,:) = vn_b(:,:) * hv_n(:,:) * e1v(:,:)        ! NB: FULL domain : put a value in last row and column
      !
      CALL dyn_cor_2d( hu_n, hv_n, un_b, vn_b, zhU, zhV,  &   ! <<== in
         &                               zu_trd, zv_trd   )   ! ==>> out
      !
      IF( .NOT.ln_linssh ) THEN                 !* surface pressure gradient   (variable volume only)
         !
         IF( ln_wd_il ) THEN                       ! W/D : limiter applied to spgspg
            CALL wad_spg( sshn, zcpx, zcpy )          ! Calculating W/D gravity filters, zcpx and zcpy
            DO jj = 2, jpjm1
               DO ji = 2, jpim1                ! SPG with the application of W/D gravity filters
                  zu_trd(ji,jj) = zu_trd(ji,jj) - grav * ( sshn(ji+1,jj  ) - sshn(ji  ,jj ) )   &
                     &                          * r1_e1u(ji,jj) * zcpx(ji,jj)  * wdrampu(ji,jj)  !jth
                  zv_trd(ji,jj) = zv_trd(ji,jj) - grav * ( sshn(ji  ,jj+1) - sshn(ji  ,jj ) )   &
                     &                          * r1_e2v(ji,jj) * zcpy(ji,jj)  * wdrampv(ji,jj)  !jth
               END DO
            END DO
         ELSE                                      ! now suface pressure gradient
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zu_trd(ji,jj) = zu_trd(ji,jj) - grav * (  sshn(ji+1,jj  ) - sshn(ji  ,jj  )  ) * r1_e1u(ji,jj)
                  zv_trd(ji,jj) = zv_trd(ji,jj) - grav * (  sshn(ji  ,jj+1) - sshn(ji  ,jj  )  ) * r1_e2v(ji,jj) 
               END DO
            END DO
         ENDIF
         !
      ENDIF
      !
      DO jj = 2, jpjm1                          ! Remove coriolis term (and possibly spg) from barotropic trend
         DO ji = fs_2, fs_jpim1
             zu_frc(ji,jj) = zu_frc(ji,jj) - zu_trd(ji,jj) * ssumask(ji,jj)
             zv_frc(ji,jj) = zv_frc(ji,jj) - zv_trd(ji,jj) * ssvmask(ji,jj)
          END DO
      END DO 
      !
      !                                   !=  Add bottom stress contribution from baroclinic velocities  =!
      !                                   !  -----------------------------------------------------------  !
      CALL dyn_drg_init( zu_frc, zv_frc,  zCdU_u, zCdU_v )      ! also provide the barotropic drag coefficients
      !
      !                                   !=  Add atmospheric pressure forcing  =!
      !                                   !  ----------------------------------  !
      IF( ln_apr_dyn ) THEN
         IF( ln_bt_fw ) THEN                          ! FORWARD integration: use kt+1/2 pressure (NOW+1/2)
            DO jj = 2, jpjm1              
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zu_frc(ji,jj) = zu_frc(ji,jj) + grav * (  ssh_ib (ji+1,jj  ) - ssh_ib (ji,jj) ) * r1_e1u(ji,jj)
                  zv_frc(ji,jj) = zv_frc(ji,jj) + grav * (  ssh_ib (ji  ,jj+1) - ssh_ib (ji,jj) ) * r1_e2v(ji,jj)
               END DO
            END DO
         ELSE                                         ! CENTRED integration: use kt-1/2 + kt+1/2 pressure (NOW)
            zztmp = grav * r1_2
            DO jj = 2, jpjm1              
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zu_frc(ji,jj) = zu_frc(ji,jj) + zztmp * (  ssh_ib (ji+1,jj  ) - ssh_ib (ji,jj)  &
                       &                                   + ssh_ibb(ji+1,jj  ) - ssh_ibb(ji,jj)  ) * r1_e1u(ji,jj)
                  zv_frc(ji,jj) = zv_frc(ji,jj) + zztmp * (  ssh_ib (ji  ,jj+1) - ssh_ib (ji,jj)  &
                       &                                   + ssh_ibb(ji  ,jj+1) - ssh_ibb(ji,jj)  ) * r1_e2v(ji,jj)
               END DO
            END DO
         ENDIF 
      ENDIF
      !
      !                                   !=  Add atmospheric pressure forcing  =!
      !                                   !  ----------------------------------  !
      IF( ln_bt_fw ) THEN                        ! Add wind forcing
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zu_frc(ji,jj) =  zu_frc(ji,jj) + r1_rau0 * utau(ji,jj) * r1_hu_n(ji,jj)
               zv_frc(ji,jj) =  zv_frc(ji,jj) + r1_rau0 * vtau(ji,jj) * r1_hv_n(ji,jj)
            END DO
         END DO
      ELSE
         zztmp = r1_rau0 * r1_2
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zu_frc(ji,jj) =  zu_frc(ji,jj) + zztmp * ( utau_b(ji,jj) + utau(ji,jj) ) * r1_hu_n(ji,jj)
               zv_frc(ji,jj) =  zv_frc(ji,jj) + zztmp * ( vtau_b(ji,jj) + vtau(ji,jj) ) * r1_hv_n(ji,jj)
            END DO
         END DO
      ENDIF  
      !
      !              !----------------!
      !              !==  sssh_frc  ==!   Right-Hand-Side of the barotropic ssh equation   (over the FULL domain)
      !              !----------------!
      !                                   !=  Net water flux forcing applied to a water column  =!
      !                                   ! ---------------------------------------------------  !
      IF (ln_bt_fw) THEN                          ! FORWARD integration: use kt+1/2 fluxes (NOW+1/2)
         zssh_frc(:,:) = r1_rau0 * ( emp(:,:)             - rnf(:,:)              + fwfisf(:,:)                  )
      ELSE                                        ! CENTRED integration: use kt-1/2 + kt+1/2 fluxes (NOW)
         zztmp = r1_rau0 * r1_2
         zssh_frc(:,:) = zztmp * (  emp(:,:) + emp_b(:,:) - rnf(:,:) - rnf_b(:,:) + fwfisf(:,:) + fwfisf_b(:,:)  )
      ENDIF
      !                                   !=  Add Stokes drift divergence  =!   (if exist)
      IF( ln_sdw ) THEN                   !  -----------------------------  !
         zssh_frc(:,:) = zssh_frc(:,:) + div_sd(:,:)
      ENDIF
      !
#if defined key_asminc
      !                                   !=  Add the IAU weighted SSH increment  =!
      !                                   !  ------------------------------------  !
      IF( lk_asminc .AND. ln_sshinc .AND. ln_asmiau ) THEN
         zssh_frc(:,:) = zssh_frc(:,:) - ssh_iau(:,:)
      ENDIF
#endif
      !                                   != Fill boundary data arrays for AGRIF
      !                                   ! ------------------------------------
#if defined key_agrif
         IF( .NOT.Agrif_Root() ) CALL agrif_dta_ts( kt )
#endif
      !
      ! -----------------------------------------------------------------------
      !  Phase 2 : Integration of the barotropic equations 
      ! -----------------------------------------------------------------------
      !
      !                                             ! ==================== !
      !                                             !    Initialisations   !
      !                                             ! ==================== !  
      ! Initialize barotropic variables:      
      IF( ll_init )THEN
         sshbb_e(:,:) = 0._wp
         ubb_e  (:,:) = 0._wp
         vbb_e  (:,:) = 0._wp
         sshb_e (:,:) = 0._wp
         ub_e   (:,:) = 0._wp
         vb_e   (:,:) = 0._wp
      ENDIF
      !
      IF( ln_linssh ) THEN    ! mid-step ocean depth is fixed (hup2_e=hu_n=hu_0)
         zhup2_e(:,:) = hu_n(:,:)
         zhvp2_e(:,:) = hv_n(:,:)
         zhtp2_e(:,:) = ht_n(:,:)
      ENDIF
      !
      IF (ln_bt_fw) THEN                  ! FORWARD integration: start from NOW fields                    
         sshn_e(:,:) =    sshn(:,:)            
         un_e  (:,:) =    un_b(:,:)            
         vn_e  (:,:) =    vn_b(:,:)
         !
         hu_e  (:,:) =    hu_n(:,:)       
         hv_e  (:,:) =    hv_n(:,:) 
         hur_e (:,:) = r1_hu_n(:,:)    
         hvr_e (:,:) = r1_hv_n(:,:)
      ELSE                                ! CENTRED integration: start from BEFORE fields
         sshn_e(:,:) =    sshb(:,:)
         un_e  (:,:) =    ub_b(:,:)         
         vn_e  (:,:) =    vb_b(:,:)
         !
         hu_e  (:,:) =    hu_b(:,:)       
         hv_e  (:,:) =    hv_b(:,:) 
         hur_e (:,:) = r1_hu_b(:,:)    
         hvr_e (:,:) = r1_hv_b(:,:)
      ENDIF
      !
      ! Initialize sums:
      ua_b  (:,:) = 0._wp       ! After barotropic velocities (or transport if flux form)          
      va_b  (:,:) = 0._wp
      ssha  (:,:) = 0._wp       ! Sum for after averaged sea level
      un_adv(:,:) = 0._wp       ! Sum for now transport issued from ts loop
      vn_adv(:,:) = 0._wp
      !
      IF( ln_wd_dl ) THEN
         zuwdmask(:,:) = 0._wp  ! set to zero for definiteness (not sure this is necessary) 
         zvwdmask(:,:) = 0._wp  ! 
         zuwdav2 (:,:) = 0._wp 
         zvwdav2 (:,:) = 0._wp   
      END IF 

      !                                             ! ==================== !
      DO jn = 1, icycle                             !  sub-time-step loop  !
         !                                          ! ==================== !
         !
         l_full_nf_update = jn == icycle   ! false: disable full North fold update (performances) for jn = 1 to icycle-1
         !
         !                    !==  Update the forcing ==! (BDY and tides)
         !
         IF( ln_bdy      .AND. ln_tide )   CALL bdy_dta_tides( kt, kit=jn, kt_offset= noffset+1 )
         IF( ln_tide_pot .AND. ln_tide )   CALL upd_tide     ( kt, kit=jn, kt_offset= noffset   )
         !
         !                    !==  extrapolation at mid-step  ==!   (jn+1/2)
         !
         !                       !* Set extrapolation coefficients for predictor step:
         IF ((jn<3).AND.ll_init) THEN      ! Forward           
           za1 = 1._wp                                          
           za2 = 0._wp                        
           za3 = 0._wp                        
         ELSE                              ! AB3-AM4 Coefficients: bet=0.281105 
           za1 =  1.781105_wp              ! za1 =   3/2 +   bet
           za2 = -1.06221_wp               ! za2 = -(1/2 + 2*bet)
           za3 =  0.281105_wp              ! za3 = bet
         ENDIF
         !
         !                       !* Extrapolate barotropic velocities at mid-step (jn+1/2)
         !--        m+1/2               m                m-1           m-2       --!
         !--       u      = (3/2+beta) u   -(1/2+2beta) u      + beta u          --!
         !-------------------------------------------------------------------------!
         ua_e(:,:) = za1 * un_e(:,:) + za2 * ub_e(:,:) + za3 * ubb_e(:,:)
         va_e(:,:) = za1 * vn_e(:,:) + za2 * vb_e(:,:) + za3 * vbb_e(:,:)

         IF( .NOT.ln_linssh ) THEN                        !* Update ocean depth (variable volume case only)
            !                                             !  ------------------
            ! Extrapolate Sea Level at step jit+0.5:
            !--         m+1/2                 m                  m-1             m-2       --!
            !--      ssh      = (3/2+beta) ssh   -(1/2+2beta) ssh      + beta ssh          --!
            !--------------------------------------------------------------------------------!
            zsshp2_e(:,:) = za1 * sshn_e(:,:)  + za2 * sshb_e(:,:) + za3 * sshbb_e(:,:)
            
            ! set wetting & drying mask at tracer points for this barotropic mid-step
            IF( ln_wd_dl )   CALL wad_tmsk( zsshp2_e, ztwdmask )
            !
            !                          ! ocean t-depth at mid-step
            zhtp2_e(:,:) = ht_0(:,:) + zsshp2_e(:,:)
            !
            !                          ! ocean u- and v-depth at mid-step   (separate DO-loops remove the need of a lbc_lnk)
            DO jj = 1, jpj
               DO ji = 1, jpim1   ! not jpi-column
                  zhup2_e(ji,jj) = hu_0(ji,jj) + r1_2 * r1_e1e2u(ji,jj)                        &
                       &                              * (  e1e2t(ji  ,jj) * zsshp2_e(ji  ,jj)  &
                       &                                 + e1e2t(ji+1,jj) * zsshp2_e(ji+1,jj)  ) * ssumask(ji,jj)
               END DO
            END DO
            DO jj = 1, jpjm1        ! not jpj-row
               DO ji = 1, jpi
                  zhvp2_e(ji,jj) = hv_0(ji,jj) + r1_2 * r1_e1e2v(ji,jj)                        &
                       &                              * (  e1e2t(ji,jj  ) * zsshp2_e(ji,jj  )  &
                       &                                 + e1e2t(ji,jj+1) * zsshp2_e(ji,jj+1)  ) * ssvmask(ji,jj)
               END DO
            END DO
            !
         ENDIF
         !
         !                    !==  after SSH  ==!   (jn+1)
         !
         !                             ! update (ua_e,va_e) to enforce volume conservation at open boundaries
         !                             ! values of zhup2_e and zhvp2_e on the halo are not needed in bdy_vol2d
         IF( ln_bdy .AND. ln_vol ) CALL bdy_vol2d( kt, jn, ua_e, va_e, zhup2_e, zhvp2_e )
         !
         !                             ! resulting flux at mid-step (not over the full domain)
         zhU(1:jpim1,1:jpj  ) = e2u(1:jpim1,1:jpj  ) * ua_e(1:jpim1,1:jpj  ) * zhup2_e(1:jpim1,1:jpj  )   ! not jpi-column
         zhV(1:jpi  ,1:jpjm1) = e1v(1:jpi  ,1:jpjm1) * va_e(1:jpi  ,1:jpjm1) * zhvp2_e(1:jpi  ,1:jpjm1)   ! not jpj-row
         !
#if defined key_agrif
         ! Set fluxes during predictor step to ensure volume conservation
         IF( .NOT.Agrif_Root() .AND. ln_bt_fw ) THEN
            IF( l_Westedge ) THEN
               DO jj = 1, jpj
                  zhU(2:nbghostcells+1,jj) = ubdy_w(1:nbghostcells,jj) * e2u(2:nbghostcells+1,jj)
                  zhV(2:nbghostcells+1,jj) = vbdy_w(1:nbghostcells,jj) * e1v(2:nbghostcells+1,jj)
               END DO
            ENDIF
            IF( l_Eastedge ) THEN
               DO jj=1,jpj
                  zhU(nlci-nbghostcells-1:nlci-2,jj) = ubdy_e(1:nbghostcells,jj) * e2u(nlci-nbghostcells-1:nlci-2,jj)
                  zhV(nlci-nbghostcells  :nlci-1,jj) = vbdy_e(1:nbghostcells,jj) * e1v(nlci-nbghostcells  :nlci-1,jj)
               END DO
            ENDIF
            IF( l_Southedge ) THEN
               DO ji=1,jpi
                  zhV(ji,2:nbghostcells+1) = vbdy_s(ji,1:nbghostcells) * e1v(ji,2:nbghostcells+1)
                  zhU(ji,2:nbghostcells+1) = ubdy_s(ji,1:nbghostcells) * e2u(ji,2:nbghostcells+1)
               END DO
            ENDIF
            IF( l_Northedge ) THEN
               DO ji=1,jpi
                  zhV(ji,nlcj-nbghostcells-1:nlcj-2) = vbdy_n(ji,1:nbghostcells) * e1v(ji,nlcj-nbghostcells-1:nlcj-2)
                  zhU(ji,nlcj-nbghostcells  :nlcj-1) = ubdy_n(ji,1:nbghostcells) * e2u(ji,nlcj-nbghostcells  :nlcj-1)
               END DO
            ENDIF
         ENDIF
#endif
         IF( ln_wd_il )   CALL wad_lmt_bt(zhU, zhV, sshn_e, zssh_frc, rdtbt)    !!gm wad_lmt_bt use of lbc_lnk on zhU, zhV

         IF( ln_wd_dl ) THEN           ! un_e and vn_e are set to zero at faces where 
            !                          ! the direction of the flow is from dry cells
            CALL wad_Umsk( ztwdmask, zhU, zhV, un_e, vn_e, zuwdmask, zvwdmask )   ! not jpi colomn for U, not jpj row for V
            !
         ENDIF    
         !
         !
         !     Compute Sea Level at step jit+1
         !--           m+1        m                               m+1/2          --!
         !--        ssh    =  ssh   - delta_t' * [ frc + div( flux      ) ]      --!
         !-------------------------------------------------------------------------!
         DO jj = 2, jpjm1        ! INNER domain                             
            DO ji = 2, jpim1
               zhdiv = (   zhU(ji,jj) - zhU(ji-1,jj) + zhV(ji,jj) - zhV(ji,jj-1)   ) * r1_e1e2t(ji,jj)
               ssha_e(ji,jj) = (  sshn_e(ji,jj) - rdtbt * ( zssh_frc(ji,jj) + zhdiv )  ) * ssmask(ji,jj)
            END DO
         END DO
         !
         CALL lbc_lnk_multi( 'dynspg_ts', ssha_e, 'T', 1._wp,  zhU, 'U', -1._wp,  zhV, 'V', -1._wp )
         !
         !                             ! Sum over sub-time-steps to compute advective velocities
         za2 = wgtbtp2(jn)             ! zhU, zhV hold fluxes extrapolated at jn+0.5
         un_adv(:,:) = un_adv(:,:) + za2 * zhU(:,:) * r1_e2u(:,:)
         vn_adv(:,:) = vn_adv(:,:) + za2 * zhV(:,:) * r1_e1v(:,:)
         ! sum over sub-time-steps to decide which baroclinic velocities to set to zero (zuwdav2 is only used when ln_wd_dl_bc=True) 
         IF ( ln_wd_dl_bc ) THEN
            zuwdav2(1:jpim1,1:jpj  ) = zuwdav2(1:jpim1,1:jpj  ) + za2 * zuwdmask(1:jpim1,1:jpj  )   ! not jpi-column
            zvwdav2(1:jpi  ,1:jpjm1) = zvwdav2(1:jpi  ,1:jpjm1) + za2 * zvwdmask(1:jpi  ,1:jpjm1)   ! not jpj-row
         END IF
         !
         ! Duplicate sea level across open boundaries (this is only cosmetic if linssh=T)
         IF( ln_bdy )   CALL bdy_ssh( ssha_e )
#if defined key_agrif
         IF( .NOT.Agrif_Root() )   CALL agrif_ssh_ts( jn )
#endif
         !  
         ! Sea Surface Height at u-,v-points (vvl case only)
         IF( .NOT.ln_linssh ) THEN                                
            DO jj = 2, jpjm1   ! INNER domain, will be extended to whole domain later
               DO ji = 2, jpim1      ! NO Vector Opt.
                  zsshu_a(ji,jj) = r1_2 * ssumask(ji,jj) * r1_e1e2u(ji,jj)    &
                     &              * ( e1e2t(ji  ,jj  )  * ssha_e(ji  ,jj  ) &
                     &              +   e1e2t(ji+1,jj  )  * ssha_e(ji+1,jj  ) )
                  zsshv_a(ji,jj) = r1_2 * ssvmask(ji,jj) * r1_e1e2v(ji,jj)    &
                     &              * ( e1e2t(ji  ,jj  )  * ssha_e(ji  ,jj  ) &
                     &              +   e1e2t(ji  ,jj+1)  * ssha_e(ji  ,jj+1) )
               END DO
            END DO
         ENDIF   
         !         
         ! Half-step back interpolation of SSH for surface pressure computation at step jit+1/2
         !--            m+1/2           m+1              m               m-1              m-2     --!
         !--        ssh'    =  za0 * ssh     +  za1 * ssh   +  za2 * ssh      +  za3 * ssh        --!
         !------------------------------------------------------------------------------------------!
         CALL ts_bck_interp( jn, ll_init, za0, za1, za2, za3 )   ! coeficients of the interpolation
         zsshp2_e(:,:) = za0 *  ssha_e(:,:) + za1 *  sshn_e (:,:)   &
            &          + za2 *  sshb_e(:,:) + za3 *  sshbb_e(:,:)
         !
         !                             ! Surface pressure gradient
         zldg = ( 1._wp - rn_scal_load ) * grav    ! local factor
         DO jj = 2, jpjm1                            
            DO ji = 2, jpim1
               zu_spg(ji,jj) = - zldg * ( zsshp2_e(ji+1,jj) - zsshp2_e(ji,jj) ) * r1_e1u(ji,jj)
               zv_spg(ji,jj) = - zldg * ( zsshp2_e(ji,jj+1) - zsshp2_e(ji,jj) ) * r1_e2v(ji,jj)
            END DO
         END DO
         IF( ln_wd_il ) THEN        ! W/D : gravity filters applied on pressure gradient
            CALL wad_spg( zsshp2_e, zcpx, zcpy )   ! Calculating W/D gravity filters
            zu_spg(2:jpim1,2:jpjm1) = zu_spg(2:jpim1,2:jpjm1) * zcpx(2:jpim1,2:jpjm1)
            zv_spg(2:jpim1,2:jpjm1) = zv_spg(2:jpim1,2:jpjm1) * zcpy(2:jpim1,2:jpjm1)
         ENDIF
         !
         ! Add Coriolis trend:
         ! zwz array below or triads normally depend on sea level with ln_linssh=F and should be updated
         ! at each time step. We however keep them constant here for optimization.
         ! Recall that zhU and zhV hold fluxes at jn+0.5 (extrapolated not backward interpolated)
         CALL dyn_cor_2d( zhup2_e, zhvp2_e, ua_e, va_e, zhU, zhV,    zu_trd, zv_trd   )
         !
         ! Add tidal astronomical forcing if defined
         IF ( ln_tide .AND. ln_tide_pot ) THEN
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zu_trd(ji,jj) = zu_trd(ji,jj) + grav * ( pot_astro(ji+1,jj) - pot_astro(ji,jj) ) * r1_e1u(ji,jj)
                  zv_trd(ji,jj) = zv_trd(ji,jj) + grav * ( pot_astro(ji,jj+1) - pot_astro(ji,jj) ) * r1_e2v(ji,jj)
               END DO
            END DO
         ENDIF
         !
         ! Add bottom stresses:
!jth do implicitly instead
         IF ( .NOT. ll_wd ) THEN ! Revert to explicit for bit comparison tests in non wad runs
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zu_trd(ji,jj) = zu_trd(ji,jj) + zCdU_u(ji,jj) * un_e(ji,jj) * hur_e(ji,jj)
                  zv_trd(ji,jj) = zv_trd(ji,jj) + zCdU_v(ji,jj) * vn_e(ji,jj) * hvr_e(ji,jj)
               END DO
            END DO
         ENDIF
         !
         ! Set next velocities:
         !     Compute barotropic speeds at step jit+1    (h : total height of the water colomn)
         !--                              VECTOR FORM
         !--   m+1                 m               /                                                       m+1/2           \    --!
         !--  u     =             u   + delta_t' * \         (1-r)*g * grad_x( ssh') -         f * k vect u      +     frc /    --!
         !--                                                                                                                    --!
         !--                             FLUX FORM                                                                              --!
         !--  m+1   __1__  /  m    m               /  m+1/2                             m+1/2              m+1/2    n      \ \  --!
         !-- u    =   m+1 |  h  * u   + delta_t' * \ h     * (1-r)*g * grad_x( ssh') - h     * f * k vect u      + h * frc /  | --!
         !--         h     \                                                                                                 /  --!
         !------------------------------------------------------------------------------------------------------------------------!
         IF( ln_dynadv_vec .OR. ln_linssh ) THEN      !* Vector form
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ua_e(ji,jj) = (                                 un_e(ji,jj)   & 
                            &     + rdtbt * (                   zu_spg(ji,jj)   &
                            &                                 + zu_trd(ji,jj)   &
                            &                                 + zu_frc(ji,jj) ) & 
                            &   ) * ssumask(ji,jj)

                  va_e(ji,jj) = (                                 vn_e(ji,jj)   &
                            &     + rdtbt * (                   zv_spg(ji,jj)   &
                            &                                 + zv_trd(ji,jj)   &
                            &                                 + zv_frc(ji,jj) ) &
                            &   ) * ssvmask(ji,jj)
               END DO
            END DO
            !
         ELSE                           !* Flux form
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  !                    ! hu_e, hv_e hold depth at jn,  zhup2_e, zhvp2_e hold extrapolated depth at jn+1/2
                  !                    ! backward interpolated depth used in spg terms at jn+1/2
                  zhu_bck = hu_0(ji,jj) + r1_2*r1_e1e2u(ji,jj) * (  e1e2t(ji  ,jj) * zsshp2_e(ji  ,jj)    &
                       &                                          + e1e2t(ji+1,jj) * zsshp2_e(ji+1,jj)  ) * ssumask(ji,jj)
                  zhv_bck = hv_0(ji,jj) + r1_2*r1_e1e2v(ji,jj) * (  e1e2t(ji,jj  ) * zsshp2_e(ji,jj  )    &
                       &                                          + e1e2t(ji,jj+1) * zsshp2_e(ji,jj+1)  ) * ssvmask(ji,jj)
                  !                    ! inverse depth at jn+1
                  z1_hu = ssumask(ji,jj) / ( hu_0(ji,jj) + zsshu_a(ji,jj) + 1._wp - ssumask(ji,jj) )
                  z1_hv = ssvmask(ji,jj) / ( hv_0(ji,jj) + zsshv_a(ji,jj) + 1._wp - ssvmask(ji,jj) )
                  !
                  ua_e(ji,jj) = (               hu_e  (ji,jj) *   un_e (ji,jj)      & 
                       &            + rdtbt * (  zhu_bck        * zu_spg (ji,jj)  &   !
                       &                       + zhup2_e(ji,jj) * zu_trd (ji,jj)  &   !
                       &                       +  hu_n  (ji,jj) * zu_frc (ji,jj)  )   ) * z1_hu
                  !
                  va_e(ji,jj) = (               hv_e  (ji,jj) *   vn_e (ji,jj)      &
                       &            + rdtbt * (  zhv_bck        * zv_spg (ji,jj)  &   !
                       &                       + zhvp2_e(ji,jj) * zv_trd (ji,jj)  &   !
                       &                       +  hv_n  (ji,jj) * zv_frc (ji,jj)  )   ) * z1_hv
               END DO
            END DO
         ENDIF
!jth implicit bottom friction:
         IF ( ll_wd ) THEN ! revert to explicit for bit comparison tests in non wad runs
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                     ua_e(ji,jj) =  ua_e(ji,jj) /(1.0 -   rdtbt * zCdU_u(ji,jj) * hur_e(ji,jj))
                     va_e(ji,jj) =  va_e(ji,jj) /(1.0 -   rdtbt * zCdU_v(ji,jj) * hvr_e(ji,jj))
               END DO
            END DO
         ENDIF
       
         IF( .NOT.ln_linssh ) THEN   !* Update ocean depth (variable volume case only)
            hu_e (2:jpim1,2:jpjm1) = hu_0(2:jpim1,2:jpjm1) + zsshu_a(2:jpim1,2:jpjm1)
            hur_e(2:jpim1,2:jpjm1) = ssumask(2:jpim1,2:jpjm1) / ( hu_e(2:jpim1,2:jpjm1) + 1._wp - ssumask(2:jpim1,2:jpjm1) )
            hv_e (2:jpim1,2:jpjm1) = hv_0(2:jpim1,2:jpjm1) + zsshv_a(2:jpim1,2:jpjm1)
            hvr_e(2:jpim1,2:jpjm1) = ssvmask(2:jpim1,2:jpjm1) / ( hv_e(2:jpim1,2:jpjm1) + 1._wp - ssvmask(2:jpim1,2:jpjm1) )
            CALL lbc_lnk_multi( 'dynspg_ts', ua_e , 'U', -1._wp, va_e , 'V', -1._wp  &
                 &                         , hu_e , 'U',  1._wp, hv_e , 'V',  1._wp  &
                 &                         , hur_e, 'U',  1._wp, hvr_e, 'V',  1._wp  )
         ELSE
            CALL lbc_lnk_multi( 'dynspg_ts', ua_e , 'U', -1._wp, va_e , 'V', -1._wp  )
         ENDIF
         !
         !
         !                                                 ! open boundaries
         IF( ln_bdy )   CALL bdy_dyn2d( jn, ua_e, va_e, un_e, vn_e, hur_e, hvr_e, ssha_e )
#if defined key_agrif                                                           
         IF( .NOT.Agrif_Root() )  CALL agrif_dyn_ts( jn )  ! Agrif
#endif
         !                                             !* Swap
         !                                             !  ----
         ubb_e  (:,:) = ub_e  (:,:)
         ub_e   (:,:) = un_e  (:,:)
         un_e   (:,:) = ua_e  (:,:)
         !
         vbb_e  (:,:) = vb_e  (:,:)
         vb_e   (:,:) = vn_e  (:,:)
         vn_e   (:,:) = va_e  (:,:)
         !
         sshbb_e(:,:) = sshb_e(:,:)
         sshb_e (:,:) = sshn_e(:,:)
         sshn_e (:,:) = ssha_e(:,:)

         !                                             !* Sum over whole bt loop
         !                                             !  ----------------------
         za1 = wgtbtp1(jn)                                    
         IF( ln_dynadv_vec .OR. ln_linssh ) THEN    ! Sum velocities
            ua_b  (:,:) = ua_b  (:,:) + za1 * ua_e  (:,:) 
            va_b  (:,:) = va_b  (:,:) + za1 * va_e  (:,:) 
         ELSE                                       ! Sum transports
            IF ( .NOT.ln_wd_dl ) THEN  
               ua_b  (:,:) = ua_b  (:,:) + za1 * ua_e  (:,:) * hu_e (:,:)
               va_b  (:,:) = va_b  (:,:) + za1 * va_e  (:,:) * hv_e (:,:)
            ELSE 
               ua_b  (:,:) = ua_b  (:,:) + za1 * ua_e  (:,:) * hu_e (:,:) * zuwdmask(:,:)
               va_b  (:,:) = va_b  (:,:) + za1 * va_e  (:,:) * hv_e (:,:) * zvwdmask(:,:)
            END IF 
         ENDIF
         !                                          ! Sum sea level
         ssha(:,:) = ssha(:,:) + za1 * ssha_e(:,:)

         !                                                 ! ==================== !
      END DO                                               !        end loop      !
      !                                                    ! ==================== !
      ! -----------------------------------------------------------------------------
      ! Phase 3. update the general trend with the barotropic trend
      ! -----------------------------------------------------------------------------
      !
      ! Set advection velocity correction:
      IF (ln_bt_fw) THEN
         IF( .NOT.( kt == nit000 .AND. neuler==0 ) ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zun_save = un_adv(ji,jj)
                  zvn_save = vn_adv(ji,jj)
                  !                          ! apply the previously computed correction 
                  un_adv(ji,jj) = r1_2 * ( ub2_b(ji,jj) + zun_save - atfp * un_bf(ji,jj) )
                  vn_adv(ji,jj) = r1_2 * ( vb2_b(ji,jj) + zvn_save - atfp * vn_bf(ji,jj) )
                  !                          ! Update corrective fluxes for next time step
                  un_bf(ji,jj)  = atfp * un_bf(ji,jj) + ( zun_save - ub2_b(ji,jj) )
                  vn_bf(ji,jj)  = atfp * vn_bf(ji,jj) + ( zvn_save - vb2_b(ji,jj) )
                  !                          ! Save integrated transport for next computation
                  ub2_b(ji,jj) = zun_save
                  vb2_b(ji,jj) = zvn_save
               END DO
            END DO
         ELSE
            un_bf(:,:) = 0._wp            ! corrective fluxes for next time step set to zero
            vn_bf(:,:) = 0._wp
            ub2_b(:,:) = un_adv(:,:)      ! Save integrated transport for next computation
            vb2_b(:,:) = vn_adv(:,:)
         END IF
      ENDIF


      !
      ! Update barotropic trend:
      IF( ln_dynadv_vec .OR. ln_linssh ) THEN
         DO jk=1,jpkm1
            ua(:,:,jk) = ua(:,:,jk) + ( ua_b(:,:) - ub_b(:,:) ) * r1_2dt_b
            va(:,:,jk) = va(:,:,jk) + ( va_b(:,:) - vb_b(:,:) ) * r1_2dt_b
         END DO
      ELSE
         ! At this stage, ssha has been corrected: compute new depths at velocity points
         DO jj = 1, jpjm1
            DO ji = 1, jpim1      ! NO Vector Opt.
               zsshu_a(ji,jj) = r1_2 * ssumask(ji,jj)  * r1_e1e2u(ji,jj) &
                  &              * ( e1e2t(ji  ,jj) * ssha(ji  ,jj)      &
                  &              +   e1e2t(ji+1,jj) * ssha(ji+1,jj) )
               zsshv_a(ji,jj) = r1_2 * ssvmask(ji,jj)  * r1_e1e2v(ji,jj) &
                  &              * ( e1e2t(ji,jj  ) * ssha(ji,jj  )      &
                  &              +   e1e2t(ji,jj+1) * ssha(ji,jj+1) )
            END DO
         END DO
         CALL lbc_lnk_multi( 'dynspg_ts', zsshu_a, 'U', 1._wp, zsshv_a, 'V', 1._wp ) ! Boundary conditions
         !
         DO jk=1,jpkm1
            ua(:,:,jk) = ua(:,:,jk) + r1_hu_n(:,:) * ( ua_b(:,:) - ub_b(:,:) * hu_b(:,:) ) * r1_2dt_b
            va(:,:,jk) = va(:,:,jk) + r1_hv_n(:,:) * ( va_b(:,:) - vb_b(:,:) * hv_b(:,:) ) * r1_2dt_b
         END DO
         ! Save barotropic velocities not transport:
         ua_b(:,:) =  ua_b(:,:) / ( hu_0(:,:) + zsshu_a(:,:) + 1._wp - ssumask(:,:) )
         va_b(:,:) =  va_b(:,:) / ( hv_0(:,:) + zsshv_a(:,:) + 1._wp - ssvmask(:,:) )
      ENDIF


      ! Correct velocities so that the barotropic velocity equals (un_adv, vn_adv) (in all cases)  
      DO jk = 1, jpkm1
         un(:,:,jk) = ( un(:,:,jk) + un_adv(:,:)*r1_hu_n(:,:) - un_b(:,:) ) * umask(:,:,jk)
         vn(:,:,jk) = ( vn(:,:,jk) + vn_adv(:,:)*r1_hv_n(:,:) - vn_b(:,:) ) * vmask(:,:,jk)
      END DO

      IF ( ln_wd_dl .and. ln_wd_dl_bc) THEN 
         ! need to set lbc here because not done prior time averaging
         CALL lbc_lnk_multi( 'dynspg_ts', zuwdav2, 'U', 1._wp, zvwdav2, 'V', 1._wp)
         DO jk = 1, jpkm1
            un(:,:,jk) = ( un_adv(:,:)*r1_hu_n(:,:) &
                       & + zuwdav2(:,:)*(un(:,:,jk) - un_adv(:,:)*r1_hu_n(:,:)) ) * umask(:,:,jk) 
            vn(:,:,jk) = ( vn_adv(:,:)*r1_hv_n(:,:) & 
                       & + zvwdav2(:,:)*(vn(:,:,jk) - vn_adv(:,:)*r1_hv_n(:,:)) ) * vmask(:,:,jk)  
         END DO
      END IF 

      
      CALL iom_put(  "ubar", un_adv(:,:)*r1_hu_n(:,:) )    ! barotropic i-current
      CALL iom_put(  "vbar", vn_adv(:,:)*r1_hv_n(:,:) )    ! barotropic i-current
      !
#if defined key_agrif
      ! Save time integrated fluxes during child grid integration
      ! (used to update coarse grid transports at next time step)
      !
      IF( .NOT.Agrif_Root() .AND. ln_bt_fw ) THEN
         IF( Agrif_NbStepint() == 0 ) THEN
            ub2_i_b(:,:) = 0._wp
            vb2_i_b(:,:) = 0._wp
         END IF
         !
         za1 = 1._wp / REAL(Agrif_rhot(), wp)
         ub2_i_b(:,:) = ub2_i_b(:,:) + za1 * ub2_b(:,:)
         vb2_i_b(:,:) = vb2_i_b(:,:) + za1 * vb2_b(:,:)
      ENDIF
#endif      
      !                                   !* write time-spliting arrays in the restart
      IF( lrst_oce .AND.ln_bt_fw )   CALL ts_rst( kt, 'WRITE' )
      !
      IF( ln_wd_il )   DEALLOCATE( zcpx, zcpy )
      IF( ln_wd_dl )   DEALLOCATE( ztwdmask, zuwdmask, zvwdmask, zuwdav2, zvwdav2 )
      !
      CALL iom_put( "baro_u" , un_b )  ! Barotropic  U Velocity
      CALL iom_put( "baro_v" , vn_b )  ! Barotropic  V Velocity
      !
   END SUBROUTINE dyn_spg_ts


   SUBROUTINE ts_wgt( ll_av, ll_fw, jpit, zwgt1, zwgt2)
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ts_wgt  ***
      !!
      !! ** Purpose : Set time-splitting weights for temporal averaging (or not)
      !!----------------------------------------------------------------------
      LOGICAL, INTENT(in) ::   ll_av      ! temporal averaging=.true.
      LOGICAL, INTENT(in) ::   ll_fw      ! forward time splitting =.true.
      INTEGER, INTENT(inout) :: jpit      ! cycle length    
      REAL(wp), DIMENSION(3*nn_baro), INTENT(inout) ::   zwgt1, & ! Primary weights
                                                         zwgt2    ! Secondary weights
      
      INTEGER ::  jic, jn, ji                      ! temporary integers
      REAL(wp) :: za1, za2
      !!----------------------------------------------------------------------

      zwgt1(:) = 0._wp
      zwgt2(:) = 0._wp

      ! Set time index when averaged value is requested
      IF (ll_fw) THEN 
         jic = nn_baro
      ELSE
         jic = 2 * nn_baro
      ENDIF

      ! Set primary weights:
      IF (ll_av) THEN
           ! Define simple boxcar window for primary weights 
           ! (width = nn_baro, centered around jic)     
         SELECT CASE ( nn_bt_flt )
              CASE( 0 )  ! No averaging
                 zwgt1(jic) = 1._wp
                 jpit = jic

              CASE( 1 )  ! Boxcar, width = nn_baro
                 DO jn = 1, 3*nn_baro
                    za1 = ABS(float(jn-jic))/float(nn_baro) 
                    IF (za1 < 0.5_wp) THEN
                      zwgt1(jn) = 1._wp
                      jpit = jn
                    ENDIF
                 ENDDO

              CASE( 2 )  ! Boxcar, width = 2 * nn_baro
                 DO jn = 1, 3*nn_baro
                    za1 = ABS(float(jn-jic))/float(nn_baro) 
                    IF (za1 < 1._wp) THEN
                      zwgt1(jn) = 1._wp
                      jpit = jn
                    ENDIF
                 ENDDO
              CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for nn_bt_flt' )
         END SELECT

      ELSE ! No time averaging
         zwgt1(jic) = 1._wp
         jpit = jic
      ENDIF
    
      ! Set secondary weights
      DO jn = 1, jpit
        DO ji = jn, jpit
             zwgt2(jn) = zwgt2(jn) + zwgt1(ji)
        END DO
      END DO

      ! Normalize weigths:
      za1 = 1._wp / SUM(zwgt1(1:jpit))
      za2 = 1._wp / SUM(zwgt2(1:jpit))
      DO jn = 1, jpit
        zwgt1(jn) = zwgt1(jn) * za1
        zwgt2(jn) = zwgt2(jn) * za2
      END DO
      !
   END SUBROUTINE ts_wgt


   SUBROUTINE ts_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ts_rst  ***
      !!
      !! ** Purpose : Read or write time-splitting arrays in restart file
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt     ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
         !                                   ! ---------------
         IF( ln_rstart .AND. ln_bt_fw .AND. (neuler/=0) ) THEN    !* Read the restart file
            CALL iom_get( numror, jpdom_autoglo, 'ub2_b'  , ub2_b  (:,:), ldxios = lrxios )   
            CALL iom_get( numror, jpdom_autoglo, 'vb2_b'  , vb2_b  (:,:), ldxios = lrxios ) 
            CALL iom_get( numror, jpdom_autoglo, 'un_bf'  , un_bf  (:,:), ldxios = lrxios )   
            CALL iom_get( numror, jpdom_autoglo, 'vn_bf'  , vn_bf  (:,:), ldxios = lrxios ) 
            IF( .NOT.ln_bt_av ) THEN
               CALL iom_get( numror, jpdom_autoglo, 'sshbb_e'  , sshbb_e(:,:), ldxios = lrxios )   
               CALL iom_get( numror, jpdom_autoglo, 'ubb_e'    ,   ubb_e(:,:), ldxios = lrxios )   
               CALL iom_get( numror, jpdom_autoglo, 'vbb_e'    ,   vbb_e(:,:), ldxios = lrxios )
               CALL iom_get( numror, jpdom_autoglo, 'sshb_e'   ,  sshb_e(:,:), ldxios = lrxios ) 
               CALL iom_get( numror, jpdom_autoglo, 'ub_e'     ,    ub_e(:,:), ldxios = lrxios )   
               CALL iom_get( numror, jpdom_autoglo, 'vb_e'     ,    vb_e(:,:), ldxios = lrxios )
            ENDIF
#if defined key_agrif
            ! Read time integrated fluxes
            IF ( .NOT.Agrif_Root() ) THEN
               CALL iom_get( numror, jpdom_autoglo, 'ub2_i_b'  , ub2_i_b(:,:), ldxios = lrxios )   
               CALL iom_get( numror, jpdom_autoglo, 'vb2_i_b'  , vb2_i_b(:,:), ldxios = lrxios )
            ENDIF
#endif
         ELSE                                   !* Start from rest
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set barotropic values to 0'
            ub2_b (:,:) = 0._wp   ;   vb2_b (:,:) = 0._wp   ! used in the 1st interpol of agrif
            un_adv(:,:) = 0._wp   ;   vn_adv(:,:) = 0._wp   ! used in the 1st interpol of agrif
            un_bf (:,:) = 0._wp   ;   vn_bf (:,:) = 0._wp   ! used in the 1st update   of agrif
#if defined key_agrif
            IF ( .NOT.Agrif_Root() ) THEN
               ub2_i_b(:,:) = 0._wp   ;   vb2_i_b(:,:) = 0._wp   ! used in the 1st update of agrif
            ENDIF
#endif
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*) '---- ts_rst ----'
         IF( lwxios ) CALL iom_swap(      cwxios_context          )
         CALL iom_rstput( kt, nitrst, numrow, 'ub2_b'   , ub2_b  (:,:), ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'vb2_b'   , vb2_b  (:,:), ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'un_bf'   , un_bf  (:,:), ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'vn_bf'   , vn_bf  (:,:), ldxios = lwxios )
         !
         IF (.NOT.ln_bt_av) THEN
            CALL iom_rstput( kt, nitrst, numrow, 'sshbb_e'  , sshbb_e(:,:), ldxios = lwxios ) 
            CALL iom_rstput( kt, nitrst, numrow, 'ubb_e'    ,   ubb_e(:,:), ldxios = lwxios )
            CALL iom_rstput( kt, nitrst, numrow, 'vbb_e'    ,   vbb_e(:,:), ldxios = lwxios )
            CALL iom_rstput( kt, nitrst, numrow, 'sshb_e'   ,  sshb_e(:,:), ldxios = lwxios )
            CALL iom_rstput( kt, nitrst, numrow, 'ub_e'     ,    ub_e(:,:), ldxios = lwxios )
            CALL iom_rstput( kt, nitrst, numrow, 'vb_e'     ,    vb_e(:,:), ldxios = lwxios )
         ENDIF
#if defined key_agrif
         ! Save time integrated fluxes
         IF ( .NOT.Agrif_Root() ) THEN
            CALL iom_rstput( kt, nitrst, numrow, 'ub2_i_b'  , ub2_i_b(:,:), ldxios = lwxios )
            CALL iom_rstput( kt, nitrst, numrow, 'vb2_i_b'  , vb2_i_b(:,:), ldxios = lwxios )
         ENDIF
#endif
         IF( lwxios ) CALL iom_swap(      cxios_context          )
      ENDIF
      !
   END SUBROUTINE ts_rst


   SUBROUTINE dyn_spg_ts_init
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_spg_ts_init  ***
      !!
      !! ** Purpose : Set time splitting options
      !!----------------------------------------------------------------------
      INTEGER  ::   ji ,jj              ! dummy loop indices
      REAL(wp) ::   zxr2, zyr2, zcmax   ! local scalar
      REAL(wp), DIMENSION(jpi,jpj) ::   zcu
      INTEGER  :: inum
      !!----------------------------------------------------------------------
      !
      ! Max courant number for ext. grav. waves
      !
      DO jj = 1, jpj
         DO ji =1, jpi
            zxr2 = r1_e1t(ji,jj) * r1_e1t(ji,jj)
            zyr2 = r1_e2t(ji,jj) * r1_e2t(ji,jj)
            zcu(ji,jj) = SQRT( grav * MAX(ht_0(ji,jj),0._wp) * (zxr2 + zyr2) )
         END DO
      END DO
      !
      zcmax = MAXVAL( zcu(:,:) )
      CALL mpp_max( 'dynspg_ts', zcmax )

      ! Estimate number of iterations to satisfy a max courant number= rn_bt_cmax
      IF( ln_bt_auto )   nn_baro = CEILING( rdt / rn_bt_cmax * zcmax)
      
      rdtbt = rdt / REAL( nn_baro , wp )
      zcmax = zcmax * rdtbt
      ! Print results
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dyn_spg_ts_init : split-explicit free surface'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      IF( ln_bt_auto ) THEN
         IF(lwp) WRITE(numout,*) '     ln_ts_auto =.true. Automatically set nn_baro '
         IF(lwp) WRITE(numout,*) '     Max. courant number allowed: ', rn_bt_cmax
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_ts_auto=.false.: Use nn_baro in namelist   nn_baro = ', nn_baro
      ENDIF

      IF(ln_bt_av) THEN
         IF(lwp) WRITE(numout,*) '     ln_bt_av =.true.  ==> Time averaging over nn_baro time steps is on '
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_bt_av =.false. => No time averaging of barotropic variables '
      ENDIF
      !
      !
      IF(ln_bt_fw) THEN
         IF(lwp) WRITE(numout,*) '     ln_bt_fw=.true.  => Forward integration of barotropic variables '
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_bt_fw =.false.=> Centred integration of barotropic variables '
      ENDIF
      !
#if defined key_agrif
      ! Restrict the use of Agrif to the forward case only
!!!      IF( .NOT.ln_bt_fw .AND. .NOT.Agrif_Root() )   CALL ctl_stop( 'AGRIF not implemented if ln_bt_fw=.FALSE.' )
#endif
      !
      IF(lwp) WRITE(numout,*)    '     Time filter choice, nn_bt_flt: ', nn_bt_flt
      SELECT CASE ( nn_bt_flt )
         CASE( 0 )      ;   IF(lwp) WRITE(numout,*) '           Dirac'
         CASE( 1 )      ;   IF(lwp) WRITE(numout,*) '           Boxcar: width = nn_baro'
         CASE( 2 )      ;   IF(lwp) WRITE(numout,*) '           Boxcar: width = 2*nn_baro' 
         CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for nn_bt_flt: should 0,1, or 2' )
      END SELECT
      !
      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) '     nn_baro = ', nn_baro
      IF(lwp) WRITE(numout,*) '     Barotropic time step [s] is :', rdtbt
      IF(lwp) WRITE(numout,*) '     Maximum Courant number is   :', zcmax
      !
      IF(lwp) WRITE(numout,*)    '     Time diffusion parameter rn_bt_alpha: ', rn_bt_alpha
      IF ((ln_bt_av.AND.nn_bt_flt/=0).AND.(rn_bt_alpha>0._wp)) THEN
         CALL ctl_stop( 'dynspg_ts ERROR: if rn_bt_alpha > 0, remove temporal averaging' )
      ENDIF
      !
      IF( .NOT.ln_bt_av .AND. .NOT.ln_bt_fw ) THEN
         CALL ctl_stop( 'dynspg_ts ERROR: No time averaging => only forward integration is possible' )
      ENDIF
      IF( zcmax>0.9_wp ) THEN
         CALL ctl_stop( 'dynspg_ts ERROR: Maximum Courant number is greater than 0.9: Inc. nn_baro !' )          
      ENDIF
      !
      !                             ! Allocate time-splitting arrays
      IF( dyn_spg_ts_alloc() /= 0    )   CALL ctl_stop('STOP', 'dyn_spg_init: failed to allocate dynspg_ts  arrays' )
      !
      !                             ! read restart when needed
      CALL ts_rst( nit000, 'READ' )
      !
      IF( lwxios ) THEN
! define variables in restart file when writing with XIOS
         CALL iom_set_rstw_var_active('ub2_b')
         CALL iom_set_rstw_var_active('vb2_b')
         CALL iom_set_rstw_var_active('un_bf')
         CALL iom_set_rstw_var_active('vn_bf')
         !
         IF (.NOT.ln_bt_av) THEN
            CALL iom_set_rstw_var_active('sshbb_e')
            CALL iom_set_rstw_var_active('ubb_e')
            CALL iom_set_rstw_var_active('vbb_e')
            CALL iom_set_rstw_var_active('sshb_e')
            CALL iom_set_rstw_var_active('ub_e')
            CALL iom_set_rstw_var_active('vb_e')
         ENDIF
#if defined key_agrif
         ! Save time integrated fluxes
         IF ( .NOT.Agrif_Root() ) THEN
            CALL iom_set_rstw_var_active('ub2_i_b')
            CALL iom_set_rstw_var_active('vb2_i_b')
         ENDIF
#endif
      ENDIF
      !
   END SUBROUTINE dyn_spg_ts_init

   
   SUBROUTINE dyn_cor_2d_init
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_cor_2d_init  ***
      !!
      !! ** Purpose : Set time splitting options
      !! Set arrays to remove/compute coriolis trend.
      !! Do it once during initialization if volume is fixed, else at each long time step.
      !! Note that these arrays are also used during barotropic loop. These are however frozen
      !! although they should be updated in the variable volume case. Not a big approximation.
      !! To remove this approximation, copy lines below inside barotropic loop
      !! and update depths at T-F points (ht and zhf resp.) at each barotropic time step
      !!
      !! Compute zwz = f / ( height of the water colomn )
      !!----------------------------------------------------------------------
      INTEGER  ::   ji ,jj, jk              ! dummy loop indices
      REAL(wp) ::   z1_ht
      REAL(wp), DIMENSION(jpi,jpj) :: zhf
      !!----------------------------------------------------------------------
      !
      SELECT CASE( nvor_scheme )
      CASE( np_EEN )                != EEN scheme using e3f (energy & enstrophy scheme)
         SELECT CASE( nn_een_e3f )              !* ff_f/e3 at F-point
         CASE ( 0 )                                   ! original formulation  (masked averaging of e3t divided by 4)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zwz(ji,jj) =   ( ht_n(ji  ,jj+1) + ht_n(ji+1,jj+1) +                    &
                       &             ht_n(ji  ,jj  ) + ht_n(ji+1,jj  )   ) * 0.25_wp  
                  IF( zwz(ji,jj) /= 0._wp )   zwz(ji,jj) = ff_f(ji,jj) / zwz(ji,jj)
               END DO
            END DO
         CASE ( 1 )                                   ! new formulation  (masked averaging of e3t divided by the sum of mask)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zwz(ji,jj) =             (  ht_n  (ji  ,jj+1) + ht_n  (ji+1,jj+1)      &
                       &                      + ht_n  (ji  ,jj  ) + ht_n  (ji+1,jj  )  )   &
                       &       / ( MAX( 1._wp,  ssmask(ji  ,jj+1) + ssmask(ji+1,jj+1)      &
                       &                      + ssmask(ji  ,jj  ) + ssmask(ji+1,jj  )  )   )
                  IF( zwz(ji,jj) /= 0._wp )   zwz(ji,jj) = ff_f(ji,jj) / zwz(ji,jj)
               END DO
            END DO
         END SELECT
         CALL lbc_lnk( 'dynspg_ts', zwz, 'F', 1._wp )
         !
         ftne(1,:) = 0._wp ; ftnw(1,:) = 0._wp ; ftse(1,:) = 0._wp ; ftsw(1,:) = 0._wp
         DO jj = 2, jpj
            DO ji = 2, jpi
               ftne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
               ftnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
               ftse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
               ftsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
            END DO
         END DO
         !
      CASE( np_EET )                  != EEN scheme using e3t (energy conserving scheme)
         ftne(1,:) = 0._wp ; ftnw(1,:) = 0._wp ; ftse(1,:) = 0._wp ; ftsw(1,:) = 0._wp
         DO jj = 2, jpj
            DO ji = 2, jpi
               z1_ht = ssmask(ji,jj) / ( ht_n(ji,jj) + 1._wp - ssmask(ji,jj) )
               ftne(ji,jj) = ( ff_f(ji-1,jj  ) + ff_f(ji  ,jj  ) + ff_f(ji  ,jj-1) ) * z1_ht
               ftnw(ji,jj) = ( ff_f(ji-1,jj-1) + ff_f(ji-1,jj  ) + ff_f(ji  ,jj  ) ) * z1_ht
               ftse(ji,jj) = ( ff_f(ji  ,jj  ) + ff_f(ji  ,jj-1) + ff_f(ji-1,jj-1) ) * z1_ht
               ftsw(ji,jj) = ( ff_f(ji  ,jj-1) + ff_f(ji-1,jj-1) + ff_f(ji-1,jj  ) ) * z1_ht
            END DO
         END DO
         !
      CASE( np_ENE, np_ENS , np_MIX )  != all other schemes (ENE, ENS, MIX) except ENT !
         !
         zwz(:,:) = 0._wp
         zhf(:,:) = 0._wp
         
         !!gm  assume 0 in both cases (which is almost surely WRONG ! ) as hvatf has been removed 
!!gm    A priori a better value should be something like :
!!gm          zhf(i,j) = masked sum of  ht(i,j) , ht(i+1,j) , ht(i,j+1) , (i+1,j+1) 
!!gm                     divided by the sum of the corresponding mask 
!!gm 
!!            
         IF( .NOT.ln_sco ) THEN
  
   !!gm  agree the JC comment  : this should be done in a much clear way
  
   ! JC: It not clear yet what should be the depth at f-points over land in z-coordinate case
   !     Set it to zero for the time being 
   !              IF( rn_hmin < 0._wp ) THEN    ;   jk = - INT( rn_hmin )                                      ! from a nb of level
   !              ELSE                          ;   jk = MINLOC( gdepw_0, mask = gdepw_0 > rn_hmin, dim = 1 )  ! from a depth
   !              ENDIF
   !              zhf(:,:) = gdepw_0(:,:,jk+1)
            !
         ELSE
            !
            !zhf(:,:) = hbatf(:,:)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zhf(ji,jj) =    (   ht_0  (ji,jj  ) + ht_0  (ji+1,jj  )          &
                       &            + ht_0  (ji,jj+1) + ht_0  (ji+1,jj+1)   )      &
                       &     / MAX(   ssmask(ji,jj  ) + ssmask(ji+1,jj  )          &
                       &            + ssmask(ji,jj+1) + ssmask(ji+1,jj+1) , 1._wp  )
               END DO
            END DO
         ENDIF
         !
         DO jj = 1, jpjm1
            zhf(:,jj) = zhf(:,jj) * (1._wp- umask(:,jj,1) * umask(:,jj+1,1))
         END DO
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               zhf(:,jj) = zhf(:,jj) + e3f_n(:,jj,jk) * umask(:,jj,jk) * umask(:,jj+1,jk)
            END DO
         END DO
         CALL lbc_lnk( 'dynspg_ts', zhf, 'F', 1._wp )
         ! JC: TBC. hf should be greater than 0 
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( zhf(ji,jj) /= 0._wp )   zwz(ji,jj) = 1._wp / zhf(ji,jj)
            END DO
         END DO
         zwz(:,:) = ff_f(:,:) * zwz(:,:)
      END SELECT
      
   END SUBROUTINE dyn_cor_2d_init



   SUBROUTINE dyn_cor_2d( hu_n, hv_n, un_b, vn_b, zhU, zhV,    zu_trd, zv_trd   )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_cor_2d  ***
      !!
      !! ** Purpose : Compute u and v coriolis trends
      !!----------------------------------------------------------------------
      INTEGER  ::   ji ,jj                             ! dummy loop indices
      REAL(wp) ::   zx1, zx2, zy1, zy2, z1_hu, z1_hv   !   -      -
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: hu_n, hv_n, un_b, vn_b, zhU, zhV
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) :: zu_trd, zv_trd
      !!----------------------------------------------------------------------
      SELECT CASE( nvor_scheme )
      CASE( np_ENT )                ! enstrophy conserving scheme (f-point)
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               z1_hu = ssumask(ji,jj) / ( hu_n(ji,jj) + 1._wp - ssumask(ji,jj) )
               z1_hv = ssvmask(ji,jj) / ( hv_n(ji,jj) + 1._wp - ssvmask(ji,jj) )
               zu_trd(ji,jj) = + r1_4 * r1_e1e2u(ji,jj) * z1_hu                    &
                  &               * (  e1e2t(ji+1,jj)*ht_n(ji+1,jj)*ff_t(ji+1,jj) * ( vn_b(ji+1,jj) + vn_b(ji+1,jj-1) )   &
                  &                  + e1e2t(ji  ,jj)*ht_n(ji  ,jj)*ff_t(ji  ,jj) * ( vn_b(ji  ,jj) + vn_b(ji  ,jj-1) )   )
                  !
               zv_trd(ji,jj) = - r1_4 * r1_e1e2v(ji,jj) * z1_hv                    &
                  &               * (  e1e2t(ji,jj+1)*ht_n(ji,jj+1)*ff_t(ji,jj+1) * ( un_b(ji,jj+1) + un_b(ji-1,jj+1) )   & 
                  &                  + e1e2t(ji,jj  )*ht_n(ji,jj  )*ff_t(ji,jj  ) * ( un_b(ji,jj  ) + un_b(ji-1,jj  ) )   ) 
            END DO  
         END DO  
         !         
      CASE( np_ENE , np_MIX )        ! energy conserving scheme (t-point) ENE or MIX
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zy1 = ( zhV(ji,jj-1) + zhV(ji+1,jj-1) ) * r1_e1u(ji,jj)
               zy2 = ( zhV(ji,jj  ) + zhV(ji+1,jj  ) ) * r1_e1u(ji,jj)
               zx1 = ( zhU(ji-1,jj) + zhU(ji-1,jj+1) ) * r1_e2v(ji,jj)
               zx2 = ( zhU(ji  ,jj) + zhU(ji  ,jj+1) ) * r1_e2v(ji,jj)
               ! energy conserving formulation for planetary vorticity term
               zu_trd(ji,jj) =   r1_4 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
               zv_trd(ji,jj) = - r1_4 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
            END DO
         END DO
         !
      CASE( np_ENS )                ! enstrophy conserving scheme (f-point)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zy1 =   r1_8 * ( zhV(ji  ,jj-1) + zhV(ji+1,jj-1) &
                 &            + zhV(ji  ,jj  ) + zhV(ji+1,jj  ) ) * r1_e1u(ji,jj)
               zx1 = - r1_8 * ( zhU(ji-1,jj  ) + zhU(ji-1,jj+1) &
                 &            + zhU(ji  ,jj  ) + zhU(ji  ,jj+1) ) * r1_e2v(ji,jj)
               zu_trd(ji,jj)  = zy1 * ( zwz(ji  ,jj-1) + zwz(ji,jj) )
               zv_trd(ji,jj)  = zx1 * ( zwz(ji-1,jj  ) + zwz(ji,jj) )
            END DO
         END DO
         !
      CASE( np_EET , np_EEN )      ! energy & enstrophy scheme (using e3t or e3f)         
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zu_trd(ji,jj) = + r1_12 * r1_e1u(ji,jj) * (  ftne(ji,jj  ) * zhV(ji  ,jj  ) &
                &                                         + ftnw(ji+1,jj) * zhV(ji+1,jj  ) &
                &                                         + ftse(ji,jj  ) * zhV(ji  ,jj-1) &
                &                                         + ftsw(ji+1,jj) * zhV(ji+1,jj-1) )
               zv_trd(ji,jj) = - r1_12 * r1_e2v(ji,jj) * (  ftsw(ji,jj+1) * zhU(ji-1,jj+1) &
                &                                         + ftse(ji,jj+1) * zhU(ji  ,jj+1) &
                &                                         + ftnw(ji,jj  ) * zhU(ji-1,jj  ) &
                &                                         + ftne(ji,jj  ) * zhU(ji  ,jj  ) )
            END DO
         END DO
         !
      END SELECT
      !
   END SUBROUTINE dyn_cor_2D


   SUBROUTINE wad_tmsk( pssh, ptmsk )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE wad_lmt  ***
      !!                    
      !! ** Purpose :   set wetting & drying mask at tracer points 
      !!              for the current barotropic sub-step 
      !!
      !! ** Method  :   ??? 
      !!
      !! ** Action  :  ptmsk : wetting & drying t-mask
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pssh    !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   ptmsk   !
      !
      INTEGER  ::   ji, jj   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( ln_wd_dl_rmp ) THEN     
         DO jj = 1, jpj
            DO ji = 1, jpi                    
               IF    ( pssh(ji,jj) + ht_0(ji,jj) >  2._wp * rn_wdmin1 ) THEN 
                  !           IF    ( pssh(ji,jj) + ht_0(ji,jj) >          rn_wdmin2 ) THEN 
                  ptmsk(ji,jj) = 1._wp
               ELSEIF( pssh(ji,jj) + ht_0(ji,jj) >          rn_wdmin1 ) THEN
                  ptmsk(ji,jj) = TANH( 50._wp*( ( pssh(ji,jj) + ht_0(ji,jj) -  rn_wdmin1 )*r_rn_wdmin1) )
               ELSE 
                  ptmsk(ji,jj) = 0._wp
               ENDIF
            END DO
         END DO
      ELSE  
         DO jj = 1, jpj
            DO ji = 1, jpi                              
               IF ( pssh(ji,jj) + ht_0(ji,jj) >  rn_wdmin1 ) THEN   ;   ptmsk(ji,jj) = 1._wp
               ELSE                                                 ;   ptmsk(ji,jj) = 0._wp
               ENDIF
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE wad_tmsk


   SUBROUTINE wad_Umsk( pTmsk, phU, phV, pu, pv, pUmsk, pVmsk )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE wad_lmt  ***
      !!                    
      !! ** Purpose :   set wetting & drying mask at tracer points 
      !!              for the current barotropic sub-step 
      !!
      !! ** Method  :   ??? 
      !!
      !! ** Action  :  ptmsk : wetting & drying t-mask
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pTmsk              ! W & D t-mask
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   phU, phV, pu, pv   ! ocean velocities and transports
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pUmsk, pVmsk       ! W & D u- and v-mask
      !
      INTEGER  ::   ji, jj   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpim1   ! not jpi-column
            IF ( phU(ji,jj) > 0._wp ) THEN   ;   pUmsk(ji,jj) = pTmsk(ji  ,jj) 
            ELSE                             ;   pUmsk(ji,jj) = pTmsk(ji+1,jj)  
            ENDIF
            phU(ji,jj) = pUmsk(ji,jj)*phU(ji,jj)
            pu (ji,jj) = pUmsk(ji,jj)*pu (ji,jj)
         END DO
      END DO
      !
      DO jj = 1, jpjm1   ! not jpj-row
         DO ji = 1, jpi
            IF ( phV(ji,jj) > 0._wp ) THEN   ;   pVmsk(ji,jj) = pTmsk(ji,jj  )
            ELSE                             ;   pVmsk(ji,jj) = pTmsk(ji,jj+1)  
            ENDIF
            phV(ji,jj) = pVmsk(ji,jj)*phV(ji,jj) 
            pv (ji,jj) = pVmsk(ji,jj)*pv (ji,jj)
         END DO
      END DO
      !
   END SUBROUTINE wad_Umsk


   SUBROUTINE wad_spg( sshn, zcpx, zcpy )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE  wad_sp  ***
      !!
      !! ** Purpose : 
      !!----------------------------------------------------------------------
      INTEGER  ::   ji ,jj               ! dummy loop indices
      LOGICAL  ::   ll_tmp1, ll_tmp2
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: sshn
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: zcpx, zcpy
      !!----------------------------------------------------------------------
      DO jj = 2, jpjm1
         DO ji = 2, jpim1 
            ll_tmp1 = MIN(  sshn(ji,jj)               ,  sshn(ji+1,jj) ) >                &
                 &      MAX( -ht_0(ji,jj)               , -ht_0(ji+1,jj) ) .AND.            &
                 &      MAX(  sshn(ji,jj) + ht_0(ji,jj) ,  sshn(ji+1,jj) + ht_0(ji+1,jj) )  &
                 &                                                         > rn_wdmin1 + rn_wdmin2
            ll_tmp2 = ( ABS( sshn(ji+1,jj)            -  sshn(ji  ,jj))  > 1.E-12 ).AND.( &
                 &      MAX(   sshn(ji,jj)              ,  sshn(ji+1,jj) ) >                &
                 &      MAX(  -ht_0(ji,jj)              , -ht_0(ji+1,jj) ) + rn_wdmin1 + rn_wdmin2 )
            IF(ll_tmp1) THEN
               zcpx(ji,jj) = 1.0_wp
            ELSEIF(ll_tmp2) THEN
               ! no worries about  sshn(ji+1,jj) -  sshn(ji  ,jj) = 0, it won't happen ! here
               zcpx(ji,jj) = ABS( (sshn(ji+1,jj) + ht_0(ji+1,jj) - sshn(ji,jj) - ht_0(ji,jj)) &
                    &    / (sshn(ji+1,jj) - sshn(ji  ,jj)) )
               zcpx(ji,jj) = max(min( zcpx(ji,jj) , 1.0_wp),0.0_wp)
            ELSE
               zcpx(ji,jj) = 0._wp
            ENDIF
            !
            ll_tmp1 = MIN(  sshn(ji,jj)               ,  sshn(ji,jj+1) ) >                &
                 &      MAX( -ht_0(ji,jj)               , -ht_0(ji,jj+1) ) .AND.            &
                 &      MAX(  sshn(ji,jj) + ht_0(ji,jj) ,  sshn(ji,jj+1) + ht_0(ji,jj+1) )  &
                 &                                                       > rn_wdmin1 + rn_wdmin2
            ll_tmp2 = ( ABS( sshn(ji,jj)              -  sshn(ji,jj+1))  > 1.E-12 ).AND.( &
                 &      MAX(   sshn(ji,jj)              ,  sshn(ji,jj+1) ) >                &
                 &      MAX(  -ht_0(ji,jj)              , -ht_0(ji,jj+1) ) + rn_wdmin1 + rn_wdmin2 )
            
            IF(ll_tmp1) THEN
               zcpy(ji,jj) = 1.0_wp
            ELSE IF(ll_tmp2) THEN
               ! no worries about  sshn(ji,jj+1) -  sshn(ji,jj  ) = 0, it won't happen ! here
               zcpy(ji,jj) = ABS( (sshn(ji,jj+1) + ht_0(ji,jj+1) - sshn(ji,jj) - ht_0(ji,jj)) &
                    &             / (sshn(ji,jj+1) - sshn(ji,jj  )) )
               zcpy(ji,jj) = MAX(  0._wp , MIN( zcpy(ji,jj) , 1.0_wp )  )
            ELSE
               zcpy(ji,jj) = 0._wp
            ENDIF
         END DO
      END DO
            
   END SUBROUTINE wad_spg
     


   SUBROUTINE dyn_drg_init( pu_RHSi, pv_RHSi, pCdU_u, pCdU_v )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_drg_init  ***
      !!                    
      !! ** Purpose : - add the baroclinic top/bottom drag contribution to 
      !!              the baroclinic part of the barotropic RHS
      !!              - compute the barotropic drag coefficients
      !!
      !! ** Method  :   computation done over the INNER domain only 
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pu_RHSi, pv_RHSi   ! baroclinic part of the barotropic RHS
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   pCdU_u , pCdU_v    ! barotropic drag coefficients
      !
      INTEGER  ::   ji, jj   ! dummy loop indices
      INTEGER  ::   ikbu, ikbv, iktu, iktv
      REAL(wp) ::   zztmp
      REAL(wp), DIMENSION(jpi,jpj) ::   zu_i, zv_i
      !!----------------------------------------------------------------------
      !
      !                    !==  Set the barotropic drag coef.  ==!
      !
      IF( ln_isfcav.OR.ln_drgice_imp ) THEN          ! top+bottom friction (ocean cavities)
         
         DO jj = 2, jpjm1
            DO ji = 2, jpim1     ! INNER domain
               pCdU_u(ji,jj) = r1_2*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) + rCdU_top(ji+1,jj)+rCdU_top(ji,jj) )
               pCdU_v(ji,jj) = r1_2*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) + rCdU_top(ji,jj+1)+rCdU_top(ji,jj) )
            END DO
         END DO
      ELSE                          ! bottom friction only
         DO jj = 2, jpjm1
            DO ji = 2, jpim1  ! INNER domain
               pCdU_u(ji,jj) = r1_2*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) )
               pCdU_v(ji,jj) = r1_2*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) )
            END DO
         END DO
      ENDIF
      !
      !                    !==  BOTTOM stress contribution from baroclinic velocities  ==!
      !
      IF( ln_bt_fw ) THEN                 ! FORWARD integration: use NOW bottom baroclinic velocities
         
         DO jj = 2, jpjm1
            DO ji = 2, jpim1  ! INNER domain
               ikbu = mbku(ji,jj)       
               ikbv = mbkv(ji,jj)    
               zu_i(ji,jj) = un(ji,jj,ikbu) - un_b(ji,jj)
               zv_i(ji,jj) = vn(ji,jj,ikbv) - vn_b(ji,jj)
            END DO
         END DO
      ELSE                                ! CENTRED integration: use BEFORE bottom baroclinic velocities
         
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! INNER domain
               ikbu = mbku(ji,jj)       
               ikbv = mbkv(ji,jj)    
               zu_i(ji,jj) = ub(ji,jj,ikbu) - ub_b(ji,jj)
               zv_i(ji,jj) = vb(ji,jj,ikbv) - vb_b(ji,jj)
            END DO
         END DO
      ENDIF
      !
      IF( ln_wd_il ) THEN      ! W/D : use the "clipped" bottom friction   !!gm   explain WHY, please !
         zztmp = -1._wp / rdtbt
         DO jj = 2, jpjm1
            DO ji = 2, jpim1    ! INNER domain
               pu_RHSi(ji,jj) = pu_RHSi(ji,jj) + zu_i(ji,jj) *  wdrampu(ji,jj) * MAX(                                 & 
                    &                              r1_hu_n(ji,jj) * r1_2*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) ) , zztmp  )
               pv_RHSi(ji,jj) = pv_RHSi(ji,jj) + zv_i(ji,jj) *  wdrampv(ji,jj) * MAX(                                 & 
                    &                              r1_hv_n(ji,jj) * r1_2*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) ) , zztmp  )
            END DO
         END DO
      ELSE                    ! use "unclipped" drag (even if explicit friction is used in 3D calculation)
         
         DO jj = 2, jpjm1
            DO ji = 2, jpim1    ! INNER domain
               pu_RHSi(ji,jj) = pu_RHSi(ji,jj) + r1_hu_n(ji,jj) * r1_2*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) ) * zu_i(ji,jj)
               pv_RHSi(ji,jj) = pv_RHSi(ji,jj) + r1_hv_n(ji,jj) * r1_2*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) ) * zv_i(ji,jj)
            END DO
         END DO
      END IF
      !
      !                    !==  TOP stress contribution from baroclinic velocities  ==!   (no W/D case)
      !
      IF( ln_isfcav.OR.ln_drgice_imp ) THEN
         !
         IF( ln_bt_fw ) THEN                ! FORWARD integration: use NOW top baroclinic velocity
            
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! INNER domain
                  iktu = miku(ji,jj)
                  iktv = mikv(ji,jj)
                  zu_i(ji,jj) = un(ji,jj,iktu) - un_b(ji,jj)
                  zv_i(ji,jj) = vn(ji,jj,iktv) - vn_b(ji,jj)
               END DO
            END DO
         ELSE                                ! CENTRED integration: use BEFORE top baroclinic velocity
            
            DO jj = 2, jpjm1
               DO ji = 2, jpim1      ! INNER domain
                  iktu = miku(ji,jj)
                  iktv = mikv(ji,jj)
                  zu_i(ji,jj) = ub(ji,jj,iktu) - ub_b(ji,jj)
                  zv_i(ji,jj) = vb(ji,jj,iktv) - vb_b(ji,jj)
               END DO
            END DO
         ENDIF
         !
         !                    ! use "unclipped" top drag (even if explicit friction is used in 3D calculation)
         
         DO jj = 2, jpjm1
            DO ji = 2, jpim1    ! INNER domain
               pu_RHSi(ji,jj) = pu_RHSi(ji,jj) + r1_hu_n(ji,jj) * r1_2*( rCdU_top(ji+1,jj)+rCdU_top(ji,jj) ) * zu_i(ji,jj)
               pv_RHSi(ji,jj) = pv_RHSi(ji,jj) + r1_hv_n(ji,jj) * r1_2*( rCdU_top(ji,jj+1)+rCdU_top(ji,jj) ) * zv_i(ji,jj)
            END DO
         END DO
         !
      ENDIF
      !
   END SUBROUTINE dyn_drg_init

   SUBROUTINE ts_bck_interp( jn, ll_init,       &   ! <== in
      &                      za0, za1, za2, za3 )   ! ==> out
      !!----------------------------------------------------------------------
      INTEGER ,INTENT(in   ) ::   jn                   ! index of sub time step
      LOGICAL ,INTENT(in   ) ::   ll_init              !
      REAL(wp),INTENT(  out) ::   za0, za1, za2, za3   ! Half-step back interpolation coefficient
      !
      REAL(wp) ::   zepsilon, zgamma                   !   -      -
      !!----------------------------------------------------------------------
      !                             ! set Half-step back interpolation coefficient
      IF    ( jn==1 .AND. ll_init ) THEN   !* Forward-backward
         za0 = 1._wp                        
         za1 = 0._wp                           
         za2 = 0._wp
         za3 = 0._wp
      ELSEIF( jn==2 .AND. ll_init ) THEN   !* AB2-AM3 Coefficients; bet=0 ; gam=-1/6 ; eps=1/12
         za0 = 1.0833333333333_wp                 ! za0 = 1-gam-eps
         za1 =-0.1666666666666_wp                 ! za1 = gam
         za2 = 0.0833333333333_wp                 ! za2 = eps
         za3 = 0._wp              
      ELSE                                 !* AB3-AM4 Coefficients; bet=0.281105 ; eps=0.013 ; gam=0.0880 
         IF( rn_bt_alpha == 0._wp ) THEN      ! Time diffusion  
            za0 = 0.614_wp                        ! za0 = 1/2 +   gam + 2*eps
            za1 = 0.285_wp                        ! za1 = 1/2 - 2*gam - 3*eps
            za2 = 0.088_wp                        ! za2 = gam
            za3 = 0.013_wp                        ! za3 = eps
         ELSE                                 ! no time diffusion
            zepsilon = 0.00976186_wp - 0.13451357_wp * rn_bt_alpha
            zgamma   = 0.08344500_wp - 0.51358400_wp * rn_bt_alpha
            za0 = 0.5_wp + zgamma + 2._wp * rn_bt_alpha + 2._wp * zepsilon
            za1 = 1._wp - za0 - zgamma - zepsilon
            za2 = zgamma
            za3 = zepsilon
         ENDIF 
      ENDIF
   END SUBROUTINE ts_bck_interp


   !!======================================================================
END MODULE dynspg_ts
