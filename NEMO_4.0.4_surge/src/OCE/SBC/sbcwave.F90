MODULE sbcwave
   !!======================================================================
   !!                       ***  MODULE  sbcwave  ***
   !! Wave module 
   !!======================================================================
   !! History :  3.3  !  2011-09  (M. Adani)  Original code: Drag Coefficient 
   !!         :  3.4  !  2012-10  (M. Adani)  Stokes Drift 
   !!            3.6  !  2014-09  (E. Clementi,P. Oddo) New Stokes Drift Computation
   !!             -   !  2016-12  (G. Madec, E. Clementi) update Stoke drift computation
   !!                                                    + add sbc_wave_ini routine
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_stokes    : calculate 3D Stokes-drift velocities
   !!   sbc_wave      : wave data from wave model in netcdf files 
   !!   sbc_wave_init : initialisation fo surface waves 
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants 
   USE oce            ! ocean variables
   USE sbc_oce	       ! Surface boundary condition: ocean fields
   USE zdf_oce,  ONLY : ln_zdfswm
   USE bdy_oce        ! open boundary condition variables
   USE domvvl         ! domain: variable volume layers
   !
   USE iom            ! I/O manager library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE fldread	       ! read input fields

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_stokes      ! routine called in sbccpl
   PUBLIC   sbc_wstress     ! routine called in sbcmod 
   PUBLIC   sbc_wave        ! routine called in sbcmod
   PUBLIC   sbc_wave_init   ! routine called in sbcmod
   
   ! Variables checking if the wave parameters are coupled (if not, they are read from file)
   LOGICAL, PUBLIC ::   cpl_hsig   = .FALSE.
   LOGICAL, PUBLIC ::   cpl_phioc  = .FALSE.
   LOGICAL, PUBLIC ::   cpl_sdrftx = .FALSE.
   LOGICAL, PUBLIC ::   cpl_sdrfty = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wper   = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wfreq  = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wnum   = .FALSE.
   LOGICAL, PUBLIC ::   cpl_tauwoc = .FALSE.
   LOGICAL, PUBLIC ::   cpl_tauw   = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wdrag  = .FALSE.

   INTEGER ::   jpfld    ! number of files to read for stokes drift
   INTEGER ::   jp_usd   ! index of stokes drift  (i-component) (m/s)    at T-point
   INTEGER ::   jp_vsd   ! index of stokes drift  (j-component) (m/s)    at T-point
   INTEGER ::   jp_hsw   ! index of significant wave hight      (m)      at T-point
   INTEGER ::   jp_wmp   ! index of mean wave period            (s)      at T-point
   INTEGER ::   jp_wfr   ! index of wave peak frequency         (1/s)    at T-point

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_cd      ! structure of input fields (file informations, fields read) Drag Coefficient
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sd      ! structure of input fields (file informations, fields read) Stokes Drift
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_wn      ! structure of input fields (file informations, fields read) wave number for Qiao
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_tauwoc  ! structure of input fields (file informations, fields read) normalized wave stress into the ocean
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_tauw    ! structure of input fields (file informations, fields read) ocean stress components from wave model

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   cdn_wave            !:
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   hsw, wmp, wnum      !: 
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   wfreq               !: 
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tauoc_wave          !:  
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tauw_x, tauw_y      !:  
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tsd2d               !: 
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   div_sd              !: barotropic stokes drift divergence
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   ut0sd, vt0sd        !: surface Stokes drift velocities at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   usd  , vsd  , wsd   !: Stokes drift velocities at u-, v- & w-points, resp.

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_stokes( )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_stokes  ***
      !!
      !! ** Purpose :   compute the 3d Stokes Drift according to Breivik et al.,
      !!                2014 (DOI: 10.1175/JPO-D-14-0020.1)
      !!
      !! ** Method  : - Calculate Stokes transport speed 
      !!              - Calculate horizontal divergence 
      !!              - Integrate the horizontal divergenze from the bottom 
      !! ** action  
      !!---------------------------------------------------------------------
      INTEGER  ::   jj, ji, jk   ! dummy loop argument
      INTEGER  ::   ik           ! local integer 
      REAL(wp) ::  ztransp, zfac, zsp0
      REAL(wp) ::  zdepth, zsqrt_depth,  zexp_depth, z_two_thirds, zsqrtpi !sqrt of pi
      REAL(wp) ::  zbot_u, zbot_v, zkb_u, zkb_v, zke3_u, zke3_v, zda_u, zda_v
      REAL(wp) ::  zstokes_psi_u_bot, zstokes_psi_v_bot
      REAL(wp) ::  zdep_u, zdep_v, zkh_u, zkh_v
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   zk_t, zk_u, zk_v, zu0_sd, zv0_sd     ! 2D workspace
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   zstokes_psi_u_top, zstokes_psi_v_top ! 2D workspace
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ze3divh                              ! 3D workspace
      !!---------------------------------------------------------------------
      !
      ALLOCATE( ze3divh(jpi,jpj,jpk) )
      ALLOCATE( zk_t(jpi,jpj), zk_u(jpi,jpj), zk_v(jpi,jpj), zu0_sd(jpi,jpj), zv0_sd(jpi,jpj) )
      !
      ! select parameterization for the calculation of vertical Stokes drift
      ! exp. wave number at t-point
      IF( ll_st_bv_li ) THEN   ! (Eq. (19) in Breivik et al. (2014) )
         zfac = 2.0_wp * rpi / 16.0_wp
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Stokes drift velocity estimated from Hs and Tmean
               ztransp = zfac * hsw(ji,jj)*hsw(ji,jj) / MAX( wmp(ji,jj), 0.0000001_wp )
               ! Stokes surface speed
               tsd2d(ji,jj) = SQRT( ut0sd(ji,jj)*ut0sd(ji,jj) + vt0sd(ji,jj)*vt0sd(ji,jj))
               ! Wavenumber scale
               zk_t(ji,jj) = ABS( tsd2d(ji,jj) ) / MAX( ABS( 5.97_wp*ztransp ), 0.0000001_wp )
            END DO
         END DO
         DO jj = 1, jpjm1              ! exp. wave number & Stokes drift velocity at u- & v-points
            DO ji = 1, jpim1
               zk_u(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji+1,jj) )
               zk_v(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji,jj+1) )
               !
               zu0_sd(ji,jj) = 0.5_wp * ( ut0sd(ji,jj) + ut0sd(ji+1,jj) )
               zv0_sd(ji,jj) = 0.5_wp * ( vt0sd(ji,jj) + vt0sd(ji,jj+1) )
            END DO
         END DO
      ELSE IF( ll_st_peakfr ) THEN    ! peak wave number calculated from the peak frequency received by the wave model
         DO jj = 1, jpj
            DO ji = 1, jpi
               zk_t(ji,jj) = ( 2.0_wp * rpi * wfreq(ji,jj) ) * ( 2.0_wp * rpi * wfreq(ji,jj) ) / grav
            END DO
         END DO
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               zk_u(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji+1,jj) )
               zk_v(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji,jj+1) )
               !
               zu0_sd(ji,jj) = 0.5_wp * ( ut0sd(ji,jj) + ut0sd(ji+1,jj) )
               zv0_sd(ji,jj) = 0.5_wp * ( vt0sd(ji,jj) + vt0sd(ji,jj+1) )
            END DO
         END DO
      ENDIF
      !
      !                       !==  horizontal Stokes Drift 3D velocity  ==!
      IF( ll_st_bv2014 ) THEN
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zdep_u = 0.5_wp * ( gdept_n(ji,jj,jk) + gdept_n(ji+1,jj,jk) )
                  zdep_v = 0.5_wp * ( gdept_n(ji,jj,jk) + gdept_n(ji,jj+1,jk) )
                  !                          
                  zkh_u = zk_u(ji,jj) * zdep_u     ! k * depth
                  zkh_v = zk_v(ji,jj) * zdep_v
                  !                                ! Depth attenuation
                  zda_u = EXP( -2.0_wp*zkh_u ) / ( 1.0_wp + 8.0_wp*zkh_u )
                  zda_v = EXP( -2.0_wp*zkh_v ) / ( 1.0_wp + 8.0_wp*zkh_v )
                  !
                  usd(ji,jj,jk) = zda_u * zu0_sd(ji,jj) * umask(ji,jj,jk)
                  vsd(ji,jj,jk) = zda_v * zv0_sd(ji,jj) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ELSE IF( ll_st_li2017 .OR. ll_st_peakfr ) THEN
         ALLOCATE( zstokes_psi_u_top(jpi,jpj), zstokes_psi_v_top(jpi,jpj) )
         DO jj = 1, jpjm1              ! exp. wave number & Stokes drift velocity at u- & v-points
            DO ji = 1, jpim1
               zstokes_psi_u_top(ji,jj) = 0._wp
               zstokes_psi_v_top(ji,jj) = 0._wp
            END DO
         END DO
         zsqrtpi = SQRT(rpi)
         z_two_thirds = 2.0_wp / 3.0_wp
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zbot_u = ( gdepw_n(ji,jj,jk+1) + gdepw_n(ji+1,jj,jk+1) )  ! 2 * bottom depth
                  zbot_v = ( gdepw_n(ji,jj,jk+1) + gdepw_n(ji,jj+1,jk+1) )  ! 2 * bottom depth
                  zkb_u  = zk_u(ji,jj) * zbot_u                             ! 2 * k * bottom depth
                  zkb_v  = zk_v(ji,jj) * zbot_v                             ! 2 * k * bottom depth
                  !
                  zke3_u = MAX(1.e-8_wp, 2.0_wp * zk_u(ji,jj) * e3u_n(ji,jj,jk))     ! 2k * thickness
                  zke3_v = MAX(1.e-8_wp, 2.0_wp * zk_v(ji,jj) * e3v_n(ji,jj,jk))     ! 2k * thickness

                  ! Depth attenuation .... do u component first..
                  zdepth      = zkb_u
                  zsqrt_depth = SQRT(zdepth)
                  zexp_depth  = EXP(-zdepth)
                  zstokes_psi_u_bot = 1.0_wp - zexp_depth  &
                       &              - z_two_thirds * ( zsqrtpi*zsqrt_depth*zdepth*ERFC(zsqrt_depth) &
                       &              + 1.0_wp - (1.0_wp + zdepth)*zexp_depth )
                  zda_u                    = ( zstokes_psi_u_bot - zstokes_psi_u_top(ji,jj) ) / zke3_u
                  zstokes_psi_u_top(ji,jj) =   zstokes_psi_u_bot

                  !         ... and then v component
                  zdepth      =zkb_v
                  zsqrt_depth = SQRT(zdepth)
                  zexp_depth  = EXP(-zdepth)
                  zstokes_psi_v_bot = 1.0_wp - zexp_depth  &
                       &              - z_two_thirds * ( zsqrtpi*zsqrt_depth*zdepth*ERFC(zsqrt_depth) &
                       &              + 1.0_wp - (1.0_wp + zdepth)*zexp_depth )
                  zda_v                    = ( zstokes_psi_v_bot - zstokes_psi_v_top(ji,jj) ) / zke3_v
                  zstokes_psi_v_top(ji,jj) =   zstokes_psi_v_bot
                  !
                  usd(ji,jj,jk) = zda_u * zu0_sd(ji,jj) * umask(ji,jj,jk)
                  vsd(ji,jj,jk) = zda_v * zv0_sd(ji,jj) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
         DEALLOCATE( zstokes_psi_u_top, zstokes_psi_v_top )
      ENDIF

      CALL lbc_lnk_multi( 'sbcwave', usd, 'U', -1., vsd, 'V', -1. )

      !
      !                       !==  vertical Stokes Drift 3D velocity  ==!
      !
      DO jk = 1, jpkm1               ! Horizontal e3*divergence
         DO jj = 2, jpj
            DO ji = fs_2, jpi
               ze3divh(ji,jj,jk) = (  e2u(ji  ,jj) * e3u_n(ji  ,jj,jk) * usd(ji  ,jj,jk)    &
                  &                 - e2u(ji-1,jj) * e3u_n(ji-1,jj,jk) * usd(ji-1,jj,jk)    &
                  &                 + e1v(ji,jj  ) * e3v_n(ji,jj  ,jk) * vsd(ji,jj  ,jk)    &
                  &                 - e1v(ji,jj-1) * e3v_n(ji,jj-1,jk) * vsd(ji,jj-1,jk)  ) * r1_e1e2t(ji,jj)
            END DO
         END DO
      END DO
      !
#if defined key_agrif
      IF( .NOT. Agrif_Root() ) THEN
         IF( l_Westedge )    ze3divh( 2:nbghostcells+1,:        ,:) = 0._wp      ! west
         IF( l_Eastedge )    ze3divh( nlci-nbghostcells:nlci-1,:,:) = 0._wp      ! east
         IF( l_Southedge )   ze3divh( :,2:nbghostcells+1        ,:) = 0._wp      ! south
         IF( l_Northedge )   ze3divh( :,nlcj-nbghostcells:nlcj-1,:) = 0._wp      ! north
      ENDIF
#endif
      !
      CALL lbc_lnk( 'sbcwave', ze3divh, 'T', 1. )
      !
      IF( ln_linssh ) THEN   ;   ik = 1   ! none zero velocity through the sea surface
      ELSE                   ;   ik = 2   ! w=0 at the surface (set one for all in sbc_wave_init)
      ENDIF
      DO jk = jpkm1, ik, -1          ! integrate from the bottom the hor. divergence (NB: at k=jpk w is always zero)
         wsd(:,:,jk) = wsd(:,:,jk+1) - ze3divh(:,:,jk)
      END DO
      !
      IF( ln_bdy ) THEN
         DO jk = 1, jpkm1
            wsd(:,:,jk) = wsd(:,:,jk) * bdytmask(:,:)
         END DO
      ENDIF
      !                       !==  Horizontal divergence of barotropic Stokes transport  ==!
      div_sd(:,:) = 0._wp
      DO jk = 1, jpkm1                                 ! 
        div_sd(:,:) = div_sd(:,:) + ze3divh(:,:,jk)
      END DO
      !
      CALL iom_put( "ustokes",  usd  )
      CALL iom_put( "vstokes",  vsd  )
      CALL iom_put( "wstokes",  wsd  )
      !
      DEALLOCATE( ze3divh )
      DEALLOCATE( zk_t, zk_u, zk_v, zu0_sd, zv0_sd )
      !
   END SUBROUTINE sbc_stokes


   SUBROUTINE sbc_wstress( )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_wstress  ***
      !!
      !! ** Purpose :   Updates the ocean momentum modified by waves
      !!
      !! ** Method  : - Calculate u,v components of stress depending on stress
      !!                model 
      !!              - Calculate the stress module
      !!              - The wind module is not modified by waves 
      !! ** action  
      !!---------------------------------------------------------------------
      INTEGER  ::   jj, ji   ! dummy loop argument
      !
      IF( ln_tauwoc ) THEN
         utau(:,:) = utau(:,:)*tauoc_wave(:,:)
         vtau(:,:) = vtau(:,:)*tauoc_wave(:,:)
         taum(:,:) = taum(:,:)*tauoc_wave(:,:)
      ENDIF
      !
      IF( ln_tauw ) THEN
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               ! Stress components at u- & v-points
               utau(ji,jj) = 0.5_wp * ( tauw_x(ji,jj) + tauw_x(ji+1,jj) )
               vtau(ji,jj) = 0.5_wp * ( tauw_y(ji,jj) + tauw_y(ji,jj+1) )
               !
               ! Stress module at t points
               taum(ji,jj) = SQRT( tauw_x(ji,jj)*tauw_x(ji,jj) + tauw_y(ji,jj)*tauw_y(ji,jj) )
            END DO
         END DO
         CALL lbc_lnk_multi( 'sbcwave', utau(:,:), 'U', -1. , vtau(:,:), 'V', -1. , taum(:,:) , 'T', -1. )
      ENDIF
      !
   END SUBROUTINE sbc_wstress


   SUBROUTINE sbc_wave( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_wave  ***
      !!
      !! ** Purpose :   read wave parameters from wave model  in netcdf files.
      !!
      !! ** Method  : - Read namelist namsbc_wave
      !!              - Read Cd_n10 fields in netcdf files 
      !!              - Read stokes drift 2d in netcdf files 
      !!              - Read wave number in netcdf files 
      !!              - Compute 3d stokes drift using Breivik et al.,2014
      !!                formulation
      !! ** action  
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt   ! ocean time step
      !!---------------------------------------------------------------------
      !
      IF( ln_cdgw .AND. .NOT. cpl_wdrag ) THEN     !==  Neutral drag coefficient  ==!
         CALL fld_read( kt, nn_fsbc, sf_cd )             ! read from external forcing
         cdn_wave(:,:) = sf_cd(1)%fnow(:,:,1) * tmask(:,:,1)
      ENDIF

      IF( ln_tauwoc .AND. .NOT. cpl_tauwoc ) THEN  !==  Wave induced stress  ==!
         CALL fld_read( kt, nn_fsbc, sf_tauwoc )         ! read wave norm stress from external forcing
         tauoc_wave(:,:) = sf_tauwoc(1)%fnow(:,:,1) * tmask(:,:,1)
      ENDIF

      IF( ln_tauw .AND. .NOT. cpl_tauw ) THEN      !==  Wave induced stress  ==!
         CALL fld_read( kt, nn_fsbc, sf_tauw )           ! read ocean stress components from external forcing (T grid)
         tauw_x(:,:) = sf_tauw(1)%fnow(:,:,1) * tmask(:,:,1)
         tauw_y(:,:) = sf_tauw(2)%fnow(:,:,1) * tmask(:,:,1)
      ENDIF

      IF( ln_sdw )  THEN                           !==  Computation of the 3d Stokes Drift  ==! 
         !
         IF( jpfld > 0 ) THEN                            ! Read from file only if the field is not coupled
            CALL fld_read( kt, nn_fsbc, sf_sd )          ! read wave parameters from external forcing
            IF( jp_hsw > 0 )   hsw  (:,:) = sf_sd(jp_hsw)%fnow(:,:,1) * tmask(:,:,1)  ! significant wave height
            IF( jp_wmp > 0 )   wmp  (:,:) = sf_sd(jp_wmp)%fnow(:,:,1) * tmask(:,:,1)  ! wave mean period
            IF( jp_wfr > 0 )   wfreq(:,:) = sf_sd(jp_wfr)%fnow(:,:,1) * tmask(:,:,1)  ! Peak wave frequency
            IF( jp_usd > 0 )   ut0sd(:,:) = sf_sd(jp_usd)%fnow(:,:,1) * tmask(:,:,1)  ! 2D zonal Stokes Drift at T point
            IF( jp_vsd > 0 )   vt0sd(:,:) = sf_sd(jp_vsd)%fnow(:,:,1) * tmask(:,:,1)  ! 2D meridional Stokes Drift at T point
         ENDIF
         !
         ! Read also wave number if needed, so that it is available in coupling routines
         IF( ln_zdfswm .AND. .NOT.cpl_wnum ) THEN
            CALL fld_read( kt, nn_fsbc, sf_wn )          ! read wave parameters from external forcing
            wnum(:,:) = sf_wn(1)%fnow(:,:,1) * tmask(:,:,1)
         ENDIF
           
         ! Calculate only if required fields have been read
         ! In coupled wave model-NEMO case the call is done after coupling
         !
         IF( ( ll_st_bv_li   .AND. jp_hsw>0 .AND. jp_wmp>0 .AND. jp_usd>0 .AND. jp_vsd>0 ) .OR. &
           & ( ll_st_peakfr  .AND. jp_wfr>0 .AND. jp_usd>0 .AND. jp_vsd>0                ) ) CALL sbc_stokes()
         !
      ENDIF
      !
   END SUBROUTINE sbc_wave


   SUBROUTINE sbc_wave_init
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_wave_init  ***
      !!
      !! ** Purpose :   read wave parameters from wave model  in netcdf files.
      !!
      !! ** Method  : - Read namelist namsbc_wave
      !!              - Read Cd_n10 fields in netcdf files 
      !!              - Read stokes drift 2d in netcdf files 
      !!              - Read wave number in netcdf files 
      !!              - Compute 3d stokes drift using Breivik et al.,2014
      !!                formulation
      !! ** action  
      !!---------------------------------------------------------------------
      INTEGER ::   ierror, ios   ! local integer
      INTEGER ::   ifpr
      !!
      CHARACTER(len=100)     ::  cn_dir                            ! Root directory for location of drag coefficient files
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::   slf_i, slf_j     ! array of namelist informations on the fields to read
      TYPE(FLD_N)            ::  sn_cdg, sn_usd, sn_vsd,  &
                             &   sn_hsw, sn_wmp, sn_wfr, sn_wnum, &
                             &   sn_tauwoc, sn_tauwx, sn_tauwy     ! informations about the fields to be read
      !
      NAMELIST/namsbc_wave/  sn_cdg, cn_dir, sn_usd, sn_vsd, sn_hsw, sn_wmp, sn_wfr, &
                          &  sn_wnum, sn_tauwoc, sn_tauwx, sn_tauwy
      !!---------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namsbc_wave in reference namelist : File for drag coeff. from wave model
      READ  ( numnam_ref, namsbc_wave, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_wave in reference namelist' )
         
      REWIND( numnam_cfg )              ! Namelist namsbc_wave in configuration namelist : File for drag coeff. from wave model
      READ  ( numnam_cfg, namsbc_wave, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_wave in configuration namelist' )
      IF(lwm) WRITE ( numond, namsbc_wave )
      !
      IF( ln_cdgw ) THEN
         IF( .NOT. cpl_wdrag ) THEN
            ALLOCATE( sf_cd(1), STAT=ierror )               !* allocate and fill sf_wave with sn_cdg
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wave structure' )
            !
                                   ALLOCATE( sf_cd(1)%fnow(jpi,jpj,1)   )
            IF( sn_cdg%ln_tint )   ALLOCATE( sf_cd(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_cd, (/ sn_cdg /), cn_dir, 'sbc_wave_init', 'Wave module ', 'namsbc_wave' )
         ENDIF
         ALLOCATE( cdn_wave(jpi,jpj) )
      ENDIF

      IF( ln_tauwoc ) THEN
         IF( .NOT. cpl_tauwoc ) THEN
            ALLOCATE( sf_tauwoc(1), STAT=ierror )           !* allocate and fill sf_wave with sn_tauwoc
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wave structure' )
            !
                                     ALLOCATE( sf_tauwoc(1)%fnow(jpi,jpj,1)   )
            IF( sn_tauwoc%ln_tint )  ALLOCATE( sf_tauwoc(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_tauwoc, (/ sn_tauwoc /), cn_dir, 'sbc_wave_init', 'Wave module', 'namsbc_wave' )
         ENDIF
         ALLOCATE( tauoc_wave(jpi,jpj) )
      ENDIF

      IF( ln_tauw ) THEN
         IF( .NOT. cpl_tauw ) THEN
            ALLOCATE( sf_tauw(2), STAT=ierror )           !* allocate and fill sf_wave with sn_tauwx/y
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_tauw structure' )
            !
            ALLOCATE( slf_j(2) )
            slf_j(1) = sn_tauwx
            slf_j(2) = sn_tauwy
                                    ALLOCATE( sf_tauw(1)%fnow(jpi,jpj,1)   )
                                    ALLOCATE( sf_tauw(2)%fnow(jpi,jpj,1)   )
            IF( slf_j(1)%ln_tint )  ALLOCATE( sf_tauw(1)%fdta(jpi,jpj,1,2) )
            IF( slf_j(2)%ln_tint )  ALLOCATE( sf_tauw(2)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_tauw, (/ slf_j /), cn_dir, 'sbc_wave_init', 'read wave input', 'namsbc_wave' )
         ENDIF
         ALLOCATE( tauw_x(jpi,jpj) )
         ALLOCATE( tauw_y(jpi,jpj) )
      ENDIF

      IF( ln_sdw ) THEN   ! Find out how many fields have to be read from file if not coupled
         jpfld=0
         jp_usd=0   ;   jp_vsd=0   ;   jp_hsw=0   ;   jp_wmp=0   ;   jp_wfr=0
         IF( .NOT. cpl_sdrftx ) THEN
            jpfld  = jpfld + 1
            jp_usd = jpfld
         ENDIF
         IF( .NOT. cpl_sdrfty ) THEN
            jpfld  = jpfld + 1
            jp_vsd = jpfld
         ENDIF
         IF( .NOT. cpl_hsig  .AND. ll_st_bv_li  ) THEN
            jpfld  = jpfld + 1
            jp_hsw = jpfld
         ENDIF
         IF( .NOT. cpl_wper  .AND. ll_st_bv_li  ) THEN
            jpfld  = jpfld + 1
            jp_wmp = jpfld
         ENDIF
         IF( .NOT. cpl_wfreq .AND. ll_st_peakfr ) THEN
            jpfld  = jpfld + 1
            jp_wfr = jpfld
         ENDIF

         ! Read from file only the non-coupled fields 
         IF( jpfld > 0 ) THEN
            ALLOCATE( slf_i(jpfld) )
            IF( jp_usd > 0 )   slf_i(jp_usd) = sn_usd
            IF( jp_vsd > 0 )   slf_i(jp_vsd) = sn_vsd
            IF( jp_hsw > 0 )   slf_i(jp_hsw) = sn_hsw
            IF( jp_wmp > 0 )   slf_i(jp_wmp) = sn_wmp
            IF( jp_wfr > 0 )   slf_i(jp_wfr) = sn_wfr

            ALLOCATE( sf_sd(jpfld), STAT=ierror )   !* allocate and fill sf_sd with stokes drift
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wave structure' )
            !
            DO ifpr= 1, jpfld
               ALLOCATE( sf_sd(ifpr)%fnow(jpi,jpj,1) )
               IF( slf_i(ifpr)%ln_tint )   ALLOCATE( sf_sd(ifpr)%fdta(jpi,jpj,1,2) )
            END DO
            !
            CALL fld_fill( sf_sd, slf_i, cn_dir, 'sbc_wave_init', 'Wave module ', 'namsbc_wave' )
         ENDIF
         ALLOCATE( usd  (jpi,jpj,jpk), vsd  (jpi,jpj,jpk), wsd(jpi,jpj,jpk) )
         ALLOCATE( hsw  (jpi,jpj)    , wmp  (jpi,jpj)     )
         ALLOCATE( wfreq(jpi,jpj) )
         ALLOCATE( ut0sd(jpi,jpj)    , vt0sd(jpi,jpj)     )
         ALLOCATE( div_sd(jpi,jpj) )
         ALLOCATE( tsd2d (jpi,jpj) )

         ut0sd(:,:) = 0._wp
         vt0sd(:,:) = 0._wp
         hsw(:,:) = 0._wp
         wmp(:,:) = 0._wp

         usd(:,:,:) = 0._wp
         vsd(:,:,:) = 0._wp
         wsd(:,:,:) = 0._wp
         ! Wave number needed only if ln_zdfswm=T
         IF( .NOT. cpl_wnum ) THEN
            ALLOCATE( sf_wn(1), STAT=ierror )           !* allocate and fill sf_wave with sn_wnum
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable toallocate sf_wave structure' )
                                   ALLOCATE( sf_wn(1)%fnow(jpi,jpj,1)   )
            IF( sn_wnum%ln_tint )  ALLOCATE( sf_wn(1)%fdta(jpi,jpj,1,2) )
            CALL fld_fill( sf_wn, (/ sn_wnum /), cn_dir, 'sbc_wave', 'Wave module', 'namsbc_wave' )
         ENDIF
         ALLOCATE( wnum(jpi,jpj) )
      ENDIF
      !
   END SUBROUTINE sbc_wave_init

   !!======================================================================
END MODULE sbcwave
