MODULE sbcisf
   !!======================================================================
   !!                       ***  MODULE  sbcisf  ***
   !! Surface module :  update surface ocean boundary condition under ice
   !!                   shelf
   !!======================================================================
   !! History :  3.2  !  2011-02  (C.Harris  ) Original code isf cav
   !!            X.X  !  2006-02  (C. Wang   ) Original code bg03
   !!            3.4  !  2013-03  (P. Mathiot) Merging + parametrization
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_isf       : update sbc under ice shelf
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE eosbn2         ! equation of state
   USE sbc_oce        ! surface boundary condition: ocean fields
   USE zdfdrg         ! vertical physics: top/bottom drag coef.
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE fldread        ! read input field at current time step
   USE lbclnk         !
   USE lib_fortran    ! glob_sum

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_isf, sbc_isf_init, sbc_isf_div, sbc_isf_alloc  ! routine called in sbcmod and divhor

   ! public in order to be able to output then 

   REAL(wp), PUBLIC ::   rn_hisf_tbl                 !: thickness of top boundary layer [m]
   INTEGER , PUBLIC ::   nn_isf                      !: flag to choose between explicit/param/specified  
   INTEGER , PUBLIC ::   nn_isfblk                   !: flag to choose the bulk formulation to compute the ice shelf melting
   INTEGER , PUBLIC ::   nn_gammablk                 !: flag to choose how the exchange coefficient is computed
   REAL(wp), PUBLIC ::   rn_gammat0                  !: temperature exchange coeficient []
   REAL(wp), PUBLIC ::   rn_gammas0                  !: salinity    exchange coeficient []

   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   misfkt   , misfkb        !: Level of ice shelf base
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   rzisf_tbl                !: depth of calving front (shallowest point) nn_isf ==2/3
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   rhisf_tbl, rhisf_tbl_0   !: thickness of tbl  [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   r1_hisf_tbl              !: 1/thickness of tbl
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ralpha                   !: proportion of bottom cell influenced by tbl 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   risfLeff                 !: effective length (Leff) BG03 nn_isf==2
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ttbl, stbl, utbl, vtbl   !: top boundary layer variable at T point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   qisf                     !: net heat flux from ice shelf      [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   risf_tsc_b, risf_tsc     !: before and now T & S isf contents [K.m/s & PSU.m/s]  

   LOGICAL, PUBLIC ::   l_isfcpl = .false.       !: isf recieved from oasis

   REAL(wp), PUBLIC, SAVE ::   rcpisf   = 2000.0_wp     !: specific heat of ice shelf             [J/kg/K]
   REAL(wp), PUBLIC, SAVE ::   rkappa   = 1.54e-6_wp    !: heat diffusivity through the ice-shelf [m2/s]
   REAL(wp), PUBLIC, SAVE ::   rhoisf   = 920.0_wp      !: volumic mass of ice shelf              [kg/m3]
   REAL(wp), PUBLIC, SAVE ::   tsurf    = -20.0_wp      !: air temperature on top of ice shelf    [C]
   REAL(wp), PUBLIC, SAVE ::   rLfusisf = 0.334e6_wp    !: latent heat of fusion of ice shelf     [J/kg]

!: Variable used in fldread to read the forcing file (nn_isf == 4 .OR. nn_isf == 3)
   CHARACTER(len=100), PUBLIC           :: cn_dirisf  = './' !: Root directory for location of ssr files
   TYPE(FLD_N)       , PUBLIC           :: sn_fwfisf         !: information about the isf melting file to be read
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_fwfisf
   TYPE(FLD_N)       , PUBLIC           :: sn_rnfisf         !: information about the isf melting param.   file to be read
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_rnfisf           
   TYPE(FLD_N)       , PUBLIC           :: sn_depmax_isf     !: information about the grounding line depth file to be read
   TYPE(FLD_N)       , PUBLIC           :: sn_depmin_isf     !: information about the calving   line depth file to be read
   TYPE(FLD_N)       , PUBLIC           :: sn_Leff_isf       !: information about the effective length     file to be read
   
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
 
  SUBROUTINE sbc_isf( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_isf  ***
      !!
      !! ** Purpose : Compute Salt and Heat fluxes related to ice_shelf 
      !!              melting and freezing 
      !!
      !! ** Method  :  4 parameterizations are available according to nn_isf 
      !!               nn_isf = 1 : Realistic ice_shelf formulation
      !!                        2 : Beckmann & Goose parameterization
      !!                        3 : Specified runoff in deptht (Mathiot & al. )
      !!                        4 : specified fwf and heat flux forcing beneath the ice shelf
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      INTEGER ::   ji, jj, jk   ! loop index
      INTEGER ::   ikt, ikb     ! local integers
      REAL(wp), DIMENSION(jpi,jpj) ::   zt_frz, zdep   ! freezing temperature (zt_frz) at depth (zdep) 
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   zqhcisf2d
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   zfwfisf3d, zqhcisf3d, zqlatisf3d
      !!---------------------------------------------------------------------
      !
      IF( MOD( kt-1, nn_fsbc) == 0 ) THEN    ! compute salt and heat flux
         !
         SELECT CASE ( nn_isf )
         CASE ( 1 )    ! realistic ice shelf formulation
            ! compute T/S/U/V for the top boundary layer
            CALL sbc_isf_tbl(tsn(:,:,:,jp_tem),ttbl(:,:),'T')
            CALL sbc_isf_tbl(tsn(:,:,:,jp_sal),stbl(:,:),'T')
            CALL sbc_isf_tbl(un(:,:,:)        ,utbl(:,:),'U')
            CALL sbc_isf_tbl(vn(:,:,:)        ,vtbl(:,:),'V')
            ! iom print
            CALL iom_put('ttbl',ttbl(:,:))
            CALL iom_put('stbl',stbl(:,:))
            CALL iom_put('utbl',utbl(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:))
            CALL iom_put('vtbl',vtbl(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:))
            ! compute fwf and heat flux
            ! compute fwf and heat flux
            IF( .NOT.l_isfcpl ) THEN    ;   CALL sbc_isf_cav (kt)
            ELSE                        ;   qisf(:,:)  = fwfisf(:,:) * rLfusisf  ! heat        flux
            ENDIF
            !
         CASE ( 2 )    ! Beckmann and Goosse parametrisation 
            stbl(:,:)   = soce
            CALL sbc_isf_bg03(kt)
            !
         CASE ( 3 )    ! specified runoff in depth (Mathiot et al., XXXX in preparation)
            ! specified runoff in depth (Mathiot et al., XXXX in preparation)
            IF( .NOT.l_isfcpl ) THEN
               CALL fld_read ( kt, nn_fsbc, sf_rnfisf   )
               fwfisf(:,:) = - sf_rnfisf(1)%fnow(:,:,1)         ! fresh water flux from the isf (fwfisf <0 mean melting) 
            ENDIF
            qisf(:,:)   = fwfisf(:,:) * rLfusisf             ! heat flux
            stbl(:,:)   = soce
            !
         CASE ( 4 )    ! specified fwf and heat flux forcing beneath the ice shelf
            !          ! specified fwf and heat flux forcing beneath the ice shelf
            IF( .NOT.l_isfcpl ) THEN
               CALL fld_read ( kt, nn_fsbc, sf_fwfisf   )
               !CALL fld_read ( kt, nn_fsbc, sf_qisf   )
               fwfisf(:,:) = -sf_fwfisf(1)%fnow(:,:,1)            ! fwf
            ENDIF
            qisf(:,:)   = fwfisf(:,:) * rLfusisf               ! heat flux
            stbl(:,:)   = soce
            !
         END SELECT

         ! compute tsc due to isf
         ! isf melting implemented as a volume flux and we assume that melt water is at 0 PSU.
         ! WARNING water add at temp = 0C, need to add a correction term (fwfisf * tfreez / rau0).
         ! compute freezing point beneath ice shelf (or top cell if nn_isf = 3)
         DO jj = 1,jpj
            DO ji = 1,jpi
               zdep(ji,jj)=gdepw_n(ji,jj,misfkt(ji,jj))
            END DO
         END DO
         CALL eos_fzp( stbl(:,:), zt_frz(:,:), zdep(:,:) )
         
         risf_tsc(:,:,jp_tem) = qisf(:,:) * r1_rau0_rcp - fwfisf(:,:) * zt_frz(:,:) * r1_rau0 !
         risf_tsc(:,:,jp_sal) = 0.0_wp

         ! lbclnk
         CALL lbc_lnk_multi( 'sbcisf', risf_tsc(:,:,jp_tem), 'T', 1., risf_tsc(:,:,jp_sal), 'T', 1., fwfisf,'T', 1., qisf, 'T', 1.)
         ! output
         IF( iom_use('iceshelf_cea') )   CALL iom_put( 'iceshelf_cea', -fwfisf(:,:)                      )   ! isf mass flux
         IF( iom_use('hflx_isf_cea') )   CALL iom_put( 'hflx_isf_cea', risf_tsc(:,:,jp_tem) * rau0 * rcp )   ! isf sensible+latent heat (W/m2)
         IF( iom_use('qlatisf' ) )       CALL iom_put( 'qlatisf'     , qisf(:,:)                         )   ! isf latent heat
         IF( iom_use('fwfisf'  ) )       CALL iom_put( 'fwfisf'      , fwfisf(:,:)                       )   ! isf mass flux (opposite sign)

         ! Diagnostics
         IF( iom_use('fwfisf3d') .OR. iom_use('qlatisf3d') .OR. iom_use('qhcisf3d') .OR. iom_use('qhcisf')) THEN
            ALLOCATE( zfwfisf3d(jpi,jpj,jpk) , zqhcisf3d(jpi,jpj,jpk) , zqlatisf3d(jpi,jpj,jpk) )
            ALLOCATE( zqhcisf2d(jpi,jpj) )
            !
            zfwfisf3d (:,:,:) = 0._wp                         ! 3d ice shelf melting (kg/m2/s)
            zqhcisf3d (:,:,:) = 0._wp                         ! 3d heat content flux (W/m2)
            zqlatisf3d(:,:,:) = 0._wp                         ! 3d ice shelf melting latent heat flux (W/m2)
            zqhcisf2d (:,:)   = fwfisf(:,:) * zt_frz * rcp    ! 2d heat content flux (W/m2)
            !
            DO jj = 1,jpj
               DO ji = 1,jpi
                  ikt = misfkt(ji,jj)
                  ikb = misfkb(ji,jj)
                  DO jk = ikt, ikb - 1
                     zfwfisf3d (ji,jj,jk) = zfwfisf3d (ji,jj,jk) + fwfisf   (ji,jj) * r1_hisf_tbl(ji,jj) * e3t_n(ji,jj,jk)
                     zqhcisf3d (ji,jj,jk) = zqhcisf3d (ji,jj,jk) + zqhcisf2d(ji,jj) * r1_hisf_tbl(ji,jj) * e3t_n(ji,jj,jk)
                     zqlatisf3d(ji,jj,jk) = zqlatisf3d(ji,jj,jk) + qisf     (ji,jj) * r1_hisf_tbl(ji,jj) * e3t_n(ji,jj,jk)
                  END DO
                  zfwfisf3d (ji,jj,jk) = zfwfisf3d (ji,jj,jk) + fwfisf   (ji,jj) * r1_hisf_tbl(ji,jj)   & 
                     &                                                                   * ralpha(ji,jj) * e3t_n(ji,jj,jk)
                  zqhcisf3d (ji,jj,jk) = zqhcisf3d (ji,jj,jk) + zqhcisf2d(ji,jj) * r1_hisf_tbl(ji,jj)   & 
                     &                                                                   * ralpha(ji,jj) * e3t_n(ji,jj,jk)
                  zqlatisf3d(ji,jj,jk) = zqlatisf3d(ji,jj,jk) + qisf     (ji,jj) * r1_hisf_tbl(ji,jj)   &  
                     &                                                                   * ralpha(ji,jj) * e3t_n(ji,jj,jk)
               END DO
            END DO
            !
            CALL iom_put('fwfisf3d' , zfwfisf3d (:,:,:))
            CALL iom_put('qlatisf3d', zqlatisf3d(:,:,:))
            CALL iom_put('qhcisf3d' , zqhcisf3d (:,:,:))
            CALL iom_put('qhcisf'   , zqhcisf2d (:,:  ))
            !
            DEALLOCATE( zfwfisf3d, zqhcisf3d, zqlatisf3d )
            DEALLOCATE( zqhcisf2d )
         ENDIF
         !
      ENDIF

      IF( kt == nit000 ) THEN                         !   set the forcing field at nit000 - 1    !
         IF( ln_rstart .AND.    &                     ! Restart: read in restart file
            &   iom_varid( numror, 'fwf_isf_b', ldstop = .FALSE. ) > 0 ) THEN
            IF(lwp) WRITE(numout,*) '          nit000-1 isf tracer content forcing fields read in the restart file'
            CALL iom_get( numror, jpdom_autoglo, 'fwf_isf_b', fwfisf_b(:,:)         , ldxios = lrxios )   ! before salt content isf_tsc trend
            CALL iom_get( numror, jpdom_autoglo, 'isf_sc_b' , risf_tsc_b(:,:,jp_sal), ldxios = lrxios )   ! before salt content isf_tsc trend
            CALL iom_get( numror, jpdom_autoglo, 'isf_hc_b' , risf_tsc_b(:,:,jp_tem), ldxios = lrxios )   ! before salt content isf_tsc trend
         ELSE
            fwfisf_b(:,:)    = fwfisf(:,:)
            risf_tsc_b(:,:,:)= risf_tsc(:,:,:)
         ENDIF
      ENDIF
      ! 
      IF( lrst_oce ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbc : isf surface tracer content forcing fields written in ocean restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         IF( lwxios ) CALL iom_swap(      cwxios_context          )
         CALL iom_rstput( kt, nitrst, numrow, 'fwf_isf_b', fwfisf(:,:)         , ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'isf_hc_b' , risf_tsc(:,:,jp_tem), ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'isf_sc_b' , risf_tsc(:,:,jp_sal), ldxios = lwxios )
         IF( lwxios ) CALL iom_swap(      cxios_context          )
      ENDIF
      !
   END SUBROUTINE sbc_isf


   INTEGER FUNCTION sbc_isf_alloc()
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION sbc_isf_rnf_alloc  ***
      !!----------------------------------------------------------------------
      sbc_isf_alloc = 0       ! set to zero if no array to be allocated
      IF( .NOT. ALLOCATED( qisf ) ) THEN
         ALLOCATE(  risf_tsc(jpi,jpj,jpts), risf_tsc_b(jpi,jpj,jpts), qisf(jpi,jpj)   , &
               &    rhisf_tbl(jpi,jpj)    , r1_hisf_tbl(jpi,jpj), rzisf_tbl(jpi,jpj)  , &
               &    ttbl(jpi,jpj)         , stbl(jpi,jpj)       , utbl(jpi,jpj)       , &
               &    vtbl(jpi, jpj)        , risfLeff(jpi,jpj)   , rhisf_tbl_0(jpi,jpj), &
               &    ralpha(jpi,jpj)       , misfkt(jpi,jpj)     , misfkb(jpi,jpj)     , &
               &    STAT= sbc_isf_alloc )
         !
         CALL mpp_sum ( 'sbcisf', sbc_isf_alloc )
         IF( sbc_isf_alloc /= 0 )   CALL ctl_stop( 'STOP', 'sbc_isf_alloc: failed to allocate arrays.' )
         !
      ENDIF
   END FUNCTION


  SUBROUTINE sbc_isf_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_isf_init  ***
      !!
      !! ** Purpose : Initialisation of variables for iceshelf fluxes formulation
      !!
      !! ** Method  :  4 parameterizations are available according to nn_isf 
      !!               nn_isf = 1 : Realistic ice_shelf formulation
      !!                        2 : Beckmann & Goose parameterization
      !!                        3 : Specified runoff in deptht (Mathiot & al. )
      !!                        4 : specified fwf and heat flux forcing beneath the ice shelf
      !!----------------------------------------------------------------------
      INTEGER               :: ji, jj, jk           ! loop index
      INTEGER               :: ik                   ! current level index
      INTEGER               :: ikt, ikb             ! top and bottom level of the isf boundary layer
      INTEGER               :: inum, ierror
      INTEGER               :: ios                  ! Local integer output status for namelist read
      REAL(wp)              :: zhk
      CHARACTER(len=256)    :: cvarzisf, cvarhisf   ! name for isf file
      CHARACTER(LEN=32 )    :: cvarLeff             ! variable name for efficient Length scale
      !!----------------------------------------------------------------------
      NAMELIST/namsbc_isf/ nn_isfblk, rn_hisf_tbl, rn_gammat0, rn_gammas0, nn_gammablk, nn_isf, &
                         & sn_fwfisf, sn_rnfisf, sn_depmax_isf, sn_depmin_isf, sn_Leff_isf
      !!----------------------------------------------------------------------

      REWIND( numnam_ref )              ! Namelist namsbc_rnf in reference namelist : Runoffs 
      READ  ( numnam_ref, namsbc_isf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_isf in reference namelist' )

      REWIND( numnam_cfg )              ! Namelist namsbc_rnf in configuration namelist : Runoffs
      READ  ( numnam_cfg, namsbc_isf, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_isf in configuration namelist' )
      IF(lwm) WRITE ( numond, namsbc_isf )

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'sbc_isf_init : heat flux of the ice shelf'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) '   Namelist namsbc_isf :'
      IF(lwp) WRITE(numout,*) '      type ice shelf melting/freezing         nn_isf      = ', nn_isf
      IF(lwp) WRITE(numout,*) '      bulk formulation (nn_isf=1 only)        nn_isfblk   = ', nn_isfblk
      IF(lwp) WRITE(numout,*) '      thickness of the top boundary layer     rn_hisf_tbl = ', rn_hisf_tbl
      IF(lwp) WRITE(numout,*) '      gamma formulation                       nn_gammablk = ', nn_gammablk 
      IF(lwp) WRITE(numout,*) '      gammat coefficient                      rn_gammat0  = ', rn_gammat0  
      IF(lwp) WRITE(numout,*) '      gammas coefficient                      rn_gammas0  = ', rn_gammas0  
      IF(lwp) WRITE(numout,*) '      top drag coef. used (from namdrg_top)   rn_Cd0      = ', r_Cdmin_top 


                           !  1 = presence of ISF    2 = bg03 parametrisation 
                           !  3 = rnf file for isf   4 = ISF fwf specified
                           !  option 1 and 4 need ln_isfcav = .true. (domzgr)
      !
      ! Allocate public variable
      IF ( sbc_isf_alloc()  /= 0 )         CALL ctl_stop( 'STOP', 'sbc_isf : unable to allocate arrays' )
      !
      ! initialisation
      qisf    (:,:)    = 0._wp   ;   fwfisf  (:,:) = 0._wp
      risf_tsc(:,:,:)  = 0._wp   ;   fwfisf_b(:,:) = 0._wp
      !
      ! define isf tbl tickness, top and bottom indice
      SELECT CASE ( nn_isf )
      CASE ( 1 ) 
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '      ==>>>   presence of under iceshelf seas (nn_isf = 1)'
         rhisf_tbl(:,:) = rn_hisf_tbl
         misfkt   (:,:) = mikt(:,:)         ! same indice for bg03 et cav => used in isfdiv
         !
      CASE ( 2 , 3 )
         IF( .NOT.l_isfcpl ) THEN
             ALLOCATE( sf_rnfisf(1), STAT=ierror )
             ALLOCATE( sf_rnfisf(1)%fnow(jpi,jpj,1), sf_rnfisf(1)%fdta(jpi,jpj,1,2) )
             CALL fld_fill( sf_rnfisf, (/ sn_rnfisf /), cn_dirisf, 'sbc_isf_init', 'read fresh water flux isf data', 'namsbc_isf' )
          ENDIF
          !  read effective lenght (BG03)
          IF( nn_isf == 2 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '      ==>>>   bg03 parametrisation (nn_isf = 2)'
            CALL iom_open( sn_Leff_isf%clname, inum )
            cvarLeff = TRIM(sn_Leff_isf%clvar)
            CALL iom_get( inum, jpdom_data, cvarLeff, risfLeff , 1)
            CALL iom_close(inum)
            !
            risfLeff = risfLeff*1000.0_wp           !: convertion in m
         ELSE
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '      ==>>>   rnf file for isf (nn_isf = 3)'
         ENDIF
         ! read depth of the top and bottom of the isf top boundary layer (in this case, isf front depth and grounding line depth)
         CALL iom_open( sn_depmax_isf%clname, inum )
         cvarhisf = TRIM(sn_depmax_isf%clvar)
         CALL iom_get( inum, jpdom_data, cvarhisf, rhisf_tbl, 1) !: depth of deepest point of the ice shelf base
         CALL iom_close(inum)
         !
         CALL iom_open( sn_depmin_isf%clname, inum )
         cvarzisf = TRIM(sn_depmin_isf%clvar)
         CALL iom_get( inum, jpdom_data, cvarzisf, rzisf_tbl, 1) !: depth of shallowest point of the ice shelves base
         CALL iom_close(inum)
         !
         rhisf_tbl(:,:) = rhisf_tbl(:,:) - rzisf_tbl(:,:)        !: tickness isf boundary layer

         !! compute first level of the top boundary layer
         DO ji = 1, jpi
            DO jj = 1, jpj
                ik = 2
!!gm potential bug: use gdepw_0 not _n
                DO WHILE ( ik <= mbkt(ji,jj) .AND. gdepw_n(ji,jj,ik) < rzisf_tbl(ji,jj) ) ;  ik = ik + 1 ;  END DO
                misfkt(ji,jj) = ik-1
            END DO
         END DO
         !
      CASE ( 4 ) 
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '      ==>>>   specified fresh water flux in ISF (nn_isf = 4)'
         ! as in nn_isf == 1
         rhisf_tbl(:,:) = rn_hisf_tbl
         misfkt   (:,:) = mikt(:,:)         ! same indice for bg03 et cav => used in isfdiv
         !
         ! load variable used in fldread (use for temporal interpolation of isf fwf forcing)
         IF( .NOT.l_isfcpl ) THEN
           ALLOCATE( sf_fwfisf(1), STAT=ierror )
           ALLOCATE( sf_fwfisf(1)%fnow(jpi,jpj,1), sf_fwfisf(1)%fdta(jpi,jpj,1,2) )
           CALL fld_fill( sf_fwfisf, (/ sn_fwfisf /), cn_dirisf, 'sbc_isf_init', 'read fresh water flux isf data', 'namsbc_isf' )
         ENDIF
         !
      CASE DEFAULT
         CALL ctl_stop( 'sbc_isf_init: wrong value of nn_isf' )
      END SELECT
         
      rhisf_tbl_0(:,:) = rhisf_tbl(:,:)

      ! compute bottom level of isf tbl and thickness of tbl below the ice shelf
      DO jj = 1,jpj
         DO ji = 1,jpi
            ikt = misfkt(ji,jj)
            ikb = misfkt(ji,jj)
            ! thickness of boundary layer at least the top level thickness
            rhisf_tbl(ji,jj) = MAX(rhisf_tbl_0(ji,jj), e3t_n(ji,jj,ikt))

            ! determine the deepest level influenced by the boundary layer
            DO jk = ikt+1, mbkt(ji,jj)
               IF( (SUM(e3t_n(ji,jj,ikt:jk-1)) < rhisf_tbl(ji,jj)) .AND. (tmask(ji,jj,jk) == 1) )   ikb = jk
            END DO
            rhisf_tbl(ji,jj) = MIN(rhisf_tbl(ji,jj), SUM(e3t_n(ji,jj,ikt:ikb))) ! limit the tbl to water thickness.
            misfkb(ji,jj) = ikb                                                   ! last wet level of the tbl
            r1_hisf_tbl(ji,jj) = 1._wp / rhisf_tbl(ji,jj)

            zhk           = SUM( e3t_n(ji, jj, ikt:ikb - 1)) * r1_hisf_tbl(ji,jj) ! proportion of tbl cover by cell from ikt to ikb - 1
            ralpha(ji,jj) = rhisf_tbl(ji,jj) * (1._wp - zhk ) / e3t_n(ji,jj,ikb)  ! proportion of bottom cell influenced by boundary layer
         END DO
      END DO

      IF( lwxios ) THEN
          CALL iom_set_rstw_var_active('fwf_isf_b')
          CALL iom_set_rstw_var_active('isf_hc_b')
          CALL iom_set_rstw_var_active('isf_sc_b')
      ENDIF


  END SUBROUTINE sbc_isf_init


  SUBROUTINE sbc_isf_bg03(kt)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_isf_bg03  ***
      !!
      !! ** Purpose : add net heat and fresh water flux from ice shelf melting
      !!          into the adjacent ocean
      !!
      !! ** Method  :   See reference
      !!
      !! ** Reference : Beckmann and Goosse (2003), "A parameterization of ice shelf-ocean
      !!         interaction for climate models", Ocean Modelling 5(2003) 157-170.
      !!         (hereafter BG)
      !! History :  06-02  (C. Wang) Original code
      !!----------------------------------------------------------------------
      INTEGER, INTENT ( in ) :: kt
      !
      INTEGER  :: ji, jj, jk ! dummy loop index
      INTEGER  :: ik         ! current level
      REAL(wp) :: zt_sum     ! sum of the temperature between 200m and 600m
      REAL(wp) :: zt_ave     ! averaged temperature between 200m and 600m
      REAL(wp) :: zt_frz     ! freezing point temperature at depth z
      REAL(wp) :: zpress     ! pressure to compute the freezing point in depth
      !!----------------------------------------------------------------------
      !
      DO ji = 1, jpi
         DO jj = 1, jpj
            ik = misfkt(ji,jj)
            !! Initialize arrays to 0 (each step)
            zt_sum = 0.e0_wp
            IF ( ik > 1 ) THEN
               ! 1. -----------the average temperature between 200m and 600m ---------------------
               DO jk = misfkt(ji,jj),misfkb(ji,jj)
                  ! Calculate freezing temperature
                  zpress = grav*rau0*gdept_n(ji,jj,ik)*1.e-04
                  CALL eos_fzp(stbl(ji,jj), zt_frz, zpress) 
                  zt_sum = zt_sum + (tsn(ji,jj,jk,jp_tem)-zt_frz) * e3t_n(ji,jj,jk) * tmask(ji,jj,jk)  ! sum temp
               END DO
               zt_ave = zt_sum/rhisf_tbl(ji,jj) ! calcul mean value
               ! 2. ------------Net heat flux and fresh water flux due to the ice shelf
               ! For those corresponding to zonal boundary    
               qisf(ji,jj) = - rau0 * rcp * rn_gammat0 * risfLeff(ji,jj) * e1t(ji,jj) * zt_ave  &
                           & * r1_e1e2t(ji,jj) * tmask(ji,jj,jk)
             
               fwfisf(ji,jj) = qisf(ji,jj) / rLfusisf          !fresh water flux kg/(m2s)                  
               fwfisf(ji,jj) = fwfisf(ji,jj) * ( soce / stbl(ji,jj) )
               !add to salinity trend
            ELSE
               qisf(ji,jj) = 0._wp   ;   fwfisf(ji,jj) = 0._wp
            END IF
         END DO
      END DO
      !
  END SUBROUTINE sbc_isf_bg03


  SUBROUTINE sbc_isf_cav( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_isf_cav  ***
      !!
      !! ** Purpose :   handle surface boundary condition under ice shelf
      !!
      !! ** Method  : -
      !!
      !! ** Action  :   utau, vtau : remain unchanged
      !!                taum, wndm : remain unchanged
      !!                qns        : update heat flux below ice shelf
      !!                emp, emps  : update freshwater flux below ice shelf
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      INTEGER  ::   nit
      LOGICAL  ::   lit
      REAL(wp) ::   zlamb1, zlamb2, zlamb3
      REAL(wp) ::   zeps1,zeps2,zeps3,zeps4,zeps6,zeps7
      REAL(wp) ::   zaqe,zbqe,zcqe,zaqer,zdis,zsfrz,zcfac
      REAL(wp) ::   zeps = 1.e-20_wp        
      REAL(wp) ::   zerr
      REAL(wp), DIMENSION(jpi,jpj) ::   zfrz
      REAL(wp), DIMENSION(jpi,jpj) ::   zgammat, zgammas 
      REAL(wp), DIMENSION(jpi,jpj) ::   zfwflx, zhtflx, zhtflx_b
      !!---------------------------------------------------------------------
      !
      ! coeficient for linearisation of potential tfreez
      ! Crude approximation for pressure (but commonly used)
      IF ( l_useCT ) THEN   ! linearisation from Jourdain et al. (2017)
         zlamb1 =-0.0564_wp
         zlamb2 = 0.0773_wp
         zlamb3 =-7.8633e-8 * grav * rau0
      ELSE                  ! linearisation from table 4 (Asay-Davis et al., 2015)
         zlamb1 =-0.0573_wp
         zlamb2 = 0.0832_wp
         zlamb3 =-7.53e-8 * grav * rau0
      ENDIF
      !
      ! initialisation
      zgammat(:,:) = rn_gammat0 ; zgammas (:,:) = rn_gammas0
      zhtflx (:,:) = 0.0_wp     ; zhtflx_b(:,:) = 0.0_wp    
      zfwflx (:,:) = 0.0_wp

      ! compute ice shelf melting
      nit = 1 ; lit = .TRUE.
      DO WHILE ( lit )    ! maybe just a constant number of iteration as in blk_core is fine
         SELECT CASE ( nn_isfblk )
         CASE ( 1 )   !  ISOMIP formulation (2 equations) for volume flux (Hunter et al., 2006)
            ! Calculate freezing temperature
            CALL eos_fzp( stbl(:,:), zfrz(:,:), risfdep(:,:) )

            ! compute gammat every where (2d)
            CALL sbc_isf_gammats(zgammat, zgammas, zhtflx, zfwflx)
            
            ! compute upward heat flux zhtflx and upward water flux zwflx
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zhtflx(ji,jj) =   zgammat(ji,jj)*rcp*rau0*(ttbl(ji,jj)-zfrz(ji,jj))
                  zfwflx(ji,jj) = - zhtflx(ji,jj)/rLfusisf
               END DO
            END DO

            ! Compute heat flux and upward fresh water flux
            qisf  (:,:) = - zhtflx(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:)
            fwfisf(:,:) =   zfwflx(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:)

         CASE ( 2 )  ! ISOMIP+ formulation (3 equations) for volume flux (Asay-Davis et al., 2015)
            ! compute gammat every where (2d)
            CALL sbc_isf_gammats(zgammat, zgammas, zhtflx, zfwflx)

            ! compute upward heat flux zhtflx and upward water flux zwflx
            ! Resolution of a 2d equation from equation 21, 22 and 23 to find Sb (Asay-Davis et al., 2015)
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! compute coeficient to solve the 2nd order equation
                  zeps1 = rcp*rau0*zgammat(ji,jj)
                  zeps2 = rLfusisf*rau0*zgammas(ji,jj)
                  zeps3 = rhoisf*rcpisf*rkappa/MAX(risfdep(ji,jj),zeps)
                  zeps4 = zlamb2+zlamb3*risfdep(ji,jj)
                  zeps6 = zeps4-ttbl(ji,jj)
                  zeps7 = zeps4-tsurf
                  zaqe  = zlamb1 * (zeps1 + zeps3)
                  zaqer = 0.5_wp/MIN(zaqe,-zeps)
                  zbqe  = zeps1*zeps6+zeps3*zeps7-zeps2
                  zcqe  = zeps2*stbl(ji,jj)
                  zdis  = zbqe*zbqe-4.0_wp*zaqe*zcqe               

                  ! Presumably zdis can never be negative because gammas is very small compared to gammat
                  ! compute s freeze
                  zsfrz=(-zbqe-SQRT(zdis))*zaqer
                  IF ( zsfrz < 0.0_wp ) zsfrz=(-zbqe+SQRT(zdis))*zaqer

                  ! compute t freeze (eq. 22)
                  zfrz(ji,jj)=zeps4+zlamb1*zsfrz
  
                  ! zfwflx is upward water flux
                  ! zhtflx is upward heat flux (out of ocean)
                  ! compute the upward water and heat flux (eq. 28 and eq. 29)
                  zfwflx(ji,jj) = rau0 * zgammas(ji,jj) * (zsfrz-stbl(ji,jj)) / MAX(zsfrz,zeps)
                  zhtflx(ji,jj) = zgammat(ji,jj) * rau0 * rcp * (ttbl(ji,jj) - zfrz(ji,jj) ) 
               END DO
            END DO

            ! compute heat and water flux
            qisf  (:,:) = - zhtflx(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:)
            fwfisf(:,:) =   zfwflx(:,:) * (1._wp - tmask(:,:,1)) * ssmask(:,:)

         END SELECT

         ! define if we need to iterate (nn_gammablk 0/1 do not need iteration)
         IF ( nn_gammablk <  2 ) THEN ; lit = .FALSE.
         ELSE                           
            ! check total number of iteration
            IF (nit >= 100) THEN ; CALL ctl_stop( 'STOP', 'sbc_isf_hol99 : too many iteration ...' )
            ELSE                 ; nit = nit + 1
            END IF

            ! compute error between 2 iterations
            ! if needed save gammat and compute zhtflx_b for next iteration
            zerr = MAXVAL(ABS(zhtflx-zhtflx_b))
            IF ( zerr <= 0.01_wp ) THEN ; lit = .FALSE.
            ELSE                        ; zhtflx_b(:,:) = zhtflx(:,:)
            END IF
         END IF
      END DO
      !
      CALL iom_put('isfgammat', zgammat)
      CALL iom_put('isfgammas', zgammas)
      ! 
   END SUBROUTINE sbc_isf_cav


   SUBROUTINE sbc_isf_gammats(pgt, pgs, pqhisf, pqwisf )
      !!----------------------------------------------------------------------
      !! ** Purpose    : compute the coefficient echange for heat flux
      !!
      !! ** Method     : gamma assume constant or depends of u* and stability
      !!
      !! ** References : Holland and Jenkins, 1999, JPO, p1787-1800, eq 14
      !!                Jenkins et al., 2010, JPO, p2298-2312
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   pgt   , pgs      ! 
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   pqhisf, pqwisf   ! 
      !
      INTEGER  :: ji, jj                     ! loop index
      INTEGER  :: ikt                        ! local integer
      REAL(wp) :: zdku, zdkv                 ! U, V shear 
      REAL(wp) :: zPr, zSc, zRc              ! Prandtl, Scmidth and Richardson number 
      REAL(wp) :: zmob, zmols                ! Monin Obukov length, coriolis factor at T point
      REAL(wp) :: zbuofdep, zhnu             ! Bouyancy length scale, sublayer tickness
      REAL(wp) :: zhmax                      ! limitation of mol
      REAL(wp) :: zetastar                   ! stability parameter
      REAL(wp) :: zgmolet, zgmoles, zgturb   ! contribution of modelecular sublayer and turbulence 
      REAL(wp) :: zcoef                      ! temporary coef
      REAL(wp) :: zdep
      REAL(wp) :: zeps = 1.0e-20_wp    
      REAL(wp), PARAMETER :: zxsiN = 0.052_wp   ! dimensionless constant
      REAL(wp), PARAMETER :: znu   = 1.95e-6_wp ! kinamatic viscosity of sea water (m2.s-1)
      REAL(wp), DIMENSION(2) :: zts, zab
      REAL(wp), DIMENSION(jpi,jpj) :: zustar   ! U, V at T point and friction velocity
      !!---------------------------------------------------------------------
      !
      SELECT CASE ( nn_gammablk )
      CASE ( 0 ) ! gamma is constant (specified in namelist)
         !! ISOMIP formulation (Hunter et al, 2006)
         pgt(:,:) = rn_gammat0
         pgs(:,:) = rn_gammas0

      CASE ( 1 ) ! gamma is assume to be proportional to u*
         !! Jenkins et al., 2010, JPO, p2298-2312
         !! Adopted by Asay-Davis et al. (2015)
         !! compute ustar (eq. 24)
!!gm  NB  use pCdU here so that it will incorporate local boost of Cd0 and log layer case :
!!         zustar(:,:) = SQRT( rCdU_top(:,:) * SQRT(utbl(:,:) * utbl(:,:) + vtbl(:,:) * vtbl(:,:) + r_ke0_top) )
!! or better :  compute ustar in zdfdrg  and use it here as well as in TKE, GLS and Co
!!
!!     ===>>>>    GM  to be done this chrismas
!!         
!!gm end
         zustar(:,:) = SQRT( r_Cdmin_top * (utbl(:,:) * utbl(:,:) + vtbl(:,:) * vtbl(:,:) + r_ke0_top) )

         !! Compute gammats
         pgt(:,:) = zustar(:,:) * rn_gammat0
         pgs(:,:) = zustar(:,:) * rn_gammas0
      
      CASE ( 2 ) ! gamma depends of stability of boundary layer
         !! Holland and Jenkins, 1999, JPO, p1787-1800, eq 14
         !! as MOL depends of flux and flux depends of MOL, best will be iteration (TO DO)
         !! compute ustar
         zustar(:,:) = SQRT( r_Cdmin_top * (utbl(:,:) * utbl(:,:) + vtbl(:,:) * vtbl(:,:) + r_ke0_top) )

         !! compute Pr and Sc number (can be improved)
         zPr =   13.8_wp
         zSc = 2432.0_wp

         !! compute gamma mole
         zgmolet = 12.5_wp * zPr ** (2.0/3.0) - 6.0_wp
         zgmoles = 12.5_wp * zSc ** (2.0/3.0) - 6.0_wp

         !! compute gamma
         DO ji = 2, jpi
            DO jj = 2, jpj
               ikt = mikt(ji,jj)

               IF( zustar(ji,jj) == 0._wp ) THEN           ! only for kt = 1 I think
                  pgt = rn_gammat0
                  pgs = rn_gammas0
               ELSE
                  !! compute Rc number (as done in zdfric.F90)
!!gm better to do it like in the new zdfric.F90   i.e. avm weighted Ri computation
!!gm moreover, use Max(rn2,0) to take care of static instabilities....
                  zcoef = 0.5_wp / e3w_n(ji,jj,ikt+1)
                  !                                            ! shear of horizontal velocity
                  zdku = zcoef * (  un(ji-1,jj  ,ikt  ) + un(ji,jj,ikt  )  &
                     &             -un(ji-1,jj  ,ikt+1) - un(ji,jj,ikt+1)  )
                  zdkv = zcoef * (  vn(ji  ,jj-1,ikt  ) + vn(ji,jj,ikt  )  &
                     &             -vn(ji  ,jj-1,ikt+1) - vn(ji,jj,ikt+1)  )
                  !                                            ! richardson number (minimum value set to zero)
                  zRc = rn2(ji,jj,ikt+1) / MAX( zdku*zdku + zdkv*zdkv, zeps )

                  !! compute bouyancy 
                  zts(jp_tem) = ttbl(ji,jj)
                  zts(jp_sal) = stbl(ji,jj)
                  zdep        = gdepw_n(ji,jj,ikt)
                  !
                  CALL eos_rab( zts, zdep, zab )
                  !
                  !! compute length scale 
                  zbuofdep = grav * ( zab(jp_tem) * pqhisf(ji,jj) - zab(jp_sal) * pqwisf(ji,jj) )  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

                  !! compute Monin Obukov Length
                  ! Maximum boundary layer depth
                  zhmax = gdept_n(ji,jj,mbkt(ji,jj)) - gdepw_n(ji,jj,mikt(ji,jj)) -0.001_wp
                  ! Compute Monin obukhov length scale at the surface and Ekman depth:
                  zmob   = zustar(ji,jj) ** 3 / (vkarmn * (zbuofdep + zeps))
                  zmols  = SIGN(1._wp, zmob) * MIN(ABS(zmob), zhmax) * tmask(ji,jj,ikt)

                  !! compute eta* (stability parameter)
                  zetastar = 1._wp / ( SQRT(1._wp + MAX(zxsiN * zustar(ji,jj) / ( ABS(ff_f(ji,jj)) * zmols * zRc ), 0._wp)))

                  !! compute the sublayer thickness
                  zhnu = 5 * znu / zustar(ji,jj)

                  !! compute gamma turb
                  zgturb = 1._wp / vkarmn * LOG(zustar(ji,jj) * zxsiN * zetastar * zetastar / ( ABS(ff_f(ji,jj)) * zhnu )) &
                  &      + 1._wp / ( 2 * zxsiN * zetastar ) - 1._wp / vkarmn

                  !! compute gammats
                  pgt(ji,jj) = zustar(ji,jj) / (zgturb + zgmolet)
                  pgs(ji,jj) = zustar(ji,jj) / (zgturb + zgmoles)
               END IF
            END DO
         END DO
         CALL lbc_lnk_multi( 'sbcisf', pgt, 'T', 1., pgs, 'T', 1.)
      END SELECT
      !
   END SUBROUTINE sbc_isf_gammats


   SUBROUTINE sbc_isf_tbl( pvarin, pvarout, cd_ptin )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE sbc_isf_tbl  ***
      !!
      !! ** Purpose : compute mean T/S/U/V in the boundary layer at T- point
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) :: pvarin
      REAL(wp), DIMENSION(:,:)  , INTENT(  out) :: pvarout
      CHARACTER(len=1),           INTENT(in   ) :: cd_ptin ! point of variable in/out
      !
      INTEGER ::   ji, jj, jk                   ! loop index
      INTEGER ::   ikt, ikb                     ! top and bottom index of the tbl
      REAL(wp) ::   ze3, zhk
      REAL(wp), DIMENSION(jpi,jpj) :: zhisf_tbl ! thickness of the tbl
      REAL(wp), DIMENSION(jpi,jpj) :: zvarout
      !!----------------------------------------------------------------------
      
      ! initialisation
      pvarout(:,:)=0._wp
   
      SELECT CASE ( cd_ptin )
      CASE ( 'U' ) ! compute U in the top boundary layer at T- point 
         !
         zvarout(:,:)=0._wp
         !
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = miku(ji,jj) ; ikb = miku(ji,jj)
               ! thickness of boundary layer at least the top level thickness
               zhisf_tbl(ji,jj) = MAX( rhisf_tbl_0(ji,jj) , e3u_n(ji,jj,ikt) )

               ! determine the deepest level influenced by the boundary layer
               DO jk = ikt+1, mbku(ji,jj)
                  IF ( (SUM(e3u_n(ji,jj,ikt:jk-1)) < zhisf_tbl(ji,jj)) .AND. (umask(ji,jj,jk) == 1) ) ikb = jk
               END DO
               zhisf_tbl(ji,jj) = MIN(zhisf_tbl(ji,jj), SUM(e3u_n(ji,jj,ikt:ikb)))  ! limit the tbl to water thickness.

               ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  ze3 = e3u_n(ji,jj,jk)
                  zvarout(ji,jj) = zvarout(ji,jj) + pvarin(ji,jj,jk) / zhisf_tbl(ji,jj) * ze3
               END DO

               ! level partially include in ice shelf boundary layer 
               zhk = SUM( e3u_n(ji, jj, ikt:ikb - 1)) / zhisf_tbl(ji,jj)
               zvarout(ji,jj) = zvarout(ji,jj) + pvarin(ji,jj,ikb) * (1._wp - zhk)
            END DO
         END DO
         DO jj = 2, jpj
            DO ji = 2, jpi
!!gm a wet-point only average should be used here !!!
               pvarout(ji,jj) = 0.5_wp * (zvarout(ji,jj) + zvarout(ji-1,jj))
            END DO
         END DO
         CALL lbc_lnk('sbcisf', pvarout,'T',-1.)
      
      CASE ( 'V' ) ! compute V in the top boundary layer at T- point 
         !
         zvarout(:,:)=0._wp
         !
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = mikv(ji,jj) ; ikb = mikv(ji,jj)
               ! thickness of boundary layer at least the top level thickness
               zhisf_tbl(ji,jj) = MAX(rhisf_tbl_0(ji,jj), e3v_n(ji,jj,ikt))

               ! determine the deepest level influenced by the boundary layer
               DO jk = ikt+1, mbkv(ji,jj)
                  IF ( (SUM(e3v_n(ji,jj,ikt:jk-1)) < zhisf_tbl(ji,jj)) .AND. (vmask(ji,jj,jk) == 1) ) ikb = jk
               END DO
               zhisf_tbl(ji,jj) = MIN(zhisf_tbl(ji,jj), SUM(e3v_n(ji,jj,ikt:ikb)))  ! limit the tbl to water thickness.

               ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  ze3 = e3v_n(ji,jj,jk)
                  zvarout(ji,jj) = zvarout(ji,jj) + pvarin(ji,jj,jk) / zhisf_tbl(ji,jj) * ze3
               END DO

               ! level partially include in ice shelf boundary layer 
               zhk = SUM( e3v_n(ji, jj, ikt:ikb - 1)) / zhisf_tbl(ji,jj)
               zvarout(ji,jj) = zvarout(ji,jj) + pvarin(ji,jj,ikb) * (1._wp - zhk)
            END DO
         END DO
         DO jj = 2, jpj
            DO ji = 2, jpi
!!gm a wet-point only average should be used here !!!
               pvarout(ji,jj) = 0.5_wp * (zvarout(ji,jj) + zvarout(ji,jj-1))
            END DO
         END DO
         CALL lbc_lnk('sbcisf', pvarout,'T',-1.)

      CASE ( 'T' ) ! compute T in the top boundary layer at T- point 
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = misfkt(ji,jj)
               ikb = misfkb(ji,jj)

               ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  ze3 = e3t_n(ji,jj,jk)
                  pvarout(ji,jj) = pvarout(ji,jj) + pvarin(ji,jj,jk) * r1_hisf_tbl(ji,jj) * ze3
               END DO

               ! level partially include in ice shelf boundary layer 
               zhk = SUM( e3t_n(ji, jj, ikt:ikb - 1)) * r1_hisf_tbl(ji,jj)
               pvarout(ji,jj) = pvarout(ji,jj) + pvarin(ji,jj,ikb) * (1._wp - zhk)
            END DO
         END DO
      END SELECT
      !
      ! mask mean tbl value
      pvarout(:,:) = pvarout(:,:) * ssmask(:,:)
      !
   END SUBROUTINE sbc_isf_tbl
      

   SUBROUTINE sbc_isf_div( phdivn )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE sbc_isf_div  ***
      !!       
      !! ** Purpose :   update the horizontal divergence with the runoff inflow
      !!
      !! ** Method  :   
      !!                CAUTION : risf_tsc(:,:,jp_sal) is negative (outflow) increase the 
      !!                          divergence and expressed in m/s
      !!
      !! ** Action  :   phdivn   decreased by the runoff inflow
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT( inout ) ::   phdivn   ! horizontal divergence
      ! 
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikt, ikb 
      REAL(wp) ::   zhk
      REAL(wp) ::   zfact     ! local scalar
      !!----------------------------------------------------------------------
      !
      zfact   = 0.5_wp
      !
      IF(.NOT.ln_linssh ) THEN     ! need to re compute level distribution of isf fresh water
         DO jj = 1,jpj
            DO ji = 1,jpi
               ikt = misfkt(ji,jj)
               ikb = misfkt(ji,jj)
               ! thickness of boundary layer at least the top level thickness
               rhisf_tbl(ji,jj) = MAX(rhisf_tbl_0(ji,jj), e3t_n(ji,jj,ikt))

               ! determine the deepest level influenced by the boundary layer
               DO jk = ikt, mbkt(ji,jj)
                  IF ( (SUM(e3t_n(ji,jj,ikt:jk-1)) .LT. rhisf_tbl(ji,jj)) .AND. (tmask(ji,jj,jk) == 1) ) ikb = jk
               END DO
               rhisf_tbl(ji,jj) = MIN(rhisf_tbl(ji,jj), SUM(e3t_n(ji,jj,ikt:ikb)))  ! limit the tbl to water thickness.
               misfkb(ji,jj) = ikb                                                  ! last wet level of the tbl
               r1_hisf_tbl(ji,jj) = 1._wp / rhisf_tbl(ji,jj)

               zhk           = SUM( e3t_n(ji, jj, ikt:ikb - 1)) * r1_hisf_tbl(ji,jj)  ! proportion of tbl cover by cell from ikt to ikb - 1
               ralpha(ji,jj) = rhisf_tbl(ji,jj) * (1._wp - zhk ) / e3t_n(ji,jj,ikb)  ! proportion of bottom cell influenced by boundary layer
            END DO
         END DO
      END IF 
      !
      !==   ice shelf melting distributed over several levels   ==!
      DO jj = 1,jpj
         DO ji = 1,jpi
               ikt = misfkt(ji,jj)
               ikb = misfkb(ji,jj)
               ! level fully include in the ice shelf boundary layer
               DO jk = ikt, ikb - 1
                  phdivn(ji,jj,jk) = phdivn(ji,jj,jk) + ( fwfisf(ji,jj) + fwfisf_b(ji,jj) ) &
                    &              * r1_hisf_tbl(ji,jj) * r1_rau0 * zfact
               END DO
               ! level partially include in ice shelf boundary layer 
               phdivn(ji,jj,ikb) = phdivn(ji,jj,ikb) + ( fwfisf(ji,jj) &
                    &            + fwfisf_b(ji,jj) ) * r1_hisf_tbl(ji,jj) * r1_rau0 * zfact * ralpha(ji,jj) 
         END DO
      END DO
      !
   END SUBROUTINE sbc_isf_div

   !!======================================================================
END MODULE sbcisf
