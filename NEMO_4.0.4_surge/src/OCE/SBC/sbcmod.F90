MODULE sbcmod
   !!======================================================================
   !!                       ***  MODULE  sbcmod  ***
   !! Surface module :  provide to the ocean its surface boundary condition
   !!======================================================================
   !! History :  3.0  ! 2006-07  (G. Madec)  Original code
   !!            3.1  ! 2008-08  (S. Masson, A. Caubel, E. Maisonnave, G. Madec) coupled interface
   !!            3.3  ! 2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!            3.3  ! 2010-10  (S. Masson)  add diurnal cycle
   !!            3.3  ! 2010-09  (D. Storkey) add ice boundary conditions (BDY)
   !!             -   ! 2010-11  (G. Madec) ice-ocean stress always computed at each ocean time-step
   !!             -   ! 2010-10  (J. Chanut, C. Bricaud, G. Madec)  add the surface pressure forcing
   !!            3.4  ! 2011-11  (C. Harris) CICE added as an option
   !!            3.5  ! 2012-11  (A. Coward, G. Madec) Rethink of heat, mass and salt surface fluxes
   !!            3.6  ! 2014-11  (P. Mathiot, C. Harris) add ice shelves melting
   !!            4.0  ! 2016-06  (L. Brodeau) new general bulk formulation
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_init      : read namsbc namelist
   !!   sbc           : surface ocean momentum, heat and freshwater boundary conditions
   !!   sbc_final     : Finalize CICE ice model (if used)
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE trc_oce        ! shared ocean-passive tracers variables
   USE sbc_ice        ! Surface boundary condition: ice fields
   USE sbcdcy         ! surface boundary condition: diurnal cycle
   USE sbcssm         ! surface boundary condition: sea-surface mean variables
   USE sbcflx         ! surface boundary condition: flux formulation
   USE sbcblk         ! surface boundary condition: bulk formulation
   USE sbcice_if      ! surface boundary condition: ice-if sea-ice model
#if defined key_si3
   USE icestp         ! surface boundary condition: SI3 sea-ice model
#endif
   USE sbcice_cice    ! surface boundary condition: CICE sea-ice model
   USE sbcisf         ! surface boundary condition: ice-shelf
   USE sbccpl         ! surface boundary condition: coupled formulation
   USE cpl_oasis3     ! OASIS routines for coupling
   USE sbcssr         ! surface boundary condition: sea surface restoring
   USE sbcrnf         ! surface boundary condition: runoffs
   USE sbcapr         ! surface boundary condition: atmo pressure 
   USE sbcisf         ! surface boundary condition: ice shelf
   USE sbcfwb         ! surface boundary condition: freshwater budget
   USE icbstp         ! Icebergs
   USE icb_oce  , ONLY : ln_passive_mode      ! iceberg interaction mode
   USE traqsr         ! active tracers: light penetration
   USE sbcwave        ! Wave module
   USE bdy_oce   , ONLY: ln_bdy
   USE usrdef_sbc     ! user defined: surface boundary condition
   USE closea         ! closed sea
   !
   USE prtctl         ! Print control                    (prt_ctl routine)
   USE iom            ! IOM library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing
   USE wet_dry
   USE diurnal_bulk, ONLY:   ln_diurnal_only   ! diurnal SST diagnostic

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc        ! routine called by step.F90
   PUBLIC   sbc_init   ! routine called by opa.F90

   INTEGER ::   nsbc   ! type of surface boundary condition (deduced from namsbc informations)

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_init
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_init ***
      !!
      !! ** Purpose :   Initialisation of the ocean surface boundary computation
      !!
      !! ** Method  :   Read the namsbc namelist and set derived parameters
      !!                Call init routines for all other SBC modules that have one
      !!
      !! ** Action  : - read namsbc parameters
      !!              - nsbc: type of sbc
      !!----------------------------------------------------------------------
      INTEGER ::   ios, icpt                         ! local integer
      LOGICAL ::   ll_purecpl, ll_opa, ll_not_nemo   ! local logical
      !!
      NAMELIST/namsbc/ nn_fsbc  ,                                                    &
         &             ln_usr   , ln_flx   , ln_blk       ,                          &
         &             ln_cpl   , ln_mixcpl, nn_components,                          &
         &             nn_ice   , ln_ice_embd,                                       &
         &             ln_traqsr, ln_dm2dc ,                                         &
         &             ln_rnf   , nn_fwb   , ln_ssr   , ln_isf    , ln_apr_dyn ,     &
         &             ln_wave  , ln_cdgw  , ln_sdw   , ln_tauwoc  , ln_stcor   ,     &
         &             ln_tauw  , nn_lsm, nn_sdrift
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_init : surface boundary condition setting'
         WRITE(numout,*) '~~~~~~~~ '
      ENDIF
      !
      !                       !**  read Surface Module namelist
      REWIND( numnam_ref )          !* Namelist namsbc in reference namelist : Surface boundary
      READ  ( numnam_ref, namsbc, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc in reference namelist' )
      REWIND( numnam_cfg )          !* Namelist namsbc in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, namsbc, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc in configuration namelist' )
      IF(lwm) WRITE( numond, namsbc )
      !
#if defined key_mpp_mpi
      ncom_fsbc = nn_fsbc    ! make nn_fsbc available for lib_mpp
#endif
      !                             !* overwrite namelist parameter using CPP key information
#if defined key_agrif
      IF( Agrif_Root() ) THEN                ! AGRIF zoom (cf r1242: possibility to run without ice in fine grid)
         IF( lk_si3  )   nn_ice      = 2
         IF( lk_cice )   nn_ice      = 3
      ENDIF
#else
      IF( lk_si3  )   nn_ice      = 2
      IF( lk_cice )   nn_ice      = 3
#endif
      !
#if ! defined key_si3
      IF( nn_ice == 2 )    nn_ice = 0  ! without key key_si3 you cannot use si3...
#endif
      !
      IF(lwp) THEN                  !* Control print
         WRITE(numout,*) '   Namelist namsbc (partly overwritten with CPP key setting)'
         WRITE(numout,*) '      frequency update of sbc (and ice)             nn_fsbc       = ', nn_fsbc
         WRITE(numout,*) '      Type of air-sea fluxes : '
         WRITE(numout,*) '         user defined formulation                   ln_usr        = ', ln_usr
         WRITE(numout,*) '         flux         formulation                   ln_flx        = ', ln_flx
         WRITE(numout,*) '         bulk         formulation                   ln_blk        = ', ln_blk
         WRITE(numout,*) '      Type of coupling (Ocean/Ice/Atmosphere) : '
         WRITE(numout,*) '         ocean-atmosphere coupled formulation       ln_cpl        = ', ln_cpl
         WRITE(numout,*) '         mixed forced-coupled     formulation       ln_mixcpl     = ', ln_mixcpl
!!gm  lk_oasis is controlled by key_oasis3  ===>>>  It shoud be removed from the namelist 
         WRITE(numout,*) '         OASIS coupling (with atm or sas)           lk_oasis      = ', lk_oasis
         WRITE(numout,*) '         components of your executable              nn_components = ', nn_components
         WRITE(numout,*) '      Sea-ice : '
         WRITE(numout,*) '         ice management in the sbc (=0/1/2/3)       nn_ice        = ', nn_ice
         WRITE(numout,*) '         ice embedded into ocean                    ln_ice_embd   = ', ln_ice_embd
         WRITE(numout,*) '      Misc. options of sbc : '
         WRITE(numout,*) '         Light penetration in temperature Eq.       ln_traqsr     = ', ln_traqsr
         WRITE(numout,*) '            daily mean to diurnal cycle qsr            ln_dm2dc   = ', ln_dm2dc
         WRITE(numout,*) '         Sea Surface Restoring on SST and/or SSS    ln_ssr        = ', ln_ssr
         WRITE(numout,*) '         FreshWater Budget control  (=0/1/2)        nn_fwb        = ', nn_fwb
         WRITE(numout,*) '         Patm gradient added in ocean & ice Eqs.    ln_apr_dyn    = ', ln_apr_dyn
         WRITE(numout,*) '         runoff / runoff mouths                     ln_rnf        = ', ln_rnf
         WRITE(numout,*) '         iceshelf formulation                       ln_isf        = ', ln_isf
         WRITE(numout,*) '         nb of iterations if land-sea-mask applied  nn_lsm        = ', nn_lsm
         WRITE(numout,*) '         surface wave                               ln_wave       = ', ln_wave
         WRITE(numout,*) '               Stokes drift corr. to vert. velocity ln_sdw        = ', ln_sdw
         WRITE(numout,*) '                  vertical parametrization          nn_sdrift     = ', nn_sdrift
         WRITE(numout,*) '               wave modified ocean stress           ln_tauwoc     = ', ln_tauwoc
         WRITE(numout,*) '               wave modified ocean stress component ln_tauw       = ', ln_tauw
         WRITE(numout,*) '               Stokes coriolis term                 ln_stcor      = ', ln_stcor
         WRITE(numout,*) '               neutral drag coefficient (CORE,NCAR) ln_cdgw       = ', ln_cdgw
      ENDIF
      !
      IF( .NOT.ln_wave ) THEN
         ln_sdw = .false. ; ln_cdgw = .false. ; ln_tauwoc = .false. ; ln_tauw = .false. ; ln_stcor = .false.
      ENDIF 
      IF( ln_sdw ) THEN
         IF( .NOT.(nn_sdrift==jp_breivik_2014 .OR. nn_sdrift==jp_li_2017 .OR. nn_sdrift==jp_peakfr) ) &
            CALL ctl_stop( 'The chosen nn_sdrift for Stokes drift vertical velocity must be 0, 1, or 2' )
      ENDIF
      ll_st_bv2014  = ( nn_sdrift==jp_breivik_2014 )
      ll_st_li2017  = ( nn_sdrift==jp_li_2017 )
      ll_st_bv_li   = ( ll_st_bv2014 .OR. ll_st_li2017 )
      ll_st_peakfr  = ( nn_sdrift==jp_peakfr )
      IF( ln_tauwoc .AND. ln_tauw ) &
         CALL ctl_stop( 'More than one method for modifying the ocean stress has been selected ', &
                                  '(ln_tauwoc=.true. and ln_tauw=.true.)' )
      IF( ln_tauwoc ) &
         CALL ctl_warn( 'You are subtracting the wave stress to the ocean (ln_tauwoc=.true.)' )
      IF( ln_tauw ) &
         CALL ctl_warn( 'The wave modified ocean stress components are used (ln_tauw=.true.) ', &
                              'This will override any other specification of the ocean stress' )
      !
      IF( .NOT.ln_usr ) THEN     ! the model calendar needs some specificities (except in user defined case)
         IF( MOD( rday , rdt ) /= 0. )   CALL ctl_stop( 'the time step must devide the number of second of in a day' )
         IF( MOD( rday , 2.  ) /= 0. )   CALL ctl_stop( 'the number of second of in a day must be an even number'    )
         IF( MOD( rdt  , 2.  ) /= 0. )   CALL ctl_stop( 'the time step (in second) must be an even number'           )
      ENDIF
      !                       !**  check option consistency
      !
      IF(lwp) WRITE(numout,*)       !* Single / Multi - executable (NEMO / OPA+SAS) 
      SELECT CASE( nn_components )
      CASE( jp_iam_nemo )
         IF(lwp) WRITE(numout,*) '   ==>>>   NEMO configured as a single executable (i.e. including both OPA and Surface module)'
      CASE( jp_iam_opa  )
         IF(lwp) WRITE(numout,*) '   ==>>>   Multi executable configuration. Here, OPA component'
         IF( .NOT.lk_oasis )   CALL ctl_stop( 'sbc_init : OPA-SAS coupled via OASIS, but key_oasis3 disabled' )
         IF( ln_cpl        )   CALL ctl_stop( 'sbc_init : OPA-SAS coupled via OASIS, but ln_cpl = T in OPA'   )
         IF( ln_mixcpl     )   CALL ctl_stop( 'sbc_init : OPA-SAS coupled via OASIS, but ln_mixcpl = T in OPA' )
      CASE( jp_iam_sas  )
         IF(lwp) WRITE(numout,*) '   ==>>>   Multi executable configuration. Here, SAS component'
         IF( .NOT.lk_oasis )   CALL ctl_stop( 'sbc_init : OPA-SAS coupled via OASIS, but key_oasis3 disabled' )
         IF( ln_mixcpl     )   CALL ctl_stop( 'sbc_init : OPA-SAS coupled via OASIS, but ln_mixcpl = T in OPA' )
      CASE DEFAULT
         CALL ctl_stop( 'sbc_init : unsupported value for nn_components' )
      END SELECT
      !                             !* coupled options
      IF( ln_cpl ) THEN
         IF( .NOT. lk_oasis )   CALL ctl_stop( 'sbc_init : coupled mode with an atmosphere model (ln_cpl=T)',   &
            &                                  '           required to defined key_oasis3' )
      ENDIF
      IF( ln_mixcpl ) THEN
         IF( .NOT. lk_oasis )   CALL ctl_stop( 'sbc_init : mixed forced-coupled mode (ln_mixcpl=T) ',   &
            &                                  '           required to defined key_oasis3' )
         IF( .NOT.ln_cpl    )   CALL ctl_stop( 'sbc_init : mixed forced-coupled mode (ln_mixcpl=T) requires ln_cpl = T' )
         IF( nn_components /= jp_iam_nemo )    &
            &                   CALL ctl_stop( 'sbc_init : the mixed forced-coupled mode (ln_mixcpl=T) ',   &
            &                                   '          not yet working with sas-opa coupling via oasis' )
      ENDIF
      !                             !* sea-ice
      SELECT CASE( nn_ice )
      CASE( 0 )                        !- no ice in the domain
      CASE( 1 )                        !- Ice-cover climatology ("Ice-if" model)  
      CASE( 2 )                        !- SI3  ice model
      CASE( 3 )                        !- CICE ice model
         IF( .NOT.( ln_blk .OR. ln_cpl ) )   CALL ctl_stop( 'sbc_init : CICE sea-ice model requires ln_blk or ln_cpl = T' )
         IF( lk_agrif                    )   CALL ctl_stop( 'sbc_init : CICE sea-ice model not currently available with AGRIF' ) 
      CASE DEFAULT                     !- not supported
      END SELECT
      !
      !                       !**  allocate and set required variables
      !
      !                             !* allocate sbc arrays
      IF( sbc_oce_alloc() /= 0 )   CALL ctl_stop( 'sbc_init : unable to allocate sbc_oce arrays' )
#if ! defined key_si3 && ! defined key_cice
      IF( sbc_ice_alloc() /= 0 )   CALL ctl_stop( 'sbc_init : unable to allocate sbc_ice arrays' )
#endif
      !
      IF( .NOT.ln_isf ) THEN        !* No ice-shelf in the domain : allocate and set to zero
         IF( sbc_isf_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_init : unable to allocate sbc_isf arrays' )
         fwfisf  (:,:)   = 0._wp   ;   risf_tsc  (:,:,:) = 0._wp
         fwfisf_b(:,:)   = 0._wp   ;   risf_tsc_b(:,:,:) = 0._wp
      END IF
      !
      IF( sbc_ssr_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_init : unable to allocate sbc_ssr arrays' )
      IF( .NOT.ln_ssr ) THEN               !* Initialize qrp and erp if no restoring 
         qrp(:,:) = 0._wp
         erp(:,:) = 0._wp
      ENDIF
      !

      IF( nn_ice == 0 ) THEN        !* No sea-ice in the domain : ice fraction is always zero
         IF( nn_components /= jp_iam_opa )   fr_i(:,:) = 0._wp    ! except for OPA in SAS-OPA coupled case
      ENDIF
      !
      sfx   (:,:) = 0._wp           !* salt flux due to freezing/melting
      fmmflx(:,:) = 0._wp           !* freezing minus melting flux

      taum(:,:) = 0._wp             !* wind stress module (needed in GLS in case of reduced restart)

      !                          ! Choice of the Surface Boudary Condition (set nsbc)
      IF( ln_dm2dc ) THEN           !* daily mean to diurnal cycle
         nday_qsr = -1   ! allow initialization at the 1st call
         IF( .NOT.( ln_flx .OR. ln_blk ) .AND. nn_components /= jp_iam_opa )   &
            &   CALL ctl_stop( 'qsr diurnal cycle from daily values requires a flux or bulk formulation' )
      ENDIF
      !                             !* Choice of the Surface Boudary Condition
      !                             (set nsbc)
      !
      ll_purecpl  = ln_cpl .AND. .NOT.ln_mixcpl
      ll_opa      = nn_components == jp_iam_opa
      ll_not_nemo = nn_components /= jp_iam_nemo
      icpt = 0
      !
      IF( ln_usr          ) THEN   ;   nsbc = jp_usr     ; icpt = icpt + 1   ;   ENDIF       ! user defined         formulation
      IF( ln_flx          ) THEN   ;   nsbc = jp_flx     ; icpt = icpt + 1   ;   ENDIF       ! flux                 formulation
      IF( ln_blk          ) THEN   ;   nsbc = jp_blk     ; icpt = icpt + 1   ;   ENDIF       ! bulk                 formulation
      IF( ll_purecpl      ) THEN   ;   nsbc = jp_purecpl ; icpt = icpt + 1   ;   ENDIF       ! Pure Coupled         formulation
      IF( ll_opa          ) THEN   ;   nsbc = jp_none    ; icpt = icpt + 1   ;   ENDIF       ! opa coupling via SAS module
      !
      IF( icpt /= 1 )    CALL ctl_stop( 'sbc_init : choose ONE and only ONE sbc option' )
      !
      IF(lwp) THEN                     !- print the choice of surface flux formulation
         WRITE(numout,*)
         SELECT CASE( nsbc )
         CASE( jp_usr     )   ;   WRITE(numout,*) '   ==>>>   user defined forcing formulation'
         CASE( jp_flx     )   ;   WRITE(numout,*) '   ==>>>   flux formulation'
         CASE( jp_blk     )   ;   WRITE(numout,*) '   ==>>>   bulk formulation'
         CASE( jp_purecpl )   ;   WRITE(numout,*) '   ==>>>   pure coupled formulation'
!!gm abusive use of jp_none ??   ===>>> need to be check and changed by adding a jp_sas parameter
         CASE( jp_none    )   ;   WRITE(numout,*) '   ==>>>   OPA coupled to SAS via oasis'
            IF( ln_mixcpl )       WRITE(numout,*) '               + forced-coupled mixed formulation'
         END SELECT
         IF( ll_not_nemo  )       WRITE(numout,*) '               + OASIS coupled SAS'
      ENDIF
      !
      !                             !* OASIS initialization
      !
      IF( lk_oasis )   CALL sbc_cpl_init( nn_ice )   ! Must be done before: (1) first time step
      !                                              !                      (2) the use of nn_fsbc
      !     nn_fsbc initialization if OPA-SAS coupling via OASIS
      !     SAS time-step has to be declared in OASIS (mandatory) -> nn_fsbc has to be modified accordingly
      IF( nn_components /= jp_iam_nemo ) THEN
         IF( nn_components == jp_iam_opa )   nn_fsbc = cpl_freq('O_SFLX') / NINT(rdt)
         IF( nn_components == jp_iam_sas )   nn_fsbc = cpl_freq('I_SFLX') / NINT(rdt)
         !
         IF(lwp)THEN
            WRITE(numout,*)
            WRITE(numout,*)"   OPA-SAS coupled via OASIS : nn_fsbc re-defined from OASIS namcouple ", nn_fsbc
            WRITE(numout,*)
         ENDIF
      ENDIF
      !
      !                             !* check consistency between model timeline and nn_fsbc
      IF( ln_rst_list .OR. nn_stock /= -1 ) THEN   ! we will do restart files
         IF( MOD( nitend - nit000 + 1, nn_fsbc) /= 0 ) THEN
            WRITE(ctmp1,*) 'sbc_init : experiment length (', nitend - nit000 + 1, ') is NOT a multiple of nn_fsbc (', nn_fsbc, ')'
            CALL ctl_stop( ctmp1, 'Impossible to properly do model restart' )
         ENDIF
         IF( .NOT. ln_rst_list .AND. MOD( nn_stock, nn_fsbc) /= 0 ) THEN   ! we don't use nn_stock if ln_rst_list
            WRITE(ctmp1,*) 'sbc_init : nn_stock (', nn_stock, ') is NOT a multiple of nn_fsbc (', nn_fsbc, ')'
            CALL ctl_stop( ctmp1, 'Impossible to properly do model restart' )
         ENDIF
      ENDIF
      !
      IF( MOD( rday, REAL(nn_fsbc, wp) * rdt ) /= 0 )   &
         &  CALL ctl_warn( 'sbc_init : nn_fsbc is NOT a multiple of the number of time steps in a day' )
      !
      IF( ln_dm2dc .AND. NINT(rday) / ( nn_fsbc * NINT(rdt) ) < 8  )   &
         &   CALL ctl_warn( 'sbc_init : diurnal cycle for qsr: the sampling of the diurnal cycle is too small...' )
      !
   
      !                       !**  associated modules : initialization
      !
                          CALL sbc_ssm_init            ! Sea-surface mean fields initialization
      !
      IF( ln_blk      )   CALL sbc_blk_init            ! bulk formulae initialization

      IF( ln_ssr      )   CALL sbc_ssr_init            ! Sea-Surface Restoring initialization
      !
      IF( ln_isf      )   CALL sbc_isf_init            ! Compute iceshelves
      !
                          CALL sbc_rnf_init            ! Runof initialization
      !
      IF( ln_apr_dyn )    CALL sbc_apr_init            ! Atmo Pressure Forcing initialization
      !
#if defined key_si3
      IF( lk_agrif .AND. nn_ice == 0 ) THEN            ! allocate ice arrays in case agrif + ice-model + no-ice in child grid
                          IF( sbc_ice_alloc() /= 0 )   CALL ctl_stop('STOP', 'sbc_ice_alloc : unable to allocate arrays' )
      ELSEIF( nn_ice == 2 ) THEN
                          CALL ice_init                ! ICE initialization
      ENDIF
#endif
      IF( nn_ice == 3 )   CALL cice_sbc_init( nsbc )   ! CICE initialization
      !
      IF( ln_wave     )   CALL sbc_wave_init           ! surface wave initialisation
      !
      IF( lwxios ) THEN
         CALL iom_set_rstw_var_active('utau_b')
         CALL iom_set_rstw_var_active('vtau_b')
         CALL iom_set_rstw_var_active('qns_b')
         ! The 3D heat content due to qsr forcing is treated in traqsr
         ! CALL iom_set_rstw_var_active('qsr_b')
         CALL iom_set_rstw_var_active('emp_b')
         CALL iom_set_rstw_var_active('sfx_b')
      ENDIF

   END SUBROUTINE sbc_init


   SUBROUTINE sbc( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc  ***
      !!
      !! ** Purpose :   provide at each time-step the ocean surface boundary
      !!                condition (momentum, heat and freshwater fluxes)
      !!
      !! ** Method  :   blah blah  to be written ?????????
      !!                CAUTION : never mask the surface stress field (tke sbc)
      !!
      !! ** Action  : - set the ocean surface boundary condition at before and now
      !!                time step, i.e.
      !!                utau_b, vtau_b, qns_b, qsr_b, emp_n, sfx_b, qrp_b, erp_b
      !!                utau  , vtau  , qns  , qsr  , emp  , sfx  , qrp  , erp
      !!              - updte the ice fraction : fr_i
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      LOGICAL ::   ll_sas, ll_opa   ! local logical
      !
      REAL(wp) ::     zthscl        ! wd  tanh scale
      REAL(wp), DIMENSION(jpi,jpj) ::  zwdht, zwght  ! wd dep over wd limit, wgt  

      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('sbc')
      !
      !                                            ! ---------------------------------------- !
      IF( kt /= nit000 ) THEN                      !          Swap of forcing fields          !
         !                                         ! ---------------------------------------- !
         utau_b(:,:) = utau(:,:)                         ! Swap the ocean forcing fields
         vtau_b(:,:) = vtau(:,:)                         ! (except at nit000 where before fields
         qns_b (:,:) = qns (:,:)                         !  are set at the end of the routine)
         emp_b (:,:) = emp (:,:)
         sfx_b (:,:) = sfx (:,:)
         IF ( ln_rnf ) THEN
            rnf_b    (:,:  ) = rnf    (:,:  )
            rnf_tsc_b(:,:,:) = rnf_tsc(:,:,:)
         ENDIF
         IF( ln_isf )  THEN
            fwfisf_b  (:,:  ) = fwfisf  (:,:  )               
            risf_tsc_b(:,:,:) = risf_tsc(:,:,:)              
         ENDIF
        !
      ENDIF
      !                                            ! ---------------------------------------- !
      !                                            !        forcing field computation         !
      !                                            ! ---------------------------------------- !
      !
      ll_sas = nn_components == jp_iam_sas               ! component flags
      ll_opa = nn_components == jp_iam_opa
      !
      IF( .NOT.ll_sas )   CALL sbc_ssm ( kt )            ! mean ocean sea surface variables (sst_m, sss_m, ssu_m, ssv_m)
      IF( ln_wave     )   CALL sbc_wave( kt )            ! surface waves

      !
      !                                            !==  sbc formulation  ==!
      !                                                   
      SELECT CASE( nsbc )                                ! Compute ocean surface boundary condition
      !                                                  ! (i.e. utau,vtau, qns, qsr, emp, sfx)
      CASE( jp_usr   )     ;   CALL usrdef_sbc_oce( kt )                    ! user defined formulation 
      CASE( jp_flx     )   ;   CALL sbc_flx       ( kt )                    ! flux formulation
      CASE( jp_blk     )
         IF( ll_sas    )       CALL sbc_cpl_rcv   ( kt, nn_fsbc, nn_ice )   ! OPA-SAS coupling: SAS receiving fields from OPA
                               CALL sbc_blk       ( kt )                    ! bulk formulation for the ocean
                               !
      CASE( jp_purecpl )   ;   CALL sbc_cpl_rcv   ( kt, nn_fsbc, nn_ice )   ! pure coupled formulation
      CASE( jp_none    )
         IF( ll_opa    )       CALL sbc_cpl_rcv   ( kt, nn_fsbc, nn_ice )   ! OPA-SAS coupling: OPA receiving fields from SAS
      END SELECT
      !
      IF( ln_mixcpl )          CALL sbc_cpl_rcv   ( kt, nn_fsbc, nn_ice )   ! forced-coupled mixed formulation after forcing
      !
      IF ( ln_wave .AND. (ln_tauwoc .OR. ln_tauw) ) CALL sbc_wstress( )      ! Wind stress provided by waves 
      !
      !                                            !==  Misc. Options  ==!
      !
      SELECT CASE( nn_ice )                                       ! Update heat and freshwater fluxes over sea-ice areas
      CASE(  1 )   ;         CALL sbc_ice_if   ( kt )             ! Ice-cover climatology ("Ice-if" model)
#if defined key_si3
      CASE(  2 )   ;         CALL ice_stp  ( kt, nsbc )           ! SI3 ice model
#endif
      CASE(  3 )   ;         CALL sbc_ice_cice ( kt, nsbc )       ! CICE ice model
      END SELECT

      IF( ln_icebergs    )   THEN
                                     CALL icb_stp( kt )           ! compute icebergs
         ! icebergs may advect into haloes during the icb step and alter emp.
         ! A lbc_lnk is necessary here to ensure restartability (#2113)
         IF( .NOT. ln_passive_mode ) CALL lbc_lnk( 'sbcmod', emp, 'T', 1. ) ! ensure restartability with icebergs
      ENDIF

      IF( ln_isf         )   CALL sbc_isf( kt )                   ! compute iceshelves

      IF( ln_rnf         )   CALL sbc_rnf( kt )                   ! add runoffs to fresh water fluxes

      IF( ln_ssr         )   CALL sbc_ssr( kt )                   ! add SST/SSS damping term

      IF( nn_fwb    /= 0 )   CALL sbc_fwb( kt, nn_fwb, nn_fsbc )  ! control the freshwater budget

      ! Special treatment of freshwater fluxes over closed seas in the model domain
      ! Should not be run if ln_diurnal_only
      IF( l_sbc_clo .AND. (.NOT. ln_diurnal_only) )   CALL sbc_clo( kt )   

!!$!RBbug do not understand why see ticket 667
!!$!clem: it looks like it is necessary for the north fold (in certain circumstances). Don't know why.
!!$      CALL lbc_lnk( 'sbcmod', emp, 'T', 1. )
      IF ( ll_wd ) THEN     ! If near WAD point limit the flux for now
         zthscl = atanh(rn_wd_sbcfra)                     ! taper frac default is .999 
         zwdht(:,:) = sshn(:,:) + ht_0(:,:) - rn_wdmin1   ! do this calc of water
                                                     ! depth above wd limit once
         WHERE( zwdht(:,:) <= 0.0 )
            taum(:,:) = 0.0
            utau(:,:) = 0.0
            vtau(:,:) = 0.0
            qns (:,:) = 0.0
            qsr (:,:) = 0.0
            emp (:,:) = min(emp(:,:),0.0) !can allow puddles to grow but not shrink
            sfx (:,:) = 0.0
         END WHERE
         zwght(:,:) = tanh(zthscl*zwdht(:,:))
         WHERE( zwdht(:,:) > 0.0  .and. zwdht(:,:) < rn_wd_sbcdep ) !  5 m hard limit here is arbitrary
            qsr  (:,:) =  qsr(:,:)  * zwght(:,:)
            qns  (:,:) =  qns(:,:)  * zwght(:,:)
            taum (:,:) =  taum(:,:) * zwght(:,:)
            utau (:,:) =  utau(:,:) * zwght(:,:)
            vtau (:,:) =  vtau(:,:) * zwght(:,:)
            sfx  (:,:) =  sfx(:,:)  * zwght(:,:)
            emp  (:,:) =  emp(:,:)  * zwght(:,:)
         END WHERE
      ENDIF
      !
      IF( kt == nit000 ) THEN                          !   set the forcing field at nit000 - 1    !
         !                                             ! ---------------------------------------- !
         IF( ln_rstart .AND.    &                               !* Restart: read in restart file
            & iom_varid( numror, 'utau_b', ldstop = .FALSE. ) > 0 ) THEN
            IF(lwp) WRITE(numout,*) '          nit000-1 surface forcing fields red in the restart file'
            CALL iom_get( numror, jpdom_autoglo, 'utau_b', utau_b, ldxios = lrxios )   ! before i-stress  (U-point)
            CALL iom_get( numror, jpdom_autoglo, 'vtau_b', vtau_b, ldxios = lrxios )   ! before j-stress  (V-point)
            CALL iom_get( numror, jpdom_autoglo, 'qns_b' , qns_b, ldxios = lrxios  )   ! before non solar heat flux (T-point)
            ! The 3D heat content due to qsr forcing is treated in traqsr
            ! CALL iom_get( numror, jpdom_autoglo, 'qsr_b' , qsr_b, ldxios = lrxios  ) ! before     solar heat flux (T-point)
            CALL iom_get( numror, jpdom_autoglo, 'emp_b', emp_b, ldxios = lrxios  )    ! before     freshwater flux (T-point)
            ! To ensure restart capability with 3.3x/3.4 restart files    !! to be removed in v3.6
            IF( iom_varid( numror, 'sfx_b', ldstop = .FALSE. ) > 0 ) THEN
               CALL iom_get( numror, jpdom_autoglo, 'sfx_b', sfx_b, ldxios = lrxios )  ! before salt flux (T-point)
            ELSE
               sfx_b (:,:) = sfx(:,:)
            ENDIF
         ELSE                                                   !* no restart: set from nit000 values
            IF(lwp) WRITE(numout,*) '          nit000-1 surface forcing fields set to nit000'
            utau_b(:,:) = utau(:,:)
            vtau_b(:,:) = vtau(:,:)
            qns_b (:,:) = qns (:,:)
            emp_b (:,:) = emp (:,:)
            sfx_b (:,:) = sfx (:,:)
         ENDIF
      ENDIF
      !                                                ! ---------------------------------------- !
      IF( lrst_oce ) THEN                              !      Write in the ocean restart file     !
         !                                             ! ---------------------------------------- !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbc : ocean surface forcing fields written in ocean restart file ',   &
            &                    'at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         IF( lwxios ) CALL iom_swap(      cwxios_context          )
         CALL iom_rstput( kt, nitrst, numrow, 'utau_b' , utau, ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'vtau_b' , vtau, ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'qns_b'  , qns, ldxios = lwxios  )
         ! The 3D heat content due to qsr forcing is treated in traqsr
         ! CALL iom_rstput( kt, nitrst, numrow, 'qsr_b'  , qsr  )
         CALL iom_rstput( kt, nitrst, numrow, 'emp_b'  , emp, ldxios = lwxios  )
         CALL iom_rstput( kt, nitrst, numrow, 'sfx_b'  , sfx, ldxios = lwxios  )
         IF( lwxios ) CALL iom_swap(      cxios_context          )
      ENDIF
      !                                                ! ---------------------------------------- !
      !                                                !        Outputs and control print         !
      !                                                ! ---------------------------------------- !
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN
         CALL iom_put( "empmr"  , emp    - rnf )                ! upward water flux
         CALL iom_put( "empbmr" , emp_b  - rnf )                ! before upward water flux ( needed to recalculate the time evolution of ssh in offline )
         CALL iom_put( "saltflx", sfx  )                        ! downward salt flux (includes virtual salt flux beneath ice in linear free surface case)
         CALL iom_put( "fmmflx", fmmflx  )                      ! Freezing-melting water flux
         CALL iom_put( "qt"    , qns  + qsr )                   ! total heat flux
         CALL iom_put( "qns"   , qns        )                   ! solar heat flux
         CALL iom_put( "qsr"   ,       qsr  )                   ! solar heat flux
         IF( nn_ice > 0 .OR. ll_opa )   CALL iom_put( "ice_cover", fr_i )   ! ice fraction
         CALL iom_put( "taum"  , taum       )                   ! wind stress module
         CALL iom_put( "wspd"  , wndm       )                   ! wind speed  module over free ocean or leads in presence of sea-ice
         CALL iom_put( "qrp", qrp )                             ! heat flux damping
         CALL iom_put( "erp", erp )                             ! freshwater flux damping
      ENDIF
      !
      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab2d_1=fr_i              , clinfo1=' fr_i    - : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=(emp-rnf + fwfisf), clinfo1=' emp-rnf - : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=(sfx-rnf + fwfisf), clinfo1=' sfx-rnf - : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=qns              , clinfo1=' qns      - : ', mask1=tmask )
         CALL prt_ctl(tab2d_1=qsr              , clinfo1=' qsr      - : ', mask1=tmask )
         CALL prt_ctl(tab3d_1=tmask            , clinfo1=' tmask    - : ', mask1=tmask, kdim=jpk )
         CALL prt_ctl(tab3d_1=tsn(:,:,:,jp_tem), clinfo1=' sst      - : ', mask1=tmask, kdim=1   )
         CALL prt_ctl(tab3d_1=tsn(:,:,:,jp_sal), clinfo1=' sss      - : ', mask1=tmask, kdim=1   )
         CALL prt_ctl(tab2d_1=utau             , clinfo1=' utau     - : ', mask1=umask,                      &
            &         tab2d_2=vtau             , clinfo2=' vtau     - : ', mask2=vmask )
      ENDIF

      IF( kt == nitend )   CALL sbc_final         ! Close down surface module if necessary
      !
      IF( ln_timing )   CALL timing_stop('sbc')
      !
   END SUBROUTINE sbc


   SUBROUTINE sbc_final
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_final  ***
      !!
      !! ** Purpose :   Finalize CICE (if used)
      !!---------------------------------------------------------------------
      !
      IF( nn_ice == 3 )   CALL cice_sbc_final
      !
   END SUBROUTINE sbc_final

   !!======================================================================
END MODULE sbcmod
