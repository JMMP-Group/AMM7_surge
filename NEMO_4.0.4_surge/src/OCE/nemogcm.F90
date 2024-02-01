MODULE nemogcm
   !!======================================================================
   !!                       ***  MODULE nemogcm   ***
   !! Ocean system   : NEMO GCM (ocean dynamics, on-line tracers, biochemistry and sea-ice)
   !!======================================================================
   !! History :  OPA  ! 1990-10  (C. Levy, G. Madec)  Original code
   !!            7.0  ! 1991-11  (M. Imbard, C. Levy, G. Madec)
   !!            7.1  ! 1993-03  (M. Imbard, C. Levy, G. Madec, O. Marti, M. Guyon, A. Lazar,
   !!                             P. Delecluse, C. Perigaud, G. Caniaux, B. Colot, C. Maes) release 7.1
   !!             -   ! 1992-06  (L.Terray)  coupling implementation
   !!             -   ! 1993-11  (M.A. Filiberti) IGLOO sea-ice
   !!            8.0  ! 1996-03  (M. Imbard, C. Levy, G. Madec, O. Marti, M. Guyon, A. Lazar,
   !!                             P. Delecluse, L.Terray, M.A. Filiberti, J. Vialar, A.M. Treguier, M. Levy) release 8.0
   !!            8.1  ! 1997-06  (M. Imbard, G. Madec)
   !!            8.2  ! 1999-11  (M. Imbard, H. Goosse)  sea-ice model
   !!                 ! 1999-12  (V. Thierry, A-M. Treguier, M. Imbard, M-A. Foujols)  OPEN-MP
   !!                 ! 2000-07  (J-M Molines, M. Imbard)  Open Boundary Conditions  (CLIPPER)
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  F90: Free form and modules
   !!             -   ! 2004-06  (R. Redler, NEC CCRLE, Germany) add OASIS[3/4] coupled interfaces
   !!             -   ! 2004-08  (C. Talandier) New trends organization
   !!             -   ! 2005-06  (C. Ethe) Add the 1D configuration possibility
   !!             -   ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!             -   ! 2006-03  (L. Debreu, C. Mazauric)  Agrif implementation
   !!             -   ! 2006-04  (G. Madec, R. Benshila)  Step reorganization
   !!             -   ! 2007-07  (J. Chanut, A. Sellar) Unstructured open boundaries (BDY)
   !!            3.2  ! 2009-08  (S. Masson)  open/write in the listing file in mpp
   !!            3.3  ! 2010-05  (K. Mogensen, A. Weaver, M. Martin, D. Lea) Assimilation interface
   !!             -   ! 2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.3.1! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!             -   ! 2011-11  (C. Harris) decomposition changes for running with CICE
   !!            3.6  ! 2012-05  (C. Calone, J. Simeon, G. Madec, C. Ethe) Add grid coarsening 
   !!             -   ! 2014-12  (G. Madec) remove KPP scheme and cross-land advection (cla)
   !!            4.0  ! 2016-10  (G. Madec, S. Flavoni)  domain configuration / user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   nemo_gcm      : solve ocean dynamics, tracer, biogeochemistry and/or sea-ice
   !!   nemo_init     : initialization of the NEMO system
   !!   nemo_ctl      : initialisation of the contol print
   !!   nemo_closefile: close remaining open files
   !!   nemo_alloc    : dynamical allocation
   !!----------------------------------------------------------------------
   USE step_oce       ! module used in the ocean time stepping module (step.F90)
   USE phycst         ! physical constant                  (par_cst routine)
   USE domain         ! domain initialization   (dom_init & dom_cfg routines)
   USE closea         ! treatment of closed seas (for ln_closea)
   USE usrdef_nam     ! user defined configuration
   USE tideini        ! tidal components initialization   (tide_ini routine)
   USE bdy_oce,  ONLY : ln_bdy
   USE bdyini         ! open boundary cond. setting       (bdy_init routine)
   USE istate         ! initial state setting          (istate_init routine)
   USE ldfdyn         ! lateral viscosity setting      (ldfdyn_init routine)
   USE ldftra         ! lateral diffusivity setting    (ldftra_init routine)
   USE trdini         ! dyn/tra trends initialization     (trd_init routine)
   USE asminc         ! assimilation increments     
   USE asmbkg         ! writing out state trajectory
   USE diaptr         ! poleward transports           (dia_ptr_init routine)
   USE diadct         ! sections transports           (dia_dct_init routine)
   USE diaobs         ! Observation diagnostics       (dia_obs_init routine)
   USE diacfl         ! CFL diagnostics               (dia_cfl_init routine)
   USE diaharm        ! tidal harmonics diagnostics  (dia_harm_init routine)
   USE step           ! NEMO time-stepping                 (stp     routine)
   USE icbini         ! handle bergs, initialisation
   USE icbstp         ! handle bergs, calving, themodynamics and transport
   USE cpl_oasis3     ! OASIS3 coupling
   USE c1d            ! 1D configuration
   USE step_c1d       ! Time stepping loop for the 1D configuration
   USE dyndmp         ! Momentum damping
   USE stopar         ! Stochastic param.: ???
   USE stopts         ! Stochastic param.: ???
   USE diurnal_bulk   ! diurnal bulk SST 
   USE step_diu       ! diurnal bulk SST timestepping (called from here if run offline)
   USE crsini         ! initialise grid coarsening utility
   USE dia25h         ! 25h mean output
   USE sbc_oce , ONLY : lk_oasis
   USE wet_dry        ! Wetting and drying setting   (wad_init routine)
   USE dom_oce , ONLY : ln_2d 
#if defined key_top
   USE trcini         ! passive tracer initialisation
#endif
#if defined key_nemocice_decomp
   USE ice_domain_size, only: nx_global, ny_global
#endif
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distributed memory computing
   USE mppini         ! shared/distributed memory setting (mpp_init routine)
   USE lbcnfd  , ONLY : isendto, nsndto, nfsloop, nfeloop   ! Setup of north fold exchanges 
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)
#if defined key_iomput
   USE xios           ! xIOserver
#endif
#if defined key_agrif
   USE agrif_all_update   ! Master Agrif update
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   nemo_gcm    ! called by model.F90
   PUBLIC   nemo_init   ! needed by AGRIF
   PUBLIC   nemo_alloc  ! needed by TAM

   CHARACTER(lc) ::   cform_aaa="( /, 'AAAAAAAA', / ) "     ! flag for output listing

#if defined key_mpp_mpi
   ! need MPI_Wtime
   INCLUDE 'mpif.h'
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE nemo_gcm
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_gcm  ***
      !!
      !! ** Purpose :   NEMO solves the primitive equations on an orthogonal
      !!              curvilinear mesh on the sphere.
      !!
      !! ** Method  : - model general initialization
      !!              - launch the time-stepping (stp routine)
      !!              - finalize the run by closing files and communications
      !!
      !! References : Madec, Delecluse, Imbard, and Levy, 1997:  internal report, IPSL.
      !!              Madec, 2008, internal report, IPSL.
      !!----------------------------------------------------------------------
      INTEGER ::   istp   ! time step index
      REAL(wp)::   zstptiming   ! elapsed time for 1 time step
      !!----------------------------------------------------------------------
      !
#if defined key_agrif
      CALL Agrif_Init_Grids()      ! AGRIF: set the meshes
#endif
      !                            !-----------------------!
      CALL nemo_init               !==  Initialisations  ==!
      !                            !-----------------------!
#if defined key_agrif
      CALL Agrif_Declare_Var_dom   ! AGRIF: set the meshes for DOM
      CALL Agrif_Declare_Var       !  "      "   "   "      "  DYN/TRA 
# if defined key_top
      CALL Agrif_Declare_Var_top   !  "      "   "   "      "  TOP
# endif
# if defined key_si3
      CALL Agrif_Declare_Var_ice   !  "      "   "   "      "  Sea ice
# endif
#endif
      ! check that all process are still there... If some process have an error,
      ! they will never enter in step and other processes will wait until the end of the cpu time!
      CALL mpp_max( 'nemogcm', nstop )

      IF(lwp) WRITE(numout,cform_aaa)   ! Flag AAAAAAA

      !                            !-----------------------!
      !                            !==   time stepping   ==!
      !                            !-----------------------!
      istp = nit000
      !
#if defined key_c1d
      DO WHILE ( istp <= nitend .AND. nstop == 0 )    !==  C1D time-stepping  ==!
         CALL stp_c1d( istp )
         istp = istp + 1
      END DO
#else
      !
# if defined key_agrif
      !                                               !==  AGRIF time-stepping  ==!
      CALL Agrif_Regrid()
      !
      ! Recursive update from highest nested level to lowest:
      CALL Agrif_step_child_adj(Agrif_Update_All)
      !
      DO WHILE( istp <= nitend .AND. nstop == 0 )
         CALL stp
         istp = istp + 1
      END DO
      !
# else
      !
      IF( .NOT.ln_diurnal_only ) THEN                 !==  Standard time-stepping  ==!
         !
         DO WHILE( istp <= nitend .AND. nstop == 0 )

            ncom_stp = istp
            IF( ln_timing ) THEN
               zstptiming = MPI_Wtime()
               IF ( istp == ( nit000 + 1 ) ) elapsed_time = zstptiming
               IF ( istp ==         nitend ) elapsed_time = zstptiming - elapsed_time
            ENDIF
            
            CALL stp        ( istp ) 
            istp = istp + 1

            IF( lwp .AND. ln_timing )   WRITE(numtime,*) 'timing step ', istp-1, ' : ', MPI_Wtime() - zstptiming

         END DO
         !
      ELSE                                            !==  diurnal SST time-steeping only  ==!
         !
         DO WHILE( istp <= nitend .AND. nstop == 0 )
            CALL stp_diurnal( istp )   ! time step only the diurnal SST 
            istp = istp + 1
         END DO
         !
      ENDIF
      !
# endif
      !
#endif
      !
      IF( ln_diaobs   )   CALL dia_obs_wri
      !
      IF( ln_icebergs )   CALL icb_end( nitend )

      !                            !------------------------!
      !                            !==  finalize the run  ==!
      !                            !------------------------!
      IF(lwp) WRITE(numout,cform_aaa)        ! Flag AAAAAAA
      !
      IF( nstop /= 0 .AND. lwp ) THEN        ! error print
         WRITE(ctmp1,*) '   ==>>>   nemo_gcm: a total of ', nstop, ' errors have been found'
         IF( ngrdstop > 0 ) THEN
            WRITE(ctmp9,'(i2)') ngrdstop
            WRITE(ctmp2,*) '           E R R O R detected in Agrif grid '//TRIM(ctmp9)
            WRITE(ctmp3,*) '           Look for "E R R O R" messages in all existing '//TRIM(ctmp9)//'_ocean_output* files'
            CALL ctl_stop( ' ', ctmp1, ' ', ctmp2, ' ', ctmp3 )
         ELSE
            WRITE(ctmp2,*) '           Look for "E R R O R" messages in all existing ocean_output* files'
            CALL ctl_stop( ' ', ctmp1, ' ', ctmp2 )
         ENDIF
      ENDIF
      !
      IF( ln_timing )   CALL timing_finalize
      !
      CALL nemo_closefile
      !
#if defined key_iomput
                                    CALL xios_finalize  ! end mpp communications with xios
      IF( lk_oasis     )            CALL cpl_finalize   ! end coupling and mpp communications with OASIS
#else
      IF    ( lk_oasis ) THEN   ;   CALL cpl_finalize   ! end coupling and mpp communications with OASIS
      ELSEIF( lk_mpp   ) THEN   ;   CALL mppstop      ! end mpp communications
      ENDIF
#endif
      !
      IF(lwm) THEN
         IF( nstop == 0 ) THEN   ;   STOP 0
         ELSE                    ;   STOP 123
         ENDIF
      ENDIF
      !
   END SUBROUTINE nemo_gcm


   SUBROUTINE nemo_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_init  ***
      !!
      !! ** Purpose :   initialization of the NEMO GCM
      !!----------------------------------------------------------------------
      INTEGER ::   ios, ilocal_comm   ! local integers
      !!
      NAMELIST/namctl/ ln_ctl   , sn_cfctl, nn_print, nn_ictls, nn_ictle,   &
         &             nn_isplt , nn_jsplt, nn_jctls, nn_jctle,             &
         &             ln_timing, ln_diacfl
      NAMELIST/namcfg/ ln_read_cfg, cn_domcfg, ln_closea, ln_write_cfg, cn_domcfg_out, ln_use_jattr
      !!----------------------------------------------------------------------
      !
      cxios_context = 'nemo'
      !
      !                             !-------------------------------------------------!
      !                             !     set communicator & select the local rank    !
      !                             !  must be done as soon as possible to get narea  !
      !                             !-------------------------------------------------!
      !
#if defined key_iomput
      IF( Agrif_Root() ) THEN
         IF( lk_oasis ) THEN
            CALL cpl_init( "oceanx", ilocal_comm )                               ! nemo local communicator given by oasis
            CALL xios_initialize( "not used"       , local_comm =ilocal_comm )   ! send nemo communicator to xios
         ELSE
            CALL xios_initialize( "for_xios_mpi_id", return_comm=ilocal_comm )   ! nemo local communicator given by xios
         ENDIF
      ENDIF
      CALL mpp_start( ilocal_comm )
#else
      IF( lk_oasis ) THEN
         IF( Agrif_Root() ) THEN
            CALL cpl_init( "oceanx", ilocal_comm )          ! nemo local communicator given by oasis
         ENDIF
         CALL mpp_start( ilocal_comm )
      ELSE
         CALL mpp_start( )
      ENDIF
#endif
      !
      narea = mpprank + 1               ! mpprank: the rank of proc (0 --> mppsize -1 )
      lwm = (narea == 1)                ! control of output namelists
      !
      !                             !---------------------------------------------------------------!
      !                             ! Open output files, reference and configuration namelist files !
      !                             !---------------------------------------------------------------!
      !
      ! open ocean.output as soon as possible to get all output prints (including errors messages)
      IF( lwm )   CALL ctl_opn(     numout,        'ocean.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )
      ! open reference and configuration namelist files
                  CALL ctl_opn( numnam_ref,        'namelist_ref',     'OLD', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )
                  CALL ctl_opn( numnam_cfg,        'namelist_cfg',     'OLD', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )
      IF( lwm )   CALL ctl_opn(     numond, 'output.namelist.dyn', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )
      ! open /dev/null file to be able to supress output write easily
      IF( Agrif_Root() ) THEN
                  CALL ctl_opn(     numnul,           '/dev/null', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )
#ifdef key_agrif
      ELSE
                  numnul = Agrif_Parent(numnul)
#endif
      ENDIF
      !
      !                             !--------------------!
      !                             ! Open listing units !  -> need ln_ctl from namctl to define lwp
      !                             !--------------------!
      !
      REWIND( numnam_ref )              ! Namelist namctl in reference namelist
      READ  ( numnam_ref, namctl, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namctl in reference namelist' )
      REWIND( numnam_cfg )              ! Namelist namctl in confguration namelist
      READ  ( numnam_cfg, namctl, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namctl in configuration namelist' )
      !
      lwp = (narea == 1) .OR. ln_ctl    ! control of all listing output print
      !
      IF(lwp) THEN                      ! open listing units
         !
         IF( .NOT. lwm )   &            ! alreay opened for narea == 1
            &            CALL ctl_opn( numout, 'ocean.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE., narea )
         !
         WRITE(numout,*)
         WRITE(numout,*) '   CNRS - NERC - Met OFFICE - MERCATOR-ocean - CMCC'
         WRITE(numout,*) '                       NEMO team'
         WRITE(numout,*) '            Ocean General Circulation Model'
         WRITE(numout,*) '                NEMO version 4.0  (2019) '
         WRITE(numout,*)
         WRITE(numout,*) "           ._      ._      ._      ._      ._    "
         WRITE(numout,*) "       _.-._)`\_.-._)`\_.-._)`\_.-._)`\_.-._)`\_ "
         WRITE(numout,*)
         WRITE(numout,*) "           o         _,           _,             "
         WRITE(numout,*) "            o      .' (        .-' /             "
         WRITE(numout,*) "           o     _/..._'.    .'   /              "
         WRITE(numout,*) "      (    o .-'`      ` '-./  _.'               "
         WRITE(numout,*) "       )    ( o)           ;= <_         (       "
         WRITE(numout,*) "      (      '-.,\\__ __.-;`\   '.        )      "
         WRITE(numout,*) "       )  )       \) |`\ \)  '.   \      (   (   "
         WRITE(numout,*) "      (  (           \_/       '-._\      )   )  "
         WRITE(numout,*) "       )  ) jgs                     `    (   (   "
         WRITE(numout,*) "     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ "
         WRITE(numout,*)
         !
         WRITE(numout,cform_aaa)                                        ! Flag AAAAAAA
         !
      ENDIF
      !
      ! finalize the definition of namctl variables
      IF( sn_cfctl%l_config ) THEN
         ! Activate finer control of report outputs
         ! optionally switch off output from selected areas (note this only
         ! applies to output which does not involve global communications)
         IF( ( narea < sn_cfctl%procmin .OR. narea > sn_cfctl%procmax  ) .OR. &
           & ( MOD( narea - sn_cfctl%procmin, sn_cfctl%procincr ) /= 0 ) )    &
           &   CALL nemo_set_cfctl( sn_cfctl, .FALSE., .FALSE. )
      ELSE
         ! Use ln_ctl to turn on or off all options.
         CALL nemo_set_cfctl( sn_cfctl, ln_ctl, .TRUE. )
      ENDIF
      !
      IF(lwm) WRITE( numond, namctl )
      !
      !                             !------------------------------------!
      !                             !  Set global domain size parameters !
      !                             !------------------------------------!
      !
      REWIND( numnam_ref )              ! Namelist namcfg in reference namelist
      READ  ( numnam_ref, namcfg, IOSTAT = ios, ERR = 903 )
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcfg in reference namelist' )
      REWIND( numnam_cfg )              ! Namelist namcfg in confguration namelist
      READ  ( numnam_cfg, namcfg, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcfg in configuration namelist' )   
      !
      IF( ln_read_cfg ) THEN            ! Read sizes in domain configuration file
         CALL domain_cfg ( cn_cfg, nn_cfg, jpiglo, jpjglo, jpkglo, jperio )
      ELSE                              ! user-defined namelist
         CALL usr_def_nam( cn_cfg, nn_cfg, jpiglo, jpjglo, jpkglo, jperio )
      ENDIF
      !
      IF(lwm)   WRITE( numond, namcfg )
      !
      !                             !-----------------------------------------!
      !                             ! mpp parameters and domain decomposition !
      !                             !-----------------------------------------!
      CALL mpp_init

      ! Now we know the dimensions of the grid and numout has been set: we can allocate arrays
      CALL nemo_alloc()

      !                             !-------------------------------!
      !                             !  NEMO general initialization  !
      !                             !-------------------------------!

      CALL nemo_ctl                          ! Control prints
      !
      !                                      ! General initialization
      IF( ln_timing    )   CALL timing_init     ! timing
      IF( ln_timing    )   CALL timing_start( 'nemo_init')
      !
                           CALL     phy_cst         ! Physical constants
                           CALL     eos_init        ! Equation of state
      IF( lk_c1d       )   CALL     c1d_init        ! 1D column configuration
                           CALL     wad_init        ! Wetting and drying options
                           CALL     dom_init("OPA") ! Domain
      IF( ln_crs       )   CALL     crs_init        ! coarsened grid: domain initialization 
      IF( ln_ctl       )   CALL prt_ctl_init        ! Print control
      
      CALL diurnal_sst_bulk_init                ! diurnal sst
      IF( ln_diurnal   )   CALL diurnal_sst_coolskin_init   ! cool skin   
      !                            
      IF( ln_diurnal_only ) THEN                   ! diurnal only: a subset of the initialisation routines
         CALL  istate_init                            ! ocean initial state (Dynamics and tracers)
         CALL     sbc_init                            ! Forcings : surface module
         CALL tra_qsr_init                            ! penetrative solar radiation qsr
         IF( ln_diaobs ) THEN                         ! Observation & model comparison
            CALL dia_obs_init                            ! Initialize observational data
            CALL dia_obs( nit000 - 1 )                   ! Observation operator for restart
         ENDIF     
         IF( lk_asminc )   CALL asm_inc_init          ! Assimilation increments
         !
         RETURN                                       ! end of initialization
      ENDIF
      
                           CALL  istate_init    ! ocean initial state (Dynamics and tracers)

      !                                      ! external forcing 
                           CALL    tide_init    ! tidal harmonics
                           CALL     sbc_init    ! surface boundary conditions (including sea-ice)
                           CALL     bdy_init    ! Open boundaries initialisation

      !                                      ! Ocean physics
                           CALL zdf_phy_init    ! Vertical physics
                                     
      !                                         ! Lateral physics
      IF (.NOT. ln_2d)     CALL ldf_tra_init      ! Lateral ocean tracer physics
                           CALL ldf_eiv_init      ! eddy induced velocity param.
                           CALL ldf_dyn_init      ! Lateral ocean momentum physics

      !                                      ! Active tracers
      IF (.NOT. ln_2d) THEN  ! No tracers in 2D mode 
          IF( ln_traqsr    )   CALL tra_qsr_init      ! penetrative solar radiation qsr
                               CALL tra_bbc_init      ! bottom heat flux
                               CALL tra_bbl_init      ! advective (and/or diffusive) bottom boundary layer scheme
                               CALL tra_dmp_init      ! internal tracer damping
                               CALL tra_adv_init      ! horizontal & vertical advection
                               CALL tra_ldf_init      ! lateral mixing
      ENDIF

      !                                      ! Dynamics
      IF( lk_c1d       )   CALL dyn_dmp_init      ! internal momentum damping
                           CALL dyn_adv_init      ! advection (vector or flux form)
                           CALL dyn_vor_init      ! vorticity term including Coriolis
                           CALL dyn_ldf_init      ! lateral mixing
                           CALL dyn_hpg_init      ! horizontal gradient of Hydrostatic pressure
                           CALL dyn_spg_init      ! surface pressure gradient

#if defined key_top
      !                                      ! Passive tracers
                           CALL     trc_init
#endif
      IF( l_ldfslp     )   CALL ldf_slp_init    ! slope of lateral mixing

      !                                      ! Icebergs
                           CALL icb_init( rdt, nit000)   ! initialise icebergs instance

      !                                      ! Misc. options
                           CALL sto_par_init    ! Stochastic parametrization
      IF( ln_sto_eos   )   CALL sto_pts_init    ! RRandom T/S fluctuations
     
      !                                      ! Diagnostics
                           CALL     flo_init    ! drifting Floats
      IF( ln_diacfl    )   CALL dia_cfl_init    ! Initialise CFL diagnostics
                           CALL dia_ptr_init    ! Poleward TRansports initialization
                           CALL dia_dct_init    ! Sections tranports
                           CALL dia_hsb_init    ! heat content, salt content and volume budgets
                           CALL     trd_init    ! Mixed-layer/Vorticity/Integral constraints trends
                           CALL dia_obs_init    ! Initialize observational data
                           CALL dia_25h_init    ! 25h mean  outputs
                           CALL dia_harm_init   ! tidal harmonics outputs
     IF( ln_diaobs    )    CALL dia_obs( nit000-1 )   ! Observation operator for restart

      !                                      ! Assimilation increments
      IF( lk_asminc    )   CALL asm_inc_init    ! Initialize assimilation increments
      !
      IF(lwp) WRITE(numout,cform_aaa)           ! Flag AAAAAAA
      !
      IF( ln_timing    )   CALL timing_stop( 'nemo_init')
      !
   END SUBROUTINE nemo_init


   SUBROUTINE nemo_ctl
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_ctl  ***
      !!
      !! ** Purpose :   control print setting
      !!
      !! ** Method  : - print namctl and namcfg information and check some consistencies
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'nemo_ctl: Control prints'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '   Namelist namctl'
         WRITE(numout,*) '      run control (for debugging)     ln_ctl     = ', ln_ctl
         WRITE(numout,*) '       finer control over o/p sn_cfctl%l_config  = ', sn_cfctl%l_config
         WRITE(numout,*) '                              sn_cfctl%l_runstat = ', sn_cfctl%l_runstat
         WRITE(numout,*) '                              sn_cfctl%l_trcstat = ', sn_cfctl%l_trcstat
         WRITE(numout,*) '                              sn_cfctl%l_oceout  = ', sn_cfctl%l_oceout
         WRITE(numout,*) '                              sn_cfctl%l_layout  = ', sn_cfctl%l_layout
         WRITE(numout,*) '                              sn_cfctl%l_mppout  = ', sn_cfctl%l_mppout
         WRITE(numout,*) '                              sn_cfctl%l_mpptop  = ', sn_cfctl%l_mpptop
         WRITE(numout,*) '                              sn_cfctl%procmin   = ', sn_cfctl%procmin  
         WRITE(numout,*) '                              sn_cfctl%procmax   = ', sn_cfctl%procmax  
         WRITE(numout,*) '                              sn_cfctl%procincr  = ', sn_cfctl%procincr 
         WRITE(numout,*) '                              sn_cfctl%ptimincr  = ', sn_cfctl%ptimincr 
         WRITE(numout,*) '      level of print                  nn_print   = ', nn_print
         WRITE(numout,*) '      Start i indice for SUM control  nn_ictls   = ', nn_ictls
         WRITE(numout,*) '      End i indice for SUM control    nn_ictle   = ', nn_ictle
         WRITE(numout,*) '      Start j indice for SUM control  nn_jctls   = ', nn_jctls
         WRITE(numout,*) '      End j indice for SUM control    nn_jctle   = ', nn_jctle
         WRITE(numout,*) '      number of proc. following i     nn_isplt   = ', nn_isplt
         WRITE(numout,*) '      number of proc. following j     nn_jsplt   = ', nn_jsplt
         WRITE(numout,*) '      timing by routine               ln_timing  = ', ln_timing
         WRITE(numout,*) '      CFL diagnostics                 ln_diacfl  = ', ln_diacfl
      ENDIF
      !
      nprint    = nn_print          ! convert DOCTOR namelist names into OLD names
      nictls    = nn_ictls
      nictle    = nn_ictle
      njctls    = nn_jctls
      njctle    = nn_jctle
      isplt     = nn_isplt
      jsplt     = nn_jsplt

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namcfg'
         WRITE(numout,*) '      read domain configuration file                ln_read_cfg      = ', ln_read_cfg
         WRITE(numout,*) '         filename to be read                           cn_domcfg     = ', TRIM(cn_domcfg)
         WRITE(numout,*) '         keep closed seas in the domain (if exist)     ln_closea     = ', ln_closea
         WRITE(numout,*) '      create a configuration definition file        ln_write_cfg     = ', ln_write_cfg
         WRITE(numout,*) '         filename to be written                        cn_domcfg_out = ', TRIM(cn_domcfg_out)
         WRITE(numout,*) '      use file attribute if exists as i/p j-start   ln_use_jattr     = ', ln_use_jattr
      ENDIF
      IF( .NOT.ln_read_cfg )   ln_closea = .false.   ! dealing possible only with a domcfg file
      !
      !                             ! Parameter control
      !
      IF( ln_ctl ) THEN                 ! sub-domain area indices for the control prints
         IF( lk_mpp .AND. jpnij > 1 ) THEN
            isplt = jpni   ;   jsplt = jpnj   ;   ijsplt = jpni*jpnj   ! the domain is forced to the real split domain
         ELSE
            IF( isplt == 1 .AND. jsplt == 1  ) THEN
               CALL ctl_warn( ' - isplt & jsplt are equal to 1',   &
                  &           ' - the print control will be done over the whole domain' )
            ENDIF
            ijsplt = isplt * jsplt            ! total number of processors ijsplt
         ENDIF
         IF(lwp) WRITE(numout,*)'          - The total number of processors over which the'
         IF(lwp) WRITE(numout,*)'            print control will be done is ijsplt : ', ijsplt
         !
         !                              ! indices used for the SUM control
         IF( nictls+nictle+njctls+njctle == 0 )   THEN    ! print control done over the default area
            lsp_area = .FALSE.
         ELSE                                             ! print control done over a specific  area
            lsp_area = .TRUE.
            IF( nictls < 1 .OR. nictls > jpiglo )   THEN
               CALL ctl_warn( '          - nictls must be 1<=nictls>=jpiglo, it is forced to 1' )
               nictls = 1
            ENDIF
            IF( nictle < 1 .OR. nictle > jpiglo )   THEN
               CALL ctl_warn( '          - nictle must be 1<=nictle>=jpiglo, it is forced to jpiglo' )
               nictle = jpiglo
            ENDIF
            IF( njctls < 1 .OR. njctls > jpjglo )   THEN
               CALL ctl_warn( '          - njctls must be 1<=njctls>=jpjglo, it is forced to 1' )
               njctls = 1
            ENDIF
            IF( njctle < 1 .OR. njctle > jpjglo )   THEN
               CALL ctl_warn( '          - njctle must be 1<=njctle>=jpjglo, it is forced to jpjglo' )
               njctle = jpjglo
            ENDIF
         ENDIF
      ENDIF
      !
      IF( 1._wp /= SIGN(1._wp,-0._wp)  )   CALL ctl_stop( 'nemo_ctl: The intrinsec SIGN function follows f2003 standard.',  &
         &                                                'Compile with key_nosignedzero enabled:',   &
         &                                                '--> add -Dkey_nosignedzero to the definition of %CPP in your arch file' )
      !
#if defined key_agrif
      IF( ln_timing )   CALL ctl_stop( 'AGRIF not implemented with ln_timing = true')
#endif
      !
   END SUBROUTINE nemo_ctl


   SUBROUTINE nemo_closefile
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_closefile  ***
      !!
      !! ** Purpose :   Close the files
      !!----------------------------------------------------------------------
      !
      IF( lk_mpp )   CALL mppsync
      !
      CALL iom_close                                 ! close all input/output files managed by iom_*
      !
      IF( numstp          /= -1 )   CLOSE( numstp          )   ! time-step file
      IF( numrun          /= -1 )   CLOSE( numrun          )   ! run statistics file
      IF( numnam_ref      /= -1 )   CLOSE( numnam_ref      )   ! oce reference namelist
      IF( numnam_cfg      /= -1 )   CLOSE( numnam_cfg      )   ! oce configuration namelist
      IF( lwm.AND.numond  /= -1 )   CLOSE( numond          )   ! oce output namelist
      IF( numnam_ice_ref  /= -1 )   CLOSE( numnam_ice_ref  )   ! ice reference namelist
      IF( numnam_ice_cfg  /= -1 )   CLOSE( numnam_ice_cfg  )   ! ice configuration namelist
      IF( lwm.AND.numoni  /= -1 )   CLOSE( numoni          )   ! ice output namelist
      IF( numevo_ice      /= -1 )   CLOSE( numevo_ice      )   ! ice variables (temp. evolution)
      IF( numout          /=  6 )   CLOSE( numout          )   ! standard model output file
      IF( numdct_vol      /= -1 )   CLOSE( numdct_vol      )   ! volume transports
      IF( numdct_heat     /= -1 )   CLOSE( numdct_heat     )   ! heat transports
      IF( numdct_salt     /= -1 )   CLOSE( numdct_salt     )   ! salt transports
      !
      numout = 6                                     ! redefine numout in case it is used after this point...
      !
   END SUBROUTINE nemo_closefile


   SUBROUTINE nemo_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of the OPA modules
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      USE diawri    , ONLY : dia_wri_alloc
      USE dom_oce   , ONLY : dom_oce_alloc
      USE trc_oce   , ONLY : trc_oce_alloc
      USE bdy_oce   , ONLY : bdy_oce_alloc
      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      ierr =        oce_alloc    ()    ! ocean 
      ierr = ierr + dia_wri_alloc()
      ierr = ierr + dom_oce_alloc()    ! ocean domain
      ierr = ierr + zdf_oce_alloc()    ! ocean vertical physics
      ierr = ierr + trc_oce_alloc()    ! shared TRC / TRA arrays
      ierr = ierr + bdy_oce_alloc()    ! bdy masks (incl. initialization)
      !
      CALL mpp_sum( 'nemogcm', ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'nemo_alloc: unable to allocate standard ocean arrays' )
      !
   END SUBROUTINE nemo_alloc

   
   SUBROUTINE nemo_set_cfctl(sn_cfctl, setto, for_all )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_set_cfctl  ***
      !!
      !! ** Purpose :   Set elements of the output control structure to setto.
      !!                for_all should be .false. unless all areas are to be
      !!                treated identically.
      !!
      !! ** Method  :   Note this routine can be used to switch on/off some
      !!                types of output for selected areas but any output types
      !!                that involve global communications (e.g. mpp_max, glob_sum)
      !!                should be protected from selective switching by the
      !!                for_all argument
      !!----------------------------------------------------------------------
      LOGICAL :: setto, for_all
      TYPE(sn_ctl) :: sn_cfctl
      !!----------------------------------------------------------------------
      IF( for_all ) THEN
         sn_cfctl%l_runstat = setto
         sn_cfctl%l_trcstat = setto
      ENDIF
      sn_cfctl%l_oceout  = setto
      sn_cfctl%l_layout  = setto
      sn_cfctl%l_mppout  = setto
      sn_cfctl%l_mpptop  = setto
   END SUBROUTINE nemo_set_cfctl

   !!======================================================================
END MODULE nemogcm

