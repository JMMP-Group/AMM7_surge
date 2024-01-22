MODULE step
   !!======================================================================
   !!                       ***  MODULE step  ***
   !! Time-stepping   : manager of the ocean, tracer and ice time stepping
   !!======================================================================
   !! History :  OPA  !  1991-03  (G. Madec)  Original code
   !!             -   !  1991-11  (G. Madec)
   !!             -   !  1992-06  (M. Imbard)  add a first output record
   !!             -   !  1996-04  (G. Madec)  introduction of dynspg
   !!             -   !  1996-04  (M.A. Foujols)  introduction of passive tracer
   !!            8.0  !  1997-06  (G. Madec)  new architecture of call
   !!            8.2  !  1997-06  (G. Madec, M. Imbard, G. Roullet)  free surface
   !!             -   !  1999-02  (G. Madec, N. Grima)  hpg implicit
   !!             -   !  2000-07  (J-M Molines, M. Imbard)  Open Bondary Conditions
   !!   NEMO     1.0  !  2002-06  (G. Madec)  free form, suppress macro-tasking
   !!             -   !  2004-08  (C. Talandier) New trends organization
   !!             -   !  2005-01  (C. Ethe) Add the KPP closure scheme
   !!             -   !  2005-11  (G. Madec)  Reorganisation of tra and dyn calls
   !!             -   !  2006-01  (L. Debreu, C. Mazauric)  Agrif implementation
   !!             -   !  2006-07  (S. Masson)  restart using iom
   !!            3.2  !  2009-02  (G. Madec, R. Benshila)  reintroduicing z*-coordinate
   !!             -   !  2009-06  (S. Masson, G. Madec)  TKE restart compatible with key_cpl
   !!            3.3  !  2010-05  (K. Mogensen, A. Weaver, M. Martin, D. Lea) Assimilation interface
   !!             -   !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase + merge TRC-TRA
   !!            3.4  !  2011-04  (G. Madec, C. Ethe) Merge of dtatem and dtasal
   !!            3.6  !  2012-07  (J. Simeon, G. Madec. C. Ethe)  Online coarsening of outputs
   !!            3.6  !  2014-04  (F. Roquet, G. Madec) New equations of state
   !!            3.6  !  2014-10  (E. Clementi, P. Oddo) Add Qiao vertical mixing in case of waves
   !!            3.7  !  2014-10  (G. Madec)  LDF simplication 
   !!             -   !  2014-12  (G. Madec) remove KPP scheme
   !!             -   !  2015-11  (J. Chanut) free surface simplification (remove filtered free surface)
   !!            4.0  !  2017-05  (G. Madec)  introduction of the vertical physics manager (zdfphy)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   stp             : OPA system time-stepping
   !!----------------------------------------------------------------------
   USE step_oce         ! time stepping definition modules
   !
   USE iom              ! xIOs server

   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp   ! called by nemogcm.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_agrif
   RECURSIVE SUBROUTINE stp( )
      INTEGER             ::   kstp   ! ocean time-step index
#else
   SUBROUTINE stp( kstp )
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index
#endif
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp  ***
      !!
      !! ** Purpose : - Time stepping of OPA  (momentum and active tracer eqs.)
      !!              - Time stepping of SI3 (dynamic and thermodynamic eqs.)
      !!              - Time stepping of TRC  (passive tracer eqs.)
      !!
      !! ** Method  : -1- Update forcings and data
      !!              -2- Update ocean physics
      !!              -3- Compute the t and s trends
      !!              -4- Update t and s
      !!              -5- Compute the momentum trends
      !!              -6- Update the horizontal velocity
      !!              -7- Compute the diagnostics variables (rd,N2, hdiv,w)
      !!              -8- Outputs and diagnostics
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk   ! dummy loop indice
!!gm kcall can be removed, I guess
      INTEGER ::   kcall        ! optional integer argument (dom_vvl_sf_nxt)
      !! ---------------------------------------------------------------------
#if defined key_agrif
      IF( nstop > 0 ) RETURN   ! avoid to go further if an error was detected during previous time step (child grid)
      kstp = nit000 + Agrif_Nb_Step()
      IF( lk_agrif_debug ) THEN
         IF( Agrif_Root() .and. lwp)   WRITE(*,*) '---'
         IF(lwp)   WRITE(*,*) 'Grid Number', Agrif_Fixed(),' time step ', kstp, 'int tstep', Agrif_NbStepint()
      ENDIF
      IF( kstp == nit000 + 1 )   lk_agrif_fstep = .FALSE.
# if defined key_iomput
      IF( Agrif_Nbstepint() == 0 )   CALL iom_swap( cxios_context )
# endif
#endif
      !
      IF( ln_timing )   CALL timing_start('stp')
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! update I/O and calendar 
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( kstp == nit000 ) THEN                       ! initialize IOM context (must be done after nemo_init for AGRIF+XIOS+OASIS)
                             CALL iom_init(      cxios_context          )  ! for model grid (including passible AGRIF zoom)
         IF( ln_crs      )   CALL iom_init( TRIM(cxios_context)//"_crs" )  ! for coarse grid
      ENDIF
      IF( kstp /= nit000 )   CALL day( kstp )         ! Calendar (day was already called at nit000 in day_init)
                             CALL iom_setkt( kstp - nit000 + 1,      cxios_context          )   ! tell IOM we are at time step kstp
      IF( ln_crs         )   CALL iom_setkt( kstp - nit000 + 1, TRIM(cxios_context)//"_crs" )   ! tell IOM we are at time step kstp

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update external forcing (tides, open boundaries, and surface boundary condition (including sea-ice)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_tide    )   CALL sbc_tide( kstp )                   ! update tide potential
      IF( ln_apr_dyn )   CALL sbc_apr ( kstp )                   ! atmospheric pressure (NB: call before bdy_dta which needs ssh_ib) 
      IF( ln_bdy     )   CALL bdy_dta ( kstp, kt_offset = +1 )   ! update dynamic & tracer data at open boundaries
                         CALL sbc     ( kstp )                   ! Sea Boundary Condition (including sea-ice)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update stochastic parameters and random T/S fluctuations
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF( ln_sto_eos ) CALL sto_par( kstp )          ! Stochastic parameters
      IF( ln_sto_eos ) CALL sto_pts( tsn  )          ! Random T/S fluctuations

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Ocean physics update
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !  THERMODYNAMICS
      IF ( .NOT. ln_2d ) THEN 
                         CALL eos_rab( tsb, rab_b )       ! before local thermal/haline expension ratio at T-points
                         CALL eos_rab( tsn, rab_n )       ! now    local thermal/haline expension ratio at T-points
                         CALL bn2    ( tsb, rab_b, rn2b ) ! before Brunt-Vaisala frequency
                         CALL bn2    ( tsn, rab_n, rn2  ) ! now    Brunt-Vaisala frequency
      ENDIF
      
      !  VERTICAL PHYSICS
                         CALL zdf_phy( kstp )         ! vertical physics update (top/bot drag, avt, avs, avm + MLD)

      !  LATERAL  PHYSICS
      !
      IF( l_ldfslp ) THEN                             ! slope of lateral mixing
                         CALL eos( tsb, rhd, gdept_0(:,:,:) )               ! before in situ density

         IF( ln_zps .AND. .NOT. ln_isfcav)                               &
            &            CALL zps_hde    ( kstp, jpts, tsb, gtsu, gtsv,  &  ! Partial steps: before horizontal gradient
            &                                          rhd, gru , grv    )  ! of t, s, rd at the last ocean level

         IF( ln_zps .AND.       ln_isfcav)                               &
            &            CALL zps_hde_isf( kstp, jpts, tsb, gtsu, gtsv, gtui, gtvi,  &  ! Partial steps for top cell (ISF)
            &                                          rhd, gru , grv , grui, grvi   )  ! of t, s, rd at the first ocean level
         IF( ln_traldf_triad ) THEN 
                         CALL ldf_slp_triad( kstp )                       ! before slope for triad operator
         ELSE     
                         CALL ldf_slp     ( kstp, rhd, rn2b )             ! before slope for standard operator
         ENDIF
      ENDIF
      !                                                                   ! eddy diffusivity coeff.
      IF( l_ldftra_time .OR. l_ldfeiv_time )   CALL ldf_tra( kstp )       !       and/or eiv coeff.
      IF( l_ldfdyn_time                    )   CALL ldf_dyn( kstp )       ! eddy viscosity coeff. 
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !  Ocean dynamics : hdiv, ssh, e3, u, v, w
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                            CALL ssh_nxt       ( kstp )  ! after ssh (includes call to div_hor)
      IF( .NOT.ln_linssh )  CALL dom_vvl_sf_nxt( kstp )  ! after vertical scale factors 
                            CALL wzv           ( kstp )  ! now cross-level velocity 
      IF( ln_zad_Aimp )     CALL wAimp         ( kstp )  ! Adaptive-implicit vertical advection partitioning
      IF( .NOT. ln_2d )     CALL eos    ( tsn, rhd, rhop, gdept_n(:,:,:) )  ! now in situ density for hpg computation
                            
                            
                         ua(:,:,:) = 0._wp            ! set dynamics trends to zero
                         va(:,:,:) = 0._wp

      IF(  lk_asminc .AND. ln_asmiau .AND. ln_dyninc )   &
               &         CALL dyn_asm_inc   ( kstp )  ! apply dynamics assimilation increment
      IF( ln_bdy     )   CALL bdy_dyn3d_dmp ( kstp )  ! bdy damping trends
#if defined key_agrif
      IF(.NOT. Agrif_Root())  & 
               &         CALL Agrif_Sponge_dyn        ! momentum sponge
#endif
                         CALL dyn_adv       ( kstp )  ! advection (vector or flux form)
                         CALL dyn_vor       ( kstp )  ! vorticity term including Coriolis
                         CALL dyn_ldf       ( kstp )  ! lateral mixing
      IF( ln_zdfosm  )   CALL dyn_osm       ( kstp )  ! OSMOSIS non-local velocity fluxes
                         CALL dyn_hpg       ( kstp )  ! horizontal gradient of Hydrostatic pressure
                         CALL dyn_spg       ( kstp )  ! surface pressure gradient

                                                      ! With split-explicit free surface, since now transports have been updated and ssha as well
      IF( ln_dynspg_ts ) THEN                         ! vertical scale factors and vertical velocity need to be updated
                            CALL div_hor    ( kstp )              ! Horizontal divergence  (2nd call in time-split case)
         IF(.NOT.ln_linssh) CALL dom_vvl_sf_nxt( kstp, kcall=2 )  ! after vertical scale factors (update depth average component)
      ENDIF
      IF( .NOT. ln_2d )     CALL dyn_zdf    ( kstp )  ! vertical diffusion
      IF( ln_dynspg_ts ) THEN                         ! vertical scale factors and vertical velocity need to be updated
                            CALL wzv        ( kstp )              ! now cross-level velocity 
         IF( ln_zad_Aimp )  CALL wAimp      ( kstp )  ! Adaptive-implicit vertical advection partitioning
      ENDIF


      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! cool skin
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
      IF ( ln_diurnal )  CALL stp_diurnal( kstp )
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! diagnostics and outputs
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_floats  )   CALL flo_stp ( kstp )        ! drifting Floats
      IF( ln_diacfl  )   CALL dia_cfl ( kstp )        ! Courant number diagnostics
                         CALL dia_hth ( kstp )        ! Thermocline depth (20 degres isotherm depth)
      IF( ln_diadct  )   CALL dia_dct ( kstp )        ! Transports
                         CALL dia_ar5 ( kstp )        ! ar5 diag
      IF( ln_diaptr  )   CALL dia_ptr                 ! Poleward adv/ldf TRansports diagnostics
      IF( ln_diaharm )   CALL dia_harm( kstp )        ! Tidal harmonic analysis
                         CALL dia_prod( kstp )        ! ocean model: product diagnostics
                         CALL dia_wri ( kstp )        ! ocean model: outputs
      !
      IF( ln_crs     )   CALL crs_fld       ( kstp )  ! ocean model: online field coarsening & output
      
#if defined key_top
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Passive Tracer Model
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL trc_stp       ( kstp )  ! time-stepping
#endif

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Active tracers                              
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF (.NOT. ln_2d) THEN  ! No tracers in 2D mode
                         tsa(:,:,:,:) = 0._wp         ! set tracer trends to zero

      IF(  lk_asminc .AND. ln_asmiau .AND. &
         & ln_trainc )   CALL tra_asm_inc   ( kstp )  ! apply tracer assimilation increment
                         CALL tra_sbc       ( kstp )  ! surface boundary condition
      IF( ln_traqsr  )   CALL tra_qsr       ( kstp )  ! penetrative solar radiation qsr
      IF( ln_trabbc  )   CALL tra_bbc       ( kstp )  ! bottom heat flux
      IF( ln_trabbl  )   CALL tra_bbl       ( kstp )  ! advective (and/or diffusive) bottom boundary layer scheme
      IF( ln_tradmp  )   CALL tra_dmp       ( kstp )  ! internal damping trends
      IF( ln_bdy     )   CALL bdy_tra_dmp   ( kstp )  ! bdy damping trends
#if defined key_agrif
      IF(.NOT. Agrif_Root())  & 
               &         CALL Agrif_Sponge_tra        ! tracers sponge
#endif
                         CALL tra_adv       ( kstp )  ! horizontal & vertical advection
      IF( ln_zdfosm  )   CALL tra_osm       ( kstp )  ! OSMOSIS non-local tracer fluxes
      IF( lrst_oce .AND. ln_zdfosm ) &
           &             CALL osm_rst( kstp, 'WRITE' )! write OSMOSIS outputs + wn (so must do here) to restarts
                         CALL tra_ldf       ( kstp )  ! lateral mixing

                         CALL tra_zdf       ( kstp )  ! vertical mixing and after tracer fields
      IF( ln_zdfnpc  )   CALL tra_npc       ( kstp )  ! update after fields by non-penetrative convection
      ENDIF ! not ln_2d 
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Set boundary conditions and Swap
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!jc1: For agrif, it would be much better to finalize tracers/momentum here (e.g. bdy conditions) and move the swap 
!!    (and time filtering) after Agrif update. Then restart would be done after and would contain updated fields. 
!!    If so: 
!!    (i) no need to call agrif update at initialization time
!!    (ii) no need to update "before" fields 
!!
!!    Apart from creating new tra_swp/dyn_swp routines, this however: 
!!    (i) makes boundary conditions at initialization time computed from updated fields which is not the case between 
!!    two restarts => restartability issue. One can circumvent this, maybe, by assuming "interface separation", 
!!    e.g. a shift of the feedback interface inside child domain. 
!!    (ii) requires that all restart outputs of updated variables by agrif (e.g. passive tracers/tke/barotropic arrays) are done at the same
!!    place.
!! 
!!jc2: dynnxt must be the latest call. e3t_b are indeed updated in that routine
      IF (.NOT. ln_2d)   CALL tra_nxt       ( kstp )  ! finalize (bcs) tracer fields at next time step and swap
                         CALL dyn_nxt       ( kstp )  ! finalize (bcs) velocities at next time step and swap (always called after tra_nxt)
                         CALL ssh_swp       ( kstp )  ! swap of sea surface height
      IF(.NOT.ln_linssh) CALL dom_vvl_sf_swp( kstp )  ! swap of vertical scale factors
      !
      IF( ln_diahsb  )   CALL dia_hsb       ( kstp )  ! - ML - global conservation diagnostics

!!gm : This does not only concern the dynamics ==>>> add a new title
!!gm2: why ouput restart before AGRIF update?
!!
!!jc: That would be better, but see comment above
!!
      IF( lrst_oce   )   CALL rst_write    ( kstp )   ! write output ocean restart file
      IF( ln_sto_eos )   CALL sto_rst_write( kstp )   ! write restart file for stochastic parameters

#if defined key_agrif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! AGRIF recursive integration
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
                         CALL Agrif_Integrate_ChildGrids( stp )  ! allows to finish all the Child Grids before updating
#endif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Control
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL stp_ctl      ( kstp )
                         
#if defined key_agrif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! AGRIF update
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
      IF( Agrif_NbStepint() == 0 .AND. nstop == 0 ) THEN
                         CALL Agrif_update_all( )                  ! Update all components
      ENDIF
#endif

      IF( ln_diaobs  )   CALL dia_obs      ( kstp )      ! obs-minus-model (assimilation) diagnostics (call after dynamics update)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! File manipulation at the end of the first time step
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<                         
      IF( kstp == nit000 ) THEN                          ! 1st time step only
                                        CALL iom_close( numror )   ! close input  ocean restart file
         IF(lwm)                        CALL FLUSH    ( numond )   ! flush output namelist oce
         IF(lwm .AND. numoni /= -1 )    CALL FLUSH    ( numoni )   ! flush output namelist ice (if exist)
      ENDIF

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Coupled mode
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!gm why lk_oasis and not lk_cpl ????
      IF( lk_oasis .AND. nstop == 0 )   CALL sbc_cpl_snd( kstp )     ! coupled mode : field exchanges
      !
#if defined key_iomput
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Finalize contextes if end of simulation or error detected
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<                         
      IF( kstp == nitend .OR. nstop > 0 ) THEN 
                      CALL iom_context_finalize(       cxios_context         ) ! needed for XIOS+AGRIF
         IF( lrxios ) CALL iom_context_finalize(      crxios_context         )
         IF( ln_crs ) CALL iom_context_finalize( trim(cxios_context)//"_crs" ) ! 
      ENDIF
#endif
      !
      IF( ln_timing )   CALL timing_stop('stp')
      !
   END SUBROUTINE stp
   !
   !!======================================================================
END MODULE step
