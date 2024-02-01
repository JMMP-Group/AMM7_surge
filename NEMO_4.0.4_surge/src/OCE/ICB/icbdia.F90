MODULE icbdia
   !!======================================================================
   !!                       ***  MODULE  icbdia  ***
   !! Icebergs:  initialise variables for iceberg budgets and diagnostics
   !!======================================================================
   !! History : 3.3 !  2010-01  (Martin, Adcroft) Original code
   !!            -  !  2011-03  (Madec)          Part conversion to NEMO form
   !!            -  !                            Removal of mapping from another grid
   !!            -  !  2011-04  (Alderson)       Split into separate modules
   !!            -  !  2011-05  (Alderson)       Budgets are now all here with lots
   !!            -  !                            of silly routines to call to get values in
   !!            -  !                            from the right points in the code
   !!----------------------------------------------------------------------
 
   !!----------------------------------------------------------------------
   !!   icb_dia_init  : initialise iceberg budgeting
   !!   icb_dia       : global iceberg diagnostics
   !!   icb_dia_step  : reset at the beginning of each timestep
   !!   icb_dia_put   : output (via iom_put) iceberg fields
   !!   icb_dia_calve : 
   !!   icb_dia_income: 
   !!   icb_dia_size  : 
   !!   icb_dia_speed : 
   !!   icb_dia_melt  : 
   !!   report_state  : 
   !!   report_consistant : 
   !!   report_budget : 
   !!   report_istate : 
   !!   report_ibudget: 
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE dom_oce        ! ocean domain
   USE in_out_manager ! nemo IO
   USE lib_mpp        ! MPP library
   USE iom            ! I/O library
   USE icb_oce        ! iceberg variables
   USE icbutl         ! iceberg utility routines

   IMPLICIT NONE
   PRIVATE

   PUBLIC   icb_dia_init      ! routine called in icbini.F90 module
   PUBLIC   icb_dia           ! routine called in icbstp.F90 module
   PUBLIC   icb_dia_step      ! routine called in icbstp.F90 module
   PUBLIC   icb_dia_put       ! routine called in icbstp.F90 module
   PUBLIC   icb_dia_melt      ! routine called in icbthm.F90 module
   PUBLIC   icb_dia_size      ! routine called in icbthm.F90 module
   PUBLIC   icb_dia_speed     ! routine called in icbdyn.F90 module
   PUBLIC   icb_dia_calve     ! routine called in icbclv.F90 module
   PUBLIC   icb_dia_income    ! routine called in icbclv.F90 module

   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, PUBLIC  ::   berg_melt       ! Melting+erosion rate of icebergs     [kg/s/m2]
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, PUBLIC  ::   berg_melt_hcflx ! Heat flux to ocean due to heat content of melting icebergs [J/s/m2]
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, PUBLIC  ::   berg_melt_qlat  ! Heat flux to ocean due to latent heat of melting icebergs [J/s/m2]
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, PUBLIC  ::   buoy_melt       ! Buoyancy component of melting rate   [kg/s/m2]
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, PUBLIC  ::   eros_melt       ! Erosion component of melting rate    [kg/s/m2]
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, PUBLIC  ::   conv_melt       ! Convective component of melting rate [kg/s/m2]
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, PUBLIC  ::   bits_src        ! Mass flux from berg erosion into bergy bits [kg/s/m2]
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, PUBLIC  ::   bits_melt       ! Melting rate of bergy bits           [kg/s/m2]
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, PUBLIC  ::   bits_mass       ! Mass distribution of bergy bits      [kg/s/m2]
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, PUBLIC  ::   virtual_area    ! Virtual surface coverage by icebergs [m2]
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, PUBLIC  ::   berg_mass       ! Mass distribution                    [kg/m2]
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC  ::   real_calving    ! Calving rate into iceberg class at
   !                                                                          ! calving locations                    [kg/s]
   
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   tmpc                     ! Temporary work space
   REAL(wp), DIMENSION(:)    , ALLOCATABLE ::   rsumbuf                  ! Temporary work space to reduce mpp exchanges
   INTEGER , DIMENSION(:)    , ALLOCATABLE ::   nsumbuf                  ! Temporary work space to reduce mpp exchanges

   REAL(wp)                      ::  berg_melt_net
   REAL(wp)                      ::  bits_src_net
   REAL(wp)                      ::  bits_melt_net
   REAL(wp)                      ::  bits_mass_start     , bits_mass_end
   REAL(wp)                      ::  floating_heat_start , floating_heat_end
   REAL(wp)                      ::  floating_mass_start , floating_mass_end
   REAL(wp)                      ::  bergs_mass_start    , bergs_mass_end
   REAL(wp)                      ::  stored_start        , stored_heat_start
   REAL(wp)                      ::  stored_end          , stored_heat_end
   REAL(wp)                      ::  calving_src_net     , calving_out_net
   REAL(wp)                      ::  calving_src_heat_net, calving_out_heat_net
   REAL(wp)                      ::  calving_src_heat_used_net
   REAL(wp)                      ::  calving_rcv_net  , calving_ret_net  , calving_used_net
   REAL(wp)                      ::  heat_to_bergs_net, heat_to_ocean_net, melt_net
   REAL(wp)                      ::  calving_to_bergs_net

   INTEGER                       ::  nbergs_start, nbergs_end, nbergs_calved
   INTEGER                       ::  nbergs_melted
   INTEGER                       ::  nspeeding_tickets
   INTEGER , DIMENSION(nclasses) ::  nbergs_calved_by_class

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE icb_dia_init( )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      !
      IF( .NOT.ln_bergdia )   RETURN

      ALLOCATE( berg_melt    (jpi,jpj)   )           ;   berg_melt   (:,:)   = 0._wp
      ALLOCATE( berg_melt_hcflx(jpi,jpj) )           ;   berg_melt_hcflx(:,:)   = 0._wp
      ALLOCATE( berg_melt_qlat(jpi,jpj)  )           ;   berg_melt_qlat(:,:)   = 0._wp
      ALLOCATE( buoy_melt    (jpi,jpj)   )           ;   buoy_melt   (:,:)   = 0._wp
      ALLOCATE( eros_melt    (jpi,jpj)   )           ;   eros_melt   (:,:)   = 0._wp
      ALLOCATE( conv_melt    (jpi,jpj)   )           ;   conv_melt   (:,:)   = 0._wp
      ALLOCATE( bits_src     (jpi,jpj)   )           ;   bits_src    (:,:)   = 0._wp
      ALLOCATE( bits_melt    (jpi,jpj)   )           ;   bits_melt   (:,:)   = 0._wp
      ALLOCATE( bits_mass    (jpi,jpj)   )           ;   bits_mass   (:,:)   = 0._wp
      ALLOCATE( virtual_area (jpi,jpj)   )           ;   virtual_area(:,:)   = 0._wp
      ALLOCATE( berg_mass    (jpi,jpj)   )           ;   berg_mass   (:,:)   = 0._wp
      ALLOCATE( real_calving (jpi,jpj,nclasses) )    ;   real_calving(:,:,:) = 0._wp
      ALLOCATE( tmpc(jpi,jpj) )                      ;   tmpc        (:,:)   = 0._wp

      nbergs_start              = 0
      nbergs_end                = 0
      stored_end                = 0._wp
      nbergs_start              = 0._wp
      stored_start              = 0._wp
      nbergs_melted             = 0
      nbergs_calved             = 0
      nbergs_calved_by_class(:) = 0
      nspeeding_tickets         = 0
      stored_heat_end           = 0._wp
      floating_heat_end         = 0._wp
      floating_mass_end         = 0._wp
      bergs_mass_end            = 0._wp
      bits_mass_end             = 0._wp
      stored_heat_start         = 0._wp
      floating_heat_start       = 0._wp
      floating_mass_start       = 0._wp
      bergs_mass_start          = 0._wp
      bits_mass_start           = 0._wp
      bits_mass_end             = 0._wp
      calving_used_net          = 0._wp
      calving_to_bergs_net      = 0._wp
      heat_to_bergs_net         = 0._wp
      heat_to_ocean_net         = 0._wp
      calving_rcv_net           = 0._wp
      calving_ret_net           = 0._wp
      calving_src_net           = 0._wp
      calving_out_net           = 0._wp
      calving_src_heat_net      = 0._wp
      calving_src_heat_used_net = 0._wp
      calving_out_heat_net      = 0._wp
      melt_net                  = 0._wp
      berg_melt_net             = 0._wp
      bits_melt_net             = 0._wp
      bits_src_net              = 0._wp

      floating_mass_start       = icb_utl_mass( first_berg )
      bergs_mass_start          = icb_utl_mass( first_berg, justbergs=.TRUE. )
      bits_mass_start           = icb_utl_mass( first_berg, justbits =.TRUE. )
      IF( lk_mpp ) THEN
         ALLOCATE( rsumbuf(23) )          ; rsumbuf(:) = 0._wp
         ALLOCATE( nsumbuf(4+nclasses) )  ; nsumbuf(:) = 0
         rsumbuf(1) = floating_mass_start
         rsumbuf(2) = bergs_mass_start
         rsumbuf(3) = bits_mass_start
         CALL mpp_sum( 'icbdia', rsumbuf(1:3), 3 )
         floating_mass_start = rsumbuf(1)
         bergs_mass_start = rsumbuf(2)
         bits_mass_start = rsumbuf(3)
      ENDIF
      !
   END SUBROUTINE icb_dia_init


   SUBROUTINE icb_dia( ld_budge )
      !!----------------------------------------------------------------------
      !! sum all the things we've accumulated so far in the current processor
      !! in MPP case then add these sums across all processors
      !! for this we pack variables into buffer so we only need one mpp_sum
      !!----------------------------------------------------------------------
      LOGICAL, INTENT(in) ::   ld_budge   !
      !
      INTEGER ::   ik
      REAL(wp)::   zunused_calving, ztmpsum, zgrdd_berg_mass, zgrdd_bits_mass
      !!----------------------------------------------------------------------
      !
      IF( .NOT.ln_bergdia )   RETURN

      zunused_calving      = SUM( berg_grid%calving(:,:) )
      ztmpsum              = SUM( berg_grid%floating_melt(:,:) * e1e2t(:,:) * tmask_i(:,:) )
      melt_net             = melt_net + ztmpsum * berg_dt
      calving_out_net      = calving_out_net + ( zunused_calving + ztmpsum ) * berg_dt
      ztmpsum              = SUM( berg_melt(:,:) * e1e2t(:,:) * tmask_i(:,:) )
      berg_melt_net        = berg_melt_net + ztmpsum * berg_dt
      ztmpsum              = SUM( bits_src(:,:) * e1e2t(:,:) * tmask_i(:,:) )
      bits_src_net         = bits_src_net + ztmpsum * berg_dt
      ztmpsum              = SUM( bits_melt(:,:) * e1e2t(:,:) * tmask_i(:,:) )
      bits_melt_net        = bits_melt_net + ztmpsum * berg_dt
      ztmpsum              = SUM( src_calving(:,:) * tmask_i(:,:) )
      calving_ret_net      = calving_ret_net + ztmpsum * berg_dt
      ztmpsum              = SUM( berg_grid%calving_hflx(:,:) * e1e2t(:,:) * tmask_i(:,:) )
      calving_out_heat_net = calving_out_heat_net + ztmpsum * berg_dt   ! Units of J
      !
      IF( ld_budge ) THEN
         stored_end        = SUM( berg_grid%stored_ice(:,:,:) )
         stored_heat_end   = SUM( berg_grid%stored_heat(:,:) )
         floating_mass_end = icb_utl_mass( first_berg )
         bergs_mass_end    = icb_utl_mass( first_berg,justbergs=.TRUE. )
         bits_mass_end     = icb_utl_mass( first_berg,justbits =.TRUE. )
         floating_heat_end = icb_utl_heat( first_berg )
         !
         nbergs_end        = icb_utl_count()
         zgrdd_berg_mass   = SUM( berg_mass(:,:)*e1e2t(:,:)*tmask_i(:,:) )
         zgrdd_bits_mass   = SUM( bits_mass(:,:)*e1e2t(:,:)*tmask_i(:,:) )
         !
         IF( lk_mpp ) THEN
            rsumbuf( 1) = stored_end
            rsumbuf( 2) = stored_heat_end
            rsumbuf( 3) = floating_mass_end
            rsumbuf( 4) = bergs_mass_end
            rsumbuf( 5) = bits_mass_end
            rsumbuf( 6) = floating_heat_end
            rsumbuf( 7) = calving_ret_net
            rsumbuf( 8) = calving_out_net
            rsumbuf( 9) = calving_rcv_net
            rsumbuf(10) = calving_src_net
            rsumbuf(11) = calving_src_heat_net
            rsumbuf(12) = calving_src_heat_used_net
            rsumbuf(13) = calving_out_heat_net
            rsumbuf(14) = calving_used_net
            rsumbuf(15) = calving_to_bergs_net
            rsumbuf(16) = heat_to_bergs_net
            rsumbuf(17) = heat_to_ocean_net
            rsumbuf(18) = melt_net
            rsumbuf(19) = berg_melt_net
            rsumbuf(20) = bits_src_net
            rsumbuf(21) = bits_melt_net
            rsumbuf(22) = zgrdd_berg_mass
            rsumbuf(23) = zgrdd_bits_mass
            !
            CALL mpp_sum( 'icbdia', rsumbuf(1:23), 23)
            !
            stored_end                = rsumbuf( 1)
            stored_heat_end           = rsumbuf( 2)
            floating_mass_end         = rsumbuf( 3)
            bergs_mass_end            = rsumbuf( 4)
            bits_mass_end             = rsumbuf( 5)
            floating_heat_end         = rsumbuf( 6)
            calving_ret_net           = rsumbuf( 7)
            calving_out_net           = rsumbuf( 8)
            calving_rcv_net           = rsumbuf( 9)
            calving_src_net           = rsumbuf(10)
            calving_src_heat_net      = rsumbuf(11)
            calving_src_heat_used_net = rsumbuf(12)
            calving_out_heat_net      = rsumbuf(13)
            calving_used_net          = rsumbuf(14)
            calving_to_bergs_net      = rsumbuf(15)
            heat_to_bergs_net         = rsumbuf(16)
            heat_to_ocean_net         = rsumbuf(17)
            melt_net                  = rsumbuf(18)
            berg_melt_net             = rsumbuf(19)
            bits_src_net              = rsumbuf(20)
            bits_melt_net             = rsumbuf(21)
            zgrdd_berg_mass           = rsumbuf(22)
            zgrdd_bits_mass           = rsumbuf(23)
            !
            nsumbuf(1) = nbergs_end
            nsumbuf(2) = nbergs_calved
            nsumbuf(3) = nbergs_melted
            nsumbuf(4) = nspeeding_tickets
            DO ik = 1, nclasses
               nsumbuf(4+ik) = nbergs_calved_by_class(ik)
            END DO
            CALL mpp_sum( 'icbdia', nsumbuf(1:nclasses+4), nclasses+4 )
            !
            nbergs_end        = nsumbuf(1)
            nbergs_calved     = nsumbuf(2)
            nbergs_melted     = nsumbuf(3)
            nspeeding_tickets = nsumbuf(4)
            DO ik = 1,nclasses
               nbergs_calved_by_class(ik)= nsumbuf(4+ik)
            END DO
            !
         ENDIF
         !
         CALL report_state  ( 'stored ice','kg','',stored_start,'',stored_end,'')
         CALL report_state  ( 'floating','kg','',floating_mass_start,'',floating_mass_end,'',nbergs_end )
         CALL report_state  ( 'icebergs','kg','',bergs_mass_start,'',bergs_mass_end,'')
         CALL report_state  ( 'bits','kg','',bits_mass_start,'',bits_mass_end,'')
         CALL report_istate ( 'berg #','',nbergs_start,'',nbergs_end,'')
         CALL report_ibudget( 'berg #','calved',nbergs_calved, &
            &                          'melted',nbergs_melted, &
            &                          '#',nbergs_start,nbergs_end)
         CALL report_budget( 'stored mass','kg','calving used',calving_used_net, &
            &                              'bergs',calving_to_bergs_net, &
            &                              'stored mass',stored_start,stored_end)
         CALL report_budget( 'floating mass','kg','calving used',calving_to_bergs_net, &
            &                                'bergs',melt_net, &
            &                                'stored mass',floating_mass_start,floating_mass_end)
         CALL report_budget( 'berg mass','kg','calving',calving_to_bergs_net, &
            &                                 'melt+eros',berg_melt_net, &
            &                                 'berg mass',bergs_mass_start,bergs_mass_end)
         CALL report_budget( 'bits mass','kg','eros used',bits_src_net, &
            &                                 'bergs',bits_melt_net, &
            &                                 'stored mass',bits_mass_start,bits_mass_end)
         CALL report_budget( 'net mass','kg','recvd',calving_rcv_net, &
            &                                'rtrnd',calving_ret_net, &
            &                                'net mass',stored_start+floating_mass_start, &
            &                                           stored_end+floating_mass_end)
         CALL report_consistant( 'iceberg mass','kg','gridded',zgrdd_berg_mass,'bergs',bergs_mass_end)
         CALL report_consistant( 'bits mass','kg','gridded',zgrdd_bits_mass,'bits',bits_mass_end)
         CALL report_state( 'net heat','J','',stored_heat_start+floating_heat_start,'', &
            &                                 stored_heat_end+floating_heat_end,'')
         CALL report_state( 'stored heat','J','',stored_heat_start,'',stored_heat_end,'')
         CALL report_state( 'floating heat','J','',floating_heat_start,'',floating_heat_end,'')
         CALL report_budget( 'net heat','J','net heat',calving_src_heat_net, &
            &                               'net heat',calving_out_heat_net, &
            &                               'net heat',stored_heat_start+floating_heat_start, &
            &                                          stored_heat_end+floating_heat_end)
         CALL report_budget( 'stored heat','J','calving used',calving_src_heat_used_net, &
            &                                  'bergs',heat_to_bergs_net, &
            &                                  'net heat',stored_heat_start,stored_heat_end)
         CALL report_budget( 'flting heat','J','calved',heat_to_bergs_net, &
            &                                  'melt',heat_to_ocean_net, &
            &                                  'net heat',floating_heat_start,floating_heat_end)
         IF (nn_verbose_level >= 1) THEN
            CALL report_consistant( 'top interface','kg','from SIS',calving_src_net, &
               &                    'received',calving_rcv_net)
            CALL report_consistant( 'bot interface','kg','sent',calving_out_net, &
               &                    'returned',calving_ret_net)
         ENDIF
         IF (nn_verbose_level > 0) THEN
            WRITE( numicb, '("calved by class = ",i6,20(",",i6))') (nbergs_calved_by_class(ik),ik=1,nclasses)
            IF( nspeeding_tickets > 0 )   WRITE( numicb, '("speeding tickets issued = ",i6)') nspeeding_tickets
         ENDIF
         !
         nbergs_start              = nbergs_end
         stored_start              = stored_end
         nbergs_melted             = 0
         nbergs_calved             = 0
         nbergs_calved_by_class(:) = 0
         nspeeding_tickets         = 0
         stored_heat_start         = stored_heat_end
         floating_heat_start       = floating_heat_end
         floating_mass_start       = floating_mass_end
         bergs_mass_start          = bergs_mass_end
         bits_mass_start           = bits_mass_end
         calving_used_net          = 0._wp
         calving_to_bergs_net      = 0._wp
         heat_to_bergs_net         = 0._wp
         heat_to_ocean_net         = 0._wp
         calving_rcv_net           = 0._wp
         calving_ret_net           = 0._wp
         calving_src_net           = 0._wp
         calving_out_net           = 0._wp
         calving_src_heat_net      = 0._wp
         calving_src_heat_used_net = 0._wp
         calving_out_heat_net      = 0._wp
         melt_net                  = 0._wp
         berg_melt_net             = 0._wp
         bits_melt_net             = 0._wp
         bits_src_net              = 0._wp
      ENDIF
      !
   END SUBROUTINE icb_dia


   SUBROUTINE icb_dia_step
      !!----------------------------------------------------------------------
      !! things to reset at the beginning of each timestep
      !!----------------------------------------------------------------------
      !
      IF( .NOT.ln_bergdia )   RETURN
      berg_melt   (:,:)   = 0._wp
      berg_melt_hcflx(:,:)   = 0._wp
      berg_melt_qlat(:,:)   = 0._wp
      buoy_melt   (:,:)   = 0._wp
      eros_melt   (:,:)   = 0._wp
      conv_melt   (:,:)   = 0._wp
      bits_src    (:,:)   = 0._wp
      bits_melt   (:,:)   = 0._wp
      bits_mass   (:,:)   = 0._wp
      berg_mass   (:,:)   = 0._wp
      virtual_area(:,:)   = 0._wp
      real_calving(:,:,:) = 0._wp
      !
   END SUBROUTINE icb_dia_step


   SUBROUTINE icb_dia_put
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      !
      IF( .NOT.ln_bergdia )   RETURN            !!gm useless iom will control whether it is output or not
      !
      CALL iom_put( "berg_melt"        , berg_melt   (:,:)   )   ! Melt rate of icebergs                     [kg/m2/s]
      !! NB. The berg_melt_hcflx field is currently always zero - see comment in icbthm.F90
      CALL iom_put( "berg_melt_hcflx"  , berg_melt_hcflx(:,:))   ! Heat flux to ocean due to heat content of melting icebergs [J/m2/s]
      CALL iom_put( "berg_melt_qlat"   , berg_melt_qlat(:,:) )   ! Heat flux to ocean due to latent heat of melting icebergs [J/m2/s]
      CALL iom_put( "berg_buoy_melt"   , buoy_melt   (:,:)   )   ! Buoyancy component of iceberg melt rate   [kg/m2/s]
      CALL iom_put( "berg_eros_melt"   , eros_melt   (:,:)   )   ! Erosion component of iceberg melt rate    [kg/m2/s]
      CALL iom_put( "berg_conv_melt"   , conv_melt   (:,:)   )   ! Convective component of iceberg melt rate [kg/m2/s]
      CALL iom_put( "berg_virtual_area", virtual_area(:,:)   )   ! Virtual coverage by icebergs              [m2]
      CALL iom_put( "bits_src"         , bits_src    (:,:)   )   ! Mass source of bergy bits                 [kg/m2/s]
      CALL iom_put( "bits_melt"        , bits_melt   (:,:)   )   ! Melt rate of bergy bits                   [kg/m2/s]
      CALL iom_put( "bits_mass"        , bits_mass   (:,:)   )   ! Bergy bit density field                   [kg/m2]
      CALL iom_put( "berg_mass"        , berg_mass   (:,:)   )   ! Iceberg density field                     [kg/m2]
      CALL iom_put( "berg_real_calving", real_calving(:,:,:) )   ! Calving into iceberg class                [kg/s]
      !
   END SUBROUTINE icb_dia_put


   SUBROUTINE icb_dia_calve( ki, kj, kn, pcalved, pheated )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in)  ::   ki, kj, kn
      REAL(wp), INTENT(in)  ::   pcalved
      REAL(wp), INTENT(in)  ::   pheated
      !!----------------------------------------------------------------------
      !
      IF( .NOT. ln_bergdia ) RETURN
      real_calving(ki,kj,kn)     = real_calving(ki,kj,kn) + pcalved / berg_dt
      nbergs_calved              = nbergs_calved              + 1
      nbergs_calved_by_class(kn) = nbergs_calved_by_class(kn) + 1
      calving_to_bergs_net       = calving_to_bergs_net + pcalved
      heat_to_bergs_net          = heat_to_bergs_net    + pheated
      !
   END SUBROUTINE icb_dia_calve


   SUBROUTINE icb_dia_income( kt,  pcalving_used, pheat_used )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER ,                 INTENT(in)  :: kt
      REAL(wp),                 INTENT(in)  :: pcalving_used
      REAL(wp), DIMENSION(:,:), INTENT(in)  :: pheat_used
      !!----------------------------------------------------------------------
      !
      IF( .NOT.ln_bergdia )   RETURN
      !
      IF( kt == nit000 ) THEN
         stored_start = SUM( berg_grid%stored_ice(:,:,:) )
         CALL mpp_sum( 'icbdia', stored_start )
         !
         stored_heat_start = SUM( berg_grid%stored_heat(:,:) )
         CALL mpp_sum( 'icbdia', stored_heat_start )
         IF (nn_verbose_level > 0) THEN
            WRITE(numicb,'(a,es13.6,a)')   'icb_dia_income: initial stored mass=',stored_start,' kg'
            WRITE(numicb,'(a,es13.6,a)')   'icb_dia_income: initial stored heat=',stored_heat_start,' J'
         ENDIF
      ENDIF
      !
      calving_rcv_net = calving_rcv_net + SUM( berg_grid%calving(:,:) ) * berg_dt
      calving_src_net = calving_rcv_net
      calving_src_heat_net = calving_src_heat_net +  &
         &                      SUM( berg_grid%calving_hflx(:,:) * e1e2t(:,:) ) * berg_dt   ! Units of J
      calving_used_net = calving_used_net + pcalving_used * berg_dt
      calving_src_heat_used_net = calving_src_heat_used_net + SUM( pheat_used(:,:) )
      !
   END SUBROUTINE icb_dia_income


   SUBROUTINE icb_dia_size(ki, kj, pWn, pLn, pAbits,   &
      &                    pmass_scale, pMnew, pnMbits, pz1_e1e2)
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in) ::   ki, kj
      REAL(wp), INTENT(in) ::   pWn, pLn, pAbits, pmass_scale, pMnew, pnMbits, pz1_e1e2
      !!----------------------------------------------------------------------
      !
      IF( .NOT.ln_bergdia )   RETURN
      virtual_area(ki,kj) = virtual_area(ki,kj) + ( pWn * pLn + pAbits ) * pmass_scale      ! m^2
      berg_mass(ki,kj)    = berg_mass(ki,kj) + pMnew * pz1_e1e2                             ! kg/m2
      bits_mass(ki,kj)    = bits_mass(ki,kj) + pnMbits * pz1_e1e2                           ! kg/m2
      !
   END SUBROUTINE icb_dia_size


   SUBROUTINE icb_dia_speed()
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      !
      IF( .NOT.ln_bergdia )   RETURN
      nspeeding_tickets = nspeeding_tickets + 1
      !
   END SUBROUTINE icb_dia_speed


   SUBROUTINE icb_dia_melt(ki, kj, pmnew, pheat_hcflux, pheat_latent, pmass_scale,     &
      &                    pdM, pdMbitsE, pdMbitsM, pdMb, pdMe,   &
      &                    pdMv, pz1_dt_e1e2 )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in) ::   ki, kj
      REAL(wp), INTENT(in) ::   pmnew, pheat_hcflux, pheat_latent, pmass_scale
      REAL(wp), INTENT(in) ::   pdM, pdMbitsE, pdMbitsM, pdMb, pdMe, pdMv, pz1_dt_e1e2
      !!----------------------------------------------------------------------
      !
      IF( .NOT.ln_bergdia )   RETURN
      !
      berg_melt (ki,kj) = berg_melt (ki,kj) + pdM      * pz1_dt_e1e2   ! kg/m2/s
      berg_melt_hcflx (ki,kj) = berg_melt_hcflx (ki,kj) + pheat_hcflux * pz1_dt_e1e2   ! J/m2/s
      berg_melt_qlat (ki,kj) = berg_melt_qlat (ki,kj) + pheat_latent * pz1_dt_e1e2   ! J/m2/s
      bits_src  (ki,kj) = bits_src  (ki,kj) + pdMbitsE * pz1_dt_e1e2   ! mass flux into bergy bitskg/m2/s
      bits_melt (ki,kj) = bits_melt (ki,kj) + pdMbitsM * pz1_dt_e1e2   ! melt rate of bergy bits kg/m2/s
      buoy_melt (ki,kj) = buoy_melt (ki,kj) + pdMb     * pz1_dt_e1e2   ! kg/m2/s
      eros_melt (ki,kj) = eros_melt (ki,kj) + pdMe     * pz1_dt_e1e2   ! erosion rate kg/m2/s
      conv_melt (ki,kj) = conv_melt (ki,kj) + pdMv     * pz1_dt_e1e2   ! kg/m2/s
      heat_to_ocean_net = heat_to_ocean_net + (pheat_hcflux + pheat_latent) * pmass_scale * berg_dt         ! J
      IF( pmnew <= 0._wp ) nbergs_melted = nbergs_melted + 1                        ! Delete the berg if completely melted
      !
   END SUBROUTINE icb_dia_melt


   SUBROUTINE report_state( cd_budgetstr, cd_budgetunits, cd_startstr, pstartval, cd_endstr,   &
      &                     pendval, cd_delstr, kbergs )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER*(*), INTENT(in)           :: cd_budgetstr, cd_budgetunits, cd_startstr, cd_endstr, cd_delstr
      REAL(wp),      INTENT(in)           :: pstartval, pendval
      INTEGER,       INTENT(in), OPTIONAL :: kbergs
      !!----------------------------------------------------------------------
      !
      IF (nn_verbose_level == 0) RETURN
      IF( PRESENT(kbergs) ) THEN
         WRITE(numicb,100) cd_budgetstr // ' state:',                                    &
            &              cd_startstr  // ' start',  pstartval,         cd_budgetunits, &
            &              cd_endstr    // ' end',    pendval,           cd_budgetunits, &
            &              'Delta '     // cd_delstr, pendval-pstartval, cd_budgetunits, &
            &              '# of bergs', kbergs
      ELSE
         WRITE(numicb,100) cd_budgetstr // ' state:',                                   &
            &              cd_startstr  // ' start', pstartval,         cd_budgetunits, &
            &              cd_endstr    // ' end',   pendval,           cd_budgetunits, &
            &              cd_delstr    // 'Delta',  pendval-pstartval, cd_budgetunits
      ENDIF
100   FORMAT(a19,3(a18,"=",es14.7,x,a2,:,","),a12,i8)
      !
   END SUBROUTINE report_state


   SUBROUTINE report_consistant( cd_budgetstr, cd_budgetunits, cd_startstr, pstartval, cd_endstr, pendval)
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER*(*), INTENT(in) :: cd_budgetstr, cd_budgetunits, cd_startstr, cd_endstr
      REAL(wp),      INTENT(in) :: pstartval, pendval
      !!----------------------------------------------------------------------
      !
      IF (nn_verbose_level == 0) RETURN
      WRITE(numicb,200) cd_budgetstr // ' check:',                 &
         &              cd_startstr,    pstartval, cd_budgetunits, &
         &              cd_endstr,      pendval,   cd_budgetunits, &
         &              'error',        (pendval-pstartval)/((pendval+pstartval)+1e-30), 'nd'
200   FORMAT(a19,10(a18,"=",es14.7,x,a2,:,","))
      !
   END SUBROUTINE report_consistant


   SUBROUTINE report_budget( cd_budgetstr, cd_budgetunits, cd_instr, pinval, cd_outstr,   &
      &                      poutval, cd_delstr, pstartval, pendval)
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER*(*), INTENT(in) :: cd_budgetstr, cd_budgetunits, cd_instr, cd_outstr, cd_delstr
      REAL(wp),      INTENT(in) :: pinval, poutval, pstartval, pendval
      !
      REAL(wp) ::   zval
      !!----------------------------------------------------------------------
      !
      IF (nn_verbose_level == 0) RETURN
      zval = ( ( pendval - pstartval ) - ( pinval - poutval ) ) /   &
         &   MAX( 1.e-30, MAX( ABS( pendval - pstartval ) , ABS( pinval - poutval ) ) )
         !
      WRITE(numicb,200) cd_budgetstr // ' budget:', &
         &              cd_instr     // ' in',      pinval,         cd_budgetunits, &
         &              cd_outstr    // ' out',     poutval,        cd_budgetunits, &
         &              'Delta '     // cd_delstr,  pinval-poutval, cd_budgetunits, &
         &              'error',        zval,                       'nd'
  200 FORMAT(a19,3(a18,"=",es14.7,x,a2,:,","),a8,"=",es10.3,x,a2)
      !
   END SUBROUTINE report_budget


   SUBROUTINE report_istate( cd_budgetstr, cd_startstr, pstartval, cd_endstr, pendval, cd_delstr)
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER*(*), INTENT(in) ::   cd_budgetstr, cd_startstr, cd_endstr, cd_delstr
      INTEGER      , INTENT(in) ::   pstartval, pendval
      !!----------------------------------------------------------------------
      !
      IF (nn_verbose_level == 0) RETURN
      WRITE(numicb,100) cd_budgetstr // ' state:',           &
         &              cd_startstr  // ' start', pstartval, &
         &              cd_endstr    // ' end',   pendval,   &
         &              cd_delstr    // 'Delta',  pendval-pstartval
  100 FORMAT(a19,3(a18,"=",i14,x,:,","))
      !
   END SUBROUTINE report_istate


   SUBROUTINE report_ibudget( cd_budgetstr, cd_instr, pinval, cd_outstr, poutval,   &
      &                       cd_delstr, pstartval, pendval)
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER*(*), INTENT(in) :: cd_budgetstr, cd_instr, cd_outstr, cd_delstr
      INTEGER,       INTENT(in) :: pinval, poutval, pstartval, pendval
      !!----------------------------------------------------------------------
      !
      IF (nn_verbose_level == 0) RETURN
      WRITE(numicb,200) cd_budgetstr // ' budget:', &
         &              cd_instr     // ' in',      pinval, &
         &              cd_outstr    // ' out',     poutval, &
         &              'Delta '     // cd_delstr,  pinval-poutval, &
         &              'error',                    ( ( pendval - pstartval ) - ( pinval - poutval ) )
200   FORMAT(a19,10(a18,"=",i14,x,:,","))
      !
   END SUBROUTINE report_ibudget

   !!======================================================================
END MODULE icbdia
