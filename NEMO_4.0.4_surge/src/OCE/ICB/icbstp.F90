MODULE icbstp
   !!======================================================================
   !!                       ***  MODULE  icbstp  ***
   !! Icebergs:  initialise variables for iceberg tracking
   !!======================================================================
   !! History : 3.3.1 !  2010-01  (Martin&Adcroft) Original code
   !!            -    !  2011-03  (Madec)          Part conversion to NEMO form
   !!            -    !                            Removal of mapping from another grid
   !!            -    !  2011-04  (Alderson)       Split into separate modules
   !!            -    !                            Move budgets to icbdia routine
   !!            -    !  2011-05  (Alderson)       Add call to copy forcing arrays
   !!            -    !                            into icb copies with haloes
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   icb_stp       : start iceberg tracking
   !!   icb_end       : end   iceberg tracking
   !!----------------------------------------------------------------------
   USE par_oce        ! nemo parameters
   USE dom_oce        ! ocean domain
   USE sbc_oce        ! ocean surface forcing
   USE phycst         ! physical constants
   !
   USE icb_oce        ! iceberg: define arrays
   USE icbini         ! iceberg: initialisation routines
   USE icbutl         ! iceberg: utility routines
   USE icbrst         ! iceberg: restart routines
   USE icbdyn         ! iceberg: dynamics (ie advection) routines
   USE icbclv         ! iceberg: calving routines
   USE icbthm         ! iceberg: thermodynamics routines
   USE icblbc         ! iceberg: lateral boundary routines (including mpp)
   USE icbtrj         ! iceberg: trajectory I/O routines
   USE icbdia         ! iceberg: budget
   !
   USE in_out_manager ! nemo IO
   USE lib_mpp        ! massively parallel library 
   USE iom            ! I/O manager
   USE fldread        ! field read
   USE timing         ! timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   icb_stp        ! routine called in sbcmod.F90 module
   PUBLIC   icb_end        ! routine called in nemogcm.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE icb_stp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_stp  ***
      !!
      !! ** Purpose :   iceberg time stepping.
      !!
      !! ** Method  : - top level routine to do things in the correct order
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time step index
      !
      LOGICAL ::   ll_sample_traj, ll_budget, ll_verbose   ! local logical
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('icb_stp')

      !                       !==  start of timestep housekeeping  ==!
      !
      nktberg = kt
      !
      IF( nn_test_icebergs < 0 .OR. ln_use_calving ) THEN !* read calving data
         !
         CALL fld_read ( kt, 1, sf_icb )
         src_calving     (:,:) = sf_icb(1)%fnow(:,:,1)    ! calving in km^3/year (water equivalent)
         src_calving_hflx(:,:) = 0._wp                    ! NO heat flux for now
         !
      ENDIF
      !
      berg_grid%floating_melt(:,:) = 0._wp
      !
      !                                   !* anything that needs to be reset to zero each timestep 
      CALL icb_dia_step()                 !  for budgets is dealt with here
      !
      !                                   !* write out time
      ll_verbose = .FALSE.
      IF( nn_verbose_write > 0 .AND. MOD( kt-1 , nn_verbose_write ) == 0 )   ll_verbose = ( nn_verbose_level > 0 )
      !
      IF( ll_verbose )   WRITE(numicb,9100) nktberg, ndastp, nsec_day
 9100 FORMAT('kt= ',i8, ' day= ',i8,' secs=',i8)
      !
      !                                   !* copy nemo forcing arrays into iceberg versions with extra halo
      CALL icb_utl_copy()                 ! only necessary for variables not on T points
      !
      !
      !                       !==  process icebergs  ==!
      !                              !
                                     CALL icb_clv_flx( kt )   ! Accumulate ice from calving
      !                              !
                                     CALL icb_clv( kt )       ! Calve excess stored ice into icebergs
      !                              !
      !
      !                       !==  For each berg, evolve  ==!
      !
      IF( ASSOCIATED(first_berg) )   CALL icb_dyn( kt )       ! ice berg dynamics

      IF( lk_mpp ) THEN   ;          CALL icb_lbc_mpp()       ! Send bergs to other PEs
      ELSE                ;          CALL icb_lbc()           ! Deal with any cyclic boundaries in non-mpp case
      ENDIF

      IF( ASSOCIATED(first_berg) )   CALL icb_thm( kt )       ! Ice berg thermodynamics (melting) + rolling
      !
      !
      !                       !==  diagnostics and output  ==!
      !
      !                                   !* For each berg, record trajectory (when needed)
      ll_sample_traj = .FALSE.
      IF( nn_sample_rate > 0 .AND. MOD(kt-1,nn_sample_rate) == 0 )   ll_sample_traj = .TRUE.
      IF( ll_sample_traj .AND. ASSOCIATED(first_berg) )   CALL icb_trj_write( kt )

      !                                   !* Gridded diagnostics
      !                                   !  To get these iom_put's and those preceding to actually do something
      !                                   !  use key_iomput in cpp file and create content for XML file
      !
      CALL iom_put( "calving"           , berg_grid%calving      (:,:)   )  ! 'calving mass input'
      CALL iom_put( "berg_floating_melt", berg_grid%floating_melt(:,:)   )  ! 'Melt rate of icebergs + bits' , 'kg/m2/s'
      CALL iom_put( "berg_stored_ice"   , berg_grid%stored_ice   (:,:,:) )  ! 'Accumulated ice mass by class', 'kg'
      !
      CALL icb_dia_put()                  !* store mean budgets
      !
      !                                   !*  Dump icebergs to screen
      IF( nn_verbose_level >= 2 )   CALL icb_utl_print( 'icb_stp, status', kt )
      !
      !                                   !* Diagnose budgets
      ll_budget = .FALSE.
      IF( nn_verbose_write > 0 .AND. MOD(kt-1,nn_verbose_write) == 0 )   ll_budget = ln_bergdia
      CALL icb_dia( ll_budget )
      !
      IF( lrst_oce ) THEN    !* restart
         CALL icb_rst_write( kt )
         IF( nn_sample_rate > 0 )   CALL icb_trj_sync()
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('icb_stp')
      !
   END SUBROUTINE icb_stp


   SUBROUTINE icb_end( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_end  ***
      !!
      !! ** Purpose :   close iceberg files
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::   kt   ! model time-step index
      !!----------------------------------------------------------------------
      !
      ! finish with trajectories if they were written
      IF( nn_sample_rate > 0 )   CALL icb_trj_end()

      IF(lwp) WRITE(numout,'(a,i6)') 'icebergs: icb_end complete', narea
      !
      IF( nn_verbose_level > 0 ) THEN
         CALL flush( numicb )
         CLOSE( numicb )
      ENDIF
      !
   END SUBROUTINE icb_end

   !!======================================================================
END MODULE icbstp
