MODULE step_diu
   !!======================================================================
   !!                       ***  MODULE stp_diu  ***
   !! Time-stepping of diurnal cycle models
   !!======================================================================
   !! History :  3.7  ! 2015-11  (J. While)  Original code

   USE diurnal_bulk    ! diurnal SST bulk routines  (diurnal_sst_takaya routine) 
   USE cool_skin       ! diurnal cool skin correction (diurnal_sst_coolskin routine)   
   USE iom
   USE sbc_oce
   USE sbcmod           ! surface boundary condition       (sbc     routine)
   USE diaobs           ! Observation operator
   USE oce
   USE daymod
   USE restart          ! ocean restart                    (rst_wri routine)
   USE timing           ! Timing
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp_diurnal   ! called by nemogcm.F90 or step.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

   CONTAINS

   SUBROUTINE stp_diurnal( kstp ) 
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index 
      !!---------------------------------------------------------------------- 
      !!                     ***  ROUTINE stp_diurnal  *** 
      !!                       
      !! ** Purpose : - Time stepping of diurnal SST model only 
      !!   
      !! ** Method  : -1- Update forcings and data   
      !!              -2- Update ocean physics   
      !!              -3- Compute the t and s trends   
      !!              -4- Update t and s   
      !!              -5- Compute the momentum trends 
      !!              -6- Update the horizontal velocity 
      !!              -7- Compute the diagnostics variables (rd,N2, div,cur,w) 
      !!              -8- Outputs and diagnostics 
      !!---------------------------------------------------------------------- 
      INTEGER ::   jk       ! dummy loop indices
      INTEGER ::   indic    ! error indicator if < 0 
      REAL(wp), DIMENSION(jpi,jpj) :: z_fvel_bkginc, z_hflux_bkginc     
      !! --------------------------------------------------------------------- 
      
      IF(ln_diurnal_only) THEN
         indic = 0                                 ! reset to no error condition 
         IF( kstp /= nit000 )   CALL day( kstp )   ! Calendar (day was already called at nit000 in day_init) 
 
         CALL iom_setkt( kstp - nit000 + 1, cxios_context )   ! tell iom we are at time step kstp
         IF( ln_crs ) THEN
            CALL iom_setkt( kstp - nit000 + 1, TRIM(cxios_context)//"_crs" ) ! tell iom we are at time step kstp
         ENDIF
       
            CALL sbc    ( kstp )                      ! Sea Boundary Conditions 
      ENDIF
     
      ! Cool skin
      IF( .NOT.ln_diurnal )   CALL ctl_stop( "stp_diurnal: ln_diurnal not set" )
         
      IF( .NOT. ln_blk    )   CALL ctl_stop( "stp_diurnal: diurnal flux processing only implemented for bulk forcing" ) 

      CALL diurnal_sst_coolskin_step( qns, taum, rhop(:,:,1), rdt)

      CALL iom_put( "sst_wl"   , x_dsst               )    ! warm layer (write out before update below).
      CALL iom_put( "sst_cs"   , x_csdsst             )    ! cool skin

      ! Diurnal warm layer model       
      CALL diurnal_sst_takaya_step( kstp, & 
      &    qsr, qns, taum, rhop(:,:,1), rdt) 

      IF( ln_diurnal_only ) THEN
         IF( ln_diaobs )         CALL dia_obs( kstp )         ! obs-minus-model (assimilation) diagnostics (call after dynamics update)
     
         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
         ! Control and restarts 
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
         IF( kstp == nit000   )   CALL iom_close( numror )     ! close input  ocean restart file 
         IF( lrst_oce         )   CALL rst_write    ( kstp )   ! write output ocean restart file
     
         IF( ln_timing .AND.  kstp == nit000  )   CALL timing_reset 
      ENDIF
       
   END SUBROUTINE stp_diurnal  
   
END MODULE step_diu
