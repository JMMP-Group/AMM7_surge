MODULE trdini
   !!======================================================================
   !!                       ***  MODULE  trdini  ***
   !! Ocean diagnostics:  ocean tracers and dynamic trends
   !!=====================================================================
   !! History :   3.5  !  2012-02  (G. Madec) add 3D trends output for T, S, U, V, PE and KE
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   trd_init      : initialization step
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean domain
   USE trd_oce        ! trends: ocean variables
   USE trdken         ! trends: 3D kinetic   energy
   USE trdpen         ! trends: 3D potential energy
   USE trdglo         ! trends: global domain averaged tracers and dynamics
   USE trdmxl         ! trends: mixed layer averaged trends (tracer only)
   USE trdvor         ! trends: vertical averaged vorticity 
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trd_init   ! called by nemogcm.F90 module

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trd_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_init  ***
      !! 
      !! ** Purpose :   Initialization of trend diagnostics
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! local integer
      !!
      NAMELIST/namtrd/ ln_dyn_trd, ln_KE_trd, ln_vor_trd, ln_dyn_mxl,   &
         &             ln_tra_trd, ln_PE_trd, ln_glo_trd, ln_tra_mxl, nn_trd 
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namtrd in reference namelist : trends diagnostic
      READ  ( numnam_ref, namtrd, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtrd in reference namelist' )
      !
      REWIND( numnam_cfg )              ! Namelist namtrd in configuration namelist : trends diagnostic
      READ  ( numnam_cfg, namtrd, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtrd in configuration namelist' )
      IF(lwm) WRITE( numond, namtrd )
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'trd_init : Momentum/Tracers trends'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '   Namelist namtrd : set trends parameters'
         WRITE(numout,*) '      global domain averaged dyn & tra trends   ln_glo_trd  = ', ln_glo_trd
         WRITE(numout,*) '      U & V trends: 3D output                   ln_dyn_trd  = ', ln_dyn_trd
         WRITE(numout,*) '      U & V trends: Mixed Layer averaged        ln_dyn_mxl  = ', ln_dyn_mxl
         WRITE(numout,*) '      T & S trends: 3D output                   ln_tra_trd  = ', ln_tra_trd
         WRITE(numout,*) '      T & S trends: Mixed Layer averaged        ln_tra_mxl  = ', ln_tra_mxl
         WRITE(numout,*) '      Kinetic   Energy trends                   ln_KE_trd   = ', ln_KE_trd
         WRITE(numout,*) '      Potential Energy trends                   ln_PE_trd   = ', ln_PE_trd
         WRITE(numout,*) '      Barotropic vorticity trends               ln_vor_trd  = ', ln_vor_trd
         !
         WRITE(numout,*) '      frequency of trends diagnostics (glo)     nn_trd      = ', nn_trd
      ENDIF
      !
      !                             ! trend extraction flags  
      l_trdtra = .FALSE.                                                       ! tracers  
      IF ( ln_tra_trd .OR. ln_PE_trd .OR. ln_tra_mxl .OR.   &
         & ln_glo_trd                                       )   l_trdtra = .TRUE. 
      !
      l_trddyn = .FALSE.                                                       ! momentum
      IF ( ln_dyn_trd .OR. ln_KE_trd .OR. ln_dyn_mxl .OR.   &
         & ln_vor_trd .OR. ln_glo_trd                       )   l_trddyn = .TRUE.
      !

!!gm check the stop below      
      IF( ln_dyn_mxl )   CALL ctl_stop( 'ML diag on momentum are not yet coded we stop' )
      !

!!gm end
      IF( ln_tra_mxl .OR. ln_vor_trd )   CALL ctl_stop( 'ML tracer and Barotropic vorticity diags are still using old IOIPSL' )
!!gm end
      !
!      IF( .NOT.ln_linssh .AND. ( l_trdtra .OR. l_trddyn ) )  CALL ctl_stop( 'trend diagnostics with variable volume not validated' )
      
!!gm  : Potential BUG : 3D output only for vector invariant form!  add a ctl_stop or code the flux form case
!!gm  : bug/pb for vertical advection of tracer in vvl case: add T.dt[eta] in the output... 

      !                             ! diagnostic initialization  
      IF( ln_glo_trd )   CALL trd_glo_init      ! global domain averaged trends
      IF( ln_tra_mxl )   CALL trd_mxl_init      ! mixed-layer          trends  
      IF( ln_vor_trd )   CALL trd_vor_init      ! barotropic vorticity trends
      IF( ln_KE_trd  )   CALL trd_ken_init      ! 3D Kinetic    energy trends
      IF( ln_PE_trd  )   CALL trd_pen_init      ! 3D Potential  energy trends
      !
   END SUBROUTINE trd_init

   !!======================================================================
END MODULE trdini
