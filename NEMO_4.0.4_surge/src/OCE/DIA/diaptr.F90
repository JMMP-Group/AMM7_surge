MODULE diaptr
   !!======================================================================
   !!                       ***  MODULE  diaptr  ***
   !! Ocean physics:  Computes meridonal transports and zonal means
   !!=====================================================================
   !! History :  1.0  ! 2003-09  (C. Talandier, G. Madec)  Original code
   !!            2.0  ! 2006-01  (A. Biastoch)  Allow sub-basins computation
   !!            3.2  ! 2010-03  (O. Marti, S. Flavoni) Add fields
   !!            3.3  ! 2010-10  (G. Madec)  dynamical allocation
   !!            3.6  ! 2014-12  (C. Ethe) use of IOM
   !!            3.6  ! 2016-06  (T. Graham) Addition of diagnostics for CMIP6
   !!            4.0  ! 2010-08  ( C. Ethe, J. Deshayes ) Improvment
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_ptr      : Poleward Transport Diagnostics module
   !!   dia_ptr_init : Initialization, namelist read
   !!   ptr_sjk      : "zonal" mean computation of a field - tracer or flux array
   !!   ptr_sj       : "zonal" and vertical sum computation of a "meridional" flux array
   !!                   (Generic interface to ptr_sj_3d, ptr_sj_2d)
   !!----------------------------------------------------------------------
   USE oce              ! ocean dynamics and active tracers
   USE dom_oce          ! ocean space and time domain
   USE phycst           ! physical constants
   !
   USE iom              ! IOM library
   USE in_out_manager   ! I/O manager
   USE lib_mpp          ! MPP library
   USE timing           ! preformance summary

   IMPLICIT NONE
   PRIVATE

   INTERFACE ptr_sj
      MODULE PROCEDURE ptr_sj_3d, ptr_sj_2d
   END INTERFACE

   PUBLIC   ptr_sj         ! call by tra_ldf & tra_adv routines
   PUBLIC   ptr_sjk        ! 
   PUBLIC   dia_ptr_init   ! call in memogcm
   PUBLIC   dia_ptr        ! call in step module
   PUBLIC   dia_ptr_hst    ! called from tra_ldf/tra_adv routines

   !                                  !!** namelist  namptr  **
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hstr_adv, hstr_ldf, hstr_eiv   !: Heat/Salt TRansports(adv, diff, Bolus.)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hstr_ove, hstr_btr, hstr_vtr   !: heat Salt TRansports(overturn, baro, merional)

   LOGICAL, PUBLIC ::   ln_diaptr   !  Poleward transport flag (T) or not (F)
   LOGICAL, PUBLIC ::   ln_subbas   !  Atlantic/Pacific/Indian basins calculation
   INTEGER, PARAMETER, PUBLIC ::   nptr = 5  ! (glo, atl, pac, ind, ipc)

   REAL(wp) ::   rc_sv    = 1.e-6_wp   ! conversion from m3/s to Sverdrup
   REAL(wp) ::   rc_pwatt = 1.e-15_wp  ! conversion from W    to PW (further x rau0 x Cp)
   REAL(wp) ::   rc_ggram = 1.e-9_wp   ! conversion from g    to Gg  (further x rau0)

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: btmsk   ! T-point basin interior masks
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: btmsk34 ! mask out Southern Ocean (=0 south of 34°S)

   REAL(wp), TARGET, ALLOCATABLE, SAVE, DIMENSION(:)   :: p_fval1d
   REAL(wp), TARGET, ALLOCATABLE, SAVE, DIMENSION(:,:) :: p_fval2d

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_ptr( pvtr )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_ptr  ***
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::   pvtr   ! j-effective transport
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::   zsfc,zvfc               ! local scalar
      REAL(wp), DIMENSION(jpi,jpj)     ::  z2d   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zmask   ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  z3d    ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts) ::  zts   ! 3D workspace
      REAL(wp), DIMENSION(jpj)      ::  zvsum, ztsum, zssum   ! 1D workspace
      !
      !overturning calculation
      REAL(wp), DIMENSION(jpj,jpk,nptr) :: sjk, r1_sjk, v_msf  ! i-mean i-k-surface and its inverse
      REAL(wp), DIMENSION(jpj,jpk,nptr) :: zt_jk, zs_jk ! i-mean T and S, j-Stream-Function

      REAL(wp), DIMENSION(jpi,jpj,jpk,nptr)  :: z4d1, z4d2
      REAL(wp), DIMENSION(jpi,jpj,nptr)      :: z3dtr ! i-mean T and S, j-Stream-Function
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dia_ptr')
      !
      IF( PRESENT( pvtr ) ) THEN
         IF( iom_use( 'zomsf' ) ) THEN    ! effective MSF
            DO jn = 1, nptr                                    ! by sub-basins
               z4d1(1,:,:,jn) =  ptr_sjk( pvtr(:,:,:), btmsk34(:,:,jn) )  ! zonal cumulative effective transport excluding closed seas
               DO jk = jpkm1, 1, -1 
                  z4d1(1,:,jk,jn) = z4d1(1,:,jk+1,jn) - z4d1(1,:,jk,jn)    ! effective j-Stream-Function (MSF)
               END DO
               DO ji = 1, jpi
                  z4d1(ji,:,:,jn) = z4d1(1,:,:,jn)
               ENDDO
            END DO
            CALL iom_put( 'zomsf', z4d1 * rc_sv )
         ENDIF
         IF(  iom_use( 'sopstove' ) .OR. iom_use( 'sophtove' ) .OR.   &
            & iom_use( 'sopstbtr' ) .OR. iom_use( 'sophtbtr' ) ) THEN
            ! define fields multiplied by scalar
            zmask(:,:,:) = 0._wp
            zts(:,:,:,:) = 0._wp
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, jpi
                     zvfc = e1v(ji,jj) * e3v_n(ji,jj,jk)
                     zmask(ji,jj,jk)      = vmask(ji,jj,jk)      * zvfc
                     zts(ji,jj,jk,jp_tem) = (tsn(ji,jj,jk,jp_tem)+tsn(ji,jj+1,jk,jp_tem)) * 0.5 * zvfc  !Tracers averaged onto V grid
                     zts(ji,jj,jk,jp_sal) = (tsn(ji,jj,jk,jp_sal)+tsn(ji,jj+1,jk,jp_sal)) * 0.5 * zvfc
                  ENDDO
               ENDDO
             ENDDO
         ENDIF
         IF( iom_use( 'sopstove' ) .OR. iom_use( 'sophtove' ) ) THEN
            DO jn = 1, nptr
               sjk(:,:,jn) = ptr_sjk( zmask(:,:,:), btmsk(:,:,jn) )
               r1_sjk(:,:,jn) = 0._wp
               WHERE( sjk(:,:,jn) /= 0._wp )   r1_sjk(:,:,jn) = 1._wp / sjk(:,:,jn)
               ! i-mean T and S, j-Stream-Function, basin
               zt_jk(:,:,jn) = ptr_sjk( zts(:,:,:,jp_tem), btmsk(:,:,jn) ) * r1_sjk(:,:,jn)
               zs_jk(:,:,jn) = ptr_sjk( zts(:,:,:,jp_sal), btmsk(:,:,jn) ) * r1_sjk(:,:,jn)
               v_msf(:,:,jn) = ptr_sjk( pvtr(:,:,:), btmsk34(:,:,jn) ) 
               hstr_ove(:,jp_tem,jn) = SUM( v_msf(:,:,jn)*zt_jk(:,:,jn), 2 )
               hstr_ove(:,jp_sal,jn) = SUM( v_msf(:,:,jn)*zs_jk(:,:,jn), 2 )
               !
            ENDDO
            DO jn = 1, nptr
               z3dtr(1,:,jn) = hstr_ove(:,jp_tem,jn) * rc_pwatt  !  (conversion in PW)
               DO ji = 1, jpi
                  z3dtr(ji,:,jn) = z3dtr(1,:,jn)
               ENDDO
            ENDDO
            CALL iom_put( 'sophtove', z3dtr )
            DO jn = 1, nptr
               z3dtr(1,:,jn) = hstr_ove(:,jp_sal,jn) * rc_ggram !  (conversion in Gg)
               DO ji = 1, jpi
                  z3dtr(ji,:,jn) = z3dtr(1,:,jn)
               ENDDO
            ENDDO
            CALL iom_put( 'sopstove', z3dtr )
         ENDIF

         IF( iom_use( 'sopstbtr' ) .OR. iom_use( 'sophtbtr' ) ) THEN
            ! Calculate barotropic heat and salt transport here 
            DO jn = 1, nptr
               sjk(:,1,jn) = ptr_sj( zmask(:,:,:), btmsk(:,:,jn) )
               r1_sjk(:,1,jn) = 0._wp
               WHERE( sjk(:,1,jn) /= 0._wp )   r1_sjk(:,1,jn) = 1._wp / sjk(:,1,jn)
               !
               zvsum(:) = ptr_sj( pvtr(:,:,:), btmsk34(:,:,jn) )
               ztsum(:) = ptr_sj( zts(:,:,:,jp_tem), btmsk(:,:,jn) )
               zssum(:) = ptr_sj( zts(:,:,:,jp_sal), btmsk(:,:,jn) )
               hstr_btr(:,jp_tem,jn) = zvsum(:) * ztsum(:) * r1_sjk(:,1,jn)
               hstr_btr(:,jp_sal,jn) = zvsum(:) * zssum(:) * r1_sjk(:,1,jn)
               !
            ENDDO
            DO jn = 1, nptr
               z3dtr(1,:,jn) = hstr_btr(:,jp_tem,jn) * rc_pwatt  !  (conversion in PW)
               DO ji = 1, jpi
                  z3dtr(ji,:,jn) = z3dtr(1,:,jn)
               ENDDO
            ENDDO
            CALL iom_put( 'sophtbtr', z3dtr )
            DO jn = 1, nptr
               z3dtr(1,:,jn) = hstr_btr(:,jp_sal,jn) * rc_ggram !  (conversion in Gg)
               DO ji = 1, jpi
                  z3dtr(ji,:,jn) = z3dtr(1,:,jn)
               ENDDO
            ENDDO
            CALL iom_put( 'sopstbtr', z3dtr )
         ENDIF 
         !
      ELSE
         !
         zmask(:,:,:) = 0._wp
         zts(:,:,:,:) = 0._wp
         IF( iom_use( 'zotem' ) .OR. iom_use( 'zosal' ) .OR. iom_use( 'zosrf' )  ) THEN    ! i-mean i-k-surface 
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zsfc = e1t(ji,jj) * e3t_n(ji,jj,jk)
                     zmask(ji,jj,jk)      = tmask(ji,jj,jk)      * zsfc
                     zts(ji,jj,jk,jp_tem) = tsn(ji,jj,jk,jp_tem) * zsfc
                     zts(ji,jj,jk,jp_sal) = tsn(ji,jj,jk,jp_sal) * zsfc
                  END DO
               END DO
            END DO
            !
            DO jn = 1, nptr
               zmask(1,:,:) = ptr_sjk( zmask(:,:,:), btmsk(:,:,jn) )
               z4d1(:,:,:,jn) = zmask(:,:,:)
            ENDDO
            CALL iom_put( 'zosrf', z4d1 )
            !
            DO jn = 1, nptr
               z4d2(1,:,:,jn) = ptr_sjk( zts(:,:,:,jp_tem), btmsk(:,:,jn) ) &
                  &            / MAX( z4d1(1,:,:,jn), 10.e-15 )
               DO ji = 1, jpi
                  z4d2(ji,:,:,jn) = z4d2(1,:,:,jn)
               ENDDO
            ENDDO
            CALL iom_put( 'zotem', z4d2 )
            !
            DO jn = 1, nptr
               z4d2(1,:,:,jn) = ptr_sjk( zts(:,:,:,jp_sal), btmsk(:,:,jn) ) &
                  &            / MAX( z4d1(1,:,:,jn), 10.e-15 )
               DO ji = 1, jpi
                  z4d2(ji,:,:,jn) = z4d2(1,:,:,jn)
               ENDDO
            ENDDO
            CALL iom_put( 'zosal', z4d2 )
            !
         ENDIF
         !
         !                                ! Advective and diffusive heat and salt transport
         IF( iom_use( 'sophtadv' ) .OR. iom_use( 'sopstadv' ) ) THEN  
            ! 
            DO jn = 1, nptr
               z3dtr(1,:,jn) = hstr_adv(:,jp_tem,jn) * rc_pwatt  !  (conversion in PW)
               DO ji = 1, jpi
                  z3dtr(ji,:,jn) = z3dtr(1,:,jn)
               ENDDO
            ENDDO
            CALL iom_put( 'sophtadv', z3dtr )
            DO jn = 1, nptr
               z3dtr(1,:,jn) = hstr_adv(:,jp_sal,jn) * rc_ggram !  (conversion in Gg)
               DO ji = 1, jpi
                  z3dtr(ji,:,jn) = z3dtr(1,:,jn)
               ENDDO
            ENDDO
            CALL iom_put( 'sopstadv', z3dtr )
         ENDIF
         !
         IF( iom_use( 'sophtldf' ) .OR. iom_use( 'sopstldf' ) ) THEN  
            ! 
            DO jn = 1, nptr
               z3dtr(1,:,jn) = hstr_ldf(:,jp_tem,jn) * rc_pwatt  !  (conversion in PW)
               DO ji = 1, jpi
                  z3dtr(ji,:,jn) = z3dtr(1,:,jn)
               ENDDO
            ENDDO
            CALL iom_put( 'sophtldf', z3dtr )
            DO jn = 1, nptr
               z3dtr(1,:,jn) = hstr_ldf(:,jp_sal,jn) * rc_ggram !  (conversion in Gg)
               DO ji = 1, jpi
                  z3dtr(ji,:,jn) = z3dtr(1,:,jn)
               ENDDO
            ENDDO
            CALL iom_put( 'sopstldf', z3dtr )
         ENDIF
         !
         IF( iom_use( 'sophteiv' ) .OR. iom_use( 'sopsteiv' ) ) THEN  
            ! 
            DO jn = 1, nptr
               z3dtr(1,:,jn) = hstr_eiv(:,jp_tem,jn) * rc_pwatt  !  (conversion in PW)
               DO ji = 1, jpi
                  z3dtr(ji,:,jn) = z3dtr(1,:,jn)
               ENDDO
            ENDDO
            CALL iom_put( 'sophteiv', z3dtr )
            DO jn = 1, nptr
               z3dtr(1,:,jn) = hstr_eiv(:,jp_sal,jn) * rc_ggram !  (conversion in Gg)
               DO ji = 1, jpi
                  z3dtr(ji,:,jn) = z3dtr(1,:,jn)
               ENDDO
            ENDDO
            CALL iom_put( 'sopsteiv', z3dtr )
         ENDIF
         !
         IF( iom_use( 'sopstvtr' ) .OR. iom_use( 'sophtvtr' ) ) THEN
            zts(:,:,:,:) = 0._wp
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, jpi
                     zvfc = e1v(ji,jj) * e3v_n(ji,jj,jk)
                     zts(ji,jj,jk,jp_tem) = (tsn(ji,jj,jk,jp_tem)+tsn(ji,jj+1,jk,jp_tem)) * 0.5 * zvfc  !Tracers averaged onto V grid
                     zts(ji,jj,jk,jp_sal) = (tsn(ji,jj,jk,jp_sal)+tsn(ji,jj+1,jk,jp_sal)) * 0.5 * zvfc
                  ENDDO
               ENDDO
             ENDDO
             CALL dia_ptr_hst( jp_tem, 'vtr', zts(:,:,:,jp_tem) )
             CALL dia_ptr_hst( jp_sal, 'vtr', zts(:,:,:,jp_sal) )
             DO jn = 1, nptr
                z3dtr(1,:,jn) = hstr_vtr(:,jp_tem,jn) * rc_pwatt  !  (conversion in PW)
                DO ji = 1, jpi
                   z3dtr(ji,:,jn) = z3dtr(1,:,jn)
                ENDDO
             ENDDO
             CALL iom_put( 'sophtvtr', z3dtr )
             DO jn = 1, nptr
               z3dtr(1,:,jn) = hstr_vtr(:,jp_sal,jn) * rc_ggram !  (conversion in Gg)
               DO ji = 1, jpi
                  z3dtr(ji,:,jn) = z3dtr(1,:,jn)
               ENDDO
            ENDDO
            CALL iom_put( 'sopstvtr', z3dtr )
         ENDIF
         !
         IF( iom_use( 'uocetr_vsum_cumul' ) ) THEN
            CALL iom_get_var(  'uocetr_vsum_op', z2d ) ! get uocetr_vsum_op from xml
            z2d(:,:) = ptr_ci_2d( z2d(:,:) )  
            CALL iom_put( 'uocetr_vsum_cumul', z2d )
         ENDIF
         !
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('dia_ptr')
      !
   END SUBROUTINE dia_ptr


   SUBROUTINE dia_ptr_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_ptr_init  ***
      !!                   
      !! ** Purpose :   Initialization, namelist read
      !!----------------------------------------------------------------------
      INTEGER ::  inum, jn, ios, ierr           ! local integers
      !!
      NAMELIST/namptr/ ln_diaptr, ln_subbas
      REAL(wp), DIMENSION(jpi,jpj) :: zmsk
      !!----------------------------------------------------------------------


      REWIND( numnam_ref )              ! Namelist namptr in reference namelist : Poleward transport
      READ  ( numnam_ref, namptr, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namptr in reference namelist' )

      REWIND( numnam_cfg )              ! Namelist namptr in configuration namelist : Poleward transport
      READ  ( numnam_cfg, namptr, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namptr in configuration namelist' )
      IF(lwm) WRITE ( numond, namptr )

      IF(lwp) THEN                     ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dia_ptr_init : poleward transport and msf initialization'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namptr : set ptr parameters'
         WRITE(numout,*) '      Poleward heat & salt transport (T) or not (F)      ln_diaptr  = ', ln_diaptr
      ENDIF

      IF( ln_diaptr ) THEN  
         !
         IF( dia_ptr_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dia_ptr_init : unable to allocate arrays' )

         rc_pwatt = rc_pwatt * rau0_rcp          ! conversion from K.s-1 to PetaWatt
         rc_ggram = rc_ggram * rau0              ! conversion from m3/s to Gg/s

         IF( lk_mpp )   CALL mpp_ini_znl( numout )     ! Define MPI communicator for zonal sum

         btmsk(:,:,1) = tmask_i(:,:)                 
         CALL iom_open( 'subbasins', inum,  ldstop = .FALSE.  )
         CALL iom_get( inum, jpdom_data, 'atlmsk', btmsk(:,:,2) )   ! Atlantic basin
         CALL iom_get( inum, jpdom_data, 'pacmsk', btmsk(:,:,3) )   ! Pacific  basin
         CALL iom_get( inum, jpdom_data, 'indmsk', btmsk(:,:,4) )   ! Indian   basin
         CALL iom_close( inum )
         btmsk(:,:,5) = MAX ( btmsk(:,:,3), btmsk(:,:,4) )          ! Indo-Pacific basin
         DO jn = 2, nptr
            btmsk(:,:,jn) = btmsk(:,:,jn) * tmask_i(:,:)               ! interior domain only
         END DO
         ! JD : modification so that overturning streamfunction is available in Atlantic at 34S to compare with observations
         WHERE( gphit(:,:)*tmask_i(:,:) < -34._wp)
           zmsk(:,:) = 0._wp      ! mask out Southern Ocean
         ELSE WHERE                  
           zmsk(:,:) = ssmask(:,:)
         END WHERE
         btmsk34(:,:,1) = btmsk(:,:,1)                 
         DO jn = 2, nptr
            btmsk34(:,:,jn) = btmsk(:,:,jn) * zmsk(:,:)               ! interior domain only
         ENDDO

         ! Initialise arrays to zero because diatpr is called before they are first calculated
         ! Note that this means diagnostics will not be exactly correct when model run is restarted.
         hstr_adv(:,:,:) = 0._wp           
         hstr_ldf(:,:,:) = 0._wp           
         hstr_eiv(:,:,:) = 0._wp           
         hstr_ove(:,:,:) = 0._wp           
         hstr_btr(:,:,:) = 0._wp           !
         hstr_vtr(:,:,:) = 0._wp           !
         !
      ENDIF 
      ! 
   END SUBROUTINE dia_ptr_init


   SUBROUTINE dia_ptr_hst( ktra, cptr, pva ) 
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_ptr_hst ***
      !!----------------------------------------------------------------------
      !! Wrapper for heat and salt transport calculations to calculate them for each basin
      !! Called from all advection and/or diffusion routines
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in )  :: ktra  ! tracer index
      CHARACTER(len=3)                , INTENT(in)   :: cptr  ! transport type  'adv'/'ldf'/'eiv'
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)   :: pva   ! 3D input array of advection/diffusion
      INTEGER                                        :: jn    !

      !
      IF( cptr == 'adv' ) THEN
         IF( ktra == jp_tem )  THEN
             DO jn = 1, nptr
                hstr_adv(:,jp_tem,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
             ENDDO
         ENDIF
         IF( ktra == jp_sal )  THEN
             DO jn = 1, nptr
                hstr_adv(:,jp_sal,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
             ENDDO
         ENDIF
      ENDIF
      !
      IF( cptr == 'ldf' ) THEN
         IF( ktra == jp_tem )  THEN
             DO jn = 1, nptr
                hstr_ldf(:,jp_tem,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
             ENDDO
         ENDIF
         IF( ktra == jp_sal )  THEN
             DO jn = 1, nptr
                hstr_ldf(:,jp_sal,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
             ENDDO
         ENDIF
      ENDIF
      !
      IF( cptr == 'eiv' ) THEN
         IF( ktra == jp_tem )  THEN
             DO jn = 1, nptr
                hstr_eiv(:,jp_tem,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
             ENDDO
         ENDIF
         IF( ktra == jp_sal )  THEN
             DO jn = 1, nptr
                hstr_eiv(:,jp_sal,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
             ENDDO
         ENDIF
      ENDIF
      !
      IF( cptr == 'vtr' ) THEN
         IF( ktra == jp_tem )  THEN
             DO jn = 1, nptr
                hstr_vtr(:,jp_tem,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
             ENDDO
         ENDIF
         IF( ktra == jp_sal )  THEN
             DO jn = 1, nptr
                hstr_vtr(:,jp_sal,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
             ENDDO
         ENDIF
      ENDIF
      !
   END SUBROUTINE dia_ptr_hst


   FUNCTION dia_ptr_alloc()
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_ptr_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER               ::   dia_ptr_alloc   ! return value
      INTEGER, DIMENSION(3) ::   ierr
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      IF( .NOT. ALLOCATED( btmsk ) ) THEN
         ALLOCATE( btmsk(jpi,jpj,nptr)    , btmsk34(jpi,jpj,nptr),   &
            &      hstr_adv(jpj,jpts,nptr), hstr_eiv(jpj,jpts,nptr), &
            &      hstr_ove(jpj,jpts,nptr), hstr_btr(jpj,jpts,nptr), &
            &      hstr_ldf(jpj,jpts,nptr), hstr_vtr(jpj,jpts,nptr), STAT=ierr(1)  )
            !
         ALLOCATE( p_fval1d(jpj), p_fval2d(jpj,jpk), Stat=ierr(2))
         !
         dia_ptr_alloc = MAXVAL( ierr )
         CALL mpp_sum( 'diaptr', dia_ptr_alloc )
      ENDIF
      !
   END FUNCTION dia_ptr_alloc


   FUNCTION ptr_sj_3d( pva, pmsk )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_sj_3d  ***
      !!
      !! ** Purpose :   i-k sum computation of a j-flux array
      !!
      !! ** Method  : - i-k sum of pva using the interior 2D vmask (vmask_i).
      !!              pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
      !!
      !! ** Action  : - p_fval: i-k-mean poleward flux of pva
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk)  ::   pva   ! mask flux array at V-point
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj)      ::   pmsk   ! Optional 2D basin mask
      !
      INTEGER                  ::   ji, jj, jk   ! dummy loop arguments
      INTEGER                  ::   ijpj         ! ???
      REAL(wp), POINTER, DIMENSION(:) :: p_fval  ! function value
      !!--------------------------------------------------------------------
      !
      p_fval => p_fval1d

      ijpj = jpj
      p_fval(:) = 0._wp
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! Vector opt.
               p_fval(jj) = p_fval(jj) + pva(ji,jj,jk) * pmsk(ji,jj) * tmask_i(ji,jj)
            END DO
         END DO
      END DO
#if defined key_mpp_mpi
      CALL mpp_sum( 'diaptr', p_fval, ijpj, ncomm_znl)
#endif
      !
   END FUNCTION ptr_sj_3d


   FUNCTION ptr_sj_2d( pva, pmsk )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_sj_2d  ***
      !!
      !! ** Purpose :   "zonal" and vertical sum computation of a i-flux array
      !!
      !! ** Method  : - i-k sum of pva using the interior 2D vmask (vmask_i).
      !!      pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
      !!
      !! ** Action  : - p_fval: i-k-mean poleward flux of pva
      !!----------------------------------------------------------------------
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj) ::   pva   ! mask flux array at V-point
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj) ::   pmsk   ! Optional 2D basin mask
      !
      INTEGER                  ::   ji,jj       ! dummy loop arguments
      INTEGER                  ::   ijpj        ! ???
      REAL(wp), POINTER, DIMENSION(:) :: p_fval ! function value
      !!--------------------------------------------------------------------
      ! 
      p_fval => p_fval1d

      ijpj = jpj
      p_fval(:) = 0._wp
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! Vector opt.
            p_fval(jj) = p_fval(jj) + pva(ji,jj) * pmsk(ji,jj) * tmask_i(ji,jj)
         END DO
      END DO
#if defined key_mpp_mpi
      CALL mpp_sum( 'diaptr', p_fval, ijpj, ncomm_znl )
#endif
      ! 
   END FUNCTION ptr_sj_2d

   FUNCTION ptr_ci_2d( pva )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_ci_2d  ***
      !!
      !! ** Purpose :   "meridional" cumulated sum computation of a j-flux array
      !!
      !! ** Method  : - j cumulated sum of pva using the interior 2D vmask (umask_i).
      !!
      !! ** Action  : - p_fval: j-cumulated sum of pva
      !!----------------------------------------------------------------------
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj)  ::   pva   ! mask flux array at V-point
      !
      INTEGER                  ::   ji,jj,jc       ! dummy loop arguments
      INTEGER                  ::   ijpj        ! ??? 
      REAL(wp), DIMENSION(jpi,jpj) :: p_fval ! function value
      !!--------------------------------------------------------------------
      ! 
      ijpj = jpj  ! ???
      p_fval(:,:) = 0._wp
      DO jc = 1, jpnj ! looping over all processors in j axis
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! Vector opt.
               p_fval(ji,jj) = p_fval(ji,jj-1) + pva(ji,jj) * tmask_i(ji,jj)
            END DO
         END DO
         CALL lbc_lnk( 'diaptr', p_fval, 'U', -1. )
      END DO
      ! 
   END FUNCTION ptr_ci_2d



   FUNCTION ptr_sjk( pta, pmsk )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_sjk  ***
      !!
      !! ** Purpose :   i-sum computation of an array
      !!
      !! ** Method  : - i-sum of pva using the interior 2D vmask (vmask_i).
      !!
      !! ** Action  : - p_fval: i-mean poleward flux of pva
      !!----------------------------------------------------------------------
      !!
      IMPLICIT none
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj,jpk) ::   pta    ! mask flux array at V-point
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj)     ::   pmsk   ! Optional 2D basin mask
      !!
      INTEGER                           :: ji, jj, jk ! dummy loop arguments
      REAL(wp), POINTER, DIMENSION(:,:) :: p_fval     ! return function value
#if defined key_mpp_mpi
      INTEGER, DIMENSION(1) ::   ish
      INTEGER, DIMENSION(2) ::   ish2
      INTEGER               ::   ijpjjpk
      REAL(wp), DIMENSION(jpj*jpk) ::   zwork    ! mask flux array at V-point
#endif
      !!--------------------------------------------------------------------
      !
      p_fval => p_fval2d

      p_fval(:,:) = 0._wp
      !
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! Vector opt.
               p_fval(jj,jk) = p_fval(jj,jk) + pta(ji,jj,jk) * pmsk(ji,jj) * tmask_i(ji,jj)
            END DO
         END DO
      END DO
      !
#if defined key_mpp_mpi
      ijpjjpk = jpj*jpk
      ish(1) = ijpjjpk  ;   ish2(1) = jpj   ;   ish2(2) = jpk
      zwork(1:ijpjjpk) = RESHAPE( p_fval, ish )
      CALL mpp_sum( 'diaptr', zwork, ijpjjpk, ncomm_znl )
      p_fval(:,:) = RESHAPE( zwork, ish2 )
#endif
      !
   END FUNCTION ptr_sjk


   !!======================================================================
END MODULE diaptr
