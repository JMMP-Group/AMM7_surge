MODULE dia25h 
   !!======================================================================
   !!                       ***  MODULE  diaharm  ***
   !! Harmonic analysis of tidal constituents 
   !!======================================================================
   !! History :  3.6  !  2014  (E O'Dea)  Original code
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE zdf_oce         ! ocean vertical physics    
   USE zdfgls   , ONLY : hmxl_n
   !
   USE in_out_manager  ! I/O units
   USE iom             ! I/0 library
   USE wet_dry

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_25h_init               ! routine called by nemogcm.F90
   PUBLIC   dia_25h                    ! routine called by diawri.F90

   LOGICAL, PUBLIC ::   ln_dia25h      !:  25h mean output

   ! variables for calculating 25-hourly means
   INTEGER , SAVE ::   cnt_25h           ! Counter for 25 hour means
   REAL(wp), SAVE ::   r1_25 = 0.04_wp   ! =1/25 
   REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   tn_25h  , sn_25h
   REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:)   ::   sshn_25h 
   REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   un_25h  , vn_25h  , wn_25h
   REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   avt_25h , avm_25h
   REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   en_25h  , rmxln_25h

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_25h_init 
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_25h_init  ***
      !!     
      !! ** Purpose: Initialization of 25h mean namelist 
      !!        
      !! ** Method : Read namelist
      !!---------------------------------------------------------------------------
      INTEGER ::   ios                 ! Local integer output status for namelist read
      INTEGER ::   ierror              ! Local integer for memory allocation
      !
      NAMELIST/nam_dia25h/ ln_dia25h
      !!----------------------------------------------------------------------
      !
      REWIND ( numnam_ref )              ! Read Namelist nam_dia25h in reference namelist : 25hour mean diagnostics
      READ   ( numnam_ref, nam_dia25h, IOSTAT=ios, ERR= 901 )
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nam_dia25h in reference namelist' )
      REWIND( numnam_cfg )              ! Namelist nam_dia25h in configuration namelist  25hour diagnostics
      READ  ( numnam_cfg, nam_dia25h, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nam_dia25h in configuration namelist' )
      IF(lwm) WRITE ( numond, nam_dia25h )

      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dia_25h_init : Output 25 hour mean diagnostics'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist nam_dia25h : set 25h outputs '
         WRITE(numout,*) '      Switch for 25h diagnostics (T) or not (F)  ln_dia25h  = ', ln_dia25h
      ENDIF
      IF( .NOT. ln_dia25h )   RETURN
      ! ------------------- !
      ! 1 - Allocate memory !
      ! ------------------- !
      !                                ! ocean arrays
      ALLOCATE( tn_25h (jpi,jpj,jpk), sn_25h (jpi,jpj,jpk), sshn_25h(jpi,jpj)  ,     &
         &      un_25h (jpi,jpj,jpk), vn_25h (jpi,jpj,jpk), wn_25h(jpi,jpj,jpk),     &
         &      avt_25h(jpi,jpj,jpk), avm_25h(jpi,jpj,jpk),                      STAT=ierror )
      IF( ierror > 0 ) THEN
         CALL ctl_stop( 'dia_25h: unable to allocate ocean arrays' )   ;   RETURN
      ENDIF
      IF( ln_zdftke ) THEN             ! TKE physics
         ALLOCATE( en_25h(jpi,jpj,jpk), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'dia_25h: unable to allocate en_25h' )   ;   RETURN
         ENDIF
      ENDIF
      IF( ln_zdfgls ) THEN             ! GLS physics
         ALLOCATE( en_25h(jpi,jpj,jpk), rmxln_25h(jpi,jpj,jpk), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'dia_25h: unable to allocate en_25h and rmxln_25h' )   ;   RETURN
         ENDIF
      ENDIF
      ! ------------------------- !
      ! 2 - Assign Initial Values !
      ! ------------------------- !
      cnt_25h = 1  ! sets the first value of sum at timestep 1 (note - should strictly be at timestep zero so before values used where possible) 
      tn_25h  (:,:,:) = tsb (:,:,:,jp_tem)
      sn_25h  (:,:,:) = tsb (:,:,:,jp_sal)
      sshn_25h(:,:)   = sshb(:,:)
      un_25h  (:,:,:) = ub  (:,:,:)
      vn_25h  (:,:,:) = vb  (:,:,:)
      avt_25h (:,:,:) = avt (:,:,:)
      avm_25h (:,:,:) = avm (:,:,:)
      IF( ln_zdftke ) THEN
         en_25h(:,:,:) = en(:,:,:)
      ENDIF
      IF( ln_zdfgls ) THEN
         en_25h   (:,:,:) = en    (:,:,:)
         rmxln_25h(:,:,:) = hmxl_n(:,:,:)
      ENDIF
#if defined key_si3
      CALL ctl_stop('STOP', 'dia_25h not setup yet to do tidemean ice')
#endif 
      !
   END SUBROUTINE dia_25h_init


   SUBROUTINE dia_25h( kt )  
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dia_25h  ***
      !!         
      !! ** Purpose :   Write diagnostics with M2/S2 tide removed
      !!
      !! ** Method  :   25hr mean outputs for shelf seas
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER ::   ji, jj, jk
      INTEGER                          ::   iyear0, nimonth0,iday0            ! start year,imonth,day
      LOGICAL ::   ll_print = .FALSE.    ! =T print and flush numout
      REAL(wp)                         ::   zsto, zout, zmax, zjulian, zmdi   ! local scalars
      INTEGER                          ::   i_steps                           ! no of timesteps per hour
      REAL(wp), DIMENSION(jpi,jpj    ) ::   zw2d, un_dm, vn_dm                ! workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zw3d                              ! workspace
      REAL(wp), DIMENSION(jpi,jpj,3)   ::   zwtmb                             ! workspace
      !!----------------------------------------------------------------------

      ! 0. Initialisation
      ! -----------------
      ! Define frequency of summing to create 25 h mean
      IF( MOD( 3600,NINT(rdt) ) == 0 ) THEN
         i_steps = 3600/NINT(rdt)
      ELSE
         CALL ctl_stop('STOP', 'dia_wri_tide: timestep must give MOD(3600,rdt) = 0 otherwise no hourly values are possible')
      ENDIF

      ! local variable for debugging
      ll_print = ll_print .AND. lwp

      ! wn_25h could not be initialised in dia_25h_init, so we do it here instead
      IF( kt == nn_it000 ) THEN
         wn_25h(:,:,:) = wn(:,:,:)
      ENDIF

      ! Sum of 25 hourly instantaneous values to give a 25h mean from 24hours every day
      IF( MOD( kt, i_steps ) == 0  .AND. kt /= nn_it000 ) THEN

         IF (lwp) THEN
              WRITE(numout,*) 'dia_wri_tide : Summing instantaneous hourly diagnostics at timestep ',kt
              WRITE(numout,*) '~~~~~~~~~~~~ '
         ENDIF

         tn_25h  (:,:,:)     = tn_25h  (:,:,:) + tsn (:,:,:,jp_tem)
         sn_25h  (:,:,:)     = sn_25h  (:,:,:) + tsn (:,:,:,jp_sal)
         sshn_25h(:,:)       = sshn_25h(:,:)   + sshn(:,:)
         un_25h  (:,:,:)     = un_25h  (:,:,:) + un  (:,:,:)
         vn_25h  (:,:,:)     = vn_25h  (:,:,:) + vn  (:,:,:)
         wn_25h  (:,:,:)     = wn_25h  (:,:,:) + wn  (:,:,:)
         avt_25h (:,:,:)     = avt_25h (:,:,:) + avt (:,:,:)
         avm_25h (:,:,:)     = avm_25h (:,:,:) + avm (:,:,:)
         IF( ln_zdftke ) THEN
            en_25h(:,:,:)    = en_25h  (:,:,:) + en(:,:,:)
         ENDIF
         IF( ln_zdfgls ) THEN
            en_25h   (:,:,:) = en_25h   (:,:,:) + en    (:,:,:)
            rmxln_25h(:,:,:) = rmxln_25h(:,:,:) + hmxl_n(:,:,:)
         ENDIF
         cnt_25h = cnt_25h + 1
         !
         IF (lwp) THEN
            WRITE(numout,*) 'dia_tide : Summed the following number of hourly values so far',cnt_25h
         ENDIF
         !
      ENDIF ! MOD( kt, i_steps ) == 0

      ! Write data for 25 hour mean output streams
      IF( cnt_25h == 25 .AND.  MOD( kt, i_steps*24) == 0 .AND. kt /= nn_it000 ) THEN
         !
         IF(lwp) THEN
            WRITE(numout,*) 'dia_wri_tide : Writing 25 hour mean tide diagnostics at timestep', kt
            WRITE(numout,*) '~~~~~~~~~~~~ '
         ENDIF
         !
         tn_25h  (:,:,:) = tn_25h  (:,:,:) * r1_25
         sn_25h  (:,:,:) = sn_25h  (:,:,:) * r1_25
         sshn_25h(:,:)   = sshn_25h(:,:)   * r1_25
         un_25h  (:,:,:) = un_25h  (:,:,:) * r1_25
         vn_25h  (:,:,:) = vn_25h  (:,:,:) * r1_25
         wn_25h  (:,:,:) = wn_25h  (:,:,:) * r1_25
         avt_25h (:,:,:) = avt_25h (:,:,:) * r1_25
         avm_25h (:,:,:) = avm_25h (:,:,:) * r1_25
         IF( ln_zdftke ) THEN
            en_25h(:,:,:) = en_25h(:,:,:) * r1_25
         ENDIF
         IF( ln_zdfgls ) THEN
            en_25h   (:,:,:) = en_25h   (:,:,:) * r1_25
            rmxln_25h(:,:,:) = rmxln_25h(:,:,:) * r1_25
         ENDIF
         !
         IF(lwp)  WRITE(numout,*) 'dia_wri_tide : Mean calculated by dividing 25 hour sums and writing output'
         zmdi=1.e+20 !missing data indicator for masking
         ! write tracers (instantaneous)
         zw3d(:,:,:) = tn_25h(:,:,:)*tmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
         CALL iom_put("temper25h", zw3d)   ! potential temperature
         zw3d(:,:,:) = sn_25h(:,:,:)*tmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
         CALL iom_put( "salin25h", zw3d  )   ! salinity
         zw2d(:,:) = sshn_25h(:,:)*tmask(:,:,1) + zmdi*(1.0-tmask(:,:,1))
         IF( ll_wd ) THEN
            CALL iom_put( "ssh25h", zw2d+ssh_ref )   ! sea surface 
         ELSE
            CALL iom_put( "ssh25h", zw2d )   ! sea surface
         ENDIF
         ! Write velocities (instantaneous)
         zw3d(:,:,:) = un_25h(:,:,:)*umask(:,:,:) + zmdi*(1.0-umask(:,:,:))
         CALL iom_put("vozocrtx25h", zw3d)    ! i-current
         zw3d(:,:,:) = vn_25h(:,:,:)*vmask(:,:,:) + zmdi*(1.0-vmask(:,:,:))
         CALL iom_put("vomecrty25h", zw3d  )   ! j-current
         zw3d(:,:,:) = wn_25h(:,:,:)*wmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
         CALL iom_put("vovecrtz25h", zw3d )   ! k-current
         ! Write vertical physics
         zw3d(:,:,:) = avt_25h(:,:,:)*wmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
         CALL iom_put("avt25h", zw3d )   ! diffusivity
         zw3d(:,:,:) = avm_25h(:,:,:)*wmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
         CALL iom_put("avm25h", zw3d)   ! viscosity
         IF( ln_zdftke ) THEN
            zw3d(:,:,:) = en_25h(:,:,:)*wmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
            CALL iom_put("tke25h", zw3d)   ! tke
         ENDIF
         IF( ln_zdfgls ) THEN
            zw3d(:,:,:) = en_25h(:,:,:)*wmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
            CALL iom_put("tke25h", zw3d)   ! tke
            zw3d(:,:,:) = rmxln_25h(:,:,:)*wmask(:,:,:) + zmdi*(1.0-tmask(:,:,:))
            CALL iom_put( "mxln25h",zw3d)
         ENDIF
         !
         ! After the write reset the values to cnt=1 and sum values equal current value 
         tn_25h  (:,:,:) = tsn (:,:,:,jp_tem)
         sn_25h  (:,:,:) = tsn (:,:,:,jp_sal)
         sshn_25h(:,:)   = sshn(:,:)
         un_25h  (:,:,:) = un  (:,:,:)
         vn_25h  (:,:,:) = vn  (:,:,:)
         wn_25h  (:,:,:) = wn  (:,:,:)
         avt_25h (:,:,:) = avt (:,:,:)
         avm_25h (:,:,:) = avm (:,:,:)
         IF( ln_zdftke ) THEN
            en_25h(:,:,:) = en(:,:,:)
         ENDIF
         IF( ln_zdfgls ) THEN
            en_25h   (:,:,:) = en    (:,:,:)
            rmxln_25h(:,:,:) = hmxl_n(:,:,:)
         ENDIF
         cnt_25h = 1
         IF(lwp)  WRITE(numout,*) 'dia_wri_tide :   &
            &    After 25hr mean write, reset sum to current value and cnt_25h to one for overlapping average', cnt_25h
      ENDIF !  cnt_25h .EQ. 25 .AND.  MOD( kt, i_steps * 24) == 0 .AND. kt .NE. nn_it000
      !
   END SUBROUTINE dia_25h 

   !!======================================================================
END MODULE dia25h
