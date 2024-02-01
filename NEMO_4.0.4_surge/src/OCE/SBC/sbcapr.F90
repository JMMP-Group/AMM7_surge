MODULE sbcapr
   !!======================================================================
   !!                       ***  MODULE  sbcapr  ***
   !! Surface module :   atmospheric pressure forcing
   !!======================================================================
   !! History :  3.3  !   2010-09  (J. Chanut, C. Bricaud, G. Madec)  Original code
   !!----------------------------------------------------------------------
   
   !!----------------------------------------------------------------------
   !!   sbc_apr        : read atmospheric pressure in netcdf files 
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition
   USE phycst          ! physical constants
   !
   USE fldread         ! read input fields
   USE in_out_manager  ! I/O manager
   USE lib_fortran     ! distribued memory computing library
   USE iom             ! IOM library
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_apr       ! routine called in sbcmod
   PUBLIC   sbc_apr_init  ! routine called in sbcmod
   
   !                                          !!* namsbc_apr namelist (Atmospheric PRessure) *
   LOGICAL, PUBLIC ::   ln_apr_obc = .false.   !: inverse barometer added to OBC ssh data 
   LOGICAL, PUBLIC ::   ln_ref_apr             !: ref. pressure: global mean Patm (F) or a constant (F)
   REAL(wp)        ::   rn_pref                !  reference atmospheric pressure   [N/m2]

   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) ::   ssh_ib    ! Inverse barometer now    sea surface height   [m]
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) ::   ssh_ibb   ! Inverse barometer before sea surface height   [m]
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) ::   apr       ! atmospheric pressure at kt                 [N/m2]
   
   REAL(wp) ::   tarea                ! whole domain mean masked ocean surface
   REAL(wp) ::   r1_grau              ! = 1.e0 / (grav * rau0)
   
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_apr   ! structure of input fields (file informations, fields read)

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_apr_init
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_apr  ***
      !!
      !! ** Purpose :   read atmospheric pressure fields in netcdf files.
      !!
      !! ** Method  : - Read namelist namsbc_apr
      !!              - Read Patm fields in netcdf files 
      !!              - Compute reference atmospheric pressure
      !!              - Compute inverse barometer ssh
      !! ** action  :   apr      : atmospheric pressure at kt
      !!                ssh_ib   : inverse barometer ssh at kt
      !!---------------------------------------------------------------------
      INTEGER            ::   ierror  ! local integer 
      INTEGER            ::   ios     ! Local integer output status for namelist read
      !!
      CHARACTER(len=100) ::  cn_dir   ! Root directory for location of ssr files
      TYPE(FLD_N)        ::  sn_apr   ! informations about the fields to be read
      LOGICAL            ::  lrxios   ! read restart using XIOS?
      !!
      NAMELIST/namsbc_apr/ cn_dir, sn_apr, ln_ref_apr, rn_pref, ln_apr_obc
      !!----------------------------------------------------------------------
      REWIND( numnam_ref )              ! Namelist namsbc_apr in reference namelist : File for atmospheric pressure forcing
      READ  ( numnam_ref, namsbc_apr, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_apr in reference namelist' )

      REWIND( numnam_cfg )              ! Namelist namsbc_apr in configuration namelist : File for atmospheric pressure forcing
      READ  ( numnam_cfg, namsbc_apr, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_apr in configuration namelist' )
      IF(lwm) WRITE ( numond, namsbc_apr )
      !
      ALLOCATE( sf_apr(1), STAT=ierror )           !* allocate and fill sf_sst (forcing structure) with sn_sst
      IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_apr: unable to allocate sf_apr structure' )
      !
      CALL fld_fill( sf_apr, (/ sn_apr /), cn_dir, 'sbc_apr', 'Atmospheric pressure ', 'namsbc_apr' )
                                ALLOCATE( sf_apr(1)%fnow(jpi,jpj,1)   )
      IF( sn_apr%ln_tint )   ALLOCATE( sf_apr(1)%fdta(jpi,jpj,1,2) )
                             ALLOCATE( ssh_ib(jpi,jpj) , ssh_ibb(jpi,jpj) )
                             ALLOCATE( apr (jpi,jpj) )
      !
      IF( lwp )THEN                                 !* control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namsbc_apr : Atmospheric PRessure as extrenal forcing'
         WRITE(numout,*) '      ref. pressure: global mean Patm (T) or a constant (F)  ln_ref_apr = ', ln_ref_apr
      ENDIF
      !
      IF( ln_ref_apr ) THEN                        !* Compute whole inner domain mean masked ocean surface
         tarea = glob_sum( 'sbcapr', e1e2t(:,:) )
         IF(lwp) WRITE(numout,*) '         Variable ref. Patm computed over a ocean surface of ', tarea*1e-6, 'km2'
      ELSE
         IF(lwp) WRITE(numout,*) '         Reference Patm used : ', rn_pref, ' N/m2'
      ENDIF
      !
      r1_grau = 1.e0 / (grav * rau0)               !* constant for optimization
      !
      !                                            !* control check
      IF ( ln_apr_obc  ) THEN
         IF(lwp) WRITE(numout,*) '         Inverse barometer added to OBC ssh data'
      ENDIF
!jc: stop below should rather be a warning 
      IF( ln_apr_obc .AND. .NOT.ln_apr_dyn   )   &
            CALL ctl_warn( 'sbc_apr: use inverse barometer ssh at open boundary ONLY requires ln_apr_dyn=T' )
      !
      IF( lwxios ) THEN
         CALL iom_set_rstw_var_active('ssh_ibb')
      ENDIF
   END SUBROUTINE sbc_apr_init

   SUBROUTINE sbc_apr( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_apr  ***
      !!
      !! ** Purpose :   read atmospheric pressure fields in netcdf files.
      !!
      !! ** Method  : - Read namelist namsbc_apr
      !!              - Read Patm fields in netcdf files 
      !!              - Compute reference atmospheric pressure
      !!              - Compute inverse barometer ssh
      !! ** action  :   apr      : atmospheric pressure at kt
      !!                ssh_ib   : inverse barometer ssh at kt
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)::   kt   ! ocean time step
      !
      !!----------------------------------------------------------------------

      !                                         ! ========================== !
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN      !    At each sbc time-step   !
         !                                      ! ===========+++============ !
         !
         IF( kt /= nit000 )   ssh_ibb(:,:) = ssh_ib(:,:)    !* Swap of ssh_ib fields
         !
         CALL fld_read( kt, nn_fsbc, sf_apr )               !* input Patm provided at kt + nn_fsbc/2
         !
         !                                                  !* update the reference atmospheric pressure (if necessary)
         IF( ln_ref_apr )   rn_pref = glob_sum( 'sbcapr', sf_apr(1)%fnow(:,:,1) * e1e2t(:,:) ) / tarea
         !
         !                                                  !* Patm related forcing at kt
         ssh_ib(:,:) = - ( sf_apr(1)%fnow(:,:,1) - rn_pref ) * r1_grau    ! equivalent ssh (inverse barometer)
         apr   (:,:) =     sf_apr(1)%fnow(:,:,1)                        ! atmospheric pressure
         !
         CALL iom_put( "ssh_ib", ssh_ib )                   !* output the inverse barometer ssh
      ENDIF

      !                                         ! ---------------------------------------- !
      IF( kt == nit000 ) THEN                   !   set the forcing field at nit000 - 1    !
         !                                      ! ---------------------------------------- !
         !                                            !* Restart: read in restart file
         IF( ln_rstart .AND. iom_varid( numror, 'ssh_ibb', ldstop = .FALSE. ) > 0 ) THEN 
            IF(lwp) WRITE(numout,*) 'sbc_apr:   ssh_ibb read in the restart file'
            CALL iom_get( numror, jpdom_autoglo, 'ssh_ibb', ssh_ibb, ldxios = lrxios )   ! before inv. barometer ssh
            !
         ELSE                                         !* no restart: set from nit000 values
            IF(lwp) WRITE(numout,*) 'sbc_apr:   ssh_ibb set to nit000 values'
            ssh_ibb(:,:) = ssh_ib(:,:)
         ENDIF
      ENDIF
      !                                         ! ---------------------------------------- !
      IF( lrst_oce ) THEN                       !      Write in the ocean restart file     !
         !                                      ! ---------------------------------------- !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbc_apr : ssh_ib written in ocean restart file at it= ', kt,' date= ', ndastp
         IF(lwp) WRITE(numout,*) '~~~~'
         IF( lwxios ) CALL iom_swap(      cwxios_context          )
         CALL iom_rstput( kt, nitrst, numrow, 'ssh_ibb' , ssh_ib, ldxios = lwxios )
         IF( lwxios ) CALL iom_swap(      cxios_context          )
      ENDIF
      !
   END SUBROUTINE sbc_apr
      
   !!======================================================================
END MODULE sbcapr
