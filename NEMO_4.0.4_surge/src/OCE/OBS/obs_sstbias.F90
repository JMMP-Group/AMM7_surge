MODULE obs_sstbias
   !!======================================================================
   !!                       ***  MODULE obs_sstbias  ***
   !! Observation diagnostics: Read the bias for SST data
   !!======================================================================
   !!----------------------------------------------------------------------
   !!   obs_app_sstbias : Driver for reading and applying the SST bias
   !!----------------------------------------------------------------------
   !! * Modules used   
   USE par_kind, ONLY : &       ! Precision variables
      & wp, &
      & dp, &
      & sp
   USE par_oce, ONLY : &        ! Domain parameters
      & jpi, &
      & jpj, &
      & jpim1
   USE in_out_manager, ONLY : & ! I/O manager
      & lwp,    &
      & numout 
   USE obs_surf_def             ! Surface observation definitions
   USE dom_oce, ONLY : &        ! Domain variables
      & tmask, &
      & tmask_i, &
      & e1t,   &
      & e2t,   &
      & gphit, &
      & glamt
   USE oce, ONLY : &           ! Model variables
      & sshn
   USE obs_inter_h2d
   USE obs_utils               ! Various observation tools
   USE obs_inter_sup
   IMPLICIT NONE
   !! * Routine accessibility
   PRIVATE
   PUBLIC obs_app_sstbias     ! Read the altimeter bias
CONTAINS
   SUBROUTINE obs_app_sstbias( sstdata, k2dint, knumtypes, &
                               cl_bias_files )
      !!---------------------------------------------------------------------
      !!
      !!                   *** ROUTINE obs_app_sstbias ***
      !!
      !! ** Purpose : Read SST bias data from files and apply correction to
      !!             observations
      !!
      !! ** Method  :
      !!
      !! ** Action  :
      !!
      !! References :
      !!
      !! History : 
      !!      ! :  2014-08 (J. While) Bias correction code for SST obs,
      !!      !                      based on obs_rea_altbias
      !!----------------------------------------------------------------------
      !! * Modules used
      USE iom
      USE netcdf
      !! * Arguments

      TYPE(obs_surf), INTENT(INOUT) :: sstdata       ! SST data
      INTEGER, INTENT(IN) :: k2dint
      INTEGER, INTENT(IN) :: knumtypes !number of bias types to read in
      CHARACTER(LEN=128), DIMENSION(knumtypes), INTENT(IN) :: &
                          cl_bias_files !List of files to read
      !! * Local declarations
      INTEGER :: jobs         ! Obs loop variable
      INTEGER :: jpisstbias   ! Number of grid point in latitude for the bias
      INTEGER :: jpjsstbias   ! Number of grid point in longitude for the bias
      INTEGER :: iico         ! Grid point indices
      INTEGER :: ijco
      INTEGER :: jt
      INTEGER :: i_nx_id      ! Index to read the NetCDF file
      INTEGER :: i_ny_id      !
      INTEGER :: i_file_id    !
      INTEGER :: i_var_id
      INTEGER, DIMENSION(knumtypes) :: &
         & ibiastypes             ! Array of the bias types in each file
      REAL(wp), DIMENSION(jpi,jpj,knumtypes) :: & 
         & z_sstbias              ! Array to store the SST bias values
      REAL(wp), DIMENSION(jpi,jpj) :: & 
         & z_sstbias_2d           ! Array to store the SST bias values   
      REAL(wp), DIMENSION(1) :: &
         & zext, &
         & zobsmask
      REAL(wp), DIMENSION(2,2,1) :: &
         & zweig
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zmask, &
         & zglam, &
         & zgphi
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zmask_tmp, &
         & zglam_tmp, &
         & zgphi_tmp   
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  zbias   
      REAL(wp) :: zlam
      REAL(wp) :: zphi
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi, &
         & igrdj
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi_tmp, &
         & igrdj_tmp   
      INTEGER :: numsstbias
      INTEGER(KIND=NF90_INT) :: ifile_source
     
      INTEGER :: incfile
      INTEGER :: jtype
      INTEGER :: iret 
      INTEGER :: inumtype
      IF(lwp)WRITE(numout,*) 
      IF(lwp)WRITE(numout,*) 'obs_rea_sstbias : '
      IF(lwp)WRITE(numout,*) '----------------- '
      IF(lwp)WRITE(numout,*) 'Read SST bias '
      ! Open and read the files
      z_sstbias(:,:,:)=0.0_wp
      DO jtype = 1, knumtypes
     
         numsstbias=0
         IF(lwp)WRITE(numout,*) 'Opening ',cl_bias_files(jtype)
         CALL iom_open( cl_bias_files(jtype), numsstbias, ldstop=.FALSE. )       
         IF (numsstbias > 0) THEN
     
            !Read the bias type from the file
            !No IOM get attribute command at time of writing, 
            !so have to use NETCDF
            !routines directly - should be upgraded in the future
            iret=NF90_OPEN(TRIM(cl_bias_files(jtype)), NF90_NOWRITE, incfile)
            iret=NF90_GET_ATT( incfile, NF90_GLOBAL, "SST_source", &
                              ifile_source )
            ibiastypes(jtype) = ifile_source                 
            iret=NF90_CLOSE(incfile)       
           
            IF ( iret /= 0  ) CALL ctl_stop( &
               'obs_rea_sstbias : Cannot read bias type from file '// &
               cl_bias_files(jtype) )
            ! Get the SST bias data
            CALL iom_get( numsstbias, jpdom_data, 'tn', z_sstbias_2d(:,:), 1 )
            z_sstbias(:,:,jtype) = z_sstbias_2d(:,:)       
            ! Close the file
            CALL iom_close(numsstbias)       
         ELSE
            CALL ctl_stop('obs_read_sstbias: File '// & 
                           TRIM( cl_bias_files(jtype) )//' Not found')
         ENDIF
      END DO
           
      ! Interpolate the bias already on the model grid at the observation point
      ALLOCATE( &
         & igrdi(2,2,sstdata%nsurf), &
         & igrdj(2,2,sstdata%nsurf), &
         & zglam(2,2,sstdata%nsurf), &
         & zgphi(2,2,sstdata%nsurf), &
         & zmask(2,2,sstdata%nsurf)  )
       
      DO jobs = 1, sstdata%nsurf 
         igrdi(1,1,jobs) = sstdata%mi(jobs)-1
         igrdj(1,1,jobs) = sstdata%mj(jobs)-1
         igrdi(1,2,jobs) = sstdata%mi(jobs)-1
         igrdj(1,2,jobs) = sstdata%mj(jobs)
         igrdi(2,1,jobs) = sstdata%mi(jobs)
         igrdj(2,1,jobs) = sstdata%mj(jobs)-1
         igrdi(2,2,jobs) = sstdata%mi(jobs)
         igrdj(2,2,jobs) = sstdata%mj(jobs)
      END DO
      CALL obs_int_comm_2d( 2, 2, sstdata%nsurf, jpi, jpj, &
         &                  igrdi, igrdj, glamt, zglam )
      CALL obs_int_comm_2d( 2, 2, sstdata%nsurf, jpi, jpj, &
         &                  igrdi, igrdj, gphit, zgphi )
      CALL obs_int_comm_2d( 2, 2, sstdata%nsurf, jpi, jpj, &
         &                  igrdi, igrdj, tmask(:,:,1), zmask )
      DO jtype = 1, knumtypes
         
         !Find the number observations of type and allocate tempory arrays
         inumtype = COUNT( sstdata%ntyp(:) == ibiastypes(jtype) )
         ALLOCATE( &
            & igrdi_tmp(2,2,inumtype), &
            & igrdj_tmp(2,2,inumtype), &
            & zglam_tmp(2,2,inumtype), &
            & zgphi_tmp(2,2,inumtype), &
            & zmask_tmp(2,2,inumtype), &
            & zbias( 2,2,inumtype ) )
         jt=1
         DO jobs = 1, sstdata%nsurf 
            IF ( sstdata%ntyp(jobs) == ibiastypes(jtype) ) THEN
               igrdi_tmp(:,:,jt) = igrdi(:,:,jobs) 
               igrdj_tmp(:,:,jt) = igrdj(:,:,jobs)
               zglam_tmp(:,:,jt) = zglam(:,:,jobs)
               zgphi_tmp(:,:,jt) = zgphi(:,:,jobs)
               zmask_tmp(:,:,jt) = zmask(:,:,jobs)
               jt = jt +1
            ENDIF
         END DO
                         
         CALL obs_int_comm_2d( 2, 2, inumtype, jpi, jpj, &
               &           igrdi_tmp(:,:,:), igrdj_tmp(:,:,:), &
               &           z_sstbias(:,:,jtype), zbias(:,:,:) )
         jt=1
         DO jobs = 1, sstdata%nsurf
            IF ( sstdata%ntyp(jobs) == ibiastypes(jtype) ) THEN
               zlam = sstdata%rlam(jobs)
               zphi = sstdata%rphi(jobs)
               iico = sstdata%mi(jobs)
               ijco = sstdata%mj(jobs)         
               CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,         &
                  &                   zglam_tmp(:,:,jt), &
                  &                   zgphi_tmp(:,:,jt), &
                  &                   zmask_tmp(:,:,jt), zweig, zobsmask )
               CALL obs_int_h2d( 1, 1, zweig, zbias(:,:,jt),  zext )
               ! adjust sst with bias field
               sstdata%robs(jobs,1) = sstdata%robs(jobs,1) - zext(1)
               jt=jt+1
            ENDIF
         END DO 
               
         !Deallocate arrays
         DEALLOCATE( &
         & igrdi_tmp, &
         & igrdj_tmp, &
         & zglam_tmp, &
         & zgphi_tmp, &
         & zmask_tmp, &
         & zbias )           
      END DO
      DEALLOCATE( &
         & igrdi, &
         & igrdj, &
         & zglam, &
         & zgphi, &
         & zmask )

      IF(lwp) THEN
         WRITE(numout,*) " "
         WRITE(numout,*) "SST bias correction applied successfully"
         WRITE(numout,*) "Obs types: ",ibiastypes(:), &
                              " Have all been bias corrected\n"
      ENDIF
   END SUBROUTINE obs_app_sstbias
 
END MODULE obs_sstbias
