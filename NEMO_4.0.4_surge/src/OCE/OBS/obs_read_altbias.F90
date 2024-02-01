MODULE obs_read_altbias
   !!======================================================================
   !!                       ***  MODULE obs_readaltbias  ***
   !! Observation diagnostics: Read the bias for SLA data
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   obs_rea_altbias : Driver for reading altimeter bias
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
      & gphit
   USE oce, ONLY : &           ! Model variables
      & sshn
   USE obs_inter_h2d
   USE obs_utils               ! Various observation tools
   USE obs_inter_sup

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE

   PUBLIC obs_rea_altbias     ! Read the altimeter bias

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE obs_rea_altbias( sladata, k2dint, bias_file )
      !!---------------------------------------------------------------------
      !!
      !!                   *** ROUTINE obs_rea_altbias ***
      !!
      !! ** Purpose : Read from file the bias data 
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!
      !! References :
      !!
      !! History :  
      !!      ! :  2008-02 (D. Lea) Initial version
      !!----------------------------------------------------------------------
      !! * Modules used
      USE iom
      !
      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) :: &
         & sladata       ! SLA data
      INTEGER, INTENT(IN) :: k2dint
      CHARACTER(LEN=128) :: bias_file

      !! * Local declarations

      CHARACTER(LEN=12), PARAMETER :: cpname = 'obs_rea_altbias'

      INTEGER :: jobs         ! Obs loop variable
      INTEGER :: jpialtbias   ! Number of grid point in latitude for the bias
      INTEGER :: jpjaltbias   ! Number of grid point in longitude for the bias
      INTEGER :: iico         ! Grid point indicies
      INTEGER :: ijco
      INTEGER :: i_nx_id      ! Index to read the NetCDF file
      INTEGER :: i_ny_id      ! 
      INTEGER :: i_file_id    ! 
      INTEGER :: i_var_id

      REAL(wp), DIMENSION(1) :: &
         & zext, &
         & zobsmask
      REAL(wp), DIMENSION(2,2,1) :: &
         & zweig
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zmask, &
         & zbias, &
         & zglam, &
         & zgphi
      REAL(wp), DIMENSION(jpi,jpj) ::   z_altbias
      REAL(wp) :: zlam
      REAL(wp) :: zphi
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi, &
         & igrdj
      INTEGER :: numaltbias

      IF(lwp)WRITE(numout,*) 
      IF(lwp)WRITE(numout,*) ' obs_rea_altbias : '
      IF(lwp)WRITE(numout,*) ' ------------- '
      IF(lwp)WRITE(numout,*) '   Read altimeter bias'

      ! Open the file

      z_altbias(:,:)=0.0_wp
      numaltbias=0

      IF(lwp)WRITE(numout,*) 'Opening ',bias_file

      CALL iom_open( bias_file, numaltbias, ldstop=.FALSE. )
      

      IF (numaltbias .GT. 0) THEN      

         ! Get the Alt bias data
         
         CALL iom_get( numaltbias, jpdom_data, 'altbias', z_altbias(:,:), 1 )
         
         ! Close the file
         
         CALL iom_close(numaltbias)     
         
      ELSE

         IF(lwp)WRITE(numout,*) 'no file found'
      
      ENDIF

      ! Intepolate the bias already on the model grid at the observation point
  
      ALLOCATE( &
         & igrdi(2,2,sladata%nsurf), &
         & igrdj(2,2,sladata%nsurf), &
         & zglam(2,2,sladata%nsurf), &
         & zgphi(2,2,sladata%nsurf), &
         & zmask(2,2,sladata%nsurf), &
         & zbias(2,2,sladata%nsurf)  &
         & )
         
      DO jobs = 1, sladata%nsurf

         igrdi(1,1,jobs) = sladata%mi(jobs)-1
         igrdj(1,1,jobs) = sladata%mj(jobs)-1
         igrdi(1,2,jobs) = sladata%mi(jobs)-1
         igrdj(1,2,jobs) = sladata%mj(jobs)
         igrdi(2,1,jobs) = sladata%mi(jobs)
         igrdj(2,1,jobs) = sladata%mj(jobs)-1
         igrdi(2,2,jobs) = sladata%mi(jobs)
         igrdj(2,2,jobs) = sladata%mj(jobs)

      END DO

      CALL obs_int_comm_2d( 2, 2, sladata%nsurf, jpi, jpj, &
         &                  igrdi, igrdj, glamt, zglam )
      CALL obs_int_comm_2d( 2, 2, sladata%nsurf, jpi, jpj, &
         &                  igrdi, igrdj, gphit, zgphi )
      CALL obs_int_comm_2d( 2, 2, sladata%nsurf, jpi, jpj, &
         &                  igrdi, igrdj, tmask(:,:,1), zmask )
      CALL obs_int_comm_2d( 2, 2, sladata%nsurf, jpi, jpj, &
         &                  igrdi, igrdj, z_altbias, zbias )

      DO jobs = 1, sladata%nsurf

         zlam = sladata%rlam(jobs)
         zphi = sladata%rphi(jobs)
         iico = sladata%mi(jobs)
         ijco = sladata%mj(jobs)
            
         CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,         &
            &                   zglam(:,:,jobs), zgphi(:,:,jobs), &
            &                   zmask(:,:,jobs), zweig, zobsmask )
            
         CALL obs_int_h2d( 1, 1,      &
            &              zweig, zbias(:,:,jobs),  zext )

         ! adjust mdt with bias field
         sladata%rext(jobs,2) = sladata%rext(jobs,2) - zext(1)
            
      END DO

      DEALLOCATE( &
         & igrdi, &
         & igrdj, &
         & zglam, &
         & zgphi, &
         & zmask, &
         & zbias  &
         & )
         
   END SUBROUTINE obs_rea_altbias


 
END MODULE obs_read_altbias
