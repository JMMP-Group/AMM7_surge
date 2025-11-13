MODULE obs_readmdt
   !!======================================================================
   !!                       ***  MODULE obs_readmdt  ***
   !! Observation diagnostics: Read the MDT for SLA data (skeleton for now)
   !!======================================================================
   !! History :      ! 2007-03 (K. Mogensen) Initial skeleton version
   !!                ! 2007-04 (E. Remy) migration and improvement from OPAVAR
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   obs_rea_mdt    : Driver for reading MDT
   !!   obs_offset_mdt : Remove the offset between the model MDT and the used one
   !!----------------------------------------------------------------------
   USE par_kind         ! Precision variables
   USE par_oce          ! Domain parameters
   USE in_out_manager   ! I/O manager
   USE obs_surf_def     ! Surface observation definitions
   USE obs_inter_sup    ! Interpolation support routines
   USE obs_inter_h2d    ! 2D interpolation
   USE obs_utils        ! Various observation tools
   USE iom_nf90         ! IOM NetCDF
   USE netcdf           ! NetCDF library
   USE lib_mpp          ! MPP library
   USE dom_oce, ONLY : &                  ! Domain variables
      &                    tmask, tmask_i, e1e2t, gphit, glamt
   USE obs_const, ONLY :   obfillflt      ! Fillvalue
   USE oce      , ONLY :   sshn           ! Model variables

   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   obs_rea_mdt     ! called by dia_obs_init
   PUBLIC   obs_offset_mdt  ! called by obs_rea_mdt

   INTEGER , PUBLIC :: nn_msshc    = 1         ! MDT correction scheme
   REAL(wp), PUBLIC :: rn_mdtcorr   = 1.61_wp  ! User specified MDT correction
   REAL(wp), PUBLIC :: rn_mdtcutoff = 65.0_wp  ! MDT cutoff for computed correction

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE obs_rea_mdt( sladata, k2dint )
      !!---------------------------------------------------------------------
      !!
      !!                   *** ROUTINE obs_rea_mdt ***
      !!
      !! ** Purpose : Read from file the MDT data (skeleton)
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!----------------------------------------------------------------------
      USE iom
      !
      TYPE(obs_surf), INTENT(inout) ::   sladata   ! SLA data
      INTEGER       , INTENT(in)    ::   k2dint    ! ?
      !
      CHARACTER(LEN=12), PARAMETER ::   cpname  = 'obs_rea_mdt'
      CHARACTER(LEN=20), PARAMETER ::   mdtname = 'slaReferenceLevel.nc'

      INTEGER ::   jobs                ! Obs loop variable
      INTEGER ::   jpimdt, jpjmdt      ! Number of grid point in lat/lon for the MDT
      INTEGER ::   iico, ijco          ! Grid point indicies
      INTEGER ::   i_nx_id, i_ny_id, i_file_id, i_var_id, i_stat
      INTEGER ::   nummdt
      !
      REAL(wp), DIMENSION(1)     ::   zext, zobsmask
      REAL(wp), DIMENSION(2,2,1) ::   zweig
      !
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   zmask, zmdtl, zglam, zgphi
      INTEGER , DIMENSION(:,:,:), ALLOCATABLE ::   igrdi, igrdj
      !
      REAL(wp), DIMENSION(jpi,jpj) ::  z_mdt, mdtmask
         
      REAL(wp) :: zlam, zphi, zfill, zinfill    ! local scalar
      !!----------------------------------------------------------------------

      IF(lwp)WRITE(numout,*) 
      IF(lwp)WRITE(numout,*) ' obs_rea_mdt : Read MDT for referencing altimeter anomalies'
      IF(lwp)WRITE(numout,*) ' ------------- '
      CALL FLUSH(numout)

      CALL iom_open( mdtname, nummdt )       ! Open the file
      !                                      ! Get the MDT data
      CALL iom_get ( nummdt, jpdom_data, 'sossheig', z_mdt(:,:), 1 )
      CALL iom_close(nummdt)                 ! Close the file
      
      ! Read in the fill value
      zinfill = 0.0
      i_stat = nf90_open( mdtname, nf90_nowrite, nummdt )
      i_stat = nf90_inq_varid( nummdt, 'sossheig', i_var_id )
      i_stat = nf90_get_att( nummdt, i_var_id, "_FillValue",zinfill)
      zfill = zinfill
      i_stat = nf90_close( nummdt )

      ! setup mask based on tmask and MDT mask
      ! set mask to 0 where the MDT is set to fillvalue
      WHERE(z_mdt(:,:) /= zfill)   ;   mdtmask(:,:) = tmask(:,:,1)
      ELSE WHERE                   ;   mdtmask(:,:) = 0
      END WHERE

      ! Remove the offset between the MDT used with the sla and the model MDT
      IF( nn_msshc == 1 .OR. nn_msshc == 2 ) &
         & CALL obs_offset_mdt( jpi, jpj, z_mdt, zfill )

      ! Intepolate the MDT already on the model grid at the observation point
  
      ALLOCATE( &
         & igrdi(2,2,sladata%nsurf), &
         & igrdj(2,2,sladata%nsurf), &
         & zglam(2,2,sladata%nsurf), &
         & zgphi(2,2,sladata%nsurf), &
         & zmask(2,2,sladata%nsurf), &
         & zmdtl(2,2,sladata%nsurf)  &
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

      CALL obs_int_comm_2d( 2, 2, sladata%nsurf, jpi, jpj, igrdi, igrdj, glamt  , zglam )
      CALL obs_int_comm_2d( 2, 2, sladata%nsurf, jpi, jpj, igrdi, igrdj, gphit  , zgphi )
      CALL obs_int_comm_2d( 2, 2, sladata%nsurf, jpi, jpj, igrdi, igrdj, mdtmask, zmask )
      CALL obs_int_comm_2d( 2, 2, sladata%nsurf, jpi, jpj, igrdi, igrdj, z_mdt  , zmdtl )

      DO jobs = 1, sladata%nsurf
            
         zlam = sladata%rlam(jobs)
         zphi = sladata%rphi(jobs)

         CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,         &
            &                   zglam(:,:,jobs), zgphi(:,:,jobs), &
            &                   zmask(:,:,jobs), zweig, zobsmask )
            
         CALL obs_int_h2d( 1, 1, zweig, zmdtl(:,:,jobs),  zext )
 
         sladata%rext(jobs,2) = zext(1)

! mark any masked data with a QC flag
         IF( zobsmask(1) == 0 )   sladata%nqc(jobs) = IBSET(sladata%nqc(jobs),15)

         END DO
         
      DEALLOCATE( &
         & igrdi, &
         & igrdj, &
         & zglam, &
         & zgphi, &
         & zmask, &
         & zmdtl  &
         & )

      IF(lwp)WRITE(numout,*) ' ------------- '
      !
   END SUBROUTINE obs_rea_mdt


   SUBROUTINE obs_offset_mdt( kpi, kpj, mdt, zfill )
      !!---------------------------------------------------------------------
      !!
      !!                   *** ROUTINE obs_offset_mdt ***
      !!
      !! ** Purpose : Compute a correction term for the MDT on the model grid
      !!             !!!!! IF it is on the model grid
      !!
      !! ** Method  : Compute the mean difference between the model and the 
      !!              used MDT and remove the offset. 
      !!
      !! ** Action  : 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN) ::  kpi, kpj
      REAL(wp), DIMENSION(kpi,kpj), INTENT(INOUT) ::   mdt     ! MDT used on the model grid
      REAL(wp)                    , INTENT(IN   ) ::   zfill 
      ! 
      INTEGER  :: ji, jj
      REAL(wp) :: zdxdy, zarea, zeta1, zeta2, zcorr_mdt, zcorr_bcketa, zcorr     ! local scalar
      REAL(wp), DIMENSION(jpi,jpj) :: zpromsk
      CHARACTER(LEN=14), PARAMETER ::   cpname = 'obs_offset_mdt'
      !!----------------------------------------------------------------------

      !  Initialize the local mask, for domain projection 
      !  Also exclude mdt points which are set to missing data

      DO ji = 1, jpi
        DO jj = 1, jpj
           zpromsk(ji,jj) = tmask_i(ji,jj)
           IF (    ( gphit(ji,jj) .GT.  rn_mdtcutoff ) &
              &.OR.( gphit(ji,jj) .LT. -rn_mdtcutoff ) &
              &.OR.( mdt(ji,jj) .EQ. zfill ) ) &
              &        zpromsk(ji,jj) = 0.0
        END DO
      END DO 

      ! Compute MSSH mean over [0,360] x [-rn_mdtcutoff,rn_mdtcutoff]

      zarea = 0.0
      zeta1 = 0.0
      zeta2 = 0.0

      DO jj = 1, jpj
         DO ji = 1, jpi
          zdxdy = e1e2t(ji,jj) * zpromsk(ji,jj)
          zarea = zarea + zdxdy
          zeta1 = zeta1 + mdt(ji,jj) * zdxdy
          zeta2 = zeta2 + sshn (ji,jj) * zdxdy
        END DO      
      END DO

      CALL mpp_sum( 'obs_readmdt', zeta1 )
      CALL mpp_sum( 'obs_readmdt', zeta2 )
      CALL mpp_sum( 'obs_readmdt', zarea )
      
      zcorr_mdt    = zeta1 / zarea
      zcorr_bcketa = zeta2 / zarea

      !  Define correction term

      zcorr = zcorr_mdt - zcorr_bcketa

      !  Correct spatial mean of the MSSH

      IF( nn_msshc == 1 )   mdt(:,:) = mdt(:,:) - zcorr  

      ! User defined value : 1.6 m for the Rio MDT compared to ORCA2 MDT

      IF( nn_msshc == 2 )   mdt(:,:) = mdt(:,:) - rn_mdtcorr

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' obs_readmdt : rn_mdtcutoff     = ', rn_mdtcutoff
         WRITE(numout,*) ' -----------   zcorr_mdt     = ', zcorr_mdt
         WRITE(numout,*) '               zcorr_bcketa  = ', zcorr_bcketa
         WRITE(numout,*) '               zcorr         = ', zcorr
         WRITE(numout,*) '               nn_msshc        = ', nn_msshc
      ENDIF

      IF ( nn_msshc == 0 ) WRITE(numout,*) '           MSSH correction is not applied'
      IF ( nn_msshc == 1 ) WRITE(numout,*) '           MSSH correction is applied'
      IF ( nn_msshc == 2 ) WRITE(numout,*) '           User defined MSSH correction' 

      !
   END SUBROUTINE obs_offset_mdt
 
   !!======================================================================
END MODULE obs_readmdt
