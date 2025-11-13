MODULE diaobs
   !!======================================================================
   !!                       ***  MODULE diaobs  ***
   !! Observation diagnostics: Computation of the misfit between data and
   !!                          their model equivalent 
   !!======================================================================
   !! History :  1.0  !  2006-03  (K. Mogensen) Original code
   !!             -   !  2006-05  (K. Mogensen, A. Weaver) Reformatted
   !!             -   !  2006-10  (A. Weaver) Cleaning and add controls
   !!             -   !  2007-03  (K. Mogensen) General handling of profiles
   !!             -   !  2007-04  (G. Smith) Generalized surface operators
   !!            2.0  !  2008-10  (M. Valdivieso) obs operator for velocity profiles
   !!            3.4  !  2014-08  (J. While) observation operator for profiles in all vertical coordinates
   !!             -   !                      Incorporated SST bias correction  
   !!            3.6  !  2015-02  (M. Martin) Simplification of namelist and code
   !!             -   !  2015-08  (M. Martin) Combined surface/profile routines.
   !!            4.0  !  2017-11  (G. Madec) style only
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_obs_init  : Reading and prepare observations
   !!   dia_obs       : Compute model equivalent to observations
   !!   dia_obs_wri   : Write observational diagnostics
   !!   calc_date     : Compute the date of timestep in YYYYMMDD.HHMMSS format
   !!   ini_date      : Compute the initial date YYYYMMDD.HHMMSS
   !!   fin_date      : Compute the final date YYYYMMDD.HHMMSS
   !!----------------------------------------------------------------------
   USE par_kind       ! Precision variables
   USE in_out_manager ! I/O manager
   USE par_oce        ! ocean parameter
   USE dom_oce        ! Ocean space and time domain variables
   USE sbc_oce        ! Sea-ice fraction
   !
   USE obs_read_prof  ! Reading and allocation of profile obs
   USE obs_read_surf  ! Reading and allocation of surface obs
   USE obs_sstbias    ! Bias correction routine for SST 
   USE obs_readmdt    ! Reading and allocation of MDT for SLA.
   USE obs_prep       ! Preparation of obs. (grid search etc).
   USE obs_oper       ! Observation operators
   USE obs_write      ! Writing of observation related diagnostics
   USE obs_grid       ! Grid searching
   USE obs_read_altbias ! Bias treatment for altimeter
   USE obs_profiles_def ! Profile data definitions
   USE obs_surf_def   ! Surface data definitions
   USE obs_types      ! Definitions for observation types
   !
   USE mpp_map        ! MPP mapping
   USE lib_mpp        ! For ctl_warn/stop

   IMPLICIT NONE
   PRIVATE

   PUBLIC dia_obs_init     ! Initialize and read observations
   PUBLIC dia_obs          ! Compute model equivalent to observations
   PUBLIC dia_obs_wri      ! Write model equivalent to observations
   PUBLIC dia_obs_dealloc  ! Deallocate dia_obs data
   PUBLIC calc_date        ! Compute the date of a timestep

   LOGICAL, PUBLIC :: ln_diaobs          !: Logical switch for the obs operator
   LOGICAL         :: ln_sstnight        !  Logical switch for night mean SST obs
   LOGICAL         :: ln_sla_fp_indegs   !  T=> SLA obs footprint size specified in degrees, F=> in metres
   LOGICAL         :: ln_sst_fp_indegs   !  T=> SST obs footprint size specified in degrees, F=> in metres
   LOGICAL         :: ln_sss_fp_indegs   !  T=> SSS obs footprint size specified in degrees, F=> in metres
   LOGICAL         :: ln_sic_fp_indegs   !  T=> sea-ice obs footprint size specified in degrees, F=> in metres

   REAL(wp) ::   rn_sla_avglamscl   ! E/W diameter of SLA observation footprint (metres)
   REAL(wp) ::   rn_sla_avgphiscl   ! N/S diameter of SLA observation footprint (metres)
   REAL(wp) ::   rn_sst_avglamscl   ! E/W diameter of SST observation footprint (metres)
   REAL(wp) ::   rn_sst_avgphiscl   ! N/S diameter of SST observation footprint (metres)
   REAL(wp) ::   rn_sss_avglamscl   ! E/W diameter of SSS observation footprint (metres)
   REAL(wp) ::   rn_sss_avgphiscl   ! N/S diameter of SSS observation footprint (metres)
   REAL(wp) ::   rn_sic_avglamscl   ! E/W diameter of sea-ice observation footprint (metres)
   REAL(wp) ::   rn_sic_avgphiscl   ! N/S diameter of sea-ice observation footprint (metres)

   INTEGER :: nn_1dint       ! Vertical interpolation method
   INTEGER :: nn_2dint       ! Default horizontal interpolation method
   INTEGER :: nn_2dint_sla   ! SLA horizontal interpolation method 
   INTEGER :: nn_2dint_sst   ! SST horizontal interpolation method 
   INTEGER :: nn_2dint_sss   ! SSS horizontal interpolation method 
   INTEGER :: nn_2dint_sic   ! Seaice horizontal interpolation method 
   INTEGER, DIMENSION(imaxavtypes) ::   nn_profdavtypes   ! Profile data types representing a daily average
   INTEGER :: nproftypes     ! Number of profile obs types
   INTEGER :: nsurftypes     ! Number of surface obs types
   INTEGER , DIMENSION(:), ALLOCATABLE ::   nvarsprof, nvarssurf   ! Number of profile & surface variables
   INTEGER , DIMENSION(:), ALLOCATABLE ::   nextrprof, nextrsurf   ! Number of profile & surface extra variables
   INTEGER , DIMENSION(:), ALLOCATABLE ::   n2dintsurf             ! Interpolation option for surface variables
   REAL(wp), DIMENSION(:), ALLOCATABLE ::   zavglamscl, zavgphiscl ! E/W & N/S diameter of averaging footprint for surface variables
   LOGICAL , DIMENSION(:), ALLOCATABLE ::   lfpindegs              ! T=> surface obs footprint size specified in degrees, F=> in metres
   LOGICAL , DIMENSION(:), ALLOCATABLE ::   llnightav              ! Logical for calculating night-time averages

   TYPE(obs_surf), PUBLIC, POINTER, DIMENSION(:) ::   surfdata     !: Initial surface data
   TYPE(obs_surf), PUBLIC, POINTER, DIMENSION(:) ::   surfdataqc   !: Surface data after quality control
   TYPE(obs_prof), PUBLIC, POINTER, DIMENSION(:) ::   profdata     !: Initial profile data
   TYPE(obs_prof), PUBLIC, POINTER, DIMENSION(:) ::   profdataqc   !: Profile data after quality control

   CHARACTER(len=6), PUBLIC, DIMENSION(:), ALLOCATABLE ::   cobstypesprof, cobstypessurf   !: Profile & surface obs types

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_obs_init
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_obs_init  ***
      !!          
      !! ** Purpose : Initialize and read observations
      !!
      !! ** Method  : Read the namelist and call reading routines
      !!
      !! ** Action  : Read the namelist and call reading routines
      !!
      !!----------------------------------------------------------------------
      INTEGER, PARAMETER ::   jpmaxnfiles = 1000    ! Maximum number of files for each obs type
      INTEGER, DIMENSION(:), ALLOCATABLE ::   ifilesprof, ifilessurf   ! Number of profile & surface files
      INTEGER :: ios             ! Local integer output status for namelist read
      INTEGER :: jtype           ! Counter for obs types
      INTEGER :: jvar            ! Counter for variables
      INTEGER :: jfile           ! Counter for files
      INTEGER :: jnumsstbias
      !
      CHARACTER(len=128), DIMENSION(jpmaxnfiles) :: &
         & cn_profbfiles, &      ! T/S profile input filenames
         & cn_sstfbfiles, &      ! Sea surface temperature input filenames
         & cn_sssfbfiles, &      ! Sea surface salinity input filenames
         & cn_slafbfiles, &      ! Sea level anomaly input filenames
         & cn_sicfbfiles, &      ! Seaice concentration input filenames
         & cn_velfbfiles, &      ! Velocity profile input filenames
         & cn_sstbiasfiles      ! SST bias input filenames
      CHARACTER(LEN=128) :: &
         & cn_altbiasfile        ! Altimeter bias input filename
      CHARACTER(len=128), DIMENSION(:,:), ALLOCATABLE :: &
         & clproffiles, &        ! Profile filenames
         & clsurffiles           ! Surface filenames
         !
      LOGICAL :: ln_t3d          ! Logical switch for temperature profiles
      LOGICAL :: ln_s3d          ! Logical switch for salinity profiles
      LOGICAL :: ln_sla          ! Logical switch for sea level anomalies 
      LOGICAL :: ln_sst          ! Logical switch for sea surface temperature
      LOGICAL :: ln_sss          ! Logical switch for sea surface salinity
      LOGICAL :: ln_sic          ! Logical switch for sea ice concentration
      LOGICAL :: ln_vel3d        ! Logical switch for velocity (u,v) obs
      LOGICAL :: ln_nea          ! Logical switch to remove obs near land
      LOGICAL :: ln_altbias      ! Logical switch for altimeter bias
      LOGICAL :: ln_sstbias      ! Logical switch for bias corection of SST 
      LOGICAL :: ln_ignmis       ! Logical switch for ignoring missing files
      LOGICAL :: ln_s_at_t       ! Logical switch to compute model S at T obs
      LOGICAL :: ln_bound_reject ! Logical to remove obs near boundaries in LAMs.
      LOGICAL :: llvar1          ! Logical for profile variable 1
      LOGICAL :: llvar2          ! Logical for profile variable 1
      LOGICAL, DIMENSION(jpmaxnfiles) :: lmask ! Used for finding number of sstbias files
      !
      REAL(dp) :: rn_dobsini     ! Obs window start date YYYYMMDD.HHMMSS
      REAL(dp) :: rn_dobsend     ! Obs window end date   YYYYMMDD.HHMMSS
      REAL(wp), DIMENSION(jpi,jpj)     ::   zglam1, zglam2   ! Model longitudes for profile variable 1 & 2
      REAL(wp), DIMENSION(jpi,jpj)     ::   zgphi1, zgphi2   ! Model latitudes  for profile variable 1 & 2
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zmask1, zmask2   ! Model land/sea mask associated with variable 1 & 2
      !!
      NAMELIST/namobs/ln_diaobs, ln_t3d, ln_s3d, ln_sla,              &
         &            ln_sst, ln_sic, ln_sss, ln_vel3d,               &
         &            ln_altbias, ln_sstbias, ln_nea,                 &
         &            ln_grid_global, ln_grid_search_lookup,          &
         &            ln_ignmis, ln_s_at_t, ln_bound_reject,          &
         &            ln_sstnight,                                    &
         &            ln_sla_fp_indegs, ln_sst_fp_indegs,             &
         &            ln_sss_fp_indegs, ln_sic_fp_indegs,             &
         &            cn_profbfiles, cn_slafbfiles,                   &
         &            cn_sstfbfiles, cn_sicfbfiles,                   &
         &            cn_velfbfiles, cn_sssfbfiles,                   &
         &            cn_sstbiasfiles, cn_altbiasfile,                &
         &            cn_gridsearchfile, rn_gridsearchres,            &
         &            rn_dobsini, rn_dobsend,                         &
         &            rn_sla_avglamscl, rn_sla_avgphiscl,             &
         &            rn_sst_avglamscl, rn_sst_avgphiscl,             &
         &            rn_sss_avglamscl, rn_sss_avgphiscl,             &
         &            rn_sic_avglamscl, rn_sic_avgphiscl,             &
         &            nn_1dint, nn_2dint,                             &
         &            nn_2dint_sla, nn_2dint_sst,                     &
         &            nn_2dint_sss, nn_2dint_sic,                     &
         &            nn_msshc, rn_mdtcorr, rn_mdtcutoff,             &
         &            nn_profdavtypes
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! Read namelist parameters
      !-----------------------------------------------------------------------
      ! Some namelist arrays need initialising
      cn_profbfiles  (:) = ''
      cn_slafbfiles  (:) = ''
      cn_sstfbfiles  (:) = ''
      cn_sicfbfiles  (:) = ''
      cn_velfbfiles  (:) = ''
      cn_sssfbfiles  (:) = ''
      cn_sstbiasfiles(:) = ''
      nn_profdavtypes(:) = -1

      CALL ini_date( rn_dobsini )
      CALL fin_date( rn_dobsend )

      ! Read namelist namobs : control observation diagnostics
      REWIND( numnam_ref )   ! Namelist namobs in reference namelist
      READ  ( numnam_ref, namobs, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namobs in reference namelist' )
      REWIND( numnam_cfg )   ! Namelist namobs in configuration namelist
      READ  ( numnam_cfg, namobs, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namobs in configuration namelist' )
      IF(lwm) WRITE ( numond, namobs )

      IF( .NOT.ln_diaobs ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dia_obs_init : NO Observation diagnostic used'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
         RETURN
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_obs_init : Observation diagnostic initialization'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namobs : set observation diagnostic parameters' 
         WRITE(numout,*) '      Logical switch for T profile observations                ln_t3d = ', ln_t3d
         WRITE(numout,*) '      Logical switch for S profile observations                ln_s3d = ', ln_s3d
         WRITE(numout,*) '      Logical switch for SLA observations                      ln_sla = ', ln_sla
         WRITE(numout,*) '      Logical switch for SST observations                      ln_sst = ', ln_sst
         WRITE(numout,*) '      Logical switch for Sea Ice observations                  ln_sic = ', ln_sic
         WRITE(numout,*) '      Logical switch for velocity observations               ln_vel3d = ', ln_vel3d
         WRITE(numout,*) '      Logical switch for SSS observations                      ln_sss = ', ln_sss
         WRITE(numout,*) '      Global distribution of observations              ln_grid_global = ', ln_grid_global
         WRITE(numout,*) '      Logical switch for obs grid search lookup ln_grid_search_lookup = ', ln_grid_search_lookup
         IF (ln_grid_search_lookup) &
            WRITE(numout,*) '      Grid search lookup file header                cn_gridsearchfile = ', cn_gridsearchfile
         WRITE(numout,*) '      Initial date in window YYYYMMDD.HHMMSS               rn_dobsini = ', rn_dobsini
         WRITE(numout,*) '      Final date in window YYYYMMDD.HHMMSS                 rn_dobsend = ', rn_dobsend
         WRITE(numout,*) '      Type of vertical interpolation method                  nn_1dint = ', nn_1dint
         WRITE(numout,*) '      Type of horizontal interpolation method                nn_2dint = ', nn_2dint
         WRITE(numout,*) '      Rejection of observations near land switch               ln_nea = ', ln_nea
         WRITE(numout,*) '      Rejection of obs near open bdys                 ln_bound_reject = ', ln_bound_reject
         WRITE(numout,*) '      MSSH correction scheme                                 nn_msshc = ', nn_msshc
         WRITE(numout,*) '      MDT  correction                                      rn_mdtcorr = ', rn_mdtcorr
         WRITE(numout,*) '      MDT cutoff for computed correction                 rn_mdtcutoff = ', rn_mdtcutoff
         WRITE(numout,*) '      Logical switch for alt bias                          ln_altbias = ', ln_altbias
         WRITE(numout,*) '      Logical switch for sst bias                          ln_sstbias = ', ln_sstbias
         WRITE(numout,*) '      Logical switch for ignoring missing files             ln_ignmis = ', ln_ignmis
         WRITE(numout,*) '      Daily average types                             nn_profdavtypes = ', nn_profdavtypes
         WRITE(numout,*) '      Logical switch for night-time SST obs               ln_sstnight = ', ln_sstnight
      ENDIF
      !-----------------------------------------------------------------------
      ! Set up list of observation types to be used
      ! and the files associated with each type
      !-----------------------------------------------------------------------

      nproftypes = COUNT( (/ln_t3d .OR. ln_s3d, ln_vel3d /) )
      nsurftypes = COUNT( (/ln_sla, ln_sst, ln_sic, ln_sss /) )

      IF( ln_sstbias ) THEN 
         lmask(:) = .FALSE. 
         WHERE( cn_sstbiasfiles(:) /= '' )   lmask(:) = .TRUE. 
         jnumsstbias = COUNT(lmask) 
         lmask(:) = .FALSE. 
      ENDIF      

      IF( nproftypes == 0 .AND. nsurftypes == 0 ) THEN
         CALL ctl_warn( 'dia_obs_init: ln_diaobs is set to true, but all obs operator logical flags',   &
            &           ' (ln_t3d, ln_s3d, ln_sla, ln_sst, ln_sic, ln_vel3d)',                          &
            &           ' are set to .FALSE. so turning off calls to dia_obs'  )
         ln_diaobs = .FALSE.
         RETURN
      ENDIF

      IF( nproftypes > 0 ) THEN
         !
         ALLOCATE( cobstypesprof(nproftypes)             )
         ALLOCATE( ifilesprof   (nproftypes)             )
         ALLOCATE( clproffiles  (nproftypes,jpmaxnfiles) )
         !
         jtype = 0
         IF( ln_t3d .OR. ln_s3d ) THEN
            jtype = jtype + 1
            CALL obs_settypefiles( nproftypes, jpmaxnfiles, jtype, 'prof  ', &
               &                   cn_profbfiles, ifilesprof, cobstypesprof, clproffiles )
         ENDIF
         IF( ln_vel3d ) THEN
            jtype = jtype + 1
            CALL obs_settypefiles( nproftypes, jpmaxnfiles, jtype, 'vel   ', &
               &                   cn_velfbfiles, ifilesprof, cobstypesprof, clproffiles )
         ENDIF
         !
      ENDIF

      IF( nsurftypes > 0 ) THEN
         !
         ALLOCATE( cobstypessurf(nsurftypes)             )
         ALLOCATE( ifilessurf   (nsurftypes)             )
         ALLOCATE( clsurffiles  (nsurftypes,jpmaxnfiles) )
         ALLOCATE( n2dintsurf   (nsurftypes)             )
         ALLOCATE( zavglamscl   (nsurftypes)             )
         ALLOCATE( zavgphiscl   (nsurftypes)             )
         ALLOCATE( lfpindegs    (nsurftypes)             )
         ALLOCATE( llnightav    (nsurftypes)             )
         !
         jtype = 0
         IF( ln_sla ) THEN
            jtype = jtype + 1
            CALL obs_settypefiles( nsurftypes, jpmaxnfiles, jtype, 'sla   ', &
               &                   cn_slafbfiles, ifilessurf, cobstypessurf, clsurffiles )
            CALL obs_setinterpopts( nsurftypes, jtype, 'sla   ',      &
               &                  nn_2dint, nn_2dint_sla,             &
               &                  rn_sla_avglamscl, rn_sla_avgphiscl, &
               &                  ln_sla_fp_indegs, .FALSE.,          &
               &                  n2dintsurf, zavglamscl, zavgphiscl, &
               &                  lfpindegs, llnightav )
         ENDIF
         IF( ln_sst ) THEN
            jtype = jtype + 1
            CALL obs_settypefiles( nsurftypes, jpmaxnfiles, jtype, 'sst   ', &
               &                   cn_sstfbfiles, ifilessurf, cobstypessurf, clsurffiles )
            CALL obs_setinterpopts( nsurftypes, jtype, 'sst   ',      &
               &                  nn_2dint, nn_2dint_sst,             &
               &                  rn_sst_avglamscl, rn_sst_avgphiscl, &
               &                  ln_sst_fp_indegs, ln_sstnight,      &
               &                  n2dintsurf, zavglamscl, zavgphiscl, &
               &                  lfpindegs, llnightav )
         ENDIF
#if defined key_si3 || defined key_cice
         IF( ln_sic ) THEN
            jtype = jtype + 1
            CALL obs_settypefiles( nsurftypes, jpmaxnfiles, jtype, 'sic   ', &
               &                   cn_sicfbfiles, ifilessurf, cobstypessurf, clsurffiles )
            CALL obs_setinterpopts( nsurftypes, jtype, 'sic   ',      &
               &                  nn_2dint, nn_2dint_sic,             &
               &                  rn_sic_avglamscl, rn_sic_avgphiscl, &
               &                  ln_sic_fp_indegs, .FALSE.,          &
               &                  n2dintsurf, zavglamscl, zavgphiscl, &
               &                  lfpindegs, llnightav )
         ENDIF
#endif
         IF( ln_sss ) THEN
            jtype = jtype + 1
            CALL obs_settypefiles( nsurftypes, jpmaxnfiles, jtype, 'sss   ', &
               &                   cn_sssfbfiles, ifilessurf, cobstypessurf, clsurffiles )
            CALL obs_setinterpopts( nsurftypes, jtype, 'sss   ',      &
               &                  nn_2dint, nn_2dint_sss,             &
               &                  rn_sss_avglamscl, rn_sss_avgphiscl, &
               &                  ln_sss_fp_indegs, .FALSE.,          &
               &                  n2dintsurf, zavglamscl, zavgphiscl, &
               &                  lfpindegs, llnightav )
         ENDIF
         !
      ENDIF


      !-----------------------------------------------------------------------
      ! Obs operator parameter checking and initialisations
      !-----------------------------------------------------------------------
      !
      IF( ln_vel3d  .AND.  .NOT.ln_grid_global ) THEN
         CALL ctl_stop( 'Velocity data only works with ln_grid_global=.true.' )
         RETURN
      ENDIF
      !
      IF( ln_grid_global ) THEN
         CALL ctl_warn( 'dia_obs_init: ln_grid_global=T may cause memory issues when used with a large number of processors' )
      ENDIF
      !
      IF( nn_1dint < 0  .OR.  nn_1dint > 1 ) THEN
         CALL ctl_stop('dia_obs_init: Choice of vertical (1D) interpolation method is not available')
      ENDIF
      !
      IF( nn_2dint < 0  .OR.  nn_2dint > 6  ) THEN
         CALL ctl_stop('dia_obs_init: Choice of horizontal (2D) interpolation method is not available')
      ENDIF
      !
      CALL obs_typ_init
      IF( ln_grid_global )   CALL mppmap_init
      !
      CALL obs_grid_setup( )

      !-----------------------------------------------------------------------
      ! Depending on switches read the various observation types
      !-----------------------------------------------------------------------
      !
      IF( nproftypes > 0 ) THEN
         !
         ALLOCATE( profdata  (nproftypes) , nvarsprof (nproftypes) )
         ALLOCATE( profdataqc(nproftypes) , nextrprof (nproftypes) )
         !
         DO jtype = 1, nproftypes
            !
            nvarsprof(jtype) = 2
            IF ( TRIM(cobstypesprof(jtype)) == 'prof' ) THEN
               nextrprof(jtype) = 1
               llvar1 = ln_t3d
               llvar2 = ln_s3d
               zglam1 = glamt
               zgphi1 = gphit
               zmask1 = tmask
               zglam2 = glamt
               zgphi2 = gphit
               zmask2 = tmask
            ENDIF
            IF ( TRIM(cobstypesprof(jtype)) == 'vel' )  THEN
               nextrprof(jtype) = 2
               llvar1 = ln_vel3d
               llvar2 = ln_vel3d
               zglam1 = glamu
               zgphi1 = gphiu
               zmask1 = umask
               zglam2 = glamv
               zgphi2 = gphiv
               zmask2 = vmask
            ENDIF
            !
            ! Read in profile or profile obs types
            CALL obs_rea_prof( profdata(jtype), ifilesprof(jtype),       &
               &               clproffiles(jtype,1:ifilesprof(jtype)), &
               &               nvarsprof(jtype), nextrprof(jtype), nitend-nit000+2, &
               &               rn_dobsini, rn_dobsend, llvar1, llvar2, &
               &               ln_ignmis, ln_s_at_t, .FALSE., &
               &               kdailyavtypes = nn_profdavtypes )
               !
            DO jvar = 1, nvarsprof(jtype)
               CALL obs_prof_staend( profdata(jtype), jvar )
            END DO
            !
            CALL obs_pre_prof( profdata(jtype), profdataqc(jtype), &
               &               llvar1, llvar2, &
               &               jpi, jpj, jpk, &
               &               zmask1, zglam1, zgphi1, zmask2, zglam2, zgphi2,  &
               &               ln_nea, ln_bound_reject, &
               &               kdailyavtypes = nn_profdavtypes )
         END DO
         !
         DEALLOCATE( ifilesprof, clproffiles )
         !
      ENDIF
      !
      IF( nsurftypes > 0 ) THEN
         !
         ALLOCATE( surfdata  (nsurftypes) , nvarssurf(nsurftypes) )
         ALLOCATE( surfdataqc(nsurftypes) , nextrsurf(nsurftypes) )
         !
         DO jtype = 1, nsurftypes
            !
            nvarssurf(jtype) = 1
            nextrsurf(jtype) = 0
            llnightav(jtype) = .FALSE.
            IF( TRIM(cobstypessurf(jtype)) == 'sla' )   nextrsurf(jtype) = 2
            IF( TRIM(cobstypessurf(jtype)) == 'sst' )   llnightav(jtype) = ln_sstnight
            !
            ! Read in surface obs types
            CALL obs_rea_surf( surfdata(jtype), ifilessurf(jtype), &
               &               clsurffiles(jtype,1:ifilessurf(jtype)), &
               &               nvarssurf(jtype), nextrsurf(jtype), nitend-nit000+2, &
               &               rn_dobsini, rn_dobsend, ln_ignmis, .FALSE., llnightav(jtype) )
               !
            CALL obs_pre_surf( surfdata(jtype), surfdataqc(jtype), ln_nea, ln_bound_reject )
            !
            IF( TRIM(cobstypessurf(jtype)) == 'sla' ) THEN
               CALL obs_rea_mdt( surfdataqc(jtype), n2dintsurf(jtype) )
               IF( ln_altbias )   &
                  & CALL obs_rea_altbias ( surfdataqc(jtype), n2dintsurf(jtype), cn_altbiasfile )
            ENDIF
            !
            IF( TRIM(cobstypessurf(jtype)) == 'sst' .AND. ln_sstbias ) THEN
               jnumsstbias = 0
               DO jfile = 1, jpmaxnfiles
                  IF( TRIM(cn_sstbiasfiles(jfile)) /= '' )   jnumsstbias = jnumsstbias + 1
               END DO
               IF( jnumsstbias == 0 )   CALL ctl_stop("ln_sstbias set but no bias files to read in")    
               !
               CALL obs_app_sstbias( surfdataqc(jtype), n2dintsurf(jtype)             ,   & 
                  &                  jnumsstbias      , cn_sstbiasfiles(1:jnumsstbias) ) 
            ENDIF
         END DO
         !
         DEALLOCATE( ifilessurf, clsurffiles )
         !
      ENDIF
      !
   END SUBROUTINE dia_obs_init


   SUBROUTINE dia_obs( kstp )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_obs  ***
      !!          
      !! ** Purpose : Call the observation operators on each time step
      !!
      !! ** Method  : Call the observation operators on each time step to
      !!              compute the model equivalent of the following data:
      !!               - Profile data, currently T/S or U/V
      !!               - Surface data, currently SST, SLA or sea-ice concentration.
      !!
      !! ** Action  :
      !!----------------------------------------------------------------------
      USE dom_oce, ONLY : gdept_n, gdept_1d   ! Ocean space and time domain variables
      USE phycst , ONLY : rday                ! Physical constants
      USE oce    , ONLY : tsn, un, vn, sshn   ! Ocean dynamics and tracers variables
      USE phycst , ONLY : rday                ! Physical constants
#if defined  key_si3
      USE ice    , ONLY : at_i                ! SI3 Ice model variables
#endif
#if defined key_cice
      USE sbc_oce, ONLY : fr_i     ! ice fraction
#endif

      IMPLICIT NONE

      !! * Arguments
      INTEGER, INTENT(IN) :: kstp  ! Current timestep
      !! * Local declarations
      INTEGER :: idaystp           ! Number of timesteps per day
      INTEGER :: jtype             ! Data loop variable
      INTEGER :: jvar              ! Variable number
      INTEGER :: ji, jj            ! Loop counters
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: &
         & zprofvar1, &            ! Model values for 1st variable in a prof ob
         & zprofvar2               ! Model values for 2nd variable in a prof ob
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: &
         & zprofmask1, &           ! Mask associated with zprofvar1
         & zprofmask2              ! Mask associated with zprofvar2
      REAL(wp), DIMENSION(jpi,jpj) :: &
         & zsurfvar, &             ! Model values equivalent to surface ob.
         & zsurfmask               ! Mask associated with surface variable
      REAL(wp), DIMENSION(jpi,jpj) :: &
         & zglam1,    &            ! Model longitudes for prof variable 1
         & zglam2,    &            ! Model longitudes for prof variable 2
         & zgphi1,    &            ! Model latitudes for prof variable 1
         & zgphi2                  ! Model latitudes for prof variable 2

      !-----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_obs : Call the observation operators', kstp
         WRITE(numout,*) '~~~~~~~'
      ENDIF

      idaystp = NINT( rday / rdt )

      !-----------------------------------------------------------------------
      ! Call the profile and surface observation operators
      !-----------------------------------------------------------------------

      IF ( nproftypes > 0 ) THEN

         DO jtype = 1, nproftypes

            SELECT CASE ( TRIM(cobstypesprof(jtype)) )
            CASE('prof')
               zprofvar1(:,:,:) = tsn(:,:,:,jp_tem)
               zprofvar2(:,:,:) = tsn(:,:,:,jp_sal)
               zprofmask1(:,:,:) = tmask(:,:,:)
               zprofmask2(:,:,:) = tmask(:,:,:)
               zglam1(:,:) = glamt(:,:)
               zglam2(:,:) = glamt(:,:)
               zgphi1(:,:) = gphit(:,:)
               zgphi2(:,:) = gphit(:,:)
            CASE('vel')
               zprofvar1(:,:,:) = un(:,:,:)
               zprofvar2(:,:,:) = vn(:,:,:)
               zprofmask1(:,:,:) = umask(:,:,:)
               zprofmask2(:,:,:) = vmask(:,:,:)
               zglam1(:,:) = glamu(:,:)
               zglam2(:,:) = glamv(:,:)
               zgphi1(:,:) = gphiu(:,:)
               zgphi2(:,:) = gphiv(:,:)
            CASE DEFAULT
               CALL ctl_stop( 'Unknown profile observation type '//TRIM(cobstypesprof(jtype))//' in dia_obs' )
            END SELECT

            CALL obs_prof_opt( profdataqc(jtype), kstp, jpi, jpj, jpk,  &
               &               nit000, idaystp,                         &
               &               zprofvar1, zprofvar2,                    &
               &               gdept_n(:,:,:), gdepw_n(:,:,:),            & 
               &               zprofmask1, zprofmask2,                  &
               &               zglam1, zglam2, zgphi1, zgphi2,          &
               &               nn_1dint, nn_2dint,                      &
               &               kdailyavtypes = nn_profdavtypes )

         END DO

      ENDIF

      IF ( nsurftypes > 0 ) THEN

         DO jtype = 1, nsurftypes

            !Defaults which might be changed
            zsurfmask(:,:) = tmask(:,:,1)

            SELECT CASE ( TRIM(cobstypessurf(jtype)) )
            CASE('sst')
               zsurfvar(:,:) = tsn(:,:,1,jp_tem)
            CASE('sla')
               zsurfvar(:,:) = sshn(:,:)
            CASE('sss')
               zsurfvar(:,:) = tsn(:,:,1,jp_sal)
            CASE('sic')
               IF ( kstp == 0 ) THEN
                  IF ( lwp .AND. surfdataqc(jtype)%nsstpmpp(1) > 0 ) THEN
                     CALL ctl_warn( 'Sea-ice not initialised on zeroth '// &
                        &           'time-step but some obs are valid then.' )
                     WRITE(numout,*)surfdataqc(jtype)%nsstpmpp(1), &
                        &           ' sea-ice obs will be missed'
                  ENDIF
                  surfdataqc(jtype)%nsurfup = surfdataqc(jtype)%nsurfup + &
                     &                        surfdataqc(jtype)%nsstp(1)
                  CYCLE
               ELSE
#if defined key_cice || defined key_si3
                  zsurfvar(:,:) = fr_i(:,:)
#else
                  CALL ctl_stop( ' Trying to run sea-ice observation operator', &
                     &           ' but no sea-ice model appears to have been defined' )
#endif
               ENDIF

            END SELECT

            CALL obs_surf_opt( surfdataqc(jtype), kstp, jpi, jpj,       &
               &               nit000, idaystp, zsurfvar, zsurfmask,    &
               &               n2dintsurf(jtype), llnightav(jtype),     &
               &               zavglamscl(jtype), zavgphiscl(jtype),     &
               &               lfpindegs(jtype) )

         END DO

      ENDIF

   END SUBROUTINE dia_obs

   SUBROUTINE dia_obs_wri
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_obs_wri  ***
      !!          
      !! ** Purpose : Call observation diagnostic output routines
      !!
      !! ** Method  : Call observation diagnostic output routines
      !!
      !! ** Action  :
      !!
      !! History :
      !!        !  06-03  (K. Mogensen) Original code
      !!        !  06-05  (K. Mogensen) Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  07-03  (K. Mogensen) General handling of profiles
      !!        !  08-09  (M. Valdivieso) Velocity component (U,V) profiles
      !!        !  15-08  (M. Martin) Combined writing for prof and surf types
      !!----------------------------------------------------------------------
      !! * Modules used
      USE obs_rot_vel          ! Rotation of velocities

      IMPLICIT NONE

      !! * Local declarations
      INTEGER :: jtype                    ! Data set loop variable
      INTEGER :: jo, jvar, jk
      REAL(wp), DIMENSION(:), ALLOCATABLE :: &
         & zu, &
         & zv

      !-----------------------------------------------------------------------
      ! Depending on switches call various observation output routines
      !-----------------------------------------------------------------------

      IF ( nproftypes > 0 ) THEN

         DO jtype = 1, nproftypes

            IF ( TRIM(cobstypesprof(jtype)) == 'vel' ) THEN

               ! For velocity data, rotate the model velocities to N/S, E/W
               ! using the compressed data structure.
               ALLOCATE( &
                  & zu(profdataqc(jtype)%nvprot(1)), &
                  & zv(profdataqc(jtype)%nvprot(2))  &
                  & )

               CALL obs_rotvel( profdataqc(jtype), nn_2dint, zu, zv )

               DO jo = 1, profdataqc(jtype)%nprof
                  DO jvar = 1, 2
                     DO jk = profdataqc(jtype)%npvsta(jo,jvar), profdataqc(jtype)%npvend(jo,jvar)

                        IF ( jvar == 1 ) THEN
                           profdataqc(jtype)%var(jvar)%vmod(jk) = zu(jk)
                        ELSE
                           profdataqc(jtype)%var(jvar)%vmod(jk) = zv(jk)
                        ENDIF

                     END DO
                  END DO
               END DO

               DEALLOCATE( zu )
               DEALLOCATE( zv )

            END IF

            CALL obs_prof_decompress( profdataqc(jtype), &
               &                      profdata(jtype), .TRUE., numout )

            CALL obs_wri_prof( profdata(jtype) )

         END DO

      ENDIF

      IF ( nsurftypes > 0 ) THEN

         DO jtype = 1, nsurftypes

            CALL obs_surf_decompress( surfdataqc(jtype), &
               &                      surfdata(jtype), .TRUE., numout )

            CALL obs_wri_surf( surfdata(jtype) )

         END DO

      ENDIF

   END SUBROUTINE dia_obs_wri

   SUBROUTINE dia_obs_dealloc
      IMPLICIT NONE
      !!----------------------------------------------------------------------
      !!                    *** ROUTINE dia_obs_dealloc ***
      !!
      !!  ** Purpose : To deallocate data to enable the obs_oper online loop.
      !!               Specifically: dia_obs_init --> dia_obs --> dia_obs_wri
      !!
      !!  ** Method : Clean up various arrays left behind by the obs_oper.
      !!
      !!  ** Action :
      !!
      !!----------------------------------------------------------------------
      ! obs_grid deallocation
      CALL obs_grid_deallocate

      ! diaobs deallocation
      IF ( nproftypes > 0 ) &
         &   DEALLOCATE( cobstypesprof, profdata, profdataqc, nvarsprof, nextrprof )

      IF ( nsurftypes > 0 ) &
         &   DEALLOCATE( cobstypessurf, surfdata, surfdataqc, nvarssurf, nextrsurf, &
         &               n2dintsurf, zavglamscl, zavgphiscl, lfpindegs, llnightav )

   END SUBROUTINE dia_obs_dealloc

   SUBROUTINE calc_date( kstp, ddobs )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE calc_date  ***
      !!          
      !! ** Purpose : Get date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Method  : Get date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Action  : Get date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Action  : Get initial date in double precision YYYYMMDD.HHMMSS format
      !!
      !! History :
      !!        !  06-03  (K. Mogensen)  Original code
      !!        !  06-05  (K. Mogensen)  Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  06-10  (G. Smith) Calculates initial date the same as method for final date
      !!        !  10-05  (D. Lea) Update to month length calculation for NEMO vn3.2
      !!        !  2014-09  (D. Lea) New generic routine now deals with arbitrary initial time of day
      !!----------------------------------------------------------------------
      USE phycst, ONLY : &            ! Physical constants
         & rday
      USE dom_oce, ONLY : &           ! Ocean space and time domain variables
         & rdt

      IMPLICIT NONE

      !! * Arguments
      REAL(KIND=dp), INTENT(OUT) :: ddobs                        ! Date in YYYYMMDD.HHMMSS
      INTEGER :: kstp

      !! * Local declarations
      INTEGER :: iyea        ! date - (year, month, day, hour, minute)
      INTEGER :: imon
      INTEGER :: iday
      INTEGER :: ihou
      INTEGER :: imin
      INTEGER :: imday       ! Number of days in month.
      REAL(wp) :: zdayfrc    ! Fraction of day

      INTEGER, DIMENSION(12) ::   imonth_len    !: length in days of the months of the current year

      !!----------------------------------------------------------------------
      !! Initial date initialization (year, month, day, hour, minute)
      !!----------------------------------------------------------------------
      iyea =   ndate0 / 10000
      imon = ( ndate0 - iyea * 10000 ) / 100
      iday =   ndate0 - iyea * 10000 - imon * 100
      ihou =   nn_time0 / 100
      imin = ( nn_time0 - ihou * 100 ) 

      !!----------------------------------------------------------------------
      !! Compute number of days + number of hours + min since initial time
      !!----------------------------------------------------------------------
      zdayfrc = kstp * rdt / rday
      zdayfrc = zdayfrc - aint(zdayfrc)
      imin = imin + int( zdayfrc * 24 * 60 ) 
      DO WHILE (imin >= 60) 
        imin=imin-60
        ihou=ihou+1
      END DO
      DO WHILE (ihou >= 24)
        ihou=ihou-24
        iday=iday+1
      END DO 
      iday = iday + kstp * rdt / rday 

      !-----------------------------------------------------------------------
      ! Convert number of days (iday) into a real date
      !----------------------------------------------------------------------

      CALL calc_month_len( iyea, imonth_len )

      DO WHILE ( iday > imonth_len(imon) )
         iday = iday - imonth_len(imon)
         imon = imon + 1 
         IF ( imon > 12 ) THEN
            imon = 1
            iyea = iyea + 1
            CALL calc_month_len( iyea, imonth_len )  ! update month lengths
         ENDIF
      END DO

      !----------------------------------------------------------------------
      ! Convert it into YYYYMMDD.HHMMSS format.
      !----------------------------------------------------------------------
      ddobs = iyea * 10000_dp + imon * 100_dp + &
         &    iday + ihou * 0.01_dp + imin * 0.0001_dp

   END SUBROUTINE calc_date

   SUBROUTINE ini_date( ddobsini )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ini_date  ***
      !!          
      !! ** Purpose : Get initial date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  06-03  (K. Mogensen)  Original code
      !!        !  06-05  (K. Mogensen)  Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  10-05  (D. Lea) Update to month length calculation for NEMO vn3.2
      !!        !  2014-09  (D. Lea) Change to call generic routine calc_date
      !!----------------------------------------------------------------------

      IMPLICIT NONE

      !! * Arguments
      REAL(KIND=dp), INTENT(OUT) :: ddobsini                   ! Initial date in YYYYMMDD.HHMMSS

      CALL calc_date( nit000 - 1, ddobsini )

   END SUBROUTINE ini_date

   SUBROUTINE fin_date( ddobsfin )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE fin_date  ***
      !!          
      !! ** Purpose : Get final date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  06-03  (K. Mogensen)  Original code
      !!        !  06-05  (K. Mogensen)  Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  10-05  (D. Lea) Update to month length calculation for NEMO vn3.2
      !!        !  2014-09  (D. Lea) Change to call generic routine calc_date
      !!----------------------------------------------------------------------

      IMPLICIT NONE

      !! * Arguments
      REAL(dp), INTENT(OUT) :: ddobsfin ! Final date in YYYYMMDD.HHMMSS

      CALL calc_date( nitend, ddobsfin )

   END SUBROUTINE fin_date
   
    SUBROUTINE obs_settypefiles( ntypes, jpmaxnfiles, jtype, ctypein, &
       &                         cfilestype, ifiles, cobstypes, cfiles )

    INTEGER, INTENT(IN) :: ntypes      ! Total number of obs types
    INTEGER, INTENT(IN) :: jpmaxnfiles ! Maximum number of files allowed for each type
    INTEGER, INTENT(IN) :: jtype       ! Index of the current type of obs
    INTEGER, DIMENSION(ntypes), INTENT(INOUT) :: &
       &                   ifiles      ! Out appended number of files for this type

    CHARACTER(len=6), INTENT(IN) :: ctypein 
    CHARACTER(len=128), DIMENSION(jpmaxnfiles), INTENT(IN) :: &
       &                   cfilestype  ! In list of files for this obs type
    CHARACTER(len=6), DIMENSION(ntypes), INTENT(INOUT) :: &
       &                   cobstypes   ! Out appended list of obs types
    CHARACTER(len=128), DIMENSION(ntypes, jpmaxnfiles), INTENT(INOUT) :: &
       &                   cfiles      ! Out appended list of files for all types

    !Local variables
    INTEGER :: jfile

    cfiles(jtype,:) = cfilestype(:)
    cobstypes(jtype) = ctypein
    ifiles(jtype) = 0
    DO jfile = 1, jpmaxnfiles
       IF ( trim(cfiles(jtype,jfile)) /= '' ) &
                 ifiles(jtype) = ifiles(jtype) + 1
    END DO

    IF ( ifiles(jtype) == 0 ) THEN
         CALL ctl_stop( 'Logical for observation type '//TRIM(ctypein)//   &
            &           ' set to true but no files available to read' )
    ENDIF

    IF(lwp) THEN    
       WRITE(numout,*) '             '//cobstypes(jtype)//' input observation file names:'
       DO jfile = 1, ifiles(jtype)
          WRITE(numout,*) '                '//TRIM(cfiles(jtype,jfile))
       END DO
    ENDIF

    END SUBROUTINE obs_settypefiles

    SUBROUTINE obs_setinterpopts( ntypes, jtype, ctypein,             &
               &                  n2dint_default, n2dint_type,        &
               &                  zavglamscl_type, zavgphiscl_type,   &
               &                  lfp_indegs_type, lavnight_type,     &
               &                  n2dint, zavglamscl, zavgphiscl,     &
               &                  lfpindegs, lavnight )

    INTEGER, INTENT(IN)  :: ntypes             ! Total number of obs types
    INTEGER, INTENT(IN)  :: jtype              ! Index of the current type of obs
    INTEGER, INTENT(IN)  :: n2dint_default     ! Default option for interpolation type
    INTEGER, INTENT(IN)  :: n2dint_type        ! Option for interpolation type
    REAL(wp), INTENT(IN) :: &
       &                    zavglamscl_type, & !E/W diameter of obs footprint for this type
       &                    zavgphiscl_type    !N/S diameter of obs footprint for this type
    LOGICAL, INTENT(IN)  :: lfp_indegs_type    !T=> footprint in degrees, F=> in metres
    LOGICAL, INTENT(IN)  :: lavnight_type      !T=> obs represent night time average
    CHARACTER(len=6), INTENT(IN) :: ctypein 

    INTEGER, DIMENSION(ntypes), INTENT(INOUT) :: &
       &                    n2dint 
    REAL(wp), DIMENSION(ntypes), INTENT(INOUT) :: &
       &                    zavglamscl, zavgphiscl
    LOGICAL, DIMENSION(ntypes), INTENT(INOUT) :: &
       &                    lfpindegs, lavnight

    lavnight(jtype) = lavnight_type

    IF ( (n2dint_type >= 1) .AND. (n2dint_type <= 6) ) THEN
       n2dint(jtype) = n2dint_type
    ELSE
       n2dint(jtype) = n2dint_default
    ENDIF

    ! For averaging observation footprints set options for size of footprint 
    IF ( (n2dint(jtype) > 4) .AND. (n2dint(jtype) <= 6) ) THEN
       IF ( zavglamscl_type > 0._wp ) THEN
          zavglamscl(jtype) = zavglamscl_type
       ELSE
          CALL ctl_stop( 'Incorrect value set for averaging footprint '// &
                         'scale (zavglamscl) for observation type '//TRIM(ctypein) )      
       ENDIF

       IF ( zavgphiscl_type > 0._wp ) THEN
          zavgphiscl(jtype) = zavgphiscl_type
       ELSE
          CALL ctl_stop( 'Incorrect value set for averaging footprint '// &
                         'scale (zavgphiscl) for observation type '//TRIM(ctypein) )      
       ENDIF

       lfpindegs(jtype) = lfp_indegs_type 

    ENDIF

    ! Write out info 
    IF(lwp) THEN
       IF ( n2dint(jtype) <= 4 ) THEN
          WRITE(numout,*) '             '//TRIM(ctypein)// &
             &            ' model counterparts will be interpolated horizontally'
       ELSE IF ( n2dint(jtype) <= 6 ) THEN
          WRITE(numout,*) '             '//TRIM(ctypein)// &
             &            ' model counterparts will be averaged horizontally'
          WRITE(numout,*) '             '//'    with E/W scale: ',zavglamscl(jtype)
          WRITE(numout,*) '             '//'    with N/S scale: ',zavgphiscl(jtype)
          IF ( lfpindegs(jtype) ) THEN
              WRITE(numout,*) '             '//'    (in degrees)'
          ELSE
              WRITE(numout,*) '             '//'    (in metres)'
          ENDIF
       ENDIF
    ENDIF

    END SUBROUTINE obs_setinterpopts

END MODULE diaobs
