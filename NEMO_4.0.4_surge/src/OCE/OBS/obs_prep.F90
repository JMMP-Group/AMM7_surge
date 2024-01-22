MODULE obs_prep
   !!=====================================================================
   !!                       ***  MODULE  obs_prep  ***
   !! Observation diagnostics: Prepare observation arrays: screening, 
   !!                          sorting, coordinate search
   !!=====================================================================

   !!---------------------------------------------------------------------
   !!   obs_pre_prof  : First level check and screening of profile observations
   !!   obs_pre_surf  : First level check and screening of surface observations
   !!   obs_scr       : Basic screening of the observations
   !!   obs_coo_tim   : Compute number of time steps to the observation time
   !!   obs_sor       : Sort the observation arrays
   !!---------------------------------------------------------------------
   USE par_kind, ONLY : wp ! Precision variables
   USE in_out_manager     ! I/O manager
   USE obs_profiles_def   ! Definitions for storage arrays for profiles
   USE obs_surf_def       ! Definitions for storage arrays for surface data
   USE obs_mpp, ONLY : &  ! MPP support routines for observation diagnostics
      & obs_mpp_sum_integer, &
      & obs_mpp_sum_integers
   USE obs_inter_sup      ! Interpolation support
   USE obs_oper           ! Observation operators
   USE lib_mpp, ONLY :   ctl_warn, ctl_stop
   USE bdy_oce, ONLY : &        ! Boundary information
      idx_bdy, nb_bdy, ln_bdy

   IMPLICIT NONE
   PRIVATE

   PUBLIC   obs_pre_prof     ! First level check and screening of profile obs
   PUBLIC   obs_pre_surf     ! First level check and screening of surface obs
   PUBLIC   calc_month_len   ! Calculate the number of days in the months of a year

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------


CONTAINS

   SUBROUTINE obs_pre_surf( surfdata, surfdataqc, ld_nea, ld_bound_reject, &
                            kqc_cutoff )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE obs_pre_sla  ***
      !!
      !! ** Purpose : First level check and screening of surface observations
      !!
      !! ** Method  : First level check and screening of surface observations
      !!
      !! ** Action  : 
      !!
      !! References :
      !!   
      !! History :
      !!        !  2007-03  (A. Weaver, K. Mogensen) Original
      !!        !  2007-06  (K. Mogensen et al) Reject obs. near land.
      !!        !  2015-02  (M. Martin) Combined routine for surface types.
      !!----------------------------------------------------------------------
      !! * Modules used
      USE par_oce             ! Ocean parameters
      USE dom_oce, ONLY       :   glamt, gphit, tmask, nproc   ! Geographical information
      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) :: surfdata    ! Full set of surface data
      TYPE(obs_surf), INTENT(INOUT) :: surfdataqc   ! Subset of surface data not failing screening
      LOGICAL, INTENT(IN) :: ld_nea         ! Switch for rejecting observation near land
      LOGICAL, INTENT(IN) :: ld_bound_reject       ! Switch for rejecting obs near the boundary
      INTEGER, INTENT(IN), OPTIONAL :: kqc_cutoff   ! cut off for QC value
      !! * Local declarations
      INTEGER :: iqc_cutoff = 255   ! cut off for QC value
      INTEGER :: iyea0        ! Initial date
      INTEGER :: imon0        !  - (year, month, day, hour, minute)
      INTEGER :: iday0    
      INTEGER :: ihou0    
      INTEGER :: imin0
      INTEGER :: icycle       ! Current assimilation cycle
                              ! Counters for observations that
      INTEGER :: iotdobs      !  - outside time domain
      INTEGER :: iosdsobs     !  - outside space domain
      INTEGER :: ilansobs     !  - within a model land cell
      INTEGER :: inlasobs     !  - close to land
      INTEGER :: igrdobs      !  - fail the grid search
      INTEGER :: ibdysobs     !  - close to open boundary
                              ! Global counters for observations that
      INTEGER :: iotdobsmpp     !  - outside time domain
      INTEGER :: iosdsobsmpp    !  - outside space domain
      INTEGER :: ilansobsmpp    !  - within a model land cell
      INTEGER :: inlasobsmpp    !  - close to land
      INTEGER :: igrdobsmpp     !  - fail the grid search
      INTEGER :: ibdysobsmpp  !  - close to open boundary
      LOGICAL, DIMENSION(:), ALLOCATABLE :: & 
         & llvalid            ! SLA data selection
      INTEGER :: jobs         ! Obs. loop variable
      INTEGER :: jstp         ! Time loop variable
      INTEGER :: inrc         ! Time index variable
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*) 'obs_pre_surf : Preparing the surface observations...'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      
      ! Initial date initialization (year, month, day, hour, minute)
      iyea0 =   ndate0 / 10000
      imon0 = ( ndate0 - iyea0 * 10000 ) / 100
      iday0 =   ndate0 - iyea0 * 10000 - imon0 * 100
      ihou0 = nn_time0 / 100
      imin0 = ( nn_time0 - ihou0 * 100 )

      icycle = nn_no     ! Assimilation cycle

      ! Diagnotics counters for various failures.

      iotdobs  = 0
      igrdobs  = 0
      iosdsobs = 0
      ilansobs = 0
      inlasobs = 0
      ibdysobs = 0 

      ! Set QC cutoff to optional value if provided
      IF ( PRESENT(kqc_cutoff) ) iqc_cutoff=kqc_cutoff

      ! -----------------------------------------------------------------------
      ! Find time coordinate for surface data
      ! -----------------------------------------------------------------------

      CALL obs_coo_tim( icycle, &
         &              iyea0,   imon0,   iday0,   ihou0,   imin0,      &
         &              surfdata%nsurf,   surfdata%nyea, surfdata%nmon, &
         &              surfdata%nday,    surfdata%nhou, surfdata%nmin, &
         &              surfdata%nqc,     surfdata%mstp, iotdobs        )

      CALL obs_mpp_sum_integer( iotdobs, iotdobsmpp )
      
      ! -----------------------------------------------------------------------
      ! Check for surface data failing the grid search
      ! -----------------------------------------------------------------------

      CALL obs_coo_grd( surfdata%nsurf,   surfdata%mi, surfdata%mj, &
         &              surfdata%nqc,     igrdobs                         )

      CALL obs_mpp_sum_integer( igrdobs, igrdobsmpp )

      ! -----------------------------------------------------------------------
      ! Check for land points. 
      ! -----------------------------------------------------------------------

      CALL obs_coo_spc_2d( surfdata%nsurf,              &
         &                 jpi,          jpj,          &
         &                 surfdata%mi,   surfdata%mj,   & 
         &                 surfdata%rlam, surfdata%rphi, &
         &                 glamt,        gphit,        &
         &                 tmask(:,:,1), surfdata%nqc,  &
         &                 iosdsobs,     ilansobs,     &
         &                 inlasobs,     ld_nea,       &
         &                 ibdysobs,     ld_bound_reject, &
         &                 iqc_cutoff                     )

      CALL obs_mpp_sum_integer( iosdsobs, iosdsobsmpp )
      CALL obs_mpp_sum_integer( ilansobs, ilansobsmpp )
      CALL obs_mpp_sum_integer( inlasobs, inlasobsmpp )
      CALL obs_mpp_sum_integer( ibdysobs, ibdysobsmpp )

      ! -----------------------------------------------------------------------
      ! Copy useful data from the surfdata data structure to
      ! the surfdataqc data structure 
      ! -----------------------------------------------------------------------

      ! Allocate the selection arrays

      ALLOCATE( llvalid(surfdata%nsurf) )
      
      ! We want all data which has qc flags <= iqc_cutoff

      llvalid(:)  = ( surfdata%nqc(:)  <= iqc_cutoff )

      ! The actual copying

      CALL obs_surf_compress( surfdata,     surfdataqc,       .TRUE.,  numout, &
         &                    lvalid=llvalid )

      ! Dellocate the selection arrays
      DEALLOCATE( llvalid )

      ! -----------------------------------------------------------------------
      ! Print information about what observations are left after qc
      ! -----------------------------------------------------------------------

      ! Update the total observation counter array
      
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' '//surfdataqc%cvars(1)//' data outside time domain                  = ', &
            &            iotdobsmpp
         WRITE(numout,*) ' Remaining '//surfdataqc%cvars(1)//' data that failed grid search    = ', &
            &            igrdobsmpp
         WRITE(numout,*) ' Remaining '//surfdataqc%cvars(1)//' data outside space domain       = ', &
            &            iosdsobsmpp
         WRITE(numout,*) ' Remaining '//surfdataqc%cvars(1)//' data at land points             = ', &
            &            ilansobsmpp
         IF (ld_nea) THEN
            WRITE(numout,*) ' Remaining '//surfdataqc%cvars(1)//' data near land points (removed) = ', &
               &            inlasobsmpp
         ELSE
            WRITE(numout,*) ' Remaining '//surfdataqc%cvars(1)//' data near land points (kept)    = ', &
               &            inlasobsmpp
         ENDIF
         WRITE(numout,*) ' Remaining '//surfdataqc%cvars(1)//' data near open boundary (removed) = ', &
            &            ibdysobsmpp  
         WRITE(numout,*) ' '//surfdataqc%cvars(1)//' data accepted                             = ', &
            &            surfdataqc%nsurfmpp

         WRITE(numout,*)
         WRITE(numout,*) ' Number of observations per time step :'
         WRITE(numout,*)
         WRITE(numout,'(10X,A,10X,A)')'Time step',surfdataqc%cvars(1)
         WRITE(numout,'(10X,A,5X,A)')'---------','-----------------'
         CALL FLUSH(numout)
      ENDIF
      
      DO jobs = 1, surfdataqc%nsurf
         inrc = surfdataqc%mstp(jobs) + 2 - nit000
         surfdataqc%nsstp(inrc)  = surfdataqc%nsstp(inrc) + 1
      END DO
      
      CALL obs_mpp_sum_integers( surfdataqc%nsstp, surfdataqc%nsstpmpp, &
         &                       nitend - nit000 + 2 )

      IF ( lwp ) THEN
         DO jstp = nit000 - 1, nitend
            inrc = jstp - nit000 + 2
            WRITE(numout,1999) jstp, surfdataqc%nsstpmpp(inrc)
            CALL FLUSH(numout)
         END DO
      ENDIF

1999  FORMAT(10X,I9,5X,I17)

   END SUBROUTINE obs_pre_surf


   SUBROUTINE obs_pre_prof( profdata, prodatqc, ld_var1, ld_var2, &
      &                     kpi, kpj, kpk, &
      &                     zmask1, pglam1, pgphi1, zmask2, pglam2, pgphi2,  &
      &                     ld_nea, ld_bound_reject, kdailyavtypes,  kqc_cutoff )

!!----------------------------------------------------------------------
      !!                    ***  ROUTINE obs_pre_prof  ***
      !!
      !! ** Purpose : First level check and screening of profiles
      !!
      !! ** Method  : First level check and screening of profiles
      !!
      !! History :
      !!        !  2007-06  (K. Mogensen) original : T and S profile data
      !!        !  2008-09  (M. Valdivieso) : TAO velocity data
      !!        !  2009-01  (K. Mogensen) : New feedback stricture
      !!        !  2015-02  (M. Martin) : Combined profile routine.
      !!
      !!----------------------------------------------------------------------
      !! * Modules used
      USE par_oce             ! Ocean parameters
      USE dom_oce, ONLY : &   ! Geographical information
         & gdept_1d,             &
         & nproc

      !! * Arguments
      TYPE(obs_prof), INTENT(INOUT) :: profdata   ! Full set of profile data
      TYPE(obs_prof), INTENT(INOUT) :: prodatqc   ! Subset of profile data not failing screening
      LOGICAL, INTENT(IN) :: ld_var1              ! Observed variables switches
      LOGICAL, INTENT(IN) :: ld_var2
      LOGICAL, INTENT(IN) :: ld_nea               ! Switch for rejecting observation near land
      LOGICAL, INTENT(IN) :: ld_bound_reject      ! Switch for rejecting observations near the boundary
      INTEGER, INTENT(IN) :: kpi, kpj, kpk        ! Local domain sizes
      INTEGER, DIMENSION(imaxavtypes), OPTIONAL :: &
         & kdailyavtypes                          ! Types for daily averages
      REAL(wp), INTENT(IN), DIMENSION(kpi,kpj,kpk) :: &
         & zmask1, &
         & zmask2
      REAL(wp), INTENT(IN), DIMENSION(kpi,kpj) :: &
         & pglam1, &
         & pglam2, &
         & pgphi1, &
         & pgphi2
      INTEGER, INTENT(IN), OPTIONAL :: kqc_cutoff   ! cut off for QC value

      !! * Local declarations
      INTEGER :: iqc_cutoff = 255   ! cut off for QC value
      INTEGER :: iyea0        ! Initial date
      INTEGER :: imon0        !  - (year, month, day, hour, minute)
      INTEGER :: iday0    
      INTEGER :: ihou0    
      INTEGER :: imin0
      INTEGER :: icycle       ! Current assimilation cycle
                              ! Counters for observations that are
      INTEGER :: iotdobs      !  - outside time domain
      INTEGER :: iosdv1obs    !  - outside space domain (variable 1)
      INTEGER :: iosdv2obs    !  - outside space domain (variable 2)
      INTEGER :: ilanv1obs    !  - within a model land cell (variable 1)
      INTEGER :: ilanv2obs    !  - within a model land cell (variable 2)
      INTEGER :: inlav1obs    !  - close to land (variable 1)
      INTEGER :: inlav2obs    !  - close to land (variable 2)
      INTEGER :: ibdyv1obs    !  - boundary (variable 1) 
      INTEGER :: ibdyv2obs    !  - boundary (variable 2)      
      INTEGER :: igrdobs      !  - fail the grid search
      INTEGER :: iuvchku      !  - reject u if v rejected and vice versa
      INTEGER :: iuvchkv      !
                              ! Global counters for observations that are
      INTEGER :: iotdobsmpp   !  - outside time domain
      INTEGER :: iosdv1obsmpp !  - outside space domain (variable 1)
      INTEGER :: iosdv2obsmpp !  - outside space domain (variable 2)
      INTEGER :: ilanv1obsmpp !  - within a model land cell (variable 1)
      INTEGER :: ilanv2obsmpp !  - within a model land cell (variable 2)
      INTEGER :: inlav1obsmpp !  - close to land (variable 1)
      INTEGER :: inlav2obsmpp !  - close to land (variable 2)
      INTEGER :: ibdyv1obsmpp !  - boundary (variable 1) 
      INTEGER :: ibdyv2obsmpp !  - boundary (variable 2)      
      INTEGER :: igrdobsmpp   !  - fail the grid search
      INTEGER :: iuvchkumpp   !  - reject var1 if var2 rejected and vice versa
      INTEGER :: iuvchkvmpp   !
      TYPE(obs_prof_valid) ::  llvalid      ! Profile selection 
      TYPE(obs_prof_valid), DIMENSION(profdata%nvar) :: &
         & llvvalid           ! var1,var2 selection 
      INTEGER :: jvar         ! Variable loop variable
      INTEGER :: jobs         ! Obs. loop variable
      INTEGER :: jstp         ! Time loop variable
      INTEGER :: inrc         ! Time index variable
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)'obs_pre_prof: Preparing the profile data...'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'

      ! Initial date initialization (year, month, day, hour, minute)
      iyea0 =   ndate0 / 10000
      imon0 = ( ndate0 - iyea0 * 10000 ) / 100
      iday0 =   ndate0 - iyea0 * 10000 - imon0 * 100
      ihou0 = nn_time0 / 100
      imin0 = ( nn_time0 - ihou0 * 100 )

      icycle = nn_no     ! Assimilation cycle

      ! Diagnotics counters for various failures.

      iotdobs   = 0
      igrdobs   = 0
      iosdv1obs = 0
      iosdv2obs = 0
      ilanv1obs = 0
      ilanv2obs = 0
      inlav1obs = 0
      inlav2obs = 0
      ibdyv1obs = 0
      ibdyv2obs = 0
      iuvchku   = 0
      iuvchkv   = 0


      ! Set QC cutoff to optional value if provided
      IF ( PRESENT(kqc_cutoff) ) iqc_cutoff=kqc_cutoff

      ! -----------------------------------------------------------------------
      ! Find time coordinate for profiles
      ! -----------------------------------------------------------------------

      IF ( PRESENT(kdailyavtypes) ) THEN
         CALL obs_coo_tim_prof( icycle, &
            &              iyea0,   imon0,   iday0,   ihou0,   imin0,      &
            &              profdata%nprof,   profdata%nyea, profdata%nmon, &
            &              profdata%nday,    profdata%nhou, profdata%nmin, &
            &              profdata%ntyp,    profdata%nqc,  profdata%mstp, &
            &              iotdobs, kdailyavtypes = kdailyavtypes,         &
            &              kqc_cutoff = iqc_cutoff )
      ELSE
         CALL obs_coo_tim_prof( icycle, &
            &              iyea0,   imon0,   iday0,   ihou0,   imin0,      &
            &              profdata%nprof,   profdata%nyea, profdata%nmon, &
            &              profdata%nday,    profdata%nhou, profdata%nmin, &
            &              profdata%ntyp,    profdata%nqc,  profdata%mstp, &
            &              iotdobs,          kqc_cutoff = iqc_cutoff )
      ENDIF

      CALL obs_mpp_sum_integer( iotdobs, iotdobsmpp )
      
      ! -----------------------------------------------------------------------
      ! Check for profiles failing the grid search
      ! -----------------------------------------------------------------------

      CALL obs_coo_grd( profdata%nprof,   profdata%mi(:,1), profdata%mj(:,1), &
         &              profdata%nqc,     igrdobs                         )
      CALL obs_coo_grd( profdata%nprof,   profdata%mi(:,2), profdata%mj(:,2), &
         &              profdata%nqc,     igrdobs                         )

      CALL obs_mpp_sum_integer( igrdobs, igrdobsmpp )

      ! -----------------------------------------------------------------------
      ! Reject all observations for profiles with nqc > iqc_cutoff
      ! -----------------------------------------------------------------------

      CALL obs_pro_rej( profdata, kqc_cutoff = iqc_cutoff )

      ! -----------------------------------------------------------------------
      ! Check for land points. This includes points below the model
      ! bathymetry so this is done for every point in the profile
      ! -----------------------------------------------------------------------

      ! Variable 1
      CALL obs_coo_spc_3d( profdata%nprof,        profdata%nvprot(1),   &
         &                 profdata%npvsta(:,1),  profdata%npvend(:,1), &
         &                 jpi,                   jpj,                  &
         &                 jpk,                                         &
         &                 profdata%mi,           profdata%mj,          &
         &                 profdata%var(1)%mvk,                         &
         &                 profdata%rlam,         profdata%rphi,        &
         &                 profdata%var(1)%vdep,                        &
         &                 pglam1,                pgphi1,               &
         &                 gdept_1d,              zmask1,               &
         &                 profdata%nqc,          profdata%var(1)%nvqc, &
         &                 iosdv1obs,             ilanv1obs,            &
         &                 inlav1obs,             ld_nea,               &
         &                 ibdyv1obs,             ld_bound_reject,      &
         &                 iqc_cutoff       )

      CALL obs_mpp_sum_integer( iosdv1obs, iosdv1obsmpp )
      CALL obs_mpp_sum_integer( ilanv1obs, ilanv1obsmpp )
      CALL obs_mpp_sum_integer( inlav1obs, inlav1obsmpp )
      CALL obs_mpp_sum_integer( ibdyv1obs, ibdyv1obsmpp )

      ! Variable 2
      CALL obs_coo_spc_3d( profdata%nprof,        profdata%nvprot(2),   &
         &                 profdata%npvsta(:,2),  profdata%npvend(:,2), &
         &                 jpi,                   jpj,                  &
         &                 jpk,                                         &
         &                 profdata%mi,           profdata%mj,          & 
         &                 profdata%var(2)%mvk,                         &
         &                 profdata%rlam,         profdata%rphi,        &
         &                 profdata%var(2)%vdep,                        &
         &                 pglam2,                pgphi2,               &
         &                 gdept_1d,              zmask2,               &
         &                 profdata%nqc,          profdata%var(2)%nvqc, &
         &                 iosdv2obs,             ilanv2obs,            &
         &                 inlav2obs,             ld_nea,               &
         &                 ibdyv2obs,             ld_bound_reject,      &
         &                 iqc_cutoff       )

      CALL obs_mpp_sum_integer( iosdv2obs, iosdv2obsmpp )
      CALL obs_mpp_sum_integer( ilanv2obs, ilanv2obsmpp )
      CALL obs_mpp_sum_integer( inlav2obs, inlav2obsmpp )
      CALL obs_mpp_sum_integer( ibdyv2obs, ibdyv2obsmpp )

      ! -----------------------------------------------------------------------
      ! Reject u if v is rejected and vice versa
      ! -----------------------------------------------------------------------

      IF ( TRIM(profdata%cvars(1)) == 'UVEL' ) THEN
         CALL obs_uv_rej( profdata, iuvchku, iuvchkv, iqc_cutoff )
         CALL obs_mpp_sum_integer( iuvchku, iuvchkumpp )
         CALL obs_mpp_sum_integer( iuvchkv, iuvchkvmpp )
      ENDIF

      ! -----------------------------------------------------------------------
      ! Copy useful data from the profdata data structure to
      ! the prodatqc data structure 
      ! -----------------------------------------------------------------------

      ! Allocate the selection arrays

      ALLOCATE( llvalid%luse(profdata%nprof) )
      DO jvar = 1,profdata%nvar
         ALLOCATE( llvvalid(jvar)%luse(profdata%nvprot(jvar)) )
      END DO

      ! We want all data which has qc flags <= iqc_cutoff

      llvalid%luse(:) = ( profdata%nqc(:)  <= iqc_cutoff )
      DO jvar = 1,profdata%nvar
         llvvalid(jvar)%luse(:) = ( profdata%var(jvar)%nvqc(:) <= iqc_cutoff )
      END DO

      ! The actual copying

      CALL obs_prof_compress( profdata,     prodatqc,       .TRUE.,  numout, &
         &                    lvalid=llvalid, lvvalid=llvvalid )

      ! Dellocate the selection arrays
      DEALLOCATE( llvalid%luse )
      DO jvar = 1,profdata%nvar
         DEALLOCATE( llvvalid(jvar)%luse )
      END DO

      ! -----------------------------------------------------------------------
      ! Print information about what observations are left after qc
      ! -----------------------------------------------------------------------

      ! Update the total observation counter array
      
      IF(lwp) THEN
      
         WRITE(numout,*)
         WRITE(numout,*) ' Profiles outside time domain                     = ', &
            &            iotdobsmpp
         WRITE(numout,*) ' Remaining profiles that failed grid search       = ', &
            &            igrdobsmpp
         WRITE(numout,*) ' Remaining '//prodatqc%cvars(1)//' data outside space domain       = ', &
            &            iosdv1obsmpp
         WRITE(numout,*) ' Remaining '//prodatqc%cvars(1)//' data at land points             = ', &
            &            ilanv1obsmpp
         IF (ld_nea) THEN
            WRITE(numout,*) ' Remaining '//prodatqc%cvars(1)//' data near land points (removed) = ',&
               &            inlav1obsmpp
         ELSE
            WRITE(numout,*) ' Remaining '//prodatqc%cvars(1)//' data near land points (kept)    = ',&
               &            inlav1obsmpp
         ENDIF
         IF ( TRIM(profdata%cvars(1)) == 'UVEL' ) THEN
            WRITE(numout,*) ' U observation rejected since V rejected     = ', &
               &            iuvchku
         ENDIF
         WRITE(numout,*) ' Remaining '//prodatqc%cvars(1)//' data near open boundary (removed) = ',&
               &            ibdyv1obsmpp
         WRITE(numout,*) ' '//prodatqc%cvars(1)//' data accepted                             = ', &
            &            prodatqc%nvprotmpp(1)
         WRITE(numout,*) ' Remaining '//prodatqc%cvars(2)//' data outside space domain       = ', &
            &            iosdv2obsmpp
         WRITE(numout,*) ' Remaining '//prodatqc%cvars(2)//' data at land points             = ', &
            &            ilanv2obsmpp
         IF (ld_nea) THEN
            WRITE(numout,*) ' Remaining '//prodatqc%cvars(2)//' data near land points (removed) = ',&
               &            inlav2obsmpp
         ELSE
            WRITE(numout,*) ' Remaining '//prodatqc%cvars(2)//' data near land points (kept)    = ',&
               &            inlav2obsmpp
         ENDIF
         IF ( TRIM(profdata%cvars(1)) == 'UVEL' ) THEN
            WRITE(numout,*) ' V observation rejected since U rejected     = ', &
               &            iuvchkv
         ENDIF
         WRITE(numout,*) ' Remaining '//prodatqc%cvars(2)//' data near open boundary (removed) = ',&
               &            ibdyv2obsmpp
         WRITE(numout,*) ' '//prodatqc%cvars(2)//' data accepted                             = ', &
            &            prodatqc%nvprotmpp(2)

         WRITE(numout,*)
         WRITE(numout,*) ' Number of observations per time step :'
         WRITE(numout,*)
         WRITE(numout,'(10X,A,5X,A,5X,A,A)')'Time step','Profiles', &
            &                               '     '//prodatqc%cvars(1)//'     ', &
            &                               '     '//prodatqc%cvars(2)//'     '
         WRITE(numout,998)
      ENDIF
      
      DO jobs = 1, prodatqc%nprof
         inrc = prodatqc%mstp(jobs) + 2 - nit000
         prodatqc%npstp(inrc)  = prodatqc%npstp(inrc) + 1
         DO jvar = 1, prodatqc%nvar
            IF ( prodatqc%npvend(jobs,jvar) > 0 ) THEN
               prodatqc%nvstp(inrc,jvar) = prodatqc%nvstp(inrc,jvar) + &
                  &                      ( prodatqc%npvend(jobs,jvar) - &
                  &                        prodatqc%npvsta(jobs,jvar) + 1 )
            ENDIF
         END DO
      END DO
      
      
      CALL obs_mpp_sum_integers( prodatqc%npstp, prodatqc%npstpmpp, &
         &                       nitend - nit000 + 2 )
      DO jvar = 1, prodatqc%nvar
         CALL obs_mpp_sum_integers( prodatqc%nvstp(:,jvar), &
            &                       prodatqc%nvstpmpp(:,jvar), &
            &                       nitend - nit000 + 2 )
      END DO

      IF ( lwp ) THEN
         DO jstp = nit000 - 1, nitend
            inrc = jstp - nit000 + 2
            WRITE(numout,999) jstp, prodatqc%npstpmpp(inrc), &
               &                    prodatqc%nvstpmpp(inrc,1), &
               &                    prodatqc%nvstpmpp(inrc,2)
         END DO
      ENDIF

998   FORMAT(10X,'---------',5X,'--------',5X,'-----------',5X,'----------------')
999   FORMAT(10X,I9,5X,I8,5X,I11,5X,I8)

   END SUBROUTINE obs_pre_prof

   SUBROUTINE obs_coo_tim( kcycle, &
      &                    kyea0,   kmon0,   kday0,   khou0,   kmin0,     &
      &                    kobsno,                                        &
      &                    kobsyea, kobsmon, kobsday, kobshou, kobsmin,   &
      &                    kobsqc,  kobsstp, kotdobs                      )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE obs_coo_tim ***
      !!
      !! ** Purpose : Compute the number of time steps to the observation time.
      !!
      !! ** Method  : For time coordinates ( yea_obs, mon_obs, day_obs, 
      !!              hou_obs, min_obs ), this routine locates the time step 
      !!              that is closest to this time.
      !!
      !! ** Action  : 
      !!
      !! References :
      !!   
      !! History :
      !!        !  1997-07  (A. Weaver) Original
      !!        !  2006-08  (A. Weaver) NEMOVAR migration
      !!        !  2006-10  (A. Weaver) Cleanup
      !!        !  2007-01  (K. Mogensen) Rewritten with loop
      !!        !  2010-05  (D. Lea) Fix in leap year calculation for NEMO vn3.2
      !!----------------------------------------------------------------------
      !! * Modules used
      USE dom_oce, ONLY : &  ! Geographical information
         & rdt
      USE phycst, ONLY : &   ! Physical constants
         & rday,  &             
         & rmmss, &             
         & rhhmm                        
      !! * Arguments
      INTEGER, INTENT(IN) :: kcycle     ! Current cycle
      INTEGER, INTENT(IN) :: kyea0      ! Initial date coordinates
      INTEGER, INTENT(IN) :: kmon0
      INTEGER, INTENT(IN) :: kday0 
      INTEGER, INTENT(IN) :: khou0
      INTEGER, INTENT(IN) :: kmin0
      INTEGER, INTENT(IN) :: kobsno     ! Number of observations
      INTEGER, INTENT(INOUT) :: kotdobs   ! Number of observations failing time check
      INTEGER, DIMENSION(kobsno), INTENT(IN ) :: &
         & kobsyea,  &      ! Observation time coordinates
         & kobsmon,  &
         & kobsday,  & 
         & kobshou,  &
         & kobsmin
      INTEGER, DIMENSION(kobsno), INTENT(INOUT) :: &
         & kobsqc           ! Quality control flag
      INTEGER, DIMENSION(kobsno), INTENT(OUT) :: &
         & kobsstp          ! Number of time steps up to the 
                            ! observation time

      !! * Local declarations
      INTEGER :: jyea
      INTEGER :: jmon
      INTEGER :: jday
      INTEGER :: jobs
      INTEGER :: iyeastr
      INTEGER :: iyeaend
      INTEGER :: imonstr
      INTEGER :: imonend
      INTEGER :: idaystr
      INTEGER :: idayend 
      INTEGER :: iskip
      INTEGER :: idaystp
      REAL(KIND=wp) :: zminstp
      REAL(KIND=wp) :: zhoustp
      REAL(KIND=wp) :: zobsstp 
      INTEGER, DIMENSION(12) ::   imonth_len    !: length in days of the months of the current year
 
      !-----------------------------------------------------------------------
      ! Initialization
      !-----------------------------------------------------------------------

      ! Intialize the number of time steps per day
      idaystp = NINT( rday / rdt )

      !---------------------------------------------------------------------
      ! Locate the model time coordinates for interpolation
      !---------------------------------------------------------------------

      DO jobs = 1, kobsno

         ! Initialize the time step counter
         kobsstp(jobs) = nit000 - 1

         ! Flag if observation date is less than the initial date

         IF ( ( kobsyea(jobs) < kyea0 )                   &
            & .OR. ( ( kobsyea(jobs) == kyea0 )           &
            &        .AND. ( kobsmon(jobs) <  kmon0 ) )   &
            & .OR. ( ( kobsyea(jobs) == kyea0 )           &
            &        .AND. ( kobsmon(jobs) == kmon0 )     &
            &        .AND. ( kobsday(jobs) <  kday0 ) )   &
            & .OR. ( ( kobsyea(jobs) == kyea0 )           &
            &        .AND. ( kobsmon(jobs) == kmon0 )     &
            &        .AND. ( kobsday(jobs) == kday0 )     &
            &        .AND. ( kobshou(jobs) <  khou0 ) )   &
            &  .OR. ( ( kobsyea(jobs) == kyea0 )          &
            &        .AND. ( kobsmon(jobs) == kmon0 )     &
            &        .AND. ( kobsday(jobs) == kday0 )          &
            &        .AND. ( kobshou(jobs) == khou0 )          &
            &        .AND. ( kobsmin(jobs) <= kmin0 ) ) ) THEN
            kobsstp(jobs) = -1
            kobsqc(jobs)  = IBSET(kobsqc(jobs),13)
            kotdobs       = kotdobs + 1
            CYCLE
         ENDIF

         ! Compute the number of time steps to the observation day
         iyeastr = kyea0
         iyeaend = kobsyea(jobs)
         
         !---------------------------------------------------------------------
         ! Year loop
         !---------------------------------------------------------------------
         DO jyea = iyeastr, iyeaend

            CALL calc_month_len( jyea, imonth_len )
            
            imonstr = 1
            IF ( jyea == kyea0         ) imonstr = kmon0
            imonend = 12
            IF ( jyea == kobsyea(jobs) ) imonend = kobsmon(jobs)
            
            ! Month loop
            DO jmon = imonstr, imonend
               
               idaystr = 1
               IF (       ( jmon == kmon0   ) &
                  & .AND. ( jyea == kyea0   ) ) idaystr = kday0
               idayend = imonth_len(jmon)
               IF (       ( jmon == kobsmon(jobs) ) &
                  & .AND. ( jyea == kobsyea(jobs) ) ) idayend = kobsday(jobs) - 1
               
               ! Day loop
               DO jday = idaystr, idayend
                  kobsstp(jobs) = kobsstp(jobs) + idaystp
               END DO
               
            END DO

         END DO

         ! Add in the number of time steps to the observation minute
         zminstp = rmmss / rdt
         zhoustp = rhhmm * zminstp

         zobsstp =   REAL( kobsmin(jobs) - kmin0, KIND=wp ) * zminstp &
            &      + REAL( kobshou(jobs) - khou0, KIND=wp ) * zhoustp
         kobsstp(jobs) = kobsstp(jobs) + NINT( zobsstp )

         ! Flag if observation step outside the time window
         IF ( ( kobsstp(jobs) < ( nit000 - 1 ) ) &
            & .OR.( kobsstp(jobs) > nitend ) ) THEN
            kobsqc(jobs) = IBSET(kobsqc(jobs),13)
            kotdobs = kotdobs + 1
            CYCLE
         ENDIF

      END DO

   END SUBROUTINE obs_coo_tim

   SUBROUTINE calc_month_len( iyear, imonth_len )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE calc_month_len  ***
      !!          
      !! ** Purpose : Compute the number of days in a months given a year.
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  10-05  (D. Lea)   New routine based on day_init 
      !!----------------------------------------------------------------------

      INTEGER, DIMENSION(12) ::   imonth_len    !: length in days of the months of the current year
      INTEGER :: iyear         !: year
      
      ! length of the month of the current year (from nleapy, read in namelist)
      IF ( nleapy < 2 ) THEN 
         imonth_len(:) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
         IF ( nleapy == 1 ) THEN   ! we are using calendar with leap years
            IF ( MOD(iyear, 4) == 0 .AND. ( MOD(iyear, 400) == 0 .OR. MOD(iyear, 100) /= 0 ) ) THEN
               imonth_len(2) = 29
            ENDIF
         ENDIF
      ELSE
         imonth_len(:) = nleapy   ! all months with nleapy days per year
      ENDIF

   END SUBROUTINE

   SUBROUTINE obs_coo_tim_prof( kcycle,                                   &
      &                    kyea0,   kmon0,   kday0,   khou0,   kmin0,     &
      &                    kobsno,                                        &
      &                    kobsyea, kobsmon, kobsday, kobshou, kobsmin,   &
      &                    ktyp,    kobsqc,  kobsstp, kotdobs, kdailyavtypes, &
      &                    kqc_cutoff )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE obs_coo_tim ***
      !!
      !! ** Purpose : Compute the number of time steps to the observation time.
      !!
      !! ** Method  : For time coordinates ( yea_obs, mon_obs, day_obs, 
      !!              hou_obs, min_obs ), this routine locates the time step 
      !!              that is closest to this time.
      !!
      !! ** Action  : 
      !!
      !! References :
      !!   
      !! History :
      !!        !  1997-07  (A. Weaver) Original
      !!        !  2006-08  (A. Weaver) NEMOVAR migration
      !!        !  2006-10  (A. Weaver) Cleanup
      !!        !  2007-01  (K. Mogensen) Rewritten with loop
      !!----------------------------------------------------------------------
      !! * Modules used
      !! * Arguments
      INTEGER, INTENT(IN) :: kcycle      ! Current cycle
      INTEGER, INTENT(IN) :: kyea0       ! Initial date coordinates
      INTEGER, INTENT(IN) :: kmon0
      INTEGER, INTENT(IN) :: kday0
      INTEGER, INTENT(IN) :: khou0
      INTEGER, INTENT(IN) :: kmin0
      INTEGER, INTENT(IN) :: kobsno      ! Number of observations
      INTEGER, INTENT(INOUT) ::  kotdobs    ! Number of observations failing time check
      INTEGER, DIMENSION(kobsno), INTENT(IN ) :: &
         & kobsyea,  &      ! Observation time coordinates
         & kobsmon,  &
         & kobsday,  & 
         & kobshou,  &
         & kobsmin,  &
         & ktyp             ! Observation type.
      INTEGER, DIMENSION(kobsno), INTENT(INOUT) :: &
         & kobsqc           ! Quality control flag
      INTEGER, DIMENSION(kobsno), INTENT(OUT) :: &
         & kobsstp          ! Number of time steps up to the 
                            ! observation time
      INTEGER, DIMENSION(imaxavtypes), OPTIONAL :: &
         & kdailyavtypes    ! Types for daily averages
      INTEGER, OPTIONAL, INTENT(IN) :: kqc_cutoff           ! QC cutoff value

      !! * Local declarations
      INTEGER :: jobs
      INTEGER :: iqc_cutoff=255

      !-----------------------------------------------------------------------
      ! Call standard obs_coo_tim
      !-----------------------------------------------------------------------

      CALL obs_coo_tim( kcycle, &
         &              kyea0,   kmon0,   kday0,   khou0,   kmin0,     &
         &              kobsno,                                        &
         &              kobsyea, kobsmon, kobsday, kobshou, kobsmin,   &
         &              kobsqc,  kobsstp, kotdobs                      )

      !------------------------------------------------------------------------
      ! Always reject daily averaged data (e.g. MRB data (820)) at initial time
      !------------------------------------------------------------------------

      IF ( PRESENT(kdailyavtypes) ) THEN
         DO jobs = 1, kobsno
            
            IF ( kobsqc(jobs) <= iqc_cutoff ) THEN
               
               IF ( ( kobsstp(jobs) == (nit000 - 1) ).AND.&
                  & ( ANY (kdailyavtypes(:) == ktyp(jobs)) ) ) THEN
                  kobsqc(jobs) = IBSET(kobsqc(jobs),13)
                  kotdobs      = kotdobs + 1
                  CYCLE
               ENDIF
               
            ENDIF
         END DO
      ENDIF


   END SUBROUTINE obs_coo_tim_prof

   SUBROUTINE obs_coo_grd( kobsno, kobsi, kobsj, kobsqc, kgrdobs )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE obs_coo_grd ***
      !!
      !! ** Purpose : Verify that the grid search has not failed
      !!
      !! ** Method  : The previously computed i,j indeces are checked  
      !!
      !! ** Action  : 
      !!
      !! References :
      !!   
      !! History :
      !!        !  2007-01  (K. Mogensen)  Original
      !!----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: kobsno        ! Number of observations
      INTEGER, DIMENSION(kobsno), INTENT(IN ) :: &
         & kobsi, &         ! i,j indeces previously computed
         & kobsj
      INTEGER, INTENT(INOUT) ::  kgrdobs   ! Number of observations failing the check
      INTEGER, DIMENSION(kobsno), INTENT(INOUT) :: &
         & kobsqc           ! Quality control flag

      !! * Local declarations
      INTEGER :: jobs       ! Loop variable

      ! Flag if the grid search failed

      DO jobs = 1, kobsno
         IF ( ( kobsi(jobs) <= 0 ) .AND. ( kobsj(jobs) <= 0 ) ) THEN
            kobsqc(jobs) = IBSET(kobsqc(jobs),12)
            kgrdobs = kgrdobs + 1
         ENDIF
      END DO
      
   END SUBROUTINE obs_coo_grd

   SUBROUTINE obs_coo_spc_2d( kobsno, kpi,     kpj,              &
      &                       kobsi,  kobsj,   pobslam, pobsphi, & 
      &                       plam,   pphi,    pmask,            &
      &                       kobsqc, kosdobs, klanobs,          &
      &                       knlaobs,ld_nea,                    &
      &                       kbdyobs,ld_bound_reject,           &
      &                       kqc_cutoff                         )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE obs_coo_spc_2d  ***
      !!
      !! ** Purpose : Check for points outside the domain and land points
      !!
      !! ** Method  : Remove the observations that are outside the model space
      !!              and time domain or located within model land cells.
      !!
      !! ** Action  : 
      !!   
      !! History :  2007-03  (A. Weaver, K. Mogensen)  Original
      !!         !  2007-06  (K. Mogensen et al) Reject obs. near land.
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in   )                     ::   kobsno            ! Total number of observations
      INTEGER , INTENT(in   )                     ::   kpi    , kpj      ! Number of grid points in (i,j)
      INTEGER , INTENT(in   ), DIMENSION(kobsno)  ::   kobsi  , kobsj    ! Observation (i,j) coordinates
      REAL(wp), INTENT(in   ), DIMENSION(kobsno)  ::   pobslam, pobsphi  ! Observation (lon,lat) coordinates
      REAL(wp), INTENT(in   ), DIMENSION(kpi,kpj) ::   plam   , pphi     ! Model (lon,lat) coordinates
      REAL(wp), INTENT(in   ), DIMENSION(kpi,kpj) ::   pmask             ! Land mask array
      INTEGER , INTENT(inout), DIMENSION(kobsno)  ::   kobsqc            ! Observation quality control
      INTEGER , INTENT(inout)                     ::   kosdobs           ! Observations outside space domain
      INTEGER , INTENT(inout)                     ::   klanobs           ! Observations within a model land cell
      INTEGER , INTENT(inout)                     ::   knlaobs           ! Observations near land
      INTEGER , INTENT(inout)                     ::   kbdyobs           ! Observations near boundary
      LOGICAL , INTENT(in   )                     ::   ld_nea            ! Flag observations near land
      LOGICAL , INTENT(in   )                     ::   ld_bound_reject   ! Flag observations near open boundary 
      INTEGER , INTENT(in   )                     ::   kqc_cutoff        ! Cutoff QC value
      !
      REAL(KIND=wp), DIMENSION(2,2,kobsno) ::   zgmsk          ! Grid mask
      REAL(KIND=wp), DIMENSION(2,2,kobsno) ::   zbmsk          ! Boundary mask
      REAL(KIND=wp), DIMENSION(jpi,jpj)    ::   zbdymask
      REAL(KIND=wp), DIMENSION(2,2,kobsno) ::   zglam, zgphi   ! Model Lon/lat at grid points
      INTEGER      , DIMENSION(2,2,kobsno) ::   igrdi, igrdj   ! Grid i,j
      LOGICAL ::   lgridobs           ! Is observation on a model grid point.
      INTEGER ::   iig, ijg           ! i,j of observation on model grid point.
      INTEGER ::   jobs, ji, jj
      !!----------------------------------------------------------------------
      
      ! Get grid point indices

      DO jobs = 1, kobsno
         
         ! For invalid points use 2,2

         IF ( kobsqc(jobs) >= kqc_cutoff ) THEN

            igrdi(1,1,jobs) = 1
            igrdj(1,1,jobs) = 1
            igrdi(1,2,jobs) = 1
            igrdj(1,2,jobs) = 2
            igrdi(2,1,jobs) = 2
            igrdj(2,1,jobs) = 1
            igrdi(2,2,jobs) = 2
            igrdj(2,2,jobs) = 2

         ELSE

            igrdi(1,1,jobs) = kobsi(jobs)-1
            igrdj(1,1,jobs) = kobsj(jobs)-1
            igrdi(1,2,jobs) = kobsi(jobs)-1
            igrdj(1,2,jobs) = kobsj(jobs)
            igrdi(2,1,jobs) = kobsi(jobs)
            igrdj(2,1,jobs) = kobsj(jobs)-1
            igrdi(2,2,jobs) = kobsi(jobs)
            igrdj(2,2,jobs) = kobsj(jobs)

         ENDIF

      END DO

      IF (ln_bdy) THEN
        ! Create a mask grid points in boundary rim
        IF (ld_bound_reject) THEN
           zbdymask(:,:) = 1.0_wp
           DO ji = 1, nb_bdy
              DO jj = 1, idx_bdy(ji)%nblen(1)
                 zbdymask(idx_bdy(ji)%nbi(jj,1),idx_bdy(ji)%nbj(jj,1)) = 0.0_wp
              ENDDO
           ENDDO

           CALL obs_int_comm_2d( 2, 2, kobsno, kpi, kpj, igrdi, igrdj, zbdymask, zbmsk )
        ENDIF
      ENDIF

      
      CALL obs_int_comm_2d( 2, 2, kobsno, kpi, kpj, igrdi, igrdj, pmask, zgmsk )
      CALL obs_int_comm_2d( 2, 2, kobsno, kpi, kpj, igrdi, igrdj, plam, zglam )
      CALL obs_int_comm_2d( 2, 2, kobsno, kpi, kpj, igrdi, igrdj, pphi, zgphi )

      DO jobs = 1, kobsno

         ! Skip bad observations
         IF ( kobsqc(jobs) >= kqc_cutoff ) CYCLE

         ! Flag if the observation falls outside the model spatial domain
         IF (       ( pobslam(jobs) < -180. ) &
            &  .OR. ( pobslam(jobs) >  180. ) &
            &  .OR. ( pobsphi(jobs) <  -90. ) &
            &  .OR. ( pobsphi(jobs) >   90. ) ) THEN
            kobsqc(jobs) = IBSET(kobsqc(jobs),11)
            kosdobs = kosdobs + 1
            CYCLE
         ENDIF

         ! Flag if the observation falls with a model land cell
         IF ( SUM( zgmsk(1:2,1:2,jobs) ) == 0.0_wp ) THEN
            kobsqc(jobs) = IBSET(kobsqc(jobs),10)
            klanobs = klanobs + 1
            CYCLE
         ENDIF

         ! Check if this observation is on a grid point

         lgridobs = .FALSE.
         iig = -1
         ijg = -1
         DO jj = 1, 2
            DO ji = 1, 2
               IF ( ( ABS( zgphi(ji,jj,jobs) - pobsphi(jobs) ) < 1.0e-6_wp ) &
                  & .AND. &
                  & ( ABS( MOD( zglam(ji,jj,jobs) - pobslam(jobs),360.0) )  &
                  & < 1.0e-6_wp ) ) THEN
                  lgridobs = .TRUE.
                  iig = ji
                  ijg = jj
               ENDIF
            END DO
         END DO
  
         ! For observations on the grid reject them if their are at
         ! a masked point
         
         IF (lgridobs) THEN
            IF (zgmsk(iig,ijg,jobs) == 0.0_wp ) THEN
               kobsqc(jobs) = IBSET(kobsqc(jobs),10)
               klanobs = klanobs + 1
               CYCLE
            ENDIF
         ENDIF
                      
         ! Flag if the observation falls is close to land
         IF ( MINVAL( zgmsk(1:2,1:2,jobs) ) == 0.0_wp) THEN
            knlaobs = knlaobs + 1
            IF (ld_nea) THEN
               kobsqc(jobs) = IBSET(kobsqc(jobs),9)
               CYCLE
            ENDIF
         ENDIF

         IF (ln_bdy) THEN
         ! Flag if the observation falls close to the boundary rim
           IF (ld_bound_reject) THEN
              IF ( MINVAL( zbmsk(1:2,1:2,jobs) ) == 0.0_wp ) THEN
                 kobsqc(jobs) = IBSET(kobsqc(jobs),8)
                 kbdyobs = kbdyobs + 1
                 CYCLE
              ENDIF
              ! for observations on the grid...
              IF (lgridobs) THEN
                 IF (zbmsk(iig,ijg,jobs) == 0.0_wp ) THEN
                    kobsqc(jobs) = IBSET(kobsqc(jobs),8)
                    kbdyobs = kbdyobs + 1
                    CYCLE
                 ENDIF
              ENDIF
            ENDIF
         ENDIF
         !
      END DO
      !
   END SUBROUTINE obs_coo_spc_2d


   SUBROUTINE obs_coo_spc_3d( kprofno, kobsno,  kpstart, kpend, &
      &                       kpi,     kpj,     kpk,            &
      &                       kobsi,   kobsj,   kobsk,          &
      &                       pobslam, pobsphi, pobsdep,        &
      &                       plam,    pphi,    pdep,    pmask, &
      &                       kpobsqc, kobsqc,  kosdobs,        &
      &                       klanobs, knlaobs, ld_nea,         &
      &                       kbdyobs, ld_bound_reject,         &
      &                       kqc_cutoff                        )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE obs_coo_spc_3d  ***
      !!
      !! ** Purpose : Check for points outside the domain and land points
      !!              Reset depth of observation above highest model level
      !!              to the value of highest model level
      !!
      !! ** Method  : Remove the observations that are outside the model space
      !!              and time domain or located within model land cells.
      !!
      !!              NB. T and S profile observations lying between the ocean
      !!              surface and the depth of the first model T point are 
      !!              assigned a depth equal to that of the first model T pt.
      !!
      !! ** Action  : 
      !!   
      !! History :
      !!        !  2007-01  (K. Mogensen) Rewrite of parts of obs_scr
      !!        !  2007-06  (K. Mogensen et al) Reject obs. near land.
      !!----------------------------------------------------------------------
      !! * Modules used
      USE dom_oce, ONLY : &       ! Geographical information
         & gdepw_1d,      &
         & gdepw_0,       &                       
         & gdepw_n,       &
         & gdept_n,       &
         & ln_zco,        &
         & ln_zps             

      !! * Arguments
      INTEGER, INTENT(IN) :: kprofno      ! Number of profiles
      INTEGER, INTENT(IN) :: kobsno       ! Total number of observations
      INTEGER, INTENT(IN) :: kpi          ! Number of grid points in (i,j,k)
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kpk    
      INTEGER, DIMENSION(kprofno), INTENT(IN) :: &
         & kpstart, &         ! Start of individual profiles
         & kpend              ! End of individual profiles
      INTEGER, DIMENSION(kprofno), INTENT(IN) :: &
         & kobsi, &           ! Observation (i,j) coordinates
         & kobsj
      INTEGER, DIMENSION(kobsno), INTENT(IN) :: &
         & kobsk              ! Observation k coordinate
      REAL(KIND=wp), DIMENSION(kprofno), INTENT(IN) :: &
         & pobslam, &         ! Observation (lon,lat) coordinates
         & pobsphi
      REAL(KIND=wp), DIMENSION(kobsno), INTENT(INOUT) :: &
         & pobsdep            ! Observation depths  
      REAL(KIND=wp), DIMENSION(kpi,kpj), INTENT(IN) :: &
         & plam, pphi         ! Model (lon,lat) coordinates
      REAL(KIND=wp), DIMENSION(kpk), INTENT(IN) :: &
         & pdep               ! Model depth coordinates
      REAL(KIND=wp), DIMENSION(kpi,kpj,kpk), INTENT(IN) :: &
         & pmask              ! Land mask array
      INTEGER, DIMENSION(kprofno), INTENT(INOUT) :: &
         & kpobsqc             ! Profile quality control
      INTEGER, DIMENSION(kobsno), INTENT(INOUT) :: &
         & kobsqc             ! Observation quality control
      INTEGER, INTENT(INOUT) :: kosdobs     ! Observations outside space domain
      INTEGER, INTENT(INOUT) :: klanobs     ! Observations within a model land cell
      INTEGER, INTENT(INOUT) :: knlaobs     ! Observations near land
      INTEGER, INTENT(INOUT) :: kbdyobs     ! Observations near boundary
      LOGICAL, INTENT(IN) :: ld_nea         ! Flag observations near land
      LOGICAL, INTENT(IN) :: ld_bound_reject  ! Flag observations near open boundary
      INTEGER, INTENT(IN) :: kqc_cutoff     ! Cutoff QC value

      !! * Local declarations
      REAL(KIND=wp), DIMENSION(2,2,kpk,kprofno) :: &
         & zgmsk              ! Grid mask
      REAL(KIND=wp), DIMENSION(2,2,kprofno) :: &
         & zbmsk              ! Boundary mask
      REAL(KIND=wp), DIMENSION(jpi,jpj) :: zbdymask
      REAL(KIND=wp), DIMENSION(2,2,kpk,kprofno) :: &
         & zgdepw
      REAL(KIND=wp), DIMENSION(2,2,kprofno) :: &
         & zglam, &           ! Model longitude at grid points
         & zgphi              ! Model latitude at grid points
      INTEGER, DIMENSION(2,2,kprofno) :: &
         & igrdi, &           ! Grid i,j
         & igrdj
      LOGICAL :: lgridobs           ! Is observation on a model grid point.
      LOGICAL :: ll_next_to_land    ! Is a profile next to land 
      INTEGER :: iig, ijg           ! i,j of observation on model grid point.
      INTEGER :: jobs, jobsp, jk, ji, jj
      !!----------------------------------------------------------------------

      ! Get grid point indices
      
      DO jobs = 1, kprofno

         ! For invalid points use 2,2

         IF ( kpobsqc(jobs) >= kqc_cutoff ) THEN
            
            igrdi(1,1,jobs) = 1
            igrdj(1,1,jobs) = 1
            igrdi(1,2,jobs) = 1
            igrdj(1,2,jobs) = 2
            igrdi(2,1,jobs) = 2
            igrdj(2,1,jobs) = 1
            igrdi(2,2,jobs) = 2
            igrdj(2,2,jobs) = 2
            
         ELSE
            
            igrdi(1,1,jobs) = kobsi(jobs)-1
            igrdj(1,1,jobs) = kobsj(jobs)-1
            igrdi(1,2,jobs) = kobsi(jobs)-1
            igrdj(1,2,jobs) = kobsj(jobs)
            igrdi(2,1,jobs) = kobsi(jobs)
            igrdj(2,1,jobs) = kobsj(jobs)-1
            igrdi(2,2,jobs) = kobsi(jobs)
            igrdj(2,2,jobs) = kobsj(jobs)
            
         ENDIF
         
      END DO

      IF (ln_bdy) THEN
        ! Create a mask grid points in boundary rim
        IF (ld_bound_reject) THEN           
           zbdymask(:,:) = 1.0_wp
           DO ji = 1, nb_bdy
              DO jj = 1, idx_bdy(ji)%nblen(1)
                 zbdymask(idx_bdy(ji)%nbi(jj,1),idx_bdy(ji)%nbj(jj,1)) = 0.0_wp
              ENDDO
           ENDDO
        ENDIF
   
        CALL obs_int_comm_2d( 2, 2, kprofno, kpi, kpj, igrdi, igrdj, zbdymask, zbmsk )
      ENDIF
      
      CALL obs_int_comm_3d( 2, 2, kprofno, kpi, kpj, kpk, igrdi, igrdj, pmask, zgmsk )
      CALL obs_int_comm_2d( 2, 2, kprofno, kpi, kpj, igrdi, igrdj, plam, zglam )
      CALL obs_int_comm_2d( 2, 2, kprofno, kpi, kpj, igrdi, igrdj, pphi, zgphi )
      CALL obs_int_comm_3d( 2, 2, kprofno, kpi, kpj, kpk, igrdi, igrdj, gdepw_n(:,:,:), &
        &                     zgdepw )

      DO jobs = 1, kprofno

         ! Skip bad profiles
         IF ( kpobsqc(jobs) >= kqc_cutoff ) CYCLE

         ! Check if this observation is on a grid point

         lgridobs = .FALSE.
         iig = -1
         ijg = -1
         DO jj = 1, 2
            DO ji = 1, 2
               IF ( ( ABS( zgphi(ji,jj,jobs) - pobsphi(jobs) ) < 1.0e-6_wp ) &
                  & .AND. &
                  & ( ABS( MOD( zglam(ji,jj,jobs) - pobslam(jobs),360.0) ) < 1.0e-6_wp ) &
                  & ) THEN
                  lgridobs = .TRUE.
                  iig = ji
                  ijg = jj
               ENDIF
            END DO
         END DO

         ! Check if next to land
         IF (  ANY( zgmsk(1:2,1:2,1,jobs) == 0.0_wp ) ) THEN
            ll_next_to_land=.TRUE.
         ELSE
            ll_next_to_land=.FALSE.
         ENDIF 

         ! Reject observations

         DO jobsp = kpstart(jobs), kpend(jobs)

            ! Flag if the observation falls outside the model spatial domain
            IF (       ( pobslam(jobs) < -180.         )       &
               &  .OR. ( pobslam(jobs) >  180.         )       &
               &  .OR. ( pobsphi(jobs) <  -90.         )       &
               &  .OR. ( pobsphi(jobs) >   90.         )       &
               &  .OR. ( pobsdep(jobsp) < 0.0          )       &
               &  .OR. ( pobsdep(jobsp) > gdepw_1d(kpk)) ) THEN
               kobsqc(jobsp) = IBSET(kobsqc(jobsp),11)
               kosdobs = kosdobs + 1
               CYCLE
            ENDIF

            ! To check if an observations falls within land:
             
            ! Flag if the observation is deeper than the bathymetry
            ! Or if it is within the mask
            IF ( ALL( zgdepw(1:2,1:2,kpk,jobs) < pobsdep(jobsp) ) &
               &     .OR. &
               &  ( SUM( zgmsk(1:2,1:2,kobsk(jobsp)-1:kobsk(jobsp),jobs) ) &
               &  == 0.0_wp) ) THEN
               kobsqc(jobsp) = IBSET(kobsqc(jobsp),10)
               klanobs = klanobs + 1
               CYCLE
            ENDIF
               
            ! Flag if the observation is close to land
            IF ( ll_next_to_land ) THEN
               knlaobs = knlaobs + 1
               IF (ld_nea) THEN   
                  kobsqc(jobsp) = IBSET(kobsqc(jobsp),10)
               ENDIF 
            ENDIF
            
            ! For observations on the grid reject them if their are at
            ! a masked point
            
            IF (lgridobs) THEN
               IF (zgmsk(iig,ijg,kobsk(jobsp)-1,jobs) == 0.0_wp ) THEN
                  kobsqc(jobsp) = IBSET(kobsqc(jobs),10)
                  klanobs = klanobs + 1
                  CYCLE
               ENDIF
            ENDIF
            
            ! Flag if the observation falls is close to land
            IF ( MINVAL( zgmsk(1:2,1:2,kobsk(jobsp)-1:kobsk(jobsp),jobs) ) == &
               &  0.0_wp) THEN
               IF (ld_nea) kobsqc(jobsp) = kobsqc(jobsp) + 14
               knlaobs = knlaobs + 1
            ENDIF

            ! Set observation depth equal to that of the first model depth
            IF ( pobsdep(jobsp) <= pdep(1) ) THEN
               pobsdep(jobsp) = pdep(1)  
            ENDIF
            
            IF (ln_bdy) THEN
               ! Flag if the observation falls close to the boundary rim
               IF (ld_bound_reject) THEN
                  IF ( MINVAL( zbmsk(1:2,1:2,jobs) ) == 0.0_wp ) THEN
                     kobsqc(jobsp) = IBSET(kobsqc(jobs),8)
                     kbdyobs = kbdyobs + 1
                     CYCLE
                  ENDIF
                  ! for observations on the grid...
                  IF (lgridobs) THEN
                     IF (zbmsk(iig,ijg,jobs) == 0.0_wp ) THEN
                        kobsqc(jobsp) = IBSET(kobsqc(jobs),8)
                        kbdyobs = kbdyobs + 1
                        CYCLE
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            !
         END DO
      END DO
      !
   END SUBROUTINE obs_coo_spc_3d


   SUBROUTINE obs_pro_rej( profdata, kqc_cutoff )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE obs_pro_rej ***
      !!
      !! ** Purpose : Reject all data within a rejected profile
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!
      !! References :
      !!   
      !! History :   2007-10  (K. Mogensen) Original code
      !!----------------------------------------------------------------------
      TYPE(obs_prof), INTENT(inout) ::   profdata     ! Profile data
      INTEGER       , INTENT(in   ) ::   kqc_cutoff   ! QC cutoff value
      !
      INTEGER :: jprof
      INTEGER :: jvar
      INTEGER :: jobs
      !!----------------------------------------------------------------------
      
      ! Loop over profiles

      DO jprof = 1, profdata%nprof

         IF ( profdata%nqc(jprof) > kqc_cutoff ) THEN
            
            DO jvar = 1, profdata%nvar

               DO jobs = profdata%npvsta(jprof,jvar),  &
                  &      profdata%npvend(jprof,jvar)
                  
                  profdata%var(jvar)%nvqc(jobs) = &
                     & IBSET(profdata%var(jvar)%nvqc(jobs),14)

               END DO

            END DO

         ENDIF

      END DO
      !
   END SUBROUTINE obs_pro_rej


   SUBROUTINE obs_uv_rej( profdata, knumu, knumv, kqc_cutoff )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE obs_uv_rej ***
      !!
      !! ** Purpose : Reject u if v is rejected and vice versa
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!
      !! References :
      !!   
      !! History :   2009-2  (K. Mogensen) Original code
      !!----------------------------------------------------------------------
      TYPE(obs_prof), INTENT(INOUT) :: profdata   ! Profile data
      INTEGER, INTENT(INOUT) :: knumu             ! Number of u rejected
      INTEGER, INTENT(INOUT) :: knumv             ! Number of v rejected
      INTEGER, INTENT(IN) :: kqc_cutoff           ! QC cutoff value
      !
      INTEGER :: jprof
      INTEGER :: jvar
      INTEGER :: jobs
      !!----------------------------------------------------------------------

      DO jprof = 1, profdata%nprof      !==  Loop over profiles  ==!
         !
         IF ( ( profdata%npvsta(jprof,1) /= profdata%npvsta(jprof,2) ) .OR. &
            & ( profdata%npvend(jprof,1) /= profdata%npvend(jprof,2) ) ) THEN
            !
            CALL ctl_stop('U,V profiles inconsistent in obs_uv_rej')
            RETURN
            !
         ENDIF
         !
         DO jobs = profdata%npvsta(jprof,1), profdata%npvend(jprof,1)
            !  
            IF ( ( profdata%var(1)%nvqc(jobs) >  kqc_cutoff ) .AND. &
               & ( profdata%var(2)%nvqc(jobs) <=  kqc_cutoff) ) THEN
               profdata%var(2)%nvqc(jobs) = IBSET(profdata%var(1)%nvqc(jobs),15)
               knumv = knumv + 1
            ENDIF
            IF ( ( profdata%var(2)%nvqc(jobs) >  kqc_cutoff ) .AND. &
               & ( profdata%var(1)%nvqc(jobs) <=  kqc_cutoff) ) THEN
               profdata%var(1)%nvqc(jobs) = IBSET(profdata%var(1)%nvqc(jobs),15)
               knumu = knumu + 1
            ENDIF
            !
         END DO
         !
      END DO
      !
   END SUBROUTINE obs_uv_rej

   !!=====================================================================
END MODULE obs_prep
