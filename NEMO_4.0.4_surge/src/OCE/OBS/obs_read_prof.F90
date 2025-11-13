MODULE obs_read_prof
   !!======================================================================
   !!                       ***  MODULE obs_read_prof  ***
   !! Observation diagnostics: Read the T and S profile observations
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   obs_rea_pro_dri : Driver for reading profile obs
   !!----------------------------------------------------------------------

   !! * Modules used   
   USE par_kind                 ! Precision variables
   USE par_oce                  ! Ocean parameters
   USE in_out_manager           ! I/O manager
   USE dom_oce                  ! Ocean space and time domain variables
   USE obs_mpp                  ! MPP support routines for observation diagnostics
   USE julian                   ! Julian date routines
   USE obs_utils                ! Observation operator utility functions
   USE obs_prep                 ! Prepare observation arrays
   USE obs_grid                 ! Grid search
   USE obs_sort                 ! Sorting observation arrays
   USE obs_profiles_def         ! Profile definitions
   USE obs_conv                 ! Various conversion routines
   USE obs_types                ! Observation type definitions
   USE netcdf                   ! NetCDF library
   USE obs_oper                 ! Observation operators
   USE lib_mpp                  ! For ctl_warn/stop
   USE obs_fbm                  ! Feedback routines

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE

   PUBLIC obs_rea_prof  ! Read the profile observations 

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE obs_rea_prof( profdata, knumfiles, cdfilenames, &
      &                     kvars, kextr, kstp, ddobsini, ddobsend, &
      &                     ldvar1, ldvar2, ldignmis, ldsatt, &
      &                     ldmod, kdailyavtypes )
      !!---------------------------------------------------------------------
      !!
      !!                   *** ROUTINE obs_rea_prof ***
      !!
      !! ** Purpose : Read from file the profile observations
      !!
      !! ** Method  : Read feedback data in and transform to NEMO internal 
      !!              profile data structure
      !!
      !! ** Action  : 
      !!
      !! References : 
      !!
      !! History :  
      !!      ! :  2009-09 (K. Mogensen) : New merged version of old routines
      !!      ! :  2015-08 (M. Martin) : Merged profile and velocity routines
      !!----------------------------------------------------------------------

      !! * Arguments
      TYPE(obs_prof), INTENT(OUT) :: &
         & profdata                     ! Profile data to be read
      INTEGER, INTENT(IN) :: knumfiles  ! Number of files to read
      CHARACTER(LEN=128), INTENT(IN) ::  &
         & cdfilenames(knumfiles)        ! File names to read in
      INTEGER, INTENT(IN) :: kvars      ! Number of variables in profdata
      INTEGER, INTENT(IN) :: kextr      ! Number of extra fields for each var
      INTEGER, INTENT(IN) :: kstp       ! Ocean time-step index
      LOGICAL, INTENT(IN) :: ldvar1     ! Observed variables switches
      LOGICAL, INTENT(IN) :: ldvar2
      LOGICAL, INTENT(IN) :: ldignmis   ! Ignore missing files
      LOGICAL, INTENT(IN) :: ldsatt     ! Compute salinity at all temperature points
      LOGICAL, INTENT(IN) :: ldmod      ! Initialize model from input data
      REAL(dp), INTENT(IN) :: ddobsini  ! Obs. ini time in YYYYMMDD.HHMMSS
      REAL(dp), INTENT(IN) :: ddobsend  ! Obs. end time in YYYYMMDD.HHMMSS
      INTEGER, DIMENSION(imaxavtypes), OPTIONAL :: &
         & kdailyavtypes                ! Types of daily average observations

      !! * Local declarations
      CHARACTER(LEN=15), PARAMETER :: cpname='obs_rea_prof'
      CHARACTER(len=8) :: clrefdate
      CHARACTER(len=8), DIMENSION(:), ALLOCATABLE :: clvars
      INTEGER :: jvar
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jk
      INTEGER :: ij
      INTEGER :: iflag
      INTEGER :: inobf
      INTEGER :: i_file_id
      INTEGER :: inowin
      INTEGER :: iyea
      INTEGER :: imon
      INTEGER :: iday
      INTEGER :: ihou
      INTEGER :: imin
      INTEGER :: isec
      INTEGER :: iprof
      INTEGER :: iproftot
      INTEGER :: ivar1t0
      INTEGER :: ivar2t0
      INTEGER :: ivar1t
      INTEGER :: ivar2t
      INTEGER :: ip3dt
      INTEGER :: ios
      INTEGER :: ioserrcount
      INTEGER :: ivar1tmpp
      INTEGER :: ivar2tmpp
      INTEGER :: ip3dtmpp
      INTEGER :: itype
      INTEGER, DIMENSION(knumfiles) :: &
         & irefdate
      INTEGER, DIMENSION(ntyp1770+1) :: &
         & itypvar1,    &
         & itypvar1mpp, &
         & itypvar2,    &
         & itypvar2mpp 
      INTEGER, DIMENSION(:), ALLOCATABLE :: &
         & iobsi1,    &
         & iobsj1,    &
         & iproc1,    &
         & iobsi2,    &
         & iobsj2,    &
         & iproc2,    &
         & iindx,    &
         & ifileidx, &
         & iprofidx
      INTEGER, DIMENSION(imaxavtypes) :: &
         & idailyavtypes
      INTEGER, DIMENSION(kvars) :: &
         & iv3dt
      REAL(wp), DIMENSION(:), ALLOCATABLE :: &
         & zphi, &
         & zlam
      REAL(wp), DIMENSION(:), ALLOCATABLE :: &
         & zdat
      REAL(wp), DIMENSION(knumfiles) :: &
         & djulini, &
         & djulend
      LOGICAL :: llvalprof
      LOGICAL :: lldavtimset
      TYPE(obfbdata), POINTER, DIMENSION(:) :: &
         & inpfiles

      ! Local initialization
      iprof = 0
      ivar1t0 = 0
      ivar2t0 = 0
      ip3dt = 0

      ! Daily average types
      lldavtimset = .FALSE.
      IF ( PRESENT(kdailyavtypes) ) THEN
         idailyavtypes(:) = kdailyavtypes(:)
         IF ( ANY (idailyavtypes(:) /= -1) ) lldavtimset = .TRUE.
      ELSE
         idailyavtypes(:) = -1
      ENDIF

      !-----------------------------------------------------------------------
      ! Count the number of files needed and allocate the obfbdata type
      !-----------------------------------------------------------------------

      inobf = knumfiles

      ALLOCATE( inpfiles(inobf) )

      prof_files : DO jj = 1, inobf

         !---------------------------------------------------------------------
         ! Prints
         !---------------------------------------------------------------------
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) ' obs_rea_pro_dri : Reading from file = ', &
               & TRIM( TRIM( cdfilenames(jj) ) )
            WRITE(numout,*) ' ~~~~~~~~~~~~~~~'
            WRITE(numout,*)
         ENDIF

         !---------------------------------------------------------------------
         !  Initialization: Open file and get dimensions only
         !---------------------------------------------------------------------

         iflag = nf90_open( TRIM( cdfilenames(jj) ), nf90_nowrite, &
            &                      i_file_id )

         IF ( iflag /= nf90_noerr ) THEN

            IF ( ldignmis ) THEN
               inpfiles(jj)%nobs = 0
               CALL ctl_warn( 'File ' // TRIM( cdfilenames(jj) ) // &
                  &           ' not found' )
            ELSE 
               CALL ctl_stop( 'File ' // TRIM( cdfilenames(jj) ) // &
                  &           ' not found' )
            ENDIF

         ELSE 

            !------------------------------------------------------------------
            !  Close the file since it is opened in read_obfbdata
            !------------------------------------------------------------------

            iflag = nf90_close( i_file_id )

            !------------------------------------------------------------------
            !  Read the profile file into inpfiles
            !------------------------------------------------------------------
            CALL init_obfbdata( inpfiles(jj) )
            CALL read_obfbdata( TRIM( cdfilenames(jj) ), inpfiles(jj), &
               &                ldgrid = .TRUE. )

            IF ( inpfiles(jj)%nvar < 2 ) THEN
               CALL ctl_stop( 'Feedback format error: ', &
                  &           ' less than 2 vars in profile file' )
            ENDIF

            IF ( ldmod .AND. ( inpfiles(jj)%nadd == 0 ) ) THEN
               CALL ctl_stop( 'Model not in input data' )
            ENDIF

            IF ( jj == 1 ) THEN
               ALLOCATE( clvars( inpfiles(jj)%nvar ) )
               DO ji = 1, inpfiles(jj)%nvar
                 clvars(ji) = inpfiles(jj)%cname(ji)
               END DO
            ELSE
               DO ji = 1, inpfiles(jj)%nvar
                  IF ( inpfiles(jj)%cname(ji) /= clvars(ji) ) THEN
                     CALL ctl_stop( 'Feedback file variables not consistent', &
                        &           ' with previous files for this type' )
                  ENDIF
               END DO
            ENDIF

            !------------------------------------------------------------------
            !  Change longitude (-180,180)
            !------------------------------------------------------------------

            DO ji = 1, inpfiles(jj)%nobs 

               IF ( inpfiles(jj)%plam(ji) < -180. ) &
                  &   inpfiles(jj)%plam(ji) = inpfiles(jj)%plam(ji) + 360.

               IF ( inpfiles(jj)%plam(ji) >  180. ) &
                  &   inpfiles(jj)%plam(ji) = inpfiles(jj)%plam(ji) - 360.

            END DO

            !------------------------------------------------------------------
            !  Calculate the date  (change eventually)
            !------------------------------------------------------------------
            clrefdate=inpfiles(jj)%cdjuldref(1:8)
            READ(clrefdate,'(I8)') irefdate(jj)

            CALL ddatetoymdhms( ddobsini, iyea, imon, iday, ihou, imin, isec )
            CALL greg2jul( isec, imin, ihou, iday, imon, iyea, djulini(jj), &
               &           krefdate = irefdate(jj) )
            CALL ddatetoymdhms( ddobsend, iyea, imon, iday, ihou, imin, isec )
            CALL greg2jul( isec, imin, ihou, iday, imon, iyea, djulend(jj), &
               &           krefdate = irefdate(jj) )

            ioserrcount=0
            IF ( lldavtimset ) THEN

               IF ( ANY ( idailyavtypes(:) /= -1 ) .AND. lwp) THEN
                  WRITE(numout,*)' Resetting time of daily averaged', &
                     &           ' observations to the end of the day'
               ENDIF

               DO ji = 1, inpfiles(jj)%nobs
                  READ( inpfiles(jj)%cdtyp(ji), '(I4)', IOSTAT = ios, ERR = 900 ) itype
900               IF ( ios /= 0 ) THEN
                     ! Set type to zero if there is a problem in the string conversion
                     itype = 0
                  ENDIF

                  IF ( ANY ( idailyavtypes(:) == itype ) ) THEN
                  !  for daily averaged data force the time
                  !  to be the last time-step of the day, but still within the day.
                     IF ( inpfiles(jj)%ptim(ji) >= 0. ) THEN
                        inpfiles(jj)%ptim(ji) = &
                           & INT(inpfiles(jj)%ptim(ji)) + 0.9999
                     ELSE
                        inpfiles(jj)%ptim(ji) = &
                           & INT(inpfiles(jj)%ptim(ji)) - 0.0001
                     ENDIF
                  ENDIF

               END DO

            ENDIF

            IF ( inpfiles(jj)%nobs > 0 ) THEN
               inpfiles(jj)%iproc(:,:) = -1
               inpfiles(jj)%iobsi(:,:) = -1
               inpfiles(jj)%iobsj(:,:) = -1
            ENDIF
            inowin = 0
            DO ji = 1, inpfiles(jj)%nobs
               IF ( BTEST(inpfiles(jj)%ioqc(ji),2 ) ) CYCLE
               IF ( BTEST(inpfiles(jj)%ivqc(ji,1),2) .AND. &
                  & BTEST(inpfiles(jj)%ivqc(ji,2),2) ) CYCLE
               IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND. &
                  & ( inpfiles(jj)%ptim(ji) <= djulend(jj) )       ) THEN
                  inowin = inowin + 1
               ENDIF
            END DO
            ALLOCATE( zlam(inowin)  )
            ALLOCATE( zphi(inowin)  )
            ALLOCATE( iobsi1(inowin) )
            ALLOCATE( iobsj1(inowin) )
            ALLOCATE( iproc1(inowin) )
            ALLOCATE( iobsi2(inowin) )
            ALLOCATE( iobsj2(inowin) )
            ALLOCATE( iproc2(inowin) )
            inowin = 0
            DO ji = 1, inpfiles(jj)%nobs
               IF ( BTEST(inpfiles(jj)%ioqc(ji),2 ) ) CYCLE
               IF ( BTEST(inpfiles(jj)%ivqc(ji,1),2) .AND. &
                  & BTEST(inpfiles(jj)%ivqc(ji,2),2) ) CYCLE
               IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND. &
                  & ( inpfiles(jj)%ptim(ji) <= djulend(jj) )       ) THEN
                  inowin = inowin + 1
                  zlam(inowin) = inpfiles(jj)%plam(ji)
                  zphi(inowin) = inpfiles(jj)%pphi(ji)
               ENDIF
            END DO

            IF ( TRIM( inpfiles(jj)%cname(1) ) == 'POTM' ) THEN
               CALL obs_grid_search( inowin, zlam, zphi, iobsi1, iobsj1, &
                  &                  iproc1, 'T' )
               iobsi2(:) = iobsi1(:)
               iobsj2(:) = iobsj1(:)
               iproc2(:) = iproc1(:)
            ELSEIF ( TRIM( inpfiles(jj)%cname(1) ) == 'UVEL' ) THEN
               CALL obs_grid_search( inowin, zlam, zphi, iobsi1, iobsj1, &
                  &                  iproc1, 'U' )
               CALL obs_grid_search( inowin, zlam, zphi, iobsi2, iobsj2, &
                  &                  iproc2, 'V' )
            ENDIF

            inowin = 0
            DO ji = 1, inpfiles(jj)%nobs
               IF ( BTEST(inpfiles(jj)%ioqc(ji),2 ) ) CYCLE
               IF ( BTEST(inpfiles(jj)%ivqc(ji,1),2) .AND. &
                  & BTEST(inpfiles(jj)%ivqc(ji,2),2) ) CYCLE
               IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND. &
                  & ( inpfiles(jj)%ptim(ji) <= djulend(jj) )       ) THEN
                  inowin = inowin + 1
                  inpfiles(jj)%iproc(ji,1) = iproc1(inowin)
                  inpfiles(jj)%iobsi(ji,1) = iobsi1(inowin)
                  inpfiles(jj)%iobsj(ji,1) = iobsj1(inowin)
                  inpfiles(jj)%iproc(ji,2) = iproc2(inowin)
                  inpfiles(jj)%iobsi(ji,2) = iobsi2(inowin)
                  inpfiles(jj)%iobsj(ji,2) = iobsj2(inowin)
                  IF ( inpfiles(jj)%iproc(ji,1) /= &
                     & inpfiles(jj)%iproc(ji,2) ) THEN
                     CALL ctl_stop( 'Error in obs_read_prof:', &
                        & 'var1 and var2 observation on different processors')
                  ENDIF
               ENDIF
            END DO
            DEALLOCATE( zlam, zphi, iobsi1, iobsj1, iproc1, iobsi2, iobsj2, iproc2 )

            DO ji = 1, inpfiles(jj)%nobs
               IF ( BTEST(inpfiles(jj)%ioqc(ji),2 ) ) CYCLE
               IF ( BTEST(inpfiles(jj)%ivqc(ji,1),2) .AND. &
                  & BTEST(inpfiles(jj)%ivqc(ji,2),2) ) CYCLE
               IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND. &
                  & ( inpfiles(jj)%ptim(ji) <= djulend(jj) )       ) THEN
                  IF ( nproc == 0 ) THEN
                     IF ( inpfiles(jj)%iproc(ji,1) >  nproc ) CYCLE
                  ELSE
                     IF ( inpfiles(jj)%iproc(ji,1) /= nproc ) CYCLE
                  ENDIF
                  llvalprof = .FALSE.
                  IF ( ldvar1 ) THEN
                     loop_t_count : DO ij = 1,inpfiles(jj)%nlev
                        IF ( inpfiles(jj)%pdep(ij,ji) >= 6000. ) &
                           & CYCLE
                        IF ( .NOT. BTEST(inpfiles(jj)%ivlqc(ij,ji,1),2) .AND. &
                           & .NOT. BTEST(inpfiles(jj)%idqc(ij,ji),2) ) THEN
                           ivar1t0 = ivar1t0 + 1
                        ENDIF
                     END DO loop_t_count
                  ENDIF
                  IF ( ldvar2 ) THEN
                     loop_s_count : DO ij = 1,inpfiles(jj)%nlev
                        IF ( inpfiles(jj)%pdep(ij,ji) >= 6000. ) &
                           & CYCLE
                        IF ( .NOT. BTEST(inpfiles(jj)%ivlqc(ij,ji,2),2) .AND. &
                           & .NOT. BTEST(inpfiles(jj)%idqc(ij,ji),2) ) THEN
                           ivar2t0 = ivar2t0 + 1
                        ENDIF
                     END DO loop_s_count
                  ENDIF
                  loop_p_count : DO ij = 1,inpfiles(jj)%nlev
                     IF ( inpfiles(jj)%pdep(ij,ji) >= 6000. ) &
                        & CYCLE
                     IF ( ( .NOT. BTEST(inpfiles(jj)%ivlqc(ij,ji,1),2) .AND. &
                        &   .NOT. BTEST(inpfiles(jj)%idqc(ij,ji),2) .AND. &
                        &    ldvar1 ) .OR. &
                        & ( .NOT. BTEST(inpfiles(jj)%ivlqc(ij,ji,2),2) .AND. &
                        &   .NOT. BTEST(inpfiles(jj)%idqc(ij,ji),2) .AND. &
                        &     ldvar2 ) ) THEN
                        ip3dt = ip3dt + 1
                        llvalprof = .TRUE.
                     ENDIF
                  END DO loop_p_count

                  IF ( llvalprof ) iprof = iprof + 1

               ENDIF
            END DO

         ENDIF

      END DO prof_files

      !-----------------------------------------------------------------------
      ! Get the time ordered indices of the input data
      !-----------------------------------------------------------------------

      !---------------------------------------------------------------------
      !  Loop over input data files to count total number of profiles
      !---------------------------------------------------------------------
      iproftot = 0
      DO jj = 1, inobf
         DO ji = 1, inpfiles(jj)%nobs
            IF ( BTEST(inpfiles(jj)%ioqc(ji),2 ) ) CYCLE
            IF ( BTEST(inpfiles(jj)%ivqc(ji,1),2) .AND. &
               & BTEST(inpfiles(jj)%ivqc(ji,2),2) ) CYCLE
            IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND. &
               & ( inpfiles(jj)%ptim(ji) <= djulend(jj) )       ) THEN
               iproftot = iproftot + 1
            ENDIF
         END DO
      END DO

      ALLOCATE( iindx(iproftot), ifileidx(iproftot), &
         &      iprofidx(iproftot), zdat(iproftot) )
      jk = 0
      DO jj = 1, inobf
         DO ji = 1, inpfiles(jj)%nobs
            IF ( BTEST(inpfiles(jj)%ioqc(ji),2 ) ) CYCLE
            IF ( BTEST(inpfiles(jj)%ivqc(ji,1),2) .AND. &
               & BTEST(inpfiles(jj)%ivqc(ji,2),2) ) CYCLE
            IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND. &
               & ( inpfiles(jj)%ptim(ji) <= djulend(jj) )       ) THEN
               jk = jk + 1
               ifileidx(jk) = jj
               iprofidx(jk) = ji
               zdat(jk)     = inpfiles(jj)%ptim(ji)
            ENDIF
         END DO
      END DO
      CALL sort_dp_indx( iproftot, &
         &               zdat,     &
         &               iindx   )

      iv3dt(:) = -1
      IF (ldsatt) THEN
         iv3dt(1) = ip3dt
         iv3dt(2) = ip3dt
      ELSE
         iv3dt(1) = ivar1t0
         iv3dt(2) = ivar2t0
      ENDIF
      CALL obs_prof_alloc( profdata, kvars, kextr, iprof, iv3dt, &
         &                 kstp, jpi, jpj, jpk )

      ! * Read obs/positions, QC, all variable and assign to profdata

      profdata%nprof     = 0
      profdata%nvprot(:) = 0
      profdata%cvars(:)  = clvars(:)
      iprof = 0

      ip3dt = 0
      ivar1t = 0
      ivar2t = 0
      itypvar1   (:) = 0
      itypvar1mpp(:) = 0

      itypvar2   (:) = 0
      itypvar2mpp(:) = 0

      ioserrcount = 0
      DO jk = 1, iproftot

         jj = ifileidx(iindx(jk))
         ji = iprofidx(iindx(jk))

            IF ( BTEST(inpfiles(jj)%ioqc(ji),2 ) ) CYCLE
            IF ( BTEST(inpfiles(jj)%ivqc(ji,1),2) .AND. &
               & BTEST(inpfiles(jj)%ivqc(ji,2),2) ) CYCLE

         IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND.  &
            & ( inpfiles(jj)%ptim(ji) <= djulend(jj) ) ) THEN

            IF ( nproc == 0 ) THEN
               IF ( inpfiles(jj)%iproc(ji,1) >  nproc ) CYCLE
            ELSE
               IF ( inpfiles(jj)%iproc(ji,1) /= nproc ) CYCLE
            ENDIF

            llvalprof = .FALSE.

            IF ( inpfiles(jj)%ioqc(ji) > 2 ) CYCLE

            IF ( BTEST(inpfiles(jj)%ioqc(ji),2 ) ) CYCLE
            IF ( BTEST(inpfiles(jj)%ivqc(ji,1),2) .AND. &
               & BTEST(inpfiles(jj)%ivqc(ji,2),2) ) CYCLE

            loop_prof : DO ij = 1, inpfiles(jj)%nlev

               IF ( inpfiles(jj)%pdep(ij,ji) >= 6000. ) &
                  & CYCLE

               IF ( .NOT. BTEST(inpfiles(jj)%ivlqc(ij,ji,1),2) .AND. &
                  & .NOT. BTEST(inpfiles(jj)%idqc(ij,ji),2) ) THEN

                  llvalprof = .TRUE. 
                  EXIT loop_prof

               ENDIF

               IF ( .NOT. BTEST(inpfiles(jj)%ivlqc(ij,ji,2),2) .AND. &
                  & .NOT. BTEST(inpfiles(jj)%idqc(ij,ji),2) ) THEN

                  llvalprof = .TRUE. 
                  EXIT loop_prof

               ENDIF

            END DO loop_prof

            ! Set profile information

            IF ( llvalprof ) THEN

               iprof = iprof + 1

               CALL jul2greg( isec,                   &
                  &           imin,                   &
                  &           ihou,                   &
                  &           iday,                   &
                  &           imon,                   &
                  &           iyea,                   &
                  &           inpfiles(jj)%ptim(ji), &
                  &           irefdate(jj) )


               ! Profile time coordinates
               profdata%nyea(iprof) = iyea
               profdata%nmon(iprof) = imon
               profdata%nday(iprof) = iday
               profdata%nhou(iprof) = ihou
               profdata%nmin(iprof) = imin

               ! Profile space coordinates
               profdata%rlam(iprof) = inpfiles(jj)%plam(ji)
               profdata%rphi(iprof) = inpfiles(jj)%pphi(ji)

               ! Coordinate search parameters
               profdata%mi  (iprof,1) = inpfiles(jj)%iobsi(ji,1)
               profdata%mj  (iprof,1) = inpfiles(jj)%iobsj(ji,1)
               profdata%mi  (iprof,2) = inpfiles(jj)%iobsi(ji,2)
               profdata%mj  (iprof,2) = inpfiles(jj)%iobsj(ji,2)

               ! Profile WMO number
               profdata%cwmo(iprof) = inpfiles(jj)%cdwmo(ji)

               ! Instrument type
               READ( inpfiles(jj)%cdtyp(ji), '(I4)', IOSTAT = ios, ERR = 901 ) itype
901            IF ( ios /= 0 ) THEN
                  IF (ioserrcount == 0) CALL ctl_warn ( 'Problem converting an instrument type to integer. Setting type to zero' )
                  ioserrcount = ioserrcount + 1
                  itype = 0
               ENDIF

               profdata%ntyp(iprof) = itype

               ! QC stuff

               profdata%nqc(iprof)     = inpfiles(jj)%ioqc(ji)
               profdata%nqcf(:,iprof)  = inpfiles(jj)%ioqcf(:,ji)
               profdata%ipqc(iprof)    = inpfiles(jj)%ipqc(ji)
               profdata%ipqcf(:,iprof) = inpfiles(jj)%ipqcf(:,ji)
               profdata%itqc(iprof)    = inpfiles(jj)%itqc(ji)
               profdata%itqcf(:,iprof) = inpfiles(jj)%itqcf(:,ji)
               profdata%ivqc(iprof,:)  = inpfiles(jj)%ivqc(ji,:)
               profdata%ivqcf(:,iprof,:) = inpfiles(jj)%ivqcf(:,ji,:)

               ! Bookkeeping data to match profiles
               profdata%npidx(iprof) = iprof
               profdata%npfil(iprof) = iindx(jk)

               ! Observation QC flag (whole profile)
               profdata%nqc(iprof)  = 0 !TODO

               loop_p : DO ij = 1, inpfiles(jj)%nlev

                  IF ( inpfiles(jj)%pdep(ij,ji) >= 6000. ) &
                     & CYCLE

                  IF (ldsatt) THEN

                     IF ( ( .NOT. BTEST(inpfiles(jj)%ivlqc(ij,ji,1),2) .AND. &
                        &   .NOT. BTEST(inpfiles(jj)%idqc(ij,ji),2) .AND. &
                        &    ldvar1 ) .OR. &
                        & ( .NOT. BTEST(inpfiles(jj)%ivlqc(ij,ji,2),2) .AND. &
                        &   .NOT. BTEST(inpfiles(jj)%idqc(ij,ji),2) .AND. &
                        &   ldvar2 ) ) THEN
                        ip3dt = ip3dt + 1
                     ELSE
                        CYCLE
                     ENDIF

                  ENDIF

                  IF ( ( .NOT. BTEST(inpfiles(jj)%ivlqc(ij,ji,1),2) .AND. &
                    &   .NOT. BTEST(inpfiles(jj)%idqc(ij,ji),2) .AND. &
                    &    ldvar1 ) .OR. ldsatt ) THEN

                     IF (ldsatt) THEN

                        ivar1t = ip3dt

                     ELSE

                        ivar1t = ivar1t + 1

                     ENDIF

                     ! Depth of var1 observation
                     profdata%var(1)%vdep(ivar1t) = &
                        &                inpfiles(jj)%pdep(ij,ji)

                     ! Depth of var1 observation QC
                     profdata%var(1)%idqc(ivar1t) = &
                        &                inpfiles(jj)%idqc(ij,ji)

                     ! Depth of var1 observation QC flags
                     profdata%var(1)%idqcf(:,ivar1t) = &
                        &                inpfiles(jj)%idqcf(:,ij,ji)

                     ! Profile index
                     profdata%var(1)%nvpidx(ivar1t) = iprof

                     ! Vertical index in original profile
                     profdata%var(1)%nvlidx(ivar1t) = ij

                     ! Profile var1 value
                     IF ( .NOT. BTEST(inpfiles(jj)%ivlqc(ij,ji,1),2) .AND. &
                        & .NOT. BTEST(inpfiles(jj)%idqc(ij,ji),2) ) THEN
                        profdata%var(1)%vobs(ivar1t) = &
                           &                inpfiles(jj)%pob(ij,ji,1)
                        IF ( ldmod ) THEN
                           profdata%var(1)%vmod(ivar1t) = &
                              &                inpfiles(jj)%padd(ij,ji,1,1)
                        ENDIF
                        ! Count number of profile var1 data as function of type
                        itypvar1( profdata%ntyp(iprof) + 1 ) = &
                           & itypvar1( profdata%ntyp(iprof) + 1 ) + 1
                     ELSE
                        profdata%var(1)%vobs(ivar1t) = fbrmdi
                     ENDIF

                     ! Profile var1 qc
                     profdata%var(1)%nvqc(ivar1t) = &
                        & inpfiles(jj)%ivlqc(ij,ji,1)

                     ! Profile var1 qc flags
                     profdata%var(1)%nvqcf(:,ivar1t) = &
                        & inpfiles(jj)%ivlqcf(:,ij,ji,1)

                     ! Profile insitu T value
                     IF ( TRIM( inpfiles(jj)%cname(1) ) == 'POTM' ) THEN
                        profdata%var(1)%vext(ivar1t,1) = &
                           &                inpfiles(jj)%pext(ij,ji,1)
                     ENDIF

                  ENDIF

                  IF ( ( .NOT. BTEST(inpfiles(jj)%ivlqc(ij,ji,2),2) .AND. &
                     &   .NOT. BTEST(inpfiles(jj)%idqc(ij,ji),2)    .AND. &
                     &   ldvar2 ) .OR. ldsatt ) THEN

                     IF (ldsatt) THEN

                        ivar2t = ip3dt

                     ELSE

                        ivar2t = ivar2t + 1

                     ENDIF

                     ! Depth of var2 observation
                     profdata%var(2)%vdep(ivar2t) = &
                        &                inpfiles(jj)%pdep(ij,ji)

                     ! Depth of var2 observation QC
                     profdata%var(2)%idqc(ivar2t) = &
                        &                inpfiles(jj)%idqc(ij,ji)

                     ! Depth of var2 observation QC flags
                     profdata%var(2)%idqcf(:,ivar2t) = &
                        &                inpfiles(jj)%idqcf(:,ij,ji)

                     ! Profile index
                     profdata%var(2)%nvpidx(ivar2t) = iprof

                     ! Vertical index in original profile
                     profdata%var(2)%nvlidx(ivar2t) = ij

                     ! Profile var2 value
                  IF (  ( .NOT. BTEST(inpfiles(jj)%ivlqc(ij,ji,2),2) ) .AND. &
                    &   ( .NOT. BTEST(inpfiles(jj)%idqc(ij,ji),2)    )  ) THEN
                        profdata%var(2)%vobs(ivar2t) = &
                           &                inpfiles(jj)%pob(ij,ji,2)
                        IF ( ldmod ) THEN
                           profdata%var(2)%vmod(ivar2t) = &
                              &                inpfiles(jj)%padd(ij,ji,1,2)
                        ENDIF
                        ! Count number of profile var2 data as function of type
                        itypvar2( profdata%ntyp(iprof) + 1 ) = &
                           & itypvar2( profdata%ntyp(iprof) + 1 ) + 1
                     ELSE
                        profdata%var(2)%vobs(ivar2t) = fbrmdi
                     ENDIF

                     ! Profile var2 qc
                     profdata%var(2)%nvqc(ivar2t) = &
                        & inpfiles(jj)%ivlqc(ij,ji,2)

                     ! Profile var2 qc flags
                     profdata%var(2)%nvqcf(:,ivar2t) = &
                        & inpfiles(jj)%ivlqcf(:,ij,ji,2)

                  ENDIF

               END DO loop_p

            ENDIF

         ENDIF

      END DO

      !-----------------------------------------------------------------------
      ! Sum up over processors
      !-----------------------------------------------------------------------

      CALL obs_mpp_sum_integer ( ivar1t0, ivar1tmpp )
      CALL obs_mpp_sum_integer ( ivar2t0, ivar2tmpp )
      CALL obs_mpp_sum_integer ( ip3dt,   ip3dtmpp  )

      CALL obs_mpp_sum_integers( itypvar1, itypvar1mpp, ntyp1770 + 1 )
      CALL obs_mpp_sum_integers( itypvar2, itypvar2mpp, ntyp1770 + 1 )

      !-----------------------------------------------------------------------
      ! Output number of observations.
      !-----------------------------------------------------------------------
      IF(lwp) THEN
         WRITE(numout,*) 
         WRITE(numout,'(A)') ' Profile data'
         WRITE(numout,'(1X,A)') '------------'
         WRITE(numout,*) 
         WRITE(numout,'(1X,A)') 'Profile data, '//TRIM( profdata%cvars(1) )
         WRITE(numout,'(1X,A)') '------------------------'
         DO ji = 0, ntyp1770
            IF ( itypvar1mpp(ji+1) > 0 ) THEN
               WRITE(numout,'(1X,A3,1X,A48,A3,I8)') ctypshort(ji), &
                  & cwmonam1770(ji)(1:52),' = ', &
                  & itypvar1mpp(ji+1)
            ENDIF
         END DO
         WRITE(numout,'(1X,A)') &
            & '---------------------------------------------------------------'
         WRITE(numout,'(1X,A55,I8)') &
            & 'Total profile data for variable '//TRIM( profdata%cvars(1) )// &
            & '             = ', ivar1tmpp
         WRITE(numout,'(1X,A)') &
            & '---------------------------------------------------------------'
         WRITE(numout,*) 
         WRITE(numout,'(1X,A)') 'Profile data, '//TRIM( profdata%cvars(2) )
         WRITE(numout,'(1X,A)') '------------------------'
         DO ji = 0, ntyp1770
            IF ( itypvar2mpp(ji+1) > 0 ) THEN
               WRITE(numout,'(1X,A3,1X,A48,A3,I8)') ctypshort(ji), &
                  & cwmonam1770(ji)(1:52),' = ', &
                  & itypvar2mpp(ji+1)
            ENDIF
         END DO
         WRITE(numout,'(1X,A)') &
            & '---------------------------------------------------------------'
         WRITE(numout,'(1X,A55,I8)') &
            & 'Total profile data for variable '//TRIM( profdata%cvars(2) )// &
            & '             = ', ivar2tmpp
         WRITE(numout,'(1X,A)') &
            & '---------------------------------------------------------------'
         WRITE(numout,*) 
      ENDIF

      IF (ldsatt) THEN
         profdata%nvprot(1)    = ip3dt
         profdata%nvprot(2)    = ip3dt
         profdata%nvprotmpp(1) = ip3dtmpp
         profdata%nvprotmpp(2) = ip3dtmpp
      ELSE
         profdata%nvprot(1)    = ivar1t
         profdata%nvprot(2)    = ivar2t
         profdata%nvprotmpp(1) = ivar1tmpp
         profdata%nvprotmpp(2) = ivar2tmpp
      ENDIF
      profdata%nprof        = iprof

      !-----------------------------------------------------------------------
      ! Model level search
      !-----------------------------------------------------------------------
      IF ( ldvar1 ) THEN
         CALL obs_level_search( jpk, gdept_1d, &
            & profdata%nvprot(1), profdata%var(1)%vdep, &
            & profdata%var(1)%mvk )
      ENDIF
      IF ( ldvar2 ) THEN
         CALL obs_level_search( jpk, gdept_1d, &
            & profdata%nvprot(2), profdata%var(2)%vdep, &
            & profdata%var(2)%mvk )
      ENDIF

      !-----------------------------------------------------------------------
      ! Set model equivalent to 99999
      !-----------------------------------------------------------------------
      IF ( .NOT. ldmod ) THEN
         DO jvar = 1, kvars
            profdata%var(jvar)%vmod(:) = fbrmdi
         END DO
      ENDIF
      !-----------------------------------------------------------------------
      ! Deallocate temporary data
      !-----------------------------------------------------------------------
      DEALLOCATE( ifileidx, iprofidx, zdat, clvars )

      !-----------------------------------------------------------------------
      ! Deallocate input data
      !-----------------------------------------------------------------------
      DO jj = 1, inobf
         CALL dealloc_obfbdata( inpfiles(jj) )
      END DO
      DEALLOCATE( inpfiles )

   END SUBROUTINE obs_rea_prof

END MODULE obs_read_prof
