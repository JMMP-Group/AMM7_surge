MODULE obs_read_surf
   !!======================================================================
   !!                       ***  MODULE obs_read_surf  ***
   !! Observation diagnostics: Read the surface data from feedback files
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   obs_rea_surf : Driver for reading surface data from feedback files
   !!----------------------------------------------------------------------

   !! * Modules used
   USE par_kind                 ! Precision variables
   USE in_out_manager           ! I/O manager
   USE dom_oce                  ! Ocean space and time domain variables
   USE obs_mpp                  ! MPP support routines for observation diagnostics
   USE julian                   ! Julian date routines
   USE obs_utils                ! Observation operator utility functions
   USE obs_grid                 ! Grid search
   USE obs_sort                 ! Sorting observation arrays
   USE obs_surf_def             ! Surface observation definitions
   USE obs_types                ! Observation type definitions
   USE obs_fbm                  ! Feedback routines
   USE netcdf                   ! NetCDF library

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE

   PUBLIC obs_rea_surf      ! Read the surface observations from the point data

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE obs_rea_surf( surfdata, knumfiles, cdfilenames, &
      &                     kvars, kextr, kstp, ddobsini, ddobsend, &
      &                     ldignmis, ldmod, ldnightav )
      !!---------------------------------------------------------------------
      !!
      !!                   *** ROUTINE obs_rea_surf ***
      !!
      !! ** Purpose : Read from file the surface data
      !!
      !! ** Method  : Read in the data from feedback format files and 
      !!              put into the NEMO internal surface data structure
      !!
      !! ** Action  : 
      !!
      !!
      !! History :  
      !!      ! :  2009-01 (K. Mogensen) Initial version based on old versions
      !!      ! :  2015-02 (M. Martin)   Unify the different surface data type reading.
      !!----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) :: &
         & surfdata                     ! Surface data to be read
      INTEGER, INTENT(IN) :: knumfiles  ! Number of corio format files to read
      CHARACTER(LEN=128), INTENT(IN) :: &
         & cdfilenames(knumfiles)       ! File names to read in
      INTEGER, INTENT(IN) :: kvars      ! Number of variables in surfdata
      INTEGER, INTENT(IN) :: kextr      ! Number of extra fields for each var
      INTEGER, INTENT(IN) :: kstp       ! Ocean time-step index
      LOGICAL, INTENT(IN) :: ldignmis   ! Ignore missing files
      LOGICAL, INTENT(IN) :: ldmod      ! Initialize model from input data
      LOGICAL, INTENT(IN) :: ldnightav  ! Observations represent a night-time average
      REAL(dp), INTENT(IN) :: ddobsini   ! Obs. ini time in YYYYMMDD.HHMMSS
      REAL(dp), INTENT(IN) :: ddobsend   ! Obs. end time in YYYYMMDD.HHMMSS

      !! * Local declarations
      CHARACTER(LEN=11), PARAMETER :: cpname='obs_rea_surf'
      CHARACTER(len=8) :: clrefdate
      CHARACTER(len=8), DIMENSION(:), ALLOCATABLE :: clvars
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jk
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
      INTEGER :: itype
      INTEGER :: iobsmpp
      INTEGER :: iobs
      INTEGER :: iobstot
      INTEGER :: ios
      INTEGER :: ioserrcount
      INTEGER, PARAMETER :: jpsurfmaxtype = 1024
      INTEGER, DIMENSION(knumfiles) :: irefdate
      INTEGER, DIMENSION(jpsurfmaxtype+1) :: &
         & ityp, &
         & itypmpp
      INTEGER, DIMENSION(:), ALLOCATABLE :: &
         & iobsi,    &
         & iobsj,    &
         & iproc,    &
         & iindx,    &
         & ifileidx, &
         & isurfidx
      REAL(wp), DIMENSION(:), ALLOCATABLE :: &
         & zphi, &
         & zlam
      REAL(wp), DIMENSION(:), ALLOCATABLE :: &
         & zdat
      REAL(wp), DIMENSION(knumfiles) :: &
         & djulini, &
         & djulend
      LOGICAL :: llvalprof
      TYPE(obfbdata), POINTER, DIMENSION(:) :: &
         & inpfiles

      ! Local initialization
      iobs = 0

      !-----------------------------------------------------------------------
      ! Count the number of files needed and allocate the obfbdata type
      !-----------------------------------------------------------------------

      inobf = knumfiles

      ALLOCATE( inpfiles(inobf) )

      surf_files : DO jj = 1, inobf

         !---------------------------------------------------------------------
         ! Prints
         !---------------------------------------------------------------------
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) ' obs_rea_surf : Reading from file = ', &
               & TRIM( TRIM( cdfilenames(jj) ) )
            WRITE(numout,*) ' ~~~~~~~~~~~'
            WRITE(numout,*)
         ENDIF

         !---------------------------------------------------------------------
         !  Initialization: Open file and get dimensions only
         !---------------------------------------------------------------------

         iflag = nf90_open( TRIM( TRIM( cdfilenames(jj) ) ), nf90_nowrite, &
            &                      i_file_id )

         IF ( iflag /= nf90_noerr ) THEN

            IF ( ldignmis ) THEN
               inpfiles(jj)%nobs = 0
               CALL ctl_warn( 'File ' // TRIM( TRIM( cdfilenames(jj) ) ) // &
                  &           ' not found' )
            ELSE 
               CALL ctl_stop( 'File ' // TRIM( TRIM( cdfilenames(jj) ) ) // &
                  &           ' not found' )
            ENDIF

         ELSE 

            !------------------------------------------------------------------
            !  Close the file since it is opened in read_obfbdata
            !------------------------------------------------------------------

            iflag = nf90_close( i_file_id )

            !------------------------------------------------------------------
            !  Read the surface file into inpfiles
            !------------------------------------------------------------------
            CALL init_obfbdata( inpfiles(jj) )
            CALL read_obfbdata( TRIM( cdfilenames(jj) ), inpfiles(jj), &
               &                ldgrid = .TRUE. )

            IF ( ldmod .AND. ( inpfiles(jj)%nadd == 0 ) ) THEN
               CALL ctl_stop( 'Model not in input data' )
               RETURN
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

            IF (lwp) WRITE(numout,*)'Observation file contains ',inpfiles(jj)%nobs,' observations'

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

            IF ( ldnightav ) THEN

               IF ( lwp ) THEN
                  WRITE(numout,*)'Resetting time of night-time averaged observations', &
                     &             ' to the end of the day'
               ENDIF

               DO ji = 1, inpfiles(jj)%nobs
                  !  for night-time averaged data force the time
                  !  to be the last time-step of the day, but still within the day.
                  IF ( inpfiles(jj)%ptim(ji) >= 0. ) THEN
                     inpfiles(jj)%ptim(ji) = &
                        & INT(inpfiles(jj)%ptim(ji)) + 0.9999
                  ELSE
                     inpfiles(jj)%ptim(ji) = &
                        & INT(inpfiles(jj)%ptim(ji)) - 0.0001
                  ENDIF
               END DO
            ENDIF

            IF ( inpfiles(jj)%nobs > 0 ) THEN
               inpfiles(jj)%iproc = -1
               inpfiles(jj)%iobsi = -1
               inpfiles(jj)%iobsj = -1
            ENDIF
            inowin = 0
            DO ji = 1, inpfiles(jj)%nobs
               IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND. &
                  & ( inpfiles(jj)%ptim(ji) <= djulend(jj) )       ) THEN
                  inowin = inowin + 1
               ENDIF
            END DO
            ALLOCATE( zlam(inowin)  )
            ALLOCATE( zphi(inowin)  )
            ALLOCATE( iobsi(inowin) )
            ALLOCATE( iobsj(inowin) )
            ALLOCATE( iproc(inowin) )
            inowin = 0
            DO ji = 1, inpfiles(jj)%nobs
               IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND. &
                  & ( inpfiles(jj)%ptim(ji) <= djulend(jj) )       ) THEN
                  inowin = inowin + 1
                  zlam(inowin) = inpfiles(jj)%plam(ji)
                  zphi(inowin) = inpfiles(jj)%pphi(ji)
               ENDIF
            END DO

            CALL obs_grid_search( inowin, zlam, zphi, iobsi, iobsj, iproc, 'T' )

            inowin = 0
            DO ji = 1, inpfiles(jj)%nobs
               IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND. &
                  & ( inpfiles(jj)%ptim(ji) <= djulend(jj) )       ) THEN
                  inowin = inowin + 1
                  inpfiles(jj)%iproc(ji,1) = iproc(inowin)
                  inpfiles(jj)%iobsi(ji,1) = iobsi(inowin)
                  inpfiles(jj)%iobsj(ji,1) = iobsj(inowin)
               ENDIF
            END DO
            DEALLOCATE( zlam, zphi, iobsi, iobsj, iproc )

            DO ji = 1, inpfiles(jj)%nobs
               IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND. &
                  & ( inpfiles(jj)%ptim(ji) <= djulend(jj) )       ) THEN
                  IF ( nproc == 0 ) THEN
                     IF ( inpfiles(jj)%iproc(ji,1) >  nproc ) CYCLE
                  ELSE
                     IF ( inpfiles(jj)%iproc(ji,1) /= nproc ) CYCLE
                  ENDIF
                  llvalprof = .FALSE.
                  IF ( .NOT. BTEST(inpfiles(jj)%ivlqc(1,ji,1),2) ) THEN
                     iobs = iobs + 1
                  ENDIF
               ENDIF
            END DO

         ENDIF

      END DO surf_files

      !-----------------------------------------------------------------------
      ! Get the time ordered indices of the input data
      !-----------------------------------------------------------------------

      !---------------------------------------------------------------------
      !  Loop over input data files to count total number of profiles
      !---------------------------------------------------------------------
      iobstot = 0
      DO jj = 1, inobf
         DO ji = 1, inpfiles(jj)%nobs
            IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND. &
               & ( inpfiles(jj)%ptim(ji) <= djulend(jj) )       ) THEN
               iobstot = iobstot + 1
            ENDIF
         END DO
      END DO

      ALLOCATE( iindx(iobstot), ifileidx(iobstot), &
         &      isurfidx(iobstot), zdat(iobstot) )
      jk = 0
      DO jj = 1, inobf
         DO ji = 1, inpfiles(jj)%nobs
            IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND. &
               & ( inpfiles(jj)%ptim(ji) <= djulend(jj) )       ) THEN
               jk = jk + 1
               ifileidx(jk) = jj
               isurfidx(jk) = ji
               zdat(jk)     = inpfiles(jj)%ptim(ji)
            ENDIF
         END DO
      END DO
      CALL sort_dp_indx( iobstot, &
         &               zdat,     &
         &               iindx   )

      CALL obs_surf_alloc( surfdata, iobs, kvars, kextr, kstp, jpi, jpj )

      ! Read obs/positions, QC, all variable and assign to surfdata

      iobs = 0

      surfdata%cvars(:)  = clvars(:)

      ityp   (:) = 0
      itypmpp(:) = 0

      ioserrcount = 0

      DO jk = 1, iobstot

         jj = ifileidx(iindx(jk))
         ji = isurfidx(iindx(jk))
         IF ( ( inpfiles(jj)%ptim(ji) >  djulini(jj) ) .AND.  &
            & ( inpfiles(jj)%ptim(ji) <= djulend(jj) ) ) THEN

            IF ( nproc == 0 ) THEN
               IF ( inpfiles(jj)%iproc(ji,1) >  nproc ) CYCLE
            ELSE
               IF ( inpfiles(jj)%iproc(ji,1) /= nproc ) CYCLE
            ENDIF

            ! Set observation information

            IF ( .NOT. BTEST(inpfiles(jj)%ivlqc(1,ji,1),2) ) THEN

               iobs = iobs + 1

               CALL jul2greg( isec,                   &
                  &           imin,                   &
                  &           ihou,                   &
                  &           iday,                   &
                  &           imon,                   &
                  &           iyea,                   &
                  &           inpfiles(jj)%ptim(ji), &
                  &           irefdate(jj) )


               ! Surface time coordinates
               surfdata%nyea(iobs) = iyea
               surfdata%nmon(iobs) = imon
               surfdata%nday(iobs) = iday
               surfdata%nhou(iobs) = ihou
               surfdata%nmin(iobs) = imin

               ! Surface space coordinates
               surfdata%rlam(iobs) = inpfiles(jj)%plam(ji)
               surfdata%rphi(iobs) = inpfiles(jj)%pphi(ji)

               ! Coordinate search parameters
               surfdata%mi  (iobs) = inpfiles(jj)%iobsi(ji,1)
               surfdata%mj  (iobs) = inpfiles(jj)%iobsj(ji,1)

               ! WMO number
               surfdata%cwmo(iobs) = inpfiles(jj)%cdwmo(ji)

               ! Instrument type
               READ( inpfiles(jj)%cdtyp(ji), '(I4)', IOSTAT = ios, ERR = 901 ) itype
901            IF ( ios /= 0 ) THEN
                  IF (ioserrcount == 0) THEN
                     CALL ctl_warn ( 'Problem converting an instrument type ', &
                        &            'to integer. Setting type to zero' )
                  ENDIF
                  ioserrcount = ioserrcount + 1
                  itype = 0
               ENDIF
               surfdata%ntyp(iobs) = itype
               IF ( itype < jpsurfmaxtype + 1 ) THEN
                  ityp(itype+1) = ityp(itype+1) + 1
               ELSE
                  IF(lwp)WRITE(numout,*)'WARNING:Increase jpsurfmaxtype in ',&
                     &                  cpname
               ENDIF

               ! Bookkeeping data to match observations
               surfdata%nsidx(iobs) = iobs
               surfdata%nsfil(iobs) = iindx(jk)

               ! QC flags
               surfdata%nqc(iobs) = inpfiles(jj)%ivqc(ji,1)

               ! Observed value
               surfdata%robs(iobs,1) = inpfiles(jj)%pob(1,ji,1)


               ! Model and MDT is set to fbrmdi unless read from file
               IF ( ldmod ) THEN
                  surfdata%rmod(iobs,1) = inpfiles(jj)%padd(1,ji,1,1)
                  IF ( TRIM(surfdata%cvars(1)) == 'SLA' ) THEN
                     surfdata%rext(iobs,1) = inpfiles(jj)%padd(1,ji,2,1)
                     surfdata%rext(iobs,2) = inpfiles(jj)%pext(1,ji,1)
                  ENDIF
                ELSE
                  surfdata%rmod(iobs,1) = fbrmdi
                  IF ( TRIM(surfdata%cvars(1)) == 'SLA' ) surfdata%rext(iobs,:) = fbrmdi
               ENDIF
            ENDIF
         ENDIF

      END DO

      !-----------------------------------------------------------------------
      ! Sum up over processors
      !-----------------------------------------------------------------------

      CALL obs_mpp_sum_integer( iobs, iobsmpp )
      CALL obs_mpp_sum_integers( ityp, itypmpp, jpsurfmaxtype + 1 )

      !-----------------------------------------------------------------------
      ! Output number of observations.
      !-----------------------------------------------------------------------
      IF (lwp) THEN

         WRITE(numout,*)
         WRITE(numout,'(1X,A)')TRIM( surfdata%cvars(1) )//' data'
         WRITE(numout,'(1X,A)')'--------------'
         DO jj = 1,8
            IF ( itypmpp(jj) > 0 ) THEN
               WRITE(numout,'(1X,A4,I4,A3,I10)')'Type ', jj,' = ',itypmpp(jj)
            ENDIF
         END DO
         WRITE(numout,'(1X,A)') &
            & '---------------------------------------------------------------'
         WRITE(numout,'(1X,A,I8)') &
            & 'Total data for variable '//TRIM( surfdata%cvars(1) )// &
            & '           = ', iobsmpp
         WRITE(numout,'(1X,A)') &
            & '---------------------------------------------------------------'
         WRITE(numout,*)

      ENDIF

      !-----------------------------------------------------------------------
      ! Deallocate temporary data
      !-----------------------------------------------------------------------
      DEALLOCATE( ifileidx, isurfidx, zdat, clvars )

      !-----------------------------------------------------------------------
      ! Deallocate input data
      !-----------------------------------------------------------------------
      DO jj = 1, inobf
         IF ( inpfiles(jj)%lalloc ) THEN
            CALL dealloc_obfbdata( inpfiles(jj) )
         ENDIF
      END DO
      DEALLOCATE( inpfiles )

   END SUBROUTINE obs_rea_surf

END MODULE obs_read_surf
