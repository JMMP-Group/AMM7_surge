MODULE fldread
   !!======================================================================
   !!                       ***  MODULE  fldread  ***
   !! Ocean forcing:  read input field for surface boundary condition
   !!=====================================================================
   !! History :  2.0  !  2006-06  (S. Masson, G. Madec)  Original code
   !!            3.0  !  2008-05  (S. Alderson)  Modified for Interpolation in memory from input grid to model grid
   !!            3.4  !  2013-10  (D. Delrosso, P. Oddo)  suppression of land point prior to interpolation
   !!                 !  12-2015  (J. Harle) Adding BDY on-the-fly interpolation
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   fld_read      : read input fields used for the computation of the surface boundary condition
   !!   fld_init      : initialization of field read
   !!   fld_rec       : determined the record(s) to be read
   !!   fld_get       : read the data
   !!   fld_map       : read global data from file and map onto local data using a general mapping (use for open boundaries)
   !!   fld_rot       : rotate the vector fields onto the local grid direction
   !!   fld_clopn     : update the data file name and close/open the files
   !!   fld_fill      : fill the data structure with the associated information read in namelist
   !!   wgt_list      : manage the weights used for interpolation
   !!   wgt_print     : print the list of known weights
   !!   fld_weight    : create a WGT structure and fill in data from file, restructuring as required
   !!   apply_seaoverland : fill land with ocean values
   !!   seaoverland   : create shifted matrices for seaoverland application
   !!   fld_interp    : apply weights to input gridded data to create data on model grid
   !!   ksec_week     : function returning the first 3 letters of the first day of the weekly file
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constant
   USE sbc_oce        ! surface boundary conditions : fields
   USE geo2ocean      ! for vector rotation on to model grid
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE ioipsl  , ONLY : ymds2ju, ju2ymds   ! for calendar
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (C1D case)
   
   IMPLICIT NONE
   PRIVATE   
 
   PUBLIC   fld_map    ! routine called by tides_init
   PUBLIC   fld_read, fld_fill   ! called by sbc... modules
   PUBLIC   fld_clopn

   TYPE, PUBLIC ::   FLD_N      !: Namelist field informations
      CHARACTER(len = 256) ::   clname      ! generic name of the NetCDF flux file
      REAL(wp)             ::   freqh       ! frequency of each flux file
      CHARACTER(len = 34)  ::   clvar       ! generic name of the variable in the NetCDF flux file
      LOGICAL              ::   ln_tint     ! time interpolation or not (T/F)
      LOGICAL              ::   ln_clim     ! climatology or not (T/F)
      CHARACTER(len = 8)   ::   cltype      ! type of data file 'daily', 'monthly' or yearly'
      CHARACTER(len = 256) ::   wname       ! generic name of a NetCDF weights file to be used, blank if not
      CHARACTER(len = 34)  ::   vcomp       ! symbolic component name if a vector that needs rotation
      !                                     ! a string starting with "U" or "V" for each component   
      !                                     ! chars 2 onwards identify which components go together  
      CHARACTER(len = 34)  ::   lname       ! generic name of a NetCDF land/sea mask file to be used, blank if not 
      !                                     ! 0=sea 1=land
   END TYPE FLD_N

   TYPE, PUBLIC ::   FLD        !: Input field related variables
      CHARACTER(len = 256)            ::   clrootname   ! generic name of the NetCDF file
      CHARACTER(len = 256)            ::   clname       ! current name of the NetCDF file
      REAL(wp)                        ::   freqh        ! frequency of each flux file
      CHARACTER(len = 34)             ::   clvar        ! generic name of the variable in the NetCDF flux file
      LOGICAL                         ::   ln_tint      ! time interpolation or not (T/F)
      LOGICAL                         ::   ln_clim      ! climatology or not (T/F)
      CHARACTER(len = 8)              ::   cltype       ! type of data file 'daily', 'monthly' or yearly'
      INTEGER                         ::   num          ! iom id of the jpfld files to be read
      INTEGER , DIMENSION(2)          ::   nrec_b       ! before record (1: index, 2: second since Jan. 1st 00h of nit000 year)
      INTEGER , DIMENSION(2)          ::   nrec_a       ! after  record (1: index, 2: second since Jan. 1st 00h of nit000 year)
      REAL(wp) , ALLOCATABLE, DIMENSION(:,:,:  ) ::   fnow   ! input fields interpolated to now time step
      REAL(wp) , ALLOCATABLE, DIMENSION(:,:,:,:) ::   fdta   ! 2 consecutive record of input fields
      CHARACTER(len = 256)            ::   wgtname      ! current name of the NetCDF weight file acting as a key
      !                                                 ! into the WGTLIST structure
      CHARACTER(len = 34)             ::   vcomp        ! symbolic name for a vector component that needs rotation
      LOGICAL, DIMENSION(2)           ::   rotn         ! flag to indicate whether before/after field has been rotated
      INTEGER                         ::   nreclast     ! last record to be read in the current file
      CHARACTER(len = 256)            ::   lsmname      ! current name of the NetCDF mask file acting as a key
      !                                                 ! 
      !                                                 ! Variables related to BDY
      INTEGER                         ::   igrd         !   grid type for bdy data
      INTEGER                         ::   ibdy         !   bdy set id number
      INTEGER, POINTER, DIMENSION(:)  ::   imap         !   Array of integer pointers to 1D arrays
      LOGICAL                         ::   ltotvel      !   total velocity or not (T/F)
      LOGICAL                         ::   lzint        !   T if it requires a vertical interpolation
   END TYPE FLD

!$AGRIF_DO_NOT_TREAT

   !! keep list of all weights variables so they're only read in once
   !! need to add AGRIF directives not to process this structure
   !! also need to force wgtname to include AGRIF nest number
   TYPE         ::   WGT        !: Input weights related variables
      CHARACTER(len = 256)                    ::   wgtname      ! current name of the NetCDF weight file
      INTEGER , DIMENSION(2)                  ::   ddims        ! shape of input grid
      INTEGER , DIMENSION(2)                  ::   botleft      ! top left corner of box in input grid containing 
      !                                                         ! current processor grid
      INTEGER , DIMENSION(2)                  ::   topright     ! top right corner of box 
      INTEGER                                 ::   jpiwgt       ! width of box on input grid
      INTEGER                                 ::   jpjwgt       ! height of box on input grid
      INTEGER                                 ::   numwgt       ! number of weights (4=bilinear, 16=bicubic)
      INTEGER                                 ::   nestid       ! for agrif, keep track of nest we're in
      INTEGER                                 ::   overlap      ! =0 when cyclic grid has no overlapping EW columns
      !                                                         ! =>1 when they have one or more overlapping columns      
      !                                                         ! =-1 not cyclic
      LOGICAL                                 ::   cyclic       ! east-west cyclic or not
      INTEGER,  DIMENSION(:,:,:), POINTER     ::   data_jpi     ! array of source integers
      INTEGER,  DIMENSION(:,:,:), POINTER     ::   data_jpj     ! array of source integers
      REAL(wp), DIMENSION(:,:,:), POINTER     ::   data_wgt     ! array of weights on model grid
      REAL(wp), DIMENSION(:,:,:), POINTER     ::   fly_dta      ! array of values on input grid
      REAL(wp), DIMENSION(:,:,:), POINTER     ::   col          ! temporary array for reading in columns
   END TYPE WGT

   INTEGER,     PARAMETER             ::   tot_wgts = 20
   TYPE( WGT ), DIMENSION(tot_wgts)   ::   ref_wgts     ! array of wgts
   INTEGER                            ::   nxt_wgt = 1  ! point to next available space in ref_wgts array
   REAL(wp), PARAMETER                ::   undeff_lsm = -999.00_wp

!$AGRIF_END_DO_NOT_TREAT

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE fld_read( kt, kn_fsbc, sd, kit, kt_offset )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_read  ***
      !!                   
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!                (momentum, heat, freshwater and runoff) 
      !!
      !! ** Method  :   READ each input fields in NetCDF files using IOM
      !!      and intepolate it to the model time-step.
      !!         Several assumptions are made on the input file:
      !!      blahblahblah....
      !!----------------------------------------------------------------------
      INTEGER  , INTENT(in   )               ::   kt        ! ocean time step
      INTEGER  , INTENT(in   )               ::   kn_fsbc   ! sbc computation period (in time step) 
      TYPE(FLD), INTENT(inout), DIMENSION(:) ::   sd        ! input field related variables
      INTEGER  , INTENT(in   ), OPTIONAL     ::   kit       ! subcycle timestep for timesplitting option
      INTEGER  , INTENT(in   ), OPTIONAL     ::   kt_offset ! provide fields at time other than "now"
      !                                                     !   kt_offset = -1 => fields at "before" time level
      !                                                     !   kt_offset = +1 => fields at "after"  time level
      !                                                     !   etc.
      !!
      INTEGER  ::   itmp         ! local variable
      INTEGER  ::   imf          ! size of the structure sd
      INTEGER  ::   jf           ! dummy indices
      INTEGER  ::   isecend      ! number of second since Jan. 1st 00h of nit000 year at nitend
      INTEGER  ::   isecsbc      ! number of seconds between Jan. 1st 00h of nit000 year and the middle of sbc time step
      INTEGER  ::   it_offset    ! local time offset variable
      LOGICAL  ::   llnxtyr      ! open next year  file?
      LOGICAL  ::   llnxtmth     ! open next month file?
      LOGICAL  ::   llstop       ! stop is the file does not exist
      LOGICAL  ::   ll_firstcall ! true if this is the first call to fld_read for this set of fields
      REAL(wp) ::   ztinta       ! ratio applied to after  records when doing time interpolation
      REAL(wp) ::   ztintb       ! ratio applied to before records when doing time interpolation
      CHARACTER(LEN=1000) ::   clfmt  ! write format
      !!---------------------------------------------------------------------
      ll_firstcall = kt == nit000
      IF( PRESENT(kit) )   ll_firstcall = ll_firstcall .and. kit == 1

      IF ( nn_components == jp_iam_sas ) THEN   ;   it_offset = nn_fsbc
      ELSE                                      ;   it_offset = 0
      ENDIF
      IF( PRESENT(kt_offset) )   it_offset = kt_offset

      ! Note that shifting time to be centrered in the middle of sbc time step impacts only nsec_* variables of the calendar 
      IF( present(kit) ) THEN   ! ignore kn_fsbc in this case
         isecsbc = nsec_year + nsec1jan000 + (kit+it_offset)*NINT( rdt/REAL(nn_baro,wp) )
      ELSE                      ! middle of sbc time step
         isecsbc = nsec_year + nsec1jan000 + NINT(0.5 * REAL(kn_fsbc - 1,wp) * rdt) + it_offset * NINT(rdt)
      ENDIF
      imf = SIZE( sd )
      !
      IF( ll_firstcall ) THEN                      ! initialization
         DO jf = 1, imf 
            IF( TRIM(sd(jf)%clrootname) == 'NOT USED' )   CYCLE
            CALL fld_init( kn_fsbc, sd(jf) )       ! read each before field (put them in after as they will be swapped)
         END DO
         IF( lwp ) CALL wgt_print()                ! control print
      ENDIF
      !                                            ! ====================================== !
      IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN         ! update field at each kn_fsbc time-step !
         !                                         ! ====================================== !
         !
         DO jf = 1, imf                            ! ---   loop over field   --- !

            IF( TRIM(sd(jf)%clrootname) == 'NOT USED' )   CYCLE
                      
            IF( isecsbc > sd(jf)%nrec_a(2) .OR. ll_firstcall ) THEN    ! read/update the after data?

               sd(jf)%nrec_b(:) = sd(jf)%nrec_a(:)                                  ! swap before record informations
               sd(jf)%rotn(1) = sd(jf)%rotn(2)                                      ! swap before rotate informations
               IF( sd(jf)%ln_tint )   sd(jf)%fdta(:,:,:,1) = sd(jf)%fdta(:,:,:,2)   ! swap before record field

               CALL fld_rec( kn_fsbc, sd(jf), kt_offset = it_offset, kit = kit )    ! update after record informations

               ! if kn_fsbc*rdt is larger than freqh (which is kind of odd),
               ! it is possible that the before value is no more the good one... we have to re-read it
               ! if before is not the last record of the file currently opened and after is the first record to be read
               ! in a new file which means after = 1 (the file to be opened corresponds to the current time)
               ! or after = nreclast + 1 (the file to be opened corresponds to a future time step)
               IF( .NOT. ll_firstcall .AND. sd(jf)%ln_tint .AND. sd(jf)%nrec_b(1) /= sd(jf)%nreclast &
                  &                   .AND. MOD( sd(jf)%nrec_a(1), sd(jf)%nreclast ) == 1 ) THEN
                  itmp = sd(jf)%nrec_a(1)                       ! temporary storage
                  sd(jf)%nrec_a(1) = sd(jf)%nreclast            ! read the last record of the file currently opened
                  CALL fld_get( sd(jf) )                        ! read after data
                  sd(jf)%fdta(:,:,:,1) = sd(jf)%fdta(:,:,:,2)   ! re-swap before record field
                  sd(jf)%nrec_b(1) = sd(jf)%nrec_a(1)           ! update before record informations
                  sd(jf)%nrec_b(2) = sd(jf)%nrec_a(2) - NINT( sd(jf)%freqh * 3600. )  ! assume freq to be in hours in this case
                  sd(jf)%rotn(1)   = sd(jf)%rotn(2)             ! update before rotate informations
                  sd(jf)%nrec_a(1) = itmp                       ! move back to after record 
               ENDIF

               CALL fld_clopn( sd(jf) )   ! Do we need to open a new year/month/week/day file?
               
               IF( sd(jf)%ln_tint ) THEN
                  
                  ! if kn_fsbc*rdt is larger than freqh (which is kind of odd),
                  ! it is possible that the before value is no more the good one... we have to re-read it
                  ! if before record is not just just before the after record...
                  IF( .NOT. ll_firstcall .AND. MOD( sd(jf)%nrec_a(1), sd(jf)%nreclast ) /= 1 &
                     &                   .AND. sd(jf)%nrec_b(1) /= sd(jf)%nrec_a(1) - 1 ) THEN   
                     sd(jf)%nrec_a(1) = sd(jf)%nrec_a(1) - 1       ! move back to before record
                     CALL fld_get( sd(jf) )                        ! read after data
                     sd(jf)%fdta(:,:,:,1) = sd(jf)%fdta(:,:,:,2)   ! re-swap before record field
                     sd(jf)%nrec_b(1) = sd(jf)%nrec_a(1)           ! update before record informations
                     sd(jf)%nrec_b(2) = sd(jf)%nrec_a(2) - NINT( sd(jf)%freqh * 3600. )  ! assume freq to be in hours in this case
                     sd(jf)%rotn(1)   = sd(jf)%rotn(2)             ! update before rotate informations
                     sd(jf)%nrec_a(1) = sd(jf)%nrec_a(1) + 1       ! move back to after record
                  ENDIF
               ENDIF ! temporal interpolation?

               ! do we have to change the year/month/week/day of the forcing field?? 
               ! if we do time interpolation we will need to open next year/month/week/day file before the end of the current
               ! one. If so, we are still before the end of the year/month/week/day when calling fld_rec so sd(jf)%nrec_a(1)
               ! will be larger than the record number that should be read for current year/month/week/day
               ! do we need next file data?
               ! This applies to both cases with or without time interpolation
               IF( sd(jf)%nrec_a(1) > sd(jf)%nreclast ) THEN
                  
                  sd(jf)%nrec_a(1) = sd(jf)%nrec_a(1) - sd(jf)%nreclast   ! 
                  
                  IF( .NOT. ( sd(jf)%ln_clim .AND. sd(jf)%cltype == 'yearly' ) ) THEN   ! close/open the current/new file
                     
                     llnxtmth = sd(jf)%cltype == 'monthly' .OR. nday == nmonth_len(nmonth)      ! open next month file?
                     llnxtyr  = sd(jf)%cltype == 'yearly'  .OR. (nmonth == 12 .AND. llnxtmth)   ! open next year  file?

                     ! if the run finishes at the end of the current year/month/week/day, we will allow next
                     ! year/month/week/day file to be not present. If the run continue further than the current
                     ! year/month/week/day, next year/month/week/day file must exist
                     isecend = nsec_year + nsec1jan000 + (nitend - kt) * NINT(rdt)   ! second at the end of the run
                     llstop = isecend > sd(jf)%nrec_a(2)                             ! read more than 1 record of next year
                     ! we suppose that the date of next file is next day (should be ok even for weekly files...)
                     CALL fld_clopn( sd(jf), nyear  + COUNT((/llnxtyr /))                                           ,         &
                        &                    nmonth + COUNT((/llnxtmth/)) - 12                 * COUNT((/llnxtyr /)),         &
                        &                    nday   + 1                   - nmonth_len(nmonth) * COUNT((/llnxtmth/)), llstop )

                     IF( sd(jf)%num <= 0 .AND. .NOT. llstop ) THEN    ! next year file does not exist
                        CALL ctl_warn('next year/month/week/day file: '//TRIM(sd(jf)%clname)//     &
                           &     ' not present -> back to current year/month/day')
                        CALL fld_clopn( sd(jf) )               ! back to the current year/month/day
                        sd(jf)%nrec_a(1) = sd(jf)%nreclast     ! force to read the last record in the current year file
                     ENDIF
                     
                  ENDIF
               ENDIF   ! open need next file?
                  
               ! read after data
               CALL fld_get( sd(jf) )
               
            ENDIF   ! read new data?
         END DO                                    ! --- end loop over field --- !

         CALL fld_rot( kt, sd )                    ! rotate vector before/now/after fields if needed

         DO jf = 1, imf                            ! ---   loop over field   --- !
            !
            IF( TRIM(sd(jf)%clrootname) == 'NOT USED' )   CYCLE
            !
            IF( sd(jf)%ln_tint ) THEN              ! temporal interpolation
               IF(lwp .AND. kt - nit000 <= 100 ) THEN 
                  clfmt = "('   fld_read: var ', a, ' kt = ', i8, ' (', f9.4,' days), Y/M/D = ', i4.4,'/', i2.2,'/', i2.2," //   &
                     &    "', records b/a: ', i6.4, '/', i6.4, ' (days ', f9.4,'/', f9.4, ')')"
                  WRITE(numout, clfmt)  TRIM( sd(jf)%clvar ), kt, REAL(isecsbc,wp)/rday, nyear, nmonth, nday,   &            
                     & sd(jf)%nrec_b(1), sd(jf)%nrec_a(1), REAL(sd(jf)%nrec_b(2),wp)/rday, REAL(sd(jf)%nrec_a(2),wp)/rday
                  WRITE(numout, *) '      it_offset is : ',it_offset
               ENDIF
               ! temporal interpolation weights
               ztinta =  REAL( isecsbc - sd(jf)%nrec_b(2), wp ) / REAL( sd(jf)%nrec_a(2) - sd(jf)%nrec_b(2), wp )
               ztintb =  1. - ztinta
               sd(jf)%fnow(:,:,:) = ztintb * sd(jf)%fdta(:,:,:,1) + ztinta * sd(jf)%fdta(:,:,:,2)
            ELSE   ! nothing to do...
               IF(lwp .AND. kt - nit000 <= 100 ) THEN
                  clfmt = "('   fld_read: var ', a, ' kt = ', i8,' (', f9.4,' days), Y/M/D = ', i4.4,'/', i2.2,'/', i2.2," //   &
                     &    "', record: ', i6.4, ' (days ', f9.4, ' <-> ', f9.4, ')')"
                  WRITE(numout, clfmt) TRIM(sd(jf)%clvar), kt, REAL(isecsbc,wp)/rday, nyear, nmonth, nday,    &
                     &                 sd(jf)%nrec_a(1), REAL(sd(jf)%nrec_b(2),wp)/rday, REAL(sd(jf)%nrec_a(2),wp)/rday
               ENDIF
            ENDIF
            !
            IF( kt == nitend - kn_fsbc + 1 )   CALL iom_close( sd(jf)%num )   ! Close the input files

         END DO                                    ! --- end loop over field --- !
         !
      ENDIF
      !
   END SUBROUTINE fld_read


   SUBROUTINE fld_init( kn_fsbc, sdjf )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_init  ***
      !!
      !! ** Purpose :  - first call to fld_rec to define before values
      !!               - if time interpolation, read before data 
      !!----------------------------------------------------------------------
      INTEGER  , INTENT(in   ) ::   kn_fsbc      ! sbc computation period (in time step) 
      TYPE(FLD), INTENT(inout) ::   sdjf         ! input field related variables
      !!
      LOGICAL :: llprevyr              ! are we reading previous year  file?
      LOGICAL :: llprevmth             ! are we reading previous month file?
      LOGICAL :: llprevweek            ! are we reading previous week  file?
      LOGICAL :: llprevday             ! are we reading previous day   file?
      LOGICAL :: llprev                ! llprevyr .OR. llprevmth .OR. llprevweek .OR. llprevday
      INTEGER :: idvar                 ! variable id 
      INTEGER :: inrec                 ! number of record existing for this variable
      INTEGER :: iyear, imonth, iday   ! first day of the current file in yyyy mm dd
      INTEGER :: isec_week             ! number of seconds since start of the weekly file
      CHARACTER(LEN=1000) ::   clfmt   ! write format
      !!---------------------------------------------------------------------
      !
      llprevyr   = .FALSE.
      llprevmth  = .FALSE.
      llprevweek = .FALSE.
      llprevday  = .FALSE.
      isec_week  = 0
      !
      ! define record informations
      CALL fld_rec( kn_fsbc, sdjf, ldbefore = .TRUE. )  ! return before values in sdjf%nrec_a (as we will swap it later)
      !
      ! Note that shifting time to be centrered in the middle of sbc time step impacts only nsec_* variables of the calendar 
      !
      IF( sdjf%ln_tint ) THEN ! we need to read the previous record and we will put it in the current record structure
         !
         IF( sdjf%nrec_a(1) == 0  ) THEN   ! we redefine record sdjf%nrec_a(1) with the last record of previous year file
            IF    ( NINT(sdjf%freqh) == -12 ) THEN   ! yearly mean
               IF( sdjf%cltype == 'yearly' ) THEN             ! yearly file
                  sdjf%nrec_a(1) = 1                                                       ! force to read the unique record
                  llprevyr  = .NOT. sdjf%ln_clim                                           ! use previous year  file?
               ELSE
                  CALL ctl_stop( "fld_init: yearly mean file must be in a yearly type of file: "//TRIM(sdjf%clrootname) )
               ENDIF
            ELSEIF( NINT(sdjf%freqh) ==  -1 ) THEN   ! monthly mean
               IF( sdjf%cltype == 'monthly' ) THEN            ! monthly file
                  sdjf%nrec_a(1) = 1                                                       ! force to read the unique record
                  llprevmth = .TRUE.                                                       ! use previous month file?
                  llprevyr  = llprevmth .AND. nmonth == 1                                  ! use previous year  file?
               ELSE                                           ! yearly file
                  sdjf%nrec_a(1) = 12                                                      ! force to read december mean
                  llprevyr = .NOT. sdjf%ln_clim                                            ! use previous year  file?
               ENDIF
            ELSE                                     ! higher frequency mean (in hours) 
               IF    ( sdjf%cltype      == 'monthly' ) THEN   ! monthly file
                  sdjf%nrec_a(1) = NINT( 24. * REAL(nmonth_len(nmonth-1),wp) / sdjf%freqh )! last record of previous month
                  llprevmth = .TRUE.                                                       ! use previous month file?
                  llprevyr  = llprevmth .AND. nmonth == 1                                  ! use previous year  file?
               ELSEIF( sdjf%cltype(1:4) == 'week'    ) THEN   ! weekly file
                  llprevweek = .TRUE.                                                      ! use previous week  file?
                  sdjf%nrec_a(1) = NINT( 24. * 7. / sdjf%freqh )                           ! last record of previous week
                  isec_week = NINT(rday) * 7                                               ! add a shift toward previous week
               ELSEIF( sdjf%cltype      == 'daily'   ) THEN   ! daily file
                  sdjf%nrec_a(1) = NINT( 24. / sdjf%freqh )                                ! last record of previous day
                  llprevday = .TRUE.                                                       ! use previous day   file?
                  llprevmth = llprevday .AND. nday   == 1                                  ! use previous month file?
                  llprevyr  = llprevmth .AND. nmonth == 1                                  ! use previous year  file?
               ELSE                                           ! yearly file
                  sdjf%nrec_a(1) = NINT( 24. * REAL(nyear_len(0),wp) / sdjf%freqh )        ! last record of previous year 
                  llprevyr = .NOT. sdjf%ln_clim                                            ! use previous year  file?
               ENDIF
            ENDIF
         ENDIF
         !
         IF ( sdjf%cltype(1:4) == 'week' ) THEN
            isec_week = isec_week + ksec_week( sdjf%cltype(6:8) )   ! second since the beginning of the week
            llprevmth = isec_week > nsec_month                      ! longer time since the beginning of the week than the month
            llprevyr  = llprevmth .AND. nmonth == 1
         ENDIF
         llprev = llprevyr .OR. llprevmth .OR. llprevweek .OR. llprevday
         !
         iyear  = nyear  - COUNT((/llprevyr /))
         imonth = nmonth - COUNT((/llprevmth/)) + 12 * COUNT((/llprevyr /))
         iday   = nday   - COUNT((/llprevday/)) + nmonth_len(nmonth-1) * COUNT((/llprevmth/)) - isec_week / NINT(rday)
         !
         CALL fld_clopn( sdjf, iyear, imonth, iday, .NOT. llprev )
         !
         ! if previous year/month/day file does not exist, we switch to the current year/month/day
         IF( llprev .AND. sdjf%num <= 0 ) THEN
            CALL ctl_warn( 'previous year/month/week/day file: '//TRIM(sdjf%clrootname)//   &
               &           ' not present -> back to current year/month/week/day' )
            ! we force to read the first record of the current year/month/day instead of last record of previous year/month/day
            llprev = .FALSE.
            sdjf%nrec_a(1) = 1
            CALL fld_clopn( sdjf )
         ENDIF
         !
         IF( llprev ) THEN   ! check if the record sdjf%nrec_a(1) exists in the file
            idvar = iom_varid( sdjf%num, sdjf%clvar )                                        ! id of the variable sdjf%clvar
            IF( idvar <= 0 )   RETURN
            inrec = iom_file( sdjf%num )%dimsz( iom_file( sdjf%num )%ndims(idvar), idvar )   ! size of the last dim of idvar
            sdjf%nrec_a(1) = MIN( sdjf%nrec_a(1), inrec )   ! make sure we select an existing record
         ENDIF
         !
         ! read before data in after arrays(as we will swap it later)
         CALL fld_get( sdjf )
         !
         clfmt = "('   fld_init : time-interpolation for ', a, ' read previous record = ', i6, ' at time = ', f7.2, ' days')"
         IF(lwp) WRITE(numout, clfmt) TRIM(sdjf%clvar), sdjf%nrec_a(1), REAL(sdjf%nrec_a(2),wp)/rday
         !
      ENDIF
      !
   END SUBROUTINE fld_init


   SUBROUTINE fld_rec( kn_fsbc, sdjf, ldbefore, kit, kt_offset )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_rec  ***
      !!
      !! ** Purpose : Compute
      !!              if sdjf%ln_tint = .TRUE.
      !!                  nrec_a: record number and its time (nrec_b is obtained from nrec_a when swapping)
      !!              if sdjf%ln_tint = .FALSE.
      !!                  nrec_a(1): record number
      !!                  nrec_b(2) and nrec_a(2): time of the beginning and end of the record
      !!----------------------------------------------------------------------
      INTEGER  , INTENT(in   )           ::   kn_fsbc   ! sbc computation period (in time step) 
      TYPE(FLD), INTENT(inout)           ::   sdjf      ! input field related variables
      LOGICAL  , INTENT(in   ), OPTIONAL ::   ldbefore  ! sent back before record values (default = .FALSE.)
      INTEGER  , INTENT(in   ), OPTIONAL ::   kit       ! index of barotropic subcycle
      !                                                 ! used only if sdjf%ln_tint = .TRUE.
      INTEGER  , INTENT(in   ), OPTIONAL ::   kt_offset ! Offset of required time level compared to "now"
      !                                                 !   time level in units of time steps.
      !
      LOGICAL  ::   llbefore    ! local definition of ldbefore
      INTEGER  ::   iendrec     ! end of this record (in seconds)
      INTEGER  ::   imth        ! month number
      INTEGER  ::   ifreq_sec   ! frequency mean (in seconds)
      INTEGER  ::   isec_week   ! number of seconds since the start of the weekly file
      INTEGER  ::   it_offset   ! local time offset variable
      REAL(wp) ::   ztmp        ! temporary variable
      !!----------------------------------------------------------------------
      !
      ! Note that shifting time to be centrered in the middle of sbc time step impacts only nsec_* variables of the calendar 
      !
      IF( PRESENT(ldbefore) ) THEN   ;   llbefore = ldbefore .AND. sdjf%ln_tint   ! needed only if sdjf%ln_tint = .TRUE.
      ELSE                           ;   llbefore = .FALSE.
      ENDIF
      !
      IF ( nn_components == jp_iam_sas ) THEN   ;   it_offset = nn_fsbc
      ELSE                                      ;   it_offset = 0
      ENDIF
      IF( PRESENT(kt_offset) )      it_offset = kt_offset
      IF( PRESENT(kit) ) THEN   ;   it_offset = ( kit + it_offset ) * NINT( rdt/REAL(nn_baro,wp) )
      ELSE                      ;   it_offset =         it_offset   * NINT(       rdt            )
      ENDIF
      !
      !                                           ! =========== !
      IF    ( NINT(sdjf%freqh) == -12 ) THEN      ! yearly mean
         !                                        ! =========== !
         !
         IF( sdjf%ln_tint ) THEN                  ! time interpolation, shift by 1/2 record
            !
            !                  INT( ztmp )
            !                     /|\
            !                    1 |    *----
            !                    0 |----(              
            !                      |----+----|--> time
            !                      0   /|\   1   (nday/nyear_len(1))
            !                           |   
            !                           |   
            !       forcing record :    1 
            !                            
            ztmp =  REAL( nsec_year, wp ) / ( REAL( nyear_len(1), wp ) * rday ) + 0.5 &
               &  + REAL( it_offset, wp ) / ( REAL( nyear_len(1), wp ) * rday )
            sdjf%nrec_a(1) = 1 + INT( ztmp ) - COUNT((/llbefore/))
            ! swap at the middle of the year
            IF( llbefore ) THEN   ;   sdjf%nrec_a(2) = nsec1jan000 - (1 - INT(ztmp)) * NINT(0.5 * rday) * nyear_len(0) + &
                                    & INT(ztmp) * NINT( 0.5 * rday) * nyear_len(1) 
            ELSE                  ;   sdjf%nrec_a(2) = nsec1jan000 + (1 - INT(ztmp)) * NINT(0.5 * rday) * nyear_len(1) + &
                                    & INT(ztmp) * INT(rday) * nyear_len(1) + INT(ztmp) * NINT( 0.5 * rday) * nyear_len(2) 
            ENDIF
         ELSE                                     ! no time interpolation
            sdjf%nrec_a(1) = 1
            sdjf%nrec_a(2) = NINT(rday) * nyear_len(1) + nsec1jan000   ! swap at the end    of the year
            sdjf%nrec_b(2) = nsec1jan000                               ! beginning of the year (only for print)
         ENDIF
         !
         !                                        ! ============ !
      ELSEIF( NINT(sdjf%freqh) ==  -1 ) THEN      ! monthly mean !
         !                                        ! ============ !
         !
         IF( sdjf%ln_tint ) THEN                  ! time interpolation, shift by 1/2 record
            !
            !                  INT( ztmp )
            !                     /|\
            !                    1 |    *----
            !                    0 |----(              
            !                      |----+----|--> time
            !                      0   /|\   1   (nday/nmonth_len(nmonth))
            !                           |   
            !                           |   
            !       forcing record :  nmonth 
            !                            
            ztmp =  REAL( nsec_month, wp ) / ( REAL( nmonth_len(nmonth), wp ) * rday ) + 0.5 &
           &      + REAL(  it_offset, wp ) / ( REAL( nmonth_len(nmonth), wp ) * rday )
            imth = nmonth + INT( ztmp ) - COUNT((/llbefore/))
            IF( sdjf%cltype == 'monthly' ) THEN   ;   sdjf%nrec_a(1) = 1 + INT( ztmp ) - COUNT((/llbefore/))
            ELSE                                  ;   sdjf%nrec_a(1) = imth
            ENDIF
            sdjf%nrec_a(2) = nmonth_half(   imth ) + nsec1jan000   ! swap at the middle of the month
         ELSE                                    ! no time interpolation
            IF( sdjf%cltype == 'monthly' ) THEN   ;   sdjf%nrec_a(1) = 1
            ELSE                                  ;   sdjf%nrec_a(1) = nmonth
            ENDIF
            sdjf%nrec_a(2) =  nmonth_end(nmonth  ) + nsec1jan000   ! swap at the end    of the month
            sdjf%nrec_b(2) =  nmonth_end(nmonth-1) + nsec1jan000   ! beginning of the month (only for print)
         ENDIF
         !
         !                                        ! ================================ !
      ELSE                                        ! higher frequency mean (in hours)
         !                                        ! ================================ !
         !
         ifreq_sec = NINT( sdjf%freqh * 3600. )                                         ! frequency mean (in seconds)
         IF( sdjf%cltype(1:4) == 'week' )   isec_week = ksec_week( sdjf%cltype(6:8) )   ! since the first day of the current week
         ! number of second since the beginning of the file
         IF(     sdjf%cltype      == 'monthly' ) THEN   ;   ztmp = REAL(nsec_month,wp)  ! since the first day of the current month
         ELSEIF( sdjf%cltype(1:4) == 'week'    ) THEN   ;   ztmp = REAL(isec_week ,wp)  ! since the first day of the current week
         ELSEIF( sdjf%cltype      == 'daily'   ) THEN   ;   ztmp = REAL(nsec_day  ,wp)  ! since 00h of the current day
         ELSE                                           ;   ztmp = REAL(nsec_year ,wp)  ! since 00h on Jan 1 of the current year
         ENDIF
         ztmp = ztmp + 0.5 * REAL(kn_fsbc - 1, wp) * rdt + REAL( it_offset, wp )        ! centrered in the middle of sbc time step
         ztmp = ztmp + 0.01 * rdt                                                       ! avoid truncation error 
         IF( sdjf%ln_tint ) THEN                 ! time interpolation, shift by 1/2 record
            !
            !          INT( ztmp/ifreq_sec + 0.5 )
            !                     /|\
            !                    2 |        *-----(
            !                    1 |  *-----(
            !                    0 |--(              
            !                      |--+--|--+--|--+--|--> time
            !                      0 /|\ 1 /|\ 2 /|\ 3    (ztmp/ifreq_sec)
            !                         |     |     |
            !                         |     |     |
            !       forcing record :  1     2     3
            !                   
            ztmp= ztmp / REAL(ifreq_sec, wp) + 0.5
         ELSE                                    ! no time interpolation
            !
            !           INT( ztmp/ifreq_sec )
            !                     /|\
            !                    2 |           *-----(
            !                    1 |     *-----(
            !                    0 |-----(              
            !                      |--+--|--+--|--+--|--> time
            !                      0 /|\ 1 /|\ 2 /|\ 3    (ztmp/ifreq_sec)
            !                         |     |     |
            !                         |     |     |
            !       forcing record :  1     2     3
            !                            
            ztmp= ztmp / REAL(ifreq_sec, wp)
         ENDIF
         sdjf%nrec_a(1) = 1 + INT( ztmp ) - COUNT((/llbefore/))   ! record number to be read

         iendrec = ifreq_sec * sdjf%nrec_a(1) + nsec1jan000       ! end of this record (in second)
         ! add the number of seconds between 00h Jan 1 and the end of previous month/week/day (ok if nmonth=1)
         IF( sdjf%cltype      == 'monthly' )   iendrec = iendrec + NINT(rday) * SUM(nmonth_len(1:nmonth -1))
         IF( sdjf%cltype(1:4) == 'week'    )   iendrec = iendrec + ( nsec_year - isec_week )
         IF( sdjf%cltype      == 'daily'   )   iendrec = iendrec + NINT(rday) * ( nday_year - 1 )
         IF( sdjf%ln_tint ) THEN
             sdjf%nrec_a(2) = iendrec - ifreq_sec / 2        ! swap at the middle of the record
         ELSE
             sdjf%nrec_a(2) = iendrec                        ! swap at the end    of the record
             sdjf%nrec_b(2) = iendrec - ifreq_sec            ! beginning of the record (only for print)
         ENDIF
         !
      ENDIF
      !
      IF( .NOT. sdjf%ln_tint ) sdjf%nrec_a(2) = sdjf%nrec_a(2) - 1   ! last second belongs to bext record : *----(
      !
   END SUBROUTINE fld_rec


   SUBROUTINE fld_get( sdjf )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_get  ***
      !!
      !! ** Purpose :   read the data
      !!----------------------------------------------------------------------
      TYPE(FLD)        , INTENT(inout) ::   sdjf   ! input field related variables
      !
      INTEGER ::   ipk      ! number of vertical levels of sdjf%fdta ( 2D: ipk=1 ; 3D: ipk=jpk )
      INTEGER ::   iw       ! index into wgts array
      INTEGER ::   ipdom    ! index of the domain
      INTEGER ::   idvar    ! variable ID
      INTEGER ::   idmspc   ! number of spatial dimensions
      LOGICAL ::   lmoor    ! C1D case: point data
      !!---------------------------------------------------------------------
      !
      ipk = SIZE( sdjf%fnow, 3 )
      !
      IF( ASSOCIATED(sdjf%imap) ) THEN
         IF( sdjf%ln_tint ) THEN   ;   CALL fld_map( sdjf%num, sdjf%clvar, sdjf%fdta(:,:,:,2), sdjf%nrec_a(1),   &
            &                                        sdjf%imap, sdjf%igrd, sdjf%ibdy, sdjf%ltotvel, sdjf%lzint )
         ELSE                      ;   CALL fld_map( sdjf%num, sdjf%clvar, sdjf%fnow(:,:,:  ), sdjf%nrec_a(1),   &
            &                                        sdjf%imap, sdjf%igrd, sdjf%ibdy, sdjf%ltotvel, sdjf%lzint )
         ENDIF
      ELSE IF( LEN(TRIM(sdjf%wgtname)) > 0 ) THEN
         CALL wgt_list( sdjf, iw )
         IF( sdjf%ln_tint ) THEN   ;   CALL fld_interp( sdjf%num, sdjf%clvar, iw, ipk, sdjf%fdta(:,:,:,2),          & 
            &                                                                          sdjf%nrec_a(1), sdjf%lsmname )
         ELSE                      ;   CALL fld_interp( sdjf%num, sdjf%clvar, iw, ipk, sdjf%fnow(:,:,:  ),          &
            &                                                                          sdjf%nrec_a(1), sdjf%lsmname )
         ENDIF
      ELSE
         IF( SIZE(sdjf%fnow, 1) == jpi ) THEN   ;   ipdom = jpdom_data
         ELSE                                   ;   ipdom = jpdom_unknown
         ENDIF
         ! C1D case: If product of spatial dimensions == ipk, then x,y are of
         ! size 1 (point/mooring data): this must be read onto the central grid point
         idvar  = iom_varid( sdjf%num, sdjf%clvar )
         idmspc = iom_file ( sdjf%num )%ndims( idvar )
         IF( iom_file( sdjf%num )%luld( idvar ) )   idmspc = idmspc - 1
         lmoor  = (  idmspc == 0 .OR. PRODUCT( iom_file( sdjf%num )%dimsz( 1:MAX(idmspc,1) ,idvar ) ) == ipk  )
         !
         SELECT CASE( ipk )
         CASE(1)
            IF( lk_c1d .AND. lmoor ) THEN
               IF( sdjf%ln_tint ) THEN
                  CALL iom_get( sdjf%num, sdjf%clvar, sdjf%fdta(2,2,1,2), sdjf%nrec_a(1) )
                  CALL lbc_lnk( 'fldread', sdjf%fdta(:,:,1,2),'Z',1. )
               ELSE
                  CALL iom_get( sdjf%num, sdjf%clvar, sdjf%fnow(2,2,1  ), sdjf%nrec_a(1) )
                  CALL lbc_lnk( 'fldread', sdjf%fnow(:,:,1  ),'Z',1. )
               ENDIF
            ELSE
               IF( sdjf%ln_tint ) THEN   ;   CALL iom_get( sdjf%num, ipdom, sdjf%clvar, sdjf%fdta(:,:,1,2), sdjf%nrec_a(1) )
               ELSE                      ;   CALL iom_get( sdjf%num, ipdom, sdjf%clvar, sdjf%fnow(:,:,1  ), sdjf%nrec_a(1) )
               ENDIF
            ENDIF
         CASE DEFAULT
            IF (lk_c1d .AND. lmoor ) THEN
               IF( sdjf%ln_tint ) THEN
                  CALL iom_get( sdjf%num, jpdom_unknown, sdjf%clvar, sdjf%fdta(2,2,:,2), sdjf%nrec_a(1) )
                  CALL lbc_lnk( 'fldread', sdjf%fdta(:,:,:,2),'Z',1. )
               ELSE
                  CALL iom_get( sdjf%num, jpdom_unknown, sdjf%clvar, sdjf%fnow(2,2,:  ), sdjf%nrec_a(1) )
                  CALL lbc_lnk( 'fldread', sdjf%fnow(:,:,:  ),'Z',1. )
               ENDIF
            ELSE
               IF( sdjf%ln_tint ) THEN   ;   CALL iom_get( sdjf%num, ipdom, sdjf%clvar, sdjf%fdta(:,:,:,2), sdjf%nrec_a(1) )
               ELSE                      ;   CALL iom_get( sdjf%num, ipdom, sdjf%clvar, sdjf%fnow(:,:,:  ), sdjf%nrec_a(1) )
               ENDIF
            ENDIF
         END SELECT
      ENDIF
      !
      sdjf%rotn(2) = .false.   ! vector not yet rotated
      !
   END SUBROUTINE fld_get

   
   SUBROUTINE fld_map( knum, cdvar, pdta, krec, kmap, kgrd, kbdy, ldtotvel, ldzint )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_map  ***
      !!
      !! ** Purpose :   read global data from file and map onto local data
      !!                using a general mapping (for open boundaries)
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   knum         ! stream number
      CHARACTER(LEN=*)          , INTENT(in   ) ::   cdvar        ! variable name
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdta         ! bdy output field on model grid
      INTEGER                   , INTENT(in   ) ::   krec         ! record number to read (ie time slice)
      INTEGER , DIMENSION(:)    , INTENT(in   ) ::   kmap         ! global-to-local bdy mapping indices
      ! optional variables used for vertical interpolation:
      INTEGER, OPTIONAL         , INTENT(in   ) ::   kgrd         ! grid type (t, u, v)
      INTEGER, OPTIONAL         , INTENT(in   ) ::   kbdy         ! bdy number
      LOGICAL, OPTIONAL         , INTENT(in   ) ::   ldtotvel     ! true if total ( = barotrop + barocline) velocity
      LOGICAL, OPTIONAL         , INTENT(in   ) ::   ldzint       ! true if 3D variable requires a vertical interpolation
      !!
      INTEGER                                   ::   ipi          ! length of boundary data on local process
      INTEGER                                   ::   ipj          ! length of dummy dimension ( = 1 )
      INTEGER                                   ::   ipk          ! number of vertical levels of pdta ( 2D: ipk=1 ; 3D: ipk=jpk )
      INTEGER                                   ::   ipkb         ! number of vertical levels in boundary data file
      INTEGER                                   ::   idvar        ! variable ID
      INTEGER                                   ::   indims       ! number of dimensions of the variable
      INTEGER, DIMENSION(4)                     ::   idimsz       ! size of variable dimensions 
      REAL(wp)                                  ::   zfv          ! fillvalue 
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::   zz_read      ! work space for global boundary data
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::   zdta_read    ! work space local data requiring vertical interpolation
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::   zdta_read_z  ! work space local data requiring vertical interpolation
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::   zdta_read_dz ! work space local data requiring vertical interpolation
      CHARACTER(LEN=1),DIMENSION(3)             ::   clgrid
      LOGICAL                                   ::   lluld        ! is the variable using the unlimited dimension
      LOGICAL                                   ::   llzint       ! local value of ldzint
      !!---------------------------------------------------------------------
      !
      clgrid = (/'t','u','v'/)
      !
      ipi = SIZE( pdta, 1 )
      ipj = SIZE( pdta, 2 )   ! must be equal to 1
      ipk = SIZE( pdta, 3 )
      !
      llzint = .FALSE.
      IF( PRESENT(ldzint) )   llzint = ldzint
      !
      idvar = iom_varid( knum, cdvar, kndims = indims, kdimsz = idimsz, lduld = lluld  )
      IF( indims == 4 .OR. ( indims == 3 .AND. .NOT. lluld ) ) THEN   ;   ipkb = idimsz(3)   ! xy(zl)t or xy(zl)
      ELSE                                                            ;   ipkb = 1           ! xy or xyt
      ENDIF
      !
      ALLOCATE( zz_read( idimsz(1), idimsz(2), ipkb ) )  ! ++++++++ !!! this can be very big...         
      !
      IF( ipk == 1 ) THEN

         IF( ipkb /= 1 ) CALL ctl_stop( 'fld_map : we must have ipkb = 1 to read surface data' )
         CALL iom_get ( knum, jpdom_unknown, cdvar, zz_read(:,:,1), krec )   ! call iom_get with a 2D file
         CALL fld_map_core( zz_read, kmap, pdta )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Do we include something here to adjust barotropic velocities !
      ! in case of a depth difference between bdy files and          !
      ! bathymetry in the case ln_totvel = .false. and ipkb>0?       !
      ! [as the enveloping and parital cells could change H]         !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ELSE
         !
         CALL iom_get ( knum, jpdom_unknown, cdvar, zz_read(:,:,:), krec )   ! call iom_get with a 3D file
         !
         IF( ipkb /= ipk .OR. llzint ) THEN   ! boundary data not on model vertical grid : vertical interpolation
            !
            IF( ipk == jpk .AND. iom_varid(knum,'gdep'//clgrid(kgrd)) /= -1 .AND. iom_varid(knum,'e3'//clgrid(kgrd)) /= -1 ) THEN
               
               ALLOCATE( zdta_read(ipi,ipj,ipkb), zdta_read_z(ipi,ipj,ipkb), zdta_read_dz(ipi,ipj,ipkb) )
                
               CALL fld_map_core( zz_read, kmap, zdta_read )
               CALL iom_get ( knum, jpdom_unknown, 'gdep'//clgrid(kgrd), zz_read )   ! read only once? Potential temporal evolution?
               CALL fld_map_core( zz_read, kmap, zdta_read_z )
               CALL iom_get ( knum, jpdom_unknown,   'e3'//clgrid(kgrd), zz_read )   ! read only once? Potential temporal evolution?
               CALL fld_map_core( zz_read, kmap, zdta_read_dz )
               
               CALL iom_getatt(knum, '_FillValue', zfv, cdvar=cdvar )
               CALL fld_bdy_interp(zdta_read, zdta_read_z, zdta_read_dz, pdta, kgrd, kbdy, zfv, ldtotvel)
               DEALLOCATE( zdta_read, zdta_read_z, zdta_read_dz )
               
            ELSE
               IF( ipk /= jpk ) CALL ctl_stop( 'fld_map : this should be an impossible case...' )
               WRITE(ctmp1,*) 'fld_map : vertical interpolation for bdy variable '//TRIM(cdvar)//' requires ' 
               IF( iom_varid(knum, 'gdep'//clgrid(kgrd)) == -1 ) CALL ctl_stop( ctmp1//'gdep'//clgrid(kgrd)//' variable' )
               IF( iom_varid(knum,   'e3'//clgrid(kgrd)) == -1 ) CALL ctl_stop( ctmp1//  'e3'//clgrid(kgrd)//' variable' )

            ENDIF
            !
         ELSE                            ! bdy data assumed to be the same levels as bdy variables
            !
            CALL fld_map_core( zz_read, kmap, pdta )
            !
         ENDIF   ! ipkb /= ipk
      ENDIF   ! ipk == 1
      
      DEALLOCATE( zz_read )

   END SUBROUTINE fld_map

     
   SUBROUTINE fld_map_core( pdta_read, kmap, pdta_bdy )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_map_core  ***
      !!
      !! ** Purpose :  inner core of fld_map
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdta_read    ! global boundary data
      INTEGER,  DIMENSION(:    ), INTENT(in   ) ::   kmap         ! global-to-local bdy mapping indices
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdta_bdy     ! bdy output field on model grid
      !!
      INTEGER,  DIMENSION(3) ::   idim_read,  idim_bdy            ! arrays dimensions
      INTEGER                ::   ji, jj, jk, jb                  ! loop counters
      INTEGER                ::   im1
      !!---------------------------------------------------------------------
      !
      idim_read = SHAPE( pdta_read )
      idim_bdy  = SHAPE( pdta_bdy  )
      !
      ! in all cases: idim_bdy(2) == 1 .AND. idim_read(1) * idim_read(2) == idim_bdy(1)
      ! structured BDY with rimwidth > 1                     : idim_read(2) == rimwidth /= 1
      ! structured BDY with rimwidth == 1 or unstructured BDY: idim_read(2) == 1
      !
      IF( idim_read(2) > 1 ) THEN    ! structured BDY with rimwidth > 1  
         DO jk = 1, idim_bdy(3)
            DO jb = 1, idim_bdy(1)
               im1 = kmap(jb) - 1
               jj = im1 / idim_read(1) + 1
               ji = MOD( im1, idim_read(1) ) + 1
               pdta_bdy(jb,1,jk) =  pdta_read(ji,jj,jk)
            END DO
         END DO
      ELSE
         DO jk = 1, idim_bdy(3)
            DO jb = 1, idim_bdy(1)   ! horizontal remap of bdy data on the local bdy 
               pdta_bdy(jb,1,jk) = pdta_read(kmap(jb),1,jk)
            END DO
         END DO
      ENDIF
      
   END SUBROUTINE fld_map_core
   
   
   SUBROUTINE fld_bdy_interp(pdta_read, pdta_read_z, pdta_read_dz, pdta, kgrd, kbdy, pfv, ldtotvel)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_bdy_interp  ***
      !!
      !! ** Purpose :   on the fly vertical interpolation to allow the use of
      !!                boundary data from non-native vertical grid
      !!----------------------------------------------------------------------
      USE bdy_oce, ONLY:  idx_bdy         ! indexing for map <-> ij transformation

      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdta_read       ! data read in bdy file
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdta_read_z     ! depth of the data read in bdy file
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdta_read_dz    ! thickness of the levels in bdy file
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdta            ! output field on model grid (2 dimensional)
      REAL(wp)                  , INTENT(in   ) ::   pfv             ! fillvalue of the data read in bdy file
      LOGICAL                   , INTENT(in   ) ::   ldtotvel        ! true if toal ( = barotrop + barocline) velocity
      INTEGER                   , INTENT(in   ) ::   kgrd            ! grid type (t, u, v)
      INTEGER                   , INTENT(in   ) ::   kbdy            ! bdy number
      !!
      INTEGER                  ::   ipi                 ! length of boundary data on local process
      INTEGER                  ::   ipkb                ! number of vertical levels in boundary data file
      INTEGER                  ::   ipkmax              ! number of vertical levels in boundary data file where no mask
      INTEGER                  ::   jb, ji, jj, jk, jkb ! loop counters
      REAL(wp)                 ::   zcoef, zi           ! 
      REAL(wp)                 ::   ztrans, ztrans_new  ! transports
      REAL(wp), DIMENSION(jpk) ::   zdepth, zdhalf      ! level and half-level depth
      !!---------------------------------------------------------------------
     
      ipi  = SIZE( pdta, 1 )
      ipkb = SIZE( pdta_read, 3 )
      
      DO jb = 1, ipi
         ji = idx_bdy(kbdy)%nbi(jb,kgrd)
         jj = idx_bdy(kbdy)%nbj(jb,kgrd)
         !
         ! --- max jk where input data /= FillValue --- !
         ipkmax = 1
         DO jkb = 2, ipkb
            IF( pdta_read(jb,1,jkb) /= pfv )   ipkmax = MAX( ipkmax, jkb )
         END DO
         !
         ! --- calculate depth at t,u,v points --- !
         SELECT CASE( kgrd )                         
         CASE(1)            ! depth of T points:
            zdepth(:) = gdept_n(ji,jj,:)
         CASE(2)            ! depth of U points: we must not use gdept_n as we don't want to do a communication
            !                 --> copy what is done for gdept_n in domvvl...
            zdhalf(1) = 0.0_wp
            zdepth(1) = 0.5_wp * e3uw_n(ji,jj,1)
            DO jk = 2, jpk                               ! vertical sum
               !    zcoef = umask - wumask    ! 0 everywhere tmask = wmask, ie everywhere expect at jk = mikt
               !                              ! 1 everywhere from mbkt to mikt + 1 or 1 (if no isf)
               !                              ! 0.5 where jk = mikt     
               !!gm ???????   BUG ?  gdept_n as well as gde3w_n  does not include the thickness of ISF ??
               zcoef = ( umask(ji,jj,jk) - wumask(ji,jj,jk) )
               zdhalf(jk) = zdhalf(jk-1) + e3u_n(ji,jj,jk-1)
               zdepth(jk) =       zcoef  * ( zdhalf(jk  ) + 0.5 * e3uw_n(ji,jj,jk))  &
                  &         + (1.-zcoef) * ( zdepth(jk-1) +       e3uw_n(ji,jj,jk))
            END DO
         CASE(3)            ! depth of V points: we must not use gdept_n as we don't want to do a communication
            !                 --> copy what is done for gdept_n in domvvl...
            zdhalf(1) = 0.0_wp
            zdepth(1) = 0.5_wp * e3vw_n(ji,jj,1)
            DO jk = 2, jpk                               ! vertical sum
               !    zcoef = vmask - wvmask    ! 0 everywhere tmask = wmask, ie everywhere expect at jk = mikt
               !                              ! 1 everywhere from mbkt to mikt + 1 or 1 (if no isf)
               !                              ! 0.5 where jk = mikt     
               !!gm ???????   BUG ?  gdept_n as well as gde3w_n  does not include the thickness of ISF ??
               zcoef = ( vmask(ji,jj,jk) - wvmask(ji,jj,jk) )
               zdhalf(jk) = zdhalf(jk-1) + e3v_n(ji,jj,jk-1)
               zdepth(jk) =       zcoef  * ( zdhalf(jk  ) + 0.5 * e3vw_n(ji,jj,jk))  &
                  &         + (1.-zcoef) * ( zdepth(jk-1) +       e3vw_n(ji,jj,jk))
            END DO
         END SELECT
         !         
         ! --- interpolate bdy data on the model grid --- !
         DO jk = 1, jpk
            IF(     zdepth(jk) <= pdta_read_z(jb,1,1)      ) THEN   ! above the first level of external data
               pdta(jb,1,jk) = pdta_read(jb,1,1)
            ELSEIF( zdepth(jk) >  pdta_read_z(jb,1,ipkmax) ) THEN   ! below the last level of external data /= FillValue
               pdta(jb,1,jk) = pdta_read(jb,1,ipkmax)
            ELSE                                                    ! inbetween: vertical interpolation between jkb & jkb+1
               DO jkb = 1, ipkmax-1
                  IF( ( ( zdepth(jk) - pdta_read_z(jb,1,jkb) ) * ( zdepth(jk) - pdta_read_z(jb,1,jkb+1) ) ) <= 0._wp ) THEN ! linear interpolation between 2 levels
                     zi = ( zdepth(jk) - pdta_read_z(jb,1,jkb) ) / ( pdta_read_z(jb,1,jkb+1) - pdta_read_z(jb,1,jkb) )
                     pdta(jb,1,jk) = pdta_read(jb,1,jkb) + zi * ( pdta_read(jb,1,jkb+1) - pdta_read(jb,1,jkb) )
                  ENDIF
               END DO
            ENDIF
         END DO
         !
      END DO   ! ipi

      ! --- mask data and adjust transport --- !
      SELECT CASE( kgrd )                         

      CASE(1)                                 ! mask data (probably unecessary)
         DO jb = 1, ipi
            ji = idx_bdy(kbdy)%nbi(jb,kgrd)
            jj = idx_bdy(kbdy)%nbj(jb,kgrd)
            DO jk = 1, jpk                      
               pdta(jb,1,jk) = pdta(jb,1,jk) * tmask(ji,jj,jk)
            END DO
         END DO
         
      CASE(2)                                 ! adjust the U-transport term
         DO jb = 1, ipi
            ji = idx_bdy(kbdy)%nbi(jb,kgrd)
            jj = idx_bdy(kbdy)%nbj(jb,kgrd)
            ztrans = 0._wp
            DO jkb = 1, ipkb                              ! calculate transport on input grid
               IF( pdta_read(jb,1,jkb) /= pfv )   ztrans = ztrans + pdta_read(jb,1,jkb) * pdta_read_dz(jb,1,jkb)
            ENDDO
            ztrans_new = 0._wp
            DO jk = 1, jpk                                ! calculate transport on model grid
               ztrans_new = ztrans_new + pdta(jb,1,jk ) * e3u_n(ji,jj,jk) * umask(ji,jj,jk)
            ENDDO
            DO jk = 1, jpk                                ! make transport correction
               IF(ldtotvel) THEN ! bdy data are total velocity so adjust bt transport term to match input data
                  pdta(jb,1,jk) = ( pdta(jb,1,jk) + ( ztrans - ztrans_new ) * r1_hu_n(ji,jj) ) * umask(ji,jj,jk)
               ELSE              ! we're just dealing with bc velocity so bt transport term should sum to zero
                  pdta(jb,1,jk) = ( pdta(jb,1,jk) + (  0._wp - ztrans_new ) * r1_hu_n(ji,jj) ) * umask(ji,jj,jk)
               ENDIF
            ENDDO
         ENDDO

      CASE(3)                                 ! adjust the V-transport term
         DO jb = 1, ipi
            ji = idx_bdy(kbdy)%nbi(jb,kgrd)
            jj = idx_bdy(kbdy)%nbj(jb,kgrd)
            ztrans = 0._wp
            DO jkb = 1, ipkb                              ! calculate transport on input grid
               IF( pdta_read(jb,1,jkb) /= pfv )   ztrans = ztrans + pdta_read(jb,1,jkb) * pdta_read_dz(jb,1,jkb)
            ENDDO
            ztrans_new = 0._wp
            DO jk = 1, jpk                                ! calculate transport on model grid
               ztrans_new = ztrans_new + pdta(jb,1,jk ) * e3v_n(ji,jj,jk) * vmask(ji,jj,jk)
            ENDDO
            DO jk = 1, jpk                                ! make transport correction
               IF(ldtotvel) THEN ! bdy data are total velocity so adjust bt transport term to match input data
                  pdta(jb,1,jk) = ( pdta(jb,1,jk) + ( ztrans - ztrans_new ) * r1_hv_n(ji,jj) ) * vmask(ji,jj,jk)
               ELSE              ! we're just dealing with bc velocity so bt transport term should sum to zero
                  pdta(jb,1,jk) = ( pdta(jb,1,jk) + (  0._wp - ztrans_new ) * r1_hv_n(ji,jj) ) * vmask(ji,jj,jk)
               ENDIF
            ENDDO
         ENDDO
      END SELECT

   END SUBROUTINE fld_bdy_interp

   
   SUBROUTINE fld_rot( kt, sd )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_rot  ***
      !!
      !! ** Purpose :   Vector fields may need to be rotated onto the local grid direction
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in   ) ::   kt   ! ocean time step
      TYPE(FLD), DIMENSION(:), INTENT(inout) ::   sd   ! input field related variables
      !
      INTEGER ::   ju, jv, jk, jn  ! loop indices
      INTEGER ::   imf             ! size of the structure sd
      INTEGER ::   ill             ! character length
      INTEGER ::   iv              ! indice of V component
      CHARACTER (LEN=100)          ::   clcomp       ! dummy weight name
      REAL(wp), DIMENSION(jpi,jpj) ::   utmp, vtmp   ! temporary arrays for vector rotation
      !!---------------------------------------------------------------------
      !
      !! (sga: following code should be modified so that pairs arent searched for each time
      !
      imf = SIZE( sd )
      DO ju = 1, imf
         IF( TRIM(sd(ju)%clrootname) == 'NOT USED' )   CYCLE
         ill = LEN_TRIM( sd(ju)%vcomp )
         DO jn = 2-COUNT((/sd(ju)%ln_tint/)), 2
            IF( ill > 0 .AND. .NOT. sd(ju)%rotn(jn) ) THEN   ! find vector rotations required             
               IF( sd(ju)%vcomp(1:1) == 'U' ) THEN      ! east-west component has symbolic name starting with 'U'
                  ! look for the north-south component which has same symbolic name but with 'U' replaced with 'V'
                  clcomp = 'V' // sd(ju)%vcomp(2:ill)   ! works even if ill == 1
                  iv = -1
                  DO jv = 1, imf
                     IF( TRIM(sd(jv)%clrootname) == 'NOT USED' )   CYCLE
                     IF( TRIM(sd(jv)%vcomp) == TRIM(clcomp) )   iv = jv
                  END DO
                  IF( iv > 0 ) THEN   ! fields ju and iv are two components which need to be rotated together
                     DO jk = 1, SIZE( sd(ju)%fnow, 3 )
                        IF( sd(ju)%ln_tint )THEN
                           CALL rot_rep( sd(ju)%fdta(:,:,jk,jn), sd(iv)%fdta(:,:,jk,jn), 'T', 'en->i', utmp(:,:) )
                           CALL rot_rep( sd(ju)%fdta(:,:,jk,jn), sd(iv)%fdta(:,:,jk,jn), 'T', 'en->j', vtmp(:,:) )
                           sd(ju)%fdta(:,:,jk,jn) = utmp(:,:)   ;   sd(iv)%fdta(:,:,jk,jn) = vtmp(:,:)
                        ELSE 
                           CALL rot_rep( sd(ju)%fnow(:,:,jk  ), sd(iv)%fnow(:,:,jk  ), 'T', 'en->i', utmp(:,:) )
                           CALL rot_rep( sd(ju)%fnow(:,:,jk  ), sd(iv)%fnow(:,:,jk  ), 'T', 'en->j', vtmp(:,:) )
                           sd(ju)%fnow(:,:,jk   ) = utmp(:,:)   ;   sd(iv)%fnow(:,:,jk   ) = vtmp(:,:)
                        ENDIF
                     END DO
                     sd(ju)%rotn(jn) = .TRUE.               ! vector was rotated 
                     IF( lwp .AND. kt == nit000 )   WRITE(numout,*)   &
                        &   'fld_read: vector pair ('//TRIM(sd(ju)%clvar)//', '//TRIM(sd(iv)%clvar)//') rotated on to model grid'
                  ENDIF
               ENDIF
            ENDIF
         END DO
       END DO
      !
   END SUBROUTINE fld_rot


   SUBROUTINE fld_clopn( sdjf, kyear, kmonth, kday, ldstop )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_clopn  ***
      !!
      !! ** Purpose :   update the file name and close/open the files
      !!----------------------------------------------------------------------
      TYPE(FLD)        , INTENT(inout) ::   sdjf     ! input field related variables
      INTEGER, OPTIONAL, INTENT(in   ) ::   kyear    ! year value
      INTEGER, OPTIONAL, INTENT(in   ) ::   kmonth   ! month value
      INTEGER, OPTIONAL, INTENT(in   ) ::   kday     ! day value
      LOGICAL, OPTIONAL, INTENT(in   ) ::   ldstop   ! stop if open to read a non-existing file (default = .TRUE.)
      !
      LOGICAL  :: llprevyr              ! are we reading previous year  file?
      LOGICAL  :: llprevmth             ! are we reading previous month file?
      INTEGER  :: iyear, imonth, iday   ! first day of the current file in yyyy mm dd
      INTEGER  :: isec_week             ! number of seconds since start of the weekly file
      INTEGER  :: indexyr               ! year undex (O/1/2: previous/current/next)
      REAL(wp) :: zyear_len, zmonth_len ! length (days) of iyear and imonth             ! 
      CHARACTER(len = 256) ::   clname  ! temporary file name
      !!----------------------------------------------------------------------
      IF( PRESENT(kyear) ) THEN                             ! use given values 
         iyear = kyear
         imonth = kmonth
         iday = kday
         IF ( sdjf%cltype(1:4) == 'week' ) THEN             ! find the day of the beginning of the week
            isec_week = ksec_week( sdjf%cltype(6:8) )- (86400 * 8 )  
            llprevmth  = isec_week > nsec_month             ! longer time since beginning of the week than the month
            llprevyr   = llprevmth .AND. nmonth == 1
            iyear  = nyear  - COUNT((/llprevyr /))
            imonth = nmonth - COUNT((/llprevmth/)) + 12 * COUNT((/llprevyr /))
            iday   = nday   + nmonth_len(nmonth-1) * COUNT((/llprevmth/)) - isec_week / NINT(rday)
         ENDIF
      ELSE                                                  ! use current day values
         IF ( sdjf%cltype(1:4) == 'week' ) THEN             ! find the day of the beginning of the week
            isec_week  = ksec_week( sdjf%cltype(6:8) )      ! second since the beginning of the week
            llprevmth  = isec_week > nsec_month             ! longer time since beginning of the week than the month
            llprevyr   = llprevmth .AND. nmonth == 1
         ELSE
            isec_week  = 0
            llprevmth  = .FALSE.
            llprevyr   = .FALSE.
         ENDIF
         iyear  = nyear  - COUNT((/llprevyr /))
         imonth = nmonth - COUNT((/llprevmth/)) + 12 * COUNT((/llprevyr /))
         iday   = nday   + nmonth_len(nmonth-1) * COUNT((/llprevmth/)) - isec_week / NINT(rday)
      ENDIF

      ! build the new filename if not climatological data
      clname=TRIM(sdjf%clrootname)
      !
      ! note that sdjf%ln_clim is is only acting on the presence of the year in the file name
      IF( .NOT. sdjf%ln_clim ) THEN   
                                         WRITE(clname, '(a,"_y",i4.4)' ) TRIM( sdjf%clrootname ), iyear    ! add year
         IF( sdjf%cltype /= 'yearly' )   WRITE(clname, '(a,"m" ,i2.2)' ) TRIM( clname          ), imonth   ! add month
      ELSE
         ! build the new filename if climatological data
         IF( sdjf%cltype /= 'yearly' )   WRITE(clname, '(a,"_m",i2.2)' ) TRIM( sdjf%clrootname ), imonth   ! add month
      ENDIF
      IF( sdjf%cltype == 'daily' .OR. sdjf%cltype(1:4) == 'week' ) &
            &                            WRITE(clname, '(a,"d" ,i2.2)' ) TRIM( clname          ), iday     ! add day
      !
      IF( TRIM(clname) /= TRIM(sdjf%clname) .OR. sdjf%num == 0 ) THEN   ! new file to be open 
         !
         sdjf%clname = TRIM(clname)
         IF( sdjf%num /= 0 )   CALL iom_close( sdjf%num )   ! close file if already open
         CALL iom_open( sdjf%clname, sdjf%num, ldstop = ldstop, ldiof =  LEN(TRIM(sdjf%wgtname)) > 0 )
         !
         ! find the last record to be read -> update sdjf%nreclast
         indexyr = iyear - nyear + 1
         zyear_len = REAL(nyear_len( indexyr ), wp)
         SELECT CASE ( indexyr )
         CASE ( 0 )   ;   zmonth_len = 31.   ! previous year -> imonth = 12
         CASE ( 1 )   ;   zmonth_len = REAL(nmonth_len(imonth), wp)
         CASE ( 2 )   ;   zmonth_len = 31.   ! next     year -> imonth = 1
         END SELECT
         !
         ! last record to be read in the current file
         IF    ( sdjf%freqh == -12. ) THEN                 ;   sdjf%nreclast = 1    !  yearly mean
         ELSEIF( sdjf%freqh ==  -1. ) THEN                                          ! monthly mean
            IF(     sdjf%cltype      == 'monthly' ) THEN   ;   sdjf%nreclast = 1
            ELSE                                           ;   sdjf%nreclast = 12
            ENDIF
         ELSE                                                                       ! higher frequency mean (in hours)
            IF(     sdjf%cltype      == 'monthly' ) THEN   ;   sdjf%nreclast = NINT( 24. * zmonth_len / sdjf%freqh )
            ELSEIF( sdjf%cltype(1:4) == 'week'    ) THEN   ;   sdjf%nreclast = NINT( 24. * 7.         / sdjf%freqh )
            ELSEIF( sdjf%cltype      == 'daily'   ) THEN   ;   sdjf%nreclast = NINT( 24.              / sdjf%freqh )
            ELSE                                           ;   sdjf%nreclast = NINT( 24. * zyear_len  / sdjf%freqh )
            ENDIF
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE fld_clopn


   SUBROUTINE fld_fill( sdf, sdf_n, cdir, cdcaller, cdtitle, cdnam, knoprint )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_fill  ***
      !!
      !! ** Purpose :   fill the data structure (sdf) with the associated information 
      !!              read in namelist (sdf_n) and control print
      !!----------------------------------------------------------------------
      TYPE(FLD)  , DIMENSION(:)          , INTENT(inout) ::   sdf        ! structure of input fields (file informations, fields read)
      TYPE(FLD_N), DIMENSION(:)          , INTENT(in   ) ::   sdf_n      ! array of namelist information structures
      CHARACTER(len=*)                   , INTENT(in   ) ::   cdir       ! Root directory for location of flx files
      CHARACTER(len=*)                   , INTENT(in   ) ::   cdcaller   ! name of the calling routine
      CHARACTER(len=*)                   , INTENT(in   ) ::   cdtitle    ! description of the calling routine 
      CHARACTER(len=*)                   , INTENT(in   ) ::   cdnam      ! name of the namelist from which sdf_n comes
      INTEGER                  , OPTIONAL, INTENT(in   ) ::   knoprint   ! no calling routine information printed
      !
      INTEGER  ::   jf   ! dummy indices
      !!---------------------------------------------------------------------
      !
      DO jf = 1, SIZE(sdf)
         sdf(jf)%clrootname = sdf_n(jf)%clname
         IF( TRIM(sdf_n(jf)%clname) /= 'NOT USED' )   sdf(jf)%clrootname = TRIM( cdir )//sdf(jf)%clrootname
         sdf(jf)%clname     = "not yet defined"
         sdf(jf)%freqh      = sdf_n(jf)%freqh
         sdf(jf)%clvar      = sdf_n(jf)%clvar
         sdf(jf)%ln_tint    = sdf_n(jf)%ln_tint
         sdf(jf)%ln_clim    = sdf_n(jf)%ln_clim
         sdf(jf)%cltype     = sdf_n(jf)%cltype
         sdf(jf)%num        = -1
         sdf(jf)%wgtname    = " "
         IF( LEN( TRIM(sdf_n(jf)%wname) ) > 0 )   sdf(jf)%wgtname = TRIM( cdir )//sdf_n(jf)%wname
         sdf(jf)%lsmname = " "
         IF( LEN( TRIM(sdf_n(jf)%lname) ) > 0 )   sdf(jf)%lsmname = TRIM( cdir )//sdf_n(jf)%lname
         sdf(jf)%vcomp      = sdf_n(jf)%vcomp
         sdf(jf)%rotn(:)    = .TRUE.   ! pretend to be rotated -> won't try to rotate data before the first call to fld_get
         IF( sdf(jf)%cltype(1:4) == 'week' .AND. nn_leapy == 0  )   &
            &   CALL ctl_stop('fld_clopn: weekly file ('//TRIM(sdf(jf)%clrootname)//') needs nn_leapy = 1')
         IF( sdf(jf)%cltype(1:4) == 'week' .AND. sdf(jf)%ln_clim )   &
            &   CALL ctl_stop('fld_clopn: weekly file ('//TRIM(sdf(jf)%clrootname)//') needs ln_clim = .FALSE.')
         sdf(jf)%nreclast   = -1 ! Set to non zero default value to avoid errors, is updated to meaningful value during fld_clopn
         sdf(jf)%igrd       = 0
         sdf(jf)%ibdy       = 0
         sdf(jf)%imap       => NULL()
         sdf(jf)%ltotvel    = .FALSE.
         sdf(jf)%lzint      = .FALSE.
      END DO
      !
      IF(lwp) THEN      ! control print
         WRITE(numout,*)
         IF( .NOT.PRESENT( knoprint) ) THEN
            WRITE(numout,*) TRIM( cdcaller )//' : '//TRIM( cdtitle )
            WRITE(numout,*) (/ ('~', jf = 1, LEN_TRIM( cdcaller ) ) /)
         ENDIF
         WRITE(numout,*) '   fld_fill : fill data structure with information from namelist '//TRIM( cdnam )
         WRITE(numout,*) '   ~~~~~~~~'
         WRITE(numout,*) '      list of files and frequency (>0: in hours ; <0 in months)'
         DO jf = 1, SIZE(sdf)
            WRITE(numout,*) '      root filename: '  , TRIM( sdf(jf)%clrootname ), '   variable name: ', TRIM( sdf(jf)%clvar )
            WRITE(numout,*) '         frequency: '      ,       sdf(jf)%freqh       ,   &
               &                  '   time interp: '    ,       sdf(jf)%ln_tint     ,   &
               &                  '   climatology: '    ,       sdf(jf)%ln_clim
            WRITE(numout,*) '         weights: '        , TRIM( sdf(jf)%wgtname    ),   &
               &                  '   pairing: '        , TRIM( sdf(jf)%vcomp      ),   &
               &                  '   data type: '      ,       sdf(jf)%cltype      ,   &
               &                  '   land/sea mask:'   , TRIM( sdf(jf)%lsmname    )
            call flush(numout)
         END DO
      ENDIF
      !
   END SUBROUTINE fld_fill


   SUBROUTINE wgt_list( sd, kwgt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE wgt_list  ***
      !!
      !! ** Purpose :   search array of WGTs and find a weights file entry,
      !!                or return a new one adding it to the end if new entry.
      !!                the weights data is read in and restructured (fld_weight)
      !!----------------------------------------------------------------------
      TYPE( FLD ), INTENT(in   ) ::   sd        ! field with name of weights file
      INTEGER    , INTENT(inout) ::   kwgt      ! index of weights
      !
      INTEGER ::   kw, nestid   ! local integer
      LOGICAL ::   found        ! local logical
      !!----------------------------------------------------------------------
      !
      !! search down linked list 
      !! weights filename is either present or we hit the end of the list
      found = .FALSE.
      !
      !! because agrif nest part of filenames are now added in iom_open
      !! to distinguish between weights files on the different grids, need to track
      !! nest number explicitly
      nestid = 0
#if defined key_agrif
      nestid = Agrif_Fixed()
#endif
      DO kw = 1, nxt_wgt-1
         IF( TRIM(ref_wgts(kw)%wgtname) == TRIM(sd%wgtname) .AND. &
             ref_wgts(kw)%nestid == nestid) THEN
            kwgt = kw
            found = .TRUE.
            EXIT
         ENDIF
      END DO
      IF( .NOT.found ) THEN
         kwgt = nxt_wgt
         CALL fld_weight( sd )
      ENDIF
      !
   END SUBROUTINE wgt_list


   SUBROUTINE wgt_print( )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE wgt_print  ***
      !!
      !! ** Purpose :   print the list of known weights
      !!----------------------------------------------------------------------
      INTEGER ::   kw   !
      !!----------------------------------------------------------------------
      !
      DO kw = 1, nxt_wgt-1
         WRITE(numout,*) 'weight file:  ',TRIM(ref_wgts(kw)%wgtname)
         WRITE(numout,*) '      ddims:  ',ref_wgts(kw)%ddims(1),ref_wgts(kw)%ddims(2)
         WRITE(numout,*) '     numwgt:  ',ref_wgts(kw)%numwgt
         WRITE(numout,*) '     jpiwgt:  ',ref_wgts(kw)%jpiwgt
         WRITE(numout,*) '     jpjwgt:  ',ref_wgts(kw)%jpjwgt
         WRITE(numout,*) '    botleft:  ',ref_wgts(kw)%botleft
         WRITE(numout,*) '   topright:  ',ref_wgts(kw)%topright
         IF( ref_wgts(kw)%cyclic ) THEN
            WRITE(numout,*) '       cyclical'
            IF( ref_wgts(kw)%overlap > 0 ) WRITE(numout,*) '              with overlap of ', ref_wgts(kw)%overlap
         ELSE
            WRITE(numout,*) '       not cyclical'
         ENDIF
         IF( ASSOCIATED(ref_wgts(kw)%data_wgt) )  WRITE(numout,*) '       allocated'
      END DO
      !
   END SUBROUTINE wgt_print


   SUBROUTINE fld_weight( sd )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_weight  ***
      !!
      !! ** Purpose :   create a new WGT structure and fill in data from file,
      !!              restructuring as required
      !!----------------------------------------------------------------------
      TYPE( FLD ), INTENT(in) ::   sd   ! field with name of weights file
      !!
      INTEGER ::   jn         ! dummy loop indices
      INTEGER ::   inum       ! local logical unit
      INTEGER ::   id         ! local variable id
      INTEGER ::   ipk        ! local vertical dimension
      INTEGER ::   zwrap      ! local integer
      LOGICAL ::   cyclical   ! 
      CHARACTER (len=5) ::   aname   !
      INTEGER , DIMENSION(:), ALLOCATABLE ::   ddims
      INTEGER,  DIMENSION(jpi,jpj) ::   data_src
      REAL(wp), DIMENSION(jpi,jpj) ::   data_tmp
      !!----------------------------------------------------------------------
      !
      IF( nxt_wgt > tot_wgts ) THEN
        CALL ctl_stop("fld_weight: weights array size exceeded, increase tot_wgts")
      ENDIF
      !
      !! new weights file entry, add in extra information
      !! a weights file represents a 2D grid of a certain shape, so we assume that the current
      !! input data file is representative of all other files to be opened and processed with the
      !! current weights file

      !! open input data file (non-model grid)
      CALL iom_open( sd%clname, inum, ldiof =  LEN(TRIM(sd%wgtname)) > 0 )

      !! get dimensions
      IF ( SIZE(sd%fnow, 3) > 1 ) THEN
         ALLOCATE( ddims(4) )
      ELSE
         ALLOCATE( ddims(3) )
      ENDIF
      id = iom_varid( inum, sd%clvar, ddims )

      !! close it
      CALL iom_close( inum )

      !! now open the weights file

      CALL iom_open ( sd%wgtname, inum )   ! interpolation weights
      IF ( inum > 0 ) THEN

         !! determine whether we have an east-west cyclic grid
         !! from global attribute called "ew_wrap" in the weights file
         !! note that if not found, iom_getatt returns -999 and cyclic with no overlap is assumed
         !! since this is the most common forcing configuration

         CALL iom_getatt(inum, 'ew_wrap', zwrap)
         IF( zwrap >= 0 ) THEN
            cyclical = .TRUE.
         ELSE IF( zwrap == -999 ) THEN
            cyclical = .TRUE.
            zwrap = 0
         ELSE
            cyclical = .FALSE.
         ENDIF

         ref_wgts(nxt_wgt)%ddims(1) = ddims(1)
         ref_wgts(nxt_wgt)%ddims(2) = ddims(2)
         ref_wgts(nxt_wgt)%wgtname = sd%wgtname
         ref_wgts(nxt_wgt)%overlap = zwrap
         ref_wgts(nxt_wgt)%cyclic = cyclical
         ref_wgts(nxt_wgt)%nestid = 0
#if defined key_agrif
         ref_wgts(nxt_wgt)%nestid = Agrif_Fixed()
#endif
         !! weights file is stored as a set of weights (wgt01->wgt04 or wgt01->wgt16)
         !! for each weight wgtNN there is an integer array srcNN which gives the point in
         !! the input data grid which is to be multiplied by the weight
         !! they are both arrays on the model grid so the result of the multiplication is
         !! added into an output array on the model grid as a running sum

         !! two possible cases: bilinear (4 weights) or bicubic (16 weights)
         id = iom_varid(inum, 'src05', ldstop=.FALSE.)
         IF( id <= 0) THEN
            ref_wgts(nxt_wgt)%numwgt = 4
         ELSE
            ref_wgts(nxt_wgt)%numwgt = 16
         ENDIF

         ALLOCATE( ref_wgts(nxt_wgt)%data_jpi(jpi,jpj,4) )
         ALLOCATE( ref_wgts(nxt_wgt)%data_jpj(jpi,jpj,4) )
         ALLOCATE( ref_wgts(nxt_wgt)%data_wgt(jpi,jpj,ref_wgts(nxt_wgt)%numwgt) )

         DO jn = 1,4
            aname = ' '
            WRITE(aname,'(a3,i2.2)') 'src',jn
            data_tmp(:,:) = 0
            CALL iom_get ( inum, jpdom_data, aname, data_tmp(:,:) )
            data_src(:,:) = INT(data_tmp(:,:))
            ref_wgts(nxt_wgt)%data_jpj(:,:,jn) = 1 + (data_src(:,:)-1) / ref_wgts(nxt_wgt)%ddims(1)
            ref_wgts(nxt_wgt)%data_jpi(:,:,jn) = data_src(:,:) - ref_wgts(nxt_wgt)%ddims(1)*(ref_wgts(nxt_wgt)%data_jpj(:,:,jn)-1)
         END DO

         DO jn = 1, ref_wgts(nxt_wgt)%numwgt
            aname = ' '
            WRITE(aname,'(a3,i2.2)') 'wgt',jn
            ref_wgts(nxt_wgt)%data_wgt(:,:,jn) = 0.0
            CALL iom_get ( inum, jpdom_data, aname, ref_wgts(nxt_wgt)%data_wgt(:,:,jn) )
         END DO
         CALL iom_close (inum)
 
         ! find min and max indices in grid
         ref_wgts(nxt_wgt)%botleft(1) = MINVAL(ref_wgts(nxt_wgt)%data_jpi(:,:,:))
         ref_wgts(nxt_wgt)%botleft(2) = MINVAL(ref_wgts(nxt_wgt)%data_jpj(:,:,:))
         ref_wgts(nxt_wgt)%topright(1) = MAXVAL(ref_wgts(nxt_wgt)%data_jpi(:,:,:))
         ref_wgts(nxt_wgt)%topright(2) = MAXVAL(ref_wgts(nxt_wgt)%data_jpj(:,:,:))

         ! and therefore dimensions of the input box
         ref_wgts(nxt_wgt)%jpiwgt = ref_wgts(nxt_wgt)%topright(1) - ref_wgts(nxt_wgt)%botleft(1) + 1
         ref_wgts(nxt_wgt)%jpjwgt = ref_wgts(nxt_wgt)%topright(2) - ref_wgts(nxt_wgt)%botleft(2) + 1

         ! shift indexing of source grid
         ref_wgts(nxt_wgt)%data_jpi(:,:,:) = ref_wgts(nxt_wgt)%data_jpi(:,:,:) - ref_wgts(nxt_wgt)%botleft(1) + 1
         ref_wgts(nxt_wgt)%data_jpj(:,:,:) = ref_wgts(nxt_wgt)%data_jpj(:,:,:) - ref_wgts(nxt_wgt)%botleft(2) + 1

         ! create input grid, give it a halo to allow gradient calculations
         ! SA: +3 stencil is a patch to avoid out-of-bound computation in some configuration. 
         ! a more robust solution will be given in next release
         ipk =  SIZE(sd%fnow, 3)
         ALLOCATE( ref_wgts(nxt_wgt)%fly_dta(ref_wgts(nxt_wgt)%jpiwgt+3, ref_wgts(nxt_wgt)%jpjwgt+3 ,ipk) )
         IF( ref_wgts(nxt_wgt)%cyclic ) ALLOCATE( ref_wgts(nxt_wgt)%col(1,ref_wgts(nxt_wgt)%jpjwgt+3,ipk) )
         !
         nxt_wgt = nxt_wgt + 1
         !
      ELSE 
         CALL ctl_stop( '    fld_weight : unable to read the file ' )
      ENDIF

      DEALLOCATE (ddims )
      !
   END SUBROUTINE fld_weight


   SUBROUTINE apply_seaoverland( clmaskfile, zfieldo, jpi1_lsm, jpi2_lsm, jpj1_lsm,   &
      &                          jpj2_lsm, itmpi, itmpj, itmpz, rec1_lsm, recn_lsm )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE apply_seaoverland  ***
      !!
      !! ** Purpose :   avoid spurious fluxes in coastal or near-coastal areas
      !!                due to the wrong usage of "land" values from the coarse
      !!                atmospheric model when spatial interpolation is required
      !!      D. Delrosso INGV          
      !!---------------------------------------------------------------------- 
      INTEGER,                   INTENT(in   ) :: itmpi,itmpj,itmpz                    ! lengths
      INTEGER,                   INTENT(in   ) :: jpi1_lsm,jpi2_lsm,jpj1_lsm,jpj2_lsm  ! temporary indices
      INTEGER, DIMENSION(3),     INTENT(in   ) :: rec1_lsm,recn_lsm                    ! temporary arrays for start and length
      REAL(wp),DIMENSION (:,:,:),INTENT(inout) :: zfieldo                              ! input/output array for seaoverland application
      CHARACTER (len=100),       INTENT(in   ) :: clmaskfile                           ! land/sea mask file name
      !
      INTEGER :: inum,jni,jnj,jnz,jc   ! local indices
      REAL(wp),DIMENSION (:,:,:),ALLOCATABLE :: zslmec1             ! local array for land point detection
      REAL(wp),DIMENSION (:,:),  ALLOCATABLE :: zfieldn   ! array of forcing field with undeff for land points
      REAL(wp),DIMENSION (:,:),  ALLOCATABLE :: zfield    ! array of forcing field
      !!---------------------------------------------------------------------
      !
      ALLOCATE ( zslmec1(itmpi,itmpj,itmpz), zfieldn(itmpi,itmpj), zfield(itmpi,itmpj) )
      !
      ! Retrieve the land sea mask data
      CALL iom_open( clmaskfile, inum )
      SELECT CASE( SIZE(zfieldo(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,:),3) )
      CASE(1)
         CALL iom_get( inum, jpdom_unknown, 'LSM', zslmec1(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,1), 1, rec1_lsm, recn_lsm)
      CASE DEFAULT
         CALL iom_get( inum, jpdom_unknown, 'LSM', zslmec1(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,:), 1, rec1_lsm, recn_lsm)
      END SELECT
      CALL iom_close( inum )
      !
      DO jnz=1,rec1_lsm(3)             !! Loop over k dimension
         !
         DO jni = 1, itmpi                               !! copy the original field into a tmp array
            DO jnj = 1, itmpj                            !! substituting undeff over land points
               zfieldn(jni,jnj) = zfieldo(jni,jnj,jnz)
               IF( zslmec1(jni,jnj,jnz) == 1. )   zfieldn(jni,jnj) = undeff_lsm
            END DO
         END DO
         !
         CALL seaoverland( zfieldn, itmpi, itmpj, zfield )
         DO jc = 1, nn_lsm
            CALL seaoverland( zfield, itmpi, itmpj, zfield )
         END DO
         !
         !   Check for Undeff and substitute original values
         IF( ANY(zfield==undeff_lsm) ) THEN
            DO jni = 1, itmpi
               DO jnj = 1, itmpj
                  IF( zfield(jni,jnj)==undeff_lsm )   zfield(jni,jnj) = zfieldo(jni,jnj,jnz)
               END DO
            END DO
         ENDIF
         !
         zfieldo(:,:,jnz) = zfield(:,:)
         !
      END DO                           !! End Loop over k dimension
      !
      DEALLOCATE ( zslmec1, zfieldn, zfield )
      !
   END SUBROUTINE apply_seaoverland 


   SUBROUTINE seaoverland( zfieldn, ileni, ilenj, zfield )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE seaoverland  ***
      !!
      !! ** Purpose :   create shifted matrices for seaoverland application  
      !!      D. Delrosso INGV
      !!---------------------------------------------------------------------- 
      INTEGER                      , INTENT(in   ) :: ileni,ilenj   ! lengths 
      REAL, DIMENSION (ileni,ilenj), INTENT(in   ) :: zfieldn       ! array of forcing field with undeff for land points
      REAL, DIMENSION (ileni,ilenj), INTENT(  out) :: zfield        ! array of forcing field
      !
      REAL   , DIMENSION (ileni,ilenj)   :: zmat1, zmat2, zmat3, zmat4  ! local arrays 
      REAL   , DIMENSION (ileni,ilenj)   :: zmat5, zmat6, zmat7, zmat8  !   -     - 
      REAL   , DIMENSION (ileni,ilenj)   :: zlsm2d                      !   -     - 
      REAL   , DIMENSION (ileni,ilenj,8) :: zlsm3d                      !   -     -
      LOGICAL, DIMENSION (ileni,ilenj,8) :: ll_msknan3d                 ! logical mask for undeff detection
      LOGICAL, DIMENSION (ileni,ilenj)   :: ll_msknan2d                 ! logical mask for undeff detection
      !!---------------------------------------------------------------------- 
      zmat8 = eoshift( zfieldn , SHIFT=-1 , BOUNDARY = (/zfieldn(:,1)/)     , DIM=2 )
      zmat1 = eoshift( zmat8   , SHIFT=-1 , BOUNDARY = (/zmat8(1,:)/)       , DIM=1 )
      zmat2 = eoshift( zfieldn , SHIFT=-1 , BOUNDARY = (/zfieldn(1,:)/)     , DIM=1 )
      zmat4 = eoshift( zfieldn , SHIFT= 1 , BOUNDARY = (/zfieldn(:,ilenj)/) , DIM=2 )
      zmat3 = eoshift( zmat4   , SHIFT=-1 , BOUNDARY = (/zmat4(1,:)/)       , DIM=1 )
      zmat5 = eoshift( zmat4   , SHIFT= 1 , BOUNDARY = (/zmat4(ileni,:)/)   , DIM=1 )
      zmat6 = eoshift( zfieldn , SHIFT= 1 , BOUNDARY = (/zfieldn(ileni,:)/) , DIM=1 )
      zmat7 = eoshift( zmat8   , SHIFT= 1 , BOUNDARY = (/zmat8(ileni,:)/)   , DIM=1 )
      !
      zlsm3d  = RESHAPE( (/ zmat1, zmat2, zmat3, zmat4, zmat5, zmat6, zmat7, zmat8 /), (/ ileni, ilenj, 8 /))
      ll_msknan3d = .NOT.( zlsm3d  == undeff_lsm )
      ll_msknan2d = .NOT.( zfieldn == undeff_lsm )  ! FALSE where is Undeff (land)
      zlsm2d = SUM( zlsm3d, 3 , ll_msknan3d ) / MAX( 1 , COUNT( ll_msknan3d , 3 ) )
      WHERE( COUNT( ll_msknan3d , 3 ) == 0._wp )   zlsm2d = undeff_lsm
      zfield = MERGE( zfieldn, zlsm2d, ll_msknan2d )
      !
   END SUBROUTINE seaoverland


   SUBROUTINE fld_interp( num, clvar, kw, kk, dta,  &
                          &         nrec, lsmfile)      
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE fld_interp  ***
      !!
      !! ** Purpose :   apply weights to input gridded data to create data
      !!                on model grid
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   num     ! stream number
      CHARACTER(LEN=*)          , INTENT(in   ) ::   clvar   ! variable name
      INTEGER                   , INTENT(in   ) ::   kw      ! weights number
      INTEGER                   , INTENT(in   ) ::   kk      ! vertical dimension of kk
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   dta     ! output field on model grid
      INTEGER                   , INTENT(in   ) ::   nrec    ! record number to read (ie time slice)
      CHARACTER(LEN=*)          , INTENT(in   ) ::   lsmfile ! land sea mask file name
      !
      INTEGER, DIMENSION(3) ::   rec1, recn           ! temporary arrays for start and length
      INTEGER, DIMENSION(3) ::   rec1_lsm, recn_lsm   ! temporary arrays for start and length in case of seaoverland
      INTEGER ::   ii_lsm1,ii_lsm2,ij_lsm1,ij_lsm2    ! temporary indices
      INTEGER ::   jk, jn, jm, jir, jjr               ! loop counters
      INTEGER ::   ni, nj                             ! lengths
      INTEGER ::   jpimin,jpiwid                      ! temporary indices
      INTEGER ::   jpimin_lsm,jpiwid_lsm              ! temporary indices
      INTEGER ::   jpjmin,jpjwid                      ! temporary indices
      INTEGER ::   jpjmin_lsm,jpjwid_lsm              ! temporary indices
      INTEGER ::   jpi1,jpi2,jpj1,jpj2                ! temporary indices
      INTEGER ::   jpi1_lsm,jpi2_lsm,jpj1_lsm,jpj2_lsm   ! temporary indices
      INTEGER ::   itmpi,itmpj,itmpz                     ! lengths
      REAL(wp),DIMENSION(:,:,:), ALLOCATABLE ::   ztmp_fly_dta                 ! local array of values on input grid     
      !!----------------------------------------------------------------------
      !
      !! for weighted interpolation we have weights at four corners of a box surrounding 
      !! a model grid point, each weight is multiplied by a grid value (bilinear case)
      !! or by a grid value and gradients at the corner point (bicubic case) 
      !! so we need to have a 4 by 4 subgrid surrounding each model point to cover both cases

      !! sub grid from non-model input grid which encloses all grid points in this nemo process
      jpimin = ref_wgts(kw)%botleft(1)
      jpjmin = ref_wgts(kw)%botleft(2)
      jpiwid = ref_wgts(kw)%jpiwgt
      jpjwid = ref_wgts(kw)%jpjwgt

      !! when reading in, expand this sub-grid by one halo point all the way round for calculating gradients
      rec1(1) = MAX( jpimin-1, 1 )
      rec1(2) = MAX( jpjmin-1, 1 )
      rec1(3) = 1
      recn(1) = MIN( jpiwid+2, ref_wgts(kw)%ddims(1)-rec1(1)+1 )
      recn(2) = MIN( jpjwid+2, ref_wgts(kw)%ddims(2)-rec1(2)+1 )
      recn(3) = kk

      !! where we need to put it in the non-nemo grid fly_dta
      !! note that jpi1 and jpj1 only differ from 1 when jpimin and jpjmin are 1
      !! (ie at the extreme west or south of the whole input grid) and similarly for jpi2 and jpj2
      jpi1 = 2 + rec1(1) - jpimin
      jpj1 = 2 + rec1(2) - jpjmin
      jpi2 = jpi1 + recn(1) - 1
      jpj2 = jpj1 + recn(2) - 1


      IF( LEN( TRIM(lsmfile) ) > 0 ) THEN
      !! indeces for ztmp_fly_dta
      ! --------------------------
         rec1_lsm(1)=MAX(rec1(1)-nn_lsm,1)  ! starting index for enlarged external data, x direction
         rec1_lsm(2)=MAX(rec1(2)-nn_lsm,1)  ! starting index for enlarged external data, y direction
         rec1_lsm(3) = 1                    ! vertical dimension
         recn_lsm(1)=MIN(rec1(1)-rec1_lsm(1)+recn(1)+nn_lsm,ref_wgts(kw)%ddims(1)-rec1_lsm(1)) ! n points in x direction
         recn_lsm(2)=MIN(rec1(2)-rec1_lsm(2)+recn(2)+nn_lsm,ref_wgts(kw)%ddims(2)-rec1_lsm(2)) ! n points in y direction
         recn_lsm(3) = kk                   ! number of vertical levels in the input file

      !  Avoid out of bound
         jpimin_lsm = MAX( rec1_lsm(1)+1, 1 )
         jpjmin_lsm = MAX( rec1_lsm(2)+1, 1 )
         jpiwid_lsm = MIN( recn_lsm(1)-2,ref_wgts(kw)%ddims(1)-rec1(1)+1)
         jpjwid_lsm = MIN( recn_lsm(2)-2,ref_wgts(kw)%ddims(2)-rec1(2)+1)

         jpi1_lsm = 2+rec1_lsm(1)-jpimin_lsm
         jpj1_lsm = 2+rec1_lsm(2)-jpjmin_lsm
         jpi2_lsm = jpi1_lsm + recn_lsm(1) - 1
         jpj2_lsm = jpj1_lsm + recn_lsm(2) - 1


         itmpi=jpi2_lsm-jpi1_lsm+1
         itmpj=jpj2_lsm-jpj1_lsm+1
         itmpz=kk
         ALLOCATE(ztmp_fly_dta(itmpi,itmpj,itmpz))
         ztmp_fly_dta(:,:,:) = 0.0
         SELECT CASE( SIZE(ztmp_fly_dta(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,:),3) )
         CASE(1)
              CALL iom_get( num, jpdom_unknown, clvar, ztmp_fly_dta(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,1),   &
                 &                                                                nrec, rec1_lsm, recn_lsm)
         CASE DEFAULT
              CALL iom_get( num, jpdom_unknown, clvar, ztmp_fly_dta(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,:),   &
                 &                                                                nrec, rec1_lsm, recn_lsm)
         END SELECT
         CALL apply_seaoverland(lsmfile,ztmp_fly_dta(jpi1_lsm:jpi2_lsm,jpj1_lsm:jpj2_lsm,:),                  &
                 &                                      jpi1_lsm,jpi2_lsm,jpj1_lsm,jpj2_lsm,                  &
                 &                                      itmpi,itmpj,itmpz,rec1_lsm,recn_lsm)


         ! Relative indeces for remapping
         ii_lsm1 = (rec1(1)-rec1_lsm(1))+1
         ii_lsm2 = (ii_lsm1+recn(1))-1
         ij_lsm1 = (rec1(2)-rec1_lsm(2))+1
         ij_lsm2 = (ij_lsm1+recn(2))-1

         ref_wgts(kw)%fly_dta(:,:,:) = 0.0
         ref_wgts(kw)%fly_dta(jpi1:jpi2,jpj1:jpj2,:) = ztmp_fly_dta(ii_lsm1:ii_lsm2,ij_lsm1:ij_lsm2,:)
         DEALLOCATE(ztmp_fly_dta)

      ELSE
         
         ref_wgts(kw)%fly_dta(:,:,:) = 0.0
         SELECT CASE( SIZE(ref_wgts(kw)%fly_dta(jpi1:jpi2,jpj1:jpj2,:),3) )
         CASE(1)
              CALL iom_get( num, jpdom_unknown, clvar, ref_wgts(kw)%fly_dta(jpi1:jpi2,jpj1:jpj2,1), nrec, rec1, recn)
         CASE DEFAULT
              CALL iom_get( num, jpdom_unknown, clvar, ref_wgts(kw)%fly_dta(jpi1:jpi2,jpj1:jpj2,:), nrec, rec1, recn)
         END SELECT 
      ENDIF
      

      !! first four weights common to both bilinear and bicubic
      !! data_jpi, data_jpj have already been shifted to (1,1) corresponding to botleft
      !! note that we have to offset by 1 into fly_dta array because of halo
      dta(:,:,:) = 0.0
      DO jk = 1,4
        DO jn = 1, jpj
          DO jm = 1,jpi
            ni = ref_wgts(kw)%data_jpi(jm,jn,jk)
            nj = ref_wgts(kw)%data_jpj(jm,jn,jk)
            dta(jm,jn,:) = dta(jm,jn,:) + ref_wgts(kw)%data_wgt(jm,jn,jk) * ref_wgts(kw)%fly_dta(ni+1,nj+1,:)
          END DO
        END DO
      END DO

      IF (ref_wgts(kw)%numwgt .EQ. 16) THEN

        !! fix up halo points that we couldnt read from file
        IF( jpi1 == 2 ) THEN
           ref_wgts(kw)%fly_dta(jpi1-1,:,:) = ref_wgts(kw)%fly_dta(jpi1,:,:)
        ENDIF
        IF( jpi2 + jpimin - 1 == ref_wgts(kw)%ddims(1)+1 ) THEN
           ref_wgts(kw)%fly_dta(jpi2+1,:,:) = ref_wgts(kw)%fly_dta(jpi2,:,:)
        ENDIF
        IF( jpj1 == 2 ) THEN
           ref_wgts(kw)%fly_dta(:,jpj1-1,:) = ref_wgts(kw)%fly_dta(:,jpj1,:)
        ENDIF
        IF( jpj2 + jpjmin - 1 == ref_wgts(kw)%ddims(2)+1 .AND. jpj2 .lt. jpjwid+2 ) THEN
           ref_wgts(kw)%fly_dta(:,jpj2+1,:) = 2.0*ref_wgts(kw)%fly_dta(:,jpj2,:) - ref_wgts(kw)%fly_dta(:,jpj2-1,:)
        ENDIF

        !! if data grid is cyclic we can do better on east-west edges
        !! but have to allow for whether first and last columns are coincident
        IF( ref_wgts(kw)%cyclic ) THEN
           rec1(2) = MAX( jpjmin-1, 1 )
           recn(1) = 1
           recn(2) = MIN( jpjwid+2, ref_wgts(kw)%ddims(2)-rec1(2)+1 )
           jpj1 = 2 + rec1(2) - jpjmin
           jpj2 = jpj1 + recn(2) - 1
           IF( jpi1 == 2 ) THEN
              rec1(1) = ref_wgts(kw)%ddims(1) - ref_wgts(kw)%overlap
              SELECT CASE( SIZE( ref_wgts(kw)%col(:,jpj1:jpj2,:),3) )
              CASE(1)
                   CALL iom_get( num, jpdom_unknown, clvar, ref_wgts(kw)%col(:,jpj1:jpj2,1), nrec, rec1, recn)
              CASE DEFAULT
                   CALL iom_get( num, jpdom_unknown, clvar, ref_wgts(kw)%col(:,jpj1:jpj2,:), nrec, rec1, recn)
              END SELECT      
              ref_wgts(kw)%fly_dta(jpi1-1,jpj1:jpj2,:) = ref_wgts(kw)%col(1,jpj1:jpj2,:)
           ENDIF
           IF( jpi2 + jpimin - 1 == ref_wgts(kw)%ddims(1)+1 ) THEN
              rec1(1) = 1 + ref_wgts(kw)%overlap
              SELECT CASE( SIZE( ref_wgts(kw)%col(:,jpj1:jpj2,:),3) )
              CASE(1)
                   CALL iom_get( num, jpdom_unknown, clvar, ref_wgts(kw)%col(:,jpj1:jpj2,1), nrec, rec1, recn)
              CASE DEFAULT
                   CALL iom_get( num, jpdom_unknown, clvar, ref_wgts(kw)%col(:,jpj1:jpj2,:), nrec, rec1, recn)
              END SELECT
              ref_wgts(kw)%fly_dta(jpi2+1,jpj1:jpj2,:) = ref_wgts(kw)%col(1,jpj1:jpj2,:)
           ENDIF
        ENDIF

        ! gradient in the i direction
        DO jk = 1,4
          DO jn = 1, jpj
            DO jm = 1,jpi
              ni = ref_wgts(kw)%data_jpi(jm,jn,jk)
              nj = ref_wgts(kw)%data_jpj(jm,jn,jk)
              dta(jm,jn,:) = dta(jm,jn,:) + ref_wgts(kw)%data_wgt(jm,jn,jk+4) * 0.5 *         &
                               (ref_wgts(kw)%fly_dta(ni+2,nj+1,:) - ref_wgts(kw)%fly_dta(ni,nj+1,:))
            END DO
          END DO
        END DO

        ! gradient in the j direction
        DO jk = 1,4
          DO jn = 1, jpj
            DO jm = 1,jpi
              ni = ref_wgts(kw)%data_jpi(jm,jn,jk)
              nj = ref_wgts(kw)%data_jpj(jm,jn,jk)
              dta(jm,jn,:) = dta(jm,jn,:) + ref_wgts(kw)%data_wgt(jm,jn,jk+8) * 0.5 *         &
                               (ref_wgts(kw)%fly_dta(ni+1,nj+2,:) - ref_wgts(kw)%fly_dta(ni+1,nj,:))
            END DO
          END DO
        END DO

         ! gradient in the ij direction
         DO jk = 1,4
            DO jn = 1, jpj
               DO jm = 1,jpi
                  ni = ref_wgts(kw)%data_jpi(jm,jn,jk)
                  nj = ref_wgts(kw)%data_jpj(jm,jn,jk)
                  dta(jm,jn,:) = dta(jm,jn,:) + ref_wgts(kw)%data_wgt(jm,jn,jk+12) * 0.25 * ( &
                               (ref_wgts(kw)%fly_dta(ni+2,nj+2,:) - ref_wgts(kw)%fly_dta(ni  ,nj+2,:)) -   &
                               (ref_wgts(kw)%fly_dta(ni+2,nj  ,:) - ref_wgts(kw)%fly_dta(ni  ,nj  ,:)))
               END DO
            END DO
         END DO
         !
      END IF
      !
   END SUBROUTINE fld_interp


   FUNCTION ksec_week( cdday )
      !!---------------------------------------------------------------------
      !!                    ***  FUNCTION kshift_week *** 
      !!
      !! ** Purpose :   return the first 3 letters of the first day of the weekly file
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in)   ::   cdday   ! first 3 letters of the first day of the weekly file
      !!
      INTEGER                        ::   ksec_week      ! output variable
      INTEGER                        ::   ijul, ishift   ! local integer
      CHARACTER(len=3),DIMENSION(7)  ::   cl_week 
      !!----------------------------------------------------------------------
      cl_week = (/"sun","sat","fri","thu","wed","tue","mon"/)
      DO ijul = 1, 7
         IF( cl_week(ijul) == TRIM(cdday) ) EXIT
      END DO
      IF( ijul .GT. 7 )   CALL ctl_stop( 'ksec_week: wrong day for sdjf%cltype(6:8): '//TRIM(cdday) )
      !
      ishift = ijul * NINT(rday)
      ! 
      ksec_week = nsec_week + ishift
      ksec_week = MOD( ksec_week, 7*NINT(rday) )
      ! 
   END FUNCTION ksec_week

   !!======================================================================
END MODULE fldread
