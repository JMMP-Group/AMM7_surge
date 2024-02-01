MODULE icbini
   !!======================================================================
   !!                       ***  MODULE  icbini  ***
   !! Icebergs:  initialise variables for iceberg tracking
   !!======================================================================
   !! History :   -   !  2010-01  (T. Martin & A. Adcroft)  Original code
   !!            3.3  !  2011-03  (G. Madec)  Part conversion to NEMO form ; Removal of mapping from another grid
   !!             -   !  2011-04  (S. Alderson)  Split into separate modules ; Restore restart routines
   !!             -   !  2011-05  (S. Alderson)  generate_test_icebergs restored ; new forcing arrays with extra halo ;
   !!             -   !                          north fold exchange arrays added
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   icb_init     : initialise icebergs
   !!   icb_ini_gen  : generate test icebergs
   !!   icb_nam      : read iceberg namelist
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean domain
   USE in_out_manager ! IO routines and numout in particular
   USE lib_mpp        ! mpi library and lk_mpp in particular
   USE sbc_oce        ! ocean  : surface boundary condition
   USE sbc_ice        ! sea-ice: surface boundary condition
   USE iom            ! IOM library
   USE fldread        ! field read
   USE lbclnk         ! lateral boundary condition - MPP link
   !
   USE icb_oce        ! define iceberg arrays
   USE icbutl         ! iceberg utility routines
   USE icbrst         ! iceberg restart routines
   USE icbtrj         ! iceberg trajectory I/O routines
   USE icbdia         ! iceberg budget routines

   IMPLICIT NONE
   PRIVATE

   PUBLIC   icb_init  ! routine called in nemogcm.F90 module

   CHARACTER(len=100)                                 ::   cn_dir = './'   !: Root directory for location of icb files
   TYPE(FLD_N)                                        ::   sn_icb          !: information about the calving file to be read
   TYPE(FLD), PUBLIC, ALLOCATABLE     , DIMENSION(:)  ::   sf_icb          !: structure: file information, fields read
                                                                           !: used in icbini and icbstp
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE icb_init( pdt, kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_init  ***
      !!
      !! ** Purpose :   iceberg initialization.
      !!
      !! ** Method  : - read the iceberg namelist
      !!              - find non-overlapping processor interior since we can only
      !!                have one instance of a particular iceberg
      !!              - calculate the destinations for north fold exchanges
      !!              - setup either test icebergs or calving file
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::   pdt   ! iceberg time-step (rdt*nn_fsbc)
      INTEGER , INTENT(in) ::   kt    ! time step number
      !
      INTEGER ::   ji, jj, jn               ! dummy loop indices
      INTEGER ::   i1, i2, i3               ! local integers
      INTEGER ::   ii, inum, ivar           !   -       -
      INTEGER ::   istat1, istat2, istat3   !   -       -
      CHARACTER(len=300) ::   cl_sdist      ! local character
      !!----------------------------------------------------------------------
      !
      CALL icb_nam               ! Read and print namelist parameters
      !
      IF( .NOT. ln_icebergs )   RETURN

      !                          ! allocate gridded fields
      IF( icb_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'icb_alloc : unable to allocate arrays' )
      !
      !                          ! initialised variable with extra haloes to zero
      uo_e(:,:) = 0._wp   ;   vo_e(:,:) = 0._wp   ;
      ua_e(:,:) = 0._wp   ;   va_e(:,:) = 0._wp   ;
      ff_e(:,:) = 0._wp   ;   tt_e(:,:) = 0._wp   ;
      fr_e(:,:) = 0._wp   ;   ss_e(:,:) = 0._wp   ;
#if defined key_si3
      hi_e(:,:) = 0._wp   ;
      ui_e(:,:) = 0._wp   ;   vi_e(:,:) = 0._wp   ;
#endif
      ssh_e(:,:) = 0._wp  ; 
      !
      !                          ! open ascii output file or files for iceberg status information
      !                          ! note that we choose to do this on all processors since we cannot
      !                          ! predict where icebergs will be ahead of time
      IF( nn_verbose_level > 0) THEN
         CALL ctl_opn( numicb, 'icebergs.stat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
      ENDIF

      ! set parameters (mostly from namelist)
      !
      berg_dt         = pdt
      first_width (:) = SQRT(  rn_initial_mass(:) / ( rn_LoW_ratio * rn_rho_bergs * rn_initial_thickness(:) )  )
      first_length(:) = rn_LoW_ratio * first_width(:)

      berg_grid%calving      (:,:)   = 0._wp
      berg_grid%calving_hflx (:,:)   = 0._wp
      berg_grid%stored_heat  (:,:)   = 0._wp
      berg_grid%floating_melt(:,:)   = 0._wp
      berg_grid%maxclass     (:,:)   = nclasses
      berg_grid%stored_ice   (:,:,:) = 0._wp
      berg_grid%tmp          (:,:)   = 0._wp
      src_calving            (:,:)   = 0._wp
      src_calving_hflx       (:,:)   = 0._wp

      !                          ! domain for icebergs
      IF( lk_mpp .AND. jpni == 1 )   CALL ctl_stop( 'icbinit: having ONE processor in x currently does not work' )
      ! NB: the issue here is simply that cyclic east-west boundary condition have not been coded in mpp case
      ! for the north fold we work out which points communicate by asking
      ! lbc_lnk to pass processor number (valid even in single processor case)
      ! borrow src_calving arrays for this
      !
      ! pack i and j together using a scaling of a power of 10
      nicbpack = 10000
      IF( jpiglo >= nicbpack )   CALL ctl_stop( 'icbini: processor index packing failure' )
      nicbfldproc(:) = -1

      DO jj = 1, jpj
         DO ji = 1, jpi
            src_calving_hflx(ji,jj) = narea
            src_calving     (ji,jj) = nicbpack * mjg(jj) + mig(ji)
         END DO
      END DO
      CALL lbc_lnk( 'icbini', src_calving_hflx, 'T', 1._wp )
      CALL lbc_lnk( 'icbini', src_calving     , 'T', 1._wp )

      ! work out interior of processor from exchange array
      ! first entry with narea for this processor is left hand interior index
      ! last  entry                               is right hand interior index
      jj = nlcj/2
      nicbdi = -1
      nicbei = -1
      DO ji = 1, jpi
         i3 = INT( src_calving(ji,jj) )
         i2 = INT( i3/nicbpack )
         i1 = i3 - i2*nicbpack
         i3 = INT( src_calving_hflx(ji,jj) )
         IF( i1 == mig(ji) .AND. i3 == narea ) THEN
            IF( nicbdi < 0 ) THEN   ;   nicbdi = ji
            ELSE                    ;   nicbei = ji
            ENDIF
         ENDIF
      END DO
      !
      ! repeat for j direction
      ji = nlci/2
      nicbdj = -1
      nicbej = -1
      DO jj = 1, jpj
         i3 = INT( src_calving(ji,jj) )
         i2 = INT( i3/nicbpack )
         i1 = i3 - i2*nicbpack
         i3 = INT( src_calving_hflx(ji,jj) )
         IF( i2 == mjg(jj) .AND. i3 == narea ) THEN
            IF( nicbdj < 0 ) THEN   ;   nicbdj = jj
            ELSE                    ;   nicbej = jj
            ENDIF
         ENDIF
      END DO
      !   
      ! special for east-west boundary exchange we save the destination index
      i1 = MAX( nicbdi-1, 1)
      i3 = INT( src_calving(i1,nlcj/2) )
      jj = INT( i3/nicbpack )
      ricb_left = REAL( i3 - nicbpack*jj, wp )
      i1 = MIN( nicbei+1, jpi )
      i3 = INT( src_calving(i1,nlcj/2) )
      jj = INT( i3/nicbpack )
      ricb_right = REAL( i3 - nicbpack*jj, wp )
      
      ! north fold
      IF( npolj > 0 ) THEN
         !
         ! icebergs in row nicbej+1 get passed across fold
         nicbfldpts(:)  = INT( src_calving(:,nicbej+1) )
         nicbflddest(:) = INT( src_calving_hflx(:,nicbej+1) )
         !
         ! work out list of unique processors to talk to
         ! pack them into a fixed size array where empty slots are marked by a -1
         DO ji = nicbdi, nicbei
            ii = nicbflddest(ji)
            IF( ii .GT. 0 ) THEN     ! Needed because land suppression can mean
                                     ! that unused points are not set in edge haloes
               DO jn = 1, jpni
                  ! work along array until we find an empty slot
                  IF( nicbfldproc(jn) == -1 ) THEN
                     nicbfldproc(jn) = ii
                     EXIT                             !!gm EXIT should be avoided: use DO WHILE expression instead
                  ENDIF
                  ! before we find an empty slot, we may find processor number is already here so we exit
                  IF( nicbfldproc(jn) == ii ) EXIT
               END DO
            ENDIF
         END DO
      ENDIF
      !
      IF( nn_verbose_level > 0) THEN
         WRITE(numicb,*) 'processor ', narea
         WRITE(numicb,*) 'jpi, jpj   ', jpi, jpj
         WRITE(numicb,*) 'nldi, nlei ', nldi, nlei
         WRITE(numicb,*) 'nldj, nlej ', nldj, nlej
         WRITE(numicb,*) 'berg i interior ', nicbdi, nicbei
         WRITE(numicb,*) 'berg j interior ', nicbdj, nicbej
         WRITE(numicb,*) 'berg left       ', ricb_left
         WRITE(numicb,*) 'berg right      ', ricb_right
         jj = nlcj/2
         WRITE(numicb,*) "central j line:"
         WRITE(numicb,*) "i processor"
         WRITE(numicb,*) (INT(src_calving_hflx(ji,jj)), ji=1,jpi)
         WRITE(numicb,*) "i point"
         WRITE(numicb,*) (INT(src_calving(ji,jj)), ji=1,jpi)
         ji = nlci/2
         WRITE(numicb,*) "central i line:"
         WRITE(numicb,*) "j processor"
         WRITE(numicb,*) (INT(src_calving_hflx(ji,jj)), jj=1,jpj)
         WRITE(numicb,*) "j point"
         WRITE(numicb,*) (INT(src_calving(ji,jj)), jj=1,jpj)
         IF( npolj > 0 ) THEN
            WRITE(numicb,*) 'north fold destination points '
            WRITE(numicb,*) nicbfldpts
            WRITE(numicb,*) 'north fold destination procs  '
            WRITE(numicb,*) nicbflddest
            WRITE(numicb,*) 'north fold destination proclist  '
            WRITE(numicb,*) nicbfldproc
         ENDIF
         CALL flush(numicb)
      ENDIF
      
      src_calving     (:,:) = 0._wp
      src_calving_hflx(:,:) = 0._wp

      ! definition of extended surface masked needed by icb_bilin_h
      tmask_e(:,:) = 0._wp   ;   tmask_e(1:jpi,1:jpj) = tmask(:,:,1)
      umask_e(:,:) = 0._wp   ;   umask_e(1:jpi,1:jpj) = umask(:,:,1)
      vmask_e(:,:) = 0._wp   ;   vmask_e(1:jpi,1:jpj) = vmask(:,:,1)
      CALL lbc_lnk_icb( 'icbini', tmask_e, 'T', +1._wp, 1, 1 )
      CALL lbc_lnk_icb( 'icbini', umask_e, 'U', +1._wp, 1, 1 )
      CALL lbc_lnk_icb( 'icbini', vmask_e, 'V', +1._wp, 1, 1 )
      !
      ! assign each new iceberg with a unique number constructed from the processor number
      ! and incremented by the total number of processors
      num_bergs(:) = 0
      num_bergs(1) = narea - jpnij

      ! when not generating test icebergs we need to setup calving file
      IF( nn_test_icebergs < 0 .OR. ln_use_calving ) THEN
         !
         ! maximum distribution class array does not change in time so read it once
         cl_sdist = TRIM( cn_dir )//TRIM( sn_icb%clname )
         CALL iom_open ( cl_sdist, inum )                              ! open file
         ivar = iom_varid( inum, 'maxclass', ldstop=.FALSE. )
         IF( ivar > 0 ) THEN
            CALL iom_get  ( inum, jpdom_data, 'maxclass', src_calving )   ! read the max distribution array
            berg_grid%maxclass(:,:) = INT( src_calving )
            src_calving(:,:) = 0._wp
         ENDIF
         CALL iom_close( inum )                                     ! close file
         !
         IF( nn_verbose_level > 0) THEN
            WRITE(numicb,*)
            WRITE(numicb,*) '          calving read in a file'
         ENDIF
         ALLOCATE( sf_icb(1), STAT=istat1 )         ! Create sf_icb structure (calving)
         ALLOCATE( sf_icb(1)%fnow(jpi,jpj,1), STAT=istat2 )
         ALLOCATE( sf_icb(1)%fdta(jpi,jpj,1,2), STAT=istat3 )
         IF( istat1+istat2+istat3 > 0 ) THEN
            CALL ctl_stop( 'sbc_icb: unable to allocate sf_icb structure' )   ;   RETURN
         ENDIF
         !                                          ! fill sf_icb with the namelist (sn_icb) and control print
         CALL fld_fill( sf_icb, (/ sn_icb /), cn_dir, 'icb_init', 'read calving data', 'namicb' )
         !
      ENDIF

      IF( .NOT.ln_rstart ) THEN
         IF( nn_test_icebergs > 0 )   CALL icb_ini_gen()
      ELSE
         IF( nn_test_icebergs > 0 ) THEN
            CALL icb_ini_gen()
         ELSE
            CALL icb_rst_read()
            l_restarted_bergs = .TRUE.
         ENDIF
      ENDIF
      !
      IF( nn_sample_rate .GT. 0 ) CALL icb_trj_init( nitend )
      !
      CALL icb_dia_init()
      !
      IF( nn_verbose_level >= 2 )   CALL icb_utl_print('icb_init, initial status', nit000-1)
      !
   END SUBROUTINE icb_init


   SUBROUTINE icb_ini_gen()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_ini_gen  ***
      !!
      !! ** Purpose :   iceberg generation
      !!
      !! ** Method  : - at each grid point of the test box supplied in the namelist
      !!                generate an iceberg in one class determined by the value of
      !!                parameter nn_test_icebergs
      !!----------------------------------------------------------------------
      INTEGER                         ::   ji, jj, ibergs
      TYPE(iceberg)                   ::   localberg ! NOT a pointer but an actual local variable
      TYPE(point)                     ::   localpt
      INTEGER                         ::   iyr, imon, iday, ihr, imin, isec
      INTEGER                         ::   iberg
      !!----------------------------------------------------------------------

      ! For convenience
      iberg = nn_test_icebergs

      ! call get_date(Time, iyr, imon, iday, ihr, imin, isec)
      ! Convert nemo time variables from dom_oce into local versions
      iyr  = nyear
      imon = nmonth
      iday = nday
      ihr = INT(nsec_day/3600)
      imin = INT((nsec_day-ihr*3600)/60)
      isec = nsec_day - ihr*3600 - imin*60

      ! no overlap for icebergs since we want only one instance of each across the whole domain
      ! so restrict area of interest
      ! use tmask here because tmask_i has been doctored on one side of the north fold line

      DO jj = nicbdj, nicbej
         DO ji = nicbdi, nicbei
            IF( tmask(ji,jj,1) > 0._wp        .AND.                                       &
                rn_test_box(1) < glamt(ji,jj) .AND. glamt(ji,jj) < rn_test_box(2) .AND.   &
                rn_test_box(3) < gphit(ji,jj) .AND. gphit(ji,jj) < rn_test_box(4) ) THEN
               localberg%mass_scaling = rn_mass_scaling(iberg)
               localpt%xi = REAL( mig(ji), wp )
               localpt%yj = REAL( mjg(jj), wp )
               localpt%lon = icb_utl_bilin(glamt, localpt%xi, localpt%yj, 'T' )
               localpt%lat = icb_utl_bilin(gphit, localpt%xi, localpt%yj, 'T' )
               localpt%mass      = rn_initial_mass     (iberg)
               localpt%thickness = rn_initial_thickness(iberg)
               localpt%width  = first_width (iberg)
               localpt%length = first_length(iberg)
               localpt%year = iyr
               localpt%day = REAL(iday,wp)+(REAL(ihr,wp)+REAL(imin,wp)/60._wp)/24._wp
               localpt%mass_of_bits = 0._wp
               localpt%heat_density = 0._wp
               localpt%uvel = 0._wp
               localpt%vvel = 0._wp
               CALL icb_utl_incr()
               localberg%number(:) = num_bergs(:)
               call icb_utl_add(localberg, localpt)
            ENDIF
         END DO
      END DO
      !
      ibergs = icb_utl_count()
      CALL mpp_sum('icbini', ibergs)
      IF( nn_verbose_level > 0) THEN
         WRITE(numicb,'(a,i6,a)') 'diamonds, icb_ini_gen: ',ibergs,' were generated'
      ENDIF
      !
   END SUBROUTINE icb_ini_gen


   SUBROUTINE icb_nam
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE icb_nam  ***
      !!
      !! ** Purpose :   read iceberg namelist and print the variables.
      !!
      !! ** input   : - namberg namelist
      !!----------------------------------------------------------------------
      INTEGER  ::   jn      ! dummy loop indices
      INTEGER  ::   ios     ! Local integer output status for namelist read
      REAL(wp) ::   zfact   ! local scalar
      !
      NAMELIST/namberg/ ln_icebergs    , ln_bergdia     , nn_sample_rate      , rn_initial_mass      ,   &
         &              rn_distribution, rn_mass_scaling, rn_initial_thickness, nn_verbose_write     ,   &
         &              rn_rho_bergs   , rn_LoW_ratio   , nn_verbose_level    , ln_operator_splitting,   &
         &              rn_bits_erosion_fraction        , rn_sicn_shift       , ln_passive_mode      ,   &
         &              nn_test_icebergs                , rn_test_box         , ln_use_calving       ,   &
         &              rn_speed_limit , cn_dir, sn_icb
      !!----------------------------------------------------------------------

#if defined key_agrif
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'icb_nam : AGRIF is not compatible with namelist namberg :  '
         WRITE(numout,*) '~~~~~~~   definition of rn_initial_mass(nclasses) with nclasses as PARAMETER '
         WRITE(numout,*)
         WRITE(numout,*) '   ==>>>   force  NO icebergs used. The namelist namberg is not read'
      ENDIF
      ln_icebergs = .false.      
      RETURN
#else
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'icb_nam : iceberg initialization through namberg namelist read'
         WRITE(numout,*) '~~~~~~~~ '
      ENDIF
#endif   
      !                             !==  read namelist  ==!
      REWIND( numnam_ref )              ! Namelist namberg in reference namelist : Iceberg parameters
      READ  ( numnam_ref, namberg, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namberg in reference namelist' )
      REWIND( numnam_cfg )              ! Namelist namberg in configuration namelist : Iceberg parameters
      READ  ( numnam_cfg, namberg, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namberg in configuration namelist' )
      IF(lwm) WRITE ( numond, namberg )
      !
      IF(lwp) WRITE(numout,*)
      IF( ln_icebergs ) THEN
         IF(lwp) WRITE(numout,*) '   ==>>>   icebergs are used'
      ELSE
         IF(lwp) WRITE(numout,*) '   ==>>>   No icebergs used'
         RETURN
      ENDIF
      !
      IF( nn_test_icebergs > nclasses ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   Resetting of nn_test_icebergs to ', nclasses
         nn_test_icebergs = nclasses
      ENDIF
      !
      IF( nn_test_icebergs < 0 .AND. .NOT. ln_use_calving ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   Resetting ln_use_calving to .true. since we are not using test icebergs'
         ln_use_calving = .true.
      ENDIF
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'icb_nam : iceberg initialization through namberg namelist read'
         WRITE(numout,*) '~~~~~~~~ '
         WRITE(numout,*) '   Calculate budgets                                            ln_bergdia       = ', ln_bergdia
         WRITE(numout,*) '   Period between sampling of position for trajectory storage   nn_sample_rate = ', nn_sample_rate
         WRITE(numout,*) '   Mass thresholds between iceberg classes (kg)                 rn_initial_mass     ='
         DO jn = 1, nclasses
            WRITE(numout,'(a,f15.2)') '                                                                ', rn_initial_mass(jn)
         ENDDO
         WRITE(numout,*) '   Fraction of calving to apply to this class (non-dim)         rn_distribution     ='
         DO jn = 1, nclasses
            WRITE(numout,'(a,f10.4)') '                                                                ', rn_distribution(jn)
         END DO
         WRITE(numout,*) '   Ratio between effective and real iceberg mass (non-dim)      rn_mass_scaling     = '
         DO jn = 1, nclasses
            WRITE(numout,'(a,f10.2)') '                                                                ', rn_mass_scaling(jn)
         END DO
         WRITE(numout,*) '   Total thickness of newly calved bergs (m)                    rn_initial_thickness = '
         DO jn = 1, nclasses
            WRITE(numout,'(a,f10.2)') '                                                                ', rn_initial_thickness(jn)
         END DO
         WRITE(numout,*) '   Timesteps between verbose messages                           nn_verbose_write    = ', nn_verbose_write

         WRITE(numout,*) '   Density of icebergs                           rn_rho_bergs  = ', rn_rho_bergs
         WRITE(numout,*) '   Initial ratio L/W for newly calved icebergs   rn_LoW_ratio  = ', rn_LoW_ratio
         WRITE(numout,*) '   Turn on more verbose output                          level  = ', nn_verbose_level
         WRITE(numout,*) '   Use first order operator splitting for thermodynamics    ',   &
            &                    'use_operator_splitting = ', ln_operator_splitting
         WRITE(numout,*) '   Fraction of erosion melt flux to divert to bergy bits    ',   &
            &                    'bits_erosion_fraction = ', rn_bits_erosion_fraction

         WRITE(numout,*) '   Shift of sea-ice concentration in erosion flux modulation ',   &
            &                    '(0<sicn_shift<1)    rn_sicn_shift  = ', rn_sicn_shift
         WRITE(numout,*) '   Do not add freshwater flux from icebergs to ocean                ',   &
            &                    '                  passive_mode            = ', ln_passive_mode
         WRITE(numout,*) '   Create icebergs in absence of a restart file   nn_test_icebergs  = ', nn_test_icebergs
         WRITE(numout,*) '                   in lon/lat box                                   = ', rn_test_box
         WRITE(numout,*) '   Use calving data even if nn_test_icebergs > 0    ln_use_calving  = ', ln_use_calving
         WRITE(numout,*) '   CFL speed limit for a berg            speed_limit                = ', rn_speed_limit
         WRITE(numout,*) '   Writing Iceberg status information to icebergs.stat file        '
      ENDIF
      !
      ! ensure that the sum of berg input distribution is equal to one
      zfact = SUM( rn_distribution )
      IF( zfact /= 1._wp .AND. 0_wp /= zfact ) THEN
         rn_distribution(:) = rn_distribution(:) / zfact
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) '      ==>>> CAUTION:    sum of berg input distribution = ', zfact
            WRITE(numout,*) '            *******     redistribution has been rescaled'
            WRITE(numout,*) '                        updated berg distribution is :'
            DO jn = 1, nclasses
               WRITE(numout,'(a,f10.4)') '                                   ',rn_distribution(jn)
            END DO
         ENDIF
      ENDIF
      IF( MINVAL( rn_distribution(:) ) < 0._wp ) THEN
         CALL ctl_stop( 'icb_nam: a negative rn_distribution value encountered ==>> change your namelist namberg' )
      ENDIF
      !
   END SUBROUTINE icb_nam

   !!======================================================================
END MODULE icbini
