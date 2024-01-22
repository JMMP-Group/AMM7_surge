MODULE obs_write
   !!======================================================================
   !!                       ***  MODULE obs_write   ***
   !! Observation diagnosticss: Write observation related diagnostics
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   obs_wri_prof   : Write profile observations in feedback format
   !!   obs_wri_surf   : Write surface observations in feedback format
   !!   obs_wri_stats  : Print basic statistics on the data being written out
   !!----------------------------------------------------------------------

   !! * Modules used
   USE par_kind, ONLY : &   ! Precision variables
      & wp
   USE in_out_manager       ! I/O manager
   USE dom_oce              ! Ocean space and time domain variables
   USE obs_types            ! Observation type integer to character translation
   USE julian, ONLY : &         ! Julian date routines
      & greg2jul
   USE obs_utils, ONLY : &  ! Observation operator utility functions
      & chkerr
   USE obs_profiles_def     ! Type definitions for profiles
   USE obs_surf_def         ! Type defintions for surface observations
   USE obs_fbm              ! Observation feedback I/O
   USE obs_grid             ! Grid tools
   USE obs_conv             ! Conversion between units
   USE obs_const
   USE obs_mpp              ! MPP support routines for observation diagnostics
   USE lib_mpp		    ! MPP routines

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE
   PUBLIC obs_wri_prof, &    ! Write profile observation files
      &   obs_wri_surf, &    ! Write surface observation files
      &   obswriinfo
   
   TYPE obswriinfo
      INTEGER :: inum
      INTEGER, POINTER, DIMENSION(:) :: ipoint
      CHARACTER(len=ilenname), POINTER, DIMENSION(:) :: cdname
      CHARACTER(len=ilenlong), POINTER, DIMENSION(:,:) :: cdlong
      CHARACTER(len=ilenunit), POINTER, DIMENSION(:,:) :: cdunit
   END TYPE obswriinfo

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE obs_wri_prof( profdata, padd, pext )
      !!-----------------------------------------------------------------------
      !!
      !!                     *** ROUTINE obs_wri_prof  ***
      !!
      !! ** Purpose : Write profile feedback files
      !!
      !! ** Method  : NetCDF
      !! 
      !! ** Action  :
      !!
      !! History :
      !!      ! 06-04  (A. Vidard) Original
      !!      ! 06-04  (A. Vidard) Reformatted
      !!      ! 06-10  (A. Weaver) Cleanup
      !!      ! 07-01  (K. Mogensen) Use profile data types
      !!      ! 07-03  (K. Mogensen) General handling of profiles
      !!      ! 09-01  (K. Mogensen) New feedback format
      !!      ! 15-02  (M. Martin) Combined routine for writing profiles
      !!-----------------------------------------------------------------------

      !! * Arguments
      TYPE(obs_prof), INTENT(INOUT) :: profdata      ! Full set of profile data
      TYPE(obswriinfo), OPTIONAL :: padd             ! Additional info for each variable
      TYPE(obswriinfo), OPTIONAL :: pext             ! Extra info

      !! * Local declarations
      TYPE(obfbdata) :: fbdata
      CHARACTER(LEN=40) :: clfname
      CHARACTER(LEN=10) :: clfiletype
      INTEGER :: ilevel
      INTEGER :: jvar
      INTEGER :: jo
      INTEGER :: jk
      INTEGER :: ik
      INTEGER :: ja
      INTEGER :: je
      INTEGER :: iadd
      INTEGER :: iext
      REAL(wp) :: zpres

      IF ( PRESENT( padd ) ) THEN
         iadd = padd%inum
      ELSE
         iadd = 0
      ENDIF

      IF ( PRESENT( pext ) ) THEN
         iext = pext%inum
      ELSE
         iext = 0
      ENDIF

      CALL init_obfbdata( fbdata )

      ! Find maximum level
      ilevel = 0
      DO jvar = 1, 2
         ilevel = MAX( ilevel, MAXVAL( profdata%var(jvar)%nvlidx(:) ) )
      END DO

      SELECT CASE ( TRIM(profdata%cvars(1)) )
      CASE('POTM')

         clfiletype='profb'
         CALL alloc_obfbdata( fbdata, 2, profdata%nprof, ilevel, &
            &                 1 + iadd, 1 + iext, .TRUE. )
         fbdata%cname(1)      = profdata%cvars(1)
         fbdata%cname(2)      = profdata%cvars(2)
         fbdata%coblong(1)    = 'Potential temperature'
         fbdata%coblong(2)    = 'Practical salinity'
         fbdata%cobunit(1)    = 'Degrees centigrade'
         fbdata%cobunit(2)    = 'PSU'
         fbdata%cextname(1)   = 'TEMP'
         fbdata%cextlong(1)   = 'Insitu temperature'
         fbdata%cextunit(1)   = 'Degrees centigrade'
         fbdata%caddlong(1,1) = 'Model interpolated potential temperature'
         fbdata%caddlong(1,2) = 'Model interpolated practical salinity'
         fbdata%caddunit(1,1) = 'Degrees centigrade'
         fbdata%caddunit(1,2) = 'PSU'
         fbdata%cgrid(:)      = 'T'
         DO je = 1, iext
            fbdata%cextname(1+je) = pext%cdname(je)
            fbdata%cextlong(1+je) = pext%cdlong(je,1)
            fbdata%cextunit(1+je) = pext%cdunit(je,1)
         END DO
         DO ja = 1, iadd
            fbdata%caddname(1+ja) = padd%cdname(ja)
            DO jvar = 1, 2
               fbdata%caddlong(1+ja,jvar) = padd%cdlong(ja,jvar)
               fbdata%caddunit(1+ja,jvar) = padd%cdunit(ja,jvar)
            END DO
         END DO

      CASE('UVEL')

         clfiletype='velfb'
         CALL alloc_obfbdata( fbdata, 2, profdata%nprof, ilevel, 1, 0, .TRUE. )
         fbdata%cname(1)      = profdata%cvars(1)
         fbdata%cname(2)      = profdata%cvars(2)
         fbdata%coblong(1)    = 'Zonal velocity'
         fbdata%coblong(2)    = 'Meridional velocity'
         fbdata%cobunit(1)    = 'm/s'
         fbdata%cobunit(2)    = 'm/s'
         DO je = 1, iext
            fbdata%cextname(je) = pext%cdname(je)
            fbdata%cextlong(je) = pext%cdlong(je,1)
            fbdata%cextunit(je) = pext%cdunit(je,1)
         END DO
         fbdata%caddlong(1,1) = 'Model interpolated zonal velocity'
         fbdata%caddlong(1,2) = 'Model interpolated meridional velocity'
         fbdata%caddunit(1,1) = 'm/s'
         fbdata%caddunit(1,2) = 'm/s'
         fbdata%cgrid(1)      = 'U' 
         fbdata%cgrid(2)      = 'V'
         DO ja = 1, iadd
            fbdata%caddname(1+ja) = padd%cdname(ja)
            fbdata%caddlong(1+ja,1) = padd%cdlong(ja,1)
            fbdata%caddunit(1+ja,1) = padd%cdunit(ja,1)
         END DO

      END SELECT

      fbdata%caddname(1)   = 'Hx'

      WRITE(clfname, FMT="(A,'_fdbk_',I4.4,'.nc')") TRIM(clfiletype), nproc

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*)'obs_wri_prof :'
         WRITE(numout,*)'~~~~~~~~~~~~~'
         WRITE(numout,*)'Writing '//TRIM(clfiletype)//' feedback file : ',TRIM(clfname)
      ENDIF

      ! Transform obs_prof data structure into obfb data structure
      fbdata%cdjuldref = '19500101000000'
      DO jo = 1, profdata%nprof
         fbdata%plam(jo)      = profdata%rlam(jo)
         fbdata%pphi(jo)      = profdata%rphi(jo)
         WRITE(fbdata%cdtyp(jo),'(I4)') profdata%ntyp(jo)
         fbdata%ivqc(jo,:)    = profdata%ivqc(jo,:)
         fbdata%ivqcf(:,jo,:) = profdata%ivqcf(:,jo,:)
         IF ( profdata%nqc(jo) > 255 ) THEN
            fbdata%ioqc(jo)    = IBSET(profdata%nqc(jo),2)
            fbdata%ioqcf(1,jo) = profdata%nqcf(1,jo)
            fbdata%ioqcf(2,jo) = profdata%nqc(jo)
         ELSE
            fbdata%ioqc(jo)    = profdata%nqc(jo)
            fbdata%ioqcf(:,jo) = profdata%nqcf(:,jo)
         ENDIF
         fbdata%ipqc(jo)      = profdata%ipqc(jo)
         fbdata%ipqcf(:,jo)   = profdata%ipqcf(:,jo)
         fbdata%itqc(jo)      = profdata%itqc(jo)
         fbdata%itqcf(:,jo)   = profdata%itqcf(:,jo)
         fbdata%cdwmo(jo)     = profdata%cwmo(jo)
         fbdata%kindex(jo)    = profdata%npfil(jo)
         DO jvar = 1, profdata%nvar
            IF (ln_grid_global) THEN
               fbdata%iobsi(jo,jvar) = profdata%mi(jo,jvar)
               fbdata%iobsj(jo,jvar) = profdata%mj(jo,jvar)
            ELSE
               fbdata%iobsi(jo,jvar) = mig(profdata%mi(jo,jvar))
               fbdata%iobsj(jo,jvar) = mjg(profdata%mj(jo,jvar))
            ENDIF
         END DO
         CALL greg2jul( 0, &
            &           profdata%nmin(jo), &
            &           profdata%nhou(jo), &
            &           profdata%nday(jo), &
            &           profdata%nmon(jo), &
            &           profdata%nyea(jo), &
            &           fbdata%ptim(jo),   &
            &           krefdate = 19500101 )
         ! Reform the profiles arrays for output
         DO jvar = 1, 2
            DO jk = profdata%npvsta(jo,jvar), profdata%npvend(jo,jvar)
               ik = profdata%var(jvar)%nvlidx(jk)
               fbdata%padd(ik,jo,1,jvar) = profdata%var(jvar)%vmod(jk)
               fbdata%pob(ik,jo,jvar)    = profdata%var(jvar)%vobs(jk)
               fbdata%pdep(ik,jo)        = profdata%var(jvar)%vdep(jk)
               fbdata%idqc(ik,jo)        = profdata%var(jvar)%idqc(jk)
               fbdata%idqcf(:,ik,jo)     = profdata%var(jvar)%idqcf(:,jk)
               IF ( profdata%var(jvar)%nvqc(jk) > 255 ) THEN
                  fbdata%ivlqc(ik,jo,jvar) = IBSET(profdata%var(jvar)%nvqc(jk),2)
                  fbdata%ivlqcf(1,ik,jo,jvar) = profdata%var(jvar)%nvqcf(1,jk)
!$AGRIF_DO_NOT_TREAT
                  fbdata%ivlqcf(2,ik,jo,jvar) = IAND(profdata%var(jvar)%nvqc(jk),b'0000000011111111')
!$AGRIF_END_DO_NOT_TREAT
               ELSE
                  fbdata%ivlqc(ik,jo,jvar) = profdata%var(jvar)%nvqc(jk)
                  fbdata%ivlqcf(:,ik,jo,jvar) = profdata%var(jvar)%nvqcf(:,jk)
               ENDIF
               fbdata%iobsk(ik,jo,jvar)  = profdata%var(jvar)%mvk(jk)
               DO ja = 1, iadd
                  fbdata%padd(ik,jo,1+ja,jvar) = &
                     & profdata%var(jvar)%vext(jk,padd%ipoint(ja))
               END DO
               DO je = 1, iext
                  fbdata%pext(ik,jo,1+je) = &
                     & profdata%var(jvar)%vext(jk,pext%ipoint(je))
               END DO
               IF ( ( jvar == 1 ) .AND. &
                  & ( TRIM(profdata%cvars(1)) == 'POTM' ) ) THEN
                  fbdata%pext(ik,jo,1) = profdata%var(jvar)%vext(jk,1)
               ENDIF 
            END DO
         END DO
      END DO

      IF ( TRIM(profdata%cvars(1)) == 'POTM' ) THEN
         ! Convert insitu temperature to potential temperature using the model
         ! salinity if no potential temperature
         DO jo = 1, fbdata%nobs
            IF ( fbdata%pphi(jo) < 9999.0 ) THEN
               DO jk = 1, fbdata%nlev
                  IF ( ( fbdata%pob(jk,jo,1) >= 9999.0 ) .AND. &
                     & ( fbdata%pdep(jk,jo) < 9999.0 ) .AND. &
                     & ( fbdata%padd(jk,jo,1,2) < 9999.0 ) .AND. &
                     & ( fbdata%pext(jk,jo,1) < 9999.0 ) ) THEN
                     zpres = dep_to_p( REAL(fbdata%pdep(jk,jo),wp), &
                        &              REAL(fbdata%pphi(jo),wp) )
                     fbdata%pob(jk,jo,1) = potemp( &
                        &                     REAL(fbdata%padd(jk,jo,1,2), wp),  &
                        &                     REAL(fbdata%pext(jk,jo,1), wp), &
                        &                     zpres, 0.0_wp )
                  ENDIF
               END DO
            ENDIF
         END DO
      ENDIF

      ! Write the obfbdata structure
      CALL write_obfbdata( clfname, fbdata )

      ! Output some basic statistics
      CALL obs_wri_stats( fbdata )

      CALL dealloc_obfbdata( fbdata )

   END SUBROUTINE obs_wri_prof

   SUBROUTINE obs_wri_surf( surfdata, padd, pext )
      !!-----------------------------------------------------------------------
      !!
      !!                     *** ROUTINE obs_wri_surf  ***
      !!
      !! ** Purpose : Write surface observation files
      !!
      !! ** Method  : NetCDF
      !! 
      !! ** Action  :
      !!
      !!      ! 07-03  (K. Mogensen) Original
      !!      ! 09-01  (K. Mogensen) New feedback format.
      !!      ! 15-02  (M. Martin) Combined surface writing routine.
      !!-----------------------------------------------------------------------

      !! * Modules used
      IMPLICIT NONE

      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) :: surfdata         ! Full set of surface data
      TYPE(obswriinfo), OPTIONAL :: padd               ! Additional info for each variable
      TYPE(obswriinfo), OPTIONAL :: pext               ! Extra info

      !! * Local declarations
      TYPE(obfbdata) :: fbdata
      CHARACTER(LEN=40) :: clfname         ! netCDF filename
      CHARACTER(LEN=10) :: clfiletype
      CHARACTER(LEN=12), PARAMETER :: cpname = 'obs_wri_surf'
      INTEGER :: jo
      INTEGER :: ja
      INTEGER :: je
      INTEGER :: iadd
      INTEGER :: iext

      IF ( PRESENT( padd ) ) THEN
         iadd = padd%inum
      ELSE
         iadd = 0
      ENDIF

      IF ( PRESENT( pext ) ) THEN
         iext = pext%inum
      ELSE
         iext = 0
      ENDIF

      CALL init_obfbdata( fbdata )

      SELECT CASE ( TRIM(surfdata%cvars(1)) )
      CASE('SLA')

         CALL alloc_obfbdata( fbdata, 1, surfdata%nsurf, 1, &
            &                 2 + iadd, 1 + iext, .TRUE. )

         clfiletype = 'slafb'
         fbdata%cname(1)      = surfdata%cvars(1)
         fbdata%coblong(1)    = 'Sea level anomaly'
         fbdata%cobunit(1)    = 'Metres'
         fbdata%cextname(1)   = 'MDT'
         fbdata%cextlong(1)   = 'Mean dynamic topography'
         fbdata%cextunit(1)   = 'Metres'
         DO je = 1, iext
            fbdata%cextname(je) = pext%cdname(je)
            fbdata%cextlong(je) = pext%cdlong(je,1)
            fbdata%cextunit(je) = pext%cdunit(je,1)
         END DO
         fbdata%caddlong(1,1) = 'Model interpolated SSH - MDT'
         fbdata%caddunit(1,1) = 'Metres' 
         fbdata%caddname(2)   = 'SSH'
         fbdata%caddlong(2,1) = 'Model Sea surface height'
         fbdata%caddunit(2,1) = 'Metres'
         fbdata%cgrid(1)      = 'T'
         DO ja = 1, iadd
            fbdata%caddname(2+ja) = padd%cdname(ja)
            fbdata%caddlong(2+ja,1) = padd%cdlong(ja,1)
            fbdata%caddunit(2+ja,1) = padd%cdunit(ja,1)
         END DO

      CASE('SST')

         CALL alloc_obfbdata( fbdata, 1, surfdata%nsurf, 1, &
            &                 1 + iadd, iext, .TRUE. )

         clfiletype = 'sstfb'
         fbdata%cname(1)      = surfdata%cvars(1)
         fbdata%coblong(1)    = 'Sea surface temperature'
         fbdata%cobunit(1)    = 'Degree centigrade'
         DO je = 1, iext
            fbdata%cextname(je) = pext%cdname(je)
            fbdata%cextlong(je) = pext%cdlong(je,1)
            fbdata%cextunit(je) = pext%cdunit(je,1)
         END DO
         fbdata%caddlong(1,1) = 'Model interpolated SST'
         fbdata%caddunit(1,1) = 'Degree centigrade'
         fbdata%cgrid(1)      = 'T'
         DO ja = 1, iadd
            fbdata%caddname(1+ja) = padd%cdname(ja)
            fbdata%caddlong(1+ja,1) = padd%cdlong(ja,1)
            fbdata%caddunit(1+ja,1) = padd%cdunit(ja,1)
         END DO

      CASE('ICECONC')

         CALL alloc_obfbdata( fbdata, 1, surfdata%nsurf, 1, &
            &                 1 + iadd, iext, .TRUE. )

         clfiletype = 'sicfb'
         fbdata%cname(1)      = surfdata%cvars(1)
         fbdata%coblong(1)    = 'Sea ice'
         fbdata%cobunit(1)    = 'Fraction'
         DO je = 1, iext
            fbdata%cextname(je) = pext%cdname(je)
            fbdata%cextlong(je) = pext%cdlong(je,1)
            fbdata%cextunit(je) = pext%cdunit(je,1)
         END DO
         fbdata%caddlong(1,1) = 'Model interpolated ICE'
         fbdata%caddunit(1,1) = 'Fraction'
         fbdata%cgrid(1)      = 'T'
         DO ja = 1, iadd
            fbdata%caddname(1+ja) = padd%cdname(ja)
            fbdata%caddlong(1+ja,1) = padd%cdlong(ja,1)
            fbdata%caddunit(1+ja,1) = padd%cdunit(ja,1)
         END DO

      CASE('SSS')

         CALL alloc_obfbdata( fbdata, 1, surfdata%nsurf, 1, &
            &                 1 + iadd, iext, .TRUE. )

         clfiletype = 'sssfb'
         fbdata%cname(1)      = surfdata%cvars(1)
         fbdata%coblong(1)    = 'Sea surface salinity'
         fbdata%cobunit(1)    = 'psu'
         DO je = 1, iext
            fbdata%cextname(je) = pext%cdname(je)
            fbdata%cextlong(je) = pext%cdlong(je,1)
            fbdata%cextunit(je) = pext%cdunit(je,1)
         END DO
         fbdata%caddlong(1,1) = 'Model interpolated SSS'
         fbdata%caddunit(1,1) = 'psu'
         fbdata%cgrid(1)      = 'T'
         DO ja = 1, iadd
            fbdata%caddname(1+ja) = padd%cdname(ja)
            fbdata%caddlong(1+ja,1) = padd%cdlong(ja,1)
            fbdata%caddunit(1+ja,1) = padd%cdunit(ja,1)
         END DO

      CASE DEFAULT

         CALL ctl_stop( 'Unknown observation type '//TRIM(surfdata%cvars(1))//' in obs_wri_surf' )

      END SELECT

      fbdata%caddname(1)   = 'Hx'

      WRITE(clfname, FMT="(A,'_fdbk_',I4.4,'.nc')") TRIM(clfiletype), nproc

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*)'obs_wri_surf :'
         WRITE(numout,*)'~~~~~~~~~~~~~'
         WRITE(numout,*)'Writing '//TRIM(surfdata%cvars(1))//' feedback file : ',TRIM(clfname)
      ENDIF

      ! Transform surf data structure into obfbdata structure
      fbdata%cdjuldref = '19500101000000'
      DO jo = 1, surfdata%nsurf
         fbdata%plam(jo)      = surfdata%rlam(jo)
         fbdata%pphi(jo)      = surfdata%rphi(jo)
         WRITE(fbdata%cdtyp(jo),'(I4)') surfdata%ntyp(jo)
         fbdata%ivqc(jo,:)    = 0
         fbdata%ivqcf(:,jo,:) = 0
         IF ( surfdata%nqc(jo) > 255 ) THEN
            fbdata%ioqc(jo)    = 4
            fbdata%ioqcf(1,jo) = 0
!$AGRIF_DO_NOT_TREAT
            fbdata%ioqcf(2,jo) = IAND(surfdata%nqc(jo),b'0000000011111111')
!$AGRIF_END_DO_NOT_TREAT
         ELSE
            fbdata%ioqc(jo)    = surfdata%nqc(jo)
            fbdata%ioqcf(:,jo) = 0
         ENDIF
         fbdata%ipqc(jo)      = 0
         fbdata%ipqcf(:,jo)   = 0
         fbdata%itqc(jo)      = 0
         fbdata%itqcf(:,jo)   = 0
         fbdata%cdwmo(jo)     = surfdata%cwmo(jo)
         fbdata%kindex(jo)    = surfdata%nsfil(jo)
         IF (ln_grid_global) THEN
            fbdata%iobsi(jo,1) = surfdata%mi(jo)
            fbdata%iobsj(jo,1) = surfdata%mj(jo)
         ELSE
            fbdata%iobsi(jo,1) = mig(surfdata%mi(jo))
            fbdata%iobsj(jo,1) = mjg(surfdata%mj(jo))
         ENDIF
         CALL greg2jul( 0, &
            &           surfdata%nmin(jo), &
            &           surfdata%nhou(jo), &
            &           surfdata%nday(jo), &
            &           surfdata%nmon(jo), &
            &           surfdata%nyea(jo), &
            &           fbdata%ptim(jo),   &
            &           krefdate = 19500101 )
         fbdata%padd(1,jo,1,1) = surfdata%rmod(jo,1)
         IF ( TRIM(surfdata%cvars(1)) == 'SLA' ) fbdata%padd(1,jo,2,1) = surfdata%rext(jo,1)
         fbdata%pob(1,jo,1)    = surfdata%robs(jo,1) 
         fbdata%pdep(1,jo)     = 0.0
         fbdata%idqc(1,jo)     = 0
         fbdata%idqcf(:,1,jo)  = 0
         IF ( surfdata%nqc(jo) > 255 ) THEN
            fbdata%ivqc(jo,1)       = 4
            fbdata%ivlqc(1,jo,1)    = 4
            fbdata%ivlqcf(1,1,jo,1) = 0
!$AGRIF_DO_NOT_TREAT
            fbdata%ivlqcf(2,1,jo,1) = IAND(surfdata%nqc(jo),b'0000000011111111')
!$AGRIF_END_DO_NOT_TREAT
         ELSE
            fbdata%ivqc(jo,1)       = surfdata%nqc(jo)
            fbdata%ivlqc(1,jo,1)    = surfdata%nqc(jo)
            fbdata%ivlqcf(:,1,jo,1) = 0
         ENDIF
         fbdata%iobsk(1,jo,1)  = 0
         IF ( TRIM(surfdata%cvars(1)) == 'SLA' ) fbdata%pext(1,jo,1) = surfdata%rext(jo,2)
         DO ja = 1, iadd
            fbdata%padd(1,jo,2+ja,1) = &
               & surfdata%rext(jo,padd%ipoint(ja))
         END DO
         DO je = 1, iext
            fbdata%pext(1,jo,1+je) = &
               & surfdata%rext(jo,pext%ipoint(je))
         END DO
      END DO

      ! Write the obfbdata structure
      CALL write_obfbdata( clfname, fbdata )

      ! Output some basic statistics
      CALL obs_wri_stats( fbdata )

      CALL dealloc_obfbdata( fbdata )

   END SUBROUTINE obs_wri_surf

   SUBROUTINE obs_wri_stats( fbdata )
      !!-----------------------------------------------------------------------
      !!
      !!                     *** ROUTINE obs_wri_stats  ***
      !!
      !! ** Purpose : Output some basic statistics of the data being written out
      !!
      !! ** Method  :
      !! 
      !! ** Action  :
      !!
      !!      ! 2014-08  (D. Lea) Initial version 
      !!-----------------------------------------------------------------------

      !! * Arguments
      TYPE(obfbdata) :: fbdata

      !! * Local declarations
      INTEGER :: jvar
      INTEGER :: jo
      INTEGER :: jk
      INTEGER :: inumgoodobs
      INTEGER :: inumgoodobsmpp
      REAL(wp) :: zsumx
      REAL(wp) :: zsumx2
      REAL(wp) :: zomb
      

      IF (lwp) THEN
         WRITE(numout,*) ''
         WRITE(numout,*) 'obs_wri_stats :'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
      ENDIF

      DO jvar = 1, fbdata%nvar
         zsumx=0.0_wp
         zsumx2=0.0_wp
         inumgoodobs=0
         DO jo = 1, fbdata%nobs
            DO jk = 1, fbdata%nlev
               IF ( ( fbdata%pob(jk,jo,jvar) < 9999.0 ) .AND. &
                  & ( fbdata%pdep(jk,jo) < 9999.0 ) .AND. &
                  & ( fbdata%padd(jk,jo,1,jvar) < 9999.0 ) ) THEN

                  zomb=fbdata%pob(jk, jo, jvar)-fbdata%padd(jk, jo, 1, jvar)
                  zsumx=zsumx+zomb
                  zsumx2=zsumx2+zomb**2
                  inumgoodobs=inumgoodobs+1
               ENDIF
            ENDDO
         ENDDO

         CALL obs_mpp_sum_integer( inumgoodobs, inumgoodobsmpp )
         CALL mpp_sum('obs_write', zsumx)
         CALL mpp_sum('obs_write', zsumx2)

         IF (lwp) THEN
            WRITE(numout,*) 'Type: ',fbdata%cname(jvar),'  Total number of good observations: ',inumgoodobsmpp 
            WRITE(numout,*) 'Overall mean obs minus model of the good observations: ',zsumx/inumgoodobsmpp
            WRITE(numout,*) 'Overall RMS obs minus model of the good observations: ',sqrt( zsumx2/inumgoodobsmpp )
            WRITE(numout,*) ''
         ENDIF

      ENDDO

   END SUBROUTINE obs_wri_stats

END MODULE obs_write
