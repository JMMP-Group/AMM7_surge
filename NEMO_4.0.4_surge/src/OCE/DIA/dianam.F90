MODULE dianam
   !!======================================================================
   !!                       ***  MODULE  dianam  ***
   !! Ocean diagnostics:  Builds output file name
   !!=====================================================================
   !! History :  OPA  ! 1999-02  (E. Guilyardi)  Creation for 30 days/month
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            3.2  ! 2009-11  (S. Masson) complete rewriting, works for all calendars...
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_nam       : Builds output file name
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE ioipsl, ONLY :  ju2ymds    ! for calendar

   IMPLICIT NONE
   PRIVATE

   PUBLIC dia_nam

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dia_nam( cdfnam, kfreq, cdsuff, ldfsec )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_nam  ***
      !!                   
      !! ** Purpose :   Builds output file name
      !!
      !! ** Method  :   File name is a function of date and output frequency
      !!      cdfnam=<cexper>_<clave>_<idtbeg>_<idtend>_<cdsuff>
      !!      <clave> = averaging frequency (DA, MO, etc...)
      !!      <idtbeg>,<idtend> date of beginning and end of run
      !!
      !!----------------------------------------------------------------------
      CHARACTER (len=*), INTENT(  out)           ::   cdfnam   ! file name
      CHARACTER (len=*), INTENT(in   )           ::   cdsuff   ! to be added at the end of the file name
      INTEGER          , INTENT(in   )           ::   kfreq    ! output frequency: > 0 in time-step (or seconds see ldfsec)
      !                                                                            < 0 in months
      !                                                                            = 0 no frequency
      LOGICAL          , INTENT(in   ), OPTIONAL ::   ldfsec   ! kfreq in second(in time-step) if .true.(.false. default)
      !
      CHARACTER (len=20) ::   clfmt, clfmt0                    ! writing format
      CHARACTER (len=20) ::   clave                            ! name for output frequency
      CHARACTER (len=20) ::   cldate1                          ! date of the beginning of run
      CHARACTER (len=20) ::   cldate2                          ! date of the end       of run
      LOGICAL            ::   llfsec                           ! local value of ldfsec
      INTEGER            ::   iyear1, imonth1, iday1           ! year, month, day of the first day of the run
      INTEGER            ::   iyear2, imonth2, iday2           ! year, month, day of the last  day of the run
      INTEGER            ::   indg                             ! number of digits needed to write a number     
      INTEGER            ::   inbsec, inbmn, inbhr             ! output frequency in seconds, minutes and hours
      INTEGER            ::   inbday, inbmo, inbyr             ! output frequency in days, months and years
      INTEGER            ::   iyyss, iddss, ihhss, immss       ! number of seconds in 1 year, 1 day, 1 hour and 1 minute
      INTEGER            ::   iyymo                            ! number of months in 1 year
      REAL(wp)           ::   zsec1, zsec2                     ! not used
      REAL(wp)           ::   zdrun, zjul                      ! temporary scalars
      !!----------------------------------------------------------------------

      ! name for output frequency

      IF( PRESENT(ldfsec) ) THEN   ;   llfsec = ldfsec
      ELSE                         ;   llfsec = .FALSE.
      ENDIF

      IF( llfsec .OR. kfreq < 0 ) THEN   ;   inbsec = kfreq                       ! output frequency already in seconds
      ELSE                               ;   inbsec = kfreq * NINT( rdt )   ! from time-step to seconds
      ENDIF
      iddss = NINT( rday          )                                         ! number of seconds in 1 day
      ihhss = NINT( rmmss * rhhmm )                                         ! number of seconds in 1 hour
      immss = NINT( rmmss         )                                         ! number of seconds in 1 minute
      iyymo = NINT( raamo         )                                         ! number of months  in 1 year
      iyyss = iddss * nyear_len(1)                                          ! seconds in 1 year (not good: multi years with leap)
      clfmt0 = "('(a,i',i1,',a)')"                                          ! format '(a,ix,a)' with x to be defined
      ! 
      IF(          inbsec == 0           ) THEN   ;   clave = ''            ! no frequency
      ELSEIF(      inbsec <  0           ) THEN        
         inbmo = -inbsec                                                    ! frequency in month
         IF( MOD( inbmo, iyymo  ) == 0 ) THEN                               ! frequency in years
            inbyr  = inbmo / iyymo
            indg   = INT(LOG10(REAL(inbyr,wp))) + 1                         ! number of digits needed to write years   frequency
            WRITE(clfmt, clfmt0) indg             ;   WRITE(clave, clfmt) '_', inbyr , 'y'
         ELSE                                                               ! frequency in month
            indg   = INT(LOG10(REAL(inbmo,wp))) + 1                         ! number of digits needed to write months  frequency
            WRITE(clfmt, clfmt0) indg             ;   WRITE(clave, clfmt) '_', inbmo, 'm'
         ENDIF
      ELSEIF( MOD( inbsec, iyyss  ) == 0 ) THEN                             ! frequency in years
         inbyr  = inbsec / iyyss
         indg   = INT(LOG10(REAL(inbyr ,wp))) + 1                           ! number of digits needed to write years   frequency
         WRITE(clfmt, clfmt0) indg                ;   WRITE(clave, clfmt) '_', inbyr , 'y'
      ELSEIF( MOD( inbsec, iddss  ) == 0 ) THEN                             ! frequency in days
         inbday = inbsec / iddss
         indg   = INT(LOG10(REAL(inbday,wp))) + 1                           ! number of digits needed to write days    frequency
         WRITE(clfmt, clfmt0) indg                ;   WRITE(clave, clfmt) '_', inbday, 'd'
         IF( inbday == nmonth_len(nmonth) )           clave = '_1m'
      ELSEIF( MOD( inbsec, ihhss ) == 0 ) THEN                              ! frequency in hours
         inbhr  = inbsec / ihhss
         indg   = INT(LOG10(REAL(inbhr ,wp))) + 1                           ! number of digits needed to write hours   frequency
         WRITE(clfmt, clfmt0) indg                ;   WRITE(clave, clfmt) '_', inbhr , 'h'
      ELSEIF( MOD( inbsec, immss ) == 0 ) THEN                              ! frequency in minutes
         inbmn  = inbsec / immss
         indg   = INT(LOG10(REAL(inbmn ,wp))) + 1                           ! number of digits needed to write minutes frequency
         WRITE(clfmt, clfmt0) indg                ;   WRITE(clave, clfmt) '_', inbmn , 'mn'
      ELSE                                                                  ! frequency in seconds
         indg   = INT(LOG10(REAL(inbsec,wp))) + 1                           ! number of digits needed to write seconds frequency
         WRITE(clfmt, clfmt0) indg                ;   WRITE(clave, clfmt) '_', inbsec, 's'
      ENDIF

      ! date of the beginning and the end of the run

      zdrun = rdt / rday * REAL( nitend - nit000, wp )                ! length of the run in days
      zjul  = fjulday - rdt / rday
      CALL ju2ymds( zjul        , iyear1, imonth1, iday1, zsec1 )           ! year/month/day of the beginning of run
      CALL ju2ymds( zjul + zdrun, iyear2, imonth2, iday2, zsec2 )           ! year/month/day of the end       of run

      IF( iyear2 < 10000 ) THEN   ;   clfmt = "(i4.4,2i2.2)"                ! format used to write the date 
      ELSE                        ;   WRITE(clfmt, "('(i',i1,',2i2.2)')") INT(LOG10(REAL(iyear2,wp))) + 1
      ENDIF

      WRITE(cldate1, clfmt) iyear1, imonth1, iday1                          ! date of the beginning of run
      WRITE(cldate2, clfmt) iyear2, imonth2, iday2                          ! date of the end       of run
 
      cdfnam = TRIM(cexper)//TRIM(clave)//"_"//TRIM(cldate1)//"_"//TRIM(cldate2)//"_"//TRIM(cdsuff)
      IF( .NOT. Agrif_Root() ) cdfnam = TRIM(Agrif_CFixed())//'_'//TRIM(cdfnam)

   END SUBROUTINE dia_nam

   !!======================================================================
END MODULE dianam
