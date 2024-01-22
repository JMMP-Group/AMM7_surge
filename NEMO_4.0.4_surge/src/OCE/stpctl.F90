MODULE stpctl
   !!======================================================================
   !!                       ***  MODULE  stpctl  ***
   !! Ocean run control :  gross check of the ocean time stepping
   !!======================================================================
   !! History :  OPA  ! 1991-03  (G. Madec) Original code
   !!            6.0  ! 1992-06  (M. Imbard)
   !!            8.0  ! 1997-06  (A.M. Treguier)
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            2.0  ! 2009-07  (G. Madec)  Add statistic for time-spliting
   !!            3.7  ! 2016-09  (G. Madec)  Remove solver
   !!            4.0  ! 2017-04  (G. Madec)  regroup global communications
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   stp_ctl      : Control the run
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE c1d             ! 1D vertical configuration
   USE diawri          ! Standard run outputs       (dia_wri_state routine)
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! distributed memory computing
   USE zdf_oce ,  ONLY : ln_zad_Aimp       ! ocean vertical physics variables
   USE wet_dry,   ONLY : ll_wd, ssh_ref    ! reference depth for negative bathy

   USE netcdf          ! NetCDF library
   IMPLICIT NONE
   PRIVATE

   PUBLIC stp_ctl           ! routine called by step.F90

   INTEGER  ::   idrun, idtime, idssh, idu, ids1, ids2, idt1, idt2, idc1, idw1, istatus
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE stp_ctl( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE stp_ctl  ***
      !!                     
      !! ** Purpose :   Control the run
      !!
      !! ** Method  : - Save the time step in numstp
      !!              - Print it each 50 time steps
      !!              - Stop the run IF problem encountered by setting nstop > 0
      !!                Problems checked: |ssh| maximum larger than 10 m
      !!                                  |U|   maximum larger than 10 m/s 
      !!                                  negative sea surface salinity
      !!
      !! ** Actions :   "time.step" file = last ocean time-step
      !!                "run.stat"  file = run statistics
      !!                 nstop indicator sheared among all local domain
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt       ! ocean time-step index
      !!
      INTEGER                ::   ji, jj, jk          ! dummy loop indices
      INTEGER,  DIMENSION(3) ::   ih, iu, is1, is2    ! min/max loc indices
      INTEGER,  DIMENSION(9) ::   iareasum, iareamin, iareamax
      REAL(wp)               ::   zzz                 ! local real 
      REAL(wp), DIMENSION(9) ::   zmax, zmaxlocal
      LOGICAL                ::   ll_wrtstp, ll_colruns, ll_wrtruns
      LOGICAL, DIMENSION(jpi,jpj,jpk) ::   llmsk
      CHARACTER(len=20) :: clname
      !!----------------------------------------------------------------------
      IF( nstop > 0 .AND. ngrdstop > -1 )   RETURN   !   stpctl was already called by a child grid
      !
      ll_wrtstp  = ( MOD( kt-nit000, sn_cfctl%ptimincr ) == 0 ) .OR. ( kt == nitend )
      ll_colruns = ll_wrtstp .AND. sn_cfctl%l_runstat .AND. jpnij > 1 
      ll_wrtruns = ( ll_colruns .OR. jpnij == 1 ) .AND. lwm
      !
      IF( kt == nit000 ) THEN
         !
         IF( lwp ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'stp_ctl : time-stepping control'
            WRITE(numout,*) '~~~~~~~'
         ENDIF
         !                                ! open time.step file
         IF( lwm ) CALL ctl_opn( numstp, 'time.step', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
         !                                ! open run.stat file(s) at start whatever
         !                                ! the value of sn_cfctl%ptimincr
         IF( ll_wrtruns ) THEN
            CALL ctl_opn( numrun, 'run.stat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
            clname = 'run.stat.nc'
            IF( .NOT. Agrif_Root() )   clname = TRIM(Agrif_CFixed())//"_"//TRIM(clname)
            istatus = NF90_CREATE( TRIM(clname), NF90_CLOBBER, idrun )
            istatus = NF90_DEF_DIM( idrun, 'time', NF90_UNLIMITED, idtime )
            istatus = NF90_DEF_VAR( idrun, 'abs_ssh_max', NF90_DOUBLE, (/ idtime /), idssh )
            istatus = NF90_DEF_VAR( idrun,   'abs_u_max', NF90_DOUBLE, (/ idtime /), idu   )
            istatus = NF90_DEF_VAR( idrun,       's_min', NF90_DOUBLE, (/ idtime /), ids1  )
            istatus = NF90_DEF_VAR( idrun,       's_max', NF90_DOUBLE, (/ idtime /), ids2  )
            istatus = NF90_DEF_VAR( idrun,       't_min', NF90_DOUBLE, (/ idtime /), idt1  )
            istatus = NF90_DEF_VAR( idrun,       't_max', NF90_DOUBLE, (/ idtime /), idt2  )
            IF( ln_zad_Aimp ) THEN
               istatus = NF90_DEF_VAR( idrun,   'abs_wi_max', NF90_DOUBLE, (/ idtime /), idw1  )
               istatus = NF90_DEF_VAR( idrun,       'Cf_max', NF90_DOUBLE, (/ idtime /), idc1  )
            ENDIF
            istatus = NF90_ENDDEF(idrun)
         ENDIF
      ENDIF
      !
      IF(lwm .AND. ll_wrtstp) THEN        !==  current time step  ==!   ("time.step" file)
         WRITE ( numstp, '(1x, i8)' )   kt
         REWIND( numstp )
      ENDIF
      !
      !                                   !==  test of extrema  ==!
      !
      ! define zmax default value. needed for land processors
      IF( ll_colruns ) THEN    ! default value: must not be kept when calling mpp_max -> must be as small as possible
         zmax(:) = -HUGE(1._wp)
      ELSE                     ! default value: must not give true for any of the tests bellow (-> avoid manipulating HUGE...)
         zmax(:) =  0._wp
         zmax(3) = -1._wp      ! avoid salinity minimum at 0.
      ENDIF
      !
      IF( ll_wd ) THEN
         zmax(1) = MAXVAL(  ABS( sshn(:,:) + ssh_ref*tmask(:,:,1) )  )        ! ssh max 
      ELSE
         zmax(1) = MAXVAL(  ABS( sshn(:,:) )  )                               ! ssh max
      ENDIF
      zmax(2) = MAXVAL(  ABS( un(:,:,:) )  )                                  ! velocity max (zonal only)
      llmsk(:,:,:) = tmask(:,:,:) == 1._wp
      IF( COUNT( llmsk(:,:,:) ) > 0 ) THEN   ! avoid huge values sent back for land processors...      
         zmax(4) = MAXVAL(  tsn(:,:,:,jp_sal)*tmask )            !       salinity max
         zmax(3) = MAXVAL( (zmax(4)-tsn(:,:,:,jp_sal))*tmask ) - zmax(4)          ! minus salinity max
         IF( ll_colruns .OR. jpnij == 1 ) THEN     ! following variables are used only in the netcdf file
            zmax(6) = MAXVAL(  tsn(:,:,:,jp_tem)*tmask )         !       temperature max
            zmax(5) = zmax(6) - MAXVAL( (zmax(6)-tsn(:,:,:,jp_tem))*tmask )        ! minus temperature max
            IF( ln_zad_Aimp ) THEN
               zmax(9) = MAXVAL(   Cu_adv(:,:,:)*tmask )         ! partitioning coeff. max
               llmsk(:,:,:) = wmask(:,:,:) == 1._wp
               IF( COUNT( llmsk(:,:,:) ) > 0 ) THEN   ! avoid huge values sent back for land processors...
                  zmax(8) = MAXVAL( ABS( wi(:,:,:) )*wmask(:,:,:) ) ! implicit vertical vel. max
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      zmax(7) = REAL( nstop , wp )                                            ! stop indicator
      !
      IF( ll_colruns ) THEN
         zmaxlocal(:) = zmax(:)
         CALL mpp_max( "stpctl", zmax )          ! max over the global domain
         nstop = NINT( zmax(7) )                 ! nstop indicator sheared among all local domains
      ENDIF
      !                                   !==  run statistics  ==!   ("run.stat" files)
      IF( ll_wrtruns ) THEN
         WRITE(numrun,9500) kt, zmax(1), zmax(2), -zmax(3), zmax(4)
         istatus = NF90_PUT_VAR( idrun, idssh, (/ zmax(1)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun,   idu, (/ zmax(2)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun,  ids1, (/-zmax(3)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun,  ids2, (/ zmax(4)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun,  idt1, (/-zmax(5)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun,  idt2, (/ zmax(6)/), (/kt/), (/1/) )
         IF( ln_zad_Aimp ) THEN
            istatus = NF90_PUT_VAR( idrun,  idw1, (/ zmax(8)/), (/kt/), (/1/) )
            istatus = NF90_PUT_VAR( idrun,  idc1, (/ zmax(9)/), (/kt/), (/1/) )
         ENDIF
         IF( kt == nitend ) istatus = NF90_CLOSE(idrun)
      END IF
      !                                   !==  error handling  ==!
      IF(   zmax(1) >   20._wp .OR.   &                    ! too large sea surface height ( > 20 m )
         &  zmax(2) >   10._wp .OR.   &                    ! too large velocity ( > 10 m/s)
         &  zmax(3) >=   0._wp .OR.   &                    ! negative or zero sea surface salinity
         &  zmax(4) >= 100._wp .OR.   &                    ! too large sea surface salinity ( > 100 )
         &  zmax(4) <    0._wp .OR.   &                    ! too large sea surface salinity (keep this line for sea-ice)
         &  ISNAN( zmax(1) + zmax(2) + zmax(3) ) .OR.  &   ! NaN encounter in the tests
         &  ABS(   zmax(1) + zmax(2) + zmax(3) ) > HUGE(1._wp) ) THEN    ! Infinity encounter in the tests
         IF( ll_colruns ) THEN
            ! first: close the netcdf file, so we can read it
            IF( lwm .AND. kt /= nitend )   istatus = NF90_CLOSE(idrun)
            CALL mpp_maxloc( 'stpctl', ABS(sshn)        , ssmask(:,:)  , zzz, ih(1:2)  )   ;   ih(3) = 0
            CALL mpp_maxloc( 'stpctl', ABS(un)          , umask (:,:,:), zzz, iu  )
            CALL mpp_minloc( 'stpctl', tsn(:,:,:,jp_sal), tmask (:,:,:), zzz, is1 )
            CALL mpp_maxloc( 'stpctl', tsn(:,:,:,jp_sal), tmask (:,:,:), zzz, is2 )
            ! find which subdomain has the max.
            iareamin(:) = jpnij+1   ;   iareamax(:) = 0   ;   iareasum(:) = 0
            DO ji = 1, 9
               IF( zmaxlocal(ji) == zmax(ji) ) THEN
                  iareamin(ji) = narea   ;   iareamax(ji) = narea   ;   iareasum(ji) = 1
               ENDIF
            END DO
            CALL mpp_min( "stpctl", iareamin )         ! min over the global domain
            CALL mpp_max( "stpctl", iareamax )         ! max over the global domain
            CALL mpp_sum( "stpctl", iareasum )         ! sum over the global domain
         ELSE
            ih(1:2)= MAXLOC( ABS( sshn(:,:)   )                              ) + (/ nimpp - 1, njmpp - 1    /)   ;   ih(3) = 0
            iu(:)  = MAXLOC( ABS( un  (:,:,:) )                              ) + (/ nimpp - 1, njmpp - 1, 0 /)
            is1(:) = MAXLOC( (zmax(4)-tsn(:,:,:,jp_sal))*tmask(:,:,:) ) + (/ nimpp - 1, njmpp - 1, 0 /)
            is2(:) = MAXLOC( tsn(:,:,:,jp_sal)*tmask(:,:,:) ) + (/ nimpp - 1, njmpp - 1, 0 /)
            iareamin(:) = narea   ;   iareamax(:) = narea   ;   iareasum(:) = 1         ! this is local information
         ENDIF
         !
         WRITE(ctmp1,*) ' stp_ctl: |ssh| > 20 m  or  |U| > 10 m/s  or  S <= 0  or  S >= 100  or  NaN encounter in the tests'
         CALL wrt_line(ctmp2, kt, ' |ssh| max ',   zmax(1), ih , iareasum(1), iareamin(1), iareamax(1) ) 
         CALL wrt_line(ctmp3, kt, ' |U|   max ',   zmax(2), iu , iareasum(2), iareamin(2), iareamax(2) ) 
         CALL wrt_line(ctmp4, kt, ' Sal   min ', - zmax(3), is1, iareasum(3), iareamin(3), iareamax(3) ) 
         CALL wrt_line(ctmp5, kt, ' Sal   max ',   zmax(4), is2, iareasum(4), iareamin(4), iareamax(4) ) 
         IF( Agrif_Root() ) THEN
            WRITE(ctmp6,*) '      ===> output of last computed fields in output.abort* files'
         ELSE
            WRITE(ctmp6,*) '      ===> output of last computed fields in '//TRIM(Agrif_CFixed())//'_output.abort* files'
         ENDIF
         !
         CALL dia_wri_state( 'output.abort' )    ! create an output.abort file
         !
         IF( ll_colruns .or. jpnij == 1 ) THEN   ! all processes synchronized -> use lwp to print in opened ocean.output files
            IF(lwp) THEN   ;   CALL ctl_stop( ctmp1, ' ', ctmp2, ctmp3, ctmp4, ctmp5, ' ', ctmp6 )
            ELSE           ;   nstop = MAX(1, nstop)   ! make sure nstop > 0 (automatically done when calling ctl_stop)
            ENDIF
         ELSE                                    ! only mpi subdomains with errors are here -> STOP now
            CALL ctl_stop( 'STOP', ctmp1, ' ', ctmp2, ctmp3, ctmp4, ctmp5, ' ', ctmp6 )
         ENDIF
         !
      ENDIF
      !
      IF( nstop > 0 ) THEN                                                  ! an error was detected and we did not abort yet...
         ngrdstop = Agrif_Fixed()                                           ! store which grid got this error
         IF( .NOT. ll_colruns .AND. jpnij > 1 )   CALL ctl_stop( 'STOP' )   ! we must abort here to avoid MPI deadlock
      ENDIF
      !
9500  FORMAT(' it :', i8, '    |ssh|_max: ', D23.16, ' |U|_max: ', D23.16,' S_min: ', D23.16,' S_max: ', D23.16)
      !
   END SUBROUTINE stp_ctl


   SUBROUTINE wrt_line( cdline, kt, cdprefix, pval, kloc, ksum, kmin, kmax )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE wrt_line  ***
      !!
      !! ** Purpose :   write information line
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=*),      INTENT(  out) ::   cdline
      CHARACTER(len=*),      INTENT(in   ) ::   cdprefix
      REAL(wp),              INTENT(in   ) ::   pval
      INTEGER, DIMENSION(3), INTENT(in   ) ::   kloc
      INTEGER,               INTENT(in   ) ::   kt, ksum, kmin, kmax
      !
      CHARACTER(len=80) ::   clsuff
      CHARACTER(len=9 ) ::   clkt, clsum, clmin, clmax
      CHARACTER(len=9 ) ::   cli, clj, clk
      CHARACTER(len=1 ) ::   clfmt
      CHARACTER(len=4 ) ::   cl4   ! needed to be able to compile with Agrif, I don't know why
      INTEGER           ::   ifmtk
      !!----------------------------------------------------------------------
      WRITE(clkt , '(i9)') kt
      
      WRITE(clfmt, '(i1)') INT(LOG10(REAL(jpnij  ,wp))) + 1     ! how many digits to we need to write ? (we decide max = 9)
      !!! WRITE(clsum, '(i'//clfmt//')') ksum                   ! this is creating a compilation error with AGRIF
      cl4 = '(i'//clfmt//')'   ;   WRITE(clsum, cl4) ksum
      WRITE(clfmt, '(i1)') INT(LOG10(REAL(MAX(1,jpnij-1),wp))) + 1    ! how many digits to we need to write ? (we decide max = 9)
      cl4 = '(i'//clfmt//')'   ;   WRITE(clmin, cl4) kmin-1
                                   WRITE(clmax, cl4) kmax-1
      !
      WRITE(clfmt, '(i1)') INT(LOG10(REAL(jpiglo,wp))) + 1      ! how many digits to we need to write jpiglo? (we decide max = 9)
      cl4 = '(i'//clfmt//')'   ;   WRITE(cli, cl4) kloc(1)      ! this is ok with AGRIF
      WRITE(clfmt, '(i1)') INT(LOG10(REAL(jpjglo,wp))) + 1      ! how many digits to we need to write jpjglo? (we decide max = 9)
      cl4 = '(i'//clfmt//')'   ;   WRITE(clj, cl4) kloc(2)      ! this is ok with AGRIF
      !
      IF( ksum == 1 ) THEN   ;   WRITE(clsuff,9100) TRIM(clmin)
      ELSE                   ;   WRITE(clsuff,9200) TRIM(clsum), TRIM(clmin), TRIM(clmax)
      ENDIF
      IF(kloc(3) == 0) THEN
         ifmtk = INT(LOG10(REAL(jpk,wp))) + 1                   ! how many digits to we need to write jpk? (we decide max = 9)
         clk = REPEAT(' ', ifmtk)                               ! create the equivalent in blank string
         WRITE(cdline,9300) TRIM(ADJUSTL(clkt)), TRIM(ADJUSTL(cdprefix)), pval, TRIM(cli), TRIM(clj), clk(1:ifmtk), TRIM(clsuff)
      ELSE
         WRITE(clfmt, '(i1)') INT(LOG10(REAL(jpk,wp))) + 1      ! how many digits to we need to write jpk? (we decide max = 9)
         !!! WRITE(clk, '(i'//clfmt//')') kloc(3)               ! this is creating a compilation error with AGRIF
         cl4 = '(i'//clfmt//')'   ;   WRITE(clk, cl4) kloc(3)   ! this is ok with AGRIF
         WRITE(cdline,9400) TRIM(ADJUSTL(clkt)), TRIM(ADJUSTL(cdprefix)), pval, TRIM(cli), TRIM(clj),    TRIM(clk), TRIM(clsuff)
      ENDIF
      !
9100  FORMAT('MPI rank ', a)
9200  FORMAT('found in ', a, ' MPI tasks, spread out among ranks ', a, ' to ', a)
9300  FORMAT('kt ', a, ' ', a, ' ', 1pg11.4, ' at i j   ', a, ' ', a, ' ', a, ' ', a)
9400  FORMAT('kt ', a, ' ', a, ' ', 1pg11.4, ' at i j k ', a, ' ', a, ' ', a, ' ', a)
      !
   END SUBROUTINE wrt_line


   !!======================================================================
END MODULE stpctl
