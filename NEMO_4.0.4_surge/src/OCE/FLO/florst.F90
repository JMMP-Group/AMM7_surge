MODULE florst
   !!======================================================================
   !!                       ***  MODULE  florst  ***
   !! Ocean floats : write floats restart files
   !!======================================================================
   !!  History :  OPA  !  1999-09  (Y. Drillet)    : Original code
   !!              -   !  2000-06  (J.-M. Molines) : Profiling floats for CLS 
   !!   NEMO      1.0  !  2002-10  (A. Bozec)  F90 : Free form and module
   !!             3.2  !  2010-08  (slaw, cbricaud): netcdf outputs and others 
   !!----------------------------------------------------------------------
   USE flo_oce         ! ocean drifting floats
   USE dom_oce         ! ocean space and time domain
   USE lib_mpp         ! distribued memory computing library
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC flo_rst         ! routine called by floats.F90
   PUBLIC flo_rst_alloc   ! routine called by floats.F90

   INTEGER, ALLOCATABLE, DIMENSION(:) :: iperproc   ! 1D workspace

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION flo_rst_alloc()
      !!-------------------------------------------------------------------
      !!                ***  FUNCTION flo_rst_alloc  ***
      !!-------------------------------------------------------------------
      ALLOCATE( iperproc(jpnij), STAT=flo_rst_alloc )
      !  
      CALL mpp_sum ( 'florst', flo_rst_alloc )
      IF( flo_rst_alloc /= 0 )   CALL ctl_stop( 'STOP', 'flo_rst_alloc: failed to allocate arrays.' )
   END FUNCTION flo_rst_alloc


   SUBROUTINE flo_rst( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE flo_rst ***
      !!             
      !! ** Purpose :  
      !!             
      !!      
      !! ** Method  :   The frequency of  ??? is nwritefl
      !!      
      !!----------------------------------------------------------------------
      INTEGER  :: kt                            ! time step
      !
      CHARACTER (len=80)       :: clname             ! restart filename
      INTEGER                  :: ic , jc , jpn ,jfl ! temporary integer
      INTEGER                  :: inum               ! temporary logical unit for restart file
      !!----------------------------------------------------------------------

      IF(  ( MOD(kt,nn_stockfl) == 0 ) .OR. ( kt == nitend )  )THEN      

         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'flo_rst : write in  restart_float file '
            WRITE(numout,*) '~~~~~~~    '
         ENDIF

         ! file is opened and closed every time it is used.

         clname = 'restart.float.'
         ic = 1
         DO jc = 1, 16
            IF( cexper(jc:jc) /= ' ' ) ic = jc
         END DO
         clname = clname(1:14)//cexper(1:ic)
         ic = 1
         DO jc = 1, 48
            IF( clname(jc:jc) /= ' ' ) ic = jc
         END DO

         inum=0
         IF( lwp )THEN
            CALL ctl_opn( inum, clname, 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
            REWIND inum
         ENDIF
         !
         DO jpn = 1, jpnij
            iperproc(jpn) = 0
         END DO
         !
         IF(lwp) THEN
            REWIND(inum)
            WRITE (inum,*) tpifl,tpjfl,tpkfl,nisobfl,ngrpfl
            CLOSE (inum)
         ENDIF
         !
         ! Compute the number of trajectories for each processor
         !
         IF( lk_mpp ) THEN
            DO jfl = 1, jpnfl
               IF( (INT(tpifl(jfl)) >= mig(nldi)) .AND.   &
                  &(INT(tpifl(jfl)) <= mig(nlei)) .AND.   &
                  &(INT(tpjfl(jfl)) >= mjg(nldj)) .AND.   &
                  &(INT(tpjfl(jfl)) <= mjg(nlej)) ) THEN
                  iperproc(narea) = iperproc(narea)+1
               ENDIF
            END DO
            CALL mpp_sum( 'florst', iperproc, jpnij )
            !
            IF(lwp) THEN
               WRITE(numout,*) 'DATE',adatrj
               DO jpn = 1, jpnij
                  IF( iperproc(jpn) /= 0 ) THEN
                     WRITE(numout,*)'PROCESSOR',jpn-1,'compute',iperproc(jpn), 'trajectories.'
                  ENDIF
               END DO
            ENDIF
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE flo_rst

   !!=======================================================================
END MODULE florst
