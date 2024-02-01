MODULE obs_utils
   !!======================================================================
   !!                       ***  MODULE obs_utils   ***
   !! Observation diagnostics: Utility functions
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   grt_cir_dis     : Great circle distance 
   !!   grt_cir_dis_saa : Great circle distance (small angle)
   !!   chkerr          : Error-message managment for NetCDF files
   !!   chkdim          : Error-message managment for NetCDF files
   !!   fatal_error     : Fatal error handling
   !!   ddatetoymdhms   : Convert YYYYMMDD.hhmmss to components
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce, ONLY : &        ! Precision variables
      & wp, &
      & dp, &
      & i8  
   USE in_out_manager           ! I/O manager
   USE lib_mpp                  ! For ctl_warn/stop

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE
   PUBLIC grt_cir_dis,     &  ! Great circle distance 
      &   grt_cir_dis_saa, &  ! Great circle distance (small angle)
      &   str_c_to_for,    &  ! Remove non-printable chars from string
      &   chkerr,          &  ! Error-message managment for NetCDF files
      &   chkdim,          &  ! Check if dimensions are correct for a variable
      &   fatal_error,     &  ! Fatal error handling
      &   warning,         &  ! Warning handling
      &   ddatetoymdhms       ! Convert YYYYMMDD.hhmmss to components
         
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS
 
#include "grt_cir_dis.h90"

#include "grt_cir_dis_saa.h90"

#include "str_c_to_for.h90"

   SUBROUTINE chkerr( kstatus, cd_name, klineno )
      !!----------------------------------------------------------------------
      !!   
      !!                    *** ROUTINE chkerr ***
      !! 
      !! ** Purpose : Error-message managment for NetCDF files.
      !! 
      !! ** Method  : 
      !!
      !! ** Action  :
      !!
      !! History
      !!      ! 02-12  (N. Daget)  hdlerr
      !!      ! 06-04  (A. Vidard) f90/nemovar migration, change name
      !!      ! 06-10  (A. Weaver) Cleanup
      !!----------------------------------------------------------------------
      !! * Modules used
      USE netcdf             ! NetCDF library
      USE dom_oce, ONLY : &  ! Ocean space and time domain variables
         & nproc

      !! * Arguments
      INTEGER :: kstatus
      INTEGER :: klineno
      CHARACTER(LEN=*) :: cd_name
      
      !! * Local declarations
      CHARACTER(len=200) :: clineno

      ! Main computation
      IF ( kstatus /= nf90_noerr ) THEN
         WRITE(clineno,'(A,I8)')' at line number ', klineno
         CALL ctl_stop( ' chkerr', ' Netcdf Error in ' // TRIM( cd_name ), &
            &           clineno, nf90_strerror( kstatus ) )
      ENDIF

   END SUBROUTINE chkerr

   SUBROUTINE chkdim( kfileid, kvarid, kndim, kdim, cd_name, klineno )
      !!----------------------------------------------------------------------
      !!   
      !!                    *** ROUTINE chkerr ***
      !! 
      !! ** Purpose : Error-message managment for NetCDF files.
      !! 
      !! ** Method  : 
      !!
      !! ** Action  :
      !!
      !! History
      !!      ! 07-03  (K. Mogenen + E. Remy) Initial version
      !!----------------------------------------------------------------------
      !! * Modules used
      USE netcdf             ! NetCDF library
      USE dom_oce, ONLY : &  ! Ocean space and time domain variables
         & nproc

      !! * Arguments
      INTEGER :: kfileid       ! NetCDF file id   
      INTEGER :: kvarid        ! NetCDF variable id   
      INTEGER :: kndim         ! Expected number of dimensions
      INTEGER, DIMENSION(kndim) :: kdim      ! Expected dimensions
      CHARACTER(LEN=*) :: cd_name            ! Calling routine name
      INTEGER ::  klineno      ! Calling line number

      !! * Local declarations
      INTEGER :: indim
      INTEGER, ALLOCATABLE, DIMENSION(:) :: &
         & idim,ilendim
      INTEGER :: ji
      LOGICAL :: llerr
      CHARACTER(len=200) :: clineno

      CALL chkerr( nf90_inquire_variable( kfileid, kvarid, ndims=indim ), &
         &         cd_name, klineno )

      ALLOCATE(idim(indim),ilendim(indim))

      CALL chkerr( nf90_inquire_variable( kfileid, kvarid, dimids=idim ), &
         &         cd_name, klineno )

      DO ji = 1, indim
         CALL chkerr( nf90_inquire_dimension( kfileid, idim(ji), &
            &                                 len=ilendim(ji) ), &
            &         cd_name, klineno )
      END DO
      
      IF ( indim /= kndim ) THEN
         WRITE(clineno,'(A,I8)')' at line number ', klineno
         CALL ctl_stop( ' chkdim',  &
            &           ' Netcdf no dim error in ' // TRIM( cd_name ), &
            &           clineno )
      ENDIF

      DO ji = 1, indim
         IF ( ilendim(ji) /= kdim(ji) ) THEN
            WRITE(clineno,'(A,I8)')' at line number ', klineno
            CALL ctl_stop( ' chkdim',  &
               &           ' Netcdf dim len error in ' // TRIM( cd_name ), &
               &           clineno )
         ENDIF
      END DO
         
      DEALLOCATE(idim,ilendim)

   END SUBROUTINE chkdim
   
  SUBROUTINE fatal_error( cd_name, klineno )
      !!----------------------------------------------------------------------
      !!
      !!                    *** ROUTINE fatal_error ***
      !!
      !! ** Purpose : Fatal error handling
      !!
      !! ** Method  :
      !!
      !! ** Action  :
      !!
      !! History
      !!----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER :: klineno
      CHARACTER(LEN=*) :: cd_name
      !! * Local declarations
      CHARACTER(len=200) :: clineno

      WRITE(clineno,'(A,I8)')' at line number ', klineno
      CALL ctl_stop( ' fatal_error', ' Error in ' // TRIM( cd_name ), &
         &           clineno)
      
   END SUBROUTINE fatal_error

   SUBROUTINE warning( cd_name, klineno )
      !!----------------------------------------------------------------------
      !!
      !!                    *** ROUTINE warning ***
      !!
      !! ** Purpose : Warning handling
      !!
      !! ** Method  :
      !!
      !! ** Action  :
      !!
      !! History
      !!----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER :: klineno
      CHARACTER(LEN=*) :: cd_name
      !! * Local declarations
      CHARACTER(len=200) :: clineno

      WRITE(clineno,'(A,I8)')' at line number ', klineno
      CALL ctl_warn( ' warning', ' Potential problem in ' // TRIM( cd_name ), &
         &           clineno)
      
   END SUBROUTINE warning

#include "ddatetoymdhms.h90"

END MODULE obs_utils
