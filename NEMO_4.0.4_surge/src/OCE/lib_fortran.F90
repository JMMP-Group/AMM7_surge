MODULE lib_fortran
   !!======================================================================
   !!                       ***  MODULE  lib_fortran  ***
   !! Fortran utilities:  includes some low levels fortran functionality
   !!======================================================================
   !! History :  3.2  !  2010-05  (M. Dunphy, R. Benshila)  Original code
   !!            3.4  !  2013-06  (C. Rousset)  add glob_min, glob_max 
   !!                                           + 3d dim. of input is fexible (jpk, jpl...) 
   !!            4.0  !  2016-06  (T. Lovato)  double precision global sum by default 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   glob_sum    : generic interface for global masked summation over
   !!                 the interior domain for 1 or 2 2D or 3D arrays
   !!                 it works only for T points
   !!   SIGN        : generic interface for SIGN to overwrite f95 behaviour
   !!                 of intrinsinc sign function
   !!----------------------------------------------------------------------
   USE par_oce         ! Ocean parameter
   USE dom_oce         ! ocean domain
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing
   USE lbclnk          ! ocean lateral boundary conditions

   IMPLICIT NONE
   PRIVATE

   PUBLIC   glob_sum      ! used in many places (masked with tmask_i)
   PUBLIC   glob_sum_full ! used in many places (masked with tmask_h, ie only over the halos)
   PUBLIC   local_sum     ! used in trcrad, local operation before glob_sum_delay
   PUBLIC   sum3x3        ! used in trcrad, do a sum over 3x3 boxes
   PUBLIC   DDPDD         ! also used in closea module
   PUBLIC   glob_min, glob_max
#if defined key_nosignedzero
   PUBLIC SIGN
#endif

   INTERFACE glob_sum
      MODULE PROCEDURE glob_sum_1d, glob_sum_2d, glob_sum_3d
   END INTERFACE
   INTERFACE glob_sum_full
      MODULE PROCEDURE glob_sum_full_2d, glob_sum_full_3d
   END INTERFACE
   INTERFACE local_sum
      MODULE PROCEDURE local_sum_2d, local_sum_3d
   END INTERFACE
   INTERFACE sum3x3
      MODULE PROCEDURE sum3x3_2d, sum3x3_3d
   END INTERFACE
   INTERFACE glob_min
      MODULE PROCEDURE glob_min_2d, glob_min_3d
   END INTERFACE
   INTERFACE glob_max
      MODULE PROCEDURE glob_max_2d, glob_max_3d
   END INTERFACE

#if defined key_nosignedzero
   INTERFACE SIGN
      MODULE PROCEDURE SIGN_SCALAR, SIGN_ARRAY_1D, SIGN_ARRAY_2D, SIGN_ARRAY_3D,   &
         &             SIGN_ARRAY_1D_A, SIGN_ARRAY_2D_A, SIGN_ARRAY_3D_A,          &
         &             SIGN_ARRAY_1D_B, SIGN_ARRAY_2D_B, SIGN_ARRAY_3D_B
   END INTERFACE
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

#  define GLOBSUM_CODE

#     define DIM_1d
#     define FUNCTION_GLOBSUM           glob_sum_1d
#     include "lib_fortran_generic.h90"
#     undef FUNCTION_GLOBSUM
#     undef DIM_1d

#     define DIM_2d
#     define OPERATION_GLOBSUM
#     define FUNCTION_GLOBSUM           glob_sum_2d
#     include "lib_fortran_generic.h90"
#     undef FUNCTION_GLOBSUM
#     undef OPERATION_GLOBSUM
#     define OPERATION_FULL_GLOBSUM
#     define FUNCTION_GLOBSUM           glob_sum_full_2d
#     include "lib_fortran_generic.h90"
#     undef FUNCTION_GLOBSUM
#     undef OPERATION_FULL_GLOBSUM
#     undef DIM_2d

#     define DIM_3d
#     define OPERATION_GLOBSUM
#     define FUNCTION_GLOBSUM           glob_sum_3d
#     include "lib_fortran_generic.h90"
#     undef FUNCTION_GLOBSUM
#     undef OPERATION_GLOBSUM
#     define OPERATION_FULL_GLOBSUM
#     define FUNCTION_GLOBSUM           glob_sum_full_3d
#     include "lib_fortran_generic.h90"
#     undef FUNCTION_GLOBSUM
#     undef OPERATION_FULL_GLOBSUM
#     undef DIM_3d

#  undef GLOBSUM_CODE


#  define GLOBMINMAX_CODE

#     define DIM_2d
#     define OPERATION_GLOBMIN
#     define FUNCTION_GLOBMINMAX           glob_min_2d
#     include "lib_fortran_generic.h90"
#     undef FUNCTION_GLOBMINMAX
#     undef OPERATION_GLOBMIN
#     define OPERATION_GLOBMAX
#     define FUNCTION_GLOBMINMAX           glob_max_2d
#     include "lib_fortran_generic.h90"
#     undef FUNCTION_GLOBMINMAX
#     undef OPERATION_GLOBMAX
#     undef DIM_2d

#     define DIM_3d
#     define OPERATION_GLOBMIN
#     define FUNCTION_GLOBMINMAX           glob_min_3d
#     include "lib_fortran_generic.h90"
#     undef FUNCTION_GLOBMINMAX
#     undef OPERATION_GLOBMIN
#     define OPERATION_GLOBMAX
#     define FUNCTION_GLOBMINMAX           glob_max_3d
#     include "lib_fortran_generic.h90"
#     undef FUNCTION_GLOBMINMAX
#     undef OPERATION_GLOBMAX
#     undef DIM_3d
#  undef GLOBMINMAX_CODE

!                          ! FUNCTION local_sum !

   FUNCTION local_sum_2d( ptab )
      !!----------------------------------------------------------------------
      REAL(wp),  INTENT(in   ) ::   ptab(:,:) ! array on which operation is applied
      COMPLEX(wp)              ::  local_sum_2d
      !
      !!-----------------------------------------------------------------------
      !
      COMPLEX(wp)::   ctmp
      REAL(wp)   ::   ztmp
      INTEGER    ::   ji, jj    ! dummy loop indices
      INTEGER    ::   ipi, ipj  ! dimensions
      !!-----------------------------------------------------------------------
      !
      ipi = SIZE(ptab,1)   ! 1st dimension
      ipj = SIZE(ptab,2)   ! 2nd dimension
      !
      ctmp = CMPLX( 0.e0, 0.e0, wp )   ! warning ctmp is cumulated

      DO jj = 1, ipj
         DO ji = 1, ipi
            ztmp =  ptab(ji,jj) * tmask_i(ji,jj)
            CALL DDPDD( CMPLX( ztmp, 0.e0, wp ), ctmp )
         END DO
      END DO
      !
      local_sum_2d = ctmp
       
   END FUNCTION local_sum_2d

   FUNCTION local_sum_3d( ptab )
      !!----------------------------------------------------------------------
      REAL(wp),  INTENT(in   ) ::   ptab(:,:,:) ! array on which operation is applied
      COMPLEX(wp)              ::  local_sum_3d
      !
      !!-----------------------------------------------------------------------
      !
      COMPLEX(wp)::   ctmp
      REAL(wp)   ::   ztmp
      INTEGER    ::   ji, jj, jk   ! dummy loop indices
      INTEGER    ::   ipi, ipj, ipk    ! dimensions
      !!-----------------------------------------------------------------------
      !
      ipi = SIZE(ptab,1)   ! 1st dimension
      ipj = SIZE(ptab,2)   ! 2nd dimension
      ipk = SIZE(ptab,3)   ! 3rd dimension
      !
      ctmp = CMPLX( 0.e0, 0.e0, wp )   ! warning ctmp is cumulated

      DO jk = 1, ipk
        DO jj = 1, ipj
          DO ji = 1, ipi
             ztmp =  ptab(ji,jj,jk) * tmask_i(ji,jj)
             CALL DDPDD( CMPLX( ztmp, 0.e0, wp ), ctmp )
          END DO
        END DO
      END DO
      !
      local_sum_3d = ctmp
       
   END FUNCTION local_sum_3d

!                          ! FUNCTION sum3x3 !

   SUBROUTINE sum3x3_2d( p2d )
      !!-----------------------------------------------------------------------
      !!                  ***  routine sum3x3_2d  ***
      !!
      !! ** Purpose : sum over 3x3 boxes
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION (:,:), INTENT(inout) ::   p2d
      !
      INTEGER ::   ji, ji2, jj, jj2     ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( SIZE(p2d,1) /= jpi ) CALL ctl_stop( 'STOP', 'wrong call of sum3x3_2d, the first dimension is not equal to jpi' ) 
      IF( SIZE(p2d,2) /= jpj ) CALL ctl_stop( 'STOP', 'wrong call of sum3x3_2d, the second dimension is not equal to jpj' ) 
      !
      DO jj = 1, jpj
         DO ji = 1, jpi 
            IF( MOD(mig(ji), 3) == 1 .AND. MOD(mjg(jj), 3) == 1 ) THEN   ! bottom left corber of a 3x3 box
               ji2 = MIN(mig(ji)+2, jpiglo) - nimpp + 1                  ! right position of the box
               jj2 = MIN(mjg(jj)+2, jpjglo) - njmpp + 1                  ! upper position of the box
               IF( ji2 <= jpi .AND. jj2 <= jpj ) THEN                    ! the box is fully included in the local mpi domain
                  p2d(ji:ji2,jj:jj2) = SUM(p2d(ji:ji2,jj:jj2))
               ENDIF
            ENDIF
         END DO
      END DO
      CALL lbc_lnk( 'lib_fortran', p2d, 'T', 1. )
      IF( nbondi /= -1 ) THEN
         IF( MOD(mig(    1), 3) == 1 )   p2d(    1,:) = p2d(    2,:)
         IF( MOD(mig(    1), 3) == 2 )   p2d(    2,:) = p2d(    1,:)
      ENDIF
      IF( nbondi /=  1 ) THEN
         IF( MOD(mig(jpi-2), 3) == 1 )   p2d(  jpi,:) = p2d(jpi-1,:)
         IF( MOD(mig(jpi-2), 3) == 0 )   p2d(jpi-1,:) = p2d(  jpi,:)
      ENDIF
      IF( nbondj /= -1 ) THEN
         IF( MOD(mjg(    1), 3) == 1 )   p2d(:,    1) = p2d(:,    2)
         IF( MOD(mjg(    1), 3) == 2 )   p2d(:,    2) = p2d(:,    1)
      ENDIF
      IF( nbondj /=  1 ) THEN
         IF( MOD(mjg(jpj-2), 3) == 1 )   p2d(:,  jpj) = p2d(:,jpj-1)
         IF( MOD(mjg(jpj-2), 3) == 0 )   p2d(:,jpj-1) = p2d(:,  jpj)
      ENDIF
      CALL lbc_lnk( 'lib_fortran', p2d, 'T', 1. )

   END SUBROUTINE sum3x3_2d

   SUBROUTINE sum3x3_3d( p3d )
      !!-----------------------------------------------------------------------
      !!                  ***  routine sum3x3_3d  ***
      !!
      !! ** Purpose : sum over 3x3 boxes
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION (:,:,:), INTENT(inout) ::   p3d
      !
      INTEGER ::   ji, ji2, jj, jj2, jn     ! dummy loop indices
      INTEGER ::   ipn                      ! Third dimension size
      !!----------------------------------------------------------------------
      !
      IF( SIZE(p3d,1) /= jpi ) CALL ctl_stop( 'STOP', 'wrong call of sum3x3_3d, the first dimension is not equal to jpi' ) 
      IF( SIZE(p3d,2) /= jpj ) CALL ctl_stop( 'STOP', 'wrong call of sum3x3_3d, the second dimension is not equal to jpj' ) 
      ipn = SIZE(p3d,3)
      !
      DO jn = 1, ipn
         DO jj = 1, jpj
            DO ji = 1, jpi 
               IF( MOD(mig(ji), 3) == 1 .AND. MOD(mjg(jj), 3) == 1 ) THEN   ! bottom left corber of a 3x3 box
                  ji2 = MIN(mig(ji)+2, jpiglo) - nimpp + 1                  ! right position of the box
                  jj2 = MIN(mjg(jj)+2, jpjglo) - njmpp + 1                  ! upper position of the box
                  IF( ji2 <= jpi .AND. jj2 <= jpj ) THEN                    ! the box is fully included in the local mpi domain
                     p3d(ji:ji2,jj:jj2,jn) = SUM(p3d(ji:ji2,jj:jj2,jn))
                  ENDIF
               ENDIF
            END DO
         END DO
      END DO
      CALL lbc_lnk( 'lib_fortran', p3d, 'T', 1. )
      IF( nbondi /= -1 ) THEN
         IF( MOD(mig(    1), 3) == 1 )   p3d(    1,:,:) = p3d(    2,:,:)
         IF( MOD(mig(    1), 3) == 2 )   p3d(    2,:,:) = p3d(    1,:,:)
      ENDIF
      IF( nbondi /=  1 ) THEN
         IF( MOD(mig(jpi-2), 3) == 1 )   p3d(  jpi,:,:) = p3d(jpi-1,:,:)
         IF( MOD(mig(jpi-2), 3) == 0 )   p3d(jpi-1,:,:) = p3d(  jpi,:,:)
      ENDIF
      IF( nbondj /= -1 ) THEN
         IF( MOD(mjg(    1), 3) == 1 )   p3d(:,    1,:) = p3d(:,    2,:)
         IF( MOD(mjg(    1), 3) == 2 )   p3d(:,    2,:) = p3d(:,    1,:)
      ENDIF
      IF( nbondj /=  1 ) THEN
         IF( MOD(mjg(jpj-2), 3) == 1 )   p3d(:,  jpj,:) = p3d(:,jpj-1,:)
         IF( MOD(mjg(jpj-2), 3) == 0 )   p3d(:,jpj-1,:) = p3d(:,  jpj,:)
      ENDIF
      CALL lbc_lnk( 'lib_fortran', p3d, 'T', 1. )

   END SUBROUTINE sum3x3_3d


   SUBROUTINE DDPDD( ydda, yddb )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE DDPDD ***
      !!
      !! ** Purpose : Add a scalar element to a sum
      !!
      !!
      !! ** Method  : The code uses the compensated summation with doublet
      !!              (sum,error) emulated useing complex numbers. ydda is the
      !!               scalar to add to the summ yddb
      !!
      !! ** Action  : This does only work for MPI.
      !!
      !! References : Using Acurate Arithmetics to Improve Numerical
      !!              Reproducibility and Sability in Parallel Applications
      !!              Yun HE and Chris H. Q. DING, Journal of Supercomputing 18, 259-277, 2001
      !!----------------------------------------------------------------------
      COMPLEX(wp), INTENT(in   ) ::   ydda
      COMPLEX(wp), INTENT(inout) ::   yddb
      !
      REAL(wp) :: zerr, zt1, zt2  ! local work variables
      !!-----------------------------------------------------------------------
      !
      ! Compute ydda + yddb using Knuth's trick.
      zt1  = REAL(ydda) + REAL(yddb)
      zerr = zt1 - REAL(ydda)
      zt2  = ( (REAL(yddb) - zerr) + (REAL(ydda) - (zt1 - zerr)) )   &
         &   + AIMAG(ydda)         + AIMAG(yddb)
      !
      ! The result is t1 + t2, after normalization.
      yddb = CMPLX( zt1 + zt2, zt2 - ((zt1 + zt2) - zt1), wp )
      !
   END SUBROUTINE DDPDD

#if defined key_nosignedzero
   !!----------------------------------------------------------------------
   !!   'key_nosignedzero'                                         F90 SIGN
   !!----------------------------------------------------------------------

   FUNCTION SIGN_SCALAR( pa, pb )
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_SCALAR  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa,pb          ! input
      REAL(wp) :: SIGN_SCALAR    ! result
      !!-----------------------------------------------------------------------
      IF ( pb >= 0.e0) THEN   ;   SIGN_SCALAR = ABS(pa)
      ELSE                    ;   SIGN_SCALAR =-ABS(pa)
      ENDIF
   END FUNCTION SIGN_SCALAR


   FUNCTION SIGN_ARRAY_1D( pa, pb )
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_1D  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa,pb(:)                   ! input
      REAL(wp) :: SIGN_ARRAY_1D(SIZE(pb,1))  ! result
      !!-----------------------------------------------------------------------
      WHERE ( pb >= 0.e0 )   ;   SIGN_ARRAY_1D = ABS(pa)
      ELSEWHERE              ;   SIGN_ARRAY_1D =-ABS(pa)
      END WHERE
   END FUNCTION SIGN_ARRAY_1D


   FUNCTION SIGN_ARRAY_2D(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_2D  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa,pb(:,:)      ! input
      REAL(wp) :: SIGN_ARRAY_2D(SIZE(pb,1),SIZE(pb,2))  ! result
      !!-----------------------------------------------------------------------
      WHERE ( pb >= 0.e0 )   ;   SIGN_ARRAY_2D = ABS(pa)
      ELSEWHERE              ;   SIGN_ARRAY_2D =-ABS(pa)
      END WHERE
   END FUNCTION SIGN_ARRAY_2D

   FUNCTION SIGN_ARRAY_3D(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_3D  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa,pb(:,:,:)      ! input
      REAL(wp) :: SIGN_ARRAY_3D(SIZE(pb,1),SIZE(pb,2),SIZE(pb,3))  ! result
      !!-----------------------------------------------------------------------
      WHERE ( pb >= 0.e0 )   ;   SIGN_ARRAY_3D = ABS(pa)
      ELSEWHERE              ;   SIGN_ARRAY_3D =-ABS(pa)
      END WHERE
   END FUNCTION SIGN_ARRAY_3D


   FUNCTION SIGN_ARRAY_1D_A(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_1D_A  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa(:),pb(:)      ! input
      REAL(wp) :: SIGN_ARRAY_1D_A(SIZE(pb,1))  ! result
      !!-----------------------------------------------------------------------
      WHERE ( pb >= 0.e0 )   ;   SIGN_ARRAY_1D_A = ABS(pa)
      ELSEWHERE              ;   SIGN_ARRAY_1D_A =-ABS(pa)
      END WHERE
   END FUNCTION SIGN_ARRAY_1D_A


   FUNCTION SIGN_ARRAY_2D_A(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_2D_A  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa(:,:),pb(:,:)      ! input
      REAL(wp) :: SIGN_ARRAY_2D_A(SIZE(pb,1),SIZE(pb,2))  ! result
      !!-----------------------------------------------------------------------
      WHERE ( pb >= 0.e0 )   ;   SIGN_ARRAY_2D_A = ABS(pa)
      ELSEWHERE              ;   SIGN_ARRAY_2D_A =-ABS(pa)
      END WHERE
   END FUNCTION SIGN_ARRAY_2D_A


   FUNCTION SIGN_ARRAY_3D_A(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_3D_A  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa(:,:,:),pb(:,:,:)  ! input
      REAL(wp) :: SIGN_ARRAY_3D_A(SIZE(pb,1),SIZE(pb,2),SIZE(pb,3)) ! result
      !!-----------------------------------------------------------------------
      WHERE ( pb >= 0.e0 )   ;   SIGN_ARRAY_3D_A = ABS(pa)
      ELSEWHERE              ;   SIGN_ARRAY_3D_A =-ABS(pa)
      END WHERE
   END FUNCTION SIGN_ARRAY_3D_A


   FUNCTION SIGN_ARRAY_1D_B(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_1D_B  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa(:),pb      ! input
      REAL(wp) :: SIGN_ARRAY_1D_B(SIZE(pa,1))  ! result
      !!-----------------------------------------------------------------------
      IF( pb >= 0.e0 ) THEN   ;   SIGN_ARRAY_1D_B = ABS(pa)
      ELSE                    ;   SIGN_ARRAY_1D_B =-ABS(pa)
      ENDIF
   END FUNCTION SIGN_ARRAY_1D_B


   FUNCTION SIGN_ARRAY_2D_B(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_2D_B  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa(:,:),pb      ! input
      REAL(wp) :: SIGN_ARRAY_2D_B(SIZE(pa,1),SIZE(pa,2))  ! result
      !!-----------------------------------------------------------------------
      IF( pb >= 0.e0 ) THEN   ;   SIGN_ARRAY_2D_B = ABS(pa)
      ELSE                    ;   SIGN_ARRAY_2D_B =-ABS(pa)
      ENDIF
   END FUNCTION SIGN_ARRAY_2D_B


   FUNCTION SIGN_ARRAY_3D_B(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_3D_B  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa(:,:,:),pb      ! input
      REAL(wp) :: SIGN_ARRAY_3D_B(SIZE(pa,1),SIZE(pa,2),SIZE(pa,3))  ! result
      !!-----------------------------------------------------------------------
      IF( pb >= 0.e0 ) THEN   ;   SIGN_ARRAY_3D_B = ABS(pa)
      ELSE                    ;   SIGN_ARRAY_3D_B =-ABS(pa)
      ENDIF
   END FUNCTION SIGN_ARRAY_3D_B
#endif

   !!======================================================================
END MODULE lib_fortran
