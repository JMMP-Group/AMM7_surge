MODULE storng
!$AGRIF_DO_NOT_TREAT
   !!======================================================================
   !!                       ***  MODULE  storng  ***
   !! Random number generator, used in NEMO stochastic parameterization
   !!
   !!=====================================================================
   !! History :  3.3  ! 2011-10 (J.-M. Brankart)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !! The module is based on (and includes) the
   !! 64-bit KISS (Keep It Simple Stupid) random number generator
   !! distributed by George Marsaglia :
   !! http://groups.google.com/group/comp.lang.fortran/
   !!        browse_thread/thread/a85bf5f2a97f5a55
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   kiss          : 64-bit KISS random number generator (period ~ 2^250)
   !!   kiss_seed     : Define seeds for KISS random number generator
   !!   kiss_state    : Get current state of KISS random number generator
   !!   kiss_save     : Save current state of KISS (for future restart)
   !!   kiss_load     : Load the saved state of KISS
   !!   kiss_reset    : Reset the default seeds
   !!   kiss_check    : Check the KISS pseudo-random sequence
   !!   kiss_uniform  : Real random numbers with uniform distribution in [0,1]
   !!   kiss_gaussian : Real random numbers with Gaussian distribution N(0,1)
   !!   kiss_gamma    : Real random numbers with Gamma distribution Gamma(k,1)
   !!   kiss_sample   : Select a random sample from a set of integers
   !!
   !!   ---CURRENTLY NOT USED IN NEMO :
   !!   kiss_save, kiss_load, kiss_check, kiss_gamma, kiss_sample
   !!----------------------------------------------------------------------
   USE par_kind
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   ! Public functions/subroutines
   PUBLIC :: kiss, kiss_seed, kiss_state, kiss_reset ! kiss_save, kiss_load, kiss_check
   PUBLIC :: kiss_uniform, kiss_gaussian, kiss_gamma, kiss_sample

   ! Default/initial seeds
   INTEGER(KIND=i8) :: x=1234567890987654321_8
   INTEGER(KIND=i8) :: y=362436362436362436_8
   INTEGER(KIND=i8) :: z=1066149217761810_8
   INTEGER(KIND=i8) :: w=123456123456123456_8

   ! Parameters to generate real random variates
   REAL(KIND=wp), PARAMETER :: huge64=9223372036854775808.0  ! +1
   REAL(KIND=wp), PARAMETER :: zero=0.0, half=0.5, one=1.0, two=2.0

   ! Variables to store 2 Gaussian random numbers with current index (ig)
   INTEGER(KIND=i8), SAVE :: ig=1
   REAL(KIND=wp), SAVE :: gran1, gran2

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION kiss()
      !! --------------------------------------------------------------------
      !!                  ***  FUNCTION kiss  ***
      !!
      !! ** Purpose :   64-bit KISS random number generator
      !!
      !! ** Method  :   combine several random number generators:
      !!                (1) Xorshift (XSH), period 2^64-1,
      !!                (2) Multiply-with-carry (MWC), period (2^121+2^63-1)
      !!                (3) Congruential generator (CNG), period 2^64.
      !!
      !!                overall period:
      !!                (2^250+2^192+2^64-2^186-2^129)/6
      !!                            ~= 2^(247.42) or 10^(74.48)
      !!
      !!                set your own seeds with 'kiss_seed'
      ! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i8) :: kiss, t

      t = ISHFT(x,58) + w
      IF (s(x).eq.s(t)) THEN
         w = ISHFT(x,-6) + s(x)
      ELSE
         w = ISHFT(x,-6) + 1 - s(x+t)
      ENDIF
      x = t + x
      y = m( m( m(y,13_8), -17_8 ), 43_8 )
      z = 6906969069_8 * z + 1234567_8

      kiss = x + y + z

      CONTAINS

         FUNCTION s(k)
            INTEGER(KIND=i8) :: s, k
            s = ISHFT(k,-63)
         END FUNCTION s

         FUNCTION m(k, n)
            INTEGER(KIND=i8) :: m, k, n
            m =  IEOR(k, ISHFT(k, n) )
         END FUNCTION m

   END FUNCTION kiss


   SUBROUTINE kiss_seed(ix, iy, iz, iw)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_seed  ***
      !!
      !! ** Purpose :   Define seeds for KISS random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i8) :: ix, iy, iz, iw

      x = ix
      y = iy
      z = iz
      w = iw

   END SUBROUTINE kiss_seed


   SUBROUTINE kiss_state(ix, iy, iz, iw)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_state  ***
      !!
      !! ** Purpose :   Get current state of KISS random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i8) :: ix, iy, iz, iw

      ix = x
      iy = y
      iz = z
      iw = w

   END SUBROUTINE kiss_state


   SUBROUTINE kiss_reset()
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_reset  ***
      !!
      !! ** Purpose :   Reset the default seeds for KISS random number generator
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE

      x=1234567890987654321_8
      y=362436362436362436_8
      z=1066149217761810_8
      w=123456123456123456_8

    END SUBROUTINE kiss_reset


   !  SUBROUTINE kiss_check(check_type)
   !    !! --------------------------------------------------------------------
   !    !!                  ***  ROUTINE kiss_check  ***
   !    !!
   !    !! ** Purpose :   Check the KISS pseudo-random sequence
   !    !!
   !    !! ** Method  :   Check that it reproduces the correct sequence
   !    !!                from the default seed
   !    !!
   !    !! --------------------------------------------------------------------
   !    IMPLICIT NONE
   !    INTEGER(KIND=i8) :: iter, niter, correct, iran
   !    CHARACTER(LEN=*) :: check_type
   !    LOGICAL :: print_success

   !    ! Save current state of KISS
   !    CALL kiss_save()
   !    ! Reset the default seed
   !    CALL kiss_reset()

   !    ! Select check type
   !    SELECT CASE(check_type)
   !    CASE('short')
   !       niter = 5_8
   !       correct = 542381058189297533
   !       print_success = .FALSE.
   !    CASE('long')
   !       niter = 100000000_8
   !       correct = 1666297717051644203 ! Check provided by G. Marsaglia
   !       print_success = .TRUE.
   !    CASE('default')
   !    CASE DEFAULT
   !       STOP 'Bad check type in kiss_check'
   !    END SELECT

   !    ! Run kiss for the required number of iterations (niter)
   !    DO iter=1,niter
   !       iran = kiss()
   !    ENDDO

   !    ! Check that last iterate is correct
   !    IF (iran.NE.correct) THEN
   !       STOP 'Check failed: KISS internal error !!'
   !    ELSE
   !       IF (print_success) PRINT *, 'Check successful: 100 million calls to KISS OK'
   !    ENDIF

   !    ! Reload the previous state of KISS
   !    CALL kiss_load()

   ! END SUBROUTINE kiss_check


   ! SUBROUTINE kiss_save
   !    !! --------------------------------------------------------------------
   !    !!                  ***  ROUTINE kiss_save  ***
   !    !!
   !    !! ** Purpose :   Save current state of KISS random number generator
   !    !!
   !    !! --------------------------------------------------------------------
   !    INTEGER :: inum     !! Local integer

   !    IMPLICIT NONE

   !    CALL ctl_opn( inum, '.kiss_restart', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )

   !    ! OPEN(UNIT=30,FILE='.kiss_restart')
   !    WRITE(inum,*) x
   !    WRITE(inum,*) y
   !    WRITE(inum,*) z
   !    WRITE(inum,*) w
   !    CALL flush(inum)

   !  END SUBROUTINE kiss_save


   !  SUBROUTINE kiss_load
   !    !! --------------------------------------------------------------------
   !    !!                  ***  ROUTINE kiss_load  ***
   !    !!
   !    !! ** Purpose :   Load the saved state of KISS random number generator
   !    !!
   !    !! --------------------------------------------------------------------
   !    IMPLICIT NONE
   !    LOGICAL :: filexists
   ! Use ctl_opn routine rather than fortran intrinsic functions
   !    INQUIRE(FILE='.kiss_restart',EXIST=filexists)
   !    IF (filexists) THEN
   !       OPEN(UNIT=30,FILE='.kiss_restart')
   !       READ(30,*) x
   !       READ(30,*) y
   !       READ(30,*) z
   !       READ(30,*) w
   !       CLOSE(30)
   !    ENDIF

   ! END SUBROUTINE kiss_load


   SUBROUTINE kiss_uniform(uran)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_uniform  ***
      !!
      !! ** Purpose :   Real random numbers with uniform distribution in [0,1]
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=wp) :: uran

      uran = half * ( one + REAL(kiss(),wp) / huge64 )

   END SUBROUTINE kiss_uniform


   SUBROUTINE kiss_gaussian(gran)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_gaussian  ***
      !!
      !! ** Purpose :   Real random numbers with Gaussian distribution N(0,1)
      !!
      !! ** Method  :   Generate 2 new Gaussian draws (gran1 and gran2)
      !!                from 2 uniform draws on [-1,1] (u1 and u2),
      !!                using the Marsaglia polar method
      !!                (see Devroye, Non-Uniform Random Variate Generation, p. 235-236)
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=wp) :: gran, u1, u2, rsq, fac

      IF (ig.EQ.1) THEN
         rsq = two
         DO WHILE ( (rsq.GE.one).OR. (rsq.EQ.zero) )
            u1 = REAL(kiss(),wp) / huge64
            u2 = REAL(kiss(),wp) / huge64
            rsq = u1*u1 + u2*u2
         ENDDO
         fac = SQRT(-two*LOG(rsq)/rsq)
         gran1 = u1 * fac
         gran2 = u2 * fac
      ENDIF

      ! Output one of the 2 draws
      IF (ig.EQ.1) THEN
         gran = gran1 ; ig = 2
      ELSE
         gran = gran2 ; ig = 1
      ENDIF

   END SUBROUTINE kiss_gaussian


   SUBROUTINE kiss_gamma(gamr,k)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_gamma  ***
      !!
      !! ** Purpose :   Real random numbers with Gamma distribution Gamma(k,1)
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=wp), PARAMETER :: p1 = 4.5_8
      REAL(KIND=wp), PARAMETER :: p2 = 2.50407739677627_8  ! 1+LOG(9/2)
      REAL(KIND=wp), PARAMETER :: p3 = 1.38629436111989_8  ! LOG(4)
      REAL(KIND=wp) :: gamr, k, u1, u2, b, c, d, xx, yy, zz, rr, ee
      LOGICAL :: accepted

      IF (k.GT.one) THEN
         ! Cheng's rejection algorithm
         ! (see Devroye, Non-Uniform Random Variate Generation, p. 413)
         b = k - p3 ; d = SQRT(two*k-one) ; c = k + d

         accepted=.FALSE.
         DO WHILE (.NOT.accepted)
            CALL kiss_uniform(u1)
            yy = LOG(u1/(one-u1)) / d  ! Mistake in Devroye: "* k" instead of "/ d"
            xx = k * EXP(yy)
            rr = b + c * yy - xx
            CALL kiss_uniform(u2)
            zz = u1 * u1 * u2

            accepted = rr .GE. (zz*p1-p2)
            IF (.NOT.accepted) accepted =  rr .GE. LOG(zz)
         ENDDO

         gamr = xx

      ELSEIF (k.LT.one) THEN
        ! Rejection from the Weibull density
        ! (see Devroye, Non-Uniform Random Variate Generation, p. 415)
        c = one/k ; d = (one-k) * EXP( (k/(one-k)) * LOG(k) )

        accepted=.FALSE.
        DO WHILE (.NOT.accepted)
           CALL kiss_uniform(u1)
           zz = -LOG(u1)
           xx = EXP( c * LOG(zz) )
           CALL kiss_uniform(u2)
           ee = -LOG(u2)

           accepted = (zz+ee) .GE. (d+xx)  ! Mistake in Devroye: "LE" instead of "GE"
        ENDDO

        gamr = xx

      ELSE
         ! Exponential distribution
         CALL kiss_uniform(u1)
         gamr = -LOG(u1)

      ENDIF

   END SUBROUTINE kiss_gamma


   SUBROUTINE kiss_sample(a,n,k)
      !! --------------------------------------------------------------------
      !!                  ***  ROUTINE kiss_sample  ***
      !!
      !! ** Purpose :   Select a random sample of size k from a set of n integers
      !!
      !! ** Method  :   The sample is output in the first k elements of a
      !!                Set k equal to n to obtain a random permutation
      !!                  of the whole set of integers
      !!
      !! --------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER(KIND=i8), DIMENSION(:) :: a
      INTEGER(KIND=i8) :: n, k, i, j, atmp
      REAL(KIND=wp) :: uran

      ! Select the sample using the swapping method
      ! (see Devroye, Non-Uniform Random Variate Generation, p. 612)
      DO i=1,k
         ! Randomly select the swapping element between i and n (inclusive)
         CALL kiss_uniform(uran)
         j = i - 1 + CEILING( REAL(n-i+1,8) * uran )
         ! Swap elements i and j
         atmp = a(i) ; a(i) = a(j) ; a(j) = atmp
      ENDDO

   END SUBROUTINE kiss_sample
!$AGRIF_END_DO_NOT_TREAT
END MODULE storng
