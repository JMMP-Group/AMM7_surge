MODULE dyncor_c1d
   !!======================================================================
   !!                     ***  MODULE  dyncor_c1d  ***
   !! Ocean Dynamics :   Coriolis term in 1D configuration
   !!=====================================================================
   !! History :  2.0  !  2004-09  (C. Ethe)  Original code
   !!            3.0  !  2008-04  (G. Madec)  style only
   !!----------------------------------------------------------------------
#if defined key_c1d
   !!----------------------------------------------------------------------
   !!   'key_c1d'                                          1D Configuration
   !!----------------------------------------------------------------------
   !!   cor_c1d       : Coriolis factor at T-point (1D configuration)
   !!   dyn_cor_c1d   : vorticity trend due to Coriolis at T-point
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control

   USE sbcwave        ! Surface Waves (add Stokes-Coriolis force)
   USE sbc_oce , ONLY : ln_stcor    ! use Stoke-Coriolis force
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   cor_c1d      ! called by nemogcm.F90
   PUBLIC   dyn_cor_c1d  ! called by step1d.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE cor_c1d
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE cor_c1d  ***
      !! 
      !! ** Purpose : set the Coriolis factor at T-point
      !!----------------------------------------------------------------------
      REAL(wp) ::   zphi0, zbeta, zf0   ! local scalars
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'cor_c1d : Coriolis factor at T-point'
      IF(lwp) WRITE(numout,*) '~~~~~~~'

      !
   END SUBROUTINE cor_c1d


   SUBROUTINE dyn_cor_c1d( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_cor_c1d  ***
      !! 
      !! ** Purpose :   Compute the now Coriolis trend and add it to 
      !!               the general trend of the momentum equation in 1D case.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_cor_c1d : total vorticity trend in 1D'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      ENDIF
      !
      IF( ln_stcor ) THEN
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ua(ji,jj,jk) = ua(ji,jj,jk) + ff_t(ji,jj) * (vn(ji,jj,jk) + vsd(ji,jj,jk))
                  va(ji,jj,jk) = va(ji,jj,jk) - ff_t(ji,jj) * (un(ji,jj,jk) + usd(ji,jj,jk))
               END DO
            END DO
         END DO
      ELSE
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ua(ji,jj,jk) = ua(ji,jj,jk) + ff_t(ji,jj) * vn(ji,jj,jk)
                  va(ji,jj,jk) = va(ji,jj,jk) - ff_t(ji,jj) * un(ji,jj,jk)
               END DO
            END DO
         END DO
      END IF
      
      !
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' cor  - Ua: ', mask1=umask,  &
         &                       tab3d_2=va, clinfo2=' Va: '       , mask2=vmask )
      !
   END SUBROUTINE dyn_cor_c1d

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO 1D Configuration
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE cor_c1d              ! Empty routine
       IMPLICIT NONE
   END SUBROUTINE cor_c1d   
   SUBROUTINE dyn_cor_c1d ( kt )      ! Empty routine
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: kt
      WRITE(*,*) 'dyn_cor_c1d: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_cor_c1d
#endif

   !!=====================================================================
END MODULE dyncor_c1d
