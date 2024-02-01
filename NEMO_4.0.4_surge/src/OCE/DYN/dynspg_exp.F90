MODULE dynspg_exp
   !!======================================================================
   !!                   ***  MODULE  dynspg_exp  ***
   !! Ocean dynamics:  surface pressure gradient trend, explicit scheme
   !!======================================================================
   !! History :  2.0  !  2005-11  (V. Garnier, G. Madec, L. Bessieres) Original code
   !!            3.2  !  2009-06  (G. Madec, M. Leclair, R. Benshila) introduce sshwzv module
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_spg_exp  : update the momentum trend with the surface 
   !!                      pressure gradient in the free surface constant  
   !!                      volume case with vector optimization
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain 
   USE sbc_oce         ! surface boundary condition: ocean
   USE phycst          ! physical constants
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE iom             ! I/O library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_spg_exp   ! called in dynspg.F90 

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_spg_exp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_exp  ***
      !!
      !! ** Purpose :   Compute the now trend due to the surface pressure
      !!              gradient in case of explicit free surface formulation and 
      !!              add it to the general trend of momentum equation.
      !!
      !! ** Method  :   Explicit free surface formulation. Add to the general
      !!              momentum trend the surface pressure gradient :
      !!                      (ua,va) = (ua,va) + (spgu,spgv)
      !!              where spgu = -1/rau0 d/dx(ps) = -g/e1u di( sshn )
      !!                    spgv = -1/rau0 d/dy(ps) = -g/e2v dj( sshn )
      !!
      !! ** Action :   (ua,va)   trend of horizontal velocity increased by 
      !!                         the surf. pressure gradient trend
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  ::   kt   ! ocean time-step index
      !!
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_exp : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   (explicit free surface)'
         !
         spgu(:,:) = 0._wp   ;   spgv(:,:) = 0._wp
         !
         IF( .NOT.ln_linssh .AND. lwp ) WRITE(numout,*) '      non linear free surface: spg is included in dynhpg'
      ENDIF

      IF( ln_linssh ) THEN          !* linear free surface : add the surface pressure gradient trend
         !
         DO jj = 2, jpjm1                    ! now surface pressure gradient
            DO ji = fs_2, fs_jpim1   ! vector opt.
               spgu(ji,jj) = - grav * ( sshn(ji+1,jj) - sshn(ji,jj) ) * r1_e1u(ji,jj)
               spgv(ji,jj) = - grav * ( sshn(ji,jj+1) - sshn(ji,jj) ) * r1_e2v(ji,jj)
            END DO 
         END DO
         !
         DO jk = 1, jpkm1                    ! Add it to the general trend
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ua(ji,jj,jk) = ua(ji,jj,jk) + spgu(ji,jj)
                  va(ji,jj,jk) = va(ji,jj,jk) + spgv(ji,jj)
               END DO
            END DO
         END DO
         !
      ENDIF
      !
   END SUBROUTINE dyn_spg_exp

   !!======================================================================
END MODULE dynspg_exp
