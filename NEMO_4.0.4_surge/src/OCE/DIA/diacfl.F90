MODULE diacfl
   !!======================================================================
   !!                       ***  MODULE  diacfl  ***
   !! Output CFL diagnostics to ascii file
   !!======================================================================
   !! History :  3.4  !  2010-03  (E. Blockley)  Original code
   !!            3.6  !  2014-06  (T. Graham)  Removed CPP key & Updated to vn3.6
   !!            4.0  !  2017-09  (G. Madec)  style + comments
   !!----------------------------------------------------------------------
   !!   dia_cfl        : Compute and output Courant numbers at each timestep
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE domvvl          ! 
   !
   USE lib_mpp         ! distribued memory computing
   USE lbclnk          ! ocean lateral boundary condition (or mpp link)
   USE in_out_manager  ! I/O manager
   USE iom             ! 
   USE timing          ! Performance output

   IMPLICIT NONE
   PRIVATE

   CHARACTER(LEN=50) :: clname="cfl_diagnostics.ascii"    ! ascii filename
   INTEGER           :: numcfl                            ! outfile unit
   !
   INTEGER, DIMENSION(3) ::   nCu_loc, nCv_loc, nCw_loc   ! U, V, and W run max locations in the global domain
   REAL(wp)              ::   rCu_max, rCv_max, rCw_max   ! associated run max Courant number 

   PUBLIC   dia_cfl       ! routine called by step.F90
   PUBLIC   dia_cfl_init  ! routine called by nemogcm

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_cfl ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_cfl  ***
      !!
      !! ** Purpose :  Compute the Courant numbers Cu=u*dt/dx and Cv=v*dt/dy
      !!               and output to ascii file 'cfl_diagnostics.ascii'
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER                          ::   ji, jj, jk                       ! dummy loop indices
      REAL(wp)                         ::   z2dt, zCu_max, zCv_max, zCw_max  ! local scalars
      INTEGER , DIMENSION(3)           ::   iloc_u , iloc_v , iloc_w , iloc  ! workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zCu_cfl, zCv_cfl, zCw_cfl        ! workspace
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dia_cfl')
      !
      !                       ! setup timestep multiplier to account for initial Eulerian timestep
      IF( neuler == 0 .AND. kt == nit000 ) THEN   ;    z2dt = rdt
      ELSE                                        ;    z2dt = rdt * 2._wp
      ENDIF
      !
      !                
      DO jk = 1, jpk       ! calculate Courant numbers
         DO jj = 1, jpj
            DO ji = 1, jpi
               zCu_cfl(ji,jj,jk) = ABS( un(ji,jj,jk) ) * z2dt / e1u  (ji,jj)      ! for i-direction
               zCv_cfl(ji,jj,jk) = ABS( vn(ji,jj,jk) ) * z2dt / e2v  (ji,jj)      ! for j-direction
               zCw_cfl(ji,jj,jk) = ABS( wn(ji,jj,jk) ) * z2dt / e3w_n(ji,jj,jk)   ! for k-direction
            END DO
         END DO         
      END DO
      !
      ! write outputs
      IF( iom_use('cfl_cu') )   CALL iom_put( 'cfl_cu', MAXVAL( zCu_cfl, dim=3 ) )
      IF( iom_use('cfl_cv') )   CALL iom_put( 'cfl_cv', MAXVAL( zCv_cfl, dim=3 ) )
      IF( iom_use('cfl_cw') )   CALL iom_put( 'cfl_cw', MAXVAL( zCw_cfl, dim=3 ) )

      !                    ! calculate maximum values and locations
      IF( lk_mpp ) THEN
         CALL mpp_maxloc( 'diacfl', zCu_cfl, umask, zCu_max, iloc_u )
         CALL mpp_maxloc( 'diacfl', zCv_cfl, vmask, zCv_max, iloc_v )
         CALL mpp_maxloc( 'diacfl', zCw_cfl, wmask, zCw_max, iloc_w )
      ELSE
         iloc = MAXLOC( ABS( zcu_cfl(:,:,:) ) )
         iloc_u(1) = iloc(1) + nimpp - 1
         iloc_u(2) = iloc(2) + njmpp - 1
         iloc_u(3) = iloc(3)
         zCu_max = zCu_cfl(iloc(1),iloc(2),iloc(3))
         !
         iloc = MAXLOC( ABS( zcv_cfl(:,:,:) ) )
         iloc_v(1) = iloc(1) + nimpp - 1
         iloc_v(2) = iloc(2) + njmpp - 1
         iloc_v(3) = iloc(3)
         zCv_max = zCv_cfl(iloc(1),iloc(2),iloc(3))
         !
         iloc = MAXLOC( ABS( zcw_cfl(:,:,:) ) )
         iloc_w(1) = iloc(1) + nimpp - 1
         iloc_w(2) = iloc(2) + njmpp - 1
         iloc_w(3) = iloc(3)
         zCw_max = zCw_cfl(iloc(1),iloc(2),iloc(3))
      ENDIF
      !
      !                    ! write out to file
      IF( lwp ) THEN
         WRITE(numcfl,FMT='(2x,i6,3x,a6,4x,f7.4,1x,i4,1x,i4,1x,i4)') kt, 'Max Cu', zCu_max, iloc_u(1), iloc_u(2), iloc_u(3)
         WRITE(numcfl,FMT='(11x,     a6,4x,f7.4,1x,i4,1x,i4,1x,i4)')     'Max Cv', zCv_max, iloc_v(1), iloc_v(2), iloc_v(3)
         WRITE(numcfl,FMT='(11x,     a6,4x,f7.4,1x,i4,1x,i4,1x,i4)')     'Max Cw', zCw_max, iloc_w(1), iloc_w(2), iloc_w(3)
      ENDIF
      !
      !                    ! update maximum Courant numbers from whole run if applicable
      IF( zCu_max > rCu_max ) THEN   ;   rCu_max = zCu_max   ;   nCu_loc(:) = iloc_u(:)   ;   ENDIF
      IF( zCv_max > rCv_max ) THEN   ;   rCv_max = zCv_max   ;   nCv_loc(:) = iloc_v(:)   ;   ENDIF
      IF( zCw_max > rCw_max ) THEN   ;   rCw_max = zCw_max   ;   nCw_loc(:) = iloc_w(:)   ;   ENDIF

      !                    ! at end of run output max Cu and Cv and close ascii file
      IF( kt == nitend .AND. lwp ) THEN
         ! to ascii file
         WRITE(numcfl,*) '******************************************'
         WRITE(numcfl,FMT='(3x,a12,6x,f7.4,1x,i4,1x,i4,1x,i4)') 'Run Max Cu', rCu_max, nCu_loc(1), nCu_loc(2), nCu_loc(3)
         WRITE(numcfl,FMT='(3x,a8,11x,f15.1)') ' => dt/C', z2dt/rCu_max
         WRITE(numcfl,*) '******************************************'
         WRITE(numcfl,FMT='(3x,a12,6x,f7.4,1x,i4,1x,i4,1x,i4)') 'Run Max Cv', rCv_max, nCv_loc(1), nCv_loc(2), nCv_loc(3)
         WRITE(numcfl,FMT='(3x,a8,11x,f15.1)') ' => dt/C', z2dt/rCv_max
         WRITE(numcfl,*) '******************************************'
         WRITE(numcfl,FMT='(3x,a12,6x,f7.4,1x,i4,1x,i4,1x,i4)') 'Run Max Cw', rCw_max, nCw_loc(1), nCw_loc(2), nCw_loc(3)
         WRITE(numcfl,FMT='(3x,a8,11x,f15.1)') ' => dt/C', z2dt/rCw_max
         CLOSE( numcfl ) 
         !
         ! to ocean output
         WRITE(numout,*)
         WRITE(numout,*) 'dia_cfl : Maximum Courant number information for the run '
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Max Cu = ', rCu_max, ' at (i,j,k) = (',nCu_loc(1),nCu_loc(2),nCu_loc(3),') => dt/C = ', z2dt/rCu_max
         WRITE(numout,*) '   Max Cv = ', rCv_max, ' at (i,j,k) = (',nCv_loc(1),nCv_loc(2),nCv_loc(3),') => dt/C = ', z2dt/rCv_max
         WRITE(numout,*) '   Max Cw = ', rCw_max, ' at (i,j,k) = (',nCw_loc(1),nCw_loc(2),nCw_loc(3),') => dt/C = ', z2dt/rCw_max
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('dia_cfl')
      !
   END SUBROUTINE dia_cfl


   SUBROUTINE dia_cfl_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_cfl_init  ***
      !!                   
      !! ** Purpose :   create output file, initialise arrays
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_cfl : Outputting CFL diagnostics to ',TRIM(clname), ' file'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*)
         !
         ! create output ascii file
         CALL ctl_opn( numcfl, clname, 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', 1, numout, lwp, 1 )
         WRITE(numcfl,*) 'Timestep  Direction  Max C     i    j    k'
         WRITE(numcfl,*) '******************************************'
      ENDIF
      !
      rCu_max = 0._wp
      rCv_max = 0._wp
      rCw_max = 0._wp
      !
   END SUBROUTINE dia_cfl_init

   !!======================================================================
END MODULE diacfl
