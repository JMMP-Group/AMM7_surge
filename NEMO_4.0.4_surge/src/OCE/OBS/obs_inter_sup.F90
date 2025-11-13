MODULE obs_inter_sup
   !!=====================================================================
   !!                       ***  MODULE obs_inter_sup  ***
   !! Observation diagnostics: Support for interpolation
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   obs_int_comm_3d : Get 3D interpolation stencil
   !!   obs_int_comm_2d : Get 2D interpolation stencil
   !!---------------------------------------------------------------------
   !! * Modules used
   USE par_kind        ! Precision variables
   USE dom_oce         ! Domain variables
   USE mpp_map         ! Map of processor points
   USE lib_mpp         ! MPP stuff
   USE obs_mpp         ! MPP stuff for observations
   USE obs_grid        ! Grid tools
   USE in_out_manager  ! I/O stuff

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE
   
   PUBLIC obs_int_comm_3d, & ! Get 3D interpolation stencil
      &   obs_int_comm_2d    ! Get 2D interpolation stencil
   
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE obs_int_comm_3d( kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, &
      &                        pval, pgval, kproc )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE obs_int_comm_3d  ***
      !!          
      !! ** Purpose : Get 3D interpolation stencil
      !!
      !! ** Method  : Either on-demand communication with 
      !!              obs_int_comm_3d_global 
      !!              or local memory with
      !!              obs_int_comm_3D_local
      !!              depending on ln_global_grid
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  08-02  (K. Mogensen)  Original code
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(IN) :: kptsi     ! Number of i horizontal points per stencil 
      INTEGER, INTENT(IN) :: kptsj     ! Number of j horizontal points per stencil
      INTEGER, INTENT(IN) :: kobs      ! Local number of observations
      INTEGER, INTENT(IN) :: kpi       ! Number of points in i direction
      INTEGER, INTENT(IN) :: kpj       ! Number of points in j direction
      INTEGER, INTENT(IN) :: kpk       ! Number of levels
      INTEGER, DIMENSION(kptsi,kptsj,kobs), INTENT(IN) :: &
         & kgrdi, &         ! i,j indicies for each stencil
         & kgrdj
      INTEGER, OPTIONAL, DIMENSION(kptsi,kptsj,kobs), INTENT(IN) :: &
         & kproc            ! Precomputed processor for each i,j,iobs points
      REAL(KIND=wp), DIMENSION(kpi,kpj,kpk), INTENT(IN) ::&
         & pval             ! Local 3D array to extract data from
      REAL(KIND=wp), DIMENSION(kptsi,kptsj,kpk,kobs), INTENT(OUT) ::&
         & pgval            ! Stencil at each point
      !! * Local declarations
      
      IF (ln_grid_global) THEN
         
         IF (PRESENT(kproc)) THEN

            CALL obs_int_comm_3d_global( kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, &
               &                         kgrdj, pval, pgval, kproc=kproc )

         ELSE

            CALL obs_int_comm_3d_global( kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, &
               &                         kgrdj, pval, pgval )

         ENDIF

      ELSE

         CALL obs_int_comm_3d_local( kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, &
            &                        pval, pgval )

      ENDIF

   END SUBROUTINE obs_int_comm_3d

   SUBROUTINE obs_int_comm_2d( kptsi, kptsj, kobs, kpi, kpj, kgrdi, kgrdj, pval, pgval, &
      &                        kproc )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE obs_int_comm_2d  ***
      !!          
      !! ** Purpose : Get 2D interpolation stencil
      !!
      !! ** Method  : Call to obs_int_comm_3d
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  08-02  (K. Mogensen)  Original code
      !!----------------------------------------------------------------------
      !!
      !! * Arguments
      INTEGER, INTENT(IN) :: kptsi        ! Number of i horizontal points per stencil
      INTEGER, INTENT(IN) :: kptsj        ! Number of j horizontal points per stencil
      INTEGER, INTENT(IN) :: kobs          ! Local number of observations
      INTEGER, INTENT(IN) :: kpi          ! Number of model grid points in i direction
      INTEGER, INTENT(IN) :: kpj          ! Number of model grid points in j direction
      INTEGER, DIMENSION(kptsi,kptsj,kobs), INTENT(IN) :: &
         & kgrdi, &         ! i,j indicies for each stencil
         & kgrdj
      INTEGER, OPTIONAL, DIMENSION(kptsi,kptsj,kobs), INTENT(IN) :: &
         & kproc            ! Precomputed processor for each i,j,iobs points
      REAL(KIND=wp), DIMENSION(kpi,kpj), INTENT(IN) ::&
         & pval             ! Local 3D array to extra data from
      REAL(KIND=wp), DIMENSION(kptsi,kptsj,kobs), INTENT(OUT) ::&
         & pgval            ! Stencil at each point
      !! * Local declarations
      REAL(KIND=wp), DIMENSION(jpi,jpj,1) ::   zval
      REAL(KIND=wp), DIMENSION(kptsi,kptsj,1,kobs) ::&
         & zgval 

      ! Set up local "3D" buffer

      zval(:,:,1) = pval(:,:)

      ! Call the 3D version

      IF (PRESENT(kproc)) THEN

         CALL obs_int_comm_3d( kptsi, kptsj, kobs, kpi, kpj, 1, kgrdi, kgrdj, zval, &
            &                  zgval, kproc=kproc )
      ELSE

         CALL obs_int_comm_3d( kptsi, kptsj, kobs, kpi, kpj, 1, kgrdi, kgrdj, zval, &
            &                  zgval )

      ENDIF

      ! Copy "3D" data back to 2D

      pgval(:,:,:) = zgval(:,:,1,:)

   END SUBROUTINE obs_int_comm_2d

   SUBROUTINE obs_int_comm_3d_global( kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, &
      &                               pval, pgval, kproc )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE obs_int_comm_3d_global  ***
      !!          
      !! ** Purpose : Get 3D interpolation stencil (global version)
      !!
      !! ** Method  : On-demand communication where each processor send its
      !!              list of (i,j) of points to all processors and receive 
      !!              the corresponding values
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  08-02  (K. Mogensen)  Original code
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(IN) :: kptsi     ! Number of i horizontal points per stencil 
      INTEGER, INTENT(IN) :: kptsj     ! Number of j horizontal points per stencil
      INTEGER, INTENT(IN) :: kobs      ! Local number of observations
      INTEGER, INTENT(IN) :: kpi       ! Number of model points in i direction
      INTEGER, INTENT(IN) :: kpj       ! Number of model points in j direction
      INTEGER, INTENT(IN) :: kpk       ! Number of levels
      INTEGER, DIMENSION(kptsi,kptsj,kobs), INTENT(IN) :: &
         & kgrdi, &         ! i,j indicies for each stencil
         & kgrdj
      INTEGER, OPTIONAL, DIMENSION(kptsi,kptsj,kobs), INTENT(IN) :: &
         & kproc            ! Precomputed processor for each i,j,iobs points
      REAL(KIND=wp), DIMENSION(kpi,kpj,kpk), INTENT(IN) ::&
         & pval             ! Local 3D array to extract data from
      REAL(KIND=wp), DIMENSION(kptsi,kptsj,kpk,kobs), INTENT(OUT) ::&
         & pgval            ! Stencil at each point
      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zsend, &
         & zrecv
      INTEGER, DIMENSION(:), ALLOCATABLE :: &
         & igrdij_send, &
         & igrdij_recv
      INTEGER, DIMENSION(kptsi,kptsj,kobs) :: &
         & iorder, &
         & iproc
      INTEGER :: nplocal(jpnij)
      INTEGER :: npglobal(jpnij)
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jk
      INTEGER :: jp
      INTEGER :: jobs
      INTEGER :: it
      INTEGER :: itot
      INTEGER :: ii
      INTEGER :: ij

      ! Check valid points

      IF ( ( MAXVAL(kgrdi) > jpiglo ) .OR. ( MINVAL(kgrdi) < 1 ) .OR. &
         & ( MAXVAL(kgrdj) > jpjglo ) .OR. ( MINVAL(kgrdj) < 1 ) ) THEN

         CALL ctl_stop( 'Error in obs_int_comm_3d_global', &
            &           'Point outside global domain' )

      ENDIF

      ! Count number of points on each processors

      nplocal(:) = 0
      IF (PRESENT(kproc)) THEN 
         iproc(:,:,:) = kproc(:,:,:)
         DO jobs = 1, kobs
            DO jj = 1, kptsj
               DO ji = 1, kptsi
                  nplocal(iproc(ji,jj,jobs)) = nplocal(iproc(ji,jj,jobs)) + 1
               END DO
            END DO
         END DO
      ELSE
         DO jobs = 1, kobs
            DO jj = 1, kptsj
               DO ji = 1, kptsi
                  iproc(ji,jj,jobs) = mppmap(kgrdi(ji,jj,jobs),&
                     &                       kgrdj(ji,jj,jobs))
                  nplocal(iproc(ji,jj,jobs)) = nplocal(iproc(ji,jj,jobs)) + 1
               END DO
            END DO
         END DO
      ENDIF

      ! Send local number of points and receive points on current domain

      CALL mpp_alltoall_int( 1, nplocal, npglobal )

      ! Allocate message parsing workspace

      itot = SUM(npglobal)

      ALLOCATE( &
         & igrdij_send(kptsi*kptsj*kobs*2), &
         & igrdij_recv(itot*2),             &
         & zsend(kpk,itot),                 &
         & zrecv(kpk,kptsi*kptsj*kobs)      &
         & )

      ! Pack buffers for list of points

      it = 0
      DO jp = 1, jpnij 
         DO jobs = 1, kobs
            DO jj = 1, kptsj
               DO ji = 1, kptsi
                  IF ( iproc(ji,jj,jobs) == jp ) THEN
                     it = it + 1
                     iorder(ji,jj,jobs) = it
                     igrdij_send(2*it-1) = kgrdi(ji,jj,jobs)
                     igrdij_send(2*it  ) = kgrdj(ji,jj,jobs)
                  ENDIF
               END DO
            END DO
         END DO
      END DO

      ! Send and recieve buffers for list of points

      CALL mpp_alltoallv_int( igrdij_send, kptsi*kptsj*kobs*2, nplocal(:)*2, &
         &                    igrdij_recv, itot*2, npglobal(:)*2 )

      ! Pack interpolation data to be sent

      DO ji = 1, itot
         ii = mi1(igrdij_recv(2*ji-1))
         ij = mj1(igrdij_recv(2*ji))
         DO jk = 1, kpk
            zsend(jk,ji) = pval(ii,ij,jk)
         END DO
      END DO

      ! Re-adjust sizes

      nplocal(:)  = kpk*nplocal(:)
      npglobal(:) = kpk*npglobal(:)


      ! Send and receive data for interpolation stencil

      CALL mpp_alltoallv_real( zsend, kpk*itot,             npglobal, &
         &                     zrecv, kpk*kptsi*kptsj*kobs, nplocal )

      ! Copy the received data into output data structure
      
      DO jobs = 1, kobs
         DO jj = 1, kptsj
            DO ji = 1, kptsi
               it = iorder(ji,jj,jobs)  
               DO jk = 1, kpk
                  pgval(ji,jj,jk,jobs) = zrecv(jk,it)
               END DO
            END DO
         END DO
      END DO

      ! Deallocate message parsing workspace

      DEALLOCATE( &
         & igrdij_send, &
         & igrdij_recv, &
         & zsend,       &
         & zrecv        &
         & )

   END SUBROUTINE obs_int_comm_3d_global
   
   SUBROUTINE obs_int_comm_3d_local( kptsi, kptsj, kobs, kpi, kpj, kpk, kgrdi, kgrdj, &
      &                              pval, pgval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE obs_int_comm_3d_global  ***
      !!          
      !! ** Purpose : Get 3D interpolation stencil (global version)
      !!
      !! ** Method  : On-demand communication where each processor send its
      !!              list of (i,j) of points to all processors and receive 
      !!              the corresponding values
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  08-02  (K. Mogensen)  Original code
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(IN) :: kptsi        ! Number of i horizontal points per stencil 
      INTEGER, INTENT(IN) :: kptsj        ! Number of j horizontal points per stencil
      INTEGER, INTENT(IN) :: kobs         ! Local number of observations
      INTEGER, INTENT(IN) :: kpi          ! Number of model points in i direction
      INTEGER, INTENT(IN) :: kpj          ! Number of model points in j direction
      INTEGER, INTENT(IN) :: kpk          ! Number of levels
      INTEGER, DIMENSION(kptsi,kptsj,kobs), INTENT(IN) :: &
         & kgrdi, &         ! i,j indicies for each stencil
         & kgrdj
      REAL(KIND=wp), DIMENSION(kpi,kpj,kpk), INTENT(IN) ::&
         & pval             ! Local 3D array to extract data from
      REAL(KIND=wp), DIMENSION(kptsi,kptsj,kpk,kobs), INTENT(OUT) ::&
         & pgval            ! Stencil at each point
      !! * Local declarations
      INTEGER ::  ji
      INTEGER ::  jj
      INTEGER ::  jk
      INTEGER ::  jobs

      ! Check valid points

      IF ( ( MAXVAL(kgrdi) > jpi ) .OR. ( MINVAL(kgrdi) < 1 ) .OR. &
         & ( MAXVAL(kgrdj) > jpj ) .OR. ( MINVAL(kgrdj) < 1 ) ) THEN
         
         CALL ctl_stop( 'Error in obs_int_comm_3d_local', &
            &           'Point outside local domain' )
 
      ENDIF

      ! Copy local data

      DO jobs = 1, kobs
         DO jj = 1, kptsj
            DO ji = 1, kptsi
               DO jk = 1, kpk
                  pgval(ji,jj,jk,jobs) = &
                     &            pval(kgrdi(ji,jj,jobs),kgrdj(ji,jj,jobs),jk)
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE obs_int_comm_3d_local

END MODULE obs_inter_sup
 
