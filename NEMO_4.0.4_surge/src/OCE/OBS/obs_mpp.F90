MODULE obs_mpp
   !!======================================================================
   !!                       ***  MODULE obs_mpp  ***
   !! Observation diagnostics: Various MPP support routines
   !!======================================================================
   !! History :  2.0  ! 2006-03  (K. Mogensen)  Original code
   !!             -   ! 2006-05  (K. Mogensen)  Reformatted
   !!             -   ! 2008-01  (K. Mogensen)  add mpp_global_max
   !!            3.6  ! 2015-01  (J. Waters) obs_mpp_find_obs_proc 
   !!                            rewritten to avoid global arrays
   !!----------------------------------------------------------------------
#  define mpivar mpi_double_precision
   !!----------------------------------------------------------------------
   !! obs_mpp_bcast_integer : Broadcast an integer array from a processor to all processors
   !! obs_mpp_max_integer   : Find maximum on all processors of each value in an integer on all processors
   !! obs_mpp_find_obs_proc : Find processors which should hold the observations, avoiding global arrays
   !! obs_mpp_sum_integers  : Sum an integer array from all processors
   !! obs_mpp_sum_integer   : Sum an integer from all processors
   !!----------------------------------------------------------------------
   USE dom_oce, ONLY :   nproc, mig, mjg   ! Ocean space and time domain variables
   USE mpp_map, ONLY :   mppmap
   USE in_out_manager
#if defined key_mpp_mpi
   USE lib_mpp, ONLY :   mpi_comm_oce      ! MPP library
#endif
   IMPLICIT NONE
   PRIVATE

   PUBLIC obs_mpp_bcast_integer, & !: Broadcast an integer array from a proc to all procs
      &   obs_mpp_max_integer,   & !: Find maximum across processors in an integer array
      &   obs_mpp_find_obs_proc, & !: Find processors which should hold the observations
      &   obs_mpp_sum_integers,  & !: Sum an integer array from all processors
      &   obs_mpp_sum_integer,   & !: Sum an integer from all processors
      &   mpp_alltoall_int,      &
      &   mpp_alltoallv_int,     &
      &   mpp_alltoallv_real,    &
      &   mpp_global_max

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE obs_mpp_bcast_integer( kvals, kno, kroot )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE obs_mpp_bcast_integer ***
      !!          
      !! ** Purpose : Send array kvals to all processors
      !!
      !! ** Method  : MPI broadcast
      !!
      !! ** Action  : This does only work for MPI. 
      !!              MPI_COMM_OCE needs to be replace for OASIS4.!
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in   ) ::   kno     ! Number of elements in array
      INTEGER                , INTENT(in   ) ::   kroot   ! Processor to send data
      INTEGER, DIMENSION(kno), INTENT(inout) ::   kvals   ! Array to send on kroot, receive for non-kroot
      !
#if defined key_mpp_mpi
      !
      INTEGER :: ierr 
      !
INCLUDE 'mpif.h'
      !!----------------------------------------------------------------------

      ! Call the MPI library to broadcast data
      CALL mpi_bcast( kvals, kno, mpi_integer,  &
         &            kroot, mpi_comm_oce, ierr )
#else
      ! no MPI: empty routine
#endif
      !
   END SUBROUTINE obs_mpp_bcast_integer

  
   SUBROUTINE obs_mpp_max_integer( kvals, kno )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE obs_mpp_bcast_integer ***
      !!          
      !! ** Purpose : Find maximum across processors in an integer array.
      !!
      !! ** Method  : MPI all reduce.
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!              MPI_COMM_OCE needs to be replace for OASIS4.!
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in   ) ::   kno     ! Number of elements in array
      INTEGER, DIMENSION(kno), INTENT(inout) ::   kvals   ! Array to send on kroot, receive for non-kroot  
      !
#if defined key_mpp_mpi
      !
      INTEGER :: ierr 
      INTEGER, DIMENSION(kno) ::   ivals
      !
INCLUDE 'mpif.h'
      !!----------------------------------------------------------------------

      ! Call the MPI library to find the maximum across processors
      CALL mpi_allreduce( kvals, ivals, kno, mpi_integer,   &
         &                mpi_max, mpi_comm_oce, ierr )
      kvals(:) = ivals(:)
#else
      ! no MPI: empty routine
#endif
   END SUBROUTINE obs_mpp_max_integer


   SUBROUTINE obs_mpp_find_obs_proc( kobsp,kno )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE obs_mpp_find_obs_proc  ***
      !!         
      !! ** Purpose : From the array kobsp containing the results of the
      !!              grid search on each processor the processor return a
      !!              decision of which processors should hold the observation.
      !!
      !! ** Method  : Synchronize the processor number for each obs using
      !!              obs_mpp_max_integer. If an observation exists on two 
      !!              processors it will be allocated to the lower numbered
      !!              processor.
      !!
      !! ** Action  : This does only work for MPI.
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in   ) ::   kno
      INTEGER, DIMENSION(kno), INTENT(inout) ::   kobsp
      !
#if defined key_mpp_mpi
      !
      !
      INTEGER :: ji, isum
      INTEGER, DIMENSION(kno) ::   iobsp
      !!
      !!

      iobsp(:)=kobsp(:)

      WHERE( iobsp(:) == -1 )
         iobsp(:) = 9999999
      END WHERE

      iobsp(:)=-1*iobsp(:)

      CALL obs_mpp_max_integer( iobsp, kno )

      kobsp(:)=-1*iobsp(:)

      isum=0
      DO ji = 1, kno
         IF ( kobsp(ji) == 9999999 ) THEN
            isum=isum+1
            kobsp(ji)=-1
         ENDIF
      ENDDO


      IF ( isum > 0 ) THEN
         IF (lwp) WRITE(numout,*) isum, ' observations failed the grid search.'
         IF (lwp) WRITE(numout,*)'If ln_grid_search_lookup=.TRUE., try reducing grid_search_res'
      ENDIF

#else
      ! no MPI: empty routine
#endif     
      
   END SUBROUTINE obs_mpp_find_obs_proc


   SUBROUTINE obs_mpp_sum_integers( kvalsin, kvalsout, kno )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE obs_mpp_sum_integers ***
      !!          
      !! ** Purpose : Sum an integer array.
      !!
      !! ** Method  : MPI all reduce.
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                , INTENT(in   ) :: kno
      INTEGER, DIMENSION(kno), INTENT(in   ) ::   kvalsin
      INTEGER, DIMENSION(kno), INTENT(  out) ::   kvalsout
      !
#if defined key_mpp_mpi
      !
      INTEGER :: ierr
      !
INCLUDE 'mpif.h'
      !!----------------------------------------------------------------------
      !
      !-----------------------------------------------------------------------
      ! Call the MPI library to find the sum across processors
      !-----------------------------------------------------------------------
      CALL mpi_allreduce( kvalsin, kvalsout, kno, mpi_integer, &
         &                mpi_sum, mpi_comm_oce, ierr )
#else
      !-----------------------------------------------------------------------
      ! For no-MPP just return input values
      !-----------------------------------------------------------------------
      kvalsout(:) = kvalsin(:)
#endif
      !
   END SUBROUTINE obs_mpp_sum_integers


   SUBROUTINE obs_mpp_sum_integer( kvalin, kvalout )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE obs_mpp_sum_integers ***
      !!          
      !! ** Purpose : Sum a single integer
      !!
      !! ** Method  : MPI all reduce.
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kvalin
      INTEGER, INTENT(  out) ::   kvalout
      !
#if defined key_mpp_mpi
      !
      INTEGER :: ierr
      !
INCLUDE 'mpif.h'
      !!----------------------------------------------------------------------
      !
      !-----------------------------------------------------------------------
      ! Call the MPI library to find the sum across processors
      !-----------------------------------------------------------------------
      CALL mpi_allreduce( kvalin, kvalout, 1, mpi_integer,   &
         &                mpi_sum, mpi_comm_oce, ierr )
#else
      !-----------------------------------------------------------------------
      ! For no-MPP just return input values
      !-----------------------------------------------------------------------
      kvalout = kvalin
#endif
      !
   END SUBROUTINE obs_mpp_sum_integer


   SUBROUTINE mpp_global_max( pval )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE mpp_global_or ***
      !!          
      !! ** Purpose : Get the maximum value across processors for a global 
      !!              real array
      !!
      !! ** Method  : MPI allreduce
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      REAL(KIND=wp), DIMENSION(jpiglo,jpjglo), INTENT(inout) ::   pval
      !
      INTEGER :: ierr
      !
#if defined key_mpp_mpi
      !
INCLUDE 'mpif.h'
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE ::   zcp
      !!----------------------------------------------------------------------

      ! Copy data for input to MPI

      ALLOCATE( &
         & zcp(jpiglo,jpjglo) &
         & )
      zcp(:,:) = pval(:,:)

      ! Call the MPI library to find the coast lines globally

      CALL mpi_allreduce( zcp, pval, jpiglo*jpjglo, mpivar, &
         &                mpi_max, mpi_comm_oce, ierr )

      DEALLOCATE( &
         & zcp &
         & )

#else
      ! no MPI: empty routine
#endif
      !
   END SUBROUTINE mpp_global_max


   SUBROUTINE mpp_alltoall_int( kno, kvalsin, kvalsout )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE mpp_allgatherv ***
      !!          
      !! ** Purpose : all to all.
      !!
      !! ** Method  : MPI alltoall
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                      , INTENT(in   ) ::   kno
      INTEGER, DIMENSION(kno*jpnij), INTENT(in   ) ::   kvalsin
      INTEGER, DIMENSION(kno*jpnij), INTENT(  out) ::   kvalsout
      !!
      INTEGER :: ierr
      !
#if defined key_mpp_mpi
      !
INCLUDE 'mpif.h'
      !-----------------------------------------------------------------------
      ! Call the MPI library to do the all to all operation of the data
      !-----------------------------------------------------------------------
      CALL mpi_alltoall( kvalsin,  kno, mpi_integer, &
         &               kvalsout, kno, mpi_integer, &
         &               mpi_comm_oce, ierr )
#else
      !-----------------------------------------------------------------------
      ! For no-MPP just return input values
      !-----------------------------------------------------------------------
      kvalsout = kvalsin
#endif
      !
   END SUBROUTINE mpp_alltoall_int


   SUBROUTINE mpp_alltoallv_int( kvalsin, knoin , kinv , kvalsout,   &
      &                                   knoout, koutv )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE mpp_alltoallv_int ***
      !!          
      !! ** Purpose : all to all (integer version).
      !!
      !! ** Method  : MPI alltoall
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in) :: knoin
      INTEGER                   , INTENT(in) :: knoout
      INTEGER, DIMENSION(jpnij)                 ::   kinv, koutv
      INTEGER, DIMENSION(knoin) , INTENT(in   ) ::   kvalsin
      INTEGER, DIMENSION(knoout), INTENT(  out) ::   kvalsout
      !!
      INTEGER :: ierr
      INTEGER :: jproc
      !
#if defined key_mpp_mpi
      !
INCLUDE 'mpif.h'
      INTEGER, DIMENSION(jpnij) ::   irdsp, isdsp
      !-----------------------------------------------------------------------
      ! Compute displacements
      !-----------------------------------------------------------------------
      irdsp(1) = 0
      isdsp(1) = 0
      DO jproc = 2, jpnij
         isdsp(jproc) = isdsp(jproc-1) + kinv(jproc-1)
         irdsp(jproc) = irdsp(jproc-1) + koutv(jproc-1)
      END DO
      !-----------------------------------------------------------------------
      ! Call the MPI library to do the all to all operation of the data
      !-----------------------------------------------------------------------
      CALL mpi_alltoallv( kvalsin,  kinv,  isdsp, mpi_integer, &
         &                kvalsout, koutv, irdsp, mpi_integer, &
         &                mpi_comm_oce, ierr )
#else
      !-----------------------------------------------------------------------
      ! For no-MPP just return input values
      !-----------------------------------------------------------------------
      kvalsout = kvalsin
#endif
      !
   END SUBROUTINE mpp_alltoallv_int


   SUBROUTINE mpp_alltoallv_real( pvalsin, knoin , kinv , pvalsout,   &
      &                                    knoout, koutv )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE mpp_alltoallv_real ***
      !!          
      !! ** Purpose : all to all (integer version).
      !!
      !! ** Method  : MPI alltoall
      !!
      !! ** Action  : This does only work for MPI. 
      !!              It does not work for SHMEM.
      !!
      !! References : http://www.mpi-forum.org
      !!----------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) :: knoin
      INTEGER                    , INTENT(in   ) :: knoout
      INTEGER , DIMENSION(jpnij)                 ::   kinv, koutv
      REAL(wp), DIMENSION(knoin) , INTENT(in   ) ::   pvalsin
      REAL(wp), DIMENSION(knoout), INTENT(  out) ::   pvalsout
      !!
      INTEGER :: ierr
      INTEGER :: jproc
      !
#if defined key_mpp_mpi
      !
INCLUDE 'mpif.h'
      INTEGER, DIMENSION(jpnij) ::   irdsp, isdsp
      !!----------------------------------------------------------------------
      !
      !-----------------------------------------------------------------------
      ! Compute displacements
      !-----------------------------------------------------------------------
      irdsp(1) = 0
      isdsp(1) = 0
      DO jproc = 2, jpnij
         isdsp(jproc) = isdsp(jproc-1) + kinv(jproc-1)
         irdsp(jproc) = irdsp(jproc-1) + koutv(jproc-1)
      END DO
      !-----------------------------------------------------------------------
      ! Call the MPI library to do the all to all operation of the data
      !-----------------------------------------------------------------------
      CALL mpi_alltoallv( pvalsin,  kinv,  isdsp, mpivar, &
         &                pvalsout, koutv, irdsp, mpivar, &
         &                mpi_comm_oce, ierr )
#else
      !-----------------------------------------------------------------------
      ! For no-MPP just return input values
      !-----------------------------------------------------------------------
      pvalsout = pvalsin
#endif
      !
   END SUBROUTINE mpp_alltoallv_real

   !!======================================================================
END MODULE obs_mpp
