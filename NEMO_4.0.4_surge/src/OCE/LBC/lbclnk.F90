MODULE lbclnk
   !!======================================================================
   !!                       ***  MODULE  lbclnk  ***
   !! NEMO        : lateral boundary conditions
   !!=====================================================================
   !! History :  OPA  ! 1997-06  (G. Madec)  Original code
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            3.2  ! 2009-03  (R. Benshila)  External north fold treatment  
   !!            3.5  ! 2012     (S.Mocavero, I. Epicoco)  optimization of BDY comm. via lbc_bdy_lnk and lbc_obc_lnk
   !!            3.4  ! 2012-12  (R. Bourdalle-Badie, G. Reffray)  add a C1D case  
   !!            3.6  ! 2015-06  (O. Tintó and M. Castrillo)  add lbc_lnk_multi  
   !!            4.0  ! 2017-03  (G. Madec) automatique allocation of array size (use with any 3rd dim size)
   !!             -   ! 2017-04  (G. Madec) remove duplicated routines (lbc_lnk_2d_9, lbc_lnk_2d_multiple, lbc_lnk_3d_gather)
   !!             -   ! 2017-05  (G. Madec) create generic.h90 files to generate all lbc and north fold routines
   !!----------------------------------------------------------------------
   !!           define the generic interfaces of lib_mpp routines
   !!----------------------------------------------------------------------
   !!   lbc_lnk       : generic interface for mpp_lnk_3d and mpp_lnk_2d routines defined in lib_mpp
   !!   lbc_bdy_lnk   : generic interface for mpp_lnk_bdy_2d and mpp_lnk_bdy_3d routines defined in lib_mpp
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE lib_mpp        ! distributed memory computing library
   USE lbcnfd         ! north fold
   USE in_out_manager ! I/O manager

   IMPLICIT NONE
   PRIVATE

   INTERFACE lbc_lnk
      MODULE PROCEDURE   mpp_lnk_2d      , mpp_lnk_3d      , mpp_lnk_4d
   END INTERFACE
   INTERFACE lbc_lnk_ptr
      MODULE PROCEDURE   mpp_lnk_2d_ptr  , mpp_lnk_3d_ptr  , mpp_lnk_4d_ptr
   END INTERFACE
   INTERFACE lbc_lnk_multi
      MODULE PROCEDURE   lbc_lnk_2d_multi, lbc_lnk_3d_multi, lbc_lnk_4d_multi
   END INTERFACE
   !
   INTERFACE lbc_lnk_icb
      MODULE PROCEDURE mpp_lnk_2d_icb
   END INTERFACE

   INTERFACE mpp_nfd
      MODULE PROCEDURE   mpp_nfd_2d    , mpp_nfd_3d    , mpp_nfd_4d
      MODULE PROCEDURE   mpp_nfd_2d_ptr, mpp_nfd_3d_ptr, mpp_nfd_4d_ptr
   END INTERFACE

   PUBLIC   lbc_lnk       ! ocean/ice lateral boundary conditions
   PUBLIC   lbc_lnk_multi ! modified ocean/ice lateral boundary conditions
   PUBLIC   lbc_lnk_icb   ! iceberg lateral boundary conditions

#if   defined key_mpp_mpi
!$AGRIF_DO_NOT_TREAT
   INCLUDE 'mpif.h'
!$AGRIF_END_DO_NOT_TREAT
#endif

   INTEGER, PUBLIC, PARAMETER ::   jpfillnothing = 1
   INTEGER, PUBLIC, PARAMETER ::   jpfillcst     = 2
   INTEGER, PUBLIC, PARAMETER ::   jpfillcopy    = 3
   INTEGER, PUBLIC, PARAMETER ::   jpfillperio   = 4
   INTEGER, PUBLIC, PARAMETER ::   jpfillmpi     = 5

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!                   ***   load_ptr_(2,3,4)d   ***
   !!
   !!   * Dummy Argument :
   !!       in    ==>   ptab       ! array to be loaded (2D, 3D or 4D)
   !!                   cd_nat     ! nature of pt2d array grid-points
   !!                   psgn       ! sign used across the north fold boundary
   !!       inout <=>   ptab_ptr   ! array of 2D, 3D or 4D pointers
   !!                   cdna_ptr   ! nature of ptab array grid-points
   !!                   psgn_ptr   ! sign used across the north fold boundary
   !!                   kfld       ! number of elements that has been attributed
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!                  ***   lbc_lnk_(2,3,4)d_multi   ***
   !!                     ***   load_ptr_(2,3,4)d   ***
   !!
   !!   * Argument : dummy argument use in lbc_lnk_multi_... routines
   !!
   !!----------------------------------------------------------------------

#  define DIM_2d
#     define ROUTINE_LOAD           load_ptr_2d
#     define ROUTINE_MULTI          lbc_lnk_2d_multi
#     include "lbc_lnk_multi_generic.h90"
#     undef ROUTINE_MULTI
#     undef ROUTINE_LOAD
#  undef DIM_2d

#  define DIM_3d
#     define ROUTINE_LOAD           load_ptr_3d
#     define ROUTINE_MULTI          lbc_lnk_3d_multi
#     include "lbc_lnk_multi_generic.h90"
#     undef ROUTINE_MULTI
#     undef ROUTINE_LOAD
#  undef DIM_3d

#  define DIM_4d
#     define ROUTINE_LOAD           load_ptr_4d
#     define ROUTINE_MULTI          lbc_lnk_4d_multi
#     include "lbc_lnk_multi_generic.h90"
#     undef ROUTINE_MULTI
#     undef ROUTINE_LOAD
#  undef DIM_4d

   !!----------------------------------------------------------------------
   !!                   ***  routine mpp_lnk_(2,3,4)d  ***
   !!
   !!   * Argument : dummy argument use in mpp_lnk_... routines
   !!                ptab      :   array or pointer of arrays on which the boundary condition is applied
   !!                cd_nat    :   nature of array grid-points
   !!                psgn      :   sign used across the north fold boundary
   !!                kfld      :   optional, number of pt3d arrays
   !!                kfillmode :   optional, method to be use to fill the halos (see jpfill* variables)
   !!                pfillval  :   optional, background value (used with jpfillcopy)
   !!----------------------------------------------------------------------
   !
   !                       !==  2D array and array of 2D pointer  ==!
   !
#  define DIM_2d
#     define ROUTINE_LNK           mpp_lnk_2d
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     define MULTI
#     define ROUTINE_LNK           mpp_lnk_2d_ptr
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     undef MULTI
#  undef DIM_2d
   !
   !                       !==  3D array and array of 3D pointer  ==!
   !
#  define DIM_3d
#     define ROUTINE_LNK           mpp_lnk_3d
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     define MULTI
#     define ROUTINE_LNK           mpp_lnk_3d_ptr
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     undef MULTI
#  undef DIM_3d
   !
   !                       !==  4D array and array of 4D pointer  ==!
   !
#  define DIM_4d
#     define ROUTINE_LNK           mpp_lnk_4d
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     define MULTI
#     define ROUTINE_LNK           mpp_lnk_4d_ptr
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     undef MULTI
#  undef DIM_4d

   !!----------------------------------------------------------------------
   !!                   ***  routine mpp_nfd_(2,3,4)d  ***
   !!
   !!   * Argument : dummy argument use in mpp_nfd_... routines
   !!                ptab      :   array or pointer of arrays on which the boundary condition is applied
   !!                cd_nat    :   nature of array grid-points
   !!                psgn      :   sign used across the north fold boundary
   !!                kfld      :   optional, number of pt3d arrays
   !!                kfillmode :   optional, method to be use to fill the halos (see jpfill* variables)
   !!                pfillval  :   optional, background value (used with jpfillcopy)
   !!----------------------------------------------------------------------
   !
   !                       !==  2D array and array of 2D pointer  ==!
   !
#  define DIM_2d
#     define ROUTINE_NFD           mpp_nfd_2d
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     define MULTI
#     define ROUTINE_NFD           mpp_nfd_2d_ptr
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     undef MULTI
#  undef DIM_2d
   !
   !                       !==  3D array and array of 3D pointer  ==!
   !
#  define DIM_3d
#     define ROUTINE_NFD           mpp_nfd_3d
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     define MULTI
#     define ROUTINE_NFD           mpp_nfd_3d_ptr
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     undef MULTI
#  undef DIM_3d
   !
   !                       !==  4D array and array of 4D pointer  ==!
   !
#  define DIM_4d
#     define ROUTINE_NFD           mpp_nfd_4d
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     define MULTI
#     define ROUTINE_NFD           mpp_nfd_4d_ptr
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     undef MULTI
#  undef DIM_4d


   !!======================================================================



   SUBROUTINE mpp_lbc_north_icb( pt2d, cd_type, psgn, kextj)
      !!---------------------------------------------------------------------
      !!                   ***  routine mpp_lbc_north_icb  ***
      !!
      !! ** Purpose :   Ensure proper north fold horizontal bondary condition
      !!              in mpp configuration in case of jpn1 > 1 and for 2d
      !!              array with outer extra halo
      !!
      !! ** Method  :   North fold condition and mpp with more than one proc
      !!              in i-direction require a specific treatment. We gather
      !!              the 4+kextj northern lines of the global domain on 1
      !!              processor and apply lbc north-fold on this sub array.
      !!              Then we scatter the north fold array back to the processors.
      !!              This routine accounts for an extra halo with icebergs
      !!              and assumes ghost rows and columns have been suppressed.
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pt2d     ! 2D array with extra halo
      CHARACTER(len=1)        , INTENT(in   ) ::   cd_type  ! nature of pt3d grid-points
      !                                                     !   = T ,  U , V , F or W -points
      REAL(wp)                , INTENT(in   ) ::   psgn     ! = -1. the sign change across the
      !!                                                    ! north fold, =  1. otherwise
      INTEGER                 , INTENT(in   ) ::   kextj    ! Extra halo width at north fold
      !
      INTEGER ::   ji, jj, jr
      INTEGER ::   ierr, itaille, ildi, ilei, iilb
      INTEGER ::   ipj, ij, iproc
      !
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE  ::  ztab_e, znorthloc_e
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE  ::  znorthgloio_e
      !!----------------------------------------------------------------------
#if defined key_mpp_mpi
      !
      ipj=4
      ALLOCATE(        ztab_e(jpiglo, 1-kextj:ipj+kextj)       ,       &
     &            znorthloc_e(jpimax, 1-kextj:ipj+kextj)       ,       &
     &          znorthgloio_e(jpimax, 1-kextj:ipj+kextj,jpni)    )
      !
      ztab_e(:,:)      = 0._wp
      znorthloc_e(:,:) = 0._wp
      !
      ij = 1 - kextj
      ! put the last ipj+2*kextj lines of pt2d into znorthloc_e 
      DO jj = jpj - ipj + 1 - kextj , jpj + kextj
         znorthloc_e(1:jpi,ij)=pt2d(1:jpi,jj)
         ij = ij + 1
      END DO
      !
      itaille = jpimax * ( ipj + 2*kextj )
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      CALL MPI_ALLGATHER( znorthloc_e(1,1-kextj)    , itaille, MPI_DOUBLE_PRECISION,    &
         &                znorthgloio_e(1,1-kextj,1), itaille, MPI_DOUBLE_PRECISION,    &
         &                ncomm_north, ierr )
      !
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      DO jr = 1, ndim_rank_north            ! recover the global north array
         iproc = nrank_north(jr) + 1
         ildi = nldit (iproc)
         ilei = nleit (iproc)
         iilb = nimppt(iproc)
         DO jj = 1-kextj, ipj+kextj
            DO ji = ildi, ilei
               ztab_e(ji+iilb-1,jj) = znorthgloio_e(ji,jj,jr)
            END DO
         END DO
      END DO

      ! 2. North-Fold boundary conditions
      ! ----------------------------------
      CALL lbc_nfd( ztab_e(:,1-kextj:ipj+kextj), cd_type, psgn, kextj )

      ij = 1 - kextj
      !! Scatter back to pt2d
      DO jj = jpj - ipj + 1 - kextj , jpj + kextj
         DO ji= 1, jpi
            pt2d(ji,jj) = ztab_e(ji+nimpp-1,ij)
         END DO
         ij  = ij +1
      END DO
      !
      DEALLOCATE( ztab_e, znorthloc_e, znorthgloio_e )
      !
#endif
   END SUBROUTINE mpp_lbc_north_icb


   SUBROUTINE mpp_lnk_2d_icb( cdname, pt2d, cd_type, psgn, kexti, kextj )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_2d_icb  ***
      !!
      !! ** Purpose :   Message passing management for 2d array (with extra halo for icebergs)
      !!                This routine receives a (1-kexti:jpi+kexti,1-kexti:jpj+kextj)
      !!                array (usually (0:jpi+1, 0:jpj+1)) from lbc_lnk_icb calls.
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    jpi    : first dimension of the local subdomain
      !!                    jpj    : second dimension of the local subdomain
      !!                    kexti  : number of columns for extra outer halo
      !!                    kextj  : number of rows for extra outer halo
      !!                    nbondi : mark for "east-west local boundary"
      !!                    nbondj : mark for "north-south local boundary"
      !!                    noea   : number for local neighboring processors
      !!                    nowe   : number for local neighboring processors
      !!                    noso   : number for local neighboring processors
      !!                    nono   : number for local neighboring processors
      !!----------------------------------------------------------------------
      CHARACTER(len=*)                                        , INTENT(in   ) ::   cdname      ! name of the calling subroutine
      REAL(wp), DIMENSION(1-kexti:jpi+kexti,1-kextj:jpj+kextj), INTENT(inout) ::   pt2d     ! 2D array with extra halo
      CHARACTER(len=1)                                        , INTENT(in   ) ::   cd_type  ! nature of ptab array grid-points
      REAL(wp)                                                , INTENT(in   ) ::   psgn     ! sign used across the north fold
      INTEGER                                                 , INTENT(in   ) ::   kexti    ! extra i-halo width
      INTEGER                                                 , INTENT(in   ) ::   kextj    ! extra j-halo width
      !
      INTEGER  ::   jl   ! dummy loop indices
      INTEGER  ::   imigr, iihom, ijhom        ! local integers
      INTEGER  ::   ipreci, iprecj             !   -       -
      INTEGER  ::   ml_req1, ml_req2, ml_err   ! for key_mpi_isend
      INTEGER, DIMENSION(MPI_STATUS_SIZE) ::   ml_stat   ! for key_mpi_isend
      !!
      REAL(wp), DIMENSION(1-kexti:jpi+kexti,nn_hls+kextj,2) ::   r2dns, r2dsn
      REAL(wp), DIMENSION(1-kextj:jpj+kextj,nn_hls+kexti,2) ::   r2dwe, r2dew
      !!----------------------------------------------------------------------

      ipreci = nn_hls + kexti      ! take into account outer extra 2D overlap area
      iprecj = nn_hls + kextj

      IF( narea == 1 .AND. numcom == -1 ) CALL mpp_report( cdname, 1, 1, 1, ld_lbc = .TRUE. )

      ! 1. standard boundary treatment
      ! ------------------------------
      ! Order matters Here !!!!
      !
      !                                      ! East-West boundaries
      !                                           !* Cyclic east-west
      IF( l_Iperio ) THEN
         pt2d(1-kexti:     1   ,:) = pt2d(jpim1-kexti: jpim1 ,:)       ! east
         pt2d(  jpi  :jpi+kexti,:) = pt2d(     2     :2+kexti,:)       ! west
         !
      ELSE                                        !* closed
         IF( .NOT. cd_type == 'F' )   pt2d(  1-kexti   :nn_hls   ,:) = 0._wp    ! east except at F-point
                                      pt2d(jpi-nn_hls+1:jpi+kexti,:) = 0._wp    ! west
      ENDIF
      !                                      ! North-South boundaries
      IF( l_Jperio ) THEN                         !* cyclic (only with no mpp j-split)
         pt2d(:,1-kextj:     1   ) = pt2d(:,jpjm1-kextj:  jpjm1)       ! north
         pt2d(:,  jpj  :jpj+kextj) = pt2d(:,     2     :2+kextj)       ! south
      ELSE                                        !* closed
         IF( .NOT. cd_type == 'F' )   pt2d(:,  1-kextj   :nn_hls   ) = 0._wp    ! north except at F-point
                                      pt2d(:,jpj-nn_hls+1:jpj+kextj) = 0._wp    ! south
      ENDIF
      !

      ! north fold treatment
      ! -----------------------
      IF( npolj /= 0 ) THEN
         !
         SELECT CASE ( jpni )
                   CASE ( 1 )     ;   CALL lbc_nfd          ( pt2d(1:jpi,1-kextj:jpj+kextj), cd_type, psgn, kextj )
                   CASE DEFAULT   ;   CALL mpp_lbc_north_icb( pt2d(1:jpi,1:jpj+kextj), cd_type, psgn, kextj )
         END SELECT
         !
      ENDIF

      ! 2. East and west directions exchange
      ! ------------------------------------
      ! we play with the neigbours AND the row number because of the periodicity
      !
      SELECT CASE ( nbondi )      ! Read Dirichlet lateral conditions
      CASE ( -1, 0, 1 )                ! all exept 2 (i.e. close case)
         iihom = jpi-nreci-kexti
         DO jl = 1, ipreci
            r2dew(:,jl,1) = pt2d(nn_hls+jl,:)
            r2dwe(:,jl,1) = pt2d(iihom +jl,:)
         END DO
      END SELECT
      !
      !                           ! Migrations
      imigr = ipreci * ( jpj + 2*kextj )
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         CALL mppsend( 2, r2dwe(1-kextj,1,1), imigr, noea, ml_req1 )
         CALL mpprecv( 1, r2dew(1-kextj,1,2), imigr, noea )
         CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 1, r2dew(1-kextj,1,1), imigr, nowe, ml_req1 )
         CALL mppsend( 2, r2dwe(1-kextj,1,1), imigr, noea, ml_req2 )
         CALL mpprecv( 1, r2dew(1-kextj,1,2), imigr, noea )
         CALL mpprecv( 2, r2dwe(1-kextj,1,2), imigr, nowe )
         CALL mpi_wait(ml_req1,ml_stat,ml_err)
         CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 1, r2dew(1-kextj,1,1), imigr, nowe, ml_req1 )
         CALL mpprecv( 2, r2dwe(1-kextj,1,2), imigr, nowe )
         CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT
      !
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      !                           ! Write Dirichlet lateral conditions
      iihom = jpi - nn_hls
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         DO jl = 1, ipreci
            pt2d(iihom+jl,:) = r2dew(:,jl,2)
         END DO
      CASE ( 0 )
         DO jl = 1, ipreci
            pt2d(jl-kexti,:) = r2dwe(:,jl,2)
            pt2d(iihom+jl,:) = r2dew(:,jl,2)
         END DO
      CASE ( 1 )
         DO jl = 1, ipreci
            pt2d(jl-kexti,:) = r2dwe(:,jl,2)
         END DO
      END SELECT


      ! 3. North and south directions
      ! -----------------------------
      ! always closed : we play only with the neigbours
      !
      IF( nbondj /= 2 ) THEN      ! Read Dirichlet lateral conditions
         ijhom = jpj-nrecj-kextj
         DO jl = 1, iprecj
            r2dsn(:,jl,1) = pt2d(:,ijhom +jl)
            r2dns(:,jl,1) = pt2d(:,nn_hls+jl)
         END DO
      ENDIF
      !
      !                           ! Migrations
      imigr = iprecj * ( jpi + 2*kexti )
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL mppsend( 4, r2dsn(1-kexti,1,1), imigr, nono, ml_req1 )
         CALL mpprecv( 3, r2dns(1-kexti,1,2), imigr, nono )
         CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 3, r2dns(1-kexti,1,1), imigr, noso, ml_req1 )
         CALL mppsend( 4, r2dsn(1-kexti,1,1), imigr, nono, ml_req2 )
         CALL mpprecv( 3, r2dns(1-kexti,1,2), imigr, nono )
         CALL mpprecv( 4, r2dsn(1-kexti,1,2), imigr, noso )
         CALL mpi_wait(ml_req1,ml_stat,ml_err)
         CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 3, r2dns(1-kexti,1,1), imigr, noso, ml_req1 )
         CALL mpprecv( 4, r2dsn(1-kexti,1,2), imigr, noso )
         CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT
      !
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      !                           ! Write Dirichlet lateral conditions
      ijhom = jpj - nn_hls
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         DO jl = 1, iprecj
            pt2d(:,ijhom+jl) = r2dns(:,jl,2)
         END DO
      CASE ( 0 )
         DO jl = 1, iprecj
            pt2d(:,jl-kextj) = r2dsn(:,jl,2)
            pt2d(:,ijhom+jl) = r2dns(:,jl,2)
         END DO
      CASE ( 1 )
         DO jl = 1, iprecj
            pt2d(:,jl-kextj) = r2dsn(:,jl,2)
         END DO
      END SELECT
      !
   END SUBROUTINE mpp_lnk_2d_icb
   
END MODULE lbclnk

