MODULE obs_surf_def
   !!=====================================================================
   !!                       ***  MODULE  obs_surf_def  ***
   !! Observation diagnostics: Storage handling for surface observation
   !!                          arrays and additional flags etc.
   !!                          This module only defines the data type and
   !!                          operations on the data type. There is no
   !!                          actual data in the module.
   !!===================================================================== 

   !!----------------------------------------------------------------------
   !!   obs_surf            : F90 type containing the surface information
   !!   obs_surf_alloc      : Allocates surface data arrays
   !!   obs_surf_dealloc    : Deallocates surface data arrays
   !!   obs_surf_compress   : Extract sub-information from a obs_surf type
   !!                         to a new obs_surf type
   !!   obs_surf_decompress : Reinsert sub-information from a obs_surf type
   !!                         into the original obs_surf type
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind, ONLY : & ! Precision variables
      & wp         
   USE obs_mpp, ONLY : &  ! MPP tools 
      obs_mpp_sum_integer

   IMPLICIT NONE

   !! * Routine/type accessibility
   PRIVATE

   PUBLIC &
      & obs_surf,           &
      & obs_surf_alloc,     &
      & obs_surf_dealloc,   &
      & obs_surf_compress,  &
      & obs_surf_decompress

   !! * Type definition for surface observation type

   TYPE obs_surf

      ! Bookkeeping

      INTEGER :: nsurf      !: Local number of surface data within window
      INTEGER :: nsurfmpp   !: Global number of surface data within window
      INTEGER :: nvar       !: Number of variables at observation points
      INTEGER :: nextra     !: Number of extra fields at observation points
      INTEGER :: nstp       !: Number of time steps
      INTEGER :: npi        !: Number of 3D grid points
      INTEGER :: npj
      INTEGER :: nsurfup    !: Observation counter used in obs_oper
      INTEGER :: nrec       !: Number of surface observation records in window

      ! Arrays with size equal to the number of surface observations

      INTEGER, POINTER, DIMENSION(:) :: &
         & mi,   &        !: i-th grid coord. for interpolating to surface observation
         & mj,   &        !: j-th grid coord. for interpolating to surface observation
         & mt,   &        !: time record number for gridded data
         & nsidx,&        !: Surface observation number
         & nsfil,&        !: Surface observation number in file
         & nyea, &        !: Year of surface observation
         & nmon, &        !: Month of surface observation 
         & nday, &        !: Day of surface observation
         & nhou, &        !: Hour of surface observation
         & nmin, &        !: Minute of surface observation
         & mstp, &        !: Time step nearest to surface observation
         & nqc,  &        !: Surface observation qc flag
         & ntyp           !: Type of surface observation product

      CHARACTER(len=8), POINTER, DIMENSION(:) :: &
         & cvars          !: Variable names

      CHARACTER(LEN=8), POINTER, DIMENSION(:) :: &
         & cwmo           !: WMO indentifier
         
      REAL(KIND=wp), POINTER, DIMENSION(:) :: &
         & rlam, &        !: Longitude coordinate of surface observation
         & rphi           !: Latitude coordinate of surface observation

      REAL(KIND=wp), POINTER, DIMENSION(:,:) :: &
         & robs, &        !: Surface observation 
         & rmod           !: Model counterpart of the surface observation vector

      REAL(KIND=wp), POINTER, DIMENSION(:,:) :: &
         & rext           !: Extra fields interpolated to observation points

      REAL(KIND=wp), POINTER, DIMENSION(:,:) :: &
         & vdmean         !: Time averaged of model field

      ! Arrays with size equal to the number of time steps in the window

      INTEGER, POINTER, DIMENSION(:) :: &
         & nsstp,     &   !: Local number of surface observations per time step
         & nsstpmpp       !: Global number of surface observations per time step

      ! Arrays with size equal to the number of observation records in the window
      INTEGER, POINTER, DIMENSION(:) :: &
         & mrecstp   ! Time step of the records

      ! Arrays used to store source indices when 
      ! compressing obs_surf derived types
      
      ! Array with size nsurf

      INTEGER, POINTER, DIMENSION(:) :: &
         & nsind          !: Source indices of surface data in compressed data

      ! Is this a gridded product?
     
      LOGICAL :: lgrid

   END TYPE obs_surf

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS
   
   SUBROUTINE obs_surf_alloc( surf, ksurf, kvar, kextra, kstp, kpi, kpj )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE obs_surf_alloc  ***
      !!                      
      !! ** Purpose : - Allocate data for surface data arrays
      !! 
      !! ** Method  : - Fortran-90 dynamic arrays
      !!
      !! History :
      !!        !  07-03  (K. Mogensen, A. Weaver, E. Remy, S. Ricci) original
      !!----------------------------------------------------------------------
      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) ::  surf      ! Surface data to be allocated
      INTEGER, INTENT(IN) :: ksurf   ! Number of surface observations
      INTEGER, INTENT(IN) :: kvar    ! Number of surface variables
      INTEGER, INTENT(IN) :: kextra  ! Number of extra fields at observation points
      INTEGER, INTENT(IN) :: kstp    ! Number of time steps
      INTEGER, INTENT(IN) :: kpi     ! Number of 3D grid points
      INTEGER, INTENT(IN) :: kpj

      !!* Local variables
      INTEGER :: ji
      INTEGER :: jvar

      ! Set bookkeeping variables

      surf%nsurf    = ksurf
      surf%nsurfmpp = 0
      surf%nextra   = kextra
      surf%nvar     = kvar
      surf%nstp     = kstp
      surf%npi      = kpi
      surf%npj      = kpj

      ! Allocate arrays of size number of variables

      ALLOCATE( &
         & surf%cvars(kvar)    &
         & )

      DO jvar = 1, kvar
         surf%cvars(jvar) = "NotSet"
      END DO
      
      ! Allocate arrays of number of surface data size

      ALLOCATE( &
         & surf%mi(ksurf),      &
         & surf%mj(ksurf),      &
         & surf%mt(ksurf),      &
         & surf%nsidx(ksurf),   &
         & surf%nsfil(ksurf),   &
         & surf%nyea(ksurf),    &
         & surf%nmon(ksurf),    &
         & surf%nday(ksurf),    &
         & surf%nhou(ksurf),    &
         & surf%nmin(ksurf),    &
         & surf%mstp(ksurf),    &
         & surf%nqc(ksurf),     &
         & surf%ntyp(ksurf),    &
         & surf%cwmo(ksurf),    &
         & surf%rlam(ksurf),    &
         & surf%rphi(ksurf),    &
         & surf%nsind(ksurf)    &
         & )

      surf%mt(:) = -1


      ! Allocate arrays of number of surface data size * number of variables

      ALLOCATE( & 
         & surf%robs(ksurf,kvar), &
         & surf%rmod(ksurf,kvar)  &
         & )   

      ! Allocate arrays of number of extra fields at observation points

      ALLOCATE( & 
         & surf%rext(ksurf,kextra) &
         & )

      surf%rext(:,:) = 0.0_wp 

      ! Allocate arrays of number of time step size

      ALLOCATE( &
         & surf%nsstp(kstp),     &
         & surf%nsstpmpp(kstp)   &
         & )

      ! Allocate arrays of size number of grid points

      ALLOCATE( &
         & surf%vdmean(kpi,kpj) &
         & )

      ! Set defaults for compression indices
      
      DO ji = 1, ksurf
         surf%nsind(ji) = ji
      END DO
      
      ! Set defaults for number of observations per time step

      surf%nsstp(:)     = 0
      surf%nsstpmpp(:)  = 0

      ! Set the observation counter used in obs_oper

      surf%nsurfup     = 0
      
      ! Not gridded by default
          
      surf%lgrid       = .FALSE.
              
   END SUBROUTINE obs_surf_alloc

   SUBROUTINE obs_surf_dealloc( surf )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE obs_surf_dealloc  ***
      !!                      
      !! ** Purpose : - Deallocate data for surface data arrays
      !! 
      !! ** Method  : - Fortran-90 dynamic arrays
      !!
      !! History :
      !!        !  07-03  (K. Mogensen, A. Weaver, E. Remy, S. Ricci) original
      !!----------------------------------------------------------------------
      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) :: &
         & surf      ! Surface data to be allocated

      !!* Local variables

      ! Deallocate arrays of number of surface data size

      DEALLOCATE( &
         & surf%mi,      &
         & surf%mj,      &
         & surf%mt,      &
         & surf%nsidx,   &
         & surf%nsfil,   &
         & surf%nyea,    &
         & surf%nmon,    &
         & surf%nday,    &
         & surf%nhou,    &
         & surf%nmin,    &
         & surf%mstp,    &
         & surf%nqc,     &
         & surf%ntyp,    &
         & surf%cwmo,    &
         & surf%rlam,    &
         & surf%rphi,    &
         & surf%nsind    &
         & )

      ! Allocate arrays of number of surface data size * number of variables

      DEALLOCATE( & 
         & surf%robs,    &
         & surf%rmod     &
         & )

      ! Deallocate arrays of number of extra fields at observation points

      DEALLOCATE( & 
         & surf%rext &
         & )

      ! Deallocate arrays of size number of grid points size times
      ! number of variables

      DEALLOCATE( &
         & surf%vdmean &
         & )

      ! Deallocate arrays of number of time step size

      DEALLOCATE( &
         & surf%nsstp,     &
         & surf%nsstpmpp   &
         & )

      ! Dellocate arrays of size number of variables

      DEALLOCATE( &
         & surf%cvars     &
         & )

   END SUBROUTINE obs_surf_dealloc

   SUBROUTINE obs_surf_compress( surf, newsurf, lallocate, kumout, lvalid )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE obs_surf_compress  ***
      !!                      
      !! ** Purpose : - Extract sub-information from a obs_surf type
      !!                into a new obs_surf type
      !! 
      !! ** Method  : - The data is copied from surf to new surf.
      !!                In the case of lvalid being present only the 
      !!                selected data will be copied.
      !!                If lallocate is true the data in the newsurf is 
      !!                allocated either with the same number of elements
      !!                as surf or with only the subset of elements defined
      !!                by the optional selection.
      !!
      !! History :
      !!        !  07-03  (K. Mogensen, A. Weaver, E. Remy, S. Ricci) original
      !!----------------------------------------------------------------------
      !! * Arguments
      TYPE(obs_surf), INTENT(IN)    :: surf      ! Original surface data
      TYPE(obs_surf), INTENT(INOUT) :: newsurf   ! New surface data with a subset of the original data
      LOGICAL :: lallocate     ! Allocate newsurf data
      INTEGER,INTENT(IN) :: kumout        ! Fortran unit for messages
      LOGICAL, OPTIONAL, INTENT(in), DIMENSION(:) :: &
         & lvalid         ! Valid of surface observations
      
      !!* Local variables
      INTEGER :: insurf
      INTEGER :: ji
      INTEGER :: jk
      LOGICAL, DIMENSION(:), ALLOCATABLE :: llvalid

      ! Count how many elements there should be in the new data structure

      IF ( PRESENT(lvalid) ) THEN
         insurf = 0
         DO ji = 1, surf%nsurf
            IF ( lvalid(ji) ) THEN
               insurf = insurf + 1
            ENDIF
         END DO
      ELSE
         insurf = surf%nsurf
      ENDIF

      ! Optionally allocate data in the new data structure

      IF ( lallocate ) THEN
         CALL obs_surf_alloc( newsurf,  insurf, surf%nvar, &
            & surf%nextra, surf%nstp, surf%npi, surf%npj )
      ENDIF

      ! Allocate temporary valid array to unify the code for both cases

      ALLOCATE( llvalid(surf%nsurf) )
      IF ( PRESENT(lvalid) ) THEN
         llvalid(:)  = lvalid(:)
      ELSE
         llvalid(:)  = .TRUE.
      ENDIF

      ! Setup bookkeeping variables

      insurf = 0

      ! Loop over source surface data

      DO ji = 1, surf%nsurf

         IF ( llvalid(ji) ) THEN

            ! Copy the header information

            insurf = insurf + 1

            newsurf%mi(insurf)    = surf%mi(ji)
            newsurf%mj(insurf)    = surf%mj(ji)
            newsurf%mt(insurf)    = surf%mt(ji)
            newsurf%nsidx(insurf) = surf%nsidx(ji)
            newsurf%nsfil(insurf) = surf%nsfil(ji)
            newsurf%nyea(insurf)  = surf%nyea(ji)
            newsurf%nmon(insurf)  = surf%nmon(ji)
            newsurf%nday(insurf)  = surf%nday(ji)
            newsurf%nhou(insurf)  = surf%nhou(ji)
            newsurf%nmin(insurf)  = surf%nmin(ji)
            newsurf%mstp(insurf)  = surf%mstp(ji)
            newsurf%nqc(insurf)   = surf%nqc(ji)
            newsurf%ntyp(insurf)  = surf%ntyp(ji)
            newsurf%cwmo(insurf)  = surf%cwmo(ji)
            newsurf%rlam(insurf)  = surf%rlam(ji)
            newsurf%rphi(insurf)  = surf%rphi(ji)

            DO jk = 1, surf%nvar

               newsurf%robs(insurf,jk)  = surf%robs(ji,jk)
               newsurf%rmod(insurf,jk)  = surf%rmod(ji,jk)
               
            END DO

            DO jk = 1, surf%nextra

               newsurf%rext(insurf,jk) = surf%rext(ji,jk)

            END DO
            
            ! nsind is the index of the original surface data
            
            newsurf%nsind(insurf) = ji

         ENDIF

      END DO

      ! Update MPP counters

      newsurf%nsurf = insurf
      CALL obs_mpp_sum_integer ( newsurf%nsurf, newsurf%nsurfmpp )

      ! Set book keeping variables which do not depend on number of obs.

      newsurf%nstp     = surf%nstp
      newsurf%cvars(:) = surf%cvars(:)
      
      ! Set gridded stuff
      
      newsurf%mt(insurf)    = surf%mt(ji)
 
      ! Deallocate temporary data

      DEALLOCATE( llvalid )
     
   END SUBROUTINE obs_surf_compress

   SUBROUTINE obs_surf_decompress( surf, oldsurf, ldeallocate, kumout )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE obs_surf_decompress  ***
      !!                      
      !! ** Purpose : - Copy back information to original surface data type
      !!
      !! ** Method  : - Reinsert updated information from a previous
      !!                copied/compressed surface data type into the original
      !!                surface data and optionally deallocate the surface
      !!                data input
      !! 
      !! History :
      !!        !  07-03  (K. Mogensen, A. Weaver, E. Remy, S. Ricci) original
      !!----------------------------------------------------------------------
      !! * Arguments
      TYPE(obs_surf),INTENT(INOUT) :: surf       ! Updated surface data
      TYPE(obs_surf),INTENT(INOUT) :: oldsurf    ! Original surface data
      LOGICAL :: ldeallocate ! Deallocate the updated data of insertion
      INTEGER,INTENT(in) :: kumout      ! Output unit
      
      !!* Local variables
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jk

      ! Copy data from surf to old surf

      DO ji = 1, surf%nsurf

         jj=surf%nsind(ji)

         oldsurf%mi(jj)    = surf%mi(ji)
         oldsurf%mj(jj)    = surf%mj(ji)
         oldsurf%mt(jj)    = surf%mt(ji)
         oldsurf%nsidx(jj) = surf%nsidx(ji)
         oldsurf%nsfil(jj) = surf%nsfil(ji)
         oldsurf%nyea(jj)  = surf%nyea(ji)
         oldsurf%nmon(jj)  = surf%nmon(ji)
         oldsurf%nday(jj)  = surf%nday(ji)
         oldsurf%nhou(jj)  = surf%nhou(ji)
         oldsurf%nmin(jj)  = surf%nmin(ji)
         oldsurf%mstp(jj)  = surf%mstp(ji)
         oldsurf%nqc(jj)   = surf%nqc(ji)
         oldsurf%ntyp(jj)  = surf%ntyp(ji)
         oldsurf%cwmo(jj)  = surf%cwmo(ji)
         oldsurf%rlam(jj)  = surf%rlam(ji)
         oldsurf%rphi(jj)  = surf%rphi(ji)

      END DO

      DO jk = 1, surf%nvar

         DO ji = 1, surf%nsurf
            
            jj=surf%nsind(ji)

            oldsurf%robs(jj,jk)  = surf%robs(ji,jk)
            oldsurf%rmod(jj,jk)  = surf%rmod(ji,jk)

         END DO

      END DO

      DO jk = 1, surf%nextra

         DO ji = 1, surf%nsurf
            
            jj=surf%nsind(ji)

            oldsurf%rext(jj,jk)  = surf%rext(ji,jk)

         END DO

      END DO

      ! Optionally deallocate the updated surface data

      IF ( ldeallocate ) CALL obs_surf_dealloc( surf )
      
   END SUBROUTINE obs_surf_decompress
   
END MODULE obs_surf_def

