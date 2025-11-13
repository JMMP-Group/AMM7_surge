MODULE obs_profiles_def
   !!=====================================================================
   !!                       ***  MODULE  obs_profiles_def  ***
   !! Observation diagnostics: Storage handling for T,S profiles
   !!                          arrays and additional flags etc.
   !!                          This module only defines the data type and
   !!                          operations on the data type. There is no
   !!                          actual data in the module.
   !!===================================================================== 

   !!----------------------------------------------------------------------
   !!   obs_prof            : F90 type containing the profile information
   !!   obs_prof_var        : F90 type containing the variable definition
   !!   obs_prof_valid       : F90 type containing the valid obs. definition
   !!   obs_prof_alloc      : Allocates profile arrays
   !!   obs_prof_dealloc    : Deallocates profile arrays
   !!   obs_prof_compress   : Extract sub-information from a obs_prof type
   !!                         to a new obs_prof type
   !!   obs_prof_decompress : Reinsert sub-information from a obs_prof type
   !!                         into the original obs_prof type
   !!   obs_prof_staend     : Set npvsta and npvend of a variable within an 
   !!                         obs_prof_var type
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind, ONLY : & ! Precision variables
      & wp         
   USE in_out_manager     ! I/O manager
   USE obs_mpp, ONLY : &  ! MPP tools 
      obs_mpp_sum_integers
   USE obs_fbm            ! Obs feedback format
   USE lib_mpp, ONLY : &
      & ctl_warn, ctl_stop

   IMPLICIT NONE

   !! * Routine/type accessibility
   PRIVATE

   PUBLIC &
      & obs_prof,           &
      & obs_prof_var,       &
      & obs_prof_valid,     &
      & obs_prof_alloc,     &
      & obs_prof_alloc_var, &
      & obs_prof_dealloc,   &
      & obs_prof_compress,  &
      & obs_prof_decompress,&
      & obs_prof_staend

   !! * Type definition for valid observations

   TYPE obs_prof_valid
      
      LOGICAL, POINTER, DIMENSION(:) :: luse

   END TYPE obs_prof_valid

   !! * Type definition for each variable

   TYPE obs_prof_var

      ! Arrays with size equal to the number of observations

      INTEGER, POINTER, DIMENSION(:) :: &
         & mvk,   &       !: k-th grid coord. for interpolating to profile data
         & nvpidx,&       !: Profile number
         & nvlidx,&       !: Level number in profile
         & nvqc,  &       !: Variable QC flags
         & idqc           !: Depth QC flag

      REAL(KIND=wp), POINTER, DIMENSION(:) :: &
         & vdep,  &       !: Depth coordinate of profile data
         & vobs,  &       !: Profile data
         & vmod           !: Model counterpart of the profile data vector

      REAL(KIND=wp), POINTER, DIMENSION(:,:) :: &
         & vext           !: Extra variables

      INTEGER, POINTER, DIMENSION(:) :: &
         & nvind          !: Source indices of temp. data in compressed data

      ! Arrays with size equal to idefnqcf times the number of observations
      INTEGER, POINTER, DIMENSION(:,:) :: &
         & idqcf,  &      !: Depth QC flags
         & nvqcf          !: Variable QC flags

   END TYPE obs_prof_var

   !! * Type definition for profile observation type

   TYPE obs_prof

      ! Bookkeeping

      INTEGER :: nvar     !: Number of variables
      INTEGER :: next     !: Number of extra fields
      INTEGER :: nprof    !: Total number of profiles within window.
      INTEGER :: nstp     !: Number of time steps
      INTEGER :: npi      !: Number of 3D grid points
      INTEGER :: npj
      INTEGER :: npk
      INTEGER :: nprofup  !: Observation counter used in obs_oper

      ! Bookkeeping arrays with sizes equal to number of variables

      CHARACTER(len=8), POINTER, DIMENSION(:) :: &
         & cvars          !: Variable names

      INTEGER, POINTER, DIMENSION(:) :: &
         & nvprot,   &    !: Local total number of profile T data
         & nvprotmpp      !: Global total number of profile T data
      
      ! Arrays with size equal to the number of profiles

      INTEGER, POINTER, DIMENSION(:) :: &
         & npidx,&        !: Profile number
         & npfil,&        !: Profile number in file
         & nyea, &        !: Year of profile
         & nmon, &        !: Month of profile
         & nday, &        !: Day of profile
         & nhou, &        !: Hour of profile
         & nmin, &        !: Minute of profile
         & mstp, &        !: Time step nearest to profile
         & nqc,  &        !: Profile QC
         & ntyp, &        !: Type of profile product (WMO table 1770)
         & ipqc, &        !: Position QC
         & itqc           !: Time QC

      REAL(KIND=wp), POINTER, DIMENSION(:) :: &
         & rlam, &        !: Longitude coordinate of profile data
         & rphi           !: Latitude coordinate of profile data

      CHARACTER(LEN=8), POINTER, DIMENSION(:) :: &
         & cwmo           !: Profile WMO indentifier
      
      ! Arrays with size equal to the number of profiles times
      ! number of variables

      INTEGER, POINTER, DIMENSION(:,:) :: &
         & npvsta, &      !: Start of each variable profile in full arrays
         & npvend, &      !: End of each variable profile in full arrays
         & mi,     &      !: i-th grid coord. for interpolating to profile T data
         & mj,     &      !: j-th grid coord. for interpolating to profile T data
         & ivqc           !: QC flags for all levels for a variable

      ! Arrays with size equal to idefnqcf
      ! the number of profiles times number of variables
      INTEGER, POINTER, DIMENSION(:,:) :: &
         & nqcf,  &       !: Observation QC flags
         & ipqcf, &       !: Position QC flags
         & itqcf          !: Time QC flags

      ! Arrays with size equal to idefnqcf
      ! the number of profiles times number of variables
      INTEGER, POINTER, DIMENSION(:,:,:) :: &
         & ivqcf

      ! Arrays of variables

      TYPE(obs_prof_var), POINTER, DIMENSION(:) :: var

      ! Arrays with size equal to the number of time steps in the window

      INTEGER, POINTER, DIMENSION(:) :: &
         & npstp,    &    !: Total number of profiles
         & npstpmpp       !: Total number of profiles

      ! Arrays with size equal to the number of time steps in the window times
      ! number of variables

      INTEGER, POINTER, DIMENSION(:,:) :: &
         & nvstp,    &    !: Local total num. of profile data each time step 
         & nvstpmpp       !: Global total num. of profile data each time step
      
      ! Arrays with size equal to the number of grid points times number of
      ! variables

      REAL(KIND=wp), POINTER, DIMENSION(:,:,:,:) :: &
         & vdmean        !: Daily averaged model field

      ! Arrays used to store source indices when 
      ! compressing obs_prof derived types
      
      ! Array with size nprof

      INTEGER, POINTER, DIMENSION(:) :: &
         & npind         !: Source indices of profile data in compressed data

   END TYPE obs_prof

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS
   
   SUBROUTINE obs_prof_alloc( prof,  kvar, kext, kprof,  &
      &                       ko3dt, kstp, kpi, kpj, kpk )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE obs_prof_alloc  ***
      !!                      
      !! ** Purpose : - Allocate data for profile arrays
      !! 
      !! ** Method  : - Fortran-90 dynamic arrays
      !!
      !! History :
      !!        !  07-01  (K. Mogensen) Original code
      !!        !  07-03  (K. Mogensen) Generalized profiles
      !!----------------------------------------------------------------------
      !! * Arguments
      TYPE(obs_prof), INTENT(INOUT) :: prof      ! Profile data to be allocated
      INTEGER, INTENT(IN) :: kprof  ! Number of profiles
      INTEGER, INTENT(IN) :: kvar   ! Number of variables
      INTEGER, INTENT(IN) :: kext   ! Number of extra fields within each variable
      INTEGER, INTENT(IN), DIMENSION(kvar) :: &
         & ko3dt     ! Number of observations per variables
      INTEGER, INTENT(IN) :: kstp   ! Number of time steps
      INTEGER, INTENT(IN) :: kpi    ! Number of 3D grid points
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kpk

      !!* Local variables
      INTEGER :: jvar
      INTEGER :: ji

      ! Set bookkeeping variables

      prof%nvar      = kvar
      prof%next      = kext
      prof%nprof     = kprof

      prof%nstp      = kstp
      prof%npi       = kpi
      prof%npj       = kpj
      prof%npk       = kpk

      ! Allocate arrays of size number of variables

      ALLOCATE( &
         & prof%cvars(kvar),    &
         & prof%nvprot(kvar),   &
         & prof%nvprotmpp(kvar) &
         )
         
      DO jvar = 1, kvar
         prof%cvars    (jvar) = "NotSet"
         prof%nvprot   (jvar) = ko3dt(jvar)
         prof%nvprotmpp(jvar) = 0
      END DO

      ! Allocate arrays of size number of profiles
      ! times number of variables
      
      ALLOCATE( &
         & prof%npvsta(kprof,kvar), &  
         & prof%npvend(kprof,kvar), &
         & prof%mi(kprof,kvar),     &
         & prof%mj(kprof,kvar),     &
         & prof%ivqc(kprof,kvar)    &
         )

      ! Allocate arrays of size iqcfdef times number of profiles
      ! times number of variables

      ALLOCATE( &
         & prof%ivqcf(idefnqcf,kprof,kvar) &
         & )

      ! Allocate arrays of size number of profiles

      ALLOCATE( &
         & prof%npidx(kprof),   &
         & prof%npfil(kprof),   &
         & prof%nyea(kprof),    &
         & prof%nmon(kprof),    &
         & prof%nday(kprof),    &
         & prof%nhou(kprof),    &
         & prof%nmin(kprof),    &
         & prof%mstp(kprof),    &
         & prof%nqc(kprof),     &
         & prof%ipqc(kprof),    &
         & prof%itqc(kprof),    &
         & prof%ntyp(kprof),    &
         & prof%rlam(kprof),    &
         & prof%rphi(kprof),    &
         & prof%cwmo(kprof),    &
         & prof%npind(kprof)    &
         & )

      ! Allocate arrays of size idefnqcf times number of profiles

      ALLOCATE( &
         & prof%nqcf(idefnqcf,kprof),  &
         & prof%ipqcf(idefnqcf,kprof), &
         & prof%itqcf(idefnqcf,kprof)  &
         & )

      ! Allocate obs_prof_var type
      ALLOCATE( &
         & prof%var(kvar) &
         & )

      ! For each variables allocate arrays of size number of observations

      DO jvar = 1, kvar

         IF ( ko3dt(jvar) >= 0 ) THEN
            CALL obs_prof_alloc_var( prof, jvar, kext, ko3dt(jvar) )
         ENDIF
         
      END DO

      ! Allocate arrays of size number of time step size

      ALLOCATE( &
         & prof%npstp(kstp),   &
         & prof%npstpmpp(kstp) &
         & )

      ! Allocate arrays of size number of time step size times 
      ! number of variables
      
      ALLOCATE( &
         & prof%nvstp(kstp,kvar),   &
         & prof%nvstpmpp(kstp,kvar) &
         & )

      ! Allocate arrays of size number of grid points size times
      ! number of variables

      ALLOCATE( &
         & prof%vdmean(kpi,kpj,kpk,kvar) &
         & )

      ! Set defaults for compression indices
      
      DO ji = 1, kprof
         prof%npind(ji) = ji
      END DO

      DO jvar = 1, kvar
         DO ji = 1, ko3dt(jvar)
            prof%var(jvar)%nvind(ji) = ji
         END DO
      END DO

      ! Set defaults for number of observations per time step

      prof%npstp(:)      = 0
      prof%npstpmpp(:)   = 0
      prof%nvstp(:,:)    = 0
      prof%nvstpmpp(:,:) = 0
      
      ! Set the observation counter used in obs_oper

      prof%nprofup     = 0

   END SUBROUTINE obs_prof_alloc

   SUBROUTINE obs_prof_dealloc( prof )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE obs_prof_dealloc  ***
      !!                      
      !! ** Purpose : - Deallocate data for profile arrays
      !! 
      !! ** Method  : - Fortran-90 dynamic arrays
      !!
      !! History :
      !!        !  07-01  (K. Mogensen) Original code
      !!----------------------------------------------------------------------
      !! * Arguments
      TYPE(obs_prof), INTENT(INOUT) :: &
         & prof      ! Profile data to be deallocated

      !!* Local variables
      INTEGER :: &
         & jvar

      ! Deallocate arrays of size number of profiles
      ! times number of variables
      
      DEALLOCATE( &
         & prof%npvsta, &  
         & prof%npvend  &
         )

      ! Dellocate arrays of size number of profiles size

      DEALLOCATE( &
         & prof%mi,      &
         & prof%mj,      &
         & prof%ivqc,    &
         & prof%ivqcf,   &
         & prof%npidx,   &
         & prof%npfil,   &
         & prof%nyea,    &
         & prof%nmon,    &
         & prof%nday,    &
         & prof%nhou,    &
         & prof%nmin,    &
         & prof%mstp,    &
         & prof%nqc,     &
         & prof%ipqc,    &
         & prof%itqc,    &
         & prof%nqcf,    &
         & prof%ipqcf,   &
         & prof%itqcf,   &
         & prof%ntyp,    &
         & prof%rlam,    &
         & prof%rphi,    &
         & prof%cwmo,    &
         & prof%npind    &
         & )

      ! For each variables allocate arrays of size number of observations

      DO jvar = 1, prof%nvar

         IF ( prof%nvprot(jvar) >= 0 ) THEN

            CALL obs_prof_dealloc_var( prof, jvar )

         ENDIF
         
      END DO

      ! Dellocate obs_prof_var type
      DEALLOCATE( &
         & prof%var &
         & )

      ! Deallocate arrays of size number of time step size

      DEALLOCATE( &
         & prof%npstp,   &
         & prof%npstpmpp &
         & )

      ! Deallocate arrays of size number of time step size times 
      ! number of variables
      
      DEALLOCATE( &
         & prof%nvstp,   &
         & prof%nvstpmpp &
         & )

      ! Deallocate arrays of size number of grid points size times
      ! number of variables

      DEALLOCATE( &
         & prof%vdmean &
         & )

      ! Dellocate arrays of size number of variables

      DEALLOCATE( &
         & prof%cvars,    &
         & prof%nvprot,   &
         & prof%nvprotmpp &
         )


   END SUBROUTINE obs_prof_dealloc


   SUBROUTINE obs_prof_alloc_var( prof, kvar, kext, kobs )

      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE obs_prof_alloc_var  ***
      !!                      
      !! ** Purpose : - Allocate data for variable data in profile arrays
      !! 
      !! ** Method  : - Fortran-90 dynamic arrays
      !!
      !! History :
      !!        !  07-03  (K. Mogensen) Original code
      !! * Arguments
      TYPE(obs_prof), INTENT(INOUT) :: prof   ! Profile data to be allocated
      INTEGER, INTENT(IN) :: kvar   ! Variable number
      INTEGER, INTENT(IN) :: kext   ! Number of extra fields within each variable
      INTEGER, INTENT(IN) :: kobs   ! Number of observations
      
      ALLOCATE( & 
         & prof%var(kvar)%mvk(kobs),       &
         & prof%var(kvar)%nvpidx(kobs),    &
         & prof%var(kvar)%nvlidx(kobs),    &
         & prof%var(kvar)%nvqc(kobs),      &
         & prof%var(kvar)%idqc(kobs),      &
         & prof%var(kvar)%vdep(kobs),      &
         & prof%var(kvar)%vobs(kobs),      &
         & prof%var(kvar)%vmod(kobs),      &
         & prof%var(kvar)%nvind(kobs)      &
         & )
      ALLOCATE( & 
         & prof%var(kvar)%idqcf(idefnqcf,kobs), &
         & prof%var(kvar)%nvqcf(idefnqcf,kobs)  &
         & )
      IF (kext>0) THEN
         ALLOCATE( & 
            & prof%var(kvar)%vext(kobs,kext) &
            & )
      ENDIF

   END SUBROUTINE obs_prof_alloc_var

   SUBROUTINE obs_prof_dealloc_var( prof, kvar )

      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE obs_prof_alloc_var  ***
      !!                      
      !! ** Purpose : - Allocate data for variable data in profile arrays
      !! 
      !! ** Method  : - Fortran-90 dynamic arrays
      !!
      !! History :
      !!        !  07-03  (K. Mogensen) Original code
      !! * Arguments
      TYPE(obs_prof), INTENT(INOUT) :: prof   ! Profile data to be allocated
      INTEGER, INTENT(IN) :: kvar      ! Variable number
      
      DEALLOCATE( & 
         & prof%var(kvar)%mvk,    &
         & prof%var(kvar)%nvpidx, &
         & prof%var(kvar)%nvlidx, &
         & prof%var(kvar)%nvqc,   &
         & prof%var(kvar)%idqc,   &
         & prof%var(kvar)%vdep,   &
         & prof%var(kvar)%vobs,   &
         & prof%var(kvar)%vmod,   &
         & prof%var(kvar)%nvind,  &
         & prof%var(kvar)%idqcf,  &
         & prof%var(kvar)%nvqcf   &
         & )
      IF (prof%next>0) THEN
         DEALLOCATE( & 
            & prof%var(kvar)%vext  &
            & )
      ENDIF

   END SUBROUTINE obs_prof_dealloc_var

   SUBROUTINE obs_prof_compress( prof,   newprof, lallocate, &
      &                          kumout, lvalid,   lvvalid )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE obs_prof_compress  ***
      !!                      
      !! ** Purpose : - Extract sub-information from a obs_prof type
      !!                into a new obs_prof type
      !! 
      !! ** Method  : - The data is copied from prof to new prof.
      !!                In the case of lvalid and lvvalid both being
      !!                present only the selected data will be copied.
      !!                If lallocate is true the data in the newprof is 
      !!                allocated either with the same number of elements
      !!                as prof or with only the subset of elements defined
      !!                by the optional selection in lvalid and lvvalid
      !!
      !! History :
      !!        !  07-01  (K. Mogensen) Original code
      !!----------------------------------------------------------------------
      !! * Arguments
      TYPE(obs_prof), INTENT(IN)    :: prof      ! Original profile
      TYPE(obs_prof), INTENT(INOUT) :: newprof   ! New profile with the copy of the data
      LOGICAL :: lallocate                ! Allocate newprof data
      INTEGER,INTENT(IN) :: kumout        ! Fortran unit for messages
      TYPE(obs_prof_valid), OPTIONAL, INTENT(in) :: &
         & lvalid        ! Valid profiles
      TYPE(obs_prof_valid), OPTIONAL, INTENT(in), DIMENSION(prof%nvar) :: &
         & lvvalid       ! Valid data within the profiles
      
      !!* Local variables
      INTEGER :: inprof
      INTEGER, DIMENSION(prof%nvar) :: &
         & invpro
      INTEGER :: jvar
      INTEGER :: jext
      INTEGER :: ji
      INTEGER :: jj 
      LOGICAL :: lfirst
      TYPE(obs_prof_valid) :: &
         & llvalid
      TYPE(obs_prof_valid), DIMENSION(prof%nvar) :: &
         & llvvalid
      LOGICAL :: lallpresent
      LOGICAL :: lnonepresent

      ! Check that either all or none of the masks are persent.

      lallpresent  = .FALSE.
      lnonepresent = .FALSE.
      IF ( PRESENT(lvalid) .AND. PRESENT(lvvalid) ) THEN
         lallpresent =  .TRUE.
      ELSEIF ( ( .NOT. PRESENT(lvalid)  ) .AND. &
         &     ( .NOT. PRESENT(lvvalid) ) ) THEN
         lnonepresent = .TRUE.
      ELSE
         CALL ctl_stop('Error in obs_prof_compress:', &
            &          'Either all selection variables should be set', &
            &          'or no selection variable should be set' )
      ENDIF
      
      ! Count how many elements there should be in the new data structure

      IF ( lallpresent ) THEN
         inprof = 0
         invpro(:) = 0
         DO ji = 1, prof%nprof
            IF ( lvalid%luse(ji) ) THEN
               inprof=inprof+1
               DO jvar = 1, prof%nvar
                  DO jj = prof%npvsta(ji,jvar), prof%npvend(ji,jvar)
                     IF ( lvvalid(jvar)%luse(jj) ) &
                        &           invpro(jvar) = invpro(jvar) +1
                  END DO
               END DO
            ENDIF
         END DO
      ELSE
         inprof    = prof%nprof
         invpro(:) = prof%nvprot(:)
      ENDIF

      ! Optionally allocate data in the new data structure

      IF ( lallocate ) THEN
         CALL obs_prof_alloc( newprof,   prof%nvar, &
            &                 prof%next,            &
            &                 inprof,    invpro,    &
            &                 prof%nstp, prof%npi,  &
            &                 prof%npj,  prof%npk )
      ENDIF

      ! Allocate temporary mask array to unify the code for both cases

      ALLOCATE( llvalid%luse(prof%nprof) )
      DO jvar = 1, prof%nvar
         ALLOCATE( llvvalid(jvar)%luse(prof%nvprot(jvar)) )
      END DO
      IF ( lallpresent ) THEN
         llvalid%luse(:) = lvalid%luse(:)
         DO jvar = 1, prof%nvar
            llvvalid(jvar)%luse(:) = lvvalid(jvar)%luse(:)
         END DO
      ELSE
         llvalid%luse(:) = .TRUE.
         DO jvar = 1, prof%nvar
            llvvalid(jvar)%luse(:) = .TRUE.
         END DO
      ENDIF

      ! Setup bookkeeping variables

      inprof    = 0
      invpro(:) = 0

      newprof%npvsta(:,:) =  0
      newprof%npvend(:,:) = -1
      
      ! Loop over source profiles

      DO ji = 1, prof%nprof

         IF ( llvalid%luse(ji) ) THEN

            ! Copy the header information

            inprof = inprof + 1

            newprof%mi(inprof,:)  = prof%mi(ji,:)
            newprof%mj(inprof,:) = prof%mj(ji,:)
            newprof%npidx(inprof) = prof%npidx(ji)
            newprof%npfil(inprof) = prof%npfil(ji)
            newprof%nyea(inprof)  = prof%nyea(ji)
            newprof%nmon(inprof)  = prof%nmon(ji)
            newprof%nday(inprof)  = prof%nday(ji)
            newprof%nhou(inprof)  = prof%nhou(ji)
            newprof%nmin(inprof)  = prof%nmin(ji)
            newprof%mstp(inprof)  = prof%mstp(ji)
            newprof%nqc(inprof)   = prof%nqc(ji)
            newprof%ipqc(inprof)  = prof%ipqc(ji)
            newprof%itqc(inprof)  = prof%itqc(ji)
            newprof%ivqc(inprof,:)= prof%ivqc(ji,:)
            newprof%ntyp(inprof)  = prof%ntyp(ji)
            newprof%rlam(inprof)  = prof%rlam(ji)
            newprof%rphi(inprof)  = prof%rphi(ji)
            newprof%cwmo(inprof)  = prof%cwmo(ji)
            
            ! QC info

            newprof%nqcf(:,inprof)    = prof%nqcf(:,ji)
            newprof%ipqcf(:,inprof)   = prof%ipqcf(:,ji)
            newprof%itqcf(:,inprof)   = prof%itqcf(:,ji)
            newprof%ivqcf(:,inprof,:) = prof%ivqcf(:,ji,:)
            
            ! npind is the index of the original profile
            
            newprof%npind(inprof) = ji

            ! Copy the variable information

            DO jvar = 1, prof%nvar

               lfirst = .TRUE.
               
               DO jj = prof%npvsta(ji,jvar), prof%npvend(ji,jvar)
                  
                  IF ( llvvalid(jvar)%luse(jj) ) THEN

                     invpro(jvar) = invpro(jvar) + 1
                  
                     ! Book keeping information
                                    
                     IF ( lfirst ) THEN
                        lfirst = .FALSE.
                        newprof%npvsta(inprof,jvar) = invpro(jvar)
                     ENDIF
                     newprof%npvend(inprof,jvar) = invpro(jvar)

                     ! Variable data
                     
                     newprof%var(jvar)%mvk(invpro(jvar))    = &
                        &                           prof%var(jvar)%mvk(jj)
                     newprof%var(jvar)%nvpidx(invpro(jvar)) = &
                        &                           prof%var(jvar)%nvpidx(jj)
                     newprof%var(jvar)%nvlidx(invpro(jvar)) = &
                        &                           prof%var(jvar)%nvlidx(jj)
                     newprof%var(jvar)%nvqc(invpro(jvar))   = &
                        &                           prof%var(jvar)%nvqc(jj)
                     newprof%var(jvar)%idqc(invpro(jvar))   = &
                        &                           prof%var(jvar)%idqc(jj)
                     newprof%var(jvar)%idqcf(:,invpro(jvar))= &
                        &                           prof%var(jvar)%idqcf(:,jj)
                     newprof%var(jvar)%nvqcf(:,invpro(jvar))= &
                        &                           prof%var(jvar)%nvqcf(:,jj)
                     newprof%var(jvar)%vdep(invpro(jvar))   = &
                        &                           prof%var(jvar)%vdep(jj)
                     newprof%var(jvar)%vobs(invpro(jvar))   = &
                        &                           prof%var(jvar)%vobs(jj)
                     newprof%var(jvar)%vmod(invpro(jvar))   = &
                        &                           prof%var(jvar)%vmod(jj)
                     DO jext = 1, prof%next
                        newprof%var(jvar)%vext(invpro(jvar),jext) = &
                           &                      prof%var(jvar)%vext(jj,jext)
                     END DO
                  
                     ! nvind is the index of the original variable data
                     
                     newprof%var(jvar)%nvind(invpro(jvar))  = jj
                     
                  ENDIF

               END DO

            END DO

         ENDIF

      END DO

      ! Update MPP counters

      DO jvar = 1, prof%nvar
         newprof%nvprot(jvar) = invpro(jvar)
      END DO
      CALL obs_mpp_sum_integers ( newprof%nvprot, newprof%nvprotmpp,&
         &                        prof%nvar )
      
      ! Set book keeping variables which do not depend on number of obs.

      newprof%nvar     = prof%nvar
      newprof%next     = prof%next
      newprof%nstp     = prof%nstp
      newprof%npi      = prof%npi
      newprof%npj      = prof%npj
      newprof%npk      = prof%npk
      newprof%cvars(:) = prof%cvars(:)
 
      ! Deallocate temporary data

      DO jvar = 1, prof%nvar
         DEALLOCATE( llvvalid(jvar)%luse )
      END DO
 
      DEALLOCATE( llvalid%luse )
     
   END SUBROUTINE obs_prof_compress

   SUBROUTINE obs_prof_decompress( prof, oldprof, ldeallocate, kumout )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE obs_prof_decompress  ***
      !!                      
      !! ** Purpose : - Copy back information to original profile type
      !!
      !! ** Method  : - Reinsert updated information from a previous
      !!                copied/compressed profile type into the original
      !!                profile data and optionally deallocate the prof
      !!                data input
      !! 
      !! History :
      !!        !  07-01  (K. Mogensen) Original code
      !!----------------------------------------------------------------------
      !! * Arguments
      TYPE(obs_prof),INTENT(INOUT) :: prof      ! Updated profile data
      TYPE(obs_prof),INTENT(INOUT) :: oldprof   ! Original profile data
      LOGICAL :: ldeallocate         ! Deallocate the updated data of insertion
      INTEGER,INTENT(in) :: kumout   ! Output unit
      
      !!* Local variables
      INTEGER :: jvar
      INTEGER :: jext
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jk
      INTEGER :: jl

      DO ji = 1, prof%nprof

         ! Copy header information
         
         jk = prof%npind(ji)
      
         oldprof%mi(jk,:)  = prof%mi(ji,:)
         oldprof%mj(jk,:)  = prof%mj(ji,:)
         oldprof%npidx(jk) = prof%npidx(ji)
         oldprof%npfil(jk) = prof%npfil(ji)
         oldprof%nyea(jk)  = prof%nyea(ji)
         oldprof%nmon(jk)  = prof%nmon(ji)
         oldprof%nday(jk)  = prof%nday(ji)
         oldprof%nhou(jk)  = prof%nhou(ji)
         oldprof%nmin(jk)  = prof%nmin(ji)
         oldprof%mstp(jk)  = prof%mstp(ji)
         oldprof%nqc(jk)   = prof%nqc(ji)
         oldprof%ipqc(jk)  = prof%ipqc(ji)
         oldprof%itqc(jk)  = prof%itqc(ji)
         oldprof%ivqc(jk,:)= prof%ivqc(ji,:)
         oldprof%ntyp(jk)  = prof%ntyp(ji)
         oldprof%rlam(jk)  = prof%rlam(ji)
         oldprof%rphi(jk)  = prof%rphi(ji)
         oldprof%cwmo(jk)  = prof%cwmo(ji)
         
         ! QC info

         oldprof%nqcf(:,jk)    = prof%nqcf(:,ji)
         oldprof%ipqcf(:,jk)   = prof%ipqcf(:,ji)
         oldprof%itqcf(:,jk)   = prof%itqcf(:,ji)
         oldprof%ivqcf(:,jk,:) = prof%ivqcf(:,ji,:)

         ! Copy the variable information

         DO jvar = 1, prof%nvar

            DO jj = prof%npvsta(ji,jvar), prof%npvend(ji,jvar)
               
               jl = prof%var(jvar)%nvind(jj)

               oldprof%var(jvar)%mvk(jl)    = prof%var(jvar)%mvk(jj)
               oldprof%var(jvar)%nvpidx(jl) = prof%var(jvar)%nvpidx(jj)
               oldprof%var(jvar)%nvlidx(jl) = prof%var(jvar)%nvlidx(jj)
               oldprof%var(jvar)%nvqc(jl)   = prof%var(jvar)%nvqc(jj)
               oldprof%var(jvar)%idqc(jl)   = prof%var(jvar)%idqc(jj)
               oldprof%var(jvar)%vdep(jl)   = prof%var(jvar)%vdep(jj)
               oldprof%var(jvar)%vobs(jl)   = prof%var(jvar)%vobs(jj)
               oldprof%var(jvar)%vmod(jl)   = prof%var(jvar)%vmod(jj)
               oldprof%var(jvar)%idqcf(:,jl) = prof%var(jvar)%idqcf(:,jj)
               oldprof%var(jvar)%nvqcf(:,jl) = prof%var(jvar)%nvqcf(:,jj)
               DO jext = 1, prof%next
                  oldprof%var(jvar)%vext(jl,jext) = &
                     &                        prof%var(jvar)%vext(jj,jext)
               END DO
               
            END DO

         END DO
         
      END DO

      ! Optionally deallocate the updated profile data

      IF ( ldeallocate ) CALL obs_prof_dealloc( prof )
      
   END SUBROUTINE obs_prof_decompress

   SUBROUTINE obs_prof_staend( prof, kvarno )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE obs_prof_decompress  ***
      !!                      
      !! ** Purpose : - Set npvsta and npvend of a variable within 
      !!                an obs_prof_var type
      !!
      !! ** Method  : - Find the start and stop of a profile by searching 
      !!                through the data
      !! 
      !! History :
      !!        !  07-04  (K. Mogensen) Original code
      !!----------------------------------------------------------------------
      !! * Arguments
      TYPE(obs_prof),INTENT(INOUT) :: prof     ! Profile data
      INTEGER,INTENT(IN) :: kvarno     ! Variable number

      !!* Local variables
      INTEGER :: ji
      INTEGER :: iprofno

      !-----------------------------------------------------------------------
      ! Compute start and end bookkeeping arrays
      !-----------------------------------------------------------------------

      prof%npvsta(:,kvarno) = prof%nvprot(kvarno) + 1
      prof%npvend(:,kvarno) = -1
      DO ji = 1, prof%nvprot(kvarno)
         iprofno = prof%var(kvarno)%nvpidx(ji)
         prof%npvsta(iprofno,kvarno) = &
            & MIN( ji, prof%npvsta(iprofno,kvarno) )
         prof%npvend(iprofno,kvarno) = &
            & MAX( ji, prof%npvend(iprofno,kvarno) )
      END DO

      DO ji = 1, prof%nprof
         IF ( prof%npvsta(ji,kvarno) == ( prof%nvprot(kvarno) + 1 ) ) &
            & prof%npvsta(ji,kvarno) = 0
      END DO

   END SUBROUTINE obs_prof_staend
   
END MODULE obs_profiles_def

