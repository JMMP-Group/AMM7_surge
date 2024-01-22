MODULE obs_oper
   !!======================================================================
   !!                       ***  MODULE obs_oper  ***
   !! Observation diagnostics: Observation operators for various observation
   !!                          types
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   obs_prof_opt :    Compute the model counterpart of profile data
   !!   obs_surf_opt :    Compute the model counterpart of surface data
   !!----------------------------------------------------------------------
   USE obs_inter_sup                                        ! Interpolation support
   USE obs_inter_h2d, ONLY : obs_int_h2d, obs_int_h2d_init  ! Horizontal interpolation to the obs pt
   USE obs_averg_h2d, ONLY : obs_avg_h2d, obs_avg_h2d_init, obs_max_fpsize    ! Horizontal averaging to the obs footprint
   USE obs_inter_z1d, ONLY : obs_int_z1d, obs_int_z1d_spl   ! Vertical interpolation to the obs pt
   USE obs_const    , ONLY : obfillflt                      ! Obs fill value
   USE dom_oce,       ONLY :   glamt, glamf, gphit, gphif   ! lat/lon of ocean grid-points
   USE lib_mpp,       ONLY :   ctl_warn, ctl_stop           ! Warning and stopping routines
   USE sbcdcy,        ONLY :   sbc_dcy, nday_qsr            ! For calculation of where it is night-time
   USE obs_grid,      ONLY :   obs_level_search     
   !
   USE par_kind     , ONLY :   wp   ! Precision variables
   USE in_out_manager               ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   obs_prof_opt   !: Compute the model counterpart of profile obs
   PUBLIC   obs_surf_opt   !: Compute the model counterpart of surface obs

   INTEGER, PARAMETER, PUBLIC ::   imaxavtypes = 20   !: Max number of daily avgd obs types

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE obs_prof_opt( prodatqc, kt, kpi, kpj, kpk,          &
      &                     kit000, kdaystp,                      &
      &                     pvar1, pvar2, pgdept, pgdepw,         &
      &                     pmask1, pmask2,                       &  
      &                     plam1, plam2, pphi1, pphi2,           &
      &                     k1dint, k2dint, kdailyavtypes )
      !!-----------------------------------------------------------------------
      !!                     ***  ROUTINE obs_pro_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of profiles
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : Linearly interpolate to each observation point using 
      !!              the model values at the corners of the surrounding grid box.
      !!
      !!    First, a vertical profile of horizontally interpolated model
      !!    now values is computed at the obs (lon, lat) point.
      !!    Several horizontal interpolation schemes are available:
      !!        - distance-weighted (great circle) (k2dint = 0)
      !!        - distance-weighted (small angle)  (k2dint = 1)
      !!        - bilinear (geographical grid)     (k2dint = 2)
      !!        - bilinear (quadrilateral grid)    (k2dint = 3)
      !!        - polynomial (quadrilateral grid)  (k2dint = 4)
      !!
      !!    Next, the vertical profile is interpolated to the
      !!    data depth points. Two vertical interpolation schemes are
      !!    available:
      !!        - linear       (k1dint = 0)
      !!        - Cubic spline (k1dint = 1)
      !!
      !!    For the cubic spline the 2nd derivative of the interpolating 
      !!    polynomial is computed before entering the vertical interpolation 
      !!    routine.
      !!
      !!    If the logical is switched on, the model equivalent is
      !!    a daily mean model temperature field. So, we first compute
      !!    the mean, then interpolate only at the end of the day.
      !!
      !!    Note: in situ temperature observations must be converted
      !!    to potential temperature (the model variable) prior to
      !!    assimilation. 
      !!
      !! ** Action  :
      !!
      !! History :
      !!      ! 97-11 (A. Weaver, S. Ricci, N. Daget)
      !!      ! 06-03 (G. Smith) NEMOVAR migration
      !!      ! 06-10 (A. Weaver) Cleanup
      !!      ! 07-01 (K. Mogensen) Merge of temperature and salinity
      !!      ! 07-03 (K. Mogensen) General handling of profiles
      !!      ! 15-02 (M. Martin) Combined routine for all profile types
      !!      ! 17-02 (M. Martin) Include generalised vertical coordinate changes
      !!-----------------------------------------------------------------------
      USE obs_profiles_def ! Definition of storage space for profile obs.

      IMPLICIT NONE

      TYPE(obs_prof), INTENT(inout) ::   prodatqc        ! Subset of profile data passing QC
      INTEGER       , INTENT(in   ) ::   kt              ! Time step
      INTEGER       , INTENT(in   ) ::   kpi, kpj, kpk   ! Model grid parameters
      INTEGER       , INTENT(in   ) ::   kit000          ! Number of the first time step (kit000-1 = restart time)
      INTEGER       , INTENT(in   ) ::   k1dint          ! Vertical interpolation type (see header)
      INTEGER       , INTENT(in   ) ::   k2dint          ! Horizontal interpolation type (see header)
      INTEGER       , INTENT(in   ) ::   kdaystp         ! Number of time steps per day
      REAL(KIND=wp) , INTENT(in   ), DIMENSION(kpi,kpj,kpk) ::   pvar1 , pvar2    ! Model field     1 and 2
      REAL(KIND=wp) , INTENT(in   ), DIMENSION(kpi,kpj,kpk) ::   pmask1, pmask2   ! Land-sea mask   1 and 2
      REAL(KIND=wp) , INTENT(in   ), DIMENSION(kpi,kpj)     ::   plam1 , plam2    ! Model longitude 1 and 2
      REAL(KIND=wp) , INTENT(in   ), DIMENSION(kpi,kpj)     ::   pphi1 , pphi2    ! Model latitudes 1 and 2
      REAL(KIND=wp) , INTENT(in   ), DIMENSION(kpi,kpj,kpk) ::   pgdept, pgdepw   ! depth of T and W levels 
      INTEGER, DIMENSION(imaxavtypes), OPTIONAL ::   kdailyavtypes             ! Types for daily averages

      !! * Local declarations
      INTEGER ::   ji
      INTEGER ::   jj
      INTEGER ::   jk
      INTEGER ::   jobs
      INTEGER ::   inrc
      INTEGER ::   ipro
      INTEGER ::   idayend
      INTEGER ::   ista
      INTEGER ::   iend
      INTEGER ::   iobs
      INTEGER ::   iin, ijn, ikn, ik   ! looping indices over interpolation nodes 
      INTEGER ::   inum_obs
      INTEGER, DIMENSION(imaxavtypes) :: &
         & idailyavtypes
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi1, &
         & igrdi2, &
         & igrdj1, &
         & igrdj2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iv_indic

      REAL(KIND=wp) :: zlam
      REAL(KIND=wp) :: zphi
      REAL(KIND=wp) :: zdaystp
      REAL(KIND=wp), DIMENSION(kpk) :: &
         & zobsmask1, &
         & zobsmask2, &
         & zobsk,    &
         & zobs2k
      REAL(KIND=wp), DIMENSION(2,2,1) :: &
         & zweig1, &
         & zweig2, &
         & zweig
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE :: &
         & zmask1, &
         & zmask2, &
         & zint1,  &
         & zint2,  &
         & zinm1,  &
         & zinm2,  &
         & zgdept, & 
         & zgdepw
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zglam1, &
         & zglam2, &
         & zgphi1, &
         & zgphi2
      REAL(KIND=wp), DIMENSION(1) :: zmsk_1, zmsk_2   
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: interp_corner

      LOGICAL :: ld_dailyav

      !------------------------------------------------------------------------
      ! Local initialization 
      !------------------------------------------------------------------------
      ! Record and data counters
      inrc = kt - kit000 + 2
      ipro = prodatqc%npstp(inrc)

      ! Daily average types
      ld_dailyav = .FALSE.
      IF ( PRESENT(kdailyavtypes) ) THEN
         idailyavtypes(:) = kdailyavtypes(:)
         IF ( ANY (idailyavtypes(:) /= -1) ) ld_dailyav = .TRUE.
      ELSE
         idailyavtypes(:) = -1
      ENDIF

      ! Daily means are calculated for values over timesteps:
      !  [1 <= kt <= kdaystp], [kdaystp+1 <= kt <= 2*kdaystp], ...
      idayend = MOD( kt - kit000 + 1, kdaystp )

      IF ( ld_dailyav ) THEN

         ! Initialize daily mean for first timestep of the day
         IF ( idayend == 1 .OR. kt == 0 ) THEN
            DO jk = 1, jpk
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     prodatqc%vdmean(ji,jj,jk,1) = 0.0
                     prodatqc%vdmean(ji,jj,jk,2) = 0.0
                  END DO
               END DO
            END DO
         ENDIF

         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! Increment field 1 for computing daily mean
                  prodatqc%vdmean(ji,jj,jk,1) = prodatqc%vdmean(ji,jj,jk,1) &
                     &                        + pvar1(ji,jj,jk)
                  ! Increment field 2 for computing daily mean
                  prodatqc%vdmean(ji,jj,jk,2) = prodatqc%vdmean(ji,jj,jk,2) &
                     &                        + pvar2(ji,jj,jk)
               END DO
            END DO
         END DO

         ! Compute the daily mean at the end of day
         zdaystp = 1.0 / REAL( kdaystp )
         IF ( idayend == 0 ) THEN
            IF (lwp) WRITE(numout,*) 'Calculating prodatqc%vdmean on time-step: ',kt
            CALL FLUSH(numout)
            DO jk = 1, jpk
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     prodatqc%vdmean(ji,jj,jk,1) = prodatqc%vdmean(ji,jj,jk,1) &
                        &                        * zdaystp
                     prodatqc%vdmean(ji,jj,jk,2) = prodatqc%vdmean(ji,jj,jk,2) &
                        &                        * zdaystp
                  END DO
               END DO
            END DO
         ENDIF

      ENDIF

      ! Get the data for interpolation
      ALLOCATE( &
         & igrdi1(2,2,ipro),      &
         & igrdi2(2,2,ipro),      &
         & igrdj1(2,2,ipro),      &
         & igrdj2(2,2,ipro),      &
         & zglam1(2,2,ipro),      &
         & zglam2(2,2,ipro),      &
         & zgphi1(2,2,ipro),      &
         & zgphi2(2,2,ipro),      &
         & zmask1(2,2,kpk,ipro),  &
         & zmask2(2,2,kpk,ipro),  &
         & zint1(2,2,kpk,ipro),   &
         & zint2(2,2,kpk,ipro),   &
         & zgdept(2,2,kpk,ipro),  & 
         & zgdepw(2,2,kpk,ipro)   & 
         & )

      DO jobs = prodatqc%nprofup + 1, prodatqc%nprofup + ipro
         iobs = jobs - prodatqc%nprofup
         igrdi1(1,1,iobs) = prodatqc%mi(jobs,1)-1
         igrdj1(1,1,iobs) = prodatqc%mj(jobs,1)-1
         igrdi1(1,2,iobs) = prodatqc%mi(jobs,1)-1
         igrdj1(1,2,iobs) = prodatqc%mj(jobs,1)
         igrdi1(2,1,iobs) = prodatqc%mi(jobs,1)
         igrdj1(2,1,iobs) = prodatqc%mj(jobs,1)-1
         igrdi1(2,2,iobs) = prodatqc%mi(jobs,1)
         igrdj1(2,2,iobs) = prodatqc%mj(jobs,1)
         igrdi2(1,1,iobs) = prodatqc%mi(jobs,2)-1
         igrdj2(1,1,iobs) = prodatqc%mj(jobs,2)-1
         igrdi2(1,2,iobs) = prodatqc%mi(jobs,2)-1
         igrdj2(1,2,iobs) = prodatqc%mj(jobs,2)
         igrdi2(2,1,iobs) = prodatqc%mi(jobs,2)
         igrdj2(2,1,iobs) = prodatqc%mj(jobs,2)-1
         igrdi2(2,2,iobs) = prodatqc%mi(jobs,2)
         igrdj2(2,2,iobs) = prodatqc%mj(jobs,2)
      END DO

      ! Initialise depth arrays
      zgdept(:,:,:,:) = 0.0
      zgdepw(:,:,:,:) = 0.0

      CALL obs_int_comm_2d( 2, 2, ipro, kpi, kpj, igrdi1, igrdj1, plam1, zglam1 )
      CALL obs_int_comm_2d( 2, 2, ipro, kpi, kpj, igrdi1, igrdj1, pphi1, zgphi1 )
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi1, igrdj1, pmask1, zmask1 )
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi1, igrdj1, pvar1,   zint1 )
      
      CALL obs_int_comm_2d( 2, 2, ipro, kpi, kpj, igrdi2, igrdj2, plam2, zglam2 )
      CALL obs_int_comm_2d( 2, 2, ipro, kpi, kpj, igrdi2, igrdj2, pphi2, zgphi2 )
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi2, igrdj2, pmask2, zmask2 )
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi2, igrdj2, pvar2,   zint2 )

      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi1, igrdj1, pgdept, zgdept ) 
      CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi1, igrdj1, pgdepw, zgdepw ) 

      ! At the end of the day also get interpolated means
      IF ( ld_dailyav .AND. idayend == 0 ) THEN

         ALLOCATE( &
            & zinm1(2,2,kpk,ipro),  &
            & zinm2(2,2,kpk,ipro)   &
            & )

         CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi1, igrdj1, &
            &                  prodatqc%vdmean(:,:,:,1), zinm1 )
         CALL obs_int_comm_3d( 2, 2, ipro, kpi, kpj, kpk, igrdi2, igrdj2, &
            &                  prodatqc%vdmean(:,:,:,2), zinm2 )

      ENDIF

      ! Return if no observations to process 
      ! Has to be done after comm commands to ensure processors 
      ! stay in sync 
      IF ( ipro == 0 ) RETURN 

      DO jobs = prodatqc%nprofup + 1, prodatqc%nprofup + ipro

         iobs = jobs - prodatqc%nprofup

         IF ( kt /= prodatqc%mstp(jobs) ) THEN

            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' E R R O R : Observation',              &
                  &            ' time step is not consistent with the', &
                  &            ' model time step'
               WRITE(numout,*) ' ========='
               WRITE(numout,*)
               WRITE(numout,*) ' Record  = ', jobs,                    &
                  &            ' kt      = ', kt,                      &
                  &            ' mstp    = ', prodatqc%mstp(jobs), &
                  &            ' ntyp    = ', prodatqc%ntyp(jobs)
            ENDIF
            CALL ctl_stop( 'obs_pro_opt', 'Inconsistent time' )
         ENDIF

         zlam = prodatqc%rlam(jobs)
         zphi = prodatqc%rphi(jobs)

         ! Horizontal weights 
         ! Masked values are calculated later.  
         IF ( prodatqc%npvend(jobs,1) > 0 ) THEN

            CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,     &
               &                   zglam1(:,:,iobs), zgphi1(:,:,iobs), &
               &                   zmask1(:,:,1,iobs), zweig1, zmsk_1 )

         ENDIF

         IF ( prodatqc%npvend(jobs,2) > 0 ) THEN

            CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,     &
               &                   zglam2(:,:,iobs), zgphi2(:,:,iobs), &
               &                   zmask2(:,:,1,iobs), zweig2, zmsk_2 )
 
         ENDIF

         IF ( prodatqc%npvend(jobs,1) > 0 ) THEN

            zobsk(:) = obfillflt

            IF ( ANY (idailyavtypes(:) == prodatqc%ntyp(jobs)) ) THEN

               IF ( idayend == 0 )  THEN
                  ! Daily averaged data

                  ! vertically interpolate all 4 corners 
                  ista = prodatqc%npvsta(jobs,1) 
                  iend = prodatqc%npvend(jobs,1) 
                  inum_obs = iend - ista + 1 
                  ALLOCATE(interp_corner(2,2,inum_obs),iv_indic(inum_obs)) 

                  DO iin=1,2 
                     DO ijn=1,2 

                        IF ( k1dint == 1 ) THEN 
                           CALL obs_int_z1d_spl( kpk, & 
                              &     zinm1(iin,ijn,:,iobs), & 
                              &     zobs2k, zgdept(iin,ijn,:,iobs), & 
                              &     zmask1(iin,ijn,:,iobs)) 
                        ENDIF 
       
                        CALL obs_level_search(kpk, & 
                           &    zgdept(iin,ijn,:,iobs), & 
                           &    inum_obs, prodatqc%var(1)%vdep(ista:iend), & 
                           &    iv_indic) 

                        CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs, & 
                           &    prodatqc%var(1)%vdep(ista:iend), & 
                           &    zinm1(iin,ijn,:,iobs), & 
                           &    zobs2k, interp_corner(iin,ijn,:), & 
                           &    zgdept(iin,ijn,:,iobs), & 
                           &    zmask1(iin,ijn,:,iobs)) 
       
                     ENDDO 
                  ENDDO 

               ENDIF !idayend

            ELSE   

               ! Point data 
     
               ! vertically interpolate all 4 corners 
               ista = prodatqc%npvsta(jobs,1) 
               iend = prodatqc%npvend(jobs,1) 
               inum_obs = iend - ista + 1 
               ALLOCATE(interp_corner(2,2,inum_obs), iv_indic(inum_obs)) 
               DO iin=1,2  
                  DO ijn=1,2 
                    
                     IF ( k1dint == 1 ) THEN 
                        CALL obs_int_z1d_spl( kpk, & 
                           &    zint1(iin,ijn,:,iobs),& 
                           &    zobs2k, zgdept(iin,ijn,:,iobs), & 
                           &    zmask1(iin,ijn,:,iobs)) 
  
                     ENDIF 
       
                     CALL obs_level_search(kpk, & 
                         &        zgdept(iin,ijn,:,iobs),& 
                         &        inum_obs, prodatqc%var(1)%vdep(ista:iend), & 
                         &        iv_indic) 

                     CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs,     & 
                         &          prodatqc%var(1)%vdep(ista:iend),     & 
                         &          zint1(iin,ijn,:,iobs),            & 
                         &          zobs2k,interp_corner(iin,ijn,:), & 
                         &          zgdept(iin,ijn,:,iobs),         & 
                         &          zmask1(iin,ijn,:,iobs) )      
         
                  ENDDO 
               ENDDO 
             
            ENDIF 

            !------------------------------------------------------------- 
            ! Compute the horizontal interpolation for every profile level 
            !------------------------------------------------------------- 
             
            DO ikn=1,inum_obs 
               iend=ista+ikn-1
                  
               zweig(:,:,1) = 0._wp 
   
               ! This code forces the horizontal weights to be  
               ! zero IF the observation is below the bottom of the  
               ! corners of the interpolation nodes, Or if it is in  
               ! the mask. This is important for observations near  
               ! steep bathymetry 
               DO iin=1,2 
                  DO ijn=1,2 
     
                     depth_loop1: DO ik=kpk,2,-1 
                        IF(zmask1(iin,ijn,ik-1,iobs ) > 0.9 )THEN   
                            
                           zweig(iin,ijn,1) = &  
                              & zweig1(iin,ijn,1) * & 
                              & MAX( SIGN(1._wp,(zgdepw(iin,ijn,ik,iobs) ) & 
                              &  - prodatqc%var(1)%vdep(iend)),0._wp) 
                            
                           EXIT depth_loop1 

                        ENDIF 

                     ENDDO depth_loop1 
     
                  ENDDO 
               ENDDO 
   
               CALL obs_int_h2d( 1, 1, zweig, interp_corner(:,:,ikn), & 
                  &              prodatqc%var(1)%vmod(iend:iend) ) 

                  ! Set QC flag for any observations found below the bottom
                  ! needed as the check here is more strict than that in obs_prep
               IF (sum(zweig) == 0.0_wp) prodatqc%var(1)%nvqc(iend:iend)=4
 
            ENDDO 
 
            DEALLOCATE(interp_corner,iv_indic) 
          
         ENDIF 

         ! For the second variable
         IF ( prodatqc%npvend(jobs,2) > 0 ) THEN

            zobsk(:) = obfillflt

            IF ( ANY (idailyavtypes(:) == prodatqc%ntyp(jobs)) ) THEN

               IF ( idayend == 0 )  THEN
                  ! Daily averaged data

                  ! vertically interpolate all 4 corners 
                  ista = prodatqc%npvsta(jobs,2) 
                  iend = prodatqc%npvend(jobs,2) 
                  inum_obs = iend - ista + 1 
                  ALLOCATE(interp_corner(2,2,inum_obs),iv_indic(inum_obs)) 

                  DO iin=1,2 
                     DO ijn=1,2 

                        IF ( k1dint == 1 ) THEN 
                           CALL obs_int_z1d_spl( kpk, & 
                              &     zinm2(iin,ijn,:,iobs), & 
                              &     zobs2k, zgdept(iin,ijn,:,iobs), & 
                              &     zmask2(iin,ijn,:,iobs)) 
                        ENDIF 
       
                        CALL obs_level_search(kpk, & 
                           &    zgdept(iin,ijn,:,iobs), & 
                           &    inum_obs, prodatqc%var(2)%vdep(ista:iend), & 
                           &    iv_indic) 

                        CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs, & 
                           &    prodatqc%var(2)%vdep(ista:iend), & 
                           &    zinm2(iin,ijn,:,iobs), & 
                           &    zobs2k, interp_corner(iin,ijn,:), & 
                           &    zgdept(iin,ijn,:,iobs), & 
                           &    zmask2(iin,ijn,:,iobs)) 
       
                     ENDDO 
                  ENDDO 

               ENDIF !idayend

            ELSE   

               ! Point data 
     
               ! vertically interpolate all 4 corners 
               ista = prodatqc%npvsta(jobs,2) 
               iend = prodatqc%npvend(jobs,2) 
               inum_obs = iend - ista + 1 
               ALLOCATE(interp_corner(2,2,inum_obs), iv_indic(inum_obs)) 
               DO iin=1,2  
                  DO ijn=1,2 
                    
                     IF ( k1dint == 1 ) THEN 
                        CALL obs_int_z1d_spl( kpk, & 
                           &    zint2(iin,ijn,:,iobs),& 
                           &    zobs2k, zgdept(iin,ijn,:,iobs), & 
                           &    zmask2(iin,ijn,:,iobs)) 
  
                     ENDIF 
       
                     CALL obs_level_search(kpk, & 
                         &        zgdept(iin,ijn,:,iobs),& 
                         &        inum_obs, prodatqc%var(2)%vdep(ista:iend), & 
                         &        iv_indic) 

                     CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs,     & 
                         &          prodatqc%var(2)%vdep(ista:iend),     & 
                         &          zint2(iin,ijn,:,iobs),            & 
                         &          zobs2k,interp_corner(iin,ijn,:), & 
                         &          zgdept(iin,ijn,:,iobs),         & 
                         &          zmask2(iin,ijn,:,iobs) )      
         
                  ENDDO 
               ENDDO 
             
            ENDIF 

            !------------------------------------------------------------- 
            ! Compute the horizontal interpolation for every profile level 
            !------------------------------------------------------------- 
             
            DO ikn=1,inum_obs 
               iend=ista+ikn-1
                  
               zweig(:,:,1) = 0._wp 
   
               ! This code forces the horizontal weights to be  
               ! zero IF the observation is below the bottom of the  
               ! corners of the interpolation nodes, Or if it is in  
               ! the mask. This is important for observations near  
               ! steep bathymetry 
               DO iin=1,2 
                  DO ijn=1,2 
     
                     depth_loop2: DO ik=kpk,2,-1 
                        IF(zmask2(iin,ijn,ik-1,iobs ) > 0.9 )THEN   
                            
                           zweig(iin,ijn,1) = &  
                              & zweig2(iin,ijn,1) * & 
                              & MAX( SIGN(1._wp,(zgdepw(iin,ijn,ik,iobs) ) & 
                              &  - prodatqc%var(2)%vdep(iend)),0._wp) 
                            
                           EXIT depth_loop2 

                        ENDIF 

                     ENDDO depth_loop2 
     
                  ENDDO 
               ENDDO 
   
               CALL obs_int_h2d( 1, 1, zweig, interp_corner(:,:,ikn), & 
                  &              prodatqc%var(2)%vmod(iend:iend) ) 

                  ! Set QC flag for any observations found below the bottom
                  ! needed as the check here is more strict than that in obs_prep
               IF (sum(zweig) == 0.0_wp) prodatqc%var(2)%nvqc(iend:iend)=4
 
            ENDDO 
 
            DEALLOCATE(interp_corner,iv_indic) 
          
         ENDIF 

      ENDDO

      ! Deallocate the data for interpolation
      DEALLOCATE( &
         & igrdi1, &
         & igrdi2, &
         & igrdj1, &
         & igrdj2, &
         & zglam1, &
         & zglam2, &
         & zgphi1, &
         & zgphi2, &
         & zmask1, &
         & zmask2, &
         & zint1,  &
         & zint2,  &
         & zgdept, &
         & zgdepw  &
         & )

      ! At the end of the day also get interpolated means
      IF ( ld_dailyav .AND. idayend == 0 ) THEN
         DEALLOCATE( &
            & zinm1,  &
            & zinm2   &
            & )
      ENDIF

      prodatqc%nprofup = prodatqc%nprofup + ipro 

   END SUBROUTINE obs_prof_opt

   SUBROUTINE obs_surf_opt( surfdataqc, kt, kpi, kpj,            &
      &                     kit000, kdaystp, psurf, psurfmask,   &
      &                     k2dint, ldnightav, plamscl, pphiscl, &
      &                     lindegrees )

      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_surf_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of surface
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : Linearly interpolate to each observation point using 
      !!              the model values at the corners of the surrounding grid box.
      !!
      !!    The new model value is first computed at the obs (lon, lat) point.
      !!
      !!    Several horizontal interpolation schemes are available:
      !!        - distance-weighted (great circle) (k2dint = 0)
      !!        - distance-weighted (small angle)  (k2dint = 1)
      !!        - bilinear (geographical grid)     (k2dint = 2)
      !!        - bilinear (quadrilateral grid)    (k2dint = 3)
      !!        - polynomial (quadrilateral grid)  (k2dint = 4)
      !!
      !!    Two horizontal averaging schemes are also available:
      !!        - weighted radial footprint        (k2dint = 5)
      !!        - weighted rectangular footprint   (k2dint = 6)
      !!
      !!
      !! ** Action  :
      !!
      !! History :
      !!      ! 07-03 (A. Weaver)
      !!      ! 15-02 (M. Martin) Combined routine for surface types
      !!      ! 17-03 (M. Martin) Added horizontal averaging options
      !!-----------------------------------------------------------------------
      USE obs_surf_def  ! Definition of storage space for surface observations

      IMPLICIT NONE

      TYPE(obs_surf), INTENT(INOUT) :: &
         & surfdataqc                  ! Subset of surface data passing QC
      INTEGER, INTENT(IN) :: kt        ! Time step
      INTEGER, INTENT(IN) :: kpi       ! Model grid parameters
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kit000    ! Number of the first time step 
                                       !   (kit000-1 = restart time)
      INTEGER, INTENT(IN) :: kdaystp   ! Number of time steps per day
      INTEGER, INTENT(IN) :: k2dint    ! Horizontal interpolation type (see header)
      REAL(wp), INTENT(IN), DIMENSION(kpi,kpj) :: &
         & psurf,  &                   ! Model surface field
         & psurfmask                   ! Land-sea mask
      LOGICAL, INTENT(IN) :: ldnightav ! Logical for averaging night-time data
      REAL(KIND=wp), INTENT(IN) :: &
         & plamscl, &                  ! Diameter in metres of obs footprint in E/W, N/S directions
         & pphiscl                     ! This is the full width (rather than half-width)
      LOGICAL, INTENT(IN) :: &
         & lindegrees                  ! T=> plamscl and pphiscl are specified in degrees, F=> in metres

      !! * Local declarations
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jobs
      INTEGER :: inrc
      INTEGER :: isurf
      INTEGER :: iobs
      INTEGER :: imaxifp, imaxjfp
      INTEGER :: imodi, imodj
      INTEGER :: idayend
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi,   &
         & igrdj,   &
         & igrdip1, &
         & igrdjp1
      INTEGER, DIMENSION(:,:), SAVE, ALLOCATABLE :: &
         & icount_night,      &
         & imask_night
      REAL(wp) :: zlam
      REAL(wp) :: zphi
      REAL(wp), DIMENSION(1) :: zext, zobsmask
      REAL(wp) :: zdaystp
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zweig,  &
         & zmask,  &
         & zsurf,  &
         & zsurfm, &
         & zsurftmp, &
         & zglam,  &
         & zgphi,  &
         & zglamf, &
         & zgphif

      REAL(wp), DIMENSION(:,:), SAVE, ALLOCATABLE :: &
         & zintmp,  &
         & zouttmp, &
         & zmeanday    ! to compute model sst in region of 24h daylight (pole)

      !------------------------------------------------------------------------
      ! Local initialization 
      !------------------------------------------------------------------------
      ! Record and data counters
      inrc = kt - kit000 + 2
      isurf = surfdataqc%nsstp(inrc)

      ! Work out the maximum footprint size for the 
      ! interpolation/averaging in model grid-points - has to be even.

      CALL obs_max_fpsize( k2dint, plamscl, pphiscl, lindegrees, psurfmask, imaxifp, imaxjfp )


      IF ( ldnightav ) THEN

      ! Initialize array for night mean
         IF ( kt == 0 ) THEN
            ALLOCATE ( icount_night(kpi,kpj) )
            ALLOCATE ( imask_night(kpi,kpj) )
            ALLOCATE ( zintmp(kpi,kpj) )
            ALLOCATE ( zouttmp(kpi,kpj) )
            ALLOCATE ( zmeanday(kpi,kpj) )
            nday_qsr = -1   ! initialisation flag for nbc_dcy
         ENDIF

         ! Night-time means are calculated for night-time values over timesteps:
         !  [1 <= kt <= kdaystp], [kdaystp+1 <= kt <= 2*kdaystp], .....
         idayend = MOD( kt - kit000 + 1, kdaystp )

         ! Initialize night-time mean for first timestep of the day
         IF ( idayend == 1 .OR. kt == 0 ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  surfdataqc%vdmean(ji,jj) = 0.0
                  zmeanday(ji,jj) = 0.0
                  icount_night(ji,jj) = 0
               END DO
            END DO
         ENDIF

         zintmp(:,:) = 0.0
         zouttmp(:,:) = sbc_dcy( zintmp(:,:), .TRUE. )
         imask_night(:,:) = INT( zouttmp(:,:) )

         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Increment the temperature field for computing night mean and counter
               surfdataqc%vdmean(ji,jj) = surfdataqc%vdmean(ji,jj)  &
                      &                    + psurf(ji,jj) * REAL( imask_night(ji,jj) )
               zmeanday(ji,jj)          = zmeanday(ji,jj) + psurf(ji,jj)
               icount_night(ji,jj)      = icount_night(ji,jj) + imask_night(ji,jj)
            END DO
         END DO

         ! Compute the night-time mean at the end of the day
         zdaystp = 1.0 / REAL( kdaystp )
         IF ( idayend == 0 ) THEN
            IF (lwp) WRITE(numout,*) 'Calculating surfdataqc%vdmean on time-step: ',kt
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! Test if "no night" point
                  IF ( icount_night(ji,jj) > 0 ) THEN
                     surfdataqc%vdmean(ji,jj) = surfdataqc%vdmean(ji,jj) &
                       &                        / REAL( icount_night(ji,jj) )
                  ELSE
                     !At locations where there is no night (e.g. poles),
                     ! calculate daily mean instead of night-time mean.
                     surfdataqc%vdmean(ji,jj) = zmeanday(ji,jj) * zdaystp
                  ENDIF
               END DO
            END DO
         ENDIF

      ENDIF

      ! Get the data for interpolation

      ALLOCATE( &
         & zweig(imaxifp,imaxjfp,1),      &
         & igrdi(imaxifp,imaxjfp,isurf), &
         & igrdj(imaxifp,imaxjfp,isurf), &
         & zglam(imaxifp,imaxjfp,isurf), &
         & zgphi(imaxifp,imaxjfp,isurf), &
         & zmask(imaxifp,imaxjfp,isurf), &
         & zsurf(imaxifp,imaxjfp,isurf), &
         & zsurftmp(imaxifp,imaxjfp,isurf),  &
         & zglamf(imaxifp+1,imaxjfp+1,isurf), &
         & zgphif(imaxifp+1,imaxjfp+1,isurf), &
         & igrdip1(imaxifp+1,imaxjfp+1,isurf), &
         & igrdjp1(imaxifp+1,imaxjfp+1,isurf) &
         & )

      DO jobs = surfdataqc%nsurfup + 1, surfdataqc%nsurfup + isurf
         iobs = jobs - surfdataqc%nsurfup
         DO ji = 0, imaxifp
            imodi = surfdataqc%mi(jobs) - int(imaxifp/2) + ji - 1
            !
            !Deal with wrap around in longitude
            IF ( imodi < 1      ) imodi = imodi + jpiglo
            IF ( imodi > jpiglo ) imodi = imodi - jpiglo
            !
            DO jj = 0, imaxjfp
               imodj = surfdataqc%mj(jobs) - int(imaxjfp/2) + jj - 1
               !If model values are out of the domain to the north/south then
               !set them to be the edge of the domain
               IF ( imodj < 1      ) imodj = 1
               IF ( imodj > jpjglo ) imodj = jpjglo
               !
               igrdip1(ji+1,jj+1,iobs) = imodi
               igrdjp1(ji+1,jj+1,iobs) = imodj
               !
               IF ( ji >= 1 .AND. jj >= 1 ) THEN
                  igrdi(ji,jj,iobs) = imodi
                  igrdj(ji,jj,iobs) = imodj
               ENDIF
               !
            END DO
         END DO
      END DO

      CALL obs_int_comm_2d( imaxifp, imaxjfp, isurf, kpi, kpj, &
         &                  igrdi, igrdj, glamt, zglam )
      CALL obs_int_comm_2d( imaxifp, imaxjfp, isurf, kpi, kpj, &
         &                  igrdi, igrdj, gphit, zgphi )
      CALL obs_int_comm_2d( imaxifp, imaxjfp, isurf, kpi, kpj, &
         &                  igrdi, igrdj, psurfmask, zmask )
      CALL obs_int_comm_2d( imaxifp, imaxjfp, isurf, kpi, kpj, &
         &                  igrdi, igrdj, psurf, zsurf )
      CALL obs_int_comm_2d( imaxifp+1, imaxjfp+1, isurf, kpi, kpj, &
         &                  igrdip1, igrdjp1, glamf, zglamf )
      CALL obs_int_comm_2d( imaxifp+1, imaxjfp+1, isurf, kpi, kpj, &
         &                  igrdip1, igrdjp1, gphif, zgphif )

      ! At the end of the day get interpolated means
      IF ( idayend == 0 .AND. ldnightav ) THEN

         ALLOCATE( &
            & zsurfm(imaxifp,imaxjfp,isurf)  &
            & )

         CALL obs_int_comm_2d( imaxifp,imaxjfp, isurf, kpi, kpj, igrdi, igrdj, &
            &               surfdataqc%vdmean(:,:), zsurfm )

      ENDIF

      ! Loop over observations
      DO jobs = surfdataqc%nsurfup + 1, surfdataqc%nsurfup + isurf

         iobs = jobs - surfdataqc%nsurfup

         IF ( kt /= surfdataqc%mstp(jobs) ) THEN

            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' E R R O R : Observation',              &
                  &            ' time step is not consistent with the', &
                  &            ' model time step'
               WRITE(numout,*) ' ========='
               WRITE(numout,*)
               WRITE(numout,*) ' Record  = ', jobs,                &
                  &            ' kt      = ', kt,                  &
                  &            ' mstp    = ', surfdataqc%mstp(jobs), &
                  &            ' ntyp    = ', surfdataqc%ntyp(jobs)
            ENDIF
            CALL ctl_stop( 'obs_surf_opt', 'Inconsistent time' )

         ENDIF

         zlam = surfdataqc%rlam(jobs)
         zphi = surfdataqc%rphi(jobs)

         IF ( ldnightav .AND. idayend == 0 ) THEN
            ! Night-time averaged data
            zsurftmp(:,:,iobs) = zsurfm(:,:,iobs)
         ELSE
            zsurftmp(:,:,iobs) = zsurf(:,:,iobs)
         ENDIF

         IF ( k2dint <= 4 ) THEN

            ! Get weights to interpolate the model value to the observation point
            CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,         &
               &                   zglam(:,:,iobs), zgphi(:,:,iobs), &
               &                   zmask(:,:,iobs), zweig, zobsmask )

            ! Interpolate the model value to the observation point 
            CALL obs_int_h2d( 1, 1, zweig, zsurftmp(:,:,iobs), zext )

         ELSE

            ! Get weights to average the model SLA to the observation footprint
            CALL obs_avg_h2d_init( 1, 1, imaxifp, imaxjfp, k2dint, zlam,  zphi, &
               &                   zglam(:,:,iobs), zgphi(:,:,iobs), &
               &                   zglamf(:,:,iobs), zgphif(:,:,iobs), &
               &                   zmask(:,:,iobs), plamscl, pphiscl, &
               &                   lindegrees, zweig, zobsmask )

            ! Average the model SST to the observation footprint
            CALL obs_avg_h2d( 1, 1, imaxifp, imaxjfp, &
               &              zweig, zsurftmp(:,:,iobs),  zext )

         ENDIF

         IF ( TRIM(surfdataqc%cvars(1)) == 'SLA' .AND. surfdataqc%nextra == 2 ) THEN
            ! ... Remove the MDT from the SSH at the observation point to get the SLA
            surfdataqc%rext(jobs,1) = zext(1)
            surfdataqc%rmod(jobs,1) = surfdataqc%rext(jobs,1) - surfdataqc%rext(jobs,2)
         ELSE
            surfdataqc%rmod(jobs,1) = zext(1)
         ENDIF
         
         IF ( zext(1) == obfillflt ) THEN
            ! If the observation value is a fill value, set QC flag to bad
            surfdataqc%nqc(jobs) = 4
         ENDIF

      END DO

      ! Deallocate the data for interpolation
      DEALLOCATE( &
         & zweig, &
         & igrdi, &
         & igrdj, &
         & zglam, &
         & zgphi, &
         & zmask, &
         & zsurf, &
         & zsurftmp, &
         & zglamf, &
         & zgphif, &
         & igrdip1,&
         & igrdjp1 &
         & )

      ! At the end of the day also deallocate night-time mean array
      IF ( idayend == 0 .AND. ldnightav ) THEN
         DEALLOCATE( &
            & zsurfm  &
            & )
      ENDIF
      !
      surfdataqc%nsurfup = surfdataqc%nsurfup + isurf
      !
   END SUBROUTINE obs_surf_opt

   !!======================================================================
END MODULE obs_oper
