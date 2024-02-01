MODULE obs_grid
   !!======================================================================
   !!                       ***  MODULE  obs_grid  ***
   !! Observation diagnostics: Various tools for grid searching etc.
   !!======================================================================
   !!----------------------------------------------------------------------
   !!   obs_grid_search   : Find i,j on the ORCA grid from lat,lon
   !!   obs_level_search  : Find level from depth
   !!   obs_zlevel_search : Find depth level from observed depth
   !!   obs_tlevel_search : Find temperature level from observed temp
   !!   obs_rlevel_search : Find density level from observed density
   !!----------------------------------------------------------------------
   !! * Modules used   
   USE par_kind, ONLY : &          ! Precision variables
      & wp 
   USE par_oce, ONLY :  &          ! Ocean parameters
      & jpk,     &
      & jpni,    &
      & jpnj,    &
      & jpnij
   USE dom_oce                     ! Ocean space and time domain variables
   USE obs_mpp, ONLY : &           ! MPP support routines for observation diagnostics
      & obs_mpp_find_obs_proc, &
      & mpp_global_max,        &
      & obs_mpp_max_integer
   USE phycst, ONLY : &            ! Physical constants
      & rad
   USE obs_utils, ONLY : &         ! Observation operator utility functions
      & grt_cir_dis, &
      & chkerr
   USE in_out_manager              ! Printing support
   USE netcdf
   USE obs_const, ONLY :      &
      & obfillflt                  ! Fillvalue   
   USE lib_mpp, ONLY :   &
      & ctl_warn, ctl_stop

   IMPLICIT NONE

   !! * Routine accessibility
   PUBLIC  obs_grid_setup,      & ! Setup grid searching
      &    obs_grid_search,     & ! Find i, j on the ORCA grid from lat, lon
      &    obs_grid_deallocate, & ! Deallocate the look up table
      &    obs_level_search       ! Find level from depth

   PRIVATE linquad,                    & ! Determine whether a point lies within a cell
      &    maxdist,                    & ! Find the maximum distance between 2 pts in a cell
      &    obs_grd_bruteforce, & ! Find i, j on the ORCA grid from lat, lon
      &    obs_grd_lookup        ! Find i, j on the ORCA grid from lat, lon quicker

   !!* Module variables

   !! Default values
   REAL, PUBLIC :: rn_gridsearchres = 0.5   ! Resolution of grid
   INTEGER, PRIVATE :: gsearch_nlons_def    ! Num of longitudes
   INTEGER, PRIVATE :: gsearch_nlats_def    ! Num of latitudes
   REAL(wp), PRIVATE :: gsearch_lonmin_def  ! Min longitude
   REAL(wp), PRIVATE :: gsearch_latmin_def  ! Min latitude
   REAL(wp), PRIVATE :: gsearch_dlon_def    ! Lon spacing
   REAL(wp), PRIVATE :: gsearch_dlat_def    ! Lat spacing
   !! Variable versions
   INTEGER, PRIVATE :: nlons     ! Num of longitudes
   INTEGER, PRIVATE :: nlats     ! Num of latitudes
   REAL(wp), PRIVATE :: lonmin ! Min longitude
   REAL(wp), PRIVATE :: latmin ! Min latitude
   REAL(wp), PRIVATE :: dlon     ! Lon spacing
   REAL(wp), PRIVATE :: dlat     ! Lat spacing
   
   INTEGER, PRIVATE :: maxxdiff, maxydiff ! Max diffs between model points
   INTEGER, PRIVATE :: limxdiff, limydiff
   
   ! Data storage
   REAL(wp), PRIVATE, DIMENSION(:,:), ALLOCATABLE :: &
      & lons,  &
      & lats
   INTEGER, PRIVATE, DIMENSION(:,:), ALLOCATABLE :: &
      & ixpos, &
      & iypos, &
      & iprocn    

   ! Switches
   LOGICAL, PUBLIC :: ln_grid_search_lookup  ! Use lookup table to speed up grid search
   LOGICAL, PUBLIC :: ln_grid_global         ! Use global distribution of observations
   CHARACTER(LEN=44), PUBLIC :: &
      & cn_gridsearchfile    ! file name head for grid search lookup 

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE obs_grid_search( kobsin, plam, pphi, kobsi, kobsj, kproc, &
      &                        cdgrid )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE obs_grid_search ***
      !!
      !! ** Purpose : Search local gridpoints to find the grid box containing
      !!              the observations calls either
      !!              obs_grd_bruteforce - the original brute force search
      !!                     or
      !!              obs_grd_lookup - uses a lookup table to do a fast 
      !!search
      !!History :
      !!        !  2007-12  (D. Lea) 
      !!------------------------------------------------------------------------

      !! * Arguments
      INTEGER :: &
         & kobsin                     ! Size of the observation arrays
      REAL(KIND=wp), DIMENSION(kobsin), INTENT(IN) :: &
         & plam, &                  ! Longitude of obsrvations 
         & pphi                     ! Latitude of observations
      INTEGER, DIMENSION(kobsin), INTENT(OUT) :: &
         & kobsi, &                 ! I-index of observations 
         & kobsj, &                 ! J-index of observations 
         & kproc                    ! Processor number of observations
      CHARACTER(LEN=1) :: &
         & cdgrid                   ! Grid to search

      IF(kobsin > 0) THEN

         IF ( ln_grid_search_lookup .AND. ( cdgrid == 'T' ) ) THEN
            CALL obs_grd_lookup( kobsin, plam, pphi, &
               &                         kobsi, kobsj, kproc )
         ELSE
            IF ( cdgrid == 'T' ) THEN
               CALL obs_grd_bruteforce( jpi, jpj, jpiglo, jpjglo, &
                  &                             1, nlci, 1, nlcj,         &
                  &                             nproc, jpnij,             &
                  &                             glamt, gphit, tmask,      &
                  &                             kobsin, plam, pphi,       &
                  &                             kobsi, kobsj, kproc )
            ELSEIF ( cdgrid == 'U' ) THEN
               CALL obs_grd_bruteforce( jpi, jpj, jpiglo, jpjglo, &
                  &                             1, nlci, 1, nlcj,         &
                  &                             nproc, jpnij,             &
                  &                             glamu, gphiu, umask,      &
                  &                             kobsin, plam, pphi,       &
                  &                             kobsi, kobsj, kproc )
            ELSEIF ( cdgrid == 'V' ) THEN
               CALL obs_grd_bruteforce( jpi, jpj, jpiglo, jpjglo, &
                  &                             1, nlci, 1, nlcj,         &
                  &                             nproc, jpnij,             &
                  &                             glamv, gphiv, vmask,      &
                  &                             kobsin, plam, pphi,       &
                  &                             kobsi, kobsj, kproc )
            ELSEIF ( cdgrid == 'F' ) THEN
               CALL obs_grd_bruteforce( jpi, jpj, jpiglo, jpjglo, &
                  &                             1, nlci, 1, nlcj,         &
                  &                             nproc, jpnij,             &
                  &                             glamf, gphif, fmask,      &
                  &                             kobsin, plam, pphi,       &
                  &                             kobsi, kobsj, kproc )
            ELSE
               CALL ctl_stop( 'Grid not supported' )
            ENDIF
         ENDIF
         
      ENDIF

   END SUBROUTINE obs_grid_search

#include "obs_grd_bruteforce.h90"
   
   SUBROUTINE obs_grd_lookup( kobs, plam, pphi, kobsi, kobsj, kproc )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE obs_grid_lookup ***
      !!
      !! ** Purpose : Search local gridpoints to find the grid box containing
      !!              the observations (much faster then obs_grd_bruteforce)
      !!
      !! ** Method  : Call to linquad
      !!
      !! ** Action  : Return kproc holding the observation and kiobsi,kobsj
      !!              valid on kproc=nproc processor only.
      !!   
      !! History :
      !!        !  2007-12 (D. Lea) new routine based on obs_grid_search
      !!! updated with fixes from new version of obs_grid_search_bruteforce
      !!! speeded up where points are not near a "difficult" region like an edge
      !!----------------------------------------------------------------------

      !! * Arguments
      INTEGER :: kobs                     ! Size of the observation arrays
      REAL(KIND=wp), DIMENSION(kobs), INTENT(IN) :: &
         & plam, &                  ! Longitude of obsrvations 
         & pphi                     ! Latitude of observations
      INTEGER, DIMENSION(kobs), INTENT(OUT) :: &
         & kobsi, &                 ! I-index of observations 
         & kobsj, &                 ! J-index of observations 
         & kproc                    ! Processor number of observations
  
      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: &
         & zplam
      REAL(wp) :: zlammax
      REAL(wp) :: zlam
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jk
      INTEGER :: jo
      INTEGER :: isx
      INTEGER :: isy
      INTEGER :: jimin
      INTEGER :: jimax
      INTEGER :: jjmin
      INTEGER :: jjmax
      INTEGER :: jojimin
      INTEGER :: jojimax
      INTEGER :: jojjmin
      INTEGER :: jojjmax
      INTEGER :: ipx1
      INTEGER :: ipy1
      INTEGER :: ip
      INTEGER :: jp
      INTEGER :: ipx
      INTEGER :: ipy
      INTEGER :: ipmx
      INTEGER :: jlon
      INTEGER :: jlat
      INTEGER :: joffset
      INTEGER :: jostride
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zlamg, &
         & zphig, &
         & zmskg, &
         & zphitmax,&
         & zphitmin,&
         & zlamtmax,&
         & zlamtmin
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: &
         & llinvalidcell
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zlamtm,  &
         & zphitm
      LOGICAL :: llfourflag
      INTEGER :: ifourflagcountt
      INTEGER :: ifourflagcountf
      INTEGER, DIMENSION(5) :: ifourflagcountr

      !-----------------------------------------------------------------------
      ! Define grid for grid search
      !-----------------------------------------------------------------------
      IF (ln_grid_global) THEN
         jlon     = jpiglo
         jlat     = jpjglo
         joffset  = nproc
         jostride = jpnij
      ELSE
         jlon     = jpi
         jlat     = jpj
         joffset  = 0
         jostride = 1
      ENDIF
      !-----------------------------------------------------------------------
      ! Set up data for grid search
      !-----------------------------------------------------------------------
      ALLOCATE( &
         & zlamg(jlon,jlat),             &
         & zphig(jlon,jlat),             &
         & zmskg(jlon,jlat),             &
         & zphitmax(jlon-1,jlat-1),      &
         & zphitmin(jlon-1,jlat-1),      &
         & zlamtmax(jlon-1,jlat-1),      &
         & zlamtmin(jlon-1,jlat-1),      &
         & llinvalidcell(jlon-1,jlat-1), &
         & zlamtm(4,jlon-1,jlat-1),      &
         & zphitm(4,jlon-1,jlat-1)       &
         & )
      !-----------------------------------------------------------------------
      ! Copy data to local arrays
      !-----------------------------------------------------------------------
      IF (ln_grid_global) THEN
         zlamg(:,:) = -1.e+10
         zphig(:,:) = -1.e+10
         zmskg(:,:) = -1.e+10
         ! Add various grids here.
         DO jj = 1, nlcj
            DO ji = 1, nlci
               zlamg(mig(ji),mjg(jj)) = glamt(ji,jj)
               zphig(mig(ji),mjg(jj)) = gphit(ji,jj)
               zmskg(mig(ji),mjg(jj)) = tmask(ji,jj,1)
            END DO
         END DO
         CALL mpp_global_max( zlamg )
         CALL mpp_global_max( zphig )
         CALL mpp_global_max( zmskg )
      ELSE
         ! Add various grids here.
         DO jj = 1, jlat
            DO ji = 1, jlon
               zlamg(ji,jj) = glamt(ji,jj)
               zphig(ji,jj) = gphit(ji,jj)
               zmskg(ji,jj) = tmask(ji,jj,1)
            END DO
         END DO
      ENDIF
      !-----------------------------------------------------------------------
      ! Copy longitudes
      !----------------------------------------------------------------------- 
      ALLOCATE( &
         & zplam(kobs) &
         & )
      DO jo = 1, kobs
         zplam(jo) = plam(jo)
      END DO
      !-----------------------------------------------------------------------
      ! Set default values for output
      !-----------------------------------------------------------------------
      kproc(:) = -1
      kobsi(:) = -1
      kobsj(:) = -1
      !-----------------------------------------------------------------------
      ! Copy grid positions to temporary arrays and renormalize to 0 to 360.
      !-----------------------------------------------------------------------
      DO jj = 1, jlat-1
         DO ji = 1, jlon-1
            zlamtm(1,ji,jj) = zlamg(ji  ,jj  )
            zphitm(1,ji,jj) = zphig(ji  ,jj  )
            zlamtm(2,ji,jj) = zlamg(ji+1,jj  )
            zphitm(2,ji,jj) = zphig(ji+1,jj  )
            zlamtm(3,ji,jj) = zlamg(ji+1,jj+1)
            zphitm(3,ji,jj) = zphig(ji+1,jj+1)
            zlamtm(4,ji,jj) = zlamg(ji  ,jj+1)
            zphitm(4,ji,jj) = zphig(ji  ,jj+1)
         END DO
      END DO
      WHERE ( zlamtm(:,:,:) < 0.0_wp )
         zlamtm(:,:,:) = zlamtm(:,:,:) + 360.0_wp
      END WHERE
      WHERE ( zlamtm(:,:,:) > 360.0_wp )
         zlamtm(:,:,:) = zlamtm(:,:,:) - 360.0_wp
      END WHERE
      !-----------------------------------------------------------------------
      ! Handle case of the wraparound; beware, not working with orca180
      !-----------------------------------------------------------------------
      DO jj = 1, jlat-1
         DO ji = 1, jlon-1
            zlammax = MAXVAL( zlamtm(:,ji,jj) )
            WHERE (zlammax - zlamtm(:, ji, jj) > 180 ) & 
               & zlamtm(:,ji,jj) = zlamtm(:,ji,jj) + 360._wp
            zphitmax(ji,jj) = MAXVAL(zphitm(:,ji,jj))
            zphitmin(ji,jj) = MINVAL(zphitm(:,ji,jj))
            zlamtmax(ji,jj) = MAXVAL(zlamtm(:,ji,jj))
            zlamtmin(ji,jj) = MINVAL(zlamtm(:,ji,jj))
         END DO
      END DO
      !-----------------------------------------------------------------------
      ! Search for boxes with only land points mark them invalid
      !-----------------------------------------------------------------------
      llinvalidcell(:,:) = .FALSE.
      DO jj = 1, jlat-1
         DO ji = 1, jlon-1
            llinvalidcell(ji,jj) =               &
               & zmskg(ji  ,jj  ) == 0.0_wp .AND. &
               & zmskg(ji+1,jj  ) == 0.0_wp .AND. &
               & zmskg(ji+1,jj+1) == 0.0_wp .AND. &
               & zmskg(ji  ,jj+1) == 0.0_wp
         END DO
      END DO

      if(lwp) WRITE(numout,*) 'obs_grid_lookup do coordinate search using lookup table'

      !-----------------------------------------------------------------------
      ! Do coordinate search using lookup table with local searches.
      ! - For land points kproc is set to number of the processor + 1000000
      !   and we continue the search.
      ! - For ocean points kproc is set to the number of the processor 
      !   and we stop the search.
      !-----------------------------------------------------------------------
      ifourflagcountt = 0
      ifourflagcountf = 0
      ifourflagcountr(:) = 0

      !------------------------------------------------------------------------
      ! Master loop for grid search
      !------------------------------------------------------------------------
         
      gpkobs: DO jo = 1+joffset, kobs, jostride
         ! Normal case
         !        specify 4 points which surround the lat lon of interest
         !                          x i,j+1  x i+1, j+1
         !
         !
         !                             * lon,lat
         !                          x i,j    x i+1,j

         ! bottom corner point
         ipx1 = INT( ( zplam(jo)  - lonmin ) / dlon + 1.0 ) 
         ipy1 = INT( ( pphi (jo)  - latmin ) / dlat + 1.0 )  
         
         ipx = ipx1 + 1
         ipy = ipy1 + 1

         ! flag for searching around four points separately
         ! default to false
         llfourflag = .FALSE.
         
         ! check for point fully outside of region
         IF ( (ipx1 > nlons) .OR. (ipy1 > nlats) .OR. &
            & (ipx < 1) .OR. (ipy < 1) ) THEN
            CYCLE
         ENDIF
         ! check wrap around
         IF ( (ipx > nlons) .OR. (ipy > nlats) .OR. &
            & (ipx1 < 1) .OR. (ipy1 < 1) ) THEN
            llfourflag=.TRUE.
            ifourflagcountr(1) = ifourflagcountr(1) + 1
         ENDIF

         IF (.NOT. llfourflag) THEN
            IF (MAXVAL(ixpos(ipx1:ipx,ipy1:ipy)) == -1) CYCLE! cycle if no lookup points found
         ENDIF
         
         jimin = 0
         jimax = 0
         jjmin = 0
         jjmax = 0
         
         IF (.NOT. llfourflag) THEN 

            ! calculate points range
            ! define a square region encompassing the four corner points
            ! do I need the -1 points?

            jojimin = MINVAL(ixpos(ipx1:ipx,ipy1:ipy)) - 1
            jojimax = MAXVAL(ixpos(ipx1:ipx,ipy1:ipy)) + 1
            jojjmin = MINVAL(iypos(ipx1:ipx,ipy1:ipy)) - 1
            jojjmax = MAXVAL(iypos(ipx1:ipx,ipy1:ipy)) + 1

            jimin = jojimin - 1
            jimax = jojimax + 1
            jjmin = jojjmin - 1
            jjmax = jojjmax + 1
            
            IF ( jojimin < 0 .OR. jojjmin < 0) THEN 
               llfourflag = .TRUE.
               ifourflagcountr(2) = ifourflagcountr(2) + 1
            ENDIF
            IF ( jojimax - jojimin > maxxdiff) THEN
               llfourflag = .TRUE.
               ifourflagcountr(3) = ifourflagcountr(3) + 1
            ENDIF
            IF ( jojjmax - jojjmin > maxydiff) THEN
               llfourflag = .TRUE.
               ifourflagcountr(4) = ifourflagcountr(4) + 1
            ENDIF
            
         ENDIF

         ipmx = 0
         IF (llfourflag) ipmx = 1

         IF (llfourflag) THEN 
            ifourflagcountt = ifourflagcountt + 1
         ELSE
            ifourflagcountf = ifourflagcountf + 1
         ENDIF

         gridpointsn : DO ip = 0, ipmx
            DO jp = 0, ipmx
               
               IF ( kproc(jo) /= -1 ) EXIT gridpointsn
        
               ipx = ipx1 + ip
               ipy = ipy1 + jp
               
               IF (llfourflag) THEN

                  ! deal with wrap around
                  IF ( ipx > nlons ) ipx = 1
                  IF ( ipy > nlats ) ipy = 1
                  IF ( ipx < 1     ) ipx = nlons
                  IF ( ipy < 1     ) ipy = nlats

                  ! get i,j 
                  isx = ixpos(ipx,ipy)
                  isy = iypos(ipx,ipy)
                  
                  ! estimate appropriate search region (use max/min values)
                  jimin = isx - maxxdiff - 1
                  jimax = isx + maxxdiff + 1
                  jjmin = isy - maxydiff - 1
                  jjmax = isy + maxydiff + 1

               ENDIF

               IF ( jimin < 1      ) jimin = 1
               IF ( jimax > jlon-1 ) jimax = jlon-1
               IF ( jjmin < 1      ) jjmin = 1
               IF ( jjmax > jlat-1 ) jjmax = jlat-1
               
               !---------------------------------------------------------------
               ! Ensure that all observation longtiudes are between 0 and 360 
               !---------------------------------------------------------------

               IF ( zplam(jo) <   0.0_wp ) zplam(jo) = zplam(jo) + 360.0_wp
               IF ( zplam(jo) > 360.0_wp ) zplam(jo) = zplam(jo) - 360.0_wp
      
               !---------------------------------------------------------------
               ! Find observations which are on within 1e-6 of a grid point
               !---------------------------------------------------------------

               gridloop: DO jj = jjmin, jjmax
                  DO ji = jimin, jimax
                     IF ( ABS( zphig(ji,jj) - pphi(jo) ) < 1e-6 )  THEN
                        zlam = zlamg(ji,jj)
                        IF ( zlam <   0.0_wp ) zlam = zlam + 360.0_wp
                        IF ( zlam > 360.0_wp ) zlam = zlam - 360.0_wp
                        IF ( ABS( zlam - zplam(jo) ) < 1e-6 ) THEN
                           IF ( llinvalidcell(ji,jj) ) THEN
                              kproc(jo) = nproc + 1000000
                              kobsi(jo) = ji + 1
                              kobsj(jo) = jj + 1
                              CYCLE
                           ELSE
                              kproc(jo) = nproc
                              kobsi(jo) = ji + 1
                              kobsj(jo) = jj + 1
                              EXIT gridloop
                           ENDIF
                        ENDIF
                     ENDIF
                  END DO
               END DO gridloop

               !---------------------------------------------------------------
               ! Ensure that all observation longtiudes are between -180/180
               !---------------------------------------------------------------

               IF ( zplam(jo) > 180 ) zplam(jo) = zplam(jo) - 360.0_wp

               IF ( kproc(jo) == -1 ) THEN
                  
                  ! Normal case 
                  gridpoints : DO jj = jjmin, jjmax
                     DO ji = jimin, jimax


                        IF ( ( zplam(jo) > zlamtmax(ji,jj) ) .OR. &
                           & ( zplam(jo) < zlamtmin(ji,jj) ) ) CYCLE
                        
                        IF ( ABS( pphi(jo) ) < 85 ) THEN
                           IF ( ( pphi(jo) > zphitmax(ji,jj) ) .OR. &
                              & ( pphi(jo) < zphitmin(ji,jj) ) ) CYCLE
                        ENDIF
                        
                        IF ( linquad( zplam(jo), pphi(jo), &
                           &          zlamtm(:,ji,jj), zphitm(:,ji,jj) ) ) THEN
                           IF ( llinvalidcell(ji,jj) ) THEN
                              kproc(jo) = nproc + 1000000
                              kobsi(jo) = ji + 1
                              kobsj(jo) = jj + 1
                              CYCLE
                           ELSE
                              kproc(jo) = nproc
                              kobsi(jo) = ji + 1
                              kobsj(jo) = jj + 1
                              EXIT gridpoints
                           ENDIF
                        ENDIF
                        
                     END DO
                  END DO gridpoints
               ENDIF

               ! In case of failure retry for obs. longtiude + 360.
               IF ( kproc(jo) == -1 ) THEN
                  gridpoints_greenwich : DO jj = jjmin, jjmax
                     DO ji = jimin, jimax
                        
                        IF ( ( zplam(jo)+360.0_wp > zlamtmax(ji,jj) ) .OR. &
                           & ( zplam(jo)+360.0_wp < zlamtmin(ji,jj) ) ) CYCLE

                        IF ( ABS( pphi(jo) ) < 85 ) THEN
                           IF ( ( pphi(jo) > zphitmax(ji,jj) ) .OR. &
                              & ( pphi(jo) < zphitmin(ji,jj) ) ) CYCLE
                        ENDIF

                        IF ( linquad( zplam(jo)+360.0_wp, pphi(jo), &
                           &          zlamtm(:,ji,jj), zphitm(:,ji,jj) ) ) THEN
                           IF ( llinvalidcell(ji,jj) ) THEN
                              kproc(jo) = nproc + 1000000
                              kobsi(jo) = ji + 1
                              kobsj(jo) = jj + 1
                              CYCLE
                           ELSE
                              kproc(jo) = nproc
                              kobsi(jo) = ji + 1
                              kobsj(jo) = jj + 1
                              EXIT gridpoints_greenwich
                           ENDIF
                        ENDIF
                        
                     END DO
                  END DO gridpoints_greenwich
                  
               ENDIF   ! kproc
               
            END DO
         END DO gridpointsn
      END DO gpkobs  ! kobs

      !----------------------------------------------------------------------
      ! Synchronize kproc on all processors
      !----------------------------------------------------------------------
      IF ( ln_grid_global ) THEN
         CALL obs_mpp_max_integer( kproc, kobs )
         CALL obs_mpp_max_integer( kobsi, kobs )
         CALL obs_mpp_max_integer( kobsj, kobs )
      ELSE
         CALL obs_mpp_find_obs_proc( kproc, kobs )
      ENDIF

      WHERE( kproc(:) >= 1000000 )
         kproc(:) = kproc(:) - 1000000
      END WHERE

      DEALLOCATE( &
         & zlamg,         &
         & zphig,         &
         & zmskg,         &
         & zphitmax,      &
         & zphitmin,      &
         & zlamtmax,      &
         & zlamtmin,      &
         & llinvalidcell, &
         & zlamtm,        &
         & zphitm,        &
         & zplam          &
         & )
      
   END SUBROUTINE obs_grd_lookup


   SUBROUTINE obs_grid_setup
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE obs_grid_setup ***
      !!
      !! ** Purpose : Setup a lookup table to reduce the searching required
      !!              for converting lat lons to grid point location
      !!              produces or reads in a preexisting file for use in 
      !!              obs_grid_search_lookup_local
      !!
      !! ** Method : calls obs_grid_search_bruteforce_local with a array 
      !!             of lats and lons
      !!
      !! History :
      !!        !  2007-12 (D. Lea) new routine
      !!----------------------------------------------------------------------
      
      !! * Local declarations
      CHARACTER(LEN=15), PARAMETER :: &
         & cpname = 'obs_grid_setup'
      CHARACTER(LEN=40) :: cfname      
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jk
      INTEGER :: jo
      INTEGER :: idfile, idny, idnx, idxpos, idypos
      INTEGER :: idlat, idlon, fileexist
      INTEGER, DIMENSION(2) :: incdim
      CHARACTER(LEN=20) :: datestr=" ",timestr=" "
      REAL(wp) :: tmpx1, tmpx2, tmpy1, tmpy2
      REAL(wp) :: meanxdiff, meanydiff
      REAL(wp) :: meanxdiff1, meanydiff1
      REAL(wp) :: meanxdiff2, meanydiff2
      INTEGER :: numx1, numx2, numy1, numy2, df
      INTEGER :: jimin, jimax, jjmin, jjmax
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
         & lonsi,     &
         & latsi
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: &  
         & ixposi,    &
         & iyposi,    & 
         & iproci    
      INTEGER, PARAMETER :: histsize=90
      INTEGER, DIMENSION(histsize) :: &
         & histx1, histx2, histy1, histy2
      REAL, DIMENSION(histsize) :: &
         & fhistx1, fhistx2, fhisty1, fhisty2
      REAL(wp) :: histtol
      
      IF (ln_grid_search_lookup) THEN
         
         WRITE(numout,*) 'Calling obs_grid_setup'
         
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'Grid search resolution : ', rn_gridsearchres
         
         gsearch_nlons_def  = NINT( 360.0_wp / rn_gridsearchres ) 
         gsearch_nlats_def  = NINT( 180.0_wp / rn_gridsearchres )
         gsearch_lonmin_def = -180.0_wp + 0.5_wp * rn_gridsearchres
         gsearch_latmin_def =  -90.0_wp + 0.5_wp * rn_gridsearchres
         gsearch_dlon_def   = rn_gridsearchres
         gsearch_dlat_def   = rn_gridsearchres
         
         IF (lwp) THEN
            WRITE(numout,*)'Grid search gsearch_nlons_def  = ',gsearch_nlons_def
            WRITE(numout,*)'Grid search gsearch_nlats_def  = ',gsearch_nlats_def
            WRITE(numout,*)'Grid search gsearch_lonmin_def = ',gsearch_lonmin_def
            WRITE(numout,*)'Grid search gsearch_latmin_def = ',gsearch_latmin_def
            WRITE(numout,*)'Grid search gsearch_dlon_def   = ',gsearch_dlon_def
            WRITE(numout,*)'Grid search gsearch_dlat_def   = ',gsearch_dlat_def
         ENDIF

         IF ( ln_grid_global ) THEN
            WRITE(cfname, FMT="(A,'_',A)") &
               &          TRIM(cn_gridsearchfile), 'global.nc'
         ELSE
            WRITE(cfname, FMT="(A,'_',I4.4,'of',I4.4,'by',I4.4,'.nc')") &
               &          TRIM(cn_gridsearchfile), nproc, jpni, jpnj
         ENDIF

         fileexist=nf90_open( TRIM( cfname ), nf90_nowrite, &
            &                  idfile )
         
         IF ( fileexist == nf90_noerr ) THEN
            
            ! read data
            ! initially assume size is as defined (to be fixed)
            
            WRITE(numout,*) 'Reading: ',cfname
            
            CALL chkerr( nf90_open( TRIM( cfname ), nf90_nowrite, idfile ), &
               &         cpname, __LINE__ )
            CALL chkerr( nf90_get_att( idfile, nf90_global, 'maxxdiff', maxxdiff ), &
               &         cpname, __LINE__ )        
            CALL chkerr( nf90_get_att( idfile, nf90_global, 'maxydiff', maxydiff ), &
               &         cpname, __LINE__ )        
            CALL chkerr( nf90_get_att( idfile, nf90_global, 'dlon', dlon ), &
            &         cpname, __LINE__ )        
            CALL chkerr( nf90_get_att( idfile, nf90_global, 'dlat', dlat ), &
               &         cpname, __LINE__ )        
            CALL chkerr( nf90_get_att( idfile, nf90_global, 'lonmin', lonmin ), &
               &         cpname, __LINE__ )        
            CALL chkerr( nf90_get_att( idfile, nf90_global, 'latmin', latmin ), &
               &         cpname, __LINE__ )        
            
            CALL chkerr( nf90_inq_dimid(idfile, 'nx'  , idnx), &
               &         cpname, __LINE__ )
            CALL chkerr( nf90_inquire_dimension( idfile, idnx, len = nlons ),     &
               &         cpname, __LINE__ ) 
            CALL chkerr( nf90_inq_dimid(idfile, 'ny'  , idny), &
               &         cpname, __LINE__ )
            CALL chkerr( nf90_inquire_dimension( idfile, idny, len = nlats ),     &
               &         cpname, __LINE__ ) 
            
            ALLOCATE( &
               & lons(nlons,nlats),  &
               & lats(nlons,nlats),  &
               & ixpos(nlons,nlats), &
               & iypos(nlons,nlats), &
               & iprocn(nlons,nlats)  &
               & )
            
            CALL chkerr( nf90_inq_varid( idfile, 'XPOS', idxpos ), & 
               &         cpname, __LINE__ )
            CALL chkerr( nf90_get_var  ( idfile, idxpos, ixpos),   &
               &         cpname, __LINE__ )
            CALL chkerr( nf90_inq_varid( idfile, 'YPOS', idypos ), & 
               &         cpname, __LINE__ )
            CALL chkerr( nf90_get_var  ( idfile, idypos, iypos),   &
               &         cpname, __LINE__ )
            
            CALL chkerr( nf90_close( idfile ), cpname, __LINE__ )
            
            ! setup arrays
            
            DO ji = 1, nlons
               DO jj = 1, nlats
                  lons(ji,jj) = lonmin + (ji-1) * dlon
                  lats(ji,jj) = latmin + (jj-1) * dlat
               END DO
            END DO
            
            ! if we are not reading the file we need to create it
            ! create new obs grid search lookup file
            
         ELSE 
            
            ! call obs_grid_search
            
            IF (lwp) THEN
               WRITE(numout,*) 'creating: ',cfname
               WRITE(numout,*) 'calling obs_grid_search: ',nlons*nlats
            ENDIF

            ! set parameters from default values
            nlons  = gsearch_nlons_def
            nlats  = gsearch_nlats_def
            lonmin = gsearch_lonmin_def
            latmin = gsearch_latmin_def
            dlon   = gsearch_dlon_def
            dlat   = gsearch_dlat_def
            
            ! setup arrays
            
            ALLOCATE( &
               & lonsi(nlons,nlats),   &
               & latsi(nlons,nlats),   &
               & ixposi(nlons,nlats),  &
               & iyposi(nlons,nlats),  &
               & iproci(nlons,nlats)   &
               & )
         
            DO ji = 1, nlons
               DO jj = 1, nlats
                  lonsi(ji,jj) = lonmin + (ji-1) * dlon
                  latsi(ji,jj) = latmin + (jj-1) * dlat
               END DO
            END DO
            
            CALL obs_grd_bruteforce( jpi, jpj, jpiglo, jpjglo,  &
               &                     1, nlci, 1, nlcj,          &
               &                     nproc, jpnij,              &
               &                     glamt, gphit, tmask,       &
               &                     nlons*nlats, lonsi, latsi, &
               &                     ixposi, iyposi, iproci )
            
            ! minimise file size by removing regions with no data from xypos file
            ! should be able to just use xpos (ypos will have the same areas of missing data)
         
            jimin=1
            jimax=nlons
            jjmin=1
            jjmax=nlats

            minlon_xpos: DO ji= 1, nlons
               IF (COUNT(ixposi(ji,:) >= 0) > 0) THEN
                  jimin=ji
                  EXIT minlon_xpos
               ENDIF
            END DO minlon_xpos

            maxlon_xpos: DO ji= nlons, 1, -1
               IF (COUNT(ixposi(ji,:) >= 0) > 0) THEN
                  jimax=ji
                  EXIT maxlon_xpos
               ENDIF
            END DO maxlon_xpos

            minlat_xpos: DO jj= 1, nlats
               IF (COUNT(ixposi(:,jj) >= 0) > 0) THEN
                  jjmin=jj
                  EXIT minlat_xpos
               ENDIF
            END DO minlat_xpos

            maxlat_xpos: DO jj= nlats, 1, -1
               IF (COUNT(ixposi(:,jj) >= 0) > 0) THEN 
                  jjmax=jj
                  EXIT maxlat_xpos
               ENDIF
            END DO maxlat_xpos

            lonmin = lonsi(jimin,jjmin)
            latmin = latsi(jimin,jjmin)
            nlons  = jimax-jimin+1
            nlats  = jjmax-jjmin+1

            ! construct new arrays

            ALLOCATE( &
               & lons(nlons,nlats),    &
               & lats(nlons,nlats),    &
               & ixpos(nlons,nlats),   &
               & iypos(nlons,nlats),   &
               & iprocn(nlons,nlats)    &
               & )

            lons(:,:) = lonsi(jimin:jimax,jjmin:jjmax)
            lats(:,:) = latsi(jimin:jimax,jjmin:jjmax)
            ixpos(:,:) = ixposi(jimin:jimax,jjmin:jjmax)
            iypos(:,:) = iyposi(jimin:jimax,jjmin:jjmax)
            iprocn(:,:) = iproci(jimin:jimax,jjmin:jjmax)

            DEALLOCATE(lonsi,latsi,ixposi,iyposi,iproci)

            ! calculate (estimate) maxxdiff, maxydiff
            ! this is to help define the search area for obs_grid_search_lookup

            maxxdiff = 1
            maxydiff = 1

            tmpx1 = 0
            tmpx2 = 0
            tmpy1 = 0
            tmpy2 = 0 

            numx1 = 0
            numx2 = 0
            numy1 = 0
            numy2 = 0

            ! calculate the mean absolute xdiff and ydiff
            ! also calculate a histogram
            ! note the reason why looking for xdiff and ydiff in both directions
            ! is to allow for rotated grids

            DO ji = 1, nlons-1
               DO jj = 1, nlats-1
                  IF ( ixpos(ji,jj) > 0 .AND. iypos(ji,jj) > 0 ) THEN
                     IF ( ixpos(ji+1,jj) > 0 ) THEN
                        df = ABS( ixpos(ji+1,jj) - ixpos(ji,jj) )
                        tmpx1 = tmpx1+df
                        numx1 = numx1+1
                        IF ( df < histsize ) histx1(df+1) = histx1(df+1) + 1
                     ENDIF
                     IF ( ixpos(ji,jj+1) > 0 ) THEN
                        df = ABS( ixpos(ji,jj+1) - ixpos(ji,jj) )
                        tmpx2 = tmpx2 + df
                        numx2 = numx2 + 1
                        IF ( df < histsize ) histx2(df+1) = histx2(df+1) + 1
                     ENDIF
                     IF (iypos(ji+1,jj) > 0) THEN
                        df = ABS( iypos(ji+1,jj) - iypos(ji,jj) )
                        tmpy1 = tmpy1 + df
                        numy1 = numy1 + 1
                        IF ( df < histsize ) histy1(df+1) = histy1(df+1) + 1
                     ENDIF
                     IF ( iypos(ji,jj+1) > 0 ) THEN
                        df = ABS( iypos(ji,jj+1) - iypos(ji,jj) )
                        tmpy2 = tmpy2 + df
                        numy2 = numy2 + 1
                        IF ( df < histsize ) histy2(df+1) = histy2(df+1) + 1
                     ENDIF
                  ENDIF
               END DO
            END DO

            IF (lwp) THEN
               WRITE(numout,*) 'histograms'
               WRITE(numout,*) '0   1   2   3   4   5   6   7   8   9   10 ...'
               WRITE(numout,*) 'histx1'
               WRITE(numout,*) histx1
               WRITE(numout,*) 'histx2'
               WRITE(numout,*) histx2
               WRITE(numout,*) 'histy1'
               WRITE(numout,*) histy1
               WRITE(numout,*) 'histy2'
               WRITE(numout,*) histy2
            ENDIF

            meanxdiff1 = tmpx1 / numx1
            meanydiff1 = tmpy1 / numy1
            meanxdiff2 = tmpx2 / numx2
            meanydiff2 = tmpy2 / numy2

            meanxdiff = MAXVAL((/ meanxdiff1, meanxdiff2 /))
            meanydiff = MAXVAL((/ meanydiff1, meanydiff2 /))

            IF (lwp) THEN
               WRITE(numout,*) tmpx1, tmpx2, tmpy1, tmpy2
               WRITE(numout,*) numx1, numx2, numy1, numy2
               WRITE(numout,*) 'meanxdiff: ',meanxdiff, meanxdiff1, meanxdiff2
               WRITE(numout,*) 'meanydiff: ',meanydiff, meanydiff1, meanydiff2
            ENDIF

            tmpx1 = 0
            tmpx2 = 0
            tmpy1 = 0
            tmpy2 = 0

            numx1 = 0
            numx2 = 0
            numy1 = 0
            numy2 = 0

            histx1(:) = 0
            histx2(:) = 0
            histy1(:) = 0
            histy2(:) = 0

            limxdiff = meanxdiff * 4! limit the difference to avoid picking up wraparound
            limydiff = meanydiff * 4

            DO ji = 1, nlons-1
               DO jj = 1, nlats-1
                  IF ( ixpos(ji,jj) > 0 .AND. iypos(ji,jj) > 0 ) THEN

                     IF ( ixpos(ji+1,jj) > 0 ) THEN
                        df = ABS( ixpos(ji+1,jj)-ixpos(ji,jj) )
                        tmpx1 = df
                        IF ( df < limxdiff ) numx1 = numx1+1
                        IF ( df < histsize ) histx1(df+1) = histx1(df+1) + 1
                     ENDIF
                     IF ( ixpos(ji,jj+1) > 0 ) THEN
                        df = ABS( ixpos(ji,jj+1) - ixpos(ji,jj) )
                        tmpx2 = df
                        IF ( df < limxdiff ) numx2 = numx2 + 1
                        IF ( df < histsize ) histx2(df+1) = histx2(df+1) + 1
                     ENDIF
                     IF (iypos(ji+1,jj) > 0) THEN
                        df = ABS( iypos(ji+1,jj) - iypos(ji,jj) )
                        tmpy1 = df
                        IF ( df < limydiff ) numy1 = numy1 + 1
                        IF ( df < histsize ) histy1(df+1) = histy1(df+1) + 1
                     ENDIF
                     IF (iypos(ji,jj+1) > 0) THEN
                        df = ABS( iypos(ji,jj+1) - iypos(ji,jj) )
                        tmpy2 = df
                        IF ( df < limydiff ) numy2 = numy2+1
                        IF ( df < histsize ) histy2(df+1) = histy2(df+1)+1
                     ENDIF

                     IF ( maxxdiff < tmpx1 .AND. tmpx1 < limxdiff ) &
                        & maxxdiff = tmpx1
                     IF ( maxxdiff < tmpx2 .AND. tmpx2 < limxdiff ) &
                        & maxxdiff = tmpx2
                     IF ( maxydiff < tmpy1 .AND. tmpy1 < limydiff ) &
                        & maxydiff = tmpy1
                     IF ( maxydiff < tmpy2 .AND. tmpy2 < limydiff ) &
                        & maxydiff = tmpy2

                  ENDIF
               END DO
            END DO

            ! cumulative histograms

            DO ji = 1, histsize - 1
               histx1(ji+1) = histx1(ji+1) + histx1(ji)
               histx2(ji+1) = histx2(ji+1) + histx2(ji)
               histy1(ji+1) = histy1(ji+1) + histy1(ji)
               histy2(ji+1) = histy2(ji+1) + histy2(ji)
            END DO

            fhistx1(:) = histx1(:) * 1.0 / numx1
            fhistx2(:) = histx2(:) * 1.0 / numx2
            fhisty1(:) = histy1(:) * 1.0 / numy1
            fhisty2(:) = histy2(:) * 1.0 / numy2

            ! output new histograms

            IF (lwp) THEN
               WRITE(numout,*) 'cumulative histograms'
               WRITE(numout,*) '0   1   2   3   4   5   6   7   8   9   10 ...'
               WRITE(numout,*) 'fhistx1'
               WRITE(numout,*) fhistx1
               WRITE(numout,*) 'fhistx2'
               WRITE(numout,*) fhistx2
               WRITE(numout,*) 'fhisty1'
               WRITE(numout,*) fhisty1
               WRITE(numout,*) 'fhisty2'
               WRITE(numout,*) fhisty2
            ENDIF

            ! calculate maxxdiff and maxydiff based on cumulative histograms
            ! where > 0.999 of points are

            ! maxval just converts 1x1 vector return from maxloc to a scalar 

            histtol = 0.999
            tmpx1 = MAXVAL( MAXLOC( fhistx1(:), mask = ( fhistx1(:) <= histtol ) ) )
            tmpx2 = MAXVAL( MAXLOC( fhistx2(:), mask = ( fhistx2(:) <= histtol ) ) )
            tmpy1 = MAXVAL( MAXLOC( fhisty1(:), mask = ( fhisty1(:) <= histtol ) ) )
            tmpy2 = MAXVAL( MAXLOC( fhisty2(:), mask = ( fhisty2(:) <= histtol ) ) )

            maxxdiff = MAXVAL( (/ tmpx1, tmpx2 /) ) + 1
            maxydiff = MAXVAL( (/ tmpy1, tmpy2 /) ) + 1

            ! Write out data

            IF ( ( .NOT. ln_grid_global ) .OR. &
               & ( ( ln_grid_global ) .AND. ( nproc==0 ) ) ) THEN

               CALL chkerr( nf90_create (TRIM(cfname), nf90_clobber, idfile), &
                  &         cpname, __LINE__ )
               CALL chkerr( nf90_put_att( idfile, nf90_global, 'title',       &
                  &         'Mapping file from lon/lat to model grid point' ),&
                  &         cpname,__LINE__ ) 
               CALL chkerr( nf90_put_att( idfile, nf90_global, 'maxxdiff',    &
                  &                       maxxdiff ),                         &
                  &         cpname,__LINE__ ) 
               CALL chkerr( nf90_put_att( idfile, nf90_global, 'maxydiff',    &
                  &                       maxydiff ),                         &
                  &         cpname,__LINE__ ) 
               CALL chkerr( nf90_put_att( idfile, nf90_global, 'dlon', dlon ),&
                  &         cpname,__LINE__ ) 
               CALL chkerr( nf90_put_att( idfile, nf90_global, 'dlat', dlat ),&
                  &         cpname,__LINE__ ) 
               CALL chkerr( nf90_put_att( idfile, nf90_global, 'lonmin',      &
                  &                       lonmin ),                           &
                  &         cpname,__LINE__ ) 
               CALL chkerr( nf90_put_att( idfile, nf90_global, 'latmin',      &
                  &                       latmin ),                           &
                  &         cpname,__LINE__ ) 

               CALL chkerr( nf90_def_dim(idfile, 'nx'  , nlons, idnx),        &
                  &         cpname,__LINE__ )
               CALL chkerr( nf90_def_dim(idfile, 'ny'  , nlats, idny),        &
                  &         cpname,__LINE__ )

               incdim(1) = idnx
               incdim(2) = idny
               
               CALL chkerr( nf90_def_var( idfile, 'LON', nf90_float, incdim,  &
                  &                       idlon ),                            &
                  &         cpname, __LINE__ )
               CALL chkerr( nf90_put_att( idfile, idlon, 'long_name',         &
                  &                       'longitude' ),                      &
                  &         cpname, __LINE__ )
               
               CALL chkerr( nf90_def_var( idfile, 'LAT', nf90_float, incdim,  &
                  &                       idlat ),                            &
                  &         cpname, __LINE__ )
               CALL chkerr( nf90_put_att( idfile, idlat, 'long_name',         &
                  &                       'latitude' ),                       &
                  &         cpname, __LINE__ )

               CALL chkerr( nf90_def_var( idfile, 'XPOS', nf90_int, incdim,   &
                  &                       idxpos ),                           &
                  &         cpname, __LINE__ )
               CALL chkerr( nf90_put_att( idfile, idxpos, 'long_name',        &
                  &                       'x position' ),                     &
                  &         cpname, __LINE__ )
               CALL chkerr( nf90_put_att( idfile, idxpos, '_FillValue', -1 ), &
                  &         cpname, __LINE__ )

               CALL chkerr( nf90_def_var( idfile, 'YPOS', nf90_int, incdim,   &
                  &                       idypos ),                           &
                  &         cpname, __LINE__ )
               CALL chkerr( nf90_put_att( idfile, idypos, 'long_name',        &
                  &                       'y position' ),                     &
                  &         cpname, __LINE__ )
               CALL chkerr( nf90_put_att( idfile, idypos, '_FillValue', -1 ), &
                  &         cpname, __LINE__ )

               CALL chkerr( nf90_enddef( idfile ), cpname, __LINE__ )
               
               CALL chkerr( nf90_put_var( idfile, idlon, lons),               &
                  &         cpname, __LINE__ )
               CALL chkerr( nf90_put_var( idfile, idlat, lats),               &
                  &         cpname, __LINE__ )
               CALL chkerr( nf90_put_var( idfile, idxpos, ixpos),             &
                  &         cpname, __LINE__ )
               CALL chkerr( nf90_put_var( idfile, idypos, iypos),             &
                  &         cpname, __LINE__ )
               
               CALL chkerr( nf90_close( idfile ), cpname, __LINE__ )
               
               ! should also output max i, max j spacing for use in 
               ! obs_grid_search_lookup
               
            ENDIF

         ENDIF

      ENDIF

   END SUBROUTINE obs_grid_setup
   
   SUBROUTINE obs_grid_deallocate( )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE obs_grid_setup ***
      !!
      !! ** Purpose : Deallocate arrays setup by obs_grid_setup
      !!
      !! History :
      !!        !  2007-12 (D. Lea) new routine
      !!-----------------------------------------------------------------------

      IF (ln_grid_search_lookup) THEN
         DEALLOCATE( lons, lats, ixpos, iypos, iprocn )
      ENDIF
      
   END SUBROUTINE obs_grid_deallocate

#include "obs_level_search.h90"

#include "linquad.h90"

#include "maxdist.h90"

#include "find_obs_proc.h90"

END MODULE obs_grid

