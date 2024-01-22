MODULE obs_averg_h2d
   !!======================================================================
   !!                       ***  MODULE obs_averg_h2d   ***
   !! Observation diagnostics: Perform the horizontal averaging
   !!                          from model grid to observation footprint
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   obs_averg_h2d     : Horizontal averaging to the observation footprint
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind, ONLY : &  ! Precision variables
      & wp
   USE par_oce, ONLY : & 
      & jpi, jpj
   USE phycst,   ONLY : &  ! Physical constants
      & rad,  &
      & ra,   &
      & rpi
   USE dom_oce,   ONLY : &
      & e1t, e2t, &
      & e1f, e2f, &
      & glamt, gphit, &
      & nproc
   USE in_out_manager
   USE obs_const, ONLY : &
      & obfillflt		! Fillvalue
   USE obs_utils           ! Utility functions
   USE lib_mpp,   ONLY : &
      & ctl_warn, ctl_stop, &
      & mpp_min, lk_mpp

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE obs_avg_h2d_rad, & ! Horizontal averaging using a radial footprint
      &    obs_avg_h2d_rec, & ! Horizontal averaging using a rectangular footprint
      &    obs_deg2dist,    & ! Conversion of distance in degrees to distance in metres
      &    obs_dist2corners   ! Distance from the centre of obs footprint to the corners of a grid box
   
   PUBLIC obs_avg_h2d,      & ! Horizontal averaging to the observation footprint
      &   obs_avg_h2d_init, & ! Set up weights for the averaging
      &   obs_max_fpsize      ! Works out the maximum number of grid points required for the averaging

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS
   SUBROUTINE obs_avg_h2d_init( kpk, kpk2, kmaxifp, kmaxjfp, k2dint, plam,  pphi, &
      &                         pglam, pgphi, pglamf, pgphif, pmask, plamscl, pphiscl, lindegrees, &
      &                         pweig, pobsmask, iminpoints )
      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_avg_h2d_init  ***
      !!
      !! ** Purpose : Computes weights for horizontal averaging to the 
      !!              observation footprint.
      !!
      !! ** Method  : Horizontal averaging to the observation footprint using
      !!              model values at a defined area.
      !!
      !!    Averaging schemes :
      !!
      !!    Two horizontal averaging schemes are available:
      !!        - weighted radial footprint        (k2dint = 5)
      !!        - weighted rectangular footprint   (k2dint = 6)
      !!
      !! History :
      !!        ! 13-10 (M. Martin)
      !!-----------------------------------------------------------------------
      !! * Modules used
      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kpk,   &             ! Parameter values for automatic arrays
         & kpk2,  &
         & kmaxifp,  &          ! Max size of model points in i-direction for obs footprint
         & kmaxjfp,  &          ! Max size of model points in j-direction for obs footprint
         & k2dint               ! Averaging scheme options - see header
      REAL(KIND=wp), INTENT(INOUT) :: &
         & plam, &              ! Geographical (lat,lon) coordinates of
         & pphi                 ! observation
      REAL(KIND=wp), DIMENSION(kmaxifp,kmaxjfp), INTENT(IN) :: &
         & pglam, &             ! Model variable lon
         & pgphi                ! Model variable lat
      REAL(KIND=wp), DIMENSION(kmaxifp+1,kmaxjfp+1), INTENT(IN) :: &
         & pglamf, &            ! Model variable lon at corners of grid-boxes
         & pgphif               ! Model variable lat at corners of grid-boxes
      REAL(KIND=wp), DIMENSION(kmaxifp,kmaxjfp,kpk2), INTENT(IN) :: &
         & pmask                ! Model variable mask
      REAL(KIND=wp), INTENT(IN) :: &
         & plamscl, &           ! Diameter (lat,lon) of obs footprint in metres
         & pphiscl              ! This is the full width (rather than half-width)
      LOGICAL, INTENT(IN) :: &
         & lindegrees           ! T=> obs footprint specified in degrees, F=> in metres
      REAL(KIND=wp), DIMENSION(kmaxifp,kmaxjfp,kpk2), INTENT(OUT) ::  &
         & pweig                ! Weights for averaging
      REAL(KIND=wp), DIMENSION(kpk2), INTENT(OUT) ::  &
         & pobsmask             ! Vertical mask for observations
      INTEGER, INTENT(IN), OPTIONAL :: &
         & iminpoints           ! Reject point which is not surrounded
                                ! by at least iminpoints sea points

      !! * Local declarations
      INTEGER :: &
         & jk
      INTEGER :: &
         & ikmax


      !------------------------------------------------------------------------
      !
      !------------------------------------------------------------------------

      !------------------------------------------------------------------------
      ! Initialize number of levels
      !------------------------------------------------------------------------
      IF ( kpk2 == 1 ) THEN
         ikmax = 1
      ELSEIF ( kpk2 == kpk) THEN
         ikmax = kpk-1
      ENDIF


         SELECT CASE (k2dint)
         CASE(5)
            CALL obs_avg_h2d_rad( kpk2, ikmax, kmaxifp, kmaxjfp, plam, pphi, &
      &                           plamscl, pphiscl, lindegrees, pmask, pglam, pgphi, pglamf, pgphif, pweig )
         CASE(6)
            CALL obs_avg_h2d_rec( kpk2, ikmax, kmaxifp, kmaxjfp, plam, pphi, &
      &                           plamscl, pphiscl, lindegrees, pmask, pglam, pgphi, pglamf, pgphif, pweig )
         END SELECT


      END SUBROUTINE obs_avg_h2d_init


      SUBROUTINE obs_avg_h2d_rad( kpk2, kmax, kmaxifp, kmaxjfp, plam, pphi, &
      &                           plamscl, pphiscl, lindegrees, pmask, pglam, pgphi, pglamf, pgphif, pweig )
      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_avg_h2d_rad  ***
      !!
      !! ** Purpose : Computes weights for horizontal averaging to the 
      !!              observation using a radial footprint.
      !!
      !! ** Method  : Calculate whether each grid box is completely or 
      !!              partially within the observation footprint.
      !!              If it is partially in the footprint then calculate
      !!              the ratio of the area inside the footprint to the total 
      !!              area of the grid box.
      !!
      !! History :
      !!        ! 14-01 (M. Martin)
      !!-----------------------------------------------------------------------
      !! * Modules used
      USE phycst,   ONLY : &  ! Physical constants
         & ra,  &
         & rpi

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kpk2, &             ! Parameter values for automatic arrays
         & kmax

      INTEGER, INTENT(IN) :: &
         & kmaxifp,   &         ! Max size of model points in i-direction for obs footprint
         & kmaxjfp              ! Max size of model points in j-direction for obs footprint

      REAL(KIND=wp), INTENT(IN) :: &
         & plam, &
         & pphi                 ! Geographical (lat,lon) coordinates of
                                ! observation
      REAL(KIND=wp), INTENT(IN) :: &
         & plamscl, &           ! Diameter (lat,lon) of obs footprint in metres or degrees (see below)
         & pphiscl              ! This is the full width (rather than half-width)
      LOGICAL, INTENT(IN) :: &
         & lindegrees           ! T=>scales specified in degrees, F=> in metres
      REAL(KIND=wp), DIMENSION(kmaxifp,kmaxjfp,kpk2), INTENT(IN) :: &
         & pmask                ! Model variable mask
      REAL(KIND=wp), DIMENSION(kmaxifp,kmaxjfp), INTENT(IN) :: &
         & pglam, &             ! Model variable lon
         & pgphi                ! Model variable lat
      REAL(KIND=wp), DIMENSION(kmaxifp+1,kmaxjfp+1), INTENT(IN) :: &
         & pglamf, &             ! Model variable lon at corners of grid boxes
         & pgphif                ! Model variable lat at corners of grid boxes
      REAL(KIND=wp), DIMENSION(kmaxifp,kmaxjfp,kpk2), INTENT(OUT) ::  &
         & pweig                ! Weights for interpolation

      !! Local declarations
      INTEGER :: ji, jj, jk
      INTEGER :: jvert, jis, jjs
      INTEGER :: jnumvert, jnumvertbig
      INTEGER, PARAMETER :: &
         & jnumsubgrid = 20     ! The number of sub grid-boxes (in x and y directions) used to approximate area of obs fp

      REAL(KIND=wp), DIMENSION(4) :: &
         & zxvert, zyvert, &    ! The lon/lat of the vertices(corners) of the grid box in m relative to centre of obs fp
         & zdist                ! Distance of each vertex to the centre of the obs footprint
      REAL(KIND=wp), DIMENSION(4) :: &
         & zxgrid, zygrid, &      ! Distance of each vertex of grid box to the centre of the grid box in x/y directions
         & zdgrid
      REAL(KIND=wp) :: &
         & zdx, zdy, &          ! The sub grid-box sizes (in metres)
         & zarea_subbox, &      ! The area of each sub grid-box (in metres squared)
         & zxpos, zypos, &      ! The x,y position (relative to centre of obs footprint) of the centre of each sub grid-box
         & zsubdist, &          ! The distance of the centre of each sub grid-box from the centre of the obs footprint
         & zarea_fp, &          ! Total area of obs footprint within the grid box
         & zareabox             ! Total area of the grid box
      REAL(KIND=wp) :: &
         & zphiscl_m, &         ! Diameter of obs footprint in metres
         & zlamscl_m
      !---------------------------------------------------------------------------------------------------
      !Initialise weights to zero.
      pweig(:,:,:) = 0.0_wp
      
      !Two footprint sizes can be specified in the namelist but this routine assumes a circular footprint. 
      !If the two sizes are different then write out a warning.
      IF ( pphiscl /= plamscl ) THEN
               CALL ctl_warn( 'obs_avg_h2d_rad:',   &
                  &           'The two components of the obs footprint size are not equal', &
                  &           'yet the radial option has been selected - using pphiscl here' )
      ENDIF
      
      DO jk = 1, kmax
         DO ji = 1, kmaxifp
            DO jj = 1, kmaxjfp
            
               IF ( pmask(ji,jj,jk) == 1.0_wp ) THEN

                  IF ( lindegrees ) THEN
                     !If the scales are specified in degrees, work out the
                     !scales (metres) in x/y directions
                     CALL obs_deg2dist( 1, 1, pglam(ji,jj), pgphi(ji,jj), &
                        &               plamscl, pphiscl, zlamscl_m, zphiscl_m )
                  ELSE
                     zphiscl_m = pphiscl
                  ENDIF


                  ! Work out the area of the grid box using distance of corners relative to centre of grid box
                  CALL obs_dist2corners(pglamf(ji,jj), pglamf(ji+1,jj), pglamf(ji,jj+1), pglamf(ji+1,jj+1), &
                     &                  pgphif(ji,jj), pgphif(ji+1,jj), pgphif(ji,jj+1), pgphif(ji+1,jj+1), &
                     &                  pglam(ji,jj), pgphi(ji,jj), zxgrid, zygrid, zdgrid)
                  zareabox = ABS( zxgrid(1) - zxgrid(2) ) * ABS( zygrid(1) - zygrid(4) )

                  !1. Determine how many of the vertices of the grid box lie within the circle
                  
                  !For each vertex, calculate its location and distance relative 
                  !to the centre of the observation footprint
                  
                  CALL obs_dist2corners(pglamf(ji,jj), pglamf(ji+1,jj), pglamf(ji,jj+1), pglamf(ji+1,jj+1), &
                     &                  pgphif(ji,jj), pgphif(ji+1,jj), pgphif(ji,jj+1), pgphif(ji+1,jj+1), &
                     &                  plam, pphi, zxvert, zyvert, zdist)

                  jnumvert = 0
                  jnumvertbig = 0
                  DO jvert = 1, 4

                     !If the distance to the center to the observation footprint is less
                     !than the radius of the footprint (half the diameter) then this
                     !vertex is within the observation footprint
                     IF ( zdist(jvert) <= ( zphiscl_m / 2.0_wp ) ) jnumvert = jnumvert + 1
                     
                     !For expediency, check if the vertices are "nearly" within the obs
                     !footprint as if none of them are close to the edge of the footprint
                     !then the footprint is unlikely to be intersecting the grid box
                     IF ( zdist(jvert) - ( 0.5_wp * zareabox ) <= ( zphiscl_m / 2.0 ) ) &
                        & jnumvertbig = jnumvertbig + 1
                     
                  END DO
                  
                  !2. If none of the vertices are even close to the edge of the obs
                  !footprint then leave weight as zero and cycle to next grid box.
                  IF ( jnumvertbig == 0 ) CYCLE

                  !3. If all the vertices of the box are within the observation footprint then the 
                  !      whole grid box is within the footprint so set the weight to one and 
                  !      move to the next grid box.
                  IF ( jnumvert == 4 ) THEN
                     pweig(ji,jj,jk) = 1.0_wp
                     CYCLE
                  ENDIF


                  !4. Use a brute force technique for calculating the area within 
                  !   the grid box covered by the obs footprint.
                  !  (alternative could be to use formulae on 
                  !         http://mathworld.wolfram.com/Circle-LineIntersection.html)
                  !   For now split the grid box into a specified number of smaller 
                  !   boxes and add up the area of those whose centre is within the obs footprint.
                  !   Order of vertices is 1=topleft, 2=topright, 3=bottomright, 4=bottomleft
                  zdx = ABS( zxvert(3) - zxvert(4) ) / REAL(jnumsubgrid, wp)
                  zdy = ABS( zyvert(1) - zyvert(4) ) / REAL(jnumsubgrid, wp)
                  zarea_subbox = zdx * zdy

                  zarea_fp = 0.0_wp
                  DO jis = 1, jnumsubgrid
                     zxpos = zxvert(4) + ( REAL(jis, wp) * zdx ) - (0.5_wp * zdx )
                     DO jjs = 1, jnumsubgrid
                        !Find the distance of the centre of this sub grid box to the
                        !centre of the obs footprint
                        zypos = zyvert(4) + ( REAL(jjs, wp) * zdy ) - ( 0.5_wp * zdy )
                        zsubdist = SQRT( (zxpos * zxpos) + (zypos * zypos) )
                        IF ( zsubdist < ( zphiscl_m / 2.0_wp ) ) &
                           &  zarea_fp = zarea_fp + zarea_subbox
                     END DO
                  END DO

                  !6. Calculate the ratio of the area of the footprint within the box 
                  !   to the total area of the grid box and use this fraction to weight 
                  !   the model value in this grid box.
                  pweig(ji,jj,jk) = MIN( zarea_fp / zareabox, 1.0_wp )

               END IF  !pmask
            END DO
         END DO
      END DO
      
      END SUBROUTINE obs_avg_h2d_rad


      SUBROUTINE obs_avg_h2d_rec( kpk2, kmax, kmaxifp, kmaxjfp, plam, pphi, &
      &                           plamscl, pphiscl, lindegrees, pmask, pglam, pgphi, pglamf, pgphif, pweig )
      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_avg_h2d_rec  ***
      !!
      !! ** Purpose : Computes weights for horizontal averaging to the
      !!              observation using a rectangular footprint which
      !!              is aligned with lines of lat/lon.
      !!
      !! ** Method  : Horizontal averaging to the observation footprint using
      !!              model values at a defined area.
      !!
      !! History :
      !!        ! 14-01 (M. Martin)
      !!-----------------------------------------------------------------------
      !! * Modules used
      USE phycst,   ONLY : &    ! Physical constants
         & ra,  &
         & rpi

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kpk2, &              ! Parameter values for automatic arrays
         & kmax

      INTEGER, INTENT(IN) :: &
         & kmaxifp,   &         ! Max size of model points in i-direction for obs footprint
         & kmaxjfp              ! Max size of model points in j-direction for obs footprint

      REAL(KIND=wp), INTENT(IN) :: &
         & plam, &
         & pphi                 ! Geographical (lat,lon) coordinates of
                                ! observation
      REAL(KIND=wp), INTENT(IN) :: &
         & plamscl, &
         & pphiscl              ! Width in x/y directions of obs footprint in metres
                                ! This is the full width (rather than half-width)
      LOGICAL, INTENT(IN) :: &
         & lindegrees           !T=> scales specified in degrees, F=> in metres
      REAL(KIND=wp), DIMENSION(kmaxifp,kmaxjfp,kpk2), INTENT(IN) :: &
         & pmask                ! Model variable mask
      REAL(KIND=wp), DIMENSION(kmaxifp,kmaxjfp), INTENT(IN) :: &
         & pglam, &             ! Model variable lat at centre of grid boxes
         & pgphi                ! Model variable lon at centre of grid boxes
      REAL(KIND=wp), DIMENSION(kmaxifp+1,kmaxjfp+1), INTENT(IN) :: &
         & pglamf, &             ! Model variable lat at corners of grid boxes
         & pgphif                ! Model variable lon at corners of grid boxes
      REAL(KIND=wp), DIMENSION(kmaxifp,kmaxjfp,kpk2), INTENT(OUT) ::  &
         & pweig                ! Weights for interpolation

      !! Local declarations
      INTEGER :: ji, jj, jk
      INTEGER :: jvert
      INTEGER, DIMENSION(4) :: &
         & jnumvert
      REAL(KIND=wp), DIMENSION(4) :: &
         & zxvert, zyvert         ! The lon/lat of the vertices(corners) of the grid box in m relative to centre of obs fp
      REAL(KIND=wp), DIMENSION(4) :: &
         & zdist                  ! Distance of each vertex to the centre of the obs footprint
      REAL(KIND=wp), DIMENSION(4) :: &
         & zxgrid, zygrid, &      ! Distance of each vertex of grid box to the centre of the grid box in x/y directions
         & zdgrid
      REAL(KIND=wp) :: &
         & zareabox, &            ! Total area of grid box
         & zarea_fp, &            ! Total area of obs footprint
         & zarea_intersect        ! Area of the intersection between the grid box and the obs footprint
      REAL(KIND=wp) :: &
         & zlamscl_m, &
         & zphiscl_m              ! Total width (lat,lon) of obs footprint in metres
      REAL(KIND=wp) :: &
         & z_awidth, z_aheight, & ! Width and height of model grid box
         & z_cwidth, z_cheight    ! Width and height of union of model grid box and obs footprint
      REAL(KIND=wp) :: &
         & zleft, zright, &       ! Distance (metres) of corners area of intersection
         & ztop, zbottom          ! between grid box and obs footprint

      !-----------------------------------------------------------------------

      !Initialise weights to zero
      pweig(:,:,:) = 0.0_wp
      
      !Loop over the grid boxes which have been identified as potentially being within the
      !observation footprint
      DO jk = 1, kmax
         DO ji = 1, kmaxifp
            DO jj = 1, kmaxjfp

               IF ( pmask(ji,jj,jk) == 1.0_wp ) THEN


                  IF ( lindegrees ) THEN
                     !If the scales are specified in degrees, work out the
                     !scales (metres) in x/y directions
                     CALL obs_deg2dist( 1, 1, pglam(ji,jj), pgphi(ji,jj), &
                        &               plamscl, pphiscl, zlamscl_m, zphiscl_m )
                  ELSE
                     zlamscl_m = plamscl
                     zphiscl_m = pphiscl
                  ENDIF

                  ! Work out the area of the grid box using distance of corners relative to centre of grid box
                  CALL obs_dist2corners(pglamf(ji,jj), pglamf(ji+1,jj), pglamf(ji,jj+1), pglamf(ji+1,jj+1), &
                     &                  pgphif(ji,jj), pgphif(ji+1,jj), pgphif(ji,jj+1), pgphif(ji+1,jj+1), &
                     &                  pglam(ji,jj), pgphi(ji,jj), zxgrid, zygrid, zdgrid)

                  !Calculate width and height of model grid box
                  z_awidth  = ABS( zxgrid(1) - zxgrid(2) )
                  z_aheight = ABS( zygrid(1) - zygrid(4) )
                  zareabox = z_awidth * z_aheight

                  ! Work out area of the observation footprint
                  zarea_fp = zlamscl_m * zphiscl_m

                  ! For each corner of the grid box, calculate its location and distance relative
                  ! to the centre of the observation footprint
                  CALL obs_dist2corners(pglamf(ji,jj), pglamf(ji+1,jj), pglamf(ji,jj+1), pglamf(ji+1,jj+1), &
                     &                  pgphif(ji,jj), pgphif(ji+1,jj), pgphif(ji,jj+1), pgphif(ji+1,jj+1), &
                     &                  plam, pphi, zxvert, zyvert, zdist)

                  !Work out maximum width and height of rectangle covered by corners of obs fp and grid box
                  z_cwidth  = MAX( zxvert(1), zxvert(2), -zlamscl_m/2.0_wp, zlamscl_m/2.0_wp ) - &
                     &       MIN( zxvert(1), zxvert(2), -zlamscl_m/2.0_wp, zlamscl_m/2.0_wp )
                     
                  z_cheight = MAX( zyvert(1), zyvert(4), zphiscl_m/2.0_wp, -zphiscl_m/2.0_wp ) - &
                     &       MIN( zyvert(1), zyvert(4), zphiscl_m/2.0_wp, -zphiscl_m/2.0_wp )
                  
                  IF ( ( z_cwidth  >= z_awidth  + zlamscl_m  ) .OR. &
                     & ( z_cheight >= z_aheight + zphiscl_m ) ) THEN
                     !The obs footprint and the model grid box don't overlap so set weight to zero
                     pweig(ji,jj,jk) = 0.0_wp
                  ELSE IF ( ( z_cwidth  == zlamscl_m ) .AND. &
                     &      ( z_cheight == zphiscl_m ) ) THEN
                     !The grid box is totally contained within the obs footprint so set weight to one
                     pweig(ji,jj,jk) = 1.0_wp
                  ELSE IF ( ( z_cwidth  == z_awidth  ) .AND. &
                     &      ( z_cheight == z_aheight ) ) THEN
                     !The obs footprint is totally contained within the grid box so set weight as ratio of the two
                     pweig(ji,jj,jk) = zarea_fp / zareabox
                  ELSE 
                     !The obs footprint and the grid box overlap so calculate the area of the intersection of the two
                     zleft   = max(zxvert(1), -zlamscl_m/2.0_wp)
                     zright  = min(zxvert(2),  zlamscl_m/2.0_wp)
                     zbottom = max(zyvert(4), -zphiscl_m/2.0_wp)
                     ztop    = min(zyvert(1),  zphiscl_m/2.0_wp)
                     
                     IF (  ( zleft < zright ) .AND. ( zbottom < ztop ) ) THEN
                        zarea_intersect = ( zright - zleft ) * ( ztop - zbottom )
                        pweig(ji,jj,jk) = zarea_intersect / zareabox
                     ENDIF
                  ENDIF

               END IF !pmask
            END DO
         END DO
      END DO

   END SUBROUTINE obs_avg_h2d_rec

   SUBROUTINE obs_avg_h2d( kpk, kpk2, kmaxifp, kmaxjfp, pweig, pmod, pobsk )

      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_int_h2d  ***
      !!
      !! ** Purpose : Horizontal averaging to the observation footprint.
      !!
      !! ** Method  : Average the model points based on the weights already calculated.
      !!
      !! ** Action  :
      !!                   
      !! References : 
      !!
      !! History :
      !!        ! 13/10. M. Martin.
      !!-----------------------------------------------------------------------
      !! * Modules used
      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kpk,   &             ! Parameter values for automatic arrays
         & kpk2
      INTEGER, INTENT(IN) :: &
         & kmaxifp,   &         ! Max size of model points in i-direction for obs footprint
         & kmaxjfp              ! Max size of model points in j-direction for obs footprint
      REAL(KIND=wp), DIMENSION(kmaxifp,kmaxjfp,kpk2), INTENT(IN) :: &
         & pweig                ! Interpolation weights
      REAL(KIND=wp), DIMENSION(kmaxifp,kmaxjfp,kpk2), INTENT(IN) :: &
         & pmod                 ! Model variable to interpolate
      REAL(KIND=wp), DIMENSION(kpk2), INTENT(OUT) ::  &
         & pobsk                ! Model profile interpolated to obs (i,j) pt

      INTEGER :: &
         & jk
      INTEGER :: &
         & ikmax
      REAL(KIND=wp) :: &
         & zsum

      !------------------------------------------------------------------------
      ! Initialize number of levels
      !------------------------------------------------------------------------
      IF ( kpk2 == 1 ) THEN
         ikmax = 1
      ELSEIF ( kpk2 == kpk) THEN
         ikmax = kpk-1
      ENDIF

      !------------------------------------------------------------------------
      ! Average model values to the observation footprint
      !------------------------------------------------------------------------
      pobsk = obfillflt

      DO jk = 1, ikmax

         zsum = SUM( pweig(:,:,jk) )

         IF ( zsum /= 0.0_wp ) THEN
            pobsk(jk) = SUM ( pweig(:,:,jk) * pmod(:,:,jk), Mask=pweig(:,:,jk) > 0.0_wp )
            pobsk(jk) = pobsk(jk) / zsum
         END IF

      END DO

   END SUBROUTINE obs_avg_h2d

   SUBROUTINE obs_max_fpsize( k2dint, plamscl, pphiscl, lindegrees, pmask, kmaxifp, kmaxjfp )
      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_max_fpsize  ***
      !!
      !! ** Purpose : Calculate maximum number of grid points which may
      !!              need to be used in the averaging in the global domain.
      !!
      !!
      !! ** Method  : Work out the minimum grid size and work out
      !!              how many of the smallest grid points would be needed
      !!              to cover the scale of the observation footprint.
      !!              This needs to be done using the max/min of the global domain
      !!              as the obs can be distributed from other parts of the grid.
      !!
      !! ** Action  :
      !!                   
      !! References : 
      !!
      !! History :
      !!        ! 14/01. M. Martin. 
      !!-----------------------------------------------------------------------
      !! * Modules used
      !! * Arguments
      INTEGER , INTENT(IN) :: &
         & k2dint                   !Type of interpolation/averaging used
      REAL(KIND=wp), INTENT(IN) :: &
         & plamscl,   &             !Total width/radius in metres of the observation footprint
         & pphiscl                  !
      LOGICAL, INTENT(IN) :: &
         & lindegrees               !T=> plamscl and pphiscl are specified in degrees
      REAL(KIND=wp), DIMENSION(jpi,jpj), INTENT(IN) :: &
         & pmask                    !Land/sea mask
                                    !F=> plamscl and pphiscl are specified in metres
      INTEGER, INTENT(OUT)  :: &
         & kmaxifp,   &             !Max number of grid points in i,j directions to use in averaging
         & kmaxjfp                  !these have to be even so that the footprint is centred

      !! * Local variables
      REAL(KIND=wp) :: &
         & ze1min,     &            !Minimum global grid-size in i,j directions
         & ze2min
      REAL(KIND=wp) :: &
         & zphiscl_m, &
         & zlamscl_m
      !------------------------------------------------------------------------

      IF ( k2dint <= 4 ) THEN
         !If interpolation is being used then only need to use a 2x2 footprint
         kmaxifp = 2
         kmaxjfp = 2

      ELSE

         IF ( lindegrees ) THEN
            !If the scales are specified in degrees, work out the max 
            !distance (metres) in x/y directions
            CALL obs_deg2dist( jpi, jpj, glamt, gphit, &
               &               plamscl, pphiscl, zlamscl_m, zphiscl_m )
         ELSE
            zlamscl_m = plamscl
            zphiscl_m = pphiscl
         ENDIF

         ze1min = MAXVAL( e1t(:,:)*pmask(:,:) )
         ze2min = MAXVAL( e2t(:,:)*pmask(:,:) )

         ze1min = ze1min - MAXVAL( (ze1min-e1t(:,:))*pmask(:,:) )
         ze2min = ze2min - MAXVAL( (ze2min-e2t(:,:))*pmask(:,:) )
         
         IF(lk_mpp) THEN
            CALL mpp_min( 'obs_averg_h2d', ze1min )
            CALL mpp_min( 'obs_averg_h2d', ze2min )
         ENDIF

         kmaxifp = ceiling(zlamscl_m/ze1min) + 1
         kmaxjfp = ceiling(zphiscl_m/ze2min) + 1
         
         !Ensure that these numbers are even
         kmaxifp = kmaxifp + MOD(kmaxifp,2)
         kmaxjfp = kmaxjfp + MOD(kmaxjfp,2)
         

      ENDIF

   END SUBROUTINE obs_max_fpsize

   SUBROUTINE obs_deg2dist( ki, kj, pglam, pgphi, plamscl_deg, pphiscl_deg, &
      &                     plamscl_max, pphiscl_max )
      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_deg2dist  ***
      !!
      !! ** Purpose : Calculate the maximum distance in m of the length scale
      !!              in degrees.
      !!
      !! ** Method  : At each lon/lat point, work out the distances in the 
      !!              zonal and meridional directions.
      !!
      !! ** Action  :
      !!                   
      !! References : 
      !!
      !! History :
      !!        ! 14/01. M. Martin. 
      !!-----------------------------------------------------------------------
      !! * Modules used
      !! * Arguments
      INTEGER , INTENT(IN) :: &
         & ki, kj                   !x/y dimensions of input lat/lon variables
      REAL(KIND=wp), INTENT(IN), DIMENSION(ki,kj) :: &
         & pglam, pgphi             !Longitude and latitudes of grid points
      REAL(KIND=wp), INTENT(IN) :: &
         & plamscl_deg,   &         !Size in degrees of the observation footprint
         & pphiscl_deg              !
      REAL(KIND=wp), INTENT(OUT) :: &
         & plamscl_max, &           !Maximum size of obs footprint in metres
         & pphiscl_max
      
      !! * Local declarations
      INTEGER :: &
         & ji, jj                   !Counters
      REAL(KIND=wp) :: &
         & zlon1, zlon2, &          !Lon values surrounding centre of grid point
         & zlat1, zlat2, &          !Lat values surrounding centre of grid point
         & zdlat, zdlon             !Distance in radians in lat/lon directions
      REAL(KIND=wp) :: &
         & za1, za2, za, zc, zd
      
      plamscl_max = -1.0_wp
      pphiscl_max = -1.0_wp

      DO ji = 1, ki
         DO jj = 1, kj

            !Calculate distance in metres in zonal(x) direction

            zlon1 = rad * ( pglam(ji,jj) + ( 0.5_wp * plamscl_deg ) )
            zlon2 = rad * ( pglam(ji,jj) - ( 0.5_wp * plamscl_deg ) )
            zlat1 = rad * pgphi(ji,jj)
            zlat2 = rad * pgphi(ji,jj)
            zdlon = zlon2 - zlon1
            zdlat = zlat2 - zlat1

            za1 = sin( zdlat/2.0_wp )
            za2 = sin( zdlon/2.0_wp )
            za = ( za1 * za1 ) + ( COS( zlat1 ) * COS( zlat2 ) * ( za2 * za2 ) )
            zc = 2.0_wp * atan2( SQRT( za ), SQRT( 1.0_wp-za ) )
            zd = ra * zc

            IF ( zd > plamscl_max ) plamscl_max = zd

            !Calculate distance in metres in meridional(y) direction

            zlon1 = rad * pglam(ji,jj)
            zlon2 = rad * pglam(ji,jj)
            zlat1 = rad * ( pgphi(ji,jj) + ( 0.5_wp * pphiscl_deg ) )
            zlat2 = rad * ( pgphi(ji,jj) - ( 0.5_wp * pphiscl_deg ) )
            zdlon = zlon2 - zlon1
            zdlat = zlat2 - zlat1

            za1 = sin( zdlat/2.0_wp )
            za2 = sin( zdlon/2.0_wp )
            za = ( za1 * za1 ) + ( COS( zlat1 ) * COS( zlat2 ) * ( za2 * za2 ) )
            zc = 2.0_wp * atan2( SQRT( za ), SQRT( 1.0_wp-za ) )
            zd = ra * zc

            IF ( zd > pphiscl_max ) pphiscl_max = zd

         END DO
      END DO
         
   END SUBROUTINE obs_deg2dist

   SUBROUTINE obs_dist2corners(pglam_bl, pglam_br, pglam_tl, pglam_tr, &
      &                        pgphi_bl, pgphi_br, pgphi_tl, pgphi_tr, &
      &                        plam, pphi, pxvert, pyvert, pdist)
      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_dist2corners  ***
      !!
      !! ** Purpose : Calculate distance from centre of obs footprint to the corners of a grid box
      !!
      !! ** Method  : Use great circle distance formulae.
      !!              Order of corners is 1=topleft, 2=topright, 3=bottomright, 4=bottomleft
      !!
      !! ** Action  :
      !!                   
      !! References : 
      !!
      !! History :
      !!        ! 14/01. M. Martin.
      !!-----------------------------------------------------------------------
      !! * Modules used
      !! * Arguments
         REAL(KIND=wp), INTENT(IN) :: &
            &  pglam_bl, pglam_br, &       !lon at corners of grid box
            &  pglam_tl, pglam_tr
         REAL(KIND=wp), INTENT(IN) :: &
            &  pgphi_bl, pgphi_br, &       !lat at corners of grid box
            &  pgphi_tl, pgphi_tr
         REAL(KIND=wp), INTENT(IN) :: &
            &  pphi, plam                  !lat/lon of centre of obs footprint
         REAL(KIND=wp), DIMENSION(4), INTENT(OUT) :: &
            &  pxvert, pyvert              !x/y location (in metres relative to centre of obs footprint) of corners
         REAL(KIND=wp), DIMENSION(4), INTENT(OUT) :: &
            &  pdist                       !distance (in metres) of each corner relative to centre of obs footprint

      !! * Local variables
         INTEGER :: &
            &  jvert                       !Counter for corners
         REAL(KIND=wp) :: &
            &  zphi, zlam                  !Local values for lon/lat of corners
         REAL(KIND=wp) :: & 
            &  za1, za2,  &                !For great circle distance calculations
            &  zb1, zb2,  &
            &  zc1, zc2
         REAL(KIND=wp) :: &
            &  zdist_centre_lat, &         !Distances in lat/lon directions (in metres)
            &  zdist_centre_lon

      !!-----------------------------------------------------------------------

         ! Work out latitudinal and longitudinal distance from centre of 
         ! obs fp to corners of grid box
         DO jvert = 1, 4
            SELECT CASE(jvert)
            CASE(1)
               zphi = pgphi_tl
               zlam = pglam_tl
            CASE(2)
               zphi = pgphi_tr
               zlam = pglam_tr
            CASE(3)
               zphi = pgphi_br
               zlam = pglam_br
            CASE(4)
               zphi = pgphi_bl
               zlam = pglam_bl
            END SELECT
            
            IF (zlam == plam ) THEN
               pxvert(jvert) = 0.0_wp
            ELSE
               za1 = SIN( zphi * rad )
               za2 = SIN( zphi * rad )
               zb1 = COS( zphi * rad ) * COS( zlam * rad )
               zb2 = COS( zphi * rad ) * COS( plam  * rad )
               zc1 = COS( zphi * rad ) * SIN( zlam * rad )
               zc2 = COS( zphi * rad ) * SIN( plam  * rad )
               pxvert(jvert) = grt_cir_dis( za1, za2, zb1, zb2, zc1, zc2 )
               pxvert(jvert) =  ra * pxvert(jvert)
               IF ( zlam < plam ) pxvert(jvert) = - pxvert(jvert)
            ENDIF
            
            IF ( zphi == pphi ) THEN
               pyvert(jvert) = 0.0_wp
            ELSE
               za1 = SIN( zphi * rad )
               za2 = SIN( pphi * rad )
               zb1 = COS( zphi * rad ) * COS( zlam * rad )
               zb2 = COS( pphi * rad ) * COS( zlam * rad )
               zc1 = COS( zphi * rad ) * SIN( zlam * rad )
               zc2 = COS( pphi * rad ) * SIN( zlam * rad )
               pyvert(jvert) = grt_cir_dis( za1, za2, zb1, zb2, zc1, zc2 )
               pyvert(jvert) =  ra * pyvert(jvert)
               IF ( zphi < pphi ) pyvert(jvert) = - pyvert(jvert)
            ENDIF

            !Calculate the distance of each vertex relative to centre of obs footprint
            pdist(jvert) = SQRT( ( pxvert(jvert) * pxvert(jvert) ) + &
               &                 ( pyvert(jvert) * pyvert(jvert) ) )

         END DO

      END SUBROUTINE obs_dist2corners

END MODULE obs_averg_h2d
