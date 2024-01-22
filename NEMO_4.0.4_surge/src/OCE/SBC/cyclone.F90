MODULE cyclone
   !!======================================================================
   !!                       ***  MODULE  cyclone  ***
   !! add the Tropical Cyclones along tracks to the surface wind forcing
   !!                 
   !!======================================================================
   !! History : 3.3  ! 2010-05  (E Vincent, G Madec, S Masson)  Original code
   !!----------------------------------------------------------------------

#if defined key_cyclone
   !!----------------------------------------------------------------------
   !!  'key_cyclone' : key option add Tropical Cyclones in the wind forcing
   !!----------------------------------------------------------------------
   !!   wnd_cyc      : 1 module subroutine
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE sbc_oce         ! surface boundary condition: ocean
   USE dom_oce         ! ocean space domain variables
   USE phycst          ! physical constant
   USE fldread         ! read input fields
   USE in_out_manager  ! I/O manager
   USE geo2ocean       ! tools for projection on ORCA grid
   USE lib_mpp       

   IMPLICIT NONE
   PRIVATE

   PUBLIC   wnd_cyc   ! routine called in sbcblk.F90 module

   INTEGER , PARAMETER ::   jp_is1  = 1   ! index of presence 1 or absence 0 of a TC record
   INTEGER , PARAMETER ::   jp_lon  = 2   ! index of longitude for present TCs
   INTEGER , PARAMETER ::   jp_lat  = 3   ! index of latitude  for present TCs
   INTEGER , PARAMETER ::   jp_vmax = 4   ! index of max wind  for present TCs
   INTEGER , PARAMETER ::   jp_pres = 5   ! index of eye-pres  for present TCs

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf   ! structure of input fields (file informations, fields read)

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE wnd_cyc( kt, pwnd_i, pwnd_j )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE wnd_cyc  ***
      !!
      !! ** Purpose :  Add cyclone winds on the ORCA grid
      !!
      !! ** Action  : - open TC data, find TCs for the current timestep
      !!              - for each potential TC, add the winds on the grid
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in)                      ::   kt       ! time step index 
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   pwnd_i   ! wind speed i-components at T-point ORCA direction
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::   pwnd_j   ! wind speed j-components at T-point ORCA direction
      ! 
      !!
      INTEGER  ::   ji, jj , jtc        ! loop arguments
      INTEGER  ::   ierror              ! loop arguments
      INTEGER  ::   vortex=1            ! vortex shape to be used: 0=Holland 1=Willoughby
      REAL(wp) ::   zrout1=1.5e6        ! distance from center where we begin to kill vortex (m)
      REAL(wp) ::   zrout2=2.5e6        ! distance from center where we bring vortex to zero (m)
      REAL(wp) ::   zb                  ! power in Holland vortex shape
      REAL(wp) ::   zA                  ! shape parameter in Willoughby vortex : A transtion between first and second outter exp
      REAL(wp) ::   zn                  ! shape parameter in Willoughby vortex : n power law in the eye
      REAL(wp) ::   zXX1                ! shape parameter in Willoughby vortex : decay length second outter exponential
      REAL(wp) ::   zXX2                ! shape parameter in Willoughby vortex : decay length first  outter exponential
      REAL(wp) ::   zztmp               ! temporary 
      REAL(wp) ::   zzrglam, zzrgphi    ! temporary 
      REAL(wp) ::   ztheta              ! azimuthal angle
      REAL(wp) ::   zdist               ! dist to the TC center
      REAL(wp) ::   zhemi               ! 1 for NH ;  -1 for SH
      REAL(wp) ::   zinfl               ! clim inflow angle in TCs
      REAL(wp) ::   zrmw                ! mean radius of Max wind of a tropical cyclone (Willoughby 2004) [m]
      REAL(wp) ::   zwnd_r, zwnd_t      ! radial and tangential components of the wind
      REAL(wp) ::   zvmax               ! timestep interpolated vmax
      REAL(wp) ::   zrlon, zrlat        ! temporary 
      REAL(wp), DIMENSION(jpi,jpj) ::   zwnd_x, zwnd_y   ! zonal and meridional components of the wind
      REAL(wp), DIMENSION(14,5)    ::   ztct                ! tropical cyclone track data at kt
      !
      CHARACTER(len=100) ::  cn_dir            ! Root directory for location of files
      TYPE(FLD_N), DIMENSION(1) ::   slf_i     ! array of namelist informations on the TC position
      TYPE(FLD_N) ::   sn_tc                   ! informations about the fields to be read
      !!--------------------------------------------------------------------

      !                                         ! ====================== !
      IF( kt == nit000 ) THEN                   !  First call kt=nit000  !
         !                                      ! ====================== !
         ! set file information (default values)
         cn_dir = './'       ! directory in which the model is executed
         !
         ! (NB: frequency positive => hours, negative => months)
         !          !    file     ! frequency !  variable ! time intep !  clim   ! 'yearly' or ! weights  ! rotation   ! land/sea mask !
         !          !    name     !  (hours)  !   name    !   (T/F)    !  (T/F)  !  'monthly'  ! filename ! pairs      ! filename      !
         sn_tc = FLD_N( 'tc_track',     6     ,  'tc'     ,  .true.    , .false. ,   'yearly'  , ''       , ''         , ''            )
         !
         !  Namelist is read in namsbc_blk
         ! set sf structure
         ALLOCATE( sf(1), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'wnd_cyc: unable to allocate sf structure' )   ;   RETURN
         ENDIF
         ALLOCATE( sf(1)%fnow(14,5,1) )
         ALLOCATE( sf(1)%fdta(14,5,1,2) )
         slf_i(1) = sn_tc
         !
         ! fill sf with slf_i and control print
         CALL fld_fill( sf, slf_i, cn_dir, 'sbc_blk_tc', 'tropical cyclone track', 'namsbc_tc' )
         !
      ENDIF


      !       Interpolation of lon lat vmax... at the current timestep
      !       ***************************************************************

      CALL fld_read( kt, nn_fsbc, sf )                   ! input fields provided at the current time-step

      ztct(:,:) = sf(1)%fnow(:,:,1)

      !       Add TC wind on the grid
      !       ***************************************************************

      zwnd_x(:,:) = 0.e0  
      zwnd_y(:,:) = 0.e0 
      
      DO jtc = 1, 14
         !
         IF( ztct(jtc,jp_is1) == 1 ) THEN                ! cyclone is defined in this slot ? yes--> begin

            zvmax =       ztct(jtc,jp_vmax)
            zrlon = rad * ztct(jtc,jp_lon )
            zrlat = rad * ztct(jtc,jp_lat )
            zhemi = SIGN( 1. , zrlat )
            zinfl = 15.* rad                             ! clim inflow angle in Tropical Cyclones
         IF ( vortex == 0 ) THEN

            ! Vortex Holland reconstruct wind at each lon-lat position
            ! ********************************************************
            zrmw = 51.6 * EXP( -0.0223*zvmax + 0.0281* ABS( ztct(jtc,jp_lat) ) ) * 1000.
            ! climatological ZRMW of cyclones as a function of wind and latitude (Willoughby 2004)             
            ! zb = 1.0036 + 0.0173 * zvmax - 0.0313 * LOG(zrmw/1000.) + 0.0087 * ABS( ztct(jtc,jp_lat) ) 
            ! fitted B parameter (Willoughby 2004)
            zb = 2.

            DO jj = 1, jpj
               DO ji = 1, jpi

                  ! calc distance between TC center and any point following great circle
                  ! source : http://www.movable-type.co.uk/scripts/latlong.html
                  zzrglam = rad * glamt(ji,jj) - zrlon
                  zzrgphi = rad * gphit(ji,jj)
                  zdist = ra * ACOS(  SIN( zrlat ) * SIN( zzrgphi )   &
                     &              + COS( zrlat ) * COS( zzrgphi ) * COS( zzrglam ) )

                 IF (zdist < zrout2) THEN ! calculation of wind only to a given max radius
                  ! shape of the wind profile
                  zztmp = ( zrmw / ( zdist + 1.e-12 ) )**zb
                  zztmp =  zvmax * SQRT( zztmp * EXP(1. - zztmp) )    

                  IF (zdist > zrout1) THEN ! bring to zero between r_out1 and r_out2
                     zztmp = zztmp * ( (zrout2-zdist)*1.e-6 )
                  ENDIF

                  ! !!! KILL EQ WINDS
                  ! IF (SIGN( 1. , zrlat ) /= zhemi) THEN
                  !    zztmp = 0.                              ! winds in other hemisphere
                  !    IF (ABS(gphit(ji,jj)) <= 5.) zztmp=0.   ! kill between 5N-5S
                  ! ENDIF
                  ! IF (ABS(gphit(ji,jj)) <= 10. .and. ABS(gphit(ji,jj)) > 5.) THEN
                  !    zztmp = zztmp * ( 1./5. * (ABS(gphit(ji,jj)) - 5.) ) 
                  !    !linear to zero between 10 and 5
                  ! ENDIF
                  ! !!! / KILL EQ

                  IF (ABS(gphit(ji,jj)) >= 55.) zztmp = 0. ! kill weak spurious winds at high latitude

                  zwnd_t =   COS( zinfl ) * zztmp    
                  zwnd_r = - SIN( zinfl ) * zztmp

                  ! Project radial-tangential components on zonal-meridional components
                  ! -------------------------------------------------------------------
                  
                  ! ztheta = azimuthal angle of the great circle between two points
                  zztmp = COS( zrlat ) * SIN( zzrgphi ) &
                     &  - SIN( zrlat ) * COS( zzrgphi ) * COS( zzrglam )
                  ztheta = ATAN2(        COS( zzrgphi ) * SIN( zzrglam ) , zztmp )

                  zwnd_x(ji,jj) = zwnd_x(ji,jj) - zhemi * COS(ztheta)*zwnd_t + SIN(ztheta)*zwnd_r
                  zwnd_y(ji,jj) = zwnd_y(ji,jj) + zhemi * SIN(ztheta)*zwnd_t + COS(ztheta)*zwnd_r
                 ENDIF
               END DO
            END DO
         
         ELSE IF ( vortex == 1 ) THEN

            ! Vortex Willoughby reconstruct wind at each lon-lat position
            ! ***********************************************************
            zrmw = 46.4 * EXP( -0.0155*zvmax + 0.0169* ABS( ztct(jtc,jp_lat) ) )*1000.
            ! climatological ZRMW of cyclones as a function of wind and latitude (Willoughby 2006)
            zXX2 = 25.*1000.                                              ! 25km fixed "near-eye" exponential decay
            zXX1 = ( 287.6  - 1.942 *zvmax + 7.799 *LOG(zrmw/1000.) + 1.819 *ABS( ztct(jtc,jp_lat) ) )*1000.    
            zn   =   2.1340 + 0.0077*zvmax - 0.4522*LOG(zrmw/1000.) - 0.0038*ABS( ztct(jtc,jp_lat) )            
            zA   =   0.5913 + 0.0029*zvmax - 0.1361*LOG(zrmw/1000.) - 0.0042*ABS( ztct(jtc,jp_lat) )  
            IF (zA < 0) THEN 
               zA=0
            ENDIF           
        
            DO jj = 1, jpj
               DO ji = 1, jpi
                                  
                  zzrglam = rad * glamt(ji,jj) - zrlon
                  zzrgphi = rad * gphit(ji,jj)
                  zdist = ra * ACOS(  SIN( zrlat ) * SIN( zzrgphi )   &
                     &              + COS( zrlat ) * COS( zzrgphi ) * COS( zzrglam ) )

                 IF (zdist < zrout2) THEN ! calculation of wind only to a given max radius
               
                  ! shape of the wind profile                     
                  IF (zdist <= zrmw) THEN     ! inside the Radius of Maximum Wind
                     zztmp  = zvmax * (zdist/zrmw)**zn
                  ELSE 
                     zztmp  = zvmax * ( (1-zA) * EXP(- (zdist-zrmw)/zXX1 ) + zA * EXP(- (zdist-zrmw)/zXX2 ) )
                  ENDIF

                  IF (zdist > zrout1) THEN ! bring to zero between r_out1 and r_out2
                     zztmp = zztmp * ( (zrout2-zdist)*1.e-6 )
                  ENDIF

                  ! !!! KILL EQ WINDS
                  ! IF (SIGN( 1. , zrlat ) /= zhemi) THEN
                  !    zztmp = 0.                              ! winds in other hemisphere
                  !    IF (ABS(gphit(ji,jj)) <= 5.) zztmp=0.   ! kill between 5N-5S
                  ! ENDIF
                  ! IF (ABS(gphit(ji,jj)) <= 10. .and. ABS(gphit(ji,jj)) > 5.) THEN
                  !    zztmp = zztmp * ( 1./5. * (ABS(gphit(ji,jj)) - 5.) ) 
                  !    !linear to zero between 10 and 5
                  ! ENDIF
                  ! !!! / KILL EQ

                  IF (ABS(gphit(ji,jj)) >= 55.) zztmp = 0. ! kill weak spurious winds at high latitude

                  zwnd_t =   COS( zinfl ) * zztmp    
                  zwnd_r = - SIN( zinfl ) * zztmp

                  ! Project radial-tangential components on zonal-meridional components
                  ! -------------------------------------------------------------------
                  
                  ! ztheta = azimuthal angle of the great circle between two points
                  zztmp = COS( zrlat ) * SIN( zzrgphi ) &
                     &  - SIN( zrlat ) * COS( zzrgphi ) * COS( zzrglam )
                  ztheta = ATAN2(        COS( zzrgphi ) * SIN( zzrglam ) , zztmp )

                  zwnd_x(ji,jj) = zwnd_x(ji,jj) - zhemi * COS(ztheta)*zwnd_t + SIN(ztheta)*zwnd_r
                  zwnd_y(ji,jj) = zwnd_y(ji,jj) + zhemi * SIN(ztheta)*zwnd_t + COS(ztheta)*zwnd_r
                  
                 ENDIF
               END DO
            END DO
         ENDIF                                         ! / vortex Holland or Wiloughby
         ENDIF                                           ! / cyclone is defined in this slot ? yes--> begin
      END DO ! / end simultaneous cyclones loop

      CALL rot_rep ( zwnd_x, zwnd_y, 'T', 'en->i', pwnd_i ) !rotation of components on ORCA grid
      CALL rot_rep ( zwnd_x, zwnd_y, 'T', 'en->j', pwnd_j ) !rotation of components on ORCA grid

   END SUBROUTINE wnd_cyc

#endif

   !!======================================================================
END MODULE cyclone
