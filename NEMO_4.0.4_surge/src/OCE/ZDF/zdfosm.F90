MODULE zdfosm
   !!======================================================================
   !!                       ***  MODULE  zdfosm  ***
   !! Ocean physics:  vertical mixing coefficient compute from the OSMOSIS
   !!                 turbulent closure parameterization
   !!=====================================================================
   !!  History : NEMO 4.0  ! A. Grant, G. Nurser
   !! 15/03/2017  Changed calculation of pycnocline thickness in unstable conditions and stable conditions AG
   !! 15/03/2017  Calculation of pycnocline gradients for stable conditions changed. Pycnocline gradients now depend on stability of the OSBL. A.G
   !! 06/06/2017  (1) Checks on sign of buoyancy jump in calculation of  OSBL depth.  A.G.
   !!             (2) Removed variable zbrad0, zbradh and zbradav since they are not used.
   !!             (3) Approximate treatment for shear turbulence.
   !!                        Minimum values for zustar and zustke.
   !!                        Add velocity scale, zvstr, that tends to zustar for large Langmuir numbers.
   !!                        Limit maximum value for Langmuir number.
   !!                        Use zvstr in definition of stability parameter zhol.
   !!             (4) Modified parametrization of entrainment flux, changing original coefficient 0.0485 for Langmuir contribution to 0.135 * zla
   !!             (5) For stable boundary layer add factor that depends on length of timestep to 'slow' collapse and growth. Make sure buoyancy jump not negative.
   !!             (6) For unstable conditions when growth is over multiple levels, limit change to maximum of one level per cycle through loop.
   !!             (7) Change lower limits for loops that calculate OSBL averages from 1 to 2. Large gradients between levels 1 and 2 can cause problems.
   !!             (8) Change upper limits from ibld-1 to ibld.
   !!             (9) Calculation of pycnocline thickness in unstable conditions. Check added to ensure that buoyancy jump is positive before calculating Ri.
   !!            (10) Thickness of interface layer at base of the stable OSBL set by Richardson number. Gives continuity in transition from unstable OSBL.
   !!            (11) Checks that buoyancy jump is poitive when calculating pycnocline profiles.
   !!            (12) Replace zwstrl with zvstr in calculation of eddy viscosity.
   !! 27/09/2017 (13) Calculate Stokes drift and Stokes penetration depth from wave information
   !!            (14) Bouyancy flux due to entrainment changed to include contribution from shear turbulence (for testing commented out).
   !! 28/09/2017 (15) Calculation of Stokes drift moved into separate do-loops to allow for different options for the determining the Stokes drift to be added.
   !!            (16) Calculation of Stokes drift from windspeed for PM spectrum (for testing, commented out)
   !!            (17) Modification to Langmuir velocity scale to include effects due to the Stokes penetration depth (for testing, commented out)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'ln_zdfosm'                                             OSMOSIS scheme
   !!----------------------------------------------------------------------
   !!   zdf_osm       : update momentum and tracer Kz from osm scheme
   !!   zdf_osm_init  : initialization, namelist read, and parameters control
   !!   osm_rst       : read (or initialize) and write osmosis restart fields
   !!   tra_osm       : compute and add to the T & S trend the non-local flux
   !!   trc_osm       : compute and add to the passive tracer trend the non-local flux (TBD)
   !!   dyn_osm       : compute and add to u & v trensd the non-local flux
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
                      ! uses wn from previous time step (which is now wb) to calculate hbl
   USE dom_oce        ! ocean space and time domain
   USE zdf_oce        ! ocean vertical physics
   USE sbc_oce        ! surface boundary condition: ocean
   USE sbcwave        ! surface wave parameters
   USE phycst         ! physical constants
   USE eosbn2         ! equation of state
   USE traqsr         ! details of solar radiation absorption
   USE zdfddm         ! double diffusion mixing (avs array)
   USE iom            ! I/O library
   USE lib_mpp        ! MPP library
   USE trd_oce        ! ocean trends definition
   USE trdtra         ! tracers trends
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_osm       ! routine called by step.F90
   PUBLIC   zdf_osm_init  ! routine called by nemogcm.F90
   PUBLIC   osm_rst       ! routine called by step.F90
   PUBLIC   tra_osm       ! routine called by step.F90
   PUBLIC   trc_osm       ! routine called by trcstp.F90
   PUBLIC   dyn_osm       ! routine called by 'step.F90'

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ghamu    !: non-local u-momentum flux
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ghamv    !: non-local v-momentum flux
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ghamt    !: non-local temperature flux (gamma/<ws>o)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ghams    !: non-local salinity flux (gamma/<ws>o)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   etmean   !: averaging operator for avt
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hbl      !: boundary layer depth
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   hbli     !: intial boundary layer depth for stable blayer
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   dstokes  !: penetration depth of the Stokes drift.

   !                      !!** Namelist  namzdf_osm  **
   LOGICAL  ::   ln_use_osm_la      ! Use namelist  rn_osm_la
   REAL(wp) ::   rn_osm_la          ! Turbulent Langmuir number
   REAL(wp) ::   rn_osm_dstokes     ! Depth scale of Stokes drift
   REAL(wp) ::   rn_osm_hbl0 = 10._wp       ! Initial value of hbl for 1D runs
   INTEGER  ::   nn_ave             ! = 0/1 flag for horizontal average on avt
   INTEGER  ::   nn_osm_wave = 0    ! = 0/1/2 flag for getting stokes drift from La# / PM wind-waves/Inputs into sbcwave
   LOGICAL  ::   ln_dia_osm         ! Use namelist  rn_osm_la


   LOGICAL  ::   ln_kpprimix  = .true.  ! Shear instability mixing
   REAL(wp) ::   rn_riinfty   = 0.7     ! local Richardson Number limit for shear instability
   REAL(wp) ::   rn_difri    =  0.005   ! maximum shear mixing at Rig = 0    (m2/s)
   LOGICAL  ::   ln_convmix  = .true.   ! Convective instability mixing
   REAL(wp) ::   rn_difconv = 1._wp     ! diffusivity when unstable below BL  (m2/s)

   !                                    !!! ** General constants  **
   REAL(wp) ::   epsln   = 1.0e-20_wp   ! a small positive number
   REAL(wp) ::   pthird  = 1._wp/3._wp  ! 1/3
   REAL(wp) ::   p2third = 2._wp/3._wp  ! 2/3

   INTEGER :: idebug = 236
   INTEGER :: jdebug = 228
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_osm_alloc()
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION zdf_osm_alloc  ***
      !!----------------------------------------------------------------------
     ALLOCATE( ghamu(jpi,jpj,jpk), ghamv(jpi,jpj,jpk), ghamt(jpi,jpj,jpk), ghams(jpi,jpj,jpk), &
          &      hbl(jpi,jpj)    ,  hbli(jpi,jpj)    , dstokes(jpi, jpj) ,                     &
          &   etmean(jpi,jpj,jpk),  STAT= zdf_osm_alloc )
     IF( zdf_osm_alloc /= 0 )   CALL ctl_warn('zdf_osm_alloc: failed to allocate zdf_osm arrays')
     CALL mpp_sum ( 'zdfosm', zdf_osm_alloc )
   END FUNCTION zdf_osm_alloc


   SUBROUTINE zdf_osm( kt, p_avm, p_avt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zdf_osm  ***
      !!
      !! ** Purpose :   Compute the vertical eddy viscosity and diffusivity
      !!      coefficients and non local mixing using the OSMOSIS scheme
      !!
      !! ** Method :   The boundary layer depth hosm is diagnosed at tracer points
      !!      from profiles of buoyancy, and shear, and the surface forcing.
      !!      Above hbl (sigma=-z/hbl <1) the mixing coefficients are computed from
      !!
      !!                      Kx =  hosm  Wx(sigma) G(sigma)
      !!
      !!             and the non local term ghamt = Cs / Ws(sigma) / hosm
      !!      Below hosm  the coefficients are the sum of mixing due to internal waves
      !!      shear instability and double diffusion.
      !!
      !!      -1- Compute the now interior vertical mixing coefficients at all depths.
      !!      -2- Diagnose the boundary layer depth.
      !!      -3- Compute the now boundary layer vertical mixing coefficients.
      !!      -4- Compute the now vertical eddy vicosity and diffusivity.
      !!      -5- Smoothing
      !!
      !!        N.B. The computation is done from jk=2 to jpkm1
      !!             Surface value of avt are set once a time to zero
      !!             in routine zdf_osm_init.
      !!
      !! ** Action  :   update the non-local terms ghamts
      !!                update avt (before vertical eddy coef.)
      !!
      !! References : Large W.G., Mc Williams J.C. and Doney S.C.
      !!         Reviews of Geophysics, 32, 4, November 1994
      !!         Comments in the code refer to this paper, particularly
      !!         the equation number. (LMD94, here after)
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   kt            ! ocean time step
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::  p_avm, p_avt   ! momentum and tracer Kz (w-points)
      !!
      INTEGER ::   ji, jj, jk                   ! dummy loop indices
      INTEGER ::   ikbot, jkmax, jkm1, jkp2     !

      REAL(wp) ::   ztx, zty, zflageos, zstabl, zbuofdep,zucube      !
      REAL(wp) ::   zbeta, zthermal                                  !
      REAL(wp) ::   zehat, zeta, zhrib, zsig, zscale, zwst, zws, zwm ! Velocity scales
      REAL(wp) ::   zwsun, zwmun, zcons, zconm, zwcons, zwconm       !
      REAL(wp) ::   zsr, zbw, ze, zb, zd, zc, zaw, za, zb1, za1, zkw, zk0, zcomp , zrhd,zrhdr,zbvzed   ! In situ density
      INTEGER  ::   jm                          ! dummy loop indices
      REAL(wp) ::   zr1, zr2, zr3, zr4, zrhop   ! Compression terms
      REAL(wp) ::   zflag, zrn2, zdep21, zdep32, zdep43
      REAL(wp) ::   zesh2, zri, zfri            ! Interior richardson mixing
      REAL(wp) ::   zdelta, zdelta2, zdzup, zdzdn, zdzh, zvath, zgat1, zdat1, zkm1m, zkm1t
      REAL(wp) :: zt,zs,zu,zv,zrh               ! variables used in constructing averages
! Scales
      REAL(wp), DIMENSION(jpi,jpj) :: zrad0     ! Surface solar temperature flux (deg m/s)
      REAL(wp), DIMENSION(jpi,jpj) :: zradh     ! Radiative flux at bl base (Buoyancy units)
      REAL(wp), DIMENSION(jpi,jpj) :: zradav    ! Radiative flux, bl average (Buoyancy Units)
      REAL(wp), DIMENSION(jpi,jpj) :: zustar    ! friction velocity
      REAL(wp), DIMENSION(jpi,jpj) :: zwstrl    ! Langmuir velocity scale
      REAL(wp), DIMENSION(jpi,jpj) :: zvstr     ! Velocity scale that ends to zustar for large Langmuir number.
      REAL(wp), DIMENSION(jpi,jpj) :: zwstrc    ! Convective velocity scale
      REAL(wp), DIMENSION(jpi,jpj) :: zuw0      ! Surface u-momentum flux
      REAL(wp), DIMENSION(jpi,jpj) :: zvw0      ! Surface v-momentum flux
      REAL(wp), DIMENSION(jpi,jpj) :: zwth0     ! Surface heat flux (Kinematic)
      REAL(wp), DIMENSION(jpi,jpj) :: zws0      ! Surface freshwater flux
      REAL(wp), DIMENSION(jpi,jpj) :: zwb0      ! Surface buoyancy flux
      REAL(wp), DIMENSION(jpi,jpj) :: zwthav    ! Heat flux - bl average
      REAL(wp), DIMENSION(jpi,jpj) :: zwsav     ! freshwater flux - bl average
      REAL(wp), DIMENSION(jpi,jpj) :: zwbav     ! Buoyancy flux - bl average
      REAL(wp), DIMENSION(jpi,jpj) :: zwb_ent   ! Buoyancy entrainment flux
      REAL(wp), DIMENSION(jpi,jpj) :: zustke    ! Surface Stokes drift
      REAL(wp), DIMENSION(jpi,jpj) :: zla       ! Trubulent Langmuir number
      REAL(wp), DIMENSION(jpi,jpj) :: zcos_wind ! Cos angle of surface stress
      REAL(wp), DIMENSION(jpi,jpj) :: zsin_wind ! Sin angle of surface stress
      REAL(wp), DIMENSION(jpi,jpj) :: zhol      ! Stability parameter for boundary layer
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: lconv ! unstable/stable bl

      ! mixed-layer variables

      INTEGER, DIMENSION(jpi,jpj) :: ibld ! level of boundary layer base
      INTEGER, DIMENSION(jpi,jpj) :: imld ! level of mixed-layer depth (pycnocline top)

      REAL(wp) :: ztgrad,zsgrad,zbgrad ! Temporary variables used to calculate pycnocline gradients
      REAL(wp) :: zugrad,zvgrad        ! temporary variables for calculating pycnocline shear

      REAL(wp), DIMENSION(jpi,jpj) :: zhbl  ! bl depth - grid
      REAL(wp), DIMENSION(jpi,jpj) :: zhml  ! ml depth - grid
      REAL(wp), DIMENSION(jpi,jpj) :: zdh   ! pycnocline depth - grid
      REAL(wp), DIMENSION(jpi,jpj) :: zdhdt ! BL depth tendency
      REAL(wp), DIMENSION(jpi,jpj) :: zt_bl,zs_bl,zu_bl,zv_bl,zrh_bl  ! averages over the depth of the blayer
      REAL(wp), DIMENSION(jpi,jpj) :: zt_ml,zs_ml,zu_ml,zv_ml,zrh_ml  ! averages over the depth of the mixed layer
      REAL(wp), DIMENSION(jpi,jpj) :: zdt_bl,zds_bl,zdu_bl,zdv_bl,zdrh_bl,zdb_bl ! difference between blayer average and parameter at base of blayer
      REAL(wp), DIMENSION(jpi,jpj) :: zdt_ml,zds_ml,zdu_ml,zdv_ml,zdrh_ml,zdb_ml ! difference between mixed layer average and parameter at base of blayer
      REAL(wp), DIMENSION(jpi,jpj) :: zwth_ent,zws_ent ! heat and salinity fluxes at the top of the pycnocline
      REAL(wp), DIMENSION(jpi,jpj) :: zuw_bse,zvw_bse  ! momentum fluxes at the top of the pycnocline
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zdtdz_pyc    ! parametrized gradient of temperature in pycnocline
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zdsdz_pyc    ! parametrised gradient of salinity in pycnocline
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zdbdz_pyc    ! parametrised gradient of buoyancy in the pycnocline
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zdudz_pyc    ! u-shear across the pycnocline
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zdvdz_pyc    ! v-shear across the pycnocline

      ! Flux-gradient relationship variables

      REAL(wp) :: zl_c,zl_l,zl_eps  ! Used to calculate turbulence length scale.

      REAL(wp), DIMENSION(jpi,jpj) :: zdifml_sc,zvisml_sc,zdifpyc_sc,zvispyc_sc,zbeta_d_sc,zbeta_v_sc ! Scales for eddy diffusivity/viscosity
      REAL(wp), DIMENSION(jpi,jpj) :: zsc_wth_1,zsc_ws_1 ! Temporary scales used to calculate scalar non-gradient terms.
      REAL(wp), DIMENSION(jpi,jpj) :: zsc_uw_1,zsc_uw_2,zsc_vw_1,zsc_vw_2 ! Temporary scales for non-gradient momentum flux terms.
      REAL(wp), DIMENSION(jpi,jpj) :: zhbl_t ! holds boundary layer depth updated by full timestep

      ! For calculating Ri#-dependent mixing
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: z3du   ! u-shear^2
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: z3dv   ! v-shear^2
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zrimix ! spatial form of ri#-induced diffusion

      ! Temporary variables
      INTEGER :: inhml
      INTEGER :: i_lconv_alloc
      REAL(wp) :: znd,znd_d,zznd_ml,zznd_pyc,zznd_d ! temporary non-dimensional depths used in various routines
      REAL(wp) :: ztemp, zari, zpert, zzdhdt, zdb   ! temporary variables
      REAL(wp) :: zthick, zz0, zz1 ! temporary variables
      REAL(wp) :: zvel_max, zhbl_s ! temporary variables
      REAL(wp) :: zfac             ! temporary variable
      REAL(wp) :: zus_x, zus_y     ! temporary Stokes drift
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zviscos ! viscosity
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zdiffut ! t-diffusivity

      ! For debugging
      INTEGER :: ikt
      !!--------------------------------------------------------------------
      !
      ALLOCATE( lconv(jpi,jpj),  STAT= i_lconv_alloc )
      IF( i_lconv_alloc /= 0 )   CALL ctl_warn('zdf_osm: failed to allocate lconv')

      ibld(:,:)   = 0     ; imld(:,:)  = 0
      zrad0(:,:)  = 0._wp ; zradh(:,:) = 0._wp ; zradav(:,:)    = 0._wp ; zustar(:,:)    = 0._wp
      zwstrl(:,:) = 0._wp ; zvstr(:,:) = 0._wp ; zwstrc(:,:)    = 0._wp ; zuw0(:,:)      = 0._wp
      zvw0(:,:)   = 0._wp ; zwth0(:,:) = 0._wp ; zws0(:,:)      = 0._wp ; zwb0(:,:)      = 0._wp
      zwthav(:,:) = 0._wp ; zwsav(:,:) = 0._wp ; zwbav(:,:)     = 0._wp ; zwb_ent(:,:)   = 0._wp
      zustke(:,:) = 0._wp ; zla(:,:)   = 0._wp ; zcos_wind(:,:) = 0._wp ; zsin_wind(:,:) = 0._wp
      zhol(:,:)   = 0._wp
      lconv(:,:)  = .FALSE.
      ! mixed layer
      ! no initialization of zhbl or zhml (or zdh?)
      zhbl(:,:)    = 1._wp ; zhml(:,:)    = 1._wp ; zdh(:,:)      = 1._wp ; zdhdt(:,:)   = 0._wp
      zt_bl(:,:)   = 0._wp ; zs_bl(:,:)   = 0._wp ; zu_bl(:,:)    = 0._wp ; zv_bl(:,:)   = 0._wp
      zrh_bl(:,:)  = 0._wp ; zt_ml(:,:)   = 0._wp ; zs_ml(:,:)    = 0._wp ; zu_ml(:,:)   = 0._wp
      zv_ml(:,:)   = 0._wp ; zrh_ml(:,:)  = 0._wp ; zdt_bl(:,:)   = 0._wp ; zds_bl(:,:)  = 0._wp
      zdu_bl(:,:)  = 0._wp ; zdv_bl(:,:)  = 0._wp ; zdrh_bl(:,:)  = 0._wp ; zdb_bl(:,:)  = 0._wp
      zdt_ml(:,:)  = 0._wp ; zds_ml(:,:)  = 0._wp ; zdu_ml(:,:)   = 0._wp ; zdv_ml(:,:)  = 0._wp
      zdrh_ml(:,:) = 0._wp ; zdb_ml(:,:)  = 0._wp ; zwth_ent(:,:) = 0._wp ; zws_ent(:,:) = 0._wp
      zuw_bse(:,:) = 0._wp ; zvw_bse(:,:) = 0._wp
      !
      zdtdz_pyc(:,:,:) = 0._wp ; zdsdz_pyc(:,:,:) = 0._wp ; zdbdz_pyc(:,:,:) = 0._wp
      zdudz_pyc(:,:,:) = 0._wp ; zdvdz_pyc(:,:,:) = 0._wp
      !
      ! Flux-Gradient arrays.
      zdifml_sc(:,:)  = 0._wp ; zvisml_sc(:,:)  = 0._wp ; zdifpyc_sc(:,:) = 0._wp
      zvispyc_sc(:,:) = 0._wp ; zbeta_d_sc(:,:) = 0._wp ; zbeta_v_sc(:,:) = 0._wp
      zsc_wth_1(:,:)  = 0._wp ; zsc_ws_1(:,:)   = 0._wp ; zsc_uw_1(:,:)   = 0._wp
      zsc_uw_2(:,:)   = 0._wp ; zsc_vw_1(:,:)   = 0._wp ; zsc_vw_2(:,:)   = 0._wp
      zhbl_t(:,:)     = 0._wp ; zdhdt(:,:)      = 0._wp

      zdiffut(:,:,:) = 0._wp ; zviscos(:,:,:) = 0._wp ; ghamt(:,:,:) = 0._wp
      ghams(:,:,:)   = 0._wp ; ghamu(:,:,:)   = 0._wp ; ghamv(:,:,:) = 0._wp

      ! hbl = MAX(hbl,epsln)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Calculate boundary layer scales
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      ! Assume two-band radiation model for depth of OSBL
     zz0 =       rn_abs       ! surface equi-partition in 2-bands
     zz1 =  1. - rn_abs
     DO jj = 2, jpjm1
        DO ji = 2, jpim1
           ! Surface downward irradiance (so always +ve)
           zrad0(ji,jj) = qsr(ji,jj) * r1_rau0_rcp
           ! Downwards irradiance at base of boundary layer
           zradh(ji,jj) = zrad0(ji,jj) * ( zz0 * EXP( -hbl(ji,jj)/rn_si0 ) + zz1 * EXP( -hbl(ji,jj)/rn_si1) )
           ! Downwards irradiance averaged over depth of the OSBL
           zradav(ji,jj) = zrad0(ji,jj) * ( zz0 * ( 1.0 - EXP( -hbl(ji,jj)/rn_si0 ) )*rn_si0 &
                 &                         + zz1 * ( 1.0 - EXP( -hbl(ji,jj)/rn_si1 ) )*rn_si1 ) / hbl(ji,jj)
        END DO
     END DO
     ! Turbulent surface fluxes and fluxes averaged over depth of the OSBL
     DO jj = 2, jpjm1
        DO ji = 2, jpim1
           zthermal = rab_n(ji,jj,1,jp_tem)
           zbeta    = rab_n(ji,jj,1,jp_sal)
           ! Upwards surface Temperature flux for non-local term
           zwth0(ji,jj) = - qns(ji,jj) * r1_rau0_rcp * tmask(ji,jj,1)
           ! Upwards surface salinity flux for non-local term
           zws0(ji,jj) = - ( ( emp(ji,jj)-rnf(ji,jj) ) * tsn(ji,jj,1,jp_sal)  + sfx(ji,jj) ) * r1_rau0 * tmask(ji,jj,1)
           ! Non radiative upwards surface buoyancy flux
           zwb0(ji,jj) = grav * zthermal * zwth0(ji,jj) -  grav * zbeta * zws0(ji,jj)
           ! turbulent heat flux averaged over depth of OSBL
           zwthav(ji,jj) = 0.5 * zwth0(ji,jj) - ( 0.5*( zrad0(ji,jj) + zradh(ji,jj) ) - zradav(ji,jj) )
           ! turbulent salinity flux averaged over depth of the OBSL
           zwsav(ji,jj) = 0.5 * zws0(ji,jj)
           ! turbulent buoyancy flux averaged over the depth of the OBSBL
           zwbav(ji,jj) = grav  * zthermal * zwthav(ji,jj) - grav  * zbeta * zwsav(ji,jj)
           ! Surface upward velocity fluxes
           zuw0(ji,jj) = -utau(ji,jj) * r1_rau0 * tmask(ji,jj,1)
           zvw0(ji,jj) = -vtau(ji,jj) * r1_rau0 * tmask(ji,jj,1)
           ! Friction velocity (zustar), at T-point : LMD94 eq. 2
           zustar(ji,jj) = MAX( SQRT( SQRT( zuw0(ji,jj) * zuw0(ji,jj) + zvw0(ji,jj) * zvw0(ji,jj) ) ), 1.0e-8 )
           zcos_wind(ji,jj) = -zuw0(ji,jj) / ( zustar(ji,jj) * zustar(ji,jj) )
           zsin_wind(ji,jj) = -zvw0(ji,jj) / ( zustar(ji,jj) * zustar(ji,jj) )
        END DO
     END DO
     ! Calculate Stokes drift in direction of wind (zustke) and Stokes penetration depth (dstokes)
     SELECT CASE (nn_osm_wave)
     ! Assume constant La#=0.3
     CASE(0)
        DO jj = 2, jpjm1
           DO ji = 2, jpim1
              zus_x = zcos_wind(ji,jj) * zustar(ji,jj) / 0.3**2
              zus_y = zsin_wind(ji,jj) * zustar(ji,jj) / 0.3**2
              zustke(ji,jj) = MAX ( SQRT( zus_x*zus_x + zus_y*zus_y), 1.0e-8 )
              ! dstokes(ji,jj) set to constant value rn_osm_dstokes from namelist in zdf_osm_init
           END DO
        END DO
     ! Assume Pierson-Moskovitz wind-wave spectrum
     CASE(1)
        DO jj = 2, jpjm1
           DO ji = 2, jpim1
              ! Use wind speed wndm included in sbc_oce module
              zustke(ji,jj) = MAX ( 0.016 * wndm(ji,jj), 1.0e-8 )
              dstokes(ji,jj) = 0.12 * wndm(ji,jj)**2 / grav
           END DO
        END DO
     ! Use ECMWF wave fields as output from SBCWAVE
     CASE(2)
        zfac =  2.0_wp * rpi / 16.0_wp
        DO jj = 2, jpjm1
           DO ji = 2, jpim1
              ! The Langmur number from the ECMWF model appears to give La<0.3 for wind-driven seas.
              !    The coefficient 0.8 gives La=0.3  in this situation.
              ! It could represent the effects of the spread of wave directions
              ! around the mean wind. The effect of this adjustment needs to be tested.
              zustke(ji,jj) = MAX ( 1.0 * ( zcos_wind(ji,jj) * ut0sd(ji,jj ) + zsin_wind(ji,jj)  * vt0sd(ji,jj) ), &
                   &                zustar(ji,jj) / ( 0.45 * 0.45 )                                                  )
              dstokes(ji,jj) = MAX(zfac * hsw(ji,jj)*hsw(ji,jj) / ( MAX(zustke(ji,jj)*wmp(ji,jj), 1.0e-7 ) ), 5.0e-1) !rn_osm_dstokes !
           END DO
        END DO
     END SELECT

     ! Langmuir velocity scale (zwstrl), La # (zla)
     ! mixed scale (zvstr), convective velocity scale (zwstrc)
     DO jj = 2, jpjm1
        DO ji = 2, jpim1
           ! Langmuir velocity scale (zwstrl), at T-point
           zwstrl(ji,jj) = ( zustar(ji,jj) * zustar(ji,jj) * zustke(ji,jj) )**pthird
           ! Modify zwstrl to allow for small and large values of dstokes/hbl.
           ! Intended as a possible test. Doesn't affect LES results for entrainment,
           !  but hasn't been shown to be correct as dstokes/h becomes large or small.
           zwstrl(ji,jj) = zwstrl(ji,jj) *  &
                & (1.12 * ( 1.0 - ( 1.0 - EXP( -hbl(ji,jj) / dstokes(ji,jj) ) ) * dstokes(ji,jj) / hbl(ji,jj) ))**pthird * &
                & ( 1.0 - EXP( -15.0 * dstokes(ji,jj) / hbl(ji,jj) ))
           ! define La this way so effects of Stokes penetration depth on velocity scale are included
           zla(ji,jj) = SQRT ( zustar(ji,jj) / zwstrl(ji,jj) )**3
           ! Velocity scale that tends to zustar for large Langmuir numbers
           zvstr(ji,jj) = ( zwstrl(ji,jj)**3  + &
                & ( 1.0 - EXP( -0.5 * zla(ji,jj)**2 ) ) * zustar(ji,jj) * zustar(ji,jj) * zustar(ji,jj) )**pthird

           ! limit maximum value of Langmuir number as approximate treatment for shear turbulence.
           ! Note zustke and zwstrl are not amended.
           IF ( zla(ji,jj) >= 0.45 ) zla(ji,jj) = 0.45
           !
           ! get convective velocity (zwstrc), stabilty scale (zhol) and logical conection flag lconv
           IF ( zwbav(ji,jj) > 0.0) THEN
              zwstrc(ji,jj) = ( 2.0 * zwbav(ji,jj) * 0.9 * hbl(ji,jj) )**pthird
              zhol(ji,jj) = -0.9 * hbl(ji,jj) * 2.0 * zwbav(ji,jj) / (zvstr(ji,jj)**3 + epsln )
              lconv(ji,jj) = .TRUE.
           ELSE
              zhol(ji,jj) = -hbl(ji,jj) *  2.0 * zwbav(ji,jj)/ (zvstr(ji,jj)**3  + epsln )
              lconv(ji,jj) = .FALSE.
           ENDIF
        END DO
     END DO

     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     ! Mixed-layer model - calculate averages over the boundary layer, and the change in the boundary layer depth
     !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     ! BL must be always 2 levels deep.
      hbl(:,:) = MAX(hbl(:,:), gdepw_n(:,:,3) )
      ibld(:,:) = 3
      DO jk = 4, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               IF ( hbl(ji,jj) >= gdepw_n(ji,jj,jk) ) THEN
                  ibld(ji,jj) = MIN(mbkt(ji,jj), jk)
               ENDIF
            END DO
         END DO
      END DO

      DO jj = 2, jpjm1                                 !  Vertical slab
         DO ji = 2, jpim1
               zthermal = rab_n(ji,jj,1,jp_tem) !ideally use ibld not 1??
               zbeta    = rab_n(ji,jj,1,jp_sal)
               zt   = 0._wp
               zs   = 0._wp
               zu   = 0._wp
               zv   = 0._wp
               ! average over depth of boundary layer
               zthick=0._wp
               DO jm = 2, ibld(ji,jj)
                  zthick=zthick+e3t_n(ji,jj,jm)
                  zt   = zt  + e3t_n(ji,jj,jm) * tsn(ji,jj,jm,jp_tem)
                  zs   = zs  + e3t_n(ji,jj,jm) * tsn(ji,jj,jm,jp_sal)
                  zu   = zu  + e3t_n(ji,jj,jm) &
                     &            * ( ub(ji,jj,jm) + ub(ji - 1,jj,jm) ) &
                     &            / MAX( 1. , umask(ji,jj,jm) + umask(ji - 1,jj,jm) )
                  zv   = zv  + e3t_n(ji,jj,jm) &
                     &            * ( vb(ji,jj,jm) + vb(ji,jj - 1,jm) ) &
                     &            / MAX( 1. , vmask(ji,jj,jm) + vmask(ji,jj - 1,jm) )
               END DO
               zt_bl(ji,jj) = zt / zthick
               zs_bl(ji,jj) = zs / zthick
               zu_bl(ji,jj) = zu / zthick
               zv_bl(ji,jj) = zv / zthick
               zdt_bl(ji,jj) = zt_bl(ji,jj) - tsn(ji,jj,ibld(ji,jj),jp_tem)
               zds_bl(ji,jj) = zs_bl(ji,jj) - tsn(ji,jj,ibld(ji,jj),jp_sal)
               zdu_bl(ji,jj) = zu_bl(ji,jj) - ( ub(ji,jj,ibld(ji,jj)) + ub(ji-1,jj,ibld(ji,jj) ) ) &
                     &    / MAX(1. , umask(ji,jj,ibld(ji,jj) ) + umask(ji-1,jj,ibld(ji,jj) ) )
               zdv_bl(ji,jj) = zv_bl(ji,jj) - ( vb(ji,jj,ibld(ji,jj)) + vb(ji,jj-1,ibld(ji,jj) ) ) &
                     &   / MAX(1. , vmask(ji,jj,ibld(ji,jj) ) + vmask(ji,jj-1,ibld(ji,jj) ) )
               zdb_bl(ji,jj) = grav * zthermal * zdt_bl(ji,jj) - grav * zbeta * zds_bl(ji,jj)
               IF ( lconv(ji,jj) ) THEN    ! Convective
                      zwb_ent(ji,jj) = -( 2.0 * 0.2 * zwbav(ji,jj) &
                           &            + 0.135 * zla(ji,jj) * zwstrl(ji,jj)**3/hbl(ji,jj) )

                      zvel_max =  - ( 1.0 + 1.0 * ( zwstrl(ji,jj)**3 + 0.5 * zwstrc(ji,jj)**3 )**pthird * rn_rdt / hbl(ji,jj) ) &
                           &   * zwb_ent(ji,jj) / ( zwstrl(ji,jj)**3 + 0.5 * zwstrc(ji,jj)**3 )**pthird
! Entrainment including component due to shear turbulence. Modified Langmuir component, but gives same result for La=0.3 For testing uncomment.
!                      zwb_ent(ji,jj) = -( 2.0 * 0.2 * zwbav(ji,jj) &
!                           &            + ( 0.15 * ( 1.0 - EXP( -0.5 * zla(ji,jj) ) ) + 0.03 / zla(ji,jj)**2 ) * zustar(ji,jj)**3/hbl(ji,jj) )

!                      zvel_max = - ( 1.0 + 1.0 * ( zvstr(ji,jj)**3 + 0.5 * zwstrc(ji,jj)**3 )**pthird * rn_rdt / zhbl(ji,jj) ) * zwb_ent(ji,jj) / &
!                           &       ( zvstr(ji,jj)**3 + 0.5 * zwstrc(ji,jj)**3 )**pthird
                      zzdhdt = - zwb_ent(ji,jj) / ( zvel_max + MAX(zdb_bl(ji,jj),0.0) )
               ELSE                        ! Stable
                      zzdhdt = 0.32 * ( hbli(ji,jj) / hbl(ji,jj) -1.0 ) * zwstrl(ji,jj)**3 / hbli(ji,jj) &
                           &   + ( ( 0.32 / 3.0 ) * exp ( -2.5 * ( hbli(ji,jj) / hbl(ji,jj) - 1.0 ) ) &
                           & - ( 0.32 / 3.0 - 0.135 * zla(ji,jj) ) * exp ( -12.5 * ( hbli(ji,jj) / hbl(ji,jj) ) ) ) &
                           &  * zwstrl(ji,jj)**3 / hbli(ji,jj)
                      zzdhdt = zzdhdt + zwbav(ji,jj)
                      IF ( zzdhdt < 0._wp ) THEN
                      ! For long timsteps factor in brackets slows the rapid collapse of the OSBL
                         zpert   = 2.0 * ( 1.0 + 2.0 * zwstrl(ji,jj) * rn_rdt / hbl(ji,jj) ) * zwstrl(ji,jj)**2 / hbl(ji,jj)
                      ELSE
                         zpert   = 2.0 * ( 1.0 + 2.0 * zwstrl(ji,jj) * rn_rdt / hbl(ji,jj) ) * zwstrl(ji,jj)**2 / hbl(ji,jj) &
                              &  + MAX( zdb_bl(ji,jj), 0.0 )
                      ENDIF
                      zzdhdt = 2.0 * zzdhdt / zpert
               ENDIF
               zdhdt(ji,jj) = zzdhdt
           END DO
      END DO

      ! Calculate averages over depth of boundary layer
      imld = ibld           ! use imld to hold previous blayer index
      ibld(:,:) = 3

      zhbl_t(:,:) = hbl(:,:) + (zdhdt(:,:) - wn(ji,jj,ibld(ji,jj)))* rn_rdt ! certainly need wb here, so subtract it
      zhbl_t(:,:) = MIN(zhbl_t(:,:), ht_n(:,:))
      zdhdt(:,:) = MIN(zdhdt(:,:), (zhbl_t(:,:) - hbl(:,:))/rn_rdt + wn(ji,jj,ibld(ji,jj))) ! adjustment to represent limiting by ocean bottom

      DO jk = 4, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               IF ( zhbl_t(ji,jj) >= gdepw_n(ji,jj,jk) ) THEN
                  ibld(ji,jj) =  MIN(mbkt(ji,jj), jk)
               ENDIF
            END DO
         END DO
      END DO

!
! Step through model levels taking account of buoyancy change to determine the effect on dhdt
!
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            IF ( ibld(ji,jj) - imld(ji,jj) > 1 ) THEN
!
! If boundary layer changes by more than one level, need to check for stable layers between initial and final depths.
!
               zhbl_s = hbl(ji,jj)
               jm = imld(ji,jj)
               zthermal = rab_n(ji,jj,1,jp_tem)
               zbeta = rab_n(ji,jj,1,jp_sal)
               IF ( lconv(ji,jj) ) THEN
!unstable
                  zvel_max =  - ( 1.0 + 1.0 * ( zvstr(ji,jj)**3 + 0.5 * zwstrc(ji,jj)**3 )**pthird * rn_rdt / hbl(ji,jj) ) &
                       &   * zwb_ent(ji,jj) / ( zvstr(ji,jj)**3 + 0.5 * zwstrc(ji,jj)**3 )**pthird

                  DO jk = imld(ji,jj), ibld(ji,jj)
                     zdb = MAX( grav * ( zthermal * ( zt_bl(ji,jj) - tsn(ji,jj,jm,jp_tem) ) &
                          & - zbeta * ( zs_bl(ji,jj) - tsn(ji,jj,jm,jp_sal) ) ), 0.0 ) + zvel_max

                     zhbl_s = zhbl_s + MIN( - zwb_ent(ji,jj) / zdb * rn_rdt / FLOAT(ibld(ji,jj)-imld(ji,jj) ), e3w_n(ji,jj,jk) )
                     zhbl_s = MIN(zhbl_s, ht_n(ji,jj))

                     IF ( zhbl_s >= gdepw_n(ji,jj,jm+1) ) jm = jm + 1
                  END DO
                  hbl(ji,jj) = zhbl_s
                  ibld(ji,jj) = jm
                  hbli(ji,jj) = hbl(ji,jj)
               ELSE
! stable
                  DO jk = imld(ji,jj), ibld(ji,jj)
                     zdb = MAX( grav * ( zthermal * ( zt_bl(ji,jj) - tsn(ji,jj,jm,jp_tem) )          &
                          &               - zbeta * ( zs_bl(ji,jj) - tsn(ji,jj,jm,jp_sal) ) ), 0.0 ) &
                          & + 2.0 * zwstrl(ji,jj)**2 / zhbl_s

                     zhbl_s = zhbl_s +  (                                                                                &
                          &                     0.32         *                         ( hbli(ji,jj) / zhbl_s -1.0 )     &
                          &               * zwstrl(ji,jj)**3 / hbli(ji,jj)                                               &
                          &               + ( ( 0.32 / 3.0 )           * EXP( -  2.5 * ( hbli(ji,jj) / zhbl_s -1.0 ) )   &
                          &               -   ( 0.32 / 3.0  - 0.0485 ) * EXP( - 12.5 * ( hbli(ji,jj) / zhbl_s      ) ) ) &
                          &          * zwstrl(ji,jj)**3 / hbli(ji,jj) ) / zdb * e3w_n(ji,jj,jk) / zdhdt(ji,jj)  ! ALMG to investigate whether need to include wn here

                     zhbl_s = MIN(zhbl_s, ht_n(ji,jj))
                     IF ( zhbl_s >= gdepw_n(ji,jj,jm) ) jm = jm + 1
                  END DO
                  hbl(ji,jj) = MAX(zhbl_s, gdepw_n(ji,jj,3) )
                  ibld(ji,jj) = MAX(jm, 3 )
                  IF ( hbl(ji,jj) > hbli(ji,jj) ) hbli(ji,jj) = hbl(ji,jj)
               ENDIF   ! IF ( lconv )
            ELSE
! change zero or one model level.
               hbl(ji,jj) = zhbl_t(ji,jj)
               IF ( lconv(ji,jj) ) THEN
                  hbli(ji,jj) = hbl(ji,jj)
               ELSE
                  hbl(ji,jj) = MAX(hbl(ji,jj), gdepw_n(ji,jj,3) )
                  IF ( hbl(ji,jj) > hbli(ji,jj) ) hbli(ji,jj) = hbl(ji,jj)
               ENDIF
            ENDIF
            zhbl(ji,jj) = gdepw_n(ji,jj,ibld(ji,jj))
         END DO
      END DO
      dstokes(:,:) = MIN ( dstokes(:,:), hbl(:,:)/3. )  !  Limit delta for shallow boundary layers for calculating flux-gradient terms.

! Recalculate averages over boundary layer after depth updated
     ! Consider later  combining this into the loop above and looking for columns
     ! where the index for base of the boundary layer have changed
      DO jj = 2, jpjm1                                 !  Vertical slab
         DO ji = 2, jpim1
               zthermal = rab_n(ji,jj,1,jp_tem) !ideally use ibld not 1??
               zbeta    = rab_n(ji,jj,1,jp_sal)
               zt   = 0._wp
               zs   = 0._wp
               zu   = 0._wp
               zv   = 0._wp
               ! average over depth of boundary layer
               zthick=0._wp
               DO jm = 2, ibld(ji,jj)
                  zthick=zthick+e3t_n(ji,jj,jm)
                  zt   = zt  + e3t_n(ji,jj,jm) * tsn(ji,jj,jm,jp_tem)
                  zs   = zs  + e3t_n(ji,jj,jm) * tsn(ji,jj,jm,jp_sal)
                  zu   = zu  + e3t_n(ji,jj,jm) &
                     &            * ( ub(ji,jj,jm) + ub(ji - 1,jj,jm) ) &
                     &            / MAX( 1. , umask(ji,jj,jm) + umask(ji - 1,jj,jm) )
                  zv   = zv  + e3t_n(ji,jj,jm) &
                     &            * ( vb(ji,jj,jm) + vb(ji,jj - 1,jm) ) &
                     &            / MAX( 1. , vmask(ji,jj,jm) + vmask(ji,jj - 1,jm) )
               END DO
               zt_bl(ji,jj) = zt / zthick
               zs_bl(ji,jj) = zs / zthick
               zu_bl(ji,jj) = zu / zthick
               zv_bl(ji,jj) = zv / zthick
               zdt_bl(ji,jj) = zt_bl(ji,jj) - tsn(ji,jj,ibld(ji,jj),jp_tem)
               zds_bl(ji,jj) = zs_bl(ji,jj) - tsn(ji,jj,ibld(ji,jj),jp_sal)
               zdu_bl(ji,jj) = zu_bl(ji,jj) - ( ub(ji,jj,ibld(ji,jj)) + ub(ji-1,jj,ibld(ji,jj) ) ) &
                      &   / MAX(1. , umask(ji,jj,ibld(ji,jj) ) + umask(ji-1,jj,ibld(ji,jj) ) )
               zdv_bl(ji,jj) = zv_bl(ji,jj) - ( vb(ji,jj,ibld(ji,jj)) + vb(ji,jj-1,ibld(ji,jj) ) ) &
                      &  / MAX(1. , vmask(ji,jj,ibld(ji,jj) ) + vmask(ji,jj-1,ibld(ji,jj) ) )
               zdb_bl(ji,jj) = grav * zthermal * zdt_bl(ji,jj) - grav * zbeta * zds_bl(ji,jj)
               zhbl(ji,jj) = gdepw_n(ji,jj,ibld(ji,jj))
               IF ( lconv(ji,jj) ) THEN
                  IF ( zdb_bl(ji,jj) > 0._wp )THEN
                     IF ( ( zwstrc(ji,jj) / zvstr(ji,jj) )**3 <= 0.5 ) THEN  ! near neutral stability
                           zari = 4.5 * ( zvstr(ji,jj)**2 ) &
                             & / ( zdb_bl(ji,jj) * zhbl(ji,jj) ) + 0.01
                     ELSE                                                     ! unstable
                           zari = 4.5 * ( zwstrc(ji,jj)**2 ) &
                             & / ( zdb_bl(ji,jj) * zhbl(ji,jj) ) + 0.01
                     ENDIF
                     IF ( zari > 0.2 ) THEN                                                ! This test checks for weakly stratified pycnocline
                        zari = 0.2
                        zwb_ent(ji,jj) = 0._wp
                     ENDIF
                     inhml = MAX( INT( zari * zhbl(ji,jj) / e3t_n(ji,jj,ibld(ji,jj)) ) , 1 )
                     imld(ji,jj) = MAX( ibld(ji,jj) - inhml, 1)
                     zhml(ji,jj) = gdepw_n(ji,jj,imld(ji,jj))
                     zdh(ji,jj) = zhbl(ji,jj) - zhml(ji,jj)
                  ELSE  ! IF (zdb_bl)
                     imld(ji,jj) = ibld(ji,jj) - 1
                     zhml(ji,jj) = gdepw_n(ji,jj,imld(ji,jj))
                     zdh(ji,jj) = zhbl(ji,jj) - zhml(ji,jj)
                  ENDIF
               ELSE   ! IF (lconv)
                  IF ( zdhdt(ji,jj) >= 0.0 ) THEN    ! probably shouldn't include wm here
                  ! boundary layer deepening
                     IF ( zdb_bl(ji,jj) > 0._wp ) THEN
                  ! pycnocline thickness set by stratification - use same relationship as for neutral conditions.
                        zari = MIN( 4.5 * ( zvstr(ji,jj)**2 ) &
                          & / ( zdb_bl(ji,jj) * zhbl(ji,jj) ) + 0.01  , 0.2 )
                        inhml = MAX( INT( zari * zhbl(ji,jj) / e3t_n(ji,jj,ibld(ji,jj)) ) , 1 )
                        imld(ji,jj) = MAX( ibld(ji,jj) - inhml, 1)
                        zhml(ji,jj) = gdepw_n(ji,jj,imld(ji,jj))
                        zdh(ji,jj) = zhbl(ji,jj) - zhml(ji,jj)
                     ELSE
                        imld(ji,jj) = ibld(ji,jj) - 1
                        zhml(ji,jj) = gdepw_n(ji,jj,imld(ji,jj))
                        zdh(ji,jj) = zhbl(ji,jj) - zhml(ji,jj)
                     ENDIF ! IF (zdb_bl > 0.0)
                  ELSE     ! IF(dhdt >= 0)
                  ! boundary layer collapsing.
                     imld(ji,jj) = ibld(ji,jj)
                     zhml(ji,jj) = zhbl(ji,jj)
                     zdh(ji,jj) = 0._wp
                  ENDIF    ! IF (dhdt >= 0)
               ENDIF       ! IF (lconv)
         END DO
      END DO

      ! Average over the depth of the mixed layer in the convective boundary layer
      ! Also calculate entrainment fluxes for temperature and salinity
      DO jj = 2, jpjm1                                 !  Vertical slab
         DO ji = 2, jpim1
            zthermal = rab_n(ji,jj,1,jp_tem) !ideally use ibld not 1??
            zbeta    = rab_n(ji,jj,1,jp_sal)
            IF ( lconv(ji,jj) ) THEN
               zt   = 0._wp
               zs   = 0._wp
               zu   = 0._wp
               zv   = 0._wp
               ! average over depth of boundary layer
               zthick=0._wp
               DO jm = 2, imld(ji,jj)
                  zthick=zthick+e3t_n(ji,jj,jm)
                  zt   = zt  + e3t_n(ji,jj,jm) * tsn(ji,jj,jm,jp_tem)
                  zs   = zs  + e3t_n(ji,jj,jm) * tsn(ji,jj,jm,jp_sal)
                  zu   = zu  + e3t_n(ji,jj,jm) &
                     &            * ( ub(ji,jj,jm) + ub(ji - 1,jj,jm) ) &
                     &            / MAX( 1. , umask(ji,jj,jm) + umask(ji - 1,jj,jm) )
                  zv   = zv  + e3t_n(ji,jj,jm) &
                     &            * ( vb(ji,jj,jm) + vb(ji,jj - 1,jm) ) &
                     &            / MAX( 1. , vmask(ji,jj,jm) + vmask(ji,jj - 1,jm) )
               END DO
               zt_ml(ji,jj) = zt / zthick
               zs_ml(ji,jj) = zs / zthick
               zu_ml(ji,jj) = zu / zthick
               zv_ml(ji,jj) = zv / zthick
               zdt_ml(ji,jj) = zt_ml(ji,jj) - tsn(ji,jj,ibld(ji,jj),jp_tem)
               zds_ml(ji,jj) = zs_ml(ji,jj) - tsn(ji,jj,ibld(ji,jj),jp_sal)
               zdu_ml(ji,jj) = zu_ml(ji,jj) - ( ub(ji,jj,ibld(ji,jj)) + ub(ji-1,jj,ibld(ji,jj) ) ) &
                     &    / MAX(1. , umask(ji,jj,ibld(ji,jj) ) + umask(ji-1,jj,ibld(ji,jj) ) )
               zdv_ml(ji,jj) = zv_ml(ji,jj) - ( vb(ji,jj,ibld(ji,jj)) + vb(ji,jj-1,ibld(ji,jj) ) ) &
                     &    / MAX(1. , vmask(ji,jj,ibld(ji,jj) ) + vmask(ji,jj-1,ibld(ji,jj) ) )
               zdb_ml(ji,jj) = grav * zthermal * zdt_ml(ji,jj) - grav * zbeta * zds_ml(ji,jj)
            ELSE
            ! stable, if entraining calulate average below interface layer.
               IF ( zdhdt(ji,jj) >= 0._wp ) THEN
                  zt   = 0._wp
                  zs   = 0._wp
                  zu   = 0._wp
                  zv   = 0._wp
                  ! average over depth of boundary layer
                  zthick=0._wp
                  DO jm = 2, imld(ji,jj)
                     zthick=zthick+e3t_n(ji,jj,jm)
                     zt   = zt  + e3t_n(ji,jj,jm) * tsn(ji,jj,jm,jp_tem)
                     zs   = zs  + e3t_n(ji,jj,jm) * tsn(ji,jj,jm,jp_sal)
                     zu   = zu  + e3t_n(ji,jj,jm) &
                        &            * ( ub(ji,jj,jm) + ub(ji - 1,jj,jm) ) &
                        &            / MAX( 1. , umask(ji,jj,jm) + umask(ji - 1,jj,jm) )
                     zv   = zv  + e3t_n(ji,jj,jm) &
                        &            * ( vb(ji,jj,jm) + vb(ji,jj - 1,jm) ) &
                        &            / MAX( 1. , vmask(ji,jj,jm) + vmask(ji,jj - 1,jm) )
                  END DO
                  zt_ml(ji,jj) = zt / zthick
                  zs_ml(ji,jj) = zs / zthick
                  zu_ml(ji,jj) = zu / zthick
                  zv_ml(ji,jj) = zv / zthick
                  zdt_ml(ji,jj) = zt_ml(ji,jj) - tsn(ji,jj,ibld(ji,jj),jp_tem)
                  zds_ml(ji,jj) = zs_ml(ji,jj) - tsn(ji,jj,ibld(ji,jj),jp_sal)
                  zdu_ml(ji,jj) = zu_ml(ji,jj) - ( ub(ji,jj,ibld(ji,jj)) + ub(ji-1,jj,ibld(ji,jj) ) ) &
                        &    / MAX(1. , umask(ji,jj,ibld(ji,jj) ) + umask(ji-1,jj,ibld(ji,jj) ) )
                  zdv_ml(ji,jj) = zv_ml(ji,jj) - ( vb(ji,jj,ibld(ji,jj)) + vb(ji,jj-1,ibld(ji,jj) ) ) &
                        &    / MAX(1. , vmask(ji,jj,ibld(ji,jj) ) + vmask(ji,jj-1,ibld(ji,jj) ) )
                  zdb_ml(ji,jj) = grav * zthermal * zdt_ml(ji,jj) - grav * zbeta * zds_ml(ji,jj)
               ENDIF
            ENDIF
         END DO
      END DO
    !
    ! rotate mean currents and changes onto wind align co-ordinates
    !

      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            ztemp = zu_ml(ji,jj)
            zu_ml(ji,jj) = zu_ml(ji,jj) * zcos_wind(ji,jj) + zv_ml(ji,jj) * zsin_wind(ji,jj)
            zv_ml(ji,jj) = zv_ml(ji,jj) * zcos_wind(ji,jj) - ztemp * zsin_wind(ji,jj)
            ztemp = zdu_ml(ji,jj)
            zdu_ml(ji,jj) = zdu_ml(ji,jj) * zcos_wind(ji,jj) + zdv_ml(ji,jj) * zsin_wind(ji,jj)
            zdv_ml(ji,jj) = zdv_ml(ji,jj) * zsin_wind(ji,jj) - ztemp * zsin_wind(ji,jj)
    !
            ztemp = zu_bl(ji,jj)
            zu_bl = zu_bl(ji,jj) * zcos_wind(ji,jj) + zv_bl(ji,jj) * zsin_wind(ji,jj)
            zv_bl(ji,jj) = zv_bl(ji,jj) * zcos_wind(ji,jj) - ztemp * zsin_wind(ji,jj)
            ztemp = zdu_bl(ji,jj)
            zdu_bl(ji,jj) = zdu_bl(ji,jj) * zcos_wind(ji,jj) + zdv_bl(ji,jj) * zsin_wind(ji,jj)
            zdv_bl(ji,jj) = zdv_bl(ji,jj) * zsin_wind(ji,jj) - ztemp * zsin_wind(ji,jj)
         END DO
      END DO

     zuw_bse = 0._wp
     zvw_bse = 0._wp
     DO jj = 2, jpjm1
        DO ji = 2, jpim1

           IF ( lconv(ji,jj) ) THEN
              IF ( zdb_bl(ji,jj) > 0._wp ) THEN
                 zwth_ent(ji,jj) = zwb_ent(ji,jj) * zdt_ml(ji,jj) / (zdb_ml(ji,jj) + epsln)
                 zws_ent(ji,jj) = zwb_ent(ji,jj) * zds_ml(ji,jj) / (zdb_ml(ji,jj) + epsln)
              ENDIF
           ELSE
              zwth_ent(ji,jj) = -2.0 * zwthav(ji,jj) * ( (1.0 - 0.8) - ( 1.0 - 0.8)**(3.0/2.0) )
              zws_ent(ji,jj) = -2.0 * zwsav(ji,jj) * ( (1.0 - 0.8 ) - ( 1.0 - 0.8 )**(3.0/2.0) )
           ENDIF
        END DO
     END DO

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !  Pycnocline gradients for scalars and velocity
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

       DO jj = 2, jpjm1
          DO ji = 2, jpim1
          !
             IF ( lconv (ji,jj) ) THEN
             ! Unstable conditions
                IF( zdb_bl(ji,jj) > 0._wp ) THEN
                ! calculate pycnocline profiles, no need if zdb_bl <= 0. since profile is zero and arrays have been initialized to zero
                   ztgrad = ( zdt_ml(ji,jj) / zdh(ji,jj) )
                   zsgrad = ( zds_ml(ji,jj) / zdh(ji,jj) )
                   zbgrad = ( zdb_ml(ji,jj) / zdh(ji,jj) )
                   DO jk = 2 , ibld(ji,jj)
                      znd = -( gdepw_n(ji,jj,jk) - zhml(ji,jj) ) / zdh(ji,jj)
                      zdtdz_pyc(ji,jj,jk) =  ztgrad * EXP( -1.75 * ( znd + 0.75 )**2 )
                      zdbdz_pyc(ji,jj,jk) =  zbgrad * EXP( -1.75 * ( znd + 0.75 )**2 )
                      zdsdz_pyc(ji,jj,jk) =  zsgrad * EXP( -1.75 * ( znd + 0.75 )**2 )
                   END DO
                ENDIF
             ELSE
             ! stable conditions
             ! if pycnocline profile only defined when depth steady of increasing.
                IF ( zdhdt(ji,jj) >= 0.0 ) THEN        ! Depth increasing, or steady.
                   IF ( zdb_bl(ji,jj) > 0._wp ) THEN
                     IF ( zhol(ji,jj) >= 0.5 ) THEN      ! Very stable - 'thick' pycnocline
                         ztgrad = zdt_bl(ji,jj) / zhbl(ji,jj)
                         zsgrad = zds_bl(ji,jj) / zhbl(ji,jj)
                         zbgrad = zdb_bl(ji,jj) / zhbl(ji,jj)
                         DO jk = 2, ibld(ji,jj)
                            znd = gdepw_n(ji,jj,jk) / zhbl(ji,jj)
                            zdtdz_pyc(ji,jj,jk) =  ztgrad * EXP( -15.0 * ( znd - 0.9 )**2 )
                            zdbdz_pyc(ji,jj,jk) =  zbgrad * EXP( -15.0 * ( znd - 0.9 )**2 )
                            zdsdz_pyc(ji,jj,jk) =  zsgrad * EXP( -15.0 * ( znd - 0.9 )**2 )
                         END DO
                     ELSE                                   ! Slightly stable - 'thin' pycnoline - needed when stable layer begins to form.
                         ztgrad = zdt_bl(ji,jj) / zdh(ji,jj)
                         zsgrad = zds_bl(ji,jj) / zdh(ji,jj)
                         zbgrad = zdb_bl(ji,jj) / zdh(ji,jj)
                         DO jk = 2, ibld(ji,jj)
                            znd = -( gdepw_n(ji,jj,jk) - zhml(ji,jj) ) / zdh(ji,jj)
                            zdtdz_pyc(ji,jj,jk) =  ztgrad * EXP( -1.75 * ( znd + 0.75 )**2 )
                            zdbdz_pyc(ji,jj,jk) =  zbgrad * EXP( -1.75 * ( znd + 0.75 )**2 )
                            zdsdz_pyc(ji,jj,jk) =  zsgrad * EXP( -1.75 * ( znd + 0.75 )**2 )
                         END DO
                      ENDIF ! IF (zhol >=0.5)
                   ENDIF    ! IF (zdb_bl> 0.)
                ENDIF       ! IF (zdhdt >= 0) zdhdt < 0 not considered since pycnocline profile is zero, profile arrays are intialized to zero
             ENDIF          ! IF (lconv)
            !
          END DO
       END DO
!
       DO jj = 2, jpjm1
          DO ji = 2, jpim1
          !
             IF ( lconv (ji,jj) ) THEN
             ! Unstable conditions
                 zugrad = ( zdu_ml(ji,jj) / zdh(ji,jj) ) + 0.275 * zustar(ji,jj)*zustar(ji,jj) / &
               & (( zwstrl(ji,jj)**3 + 0.5 * zwstrc(ji,jj)**3 )**pthird * zhml(ji,jj) ) / zla(ji,jj)**(8.0/3.0)
                zvgrad = ( zdv_ml(ji,jj) / zdh(ji,jj) ) + 3.5 * ff_t(ji,jj) * zustke(ji,jj) / &
              & ( zwstrl(ji,jj)**3 + 0.5 * zwstrc(ji,jj)**3 )**pthird
                DO jk = 2 , ibld(ji,jj)-1
                   znd = -( gdepw_n(ji,jj,jk) - zhml(ji,jj) ) / zdh(ji,jj)
                   zdudz_pyc(ji,jj,jk) =  zugrad * EXP( -1.75 * ( znd + 0.75 )**2 )
                   zdvdz_pyc(ji,jj,jk) = zvgrad * EXP( -1.75 * ( znd + 0.75 )**2 )
                END DO
             ELSE
             ! stable conditions
                zugrad = 3.25 * zdu_bl(ji,jj) / zhbl(ji,jj)
                zvgrad = 2.75 * zdv_bl(ji,jj) / zhbl(ji,jj)
                DO jk = 2, ibld(ji,jj)
                   znd = gdepw_n(ji,jj,jk) / zhbl(ji,jj)
                   IF ( znd < 1.0 ) THEN
                      zdudz_pyc(ji,jj,jk) = zugrad * EXP( -40.0 * ( znd - 1.0 )**2 )
                   ELSE
                      zdudz_pyc(ji,jj,jk) = zugrad * EXP( -20.0 * ( znd - 1.0 )**2 )
                   ENDIF
                   zdvdz_pyc(ji,jj,jk) = zvgrad * EXP( -20.0 * ( znd - 0.85 )**2 )
                END DO
             ENDIF
            !
          END DO
       END DO
       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       ! Eddy viscosity/diffusivity and non-gradient terms in the flux-gradient relationship
       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      ! WHERE ( lconv )
      !     zdifml_sc = zhml * ( zwstrl**3 + 0.5 * zwstrc**3 )**pthird
      !     zvisml_sc = zdifml_sc
      !     zdifpyc_sc = 0.165 * ( zwstrl**3 + zwstrc**3 )**pthird * ( zhbl - zhml )
      !     zvispyc_sc = 0.142 * ( zwstrl**3 + 0.5 * zwstrc**3 )**pthird * ( zhbl - zhml )
      !     zbeta_d_sc = 1.0 - (0.165 / 0.8 * ( zhbl - zhml ) / zhbl )**p2third
      !     zbeta_v_sc = 1.0 -  2.0 * (0.142 /0.375) * (zhbl - zhml ) / zhml
      !  ELSEWHERE
      !     zdifml_sc = zwstrl * zhbl * EXP ( -( zhol / 0.183_wp )**2 )
      !     zvisml_sc = zwstrl * zhbl * EXP ( -( zhol / 0.183_wp )**2 )
      !  ENDWHERE
       DO jj = 2, jpjm1
          DO ji = 2, jpim1
             IF ( lconv(ji,jj) ) THEN
               zdifml_sc(ji,jj) = zhml(ji,jj) * ( zvstr(ji,jj)**3 + 0.5 * zwstrc(ji,jj)**3 )**pthird
               zvisml_sc(ji,jj) = zdifml_sc(ji,jj)
               zdifpyc_sc(ji,jj) = 0.165 * ( zvstr(ji,jj)**3 + 0.5 *zwstrc(ji,jj)**3 )**pthird * zdh(ji,jj)
               zvispyc_sc(ji,jj) = 0.142 * ( zvstr(ji,jj)**3 + 0.5 * zwstrc(ji,jj)**3 )**pthird * zdh(ji,jj)
               zbeta_d_sc(ji,jj) = 1.0 - (0.165 / 0.8 * zdh(ji,jj) / zhbl(ji,jj) )**p2third
               zbeta_v_sc(ji,jj) = 1.0 -  2.0 * (0.142 /0.375) * zdh(ji,jj) / zhml(ji,jj)
             ELSE
               zdifml_sc(ji,jj) = zvstr(ji,jj) * zhbl(ji,jj) * EXP ( -( zhol(ji,jj) / 0.6_wp )**2 )
               zvisml_sc(ji,jj) = zvstr(ji,jj) * zhbl(ji,jj) * EXP ( -( zhol(ji,jj) / 0.6_wp )**2 )
            END IF
        END DO
    END DO
!
       DO jj = 2, jpjm1
          DO ji = 2, jpim1
             IF ( lconv(ji,jj) ) THEN
                DO jk = 2, imld(ji,jj)   ! mixed layer diffusivity
                    zznd_ml = gdepw_n(ji,jj,jk) / zhml(ji,jj)
                    !
                    zdiffut(ji,jj,jk) = 0.8   * zdifml_sc(ji,jj) * zznd_ml * ( 1.0 - zbeta_d_sc(ji,jj) * zznd_ml    )**1.5
                    !
                    zviscos(ji,jj,jk) = 0.375 * zvisml_sc(ji,jj) * zznd_ml * ( 1.0 - zbeta_v_sc(ji,jj) * zznd_ml    ) &
                         &            *                                      ( 1.0 -               0.5 * zznd_ml**2 )
                END DO
                ! pycnocline - if present linear profile
                IF ( zdh(ji,jj) > 0._wp ) THEN
                   DO jk = imld(ji,jj)+1 , ibld(ji,jj)
                       zznd_pyc = -( gdepw_n(ji,jj,jk) - zhml(ji,jj) ) / zdh(ji,jj)
                       !
                       zdiffut(ji,jj,jk) = zdifpyc_sc(ji,jj) * ( 1.0 + zznd_pyc )
                       !
                       zviscos(ji,jj,jk) = zvispyc_sc(ji,jj) * ( 1.0 + zznd_pyc )
                   END DO
                ENDIF
                ! Temporay fix to ensure zdiffut is +ve; won't be necessary with wn taken out
                zdiffut(ji,jj,ibld(ji,jj)) = zdhdt(ji,jj)* e3t_n(ji,jj,ibld(ji,jj))
                ! could be taken out, take account of entrainment represents as a diffusivity
                ! should remove w from here, represents entrainment
             ELSE
             ! stable conditions
                DO jk = 2, ibld(ji,jj)
                   zznd_ml = gdepw_n(ji,jj,jk) / zhbl(ji,jj)
                   zdiffut(ji,jj,jk) = 0.75 * zdifml_sc(ji,jj) * zznd_ml * ( 1.0 - zznd_ml )**1.5
                   zviscos(ji,jj,jk) = 0.375 * zvisml_sc(ji,jj) * zznd_ml * (1.0 - zznd_ml) * ( 1.0 - zznd_ml**2 )
                END DO
             ENDIF   ! end if ( lconv )
!
          END DO  ! end of ji loop
       END DO     ! end of jj loop

       !
       ! calculate non-gradient components of the flux-gradient relationships
       !
! Stokes term in scalar flux, flux-gradient relationship
       WHERE ( lconv )
          zsc_wth_1 = zwstrl**3 * zwth0 / ( zvstr**3 + 0.5 * zwstrc**3 + epsln)
          !
          zsc_ws_1 = zwstrl**3 * zws0 / ( zvstr**3 + 0.5 * zwstrc**3 + epsln )
       ELSEWHERE
          zsc_wth_1 = 2.0 * zwthav
          !
          zsc_ws_1 = 2.0 * zwsav
       ENDWHERE


       DO jj = 2, jpjm1
          DO ji = 2, jpim1
            IF ( lconv(ji,jj) ) THEN
              DO jk = 2, imld(ji,jj)
                 zznd_d = gdepw_n(ji,jj,jk) / dstokes(ji,jj)
                 ghamt(ji,jj,jk) = ghamt(ji,jj,jk) + 1.35 * EXP ( -zznd_d ) * ( 1.0 - EXP ( -2.0 * zznd_d ) ) * zsc_wth_1(ji,jj)
                 !
                 ghams(ji,jj,jk) = ghams(ji,jj,jk) + 1.35 * EXP ( -zznd_d ) * ( 1.0 - EXP ( -2.0 * zznd_d ) ) *  zsc_ws_1(ji,jj)
              END DO ! end jk loop
            ELSE     ! else for if (lconv)
 ! Stable conditions
               DO jk = 2, ibld(ji,jj)
                  zznd_d=gdepw_n(ji,jj,jk) / dstokes(ji,jj)
                  ghamt(ji,jj,jk) = ghamt(ji,jj,jk) + 1.5 * EXP ( -0.9 * zznd_d ) &
                       &          *                 ( 1.0 - EXP ( -4.0 * zznd_d ) ) * zsc_wth_1(ji,jj)
                  !
                  ghams(ji,jj,jk) = ghams(ji,jj,jk) + 1.5 * EXP ( -0.9 * zznd_d ) &
                       &          *                 ( 1.0 - EXP ( -4.0 * zznd_d ) ) *  zsc_ws_1(ji,jj)
               END DO
            ENDIF               ! endif for check on lconv

          END DO  ! end of ji loop
       END DO     ! end of jj loop


! Stokes term in flux-gradient relationship (note in zsc_uw_n don't use zvstr since term needs to go to zero as zwstrl goes to zero)
       WHERE ( lconv )
          zsc_uw_1 = ( zwstrl**3 + 0.5 * zwstrc**3 )**pthird * zustke /( 1.0 - 1.0 * 6.5 * zla**(8.0/3.0) )
          zsc_uw_2 = ( zwstrl**3 + 0.5 * zwstrc**3 )**pthird * zustke / ( zla**(8.0/3.0) + epsln )
          zsc_vw_1 = ff_t * zhml * zustke**3 * zla**(8.0/3.0) / ( ( zvstr**3 + 0.5 * zwstrc**3 )**(2.0/3.0) + epsln )
       ELSEWHERE
          zsc_uw_1 = zustar**2
          zsc_vw_1 = ff_t * zhbl * zustke**3 * zla**(8.0/3.0) / (zvstr**2 + epsln)
       ENDWHERE

       DO jj = 2, jpjm1
          DO ji = 2, jpim1
             IF ( lconv(ji,jj) ) THEN
                DO jk = 2, imld(ji,jj)
                   zznd_d = gdepw_n(ji,jj,jk) / dstokes(ji,jj)
                   ghamu(ji,jj,jk) = ghamu(ji,jj,jk) +      ( -0.05 * EXP ( -0.4 * zznd_d )   * zsc_uw_1(ji,jj)   &
                        &          +                        0.00125 * EXP (      - zznd_d )   * zsc_uw_2(ji,jj) ) &
                        &          *                          ( 1.0 - EXP ( -2.0 * zznd_d ) )
!
                   ghamv(ji,jj,jk) = ghamv(ji,jj,jk) - 0.65 *  0.15 * EXP (      - zznd_d )                       &
                        &          *                          ( 1.0 - EXP ( -2.0 * zznd_d ) ) * zsc_vw_1(ji,jj)
                END DO   ! end jk loop
             ELSE
! Stable conditions
                DO jk = 2, ibld(ji,jj) ! corrected to ibld
                   zznd_d = gdepw_n(ji,jj,jk) / dstokes(ji,jj)
                   ghamu(ji,jj,jk) = ghamu(ji,jj,jk) - 0.75 *   1.3 * EXP ( -0.5 * zznd_d )                       &
                        &                                   * ( 1.0 - EXP ( -4.0 * zznd_d ) ) * zsc_uw_1(ji,jj)
                   ghamv(ji,jj,jk) = ghamv(ji,jj,jk) + 0._wp
                END DO   ! end jk loop
             ENDIF
          END DO  ! ji loop
       END DO     ! jj loo

! Buoyancy term in flux-gradient relationship [note : includes ROI ratio (X0.3) and pressure (X0.5)]

       WHERE ( lconv )
          zsc_wth_1 = zwbav * zwth0 * ( 1.0 + EXP ( 0.2 * zhol ) ) / ( zvstr**3 + 0.5 * zwstrc**3 + epsln )
          zsc_ws_1  = zwbav * zws0  * ( 1.0 + EXP ( 0.2 * zhol ) ) / ( zvstr**3 + 0.5 * zwstrc**3 + epsln )
       ELSEWHERE
          zsc_wth_1 = 0._wp
          zsc_ws_1 = 0._wp
       ENDWHERE

       DO jj = 2, jpjm1
          DO ji = 2, jpim1
             IF (lconv(ji,jj) ) THEN
                DO jk = 2, imld(ji,jj)
                   zznd_ml = gdepw_n(ji,jj,jk) / zhml(ji,jj)
                   ! calculate turbulent length scale
                   zl_c = 0.9 * ( 1.0 - EXP ( - 7.0 * ( zznd_ml - zznd_ml**3 / 3.0 ) ) )                                           &
                        &     * ( 1.0 - EXP ( -15.0 * (     1.1 - zznd_ml          ) ) )
                   zl_l = 2.0 * ( 1.0 - EXP ( - 2.0 * ( zznd_ml - zznd_ml**3 / 3.0 ) ) )                                           &
                        &     * ( 1.0 - EXP ( - 5.0 * (     1.0 - zznd_ml          ) ) ) * ( 1.0 + dstokes(ji,jj) / zhml (ji,jj) )
                   zl_eps = zl_l + ( zl_c - zl_l ) / ( 1.0 + EXP ( 3.0 * LOG10 ( - zhol(ji,jj) ) ) ) ** (3.0/2.0)
                   ! non-gradient buoyancy terms
                   ghamt(ji,jj,jk) = ghamt(ji,jj,jk) + 0.3 * 0.5 * zsc_wth_1(ji,jj) * zl_eps * zhml(ji,jj) / ( 0.15 + zznd_ml )
                   ghams(ji,jj,jk) = ghams(ji,jj,jk) + 0.3 * 0.5 *  zsc_ws_1(ji,jj) * zl_eps * zhml(ji,jj) / ( 0.15 + zznd_ml )
                END DO
             ELSE
                DO jk = 2, ibld(ji,jj)
                   ghamt(ji,jj,jk) = ghamt(ji,jj,jk) + zsc_wth_1(ji,jj)
                   ghams(ji,jj,jk) = ghams(ji,jj,jk) +  zsc_ws_1(ji,jj)
                END DO
             ENDIF
          END DO   ! ji loop
       END DO      ! jj loop


       WHERE ( lconv )
          zsc_uw_1 = -zwb0 * zustar**2 * zhml / ( zvstr**3 + 0.5 * zwstrc**3 + epsln )
          zsc_uw_2 =  zwb0 * zustke    * zhml / ( zvstr**3 + 0.5 * zwstrc**3 + epsln )**(2.0/3.0)
          zsc_vw_1 = 0._wp
       ELSEWHERE
         zsc_uw_1 = 0._wp
         zsc_vw_1 = 0._wp
       ENDWHERE

       DO jj = 2, jpjm1
          DO ji = 2, jpim1
             IF ( lconv(ji,jj) ) THEN
                DO jk = 2 , imld(ji,jj)
                   zznd_d = gdepw_n(ji,jj,jk) / dstokes(ji,jj)
                   ghamu(ji,jj,jk) = ghamu(ji,jj,jk) + 0.3 * 0.5 * ( zsc_uw_1(ji,jj) +   0.125 * EXP( -0.5 * zznd_d )     &
                        &                                                            * (   1.0 - EXP( -0.5 * zznd_d ) )   &
                        &                                          * zsc_uw_2(ji,jj)                                    )
                   ghamv(ji,jj,jk) = ghamv(ji,jj,jk) + zsc_vw_1(ji,jj)
                END DO  ! jk loop
             ELSE
             ! stable conditions
                DO jk = 2, ibld(ji,jj)
                   ghamu(ji,jj,jk) = ghamu(ji,jj,jk) + zsc_uw_1(ji,jj)
                   ghamv(ji,jj,jk) = ghamv(ji,jj,jk) + zsc_vw_1(ji,jj)
                END DO
             ENDIF
          END DO        ! ji loop
       END DO           ! jj loop

! Transport term in flux-gradient relationship [note : includes ROI ratio (X0.3) ]

       WHERE ( lconv )
          zsc_wth_1 = zwth0
          zsc_ws_1 = zws0
       ELSEWHERE
          zsc_wth_1 = 2.0 * zwthav
          zsc_ws_1 = zws0
       ENDWHERE

       DO jj = 2, jpjm1
          DO ji = 2, jpim1
            IF ( lconv(ji,jj) ) THEN
               DO jk = 2, imld(ji,jj)
                  zznd_ml=gdepw_n(ji,jj,jk) / zhml(ji,jj)
                  ghamt(ji,jj,jk) = ghamt(ji,jj,jk) + 0.3 * zsc_wth_1(ji,jj)                &
                       &          * ( -2.0 + 2.75 * (       ( 1.0 + 0.6 * zznd_ml**4 )      &
                       &                               - EXP(     - 6.0 * zznd_ml    ) ) )  &
                       &          * ( 1.0 - EXP( - 15.0 * (         1.0 - zznd_ml    ) ) )
                  !
                  ghams(ji,jj,jk) = ghams(ji,jj,jk) + 0.3 * zsc_ws_1(ji,jj)  &
                       &          * ( -2.0 + 2.75 * (       ( 1.0 + 0.6 * zznd_ml**4 )      &
                       &                               - EXP(     - 6.0 * zznd_ml    ) ) )  &
                       &          * ( 1.0 - EXP ( -15.0 * (         1.0 - zznd_ml    ) ) )
               END DO
            ELSE
               DO jk = 2, ibld(ji,jj)
                  zznd_d = gdepw_n(ji,jj,jk) / dstokes(ji,jj)
                  znd = gdepw_n(ji,jj,jk) / zhbl(ji,jj)
                  ghamt(ji,jj,jk) = ghamt(ji,jj,jk) + 0.3 * ( -4.06 * EXP( -2.0 * zznd_d ) * (1.0 - EXP( -4.0 * zznd_d ) ) + &
               &  7.5 * EXP ( -10.0 * ( 0.95 - znd )**2 ) * ( 1.0 - znd ) ) * zsc_wth_1(ji,jj)
                  ghams(ji,jj,jk) = ghams(ji,jj,jk) + 0.3 * ( -4.06 * EXP( -2.0 * zznd_d ) * (1.0 - EXP( -4.0 * zznd_d ) ) + &
               &  7.5 * EXP ( -10.0 * ( 0.95 - znd )**2 ) * ( 1.0 - znd ) ) * zsc_ws_1(ji,jj)
               END DO
            ENDIF
          ENDDO    ! ji loop
       END DO      ! jj loop


       WHERE ( lconv )
          zsc_uw_1 = zustar**2
          zsc_vw_1 = ff_t * zustke * zhml
       ELSEWHERE
          zsc_uw_1 = zustar**2
          zsc_uw_2 = (2.25 - 3.0 * ( 1.0 - EXP( -1.25 * 2.0 ) ) ) * ( 1.0 - EXP( -4.0 * 2.0 ) ) * zsc_uw_1
          zsc_vw_1 = ff_t * zustke * zhbl
          zsc_vw_2 = -0.11 * SIN( 3.14159 * ( 2.0 + 0.4 ) ) * EXP(-( 1.5 + 2.0 )**2 ) * zsc_vw_1
       ENDWHERE

       DO jj = 2, jpjm1
          DO ji = 2, jpim1
             IF ( lconv(ji,jj) ) THEN
               DO jk = 2, imld(ji,jj)
                  zznd_ml = gdepw_n(ji,jj,jk) / zhml(ji,jj)
                  zznd_d = gdepw_n(ji,jj,jk) / dstokes(ji,jj)
                  ghamu(ji,jj,jk) = ghamu(ji,jj,jk)&
                       & + 0.3 * ( -2.0 + 2.5 * ( 1.0 + 0.1 * zznd_ml**4 ) - EXP ( -8.0 * zznd_ml ) ) * zsc_uw_1(ji,jj)
                  !
                  ghamv(ji,jj,jk) = ghamv(ji,jj,jk)&
                       & + 0.3 * 0.1 * ( EXP( -zznd_d ) + EXP( -5.0 * ( 1.0 - zznd_ml ) ) ) * zsc_vw_1(ji,jj)
               END DO
             ELSE
               DO jk = 2, ibld(ji,jj)
                  znd = gdepw_n(ji,jj,jk) / zhbl(ji,jj)
                  zznd_d = gdepw_n(ji,jj,jk) / dstokes(ji,jj)
                  IF ( zznd_d <= 2.0 ) THEN
                     ghamu(ji,jj,jk) = ghamu(ji,jj,jk) + 0.5 * 0.3 &
                          &*  ( 2.25 - 3.0  * ( 1.0 - EXP( - 1.25 * zznd_d ) ) * ( 1.0 - EXP( -2.0 * zznd_d ) ) ) * zsc_uw_1(ji,jj)
                     !
                  ELSE
                     ghamu(ji,jj,jk) = ghamu(ji,jj,jk)&
                          & + 0.5 * 0.3 * ( 1.0 - EXP( -5.0 * ( 1.0 - znd ) ) ) * zsc_uw_2(ji,jj)
                     !
                  ENDIF

                  ghamv(ji,jj,jk) = ghamv(ji,jj,jk)&
                       & + 0.3 * 0.15 * SIN( 3.14159 * ( 0.65 * zznd_d ) ) * EXP( -0.25 * zznd_d**2 ) * zsc_vw_1(ji,jj)
                  ghamv(ji,jj,jk) = ghamv(ji,jj,jk)&
                       & + 0.3 * 0.15 * EXP( -5.0 * ( 1.0 - znd ) ) * ( 1.0 - EXP( -20.0 * ( 1.0 - znd ) ) ) * zsc_vw_2(ji,jj)
               END DO
             ENDIF
          END DO
       END DO
!
! Make surface forced velocity non-gradient terms go to zero at the base of the mixed layer.

      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            IF ( lconv(ji,jj) ) THEN
               DO jk = 2, ibld(ji,jj)
                  znd = ( gdepw_n(ji,jj,jk) - zhml(ji,jj) ) / zhml(ji,jj) !ALMG to think about
                  IF ( znd >= 0.0 ) THEN
                     ghamu(ji,jj,jk) = ghamu(ji,jj,jk) * ( 1.0 - EXP( -30.0 * znd**2 ) )
                     ghamv(ji,jj,jk) = ghamv(ji,jj,jk) * ( 1.0 - EXP( -30.0 * znd**2 ) )
                  ELSE
                     ghamu(ji,jj,jk) = 0._wp
                     ghamv(ji,jj,jk) = 0._wp
                  ENDIF
               END DO
            ELSE
               DO jk = 2, ibld(ji,jj)
                  znd = ( gdepw_n(ji,jj,jk) - zhml(ji,jj) ) / zhml(ji,jj) !ALMG to think about
                  IF ( znd >= 0.0 ) THEN
                     ghamu(ji,jj,jk) = ghamu(ji,jj,jk) * ( 1.0 - EXP( -10.0 * znd**2 ) )
                     ghamv(ji,jj,jk) = ghamv(ji,jj,jk) * ( 1.0 - EXP( -10.0 * znd**2 ) )
                  ELSE
                     ghamu(ji,jj,jk) = 0._wp
                     ghamv(ji,jj,jk) = 0._wp
                  ENDIF
               END DO
            ENDIF
         END DO
      END DO

      ! pynocline contributions
       ! Temporary fix to avoid instabilities when zdb_bl becomes very very small
       zsc_uw_1 = 0._wp ! 50.0 * zla**(8.0/3.0) * zustar**2 * zhbl / ( zdb_bl + epsln )
       DO jj = 2, jpjm1
          DO ji = 2, jpim1
             DO jk= 2, ibld(ji,jj)
                znd = gdepw_n(ji,jj,jk) / zhbl(ji,jj)
                ghamt(ji,jj,jk) = ghamt(ji,jj,jk) + zdiffut(ji,jj,jk) * zdtdz_pyc(ji,jj,jk)
                ghams(ji,jj,jk) = ghams(ji,jj,jk) + zdiffut(ji,jj,jk) * zdsdz_pyc(ji,jj,jk)
                ghamu(ji,jj,jk) = ghamu(ji,jj,jk) + zviscos(ji,jj,jk) * zdudz_pyc(ji,jj,jk)
                ghamu(ji,jj,jk) = ghamu(ji,jj,jk) + zsc_uw_1(ji,jj) * ( 1.0 - znd )**(7.0/4.0) * zdbdz_pyc(ji,jj,jk)
                ghamv(ji,jj,jk) = ghamv(ji,jj,jk) + zviscos(ji,jj,jk) * zdvdz_pyc(ji,jj,jk)
             END DO
           END DO
       END DO

! Entrainment contribution.

       DO jj=2, jpjm1
          DO ji = 2, jpim1
             IF ( lconv(ji,jj) ) THEN
               DO jk = 1, imld(ji,jj) - 1
                  znd=gdepw_n(ji,jj,jk) / zhml(ji,jj)
                  ghamt(ji,jj,jk) = ghamt(ji,jj,jk) + zwth_ent(ji,jj) * znd
                  ghams(ji,jj,jk) = ghams(ji,jj,jk) + zws_ent(ji,jj) * znd
                  ghamu(ji,jj,jk) = ghamu(ji,jj,jk) + zuw_bse(ji,jj) * znd
                  ghamv(ji,jj,jk) = ghamv(ji,jj,jk) + zvw_bse(ji,jj) * znd
               END DO
               DO jk = imld(ji,jj), ibld(ji,jj)
                  znd = -( gdepw_n(ji,jj,jk) - zhml(ji,jj) ) / zdh(ji,jj)
                  ghamt(ji,jj,jk) = ghamt(ji,jj,jk) + zwth_ent(ji,jj) * ( 1.0 + znd )
                  ghams(ji,jj,jk) = ghams(ji,jj,jk) + zws_ent(ji,jj) * ( 1.0 + znd )
                  ghamu(ji,jj,jk) = ghamu(ji,jj,jk) + zuw_bse(ji,jj) * ( 1.0 + znd )
                  ghamv(ji,jj,jk) = ghamv(ji,jj,jk) + zvw_bse(ji,jj) * ( 1.0 + znd )
                END DO
             ENDIF
             ghamt(ji,jj,ibld(ji,jj)) = 0._wp
             ghams(ji,jj,ibld(ji,jj)) = 0._wp
             ghamu(ji,jj,ibld(ji,jj)) = 0._wp
             ghamv(ji,jj,ibld(ji,jj)) = 0._wp
          END DO       ! ji loop
       END DO          ! jj loop


       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       ! Need to put in code for contributions that are applied explicitly to
       ! the prognostic variables
       !  1. Entrainment flux
       !
       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



       ! rotate non-gradient velocity terms back to model reference frame

       DO jj = 2, jpjm1
          DO ji = 2, jpim1
             DO jk = 2, ibld(ji,jj)
                ztemp = ghamu(ji,jj,jk)
                ghamu(ji,jj,jk) = ghamu(ji,jj,jk) * zcos_wind(ji,jj) - ghamv(ji,jj,jk) * zsin_wind(ji,jj)
                ghamv(ji,jj,jk) = ghamv(ji,jj,jk) * zcos_wind(ji,jj) + ztemp * zsin_wind(ji,jj)
             END DO
          END DO
       END DO

       IF(ln_dia_osm) THEN
          IF ( iom_use("zdtdz_pyc") ) CALL iom_put( "zdtdz_pyc", wmask*zdtdz_pyc )
       END IF

! KPP-style Ri# mixing
       IF( ln_kpprimix) THEN
          DO jk = 2, jpkm1           !* Shear production at uw- and vw-points (energy conserving form)
             DO jj = 1, jpjm1
                DO ji = 1, jpim1   ! vector opt.
                   z3du(ji,jj,jk) = 0.5 * (  un(ji,jj,jk-1) -  un(ji  ,jj,jk) )   &
                        &                 * (  ub(ji,jj,jk-1) -  ub(ji  ,jj,jk) ) * wumask(ji,jj,jk) &
                        &                 / (  e3uw_n(ji,jj,jk) * e3uw_b(ji,jj,jk) )
                   z3dv(ji,jj,jk) = 0.5 * (  vn(ji,jj,jk-1) -  vn(ji,jj  ,jk) )   &
                        &                 * (  vb(ji,jj,jk-1) -  vb(ji,jj  ,jk) ) * wvmask(ji,jj,jk) &
                        &                 / (  e3vw_n(ji,jj,jk) * e3vw_b(ji,jj,jk) )
                END DO
             END DO
          END DO
      !
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  !                                          ! shear prod. at w-point weightened by mask
                  zesh2  =  ( z3du(ji-1,jj,jk) + z3du(ji,jj,jk) ) / MAX( 1._wp , umask(ji-1,jj,jk) + umask(ji,jj,jk) )   &
                     &    + ( z3dv(ji,jj-1,jk) + z3dv(ji,jj,jk) ) / MAX( 1._wp , vmask(ji,jj-1,jk) + vmask(ji,jj,jk) )
                  !                                          ! local Richardson number
                  zri   = MAX( rn2b(ji,jj,jk), 0._wp ) / MAX(zesh2, epsln)
                  zfri =  MIN( zri / rn_riinfty , 1.0_wp )
                  zfri  = ( 1.0_wp - zfri * zfri )
                  zrimix(ji,jj,jk)  =  zfri * zfri  * zfri * wmask(ji, jj, jk)
                END DO
             END DO
          END DO

          DO jj = 2, jpjm1
             DO ji = 2, jpim1
                DO jk = ibld(ji,jj) + 1, jpkm1
                   zdiffut(ji,jj,jk) = zrimix(ji,jj,jk)*rn_difri
                   zviscos(ji,jj,jk) = zrimix(ji,jj,jk)*rn_difri
                END DO
             END DO
          END DO

       END IF ! ln_kpprimix = .true.

! KPP-style set diffusivity large if unstable below BL
       IF( ln_convmix) THEN
          DO jj = 2, jpjm1
             DO ji = 2, jpim1
                DO jk = ibld(ji,jj) + 1, jpkm1
                  IF(  MIN( rn2(ji,jj,jk), rn2b(ji,jj,jk) ) <= -1.e-12 ) zdiffut(ji,jj,jk) = rn_difconv
                END DO
             END DO
          END DO
       END IF ! ln_convmix = .true.

       ! Lateral boundary conditions on zvicos (sign unchanged), needed to caclulate viscosities on u and v grids
       CALL lbc_lnk( 'zdfosm', zviscos(:,:,:), 'W', 1. )

       ! GN 25/8: need to change tmask --> wmask

     DO jk = 2, jpkm1
         DO jj = 2, jpjm1
             DO ji = 2, jpim1
                p_avt(ji,jj,jk) = MAX( zdiffut(ji,jj,jk), avtb(jk) ) * tmask(ji,jj,jk)
                p_avm(ji,jj,jk) = MAX( zviscos(ji,jj,jk), avmb(jk) ) * tmask(ji,jj,jk)
             END DO
         END DO
     END DO
      ! Lateral boundary conditions on ghamu and ghamv, currently on W-grid  (sign unchanged), needed to caclulate gham[uv] on u and v grids
     CALL lbc_lnk_multi( 'zdfosm', p_avt, 'W', 1. , p_avm, 'W', 1.,   &
      &                  ghamu, 'W', 1. , ghamv, 'W', 1. )
       DO jk = 2, jpkm1
           DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ghamu(ji,jj,jk) = ( ghamu(ji,jj,jk) + ghamu(ji+1,jj,jk) ) &
                     &  / MAX( 1., tmask(ji,jj,jk) + tmask (ji + 1,jj,jk) ) * umask(ji,jj,jk)

                  ghamv(ji,jj,jk) = ( ghamv(ji,jj,jk) + ghamv(ji,jj+1,jk) ) &
                      &  / MAX( 1., tmask(ji,jj,jk) + tmask (ji,jj+1,jk) ) * vmask(ji,jj,jk)

                  ghamt(ji,jj,jk) =  ghamt(ji,jj,jk) * tmask(ji,jj,jk)
                  ghams(ji,jj,jk) =  ghams(ji,jj,jk) * tmask(ji,jj,jk)
               END DO
           END DO
        END DO
        ! Lateral boundary conditions on final outputs for gham[ts],  on W-grid  (sign unchanged)
        ! Lateral boundary conditions on final outputs for gham[uv],  on [UV]-grid  (sign unchanged)
        CALL lbc_lnk_multi( 'zdfosm', ghamt, 'W', 1. , ghams, 'W', 1.,   &
         &                  ghamu, 'U', 1. , ghamv, 'V', 1. )

       IF(ln_dia_osm) THEN
         SELECT CASE (nn_osm_wave)
         ! Stokes drift set by assumimg onstant La#=0.3(=0)  or Pierson-Moskovitz spectrum (=1).
         CASE(0:1)
            IF ( iom_use("us_x") ) CALL iom_put( "us_x", tmask(:,:,1)*zustke*zcos_wind )   ! x surface Stokes drift
            IF ( iom_use("us_y") ) CALL iom_put( "us_y", tmask(:,:,1)*zustke*zsin_wind )  ! y surface Stokes drift
            IF ( iom_use("wind_wave_abs_power") ) CALL iom_put( "wind_wave_abs_power", 1000.*rau0*tmask(:,:,1)*zustar**2*zustke )
         ! Stokes drift read in from sbcwave  (=2).
         CASE(2)
            IF ( iom_use("us_x") ) CALL iom_put( "us_x", ut0sd )               ! x surface Stokes drift
            IF ( iom_use("us_y") ) CALL iom_put( "us_y", vt0sd )               ! y surface Stokes drift
            IF ( iom_use("wind_wave_abs_power") ) CALL iom_put( "wind_wave_abs_power", 1000.*rau0*tmask(:,:,1)*zustar**2* &
                 & SQRT(ut0sd**2 + vt0sd**2 ) )
         END SELECT
         IF ( iom_use("ghamt") ) CALL iom_put( "ghamt", tmask*ghamt )            ! <Tw_NL>
         IF ( iom_use("ghams") ) CALL iom_put( "ghams", tmask*ghams )            ! <Sw_NL>
         IF ( iom_use("ghamu") ) CALL iom_put( "ghamu", umask*ghamu )            ! <uw_NL>
         IF ( iom_use("ghamv") ) CALL iom_put( "ghamv", vmask*ghamv )            ! <vw_NL>
         IF ( iom_use("zwth0") ) CALL iom_put( "zwth0", tmask(:,:,1)*zwth0 )            ! <Tw_0>
         IF ( iom_use("zws0") ) CALL iom_put( "zws0", tmask(:,:,1)*zws0 )                ! <Sw_0>
         IF ( iom_use("hbl") ) CALL iom_put( "hbl", tmask(:,:,1)*hbl )                  ! boundary-layer depth
         IF ( iom_use("hbli") ) CALL iom_put( "hbli", tmask(:,:,1)*hbli )               ! Initial boundary-layer depth
         IF ( iom_use("dstokes") ) CALL iom_put( "dstokes", tmask(:,:,1)*dstokes )      ! Stokes drift penetration depth
         IF ( iom_use("zustke") ) CALL iom_put( "zustke", tmask(:,:,1)*zustke )            ! Stokes drift magnitude at T-points
         IF ( iom_use("zwstrc") ) CALL iom_put( "zwstrc", tmask(:,:,1)*zwstrc )         ! convective velocity scale
         IF ( iom_use("zwstrl") ) CALL iom_put( "zwstrl", tmask(:,:,1)*zwstrl )         ! Langmuir velocity scale
         IF ( iom_use("zustar") ) CALL iom_put( "zustar", tmask(:,:,1)*zustar )         ! friction velocity scale
         IF ( iom_use("wind_power") ) CALL iom_put( "wind_power", 1000.*rau0*tmask(:,:,1)*zustar**3 ) ! BL depth internal to zdf_osm routine
         IF ( iom_use("wind_wave_power") ) CALL iom_put( "wind_wave_power", 1000.*rau0*tmask(:,:,1)*zustar**2*zustke )
         IF ( iom_use("zhbl") ) CALL iom_put( "zhbl", tmask(:,:,1)*zhbl )               ! BL depth internal to zdf_osm routine
         IF ( iom_use("zhml") ) CALL iom_put( "zhml", tmask(:,:,1)*zhml )               ! ML depth internal to zdf_osm routine
         IF ( iom_use("zdh") ) CALL iom_put( "zdh", tmask(:,:,1)*zdh )               ! ML depth internal to zdf_osm routine
         IF ( iom_use("zhol") ) CALL iom_put( "zhol", tmask(:,:,1)*zhol )               ! ML depth internal to zdf_osm routine
         IF ( iom_use("zwthav") ) CALL iom_put( "zwthav", tmask(:,:,1)*zwthav )               ! ML depth internal to zdf_osm routine
         IF ( iom_use("zwth_ent") ) CALL iom_put( "zwth_ent", tmask(:,:,1)*zwth_ent )               ! ML depth internal to zdf_osm routine
         IF ( iom_use("zt_ml") ) CALL iom_put( "zt_ml", tmask(:,:,1)*zt_ml )               ! average T in ML
      END IF
      ! Lateral boundary conditions on p_avt  (sign unchanged)
      CALL lbc_lnk( 'zdfosm', p_avt(:,:,:), 'W', 1. )
      !
   END SUBROUTINE zdf_osm


   SUBROUTINE zdf_osm_init
     !!----------------------------------------------------------------------
     !!                  ***  ROUTINE zdf_osm_init  ***
     !!
     !! ** Purpose :   Initialization of the vertical eddy diffivity and
     !!      viscosity when using a osm turbulent closure scheme
     !!
     !! ** Method  :   Read the namosm namelist and check the parameters
     !!      called at the first timestep (nit000)
     !!
     !! ** input   :   Namlist namosm
     !!----------------------------------------------------------------------
     INTEGER  ::   ios            ! local integer
     INTEGER  ::   ji, jj, jk     ! dummy loop indices
     !!
     NAMELIST/namzdf_osm/ ln_use_osm_la, rn_osm_la, rn_osm_dstokes, nn_ave &
          & ,nn_osm_wave, ln_dia_osm, rn_osm_hbl0 &
          & ,ln_kpprimix, rn_riinfty, rn_difri, ln_convmix, rn_difconv
     !!----------------------------------------------------------------------
     !
     REWIND( numnam_ref )              ! Namelist namzdf_osm in reference namelist : Osmosis ML model
     READ  ( numnam_ref, namzdf_osm, IOSTAT = ios, ERR = 901)
901  IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_osm in reference namelist' )

     REWIND( numnam_cfg )              ! Namelist namzdf_tke in configuration namelist : Turbulent Kinetic Energy
     READ  ( numnam_cfg, namzdf_osm, IOSTAT = ios, ERR = 902 )
902  IF( ios >  0 ) CALL ctl_nam ( ios , 'namzdf_osm in configuration namelist' )
     IF(lwm) WRITE ( numond, namzdf_osm )

     IF(lwp) THEN                    ! Control print
        WRITE(numout,*)
        WRITE(numout,*) 'zdf_osm_init : OSMOSIS Parameterisation'
        WRITE(numout,*) '~~~~~~~~~~~~'
        WRITE(numout,*) '   Namelist namzdf_osm : set tke mixing parameters'
        WRITE(numout,*) '     Use namelist  rn_osm_la                     ln_use_osm_la = ', ln_use_osm_la
        WRITE(numout,*) '     Turbulent Langmuir number                     rn_osm_la   = ', rn_osm_la
        WRITE(numout,*) '     Initial hbl for 1D runs                       rn_osm_hbl0   = ', rn_osm_hbl0
        WRITE(numout,*) '     Depth scale of Stokes drift                rn_osm_dstokes = ', rn_osm_dstokes
        WRITE(numout,*) '     horizontal average flag                       nn_ave      = ', nn_ave
        WRITE(numout,*) '     Stokes drift                                  nn_osm_wave = ', nn_osm_wave
        SELECT CASE (nn_osm_wave)
        CASE(0)
           WRITE(numout,*) '     calculated assuming constant La#=0.3'
        CASE(1)
           WRITE(numout,*) '     calculated from Pierson Moskowitz wind-waves'
        CASE(2)
           WRITE(numout,*) '     calculated from ECMWF wave fields'
        END SELECT
        WRITE(numout,*) '     Output osm diagnostics                       ln_dia_osm  = ',  ln_dia_osm
        WRITE(numout,*) '     Use KPP-style shear instability mixing       ln_kpprimix = ', ln_kpprimix
        WRITE(numout,*) '     local Richardson Number limit for shear instability rn_riinfty = ', rn_riinfty
        WRITE(numout,*) '     maximum shear diffusivity at Rig = 0    (m2/s) rn_difri = ', rn_difri
        WRITE(numout,*) '     Use large mixing below BL when unstable       ln_convmix = ', ln_convmix
        WRITE(numout,*) '     diffusivity when unstable below BL     (m2/s) rn_difconv = ', rn_difconv
     ENDIF


     !                              ! Check wave coupling settings !
     !                         ! Further work needed - see ticket #2447 !
     IF( nn_osm_wave == 2 ) THEN
        IF (.NOT. ( ln_wave .AND. ln_sdw )) &
           & CALL ctl_stop( 'zdf_osm_init : ln_zdfosm and nn_osm_wave=2, ln_wave and ln_sdw must be true' )
     END IF

     !                              ! allocate zdfosm arrays
     IF( zdf_osm_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'zdf_osm_init : unable to allocate arrays' )

     call osm_rst( nit000, 'READ' ) !* read or initialize hbl

     IF( ln_zdfddm) THEN
        IF(lwp) THEN
           WRITE(numout,*)
           WRITE(numout,*) '    Double diffusion mixing on temperature and salinity '
           WRITE(numout,*) '    CAUTION : done in routine zdfosm, not in routine zdfddm '
        ENDIF
     ENDIF


     !set constants not in namelist
     !-----------------------------

     IF(lwp) THEN
        WRITE(numout,*)
     ENDIF

     IF (nn_osm_wave == 0) THEN
        dstokes(:,:) = rn_osm_dstokes
     END IF

     ! Horizontal average : initialization of weighting arrays
     ! -------------------

     SELECT CASE ( nn_ave )

     CASE ( 0 )                ! no horizontal average
        IF(lwp) WRITE(numout,*) '          no horizontal average on avt'
        IF(lwp) WRITE(numout,*) '          only in very high horizontal resolution !'
        ! weighting mean arrays etmean
        !           ( 1  1 )
        ! avt = 1/4 ( 1  1 )
        !
        etmean(:,:,:) = 0.e0

        DO jk = 1, jpkm1
           DO jj = 2, jpjm1
              DO ji = 2, jpim1   ! vector opt.
                 etmean(ji,jj,jk) = tmask(ji,jj,jk)                     &
                      &  / MAX( 1.,  umask(ji-1,jj  ,jk) + umask(ji,jj,jk)   &
                      &            + vmask(ji  ,jj-1,jk) + vmask(ji,jj,jk)  )
              END DO
           END DO
        END DO

     CASE ( 1 )                ! horizontal average
        IF(lwp) WRITE(numout,*) '          horizontal average on avt'
        ! weighting mean arrays etmean
        !           ( 1/2  1  1/2 )
        ! avt = 1/8 ( 1    2  1   )
        !           ( 1/2  1  1/2 )
        etmean(:,:,:) = 0.e0

        DO jk = 1, jpkm1
           DO jj = 2, jpjm1
              DO ji = 2, jpim1   ! vector opt.
                 etmean(ji,jj,jk) = tmask(ji, jj,jk)                           &
                      & / MAX( 1., 2.* tmask(ji,jj,jk)                           &
                      &      +.5 * ( tmask(ji-1,jj+1,jk) + tmask(ji-1,jj-1,jk)   &
                      &             +tmask(ji+1,jj+1,jk) + tmask(ji+1,jj-1,jk) ) &
                      &      +1. * ( tmask(ji-1,jj  ,jk) + tmask(ji  ,jj+1,jk)   &
                      &             +tmask(ji  ,jj-1,jk) + tmask(ji+1,jj  ,jk) ) )
              END DO
           END DO
        END DO

     CASE DEFAULT
        WRITE(ctmp1,*) '          bad flag value for nn_ave = ', nn_ave
        CALL ctl_stop( ctmp1 )

     END SELECT

     ! Initialization of vertical eddy coef. to the background value
     ! -------------------------------------------------------------
     DO jk = 1, jpk
        avt (:,:,jk) = avtb(jk) * tmask(:,:,jk)
     END DO

     ! zero the surface flux for non local term and osm mixed layer depth
     ! ------------------------------------------------------------------
     ghamt(:,:,:) = 0.
     ghams(:,:,:) = 0.
     ghamu(:,:,:) = 0.
     ghamv(:,:,:) = 0.
     !
     IF( lwxios ) THEN
        CALL iom_set_rstw_var_active('wn')
        CALL iom_set_rstw_var_active('hbl')
        CALL iom_set_rstw_var_active('hbli')
     ENDIF
   END SUBROUTINE zdf_osm_init


   SUBROUTINE osm_rst( kt, cdrw )
     !!---------------------------------------------------------------------
     !!                   ***  ROUTINE osm_rst  ***
     !!
     !! ** Purpose :   Read or write BL fields in restart file
     !!
     !! ** Method  :   use of IOM library. If the restart does not contain
     !!                required fields, they are recomputed from stratification
     !!----------------------------------------------------------------------

     INTEGER, INTENT(in) :: kt
     CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag

     INTEGER ::   id1, id2   ! iom enquiry index
     INTEGER  ::   ji, jj, jk     ! dummy loop indices
     INTEGER  ::   iiki, ikt ! local integer
     REAL(wp) ::   zhbf           ! tempory scalars
     REAL(wp) ::   zN2_c           ! local scalar
     REAL(wp) ::   rho_c = 0.01_wp    !: density criterion for mixed layer depth
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: imld_rst ! level of mixed-layer depth (pycnocline top)
     !!----------------------------------------------------------------------
     !
     !!-----------------------------------------------------------------------------
     ! If READ/WRITE Flag is 'READ', try to get hbl from restart file. If successful then return
     !!-----------------------------------------------------------------------------
     IF( TRIM(cdrw) == 'READ'.AND. ln_rstart) THEN
        id1 = iom_varid( numror, 'wn'   , ldstop = .FALSE. )
        IF( id1 > 0 ) THEN                       ! 'wn' exists; read
           CALL iom_get( numror, jpdom_autoglo, 'wn', wn, ldxios = lrxios )
           WRITE(numout,*) ' ===>>>> :  wn read from restart file'
        ELSE
           wn(:,:,:) = 0._wp
           WRITE(numout,*) ' ===>>>> :  wn not in restart file, set to zero initially'
        END IF
        id1 = iom_varid( numror, 'hbl'   , ldstop = .FALSE. )
        id2 = iom_varid( numror, 'hbli'   , ldstop = .FALSE. )
        IF( id1 > 0 .AND. id2 > 0) THEN                       ! 'hbl' exists; read and return
           CALL iom_get( numror, jpdom_autoglo, 'hbl' , hbl , ldxios = lrxios )
           CALL iom_get( numror, jpdom_autoglo, 'hbli', hbli, ldxios = lrxios  )
           WRITE(numout,*) ' ===>>>> :  hbl & hbli read from restart file'
           RETURN
        ELSE                      ! 'hbl' & 'hbli' not in restart file, recalculate
           WRITE(numout,*) ' ===>>>> : previous run without osmosis scheme, hbl computed from stratification'
        END IF
     END IF

     !!-----------------------------------------------------------------------------
     ! If READ/WRITE Flag is 'WRITE', write hbl into the restart file, then return
     !!-----------------------------------------------------------------------------
     IF( TRIM(cdrw) == 'WRITE') THEN     !* Write hbli into the restart file, then return
        IF(lwp) WRITE(numout,*) '---- osm-rst ----'
         CALL iom_rstput( kt, nitrst, numrow, 'wn'     , wn  , ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'hbl'    , hbl , ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'hbli'   , hbli, ldxios = lwxios )
        RETURN
     END IF

     !!-----------------------------------------------------------------------------
     ! Getting hbl, no restart file with hbl, so calculate from surface stratification
     !!-----------------------------------------------------------------------------
     IF( lwp ) WRITE(numout,*) ' ===>>>> : calculating hbl computed from stratification'
     ALLOCATE( imld_rst(jpi,jpj) )
     ! w-level of the mixing and mixed layers
     CALL eos_rab( tsn, rab_n )
     CALL bn2(tsn, rab_n, rn2)
     imld_rst(:,:)  = nlb10         ! Initialization to the number of w ocean point
     hbl(:,:)  = 0._wp              ! here hbl used as a dummy variable, integrating vertically N^2
     zN2_c = grav * rho_c * r1_rau0 ! convert density criteria into N^2 criteria
     !
     hbl(:,:)  = 0._wp              ! here hbl used as a dummy variable, integrating vertically N^2
     DO jk = 1, jpkm1
        DO jj = 1, jpj              ! Mixed layer level: w-level
           DO ji = 1, jpi
              ikt = mbkt(ji,jj)
              hbl(ji,jj) = hbl(ji,jj) + MAX( rn2(ji,jj,jk) , 0._wp ) * e3w_n(ji,jj,jk)
              IF( hbl(ji,jj) < zN2_c )   imld_rst(ji,jj) = MIN( jk , ikt ) + 1   ! Mixed layer level
           END DO
        END DO
     END DO
     !
     DO jj = 1, jpj
        DO ji = 1, jpi
           iiki = imld_rst(ji,jj)
           hbl (ji,jj) = gdepw_n(ji,jj,iiki  ) * ssmask(ji,jj)    ! Turbocline depth
        END DO
     END DO
     hbl = MAX(hbl,epsln)
     hbli(:,:) = hbl(:,:)
     DEALLOCATE( imld_rst )
     WRITE(numout,*) ' ===>>>> : hbl computed from stratification'
   END SUBROUTINE osm_rst


   SUBROUTINE tra_osm( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_osm  ***
      !!
      !! ** Purpose :   compute and add to the tracer trend the non-local tracer flux
      !!
      !! ** Method  :   ???
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdt, ztrds   ! 3D workspace
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt
      INTEGER :: ji, jj, jk
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_osm : OSM non-local tracer fluxes'
         IF(lwp) WRITE(numout,*) '~~~~~~~   '
      ENDIF

      IF( l_trdtra )   THEN                    !* Save ta and sa trends
         ALLOCATE( ztrdt(jpi,jpj,jpk) )   ;    ztrdt(:,:,:) = tsa(:,:,:,jp_tem)
         ALLOCATE( ztrds(jpi,jpj,jpk) )   ;    ztrds(:,:,:) = tsa(:,:,:,jp_sal)
      ENDIF

      ! add non-local temperature and salinity flux
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               tsa(ji,jj,jk,jp_tem) =  tsa(ji,jj,jk,jp_tem)                      &
                  &                 - (  ghamt(ji,jj,jk  )  &
                  &                    - ghamt(ji,jj,jk+1) ) /e3t_n(ji,jj,jk)
               tsa(ji,jj,jk,jp_sal) =  tsa(ji,jj,jk,jp_sal)                      &
                  &                 - (  ghams(ji,jj,jk  )  &
                  &                    - ghams(ji,jj,jk+1) ) / e3t_n(ji,jj,jk)
            END DO
         END DO
      END DO


      ! save the non-local tracer flux trends for diagnostic
      IF( l_trdtra )   THEN
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem) - ztrdt(:,:,:)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal) - ztrds(:,:,:)
!!bug gm jpttdzdf ==> jpttosm
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_zdf, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_zdf, ztrds )
         DEALLOCATE( ztrdt )      ;     DEALLOCATE( ztrds )
      ENDIF

      IF(ln_ctl) THEN
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' osm  - Ta: ', mask1=tmask,   &
         &             tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      ENDIF
      !
   END SUBROUTINE tra_osm


   SUBROUTINE trc_osm( kt )          ! Dummy routine
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_osm  ***
      !!
      !! ** Purpose :   compute and add to the passive tracer trend the non-local
      !!                 passive tracer flux
      !!
      !!
      !! ** Method  :   ???
      !!----------------------------------------------------------------------
      !
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_osm: Not written yet', kt
   END SUBROUTINE trc_osm


   SUBROUTINE dyn_osm( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_osm  ***
      !!
      !! ** Purpose :   compute and add to the velocity trend the non-local flux
      !! copied/modified from tra_osm
      !!
      !! ** Method  :   ???
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   !
      !
      INTEGER :: ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_osm : OSM non-local velocity'
         IF(lwp) WRITE(numout,*) '~~~~~~~   '
      ENDIF
      !code saving tracer trends removed, replace with trdmxl_oce

      DO jk = 1, jpkm1           ! add non-local u and v fluxes
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ua(ji,jj,jk) =  ua(ji,jj,jk)                      &
                  &                 - (  ghamu(ji,jj,jk  )  &
                  &                    - ghamu(ji,jj,jk+1) ) / e3u_n(ji,jj,jk)
               va(ji,jj,jk) =  va(ji,jj,jk)                      &
                  &                 - (  ghamv(ji,jj,jk  )  &
                  &                    - ghamv(ji,jj,jk+1) ) / e3v_n(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      ! code for saving tracer trends removed
      !
   END SUBROUTINE dyn_osm

   !!======================================================================
END MODULE zdfosm
