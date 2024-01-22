MODULE sbcblk_algo_ncar
   !!======================================================================
   !!                   ***  MODULE  sbcblk_algo_ncar  ***
   !! Computes:
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m U_blk
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of Large & Yeager 2008
   !!
   !!       Routine turb_ncar maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!                         L. Brodeau, 2020
   !!=====================================================================
   !! History :  4.0  !  2020-06  (L.Brodeau) successor of old turb_ncar of former sbcblk_core.F90
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   turb_ncar  : computes the bulk turbulent transfer coefficients
   !!                   adjusts t_air and q_air from zt to zu m
   !!                   returns the effective bulk wind speed at 10m
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbcwave, ONLY   :  cdn_wave ! wave module
#if defined key_si3 || defined key_cice
   USE sbc_ice         ! Surface boundary condition: ice fields
#endif
   !
   USE iom             ! I/O manager library
   USE lib_mpp         ! distribued memory computing library
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE lib_fortran     ! to use key_nosignedzero


   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TURB_NCAR   ! called by sbcblk.F90

   !                              ! NCAR own values for given constants:
   REAL(wp), PARAMETER ::   rctv0 = 0.608   ! constant to obtain virtual temperature...

   INTEGER , PARAMETER ::   nb_itt = 4       ! number of itterations

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ncar( zt, zu, sst, t_zt, ssq, q_zt, U_zu, &
      &                  Cd, Ch, Ce, t_zu, q_zu, Ub,         &
      &                  CdN, ChN, CeN                       )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ncar  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Large & Yeager (2004) and Large & Yeager (2008)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at zu to be used in the bulk formulas
      !!
      !! INPUT :
      !! -------
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  sst  : bulk SST                                                [K]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  ssq  : specific humidity at saturation at SST                  [kg/kg]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd     : drag coefficient
      !!    *  Ch     : sensible heat coefficient
      !!    *  Ce     : evaporation coefficient
      !!    *  t_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  q_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  Ub  : bulk wind speed at zu                                 [m/s]
      !!
      !! ** Author: L. Brodeau, June 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in   )                     ::   zt       ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for U_zu                             [m]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   sst      ! sea surface temperature                [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   ssq      ! sea surface specific humidity           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity at zt             [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ub    ! bulk wind speed at zu                     [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   CdN, ChN, CeN ! neutral transfer coefficients
      !
      INTEGER :: j_itt
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   Cx_n10        ! 10m neutral latent/sensible coefficient
      REAL(wp), DIMENSION(jpi,jpj) ::   sqrtCdN10   ! root square of Cd_n10
      REAL(wp), DIMENSION(jpi,jpj) ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(jpi,jpj) ::   zpsi_h_u
      REAL(wp), DIMENSION(jpi,jpj) ::   ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(jpi,jpj) ::   sqrtCd
      !!----------------------------------------------------------------------------------

      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp )

      Ub = MAX( 0.5_wp , U_zu )   !  relative wind speed at zu (normally 10m), we don't want to fall under 0.5 m/s

      !! Neutral drag coefficient at zu:
      IF( ln_cdgw ) THEN      ! wave drag case
         CdN = MAX( cdn_wave(:,:) + rsmall * ( 1._wp - tmask(:,:,1) ) , 0.1E-3_wp )
      ELSE
         CdN = CD_N10_NCAR( Ub )
      ENDIF
      sqrtCdN10 = SQRT( CdN )

      !! Initializing transf. coeff. with their first guess neutral equivalents :
      Cd = CdN
      Ce = CE_N10_NCAR( sqrtCdN10 )
      ztmp0 = 0.5_wp + SIGN(0.5_wp, virt_temp(t_zt, q_zt) - virt_temp(sst, ssq)) ! we guess stability based on delta of virt. pot. temp.
      Ch = CH_N10_NCAR( sqrtCdN10 , ztmp0 )
      sqrtCd = sqrtCdN10

      !! Initializing values at z_u with z_t values:
      t_zu = t_zt
      q_zu = q_zt

      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt
         !
         ztmp1 = t_zu - sst   ! Updating air/sea differences
         ztmp2 = q_zu - ssq

         ! Updating turbulent scales :   (L&Y 2004 Eq. (7))
         ztmp0 = sqrtCd*Ub          ! u*
         ztmp1 = Ch/sqrtCd*ztmp1    ! theta*
         ztmp2 = Ce/sqrtCd*ztmp2    ! q*

         ! Estimate the inverse of Obukov length (1/L) at height zu:
         ztmp0 = One_on_L( t_zu, q_zu, ztmp0, ztmp1, ztmp2 )

         !! Stability parameters :
         zeta_u   = zu*ztmp0
         zeta_u   = sign( min(abs(zeta_u),10._wp), zeta_u )

         !! Shifting temperature and humidity at zu (L&Y 2004 Eq. (9b-9c))
         IF( .NOT. l_zt_equal_zu ) THEN
            ztmp0 = zt*ztmp0 ! zeta_t !
            ztmp0 = SIGN( MIN(ABS(ztmp0),10._wp), ztmp0 )  ! Temporaty array ztmp0 == zeta_t !!!
            ztmp0 = LOG(zt/zu) + psi_h_ncar(zeta_u) - psi_h_ncar(ztmp0)                   ! ztmp0 just used as temp array again!
            t_zu = t_zt - ztmp1/vkarmn*ztmp0    ! ztmp1 is still theta*  L&Y 2004 Eq. (9b)
            !!
            q_zu = q_zt - ztmp2/vkarmn*ztmp0    ! ztmp2 is still q*      L&Y 2004 Eq. (9c)
            q_zu = MAX(0._wp, q_zu)
         END IF

         ! Update neutral wind speed at 10m and neutral Cd at 10m (L&Y 2004 Eq. 9a)...
         !   In very rare low-wind conditions, the old way of estimating the
         !   neutral wind speed at 10m leads to a negative value that causes the code
         !   to crash. To prevent this a threshold of 0.25m/s is imposed.
         ztmp2 = psi_m_ncar(zeta_u)
         ztmp0 = MAX( 0.25_wp , UN10_from_CD(zu, Ub, Cd, ppsi=ztmp2) ) ! U_n10 (ztmp2 == psi_m_ncar(zeta_u))

         IF( ln_cdgw ) THEN      ! wave drag case
            CdN = MAX( cdn_wave(:,:) + rsmall * ( 1._wp - tmask(:,:,1) ) , 0.1E-3_wp )
         ELSE
            CdN   = CD_N10_NCAR(ztmp0)                                       ! Cd_n10
         END IF
         sqrtCdN10 = SQRT(CdN)

         !! Update of transfer coefficients:
         ztmp1  = 1._wp + sqrtCdN10/vkarmn*(LOG(zu/10._wp) - ztmp2)   ! L&Y 2004 Eq. (10a) (ztmp2 == psi_m(zeta_u))
         Cd     = MAX( CdN / ( ztmp1*ztmp1 ) , 0.1E-3_wp )
         sqrtCd = SQRT( Cd )

         ztmp0  = ( LOG(zu/10._wp) - psi_h_ncar(zeta_u) ) / vkarmn / sqrtCdN10
         ztmp2  = sqrtCd / sqrtCdN10

         ztmp1  = 0.5_wp + sign(0.5_wp,zeta_u) ! stability flag
         ChN    = CH_N10_NCAR( sqrtCdN10 , ztmp1 )
         ztmp1  = 1._wp + ChN*ztmp0
         Ch     = MAX( ChN*ztmp2 / ztmp1 , 0.1E-3_wp )   ! L&Y 2004 Eq. (10b)

         CeN = CE_N10_NCAR( sqrtCdN10 )
         ztmp1  = 1._wp + CeN*ztmp0
         Ce     = MAX( CeN*ztmp2 / ztmp1 , 0.1E-3_wp )  ! L&Y 2004 Eq. (10c)

      END DO !DO j_itt = 1, nb_itt

   END SUBROUTINE turb_ncar


   FUNCTION CD_N10_NCAR( pw10 )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral drag coefficient at 10m as a function
      !! of neutral wind  speed at 10m
      !!
      !! Origin: Large & Yeager 2008, Eq. (11)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pw10           ! scalar wind speed at 10m (m/s)
      REAL(wp), DIMENSION(jpi,jpj)             :: CD_N10_NCAR
      !
      INTEGER  ::     ji, jj     ! dummy loop indices
      REAL(wp) :: zgt33, zw, zw6 ! local scalars
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zw  = pw10(ji,jj)
            zw6 = zw*zw*zw
            zw6 = zw6*zw6
            !
            ! When wind speed > 33 m/s => Cyclone conditions => special treatment
            zgt33 = 0.5_wp + SIGN( 0.5_wp, (zw - 33._wp) )   ! If pw10 < 33. => 0, else => 1
            !
            CD_N10_NCAR(ji,jj) = 1.e-3_wp * ( &
               &       (1._wp - zgt33)*( 2.7_wp/zw + 0.142_wp + zw/13.09_wp - 3.14807E-10_wp*zw6) & ! wind <  33 m/s
               &      +    zgt33   *      2.34_wp )                                                 ! wind >= 33 m/s
            !
            CD_N10_NCAR(ji,jj) = MAX( CD_N10_NCAR(ji,jj), 0.1E-3_wp )
            !
         END DO
      END DO
      !
   END FUNCTION CD_N10_NCAR



   FUNCTION CH_N10_NCAR( psqrtcdn10 , pstab )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral heat transfer coefficient at 10m      !!
      !! Origin: Large & Yeager 2008, Eq. (9) and (12)

      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: ch_n10_ncar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: psqrtcdn10 ! sqrt( CdN10 )
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pstab      ! stable ABL => 1 / unstable ABL => 0
      !!----------------------------------------------------------------------------------
      !
      ch_n10_ncar = MAX( 1.e-3_wp * psqrtcdn10*( 18._wp*pstab + 32.7_wp*(1._wp - pstab) )  , 0.1E-3_wp )   ! Eq. (9) & (12) Large & Yeager, 2008
      !
   END FUNCTION CH_N10_NCAR

   FUNCTION CE_N10_NCAR( psqrtcdn10 )
      !!----------------------------------------------------------------------------------
      !! Estimate of the neutral heat transfer coefficient at 10m      !!
      !! Origin: Large & Yeager 2008, Eq. (9) and (13)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: ce_n10_ncar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: psqrtcdn10 ! sqrt( CdN10 )
      !!----------------------------------------------------------------------------------
      ce_n10_ncar = MAX( 1.e-3_wp * ( 34.6_wp * psqrtcdn10 ) , 0.1E-3_wp )
      !
   END FUNCTION CE_N10_NCAR


   FUNCTION psi_m_ncar( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for momentum
      !!    !! Psis, L&Y 2004, Eq. (8c), (8d), (8e)
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_m_ncar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zzeta, zx2, zx, zpsi_unst, zpsi_stab,  zstab   ! local scalars
      !!----------------------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi

            zzeta = pzeta(ji,jj)
            !
            zx2 = SQRT( ABS(1._wp - 16._wp*zzeta) )  ! (1 - 16z)^0.5
            zx2 = MAX( zx2 , 1._wp )
            zx  = SQRT(zx2)                          ! (1 - 16z)^0.25
            zpsi_unst = 2._wp*LOG( (1._wp + zx )*0.5_wp )   &
               &            + LOG( (1._wp + zx2)*0.5_wp )   &
               &          - 2._wp*ATAN(zx) + rpi*0.5_wp
            !
            zpsi_stab = -5._wp*zzeta
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zzeta) ! zzeta > 0 => zstab = 1
            !
            psi_m_ncar(ji,jj) =          zstab  * zpsi_stab &  ! (zzeta > 0) Stable
               &              + (1._wp - zstab) * zpsi_unst    ! (zzeta < 0) Unstable
            !
         END DO
      END DO
   END FUNCTION psi_m_ncar


   FUNCTION psi_h_ncar( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !!    !! Psis, L&Y 2004, Eq. (8c), (8d), (8e)
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_h_ncar
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zzeta, zx2, zpsi_unst, zpsi_stab, zstab  ! local scalars
      !!----------------------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zzeta = pzeta(ji,jj)
            !
            zx2 = SQRT( ABS(1._wp - 16._wp*zzeta) )  ! (1 -16z)^0.5
            zx2 = MAX( zx2 , 1._wp )
            zpsi_unst = 2._wp*LOG( 0.5_wp*(1._wp + zx2) )
            !
            zpsi_stab = -5._wp*zzeta
            !
            zstab = 0.5_wp + SIGN(0.5_wp, zzeta) ! zzeta > 0 => zstab = 1
            !
            psi_h_ncar(ji,jj) =          zstab  * zpsi_stab &  ! (zzeta > 0) Stable
               &              + (1._wp - zstab) * zpsi_unst    ! (zzeta < 0) Unstable
            !
         END DO
      END DO
   END FUNCTION psi_h_ncar




   FUNCTION UN10_from_CD( pzu, pUb, pCd, ppsi )
      !!----------------------------------------------------------------------------------
      !!  Provides the neutral-stability wind speed at 10 m
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: UN10_from_CD  !: [m/s]
      REAL(wp),                     INTENT(in) :: pzu  !: measurement heigh of bulk wind speed
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pUb  !: bulk wind speed at height pzu m   [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pCd  !: drag coefficient
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ppsi !: "Psi_m(pzu/L)" stability correction profile for momentum []
      !!----------------------------------------------------------------------------------
      !! Reminder: UN10 = u*/vkarmn * log(10/z0)
      !!     and: u* = sqrt(Cd) * Ub
      !!                                  u*/vkarmn * log(   10   /       z0    )
      UN10_from_CD(:,:) = SQRT(pCd(:,:))*pUb/vkarmn * LOG( 10._wp / z0_from_Cd( pzu, pCd(:,:), ppsi=ppsi(:,:) ) )
      !!
   END FUNCTION UN10_from_CD


   FUNCTION One_on_L( ptha, pqa, pus, pts, pqs )
      !!------------------------------------------------------------------------
      !!
      !! Evaluates the 1./(Obukhov length) from air temperature,
      !! air specific humidity, and frictional scales u*, t* and q*
      !!
      !! Author: L. Brodeau, June 2019 / AeroBulk
      !!         (https://github.com/brodeau/aerobulk/)
      !!------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: One_on_L     !: 1./(Obukhov length) [m^-1]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ptha         !: reference potential temperature of air [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa          !: reference specific humidity of air   [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pus          !: u*: friction velocity [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pts, pqs     !: \theta* and q* friction aka turb. scales for temp. and spec. hum.
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) ::     zqa          ! local scalar
      !!-------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zqa = (1._wp + rctv0*pqa(ji,jj))
            !
            ! The main concern is to know whether, the vertical turbulent flux of virtual temperature, < u' theta_v' > is estimated with:
            !  a/  -u* [ theta* (1 + 0.61 q) + 0.61 theta q* ] => this is the one that seems correct! chose this one!
            !                      or
            !  b/  -u* [ theta*              + 0.61 theta q* ]
            !
            One_on_L(ji,jj) = grav*vkarmn*( pts(ji,jj)*zqa + rctv0*ptha(ji,jj)*pqs(ji,jj) ) &
               &               / MAX( pus(ji,jj)*pus(ji,jj)*ptha(ji,jj)*zqa , 1.E-9_wp )
            !
         END DO
      END DO
      !
      One_on_L = SIGN( MIN(ABS(One_on_L),200._wp), One_on_L ) ! (prevent FPE from stupid values over masked regions...)
      !
   END FUNCTION One_on_L


   FUNCTION z0_from_Cd( pzu, pCd,  ppsi )
      REAL(wp), DIMENSION(jpi,jpj) :: z0_from_Cd        !: roughness length [m]
      REAL(wp)                    , INTENT(in) :: pzu   !: reference height zu [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pCd   !: (neutral or non-neutral) drag coefficient []
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in), OPTIONAL :: ppsi !: "Psi_m(pzu/L)" stability correction profile for momentum []
      !!
      !! If pCd is the NEUTRAL-STABILITY drag coefficient then ppsi must be 0 or not given
      !! If pCd is the drag coefficient (in stable or unstable conditions) then pssi must be provided
      !!----------------------------------------------------------------------------------
      IF ( PRESENT(ppsi) ) THEN
         !! Cd provided is the actual Cd (not the neutral-stability CdN) :
         z0_from_Cd = pzu * EXP( - ( vkarmn/SQRT(pCd(:,:)) + ppsi(:,:) ) ) !LB: ok, double-checked!
      ELSE
         !! Cd provided is the neutral-stability Cd, aka CdN :
         z0_from_Cd = pzu * EXP( - vkarmn/SQRT(pCd(:,:)) )            !LB: ok, double-checked!
      END IF
   END FUNCTION z0_from_Cd

   FUNCTION virt_temp( pta, pqa )
      REAL(wp), DIMENSION(jpi,jpj)             :: virt_temp !: virtual temperature [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pta !: absolute or potential air temperature [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa !: specific humidity of air   [kg/kg]
      virt_temp(:,:) = pta(:,:) * (1._wp + rctv0*pqa(:,:))
   END FUNCTION virt_temp

   !!======================================================================
END MODULE sbcblk_algo_ncar
