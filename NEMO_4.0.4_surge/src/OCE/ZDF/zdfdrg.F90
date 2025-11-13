MODULE zdfdrg
   !!======================================================================
   !!                       ***  MODULE  zdfdrg  ***
   !! Ocean physics: top and/or Bottom friction 
   !!======================================================================
   !! History :  OPA  ! 1997-06  (G. Madec, A.-M. Treguier)  Original code
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            3.2  ! 2009-09  (A.C.Coward)  Correction to include barotropic contribution
   !!            3.3  ! 2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.4  ! 2011-11  (H. Liu) implementation of semi-implicit bottom friction option
   !!                 ! 2012-06  (H. Liu) implementation of Log Layer bottom friction option
   !!            4.0  ! 2017-05  (G. Madec) zdfbfr becomes zdfdrg + variable names change
   !!                                     + drag defined at t-point + new user interface + top drag (ocean cavities)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_drg       : update bottom friction coefficient (non-linear bottom friction only)
   !!   zdf_drg_exp   : compute the top & bottom friction in explicit case
   !!   zdf_drg_init  : read in namdrg namelist and control the bottom friction parameters.
   !!       drg_init  :
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE phycst  , ONLY : vkarmn
   USE dom_oce        ! ocean space and time domain variables
   USE zdf_oce        ! ocean vertical physics variables
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O module
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing
   USE prtctl         ! Print control
   USE sbc_oce , ONLY : nn_ice 

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_drg         ! called by zdf_phy
   PUBLIC   zdf_drg_exp     ! called by dyn_zdf
   PUBLIC   zdf_drg_init    ! called by zdf_phy_init

   !                                 !!* Namelist namdrg: nature of drag coefficient namelist *
   LOGICAL , PUBLIC ::   ln_drg_OFF   ! free-slip       : Cd = 0
   LOGICAL          ::   ln_lin       !     linear  drag: Cd = Cd0_lin
   LOGICAL          ::   ln_non_lin   ! non-linear  drag: Cd = Cd0_nl |U|
   LOGICAL          ::   ln_loglayer  ! logarithmic drag: Cd = vkarmn/log(z/z0)
   LOGICAL , PUBLIC ::   ln_drgimp    ! implicit top/bottom friction flag
   LOGICAL , PUBLIC ::   ln_drgice_imp ! implicit ice-ocean drag 
   !                                 !!* Namelist namdrg_top & _bot: TOP or BOTTOM coefficient namelist *
   REAL(wp)         ::   rn_Cd0       !: drag coefficient                                           [ - ]
   REAL(wp)         ::   rn_Uc0       !: characteristic velocity (linear case: tau=rho*Cd0*Uc0*u)   [m/s]
   REAL(wp)         ::   rn_Cdmax     !: drag value maximum (ln_loglayer=T)                         [ - ]
   REAL(wp)         ::   rn_z0        !: roughness          (ln_loglayer=T)                         [ m ]
   REAL(wp)         ::   rn_ke0       !: background kinetic energy (non-linear case)                [m2/s2]
   LOGICAL          ::   ln_boost     !: =T regional boost of Cd0 ; =F Cd0 horizontally uniform
   REAL(wp)         ::   rn_boost     !: local boost factor                                         [ - ]

   REAL(wp), PUBLIC ::   r_Cdmin_top, r_Cdmax_top, r_z0_top, r_ke0_top   ! set from namdrg_top namelist values
   REAL(wp), PUBLIC ::   r_Cdmin_bot, r_Cdmax_bot, r_z0_bot, r_ke0_bot   !  -    -  namdrg_bot    -       -

   INTEGER ::              ndrg       ! choice of the type of drag coefficient
   !                                  ! associated indices:
   INTEGER, PARAMETER ::   np_OFF      = 0   ! free-slip: drag set to zero
   INTEGER, PARAMETER ::   np_lin      = 1   !     linear drag: Cd = Cd0_lin
   INTEGER, PARAMETER ::   np_non_lin  = 2   ! non-linear drag: Cd = Cd0_nl |U|
   INTEGER, PARAMETER ::   np_loglayer = 3   ! non linear drag (logarithmic formulation): Cd = vkarmn/log(z/z0)

   LOGICAL , PUBLIC ::   l_zdfdrg           !: flag to update at each time step the top/bottom Cd
   LOGICAL          ::   l_log_not_linssh   !: flag to update at each time step the position ot the velocity point 
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC ::   rCd0_top, rCd0_bot   !: precomputed top/bottom drag coeff. at t-point (>0)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC ::   rCdU_top, rCdU_bot   !: top/bottom drag coeff. at t-point (<0)  [m/s]

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE zdf_drg( kt, k_mk, pCdmin, pCdmax, pz0, pke0, pCd0,   &   ! <<== in 
      &                                                     pCdU )      ! ==>> out : bottom drag [m/s]
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zdf_drg  ***
      !!
      !! ** Purpose :   update the top/bottom drag coefficient (non-linear case only)
      !!
      !! ** Method  :   In non linear friction case, the drag coeficient is
      !!              a function of the velocity:
      !!                          Cd = cd0 * |U+Ut|   
      !!              where U is the top or bottom velocity and
      !!                    Ut a tidal velocity (Ut^2 = Tidal kinetic energy 
      !!                       assumed here here to be constant) 
      !!              Depending on the input variable, the top- or bottom drag is compted
      !!
      !! ** Action  :   p_Cd   drag coefficient at t-point
      !!----------------------------------------------------------------------
      INTEGER                 , INTENT(in   ) ::   kt       ! ocean time-step index
      !                       !               !!         !==  top or bottom variables  ==!
      INTEGER , DIMENSION(:,:), INTENT(in   ) ::   k_mk     ! wet level (1st or last)
      REAL(wp)                , INTENT(in   ) ::   pCdmin   ! min drag value
      REAL(wp)                , INTENT(in   ) ::   pCdmax   ! max drag value
      REAL(wp)                , INTENT(in   ) ::   pz0      ! roughness
      REAL(wp)                , INTENT(in   ) ::   pke0     ! background tidal KE
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   pCd0     ! masked precomputed part of Cd0
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   pCdU     ! = - Cd*|U|   (t-points) [m/s]
      !!
      INTEGER ::   ji, jj   ! dummy loop indices
      INTEGER ::   imk      ! local integers
      REAL(wp)::   zzz, zut, zvt, zcd   ! local scalars
      !!----------------------------------------------------------------------
      !
      IF( l_log_not_linssh ) THEN     !==  "log layer"  ==!   compute Cd and -Cd*|U|
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               imk = k_mk(ji,jj)          ! ocean bottom level at t-points
               zut = un(ji,jj,imk) + un(ji-1,jj,imk)     ! 2 x velocity at t-point
               zvt = vn(ji,jj,imk) + vn(ji,jj-1,imk)
               zzz = 0.5_wp * e3t_n(ji,jj,imk)           ! altitude below/above (top/bottom) the boundary
               !
!!JC: possible WAD implementation should modify line below if layers vanish
               zcd = (  vkarmn / LOG( zzz / pz0 )  )**2
               zcd = pCd0(ji,jj) * MIN(  MAX( pCdmin , zcd ) , pCdmax  )   ! here pCd0 = mask*boost
               pCdU(ji,jj) = - zcd * SQRT(  0.25 * ( zut*zut + zvt*zvt ) + pke0  )
            END DO
         END DO
      ELSE                                            !==  standard Cd  ==!
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               imk = k_mk(ji,jj)    ! ocean bottom level at t-points
               zut = un(ji,jj,imk) + un(ji-1,jj,imk)     ! 2 x velocity at t-point
               zvt = vn(ji,jj,imk) + vn(ji,jj-1,imk)
               !                                                           ! here pCd0 = mask*boost * drag
               pCdU(ji,jj) = - pCd0(ji,jj) * SQRT(  0.25 * ( zut*zut + zvt*zvt ) + pke0  )
            END DO
         END DO
      ENDIF
      !
      IF(ln_ctl)   CALL prt_ctl( tab2d_1=pCdU, clinfo1=' Cd*U ')
      !
   END SUBROUTINE zdf_drg


   SUBROUTINE zdf_drg_exp( kt, pub, pvb, pua, pva )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_drg_exp  ***
      !!
      !! ** Purpose :   compute and add the explicit top and bottom frictions.
      !!
      !! ** Method  :   in explicit case, 
      !!
      !!              NB: in implicit case the calculation is performed in dynzdf.F90
      !!
      !! ** Action  :   (pua,pva)   momentum trend increased by top & bottom friction trend
      !!---------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt         ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pub, pvb   ! the two components of the before velocity
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pua, pva   ! the two components of the velocity tendency
      !! 
      INTEGER  ::   ji, jj       ! dummy loop indexes
      INTEGER  ::   ikbu, ikbv   ! local integers
      REAL(wp) ::   zm1_2dt      ! local scalar
      REAL(wp) ::   zCdu, zCdv   !   -      -
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdu, ztrdv
      !!---------------------------------------------------------------------
      !
!!gm bug : time step is only rdt (not 2 rdt if euler start !)
      zm1_2dt = - 1._wp / ( 2._wp * rdt )

      IF( l_trddyn ) THEN      ! trends: store the input trends
         ALLOCATE( ztrdu(jpi,jpj,jpk) , ztrdv(jpi,jpj,jpk) )
         ztrdu(:,:,:) = pua(:,:,:)
         ztrdv(:,:,:) = pva(:,:,:)
      ENDIF

      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            ikbu = mbku(ji,jj)          ! deepest wet ocean u- & v-levels
            ikbv = mbkv(ji,jj)
            !
            ! Apply stability criteria on absolute value  : abs(bfr/e3) < 1/(2dt) => bfr/e3 > -1/(2dt)
            zCdu = 0.5*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) ) / e3u_n(ji,jj,ikbu)
            zCdv = 0.5*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) ) / e3v_n(ji,jj,ikbv)
            !
            pua(ji,jj,ikbu) = pua(ji,jj,ikbu) + MAX(  zCdu , zm1_2dt  ) * pub(ji,jj,ikbu)
            pva(ji,jj,ikbv) = pva(ji,jj,ikbv) + MAX(  zCdv , zm1_2dt  ) * pvb(ji,jj,ikbv)
         END DO
      END DO
      !
      IF( ln_isfcav ) THEN        ! ocean cavities
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ikbu = miku(ji,jj)          ! first wet ocean u- & v-levels
               ikbv = mikv(ji,jj)
               !
               ! Apply stability criteria on absolute value  : abs(bfr/e3) < 1/(2dt) => bfr/e3 > -1/(2dt)
               zCdu = 0.5*( rCdU_top(ji+1,jj)+rCdU_top(ji,jj) ) / e3u_n(ji,jj,ikbu)    ! NB: Cdtop masked
               zCdv = 0.5*( rCdU_top(ji,jj+1)+rCdU_top(ji,jj) ) / e3v_n(ji,jj,ikbv)
               !
               pua(ji,jj,ikbu) = pua(ji,jj,ikbu) + MAX(  zCdu , zm1_2dt  ) * pub(ji,jj,ikbu)
               pva(ji,jj,ikbv) = pva(ji,jj,ikbv) + MAX(  zCdv , zm1_2dt  ) * pvb(ji,jj,ikbv)
           END DO
         END DO
      ENDIF
      !
      IF( l_trddyn ) THEN      ! trends: send trends to trddyn for further diagnostics
         ztrdu(:,:,:) = pua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = pva(:,:,:) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu(:,:,:), ztrdv(:,:,:), jpdyn_bfr, kt )
         DEALLOCATE( ztrdu, ztrdv )
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=pua, clinfo1=' bfr  - Ua: ', mask1=umask,               &
         &                       tab3d_2=pva, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
   END SUBROUTINE zdf_drg_exp


   SUBROUTINE zdf_drg_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_brg_init  ***
      !!
      !! ** Purpose :   Initialization of the bottom friction
      !!
      !! ** Method  :   Read the namdrg namelist and check their consistency
      !!                called at the first timestep (nit000)
      !!----------------------------------------------------------------------
      INTEGER   ::   ji, jj      ! dummy loop indexes
      INTEGER   ::   ios, ioptio   ! local integers
      !!
      NAMELIST/namdrg/ ln_drg_OFF, ln_lin, ln_non_lin, ln_loglayer, ln_drgimp, ln_drgice_imp
      !!----------------------------------------------------------------------
      !
      !                     !==  drag nature  ==!
      !
      REWIND( numnam_ref )                   ! Namelist namdrg in reference namelist
      READ  ( numnam_ref, namdrg, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam( ios , 'namdrg in reference namelist' )
      REWIND( numnam_cfg )                   ! Namelist namdrg in configuration namelist
      READ  ( numnam_cfg, namdrg, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam( ios , 'namdrg in configuration namelist' )
      IF(lwm) WRITE ( numond, namdrg )
      !
      IF( ln_drgice_imp .AND. nn_ice /= 2 )   ln_drgice_imp = .FALSE.
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_drg_init : top and/or bottom drag setting'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdrg : top/bottom friction choices'
         WRITE(numout,*) '      free-slip       : Cd = 0                  ln_drg_OFF  = ', ln_drg_OFF 
         WRITE(numout,*) '      linear  drag    : Cd = Cd0                ln_lin      = ', ln_lin
         WRITE(numout,*) '      non-linear  drag: Cd = Cd0_nl |U|         ln_non_lin  = ', ln_non_lin
         WRITE(numout,*) '      logarithmic drag: Cd = vkarmn/log(z/z0)   ln_loglayer = ', ln_loglayer
         WRITE(numout,*) '      implicit friction                         ln_drgimp   = ', ln_drgimp
         WRITE(numout,*) '      implicit ice-ocean drag                   ln_drgice_imp  =', ln_drgice_imp
      ENDIF
      !
      ioptio = 0                       ! set ndrg and control check
      IF( ln_drg_OFF  ) THEN   ;   ndrg = np_OFF        ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ln_lin      ) THEN   ;   ndrg = np_lin        ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ln_non_lin  ) THEN   ;   ndrg = np_non_lin    ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ln_loglayer ) THEN   ;   ndrg = np_loglayer   ;   ioptio = ioptio + 1   ;   ENDIF
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'zdf_drg_init: Choose ONE type of drag coef in namdrg' )
      !
      IF ( ln_drgice_imp.AND.(.NOT.ln_drgimp) ) & 
         &                CALL ctl_stop( 'zdf_drg_init: ln_drgice_imp=T requires ln_drgimp=T' )
      !
      !                     !==  BOTTOM drag setting  ==!   (applied at seafloor)
      !
      ALLOCATE( rCd0_bot(jpi,jpj), rCdU_bot(jpi,jpj) )
      CALL drg_init( 'BOTTOM'   , mbkt       ,                                         &   ! <== in
         &           r_Cdmin_bot, r_Cdmax_bot, r_z0_bot, r_ke0_bot, rCd0_bot, rCdU_bot )   ! ==> out
      !
      !                     !==  TOP drag setting  ==!   (applied at the top of ocean cavities)
      !
      IF( ln_isfcav.OR.ln_drgice_imp ) THEN              ! Ocean cavities: top friction setting
         ALLOCATE( rCdU_top(jpi,jpj) )
      ENDIF
      !
      IF( ln_isfcav ) THEN
         ALLOCATE( rCd0_top(jpi,jpj))
         CALL drg_init( 'TOP   '   , mikt       ,                                         &   ! <== in
            &           r_Cdmin_top, r_Cdmax_top, r_z0_top, r_ke0_top, rCd0_top, rCdU_top )   ! ==> out
      ENDIF
      !
   END SUBROUTINE zdf_drg_init


   SUBROUTINE drg_init( cd_topbot, k_mk,  &
      &                 pCdmin, pCdmax, pz0, pke0, pCd0, pCdU ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE drg_init  ***
      !!
      !! ** Purpose :   Initialization of the top/bottom friction CdO and Cd
      !!              from namelist parameters
      !!----------------------------------------------------------------------
      CHARACTER(len=6)        , INTENT(in   ) ::   cd_topbot       ! top/ bot indicator
      INTEGER , DIMENSION(:,:), INTENT(in   ) ::   k_mk            ! 1st/last  wet level 
      REAL(wp)                , INTENT(  out) ::   pCdmin, pCdmax  ! min and max drag coef. [-]
      REAL(wp)                , INTENT(  out) ::   pz0             ! roughness              [m]
      REAL(wp)                , INTENT(  out) ::   pke0            ! background KE          [m2/s2]
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   pCd0            ! masked precomputed part of the non-linear drag coefficient
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   pCdU            ! minus linear drag*|U| at t-points  [m/s]
      !!
      CHARACTER(len=40) ::   cl_namdrg, cl_file, cl_varname, cl_namref, cl_namcfg  ! local names 
      INTEGER ::   ji, jj              ! dummy loop indexes
      LOGICAL ::   ll_top, ll_bot      ! local logical
      INTEGER ::   ios, inum, imk      ! local integers
      REAL(wp)::   zmsk, zzz, zcd      ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   zmsk_boost   ! 2D workspace
      !!
      NAMELIST/namdrg_top/ rn_Cd0, rn_Uc0, rn_Cdmax, rn_ke0, rn_z0, ln_boost, rn_boost
      NAMELIST/namdrg_bot/ rn_Cd0, rn_Uc0, rn_Cdmax, rn_ke0, rn_z0, ln_boost, rn_boost
      !!----------------------------------------------------------------------
      !
      !                          !==  set TOP / BOTTOM specificities  ==!
      ll_top = .FALSE.
      ll_bot = .FALSE.
      !
      SELECT CASE (cd_topbot)
      CASE( 'TOP   ' )
         ll_top = .TRUE.
         cl_namdrg  = 'namdrg_top'
         cl_namref  = 'namdrg_top in reference     namelist'
         cl_namcfg  = 'namdrg_top in configuration namelist'
         cl_file    = 'tfr_coef.nc'
         cl_varname = 'tfr_coef'
      CASE( 'BOTTOM' )
         ll_bot = .TRUE.
         cl_namdrg  = 'namdrg_bot'
         cl_namref  = 'namdrg_bot  in reference     namelist'
         cl_namcfg  = 'namdrg_bot  in configuration namelist'
         cl_file    = 'bfr_coef.nc'
         cl_varname = 'bfr_coef'
      CASE DEFAULT
         CALL ctl_stop( 'drg_init: bad value for cd_topbot ' )
      END SELECT
      !
      !                          !==  read namlist  ==!
      !
      REWIND( numnam_ref )                   ! Namelist cl_namdrg in reference namelist
      IF(ll_top)   READ  ( numnam_ref, namdrg_top, IOSTAT = ios, ERR = 901)
      IF(ll_bot)   READ  ( numnam_ref, namdrg_bot, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam( ios , TRIM(cl_namref) )
      REWIND( numnam_cfg )                   ! Namelist cd_namdrg in configuration namelist
      IF(ll_top)   READ  ( numnam_cfg, namdrg_top, IOSTAT = ios, ERR = 902 )
      IF(ll_bot)   READ  ( numnam_cfg, namdrg_bot, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam( ios , TRIM(cl_namcfg) )
      IF(lwm .AND. ll_top)   WRITE ( numond, namdrg_top )
      IF(lwm .AND. ll_bot)   WRITE ( numond, namdrg_bot )
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist ',TRIM(cl_namdrg),' : set ',TRIM(cd_topbot),' friction parameters'
         WRITE(numout,*) '      drag coefficient                        rn_Cd0   = ', rn_Cd0
         WRITE(numout,*) '      characteristic velocity (linear case)   rn_Uc0   = ', rn_Uc0, ' m/s'
         WRITE(numout,*) '      non-linear drag maximum                 rn_Cdmax = ', rn_Cdmax
         WRITE(numout,*) '      background kinetic energy  (n-l case)   rn_ke0   = ', rn_ke0
         WRITE(numout,*) '      bottom roughness           (n-l case)   rn_z0    = ', rn_z0
         WRITE(numout,*) '      set a regional boost of Cd0             ln_boost = ', ln_boost
         WRITE(numout,*) '         associated boost factor              rn_boost = ', rn_boost
      ENDIF
      !
      !                          !==  return some namelist parametres  ==!   (used in non_lin and loglayer cases)
      pCdmin = rn_Cd0
      pCdmax = rn_Cdmax
      pz0    = rn_z0
      pke0   = rn_ke0
      !
      !                          !==  mask * boost factor  ==!
      !
      IF( ln_boost ) THEN           !* regional boost:   boost factor = 1 + regional boost
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   use a regional boost read in ', TRIM(cl_file), ' file'
         IF(lwp) WRITE(numout,*) '           using enhancement factor of ', rn_boost
         ! cl_varname is a coefficient in [0,1] giving where to apply the regional boost
         CALL iom_open ( TRIM(cl_file), inum )
         CALL iom_get  ( inum, jpdom_data, TRIM(cl_varname), zmsk_boost, 1 )
         CALL iom_close( inum)
         zmsk_boost(:,:) = 1._wp + rn_boost * zmsk_boost(:,:)
         !
      ELSE                          !* no boost:   boost factor = 1
         zmsk_boost(:,:) = 1._wp
      ENDIF
      !                             !* mask outside ocean cavities area (top) or land area (bot)
      IF(ll_top)   zmsk_boost(:,:) = zmsk_boost(:,:) * ssmask(:,:) * (1. - tmask(:,:,1) )  ! none zero in ocean cavities only
      IF(ll_bot)   zmsk_boost(:,:) = zmsk_boost(:,:) * ssmask(:,:)                         ! x seafloor mask
      !
      !
      SELECT CASE( ndrg )
      !
      CASE( np_OFF  )            !==  No top/bottom friction  ==!   (pCdU = 0)
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   ',TRIM(cd_topbot),' free-slip, friction set to zero'
         !
         l_zdfdrg = .FALSE.         ! no time variation of the drag: set it one for all
         !
         pCdU(:,:) = 0._wp
         pCd0(:,:) = 0._wp
         !
      CASE( np_lin )             !==  linear friction  ==!   (pCdU = Cd0 * Uc0)
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   linear ',TRIM(cd_topbot),' friction (constant coef = Cd0*Uc0 = ', rn_Cd0*rn_Uc0, ')'
         !
         l_zdfdrg = .FALSE.         ! no time variation of the Cd*|U| : set it one for all
         !                      
         pCd0(:,:) = rn_Cd0 * zmsk_boost(:,:)  !* constant in time drag coefficient (= mask (and boost) Cd0)
         pCdU(:,:) = - pCd0(:,:) * rn_Uc0      !  using a constant velocity
         !
      CASE( np_non_lin )         !== non-linear friction  ==!   (pCd0 = Cd0 )
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   quadratic ',TRIM(cd_topbot),' friction (propotional to module of the velocity)'
         IF(lwp) WRITE(numout,*) '   with    a drag coefficient Cd0 = ', rn_Cd0, ', and'
         IF(lwp) WRITE(numout,*) '           a background velocity module of (rn_ke0)^1/2 = ', SQRT(rn_ke0), 'm/s)'
         !
         l_zdfdrg = .TRUE.          !* Cd*|U| updated at each time-step (it depends on ocean velocity)
         !
         pCd0(:,:) = rn_Cd0 * zmsk_boost(:,:)  !* constant in time proportionality coefficient (= mask (and boost) Cd0)
         pCdU(:,:) = 0._wp                     !  
         !
      CASE( np_loglayer )       !== logarithmic layer formulation of friction  ==!   (CdU = (vkarman log(z/z0))^2 |U| )
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   quadratic ',TRIM(cd_topbot),' drag (propotional to module of the velocity)'
         IF(lwp) WRITE(numout,*) '   with   a logarithmic Cd0 formulation Cd0 = ( vkarman log(z/z0) )^2 ,'
         IF(lwp) WRITE(numout,*) '          a background velocity module of (rn_ke0)^1/2 = ', SQRT(pke0), 'm/s), '
         IF(lwp) WRITE(numout,*) '          a logarithmic formulation: a roughness of ', pz0, ' meters,   and '
         IF(lwp) WRITE(numout,*) '          a proportionality factor bounded by min/max values of ', pCdmin, pCdmax
         !
         l_zdfdrg = .TRUE.          !* Cd*|U| updated at each time-step (it depends on ocean velocity)
         !
         IF( ln_linssh ) THEN       !* pCd0 = (v log(z/z0))^2   as velocity points have a fixed z position
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   N.B.   linear free surface case, Cd0 computed one for all'
            !
            l_log_not_linssh = .FALSE.    !- don't update Cd at each time step
            !
            DO jj = 1, jpj                   ! pCd0 = mask (and boosted) logarithmic drag coef. 
               DO ji = 1, jpi
                  zzz =  0.5_wp * e3t_0(ji,jj,k_mk(ji,jj))
                  zcd = (  vkarmn / LOG( zzz / rn_z0 )  )**2
                  pCd0(ji,jj) = zmsk_boost(ji,jj) * MIN(  MAX( rn_Cd0 , zcd ) , rn_Cdmax  )  ! rn_Cd0 < Cd0 < rn_Cdmax
               END DO
            END DO
         ELSE                       !* Cd updated at each time-step ==> pCd0 = mask * boost
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   N.B.   non-linear free surface case, Cd0 updated at each time-step '
            !
            l_log_not_linssh = .TRUE.     ! compute the drag coef. at each time-step 
            !
            pCd0(:,:) = zmsk_boost(:,:)
         ENDIF
         pCdU(:,:) = 0._wp          ! initialisation to zero (will be updated at each time step)
         !
      CASE DEFAULT
         CALL ctl_stop( 'drg_init: bad flag value for ndrg ' )
      END SELECT
      !
   END SUBROUTINE drg_init

   !!======================================================================
END MODULE zdfdrg
