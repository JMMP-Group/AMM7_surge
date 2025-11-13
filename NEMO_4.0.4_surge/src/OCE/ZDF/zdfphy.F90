MODULE zdfphy
   !!======================================================================
   !!                      ***  MODULE  zdfphy  ***
   !! Vertical ocean physics :   manager of all vertical physics packages
   !!======================================================================
   !! History :  4.0  !  2017-04  (G. Madec)  original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_phy_init  : initialization of all vertical physics packages
   !!   zdf_phy       : upadate at each time-step the vertical mixing coeff. 
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE zdf_oce        ! vertical physics: shared variables         
   USE zdfdrg         ! vertical physics: top/bottom drag coef.
   USE zdfsh2         ! vertical physics: shear production term of TKE
   USE zdfric         ! vertical physics: RIChardson dependent vertical mixing   
   USE zdftke         ! vertical physics: TKE vertical mixing
   USE zdfgls         ! vertical physics: GLS vertical mixing
   USE zdfosm         ! vertical physics: OSMOSIS vertical mixing
   USE zdfddm         ! vertical physics: double diffusion mixing      
   USE zdfevd         ! vertical physics: convection via enhanced vertical diffusion  
   USE zdfiwm         ! vertical physics: internal wave-induced mixing  
   USE zdfswm         ! vertical physics: surface  wave-induced mixing
   USE zdfmxl         ! vertical physics: mixed layer
   USE tranpc         ! convection: non penetrative adjustment
   USE trc_oce        ! variables shared between passive tracer & ocean           
   USE sbc_oce        ! surface module (only for nn_isf in the option compatibility test)
   USE sbcrnf         ! surface boundary condition: runoff variables
   USE sbc_ice        ! sea ice drag
#if defined key_agrif
   USE agrif_oce_interp   ! interpavm
#endif
   !
   USE in_out_manager ! I/O manager
   USE iom            ! IOM library
   USE lbclnk         ! lateral boundary conditions
   USE lib_mpp        ! distribued memory computing
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_phy_init  ! called by nemogcm.F90
   PUBLIC   zdf_phy       ! called by step.F90

   INTEGER ::   nzdf_phy   ! type of vertical closure used 
   !                       ! associated indicators
   INTEGER, PARAMETER ::   np_CST = 1   ! Constant Kz
   INTEGER, PARAMETER ::   np_RIC = 2   ! Richardson number dependent Kz
   INTEGER, PARAMETER ::   np_TKE = 3   ! Turbulente Kinetic Eenergy closure scheme for Kz
   INTEGER, PARAMETER ::   np_GLS = 4   ! Generic Length Scale closure scheme for Kz
   INTEGER, PARAMETER ::   np_OSM = 5   ! OSMOSIS-OBL closure scheme for Kz

   LOGICAL ::   l_zdfsh2   ! shear production term flag (=F for CST, =T otherwise (i.e. TKE, GLS, RIC))

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE zdf_phy_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_phy_init  ***
      !! 
      !! ** Purpose :   initializations of the vertical ocean physics
      !!
      !! ** Method  :   Read namelist namzdf, control logicals 
      !!                set horizontal shape and vertical profile of background mixing coef.
      !!----------------------------------------------------------------------
      INTEGER ::   jk            ! dummy loop indices
      INTEGER ::   ioptio, ios   ! local integers
      !!
      NAMELIST/namzdf/ ln_zdfcst, ln_zdfric, ln_zdftke, ln_zdfgls,   &     ! type of closure scheme
         &             ln_zdfosm,                                    &     ! type of closure scheme
         &             ln_zdfevd, nn_evdm, rn_evd ,                  &     ! convection : evd
         &             ln_zdfnpc, nn_npc , nn_npcp,                  &     ! convection : npc
         &             ln_zdfddm, rn_avts, rn_hsbfr,                 &     ! double diffusion
         &             ln_zdfswm,                                    &     ! surface  wave-induced mixing
         &             ln_zdfiwm,                                    &     ! internal  -      -      -
         &             ln_zad_Aimp,                                  &     ! apdative-implicit vertical advection
         &             rn_avm0, rn_avt0, nn_avb, nn_havtb                  ! coefficients
      !!----------------------------------------------------------------------
      !
      IF( .NOT. ln_2d ) THEN

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_phy_init: ocean vertical physics'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      !                           !==  Namelist  ==!
      REWIND( numnam_ref )              ! Namelist namzdf in reference namelist : Vertical mixing parameters
      READ  ( numnam_ref, namzdf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namzdf in reference namelist' )
      !
      REWIND( numnam_cfg )              ! Namelist namzdf in reference namelist : Vertical mixing parameters
      READ  ( numnam_cfg, namzdf, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namzdf in configuration namelist' )
      IF(lwm)   WRITE ( numond, namzdf )
      !
      IF(lwp) THEN                      ! Parameter print
         WRITE(numout,*) '   Namelist namzdf : set vertical mixing mixing parameters'
         WRITE(numout,*) '      adaptive-implicit vertical advection'
         WRITE(numout,*) '         Courant number targeted application   ln_zad_Aimp = ', ln_zad_Aimp
         WRITE(numout,*) '      vertical closure scheme'
         WRITE(numout,*) '         constant vertical mixing coefficient    ln_zdfcst = ', ln_zdfcst
         WRITE(numout,*) '         Richardson number dependent closure     ln_zdfric = ', ln_zdfric
         WRITE(numout,*) '         Turbulent Kinetic Energy closure (TKE)  ln_zdftke = ', ln_zdftke
         WRITE(numout,*) '         Generic Length Scale closure (GLS)      ln_zdfgls = ', ln_zdfgls
         WRITE(numout,*) '         OSMOSIS-OBL closure (OSM)               ln_zdfosm = ', ln_zdfosm
         WRITE(numout,*) '      convection: '
         WRITE(numout,*) '         enhanced vertical diffusion             ln_zdfevd = ', ln_zdfevd
         WRITE(numout,*) '            applied on momentum (=1/0)             nn_evdm = ', nn_evdm
         WRITE(numout,*) '            vertical coefficient for evd           rn_evd  = ', rn_evd
         WRITE(numout,*) '         non-penetrative convection (npc)        ln_zdfnpc = ', ln_zdfnpc
         WRITE(numout,*) '            npc call  frequency                    nn_npc  = ', nn_npc
         WRITE(numout,*) '            npc print frequency                    nn_npcp = ', nn_npcp
         WRITE(numout,*) '      double diffusive mixing                    ln_zdfddm = ', ln_zdfddm
         WRITE(numout,*) '         maximum avs for dd mixing                 rn_avts = ', rn_avts
         WRITE(numout,*) '         heat/salt buoyancy flux ratio             rn_hsbfr= ', rn_hsbfr
         WRITE(numout,*) '      gravity wave-induced mixing'
         WRITE(numout,*) '         surface  wave (Qiao et al 2010)         ln_zdfswm = ', ln_zdfswm                                          ! surface wave induced mixing
         WRITE(numout,*) '         internal wave (de Lavergne et al 2017)  ln_zdfiwm = ', ln_zdfiwm
         WRITE(numout,*) '      coefficients : '
         WRITE(numout,*) '         vertical eddy viscosity                 rn_avm0   = ', rn_avm0
         WRITE(numout,*) '         vertical eddy diffusivity               rn_avt0   = ', rn_avt0
         WRITE(numout,*) '         constant background or profile          nn_avb    = ', nn_avb
         WRITE(numout,*) '         horizontal variation for avtb           nn_havtb  = ', nn_havtb
      ENDIF

      IF( ln_zad_Aimp ) THEN
         IF( zdf_phy_alloc() /= 0 )   &
            &       CALL ctl_stop( 'STOP', 'zdf_phy_init : unable to allocate adaptive-implicit z-advection arrays' )
         Cu_adv(:,:,:) = 0._wp
         wi    (:,:,:) = 0._wp
      ENDIF
      !                          !==  Background eddy viscosity and diffusivity  ==!
      IF( nn_avb == 0 ) THEN             ! Define avmb, avtb from namelist parameter
         avmb(:) = rn_avm0
         avtb(:) = rn_avt0                     
      ELSE                               ! Background profile of avt (fit a theoretical/observational profile (Krauss 1990)
         avmb(:) = rn_avm0
         avtb(:) = rn_avt0 + ( 3.e-4_wp - 2._wp * rn_avt0 ) * 1.e-4_wp * gdepw_1d(:)   ! m2/s
         IF(ln_sco .AND. lwp)   CALL ctl_warn( 'avtb profile not valid in sco' )
      ENDIF
      !                                  ! 2D shape of the avtb
      avtb_2d(:,:) = 1._wp                   ! uniform 
      !
      IF( nn_havtb == 1 ) THEN               ! decrease avtb by a factor of ten in the equatorial band
           !                                 !   -15S -5S : linear decrease from avt0 to avt0/10.
           !                                 !   -5S  +5N : cst value avt0/10.
           !                                 !    5N  15N : linear increase from avt0/10, to avt0
           WHERE(-15. <= gphit .AND. gphit < -5 )   avtb_2d = (1.  - 0.09 * (gphit + 15.))
           WHERE( -5. <= gphit .AND. gphit <  5 )   avtb_2d =  0.1
           WHERE(  5. <= gphit .AND. gphit < 15 )   avtb_2d = (0.1 + 0.09 * (gphit -  5.))
      ENDIF
      !
      DO jk = 1, jpk                      ! set turbulent closure Kz to the background value (avt_k, avm_k)
         avt_k(:,:,jk) = avtb_2d(:,:) * avtb(jk) * wmask (:,:,jk)
         avm_k(:,:,jk) =                avmb(jk) * wmask (:,:,jk)
      END DO
!!gm  to be tested only the 1st & last levels
!      avt  (:,:, 1 ) = 0._wp   ;   avs(:,:, 1 ) = 0._wp   ;   avm  (:,:, 1 ) = 0._wp
!      avt  (:,:,jpk) = 0._wp   ;   avs(:,:,jpk) = 0._wp   ;   avm  (:,:,jpk) = 0._wp
!!gm
      avt  (:,:,:) = 0._wp   ;   avs(:,:,:) = 0._wp   ;   avm  (:,:,:) = 0._wp

      !                          !==  Convection  ==!
      !
      IF( ln_zdfnpc .AND. ln_zdfevd )   CALL ctl_stop( 'zdf_phy_init: chose between ln_zdfnpc and ln_zdfevd' )
      IF( ln_zdfosm .AND. ln_zdfevd )   CALL ctl_stop( 'zdf_phy_init: chose between ln_zdfosm and ln_zdfevd' )
      IF( lk_top    .AND. ln_zdfnpc )   CALL ctl_stop( 'zdf_phy_init: npc scheme is not working with key_top' )
      IF( lk_top    .AND. ln_zdfosm )   CALL ctl_stop( 'zdf_phy_init: osmosis scheme is not working with key_top' )
      IF(lwp) THEN
         WRITE(numout,*)
         IF    ( ln_zdfnpc ) THEN  ;   WRITE(numout,*) '   ==>>>   convection: use non penetrative convective scheme'
         ELSEIF( ln_zdfevd ) THEN  ;   WRITE(numout,*) '   ==>>>   convection: use enhanced vertical diffusion scheme'
         ELSE                      ;   WRITE(numout,*) '   ==>>>   convection: no specific scheme used'
         ENDIF
      ENDIF

      IF(lwp) THEN               !==  Double Diffusion Mixing parameterization  ==!   (ddm)
         WRITE(numout,*)
         IF( ln_zdfddm ) THEN   ;   WRITE(numout,*) '   ==>>>   use double diffusive mixing: avs /= avt'
         ELSE                   ;   WRITE(numout,*) '   ==>>>   No  double diffusive mixing: avs = avt'
         ENDIF
      ENDIF

      !                          !==  type of vertical turbulent closure  ==!    (set nzdf_phy)
      ioptio = 0 
      IF( ln_zdfcst ) THEN   ;   ioptio = ioptio + 1   ;    nzdf_phy = np_CST   ;   ENDIF
      IF( ln_zdfric ) THEN   ;   ioptio = ioptio + 1   ;    nzdf_phy = np_RIC   ;   CALL zdf_ric_init   ;   ENDIF
      IF( ln_zdftke ) THEN   ;   ioptio = ioptio + 1   ;    nzdf_phy = np_TKE   ;   CALL zdf_tke_init   ;   ENDIF
      IF( ln_zdfgls ) THEN   ;   ioptio = ioptio + 1   ;    nzdf_phy = np_GLS   ;   CALL zdf_gls_init   ;   ENDIF
      IF( ln_zdfosm ) THEN   ;   ioptio = ioptio + 1   ;    nzdf_phy = np_OSM   ;   CALL zdf_osm_init   ;   ENDIF
      !
      IF( ioptio /= 1 )    CALL ctl_stop( 'zdf_phy_init: one and only one vertical diffusion option has to be defined ' )
      IF( ln_isfcav ) THEN
      IF( ln_zdfric .OR. ln_zdfgls )    CALL ctl_stop( 'zdf_phy_init: zdfric and zdfgls never tested with ice shelves cavities ' )
      ENDIF
      !                                ! shear production term flag
      IF( ln_zdfcst ) THEN   ;   l_zdfsh2 = .FALSE.
      ELSE                   ;   l_zdfsh2 = .TRUE.
      ENDIF

      !                          !== gravity wave-driven mixing  ==!
      IF( ln_zdfiwm )   CALL zdf_iwm_init       ! internal wave-driven mixing
      IF( ln_zdfswm )   CALL zdf_swm_init       ! surface  wave-driven mixing

      ENDIF  ! .NOT. ln_2d

      !                          !== top/bottom friction  ==!
      CALL zdf_drg_init
      !
      !                          !== time-stepping  ==!
      ! Check/update of time stepping done in dynzdf_init/trazdf_init
      !!gm move it here ?
      !
   END SUBROUTINE zdf_phy_init


   SUBROUTINE zdf_phy( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE zdf_phy  ***
      !!
      !! ** Purpose :  Update ocean physics at each time-step
      !!
      !! ** Method  : 
      !!
      !! ** Action  :   avm, avt vertical eddy viscosity and diffusivity at w-points
      !!                nmld ??? mixed layer depth in level and meters   <<<<====verifier !
      !!                bottom stress.....                               <<<<====verifier !
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indice
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zsh2   ! shear production
      !! ---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('zdf_phy')
      !
      IF( l_zdfdrg ) THEN     !==  update top/bottom drag  ==!   (non-linear cases)
         !
         !                       !* bottom drag
         CALL zdf_drg( kt, mbkt    , r_Cdmin_bot, r_Cdmax_bot,   &   ! <<== in 
            &              r_z0_bot,   r_ke0_bot,    rCd0_bot,   &
            &                                        rCdU_bot  )     ! ==>> out : bottom drag [m/s]
         IF( ln_isfcav ) THEN    !* top drag   (ocean cavities)
            CALL zdf_drg( kt, mikt    , r_Cdmin_top, r_Cdmax_top,   &   ! <<== in 
               &              r_z0_top,   r_ke0_top,    rCd0_top,   &
               &                                        rCdU_top  )     ! ==>> out : bottom drag [m/s]
         ENDIF
      ENDIF
      
      IF( .NOT. ln_2d ) THEN ! 2D case only uses bottom friction from this routine     
      !
#if defined key_si3
      IF ( ln_drgice_imp) THEN
         IF ( ln_isfcav ) THEN
            rCdU_top(:,:) = rCdU_top(:,:) + ssmask(:,:) * tmask(:,:,1) * rCdU_ice(:,:)
         ELSE
            rCdU_top(:,:) = rCdU_ice(:,:)
         ENDIF
      ENDIF
#endif
      ! 
      !                       !==  Kz from chosen turbulent closure  ==!   (avm_k, avt_k)
      !
      IF( l_zdfsh2 )   &         !* shear production at w-points (energy conserving form)
         CALL zdf_sh2( ub, vb, un, vn, avm_k,   &     ! <<== in
            &                           zsh2    )     ! ==>> out : shear production
      !
      SELECT CASE ( nzdf_phy )                  !* Vertical eddy viscosity and diffusivity coefficients at w-points
      CASE( np_RIC )   ;   CALL zdf_ric( kt, gdept_n, zsh2, avm_k, avt_k )    ! Richardson number dependent Kz
      CASE( np_TKE )   ;   CALL zdf_tke( kt         , zsh2, avm_k, avt_k )    ! TKE closure scheme for Kz
      CASE( np_GLS )   ;   CALL zdf_gls( kt         , zsh2, avm_k, avt_k )    ! GLS closure scheme for Kz
      CASE( np_OSM )   ;   CALL zdf_osm( kt               , avm_k, avt_k )    ! OSMOSIS closure scheme for Kz
!     CASE( np_CST )                                  ! Constant Kz (reset avt, avm to the background value)
!         ! avt_k and avm_k set one for all at initialisation phase
!!gm         avt(2:jpim1,2:jpjm1,1:jpkm1) = rn_avt0 * wmask(2:jpim1,2:jpjm1,1:jpkm1)
!!gm         avm(2:jpim1,2:jpjm1,1:jpkm1) = rn_avm0 * wmask(2:jpim1,2:jpjm1,1:jpkm1)
      END SELECT
      !  
      !                          !==  ocean Kz  ==!   (avt, avs, avm)
      !
      !                                         !* start from turbulent closure values
      avt(:,:,2:jpkm1) = avt_k(:,:,2:jpkm1)
      avm(:,:,2:jpkm1) = avm_k(:,:,2:jpkm1)
      !
      IF( ln_rnf_mouth ) THEN                   !* increase diffusivity at rivers mouths
         DO jk = 2, nkrnf
            avt(:,:,jk) = avt(:,:,jk) + 2._wp * rn_avt_rnf * rnfmsk(:,:) * wmask(:,:,jk)
         END DO
      ENDIF
      !
      IF( ln_zdfevd )   CALL zdf_evd( kt, avm, avt )  !* convection: enhanced vertical eddy diffusivity
      !
      !                                         !* double diffusive mixing
      IF( ln_zdfddm ) THEN                            ! update avt and compute avs
                        CALL zdf_ddm( kt, avm, avt, avs )
      ELSE                                            ! same mixing on all tracers
         avs(2:jpim1,2:jpjm1,1:jpkm1) = avt(2:jpim1,2:jpjm1,1:jpkm1)
      ENDIF
      !
      !                                         !* wave-induced mixing 
      IF( ln_zdfswm )   CALL zdf_swm( kt, avm, avt, avs )   ! surface  wave (Qiao et al. 2004) 
      IF( ln_zdfiwm )   CALL zdf_iwm( kt, avm, avt, avs )   ! internal wave (de Lavergne et al 2017)

#if defined key_agrif 
      ! interpolation parent grid => child grid for avm_k ( ex : at west border: update column 1 and 2)
      IF( l_zdfsh2 )   CALL Agrif_avm
#endif

      !                                         !* Lateral boundary conditions (sign unchanged)
      IF( l_zdfsh2 ) THEN
         CALL lbc_lnk_multi( 'zdfphy', avm_k, 'W', 1. , avt_k, 'W', 1.,   &
            &                avm  , 'W', 1. , avt  , 'W', 1. , avs , 'W', 1. )
      ELSE
         CALL lbc_lnk_multi( 'zdfphy', avm  , 'W', 1. , avt  , 'W', 1. , avs , 'W', 1. )
      ENDIF
      !
      IF( l_zdfdrg ) THEN     ! drag  have been updated (non-linear cases)
         IF( ln_isfcav ) THEN   ;  CALL lbc_lnk_multi( 'zdfphy', rCdU_top, 'T', 1. , rCdU_bot, 'T', 1. )   ! top & bot drag
         ELSE                   ;  CALL lbc_lnk      ( 'zdfphy', rCdU_bot, 'T', 1. )                       ! bottom drag only
         ENDIF
      ENDIF
      !
      CALL zdf_mxl( kt )                        !* mixed layer depth, and level
      !
      IF( lrst_oce ) THEN                       !* write TKE, GLS or RIC fields in the restart file
         IF( ln_zdftke )   CALL tke_rst( kt, 'WRITE' )
         IF( ln_zdfgls )   CALL gls_rst( kt, 'WRITE' )
         IF( ln_zdfric )   CALL ric_rst( kt, 'WRITE' ) 
         ! NB. OSMOSIS restart (osm_rst) will be called in step.F90 after wn has been updated
      ENDIF
      
      ENDIF  ! .NOT. ln_2d
      !
      IF( ln_timing )   CALL timing_stop('zdf_phy')
      !
   END SUBROUTINE zdf_phy
   INTEGER FUNCTION zdf_phy_alloc()
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION zdf_phy_alloc  ***
      !!----------------------------------------------------------------------
     ! Allocate wi array (declared in oce.F90) for use with the adaptive-implicit vertical velocity option
     ALLOCATE(     wi(jpi,jpj,jpk), Cu_adv(jpi,jpj,jpk),  STAT= zdf_phy_alloc )
     IF( zdf_phy_alloc /= 0 )   CALL ctl_warn('zdf_phy_alloc: failed to allocate ln_zad_Aimp=T required arrays')
     CALL mpp_sum ( 'zdfphy', zdf_phy_alloc )
   END FUNCTION zdf_phy_alloc

   !!======================================================================
END MODULE zdfphy
