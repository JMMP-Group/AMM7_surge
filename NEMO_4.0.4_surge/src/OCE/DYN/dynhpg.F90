MODULE dynhpg
   !!======================================================================
   !!                       ***  MODULE  dynhpg  ***
   !! Ocean dynamics:  hydrostatic pressure gradient trend
   !!======================================================================
   !! History :  OPA  !  1987-09  (P. Andrich, M.-A. Foujols)  hpg_zco: Original code
   !!            5.0  !  1991-11  (G. Madec)
   !!            7.0  !  1996-01  (G. Madec)  hpg_sco: Original code for s-coordinates
   !!            8.0  !  1997-05  (G. Madec)  split dynber into dynkeg and dynhpg
   !!            8.5  !  2002-07  (G. Madec)  F90: Free form and module
   !!            8.5  !  2002-08  (A. Bozec)  hpg_zps: Original code
   !!   NEMO     1.0  !  2005-10  (A. Beckmann, B.W. An)  various s-coordinate options
   !!                 !         Original code for hpg_ctl, hpg_hel hpg_wdj, hpg_djc, hpg_rot
   !!             -   !  2005-11  (G. Madec) style & small optimisation
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.4  !  2011-11  (H. Liu) hpg_prj: Original code for s-coordinates
   !!                 !           (A. Coward) suppression of hel, wdj and rot options
   !!            3.6  !  2014-11  (P. Mathiot) hpg_isf: original code for ice shelf cavity
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_hpg      : update the momentum trend with the now horizontal
   !!                  gradient of the hydrostatic pressure
   !!   dyn_hpg_init : initialisation and control of options
   !!       hpg_zco  : z-coordinate scheme
   !!       hpg_zps  : z-coordinate plus partial steps (interpolation)
   !!       hpg_sco  : s-coordinate (standard jacobian formulation)
   !!       hpg_isf  : s-coordinate (sco formulation) adapted to ice shelf
   !!       hpg_djc  : s-coordinate (Density Jacobian with Cubic polynomial)
   !!       hpg_prj  : s-coordinate (Pressure Jacobian with Cubic polynomial)
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE sbc_oce         ! surface variable (only for the flag with ice shelf)
   USE dom_oce         ! ocean space and time domain
   USE wet_dry         ! wetting and drying
   USE phycst          ! physical constants
   USE trd_oce         ! trends: ocean variables
   USE trddyn          ! trend manager: dynamics
   USE zpshde          ! partial step: hor. derivative     (zps_hde routine)
   !
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE lbclnk          ! lateral boundary condition 
   USE lib_mpp         ! MPP library
   USE eosbn2          ! compute density
   USE timing          ! Timing
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_hpg        ! routine called by step module
   PUBLIC   dyn_hpg_init   ! routine called by opa module

   !                                !!* Namelist namdyn_hpg : hydrostatic pressure gradient
   LOGICAL, PUBLIC ::   ln_hpg_zco   !: z-coordinate - full steps
   LOGICAL, PUBLIC ::   ln_hpg_zps   !: z-coordinate - partial steps (interpolation)
   LOGICAL, PUBLIC ::   ln_hpg_sco   !: s-coordinate (standard jacobian formulation)
   LOGICAL, PUBLIC ::   ln_hpg_djc   !: s-coordinate (Density Jacobian with Cubic polynomial)
   LOGICAL, PUBLIC ::   ln_hpg_prj   !: s-coordinate (Pressure Jacobian scheme)
   LOGICAL, PUBLIC ::   ln_hpg_isf   !: s-coordinate similar to sco modify for isf

   !                                !! Flag to control the type of hydrostatic pressure gradient
   INTEGER, PARAMETER ::   np_ERROR  =-10   ! error in specification of lateral diffusion
   INTEGER, PARAMETER ::   np_zco    =  0   ! z-coordinate - full steps
   INTEGER, PARAMETER ::   np_zps    =  1   ! z-coordinate - partial steps (interpolation)
   INTEGER, PARAMETER ::   np_sco    =  2   ! s-coordinate (standard jacobian formulation)
   INTEGER, PARAMETER ::   np_djc    =  3   ! s-coordinate (Density Jacobian with Cubic polynomial)
   INTEGER, PARAMETER ::   np_prj    =  4   ! s-coordinate (Pressure Jacobian scheme)
   INTEGER, PARAMETER ::   np_isf    =  5   ! s-coordinate similar to sco modify for isf
   !
   INTEGER, PUBLIC ::   nhpg         !: type of pressure gradient scheme used ! (deduced from ln_hpg_... flags) (PUBLIC for TAM)

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_hpg( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_hpg  ***
      !!
      !! ** Method  :   Call the hydrostatic pressure gradient routine
      !!              using the scheme defined in the namelist
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!             - send trends to trd_dyn for futher diagnostics (l_trddyn=T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrdu, ztrdv
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_hpg')
      !
      IF( l_trddyn ) THEN                    ! Temporary saving of ua and va trends (l_trddyn)
         ALLOCATE( ztrdu(jpi,jpj,jpk) , ztrdv(jpi,jpj,jpk) )
         ztrdu(:,:,:) = ua(:,:,:)
         ztrdv(:,:,:) = va(:,:,:)
      ENDIF
      !
      SELECT CASE ( nhpg )      ! Hydrostatic pressure gradient computation
      CASE ( np_zco )   ;   CALL hpg_zco    ( kt )      ! z-coordinate
      CASE ( np_zps )   ;   CALL hpg_zps    ( kt )      ! z-coordinate plus partial steps (interpolation)
      CASE ( np_sco )   ;   CALL hpg_sco    ( kt )      ! s-coordinate (standard jacobian formulation)
      CASE ( np_djc )   ;   CALL hpg_djc    ( kt )      ! s-coordinate (Density Jacobian with Cubic polynomial)
      CASE ( np_prj )   ;   CALL hpg_prj    ( kt )      ! s-coordinate (Pressure Jacobian scheme)
      CASE ( np_isf )   ;   CALL hpg_isf    ( kt )      ! s-coordinate similar to sco modify for ice shelf
      END SELECT
      !
      IF( l_trddyn ) THEN      ! save the hydrostatic pressure gradient trends for momentum trend diagnostics
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_hpg, kt )
         DEALLOCATE( ztrdu , ztrdv )
      ENDIF
      !
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' hpg  - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_hpg')
      !
   END SUBROUTINE dyn_hpg


   SUBROUTINE dyn_hpg_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dyn_hpg_init  ***
      !!
      !! ** Purpose :   initializations for the hydrostatic pressure gradient
      !!              computation and consistency control
      !!
      !! ** Action  :   Read the namelist namdyn_hpg and check the consistency
      !!      with the type of vertical coordinate used (zco, zps, sco)
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio = 0      ! temporary integer
      INTEGER ::   ios             ! Local integer output status for namelist read
      !!
      INTEGER  ::   ji, jj, jk, ikt    ! dummy loop indices      ISF
      REAL(wp) ::   znad
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  zts_top, zrhd   ! hypothesys on isf density
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::  zrhdtop_isf    ! density at bottom of ISF
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::  ziceload       ! density at bottom of ISF
      !!
      NAMELIST/namdyn_hpg/ ln_hpg_zco, ln_hpg_zps, ln_hpg_sco,     &
         &                 ln_hpg_djc, ln_hpg_prj, ln_hpg_isf
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namdyn_hpg in reference namelist : Hydrostatic pressure gradient
      READ  ( numnam_ref, namdyn_hpg, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdyn_hpg in reference namelist' )
      !
      REWIND( numnam_cfg )              ! Namelist namdyn_hpg in configuration namelist : Hydrostatic pressure gradient
      READ  ( numnam_cfg, namdyn_hpg, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdyn_hpg in configuration namelist' )
      IF(lwm) WRITE ( numond, namdyn_hpg )
      !
      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_hpg_init : hydrostatic pressure gradient initialisation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_hpg : choice of hpg scheme'
         WRITE(numout,*) '      z-coord. - full steps                             ln_hpg_zco    = ', ln_hpg_zco
         WRITE(numout,*) '      z-coord. - partial steps (interpolation)          ln_hpg_zps    = ', ln_hpg_zps
         WRITE(numout,*) '      s-coord. (standard jacobian formulation)          ln_hpg_sco    = ', ln_hpg_sco
         WRITE(numout,*) '      s-coord. (standard jacobian formulation) for isf  ln_hpg_isf    = ', ln_hpg_isf
         WRITE(numout,*) '      s-coord. (Density Jacobian: Cubic polynomial)     ln_hpg_djc    = ', ln_hpg_djc
         WRITE(numout,*) '      s-coord. (Pressure Jacobian: Cubic polynomial)    ln_hpg_prj    = ', ln_hpg_prj
      ENDIF
      !
      IF( ln_hpg_djc )   &
         &   CALL ctl_stop('dyn_hpg_init : Density Jacobian: Cubic polynominal method',   &
         &                 '   currently disabled (bugs under investigation).'        ,   &
         &                 '   Please select either  ln_hpg_sco or  ln_hpg_prj instead' )
         !
      IF( .NOT.ln_linssh .AND. .NOT.(ln_hpg_sco.OR.ln_hpg_prj.OR.ln_hpg_isf) )          &
         &   CALL ctl_stop('dyn_hpg_init : non-linear free surface requires either ',   &
         &                 '   the standard jacobian formulation hpg_sco    or '    ,   &
         &                 '   the pressure jacobian formulation hpg_prj'             )
         !
      IF( ln_hpg_isf ) THEN
         IF( .NOT. ln_isfcav )   CALL ctl_stop( ' hpg_isf not available if ln_isfcav = false ' )
       ELSE
         IF(       ln_isfcav )   CALL ctl_stop( 'Only hpg_isf has been corrected to work with ice shelf cavity.' )
      ENDIF
      !
      !                               ! Set nhpg from ln_hpg_... flags & consistency check
      nhpg   = np_ERROR
      ioptio = 0
      IF( ln_hpg_zco ) THEN   ;   nhpg = np_zco   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_zps ) THEN   ;   nhpg = np_zps   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_sco ) THEN   ;   nhpg = np_sco   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_djc ) THEN   ;   nhpg = np_djc   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_prj ) THEN   ;   nhpg = np_prj   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_isf ) THEN   ;   nhpg = np_isf   ;   ioptio = ioptio +1   ;   ENDIF
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'NO or several hydrostatic pressure gradient options used' )
      ! 
      IF(lwp) THEN
         WRITE(numout,*)
         SELECT CASE( nhpg )
         CASE( np_zco )   ;   WRITE(numout,*) '   ==>>>   z-coord. - full steps '
         CASE( np_zps )   ;   WRITE(numout,*) '   ==>>>   z-coord. - partial steps (interpolation)'
         CASE( np_sco )   ;   WRITE(numout,*) '   ==>>>   s-coord. (standard jacobian formulation)'
         CASE( np_djc )   ;   WRITE(numout,*) '   ==>>>   s-coord. (Density Jacobian: Cubic polynomial)'
         CASE( np_prj )   ;   WRITE(numout,*) '   ==>>>   s-coord. (Pressure Jacobian: Cubic polynomial)'
         CASE( np_isf )   ;   WRITE(numout,*) '   ==>>>   s-coord. (standard jacobian formulation) for isf'
         END SELECT
         WRITE(numout,*)
      ENDIF
      !                          
      IF ( .NOT. ln_isfcav ) THEN     !--- no ice shelf load
         riceload(:,:) = 0._wp
         !
      ELSE                            !--- set an ice shelf load
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ice shelf case: set the ice-shelf load'
         ALLOCATE( zts_top(jpi,jpj,jpts) , zrhd(jpi,jpj,jpk) , zrhdtop_isf(jpi,jpj) , ziceload(jpi,jpj) ) 
         !
         znad = 1._wp                     !- To use density and not density anomaly
         !
         !                                !- assume water displaced by the ice shelf is at T=-1.9 and S=34.4 (rude)
         zts_top(:,:,jp_tem) = -1.9_wp   ;   zts_top(:,:,jp_sal) = 34.4_wp
         !
         DO jk = 1, jpk                   !- compute density of the water displaced by the ice shelf 
            CALL eos( zts_top(:,:,:), gdept_n(:,:,jk), zrhd(:,:,jk) )
         END DO
         !
         !                                !- compute rhd at the ice/oce interface (ice shelf side)
         CALL eos( zts_top , risfdep, zrhdtop_isf )
         !
         !                                !- Surface value + ice shelf gradient
         ziceload = 0._wp                       ! compute pressure due to ice shelf load 
         DO jj = 1, jpj                         ! (used to compute hpgi/j for all the level from 1 to miku/v)
            DO ji = 1, jpi                      ! divided by 2 later
               ikt = mikt(ji,jj)
               ziceload(ji,jj) = ziceload(ji,jj) + (znad + zrhd(ji,jj,1) ) * e3w_n(ji,jj,1) * (1._wp - tmask(ji,jj,1))
               DO jk = 2, ikt-1
                  ziceload(ji,jj) = ziceload(ji,jj) + (2._wp * znad + zrhd(ji,jj,jk-1) + zrhd(ji,jj,jk)) * e3w_n(ji,jj,jk) &
                     &                              * (1._wp - tmask(ji,jj,jk))
               END DO
               IF (ikt  >=  2) ziceload(ji,jj) = ziceload(ji,jj) + (2._wp * znad + zrhdtop_isf(ji,jj) + zrhd(ji,jj,ikt-1)) &
                  &                                              * ( risfdep(ji,jj) - gdept_n(ji,jj,ikt-1) )
            END DO
         END DO
         riceload(:,:) = ziceload(:,:)  ! need to be saved for diaar5
         !
         DEALLOCATE( zts_top , zrhd , zrhdtop_isf , ziceload ) 
      ENDIF
      !
   END SUBROUTINE dyn_hpg_init


   SUBROUTINE hpg_zco( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_zco  ***
      !!
      !! ** Method  :   z-coordinate case, levels are horizontal surfaces.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level:    zhpi = grav .....
      !!                zhpj = grav .....
      !!      add it to the general momentum trend (ua,va).
      !!            ua = ua - 1/e1u * zhpi
      !!            va = va - 1/e2v * zhpj
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk       ! dummy loop indices
      REAL(wp) ::   zcoef0, zcoef1   ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zhpi, zhpj
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_zco : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   z-coordinate case '
      ENDIF

      zcoef0 = - grav * 0.5_wp      ! Local constant initialization

      ! Surface value
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zcoef1 = zcoef0 * e3w_n(ji,jj,1)
            ! hydrostatic pressure gradient
            zhpi(ji,jj,1) = zcoef1 * ( rhd(ji+1,jj,1) - rhd(ji,jj,1) ) * r1_e1u(ji,jj)
            zhpj(ji,jj,1) = zcoef1 * ( rhd(ji,jj+1,1) - rhd(ji,jj,1) ) * r1_e2v(ji,jj)
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
            va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)
         END DO
      END DO

      !
      ! interior value (2=<jk=<jpkm1)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zcoef1 = zcoef0 * e3w_n(ji,jj,jk)
               ! hydrostatic pressure gradient
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd(ji+1,jj,jk)+rhd(ji+1,jj,jk-1) )    &
                  &                       - ( rhd(ji  ,jj,jk)+rhd(ji  ,jj,jk-1) )  ) * r1_e1u(ji,jj)

               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd(ji,jj+1,jk)+rhd(ji,jj+1,jk-1) )    &
                  &                       - ( rhd(ji,jj,  jk)+rhd(ji,jj  ,jk-1) )  ) * r1_e2v(ji,jj)
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
            END DO
         END DO
      END DO
      !
   END SUBROUTINE hpg_zco


   SUBROUTINE hpg_zps( kt )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE hpg_zps  ***
      !!
      !! ** Method  :   z-coordinate plus partial steps case.  blahblah...
      !!
      !! ** Action  : - Update (ua,va) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk                       ! dummy loop indices
      INTEGER  ::   iku, ikv                         ! temporary integers
      REAL(wp) ::   zcoef0, zcoef1, zcoef2, zcoef3   ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zhpi, zhpj
      REAL(wp), DIMENSION(jpi,jpj) :: zgtsu, zgtsv, zgru, zgrv
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_zps : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   z-coordinate with partial steps - vector optimization'
      ENDIF

      ! Partial steps: Compute NOW horizontal gradient of t, s, rd at the last ocean level
      CALL zps_hde( kt, jpts, tsn, zgtsu, zgtsv, rhd, zgru , zgrv )

      ! Local constant initialization
      zcoef0 = - grav * 0.5_wp

      !  Surface value (also valid in partial step case)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zcoef1 = zcoef0 * e3w_n(ji,jj,1)
            ! hydrostatic pressure gradient
            zhpi(ji,jj,1) = zcoef1 * ( rhd(ji+1,jj  ,1) - rhd(ji,jj,1) ) * r1_e1u(ji,jj)
            zhpj(ji,jj,1) = zcoef1 * ( rhd(ji  ,jj+1,1) - rhd(ji,jj,1) ) * r1_e2v(ji,jj)
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
            va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)
         END DO
      END DO

      ! interior value (2=<jk=<jpkm1)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zcoef1 = zcoef0 * e3w_n(ji,jj,jk)
               ! hydrostatic pressure gradient
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) )   &
                  &                       - ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) )  ) * r1_e1u(ji,jj)

               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) )   &
                  &                       - ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) )  ) * r1_e2v(ji,jj)
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
            END DO
         END DO
      END DO

      ! partial steps correction at the last level  (use zgru & zgrv computed in zpshde.F90)
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            iku = mbku(ji,jj)
            ikv = mbkv(ji,jj)
            zcoef2 = zcoef0 * MIN( e3w_n(ji,jj,iku), e3w_n(ji+1,jj  ,iku) )
            zcoef3 = zcoef0 * MIN( e3w_n(ji,jj,ikv), e3w_n(ji  ,jj+1,ikv) )
            IF( iku > 1 ) THEN            ! on i-direction (level 2 or more)
               ua  (ji,jj,iku) = ua(ji,jj,iku) - zhpi(ji,jj,iku)         ! subtract old value
               zhpi(ji,jj,iku) = zhpi(ji,jj,iku-1)                   &   ! compute the new one
                  &            + zcoef2 * ( rhd(ji+1,jj,iku-1) - rhd(ji,jj,iku-1) + zgru(ji,jj) ) * r1_e1u(ji,jj)
               ua  (ji,jj,iku) = ua(ji,jj,iku) + zhpi(ji,jj,iku)         ! add the new one to the general momentum trend
            ENDIF
            IF( ikv > 1 ) THEN            ! on j-direction (level 2 or more)
               va  (ji,jj,ikv) = va(ji,jj,ikv) - zhpj(ji,jj,ikv)         ! subtract old value
               zhpj(ji,jj,ikv) = zhpj(ji,jj,ikv-1)                   &   ! compute the new one
                  &            + zcoef3 * ( rhd(ji,jj+1,ikv-1) - rhd(ji,jj,ikv-1) + zgrv(ji,jj) ) * r1_e2v(ji,jj)
               va  (ji,jj,ikv) = va(ji,jj,ikv) + zhpj(ji,jj,ikv)         ! add the new one to the general momentum trend
            ENDIF
         END DO
      END DO
      !
   END SUBROUTINE hpg_zps


   SUBROUTINE hpg_sco( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_sco  ***
      !!
      !! ** Method  :   s-coordinate case. Jacobian scheme.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level. s-coordinates (ln_sco): a corrective term is added
      !!      to the horizontal pressure gradient :
      !!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
      !!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
      !!      add it to the general momentum trend (ua,va).
      !!         ua = ua - 1/e1u * zhpi
      !!         va = va - 1/e2v * zhpj
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk, jii, jjj                 ! dummy loop indices
      REAL(wp) ::   zcoef0, zuap, zvap, znad, ztmp       ! temporary scalars
      LOGICAL  ::   ll_tmp1, ll_tmp2                     ! local logical variables
      REAL(wp), DIMENSION(jpi,jpj,jpk)      ::   zhpi, zhpj
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zcpx, zcpy   !W/D pressure filter
      !!----------------------------------------------------------------------
      !
      IF( ln_wd_il ) ALLOCATE(zcpx(jpi,jpj), zcpy(jpi,jpj))
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_sco : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, OPA original scheme used'
      ENDIF
      !
      zcoef0 = - grav * 0.5_wp
      IF ( ln_linssh ) THEN   ;   znad = 0._wp         ! Fixed    volume: density anomaly
      ELSE                    ;   znad = 1._wp         ! Variable volume: density
      ENDIF
      !
      IF( ln_wd_il ) THEN
        DO jj = 2, jpjm1
           DO ji = 2, jpim1 
             ll_tmp1 = MIN(  sshn(ji,jj)               ,  sshn(ji+1,jj) ) >                &
                  &    MAX( -ht_0(ji,jj)               , -ht_0(ji+1,jj) ) .AND.            &
                  &    MAX(  sshn(ji,jj) +  ht_0(ji,jj),  sshn(ji+1,jj) + ht_0(ji+1,jj) )  &
                  &                                                       > rn_wdmin1 + rn_wdmin2
             ll_tmp2 = ( ABS( sshn(ji,jj)              -  sshn(ji+1,jj) ) > 1.E-12 ) .AND. (       &
                  &    MAX(   sshn(ji,jj)              ,  sshn(ji+1,jj) ) >                &
                  &    MAX(  -ht_0(ji,jj)              , -ht_0(ji+1,jj) ) + rn_wdmin1 + rn_wdmin2 )

             IF(ll_tmp1) THEN
               zcpx(ji,jj) = 1.0_wp
             ELSE IF(ll_tmp2) THEN
               ! no worries about  sshn(ji+1,jj) -  sshn(ji  ,jj) = 0, it won't happen ! here
               zcpx(ji,jj) = ABS( (sshn(ji+1,jj) + ht_0(ji+1,jj) - sshn(ji,jj) - ht_0(ji,jj)) &
                           &    / (sshn(ji+1,jj) - sshn(ji  ,jj)) )
             ELSE
               zcpx(ji,jj) = 0._wp
             END IF
      
             ll_tmp1 = MIN(  sshn(ji,jj)              ,  sshn(ji,jj+1) ) >                &
                  &    MAX( -ht_0(ji,jj)              , -ht_0(ji,jj+1) ) .AND.            &
                  &    MAX(  sshn(ji,jj) + ht_0(ji,jj),  sshn(ji,jj+1) + ht_0(ji,jj+1) )  &
                  &                                                      > rn_wdmin1 + rn_wdmin2
             ll_tmp2 = ( ABS( sshn(ji,jj)             -  sshn(ji,jj+1) ) > 1.E-12 ) .AND. (        &
                  &    MAX(   sshn(ji,jj)             ,  sshn(ji,jj+1) ) >                &
                  &    MAX(  -ht_0(ji,jj)             , -ht_0(ji,jj+1) ) + rn_wdmin1 + rn_wdmin2 )

             IF(ll_tmp1) THEN
               zcpy(ji,jj) = 1.0_wp
             ELSE IF(ll_tmp2) THEN
               ! no worries about  sshn(ji,jj+1) -  sshn(ji,jj  ) = 0, it won't happen ! here
               zcpy(ji,jj) = ABS( (sshn(ji,jj+1) + ht_0(ji,jj+1) - sshn(ji,jj) - ht_0(ji,jj)) &
                           &    / (sshn(ji,jj+1) - sshn(ji,jj  )) )
             ELSE
               zcpy(ji,jj) = 0._wp
             END IF
           END DO
        END DO
        CALL lbc_lnk_multi( 'dynhpg', zcpx, 'U', 1., zcpy, 'V', 1. )
      END IF

    IF (ln_2d) THEN 
        ! Surface value 
        DO jj = 2, jpjm1 
           DO ji = fs_2, fs_jpim1   ! vector opt. 
              ! hydrostatic pressure gradient along s-surfaces 
              zhpi(ji,jj,1) = zcoef0 * (  e3w_n(ji+1,jj  ,1) * znad   & 
                 &                      - e3w_n(ji  ,jj  ,1) * znad ) * r1_e1u(ji,jj) 
              zhpj(ji,jj,1) = zcoef0 * (  e3w_n(ji  ,jj+1,1) * znad   & 
                 &                      - e3w_n(ji  ,jj  ,1) * znad ) * r1_e2v(ji,jj) 
              ! s-coordinate pressure gradient correction 
              zuap = -zcoef0 * ( 2._wp * znad )   & 
                 &           * ( gde3w_n(ji+1,jj,1) - gde3w_n(ji,jj,1) ) * r1_e1u(ji,jj) 
              zvap = -zcoef0 * ( 2._wp * znad )   & 
                 &           * ( gde3w_n(ji,jj+1,1) - gde3w_n(ji,jj,1) ) * r1_e2v(ji,jj) 
 
 
              IF( ln_wd_il ) THEN 
 
                zhpi(ji,jj,1) = zhpi(ji,jj,1) * zcpx(ji,jj) 
                zhpj(ji,jj,1) = zhpj(ji,jj,1) * zcpy(ji,jj)  
                zuap = zuap * zcpx(ji,jj) 
                zvap = zvap * zcpy(ji,jj) 
              ENDIF 
 
              ! add to the general momentum trend 
              ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1) + zuap * umask(ji,jj,1) 
              va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1) + zvap * vmask(ji,jj,1) 
           END DO 
        END DO 
 
        ! interior value (2=<jk=<jpkm1) 
        DO jk = 2, jpkm1 
           DO jj = 2, jpjm1 
              DO ji = fs_2, fs_jpim1   ! vector opt. 
                 ! hydrostatic pressure gradient along s-surfaces 
                 zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) + zcoef0 * r1_e1u(ji,jj)   & 
                    &           * (  e3w_n(ji+1,jj,jk) * ( 2*znad )   & 
                    &              - e3w_n(ji  ,jj,jk) * ( 2*znad )  ) 
                 zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1) + zcoef0 * r1_e2v(ji,jj)   & 
                    &           * (  e3w_n(ji,jj+1,jk) * ( 2*znad )   & 
                    &              - e3w_n(ji,jj  ,jk) * ( 2*znad )  ) 
                 ! s-coordinate pressure gradient correction 
                 zuap = -zcoef0 * ( 2._wp * znad )   & 
                    &           * ( gde3w_n(ji+1,jj  ,jk) - gde3w_n(ji,jj,jk) ) * r1_e1u(ji,jj) 
                 zvap = -zcoef0 * ( 2._wp * znad )   & 
                    &           * ( gde3w_n(ji  ,jj+1,jk) - gde3w_n(ji,jj,jk) ) * r1_e2v(ji,jj) 
 
                 IF( ln_wd_il ) THEN 
                   zhpi(ji,jj,jk) = zhpi(ji,jj,jk) * zcpx(ji,jj) 
                   zhpj(ji,jj,jk) = zhpj(ji,jj,jk) * zcpy(ji,jj)  
                   zuap = zuap * zcpx(ji,jj) 
                   zvap = zvap * zcpy(ji,jj) 
                 ENDIF 
 
                 ! add to the general momentum trend 
                 ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk) + zuap * umask(ji,jj,1) 
                 va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk) + zvap * vmask(ji,jj,1) 
              END DO 
           END DO 
        END DO 
 
    ELSE 

      ! Surface value
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ! hydrostatic pressure gradient along s-surfaces
            zhpi(ji,jj,1) = zcoef0 * (  e3w_n(ji+1,jj  ,1) * ( znad + rhd(ji+1,jj  ,1) )    &
               &                      - e3w_n(ji  ,jj  ,1) * ( znad + rhd(ji  ,jj  ,1) )  ) * r1_e1u(ji,jj)
            zhpj(ji,jj,1) = zcoef0 * (  e3w_n(ji  ,jj+1,1) * ( znad + rhd(ji  ,jj+1,1) )    &
               &                      - e3w_n(ji  ,jj  ,1) * ( znad + rhd(ji  ,jj  ,1) )  ) * r1_e2v(ji,jj)
            ! s-coordinate pressure gradient correction
            zuap = -zcoef0 * ( rhd    (ji+1,jj,1) + rhd    (ji,jj,1) + 2._wp * znad )   &
               &           * ( gde3w_n(ji+1,jj,1) - gde3w_n(ji,jj,1) ) * r1_e1u(ji,jj)
            zvap = -zcoef0 * ( rhd    (ji,jj+1,1) + rhd    (ji,jj,1) + 2._wp * znad )   &
               &           * ( gde3w_n(ji,jj+1,1) - gde3w_n(ji,jj,1) ) * r1_e2v(ji,jj)
            !
            IF( ln_wd_il ) THEN
               zhpi(ji,jj,1) = zhpi(ji,jj,1) * zcpx(ji,jj)
               zhpj(ji,jj,1) = zhpj(ji,jj,1) * zcpy(ji,jj) 
               zuap = zuap * zcpx(ji,jj)
               zvap = zvap * zcpy(ji,jj)
            ENDIF
            !
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1) + zuap
            va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1) + zvap
         END DO
      END DO

      ! interior value (2=<jk=<jpkm1)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! hydrostatic pressure gradient along s-surfaces
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) + zcoef0 * r1_e1u(ji,jj)   &
                  &           * (  e3w_n(ji+1,jj,jk) * ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) + 2*znad )   &
                  &              - e3w_n(ji  ,jj,jk) * ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) + 2*znad )  )
               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1) + zcoef0 * r1_e2v(ji,jj)   &
                  &           * (  e3w_n(ji,jj+1,jk) * ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) + 2*znad )   &
                  &              - e3w_n(ji,jj  ,jk) * ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) + 2*znad )  )
               ! s-coordinate pressure gradient correction
               zuap = -zcoef0 * ( rhd    (ji+1,jj  ,jk) + rhd    (ji,jj,jk) + 2._wp * znad )   &
                  &           * ( gde3w_n(ji+1,jj  ,jk) - gde3w_n(ji,jj,jk) ) * r1_e1u(ji,jj)
               zvap = -zcoef0 * ( rhd    (ji  ,jj+1,jk) + rhd    (ji,jj,jk) + 2._wp * znad )   &
                  &           * ( gde3w_n(ji  ,jj+1,jk) - gde3w_n(ji,jj,jk) ) * r1_e2v(ji,jj)
               !
               IF( ln_wd_il ) THEN
                  zhpi(ji,jj,jk) = zhpi(ji,jj,jk) * zcpx(ji,jj)
                  zhpj(ji,jj,jk) = zhpj(ji,jj,jk) * zcpy(ji,jj) 
                  zuap = zuap * zcpx(ji,jj)
                  zvap = zvap * zcpy(ji,jj)
               ENDIF
               !
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk) + zuap
               va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk) + zvap
            END DO
         END DO
      END DO

    ENDIF !ln_2d
      !
      IF( ln_wd_il )  DEALLOCATE( zcpx , zcpy )
      !
   END SUBROUTINE hpg_sco


   SUBROUTINE hpg_isf( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_isf  ***
      !!
      !! ** Method  :   s-coordinate case. Jacobian scheme.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level. s-coordinates (ln_sco): a corrective term is added
      !!      to the horizontal pressure gradient :
      !!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
      !!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
      !!      add it to the general momentum trend (ua,va).
      !!         ua = ua - 1/e1u * zhpi
      !!         va = va - 1/e2v * zhpj
      !!      iceload is added and partial cell case are added to the top and bottom
      !!      
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk, ikt, iktp1i, iktp1j   ! dummy loop indices
      REAL(wp) ::   zcoef0, zuap, zvap, znad          ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk ) ::  zhpi, zhpj
      REAL(wp), DIMENSION(jpi,jpj,jpts) ::  zts_top
      REAL(wp), DIMENSION(jpi,jpj)      ::  zrhdtop_oce
      !!----------------------------------------------------------------------
      !
      zcoef0 = - grav * 0.5_wp   ! Local constant initialization
      !
      znad=1._wp                 ! To use density and not density anomaly
      !
      !                          ! iniitialised to 0. zhpi zhpi 
      zhpi(:,:,:) = 0._wp   ;   zhpj(:,:,:) = 0._wp

      ! compute rhd at the ice/oce interface (ocean side)
      ! usefull to reduce residual current in the test case ISOMIP with no melting
      DO ji = 1, jpi
        DO jj = 1, jpj
          ikt = mikt(ji,jj)
          zts_top(ji,jj,1) = tsn(ji,jj,ikt,1)
          zts_top(ji,jj,2) = tsn(ji,jj,ikt,2)
        END DO
      END DO
      CALL eos( zts_top, risfdep, zrhdtop_oce )

!==================================================================================     
!===== Compute surface value ===================================================== 
!==================================================================================
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ikt    = mikt(ji,jj)
            iktp1i = mikt(ji+1,jj)
            iktp1j = mikt(ji,jj+1)
            ! hydrostatic pressure gradient along s-surfaces and ice shelf pressure
            ! we assume ISF is in isostatic equilibrium
            zhpi(ji,jj,1) = zcoef0 / e1u(ji,jj) * ( 0.5_wp * e3w_n(ji+1,jj,iktp1i)                                    &
               &                                    * ( 2._wp * znad + rhd(ji+1,jj,iktp1i) + zrhdtop_oce(ji+1,jj) )   &
               &                                  - 0.5_wp * e3w_n(ji,jj,ikt)                                         &
               &                                    * ( 2._wp * znad + rhd(ji,jj,ikt) + zrhdtop_oce(ji,jj) )          &
               &                                  + ( riceload(ji+1,jj) - riceload(ji,jj))                            ) 
            zhpj(ji,jj,1) = zcoef0 / e2v(ji,jj) * ( 0.5_wp * e3w_n(ji,jj+1,iktp1j)                                    &
               &                                    * ( 2._wp * znad + rhd(ji,jj+1,iktp1j) + zrhdtop_oce(ji,jj+1) )   &
               &                                  - 0.5_wp * e3w_n(ji,jj,ikt)                                         & 
               &                                    * ( 2._wp * znad + rhd(ji,jj,ikt) + zrhdtop_oce(ji,jj) )          &
               &                                  + ( riceload(ji,jj+1) - riceload(ji,jj))                            ) 
            ! s-coordinate pressure gradient correction (=0 if z coordinate)
            zuap = -zcoef0 * ( rhd    (ji+1,jj,1) + rhd    (ji,jj,1) + 2._wp * znad )   &
               &           * ( gde3w_n(ji+1,jj,1) - gde3w_n(ji,jj,1) ) * r1_e1u(ji,jj)
            zvap = -zcoef0 * ( rhd    (ji,jj+1,1) + rhd    (ji,jj,1) + 2._wp * znad )   &
               &           * ( gde3w_n(ji,jj+1,1) - gde3w_n(ji,jj,1) ) * r1_e2v(ji,jj)
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + (zhpi(ji,jj,1) + zuap) * umask(ji,jj,1)
            va(ji,jj,1) = va(ji,jj,1) + (zhpj(ji,jj,1) + zvap) * vmask(ji,jj,1)
         END DO
      END DO
!==================================================================================     
!===== Compute interior value ===================================================== 
!==================================================================================
      ! interior value (2=<jk=<jpkm1)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! hydrostatic pressure gradient along s-surfaces
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) + zcoef0 / e1u(ji,jj)   &
                  &           * (  e3w_n(ji+1,jj,jk) * ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) + 2*znad ) * wmask(ji+1,jj,jk)   &
                  &              - e3w_n(ji  ,jj,jk) * ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) + 2*znad ) * wmask(ji  ,jj,jk)   )
               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1) + zcoef0 / e2v(ji,jj)   &
                  &           * (  e3w_n(ji,jj+1,jk) * ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) + 2*znad ) * wmask(ji,jj+1,jk)   &
                  &              - e3w_n(ji,jj  ,jk) * ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) + 2*znad ) * wmask(ji,jj  ,jk)   )
               ! s-coordinate pressure gradient correction
               zuap = -zcoef0 * ( rhd   (ji+1,jj  ,jk) + rhd   (ji,jj,jk) + 2._wp * znad )   &
                  &           * ( gde3w_n(ji+1,jj  ,jk) - gde3w_n(ji,jj,jk) ) / e1u(ji,jj)
               zvap = -zcoef0 * ( rhd   (ji  ,jj+1,jk) + rhd   (ji,jj,jk) + 2._wp * znad )   &
                  &           * ( gde3w_n(ji  ,jj+1,jk) - gde3w_n(ji,jj,jk) ) / e2v(ji,jj)
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + (zhpi(ji,jj,jk) + zuap) * umask(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) + (zhpj(ji,jj,jk) + zvap) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
   END SUBROUTINE hpg_isf


   SUBROUTINE hpg_djc( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_djc  ***
      !!
      !! ** Method  :   Density Jacobian with Cubic polynomial scheme
      !!
      !! Reference: Shchepetkin and McWilliams, J. Geophys. Res., 108(C3), 3090, 2003
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk          ! dummy loop indices
      REAL(wp) ::   zcoef0, zep, cffw   ! temporary scalars
      REAL(wp) ::   z1_10, cffu, cffx   !    "         "
      REAL(wp) ::   z1_12, cffv, cffy   !    "         "
      LOGICAL  ::   ll_tmp1, ll_tmp2    ! local logical variables
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zhpi, zhpj
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   dzx, dzy, dzz, dzu, dzv, dzw
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   drhox, drhoy, drhoz, drhou, drhov, drhow
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   rho_i, rho_j, rho_k
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zcpx, zcpy   !W/D pressure filter
      !!----------------------------------------------------------------------
      !
      IF( ln_wd_il ) THEN
         ALLOCATE( zcpx(jpi,jpj) , zcpy(jpi,jpj) )
        DO jj = 2, jpjm1
           DO ji = 2, jpim1 
             ll_tmp1 = MIN(  sshn(ji,jj)              ,  sshn(ji+1,jj) ) >                &
                  &    MAX( -ht_0(ji,jj)              , -ht_0(ji+1,jj) ) .AND.            &
                  &    MAX(  sshn(ji,jj) + ht_0(ji,jj),  sshn(ji+1,jj) + ht_0(ji+1,jj) )  &
                  &                                                      > rn_wdmin1 + rn_wdmin2
             ll_tmp2 = ( ABS( sshn(ji,jj)             -  sshn(ji+1,jj) ) > 1.E-12 ) .AND. (        &
                  &    MAX(   sshn(ji,jj)             ,  sshn(ji+1,jj) ) >                &
                  &    MAX(  -ht_0(ji,jj)             , -ht_0(ji+1,jj) ) + rn_wdmin1 + rn_wdmin2 )
             IF(ll_tmp1) THEN
               zcpx(ji,jj) = 1.0_wp
             ELSE IF(ll_tmp2) THEN
               ! no worries about  sshn(ji+1,jj) -  sshn(ji  ,jj) = 0, it won't happen ! here
               zcpx(ji,jj) = ABS( (sshn(ji+1,jj) + ht_0(ji+1,jj) - sshn(ji,jj) - ht_0(ji,jj)) &
                           &    / (sshn(ji+1,jj) - sshn(ji  ,jj)) )
             ELSE
               zcpx(ji,jj) = 0._wp
             END IF
      
             ll_tmp1 = MIN(  sshn(ji,jj)              ,  sshn(ji,jj+1) ) >                &
                  &    MAX( -ht_0(ji,jj)              , -ht_0(ji,jj+1) ) .AND.            &
                  &    MAX(  sshn(ji,jj) + ht_0(ji,jj),  sshn(ji,jj+1) + ht_0(ji,jj+1) )  &
                  &                                                      > rn_wdmin1 + rn_wdmin2
             ll_tmp2 = ( ABS( sshn(ji,jj)             -  sshn(ji,jj+1) ) > 1.E-12 ) .AND. (        &
                  &    MAX(   sshn(ji,jj)             ,  sshn(ji,jj+1) ) >                &
                  &    MAX(  -ht_0(ji,jj)             , -ht_0(ji,jj+1) ) + rn_wdmin1 + rn_wdmin2 )

             IF(ll_tmp1) THEN
               zcpy(ji,jj) = 1.0_wp
             ELSE IF(ll_tmp2) THEN
               ! no worries about  sshn(ji,jj+1) -  sshn(ji,jj  ) = 0, it won't happen ! here
               zcpy(ji,jj) = ABS( (sshn(ji,jj+1) + ht_0(ji,jj+1) - sshn(ji,jj) - ht_0(ji,jj)) &
                           &    / (sshn(ji,jj+1) - sshn(ji,jj  )) )
             ELSE
               zcpy(ji,jj) = 0._wp
             END IF
           END DO
        END DO
        CALL lbc_lnk_multi( 'dynhpg', zcpx, 'U', 1., zcpy, 'V', 1. )
      END IF

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_djc : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, density Jacobian with cubic polynomial scheme'
      ENDIF

      ! Local constant initialization
      zcoef0 = - grav * 0.5_wp
      z1_10  = 1._wp / 10._wp
      z1_12  = 1._wp / 12._wp

      !----------------------------------------------------------------------------------------
      !  compute and store in provisional arrays elementary vertical and horizontal differences
      !----------------------------------------------------------------------------------------

!!bug gm   Not a true bug, but... dzz=e3w  for dzx, dzy verify what it is really

      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               drhoz(ji,jj,jk) = rhd    (ji  ,jj  ,jk) - rhd    (ji,jj,jk-1)
               dzz  (ji,jj,jk) = gde3w_n(ji  ,jj  ,jk) - gde3w_n(ji,jj,jk-1)
               drhox(ji,jj,jk) = rhd    (ji+1,jj  ,jk) - rhd    (ji,jj,jk  )
               dzx  (ji,jj,jk) = gde3w_n(ji+1,jj  ,jk) - gde3w_n(ji,jj,jk  )
               drhoy(ji,jj,jk) = rhd    (ji  ,jj+1,jk) - rhd    (ji,jj,jk  )
               dzy  (ji,jj,jk) = gde3w_n(ji  ,jj+1,jk) - gde3w_n(ji,jj,jk  )
            END DO
         END DO
      END DO

      !-------------------------------------------------------------------------
      ! compute harmonic averages using eq. 5.18
      !-------------------------------------------------------------------------
      zep = 1.e-15

!!bug  gm  drhoz not defined at level 1 and used (jk-1 with jk=2)
!!bug  gm  idem for drhox, drhoy et ji=jpi and jj=jpj

      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               cffw = 2._wp * drhoz(ji  ,jj  ,jk) * drhoz(ji,jj,jk-1)

               cffu = 2._wp * drhox(ji+1,jj  ,jk) * drhox(ji,jj,jk  )
               cffx = 2._wp * dzx  (ji+1,jj  ,jk) * dzx  (ji,jj,jk  )

               cffv = 2._wp * drhoy(ji  ,jj+1,jk) * drhoy(ji,jj,jk  )
               cffy = 2._wp * dzy  (ji  ,jj+1,jk) * dzy  (ji,jj,jk  )

               IF( cffw > zep) THEN
                  drhow(ji,jj,jk) = 2._wp *   drhoz(ji,jj,jk) * drhoz(ji,jj,jk-1)   &
                     &                    / ( drhoz(ji,jj,jk) + drhoz(ji,jj,jk-1) )
               ELSE
                  drhow(ji,jj,jk) = 0._wp
               ENDIF

               dzw(ji,jj,jk) = 2._wp *   dzz(ji,jj,jk) * dzz(ji,jj,jk-1)   &
                  &                  / ( dzz(ji,jj,jk) + dzz(ji,jj,jk-1) )

               IF( cffu > zep ) THEN
                  drhou(ji,jj,jk) = 2._wp *   drhox(ji+1,jj,jk) * drhox(ji,jj,jk)   &
                     &                    / ( drhox(ji+1,jj,jk) + drhox(ji,jj,jk) )
               ELSE
                  drhou(ji,jj,jk ) = 0._wp
               ENDIF

               IF( cffx > zep ) THEN
                  dzu(ji,jj,jk) = 2._wp *   dzx(ji+1,jj,jk) * dzx(ji,jj,jk)   &
                     &                  / ( dzx(ji+1,jj,jk) + dzx(ji,jj,jk) )
               ELSE
                  dzu(ji,jj,jk) = 0._wp
               ENDIF

               IF( cffv > zep ) THEN
                  drhov(ji,jj,jk) = 2._wp *   drhoy(ji,jj+1,jk) * drhoy(ji,jj,jk)   &
                     &                    / ( drhoy(ji,jj+1,jk) + drhoy(ji,jj,jk) )
               ELSE
                  drhov(ji,jj,jk) = 0._wp
               ENDIF

               IF( cffy > zep ) THEN
                  dzv(ji,jj,jk) = 2._wp *   dzy(ji,jj+1,jk) * dzy(ji,jj,jk)   &
                     &                  / ( dzy(ji,jj+1,jk) + dzy(ji,jj,jk) )
               ELSE
                  dzv(ji,jj,jk) = 0._wp
               ENDIF

            END DO
         END DO
      END DO

      !----------------------------------------------------------------------------------
      ! apply boundary conditions at top and bottom using 5.36-5.37
      !----------------------------------------------------------------------------------
      drhow(:,:, 1 ) = 1.5_wp * ( drhoz(:,:, 2 ) - drhoz(:,:,  1  ) ) - 0.5_wp * drhow(:,:,  2  )
      drhou(:,:, 1 ) = 1.5_wp * ( drhox(:,:, 2 ) - drhox(:,:,  1  ) ) - 0.5_wp * drhou(:,:,  2  )
      drhov(:,:, 1 ) = 1.5_wp * ( drhoy(:,:, 2 ) - drhoy(:,:,  1  ) ) - 0.5_wp * drhov(:,:,  2  )

      drhow(:,:,jpk) = 1.5_wp * ( drhoz(:,:,jpk) - drhoz(:,:,jpkm1) ) - 0.5_wp * drhow(:,:,jpkm1)
      drhou(:,:,jpk) = 1.5_wp * ( drhox(:,:,jpk) - drhox(:,:,jpkm1) ) - 0.5_wp * drhou(:,:,jpkm1)
      drhov(:,:,jpk) = 1.5_wp * ( drhoy(:,:,jpk) - drhoy(:,:,jpkm1) ) - 0.5_wp * drhov(:,:,jpkm1)


      !--------------------------------------------------------------
      ! Upper half of top-most grid box, compute and store
      !-------------------------------------------------------------

!!bug gm   :  e3w-gde3w = 0.5*e3w  ....  and gde3w(2)-gde3w(1)=e3w(2) ....   to be verified
!          true if gde3w is really defined as the sum of the e3w scale factors as, it seems to me, it should be

      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            rho_k(ji,jj,1) = -grav * ( e3w_n(ji,jj,1) - gde3w_n(ji,jj,1) )               &
               &                   * (  rhd(ji,jj,1)                                     &
               &                     + 0.5_wp * ( rhd    (ji,jj,2) - rhd    (ji,jj,1) )  &
               &                              * ( e3w_n  (ji,jj,1) - gde3w_n(ji,jj,1) )  &
               &                              / ( gde3w_n(ji,jj,2) - gde3w_n(ji,jj,1) )  )
         END DO
      END DO

!!bug gm    : here also, simplification is possible
!!bug gm    : optimisation: 1/10 and 1/12 the division should be done before the loop

      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.

               rho_k(ji,jj,jk) = zcoef0 * ( rhd    (ji,jj,jk) + rhd    (ji,jj,jk-1) )                                   &
                  &                     * ( gde3w_n(ji,jj,jk) - gde3w_n(ji,jj,jk-1) )                                   &
                  &            - grav * z1_10 * (                                                                   &
                  &     ( drhow  (ji,jj,jk) - drhow  (ji,jj,jk-1) )                                                     &
                  &   * ( gde3w_n(ji,jj,jk) - gde3w_n(ji,jj,jk-1) - z1_12 * ( dzw  (ji,jj,jk) + dzw  (ji,jj,jk-1) ) )   &
                  &   - ( dzw    (ji,jj,jk) - dzw    (ji,jj,jk-1) )                                                     &
                  &   * ( rhd    (ji,jj,jk) - rhd    (ji,jj,jk-1) - z1_12 * ( drhow(ji,jj,jk) + drhow(ji,jj,jk-1) ) )   &
                  &                             )

               rho_i(ji,jj,jk) = zcoef0 * ( rhd    (ji+1,jj,jk) + rhd    (ji,jj,jk) )                                   &
                  &                     * ( gde3w_n(ji+1,jj,jk) - gde3w_n(ji,jj,jk) )                                   &
                  &            - grav* z1_10 * (                                                                    &
                  &     ( drhou  (ji+1,jj,jk) - drhou  (ji,jj,jk) )                                                     &
                  &   * ( gde3w_n(ji+1,jj,jk) - gde3w_n(ji,jj,jk) - z1_12 * ( dzu  (ji+1,jj,jk) + dzu  (ji,jj,jk) ) )   &
                  &   - ( dzu    (ji+1,jj,jk) - dzu    (ji,jj,jk) )                                                     &
                  &   * ( rhd    (ji+1,jj,jk) - rhd    (ji,jj,jk) - z1_12 * ( drhou(ji+1,jj,jk) + drhou(ji,jj,jk) ) )   &
                  &                            )

               rho_j(ji,jj,jk) = zcoef0 * ( rhd    (ji,jj+1,jk) + rhd    (ji,jj,jk) )                                 &
                  &                     * ( gde3w_n(ji,jj+1,jk) - gde3w_n(ji,jj,jk) )                                   &
                  &            - grav* z1_10 * (                                                                    &
                  &     ( drhov  (ji,jj+1,jk) - drhov  (ji,jj,jk) )                                                     &
                  &   * ( gde3w_n(ji,jj+1,jk) - gde3w_n(ji,jj,jk) - z1_12 * ( dzv  (ji,jj+1,jk) + dzv  (ji,jj,jk) ) )   &
                  &   - ( dzv    (ji,jj+1,jk) - dzv    (ji,jj,jk) )                                                     &
                  &   * ( rhd    (ji,jj+1,jk) - rhd    (ji,jj,jk) - z1_12 * ( drhov(ji,jj+1,jk) + drhov(ji,jj,jk) ) )   &
                  &                            )

            END DO
         END DO
      END DO
      CALL lbc_lnk_multi( 'dynhpg', rho_k, 'W', 1., rho_i, 'U', 1., rho_j, 'V', 1. )

      ! ---------------
      !  Surface value
      ! ---------------
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zhpi(ji,jj,1) = ( rho_k(ji+1,jj  ,1) - rho_k(ji,jj,1) - rho_i(ji,jj,1) ) * r1_e1u(ji,jj)
            zhpj(ji,jj,1) = ( rho_k(ji  ,jj+1,1) - rho_k(ji,jj,1) - rho_j(ji,jj,1) ) * r1_e2v(ji,jj)
            IF( ln_wd_il ) THEN
              zhpi(ji,jj,1) = zhpi(ji,jj,1) * zcpx(ji,jj)
              zhpj(ji,jj,1) = zhpj(ji,jj,1) * zcpy(ji,jj) 
            ENDIF
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
            va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)
         END DO
      END DO

      ! ----------------
      !  interior value   (2=<jk=<jpkm1)
      ! ----------------
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! hydrostatic pressure gradient along s-surfaces
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)                                &
                  &           + (  ( rho_k(ji+1,jj,jk) - rho_k(ji,jj,jk  ) )    &
                  &              - ( rho_i(ji  ,jj,jk) - rho_i(ji,jj,jk-1) )  ) * r1_e1u(ji,jj)
               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)                                &
                  &           + (  ( rho_k(ji,jj+1,jk) - rho_k(ji,jj,jk  ) )    &
                  &               -( rho_j(ji,jj  ,jk) - rho_j(ji,jj,jk-1) )  ) * r1_e2v(ji,jj)
               IF( ln_wd_il ) THEN
                 zhpi(ji,jj,jk) = zhpi(ji,jj,jk) * zcpx(ji,jj)
                 zhpj(ji,jj,jk) = zhpj(ji,jj,jk) * zcpy(ji,jj) 
               ENDIF
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      IF( ln_wd_il )   DEALLOCATE( zcpx, zcpy )
      !
   END SUBROUTINE hpg_djc


   SUBROUTINE hpg_prj( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_prj  ***
      !!
      !! ** Method  :   s-coordinate case.
      !!      A Pressure-Jacobian horizontal pressure gradient method
      !!      based on the constrained cubic-spline interpolation for
      !!      all vertical coordinate systems
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      INTEGER, PARAMETER  :: polynomial_type = 1    ! 1: cubic spline, 2: linear
      INTEGER, INTENT(in) ::   kt                   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk, jkk                 ! dummy loop indices
      REAL(wp) ::   zcoef0, znad                    ! local scalars
      !
      !! The local variables for the correction term
      INTEGER  :: jk1, jis, jid, jjs, jjd
      LOGICAL  :: ll_tmp1, ll_tmp2                  ! local logical variables
      REAL(wp) :: zuijk, zvijk, zpwes, zpwed, zpnss, zpnsd, zdeps
      REAL(wp) :: zrhdt1
      REAL(wp) :: zdpdx1, zdpdx2, zdpdy1, zdpdy2
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zdept, zrhh
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zhpi, zu, zv, fsp, xsp, asp, bsp, csp, dsp
      REAL(wp), DIMENSION(jpi,jpj)   ::   zsshu_n, zsshv_n
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zcpx, zcpy   !W/D pressure filter
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_prj : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, cubic spline pressure Jacobian'
      ENDIF

      ! Local constant initialization
      zcoef0 = - grav
      znad = 1._wp
      IF( ln_linssh )   znad = 0._wp

      IF( ln_wd_il ) THEN
         ALLOCATE( zcpx(jpi,jpj) , zcpy(jpi,jpj) )
         DO jj = 2, jpjm1
           DO ji = 2, jpim1 
             ll_tmp1 = MIN(  sshn(ji,jj)              ,  sshn(ji+1,jj) ) >                &
                  &    MAX( -ht_0(ji,jj)              , -ht_0(ji+1,jj) ) .AND.            &
                  &    MAX(  sshn(ji,jj) + ht_0(ji,jj),  sshn(ji+1,jj) + ht_0(ji+1,jj) )  &
                  &                                                      > rn_wdmin1 + rn_wdmin2
             ll_tmp2 = ( ABS( sshn(ji,jj)             -  sshn(ji+1,jj) ) > 1.E-12 ) .AND. (         &
                  &    MAX(   sshn(ji,jj)             ,  sshn(ji+1,jj) ) >                &
                  &    MAX(  -ht_0(ji,jj)             , -ht_0(ji+1,jj) ) + rn_wdmin1 + rn_wdmin2 )

             IF(ll_tmp1) THEN
               zcpx(ji,jj) = 1.0_wp
             ELSE IF(ll_tmp2) THEN
               ! no worries about  sshn(ji+1,jj) -  sshn(ji  ,jj) = 0, it won't happen ! here
               zcpx(ji,jj) = ABS( (sshn(ji+1,jj) + ht_0(ji+1,jj) - sshn(ji,jj) - ht_0(ji,jj)) &
                           &    / (sshn(ji+1,jj) -  sshn(ji  ,jj)) )
              
                zcpx(ji,jj) = max(min( zcpx(ji,jj) , 1.0_wp),0.0_wp)
             ELSE
               zcpx(ji,jj) = 0._wp
             END IF
      
             ll_tmp1 = MIN(  sshn(ji,jj)              ,  sshn(ji,jj+1) ) >                &
                  &    MAX( -ht_0(ji,jj)              , -ht_0(ji,jj+1) ) .AND.            &
                  &    MAX(  sshn(ji,jj) + ht_0(ji,jj),  sshn(ji,jj+1) + ht_0(ji,jj+1) )  &
                  &                                                      > rn_wdmin1 + rn_wdmin2
             ll_tmp2 = ( ABS( sshn(ji,jj)             -  sshn(ji,jj+1) ) > 1.E-12 ) .AND. (      &
                  &    MAX(   sshn(ji,jj)             ,  sshn(ji,jj+1) ) >                &
                  &    MAX(  -ht_0(ji,jj)             , -ht_0(ji,jj+1) ) + rn_wdmin1 + rn_wdmin2 )

             IF(ll_tmp1) THEN
               zcpy(ji,jj) = 1.0_wp
             ELSE IF(ll_tmp2) THEN
               ! no worries about  sshn(ji,jj+1) -  sshn(ji,jj  ) = 0, it won't happen ! here
               zcpy(ji,jj) = ABS( (sshn(ji,jj+1) + ht_0(ji,jj+1) - sshn(ji,jj) - ht_0(ji,jj)) &
                           &    / (sshn(ji,jj+1) - sshn(ji,jj  )) )
                zcpy(ji,jj) = max(min( zcpy(ji,jj) , 1.0_wp),0.0_wp)

               ELSE
                  zcpy(ji,jj) = 0._wp
               ENDIF
            END DO
         END DO
         CALL lbc_lnk_multi( 'dynhpg', zcpx, 'U', 1., zcpy, 'V', 1. )
      ENDIF

      ! Clean 3-D work arrays
      zhpi(:,:,:) = 0._wp
      zrhh(:,:,:) = rhd(:,:,:)

      ! Preparing vertical density profile "zrhh(:,:,:)" for hybrid-sco coordinate
      DO jj = 1, jpj
        DO ji = 1, jpi
          jk = mbkt(ji,jj)
          IF(     jk <=  1   ) THEN   ;   zrhh(ji,jj,    :   ) = 0._wp
          ELSEIF( jk ==  2   ) THEN   ;   zrhh(ji,jj,jk+1:jpk) = rhd(ji,jj,jk)
          ELSEIF( jk < jpkm1 ) THEN
             DO jkk = jk+1, jpk
                zrhh(ji,jj,jkk) = interp1(gde3w_n(ji,jj,jkk  ), gde3w_n(ji,jj,jkk-1),   &
                   &                      gde3w_n(ji,jj,jkk-2), zrhh    (ji,jj,jkk-1), zrhh(ji,jj,jkk-2))
             END DO
          ENDIF
        END DO
      END DO

      ! Transfer the depth of "T(:,:,:)" to vertical coordinate "zdept(:,:,:)"
      DO jj = 1, jpj
         DO ji = 1, jpi
            zdept(ji,jj,1) = 0.5_wp * e3w_n(ji,jj,1) - sshn(ji,jj) * znad
         END DO
      END DO

      DO jk = 2, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdept(ji,jj,jk) = zdept(ji,jj,jk-1) + e3w_n(ji,jj,jk)
            END DO
         END DO
      END DO

      fsp(:,:,:) = zrhh (:,:,:)
      xsp(:,:,:) = zdept(:,:,:)

      ! Construct the vertical density profile with the
      ! constrained cubic spline interpolation
      ! rho(z) = asp + bsp*z + csp*z^2 + dsp*z^3
      CALL cspline( fsp, xsp, asp, bsp, csp, dsp, polynomial_type )

      ! Integrate the hydrostatic pressure "zhpi(:,:,:)" at "T(ji,jj,1)"
      DO jj = 2, jpj
        DO ji = 2, jpi
          zrhdt1 = zrhh(ji,jj,1) - interp3( zdept(ji,jj,1), asp(ji,jj,1), bsp(ji,jj,1),  &
             &                                              csp(ji,jj,1), dsp(ji,jj,1) ) * 0.25_wp * e3w_n(ji,jj,1)

          ! assuming linear profile across the top half surface layer
          zhpi(ji,jj,1) =  0.5_wp * e3w_n(ji,jj,1) * zrhdt1
        END DO
      END DO

      ! Calculate the pressure "zhpi(:,:,:)" at "T(ji,jj,2:jpkm1)"
      DO jk = 2, jpkm1
        DO jj = 2, jpj
          DO ji = 2, jpi
            zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) +                                  &
               &             integ_spline( zdept(ji,jj,jk-1), zdept(ji,jj,jk),   &
               &                           asp  (ji,jj,jk-1), bsp  (ji,jj,jk-1), &
               &                           csp  (ji,jj,jk-1), dsp  (ji,jj,jk-1)  )
          END DO
        END DO
      END DO

      ! Z coordinate of U(ji,jj,1:jpkm1) and V(ji,jj,1:jpkm1)

      ! Prepare zsshu_n and zsshv_n
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
!!gm BUG ?    if it is ssh at u- & v-point then it should be:
!          zsshu_n(ji,jj) = (e1e2t(ji,jj) * sshn(ji,jj) + e1e2t(ji+1,jj) * sshn(ji+1,jj)) * &
!                         & r1_e1e2u(ji,jj) * umask(ji,jj,1) * 0.5_wp 
!          zsshv_n(ji,jj) = (e1e2t(ji,jj) * sshn(ji,jj) + e1e2t(ji,jj+1) * sshn(ji,jj+1)) * &
!                         & r1_e1e2v(ji,jj) * vmask(ji,jj,1) * 0.5_wp 
!!gm not this:
          zsshu_n(ji,jj) = (e1e2u(ji,jj) * sshn(ji,jj) + e1e2u(ji+1, jj) * sshn(ji+1,jj)) * &
                         & r1_e1e2u(ji,jj) * umask(ji,jj,1) * 0.5_wp 
          zsshv_n(ji,jj) = (e1e2v(ji,jj) * sshn(ji,jj) + e1e2v(ji+1, jj) * sshn(ji,jj+1)) * &
                         & r1_e1e2v(ji,jj) * vmask(ji,jj,1) * 0.5_wp 
        END DO
      END DO

      CALL lbc_lnk_multi ('dynhpg', zsshu_n, 'U', 1., zsshv_n, 'V', 1. )

      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zu(ji,jj,1) = - ( e3u_n(ji,jj,1) - zsshu_n(ji,jj) * znad) 
          zv(ji,jj,1) = - ( e3v_n(ji,jj,1) - zsshv_n(ji,jj) * znad)
        END DO
      END DO

      DO jk = 2, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zu(ji,jj,jk) = zu(ji,jj,jk-1) - e3u_n(ji,jj,jk)
            zv(ji,jj,jk) = zv(ji,jj,jk-1) - e3v_n(ji,jj,jk)
          END DO
        END DO
      END DO

      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zu(ji,jj,jk) = zu(ji,jj,jk) + 0.5_wp * e3u_n(ji,jj,jk)
            zv(ji,jj,jk) = zv(ji,jj,jk) + 0.5_wp * e3v_n(ji,jj,jk)
          END DO
        END DO
      END DO

      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zu(ji,jj,jk) = MIN(  zu(ji,jj,jk) , MAX( -zdept(ji,jj,jk) , -zdept(ji+1,jj,jk) )  )
            zu(ji,jj,jk) = MAX(  zu(ji,jj,jk) , MIN( -zdept(ji,jj,jk) , -zdept(ji+1,jj,jk) )  )
            zv(ji,jj,jk) = MIN(  zv(ji,jj,jk) , MAX( -zdept(ji,jj,jk) , -zdept(ji,jj+1,jk) )  )
            zv(ji,jj,jk) = MAX(  zv(ji,jj,jk) , MIN( -zdept(ji,jj,jk) , -zdept(ji,jj+1,jk) )  )
          END DO
        END DO
      END DO


      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zpwes = 0._wp; zpwed = 0._wp
            zpnss = 0._wp; zpnsd = 0._wp
            zuijk = zu(ji,jj,jk)
            zvijk = zv(ji,jj,jk)

            !!!!!     for u equation
            IF( jk <= mbku(ji,jj) ) THEN
               IF( -zdept(ji+1,jj,jk) >= -zdept(ji,jj,jk) ) THEN
                 jis = ji + 1; jid = ji
               ELSE
                 jis = ji;     jid = ji +1
               ENDIF

               ! integrate the pressure on the shallow side
               jk1 = jk
               DO WHILE ( -zdept(jis,jj,jk1) > zuijk )
                 IF( jk1 == mbku(ji,jj) ) THEN
                   zuijk = -zdept(jis,jj,jk1)
                   EXIT
                 ENDIF
                 zdeps = MIN(zdept(jis,jj,jk1+1), -zuijk)
                 zpwes = zpwes +                                    &
                      integ_spline(zdept(jis,jj,jk1), zdeps,            &
                             asp(jis,jj,jk1),    bsp(jis,jj,jk1), &
                             csp(jis,jj,jk1),    dsp(jis,jj,jk1))
                 jk1 = jk1 + 1
               END DO

               ! integrate the pressure on the deep side
               jk1 = jk
               DO WHILE ( -zdept(jid,jj,jk1) < zuijk )
                 IF( jk1 == 1 ) THEN
                   zdeps = zdept(jid,jj,1) + MIN(zuijk, sshn(jid,jj)*znad)
                   zrhdt1 = zrhh(jid,jj,1) - interp3(zdept(jid,jj,1), asp(jid,jj,1), &
                                                     bsp(jid,jj,1),   csp(jid,jj,1), &
                                                     dsp(jid,jj,1)) * zdeps
                   zpwed  = zpwed + 0.5_wp * (zrhh(jid,jj,1) + zrhdt1) * zdeps
                   EXIT
                 ENDIF
                 zdeps = MAX(zdept(jid,jj,jk1-1), -zuijk)
                 zpwed = zpwed +                                        &
                        integ_spline(zdeps,              zdept(jid,jj,jk1), &
                               asp(jid,jj,jk1-1), bsp(jid,jj,jk1-1),  &
                               csp(jid,jj,jk1-1), dsp(jid,jj,jk1-1) )
                 jk1 = jk1 - 1
               END DO

               ! update the momentum trends in u direction

               zdpdx1 = zcoef0 * r1_e1u(ji,jj) * ( zhpi(ji+1,jj,jk) - zhpi(ji,jj,jk) )
               IF( .NOT.ln_linssh ) THEN
                 zdpdx2 = zcoef0 * r1_e1u(ji,jj) * &
                    &    ( REAL(jis-jid, wp) * (zpwes + zpwed) + (sshn(ji+1,jj)-sshn(ji,jj)) )
                ELSE
                 zdpdx2 = zcoef0 * r1_e1u(ji,jj) * REAL(jis-jid, wp) * (zpwes + zpwed)
               ENDIF
               IF( ln_wd_il ) THEN
                  zdpdx1 = zdpdx1 * zcpx(ji,jj) * wdrampu(ji,jj)
                  zdpdx2 = zdpdx2 * zcpx(ji,jj) * wdrampu(ji,jj)
               ENDIF
               ua(ji,jj,jk) = ua(ji,jj,jk) + (zdpdx1 + zdpdx2) * umask(ji,jj,jk) 
            ENDIF

            !!!!!     for v equation
            IF( jk <= mbkv(ji,jj) ) THEN
               IF( -zdept(ji,jj+1,jk) >= -zdept(ji,jj,jk) ) THEN
                 jjs = jj + 1; jjd = jj
               ELSE
                 jjs = jj    ; jjd = jj + 1
               ENDIF

               ! integrate the pressure on the shallow side
               jk1 = jk
               DO WHILE ( -zdept(ji,jjs,jk1) > zvijk )
                 IF( jk1 == mbkv(ji,jj) ) THEN
                   zvijk = -zdept(ji,jjs,jk1)
                   EXIT
                 ENDIF
                 zdeps = MIN(zdept(ji,jjs,jk1+1), -zvijk)
                 zpnss = zpnss +                                      &
                        integ_spline(zdept(ji,jjs,jk1), zdeps,            &
                               asp(ji,jjs,jk1),    bsp(ji,jjs,jk1), &
                               csp(ji,jjs,jk1),    dsp(ji,jjs,jk1) )
                 jk1 = jk1 + 1
               END DO

               ! integrate the pressure on the deep side
               jk1 = jk
               DO WHILE ( -zdept(ji,jjd,jk1) < zvijk )
                 IF( jk1 == 1 ) THEN
                   zdeps = zdept(ji,jjd,1) + MIN(zvijk, sshn(ji,jjd)*znad)
                   zrhdt1 = zrhh(ji,jjd,1) - interp3(zdept(ji,jjd,1), asp(ji,jjd,1), &
                                                     bsp(ji,jjd,1),   csp(ji,jjd,1), &
                                                     dsp(ji,jjd,1) ) * zdeps
                   zpnsd  = zpnsd + 0.5_wp * (zrhh(ji,jjd,1) + zrhdt1) * zdeps
                   EXIT
                 ENDIF
                 zdeps = MAX(zdept(ji,jjd,jk1-1), -zvijk)
                 zpnsd = zpnsd +                                        &
                        integ_spline(zdeps,              zdept(ji,jjd,jk1), &
                               asp(ji,jjd,jk1-1), bsp(ji,jjd,jk1-1), &
                               csp(ji,jjd,jk1-1), dsp(ji,jjd,jk1-1) )
                 jk1 = jk1 - 1
               END DO


               ! update the momentum trends in v direction

               zdpdy1 = zcoef0 * r1_e2v(ji,jj) * ( zhpi(ji,jj+1,jk) - zhpi(ji,jj,jk) )
               IF( .NOT.ln_linssh ) THEN
                  zdpdy2 = zcoef0 * r1_e2v(ji,jj) * &
                          ( REAL(jjs-jjd, wp) * (zpnss + zpnsd) + (sshn(ji,jj+1)-sshn(ji,jj)) )
               ELSE
                  zdpdy2 = zcoef0 * r1_e2v(ji,jj) * REAL(jjs-jjd, wp) * (zpnss + zpnsd )
               ENDIF
               IF( ln_wd_il ) THEN
                  zdpdy1 = zdpdy1 * zcpy(ji,jj) * wdrampv(ji,jj) 
                  zdpdy2 = zdpdy2 * zcpy(ji,jj) * wdrampv(ji,jj) 
               ENDIF

               va(ji,jj,jk) = va(ji,jj,jk) + (zdpdy1 + zdpdy2) * vmask(ji,jj,jk)
            ENDIF
               !
            END DO
         END DO
      END DO
      !
      IF( ln_wd_il )   DEALLOCATE( zcpx, zcpy )
      !
   END SUBROUTINE hpg_prj


   SUBROUTINE cspline( fsp, xsp, asp, bsp, csp, dsp, polynomial_type )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE cspline  ***
      !!
      !! ** Purpose :   constrained cubic spline interpolation
      !!
      !! ** Method  :   f(x) = asp + bsp*x + csp*x^2 + dsp*x^3
      !!
      !! Reference: CJC Kruger, Constrained Cubic Spline Interpoltation
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   fsp, xsp           ! value and coordinate
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   asp, bsp, csp, dsp ! coefficients of the interpoated function
      INTEGER                   , INTENT(in   ) ::   polynomial_type    ! 1: cubic spline   ;   2: Linear
      !
      INTEGER  ::   ji, jj, jk                 ! dummy loop indices
      INTEGER  ::   jpi, jpj, jpkm1
      REAL(wp) ::   zdf1, zdf2, zddf1, zddf2, ztmp1, ztmp2, zdxtmp
      REAL(wp) ::   zdxtmp1, zdxtmp2, zalpha
      REAL(wp) ::   zdf(size(fsp,3))
      !!----------------------------------------------------------------------
      !
!!gm  WHAT !!!!!   THIS IS VERY DANGEROUS !!!!!  
      jpi   = size(fsp,1)
      jpj   = size(fsp,2)
      jpkm1 = MAX( 1, size(fsp,3) - 1 )
      !
      IF (polynomial_type == 1) THEN     ! Constrained Cubic Spline
         DO ji = 1, jpi
            DO jj = 1, jpj
           !!Fritsch&Butland's method, 1984 (preferred, but more computation)
           !    DO jk = 2, jpkm1-1
           !       zdxtmp1 = xsp(ji,jj,jk)   - xsp(ji,jj,jk-1)
           !       zdxtmp2 = xsp(ji,jj,jk+1) - xsp(ji,jj,jk)
           !       zdf1    = ( fsp(ji,jj,jk)   - fsp(ji,jj,jk-1) ) / zdxtmp1
           !       zdf2    = ( fsp(ji,jj,jk+1) - fsp(ji,jj,jk)   ) / zdxtmp2
           !
           !       zalpha = ( zdxtmp1 + 2._wp * zdxtmp2 ) / ( zdxtmp1 + zdxtmp2 ) / 3._wp
           !
           !       IF(zdf1 * zdf2 <= 0._wp) THEN
           !           zdf(jk) = 0._wp
           !       ELSE
           !         zdf(jk) = zdf1 * zdf2 / ( ( 1._wp - zalpha ) * zdf1 + zalpha * zdf2 )
           !       ENDIF
           !    END DO

           !!Simply geometric average
               DO jk = 2, jpkm1-1
                  zdf1 = (fsp(ji,jj,jk  ) - fsp(ji,jj,jk-1)) / (xsp(ji,jj,jk  ) - xsp(ji,jj,jk-1))
                  zdf2 = (fsp(ji,jj,jk+1) - fsp(ji,jj,jk  )) / (xsp(ji,jj,jk+1) - xsp(ji,jj,jk  ))

                  IF(zdf1 * zdf2 <= 0._wp) THEN
                     zdf(jk) = 0._wp
                  ELSE
                     zdf(jk) = 2._wp * zdf1 * zdf2 / (zdf1 + zdf2)
                  ENDIF
               END DO

               zdf(1)     = 1.5_wp * ( fsp(ji,jj,2) - fsp(ji,jj,1) ) / &
                          &          ( xsp(ji,jj,2) - xsp(ji,jj,1) )           -  0.5_wp * zdf(2)
               zdf(jpkm1) = 1.5_wp * ( fsp(ji,jj,jpkm1) - fsp(ji,jj,jpkm1-1) ) / &
                          &          ( xsp(ji,jj,jpkm1) - xsp(ji,jj,jpkm1-1) ) - 0.5_wp * zdf(jpkm1 - 1)

               DO jk = 1, jpkm1 - 1
                 zdxtmp = xsp(ji,jj,jk+1) - xsp(ji,jj,jk)
                 ztmp1  = (zdf(jk+1) + 2._wp * zdf(jk)) / zdxtmp
                 ztmp2  =  6._wp * (fsp(ji,jj,jk+1) - fsp(ji,jj,jk)) / zdxtmp / zdxtmp
                 zddf1  = -2._wp * ztmp1 + ztmp2
                 ztmp1  = (2._wp * zdf(jk+1) + zdf(jk)) / zdxtmp
                 zddf2  =  2._wp * ztmp1 - ztmp2

                 dsp(ji,jj,jk) = (zddf2 - zddf1) / 6._wp / zdxtmp
                 csp(ji,jj,jk) = ( xsp(ji,jj,jk+1) * zddf1 - xsp(ji,jj,jk)*zddf2 ) / 2._wp / zdxtmp
                 bsp(ji,jj,jk) = ( fsp(ji,jj,jk+1) - fsp(ji,jj,jk) ) / zdxtmp - &
                               & csp(ji,jj,jk) * ( xsp(ji,jj,jk+1) + xsp(ji,jj,jk) ) - &
                               & dsp(ji,jj,jk) * ((xsp(ji,jj,jk+1) + xsp(ji,jj,jk))**2 - &
                               &                   xsp(ji,jj,jk+1) * xsp(ji,jj,jk))
                 asp(ji,jj,jk) = fsp(ji,jj,jk) - xsp(ji,jj,jk) * (bsp(ji,jj,jk) + &
                               &                (xsp(ji,jj,jk) * (csp(ji,jj,jk) + &
                               &                 dsp(ji,jj,jk) * xsp(ji,jj,jk))))
               END DO
            END DO
         END DO

      ELSEIF ( polynomial_type == 2 ) THEN     ! Linear
         DO ji = 1, jpi
            DO jj = 1, jpj
               DO jk = 1, jpkm1-1
                  zdxtmp =xsp(ji,jj,jk+1) - xsp(ji,jj,jk)
                  ztmp1 = fsp(ji,jj,jk+1) - fsp(ji,jj,jk)

                  dsp(ji,jj,jk) = 0._wp
                  csp(ji,jj,jk) = 0._wp
                  bsp(ji,jj,jk) = ztmp1 / zdxtmp
                  asp(ji,jj,jk) = fsp(ji,jj,jk) - bsp(ji,jj,jk) * xsp(ji,jj,jk)
               END DO
            END DO
         END DO
         !
      ELSE
         CALL ctl_stop( 'invalid polynomial type in cspline' )
      ENDIF
      !
   END SUBROUTINE cspline


   FUNCTION interp1(x, xl, xr, fl, fr)  RESULT(f)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE interp1  ***
      !!
      !! ** Purpose :   1-d linear interpolation
      !!
      !! ** Method  :   interpolation is straight forward
      !!                extrapolation is also permitted (no value limit)
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::  x, xl, xr, fl, fr
      REAL(wp)             ::  f ! result of the interpolation (extrapolation)
      REAL(wp)             ::  zdeltx
      !!----------------------------------------------------------------------
      !
      zdeltx = xr - xl
      IF( abs(zdeltx) <= 10._wp * EPSILON(x) ) THEN
         f = 0.5_wp * (fl + fr)
      ELSE
         f = ( (x - xl ) * fr - ( x - xr ) * fl ) / zdeltx
      ENDIF
      !
   END FUNCTION interp1


   FUNCTION interp2( x, a, b, c, d )  RESULT(f)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE interp1  ***
      !!
      !! ** Purpose :   1-d constrained cubic spline interpolation
      !!
      !! ** Method  :  cubic spline interpolation
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::  x, a, b, c, d
      REAL(wp)             ::  f ! value from the interpolation
      !!----------------------------------------------------------------------
      !
      f = a + x* ( b + x * ( c + d * x ) )
      !
   END FUNCTION interp2


   FUNCTION interp3( x, a, b, c, d )  RESULT(f)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE interp1  ***
      !!
      !! ** Purpose :   Calculate the first order of derivative of
      !!                a cubic spline function y=a+b*x+c*x^2+d*x^3
      !!
      !! ** Method  :   f=dy/dx=b+2*c*x+3*d*x^2
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::  x, a, b, c, d
      REAL(wp)             ::  f ! value from the interpolation
      !!----------------------------------------------------------------------
      !
      f = b + x * ( 2._wp * c + 3._wp * d * x)
      !
   END FUNCTION interp3


   FUNCTION integ_spline( xl, xr, a, b, c, d )  RESULT(f)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE interp1  ***
      !!
      !! ** Purpose :   1-d constrained cubic spline integration
      !!
      !! ** Method  :  integrate polynomial a+bx+cx^2+dx^3 from xl to xr
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::  xl, xr, a, b, c, d
      REAL(wp)             ::  za1, za2, za3
      REAL(wp)             ::  f                   ! integration result
      !!----------------------------------------------------------------------
      !
      za1 = 0.5_wp * b
      za2 = c / 3.0_wp
      za3 = 0.25_wp * d
      !
      f  = xr * ( a + xr * ( za1 + xr * ( za2 + za3 * xr ) ) ) - &
         & xl * ( a + xl * ( za1 + xl * ( za2 + za3 * xl ) ) )
      !
   END FUNCTION integ_spline

   !!======================================================================
END MODULE dynhpg

