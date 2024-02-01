MODULE traadv
   !!==============================================================================
   !!                       ***  MODULE  traadv  ***
   !! Ocean active tracers:  advection trend 
   !!==============================================================================
   !! History :  2.0  !  2005-11  (G. Madec)  Original code
   !!            3.3  !  2010-09  (C. Ethe, G. Madec)  merge TRC-TRA + switch from velocity to transport
   !!            3.6  !  2011-06  (G. Madec)  Addition of Mixed Layer Eddy parameterisation
   !!            3.7  !  2014-05  (G. Madec)  Add 2nd/4th order cases for CEN and FCT schemes 
   !!             -   !  2014-12  (G. Madec) suppression of cross land advection option
   !!            3.6  !  2015-06  (E. Clementi) Addition of Stokes drift in case of wave coupling
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv       : compute ocean tracer advection trend
   !!   tra_adv_init  : control the different options of advection scheme
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE dom_oce        ! ocean space and time domain
   USE domvvl         ! variable vertical scale factors
   USE sbcwave        ! wave module
   USE sbc_oce        ! surface boundary condition: ocean
   USE traadv_cen     ! centered scheme            (tra_adv_cen  routine)
   USE traadv_fct     ! FCT      scheme            (tra_adv_fct  routine)
   USE traadv_mus     ! MUSCL    scheme            (tra_adv_mus  routine)
   USE traadv_ubs     ! UBS      scheme            (tra_adv_ubs  routine)
   USE traadv_qck     ! QUICKEST scheme            (tra_adv_qck  routine)
   USE tramle         ! Mixed Layer Eddy transport (tra_mle_trp  routine)
   USE ldftra         ! Eddy Induced transport     (ldf_eiv_trp  routine)
   USE ldfslp         ! Lateral diffusion: slopes of neutral surfaces
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers 
   USE diaptr         ! Poleward heat transport 
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O module
   USE prtctl         ! Print control
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv        ! called by step.F90
   PUBLIC   tra_adv_init   ! called by nemogcm.F90

   !                            !!* Namelist namtra_adv *
   LOGICAL ::   ln_traadv_OFF    ! no advection on T and S
   LOGICAL ::   ln_traadv_cen    ! centered scheme flag
   INTEGER ::      nn_cen_h, nn_cen_v   ! =2/4 : horizontal and vertical choices of the order of CEN scheme
   LOGICAL ::   ln_traadv_fct    ! FCT scheme flag
   INTEGER ::      nn_fct_h, nn_fct_v   ! =2/4 : horizontal and vertical choices of the order of FCT scheme
   LOGICAL ::   ln_traadv_mus    ! MUSCL scheme flag
   LOGICAL ::      ln_mus_ups           ! use upstream scheme in vivcinity of river mouths
   LOGICAL ::   ln_traadv_ubs    ! UBS scheme flag
   INTEGER ::      nn_ubs_v             ! =2/4 : vertical choice of the order of UBS scheme
   LOGICAL ::   ln_traadv_qck    ! QUICKEST scheme flag

   INTEGER ::   nadv             ! choice of the type of advection scheme
   !                             ! associated indices:
   INTEGER, PARAMETER ::   np_NO_adv  = 0   ! no T-S advection
   INTEGER, PARAMETER ::   np_CEN     = 1   ! 2nd/4th order centered scheme
   INTEGER, PARAMETER ::   np_FCT     = 2   ! 2nd/4th order Flux Corrected Transport scheme
   INTEGER, PARAMETER ::   np_MUS     = 3   ! MUSCL scheme
   INTEGER, PARAMETER ::   np_UBS     = 4   ! 3rd order Upstream Biased Scheme
   INTEGER, PARAMETER ::   np_QCK     = 5   ! QUICK scheme
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv  ***
      !!
      !! ** Purpose :   compute the ocean tracer advection trend.
      !!
      !! ** Method  : - Update (ua,va) with the advection term following nadv
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER ::   jk   ! dummy loop index
      REAL(wp), DIMENSION(jpi,jpj,jpk)        :: zun, zvn, zwn   ! 3D workspace
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdt, ztrds
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tra_adv')
      !
      !                                          ! set time step
      IF( neuler == 0 .AND. kt == nit000 ) THEN   ;   r2dt =         rdt   ! at nit000             (Euler)
      ELSEIF( kt <= nit000 + 1 )           THEN   ;   r2dt = 2._wp * rdt   ! at nit000 or nit000+1 (Leapfrog)
      ENDIF
      !
      !                                         !==  effective transport  ==!
      zun(:,:,jpk) = 0._wp
      zvn(:,:,jpk) = 0._wp
      zwn(:,:,jpk) = 0._wp
      IF( ln_wave .AND. ln_sdw )  THEN
         DO jk = 1, jpkm1                                                       ! eulerian transport + Stokes Drift
            zun(:,:,jk) = e2u  (:,:) * e3u_n(:,:,jk) * ( un(:,:,jk) + usd(:,:,jk) )
            zvn(:,:,jk) = e1v  (:,:) * e3v_n(:,:,jk) * ( vn(:,:,jk) + vsd(:,:,jk) )
            zwn(:,:,jk) = e1e2t(:,:)                 * ( wn(:,:,jk) + wsd(:,:,jk) )
         END DO
      ELSE
         DO jk = 1, jpkm1
            zun(:,:,jk) = e2u  (:,:) * e3u_n(:,:,jk) * un(:,:,jk)               ! eulerian transport only
            zvn(:,:,jk) = e1v  (:,:) * e3v_n(:,:,jk) * vn(:,:,jk)
            zwn(:,:,jk) = e1e2t(:,:)                 * wn(:,:,jk)
         END DO
      ENDIF
      !
      IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN                                ! add z-tilde and/or vvl corrections
         zun(:,:,:) = zun(:,:,:) + un_td(:,:,:)
         zvn(:,:,:) = zvn(:,:,:) + vn_td(:,:,:)
      ENDIF
      !
      zun(:,:,jpk) = 0._wp                                                      ! no transport trough the bottom
      zvn(:,:,jpk) = 0._wp
      zwn(:,:,jpk) = 0._wp
      !
      IF( ln_ldfeiv .AND. .NOT. ln_traldf_triad )   &
         &              CALL ldf_eiv_trp( kt, nit000, zun, zvn, zwn, 'TRA' )   ! add the eiv transport (if necessary)
      !
      IF( ln_mle    )   CALL tra_mle_trp( kt, nit000, zun, zvn, zwn, 'TRA' )   ! add the mle transport (if necessary)
      !
      CALL iom_put( "uocetr_eff", zun )                                        ! output effective transport      
      CALL iom_put( "vocetr_eff", zvn )
      CALL iom_put( "wocetr_eff", zwn )
      !
!!gm ???
      IF( ln_diaptr )   CALL dia_ptr( zvn )                                    ! diagnose the effective MSF 
!!gm ???
      !
      IF( l_trdtra )   THEN                    !* Save ta and sa trends
         ALLOCATE( ztrdt(jpi,jpj,jpk), ztrds(jpi,jpj,jpk) )
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal)
      ENDIF
      !
      SELECT CASE ( nadv )                      !==  compute advection trend and add it to general trend  ==!
      !
      CASE ( np_CEN )                                 ! Centered scheme : 2nd / 4th order
         CALL tra_adv_cen    ( kt, nit000, 'TRA',         zun, zvn, zwn     , tsn, tsa, jpts, nn_cen_h, nn_cen_v )
      CASE ( np_FCT )                                 ! FCT scheme      : 2nd / 4th order
         CALL tra_adv_fct    ( kt, nit000, 'TRA', r2dt, zun, zvn, zwn, tsb, tsn, tsa, jpts, nn_fct_h, nn_fct_v )
      CASE ( np_MUS )                                 ! MUSCL
         CALL tra_adv_mus    ( kt, nit000, 'TRA', r2dt, zun, zvn, zwn, tsb,      tsa, jpts        , ln_mus_ups ) 
      CASE ( np_UBS )                                 ! UBS
         CALL tra_adv_ubs    ( kt, nit000, 'TRA', r2dt, zun, zvn, zwn, tsb, tsn, tsa, jpts        , nn_ubs_v   )
      CASE ( np_QCK )                                 ! QUICKEST
         CALL tra_adv_qck    ( kt, nit000, 'TRA', r2dt, zun, zvn, zwn, tsb, tsn, tsa, jpts                     )
      !
      END SELECT
      !
      IF( l_trdtra )   THEN                      ! save the advective trends for further diagnostics
         DO jk = 1, jpkm1
            ztrdt(:,:,jk) = tsa(:,:,jk,jp_tem) - ztrdt(:,:,jk)
            ztrds(:,:,jk) = tsa(:,:,jk,jp_sal) - ztrds(:,:,jk)
         END DO
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_totad, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_totad, ztrds )
         DEALLOCATE( ztrdt, ztrds )
      ENDIF
      !                                              ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' adv  - Ta: ', mask1=tmask,               &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( ln_timing )   CALL timing_stop( 'tra_adv' )
      !
   END SUBROUTINE tra_adv


   SUBROUTINE tra_adv_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_init  ***
      !!                
      !! ** Purpose :   Control the consistency between namelist options for 
      !!              tracer advection schemes and set nadv
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ios   ! Local integers
      !
      NAMELIST/namtra_adv/ ln_traadv_OFF,                        &   ! No advection
         &                 ln_traadv_cen , nn_cen_h, nn_cen_v,   &   ! CEN
         &                 ln_traadv_fct , nn_fct_h, nn_fct_v,   &   ! FCT
         &                 ln_traadv_mus , ln_mus_ups,           &   ! MUSCL
         &                 ln_traadv_ubs ,           nn_ubs_v,   &   ! UBS
         &                 ln_traadv_qck                             ! QCK
      !!----------------------------------------------------------------------
      !
      !                                !==  Namelist  ==!
      REWIND( numnam_ref )                   ! Namelist namtra_adv in reference namelist : Tracer advection scheme
      READ  ( numnam_ref, namtra_adv, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtra_adv in reference namelist' )
      !
      REWIND( numnam_cfg )                   ! Namelist namtra_adv in configuration namelist : Tracer advection scheme
      READ  ( numnam_cfg, namtra_adv, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtra_adv in configuration namelist' )
      IF(lwm) WRITE( numond, namtra_adv )
      !
      IF(lwp) THEN                           ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_adv_init : choice/control of the tracer advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_adv : chose a advection scheme for tracers'
         WRITE(numout,*) '      No advection on T & S                     ln_traadv_OFF = ', ln_traadv_OFF
         WRITE(numout,*) '      centered scheme                           ln_traadv_cen = ', ln_traadv_cen
         WRITE(numout,*) '            horizontal 2nd/4th order               nn_cen_h   = ', nn_fct_h
         WRITE(numout,*) '            vertical   2nd/4th order               nn_cen_v   = ', nn_fct_v
         WRITE(numout,*) '      Flux Corrected Transport scheme           ln_traadv_fct = ', ln_traadv_fct
         WRITE(numout,*) '            horizontal 2nd/4th order               nn_fct_h   = ', nn_fct_h
         WRITE(numout,*) '            vertical   2nd/4th order               nn_fct_v   = ', nn_fct_v
         WRITE(numout,*) '      MUSCL scheme                              ln_traadv_mus = ', ln_traadv_mus
         WRITE(numout,*) '            + upstream scheme near river mouths    ln_mus_ups = ', ln_mus_ups
         WRITE(numout,*) '      UBS scheme                                ln_traadv_ubs = ', ln_traadv_ubs
         WRITE(numout,*) '            vertical   2nd/4th order               nn_ubs_v   = ', nn_ubs_v
         WRITE(numout,*) '      QUICKEST scheme                           ln_traadv_qck = ', ln_traadv_qck
      ENDIF
      !
      !                                !==  Parameter control & set nadv ==!
      ioptio = 0                       
      IF( ln_traadv_OFF ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_NO_adv   ;   ENDIF
      IF( ln_traadv_cen ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_CEN      ;   ENDIF
      IF( ln_traadv_fct ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_FCT      ;   ENDIF
      IF( ln_traadv_mus ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_MUS      ;   ENDIF
      IF( ln_traadv_ubs ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_UBS      ;   ENDIF
      IF( ln_traadv_qck ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_QCK      ;   ENDIF
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'tra_adv_init: Choose ONE advection option in namelist namtra_adv' )
      !
      IF( ln_traadv_cen .AND. ( nn_cen_h /= 2 .AND. nn_cen_h /= 4 )   &          ! Centered
                        .AND. ( nn_cen_v /= 2 .AND. nn_cen_v /= 4 )   ) THEN
        CALL ctl_stop( 'tra_adv_init: CEN scheme, choose 2nd or 4th order' )
      ENDIF
      IF( ln_traadv_fct .AND. ( nn_fct_h /= 2 .AND. nn_fct_h /= 4 )   &          ! FCT
                        .AND. ( nn_fct_v /= 2 .AND. nn_fct_v /= 4 )   ) THEN
        CALL ctl_stop( 'tra_adv_init: FCT scheme, choose 2nd or 4th order' )
      ENDIF
      IF( ln_traadv_ubs .AND. ( nn_ubs_v /= 2 .AND. nn_ubs_v /= 4 )   ) THEN     ! UBS
        CALL ctl_stop( 'tra_adv_init: UBS scheme, choose 2nd or 4th order' )
      ENDIF
      IF( ln_traadv_ubs .AND. nn_ubs_v == 4 ) THEN
         CALL ctl_warn( 'tra_adv_init: UBS scheme, only 2nd FCT scheme available on the vertical. It will be used' )
      ENDIF
      IF( ln_isfcav ) THEN                                                       ! ice-shelf cavities
         IF(  ln_traadv_cen .AND. nn_cen_v == 4    .OR.   &                            ! NO 4th order with ISF
            & ln_traadv_fct .AND. nn_fct_v == 4   )   CALL ctl_stop( 'tra_adv_init: 4th order COMPACT scheme not allowed with ISF' )
      ENDIF
      !
      !                                !==  Print the choice  ==!  
      IF(lwp) THEN
         WRITE(numout,*)
         SELECT CASE ( nadv )
         CASE( np_NO_adv  )   ;   WRITE(numout,*) '   ==>>>   NO T-S advection'
         CASE( np_CEN     )   ;   WRITE(numout,*) '   ==>>>   CEN      scheme is used. Horizontal order: ', nn_cen_h,   &
            &                                                                        ' Vertical   order: ', nn_cen_v
         CASE( np_FCT     )   ;   WRITE(numout,*) '   ==>>>   FCT      scheme is used. Horizontal order: ', nn_fct_h,   &
            &                                                                        ' Vertical   order: ', nn_fct_v
         CASE( np_MUS     )   ;   WRITE(numout,*) '   ==>>>   MUSCL    scheme is used'
         CASE( np_UBS     )   ;   WRITE(numout,*) '   ==>>>   UBS      scheme is used'
         CASE( np_QCK     )   ;   WRITE(numout,*) '   ==>>>   QUICKEST scheme is used'
         END SELECT
      ENDIF
      !
      CALL tra_mle_init            !== initialisation of the Mixed Layer Eddy parametrisation (MLE)  ==!
      !
   END SUBROUTINE tra_adv_init

  !!======================================================================
END MODULE traadv
