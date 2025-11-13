MODULE dynnxt
   !!=========================================================================
   !!                       ***  MODULE  dynnxt  ***
   !! Ocean dynamics: time stepping
   !!=========================================================================
   !! History :  OPA  !  1987-02  (P. Andrich, D. L Hostis)  Original code
   !!                 !  1990-10  (C. Levy, G. Madec)
   !!            7.0  !  1993-03  (M. Guyon)  symetrical conditions
   !!            8.0  !  1997-02  (G. Madec & M. Imbard)  opa, release 8.0
   !!            8.2  !  1997-04  (A. Weaver)  Euler forward step
   !!             -   !  1997-06  (G. Madec)  lateral boudary cond., lbc routine
   !!    NEMO    1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2002-10  (C. Talandier, A-M. Treguier) Open boundary cond.
   !!            2.0  !  2005-11  (V. Garnier) Surface pressure gradient organization
   !!            2.3  !  2007-07  (D. Storkey) Calls to BDY routines. 
   !!            3.2  !  2009-06  (G. Madec, R.Benshila)  re-introduce the vvl option
   !!            3.3  !  2010-09  (D. Storkey, E.O'Dea) Bug fix for BDY module
   !!            3.3  !  2011-03  (P. Oddo) Bug fix for time-splitting+(BDY-OBC) and not VVL
   !!            3.5  !  2013-07  (J. Chanut) Compliant with time splitting changes
   !!            3.6  !  2014-04  (G. Madec) add the diagnostic of the time filter trends
   !!            3.7  !  2015-11  (J. Chanut) Free surface simplification
   !!-------------------------------------------------------------------------
  
   !!-------------------------------------------------------------------------
   !!   dyn_nxt       : obtain the next (after) horizontal velocity
   !!-------------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbcrnf         ! river runoffs
   USE sbcisf         ! ice shelf
   USE phycst         ! physical constants
   USE dynadv         ! dynamics: vector invariant versus flux form
   USE dynspg_ts      ! surface pressure gradient: split-explicit scheme
   USE domvvl         ! variable volume
   USE bdy_oce , ONLY : ln_bdy
   USE bdydta         ! ocean open boundary conditions
   USE bdydyn         ! ocean open boundary conditions
   USE bdyvol         ! ocean open boundary condition (bdy_vol routines)
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   USE trdken         ! trend manager: kinetic energy
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lbclnk         ! lateral boundary condition (or mpp link)
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control
   USE timing         ! Timing
   USE zdfdrg ,  ONLY : ln_drgice_imp, rCdU_top
#if defined key_agrif
   USE agrif_oce_interp
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC    dyn_nxt   ! routine called by step.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_nxt ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_nxt  ***
      !!                   
      !! ** Purpose :   Finalize after horizontal velocity. Apply the boundary 
      !!             condition on the after velocity, achieve the time stepping 
      !!             by applying the Asselin filter on now fields and swapping 
      !!             the fields.
      !!
      !! ** Method  : * Ensure after velocities transport matches time splitting
      !!             estimate (ln_dynspg_ts=T)
      !!
      !!              * Apply lateral boundary conditions on after velocity 
      !!             at the local domain boundaries through lbc_lnk call,
      !!             at the one-way open boundaries (ln_bdy=T),
      !!             at the AGRIF zoom   boundaries (lk_agrif=T)
      !!
      !!              * Apply the time filter applied and swap of the dynamics
      !!             arrays to start the next time step:
      !!                (ub,vb) = (un,vn) + atfp [ (ub,vb) + (ua,va) - 2 (un,vn) ]
      !!                (un,vn) = (ua,va).
      !!             Note that with flux form advection and non linear free surface,
      !!             the time filter is applied on thickness weighted velocity.
      !!             As a result, dyn_nxt MUST be called after tra_nxt.
      !!
      !! ** Action :   ub,vb   filtered before horizontal velocity of next time-step
      !!               un,vn   now horizontal velocity of next time-step
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikt          ! local integers
      REAL(wp) ::   zue3a, zue3n, zue3b, zuf, zcoef    ! local scalars
      REAL(wp) ::   zve3a, zve3n, zve3b, zvf, z1_2dt   !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   zue, zve
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   zutau, zvtau
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ze3u_f, ze3v_f, zua, zva 
      !!----------------------------------------------------------------------
      !
      IF( ln_timing    )   CALL timing_start('dyn_nxt')
      IF( ln_dynspg_ts )   ALLOCATE( zue(jpi,jpj)     , zve(jpi,jpj)     )
      IF( l_trddyn     )   ALLOCATE( zua(jpi,jpj,jpk) , zva(jpi,jpj,jpk) )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_nxt : time stepping'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

      IF ( ln_dynspg_ts ) THEN
         ! Ensure below that barotropic velocities match time splitting estimate
         ! Compute actual transport and replace it with ts estimate at "after" time step
         zue(:,:) = e3u_a(:,:,1) * ua(:,:,1) * umask(:,:,1)
         zve(:,:) = e3v_a(:,:,1) * va(:,:,1) * vmask(:,:,1)
         DO jk = 2, jpkm1
            zue(:,:) = zue(:,:) + e3u_a(:,:,jk) * ua(:,:,jk) * umask(:,:,jk)
            zve(:,:) = zve(:,:) + e3v_a(:,:,jk) * va(:,:,jk) * vmask(:,:,jk)
         END DO
         DO jk = 1, jpkm1
            ua(:,:,jk) = ( ua(:,:,jk) - zue(:,:) * r1_hu_a(:,:) + ua_b(:,:) ) * umask(:,:,jk)
            va(:,:,jk) = ( va(:,:,jk) - zve(:,:) * r1_hv_a(:,:) + va_b(:,:) ) * vmask(:,:,jk)
         END DO
         !
         IF( .NOT.ln_bt_fw ) THEN
            ! Remove advective velocity from "now velocities" 
            ! prior to asselin filtering     
            ! In the forward case, this is done below after asselin filtering   
            ! so that asselin contribution is removed at the same time 
            DO jk = 1, jpkm1
               un(:,:,jk) = ( un(:,:,jk) - un_adv(:,:)*r1_hu_n(:,:) + un_b(:,:) )*umask(:,:,jk)
               vn(:,:,jk) = ( vn(:,:,jk) - vn_adv(:,:)*r1_hv_n(:,:) + vn_b(:,:) )*vmask(:,:,jk)
            END DO  
         ENDIF
      ENDIF

      ! Update after velocity on domain lateral boundaries
      ! --------------------------------------------------      
# if defined key_agrif
      CALL Agrif_dyn( kt )             !* AGRIF zoom boundaries
# endif
      !
      CALL lbc_lnk_multi( 'dynnxt', ua, 'U', -1., va, 'V', -1. )     !* local domain boundaries
      !
      !                                !* BDY open boundaries
      IF( ln_bdy .AND. ln_dynspg_exp )   CALL bdy_dyn( kt )
      IF( ln_bdy .AND. ln_dynspg_ts  )   CALL bdy_dyn( kt, dyn3d_only=.true. )

!!$   Do we need a call to bdy_vol here??
      !
      IF( l_trddyn ) THEN             ! prepare the atf trend computation + some diagnostics
         z1_2dt = 1._wp / (2. * rdt)        ! Euler or leap-frog time step 
         IF( neuler == 0 .AND. kt == nit000 )   z1_2dt = 1._wp / rdt
         !
         !                                  ! Kinetic energy and Conversion
         IF( ln_KE_trd  )   CALL trd_dyn( ua, va, jpdyn_ken, kt )
         !
         IF( ln_dyn_trd ) THEN              ! 3D output: total momentum trends
            zua(:,:,:) = ( ua(:,:,:) - ub(:,:,:) ) * z1_2dt
            zva(:,:,:) = ( va(:,:,:) - vb(:,:,:) ) * z1_2dt
            CALL iom_put( "utrd_tot", zua )        ! total momentum trends, except the asselin time filter
            CALL iom_put( "vtrd_tot", zva )
         ENDIF
         !
         zua(:,:,:) = un(:,:,:)             ! save the now velocity before the asselin filter
         zva(:,:,:) = vn(:,:,:)             ! (caution: there will be a shift by 1 timestep in the
         !                                  !  computation of the asselin filter trends)
      ENDIF

      ! Time filter and swap of dynamics arrays
      ! ------------------------------------------
      IF( neuler == 0 .AND. kt == nit000 ) THEN        !* Euler at first time-step: only swap
         DO jk = 1, jpkm1
            ub(:,:,jk) = un(:,:,jk)                         ! ub <-- un
            vb(:,:,jk) = vn(:,:,jk)
            un(:,:,jk) = ua(:,:,jk)                         ! un <-- ua
            vn(:,:,jk) = va(:,:,jk)
         END DO
         IF( .NOT.ln_linssh ) THEN                          ! e3._b <-- e3._n
!!gm BUG ????    I don't understand why it is not : e3._n <-- e3._a  
            DO jk = 1, jpkm1
!               e3t_b(:,:,jk) = e3t_n(:,:,jk)
!               e3u_b(:,:,jk) = e3u_n(:,:,jk)
!               e3v_b(:,:,jk) = e3v_n(:,:,jk)
               !
               e3t_n(:,:,jk) = e3t_a(:,:,jk)
               e3u_n(:,:,jk) = e3u_a(:,:,jk)
               e3v_n(:,:,jk) = e3v_a(:,:,jk)
            END DO
!!gm BUG end
         ENDIF
                                                            ! 
         
      ELSE                                             !* Leap-Frog : Asselin filter and swap
         !                                ! =============!
         IF( ln_linssh ) THEN             ! Fixed volume !
            !                             ! =============!
            DO jk = 1, jpkm1                              
               DO jj = 1, jpj
                  DO ji = 1, jpi    
                     zuf = un(ji,jj,jk) + atfp * ( ub(ji,jj,jk) - 2._wp * un(ji,jj,jk) + ua(ji,jj,jk) )
                     zvf = vn(ji,jj,jk) + atfp * ( vb(ji,jj,jk) - 2._wp * vn(ji,jj,jk) + va(ji,jj,jk) )
                     !
                     ub(ji,jj,jk) = zuf                      ! ub <-- filtered velocity
                     vb(ji,jj,jk) = zvf
                     un(ji,jj,jk) = ua(ji,jj,jk)             ! un <-- ua
                     vn(ji,jj,jk) = va(ji,jj,jk)
                  END DO
               END DO
            END DO
            !                             ! ================!
         ELSE                             ! Variable volume !
            !                             ! ================!
            ! Before scale factor at t-points
            ! (used as a now filtered scale factor until the swap)
            ! ----------------------------------------------------
            DO jk = 1, jpkm1
               e3t_b(:,:,jk) = e3t_n(:,:,jk) + atfp * ( e3t_b(:,:,jk) - 2._wp * e3t_n(:,:,jk) + e3t_a(:,:,jk) )
            END DO
            ! Add volume filter correction: compatibility with tracer advection scheme
            ! => time filter + conservation correction (only at the first level)
            zcoef = atfp * rdt * r1_rau0

            DO jk = 1, jpkm1
               e3t_b(:,:,jk) = e3t_b(:,:,jk) - zcoef * ( emp_b(:,:) - emp(:,:) ) * tmask(:,:,jk) & 
                             &                       * e3t_n(:,:,jk) /  ( ht_n(:,:) + 1._wp - ssmask(:,:) )
            END DO

            IF ( ln_rnf ) THEN
               DO jk = 1, jpkm1
                  e3t_b(:,:,jk) = e3t_b(:,:,jk) + zcoef * ( rnf_b(:,:) - rnf(:,:) ) * tmask(:,:,jk) & 
                                &                       * e3t_n(:,:,jk) /  ( ht_n(:,:) + 1._wp - ssmask(:,:) )
               END DO
            ENDIF

            IF ( ln_isf ) THEN
               DO jk = 1, jpkm1
                  e3t_b(:,:,jk) = e3t_b(:,:,jk) - zcoef * ( fwfisf_b(:,:) - fwfisf(:,:) ) * tmask(:,:,jk) & 
                                &                       * e3t_n(:,:,jk) /  ( ht_n(:,:) + 1._wp - ssmask(:,:) )
               END DO
            ENDIF
            !
            IF( ln_dynadv_vec ) THEN      ! Asselin filter applied on velocity
               ! Before filtered scale factor at (u/v)-points
               CALL dom_vvl_interpol( e3t_b(:,:,:), e3u_b(:,:,:), 'U' )
               CALL dom_vvl_interpol( e3t_b(:,:,:), e3v_b(:,:,:), 'V' )
               DO jk = 1, jpkm1
                  DO jj = 1, jpj
                     DO ji = 1, jpi
                        zuf = un(ji,jj,jk) + atfp * ( ub(ji,jj,jk) - 2._wp * un(ji,jj,jk) + ua(ji,jj,jk) )
                        zvf = vn(ji,jj,jk) + atfp * ( vb(ji,jj,jk) - 2._wp * vn(ji,jj,jk) + va(ji,jj,jk) )
                        !
                        ub(ji,jj,jk) = zuf                      ! ub <-- filtered velocity
                        vb(ji,jj,jk) = zvf
                        un(ji,jj,jk) = ua(ji,jj,jk)             ! un <-- ua
                        vn(ji,jj,jk) = va(ji,jj,jk)
                     END DO
                  END DO
               END DO
               !
            ELSE                          ! Asselin filter applied on thickness weighted velocity
               !
               ALLOCATE( ze3u_f(jpi,jpj,jpk) , ze3v_f(jpi,jpj,jpk) )
               ! Before filtered scale factor at (u/v)-points stored in ze3u_f, ze3v_f
               CALL dom_vvl_interpol( e3t_b(:,:,:), ze3u_f, 'U' )
               CALL dom_vvl_interpol( e3t_b(:,:,:), ze3v_f, 'V' )
               DO jk = 1, jpkm1
                  DO jj = 1, jpj
                     DO ji = 1, jpi                  
                        zue3a = e3u_a(ji,jj,jk) * ua(ji,jj,jk)
                        zve3a = e3v_a(ji,jj,jk) * va(ji,jj,jk)
                        zue3n = e3u_n(ji,jj,jk) * un(ji,jj,jk)
                        zve3n = e3v_n(ji,jj,jk) * vn(ji,jj,jk)
                        zue3b = e3u_b(ji,jj,jk) * ub(ji,jj,jk)
                        zve3b = e3v_b(ji,jj,jk) * vb(ji,jj,jk)
                        !
                        zuf = ( zue3n + atfp * ( zue3b - 2._wp * zue3n  + zue3a ) ) / ze3u_f(ji,jj,jk)
                        zvf = ( zve3n + atfp * ( zve3b - 2._wp * zve3n  + zve3a ) ) / ze3v_f(ji,jj,jk)
                        !
                        ub(ji,jj,jk) = zuf                     ! ub <-- filtered velocity
                        vb(ji,jj,jk) = zvf
                        un(ji,jj,jk) = ua(ji,jj,jk)            ! un <-- ua
                        vn(ji,jj,jk) = va(ji,jj,jk)
                     END DO
                  END DO
               END DO
               e3u_b(:,:,1:jpkm1) = ze3u_f(:,:,1:jpkm1)        ! e3u_b <-- filtered scale factor
               e3v_b(:,:,1:jpkm1) = ze3v_f(:,:,1:jpkm1)
               !
               DEALLOCATE( ze3u_f , ze3v_f )
            ENDIF
            !
         ENDIF
         !
         IF( ln_dynspg_ts .AND. ln_bt_fw ) THEN
            ! Revert "before" velocities to time split estimate
            ! Doing it here also means that asselin filter contribution is removed  
            zue(:,:) = e3u_b(:,:,1) * ub(:,:,1) * umask(:,:,1)
            zve(:,:) = e3v_b(:,:,1) * vb(:,:,1) * vmask(:,:,1)    
            DO jk = 2, jpkm1
               zue(:,:) = zue(:,:) + e3u_b(:,:,jk) * ub(:,:,jk) * umask(:,:,jk)
               zve(:,:) = zve(:,:) + e3v_b(:,:,jk) * vb(:,:,jk) * vmask(:,:,jk)    
            END DO
            DO jk = 1, jpkm1
               ub(:,:,jk) = ub(:,:,jk) - (zue(:,:) * r1_hu_n(:,:) - un_b(:,:)) * umask(:,:,jk)
               vb(:,:,jk) = vb(:,:,jk) - (zve(:,:) * r1_hv_n(:,:) - vn_b(:,:)) * vmask(:,:,jk)
            END DO
         ENDIF
         !
      ENDIF ! neuler =/0
      !
      ! Set "now" and "before" barotropic velocities for next time step:
      ! JC: Would be more clever to swap variables than to make a full vertical
      ! integration
      !
      !
      IF(.NOT.ln_linssh ) THEN
         hu_b(:,:) = e3u_b(:,:,1) * umask(:,:,1)
         hv_b(:,:) = e3v_b(:,:,1) * vmask(:,:,1)
         DO jk = 2, jpkm1
            hu_b(:,:) = hu_b(:,:) + e3u_b(:,:,jk) * umask(:,:,jk)
            hv_b(:,:) = hv_b(:,:) + e3v_b(:,:,jk) * vmask(:,:,jk)
         END DO
         r1_hu_b(:,:) = ssumask(:,:) / ( hu_b(:,:) + 1._wp - ssumask(:,:) )
         r1_hv_b(:,:) = ssvmask(:,:) / ( hv_b(:,:) + 1._wp - ssvmask(:,:) )
      ENDIF
      !
      un_b(:,:) = e3u_a(:,:,1) * un(:,:,1) * umask(:,:,1)
      ub_b(:,:) = e3u_b(:,:,1) * ub(:,:,1) * umask(:,:,1)
      vn_b(:,:) = e3v_a(:,:,1) * vn(:,:,1) * vmask(:,:,1)
      vb_b(:,:) = e3v_b(:,:,1) * vb(:,:,1) * vmask(:,:,1)
      DO jk = 2, jpkm1
         un_b(:,:) = un_b(:,:) + e3u_a(:,:,jk) * un(:,:,jk) * umask(:,:,jk)
         ub_b(:,:) = ub_b(:,:) + e3u_b(:,:,jk) * ub(:,:,jk) * umask(:,:,jk)
         vn_b(:,:) = vn_b(:,:) + e3v_a(:,:,jk) * vn(:,:,jk) * vmask(:,:,jk)
         vb_b(:,:) = vb_b(:,:) + e3v_b(:,:,jk) * vb(:,:,jk) * vmask(:,:,jk)
      END DO
      un_b(:,:) = un_b(:,:) * r1_hu_a(:,:)
      vn_b(:,:) = vn_b(:,:) * r1_hv_a(:,:)
      ub_b(:,:) = ub_b(:,:) * r1_hu_b(:,:)
      vb_b(:,:) = vb_b(:,:) * r1_hv_b(:,:)
      !
      IF( .NOT.ln_dynspg_ts ) THEN        ! output the barotropic currents
         CALL iom_put(  "ubar", un_b(:,:) )
         CALL iom_put(  "vbar", vn_b(:,:) )
      ENDIF
      IF( l_trddyn ) THEN                ! 3D output: asselin filter trends on momentum
         zua(:,:,:) = ( ub(:,:,:) - zua(:,:,:) ) * z1_2dt
         zva(:,:,:) = ( vb(:,:,:) - zva(:,:,:) ) * z1_2dt
         CALL trd_dyn( zua, zva, jpdyn_atf, kt )
      ENDIF
      !
      IF ( iom_use("utau") ) THEN
         IF ( ln_drgice_imp.OR.ln_isfcav ) THEN
            ALLOCATE(zutau(jpi,jpj)) 
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  jk = miku(ji,jj) 
                  zutau(ji,jj) = utau(ji,jj) & 
                  &  + 0.5_wp * rau0 * (rCdU_top(ji+1,jj)+rCdU_top(ji,jj)) * ua(ji,jj,jk) 
               END DO
            END DO
            CALL lbc_lnk( 'dynnxt' , zutau, 'U', -1.)
            CALL iom_put(  "utau", zutau(:,:) )
            DEALLOCATE(zutau)
         ELSE
            CALL iom_put(  "utau", utau(:,:) )
         ENDIF
      ENDIF
      !
      IF ( iom_use("vtau") ) THEN
         IF ( ln_drgice_imp.OR.ln_isfcav ) THEN
            ALLOCATE(zvtau(jpi,jpj))
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  jk = mikv(ji,jj)
                  zvtau(ji,jj) = vtau(ji,jj) &
                  &  + 0.5_wp * rau0 * (rCdU_top(ji,jj+1)+rCdU_top(ji,jj)) * va(ji,jj,jk)
               END DO
            END DO
            CALL lbc_lnk( 'dynnxt' , zvtau, 'V', -1.)
            CALL iom_put(  "vtau", zvtau(:,:) )
            DEALLOCATE(zvtau)
         ELSE
            CALL iom_put(  "vtau", vtau(:,:) )
         ENDIF
      ENDIF
      !
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=un, clinfo1=' nxt  - Un: ', mask1=umask,   &
         &                       tab3d_2=vn, clinfo2=' Vn: '       , mask2=vmask )
      ! 
      IF( ln_dynspg_ts )   DEALLOCATE( zue, zve )
      IF( l_trddyn     )   DEALLOCATE( zua, zva )
      IF( ln_timing    )   CALL timing_stop('dyn_nxt')
      !
   END SUBROUTINE dyn_nxt

   !!=========================================================================
END MODULE dynnxt
