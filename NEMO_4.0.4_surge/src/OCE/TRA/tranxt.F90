MODULE tranxt
   !!======================================================================
   !!                       ***  MODULE  tranxt  ***
   !! Ocean active tracers:  time stepping on temperature and salinity
   !!======================================================================
   !! History :  OPA  !  1991-11  (G. Madec)  Original code
   !!            7.0  !  1993-03  (M. Guyon)  symetrical conditions
   !!            8.0  !  1996-02  (G. Madec & M. Imbard)  opa release 8.0
   !!             -   !  1996-04  (A. Weaver)  Euler forward step
   !!            8.2  !  1999-02  (G. Madec, N. Grima)  semi-implicit pressure grad.
   !!  NEMO      1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2002-11  (C. Talandier, A-M Treguier) Open boundaries
   !!             -   !  2005-04  (C. Deltel) Add Asselin trend in the ML budget
   !!            2.0  !  2006-02  (L. Debreu, C. Mazauric) Agrif implementation
   !!            3.0  !  2008-06  (G. Madec)  time stepping always done in trazdf
   !!            3.1  !  2009-02  (G. Madec, R. Benshila)  re-introduce the vvl option
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  semi-implicit hpg with asselin filter + modified LF-RA
   !!             -   !  2010-05  (C. Ethe, G. Madec)  merge TRC-TRA
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_nxt       : time stepping on tracers
   !!   tra_nxt_fix   : time stepping on tracers : fixed    volume case
   !!   tra_nxt_vvl   : time stepping on tracers : variable volume case
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE sbc_oce         ! surface boundary condition: ocean
   USE sbcrnf          ! river runoffs
   USE sbcisf          ! ice shelf melting
   USE zdf_oce         ! ocean vertical mixing
   USE domvvl          ! variable volume
   USE trd_oce         ! trends: ocean variables
   USE trdtra          ! trends manager: tracers 
   USE traqsr          ! penetrative solar radiation (needed for nksr)
   USE phycst          ! physical constant
   USE ldftra          ! lateral physics : tracers
   USE ldfslp          ! lateral physics : slopes
   USE bdy_oce  , ONLY : ln_bdy
   USE bdytra          ! open boundary condition (bdy_tra routine)
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE timing          ! Timing
#if defined key_agrif
   USE agrif_oce_interp
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_nxt       ! routine called by step.F90
   PUBLIC   tra_nxt_fix   ! to be used in trcnxt
   PUBLIC   tra_nxt_vvl   ! to be used in trcnxt

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_nxt( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tranxt  ***
      !!
      !! ** Purpose :   Apply the boundary condition on the after temperature  
      !!             and salinity fields, achieved the time stepping by adding
      !!             the Asselin filter on now fields and swapping the fields.
      !! 
      !! ** Method  :   At this stage of the computation, ta and sa are the 
      !!             after temperature and salinity as the time stepping has
      !!             been performed in trazdf_imp or trazdf_exp module.
      !!
      !!              - Apply lateral boundary conditions on (ta,sa) 
      !!             at the local domain   boundaries through lbc_lnk call, 
      !!             at the one-way open boundaries (ln_bdy=T), 
      !!             at the AGRIF zoom   boundaries (lk_agrif=T)
      !!
      !!              - Update lateral boundary conditions on AGRIF children
      !!             domains (lk_agrif=T)
      !!
      !! ** Action  : - tsb & tsn ready for the next time step
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::   zfact            ! local scalars
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrdt, ztrds
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start( 'tra_nxt')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt : achieve the time stepping by Asselin filter and array swap'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

      ! Update after tracer on domain lateral boundaries
      ! 
#if defined key_agrif
      CALL Agrif_tra                     ! AGRIF zoom boundaries
#endif
      !                                              ! local domain boundaries  (T-point, unchanged sign)
      CALL lbc_lnk_multi( 'tranxt', tsa(:,:,:,jp_tem), 'T', 1., tsa(:,:,:,jp_sal), 'T', 1. )
      !
      IF( ln_bdy )   CALL bdy_tra( kt )  ! BDY open boundaries
 
      ! set time step size (Euler/Leapfrog)
      IF( neuler == 0 .AND. kt == nit000 ) THEN   ;   r2dt =        rdt   ! at nit000             (Euler)
      ELSEIF( kt <= nit000 + 1 )           THEN   ;   r2dt = 2._wp* rdt   ! at nit000 or nit000+1 (Leapfrog)
      ENDIF

      ! trends computation initialisation
      IF( l_trdtra )   THEN                    
         ALLOCATE( ztrdt(jpi,jpj,jpk) , ztrds(jpi,jpj,jpk) )
         ztrdt(:,:,jpk) = 0._wp
         ztrds(:,:,jpk) = 0._wp
         IF( ln_traldf_iso ) THEN              ! diagnose the "pure" Kz diffusive trend 
            CALL trd_tra( kt, 'TRA', jp_tem, jptra_zdfp, ztrdt )
            CALL trd_tra( kt, 'TRA', jp_sal, jptra_zdfp, ztrds )
         ENDIF
         ! total trend for the non-time-filtered variables. 
         zfact = 1.0 / rdt
         ! G Nurser 23 Mar 2017. Recalculate trend as Delta(e3t*T)/e3tn; e3tn cancel from tsn terms
         DO jk = 1, jpkm1
            ztrdt(:,:,jk) = ( tsa(:,:,jk,jp_tem)*e3t_a(:,:,jk) / e3t_n(:,:,jk) - tsn(:,:,jk,jp_tem)) * zfact
            ztrds(:,:,jk) = ( tsa(:,:,jk,jp_sal)*e3t_a(:,:,jk) / e3t_n(:,:,jk) - tsn(:,:,jk,jp_sal)) * zfact
         END DO
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_tot, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_tot, ztrds )
         IF( ln_linssh ) THEN       ! linear sea surface height only
            ! Store now fields before applying the Asselin filter 
            ! in order to calculate Asselin filter trend later.
            ztrdt(:,:,:) = tsn(:,:,:,jp_tem) 
            ztrds(:,:,:) = tsn(:,:,:,jp_sal)
         ENDIF
      ENDIF

      IF( neuler == 0 .AND. kt == nit000 ) THEN       ! Euler time-stepping at first time-step (only swap)
         DO jn = 1, jpts
            DO jk = 1, jpkm1
               tsn(:,:,jk,jn) = tsa(:,:,jk,jn)    
            END DO
         END DO
         IF (l_trdtra .AND. .NOT. ln_linssh ) THEN   ! Zero Asselin filter contribution must be explicitly written out since for vvl
            !                                        ! Asselin filter is output by tra_nxt_vvl that is not called on this time step
            ztrdt(:,:,:) = 0._wp
            ztrds(:,:,:) = 0._wp
            CALL trd_tra( kt, 'TRA', jp_tem, jptra_atf, ztrdt )
            CALL trd_tra( kt, 'TRA', jp_sal, jptra_atf, ztrds )
         END IF
         !
      ELSE                                            ! Leap-Frog + Asselin filter time stepping
         !
         IF( ln_linssh ) THEN   ;   CALL tra_nxt_fix( kt, nit000,      'TRA', tsb, tsn, tsa, jpts )  ! linear free surface 
         ELSE                   ;   CALL tra_nxt_vvl( kt, nit000, rdt, 'TRA', tsb, tsn, tsa,   &
           &                                                                sbc_tsc, sbc_tsc_b, jpts )  ! non-linear free surface
         ENDIF
         !
         CALL lbc_lnk_multi( 'tranxt', tsb(:,:,:,jp_tem), 'T', 1., tsb(:,:,:,jp_sal), 'T', 1., &
                  &          tsn(:,:,:,jp_tem), 'T', 1., tsn(:,:,:,jp_sal), 'T', 1., &
                  &          tsa(:,:,:,jp_tem), 'T', 1., tsa(:,:,:,jp_sal), 'T', 1.  )
         !
      ENDIF     
      !
      IF( l_trdtra .AND. ln_linssh ) THEN      ! trend of the Asselin filter (tb filtered - tb)/dt     
         zfact = 1._wp / r2dt             
         DO jk = 1, jpkm1
            ztrdt(:,:,jk) = ( tsb(:,:,jk,jp_tem) - ztrdt(:,:,jk) ) * zfact
            ztrds(:,:,jk) = ( tsb(:,:,jk,jp_sal) - ztrds(:,:,jk) ) * zfact
         END DO
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_atf, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_atf, ztrds )
      END IF
      IF( l_trdtra )   DEALLOCATE( ztrdt , ztrds )
      !
      !                        ! control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsn(:,:,:,jp_tem), clinfo1=' nxt  - Tn: ', mask1=tmask,   &
         &                       tab3d_2=tsn(:,:,:,jp_sal), clinfo2=       ' Sn: ', mask2=tmask )
      !
      IF( ln_timing )   CALL timing_stop('tra_nxt')
      !
   END SUBROUTINE tra_nxt


   SUBROUTINE tra_nxt_fix( kt, kit000, cdtype, ptb, ptn, pta, kjpt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_nxt_fix  ***
      !!
      !! ** Purpose :   fixed volume: apply the Asselin time filter and 
      !!                swap the tracer fields.
      !! 
      !! ** Method  : - Apply a Asselin time filter on now fields.
      !!              - swap tracer fields to prepare the next time_step.
      !!
      !! ** Action  : - tsb & tsn ready for the next time step
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::  kt        ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::  kit000    ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::  cdtype    ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::  kjpt      ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::  ptb       ! before tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::  ptn       ! now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::  pta       ! tracer trend
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::   ztn, ztd         ! local scalars
      !!----------------------------------------------------------------------
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt_fix : time stepping', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      DO jn = 1, kjpt
         !
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  ztn = ptn(ji,jj,jk,jn)                                    
                  ztd = pta(ji,jj,jk,jn) - 2._wp * ztn + ptb(ji,jj,jk,jn)  ! time laplacian on tracers
                  !
                  ptb(ji,jj,jk,jn) = ztn + atfp * ztd                      ! ptb <-- filtered ptn 
                  ptn(ji,jj,jk,jn) = pta(ji,jj,jk,jn)                      ! ptn <-- pta
               END DO
           END DO
         END DO
         !
      END DO
      !
   END SUBROUTINE tra_nxt_fix


   SUBROUTINE tra_nxt_vvl( kt, kit000, p2dt, cdtype, ptb, ptn, pta, psbc_tc, psbc_tc_b, kjpt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_nxt_vvl  ***
      !!
      !! ** Purpose :   Time varying volume: apply the Asselin time filter  
      !!                and swap the tracer fields.
      !! 
      !! ** Method  : - Apply a thickness weighted Asselin time filter on now fields.
      !!              - swap tracer fields to prepare the next time_step.
      !!             tb  = ( e3t_n*tn + atfp*[ e3t_b*tb - 2 e3t_n*tn + e3t_a*ta ] )
      !!                  /( e3t_n    + atfp*[ e3t_b    - 2 e3t_n    + e3t_a    ] )
      !!             tn  = ta 
      !!
      !! ** Action  : - tsb & tsn ready for the next time step
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::  kt        ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::  kit000    ! first time step index
      REAL(wp)                             , INTENT(in   ) ::  p2dt      ! time-step
      CHARACTER(len=3)                     , INTENT(in   ) ::  cdtype    ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::  kjpt      ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::  ptb       ! before tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::  ptn       ! now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::  pta       ! tracer trend
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(in   ) ::  psbc_tc   ! surface tracer content
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(in   ) ::  psbc_tc_b ! before surface tracer content
      !
      LOGICAL  ::   ll_traqsr, ll_rnf, ll_isf   ! local logical
      INTEGER  ::   ji, jj, jk, jn              ! dummy loop indices
      REAL(wp) ::   zfact, zfact1, ztc_a , ztc_n , ztc_b , ztc_f , ztc_d    ! local scalar
      REAL(wp) ::   zfact2, ze3t_b, ze3t_n, ze3t_a, ze3t_f, ze3t_d, zscale  !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   ztrd_atf
      !!----------------------------------------------------------------------
      !
      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt_vvl : time stepping', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      IF( cdtype == 'TRA' )  THEN   
         ll_traqsr  = ln_traqsr        ! active  tracers case  and  solar penetration
         ll_rnf     = ln_rnf           ! active  tracers case  and  river runoffs
         ll_isf     = ln_isf           ! active  tracers case  and  ice shelf melting
      ELSE                          ! passive tracers case
         ll_traqsr  = .FALSE.          ! NO solar penetration
         ll_rnf     = .FALSE.          ! NO river runoffs ????          !!gm BUG ?  
         ll_isf     = .FALSE.          ! NO ice shelf melting/freezing  !!gm BUG ?? 
      ENDIF
      !
      IF( ( l_trdtra .AND. cdtype == 'TRA' ) .OR. ( l_trdtrc .AND. cdtype == 'TRC' ) )   THEN
         ALLOCATE( ztrd_atf(jpi,jpj,jpk,kjpt) )
         ztrd_atf(:,:,:,:) = 0.0_wp
      ENDIF
      zfact = 1._wp / p2dt
      zfact1 = atfp * p2dt
      zfact2 = zfact1 * r1_rau0
      DO jn = 1, kjpt      
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  ze3t_b = e3t_b(ji,jj,jk)
                  ze3t_n = e3t_n(ji,jj,jk)
                  ze3t_a = e3t_a(ji,jj,jk)
                  !                                         ! tracer content at Before, now and after
                  ztc_b  = ptb(ji,jj,jk,jn) * ze3t_b
                  ztc_n  = ptn(ji,jj,jk,jn) * ze3t_n
                  ztc_a  = pta(ji,jj,jk,jn) * ze3t_a
                  !
                  ze3t_d = ze3t_a - 2. * ze3t_n + ze3t_b
                  ztc_d  = ztc_a  - 2. * ztc_n  + ztc_b
                  !
                  ze3t_f = ze3t_n + atfp * ze3t_d
                  ztc_f  = ztc_n  + atfp * ztc_d
                  !
                  zscale = zfact2 * e3t_n(ji,jj,jk) / ( ht_n(ji,jj) + 1._wp - ssmask(ji,jj) )
                  ze3t_f = ze3t_f - zscale * ( emp_b(ji,jj) - emp(ji,jj) )
                  IF ( ll_rnf ) ze3t_f = ze3t_f + zscale * (    rnf_b(ji,jj) -    rnf(ji,jj) )
                  IF ( ll_isf ) ze3t_f = ze3t_f - zscale * ( fwfisf_b(ji,jj) - fwfisf(ji,jj) )

                  IF( jk == mikt(ji,jj) ) THEN           ! first level 
                     ztc_f  = ztc_f  - zfact1 * ( psbc_tc(ji,jj,jn) - psbc_tc_b(ji,jj,jn) )
                  ENDIF
                  !
                  ! solar penetration (temperature only)
                  IF( ll_traqsr .AND. jn == jp_tem .AND. jk <= nksr )                            & 
                     &     ztc_f  = ztc_f  - zfact1 * ( qsr_hc(ji,jj,jk) - qsr_hc_b(ji,jj,jk) ) 
                     !
                  ! river runoff
                  IF( ll_rnf .AND. jk <= nk_rnf(ji,jj) )                                          &
                     &     ztc_f  = ztc_f  - zfact1 * ( rnf_tsc(ji,jj,jn) - rnf_tsc_b(ji,jj,jn) ) & 
                     &                              * e3t_n(ji,jj,jk) / h_rnf(ji,jj)
                     !
                  ! ice shelf
                  IF( ll_isf ) THEN
                     ! level fully include in the Losch_2008 ice shelf boundary layer
                     IF ( jk >= misfkt(ji,jj) .AND. jk < misfkb(ji,jj) )                          &
                        ztc_f  = ztc_f  - zfact1 * ( risf_tsc(ji,jj,jn) - risf_tsc_b(ji,jj,jn) )  &
                               &                 * e3t_n(ji,jj,jk) * r1_hisf_tbl (ji,jj)
                     ! level partially include in Losch_2008 ice shelf boundary layer 
                     IF ( jk == misfkb(ji,jj) )                                                   &
                        ztc_f  = ztc_f  - zfact1 * ( risf_tsc(ji,jj,jn) - risf_tsc_b(ji,jj,jn) )  &
                               &                 * e3t_n(ji,jj,jk) * r1_hisf_tbl (ji,jj) * ralpha(ji,jj)
                  END IF
                  !
                  ze3t_f = 1.e0 / ze3t_f
                  ptb(ji,jj,jk,jn) = ztc_f * ze3t_f       ! ptb <-- ptn filtered
                  ptn(ji,jj,jk,jn) = pta(ji,jj,jk,jn)     ! ptn <-- pta
                  !
                  IF( ( l_trdtra .and. cdtype == 'TRA' ) .OR. ( l_trdtrc .and. cdtype == 'TRC' ) ) THEN
                     ztrd_atf(ji,jj,jk,jn) = (ztc_f - ztc_n) * zfact/ze3t_n
                  ENDIF
                  !
               END DO
            END DO
         END DO
         ! 
      END DO
      !
      IF( ( l_trdtra .AND. cdtype == 'TRA' ) .OR. ( l_trdtrc .AND. cdtype == 'TRC' ) )   THEN
         IF( l_trdtra .AND. cdtype == 'TRA' ) THEN 
            CALL trd_tra( kt, cdtype, jp_tem, jptra_atf, ztrd_atf(:,:,:,jp_tem) )
            CALL trd_tra( kt, cdtype, jp_sal, jptra_atf, ztrd_atf(:,:,:,jp_sal) )
         ENDIF
         IF( l_trdtrc .AND. cdtype == 'TRC' ) THEN
            DO jn = 1, kjpt
               CALL trd_tra( kt, cdtype, jn, jptra_atf, ztrd_atf(:,:,:,jn) )
            END DO
         ENDIF
         DEALLOCATE( ztrd_atf )
      ENDIF
      !
   END SUBROUTINE tra_nxt_vvl

   !!======================================================================
END MODULE tranxt
