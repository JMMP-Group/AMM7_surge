MODULE traldf_lap_blp
   !!==============================================================================
   !!                       ***  MODULE  traldf_lap_blp  ***
   !! Ocean tracers:  lateral diffusivity trend  (laplacian and bilaplacian)
   !!==============================================================================
   !! History :  3.7  ! 2014-01  (G. Madec, S. Masson)  Original code, re-entrant laplacian 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_ldf_lap   : tracer trend update with iso-level laplacian diffusive operator
   !!   tra_ldf_blp   : tracer trend update with iso-level or iso-neutral bilaplacian operator
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE dom_oce        ! ocean space and time domain
   USE ldftra         ! lateral physics: eddy diffusivity
   USE traldf_iso     ! iso-neutral lateral diffusion (standard operator)     (tra_ldf_iso   routine)
   USE traldf_triad   ! iso-neutral lateral diffusion (triad    operator)     (tra_ldf_triad routine)
   USE diaptr         ! poleward transport diagnostics
   USE diaar5         ! AR5 diagnostics
   USE trc_oce        ! share passive tracers/Ocean variables
   USE zpshde         ! partial step: hor. derivative     (zps_hde routine)
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distribued memory computing library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf_lap   ! called by traldf.F90
   PUBLIC   tra_ldf_blp   ! called by traldf.F90

   LOGICAL  ::   l_ptr   ! flag to compute poleward transport
   LOGICAL  ::   l_hst   ! flag to compute heat transport

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_ldf_lap( kt, kit000, cdtype, pahu, pahv, pgu , pgv ,   &
      &                                                    pgui, pgvi,   &
      &                                        ptb , pta , kjpt, kpass ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_lap  ***
      !!                   
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   Second order diffusive operator evaluated using before
      !!      fields (forward time scheme). The horizontal diffusive trends of 
      !!      the tracer is given by:
      !!          difft = 1/(e1e2t*e3t) {  di-1[ pahu e2u*e3u/e1u di(tb) ]
      !!                                 + dj-1[ pahv e1v*e3v/e2v dj(tb) ] }
      !!      Add this trend to the general tracer trend pta :
      !!          pta = pta + difft
      !!
      !! ** Action  : - Update pta arrays with the before iso-level 
      !!                harmonic mixing trend.
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000     ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      INTEGER                              , INTENT(in   ) ::   kpass      ! =1/2 first or second passage
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pahu, pahv ! eddy diffusivity at u- and v-points  [m2/s]
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(in   ) ::   pgu, pgv   ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,    kjpt), INTENT(in   ) ::   pgui, pgvi ! tracer gradient at top   levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta        ! tracer trend 
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::   zsign            ! local scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ztu, ztv, zaheeu, zaheev
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp )  THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_ldf_lap : iso-level laplacian diffusion on ', cdtype, ', pass=', kpass
         WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF
      !
      l_hst = .FALSE.
      l_ptr = .FALSE.
      IF( cdtype == 'TRA' .AND. ln_diaptr )                                                l_ptr = .TRUE. 
      IF( cdtype == 'TRA' .AND. ( iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. &
         &                        iom_use("uadv_salttr") .OR. iom_use("vadv_salttr")  ) )  l_hst = .TRUE.
      !
      !                                !==  Initialization of metric arrays used for all tracers  ==!
      IF( kpass == 1 ) THEN   ;   zsign =  1._wp      ! bilaplacian operator require a minus sign (eddy diffusivity >0)
      ELSE                    ;   zsign = -1._wp
      ENDIF
      DO jk = 1, jpkm1
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zaheeu(ji,jj,jk) = zsign * pahu(ji,jj,jk) * e2_e1u(ji,jj) * e3u_n(ji,jj,jk)   !!gm   * umask(ji,jj,jk) pah masked!
               zaheev(ji,jj,jk) = zsign * pahv(ji,jj,jk) * e1_e2v(ji,jj) * e3v_n(ji,jj,jk)   !!gm   * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      !                             ! =========== !
      DO jn = 1, kjpt               ! tracer loop !
         !                          ! =========== !    
         !                               
         DO jk = 1, jpkm1              !== First derivative (gradient)  ==!
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  ztu(ji,jj,jk) = zaheeu(ji,jj,jk) * ( ptb(ji+1,jj  ,jk,jn) - ptb(ji,jj,jk,jn) )
                  ztv(ji,jj,jk) = zaheev(ji,jj,jk) * ( ptb(ji  ,jj+1,jk,jn) - ptb(ji,jj,jk,jn) )
               END DO
            END DO
         END DO  
         IF( ln_zps ) THEN                ! set gradient at bottom/top ocean level
            DO jj = 1, jpjm1                    ! bottom
               DO ji = 1, fs_jpim1
                  ztu(ji,jj,mbku(ji,jj)) = zaheeu(ji,jj,mbku(ji,jj)) * pgu(ji,jj,jn)
                  ztv(ji,jj,mbkv(ji,jj)) = zaheev(ji,jj,mbkv(ji,jj)) * pgv(ji,jj,jn)
               END DO
            END DO  
            IF( ln_isfcav ) THEN                ! top in ocean cavities only
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.
                     IF( miku(ji,jj) > 1 )   ztu(ji,jj,miku(ji,jj)) = zaheeu(ji,jj,miku(ji,jj)) * pgui(ji,jj,jn) 
                     IF( mikv(ji,jj) > 1 )   ztv(ji,jj,mikv(ji,jj)) = zaheev(ji,jj,mikv(ji,jj)) * pgvi(ji,jj,jn) 
                  END DO
               END DO
            ENDIF
         ENDIF
         !
         DO jk = 1, jpkm1              !== Second derivative (divergence) added to the general tracer trends  ==!
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk)     &
                     &                                   + ztv(ji,jj,jk) - ztv(ji,jj-1,jk) )   &
                     &                                / ( e1e2t(ji,jj) * e3t_n(ji,jj,jk) )
               END DO
            END DO
         END DO  
         !
         !                             !== "Poleward" diffusive heat or salt transports  ==!
         IF( ( kpass == 1 .AND. .NOT.ln_traldf_blp ) .OR.  &     !==  first pass only (  laplacian)  ==!
             ( kpass == 2 .AND.      ln_traldf_blp ) ) THEN      !==  2nd   pass only (bilaplacian)  ==!

            IF( l_ptr )  CALL dia_ptr_hst( jn, 'ldf', -ztv(:,:,:)  )
            IF( l_hst )  CALL dia_ar5_hst( jn, 'ldf', -ztu(:,:,:), -ztv(:,:,:) )
         ENDIF
         !                          ! ==================
      END DO                        ! end of tracer loop
      !                             ! ==================
      !
   END SUBROUTINE tra_ldf_lap
   

   SUBROUTINE tra_ldf_blp( kt, kit000, cdtype, pahu, pahv, pgu , pgv ,   &
      &                                                    pgui, pgvi,   &
      &                                                    ptb , pta , kjpt, kldf )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tra_ldf_blp  ***
      !!                    
      !! ** Purpose :   Compute the before lateral tracer diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   The lateral diffusive trends is provided by a bilaplacian
      !!      operator applied to before field (forward in time).
      !!      It is computed by two successive calls to laplacian routine
      !!
      !! ** Action :   pta   updated with the before rotated bilaplacian diffusion
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000     ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      INTEGER                              , INTENT(in   ) ::   kldf       ! type of operator used
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pahu, pahv ! eddy diffusivity at u- and v-points  [m2/s]
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(in   ) ::   pgu, pgv   ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,    kjpt), INTENT(in   ) ::   pgui, pgvi ! tracer gradient at top levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta        ! tracer trend
      !
      INTEGER ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt) :: zlap         ! laplacian at t-point
      REAL(wp), DIMENSION(jpi,jpj,    kjpt) :: zglu, zglv   ! bottom GRADh of the laplacian (u- and v-points)
      REAL(wp), DIMENSION(jpi,jpj,    kjpt) :: zgui, zgvi   ! top    GRADh of the laplacian (u- and v-points)
      !!---------------------------------------------------------------------
      !
      IF( kt == kit000 .AND. lwp )  THEN
         WRITE(numout,*)
         SELECT CASE ( kldf )
         CASE ( np_blp    )   ;   WRITE(numout,*) 'tra_ldf_blp : iso-level   bilaplacian operator on ', cdtype
         CASE ( np_blp_i  )   ;   WRITE(numout,*) 'tra_ldf_blp : iso-neutral bilaplacian operator on ', cdtype, ' (Standard)'
         CASE ( np_blp_it )   ;   WRITE(numout,*) 'tra_ldf_blp : iso-neutral bilaplacian operator on ', cdtype, ' (triad)'
         END SELECT
         WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      zlap(:,:,:,:) = 0._wp
      !
      SELECT CASE ( kldf )       !==  1st laplacian applied to ptb (output in zlap)  ==!
      !
      CASE ( np_blp    )               ! iso-level bilaplacian
         CALL tra_ldf_lap  ( kt, kit000, cdtype, pahu, pahv, pgu, pgv, pgui, pgvi, ptb,      zlap, kjpt, 1 )
      CASE ( np_blp_i  )               ! rotated   bilaplacian : standard operator (Madec)
         CALL tra_ldf_iso  ( kt, kit000, cdtype, pahu, pahv, pgu, pgv, pgui, pgvi, ptb, ptb, zlap, kjpt, 1 )
      CASE ( np_blp_it )               ! rotated  bilaplacian : triad operator (griffies)
         CALL tra_ldf_triad( kt, kit000, cdtype, pahu, pahv, pgu, pgv, pgui, pgvi, ptb, ptb, zlap, kjpt, 1 )
      END SELECT
      !
      CALL lbc_lnk( 'traldf_lap_blp', zlap(:,:,:,:) , 'T', 1. )     ! Lateral boundary conditions (unchanged sign)
      !                                               ! Partial top/bottom cell: GRADh( zlap )  
      IF( ln_isfcav .AND. ln_zps ) THEN   ;   CALL zps_hde_isf( kt, kjpt, zlap, zglu, zglv, zgui, zgvi )  ! both top & bottom
      ELSEIF(             ln_zps ) THEN   ;   CALL zps_hde    ( kt, kjpt, zlap, zglu, zglv )              ! only bottom 
      ENDIF
      !
      SELECT CASE ( kldf )       !==  2nd laplacian applied to zlap (output in pta)  ==!
      !
      CASE ( np_blp    )               ! iso-level bilaplacian
         CALL tra_ldf_lap  ( kt, kit000, cdtype, pahu, pahv, zglu, zglv, zgui, zgvi, zlap, pta,      kjpt, 2 )
      CASE ( np_blp_i  )               ! rotated   bilaplacian : standard operator (Madec)
         CALL tra_ldf_iso  ( kt, kit000, cdtype, pahu, pahv, zglu, zglv, zgui, zgvi, zlap, ptb, pta, kjpt, 2 )
      CASE ( np_blp_it )               ! rotated  bilaplacian : triad operator (griffies)
         CALL tra_ldf_triad( kt, kit000, cdtype, pahu, pahv, zglu, zglv, zgui, zgvi, zlap, ptb, pta, kjpt, 2 )
      END SELECT
      !
   END SUBROUTINE tra_ldf_blp

   !!==============================================================================
END MODULE traldf_lap_blp
