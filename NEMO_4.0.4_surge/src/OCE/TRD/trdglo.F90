MODULE trdglo
   !!======================================================================
   !!                       ***  MODULE  trdglo  ***
   !! Ocean diagnostics:  global domain averaged tracer and momentum trends
   !!=====================================================================
   !! History :  1.0  !  2004-08  (C. Talandier) New trends organization
   !!            3.5  !  2012-02  (G. Madec)  add 3D tracer zdf trend output using iom
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   trd_glo       : domain averaged budget of trends (including kinetic energy and T^2 trends)
   !!   glo_dyn_wri   : print dynamic trends in ocean.output file
   !!   glo_tra_wri   : print global T & T^2 trends in ocean.output file
   !!   trd_glo_init  : initialization step
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE sbc_oce        ! surface boundary condition: ocean
   USE trd_oce        ! trends: ocean variables
   USE phycst         ! physical constants
   USE ldftra         ! lateral diffusion: eddy diffusivity & EIV coeff.
   USE ldfdyn         ! ocean dynamics: lateral physics
   USE zdf_oce        ! ocean vertical physics
!!gm   USE zdfdrg         ! ocean vertical physics: bottom friction
   USE zdfddm         ! ocean vertical physics: double diffusion
   USE eosbn2         ! equation of state
   USE phycst         ! physical constants
   !
   USE lib_mpp        ! distibuted memory computing library
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trd_glo       ! called by trdtra and trddyn modules
   PUBLIC   trd_glo_init  ! called by trdini module

   !                     !!! Variables used for diagnostics
   REAL(wp) ::   tvolt    ! volume of the whole ocean computed at t-points
   REAL(wp) ::   tvolu    ! volume of the whole ocean computed at u-points
   REAL(wp) ::   tvolv    ! volume of the whole ocean computed at v-points
   REAL(wp) ::   rpktrd   ! potential to kinetic energy conversion
   REAL(wp) ::   peke     ! conversion potential energy - kinetic energy trend

   !                     !!! domain averaged trends
   REAL(wp), DIMENSION(jptot_tra) ::   tmo, smo   ! temperature and salinity trends 
   REAL(wp), DIMENSION(jptot_tra) ::   t2 , s2    ! T^2 and S^2 trends 
   REAL(wp), DIMENSION(jptot_dyn) ::   umo, vmo   ! momentum trends 
   REAL(wp), DIMENSION(jptot_dyn) ::   hke        ! kinetic energy trends (u^2+v^2) 

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trd_glo( ptrdx, ptrdy, ktrd, ctype, kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_glo  ***
      !! 
      !! ** Purpose :   compute and print global domain averaged trends for 
      !!              T, T^2, momentum, KE, and KE<->PE
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptrdx   ! Temperature or U trend 
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptrdy   ! Salinity    or V trend
      INTEGER                   , INTENT(in   ) ::   ktrd    ! tracer trend index
      CHARACTER(len=3)          , INTENT(in   ) ::   ctype   ! momentum or tracers trends type (='DYN'/'TRA')
      INTEGER                   , INTENT(in   ) ::   kt      ! time step
      !!
      INTEGER ::   ji, jj, jk      ! dummy loop indices
      INTEGER ::   ikbu, ikbv      ! local integers
      REAL(wp)::   zvm, zvt, zvs, z1_2rau0   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj)  :: ztswu, ztswv, z2dx, z2dy   ! 2D workspace 
      !!----------------------------------------------------------------------
      !
      IF( MOD(kt,nn_trd) == 0 .OR. kt == nit000 .OR. kt == nitend ) THEN
         !
      	SELECT CASE( ctype )
      	!
      	CASE( 'TRA' )          !==  Tracers (T & S)  ==!
      	   DO jk = 1, jpkm1       ! global sum of mask volume trend and trend*T (including interior mask)
      	      DO jj = 1, jpj
      	         DO ji = 1, jpi        
      	            zvm = e1e2t(ji,jj) * e3t_n(ji,jj,jk) * tmask(ji,jj,jk) * tmask_i(ji,jj)
      	            zvt = ptrdx(ji,jj,jk) * zvm
      	            zvs = ptrdy(ji,jj,jk) * zvm
         	         tmo(ktrd) = tmo(ktrd) + zvt   
         	         smo(ktrd) = smo(ktrd) + zvs
         	         t2 (ktrd) = t2(ktrd)  + zvt * tsn(ji,jj,jk,jp_tem)
         	         s2 (ktrd) = s2(ktrd)  + zvs * tsn(ji,jj,jk,jp_sal)
                  END DO
               END DO
            END DO
            !                       ! linear free surface: diagnose advective flux trough the fixed k=1 w-surface
            IF( ln_linssh .AND. ktrd == jptra_zad ) THEN  
               tmo(jptra_sad) = SUM( wn(:,:,1) * tsn(:,:,1,jp_tem) * e1e2t(:,:) * tmask_i(:,:) )
               smo(jptra_sad) = SUM( wn(:,:,1) * tsn(:,:,1,jp_sal) * e1e2t(:,:) * tmask_i(:,:)  )
               t2 (jptra_sad) = SUM( wn(:,:,1) * tsn(:,:,1,jp_tem) * tsn(:,:,1,jp_tem) * e1e2t(:,:) * tmask_i(:,:)  )
               s2 (jptra_sad) = SUM( wn(:,:,1) * tsn(:,:,1,jp_sal) * tsn(:,:,1,jp_sal) * e1e2t(:,:) * tmask_i(:,:)  )
            ENDIF
            !
            IF( ktrd == jptra_atf ) THEN     ! last trend (asselin time filter)
               ! 
               CALL glo_tra_wri( kt )             ! print the results in ocean.output
               !                
	            tmo(:) = 0._wp                     ! prepare the next time step (domain averaged array reset to zero)
   	         smo(:) = 0._wp
               t2 (:) = 0._wp
               s2 (:) = 0._wp
               !
            ENDIF
            !
         CASE( 'DYN' )          !==  Momentum and KE  ==!        
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     zvt = ptrdx(ji,jj,jk) * tmask_i(ji+1,jj) * tmask_i(ji,jj) * umask(ji,jj,jk)   &
                        &                                     * e1e2u  (ji,jj) * e3u_n(ji,jj,jk)
                     zvs = ptrdy(ji,jj,jk) * tmask_i(ji,jj+1) * tmask_i(ji,jj) * vmask(ji,jj,jk)   &
                        &                                     * e1e2v  (ji,jj) * e3u_n(ji,jj,jk)
                     umo(ktrd) = umo(ktrd) + zvt
                     vmo(ktrd) = vmo(ktrd) + zvs
                     hke(ktrd) = hke(ktrd) + un(ji,jj,jk) * zvt + vn(ji,jj,jk) * zvs
                  END DO
               END DO
            END DO
            !                 
            IF( ktrd == jpdyn_zdf ) THEN      ! zdf trend: compute separately the surface forcing trend
               z1_2rau0 = 0.5_wp / rau0
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     zvt = ( utau_b(ji,jj) + utau(ji,jj) ) * tmask_i(ji+1,jj) * tmask_i(ji,jj) * umask(ji,jj,jk)   &
                        &                                                     * z1_2rau0       * e1e2u(ji,jj)
                     zvs = ( vtau_b(ji,jj) + vtau(ji,jj) ) * tmask_i(ji,jj+1) * tmask_i(ji,jj) * vmask(ji,jj,jk)   &
                        &                                                     * z1_2rau0       * e1e2v(ji,jj)
                     umo(jpdyn_tau) = umo(jpdyn_tau) + zvt
                     vmo(jpdyn_tau) = vmo(jpdyn_tau) + zvs
                     hke(jpdyn_tau) = hke(jpdyn_tau) + un(ji,jj,1) * zvt + vn(ji,jj,1) * zvs
                  END DO
               END DO
            ENDIF
            !                       
!!gm  miss placed calculation   ===>>>> to be done in dynzdf.F90
!            IF( ktrd == jpdyn_atf ) THEN     ! last trend (asselin time filter)
!               !
!               IF( ln_drgimp ) THEN                   ! implicit drag case: compute separately the bottom friction 
!                  z1_2rau0 = 0.5_wp / rau0
!                  DO jj = 1, jpjm1
!                     DO ji = 1, jpim1
!                        ikbu = mbku(ji,jj)                  ! deepest ocean u- & v-levels
!                        ikbv = mbkv(ji,jj)
!                        zvt = 0.5*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) ) * un(ji,jj,ikbu) * e1e2u(ji,jj)
!                        zvs = 0.5*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) ) * vn(ji,jj,ikbv) * e1e2v(ji,jj)
!                        umo(jpdyn_bfri) = umo(jpdyn_bfri) + zvt
!                        vmo(jpdyn_bfri) = vmo(jpdyn_bfri) + zvs
!                        hke(jpdyn_bfri) = hke(jpdyn_bfri) + un(ji,jj,ikbu) * zvt + vn(ji,jj,ikbv) * zvs
!                     END DO
!                  END DO
!               ENDIF
!
!!gm top drag case is missing 
!
!               ! 
!               CALL glo_dyn_wri( kt )                 ! print the results in ocean.output
!               !                
!               umo(:) = 0._wp                         ! reset for the next time step
!               vmo(:) = 0._wp
!               hke(:) = 0._wp
!               !
!            ENDIF
!!gm end
            !
         END SELECT
         !
      ENDIF
      !
   END SUBROUTINE trd_glo


   SUBROUTINE glo_dyn_wri( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE glo_dyn_wri  ***
      !! 
      !! ** Purpose :  write global averaged U, KE, PE<->KE trends in ocean.output 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zcof         ! local scalar
      REAL(wp), DIMENSION(jpi,jpj,jpk)  ::  zkx, zky, zkz, zkepe  
      !!----------------------------------------------------------------------

      ! I. Momentum trends
      ! -------------------

      IF( MOD( kt, nn_trd ) == 0 .OR. kt == nit000 .OR. kt == nitend ) THEN

         ! I.1 Conversion potential energy - kinetic energy
         ! --------------------------------------------------
         ! c a u t i o n here, trends are computed at kt+1 (now , but after the swap)
         zkx  (:,:,:) = 0._wp
         zky  (:,:,:) = 0._wp
         zkz  (:,:,:) = 0._wp
         zkepe(:,:,:) = 0._wp
   
         CALL eos( tsn, rhd, rhop )       ! now potential density

         zcof = 0.5_wp / rau0             ! Density flux at w-point
         zkz(:,:,1) = 0._wp
         DO jk = 2, jpk
            zkz(:,:,jk) = zcof * e1e2t(:,:) * wn(:,:,jk) * ( rhop(:,:,jk) + rhop(:,:,jk-1) ) * tmask_i(:,:)
         END DO
         
         zcof   = 0.5_wp / rau0           ! Density flux at u and v-points
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zkx(ji,jj,jk) = zcof * e2u(ji,jj) * e3u_n(ji,jj,jk) * un(ji,jj,jk) * ( rhop(ji,jj,jk) + rhop(ji+1,jj,jk) )
                  zky(ji,jj,jk) = zcof * e1v(ji,jj) * e3v_n(ji,jj,jk) * vn(ji,jj,jk) * ( rhop(ji,jj,jk) + rhop(ji,jj+1,jk) )
               END DO
            END DO
         END DO
         
         DO jk = 1, jpkm1                 ! Density flux divergence at t-point
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zkepe(ji,jj,jk) = - (  zkz(ji,jj,jk) - zkz(ji  ,jj  ,jk+1)               &
                     &                 + zkx(ji,jj,jk) - zkx(ji-1,jj  ,jk  )               &
                     &                 + zky(ji,jj,jk) - zky(ji  ,jj-1,jk  )   )           &
                     &              / ( e1e2t(ji,jj) * e3t_n(ji,jj,jk) ) * tmask(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO

         ! I.2 Basin averaged kinetic energy trend
         ! ----------------------------------------
         peke = 0._wp
         DO jk = 1, jpkm1
            peke = peke + SUM( zkepe(:,:,jk) * gdept_n(:,:,jk) * e1e2t(:,:) * e3t_n(:,:,jk) )
         END DO
         peke = grav * peke

         ! I.3 Sums over the global domain
         ! ---------------------------------
         IF( lk_mpp ) THEN
            CALL mpp_sum( 'trdglo', peke )
            CALL mpp_sum( 'trdglo', umo , jptot_dyn )
            CALL mpp_sum( 'trdglo', vmo , jptot_dyn )
            CALL mpp_sum( 'trdglo', hke , jptot_dyn )
         ENDIF

         ! I.2 Print dynamic trends in the ocean.output file
         ! --------------------------------------------------

         IF(lwp) THEN
            WRITE (numout,*)
            WRITE (numout,*)
            WRITE (numout,9500) kt
            WRITE (numout,9501) umo(jpdyn_hpg) / tvolu, vmo(jpdyn_hpg) / tvolv
            WRITE (numout,9502) umo(jpdyn_keg) / tvolu, vmo(jpdyn_keg) / tvolv
            WRITE (numout,9503) umo(jpdyn_rvo) / tvolu, vmo(jpdyn_rvo) / tvolv
            WRITE (numout,9504) umo(jpdyn_pvo) / tvolu, vmo(jpdyn_pvo) / tvolv
            WRITE (numout,9505) umo(jpdyn_zad) / tvolu, vmo(jpdyn_zad) / tvolv
            WRITE (numout,9506) umo(jpdyn_ldf) / tvolu, vmo(jpdyn_ldf) / tvolv
            WRITE (numout,9507) umo(jpdyn_zdf) / tvolu, vmo(jpdyn_zdf) / tvolv
            WRITE (numout,9508) umo(jpdyn_spg) / tvolu, vmo(jpdyn_spg) / tvolv
            WRITE (numout,9509) umo(jpdyn_bfr) / tvolu, vmo(jpdyn_bfr) / tvolv
            WRITE (numout,9510) umo(jpdyn_atf) / tvolu, vmo(jpdyn_atf) / tvolv
            WRITE (numout,9511)
            WRITE (numout,9512)                                                 &
            &     (  umo(jpdyn_hpg) + umo(jpdyn_keg) + umo(jpdyn_rvo) + umo(jpdyn_pvo)   &
            &      + umo(jpdyn_zad) + umo(jpdyn_ldf) + umo(jpdyn_zdf) + umo(jpdyn_spg)   &
            &      + umo(jpdyn_bfr) + umo(jpdyn_atf) ) / tvolu,   &
            &     (  vmo(jpdyn_hpg) + vmo(jpdyn_keg) + vmo(jpdyn_rvo) + vmo(jpdyn_pvo)   &
            &      + vmo(jpdyn_zad) + vmo(jpdyn_ldf) + vmo(jpdyn_zdf) + vmo(jpdyn_spg)   &
            &      + vmo(jpdyn_bfr) + vmo(jpdyn_atf) ) / tvolv
            WRITE (numout,9513) umo(jpdyn_tau) / tvolu, vmo(jpdyn_tau) / tvolv
!!gm            IF( ln_drgimp )   WRITE (numout,9514) umo(jpdyn_bfri) / tvolu, vmo(jpdyn_bfri) / tvolv
         ENDIF

 9500    FORMAT(' momentum trend at it= ', i6, ' :', /' ==============================')
 9501    FORMAT(' hydro pressure gradient    u= ', e20.13, '    v= ', e20.13)
 9502    FORMAT(' ke gradient                u= ', e20.13, '    v= ', e20.13)
 9503    FORMAT(' relative vorticity term    u= ', e20.13, '    v= ', e20.13)
 9504    FORMAT(' planetary vorticity term   u= ', e20.13, '    v= ', e20.13)
 9505    FORMAT(' vertical advection         u= ', e20.13, '    v= ', e20.13)
 9506    FORMAT(' horizontal diffusion       u= ', e20.13, '    v= ', e20.13)
 9507    FORMAT(' vertical diffusion         u= ', e20.13, '    v= ', e20.13)
 9508    FORMAT(' surface pressure gradient  u= ', e20.13, '    v= ', e20.13)
 9509    FORMAT(' explicit bottom friction   u= ', e20.13, '    v= ', e20.13)
 9510    FORMAT(' Asselin time filter        u= ', e20.13, '    v= ', e20.13)
 9511    FORMAT(' -----------------------------------------------------------------------------')
 9512    FORMAT(' total trend                u= ', e20.13, '    v= ', e20.13)
 9513    FORMAT(' incl. surface wind stress  u= ', e20.13, '    v= ', e20.13)
 9514    FORMAT('       bottom stress        u= ', e20.13, '    v= ', e20.13)

         IF(lwp) THEN
            WRITE (numout,*)
            WRITE (numout,*)
            WRITE (numout,9520) kt
            WRITE (numout,9521) hke(jpdyn_hpg) / tvolt
            WRITE (numout,9522) hke(jpdyn_keg) / tvolt
            WRITE (numout,9523) hke(jpdyn_rvo) / tvolt
            WRITE (numout,9524) hke(jpdyn_pvo) / tvolt
            WRITE (numout,9525) hke(jpdyn_zad) / tvolt
            WRITE (numout,9526) hke(jpdyn_ldf) / tvolt
            WRITE (numout,9527) hke(jpdyn_zdf) / tvolt
            WRITE (numout,9528) hke(jpdyn_spg) / tvolt
            WRITE (numout,9529) hke(jpdyn_bfr) / tvolt
            WRITE (numout,9530) hke(jpdyn_atf) / tvolt
            WRITE (numout,9531)
            WRITE (numout,9532)   &
            &     (  hke(jpdyn_hpg) + hke(jpdyn_keg) + hke(jpdyn_rvo) + hke(jpdyn_pvo)   &
            &      + hke(jpdyn_zad) + hke(jpdyn_ldf) + hke(jpdyn_zdf) + hke(jpdyn_spg)   &
            &      + hke(jpdyn_bfr) + hke(jpdyn_atf) ) / tvolt
            WRITE (numout,9533) hke(jpdyn_tau) / tvolt
!!gm            IF( ln_drgimp )   WRITE (numout,9534) hke(jpdyn_bfri) / tvolt
         ENDIF

 9520    FORMAT(' kinetic energy trend at it= ', i6, ' :', /' ====================================')
 9521    FORMAT(' hydro pressure gradient   u2= ', e20.13)
 9522    FORMAT(' ke gradient               u2= ', e20.13)
 9523    FORMAT(' relative vorticity term   u2= ', e20.13)
 9524    FORMAT(' planetary vorticity term  u2= ', e20.13)
 9525    FORMAT(' vertical advection        u2= ', e20.13)
 9526    FORMAT(' horizontal diffusion      u2= ', e20.13)
 9527    FORMAT(' vertical diffusion        u2= ', e20.13)
 9528    FORMAT(' surface pressure gradient u2= ', e20.13)
 9529    FORMAT(' explicit bottom friction  u2= ', e20.13)
 9530    FORMAT(' Asselin time filter       u2= ', e20.13)
 9531    FORMAT(' --------------------------------------------------')
 9532    FORMAT(' total trend               u2= ', e20.13)
 9533    FORMAT(' incl. surface wind stress u2= ', e20.13)
 9534    FORMAT('       bottom stress       u2= ', e20.13)

         IF(lwp) THEN
            WRITE (numout,*)
            WRITE (numout,*)
            WRITE (numout,9540) kt
            WRITE (numout,9541) ( hke(jpdyn_keg) + hke(jpdyn_rvo) + hke(jpdyn_zad) ) / tvolt
            WRITE (numout,9542) ( hke(jpdyn_keg) + hke(jpdyn_zad) ) / tvolt
            WRITE (numout,9543) ( hke(jpdyn_pvo) ) / tvolt
            WRITE (numout,9544) ( hke(jpdyn_rvo) ) / tvolt
            WRITE (numout,9545) ( hke(jpdyn_spg) ) / tvolt
            WRITE (numout,9546) ( hke(jpdyn_ldf) ) / tvolt
            WRITE (numout,9547) ( hke(jpdyn_zdf) ) / tvolt
            WRITE (numout,9548) ( hke(jpdyn_hpg) ) / tvolt, rpktrd / tvolt
            WRITE (numout,*)
            WRITE (numout,*)
         ENDIF

 9540    FORMAT(' energetic consistency at it= ', i6, ' :', /' =========================================')
 9541    FORMAT(' 0 = non linear term (true if KE conserved)                : ', e20.13)
 9542    FORMAT(' 0 = ke gradient + vertical advection                      : ', e20.13)
 9543    FORMAT(' 0 = coriolis term  (true if KE conserving scheme)         : ', e20.13)
 9544    FORMAT(' 0 = vorticity term (true if KE conserving scheme)         : ', e20.13)
 9545    FORMAT(' 0 = surface pressure gradient  ???                        : ', e20.13)
 9546    FORMAT(' 0 < horizontal diffusion                                  : ', e20.13)
 9547    FORMAT(' 0 < vertical diffusion                                    : ', e20.13)
 9548    FORMAT(' pressure gradient u2 = - 1/rau0 u.dz(rhop)                : ', e20.13, '  u.dz(rhop) =', e20.13)
         !
         ! Save potential to kinetic energy conversion for next time step
         rpktrd = peke
         !
      ENDIF
      !
   END SUBROUTINE glo_dyn_wri


   SUBROUTINE glo_tra_wri( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE glo_tra_wri  ***
      !! 
      !! ** Purpose :  write global domain averaged of T and T^2 trends in ocean.output 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   jk   ! loop indices
      !!----------------------------------------------------------------------

      ! I. Tracers trends
      ! -----------------

      IF( MOD(kt,nn_trd) == 0 .OR. kt == nit000 .OR. kt == nitend ) THEN

         ! I.1 Sums over the global domain
         ! -------------------------------
         IF( lk_mpp ) THEN
            CALL mpp_sum( 'trdglo', tmo, jptot_tra )   
            CALL mpp_sum( 'trdglo', smo, jptot_tra )
            CALL mpp_sum( 'trdglo', t2 , jptot_tra )
            CALL mpp_sum( 'trdglo', s2 , jptot_tra )
         ENDIF

         ! I.2 Print tracers trends in the ocean.output file
         ! --------------------------------------------------
         
         IF(lwp) THEN
            WRITE (numout,*)
            WRITE (numout,*)
            WRITE (numout,9400) kt
            WRITE (numout,9401)  tmo(jptra_xad) / tvolt, smo(jptra_xad) / tvolt
            WRITE (numout,9411)  tmo(jptra_yad) / tvolt, smo(jptra_yad) / tvolt
            WRITE (numout,9402)  tmo(jptra_zad) / tvolt, smo(jptra_zad) / tvolt
            WRITE (numout,9403)  tmo(jptra_ldf) / tvolt, smo(jptra_ldf) / tvolt
            WRITE (numout,9404)  tmo(jptra_zdf) / tvolt, smo(jptra_zdf) / tvolt
            WRITE (numout,9405)  tmo(jptra_npc) / tvolt, smo(jptra_npc) / tvolt
            WRITE (numout,9406)  tmo(jptra_dmp) / tvolt, smo(jptra_dmp) / tvolt
            WRITE (numout,9407)  tmo(jptra_qsr) / tvolt
            WRITE (numout,9408)  tmo(jptra_nsr) / tvolt, smo(jptra_nsr) / tvolt
            WRITE (numout,9409) 
            WRITE (numout,9410) (  tmo(jptra_xad) + tmo(jptra_yad) + tmo(jptra_zad) + tmo(jptra_ldf) + tmo(jptra_zdf)   &
            &                    + tmo(jptra_npc) + tmo(jptra_dmp) + tmo(jptra_qsr) + tmo(jptra_nsr) ) / tvolt,   &
            &                   (  smo(jptra_xad) + smo(jptra_yad) + smo(jptra_zad) + smo(jptra_ldf) + smo(jptra_zdf)   &
            &                    + smo(jptra_npc) + smo(jptra_dmp)                   + smo(jptra_nsr) ) / tvolt
         ENDIF

9400     FORMAT(' tracer trend at it= ',i6,' :     temperature',   &
              '              salinity',/' ============================')
9401     FORMAT(' zonal      advection        ',e20.13,'     ',e20.13)
9411     FORMAT(' meridional advection        ',e20.13,'     ',e20.13)
9402     FORMAT(' vertical advection          ',e20.13,'     ',e20.13)
9403     FORMAT(' horizontal diffusion        ',e20.13,'     ',e20.13)
9404     FORMAT(' vertical diffusion          ',e20.13,'     ',e20.13)
9405     FORMAT(' static instability mixing   ',e20.13,'     ',e20.13)
9406     FORMAT(' damping term                ',e20.13,'     ',e20.13)
9407     FORMAT(' penetrative qsr             ',e20.13)
9408     FORMAT(' non solar radiation         ',e20.13,'     ',e20.13)
9409     FORMAT(' -------------------------------------------------------------------------')
9410     FORMAT(' total trend                 ',e20.13,'     ',e20.13)


         IF(lwp) THEN
            WRITE (numout,*)
            WRITE (numout,*)
            WRITE (numout,9420) kt
            WRITE (numout,9421)   t2(jptra_xad) / tvolt, s2(jptra_xad) / tvolt
            WRITE (numout,9431)   t2(jptra_yad) / tvolt, s2(jptra_yad) / tvolt
            WRITE (numout,9422)   t2(jptra_zad) / tvolt, s2(jptra_zad) / tvolt
            WRITE (numout,9423)   t2(jptra_ldf) / tvolt, s2(jptra_ldf) / tvolt
            WRITE (numout,9424)   t2(jptra_zdf) / tvolt, s2(jptra_zdf) / tvolt
            WRITE (numout,9425)   t2(jptra_npc) / tvolt, s2(jptra_npc) / tvolt
            WRITE (numout,9426)   t2(jptra_dmp) / tvolt, s2(jptra_dmp) / tvolt
            WRITE (numout,9427)   t2(jptra_qsr) / tvolt
            WRITE (numout,9428)   t2(jptra_nsr) / tvolt, s2(jptra_nsr) / tvolt
            WRITE (numout,9429)
            WRITE (numout,9430) (  t2(jptra_xad) + t2(jptra_yad) + t2(jptra_zad) + t2(jptra_ldf) + t2(jptra_zdf)   &
            &                    + t2(jptra_npc) + t2(jptra_dmp) + t2(jptra_qsr) + t2(jptra_nsr) ) / tvolt,   &
            &                   (  s2(jptra_xad) + s2(jptra_yad) + s2(jptra_zad) + s2(jptra_ldf) + s2(jptra_zdf)   &
            &                    + s2(jptra_npc) + s2(jptra_dmp)                  + s2(jptra_nsr) ) / tvolt
         ENDIF

9420     FORMAT(' tracer**2 trend at it= ', i6, ' :      temperature',   &
            '               salinity', /, ' ===============================')
9421     FORMAT(' zonal      advection      * t   ', e20.13, '     ', e20.13)
9431     FORMAT(' meridional advection      * t   ', e20.13, '     ', e20.13)
9422     FORMAT(' vertical advection        * t   ', e20.13, '     ', e20.13)
9423     FORMAT(' horizontal diffusion      * t   ', e20.13, '     ', e20.13)
9424     FORMAT(' vertical diffusion        * t   ', e20.13, '     ', e20.13)
9425     FORMAT(' static instability mixing * t   ', e20.13, '     ', e20.13)
9426     FORMAT(' damping term              * t   ', e20.13, '     ', e20.13)
9427     FORMAT(' penetrative qsr           * t   ', e20.13)
9428     FORMAT(' non solar radiation       * t   ', e20.13, '     ', e20.13)
9429     FORMAT(' -----------------------------------------------------------------------------')
9430     FORMAT(' total trend                *t = ', e20.13, '  *s = ', e20.13)


         IF(lwp) THEN
            WRITE (numout,*)
            WRITE (numout,*)
            WRITE (numout,9440) kt
            WRITE (numout,9441) ( tmo(jptra_xad)+tmo(jptra_yad)+tmo(jptra_zad) )/tvolt,   &
            &                   ( smo(jptra_xad)+smo(jptra_yad)+smo(jptra_zad) )/tvolt
            WRITE (numout,9442)   tmo(jptra_sad)/tvolt,  smo(jptra_sad)/tvolt
            WRITE (numout,9443)   tmo(jptra_ldf)/tvolt,  smo(jptra_ldf)/tvolt
            WRITE (numout,9444)   tmo(jptra_zdf)/tvolt,  smo(jptra_zdf)/tvolt
            WRITE (numout,9445)   tmo(jptra_npc)/tvolt,  smo(jptra_npc)/tvolt
            WRITE (numout,9446) ( t2(jptra_xad)+t2(jptra_yad)+t2(jptra_zad) )/tvolt,    &
            &                   ( s2(jptra_xad)+s2(jptra_yad)+s2(jptra_zad) )/tvolt
            WRITE (numout,9447)   t2(jptra_ldf)/tvolt,   s2(jptra_ldf)/tvolt
            WRITE (numout,9448)   t2(jptra_zdf)/tvolt,   s2(jptra_zdf)/tvolt
            WRITE (numout,9449)   t2(jptra_npc)/tvolt,   s2(jptra_npc)/tvolt
         ENDIF

9440     FORMAT(' tracer consistency at it= ',i6,   &
            ' :         temperature','                salinity',/,   &
            ' ==================================')
9441     FORMAT(' 0 = horizontal+vertical advection +    ',e20.13,'       ',e20.13)
9442     FORMAT('     1st lev vertical advection         ',e20.13,'       ',e20.13)
9443     FORMAT(' 0 = horizontal diffusion               ',e20.13,'       ',e20.13)
9444     FORMAT(' 0 = vertical diffusion                 ',e20.13,'       ',e20.13)
9445     FORMAT(' 0 = static instability mixing          ',e20.13,'       ',e20.13)
9446     FORMAT(' 0 = horizontal+vertical advection * t  ',e20.13,'       ',e20.13)
9447     FORMAT(' 0 > horizontal diffusion          * t  ',e20.13,'       ',e20.13)
9448     FORMAT(' 0 > vertical diffusion            * t  ',e20.13,'       ',e20.13)
9449     FORMAT(' 0 > static instability mixing     * t  ',e20.13,'       ',e20.13)
         !
      ENDIF
      !
   END SUBROUTINE glo_tra_wri


   SUBROUTINE trd_glo_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_glo_init  ***
      !! 
      !! ** Purpose :   Read the namtrd namelist
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trd_glo_init : integral constraints properties trends'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF

      ! Total volume at t-points:
      tvolt = 0._wp
      DO jk = 1, jpkm1
         tvolt = tvolt + SUM( e1e2t(:,:) * e3t_n(:,:,jk) * tmask(:,:,jk) * tmask_i(:,:) )
      END DO
      CALL mpp_sum( 'trdglo', tvolt )   ! sum over the global domain

      IF(lwp) WRITE(numout,*) '                total ocean volume at T-point   tvolt = ',tvolt

      ! Initialization of potential to kinetic energy conversion
      rpktrd = 0._wp

      ! Total volume at u-, v- points:
!!gm :  bug?  je suis quasi sur que le produit des tmask_i ne correspond pas exactement au umask_i et vmask_i !
      tvolu = 0._wp
      tvolv = 0._wp

      DO jk = 1, jpk
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               tvolu = tvolu + e1u(ji,jj) * e2u(ji,jj) * e3u_n(ji,jj,jk) * tmask_i(ji+1,jj  ) * tmask_i(ji,jj) * umask(ji,jj,jk)
               tvolv = tvolv + e1v(ji,jj) * e2v(ji,jj) * e3v_n(ji,jj,jk) * tmask_i(ji  ,jj+1) * tmask_i(ji,jj) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL mpp_sum( 'trdglo', tvolu )   ! sums over the global domain
      CALL mpp_sum( 'trdglo', tvolv )

      IF(lwp) THEN
         WRITE(numout,*) '                total ocean volume at U-point   tvolu = ',tvolu
         WRITE(numout,*) '                total ocean volume at V-point   tvolv = ',tvolv
      ENDIF
      !
   END SUBROUTINE trd_glo_init

   !!======================================================================
END MODULE trdglo
