MODULE crsfld
   !!======================================================================
   !!                     ***  MODULE  crsdfld  ***
   !!  Ocean coarsening :  coarse ocean fields
   !!=====================================================================
   !!   2012-07  (J. Simeon, C. Calone, G. Madec, C. Ethe)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   crs_fld       : create the standard output files for coarse grid and prep
   !!                       other variables needed to be passed to TOP
   !!----------------------------------------------------------------------
   USE crs
   USE crsdom
   USE crslbclnk
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE zdf_oce         ! vertical  physics: ocean fields
   USE ldftra          ! ocean active tracers: lateral diffusivity & EIV coefficients
   USE zdfddm          ! vertical  physics: double diffusion
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! 
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE timing          ! preformance summary

   IMPLICIT NONE
   PRIVATE

   PUBLIC   crs_fld                 ! routines called by step.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE crs_fld( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE crs_fld  ***
      !!                   
      !! ** Purpose :   Basic output of coarsened dynamics and tracer fields 
      !!      NETCDF format is used by default 
      !!      1. Accumulate in time the dimensionally-weighted fields
      !!      2. At time of output, rescale [1] by dimension and time
      !!         to yield the spatial and temporal average.
      !!  See. sbcmod.F90
      !!
      !! ** Method  :  
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk        ! dummy loop indices
      REAL(wp) ::   z2dcrsu, z2dcrsv  ! local scalars
      REAL(wp) ::   zztmp             !   -      -
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ze3t, ze3u, ze3v, ze3w   ! 3D workspace for e3
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zt  , zs  , z3d
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk) ::   zt_crs, zs_crs  
      !!----------------------------------------------------------------------
      ! 
      IF( ln_timing )   CALL timing_start('crs_fld')

      ! Depth work arrrays
      ze3t(:,:,:) = e3t_n(:,:,:)
      ze3u(:,:,:) = e3u_n(:,:,:)
      ze3v(:,:,:) = e3v_n(:,:,:)
      ze3w(:,:,:) = e3w_n(:,:,:)

      IF( kt == nit000  ) THEN
         tsn_crs  (:,:,:,:) = 0._wp    ! temp/sal  array, now 
         un_crs   (:,:,:  ) = 0._wp    ! u-velocity
         vn_crs   (:,:,:  ) = 0._wp    ! v-velocity
         wn_crs   (:,:,:  ) = 0._wp    ! w
         avs_crs  (:,:,:  ) = 0._wp    ! avt
         hdivn_crs(:,:,:  ) = 0._wp    ! hdiv
         sshn_crs (:,:    ) = 0._wp    ! ssh
         utau_crs (:,:    ) = 0._wp    ! taux
         vtau_crs (:,:    ) = 0._wp    ! tauy
         wndm_crs (:,:    ) = 0._wp    ! wind speed
         qsr_crs  (:,:    ) = 0._wp    ! qsr
         emp_crs  (:,:    ) = 0._wp    ! emp
         emp_b_crs(:,:    ) = 0._wp    ! emp
         rnf_crs  (:,:    ) = 0._wp    ! runoff
         fr_i_crs (:,:    ) = 0._wp    ! ice cover
      ENDIF

      CALL iom_swap( "nemo_crs" )    ! swap on the coarse grid

      ! 2. Coarsen fields at each time step
      ! --------------------------------------------------------

      !  Temperature
      zt(:,:,:) = tsn(:,:,:,jp_tem)  ;      zt_crs(:,:,:) = 0._wp
      CALL crs_dom_ope( zt, 'VOL', 'T', tmask, zt_crs, p_e12=e1e2t, p_e3=ze3t, psgn=1.0 )
      tsn_crs(:,:,:,jp_tem) = zt_crs(:,:,:)

      CALL iom_put( "toce", tsn_crs(:,:,:,jp_tem) )    ! temp
      CALL iom_put( "sst" , tsn_crs(:,:,1,jp_tem) )    ! sst

      
      !  Salinity
      zs(:,:,:) = tsn(:,:,:,jp_sal)  ;      zs_crs(:,:,:) = 0._wp
      CALL crs_dom_ope( zs, 'VOL', 'T', tmask, zs_crs, p_e12=e1e2t, p_e3=ze3t, psgn=1.0 )
      tsn_crs(:,:,:,jp_sal) = zt_crs(:,:,:)

      CALL iom_put( "soce" , tsn_crs(:,:,:,jp_sal) )    ! sal
      CALL iom_put( "sss"  , tsn_crs(:,:,1,jp_sal) )    ! sss

      !  U-velocity
      CALL crs_dom_ope( un, 'SUM', 'U', umask, un_crs, p_e12=e2u, p_e3=ze3u, p_surf_crs=e2e3u_msk, psgn=-1.0 )
      !
      zt(:,:,:) = 0._wp     ;    zs(:,:,:) = 0._wp  ;   zt_crs(:,:,:) = 0._wp   ;    zs_crs(:,:,:) = 0._wp
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   
               zt(ji,jj,jk)  = un(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_tem) + tsn(ji+1,jj,jk,jp_tem) ) 
               zs(ji,jj,jk)  = un(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_sal) + tsn(ji+1,jj,jk,jp_sal) ) 
            END DO
         END DO
      END DO
      CALL crs_dom_ope( zt, 'SUM', 'U', umask, zt_crs, p_e12=e2u, p_e3=ze3u, p_surf_crs=e2e3u_msk, psgn=-1.0 )
      CALL crs_dom_ope( zs, 'SUM', 'U', umask, zs_crs, p_e12=e2u, p_e3=ze3u, p_surf_crs=e2e3u_msk, psgn=-1.0 )

      CALL iom_put( "uoce"  , un_crs )   ! i-current 
      CALL iom_put( "uocet" , zt_crs )   ! uT
      CALL iom_put( "uoces" , zs_crs )   ! uS

      !  V-velocity
      CALL crs_dom_ope( vn, 'SUM', 'V', vmask, vn_crs, p_e12=e1v, p_e3=ze3v, p_surf_crs=e1e3v_msk, psgn=-1.0 )
      !                                                                                 
      zt(:,:,:) = 0._wp     ;    zs(:,:,:) = 0._wp  ;   zt_crs(:,:,:) = 0._wp   ;    zs_crs(:,:,:) = 0._wp
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   
               zt(ji,jj,jk)  = vn(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_tem) + tsn(ji,jj+1,jk,jp_tem) ) 
               zs(ji,jj,jk)  = vn(ji,jj,jk) * 0.5 * ( tsn(ji,jj,jk,jp_sal) + tsn(ji,jj+1,jk,jp_sal) ) 
            END DO
         END DO
      END DO
      CALL crs_dom_ope( zt, 'SUM', 'V', vmask, zt_crs, p_e12=e1v, p_e3=ze3v, p_surf_crs=e1e3v_msk, psgn=-1.0 )
      CALL crs_dom_ope( zs, 'SUM', 'V', vmask, zs_crs, p_e12=e1v, p_e3=ze3v, p_surf_crs=e1e3v_msk, psgn=-1.0 )
 
      CALL iom_put( "voce"  , vn_crs )   ! i-current 
      CALL iom_put( "vocet" , zt_crs )   ! vT
      CALL iom_put( "voces" , zs_crs )   ! vS

      IF( iom_use( "eken") ) THEN     !      kinetic energy
         z3d(:,:,jk) = 0._wp 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zztmp  = r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk)
                  z3d(ji,jj,jk) = 0.25_wp * zztmp * (                                    &
                     &            un(ji-1,jj,jk)**2 * e2u(ji-1,jj) * e3u_n(ji-1,jj,jk)   &
                     &          + un(ji  ,jj,jk)**2 * e2u(ji  ,jj) * e3u_n(ji  ,jj,jk)   &
                     &          + vn(ji,jj-1,jk)**2 * e1v(ji,jj-1) * e3v_n(ji,jj-1,jk)   &
                     &          + vn(ji,jj  ,jk)**2 * e1v(ji,jj  ) * e3v_n(ji,jj  ,jk)   )
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'crsfld', z3d, 'T', 1. )
         !
         CALL crs_dom_ope( z3d, 'VOL', 'T', tmask, zt_crs, p_e12=e1e2t, p_e3=ze3t, psgn=1.0 )
         CALL iom_put( "eken", zt_crs )
      ENDIF
      !  Horizontal divergence ( following OCE/DYN/divhor.F90 ) 
      DO jk = 1, jpkm1
         DO ji = 2, jpi_crsm1
            DO jj = 2, jpj_crsm1
               IF( tmask_crs(ji,jj,jk ) > 0 ) THEN
                   z2dcrsu =  ( un_crs(ji  ,jj  ,jk) * crs_surfu_wgt(ji  ,jj  ,jk) ) &
                      &     - ( un_crs(ji-1,jj  ,jk) * crs_surfu_wgt(ji-1,jj  ,jk) )
                   z2dcrsv =  ( vn_crs(ji  ,jj  ,jk) * crs_surfv_wgt(ji  ,jj  ,jk) ) &
                      &     - ( vn_crs(ji  ,jj-1,jk) * crs_surfv_wgt(ji  ,jj-1,jk) )
                   !
                   hdivn_crs(ji,jj,jk) = ( z2dcrsu + z2dcrsv ) / crs_volt_wgt(ji,jj,jk) 
               ENDIF
            END DO
         END DO
      END DO
      CALL crs_lbc_lnk( hdivn_crs, 'T', 1.0 )
      !
      CALL iom_put( "hdiv", hdivn_crs )  


      !  W-velocity
      IF( ln_crs_wn ) THEN
         CALL crs_dom_ope( wn, 'SUM', 'W', tmask, wn_crs, p_e12=e1e2t, p_surf_crs=e1e2w_msk, psgn=1.0 )
       !  CALL crs_dom_ope( wn, 'VOL', 'W', tmask, wn_crs, p_e12=e1e2t, p_e3=ze3w )
      ELSE
        wn_crs(:,:,jpk) = 0._wp
        DO jk = jpkm1, 1, -1
           wn_crs(:,:,jk) = wn_crs(:,:,jk+1) - e3t_crs(:,:,jk) * hdivn_crs(:,:,jk)
        ENDDO
      ENDIF
      CALL iom_put( "woce", wn_crs  )   ! vertical velocity
      !  free memory

      !  avs
      SELECT CASE ( nn_crs_kz )
         CASE ( 0 )
            CALL crs_dom_ope( avt, 'VOL', 'W', tmask, avt_crs, p_e12=e1e2t, p_e3=ze3w, psgn=1.0 )
            CALL crs_dom_ope( avs, 'VOL', 'W', tmask, avs_crs, p_e12=e1e2t, p_e3=ze3w, psgn=1.0 )
         CASE ( 1 )
            CALL crs_dom_ope( avt, 'MAX', 'W', tmask, avt_crs, p_e12=e1e2t, p_e3=ze3w, psgn=1.0 )
            CALL crs_dom_ope( avs, 'MAX', 'W', tmask, avs_crs, p_e12=e1e2t, p_e3=ze3w, psgn=1.0 )
         CASE ( 2 )
            CALL crs_dom_ope( avt, 'MIN', 'W', tmask, avt_crs, p_e12=e1e2t, p_e3=ze3w, psgn=1.0 )
            CALL crs_dom_ope( avs, 'MIN', 'W', tmask, avs_crs, p_e12=e1e2t, p_e3=ze3w, psgn=1.0 )
      END SELECT
      !
      CALL iom_put( "avt", avt_crs )   !  Kz on T
      CALL iom_put( "avs", avs_crs )   !  Kz on S
      
      !  sbc fields  
      CALL crs_dom_ope( sshn , 'VOL', 'T', tmask, sshn_crs , p_e12=e1e2t, p_e3=ze3t           , psgn=1.0 )  
      CALL crs_dom_ope( utau , 'SUM', 'U', umask, utau_crs , p_e12=e2u  , p_surf_crs=e2u_crs  , psgn=1.0 )
      CALL crs_dom_ope( vtau , 'SUM', 'V', vmask, vtau_crs , p_e12=e1v  , p_surf_crs=e1v_crs  , psgn=1.0 )
      CALL crs_dom_ope( wndm , 'SUM', 'T', tmask, wndm_crs , p_e12=e1e2t, p_surf_crs=e1e2t_crs, psgn=1.0 )
      CALL crs_dom_ope( rnf  , 'MAX', 'T', tmask, rnf_crs                                     , psgn=1.0 )
      CALL crs_dom_ope( qsr  , 'SUM', 'T', tmask, qsr_crs  , p_e12=e1e2t, p_surf_crs=e1e2t_crs, psgn=1.0 )
      CALL crs_dom_ope( emp_b, 'SUM', 'T', tmask, emp_b_crs, p_e12=e1e2t, p_surf_crs=e1e2t_crs, psgn=1.0 )
      CALL crs_dom_ope( emp  , 'SUM', 'T', tmask, emp_crs  , p_e12=e1e2t, p_surf_crs=e1e2t_crs, psgn=1.0 )
      CALL crs_dom_ope( sfx  , 'SUM', 'T', tmask, sfx_crs  , p_e12=e1e2t, p_surf_crs=e1e2t_crs, psgn=1.0 )
      CALL crs_dom_ope( fr_i , 'SUM', 'T', tmask, fr_i_crs , p_e12=e1e2t, p_surf_crs=e1e2t_crs, psgn=1.0 )

      CALL iom_put( "ssh"      , sshn_crs )   ! ssh output 
      CALL iom_put( "utau"     , utau_crs )   ! i-tau output 
      CALL iom_put( "vtau"     , vtau_crs )   ! j-tau output 
      CALL iom_put( "wspd"     , wndm_crs )   ! wind speed output 
      CALL iom_put( "runoffs"  , rnf_crs  )   ! runoff output 
      CALL iom_put( "qsr"      , qsr_crs  )   ! qsr output 
      CALL iom_put( "empmr"    , emp_crs  )   ! water flux output 
      CALL iom_put( "saltflx"  , sfx_crs  )   ! salt flux output 
      CALL iom_put( "ice_cover", fr_i_crs )   ! ice cover output 

      !
      CALL iom_swap( "nemo" )     ! return back on high-resolution grid
      !
      IF( ln_timing )   CALL timing_stop('crs_fld')
      !
   END SUBROUTINE crs_fld

   !!======================================================================
END MODULE crsfld
