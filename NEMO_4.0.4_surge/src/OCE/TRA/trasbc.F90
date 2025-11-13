MODULE trasbc
   !!==============================================================================
   !!                       ***  MODULE  trasbc  ***
   !! Ocean active tracers:  surface boundary condition
   !!==============================================================================
   !! History :  OPA  !  1998-10  (G. Madec, G. Roullet, M. Imbard)  Original code
   !!            8.2  !  2001-02  (D. Ludicone)  sea ice and free surface
   !!  NEMO      1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!             -   !  2010-09  (C. Ethe, G. Madec) Merge TRA-TRC
   !!            3.6  !  2014-11  (P. Mathiot) isf melting forcing 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_sbc       : update the tracer trend at ocean surface
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE sbc_oce        ! surface boundary condition: ocean
   USE dom_oce        ! ocean space domain variables
   USE phycst         ! physical constant
   USE eosbn2         ! Equation Of State
   USE sbcmod         ! ln_rnf  
   USE sbcrnf         ! River runoff  
   USE sbcisf         ! Ice shelf   
   USE iscplini       ! Ice sheet coupling
   USE traqsr         ! solar radiation penetration
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers 
#if defined key_asminc   
   USE asminc         ! Assimilation increment
#endif
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE iom            ! xIOS server
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_sbc   ! routine called by step.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_sbc ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_sbc  ***
      !!                   
      !! ** Purpose :   Compute the tracer surface boundary condition trend of
      !!      (flux through the interface, concentration/dilution effect)
      !!      and add it to the general trend of tracer equations.
      !!
      !! ** Method :   The (air+ice)-sea flux has two components: 
      !!      (1) Fext, external forcing (i.e. flux through the (air+ice)-sea interface); 
      !!      (2) Fwe , tracer carried with the water that is exchanged with air+ice. 
      !!               The input forcing fields (emp, rnf, sfx, isf) contain Fext+Fwe,
      !!             they are simply added to the tracer trend (tsa).
      !!               In linear free surface case (ln_linssh=T), the volume of the
      !!             ocean does not change with the water exchanges at the (air+ice)-sea
      !!             interface. Therefore another term has to be added, to mimic the
      !!             concentration/dilution effect associated with water exchanges.
      !!
      !! ** Action  : - Update tsa with the surface boundary condition trend 
      !!              - send trends to trdtra module for further diagnostics(l_trdtra=T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk, jn              ! dummy loop indices  
      INTEGER  ::   ikt, ikb                    ! local integers
      REAL(wp) ::   zfact, z1_e3t, zdep, ztim   ! local scalar
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  ztrdt, ztrds
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tra_sbc')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_sbc : TRAcer Surface Boundary Condition'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !
      IF( l_trdtra ) THEN                    !* Save ta and sa trends
         ALLOCATE( ztrdt(jpi,jpj,jpk) , ztrds(jpi,jpj,jpk) ) 
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal)
      ENDIF
      !
!!gm  This should be moved into sbcmod.F90 module ? (especially now that ln_traqsr is read in namsbc namelist)
      IF( .NOT.ln_traqsr ) THEN     ! no solar radiation penetration
         qns(:,:) = qns(:,:) + qsr(:,:)      ! total heat flux in qns
         qsr(:,:) = 0._wp                     ! qsr set to zero
      ENDIF

      !----------------------------------------
      !        EMP, SFX and QNS effects
      !----------------------------------------
      !                             !==  Set before sbc tracer content fields  ==!
      IF( kt == nit000 ) THEN             !* 1st time-step
         IF( ln_rstart .AND.    &               ! Restart: read in restart file
              & iom_varid( numror, 'sbc_hc_b', ldstop = .FALSE. ) > 0 ) THEN
            IF(lwp) WRITE(numout,*) '          nit000-1 sbc tracer content field read in the restart file'
            zfact = 0.5_wp
            sbc_tsc(:,:,:) = 0._wp
            CALL iom_get( numror, jpdom_autoglo, 'sbc_hc_b', sbc_tsc_b(:,:,jp_tem), ldxios = lrxios )   ! before heat content sbc trend
            CALL iom_get( numror, jpdom_autoglo, 'sbc_sc_b', sbc_tsc_b(:,:,jp_sal), ldxios = lrxios )   ! before salt content sbc trend
         ELSE                                   ! No restart or restart not found: Euler forward time stepping
            zfact = 1._wp
            sbc_tsc(:,:,:) = 0._wp
            sbc_tsc_b(:,:,:) = 0._wp
         ENDIF
      ELSE                                !* other time-steps: swap of forcing fields
         zfact = 0.5_wp
         sbc_tsc_b(:,:,:) = sbc_tsc(:,:,:)
      ENDIF
      !                             !==  Now sbc tracer content fields  ==!
      DO jj = 2, jpj
         DO ji = fs_2, fs_jpim1   ! vector opt.
            sbc_tsc(ji,jj,jp_tem) = r1_rau0_rcp * qns(ji,jj)   ! non solar heat flux
            sbc_tsc(ji,jj,jp_sal) = r1_rau0     * sfx(ji,jj)   ! salt flux due to freezing/melting
         END DO
      END DO
      IF( ln_linssh ) THEN                !* linear free surface  
         DO jj = 2, jpj                         !==>> add concentration/dilution effect due to constant volume cell
            DO ji = fs_2, fs_jpim1   ! vector opt.
               sbc_tsc(ji,jj,jp_tem) = sbc_tsc(ji,jj,jp_tem) + r1_rau0 * emp(ji,jj) * tsn(ji,jj,1,jp_tem)
               sbc_tsc(ji,jj,jp_sal) = sbc_tsc(ji,jj,jp_sal) + r1_rau0 * emp(ji,jj) * tsn(ji,jj,1,jp_sal)
            END DO
         END DO                                 !==>> output c./d. term
         IF( iom_use('emp_x_sst') )   CALL iom_put( "emp_x_sst", emp (:,:) * tsn(:,:,1,jp_tem) )
         IF( iom_use('emp_x_sss') )   CALL iom_put( "emp_x_sss", emp (:,:) * tsn(:,:,1,jp_sal) )
      ENDIF
      !
      DO jn = 1, jpts               !==  update tracer trend  ==!
         DO jj = 2, jpj
            DO ji = fs_2, fs_jpim1   ! vector opt.  
               tsa(ji,jj,1,jn) = tsa(ji,jj,1,jn) + zfact * ( sbc_tsc_b(ji,jj,jn) + sbc_tsc(ji,jj,jn) ) / e3t_n(ji,jj,1)
            END DO
         END DO
      END DO
      !                  
      IF( lrst_oce ) THEN           !==  write sbc_tsc in the ocean restart file  ==!
         IF( lwxios ) CALL iom_swap(      cwxios_context          )
         CALL iom_rstput( kt, nitrst, numrow, 'sbc_hc_b', sbc_tsc(:,:,jp_tem), ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'sbc_sc_b', sbc_tsc(:,:,jp_sal), ldxios = lwxios )
         IF( lwxios ) CALL iom_swap(      cxios_context          )
      ENDIF
      !
      !----------------------------------------
      !       Ice Shelf effects (ISF)
      !     tbl treated as in Losh (2008) JGR
      !----------------------------------------
      !
!!gm BUG ?   Why no differences between non-linear and linear free surface ?
!!gm         probably taken into account in r1_hisf_tbl : to be verified
      IF( ln_isf ) THEN
         zfact = 0.5_wp
         DO jj = 2, jpj
            DO ji = fs_2, fs_jpim1
               !
               ikt = misfkt(ji,jj)
               ikb = misfkb(ji,jj)
               !
               ! level fully include in the ice shelf boundary layer
               ! sign - because fwf sign of evapo (rnf sign of precip)
               DO jk = ikt, ikb - 1
               ! compute trend
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem)                                                &
                     &           + zfact * ( risf_tsc_b(ji,jj,jp_tem) + risf_tsc(ji,jj,jp_tem) )             &
                     &           * r1_hisf_tbl(ji,jj)
               END DO
   
               ! level partially include in ice shelf boundary layer 
               ! compute trend
               tsa(ji,jj,ikb,jp_tem) = tsa(ji,jj,ikb,jp_tem)                                                 &
                  &              + zfact * ( risf_tsc_b(ji,jj,jp_tem) + risf_tsc(ji,jj,jp_tem) )             &
                  &              * r1_hisf_tbl(ji,jj) * ralpha(ji,jj)

            END DO
         END DO
      END IF
      !
      !----------------------------------------
      !        River Runoff effects
      !----------------------------------------
      !
      IF( ln_rnf ) THEN         ! input of heat and salt due to river runoff 
         zfact = 0.5_wp
         DO jj = 2, jpj 
            DO ji = fs_2, fs_jpim1
               IF( rnf(ji,jj) /= 0._wp ) THEN
                  zdep = zfact / h_rnf(ji,jj)
                  DO jk = 1, nk_rnf(ji,jj)
                                        tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem)                                 &
                                           &                 +  ( rnf_tsc_b(ji,jj,jp_tem) + rnf_tsc(ji,jj,jp_tem) ) * zdep
                     IF( ln_rnf_sal )   tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal)                                 &
                                           &                 +  ( rnf_tsc_b(ji,jj,jp_sal) + rnf_tsc(ji,jj,jp_sal) ) * zdep 
                  END DO
               ENDIF
            END DO  
         END DO  
      ENDIF

      IF( iom_use('rnf_x_sst') )   CALL iom_put( "rnf_x_sst", rnf*tsn(:,:,1,jp_tem) )   ! runoff term on sst
      IF( iom_use('rnf_x_sss') )   CALL iom_put( "rnf_x_sss", rnf*tsn(:,:,1,jp_sal) )   ! runoff term on sss

#if defined key_asminc
      !
      !----------------------------------------
      !        Assmilation effects
      !----------------------------------------
      !
      IF( ln_sshinc ) THEN         ! input of heat and salt due to assimilation
      	 !
         IF( ln_linssh ) THEN 
            DO jj = 2, jpj 
               DO ji = fs_2, fs_jpim1
                  ztim = ssh_iau(ji,jj) / e3t_n(ji,jj,1)
                  tsa(ji,jj,1,jp_tem) = tsa(ji,jj,1,jp_tem) + tsn(ji,jj,1,jp_tem) * ztim
                  tsa(ji,jj,1,jp_sal) = tsa(ji,jj,1,jp_sal) + tsn(ji,jj,1,jp_sal) * ztim
               END DO
            END DO
         ELSE
            DO jj = 2, jpj 
               DO ji = fs_2, fs_jpim1
                  ztim = ssh_iau(ji,jj) / ( ht_n(ji,jj) + 1. - ssmask(ji, jj) )
                  tsa(ji,jj,:,jp_tem) = tsa(ji,jj,:,jp_tem) + tsn(ji,jj,:,jp_tem) * ztim
                  tsa(ji,jj,:,jp_sal) = tsa(ji,jj,:,jp_sal) + tsn(ji,jj,:,jp_sal) * ztim
               END DO  
            END DO  
         ENDIF
         !
      ENDIF
      !
#endif
      !
      !----------------------------------------
      !        Ice Sheet coupling imbalance correction to have conservation
      !----------------------------------------
      !
      IF( ln_iscpl .AND. ln_hsb) THEN         ! input of heat and salt due to river runoff 
         DO jk = 1,jpk
            DO jj = 2, jpj 
               DO ji = fs_2, fs_jpim1
                  zdep = 1._wp / e3t_n(ji,jj,jk) 
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) - htsc_iscpl(ji,jj,jk,jp_tem) * zdep
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) - htsc_iscpl(ji,jj,jk,jp_sal) * zdep  
               END DO  
            END DO  
         END DO
      ENDIF

      IF( l_trdtra )   THEN                      ! save the horizontal diffusive trends for further diagnostics
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem) - ztrdt(:,:,:)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal) - ztrds(:,:,:)
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_nsr, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_nsr, ztrds )
         DEALLOCATE( ztrdt , ztrds ) 
      ENDIF
      !
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' sbc  - Ta: ', mask1=tmask,   &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( ln_timing )   CALL timing_stop('tra_sbc')
      !
   END SUBROUTINE tra_sbc

   !!======================================================================
END MODULE trasbc
