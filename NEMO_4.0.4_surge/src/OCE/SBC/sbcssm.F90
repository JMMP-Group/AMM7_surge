MODULE sbcssm
   !!======================================================================
   !!                       ***  MODULE  sbcssm  ***
   !! Surface module :  provide time-mean ocean surface variables
   !!======================================================================
   !! History :  9.0  ! 2006-07  (G. Madec)  Original code
   !!            3.3  ! 2010-10  (C. Bricaud, G. Madec)  add the Patm forcing for sea-ice
   !!            3.7  ! 2015-11  (G. Madec)  non linear free surface by default: e3t_m always computed
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_ssm       : calculate sea surface mean currents, temperature,  
   !!                   and salinity over nn_fsbc time-step
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! surface boundary condition: ocean fields
   USE sbcapr         ! surface boundary condition: atmospheric pressure
   USE eosbn2         ! equation of state and related derivatives
   USE traqsr, ONLY: ln_traqsr
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE iom            ! IOM library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ssm        ! routine called by step.F90
   PUBLIC   sbc_ssm_init   ! routine called by sbcmod.F90

   LOGICAL, SAVE ::   l_ssm_mean = .FALSE.   ! keep track of whether means have been read from restart file
   
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_ssm( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_oce  ***
      !!                     
      !! ** Purpose :   provide ocean surface variable to sea-surface boundary
      !!                condition computation 
      !!                
      !! ** Method  :   compute mean surface velocity (2 components at U and 
      !!      V-points) [m/s], temperature [Celsius] and salinity [psu] over
      !!      the periode (kt - nn_fsbc) to kt
      !!         Note that the inverse barometer ssh (i.e. ssh associated with Patm)
      !!      is add to ssh_m when ln_apr_dyn = T. Required for sea-ice dynamics.
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      INTEGER  ::   ji, jj               ! loop index
      REAL(wp) ::   zcoef, zf_sbc       ! local scalar
      REAL(wp), DIMENSION(jpi,jpj,jpts) :: zts
      CHARACTER(len=4),SAVE :: stype
      !!---------------------------------------------------------------------
      IF( kt == nit000 ) THEN
         IF( ln_TEOS10 ) THEN
            stype='abs'   ! teos-10: using absolute salinity (sst is converted to potential temperature for the surface module)
         ELSE IF( ln_EOS80  ) THEN
            stype='pra'   ! eos-80: using practical salinity
         ELSE IF ( ln_SEOS) THEN
            stype='seos' ! seos using Simplified Equation of state (sst is converted to potential temperature for the surface module)
         ENDIF
      ENDIF
      !
      !                                        !* surface T-, U-, V- ocean level variables (T, S, depth, velocity)
      DO jj = 1, jpj
         DO ji = 1, jpi
            zts(ji,jj,jp_tem) = tsn(ji,jj,mikt(ji,jj),jp_tem)
            zts(ji,jj,jp_sal) = tsn(ji,jj,mikt(ji,jj),jp_sal)
         END DO
      END DO
      !
      IF( nn_fsbc == 1 ) THEN                             !   Instantaneous surface fields        !
         !                                                ! ---------------------------------------- !
         ssu_m(:,:) = ub(:,:,1)
         ssv_m(:,:) = vb(:,:,1)
         IF( l_useCT )  THEN    ;   sst_m(:,:) = eos_pt_from_ct( zts(:,:,jp_tem), zts(:,:,jp_sal) )
         ELSE                    ;   sst_m(:,:) = zts(:,:,jp_tem)
         ENDIF
         sss_m(:,:) = zts(:,:,jp_sal)
         !                          ! removed inverse barometer ssh when Patm forcing is used (for sea-ice dynamics)
         IF( ln_apr_dyn ) THEN   ;   ssh_m(:,:) = sshn(:,:) - 0.5 * ( ssh_ib(:,:) + ssh_ibb(:,:) )
         ELSE                    ;   ssh_m(:,:) = sshn(:,:)
         ENDIF
         !
         e3t_m(:,:) = e3t_n(:,:,1)
         !
         frq_m(:,:) = fraqsr_1lev(:,:)
         !
      ELSE
         !                                                ! ----------------------------------------------- !
         IF( kt == nit000 .AND. .NOT. l_ssm_mean ) THEN   !   Initialisation: 1st time-step, no input means !
            !                                             ! ----------------------------------------------- !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'sbc_ssm : mean fields initialised to instantaneous values'
            IF(lwp) WRITE(numout,*) '~~~~~~~   '
            zcoef = REAL( nn_fsbc - 1, wp )
            ssu_m(:,:) = zcoef * ub(:,:,1)
            ssv_m(:,:) = zcoef * vb(:,:,1)
            IF( l_useCT )  THEN    ;   sst_m(:,:) = zcoef * eos_pt_from_ct( zts(:,:,jp_tem), zts(:,:,jp_sal) )
            ELSE                    ;   sst_m(:,:) = zcoef * zts(:,:,jp_tem)
            ENDIF
            sss_m(:,:) = zcoef * zts(:,:,jp_sal)
            !                          ! removed inverse barometer ssh when Patm forcing is used (for sea-ice dynamics)
            IF( ln_apr_dyn ) THEN   ;   ssh_m(:,:) = zcoef * ( sshn(:,:) - 0.5 * ( ssh_ib(:,:) + ssh_ibb(:,:) ) )
            ELSE                    ;   ssh_m(:,:) = zcoef * sshn(:,:)
            ENDIF
            !
            e3t_m(:,:) = zcoef * e3t_n(:,:,1)
            !
            frq_m(:,:) = zcoef * fraqsr_1lev(:,:)
            !                                             ! ---------------------------------------- !
         ELSEIF( MOD( kt - 2 , nn_fsbc ) == 0 ) THEN      !   Initialisation: New mean computation   !
            !                                             ! ---------------------------------------- !
            ssu_m(:,:) = 0._wp     ! reset to zero ocean mean sbc fields
            ssv_m(:,:) = 0._wp
            sst_m(:,:) = 0._wp
            sss_m(:,:) = 0._wp
            ssh_m(:,:) = 0._wp
            e3t_m(:,:) = 0._wp
            frq_m(:,:) = 0._wp
         ENDIF
         !                                                ! ---------------------------------------- !
         !                                                !        Cumulate at each time step        !
         !                                                ! ---------------------------------------- !
         ssu_m(:,:) = ssu_m(:,:) + ub(:,:,1)
         ssv_m(:,:) = ssv_m(:,:) + vb(:,:,1)
         IF( l_useCT )  THEN     ;   sst_m(:,:) = sst_m(:,:) + eos_pt_from_ct( zts(:,:,jp_tem), zts(:,:,jp_sal) )
         ELSE                    ;   sst_m(:,:) = sst_m(:,:) + zts(:,:,jp_tem)
         ENDIF
         sss_m(:,:) = sss_m(:,:) + zts(:,:,jp_sal)
         !                          ! removed inverse barometer ssh when Patm forcing is used (for sea-ice dynamics)
         IF( ln_apr_dyn ) THEN   ;   ssh_m(:,:) = ssh_m(:,:) + sshn(:,:) - 0.5 * ( ssh_ib(:,:) + ssh_ibb(:,:) )
         ELSE                    ;   ssh_m(:,:) = ssh_m(:,:) + sshn(:,:)
         ENDIF
         !
         e3t_m(:,:) = e3t_m(:,:) + e3t_n(:,:,1)
         !
         frq_m(:,:) = frq_m(:,:) + fraqsr_1lev(:,:)

         !                                                ! ---------------------------------------- !
         IF( MOD( kt - 1 , nn_fsbc ) == 0 ) THEN          !   Mean value at each nn_fsbc time-step   !
            !                                             ! ---------------------------------------- !
            zcoef = 1. / REAL( nn_fsbc, wp )
            sst_m(:,:) = sst_m(:,:) * zcoef     ! mean SST             [Celsius]
            sss_m(:,:) = sss_m(:,:) * zcoef     ! mean SSS             [psu]
            ssu_m(:,:) = ssu_m(:,:) * zcoef     ! mean suface current  [m/s]
            ssv_m(:,:) = ssv_m(:,:) * zcoef     !
            ssh_m(:,:) = ssh_m(:,:) * zcoef     ! mean SSH             [m]
            e3t_m(:,:) = e3t_m(:,:) * zcoef     ! mean vertical scale factor [m]
            frq_m(:,:) = frq_m(:,:) * zcoef     ! mean fraction of solar net radiation absorbed in the 1st T level [-]
            !
         ENDIF
         !                                                ! ---------------------------------------- !
         IF( lrst_oce ) THEN                              !      Write in the ocean restart file     !
            !                                             ! ---------------------------------------- !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'sbc_ssm : sea surface mean fields written in ocean restart file ',   &
               &                    'at it= ', kt,' date= ', ndastp
            IF(lwp) WRITE(numout,*) '~~~~~~~'
            zf_sbc = REAL( nn_fsbc, wp )
            IF( lwxios ) CALL iom_swap(      cwxios_context          )
            CALL iom_rstput( kt, nitrst, numrow, 'nn_fsbc', zf_sbc, ldxios = lwxios )    ! sbc frequency
            CALL iom_rstput( kt, nitrst, numrow, 'ssu_m'  , ssu_m, ldxios = lwxios  )    ! sea surface mean fields
            CALL iom_rstput( kt, nitrst, numrow, 'ssv_m'  , ssv_m, ldxios = lwxios  )
            CALL iom_rstput( kt, nitrst, numrow, 'sst_m'  , sst_m, ldxios = lwxios  )
            CALL iom_rstput( kt, nitrst, numrow, 'sss_m'  , sss_m, ldxios = lwxios  )
            CALL iom_rstput( kt, nitrst, numrow, 'ssh_m'  , ssh_m, ldxios = lwxios  )
            CALL iom_rstput( kt, nitrst, numrow, 'e3t_m'  , e3t_m, ldxios = lwxios  )
            CALL iom_rstput( kt, nitrst, numrow, 'frq_m'  , frq_m, ldxios = lwxios  )
            !
            IF( lwxios ) CALL iom_swap(      cxios_context          )
         ENDIF
         !
      ENDIF
      !
      IF( MOD( kt - 1 , nn_fsbc ) == 0 ) THEN          !   Mean value at each nn_fsbc time-step   !
         CALL iom_put( 'ssu_m', ssu_m )
         CALL iom_put( 'ssv_m', ssv_m )
         CALL iom_put( 'sst_m_pot', sst_m )
         CALL iom_put( 'sss_m_'//stype, sss_m )
         CALL iom_put( 'ssh_m', ssh_m )
         CALL iom_put( 'e3t_m', e3t_m )
         CALL iom_put( 'frq_m', frq_m )
      ENDIF
      !
   END SUBROUTINE sbc_ssm


   SUBROUTINE sbc_ssm_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_ssm_init  ***
      !!
      !! ** Purpose :   Initialisation of the sbc data
      !!
      !! ** Action  : - read parameters
      !!----------------------------------------------------------------------
      REAL(wp) ::   zcoef, zf_sbc   ! local scalar
      !!----------------------------------------------------------------------
      !
      IF( nn_fsbc == 1 ) THEN
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbc_ssm_init : sea surface mean fields, nn_fsbc=1 : instantaneous values'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
         !
      ELSE
         !               
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'sbc_ssm_init : sea surface mean fields'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~ '
         !
         IF( ln_rstart .AND. iom_varid( numror, 'nn_fsbc', ldstop = .FALSE. ) > 0 ) THEN
            l_ssm_mean = .TRUE.
            CALL iom_get( numror               , 'nn_fsbc', zf_sbc, ldxios = lrxios )    ! sbc frequency of previous run
            CALL iom_get( numror, jpdom_autoglo, 'ssu_m'  , ssu_m, ldxios = lrxios  )    ! sea surface mean velocity    (U-point)
            CALL iom_get( numror, jpdom_autoglo, 'ssv_m'  , ssv_m, ldxios = lrxios  )    !   "         "    velocity    (V-point)
            CALL iom_get( numror, jpdom_autoglo, 'sst_m'  , sst_m, ldxios = lrxios  )    !   "         "    temperature (T-point)
            CALL iom_get( numror, jpdom_autoglo, 'sss_m'  , sss_m, ldxios = lrxios  )    !   "         "    salinity    (T-point)
            CALL iom_get( numror, jpdom_autoglo, 'ssh_m'  , ssh_m, ldxios = lrxios  )    !   "         "    height      (T-point)
            CALL iom_get( numror, jpdom_autoglo, 'e3t_m'  , e3t_m, ldxios = lrxios  )    ! 1st level thickness          (T-point)
            ! fraction of solar net radiation absorbed in 1st T level
            IF( iom_varid( numror, 'frq_m', ldstop = .FALSE. ) > 0 ) THEN
               CALL iom_get( numror, jpdom_autoglo, 'frq_m'  , frq_m, ldxios = lrxios  )
            ELSE
               frq_m(:,:) = 1._wp   ! default definition
            ENDIF
            !
            IF( zf_sbc /= REAL( nn_fsbc, wp ) ) THEN      ! nn_fsbc has changed between 2 runs
               IF(lwp) WRITE(numout,*) '   restart with a change in the frequency of mean from ', zf_sbc, ' to ', nn_fsbc 
               zcoef = REAL( nn_fsbc - 1, wp ) / zf_sbc 
               ssu_m(:,:) = zcoef * ssu_m(:,:) 
               ssv_m(:,:) = zcoef * ssv_m(:,:)
               sst_m(:,:) = zcoef * sst_m(:,:)
               sss_m(:,:) = zcoef * sss_m(:,:)
               ssh_m(:,:) = zcoef * ssh_m(:,:)
               e3t_m(:,:) = zcoef * e3t_m(:,:)
               frq_m(:,:) = zcoef * frq_m(:,:)
            ELSE
               IF(lwp) WRITE(numout,*) '   mean fields read in the ocean restart file'
            ENDIF
         ENDIF
      ENDIF
      !
      IF( .NOT.l_ssm_mean ) THEN   ! default initialisation. needed by iceistate
         !
         IF(lwp) WRITE(numout,*) '   default initialisation of ss._m arrays'
         ssu_m(:,:) = ub(:,:,1)
         ssv_m(:,:) = vb(:,:,1)
         IF( l_useCT )  THEN    ;   sst_m(:,:) = eos_pt_from_ct( tsn(:,:,1,jp_tem), tsn(:,:,1,jp_sal) )
         ELSE                   ;   sst_m(:,:) = tsn(:,:,1,jp_tem)
         ENDIF
         sss_m(:,:) = tsn  (:,:,1,jp_sal)
         ssh_m(:,:) = sshn (:,:)
         e3t_m(:,:) = e3t_n(:,:,1)
         frq_m(:,:) = 1._wp
         !
      ENDIF
      !
      IF( .NOT. ln_traqsr )   fraqsr_1lev(:,:) = 1._wp   ! default definition: qsr 100% in the fisrt level 
      !
      IF( lwxios.AND.nn_fsbc > 1 ) THEN
         CALL iom_set_rstw_var_active('nn_fsbc')
         CALL iom_set_rstw_var_active('ssu_m')
         CALL iom_set_rstw_var_active('ssv_m')
         CALL iom_set_rstw_var_active('sst_m')
         CALL iom_set_rstw_var_active('sss_m')
         CALL iom_set_rstw_var_active('ssh_m')
         CALL iom_set_rstw_var_active('e3t_m')
         CALL iom_set_rstw_var_active('frq_m')
      ENDIF

   END SUBROUTINE sbc_ssm_init

   !!======================================================================
END MODULE sbcssm
