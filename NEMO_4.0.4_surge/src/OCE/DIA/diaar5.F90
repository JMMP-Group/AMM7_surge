MODULE diaar5
   !!======================================================================
   !!                       ***  MODULE  diaar5  ***
   !! AR5 diagnostics
   !!======================================================================
   !! History :  3.2  !  2009-11  (S. Masson)  Original code
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase + merge TRC-TRA
   !!----------------------------------------------------------------------
   !!   dia_ar5       : AR5 diagnostics
   !!   dia_ar5_init  : initialisation of AR5 diagnostics
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers 
   USE dom_oce        ! ocean space and time domain
   USE eosbn2         ! equation of state                (eos_bn2 routine)
   USE phycst         ! physical constant
   USE in_out_manager  ! I/O manager
   USE zdfddm
   USE zdf_oce
   !
   USE lib_mpp        ! distribued memory computing library
   USE iom            ! I/O manager library
   USE fldread        ! type FLD_N
   USE timing         ! preformance summary

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_ar5        ! routine called in step.F90 module
   PUBLIC   dia_ar5_alloc  ! routine called in nemogcm.F90 module
   PUBLIC   dia_ar5_hst    ! heat/salt transport

   REAL(wp)                         ::   vol0         ! ocean volume (interior domain)
   REAL(wp)                         ::   area_tot     ! total ocean surface (interior domain)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::   thick0       ! ocean thickness (interior domain)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sn0          ! initial salinity

   LOGICAL  :: l_ar5
      
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION dia_ar5_alloc()
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_ar5_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: dia_ar5_alloc
      !!----------------------------------------------------------------------
      !
      ALLOCATE( thick0(jpi,jpj) , sn0(jpi,jpj,jpk) , STAT=dia_ar5_alloc )
      !
      CALL mpp_sum ( 'diaar5', dia_ar5_alloc )
      IF( dia_ar5_alloc /= 0 )   CALL ctl_stop( 'STOP', 'dia_ar5_alloc: failed to allocate arrays' )
      !
   END FUNCTION dia_ar5_alloc


   SUBROUTINE dia_ar5( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_ar5  ***
      !!
      !! ** Purpose :   compute and output some AR5 diagnostics
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk, iks, ikb                      ! dummy loop arguments
      REAL(wp) ::   zvolssh, zvol, zssh_steric, zztmp, zarho, ztemp, zsal, zmass, zsst
      REAL(wp) ::   zaw, zbw, zrw
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)     :: zarea_ssh , zbotpres       ! 2D workspace 
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)     :: zpe, z2d                   ! 2D workspace 
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   :: zrhd , ztpot               ! 3D workspace
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: ztsn                       ! 4D workspace

      !!--------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('dia_ar5')
 
      IF( kt == nit000 )     CALL dia_ar5_init

      IF( l_ar5 ) THEN 
         ALLOCATE( zarea_ssh(jpi,jpj), zbotpres(jpi,jpj), z2d(jpi,jpj) )
         ALLOCATE( zrhd(jpi,jpj,jpk) )
         ALLOCATE( ztsn(jpi,jpj,jpk,jpts) )
         zarea_ssh(:,:) = e1e2t(:,:) * sshn(:,:)
      ENDIF
      !
      CALL iom_put( 'e2u'      , e2u  (:,:) )
      CALL iom_put( 'e1v'      , e1v  (:,:) )
      CALL iom_put( 'areacello', e1e2t(:,:) )
      !
      IF( iom_use( 'volcello' ) .OR. iom_use( 'masscello' )  ) THEN  
         zrhd(:,:,jpk) = 0._wp        ! ocean volume ; rhd is used as workspace
         DO jk = 1, jpkm1
            zrhd(:,:,jk) = e1e2t(:,:) * e3t_n(:,:,jk) * tmask(:,:,jk)
         END DO
         CALL iom_put( 'volcello'  , zrhd(:,:,:)  )  ! WARNING not consistent with CMIP DR where volcello is at ca. 2000
         CALL iom_put( 'masscello' , rau0 * e3t_n(:,:,:) * tmask(:,:,:) )  ! ocean mass
      ENDIF 
      !
      IF( iom_use( 'e3tb' ) )  THEN    ! bottom layer thickness
         DO jj = 1, jpj
            DO ji = 1, jpi
               ikb = mbkt(ji,jj)
               z2d(ji,jj) = e3t_n(ji,jj,ikb)
            END DO
         END DO
         CALL iom_put( 'e3tb', z2d )
      ENDIF 
      !
      IF( iom_use( 'voltot' ) .OR. iom_use( 'sshtot' )  .OR. iom_use( 'sshdyn' )  ) THEN    
         !                                         ! total volume of liquid seawater
         zvolssh = glob_sum( 'diaar5', zarea_ssh(:,:) ) 
         zvol    = vol0 + zvolssh
      
         CALL iom_put( 'voltot', zvol               )
         CALL iom_put( 'sshtot', zvolssh / area_tot )
         CALL iom_put( 'sshdyn', sshn(:,:) - (zvolssh / area_tot) )
         !
      ENDIF

      IF( iom_use( 'botpres' ) .OR. iom_use( 'sshthster' )  .OR. iom_use( 'sshsteric' )  ) THEN    
         !                     
         ztsn(:,:,:,jp_tem) = tsn(:,:,:,jp_tem)                    ! thermosteric ssh
         ztsn(:,:,:,jp_sal) = sn0(:,:,:)
         CALL eos( ztsn, zrhd, gdept_n(:,:,:) )                       ! now in situ density using initial salinity
         !
         zbotpres(:,:) = 0._wp                        ! no atmospheric surface pressure, levitating sea-ice
         DO jk = 1, jpkm1
            zbotpres(:,:) = zbotpres(:,:) + e3t_n(:,:,jk) * zrhd(:,:,jk)
         END DO
         IF( ln_linssh ) THEN
            IF( ln_isfcav ) THEN
               DO ji = 1, jpi
                  DO jj = 1, jpj
                     iks = mikt(ji,jj)
                     zbotpres(ji,jj) = zbotpres(ji,jj) + sshn(ji,jj) * zrhd(ji,jj,iks) + riceload(ji,jj)
                  END DO
               END DO
            ELSE
               zbotpres(:,:) = zbotpres(:,:) + sshn(:,:) * zrhd(:,:,1)
            END IF
!!gm
!!gm   riceload should be added in both ln_linssh=T or F, no?
!!gm
         END IF
         !                                         
         zarho = glob_sum( 'diaar5', e1e2t(:,:) * zbotpres(:,:) ) 
         zssh_steric = - zarho / area_tot
         CALL iom_put( 'sshthster', zssh_steric )
      
         !                                         ! steric sea surface height
         zbotpres(:,:) = 0._wp                        ! no atmospheric surface pressure, levitating sea-ice
         DO jk = 1, jpkm1
            zbotpres(:,:) = zbotpres(:,:) + e3t_n(:,:,jk) * rhd(:,:,jk)
         END DO
         IF( ln_linssh ) THEN
            IF ( ln_isfcav ) THEN
               DO ji = 1,jpi
                  DO jj = 1,jpj
                     iks = mikt(ji,jj)
                     zbotpres(ji,jj) = zbotpres(ji,jj) + sshn(ji,jj) * rhd(ji,jj,iks) + riceload(ji,jj)
                  END DO
               END DO
            ELSE
               zbotpres(:,:) = zbotpres(:,:) + sshn(:,:) * rhd(:,:,1)
            END IF
         END IF
         !    
         zarho = glob_sum( 'diaar5', e1e2t(:,:) * zbotpres(:,:) ) 
         zssh_steric = - zarho / area_tot
         CALL iom_put( 'sshsteric', zssh_steric )
         !                                         ! ocean bottom pressure
         zztmp = rau0 * grav * 1.e-4_wp               ! recover pressure from pressure anomaly and cover to dbar = 1.e4 Pa
         zbotpres(:,:) = zztmp * ( zbotpres(:,:) + sshn(:,:) + thick0(:,:) )
         CALL iom_put( 'botpres', zbotpres )
         !
      ENDIF

      IF( iom_use( 'masstot' ) .OR. iom_use( 'temptot' )  .OR. iom_use( 'saltot' )  ) THEN    
          !                                         ! Mean density anomalie, temperature and salinity
          ztsn(:,:,:,:) = 0._wp                    ! ztsn(:,:,1,jp_tem/sal) is used here as 2D Workspace for temperature & salinity
          DO jk = 1, jpkm1
             DO jj = 1, jpj
                DO ji = 1, jpi
                   zztmp = e1e2t(ji,jj) * e3t_n(ji,jj,jk)
                   ztsn(ji,jj,1,jp_tem) = ztsn(ji,jj,1,jp_tem) + zztmp * tsn(ji,jj,jk,jp_tem)
                   ztsn(ji,jj,1,jp_sal) = ztsn(ji,jj,1,jp_sal) + zztmp * tsn(ji,jj,jk,jp_sal)
                ENDDO
             ENDDO
          ENDDO

          IF( ln_linssh ) THEN
            IF( ln_isfcav ) THEN
               DO ji = 1, jpi
                  DO jj = 1, jpj
                     iks = mikt(ji,jj)
                     ztsn(ji,jj,1,jp_tem) = ztsn(ji,jj,1,jp_tem) + zarea_ssh(ji,jj) * tsn(ji,jj,iks,jp_tem) 
                     ztsn(ji,jj,1,jp_sal) = ztsn(ji,jj,1,jp_sal) + zarea_ssh(ji,jj) * tsn(ji,jj,iks,jp_sal) 
                  END DO
               END DO
            ELSE
               ztsn(:,:,1,jp_tem) = ztsn(:,:,1,jp_tem) + zarea_ssh(:,:) * tsn(:,:,1,jp_tem) 
               ztsn(:,:,1,jp_sal) = ztsn(:,:,1,jp_sal) + zarea_ssh(:,:) * tsn(:,:,1,jp_sal) 
            END IF
         ENDIF
         !
         ztemp = glob_sum( 'diaar5', ztsn(:,:,1,jp_tem) )
         zsal  = glob_sum( 'diaar5', ztsn(:,:,1,jp_sal) )
         zmass = rau0 * ( zarho + zvol )      
         !
         CALL iom_put( 'masstot', zmass )
         CALL iom_put( 'temptot', ztemp / zvol )
         CALL iom_put( 'saltot' , zsal  / zvol )
         !
      ENDIF     

      IF( ln_teos10 ) THEN        ! ! potential temperature (TEOS-10 case)
         IF( iom_use( 'toce_pot') .OR. iom_use( 'temptot_pot' ) .OR. iom_use( 'sst_pot' )  &
                                  .OR. iom_use( 'ssttot' ) .OR.  iom_use( 'tosmint_pot' ) ) THEN
            !
            ALLOCATE( ztpot(jpi,jpj,jpk) )
            ztpot(:,:,jpk) = 0._wp
            DO jk = 1, jpkm1
               ztpot(:,:,jk) = eos_pt_from_ct( tsn(:,:,jk,jp_tem), tsn(:,:,jk,jp_sal) )
            END DO
            !
            CALL iom_put( 'toce_pot', ztpot(:,:,:) )  ! potential temperature (TEOS-10 case)
            CALL iom_put( 'sst_pot' , ztpot(:,:,1) )  ! surface temperature
            !
            IF( iom_use( 'temptot_pot' ) ) THEN   ! Output potential temperature in case we use TEOS-10
               z2d(:,:) = 0._wp
               DO jk = 1, jpkm1
                 z2d(:,:) = z2d(:,:) + e1e2t(:,:) * e3t_n(:,:,jk) * ztpot(:,:,jk)
               END DO
               ztemp = glob_sum( 'diaar5', z2d(:,:)  ) 
               CALL iom_put( 'temptot_pot', ztemp / zvol )
             ENDIF
             !
             IF( iom_use( 'ssttot' ) ) THEN   ! Output potential temperature in case we use TEOS-10
               zsst = glob_sum( 'diaar5',  e1e2t(:,:) * ztpot(:,:,1)  ) 
               CALL iom_put( 'ssttot', zsst / area_tot )
             ENDIF
             ! Vertical integral of temperature
             IF( iom_use( 'tosmint_pot') ) THEN
               z2d(:,:) = 0._wp
               DO jk = 1, jpkm1
                  DO jj = 1, jpj
                     DO ji = 1, jpi   ! vector opt.
                        z2d(ji,jj) = z2d(ji,jj) + rau0 * e3t_n(ji,jj,jk) *  ztpot(ji,jj,jk)
                     END DO
                  END DO
               END DO
               CALL iom_put( 'tosmint_pot', z2d ) 
            ENDIF
            DEALLOCATE( ztpot )
        ENDIF
      ELSE       
         IF( iom_use('ssttot') ) THEN   ! Output sst in case we use EOS-80
            zsst  = glob_sum( 'diaar5', e1e2t(:,:) * tsn(:,:,1,jp_tem) )
            CALL iom_put('ssttot', zsst / area_tot )
         ENDIF
      ENDIF

      IF( iom_use( 'tnpeo' )) THEN    
        ! Work done against stratification by vertical mixing
        ! Exclude points where rn2 is negative as convection kicks in here and
        ! work is not being done against stratification
         ALLOCATE( zpe(jpi,jpj) )
         zpe(:,:) = 0._wp
         IF( ln_zdfddm ) THEN
            DO jk = 2, jpk
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     IF( rn2(ji,jj,jk) > 0._wp ) THEN
                        zrw = ( gdept_n(ji,jj,jk) - gdepw_n(ji,jj,jk) ) / e3w_n(ji,jj,jk)
                        !
                        zaw = rab_n(ji,jj,jk,jp_tem) * (1. - zrw) + rab_n(ji,jj,jk-1,jp_tem)* zrw
                        zbw = rab_n(ji,jj,jk,jp_sal) * (1. - zrw) + rab_n(ji,jj,jk-1,jp_sal)* zrw
                        !
                        zpe(ji, jj) = zpe(ji,jj)   &
                           &        -  grav * (  avt(ji,jj,jk) * zaw * (tsn(ji,jj,jk-1,jp_tem) - tsn(ji,jj,jk,jp_tem) )  &
                           &                   - avs(ji,jj,jk) * zbw * (tsn(ji,jj,jk-1,jp_sal) - tsn(ji,jj,jk,jp_sal) ) )
                     ENDIF
                  END DO
               END DO
             END DO
          ELSE
            DO jk = 1, jpk
               DO ji = 1, jpi
                  DO jj = 1, jpj
                     zpe(ji,jj) = zpe(ji,jj) + avt(ji,jj,jk) * MIN(0._wp,rn2(ji,jj,jk)) * rau0 * e3w_n(ji,jj,jk)
                  END DO
               END DO
            END DO
         ENDIF
          CALL iom_put( 'tnpeo', zpe )
          DEALLOCATE( zpe )
      ENDIF

      IF( l_ar5 ) THEN
        DEALLOCATE( zarea_ssh , zbotpres, z2d )
        DEALLOCATE( zrhd      )
        DEALLOCATE( ztsn      )
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('dia_ar5')
      !
   END SUBROUTINE dia_ar5


   SUBROUTINE dia_ar5_hst( ktra, cptr, pua, pva ) 
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_ar5_htr ***
      !!----------------------------------------------------------------------
      !! Wrapper for heat transport calculations
      !! Called from all advection and/or diffusion routines
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in )  :: ktra  ! tracer index
      CHARACTER(len=3)                , INTENT(in)   :: cptr  ! transport type  'adv'/'ldf'
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)   :: pua   ! 3D input array of advection/diffusion
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)   :: pva   ! 3D input array of advection/diffusion
      !
      INTEGER    ::  ji, jj, jk
      REAL(wp), DIMENSION(jpi,jpj)  :: z2d

    
      z2d(:,:) = pua(:,:,1) 
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               z2d(ji,jj) = z2d(ji,jj) + pua(ji,jj,jk) 
            END DO
         END DO
       END DO
       CALL lbc_lnk( 'diaar5', z2d, 'U', -1. )
       IF( cptr == 'adv' ) THEN
          IF( ktra == jp_tem ) CALL iom_put( 'uadv_heattr' , rau0_rcp * z2d )  ! advective heat transport in i-direction
          IF( ktra == jp_sal ) CALL iom_put( 'uadv_salttr' , rau0     * z2d )  ! advective salt transport in i-direction
       ENDIF
       IF( cptr == 'ldf' ) THEN
          IF( ktra == jp_tem ) CALL iom_put( 'udiff_heattr' , rau0_rcp * z2d ) ! diffusive heat transport in i-direction
          IF( ktra == jp_sal ) CALL iom_put( 'udiff_salttr' , rau0     * z2d ) ! diffusive salt transport in i-direction
       ENDIF
       !
       z2d(:,:) = pva(:,:,1) 
       DO jk = 1, jpkm1
          DO jj = 2, jpjm1
             DO ji = fs_2, fs_jpim1   ! vector opt.
                z2d(ji,jj) = z2d(ji,jj) + pva(ji,jj,jk) 
             END DO
          END DO
       END DO
       CALL lbc_lnk( 'diaar5', z2d, 'V', -1. )
       IF( cptr == 'adv' ) THEN
          IF( ktra == jp_tem ) CALL iom_put( 'vadv_heattr' , rau0_rcp * z2d )  ! advective heat transport in j-direction
          IF( ktra == jp_sal ) CALL iom_put( 'vadv_salttr' , rau0     * z2d )  ! advective salt transport in j-direction
       ENDIF
       IF( cptr == 'ldf' ) THEN
          IF( ktra == jp_tem ) CALL iom_put( 'vdiff_heattr' , rau0_rcp * z2d ) ! diffusive heat transport in j-direction
          IF( ktra == jp_sal ) CALL iom_put( 'vdiff_salttr' , rau0     * z2d ) ! diffusive salt transport in j-direction
       ENDIF
          
   END SUBROUTINE dia_ar5_hst


   SUBROUTINE dia_ar5_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_ar5_init  ***
      !!                   
      !! ** Purpose :   initialization for AR5 diagnostic computation
      !!----------------------------------------------------------------------
      INTEGER  ::   inum
      INTEGER  ::   ik, idep
      INTEGER  ::   ji, jj, jk  ! dummy loop indices
      REAL(wp) ::   zztmp  
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   zsaldta   ! Jan/Dec levitus salinity
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)     ::   zvol0     
      !
      !!----------------------------------------------------------------------
      !
      l_ar5 = .FALSE.
      IF(   iom_use( 'voltot'  ) .OR. iom_use( 'sshtot'    )  .OR. iom_use( 'sshdyn' )  .OR.  & 
         &  iom_use( 'masstot' ) .OR. iom_use( 'temptot'   )  .OR. iom_use( 'saltot' ) .OR.  &    
         &  iom_use( 'botpres' ) .OR. iom_use( 'sshthster' )  .OR. iom_use( 'sshsteric' ) .OR. &
         &  iom_use( 'rhop' )  ) L_ar5 = .TRUE.
  
      IF( l_ar5 ) THEN
         !
         !                                      ! allocate dia_ar5 arrays
         IF( dia_ar5_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dia_ar5_init : unable to allocate arrays' )

         area_tot  = glob_sum( 'diaar5', e1e2t(:,:) )

         ALLOCATE( zvol0(jpi,jpj) )
         zvol0 (:,:) = 0._wp
         thick0(:,:) = 0._wp
         DO jk = 1, jpkm1
            DO jj = 1, jpj               ! interpolation of salinity at the last ocean level (i.e. the partial step)
               DO ji = 1, jpi
                  idep = tmask(ji,jj,jk) * e3t_0(ji,jj,jk)
                  zvol0 (ji,jj) = zvol0 (ji,jj) +  idep * e1e2t(ji,jj)
                  thick0(ji,jj) = thick0(ji,jj) +  idep    
               END DO
            END DO
         END DO
         vol0 = glob_sum( 'diaar5', zvol0 )
         DEALLOCATE( zvol0 )

         IF( iom_use( 'sshthster' ) ) THEN
            ALLOCATE( zsaldta(jpi,jpj,jpk,jpts) )
            CALL iom_open ( 'sali_ref_clim_monthly', inum )
            CALL iom_get  ( inum, jpdom_data, 'vosaline' , zsaldta(:,:,:,1), 1  )
            CALL iom_get  ( inum, jpdom_data, 'vosaline' , zsaldta(:,:,:,2), 12 )
            CALL iom_close( inum )

            sn0(:,:,:) = 0.5_wp * ( zsaldta(:,:,:,1) + zsaldta(:,:,:,2) )        
            sn0(:,:,:) = sn0(:,:,:) * tmask(:,:,:)
            IF( ln_zps ) THEN               ! z-coord. partial steps
               DO jj = 1, jpj               ! interpolation of salinity at the last ocean level (i.e. the partial step)
                  DO ji = 1, jpi
                     ik = mbkt(ji,jj)
                     IF( ik > 1 ) THEN
                        zztmp = ( gdept_1d(ik) - gdept_0(ji,jj,ik) ) / ( gdept_1d(ik) - gdept_1d(ik-1) )
                        sn0(ji,jj,ik) = ( 1._wp - zztmp ) * sn0(ji,jj,ik) + zztmp * sn0(ji,jj,ik-1)
                     ENDIF
                  END DO
               END DO
            ENDIF
            !
            DEALLOCATE( zsaldta )
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE dia_ar5_init

   !!======================================================================
END MODULE diaar5
