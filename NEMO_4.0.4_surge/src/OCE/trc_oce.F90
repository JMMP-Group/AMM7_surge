MODULE trc_oce
   !!======================================================================
   !!                      ***  MODULE  trc_oce  ***
   !! Ocean passive tracer  :  share SMS/Ocean variables
   !!======================================================================
   !! History :  1.0  !  2004-03  (C. Ethe)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   trc_oce_rgb   : tabulated attenuation coefficients for RGB light penetration         
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE dom_oce        ! ocean space and time domain
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_oce_rgb        ! routine called by traqsr.F90
   PUBLIC   trc_oce_rgb_read   ! routine called by traqsr.F90
   PUBLIC   trc_oce_ext_lev    ! function called by traqsr.F90 at least
   PUBLIC   trc_oce_alloc      ! function called by nemogcm.F90

   LOGICAL , PUBLIC ::   l_co2cpl  = .false.   !: atmospheric pco2 recieved from oasis
   LOGICAL , PUBLIC ::   l_offline = .false.   !: offline passive tracers flag
   INTEGER , PUBLIC ::   nn_dttrc              !: frequency of step on passive tracers
   REAL(wp), PUBLIC ::   r_si2                 !: largest depth of extinction (blue & 0.01 mg.m-3)  (RGB)
   !
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   etot3     !: light absortion coefficient
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)   ::   oce_co2   !: ocean carbon flux

#if defined key_top 
   !!----------------------------------------------------------------------
   !!   'key_top'                                                 bio-model          
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_top     = .TRUE.   !: TOP model
#else
   !!----------------------------------------------------------------------
   !! Default option                          No bio-model light absorption      
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_top     = .FALSE.   !: TOP model
#endif
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trc_oce_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  trc_oce_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( etot3(jpi,jpj,jpk), oce_co2(jpi,jpj), STAT=trc_oce_alloc )
      IF( trc_oce_alloc /= 0 )   CALL ctl_warn('trc_oce_alloc: failed to allocate etot3 array')
      !
   END FUNCTION trc_oce_alloc


   SUBROUTINE trc_oce_rgb( prgb )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of of the optical scheme
      !!
      !! ** Method  :   Set a look up table for the optical coefficients
      !!                i.e. the attenuation coefficient for R-G-B light 
      !!                tabulated in Chlorophyll class (from JM Andre)
      !!
      !! ** Action  :   prgb(3,61) tabulated R-G-B attenuation coef. 
      !!
      !! Reference  : Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(3,61), INTENT(out) ::   prgb   ! tabulated attenuation coefficient
      !
      INTEGER  ::   jc     ! dummy loop indice
      INTEGER  ::   irgb   ! temporary integer
      REAL(wp) ::   zchl   ! temporary scalar
      REAL(wp), DIMENSION(4,61) ::   zrgb   ! tabulated attenuation coefficient (formerly read in 'kRGB61.txt')
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   trc_oce_rgb : Initialisation of the optical look-up table'
         WRITE(numout,*) '   ~~~~~~~~~~~ '
      ENDIF
      !
      !  Chlorophyll        !     Blue attenuation     !     Green attenuation    !     Red attenuation      !
      zrgb(1, 1) =  0.010   ;   zrgb(2, 1) = 0.01618   ;   zrgb(3, 1) = 0.07464   ;   zrgb(4, 1) = 0.37807
      zrgb(1, 2) =  0.011   ;   zrgb(2, 2) = 0.01654   ;   zrgb(3, 2) = 0.07480   ;   zrgb(4, 2) = 0.37823
      zrgb(1, 3) =  0.013   ;   zrgb(2, 3) = 0.01693   ;   zrgb(3, 3) = 0.07499   ;   zrgb(4, 3) = 0.37840
      zrgb(1, 4) =  0.014   ;   zrgb(2, 4) = 0.01736   ;   zrgb(3, 4) = 0.07518   ;   zrgb(4, 4) = 0.37859
      zrgb(1, 5) =  0.016   ;   zrgb(2, 5) = 0.01782   ;   zrgb(3, 5) = 0.07539   ;   zrgb(4, 5) = 0.37879
      zrgb(1, 6) =  0.018   ;   zrgb(2, 6) = 0.01831   ;   zrgb(3, 6) = 0.07562   ;   zrgb(4, 6) = 0.37900
      zrgb(1, 7) =  0.020   ;   zrgb(2, 7) = 0.01885   ;   zrgb(3, 7) = 0.07586   ;   zrgb(4, 7) = 0.37923
      zrgb(1, 8) =  0.022   ;   zrgb(2, 8) = 0.01943   ;   zrgb(3, 8) = 0.07613   ;   zrgb(4, 8) = 0.37948
      zrgb(1, 9) =  0.025   ;   zrgb(2, 9) = 0.02005   ;   zrgb(3, 9) = 0.07641   ;   zrgb(4, 9) = 0.37976
      zrgb(1,10) =  0.028   ;   zrgb(2,10) = 0.02073   ;   zrgb(3,10) = 0.07672   ;   zrgb(4,10) = 0.38005
      zrgb(1,11) =  0.032   ;   zrgb(2,11) = 0.02146   ;   zrgb(3,11) = 0.07705   ;   zrgb(4,11) = 0.38036
      zrgb(1,12) =  0.035   ;   zrgb(2,12) = 0.02224   ;   zrgb(3,12) = 0.07741   ;   zrgb(4,12) = 0.38070
      zrgb(1,13) =  0.040   ;   zrgb(2,13) = 0.02310   ;   zrgb(3,13) = 0.07780   ;   zrgb(4,13) = 0.38107
      zrgb(1,14) =  0.045   ;   zrgb(2,14) = 0.02402   ;   zrgb(3,14) = 0.07821   ;   zrgb(4,14) = 0.38146
      zrgb(1,15) =  0.050   ;   zrgb(2,15) = 0.02501   ;   zrgb(3,15) = 0.07866   ;   zrgb(4,15) = 0.38189
      zrgb(1,16) =  0.056   ;   zrgb(2,16) = 0.02608   ;   zrgb(3,16) = 0.07914   ;   zrgb(4,16) = 0.38235
      zrgb(1,17) =  0.063   ;   zrgb(2,17) = 0.02724   ;   zrgb(3,17) = 0.07967   ;   zrgb(4,17) = 0.38285
      zrgb(1,18) =  0.071   ;   zrgb(2,18) = 0.02849   ;   zrgb(3,18) = 0.08023   ;   zrgb(4,18) = 0.38338
      zrgb(1,19) =  0.079   ;   zrgb(2,19) = 0.02984   ;   zrgb(3,19) = 0.08083   ;   zrgb(4,19) = 0.38396
      zrgb(1,20) =  0.089   ;   zrgb(2,20) = 0.03131   ;   zrgb(3,20) = 0.08149   ;   zrgb(4,20) = 0.38458
      zrgb(1,21) =  0.100   ;   zrgb(2,21) = 0.03288   ;   zrgb(3,21) = 0.08219   ;   zrgb(4,21) = 0.38526
      zrgb(1,22) =  0.112   ;   zrgb(2,22) = 0.03459   ;   zrgb(3,22) = 0.08295   ;   zrgb(4,22) = 0.38598
      zrgb(1,23) =  0.126   ;   zrgb(2,23) = 0.03643   ;   zrgb(3,23) = 0.08377   ;   zrgb(4,23) = 0.38676
      zrgb(1,24) =  0.141   ;   zrgb(2,24) = 0.03842   ;   zrgb(3,24) = 0.08466   ;   zrgb(4,24) = 0.38761
      zrgb(1,25) =  0.158   ;   zrgb(2,25) = 0.04057   ;   zrgb(3,25) = 0.08561   ;   zrgb(4,25) = 0.38852
      zrgb(1,26) =  0.178   ;   zrgb(2,26) = 0.04289   ;   zrgb(3,26) = 0.08664   ;   zrgb(4,26) = 0.38950
      zrgb(1,27) =  0.200   ;   zrgb(2,27) = 0.04540   ;   zrgb(3,27) = 0.08775   ;   zrgb(4,27) = 0.39056
      zrgb(1,28) =  0.224   ;   zrgb(2,28) = 0.04811   ;   zrgb(3,28) = 0.08894   ;   zrgb(4,28) = 0.39171
      zrgb(1,29) =  0.251   ;   zrgb(2,29) = 0.05103   ;   zrgb(3,29) = 0.09023   ;   zrgb(4,29) = 0.39294
      zrgb(1,30) =  0.282   ;   zrgb(2,30) = 0.05420   ;   zrgb(3,30) = 0.09162   ;   zrgb(4,30) = 0.39428
      zrgb(1,31) =  0.316   ;   zrgb(2,31) = 0.05761   ;   zrgb(3,31) = 0.09312   ;   zrgb(4,31) = 0.39572
      zrgb(1,32) =  0.355   ;   zrgb(2,32) = 0.06130   ;   zrgb(3,32) = 0.09474   ;   zrgb(4,32) = 0.39727
      zrgb(1,33) =  0.398   ;   zrgb(2,33) = 0.06529   ;   zrgb(3,33) = 0.09649   ;   zrgb(4,33) = 0.39894
      zrgb(1,34) =  0.447   ;   zrgb(2,34) = 0.06959   ;   zrgb(3,34) = 0.09837   ;   zrgb(4,34) = 0.40075
      zrgb(1,35) =  0.501   ;   zrgb(2,35) = 0.07424   ;   zrgb(3,35) = 0.10040   ;   zrgb(4,35) = 0.40270
      zrgb(1,36) =  0.562   ;   zrgb(2,36) = 0.07927   ;   zrgb(3,36) = 0.10259   ;   zrgb(4,36) = 0.40480
      zrgb(1,37) =  0.631   ;   zrgb(2,37) = 0.08470   ;   zrgb(3,37) = 0.10495   ;   zrgb(4,37) = 0.40707
      zrgb(1,38) =  0.708   ;   zrgb(2,38) = 0.09056   ;   zrgb(3,38) = 0.10749   ;   zrgb(4,38) = 0.40952
      zrgb(1,39) =  0.794   ;   zrgb(2,39) = 0.09690   ;   zrgb(3,39) = 0.11024   ;   zrgb(4,39) = 0.41216
      zrgb(1,40) =  0.891   ;   zrgb(2,40) = 0.10374   ;   zrgb(3,40) = 0.11320   ;   zrgb(4,40) = 0.41502
      zrgb(1,41) =  1.000   ;   zrgb(2,41) = 0.11114   ;   zrgb(3,41) = 0.11639   ;   zrgb(4,41) = 0.41809
      zrgb(1,42) =  1.122   ;   zrgb(2,42) = 0.11912   ;   zrgb(3,42) = 0.11984   ;   zrgb(4,42) = 0.42142
      zrgb(1,43) =  1.259   ;   zrgb(2,43) = 0.12775   ;   zrgb(3,43) = 0.12356   ;   zrgb(4,43) = 0.42500
      zrgb(1,44) =  1.413   ;   zrgb(2,44) = 0.13707   ;   zrgb(3,44) = 0.12757   ;   zrgb(4,44) = 0.42887
      zrgb(1,45) =  1.585   ;   zrgb(2,45) = 0.14715   ;   zrgb(3,45) = 0.13189   ;   zrgb(4,45) = 0.43304
      zrgb(1,46) =  1.778   ;   zrgb(2,46) = 0.15803   ;   zrgb(3,46) = 0.13655   ;   zrgb(4,46) = 0.43754
      zrgb(1,47) =  1.995   ;   zrgb(2,47) = 0.16978   ;   zrgb(3,47) = 0.14158   ;   zrgb(4,47) = 0.44240
      zrgb(1,48) =  2.239   ;   zrgb(2,48) = 0.18248   ;   zrgb(3,48) = 0.14701   ;   zrgb(4,48) = 0.44765
      zrgb(1,49) =  2.512   ;   zrgb(2,49) = 0.19620   ;   zrgb(3,49) = 0.15286   ;   zrgb(4,49) = 0.45331
      zrgb(1,50) =  2.818   ;   zrgb(2,50) = 0.21102   ;   zrgb(3,50) = 0.15918   ;   zrgb(4,50) = 0.45942
      zrgb(1,51) =  3.162   ;   zrgb(2,51) = 0.22703   ;   zrgb(3,51) = 0.16599   ;   zrgb(4,51) = 0.46601
      zrgb(1,52) =  3.548   ;   zrgb(2,52) = 0.24433   ;   zrgb(3,52) = 0.17334   ;   zrgb(4,52) = 0.47313
      zrgb(1,53) =  3.981   ;   zrgb(2,53) = 0.26301   ;   zrgb(3,53) = 0.18126   ;   zrgb(4,53) = 0.48080
      zrgb(1,54) =  4.467   ;   zrgb(2,54) = 0.28320   ;   zrgb(3,54) = 0.18981   ;   zrgb(4,54) = 0.48909
      zrgb(1,55) =  5.012   ;   zrgb(2,55) = 0.30502   ;   zrgb(3,55) = 0.19903   ;   zrgb(4,55) = 0.49803
      zrgb(1,56) =  5.623   ;   zrgb(2,56) = 0.32858   ;   zrgb(3,56) = 0.20898   ;   zrgb(4,56) = 0.50768
      zrgb(1,57) =  6.310   ;   zrgb(2,57) = 0.35404   ;   zrgb(3,57) = 0.21971   ;   zrgb(4,57) = 0.51810
      zrgb(1,58) =  7.079   ;   zrgb(2,58) = 0.38154   ;   zrgb(3,58) = 0.23129   ;   zrgb(4,58) = 0.52934
      zrgb(1,59) =  7.943   ;   zrgb(2,59) = 0.41125   ;   zrgb(3,59) = 0.24378   ;   zrgb(4,59) = 0.54147
      zrgb(1,60) =  8.912   ;   zrgb(2,60) = 0.44336   ;   zrgb(3,60) = 0.25725   ;   zrgb(4,60) = 0.55457
      zrgb(1,61) = 10.000   ;   zrgb(2,61) = 0.47804   ;   zrgb(3,61) = 0.27178   ;   zrgb(4,61) = 0.56870
      !
      prgb(:,:) = zrgb(2:4,:)
      !
      r_si2 = 1.e0 / zrgb(2, 1)        ! blue with the smallest chlorophyll concentration)
      IF(lwp) WRITE(numout,*) '      RGB longest depth of extinction    r_si2 = ', r_si2
      !
      DO jc = 1, 61                         ! check
         zchl = zrgb(1,jc)
         irgb = NINT( 41 + 20.* LOG10( zchl ) + 1.e-15 )
         IF(lwp .AND. nn_print >= 1 ) WRITE(numout,*) '    jc =', jc, '  Chl = ', zchl, '  irgb = ', irgb
         IF( irgb /= jc ) THEN
            IF(lwp) WRITE(numout,*) '    jc =', jc, '  Chl = ', zchl, '  Chl class = ', irgb
            CALL ctl_stop( 'trc_oce_rgb : inconsistency in Chl tabulated attenuation coeff.' )
         ENDIF
      END DO
      !
   END SUBROUTINE trc_oce_rgb


   SUBROUTINE trc_oce_rgb_read( prgb )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of of the optical scheme
      !!
      !! ** Method  :   read the look up table for the optical coefficients
      !!
      !! ** input   :   xkrgb(61) precomputed array corresponding to the  
      !!                          attenuation coefficient (from JM Andre)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(3,61), INTENT(out) ::   prgb   ! tabulated attenuation coefficient
      !
      INTEGER  ::   jc, jb ! dummy loop indice
      INTEGER  ::   irgb   ! temporary integer
      REAL(wp) ::   zchl   ! temporary scalar
      INTEGER  ::   numlight
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) ' trc_oce_rgb_read : optical look-up table read in kRGB61.txt file'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~'
         WRITE(numout,*) 
      ENDIF
      !
      CALL ctl_opn( numlight, 'kRGB61.txt', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )
      DO jc = 1, 61
         READ(numlight,*) zchl, ( prgb(jb,jc), jb = 1, 3 )
         irgb = NINT( 41 + 20.* LOG10( zchl ) + 1.e-15 )   
         IF(lwp) WRITE(numout,*) '    jc =', jc, '  Chl = ', zchl, '  irgb = ', irgb  
         IF( irgb /= jc ) THEN  
            IF(lwp) WRITE(numout,*) '    jc =', jc, '  Chl = ', zchl, '  Chl class = ', irgb
            CALL ctl_stop( 'trc_oce_rgb_read : inconsistency in Chl tabulated attenuation coeff.' )
         ENDIF
      END DO
      CLOSE( numlight )
      !
      r_si2 = 1.e0 / prgb(1, 1)      ! blue with the smallest chlorophyll concentration)
      IF(lwp) WRITE(numout,*) '      RGB longest depth of extinction    r_si2 = ', r_si2
      !
   END SUBROUTINE trc_oce_rgb_read


   FUNCTION trc_oce_ext_lev( prldex, pqsr_frc ) RESULT( pjl )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE trc_oce_ext_lev  ***
      !!       
      !! ** Purpose :   compute max. level for light penetration
      !!          
      !! ** Method  :   the function provides the level at which irradiance 
      !!                becomes negligible (i.e. = 1.e-15 W/m2) for 3 or 2 bands light
      !!                penetration: I(z) = pqsr_frc * EXP(hext/prldex) = 1.e-15 W/m2
      !!                # prldex is the longest depth of extinction:
      !!                   - prldex = 23 m (2 bands case)
      !!                   - prldex = 62 m (3 bands case: blue waveband & 0.01 mg/m2 for the chlorophyll)
      !!                # pqsr_frc is the fraction of solar radiation which penetrates,
      !!                considering Qsr=240 W/m2 and rn_abs = 0.58:
      !!                   - pqsr_frc = Qsr * (1-rn_abs)   = 1.00e2 W/m2 (2 bands case)
      !!                   - pqsr_frc = Qsr * (1-rn_abs)/3 = 0.33e2 W/m2 (3 bands case & equi-partition)
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::   prldex    ! longest depth of extinction
      REAL(wp), INTENT(in) ::   pqsr_frc  ! frac. solar radiation which penetrates 
      !
      INTEGER  ::   jk, pjl            ! levels
      REAL(wp) ::   zhext              ! deepest level till which light penetrates
      REAL(wp) ::   zprec = 15._wp     ! precision to reach -LOG10(1.e-15)
      REAL(wp) ::   zem                ! temporary scalar 
      !!----------------------------------------------------------------------
      !
      ! It is not necessary to compute anything below the following depth
      zhext = prldex * ( LOG(10._wp) * zprec + LOG(pqsr_frc) )
      !
      ! Level of light extinction
      pjl = jpkm1
      DO jk = jpkm1, 1, -1
         IF(SUM(tmask(:,:,jk)) > 0 ) THEN
            zem = MAXVAL( gdepw_0(:,:,jk+1) * tmask(:,:,jk) )
            IF( zem >= zhext )   pjl = jk                       ! last T-level reached by Qsr
         ELSE
            pjl = jk                                            ! or regional sea-bed depth 
         ENDIF
      END DO
      !
   END FUNCTION trc_oce_ext_lev

   !!======================================================================
END MODULE trc_oce
