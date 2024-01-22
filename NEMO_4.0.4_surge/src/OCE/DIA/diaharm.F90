MODULE diaharm 
   !!======================================================================
   !!                       ***  MODULE  diaharm  ***
   !! Harmonic analysis of tidal constituents 
   !!======================================================================
   !! History :  3.1  !  2007  (O. Le Galloudec, J. Chanut)  Original code
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE phycst
   USE daymod
   USE tide_mod
   USE sbctide         ! Tidal forcing or not
   !
   USE in_out_manager  ! I/O units
   USE iom             ! I/0 library
   USE ioipsl          ! NetCDF IPSL library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE timing          ! preformance summary
   USE lib_mpp           ! MPP library

   IMPLICIT NONE
   PRIVATE
   
   INTEGER, PARAMETER :: jpincomax    = 2.*jpmax_harmo
   INTEGER, PARAMETER :: jpdimsparse  = jpincomax*366*24*2 ! 30min for a 1yr-long run

   !                         !!** namelist variables **
   LOGICAL, PUBLIC ::   ln_diaharm    ! Choose tidal harmonic output or not
   INTEGER         ::   nit000_han    ! First time step used for harmonic analysis
   INTEGER         ::   nitend_han    ! Last time step used for harmonic analysis
   INTEGER         ::   nstep_han     ! Time step frequency for harmonic analysis
   INTEGER         ::   nb_ana        ! Number of harmonics to analyse

   INTEGER , ALLOCATABLE, DIMENSION(:)       ::   name
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   ana_temp
   REAL(wp), ALLOCATABLE, DIMENSION(:)       ::   ana_freq, ut, vt, ft
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)   ::   out_eta, out_u, out_v

   INTEGER  ::   ninco, nsparse
   REAL(wp) ::   z1_tmp3
   INTEGER ,       DIMENSION(jpdimsparse)         ::   njsparse, nisparse
   INTEGER , SAVE, DIMENSION(jpincomax)           ::   ipos1
   REAL(wp),       DIMENSION(jpdimsparse)         ::   valuesparse
   REAL(wp),       DIMENSION(jpincomax)           ::   ztmp4 , ztmp7, z1_pivot
   REAL(wp), SAVE, DIMENSION(jpincomax,jpincomax) ::   ztmp3 , zpilier

   CHARACTER (LEN=4), DIMENSION(jpmax_harmo) ::   tname   ! Names of tidal constituents ('M2', 'K1',...)

   PUBLIC   dia_harm        ! routine called by step.F90
   PUBLIC   dia_harm_init   ! routine called by nemogcm.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_harm_init 
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dia_harm_init  ***
      !!         
      !! ** Purpose :   Initialization of tidal harmonic analysis
      !!
      !! ** Method  :   Initialize frequency array and  nodal factor for nit000_han
      !!
      !!--------------------------------------------------------------------
      INTEGER ::   jh, nhan, ji
      INTEGER ::   ios                 ! Local integer output status for namelist read

      NAMELIST/nam_diaharm/ ln_diaharm, nit000_han, nitend_han, nstep_han, tname
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_harm_init: Tidal harmonic analysis initialization'
         WRITE(numout,*) '~~~~~~~~~~~~~ '
      ENDIF
      !
      REWIND( numnam_ref )              ! Namelist nam_diaharm in reference namelist : Tidal harmonic analysis
      READ  ( numnam_ref, nam_diaharm, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nam_diaharm in reference namelist' )
      REWIND( numnam_cfg )              ! Namelist nam_diaharm in configuration namelist : Tidal harmonic analysis
      READ  ( numnam_cfg, nam_diaharm, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nam_diaharm in configuration namelist' )
      IF(lwm) WRITE ( numond, nam_diaharm )
      !
      IF(lwp) THEN
         WRITE(numout,*) 'Tidal diagnostics = ', ln_diaharm
         WRITE(numout,*) '   First time step used for analysis:         nit000_han= ', nit000_han
         WRITE(numout,*) '   Last  time step used for analysis:         nitend_han= ', nitend_han
         WRITE(numout,*) '   Time step frequency for harmonic analysis: nstep_han = ', nstep_han
      ENDIF

      IF( ln_diaharm .AND. .NOT.ln_tide )   CALL ctl_stop( 'dia_harm_init : ln_tide must be true for harmonic analysis')

      IF( ln_diaharm ) THEN

         CALL tide_init_Wave
         !
         ! Basic checks on harmonic analysis time window:
         ! ----------------------------------------------
         IF( nit000 > nit000_han )   CALL ctl_stop( 'dia_harm_init : nit000_han must be greater than nit000',   &
            &                                       ' restart capability not implemented' )
         IF( nitend < nitend_han )   CALL ctl_stop( 'dia_harm_init : nitend_han must be lower than nitend',   &
            &                                       'restart capability not implemented' )

         IF( MOD( nitend_han-nit000_han+1 , nstep_han ) /= 0 )   &
            &                        CALL ctl_stop( 'dia_harm_init : analysis time span must be a multiple of nstep_han' )
         !
         nb_ana = 0
         DO jh=1,jpmax_harmo
            DO ji=1,jpmax_harmo
               IF(TRIM(tname(jh)) == Wave(ji)%cname_tide) THEN
                  nb_ana=nb_ana+1
               ENDIF
            END DO
         END DO
         !
         IF(lwp) THEN
            WRITE(numout,*) '        Namelist nam_diaharm'
            WRITE(numout,*) '        nb_ana    = ', nb_ana
            CALL flush(numout)
         ENDIF
         !
         IF (nb_ana > jpmax_harmo) THEN
            WRITE(ctmp1,*) ' nb_ana must be lower than jpmax_harmo'
            WRITE(ctmp2,*) ' jpmax_harmo= ', jpmax_harmo
            CALL ctl_stop( 'dia_harm_init', ctmp1, ctmp2 )
         ENDIF

         ALLOCATE(name(nb_ana))
         DO jh=1,nb_ana
            DO ji=1,jpmax_harmo
               IF (TRIM(tname(jh)) ==  Wave(ji)%cname_tide) THEN
                  name(jh) = ji
                  EXIT
               END IF
            END DO
         END DO

         ! Initialize frequency array:
         ! ---------------------------
         ALLOCATE( ana_freq(nb_ana), ut(nb_ana), vt(nb_ana), ft(nb_ana) )

         CALL tide_harmo( ana_freq, vt, ut, ft, name, nb_ana )

         IF(lwp) WRITE(numout,*) 'Analysed frequency  : ',nb_ana ,'Frequency '

         DO jh = 1, nb_ana
            IF(lwp) WRITE(numout,*) '                    : ',tname(jh),' ',ana_freq(jh)
         END DO

         ! Initialize temporary arrays:
         ! ----------------------------
         ALLOCATE( ana_temp(jpi,jpj,2*nb_ana,3) )
         ana_temp(:,:,:,:) = 0._wp

      ENDIF

   END SUBROUTINE dia_harm_init


   SUBROUTINE dia_harm( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dia_harm  ***
      !!         
      !! ** Purpose :   Tidal harmonic analysis main routine
      !!
      !! ** Action  :   Sums ssh/u/v over time analysis [nit000_han,nitend_han]
      !!
      !!--------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt
      !
      INTEGER  ::   ji, jj, jh, jc, nhc
      REAL(wp) ::   ztime, ztemp
      !!--------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('dia_harm')
      !
      IF( kt >= nit000_han .AND. kt <= nitend_han .AND. MOD(kt,nstep_han) == 0 ) THEN
         !
         ztime = (kt-nit000+1) * rdt 
         !
         nhc = 0
         DO jh = 1, nb_ana
            DO jc = 1, 2
               nhc = nhc+1
               ztemp =  (   MOD(jc,2) * ft(jh) *COS(ana_freq(jh)*ztime + vt(jh) + ut(jh))  &
                  &    +(1.-MOD(jc,2))* ft(jh) *SIN(ana_freq(jh)*ztime + vt(jh) + ut(jh)))
                  !
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1
                     ana_temp(ji,jj,nhc,1) = ana_temp(ji,jj,nhc,1) + ztemp * sshn(ji,jj) * ssmask (ji,jj) ! elevation      
                     ana_temp(ji,jj,nhc,2) = ana_temp(ji,jj,nhc,2) + ztemp * un_b(ji,jj) * ssumask(ji,jj) ! u-vel
                     ana_temp(ji,jj,nhc,3) = ana_temp(ji,jj,nhc,3) + ztemp * vn_b(ji,jj) * ssvmask(ji,jj) ! v-vel
                  END DO
               END DO
            END DO
         END DO
      END IF
      !
      IF( kt == nitend_han )   CALL dia_harm_end
      !
      IF( ln_timing )   CALL timing_stop('dia_harm')
      !
   END SUBROUTINE dia_harm


   SUBROUTINE dia_harm_end
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE diaharm_end  ***
      !!         
      !! ** Purpose :  Compute the Real and Imaginary part of tidal constituents
      !!
      !! ** Action  :  Decompose the signal on the harmonic constituents 
      !!
      !!--------------------------------------------------------------------
      INTEGER  ::   ji, jj, jh, jc, jn, nhan
      INTEGER  ::   ksp, kun, keq
      REAL(wp) ::   ztime, ztime_ini, ztime_end, z1_han
      !!--------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dia_harm_end: kt=nitend_han: Perform harmonic analysis'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      
      ALLOCATE( out_eta(jpi,jpj,2*nb_ana), out_u(jpi,jpj,2*nb_ana), out_v(jpi,jpj,2*nb_ana) )

      ztime_ini = nit000_han*rdt                 ! Initial time in seconds at the beginning of analysis
      ztime_end = nitend_han*rdt                 ! Final time in seconds at the end of analysis
      nhan = (nitend_han-nit000_han+1)/nstep_han ! Number of dumps used for analysis
      z1_han = 1._wp / REAL(nhan-1) 
      
      ninco = 2*nb_ana

      ksp = 0
      keq = 0        
      DO jn = 1, nhan
         ztime=( (nhan-jn)*ztime_ini + (jn-1)*ztime_end ) * z1_han
         keq = keq + 1
         kun = 0
         DO jh = 1, nb_ana
            DO jc = 1, 2
               kun = kun + 1
               ksp = ksp + 1
               nisparse(ksp) = keq
               njsparse(ksp) = kun
               valuesparse(ksp) = (   MOD(jc,2) * ft(jh) * COS(ana_freq(jh)*ztime + vt(jh) + ut(jh))   &
                  &             + (1.-MOD(jc,2))* ft(jh) * SIN(ana_freq(jh)*ztime + vt(jh) + ut(jh)) )
            END DO
         END DO
      END DO

      nsparse = ksp

      IF( nsparse > jpdimsparse )   CALL ctl_stop( 'STOP', 'SUR_DETERMINE : nsparse .GT. jpdimsparse')
      IF( ninco   > jpincomax   )   CALL ctl_stop( 'STOP', 'SUR_DETERMINE : ninco .GT. jpincomax')

      CALL SUR_DETERMINE_INIT

      ! Elevation:
      DO jj = 2, jpjm1
         DO ji = 2, jpim1

            ! Fill input array
            ztmp4(1:nb_ana*2) = ana_temp(ji,jj,1:nb_ana*2,1)
            CALL SUR_DETERMINE
            
            ! Fill output array
            DO jh = 1, nb_ana
               out_eta(ji,jj,jh       ) =  ztmp7((jh-1)*2+1) * ssmask(ji,jj)
               out_eta(ji,jj,jh+nb_ana) = -ztmp7((jh-1)*2+2) * ssmask(ji,jj)
            END DO
         END DO
      END DO

      ! ubar:
      DO jj = 2, jpjm1
         DO ji = 2, jpim1

            ! Fill input array
            ztmp4(1:nb_ana*2) = ana_temp(ji,jj,1:nb_ana*2,2)
            CALL SUR_DETERMINE

            ! Fill output array
            DO jh = 1, nb_ana
               out_u(ji,jj,       jh) =  ztmp7((jh-1)*2+1) * ssumask(ji,jj)
               out_u(ji,jj,nb_ana+jh) = -ztmp7((jh-1)*2+2) * ssumask(ji,jj)
            END DO

        END DO
      END DO

      ! vbar:
      DO jj = 2, jpjm1
         DO ji = 2, jpim1

            ! Fill input array
            ztmp4(1:nb_ana*2) = ana_temp(ji,jj,1:nb_ana*2,3)
            CALL SUR_DETERMINE

            ! Fill output array
            DO jh = 1, nb_ana
               out_v(ji,jj,       jh) =  ztmp7((jh-1)*2+1) * ssvmask(ji,jj)
               out_v(ji,jj,nb_ana+jh) = -ztmp7((jh-1)*2+2) * ssvmask(ji,jj)
            END DO

         END DO
      END DO
      !
      ! clem: we could avoid this call if all the loops were from 1:jpi and 1:jpj
      !       but I think this is the most efficient
      CALL lbc_lnk_multi( 'dia_harm_end', out_eta, 'T', 1., out_u, 'U', -1. , out_v, 'V', -1. )
      !
      CALL dia_wri_harm ! Write results in files
      !
      DEALLOCATE( out_eta, out_u, out_v )
      !
   END SUBROUTINE dia_harm_end


   SUBROUTINE dia_wri_harm
      !!--------------------------------------------------------------------
      !!                 ***  ROUTINE dia_wri_harm  ***
      !!         
      !! ** Purpose : Write tidal harmonic analysis results in a netcdf file
      !!--------------------------------------------------------------------
      INTEGER  ::   jh
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*) '  '
      IF(lwp) WRITE(numout,*) 'dia_wri_harm : Write harmonic analysis results'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'

      ! A) Elevation
      !/////////////
      DO jh = 1, nb_ana
      CALL iom_put( TRIM(tname(jh))//'x', out_eta(:,:,jh) )
      CALL iom_put( TRIM(tname(jh))//'y', out_eta(:,:,jh+nb_ana) )
      END DO

      ! B) ubar
      !/////////
      DO jh = 1, nb_ana
      CALL iom_put( TRIM(tname(jh))//'x_u', out_u(:,:,jh) )
      CALL iom_put( TRIM(tname(jh))//'y_u', out_u(:,:,jh+nb_ana) )
      END DO

      ! C) vbar
      !/////////
      DO jh = 1, nb_ana
         CALL iom_put( TRIM(tname(jh))//'x_v', out_v(:,:,jh       ) )
         CALL iom_put( TRIM(tname(jh))//'y_v', out_v(:,:,jh+nb_ana) )
      END DO
      !
   END SUBROUTINE dia_wri_harm


   SUBROUTINE SUR_DETERMINE_INIT
      !!---------------------------------------------------------------------------------
      !!                      *** ROUTINE SUR_DETERMINE_INIT ***
      !!       
      !!---------------------------------------------------------------------------------
      INTEGER                        :: ji_sd, jj_sd, ji1_sd, ji2_sd, jh1_sd, jh2_sd
      INTEGER                        :: ipivot
      REAL(wp)                       :: zval1, zval2, zcol1, zcol2
      INTEGER , DIMENSION(jpincomax) :: ipos2
      !!---------------------------------------------------------------------------------
      !            
      !
      ztmp3(:,:) = 0._wp
      !
      DO jh1_sd = 1, nsparse
         DO jh2_sd = 1, nsparse
            IF( nisparse(jh2_sd) == nisparse(jh1_sd) ) THEN
               ztmp3(njsparse(jh1_sd),njsparse(jh2_sd)) = ztmp3(njsparse(jh1_sd),njsparse(jh2_sd))  &
                  &                                     + valuesparse(jh1_sd)*valuesparse(jh2_sd)
            ENDIF
         END DO
      END DO
      !
      DO jj_sd = 1, ninco
         ipos1(jj_sd) = jj_sd
         ipos2(jj_sd) = jj_sd
      END DO
      !
      DO ji_sd = 1, ninco
         !
         !find greatest non-zero pivot:
         zval1 = ABS(ztmp3(ji_sd,ji_sd))
         !
         ipivot = ji_sd
         DO jj_sd = ji_sd, ninco
            zval2 = ABS(ztmp3(ji_sd,jj_sd))
            IF( zval2 >= zval1 )THEN
               ipivot = jj_sd
               zval1  = zval2
            ENDIF
         END DO
         !
         DO ji1_sd = 1, ninco
            zcol1                = ztmp3(ji1_sd,ji_sd)
            zcol2                = ztmp3(ji1_sd,ipivot)
            ztmp3(ji1_sd,ji_sd)  = zcol2
            ztmp3(ji1_sd,ipivot) = zcol1
         END DO
         !
         ipos2(ji_sd)    = ipos1(ipivot)
         ipos2(ipivot)   = ipos1(ji_sd)
         ipos1(ji_sd)    = ipos2(ji_sd)
         ipos1(ipivot)   = ipos2(ipivot)
         z1_pivot(ji_sd) = 1._wp / ztmp3(ji_sd,ji_sd)
         DO jj_sd = 1, ninco
            ztmp3(ji_sd,jj_sd) = ztmp3(ji_sd,jj_sd) * z1_pivot(ji_sd)
         END DO
         !
         DO ji2_sd = ji_sd+1, ninco
            zpilier(ji2_sd,ji_sd) = ztmp3(ji2_sd,ji_sd)
            DO jj_sd=1,ninco
               ztmp3(ji2_sd,jj_sd) = ztmp3(ji2_sd,jj_sd) - ztmp3(ji_sd,jj_sd) * zpilier(ji2_sd,ji_sd)
            END DO
         END DO
         !
      END DO
      !
      z1_tmp3 = 1._wp / ztmp3(ninco,ninco)
      !
   END SUBROUTINE SUR_DETERMINE_INIT

   
   SUBROUTINE SUR_DETERMINE
      !!---------------------------------------------------------------------------------
      !!                      *** ROUTINE SUR_DETERMINE ***
      !!    
      !!---------------------------------------------------------------------------------
      INTEGER                        :: ji_sd, jj_sd, ji1_sd, ji2_sd
      REAL(wp)                       :: zx1
      REAL(wp), DIMENSION(jpincomax) :: ztmpx
      !!---------------------------------------------------------------------------------
      !            
      DO ji_sd = 1, ninco
         ztmp4(ji_sd) = ztmp4(ji_sd) * z1_pivot(ji_sd)
         DO ji2_sd = ji_sd+1, ninco
            ztmp4(ji2_sd) = ztmp4(ji2_sd) - ztmp4(ji_sd) * zpilier(ji2_sd,ji_sd)
         END DO
      END DO

      !system solving: 
      ztmpx(ninco) = ztmp4(ninco) * z1_tmp3
      DO ji_sd = ninco-1, 1, -1
         zx1 = 0._wp
         DO jj_sd = ji_sd+1, ninco
            zx1 = zx1 + ztmpx(jj_sd) * ztmp3(ji_sd,jj_sd)
         END DO
         ztmpx(ji_sd) = ztmp4(ji_sd) - zx1
      END DO

      DO jj_sd = 1, ninco
         ztmp7(ipos1(jj_sd)) = ztmpx(jj_sd)
      END DO
      !
   END SUBROUTINE SUR_DETERMINE

   
   !!======================================================================
END MODULE diaharm
