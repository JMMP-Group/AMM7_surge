MODULE sbcdcy
   !!======================================================================
   !!                    ***  MODULE  sbcdcy  ***
   !! Ocean forcing:  compute the diurnal cycle
   !!======================================================================
   !! History : OPA  !  2005-02  (D. Bernie)  Original code
   !!   NEMO    2.0  !  2006-02  (S. Masson, G. Madec)  adaptation to NEMO
   !!           3.1  !  2009-07  (J.M. Molines)  adaptation to v3.1
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  sbc_dcy : solar flux at kt from daily mean, taking diurnal cycle into account
   !!----------------------------------------------------------------------
   USE oce              ! ocean dynamics and tracers
   USE phycst           ! ocean physics
   USE dom_oce          ! ocean space and time domain
   USE sbc_oce          ! Surface boundary condition: ocean fields
   !
   USE in_out_manager   ! I/O manager
   USE lib_mpp          ! MPP library

   IMPLICIT NONE
   PRIVATE
   
   INTEGER, PUBLIC ::   nday_qsr   !: day when parameters were computed
   
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   raa , rbb  , rcc  , rab     ! diurnal cycle parameters
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   rtmd, rdawn, rdusk, rscal   !    -      -       -
  
   PUBLIC   sbc_dcy        ! routine called by sbc

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

      INTEGER FUNCTION sbc_dcy_alloc()
         !!----------------------------------------------------------------------
         !!                ***  FUNCTION sbc_dcy_alloc  ***
         !!----------------------------------------------------------------------
         ALLOCATE( raa (jpi,jpj) , rbb  (jpi,jpj) , rcc  (jpi,jpj) , rab  (jpi,jpj) ,     &
            &      rtmd(jpi,jpj) , rdawn(jpi,jpj) , rdusk(jpi,jpj) , rscal(jpi,jpj) , STAT=sbc_dcy_alloc )
            !
         CALL mpp_sum ( 'sbcdcy', sbc_dcy_alloc )
         IF( sbc_dcy_alloc /= 0 )   CALL ctl_stop( 'STOP', 'sbc_dcy_alloc: failed to allocate arrays' )
      END FUNCTION sbc_dcy_alloc


   FUNCTION sbc_dcy( pqsrin, l_mask ) RESULT( zqsrout )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_dcy  ***
      !!
      !! ** Purpose : introduce a diurnal cycle of qsr from daily values
      !!
      !! ** Method  : see Appendix A of Bernie et al. 2007.
      !!
      !! ** Action  : redistribute daily QSR on each time step following the diurnal cycle
      !!
      !! reference  : Bernie, DJ, E Guilyardi, G Madec, JM Slingo, and SJ Woolnough, 2007
      !!              Impact of resolving the diurnal cycle in an ocean--atmosphere GCM. 
      !!              Part 1: a diurnally forced OGCM. Climate Dynamics 29:6, 575-590.
      !!----------------------------------------------------------------------
      LOGICAL , OPTIONAL          , INTENT(in) ::   l_mask    ! use the routine for night mask computation
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pqsrin    ! input daily QSR flux 
      REAL(wp), DIMENSION(jpi,jpj)             ::   zqsrout   ! output QSR flux with diurnal cycle
      !!
      INTEGER  ::   ji, jj                                       ! dummy loop indices
      INTEGER, DIMENSION(jpi,jpj) :: imask_night ! night mask
      REAL(wp) ::   ztwopi, zinvtwopi, zconvrad 
      REAL(wp) ::   zlo, zup, zlousd, zupusd
      REAL(wp) ::   zdsws, zdecrad, ztx, zsin, zcos
      REAL(wp) ::   ztmp, ztmp1, ztmp2, ztest
      REAL(wp) ::   ztmpm, ztmpm1, ztmpm2
      !---------------------------statement functions------------------------
      REAL(wp) ::   fintegral, pt1, pt2, paaa, pbbb, pccc        ! dummy statement function arguments
      fintegral( pt1, pt2, paaa, pbbb, pccc ) =                         &
         &   paaa * pt2 + zinvtwopi * pbbb * SIN(pccc + ztwopi * pt2)   &
         & - paaa * pt1 - zinvtwopi * pbbb * SIN(pccc + ztwopi * pt1)
      !!---------------------------------------------------------------------
      !
      ! Initialization
      ! --------------
      ztwopi    = 2._wp * rpi
      zinvtwopi = 1._wp / ztwopi
      zconvrad  = ztwopi / 360._wp

      ! When are we during the day (from 0 to 1)
      zlo = ( REAL(nsec_day, wp) - 0.5_wp * rdt ) / rday
      zup = zlo + ( REAL(nn_fsbc, wp)     * rdt ) / rday
      !                                          
      IF( nday_qsr == -1 ) THEN       ! first time step only  
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'sbc_dcy : introduce diurnal cycle from daily mean qsr'
            WRITE(numout,*) '~~~~~~~'
            WRITE(numout,*)
         ENDIF
         ! allocate sbcdcy arrays
         IF( sbc_dcy_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_dcy_alloc : unable to allocate arrays' )
         ! Compute rcc needed to compute the time integral of the diurnal cycle
         rcc(:,:) = zconvrad * glamt(:,:) - rpi
         ! time of midday
         rtmd(:,:) = 0.5_wp - glamt(:,:) / 360._wp
         rtmd(:,:) = MOD( (rtmd(:,:) + 1._wp) , 1._wp)
      ENDIF

      ! If this is a new day, we have to update the dawn, dusk and scaling function  
      !----------------------
    
      !     2.1 dawn and dusk  

      ! nday is the number of days since the beginning of the current month 
      IF( nday_qsr /= nday ) THEN 
         ! save the day of the year and the daily mean of qsr
         nday_qsr = nday 
         ! number of days since the previous winter solstice (supposed to be always 21 December)         
         zdsws = REAL(11 + nday_year, wp)
         ! declination of the earths orbit
         zdecrad = (-23.5_wp * zconvrad) * COS( zdsws * ztwopi / REAL(nyear_len(1),wp) )
         ! Compute A and B needed to compute the time integral of the diurnal cycle

         zsin = SIN( zdecrad )   ;   zcos = COS( zdecrad )
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztmp = zconvrad * gphit(ji,jj)
               raa(ji,jj) = SIN( ztmp ) * zsin
               rbb(ji,jj) = COS( ztmp ) * zcos
            END DO  
         END DO  
         ! Compute the time of dawn and dusk

         ! rab to test if the day time is equal to 0, less than 24h of full day        
         rab(:,:) = -raa(:,:) / rbb(:,:)
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF ( ABS(rab(ji,jj)) < 1._wp ) THEN         ! day duration is less than 24h
         ! When is it night?
                  ztx = zinvtwopi * (ACOS(rab(ji,jj)) - rcc(ji,jj))
                  ztest = -rbb(ji,jj) * SIN( rcc(ji,jj) + ztwopi * ztx )
         ! is it dawn or dusk?
                  IF ( ztest > 0._wp ) THEN
                     rdawn(ji,jj) = ztx
                     rdusk(ji,jj) = rtmd(ji,jj) + ( rtmd(ji,jj) - rdawn(ji,jj) )
                  ELSE
                     rdusk(ji,jj) = ztx
                     rdawn(ji,jj) = rtmd(ji,jj) - ( rdusk(ji,jj) - rtmd(ji,jj) )
                  ENDIF
               ELSE
                  rdawn(ji,jj) = rtmd(ji,jj) + 0.5_wp
                  rdusk(ji,jj) = rdawn(ji,jj)
               ENDIF
             END DO  
         END DO  
         rdawn(:,:) = MOD( (rdawn(:,:) + 1._wp), 1._wp )
         rdusk(:,:) = MOD( (rdusk(:,:) + 1._wp), 1._wp )
         !     2.2 Compute the scaling function:
         !         S* = the inverse of the time integral of the diurnal cycle from dawn to dusk
         !         Avoid possible infinite scaling factor, associated with very short daylight
         !         periods, by ignoring periods less than 1/1000th of a day (ticket #1040)
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF ( ABS(rab(ji,jj)) < 1._wp ) THEN         ! day duration is less than 24h
                  rscal(ji,jj) = 0.0_wp
                  IF ( rdawn(ji,jj) < rdusk(ji,jj) ) THEN      ! day time in one part
                     IF( (rdusk(ji,jj) - rdawn(ji,jj) ) .ge. 0.001_wp ) THEN
                       rscal(ji,jj) = fintegral(rdawn(ji,jj), rdusk(ji,jj), raa(ji,jj), rbb(ji,jj), rcc(ji,jj)) 
                       rscal(ji,jj) = 1._wp / rscal(ji,jj)
                     ENDIF
                  ELSE                                         ! day time in two parts
                     IF( (rdusk(ji,jj) + (1._wp - rdawn(ji,jj)) ) .ge. 0.001_wp ) THEN
                       rscal(ji,jj) = fintegral(0._wp, rdusk(ji,jj), raa(ji,jj), rbb(ji,jj), rcc(ji,jj))   &
                          &         + fintegral(rdawn(ji,jj), 1._wp, raa(ji,jj), rbb(ji,jj), rcc(ji,jj)) 
                       rscal(ji,jj) = 1. / rscal(ji,jj)
                     ENDIF
                  ENDIF
               ELSE
                  IF ( raa(ji,jj) > rbb(ji,jj) ) THEN         ! 24h day
                     rscal(ji,jj) = fintegral(0._wp, 1._wp, raa(ji,jj), rbb(ji,jj), rcc(ji,jj)) 
                     rscal(ji,jj) = 1._wp / rscal(ji,jj)
                  ELSE                                          ! No day
                     rscal(ji,jj) = 0.0_wp
                  ENDIF
               ENDIF
            END DO  
         END DO  
         !
         ztmp = rday / ( rdt * REAL(nn_fsbc, wp) )
         rscal(:,:) = rscal(:,:) * ztmp
         !
      ENDIF 
         !     3. update qsr with the diurnal cycle
         !     ------------------------------------

      imask_night(:,:) = 0
      DO jj = 1, jpj
         DO ji = 1, jpi
            ztmpm = 0._wp
            IF( ABS(rab(ji,jj)) < 1. ) THEN         ! day duration is less than 24h
               !
               IF( rdawn(ji,jj) < rdusk(ji,jj) ) THEN       ! day time in one part
                  zlousd = MAX(zlo, rdawn(ji,jj))
                  zlousd = MIN(zlousd, zup)
                  zupusd = MIN(zup, rdusk(ji,jj))
                  zupusd = MAX(zupusd, zlo)
                  ztmp = fintegral(zlousd, zupusd, raa(ji,jj), rbb(ji,jj), rcc(ji,jj)) 
                  zqsrout(ji,jj) = pqsrin(ji,jj) * ztmp * rscal(ji,jj)
                  ztmpm = zupusd - zlousd
                  IF ( ztmpm .EQ. 0 ) imask_night(ji,jj) = 1
                  !
               ELSE                                         ! day time in two parts
                  zlousd = MIN(zlo, rdusk(ji,jj))
                  zupusd = MIN(zup, rdusk(ji,jj))
                  ztmp1 = fintegral(zlousd, zupusd, raa(ji,jj), rbb(ji,jj), rcc(ji,jj)) 
                  ztmpm1=zupusd-zlousd
                  zlousd = MAX(zlo, rdawn(ji,jj))
                  zupusd = MAX(zup, rdawn(ji,jj))
                  ztmp2 = fintegral(zlousd, zupusd, raa(ji,jj), rbb(ji,jj), rcc(ji,jj)) 
                  ztmpm2 =zupusd-zlousd
                  ztmp = ztmp1 + ztmp2
                  ztmpm = ztmpm1 + ztmpm2
                  zqsrout(ji,jj) = pqsrin(ji,jj) * ztmp * rscal(ji,jj)
                  IF (ztmpm .EQ. 0.) imask_night(ji,jj) = 1
               ENDIF
            ELSE                                   ! 24h light or 24h night
               !
               IF( raa(ji,jj) > rbb(ji,jj) ) THEN           ! 24h day
                  ztmp = fintegral(zlo, zup, raa(ji,jj), rbb(ji,jj), rcc(ji,jj)) 
                  zqsrout(ji,jj) = pqsrin(ji,jj) * ztmp * rscal(ji,jj)
                  imask_night(ji,jj) = 0
                  !
               ELSE                                         ! No day
                  zqsrout(ji,jj) = 0.0_wp
                  imask_night(ji,jj) = 1
               ENDIF
            ENDIF
         END DO  
      END DO  
      !
      IF( PRESENT(l_mask) .AND. l_mask ) THEN
         zqsrout(:,:) = float(imask_night(:,:))
      ENDIF
      !
   END FUNCTION sbc_dcy

   !!======================================================================
END MODULE sbcdcy
