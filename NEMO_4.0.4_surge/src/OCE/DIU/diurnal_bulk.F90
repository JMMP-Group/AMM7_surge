MODULE diurnal_bulk
   !!======================================================================
   !!                    ***  MODULE  diurnal_bulk  ***
   !!     Takaya model of diurnal warming (Takaya, 2010)
   !!=====================================================================
   !! History :        !  11-10  (J. While)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   diurnal_sst_bulk_init  : initialise diurnal model
   !!   diurnal_sst_bulk_step  : time-step the diurnal model
   !!----------------------------------------------------------------------
   USE par_kind
   USE phycst
   USE dom_oce
   USE lib_mpp
   USE solfrac_mod
   USE in_out_manager
   
   IMPLICIT NONE
   PRIVATE
   
   ! Namelist parameters
   LOGICAL, PUBLIC :: ln_diurnal
   LOGICAL, PUBLIC :: ln_diurnal_only

   ! Parameters
   REAL(wp), PRIVATE, PARAMETER :: pp_alpha = 2.0e-4_wp
   REAL(wp), PRIVATE, PARAMETER :: pp_veltol = 0._wp
   REAL(wp), PRIVATE, PARAMETER :: pp_min_fvel = 1.e-10_wp 
   
   ! Key variables
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: x_dsst     ! Delta SST
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: x_solfrac  ! Fraction of 
   !                                                           ! absorbed radiation

   PUBLIC diurnal_sst_bulk_init, diurnal_sst_takaya_step
   
   !!----------------------------------------------------------------------
CONTAINS 
   
   SUBROUTINE diurnal_sst_bulk_init
      !!----------------------------------------------------------------------
      !! *** ROUTINE diurnal_sst_init ***
      !!
      !! ** Purpose : Initialise the Takaya diurnal model   
      !!----------------------------------------------------------------------      
      INTEGER ::   ios   ! local integer
      !!
      NAMELIST /namdiu/ ln_diurnal, ln_diurnal_only
      !!----------------------------------------------------------------------      

      ! Read the namelist
      REWIND( numnam_ref )
      READ  ( numnam_ref, namdiu, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdiu in reference namelist' )
      REWIND( numnam_cfg )
      READ  ( numnam_cfg, namdiu, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdiu in configuration namelist' )      
      !
      IF( ln_diurnal_only .AND. ( .NOT. ln_diurnal ) ) THEN
         CALL ctl_stop( "ln_diurnal_only set, but ln_diurnal = FALSE !" )
      ENDIF
      
      IF( ln_diurnal ) THEN      
         !         
         ALLOCATE( x_dsst(jpi,jpj), x_solfrac(jpi,jpj) )
         !
         x_solfrac = 0._wp         ! Initialise the solar fraction
         x_dsst    = 0._wp
         !
         IF( ln_diurnal_only ) THEN
            CALL ctl_warn( "ln_diurnal_only set; only the diurnal component of SST will be calculated" )
         ENDIF
      ENDIF
      
   END SUBROUTINE diurnal_sst_bulk_init


   SUBROUTINE diurnal_sst_takaya_step(kt, psolflux, pqflux, ptauflux, prho, p_rdt,   &
            &                  pla, pthick, pcoolthick, pmu, &
            &                  p_fvel_bkginc, p_hflux_bkginc)
      !!----------------------------------------------------------------------
      !! *** ROUTINE diurnal_sst_takaya_step ***
      !!
      !! ** Purpose :   Time-step the Takaya diurnal model
      !!
      !! ** Method :    1) Calculate the Obukhov length
      !!                2) Calculate the Similarity function
      !!                2) Calculate the increment to dsst
      !!                3) Apply the increment
      !! ** Reference : Refinements to a prognostic scheme of skin sea surface
      !!                temperature, Takaya et al, JGR, 2010 
      !!----------------------------------------------------------------------
      INTEGER                               , INTENT(in) ::   kt             ! time step
      REAL(wp), DIMENSION(jpi,jpj)          , INTENT(in) ::   psolflux       ! solar flux (Watts)
      REAL(wp), DIMENSION(jpi,jpj)          , INTENT(in) ::   pqflux         ! heat (non-solar) flux (Watts)
      REAL(wp), DIMENSION(jpi,jpj)          , INTENT(in) ::   ptauflux       ! wind stress  (kg/ m s^2)
      REAL(wp), DIMENSION(jpi,jpj)          , INTENT(in) ::   prho           ! water density  (kg/m^3)
      REAL(wp)                              , INTENT(in) ::   p_rdt          ! time-step
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL, INTENT(in) ::   pLa            ! Langmuir number
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL, INTENT(in) ::   pthick         ! warm layer thickness (m)
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL, INTENT(in) ::   pcoolthick     ! cool skin thickness (m)
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL, INTENT(in) ::   pmu            ! mu parameter
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL, INTENT(in) ::   p_hflux_bkginc ! increment to the heat flux
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL, INTENT(in) ::   p_fvel_bkginc  ! increment to the friction velocity
      !
      INTEGER :: ji,jj
      LOGICAL  :: ll_calcfrac
      REAL(wp), DIMENSION(jpi,jpj) :: z_fvel              ! friction velocity     
      REAL(wp), DIMENSION(jpi,jpj) :: zthick, zcoolthick, zmu, zla
      REAL(wp), DIMENSION(jpi,jpj) :: z_abflux            ! absorbed flux           
      REAL(wp), DIMENSION(jpi,jpj) :: z_fla               ! Langmuir function value 
      !!----------------------------------------------------------------------

      ! Set optional arguments to their defaults
      IF( .NOT. PRESENT( pthick )   ) THEN   ;   zthick(:,:) = 3._wp
      ELSE                                   ;   zthick(:,:) = pthick(:,:)
      ENDIF
      IF( .NOT. PRESENT(pcoolthick) ) THEN   ;   zcoolthick(:,:) = 0._wp
      ELSE                                   ;   zcoolthick(:,:) = pcoolthick(:,:)
      ENDIF
      IF( .NOT. PRESENT( pmu )      ) THEN   ;   zmu(:,:) = 0.3_wp
      ELSE                                   ;   zmu(:,:) = pmu(:,:)
      ENDIF
      IF( .NOT. PRESENT(pla) ) THEN          ;   zla(:,:) = 0.3_wp
      ELSE                                   ;   zla(:,:) = pla(:,:)
      ENDIF
      
      ! If not done already, calculate the solar fraction
      IF ( kt==nit000 ) THEN
         DO jj = 1,jpj
            DO ji = 1, jpi
               IF(  ( x_solfrac(ji,jj) == 0._wp ) .AND. ( tmask(ji,jj,1) == 1._wp ) ) &
                  &   x_solfrac(ji,jj) = solfrac( zcoolthick(ji,jj),zthick(ji,jj) )
            END DO
         END DO   
      ENDIF   

      ! convert solar flux and heat flux to absorbed flux   
      WHERE ( tmask(:,:,1) == 1._wp) 
         z_abflux(:,:) = (  x_solfrac(:,:) *  psolflux (:,:)) + pqflux(:,:)     
      ELSEWHERE
         z_abflux(:,:) = 0._wp
      ENDWHERE
      IF( PRESENT(p_hflux_bkginc) ) z_abflux(:,:) = z_abflux(:,:) + p_hflux_bkginc   ! Optional increment
      WHERE ( ABS( z_abflux(:,:) ) < rsmall )
         z_abflux(:,:) = rsmall
      ENDWHERE 
      
      ! Calculate the friction velocity
      WHERE ( (ptauflux /= 0) .AND. ( tmask(:,:,1) == 1.) )   
         z_fvel(:,:) = SQRT( ptauflux(:,:) / prho(:,:) )
      ELSEWHERE
         z_fvel(:,:) = 0._wp
      ENDWHERE
      IF( PRESENT(p_fvel_bkginc) ) z_fvel(:,:) = z_fvel(:,:) + p_fvel_bkginc   ! Optional increment
      
       
       
      ! Calculate the Langmuir function value
      WHERE ( tmask(:,:,1) == 1.)
         z_fla(:,:) = MAX( 1._wp, zla(:,:)**( -2._wp / 3._wp ) )  
      ELSEWHERE
         z_fla(:,:) = 0._wp
      ENDWHERE     
      
      ! Increment the temperature using the implicit solution
      x_dsst(:,:) = t_imp( x_dsst(:,:), p_rdt, z_abflux(:,:), z_fvel(:,:),   &
         &                       z_fla(:,:), zmu(:,:), zthick(:,:), prho(:,:) )
      !
   END SUBROUTINE diurnal_sst_takaya_step

   
   FUNCTION t_imp(p_dsst, p_rdt, p_abflux, p_fvel, &
                          p_fla, pmu, pthick, prho )
                          
      IMPLICIT NONE
      
      ! Function definition
      REAL(wp), DIMENSION(jpi,jpj) :: t_imp
      ! Dummy variables
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN) :: p_dsst     ! Delta SST
      REAL(wp), INTENT(IN)                     :: p_rdt      ! Time-step
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN) :: p_abflux   ! Heat forcing
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN) :: p_fvel     ! Friction velocity
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN) :: p_fla      ! Langmuir number
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN) :: pmu        ! Structure parameter
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN) :: pthick     ! Layer thickness
      REAL(wp), DIMENSION(jpi,jpj), INTENT(IN) :: prho       ! Water density
   
      ! Local variables
      REAL(wp) :: z_olength          ! Obukhov length
      REAL(wp) :: z_sigma, z_sigma2
      REAL(wp) :: z_term1, z_term2     
      REAL(wp) :: z_stabfunc         ! stability function value
      REAL(wp) :: z_fvel      
      
      CHARACTER(200) :: warn_string
       
      INTEGER :: ji,jj
                     
      DO jj = 1, jpj
         DO ji = 1, jpi   
            
            ! Only calculate outside tmask
            IF ( tmask(ji,jj,1) /= 1._wp ) THEN
               t_imp(ji,jj) = 0._wp
               CYCLE   
            END IF
            
            IF (p_fvel(ji,jj) < pp_min_fvel) THEN
               z_fvel = pp_min_fvel
               WRITE(warn_string,*) "diurnal_sst_takaya step: "&
               &//"friction velocity < minimum\n" &
               &//"Setting friction velocity =",pp_min_fvel
               CALL ctl_warn(warn_string)
               
            ELSE
               z_fvel = p_fvel(ji,jj)
            ENDIF
                 
            ! Calculate the Obukhov length
            IF ( (z_fvel < pp_veltol ) .AND. &
            &    (p_dsst(ji,jj) > 0._wp) ) THEN
               z_olength =  z_fvel  / &
               &     SQRT( p_dsst(ji,jj) * vkarmn * grav * &
               &             pp_alpha / ( 5._wp * pthick(ji,jj) ) )
            ELSE
               z_olength = &
               &   ( prho(ji,jj) * rcp * z_fvel**3._wp ) / &
               &   ( vkarmn * grav * pp_alpha *&
               &     p_abflux(ji,jj) )
            ENDIF
             
            ! Calculate the stability function          
            z_sigma = pthick(ji,jj) / z_olength
            z_sigma2 = z_sigma * z_sigma
  
            IF ( z_sigma >= 0. ) THEN
               z_stabfunc = 1._wp + &
               &  ( ( 5._wp * z_sigma + 4._wp * z_sigma2 ) / &
               &  ( 1._wp + 3._wp * z_sigma + 0.25_wp * &
               &    z_sigma2 ) )
            ELSE
               z_stabfunc = 1._wp / &
               &                   SQRT( 1._wp - 16._wp * z_sigma )
            ENDIF

            ! Calculate the T increment
            z_term1 = ( p_abflux(ji,jj) * ( pmu(ji,jj) + 1._wp)  / &
            & ( pmu(ji,jj) * pthick(ji,jj) * prho(ji,jj) * rcp ) )
            
             
            z_term2 = -( ( pmu(ji,jj) + 1._wp) * &
            &                       ( vkarmn * z_fvel * p_fla(ji,jj) ) / &
            &      ( pthick(ji,jj) * z_stabfunc ) )     
          
            t_imp(ji,jj) = ( p_dsst(ji,jj) + p_rdt * z_term1 ) / &
                           ( 1._wp - p_rdt * z_term2 )

         END DO
      END DO
      
      END FUNCTION t_imp

END MODULE diurnal_bulk
