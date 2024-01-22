MODULE icbthm
   !!======================================================================
   !!                       ***  MODULE  icbthm  ***
   !! Icebergs:  thermodynamics routines for icebergs
   !!======================================================================
   !! History : 3.3.1 !  2010-01  (Martin&Adcroft) Original code
   !!            -    !  2011-03  (Madec)          Part conversion to NEMO form
   !!            -    !                            Removal of mapping from another grid
   !!            -    !  2011-04  (Alderson)       Split into separate modules
   !!            -    !  2011-05  (Alderson)       Use tmask instead of tmask_i
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   icb_thm : initialise
   !!             reference for equations - M = Martin + Adcroft, OM 34, 2010
   !!----------------------------------------------------------------------
   USE par_oce        ! NEMO parameters
   USE dom_oce        ! NEMO domain
   USE in_out_manager ! NEMO IO routines, numout in particular
   USE lib_mpp        ! NEMO MPI routines, ctl_stop in particular
   USE phycst         ! NEMO physical constants
   USE sbc_oce
   USE eosbn2         ! equation of state
   USE lib_fortran, ONLY : DDPDD

   USE icb_oce        ! define iceberg arrays
   USE icbutl         ! iceberg utility routines
   USE icbdia         ! iceberg budget routines

   IMPLICIT NONE
   PRIVATE

   PUBLIC   icb_thm ! routine called in icbstp.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE icb_thm( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_thm  ***
      !!
      !! ** Purpose :   compute the iceberg thermodynamics.
      !!
      !! ** Method  : - See Martin & Adcroft, Ocean Modelling 34, 2010
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! timestep number, just passed to icb_utl_print_berg
      !
      INTEGER  ::   ii, ij
      REAL(wp) ::   zM, zT, zW, zL, zSST, zVol, zLn, zWn, zTn, znVol, zIC, zDn
      REAL(wp) ::   zSSS, zfzpt
      REAL(wp) ::   zMv, zMe, zMb, zmelt, zdvo, zdva, zdM, zSs, zdMe, zdMb, zdMv
      REAL(wp) ::   zMnew, zMnew1, zMnew2, zheat_hcflux, zheat_latent, z1_12
      REAL(wp) ::   zMbits, znMbits, zdMbitsE, zdMbitsM, zLbits, zAbits, zMbb
      REAL(wp) ::   zxi, zyj, zff, z1_rday, z1_e1e2, zdt, z1_dt, z1_dt_e1e2
      TYPE(iceberg), POINTER ::   this, next
      TYPE(point)  , POINTER ::   pt
      !
      COMPLEX(wp), DIMENSION(jpi,jpj) :: cicb_melt, cicb_hflx
      !!----------------------------------------------------------------------
      !
      !! initialiaze cicb_melt and cicb_heat
      cicb_melt = CMPLX( 0.e0, 0.e0, wp ) 
      cicb_hflx = CMPLX( 0.e0, 0.e0, wp ) 
      !
      z1_rday = 1._wp / rday
      z1_12   = 1._wp / 12._wp
      zdt     = berg_dt
      z1_dt   = 1._wp / zdt
      !
      ! we're either going to ignore berg fresh water melt flux and associated heat
      ! or we pass it into the ocean, so at this point we set them both to zero,
      ! accumulate the contributions to them from each iceberg in the while loop following
      ! and then pass them (or not) to the ocean
      !
      berg_grid%floating_melt(:,:) = 0._wp
      ! calving_hflx re-used here as temporary workspace for the heat flux associated with melting
      berg_grid%calving_hflx(:,:)  = 0._wp
      !
      this => first_berg
      DO WHILE( ASSOCIATED(this) )
         !
         pt => this%current_point
         nknberg = this%number(1)
         CALL icb_utl_interp( pt%xi, pt%e1, pt%uo, pt%ui, pt%ua, pt%ssh_x,   &
            &                 pt%yj, pt%e2, pt%vo, pt%vi, pt%va, pt%ssh_y,   &
            &                 pt%sst, pt%cn, pt%hi, zff, pt%sss )
         !
         zSST = pt%sst
         zSSS = pt%sss
         CALL eos_fzp(zSSS,zfzpt)                       ! freezing point
         zIC  = MIN( 1._wp, pt%cn + rn_sicn_shift )     ! Shift sea-ice concentration       !!gm ???
         zM   = pt%mass
         zT   = pt%thickness                               ! total thickness
       ! D   = (rn_rho_bergs/pp_rho_seawater)*zT ! draught (keel depth)
       ! F   = zT - D ! freeboard
         zW   = pt%width
         zL   = pt%length
         zxi  = pt%xi                                      ! position in (i,j) referential
         zyj  = pt%yj
         ii  = INT( zxi + 0.5 )                            ! T-cell of the berg
         ii  = mi1( ii )
         ij  = INT( zyj + 0.5 )              
         ij  = mj1( ij )
         zVol = zT * zW * zL

         ! Environment
         zdvo = SQRT( (pt%uvel-pt%uo)**2 + (pt%vvel-pt%vo)**2 )
         zdva = SQRT( (pt%ua  -pt%uo)**2 + (pt%va  -pt%vo)**2 )
         zSs  = 1.5_wp * SQRT( zdva ) + 0.1_wp * zdva                ! Sea state      (eqn M.A9)

         ! Melt rates in m/s (i.e. division by rday)
         zMv = MAX( 7.62d-3*zSST+1.29d-3*(zSST**2)                    , 0._wp ) * z1_rday      ! Buoyant convection at sides (eqn M.A10)
         IF ( zSST > zfzpt ) THEN                                                              ! Calculate basal melting only if SST above freezing point  
            zMb = MAX( 0.58_wp*(zdvo**0.8_wp)*(zSST+4.0_wp)/(zL**0.2_wp) , 0._wp ) * z1_rday   ! Basal turbulent melting     (eqn M.A7 )
         ELSE
            zMb = 0._wp                                                                        ! No basal melting if SST below freezing point     
         ENDIF
         zMe = MAX( z1_12*(zSST+2.)*zSs*(1._wp+COS(rpi*(zIC**3)))     , 0._wp ) * z1_rday      ! Wave erosion                (eqn M.A8 )

         IF( ln_operator_splitting ) THEN      ! Operator split update of volume/mass
            zTn    = MAX( zT - zMb*zdt , 0._wp )         ! new total thickness (m)
            znVol  = zTn * zW * zL                       ! new volume (m^3)
            zMnew1 = ( znVol / zVol ) * zM               ! new mass (kg)
            zdMb   = zM - zMnew1                         ! mass lost to basal melting (>0) (kg)
            !
            zLn    = MAX( zL - zMv*zdt , 0._wp )         ! new length (m)
            zWn    = MAX( zW - zMv*zdt , 0._wp )         ! new width (m)
            znVol  = zTn * zWn * zLn                     ! new volume (m^3)
            zMnew2 = ( znVol / zVol ) * zM               ! new mass (kg)
            zdMv   = zMnew1 - zMnew2                     ! mass lost to buoyant convection (>0) (kg)
            !
            zLn    = MAX( zLn - zMe*zdt , 0._wp )        ! new length (m)
            zWn    = MAX( zWn - zMe*zdt , 0._wp )        ! new width (m)
            znVol  = zTn * zWn * zLn                     ! new volume (m^3)
            zMnew  = ( znVol / zVol ) * zM               ! new mass (kg)
            zdMe   = zMnew2 - zMnew                      ! mass lost to erosion (>0) (kg)
            zdM    = zM - zMnew                          ! mass lost to all erosion and melting (>0) (kg)
            !
         ELSE                                         ! Update dimensions of berg
            zLn = MAX( zL -(zMv+zMe)*zdt ,0._wp )        ! (m)
            zWn = MAX( zW -(zMv+zMe)*zdt ,0._wp )        ! (m)
            zTn = MAX( zT - zMb    *zdt ,0._wp )         ! (m)
            ! Update volume and mass of berg
            znVol = zTn*zWn*zLn                          ! (m^3)
            zMnew = (znVol/zVol)*zM                      ! (kg)
            zdM   = zM - zMnew                           ! (kg)
            zdMb = (zM/zVol) * (zW*   zL ) *zMb*zdt      ! approx. mass loss to basal melting (kg)
            zdMe = (zM/zVol) * (zT*(zW+zL)) *zMe*zdt     ! approx. mass lost to erosion (kg)
            zdMv = (zM/zVol) * (zT*(zW+zL)) *zMv*zdt     ! approx. mass loss to buoyant convection (kg)
         ENDIF

         IF( rn_bits_erosion_fraction > 0._wp ) THEN     ! Bergy bits
            !
            zMbits   = pt%mass_of_bits                                               ! mass of bergy bits (kg)
            zdMbitsE = rn_bits_erosion_fraction * zdMe                               ! change in mass of bits (kg)
            znMbits  = zMbits + zdMbitsE                                             ! add new bergy bits to mass (kg)
            zLbits   = MIN( zL, zW, zT, 40._wp )                                     ! assume bergy bits are smallest dimension or 40 meters
            zAbits   = ( zMbits / rn_rho_bergs ) / zLbits                            ! Effective bottom area (assuming T=Lbits)
            zMbb     = MAX( 0.58_wp*(zdvo**0.8_wp)*(zSST+2._wp) /   &
               &                              ( zLbits**0.2_wp ) , 0._wp ) * z1_rday ! Basal turbulent melting (for bits)
            zMbb     = rn_rho_bergs * zAbits * zMbb                                  ! in kg/s
            zdMbitsM = MIN( zMbb*zdt , znMbits )                                     ! bergy bits mass lost to melting (kg)
            znMbits  = znMbits-zdMbitsM                                              ! remove mass lost to bergy bits melt
            IF( zMnew == 0._wp ) THEN                                                ! if parent berg has completely melted then
               zdMbitsM = zdMbitsM + znMbits                                         ! instantly melt all the bergy bits
               znMbits  = 0._wp
            ENDIF
         ELSE                                                     ! No bergy bits
            zAbits   = 0._wp
            zdMbitsE = 0._wp
            zdMbitsM = 0._wp
            znMbits  = pt%mass_of_bits                             ! retain previous value incase non-zero
         ENDIF

         ! use tmask rather than tmask_i when dealing with icebergs
         IF( tmask(ii,ij,1) /= 0._wp ) THEN    ! Add melting to the grid and field diagnostics
            z1_e1e2    = r1_e1e2t(ii,ij) * this%mass_scaling
            z1_dt_e1e2 = z1_dt * z1_e1e2
            !
            ! iceberg melt
            !! the use of DDPDD function for the cumulative sum is needed for reproducibility
            zmelt    = ( zdM - ( zdMbitsE - zdMbitsM ) ) * z1_dt   ! kg/s
            CALL DDPDD( CMPLX( zmelt * z1_e1e2, 0.e0, wp ), cicb_melt(ii,ij) )
            !
            ! iceberg heat flux
            !! the use of DDPDD function for the cumulative sum is needed for reproducibility
            !! NB. The src_calving_hflx field is currently hardwired to zero in icb_stp, which means that the
            !!     heat density of the icebergs is zero and the heat content flux to the ocean from iceberg
            !!     melting is always zero. Leaving the term in the code until such a time as this is fixed. DS.
            zheat_hcflux = zmelt * pt%heat_density       ! heat content flux : kg/s x J/kg = J/s
            zheat_latent = - zmelt * rLfus               ! latent heat flux:  kg/s x J/kg = J/s
            CALL DDPDD( CMPLX( ( zheat_hcflux + zheat_latent ) * z1_e1e2, 0.e0, wp ), cicb_hflx(ii,ij) )
            !
            ! diagnostics
            CALL icb_dia_melt( ii, ij, zMnew, zheat_hcflux, zheat_latent, this%mass_scaling,       &
               &                       zdM, zdMbitsE, zdMbitsM, zdMb, zdMe,   &
               &                       zdMv, z1_dt_e1e2 )
         ELSE
            WRITE(numout,*) 'icb_thm: berg ',this%number(:),' appears to have grounded  at ',narea,ii,ij
            CALL icb_utl_print_berg( this, kt )
            WRITE(numout,*) 'msk=',tmask(ii,ij,1), e1e2t(ii,ij)
            CALL ctl_stop('icb_thm', 'berg appears to have grounded!')
         ENDIF

         ! Rolling
         zDn = ( rn_rho_bergs / pp_rho_seawater ) * zTn       ! draught (keel depth)
         IF( zDn > 0._wp .AND. MAX(zWn,zLn) < SQRT( 0.92*(zDn**2) + 58.32*zDn ) ) THEN
            zT  = zTn
            zTn = zWn
            zWn = zT
         ENDIF

         ! Store the new state of iceberg (with L>W)
         pt%mass         = zMnew
         pt%mass_of_bits = znMbits
         pt%thickness    = zTn
         pt%width        = MIN( zWn , zLn )
         pt%length       = MAX( zWn , zLn )

         next=>this%next

!!gm  add a test to avoid over melting ?

         IF( zMnew <= 0._wp ) THEN       ! Delete the berg if completely melted
            CALL icb_utl_delete( first_berg, this )
            !
         ELSE                            ! Diagnose mass distribution on grid
            z1_e1e2 = r1_e1e2t(ii,ij) * this%mass_scaling
            CALL icb_dia_size( ii, ij, zWn, zLn, zAbits,   &
               &               this%mass_scaling, zMnew, znMbits, z1_e1e2 )
         ENDIF
         !
         this=>next
         !
      END DO
      !
      berg_grid%floating_melt = REAL(cicb_melt,wp)    ! kg/m2/s
      berg_grid%calving_hflx  = REAL(cicb_hflx,wp)
      !
      ! now use melt and associated heat flux in ocean (or not)
      !
      IF(.NOT. ln_passive_mode ) THEN
         emp (:,:) = emp (:,:) - berg_grid%floating_melt(:,:)
         qns (:,:) = qns (:,:) + berg_grid%calving_hflx (:,:)  
      ENDIF
      !
   END SUBROUTINE icb_thm

   !!======================================================================
END MODULE icbthm
