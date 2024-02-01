MODULE updtide
   !!======================================================================
   !!                       ***  MODULE  updtide  ***
   !! Initialization of tidal forcing
   !!======================================================================
   !! History :  9.0  !  07  (O. Le Galloudec)  Original code
   !!----------------------------------------------------------------------
   !!   upd_tide       : update tidal potential
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O units
   USE phycst          ! physical constant
   USE sbctide         ! tide potential variable
   USE tideini, ONLY: ln_tide_ramp, rdttideramp

   IMPLICIT NONE
   PUBLIC

   PUBLIC   upd_tide   ! called in dynspg_... modules
  
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE upd_tide( kt, kit, kt_offset )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE upd_tide  ***
      !!
      !! ** Purpose :   provide at each time step the astronomical potential
      !!
      !! ** Method  :   computed from pulsation and amplitude of all tide components
      !!
      !! ** Action  :   pot_astro   actronomical potential
      !!----------------------------------------------------------------------      
      INTEGER, INTENT(in)           ::   kt      ! ocean time-step index
      INTEGER, INTENT(in), OPTIONAL ::   kit     ! external mode sub-time-step index (lk_dynspg_ts=T)
      INTEGER, INTENT(in), OPTIONAL ::   kt_offset ! time offset in number 
                                                     ! of internal steps             (lk_dynspg_ts=F)
                                                     ! of external steps             (lk_dynspg_ts=T)
      !
      INTEGER  ::   ioffset      ! local integer
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zt, zramp    ! local scalar
      REAL(wp), DIMENSION(nb_harmo) ::   zwt 
      !!----------------------------------------------------------------------      
      !
      !                               ! tide pulsation at model time step (or sub-time-step)
      zt = ( kt - kt_tide ) * rdt
      !
      ioffset = 0
      IF( PRESENT( kt_offset ) )   ioffset = kt_offset
      !
      IF( PRESENT( kit ) )   THEN
         zt = zt + ( kit +  ioffset - 1 ) * rdt / REAL( nn_baro, wp )
      ELSE
         zt = zt + ioffset * rdt
      ENDIF
      !
      zwt(:) = omega_tide(:) * zt

      pot_astro(:,:) = 0._wp          ! update tidal potential (sum of all harmonics)
      DO jk = 1, nb_harmo   
         pot_astro(:,:) = pot_astro(:,:) + amp_pot(:,:,jk) * COS( zwt(jk) + phi_pot(:,:,jk) )      
      END DO
      !
      IF( ln_tide_ramp ) THEN         ! linear increase if asked
         zt = ( kt - nit000 ) * rdt
         IF( PRESENT( kit ) )   zt = zt + ( kit + ioffset -1) * rdt / REAL( nn_baro, wp )
         zramp = MIN(  MAX( zt / (rdttideramp*rday) , 0._wp ) , 1._wp  )
         pot_astro(:,:) = zramp * pot_astro(:,:)
      ENDIF
      !
   END SUBROUTINE upd_tide

  !!======================================================================

END MODULE updtide
