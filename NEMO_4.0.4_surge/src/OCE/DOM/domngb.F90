MODULE domngb
   !!======================================================================
   !!                       ***  MODULE  domngb  ***
   !! Grid search:  find the closest grid point from a given on/lat position
   !!======================================================================
   !! History : 3.2  !  2009-11  (S. Masson)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_ngb       : find the closest grid point from a given lon/lat position
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! for mppsum

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_ngb   ! routine called in iom.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_ngb( plon, plat, kii, kjj, cdgrid, kkk )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dom_ngb  ***
      !!
      !! ** Purpose :   find the closest grid point from a given lon/lat position
      !!
      !! ** Method  :   look for minimum distance in cylindrical projection 
      !!                -> not good if located at too high latitude...
      !!----------------------------------------------------------------------
      REAL(wp)        , INTENT(in   ) ::   plon, plat   ! longitude,latitude of the point
      INTEGER         , INTENT(  out) ::   kii, kjj     ! i-,j-index of the closes grid point
      INTEGER         , INTENT(in   ), OPTIONAL :: kkk  ! k-index of the mask level used
      CHARACTER(len=1), INTENT(in   ) ::   cdgrid       ! grid name 'T', 'U', 'V', 'W'
      !
      INTEGER :: ik         ! working level
      INTEGER , DIMENSION(2) ::   iloc
      REAL(wp)               ::   zlon, zmini
      REAL(wp), DIMENSION(jpi,jpj) ::   zglam, zgphi, zmask, zdist
      !!--------------------------------------------------------------------
      !
      zmask(:,:) = 0._wp
      ik = 1
      IF ( PRESENT(kkk) ) ik=kkk
      SELECT CASE( cdgrid )
      CASE( 'U' )  ; zglam(:,:) = glamu(:,:) ; zgphi(:,:) = gphiu(:,:) ; zmask(nldi:nlei,nldj:nlej) = umask(nldi:nlei,nldj:nlej,ik)
      CASE( 'V' )  ; zglam(:,:) = glamv(:,:) ; zgphi(:,:) = gphiv(:,:) ; zmask(nldi:nlei,nldj:nlej) = vmask(nldi:nlei,nldj:nlej,ik)
      CASE( 'F' )  ; zglam(:,:) = glamf(:,:) ; zgphi(:,:) = gphif(:,:) ; zmask(nldi:nlei,nldj:nlej) = fmask(nldi:nlei,nldj:nlej,ik)
      CASE DEFAULT ; zglam(:,:) = glamt(:,:) ; zgphi(:,:) = gphit(:,:) ; zmask(nldi:nlei,nldj:nlej) = tmask(nldi:nlei,nldj:nlej,ik)
      END SELECT

      zlon       = MOD( plon       + 720., 360. )                                     ! plon between    0 and 360
      zglam(:,:) = MOD( zglam(:,:) + 720., 360. )                                     ! glam between    0 and 360
      IF( zlon > 270. )   zlon = zlon - 360.                                          ! zlon between  -90 and 270
      IF( zlon <  90. )   WHERE( zglam(:,:) > 180. ) zglam(:,:) = zglam(:,:) - 360.   ! glam between -180 and 180
      zglam(:,:) = zglam(:,:) - zlon

      zgphi(:,:) = zgphi(:,:) - plat
      zdist(:,:) = zglam(:,:) * zglam(:,:) + zgphi(:,:) * zgphi(:,:)
      
      IF( lk_mpp ) THEN  
         CALL mpp_minloc( 'domngb', zdist(:,:), zmask, zmini, iloc)
         kii = iloc(1) ; kjj = iloc(2)
      ELSE
         zmini = MAXVAL( zdist(:,:)*zmask(:,:) )
         iloc(:) = MAXLOC( (zmini-zdist(:,:))*zmask(:,:) )
         kii = iloc(1) + nimpp - 1
         kjj = iloc(2) + njmpp - 1
      ENDIF
      !
   END SUBROUTINE dom_ngb

   !!======================================================================
END MODULE domngb
