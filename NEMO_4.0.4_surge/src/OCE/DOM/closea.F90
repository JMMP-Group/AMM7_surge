MODULE closea
   !!======================================================================
   !!                   ***  MODULE  closea  ***
   !!
   !! User define : specific treatments associated with closed seas
   !!======================================================================
   !! History :   8.2  !  2000-05  (O. Marti)  Original code
   !!   NEMO      1.0  !  2002-06  (E. Durand, G. Madec)  F90
   !!             3.0  !  2006-07  (G. Madec)  add clo_rnf, clo_ups, clo_bat
   !!             3.4  !  2014-12  (P.G. Fogli) sbc_clo bug fix & mpp reproducibility
   !!             4.0  !  2016-06  (G. Madec)  move to usrdef_closea, remove clo_ups
   !!             4.0  !  2017-12  (D. Storkey) new formulation based on masks read from file
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_clo    : read in masks which define closed seas and runoff areas
   !!   sbc_clo    : Special handling of freshwater fluxes over closed seas
   !!   clo_rnf    : set close sea outflows as river mouths (see sbcrnf)
   !!   clo_bat    : set to zero a field over closed sea (see domzgr)
   !!----------------------------------------------------------------------
   USE oce             ! dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! ocean surface boundary conditions
   USE iom             ! I/O routines
   !
   USE in_out_manager  ! I/O manager
   USE lib_fortran,    ONLY: glob_sum
   USE lbclnk          ! lateral boundary condition - MPP exchanges
   USE lib_mpp         ! MPP library
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC dom_clo      ! called by domain module
   PUBLIC sbc_clo      ! called by sbcmod module
   PUBLIC clo_rnf      ! called by sbcrnf module
   PUBLIC clo_bat      ! called in domzgr module

   LOGICAL, PUBLIC :: ln_closea  !:  T => keep closed seas (defined by closea_mask field) in the domain and apply
                                 !:       special treatment of freshwater fluxes.
                                 !:  F => suppress closed seas (defined by closea_mask field) from the bathymetry
                                 !:       at runtime.
                                 !:  If there is no closea_mask field in the domain_cfg file or we do not use
                                 !:  a domain_cfg file then this logical does nothing.
                                 !:
   LOGICAL, PUBLIC :: l_sbc_clo  !: T => Closed seas defined, apply special treatment of freshwater fluxes.
                                 !: F => No closed seas defined (closea_mask field not found).
   LOGICAL, PUBLIC :: l_clo_rnf  !: T => Some closed seas output freshwater (RNF or EMPMR) to specified runoff points.
   INTEGER, PUBLIC :: jncs       !: number of closed seas (inferred from closea_mask field)
   INTEGER, PUBLIC :: jncsr      !: number of closed seas rnf mappings (inferred from closea_mask_rnf field)
   INTEGER, PUBLIC :: jncse      !: number of closed seas empmr mappings (inferred from closea_mask_empmr field)
   
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::  closea_mask       !: mask of integers defining closed seas
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::  closea_mask_rnf   !: mask of integers defining closed seas rnf mappings
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::  closea_mask_empmr !: mask of integers defining closed seas empmr mappings
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:)  ::   surf         !: closed sea surface areas 
                                                                  !: (and residual global surface area) 
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:)  ::   surfr        !: closed sea target rnf surface areas 
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:)  ::   surfe        !: closed sea target empmr surface areas 

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_clo()
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dom_clo  ***
      !!        
      !! ** Purpose :   Closed sea domain initialization
      !!
      !! ** Method  :   if a closed sea is located only in a model grid point
      !!                just the thermodynamic processes are applied.
      !!
      !! ** Action  :   Read closea_mask* fields (if they exist) from domain_cfg file and infer
      !!                number of closed seas from closea_mask field.
      !!                closea_mask       : integer values defining closed seas (or groups of closed seas)
      !!                closea_mask_rnf   : integer values defining mappings from closed seas or groups of
      !!                                    closed seas to a runoff area for downwards flux only.
      !!                closea_mask_empmr : integer values defining mappings from closed seas or groups of
      !!                                    closed seas to a runoff area for net fluxes.
      !!
      !!                Python code to generate the closea_masks* fields from the old-style indices
      !!                definitions is available at TOOLS/DOMAINcfg/make_closea_masks.py
      !!----------------------------------------------------------------------
      INTEGER ::   inum    ! input file identifier
      INTEGER ::   ierr    ! error code
      INTEGER ::   id      ! netcdf variable ID

      REAL(wp), DIMENSION(jpi,jpj) :: zdata_in ! temporary real array for input
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*)'dom_clo : read in masks to define closed seas '
      IF(lwp) WRITE(numout,*)'~~~~~~~'
      !
      ! read the closed seas masks (if they exist) from domain_cfg file (if it exists)
      ! ------------------------------------------------------------------------------
      !
      IF( ln_read_cfg) THEN
         !
         CALL iom_open( cn_domcfg, inum )
         !
         id = iom_varid(inum, 'closea_mask', ldstop = .false.)
         IF( id > 0 ) THEN 
            l_sbc_clo = .true.
            ALLOCATE( closea_mask(jpi,jpj) , STAT=ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'dom_clo: failed to allocate closea_mask array')
            zdata_in(:,:) = 0.0
            CALL iom_get ( inum, jpdom_data, 'closea_mask', zdata_in )
            closea_mask(:,:) = NINT(zdata_in(:,:)) * tmask(:,:,1)
            ! number of closed seas = global maximum value in closea_mask field
            jncs = maxval(closea_mask(:,:))
            CALL mpp_max('closea', jncs)
            IF( jncs > 0 ) THEN
               IF( lwp ) WRITE(numout,*) 'Number of closed seas : ',jncs
            ELSE
               CALL ctl_stop( 'Problem with closea_mask field in domain_cfg file. Has no values > 0 so no closed seas defined.')
            ENDIF
         ELSE 
            IF( lwp ) WRITE(numout,*)
            IF( lwp ) WRITE(numout,*) '   ==>>>   closea_mask field not found in domain_cfg file.'
            IF( lwp ) WRITE(numout,*) '           No closed seas defined.'
            IF( lwp ) WRITE(numout,*)
            l_sbc_clo = .false.
            jncs = 0 
         ENDIF

         l_clo_rnf = .false.

         IF( l_sbc_clo ) THEN ! No point reading in closea_mask_rnf or closea_mask_empmr fields if no closed seas defined.

            id = iom_varid(inum, 'closea_mask_rnf', ldstop = .false.)
            IF( id > 0 ) THEN 
               l_clo_rnf = .true.            
               ALLOCATE( closea_mask_rnf(jpi,jpj) , STAT=ierr )
               IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'dom_clo: failed to allocate closea_mask_rnf array')
               CALL iom_get ( inum, jpdom_data, 'closea_mask_rnf', zdata_in )
               closea_mask_rnf(:,:) = NINT(zdata_in(:,:)) * tmask(:,:,1)
               ! number of closed seas rnf mappings = global maximum in closea_mask_rnf field
               jncsr = maxval(closea_mask_rnf(:,:))
               CALL mpp_max('closea', jncsr)
               IF( jncsr > 0 ) THEN
                  IF( lwp ) WRITE(numout,*) 'Number of closed seas rnf mappings : ',jncsr
               ELSE
                  CALL ctl_stop( 'Problem with closea_mask_rnf field in domain_cfg file. Has no values > 0 so no closed seas rnf mappings defined.')
               ENDIF
            ELSE 
               IF( lwp ) WRITE(numout,*) 'closea_mask_rnf field not found in domain_cfg file. No closed seas rnf mappings defined.'
               jncsr = 0
            ENDIF
 
            id = iom_varid(inum, 'closea_mask_empmr', ldstop = .false.)
            IF( id > 0 ) THEN 
               l_clo_rnf = .true.            
               ALLOCATE( closea_mask_empmr(jpi,jpj) , STAT=ierr )
               IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'dom_clo: failed to allocate closea_mask_empmr array')
               CALL iom_get ( inum, jpdom_data, 'closea_mask_empmr', zdata_in )
               closea_mask_empmr(:,:) = NINT(zdata_in(:,:)) * tmask(:,:,1)
               ! number of closed seas empmr mappings = global maximum value in closea_mask_empmr field
               jncse = maxval(closea_mask_empmr(:,:))
               CALL mpp_max('closea', jncse)
               IF( jncse > 0 ) THEN 
                  IF( lwp ) WRITE(numout,*) 'Number of closed seas empmr mappings : ',jncse
               ELSE
                  CALL ctl_stop( 'Problem with closea_mask_empmr field in domain_cfg file. Has no values > 0 so no closed seas empmr mappings defined.')
               ENDIF
            ELSE 
               IF( lwp ) WRITE(numout,*) 'closea_mask_empmr field not found in domain_cfg file. No closed seas empmr mappings defined.'
               jncse = 0
            ENDIF

         ENDIF ! l_sbc_clo
         !
         CALL iom_close( inum )
         !
      ELSE ! ln_read_cfg = .false. so no domain_cfg file
         IF( lwp ) WRITE(numout,*) 'No domain_cfg file so no closed seas defined.'
         l_sbc_clo = .false.
         l_clo_rnf = .false.
      ENDIF
      !
   END SUBROUTINE dom_clo


   SUBROUTINE sbc_clo( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_clo  ***
      !!                    
      !! ** Purpose :   Special handling of closed seas
      !!
      !! ** Method  :   Water flux is forced to zero over closed sea
      !!      Excess is shared between remaining ocean, or
      !!      put as run-off in open ocean.
      !!
      !! ** Action  :   emp updated surface freshwater fluxes and associated heat content at kt
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in   ) ::   kt       ! ocean model time step
      !
      INTEGER             ::   ierr
      INTEGER             ::   jc, jcr, jce   ! dummy loop indices
      REAL(wp), PARAMETER ::   rsmall = 1.e-20_wp    ! Closed sea correction epsilon
      REAL(wp)            ::   zfwf_total, zcoef, zcoef1         ! 
      REAL(wp), DIMENSION(jncs)    ::   zfwf      !:
      REAL(wp), DIMENSION(jncsr+1) ::   zfwfr     !: freshwater fluxes over closed seas
      REAL(wp), DIMENSION(jncse+1) ::   zfwfe     !: 
      REAL(wp), DIMENSION(jpi,jpj) ::   ztmp2d   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('sbc_clo')
      !
      !                                                   !------------------! 
      IF( kt == nit000 ) THEN                             !  Initialisation  !
         !                                                !------------------!
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'sbc_clo : closed seas '
         IF(lwp) WRITE(numout,*)'~~~~~~~'

         ALLOCATE( surf(jncs+1) , STAT=ierr )
         IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'sbc_clo: failed to allocate surf array')
         surf(:) = 0.e0_wp
         !
         ! jncsr can be zero so add 1 to avoid allocating zero-length array
         ALLOCATE( surfr(jncsr+1) , STAT=ierr )
         IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'sbc_clo: failed to allocate surfr array')
         surfr(:) = 0.e0_wp
         !
         ! jncse can be zero so add 1 to avoid allocating zero-length array
         ALLOCATE( surfe(jncse+1) , STAT=ierr )
         IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'sbc_clo: failed to allocate surfe array')
         surfe(:) = 0.e0_wp
         !
         surf(jncs+1) = glob_sum( 'closea', e1e2t(:,:) )   ! surface of the global ocean
         !
         !                                        ! surface areas of closed seas 
         DO jc = 1, jncs
            ztmp2d(:,:) = 0.e0_wp
            WHERE( closea_mask(:,:) == jc ) ztmp2d(:,:) = e1e2t(:,:) * tmask_i(:,:)
            surf(jc) = glob_sum( 'closea', ztmp2d(:,:) )
         END DO
         !
         ! jncs+1 : surface area of global ocean, closed seas excluded
         surf(jncs+1) = surf(jncs+1) - SUM(surf(1:jncs))
         !
         !                                        ! surface areas of rnf target areas
         IF( jncsr > 0 ) THEN
            DO jcr = 1, jncsr
               ztmp2d(:,:) = 0.e0_wp
               WHERE( closea_mask_rnf(:,:) == jcr .and. closea_mask(:,:) == 0 ) ztmp2d(:,:) = e1e2t(:,:) * tmask_i(:,:)
               surfr(jcr) = glob_sum( 'closea', ztmp2d(:,:) )
            END DO
         ENDIF
         !
         !                                        ! surface areas of empmr target areas
         IF( jncse > 0 ) THEN
            DO jce = 1, jncse
               ztmp2d(:,:) = 0.e0_wp
               WHERE( closea_mask_empmr(:,:) == jce .and. closea_mask(:,:) == 0 ) ztmp2d(:,:) = e1e2t(:,:) * tmask_i(:,:)
               surfe(jce) = glob_sum( 'closea', ztmp2d(:,:) )
            END DO
         ENDIF
         !
         IF(lwp) WRITE(numout,*)'     Closed sea surface areas (km2)'
         DO jc = 1, jncs
            IF(lwp) WRITE(numout,FMT='(1I3,5X,ES12.2)') jc, surf(jc) * 1.0e-6
         END DO
         IF(lwp) WRITE(numout,FMT='(A,ES12.2)') 'Global surface area excluding closed seas (km2): ', surf(jncs+1) * 1.0e-6
         !
         IF(jncsr > 0) THEN
            IF(lwp) WRITE(numout,*)'     Closed sea target rnf surface areas (km2)'
            DO jcr = 1, jncsr
               IF(lwp) WRITE(numout,FMT='(1I3,5X,ES12.2)') jcr, surfr(jcr) * 1.0e-6
            END DO
         ENDIF
         !
         IF(jncse > 0) THEN
            IF(lwp) WRITE(numout,*)'     Closed sea target empmr surface areas (km2)'
            DO jce = 1, jncse
               IF(lwp) WRITE(numout,FMT='(1I3,5X,ES12.2)') jce, surfe(jce) * 1.0e-6
            END DO
         ENDIF
      ENDIF
      !
      !                                                      !--------------------!
      !                                                      !  update emp        !
      !                                                      !--------------------!

      zfwf_total = 0._wp

      !
      ! 1. Work out total freshwater fluxes over closed seas from EMP - RNF.
      !
      zfwf(:) = 0.e0_wp           
      DO jc = 1, jncs
         ztmp2d(:,:) = 0.e0_wp
         WHERE( closea_mask(:,:) == jc ) ztmp2d(:,:) = e1e2t(:,:) * ( emp(:,:)-rnf(:,:) ) * tmask_i(:,:)
         zfwf(jc) = glob_sum( 'closea', ztmp2d(:,:) )
      END DO
      zfwf_total = SUM(zfwf)

      zfwfr(:) = 0.e0_wp           
      IF( jncsr > 0 ) THEN
         !
         ! 2. Work out total FW fluxes over rnf source areas and add to rnf target areas. 
         !    Where zfwf is negative add flux at specified runoff points and subtract from fluxes for global redistribution.
         !    Where positive leave in global redistribution total.
         !
         DO jcr = 1, jncsr
            !
            ztmp2d(:,:) = 0.e0_wp
            WHERE( closea_mask_rnf(:,:) == jcr .and. closea_mask(:,:) > 0 ) ztmp2d(:,:) = e1e2t(:,:) * ( emp(:,:)-rnf(:,:) ) * tmask_i(:,:)
            zfwfr(jcr) = glob_sum( 'closea', ztmp2d(:,:) )
            !
            ! The following if avoids the redistribution of the round off
            IF ( ABS(zfwfr(jcr) / surf(jncs+1) ) > rsmall) THEN
               !
               ! Add residuals to target runoff points if negative and subtract from total to be added globally
               IF( zfwfr(jcr) < 0.0 ) THEN 
                  zfwf_total = zfwf_total - zfwfr(jcr)
                  zcoef    = zfwfr(jcr) / surfr(jcr)
                  zcoef1   = rcp * zcoef
                  WHERE( closea_mask_rnf(:,:) == jcr .and. closea_mask(:,:) == 0.0)
                     emp(:,:) = emp(:,:) + zcoef
                     qns(:,:) = qns(:,:) - zcoef1 * sst_m(:,:)
                  ENDWHERE
               ENDIF
               !
            ENDIF
         END DO
      ENDIF  ! jncsr > 0    
      !
      zfwfe(:) = 0.e0_wp           
      IF( jncse > 0 ) THEN
         !
         ! 3. Work out total fluxes over empmr source areas and add to empmr target areas. 
         !
         DO jce = 1, jncse
            !
            ztmp2d(:,:) = 0.e0_wp
            WHERE( closea_mask_empmr(:,:) == jce .and. closea_mask(:,:) > 0 ) ztmp2d(:,:) = e1e2t(:,:) * ( emp(:,:)-rnf(:,:) ) * tmask_i(:,:)
            zfwfe(jce) = glob_sum( 'closea', ztmp2d(:,:) )
            !
            ! The following if avoids the redistribution of the round off
            IF ( ABS( zfwfe(jce) / surf(jncs+1) ) > rsmall ) THEN
               !
               ! Add residuals to runoff points and subtract from total to be added globally
               zfwf_total = zfwf_total - zfwfe(jce)
               zcoef    = zfwfe(jce) / surfe(jce)
               zcoef1   = rcp * zcoef
               WHERE( closea_mask_empmr(:,:) == jce .and. closea_mask(:,:) == 0.0)
                  emp(:,:) = emp(:,:) + zcoef
                  qns(:,:) = qns(:,:) - zcoef1 * sst_m(:,:)
               ENDWHERE
               !
            ENDIF
         END DO
      ENDIF ! jncse > 0    

      !
      ! 4. Spread residual flux over global ocean. 
      !
      ! The following if avoids the redistribution of the round off
      IF ( ABS(zfwf_total / surf(jncs+1) ) > rsmall) THEN
         zcoef    = zfwf_total / surf(jncs+1)
         zcoef1   = rcp * zcoef
         WHERE( closea_mask(:,:) == 0 )
            emp(:,:) = emp(:,:) + zcoef
            qns(:,:) = qns(:,:) - zcoef1 * sst_m(:,:)
         ENDWHERE
      ENDIF

      !
      ! 5. Subtract area means from emp (and qns) over closed seas to give zero mean FW flux over each sea.
      !
      DO jc = 1, jncs
         ! The following if avoids the redistribution of the round off
         IF ( ABS(zfwf(jc) / surf(jncs+1) ) > rsmall) THEN
            !
            ! Subtract residuals from fluxes over closed sea
            zcoef    = zfwf(jc) / surf(jc)
            zcoef1   = rcp * zcoef
            WHERE( closea_mask(:,:) == jc )
               emp(:,:) = emp(:,:) - zcoef
               qns(:,:) = qns(:,:) + zcoef1 * sst_m(:,:)
            ENDWHERE
            !
         ENDIF
      END DO
      !
      emp (:,:) = emp (:,:) * tmask(:,:,1)
      !
      CALL lbc_lnk( 'closea', emp , 'T', 1._wp )
      !
   END SUBROUTINE sbc_clo

   SUBROUTINE clo_rnf( p_rnfmsk )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf  ***
      !!                    
      !! ** Purpose :   allow the treatment of closed sea outflow grid-points
      !!                to be the same as river mouth grid-points
      !!
      !! ** Method  :   set to 1 the runoff mask (mskrnf, see sbcrnf module)
      !!                at the closed sea outflow grid-point.
      !!
      !! ** Action  :   update (p_)mskrnf (set 1 at closed sea outflow)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   p_rnfmsk   ! river runoff mask (rnfmsk array)
      !!----------------------------------------------------------------------
      !
      IF( jncsr > 0 ) THEN
         WHERE( closea_mask_rnf(:,:) > 0 .and. closea_mask(:,:) == 0 )
            p_rnfmsk(:,:) = MAX( p_rnfmsk(:,:), 1.0_wp )
         ENDWHERE
      ENDIF
      !
      IF( jncse > 0 ) THEN
         WHERE( closea_mask_empmr(:,:) > 0 .and. closea_mask(:,:) == 0 )
            p_rnfmsk(:,:) = MAX( p_rnfmsk(:,:), 1.0_wp )
         ENDWHERE
      ENDIF
      !
   END SUBROUTINE clo_rnf
   
      
   SUBROUTINE clo_bat( k_top, k_bot )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE clo_bat  ***
      !!                    
      !! ** Purpose :   Suppress closed sea from the domain
      !!
      !! ** Method  :   Read in closea_mask field (if it exists) from domain_cfg file.
      !!                Where closea_mask > 0 set first and last ocean level to 0
      !!                (As currently coded you can't define a closea_mask field in 
      !!                usr_def_zgr).
      !!
      !! ** Action  :   set k_top=0 and k_bot=0 over closed seas
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(:,:), INTENT(inout) ::   k_top, k_bot   ! ocean first and last level indices
      INTEGER                           :: inum, id
      INTEGER,  DIMENSION(jpi,jpj) :: closea_mask ! closea_mask field
      REAL(wp), DIMENSION(jpi,jpj) :: zdata_in ! temporary real array for input
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                     ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'clo_bat : suppression of closed seas'
         WRITE(numout,*) '~~~~~~~'
      ENDIF
      !
      IF( ln_read_cfg ) THEN
         !
         CALL iom_open( cn_domcfg, inum )
         !
         id = iom_varid(inum, 'closea_mask', ldstop = .false.)      
         IF( id > 0 ) THEN
            IF( lwp ) WRITE(numout,*) 'Suppressing closed seas in bathymetry based on closea_mask field,'
            CALL iom_get ( inum, jpdom_data, 'closea_mask', zdata_in )
            closea_mask(:,:) = NINT(zdata_in(:,:))
            WHERE( closea_mask(:,:) > 0 )
               k_top(:,:) = 0   
               k_bot(:,:) = 0   
            ENDWHERE
         ELSE
            IF( lwp ) WRITE(numout,*) 'No closea_mask field found in domain_cfg file. No suppression of closed seas.'
         ENDIF
         !
         CALL iom_close(inum)
         !
      ELSE
         IF( lwp ) WRITE(numout,*) 'No domain_cfg file => no suppression of closed seas.'
      ENDIF
      !
      ! Initialise l_sbc_clo and l_clo_rnf for this case (ln_closea=.false.)
      l_sbc_clo = .false.
      l_clo_rnf = .false.
      !
   END SUBROUTINE clo_bat

   !!======================================================================
END MODULE closea

