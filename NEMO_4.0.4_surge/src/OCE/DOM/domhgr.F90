MODULE domhgr
   !!==============================================================================
   !!                       ***  MODULE domhgr   ***
   !! Ocean initialization : domain initialization
   !!==============================================================================
   !! History :  OPA  ! 1988-03  (G. Madec) Original code
   !!            7.0  ! 1996-01  (G. Madec)  terrain following coordinates
   !!            8.0  ! 1997-02  (G. Madec)  print mesh informations
   !!            8.1  ! 1999-11  (M. Imbard) NetCDF format with IO-IPSL
   !!            8.2  ! 2000-08  (D. Ludicone) Reduced section at Bab el Mandeb
   !!             -   ! 2001-09  (M. Levy)  eel config: grid in km, beta-plane
   !!  NEMO      1.0  ! 2002-08  (G. Madec)  F90: Free form and module, namelist
   !!             -   ! 2004-01  (A.M. Treguier, J.M. Molines) Case 4 (Mercator mesh)
   !!                            use of parameters in par_CONFIG-Rxx.h90, not in namelist
   !!             -   ! 2004-05  (A. Koch-Larrouy) Add Gyre configuration 
   !!            3.7  ! 2015-09  (G. Madec, S. Flavoni) add cell surface and their inverse
   !!                                       add optional read of e1e2u & e1e2v
   !!             -   ! 2016-04  (S. Flavoni, G. Madec) new configuration interface: read or usrdef.F90
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_hgr       : initialize the horizontal mesh 
   !!   hgr_read      : read horizontal information in the domain configuration file 
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE usrdef_hgr     ! User defined routine
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_hgr   ! called by domain.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_hgr
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_hgr  ***
      !!
      !! ** Purpose :   Read or compute the geographical position (in degrees)  
      !!      of the model grid-points, the horizontal scale factors (in meters), 
      !!      the associated horizontal metrics, and the Coriolis factor (in s-1).
      !!
      !! ** Method  :   Controlled by ln_read_cfg logical
      !!              =T : all needed arrays are read in mesh_mask.nc file 
      !!              =F : user-defined configuration, all needed arrays 
      !!                   are computed in usr-def_hgr subroutine 
      !!
      !!                If Coriolis factor is neither read nor computed (iff=0)
      !!              it is computed from gphit assuming that the mesh is
      !!              defined on the sphere :
      !!                   ff = 2.*omega*sin(gphif)      (in s-1)
      !!    
      !!                If u- & v-surfaces are neither read nor computed (ie1e2u_v=0)
      !!              (i.e. no use of reduced scale factors in some straits)
      !!              they are computed from e1u, e2u, e1v and e2v as:
      !!                   e1e2u = e1u*e2u   and   e1e2v = e1v*e2v  
      !!    
      !! ** Action  : - define longitude & latitude of t-, u-, v- and f-points (in degrees)
      !!              - define Coriolis parameter at f-point                   (in 1/s)
      !!              - define i- & j-scale factors at t-, u-, v- and f-points (in meters)
      !!              - define associated horizontal metrics at t-, u-, v- and f-points
      !!                (inverse of scale factors 1/e1 & 1/e2, surface e1*e2, ratios e1/e2 & e2/e1)
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj     ! dummy loop indices
      INTEGER ::   ie1e2u_v   ! flag for u- & v-surfaces 
      INTEGER ::   iff        ! flag for Coriolis parameter
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dom_hgr')
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_hgr : define the horizontal mesh from ithe following par_oce parameters '
         WRITE(numout,*) '~~~~~~~   '
         WRITE(numout,*) '   namcfg : read (=T) or user defined (=F) configuration    ln_read_cfg  = ', ln_read_cfg
      ENDIF
      !
      !
      IF( ln_read_cfg ) THEN        !==  read in mesh_mask.nc file  ==!
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ==>>>   read horizontal mesh in ', TRIM( cn_domcfg ), ' file'
         !
         CALL hgr_read   ( glamt , glamu , glamv , glamf ,   &    ! geographic position (required)
            &              gphit , gphiu , gphiv , gphif ,   &    !     -        -
            &              iff   , ff_f  , ff_t  ,           &    ! Coriolis parameter (if not on the sphere)
            &              e1t   , e1u   , e1v   , e1f   ,   &    ! scale factors (required)
            &              e2t   , e2u   , e2v   , e2f   ,   &    !    -     -        -
            &              ie1e2u_v      , e1e2u , e1e2v     )    ! u- & v-surfaces (if gridsize reduction in some straits)
         !
      ELSE                          !==  User defined configuration  ==! 
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          User defined horizontal mesh (usr_def_hgr)'
         !
         CALL usr_def_hgr( glamt , glamu , glamv , glamf ,   &    ! geographic position (required)
            &              gphit , gphiu , gphiv , gphif ,   &    !
            &              iff   , ff_f  , ff_t  ,           &    ! Coriolis parameter  (if domain not on the sphere)
            &              e1t   , e1u   , e1v   , e1f   ,   &    ! scale factors       (required)
            &              e2t   , e2u   , e2v   , e2f   ,   &    !
            &              ie1e2u_v      , e1e2u , e1e2v     )    ! u- & v-surfaces (if gridsize reduction is used in strait(s))
         !
      ENDIF
      !
      !                             !==  Coriolis parameter  ==!   (if necessary)
      !
      IF( iff == 0 ) THEN                 ! Coriolis parameter has not been defined 
         IF(lwp) WRITE(numout,*) '          Coriolis parameter calculated on the sphere from gphif & gphit'
         ff_f(:,:) = 2. * omega * SIN( rad * gphif(:,:) )     ! compute it on the sphere at f-point
         ff_t(:,:) = 2. * omega * SIN( rad * gphit(:,:) )     !    -        -       -    at t-point
      ELSE
         IF( ln_read_cfg ) THEN
            IF(lwp) WRITE(numout,*) '          Coriolis parameter have been read in ', TRIM( cn_domcfg ), ' file'
         ELSE
            IF(lwp) WRITE(numout,*) '          Coriolis parameter have been set in usr_def_hgr routine'
         ENDIF
      ENDIF

      !
      !                             !==  associated horizontal metrics  ==!
      !
      r1_e1t(:,:) = 1._wp / e1t(:,:)   ;   r1_e2t (:,:) = 1._wp / e2t(:,:)
      r1_e1u(:,:) = 1._wp / e1u(:,:)   ;   r1_e2u (:,:) = 1._wp / e2u(:,:)
      r1_e1v(:,:) = 1._wp / e1v(:,:)   ;   r1_e2v (:,:) = 1._wp / e2v(:,:)
      r1_e1f(:,:) = 1._wp / e1f(:,:)   ;   r1_e2f (:,:) = 1._wp / e2f(:,:)
      !
      e1e2t (:,:) = e1t(:,:) * e2t(:,:)   ;   r1_e1e2t(:,:) = 1._wp / e1e2t(:,:)
      e1e2f (:,:) = e1f(:,:) * e2f(:,:)   ;   r1_e1e2f(:,:) = 1._wp / e1e2f(:,:)
      IF( ie1e2u_v == 0 ) THEN               ! u- & v-surfaces have not been defined
         IF(lwp) WRITE(numout,*) '          u- & v-surfaces calculated as e1 e2 product'
         e1e2u (:,:) = e1u(:,:) * e2u(:,:)         ! compute them
         e1e2v (:,:) = e1v(:,:) * e2v(:,:) 
      ELSE
         IF(lwp) WRITE(numout,*) '          u- & v-surfaces have been read in "mesh_mask" file:'
         IF(lwp) WRITE(numout,*) '                     grid size reduction in strait(s) is used'
      ENDIF
      r1_e1e2u(:,:) = 1._wp / e1e2u(:,:)     ! compute their invert in any cases
      r1_e1e2v(:,:) = 1._wp / e1e2v(:,:)
      !   
      e2_e1u(:,:) = e2u(:,:) / e1u(:,:)
      e1_e2v(:,:) = e1v(:,:) / e2v(:,:)
      !
      !
      IF( ln_timing )   CALL timing_stop('dom_hgr')
      !
   END SUBROUTINE dom_hgr


   SUBROUTINE hgr_read( plamt , plamu , plamv  , plamf  ,   &    ! gridpoints position (required)
      &                 pphit , pphiu , pphiv  , pphif  ,   &     
      &                 kff   , pff_f , pff_t  ,            &    ! Coriolis parameter  (if not on the sphere)
      &                 pe1t  , pe1u  , pe1v   , pe1f   ,   &    ! scale factors       (required)
      &                 pe2t  , pe2u  , pe2v   , pe2f   ,   &
      &                 ke1e2u_v      , pe1e2u , pe1e2v     )    ! u- & v-surfaces (if gridsize reduction in some straits)
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE hgr_read  ***
      !!
      !! ** Purpose :   Read a mesh_mask file in NetCDF format using IOM
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   plamt, plamu, plamv, plamf   ! longitude outputs 
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pphit, pphiu, pphiv, pphif   ! latitude outputs
      INTEGER                 , INTENT(out) ::   kff                          ! =1 Coriolis parameter read here, =0 otherwise
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pff_f, pff_t                 ! Coriolis factor at f-point (if found in file)
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1t, pe1u, pe1v, pe1f       ! i-scale factors 
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe2t, pe2u, pe2v, pe2f       ! j-scale factors
      INTEGER                 , INTENT(out) ::   ke1e2u_v                     ! =1 u- & v-surfaces read here, =0 otherwise 
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1e2u, pe1e2v              ! u- & v-surfaces (if found in file)
      !
      INTEGER  ::   inum                  ! logical unit
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   hgr_read : read the horizontal coordinates in mesh_mask'
         WRITE(numout,*) '   ~~~~~~~~      jpiglo = ', jpiglo, ' jpjglo = ', jpjglo, ' jpk = ', jpk
      ENDIF
      !
      CALL iom_open( cn_domcfg, inum )
      !
      CALL iom_get( inum, jpdom_data, 'glamt', plamt, lrowattr=ln_use_jattr )
      CALL iom_get( inum, jpdom_data, 'glamu', plamu, lrowattr=ln_use_jattr )
      CALL iom_get( inum, jpdom_data, 'glamv', plamv, lrowattr=ln_use_jattr )
      CALL iom_get( inum, jpdom_data, 'glamf', plamf, lrowattr=ln_use_jattr )
      !
      CALL iom_get( inum, jpdom_data, 'gphit', pphit, lrowattr=ln_use_jattr )
      CALL iom_get( inum, jpdom_data, 'gphiu', pphiu, lrowattr=ln_use_jattr )
      CALL iom_get( inum, jpdom_data, 'gphiv', pphiv, lrowattr=ln_use_jattr )
      CALL iom_get( inum, jpdom_data, 'gphif', pphif, lrowattr=ln_use_jattr )
      !
      CALL iom_get( inum, jpdom_data, 'e1t'  , pe1t  , lrowattr=ln_use_jattr )
      CALL iom_get( inum, jpdom_data, 'e1u'  , pe1u  , lrowattr=ln_use_jattr )
      CALL iom_get( inum, jpdom_data, 'e1v'  , pe1v  , lrowattr=ln_use_jattr )
      CALL iom_get( inum, jpdom_data, 'e1f'  , pe1f  , lrowattr=ln_use_jattr )
      !
      CALL iom_get( inum, jpdom_data, 'e2t'  , pe2t  , lrowattr=ln_use_jattr )
      CALL iom_get( inum, jpdom_data, 'e2u'  , pe2u  , lrowattr=ln_use_jattr )
      CALL iom_get( inum, jpdom_data, 'e2v'  , pe2v  , lrowattr=ln_use_jattr )
      CALL iom_get( inum, jpdom_data, 'e2f'  , pe2f  , lrowattr=ln_use_jattr )
      !
      IF(  iom_varid( inum, 'ff_f', ldstop = .FALSE. ) > 0  .AND.  &
         & iom_varid( inum, 'ff_t', ldstop = .FALSE. ) > 0    ) THEN
         IF(lwp) WRITE(numout,*) '           Coriolis factor at f- and t-points read in ', TRIM( cn_domcfg ), ' file'
         CALL iom_get( inum, jpdom_data, 'ff_f'  , pff_f  , lrowattr=ln_use_jattr )
         CALL iom_get( inum, jpdom_data, 'ff_t'  , pff_t  , lrowattr=ln_use_jattr )
         kff = 1
      ELSE
         kff = 0
      ENDIF
      !
      IF( iom_varid( inum, 'e1e2u', ldstop = .FALSE. ) > 0 ) THEN
         IF(lwp) WRITE(numout,*) '           e1e2u & e1e2v read in ', TRIM( cn_domcfg ), ' file'
         CALL iom_get( inum, jpdom_data, 'e1e2u'  , pe1e2u  , lrowattr=ln_use_jattr )
         CALL iom_get( inum, jpdom_data, 'e1e2v'  , pe1e2v  , lrowattr=ln_use_jattr )
         ke1e2u_v = 1
      ELSE
         ke1e2u_v = 0
      ENDIF
      !
      CALL iom_close( inum )
      !
   END SUBROUTINE hgr_read
    
   !!======================================================================
END MODULE domhgr
