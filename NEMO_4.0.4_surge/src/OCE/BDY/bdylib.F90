MODULE bdylib
   !!======================================================================
   !!                       ***  MODULE  bdylib  ***
   !! Unstructured Open Boundary Cond. :  Library module of generic boundary algorithms.
   !!======================================================================
   !! History :  3.6  !  2013     (D. Storkey) original code
   !!            4.0  !  2014     (T. Lovato) Generalize OBC structure
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   bdy_orlanski_2d
   !!   bdy_orlanski_3d
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers 
   USE dom_oce        ! ocean space and time domain
   USE bdy_oce        ! ocean open boundary conditions
   USE phycst         ! physical constants
   USE bdyini
   !
   USE in_out_manager !
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp, ONLY: ctl_stop

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_frs, bdy_spe, bdy_nmn, bdy_orl
   PUBLIC   bdy_orlanski_2d
   PUBLIC   bdy_orlanski_3d

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_frs( idx, pta, dta )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_frs  ***
      !!
      !! ** Purpose : Apply the Flow Relaxation Scheme for tracers at open boundaries.
      !!
      !! Reference : Engedahl H., 1995, Tellus, 365-382.
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),                     INTENT(in) ::   idx  ! OBC indices
      REAL(wp), DIMENSION(:,:),            INTENT(in) ::   dta  ! OBC external data
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pta  ! tracer trend
      !!
      REAL(wp) ::   zwgt           ! boundary weight
      INTEGER  ::   ib, ik, igrd   ! dummy loop indices
      INTEGER  ::   ii, ij         ! 2D addresses
      !!----------------------------------------------------------------------
      !
      igrd = 1                       ! Everything is at T-points here
      DO ib = 1, idx%nblen(igrd)
         DO ik = 1, jpkm1
            ii = idx%nbi(ib,igrd) 
            ij = idx%nbj(ib,igrd)
            zwgt = idx%nbw(ib,igrd)
            pta(ii,ij,ik) = ( pta(ii,ij,ik) + zwgt * (dta(ib,ik) - pta(ii,ij,ik) ) ) * tmask(ii,ij,ik)
         END DO
      END DO
      !
   END SUBROUTINE bdy_frs


   SUBROUTINE bdy_spe( idx, pta, dta )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_spe  ***
      !!
      !! ** Purpose : Apply a specified value for tracers at open boundaries.
      !!
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),                     INTENT(in) ::   idx  ! OBC indices
      REAL(wp), DIMENSION(:,:),            INTENT(in) ::   dta  ! OBC external data
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pta  ! tracer trend
      !!
      INTEGER  ::   ib, ik, igrd   ! dummy loop indices
      INTEGER  ::   ii, ij         ! 2D addresses
      !!----------------------------------------------------------------------
      !
      igrd = 1                       ! Everything is at T-points here
      DO ib = 1, idx%nblenrim(igrd)
         ii = idx%nbi(ib,igrd)
         ij = idx%nbj(ib,igrd)
         DO ik = 1, jpkm1
            pta(ii,ij,ik) = dta(ib,ik) * tmask(ii,ij,ik)
         END DO
      END DO
      !
   END SUBROUTINE bdy_spe


   SUBROUTINE bdy_orl( idx, ptb, pta, dta, lrim0, ll_npo )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_orl  ***
      !!
      !! ** Purpose : Apply Orlanski radiation for tracers at open boundaries.
      !!              This is a wrapper routine for bdy_orlanski_3d below
      !!
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),                     INTENT(in) ::   idx  ! OBC indices
      REAL(wp), DIMENSION(:,:),            INTENT(in) ::   dta  ! OBC external data
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   ptb  ! before tracer field
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pta  ! tracer trend
      LOGICAL                 , OPTIONAL,  INTENT(in) ::   lrim0   ! indicate if rim 0 is treated
      LOGICAL,                             INTENT(in) ::   ll_npo  ! switch for NPO version
      !!
      INTEGER  ::   igrd                                    ! grid index
      !!----------------------------------------------------------------------
      !
      igrd = 1                       ! Everything is at T-points here
      !
      CALL bdy_orlanski_3d( idx, igrd, ptb(:,:,:), pta(:,:,:), dta, lrim0, ll_npo )
      !
   END SUBROUTINE bdy_orl


   SUBROUTINE bdy_orlanski_2d( idx, igrd, phib, phia, phi_ext, lrim0, ll_npo )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_orlanski_2d  ***
      !!             
      !!              - Apply Orlanski radiation condition adaptively to 2D fields:
      !!                  - radiation plus weak nudging at outflow points
      !!                  - no radiation and strong nudging at inflow points
      !! 
      !!
      !! References:  Marchesiello, McWilliams and Shchepetkin, Ocean Modelling vol. 3 (2001)    
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),          INTENT(in   ) ::   idx      ! BDY indices
      INTEGER ,                 INTENT(in   ) ::   igrd     ! grid index
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   phib     ! model before 2D field
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   phia     ! model after 2D field (to be updated)
      REAL(wp), DIMENSION(:)  , INTENT(in   ) ::   phi_ext  ! external forcing data
      LOGICAL, OPTIONAL,        INTENT(in   ) ::   lrim0    ! indicate if rim 0 is treated
      LOGICAL ,                 INTENT(in   ) ::   ll_npo   ! switch for NPO version
      !
      INTEGER  ::   jb                                     ! dummy loop indices
      INTEGER  ::   ii, ij, iibm1, iibm2, ijbm1, ijbm2     ! 2D addresses
      INTEGER  ::   iijm1, iijp1, ijjm1, ijjp1             ! 2D addresses
      INTEGER  ::   iibm1jp1, iibm1jm1, ijbm1jp1, ijbm1jm1 ! 2D addresses
      INTEGER  ::   ii_offset, ij_offset                   ! offsets for mask indices
      INTEGER  ::   flagu, flagv                           ! short cuts
      INTEGER  ::   ibeg, iend                             ! length of rim to be treated (rim 0 or rim 1 or both)
      REAL(wp) ::   zmask_x, zmask_y1, zmask_y2
      REAL(wp) ::   zex1, zex2, zey, zey1, zey2
      REAL(wp) ::   zdt, zdx, zdy, znor2, zrx, zry         ! intermediate calculations
      REAL(wp) ::   zout, zwgt, zdy_centred
      REAL(wp) ::   zdy_1, zdy_2, zsign_ups
      REAL(wp), PARAMETER :: zepsilon = 1.e-30                 ! local small value
      REAL(wp), POINTER, DIMENSION(:,:)          :: zmask      ! land/sea mask for field
      REAL(wp), POINTER, DIMENSION(:,:)          :: zmask_xdif ! land/sea mask for x-derivatives
      REAL(wp), POINTER, DIMENSION(:,:)          :: zmask_ydif ! land/sea mask for y-derivatives
      REAL(wp), POINTER, DIMENSION(:,:)          :: pe_xdif    ! scale factors for x-derivatives
      REAL(wp), POINTER, DIMENSION(:,:)          :: pe_ydif    ! scale factors for y-derivatives
      !!----------------------------------------------------------------------
      !
      ! ----------------------------------!
      ! Orlanski boundary conditions     :!
      ! ----------------------------------! 
     
      SELECT CASE(igrd)
         CASE(1)
            zmask      => tmask(:,:,1)
            zmask_xdif => umask(:,:,1)
            zmask_ydif => vmask(:,:,1)
            pe_xdif    => e1u(:,:)
            pe_ydif    => e2v(:,:)
            ii_offset = 0
            ij_offset = 0
         CASE(2)
            zmask      => umask(:,:,1)
            zmask_xdif => tmask(:,:,1)
            zmask_ydif => fmask(:,:,1)
            pe_xdif    => e1t(:,:)
            pe_ydif    => e2f(:,:)
            ii_offset = 1
            ij_offset = 0
         CASE(3)
            zmask      => vmask(:,:,1)
            zmask_xdif => fmask(:,:,1)
            zmask_ydif => tmask(:,:,1)
            pe_xdif    => e1f(:,:)
            pe_ydif    => e2t(:,:)
            ii_offset = 0
            ij_offset = 1
         CASE DEFAULT ;   CALL ctl_stop( 'unrecognised value for igrd in bdy_orlanksi_2d' )
      END SELECT
      !
      IF( PRESENT(lrim0) ) THEN
         IF( lrim0 ) THEN   ;   ibeg = 1                       ;   iend = idx%nblenrim0(igrd)   ! rim 0
         ELSE               ;   ibeg = idx%nblenrim0(igrd)+1   ;   iend = idx%nblenrim(igrd)    ! rim 1
         END IF
      ELSE                  ;   ibeg = 1                       ;   iend = idx%nblenrim(igrd)    ! both
      END IF
      !
      DO jb = ibeg, iend
         ii  = idx%nbi(jb,igrd)
         ij  = idx%nbj(jb,igrd) 
         IF( ii == 1 .OR. ii == jpi .OR. ij == 1 .OR. ij == jpj )   CYCLE
         flagu = int( idx%flagu(jb,igrd) )
         flagv = int( idx%flagv(jb,igrd) )
         !
         ! Calculate positions of b-1 and b-2 points for this rim point
         ! also (b-1,j-1) and (b-1,j+1) points
         iibm1 = ii + flagu ; iibm2 = ii + 2*flagu 
         ijbm1 = ij + flagv ; ijbm2 = ij + 2*flagv
          !
         iijm1 = ii - abs(flagv) ; iijp1 = ii + abs(flagv) 
         ijjm1 = ij - abs(flagu) ; ijjp1 = ij + abs(flagu)
         !
         iibm1jm1 = ii + flagu - abs(flagv) ; iibm1jp1 = ii + flagu + abs(flagv) 
         ijbm1jm1 = ij + flagv - abs(flagu) ; ijbm1jp1 = ij + flagv + abs(flagu) 
         !
         ! Calculate scale factors for calculation of spatial derivatives.        
         zex1 = ( abs(iibm1-iibm2) * pe_xdif(iibm1   +ii_offset,ijbm1             )   &
        &       + abs(ijbm1-ijbm2) * pe_ydif(iibm1             ,ijbm1   +ij_offset) ) 
         zex2 = ( abs(iibm1-iibm2) * pe_xdif(iibm2   +ii_offset,ijbm2             )   &
        &       + abs(ijbm1-ijbm2) * pe_ydif(iibm2             ,ijbm2   +ij_offset) ) 
         zey1 = ( (iibm1-iibm1jm1) * pe_xdif(iibm1jm1+ii_offset,ijbm1jm1          )   & 
        &      +  (ijbm1-ijbm1jm1) * pe_ydif(iibm1jm1          ,ijbm1jm1+ij_offset) ) 
         zey2 = ( (iibm1jp1-iibm1) * pe_xdif(iibm1   +ii_offset,ijbm1             )   &
        &      +  (ijbm1jp1-ijbm1) * pe_ydif(iibm1             ,ijbm1   +ij_offset) ) 
         ! make sure scale factors are nonzero
         if( zey1 .lt. rsmall ) zey1 = zey2
         if( zey2 .lt. rsmall ) zey2 = zey1
         zex1 = max(zex1,rsmall); zex2 = max(zex2,rsmall)
         zey1 = max(zey1,rsmall); zey2 = max(zey2,rsmall); 
         !
         ! Calculate masks for calculation of spatial derivatives.
         zmask_x  = ( abs(iibm1-iibm2) * zmask_xdif(iibm2   +ii_offset,ijbm2               )   &
        &           + abs(ijbm1-ijbm2) * zmask_ydif(iibm2             ,ijbm2   +ij_offset) ) 
         zmask_y1 = ( (iibm1-iibm1jm1) * zmask_xdif(iibm1jm1+ii_offset,ijbm1jm1            )   & 
        &          +  (ijbm1-ijbm1jm1) * zmask_ydif(iibm1jm1          ,ijbm1jm1+ij_offset) ) 
         zmask_y2 = ( (iibm1jp1-iibm1) * zmask_xdif(iibm1   +ii_offset,ijbm1               )   &
        &          +  (ijbm1jp1-ijbm1) * zmask_ydif(iibm1             ,ijbm1   +ij_offset) ) 

         ! Calculation of terms required for both versions of the scheme. 
         ! Mask derivatives to ensure correct land boundary conditions for each variable.
         ! Centred derivative is calculated as average of "left" and "right" derivatives for 
         ! this reason. 
         ! Note no rdt factor in expression for zdt because it cancels in the expressions for 
         ! zrx and zry.
         zdt   =     phia(iibm1   ,ijbm1   ) - phib(iibm1   ,ijbm1   )
         zdx   = ( ( phia(iibm1   ,ijbm1   ) - phia(iibm2   ,ijbm2   ) ) / zex2 ) * zmask_x 
         zdy_1 = ( ( phib(iibm1   ,ijbm1   ) - phib(iibm1jm1,ijbm1jm1) ) / zey1 ) * zmask_y1    
         zdy_2 = ( ( phib(iibm1jp1,ijbm1jp1) - phib(iibm1   ,ijbm1   ) ) / zey2 ) * zmask_y2 
         zdy_centred = 0.5 * ( zdy_1 + zdy_2 )
!!$         zdy_centred = phib(iibm1jp1,ijbm1jp1) - phib(iibm1jm1,ijbm1jm1)
         ! upstream differencing for tangential derivatives
         zsign_ups = sign( 1., zdt * zdy_centred )
         zsign_ups = 0.5*( zsign_ups + abs(zsign_ups) )
         zdy = zsign_ups * zdy_1 + (1. - zsign_ups) * zdy_2
         znor2 = zdx * zdx + zdy * zdy
         znor2 = max(znor2,zepsilon)
         !
         zrx = zdt * zdx / ( zex1 * znor2 ) 
!!$         zrx = min(zrx,2.0_wp)
         zout = sign( 1., zrx )
         zout = 0.5*( zout + abs(zout) )
         zwgt = 2.*rdt*( (1.-zout) * idx%nbd(jb,igrd) + zout * idx%nbdout(jb,igrd) )
         ! only apply radiation on outflow points 
         if( ll_npo ) then     !! NPO version !!
            phia(ii,ij) = (1.-zout) * ( phib(ii,ij) + zwgt * ( phi_ext(jb) - phib(ii,ij) ) )        &
           &            + zout      * ( phib(ii,ij) + zrx*phia(iibm1,ijbm1)                         &
           &                            + zwgt * ( phi_ext(jb) - phib(ii,ij) ) ) / ( 1. + zrx ) 
         else                  !! full oblique radiation !!
            zsign_ups = sign( 1., zdt * zdy )
            zsign_ups = 0.5*( zsign_ups + abs(zsign_ups) )
            zey = zsign_ups * zey1 + (1.-zsign_ups) * zey2 
            zry = zdt * zdy / ( zey * znor2 ) 
            phia(ii,ij) = (1.-zout) * ( phib(ii,ij) + zwgt * ( phi_ext(jb) - phib(ii,ij) ) )        &
           &            + zout      * ( phib(ii,ij) + zrx*phia(iibm1,ijbm1)                         &
           &                    - zsign_ups      * zry * ( phib(ii   ,ij   ) - phib(iijm1,ijjm1 ) ) &
           &                    - (1.-zsign_ups) * zry * ( phib(iijp1,ijjp1) - phib(ii   ,ij    ) ) &
           &                    + zwgt * ( phi_ext(jb) - phib(ii,ij) ) ) / ( 1. + zrx ) 
         end if
         phia(ii,ij) = phia(ii,ij) * zmask(ii,ij)
      END DO
      !
   END SUBROUTINE bdy_orlanski_2d


   SUBROUTINE bdy_orlanski_3d( idx, igrd, phib, phia, phi_ext, lrim0, ll_npo )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_orlanski_3d  ***
      !!             
      !!              - Apply Orlanski radiation condition adaptively to 3D fields:
      !!                  - radiation plus weak nudging at outflow points
      !!                  - no radiation and strong nudging at inflow points
      !! 
      !!
      !! References:  Marchesiello, McWilliams and Shchepetkin, Ocean Modelling vol. 3 (2001)    
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),            INTENT(in   ) ::   idx      ! BDY indices
      INTEGER ,                   INTENT(in   ) ::   igrd     ! grid index
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   phib     ! model before 3D field
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   phia     ! model after 3D field (to be updated)
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   phi_ext  ! external forcing data
      LOGICAL, OPTIONAL,          INTENT(in   ) ::   lrim0    ! indicate if rim 0 is treated
      LOGICAL ,                   INTENT(in   ) ::   ll_npo   ! switch for NPO version
      !
      INTEGER  ::   jb, jk                                 ! dummy loop indices
      INTEGER  ::   ii, ij, iibm1, iibm2, ijbm1, ijbm2     ! 2D addresses
      INTEGER  ::   iijm1, iijp1, ijjm1, ijjp1             ! 2D addresses
      INTEGER  ::   iibm1jp1, iibm1jm1, ijbm1jp1, ijbm1jm1 ! 2D addresses
      INTEGER  ::   ii_offset, ij_offset                   ! offsets for mask indices
      INTEGER  ::   flagu, flagv                           ! short cuts
      INTEGER  ::   ibeg, iend                             ! length of rim to be treated (rim 0 or rim 1 or both)
      REAL(wp) ::   zmask_x, zmask_y1, zmask_y2
      REAL(wp) ::   zex1, zex2, zey, zey1, zey2
      REAL(wp) ::   zdt, zdx, zdy, znor2, zrx, zry         ! intermediate calculations
      REAL(wp) ::   zout, zwgt, zdy_centred
      REAL(wp) ::   zdy_1, zdy_2,  zsign_ups
      REAL(wp), PARAMETER :: zepsilon = 1.e-30                 ! local small value
      REAL(wp), POINTER, DIMENSION(:,:,:)        :: zmask      ! land/sea mask for field
      REAL(wp), POINTER, DIMENSION(:,:,:)        :: zmask_xdif ! land/sea mask for x-derivatives
      REAL(wp), POINTER, DIMENSION(:,:,:)        :: zmask_ydif ! land/sea mask for y-derivatives
      REAL(wp), POINTER, DIMENSION(:,:)          :: pe_xdif    ! scale factors for x-derivatives
      REAL(wp), POINTER, DIMENSION(:,:)          :: pe_ydif    ! scale factors for y-derivatives
      !!----------------------------------------------------------------------
      !
      ! ----------------------------------!
      ! Orlanski boundary conditions     :!
      ! ----------------------------------! 
      !
      SELECT CASE(igrd)
         CASE(1)
            zmask      => tmask(:,:,:)
            zmask_xdif => umask(:,:,:)
            zmask_ydif => vmask(:,:,:)
            pe_xdif    => e1u(:,:)
            pe_ydif    => e2v(:,:)
            ii_offset = 0
            ij_offset = 0
         CASE(2)
            zmask      => umask(:,:,:)
            zmask_xdif => tmask(:,:,:)
            zmask_ydif => fmask(:,:,:)
            pe_xdif    => e1t(:,:)
            pe_ydif    => e2f(:,:)
            ii_offset = 1
            ij_offset = 0
         CASE(3)
            zmask      => vmask(:,:,:)
            zmask_xdif => fmask(:,:,:)
            zmask_ydif => tmask(:,:,:)
            pe_xdif    => e1f(:,:)
            pe_ydif    => e2t(:,:)
            ii_offset = 0
            ij_offset = 1
         CASE DEFAULT ;   CALL ctl_stop( 'unrecognised value for igrd in bdy_orlanksi_2d' )
      END SELECT
      !
      IF( PRESENT(lrim0) ) THEN
         IF( lrim0 ) THEN   ;   ibeg = 1                       ;   iend = idx%nblenrim0(igrd)   ! rim 0
         ELSE               ;   ibeg = idx%nblenrim0(igrd)+1   ;   iend = idx%nblenrim(igrd)    ! rim 1
         END IF
      ELSE                  ;   ibeg = 1                       ;   iend = idx%nblenrim(igrd)    ! both
      END IF
      !
      DO jk = 1, jpk
         !            
         DO jb = ibeg, iend
            ii  = idx%nbi(jb,igrd)
            ij  = idx%nbj(jb,igrd) 
            IF( ii == 1 .OR. ii == jpi .OR. ij == 1 .OR. ij == jpj )   CYCLE
            flagu = int( idx%flagu(jb,igrd) )
            flagv = int( idx%flagv(jb,igrd) )
            !
            ! calculate positions of b-1 and b-2 points for this rim point
            ! also (b-1,j-1) and (b-1,j+1) points
            iibm1 = ii + flagu ; iibm2 = ii + 2*flagu 
            ijbm1 = ij + flagv ; ijbm2 = ij + 2*flagv
            !
            iijm1 = ii - abs(flagv) ; iijp1 = ii + abs(flagv) 
            ijjm1 = ij - abs(flagu) ; ijjp1 = ij + abs(flagu)
            !
            iibm1jm1 = ii + flagu - abs(flagv) ; iibm1jp1 = ii + flagu + abs(flagv) 
            ijbm1jm1 = ij + flagv - abs(flagu) ; ijbm1jp1 = ij + flagv + abs(flagu) 
            !
            ! Calculate scale factors for calculation of spatial derivatives.        
            zex1 = ( abs(iibm1-iibm2) * pe_xdif(iibm1   +ii_offset,ijbm1             )   &
           &       + abs(ijbm1-ijbm2) * pe_ydif(iibm1             ,ijbm1+ij_offset   ) ) 
            zex2 = ( abs(iibm1-iibm2) * pe_xdif(iibm2   +ii_offset,ijbm2             )   &
           &       + abs(ijbm1-ijbm2) * pe_ydif(iibm2             ,ijbm2+ij_offset   ) ) 
            zey1 = ( (iibm1-iibm1jm1) * pe_xdif(iibm1jm1+ii_offset,ijbm1jm1          )   & 
           &      +  (ijbm1-ijbm1jm1) * pe_ydif(iibm1jm1          ,ijbm1jm1+ij_offset) ) 
            zey2 = ( (iibm1jp1-iibm1) * pe_xdif(iibm1   +ii_offset,ijbm1             )   &
           &      +  (ijbm1jp1-ijbm1) * pe_ydif(iibm1             ,ijbm1+ij_offset   ) ) 
            ! make sure scale factors are nonzero
            if( zey1 .lt. rsmall ) zey1 = zey2
            if( zey2 .lt. rsmall ) zey2 = zey1
            zex1 = max(zex1,rsmall); zex2 = max(zex2,rsmall); 
            zey1 = max(zey1,rsmall); zey2 = max(zey2,rsmall); 
            !
            ! Calculate masks for calculation of spatial derivatives.        
            zmask_x  = ( abs(iibm1-iibm2) * zmask_xdif(iibm2   +ii_offset,ijbm2             ,jk)   &
           &           + abs(ijbm1-ijbm2) * zmask_ydif(iibm2             ,ijbm2   +ij_offset,jk) ) 
            zmask_y1 = ( (iibm1-iibm1jm1) * zmask_xdif(iibm1jm1+ii_offset,ijbm1jm1          ,jk)   & 
           &          +  (ijbm1-ijbm1jm1) * zmask_ydif(iibm1jm1          ,ijbm1jm1+ij_offset,jk) ) 
            zmask_y2 = ( (iibm1jp1-iibm1) * zmask_xdif(iibm1   +ii_offset,ijbm1             ,jk)   &
           &          +  (ijbm1jp1-ijbm1) * zmask_ydif(iibm1             ,ijbm1   +ij_offset,jk) ) 
            !
            ! Calculate normal (zrx) and tangential (zry) components of radiation velocities.
            ! Mask derivatives to ensure correct land boundary conditions for each variable.
            ! Centred derivative is calculated as average of "left" and "right" derivatives for 
            ! this reason. 
            zdt   =     phia(iibm1   ,ijbm1   ,jk) - phib(iibm1   ,ijbm1   ,jk)
            zdx   = ( ( phia(iibm1   ,ijbm1   ,jk) - phia(iibm2   ,ijbm2   ,jk) ) / zex2 ) * zmask_x                  
            zdy_1 = ( ( phib(iibm1   ,ijbm1   ,jk) - phib(iibm1jm1,ijbm1jm1,jk) ) / zey1 ) * zmask_y1  
            zdy_2 = ( ( phib(iibm1jp1,ijbm1jp1,jk) - phib(iibm1   ,ijbm1   ,jk) ) / zey2 ) * zmask_y2      
            zdy_centred = 0.5 * ( zdy_1 + zdy_2 )
!!$            zdy_centred = phib(iibm1jp1,ijbm1jp1,jk) - phib(iibm1jm1,ijbm1jm1,jk)
            ! upstream differencing for tangential derivatives
            zsign_ups = sign( 1., zdt * zdy_centred )
            zsign_ups = 0.5*( zsign_ups + abs(zsign_ups) )
            zdy = zsign_ups * zdy_1 + (1. - zsign_ups) * zdy_2
            znor2 = zdx * zdx + zdy * zdy
            znor2 = max(znor2,zepsilon)
            !
            ! update boundary value:
            zrx = zdt * zdx / ( zex1 * znor2 )
!!$            zrx = min(zrx,2.0_wp)
            zout = sign( 1., zrx )
            zout = 0.5*( zout + abs(zout) )
            zwgt = 2.*rdt*( (1.-zout) * idx%nbd(jb,igrd) + zout * idx%nbdout(jb,igrd) )
            ! only apply radiation on outflow points 
            if( ll_npo ) then     !! NPO version !!
               phia(ii,ij,jk) = (1.-zout) * ( phib(ii,ij,jk) + zwgt * ( phi_ext(jb,jk) - phib(ii,ij,jk) ) ) &
              &               + zout      * ( phib(ii,ij,jk) + zrx*phia(iibm1,ijbm1,jk)                     &
              &                            + zwgt * ( phi_ext(jb,jk) - phib(ii,ij,jk) ) ) / ( 1. + zrx ) 
            else                  !! full oblique radiation !!
               zsign_ups = sign( 1., zdt * zdy )
               zsign_ups = 0.5*( zsign_ups + abs(zsign_ups) )
               zey = zsign_ups * zey1 + (1.-zsign_ups) * zey2 
               zry = zdt * zdy / ( zey * znor2 ) 
               phia(ii,ij,jk) = (1.-zout) * ( phib(ii,ij,jk) + zwgt * ( phi_ext(jb,jk) - phib(ii,ij,jk) ) )    &
              &               + zout      * ( phib(ii,ij,jk) + zrx*phia(iibm1,ijbm1,jk)                        &
              &                       - zsign_ups      * zry * ( phib(ii   ,ij   ,jk) - phib(iijm1,ijjm1,jk) ) &
              &                       - (1.-zsign_ups) * zry * ( phib(iijp1,ijjp1,jk) - phib(ii   ,ij   ,jk) ) &
              &                       + zwgt * ( phi_ext(jb,jk) - phib(ii,ij,jk) ) ) / ( 1. + zrx ) 
            end if
            phia(ii,ij,jk) = phia(ii,ij,jk) * zmask(ii,ij,jk)
         END DO
         !
      END DO
      !
   END SUBROUTINE bdy_orlanski_3d

   SUBROUTINE bdy_nmn( idx, igrd, phia, lrim0 )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_nmn  ***
      !!                    
      !! ** Purpose : Duplicate the value at open boundaries, zero gradient.
      !! 
      !!
      !! ** Method  : - take the average of free ocean neighbours
      !!
      !!      ___   !   |_____|   !   ___|    !   __|x o   !   |_   _|     ! |      
      !!   __|x     !      x      !     x o   !      o     !     |_|       ! |x o   
      !!      o     !      o      !     o     !            !    o x o      ! |x_x_ 
      !!                                                   !      o      
      !!----------------------------------------------------------------------
      INTEGER,                    INTENT(in   )  ::   igrd     ! grid index
      REAL(wp), DIMENSION(:,:,:), INTENT(inout)  ::   phia     ! model after 3D field (to be updated), must be masked
      TYPE(OBC_INDEX),            INTENT(in   )  ::   idx      ! OBC indices
      LOGICAL, OPTIONAL,          INTENT(in   )  ::   lrim0    ! indicate if rim 0 is treated
      !! 
      REAL(wp) ::   zweight
      REAL(wp), POINTER, DIMENSION(:,:,:)      :: zmask         ! land/sea mask for field
      INTEGER  ::   ib, ik   ! dummy loop indices
      INTEGER  ::   ii, ij   ! 2D addresses
      INTEGER  ::   ipkm1    ! size of phia third dimension minus 1
      INTEGER  ::   ibeg, iend                          ! length of rim to be treated (rim 0 or rim 1 or both)
      INTEGER  ::   ii1, ii2, ii3, ij1, ij2, ij3, itreat
      !!----------------------------------------------------------------------
      !
      ipkm1 = MAX( SIZE(phia,3) - 1, 1 ) 
      !
      SELECT CASE(igrd)
         CASE(1)   ;   zmask => tmask(:,:,:)
         CASE(2)   ;   zmask => umask(:,:,:)
         CASE(3)   ;   zmask => vmask(:,:,:)
         CASE DEFAULT ;   CALL ctl_stop( 'unrecognised value for igrd in bdy_nmn' )
      END SELECT
      !
      IF( PRESENT(lrim0) ) THEN
         IF( lrim0 ) THEN   ;   ibeg = 1                       ;   iend = idx%nblenrim0(igrd)   ! rim 0
         ELSE               ;   ibeg = idx%nblenrim0(igrd)+1   ;   iend = idx%nblenrim(igrd)    ! rim 1
         END IF
      ELSE                  ;   ibeg = 1                       ;   iend = idx%nblenrim(igrd)    ! both
      END IF
      !
      DO ib = ibeg, iend
         ii = idx%nbi(ib,igrd)
         ij = idx%nbj(ib,igrd)
         itreat = idx%ntreat(ib,igrd)
         CALL find_neib( ii, ij, itreat, ii1, ij1, ii2, ij2, ii3, ij3 )   ! find free ocean neighbours
         SELECT CASE( itreat )
         CASE( 1:8 )
            IF( ii1 < 1 .OR. ii1 > jpi .OR. ij1 < 1 .OR. ij1 > jpj )   CYCLE
            DO ik = 1, ipkm1
               IF( zmask(ii1,ij1,ik) /= 0. )   phia(ii,ij,ik) = phia(ii1,ij1,ik)  
            END DO
         CASE( 9:12 )
            IF( ii1 < 1 .OR. ii1 > jpi .OR. ij1 < 1 .OR. ij1 > jpj )   CYCLE
            IF( ii2 < 1 .OR. ii2 > jpi .OR. ij2 < 1 .OR. ij2 > jpj )   CYCLE
            DO ik = 1, ipkm1
               zweight = zmask(ii1,ij1,ik) + zmask(ii2,ij2,ik)
               IF( zweight /= 0. )   phia(ii,ij,ik) = ( phia(ii1,ij1,ik) + phia(ii2,ij2,ik) ) / zweight
            END DO
         CASE( 13:16 )
            IF( ii1 < 1 .OR. ii1 > jpi .OR. ij1 < 1 .OR. ij1 > jpj )   CYCLE
            IF( ii2 < 1 .OR. ii2 > jpi .OR. ij2 < 1 .OR. ij2 > jpj )   CYCLE
            IF( ii3 < 1 .OR. ii3 > jpi .OR. ij3 < 1 .OR. ij3 > jpj )   CYCLE
            DO ik = 1, ipkm1
               zweight = zmask(ii1,ij1,ik) + zmask(ii2,ij2,ik) + zmask(ii3,ij3,ik)
               IF( zweight /= 0. )   phia(ii,ij,ik) = ( phia(ii1,ij1,ik) + phia(ii2,ij2,ik) + phia(ii3,ij3,ik) ) / zweight
            END DO
         END SELECT
      END DO
      !
   END SUBROUTINE bdy_nmn

   !!======================================================================
END MODULE bdylib
