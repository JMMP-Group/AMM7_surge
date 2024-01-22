MODULE bdydyn2d
   !!======================================================================
   !!                       ***  MODULE  bdydyn  ***
   !! Unstructured Open Boundary Cond. :   Apply boundary conditions to barotropic solution
   !!======================================================================
   !! History :  3.4  !  2011     (D. Storkey) new module as part of BDY rewrite
   !!            3.5  !  2012     (S. Mocavero, I. Epicoco) Optimization of BDY communications
   !!            3.5  !  2013-07  (J. Chanut) Compliant with time splitting changes
   !!----------------------------------------------------------------------
   !!   bdy_dyn2d          : Apply open boundary conditions to barotropic variables.
   !!   bdy_dyn2d_frs      : Apply Flow Relaxation Scheme 
   !!   bdy_dyn2d_fla      : Apply Flather condition
   !!   bdy_dyn2d_orlanski : Orlanski Radiation
   !!   bdy_ssh            : Duplicate sea level across open boundaries
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE bdy_oce         ! ocean open boundary conditions
   USE bdylib          ! BDY library routines
   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE wet_dry         ! Use wet dry to get reference ssh level
   USE in_out_manager  !
   USE lib_mpp, ONLY: ctl_stop

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dyn2d   ! routine called in dynspg_ts and bdy_dyn
   PUBLIC   bdy_ssh       ! routine called in dynspg_ts or sshwzv

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_dyn2d( kt, pua2d, pva2d, pub2d, pvb2d, phur, phvr, pssh  )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn2d  ***
      !!
      !! ** Purpose : - Apply open boundary conditions for barotropic variables
      !!
      !!----------------------------------------------------------------------
      INTEGER,                      INTENT(in) ::   kt   ! Main time step counter
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pua2d, pva2d 
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pub2d, pvb2d
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: phur, phvr
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pssh
      !!
      INTEGER  ::   ib_bdy, ir     ! BDY set index, rim index
      LOGICAL  ::   llrim0         ! indicate if rim 0 is treated
      LOGICAL, DIMENSION(4) :: llsend2, llrecv2, llsend3, llrecv3  ! indicate how communications are to be carried out
      
      llsend2(:) = .false.   ;   llrecv2(:) = .false.
      llsend3(:) = .false.   ;   llrecv3(:) = .false.
      DO ir = 1, 0, -1   ! treat rim 1 before rim 0
         IF( ir == 0 ) THEN   ;   llrim0 = .TRUE.
         ELSE                 ;   llrim0 = .FALSE.
         END IF
         DO ib_bdy=1, nb_bdy
            SELECT CASE( cn_dyn2d(ib_bdy) )
            CASE('none')
               CYCLE
            CASE('frs')   ! treat the whole boundary at once
               IF( llrim0 )   CALL bdy_dyn2d_frs( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, pua2d, pva2d )
            CASE('flather')
               CALL bdy_dyn2d_fla( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, pua2d, pva2d, pssh, phur, phvr, llrim0 )
            CASE('orlanski')
               CALL bdy_dyn2d_orlanski( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, &
                    & pua2d, pva2d, pub2d, pvb2d, llrim0, ll_npo=.false. )
            CASE('orlanski_npo')
               CALL bdy_dyn2d_orlanski( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, &
                    & pua2d, pva2d, pub2d, pvb2d, llrim0, ll_npo=.true.  )
            CASE DEFAULT
               CALL ctl_stop( 'bdy_dyn2d : unrecognised option for open boundaries for barotropic variables' )
            END SELECT
         ENDDO
         !
         IF( nn_hls > 1 .AND. ir == 1 ) CYCLE   ! at least 2 halos will be corrected -> no need to correct rim 1 before rim 0
         IF( nn_hls == 1 ) THEN
            llsend2(:) = .false.   ;   llrecv2(:) = .false.
            llsend3(:) = .false.   ;   llrecv3(:) = .false.
         END IF
         DO ib_bdy=1, nb_bdy
            SELECT CASE( cn_dyn2d(ib_bdy) )
            CASE('flather')
               llsend2(1:2) = llsend2(1:2) .OR. lsend_bdyint(ib_bdy,2,1:2,ir)   ! west/east, U points
               llsend2(1)   = llsend2(1)   .OR. lsend_bdyext(ib_bdy,2,1,ir)     ! neighbour might search point towards its east bdy
               llrecv2(1:2) = llrecv2(1:2) .OR. lrecv_bdyint(ib_bdy,2,1:2,ir)   ! west/east, U points
               llrecv2(2)   = llrecv2(2)   .OR. lrecv_bdyext(ib_bdy,2,2,ir)     ! might search point towards bdy on the east
               llsend3(3:4) = llsend3(3:4) .OR. lsend_bdyint(ib_bdy,3,3:4,ir)   ! north/south, V points
               llsend3(3)   = llsend3(3)   .OR. lsend_bdyext(ib_bdy,3,3,ir)     ! neighbour might search point towards its north bdy 
               llrecv3(3:4) = llrecv3(3:4) .OR. lrecv_bdyint(ib_bdy,3,3:4,ir)   ! north/south, V points
               llrecv3(4)   = llrecv3(4)   .OR. lrecv_bdyext(ib_bdy,3,4,ir)     ! might search point towards bdy on the north
            CASE('orlanski', 'orlanski_npo')
               llsend2(:) = llsend2(:) .OR. lsend_bdy(ib_bdy,2,:,ir)   ! possibly every direction, U points
               llrecv2(:) = llrecv2(:) .OR. lrecv_bdy(ib_bdy,2,:,ir)   ! possibly every direction, U points
               llsend3(:) = llsend3(:) .OR. lsend_bdy(ib_bdy,3,:,ir)   ! possibly every direction, V points
               llrecv3(:) = llrecv3(:) .OR. lrecv_bdy(ib_bdy,3,:,ir)   ! possibly every direction, V points
            END SELECT
         END DO
         IF( ANY(llsend2) .OR. ANY(llrecv2) ) THEN   ! if need to send/recv in at least one direction
            CALL lbc_lnk( 'bdydyn2d', pua2d, 'U', -1., kfillmode=jpfillnothing ,lsend=llsend2, lrecv=llrecv2 )
         END IF
         IF( ANY(llsend3) .OR. ANY(llrecv3) ) THEN   ! if need to send/recv in at least one direction
            CALL lbc_lnk( 'bdydyn2d', pva2d, 'V', -1., kfillmode=jpfillnothing ,lsend=llsend3, lrecv=llrecv3 )
         END IF
         !
      END DO   ! ir
      !
   END SUBROUTINE bdy_dyn2d

   SUBROUTINE bdy_dyn2d_frs( idx, dta, ib_bdy, pua2d, pva2d )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn2d_frs  ***
      !!
      !! ** Purpose : - Apply the Flow Relaxation Scheme for barotropic velocities
      !!                at open boundaries.
      !!
      !! References :- Engedahl H., 1995: Use of the flow relaxation scheme in 
      !!               a three-dimensional baroclinic ocean model with realistic
      !!               topography. Tellus, 365-382.
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pua2d, pva2d 
      !!
      INTEGER  ::   jb             ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      igrd = 2                      ! Relaxation of zonal velocity
      DO jb = 1, idx%nblen(igrd)
         ii   = idx%nbi(jb,igrd)
         ij   = idx%nbj(jb,igrd)
         zwgt = idx%nbw(jb,igrd)
         pua2d(ii,ij) = ( pua2d(ii,ij) + zwgt * ( dta%u2d(jb) - pua2d(ii,ij) ) ) * umask(ii,ij,1)
      END DO
      !
      igrd = 3                      ! Relaxation of meridional velocity
      DO jb = 1, idx%nblen(igrd)
         ii   = idx%nbi(jb,igrd)
         ij   = idx%nbj(jb,igrd)
         zwgt = idx%nbw(jb,igrd)
         pva2d(ii,ij) = ( pva2d(ii,ij) + zwgt * ( dta%v2d(jb) - pva2d(ii,ij) ) ) * vmask(ii,ij,1)
      END DO 
      !
   END SUBROUTINE bdy_dyn2d_frs


   SUBROUTINE bdy_dyn2d_fla( idx, dta, ib_bdy, pua2d, pva2d, pssh, phur, phvr, llrim0 )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_dyn2d_fla  ***
      !!             
      !!              - Apply Flather boundary conditions on normal barotropic velocities 
      !!
      !! ** WARNINGS about FLATHER implementation:
      !!1. According to Palma and Matano, 1998 "after ssh" is used. 
      !!   In ROMS and POM implementations, it is "now ssh". In the current 
      !!   implementation (tested only in the EEL-R5 conf.), both cases were unstable. 
      !!   So I use "before ssh" in the following.
      !!
      !!2. We assume that the normal ssh gradient at the bdy is zero. As a matter of 
      !!   fact, the model ssh just inside the dynamical boundary is used (the outside  
      !!   ssh in the code is not updated).
      !!
      !! References:  Flather, R. A., 1976: A tidal model of the northwest European
      !!              continental shelf. Mem. Soc. R. Sci. Liege, Ser. 6,10, 141-164.     
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),              INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),               INTENT(in) ::   dta  ! OBC external data
      INTEGER,                      INTENT(in) ::   ib_bdy  ! BDY set index
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pua2d, pva2d
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pssh, phur, phvr
      LOGICAL                     , INTENT(in) ::   llrim0   ! indicate if rim 0 is treated
      INTEGER  ::   ibeg, iend                       ! length of rim to be treated (rim 0 or rim 1)
      INTEGER  ::   jb, igrd                         ! dummy loop indices
      INTEGER  ::   ii, ij                           ! 2D addresses
      INTEGER  ::   iiTrim, ijTrim                   ! T pts i/j-indice on the rim
      INTEGER  ::   iiToce, ijToce, iiUoce, ijVoce   ! T, U and V pts i/j-indice of the ocean next to the rim
      REAL(wp) ::   flagu, flagv                     ! short cuts
      REAL(wp) ::   zfla                             ! Flather correction
      REAL(wp) ::   z1_2                             ! 
      REAL(wp), DIMENSION(jpi,jpj) ::   sshdta       ! 2D version of dta%ssh
      !!----------------------------------------------------------------------

      z1_2 = 0.5_wp

      ! ---------------------------------!
      ! Flather boundary conditions     :!
      ! ---------------------------------!

      ! Fill temporary array with ssh data (here we use spgu with the alias sshdta):
      igrd = 1
      IF( llrim0 ) THEN   ;   ibeg = 1                       ;   iend = idx%nblenrim0(igrd)
      ELSE                ;   ibeg = idx%nblenrim0(igrd)+1   ;   iend = idx%nblenrim(igrd)
      END IF
      !
      DO jb = ibeg, iend
         ii = idx%nbi(jb,igrd)
         ij = idx%nbj(jb,igrd)
         IF( ll_wd ) THEN   ;   sshdta(ii, ij) = dta%ssh(jb) - ssh_ref 
         ELSE               ;   sshdta(ii, ij) = dta%ssh(jb)
         ENDIF
      END DO
      !
      igrd = 2      ! Flather bc on u-velocity
      !             ! remember that flagu=-1 if normal velocity direction is outward
      !             ! I think we should rather use after ssh ?
      IF( llrim0 ) THEN   ;   ibeg = 1                       ;   iend = idx%nblenrim0(igrd)
      ELSE                ;   ibeg = idx%nblenrim0(igrd)+1   ;   iend = idx%nblenrim(igrd)
      END IF
      DO jb = ibeg, iend
         ii    = idx%nbi(jb,igrd)
         ij    = idx%nbj(jb,igrd)
         flagu = idx%flagu(jb,igrd)
         IF( flagu == 0. ) THEN
            pua2d(ii,ij) = dta%u2d(jb)
         ELSE      ! T pts j-indice       on the rim          on the ocean next to the rim on T and U points
            IF( flagu ==  1. ) THEN   ;   iiTrim = ii     ;   iiToce = ii+1   ;   iiUoce = ii+1   ;   ENDIF
            IF( flagu == -1. ) THEN   ;   iiTrim = ii+1   ;   iiToce = ii     ;   iiUoce = ii-1   ;   ENDIF
            !
            ! Rare case : rim is parallel to the mpi subdomain border and on the halo : point will be received
            IF( iiTrim > jpi .OR. iiToce > jpi .OR. iiUoce > jpi .OR. iiUoce < 1 )   CYCLE   
            !
            zfla = dta%u2d(jb) - flagu * SQRT( grav * phur(ii, ij) ) * ( pssh(iiToce,ij) - sshdta(iiTrim,ij) )
            !
            ! jchanut tschanges, use characteristics method (Blayo et Debreu, 2005) :
            ! mix Flather scheme with velocity of the ocean next to the rim
            pua2d(ii,ij) =  z1_2 * ( pua2d(iiUoce,ij) + zfla )
         END IF
      END DO
      !
      igrd = 3      ! Flather bc on v-velocity
      !             ! remember that flagv=-1 if normal velocity direction is outward
      IF( llrim0 ) THEN   ;   ibeg = 1                       ;   iend = idx%nblenrim0(igrd)
      ELSE                ;   ibeg = idx%nblenrim0(igrd)+1   ;   iend = idx%nblenrim(igrd)
      END IF
      DO jb = ibeg, iend
         ii    = idx%nbi(jb,igrd)
         ij    = idx%nbj(jb,igrd)
         flagv = idx%flagv(jb,igrd)
         IF( flagv == 0. ) THEN
            pva2d(ii,ij) = dta%v2d(jb)
         ELSE      ! T pts j-indice       on the rim          on the ocean next to the rim on T and V points
            IF( flagv ==  1. ) THEN   ;   ijTrim = ij     ;   ijToce = ij+1   ;   ijVoce = ij+1   ;   ENDIF
            IF( flagv == -1. ) THEN   ;   ijTrim = ij+1   ;   ijToce = ij     ;   ijVoce = ij-1   ;   ENDIF
            !
            ! Rare case : rim is parallel to the mpi subdomain border and on the halo : point will be received
            IF( ijTrim > jpj .OR. ijToce > jpj .OR. ijVoce > jpj .OR. ijVoce < 1 )   CYCLE
            !
            zfla = dta%v2d(jb) - flagv * SQRT( grav * phvr(ii, ij) ) * ( pssh(ii,ijToce) - sshdta(ii,ijTrim) )
            !
            ! jchanut tschanges, use characteristics method (Blayo et Debreu, 2005) :
            ! mix Flather scheme with velocity of the ocean next to the rim
            pva2d(ii,ij) =  z1_2 * ( pva2d(ii,ijVoce) + zfla )
         END IF
      END DO
      !
   END SUBROUTINE bdy_dyn2d_fla


   SUBROUTINE bdy_dyn2d_orlanski( idx, dta, ib_bdy, pua2d, pva2d, pub2d, pvb2d, llrim0, ll_npo )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_dyn2d_orlanski  ***
      !!             
      !!              - Apply Orlanski radiation condition adaptively:
      !!                  - radiation plus weak nudging at outflow points
      !!                  - no radiation and strong nudging at inflow points
      !! 
      !!
      !! References:  Marchesiello, McWilliams and Shchepetkin, Ocean Modelling vol. 3 (2001)    
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),              INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),               INTENT(in) ::   dta  ! OBC external data
      INTEGER,                      INTENT(in) ::   ib_bdy  ! number of current open boundary set
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pua2d, pva2d
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pub2d, pvb2d 
      LOGICAL,                      INTENT(in) ::   ll_npo  ! flag for NPO version
      LOGICAL,                      INTENT(in) ::   llrim0   ! indicate if rim 0 is treated
      INTEGER  ::   ib, igrd                               ! dummy loop indices
      INTEGER  ::   ii, ij, iibm1, ijbm1                   ! indices
      !!----------------------------------------------------------------------
      !
      igrd = 2      ! Orlanski bc on u-velocity; 
      !            
      CALL bdy_orlanski_2d( idx, igrd, pub2d, pua2d, dta%u2d, llrim0, ll_npo )

      igrd = 3      ! Orlanski bc on v-velocity
      !  
      CALL bdy_orlanski_2d( idx, igrd, pvb2d, pva2d, dta%v2d, llrim0, ll_npo )
      !
   END SUBROUTINE bdy_dyn2d_orlanski


   SUBROUTINE bdy_ssh( zssh )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_ssh  ***
      !!
      !! ** Purpose : Duplicate sea level across open boundaries
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,1), INTENT(inout) ::   zssh ! Sea level, need 3 dimensions to be used by bdy_nmn
      !!
      INTEGER ::   ib_bdy, ir      ! bdy index, rim index
      INTEGER ::   ibeg, iend      ! length of rim to be treated (rim 0 or rim 1)
      LOGICAL ::   llrim0          ! indicate if rim 0 is treated
      LOGICAL, DIMENSION(4) :: llsend1, llrecv1  ! indicate how communications are to be carried out
      !!----------------------------------------------------------------------
      llsend1(:) = .false.   ;   llrecv1(:) = .false.
      DO ir = 1, 0, -1   ! treat rim 1 before rim 0
         IF( nn_hls == 1 ) THEN   ;   llsend1(:) = .false.   ;   llrecv1(:) = .false.   ;   END IF
         IF( ir == 0 ) THEN   ;   llrim0 = .TRUE.
         ELSE                 ;   llrim0 = .FALSE.
         END IF
         DO ib_bdy = 1, nb_bdy
            CALL bdy_nmn( idx_bdy(ib_bdy), 1, zssh, llrim0 )   ! zssh is masked
            llsend1(:) = llsend1(:) .OR. lsend_bdyint(ib_bdy,1,:,ir)   ! possibly every direction, T points
            llrecv1(:) = llrecv1(:) .OR. lrecv_bdyint(ib_bdy,1,:,ir)   ! possibly every direction, T points
         END DO
         IF( nn_hls > 1 .AND. ir == 1 ) CYCLE   ! at least 2 halos will be corrected -> no need to correct rim 1 before rim 0
         IF( ANY(llsend1) .OR. ANY(llrecv1) ) THEN   ! if need to send/recv in at least one direction
            CALL lbc_lnk( 'bdydyn2d', zssh(:,:,1), 'T',  1., kfillmode=jpfillnothing ,lsend=llsend1, lrecv=llrecv1 )
         END IF
      END DO
      !
   END SUBROUTINE bdy_ssh

   !!======================================================================
END MODULE bdydyn2d

