MODULE bdytra
   !!======================================================================
   !!                       ***  MODULE  bdytra  ***
   !! Ocean tracers:   Apply boundary conditions for tracers
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!            3.5  !  2012     (S. Mocavero, I. Epicoco) Optimization of BDY communications
   !!            4.0  !  2016     (T. Lovato) Generalize OBC structure
   !!----------------------------------------------------------------------
   !!   bdy_tra       : Apply open boundary conditions & damping to T and S
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables 
   USE bdy_oce        ! ocean open boundary conditions
   USE bdylib         ! for orlanski library routines
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp, ONLY: ctl_stop
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   ! Local structure to rearrange tracers data
   TYPE, PUBLIC ::   ztrabdy
      REAL(wp), POINTER, DIMENSION(:,:) ::  tra
   END TYPE

   PUBLIC   bdy_tra      ! called in tranxt.F90 
   PUBLIC   bdy_tra_dmp  ! called in step.F90 

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_tra( kt )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_tra  ***
      !!
      !! ** Purpose : - Apply open boundary conditions for temperature and salinity
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! Main time step counter
      !
      INTEGER                        :: ib_bdy, jn, igrd, ir   ! Loop indeces
      TYPE(ztrabdy), DIMENSION(jpts) :: zdta                   ! Temporary data structure
      LOGICAL                        :: llrim0                 ! indicate if rim 0 is treated
      LOGICAL, DIMENSION(4)          :: llsend1, llrecv1       ! indicate how communications are to be carried out
      !!----------------------------------------------------------------------
      igrd = 1 
      llsend1(:) = .false.  ;   llrecv1(:) = .false.
      DO ir = 1, 0, -1   ! treat rim 1 before rim 0
         IF( ir == 0 ) THEN   ;   llrim0 = .TRUE.
         ELSE                 ;   llrim0 = .FALSE.
         END IF
         DO ib_bdy=1, nb_bdy
            !
            zdta(1)%tra => dta_bdy(ib_bdy)%tem
            zdta(2)%tra => dta_bdy(ib_bdy)%sal
            !
            DO jn = 1, jpts
               !
               SELECT CASE( TRIM(cn_tra(ib_bdy)) )
               CASE('none'        )   ;   CYCLE
               CASE('frs'         )   ! treat the whole boundary at once
                  IF( ir == 0 ) CALL bdy_frs ( idx_bdy(ib_bdy),                tsa(:,:,:,jn), zdta(jn)%tra )
               CASE('specified'   )   ! treat the whole rim      at once
                  IF( ir == 0 ) CALL bdy_spe ( idx_bdy(ib_bdy),                tsa(:,:,:,jn), zdta(jn)%tra )
               CASE('neumann'     )   ;   CALL bdy_nmn ( idx_bdy(ib_bdy), igrd         , tsa(:,:,:,jn), llrim0 )   ! tsa masked
               CASE('orlanski'    )   ;   CALL bdy_orl ( idx_bdy(ib_bdy), tsb(:,:,:,jn), tsa(:,:,:,jn), &
                    & zdta(jn)%tra, llrim0, ll_npo=.false. )
               CASE('orlanski_npo')   ;   CALL bdy_orl ( idx_bdy(ib_bdy), tsb(:,:,:,jn), tsa(:,:,:,jn), &
                    & zdta(jn)%tra, llrim0, ll_npo=.true.  )
               CASE('runoff'      )   ;   CALL bdy_rnf ( idx_bdy(ib_bdy),                tsa(:,:,:,jn), jn, llrim0 )
               CASE DEFAULT           ;   CALL ctl_stop( 'bdy_tra : unrecognised option for open boundaries for T and S' )
               END SELECT
               ! 
            END DO
         END DO
         !
         IF( nn_hls > 1 .AND. ir == 1 ) CYCLE   ! at least 2 halos will be corrected -> no need to correct rim 1 before rim 0
         IF( nn_hls == 1 ) THEN   ;   llsend1(:) = .false.   ;   llrecv1(:) = .false.   ;   END IF
         DO ib_bdy=1, nb_bdy
            SELECT CASE( TRIM(cn_tra(ib_bdy)) )
            CASE('neumann','runoff')
               llsend1(:) = llsend1(:) .OR. lsend_bdyint(ib_bdy,1,:,ir)   ! possibly every direction, T points
               llrecv1(:) = llrecv1(:) .OR. lrecv_bdyint(ib_bdy,1,:,ir)   ! possibly every direction, T points
            CASE('orlanski', 'orlanski_npo')
               llsend1(:) = llsend1(:) .OR. lsend_bdy(ib_bdy,1,:,ir)   ! possibly every direction, T points
               llrecv1(:) = llrecv1(:) .OR. lrecv_bdy(ib_bdy,1,:,ir)   ! possibly every direction, T points
            END SELECT
         END DO
         IF( ANY(llsend1) .OR. ANY(llrecv1) ) THEN   ! if need to send/recv in at least one direction
            CALL lbc_lnk( 'bdytra', tsa, 'T',  1., kfillmode=jpfillnothing ,lsend=llsend1, lrecv=llrecv1 )
         END IF
         !
      END DO   ! ir
      !
   END SUBROUTINE bdy_tra


   SUBROUTINE bdy_rnf( idx, pta, jpa, llrim0 )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_rnf  ***
      !!                    
      !! ** Purpose : Specialized routine to apply TRA runoff values at OBs:
      !!                  - duplicate the neighbour value for the temperature
      !!                  - specified to 0.1 PSU for the salinity
      !! 
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),                     INTENT(in) ::   idx      ! OBC indices
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pta      ! tracer trend
      INTEGER,                             INTENT(in) ::   jpa      ! TRA index
      LOGICAL,                             INTENT(in) ::   llrim0   ! indicate if rim 0 is treated
      !
      INTEGER  ::   ib, ii, ij, igrd   ! dummy loop indices
      INTEGER  ::   ik, ip, jp ! 2D addresses
      !!----------------------------------------------------------------------
      !
      igrd = 1                       ! Everything is at T-points here
      IF(      jpa == jp_tem ) THEN
         CALL bdy_nmn( idx, igrd, pta, llrim0 )
      ELSE IF( jpa == jp_sal ) THEN
         IF( .NOT. llrim0 )   RETURN
         DO ib = 1, idx%nblenrim(igrd)   ! if llrim0 then treat the whole rim
            ii = idx%nbi(ib,igrd)
            ij = idx%nbj(ib,igrd)
            pta(ii,ij,1:jpkm1) = 0.1 * tmask(ii,ij,1:jpkm1)
         END DO
      END IF
      !
   END SUBROUTINE bdy_rnf


   SUBROUTINE bdy_tra_dmp( kt )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_tra_dmp  ***
      !!                    
      !! ** Purpose : Apply damping for tracers at open boundaries.
      !! 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   !
      !
      REAL(wp) ::   zwgt           ! boundary weight
      REAL(wp) ::   zta, zsa, ztime
      INTEGER  ::   ib, ik, igrd   ! dummy loop indices
      INTEGER  ::   ii, ij         ! 2D addresses
      INTEGER  ::   ib_bdy         ! Loop index
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('bdy_tra_dmp')
      !
      DO ib_bdy = 1, nb_bdy
         IF( ln_tra_dmp(ib_bdy) ) THEN
            igrd = 1                       ! Everything is at T-points here
            DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               zwgt = idx_bdy(ib_bdy)%nbd(ib,igrd)
               DO ik = 1, jpkm1
                  zta = zwgt * ( dta_bdy(ib_bdy)%tem(ib,ik) - tsb(ii,ij,ik,jp_tem) ) * tmask(ii,ij,ik)
                  zsa = zwgt * ( dta_bdy(ib_bdy)%sal(ib,ik) - tsb(ii,ij,ik,jp_sal) ) * tmask(ii,ij,ik)
                  tsa(ii,ij,ik,jp_tem) = tsa(ii,ij,ik,jp_tem) + zta
                  tsa(ii,ij,ik,jp_sal) = tsa(ii,ij,ik,jp_sal) + zsa
               END DO
            END DO
         ENDIF
      END DO
      !
      IF( ln_timing )   CALL timing_stop('bdy_tra_dmp')
      !
   END SUBROUTINE bdy_tra_dmp
 
   !!======================================================================
END MODULE bdytra
