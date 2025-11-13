MODULE icblbc
   !!======================================================================
   !!                       ***  MODULE  icblbc  ***
   !! Ocean physics:  routines to handle boundary exchanges for icebergs
   !!======================================================================
   !! History :  3.3  !  2010-01  (Martin&Adcroft) Original code
   !!             -   !  2011-03  (Madec)          Part conversion to NEMO form
   !!             -   !                            Removal of mapping from another grid
   !!             -   !  2011-04  (Alderson)       Split into separate modules
   !!             -   !  2011-05  (Alderson)       MPP exchanges written based on lib_mpp
   !!             -   !  2011-05  (Alderson)       MPP and single processor boundary conditions added
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   icb_lbc       : -  Pass icebergs across cyclic boundaries
   !!   icb_lbc_mpp   : -  In MPP pass icebergs from linked list between processors
   !!                      as they advect around
   !!                   -  Lagrangian processes cannot be handled by existing NEMO MPP
   !!                      routines because they do not lie on regular jpi,jpj grids
   !!                   -  Processor exchanges are handled as in lib_mpp whenever icebergs step 
   !!                      across boundary of interior domain (nicbdi-nicbei, nicbdj-nicbej)
   !!                      so that iceberg does not exist in more than one processor
   !!                   -  North fold exchanges controlled by three arrays:
   !!                         nicbflddest - unique processor numbers that current one exchanges with
   !!                         nicbfldproc - processor number that current grid point exchanges with
   !!                         nicbfldpts  - packed i,j point in exchanging processor
   !!----------------------------------------------------------------------
   USE par_oce                             ! ocean parameters
   USE dom_oce                             ! ocean domain
   USE in_out_manager                      ! IO parameters
   USE lib_mpp                             ! MPI code and lk_mpp in particular
   USE icb_oce                             ! define iceberg arrays
   USE icbutl                              ! iceberg utility routines

   IMPLICIT NONE
   PRIVATE

#if defined key_mpp_mpi

!$AGRIF_DO_NOT_TREAT
   INCLUDE 'mpif.h'
!$AGRIF_END_DO_NOT_TREAT

   TYPE, PUBLIC :: buffer
      INTEGER :: size = 0
      REAL(wp), DIMENSION(:,:), POINTER ::   data
   END TYPE buffer

   TYPE(buffer), POINTER       ::   obuffer_n=>NULL() , ibuffer_n=>NULL()
   TYPE(buffer), POINTER       ::   obuffer_s=>NULL() , ibuffer_s=>NULL()
   TYPE(buffer), POINTER       ::   obuffer_e=>NULL() , ibuffer_e=>NULL()
   TYPE(buffer), POINTER       ::   obuffer_w=>NULL() , ibuffer_w=>NULL()

   ! north fold exchange buffers
   TYPE(buffer), POINTER       ::   obuffer_f=>NULL() , ibuffer_f=>NULL()

   INTEGER, PARAMETER, PRIVATE ::   jp_delta_buf = 25             ! Size by which to increment buffers
   INTEGER, PARAMETER, PRIVATE ::   jp_buffer_width = 15+nkounts  ! items to store for each berg

#endif

   PUBLIC   icb_lbc
   PUBLIC   icb_lbc_mpp

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE icb_lbc()
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE icb_lbc  ***
      !!
      !! ** Purpose :   in non-mpp case need to deal with cyclic conditions
      !!                including north-fold
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER ::   this
      TYPE(point)  , POINTER ::   pt
      INTEGER                ::   iine
      !!----------------------------------------------------------------------

      !! periodic east/west boundaries
      !! =============================

      IF( l_Iperio ) THEN

         this => first_berg
         DO WHILE( ASSOCIATED(this) )
            pt => this%current_point
            iine = INT( pt%xi + 0.5 )
            IF( iine > mig(nicbei) ) THEN
               pt%xi = ricb_right + MOD(pt%xi, 1._wp ) - 1._wp
            ELSE IF( iine < mig(nicbdi) ) THEN
               pt%xi = ricb_left + MOD(pt%xi, 1._wp )
            ENDIF
            this => this%next
         END DO
         !
      ENDIF

      !! north/south boundaries
      !! ======================
      IF( l_Jperio)      CALL ctl_stop(' north-south periodicity not implemented for icebergs')
      ! north fold
      IF( npolj /= 0 )   CALL icb_lbc_nfld()
      !
   END SUBROUTINE icb_lbc


   SUBROUTINE icb_lbc_nfld()
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE icb_lbc_nfld  ***
      !!
      !! ** Purpose :   single processor north fold exchange
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER ::   this
      TYPE(point)  , POINTER ::   pt
      INTEGER                ::   iine, ijne, ipts
      INTEGER                ::   iiglo, ijglo
      !!----------------------------------------------------------------------
      !
      this => first_berg
      DO WHILE( ASSOCIATED(this) )
         pt => this%current_point
         ijne = INT( pt%yj + 0.5 )
         IF( ijne .GT. mjg(nicbej) ) THEN
            !
            iine = INT( pt%xi + 0.5 )
            ipts  = nicbfldpts (mi1(iine))
            !
            ! moving across the cut line means both position and
            ! velocity must change
            ijglo = INT( ipts/nicbpack )
            iiglo = ipts - nicbpack*ijglo
            pt%xi = iiglo - ( pt%xi - REAL(iine,wp) )
            pt%yj = ijglo - ( pt%yj - REAL(ijne,wp) )
            pt%uvel = -1._wp * pt%uvel
            pt%vvel = -1._wp * pt%vvel
         ENDIF
         this => this%next
      END DO
      !
   END SUBROUTINE icb_lbc_nfld

#if defined key_mpp_mpi
   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'             MPI massively parallel processing library
   !!----------------------------------------------------------------------

   SUBROUTINE icb_lbc_mpp()
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE icb_lbc_mpp  ***
      !!
      !! ** Purpose :   multi processor exchange
      !!
      !! ** Method  :   identify direction for exchange, pack into a buffer
      !!                which is basically a real array and delete from linked list
      !!                length of buffer is exchanged first with receiving processor
      !!                then buffer is sent if necessary
      !!----------------------------------------------------------------------
      TYPE(iceberg)         , POINTER     ::   tmpberg, this
      TYPE(point)           , POINTER     ::   pt
      INTEGER                             ::   ibergs_to_send_e, ibergs_to_send_w
      INTEGER                             ::   ibergs_to_send_n, ibergs_to_send_s
      INTEGER                             ::   ibergs_rcvd_from_e, ibergs_rcvd_from_w
      INTEGER                             ::   ibergs_rcvd_from_n, ibergs_rcvd_from_s
      INTEGER                             ::   i, ibergs_start, ibergs_end
      INTEGER                             ::   iine, ijne
      INTEGER                             ::   ipe_N, ipe_S, ipe_W, ipe_E
      REAL(wp), DIMENSION(2)              ::   zewbergs, zwebergs, znsbergs, zsnbergs
      INTEGER                             ::   iml_req1, iml_req2, iml_req3, iml_req4
      INTEGER                             ::   iml_req5, iml_req6, iml_req7, iml_req8, iml_err
      INTEGER, DIMENSION(MPI_STATUS_SIZE) ::   iml_stat

      ! set up indices of neighbouring processors
      ipe_N = -1
      ipe_S = -1
      ipe_W = -1
      ipe_E = -1
      IF( nbondi .EQ.  0 .OR. nbondi .EQ. 1) ipe_W = nowe
      IF( nbondi .EQ. -1 .OR. nbondi .EQ. 0) ipe_E = noea
      IF( nbondj .EQ.  0 .OR. nbondj .EQ. 1) ipe_S = noso
      IF( nbondj .EQ. -1 .OR. nbondj .EQ. 0) ipe_N = nono
      !
      ! at northern line of processors with north fold handle bergs differently
      IF( npolj > 0 ) ipe_N = -1

      ! if there's only one processor in x direction then don't let mpp try to handle periodicity
      IF( jpni == 1 ) THEN
         ipe_E = -1
         ipe_W = -1
      ENDIF

      IF( nn_verbose_level >= 2 ) THEN
         WRITE(numicb,*) 'processor west  : ', ipe_W
         WRITE(numicb,*) 'processor east  : ', ipe_E
         WRITE(numicb,*) 'processor north : ', ipe_N
         WRITE(numicb,*) 'processor south : ', ipe_S
         WRITE(numicb,*) 'processor nimpp : ', nimpp
         WRITE(numicb,*) 'processor njmpp : ', njmpp
         WRITE(numicb,*) 'processor nbondi: ', nbondi
         WRITE(numicb,*) 'processor nbondj: ', nbondj
         CALL flush( numicb )
      ENDIF

      ! periodicity is handled here when using mpp when there is more than one processor in
      ! the i direction, but it also has to happen when jpni=1 case so this is dealt with
      ! in icb_lbc and called here

      IF( jpni == 1 ) CALL icb_lbc()

      ! Note that xi is adjusted when swapping because of periodic condition

      IF( nn_verbose_level > 0 ) THEN
         ! store the number of icebergs on this processor at start
         ibergs_start = icb_utl_count()
      ENDIF

      ibergs_to_send_e   = 0
      ibergs_to_send_w   = 0
      ibergs_to_send_n   = 0
      ibergs_to_send_s   = 0
      ibergs_rcvd_from_e = 0
      ibergs_rcvd_from_w = 0
      ibergs_rcvd_from_n = 0
      ibergs_rcvd_from_s = 0

      IF( ASSOCIATED(first_berg) ) THEN      ! Find number of bergs that headed east/west
         this => first_berg
         DO WHILE (ASSOCIATED(this))
            pt => this%current_point
            iine = INT( pt%xi + 0.5 )
            IF( ipe_E >= 0 .AND. iine > mig(nicbei) ) THEN
               tmpberg => this
               this => this%next
               ibergs_to_send_e = ibergs_to_send_e + 1
               IF( nn_verbose_level >= 4 ) THEN
                  WRITE(numicb,*) 'bergstep ',nktberg,' packing berg ',tmpberg%number(:),' for transfer to east'
                  CALL flush( numicb )
               ENDIF
               ! deal with periodic case
               tmpberg%current_point%xi = ricb_right + MOD(tmpberg%current_point%xi, 1._wp ) - 1._wp
               ! now pack it into buffer and delete from list
               CALL icb_pack_into_buffer( tmpberg, obuffer_e, ibergs_to_send_e)
               CALL icb_utl_delete(first_berg, tmpberg)
            ELSE IF( ipe_W >= 0 .AND. iine < mig(nicbdi) ) THEN
               tmpberg => this
               this => this%next
               ibergs_to_send_w = ibergs_to_send_w + 1
               IF( nn_verbose_level >= 4 ) THEN
                  WRITE(numicb,*) 'bergstep ',nktberg,' packing berg ',tmpberg%number(:),' for transfer to west'
                  CALL flush( numicb )
               ENDIF
               ! deal with periodic case
               tmpberg%current_point%xi = ricb_left + MOD(tmpberg%current_point%xi, 1._wp )
               ! now pack it into buffer and delete from list
               CALL icb_pack_into_buffer( tmpberg, obuffer_w, ibergs_to_send_w)
               CALL icb_utl_delete(first_berg, tmpberg)
            ELSE
               this => this%next
            ENDIF
         END DO
      ENDIF
      IF( nn_verbose_level >= 3) THEN
         WRITE(numicb,*) 'bergstep ',nktberg,' send ew: ', ibergs_to_send_e, ibergs_to_send_w
         CALL flush(numicb)
      ENDIF

      ! send bergs east and receive bergs from west (ie ones that were sent east) and vice versa

      ! pattern here is copied from lib_mpp code

      SELECT CASE ( nbondi )
      CASE( -1 )
         zwebergs(1) = ibergs_to_send_e
         CALL mppsend( 12, zwebergs(1), 1, ipe_E, iml_req1)
         CALL mpprecv( 11, zewbergs(2), 1, ipe_E )
         CALL mpi_wait( iml_req1, iml_stat, iml_err )
         ibergs_rcvd_from_e = INT( zewbergs(2) )
      CASE(  0 )
         zewbergs(1) = ibergs_to_send_w
         zwebergs(1) = ibergs_to_send_e
         CALL mppsend( 11, zewbergs(1), 1, ipe_W, iml_req2)
         CALL mppsend( 12, zwebergs(1), 1, ipe_E, iml_req3)
         CALL mpprecv( 11, zewbergs(2), 1, ipe_E )
         CALL mpprecv( 12, zwebergs(2), 1, ipe_W )
         CALL mpi_wait( iml_req2, iml_stat, iml_err )
         CALL mpi_wait( iml_req3, iml_stat, iml_err )
         ibergs_rcvd_from_e = INT( zewbergs(2) )
         ibergs_rcvd_from_w = INT( zwebergs(2) )
      CASE(  1 )
         zewbergs(1) = ibergs_to_send_w
         CALL mppsend( 11, zewbergs(1), 1, ipe_W, iml_req4)
         CALL mpprecv( 12, zwebergs(2), 1, ipe_W )
         CALL mpi_wait( iml_req4, iml_stat, iml_err )
         ibergs_rcvd_from_w = INT( zwebergs(2) )
      END SELECT
      IF( nn_verbose_level >= 3) THEN
         WRITE(numicb,*) 'bergstep ',nktberg,' recv ew: ', ibergs_rcvd_from_w, ibergs_rcvd_from_e
         CALL flush(numicb)
      ENDIF

      SELECT CASE ( nbondi )
      CASE( -1 )
         IF( ibergs_to_send_e > 0 ) CALL mppsend( 14, obuffer_e%data, ibergs_to_send_e*jp_buffer_width, ipe_E, iml_req1 )
         IF( ibergs_rcvd_from_e > 0 ) THEN
            CALL icb_increase_ibuffer(ibuffer_e, ibergs_rcvd_from_e)
            CALL mpprecv( 13, ibuffer_e%data, ibergs_rcvd_from_e*jp_buffer_width )
         ENDIF
         IF( ibergs_to_send_e > 0 ) CALL mpi_wait( iml_req1, iml_stat, iml_err )
         DO i = 1, ibergs_rcvd_from_e
            IF( nn_verbose_level >= 4 ) THEN
               WRITE(numicb,*) 'bergstep ',nktberg,' unpacking berg ',INT(ibuffer_e%data(16,i)),' from east'
               CALL flush( numicb )
            ENDIF
            CALL icb_unpack_from_buffer(first_berg, ibuffer_e, i)
         ENDDO
      CASE(  0 )
         IF( ibergs_to_send_w > 0 ) CALL mppsend( 13, obuffer_w%data, ibergs_to_send_w*jp_buffer_width, ipe_W, iml_req2 )
         IF( ibergs_to_send_e > 0 ) CALL mppsend( 14, obuffer_e%data, ibergs_to_send_e*jp_buffer_width, ipe_E, iml_req3 )
         IF( ibergs_rcvd_from_e > 0 ) THEN
            CALL icb_increase_ibuffer(ibuffer_e, ibergs_rcvd_from_e)
            CALL mpprecv( 13, ibuffer_e%data, ibergs_rcvd_from_e*jp_buffer_width )
         ENDIF
         IF( ibergs_rcvd_from_w > 0 ) THEN
            CALL icb_increase_ibuffer(ibuffer_w, ibergs_rcvd_from_w)
            CALL mpprecv( 14, ibuffer_w%data, ibergs_rcvd_from_w*jp_buffer_width )
         ENDIF
         IF( ibergs_to_send_w > 0 ) CALL mpi_wait( iml_req2, iml_stat, iml_err )
         IF( ibergs_to_send_e > 0 ) CALL mpi_wait( iml_req3, iml_stat, iml_err )
         DO i = 1, ibergs_rcvd_from_e
            IF( nn_verbose_level >= 4 ) THEN
               WRITE(numicb,*) 'bergstep ',nktberg,' unpacking berg ',INT(ibuffer_e%data(16,i)),' from east'
               CALL flush( numicb )
            ENDIF
            CALL icb_unpack_from_buffer(first_berg, ibuffer_e, i)
         END DO
         DO i = 1, ibergs_rcvd_from_w
            IF( nn_verbose_level >= 4 ) THEN
               WRITE(numicb,*) 'bergstep ',nktberg,' unpacking berg ',INT(ibuffer_w%data(16,i)),' from west'
               CALL flush( numicb )
            ENDIF
            CALL icb_unpack_from_buffer(first_berg, ibuffer_w, i)
         ENDDO
      CASE(  1 )
         IF( ibergs_to_send_w > 0 ) CALL mppsend( 13, obuffer_w%data, ibergs_to_send_w*jp_buffer_width, ipe_W, iml_req4 )
         IF( ibergs_rcvd_from_w > 0 ) THEN
            CALL icb_increase_ibuffer(ibuffer_w, ibergs_rcvd_from_w)
            CALL mpprecv( 14, ibuffer_w%data, ibergs_rcvd_from_w*jp_buffer_width )
         ENDIF
         IF( ibergs_to_send_w > 0 ) CALL mpi_wait( iml_req4, iml_stat, iml_err )
         DO i = 1, ibergs_rcvd_from_w
            IF( nn_verbose_level >= 4 ) THEN
               WRITE(numicb,*) 'bergstep ',nktberg,' unpacking berg ',INT(ibuffer_w%data(16,i)),' from west'
               CALL flush( numicb )
            ENDIF
            CALL icb_unpack_from_buffer(first_berg, ibuffer_w, i)
         END DO
      END SELECT

      ! Find number of bergs that headed north/south
      ! (note: this block should technically go ahead of the E/W recv block above
      !  to handle arbitrary orientation of PEs. But for simplicity, it is
      !  here to accomodate diagonal transfer of bergs between PEs -AJA)

      IF( ASSOCIATED(first_berg) ) THEN
         this => first_berg
         DO WHILE (ASSOCIATED(this))
            pt => this%current_point
            ijne = INT( pt%yj + 0.5 )
            IF( ipe_N >= 0 .AND. ijne .GT. mjg(nicbej) ) THEN
               tmpberg => this
               this => this%next
               ibergs_to_send_n = ibergs_to_send_n + 1
               IF( nn_verbose_level >= 4 ) THEN
                  WRITE(numicb,*) 'bergstep ',nktberg,' packing berg ',tmpberg%number(:),' for transfer to north'
                  CALL flush( numicb )
               ENDIF
               CALL icb_pack_into_buffer( tmpberg, obuffer_n, ibergs_to_send_n)
               CALL icb_utl_delete(first_berg, tmpberg)
            ELSE IF( ipe_S >= 0 .AND. ijne .LT. mjg(nicbdj) ) THEN
               tmpberg => this
               this => this%next
               ibergs_to_send_s = ibergs_to_send_s + 1
               IF( nn_verbose_level >= 4 ) THEN
                  WRITE(numicb,*) 'bergstep ',nktberg,' packing berg ',tmpberg%number(:),' for transfer to south'
                  CALL flush( numicb )
               ENDIF
               CALL icb_pack_into_buffer( tmpberg, obuffer_s, ibergs_to_send_s)
               CALL icb_utl_delete(first_berg, tmpberg)
            ELSE
               this => this%next
            ENDIF
         END DO
      ENDIF
      if( nn_verbose_level >= 3) then
         write(numicb,*) 'bergstep ',nktberg,' send ns: ', ibergs_to_send_n, ibergs_to_send_s
         call flush(numicb)
      endif

      ! send bergs north
      ! and receive bergs from south (ie ones sent north)

      SELECT CASE ( nbondj )
      CASE( -1 )
         zsnbergs(1) = ibergs_to_send_n
         CALL mppsend( 16, zsnbergs(1), 1, ipe_N, iml_req1)
         CALL mpprecv( 15, znsbergs(2), 1, ipe_N )
         CALL mpi_wait( iml_req1, iml_stat, iml_err )
         ibergs_rcvd_from_n = INT( znsbergs(2) )
      CASE(  0 )
         znsbergs(1) = ibergs_to_send_s
         zsnbergs(1) = ibergs_to_send_n
         CALL mppsend( 15, znsbergs(1), 1, ipe_S, iml_req2)
         CALL mppsend( 16, zsnbergs(1), 1, ipe_N, iml_req3)
         CALL mpprecv( 15, znsbergs(2), 1, ipe_N )
         CALL mpprecv( 16, zsnbergs(2), 1, ipe_S )
         CALL mpi_wait( iml_req2, iml_stat, iml_err )
         CALL mpi_wait( iml_req3, iml_stat, iml_err )
         ibergs_rcvd_from_n = INT( znsbergs(2) )
         ibergs_rcvd_from_s = INT( zsnbergs(2) )
      CASE(  1 )
         znsbergs(1) = ibergs_to_send_s
         CALL mppsend( 15, znsbergs(1), 1, ipe_S, iml_req4)
         CALL mpprecv( 16, zsnbergs(2), 1, ipe_S )
         CALL mpi_wait( iml_req4, iml_stat, iml_err )
         ibergs_rcvd_from_s = INT( zsnbergs(2) )
      END SELECT
      if( nn_verbose_level >= 3) then
         write(numicb,*) 'bergstep ',nktberg,' recv ns: ', ibergs_rcvd_from_s, ibergs_rcvd_from_n
         call flush(numicb)
      endif

      SELECT CASE ( nbondj )
      CASE( -1 )
         IF( ibergs_to_send_n > 0 ) CALL mppsend( 18, obuffer_n%data, ibergs_to_send_n*jp_buffer_width, ipe_N, iml_req1 )
         IF( ibergs_rcvd_from_n > 0 ) THEN
            CALL icb_increase_ibuffer(ibuffer_n, ibergs_rcvd_from_n)
            CALL mpprecv( 17, ibuffer_n%data, ibergs_rcvd_from_n*jp_buffer_width )
         ENDIF
         IF( ibergs_to_send_n > 0 ) CALL mpi_wait( iml_req1, iml_stat, iml_err )
         DO i = 1, ibergs_rcvd_from_n
            IF( nn_verbose_level >= 4 ) THEN
               WRITE(numicb,*) 'bergstep ',nktberg,' unpacking berg ',INT(ibuffer_n%data(16,i)),' from north'
               CALL flush( numicb )
            ENDIF
            CALL icb_unpack_from_buffer(first_berg, ibuffer_n, i)
         END DO
      CASE(  0 )
         IF( ibergs_to_send_s > 0 ) CALL mppsend( 17, obuffer_s%data, ibergs_to_send_s*jp_buffer_width, ipe_S, iml_req2 )
         IF( ibergs_to_send_n > 0 ) CALL mppsend( 18, obuffer_n%data, ibergs_to_send_n*jp_buffer_width, ipe_N, iml_req3 )
         IF( ibergs_rcvd_from_n > 0 ) THEN
            CALL icb_increase_ibuffer(ibuffer_n, ibergs_rcvd_from_n)
            CALL mpprecv( 17, ibuffer_n%data, ibergs_rcvd_from_n*jp_buffer_width )
         ENDIF
         IF( ibergs_rcvd_from_s > 0 ) THEN
            CALL icb_increase_ibuffer(ibuffer_s, ibergs_rcvd_from_s)
            CALL mpprecv( 18, ibuffer_s%data, ibergs_rcvd_from_s*jp_buffer_width )
         ENDIF
         IF( ibergs_to_send_s > 0 ) CALL mpi_wait( iml_req2, iml_stat, iml_err )
         IF( ibergs_to_send_n > 0 ) CALL mpi_wait( iml_req3, iml_stat, iml_err )
         DO i = 1, ibergs_rcvd_from_n
            IF( nn_verbose_level >= 4 ) THEN
               WRITE(numicb,*) 'bergstep ',nktberg,' unpacking berg ',INT(ibuffer_n%data(16,i)),' from north'
               CALL flush( numicb )
            ENDIF
            CALL icb_unpack_from_buffer(first_berg, ibuffer_n, i)
         END DO
         DO i = 1, ibergs_rcvd_from_s
            IF( nn_verbose_level >= 4 ) THEN
               WRITE(numicb,*) 'bergstep ',nktberg,' unpacking berg ',INT(ibuffer_s%data(16,i)),' from south'
               CALL flush( numicb )
            ENDIF
            CALL icb_unpack_from_buffer(first_berg, ibuffer_s, i)
         ENDDO
      CASE(  1 )
         IF( ibergs_to_send_s > 0 ) CALL mppsend( 17, obuffer_s%data, ibergs_to_send_s*jp_buffer_width, ipe_S, iml_req4 )
         IF( ibergs_rcvd_from_s > 0 ) THEN
            CALL icb_increase_ibuffer(ibuffer_s, ibergs_rcvd_from_s)
            CALL mpprecv( 18, ibuffer_s%data, ibergs_rcvd_from_s*jp_buffer_width )
         ENDIF
         IF( ibergs_to_send_s > 0 ) CALL mpi_wait( iml_req4, iml_stat, iml_err )
         DO i = 1, ibergs_rcvd_from_s
            IF( nn_verbose_level >= 4 ) THEN
               WRITE(numicb,*) 'bergstep ',nktberg,' unpacking berg ',INT(ibuffer_s%data(16,i)),' from south'
               CALL flush( numicb )
            ENDIF
            CALL icb_unpack_from_buffer(first_berg, ibuffer_s, i)
         END DO
      END SELECT

      IF( nn_verbose_level > 0 ) THEN
         ! compare the number of icebergs on this processor from the start to the end
         ibergs_end = icb_utl_count()
         i = ( ibergs_rcvd_from_n + ibergs_rcvd_from_s + ibergs_rcvd_from_e + ibergs_rcvd_from_w ) - &
             ( ibergs_to_send_n + ibergs_to_send_s + ibergs_to_send_e + ibergs_to_send_w )
         IF( ibergs_end-(ibergs_start+i) .NE. 0 ) THEN
            WRITE( numicb,*   ) 'send_bergs_to_other_pes: net change in number of icebergs'
            WRITE( numicb,1000) 'send_bergs_to_other_pes: ibergs_end=', &
                                ibergs_end,' on PE',narea
            WRITE( numicb,1000) 'send_bergs_to_other_pes: ibergs_start=', &
                                ibergs_start,' on PE',narea
            WRITE( numicb,1000) 'send_bergs_to_other_pes: delta=', &
                                i,' on PE',narea
            WRITE( numicb,1000) 'send_bergs_to_other_pes: error=', &
                                ibergs_end-(ibergs_start+i),' on PE',narea
            WRITE( numicb,1000) 'send_bergs_to_other_pes: ibergs_to_send_n=', &
                                ibergs_to_send_n,' on PE',narea
            WRITE( numicb,1000) 'send_bergs_to_other_pes: ibergs_to_send_s=', &
                                ibergs_to_send_s,' on PE',narea
            WRITE( numicb,1000) 'send_bergs_to_other_pes: ibergs_to_send_e=', &
                                ibergs_to_send_e,' on PE',narea
            WRITE( numicb,1000) 'send_bergs_to_other_pes: ibergs_to_send_w=', &
                                ibergs_to_send_w,' on PE',narea
            WRITE( numicb,1000) 'send_bergs_to_other_pes: ibergs_rcvd_from_n=', &
                                ibergs_rcvd_from_n,' on PE',narea
            WRITE( numicb,1000) 'send_bergs_to_other_pes: ibergs_rcvd_from_s=', &
                                ibergs_rcvd_from_s,' on PE',narea
            WRITE( numicb,1000) 'send_bergs_to_other_pes: ibergs_rcvd_from_e=', &
                                ibergs_rcvd_from_e,' on PE',narea
            WRITE( numicb,1000) 'send_bergs_to_other_pes: ibergs_rcvd_from_w=', &
                                ibergs_rcvd_from_w,' on PE',narea
  1000      FORMAT(a,i5,a,i4)
            CALL ctl_stop('send_bergs_to_other_pes: lost or gained an iceberg or two')
         ENDIF
      ENDIF

      ! deal with north fold if we necessary when there is more than one top row processor
      ! note that for jpni=1 north fold has been dealt with above in call to icb_lbc
      IF( npolj /= 0 .AND. jpni > 1 ) CALL icb_lbc_mpp_nfld( )

      IF( nn_verbose_level > 0 ) THEN
         i = 0
         this => first_berg
         DO WHILE (ASSOCIATED(this))
            pt => this%current_point
            iine = INT( pt%xi + 0.5 )
            ijne = INT( pt%yj + 0.5 )
            IF( iine .LT. mig(nicbdi) .OR. &
                iine .GT. mig(nicbei) .OR. &
                ijne .LT. mjg(nicbdj) .OR. &
                ijne .GT. mjg(nicbej)) THEN
               i = i + 1
               WRITE(numicb,*) 'berg lost in halo: ', this%number(:),iine,ijne
               WRITE(numicb,*) '                   ', nimpp, njmpp
               WRITE(numicb,*) '                   ', nicbdi, nicbei, nicbdj, nicbej
               CALL flush( numicb )
            ENDIF
            this => this%next
         ENDDO ! WHILE
         CALL mpp_sum('icblbc', i)
         IF( i .GT. 0 ) THEN
            WRITE( numicb,'(a,i4)') 'send_bergs_to_other_pes: # of bergs outside computational domain = ',i
            CALL ctl_stop('send_bergs_to_other_pes:  there are bergs still in halos!')
         ENDIF ! root_pe
      ENDIF ! debug
      !
      CALL mppsync()
      !
   END SUBROUTINE icb_lbc_mpp


   SUBROUTINE icb_lbc_mpp_nfld()
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE icb_lbc_mpp_nfld  ***
      !!
      !! ** Purpose :   north fold treatment in multi processor exchange
      !!
      !! ** Method  :   
      !!----------------------------------------------------------------------
      TYPE(iceberg)         , POINTER     :: tmpberg, this
      TYPE(point)           , POINTER     :: pt
      INTEGER                             :: ibergs_to_send
      INTEGER                             :: ibergs_to_rcv
      INTEGER                             :: iiglo, ijglo, jk, jn
      INTEGER                             :: ifldproc, iproc, ipts
      INTEGER                             :: iine, ijne
      INTEGER                             :: jjn
      REAL(wp), DIMENSION(0:3)            :: zsbergs, znbergs
      INTEGER                             :: iml_req1, iml_req2, iml_err
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: iml_stat

      ! set up indices of neighbouring processors

      ! nicbfldproc is a list of unique processor numbers that this processor
      ! exchanges with (including itself), so we loop over this array; since
      ! its of fixed size, the first -1 marks end of list of processors
      !
      nicbfldnsend(:) = 0
      nicbfldexpect(:) = 0
      nicbfldreq(:) = 0
      !
      ! Since each processor may be communicating with more than one northern
      ! neighbour, cycle through the sends so that the receive order can be
      ! controlled.
      !
      ! First compute how many icebergs each active neighbour should expect
      DO jn = 1, jpni
         IF( nicbfldproc(jn) /= -1 ) THEN
            ifldproc = nicbfldproc(jn)
            nicbfldnsend(jn) = 0

            ! Find number of bergs that need to be exchanged
            ! Pick out exchanges with processor ifldproc
            ! if ifldproc is this processor then don't send
            !
            IF( ASSOCIATED(first_berg) ) THEN
               this => first_berg
               DO WHILE (ASSOCIATED(this))
                  pt => this%current_point
                  iine = INT( pt%xi + 0.5 )
                  ijne = INT( pt%yj + 0.5 )
                  iproc = nicbflddest(mi1(iine))
                  IF( ijne .GT. mjg(nicbej) ) THEN
                     IF( iproc == ifldproc ) THEN
                        !
                        IF( iproc /= narea ) THEN
                           tmpberg => this
                           nicbfldnsend(jn) = nicbfldnsend(jn) + 1
                        ENDIF
                        !
                     ENDIF
                  ENDIF
                  this => this%next
               END DO
            ENDIF
            !
         ENDIF
         !
      END DO
      !
      ! Now tell each active neighbour how many icebergs to expect
      DO jn = 1, jpni
         IF( nicbfldproc(jn) /= -1 ) THEN
            ifldproc = nicbfldproc(jn)
            IF( ifldproc == narea ) CYCLE
   
            zsbergs(0) = narea
            zsbergs(1) = nicbfldnsend(jn)
            !IF ( nicbfldnsend(jn) .GT. 0 .AND. nn_verbose_level > 0 ) write(numicb,*) 'ICB sending ',nicbfldnsend(jn),' to ', ifldproc
            CALL mppsend( 21, zsbergs(0:1), 2, ifldproc-1, nicbfldreq(jn))
         ENDIF
         !
      END DO
      !
      ! and receive the heads-up from active neighbours preparing to send
      DO jn = 1, jpni
         IF( nicbfldproc(jn) /= -1 ) THEN
            ifldproc = nicbfldproc(jn)
            IF( ifldproc == narea ) CYCLE

            CALL mpprecv( 21, znbergs(1:2), 2 )
            DO jjn = 1,jpni
             IF( nicbfldproc(jjn) .eq. INT(znbergs(1)) ) EXIT
            END DO
            IF( jjn .GT. jpni .AND. nn_verbose_level > 0 ) write(numicb,*) 'ICB ERROR'
            nicbfldexpect(jjn) = INT( znbergs(2) )
            !IF ( nicbfldexpect(jjn) .GT. 0 .AND. nn_verbose_level > 0 ) write(numicb,*) 'ICB expecting ',nicbfldexpect(jjn),' from ', nicbfldproc(jjn)
            !IF (nn_verbose_level > 0) CALL FLUSH(numicb)
         ENDIF
         !
      END DO
      !
      ! post the mpi waits if using immediate send protocol
      DO jn = 1, jpni
         IF( nicbfldproc(jn) /= -1 ) THEN
            ifldproc = nicbfldproc(jn)
            IF( ifldproc == narea ) CYCLE
            CALL mpi_wait( nicbfldreq(jn), iml_stat, iml_err )
         ENDIF
         !
      END DO
   
         !
         ! Cycle through the icebergs again, this time packing and sending any
         ! going through the north fold. They will be expected.
      DO jn = 1, jpni
         IF( nicbfldproc(jn) /= -1 ) THEN
            ifldproc = nicbfldproc(jn)
            ibergs_to_send = 0
   
            ! Find number of bergs that need to be exchanged
            ! Pick out exchanges with processor ifldproc
            ! if ifldproc is this processor then don't send
            !
            IF( ASSOCIATED(first_berg) ) THEN
               this => first_berg
               DO WHILE (ASSOCIATED(this))
                  pt => this%current_point
                  iine = INT( pt%xi + 0.5 )
                  ijne = INT( pt%yj + 0.5 )
                  ipts  = nicbfldpts (mi1(iine))
                  iproc = nicbflddest(mi1(iine))
                  IF( ijne .GT. mjg(nicbej) ) THEN
                     IF( iproc == ifldproc ) THEN
                        !
                        ! moving across the cut line means both position and
                        ! velocity must change
                        ijglo = INT( ipts/nicbpack )
                        iiglo = ipts - nicbpack*ijglo
                        pt%xi = iiglo - ( pt%xi - REAL(iine,wp) )
                        pt%yj = ijglo - ( pt%yj - REAL(ijne,wp) )
                        pt%uvel = -1._wp * pt%uvel
                        pt%vvel = -1._wp * pt%vvel
                        !
                        ! now remove berg from list and pack it into a buffer
                        IF( iproc /= narea ) THEN
                           tmpberg => this
                           ibergs_to_send = ibergs_to_send + 1
                           IF( nn_verbose_level >= 4 ) THEN
                              WRITE(numicb,*) 'bergstep ',nktberg,' packing berg ',tmpberg%number(:),' for north fold'
                              CALL flush( numicb )
                           ENDIF
                           CALL icb_pack_into_buffer( tmpberg, obuffer_f, ibergs_to_send)
                           CALL icb_utl_delete(first_berg, tmpberg)
                        ENDIF
                        !
                     ENDIF
                  ENDIF
                  this => this%next
               END DO
            ENDIF
            if( nn_verbose_level >= 3) then
               write(numicb,*) 'bergstep ',nktberg,' send nfld: ', ibergs_to_send
               call flush(numicb)
            endif
            !
            ! if we're in this processor, then we've done everything we need to
            ! so go on to next element of loop
            IF( ifldproc == narea ) CYCLE
   
            ! send bergs
   
            IF( ibergs_to_send > 0 )  &
                CALL mppsend( 12, obuffer_f%data, ibergs_to_send*jp_buffer_width, ifldproc-1, nicbfldreq(jn) )
            !
         ENDIF
         !
      END DO
      !
      ! Now receive the expected number of bergs from the active neighbours
      DO jn = 1, jpni
         IF( nicbfldproc(jn) /= -1 ) THEN
            ifldproc = nicbfldproc(jn)
            IF( ifldproc == narea ) CYCLE
            ibergs_to_rcv = nicbfldexpect(jn)

            IF( ibergs_to_rcv  > 0 ) THEN
               CALL icb_increase_ibuffer(ibuffer_f, ibergs_to_rcv)
               CALL mpprecv( 12, ibuffer_f%data, ibergs_to_rcv*jp_buffer_width, ifldproc-1 )
            ENDIF
            !
            DO jk = 1, ibergs_to_rcv
               IF( nn_verbose_level >= 4 ) THEN
                  WRITE(numicb,*) 'bergstep ',nktberg,' unpacking berg ',INT(ibuffer_f%data(16,jk)),' from north fold'
                  CALL flush( numicb )
               ENDIF
               CALL icb_unpack_from_buffer(first_berg, ibuffer_f, jk )
            END DO
         ENDIF
         !
      END DO
      !
      ! Finally post the mpi waits if using immediate send protocol
      DO jn = 1, jpni
         IF( nicbfldproc(jn) /= -1 ) THEN
            ifldproc = nicbfldproc(jn)
            IF( ifldproc == narea ) CYCLE
            CALL mpi_wait( nicbfldreq(jn), iml_stat, iml_err )
         ENDIF
         !
      END DO
      !
   END SUBROUTINE icb_lbc_mpp_nfld


   SUBROUTINE icb_pack_into_buffer( berg, pbuff, kb )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER :: berg
      TYPE(buffer) , POINTER :: pbuff
      INTEGER               , INTENT(in) :: kb
      ! 
      INTEGER ::   k   ! local integer
      !!----------------------------------------------------------------------
      !
      IF( .NOT. ASSOCIATED(pbuff) ) CALL icb_increase_buffer( pbuff, jp_delta_buf )
      IF( kb .GT. pbuff%size ) CALL icb_increase_buffer( pbuff, jp_delta_buf )

      !! pack points into buffer

      pbuff%data( 1,kb) = berg%current_point%lon
      pbuff%data( 2,kb) = berg%current_point%lat
      pbuff%data( 3,kb) = berg%current_point%uvel
      pbuff%data( 4,kb) = berg%current_point%vvel
      pbuff%data( 5,kb) = berg%current_point%xi
      pbuff%data( 6,kb) = berg%current_point%yj
      pbuff%data( 7,kb) = float(berg%current_point%year)
      pbuff%data( 8,kb) = berg%current_point%day
      pbuff%data( 9,kb) = berg%current_point%mass
      pbuff%data(10,kb) = berg%current_point%thickness
      pbuff%data(11,kb) = berg%current_point%width
      pbuff%data(12,kb) = berg%current_point%length
      pbuff%data(13,kb) = berg%current_point%mass_of_bits
      pbuff%data(14,kb) = berg%current_point%heat_density

      pbuff%data(15,kb) = berg%mass_scaling
      DO k=1,nkounts
         pbuff%data(15+k,kb) = REAL( berg%number(k), wp )
      END DO
      !
   END SUBROUTINE icb_pack_into_buffer


   SUBROUTINE icb_unpack_from_buffer(first, pbuff, kb)
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      TYPE(iceberg),             POINTER :: first
      TYPE(buffer) ,             POINTER :: pbuff
      INTEGER      , INTENT(in)          :: kb
      ! 
      TYPE(iceberg)                      :: currentberg
      TYPE(point)                        :: pt
      INTEGER                            :: ik
      !!----------------------------------------------------------------------
      !
      pt%lon            =      pbuff%data( 1,kb)
      pt%lat            =      pbuff%data( 2,kb)
      pt%uvel           =      pbuff%data( 3,kb)
      pt%vvel           =      pbuff%data( 4,kb)
      pt%xi             =      pbuff%data( 5,kb)
      pt%yj             =      pbuff%data( 6,kb)
      pt%year           = INT( pbuff%data( 7,kb) )
      pt%day            =      pbuff%data( 8,kb)
      pt%mass           =      pbuff%data( 9,kb)
      pt%thickness      =      pbuff%data(10,kb)
      pt%width          =      pbuff%data(11,kb)
      pt%length         =      pbuff%data(12,kb)
      pt%mass_of_bits   =      pbuff%data(13,kb)
      pt%heat_density   =      pbuff%data(14,kb)

      currentberg%mass_scaling =      pbuff%data(15,kb)
      DO ik = 1, nkounts
         currentberg%number(ik) = INT( pbuff%data(15+ik,kb) )
      END DO
      !
      CALL icb_utl_add(currentberg, pt )
      !
   END SUBROUTINE icb_unpack_from_buffer


   SUBROUTINE icb_increase_buffer(old,kdelta)
      !!----------------------------------------------------------------------
      TYPE(buffer), POINTER    :: old
      INTEGER     , INTENT(in) :: kdelta
      ! 
      TYPE(buffer), POINTER ::   new
      INTEGER ::   inew_size
      !!----------------------------------------------------------------------
      !
      IF( .NOT. ASSOCIATED(old) ) THEN   ;   inew_size = kdelta
      ELSE                               ;   inew_size = old%size + kdelta
      ENDIF
      ALLOCATE( new )
      ALLOCATE( new%data( jp_buffer_width, inew_size) )
      new%size = inew_size
      IF( ASSOCIATED(old) ) THEN
         new%data(:,1:old%size) = old%data(:,1:old%size)
         DEALLOCATE(old%data)
         DEALLOCATE(old)
      ENDIF
      old => new
      !
   END SUBROUTINE icb_increase_buffer


   SUBROUTINE icb_increase_ibuffer(old,kdelta)
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      TYPE(buffer),            POINTER :: old
      INTEGER     , INTENT(in)         :: kdelta
      !
      TYPE(buffer),            POINTER :: new
      INTEGER                          :: inew_size, iold_size
      !!----------------------------------------------------------------------

      IF( .NOT. ASSOCIATED(old) ) THEN
         inew_size = kdelta + jp_delta_buf
         iold_size = 0
      ELSE
         iold_size = old%size
         IF( kdelta .LT. old%size ) THEN
            inew_size = old%size + kdelta
         ELSE
            inew_size = kdelta + jp_delta_buf
         ENDIF
      ENDIF

      IF( iold_size .NE. inew_size ) THEN
         ALLOCATE( new )
         ALLOCATE( new%data( jp_buffer_width, inew_size) )
         new%size = inew_size
         IF( ASSOCIATED(old) ) THEN
            new%data(:,1:old%size) = old%data(:,1:old%size)
            DEALLOCATE(old%data)
            DEALLOCATE(old)
         ENDIF
         old => new
         !IF (nn_verbose_level > 0) WRITE( numicb,*) 'icb_increase_ibuffer',narea,' increased to',inew_size
      ENDIF
      !
   END SUBROUTINE icb_increase_ibuffer

#else
   !!----------------------------------------------------------------------
   !!   Default case:            Dummy module        share memory computing
   !!----------------------------------------------------------------------
   SUBROUTINE icb_lbc_mpp()
      WRITE(numout,*) 'icb_lbc_mpp: You should not have seen this message!!'
   END SUBROUTINE icb_lbc_mpp
#endif

   !!======================================================================
END MODULE icblbc
