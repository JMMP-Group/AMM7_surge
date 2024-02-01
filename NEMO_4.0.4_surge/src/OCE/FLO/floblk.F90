MODULE floblk
   !!======================================================================
   !!                     ***  MODULE  floblk  ***
   !! Ocean floats :   trajectory computation
   !!======================================================================
   !!
   !!----------------------------------------------------------------------
   !!    flotblk     : compute float trajectories with Blanke algorithme
   !!----------------------------------------------------------------------
   USE flo_oce         ! ocean drifting floats
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   flo_blk    ! routine called by floats.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE flo_blk( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE flo_blk  ***
      !!           
      !! ** Purpose :   Compute the geographical position,latitude, longitude
      !!      and depth of each float at each time step.
      !! 
      !! ** Method  :   The position of a float is computed with Bruno Blanke
      !!      algorithm. We need to know the velocity field, the old positions
      !!      of the floats and the grid defined on the domain.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in  ) ::   kt ! ocean time step
      !!
      INTEGER :: jfl              ! dummy loop arguments
      INTEGER :: ind, ifin, iloop
      REAL(wp)   ::       &
         zuinfl,zvinfl,zwinfl,      &     ! transport across the input face
         zuoutfl,zvoutfl,zwoutfl,   &     ! transport across the ouput face
         zvol,                      &     ! volume of the mesh
         zsurfz,                    &     ! surface of the face of the mesh 
         zind

      REAL(wp), DIMENSION ( 2 )  ::   zsurfx, zsurfy   ! surface of the face of the mesh

      INTEGER  , DIMENSION ( jpnfl )  ::   iil, ijl, ikl                   ! index of nearest mesh
      INTEGER  , DIMENSION ( jpnfl )  ::   iiloc , ijloc              
      INTEGER  , DIMENSION ( jpnfl )  ::   iiinfl, ijinfl, ikinfl          ! index of input mesh of the float.
      INTEGER  , DIMENSION ( jpnfl )  ::   iioutfl, ijoutfl, ikoutfl       ! index of output mesh of the float.
      REAL(wp) , DIMENSION ( jpnfl )  ::   zgifl, zgjfl, zgkfl             ! position of floats, index on 
      !                                                                         ! velocity mesh.
      REAL(wp) , DIMENSION ( jpnfl )  ::    ztxfl, ztyfl, ztzfl            ! time for a float to quit the mesh
      !                                                                         ! across one of the face x,y and z 
      REAL(wp) , DIMENSION ( jpnfl )  ::    zttfl                          ! time for a float to quit the mesh 
      REAL(wp) , DIMENSION ( jpnfl )  ::    zagefl                         ! time during which, trajectorie of 
      !                                                                         ! the float has been computed
      REAL(wp) , DIMENSION ( jpnfl )  ::   zagenewfl                       ! new age of float after calculation 
      !                                                                         ! of new position
      REAL(wp) , DIMENSION ( jpnfl )  ::   zufl, zvfl, zwfl                ! interpolated vel. at float position
      REAL(wp) , DIMENSION ( jpnfl )  ::   zudfl, zvdfl, zwdfl             ! velocity diff input/output of mesh
      REAL(wp) , DIMENSION ( jpnfl )  ::   zgidfl, zgjdfl, zgkdfl          ! direction index of float
      !!---------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'flo_blk : compute Blanke trajectories for floats '
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF

      ! Initialisation of parameters
      
      DO jfl = 1, jpnfl
         ! ages of floats are put at zero
         zagefl(jfl) = 0.
         ! index on the velocity grid 
         ! We considere k coordinate negative, with this transformation 
         ! the computation in the 3 direction is the same. 
         zgifl(jfl) = tpifl(jfl) - 0.5
         zgjfl(jfl) = tpjfl(jfl) - 0.5
         zgkfl(jfl) = MIN(-1.,-(tpkfl(jfl)))
         ! surface drift every 10 days 
         IF( ln_argo ) THEN
            IF( MOD(kt,150) >= 146 .OR. MOD(kt,150) == 0 )  zgkfl(jfl) = -1.
         ENDIF
         ! index of T mesh
         iil(jfl) = 1 + INT(zgifl(jfl))
         ijl(jfl) = 1 + INT(zgjfl(jfl))
         ikl(jfl) =     INT(zgkfl(jfl))
      END DO
       
      iloop = 0
222   DO jfl = 1, jpnfl
# if   defined key_mpp_mpi
         IF( iil(jfl) >= mig(nldi) .AND. iil(jfl) <= mig(nlei) .AND.   &
             ijl(jfl) >= mjg(nldj) .AND. ijl(jfl) <= mjg(nlej)   ) THEN
            iiloc(jfl) = iil(jfl) - mig(1) + 1
            ijloc(jfl) = ijl(jfl) - mjg(1) + 1
# else 
            iiloc(jfl) = iil(jfl)
            ijloc(jfl) = ijl(jfl)
# endif
            
            ! compute the transport across the mesh where the float is.            
!!bug (gm) change e3t into e3. but never checked 
            zsurfx(1) = e2u(iiloc(jfl)-1,ijloc(jfl)  ) * e3u_n(iiloc(jfl)-1,ijloc(jfl)  ,-ikl(jfl))
            zsurfx(2) = e2u(iiloc(jfl)  ,ijloc(jfl)  ) * e3u_n(iiloc(jfl)  ,ijloc(jfl)  ,-ikl(jfl))
            zsurfy(1) = e1v(iiloc(jfl)  ,ijloc(jfl)-1) * e3v_n(iiloc(jfl)  ,ijloc(jfl)-1,-ikl(jfl))
            zsurfy(2) = e1v(iiloc(jfl)  ,ijloc(jfl)  ) * e3v_n(iiloc(jfl)  ,ijloc(jfl)  ,-ikl(jfl))

            ! for a isobar float zsurfz is put to zero. The vertical velocity will be zero too.
            zsurfz =          e1e2t(iiloc(jfl),ijloc(jfl))
            zvol   = zsurfz * e3t_n(iiloc(jfl),ijloc(jfl),-ikl(jfl))

            !
            zuinfl =( ub(iiloc(jfl)-1,ijloc(jfl),-ikl(jfl)) + un(iiloc(jfl)-1,ijloc(jfl),-ikl(jfl)) )/2.*zsurfx(1)
            zuoutfl=( ub(iiloc(jfl)  ,ijloc(jfl),-ikl(jfl)) + un(iiloc(jfl)  ,ijloc(jfl),-ikl(jfl)) )/2.*zsurfx(2)
            zvinfl =( vb(iiloc(jfl),ijloc(jfl)-1,-ikl(jfl)) + vn(iiloc(jfl),ijloc(jfl)-1,-ikl(jfl)) )/2.*zsurfy(1)
            zvoutfl=( vb(iiloc(jfl),ijloc(jfl)  ,-ikl(jfl)) + vn(iiloc(jfl),ijloc(jfl)  ,-ikl(jfl)) )/2.*zsurfy(2)
            zwinfl =-(wb(iiloc(jfl),ijloc(jfl),-(ikl(jfl)-1))    &
               &   +  wn(iiloc(jfl),ijloc(jfl),-(ikl(jfl)-1)) )/2. *  zsurfz*nisobfl(jfl)
            zwoutfl=-(wb(iiloc(jfl),ijloc(jfl),- ikl(jfl)   )   &
               &   +  wn(iiloc(jfl),ijloc(jfl),- ikl(jfl)   ) )/2. *  zsurfz*nisobfl(jfl)
            
            ! interpolation of velocity field on the float initial position            
            zufl(jfl)=  zuinfl  + ( zgifl(jfl) - float(iil(jfl)-1) ) * ( zuoutfl - zuinfl)
            zvfl(jfl)=  zvinfl  + ( zgjfl(jfl) - float(ijl(jfl)-1) ) * ( zvoutfl - zvinfl)
            zwfl(jfl)=  zwinfl  + ( zgkfl(jfl) - float(ikl(jfl)-1) ) * ( zwoutfl - zwinfl)
            
            ! faces of input and output
            ! u-direction
            IF( zufl(jfl) < 0. ) THEN
               iioutfl(jfl) = iil(jfl) - 1.
               iiinfl (jfl) = iil(jfl)
               zind   = zuinfl
               zuinfl = zuoutfl
               zuoutfl= zind
            ELSE
               iioutfl(jfl) = iil(jfl)
               iiinfl (jfl) = iil(jfl) - 1
            ENDIF
            ! v-direction       
            IF( zvfl(jfl) < 0. ) THEN
               ijoutfl(jfl) = ijl(jfl) - 1.
               ijinfl (jfl) = ijl(jfl)
               zind    = zvinfl
               zvinfl  = zvoutfl
               zvoutfl = zind
            ELSE
               ijoutfl(jfl) = ijl(jfl)
               ijinfl (jfl) = ijl(jfl) - 1.
            ENDIF
            ! w-direction
            IF( zwfl(jfl) < 0. ) THEN
               ikoutfl(jfl) = ikl(jfl) - 1.
               ikinfl (jfl) = ikl(jfl)
               zind    = zwinfl
               zwinfl  = zwoutfl
               zwoutfl = zind
            ELSE
               ikoutfl(jfl) = ikl(jfl)
               ikinfl (jfl) = ikl(jfl) - 1.
            ENDIF
            
            ! compute the time to go out the mesh across a face
            ! u-direction
            zudfl (jfl) = zuoutfl - zuinfl
            zgidfl(jfl) = float(iioutfl(jfl) - iiinfl(jfl))
            IF( zufl(jfl)*zuoutfl <= 0. ) THEN
               ztxfl(jfl) = 1.E99
            ELSE
               IF( ABS(zudfl(jfl)) >= 1.E-5 ) THEN
                  ztxfl(jfl)= zgidfl(jfl)/zudfl(jfl) * LOG(zuoutfl/zufl (jfl))
               ELSE
                  ztxfl(jfl)=(float(iioutfl(jfl))-zgifl(jfl))/zufl(jfl)
               ENDIF
               IF( (ABS(zgifl(jfl)-float(iiinfl (jfl))) <=  1.E-7) .OR.   &
                   (ABS(zgifl(jfl)-float(iioutfl(jfl))) <=  1.E-7) ) THEN
                  ztxfl(jfl)=(zgidfl(jfl))/zufl(jfl)
               ENDIF
            ENDIF
            ! v-direction
            zvdfl (jfl) = zvoutfl - zvinfl
            zgjdfl(jfl) = float(ijoutfl(jfl)-ijinfl(jfl))
            IF( zvfl(jfl)*zvoutfl <= 0. ) THEN
               ztyfl(jfl) = 1.E99
            ELSE
               IF( ABS(zvdfl(jfl)) >= 1.E-5 ) THEN
                  ztyfl(jfl) = zgjdfl(jfl)/zvdfl(jfl) * LOG(zvoutfl/zvfl (jfl))
               ELSE
                  ztyfl(jfl) = (float(ijoutfl(jfl)) - zgjfl(jfl))/zvfl(jfl)
               ENDIF
               IF( (ABS(zgjfl(jfl)-float(ijinfl (jfl))) <= 1.E-7) .OR.   &
                   (ABS(zgjfl(jfl)-float(ijoutfl(jfl))) <=  1.E-7) ) THEN
                  ztyfl(jfl) = (zgjdfl(jfl)) / zvfl(jfl)
               ENDIF
            ENDIF
            ! w-direction        
            IF( nisobfl(jfl) == 1. ) THEN 
               zwdfl (jfl) = zwoutfl - zwinfl
               zgkdfl(jfl) = float(ikoutfl(jfl) - ikinfl(jfl))
               IF( zwfl(jfl)*zwoutfl <= 0. ) THEN
                  ztzfl(jfl) = 1.E99
               ELSE
                  IF( ABS(zwdfl(jfl)) >= 1.E-5 ) THEN
                     ztzfl(jfl) = zgkdfl(jfl)/zwdfl(jfl) * LOG(zwoutfl/zwfl (jfl))
                  ELSE
                     ztzfl(jfl) = (float(ikoutfl(jfl)) - zgkfl(jfl))/zwfl(jfl)
                  ENDIF
                  IF( (ABS(zgkfl(jfl)-float(ikinfl (jfl))) <=  1.E-7) .OR.   &
                      (ABS(zgkfl(jfl)-float(ikoutfl(jfl))) <= 1.E-7) ) THEN
                     ztzfl(jfl) = (zgkdfl(jfl)) / zwfl(jfl)
                  ENDIF
               ENDIF
            ENDIF
            
            ! the time to go leave the mesh is the smallest time
                   
            IF( nisobfl(jfl) == 1. ) THEN 
               zttfl(jfl) = MIN(ztxfl(jfl),ztyfl(jfl),ztzfl(jfl))
            ELSE
               zttfl(jfl) = MIN(ztxfl(jfl),ztyfl(jfl))
            ENDIF
            ! new age of the FLOAT
            zagenewfl(jfl) = zagefl(jfl) + zttfl(jfl)*zvol
            ! test to know if the "age" of the float is not bigger than the 
            ! time step
            IF( zagenewfl(jfl) > rdt ) THEN
               zttfl(jfl) = (rdt-zagefl(jfl)) / zvol
               zagenewfl(jfl) = rdt
            ENDIF
            
            ! In the "minimal" direction we compute the index of new mesh
            ! on i-direction
            IF( ztxfl(jfl) <=  zttfl(jfl) ) THEN
               zgifl(jfl) = float(iioutfl(jfl))
               ind = iioutfl(jfl)
               IF( iioutfl(jfl) >= iiinfl(jfl) ) THEN
                  iioutfl(jfl) = iioutfl(jfl) + 1
               ELSE
                  iioutfl(jfl) = iioutfl(jfl) - 1
               ENDIF
               iiinfl(jfl) = ind
            ELSE
               IF( ABS(zudfl(jfl)) >= 1.E-5 ) THEN 
                  zgifl(jfl) = zgifl(jfl) + zgidfl(jfl)*zufl(jfl)    &
                     &       * ( EXP( zudfl(jfl)/zgidfl(jfl)*zttfl(jfl) ) - 1. ) /  zudfl(jfl)
               ELSE
                  zgifl(jfl) = zgifl(jfl) + zufl(jfl) * zttfl(jfl)
               ENDIF
            ENDIF
            ! on j-direction
            IF( ztyfl(jfl) <= zttfl(jfl) ) THEN
               zgjfl(jfl) = float(ijoutfl(jfl))
               ind = ijoutfl(jfl)
               IF( ijoutfl(jfl) >= ijinfl(jfl) ) THEN
                  ijoutfl(jfl) = ijoutfl(jfl) + 1
               ELSE
                  ijoutfl(jfl) = ijoutfl(jfl) - 1
               ENDIF
               ijinfl(jfl) = ind
            ELSE
               IF( ABS(zvdfl(jfl)) >= 1.E-5 ) THEN 
                  zgjfl(jfl) = zgjfl(jfl)+zgjdfl(jfl)*zvfl(jfl)   &
                     &       * ( EXP(zvdfl(jfl)/zgjdfl(jfl)*zttfl(jfl)) - 1. ) /  zvdfl(jfl)
               ELSE
                  zgjfl(jfl) = zgjfl(jfl)+zvfl(jfl)*zttfl(jfl)
               ENDIF
            ENDIF
            ! on k-direction
            IF( nisobfl(jfl) == 1. ) THEN
               IF( ztzfl(jfl) <= zttfl(jfl) ) THEN
                  zgkfl(jfl) = float(ikoutfl(jfl))
                  ind = ikoutfl(jfl)
                  IF( ikoutfl(jfl) >= ikinfl(jfl) ) THEN
                     ikoutfl(jfl) = ikoutfl(jfl)+1
                  ELSE
                     ikoutfl(jfl) = ikoutfl(jfl)-1
                  ENDIF
                  ikinfl(jfl) = ind
               ELSE
                  IF( ABS(zwdfl(jfl)) >= 1.E-5 ) THEN 
                     zgkfl(jfl) = zgkfl(jfl)+zgkdfl(jfl)*zwfl(jfl)    &
                        &       * ( EXP(zwdfl(jfl)/zgkdfl(jfl)*zttfl(jfl)) - 1. ) /  zwdfl(jfl)
                  ELSE
                     zgkfl(jfl) = zgkfl(jfl)+zwfl(jfl)*zttfl(jfl)
                  ENDIF
               ENDIF
            ENDIF
            
            ! coordinate of the new point on the temperature grid
            
            iil(jfl) = MAX(iiinfl(jfl),iioutfl(jfl))
            ijl(jfl) = MAX(ijinfl(jfl),ijoutfl(jfl))
            IF( nisobfl(jfl) ==  1 ) ikl(jfl) = MAX(ikinfl(jfl),ikoutfl(jfl))
!!Alexcadm	 write(*,*)'PE ',narea,
!!Alexcadm     .    iiinfl(jfl),iioutfl(jfl),ijinfl(jfl)
!!Alexcadm     .		,ijoutfl(jfl),ikinfl(jfl),
!!Alexcadm     .	  ikoutfl(jfl),ztxfl(jfl),ztyfl(jfl)
!!Alexcadm     .		,ztzfl(jfl),zgifl(jfl),
!!Alexcadm     .  zgjfl(jfl) 
!!Alexcadm	IF (jfl == 910) write(*,*)'Flotteur 910',
!!Alexcadm     .    iiinfl(jfl),iioutfl(jfl),ijinfl(jfl)
!!Alexcadm     .		,ijoutfl(jfl),ikinfl(jfl),
!!Alexcadm     .	  ikoutfl(jfl),ztxfl(jfl),ztyfl(jfl)
!!Alexcadm     .		,ztzfl(jfl),zgifl(jfl),
!!Alexcadm     .  zgjfl(jfl) 
            ! reinitialisation of the age of FLOAT
            zagefl(jfl) = zagenewfl(jfl)
# if   defined key_mpp_mpi
         ELSE
            ! we put zgifl, zgjfl, zgkfl, zagefl
            zgifl (jfl) = 0.
            zgjfl (jfl) = 0.
            zgkfl (jfl) = 0.
            zagefl(jfl) = 0.
            iil(jfl) = 0
            ijl(jfl) = 0
         ENDIF
# endif
      END DO
      
      ! synchronisation
      CALL mpp_sum( 'floblk', zgifl , jpnfl )   ! sums over the global domain
      CALL mpp_sum( 'floblk', zgjfl , jpnfl )
      CALL mpp_sum( 'floblk', zgkfl , jpnfl )
      CALL mpp_sum( 'floblk', zagefl, jpnfl )
      CALL mpp_sum( 'floblk', iil   , jpnfl )
      CALL mpp_sum( 'floblk', ijl   , jpnfl )
      
      ! Test to know if a  float hasn't integrated enought time
      IF( ln_argo ) THEN
         ifin = 1
         DO jfl = 1, jpnfl
            IF( zagefl(jfl) < rdt )   ifin = 0
            tpifl(jfl) = zgifl(jfl) + 0.5
            tpjfl(jfl) = zgjfl(jfl) + 0.5
         END DO
      ELSE
         ifin = 1
         DO jfl = 1, jpnfl
            IF( zagefl(jfl) < rdt )   ifin = 0
            tpifl(jfl) = zgifl(jfl) + 0.5
            tpjfl(jfl) = zgjfl(jfl) + 0.5
            IF( nisobfl(jfl) == 1 ) tpkfl(jfl) = -(zgkfl(jfl))
         END DO
      ENDIF
!!Alexcadm	IF (lwp) write(numout,*) '---------'
!!Alexcadm	IF (lwp) write(numout,*) 'before Erika:',tpifl(880),tpjfl(880),
!!Alexcadm     .       tpkfl(880),zufl(880),zvfl(880),zwfl(880)
!!Alexcadm	IF (lwp) write(numout,*) 'first Erika:',tpifl(900),tpjfl(900),
!!Alexcadm     .       tpkfl(900),zufl(900),zvfl(900),zwfl(900)
!!Alexcadm	IF (lwp) write(numout,*) 'last Erika:',tpifl(jpnfl),tpjfl(jpnfl),
!!Alexcadm     .       tpkfl(jpnfl),zufl(jpnfl),zvfl(jpnfl),zwfl(jpnfl)
      IF( ifin == 0 ) THEN
         iloop = iloop + 1 
         GO TO 222
      ENDIF
      !
      !
   END SUBROUTINE flo_blk

   !!======================================================================
END MODULE floblk 
