MODULE flo4rk
   !!======================================================================
   !!                    ***  MODULE  flo4rk  ***
   !! Ocean floats :   trajectory computation using a 4th order Runge-Kutta
   !!======================================================================
   !!
   !!----------------------------------------------------------------------
   !!   flo_4rk        : Compute the geographical position of floats
   !!   flo_interp     : interpolation
   !!----------------------------------------------------------------------
   USE flo_oce         ! ocean drifting floats
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   flo_4rk    ! routine called by floats.F90

   !                                   ! RK4 and Lagrange interpolation coefficients
   REAL(wp), DIMENSION (4) ::   tcoef1 = (/  1.0  ,  0.5  ,  0.5  ,  0.0  /)   ! 
   REAL(wp), DIMENSION (4) ::   tcoef2 = (/  0.0  ,  0.5  ,  0.5  ,  1.0  /)   !
   REAL(wp), DIMENSION (4) ::   scoef2 = (/  1.0  ,  2.0  ,  2.0  ,  1.0  /)   !
   REAL(wp), DIMENSION (4) ::   rcoef  = (/-1./6. , 1./2. ,-1./2. , 1./6. /)   !
   REAL(wp), DIMENSION (3) ::   scoef1 = (/  0.5  ,  0.5  ,  1.0  /)           !

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE flo_4rk( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE flo_4rk  ***
      !!
      !!  ** Purpose :   Compute the geographical position (lat,lon,depth)
      !!       of each float at each time step.
      !! 
      !!  ** Method  :   The position of a float is computed with a 4th order
      !!       Runge-Kutta scheme and and Lagrange interpolation.
      !!         We need to know the velocity field, the old positions of the
      !!       floats and the grid defined on the domain.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER ::  jfl, jind           ! dummy loop indices
      INTEGER ::  ierror              ! error value

      REAL(wp), DIMENSION(jpnfl)   ::   zgifl , zgjfl , zgkfl    ! index RK  positions
      REAL(wp), DIMENSION(jpnfl)   ::   zufl  , zvfl  , zwfl     ! interpolated velocity at the float position 
      REAL(wp), DIMENSION(jpnfl,4) ::   zrkxfl, zrkyfl, zrkzfl   ! RK coefficients
      !!---------------------------------------------------------------------
      !
      IF( ierror /= 0 ) THEN
         WRITE(numout,*) 'flo_4rk: allocation of workspace arrays failed'
      ENDIF

    
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'flo_4rk : compute Runge Kutta trajectories for floats '
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

      ! Verification of the floats positions. If one of them leave the domain
      ! domain we replace the float near the border.
      DO jfl = 1, jpnfl
         ! i-direction
         IF( tpifl(jfl) <= 1.5 ) THEN
            IF(lwp)WRITE(numout,*)'!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!'
            IF(lwp)WRITE(numout,*)'The float',jfl,'is out of the domain at the WEST border.'
            tpifl(jfl) = tpifl(jfl) + 1.
            IF(lwp)WRITE(numout,*)'New initialisation for this float at i=',tpifl(jfl)
         ENDIF
          
         IF( tpifl(jfl) >= jpi-.5 ) THEN
            IF(lwp)WRITE(numout,*)'!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!'
            IF(lwp)WRITE(numout,*)'The float',jfl,'is out of the domain at the EAST border.'
            tpifl(jfl) = tpifl(jfl) - 1.
            IF(lwp)WRITE(numout,*)'New initialisation for this float at i=', tpifl(jfl)
         ENDIF
         ! j-direction
         IF( tpjfl(jfl) <= 1.5 ) THEN
            IF(lwp)WRITE(numout,*)'!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!'
            IF(lwp)WRITE(numout,*)'The float',jfl,'is out of the domain at the SOUTH border.'
            tpjfl(jfl) = tpjfl(jfl) + 1.
            IF(lwp)WRITE(numout,*)'New initialisation for this float at j=', tpjfl(jfl)
         ENDIF
           
         IF( tpjfl(jfl) >= jpj-.5 ) THEN
            IF(lwp)WRITE(numout,*)'!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!'
            IF(lwp)WRITE(numout,*)'The float',jfl,'is out of the domain at the NORTH border.'
            tpjfl(jfl) = tpjfl(jfl) - 1.
            IF(lwp)WRITE(numout,*)'New initialisation for this float at j=', tpjfl(jfl)
         ENDIF
         ! k-direction
         IF( tpkfl(jfl) <= .5 ) THEN
            IF(lwp)WRITE(numout,*)'!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!'
            IF(lwp)WRITE(numout,*)'The float',jfl,'is out of the domain at the TOP border.'
            tpkfl(jfl) = tpkfl(jfl) + 1.
            IF(lwp)WRITE(numout,*)'New initialisation for this float at k=', tpkfl(jfl)
         ENDIF
         
         IF( tpkfl(jfl) >= jpk-.5 )  THEN
            IF(lwp)WRITE(numout,*)'!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!'
            IF(lwp)WRITE(numout,*)'The float',jfl,'is out of the domain at the BOTTOM border.'
            tpkfl(jfl) = tpkfl(jfl) - 1.
            IF(lwp)WRITE(numout,*)'New initialisation for this float at k=', tpkfl(jfl)
         ENDIF
      END DO
      
      ! 4 steps of Runge-Kutta algorithme
      ! initialisation of the positions 
      
      DO jfl = 1, jpnfl
         zgifl(jfl) = tpifl(jfl)
         zgjfl(jfl) = tpjfl(jfl)
         zgkfl(jfl) = tpkfl(jfl)
      END DO
       
      DO  jind = 1, 4         
      
         ! for each step we compute the compute the velocity with Lagrange interpolation
         CALL flo_interp( zgifl, zgjfl, zgkfl, zufl, zvfl, zwfl, jind )
         
         ! computation of Runge-Kutta factor
         DO jfl = 1, jpnfl
            zrkxfl(jfl,jind) = rdt*zufl(jfl)
            zrkyfl(jfl,jind) = rdt*zvfl(jfl)
            zrkzfl(jfl,jind) = rdt*zwfl(jfl)
         END DO
         IF( jind /= 4 ) THEN
            DO jfl = 1, jpnfl
               zgifl(jfl) = (tpifl(jfl)) + scoef1(jind)*zrkxfl(jfl,jind)
               zgjfl(jfl) = (tpjfl(jfl)) + scoef1(jind)*zrkyfl(jfl,jind)
               zgkfl(jfl) = (tpkfl(jfl)) + scoef1(jind)*zrkzfl(jfl,jind)
            END DO
         ENDIF
      END DO
      DO jind = 1, 4
         DO jfl = 1, jpnfl
            tpifl(jfl) = tpifl(jfl) + scoef2(jind)*zrkxfl(jfl,jind)/6.
            tpjfl(jfl) = tpjfl(jfl) + scoef2(jind)*zrkyfl(jfl,jind)/6.
            tpkfl(jfl) = tpkfl(jfl) + scoef2(jind)*zrkzfl(jfl,jind)/6.
         END DO
      END DO
      !
      !
   END SUBROUTINE flo_4rk


   SUBROUTINE flo_interp( pxt , pyt , pzt ,      &
      &                   pufl, pvfl, pwfl, ki )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE flointerp  ***
      !!
      !! ** Purpose :   Interpolation of the velocity on the float position
      !! 
      !! ** Method  :   Lagrange interpolation with the 64 neighboring
      !!      points. This routine is call 4 time at each time step to
      !!      compute velocity at the date and the position we need to
      !!      integrated with RK method.
      !!----------------------------------------------------------------------
      REAL(wp) , DIMENSION(jpnfl), INTENT(in   ) ::   pxt , pyt , pzt    ! position of the float
      REAL(wp) , DIMENSION(jpnfl), INTENT(  out) ::   pufl, pvfl, pwfl   ! velocity at this position
      INTEGER                    , INTENT(in   ) ::   ki                 !
      !!
      INTEGER  ::   jfl, jind1, jind2, jind3   ! dummy loop indices
      REAL(wp) ::   zsumu, zsumv, zsumw        ! local scalar
      INTEGER  , DIMENSION(jpnfl)       ::   iilu, ijlu, iklu   ! nearest neighbour INDEX-u
      INTEGER  , DIMENSION(jpnfl)       ::   iilv, ijlv, iklv   ! nearest neighbour INDEX-v
      INTEGER  , DIMENSION(jpnfl)       ::   iilw, ijlw, iklw   ! nearest neighbour INDEX-w
      INTEGER  , DIMENSION(jpnfl,4)     ::   iidu, ijdu, ikdu   ! 64 nearest neighbour INDEX-u
      INTEGER  , DIMENSION(jpnfl,4)     ::   iidv, ijdv, ikdv   ! 64 nearest neighbour INDEX-v
      INTEGER  , DIMENSION(jpnfl,4)     ::   iidw, ijdw, ikdw   ! 64 nearest neighbour INDEX-w
      REAL(wp) , DIMENSION(jpnfl,4)     ::   zlagxu, zlagyu, zlagzu   ! Lagrange  coefficients
      REAL(wp) , DIMENSION(jpnfl,4)     ::   zlagxv, zlagyv, zlagzv   !    -           -
      REAL(wp) , DIMENSION(jpnfl,4)     ::   zlagxw, zlagyw, zlagzw   !    -           -
      REAL(wp) , DIMENSION(jpnfl,4,4,4) ::   ztufl , ztvfl , ztwfl   ! velocity at choosen time step
      !!---------------------------------------------------------------------
 
      ! Interpolation of U velocity

      ! nearest neightboring point for computation of u       
      DO jfl = 1, jpnfl
         iilu(jfl) = INT(pxt(jfl)-.5)
         ijlu(jfl) = INT(pyt(jfl)-.5)
         iklu(jfl) = INT(pzt(jfl))
      END DO
      
      !  64 neightboring points for computation of u 
      DO jind1 = 1, 4
         DO jfl = 1, jpnfl
            !  i-direction
            IF( iilu(jfl) <= 2 ) THEN          ;   iidu(jfl,jind1) = jind1
            ELSE
               IF( iilu(jfl) >= jpi-1 ) THEN   ;   iidu(jfl,jind1) = jpi       + jind1 - 4
               ELSE                            ;   iidu(jfl,jind1) = iilu(jfl) + jind1 - 2
               ENDIF
            ENDIF
            !  j-direction
            IF( ijlu(jfl) <= 2 ) THEN          ;   ijdu(jfl,jind1) = jind1
            ELSE
               IF( ijlu(jfl) >= jpj-1 ) THEN   ;   ijdu(jfl,jind1) = jpj       + jind1 - 4
               ELSE                            ;   ijdu(jfl,jind1) = ijlu(jfl) + jind1 - 2
               ENDIF
            ENDIF
            ! k-direction
            IF( iklu(jfl) <= 2 ) THEN          ;   ikdu(jfl,jind1) = jind1
            ELSE
               IF( iklu(jfl) >= jpk-1 ) THEN   ;   ikdu(jfl,jind1) = jpk + jind1 - 4
               ELSE                            ;   ikdu(jfl,jind1) = iklu(jfl) + jind1 - 2
               ENDIF
            ENDIF
         END DO
      END DO
      
      ! Lagrange coefficients
      DO jfl = 1, jpnfl
         DO jind1 = 1, 4
            zlagxu(jfl,jind1) = 1.
            zlagyu(jfl,jind1) = 1.
            zlagzu(jfl,jind1) = 1.
         END DO
      END DO
      DO jind1 = 1, 4
         DO jind2 = 1, 4
            DO jfl= 1, jpnfl
               IF( jind1 /= jind2 ) THEN
                  zlagxu(jfl,jind1) = zlagxu(jfl,jind1) * ( pxt(jfl)-(float(iidu(jfl,jind2))+.5) )
                  zlagyu(jfl,jind1) = zlagyu(jfl,jind1) * ( pyt(jfl)-(float(ijdu(jfl,jind2))) )
                  zlagzu(jfl,jind1) = zlagzu(jfl,jind1) * ( pzt(jfl)-(float(ikdu(jfl,jind2))) )
               ENDIF
            END DO
         END DO
      END DO
      
      ! velocity when we compute at middle time step
      
      DO jfl = 1, jpnfl
         DO jind1 = 1, 4
            DO jind2 = 1, 4
               DO jind3 = 1, 4
                  ztufl(jfl,jind1,jind2,jind3) =   &
                     &   (  tcoef1(ki) * ub(iidu(jfl,jind1),ijdu(jfl,jind2),ikdu(jfl,jind3)) +   &
                     &      tcoef2(ki) * un(iidu(jfl,jind1),ijdu(jfl,jind2),ikdu(jfl,jind3)) )   &
                     &      / e1u(iidu(jfl,jind1),ijdu(jfl,jind2)) 
               END DO
            END DO
         END DO
         
         zsumu = 0.
         DO jind1 = 1, 4
            DO jind2 = 1, 4
               DO jind3 = 1, 4
                  zsumu = zsumu + ztufl(jfl,jind1,jind2,jind3) * zlagxu(jfl,jind1) * zlagyu(jfl,jind2)   &
                     &  * zlagzu(jfl,jind3) * rcoef(jind1)*rcoef(jind2)*rcoef(jind3) 
               END DO
            END DO
         END DO
         pufl(jfl) = zsumu
      END DO
      
      ! Interpolation of V velocity 

      ! nearest neightboring point for computation of v 
      DO jfl = 1, jpnfl
         iilv(jfl) = INT(pxt(jfl)-.5)
         ijlv(jfl) = INT(pyt(jfl)-.5)
         iklv(jfl) = INT(pzt(jfl))
      END DO
      
      ! 64 neightboring points for computation of v 
      DO jind1 = 1, 4
         DO jfl = 1, jpnfl
            ! i-direction
            IF( iilv(jfl) <= 2 ) THEN          ;   iidv(jfl,jind1) = jind1
            ELSE
               IF( iilv(jfl) >= jpi-1 ) THEN   ;   iidv(jfl,jind1) = jpi       + jind1 - 4
               ELSE                            ;   iidv(jfl,jind1) = iilv(jfl) + jind1 - 2
               ENDIF
            ENDIF
            ! j-direction
            IF( ijlv(jfl) <= 2 ) THEN          ;   ijdv(jfl,jind1) = jind1
            ELSE
               IF( ijlv(jfl) >= jpj-1 ) THEN   ;   ijdv(jfl,jind1) = jpj       + jind1 - 4
               ELSE                            ;   ijdv(jfl,jind1) = ijlv(jfl) + jind1 - 2
               ENDIF
            ENDIF
            ! k-direction
            IF( iklv(jfl) <= 2 ) THEN          ;   ikdv(jfl,jind1) = jind1
            ELSE
               IF( iklv(jfl) >= jpk-1 ) THEN   ;   ikdv(jfl,jind1) = jpk + jind1 - 4
               ELSE                            ;   ikdv(jfl,jind1) = iklv(jfl) + jind1 - 2
               ENDIF
            ENDIF
         END DO
      END DO
      
      ! Lagrange coefficients
      
      DO jfl = 1, jpnfl
         DO jind1 = 1, 4
            zlagxv(jfl,jind1) = 1.
            zlagyv(jfl,jind1) = 1.
            zlagzv(jfl,jind1) = 1.
         END DO
      END DO
      
      DO jind1 = 1, 4
         DO jind2 = 1, 4
            DO jfl = 1, jpnfl
               IF( jind1 /= jind2 ) THEN
                  zlagxv(jfl,jind1)= zlagxv(jfl,jind1)*(pxt(jfl) - (float(iidv(jfl,jind2))   ) )
                  zlagyv(jfl,jind1)= zlagyv(jfl,jind1)*(pyt(jfl) - (float(ijdv(jfl,jind2))+.5) )
                  zlagzv(jfl,jind1)= zlagzv(jfl,jind1)*(pzt(jfl) - (float(ikdv(jfl,jind2))   ) )
               ENDIF
            END DO
         END DO
      END DO
      
      ! velocity when we compute at middle time step
      
      DO jfl = 1, jpnfl
         DO jind1 = 1, 4
            DO jind2 = 1, 4
               DO jind3 = 1 ,4
                  ztvfl(jfl,jind1,jind2,jind3)=   &
                     &   ( tcoef1(ki) * vb(iidv(jfl,jind1),ijdv(jfl,jind2),ikdv(jfl,jind3))  +   &
                     &     tcoef2(ki) * vn(iidv(jfl,jind1),ijdv(jfl,jind2),ikdv(jfl,jind3)) )    & 
                     &     / e2v(iidv(jfl,jind1),ijdv(jfl,jind2))
               END DO
            END DO
         END DO
         
         zsumv=0.
         DO jind1 = 1, 4
            DO jind2 = 1, 4
               DO jind3 = 1, 4
                  zsumv = zsumv + ztvfl(jfl,jind1,jind2,jind3) * zlagxv(jfl,jind1) * zlagyv(jfl,jind2)   &
                     &  * zlagzv(jfl,jind3) * rcoef(jind1)*rcoef(jind2)*rcoef(jind3)
               END DO
            END DO
         END DO
         pvfl(jfl) = zsumv
      END DO
      
      ! Interpolation of W velocity

      ! nearest neightboring point for computation of w 
      DO jfl = 1, jpnfl
         iilw(jfl) = INT( pxt(jfl)   )
         ijlw(jfl) = INT( pyt(jfl)   )
         iklw(jfl) = INT( pzt(jfl)+.5)
      END DO
      
      ! 64 neightboring points for computation of w 
      DO jind1 = 1, 4
         DO jfl = 1, jpnfl
            ! i-direction
            IF( iilw(jfl) <= 2 ) THEN          ;   iidw(jfl,jind1) = jind1
            ELSE
               IF( iilw(jfl) >= jpi-1 ) THEN   ;   iidw(jfl,jind1) = jpi       + jind1 - 4
               ELSE                            ;   iidw(jfl,jind1) = iilw(jfl) + jind1 - 2
               ENDIF
            ENDIF
            ! j-direction
            IF( ijlw(jfl) <= 2 ) THEN          ;   ijdw(jfl,jind1) = jind1
            ELSE
               IF( ijlw(jfl) >= jpj-1 ) THEN   ;   ijdw(jfl,jind1) = jpj       + jind1 - 4
               ELSE                            ;   ijdw(jfl,jind1) = ijlw(jfl) + jind1 - 2
               ENDIF
            ENDIF
            ! k-direction
            IF( iklw(jfl) <= 2 ) THEN          ;   ikdw(jfl,jind1) = jind1
            ELSE
               IF( iklw(jfl) >= jpk-1 ) THEN   ;   ikdw(jfl,jind1) = jpk       + jind1 - 4
               ELSE                            ;   ikdw(jfl,jind1) = iklw(jfl) + jind1 - 2
               ENDIF
            ENDIF
         END DO
      END DO
      DO jind1 = 1, 4
         DO jfl = 1, jpnfl
            IF( iklw(jfl) <= 2 ) THEN          ;   ikdw(jfl,jind1) = jind1
            ELSE
               IF( iklw(jfl) >= jpk-1 ) THEN   ;   ikdw(jfl,jind1) = jpk       + jind1 - 4
               ELSE                            ;   ikdw(jfl,jind1) = iklw(jfl) + jind1 - 2
               ENDIF
            ENDIF
         END DO
      END DO
      
      ! Lagrange coefficients  for w interpolation
      DO jfl = 1, jpnfl
         DO jind1 = 1, 4
            zlagxw(jfl,jind1) = 1.
            zlagyw(jfl,jind1) = 1.
            zlagzw(jfl,jind1) = 1.
         END DO
      END DO
      DO jind1 = 1, 4
         DO jind2 = 1, 4
            DO jfl = 1, jpnfl
               IF( jind1 /= jind2 ) THEN
                  zlagxw(jfl,jind1) = zlagxw(jfl,jind1) * (pxt(jfl) - (float(iidw(jfl,jind2))   ) )
                  zlagyw(jfl,jind1) = zlagyw(jfl,jind1) * (pyt(jfl) - (float(ijdw(jfl,jind2))   ) )
                  zlagzw(jfl,jind1) = zlagzw(jfl,jind1) * (pzt(jfl) - (float(ikdw(jfl,jind2))-.5) )
               ENDIF
            END DO
         END DO
      END DO
      
      ! velocity w  when we compute at middle time step
      DO jfl = 1, jpnfl
         DO jind1 = 1, 4
            DO jind2 = 1, 4
               DO jind3 = 1, 4
                  ztwfl(jfl,jind1,jind2,jind3)=   &
                     &   ( tcoef1(ki) * wb(iidw(jfl,jind1),ijdw(jfl,jind2),ikdw(jfl,jind3))+   &
                     &     tcoef2(ki) * wn(iidw(jfl,jind1),ijdw(jfl,jind2),ikdw(jfl,jind3)) )  &
                     &   / e3w_n(iidw(jfl,jind1),ijdw(jfl,jind2),ikdw(jfl,jind3))
               END DO
            END DO
         END DO
         
         zsumw = 0.e0
         DO jind1 = 1, 4
            DO jind2 = 1, 4
               DO jind3 = 1, 4
                  zsumw = zsumw + ztwfl(jfl,jind1,jind2,jind3) * zlagxw(jfl,jind1) * zlagyw(jfl,jind2)   &
                     &  * zlagzw(jfl,jind3) * rcoef(jind1)*rcoef(jind2)*rcoef(jind3)
               END DO
            END DO
         END DO
         pwfl(jfl) = zsumw
      END DO
      !   
      !   
   END SUBROUTINE flo_interp

   !!======================================================================
END MODULE flo4rk
