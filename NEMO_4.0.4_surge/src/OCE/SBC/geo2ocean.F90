MODULE geo2ocean
   !!======================================================================
   !!                     ***  MODULE  geo2ocean  ***
   !! Ocean mesh    :  ???
   !!======================================================================
   !! History :  OPA  !  07-1996  (O. Marti)  Original code
   !!   NEMO     1.0  !  06-2006  (G. Madec )  Free form, F90 + opt.
   !!                 !  04-2007  (S. Masson)  angle: Add T, F points and bugfix in cos lateral boundary
   !!            3.0  !  07-2008  (G. Madec)  geo2oce suppress lon/lat agruments
   !!            3.7  !  11-2015  (G. Madec)  remove the unused repere and repcmo routines
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   rot_rep       : Rotate the Repere: geographic grid <==> stretched coordinates grid
   !!   angle         :
   !!   geo2oce       :
   !!   oce2geo       :
   !!----------------------------------------------------------------------
   USE dom_oce        ! mesh and scale factors
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   rot_rep   ! called in sbccpl, fldread, and cyclone
   PUBLIC   geo2oce   ! called in sbccpl
   PUBLIC   oce2geo   ! called in sbccpl
   PUBLIC   obs_rot   ! called in obs_rot_vel and obs_write

   !                                         ! cos/sin between model grid lines and NP direction
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   gsint, gcost   ! at T point
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   gsinu, gcosu   ! at U point
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   gsinv, gcosv   ! at V point
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   gsinf, gcosf   ! at F point

   LOGICAL ,              SAVE, DIMENSION(4)     ::   linit = .FALSE.
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   gsinlon, gcoslon, gsinlat, gcoslat

   LOGICAL ::   lmust_init = .TRUE.        !: used to initialize the cos/sin variables (see above)

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE rot_rep ( pxin, pyin, cd_type, cdtodo, prot )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE rot_rep  ***
      !!
      !! ** Purpose :   Rotate the Repere: Change vector componantes between
      !!                geographic grid <--> stretched coordinates grid.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pxin, pyin   ! vector componantes
      CHARACTER(len=1),             INTENT(in   ) ::   cd_type      ! define the nature of pt2d array grid-points
      CHARACTER(len=5),             INTENT(in   ) ::   cdtodo       ! type of transpormation:
      !                                                             ! 'en->i' = east-north to i-component
      !                                                             ! 'en->j' = east-north to j-component
      !                                                             ! 'ij->e' = (i,j) components to east
      !                                                             ! 'ij->n' = (i,j) components to north
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   prot      
      !!----------------------------------------------------------------------
      !
      IF( lmust_init ) THEN      ! at 1st call only: set  gsin. & gcos.
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' rot_rep: coordinate transformation : geographic <==> model (i,j)-components'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~    '
         !
         CALL angle( glamt, gphit, glamu, gphiu, glamv, gphiv, glamf, gphif )       ! initialization of the transformation
         lmust_init = .FALSE.
      ENDIF
      !
      SELECT CASE( cdtodo )      ! type of rotation
      !
      CASE( 'en->i' )                  ! east-north to i-component
         SELECT CASE (cd_type)
         CASE ('T')   ;   prot(:,:) = pxin(:,:) * gcost(:,:) + pyin(:,:) * gsint(:,:)
         CASE ('U')   ;   prot(:,:) = pxin(:,:) * gcosu(:,:) + pyin(:,:) * gsinu(:,:)
         CASE ('V')   ;   prot(:,:) = pxin(:,:) * gcosv(:,:) + pyin(:,:) * gsinv(:,:)
         CASE ('F')   ;   prot(:,:) = pxin(:,:) * gcosf(:,:) + pyin(:,:) * gsinf(:,:)
         CASE DEFAULT   ;   CALL ctl_stop( 'Only T, U, V and F grid points are coded' )
         END SELECT
      CASE ('en->j')                   ! east-north to j-component
         SELECT CASE (cd_type)
         CASE ('T')   ;   prot(:,:) = pyin(:,:) * gcost(:,:) - pxin(:,:) * gsint(:,:)
         CASE ('U')   ;   prot(:,:) = pyin(:,:) * gcosu(:,:) - pxin(:,:) * gsinu(:,:)
         CASE ('V')   ;   prot(:,:) = pyin(:,:) * gcosv(:,:) - pxin(:,:) * gsinv(:,:)   
         CASE ('F')   ;   prot(:,:) = pyin(:,:) * gcosf(:,:) - pxin(:,:) * gsinf(:,:)   
         CASE DEFAULT   ;   CALL ctl_stop( 'Only T, U, V and F grid points are coded' )
         END SELECT
      CASE ('ij->e')                   ! (i,j)-components to east
         SELECT CASE (cd_type)
         CASE ('T')   ;   prot(:,:) = pxin(:,:) * gcost(:,:) - pyin(:,:) * gsint(:,:)
         CASE ('U')   ;   prot(:,:) = pxin(:,:) * gcosu(:,:) - pyin(:,:) * gsinu(:,:)
         CASE ('V')   ;   prot(:,:) = pxin(:,:) * gcosv(:,:) - pyin(:,:) * gsinv(:,:)
         CASE ('F')   ;   prot(:,:) = pxin(:,:) * gcosf(:,:) - pyin(:,:) * gsinf(:,:)
         CASE DEFAULT   ;   CALL ctl_stop( 'Only T, U, V and F grid points are coded' )
         END SELECT
      CASE ('ij->n')                   ! (i,j)-components to north 
         SELECT CASE (cd_type)
         CASE ('T')   ;   prot(:,:) = pyin(:,:) * gcost(:,:) + pxin(:,:) * gsint(:,:)
         CASE ('U')   ;   prot(:,:) = pyin(:,:) * gcosu(:,:) + pxin(:,:) * gsinu(:,:)
         CASE ('V')   ;   prot(:,:) = pyin(:,:) * gcosv(:,:) + pxin(:,:) * gsinv(:,:)
         CASE ('F')   ;   prot(:,:) = pyin(:,:) * gcosf(:,:) + pxin(:,:) * gsinf(:,:)
         CASE DEFAULT   ;   CALL ctl_stop( 'Only T, U, V and F grid points are coded' )
         END SELECT
      CASE DEFAULT   ;   CALL ctl_stop( 'rot_rep: Syntax Error in the definition of cdtodo' )
      !
      END SELECT
      !
   END SUBROUTINE rot_rep


   SUBROUTINE angle( plamt, pphit, plamu, pphiu, plamv, pphiv, plamf, pphif )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE angle  ***
      !! 
      !! ** Purpose :   Compute angles between model grid lines and the North direction
      !!
      !! ** Method  :   sinus and cosinus of the angle between the north-south axe 
      !!              and the j-direction at t, u, v and f-points
      !!                dot and cross products are used to obtain cos and sin, resp.
      !!
      !! ** Action  : - gsint, gcost, gsinu, gcosu, gsinv, gcosv, gsinf, gcosf
      !!----------------------------------------------------------------------
      ! WARNING: for an unexplained reason, we need to pass all glam, gphi arrays as input parameters in
      !          order to get AGRIF working with -03 compilation option
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: plamt, pphit, plamu, pphiu, plamv, pphiv, plamf, pphif  
      !
      INTEGER  ::   ji, jj   ! dummy loop indices
      INTEGER  ::   ierr     ! local integer
      REAL(wp) ::   zlam, zphi            ! local scalars
      REAL(wp) ::   zlan, zphh            !   -      -
      REAL(wp) ::   zxnpt, zynpt, znnpt   ! x,y components and norm of the vector: T point to North Pole
      REAL(wp) ::   zxnpu, zynpu, znnpu   ! x,y components and norm of the vector: U point to North Pole
      REAL(wp) ::   zxnpv, zynpv, znnpv   ! x,y components and norm of the vector: V point to North Pole
      REAL(wp) ::   zxnpf, zynpf, znnpf   ! x,y components and norm of the vector: F point to North Pole
      REAL(wp) ::   zxvvt, zyvvt, znvvt   ! x,y components and norm of the vector: between V points below and above a T point
      REAL(wp) ::   zxffu, zyffu, znffu   ! x,y components and norm of the vector: between F points below and above a U point
      REAL(wp) ::   zxffv, zyffv, znffv   ! x,y components and norm of the vector: between F points left  and right a V point
      REAL(wp) ::   zxuuf, zyuuf, znuuf   ! x,y components and norm of the vector: between U points below and above a F point
      !!----------------------------------------------------------------------
      !
      ALLOCATE( gsint(jpi,jpj), gcost(jpi,jpj),   & 
         &      gsinu(jpi,jpj), gcosu(jpi,jpj),   & 
         &      gsinv(jpi,jpj), gcosv(jpi,jpj),   &  
         &      gsinf(jpi,jpj), gcosf(jpi,jpj), STAT=ierr )
      CALL mpp_sum( 'geo2ocean', ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'angle: unable to allocate arrays' )
      !
      ! ============================= !
      ! Compute the cosinus and sinus !
      ! ============================= !
      ! (computation done on the north stereographic polar plane)
      !
      DO jj = 2, jpjm1
         DO ji = fs_2, jpi   ! vector opt.
            !                  
            zlam = plamt(ji,jj)     ! north pole direction & modulous (at t-point)
            zphi = pphit(ji,jj)
            zxnpt = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            zynpt = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            znnpt = zxnpt*zxnpt + zynpt*zynpt
            !
            zlam = plamu(ji,jj)     ! north pole direction & modulous (at u-point)
            zphi = pphiu(ji,jj)
            zxnpu = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            zynpu = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            znnpu = zxnpu*zxnpu + zynpu*zynpu
            !
            zlam = plamv(ji,jj)     ! north pole direction & modulous (at v-point)
            zphi = pphiv(ji,jj)
            zxnpv = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            zynpv = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            znnpv = zxnpv*zxnpv + zynpv*zynpv
            !
            zlam = plamf(ji,jj)     ! north pole direction & modulous (at f-point)
            zphi = pphif(ji,jj)
            zxnpf = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            zynpf = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            znnpf = zxnpf*zxnpf + zynpf*zynpf
            !
            zlam = plamv(ji,jj  )   ! j-direction: v-point segment direction (around t-point)
            zphi = pphiv(ji,jj  )
            zlan = plamv(ji,jj-1)
            zphh = pphiv(ji,jj-1)
            zxvvt =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            zyvvt =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            znvvt = SQRT( znnpt * ( zxvvt*zxvvt + zyvvt*zyvvt )  )
            znvvt = MAX( znvvt, 1.e-14 )
            !
            zlam = plamf(ji,jj  )   ! j-direction: f-point segment direction (around u-point)
            zphi = pphif(ji,jj  )
            zlan = plamf(ji,jj-1)
            zphh = pphif(ji,jj-1)
            zxffu =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            zyffu =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            znffu = SQRT( znnpu * ( zxffu*zxffu + zyffu*zyffu )  )
            znffu = MAX( znffu, 1.e-14 )
            !
            zlam = plamf(ji  ,jj)   ! i-direction: f-point segment direction (around v-point)
            zphi = pphif(ji  ,jj)
            zlan = plamf(ji-1,jj)
            zphh = pphif(ji-1,jj)
            zxffv =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            zyffv =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            znffv = SQRT( znnpv * ( zxffv*zxffv + zyffv*zyffv )  )
            znffv = MAX( znffv, 1.e-14 )
            !
            zlam = plamu(ji,jj+1)   ! j-direction: u-point segment direction (around f-point)
            zphi = pphiu(ji,jj+1)
            zlan = plamu(ji,jj  )
            zphh = pphiu(ji,jj  )
            zxuuf =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            zyuuf =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            znuuf = SQRT( znnpf * ( zxuuf*zxuuf + zyuuf*zyuuf )  )
            znuuf = MAX( znuuf, 1.e-14 )
            !
            !                       ! cosinus and sinus using dot and cross products
            gsint(ji,jj) = ( zxnpt*zyvvt - zynpt*zxvvt ) / znvvt
            gcost(ji,jj) = ( zxnpt*zxvvt + zynpt*zyvvt ) / znvvt
            !
            gsinu(ji,jj) = ( zxnpu*zyffu - zynpu*zxffu ) / znffu
            gcosu(ji,jj) = ( zxnpu*zxffu + zynpu*zyffu ) / znffu
            !
            gsinf(ji,jj) = ( zxnpf*zyuuf - zynpf*zxuuf ) / znuuf
            gcosf(ji,jj) = ( zxnpf*zxuuf + zynpf*zyuuf ) / znuuf
            !
            gsinv(ji,jj) = ( zxnpv*zxffv + zynpv*zyffv ) / znffv
            gcosv(ji,jj) =-( zxnpv*zyffv - zynpv*zxffv ) / znffv     ! (caution, rotation of 90 degres)
            !
         END DO
      END DO

      ! =============== !
      ! Geographic mesh !
      ! =============== !

      DO jj = 2, jpjm1
         DO ji = fs_2, jpi   ! vector opt.
            IF( MOD( ABS( plamv(ji,jj) - plamv(ji,jj-1) ), 360. ) < 1.e-8 ) THEN
               gsint(ji,jj) = 0.
               gcost(ji,jj) = 1.
            ENDIF
            IF( MOD( ABS( plamf(ji,jj) - plamf(ji,jj-1) ), 360. ) < 1.e-8 ) THEN
               gsinu(ji,jj) = 0.
               gcosu(ji,jj) = 1.
            ENDIF
            IF(      ABS( pphif(ji,jj) - pphif(ji-1,jj) )         < 1.e-8 ) THEN
               gsinv(ji,jj) = 0.
               gcosv(ji,jj) = 1.
            ENDIF
            IF( MOD( ABS( plamu(ji,jj) - plamu(ji,jj+1) ), 360. ) < 1.e-8 ) THEN
               gsinf(ji,jj) = 0.
               gcosf(ji,jj) = 1.
            ENDIF
         END DO
      END DO

      ! =========================== !
      ! Lateral boundary conditions !
      ! =========================== !
      !           ! lateral boundary cond.: T-, U-, V-, F-pts, sgn
      CALL lbc_lnk_multi( 'geo2ocean', gcost, 'T', -1., gsint, 'T', -1., gcosu, 'U', -1., gsinu, 'U', -1., & 
                      &   gcosv, 'V', -1., gsinv, 'V', -1., gcosf, 'F', -1., gsinf, 'F', -1.  )
      !
   END SUBROUTINE angle


   SUBROUTINE geo2oce ( pxx, pyy, pzz, cgrid, pte, ptn )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE geo2oce  ***
      !!      
      !! ** Purpose :
      !!
      !! ** Method  :   Change a vector from geocentric to east/north 
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::  pxx, pyy, pzz
      CHARACTER(len=1)            , INTENT(in   ) ::  cgrid
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::  pte, ptn
      !
      REAL(wp), PARAMETER :: rpi = 3.141592653e0
      REAL(wp), PARAMETER :: rad = rpi / 180.e0
      INTEGER ::   ig     !
      INTEGER ::   ierr   ! local integer
      !!----------------------------------------------------------------------
      !
      IF( .NOT. ALLOCATED( gsinlon ) ) THEN
         ALLOCATE( gsinlon(jpi,jpj,4) , gcoslon(jpi,jpj,4) ,   &
            &      gsinlat(jpi,jpj,4) , gcoslat(jpi,jpj,4) , STAT=ierr )
         CALL mpp_sum( 'geo2ocean', ierr )
         IF( ierr /= 0 )   CALL ctl_stop('geo2oce: unable to allocate arrays' )
      ENDIF
      !
      SELECT CASE( cgrid)
      CASE ( 'T' )   
         ig = 1
         IF( .NOT. linit(ig) ) THEN 
            gsinlon(:,:,ig) = SIN( rad * glamt(:,:) )
            gcoslon(:,:,ig) = COS( rad * glamt(:,:) )
            gsinlat(:,:,ig) = SIN( rad * gphit(:,:) )
            gcoslat(:,:,ig) = COS( rad * gphit(:,:) )
            linit(ig) = .TRUE.
         ENDIF
      CASE ( 'U' )   
         ig = 2
         IF( .NOT. linit(ig) ) THEN 
            gsinlon(:,:,ig) = SIN( rad * glamu(:,:) )
            gcoslon(:,:,ig) = COS( rad * glamu(:,:) )
            gsinlat(:,:,ig) = SIN( rad * gphiu(:,:) )
            gcoslat(:,:,ig) = COS( rad * gphiu(:,:) )
            linit(ig) = .TRUE.
         ENDIF
      CASE ( 'V' )   
         ig = 3
         IF( .NOT. linit(ig) ) THEN 
            gsinlon(:,:,ig) = SIN( rad * glamv(:,:) )
            gcoslon(:,:,ig) = COS( rad * glamv(:,:) )
            gsinlat(:,:,ig) = SIN( rad * gphiv(:,:) )
            gcoslat(:,:,ig) = COS( rad * gphiv(:,:) )
            linit(ig) = .TRUE.
         ENDIF
      CASE ( 'F' )   
         ig = 4
         IF( .NOT. linit(ig) ) THEN 
            gsinlon(:,:,ig) = SIN( rad * glamf(:,:) )
            gcoslon(:,:,ig) = COS( rad * glamf(:,:) )
            gsinlat(:,:,ig) = SIN( rad * gphif(:,:) )
            gcoslat(:,:,ig) = COS( rad * gphif(:,:) )
            linit(ig) = .TRUE.
         ENDIF
      CASE default   
         WRITE(ctmp1,*) 'geo2oce : bad grid argument : ', cgrid
         CALL ctl_stop( ctmp1 )
      END SELECT
      !
      pte = - gsinlon(:,:,ig) * pxx + gcoslon(:,:,ig) * pyy
      ptn = - gcoslon(:,:,ig) * gsinlat(:,:,ig) * pxx    &
         &  - gsinlon(:,:,ig) * gsinlat(:,:,ig) * pyy    &
         &  + gcoslat(:,:,ig) * pzz
      !
   END SUBROUTINE geo2oce


   SUBROUTINE oce2geo ( pte, ptn, cgrid, pxx , pyy , pzz )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE oce2geo  ***
      !!      
      !! ** Purpose :
      !!
      !! ** Method  :   Change vector from east/north to geocentric
      !!
      !! History :     ! (A. Caubel)  oce2geo - Original code
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT( IN    ) ::  pte, ptn
      CHARACTER(len=1)            , INTENT( IN    ) ::  cgrid
      REAL(wp), DIMENSION(jpi,jpj), INTENT(   OUT ) ::  pxx , pyy , pzz
      !!
      REAL(wp), PARAMETER :: rpi = 3.141592653E0
      REAL(wp), PARAMETER :: rad = rpi / 180.e0
      INTEGER ::   ig     !
      INTEGER ::   ierr   ! local integer
      !!----------------------------------------------------------------------

      IF( .NOT. ALLOCATED( gsinlon ) ) THEN
         ALLOCATE( gsinlon(jpi,jpj,4) , gcoslon(jpi,jpj,4) ,   &
            &      gsinlat(jpi,jpj,4) , gcoslat(jpi,jpj,4) , STAT=ierr )
         CALL mpp_sum( 'geo2ocean', ierr )
         IF( ierr /= 0 )   CALL ctl_stop('oce2geo: unable to allocate arrays' )
      ENDIF

      SELECT CASE( cgrid)
         CASE ( 'T' )   
            ig = 1
            IF( .NOT. linit(ig) ) THEN 
               gsinlon(:,:,ig) = SIN( rad * glamt(:,:) )
               gcoslon(:,:,ig) = COS( rad * glamt(:,:) )
               gsinlat(:,:,ig) = SIN( rad * gphit(:,:) )
               gcoslat(:,:,ig) = COS( rad * gphit(:,:) )
               linit(ig) = .TRUE.
            ENDIF
         CASE ( 'U' )   
            ig = 2
            IF( .NOT. linit(ig) ) THEN 
               gsinlon(:,:,ig) = SIN( rad * glamu(:,:) )
               gcoslon(:,:,ig) = COS( rad * glamu(:,:) )
               gsinlat(:,:,ig) = SIN( rad * gphiu(:,:) )
               gcoslat(:,:,ig) = COS( rad * gphiu(:,:) )
               linit(ig) = .TRUE.
            ENDIF
         CASE ( 'V' )   
            ig = 3
            IF( .NOT. linit(ig) ) THEN 
               gsinlon(:,:,ig) = SIN( rad * glamv(:,:) )
               gcoslon(:,:,ig) = COS( rad * glamv(:,:) )
               gsinlat(:,:,ig) = SIN( rad * gphiv(:,:) )
               gcoslat(:,:,ig) = COS( rad * gphiv(:,:) )
               linit(ig) = .TRUE.
            ENDIF
         CASE ( 'F' )   
            ig = 4
            IF( .NOT. linit(ig) ) THEN 
               gsinlon(:,:,ig) = SIN( rad * glamf(:,:) )
               gcoslon(:,:,ig) = COS( rad * glamf(:,:) )
               gsinlat(:,:,ig) = SIN( rad * gphif(:,:) )
               gcoslat(:,:,ig) = COS( rad * gphif(:,:) )
               linit(ig) = .TRUE.
            ENDIF
         CASE default   
            WRITE(ctmp1,*) 'geo2oce : bad grid argument : ', cgrid
            CALL ctl_stop( ctmp1 )
      END SELECT
      !
      pxx = - gsinlon(:,:,ig) * pte - gcoslon(:,:,ig) * gsinlat(:,:,ig) * ptn 
      pyy =   gcoslon(:,:,ig) * pte - gsinlon(:,:,ig) * gsinlat(:,:,ig) * ptn
      pzz =   gcoslat(:,:,ig) * ptn
      !
   END SUBROUTINE oce2geo


   SUBROUTINE obs_rot( psinu, pcosu, psinv, pcosv )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE obs_rot  ***
      !!
      !! ** Purpose :   Copy gsinu, gcosu, gsinv and gsinv
      !!                to input data for rotations of
      !!                current at observation points
      !!
      !! History :  9.2  !  09-02  (K. Mogensen)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT( OUT )::   psinu, pcosu, psinv, pcosv   ! copy of data
      !!----------------------------------------------------------------------
      !
      ! Initialization of gsin* and gcos* at first call
      ! -----------------------------------------------
      IF( lmust_init ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' obs_rot : geographic <--> stretched'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~   coordinate transformation'
         CALL angle( glamt, gphit, glamu, gphiu, glamv, gphiv, glamf, gphif )       ! initialization of the transformation
         lmust_init = .FALSE.
      ENDIF
      !
      psinu(:,:) = gsinu(:,:)
      pcosu(:,:) = gcosu(:,:)
      psinv(:,:) = gsinv(:,:)
      pcosv(:,:) = gcosv(:,:)
      !
   END SUBROUTINE obs_rot

  !!======================================================================
END MODULE geo2ocean
