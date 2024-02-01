MODULE ldfc1d_c2d
   !!======================================================================
   !!                    ***  MODULE  ldfc1d_c2d  ***
   !! Ocean physics:  profile and horizontal shape of lateral eddy coefficients 
   !!=====================================================================
   !! History :  3.7  ! 2013-12  (G. Madec)  restructuration/simplification of aht/aeiv specification,
   !!                 !                      add velocity dependent coefficient and optional read in file
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_c1d       : ah reduced by 1/4 on the vertical (tanh profile, inflection at 300m) 
   !!   ldf_c2d       : ah = F(e1,e2) (laplacian or = F(e1^3,e2^3) (bilaplacian)
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_c1d   ! called by ldftra and ldfdyn modules
   PUBLIC   ldf_c2d   ! called by ldftra and ldfdyn modules

   REAL(wp) ::   r1_2  = 0.5_wp           ! =1/2
   REAL(wp) ::   r1_4  = 0.25_wp          ! =1/4
   REAL(wp) ::   r1_12 = 1._wp / 12._wp   ! =1/12
 
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_c1d( cd_type, pahs1, pahs2, pah1, pah2 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_c1d  ***
      !!              
      !! ** Purpose :   1D eddy diffusivity/viscosity coefficients
      !!
      !! ** Method  :   1D eddy diffusivity coefficients F( depth )
      !!                Reduction by zratio from surface to bottom 
      !!                hyperbolic tangent profile with inflection point 
      !!                at zh=500m and a width of zw=200m
      !!
      !!   cd_type = TRA      pah1, pah2 defined at U- and V-points
      !!             DYN      pah1, pah2 defined at T- and F-points
      !!----------------------------------------------------------------------
      CHARACTER(len=3)                , INTENT(in   ) ::   cd_type        ! DYNamique or TRAcers
      REAL(wp), DIMENSION(jpi,jpj)    , INTENT(in   ) ::   pahs1, pahs2   ! surface value of eddy coefficient   [m2/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pah1 , pah2    ! eddy coefficient                    [m2/s]
      !
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   zh, zc, zdep1   ! local scalars
      REAL(wp) ::   zw    , zdep2   !   -      -
      REAL(wp) ::   zratio          !   -      -
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '   ldf_c1d : set a given profile to eddy mixing coefficients'
      !
      ! initialization of the profile
      zratio = 0.25_wp           ! surface/bottom ratio
      zh =  500._wp              ! depth    of the inflection point [m]
      zw =  1._wp / 200._wp      ! width^-1     -        -      -   [1/m]
      !                          ! associated coefficient           [-]
      zc = ( 1._wp - zratio ) / ( 1._wp + TANH( zh * zw) )
      !
      !
      SELECT CASE( cd_type )        ! point of calculation
      !
      CASE( 'DYN' )                     ! T- and F-points
         DO jk = jpkm1, 1, -1                ! pah1 at T-point
            pah1(:,:,jk) = pahs1(:,:) * (  zratio + zc * ( 1._wp + TANH( - ( gdept_0(:,:,jk) - zh ) * zw) )  )
         END DO
         DO jk = jpkm1, 1, -1                ! pah2 at F-point (zdep2 is an approximation in zps-coord.)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zdep2 = (  gdept_0(ji,jj+1,jk) + gdept_0(ji+1,jj+1,jk)   &
                     &     + gdept_0(ji,jj  ,jk) + gdept_0(ji+1,jj  ,jk)  ) * r1_4
                  pah2(ji,jj,jk) = pahs2(ji,jj) * (  zratio + zc * ( 1._wp + TANH( - ( zdep2 - zh ) * zw) )  )
               END DO
            END DO
         END DO
         CALL lbc_lnk( 'ldfc1d_c2d', pah2, 'F', 1. )   ! Lateral boundary conditions
         !
      CASE( 'TRA' )                     ! U- and V-points (zdep1 & 2 are an approximation in zps-coord.)
         DO jk = jpkm1, 1, -1
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zdep1 = (  gdept_0(ji,jj,jk) + gdept_0(ji+1,jj,jk)  ) * 0.5_wp
                  zdep2 = (  gdept_0(ji,jj,jk) + gdept_0(ji,jj+1,jk)  ) * 0.5_wp
                  pah1(ji,jj,jk) = pahs1(ji,jj) * (  zratio + zc * ( 1._wp + TANH( - ( zdep1 - zh ) * zw) )  )
                  pah2(ji,jj,jk) = pahs2(ji,jj) * (  zratio + zc * ( 1._wp + TANH( - ( zdep2 - zh ) * zw) )  )
               END DO
            END DO
         END DO
         ! Lateral boundary conditions
         CALL lbc_lnk_multi( 'ldfc1d_c2d', pah1, 'U', 1. , pah2, 'V', 1. )   
         !
      CASE DEFAULT                        ! error
         CALL ctl_stop( 'ldf_c1d: ', cd_type, ' Unknown, i.e. /= DYN or TRA' )
      END SELECT
      !
   END SUBROUTINE ldf_c1d


   SUBROUTINE ldf_c2d( cd_type, pUfac, knn, pah1, pah2 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_c2d  ***
      !!              
      !! ** Purpose :   2D eddy diffusivity/viscosity coefficients
      !!
      !! ** Method  :   2D eddy diffusivity coefficients F( e1 , e2 )
      !!       laplacian   operator :   ah proportional to the scale factor      [m2/s]
      !!       bilaplacian operator :   ah proportional to the (scale factor)^3  [m4/s]
      !!       In both cases, pah0 is the maximum value reached by the coefficient 
      !!       at the Equator in case of e1=ra*rad= ~111km, not over the whole domain.
      !!
      !!   cd_type = TRA      pah1, pah2 defined at U- and V-points
      !!             DYN      pah1, pah2 defined at T- and F-points
      !!----------------------------------------------------------------------
      CHARACTER(len=3)          , INTENT(in   ) ::   cd_type      ! DYNamique or TRAcers
      REAL(wp)                  , INTENT(in   ) ::   pUfac        ! =1/2*Uc LAPlacian BiLaPlacian
      INTEGER                   , INTENT(in   ) ::   knn          ! characteristic velocity   [m/s]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pah1, pah2   ! eddy coefficients         [m2/s or m4/s]
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   inn          ! local integer
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '   ldf_c2d :   aht = Ufac * max(e1,e2)   with Ufac = ', pUfac, ' m/s'
      !
      !
      SELECT CASE( cd_type )        !==  surface values  ==!  (chosen grid point function of DYN or TRA)
      !
      CASE( 'DYN' )                       ! T- and F-points
         DO jj = 1, jpj
            DO ji = 1, jpi 
               pah1(ji,jj,1) = pUfac * MAX( e1t(ji,jj) , e2t(ji,jj) )**knn
               pah2(ji,jj,1) = pUfac * MAX( e1f(ji,jj) , e2f(ji,jj) )**knn
            END DO
         END DO
      CASE( 'TRA' )                       ! U- and V-points
         DO jj = 1, jpj 
            DO ji = 1, jpi 
               pah1(ji,jj,1) = pUfac * MAX( e1u(ji,jj), e2u(ji,jj) )**knn
               pah2(ji,jj,1) = pUfac * MAX( e1v(ji,jj), e2v(ji,jj) )**knn
            END DO
         END DO
      CASE DEFAULT                        ! error
         CALL ctl_stop( 'ldf_c2d: ', cd_type, ' Unknown, i.e. /= DYN or TRA' )
      END SELECT
      !                             !==  deeper values = surface one  ==!  (except jpk)
      DO jk = 2, jpkm1
         pah1(:,:,jk) = pah1(:,:,1)
         pah2(:,:,jk) = pah2(:,:,1)
      END DO
      !
   END SUBROUTINE ldf_c2d

   !!======================================================================
END MODULE ldfc1d_c2d
