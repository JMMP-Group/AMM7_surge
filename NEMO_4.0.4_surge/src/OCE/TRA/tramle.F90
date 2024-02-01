MODULE tramle
   !!======================================================================
   !!                    ***  MODULE  tramle  ***
   !! Ocean tracers: Mixed Layer Eddy induced transport
   !!======================================================================
   !! History :  3.3  !  2010-08  (G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_mle_trp   : update the effective transport with the Mixed Layer Eddy induced transport
   !!   tra_mle_init  : initialisation of the Mixed Layer Eddy induced transport computation
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE phycst         ! physical constant
   USE zdfmxl         ! mixed layer depth
   !
   USE in_out_manager ! I/O manager
   USE iom            ! IOM library
   USE lib_mpp        ! MPP library
   USE lbclnk         ! lateral boundary condition / mpp link

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_mle_trp        ! routine called in traadv.F90
   PUBLIC   tra_mle_init   ! routine called in traadv.F90

   !                                    !!* namelist namtra_mle *
   LOGICAL, PUBLIC ::   ln_mle           !: flag to activate the Mixed Layer Eddy (MLE) parameterisation
   INTEGER         ::      nn_mle           ! MLE type: =0 standard Fox-Kemper ; =1 new formulation
   INTEGER         ::      nn_mld_uv        ! space interpolation of MLD at u- & v-pts (0=min,1=averaged,2=max)
   INTEGER         ::      nn_conv          ! =1 no MLE in case of convection ; =0 always MLE
   REAL(wp)        ::      rn_ce            ! MLE coefficient
   !                                        ! parameters used in nn_mle = 0 case
   REAL(wp)        ::      rn_lf               ! typical scale of mixed layer front
   REAL(wp)        ::      rn_time             ! time scale for mixing momentum across the mixed layer
   !                                        ! parameters used in nn_mle = 1 case
   REAL(wp)        ::      rn_lat              ! reference latitude for a 5 km scale of ML front
   REAL(wp)        ::      rn_rho_c_mle        ! Density criterion for definition of MLD used by FK

   REAL(wp) ::   r5_21 = 5.e0 / 21.e0   ! factor used in mle streamfunction computation
   REAL(wp) ::   rb_c                   ! ML buoyancy criteria = g rho_c /rau0 where rho_c is defined in zdfmld
   REAL(wp) ::   rc_f                   ! MLE coefficient (= rn_ce / (5 km * fo) ) in nn_mle=1 case

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   rfu, rfv   ! modified Coriolis parameter (f+tau) at u- & v-pts
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   r1_ft      ! inverse of the modified Coriolis parameter at t-pts

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_mle_trp( kt, kit000, pu, pv, pw, cdtype )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_mle_trp  ***
      !!
      !! ** Purpose :   Add to the transport the Mixed Layer Eddy induced transport
      !!
      !! ** Method  :   The 3 components of the Mixed Layer Eddy (MLE) induced
      !!              transport are computed as follows :
      !!                zu_mle = dk[ zpsi_uw ]
      !!                zv_mle = dk[ zpsi_vw ]
      !!                zw_mle = - di[ zpsi_uw ] - dj[ zpsi_vw ]
      !!                where zpsi is the MLE streamfunction at uw and vw points (see the doc)
      !!              and added to the input velocity :
      !!                p.n = p.n + z._mle
      !!
      !! ** Action  : - (pun,pvn,pwn) increased by the mle transport
      !!                CAUTION, the transport is not updated at the last line/raw
      !!                         this may be a problem for some advection schemes
      !!
      !! References: Fox-Kemper et al., JPO, 38, 1145-1165, 2008
      !!             Fox-Kemper and Ferrari, JPO, 38, 1166-1179, 2008
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   kit000     ! first time step index
      CHARACTER(len=3)                , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu         ! in : 3 ocean transport components
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pv         ! out: same 3  transport components
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pw         !   increased by the MLE induced transport
      !
      INTEGER  ::   ji, jj, jk          ! dummy loop indices
      INTEGER  ::   ii, ij, ik, ikmax   ! local integers
      REAL(wp) ::   zcuw, zmuw, zc      ! local scalar
      REAL(wp) ::   zcvw, zmvw          !   -      -
      INTEGER , DIMENSION(jpi,jpj)     :: inml_mle
      REAL(wp), DIMENSION(jpi,jpj)     :: zpsim_u, zpsim_v, zmld, zbm, zhu, zhv, zn2, zLf_NH, zLf_MH
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpsi_uw, zpsi_vw
      !!----------------------------------------------------------------------
      !
      !                                      !==  MLD used for MLE  ==!
      !                                                ! compute from the 10m density to deal with the diurnal cycle
      inml_mle(:,:) = mbkt(:,:) + 1                    ! init. to number of ocean w-level (T-level + 1)
      IF ( nla10 > 0 ) THEN                            ! avoid case where first level is thicker than 10m
         DO jk = jpkm1, nlb10, -1                      ! from the bottom to nlb10 (10m)
            DO jj = 1, jpj
               DO ji = 1, jpi                          ! index of the w-level at the ML based
                  IF( rhop(ji,jj,jk) > rhop(ji,jj,nla10) + rn_rho_c_mle )   inml_mle(ji,jj) = jk      ! Mixed layer
               END DO
            END DO
         END DO
      ENDIF
      ikmax = MIN( MAXVAL( inml_mle(:,:) ), jpkm1 )                  ! max level of the computation
      !
      !
      zmld(:,:) = 0._wp                      !==   Horizontal shape of the MLE  ==!
      zbm (:,:) = 0._wp
      zn2 (:,:) = 0._wp
      DO jk = 1, ikmax                                 ! MLD and mean buoyancy and N2 over the mixed layer
         DO jj = 1, jpj
            DO ji = 1, jpi
               zc = e3t_n(ji,jj,jk) * REAL( MIN( MAX( 0, inml_mle(ji,jj)-jk ) , 1  )  )    ! zc being 0 outside the ML t-points
               zmld(ji,jj) = zmld(ji,jj) + zc
               zbm (ji,jj) = zbm (ji,jj) + zc * (rau0 - rhop(ji,jj,jk) ) * r1_rau0
               zn2 (ji,jj) = zn2 (ji,jj) + zc * (rn2(ji,jj,jk)+rn2(ji,jj,jk+1))*0.5_wp
            END DO
         END DO
      END DO

      SELECT CASE( nn_mld_uv )                         ! MLD at u- & v-pts
      CASE ( 0 )                                               != min of the 2 neighbour MLDs
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zhu(ji,jj) = MIN( zmld(ji+1,jj), zmld(ji,jj) )
               zhv(ji,jj) = MIN( zmld(ji,jj+1), zmld(ji,jj) )
            END DO
         END DO
      CASE ( 1 )                                               != average of the 2 neighbour MLDs
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zhu(ji,jj) = ( zmld(ji+1,jj) + zmld(ji,jj) ) * 0.5_wp
               zhv(ji,jj) = ( zmld(ji,jj+1) + zmld(ji,jj) ) * 0.5_wp
            END DO
         END DO
      CASE ( 2 )                                               != max of the 2 neighbour MLDs
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zhu(ji,jj) = MAX( zmld(ji+1,jj), zmld(ji,jj) )
               zhv(ji,jj) = MAX( zmld(ji,jj+1), zmld(ji,jj) )
            END DO
         END DO
      END SELECT
      !                                                ! convert density into buoyancy
      zbm(:,:) = + grav * zbm(:,:) / MAX( e3t_n(:,:,1), zmld(:,:) )
      !
      !
      !                                      !==  Magnitude of the MLE stream function  ==!
      !
      !                 di[bm]  Ds
      ! Psi = Ce  H^2 ---------------- e2u  mu(z)   where fu Lf = MAX( fu*rn_fl , (Db H)^1/2 )
      !                  e1u   Lf fu                      and the e2u for the "transport"
      !                                                      (not *e3u as divided by e3u at the end)
      !
      IF( nn_mle == 0 ) THEN           ! Fox-Kemper et al. 2010 formulation
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zpsim_u(ji,jj) = rn_ce * zhu(ji,jj) * zhu(ji,jj)  * e2_e1u(ji,jj)                                            &
                  &           * ( zbm(ji+1,jj) - zbm(ji,jj) ) * MIN( 111.e3_wp , e1u(ji,jj) )   &
                  &           / (  MAX( rn_lf * rfu(ji,jj) , SQRT( rb_c * zhu(ji,jj) ) )   )
                  !
               zpsim_v(ji,jj) = rn_ce * zhv(ji,jj) * zhv(ji,jj)  * e1_e2v(ji,jj)                                            &
                  &           * ( zbm(ji,jj+1) - zbm(ji,jj) ) * MIN( 111.e3_wp , e2v(ji,jj) )   &
                  &           / (  MAX( rn_lf * rfv(ji,jj) , SQRT( rb_c * zhv(ji,jj) ) )   )
            END DO
         END DO
         !
      ELSEIF( nn_mle == 1 ) THEN       ! New formulation (Lf = 5km fo/ff with fo=Coriolis parameter at latitude rn_lat)
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zpsim_u(ji,jj) = rc_f *   zhu(ji,jj)   * zhu(ji,jj)   * e2_e1u(ji,jj)               &
                  &                  * ( zbm(ji+1,jj) - zbm(ji,jj) ) * MIN( 111.e3_wp , e1u(ji,jj) )
                  !
               zpsim_v(ji,jj) = rc_f *   zhv(ji,jj)   * zhv(ji,jj)   * e1_e2v(ji,jj)               &
                  &                  * ( zbm(ji,jj+1) - zbm(ji,jj) ) * MIN( 111.e3_wp , e2v(ji,jj) )
            END DO
         END DO
      ENDIF
      !
      IF( nn_conv == 1 ) THEN              ! No MLE in case of convection
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               IF( MIN( zn2(ji,jj) , zn2(ji+1,jj) ) < 0._wp )   zpsim_u(ji,jj) = 0._wp
               IF( MIN( zn2(ji,jj) , zn2(ji,jj+1) ) < 0._wp )   zpsim_v(ji,jj) = 0._wp
            END DO
         END DO
      ENDIF
      !
      !                                      !==  structure function value at uw- and vw-points  ==!
      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1   ! vector opt.
            zhu(ji,jj) = 1._wp / zhu(ji,jj)                   ! hu --> 1/hu
            zhv(ji,jj) = 1._wp / zhv(ji,jj)
         END DO
      END DO
      !
      zpsi_uw(:,:,:) = 0._wp
      zpsi_vw(:,:,:) = 0._wp
      !
      DO jk = 2, ikmax                                ! start from 2 : surface value = 0
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zcuw = 1._wp - ( gdepw_n(ji+1,jj,jk) + gdepw_n(ji,jj,jk) ) * zhu(ji,jj)
               zcvw = 1._wp - ( gdepw_n(ji,jj+1,jk) + gdepw_n(ji,jj,jk) ) * zhv(ji,jj)
               zcuw = zcuw * zcuw
               zcvw = zcvw * zcvw
               zmuw = MAX(  0._wp , ( 1._wp - zcuw ) * ( 1._wp + r5_21 * zcuw )  )
               zmvw = MAX(  0._wp , ( 1._wp - zcvw ) * ( 1._wp + r5_21 * zcvw )  )
               !
               zpsi_uw(ji,jj,jk) = zpsim_u(ji,jj) * zmuw * umask(ji,jj,jk)
               zpsi_vw(ji,jj,jk) = zpsim_v(ji,jj) * zmvw * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      !                                      !==  transport increased by the MLE induced transport ==!
      DO jk = 1, ikmax
         DO jj = 1, jpjm1                          ! CAUTION pu,pv must be defined at row/column i=1 / j=1
            DO ji = 1, fs_jpim1   ! vector opt.
               pu(ji,jj,jk) = pu(ji,jj,jk) + ( zpsi_uw(ji,jj,jk) - zpsi_uw(ji,jj,jk+1) )
               pv(ji,jj,jk) = pv(ji,jj,jk) + ( zpsi_vw(ji,jj,jk) - zpsi_vw(ji,jj,jk+1) )
            END DO
         END DO
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               pw(ji,jj,jk) = pw(ji,jj,jk) - ( zpsi_uw(ji,jj,jk) - zpsi_uw(ji-1,jj,jk)   &
                  &                          + zpsi_vw(ji,jj,jk) - zpsi_vw(ji,jj-1,jk) )
            END DO
         END DO
      END DO

      IF( cdtype == 'TRA') THEN              !==  outputs  ==!
         !
         zLf_NH(:,:) = SQRT( rb_c * zmld(:,:) ) * r1_ft(:,:)      ! Lf = N H / f
         CALL iom_put( "Lf_NHpf" , zLf_NH  )    ! Lf = N H / f
         !
         ! divide by cross distance to give streamfunction with dimensions m^2/s
         DO jk = 1, ikmax+1
            zpsi_uw(:,:,jk) = zpsi_uw(:,:,jk) * r1_e2u(:,:)
            zpsi_vw(:,:,jk) = zpsi_vw(:,:,jk) * r1_e1v(:,:)
         END DO
         CALL iom_put( "psiu_mle", zpsi_uw )    ! i-mle streamfunction
         CALL iom_put( "psiv_mle", zpsi_vw )    ! j-mle streamfunction
      ENDIF
      !
   END SUBROUTINE tra_mle_trp


   SUBROUTINE tra_mle_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE tra_mle_init  ***
      !!
      !! ** Purpose :   Control the consistency between namelist options for
      !!              tracer advection schemes and set nadv
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ierr
      INTEGER ::    ios                 ! Local integer output status for namelist read
      REAL(wp) ::   z1_t2, zfu, zfv                                !    -         -
      !
      NAMELIST/namtra_mle/ ln_mle , nn_mle, rn_ce, rn_lf, rn_time, rn_lat, nn_mld_uv, nn_conv, rn_rho_c_mle
      !!----------------------------------------------------------------------

      REWIND( numnam_ref )              ! Namelist namtra_mle in reference namelist : Tracer advection scheme
      READ  ( numnam_ref, namtra_mle, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namtra_mle in reference namelist' )

      REWIND( numnam_cfg )              ! Namelist namtra_mle in configuration namelist : Tracer advection scheme
      READ  ( numnam_cfg, namtra_mle, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namtra_mle in configuration namelist' )
      IF(lwm) WRITE ( numond, namtra_mle )

      IF(lwp) THEN                     ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_mle_init : mixed layer eddy (MLE) advection acting on tracers'
         WRITE(numout,*) '~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_mle : mixed layer eddy advection applied on tracers'
         WRITE(numout,*) '      use mixed layer eddy (MLE, i.e. Fox-Kemper param) (T/F)      ln_mle       = ', ln_mle
         WRITE(numout,*) '         MLE type: =0 standard Fox-Kemper ; =1 new formulation        nn_mle    = ', nn_mle
         WRITE(numout,*) '         magnitude of the MLE (typical value: 0.06 to 0.08)           rn_ce     = ', rn_ce
         WRITE(numout,*) '         scale of ML front (ML radius of deformation) (rn_mle=0)      rn_lf     = ', rn_lf, 'm'
         WRITE(numout,*) '         maximum time scale of MLE                    (rn_mle=0)      rn_time   = ', rn_time, 's'
         WRITE(numout,*) '         reference latitude (degrees) of MLE coef.    (rn_mle=1)      rn_lat    = ', rn_lat, 'deg'
         WRITE(numout,*) '         space interp. of MLD at u-(v-)pts (0=min,1=averaged,2=max)   nn_mld_uv = ', nn_mld_uv
         WRITE(numout,*) '         =1 no MLE in case of convection ; =0 always MLE              nn_conv   = ', nn_conv
         WRITE(numout,*) '         Density difference used to define ML for FK              rn_rho_c_mle  = ', rn_rho_c_mle
      ENDIF
      !
      IF(lwp) THEN
         WRITE(numout,*)
         IF( ln_mle ) THEN
            WRITE(numout,*) '   ==>>>   Mixed Layer Eddy induced transport added to tracer advection'
            IF( nn_mle == 0 )   WRITE(numout,*) '              Fox-Kemper et al 2010 formulation'
            IF( nn_mle == 1 )   WRITE(numout,*) '              New formulation'
         ELSE
            WRITE(numout,*) '   ==>>>   Mixed Layer Eddy parametrisation NOT used'
         ENDIF
      ENDIF
      !
      IF( ln_mle ) THEN                ! MLE initialisation
         !
         rb_c = grav * rn_rho_c_mle /rau0        ! Mixed Layer buoyancy criteria
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '      ML buoyancy criteria = ', rb_c, ' m/s2 '
         IF(lwp) WRITE(numout,*) '      associated ML density criteria defined in zdfmxl = ', rho_c, 'kg/m3'
         !
         IF( nn_mle == 0 ) THEN           ! MLE array allocation & initialisation
            ALLOCATE( rfu(jpi,jpj) , rfv(jpi,jpj) , STAT= ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'tra_adv_mle_init: failed to allocate arrays' )
            z1_t2 = 1._wp / ( rn_time * rn_time )
            DO jj = 2, jpj                           ! "coriolis+ time^-1" at u- & v-points
               DO ji = fs_2, jpi   ! vector opt.
                  zfu = ( ff_f(ji,jj) + ff_f(ji,jj-1) ) * 0.5_wp
                  zfv = ( ff_f(ji,jj) + ff_f(ji-1,jj) ) * 0.5_wp
                  rfu(ji,jj) = SQRT(  zfu * zfu + z1_t2 )
                  rfv(ji,jj) = SQRT(  zfv * zfv + z1_t2 )
               END DO
            END DO
            CALL lbc_lnk_multi( 'tramle', rfu, 'U', 1. , rfv, 'V', 1. )
            !
         ELSEIF( nn_mle == 1 ) THEN           ! MLE array allocation & initialisation
            rc_f = rn_ce / (  5.e3_wp * 2._wp * omega * SIN( rad * rn_lat )  )
            !
         ENDIF
         !
         !                                ! 1/(f^2+tau^2)^1/2 at t-point (needed in both nn_mle case)
         ALLOCATE( r1_ft(jpi,jpj) , STAT= ierr )
         IF( ierr /= 0 )   CALL ctl_stop( 'tra_adv_mle_init: failed to allocate r1_ft array' )
         !
         z1_t2 = 1._wp / ( rn_time * rn_time )
         r1_ft(:,:) = 1._wp / SQRT(  ff_t(:,:) * ff_t(:,:) + z1_t2  )
         !
      ENDIF
      !
   END SUBROUTINE tra_mle_init

   !!==============================================================================
END MODULE tramle
