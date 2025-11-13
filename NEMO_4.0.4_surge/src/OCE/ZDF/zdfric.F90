MODULE zdfric
   !!======================================================================
   !!                       ***  MODULE  zdfric  ***
   !! Ocean physics:  vertical mixing coefficient compute from the local
   !!                 Richardson number dependent formulation
   !!======================================================================
   !! History :  OPA  !  1987-09  (P. Andrich)  Original code
   !!            4.0  !  1991-11  (G. Madec)
   !!            7.0  !  1996-01  (G. Madec)  complete rewriting of multitasking suppression of common work arrays
   !!            8.0  !  1997-06  (G. Madec)  complete rewriting of zdfmix
   !!   NEMO     1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.3.1!  2011-09  (P. Oddo) Mixed layer depth parameterization
   !!            4.0  !  2017-04  (G. Madec)  remove CPP ddm key & avm at t-point only 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_ric_init  : initialization, namelist read, & parameters control
   !!   zdf_ric       : update momentum and tracer Kz from the Richardson number
   !!   ric_rst       : read/write RIC restart in ocean restart file
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE zdf_oce        ! vertical physics: variables
   USE phycst         ! physical constants
   USE sbc_oce,  ONLY :   taum
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)  


   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_ric         ! called by zdfphy.F90
   PUBLIC   ric_rst         ! called by zdfphy.F90
   PUBLIC   zdf_ric_init    ! called by nemogcm.F90

   !                        !!* Namelist namzdf_ric : Richardson number dependent Kz *
   INTEGER  ::   nn_ric      ! coefficient of the parameterization
   REAL(wp) ::   rn_avmri    ! maximum value of the vertical eddy viscosity
   REAL(wp) ::   rn_alp      ! coefficient of the parameterization
   REAL(wp) ::   rn_ekmfc    ! Ekman Factor Coeff
   REAL(wp) ::   rn_mldmin   ! minimum mixed layer (ML) depth    
   REAL(wp) ::   rn_mldmax   ! maximum mixed layer depth
   REAL(wp) ::   rn_wtmix    ! Vertical eddy Diff. in the ML
   REAL(wp) ::   rn_wvmix    ! Vertical eddy Visc. in the ML
   LOGICAL  ::   ln_mldw     ! Use or not the MLD parameters

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE zdf_ric_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE zdf_ric_init  ***
      !!                    
      !! ** Purpose :   Initialization of the vertical eddy diffusivity and
      !!      viscosity coef. for the Richardson number dependent formulation.
      !!
      !! ** Method  :   Read the namzdf_ric namelist and check the parameter values
      !!
      !! ** input   :   Namelist namzdf_ric
      !!
      !! ** Action  :   increase by 1 the nstop flag is setting problem encounter
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      INTEGER ::   ios          ! Local integer output status for namelist read
      !!
      NAMELIST/namzdf_ric/ rn_avmri, rn_alp   , nn_ric  , rn_ekmfc,  &
         &                rn_mldmin, rn_mldmax, rn_wtmix, rn_wvmix, ln_mldw
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namzdf_ric in reference namelist : Vertical diffusion Kz depends on Richardson number
      READ  ( numnam_ref, namzdf_ric, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_ric in reference namelist' )

      REWIND( numnam_cfg )              ! Namelist namzdf_ric in configuration namelist : Vertical diffusion Kz depends on Richardson number
      READ  ( numnam_cfg, namzdf_ric, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namzdf_ric in configuration namelist' )
      IF(lwm) WRITE ( numond, namzdf_ric )
      !
      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_ric_init : Ri depend vertical mixing scheme'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namzdf_ric : set Kz=F(Ri) parameters'
         WRITE(numout,*) '      maximum vertical viscosity        rn_avmri  = ', rn_avmri
         WRITE(numout,*) '      coefficient                       rn_alp    = ', rn_alp
         WRITE(numout,*) '      exponent                          nn_ric    = ', nn_ric
         WRITE(numout,*) '      Ekman layer enhanced mixing       ln_mldw   = ', ln_mldw
         WRITE(numout,*) '         Ekman Factor Coeff             rn_ekmfc  = ', rn_ekmfc
         WRITE(numout,*) '         minimum mixed layer depth      rn_mldmin = ', rn_mldmin
         WRITE(numout,*) '         maximum mixed layer depth      rn_mldmax = ', rn_mldmax
         WRITE(numout,*) '         Vertical eddy Diff. in the ML  rn_wtmix  = ', rn_wtmix
         WRITE(numout,*) '         Vertical eddy Visc. in the ML  rn_wvmix  = ', rn_wvmix
      ENDIF
      !
      CALL ric_rst( nit000, 'READ' )  !* read or initialize all required files
      !
      IF( lwxios ) THEN
         CALL iom_set_rstw_var_active('avt_k')
         CALL iom_set_rstw_var_active('avm_k')
      ENDIF
   END SUBROUTINE zdf_ric_init


   SUBROUTINE zdf_ric( kt, pdept, p_sh2, p_avm, p_avt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE zdfric  ***
      !!                    
      !! ** Purpose :   Compute the before eddy viscosity and diffusivity as
      !!                a function of the local richardson number.
      !!
      !! ** Method  :   Local richardson number dependent formulation of the 
      !!                vertical eddy viscosity and diffusivity coefficients. 
      !!                The eddy coefficients are given by:
      !!                    avm = avm0 + avmb
      !!                    avt = avm0 / (1 + rn_alp*ri)
      !!                with ri  = N^2 / dz(u)**2
      !!                         = e3w**2 * rn2/[ mi( dk(ub) )+mj( dk(vb) ) ]
      !!                    avm0= rn_avmri / (1 + rn_alp*Ri)**nn_ric
      !!                where ri is the before local Richardson number,
      !!                rn_avmri is the maximum value reaches by avm and avt 
      !!                and rn_alp, nn_ric are adjustable parameters.
      !!                Typical values : rn_alp=5. and nn_ric=2.
      !!
      !!      As second step compute Ekman depth from wind stress forcing
      !!      and apply namelist provided vertical coeff within this depth.
      !!      The Ekman depth is:
      !!              Ustar = SQRT(Taum/rho0)
      !!              ekd= rn_ekmfc * Ustar / f0
      !!      Large et al. (1994, eq.24) suggest rn_ekmfc=0.7; however, the derivation
      !!      of the above equation indicates the value is somewhat arbitrary; therefore
      !!      we allow the freedom to increase or decrease this value, if the
      !!      Ekman depth estimate appears too shallow or too deep, respectively.
      !!      Ekd is then limited by rn_mldmin and rn_mldmax provided in the
      !!      namelist
      !!        N.B. the mask are required for implicit scheme, and surface
      !!      and bottom value already set in zdfphy.F90
      !!
      !! ** Action  :   avm, avt  mixing coeff (inner domain values only)
      !!
      !! References : Pacanowski & Philander 1981, JPO, 1441-1451.
      !!              PFJ Lermusiaux 2001.
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   kt             ! ocean time-step
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdept          ! depth of t-point  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   p_sh2          ! shear production term
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   p_avm, p_avt   ! momentum and tracer Kz (w-points)
      !!
      INTEGER  ::   ji, jj, jk                  ! dummy loop indices
      REAL(wp) ::   zcfRi, zav, zustar, zhek    ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   zh_ekm  ! 2D workspace
      !!----------------------------------------------------------------------
      !
      !                       !==  avm and avt = F(Richardson number)  ==!
      DO jk = 2, jpkm1
         DO jj = 1, jpjm1
            DO ji = 1, jpim1              ! coefficient = F(richardson number) (avm-weighted Ri)
               zcfRi = 1._wp / (  1._wp + rn_alp * MAX(  0._wp , avm(ji,jj,jk) * rn2(ji,jj,jk) / ( p_sh2(ji,jj,jk) + 1.e-20 ) )  )
               zav   = rn_avmri * zcfRi**nn_ric
               !                          ! avm and avt coefficients
               p_avm(ji,jj,jk) = MAX(  zav         , avmb(jk)  ) * wmask(ji,jj,jk)
               p_avt(ji,jj,jk) = MAX(  zav * zcfRi , avtb(jk)  ) * wmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
!!gm BUG <<<<====  This param can't work at low latitude 
!!gm               it provides there much to thick mixed layer ( summer 150m in GYRE configuration !!! )
      !
      IF( ln_mldw ) THEN      !==  set a minimum value in the Ekman layer  ==!
         !
         DO jj = 2, jpjm1        !* Ekman depth
            DO ji = 2, jpim1
               zustar = SQRT( taum(ji,jj) * r1_rau0 )
               zhek   = rn_ekmfc * zustar / ( ABS( ff_t(ji,jj) ) + rsmall )   ! Ekman depth
               zh_ekm(ji,jj) = MAX(  rn_mldmin , MIN( zhek , rn_mldmax )  )   ! set allowed range
            END DO
         END DO
         DO jk = 2, jpkm1        !* minimum mixing coeff. within the Ekman layer
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  IF( pdept(ji,jj,jk) < zh_ekm(ji,jj) ) THEN
                     p_avm(ji,jj,jk) = MAX(  p_avm(ji,jj,jk), rn_wvmix  ) * wmask(ji,jj,jk)
                     p_avt(ji,jj,jk) = MAX(  p_avt(ji,jj,jk), rn_wtmix  ) * wmask(ji,jj,jk)
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE zdf_ric


   SUBROUTINE ric_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ric_rst  ***
      !!                     
      !! ** Purpose :   Read or write TKE file (en) in restart file
      !!
      !! ** Method  :   use of IOM library
      !!                if the restart does not contain TKE, en is either 
      !!                set to rn_emin or recomputed 
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt     ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      !
      INTEGER ::   jit, jk    ! dummy loop indices
      INTEGER ::   id1, id2   ! local integers
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
         !                                   ! ---------------
         !           !* Read the restart file
         IF( ln_rstart ) THEN
            id1 = iom_varid( numror, 'avt_k', ldstop = .FALSE. )
            id2 = iom_varid( numror, 'avm_k', ldstop = .FALSE. )
            !
            IF( MIN( id1, id2 ) > 0 ) THEN         ! restart exists => read it
               CALL iom_get( numror, jpdom_autoglo, 'avt_k', avt_k, ldxios = lrxios )
               CALL iom_get( numror, jpdom_autoglo, 'avm_k', avm_k, ldxios = lrxios )
            ENDIF
         ENDIF
         !           !* otherwise Kz already set to the background value in zdf_phy_init
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*) '---- ric-rst ----'
         IF( lwxios ) CALL iom_swap(      cwxios_context          )
         CALL iom_rstput( kt, nitrst, numrow, 'avt_k', avt_k, ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'avm_k', avm_k, ldxios = lwxios)
         IF( lwxios ) CALL iom_swap(      cxios_context          )
         !
      ENDIF
      !
   END SUBROUTINE ric_rst

   !!======================================================================
END MODULE zdfric
