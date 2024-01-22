MODULE ldfdyn
   !!======================================================================
   !!                       ***  MODULE  ldfdyn  ***
   !! Ocean physics:  lateral viscosity coefficient 
   !!=====================================================================
   !! History :  OPA  ! 1997-07  (G. Madec)  multi dimensional coefficients
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            3.7  ! 2014-01  (F. Lemarie, G. Madec)  restructuration/simplification of ahm specification,
   !!                 !                                  add velocity dependent coefficient and optional read in file
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_dyn_init  : initialization, namelist read, and parameters control
   !!   ldf_dyn       : update lateral eddy viscosity coefficients at each time step 
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers   
   USE dom_oce         ! ocean space and time domain 
   USE phycst          ! physical constants
   USE ldfslp          ! lateral diffusion: slopes of mixing orientation
   USE ldfc1d_c2d      ! lateral diffusion: 1D and 2D cases
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O module for ehanced bottom friction file
   USE timing          ! Timing
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_dyn_init   ! called by nemogcm.F90
   PUBLIC   ldf_dyn        ! called by step.F90

   !                                    !!* Namelist namdyn_ldf : lateral mixing on momentum *
   LOGICAL , PUBLIC ::   ln_dynldf_OFF   !: No operator (i.e. no explicit diffusion)
   LOGICAL , PUBLIC ::   ln_dynldf_lap   !: laplacian operator
   LOGICAL , PUBLIC ::   ln_dynldf_blp   !: bilaplacian operator
   LOGICAL , PUBLIC ::   ln_dynldf_lev   !: iso-level direction
   LOGICAL , PUBLIC ::   ln_dynldf_hor   !: horizontal (geopotential) direction
!  LOGICAL , PUBLIC ::   ln_dynldf_iso   !: iso-neutral direction                        (see ldfslp)
   INTEGER , PUBLIC ::   nn_ahm_ijk_t    !: choice of time & space variations of the lateral eddy viscosity coef.
   !                                        !  time invariant coefficients:  aht = 1/2  Ud*Ld   (lap case) 
   !                                           !                             bht = 1/12 Ud*Ld^3 (blp case)
   REAL(wp), PUBLIC ::   rn_Uv                 !: lateral viscous velocity  [m/s]
   REAL(wp), PUBLIC ::   rn_Lv                 !: lateral viscous length    [m]
   !                                        ! Smagorinsky viscosity  (nn_ahm_ijk_t = 32) 
   REAL(wp), PUBLIC ::   rn_csmc               !: Smagorinsky constant of proportionality 
   REAL(wp), PUBLIC ::   rn_minfac             !: Multiplicative factor of theorectical minimum Smagorinsky viscosity
   REAL(wp), PUBLIC ::   rn_maxfac             !: Multiplicative factor of theorectical maximum Smagorinsky viscosity
   !                                        ! iso-neutral laplacian (ln_dynldf_lap=ln_dynldf_iso=T)
   REAL(wp), PUBLIC ::   rn_ahm_b              !: lateral laplacian background eddy viscosity  [m2/s]

   !                                    !!* Parameter to control the type of lateral viscous operator
   INTEGER, PARAMETER, PUBLIC ::   np_ERROR  =-10                       !: error in setting the operator
   INTEGER, PARAMETER, PUBLIC ::   np_no_ldf = 00                       !: without operator (i.e. no lateral viscous trend)
   !                          !!      laplacian     !    bilaplacian    !
   INTEGER, PARAMETER, PUBLIC ::   np_lap    = 10   ,   np_blp    = 20  !: iso-level operator
   INTEGER, PARAMETER, PUBLIC ::   np_lap_i  = 11                       !: iso-neutral or geopotential operator
   !
   INTEGER           , PUBLIC ::   nldf_dyn         !: type of lateral diffusion used defined from ln_dynldf_... (namlist logicals)
   LOGICAL           , PUBLIC ::   l_ldfdyn_time    !: flag for time variation of the lateral eddy viscosity coef.

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ahmt, ahmf   !: eddy viscosity coef. at T- and F-points [m2/s or m4/s]
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dtensq       !: horizontal tension squared         (Smagorinsky only)
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dshesq       !: horizontal shearing strain squared (Smagorinsky only)
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   esqt, esqf   !: Square of the local gridscale (e1e2/(e1+e2))**2           

   REAL(wp) ::   r1_2    = 0.5_wp            ! =1/2
   REAL(wp) ::   r1_4    = 0.25_wp           ! =1/4
   REAL(wp) ::   r1_8    = 0.125_wp          ! =1/8
   REAL(wp) ::   r1_12   = 1._wp / 12._wp    ! =1/12
   REAL(wp) ::   r1_288  = 1._wp / 288._wp   ! =1/( 12^2 * 2 )

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_dyn_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_dyn_init  ***
      !!                   
      !! ** Purpose :   set the horizontal ocean dynamics physics
      !!
      !! ** Method  :   the eddy viscosity coef. specification depends on:
      !!              - the operator:
      !!             ln_dynldf_lap = T     laplacian operator
      !!             ln_dynldf_blp = T   bilaplacian operator
      !!              - the parameter nn_ahm_ijk_t:
      !!    nn_ahm_ijk_t  =  0 => = constant
      !!                  = 10 => = F(z) :     = constant with a reduction of 1/4 with depth 
      !!                  =-20 => = F(i,j)     = shape read in 'eddy_viscosity.nc' file
      !!                  = 20    = F(i,j)     = F(e1,e2) or F(e1^3,e2^3) (lap or bilap case)
      !!                  =-30 => = F(i,j,k)   = shape read in 'eddy_viscosity.nc'  file
      !!                  = 30    = F(i,j,k)   = 2D (case 20) + decrease with depth (case 10)
      !!                  = 31    = F(i,j,k,t) = F(local velocity) (  |u|e  /12   laplacian operator
      !!                                                           or |u|e^3/12 bilaplacian operator )
      !!                  = 32    = F(i,j,k,t) = F(local deformation rate and gridscale) (D and L) (Smagorinsky)  
      !!                                                           (   L^2|D|      laplacian operator
      !!                                                           or  L^4|D|/8  bilaplacian operator )
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk                     ! dummy loop indices
      INTEGER  ::   ioptio, ierr, inum, ios, inn   ! local integer
      REAL(wp) ::   zah0, zah_max, zUfac           ! local scalar
      CHARACTER(len=5) ::   cl_Units               ! units (m2/s or m4/s)
      !!
      NAMELIST/namdyn_ldf/ ln_dynldf_OFF, ln_dynldf_lap, ln_dynldf_blp,   &   ! type of operator
         &                 ln_dynldf_lev, ln_dynldf_hor, ln_dynldf_iso,   &   ! acting direction of the operator
         &                 nn_ahm_ijk_t , rn_Uv    , rn_Lv,   rn_ahm_b,   &   ! lateral eddy coefficient
         &                 rn_csmc      , rn_minfac    , rn_maxfac            ! Smagorinsky settings
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )
      READ  ( numnam_ref, namdyn_ldf, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdyn_ldf in reference namelist' )

      REWIND( numnam_cfg )
      READ  ( numnam_cfg, namdyn_ldf, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdyn_ldf in configuration namelist' )
      IF(lwm) WRITE ( numond, namdyn_ldf )

      IF(lwp) THEN                      ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_dyn : lateral momentum physics'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_ldf : set lateral mixing parameters'
         !
         WRITE(numout,*) '      type :'
         WRITE(numout,*) '         no explicit diffusion                ln_dynldf_OFF = ', ln_dynldf_OFF
         WRITE(numout,*) '         laplacian operator                   ln_dynldf_lap = ', ln_dynldf_lap
         WRITE(numout,*) '         bilaplacian operator                 ln_dynldf_blp = ', ln_dynldf_blp
         !
         WRITE(numout,*) '      direction of action :'
         WRITE(numout,*) '         iso-level                            ln_dynldf_lev = ', ln_dynldf_lev
         WRITE(numout,*) '         horizontal (geopotential)            ln_dynldf_hor = ', ln_dynldf_hor
         WRITE(numout,*) '         iso-neutral                          ln_dynldf_iso = ', ln_dynldf_iso
         !
         WRITE(numout,*) '      coefficients :'
         WRITE(numout,*) '         type of time-space variation         nn_ahm_ijk_t  = ', nn_ahm_ijk_t
         WRITE(numout,*) '         lateral viscous velocity  (if cst)      rn_Uv      = ', rn_Uv, ' m/s'
         WRITE(numout,*) '         lateral viscous length    (if cst)      rn_Lv      = ', rn_Lv, ' m'
         WRITE(numout,*) '         background viscosity (iso-lap case)     rn_ahm_b   = ', rn_ahm_b, ' m2/s'
         !
         WRITE(numout,*) '      Smagorinsky settings (nn_ahm_ijk_t  = 32) :'
         WRITE(numout,*) '         Smagorinsky coefficient              rn_csmc       = ', rn_csmc
         WRITE(numout,*) '         factor multiplier for eddy visc.'
         WRITE(numout,*) '            lower limit (default 1.0)         rn_minfac    = ', rn_minfac
         WRITE(numout,*) '            upper limit (default 1.0)         rn_maxfac    = ', rn_maxfac
      ENDIF

      !
      !           !==  type of lateral operator used  ==!   (set nldf_dyn)
      !           !=====================================!
      !
      nldf_dyn = np_ERROR
      ioptio = 0
      IF( ln_dynldf_OFF ) THEN   ;   nldf_dyn = np_no_ldf   ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ln_dynldf_lap ) THEN   ;                              ioptio = ioptio + 1   ;   ENDIF
      IF( ln_dynldf_blp ) THEN   ;                              ioptio = ioptio + 1   ;   ENDIF
      IF( ioptio /= 1   )   CALL ctl_stop( 'dyn_ldf_init: use ONE of the 3 operator options (NONE/lap/blp)' )
      !
      IF(.NOT.ln_dynldf_OFF ) THEN     !==  direction ==>> type of operator  ==!
         ioptio = 0
         IF( ln_dynldf_lev )   ioptio = ioptio + 1
         IF( ln_dynldf_hor )   ioptio = ioptio + 1
         IF( ln_dynldf_iso )   ioptio = ioptio + 1
         IF( ioptio /= 1   )   CALL ctl_stop( 'dyn_ldf_init: use ONE of the 3 direction options (level/hor/iso)' )
         !
         !                             ! Set nldf_dyn, the type of lateral diffusion, from ln_dynldf_... logicals
         ierr = 0
         IF( ln_dynldf_lap ) THEN         ! laplacian operator
            IF( ln_zco ) THEN                   ! z-coordinate
               IF ( ln_dynldf_lev )   nldf_dyn = np_lap     ! iso-level = horizontal (no rotation)
               IF ( ln_dynldf_hor )   nldf_dyn = np_lap     ! iso-level = horizontal (no rotation)
               IF ( ln_dynldf_iso )   nldf_dyn = np_lap_i   ! iso-neutral            (   rotation)
            ENDIF
            IF( ln_zps ) THEN                   ! z-coordinate with partial step
               IF ( ln_dynldf_lev )   nldf_dyn = np_lap     ! iso-level              (no rotation)
               IF ( ln_dynldf_hor )   nldf_dyn = np_lap     ! iso-level              (no rotation)
               IF ( ln_dynldf_iso )   nldf_dyn = np_lap_i   ! iso-neutral            (   rotation)
            ENDIF
            IF( ln_sco ) THEN                   ! s-coordinate
               IF ( ln_dynldf_lev )   nldf_dyn = np_lap     ! iso-level = horizontal (no rotation)
               IF ( ln_dynldf_hor )   nldf_dyn = np_lap_i   ! horizontal             (   rotation)
               IF ( ln_dynldf_iso )   nldf_dyn = np_lap_i   ! iso-neutral            (   rotation)
            ENDIF
         ENDIF
         !
         IF( ln_dynldf_blp ) THEN         ! bilaplacian operator
            IF( ln_zco ) THEN                   ! z-coordinate
               IF( ln_dynldf_lev )   nldf_dyn = np_blp   ! iso-level = horizontal (no rotation)
               IF( ln_dynldf_hor )   nldf_dyn = np_blp   ! iso-level = horizontal (no rotation)
               IF( ln_dynldf_iso )   ierr = 2            ! iso-neutral            (   rotation)
            ENDIF
            IF( ln_zps ) THEN                   ! z-coordinate with partial step
               IF( ln_dynldf_lev )   nldf_dyn = np_blp   ! iso-level              (no rotation)
               IF( ln_dynldf_hor )   nldf_dyn = np_blp   ! iso-level              (no rotation)
               IF( ln_dynldf_iso )   ierr = 2            ! iso-neutral            (   rotation)
            ENDIF
            IF( ln_sco ) THEN                   ! s-coordinate
               IF( ln_dynldf_lev )   nldf_dyn = np_blp   ! iso-level              (no rotation)
               IF( ln_dynldf_hor )   ierr = 2            ! horizontal             (   rotation)
               IF( ln_dynldf_iso )   ierr = 2            ! iso-neutral            (   rotation)
            ENDIF
         ENDIF
         !
         IF( ierr == 2 )   CALL ctl_stop( 'rotated bi-laplacian operator does not exist' )
         !
         IF( nldf_dyn == np_lap_i )   l_ldfslp = .TRUE.  ! rotation require the computation of the slopes
         !
      ENDIF
      !
      IF(lwp) THEN
         WRITE(numout,*)
         SELECT CASE( nldf_dyn )
         CASE( np_no_ldf )   ;   WRITE(numout,*) '   ==>>>   NO lateral viscosity'
         CASE( np_lap    )   ;   WRITE(numout,*) '   ==>>>   iso-level laplacian operator'
         CASE( np_lap_i  )   ;   WRITE(numout,*) '   ==>>>   rotated laplacian operator with iso-level background'
         CASE( np_blp    )   ;   WRITE(numout,*) '   ==>>>   iso-level bi-laplacian operator'
         END SELECT
         WRITE(numout,*)
      ENDIF
      
      !
      !           !==  Space/time variation of eddy coefficients  ==!
      !           !=================================================!
      !
      l_ldfdyn_time = .FALSE.                ! no time variation except in case defined below
      !
      IF( ln_dynldf_OFF ) THEN
         IF(lwp) WRITE(numout,*) '   ==>>>   No viscous operator selected. ahmt and ahmf are not allocated'
         RETURN
         !
      ELSE                                   !==  a lateral diffusion operator is used  ==!
         !
         !                                         ! allocate the ahm arrays
         ALLOCATE( ahmt(jpi,jpj,jpk) , ahmf(jpi,jpj,jpk) , STAT=ierr )
         IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'ldf_dyn_init: failed to allocate arrays')
         !
         ahmt(:,:,:) = 0._wp                       ! init to 0 needed 
         ahmf(:,:,:) = 0._wp
         !
         !                                         ! value of lap/blp eddy mixing coef.
         IF(     ln_dynldf_lap ) THEN   ;   zUfac = r1_2 *rn_Uv   ;   inn = 1   ;   cl_Units = ' m2/s'   !   laplacian
         ELSEIF( ln_dynldf_blp ) THEN   ;   zUfac = r1_12*rn_Uv   ;   inn = 3   ;   cl_Units = ' m4/s'   ! bilaplacian
         ENDIF
         zah0    = zUfac *    rn_Lv**inn              ! mixing coefficient
         zah_max = zUfac * (ra*rad)**inn              ! maximum reachable coefficient (value at the Equator)
         !
         SELECT CASE(  nn_ahm_ijk_t  )             !* Specification of space-time variations of ahmt, ahmf
         !
         CASE(   0  )      !==  constant  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity. = constant = ', zah0, cl_Units
            ahmt(:,:,1:jpkm1) = zah0
            ahmf(:,:,1:jpkm1) = zah0
            !
         CASE(  10  )      !==  fixed profile  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F( depth )'
            IF(lwp) WRITE(numout,*) '           surface viscous coef. = constant = ', zah0, cl_Units
            ahmt(:,:,1) = zah0                        ! constant surface value
            ahmf(:,:,1) = zah0
            CALL ldf_c1d( 'DYN', ahmt(:,:,1), ahmf(:,:,1), ahmt, ahmf )
            !
         CASE ( -20 )      !== fixed horizontal shape read in file  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F(i,j) read in eddy_viscosity.nc file'
            CALL iom_open( 'eddy_viscosity_2D.nc', inum )
            CALL iom_get ( inum, jpdom_data, 'ahmt_2d', ahmt(:,:,1) )
            CALL iom_get ( inum, jpdom_data, 'ahmf_2d', ahmf(:,:,1) )
            CALL iom_close( inum )
            DO jk = 2, jpkm1
               ahmt(:,:,jk) = ahmt(:,:,1)
               ahmf(:,:,jk) = ahmf(:,:,1)
            END DO
            !
         CASE(  20  )      !== fixed horizontal shape  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F( e1, e2 ) or F( e1^3, e2^3 ) (lap. or blp. case)'
            IF(lwp) WRITE(numout,*) '           using a fixed viscous velocity = ', rn_Uv  ,' m/s   and   Lv = Max(e1,e2)'
            IF(lwp) WRITE(numout,*) '           maximum reachable coefficient (at the Equator) = ', zah_max, cl_Units, '  for e1=1°)'
            CALL ldf_c2d( 'DYN', zUfac      , inn        , ahmt, ahmf )         ! surface value proportional to scale factor^inn
            !
         CASE( -30  )      !== fixed 3D shape read in file  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F(i,j,k) read in eddy_viscosity_3D.nc file'
            CALL iom_open( 'eddy_viscosity_3D.nc', inum )
            CALL iom_get ( inum, jpdom_data, 'ahmt_3d', ahmt )
            CALL iom_get ( inum, jpdom_data, 'ahmf_3d', ahmf )
            CALL iom_close( inum )
            !
         CASE(  30  )       !==  fixed 3D shape  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F( latitude, longitude, depth )'
            IF(lwp) WRITE(numout,*) '           using a fixed viscous velocity = ', rn_Uv  ,' m/s   and   Ld = Max(e1,e2)'
            IF(lwp) WRITE(numout,*) '           maximum reachable coefficient (at the Equator) = ', zah_max, cl_Units, '  for e1=1°)'
            CALL ldf_c2d( 'DYN', zUfac      , inn        , ahmt, ahmf )         ! surface value proportional to scale factor^inn
            CALL ldf_c1d( 'DYN', ahmt(:,:,1), ahmf(:,:,1), ahmt, ahmf )  ! reduction with depth
            !
         CASE(  31  )       !==  time varying 3D field  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F( latitude, longitude, depth , time )'
            IF(lwp) WRITE(numout,*) '           proportional to the local velocity : 1/2 |u|e (lap) or 1/12 |u|e^3 (blp)'
            !
            l_ldfdyn_time = .TRUE.     ! will be calculated by call to ldf_dyn routine in step.F90
            !
         CASE(  32  )       !==  time varying 3D field  ==!
            IF(lwp) WRITE(numout,*) '   ==>>>   eddy viscosity = F( latitude, longitude, depth , time )'
            IF(lwp) WRITE(numout,*) '           proportional to the local deformation rate and gridscale (Smagorinsky)'
            !
            l_ldfdyn_time = .TRUE.     ! will be calculated by call to ldf_dyn routine in step.F90
            !
            !                          ! allocate arrays used in ldf_dyn. 
            ALLOCATE( dtensq(jpi,jpj,jpk) , dshesq(jpi,jpj,jpk) , esqt(jpi,jpj) , esqf(jpi,jpj) , STAT=ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'ldf_dyn_init: failed to allocate Smagorinsky arrays')
            !
            DO jj = 1, jpj             ! Set local gridscale values
               DO ji = 1, jpi
                  esqt(ji,jj) = ( 2._wp * e1e2t(ji,jj) / ( e1t(ji,jj) + e2t(ji,jj) ) )**2 
                  esqf(ji,jj) = ( 2._wp * e1e2f(ji,jj) / ( e1f(ji,jj) + e2f(ji,jj) ) )**2 
               END DO
            END DO
            !
         CASE DEFAULT
            CALL ctl_stop('ldf_dyn_init: wrong choice for nn_ahm_ijk_t, the type of space-time variation of ahm')
         END SELECT
         !
         IF( .NOT.l_ldfdyn_time ) THEN             !* No time variation 
            IF(     ln_dynldf_lap ) THEN                 !   laplacian operator (mask only)
               ahmt(:,:,1:jpkm1) =       ahmt(:,:,1:jpkm1)   * tmask(:,:,1:jpkm1)
               ahmf(:,:,1:jpkm1) =       ahmf(:,:,1:jpkm1)   * fmask(:,:,1:jpkm1)
            ELSEIF( ln_dynldf_blp ) THEN                 ! bilaplacian operator (square root + mask)
               ahmt(:,:,1:jpkm1) = SQRT( ahmt(:,:,1:jpkm1) ) * tmask(:,:,1:jpkm1)
               ahmf(:,:,1:jpkm1) = SQRT( ahmf(:,:,1:jpkm1) ) * fmask(:,:,1:jpkm1)
            ENDIF
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE ldf_dyn_init


   SUBROUTINE ldf_dyn( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_dyn  ***
      !! 
      !! ** Purpose :   update at kt the momentum lateral mixing coeff. (ahmt and ahmf)
      !!
      !! ** Method  :   time varying eddy viscosity coefficients:
      !!
      !!    nn_ahm_ijk_t = 31    ahmt, ahmf = F(i,j,k,t) = F(local velocity) 
      !!                         ( |u|e /12  or  |u|e^3/12 for laplacian or bilaplacian operator )
      !!
      !!    nn_ahm_ijk_t = 32    ahmt, ahmf = F(i,j,k,t) = F(local deformation rate and gridscale) (D and L) (Smagorinsky)  
      !!                         ( L^2|D|    or  L^4|D|/8  for laplacian or bilaplacian operator )
      !!
      !! ** note    :    in BLP cases the sqrt of the eddy coef is returned, since bilaplacian is en re-entrant laplacian
      !! ** action  :    ahmt, ahmf   updated at each time step
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time step index
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zu2pv2_ij_p1, zu2pv2_ij, zu2pv2_ij_m1, zemax   ! local scalar (option 31)
      REAL(wp) ::   zcmsmag, zstabf_lo, zstabf_up, zdelta, zdb     ! local scalar (option 32)
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('ldf_dyn')
      !
      SELECT CASE(  nn_ahm_ijk_t  )       !== Eddy vicosity coefficients ==!
      !
      CASE(  31  )       !==  time varying 3D field  ==!   = F( local velocity )
         !
         IF( ln_dynldf_lap   ) THEN        ! laplacian operator : |u| e /12 = |u/144| e
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     zu2pv2_ij_m1 = ub(ji-1,jj  ,jk) * ub(ji-1,jj  ,jk) + vb(ji  ,jj-1,jk) * vb(ji  ,jj-1,jk)
                     zemax = MAX( e1t(ji,jj) , e2t(ji,jj) )
                     ahmt(ji,jj,jk) = SQRT( (zu2pv2_ij + zu2pv2_ij_m1) * r1_288 ) * zemax * tmask(ji,jj,jk)      ! 288= 12*12 * 2
                  END DO
               END DO
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     zu2pv2_ij_p1 = ub(ji  ,jj+1,jk) * ub(ji  ,jj+1,jk) + vb(ji+1,jj  ,jk) * vb(ji+1,jj  ,jk)
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     zemax = MAX( e1f(ji,jj) , e2f(ji,jj) )
                     ahmf(ji,jj,jk) = SQRT( (zu2pv2_ij + zu2pv2_ij_p1) * r1_288 ) * zemax * fmask(ji,jj,jk)      ! 288= 12*12 * 2
                  END DO
               END DO
            END DO
         ELSEIF( ln_dynldf_blp ) THEN      ! bilaplacian operator : sqrt( |u| e^3 /12 ) = sqrt( |u/144| e ) * e
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     zu2pv2_ij_m1 = ub(ji-1,jj  ,jk) * ub(ji-1,jj  ,jk) + vb(ji  ,jj-1,jk) * vb(ji  ,jj-1,jk)
                     zemax = MAX( e1t(ji,jj) , e2t(ji,jj) )
                     ahmt(ji,jj,jk) = SQRT(  SQRT( (zu2pv2_ij + zu2pv2_ij_m1) * r1_288 ) * zemax  ) * zemax * tmask(ji,jj,jk)
                  END DO
               END DO
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     zu2pv2_ij_p1 = ub(ji  ,jj+1,jk) * ub(ji  ,jj+1,jk) + vb(ji+1,jj  ,jk) * vb(ji+1,jj  ,jk)
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     zemax = MAX( e1f(ji,jj) , e2f(ji,jj) )
                     ahmf(ji,jj,jk) = SQRT(  SQRT( (zu2pv2_ij + zu2pv2_ij_p1) * r1_288 ) * zemax  ) * zemax * fmask(ji,jj,jk)
                  END DO
               END DO
            END DO
         ENDIF
         !
         CALL lbc_lnk_multi( 'ldfdyn', ahmt, 'T', 1.,  ahmf, 'F', 1. )
         !
         !
      CASE(  32  )       !==  time varying 3D field  ==!   = F( local deformation rate and gridscale ) (Smagorinsky)
         !
         IF( ln_dynldf_lap .OR. ln_dynldf_blp  ) THEN        ! laplacian operator : (C_smag/pi)^2 L^2 |D|
            !
            zcmsmag   = (rn_csmc/rpi)**2                                            ! (C_smag/pi)^2
            zstabf_lo = rn_minfac * rn_minfac / ( 2._wp * 12._wp * 12._wp * zcmsmag ) ! lower limit stability factor scaling
            zstabf_up = rn_maxfac / ( 4._wp * zcmsmag * 2._wp * rdt )               ! upper limit stability factor scaling
            IF( ln_dynldf_blp ) zstabf_lo = ( 16._wp / 9._wp ) * zstabf_lo          ! provide |U|L^3/12 lower limit instead 
            !                                                                       ! of |U|L^3/16 in blp case
            DO jk = 1, jpkm1
               !
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1
                     zdb =   ( ub(ji,jj,jk) * r1_e2u(ji,jj) - ub(ji-1,jj,jk) * r1_e2u(ji-1,jj) ) * r1_e1t(ji,jj) * e2t(ji,jj)  &
                        &  - ( vb(ji,jj,jk) * r1_e1v(ji,jj) - vb(ji,jj-1,jk) * r1_e1v(ji,jj-1) ) * r1_e2t(ji,jj) * e1t(ji,jj)
                     dtensq(ji,jj,jk) = zdb * zdb * tmask(ji,jj,jk)
                  END DO
               END DO
               !
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     zdb =   ( ub(ji,jj+1,jk) * r1_e1u(ji,jj+1) - ub(ji,jj,jk) * r1_e1u(ji,jj) ) * r1_e2f(ji,jj) * e1f(ji,jj)  &
                        &  + ( vb(ji+1,jj,jk) * r1_e2v(ji+1,jj) - vb(ji,jj,jk) * r1_e2v(ji,jj) ) * r1_e1f(ji,jj) * e2f(ji,jj)
                     dshesq(ji,jj,jk) = zdb * zdb * fmask(ji,jj,jk)
                  END DO
               END DO
               !
            END DO
            !
            CALL lbc_lnk_multi( 'ldfdyn', dtensq, 'T', 1. )  ! lbc_lnk on dshesq not needed
            !
            DO jk = 1, jpkm1
              !
               DO jj = 2, jpjm1                                ! T-point value
                  DO ji = fs_2, fs_jpim1
                     !
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     zu2pv2_ij_m1 = ub(ji-1,jj  ,jk) * ub(ji-1,jj  ,jk) + vb(ji  ,jj-1,jk) * vb(ji  ,jj-1,jk)
                     !
                     zdelta         = zcmsmag * esqt(ji,jj)                                        ! L^2 * (C_smag/pi)^2
                     ahmt(ji,jj,jk) = zdelta * SQRT(          dtensq(ji  ,jj,jk) +                         &
                        &                            r1_4 * ( dshesq(ji  ,jj,jk) + dshesq(ji  ,jj-1,jk) +  &
                        &                                     dshesq(ji-1,jj,jk) + dshesq(ji-1,jj-1,jk) ) )
                     ahmt(ji,jj,jk) = MAX( ahmt(ji,jj,jk), SQRT( (zu2pv2_ij + zu2pv2_ij_m1) * zdelta * zstabf_lo ) ) ! Impose lower limit == minfac  * |U|L/2
                     ahmt(ji,jj,jk) = MIN( ahmt(ji,jj,jk),                                    zdelta * zstabf_up )   ! Impose upper limit == maxfac  * L^2/(4*2dt)
                     !
                  END DO
               END DO
               !
               DO jj = 1, jpjm1                                ! F-point value
                  DO ji = 1, fs_jpim1
                     !
                     zu2pv2_ij_p1 = ub(ji  ,jj+1,jk) * ub(ji  ,jj+1,jk) + vb(ji+1,jj  ,jk) * vb(ji+1,jj  ,jk)
                     zu2pv2_ij    = ub(ji  ,jj  ,jk) * ub(ji  ,jj  ,jk) + vb(ji  ,jj  ,jk) * vb(ji  ,jj  ,jk)
                     !
                     zdelta         = zcmsmag * esqf(ji,jj)                                        ! L^2 * (C_smag/pi)^2
                     ahmf(ji,jj,jk) = zdelta * SQRT(          dshesq(ji  ,jj,jk) +                         &
                        &                            r1_4 * ( dtensq(ji  ,jj,jk) + dtensq(ji  ,jj+1,jk) +  &
                        &                                     dtensq(ji+1,jj,jk) + dtensq(ji+1,jj+1,jk) ) )
                     ahmf(ji,jj,jk) = MAX( ahmf(ji,jj,jk), SQRT( (zu2pv2_ij + zu2pv2_ij_p1) * zdelta * zstabf_lo ) ) ! Impose lower limit == minfac  * |U|L/2
                     ahmf(ji,jj,jk) = MIN( ahmf(ji,jj,jk),                                    zdelta * zstabf_up )   ! Impose upper limit == maxfac  * L^2/(4*2dt)
                     !
                  END DO
               END DO
               !
            END DO
            !
         ENDIF
         !
         IF( ln_dynldf_blp ) THEN      ! bilaplacian operator : sqrt( (C_smag/pi)^2 L^4 |D|/8)
            !                          !                      = sqrt( A_lap_smag L^2/8 )
            !                          ! stability limits already applied to laplacian values
            !                          ! effective default limits are 1/12 |U|L^3 < B_hm < 1//(32*2dt) L^4
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     ahmt(ji,jj,jk) = SQRT( r1_8 * esqt(ji,jj) * ahmt(ji,jj,jk) )
                  END DO
               END DO
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1
                     ahmf(ji,jj,jk) = SQRT( r1_8 * esqf(ji,jj) * ahmf(ji,jj,jk) )
                  END DO
               END DO
            END DO
            !
         ENDIF
         !
         CALL lbc_lnk_multi( 'ldfdyn', ahmt, 'T', 1. , ahmf, 'F', 1. )
         !
      END SELECT
      !
      CALL iom_put( "ahmt_2d", ahmt(:,:,1) )   ! surface u-eddy diffusivity coeff.
      CALL iom_put( "ahmf_2d", ahmf(:,:,1) )   ! surface v-eddy diffusivity coeff.
      CALL iom_put( "ahmt_3d", ahmt(:,:,:) )   ! 3D      u-eddy diffusivity coeff.
      CALL iom_put( "ahmf_3d", ahmf(:,:,:) )   ! 3D      v-eddy diffusivity coeff.
      !
      IF( ln_timing )   CALL timing_stop('ldf_dyn')
      !
   END SUBROUTINE ldf_dyn

   !!======================================================================
END MODULE ldfdyn
