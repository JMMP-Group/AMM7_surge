MODULE trdtra
   !!======================================================================
   !!                       ***  MODULE  trdtra  ***
   !! Ocean diagnostics:  ocean tracers trends pre-processing
   !!=====================================================================
   !! History :  3.3  !  2010-06  (C. Ethe) creation for the TRA/TRC merge
   !!            3.5  !  2012-02  (G. Madec) update the comments 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   trd_tra       : pre-process the tracer trends
   !!   trd_tra_adv   : transform a div(U.T) trend into a U.grad(T) trend
   !!   trd_tra_mng   : tracer trend manager: dispatch to the diagnostic modules
   !!   trd_tra_iom   : output 3D tracer trends using IOM
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean domain 
   USE sbc_oce        ! surface boundary condition: ocean
   USE zdf_oce        ! ocean vertical physics
   USE trd_oce        ! trends: ocean variables
   USE trdtrc         ! ocean passive mixed layer tracers trends 
   USE trdglo         ! trends: global domain averaged
   USE trdpen         ! trends: Potential ENergy
   USE trdmxl         ! ocean active mixed layer tracers trends 
   USE ldftra         ! ocean active tracers lateral physics
   USE ldfslp
   USE zdfddm         ! vertical physics: double diffusion
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trd_tra   ! called by all tra_... modules

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   trdtx, trdty, trdt   ! use to store the temperature trends
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   avt_evd  ! store avt_evd to calculate EVD trend

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trd_tra_alloc()
      !!---------------------------------------------------------------------
      !!                  ***  FUNCTION trd_tra_alloc  ***
      !!---------------------------------------------------------------------
      ALLOCATE( trdtx(jpi,jpj,jpk) , trdty(jpi,jpj,jpk) , trdt(jpi,jpj,jpk) , avt_evd(jpi,jpj,jpk), STAT= trd_tra_alloc )
      !
      CALL mpp_sum ( 'trdtra', trd_tra_alloc )
      IF( trd_tra_alloc /= 0 )   CALL ctl_stop( 'STOP', 'trd_tra_alloc: failed to allocate arrays' )
   END FUNCTION trd_tra_alloc


   SUBROUTINE trd_tra( kt, ctype, ktra, ktrd, ptrd, pun, ptra )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra  ***
      !! 
      !! ** Purpose : pre-process tracer trends
      !!
      !! ** Method  : - mask the trend
      !!              - advection (ptra present) converte the incoming flux (U.T) 
      !!              into trend (U.T => -U.grat(T)=div(U.T)-T.div(U)) through a 
      !!              call to trd_tra_adv
      !!              - 'TRA' case : regroup T & S trends
      !!              - send the trends to trd_tra_mng (trdtrc) for further processing
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in)           ::   kt      ! time step
      CHARACTER(len=3)                , INTENT(in)           ::   ctype   ! tracers trends type 'TRA'/'TRC'
      INTEGER                         , INTENT(in)           ::   ktra    ! tracer index
      INTEGER                         , INTENT(in)           ::   ktrd    ! tracer trend index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)           ::   ptrd    ! tracer trend  or flux
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::   pun     ! now velocity 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::   ptra    ! now tracer variable
      !
      INTEGER ::   jk    ! loop indices
      INTEGER ::   i01   ! 0 or 1
      REAL(wp),        DIMENSION(jpi,jpj,jpk) ::   ztrds             ! 3D workspace
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zwt, zws, ztrdt   ! 3D workspace
      !!----------------------------------------------------------------------
      !      
      IF( .NOT. ALLOCATED( trdtx ) ) THEN      ! allocate trdtra arrays
         IF( trd_tra_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trd_tra : unable to allocate arrays' )
      ENDIF
      !
      i01 = COUNT( (/ PRESENT(pun) .OR. ( ktrd /= jptra_xad .AND. ktrd /= jptra_yad .AND. ktrd /= jptra_zad ) /) )
      !
      IF( ctype == 'TRA' .AND. ktra == jp_tem ) THEN   !==  Temperature trend  ==!
         !
         SELECT CASE( ktrd*i01 )
         !                            ! advection: transform the advective flux into a trend
         CASE( jptra_xad )   ;   CALL trd_tra_adv( ptrd, pun, ptra, 'X', trdtx ) 
         CASE( jptra_yad )   ;   CALL trd_tra_adv( ptrd, pun, ptra, 'Y', trdty ) 
         CASE( jptra_zad )   ;   CALL trd_tra_adv( ptrd, pun, ptra, 'Z', trdt  ) 
         CASE( jptra_bbc,    &        ! qsr, bbc: on temperature only, send to trd_tra_mng
            &  jptra_qsr )   ;   trdt(:,:,:) = ptrd(:,:,:) * tmask(:,:,:)
                                 ztrds(:,:,:) = 0._wp
                                 CALL trd_tra_mng( trdt, ztrds, ktrd, kt )
 !!gm Gurvan, verify the jptra_evd trend please !
         CASE( jptra_evd )   ;   avt_evd(:,:,:) = ptrd(:,:,:) * tmask(:,:,:)
         CASE DEFAULT                 ! other trends: masked trends
            trdt(:,:,:) = ptrd(:,:,:) * tmask(:,:,:)              ! mask & store
         END SELECT
         !
      ENDIF

      IF( ctype == 'TRA' .AND. ktra == jp_sal ) THEN      !==  Salinity trends  ==!
         !
         SELECT CASE( ktrd*i01 )
         !                            ! advection: transform the advective flux into a trend
         !                            !            and send T & S trends to trd_tra_mng
         CASE( jptra_xad  )   ;   CALL trd_tra_adv( ptrd , pun  , ptra, 'X'  , ztrds ) 
                                  CALL trd_tra_mng( trdtx, ztrds, ktrd, kt   )
         CASE( jptra_yad  )   ;   CALL trd_tra_adv( ptrd , pun  , ptra, 'Y'  , ztrds ) 
                                  CALL trd_tra_mng( trdty, ztrds, ktrd, kt   )
         CASE( jptra_zad  )   ;   CALL trd_tra_adv( ptrd , pun  , ptra, 'Z'  , ztrds ) 
                                  CALL trd_tra_mng( trdt , ztrds, ktrd, kt   )
         CASE( jptra_zdfp )           ! diagnose the "PURE" Kz trend (here: just before the swap)
            !                         ! iso-neutral diffusion case otherwise jptra_zdf is "PURE"
            ALLOCATE( zwt(jpi,jpj,jpk), zws(jpi,jpj,jpk), ztrdt(jpi,jpj,jpk) )
            !
            zwt(:,:, 1 ) = 0._wp   ;   zws(:,:, 1 ) = 0._wp            ! vertical diffusive fluxes
            zwt(:,:,jpk) = 0._wp   ;   zws(:,:,jpk) = 0._wp
            DO jk = 2, jpk
               zwt(:,:,jk) = avt(:,:,jk) * ( tsa(:,:,jk-1,jp_tem) - tsa(:,:,jk,jp_tem) ) / e3w_n(:,:,jk) * tmask(:,:,jk)
               zws(:,:,jk) = avs(:,:,jk) * ( tsa(:,:,jk-1,jp_sal) - tsa(:,:,jk,jp_sal) ) / e3w_n(:,:,jk) * tmask(:,:,jk)
            END DO
            !
            ztrdt(:,:,jpk) = 0._wp   ;   ztrds(:,:,jpk) = 0._wp
            DO jk = 1, jpkm1
               ztrdt(:,:,jk) = ( zwt(:,:,jk) - zwt(:,:,jk+1) ) / e3t_n(:,:,jk)
               ztrds(:,:,jk) = ( zws(:,:,jk) - zws(:,:,jk+1) ) / e3t_n(:,:,jk) 
            END DO
            CALL trd_tra_mng( ztrdt, ztrds, jptra_zdfp, kt )  
            !
            !                         ! Also calculate EVD trend at this point. 
            zwt(:,:,:) = 0._wp   ;   zws(:,:,:) = 0._wp            ! vertical diffusive fluxes
            DO jk = 2, jpk
               zwt(:,:,jk) = avt_evd(:,:,jk) * ( tsa(:,:,jk-1,jp_tem) - tsa(:,:,jk,jp_tem) ) / e3w_n(:,:,jk) * tmask(:,:,jk)
               zws(:,:,jk) = avt_evd(:,:,jk) * ( tsa(:,:,jk-1,jp_sal) - tsa(:,:,jk,jp_sal) ) / e3w_n(:,:,jk) * tmask(:,:,jk)
            END DO
            !
            ztrdt(:,:,jpk) = 0._wp   ;   ztrds(:,:,jpk) = 0._wp
            DO jk = 1, jpkm1
               ztrdt(:,:,jk) = ( zwt(:,:,jk) - zwt(:,:,jk+1) ) / e3t_n(:,:,jk)
               ztrds(:,:,jk) = ( zws(:,:,jk) - zws(:,:,jk+1) ) / e3t_n(:,:,jk) 
            END DO
            CALL trd_tra_mng( ztrdt, ztrds, jptra_evd, kt )  
            !
            DEALLOCATE( zwt, zws, ztrdt )
            !
         CASE DEFAULT                 ! other trends: mask and send T & S trends to trd_tra_mng
            ztrds(:,:,:) = ptrd(:,:,:) * tmask(:,:,:)
            CALL trd_tra_mng( trdt, ztrds, ktrd, kt )  
         END SELECT
      ENDIF

      IF( ctype == 'TRC' ) THEN                           !==  passive tracer trend  ==!
         !
         SELECT CASE( ktrd*i01 )
         !                            ! advection: transform the advective flux into a masked trend
         CASE( jptra_xad )   ;   CALL trd_tra_adv( ptrd , pun , ptra, 'X', ztrds ) 
         CASE( jptra_yad )   ;   CALL trd_tra_adv( ptrd , pun , ptra, 'Y', ztrds ) 
         CASE( jptra_zad )   ;   CALL trd_tra_adv( ptrd , pun , ptra, 'Z', ztrds ) 
         CASE DEFAULT                 ! other trends: just masked 
                                 ztrds(:,:,:) = ptrd(:,:,:) * tmask(:,:,:)
         END SELECT
         !                            ! send trend to trd_trc
         CALL trd_trc( ztrds, ktra, ktrd, kt ) 
         !
      ENDIF
      !
   END SUBROUTINE trd_tra


   SUBROUTINE trd_tra_adv( pf, pun, ptn, cdir, ptrd )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra_adv  ***
      !! 
      !! ** Purpose :   transformed a advective flux into a masked advective trends
      !!
      !! ** Method  :   use the following transformation: -div(U.T) = - U grad(T) + T.div(U)
      !!       i-advective trends = -un. di-1[T] = -( di-1[fi] - tn di-1[un] )
      !!       j-advective trends = -un. di-1[T] = -( dj-1[fi] - tn dj-1[un] )
      !!       k-advective trends = -un. di+1[T] = -( dk+1[fi] - tn dk+1[un] )
      !!                where fi is the incoming advective flux.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pf      ! advective flux in one direction
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pun     ! now velocity   in one direction
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   ptn     ! now or before tracer 
      CHARACTER(len=1)                , INTENT(in   ) ::   cdir    ! X/Y/Z direction
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(  out) ::   ptrd    ! advective trend in one direction
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ii, ij, ik   ! index shift as function of the direction
      !!----------------------------------------------------------------------
      !
      SELECT CASE( cdir )      ! shift depending on the direction
      CASE( 'X' )   ;   ii = 1   ;   ij = 0   ;   ik = 0      ! i-trend
      CASE( 'Y' )   ;   ii = 0   ;   ij = 1   ;   ik = 0      ! j-trend
      CASE( 'Z' )   ;   ii = 0   ;   ij = 0   ;   ik =-1      ! k-trend
      END SELECT
      !
      !                        ! set to zero uncomputed values
      ptrd(jpi,:,:) = 0._wp   ;   ptrd(1,:,:) = 0._wp
      ptrd(:,jpj,:) = 0._wp   ;   ptrd(:,1,:) = 0._wp
      ptrd(:,:,jpk) = 0._wp
      !
      DO jk = 1, jpkm1         ! advective trend
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ptrd(ji,jj,jk) = - (     pf (ji,jj,jk) - pf (ji-ii,jj-ij,jk-ik)                        &
                 &                  - ( pun(ji,jj,jk) - pun(ji-ii,jj-ij,jk-ik) ) * ptn(ji,jj,jk)  )   &
                 &              * r1_e1e2t(ji,jj) / e3t_n(ji,jj,jk) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO
      !
   END SUBROUTINE trd_tra_adv


   SUBROUTINE trd_tra_mng( ptrdx, ptrdy, ktrd, kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra_mng  ***
      !! 
      !! ** Purpose :   Dispatch all tracer trends computation, e.g. 3D output,
      !!                integral constraints, potential energy, and/or 
      !!                mixed layer budget.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptrdx   ! Temperature or U trend 
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptrdy   ! Salinity    or V trend
      INTEGER                   , INTENT(in   ) ::   ktrd    ! tracer trend index
      INTEGER                   , INTENT(in   ) ::   kt      ! time step
      !!----------------------------------------------------------------------

      IF( neuler == 0 .AND. kt == nit000    ) THEN   ;   r2dt =      rdt      ! = rdt (restart with Euler time stepping)
      ELSEIF(               kt <= nit000 + 1) THEN   ;   r2dt = 2. * rdt      ! = 2 rdt (leapfrog)
      ENDIF

      !                   ! 3D output of tracers trends using IOM interface
      IF( ln_tra_trd )   CALL trd_tra_iom ( ptrdx, ptrdy, ktrd, kt )

      !                   ! Integral Constraints Properties for tracers trends                                       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_glo_trd )   CALL trd_glo( ptrdx, ptrdy, ktrd, 'TRA', kt )

      !                   ! Potential ENergy trends
      IF( ln_PE_trd  )   CALL trd_pen( ptrdx, ptrdy, ktrd, kt, r2dt )

      !                   ! Mixed layer trends for active tracers
      IF( ln_tra_mxl )   THEN   
         !-----------------------------------------------------------------------------------------------
         ! W.A.R.N.I.N.G :
         ! jptra_ldf : called by traldf.F90
         !                 at this stage we store:
         !                  - the lateral geopotential diffusion (here, lateral = horizontal)
         !                  - and the iso-neutral diffusion if activated 
         ! jptra_zdf : called by trazdf.F90
         !                 * in case of iso-neutral diffusion we store the vertical diffusion component in the 
         !                   lateral trend including the K_z contrib, which will be removed later (see trd_mxl)
         !-----------------------------------------------------------------------------------------------

         SELECT CASE ( ktrd )
         CASE ( jptra_xad )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_xad, '3D' )   ! zonal    advection
         CASE ( jptra_yad )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_yad, '3D' )   ! merid.   advection
         CASE ( jptra_zad )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_zad, '3D' )   ! vertical advection
         CASE ( jptra_ldf )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_ldf, '3D' )   ! lateral  diffusion
         CASE ( jptra_bbl )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_bbl, '3D' )   ! bottom boundary layer
         CASE ( jptra_zdf )
            IF( ln_traldf_iso ) THEN ; CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_ldf, '3D' )   ! lateral  diffusion (K_z)
            ELSE                   ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_zdf, '3D' )   ! vertical diffusion (K_z)
            ENDIF
         CASE ( jptra_dmp )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_dmp, '3D' )   ! internal 3D restoring (tradmp)
         CASE ( jptra_qsr )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_for, '3D' )   ! air-sea : penetrative sol radiat
         CASE ( jptra_nsr )        ;   ptrdx(:,:,2:jpk) = 0._wp   ;   ptrdy(:,:,2:jpk) = 0._wp
                                       CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_for, '2D' )   ! air-sea : non penetr sol radiation
         CASE ( jptra_bbc )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_bbc, '3D' )   ! bottom bound cond (geoth flux)
         CASE ( jptra_npc )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_npc, '3D' )   ! non penetr convect adjustment
         CASE ( jptra_atf )        ;   CALL trd_mxl_zint( ptrdx, ptrdy, jpmxl_atf, '3D' )   ! asselin time filter (last trend)
                                   !
                                       CALL trd_mxl( kt, r2dt )                             ! trends: Mixed-layer (output)
         END SELECT
         !
      ENDIF
      !
   END SUBROUTINE trd_tra_mng


   SUBROUTINE trd_tra_iom( ptrdx, ptrdy, ktrd, kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_tra_iom  ***
      !! 
      !! ** Purpose :   output 3D tracer trends using IOM
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptrdx   ! Temperature or U trend 
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptrdy   ! Salinity    or V trend
      INTEGER                   , INTENT(in   ) ::   ktrd    ! tracer trend index
      INTEGER                   , INTENT(in   ) ::   kt      ! time step
      !!
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      INTEGER ::   ikbu, ikbv   ! local integers
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   z2dx, z2dy   ! 2D workspace 
      !!----------------------------------------------------------------------
      !
!!gm Rq: mask the trends already masked in trd_tra, but lbc_lnk should probably be added
      !
      ! Trends evaluated every time step that could go to the standard T file and can be output every ts into a 1ts file if 1ts output is selected
      SELECT CASE( ktrd )
      ! This total trend is done every time step
      CASE( jptra_tot  )   ;   CALL iom_put( "ttrd_tot" , ptrdx )           ! model total trend
                               CALL iom_put( "strd_tot" , ptrdy )
      END SELECT
      !
      ! These trends are done every second time step. When 1ts output is selected must go different (2ts) file from standard T-file
      IF( MOD( kt, 2 ) == 0 ) THEN
         SELECT CASE( ktrd )
         CASE( jptra_xad  )   ;   CALL iom_put( "ttrd_xad"  , ptrdx )        ! x- horizontal advection
                                  CALL iom_put( "strd_xad"  , ptrdy )
         CASE( jptra_yad  )   ;   CALL iom_put( "ttrd_yad"  , ptrdx )        ! y- horizontal advection
                                  CALL iom_put( "strd_yad"  , ptrdy )
         CASE( jptra_zad  )   ;   CALL iom_put( "ttrd_zad"  , ptrdx )        ! z- vertical   advection
                                  CALL iom_put( "strd_zad"  , ptrdy )
                                  IF( ln_linssh ) THEN                   ! cst volume : adv flux through z=0 surface
                                     ALLOCATE( z2dx(jpi,jpj), z2dy(jpi,jpj) )
                                     z2dx(:,:) = wn(:,:,1) * tsn(:,:,1,jp_tem) / e3t_n(:,:,1)
                                     z2dy(:,:) = wn(:,:,1) * tsn(:,:,1,jp_sal) / e3t_n(:,:,1)
                                     CALL iom_put( "ttrd_sad", z2dx )
                                     CALL iom_put( "strd_sad", z2dy )
                                     DEALLOCATE( z2dx, z2dy )
                                  ENDIF
         CASE( jptra_totad  ) ;   CALL iom_put( "ttrd_totad", ptrdx )        ! total   advection
                                  CALL iom_put( "strd_totad", ptrdy )
         CASE( jptra_ldf  )   ;   CALL iom_put( "ttrd_ldf"  , ptrdx )        ! lateral diffusion
                                  CALL iom_put( "strd_ldf"  , ptrdy )
         CASE( jptra_zdf  )   ;   CALL iom_put( "ttrd_zdf"  , ptrdx )        ! vertical diffusion (including Kz contribution)
                                  CALL iom_put( "strd_zdf"  , ptrdy )
         CASE( jptra_zdfp )   ;   CALL iom_put( "ttrd_zdfp" , ptrdx )        ! PURE vertical diffusion (no isoneutral contribution)
                                  CALL iom_put( "strd_zdfp" , ptrdy )
         CASE( jptra_evd )    ;   CALL iom_put( "ttrd_evd"  , ptrdx )        ! EVD trend (convection)
                                  CALL iom_put( "strd_evd"  , ptrdy )
         CASE( jptra_dmp  )   ;   CALL iom_put( "ttrd_dmp"  , ptrdx )        ! internal restoring (damping)
                                  CALL iom_put( "strd_dmp"  , ptrdy )
         CASE( jptra_bbl  )   ;   CALL iom_put( "ttrd_bbl"  , ptrdx )        ! bottom boundary layer
                                  CALL iom_put( "strd_bbl"  , ptrdy )
         CASE( jptra_npc  )   ;   CALL iom_put( "ttrd_npc"  , ptrdx )        ! static instability mixing
                                  CALL iom_put( "strd_npc"  , ptrdy )
         CASE( jptra_bbc  )   ;   CALL iom_put( "ttrd_bbc"  , ptrdx )        ! geothermal heating   (only on temperature)
         CASE( jptra_nsr  )   ;   CALL iom_put( "ttrd_qns"  , ptrdx(:,:,1) ) ! surface forcing + runoff (ln_rnf=T)
                                  CALL iom_put( "strd_cdt"  , ptrdy(:,:,1) )        ! output as 2D surface fields
         CASE( jptra_qsr  )   ;   CALL iom_put( "ttrd_qsr"  , ptrdx )        ! penetrative solar radiat. (only on temperature)
         END SELECT
         ! the Asselin filter trend  is also every other time step but needs to be lagged one time step
         ! Even when 1ts output is selected can go to the same (2ts) file as the trends plotted every even time step.
      ELSE IF( MOD( kt, 2 ) == 1 ) THEN
         SELECT CASE( ktrd )
         CASE( jptra_atf  )   ;   CALL iom_put( "ttrd_atf" , ptrdx )        ! asselin time Filter
                                  CALL iom_put( "strd_atf" , ptrdy )
         END SELECT
      END IF
      !
   END SUBROUTINE trd_tra_iom

   !!======================================================================
END MODULE trdtra
