MODULE icbdyn
   !!======================================================================
   !!                       ***  MODULE  icbdyn  ***
   !! Iceberg:  time stepping routine for iceberg tracking
   !!======================================================================
   !! History :  3.3  !  2010-01  (Martin&Adcroft)  Original code
   !!             -   !  2011-03  (Madec)  Part conversion to NEMO form
   !!             -   !                    Removal of mapping from another grid
   !!             -   !  2011-04  (Alderson)  Split into separate modules
   !!             -   !  2011-05  (Alderson)  Replace broken grounding routine with one of
   !!             -   !                       Gurvan's suggestions (just like the broken one)
   !!----------------------------------------------------------------------
   USE par_oce        ! NEMO parameters
   USE dom_oce        ! NEMO ocean domain
   USE phycst         ! NEMO physical constants
   !
   USE icb_oce        ! define iceberg arrays
   USE icbutl         ! iceberg utility routines
   USE icbdia         ! iceberg budget routines

   IMPLICIT NONE
   PRIVATE

   PUBLIC   icb_dyn  ! routine called in icbstp.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE icb_dyn( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_dyn  ***
      !!
      !! ** Purpose :   iceberg evolution.
      !!
      !! ** Method  : - See Martin & Adcroft, Ocean Modelling 34, 2010
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   !
      !
      LOGICAL  ::   ll_bounced
      REAL(wp) ::   zuvel1 , zvvel1 , zu1, zv1, zax1, zay1, zxi1 , zyj1
      REAL(wp) ::   zuvel2 , zvvel2 , zu2, zv2, zax2, zay2, zxi2 , zyj2
      REAL(wp) ::   zuvel3 , zvvel3 , zu3, zv3, zax3, zay3, zxi3 , zyj3
      REAL(wp) ::   zuvel4 , zvvel4 , zu4, zv4, zax4, zay4, zxi4 , zyj4
      REAL(wp) ::   zuvel_n, zvvel_n, zxi_n   , zyj_n
      REAL(wp) ::   zdt, zdt_2, zdt_6, ze1, ze2
      TYPE(iceberg), POINTER ::   berg
      TYPE(point)  , POINTER ::   pt
      !!----------------------------------------------------------------------
      !
      ! 4th order Runge-Kutta to solve:   d/dt X = V,  d/dt V = A
      !                    with I.C.'s:   X=X1 and V=V1
      !
      !                                    ; A1=A(X1,V1)
      !  X2 = X1+dt/2*V1 ; V2 = V1+dt/2*A1 ; A2=A(X2,V2)
      !  X3 = X1+dt/2*V2 ; V3 = V1+dt/2*A2 ; A3=A(X3,V3)
      !  X4 = X1+  dt*V3 ; V4 = V1+  dt*A3 ; A4=A(X4,V4)
      !
      !  Xn = X1+dt*(V1+2*V2+2*V3+V4)/6
      !  Vn = V1+dt*(A1+2*A2+2*A3+A4)/6

      ! time steps
      zdt   = berg_dt
      zdt_2 = zdt * 0.5_wp
      zdt_6 = zdt / 6._wp

      berg => first_berg                    ! start from the first berg
      !
      DO WHILE ( ASSOCIATED(berg) )          !==  loop over all bergs  ==!
         !
         pt => berg%current_point

         ll_bounced = .FALSE.


         ! STEP 1 !
         ! ====== !
         zxi1 = pt%xi   ;   zuvel1 = pt%uvel     !**   X1 in (i,j)  ;  V1 in m/s
         zyj1 = pt%yj   ;   zvvel1 = pt%vvel


         !                                         !**   A1 = A(X1,V1)
         CALL icb_accel( berg , zxi1, ze1, zuvel1, zuvel1, zax1,     &
            &                   zyj1, ze2, zvvel1, zvvel1, zay1, zdt_2 )
         !
         zu1 = zuvel1 / ze1                           !**   V1 in d(i,j)/dt
         zv1 = zvvel1 / ze2

         ! STEP 2 !
         ! ====== !
         !                                         !**   X2 = X1+dt/2*V1   ;   V2 = V1+dt/2*A1
         ! position using di/dt & djdt   !   V2  in m/s
         zxi2 = zxi1 + zdt_2 * zu1          ;   zuvel2 = zuvel1 + zdt_2 * zax1
         zyj2 = zyj1 + zdt_2 * zv1          ;   zvvel2 = zvvel1 + zdt_2 * zay1
         !
         CALL icb_ground( zxi2, zxi1, zu1,   &
            &             zyj2, zyj1, zv1, ll_bounced )

         !                                         !**   A2 = A(X2,V2)
         CALL icb_accel( berg , zxi2, ze1, zuvel2, zuvel1, zax2,    &
            &                   zyj2, ze2, zvvel2, zvvel1, zay2, zdt_2 )
         !
         zu2 = zuvel2 / ze1                           !**   V2 in d(i,j)/dt
         zv2 = zvvel2 / ze2
         !
         ! STEP 3 !
         ! ====== !
         !                                         !**  X3 = X1+dt/2*V2  ;   V3 = V1+dt/2*A2; A3=A(X3)
         zxi3  = zxi1  + zdt_2 * zu2   ;   zuvel3 = zuvel1 + zdt_2 * zax2
         zyj3  = zyj1  + zdt_2 * zv2   ;   zvvel3 = zvvel1 + zdt_2 * zay2
         !
         CALL icb_ground( zxi3, zxi1, zu3,   &
            &                zyj3, zyj1, zv3, ll_bounced )

         !                                         !**   A3 = A(X3,V3)
         CALL icb_accel( berg , zxi3, ze1, zuvel3, zuvel1, zax3,    &
            &                   zyj3, ze2, zvvel3, zvvel1, zay3, zdt )
         !
         zu3 = zuvel3 / ze1                           !**   V3 in d(i,j)/dt
         zv3 = zvvel3 / ze2

         ! STEP 4 !
         ! ====== !
         !                                         !**   X4 = X1+dt*V3   ;   V4 = V1+dt*A3
         zxi4 = zxi1 + zdt * zu3   ;   zuvel4 = zuvel1 + zdt * zax3
         zyj4 = zyj1 + zdt * zv3   ;   zvvel4 = zvvel1 + zdt * zay3

         CALL icb_ground( zxi4, zxi1, zu4,   &
            &             zyj4, zyj1, zv4, ll_bounced )

         !                                         !**   A4 = A(X4,V4)
         CALL icb_accel( berg , zxi4, ze1, zuvel4, zuvel1, zax4,    &
            &                   zyj4, ze2, zvvel4, zvvel1, zay4, zdt )

         zu4 = zuvel4 / ze1                           !**   V4 in d(i,j)/dt
         zv4 = zvvel4 / ze2

         ! FINAL STEP !
         ! ========== !
         !                                         !**   Xn = X1+dt*(V1+2*V2+2*V3+V4)/6
         !                                         !**   Vn = V1+dt*(A1+2*A2+2*A3+A4)/6
         zxi_n   = pt%xi   + zdt_6 * (  zu1  + 2.*(zu2  + zu3 ) + zu4  )
         zyj_n   = pt%yj   + zdt_6 * (  zv1  + 2.*(zv2  + zv3 ) + zv4  )
         zuvel_n = pt%uvel + zdt_6 * (  zax1 + 2.*(zax2 + zax3) + zax4 )
         zvvel_n = pt%vvel + zdt_6 * (  zay1 + 2.*(zay2 + zay3) + zay4 )

         CALL icb_ground( zxi_n, zxi1, zuvel_n,   &
            &             zyj_n, zyj1, zvvel_n, ll_bounced )

         pt%uvel = zuvel_n                        !** save in berg structure
         pt%vvel = zvvel_n
         pt%xi   = zxi_n
         pt%yj   = zyj_n

         ! update actual position
         pt%lon  = icb_utl_bilin_x(glamt, pt%xi, pt%yj )
         pt%lat  = icb_utl_bilin(gphit, pt%xi, pt%yj, 'T' )

         berg => berg%next                         ! switch to the next berg
         !
      END DO                                  !==  end loop over all bergs  ==!
      !
   END SUBROUTINE icb_dyn


   SUBROUTINE icb_ground( pi, pi0, pu,   &
      &                   pj, pj0, pv, ld_bounced )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_ground  ***
      !!
      !! ** Purpose :   iceberg grounding.
      !!
      !! ** Method  : - adjust velocity and then put iceberg back to start position
      !!                NB two possibilities available one of which is hard-coded here
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(inout) ::   pi , pj      ! current iceberg position
      REAL(wp), INTENT(in   ) ::   pi0, pj0     ! previous iceberg position
      REAL(wp), INTENT(inout) ::   pu  , pv     ! current iceberg velocities
      LOGICAL , INTENT(  out) ::   ld_bounced   ! bounced indicator
      !
      INTEGER  ::   ii, ii0
      INTEGER  ::   ij, ij0
      INTEGER  ::   ibounce_method
      !!----------------------------------------------------------------------
      !
      ld_bounced = .FALSE.
      !
      ii0 = INT( pi0+0.5 )   ;   ij0 = INT( pj0+0.5 )       ! initial gridpoint position (T-cell)
      ii  = INT( pi +0.5 )   ;   ij  = INT( pj +0.5 )       ! current     -         -
      !
      IF( ii == ii0  .AND.  ij == ij0  )   RETURN           ! berg remains in the same cell
      !
      ! map into current processor
      ii0 = mi1( ii0 )
      ij0 = mj1( ij0 )
      ii  = mi1( ii  )
      ij  = mj1( ij  )
      !
      IF(  tmask(ii,ij,1)  /=   0._wp  )   RETURN           ! berg reach a new t-cell, but an ocean one
      !
      ! From here, berg have reach land: treat grounding/bouncing
      ! -------------------------------
      ld_bounced = .TRUE.

      !! not obvious what should happen now
      !! if berg tries to enter a land box, the only location we can return it to is the start 
      !! position (pi0,pj0), since it has to be in a wet box to do any melting;
      !! first option is simply to set whole velocity to zero and move back to start point
      !! second option (suggested by gm) is only to set the velocity component in the (i,j) direction
      !! of travel to zero; at a coastal boundary this has the effect of sliding the berg along the coast

      ibounce_method = 2
      SELECT CASE ( ibounce_method )
      CASE ( 1 )
         pi = pi0
         pj = pj0
         pu = 0._wp
         pv = 0._wp
      CASE ( 2 )
         IF( ii0 /= ii ) THEN
            pi = pi0                   ! return back to the initial position
            pu = 0._wp                 ! zeroing of velocity in the direction of the grounding
         ENDIF
         IF( ij0 /= ij ) THEN
            pj = pj0                   ! return back to the initial position
            pv = 0._wp                 ! zeroing of velocity in the direction of the grounding
         ENDIF
      END SELECT
      !
   END SUBROUTINE icb_ground


   SUBROUTINE icb_accel( berg , pxi, pe1, puvel, puvel0, pax,                &
      &                         pyj, pe2, pvvel, pvvel0, pay, pdt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_accel  ***
      !!
      !! ** Purpose :   compute the iceberg acceleration.
      !!
      !! ** Method  : - sum the terms in the momentum budget
      !!----------------------------------------------------------------------
      TYPE(iceberg ), POINTER, INTENT(in   ) ::   berg             ! berg
      REAL(wp)               , INTENT(in   ) ::   pxi   , pyj      ! berg position in (i,j) referential
      REAL(wp)               , INTENT(in   ) ::   puvel , pvvel    ! berg velocity [m/s]
      REAL(wp)               , INTENT(in   ) ::   puvel0, pvvel0   ! initial berg velocity [m/s]
      REAL(wp)               , INTENT(  out) ::   pe1, pe2         ! horizontal scale factor at (xi,yj)
      REAL(wp)               , INTENT(inout) ::   pax, pay         ! berg acceleration
      REAL(wp)               , INTENT(in   ) ::   pdt              ! berg time step
      !
      REAL(wp), PARAMETER ::   pp_alpha     = 0._wp      !
      REAL(wp), PARAMETER ::   pp_beta      = 1._wp      !
      REAL(wp), PARAMETER ::   pp_vel_lim   =15._wp      ! max allowed berg speed
      REAL(wp), PARAMETER ::   pp_accel_lim = 1.e-2_wp   ! max allowed berg acceleration
      REAL(wp), PARAMETER ::   pp_Cr0       = 0.06_wp    !
      !
      INTEGER  ::   itloop
      REAL(wp) ::   zuo, zui, zua, zuwave, zssh_x, zsst, zcn, zhi, zsss
      REAL(wp) ::   zvo, zvi, zva, zvwave, zssh_y
      REAL(wp) ::   zff, zT, zD, zW, zL, zM, zF
      REAL(wp) ::   zdrag_ocn, zdrag_atm, zdrag_ice, zwave_rad
      REAL(wp) ::   z_ocn, z_atm, z_ice
      REAL(wp) ::   zampl, zwmod, zCr, zLwavelength, zLcutoff, zLtop
      REAL(wp) ::   zlambda, zdetA, zA11, zA12, zaxe, zaye, zD_hi
      REAL(wp) ::   zuveln, zvveln, zus, zvs, zspeed, zloc_dx, zspeed_new
      !!----------------------------------------------------------------------

      ! Interpolate gridded fields to berg
      nknberg = berg%number(1)
      CALL icb_utl_interp( pxi, pe1, zuo, zui, zua, zssh_x,                     &
         &                 pyj, pe2, zvo, zvi, zva, zssh_y, zsst, zcn, zhi, zff, zsss )

      zM = berg%current_point%mass
      zT = berg%current_point%thickness               ! total thickness
      zD = ( rn_rho_bergs / pp_rho_seawater ) * zT    ! draught (keel depth)
      zF = zT - zD                                    ! freeboard
      zW = berg%current_point%width
      zL = berg%current_point%length

      zhi   = MIN( zhi   , zD    )
      zD_hi = MAX( 0._wp, zD-zhi )

      ! Wave radiation
      zuwave = zua - zuo   ;   zvwave = zva - zvo     ! Use wind speed rel. to ocean for wave model
      zwmod  = zuwave*zuwave + zvwave*zvwave          ! The wave amplitude and length depend on the  current;
      !                                               ! wind speed relative to the ocean. Actually wmod is wmod**2 here.
      zampl        = 0.5 * 0.02025 * zwmod            ! This is "a", the wave amplitude
      zLwavelength =       0.32    * zwmod            ! Surface wave length fitted to data in table at
      !                                               ! http://www4.ncsu.edu/eos/users/c/ceknowle/public/chapter10/part2.html
      zLcutoff     = 0.125 * zLwavelength
      zLtop        = 0.25  * zLwavelength
      zCr          = pp_Cr0 * MIN(  MAX( 0., (zL-zLcutoff) / ((zLtop-zLcutoff)+1.e-30)) , 1.)  ! Wave radiation coefficient
      !                                               ! fitted to graph from Carrieres et al.,  POAC Drift Model.
      zwave_rad    = 0.5 * pp_rho_seawater / zM * zCr * grav * zampl * MIN( zampl,zF ) * (2.*zW*zL) / (zW+zL)
      zwmod        = SQRT( zua*zua + zva*zva )        ! Wind speed
      IF( zwmod /= 0._wp ) THEN
         zuwave = zua/zwmod   ! Wave radiation force acts in wind direction ...       !!gm  this should be the wind rel. to ocean ?
         zvwave = zva/zwmod
      ELSE
         zuwave = 0.   ;    zvwave=0.   ;    zwave_rad=0. ! ... and only when wind is present.     !!gm  wave_rad=0. is useless
      ENDIF

      ! Weighted drag coefficients
      z_ocn = pp_rho_seawater / zM * (0.5*pp_Cd_wv*zW*(zD_hi)+pp_Cd_wh*zW*zL)
      z_atm = pp_rho_air      / zM * (0.5*pp_Cd_av*zW*zF     +pp_Cd_ah*zW*zL)
      z_ice = pp_rho_ice      / zM * (0.5*pp_Cd_iv*zW*zhi              )
      IF( abs(zui) + abs(zvi) == 0._wp )   z_ice = 0._wp

      zuveln = puvel   ;   zvveln = pvvel ! Copy starting uvel, vvel
      !
      DO itloop = 1, 2  ! Iterate on drag coefficients
         !
         zus = 0.5 * ( zuveln + puvel )
         zvs = 0.5 * ( zvveln + pvvel )
         zdrag_ocn = z_ocn * SQRT( (zus-zuo)*(zus-zuo) + (zvs-zvo)*(zvs-zvo) )
         zdrag_atm = z_atm * SQRT( (zus-zua)*(zus-zua) + (zvs-zva)*(zvs-zva) )
         zdrag_ice = z_ice * SQRT( (zus-zui)*(zus-zui) + (zvs-zvi)*(zvs-zvi) )
         !
         ! Explicit accelerations
         !zaxe= zff*pvvel -grav*zssh_x +zwave_rad*zuwave &
         !    -zdrag_ocn*(puvel-zuo) -zdrag_atm*(puvel-zua) -zdrag_ice*(puvel-zui)
         !zaye=-zff*puvel -grav*zssh_y +zwave_rad*zvwave &
         !    -zdrag_ocn*(pvvel-zvo) -zdrag_atm*(pvvel-zva) -zdrag_ice*(pvvel-zvi)
         zaxe = -grav * zssh_x + zwave_rad * zuwave
         zaye = -grav * zssh_y + zwave_rad * zvwave
         IF( pp_alpha > 0._wp ) THEN   ! If implicit, use time-level (n) rather than RK4 latest
            zaxe = zaxe + zff*pvvel0
            zaye = zaye - zff*puvel0
         ELSE
            zaxe = zaxe + zff*pvvel
            zaye = zaye - zff*puvel
         ENDIF
         IF( pp_beta > 0._wp ) THEN    ! If implicit, use time-level (n) rather than RK4 latest
            zaxe = zaxe - zdrag_ocn*(puvel0-zuo) - zdrag_atm*(puvel0-zua) -zdrag_ice*(puvel0-zui)
            zaye = zaye - zdrag_ocn*(pvvel0-zvo) - zdrag_atm*(pvvel0-zva) -zdrag_ice*(pvvel0-zvi)
         ELSE
            zaxe = zaxe - zdrag_ocn*(puvel -zuo) - zdrag_atm*(puvel -zua) -zdrag_ice*(puvel -zui)
            zaye = zaye - zdrag_ocn*(pvvel -zvo) - zdrag_atm*(pvvel -zva) -zdrag_ice*(pvvel -zvi)
         ENDIF

         ! Solve for implicit accelerations
         IF( pp_alpha + pp_beta > 0._wp ) THEN
            zlambda = zdrag_ocn + zdrag_atm + zdrag_ice
            zA11    = 1._wp + pp_beta *pdt*zlambda
            zA12    =         pp_alpha*pdt*zff
            zdetA   = 1._wp / ( zA11*zA11 + zA12*zA12 )
            pax     = zdetA * ( zA11*zaxe + zA12*zaye )
            pay     = zdetA * ( zA11*zaye - zA12*zaxe )
         ELSE
            pax = zaxe   ;   pay = zaye
         ENDIF

         zuveln = puvel0 + pdt*pax
         zvveln = pvvel0 + pdt*pay
         !
      END DO      ! itloop

      IF( rn_speed_limit > 0._wp ) THEN       ! Limit speed of bergs based on a CFL criteria (if asked)
         zspeed = SQRT( zuveln*zuveln + zvveln*zvveln )    ! Speed of berg
         IF( zspeed > 0._wp ) THEN
            zloc_dx = MIN( pe1, pe2 )                          ! minimum grid spacing
            zspeed_new = zloc_dx / pdt * rn_speed_limit        ! Speed limit as a factor of dx / dt
            IF( zspeed_new < zspeed ) THEN
               zuveln = zuveln * ( zspeed_new / zspeed )        ! Scale velocity to reduce speed
               zvveln = zvveln * ( zspeed_new / zspeed )        ! without changing the direction
               CALL icb_dia_speed()
            ENDIF
         ENDIF
      ENDIF
      !                                      ! check the speed and acceleration limits
      IF (nn_verbose_level > 0) THEN
         IF( ABS( zuveln ) > pp_vel_lim   .OR. ABS( zvveln ) > pp_vel_lim   )   &
            WRITE(numicb,'("pe=",i3,x,a)') narea,'Dump triggered by excessive velocity'
         IF( ABS( pax    ) > pp_accel_lim .OR. ABS( pay    ) > pp_accel_lim )   &
            WRITE(numicb,'("pe=",i3,x,a)') narea,'Dump triggered by excessive acceleration'
      ENDIF
      !
   END SUBROUTINE icb_accel

   !!======================================================================
END MODULE icbdyn
