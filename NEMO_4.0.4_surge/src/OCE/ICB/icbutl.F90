MODULE icbutl
   !!======================================================================
   !!                       ***  MODULE  icbutl  ***
   !! Icebergs:  various iceberg utility routines
   !!======================================================================
   !! History : 3.3.1 !  2010-01  (Martin&Adcroft) Original code
   !!            -    !  2011-03  (Madec)          Part conversion to NEMO form
   !!            -    !                            Removal of mapping from another grid
   !!            -    !  2011-04  (Alderson)       Split into separate modules
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   icb_utl_interp   :
   !!   icb_utl_bilin    :
   !!   icb_utl_bilin_e  :
   !!----------------------------------------------------------------------
   USE par_oce                             ! ocean parameters
   USE dom_oce                             ! ocean domain
   USE in_out_manager                      ! IO parameters
   USE lbclnk                              ! lateral boundary condition
   USE lib_mpp                             ! MPI code and lk_mpp in particular
   USE icb_oce                             ! define iceberg arrays
   USE sbc_oce                             ! ocean surface boundary conditions
#if defined key_si3
   USE ice,    ONLY: u_ice, v_ice, hm_i    ! SI3 variables
   USE icevar                              ! ice_var_sshdyn
   USE sbc_ice, ONLY: snwice_mass, snwice_mass_b
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   icb_utl_copy          ! routine called in icbstp module
   PUBLIC   icb_utl_interp        ! routine called in icbdyn, icbthm modules
   PUBLIC   icb_utl_bilin         ! routine called in icbini, icbdyn modules
   PUBLIC   icb_utl_bilin_x       ! routine called in icbdyn module
   PUBLIC   icb_utl_add           ! routine called in icbini.F90, icbclv, icblbc and icbrst modules
   PUBLIC   icb_utl_delete        ! routine called in icblbc, icbthm modules
   PUBLIC   icb_utl_destroy       ! routine called in icbstp module
   PUBLIC   icb_utl_track         ! routine not currently used, retain just in case
   PUBLIC   icb_utl_print_berg    ! routine called in icbthm module
   PUBLIC   icb_utl_print         ! routine called in icbini, icbstp module
   PUBLIC   icb_utl_count         ! routine called in icbdia, icbini, icblbc, icbrst modules
   PUBLIC   icb_utl_incr          ! routine called in icbini, icbclv modules
   PUBLIC   icb_utl_yearday       ! routine called in icbclv, icbstp module
   PUBLIC   icb_utl_mass          ! routine called in icbdia module
   PUBLIC   icb_utl_heat          ! routine called in icbdia module

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE icb_utl_copy()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_utl_copy  ***
      !!
      !! ** Purpose :   iceberg initialization.
      !!
      !! ** Method  : - blah blah
      !!----------------------------------------------------------------------
#if defined key_si3
      REAL(wp), DIMENSION(jpi,jpj) :: zssh_lead_m    !    ocean surface (ssh_m) if ice is not embedded
      !                                              !    ocean surface in leads if ice is embedded   
#endif
      ! copy nemo forcing arrays into iceberg versions with extra halo
      ! only necessary for variables not on T points
      ! and ssh which is used to calculate gradients

      uo_e(1:jpi,1:jpj) = ssu_m(:,:) * umask(:,:,1)
      vo_e(1:jpi,1:jpj) = ssv_m(:,:) * vmask(:,:,1)
      ff_e(1:jpi,1:jpj) = ff_f (:,:) 
      tt_e(1:jpi,1:jpj) = sst_m(:,:)
      ss_e(1:jpi,1:jpj) = sss_m(:,:)
      fr_e(1:jpi,1:jpj) = fr_i (:,:)
      ua_e(1:jpi,1:jpj) = utau (:,:) * umask(:,:,1) ! maybe mask useless because mask applied in sbcblk
      va_e(1:jpi,1:jpj) = vtau (:,:) * vmask(:,:,1) ! maybe mask useless because mask applied in sbcblk
      !
      CALL lbc_lnk_icb( 'icbutl', uo_e, 'U', -1._wp, 1, 1 )
      CALL lbc_lnk_icb( 'icbutl', vo_e, 'V', -1._wp, 1, 1 )
      CALL lbc_lnk_icb( 'icbutl', ff_e, 'F', +1._wp, 1, 1 )
      CALL lbc_lnk_icb( 'icbutl', ua_e, 'U', -1._wp, 1, 1 )
      CALL lbc_lnk_icb( 'icbutl', va_e, 'V', -1._wp, 1, 1 )
      CALL lbc_lnk_icb( 'icbutl', fr_e, 'T', +1._wp, 1, 1 )
      CALL lbc_lnk_icb( 'icbutl', tt_e, 'T', +1._wp, 1, 1 )
      CALL lbc_lnk_icb( 'icbutl', ss_e, 'T', +1._wp, 1, 1 )
#if defined key_si3
      hi_e(1:jpi, 1:jpj) = hm_i (:,:)  
      ui_e(1:jpi, 1:jpj) = u_ice(:,:)
      vi_e(1:jpi, 1:jpj) = v_ice(:,:)
      !      
      ! compute ssh slope using ssh_lead if embedded
      zssh_lead_m(:,:) = ice_var_sshdyn(ssh_m, snwice_mass, snwice_mass_b)
      ssh_e(1:jpi, 1:jpj) = zssh_lead_m(:,:) * tmask(:,:,1)
      !
      CALL lbc_lnk_icb( 'icbutl', hi_e , 'T', +1._wp, 1, 1 )
      CALL lbc_lnk_icb( 'icbutl', ui_e , 'U', -1._wp, 1, 1 )
      CALL lbc_lnk_icb( 'icbutl', vi_e , 'V', -1._wp, 1, 1 )
#else
      ssh_e(1:jpi, 1:jpj) = ssh_m(:,:) * tmask(:,:,1)
#endif
      CALL lbc_lnk_icb( 'icbutl', ssh_e, 'T', +1._wp, 1, 1 )
      !
   END SUBROUTINE icb_utl_copy


   SUBROUTINE icb_utl_interp( pi, pe1, puo, pui, pua, pssh_i,   &
      &                       pj, pe2, pvo, pvi, pva, pssh_j,   &
      &                       psst, pcn, phi, pff, psss        )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE icb_utl_interp  ***
      !!
      !! ** Purpose :   interpolation
      !!
      !! ** Method  : - interpolate from various ocean arrays onto iceberg position
      !!
      !!       !!gm  CAUTION here I do not care of the slip/no-slip conditions
      !!             this can be done later (not that easy to do...)
      !!             right now, U is 0 in land so that the coastal value of velocity parallel to the coast
      !!             is half the off shore value, wile the normal-to-the-coast value is zero.
      !!             This is OK as a starting point.
      !!       !!pm  HARD CODED: - rho_air now computed in sbcblk (what are the effect ?)
      !!                         - drag coefficient (should it be namelist parameter ?)
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   ) ::   pi , pj                        ! position in (i,j) referential
      REAL(wp), INTENT(  out) ::   pe1, pe2                       ! i- and j scale factors
      REAL(wp), INTENT(  out) ::   puo, pvo, pui, pvi, pua, pva   ! ocean, ice and wind speeds
      REAL(wp), INTENT(  out) ::   pssh_i, pssh_j                 ! ssh i- & j-gradients
      REAL(wp), INTENT(  out) ::   psst, pcn, phi, pff, psss      ! SST, ice concentration, ice thickness, Coriolis, SSS
      !
      REAL(wp) ::   zcd, zmod       ! local scalars
      !!----------------------------------------------------------------------

      pe1 = icb_utl_bilin_e( e1t, e1u, e1v, e1f, pi, pj )      ! scale factors
      pe2 = icb_utl_bilin_e( e2t, e2u, e2v, e2f, pi, pj )
      !
      puo  = icb_utl_bilin_h( uo_e, pi, pj, 'U', .false.  )    ! ocean velocities
      pvo  = icb_utl_bilin_h( vo_e, pi, pj, 'V', .false.  )
      psst = icb_utl_bilin_h( tt_e, pi, pj, 'T', .true.   )    ! SST
      psss = icb_utl_bilin_h( ss_e, pi, pj, 'T', .true.   )    ! SSS
      pcn  = icb_utl_bilin_h( fr_e, pi, pj, 'T', .true.   )    ! ice concentration
      pff  = icb_utl_bilin_h( ff_e, pi, pj, 'F', .false.  )    ! Coriolis parameter
      !
      pua  = icb_utl_bilin_h( ua_e, pi, pj, 'U', .true.   )    ! 10m wind
      pva  = icb_utl_bilin_h( va_e, pi, pj, 'V', .true.   )    ! here (ua,va) are stress => rough conversion from stress to speed
      zcd  = 1.22_wp * 1.5e-3_wp                               ! air density * drag coefficient 
      zmod = 1._wp / MAX(  1.e-20, SQRT(  zcd * SQRT( pua*pua + pva*pva)  )  )
      pua  = pua * zmod                                       ! note: stress module=0 necessarly implies ua=va=0
      pva  = pva * zmod

#if defined key_si3
      pui = icb_utl_bilin_h( ui_e , pi, pj, 'U', .false. )    ! sea-ice velocities
      pvi = icb_utl_bilin_h( vi_e , pi, pj, 'V', .false. )
      phi = icb_utl_bilin_h( hi_e , pi, pj, 'T', .true.  )    ! ice thickness
#else
      pui = 0._wp
      pvi = 0._wp
      phi = 0._wp
#endif

      ! Estimate SSH gradient in i- and j-direction (centred evaluation)
      pssh_i = ( icb_utl_bilin_h( ssh_e, pi+0.1_wp, pj, 'T', .true. ) -   &
         &       icb_utl_bilin_h( ssh_e, pi-0.1_wp, pj, 'T', .true. )  ) / ( 0.2_wp * pe1 )
      pssh_j = ( icb_utl_bilin_h( ssh_e, pi, pj+0.1_wp, 'T', .true. ) -   &
         &       icb_utl_bilin_h( ssh_e, pi, pj-0.1_wp, 'T', .true. )  ) / ( 0.2_wp * pe2 )
      !
   END SUBROUTINE icb_utl_interp


   REAL(wp) FUNCTION icb_utl_bilin_h( pfld, pi, pj, cd_type, plmask )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION icb_utl_bilin  ***
      !!
      !! ** Purpose :   bilinear interpolation at berg location depending on the grid-point type
      !!                this version deals with extra halo points
      !!
      !!       !!gm  CAUTION an optional argument should be added to handle
      !!             the slip/no-slip conditions  ==>>> to be done later
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(0:jpi+1,0:jpj+1), INTENT(in) ::   pfld      ! field to be interpolated
      REAL(wp)                            , INTENT(in) ::   pi, pj    ! targeted coordinates in (i,j) referential
      CHARACTER(len=1)                    , INTENT(in) ::   cd_type   ! type of pfld array grid-points: = T , U , V or F points
      LOGICAL                             , INTENT(in) ::   plmask    ! special treatment of mask point
      !
      INTEGER  ::   ii, ij   ! local integer
      REAL(wp) ::   zi, zj   ! local real
      REAL(wp) :: zw1, zw2, zw3, zw4
      REAL(wp), DIMENSION(4) :: zmask
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( cd_type )
      CASE ( 'T' )
         ! note that here there is no +0.5 added
         ! since we're looking for four T points containing quadrant we're in of 
         ! current T cell
         ii = MAX(0, INT( pi     ))
         ij = MAX(0, INT( pj     ))    ! T-point
         zi = pi - REAL(ii,wp)
         zj = pj - REAL(ij,wp)
      CASE ( 'U' )
         ii = MAX(0, INT( pi-0.5_wp ))
         ij = MAX(0, INT( pj     ))    ! U-point
         zi = pi - 0.5_wp - REAL(ii,wp)
         zj = pj - REAL(ij,wp)
      CASE ( 'V' )
         ii = MAX(0, INT( pi     ))
         ij = MAX(0, INT( pj-0.5_wp ))    ! V-point
         zi = pi - REAL(ii,wp)
         zj = pj - 0.5_wp - REAL(ij,wp)
      CASE ( 'F' )
         ii = MAX(0, INT( pi-0.5_wp ))
         ij = MAX(0, INT( pj-0.5_wp ))    ! F-point
         zi = pi - 0.5_wp - REAL(ii,wp)
         zj = pj - 0.5_wp - REAL(ij,wp)
      END SELECT
      !
      ! find position in this processor. Prevent near edge problems (see #1389)
      ! (PM) will be useless if extra halo is used in NEMO
      !
      IF    ( ii <= mig(1)-1 ) THEN   ;   ii = 0
      ELSEIF( ii  > mig(jpi) ) THEN   ;   ii = jpi
      ELSE                            ;   ii = mi1(ii)
      ENDIF
      IF    ( ij <= mjg(1)-1 ) THEN   ;   ij = 0
      ELSEIF( ij  > mjg(jpj) ) THEN   ;   ij = jpj
      ELSE                            ;   ij = mj1(ij)
      ENDIF
      !
      ! define mask array 
      IF (plmask) THEN
         ! land value is not used in the interpolation
         SELECT CASE ( cd_type )
         CASE ( 'T' )
            zmask = (/tmask_e(ii,ij), tmask_e(ii+1,ij), tmask_e(ii,ij+1), tmask_e(ii+1,ij+1)/)
         CASE ( 'U' )
            zmask = (/umask_e(ii,ij), umask_e(ii+1,ij), umask_e(ii,ij+1), umask_e(ii+1,ij+1)/)
         CASE ( 'V' )
            zmask = (/vmask_e(ii,ij), vmask_e(ii+1,ij), vmask_e(ii,ij+1), vmask_e(ii+1,ij+1)/)
         CASE ( 'F' )
            ! F case only used for coriolis, ff_f is not mask so zmask = 1
            zmask = 1.
         END SELECT
      ELSE
         ! land value is used during interpolation
         zmask = 1.
      END iF
      !
      ! compute weight
      zw1 = zmask(1) * (1._wp-zi) * (1._wp-zj)
      zw2 = zmask(2) *        zi  * (1._wp-zj)
      zw3 = zmask(3) * (1._wp-zi) *        zj
      zw4 = zmask(4) *        zi  *        zj
      !
      ! compute interpolated value
      icb_utl_bilin_h = ( pfld(ii,ij)*zw1 + pfld(ii+1,ij)*zw2 + pfld(ii,ij+1)*zw3 + pfld(ii+1,ij+1)*zw4 ) / MAX(1.e-20, zw1+zw2+zw3+zw4) 
      !
   END FUNCTION icb_utl_bilin_h


   REAL(wp) FUNCTION icb_utl_bilin( pfld, pi, pj, cd_type )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION icb_utl_bilin  ***
      !!
      !! ** Purpose :   bilinear interpolation at berg location depending on the grid-point type
      !!
      !!       !!gm  CAUTION an optional argument should be added to handle
      !!             the slip/no-slip conditions  ==>>> to be done later
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pfld      ! field to be interpolated
      REAL(wp)                    , INTENT(in) ::   pi, pj    ! targeted coordinates in (i,j) referential
      CHARACTER(len=1)            , INTENT(in) ::   cd_type   ! type of pfld array grid-points: = T , U , V or F points
      !
      INTEGER  ::   ii, ij   ! local integer
      REAL(wp) ::   zi, zj   ! local real
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( cd_type )
         CASE ( 'T' )
            ! note that here there is no +0.5 added
            ! since we're looking for four T points containing quadrant we're in of 
            ! current T cell
            ii = MAX(1, INT( pi     ))
            ij = MAX(1, INT( pj     ))    ! T-point
            zi = pi - REAL(ii,wp)
            zj = pj - REAL(ij,wp)
         CASE ( 'U' )
            ii = MAX(1, INT( pi-0.5 ))
            ij = MAX(1, INT( pj     ))    ! U-point
            zi = pi - 0.5 - REAL(ii,wp)
            zj = pj - REAL(ij,wp)
         CASE ( 'V' )
            ii = MAX(1, INT( pi     ))
            ij = MAX(1, INT( pj-0.5 ))    ! V-point
            zi = pi - REAL(ii,wp)
            zj = pj - 0.5 - REAL(ij,wp)
         CASE ( 'F' )
            ii = MAX(1, INT( pi-0.5 ))
            ij = MAX(1, INT( pj-0.5 ))    ! F-point
            zi = pi - 0.5 - REAL(ii,wp)
            zj = pj - 0.5 - REAL(ij,wp)
      END SELECT
      !
      ! find position in this processor. Prevent near edge problems (see #1389)
      IF    ( ii < mig( 1 ) ) THEN   ;   ii = 1
      ELSEIF( ii > mig(jpi) ) THEN   ;   ii = jpi
      ELSE                           ;   ii = mi1(ii)
      ENDIF
      IF    ( ij < mjg( 1 ) ) THEN   ;   ij = 1
      ELSEIF( ij > mjg(jpj) ) THEN   ;   ij = jpj
      ELSE                           ;   ij  = mj1(ij)
      ENDIF
      !
      IF( ii == jpi )   ii = ii-1      
      IF( ij == jpj )   ij = ij-1
      !
      icb_utl_bilin = ( pfld(ii,ij  ) * (1.-zi) + pfld(ii+1,ij  ) * zi ) * (1.-zj)   &
         &          + ( pfld(ii,ij+1) * (1.-zi) + pfld(ii+1,ij+1) * zi ) *     zj
      !
   END FUNCTION icb_utl_bilin


   REAL(wp) FUNCTION icb_utl_bilin_x( pfld, pi, pj )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION icb_utl_bilin_x  ***
      !!
      !! ** Purpose :   bilinear interpolation at berg location depending on the grid-point type
      !!                Special case for interpolating longitude
      !!
      !!       !!gm  CAUTION an optional argument should be added to handle
      !!             the slip/no-slip conditions  ==>>> to be done later
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pfld      ! field to be interpolated
      REAL(wp)                    , INTENT(in) ::   pi, pj    ! targeted coordinates in (i,j) referential
      !
      INTEGER                                  ::   ii, ij   ! local integer
      REAL(wp)                                 ::   zi, zj   ! local real
      REAL(wp)                                 ::   zret     ! local real
      REAL(wp), DIMENSION(4)                   ::   z4
      !!----------------------------------------------------------------------
      !
      ! note that here there is no +0.5 added
      ! since we're looking for four T points containing quadrant we're in of 
      ! current T cell
      ii = MAX(1, INT( pi     ))
      ij = MAX(1, INT( pj     ))    ! T-point
      zi = pi - REAL(ii,wp)
      zj = pj - REAL(ij,wp)
      !
      ! find position in this processor. Prevent near edge problems (see #1389)
      IF    ( ii < mig( 1 ) ) THEN   ;   ii = 1
      ELSEIF( ii > mig(jpi) ) THEN   ;   ii = jpi
      ELSE                           ;   ii = mi1(ii)
      ENDIF
      IF    ( ij < mjg( 1 ) ) THEN   ;   ij = 1
      ELSEIF( ij > mjg(jpj) ) THEN   ;   ij = jpj
      ELSE                           ;   ij  = mj1(ij)
      ENDIF
      !
      IF( ii == jpi )   ii = ii-1      
      IF( ij == jpj )   ij = ij-1
      !
      z4(1) = pfld(ii  ,ij  )
      z4(2) = pfld(ii+1,ij  )
      z4(3) = pfld(ii  ,ij+1)
      z4(4) = pfld(ii+1,ij+1)
      IF( MAXVAL(z4) - MINVAL(z4) > 90._wp ) THEN
         WHERE( z4 < 0._wp ) z4 = z4 + 360._wp
      ENDIF
      !
      zret = (z4(1) * (1.-zi) + z4(2) * zi) * (1.-zj) + (z4(3) * (1.-zi) + z4(4) * zi) * zj
      IF( zret > 180._wp ) zret = zret - 360._wp
      icb_utl_bilin_x = zret
      !
   END FUNCTION icb_utl_bilin_x


   REAL(wp) FUNCTION icb_utl_bilin_e( pet, peu, pev, pef, pi, pj )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION dom_init  ***
      !!
      !! ** Purpose :   bilinear interpolation at berg location of horizontal scale factor
      !! ** Method  :   interpolation done using the 4 nearest grid points among
      !!                t-, u-, v-, and f-points.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in) ::   pet, peu, pev, pef   ! horizontal scale factor to be interpolated at t-,u-,v- & f-pts
      REAL(wp)                , INTENT(in) ::   pi, pj               ! targeted coordinates in (i,j) referential
      !
      INTEGER  ::   ii, ij, icase, ierr   ! local integer
      !
      ! weights corresponding to corner points of a T cell quadrant
      REAL(wp) ::   zi, zj          ! local real
      !
      ! values at corner points of a T cell quadrant
      ! 00 = bottom left, 10 = bottom right, 01 = top left, 11 = top right
      REAL(wp) ::   ze00, ze10, ze01, ze11
      !!----------------------------------------------------------------------
      !
      ii = MAX(1, INT( pi ))   ;   ij = MAX(1, INT( pj ))            ! left bottom T-point (i,j) indices

      ! fractional box spacing
      ! 0   <= zi < 0.5, 0   <= zj < 0.5   -->  NW quadrant of current T cell
      ! 0.5 <= zi < 1  , 0   <= zj < 0.5   -->  NE quadrant
      ! 0   <= zi < 0.5, 0.5 <= zj < 1     -->  SE quadrant
      ! 0.5 <= zi < 1  , 0.5 <= zj < 1     -->  SW quadrant

      zi = pi - REAL(ii,wp)          !!gm use here mig, mjg arrays
      zj = pj - REAL(ij,wp)

      ! find position in this processor. Prevent near edge problems (see #1389)
      !
      ierr = 0
      IF    ( ii < mig( 1 ) ) THEN   ;   ii = 1       ; ierr = ierr + 1
      ELSEIF( ii > mig(jpi) ) THEN   ;   ii = jpi     ; ierr = ierr + 1
      ELSE                           ;   ii = mi1(ii)
      ENDIF
      IF    ( ij < mjg( 1 ) ) THEN   ;   ij = 1       ; ierr = ierr + 1
      ELSEIF( ij > mjg(jpj) ) THEN   ;   ij = jpj     ; ierr = ierr + 1
      ELSE                           ;   ij  = mj1(ij)
      ENDIF
      !
      IF( ii == jpi ) THEN ; ii = ii-1 ; ierr = ierr + 1 ; END IF     
      IF( ij == jpj ) THEN ; ij = ij-1 ; ierr = ierr + 1 ; END IF
      !
      IF ( ierr > 0 ) CALL ctl_stop('STOP','icb_utl_bilin_e: an icebergs coordinates is out of valid range (out of bound error)')
      !
      IF(    0.0_wp <= zi .AND. zi < 0.5_wp   ) THEN
         IF( 0.0_wp <= zj .AND. zj < 0.5_wp        )   THEN        !  NE quadrant
            !                                                      !             i=I       i=I+1/2
            ze01 = pev(ii  ,ij  )   ;   ze11 = pef(ii  ,ij  )      !   j=J+1/2    V ------- F
            ze00 = pet(ii  ,ij  )   ;   ze10 = peu(ii  ,ij  )      !   j=J        T ------- U
            zi = 2._wp * zi
            zj = 2._wp * zj
         ELSE                                                      !  SE quadrant
            !                                                                    !             i=I       i=I+1/2
            ze01 = pet(ii  ,ij+1)   ;   ze11 = peu(ii  ,ij+1)      !   j=J+1      T ------- U
            ze00 = pev(ii  ,ij  )   ;   ze10 = pef(ii  ,ij  )      !   j=J+1/2    V ------- F
            zi = 2._wp *  zi
            zj = 2._wp * (zj-0.5_wp)
         ENDIF
      ELSE
         IF( 0.0_wp <= zj .AND. zj < 0.5_wp        )   THEN        !  NW quadrant
            !                                                                    !             i=I       i=I+1/2
            ze01 = pef(ii  ,ij  )   ;   ze11 = pev(ii+1,ij)        !   j=J+1/2    F ------- V
            ze00 = peu(ii  ,ij  )   ;   ze10 = pet(ii+1,ij)        !   j=J        U ------- T
            zi = 2._wp * (zi-0.5_wp)
            zj = 2._wp *  zj
         ELSE                                                      !  SW quadrant
            !                                                                    !             i=I+1/2   i=I+1
            ze01 = peu(ii  ,ij+1)   ;   ze11 = pet(ii+1,ij+1)      !   j=J+1      U ------- T
            ze00 = pef(ii  ,ij  )   ;   ze10 = pev(ii+1,ij  )      !   j=J+1/2    F ------- V
            zi = 2._wp * (zi-0.5_wp)
            zj = 2._wp * (zj-0.5_wp)
         ENDIF
      ENDIF
      !
      icb_utl_bilin_e = ( ze01 * (1._wp-zi) + ze11 * zi ) *        zj    &
         &            + ( ze00 * (1._wp-zi) + ze10 * zi ) * (1._wp-zj)
      !
   END FUNCTION icb_utl_bilin_e


   SUBROUTINE icb_utl_add( bergvals, ptvals )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE icb_utl_add           ***
      !!
      !! ** Purpose :   add a new berg to the iceberg list
      !!
      !!----------------------------------------------------------------------
      TYPE(iceberg), INTENT(in)           ::   bergvals
      TYPE(point)  , INTENT(in)           ::   ptvals
      !
      TYPE(iceberg), POINTER ::   new => NULL()
      !!----------------------------------------------------------------------
      !
      new => NULL()
      CALL icb_utl_create( new, bergvals, ptvals )
      CALL icb_utl_insert( new )
      new => NULL()     ! Clear new
      !
   END SUBROUTINE icb_utl_add         


   SUBROUTINE icb_utl_create( berg, bergvals, ptvals )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE icb_utl_create  ***
      !!
      !! ** Purpose :   add a new berg to the iceberg list
      !!
      !!----------------------------------------------------------------------
      TYPE(iceberg), INTENT(in) ::   bergvals
      TYPE(point)  , INTENT(in) ::   ptvals
      TYPE(iceberg), POINTER    ::   berg
      !
      TYPE(point)  , POINTER    ::   pt
      INTEGER                   ::   istat
      !!----------------------------------------------------------------------
      !
      IF( ASSOCIATED(berg) )   CALL ctl_stop( 'icebergs, icb_utl_create: berg already associated' )
      ALLOCATE(berg, STAT=istat)
      IF( istat /= 0 ) CALL ctl_stop( 'failed to allocate iceberg' )
      berg%number(:) = bergvals%number(:)
      berg%mass_scaling = bergvals%mass_scaling
      berg%prev => NULL()
      berg%next => NULL()
      !
      ALLOCATE(pt, STAT=istat)
      IF( istat /= 0 ) CALL ctl_stop( 'failed to allocate first iceberg point' )
      pt = ptvals
      berg%current_point => pt
      !
   END SUBROUTINE icb_utl_create


   SUBROUTINE icb_utl_insert( newberg )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_insert  ***
      !!
      !! ** Purpose :   add a new berg to the iceberg list
      !!
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER  ::   newberg
      !
      TYPE(iceberg), POINTER  ::   this, prev, last
      !!----------------------------------------------------------------------
      !
      IF( ASSOCIATED( first_berg ) ) THEN
         last => first_berg
         DO WHILE (ASSOCIATED(last%next))
            last => last%next
         ENDDO
         newberg%prev => last
         last%next    => newberg
         last         => newberg
      ELSE                       ! list is empty so create it
         first_berg => newberg
      ENDIF
      !
   END SUBROUTINE icb_utl_insert


   REAL(wp) FUNCTION icb_utl_yearday(kmon, kday, khr, kmin, ksec)
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION icb_utl_yearday  ***
      !!
      !! ** Purpose :   
      !!
      ! sga - improved but still only applies to 365 day year, need to do this properly
      !
      !!gm  all these info are already known in daymod, no???
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)     :: kmon, kday, khr, kmin, ksec
      ! 
      INTEGER, DIMENSION(12)  :: imonths = (/ 0,31,28,31,30,31,30,31,31,30,31,30 /)
      !!----------------------------------------------------------------------
      !
      icb_utl_yearday = REAL( SUM( imonths(1:kmon) ), wp )
      icb_utl_yearday = icb_utl_yearday + REAL(kday-1,wp) + (REAL(khr,wp) + (REAL(kmin,wp) + REAL(ksec,wp)/60.)/60.)/24.
      !
   END FUNCTION icb_utl_yearday

   !!-------------------------------------------------------------------------

   SUBROUTINE icb_utl_delete( first, berg )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_delete  ***
      !!
      !! ** Purpose :   
      !!
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER :: first, berg
      !!----------------------------------------------------------------------
      ! Connect neighbors to each other
      IF ( ASSOCIATED(berg%prev) ) THEN
        berg%prev%next => berg%next
      ELSE
        first => berg%next
      ENDIF
      IF (ASSOCIATED(berg%next)) berg%next%prev => berg%prev
      !
      CALL icb_utl_destroy(berg)
      !
   END SUBROUTINE icb_utl_delete


   SUBROUTINE icb_utl_destroy( berg )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_destroy  ***
      !!
      !! ** Purpose :   remove a single iceberg instance
      !!
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER :: berg
      !!----------------------------------------------------------------------
      !
      ! Remove any points
      IF( ASSOCIATED( berg%current_point ) )   DEALLOCATE( berg%current_point )
      !
      DEALLOCATE(berg)
      !
   END SUBROUTINE icb_utl_destroy


   SUBROUTINE icb_utl_track( knum, cd_label, kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_track  ***
      !!
      !! ** Purpose :   
      !!
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(nkounts)    :: knum       ! iceberg number
      CHARACTER(len=*)               :: cd_label   ! 
      INTEGER                        :: kt         ! timestep number
      ! 
      TYPE(iceberg), POINTER         :: this
      LOGICAL                        :: match
      INTEGER                        :: k
      !!----------------------------------------------------------------------
      !
      this => first_berg
      DO WHILE( ASSOCIATED(this) )
         match = .TRUE.
         DO k = 1, nkounts
            IF( this%number(k) /= knum(k) ) match = .FALSE.
         END DO
         IF( match )   CALL icb_utl_print_berg(this, kt)
         this => this%next
      END DO
      !
   END SUBROUTINE icb_utl_track


   SUBROUTINE icb_utl_print_berg( berg, kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_print_berg  ***
      !!
      !! ** Purpose :   print one
      !!
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER :: berg
      TYPE(point)  , POINTER :: pt
      INTEGER                :: kt      ! timestep number
      !!----------------------------------------------------------------------
      !
      IF (nn_verbose_level == 0) RETURN
      pt => berg%current_point
      WRITE(numicb, 9200) kt, berg%number(1), &
                   pt%xi, pt%yj, pt%lon, pt%lat, pt%uvel, pt%vvel,  &
                   pt%uo, pt%vo, pt%ua, pt%va, pt%ui, pt%vi
      CALL flush( numicb )
 9200 FORMAT(5x,i5,2x,i10,6(2x,2f10.4))
      !
   END SUBROUTINE icb_utl_print_berg


   SUBROUTINE icb_utl_print( cd_label, kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_print  ***
      !!
      !! ** Purpose :   print many
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=*)       :: cd_label
      INTEGER                :: kt             ! timestep number
      ! 
      INTEGER                :: ibergs, inbergs
      TYPE(iceberg), POINTER :: this
      !!----------------------------------------------------------------------
      !
      IF (nn_verbose_level == 0) RETURN
      this => first_berg
      IF( ASSOCIATED(this) ) THEN
         WRITE(numicb,'(a," pe=(",i3,")")' ) cd_label, narea
         WRITE(numicb,'(a8,4x,a6,12x,a5,15x,a7,19x,a3,17x,a5,17x,a5,17x,a5)' )   &
            &         'timestep', 'number', 'xi,yj','lon,lat','u,v','uo,vo','ua,va','ui,vi'
      ENDIF
      DO WHILE( ASSOCIATED(this) )
        CALL icb_utl_print_berg(this, kt)
        this => this%next
      END DO
      ibergs = icb_utl_count()
      inbergs = ibergs
      CALL mpp_sum('icbutl', inbergs)
      IF( ibergs > 0 )   WRITE(numicb,'(a," there are",i5," bergs out of",i6," on PE ",i4)')   &
         &                                  cd_label, ibergs, inbergs, narea
      !
   END SUBROUTINE icb_utl_print


   SUBROUTINE icb_utl_incr()
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE icb_utl_incr  ***
      !!
      !! ** Purpose :   
      !!
      ! Small routine for coping with very large integer values labelling icebergs
      ! num_bergs is a array of integers
      ! the first member is incremented in steps of jpnij starting from narea
      ! this means each iceberg is labelled with a unique number
      ! when this gets to the maximum allowed integer the second and subsequent members are 
      ! used to count how many times the member before cycles
      !!----------------------------------------------------------------------
      INTEGER ::   ii, ibig
      !!----------------------------------------------------------------------

      ibig = HUGE(num_bergs(1))
      IF( ibig-jpnij < num_bergs(1) ) THEN
         num_bergs(1) = narea
         DO ii = 2,nkounts
            IF( num_bergs(ii) == ibig ) THEN
               num_bergs(ii) = 0
               IF( ii == nkounts ) CALL ctl_stop('Sorry, run out of iceberg number space')
            ELSE
               num_bergs(ii) = num_bergs(ii) + 1
               EXIT
            ENDIF
         END DO
      ELSE
         num_bergs(1) = num_bergs(1) + jpnij
      ENDIF
      !
   END SUBROUTINE icb_utl_incr


   INTEGER FUNCTION icb_utl_count()
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION icb_utl_count  ***
      !!
      !! ** Purpose :   
      !!----------------------------------------------------------------------
      TYPE(iceberg), POINTER :: this
      !!----------------------------------------------------------------------
      !
      icb_utl_count = 0
      this => first_berg
      DO WHILE( ASSOCIATED(this) )
         icb_utl_count = icb_utl_count+1
         this => this%next
      END DO
      !
   END FUNCTION icb_utl_count


   REAL(wp) FUNCTION icb_utl_mass( first, justbits, justbergs )
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION icb_utl_mass  ***
      !!
      !! ** Purpose :   compute the mass all iceberg, all berg bits or all bergs.
      !!----------------------------------------------------------------------
      TYPE(iceberg)      , POINTER  ::   first
      TYPE(point)        , POINTER  ::   pt
      LOGICAL, INTENT(in), OPTIONAL ::   justbits, justbergs
      !
      TYPE(iceberg), POINTER ::   this
      !!----------------------------------------------------------------------
      icb_utl_mass = 0._wp
      this => first
      !
      IF( PRESENT( justbergs  ) ) THEN
         DO WHILE( ASSOCIATED( this ) )
            pt => this%current_point
            icb_utl_mass = icb_utl_mass + pt%mass         * this%mass_scaling
            this => this%next
         END DO
      ELSEIF( PRESENT(justbits) ) THEN
         DO WHILE( ASSOCIATED( this ) )
            pt => this%current_point
            icb_utl_mass = icb_utl_mass + pt%mass_of_bits * this%mass_scaling
            this => this%next
         END DO
      ELSE
         DO WHILE( ASSOCIATED( this ) )
            pt => this%current_point
            icb_utl_mass = icb_utl_mass + ( pt%mass + pt%mass_of_bits ) * this%mass_scaling
            this => this%next
         END DO
      ENDIF
      !
   END FUNCTION icb_utl_mass


   REAL(wp) FUNCTION icb_utl_heat( first, justbits, justbergs )
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION icb_utl_heat  ***
      !!
      !! ** Purpose :   compute the heat in all iceberg, all bergies or all bergs.
      !!----------------------------------------------------------------------
      TYPE(iceberg)      , POINTER  ::   first
      LOGICAL, INTENT(in), OPTIONAL ::   justbits, justbergs
      !
      TYPE(iceberg)      , POINTER  ::   this
      TYPE(point)        , POINTER  ::   pt
      !!----------------------------------------------------------------------
      icb_utl_heat = 0._wp
      this => first
      !
      IF( PRESENT( justbergs  ) ) THEN
         DO WHILE( ASSOCIATED( this ) )
            pt => this%current_point
            icb_utl_heat = icb_utl_heat + pt%mass         * this%mass_scaling * pt%heat_density
            this => this%next
         END DO
      ELSEIF( PRESENT(justbits) ) THEN
         DO WHILE( ASSOCIATED( this ) )
            pt => this%current_point
            icb_utl_heat = icb_utl_heat + pt%mass_of_bits * this%mass_scaling * pt%heat_density
            this => this%next
         END DO
      ELSE
         DO WHILE( ASSOCIATED( this ) )
            pt => this%current_point
            icb_utl_heat = icb_utl_heat + ( pt%mass + pt%mass_of_bits ) * this%mass_scaling * pt%heat_density
            this => this%next
         END DO
      ENDIF
      !
   END FUNCTION icb_utl_heat

   !!======================================================================
END MODULE icbutl
