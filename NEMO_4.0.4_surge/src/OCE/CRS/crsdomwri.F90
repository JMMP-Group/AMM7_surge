MODULE crsdomwri
   !!======================================================================
   !! Coarse Ocean initialization : write the coarse ocean domain mesh and mask files
   !!======================================================================
   !! History :  3.6   ! 2012-06  (J. Simeon, C. Calone, C Ethe )  from domwri, reduced and modified for coarse grid
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   crs_dom_wri    : create and write mesh and mask file(s)
   !!----------------------------------------------------------------------
   USE timing          ! Timing
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE par_kind, ONLY: wp
   USE lib_mpp         ! MPP library
   USE iom_def
   USE iom
   USE crs         ! coarse grid domain
   USE crsdom         ! coarse grid domain
   USE crslbclnk       ! crs mediator to lbclnk

   IMPLICIT NONE
   PRIVATE

   PUBLIC crs_dom_wri        ! routine called by crsini.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE crs_dom_wri
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE crs_dom_wri  ***
      !!
      !! ** Purpose :   Create the NetCDF file(s) which contain(s) all the
      !!      ocean domain informations (mesh and mask arrays). This (these)
      !!      file(s) is (are) used for visualisation (SAXO software) and
      !!      diagnostic computation.
      !!
      !! ** Method  :   Write in a file all the arrays generated in routines
      !!      crsini for meshes and mask. In three separate files: 
      !!      domain size, horizontal grid-point position,
      !!      masks, depth and vertical scale factors
      !!      
      !! ** Output files :   mesh_hgr_crs.nc, mesh_zgr_crs.nc, mesh_mask.nc
      !!----------------------------------------------------------------------
      INTEGER           ::   ji, jj, jk   ! dummy loop indices
      INTEGER           ::   inum         ! local units for 'mesh_mask.nc' file
      INTEGER           ::   iif, iil, ijf, ijl
      CHARACTER(len=21) ::   clnam        ! filename (mesh and mask informations)
      !                                   !  workspace
      REAL(wp), DIMENSION(jpi_crs,jpj_crs    ) ::   zprt, zprw 
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk) ::   zdepu, zdepv
      !!----------------------------------------------------------------------
      !
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'crs_dom_wri : create NetCDF mesh and mask file'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      
      clnam = 'mesh_mask_crs'  ! filename (mesh and mask informations)
      

      !                            ! ============================
      !                            !  create 'mesh_mask.nc' file
      !                            ! ============================
      !
      CALL iom_open( TRIM(clnam), inum, ldwrt = .TRUE. )
 
      CALL iom_rstput( 0, 0, inum, 'tmask', tmask_crs, ktype = jp_i1 )    ! land-sea mask
      CALL iom_rstput( 0, 0, inum, 'umask', umask_crs, ktype = jp_i1 )
      CALL iom_rstput( 0, 0, inum, 'vmask', vmask_crs, ktype = jp_i1 )
      CALL iom_rstput( 0, 0, inum, 'fmask', fmask_crs, ktype = jp_i1 )
      
      
      tmask_i_crs(:,:) = tmask_crs(:,:,1)
      iif = nn_hls
      iil = nlci_crs - nn_hls + 1
      ijf = nn_hls
      ijl = nlcj_crs - nn_hls + 1
     
      tmask_i_crs( 1:iif ,    :  ) = 0._wp
      tmask_i_crs(iil:jpi_crs,    :  ) = 0._wp
      tmask_i_crs(   :   , 1:ijf ) = 0._wp
      tmask_i_crs(   :   ,ijl:jpj_crs) = 0._wp
      
      
      tpol_crs(1:jpiglo_crs,:) = 1._wp
      fpol_crs(1:jpiglo_crs,:) = 1._wp
      IF( jperio == 3 .OR. jperio == 4 ) THEN
         tpol_crs(jpiglo_crs/2+1:jpiglo_crs,:) = 0._wp
         fpol_crs(       1      :jpiglo_crs,:) = 0._wp
         IF( mjg_crs(nlej_crs) == jpiglo_crs ) THEN
            DO ji = iif+1, iil-1
               tmask_i_crs(ji,nlej_crs-1) = tmask_i_crs(ji,nlej_crs-1) &
               & * tpol_crs(mig_crs(ji),1)
            ENDDO
         ENDIF
      ENDIF
      IF( jperio == 5 .OR. jperio == 6 ) THEN
         tpol_crs(      1       :jpiglo_crs,:)=0._wp
         fpol_crs(jpiglo_crs/2+1:jpiglo_crs,:)=0._wp
      ENDIF
      
      CALL iom_rstput( 0, 0, inum, 'tmaskutil', tmask_i_crs, ktype = jp_i1 )
                                   !    ! unique point mask
      CALL dom_uniq_crs( zprw, 'U' )
      zprt = umask_crs(:,:,1) * zprw
      CALL iom_rstput( 0, 0, inum, 'umaskutil', zprt, ktype = jp_i1 )  
      CALL dom_uniq_crs( zprw, 'V' )
      zprt = vmask_crs(:,:,1) * zprw
      CALL iom_rstput( 0, 0, inum, 'vmaskutil', zprt, ktype = jp_i1 )  
      CALL dom_uniq_crs( zprw, 'F' )
      zprt = fmask_crs(:,:,1) * zprw
      CALL iom_rstput( 0, 0, inum, 'fmaskutil', zprt, ktype = jp_i1 )  
      !========================================================
      !                                                         ! horizontal mesh
      CALL iom_rstput( 0, 0, inum, 'glamt', glamt_crs, ktype = jp_r4 )     !    ! latitude
      CALL iom_rstput( 0, 0, inum, 'glamu', glamu_crs, ktype = jp_r4 )
      CALL iom_rstput( 0, 0, inum, 'glamv', glamv_crs, ktype = jp_r4 )
      CALL iom_rstput( 0, 0, inum, 'glamf', glamf_crs, ktype = jp_r4 )
      
      CALL iom_rstput( 0, 0, inum, 'gphit', gphit_crs, ktype = jp_r4 )     !    ! longitude
      CALL iom_rstput( 0, 0, inum, 'gphiu', gphiu_crs, ktype = jp_r4 )
      CALL iom_rstput( 0, 0, inum, 'gphiv', gphiv_crs, ktype = jp_r4 )
      CALL iom_rstput( 0, 0, inum, 'gphif', gphif_crs, ktype = jp_r4 )
      
      CALL iom_rstput( 0, 0, inum, 'e1t', e1t_crs, ktype = jp_r8 )         !    ! e1 scale factors
      CALL iom_rstput( 0, 0, inum, 'e1u', e1u_crs, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e1v', e1v_crs, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e1f', e1f_crs, ktype = jp_r8 )
      
      CALL iom_rstput( 0, 0, inum, 'e2t', e2t_crs, ktype = jp_r8 )         !    ! e2 scale factors
      CALL iom_rstput( 0, 0, inum, 'e2u', e2u_crs, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e2v', e2v_crs, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e2f', e2f_crs, ktype = jp_r8 )
      
      CALL iom_rstput( 0, 0, inum, 'ff', ff_crs, ktype = jp_r8 )           !    ! coriolis factor

      !========================================================
      !                                                         ! vertical mesh
!     ! note that mbkt is set to 1 over land ==> use surface tmask_crs
      zprt(:,:) = tmask_crs(:,:,1) * REAL( mbkt_crs(:,:) , wp )
      CALL iom_rstput( 0, 0, inum, 'mbathy', zprt, ktype = jp_i2 )     !    ! nb of ocean T-points
      !
      CALL iom_rstput( 0, 0, inum, 'e3t', e3t_crs )      
      CALL iom_rstput( 0, 0, inum, 'e3w', e3w_crs )      
      CALL iom_rstput( 0, 0, inum, 'e3u', e3u_crs )      
      CALL iom_rstput( 0, 0, inum, 'e3v', e3v_crs )      
      !
      CALL iom_rstput( 0, 0, inum, 'gdept', gdept_crs, ktype = jp_r4 ) 
      DO jk = 1,jpk   
         DO jj = 1, jpj_crsm1   
            DO ji = 1, jpi_crsm1  ! jes what to do for fs_jpim1??vector opt.
               zdepu(ji,jj,jk) = MIN( gdept_crs(ji,jj,jk) , gdept_crs(ji+1,jj  ,jk) ) * umask_crs(ji,jj,jk)
               zdepv(ji,jj,jk) = MIN( gdept_crs(ji,jj,jk) , gdept_crs(ji  ,jj+1,jk) ) * vmask_crs(ji,jj,jk)
            END DO   
         END DO   
      END DO
      CALL crs_lbc_lnk( zdepu,'U', 1. )   ;   CALL crs_lbc_lnk( zdepv,'V', 1. ) 
      !
      CALL iom_rstput( 0, 0, inum, 'gdepu', zdepu, ktype = jp_r4 )
      CALL iom_rstput( 0, 0, inum, 'gdepv', zdepv, ktype = jp_r4 )
      CALL iom_rstput( 0, 0, inum, 'gdepw', gdepw_crs, ktype = jp_r4 )
      !
      CALL iom_rstput( 0, 0, inum, 'gdept_1d', gdept_1d )     !    ! reference z-coord.
      CALL iom_rstput( 0, 0, inum, 'gdepw_1d', gdepw_1d )
      CALL iom_rstput( 0, 0, inum, 'e3t_1d'  , e3t_1d   )
      CALL iom_rstput( 0, 0, inum, 'e3w_1d'  , e3w_1d   )
      !
      CALL iom_rstput(  0, 0, inum, 'ocean_volume_t', ocean_volume_crs_t ) 
      CALL iom_rstput(  0, 0, inum, 'facvol_t' , facvol_t  ) 
      CALL iom_rstput(  0, 0, inum, 'facvol_w' , facvol_w  ) 
      CALL iom_rstput(  0, 0, inum, 'facsurfu' , facsurfu  ) 
      CALL iom_rstput(  0, 0, inum, 'facsurfv' , facsurfv  ) 
      CALL iom_rstput(  0, 0, inum, 'e1e2w_msk', e1e2w_msk ) 
      CALL iom_rstput(  0, 0, inum, 'e2e3u_msk', e2e3u_msk ) 
      CALL iom_rstput(  0, 0, inum, 'e1e3v_msk', e1e3v_msk )
      CALL iom_rstput(  0, 0, inum, 'e1e2w'    , e1e2w_crs ) 
      CALL iom_rstput(  0, 0, inum, 'e2e3u'    , e2e3u_crs ) 
      CALL iom_rstput(  0, 0, inum, 'e1e3v'    , e1e3v_crs )
      CALL iom_rstput(  0, 0, inum, 'bt'       , bt_crs    )
      CALL iom_rstput(  0, 0, inum, 'r1_bt'    , r1_bt_crs )
      !
      CALL iom_rstput(  0, 0, inum, 'crs_surfu_wgt', crs_surfu_wgt ) 
      CALL iom_rstput(  0, 0, inum, 'crs_surfv_wgt', crs_surfv_wgt ) 
      CALL iom_rstput(  0, 0, inum, 'crs_volt_wgt' , crs_volt_wgt  ) 
      !                                     ! ============================
      !                                     !        close the files 
      !                                     ! ============================
      CALL iom_close( inum )
      !
   END SUBROUTINE crs_dom_wri


   SUBROUTINE dom_uniq_crs( puniq, cdgrd )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE crs_dom_uniq_crs  ***
      !!                   
      !! ** Purpose :   identify unique point of a grid (TUVF)
      !!
      !! ** Method  :   1) apply crs_lbc_lnk on an array with different values for each element
      !!                2) check which elements have been changed
      !!----------------------------------------------------------------------
      CHARACTER(len=1)        , INTENT(in   ) ::   cdgrd   ! 
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   puniq   ! 
      !
      REAL(wp) ::  zshift   ! shift value link to the process number
      INTEGER  ::  ji       ! dummy loop indices
      LOGICAL, DIMENSION(SIZE(puniq,1),SIZE(puniq,2),1) ::  lldbl  ! store whether each point is unique or not
      REAL(wp), DIMENSION(jpi_crs,jpj_crs) :: ztstref
      !!----------------------------------------------------------------------
      !
      ! build an array with different values for each element 
      ! in mpp: make sure that these values are different even between process
      ! -> apply a shift value according to the process number
      zshift = jpi_crs * jpj_crs * ( narea - 1 )
      ztstref(:,:) = RESHAPE( (/ (zshift + REAL(ji,wp), ji = 1, jpi_crs*jpj_crs) /), (/ jpi_crs, jpj_crs /) )
      !
      puniq(:,:) = ztstref(:,:)                   ! default definition
      CALL crs_lbc_lnk( puniq,cdgrd, 1. )            ! apply boundary conditions
      lldbl(:,:,1) = puniq(:,:) == ztstref(:,:)   ! check which values have been changed 
      !
      puniq(:,:) = 1.                             ! default definition
      ! fill only the inner part of the cpu with llbl converted into real 
      puniq(nldi_crs:nlei_crs,nldj_crs:nlej_crs) = REAL( COUNT( lldbl(nldi_crs:nlei_crs,nldj_crs:nlej_crs,:), dim = 3 ) , wp )
      !
   END SUBROUTINE dom_uniq_crs

   !!======================================================================

END MODULE crsdomwri


