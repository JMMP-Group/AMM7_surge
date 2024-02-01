MODULE crsini   
   !!======================================================================
   !!                         ***  MODULE crsini  ***
   !!            Manage the grid coarsening module initialization
   !!======================================================================
   !!  History     2012-05   (J. Simeon, G. Madec, C. Ethe, C. Calone) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  crs_init    : 
   !!----------------------------------------------------------------------
   USE par_kind, ONLY: wp
   USE par_oce                  ! For parameter jpi,jpj
   USE dom_oce                  ! For parameters in par_oce
   USE crs                      ! Coarse grid domain
   USE phycst, ONLY: omega, rad ! physical constants
   USE crsdom
   USE crsdomwri
   USE crslbclnk
   !
   USE iom
   USE in_out_manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC   crs_init   ! called by nemogcm.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE crs_init 
      !!-------------------------------------------------------------------
      !!                     *** SUBROUTINE crs_init
      !!  ** Purpose : Initialization of the grid coarsening module  
      !!               1. Read namelist
      !!               X2. MOVE TO crs_dom.F90 Set the domain definitions for coarse grid
      !!                 a. Define the coarse grid starting/ending indices on parent grid
      !!                    Here is where the T-pivot or F-pivot grids are discerned
      !!                 b. TODO.  Include option for north-centric or equator-centric binning.
      !!                 (centered only for odd reduction factors; even reduction bins bias equator to the south)
      !!               3. Mask and mesh creation. => calls to crsfun
      !!                  a. Use crsfun_mask to generate tmask,umask, vmask, fmask.
      !!                  b. Use crsfun_coordinates to get coordinates
      !!                  c. Use crsfun_UV to get horizontal scale factors
      !!                  d. Use crsfun_TW to get initial vertical scale factors   
      !!               4. Volumes and weights jes.... TODO. Updates for vvl? Where to do this? crsstp.F90?
      !!                  a. Calculate initial coarse grid box volumes: pvol_T, pvol_W
      !!                  b. Calculate initial coarse grid surface-averaging weights
      !!                  c. Calculate initial coarse grid volume-averaging weights
      !!                  
      !!               X5. MOVE TO crs_dom_wri.F90 Using iom_rstput output the initial meshmask.
      !!               ?. Another set of "masks" to generate
      !!                  are the u- and v- surface areas for U- and V- area-weighted means.
      !!                  Need to put this somewhere in section 3?
      !! jes. What do to about the vvl?  GM.  could separate the weighting (denominator), so
      !! output C*dA or C*dV as summation not mran, then do mean (division) at moment of output.
      !! As is, crsfun takes into account vvl.   
      !!      Talked about pre-setting the surface array to avoid IF/ENDIF and division.
      !!      But have then to make that preset array here and elsewhere.
      !!      that is called every timestep...
      !! 
      !!               - Read in pertinent data ?
      !!-------------------------------------------------------------------
      INTEGER  :: ji,jj,jk      ! dummy indices
      INTEGER  :: ierr                                ! allocation error status
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ze3t, ze3u, ze3v, ze3w

      NAMELIST/namcrs/ nn_factx, nn_facty, nn_binref, ln_msh_crs, nn_crs_kz, ln_crs_wn
      !!----------------------------------------------------------------------
      !
     !---------------------------------------------------------
     ! 1. Read Namelist file
     !---------------------------------------------------------
     !
      REWIND( numnam_ref )              ! Namelist namrun in reference namelist : Parameters of the run
      READ  ( numnam_ref, namcrs, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcrs in reference namelist' )
      REWIND( numnam_cfg )              ! Namelist namrun in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, namcrs, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcrs in configuration namelist' )
      IF(lwm) WRITE ( numond, namcrs )

     IF(lwp) THEN
        WRITE(numout,*)
        WRITE(numout,*) 'crs_init : Initializing the grid coarsening module'
        WRITE(numout,*) '~~~~~~~~'
        WRITE(numout,*) '   Namelist namcrs '
        WRITE(numout,*) '      coarsening factor in i-direction      nn_factx   = ', nn_factx
        WRITE(numout,*) '      coarsening factor in j-direction      nn_facty   = ', nn_facty
        WRITE(numout,*) '      bin centering preference              nn_binref  = ', nn_binref
        WRITE(numout,*) '      create a mesh file (=T)               ln_msh_crs = ', ln_msh_crs
        WRITE(numout,*) '      type of Kz coarsening (0,1,2)         nn_crs_kz  = ', nn_crs_kz
        WRITE(numout,*) '      wn coarsened or computed using hdivn  ln_crs_wn  = ', ln_crs_wn
     ENDIF
              
     rfactx_r = 1. / nn_factx
     rfacty_r = 1. / nn_facty

     !---------------------------------------------------------
     ! 2. Define Global Dimensions of the coarsened grid
     !---------------------------------------------------------
     CALL crs_dom_def      

     !---------------------------------------------------------
     ! 3. Mask and Mesh
     !---------------------------------------------------------

     !     Set up the masks and meshes     

     !  3.a. Get the masks   
 
     CALL crs_dom_msk


     !  3.b. Get the coordinates
     !      Odd-numbered reduction factor, center coordinate on T-cell
     !      Even-numbered reduction factor, center coordinate on U-,V- faces or f-corner.
     !      
     IF ( nresty /= 0 .AND. nrestx /= 0 ) THEN
        CALL crs_dom_coordinates( gphit, glamt, 'T', gphit_crs, glamt_crs ) 
        CALL crs_dom_coordinates( gphiu, glamu, 'U', gphiu_crs, glamu_crs )       
        CALL crs_dom_coordinates( gphiv, glamv, 'V', gphiv_crs, glamv_crs ) 
        CALL crs_dom_coordinates( gphif, glamf, 'F', gphif_crs, glamf_crs ) 
     ELSEIF ( nresty /= 0 .AND. nrestx == 0 ) THEN
        CALL crs_dom_coordinates( gphiu, glamu, 'T', gphit_crs, glamt_crs )
        CALL crs_dom_coordinates( gphiu, glamu, 'U', gphiu_crs, glamu_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'V', gphiv_crs, glamv_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'F', gphif_crs, glamf_crs )
     ELSEIF ( nresty == 0 .AND. nrestx /= 0 ) THEN
        CALL crs_dom_coordinates( gphiv, glamv, 'T', gphit_crs, glamt_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'U', gphiu_crs, glamu_crs )
        CALL crs_dom_coordinates( gphiv, glamv, 'V', gphiv_crs, glamv_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'F', gphif_crs, glamf_crs )
     ELSE 
        CALL crs_dom_coordinates( gphif, glamf, 'T', gphit_crs, glamt_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'U', gphiu_crs, glamu_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'V', gphiv_crs, glamv_crs )
        CALL crs_dom_coordinates( gphif, glamf, 'F', gphif_crs, glamf_crs )
     ENDIF


     !  3.c. Get the horizontal mesh

     !      3.c.1 Horizontal scale factors

     CALL crs_dom_hgr( e1t, e2t, 'T', e1t_crs, e2t_crs )
     CALL crs_dom_hgr( e1u, e2u, 'U', e1u_crs, e2u_crs )
     CALL crs_dom_hgr( e1v, e2v, 'V', e1v_crs, e2v_crs )
     CALL crs_dom_hgr( e1f, e2f, 'F', e1f_crs, e2f_crs )

     e1e2t_crs(:,:) = e1t_crs(:,:) * e2t_crs(:,:)
     
     
     !      3.c.2 Coriolis factor  

!!gm  Not sure CRS needs Coriolis parameter....
!!gm  If needed, then update this to have Coriolis at both f- and t-points

      ff_crs(:,:) = 2. * omega * SIN( rad * gphif_crs(:,:) )

      CALL ctl_warn( 'crsini: CAUTION, CRS only designed for Coriolis defined on the sphere' ) 
 

     !    3.d.1 mbathy ( vertical k-levels of bathymetry )     

     CALL crs_dom_bat
     
     !
     ze3t(:,:,:) = e3t_n(:,:,:)
     ze3u(:,:,:) = e3u_n(:,:,:)
     ze3v(:,:,:) = e3v_n(:,:,:)
     ze3w(:,:,:) = e3w_n(:,:,:)

     !    3.d.2   Surfaces 
     CALL crs_dom_sfc( tmask, 'W', e1e2w_crs, e1e2w_msk, p_e1=e1t, p_e2=e2t  )
     CALL crs_dom_sfc( umask, 'U', e2e3u_crs, e2e3u_msk, p_e2=e2u, p_e3=ze3u )
     CALL crs_dom_sfc( vmask, 'V', e1e3v_crs, e1e3v_msk, p_e1=e1v, p_e3=ze3v )
   
     facsurfu(:,:,:) = umask_crs(:,:,:) * e2e3u_msk(:,:,:) / e2e3u_crs(:,:,:)
     facsurfv(:,:,:) = vmask_crs(:,:,:) * e1e3v_msk(:,:,:) / e1e3v_crs(:,:,:)

     !    3.d.3   Vertical scale factors
     !
     CALL crs_dom_e3( e1t, e2t, ze3t, e1e2w_crs, 'T', tmask, e3t_crs, e3t_max_crs)
     CALL crs_dom_e3( e1u, e2u, ze3u, e2e3u_crs, 'U', umask, e3u_crs, e3u_max_crs)
     CALL crs_dom_e3( e1v, e2v, ze3v, e1e3v_crs, 'V', vmask, e3v_crs, e3v_max_crs)
     CALL crs_dom_e3( e1t, e2t, ze3w, e1e2w_crs, 'W', tmask, e3w_crs, e3w_max_crs)

     ! Replace 0 by e3t_0 or e3w_0
     DO jk = 1, jpk
        DO ji = 1, jpi_crs
           DO jj = 1, jpj_crs
              IF( e3t_crs(ji,jj,jk) == 0._wp ) e3t_crs(ji,jj,jk) = e3t_1d(jk)
              IF( e3w_crs(ji,jj,jk) == 0._wp ) e3w_crs(ji,jj,jk) = e3w_1d(jk)
              IF( e3u_crs(ji,jj,jk) == 0._wp ) e3u_crs(ji,jj,jk) = e3t_1d(jk)
              IF( e3v_crs(ji,jj,jk) == 0._wp ) e3v_crs(ji,jj,jk) = e3t_1d(jk)
           ENDDO
        ENDDO
     ENDDO

     !    3.d.3   Vertical depth (meters)
     CALL crs_dom_ope( gdept_0, 'MAX', 'T', tmask, gdept_crs, p_e3=ze3t, psgn=1.0 ) 
     CALL crs_dom_ope( gdepw_0, 'MAX', 'W', tmask, gdepw_crs, p_e3=ze3w, psgn=1.0 )


     !---------------------------------------------------------
     ! 4.  Coarse grid ocean volume and averaging weights
     !---------------------------------------------------------
     ! 4.a. Ocean volume or area unmasked and masked
     CALL crs_dom_facvol( tmask, 'T', e1t, e2t, ze3t, ocean_volume_crs_t, facvol_t )
     !
     bt_crs(:,:,:) = ocean_volume_crs_t(:,:,:) * facvol_t(:,:,:)
     !
     r1_bt_crs(:,:,:) = 0._wp 
     WHERE( bt_crs /= 0._wp ) r1_bt_crs(:,:,:) = 1._wp / bt_crs(:,:,:)

     CALL crs_dom_facvol( tmask, 'W', e1t, e2t, ze3w, ocean_volume_crs_w, facvol_w )
     !
     !---------------------------------------------------------
     ! 5.  Write out coarse meshmask  (see OCE/DOM/domwri.F90 for ideas later)
     !---------------------------------------------------------

     IF( ln_msh_crs ) THEN 
        CALL dom_grid_crs   ! Save the parent grid information  & Switch to coarse grid domain
        CALL crs_dom_wri     
        CALL dom_grid_glo   ! Return to parent grid domain
     ENDIF
     
      !---------------------------------------------------------
      ! 7. Finish and clean-up
      !---------------------------------------------------------
      !
   END SUBROUTINE crs_init
    
   !!======================================================================
END MODULE crsini
