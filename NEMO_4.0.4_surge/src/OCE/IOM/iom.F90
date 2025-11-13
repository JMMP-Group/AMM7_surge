MODULE iom
   !!======================================================================
   !!                    ***  MODULE  iom ***
   !! Input/Output manager :  Library to read input files
   !!======================================================================
   !! History :  2.0  ! 2005-12  (J. Belier) Original code
   !!            2.0  ! 2006-02  (S. Masson) Adaptation to NEMO
   !!            3.0  ! 2007-07  (D. Storkey) Changes to iom_gettime
   !!            3.4  ! 2012-12  (R. Bourdalle-Badie and G. Reffray)  add C1D case  
   !!            3.6  ! 2014-15  DIMG format removed
   !!            3.6  ! 2015-15  (J. Harle) Added procedure to read REAL attributes
   !!            4.0  ! 2017-11  (M. Andrejczuk) Extend IOM interface to write any 3D fields
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   iom_open       : open a file read only
   !!   iom_close      : close a file or all files opened by iom
   !!   iom_get        : read a field (interfaced to several routines)
   !!   iom_varid      : get the id of a variable in a file
   !!   iom_rstput     : write a field in a restart file (interfaced to several routines)
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE c1d             ! 1D vertical configuration
   USE flo_oce         ! floats module declarations
   USE lbclnk          ! lateal boundary condition / mpp exchanges
   USE iom_def         ! iom variables definitions
   USE iom_nf90        ! NetCDF format with native NetCDF library
   USE in_out_manager  ! I/O manager
   USE lib_mpp           ! MPP library
#if defined key_iomput
   USE sbc_oce  , ONLY :   nn_fsbc         ! ocean space and time domain
   USE trc_oce  , ONLY :   nn_dttrc        !  !: frequency of step on passive tracers
   USE icb_oce  , ONLY :   nclasses, class_num       !  !: iceberg classes
#if defined key_si3
   USE ice      , ONLY :   jpl
#endif
   USE domngb          ! ocean space and time domain
   USE phycst          ! physical constants
   USE dianam          ! build name of file
   USE xios
# endif
   USE ioipsl, ONLY :  ju2ymds    ! for calendar
   USE crs             ! Grid coarsening
#if defined key_top
   USE trc, ONLY    :  profsed
#endif
   USE lib_fortran 
   USE diurnal_bulk, ONLY : ln_diurnal_only, ln_diurnal

   IMPLICIT NONE
   PUBLIC   !   must be public to be able to access iom_def through iom
   
#if defined key_iomput
   LOGICAL, PUBLIC, PARAMETER ::   lk_iomput = .TRUE.        !: iom_put flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_iomput = .FALSE.       !: iom_put flag
#endif
   PUBLIC iom_init, iom_swap, iom_open, iom_close, iom_setkt, iom_varid, iom_get, iom_get_var
   PUBLIC iom_chkatt, iom_getatt, iom_putatt, iom_getszuld, iom_rstput, iom_delay_rst, iom_put
   PUBLIC iom_use, iom_context_finalize, iom_miss_val

   PRIVATE iom_rp0d, iom_rp1d, iom_rp2d, iom_rp3d
   PRIVATE iom_g0d, iom_g1d, iom_g2d, iom_g3d, iom_get_123d
   PRIVATE iom_p1d, iom_p2d, iom_p3d, iom_p4d
#if defined key_iomput
   PRIVATE iom_set_domain_attr, iom_set_axis_attr, iom_set_field_attr, iom_set_file_attr, iom_get_file_attr, iom_set_grid_attr
   PRIVATE set_grid, set_grid_bounds, set_scalar, set_xmlatt, set_mooring, iom_update_file_name, iom_sdate
   PRIVATE iom_set_rst_context, iom_set_rstw_active, iom_set_rstr_active
# endif
   PUBLIC iom_set_rstw_var_active, iom_set_rstw_core, iom_set_rst_vars

   INTERFACE iom_get
      MODULE PROCEDURE iom_g0d, iom_g1d, iom_g2d, iom_g3d
   END INTERFACE
   INTERFACE iom_getatt
      MODULE PROCEDURE iom_g0d_iatt, iom_g1d_iatt, iom_g0d_ratt, iom_g1d_ratt, iom_g0d_catt
   END INTERFACE
   INTERFACE iom_putatt
      MODULE PROCEDURE iom_p0d_iatt, iom_p1d_iatt, iom_p0d_ratt, iom_p1d_ratt, iom_p0d_catt
   END INTERFACE
   INTERFACE iom_rstput
      MODULE PROCEDURE iom_rp0d, iom_rp1d, iom_rp2d, iom_rp3d
   END INTERFACE
   INTERFACE iom_put
      MODULE PROCEDURE iom_p0d, iom_p1d, iom_p2d, iom_p3d, iom_p4d
   END INTERFACE iom_put
  
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE iom_init( cdname, fname, ld_tmppatch ) 
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE   ***
      !!
      !! ** Purpose :   
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=*),           INTENT(in)  :: cdname
      CHARACTER(len=*), OPTIONAL, INTENT(in)  :: fname
      LOGICAL         , OPTIONAL, INTENT(in)  :: ld_tmppatch
#if defined key_iomput
      !
      TYPE(xios_duration) :: dtime    = xios_duration(0, 0, 0, 0, 0, 0)
      TYPE(xios_date)     :: start_date
      CHARACTER(len=lc) :: clname
      INTEGER             :: irefyear, irefmonth, irefday
      INTEGER           :: ji
      LOGICAL :: llrst_context              ! is context related to restart
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zt_bnds, zw_bnds
      LOGICAL ::   ll_tmppatch = .TRUE.    !: seb: patch before we remove periodicity
      INTEGER ::   nldi_save, nlei_save    !:      and close boundaries in output files
      INTEGER ::   nldj_save, nlej_save    !:
      !!----------------------------------------------------------------------
      !
      ! seb: patch before we remove periodicity and close boundaries in output files
      IF( PRESENT(ld_tmppatch) ) THEN   ;   ll_tmppatch = ld_tmppatch
      ELSE                              ;   ll_tmppatch = .TRUE.
      ENDIF
      IF ( ll_tmppatch ) THEN
         nldi_save = nldi   ;   nlei_save = nlei
         nldj_save = nldj   ;   nlej_save = nlej
         IF( nimpp           ==      1 ) nldi = 1
         IF( nimpp + jpi - 1 == jpiglo ) nlei = jpi
         IF( njmpp           ==      1 ) nldj = 1
         IF( njmpp + jpj - 1 == jpjglo ) nlej = jpj
      ENDIF
      !
      ALLOCATE( zt_bnds(2,jpk), zw_bnds(2,jpk) )
      !
      clname = cdname
      IF( TRIM(Agrif_CFixed()) /= '0' )   clname = TRIM(Agrif_CFixed())//"_"//TRIM(cdname)
      CALL xios_context_initialize(TRIM(clname), mpi_comm_oce)
      CALL iom_swap( cdname )
      llrst_context =  (TRIM(cdname) == TRIM(crxios_context) .OR. TRIM(cdname) == TRIM(cwxios_context))

      ! Calendar type is now defined in xml file 
      IF (.NOT.(xios_getvar('ref_year' ,irefyear ))) irefyear  = 1900
      IF (.NOT.(xios_getvar('ref_month',irefmonth))) irefmonth = 01
      IF (.NOT.(xios_getvar('ref_day'  ,irefday  ))) irefday   = 01

      SELECT CASE ( nleapy )        ! Choose calendar for IOIPSL
      CASE ( 1)   ; CALL xios_define_calendar( TYPE = "Gregorian", time_origin = xios_date(irefyear,irefmonth,irefday,00,00,00), &
          &                                    start_date = xios_date(nyear,nmonth,nday,nhour,nminute,0) )
      CASE ( 0)   ; CALL xios_define_calendar( TYPE = "NoLeap"   , time_origin = xios_date(irefyear,irefmonth,irefday,00,00,00), &
          &                                    start_date = xios_date(nyear,nmonth,nday,nhour,nminute,0) )
      CASE (30)   ; CALL xios_define_calendar( TYPE = "D360"     , time_origin = xios_date(irefyear,irefmonth,irefday,00,00,00), &
          &                                    start_date = xios_date(nyear,nmonth,nday,nhour,nminute,0) )
      END SELECT

      ! horizontal grid definition
      IF(.NOT.llrst_context) CALL set_scalar
      !
      IF( TRIM(cdname) == TRIM(cxios_context) ) THEN  
         CALL set_grid( "T", glamt, gphit, .FALSE., .FALSE. ) 
         CALL set_grid( "U", glamu, gphiu, .FALSE., .FALSE. )
         CALL set_grid( "V", glamv, gphiv, .FALSE., .FALSE. )
         CALL set_grid( "W", glamt, gphit, .FALSE., .FALSE. )
         CALL set_grid_znl( gphit )
         !
         IF( ln_cfmeta ) THEN   ! Add additional grid metadata
            CALL iom_set_domain_attr("grid_T", area = e1e2t(nldi:nlei, nldj:nlej))
            CALL iom_set_domain_attr("grid_U", area = e1e2u(nldi:nlei, nldj:nlej))
            CALL iom_set_domain_attr("grid_V", area = e1e2v(nldi:nlei, nldj:nlej))
            CALL iom_set_domain_attr("grid_W", area = e1e2t(nldi:nlei, nldj:nlej))
            CALL set_grid_bounds( "T", glamf, gphif, glamt, gphit )
            CALL set_grid_bounds( "U", glamv, gphiv, glamu, gphiu )
            CALL set_grid_bounds( "V", glamu, gphiu, glamv, gphiv )
            CALL set_grid_bounds( "W", glamf, gphif, glamt, gphit )
         ENDIF
      ENDIF
      !
      IF( TRIM(cdname) == TRIM(cxios_context)//"_crs" ) THEN  
         CALL dom_grid_crs   ! Save the parent grid information  & Switch to coarse grid domain
         !
         CALL set_grid( "T", glamt_crs, gphit_crs, .FALSE., .FALSE. ) 
         CALL set_grid( "U", glamu_crs, gphiu_crs, .FALSE., .FALSE. ) 
         CALL set_grid( "V", glamv_crs, gphiv_crs, .FALSE., .FALSE. ) 
         CALL set_grid( "W", glamt_crs, gphit_crs, .FALSE., .FALSE. ) 
         CALL set_grid_znl( gphit_crs )
          !
         CALL dom_grid_glo   ! Return to parent grid domain
         !
         IF( ln_cfmeta .AND. .NOT. llrst_context) THEN   ! Add additional grid metadata
            CALL iom_set_domain_attr("grid_T", area = e1e2t_crs(nldi:nlei, nldj:nlej))
            CALL iom_set_domain_attr("grid_U", area = e1u_crs(nldi:nlei, nldj:nlej) * e2u_crs(nldi:nlei, nldj:nlej))
            CALL iom_set_domain_attr("grid_V", area = e1v_crs(nldi:nlei, nldj:nlej) * e2v_crs(nldi:nlei, nldj:nlej))
            CALL iom_set_domain_attr("grid_W", area = e1e2t_crs(nldi:nlei, nldj:nlej))
            CALL set_grid_bounds( "T", glamf_crs, gphif_crs, glamt_crs, gphit_crs )
            CALL set_grid_bounds( "U", glamv_crs, gphiv_crs, glamu_crs, gphiu_crs )
            CALL set_grid_bounds( "V", glamu_crs, gphiu_crs, glamv_crs, gphiv_crs )
            CALL set_grid_bounds( "W", glamf_crs, gphif_crs, glamt_crs, gphit_crs )
         ENDIF
      ENDIF
      !
      ! vertical grid definition
      IF(.NOT.llrst_context) THEN
          CALL iom_set_axis_attr( "deptht",  paxis = gdept_1d )
          CALL iom_set_axis_attr( "depthu",  paxis = gdept_1d )
          CALL iom_set_axis_attr( "depthv",  paxis = gdept_1d )
          CALL iom_set_axis_attr( "depthw",  paxis = gdepw_1d )

          ! Add vertical grid bounds
          zt_bnds(2,:        ) = gdept_1d(:)
          zt_bnds(1,2:jpk    ) = gdept_1d(1:jpkm1)
          zt_bnds(1,1        ) = gdept_1d(1) - e3w_1d(1)
          zw_bnds(1,:        ) = gdepw_1d(:)
          zw_bnds(2,1:jpkm1  ) = gdepw_1d(2:jpk)
          zw_bnds(2,jpk:     ) = gdepw_1d(jpk) + e3t_1d(jpk)
          CALL iom_set_axis_attr( "deptht", bounds=zw_bnds )
          CALL iom_set_axis_attr( "depthu", bounds=zw_bnds )
          CALL iom_set_axis_attr( "depthv", bounds=zw_bnds )
          CALL iom_set_axis_attr( "depthw", bounds=zt_bnds )
          CALL iom_set_axis_attr( "nfloat", (/ (REAL(ji,wp), ji=1,jpnfl) /) )
# if defined key_si3
          CALL iom_set_axis_attr( "ncatice", (/ (REAL(ji,wp), ji=1,jpl) /) )
          ! SIMIP diagnostics (4 main arctic straits)
          CALL iom_set_axis_attr( "nstrait", (/ (REAL(ji,wp), ji=1,4) /) )
# endif
#if defined key_top
          IF( ALLOCATED(profsed) ) CALL iom_set_axis_attr( "profsed", paxis = profsed )
#endif
          CALL iom_set_axis_attr( "icbcla", class_num )
          CALL iom_set_axis_attr( "iax_20C", (/ REAL(20,wp) /) )   ! strange syntaxe and idea...
          CALL iom_set_axis_attr( "iax_26C", (/ REAL(26,wp) /) )   ! strange syntaxe and idea...
          CALL iom_set_axis_attr( "iax_28C", (/ REAL(28,wp) /) )   ! strange syntaxe and idea...
          CALL iom_set_axis_attr( "basin"  , (/ (REAL(ji,wp), ji=1,5) /) )
      ENDIF
      !
      ! automatic definitions of some of the xml attributs
      IF( TRIM(cdname) == TRIM(crxios_context) ) THEN
!set names of the fields in restart file IF using XIOS to read data
          CALL iom_set_rst_context(.TRUE.)
          CALL iom_set_rst_vars(rst_rfields)
!set which fields are to be read from restart file
          CALL iom_set_rstr_active()
      ELSE IF( TRIM(cdname) == TRIM(cwxios_context) ) THEN
!set names of the fields in restart file IF using XIOS to write data
          CALL iom_set_rst_context(.FALSE.)
          CALL iom_set_rst_vars(rst_wfields)
!set which fields are to be written to a restart file
          CALL iom_set_rstw_active(fname)
      ELSE
          CALL set_xmlatt
      ENDIF
      !
      ! end file definition
      dtime%second = rdt
      CALL xios_set_timestep( dtime )
      CALL xios_close_context_definition()
      CALL xios_update_calendar( 0 )
      !
      DEALLOCATE( zt_bnds, zw_bnds )
      !
      IF ( ll_tmppatch ) THEN
         nldi = nldi_save   ;   nlei = nlei_save
         nldj = nldj_save   ;   nlej = nlej_save
      ENDIF
#endif
      !
   END SUBROUTINE iom_init

   SUBROUTINE iom_set_rstw_var_active(field)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_set_rstw_var_active  ***
      !!
      !! ** Purpose :  enable variable in restart file when writing with XIOS 
      !!---------------------------------------------------------------------
   CHARACTER(len = *), INTENT(IN) :: field
   INTEGER :: i
   LOGICAL :: llis_set
   CHARACTER(LEN=256) :: clinfo    ! info character

#if defined key_iomput
   llis_set = .FALSE.

   DO i = 1, max_rst_fields
       IF(TRIM(rst_wfields(i)%vname) == field) THEN 
          rst_wfields(i)%active = .TRUE.
          llis_set = .TRUE.
          EXIT
       ENDIF
   ENDDO
!Warn if variable is not in defined in rst_wfields
   IF(.NOT.llis_set) THEN
      WRITE(ctmp1,*) 'iom_set_rstw_var_active: variable ', field ,' is available for writing but not defined' 
      CALL ctl_stop( 'iom_set_rstw_var_active:', ctmp1 )
   ENDIF
#else
        clinfo = 'iom_set_rstw_var_active: key_iomput is needed to use XIOS restart read/write functionality'
        CALL ctl_stop('STOP', TRIM(clinfo))
#endif

   END SUBROUTINE iom_set_rstw_var_active

   SUBROUTINE iom_set_rstr_active()
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_set_rstr_active  ***
      !!
      !! ** Purpose :  define file name in XIOS context for reading restart file,
      !!               enable variables present in restart file for reading with XIOS 
      !!---------------------------------------------------------------------

!sets enabled = .TRUE. for each field in restart file
   CHARACTER(len=256) :: rst_file

#if defined key_iomput
   TYPE(xios_field) :: field_hdl
   TYPE(xios_file) :: file_hdl
   TYPE(xios_filegroup) :: filegroup_hdl
   INTEGER :: i
   CHARACTER(lc)  ::   clpath

        clpath = TRIM(cn_ocerst_indir)
        IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'
        IF( TRIM(Agrif_CFixed()) == '0' ) THEN
           rst_file = TRIM(clpath)//TRIM(cn_ocerst_in)
        ELSE
           rst_file = TRIM(clpath)//TRIM(Agrif_CFixed())//'_'//TRIM(cn_ocerst_in)
        ENDIF
!set name of the restart file and enable available fields
        if(lwp) WRITE(numout,*) 'Setting restart filename (for XIOS) to: ',rst_file
        CALL xios_get_handle("file_definition", filegroup_hdl )
        CALL xios_add_child(filegroup_hdl, file_hdl, 'rrestart')
        CALL xios_set_file_attr( "rrestart", name=trim(rst_file), type="one_file", &
             par_access="collective", enabled=.TRUE., mode="read",                 &
             output_freq=xios_timestep)
!define variables for restart context
        DO i = 1, max_rst_fields
         IF( TRIM(rst_rfields(i)%vname) /= "NO_NAME") THEN
           IF( iom_varid( numror, TRIM(rst_rfields(i)%vname), ldstop = .FALSE. ) > 0 ) THEN
                CALL xios_add_child(file_hdl, field_hdl, TRIM(rst_rfields(i)%vname))
                SELECT CASE (TRIM(rst_rfields(i)%grid))
                 CASE ("grid_N_3D")
                    CALL xios_set_attr (field_hdl, enabled = .TRUE., name = TRIM(rst_rfields(i)%vname), &
                        domain_ref="grid_N", axis_ref="nav_lev", operation = "instant")
                 CASE ("grid_N")
                    CALL xios_set_attr (field_hdl, enabled = .TRUE., name = TRIM(rst_rfields(i)%vname), &
                        domain_ref="grid_N", operation = "instant") 
                CASE ("grid_vector")
                    CALL xios_set_attr (field_hdl, enabled = .TRUE., name = TRIM(rst_rfields(i)%vname), &
                         axis_ref="nav_lev", operation = "instant")
                 CASE ("grid_scalar")
                    CALL xios_set_attr (field_hdl, enabled = .TRUE., name = TRIM(rst_rfields(i)%vname), &
                        scalar_ref = "grid_scalar", operation = "instant")
                END SELECT
                IF(lwp) WRITE(numout,*) 'XIOS read: ', TRIM(rst_rfields(i)%vname), ' enabled in ', TRIM(rst_file)
           ENDIF
         ENDIF
        END DO
#endif
   END SUBROUTINE iom_set_rstr_active

   SUBROUTINE iom_set_rstw_core(cdmdl)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_set_rstw_core  ***
      !!
      !! ** Purpose :  set variables which are always in restart file 
      !!---------------------------------------------------------------------
   CHARACTER (len=*), INTENT (IN) :: cdmdl ! model OPA or SAS
   CHARACTER(LEN=256)             :: clinfo    ! info character
#if defined key_iomput
   IF(cdmdl == "OPA") THEN
!from restart.F90
   CALL iom_set_rstw_var_active("rdt")
   CALL iom_set_rstw_var_active("neos")

   IF ( .NOT. ln_diurnal_only ) THEN
        CALL iom_set_rstw_var_active('ub'  )
        CALL iom_set_rstw_var_active('vb'  )
        CALL iom_set_rstw_var_active('tb'  )
        CALL iom_set_rstw_var_active('sb'  )
        CALL iom_set_rstw_var_active('sshb')
        !
        CALL iom_set_rstw_var_active('un'  )
        CALL iom_set_rstw_var_active('vn'  )
        CALL iom_set_rstw_var_active('tn'  )
        CALL iom_set_rstw_var_active('sn'  )
        CALL iom_set_rstw_var_active('sshn')
        CALL iom_set_rstw_var_active('rhop')
     ! extra variable needed for the ice sheet coupling
        IF ( ln_iscpl ) THEN
             CALL iom_set_rstw_var_active('tmask')
             CALL iom_set_rstw_var_active('umask')
             CALL iom_set_rstw_var_active('vmask')
             CALL iom_set_rstw_var_active('smask')
             CALL iom_set_rstw_var_active('e3t_n')
             CALL iom_set_rstw_var_active('e3u_n')
             CALL iom_set_rstw_var_active('e3v_n')
             CALL iom_set_rstw_var_active('gdepw_n')
        END IF
      ENDIF
      IF(ln_diurnal) CALL iom_set_rstw_var_active('Dsst')
!from trasbc.F90
         CALL iom_set_rstw_var_active('sbc_hc_b')
         CALL iom_set_rstw_var_active('sbc_sc_b')
   ENDIF
#else
        clinfo = 'iom_set_rstw_core: key_iomput is needed to use XIOS restart read/write functionality'
        CALL ctl_stop('STOP', TRIM(clinfo))
#endif
   END SUBROUTINE iom_set_rstw_core

   SUBROUTINE iom_set_rst_vars(fields)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE iom_set_rst_vars   ***
      !!
      !! ** Purpose :  Fill array fields with the information about all 
      !!               possible variables and corresponding grids definition 
      !!               for reading/writing restart with XIOS
      !!---------------------------------------------------------------------
   TYPE(RST_FIELD), INTENT(INOUT) :: fields(max_rst_fields)
   INTEGER :: i

        i = 0
        i = i + 1; fields(i)%vname="rdt";            fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="neos";           fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="un";             fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="ub";             fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="vn";             fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="vb";             fields(i)%grid="grid_N_3D"  
        i = i + 1; fields(i)%vname="tn";             fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="tb";             fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="sn";             fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="sb";             fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="sshn";           fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="sshb";           fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="rhop";           fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="kt";             fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="ndastp";         fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="adatrj";         fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="utau_b";         fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="vtau_b";         fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="qns_b";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="emp_b";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="sfx_b";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="en" ;            fields(i)%grid="grid_N_3D" 
        i = i + 1; fields(i)%vname="avt_k";            fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="avm_k";            fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="dissl";          fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="sbc_hc_b";       fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="sbc_sc_b";       fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="qsr_hc_b";       fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="fraqsr_1lev";    fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="greenland_icesheet_mass"
                                               fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="greenland_icesheet_timelapsed"
                                               fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="greenland_icesheet_mass_roc"
                                               fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="antarctica_icesheet_mass"
                                               fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="antarctica_icesheet_timelapsed"
                                               fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="antarctica_icesheet_mass_roc"
                                               fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="frc_v";          fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="frc_t";          fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="frc_s";          fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="frc_wn_t";       fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="frc_wn_s";       fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="ssh_ini";        fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="e3t_ini";        fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="hc_loc_ini";     fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="sc_loc_ini";     fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="ssh_hc_loc_ini"; fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="ssh_sc_loc_ini"; fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="tilde_e3t_b";    fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="tilde_e3t_n";    fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="hdiv_lf";        fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="ub2_b";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="vb2_b";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="sshbb_e";        fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="ubb_e";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="vbb_e";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="sshb_e";         fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="ub_e";           fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="vb_e";           fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="fwf_isf_b";      fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="isf_sc_b";       fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="isf_hc_b";       fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="ssh_ibb";        fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="rnf_b";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="rnf_hc_b";       fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="rnf_sc_b";       fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="nn_fsbc";        fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="ssu_m";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="ssv_m";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="sst_m";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="sss_m";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="ssh_m";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="e3t_m";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="frq_m";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="avmb";           fields(i)%grid="grid_vector"
        i = i + 1; fields(i)%vname="avtb";           fields(i)%grid="grid_vector"
        i = i + 1; fields(i)%vname="ub2_i_b";        fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="vb2_i_b";        fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="ntime";          fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="Dsst";           fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="tmask";          fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="umask";          fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="vmask";          fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="smask";          fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="gdepw_n";        fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="e3t_n";          fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="e3u_n";          fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="e3v_n";          fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="surf_ini";       fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="e3t_b";          fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="hmxl_n";         fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="un_bf";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="vn_bf";          fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="hbl";            fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="hbli";           fields(i)%grid="grid_N"
        i = i + 1; fields(i)%vname="wn";             fields(i)%grid="grid_N_3D"
        i = i + 1; fields(i)%vname="a_fwb_b";        fields(i)%grid="grid_scalar"
        i = i + 1; fields(i)%vname="a_fwb";          fields(i)%grid="grid_scalar"

        IF( i-1 > max_rst_fields) THEN
           WRITE(ctmp1,*) 'E R R O R : iom_set_rst_vars SIZE of RST_FIELD array is too small'
           CALL ctl_stop( 'iom_set_rst_vars:', ctmp1 )
        ENDIF
   END SUBROUTINE iom_set_rst_vars


   SUBROUTINE iom_set_rstw_active(cdrst_file)
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE iom_set_rstw_active   ***
      !!
      !! ** Purpose :  define file name in XIOS context for writing restart
      !!               enable variables present in restart file for writing
      !!---------------------------------------------------------------------
!sets enabled = .TRUE. for each field in restart file
   CHARACTER(len=*) :: cdrst_file
#if defined key_iomput
   TYPE(xios_field) :: field_hdl
   TYPE(xios_file) :: file_hdl
   TYPE(xios_filegroup) :: filegroup_hdl
   INTEGER :: i
   CHARACTER(lc)  ::   clpath

!set name of the restart file and enable available fields
        IF(lwp) WRITE(numout,*) 'Setting restart filename (for XIOS write) to: ',cdrst_file
        CALL xios_get_handle("file_definition", filegroup_hdl )
        CALL xios_add_child(filegroup_hdl, file_hdl, 'wrestart')
        IF(nxioso.eq.1) THEN 
           CALL xios_set_file_attr( "wrestart", type="one_file", enabled=.TRUE.,& 
                                    mode="write", output_freq=xios_timestep) 
           if(lwp) write(numout,*) 'OPEN ', trim(cdrst_file), ' in one_file mode' 
        ELSE  
           CALL xios_set_file_attr( "wrestart", type="multiple_file", enabled=.TRUE.,& 
                                    mode="write", output_freq=xios_timestep) 
           if(lwp) write(numout,*) 'OPEN ', trim(cdrst_file), ' in multiple_file mode' 
        ENDIF 
        CALL xios_set_file_attr( "wrestart", name=trim(cdrst_file))
!define fields for restart context
        DO i = 1, max_rst_fields
         IF( rst_wfields(i)%active ) THEN
                CALL xios_add_child(file_hdl, field_hdl, TRIM(rst_wfields(i)%vname))
                SELECT CASE (TRIM(rst_wfields(i)%grid))
                 CASE ("grid_N_3D")
                    CALL xios_set_attr (field_hdl, enabled = .TRUE., name = TRIM(rst_wfields(i)%vname), &
                        domain_ref="grid_N", axis_ref="nav_lev", prec = 8, operation = "instant")
                 CASE ("grid_N")
                    CALL xios_set_attr (field_hdl, enabled = .TRUE., name = TRIM(rst_wfields(i)%vname), &
                        domain_ref="grid_N", prec = 8, operation = "instant") 
                 CASE ("grid_vector")
                    CALL xios_set_attr (field_hdl, enabled = .TRUE., name = TRIM(rst_wfields(i)%vname), &
                         axis_ref="nav_lev", prec = 8, operation = "instant")
                 CASE ("grid_scalar")
                    CALL xios_set_attr (field_hdl, enabled = .TRUE., name = TRIM(rst_wfields(i)%vname), &
                        scalar_ref = "grid_scalar", prec = 8, operation = "instant")
                END SELECT
         ENDIF
        END DO
#endif
   END SUBROUTINE iom_set_rstw_active

   SUBROUTINE iom_set_rst_context(ld_rstr) 
     !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_set_rst_context  ***
      !!
      !! ** Purpose : Define domain, axis and grid for restart (read/write) 
      !!              context 
      !!               
      !!---------------------------------------------------------------------
   LOGICAL, INTENT(IN)               :: ld_rstr
!ld_rstr is true for restart context. There is no need to define grid for 
!restart read, because it's read from file
#if defined key_iomput
   TYPE(xios_domaingroup)            :: domaingroup_hdl 
   TYPE(xios_domain)                 :: domain_hdl 
   TYPE(xios_axisgroup)              :: axisgroup_hdl 
   TYPE(xios_axis)                   :: axis_hdl 
   TYPE(xios_scalar)                 :: scalar_hdl 
   TYPE(xios_scalargroup)            :: scalargroup_hdl 

     CALL xios_get_handle("domain_definition",domaingroup_hdl) 
     CALL xios_add_child(domaingroup_hdl, domain_hdl, "grid_N") 
     CALL set_grid("N", glamt, gphit, .TRUE., ld_rstr) 
 
     CALL xios_get_handle("axis_definition",axisgroup_hdl) 
     CALL xios_add_child(axisgroup_hdl, axis_hdl, "nav_lev") 
!AGRIF fails to compile when unit= is in call to xios_set_axis_attr
!    CALL xios_set_axis_attr( "nav_lev", long_name="Vertical levels",  unit="m", positive="down") 
     CALL xios_set_axis_attr( "nav_lev", long_name="Vertical levels in meters", positive="down")
     CALL iom_set_axis_attr( "nav_lev", paxis = gdept_1d ) 

     CALL xios_get_handle("scalar_definition", scalargroup_hdl) 
     CALL xios_add_child(scalargroup_hdl, scalar_hdl, "grid_scalar") 
#endif
   END SUBROUTINE iom_set_rst_context

   SUBROUTINE iom_swap( cdname )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_swap  ***
      !!
      !! ** Purpose :  swap context between different agrif grid for xmlio_server
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in) :: cdname
#if defined key_iomput
      TYPE(xios_context) :: nemo_hdl

      IF( TRIM(Agrif_CFixed()) == '0' ) THEN
        CALL xios_get_handle(TRIM(cdname),nemo_hdl)
      ELSE
        CALL xios_get_handle(TRIM(Agrif_CFixed())//"_"//TRIM(cdname),nemo_hdl)
      ENDIF
      !
      CALL xios_set_current_context(nemo_hdl)
#endif
      !
   END SUBROUTINE iom_swap


   SUBROUTINE iom_open( cdname, kiomid, ldwrt, kdom, ldstop, ldiof, kdlev )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_open  ***
      !!
      !! ** Purpose :  open an input file (return 0 if not found)
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in   )           ::   cdname   ! File name
      INTEGER         , INTENT(  out)           ::   kiomid   ! iom identifier of the opened file
      LOGICAL         , INTENT(in   ), OPTIONAL ::   ldwrt    ! open in write modeb          (default = .FALSE.)
      INTEGER         , INTENT(in   ), OPTIONAL ::   kdom     ! Type of domain to be written (default = jpdom_local_noovlap)
      LOGICAL         , INTENT(in   ), OPTIONAL ::   ldstop   ! stop if open to read a non-existing file (default = .TRUE.)
      LOGICAL         , INTENT(in   ), OPTIONAL ::   ldiof    ! Interp On the Fly, needed for AGRIF (default = .FALSE.)
      INTEGER         , INTENT(in   ), OPTIONAL ::   kdlev    ! number of vertical levels
      !
      CHARACTER(LEN=256)    ::   clname    ! the name of the file based on cdname [[+clcpu]+clcpu]
      CHARACTER(LEN=256)    ::   cltmpn    ! tempory name to store clname (in writting mode)
      CHARACTER(LEN=10)     ::   clsuffix  ! ".nc" 
      CHARACTER(LEN=15)     ::   clcpu     ! the cpu number (max jpmax_digits digits)
      CHARACTER(LEN=256)    ::   clinfo    ! info character
      LOGICAL               ::   llok      ! check the existence 
      LOGICAL               ::   llwrt     ! local definition of ldwrt
      LOGICAL               ::   llnoov    ! local definition to read overlap
      LOGICAL               ::   llstop    ! local definition of ldstop
      LOGICAL               ::   lliof     ! local definition of ldiof
      INTEGER               ::   icnt      ! counter for digits in clcpu (max = jpmax_digits)
      INTEGER               ::   iln, ils  ! lengths of character
      INTEGER               ::   idom      ! type of domain
      INTEGER               ::   istop     ! 
      INTEGER, DIMENSION(2,5) ::   idompar ! domain parameters: 
      ! local number of points for x,y dimensions
      ! position of first local point for x,y dimensions
      ! position of last local point for x,y dimensions
      ! start halo size for x,y dimensions
      ! end halo size for x,y dimensions
      !---------------------------------------------------------------------
      ! Initializations and control
      ! =============
      kiomid = -1
      clinfo = '                    iom_open ~~~  '
      istop = nstop
      ! if iom_open is called for the first time: initialize iom_file(:)%nfid to 0
      ! (could be done when defining iom_file in f95 but not in f90)
      IF( Agrif_Root() ) THEN
         IF( iom_open_init == 0 ) THEN
            iom_file(:)%nfid = 0
            iom_open_init = 1
         ENDIF
      ENDIF
      ! do we read or write the file?
      IF( PRESENT(ldwrt) ) THEN   ;   llwrt = ldwrt
      ELSE                        ;   llwrt = .FALSE.
      ENDIF
      ! do we call ctl_stop if we try to open a non-existing file in read mode?
      IF( PRESENT(ldstop) ) THEN   ;   llstop = ldstop
      ELSE                         ;   llstop = .TRUE.
      ENDIF
      ! are we using interpolation on the fly?
      IF( PRESENT(ldiof) ) THEN   ;   lliof = ldiof
      ELSE                        ;   lliof = .FALSE.
      ENDIF
      ! do we read the overlap 
      ! ugly patch SM+JMM+RB to overwrite global definition in some cases
      llnoov = (jpni * jpnj ) == jpnij .AND. .NOT. lk_agrif
      ! create the file name by added, if needed, TRIM(Agrif_CFixed()) and TRIM(clsuffix)
      ! =============
      clname   = trim(cdname)
      IF ( .NOT. Agrif_Root() .AND. .NOT. lliof ) THEN
!FUS         iln    = INDEX(clname,'/') 
         iln    = INDEX(clname,'/',BACK=.true.)  ! FUS: to insert the nest index at the right location within the string, the last / has to be found (search from the right to left)
         cltmpn = clname(1:iln)
         clname = clname(iln+1:LEN_TRIM(clname))
         clname=TRIM(cltmpn)//TRIM(Agrif_CFixed())//'_'//TRIM(clname)
      ENDIF
      ! which suffix should we use?
      clsuffix = '.nc'
      ! Add the suffix if needed
      iln = LEN_TRIM(clname)
      ils = LEN_TRIM(clsuffix)
      IF( iln <= ils .OR. INDEX( TRIM(clname), TRIM(clsuffix), back = .TRUE. ) /= iln - ils + 1 )   &
         &   clname = TRIM(clname)//TRIM(clsuffix)
      cltmpn = clname   ! store this name
      ! try to find if the file to be opened already exist
      ! =============
      INQUIRE( FILE = clname, EXIST = llok )
      IF( .NOT.llok ) THEN
         ! we try to add the cpu number to the name
         WRITE(clcpu,*) narea-1

         clcpu  = TRIM(ADJUSTL(clcpu))
         iln = INDEX(clname,TRIM(clsuffix), back = .TRUE.)
         clname = clname(1:iln-1)//'_'//TRIM(clcpu)//TRIM(clsuffix)
         icnt = 0
         INQUIRE( FILE = clname, EXIST = llok ) 
         ! we try different formats for the cpu number by adding 0
         DO WHILE( .NOT.llok .AND. icnt < jpmax_digits )
            clcpu  = "0"//trim(clcpu)
            clname = clname(1:iln-1)//'_'//TRIM(clcpu)//TRIM(clsuffix)
            INQUIRE( FILE = clname, EXIST = llok )
            icnt = icnt + 1
         END DO
      ELSE
         lxios_sini = .TRUE.
      ENDIF
      IF( llwrt ) THEN
         ! check the domain definition
! JMM + SM: ugly patch before getting the new version of lib_mpp)
!         idom = jpdom_local_noovlap   ! default definition
         IF( llnoov ) THEN   ;   idom = jpdom_local_noovlap   ! default definition
         ELSE                ;   idom = jpdom_local_full      ! default definition
         ENDIF
         IF( PRESENT(kdom) )   idom = kdom
         ! create the domain informations
         ! =============
         SELECT CASE (idom)
         CASE (jpdom_local_full)
            idompar(:,1) = (/ jpi             , jpj              /)
            idompar(:,2) = (/ nimpp           , njmpp            /)
            idompar(:,3) = (/ nimpp + jpi - 1 , njmpp + jpj - 1  /)
            idompar(:,4) = (/ nldi - 1        , nldj - 1         /)
            idompar(:,5) = (/ jpi - nlei      , jpj - nlej       /)
         CASE (jpdom_local_noextra)
            idompar(:,1) = (/ nlci            , nlcj             /)
            idompar(:,2) = (/ nimpp           , njmpp            /)
            idompar(:,3) = (/ nimpp + nlci - 1, njmpp + nlcj - 1 /)
            idompar(:,4) = (/ nldi - 1        , nldj - 1         /)
            idompar(:,5) = (/ nlci - nlei     , nlcj - nlej      /)
         CASE (jpdom_local_noovlap)
            idompar(:,1) = (/ nlei  - nldi + 1, nlej  - nldj + 1 /)
            idompar(:,2) = (/ nimpp + nldi - 1, njmpp + nldj - 1 /)
            idompar(:,3) = (/ nimpp + nlei - 1, njmpp + nlej - 1 /)
            idompar(:,4) = (/ 0               , 0                /)
            idompar(:,5) = (/ 0               , 0                /)
         CASE DEFAULT
            CALL ctl_stop( TRIM(clinfo), 'wrong value of kdom, only jpdom_local* cases are accepted' )
         END SELECT
      ENDIF
      ! Open the NetCDF file
      ! =============
      ! do we have some free file identifier?
      IF( MINVAL(iom_file(:)%nfid) /= 0 )   &
         &   CALL ctl_stop( TRIM(clinfo), 'No more free file identifier', 'increase jpmax_files in iom_def' )
      ! if no file was found...
      IF( .NOT. llok ) THEN
         IF( .NOT. llwrt ) THEN   ! we are in read mode 
            IF( llstop ) THEN   ;   CALL ctl_stop( TRIM(clinfo), 'File '//TRIM(cltmpn)//'* not found' )
            ELSE                ;   istop = nstop + 1   ! make sure that istop /= nstop so we don't open the file
            ENDIF
         ELSE                     ! we are in write mode so we 
            clname = cltmpn       ! get back the file name without the cpu number
         ENDIF
      ELSE
         IF( llwrt .AND. .NOT. ln_clobber ) THEN   ! we stop as we want to write in a new file 
            CALL ctl_stop( TRIM(clinfo), 'We want to write in a new file but '//TRIM(clname)//' already exists...' )
            istop = nstop + 1                      ! make sure that istop /= nstop so we don't open the file
         ELSEIF( llwrt ) THEN     ! the file exists and we are in write mode with permission to 
            clname = cltmpn       ! overwrite so get back the file name without the cpu number
         ENDIF
      ENDIF
      IF( istop == nstop ) THEN   ! no error within this routine
         CALL iom_nf90_open( clname, kiomid, llwrt, llok, idompar, kdlev = kdlev )
      ENDIF
      !
   END SUBROUTINE iom_open


   SUBROUTINE iom_close( kiomid )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_close  ***
      !!
      !! ** Purpose : close an input file, or all files opened by iom
      !!--------------------------------------------------------------------
      INTEGER, INTENT(inout), OPTIONAL ::   kiomid   ! iom identifier of the file to be closed
      !                                              ! return 0 when file is properly closed
      !                                              ! No argument: all files opened by iom are closed

      INTEGER ::   jf         ! dummy loop indices
      INTEGER ::   i_s, i_e   ! temporary integer
      CHARACTER(LEN=100)    ::   clinfo    ! info character
      !---------------------------------------------------------------------
      !
      IF( iom_open_init == 0 )   RETURN   ! avoid to use iom_file(jf)%nfid that us not yet initialized
      !
      clinfo = '                    iom_close ~~~  '
      IF( PRESENT(kiomid) ) THEN
         i_s = kiomid
         i_e = kiomid
      ELSE
         i_s = 1
         i_e = jpmax_files
      ENDIF

      IF( i_s > 0 ) THEN
         DO jf = i_s, i_e
            IF( iom_file(jf)%nfid > 0 ) THEN
               CALL iom_nf90_close( jf )
               iom_file(jf)%nfid       = 0          ! free the id 
               IF( PRESENT(kiomid) )   kiomid = 0   ! return 0 as id to specify that the file was closed
               IF(lwp) WRITE(numout,*) TRIM(clinfo)//' close file: '//TRIM(iom_file(jf)%name)//' ok'
            ELSEIF( PRESENT(kiomid) ) THEN
               WRITE(ctmp1,*) '--->',  kiomid
               CALL ctl_stop( TRIM(clinfo)//' Invalid file identifier', ctmp1 )
            ENDIF
         END DO
      ENDIF
      !    
   END SUBROUTINE iom_close


   FUNCTION iom_varid ( kiomid, cdvar, kdimsz, kndims, lduld, ldstop )  
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  iom_varid  ***
      !!
      !! ** Purpose : get the id of a variable in a file (return 0 if not found)
      !!-----------------------------------------------------------------------
      INTEGER              , INTENT(in   )           ::   kiomid   ! file Identifier
      CHARACTER(len=*)     , INTENT(in   )           ::   cdvar    ! name of the variable
      INTEGER, DIMENSION(:), INTENT(  out), OPTIONAL ::   kdimsz   ! size of each dimension
      INTEGER              , INTENT(  out), OPTIONAL ::   kndims   ! number of dimensions
      LOGICAL              , INTENT(  out), OPTIONAL ::   lduld    ! true if the last dimension is unlimited (time)
      LOGICAL              , INTENT(in   ), OPTIONAL ::   ldstop   ! stop if looking for non-existing variable (default = .TRUE.)
      !
      INTEGER                        ::   iom_varid, iiv, i_nvd
      LOGICAL                        ::   ll_fnd
      CHARACTER(LEN=100)             ::   clinfo                   ! info character
      LOGICAL                        ::   llstop                   ! local definition of ldstop
      !!-----------------------------------------------------------------------
      iom_varid = 0                         ! default definition
      ! do we call ctl_stop if we look for non-existing variable?
      IF( PRESENT(ldstop) ) THEN   ;   llstop = ldstop
      ELSE                         ;   llstop = .TRUE.
      ENDIF
      !
      IF( kiomid > 0 ) THEN
         clinfo = 'iom_varid, file: '//trim(iom_file(kiomid)%name)//', var: '//trim(cdvar)
         IF( iom_file(kiomid)%nfid == 0 ) THEN 
            CALL ctl_stop( trim(clinfo), 'the file is not open' )
         ELSE
            ll_fnd  = .FALSE.
            iiv = 0
            !
            DO WHILE ( .NOT.ll_fnd .AND. iiv < iom_file(kiomid)%nvars )
               iiv = iiv + 1
               ll_fnd  = ( TRIM(cdvar) == TRIM(iom_file(kiomid)%cn_var(iiv)) )
            END DO
            !
            IF( .NOT.ll_fnd ) THEN
               iiv = iiv + 1
               IF( iiv <= jpmax_vars ) THEN
                  iom_varid = iom_nf90_varid( kiomid, cdvar, iiv, kdimsz, kndims, lduld )
               ELSE
                  CALL ctl_stop( trim(clinfo), 'Too many variables in the file '//iom_file(kiomid)%name,   &
                        &                      'increase the parameter jpmax_vars')
               ENDIF
               IF( llstop .AND. iom_varid == -1 )   CALL ctl_stop( TRIM(clinfo)//' not found' ) 
            ELSE
               iom_varid = iiv
               IF( PRESENT(kdimsz) ) THEN 
                  i_nvd = iom_file(kiomid)%ndims(iiv)
                  IF( i_nvd <= size(kdimsz) ) THEN
                     kdimsz(1:i_nvd) = iom_file(kiomid)%dimsz(1:i_nvd,iiv)
                  ELSE
                     WRITE(ctmp1,*) i_nvd, size(kdimsz)
                     CALL ctl_stop( trim(clinfo), 'error in kdimsz size'//trim(ctmp1) )
                  ENDIF
               ENDIF
               IF( PRESENT(kndims) )  kndims = iom_file(kiomid)%ndims(iiv)
               IF( PRESENT( lduld) )  lduld  = iom_file(kiomid)%luld( iiv)
            ENDIF
         ENDIF
      ENDIF
      !
   END FUNCTION iom_varid


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_get
   !!----------------------------------------------------------------------
   SUBROUTINE iom_g0d( kiomid, cdvar, pvar, ktime, ldxios )
      INTEGER         , INTENT(in   )                 ::   kiomid    ! Identifier of the file
      CHARACTER(len=*), INTENT(in   )                 ::   cdvar     ! Name of the variable
      REAL(wp)        , INTENT(  out)                 ::   pvar      ! read field
      INTEGER         , INTENT(in   ),     OPTIONAL   ::   ktime     ! record number
      LOGICAL         , INTENT(in   ),     OPTIONAL   ::   ldxios    ! use xios to read restart
      !
      INTEGER                                         ::   idvar     ! variable id
      INTEGER                                         ::   idmspc    ! number of spatial dimensions
      INTEGER         , DIMENSION(1)                  ::   itime     ! record number
      CHARACTER(LEN=100)                              ::   clinfo    ! info character
      CHARACTER(LEN=100)                              ::   clname    ! file name
      CHARACTER(LEN=1)                                ::   cldmspc   !
      LOGICAL                                         ::   llxios
      !
      llxios = .FALSE.
      IF( PRESENT(ldxios) ) llxios = ldxios

      IF(.NOT.llxios) THEN  ! read data using default library
         itime = 1
         IF( PRESENT(ktime) ) itime = ktime
         !
         clname = iom_file(kiomid)%name
         clinfo = '          iom_g0d, file: '//trim(clname)//', var: '//trim(cdvar)
         !
         IF( kiomid > 0 ) THEN
            idvar = iom_varid( kiomid, cdvar )
            IF( iom_file(kiomid)%nfid > 0 .AND. idvar > 0 ) THEN
               idmspc = iom_file ( kiomid )%ndims( idvar )
               IF( iom_file(kiomid)%luld(idvar) )  idmspc = idmspc - 1
               WRITE(cldmspc , fmt='(i1)') idmspc
               IF( idmspc > 0 )  CALL ctl_stop( TRIM(clinfo), 'When reading to a 0D array, we do not accept data', &
                                    &                         'with 1 or more spatial dimensions: '//cldmspc//' were found.' , &
                                    &                         'Use ncwa -a to suppress the unnecessary dimensions' )
               CALL iom_nf90_get( kiomid, idvar, pvar, itime )
            ENDIF
         ENDIF
      ELSE
#if defined key_iomput
         IF(lwp) WRITE(numout,*) 'XIOS RST READ (0D): ', trim(cdvar)
         CALL iom_swap( TRIM(crxios_context) )
         CALL xios_recv_field( trim(cdvar), pvar)
         CALL iom_swap( TRIM(cxios_context) )
#else
         WRITE(ctmp1,*) 'Can not use XIOS in iom_g0d, file: '//trim(clname)//', var:'//trim(cdvar)
         CALL ctl_stop( 'iom_g0d', ctmp1 )
#endif
      ENDIF
   END SUBROUTINE iom_g0d

   SUBROUTINE iom_g1d( kiomid, kdom, cdvar, pvar, ktime, kstart, kcount, ldxios )
      INTEGER         , INTENT(in   )                         ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                         ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                         ::   cdvar     ! Name of the variable
      REAL(wp)        , INTENT(  out), DIMENSION(:)           ::   pvar      ! read field
      INTEGER         , INTENT(in   )              , OPTIONAL ::   ktime     ! record number
      INTEGER         , INTENT(in   ), DIMENSION(1), OPTIONAL ::   kstart    ! start axis position of the reading 
      INTEGER         , INTENT(in   ), DIMENSION(1), OPTIONAL ::   kcount    ! number of points in each axis
      LOGICAL         , INTENT(in   ),               OPTIONAL ::   ldxios    ! read data using XIOS
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom       , cdvar        , pv_r1d=pvar,   &
              &                                                     ktime=ktime, kstart=kstart, kcount=kcount, &
              &                                                     ldxios=ldxios )
      ENDIF
   END SUBROUTINE iom_g1d

   SUBROUTINE iom_g2d( kiomid, kdom, cdvar, pvar, ktime, kstart, kcount, lrowattr, ldxios)
      INTEGER         , INTENT(in   )                           ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                           ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                           ::   cdvar     ! Name of the variable
      REAL(wp)        , INTENT(  out), DIMENSION(:,:)           ::   pvar      ! read field
      INTEGER         , INTENT(in   )                , OPTIONAL ::   ktime     ! record number
      INTEGER         , INTENT(in   ), DIMENSION(2)  , OPTIONAL ::   kstart    ! start axis position of the reading 
      INTEGER         , INTENT(in   ), DIMENSION(2)  , OPTIONAL ::   kcount    ! number of points in each axis
      LOGICAL         , INTENT(in   )                , OPTIONAL ::   lrowattr  ! logical flag telling iom_get to
                                                                               ! look for and use a file attribute
                                                                               ! called open_ocean_jstart to set the start
                                                                               ! value for the 2nd dimension (netcdf only)
      LOGICAL         , INTENT(in   ),               OPTIONAL ::   ldxios      ! read data using XIOS
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom       , cdvar        , pv_r2d=pvar,   &
              &                                                     ktime=ktime, kstart=kstart, kcount=kcount, &
              &                                                     lrowattr=lrowattr,  ldxios=ldxios)
      ENDIF
   END SUBROUTINE iom_g2d

   SUBROUTINE iom_g3d( kiomid, kdom, cdvar, pvar, ktime, kstart, kcount, lrowattr, ldxios )
      INTEGER         , INTENT(in   )                             ::   kiomid    ! Identifier of the file
      INTEGER         , INTENT(in   )                             ::   kdom      ! Type of domain to be read
      CHARACTER(len=*), INTENT(in   )                             ::   cdvar     ! Name of the variable
      REAL(wp)        , INTENT(  out), DIMENSION(:,:,:)           ::   pvar      ! read field
      INTEGER         , INTENT(in   )                  , OPTIONAL ::   ktime     ! record number
      INTEGER         , INTENT(in   ), DIMENSION(3)    , OPTIONAL ::   kstart    ! start axis position of the reading 
      INTEGER         , INTENT(in   ), DIMENSION(3)    , OPTIONAL ::   kcount    ! number of points in each axis
      LOGICAL         , INTENT(in   )                  , OPTIONAL ::   lrowattr  ! logical flag telling iom_get to
                                                                                 ! look for and use a file attribute
                                                                                 ! called open_ocean_jstart to set the start
                                                                                 ! value for the 2nd dimension (netcdf only)
      LOGICAL         , INTENT(in   ),               OPTIONAL ::   ldxios        ! read data using XIOS
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 ) CALL iom_get_123d( kiomid, kdom       , cdvar        , pv_r3d=pvar,   &
              &                                                     ktime=ktime, kstart=kstart, kcount=kcount, &
              &                                                     lrowattr=lrowattr, ldxios=ldxios )
      ENDIF
   END SUBROUTINE iom_g3d
   !!----------------------------------------------------------------------

   SUBROUTINE iom_get_123d( kiomid, kdom  , cdvar ,   &
         &                  pv_r1d, pv_r2d, pv_r3d,   &
         &                  ktime , kstart, kcount,   &
         &                  lrowattr, ldxios        )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_get_123d  ***
      !!
      !! ** Purpose : read a 1D/2D/3D variable
      !!
      !! ** Method : read ONE record at each CALL
      !!-----------------------------------------------------------------------
      INTEGER                    , INTENT(in   )           ::   kiomid     ! Identifier of the file
      INTEGER                    , INTENT(in   )           ::   kdom       ! Type of domain to be read
      CHARACTER(len=*)           , INTENT(in   )           ::   cdvar      ! Name of the variable
      REAL(wp), DIMENSION(:)     , INTENT(  out), OPTIONAL ::   pv_r1d     ! read field (1D case)
      REAL(wp), DIMENSION(:,:)   , INTENT(  out), OPTIONAL ::   pv_r2d     ! read field (2D case)
      REAL(wp), DIMENSION(:,:,:) , INTENT(  out), OPTIONAL ::   pv_r3d     ! read field (3D case)
      INTEGER                    , INTENT(in   ), OPTIONAL ::   ktime      ! record number
      INTEGER , DIMENSION(:)     , INTENT(in   ), OPTIONAL ::   kstart     ! start position of the reading in each axis 
      INTEGER , DIMENSION(:)     , INTENT(in   ), OPTIONAL ::   kcount     ! number of points to be read in each axis
      LOGICAL                    , INTENT(in   ), OPTIONAL ::   lrowattr   ! logical flag telling iom_get to
                                                                           ! look for and use a file attribute
                                                                           ! called open_ocean_jstart to set the start
                                                                           ! value for the 2nd dimension (netcdf only)
      LOGICAL                    , INTENT(in   ), OPTIONAL ::   ldxios     ! use XIOS to read restart
      !
      LOGICAL                        ::   llxios       ! local definition for XIOS read
      LOGICAL                        ::   llnoov      ! local definition to read overlap
      LOGICAL                        ::   luse_jattr  ! local definition to read open_ocean_jstart file attribute
      INTEGER                        ::   jstartrow   ! start point for 2nd dimension optionally set by file attribute
      INTEGER                        ::   jl          ! loop on number of dimension 
      INTEGER                        ::   idom        ! type of domain
      INTEGER                        ::   idvar       ! id of the variable
      INTEGER                        ::   inbdim      ! number of dimensions of the variable
      INTEGER                        ::   idmspc      ! number of spatial dimensions 
      INTEGER                        ::   itime       ! record number
      INTEGER                        ::   istop       ! temporary value of nstop
      INTEGER                        ::   ix1, ix2, iy1, iy2   ! subdomain indexes
      INTEGER                        ::   ji, jj      ! loop counters
      INTEGER                        ::   irankpv     ! 
      INTEGER                        ::   ind1, ind2  ! substring index
      INTEGER, DIMENSION(jpmax_dims) ::   istart      ! starting point to read for each axis
      INTEGER, DIMENSION(jpmax_dims) ::   icnt        ! number of value to read along each axis 
      INTEGER, DIMENSION(jpmax_dims) ::   idimsz      ! size of the dimensions of the variable
      INTEGER, DIMENSION(jpmax_dims) ::   ishape      ! size of the dimensions of the variable
      REAL(wp)                       ::   zscf, zofs  ! sacle_factor and add_offset
      INTEGER                        ::   itmp        ! temporary integer
      CHARACTER(LEN=256)             ::   clinfo      ! info character
      CHARACTER(LEN=256)             ::   clname      ! file name
      CHARACTER(LEN=1)               ::   clrankpv, cldmspc      ! 
      LOGICAL                        ::   ll_depth_spec ! T => if kstart, kcount present then *only* use values for 3rd spatial dimension.
      INTEGER                        ::   inlev       ! number of levels for 3D data
      REAL(wp)                       ::   gma, gmi
      !---------------------------------------------------------------------
      !
      inlev = -1
      IF( PRESENT(pv_r3d) )   inlev = SIZE(pv_r3d, 3)
      !
      llxios = .FALSE.
      if(PRESENT(ldxios)) llxios = ldxios
      idvar = iom_varid( kiomid, cdvar ) 
      idom = kdom
      !
      IF(.NOT.llxios) THEN
         clname = iom_file(kiomid)%name   !   esier to read
         clinfo = '          iom_get_123d, file: '//trim(clname)//', var: '//trim(cdvar)
         ! local definition of the domain ?
         ! do we read the overlap 
         ! ugly patch SM+JMM+RB to overwrite global definition in some cases
         llnoov = (jpni * jpnj ) == jpnij .AND. .NOT. lk_agrif 
         ! check kcount and kstart optionals parameters...
         IF( PRESENT(kcount) .AND. (.NOT. PRESENT(kstart)) ) CALL ctl_stop(trim(clinfo), 'kcount present needs kstart present')
         IF( PRESENT(kstart) .AND. (.NOT. PRESENT(kcount)) ) CALL ctl_stop(trim(clinfo), 'kstart present needs kcount present')
         IF( PRESENT(kstart) .AND. idom /= jpdom_unknown .AND.  idom /= jpdom_autoglo_xy  ) &
     &          CALL ctl_stop(trim(clinfo), 'kstart present needs kdom = jpdom_unknown or kdom = jpdom_autoglo_xy')

         luse_jattr = .false.
         IF( PRESENT(lrowattr) ) THEN
            IF( lrowattr .AND. idom /= jpdom_data   ) CALL ctl_stop(trim(clinfo), 'lrowattr present and true needs kdom = jpdom_data')
            IF( lrowattr .AND. idom == jpdom_data   ) luse_jattr = .true.
         ENDIF

         ! Search for the variable in the data base (eventually actualize data)
         istop = nstop
         !
         IF( idvar > 0 ) THEN
            ! to write iom_file(kiomid)%dimsz in a shorter way !
            idimsz(:) = iom_file(kiomid)%dimsz(:, idvar) 
            inbdim = iom_file(kiomid)%ndims(idvar)            ! number of dimensions in the file
            idmspc = inbdim                                   ! number of spatial dimensions in the file
            IF( iom_file(kiomid)%luld(idvar) )   idmspc = inbdim - 1
            IF( idmspc > 3 )   CALL ctl_stop(trim(clinfo), 'the file has more than 3 spatial dimensions this case is not coded...') 
            !
            ! update idom definition...
            ! Identify the domain in case of jpdom_auto(glo/dta) definition
            IF( idom == jpdom_autoglo_xy ) THEN
               ll_depth_spec = .TRUE.
               idom = jpdom_autoglo
            ELSE
               ll_depth_spec = .FALSE.
            ENDIF
            IF( idom == jpdom_autoglo .OR. idom == jpdom_autodta ) THEN            
               IF( idom == jpdom_autoglo ) THEN   ;   idom = jpdom_global 
               ELSE                               ;   idom = jpdom_data
               ENDIF
               ind1 = INDEX( clname, '_', back = .TRUE. ) + 1
               ind2 = INDEX( clname, '.', back = .TRUE. ) - 1
               IF( ind2 > ind1 ) THEN   ;   IF( VERIFY( clname(ind1:ind2), '0123456789' ) == 0 )   idom = jpdom_local   ;   ENDIF
            ENDIF
            ! Identify the domain in case of jpdom_local definition
            IF( idom == jpdom_local ) THEN
               IF(     idimsz(1) == jpi               .AND. idimsz(2) == jpj               ) THEN   ;   idom = jpdom_local_full
               ELSEIF( idimsz(1) == nlci              .AND. idimsz(2) == nlcj              ) THEN   ;   idom = jpdom_local_noextra
               ELSEIF( idimsz(1) == (nlei - nldi + 1) .AND. idimsz(2) == (nlej - nldj + 1) ) THEN   ;   idom = jpdom_local_noovlap
               ELSE   ;   CALL ctl_stop( trim(clinfo), 'impossible to identify the local domain' )
               ENDIF
            ENDIF
            !
            ! check the consistency between input array and data rank in the file
            !
            ! initializations
            itime = 1
            IF( PRESENT(ktime) ) itime = ktime
            !
            irankpv = 1 * COUNT( (/PRESENT(pv_r1d)/) ) + 2 * COUNT( (/PRESENT(pv_r2d)/) ) + 3 * COUNT( (/PRESENT(pv_r3d)/) )
            WRITE(clrankpv, fmt='(i1)') irankpv
            WRITE(cldmspc , fmt='(i1)') idmspc
            !
            IF(     idmspc <  irankpv ) THEN 
               CALL ctl_stop( TRIM(clinfo), 'The file has only '//cldmspc//' spatial dimension',   &
                  &                         'it is impossible to read a '//clrankpv//'D array from this file...' )
            ELSEIF( idmspc == irankpv ) THEN
               IF( PRESENT(pv_r1d) .AND. idom /= jpdom_unknown )   &
                  &   CALL ctl_stop( TRIM(clinfo), 'case not coded...You must use jpdom_unknown' )
            ELSEIF( idmspc >  irankpv ) THEN
                  IF( PRESENT(pv_r2d) .AND. itime == 1 .AND. idimsz(3) == 1 .AND. idmspc == 3 ) THEN
                     CALL ctl_warn( trim(clinfo), '2D array but 3 spatial dimensions for the data...'              ,   &
                           &         'As the size of the z dimension is 1 and as we try to read the first record, ',   &
                           &         'we accept this case, even if there is a possible mix-up between z and time dimension' )   
                     idmspc = idmspc - 1
                  ELSE
                     CALL ctl_stop( TRIM(clinfo), 'To keep iom lisibility, when reading a '//clrankpv//'D array,'         ,   &
                        &                         'we do not accept data with '//cldmspc//' spatial dimensions',   &
                        &                         'Use ncwa -a to suppress the unnecessary dimensions' )
                  ENDIF
            ENDIF
            !
            ! definition of istart and icnt
            !
            icnt  (:) = 1
            istart(:) = 1
            istart(idmspc+1) = itime
   
            IF( PRESENT(kstart) .AND. .NOT. ll_depth_spec ) THEN 
               istart(1:idmspc) = kstart(1:idmspc) 
               icnt  (1:idmspc) = kcount(1:idmspc)
            ELSE
               IF(idom == jpdom_unknown ) THEN
                  icnt(1:idmspc) = idimsz(1:idmspc)
               ELSE 
                  IF( .NOT. PRESENT(pv_r1d) ) THEN   !   not a 1D array
                     IF(     idom == jpdom_data    ) THEN
                        jstartrow = 1
                        IF( luse_jattr ) THEN
                           CALL iom_getatt(kiomid, 'open_ocean_jstart', jstartrow ) ! -999 is returned if the attribute is not found
                           jstartrow = MAX(1,jstartrow)
                        ENDIF
                        istart(1:2) = (/ mig(1), mjg(1) + jstartrow - 1 /)  ! icnt(1:2) done below
                     ELSEIF( idom == jpdom_global  ) THEN ; istart(1:2) = (/ nimpp , njmpp  /)  ! icnt(1:2) done below
                     ENDIF
                     ! we do not read the overlap                     -> we start to read at nldi, nldj
! JMM + SM: ugly patch before getting the new version of lib_mpp)
!                  IF( idom /= jpdom_local_noovlap )   istart(1:2) = istart(1:2) + (/ nldi - 1, nldj - 1 /)
                     IF( llnoov .AND. idom /= jpdom_local_noovlap ) istart(1:2) = istart(1:2) + (/ nldi - 1, nldj - 1 /)
                  ! we do not read the overlap and the extra-halos -> from nldi to nlei and from nldj to nlej 
! JMM + SM: ugly patch before getting the new version of lib_mpp)
!                  icnt(1:2) = (/ nlei - nldi + 1, nlej - nldj + 1 /)
                     IF( llnoov ) THEN   ;   icnt(1:2) = (/ nlei - nldi + 1, nlej - nldj + 1 /)
                     ELSE                ;   icnt(1:2) = (/ nlci           , nlcj            /)
                     ENDIF
                     IF( PRESENT(pv_r3d) ) THEN
                        IF( idom == jpdom_data ) THEN                        ;                               icnt(3) = inlev
                        ELSEIF( ll_depth_spec .AND. PRESENT(kstart) ) THEN   ;   istart(3) = kstart(3)   ;   icnt(3) = kcount(3)
                        ELSE                                                 ;                               icnt(3) = inlev
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF

            ! check that istart and icnt can be used with this file
            !-
            DO jl = 1, jpmax_dims
               itmp = istart(jl)+icnt(jl)-1
               IF( itmp > idimsz(jl) .AND. idimsz(jl) /= 0 ) THEN
                  WRITE( ctmp1, FMT="('(istart(', i1, ') + icnt(', i1, ') - 1) = ', i5)" ) jl, jl, itmp
                  WRITE( ctmp2, FMT="(' is larger than idimsz(', i1,') = ', i5)"         ) jl, idimsz(jl)
                  CALL ctl_stop( trim(clinfo), 'start and count too big regarding to the size of the data, ', ctmp1, ctmp2 )     
               ENDIF
            END DO

            ! check that icnt matches the input array
            !-     
            IF( idom == jpdom_unknown ) THEN
               IF( irankpv == 1 )        ishape(1:1) = SHAPE(pv_r1d)
               IF( irankpv == 2 )        ishape(1:2) = SHAPE(pv_r2d)
               IF( irankpv == 3 )        ishape(1:3) = SHAPE(pv_r3d)
               ctmp1 = 'd'
            ELSE
               IF( irankpv == 2 ) THEN
! JMM + SM: ugly patch before getting the new version of lib_mpp)
!               ishape(1:2) = SHAPE(pv_r2d(nldi:nlei,nldj:nlej  ))   ;   ctmp1 = 'd(nldi:nlei,nldj:nlej)'
                  IF( llnoov ) THEN ; ishape(1:2)=SHAPE(pv_r2d(nldi:nlei,nldj:nlej  )) ; ctmp1='d(nldi:nlei,nldj:nlej)'
                  ELSE              ; ishape(1:2)=SHAPE(pv_r2d(1   :nlci,1   :nlcj  )) ; ctmp1='d(1:nlci,1:nlcj)'
                  ENDIF
               ENDIF
               IF( irankpv == 3 ) THEN 
! JMM + SM: ugly patch before getting the new version of lib_mpp)
!               ishape(1:3) = SHAPE(pv_r3d(nldi:nlei,nldj:nlej,:))   ;   ctmp1 = 'd(nldi:nlei,nldj:nlej,:)'
                  IF( llnoov ) THEN ; ishape(1:3)=SHAPE(pv_r3d(nldi:nlei,nldj:nlej,:)) ; ctmp1='d(nldi:nlei,nldj:nlej,:)'
                  ELSE              ; ishape(1:3)=SHAPE(pv_r3d(1   :nlci,1   :nlcj,:)) ; ctmp1='d(1:nlci,1:nlcj,:)'
                  ENDIF
               ENDIF
            ENDIF
         
            DO jl = 1, irankpv
               WRITE( ctmp2, FMT="(', ', i1,'): ', i5,' /= icnt(', i1,'):', i5)" ) jl, ishape(jl), jl, icnt(jl)
               IF( ishape(jl) /= icnt(jl) )   CALL ctl_stop( TRIM(clinfo), 'size(pv_r'//clrankpv//TRIM(ctmp1)//TRIM(ctmp2) )
            END DO

         ENDIF

         ! read the data
         !-     
         IF( idvar > 0 .AND. istop == nstop ) THEN   ! no additional errors until this point...
            !
         ! find the right index of the array to be read
! JMM + SM: ugly patch before getting the new version of lib_mpp)
!         IF( idom /= jpdom_unknown ) THEN   ;   ix1 = nldi   ;   ix2 = nlei      ;   iy1 = nldj   ;   iy2 = nlej
!         ELSE                               ;   ix1 = 1      ;   ix2 = icnt(1)   ;   iy1 = 1      ;   iy2 = icnt(2)
!         ENDIF
            IF( llnoov ) THEN
               IF( idom /= jpdom_unknown ) THEN   ;   ix1 = nldi   ;   ix2 = nlei      ;   iy1 = nldj   ;   iy2 = nlej
               ELSE                               ;   ix1 = 1      ;   ix2 = icnt(1)   ;   iy1 = 1      ;   iy2 = icnt(2)
               ENDIF
            ELSE
               IF( idom /= jpdom_unknown ) THEN   ;   ix1 = 1      ;   ix2 = nlci      ;   iy1 = 1      ;   iy2 = nlcj
               ELSE                               ;   ix1 = 1      ;   ix2 = icnt(1)   ;   iy1 = 1      ;   iy2 = icnt(2)
               ENDIF
            ENDIF
      
            CALL iom_nf90_get( kiomid, idvar, inbdim, istart, icnt, ix1, ix2, iy1, iy2, pv_r1d, pv_r2d, pv_r3d )

            IF( istop == nstop ) THEN   ! no additional errors until this point...
               IF(lwp) WRITE(numout,"(10x,' read ',a,' (rec: ',i6,') in ',a,' ok')") TRIM(cdvar), itime, TRIM(iom_file(kiomid)%name)
             
               !--- overlap areas and extra hallows (mpp)
               IF(     PRESENT(pv_r2d) .AND. idom /= jpdom_unknown ) THEN
                  CALL lbc_lnk( 'iom', pv_r2d,'Z', -999., kfillmode = jpfillnothing )
               ELSEIF( PRESENT(pv_r3d) .AND. idom /= jpdom_unknown ) THEN
                  ! this if could be simplified with the new lbc_lnk that works with any size of the 3rd dimension
                  IF( icnt(3) == inlev ) THEN
                     CALL lbc_lnk( 'iom', pv_r3d,'Z', -999., kfillmode = jpfillnothing )
                  ELSE   ! put some arbitrary value (a call to lbc_lnk will be done later...)
                     DO jj = nlcj+1, jpj   ;   pv_r3d(1:nlci, jj, :) = pv_r3d(1:nlci, nlej, :)   ;   END DO
                     DO ji = nlci+1, jpi   ;   pv_r3d(ji    , : , :) = pv_r3d(nlei  , :   , :)   ;   END DO
                  ENDIF
               ENDIF
               !
            ELSE
               ! return if istop == nstop is false
               RETURN
            ENDIF
         ELSE
            ! return if statment idvar > 0 .AND. istop == nstop is false
            RETURN
         ENDIF
         !
      ELSE        ! read using XIOS. Only if KEY_IOMPUT is defined
#if defined key_iomput
!would be good to be able to check which context is active and swap only if current is not restart
         CALL iom_swap( TRIM(crxios_context) ) 
         IF( PRESENT(pv_r3d) ) THEN
            pv_r3d(:, :, :) = 0.
            if(lwp) write(numout,*) 'XIOS RST READ (3D): ',trim(cdvar)
            CALL xios_recv_field( trim(cdvar), pv_r3d)
            IF(idom /= jpdom_unknown ) then
                CALL lbc_lnk( 'iom', pv_r3d,'Z', -999., kfillmode = jpfillnothing)
            ENDIF
         ELSEIF( PRESENT(pv_r2d) ) THEN
            pv_r2d(:, :) = 0.
            if(lwp) write(numout,*) 'XIOS RST READ (2D): ', trim(cdvar)
            CALL xios_recv_field( trim(cdvar), pv_r2d)
            IF(idom /= jpdom_unknown ) THEN
                CALL lbc_lnk('iom', pv_r2d,'Z',-999., kfillmode = jpfillnothing)
            ENDIF
         ELSEIF( PRESENT(pv_r1d) ) THEN
            pv_r1d(:) = 0.
            if(lwp) write(numout,*) 'XIOS RST READ (1D): ', trim(cdvar)
            CALL xios_recv_field( trim(cdvar), pv_r1d)
         ENDIF
         CALL iom_swap( TRIM(cxios_context) )
#else
         istop = istop + 1 
         clinfo = 'Can not use XIOS in iom_get_123d, file: '//trim(clname)//', var:'//trim(cdvar)
#endif
      ENDIF
!some final adjustments
      ! C1D case : always call lbc_lnk to replicate the central value over the whole 3X3 domain
      IF( lk_c1d .AND. PRESENT(pv_r2d) )   CALL lbc_lnk( 'iom', pv_r2d,'Z',1. )
      IF( lk_c1d .AND. PRESENT(pv_r3d) )   CALL lbc_lnk( 'iom', pv_r3d,'Z',1. )

      !--- Apply scale_factor and offset
      zscf = iom_file(kiomid)%scf(idvar)      ! scale factor
      zofs = iom_file(kiomid)%ofs(idvar)      ! offset
      IF(     PRESENT(pv_r1d) ) THEN
         IF( zscf /= 1. )   pv_r1d(:) = pv_r1d(:) * zscf 
         IF( zofs /= 0. )   pv_r1d(:) = pv_r1d(:) + zofs
      ELSEIF( PRESENT(pv_r2d) ) THEN
         IF( zscf /= 1.)   pv_r2d(:,:) = pv_r2d(:,:) * zscf
         IF( zofs /= 0.)   pv_r2d(:,:) = pv_r2d(:,:) + zofs
      ELSEIF( PRESENT(pv_r3d) ) THEN
         IF( zscf /= 1.)   pv_r3d(:,:,:) = pv_r3d(:,:,:) * zscf
         IF( zofs /= 0.)   pv_r3d(:,:,:) = pv_r3d(:,:,:) + zofs
      ENDIF
      !
   END SUBROUTINE iom_get_123d

   SUBROUTINE iom_get_var( cdname, z2d)
      CHARACTER(LEN=*), INTENT(in ) ::   cdname
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d 
#if defined key_iomput
      IF( xios_field_is_active( cdname, at_current_timestep_arg = .TRUE. ) ) THEN
         z2d(:,:) = 0._wp
         CALL xios_recv_field( cdname, z2d)
      ENDIF
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, z2d ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_get_var


   FUNCTION iom_getszuld ( kiomid )  
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  iom_getszuld  ***
      !!
      !! ** Purpose : get the size of the unlimited dimension in a file
      !!              (return -1 if not found)
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kiomid   ! file Identifier
      !
      INTEGER                ::   iom_getszuld
      !!-----------------------------------------------------------------------
      iom_getszuld = -1
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%iduld > 0 )   iom_getszuld = iom_file(kiomid)%lenuld
      ENDIF
   END FUNCTION iom_getszuld
   

   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_chkatt
   !!----------------------------------------------------------------------
   SUBROUTINE iom_chkatt( kiomid, cdatt, llok, ksize, cdvar )
      INTEGER         , INTENT(in   )                 ::   kiomid    ! Identifier of the file
      CHARACTER(len=*), INTENT(in   )                 ::   cdatt     ! Name of the attribute
      LOGICAL         , INTENT(  out)                 ::   llok      ! Error code
      INTEGER         , INTENT(  out), OPTIONAL       ::   ksize     ! Size of the attribute array
      CHARACTER(len=*), INTENT(in   ), OPTIONAL       ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_chkatt( kiomid, cdatt, llok, ksize=ksize, cdvar=cdvar )
      ENDIF
      !
   END SUBROUTINE iom_chkatt

   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_getatt
   !!----------------------------------------------------------------------
   SUBROUTINE iom_g0d_iatt( kiomid, cdatt, katt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      INTEGER               , INTENT(  out)           ::   katt0d    ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt,  katt0d =  katt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g0d_iatt

   SUBROUTINE iom_g1d_iatt( kiomid, cdatt, katt1d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      INTEGER, DIMENSION(:) , INTENT(  out)           ::   katt1d    ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt,  katt1d =  katt1d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g1d_iatt

   SUBROUTINE iom_g0d_ratt( kiomid, cdatt, patt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      REAL(wp)              , INTENT(  out)           ::   patt0d    ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt,  patt0d =  patt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g0d_ratt

   SUBROUTINE iom_g1d_ratt( kiomid, cdatt, patt1d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      REAL(wp), DIMENSION(:), INTENT(  out)           ::   patt1d    ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt,  patt1d =  patt1d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g1d_ratt
   
   SUBROUTINE iom_g0d_catt( kiomid, cdatt, cdatt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      CHARACTER(len=*)      , INTENT(  out)           ::   cdatt0d   ! read field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_getatt( kiomid, cdatt, cdatt0d = cdatt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_g0d_catt


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_putatt
   !!----------------------------------------------------------------------
   SUBROUTINE iom_p0d_iatt( kiomid, cdatt, katt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      INTEGER               , INTENT(in   )           ::   katt0d    ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt,  katt0d =  katt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p0d_iatt

   SUBROUTINE iom_p1d_iatt( kiomid, cdatt, katt1d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      INTEGER, DIMENSION(:) , INTENT(in   )           ::   katt1d    ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt,  katt1d =  katt1d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p1d_iatt

   SUBROUTINE iom_p0d_ratt( kiomid, cdatt, patt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      REAL(wp)              , INTENT(in   )           ::   patt0d    ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt,  patt0d =  patt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p0d_ratt

   SUBROUTINE iom_p1d_ratt( kiomid, cdatt, patt1d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      REAL(wp), DIMENSION(:), INTENT(in   )           ::   patt1d    ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt,  patt1d =  patt1d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p1d_ratt
   
   SUBROUTINE iom_p0d_catt( kiomid, cdatt, cdatt0d, cdvar )
      INTEGER               , INTENT(in   )           ::   kiomid    ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt     ! Name of the attribute
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt0d   ! written field
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar     ! Name of the variable
      !
      IF( kiomid > 0 ) THEN
         IF( iom_file(kiomid)%nfid > 0 )   CALL iom_nf90_putatt( kiomid, cdatt, cdatt0d = cdatt0d, cdvar=cdvar )
      ENDIF
   END SUBROUTINE iom_p0d_catt


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_rstput
   !!----------------------------------------------------------------------
   SUBROUTINE iom_rp0d( kt, kwrite, kiomid, cdvar, pvar, ktype, ldxios )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file 
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(wp)        , INTENT(in)                         ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      LOGICAL, OPTIONAL :: ldxios   ! xios write flag
      LOGICAL :: llx                ! local xios write flag
      INTEGER :: ivid   ! variable id

      llx = .FALSE.
      IF(PRESENT(ldxios)) llx = ldxios
      IF( llx ) THEN
#ifdef key_iomput
      IF( kt == kwrite ) THEN
          IF(lwp) write(numout,*) 'RESTART: write (XIOS 0D) ',trim(cdvar)
          CALL xios_send_field(trim(cdvar), pvar)
      ENDIF
#endif
      ELSE
         IF( kiomid > 0 ) THEN
            IF( iom_file(kiomid)%nfid > 0 ) THEN
               ivid = iom_varid( kiomid, cdvar, ldstop = .FALSE. )
               CALL iom_nf90_rstput( kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r0d = pvar )
            ENDIF
         ENDIF
      ENDIF
   END SUBROUTINE iom_rp0d

   SUBROUTINE iom_rp1d( kt, kwrite, kiomid, cdvar, pvar, ktype, ldxios )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file 
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(wp)        , INTENT(in), DIMENSION(          :) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      LOGICAL, OPTIONAL                                    ::   ldxios   ! xios write flag
      LOGICAL :: llx                ! local xios write flag
      INTEGER :: ivid   ! variable id

      llx = .FALSE.
      IF(PRESENT(ldxios)) llx = ldxios
      IF( llx ) THEN
#ifdef key_iomput
      IF( kt == kwrite ) THEN
         IF(lwp) write(numout,*) 'RESTART: write (XIOS 1D) ',trim(cdvar)
         CALL xios_send_field(trim(cdvar), pvar)
      ENDIF
#endif
      ELSE
         IF( kiomid > 0 ) THEN
            IF( iom_file(kiomid)%nfid > 0 ) THEN
               ivid = iom_varid( kiomid, cdvar, ldstop = .FALSE. )
               CALL iom_nf90_rstput( kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r1d = pvar )
            ENDIF
         ENDIF
      ENDIF
   END SUBROUTINE iom_rp1d

   SUBROUTINE iom_rp2d( kt, kwrite, kiomid, cdvar, pvar, ktype, ldxios )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file 
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(wp)        , INTENT(in), DIMENSION(:,    :    ) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      LOGICAL, OPTIONAL :: ldxios   ! xios write flag
      LOGICAL :: llx
      INTEGER :: ivid   ! variable id

      llx = .FALSE.
      IF(PRESENT(ldxios)) llx = ldxios
      IF( llx ) THEN
#ifdef key_iomput
      IF( kt == kwrite ) THEN
         IF(lwp) write(numout,*) 'RESTART: write (XIOS 2D) ',trim(cdvar)
         CALL xios_send_field(trim(cdvar), pvar)
      ENDIF
#endif
      ELSE
         IF( kiomid > 0 ) THEN
            IF( iom_file(kiomid)%nfid > 0 ) THEN
               ivid = iom_varid( kiomid, cdvar, ldstop = .FALSE. )
               CALL iom_nf90_rstput( kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r2d = pvar )
            ENDIF
         ENDIF
      ENDIF
   END SUBROUTINE iom_rp2d

   SUBROUTINE iom_rp3d( kt, kwrite, kiomid, cdvar, pvar, ktype, ldxios )
      INTEGER         , INTENT(in)                         ::   kt       ! ocean time-step
      INTEGER         , INTENT(in)                         ::   kwrite   ! writing time-step
      INTEGER         , INTENT(in)                         ::   kiomid   ! Identifier of the file 
      CHARACTER(len=*), INTENT(in)                         ::   cdvar    ! time axis name
      REAL(wp)        , INTENT(in),       DIMENSION(:,:,:) ::   pvar     ! written field
      INTEGER         , INTENT(in), OPTIONAL               ::   ktype    ! variable external type
      LOGICAL, OPTIONAL :: ldxios   ! xios write flag
      LOGICAL :: llx                 ! local xios write flag
      INTEGER :: ivid   ! variable id

      llx = .FALSE.
      IF(PRESENT(ldxios)) llx = ldxios
      IF( llx ) THEN
#ifdef key_iomput
      IF( kt == kwrite ) THEN
         IF(lwp) write(numout,*) 'RESTART: write (XIOS 3D) ',trim(cdvar)
         CALL xios_send_field(trim(cdvar), pvar)
      ENDIF
#endif
      ELSE
         IF( kiomid > 0 ) THEN
            IF( iom_file(kiomid)%nfid > 0 ) THEN
               ivid = iom_varid( kiomid, cdvar, ldstop = .FALSE. )
               CALL iom_nf90_rstput( kt, kwrite, kiomid, cdvar, ivid, ktype, pv_r3d = pvar )
            ENDIF
         ENDIF
      ENDIF
   END SUBROUTINE iom_rp3d


  SUBROUTINE iom_delay_rst( cdaction, cdcpnt, kncid )
      !!---------------------------------------------------------------------
      !!   Routine iom_delay_rst: used read/write restart related to mpp_delay
      !!
      !!---------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in   ) ::   cdaction        !
      CHARACTER(len=*), INTENT(in   ) ::   cdcpnt
      INTEGER         , INTENT(in   ) ::   kncid
      !
      INTEGER  :: ji
      INTEGER  :: indim
      LOGICAL  :: llattexist
      REAL(wp), ALLOCATABLE, DIMENSION(:) ::   zreal1d
      !!---------------------------------------------------------------------
      !
      !                                      ===================================
      IF( TRIM(cdaction) == 'READ' ) THEN   ! read restart related to mpp_delay !
         !                                   ===================================
         DO ji = 1, nbdelay
            IF ( c_delaycpnt(ji) == cdcpnt ) THEN
               CALL iom_chkatt( kncid, 'DELAY_'//c_delaylist(ji), llattexist, indim )
               IF( llattexist )  THEN
                  ALLOCATE( todelay(ji)%z1d(indim) )
                  CALL iom_getatt( kncid, 'DELAY_'//c_delaylist(ji), todelay(ji)%z1d(:) )
                  ndelayid(ji) = 0   ! set to 0 to specify that the value was read in the restart
               ENDIF
           ENDIF
         END DO
         !                                   ====================================
      ELSE                                  ! write restart related to mpp_delay !
         !                                   ====================================
         DO ji = 1, nbdelay   ! save only ocean delayed global communication variables
            IF ( c_delaycpnt(ji) == cdcpnt ) THEN
               IF( ASSOCIATED(todelay(ji)%z1d) ) THEN
                  CALL mpp_delay_rcv(ji)   ! make sure %z1d is received
                  CALL iom_putatt( kncid, 'DELAY_'//c_delaylist(ji), todelay(ji)%z1d(:) )
               ENDIF
            ENDIF
         END DO
         !
      ENDIF
      
   END SUBROUTINE iom_delay_rst
  
   

   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_put
   !!----------------------------------------------------------------------
   SUBROUTINE iom_p0d( cdname, pfield0d )
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      REAL(wp)        , INTENT(in) ::   pfield0d
!!      REAL(wp)        , DIMENSION(jpi,jpj) ::   zz     ! masson
#if defined key_iomput
!!clem      zz(:,:)=pfield0d
!!clem      CALL xios_send_field(cdname, zz)
      CALL xios_send_field(cdname, (/pfield0d/)) 
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield0d   ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_p0d

   SUBROUTINE iom_p1d( cdname, pfield1d )
      CHARACTER(LEN=*)          , INTENT(in) ::   cdname
      REAL(wp),     DIMENSION(:), INTENT(in) ::   pfield1d
#if defined key_iomput
      CALL xios_send_field( cdname, RESHAPE( (/pfield1d/), (/1,1,SIZE(pfield1d)/) ) )
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield1d   ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_p1d

   SUBROUTINE iom_p2d( cdname, pfield2d )
      CHARACTER(LEN=*)            , INTENT(in) ::   cdname
      REAL(wp),     DIMENSION(:,:), INTENT(in) ::   pfield2d
#if defined key_iomput
      CALL xios_send_field(cdname, pfield2d)
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield2d   ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_p2d

   SUBROUTINE iom_p3d( cdname, pfield3d )
      CHARACTER(LEN=*)                , INTENT(in) ::   cdname
      REAL(wp),       DIMENSION(:,:,:), INTENT(in) ::   pfield3d
#if defined key_iomput
      CALL xios_send_field( cdname, pfield3d )
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield3d   ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_p3d

   SUBROUTINE iom_p4d( cdname, pfield4d )
      CHARACTER(LEN=*)                  , INTENT(in) ::   cdname
      REAL(wp),       DIMENSION(:,:,:,:), INTENT(in) ::   pfield4d
#if defined key_iomput
      CALL xios_send_field(cdname, pfield4d)
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pfield4d   ! useless test to avoid compilation warnings
#endif
   END SUBROUTINE iom_p4d


#if defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_iomput'                                         XIOS interface
   !!----------------------------------------------------------------------

   SUBROUTINE iom_set_domain_attr( cdid, ni_glo, nj_glo, ibegin, jbegin, ni, nj,                                               &
      &                                    data_dim, data_ibegin, data_ni, data_jbegin, data_nj, lonvalue, latvalue, mask,     &
      &                                    nvertex, bounds_lon, bounds_lat, area )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)                  , INTENT(in) ::   cdid
      INTEGER                 , OPTIONAL, INTENT(in) ::   ni_glo, nj_glo, ibegin, jbegin, ni, nj
      INTEGER                 , OPTIONAL, INTENT(in) ::   data_dim, data_ibegin, data_ni, data_jbegin, data_nj
      INTEGER                 , OPTIONAL, INTENT(in) ::   nvertex
      REAL(wp), DIMENSION(:)  , OPTIONAL, INTENT(in) ::   lonvalue, latvalue
      REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(in) ::   bounds_lon, bounds_lat, area
      LOGICAL , DIMENSION(:)  , OPTIONAL, INTENT(in) ::   mask
      !!----------------------------------------------------------------------
      !
      IF( xios_is_valid_domain     (cdid) ) THEN
         CALL xios_set_domain_attr     ( cdid, ni_glo=ni_glo, nj_glo=nj_glo, ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj,   &
            &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni, data_jbegin=data_jbegin, data_nj=data_nj ,   &
            &    lonvalue_1D=lonvalue, latvalue_1D=latvalue, mask_1D=mask, nvertex=nvertex, bounds_lon_1D=bounds_lon,      &
            &    bounds_lat_1D=bounds_lat, area=area, type='curvilinear')
      ENDIF
      IF( xios_is_valid_domaingroup(cdid) ) THEN
         CALL xios_set_domaingroup_attr( cdid, ni_glo=ni_glo, nj_glo=nj_glo, ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj,   &
            &    data_dim=data_dim, data_ibegin=data_ibegin, data_ni=data_ni, data_jbegin=data_jbegin, data_nj=data_nj ,   &
            &    lonvalue_1D=lonvalue, latvalue_1D=latvalue, mask_1D=mask, nvertex=nvertex, bounds_lon_1D=bounds_lon,      &
            &    bounds_lat_1D=bounds_lat, area=area, type='curvilinear' )
      ENDIF
      !
      CALL xios_solve_inheritance()
      !
   END SUBROUTINE iom_set_domain_attr


   SUBROUTINE iom_set_zoom_domain_attr( cdid, ibegin, jbegin, ni, nj )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) ::   cdid
      INTEGER         , INTENT(in) ::   ibegin, jbegin, ni, nj
      !
      TYPE(xios_gridgroup) :: gridgroup_hdl
      TYPE(xios_grid)      :: grid_hdl
      TYPE(xios_domain)    :: domain_hdl 
      TYPE(xios_axis)      :: axis_hdl 
      CHARACTER(LEN=64)    :: cldomrefid   ! domain_ref name
      CHARACTER(len=1)     :: cl1          ! last character of this name
      !!----------------------------------------------------------------------
      !
      IF( xios_is_valid_zoom_domain(cdid) ) THEN
         ! define the zoom_domain attributs
         CALL xios_set_zoom_domain_attr( cdid, ibegin=ibegin, jbegin=jbegin, ni=ni, nj=nj )
         ! define a new 2D grid with this new domain
         CALL xios_get_handle("grid_definition", gridgroup_hdl )
         CALL xios_add_child(gridgroup_hdl, grid_hdl, TRIM(cdid)//'_2D' )   ! add a new 2D grid to grid_definition
         CALL xios_add_child(grid_hdl, domain_hdl, TRIM(cdid) )             ! add its domain
         ! define a new 3D grid with this new domain
         CALL xios_add_child(gridgroup_hdl, grid_hdl, TRIM(cdid)//'_3D' )   ! add a new 3D grid to grid_definition
         CALL xios_add_child(grid_hdl, domain_hdl, TRIM(cdid) )             ! add its domain
         ! vertical axis
         cl1 = cdid(LEN_TRIM(cdid):)                                        ! last letter of cdid
         cl1 = CHAR(ICHAR(cl1)+32)                                          ! from upper to lower case
         CALL xios_add_child(grid_hdl, axis_hdl, 'depth'//cl1)              ! add its axis
      ENDIF
      !      
   END SUBROUTINE iom_set_zoom_domain_attr


   SUBROUTINE iom_set_axis_attr( cdid, paxis, bounds )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)      , INTENT(in) ::   cdid
      REAL(wp), DIMENSION(:)  , OPTIONAL, INTENT(in) ::   paxis
      REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(in) ::   bounds
      !!----------------------------------------------------------------------
      IF( PRESENT(paxis) ) THEN
         IF( xios_is_valid_axis     (cdid) )   CALL xios_set_axis_attr     ( cdid, n_glo=SIZE(paxis), value=paxis )
         IF( xios_is_valid_axisgroup(cdid) )   CALL xios_set_axisgroup_attr( cdid, n_glo=SIZE(paxis), value=paxis )
      ENDIF
      IF( xios_is_valid_axis     (cdid) )   CALL xios_set_axis_attr     ( cdid, bounds=bounds )
      IF( xios_is_valid_axisgroup(cdid) )   CALL xios_set_axisgroup_attr( cdid, bounds=bounds )
      CALL xios_solve_inheritance()
   END SUBROUTINE iom_set_axis_attr


   SUBROUTINE iom_set_field_attr( cdid, freq_op, freq_offset )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)             , INTENT(in) ::   cdid
      TYPE(xios_duration), OPTIONAL, INTENT(in) ::   freq_op
      TYPE(xios_duration), OPTIONAL, INTENT(in) ::   freq_offset
      !!----------------------------------------------------------------------
      IF( xios_is_valid_field     (cdid) )   CALL xios_set_field_attr     ( cdid, freq_op=freq_op, freq_offset=freq_offset )
      IF( xios_is_valid_fieldgroup(cdid) )   CALL xios_set_fieldgroup_attr( cdid, freq_op=freq_op, freq_offset=freq_offset )
      CALL xios_solve_inheritance()
   END SUBROUTINE iom_set_field_attr


   SUBROUTINE iom_set_file_attr( cdid, name, name_suffix )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      CHARACTER(LEN=*),OPTIONAL , INTENT(in) ::   name, name_suffix
      !!----------------------------------------------------------------------
      IF( xios_is_valid_file     (cdid) )   CALL xios_set_file_attr     ( cdid, name=name, name_suffix=name_suffix )
      IF( xios_is_valid_filegroup(cdid) )   CALL xios_set_filegroup_attr( cdid, name=name, name_suffix=name_suffix )
      CALL xios_solve_inheritance()
   END SUBROUTINE iom_set_file_attr


   SUBROUTINE iom_get_file_attr( cdid, name, name_suffix, output_freq )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in ) ::   cdid
      CHARACTER(LEN=*),OPTIONAL , INTENT(out) ::   name, name_suffix
      TYPE(xios_duration), OPTIONAL , INTENT(out) :: output_freq
      LOGICAL                                 ::   llexist1,llexist2,llexist3
      !---------------------------------------------------------------------
      IF( PRESENT( name        ) )   name = ''          ! default values
      IF( PRESENT( name_suffix ) )   name_suffix = ''
      IF( PRESENT( output_freq ) )   output_freq = xios_duration(0,0,0,0,0,0)
      IF( xios_is_valid_file     (cdid) ) THEN
         CALL xios_solve_inheritance()
         CALL xios_is_defined_file_attr     ( cdid, name = llexist1, name_suffix = llexist2, output_freq = llexist3)
         IF(llexist1)   CALL xios_get_file_attr     ( cdid, name = name )
         IF(llexist2)   CALL xios_get_file_attr     ( cdid, name_suffix = name_suffix )
         IF(llexist3)   CALL xios_get_file_attr     ( cdid, output_freq = output_freq )
      ENDIF
      IF( xios_is_valid_filegroup(cdid) ) THEN
         CALL xios_solve_inheritance()
         CALL xios_is_defined_filegroup_attr( cdid, name = llexist1, name_suffix = llexist2, output_freq = llexist3)
         IF(llexist1)   CALL xios_get_filegroup_attr( cdid, name = name )
         IF(llexist2)   CALL xios_get_filegroup_attr( cdid, name_suffix = name_suffix )
         IF(llexist3)   CALL xios_get_filegroup_attr( cdid, output_freq = output_freq )
      ENDIF
   END SUBROUTINE iom_get_file_attr


   SUBROUTINE iom_set_grid_attr( cdid, mask )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)                   , INTENT(in) ::   cdid
      LOGICAL, DIMENSION(:,:,:), OPTIONAL, INTENT(in) ::   mask
      !!----------------------------------------------------------------------
      IF( xios_is_valid_grid     (cdid) )   CALL xios_set_grid_attr     ( cdid, mask_3D=mask )
      IF( xios_is_valid_gridgroup(cdid) )   CALL xios_set_gridgroup_attr( cdid, mask_3D=mask )
      CALL xios_solve_inheritance()
   END SUBROUTINE iom_set_grid_attr

   SUBROUTINE iom_setkt( kt, cdname )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt 
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      !!----------------------------------------------------------------------
      CALL iom_swap( cdname )   ! swap to cdname context
      CALL xios_update_calendar(kt)
      IF( cdname /= TRIM(cxios_context) )   CALL iom_swap( TRIM(cxios_context) )   ! return back to nemo context
   END SUBROUTINE iom_setkt

   SUBROUTINE iom_context_finalize( cdname )
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: cdname
      CHARACTER(LEN=120)           :: clname
      !!----------------------------------------------------------------------
      clname = cdname
      IF( TRIM(Agrif_CFixed()) .NE. '0' ) clname = TRIM(Agrif_CFixed())//"_"//clname 
      IF( xios_is_valid_context(clname) ) THEN
         CALL iom_swap( cdname )   ! swap to cdname context
         CALL xios_context_finalize() ! finalize the context
         IF( cdname /= TRIM(cxios_context) ) CALL iom_swap( TRIM(cxios_context) )   ! return back to nemo context
      ENDIF
      !
   END SUBROUTINE iom_context_finalize


   SUBROUTINE set_grid( cdgrd, plon, plat, ldxios, ldrxios )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_grid  ***
      !!
      !! ** Purpose :   define horizontal grids
      !!----------------------------------------------------------------------
      CHARACTER(LEN=1)            , INTENT(in) ::   cdgrd
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   plon
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   plat
      !
      INTEGER  :: ni, nj
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zmask
      LOGICAL, INTENT(IN) :: ldxios, ldrxios
      !!----------------------------------------------------------------------
      !
      ni = nlei-nldi+1
      nj = nlej-nldj+1
      !
      CALL iom_set_domain_attr("grid_"//cdgrd, ni_glo=jpiglo, nj_glo=jpjglo, ibegin=nimpp+nldi-2, jbegin=njmpp+nldj-2, ni=ni, nj=nj)
      CALL iom_set_domain_attr("grid_"//cdgrd, data_dim=2, data_ibegin = 1-nldi, data_ni = jpi, data_jbegin = 1-nldj, data_nj = jpj)
!don't define lon and lat for restart reading context. 
      IF ( .NOT.ldrxios ) &
         CALL iom_set_domain_attr("grid_"//cdgrd, lonvalue = RESHAPE(plon(nldi:nlei, nldj:nlej),(/ ni*nj /)),   &
         &                                     latvalue = RESHAPE(plat(nldi:nlei, nldj:nlej),(/ ni*nj /)))  
      !
      IF ( ln_mskland .AND. (.NOT.ldxios) ) THEN
         ! mask land points, keep values on coast line -> specific mask for U, V and W points
         SELECT CASE ( cdgrd )
         CASE('T')   ;   zmask(:,:,:)       = tmask(:,:,:)
         CASE('U')   ;   zmask(2:jpim1,:,:) = tmask(2:jpim1,:,:) + tmask(3:jpi,:,:)   ;   CALL lbc_lnk( 'iom', zmask, 'U', 1. )
         CASE('V')   ;   zmask(:,2:jpjm1,:) = tmask(:,2:jpjm1,:) + tmask(:,3:jpj,:)   ;   CALL lbc_lnk( 'iom', zmask, 'V', 1. )
         CASE('W')   ;   zmask(:,:,2:jpk  ) = tmask(:,:,1:jpkm1) + tmask(:,:,2:jpk)   ;   zmask(:,:,1) = tmask(:,:,1)
         END SELECT
         !
         CALL iom_set_domain_attr( "grid_"//cdgrd       , mask = RESHAPE(zmask(nldi:nlei,nldj:nlej,1),(/ni*nj    /)) /= 0. )
         CALL iom_set_grid_attr  ( "grid_"//cdgrd//"_3D", mask = RESHAPE(zmask(nldi:nlei,nldj:nlej,:),(/ni,nj,jpk/)) /= 0. )
      ENDIF
      !
   END SUBROUTINE set_grid


   SUBROUTINE set_grid_bounds( cdgrd, plon_cnr, plat_cnr, plon_pnt, plat_pnt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE set_grid_bounds  ***
      !!
      !! ** Purpose :   define horizontal grid corners
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=1)                      , INTENT(in) :: cdgrd
      REAL(wp), DIMENSION(jpi,jpj)          , INTENT(in) :: plon_cnr, plat_cnr  ! Lat/lon coord. of a contiguous vertex of cell (i,j)
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL, INTENT(in) :: plon_pnt, plat_pnt  ! Lat/lon coord. of the point of cell (i,j)
      !
      INTEGER :: ji, jj, jn, ni, nj
      INTEGER :: icnr, jcnr                                    ! Offset such that the vertex coordinate (i+icnr,j+jcnr)
      !                                                        ! represents the bottom-left corner of cell (i,j)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: z_bnds      ! Lat/lon coordinates of the vertices of cell (i,j)
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)     :: z_fld       ! Working array to determine where to rotate cells
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)     :: z_rot       ! Lat/lon working array for rotation of cells
      !!----------------------------------------------------------------------
      !
      ALLOCATE( z_bnds(4,jpi,jpj,2), z_fld(jpi,jpj), z_rot(4,2)  )
      !
      ! Offset of coordinate representing bottom-left corner
      SELECT CASE ( TRIM(cdgrd) )
      CASE ('T', 'W')   ;   icnr = -1   ;   jcnr = -1
      CASE ('U')        ;   icnr =  0   ;   jcnr = -1
      CASE ('V')        ;   icnr = -1   ;   jcnr =  0
      END SELECT
      !
      ni = nlei-nldi+1   ! Dimensions of subdomain interior
      nj = nlej-nldj+1
      !
      z_fld(:,:) = 1._wp
      CALL lbc_lnk( 'iom', z_fld, cdgrd, -1. )    ! Working array for location of northfold
      !
      ! Cell vertices that can be defined
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            z_bnds(1,ji,jj,1) = plat_cnr(ji+icnr,  jj+jcnr  ) ! Bottom-left
            z_bnds(2,ji,jj,1) = plat_cnr(ji+icnr+1,jj+jcnr  ) ! Bottom-right
            z_bnds(3,ji,jj,1) = plat_cnr(ji+icnr+1,jj+jcnr+1) ! Top-right
            z_bnds(4,ji,jj,1) = plat_cnr(ji+icnr,  jj+jcnr+1) ! Top-left
            z_bnds(1,ji,jj,2) = plon_cnr(ji+icnr,  jj+jcnr  ) ! Bottom-left
            z_bnds(2,ji,jj,2) = plon_cnr(ji+icnr+1,jj+jcnr  ) ! Bottom-right
            z_bnds(3,ji,jj,2) = plon_cnr(ji+icnr+1,jj+jcnr+1) ! Top-right
            z_bnds(4,ji,jj,2) = plon_cnr(ji+icnr,  jj+jcnr+1) ! Top-left
         END DO
      END DO
      !
      ! Cell vertices on boundries
      DO jn = 1, 4
         CALL lbc_lnk( 'iom', z_bnds(jn,:,:,1), cdgrd, 1., pfillval=999._wp )
         CALL lbc_lnk( 'iom', z_bnds(jn,:,:,2), cdgrd, 1., pfillval=999._wp )
      END DO
      !
      ! Zero-size cells at closed boundaries if cell points provided,
      ! otherwise they are closed cells with unrealistic bounds
      IF( PRESENT(plon_pnt) .AND. PRESENT(plat_pnt) ) THEN
         IF( (nbondi == -1 .OR. nbondi == 2) .AND. .NOT. (jperio == 1 .OR. jperio == 4 .OR. jperio == 6) ) THEN
            DO jn = 1, 4        ! (West or jpni = 1), closed E-W
               z_bnds(jn,1,:,1) = plat_pnt(1,:)  ;  z_bnds(jn,1,:,2) = plon_pnt(1,:)
            END DO
         ENDIF
         IF( (nbondi == 1 .OR. nbondi == 2) .AND. .NOT. (jperio == 1 .OR. jperio == 4 .OR. jperio == 6) ) THEN
            DO jn = 1, 4        ! (East or jpni = 1), closed E-W
               z_bnds(jn,nlci,:,1) = plat_pnt(nlci,:)  ;  z_bnds(jn,nlci,:,2) = plon_pnt(nlci,:)
            END DO
         ENDIF
         IF( nbondj == -1 .OR. (nbondj == 2 .AND. jperio /= 2) ) THEN
            DO jn = 1, 4        ! South or (jpnj = 1, not symmetric)
               z_bnds(jn,:,1,1) = plat_pnt(:,1)  ;  z_bnds(jn,:,1,2) = plon_pnt(:,1)
            END DO
         ENDIF
         IF( (nbondj == 1 .OR. nbondj == 2) .AND. jperio  < 3 ) THEN
            DO jn = 1, 4        ! (North or jpnj = 1), no north fold
               z_bnds(jn,:,nlcj,1) = plat_pnt(:,nlcj)  ;  z_bnds(jn,:,nlcj,2) = plon_pnt(:,nlcj)
            END DO
         ENDIF
      ENDIF
      !
      IF( (nbondj == 1 .OR. nbondj == 2) .AND. jperio >= 3 ) THEN    ! Rotate cells at the north fold
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( z_fld(ji,jj) == -1. ) THEN
                  z_rot(1,:) = z_bnds(3,ji,jj,:) ; z_rot(2,:) = z_bnds(4,ji,jj,:)
                  z_rot(3,:) = z_bnds(1,ji,jj,:) ; z_rot(4,:) = z_bnds(2,ji,jj,:)
                  z_bnds(:,ji,jj,:) = z_rot(:,:)
               ENDIF
            END DO
         END DO
      ELSE IF( nbondj == 2 .AND. jperio == 2 ) THEN                  ! Invert cells at the symmetric equator
         DO ji = 1, jpi
            z_rot(1:2,:) = z_bnds(3:4,ji,1,:)
            z_rot(3:4,:) = z_bnds(1:2,ji,1,:)
            z_bnds(:,ji,1,:) = z_rot(:,:)
         END DO
      ENDIF
      !
      CALL iom_set_domain_attr("grid_"//cdgrd, bounds_lat = RESHAPE(z_bnds(:,nldi:nlei,nldj:nlej,1),(/ 4,ni*nj /)),           &
          &                                    bounds_lon = RESHAPE(z_bnds(:,nldi:nlei,nldj:nlej,2),(/ 4,ni*nj /)), nvertex=4 )
      !
      DEALLOCATE( z_bnds, z_fld, z_rot ) 
      !
   END SUBROUTINE set_grid_bounds


   SUBROUTINE set_grid_znl( plat )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_grid_znl  ***
      !!
      !! ** Purpose :   define grids for zonal mean
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   plat
      !
      INTEGER  :: ni, nj, ix, iy
      REAL(wp), DIMENSION(:), ALLOCATABLE  ::   zlon
      !!----------------------------------------------------------------------
      !
      ni=nlei-nldi+1       ! define zonal mean domain (jpj*jpk)
      nj=nlej-nldj+1
      ALLOCATE( zlon(ni*nj) )       ;       zlon(:) = 0._wp
      !
!      CALL dom_ngb( -168.53, 65.03, ix, iy, 'T' ) !  i-line that passes through Bering Strait: Reference latitude (used in plots)
      CALL dom_ngb( 180., 90., ix, iy, 'T' ) !  i-line that passes near the North Pole : Reference latitude (used in plots)
      CALL iom_set_domain_attr("gznl", ni_glo=jpiglo, nj_glo=jpjglo, ibegin=nimpp+nldi-2, jbegin=njmpp+nldj-2, ni=ni, nj=nj)
      CALL iom_set_domain_attr("gznl", data_dim=2, data_ibegin = 1-nldi, data_ni = jpi, data_jbegin = 1-nldj, data_nj = jpj)
      CALL iom_set_domain_attr("gznl", lonvalue = zlon,   &
         &                             latvalue = RESHAPE(plat(nldi:nlei, nldj:nlej),(/ ni*nj /)))  
      CALL iom_set_zoom_domain_attr("ptr", ibegin=ix-1, jbegin=0, ni=1, nj=jpjglo)
      !
      CALL iom_update_file_name('ptr')
      !
   END SUBROUTINE set_grid_znl


   SUBROUTINE set_scalar
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_scalar  ***
      !!
      !! ** Purpose :   define fake grids for scalar point
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(1)   ::   zz = 1.
      !!----------------------------------------------------------------------
      !
      CALL iom_set_domain_attr('scalarpoint', ni_glo=jpnij, nj_glo=1, ibegin=narea-1, jbegin=0, ni=1, nj=1)
      CALL iom_set_domain_attr('scalarpoint', data_dim=2, data_ibegin = 1, data_ni = 1, data_jbegin = 1, data_nj = 1)
      !
      zz = REAL( narea, wp )
      CALL iom_set_domain_attr('scalarpoint', lonvalue=zz, latvalue=zz)
      !
   END SUBROUTINE set_scalar


   SUBROUTINE set_xmlatt
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_xmlatt  ***
      !!
      !! ** Purpose :   automatic definitions of some of the xml attributs...
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=1),DIMENSION( 3) ::   clgrd                    ! suffix name
      CHARACTER(len=256)             ::   clsuff                   ! suffix name
      CHARACTER(len=1)               ::   cl1                      ! 1 character
      CHARACTER(len=2)               ::   cl2                      ! 2 characters
      CHARACTER(len=3)               ::   cl3                      ! 3 characters
      INTEGER                        ::   ji, jg                   ! loop counters
      INTEGER                        ::   ix, iy                   ! i-,j- index
      REAL(wp)        ,DIMENSION(11) ::   zlontao                  ! longitudes of tao    moorings
      REAL(wp)        ,DIMENSION( 7) ::   zlattao                  ! latitudes  of tao    moorings
      REAL(wp)        ,DIMENSION( 4) ::   zlonrama                 ! longitudes of rama   moorings
      REAL(wp)        ,DIMENSION(11) ::   zlatrama                 ! latitudes  of rama   moorings
      REAL(wp)        ,DIMENSION( 3) ::   zlonpira                 ! longitudes of pirata moorings
      REAL(wp)        ,DIMENSION( 9) ::   zlatpira                 ! latitudes  of pirata moorings
      TYPE(xios_duration)            ::   f_op, f_of
      !!----------------------------------------------------------------------
      ! 
      ! frequency of the call of iom_put (attribut: freq_op)
      f_op%timestep = 1        ;  f_of%timestep =  0  ; CALL iom_set_field_attr('field_definition', freq_op=f_op, freq_offset=f_of)
      f_op%timestep = 2        ;  f_of%timestep =  0  ; CALL iom_set_field_attr('trendT_even'     , freq_op=f_op, freq_offset=f_of)
      f_op%timestep = 2        ;  f_of%timestep = -1  ; CALL iom_set_field_attr('trendT_odd'      , freq_op=f_op, freq_offset=f_of)
      f_op%timestep = nn_fsbc  ;  f_of%timestep =  0  ; CALL iom_set_field_attr('SBC'             , freq_op=f_op, freq_offset=f_of)
      f_op%timestep = nn_fsbc  ;  f_of%timestep =  0  ; CALL iom_set_field_attr('SBC_scalar'      , freq_op=f_op, freq_offset=f_of)
      f_op%timestep = nn_dttrc ;  f_of%timestep =  0  ; CALL iom_set_field_attr('ptrc_T'          , freq_op=f_op, freq_offset=f_of)
      f_op%timestep = nn_dttrc ;  f_of%timestep =  0  ; CALL iom_set_field_attr('diad_T'          , freq_op=f_op, freq_offset=f_of)

      ! output file names (attribut: name)
      DO ji = 1, 9
         WRITE(cl1,'(i1)') ji 
         CALL iom_update_file_name('file'//cl1)
      END DO
      DO ji = 1, 99
         WRITE(cl2,'(i2.2)') ji 
         CALL iom_update_file_name('file'//cl2)
      END DO
      DO ji = 1, 999
         WRITE(cl3,'(i3.3)') ji 
         CALL iom_update_file_name('file'//cl3)
      END DO

      ! Zooms...
      clgrd = (/ 'T', 'U', 'W' /) 
      DO jg = 1, SIZE(clgrd)                                                                   ! grid type
         cl1 = clgrd(jg)
         ! Equatorial section (attributs: jbegin, ni, name_suffix)
         CALL dom_ngb( 0., 0., ix, iy, cl1 )
         CALL iom_set_zoom_domain_attr('Eq'//cl1, ibegin=0, jbegin=iy-1, ni=jpiglo, nj=1 )
         CALL iom_get_file_attr   ('Eq'//cl1, name_suffix = clsuff             )
         CALL iom_set_file_attr   ('Eq'//cl1, name_suffix = TRIM(clsuff)//'_Eq')
         CALL iom_update_file_name('Eq'//cl1)
      END DO
      ! TAO moorings (attributs: ibegin, jbegin, name_suffix)
      zlontao = (/ 137.0, 147.0, 156.0, 165.0, -180.0, -170.0, -155.0, -140.0, -125.0, -110.0, -95.0 /)
      zlattao = (/  -8.0,  -5.0,  -2.0,   0.0,    2.0,    5.0,    8.0 /)
      CALL set_mooring( zlontao, zlattao )
      ! RAMA moorings (attributs: ibegin, jbegin, name_suffix)
      zlonrama = (/  55.0,  67.0, 80.5, 90.0 /)
      zlatrama = (/ -16.0, -12.0, -8.0, -4.0, -1.5, 0.0, 1.5, 4.0, 8.0, 12.0, 15.0 /)
      CALL set_mooring( zlonrama, zlatrama )
      ! PIRATA moorings (attributs: ibegin, jbegin, name_suffix)
      zlonpira = (/ -38.0, -23.0, -10.0 /)
      zlatpira = (/ -19.0, -14.0,  -8.0, 0.0, 4.0, 8.0, 12.0, 15.0, 20.0 /)
      CALL set_mooring( zlonpira, zlatpira )
      !
   END SUBROUTINE set_xmlatt


   SUBROUTINE set_mooring( plon, plat )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE set_mooring  ***
      !!
      !! ** Purpose :   automatic definitions of moorings xml attributs...
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(in) ::   plon, plat   ! longitudes/latitudes oft the mooring
      !
!!$      CHARACTER(len=1),DIMENSION(4) ::   clgrd = (/ 'T', 'U', 'V', 'W' /)   ! suffix name
      CHARACTER(len=1),DIMENSION(1) ::   clgrd = (/ 'T' /)        ! suffix name
      CHARACTER(len=256)            ::   clname                   ! file name
      CHARACTER(len=256)            ::   clsuff                   ! suffix name
      CHARACTER(len=1)              ::   cl1                      ! 1 character
      CHARACTER(len=6)              ::   clon,clat                ! name of longitude, latitude
      INTEGER                       ::   ji, jj, jg               ! loop counters
      INTEGER                       ::   ix, iy                   ! i-,j- index
      REAL(wp)                      ::   zlon, zlat
      !!----------------------------------------------------------------------
      DO jg = 1, SIZE(clgrd)
         cl1 = clgrd(jg)
         DO ji = 1, SIZE(plon)
            DO jj = 1, SIZE(plat)
               zlon = plon(ji)
               zlat = plat(jj)
               ! modifications for RAMA moorings
               IF( zlon ==  67. .AND. zlat ==  15. )   zlon =  65.
               IF( zlon ==  90. .AND. zlat <=  -4. )   zlon =  95.
               IF( zlon ==  95. .AND. zlat ==  -4. )   zlat =  -5.
               ! modifications for PIRATA moorings
               IF( zlon == -38. .AND. zlat == -19. )   zlon = -34.
               IF( zlon == -38. .AND. zlat == -14. )   zlon = -32.
               IF( zlon == -38. .AND. zlat ==  -8. )   zlon = -30.
               IF( zlon == -38. .AND. zlat ==   0. )   zlon = -35.
               IF( zlon == -23. .AND. zlat ==  20. )   zlat =  21.
               IF( zlon == -10. .AND. zlat == -14. )   zlat = -10.
               IF( zlon == -10. .AND. zlat ==  -8. )   zlat =  -6.
               IF( zlon == -10. .AND. zlat ==   4. ) THEN   ;   zlon = 0.   ;   zlat = 0.   ;   ENDIF
               CALL dom_ngb( zlon, zlat, ix, iy, cl1 )
               IF( zlon >= 0. ) THEN  
                  IF( zlon == REAL(NINT(zlon), wp) ) THEN   ;   WRITE(clon, '(i3,  a)') NINT( zlon), 'e'
                  ELSE                                      ;   WRITE(clon, '(f5.1,a)')       zlon , 'e'
                  ENDIF
               ELSE             
                  IF( zlon == REAL(NINT(zlon), wp) ) THEN   ;   WRITE(clon, '(i3,  a)') NINT(-zlon), 'w'
                  ELSE                                      ;   WRITE(clon, '(f5.1,a)')      -zlon , 'w'
                  ENDIF
               ENDIF
               IF( zlat >= 0. ) THEN  
                  IF( zlat == REAL(NINT(zlat), wp) ) THEN   ;   WRITE(clat, '(i2,  a)') NINT( zlat), 'n'
                  ELSE                                      ;   WRITE(clat, '(f4.1,a)')       zlat , 'n'
                  ENDIF
               ELSE             
                  IF( zlat == REAL(NINT(zlat), wp) ) THEN   ;   WRITE(clat, '(i2,  a)') NINT(-zlat), 's'
                  ELSE                                      ;   WRITE(clat, '(f4.1,a)')      -zlat , 's'
                  ENDIF
               ENDIF
               clname = TRIM(ADJUSTL(clat))//TRIM(ADJUSTL(clon))
               CALL iom_set_zoom_domain_attr(TRIM(clname)//cl1, ibegin= ix-1, jbegin= iy-1, ni=1, nj=1)

               CALL iom_get_file_attr   (TRIM(clname)//cl1, name_suffix = clsuff                         )
               CALL iom_set_file_attr   (TRIM(clname)//cl1, name_suffix = TRIM(clsuff)//'_'//TRIM(clname))
               CALL iom_update_file_name(TRIM(clname)//cl1)
            END DO
         END DO
      END DO
      
   END SUBROUTINE set_mooring

   
   SUBROUTINE iom_update_file_name( cdid )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_update_file_name  ***
      !!
      !! ** Purpose :   
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*)          , INTENT(in) ::   cdid
      !
      CHARACTER(LEN=256) ::   clname
      CHARACTER(LEN=20)  ::   clfreq
      CHARACTER(LEN=20)  ::   cldate
      CHARACTER(LEN=256) ::   cltmpn                 !FUS needed for correct path with AGRIF
      INTEGER            ::   iln                    !FUS needed for correct path with AGRIF
      INTEGER            ::   idx
      INTEGER            ::   jn
      INTEGER            ::   itrlen
      INTEGER            ::   iyear, imonth, iday, isec
      REAL(wp)           ::   zsec
      LOGICAL            ::   llexist
      TYPE(xios_duration)   ::   output_freq 
      !!----------------------------------------------------------------------
      !
      DO jn = 1, 2
         !
         output_freq = xios_duration(0,0,0,0,0,0)
         IF( jn == 1 )   CALL iom_get_file_attr( cdid, name        = clname, output_freq = output_freq )
         IF( jn == 2 )   CALL iom_get_file_attr( cdid, name_suffix = clname )
         !
         IF ( TRIM(clname) /= '' ) THEN 
            !
            idx = INDEX(clname,'@expname@') + INDEX(clname,'@EXPNAME@')
            DO WHILE ( idx /= 0 ) 
               clname = clname(1:idx-1)//TRIM(cexper)//clname(idx+9:LEN_TRIM(clname))
               idx = INDEX(clname,'@expname@') + INDEX(clname,'@EXPNAME@')
            END DO
            !
            idx = INDEX(clname,'@freq@') + INDEX(clname,'@FREQ@')
            DO WHILE ( idx /= 0 ) 
              IF ( output_freq%timestep /= 0) THEN
                  WRITE(clfreq,'(I18,A2)')INT(output_freq%timestep),'ts' 
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%second /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%second),'s' 
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%minute /= 0 ) THEN
                  WRITE(clfreq,'(I18,A2)')INT(output_freq%minute),'mi' 
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%hour /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%hour),'h' 
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%day /= 0 ) THEN
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%day),'d' 
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%month /= 0 ) THEN   
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%month),'m' 
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE IF ( output_freq%year /= 0 ) THEN   
                  WRITE(clfreq,'(I19,A1)')INT(output_freq%year),'y' 
                  itrlen = LEN_TRIM(ADJUSTL(clfreq))
              ELSE
                  CALL ctl_stop('error in the name of file id '//TRIM(cdid),   &
                     & ' attribute output_freq is undefined -> cannot replace @freq@ in '//TRIM(clname) )
              ENDIF
              clname = clname(1:idx-1)//TRIM(ADJUSTL(clfreq))//clname(idx+6:LEN_TRIM(clname))
              idx = INDEX(clname,'@freq@') + INDEX(clname,'@FREQ@')
            END DO
            !
            idx = INDEX(clname,'@startdate@') + INDEX(clname,'@STARTDATE@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( fjulday - rdt / rday )
               clname = clname(1:idx-1)//TRIM(cldate)//clname(idx+11:LEN_TRIM(clname))
               idx = INDEX(clname,'@startdate@') + INDEX(clname,'@STARTDATE@')
            END DO
            !
            idx = INDEX(clname,'@startdatefull@') + INDEX(clname,'@STARTDATEFULL@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( fjulday - rdt / rday, ldfull = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)//clname(idx+15:LEN_TRIM(clname))
               idx = INDEX(clname,'@startdatefull@') + INDEX(clname,'@STARTDATEFULL@')
            END DO
            !
            idx = INDEX(clname,'@enddate@') + INDEX(clname,'@ENDDATE@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( fjulday + rdt / rday * REAL( nitend - nit000, wp ), ld24 = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)//clname(idx+9:LEN_TRIM(clname))
               idx = INDEX(clname,'@enddate@') + INDEX(clname,'@ENDDATE@')
            END DO
            !
            idx = INDEX(clname,'@enddatefull@') + INDEX(clname,'@ENDDATEFULL@')
            DO WHILE ( idx /= 0 ) 
               cldate = iom_sdate( fjulday + rdt / rday * REAL( nitend - nit000, wp ), ld24 = .TRUE., ldfull = .TRUE. )
               clname = clname(1:idx-1)//TRIM(cldate)//clname(idx+13:LEN_TRIM(clname))
               idx = INDEX(clname,'@enddatefull@') + INDEX(clname,'@ENDDATEFULL@')
            END DO
            !
!FUS            IF( jn == 1 .AND. TRIM(Agrif_CFixed()) /= '0' )   clname = TRIM(Agrif_CFixed())//"_"//TRIM(clname)
!FUS see comment line 700 
            IF( jn == 1 .AND. TRIM(Agrif_CFixed()) /= '0' ) THEN
             iln    = INDEX(clname,'/',BACK=.true.)
             cltmpn = clname(1:iln)
             clname = clname(iln+1:LEN_TRIM(clname))
             clname = TRIM(cltmpn)//TRIM(Agrif_CFixed())//'_'//TRIM(clname)
            ENDIF
!FUS 
            IF( jn == 1 )   CALL iom_set_file_attr( cdid, name        = clname )
            IF( jn == 2 )   CALL iom_set_file_attr( cdid, name_suffix = clname )
            !
         ENDIF
         !
      END DO
      !
   END SUBROUTINE iom_update_file_name


   FUNCTION iom_sdate( pjday, ld24, ldfull )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE iom_sdate  ***
      !!
      !! ** Purpose :   send back the date corresponding to the given julian day
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )           ::   pjday    ! julian day
      LOGICAL , INTENT(in   ), OPTIONAL ::   ld24     ! true to force 24:00 instead of 00:00
      LOGICAL , INTENT(in   ), OPTIONAL ::   ldfull   ! true to get the compleate date: yyyymmdd_hh:mm:ss
      !
      CHARACTER(LEN=20) ::   iom_sdate
      CHARACTER(LEN=50) ::   clfmt                         !  format used to write the date 
      INTEGER           ::   iyear, imonth, iday, ihour, iminute, isec
      REAL(wp)          ::   zsec
      LOGICAL           ::   ll24, llfull
      !!----------------------------------------------------------------------
      !
      IF( PRESENT(ld24) ) THEN   ;   ll24 = ld24
      ELSE                       ;   ll24 = .FALSE.
      ENDIF
      !
      IF( PRESENT(ldfull) ) THEN   ;   llfull = ldfull
      ELSE                         ;   llfull = .FALSE.
      ENDIF
      !
      CALL ju2ymds( pjday, iyear, imonth, iday, zsec )
      isec = NINT(zsec)
      !
      IF ( ll24 .AND. isec == 0 ) THEN   ! 00:00 of the next day -> move to 24:00 of the current day
         CALL ju2ymds( pjday - 1., iyear, imonth, iday, zsec )
         isec = 86400
      ENDIF
      !
      IF( iyear < 10000 ) THEN   ;   clfmt = "i4.4,2i2.2"                ! format used to write the date 
      ELSE                       ;   WRITE(clfmt, "('i',i1,',2i2.2')") INT(LOG10(REAL(iyear,wp))) + 1
      ENDIF
      !
!$AGRIF_DO_NOT_TREAT      
      ! needed in the conv
      IF( llfull ) THEN 
         clfmt = TRIM(clfmt)//",'_',i2.2,':',i2.2,':',i2.2"
         ihour   = isec / 3600
         isec    = MOD(isec, 3600)
         iminute = isec / 60
         isec    = MOD(isec, 60)
         WRITE(iom_sdate, '('//TRIM(clfmt)//')') iyear, imonth, iday, ihour, iminute, isec    ! date of the end of run
      ELSE
         WRITE(iom_sdate, '('//TRIM(clfmt)//')') iyear, imonth, iday                          ! date of the end of run
      ENDIF
!$AGRIF_END_DO_NOT_TREAT      
      !
   END FUNCTION iom_sdate

#else
   !!----------------------------------------------------------------------
   !!   NOT 'key_iomput'                               a few dummy routines
   !!----------------------------------------------------------------------
   SUBROUTINE iom_setkt( kt, cdname )
      INTEGER         , INTENT(in)::   kt 
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      IF( .FALSE. )   WRITE(numout,*) kt, cdname   ! useless test to avoid compilation warnings
   END SUBROUTINE iom_setkt

   SUBROUTINE iom_context_finalize( cdname )
      CHARACTER(LEN=*), INTENT(in) ::   cdname
      IF( .FALSE. )   WRITE(numout,*)  cdname   ! useless test to avoid compilation warnings
   END SUBROUTINE iom_context_finalize

#endif

   LOGICAL FUNCTION iom_use( cdname )
      CHARACTER(LEN=*), INTENT(in) ::   cdname
#if defined key_iomput
      iom_use = xios_field_is_active( cdname )
#else
      iom_use = .FALSE.
#endif
   END FUNCTION iom_use

   SUBROUTINE iom_miss_val( cdname, pmiss_val )
      CHARACTER(LEN=*), INTENT(in ) ::   cdname
      REAL(wp)        , INTENT(out) ::   pmiss_val   
#if defined key_iomput
      ! get missing value
      CALL xios_get_field_attr( cdname, default_value = pmiss_val )
#else
      IF( .FALSE. )   WRITE(numout,*) cdname, pmiss_val   ! useless test to avoid compilation warnings
      IF( .FALSE. )   pmiss_val = 0._wp                   ! useless assignment to avoid compilation warnings
#endif
   END SUBROUTINE iom_miss_val
  
   !!======================================================================
END MODULE iom
