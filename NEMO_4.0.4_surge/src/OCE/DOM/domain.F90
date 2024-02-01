MODULE domain
   !!==============================================================================
   !!                       ***  MODULE domain   ***
   !! Ocean initialization : domain initialization
   !!==============================================================================
   !! History :  OPA  !  1990-10  (C. Levy - G. Madec)  Original code
   !!                 !  1992-01  (M. Imbard) insert time step initialization
   !!                 !  1996-06  (G. Madec) generalized vertical coordinate 
   !!                 !  1997-02  (G. Madec) creation of domwri.F
   !!                 !  2001-05  (E.Durand - G. Madec) insert closed sea
   !!   NEMO     1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!            2.0  !  2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.3  !  2010-11  (G. Madec)  initialisation in C1D configuration
   !!            3.6  !  2013     ( J. Simeon, C. Calone, G. Madec, C. Ethe ) Online coarsening of outputs
   !!            3.7  !  2015-11  (G. Madec, A. Coward)  time varying zgr by default
   !!            4.0  !  2016-10  (G. Madec, S. Flavoni)  domain configuration / user defined interface
   !!----------------------------------------------------------------------
   
   !!----------------------------------------------------------------------
   !!   dom_init      : initialize the space and time domain
   !!   dom_glo       : initialize global domain <--> local domain indices
   !!   dom_nam       : read and contral domain namelists
   !!   dom_ctl       : control print for the ocean domain
   !!   domain_cfg    : read the global domain size in domain configuration file
   !!   cfg_write     : create the domain configuration file
   !!----------------------------------------------------------------------
   USE oce            ! ocean variables
   USE dom_oce        ! domain: ocean
   USE sbc_oce        ! surface boundary condition: ocean
   USE trc_oce        ! shared ocean & passive tracers variab
   USE phycst         ! physical constants
   USE closea         ! closed seas
   USE domhgr         ! domain: set the horizontal mesh
   USE domzgr         ! domain: set the vertical mesh
   USE dommsk         ! domain: set the mask system
   USE domwri         ! domain: write the meshmask file
   USE domvvl         ! variable volume
   USE c1d            ! 1D configuration
   USE dyncor_c1d     ! 1D configuration: Coriolis term    (cor_c1d routine)
   USE wet_dry,  ONLY : ll_wd
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library
   USE lbclnk         ! ocean lateral boundary condition (or mpp link)
   USE lib_mpp        ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_init     ! called by nemogcm.F90
   PUBLIC   domain_cfg   ! called by nemogcm.F90

   !!-------------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!-------------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_init(cdstr)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_init  ***
      !!                    
      !! ** Purpose :   Domain initialization. Call the routines that are 
      !!              required to create the arrays which define the space 
      !!              and time domain of the ocean model.
      !!
      !! ** Method  : - dom_msk: compute the masks from the bathymetry file
      !!              - dom_hgr: compute or read the horizontal grid-point position
      !!                         and scale factors, and the coriolis factor
      !!              - dom_zgr: define the vertical coordinate and the bathymetry
      !!              - dom_wri: create the meshmask file (ln_meshmask=T)
      !!              - 1D configuration, move Coriolis, u and v at T-point
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk, ik   ! dummy loop indices
      INTEGER ::   iconf = 0    ! local integers
      CHARACTER (len=64) ::   cform = "(A12, 3(A13, I7))" 
      CHARACTER (len=*), INTENT(IN) :: cdstr                  ! model: NEMO or SAS. Determines core restart variables
      INTEGER , DIMENSION(jpi,jpj) ::   ik_top , ik_bot       ! top and bottom ocean level
      REAL(wp), DIMENSION(jpi,jpj) ::   z1_hu_0, z1_hv_0
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN         ! Ocean domain Parameters (control print)
         WRITE(numout,*)
         WRITE(numout,*) 'dom_init : domain initialization'
         WRITE(numout,*) '~~~~~~~~'
         !
         WRITE(numout,*)     '   Domain info'
         WRITE(numout,*)     '      dimension of model:'
         WRITE(numout,*)     '             Local domain      Global domain       Data domain '
         WRITE(numout,cform) '        ','   jpi     : ', jpi, '   jpiglo  : ', jpiglo
         WRITE(numout,cform) '        ','   jpj     : ', jpj, '   jpjglo  : ', jpjglo
         WRITE(numout,cform) '        ','   jpk     : ', jpk, '   jpkglo  : ', jpkglo
         WRITE(numout,cform) '       ' ,'   jpij    : ', jpij
         WRITE(numout,*)     '      mpp local domain info (mpp):'
         WRITE(numout,*)     '              jpni    : ', jpni, '   nn_hls  : ', nn_hls
         WRITE(numout,*)     '              jpnj    : ', jpnj, '   nn_hls  : ', nn_hls
         WRITE(numout,*)     '              jpnij   : ', jpnij
         WRITE(numout,*)     '      lateral boundary of the Global domain : jperio  = ', jperio
         SELECT CASE ( jperio )
         CASE( 0 )   ;   WRITE(numout,*) '         (i.e. closed)'
         CASE( 1 )   ;   WRITE(numout,*) '         (i.e. cyclic east-west)'
         CASE( 2 )   ;   WRITE(numout,*) '         (i.e. cyclic north-south)'
         CASE( 3 )   ;   WRITE(numout,*) '         (i.e. north fold with T-point pivot)'
         CASE( 4 )   ;   WRITE(numout,*) '         (i.e. cyclic east-west and north fold with T-point pivot)'
         CASE( 5 )   ;   WRITE(numout,*) '         (i.e. north fold with F-point pivot)'
         CASE( 6 )   ;   WRITE(numout,*) '         (i.e. cyclic east-west and north fold with F-point pivot)'
         CASE( 7 )   ;   WRITE(numout,*) '         (i.e. cyclic east-west and north-south)'
         CASE DEFAULT
            CALL ctl_stop( 'jperio is out of range' )
         END SELECT
         WRITE(numout,*)     '      Ocean model configuration used:'
         WRITE(numout,*)     '         cn_cfg = ', TRIM( cn_cfg ), '   nn_cfg = ', nn_cfg
      ENDIF
      nn_wxios = 0
      ln_xios_read = .FALSE.
      !
      !           !==  Reference coordinate system  ==!
      !
      CALL dom_glo                     ! global domain versus local domain
      CALL dom_nam                     ! read namelist ( namrun, namdom )
      !
      IF( lwxios ) THEN
!define names for restart write and set core output (restart.F90)
         CALL iom_set_rst_vars(rst_wfields)
         CALL iom_set_rstw_core(cdstr)
      ENDIF
!reset namelist for SAS
      IF(cdstr == 'SAS') THEN
         IF(lrxios) THEN
               IF(lwp) write(numout,*) 'Disable reading restart file using XIOS for SAS'
               lrxios = .FALSE.
         ENDIF
      ENDIF
      !
      CALL dom_hgr                     ! Horizontal mesh
      CALL dom_zgr( ik_top, ik_bot )   ! Vertical mesh and bathymetry
      CALL dom_msk( ik_top, ik_bot )   ! Masks
      IF( ln_closea )   CALL dom_clo   ! ln_closea=T : closed seas included in the simulation
                                       ! Read in masks to define closed seas and lakes 
      !
      DO jj = 1, jpj                   ! depth of the iceshelves
         DO ji = 1, jpi
            ik = mikt(ji,jj)
            risfdep(ji,jj) = gdepw_0(ji,jj,ik)
         END DO
      END DO
      !
      ht_0(:,:) = 0._wp  ! Reference ocean thickness
      hu_0(:,:) = 0._wp
      hv_0(:,:) = 0._wp
      DO jk = 1, jpk
         ht_0(:,:) = ht_0(:,:) + e3t_0(:,:,jk) * tmask(:,:,jk)
         hu_0(:,:) = hu_0(:,:) + e3u_0(:,:,jk) * umask(:,:,jk)
         hv_0(:,:) = hv_0(:,:) + e3v_0(:,:,jk) * vmask(:,:,jk)
      END DO
      !
      !           !==  time varying part of coordinate system  ==!
      !
      IF( ln_linssh ) THEN       != Fix in time : set to the reference one for all
      !
         !       before        !          now          !       after         !
            gdept_b = gdept_0  ;   gdept_n = gdept_0   !        ---          ! depth of grid-points
            gdepw_b = gdepw_0  ;   gdepw_n = gdepw_0   !        ---          !
                                   gde3w_n = gde3w_0   !        ---          !
         !                                                                  
              e3t_b =   e3t_0  ;     e3t_n =   e3t_0   ;   e3t_a =  e3t_0    ! scale factors
              e3u_b =   e3u_0  ;     e3u_n =   e3u_0   ;   e3u_a =  e3u_0    !
              e3v_b =   e3v_0  ;     e3v_n =   e3v_0   ;   e3v_a =  e3v_0    !
                                     e3f_n =   e3f_0   !        ---          !
              e3w_b =   e3w_0  ;     e3w_n =   e3w_0   !        ---          !
             e3uw_b =  e3uw_0  ;    e3uw_n =  e3uw_0   !        ---          !
             e3vw_b =  e3vw_0  ;    e3vw_n =  e3vw_0   !        ---          !
         !
         z1_hu_0(:,:) = ssumask(:,:) / ( hu_0(:,:) + 1._wp - ssumask(:,:) )     ! _i mask due to ISF
         z1_hv_0(:,:) = ssvmask(:,:) / ( hv_0(:,:) + 1._wp - ssvmask(:,:) )
         !
         !        before       !          now          !       after         !
                                      ht_n =    ht_0   !                     ! water column thickness
               hu_b =    hu_0  ;      hu_n =    hu_0   ;    hu_a =    hu_0   ! 
               hv_b =    hv_0  ;      hv_n =    hv_0   ;    hv_a =    hv_0   !
            r1_hu_b = z1_hu_0  ;   r1_hu_n = z1_hu_0   ; r1_hu_a = z1_hu_0   ! inverse of water column thickness
            r1_hv_b = z1_hv_0  ;   r1_hv_n = z1_hv_0   ; r1_hv_a = z1_hv_0   !
         !
         !
      ELSE                       != time varying : initialize before/now/after variables
         !
         IF( .NOT.l_offline )  CALL dom_vvl_init 
         !
      ENDIF
      !
      IF( lk_c1d         )   CALL cor_c1d       ! 1D configuration: Coriolis set at T-point
      !
      IF( ln_meshmask .AND. .NOT.ln_iscpl )                        CALL dom_wri     ! Create a domain file
      IF( ln_meshmask .AND.      ln_iscpl .AND. .NOT.ln_rstart )   CALL dom_wri     ! Create a domain file
      IF(                                       .NOT.ln_rstart )   CALL dom_ctl     ! Domain control
      !
      IF( ln_write_cfg )   CALL cfg_write         ! create the configuration file
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_init :   ==>>>   END of domain initialization'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) 
      ENDIF
      !
   END SUBROUTINE dom_init


   SUBROUTINE dom_glo
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_glo  ***
      !!
      !! ** Purpose :   initialization of global domain <--> local domain indices
      !!
      !! ** Method  :   
      !!
      !! ** Action  : - mig , mjg : local  domain indices ==> global domain indices
      !!              - mi0 , mi1 : global domain indices ==> local  domain indices
      !!              - mj0,, mj1   (global point not in the local domain ==> mi0>mi1 and/or mj0>mj1)
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj   ! dummy loop argument
      !!----------------------------------------------------------------------
      !
      DO ji = 1, jpi                 ! local domain indices ==> global domain indices
        mig(ji) = ji + nimpp - 1
      END DO
      DO jj = 1, jpj
        mjg(jj) = jj + njmpp - 1
      END DO
      !                              ! global domain indices ==> local domain indices
      !                                   ! (return (m.0,m.1)=(1,0) if data domain gridpoint is to the west/south of the 
      !                                   ! local domain, or (m.0,m.1)=(jp.+1,jp.) to the east/north of local domain. 
      DO ji = 1, jpiglo
        mi0(ji) = MAX( 1 , MIN( ji - nimpp + 1, jpi+1 ) )
        mi1(ji) = MAX( 0 , MIN( ji - nimpp + 1, jpi   ) )
      END DO
      DO jj = 1, jpjglo
        mj0(jj) = MAX( 1 , MIN( jj - njmpp + 1, jpj+1 ) )
        mj1(jj) = MAX( 0 , MIN( jj - njmpp + 1, jpj   ) )
      END DO
      IF(lwp) THEN                   ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dom_glo : domain: global <<==>> local '
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   global domain:   jpiglo = ', jpiglo, ' jpjglo = ', jpjglo, ' jpkglo = ', jpkglo
         WRITE(numout,*) '   local  domain:   jpi    = ', jpi   , ' jpj    = ', jpj   , ' jpk    = ', jpk
         WRITE(numout,*)
         WRITE(numout,*) '   conversion from local to global domain indices (and vise versa) done'
         IF( nn_print >= 1 ) THEN
            WRITE(numout,*)
            WRITE(numout,*) '          conversion local  ==> global i-index domain (mig)'
            WRITE(numout,25)              (mig(ji),ji = 1,jpi)
            WRITE(numout,*)
            WRITE(numout,*) '          conversion global ==> local  i-index domain'
            WRITE(numout,*) '             starting index (mi0)'
            WRITE(numout,25)              (mi0(ji),ji = 1,jpiglo)
            WRITE(numout,*) '             ending index (mi1)'
            WRITE(numout,25)              (mi1(ji),ji = 1,jpiglo)
            WRITE(numout,*)
            WRITE(numout,*) '          conversion local  ==> global j-index domain (mjg)'
            WRITE(numout,25)              (mjg(jj),jj = 1,jpj)
            WRITE(numout,*)
            WRITE(numout,*) '          conversion global ==> local  j-index domain'
            WRITE(numout,*) '             starting index (mj0)'
            WRITE(numout,25)              (mj0(jj),jj = 1,jpjglo)
            WRITE(numout,*) '             ending index (mj1)'
            WRITE(numout,25)              (mj1(jj),jj = 1,jpjglo)
         ENDIF
      ENDIF
 25   FORMAT( 100(10x,19i4,/) )
      !
   END SUBROUTINE dom_glo


   SUBROUTINE dom_nam
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!                    
      !! ** Purpose :   read domaine namelists and print the variables.
      !!
      !! ** input   : - namrun namelist
      !!              - namdom namelist
      !!              - namnc4 namelist   ! "key_netcdf4" only
      !!----------------------------------------------------------------------
      USE ioipsl
      !!
      INTEGER  ::   ios   ! Local integer
      !
      NAMELIST/namrun/ cn_ocerst_indir, cn_ocerst_outdir, nn_stocklist, ln_rst_list,                 &
         &             nn_no   , cn_exp   , cn_ocerst_in, cn_ocerst_out, ln_rstart , nn_rstctl ,     &
         &             nn_it000, nn_itend , nn_date0    , nn_time0     , nn_leapy  , nn_istate ,     &
         &             nn_stock, nn_write , ln_mskland  , ln_clobber   , nn_chunksz, nn_euler  ,     &
         &             ln_cfmeta, ln_iscpl, ln_xios_read, nn_wxios, ln_rstdate, ln_rst_eos

      NAMELIST/namdom/ ln_linssh, rn_isfhmin, rn_rdt, rn_atfp, ln_crs, ln_meshmask, ln_2d
#if defined key_netcdf4
      NAMELIST/namnc4/ nn_nchunks_i, nn_nchunks_j, nn_nchunks_k, ln_nc4zip
#endif
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_nam : domain initialization through namelist read'
         WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !
      !
      REWIND( numnam_ref )              ! Namelist namrun in reference namelist : Parameters of the run
      READ  ( numnam_ref, namrun, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namrun in reference namelist' )
      REWIND( numnam_cfg )              ! Namelist namrun in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, namrun, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namrun in configuration namelist' )
      IF(lwm) WRITE ( numond, namrun )
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*) '   Namelist : namrun   ---   run parameters'
         WRITE(numout,*) '      Assimilation cycle              nn_no           = ', nn_no
         WRITE(numout,*) '      experiment name for output      cn_exp          = ', TRIM( cn_exp           )
         WRITE(numout,*) '      file prefix restart input       cn_ocerst_in    = ', TRIM( cn_ocerst_in     )
         WRITE(numout,*) '      restart input directory         cn_ocerst_indir = ', TRIM( cn_ocerst_indir  )
         WRITE(numout,*) '      file prefix restart output      cn_ocerst_out   = ', TRIM( cn_ocerst_out    )
         WRITE(numout,*) '      restart output directory        cn_ocerst_outdir= ', TRIM( cn_ocerst_outdir )
         WRITE(numout,*) '      restart logical                 ln_rstart       = ', ln_rstart
         WRITE(numout,*) '      start with forward time step    nn_euler        = ', nn_euler
         WRITE(numout,*) '      control of time step            nn_rstctl       = ', nn_rstctl
         WRITE(numout,*) '      number of the first time step   nn_it000        = ', nn_it000
         WRITE(numout,*) '      number of the last time step    nn_itend        = ', nn_itend
         WRITE(numout,*) '      initial calendar date aammjj    nn_date0        = ', nn_date0
         WRITE(numout,*) '      initial time of day in hhmm     nn_time0        = ', nn_time0
         WRITE(numout,*) '      leap year calendar (0/1)        nn_leapy        = ', nn_leapy
         WRITE(numout,*) '      initial state output            nn_istate       = ', nn_istate
         IF( ln_rst_list ) THEN
            WRITE(numout,*) '      list of restart dump times      nn_stocklist    =', nn_stocklist
         ELSE
            WRITE(numout,*) '      frequency of restart file       nn_stock        = ', nn_stock
         ENDIF
#if ! defined key_iomput
         WRITE(numout,*) '      frequency of output file        nn_write        = ', nn_write
#endif
         WRITE(numout,*) '      mask land points                ln_mskland      = ', ln_mskland
         WRITE(numout,*) '      date-stamp restart files        ln_rstdate = ', ln_rstdate
         WRITE(numout,*) '      additional CF standard metadata ln_cfmeta       = ', ln_cfmeta
         WRITE(numout,*) '      overwrite an existing file      ln_clobber      = ', ln_clobber
         WRITE(numout,*) '      NetCDF chunksize (bytes)        nn_chunksz      = ', nn_chunksz
         WRITE(numout,*) '      IS coupling at the restart step ln_iscpl        = ', ln_iscpl
         WRITE(numout,*) '      check restart equation of state ln_rst_eos      = ', ln_rst_eos

         IF( TRIM(Agrif_CFixed()) == '0' ) THEN
            WRITE(numout,*) '      READ restart for a single file using XIOS ln_xios_read =', ln_xios_read
            WRITE(numout,*) '      Write restart using XIOS        nn_wxios   = ', nn_wxios
         ELSE
            WRITE(numout,*) "      AGRIF: nn_wxios will be ingored. See setting for parent"
            WRITE(numout,*) "      AGRIF: ln_xios_read will be ingored. See setting for parent"
         ENDIF
      ENDIF

      cexper = cn_exp         ! conversion DOCTOR names into model names (this should disappear soon)
      nrstdt = nn_rstctl
      nit000 = nn_it000
      nitend = nn_itend
      ndate0 = nn_date0
      nleapy = nn_leapy
      ninist = nn_istate
      neuler = nn_euler
      IF( neuler == 1 .AND. .NOT. ln_rstart ) THEN
         IF(lwp) WRITE(numout,*)  
         IF(lwp) WRITE(numout,*)'   ==>>>   Start from rest (ln_rstart=F)'
         IF(lwp) WRITE(numout,*)'           an Euler initial time step is used : nn_euler is forced to 0 '   
         neuler = 0
      ENDIF
      !                             ! control of output frequency
      IF( .NOT. ln_rst_list ) THEN     ! we use nn_stock
         IF( nn_stock == -1 )   CALL ctl_warn( 'nn_stock = -1 --> no restart will be done' )
         IF( nn_stock == 0 .OR. nn_stock > nitend ) THEN
            WRITE(ctmp1,*) 'nn_stock = ', nn_stock, ' it is forced to ', nitend
            CALL ctl_warn( ctmp1 )
            nn_stock = nitend
         ENDIF
      ENDIF
#if ! defined key_iomput
      IF( nn_write == -1 )   CALL ctl_warn( 'nn_write = -1 --> no output files will be done' )
      IF ( nn_write == 0 ) THEN
         WRITE(ctmp1,*) 'nn_write = ', nn_write, ' it is forced to ', nitend
         CALL ctl_warn( ctmp1 )
         nn_write = nitend
      ENDIF
#endif

#if defined key_agrif
      IF( Agrif_Root() ) THEN
#endif
      IF(lwp) WRITE(numout,*)
      SELECT CASE ( nleapy )        ! Choose calendar for IOIPSL
      CASE (  1 ) 
         CALL ioconf_calendar('gregorian')
         IF(lwp) WRITE(numout,*) '   ==>>>   The IOIPSL calendar is "gregorian", i.e. leap year'
      CASE (  0 )
         CALL ioconf_calendar('noleap')
         IF(lwp) WRITE(numout,*) '   ==>>>   The IOIPSL calendar is "noleap", i.e. no leap year'
      CASE ( 30 )
         CALL ioconf_calendar('360d')
         IF(lwp) WRITE(numout,*) '   ==>>>   The IOIPSL calendar is "360d", i.e. 360 days in a year'
      END SELECT
#if defined key_agrif
      ENDIF
#endif

      REWIND( numnam_ref )              ! Namelist namdom in reference namelist : space & time domain (bathymetry, mesh, timestep)
      READ  ( numnam_ref, namdom, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdom in reference namelist' )
      REWIND( numnam_cfg )              ! Namelist namdom in configuration namelist : space & time domain (bathymetry, mesh, timestep)
      READ  ( numnam_cfg, namdom, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdom in configuration namelist' )
      IF(lwm) WRITE( numond, namdom )
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist : namdom   ---   space & time domain'
         WRITE(numout,*) '      linear free surface (=T)                ln_linssh   = ', ln_linssh
         WRITE(numout,*) '      create mesh/mask file                   ln_meshmask = ', ln_meshmask
         WRITE(numout,*) '      treshold to open the isf cavity         rn_isfhmin  = ', rn_isfhmin, ' [m]'
         WRITE(numout,*) '      ocean time step                         rn_rdt      = ', rn_rdt
         WRITE(numout,*) '      asselin time filter parameter           rn_atfp     = ', rn_atfp
         WRITE(numout,*) '      online coarsening of dynamical fields   ln_crs      = ', ln_crs
         WRITE(numout,*) '      2D mode                                 ln_2d      = ', ln_2d 
         IF(ln_2d) WRITE(numout,*) '   2D mode active: All tracer processes and vertical diffusion turned off'      
      ENDIF
      !
      !          ! conversion DOCTOR names into model names (this should disappear soon)
      atfp = rn_atfp
      rdt  = rn_rdt

      IF( TRIM(Agrif_CFixed()) == '0' ) THEN
         lrxios = ln_xios_read.AND.ln_rstart
!set output file type for XIOS based on NEMO namelist 
         IF (nn_wxios > 0) lwxios = .TRUE. 
         nxioso = nn_wxios
      ENDIF

#if defined key_netcdf4
      !                             ! NetCDF 4 case   ("key_netcdf4" defined)
      REWIND( numnam_ref )              ! Namelist namnc4 in reference namelist : NETCDF
      READ  ( numnam_ref, namnc4, IOSTAT = ios, ERR = 907)
907   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namnc4 in reference namelist' )
      REWIND( numnam_cfg )              ! Namelist namnc4 in configuration namelist : NETCDF
      READ  ( numnam_cfg, namnc4, IOSTAT = ios, ERR = 908 )
908   IF( ios >  0 )   CALL ctl_nam ( ios , 'namnc4 in configuration namelist' )
      IF(lwm) WRITE( numond, namnc4 )

      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namnc4 - Netcdf4 chunking parameters'
         WRITE(numout,*) '      number of chunks in i-dimension             nn_nchunks_i = ', nn_nchunks_i
         WRITE(numout,*) '      number of chunks in j-dimension             nn_nchunks_j = ', nn_nchunks_j
         WRITE(numout,*) '      number of chunks in k-dimension             nn_nchunks_k = ', nn_nchunks_k
         WRITE(numout,*) '      apply netcdf4/hdf5 chunking & compression   ln_nc4zip    = ', ln_nc4zip
      ENDIF

      ! Put the netcdf4 settings into a simple structure (snc4set, defined in in_out_manager module)
      ! Note the chunk size in the unlimited (time) dimension will be fixed at 1
      snc4set%ni   = nn_nchunks_i
      snc4set%nj   = nn_nchunks_j
      snc4set%nk   = nn_nchunks_k
      snc4set%luse = ln_nc4zip
#else
      snc4set%luse = .FALSE.        ! No NetCDF 4 case
#endif
      !
   END SUBROUTINE dom_nam


   SUBROUTINE dom_ctl
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_ctl  ***
      !!
      !! ** Purpose :   Domain control.
      !!
      !! ** Method  :   compute and print extrema of masked scale factors
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(2) ::   imi1, imi2, ima1, ima2
      INTEGER, DIMENSION(2) ::   iloc   ! 
      REAL(wp) ::   ze1min, ze1max, ze2min, ze2max
      !!----------------------------------------------------------------------
      !
      IF(lk_mpp) THEN
         CALL mpp_minloc( 'domain', e1t(:,:), tmask_i(:,:), ze1min, imi1 )
         CALL mpp_minloc( 'domain', e2t(:,:), tmask_i(:,:), ze2min, imi2 )
         CALL mpp_maxloc( 'domain', e1t(:,:), tmask_i(:,:), ze1max, ima1 )
         CALL mpp_maxloc( 'domain', e2t(:,:), tmask_i(:,:), ze2max, ima2 )
      ELSE
         ze1max = MAXVAL( e1t(:,:)*tmask_i(:,:) )    
         ze2max = MAXVAL( e2t(:,:)*tmask_i(:,:) )    
         ze1min = ze1max - MAXVAL( (ze1max-e1t(:,:))*tmask_i(:,:) )    
         ze2min = ze2max - MAXVAL( (ze2max-e2t(:,:))*tmask_i(:,:) )    
         !
         iloc  = MAXLOC( (ze1max-e1t(:,:))*tmask_i(:,:) )
         imi1(1) = iloc(1) + nimpp - 1
         imi1(2) = iloc(2) + njmpp - 1
         iloc  = MAXLOC( (ze2max-e2t(:,:))*tmask_i(:,:) )
         imi2(1) = iloc(1) + nimpp - 1
         imi2(2) = iloc(2) + njmpp - 1
         iloc  = MAXLOC( e1t(:,:)*tmask_i(:,:) )
         ima1(1) = iloc(1) + nimpp - 1
         ima1(2) = iloc(2) + njmpp - 1
         iloc  = MAXLOC( e2t(:,:)*tmask_i(:,:) )
         ima2(1) = iloc(1) + nimpp - 1
         ima2(2) = iloc(2) + njmpp - 1
      ENDIF
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_ctl : extrema of the masked scale factors'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,"(14x,'e1t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1max, ima1(1), ima1(2)
         WRITE(numout,"(14x,'e1t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1min, imi1(1), imi1(2)
         WRITE(numout,"(14x,'e2t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2max, ima2(1), ima2(2)
         WRITE(numout,"(14x,'e2t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2min, imi2(1), imi2(2)
      ENDIF
      !
   END SUBROUTINE dom_ctl


   SUBROUTINE domain_cfg( cd_cfg, kk_cfg, kpi, kpj, kpk, kperio )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!                    
      !! ** Purpose :   read the domain size in domain configuration file
      !!
      !! ** Method  :   read the cn_domcfg NetCDF file
      !!----------------------------------------------------------------------
      CHARACTER(len=*)              , INTENT(out) ::   cd_cfg          ! configuration name
      INTEGER                       , INTENT(out) ::   kk_cfg          ! configuration resolution
      INTEGER                       , INTENT(out) ::   kpi, kpj, kpk   ! global domain sizes 
      INTEGER                       , INTENT(out) ::   kperio          ! lateral global domain b.c. 
      !
      INTEGER ::   inum   ! local integer
      REAL(wp) ::   zorca_res                     ! local scalars
      REAL(wp) ::   zperio                        !   -      -
      INTEGER, DIMENSION(4) ::   idvar, idimsz    ! size   of dimensions
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*) '           '
         WRITE(numout,*) 'domain_cfg : domain size read in ', TRIM( cn_domcfg ), ' file'
         WRITE(numout,*) '~~~~~~~~~~ '
      ENDIF
      !
      CALL iom_open( cn_domcfg, inum )
      !
      !                                   !- ORCA family specificity
      IF(  iom_varid( inum, 'ORCA'       , ldstop = .FALSE. ) > 0  .AND.  &
         & iom_varid( inum, 'ORCA_index' , ldstop = .FALSE. ) > 0    ) THEN
         !
         cd_cfg = 'ORCA'
         CALL iom_get( inum, 'ORCA_index', zorca_res )   ;   kk_cfg = NINT( zorca_res )
         !
         IF(lwp) THEN
            WRITE(numout,*) '   .'
            WRITE(numout,*) '   ==>>>   ORCA configuration '
            WRITE(numout,*) '   .'
         ENDIF
         !
      ELSE                                !- cd_cfg & k_cfg are not used
         cd_cfg = 'UNKNOWN'
         kk_cfg = -9999999
                                          !- or they may be present as global attributes 
                                          !- (netcdf only)  
         CALL iom_getatt( inum, 'cn_cfg', cd_cfg )  ! returns   !  if not found
         CALL iom_getatt( inum, 'nn_cfg', kk_cfg )  ! returns -999 if not found
         IF( TRIM(cd_cfg) == '!') cd_cfg = 'UNKNOWN'
         IF( kk_cfg == -999     ) kk_cfg = -9999999
         !
      ENDIF
       !
      idvar = iom_varid( inum, 'e3t_0', kdimsz = idimsz )   ! use e3t_0, that must exist, to get jp(ijk)glo
      kpi = idimsz(1)
      kpj = idimsz(2)
      kpk = idimsz(3)
      CALL iom_get( inum, 'jperio', zperio )   ;   kperio = NINT( zperio )
      CALL iom_close( inum )
      !
      IF(lwp) THEN
         WRITE(numout,*) '      cn_cfg = ', TRIM(cd_cfg), '   nn_cfg = ', kk_cfg
         WRITE(numout,*) '      jpiglo = ', kpi
         WRITE(numout,*) '      jpjglo = ', kpj
         WRITE(numout,*) '      jpkglo = ', kpk
         WRITE(numout,*) '      type of global domain lateral boundary   jperio = ', kperio
      ENDIF
      !        
   END SUBROUTINE domain_cfg
   
   
   SUBROUTINE cfg_write
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE cfg_write  ***
      !!                   
      !! ** Purpose :   Create the "cn_domcfg_out" file, a NetCDF file which 
      !!              contains all the ocean domain informations required to 
      !!              define an ocean configuration.
      !!
      !! ** Method  :   Write in a file all the arrays required to set up an
      !!              ocean configuration.
      !!
      !! ** output file :   domcfg_out.nc : domain size, characteristics, horizontal 
      !!                       mesh, Coriolis parameter, and vertical scale factors
      !!                    NB: also contain ORCA family information
      !!----------------------------------------------------------------------
      INTEGER           ::   ji, jj, jk   ! dummy loop indices
      INTEGER           ::   izco, izps, isco, icav
      INTEGER           ::   inum     ! local units
      CHARACTER(len=21) ::   clnam    ! filename (mesh and mask informations)
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! workspace
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'cfg_write : create the domain configuration file (', TRIM(cn_domcfg_out),'.nc)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~'
      !
      !                       ! ============================= !
      !                       !  create 'domcfg_out.nc' file  !
      !                       ! ============================= !
      !         
      clnam = cn_domcfg_out  ! filename (configuration information)
      CALL iom_open( TRIM(clnam), inum, ldwrt = .TRUE. )
      
      !
      !                             !==  ORCA family specificities  ==!
      IF( TRIM(cn_cfg) == "orca" .OR. TRIM(cn_cfg) == "ORCA" ) THEN
         CALL iom_rstput( 0, 0, inum, 'ORCA'      , 1._wp            , ktype = jp_i4 )
         CALL iom_rstput( 0, 0, inum, 'ORCA_index', REAL( nn_cfg, wp), ktype = jp_i4 )         
      ENDIF
      !
      !                             !==  global domain size  ==!
      !
      CALL iom_rstput( 0, 0, inum, 'jpiglo', REAL( jpiglo, wp), ktype = jp_i4 )
      CALL iom_rstput( 0, 0, inum, 'jpjglo', REAL( jpjglo, wp), ktype = jp_i4 )
      CALL iom_rstput( 0, 0, inum, 'jpkglo', REAL( jpk   , wp), ktype = jp_i4 )
      !
      !                             !==  domain characteristics  ==!
      !
      !                                   ! lateral boundary of the global domain
      CALL iom_rstput( 0, 0, inum, 'jperio', REAL( jperio, wp), ktype = jp_i4 )
      !
      !                                   ! type of vertical coordinate
      IF( ln_zco    ) THEN   ;   izco = 1   ;   ELSE   ;   izco = 0   ;   ENDIF
      IF( ln_zps    ) THEN   ;   izps = 1   ;   ELSE   ;   izps = 0   ;   ENDIF
      IF( ln_sco    ) THEN   ;   isco = 1   ;   ELSE   ;   isco = 0   ;   ENDIF
      CALL iom_rstput( 0, 0, inum, 'ln_zco'   , REAL( izco, wp), ktype = jp_i4 )
      CALL iom_rstput( 0, 0, inum, 'ln_zps'   , REAL( izps, wp), ktype = jp_i4 )
      CALL iom_rstput( 0, 0, inum, 'ln_sco'   , REAL( isco, wp), ktype = jp_i4 )
      !
      !                                   ! ocean cavities under iceshelves
      IF( ln_isfcav ) THEN   ;   icav = 1   ;   ELSE   ;   icav = 0   ;   ENDIF
      CALL iom_rstput( 0, 0, inum, 'ln_isfcav', REAL( icav, wp), ktype = jp_i4 )
      !
      !                             !==  horizontal mesh  !
      !
      CALL iom_rstput( 0, 0, inum, 'glamt', glamt, ktype = jp_r8 )   ! latitude
      CALL iom_rstput( 0, 0, inum, 'glamu', glamu, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'glamv', glamv, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'glamf', glamf, ktype = jp_r8 )
      !                                
      CALL iom_rstput( 0, 0, inum, 'gphit', gphit, ktype = jp_r8 )   ! longitude
      CALL iom_rstput( 0, 0, inum, 'gphiu', gphiu, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'gphiv', gphiv, ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'gphif', gphif, ktype = jp_r8 )
      !                                
      CALL iom_rstput( 0, 0, inum, 'e1t'  , e1t  , ktype = jp_r8 )   ! i-scale factors (e1.)
      CALL iom_rstput( 0, 0, inum, 'e1u'  , e1u  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e1v'  , e1v  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e1f'  , e1f  , ktype = jp_r8 )
      !
      CALL iom_rstput( 0, 0, inum, 'e2t'  , e2t  , ktype = jp_r8 )   ! j-scale factors (e2.)
      CALL iom_rstput( 0, 0, inum, 'e2u'  , e2u  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e2v'  , e2v  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e2f'  , e2f  , ktype = jp_r8 )
      !
      CALL iom_rstput( 0, 0, inum, 'ff_f' , ff_f , ktype = jp_r8 )   ! coriolis factor
      CALL iom_rstput( 0, 0, inum, 'ff_t' , ff_t , ktype = jp_r8 )
      !
      !                             !==  vertical mesh  ==!
      !                                                     
      CALL iom_rstput( 0, 0, inum, 'e3t_1d'  , e3t_1d , ktype = jp_r8 )   ! reference 1D-coordinate
      CALL iom_rstput( 0, 0, inum, 'e3w_1d'  , e3w_1d , ktype = jp_r8 )
      !
      CALL iom_rstput( 0, 0, inum, 'e3t_0'   , e3t_0  , ktype = jp_r8 )   ! vertical scale factors
      CALL iom_rstput( 0, 0, inum, 'e3u_0'   , e3u_0  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3v_0'   , e3v_0  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3f_0'   , e3f_0  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3w_0'   , e3w_0  , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3uw_0'  , e3uw_0 , ktype = jp_r8 )
      CALL iom_rstput( 0, 0, inum, 'e3vw_0'  , e3vw_0 , ktype = jp_r8 )
      !                                         
      !                             !==  wet top and bottom level  ==!   (caution: multiplied by ssmask)
      !
      CALL iom_rstput( 0, 0, inum, 'top_level'    , REAL( mikt, wp )*ssmask , ktype = jp_i4 )   ! nb of ocean T-points (ISF)
      CALL iom_rstput( 0, 0, inum, 'bottom_level' , REAL( mbkt, wp )*ssmask , ktype = jp_i4 )   ! nb of ocean T-points
      !
      IF( ln_sco ) THEN             ! s-coordinate: store grid stiffness ratio  (Not required anyway)
         CALL dom_stiff( z2d )
         CALL iom_rstput( 0, 0, inum, 'stiffness', z2d )        !    ! Max. grid stiffness ratio
      ENDIF
      !
      IF( ll_wd ) THEN              ! wetting and drying domain
         CALL iom_rstput( 0, 0, inum, 'ht_0'   , ht_0   , ktype = jp_r8 )
      ENDIF
      !
      ! Add some global attributes ( netcdf only )
      CALL iom_putatt( inum, 'nn_cfg', nn_cfg )
      CALL iom_putatt( inum, 'cn_cfg', TRIM(cn_cfg) )
      !
      !                                ! ============================
      !                                !        close the files 
      !                                ! ============================
      CALL iom_close( inum )
      !
   END SUBROUTINE cfg_write

   !!======================================================================
END MODULE domain
