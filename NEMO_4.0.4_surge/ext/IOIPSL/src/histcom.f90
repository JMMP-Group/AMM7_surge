MODULE histcom
!-
!$Id: histcom.f90 2368 2010-11-09 15:38:45Z acc $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!-
  USE netcdf
  USE nc4interface ! needed to allow compilation with netcdf3 libraries
!-
  USE stringop, ONLY : nocomma,cmpblank,findpos,find_str,strlowercase
  USE mathelp,  ONLY : mathop,moycum,buildop
  USE fliocom,  ONLY : flio_dom_file,flio_dom_att
  USE calendar
  USE errioipsl, ONLY : ipslerr,ipsldbg
!-
  IMPLICIT NONE
!-
  PRIVATE
  PUBLIC :: histbeg,histdef,histhori,histvert,histend, &
 &          histwrite,histclo,histsync,ioconf_modname
!---------------------------------------------------------------------
!- Some confusing vocabulary in this code !
!- =========================================
!-
!- A REGULAR grid is a grid which is i,j indexes
!- and thus it is stored in a 2D matrix.
!- This is opposed to a IRREGULAR grid which is only in a vector
!- and where we do not know which neighbors we have.
!- As a consequence we need the bounds for each grid-cell.
!-
!- A RECTILINEAR grid is a special case of a regular grid
!- in which all longitudes for i constant are equal
!- and all latitudes for j constant.
!- In other words we do not need the full 2D matrix
!- to describe the grid, just two vectors.
!---------------------------------------------------------------------
!-
  INTERFACE histbeg
    MODULE PROCEDURE histb_reg1d,histb_reg2d,histb_irreg
  END INTERFACE
!-
  INTERFACE histhori
    MODULE PROCEDURE histh_reg1d,histh_reg2d,histh_irreg
  END INTERFACE
!-
  INTERFACE histwrite
!---------------------------------------------------------------------
!- The "histwrite" routines will give the data to the I/O system.
!- It will trigger the operations to be performed,
!- and the writting to the file if needed
!-
!- We test for the work to be done at this time here so that at a
!- later stage we can call different operation and write subroutine
!- for the REAL and INTEGER interfaces
!-
!- INPUT
!- idf      : The ID of the file on which this variable is to be,
!-            written. The variable should have been defined in
!-            this file before.
!- pvarname : The short name of the variable
!- pitau    : Current timestep
!- pdata    : The variable, I mean the real data !
!- nbindex  : The number of indexes provided. If it is equal to
!-            the size of the full field as provided in histdef
!-            then nothing is done.
!- nindex   : The indices used to expand the variable (pdata)
!-            onto the full field.
!---------------------------------------------------------------------
!- histwrite - we have to prepare different type of fields :
!-             real and integer, 1,2 or 3D
    MODULE PROCEDURE histwrite_r1d,histwrite_r2d,histwrite_r3d
  END INTERFACE
!-
! Fixed parameter
!-
  INTEGER,PARAMETER :: nb_files_max=20,nb_var_max=400, &
 &                     nb_hax_max=5,nb_zax_max=10,nbopp_max=10
  REAL,PARAMETER :: missing_val=nf90_fill_real
  INTEGER,PARAMETER,PUBLIC :: &
 &  hist_r4=nf90_real4, hist_r8=nf90_real8
!-
! Variable derived type
!-
TYPE T_D_V
  INTEGER :: ncvid
  INTEGER :: nbopp
  CHARACTER(LEN=20)  :: v_name,unit_name
  CHARACTER(LEN=256) :: title,std_name
  CHARACTER(LEN=80)  :: fullop
  CHARACTER(LEN=7)   :: topp
  CHARACTER(LEN=7),DIMENSION(nbopp_max) :: sopp
  REAL,DIMENSION(nbopp_max) :: scal
!-External type (for R4/R8)
  INTEGER :: v_typ
!-Sizes of the associated grid and zommed area
  INTEGER,DIMENSION(3) :: scsize,zorig,zsize
!-Sizes for the data as it goes through the various math operations
  INTEGER,DIMENSION(3) :: datasz_in = -1
  INTEGER :: datasz_max = -1
!-
  INTEGER :: h_axid,z_axid,t_axid
!-
  REAL,DIMENSION(2) :: hist_minmax
  LOGICAL :: hist_calc_rng=.FALSE.,hist_wrt_rng=.FALSE.
!-Book keeping of the axes
  INTEGER :: tdimid,tbndid=-1,tax_last
  LOGICAL :: l_bnd
  CHARACTER(LEN=40) :: tax_name
!-
  REAL :: freq_opp,freq_wrt
  INTEGER :: &
 &  last_opp,last_wrt,last_opp_chk,last_wrt_chk,nb_opp,nb_wrt
!- For future optimization
  REAL,POINTER,DIMENSION(:) :: t_bf
!#  REAL,ALLOCATABLE,DIMENSION(:) :: V_1_D
!#  REAL,ALLOCATABLE,DIMENSION(:,:) :: V_2_D
!#  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: V_3_D
END TYPE T_D_V
!-
! File derived type
!-
TYPE :: T_D_F
!-NETCDF IDs for file
  INTEGER :: ncfid=-1
!-Time variables
  INTEGER :: itau0=0
  REAL :: date0,deltat
!-Counter of elements (variables, time-horizontal-vertical axis
  INTEGER :: n_var=0,n_tax=0,n_hax=0,n_zax=0
!-NETCDF dimension IDs for time-[time_bounds]-longitude-latitude
  INTEGER :: tid,bid,xid,yid
!-General definitions in the NETCDF file
  INTEGER,DIMENSION(2) :: full_size=0,slab_ori,slab_siz
!-The horizontal axes
  CHARACTER(LEN=25),DIMENSION(nb_hax_max,2) :: hax_name
!-The vertical axes
  INTEGER,DIMENSION(nb_zax_max) :: zax_size,zax_ids
  CHARACTER(LEN=20),DIMENSION(nb_zax_max) :: zax_name
!-
  LOGICAL :: regular=.TRUE.
!-DOMAIN ID
  INTEGER :: dom_id_svg=-1
!-
  TYPE(T_D_V),DIMENSION(nb_var_max) :: W_V
END TYPE T_D_F
!-
TYPE(T_D_F),DIMENSION(nb_files_max),SAVE :: W_F
!-
! A list of functions which require special action
! (Needs to be updated when functions are added
!  but they are well located here)
!-
  CHARACTER(LEN=30),SAVE :: fuchnbout = 'scatter, fill'
!- Some configurable variables with locks
  CHARACTER(LEN=80),SAVE :: model_name='An IPSL model'
  LOGICAL,SAVE :: lock_modname=.FALSE.
!-
!===
CONTAINS
!===
!-
SUBROUTINE histb_reg1d                 &
 & (pfilename,pim,plon,pjm,plat,       &
 &  par_orix,par_szx,par_oriy,par_szy, &
 &  pitau0,pdate0,pdeltat,phoriid,idf,domain_id,mode,snc4chunks)
!---------------------------------------------------------------------
!- histbeg for 1D regular horizontal coordinates (see histb_all)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: pfilename
  INTEGER,INTENT(IN) :: pim,pjm
  REAL,DIMENSION(pim),INTENT(IN) :: plon
  REAL,DIMENSION(pjm),INTENT(IN) :: plat
  INTEGER,INTENT(IN):: par_orix,par_szx,par_oriy,par_szy
  INTEGER,INTENT(IN) :: pitau0
  REAL,INTENT(IN) :: pdate0,pdeltat
  INTEGER,INTENT(OUT) :: idf,phoriid
  INTEGER,INTENT(IN),OPTIONAL :: domain_id
  CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: mode
  TYPE(snc4_ctl), OPTIONAL,INTENT(IN) :: snc4chunks
!---------------------------------------------------------------------
  CALL histb_all &
 & (1,pfilename,pim,pjm,pitau0,pdate0,pdeltat,phoriid,idf, &
 &  x_1d=plon,y_1d=plat, &
 &  k_orx=par_orix,k_szx=par_szx,k_ory=par_oriy,k_szy=par_szy, &
 &  domain_id=domain_id,mode=mode,snc4chunks=snc4chunks)
!-------------------------
END SUBROUTINE histb_reg1d
!===
SUBROUTINE histb_reg2d &
 & (pfilename,pim,plon,pjm,plat,       &
 &  par_orix,par_szx,par_oriy,par_szy, &
 &  pitau0,pdate0,pdeltat,phoriid,idf,domain_id,mode,snc4chunks)
!---------------------------------------------------------------------
!- histbeg for 2D regular horizontal coordinates (see histb_all)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: pfilename
  INTEGER,INTENT(IN) :: pim,pjm
  REAL,DIMENSION(pim,pjm),INTENT(IN) :: plon,plat
  INTEGER,INTENT(IN):: par_orix,par_szx,par_oriy,par_szy
  INTEGER,INTENT(IN) :: pitau0
  REAL,INTENT(IN) :: pdate0,pdeltat
  INTEGER,INTENT(OUT) :: idf,phoriid
  INTEGER,INTENT(IN),OPTIONAL :: domain_id
  CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: mode
  TYPE(snc4_ctl), OPTIONAL,INTENT(IN) :: snc4chunks
!---------------------------------------------------------------------
  CALL histb_all &
 & (2,pfilename,pim,pjm,pitau0,pdate0,pdeltat,phoriid,idf,  &
 &  x_2d=plon,y_2d=plat, &
 &  k_orx=par_orix,k_szx=par_szx,k_ory=par_oriy,k_szy=par_szy,    &
 &  domain_id=domain_id,mode=mode,snc4chunks=snc4chunks)
!-------------------------
END SUBROUTINE histb_reg2d
!===
SUBROUTINE histb_irreg &
 &  (pfilename,pim,plon,plon_bounds,plat,plat_bounds, &
 &   pitau0,pdate0,pdeltat,phoriid,idf,domain_id,mode,snc4chunks)
!---------------------------------------------------------------------
!- histbeg for irregular horizontal coordinates (see histb_all)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: pfilename
  INTEGER,INTENT(IN) :: pim
  REAL,DIMENSION(pim),INTENT(IN) :: plon,plat
  REAL,DIMENSION(:,:),INTENT(IN) :: plon_bounds,plat_bounds
  INTEGER,INTENT(IN) :: pitau0
  REAL,INTENT(IN) :: pdate0,pdeltat
  INTEGER,INTENT(OUT) :: idf,phoriid
  INTEGER,INTENT(IN),OPTIONAL :: domain_id
  CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: mode
  TYPE(snc4_ctl), OPTIONAL,INTENT(IN) :: snc4chunks
!---------------------------------------------------------------------
  CALL histb_all &
 & (3,pfilename,pim,pim,pitau0,pdate0,pdeltat,phoriid,idf,  &
 &  x_1d=plon,y_1d=plat,x_bnds=plon_bounds,y_bnds=plat_bounds, &
 &  domain_id=domain_id,mode=mode,snc4chunks=snc4chunks)
!-------------------------
END SUBROUTINE histb_irreg
!===
SUBROUTINE histb_all &
 & (k_typ,nc_name,pim,pjm,pitau0,pdate0,pdeltat,phoriid,idf, &
 &  x_1d,y_1d,x_2d,y_2d,k_orx,k_szx,k_ory,k_szy, &
 &  x_bnds,y_bnds,domain_id,mode,snc4chunks)
!---------------------------------------------------------------------
!- General interface for horizontal grids.
!- This subroutine initializes a netcdf file and returns the ID.
!- It will set up the geographical space on which the data will be
!- stored and offers the possibility of seting a zoom.
!- In the case of irregular grids, all the data comes in as vectors
!- and for the grid we have the coordinates of the 4 corners.
!- It also gets the global parameters into the I/O subsystem.
!-
!- INPUT
!-
!- k_typ    : Type of the grid (1 rectilinear, 2 regular, 3 irregular)
!- nc_name  : Name of the netcdf file to be created
!- pim      : Size of arrays in longitude direction
!- pjm      : Size of arrays in latitude direction (pjm=pim for type 3)
!-
!- pitau0   : time step at which the history tape starts
!- pdate0   : The Julian date at which the itau was equal to 0
!- pdeltat  : Time step, in seconds, of the counter itau
!-            used in histwrite for instance
!-
!- OUTPUT
!-
!- phoriid  : Identifier of the horizontal grid
!- idf      : Identifier of the file
!-
!- Optional INPUT arguments
!-
!- For rectilinear or irregular grid
!- x_1d  : The longitudes
!- y_1d  : The latitudes
!- For regular grid
!- x_2d  : The longitudes
!- y_2d  : The latitudes
!- One pair (x_1d,y_1d) or (x_2d,y_2d) must be supplied.
!-
!- For regular grid (reg1d or reg2d),
!- the next 4 arguments allow to define a horizontal zoom
!- for this file. It is assumed that all variables to come
!- have the same index space. This can not be assumed for
!- the z axis and thus we define the zoom in histdef.
!- k_orx  : Origin of the slab of data within the X axis (pim)
!- k_szx  : Size of the slab of data in X
!- k_ory  : Origin of the slab of data within the Y axis (pjm)
!- k_szy  : Size of the slab of data in Y
!-
!- For irregular grid.
!- x_bnds : The boundaries of the grid in longitude
!- y_bnds : The boundaries of the grid in latitude
!-
!- For all grids.
!-
!- domain_id  : Domain identifier
!-
!- mode       : String of (case insensitive) blank-separated words
!-              defining the mode used to create the file.
!-              Supported keywords : 32, 64
!-              "32/64" defines the offset mode.
!-              The default offset mode is 64 bits.
!-              Keywords "NETCDF4" and "CLASSIC" are reserved
!-              for future use.
!-
!- snc4chunks    : Structure containing chunk partitioning parameters
!-              for 4-D variables and a logical switch to toggle
!-              between netcdf3 o/p (false) and netcdf4 chunked
!-              and compressed o/p (true)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: k_typ
  CHARACTER(LEN=*),INTENT(IN) :: nc_name
  INTEGER,INTENT(IN) :: pim,pjm
  INTEGER,INTENT(IN) :: pitau0
  REAL,INTENT(IN) :: pdate0,pdeltat
  INTEGER,INTENT(OUT) :: idf,phoriid
  REAL,DIMENSION(:),INTENT(IN),OPTIONAL :: x_1d,y_1d
  REAL,DIMENSION(:,:),INTENT(IN),OPTIONAL :: x_2d,y_2d
  INTEGER,INTENT(IN),OPTIONAL :: k_orx,k_szx,k_ory,k_szy
  REAL,DIMENSION(:,:),INTENT(IN),OPTIONAL :: x_bnds,y_bnds
  INTEGER,INTENT(IN),OPTIONAL :: domain_id
  CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: mode
  TYPE(snc4_ctl), OPTIONAL,INTENT(IN) :: snc4chunks
!-
  INTEGER :: nfid,iret,m_c
  CHARACTER(LEN=120) :: file
  CHARACTER(LEN=30) :: timenow
  CHARACTER(LEN=11) :: c_nam
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF     (k_typ == 1) THEN
    c_nam = 'histb_reg1d'
  ELSEIF (k_typ == 2) THEN
    c_nam = 'histb_reg2d'
  ELSEIF (k_typ == 3) THEN
    c_nam = 'histb_irreg'
  ELSE
    CALL ipslerr (3,"histbeg", &
 &    'Illegal value of k_typ argument','in internal interface','?')
  ENDIF
!-
  IF (l_dbg) WRITE(*,*) c_nam//" 0.0"
!-
! Search for a free index
!-
  idf = -1
  DO nfid=1,nb_files_max
    IF (W_F(nfid)%ncfid < 0) THEN
      idf = nfid; EXIT;
    ENDIF
  ENDDO
  IF (idf < 0) THEN
    CALL ipslerr (3,"histbeg", &
   &  'Table of files too small. You should increase nb_files_max', &
   &  'in histcom.f90 in order to accomodate all these files',' ')
  ENDIF
!-
! 1.0 Transfering into the common for future use
!-
  IF (l_dbg) WRITE(*,*) c_nam//" 1.0"
!-
  W_F(idf)%itau0  = pitau0
  W_F(idf)%date0  = pdate0
  W_F(idf)%deltat = pdeltat
!-
! 2.0 Initializes all variables for this file
!-
  IF (l_dbg) WRITE(*,*) c_nam//" 2.0"
!-
  W_F(idf)%n_var = 0
  W_F(idf)%n_tax = 0
  W_F(idf)%n_hax = 0
  W_F(idf)%n_zax = 0
!-
  IF ( (k_typ == 1).OR.(k_typ == 2) ) THEN
    W_F(idf)%slab_ori(1:2) = (/ k_orx,k_ory /)
    W_F(idf)%slab_siz(1:2)  = (/ k_szx,k_szy /)
  ELSE
    W_F(idf)%slab_ori(1:2) = (/ 1,1 /)
    W_F(idf)%slab_siz(1:2) = (/ pim,1 /)
  ENDIF
!-
! 3.0 Opening netcdf file and defining dimensions
!-
  IF (l_dbg) WRITE(*,*) c_nam//" 3.0"
!-
! Add DOMAIN number and ".nc" suffix in file name if needed
!-
  file = nc_name
  CALL flio_dom_file (file,domain_id)
!-
! Check the mode
!? See fliocom for HDF4 ????????????????????????????????????????????????
!-
  IF (PRESENT(mode)) THEN
    SELECT CASE (TRIM(mode))
    CASE('32')
      m_c = NF90_CLOBBER
    CASE('64')
      m_c = IOR(NF90_CLOBBER,NF90_64BIT_OFFSET)
    CASE DEFAULT
      CALL ipslerr (3,"histbeg", &
 &      'Invalid argument mode for file :',TRIM(file), &
 &      'Supported values are 32 or 64')
    END SELECT
  ELSE
    m_c = IOR(NF90_CLOBBER,NF90_64BIT_OFFSET)
  ENDIF
!-
  IF (PRESENT(snc4chunks)) THEN
    IF (snc4chunks%luse) CALL get_nf90_symbol("NF90_HDF5", m_c)
  ENDIF
!-
! Create file
!-
  iret = NF90_CREATE(file,m_c,nfid)
!-
  IF     (k_typ == 1) THEN
    iret = NF90_DEF_DIM(nfid,'lon',k_szx,W_F(idf)%xid)
    iret = NF90_DEF_DIM(nfid,'lat',k_szy,W_F(idf)%yid)
  ELSEIF (k_typ == 2) THEN
    iret = NF90_DEF_DIM(nfid,'x',k_szx,W_F(idf)%xid)
    iret = NF90_DEF_DIM(nfid,'y',k_szy,W_F(idf)%yid)
  ELSEIF (k_typ == 3) THEN
    iret = NF90_DEF_DIM(nfid,'x',pim,W_F(idf)%xid)
    W_F(idf)%yid = W_F(idf)%xid
  ENDIF
!-
! 4.0 Declaring the geographical coordinates and other attributes
!-
  IF (l_dbg) WRITE(*,*) c_nam//" 4.0"
!-
  iret = NF90_PUT_ATT(nfid,NF90_GLOBAL,'Conventions','CF-1.1')
  iret = NF90_PUT_ATT(nfid,NF90_GLOBAL,'file_name',TRIM(file))
  iret = NF90_PUT_ATT(nfid,NF90_GLOBAL,'production',TRIM(model_name))
  lock_modname = .TRUE.
  CALL ioget_timestamp (timenow)
  iret = NF90_PUT_ATT(nfid,NF90_GLOBAL,'TimeStamp',TRIM(timenow))
!-
! 5.0 Saving some important information on this file in the common
!-
  IF (l_dbg) WRITE(*,*) c_nam//" 5.0"
!-
  IF (PRESENT(domain_id)) THEN
    W_F(idf)%dom_id_svg = domain_id
  ENDIF
  W_F(idf)%ncfid = nfid
  IF ( (k_typ == 1).OR.(k_typ == 2) ) THEN
    W_F(idf)%full_size(1:2) = (/ pim,pjm /)
    W_F(idf)%regular=.TRUE.
  ELSEIF (k_typ == 3) THEN
    W_F(idf)%full_size(1:2) = (/ pim,1 /)
    W_F(idf)%regular=.FALSE.
  ENDIF
!-
! 6.0 storing the geographical coordinates
!-
  IF     (k_typ == 1) THEN
    CALL histh_all &
 &   (k_typ,idf,pim,pjm,' ','Default grid',phoriid, &
 &    x_1d=x_1d,y_1d=y_1d)
  ELSEIF (k_typ == 2) THEN
    CALL histh_all &
 &   (k_typ,idf,pim,pjm,' ','Default grid',phoriid, &
 &    x_2d=x_2d,y_2d=y_2d)
  ELSEIF (k_typ == 3) THEN
    CALL histh_all &
 &   (k_typ,idf,pim,pim,' ','Default grid',phoriid, &
 &    x_1d=x_1d,y_1d=y_1d,x_bnds=x_bnds,y_bnds=y_bnds)
  ENDIF
!-----------------------
END SUBROUTINE histb_all
!===
SUBROUTINE histh_reg1d &
 &  (idf,pim,plon,pjm,plat,phname,phtitle,phid)
!---------------------------------------------------------------------
!- histhori for 1d regular grid (see histh_all)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: idf,pim,pjm
  REAL,INTENT(IN),DIMENSION(:) :: plon,plat
  CHARACTER(LEN=*),INTENT(IN) :: phname,phtitle
  INTEGER,INTENT(OUT) :: phid
!---------------------------------------------------------------------
  CALL histh_all &
 & (1,idf,pim,pjm,phname,phtitle,phid,x_1d=plon,y_1d=plat)
!-------------------------
END SUBROUTINE histh_reg1d
!===
SUBROUTINE histh_reg2d &
 & (idf,pim,plon,pjm,plat,phname,phtitle,phid)
!---------------------------------------------------------------------
!- histhori for 2d regular grid (see histh_all)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: idf,pim,pjm
  REAL,INTENT(IN),DIMENSION(:,:) :: plon,plat
  CHARACTER(LEN=*),INTENT(IN) :: phname,phtitle
  INTEGER,INTENT(OUT) :: phid
!---------------------------------------------------------------------
  CALL histh_all &
 & (2,idf,pim,pjm,phname,phtitle,phid,x_2d=plon,y_2d=plat)
!-------------------------
END SUBROUTINE histh_reg2d
!===
SUBROUTINE histh_irreg &
 & (idf,pim,plon,plon_bounds,plat,plat_bounds,phname,phtitle,phid)
!---------------------------------------------------------------------
!- histhori for irregular grid (see histh_all)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: idf,pim
  REAL,DIMENSION(:),INTENT(IN) :: plon,plat
  REAL,DIMENSION(:,:),INTENT(IN) :: plon_bounds,plat_bounds
  CHARACTER(LEN=*),INTENT(IN) :: phname,phtitle
  INTEGER,INTENT(OUT) :: phid
!---------------------------------------------------------------------
  CALL histh_all &
 & (3,idf,pim,pim,phname,phtitle,phid, &
 &  x_1d=plon,y_1d=plat,x_bnds=plon_bounds,y_bnds=plat_bounds)
!-------------------------
END SUBROUTINE histh_irreg
!===
SUBROUTINE histh_all &
 & (k_typ,idf,pim,pjm,phname,phtitle,phid, &
 &  x_1d,y_1d,x_2d,y_2d,x_bnds,y_bnds)
!---------------------------------------------------------------------
!- General interface for horizontal grids.
!- This subroutine is made to declare a new horizontal grid.
!- It has to have the same number of points as
!- the original and thus in this routine we will only
!- add two variable (longitude and latitude).
!- Any variable in the file can thus point to this pair
!- through an attribute. This routine is very usefull
!- to allow staggered grids.
!-
!- INPUT
!-
!- k_typ   : Type of the grid (1 rectilinear, 2 regular, 3 irregular)
!- idf     : The id of the file to which the grid should be added
!- pim     : Size in the longitude direction
!- pjm     : Size in the latitude direction (pjm=pim for type 3)
!- phname  : The name of grid
!- phtitle : The title of the grid
!-
!- OUTPUT
!-
!- phid    : Id of the created grid
!-
!- Optional INPUT arguments
!-
!- For rectilinear or irregular grid
!- x_1d  : The longitudes
!- y_1d  : The latitudes
!- For regular grid
!- x_2d  : The longitudes
!- y_2d  : The latitudes
!- One pair (x_1d,y_1d) or (x_2d,y_2d) must be supplied.
!-
!- For irregular grid.
!- x_bnds : The boundaries of the grid in longitude
!- y_bnds : The boundaries of the grid in latitude
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: k_typ
  INTEGER,INTENT(IN) :: idf,pim,pjm
  CHARACTER(LEN=*),INTENT(IN) :: phname,phtitle
  INTEGER,INTENT(OUT) :: phid
  REAL,DIMENSION(:),INTENT(IN),OPTIONAL :: x_1d,y_1d
  REAL,DIMENSION(:,:),INTENT(IN),OPTIONAL :: x_2d,y_2d
  REAL,DIMENSION(:,:),INTENT(IN),OPTIONAL :: x_bnds,y_bnds
!-
  CHARACTER(LEN=25) :: lon_name,lat_name
  CHARACTER(LEN=30) :: lonbound_name,latbound_name
  INTEGER :: i_s,i_e
  INTEGER,DIMENSION(2) :: dims,dims_b
  INTEGER :: nbbounds
  INTEGER :: nlonidb,nlatidb,twoid
  LOGICAL :: transp = .FALSE.
  REAL,ALLOCATABLE,DIMENSION(:,:) :: bounds_trans
  REAL :: wmn,wmx
  INTEGER :: nlonid,nlatid
  INTEGER :: o_x,o_y,s_x,s_y
  INTEGER :: iret,nfid
  CHARACTER(LEN=11) :: c_nam
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF     (k_typ == 1) THEN
    c_nam = 'histh_reg1d'
  ELSEIF (k_typ == 2) THEN
    c_nam = 'histh_reg2d'
  ELSEIF (k_typ == 3) THEN
    c_nam = 'histh_irreg'
  ELSE
    CALL ipslerr (3,"histhori", &
 &    'Illegal value of k_typ argument','in internal interface','?')
  ENDIF
!-
! 1.0 Check that all fits in the buffers
!-
  IF (    (pim /= W_F(idf)%full_size(1)) &
 &    .OR.(W_F(idf)%regular.AND.(pjm /= W_F(idf)%full_size(2)))  &
 &    .OR.(.NOT.W_F(idf)%regular.AND.(W_F(idf)%full_size(2) /= 1)) ) THEN
    CALL ipslerr (3,"histhori", &
 &   'The new horizontal grid does not have the same size', &
 &   'as the one provided to histbeg. This is not yet ', &
 &   'possible in the hist package.')
  ENDIF
!-
! 1.1 Create all the variables needed
!-
  IF (l_dbg) WRITE(*,*) c_nam//" 1.0"
!-
  nfid = W_F(idf)%ncfid
!-
  IF (k_typ == 3) THEN
    IF     (SIZE(x_bnds,DIM=1) == pim) THEN
      nbbounds = SIZE(x_bnds,DIM=2)
      transp = .TRUE.
    ELSEIF (SIZE(x_bnds,DIM=2) == pim) THEN
      nbbounds = SIZE(x_bnds,DIM=1)
      transp = .FALSE.
    ELSE
      CALL ipslerr (3,"histhori", &
 &     'The boundary variable does not have any axis corresponding', &
 &     'to the size of the longitude or latitude variable','.')
    ENDIF
    ALLOCATE(bounds_trans(nbbounds,pim))
    iret = NF90_DEF_DIM(nfid,'nbnd',nbbounds,twoid)
    dims_b(1:2) = (/ twoid,W_F(idf)%xid /)
  ENDIF
!-
  dims(1:2) = (/ W_F(idf)%xid,W_F(idf)%yid /)
!-
  IF     (k_typ == 1) THEN
    IF (W_F(idf)%n_hax == 0) THEN
      lon_name = 'lon'
      lat_name = 'lat'
    ELSE
      lon_name = 'lon_'//TRIM(phname)
      lat_name = 'lat_'//TRIM(phname)
    ENDIF
  ELSEIF (k_typ == 2) THEN
    IF (W_F(idf)%n_hax == 0) THEN
      lon_name = 'nav_lon'
      lat_name = 'nav_lat'
    ELSE
      lon_name = 'nav_lon_'//TRIM(phname)
      lat_name = 'nav_lat_'//TRIM(phname)
    ENDIF
  ELSEIF (k_typ == 3) THEN
    IF (W_F(idf)%n_hax == 0) THEN
      lon_name = 'nav_lon'
      lat_name = 'nav_lat'
    ELSE
      lon_name = 'nav_lon_'//TRIM(phname)
      lat_name = 'nav_lat_'//TRIM(phname)
    ENDIF
    lonbound_name = TRIM(lon_name)//'_bounds'
    latbound_name = TRIM(lat_name)//'_bounds'
  ENDIF
!-
! 1.2 Save the informations
!-
  phid = W_F(idf)%n_hax+1
  W_F(idf)%n_hax = phid
  W_F(idf)%hax_name(phid,1:2) = (/ lon_name,lat_name /)
!-
! 2.0 Longitude
!-
  IF (l_dbg) WRITE(*,*) c_nam//" 2.0"
!-
  i_s = 1;
  IF ( (k_typ == 1).OR.(k_typ == 3) ) THEN
    i_e = 1; wmn = MINVAL(x_1d); wmx = MAXVAL(x_1d);
  ELSEIF (k_typ == 2) THEN
    i_e = 2; wmn = MINVAL(x_2d); wmx = MAXVAL(x_2d);
  ENDIF
  iret = NF90_DEF_VAR(nfid,lon_name,NF90_REAL4,dims(i_s:i_e),nlonid)
  IF (k_typ == 1) THEN
    iret = NF90_PUT_ATT(nfid,nlonid,'axis',"X")
  ENDIF
  iret = NF90_PUT_ATT(nfid,nlonid,'standard_name',"longitude")
  iret = NF90_PUT_ATT(nfid,nlonid,'units',"degrees_east")
  iret = NF90_PUT_ATT(nfid,nlonid,'valid_min',REAL(wmn,KIND=4))
  iret = NF90_PUT_ATT(nfid,nlonid,'valid_max',REAL(wmx,KIND=4))
  iret = NF90_PUT_ATT(nfid,nlonid,'long_name',"Longitude")
  iret = NF90_PUT_ATT(nfid,nlonid,'nav_model',TRIM(phtitle))
!-
  IF (k_typ == 3) THEN
!---
!-- 2.1 Longitude bounds
!---
    iret = NF90_PUT_ATT(nfid,nlonid,'bounds',TRIM(lonbound_name))
    iret = NF90_DEF_VAR(nfid,lonbound_name,NF90_REAL4,dims_b(1:2),nlonidb)
    iret = NF90_PUT_ATT(nfid,nlonidb,'long_name', &
 &          'Boundaries for coordinate variable '//TRIM(lon_name))
  ENDIF
!-
! 3.0 Latitude
!-
  IF (l_dbg) WRITE(*,*) c_nam//" 3.0"
!-
  i_e = 2;
  IF ( (k_typ == 1).OR.(k_typ == 3) ) THEN
    i_s = 2; wmn = MINVAL(y_1d); wmx = MAXVAL(y_1d);
  ELSEIF (k_typ == 2) THEN
    i_s = 1; wmn = MINVAL(y_2d); wmx = MAXVAL(y_2d);
  ENDIF
  iret = NF90_DEF_VAR(nfid,lat_name,NF90_REAL4,dims(i_s:i_e),nlatid)
  IF (k_typ == 1) THEN
    iret = NF90_PUT_ATT(nfid,nlatid,'axis',"Y")
  ENDIF
!-
  iret = NF90_PUT_ATT(nfid,nlatid,'standard_name',"latitude")
  iret = NF90_PUT_ATT(nfid,nlatid,'units',"degrees_north")
  iret = NF90_PUT_ATT(nfid,nlatid,'valid_min',REAL(wmn,KIND=4))
  iret = NF90_PUT_ATT(nfid,nlatid,'valid_max',REAL(wmx,KIND=4))
  iret = NF90_PUT_ATT(nfid,nlatid,'long_name',"Latitude")
  iret = NF90_PUT_ATT(nfid,nlatid,'nav_model',TRIM(phtitle))
!-
  IF (k_typ == 3) THEN
!---
!-- 3.1 Latitude bounds
!---
    iret = NF90_PUT_ATT(nfid,nlatid,'bounds',TRIM(latbound_name))
    iret = NF90_DEF_VAR(nfid,latbound_name,NF90_REAL4,dims_b(1:2),nlatidb)
    iret = NF90_PUT_ATT(nfid,nlatidb,'long_name', &
 &          'Boundaries for coordinate variable '//TRIM(lat_name))
  ENDIF
!-
  iret = NF90_ENDDEF(nfid)
!-
! 4.0 storing the geographical coordinates
!-
  IF (l_dbg) WRITE(*,*) c_nam//" 4.0"
!-
  IF ( (k_typ == 1).OR.(k_typ == 2) ) THEN
    o_x = W_F(idf)%slab_ori(1)
    o_y = W_F(idf)%slab_ori(2)
    s_x = W_F(idf)%slab_siz(1)
    s_y = W_F(idf)%slab_siz(2)
!---
!-- Transfer the longitude and  the latitude
!---
    IF     (k_typ == 1) THEN
      iret = NF90_PUT_VAR(nfid,nlonid,x_1d(o_x:o_x+s_x-1))
      iret = NF90_PUT_VAR(nfid,nlatid,y_1d(o_y:o_y+s_y-1))
    ELSEIF (k_typ == 2) THEN
      iret = NF90_PUT_VAR(nfid,nlonid, &
 &            x_2d(o_x:o_x+s_x-1,o_y:o_y+s_y-1))
      iret = NF90_PUT_VAR(nfid,nlatid, &
 &            y_2d(o_x:o_x+s_x-1,o_y:o_y+s_y-1))
    ENDIF
  ELSEIF (k_typ == 3) THEN
!---
!-- Transfer the longitude and the longitude bounds
!---
    iret = NF90_PUT_VAR(nfid,nlonid,x_1d(1:pim))
!---
    IF (transp) THEN
      bounds_trans = TRANSPOSE(x_bnds)
    ELSE
      bounds_trans = x_bnds
    ENDIF
    iret = NF90_PUT_VAR(nfid,nlonidb,bounds_trans(1:nbbounds,1:pim))
!---
!-- Transfer the latitude and the latitude bounds
!---
    iret = NF90_PUT_VAR(nfid,nlatid,y_1d(1:pim))
!---
    IF (transp) THEN
      bounds_trans = TRANSPOSE(y_bnds)
    ELSE
      bounds_trans = y_bnds
    ENDIF
    iret = NF90_PUT_VAR(nfid,nlatidb,bounds_trans(1:nbbounds,1:pim))
!---
    DEALLOCATE(bounds_trans)
  ENDIF
!-
  iret = NF90_REDEF(nfid)
!-----------------------
END SUBROUTINE histh_all
!===
SUBROUTINE histvert (idf,pzaxname,pzaxtitle,pzaxunit, &
 &                   pzsize,pzvalues,pzaxid,pdirect)
!---------------------------------------------------------------------
!- This subroutine defines a vertical axis and returns it s id.
!- It gives the user the possibility to the user to define many
!- different vertical axes. For each variable defined with histdef a
!- vertical axis can be specified with by it s ID.
!-
!- INPUT
!-
!- idf      : ID of the file the variable should be archived in
!- pzaxname : Name of the vertical axis
!- pzaxtitle: title of the vertical axis
!- pzaxunit : Units of the vertical axis (no units if blank string)
!- pzsize   : size of the vertical axis
!- pzvalues : Coordinate values of the vetical axis
!-
!- pdirect  : is an optional argument which allows to specify the
!-            the positive direction of the axis : up or down.
!- OUTPUT
!-
!- pzaxid   : Returns the ID of the axis.
!-            Note that this is not the netCDF ID !
!-
!- VERSION
!-
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: idf,pzsize
  CHARACTER(LEN=*),INTENT(IN) :: pzaxname,pzaxunit,pzaxtitle
  REAL,INTENT(IN) :: pzvalues(pzsize)
  INTEGER,INTENT(OUT) :: pzaxid
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: pdirect
!-
  INTEGER :: pos,iv,zdimid,zaxid_tmp
  CHARACTER(LEN=70) :: str71
  CHARACTER(LEN=20) :: direction
  INTEGER :: iret,leng,nfid
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! 1.0 Verifications :
!    Do we have enough space for an extra axis ?
!    Is the name already in use ?
!-
  IF (l_dbg) WRITE(*,*) "histvert : 1.0 Verifications", &
 &                      pzaxname,'---',pzaxunit,'---',pzaxtitle
!-
! Direction of the vertical axis. Can we get if from the user.
!-
  IF (PRESENT(pdirect)) THEN
    direction = TRIM(pdirect)
    CALL strlowercase (direction)
  ELSE
    direction = 'unknown'
  ENDIF
!-
! Check the consistency of the attribute
!-
  IF (     PRESENT(pdirect)    &
 &    .AND.(direction /= 'up') &
 &    .AND.(direction /= 'down') ) THEN
    direction = 'unknown'
    CALL ipslerr (2,"histvert",&
 & "The specified positive direction for the vertical axis is invalid.",&
 & "The value must be up or down.","The attribute will not be written.")
  ENDIF
!-
  IF (W_F(idf)%n_zax+1 > nb_zax_max) THEN
    CALL ipslerr (3,"histvert", &
   &  'Table of vertical axes too small. You should increase ',&
   &  'nb_zax_max in histcom.f90 in order to accomodate all ', &
   &  'these variables ')
  ENDIF
!-
  iv = W_F(idf)%n_zax
  IF (iv > 1) THEN
    CALL find_str (W_F(idf)%zax_name(1:iv-1),pzaxname,pos)
  ELSE
    pos = 0
  ENDIF
!-
  IF (pos > 0) THEN
    WRITE(str71,'("Check variable ",A," in file",I3)') &
 &    TRIM(pzaxname),idf
    CALL ipslerr (3,"histvert", &
 &    "Vertical axis already exists",TRIM(str71), &
 &    "Can also be a wrong file ID in another declaration")
  ENDIF
!-
  iv = W_F(idf)%n_zax+1
!-
! 2.0 Add the information to the file
!-
  IF (l_dbg) &
 &  WRITE(*,*) "histvert : 2.0 Add the information to the file"
!-
  nfid = W_F(idf)%ncfid
!-
  leng = MIN(LEN_TRIM(pzaxname),20)
  iret = NF90_DEF_DIM (nfid,pzaxname(1:leng),pzsize,zaxid_tmp)
  iret = NF90_DEF_VAR (nfid,pzaxname(1:leng),NF90_REAL4, &
 &                     zaxid_tmp,zdimid)
  iret = NF90_PUT_ATT (nfid,zdimid,'axis',"Z")
  iret = NF90_PUT_ATT (nfid,zdimid,'standard_name',"model_level_number")
  leng = MIN(LEN_TRIM(pzaxunit),20)
  IF (leng > 0) THEN
    iret = NF90_PUT_ATT (nfid,zdimid,'units',pzaxunit(1:leng))
  ENDIF
  IF (direction /= 'unknown') THEN
    iret = NF90_PUT_ATT (nfid,zdimid,'positive',TRIM(direction))
  ENDIF
  iret = NF90_PUT_ATT (nfid,zdimid,'valid_min', &
 &                     REAL(MINVAL(pzvalues(1:pzsize)),KIND=4))
  iret = NF90_PUT_ATT (nfid,zdimid,'valid_max', &
 &                     REAL(MAXVAL(pzvalues(1:pzsize)),KIND=4))
  leng = MIN(LEN_TRIM(pzaxname),20)
  iret = NF90_PUT_ATT (nfid,zdimid,'title',pzaxname(1:leng))
  leng = MIN(LEN_TRIM(pzaxtitle),80)
  iret = NF90_PUT_ATT (nfid,zdimid,'long_name',pzaxtitle(1:leng))
!-
  iret = NF90_ENDDEF (nfid)
!-
  iret = NF90_PUT_VAR (nfid,zdimid,pzvalues(1:pzsize))
!-
  iret = NF90_REDEF (nfid)
!-
!- 3.0 add the information to the common
!-
  IF (l_dbg) &
  &  WRITE(*,*) "histvert : 3.0 add the information to the common"
!-
  W_F(idf)%n_zax = iv
  W_F(idf)%zax_size(iv) = pzsize
  W_F(idf)%zax_name(iv) = pzaxname
  W_F(idf)%zax_ids(iv) = zaxid_tmp
  pzaxid = iv
!----------------------
END SUBROUTINE histvert
!===
SUBROUTINE histdef &
 &  (idf,pvarname,ptitle,punit, &
 &   pxsize,pysize,phoriid,pzsize,par_oriz,par_szz,pzid, &
 &   xtype,popp,pfreq_opp,pfreq_wrt,var_range,standard_name)
!---------------------------------------------------------------------
!- With this subroutine each variable to be archived on the history
!- tape should be declared.
!-
!- It gives the user the choise of operation
!- to be performed on the variables, the frequency of this operation
!- and finaly the frequency of the archiving.
!-
!- INPUT
!-
!- idf      : ID of the file the variable should be archived in
!- pvarname : Name of the variable, short and easy to remember
!- ptitle   : Full name of the variable
!- punit    : Units of the variable (no units if blank string)
!-
!- The next 3 arguments give the size of that data
!- that will be passed to histwrite. The zoom will be
!- done there with the horizontal information obtained
!- in histbeg and the vertical information to follow.
!-
!- pxsize   : Size in X direction (size of the data that will be
!-            given to histwrite)
!- pysize   : Size in Y direction
!- phoriid  : ID of the horizontal axis
!-
!- The next two arguments give the vertical zoom to use.
!-
!- pzsize   : Size in Z direction (If 1 then no axis is declared
!-            for this variable and pzid is not used)
!- par_oriz : Off set of the zoom
!- par_szz  : Size of the zoom
!-
!- pzid     : ID of the vertical axis to use. It has to have
!-            the size of the zoom.
!- xtype    : External netCDF type (hist_r4/hist_r8)
!- popp     : Operation to be performed. The following options
!-            exist today :
!- inst : keeps instantaneous values for writting
!- ave  : Computes the average from call between writes
!- pfreq_opp: Frequency of this operation (in seconds)
!- pfreq_wrt: Frequency at which the variable should be
!-            written (in seconds)
!- var_range: Range of the variable.
!-            If the minimum is greater than the maximum,
!-            the values will be calculated.
!-
!- VERSION
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: idf,pxsize,pysize,pzsize,pzid
  INTEGER,INTENT(IN) :: par_oriz,par_szz,xtype,phoriid
  CHARACTER(LEN=*),INTENT(IN) :: pvarname,punit,popp,ptitle
  REAL,INTENT(IN) :: pfreq_opp,pfreq_wrt
  REAL,DIMENSION(2),OPTIONAL,INTENT(IN) :: var_range
  CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: standard_name
!-
  INTEGER :: iv
  CHARACTER(LEN=70) :: str70,str71,str72
  CHARACTER(LEN=20) :: tmp_name
  CHARACTER(LEN=40) :: str40
  CHARACTER(LEN=10) :: str10
  CHARACTER(LEN=120) :: ex_topps
  REAL :: un_an,un_jour,test_fopp,test_fwrt
  INTEGER :: pos,buff_sz
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  ex_topps = 'ave, inst, t_min, t_max, t_sum, once, never, l_max, l_min'
!-
  W_F(idf)%n_var = W_F(idf)%n_var+1
  iv = W_F(idf)%n_var
!-
  IF (iv > nb_var_max) THEN
    CALL ipslerr (3,"histdef", &
   &  'Table of variables too small. You should increase nb_var_max',&
   &  'in histcom.f90 in order to accomodate all these variables', &
   &  ' ')
  ENDIF
!-
! 1.0 Transfer informations on the variable to the common
!     and verify that it does not already exist
!-
  IF (l_dbg) WRITE(*,*) "histdef : 1.0"
!-
  IF (iv > 1) THEN
    CALL find_str (W_F(idf)%W_V(1:iv-1)%v_name,pvarname,pos)
  ELSE
    pos = 0
  ENDIF
!-
  IF (pos > 0) THEN
    str70 = "Variable already exists"
    WRITE(str71,'("Check variable  ",a," in file",I3)') &
 &    TRIM(pvarname),idf
    str72 = "Can also be a wrong file ID in another declaration"
    CALL ipslerr (3,"histdef",str70,str71,str72)
  ENDIF
!-
  W_F(idf)%W_V(iv)%v_name = pvarname
  W_F(idf)%W_V(iv)%title = ptitle
  W_F(idf)%W_V(iv)%unit_name = punit
  IF (PRESENT(standard_name)) THEN
    W_F(idf)%W_V(iv)%std_name = standard_name
  ELSE
    W_F(idf)%W_V(iv)%std_name = ptitle
  ENDIF
  tmp_name = W_F(idf)%W_V(iv)%v_name
!-
! 1.1 decode the operations
!-
  W_F(idf)%W_V(iv)%fullop = popp
  CALL buildop &
 &  (TRIM(popp),ex_topps,W_F(idf)%W_V(iv)%topp,missing_val, &
 &   W_F(idf)%W_V(iv)%sopp,W_F(idf)%W_V(iv)%scal, &
 &   W_F(idf)%W_V(iv)%nbopp)
!-
! 1.2 If we have an even number of operations
!     then we need to add identity
!-
  IF ( MOD(W_F(idf)%W_V(iv)%nbopp,2) == 0) THEN
    W_F(idf)%W_V(iv)%nbopp = W_F(idf)%W_V(iv)%nbopp+1
    W_F(idf)%W_V(iv)%sopp(W_F(idf)%W_V(iv)%nbopp) = 'ident'
    W_F(idf)%W_V(iv)%scal(W_F(idf)%W_V(iv)%nbopp) = missing_val
  ENDIF
!-
! 1.3 External type of the variable
!-
  IF (xtype == hist_r8) THEN
    W_F(idf)%W_V(iv)%v_typ = hist_r8
  ELSE
    W_F(idf)%W_V(iv)%v_typ = hist_r4
  ENDIF
!-
! 2.0 Put the size of the variable in the common and check
!-
  IF (l_dbg) THEN
    WRITE(*,*) "histdef : 2.0",idf,iv,W_F(idf)%W_V(iv)%nbopp, &
 &    W_F(idf)%W_V(iv)%sopp(1:W_F(idf)%W_V(iv)%nbopp), &
 &    W_F(idf)%W_V(iv)%scal(1:W_F(idf)%W_V(iv)%nbopp)
  ENDIF
!-
  W_F(idf)%W_V(iv)%scsize(1:3) = (/ pxsize,pysize,pzsize /)
  W_F(idf)%W_V(iv)%zorig(1:3) = &
 &  (/ W_F(idf)%slab_ori(1),W_F(idf)%slab_ori(2),par_oriz /)
  W_F(idf)%W_V(iv)%zsize(1:3) = &
 &  (/ W_F(idf)%slab_siz(1),W_F(idf)%slab_siz(2),par_szz /)
!-
! Is the size of the full array the same as that of the coordinates ?
!-
  IF (    (pxsize > W_F(idf)%full_size(1)) &
 &    .OR.(pysize > W_F(idf)%full_size(2)) ) THEN
!-
    str70 = "The size of the variable is different "// &
 &          "from the one of the coordinates"
    WRITE(str71,'("Size of coordinates :",2I4)') &
 &   W_F(idf)%full_size(1),W_F(idf)%full_size(2)
    WRITE(str72,'("Size declared for variable ",a," :",2I4)') &
 &   TRIM(tmp_name),pxsize,pysize
    CALL ipslerr (3,"histdef",str70,str71,str72)
  ENDIF
!-
! Is the size of the zoom smaller than the coordinates ?
!-
  IF (    (W_F(idf)%full_size(1) < W_F(idf)%slab_siz(1)) &
 &    .OR.(W_F(idf)%full_size(2) < W_F(idf)%slab_siz(2)) ) THEN
    str70 = &
 &   "Size of variable should be greater or equal to those of the zoom"
    WRITE(str71,'("Size of XY zoom :",2I4)') &
 &   W_F(idf)%slab_siz(1),W_F(idf)%slab_siz(2)
    WRITE(str72,'("Size declared for variable ",A," :",2I4)') &
 &   TRIM(tmp_name),pxsize,pysize
    CALL ipslerr (3,"histdef",str70,str71,str72)
  ENDIF
!-
! 2.1 We store the horizontal grid information with minimal
!     and a fall back onto the default grid
!-
  IF ( (phoriid > 0).AND.(phoriid <= W_F(idf)%n_hax) ) THEN
    W_F(idf)%W_V(iv)%h_axid = phoriid
  ELSE
    W_F(idf)%W_V(iv)%h_axid = 1
    CALL ipslerr (2,"histdef", &
   &  'We use the default grid for variable as an invalide',&
   &  'ID was provided for variable : ',TRIM(pvarname))
  ENDIF
!-
! 2.2 Check the vertical coordinates if needed
!-
  IF (par_szz > 1) THEN
!-
!-- Does the vertical coordinate exist ?
!-
    IF (pzid > W_F(idf)%n_zax) THEN
      WRITE(str70, &
 &    '("The vertical coordinate chosen for variable ",A)') &
 &     TRIM(tmp_name)
      str71 = " Does not exist."
      CALL ipslerr (3,"histdef",str70,str71," ")
    ENDIF
!-
!-- Is the vertical size of the variable equal to that of the axis ?
!-
    IF (par_szz /= W_F(idf)%zax_size(pzid)) THEN
      str70 = "The size of the zoom does not correspond "// &
 &            "to the size of the chosen vertical axis"
      WRITE(str71,'("Size of zoom in z :",I4)') par_szz
      WRITE(str72,'("Size declared for axis ",A," :",I4)') &
 &     TRIM(W_F(idf)%zax_name(pzid)),W_F(idf)%zax_size(pzid)
      CALL ipslerr (3,"histdef",str70,str71,str72)
    ENDIF
!-
!-- Is the zoom smaller that the total size of the variable ?
!-
    IF (pzsize < par_szz) THEN
      str70 = "The vertical size of variable "// &
 &            "is smaller than that of the zoom."
      WRITE(str71,'("Declared vertical size of data :",I5)') pzsize
      WRITE(str72,'("Size of zoom for variable ",a," = ",I5)') &
 &     TRIM(tmp_name),par_szz
      CALL ipslerr (3,"histdef",str70,str71,str72)
    ENDIF
    W_F(idf)%W_V(iv)%z_axid = pzid
  ELSE
    W_F(idf)%W_V(iv)%z_axid = -99
  ENDIF
!-
! 3.0 We get the size of the arrays histwrite will get
!     and eventually allocate the time_buffer
!-
  IF (l_dbg) THEN
    WRITE(*,*) "histdef : 3.0"
  ENDIF
!-
  buff_sz = W_F(idf)%W_V(iv)%zsize(1) &
 &         *W_F(idf)%W_V(iv)%zsize(2) &
 &         *W_F(idf)%W_V(iv)%zsize(3)
!-
  IF (     (TRIM(W_F(idf)%W_V(iv)%topp) /= "inst") &
 &    .AND.(TRIM(W_F(idf)%W_V(iv)%topp) /= "once") &
 &    .AND.(TRIM(W_F(idf)%W_V(iv)%topp) /= "never") )THEN
    ALLOCATE(W_F(idf)%W_V(iv)%t_bf(buff_sz))
    W_F(idf)%W_V(iv)%t_bf(:) = 0.
    IF (l_dbg) THEN
      WRITE(*,*) "histdef : 3.0 allocating time_buffer for", &
 &      " idf = ",idf," iv = ",iv," size = ",buff_sz
    ENDIF
  ENDIF
!-
! 4.0 Transfer the frequency of the operations and check
!     for validity. We have to pay attention to negative values
!     of the frequency which indicate monthly time-steps.
!     The strategy is to bring it back to seconds for the tests
!-
  IF (l_dbg) WRITE(*,*) "histdef : 4.0"
!-
  W_F(idf)%W_V(iv)%freq_opp = pfreq_opp
  W_F(idf)%W_V(iv)%freq_wrt = pfreq_wrt
!-
  CALL ioget_calendar(un_an,un_jour)
  IF (pfreq_opp < 0) THEN
    CALL ioget_calendar(un_an)
    test_fopp = pfreq_opp*(-1.)*un_an/12.*un_jour
  ELSE
    test_fopp = pfreq_opp
  ENDIF
  IF (pfreq_wrt < 0) THEN
    CALL ioget_calendar(un_an)
    test_fwrt = pfreq_wrt*(-1.)*un_an/12.*un_jour
  ELSE
    test_fwrt = pfreq_wrt
  ENDIF
!-
! 4.1 Frequency of operations and output should be larger than deltat !
!-
  IF (test_fopp < W_F(idf)%deltat) THEN
    str70 = 'Frequency of operations should be larger than deltat'
    WRITE(str71,'("It is not the case for variable ",a," :",F10.4)') &
 &   TRIM(tmp_name),pfreq_opp
    str72 = "PATCH : frequency set to deltat"
!-
    CALL ipslerr (2,"histdef",str70,str71,str72)
!-
    W_F(idf)%W_V(iv)%freq_opp = W_F(idf)%deltat
  ENDIF
!-
  IF (test_fwrt < W_F(idf)%deltat) THEN
    str70 = 'Frequency of output should be larger than deltat'
    WRITE(str71,'("It is not the case for variable ",a," :",F10.4)') &
 &   TRIM(tmp_name),pfreq_wrt
    str72 = "PATCH : frequency set to deltat"
!-
    CALL ipslerr (2,"histdef",str70,str71,str72)
!-
    W_F(idf)%W_V(iv)%freq_wrt = W_F(idf)%deltat
  ENDIF
!-
! 4.2 First the existence of the operation is tested and then
!     its compaticility with the choice of frequencies
!-
  IF (TRIM(W_F(idf)%W_V(iv)%topp) == "inst") THEN
    IF (test_fopp /= test_fwrt) THEN
      str70 = 'For instantaneous output the frequency '// &
 &            'of operations and output'
      WRITE(str71, &
 &     '("should be the same, this was not case for variable ",a)') &
 &      TRIM(tmp_name)
      str72 = "PATCH : The smalest frequency of both is used"
      CALL ipslerr (2,"histdef",str70,str71,str72)
      IF (test_fopp < test_fwrt) THEN
        W_F(idf)%W_V(iv)%freq_opp = pfreq_opp
        W_F(idf)%W_V(iv)%freq_wrt = pfreq_opp
      ELSE
        W_F(idf)%W_V(iv)%freq_opp = pfreq_wrt
        W_F(idf)%W_V(iv)%freq_wrt = pfreq_wrt
      ENDIF
    ENDIF
  ELSE IF (INDEX(ex_topps,TRIM(W_F(idf)%W_V(iv)%topp)) > 0) THEN
    IF (test_fopp > test_fwrt) THEN
      str70 = 'For averages the frequency of operations '// &
 &            'should be smaller or equal'
      WRITE(str71, &
 &     '("to that of output. It is not the case for variable ",a)') &
 &     TRIM(tmp_name)
      str72 = 'PATCH : The output frequency is used for both'
      CALL ipslerr (2,"histdef",str70,str71,str72)
      W_F(idf)%W_V(iv)%freq_opp = pfreq_wrt
    ENDIF
  ELSE
    WRITE (str70,'("Operation on variable ",A," is unknown")') &
 &   TRIM(tmp_name)
    WRITE (str71,'("operation requested is :",A)') &
 &   W_F(idf)%W_V(iv)%topp
    WRITE (str72,'("File ID :",I3)') idf
    CALL ipslerr (3,"histdef",str70,str71,str72)
  ENDIF
!-
! 5.0 Initialize other variables of the common
!-
  IF (l_dbg) WRITE(*,*) "histdef : 5.0"
!-
  W_F(idf)%W_V(iv)%hist_wrt_rng = (PRESENT(var_range))
  IF (W_F(idf)%W_V(iv)%hist_wrt_rng) THEN
    W_F(idf)%W_V(iv)%hist_calc_rng = (var_range(1) > var_range(2))
    IF (W_F(idf)%W_V(iv)%hist_calc_rng) THEN
      W_F(idf)%W_V(iv)%hist_minmax(1:2) = &
 &      (/ ABS(missing_val),-ABS(missing_val) /)
    ELSE
      W_F(idf)%W_V(iv)%hist_minmax(1:2) = var_range(1:2)
    ENDIF
  ENDIF
!-
! - freq_opp(idf,iv)/2./deltat(idf)
  W_F(idf)%W_V(iv)%last_opp = W_F(idf)%itau0
! - freq_wrt(idf,iv)/2./deltat(idf)
  W_F(idf)%W_V(iv)%last_wrt = W_F(idf)%itau0
! - freq_opp(idf,iv)/2./deltat(idf)
  W_F(idf)%W_V(iv)%last_opp_chk = W_F(idf)%itau0
! - freq_wrt(idf,iv)/2./deltat(idf)
  W_F(idf)%W_V(iv)%last_wrt_chk = W_F(idf)%itau0
  W_F(idf)%W_V(iv)%nb_opp = 0
  W_F(idf)%W_V(iv)%nb_wrt = 0
!-
! 6.0 Get the time axis for this variable
!-
  IF (l_dbg) WRITE(*,*) "histdef : 6.0"
!-
! No time axis for once, l_max, l_min or never operation
!-
  IF (     (TRIM(W_F(idf)%W_V(iv)%topp) /= 'once')  &
 &    .AND.(TRIM(W_F(idf)%W_V(iv)%topp) /= 'never') &
 &    .AND.(TRIM(W_F(idf)%W_V(iv)%topp) /= 'l_max') &
 &    .AND.(TRIM(W_F(idf)%W_V(iv)%topp) /= 'l_min') ) THEN
    IF (TRIM(W_F(idf)%W_V(iv)%topp) == 'inst') THEN
      str10 = 't_inst_'
    ELSE
      str10 = 't_op_'
    ENDIF
    IF (W_F(idf)%W_V(iv)%freq_wrt > 0) THEN
      WRITE (UNIT=str40,FMT='(A,I8.8)') &
&      TRIM(str10),INT(W_F(idf)%W_V(iv)%freq_wrt)
    ELSE
      WRITE (UNIT=str40,FMT='(A,I2.2,"month")') &
&      TRIM(str10),ABS(INT(W_F(idf)%W_V(iv)%freq_wrt))
    ENDIF
    CALL find_str (W_F(idf)%W_V(1:W_F(idf)%n_tax)%tax_name,str40,pos)
    IF (pos < 0) THEN
      W_F(idf)%n_tax = W_F(idf)%n_tax+1
      W_F(idf)%W_V(iv)%l_bnd = &
 &      (TRIM(W_F(idf)%W_V(iv)%topp) /= 'inst')
      W_F(idf)%W_V(W_F(idf)%n_tax)%tax_name = str40
      W_F(idf)%W_V(W_F(idf)%n_tax)%tax_last = 0
      W_F(idf)%W_V(iv)%t_axid = W_F(idf)%n_tax
    ELSE
      W_F(idf)%W_V(iv)%t_axid = pos
    ENDIF
  ELSE
    IF (l_dbg) THEN
      WRITE(*,*) "histdef : 7.0 ",TRIM(W_F(idf)%W_V(iv)%topp),'----'
    ENDIF
    W_F(idf)%W_V(iv)%t_axid = -99
  ENDIF
!-
! 7.0 prepare frequence of writing and operation
!     for never or once operation
!-
  IF (    (TRIM(W_F(idf)%W_V(iv)%topp) == 'once')  &
 &    .OR.(TRIM(W_F(idf)%W_V(iv)%topp) == 'never') ) THEN
    W_F(idf)%W_V(iv)%freq_opp = 0.
    W_F(idf)%W_V(iv)%freq_wrt = 0.
  ENDIF
!---------------------
END SUBROUTINE histdef
!===
SUBROUTINE histend (idf, snc4chunks)
!---------------------------------------------------------------------
!- This subroutine end the decalaration of variables and sets the
!- time axes in the netcdf file and puts it into the write mode.
!-
!- INPUT
!-
!- idf : ID of the file to be worked on
!-
!- VERSION
!-
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: idf
  TYPE(snc4_ctl), OPTIONAL,INTENT(IN) :: snc4chunks
!-
  INTEGER :: nfid,nvid,iret,ndim,iv,itx,ziv,itax,dim_cnt
  INTEGER,DIMENSION(4) :: dims
  INTEGER :: year,month,day,hours,minutes
  REAL :: sec
  REAL :: rtime0
  CHARACTER(LEN=30) :: str30
  CHARACTER(LEN=35) :: str35
  CHARACTER(LEN=120) :: assoc
  CHARACTER(LEN=70) :: str70
  CHARACTER(LEN=3),DIMENSION(12) :: cal =   &
 &  (/ 'JAN','FEB','MAR','APR','MAY','JUN', &
 &     'JUL','AUG','SEP','OCT','NOV','DEC' /)
  CHARACTER(LEN=7) :: tmp_opp
  LOGICAL :: l_b
  LOGICAL :: l_dbg
  INTEGER, DIMENSION(4) :: ichunksz    ! NETCDF4 chunk sizes
  INTEGER :: ichunkalg, ishuffle,&
             ideflate, ideflate_level
  LOGICAL :: lchunk = .FALSE.   ! logical switch to activate chunking when appropriate
!-
  ! NetCDF4 chunking and compression parameters
  ichunkalg = 0
  ishuffle = 1
  ideflate = 1
  ideflate_level = 1
  !
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  nfid = W_F(idf)%ncfid
!-
! 1.0 Create the time axes
!-
  IF (l_dbg) WRITE(*,*) "histend : 1.0"
!-
! 1.1 Define the time dimensions needed for this file
!-
  iret = NF90_DEF_DIM (nfid,'time_counter', &
 &                     NF90_UNLIMITED,W_F(idf)%tid)
  DO iv=1,W_F(idf)%n_var
    IF (W_F(idf)%W_V(iv)%l_bnd) THEN
      iret = NF90_DEF_DIM (nfid,'tbnds',2,W_F(idf)%bid)
      EXIT
    ENDIF
  ENDDO
!-
! 1.2 Define all the time axes needed for this file
!-
  DO itx=1,W_F(idf)%n_tax
    dims(1) = W_F(idf)%tid
    l_b = (INDEX(W_F(idf)%W_V(itx)%tax_name,"t_op_") == 1)
    IF (itx > 1) THEN
      str30 = W_F(idf)%W_V(itx)%tax_name
    ELSE
      str30 = "time_counter"
    ENDIF
    IF (l_b) THEN
      str35 = TRIM(str30)//'_bnds'
    ENDIF
    iret = NF90_DEF_VAR (nfid,TRIM(str30),NF90_REAL8, &
 &          dims(1),W_F(idf)%W_V(itx)%tdimid)
    IF (itx <= 1) THEN
      iret = NF90_PUT_ATT (nfid,W_F(idf)%W_V(itx)%tdimid,'axis',"T")
    ENDIF
    iret = NF90_PUT_ATT (nfid,W_F(idf)%W_V(itx)%tdimid, &
 &          'standard_name',"time")
!---
!   To transform the current itau into a real date and take it
!   as the origin of the file requires the time counter to change.
!   Thus it is an operation the user has to ask for.
!   This function should thus only be re-instated
!   if there is a ioconf routine to control it.
!---
!-- rtime0 = itau2date(itau0(idf),date0(idf),deltat(idf))
    rtime0 = W_F(idf)%date0
!-
    CALL ju2ymds(rtime0,year,month,day,sec)
!---
!   Catch any error induced by a change in calendar !
!---
    IF (year < 0) THEN
      year = 2000+year
    ENDIF
!-
    hours = INT(sec/(60.*60.))
    minutes = INT((sec-hours*60.*60.)/60.)
    sec = sec-(hours*60.*60.+minutes*60.)
!-
    WRITE (UNIT=str70, &
 &   FMT='(A,I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
 &    'seconds since ',year,month,day,hours,minutes,INT(sec)
    iret = NF90_PUT_ATT (nfid,W_F(idf)%W_V(itx)%tdimid, &
 &           'units',TRIM(str70))
!-
    CALL ioget_calendar (str30)
    iret = NF90_PUT_ATT (nfid,W_F(idf)%W_V(itx)%tdimid, &
 &           'calendar',TRIM(str30))
!-
    iret = NF90_PUT_ATT (nfid,W_F(idf)%W_V(itx)%tdimid, &
 &           'title','Time')
!-
    iret = NF90_PUT_ATT (nfid,W_F(idf)%W_V(itx)%tdimid, &
 &           'long_name','Time axis')
!-
    WRITE (UNIT=str70, &
 &   FMT='(" ",I4.4,"-",A3,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
 &    year,cal(month),day,hours,minutes,INT(sec)
    iret = NF90_PUT_ATT (nfid,W_F(idf)%W_V(itx)%tdimid, &
 &           'time_origin',TRIM(str70))
!---
    IF (l_b) THEN
      iret = NF90_PUT_ATT (nfid,W_F(idf)%W_V(itx)%tdimid, &
 &             'bounds',TRIM(str35))
      dims(1:2) = (/ W_F(idf)%bid,W_F(idf)%tid /)
      iret = NF90_DEF_VAR (nfid,TRIM(str35),NF90_REAL8, &
 &            dims(1:2),W_F(idf)%W_V(itx)%tbndid)
    ENDIF
  ENDDO
!-
! 2.0 declare the variables
!-
  IF (l_dbg) WRITE(*,*) "histend : 2.0"
!-
  DO iv=1,W_F(idf)%n_var
!---
    itax = W_F(idf)%W_V(iv)%t_axid
!---
    IF (W_F(idf)%regular) THEN
      dims(1:2) = (/ W_F(idf)%xid,W_F(idf)%yid /)
      dim_cnt = 2
    ELSE
      dims(1) = W_F(idf)%xid
      dim_cnt = 1
    ENDIF
!---
    tmp_opp = W_F(idf)%W_V(iv)%topp
    ziv = W_F(idf)%W_V(iv)%z_axid
!---
!   2.1 dimension of field
!---
    IF ((TRIM(tmp_opp) /= 'never')) THEN
      IF (     (TRIM(tmp_opp) /= 'once')  &
     &    .AND.(TRIM(tmp_opp) /= 'l_max') &
     &    .AND.(TRIM(tmp_opp) /= 'l_min') ) THEN
        IF (ziv == -99) THEN
          ndim = dim_cnt+1
          dims(dim_cnt+1:dim_cnt+2) = (/ W_F(idf)%tid,0 /)
        ELSE
          ndim = dim_cnt+2
          dims(dim_cnt+1:dim_cnt+2) = &
 &         (/ W_F(idf)%zax_ids(ziv),W_F(idf)%tid /)
        ENDIF
      ELSE
        IF (ziv == -99) THEN
          ndim = dim_cnt
          dims(dim_cnt+1:dim_cnt+2) = (/ 0,0 /)
        ELSE
          ndim = dim_cnt+1
          dims(dim_cnt+1:dim_cnt+2) = (/ W_F(idf)%zax_ids(ziv),0 /)
        ENDIF
      ENDIF
!-
      iret = NF90_DEF_VAR (nfid,TRIM(W_F(idf)%W_V(iv)%v_name), &
 &             W_F(idf)%W_V(iv)%v_typ,dims(1:ABS(ndim)),nvid)
!-
      IF( ndim == 4 ) THEN
        IF( PRESENT( snc4chunks ) ) THEN
          IF( snc4chunks%luse ) THEN
           ichunksz = 1
           iret = NF90_INQUIRE_DIMENSION( nfid, W_F(idf)%xid, len = ichunksz(1) )
           iret = NF90_INQUIRE_DIMENSION( nfid, W_F(idf)%yid, len = ichunksz(2) )
           IF ( ziv .NE. -99 ) &
              iret = NF90_INQUIRE_DIMENSION( nfid, W_F(idf)%zax_ids(ziv), len = ichunksz(3) )
           ichunksz(1) = MIN(ichunksz(1), MAX((ichunksz(1)-1)/snc4chunks%ni + 1,16))
           ichunksz(2) = MIN(ichunksz(2), MAX((ichunksz(2)-1)/snc4chunks%nj + 1,16))
           ichunksz(3) = MIN(ichunksz(3), MAX((ichunksz(3)-1)/snc4chunks%nk + 1, 1))
           ! Always use a chunk size of 1 for the unlimited dimension
           iret = SET_NF90_DEF_VAR_CHUNKING(nfid, nvid, ichunkalg, ichunksz)
           iret = SET_NF90_DEF_VAR_DEFLATE(nfid, nvid, ishuffle, ideflate, ideflate_level)
          ENDIF
        ENDIF
      ENDIF
      W_F(idf)%W_V(iv)%ncvid = nvid
!-
      IF (LEN_TRIM(W_F(idf)%W_V(iv)%unit_name) > 0) THEN
        iret = NF90_PUT_ATT (nfid,nvid,'units', &
 &                           TRIM(W_F(idf)%W_V(iv)%unit_name))
      ENDIF
      iret = NF90_PUT_ATT (nfid,nvid,'standard_name', &
 &                         TRIM(W_F(idf)%W_V(iv)%std_name))
!-
      IF (W_F(idf)%W_V(iv)%v_typ == hist_r8) THEN
        iret = NF90_PUT_ATT (nfid,nvid,'_FillValue',NF90_FILL_REAL8)
      ELSE
        iret = NF90_PUT_ATT (nfid,nvid,'_FillValue',NF90_FILL_REAL4)
      ENDIF
      IF (W_F(idf)%W_V(iv)%hist_wrt_rng) THEN
        IF (W_F(idf)%W_V(iv)%v_typ == hist_r8) THEN
          iret = NF90_PUT_ATT (nfid,nvid,'valid_min', &
 &                 REAL(W_F(idf)%W_V(iv)%hist_minmax(1),KIND=8))
          iret = NF90_PUT_ATT (nfid,nvid,'valid_max', &
 &                 REAL(W_F(idf)%W_V(iv)%hist_minmax(2),KIND=8))
        ELSE
          iret = NF90_PUT_ATT (nfid,nvid,'valid_min', &
 &                 REAL(W_F(idf)%W_V(iv)%hist_minmax(1),KIND=4))
          iret = NF90_PUT_ATT (nfid,nvid,'valid_max', &
 &                 REAL(W_F(idf)%W_V(iv)%hist_minmax(2),KIND=4))
        ENDIF
      ENDIF
      iret = NF90_PUT_ATT (nfid,nvid,'long_name', &
 &                         TRIM(W_F(idf)%W_V(iv)%title))
      iret = NF90_PUT_ATT (nfid,nvid,'online_operation', &
 &                         TRIM(W_F(idf)%W_V(iv)%fullop))
!-
      SELECT CASE(ndim)
      CASE(-3,2:4)
      CASE DEFAULT
        CALL ipslerr (3,"histend", &
       &  'less than 2 or more than 4 dimensions are not', &
       &  'allowed at this stage',' ')
      END SELECT
!-
      assoc=TRIM(W_F(idf)%hax_name(W_F(idf)%W_V(iv)%h_axid,2)) &
 &   //' '//TRIM(W_F(idf)%hax_name(W_F(idf)%W_V(iv)%h_axid,1))
!-
      ziv = W_F(idf)%W_V(iv)%z_axid
      IF (ziv > 0) THEN
        str30 = W_F(idf)%zax_name(ziv)
        assoc = TRIM(str30)//' '//TRIM(assoc)
      ENDIF
!-
      IF (itax > 0) THEN
        IF (itax > 1) THEN
          str30 = W_F(idf)%W_V(itax)%tax_name
        ELSE
          str30 = "time_counter"
        ENDIF
        assoc = TRIM(str30)//' '//TRIM(assoc)
!-
        IF (l_dbg) THEN
          WRITE(*,*) "histend : 2.0.n, freq_opp, freq_wrt", &
 &          W_F(idf)%W_V(iv)%freq_opp,W_F(idf)%W_V(iv)%freq_wrt
        ENDIF
!-
        iret = NF90_PUT_ATT (nfid,nvid,'interval_operation', &
 &                           REAL(W_F(idf)%W_V(iv)%freq_opp,KIND=4))
        iret = NF90_PUT_ATT (nfid,nvid,'interval_write', &
 &                           REAL(W_F(idf)%W_V(iv)%freq_wrt,KIND=4))
      ENDIF
      iret = NF90_PUT_ATT (nfid,nvid,'coordinates',TRIM(assoc))
    ENDIF
  ENDDO
!-
! 2.2 Add DOMAIN attributes if needed
!-
  IF (W_F(idf)%dom_id_svg >= 0) THEN
    CALL flio_dom_att (nfid,W_F(idf)%dom_id_svg)
  ENDIF
!-
! 3.0 Put the netcdf file into write mode
!-
  IF (l_dbg) WRITE(*,*) "histend : 3.0"
!-
  iret = NF90_ENDDEF (nfid)
!-
! 4.0 Give some informations to the user
!-
  IF (l_dbg) WRITE(*,*) "histend : 4.0"
!-
!!$  WRITE(str70,'("All variables have been initialized on file :",I3)') idf
!!$  CALL ipslerr (1,'histend',str70,'',' ')
!---------------------
END SUBROUTINE histend
!===
SUBROUTINE histwrite_r1d (idf,pvarname,pitau,pdata,nbindex,nindex)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: idf,pitau,nbindex
  INTEGER,DIMENSION(nbindex),INTENT(IN) :: nindex
  REAL,DIMENSION(:),INTENT(IN) :: pdata
  CHARACTER(LEN=*),INTENT(IN) :: pvarname
!---------------------------------------------------------------------
  CALL histw_rnd (idf,pvarname,pitau,nbindex,nindex,pdata_1d=pdata)
!---------------------------
END SUBROUTINE histwrite_r1d
!===
SUBROUTINE histwrite_r2d (idf,pvarname,pitau,pdata,nbindex,nindex)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: idf,pitau,nbindex
  INTEGER,DIMENSION(nbindex),INTENT(IN) :: nindex
  REAL,DIMENSION(:,:),INTENT(IN) :: pdata
  CHARACTER(LEN=*),INTENT(IN) :: pvarname
!---------------------------------------------------------------------
  CALL histw_rnd (idf,pvarname,pitau,nbindex,nindex,pdata_2d=pdata)
!---------------------------
END SUBROUTINE histwrite_r2d
!===
SUBROUTINE histwrite_r3d (idf,pvarname,pitau,pdata,nbindex,nindex)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: idf,pitau,nbindex
  INTEGER,DIMENSION(nbindex),INTENT(IN) :: nindex
  REAL,DIMENSION(:,:,:),INTENT(IN) :: pdata
  CHARACTER(LEN=*),INTENT(IN) :: pvarname
!---------------------------------------------------------------------
  CALL histw_rnd (idf,pvarname,pitau,nbindex,nindex,pdata_3d=pdata)
!---------------------------
END SUBROUTINE histwrite_r3d
!===
SUBROUTINE histw_rnd (idf,pvarname,pitau,nbindex,nindex, &
  &                   pdata_1d,pdata_2d,pdata_3d)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: idf,pitau,nbindex
  INTEGER,DIMENSION(nbindex),INTENT(IN) :: nindex
  CHARACTER(LEN=*),INTENT(IN) :: pvarname
  REAL,DIMENSION(:),INTENT(IN),OPTIONAL     :: pdata_1d
  REAL,DIMENSION(:,:),INTENT(IN),OPTIONAL   :: pdata_2d
  REAL,DIMENSION(:,:,:),INTENT(IN),OPTIONAL :: pdata_3d
!-
  LOGICAL :: do_oper,do_write,largebuf,l1d,l2d,l3d
  INTEGER :: iv,io,nbpt_out
  INTEGER              :: nbpt_in1
  INTEGER,DIMENSION(2) :: nbpt_in2
  INTEGER,DIMENSION(3) :: nbpt_in3
  REAL,ALLOCATABLE,DIMENSION(:),SAVE :: tbf_1
  CHARACTER(LEN=7) :: tmp_opp
  CHARACTER(LEN=13) :: c_nam
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  l1d=PRESENT(pdata_1d); l2d=PRESENT(pdata_2d); l3d=PRESENT(pdata_3d);
  IF      (l1d) THEN
    c_nam = 'histwrite_r1d'
  ELSE IF (l2d) THEN
    c_nam = 'histwrite_r2d'
  ELSE IF (l3d) THEN
    c_nam = 'histwrite_r3d'
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "histwrite : ",c_nam
  ENDIF
!-
! 1.0 Try to catch errors like specifying the wrong file ID.
!     Thanks Marine for showing us what errors users can make !
!-
  IF ( (idf < 1).OR.(idf > nb_files_max) ) THEN
    CALL ipslerr (3,"histwrite", &
 &    'Illegal file ID in the histwrite of variable',pvarname,' ')
  ENDIF
!-
! 1.1 Find the id of the variable to be written and the real time
!-
  CALL histvar_seq (idf,pvarname,iv)
!-
! 2.0 do nothing for never operation
!-
  tmp_opp = W_F(idf)%W_V(iv)%topp
!-
  IF (TRIM(tmp_opp) == "never") THEN
    W_F(idf)%W_V(iv)%last_opp_chk = -99
    W_F(idf)%W_V(iv)%last_wrt_chk = -99
  ENDIF
!-
! 3.0 We check if we need to do an operation
!-
  IF (W_F(idf)%W_V(iv)%last_opp_chk == pitau) THEN
    CALL ipslerr (3,"histwrite", &
 &    'This variable has already been analysed at the present', &
 &    'time step',TRIM(pvarname))
  ENDIF
!-
  CALL isittime &
 &  (pitau,W_F(idf)%date0,W_F(idf)%deltat, &
 &   W_F(idf)%W_V(iv)%freq_opp, &
 &   W_F(idf)%W_V(iv)%last_opp, &
 &   W_F(idf)%W_V(iv)%last_opp_chk,do_oper)
!-
! 4.0 We check if we need to write the data
!-
  IF (W_F(idf)%W_V(iv)%last_wrt_chk == pitau) THEN
    CALL ipslerr (3,"histwrite", &
 &    'This variable as already been written for the present', &
 &    'time step',' ')
  ENDIF
!-
  CALL isittime &
 &  (pitau,W_F(idf)%date0,W_F(idf)%deltat, &
 &   W_F(idf)%W_V(iv)%freq_wrt, &
 &   W_F(idf)%W_V(iv)%last_wrt, &
 &   W_F(idf)%W_V(iv)%last_wrt_chk,do_write)
!-
! 5.0 histwrite called
!-
  IF (do_oper.OR.do_write) THEN
!-
!-- 5.1 Get the sizes of the data we will handle
!-
    IF (W_F(idf)%W_V(iv)%datasz_in(1) <= 0) THEN
!---- There is the risk here that the user has over-sized the array.
!---- But how can we catch this ?
!---- In the worst case we will do impossible operations
!---- on part of the data !
      W_F(idf)%W_V(iv)%datasz_in(1:3) = -1
      IF      (l1d) THEN
        W_F(idf)%W_V(iv)%datasz_in(1) = SIZE(pdata_1d)
      ELSE IF (l2d) THEN
        W_F(idf)%W_V(iv)%datasz_in(1) = SIZE(pdata_2d,DIM=1)
        W_F(idf)%W_V(iv)%datasz_in(2) = SIZE(pdata_2d,DIM=2)
      ELSE IF (l3d) THEN
        W_F(idf)%W_V(iv)%datasz_in(1) = SIZE(pdata_3d,DIM=1)
        W_F(idf)%W_V(iv)%datasz_in(2) = SIZE(pdata_3d,DIM=2)
        W_F(idf)%W_V(iv)%datasz_in(3) = SIZE(pdata_3d,DIM=3)
      ENDIF
    ENDIF
!-
!-- 5.2 The maximum size of the data will give the size of the buffer
!-
    IF (W_F(idf)%W_V(iv)%datasz_max <= 0) THEN
      largebuf = .FALSE.
      DO io=1,W_F(idf)%W_V(iv)%nbopp
        IF (INDEX(fuchnbout,W_F(idf)%W_V(iv)%sopp(io)) > 0) THEN
          largebuf = .TRUE.
        ENDIF
      ENDDO
      IF (largebuf) THEN
        W_F(idf)%W_V(iv)%datasz_max = &
 &        W_F(idf)%W_V(iv)%scsize(1) &
 &       *W_F(idf)%W_V(iv)%scsize(2) &
 &       *W_F(idf)%W_V(iv)%scsize(3)
      ELSE
        IF      (l1d) THEN
          W_F(idf)%W_V(iv)%datasz_max = &
 &          W_F(idf)%W_V(iv)%datasz_in(1)
        ELSE IF (l2d) THEN
          W_F(idf)%W_V(iv)%datasz_max = &
 &          W_F(idf)%W_V(iv)%datasz_in(1) &
 &         *W_F(idf)%W_V(iv)%datasz_in(2)
        ELSE IF (l3d) THEN
          W_F(idf)%W_V(iv)%datasz_max = &
 &          W_F(idf)%W_V(iv)%datasz_in(1) &
 &         *W_F(idf)%W_V(iv)%datasz_in(2) &
 &         *W_F(idf)%W_V(iv)%datasz_in(3)
        ENDIF
      ENDIF
    ENDIF
!-
    IF (.NOT.ALLOCATED(tbf_1)) THEN
      IF (l_dbg) THEN
        WRITE(*,*) &
 &       c_nam//" : allocate tbf_1 for size = ", &
 &       W_F(idf)%W_V(iv)%datasz_max
      ENDIF
      ALLOCATE(tbf_1(W_F(idf)%W_V(iv)%datasz_max))
    ELSE IF (W_F(idf)%W_V(iv)%datasz_max > SIZE(tbf_1)) THEN
      IF (l_dbg) THEN
        WRITE(*,*) &
 &       c_nam//" : re-allocate tbf_1 for size = ", &
 &       W_F(idf)%W_V(iv)%datasz_max
      ENDIF
      DEALLOCATE(tbf_1)
      ALLOCATE(tbf_1(W_F(idf)%W_V(iv)%datasz_max))
    ENDIF
!-
!-- We have to do the first operation anyway.
!-- Thus we do it here and change the ranke
!-- of the data at the same time. This should speed up things.
!-
    nbpt_out = W_F(idf)%W_V(iv)%datasz_max
    IF      (l1d) THEN
      nbpt_in1 = W_F(idf)%W_V(iv)%datasz_in(1)
      CALL mathop (W_F(idf)%W_V(iv)%sopp(1),nbpt_in1,pdata_1d, &
 &                 missing_val,nbindex,nindex, &
 &                 W_F(idf)%W_V(iv)%scal(1),nbpt_out,tbf_1)
    ELSE IF (l2d) THEN
      nbpt_in2(1:2) = W_F(idf)%W_V(iv)%datasz_in(1:2)
      CALL mathop (W_F(idf)%W_V(iv)%sopp(1),nbpt_in2,pdata_2d, &
 &                 missing_val,nbindex,nindex, &
 &                 W_F(idf)%W_V(iv)%scal(1),nbpt_out,tbf_1)
    ELSE IF (l3d) THEN
      nbpt_in3(1:3) = W_F(idf)%W_V(iv)%datasz_in(1:3)
      CALL mathop (W_F(idf)%W_V(iv)%sopp(1),nbpt_in3,pdata_3d, &
 &                 missing_val,nbindex,nindex, &
 &                 W_F(idf)%W_V(iv)%scal(1),nbpt_out,tbf_1)
    ENDIF
    CALL histwrite_real (idf,iv,pitau,nbpt_out, &
 &            tbf_1,nbindex,nindex,do_oper,do_write)
  ENDIF
!-
! 6.0 Manage time steps
!-
  IF ((TRIM(tmp_opp) /= "once").AND.(TRIM(tmp_opp) /= "never")) THEN
    W_F(idf)%W_V(iv)%last_opp_chk = pitau
    W_F(idf)%W_V(iv)%last_wrt_chk = pitau
  ELSE
    W_F(idf)%W_V(iv)%last_opp_chk = -99
    W_F(idf)%W_V(iv)%last_wrt_chk = -99
  ENDIF
!-----------------------
END SUBROUTINE histw_rnd
!===
SUBROUTINE histwrite_real &
 & (idf,iv,pitau,nbdpt,tbf_1,nbindex,nindex,do_oper,do_write)
!---------------------------------------------------------------------
!- This subroutine is internal and does the calculations and writing
!- if needed. At a later stage it should be split into an operation
!- and writing subroutines.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: idf,pitau,iv, &
 &                      nbindex,nindex(nbindex),nbdpt
  REAL,DIMENSION(:)  :: tbf_1
  LOGICAL,INTENT(IN) :: do_oper,do_write
!-
  INTEGER :: tsz,nfid,nvid,iret,itax,io,nbin,nbout
  INTEGER :: nx,ny,nz,ky,kz,kt,kc
  INTEGER,DIMENSION(4) :: corner,edges
  INTEGER :: itime
!-
  REAL :: rtime
  REAL,DIMENSION(2) :: t_bnd
  CHARACTER(LEN=7) :: tmp_opp
  REAL,ALLOCATABLE,DIMENSION(:),SAVE :: tbf_2
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "histwrite 0.0 : VAR : ",W_F(idf)%W_V(iv)%v_name
    WRITE(*,*) "histwrite 0.0 : nbindex :",nbindex
    WRITE(*,*) "histwrite 0.0 : nindex  :",nindex(1:MIN(3,nbindex)),'...'
  ENDIF
!-
! The sizes which can be encoutered
!-
  tsz =  W_F(idf)%W_V(iv)%zsize(1) &
 &      *W_F(idf)%W_V(iv)%zsize(2) &
 &      *W_F(idf)%W_V(iv)%zsize(3)
!-
! 1.0 We allocate and the temporary space needed for operations.
!     The buffers are only deallocated when more space is needed.
!     This reduces the umber of allocates but increases memory needs.
!-
  IF (.NOT.ALLOCATED(tbf_2)) THEN
    IF (l_dbg) THEN
      WRITE(*,*) "histwrite_real 1.1 allocate tbf_2 ",SIZE(tbf_1)
    ENDIF
    ALLOCATE(tbf_2(W_F(idf)%W_V(iv)%datasz_max))
  ELSE IF (W_F(idf)%W_V(iv)%datasz_max > SIZE(tbf_2)) THEN
    IF (l_dbg) THEN
      WRITE(*,*) "histwrite_real 1.2 re-allocate tbf_2 : ", &
     & SIZE(tbf_1)," instead of ",SIZE(tbf_2)
    ENDIF
    DEALLOCATE(tbf_2)
    ALLOCATE(tbf_2(W_F(idf)%W_V(iv)%datasz_max))
  ENDIF
!-
  rtime = pitau*W_F(idf)%deltat
  tmp_opp = W_F(idf)%W_V(iv)%topp
!-
! 3.0 Do the operations or transfer the slab of data into tbf_1
!-
  IF (l_dbg) THEN
    WRITE(*,*) "histwrite: 3.0",idf
  ENDIF
!-
! 3.1 DO the Operations only if needed
!-
  IF (do_oper) THEN
    nbout = nbdpt
!-
!-- 3.4 We continue the sequence of operations
!--     we started in the interface routine
!-
    DO io=2,W_F(idf)%W_V(iv)%nbopp,2
      nbin = nbout
      nbout = W_F(idf)%W_V(iv)%datasz_max
      CALL mathop(W_F(idf)%W_V(iv)%sopp(io),nbin,tbf_1, &
 &      missing_val,nbindex,nindex,W_F(idf)%W_V(iv)%scal(io), &
 &      nbout,tbf_2)
      IF (l_dbg) THEN
        WRITE(*,*) &
 &       "histwrite: 3.4a nbout : ",nbin,nbout,W_F(idf)%W_V(iv)%sopp(io)
      ENDIF
!-
      nbin = nbout
      nbout = W_F(idf)%W_V(iv)%datasz_max
      CALL mathop(W_F(idf)%W_V(iv)%sopp(io+1),nbin,tbf_2, &
 &      missing_val,nbindex,nindex,W_F(idf)%W_V(iv)%scal(io+1), &
 &      nbout,tbf_1)
      IF (l_dbg) THEN
        WRITE(*,*) &
 &     "histwrite: 3.4b nbout : ",nbin,nbout,W_F(idf)%W_V(iv)%sopp(io+1)
      ENDIF
    ENDDO
!-
!   3.5 Zoom into the data
!-
    IF (l_dbg) THEN
      WRITE(*,*) &
 &     "histwrite: 3.5 size(tbf_1) : ",SIZE(tbf_1)
      WRITE(*,*) &
 &     "histwrite: 3.5 slab in X :", &
 &     W_F(idf)%W_V(iv)%zorig(1),W_F(idf)%W_V(iv)%zsize(1)
      WRITE(*,*) &
 &     "histwrite: 3.5 slab in Y :", &
 &     W_F(idf)%W_V(iv)%zorig(2),W_F(idf)%W_V(iv)%zsize(2)
      WRITE(*,*) &
 &     "histwrite: 3.5 slab in Z :", &
 &     W_F(idf)%W_V(iv)%zorig(3),W_F(idf)%W_V(iv)%zsize(3)
      WRITE(*,*) &
 &     "histwrite: 3.5 slab of input:", &
 &     W_F(idf)%W_V(iv)%scsize(1), &
 &     W_F(idf)%W_V(iv)%scsize(2), &
 &     W_F(idf)%W_V(iv)%scsize(3)
    ENDIF
!---
!-- We have to consider blocks of contiguous data
!---
    nx=MAX(W_F(idf)%W_V(iv)%zsize(1),1)
    ny=MAX(W_F(idf)%W_V(iv)%zsize(2),1)
    nz=MAX(W_F(idf)%W_V(iv)%zsize(3),1)
    IF     (     (W_F(idf)%W_V(iv)%zorig(1) == 1) &
 &          .AND.(   W_F(idf)%W_V(iv)%zsize(1) &
 &                == W_F(idf)%W_V(iv)%scsize(1)) &
 &          .AND.(W_F(idf)%W_V(iv)%zorig(2) == 1) &
 &          .AND.(   W_F(idf)%W_V(iv)%zsize(2) &
 &                == W_F(idf)%W_V(iv)%scsize(2))) THEN
      kt = (W_F(idf)%W_V(iv)%zorig(3)-1)*nx*ny
      tbf_2(1:nx*ny*nz) = tbf_1(kt+1:kt+nx*ny*nz)
    ELSEIF (     (W_F(idf)%W_V(iv)%zorig(1) == 1) &
 &          .AND.(   W_F(idf)%W_V(iv)%zsize(1) &
 &                == W_F(idf)%W_V(iv)%scsize(1))) THEN
      kc = -nx*ny
      DO kz=W_F(idf)%W_V(iv)%zorig(3),W_F(idf)%W_V(iv)%zorig(3)+nz-1
        kc = kc+nx*ny
        kt = ( (kz-1)*W_F(idf)%W_V(iv)%scsize(2) &
 &            +W_F(idf)%W_V(iv)%zorig(2)-1)*nx
        tbf_2(kc+1:kc+nx*ny) = tbf_1(kt+1:kt+nx*ny)
      ENDDO
    ELSE
      kc = -nx
      DO kz=W_F(idf)%W_V(iv)%zorig(3),W_F(idf)%W_V(iv)%zorig(3)+nz-1
        DO ky=W_F(idf)%W_V(iv)%zorig(2),W_F(idf)%W_V(iv)%zorig(2)+ny-1
          kc = kc+nx
          kt = ((kz-1)*W_F(idf)%W_V(iv)%scsize(2)+ky-1) &
 &            *W_F(idf)%W_V(iv)%scsize(1) &
 &            +W_F(idf)%W_V(iv)%zorig(1)-1
          tbf_2(kc+1:kc+nx) = tbf_1(kt+1:kt+nx)
        ENDDO
      ENDDO
    ENDIF
!-
!-- 4.0 Get the min and max of the field
!-
    IF (l_dbg) THEN
      WRITE(*,*) "histwrite: 4.0 tbf_1",idf,iv, &
 &      TRIM(tmp_opp),' ---- ',LEN_TRIM(tmp_opp),nbindex
    ENDIF
!-
    IF (W_F(idf)%W_V(iv)%hist_calc_rng) THEN
      W_F(idf)%W_V(iv)%hist_minmax(1) = &
 &      MIN(W_F(idf)%W_V(iv)%hist_minmax(1), &
 &      MINVAL(tbf_2(1:tsz),MASK=tbf_2(1:tsz) /= missing_val))
      W_F(idf)%W_V(iv)%hist_minmax(2) = &
 &      MAX(W_F(idf)%W_V(iv)%hist_minmax(2), &
 &      MAXVAL(tbf_2(1:tsz),MASK=tbf_2(1:tsz) /= missing_val))
    ENDIF
!-
!-- 5.0 Do the operations if needed. In the case of instantaneous
!--     output we do not transfer to the time_buffer.
!-
    IF (l_dbg) THEN
      WRITE(*,*) "histwrite: 5.0 idf : ",idf," iv : ",iv," tsz : ",tsz
    ENDIF
!-
    IF (     (TRIM(tmp_opp) /= "inst") &
 &      .AND.(TRIM(tmp_opp) /= "once") ) THEN
      CALL moycum(tmp_opp,tsz,W_F(idf)%W_V(iv)%t_bf, &
 &           tbf_2,W_F(idf)%W_V(iv)%nb_opp)
    ENDIF
!-
    W_F(idf)%W_V(iv)%last_opp = pitau
    W_F(idf)%W_V(iv)%nb_opp = W_F(idf)%W_V(iv)%nb_opp+1
!-
  ENDIF
!-
! 6.0 Write to file if needed
!-
  IF (l_dbg) WRITE(*,*) "histwrite: 6.0",idf
!-
  IF (do_write) THEN
!-
    nfid = W_F(idf)%ncfid
    nvid = W_F(idf)%W_V(iv)%ncvid
!-
!-- 6.1 Do the operations that are needed before writting
!-
    IF (l_dbg) WRITE(*,*) "histwrite: 6.1",idf
!-
    IF (     (TRIM(tmp_opp) /= "inst") &
 &      .AND.(TRIM(tmp_opp) /= "once") ) THEN
      t_bnd(1:2) = (/ W_F(idf)%W_V(iv)%last_wrt*W_F(idf)%deltat,rtime /)
      rtime = (t_bnd(1)+t_bnd(2))/2.0
    ENDIF
!-
!-- 6.2 Add a value to the time axis of this variable if needed
!-
    IF (     (TRIM(tmp_opp) /= "l_max") &
 &      .AND.(TRIM(tmp_opp) /= "l_min") &
 &      .AND.(TRIM(tmp_opp) /= "once") ) THEN
!-
      IF (l_dbg) WRITE(*,*) "histwrite: 6.2",idf
!-
      itax  = W_F(idf)%W_V(iv)%t_axid
      itime = W_F(idf)%W_V(iv)%nb_wrt+1
!-
      IF (W_F(idf)%W_V(itax)%tax_last < itime) THEN
        iret = NF90_PUT_VAR (nfid,W_F(idf)%W_V(itax)%tdimid, &
 &               (/ rtime /),start=(/ itime /),count=(/ 1 /))
        IF (W_F(idf)%W_V(itax)%tbndid > 0) THEN
          iret = NF90_PUT_VAR (nfid,W_F(idf)%W_V(itax)%tbndid, &
 &                 t_bnd,start=(/ 1,itime /),count=(/ 2,1 /))
        ENDIF
        W_F(idf)%W_V(itax)%tax_last = itime
      ENDIF
    ELSE
      itime=1
    ENDIF
!-
!-- 6.3 Write the data. Only in the case of instantaneous output
!       we do not write the buffer.
!-
    IF (l_dbg) THEN
      WRITE(*,*) "histwrite: 6.3",idf,nfid,nvid,iv,itime
    ENDIF
!-
    IF (W_F(idf)%W_V(iv)%scsize(3) == 1) THEN
      IF (W_F(idf)%regular) THEN
        corner(1:4) = (/ 1,1,itime,0 /)
        edges(1:4) = (/ W_F(idf)%W_V(iv)%zsize(1), &
 &                      W_F(idf)%W_V(iv)%zsize(2),1,0 /)
      ELSE
        corner(1:4) = (/ 1,itime,0,0 /)
        edges(1:4) = (/ W_F(idf)%W_V(iv)%zsize(1),1,0,0 /)
      ENDIF
    ELSE
      IF (W_F(idf)%regular) THEN
        corner(1:4) = (/ 1,1,1,itime /)
        edges(1:4) = (/ W_F(idf)%W_V(iv)%zsize(1), &
 &                      W_F(idf)%W_V(iv)%zsize(2), &
 &                      W_F(idf)%W_V(iv)%zsize(3),1 /)
      ELSE
        corner(1:4) = (/ 1,1,itime,0 /)
        edges(1:4) = (/ W_F(idf)%W_V(iv)%zsize(1), &
 &                      W_F(idf)%W_V(iv)%zsize(3),1,0 /)
      ENDIF
    ENDIF
!-
    IF (     (TRIM(tmp_opp) /= "inst") &
 &      .AND.(TRIM(tmp_opp) /= "once") ) THEN
      iret = NF90_PUT_VAR (nfid,nvid,W_F(idf)%W_V(iv)%t_bf, &
 &                         start=corner(1:4),count=edges(1:4))
    ELSE
      iret = NF90_PUT_VAR (nfid,nvid,tbf_2, &
 &                         start=corner(1:4),count=edges(1:4))
    ENDIF
!-
    W_F(idf)%W_V(iv)%last_wrt = pitau
    W_F(idf)%W_V(iv)%nb_wrt = W_F(idf)%W_V(iv)%nb_wrt+1
    W_F(idf)%W_V(iv)%nb_opp = 0
!---
!   After the write the file can be synchronized so that no data is
!   lost in case of a crash. This feature gives up on the benefits of
!   buffering and should only be used in debuging mode. A flag is
!   needed here to switch to this mode.
!---
!   iret = NF90_SYNC (nfid)
!-
  ENDIF
!----------------------------
END SUBROUTINE histwrite_real
!===
SUBROUTINE histvar_seq (idf,pvarname,idv)
!---------------------------------------------------------------------
!- This subroutine optimize the search for the variable in the table.
!- In a first phase it will learn the succession of the variables
!- called and then it will use the table to guess what comes next.
!- It is the best solution to avoid lengthy searches through array
!- vectors.
!-
!- ARGUMENTS :
!-
!- idf      : id of the file on which we work
!- pvarname : The name of the variable we are looking for
!- idv      : The var id we found
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(in)  :: idf
  CHARACTER(LEN=*),INTENT(IN) :: pvarname
  INTEGER,INTENT(out) :: idv
!-
  LOGICAL,SAVE :: learning(nb_files_max)=.TRUE.
  INTEGER,SAVE :: overlap(nb_files_max) = -1
  INTEGER,SAVE :: varseq(nb_files_max,nb_var_max*3)
  INTEGER,SAVE :: varseq_len(nb_files_max) = 0
  INTEGER,SAVE :: varseq_pos(nb_files_max)
  INTEGER,SAVE :: varseq_err(nb_files_max) = 0
  INTEGER      :: ib,sp,nn,pos
  CHARACTER(LEN=70) :: str70
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) 'histvar_seq, start of the subroutine :',learning(idf)
  ENDIF
!-
  IF (learning(idf)) THEN
!-
!-- 1.0 We compute the length over which we are going
!--     to check the overlap
!-
    IF (overlap(idf) <= 0) THEN
      IF (W_F(idf)%n_var > 6) THEN
        overlap(idf) = W_F(idf)%n_var/3*2
      ELSE
        overlap(idf) = W_F(idf)%n_var
      ENDIF
    ENDIF
!-
!-- 1.1 Find the position of this string
!-
    CALL find_str (W_F(idf)%W_V(1:W_F(idf)%n_var)%v_name,pvarname,pos)
    IF (pos > 0) THEN
      idv = pos
    ELSE
      CALL ipslerr (3,"histvar_seq", &
 &      'The name of the variable you gave has not been declared', &
 &      'You should use subroutine histdef for declaring variable', &
 &      TRIM(pvarname))
    ENDIF
!-
!-- 1.2 If we have not given up we store the position
!--     in the sequence of calls
!-
    IF (varseq_err(idf) >= 0) THEN
      sp = varseq_len(idf)+1
      IF (sp <= nb_var_max*3) THEN
        varseq(idf,sp) = idv
        varseq_len(idf) = sp
      ELSE
        CALL ipslerr (2,"histvar_seq",&
 &       'The learning process has failed and we give up. '// &
 &       'Either you sequence is',&
 &       'too complex or I am too dumb. '// &
 &       'This will only affect the efficiency',&
 &       'of your code. Thus if you wish to save time'// &
 &       ' contact the IOIPSL team. ')
        WRITE(*,*) 'The sequence we have found up to now :'
        WRITE(*,*) varseq(idf,1:sp-1)
        varseq_err(idf) = -1
      ENDIF
!-
!---- 1.3 Check if we have found the right overlap
!-
      IF (varseq_len(idf) >= overlap(idf)*2) THEN
!-
!------ We skip a few variables if needed as they could come
!------ from the initialisation of the model.
!-
        DO ib = 0,sp-overlap(idf)*2
          IF ( learning(idf) .AND.&
            & SUM(ABS(varseq(idf,ib+1:ib+overlap(idf)) -&
            & varseq(idf,sp-overlap(idf)+1:sp))) == 0 ) THEN
            learning(idf) = .FALSE.
            varseq_len(idf) = sp-overlap(idf)-ib
            varseq_pos(idf) = overlap(idf)+ib
            varseq(idf,1:varseq_len(idf)) = &
 &            varseq(idf,ib+1:ib+varseq_len(idf))
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ELSE
!-
!-- 2.0 Now we know how the calls to histwrite are sequenced
!--     and we can get a guess at the var ID
!-
    nn = varseq_pos(idf)+1
    IF (nn > varseq_len(idf)) nn = 1
!-
    idv = varseq(idf,nn)
!-
    IF (TRIM(W_F(idf)%W_V(idv)%v_name) /= TRIM(pvarname)) THEN
      CALL find_str (W_F(idf)%W_V(1:W_F(idf)%n_var)%v_name,pvarname,pos)
      IF (pos > 0) THEN
        idv = pos
      ELSE
        CALL ipslerr (3,"histvar_seq", &
 &  'The name of the variable you gave has not been declared',&
 &  'You should use subroutine histdef for declaring variable', &
 &  TRIM(pvarname))
      ENDIF
      varseq_err(idf) = varseq_err(idf)+1
    ELSE
!-
!---- We only keep the new position if we have found the variable
!---- this way. This way an out of sequence call to histwrite does
!---- not defeat the process.
!-
      varseq_pos(idf) = nn
    ENDIF
!-
!!$    IF (varseq_err(idf) >= 10) THEN
!!$      WRITE(str70,'("for file ",I3)') idf
!!$      CALL ipslerr (2,"histvar_seq", &
!!$ &  'There were 10 errors in the learned sequence of variables',&
!!$ &  str70,'This looks like a bug, please report it.')
!!$         varseq_err(idf) = 0
!!$    ENDIF
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) &
 &   'histvar_seq, end of the subroutine :',TRIM(pvarname),idv
  ENDIF
!-------------------------
END SUBROUTINE histvar_seq
!===
SUBROUTINE histsync (idf)
!---------------------------------------------------------------------
!- This subroutine will synchronise all
!- (or one if defined) opened files.
!-
!- VERSION
!-
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
! idf  : optional argument for fileid
  INTEGER,INTENT(in),OPTIONAL :: idf
!-
  INTEGER :: ifile,iret,i_s,i_e
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->histsync"
  ENDIF
!-
  IF (PRESENT(idf)) THEN
    IF ( (idf >= 1).AND.(idf <= nb_files_max) ) THEN
      IF (W_F(idf)%ncfid > 0) THEN
        i_s = idf
        i_e = idf
      ELSE
        i_s = 1
        i_e = 0
        CALL ipslerr (2,'histsync', &
 &       'Unable to synchronise the file :','probably','not opened')
      ENDIF
    ELSE
      CALL ipslerr (3,'histsync','Invalid file identifier',' ',' ')
    ENDIF
  ELSE
    i_s = 1
    i_e = nb_files_max
  ENDIF
!-
  DO ifile=i_s,i_e
    IF (W_F(ifile)%ncfid > 0) THEN
      IF (l_dbg) THEN
        WRITE(*,*) '  histsync - synchronising file number ',ifile
      ENDIF
      iret = NF90_SYNC(W_F(ifile)%ncfid)
    ENDIF
  ENDDO
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-histsync"
  ENDIF
!----------------------
END SUBROUTINE histsync
!===
SUBROUTINE histclo (idf)
!---------------------------------------------------------------------
!- This subroutine will close all (or one if defined) opened files
!-
!- VERSION
!-
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
! idf  : optional argument for fileid
  INTEGER,INTENT(in),OPTIONAL :: idf
!-
  INTEGER :: ifile,nfid,nvid,iret,iv,i_s,i_e
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->histclo"
  ENDIF
!-
  IF (PRESENT(idf)) THEN
    IF ( (idf >= 1).AND.(idf <= nb_files_max) ) THEN
      IF (W_F(idf)%ncfid > 0) THEN
        i_s = idf
        i_e = idf
      ELSE
        i_s = 1
        i_e = 0
        CALL ipslerr (2,'histclo', &
 &       'Unable to close the file :','probably','not opened')
      ENDIF
    ELSE
      CALL ipslerr (3,'histclo','Invalid file identifier',' ',' ')
    ENDIF
  ELSE
    i_s = 1
    i_e = nb_files_max
  ENDIF
!-
  DO ifile=i_s,i_e
    IF (W_F(ifile)%ncfid > 0) THEN
      IF (l_dbg) THEN
        WRITE(*,*) '  histclo - closing specified file number :',ifile
      ENDIF
      nfid = W_F(ifile)%ncfid
      iret = NF90_REDEF(nfid)
!-----
!---- 1. Loop on the number of variables to add some final information
!-----
      IF (l_dbg) THEN
        WRITE(*,*) '  Entering loop on vars : ',W_F(ifile)%n_var
      ENDIF
      DO iv=1,W_F(ifile)%n_var
!------ Extrema
        IF (W_F(ifile)%W_V(iv)%hist_wrt_rng) THEN
          IF (l_dbg) THEN
            WRITE(*,*) 'min value for file :',ifile,' var n. :',iv, &
 &                     ' is : ',W_F(ifile)%W_V(iv)%hist_minmax(1)
            WRITE(*,*) 'max value for file :',ifile,' var n. :',iv, &
 &                     ' is : ',W_F(ifile)%W_V(iv)%hist_minmax(2)
          ENDIF
          IF (W_F(ifile)%W_V(iv)%hist_calc_rng) THEN
!---------- Put the min and max values on the file
            nvid = W_F(ifile)%W_V(iv)%ncvid
            IF (W_F(ifile)%W_V(iv)%v_typ == hist_r8) THEN
              iret = NF90_PUT_ATT(nfid,nvid,'valid_min', &
 &                     REAL(W_F(ifile)%W_V(iv)%hist_minmax(1),KIND=8))
              iret = NF90_PUT_ATT(nfid,nvid,'valid_max', &
 &                     REAL(W_F(ifile)%W_V(iv)%hist_minmax(2),KIND=8))
            ELSE
              iret = NF90_PUT_ATT(nfid,nvid,'valid_min', &
 &                     REAL(W_F(ifile)%W_V(iv)%hist_minmax(1),KIND=4))
              iret = NF90_PUT_ATT(nfid,nvid,'valid_max', &
 &                     REAL(W_F(ifile)%W_V(iv)%hist_minmax(2),KIND=4))
            ENDIF
          ENDIF
        ENDIF
!------ Time-Buffers
        IF (ASSOCIATED(W_F(ifile)%W_V(iv)%t_bf)) THEN
          DEALLOCATE(W_F(ifile)%W_V(iv)%t_bf)
        ENDIF
!------ Reinitialize the sizes
        W_F(ifile)%W_V(iv)%datasz_in(:) = -1
        W_F(ifile)%W_V(iv)%datasz_max = -1
      ENDDO
!-----
!---- 2. Close the file
!-----
      IF (l_dbg) WRITE(*,*) '  close file :',nfid
      iret = NF90_CLOSE(nfid)
      W_F(ifile)%ncfid = -1
      W_F(ifile)%dom_id_svg = -1
    ENDIF
  ENDDO
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-histclo"
  ENDIF
!---------------------
END SUBROUTINE histclo
!===
SUBROUTINE ioconf_modname (str)
!---------------------------------------------------------------------
!- This subroutine allows to configure the name
!- of the model written into the file
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(IN) :: str
!---------------------------------------------------------------------
  IF (.NOT.lock_modname) THEN
    model_name = str(1:MIN(LEN_TRIM(str),80))
    lock_modname = .TRUE.
  ELSE
    CALL ipslerr (2,"ioconf_modname", &
   &  'The model name can only be changed once and only', &
   &  'before it is used. It is now set to :',model_name)
  ENDIF
!----------------------------
END SUBROUTINE ioconf_modname
!-
!===
!-
!-----------------
END MODULE histcom
