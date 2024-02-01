MODULE restcom
!-
!$Id: restcom.f90 2281 2010-10-15 14:21:13Z smasson $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!-
USE netcdf
!-
USE errioipsl, ONLY : ipslerr,ipsldbg
USE stringop
USE calendar
USE mathelp
USE fliocom,   ONLY : flio_dom_file,flio_dom_att
!-
IMPLICIT NONE
!-
PRIVATE
!-
PUBLIC :: &
 &  restini, restget, restput, restclo, &
 &  ioconf_setatt, ioget_vname, ioconf_expval, &
 &  ioget_expval, ioget_vdim
!-
INTERFACE restput
  MODULE PROCEDURE &
 &  restput_r3d, restput_r2d, restput_r1d, &
 &  restput_opp_r2d, restput_opp_r1d
END INTERFACE
!-
INTERFACE restget
  MODULE PROCEDURE &
 &  restget_r3d,restget_r2d,restget_r1d, &
 &  restget_opp_r2d,restget_opp_r1d
END INTERFACE
!-
! We do not use allocatable arrays because these sizes are safe
! and we do not know from start how many variables will be in
! the out file.
!-
  INTEGER,PARAMETER :: &
 &  max_var=500, max_file=50, max_dim=NF90_MAX_VAR_DIMS
!-
  CHARACTER(LEN=9),SAVE :: calend_str='unknown'
!-
! The IDs of the netCDF files are going in pairs.
! The input one (netcdf_id(?,1)) and the output one (netcdf_id(?,2))
!-
  INTEGER,SAVE :: nb_fi = 0
  INTEGER,DIMENSION(max_file,2),SAVE :: netcdf_id = -1
!-
! Description of the content of the 'in' files and the 'out' files.
!   Number of variables   : nbvar_*
!   Number of dimensions  : nbdim_*
!   ID of the time axis   : tdimid_*
!-
  INTEGER,SAVE :: nbvar_in(max_file), nbvar_out(max_file)
  INTEGER,SAVE :: tdimid_in(max_file), tdimid_out(max_file)
!-
! Variables for one or the other file
!-
! Number of dimensions in the input file               : nbdim_in
! Number of variables read so far from the input file  : nbvar_read
! Type of variable read from the input file            : vartyp_in
!   (Could be used later to test if we have a restart file)
!-
  INTEGER,SAVE :: nbdim_in(max_file), nbvar_read(max_file)
  INTEGER,SAVE :: vartyp_in(max_file, max_var)
!-
! Time step and time origine in the input file.
!-
  REAL,DIMENSION(max_file),SAVE :: deltat,timeorig
!-
! Description of the axes in the output file
!-
!   tstp_out  : Index on the tie axis currently beeing written
!   itau_out  : Time step which is written on this index of the file
!-
  INTEGER,DIMENSION(max_file),SAVE :: tstp_out,itau_out
!-
! Description of the axes in the output file
!-
! For the ?ax_infs variable the following order is used :
!   ?ax_infs (if,in,1) = size of axis
!   ?ax_infs (if,in,2) = id of dimension
! Number of x,y and z axes in the output file :
!   ?ax_nb(if)
!-
  INTEGER,DIMENSION(max_file,max_dim,2),SAVE :: &
 &  xax_infs,yax_infs,zax_infs
  INTEGER,DIMENSION(max_file),SAVE :: &
 &  xax_nb=0,yax_nb=0,zax_nb=0
!-
! Description of the time axes in the input and output files
!-
!   ID of the variable which contains the itaus :
!     tind_varid_*
!   ID of the variables which contains the seconds since date :
!     tax_varid_*
!   Size of the time axis in the input file :
!     tax_size_in
!-
  INTEGER,SAVE :: tind_varid_in(max_file),  tax_varid_in(max_file), &
 &                tind_varid_out(max_file), tax_varid_out(max_file)
  INTEGER,SAVE :: tax_size_in(max_file)=1
!-
! The two time axes we have in the input file :
!   t_index   : dates in itaus
!               (thus the variable has a tstep_sec attribute)
!   t_julian  : Julian days of the time axis
!-
  INTEGER,SAVE,ALLOCATABLE :: t_index(:,:)
  REAL,SAVE,ALLOCATABLE :: t_julian(:,:)
!-
! Here we save a number of informations on the variables
! in the files we are handling
!-
! Name of variables :                                    varname_*
! ID of the variables :                                  varid_*
! Number of dimensions of the variable :                 varnbdim_*
! Dimensions which are used for the variable :           vardims_*
! Number of attributes for a variables :                 varatt_*
! A flag which markes the variables we have worked on :  touched_*
!-
  CHARACTER(LEN=20),DIMENSION(max_file,max_var),SAVE :: &
 &  varname_in,varname_out
  INTEGER,DIMENSION(max_file,max_var),SAVE :: &
 &  varid_in,varid_out,varnbdim_in,varatt_in
  INTEGER,DIMENSION(max_file,max_var,max_dim),SAVE :: &
 &  vardims_in
  LOGICAL,DIMENSION(max_file,max_var),SAVE :: &
 &  touched_in,touched_out
!-
  CHARACTER(LEN=120),SAVE :: indchfun= 'scatter, fill, gather, coll'
  REAL,PARAMETER :: missing_val=1.e20
!                or HUGE(1.0) (maximum real number)
!-
! The default value we will use for variables
! which are not present in the restart file
!-
  REAL,SAVE :: val_exp =  999999.
  LOGICAL,SAVE :: lock_valexp = .FALSE.
!-
! Temporary variables in which we store the attributed which are going
! to be given to a new variable which is going to be defined.
!-
  CHARACTER(LEN=80),SAVE :: rest_units='XXXXX',rest_lname='XXXXX'
!-
! For allocations
!-
  REAL,ALLOCATABLE,DIMENSION(:),SAVE :: buff_tmp1,buff_tmp2
!-
!===
CONTAINS
!===
!-
SUBROUTINE restini &
 & (fnamein,iim,jjm,lon,lat,llm,lev, &
 &  fnameout,itau,date0,dt,fid,owrite_time_in,domain_id)
!---------------------------------------------------------------------
!- This subroutine sets up all the restart process.
!- It will call the subroutine which opens the input
!- and output files.
!- The time step (itau), date of origine (date0) and time step are
!- READ from the input file.
!- A file ID, which is common to the input and output file is returned
!-
!- If fnamein = fnameout then the same file is used for the reading
!- the restart conditions and writing the new restart.
!-
!- A special mode can be switched in with filename='NONE'.
!- This means that no restart file is present.
!- Usefull for creating the first restart file
!- or to get elements in a file without creating an output file.
!-
!- A mode needs to be written in which itau, date0 and dt
!- are given to the restart process and thus
!- written into the output restart file.
!-
!- INPUT
!-
!- fnamein  : name of the file for the restart
!- iim      : Dimension in x
!- jjm      : Dimension in y
!- lon      : Longitude in the x,y domain
!- lat      : Latitude in the x,y domain
!- llm      : Dimension in the vertical
!- lev      : Positions of the levels
!- fnameout :
!-
!- OUTPUT
!-
!- itau     : Time step of the restart file and at which the model
!-            should restart
!- date0    : Time at which itau = 0
!- dt       : time step in seconds between two succesiv itaus
!- fid      : File identification of the restart file
!-
!- Optional INPUT arguments
!-
!- owrite_time_in : logical  argument which allows to
!-                  overwrite the time in the restart file
!- domain_id      : Domain identifier
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(IN) :: fnamein,fnameout
  INTEGER :: iim,jjm,llm,fid,itau
  REAL :: lon(iim,jjm),lat(iim,jjm),lev(llm)
  REAL :: date0,dt
  LOGICAL,OPTIONAL :: owrite_time_in
  INTEGER,INTENT(IN),OPTIONAL :: domain_id
!-
  INTEGER :: ncfid
  REAL :: dt_tmp,date0_tmp
  LOGICAL :: l_fi,l_fo,l_rw
  LOGICAL :: overwrite_time
  CHARACTER(LEN=120) :: fname
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! 0.0 Prepare the configuration before opening any files
!-
  IF (.NOT.PRESENT(owrite_time_in)) THEN
    overwrite_time = .FALSE.
  ELSE
    overwrite_time = owrite_time_in
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) 'restini 0.0 : ',TRIM(fnamein),' , ',TRIM(fnameout)
  ENDIF
!-
  nb_fi = nb_fi+1
!-
  IF (nb_fi > max_file) THEN
    CALL ipslerr (3,'restini',&
 &   'Too many restart files are used. The problem can be',&
 &   'solved by increasing max_file in restcom.f90 ',&
 &   'and recompiling ioipsl.')
  ENDIF
!-
! 0.1 Define the open flags
!-
  l_fi = (TRIM(fnamein)  /= 'NONE')
  l_fo = (TRIM(fnameout) /= 'NONE')
  IF ((.NOT.l_fi).AND.(.NOT.l_fo)) THEN
    CALL ipslerr (3,'restini',&
 &   'Input and output file names are both to NONE.',&
 &   'It is probably an error.','Verify your logic.')
  ENDIF
  l_rw = l_fi.AND.l_fo.AND.(TRIM(fnamein) == TRIM(fnameout))
!-
  IF (l_dbg) THEN
    WRITE(*,*) 'restini 0.1 l_fi, l_fo, l_rw : ',l_fi,l_fo,l_rw
  ENDIF
!-
! 1.0 Open the input file.
!-
  IF (l_fi) THEN
!---
    IF (l_dbg) WRITE(*,*) 'restini 1.0 : Open input file'
!-- Add DOMAIN number and ".nc" suffix in file names if needed
    fname = fnamein
    CALL flio_dom_file (fname,domain_id)
!-- Open the file
    CALL restopenin (nb_fi,fname,l_rw,iim,jjm,lon,lat,llm,lev,ncfid)
    netcdf_id(nb_fi,1) = ncfid
!---
!-- 1.3 Extract the time information
!---
    IF (overwrite_time) THEN
       date0_tmp = date0
    ENDIF
    CALL restsett (dt_tmp,date0_tmp,itau,overwrite_time)
    IF (.NOT.overwrite_time) THEN
      dt = dt_tmp
      date0 = date0_tmp
    ENDIF
!---
  ELSE
!---
!-- 2.0 The case of a missing restart file is dealt with
!---
    IF (l_dbg) WRITE(*,*) 'restini 2.0'
!---
    IF (     (ALL(MINLOC(lon(:iim,:jjm)) == MAXLOC(lon(:iim,:jjm)))) &
        .AND.(iim > 1) ) THEN
      CALL ipslerr (3,'restini',&
        & 'For creating a restart file the longitudes of the',&
        & 'grid need to be provided to restini. This ',&
        & 'information is needed for the restart files')
    ENDIF
    IF (     (ALL(MINLOC(lat(:iim,:jjm)) == MAXLOC(lat(:iim,:jjm)))) &
        .AND.(jjm > 1) ) THEN
      CALL ipslerr (3,'restini',&
        & 'For creating a restart file the latitudes of the',&
        & 'grid need to be provided to restini. This ',&
        & 'information is needed for the restart files')
    ENDIF
    IF (     (ALL(MINLOC(lev(:llm)) == MAXLOC(lev(:llm)))) &
        .AND.(llm > 1) ) THEN
      CALL ipslerr (3,'restini',&
        & 'For creating a restart file the levels of the',&
        & 'grid need to be provided to restini. This',&
        & 'information is needed for the restart files')
    ENDIF
!---
!-- 2.2 Allocate the time axes and write the inputed variables
!---
    tax_size_in(nb_fi) = 1
    CALL rest_atim (l_dbg,'restini')
    t_index(nb_fi,1) = itau
    t_julian(nb_fi,1) = date0
  ENDIF
!-
  IF (l_fo.AND.(.NOT.l_rw)) THEN
!-- Add DOMAIN number and ".nc" suffix in file names if needed
    fname = fnameout
    CALL flio_dom_file (fname,domain_id)
!-- Open the file
    CALL restopenout &
      (nb_fi,fname,iim,jjm,lon,lat,llm,lev,dt,date0,ncfid,domain_id)
    netcdf_id(nb_fi,2) = ncfid
  ELSE IF (l_fi.AND.l_fo) THEN
    netcdf_id(nb_fi,2) = netcdf_id(nb_fi,1)
    varname_out(nb_fi,:) = varname_in(nb_fi,:)
    nbvar_out(nb_fi) = nbvar_in(nb_fi)
    tind_varid_out(nb_fi) = tind_varid_in(nb_fi)
    tax_varid_out(nb_fi) = tax_varid_in(nb_fi)
    varid_out(nb_fi,:) = varid_in(nb_fi,:)
    touched_out(nb_fi,:) = .TRUE.
  ENDIF
!-
! 2.3 Set the calendar for the run.
!     This should not produce any error message if
!     This does not mean any change in calendar
!     (to be modified in ioconf_calendar)
!-
  IF (l_dbg) THEN
    WRITE(*,*) 'restini 2.3 : Configure calendar if needed : ', &
                calend_str
  ENDIF
!-
  IF (INDEX(calend_str,'unknown') < 1) THEN
    CALL ioconf_calendar (calend_str)
    IF (l_dbg) THEN
      WRITE(*,*) 'restini 2.3b : new calendar : ',calend_str
    ENDIF
  ENDIF
!-
! Save some data in the module
!-
  deltat(nb_fi) = dt
!-
! Prepare the variables which will be returned
!-
  fid = nb_fi
  IF (l_dbg) THEN
    WRITE(*,*) 'SIZE of t_index :',SIZE(t_index), &
               SIZE(t_index,dim=1),SIZE(t_index,dim=2)
    WRITE(*,*) 't_index = ',t_index(fid,:)
  ENDIF
  itau = t_index(fid,1)
!-
  IF (l_dbg) WRITE(*,*) 'restini END'
!---------------------
END SUBROUTINE restini
!===
SUBROUTINE restopenin &
  (fid,fname,l_rw,iim,jjm,lon,lat,llm,lev,ncfid)
!---------------------------------------------------------------------
!- Opens the restart file and checks that it belongsd to the model.
!- This means that the coordinates of the model are compared to the
!- ones in the file.
!-
!- The number and name of variable in the file are exctracted. Also
!- the time details.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: fid,iim,jjm,llm
  CHARACTER(LEN=*),INTENT(IN) :: fname
  REAL :: lon(iim,jjm),lat(iim,jjm),lev(llm)
  LOGICAL,INTENT(IN) :: l_rw
  INTEGER,INTENT(OUT) :: ncfid
!-
  INTEGER,DIMENSION(max_dim) :: var_dims,dimlen
  INTEGER :: nb_dim,nb_var,id_unl,id,iv
  INTEGER :: iread,jread,lread,iret
  INTEGER :: lon_vid,lat_vid
  REAL :: lon_read(iim,jjm),lat_read(iim,jjm)
  REAL :: lev_read(llm)
  REAL :: mdlon,mdlat
  CHARACTER(LEN=80) :: units
  CHARACTER(LEN=NF90_max_name),DIMENSION(max_dim) :: dimname
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! If we reuse the same file for input and output
! then we open it in write mode
!-
  IF (l_rw) THEN; id = NF90_WRITE; ELSE; id = NF90_NOWRITE; ENDIF
  iret = NF90_OPEN(fname,id,ncfid)
  IF (iret /= NF90_NOERR) THEN
    CALL ipslerr (3,'restopenin','Could not open file :',fname,' ')
  ENDIF
!-
  IF (l_dbg) WRITE (*,*) "restopenin 0.0 ",TRIM(fname)
  iret = NF90_INQUIRE(ncfid,nDimensions=nb_dim, &
 &         nVariables=nb_var,unlimitedDimId=id_unl)
  tdimid_in(fid) = id_unl
!-
  IF (nb_dim > max_dim) THEN
    CALL ipslerr (3,'restopenin',&
      & 'More dimensions present in file that can be store',&
      & 'Please increase max_dim in the global variables ',&
      & 'in restcom.F90')
  ENDIF
  IF (nb_var > max_var) THEN
    CALL ipslerr (3,'restopenin',&
      & 'More variables present in file that can be store',&
      & 'Please increase max_var in the global variables ',&
      & 'in restcom.F90')
  ENDIF
!-
  nbvar_in(fid) = nb_var
  nbdim_in(fid) = nb_dim
  iread = -1; jread = -1; lread = -1;
  DO id=1,nb_dim
    iret = NF90_INQUIRE_DIMENSION(ncfid,id, &
 &           len=dimlen(id),name=dimname(id))
    IF (l_dbg) THEN
      WRITE (*,*) "restopenin 0.0 dimname",id,TRIM(dimname(id))
    ENDIF
    IF      (TRIM(dimname(id)) == 'x') THEN
      iread = dimlen(id)
      IF (l_dbg) WRITE (*,*) "iread",iread
    ELSE IF (TRIM(dimname(id)) == 'y') THEN
      jread = dimlen(id)
      IF (l_dbg) WRITE (*,*) "jread",jread
    ELSE IF (TRIM(dimname(id)) == 'z') THEN
      lread = dimlen(id)
      IF (l_dbg) WRITE (*,*) "lread",lread
    ENDIF
  ENDDO
!-
  IF (id_unl > 0) THEN
!---
!-- 0.1 If we are going to add values to this file
!--     we need to know where it ends
!--     We also need to have all the dimensions in the file
!---
    IF (l_rw) THEN
      tstp_out(fid) = dimlen(id_unl)
      itau_out(fid) = -1
      tdimid_out(fid) =  tdimid_in(fid)
      IF (l_dbg) THEN
        WRITE (*,*) &
 &       "restopenin 0.0 unlimited axis dimname", &
 &       dimname(id_unl),tstp_out(fid)
      ENDIF
!-----
      xax_nb(fid) = 0
      yax_nb(fid) = 0
      zax_nb(fid) = 0
!-----
      DO id=1,nb_dim
        IF      (dimname(id)(1:1) == 'x') THEN
          xax_nb(fid) = xax_nb(fid)+1
          xax_infs(fid,xax_nb(fid),1) = dimlen(id)
          xax_infs(fid,xax_nb(fid),2) = id
        ELSE IF (dimname(id)(1:1) == 'y') THEN
          yax_nb(fid) = yax_nb(fid)+1
          yax_infs(fid,yax_nb(fid),1) = dimlen(id)
          yax_infs(fid,yax_nb(fid),2) = id
        ELSE IF (dimname(id)(1:1) == 'z') THEN
          zax_nb(fid) = zax_nb(fid)+1
          zax_infs(fid,zax_nb(fid),1) = dimlen(id)
          zax_infs(fid,zax_nb(fid),2) = id
        ENDIF
      ENDDO
    ENDIF
  ELSE
!---
!-- Still need to find a method for dealing with this
!---
!  CALL ipslerr (3,'restopenin',&
!    & ' We do not deal yet with files without time axis.',' ',' ')
  ENDIF
!-
! 1.0 First let us check that we have the righ restart file
!-
  IF ((iread /= iim).OR.(jread /= jjm).OR.(lread /= llm)) THEN
    CALL ipslerr (3,'restopenin',&
 &    'The grid of the restart file does not correspond',&
 &    'to that of the model',' ')
  ENDIF
!-
! 2.0 Get the list of variables
!-
  IF (l_dbg) WRITE(*,*) 'restopenin 1.2'
!-
  lat_vid = -1
  lon_vid = -1
  tind_varid_in(fid) = -1
  tax_varid_in(fid) = -1
!-
  DO iv=1,nb_var
!---
    varid_in(fid,iv) = iv
    var_dims(:) = 0
    iret = NF90_INQUIRE_VARIABLE(ncfid,iv, &
 &           name=varname_in(fid,iv),xtype=vartyp_in(fid,iv), &
 &           ndims=varnbdim_in(fid,iv),dimids=var_dims, &
 &           nAtts=varatt_in(fid,iv))
!---
    DO id=1,varnbdim_in(fid,iv)
      iret = NF90_INQUIRE_DIMENSION &
 &             (ncfid,var_dims(id),len=vardims_in(fid,iv,id))
    ENDDO
!---
!-- 2.1 Read the units of the variable
!---
    units=''
    iret = NF90_GET_ATT(ncfid,iv,'units',units)
    CALL strlowercase (units)
    CALL cmpblank (units)
!---
!-- 2.2 Catch the time variables
!---
    IF (varnbdim_in(fid,iv) == 1) THEN
      IF (     (INDEX(units,'timesteps since') > 0) &
          .AND.(tind_varid_in(fid) < 0) ) THEN
        tind_varid_in(fid) = iv
        tax_size_in(fid) = vardims_in(fid,iv,1)
      ENDIF
      IF (     (INDEX(units,'seconds since') > 0) &
          .AND.(tax_varid_in(fid) < 0) ) THEN
        tax_varid_in(fid) = iv
        tax_size_in(fid) = vardims_in(fid,iv,1)
      ENDIF
    ENDIF
!---
!-- 2.3 Catch longitude and latitude variables
!---
    IF      (INDEX(units,'degrees_nort') > 0) THEN
      lat_vid = iv
    ELSE IF (INDEX(units,'degrees_east') > 0) THEN
      lon_vid = iv
    ENDIF
!---
  ENDDO
!-
! 2.4 None of the variables was yet read
!-
  nbvar_read(fid) = 0
  touched_in(fid,:) = .FALSE.
!-
! 3.0 Reading the coordinates from the input restart file
!-
  lon_read = missing_val
  lat_read = missing_val
!-
  IF (lon_vid < 0 .OR. lat_vid < 0) THEN
    CALL ipslerr (3,'restopenin',&
      & ' No variables containing longitude or latitude were ',&
      & ' found in the restart file.',' ')
  ELSE
    iret = NF90_GET_VAR(ncfid,lon_vid,lon_read)
    iret = NF90_GET_VAR(ncfid,lat_vid,lat_read)
!---
    IF (  (ABS( MAXVAL(lon(:,:)) &
 &             -MINVAL(lon(:,:))) < EPSILON(MAXVAL(lon(:,:)))) &
 &   .AND.(ABS( MAXVAL(lat(:,:)) &
 &             -MINVAL(lat(:,:))) < EPSILON(MAXVAL(lat(:,:)))) ) THEN
!-----
!---- 3.1 No longitude nor latitude are provided thus
!---- they are taken from the restart file
!-----
      lon(:,:) = lon_read(:,:)
      lat(:,:) = lat_read(:,:)
    ELSE
!-----
!---- 3.2 We check that the longitudes and latitudes
!----     in the file and the model are the same
!-----
      mdlon = MAXVAL(ABS(lon_read-lon))
      mdlat = MAXVAL(ABS(lat_read-lat))
!-----
!---- We can not test against epsilon here as the longitude
!---- can be stored at another precision in the netCDF file.
!---- The test here does not need to be very precise.
!-----
      IF (mdlon > 1.e-4 .OR. mdlat > 1.e-4) THEN
        CALL ipslerr (3,'restopenin',&
          & ' The longitude or latitude found in the restart ',&
          & ' file are not the same as the ones used in the model.',&
          & ' ')
      ENDIF
    ENDIF
  ENDIF
!------------------------
END SUBROUTINE restopenin
!===
SUBROUTINE restsett (timestep,date0,itau,owrite_time_in)
!---------------------------------------------------------------------
!- Here we get all the time information from the file.
!-
!- The time information can come in three forms :
!- -global attributes which give the time origine and the
!-  time step is taken from the input to restinit
!- -A physical time exists and thus the julian date from the
!-  input is used for positioning using the itau as input
!- -A time-step axis exists and itau is positioned on it.
!-
!- What takes precedence : the model
!-
!- itau     : Time step of the model
!-
!- Optional INPUT arguments
!-
!- owrite_time_in : logical  argument which allows to
!-                  overwrite the time in the restart file
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  REAL :: date0,timestep
  INTEGER :: itau
  LOGICAL,OPTIONAL :: owrite_time_in
!-
  INTEGER :: ncfid,iret,it,iax,iv
  CHARACTER(LEN=80) :: itau_orig,tax_orig,calendar
  CHARACTER(LEN=9) :: tmp_cal
  INTEGER :: year0,month0,day0,hours0,minutes0,seci
  REAL :: sec0,one_day,one_year,date0_ju,ttmp
  CHARACTER :: strc
  LOGICAL :: ow_time
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (PRESENT(owrite_time_in)) THEN
    ow_time = owrite_time_in
  ELSE
    ow_time = .FALSE.
  ENDIF
!-
  ncfid = netcdf_id(nb_fi,1)
!-
! Allocate the space we need for the time axes
!-
  CALL rest_atim (l_dbg,'restsett')
!-
! Get the calendar if possible. Else it will be gregorian.
!-
  IF (tax_size_in(nb_fi) > 0) THEN
    calendar = ' '
    iret = NF90_GET_ATT(ncfid,tax_varid_in(nb_fi),'calendar',calendar)
    IF (iret == NF90_NOERR) THEN
      CALL ioconf_calendar (calendar)
      IF (l_dbg) THEN
        WRITE(*,*) 'restsett : calendar of the restart ',calendar
      ENDIF
    ENDIF
  ENDIF
  CALL ioget_calendar (one_year,one_day)
  IF (l_dbg) THEN
    WRITE(*,*) 'one_year,one_day = ',one_year,one_day
  ENDIF
!-
  itau_orig = 'XXXXX'
  tax_orig  = 'XXXXX'
!-
! Get the time steps of the time axis if available on the restart file
!-
  IF (tind_varid_in(nb_fi) > 0) THEN
    IF (ow_time) THEN
      t_index(nb_fi,:) = itau
      IF (l_dbg) THEN
        WRITE(*,*) "nb_fi,t_index",nb_fi,t_index(nb_fi,:)
      ENDIF
      CALL ju2ymds (date0,year0,month0,day0,sec0)
      hours0 = NINT(sec0/3600)
      sec0 = sec0 - 3600 * hours0
      minutes0 = NINT(sec0 / 60)
      sec0 = sec0 - 60 * minutes0
      seci = NINT(sec0)
      strc=':'
      IF (l_dbg) THEN
        WRITE(*,*) date0
        WRITE (UNIT=itau_orig,FMT='(I4.4,5(A,I2.2))') &
 &       year0,'-',month0,'-',day0,' ',hours0,':',minutes0,':',seci
        WRITE(*,*) "itau_orig : ",itau_orig
      ENDIF
    ELSE
      iret = NF90_GET_VAR(ncfid,tind_varid_in(nb_fi),t_index(nb_fi,:))
      IF (l_dbg) THEN
        WRITE(*,*) "restsett, time axis : ",t_index(nb_fi,:)
      ENDIF
      iret = NF90_GET_ATT(ncfid,tind_varid_in(nb_fi),'units',itau_orig)
      itau_orig = &
 &      itau_orig(INDEX(itau_orig,'since')+6:LEN_TRIM(itau_orig))
      iret = &
 &      NF90_GET_ATT(ncfid,tind_varid_in(nb_fi),'tstep_sec',timestep)
!-----
!---- This time origin will dominate as it is linked to the time steps.
!-----
      READ (UNIT=itau_orig,FMT='(I4.4,5(A,I2.2))') &
 &      year0,strc,month0,strc,day0,strc, &
 &      hours0,strc,minutes0,strc,seci
      sec0 = REAL(seci)
      sec0 = hours0*3600.+minutes0*60.+sec0
      CALL ymds2ju (year0,month0,day0,sec0,date0)
    ENDIF
  ENDIF
!-
! If a julian day time axis is available then we get it
!-
  IF (tax_varid_in(nb_fi) > 0) THEN
    iret = NF90_GET_VAR(ncfid,tax_varid_in(nb_fi),t_julian(nb_fi,:))
    iret = NF90_GET_ATT(ncfid,tax_varid_in(nb_fi),'units',tax_orig)
    tax_orig = tax_orig(INDEX(tax_orig,'since')+6:LEN_TRIM(tax_orig))
    tmp_cal = ' '
    iret = NF90_GET_ATT(ncfid,tax_varid_in(nb_fi),'calendar',tmp_cal)
    IF (l_dbg) THEN
      WRITE(*,*) 'restsett : tmp_calendar of the restart ',tmp_cal
    ENDIF
!---
    CALL strlowercase (tmp_cal)
    IF (INDEX(calend_str,tmp_cal) < 1) THEN
      IF (INDEX(calend_str,'unknown') > 0) THEN
        calend_str = tmp_cal
      ELSE
        CALL ipslerr (2,'restsett', &
 &       ' In the restart files two different calendars were found.', &
 &       ' Please check the files you have used.',' ')
      ENDIF
    ENDIF
!---
!-- We need to transform that into julian days
!-- to get ride of the intial date.
!---
    IF (l_dbg) WRITE(*,*) 'tax_orig : ',TRIM(tax_orig)
    READ (UNIT=tax_orig,FMT='(I4.4,5(a,I2.2))') &
      year0,strc,month0,strc,day0,strc, &
      hours0,strc,minutes0,strc,seci
    sec0 = REAL(seci)
    sec0 = hours0*3600.+minutes0*60.+sec0
    CALL ymds2ju (year0,month0,day0,sec0,date0_ju)
    t_julian(nb_fi,:) = t_julian(nb_fi,:)/one_day+date0_ju
  ENDIF
!-
  IF (     (INDEX(itau_orig,'XXXXX') > 0) &
      .AND.(INDEX(tax_orig,'XXXXX')  < 1) ) THEN
!!- Compute the t_itau from the date read and the timestep in the input
  ENDIF
!-
  IF (     (INDEX(tax_orig,'XXXXX')  > 0) &
      .AND.(INDEX(itau_orig,'XXXXX') < 1) ) THEN
    DO it=1,tax_size_in(nb_fi)
      t_julian(nb_fi,it) = itau2date(t_index(nb_fi,it),date0,timestep)
    ENDDO
  ENDIF
!-
! If neither the indices or time is present then get global attributes
! This is for compatibility reasons and should not be used.
!-
  IF ((tax_varid_in(nb_fi) < 0).AND.(tind_varid_in(nb_fi) < 0)) THEN
    iax = -1
    DO iv=1,nbvar_in(nb_fi)
      IF (    (INDEX(varname_in(nb_fi,iv),'tsteps') > 0) &
 &        .OR.(INDEX(varname_in(nb_fi,iv),'time_steps') > 0)) THEN
        iax = iv
      ENDIF
    ENDDO
!---
    IF (iax < 0) THEN
      CALL ipslerr (3,'restsett',&
        & 'No time axis was found in the restart file. Please check',&
        & 'that it corresponds to the convention used in restsett',&
        & ' ')
    ELSE
      iret = NF90_GET_VAR(ncfid,tind_varid_in(nb_fi),t_index(nb_fi,:))
      iret = NF90_GET_ATT(ncfid,NF90_GLOBAL,'delta_tstep_sec',timestep)
      iret = NF90_GET_ATT(ncfid,NF90_GLOBAL,'year0',ttmp)
      year0 = NINT(ttmp)
      iret = NF90_GET_ATT(ncfid,NF90_GLOBAL,'month0',ttmp)
      month0 = NINT(ttmp)
      iret = NF90_GET_ATT(ncfid,NF90_GLOBAL,'day0',ttmp)
      day0 = NINT(ttmp)
      iret = NF90_GET_ATT(ncfid,NF90_GLOBAL,'sec0',sec0)
!---
      CALL ymds2ju (year0,month0,day0,sec0,date0)
      t_julian(nb_fi,1) = itau2date(t_index(nb_fi,1),date0,timestep)
    ENDIF
  ENDIF
!----------------------
END SUBROUTINE restsett
!===
SUBROUTINE restopenout &
  (fid,fname,iim,jjm, &
   lon,lat,llm,lev,timestep,date,ncfid,domain_id)
!---------------------------------------------------------------------
!- Opens the restart file for output.
!- The longitude and time variables are written.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: fid,iim,jjm,llm
  CHARACTER(LEN=*) :: fname
  REAL :: date,timestep
  REAL :: lon(iim,jjm),lat(iim,jjm),lev(llm)
  INTEGER,INTENT(OUT) :: ncfid
  INTEGER,INTENT(IN),OPTIONAL :: domain_id
!-
  INTEGER :: iret
  CHARACTER(LEN=70) :: str_t
  INTEGER :: x_id,y_id,z_id,itauid
  INTEGER :: nlonid,nlatid,nlevid,timeid
  INTEGER :: year,month,day,hours,minutes
  REAL :: sec
  CHARACTER(LEN=3),DIMENSION(12) :: &
    cal = (/'JAN','FEB','MAR','APR','MAY','JUN', &
            'JUL','AUG','SEP','OCT','NOV','DEC'/)
  CHARACTER(LEN=30) :: timenow
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) WRITE(*,*) "restopenout 0.0 ",TRIM(fname)
!-
!  If we use the same file for input and output
!- we will not even call restopenout
!-
  iret = NF90_CREATE(fname,NF90_NOCLOBBER,ncfid)
  IF (iret == -35) THEN
    CALL ipslerr (3,'restopenout',&
      & ' The restart file aready exists on the disc. IOIPSL ',&
      & ' will not overwrite it. You should remove the old one or ',&
      & ' generate the new one with another name')
  ENDIF
!-
  iret = NF90_DEF_DIM(ncfid,'x',iim,x_id)
  xax_nb(fid) = xax_nb(fid)+1
  xax_infs(fid,xax_nb(fid),1) = iim
  xax_infs(fid,xax_nb(fid),2) = x_id
!-
  iret = NF90_DEF_DIM(ncfid,'y',jjm,y_id)
  yax_nb(fid)  = yax_nb(fid)+1
  yax_infs(fid,yax_nb(fid),1) = jjm
  yax_infs(fid,yax_nb(fid),2) = y_id
!-
  iret = NF90_DEF_DIM(ncfid,'z',llm,z_id)
  zax_nb(fid) = zax_nb(fid)+1
  zax_infs(fid,zax_nb(fid),1) = llm
  zax_infs(fid,zax_nb(fid),2) = z_id
!-
  iret = NF90_DEF_DIM(ncfid,'time',NF90_UNLIMITED,tdimid_out(fid))
!-
! 1.0 Longitude
!-
  IF (l_dbg) WRITE(*,*) "restopenout 1.0"
!-
  iret = NF90_DEF_VAR(ncfid,"nav_lon",NF90_FLOAT,(/x_id,y_id/),nlonid)
  iret = NF90_PUT_ATT(ncfid,nlonid,'units',"degrees_east")
  iret = NF90_PUT_ATT(ncfid,nlonid,'valid_min',REAL(-180.,KIND=4))
  iret = NF90_PUT_ATT(ncfid,nlonid,'valid_max',REAL( 180.,KIND=4))
  iret = NF90_PUT_ATT(ncfid,nlonid,'long_name',"Longitude")
!-
! 2.0 Latitude
!-
  IF (l_dbg) WRITE(*,*) "restopenout 2.0"
!-
  iret = NF90_DEF_VAR(ncfid,"nav_lat",NF90_FLOAT,(/x_id,y_id/),nlatid)
  iret = NF90_PUT_ATT(ncfid,nlatid,'units',"degrees_north")
  iret = NF90_PUT_ATT(ncfid,nlatid,'valid_min',REAL(-90.,KIND=4))
  iret = NF90_PUT_ATT(ncfid,nlatid,'valid_max',REAL( 90.,KIND=4))
  iret = NF90_PUT_ATT(ncfid,nlatid,'long_name',"Latitude")
!-
! 3.0 Levels
!-
  IF (l_dbg) WRITE(*,*) "restopenout 3.0"
!-
  iret = NF90_DEF_VAR(ncfid,"nav_lev",NF90_FLOAT,z_id,nlevid)
  iret = NF90_PUT_ATT(ncfid,nlevid,'units',"model_levels")
  iret = NF90_PUT_ATT(ncfid,nlevid,'valid_min', &
 &                     REAL(MINVAL(lev),KIND=4))
  iret = NF90_PUT_ATT(ncfid,nlevid,'valid_max', &
 &                     REAL(MAXVAL(lev),KIND=4))
  iret = NF90_PUT_ATT(ncfid,nlevid,'long_name',"Model levels")
!-
! 4.0 Time axis, this is the seconds since axis
!-
  IF (l_dbg) WRITE(*,*) "restopenout 4.0"
!-
  iret = NF90_DEF_VAR(ncfid,"time",NF90_FLOAT, &
                       tdimid_out(fid),timeid)
  tax_varid_out(fid) = timeid
!-
  timeorig(fid) = date
  CALL ju2ymds (date,year,month,day,sec)
  hours   = INT(sec/(60.*60.))
  minutes = INT((sec-hours*60.*60.)/60.)
  sec     = sec-(hours*60.*60.+minutes*60.)
  WRITE (UNIT=str_t, &
   FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
 &  year,month,day,hours,minutes,INT(sec)
  iret = NF90_PUT_ATT(ncfid,timeid,'units',TRIM(str_t))
!-
  CALL ioget_calendar (str_t)
  iret = NF90_PUT_ATT(ncfid,timeid,'calendar',TRIM(str_t))
  iret = NF90_PUT_ATT(ncfid,timeid,'title','Time')
  iret = NF90_PUT_ATT(ncfid,timeid,'long_name','Time axis')
!-
  WRITE(UNIT=str_t, &
   FMT='(" ",I4.4,"-",A3,"-",I2.2," ",I2.2,2(":",I2.2))') &
 &  year,cal(month),day,hours,minutes,INT(sec)
  iret = NF90_PUT_ATT(ncfid,timeid,'time_origin',TRIM(str_t))
!-
! 5.0 Time axis, this is the time steps since axis
!-
  IF (l_dbg) WRITE(*,*) "restopenout 5.0"
!-
  iret = NF90_DEF_VAR(ncfid,"time_steps",NF90_INT, &
 &                    tdimid_out(fid),itauid)
  tind_varid_out(fid) = itauid
!-
  CALL ju2ymds (date,year,month,day,sec)
!-
  hours   = INT(sec/(60.*60.))
  minutes = INT((sec-hours*60.*60.)/60.)
  sec     = sec-(hours*60.*60.+minutes*60.)
!-
  WRITE (UNIT=str_t, &
   FMT='("timesteps since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
 &  year,month,day,hours,minutes,INT(sec)
!-
  iret = NF90_PUT_ATT(ncfid,itauid,'units',TRIM(str_t))
  iret = NF90_PUT_ATT(ncfid,itauid,'title','Time steps')
  iret = NF90_PUT_ATT(ncfid,itauid,'tstep_sec',REAL(timestep,KIND=4))
  iret = NF90_PUT_ATT(ncfid,itauid,'long_name','Time step axis')
!-
  WRITE(UNIT=str_t, &
   FMT='(" ",I4.4,"-",A3,"-",I2.2," ",I2.2,2(":",I2.2))') &
 &  year,cal(month),day,hours,minutes,INT(sec)
  iret = NF90_PUT_ATT(ncfid,itauid,'time_origin',TRIM(str_t))
!-
!  5.2 Write global attributes
!-
  iret = NF90_PUT_ATT(ncfid,NF90_GLOBAL,'Conventions',"CF-1.1")
  iret = NF90_PUT_ATT(ncfid,NF90_GLOBAL,'file_name',TRIM(fname))
!!  TO BE DONE LATER
!!   iret = NF90_PUT_ATT(ncfid,NF90_GLOBAL, &
!!                       'production',TRIM(model_name))
!!   lock_modname = .TRUE.
  CALL ioget_timestamp (timenow)
  iret = NF90_PUT_ATT(ncfid,NF90_GLOBAL,'TimeStamp',TRIM(timenow))
!-
! Add DOMAIN attributes if needed
!-
  CALL flio_dom_att (ncfid,domain_id)
!-
! 6.0 The coordinates are written to the file
!-
  iret = NF90_ENDDEF(ncfid)
!-
  iret = NF90_PUT_VAR(ncfid,nlonid,lon)
  iret = NF90_PUT_VAR(ncfid,nlatid,lat)
  iret = NF90_PUT_VAR(ncfid,nlevid,lev)
!-
! 7.0 Set a few variables related to the out file
!-
  nbvar_out(fid) = 0
  itau_out(fid) = -1
  tstp_out(fid) = 0
  touched_out(fid,:) = .FALSE.
!-
! 7.1 The file is put back in define mode.
!     This will last until itau_out >= 0
!-
  iret = NF90_REDEF(ncfid)
!-
  IF (l_dbg) WRITE(*,*) "restopenout END"
!-------------------------
END SUBROUTINE restopenout
!===
SUBROUTINE restget_opp_r1d &
 & (fid,vname_q,iim,jjm,llm,itau,def_beha, &
 &  var,MY_OPERATOR,nbindex,ijndex)
!---------------------------------------------------------------------
!- This subroutine serves as an interface to restget_real
!-
!- Should work as restput_opp_r1d but the other way around !
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: fid
  CHARACTER(LEN=*) :: vname_q
  INTEGER :: iim,jjm,llm,itau
  LOGICAL def_beha
  REAL :: var(:)
  CHARACTER(LEN=*) :: MY_OPERATOR
  INTEGER :: nbindex,ijndex(nbindex)
!-
  INTEGER :: req_sz,siz1
  REAL :: scal
  CHARACTER(LEN=7) :: topp
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! 0.0 What size should be the data in the file
!-
  req_sz = 1
  IF (nbindex == iim .AND. jjm <= 1 .AND. llm <= 1) THEN
    IF (xax_infs(fid,1,1) > 0) req_sz = req_sz*xax_infs(fid,1,1)
    IF (yax_infs(fid,1,1) > 0) req_sz = req_sz*yax_infs(fid,1,1)
    IF (zax_infs(fid,1,1) > 0) req_sz = req_sz*zax_infs(fid,1,1)
  ELSE
    CALL ipslerr (3,'resget_opp_r1d', &
      'Unable to performe an operation on this variable as it has',&
      'a second and third dimension',vname_q)
  ENDIF
!-
! 1.0 Allocate the temporary buffer we need
!     to put the variable in right dimension
!-
  siz1 = SIZE(var)
  CALL rest_alloc (1,siz1,l_dbg,'restget_opp_r1d')
  CALL rest_alloc (2,req_sz,l_dbg,'restget_opp_r1d')
!-
! 2.0 Here we get the variable from the restart file
!-
  CALL restget_real &
    (fid,vname_q,xax_infs(fid,1,1),yax_infs(fid,1,1), &
     zax_infs(fid,1,1),itau,def_beha,buff_tmp2)
!-
! 4.0 Transfer the buffer obtained from the restart file
!     into the variable the model expects
!-
  topp = MY_OPERATOR(1:MIN(LEN_TRIM(MY_OPERATOR),7))
!-
  IF (INDEX(indchfun,topp(:LEN_TRIM(topp))) > 0) THEN
    scal = missing_val
    CALL mathop (topp,req_sz,buff_tmp2,missing_val, &
 &               nbindex,ijndex,scal,siz1,buff_tmp1)
    var(:) = buff_tmp1(1:siz1)
  ELSE
    CALL ipslerr (3,'resget_opp_r1d', &
      'The operation you wish to do on the variable for the ',&
      'restart file is not allowed.',topp)
  ENDIF
!-----------------------------
END SUBROUTINE restget_opp_r1d
!===
SUBROUTINE restget_opp_r2d &
 & (fid,vname_q,iim,jjm,llm,itau,def_beha, &
 &  var,MY_OPERATOR,nbindex,ijndex)
!---------------------------------------------------------------------
!- This subroutine serves as an interface to restget_real
!-
!- Should work as restput_opp_r2d but the other way around !
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: fid
  CHARACTER(LEN=*) :: vname_q
  INTEGER :: iim,jjm,llm,itau
  LOGICAL def_beha
  REAL :: var(:,:)
  CHARACTER(LEN=*) :: MY_OPERATOR
  INTEGER :: nbindex,ijndex(nbindex)
!-
  INTEGER :: jj,req_sz,ist,var_sz,siz1
  REAL :: scal
  CHARACTER(LEN=7) :: topp
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! 0.0 What size should be the data in the file
!-
  req_sz = 1
  IF (nbindex == iim  .AND. llm <= 1) THEN
    IF (xax_infs(fid,1,1) > 0) req_sz = req_sz*xax_infs(fid,1,1)
    IF (yax_infs(fid,1,1) > 0) req_sz = req_sz*yax_infs(fid,1,1)
  ELSE
    CALL ipslerr (3,'resget_opp_r2d', &
      'Unable to performe an operation on this variable as it has', &
      'a second and third dimension',vname_q)
  ENDIF
!-
  IF (jjm < 1) THEN
    CALL ipslerr (3,'resget_opp_r2d', &
      'Please specify a second dimension which is the', &
      'layer on which the operations are performed',vname_q)
  ENDIF
!-
! 1.0 Allocate the temporary buffer we need
!     to put the variable in right dimension
!-
  siz1 = SIZE(var,1)
  CALL rest_alloc (1,siz1,l_dbg,'restget_opp_r2d')
  CALL rest_alloc (2,req_sz*jjm,l_dbg,'restget_opp_r2d')
!-
! 2.0 Here we get the full variable from the restart file
!-
  CALL restget_real &
 & (fid,vname_q,xax_infs(fid,1,1),yax_infs(fid,1,1), &
 &  jjm,itau,def_beha,buff_tmp2)
!-
! 4.0 Transfer the buffer obtained from the restart file
!     into the variable the model expects
!-
  topp = MY_OPERATOR(1:MIN(LEN_TRIM(MY_OPERATOR),7))
!-
  IF (INDEX(indchfun,topp(:LEN_TRIM(topp))) > 0) THEN
    scal = missing_val
    var_sz = siz1
    DO jj = 1,jjm
      ist = (jj-1)*req_sz+1
      CALL mathop (topp,req_sz,buff_tmp2(ist:ist+req_sz-1), &
 &      missing_val,nbindex,ijndex,scal,var_sz,buff_tmp1)
      var(:,jj) = buff_tmp1(1:siz1)
    ENDDO
  ELSE
    CALL ipslerr (3,'resget_opp_r2d', &
      'The operation you wish to do on the variable for the ',&
      'restart file is not allowed.',topp)
  ENDIF
!-----------------------------
END SUBROUTINE restget_opp_r2d
!===
SUBROUTINE restget_r1d &
 & (fid,vname_q,iim,jjm,llm,itau,def_beha,var)
!---------------------------------------------------------------------
!- This subroutine serves as an interface to restget_real
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: fid
  CHARACTER(LEN=*) :: vname_q
  INTEGER :: iim,jjm,llm,itau
  LOGICAL :: def_beha
  REAL :: var(:)
!-
  INTEGER :: ji,jl,req_sz,var_sz,siz1
  CHARACTER(LEN=70) :: str,str2
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! 1.0 Allocate the temporary buffer we need
!     to put the variable in right dimension
!-
  siz1 = SIZE(var)
  var_sz = siz1
  CALL rest_alloc (1,var_sz,l_dbg,'restget_r1d')
!-
! 2.0 Here we could check if the sizes specified agree
!     with the size of the variable provided
!-
  req_sz = 1
  IF (iim > 0) req_sz = req_sz*iim
  IF (jjm > 0) req_sz = req_sz*jjm
  IF (llm > 0) req_sz = req_sz*llm
  IF (req_sz > var_sz) THEN
    WRITE(str, &
 &    '("Size of variable requested from file should be ",I6)') req_sz
    WRITE(str2, &
 &    '("but the provided variable can only hold ",I6)') var_sz
    CALL ipslerr (3,'restget_r1d',str,str2,' ')
  ENDIF
  IF (req_sz < var_sz) THEN
    WRITE(str, &
 &    '("the size of variable requested from file is ",I6)') req_sz
    WRITE(str2, &
 &    '("but the provided variable can hold ",I6)') var_sz
    CALL ipslerr (2,'restget_r1d', &
      'There could be a problem here :',str,str2)
  ENDIF
!-
  CALL restget_real &
 & (fid,vname_q,iim,jjm,llm,itau,def_beha,buff_tmp1)
!-
! 4.0 Transfer the buffer obtained from the restart file
!     into the variable the model expects
!-
  jl=0
  DO ji=1,siz1
    jl=jl+1
    var(ji) = buff_tmp1(jl)
  ENDDO
!-------------------------
END SUBROUTINE restget_r1d
!===
SUBROUTINE restget_r2d &
 & (fid,vname_q,iim,jjm,llm,itau,def_beha,var)
!---------------------------------------------------------------------
!- This subroutine serves as an interface to restget_real
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: fid
  CHARACTER(LEN=*) :: vname_q
  INTEGER :: iim,jjm,llm,itau
  LOGICAL :: def_beha
  REAL :: var(:,:)
!-
  INTEGER :: ji,jj,jl,req_sz,var_sz,siz1,siz2
  CHARACTER(LEN=70) :: str,str2
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! 1.0 Allocate the temporary buffer we need
!     to put the variable in right dimension
!-
  siz1 = SIZE(var,1)
  siz2 = SIZE(var,2)
  var_sz = siz1*siz2
  CALL rest_alloc (1,var_sz,l_dbg,'restget_r2d')
!-
! 2.0 Here we check if the sizes specified agree
!     with the size of the variable provided
!-
  req_sz = 1
  IF (iim > 0) req_sz = req_sz*iim
  IF (jjm > 0) req_sz = req_sz*jjm
  IF (llm > 0) req_sz = req_sz*llm
  IF (req_sz > var_sz) THEN
    WRITE(str, &
 &    '("Size of variable ",A, &
 &      //" requested from file should be ",I6)') TRIM(vname_q),req_sz
    WRITE(str2, &
 &    '("but the provided variable can only hold ",I6)') var_sz
    CALL ipslerr (3,'restget_r2d',str,str2,' ')
  ENDIF
  IF (req_sz < var_sz) THEN
    WRITE(str, &
 &    '("Size of variable ",A, &
 &      //" requested from file is ",I6)') TRIM(vname_q),req_sz
    WRITE(str2,'("but the provided variable can hold ",I6)') var_sz
    CALL ipslerr (2,'restget_r2d', &
      'There could be a problem here :',str,str2)
  ENDIF
!-
  CALL restget_real &
 & (fid,vname_q,iim,jjm,llm,itau,def_beha,buff_tmp1)
!-
! 4.0 Transfer the buffer obtained from the restart file
!     into the variable the model expects
!-
  jl=0
  DO jj=1,siz2
    DO ji=1,siz1
      jl=jl+1
      var(ji,jj) = buff_tmp1(jl)
    ENDDO
  ENDDO
!-------------------------
END SUBROUTINE restget_r2d
!===
SUBROUTINE restget_r3d &
  (fid,vname_q,iim,jjm,llm,itau,def_beha,var)
!---------------------------------------------------------------------
!- This subroutine serves as an interface to restget_real
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: fid
  CHARACTER(LEN=*) :: vname_q
  INTEGER :: iim,jjm,llm,itau
  LOGICAL def_beha
  REAL :: var(:,:,:)
!-
  INTEGER :: ji,jj,jk,jl,req_sz,var_sz,siz1,siz2,siz3
  CHARACTER(LEN=70) :: str,str2
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! 1.0 Allocate the temporary buffer we need
!     to put the variable in right dimension
!-
  siz1 = SIZE(var,1)
  siz2 = SIZE(var,2)
  siz3 = SIZE(var,3)
  var_sz = siz1*siz2*siz3
  CALL rest_alloc (1,var_sz,l_dbg,'restget_r3d')
!-
! 2.0 Here we check if the sizes specified agree
!     with the size of the variable provided
!-
  req_sz = 1
  IF (iim > 0) req_sz = req_sz*iim
  IF (jjm > 0) req_sz = req_sz*jjm
  IF (llm > 0) req_sz = req_sz*llm
  IF (req_sz > var_sz) THEN
    WRITE(str, &
 &    '("Size of variable ",A, &
 &      //" requested from file should be ",I6)') TRIM(vname_q),req_sz
    WRITE(str2, &
 &    '("but the provided variable can only hold ",I6)') var_sz
    CALL ipslerr (3,'restget_r3d',str,str2,' ')
  ENDIF
  IF (req_sz < var_sz) THEN
    WRITE(str, &
 &    '("Size of variable ",A, &
 &      //" requested from file is ",I6)') TRIM(vname_q),req_sz
    WRITE(str2,'("but the provided variable can hold ",I6)') var_sz
    CALL ipslerr (2,'restget_r3d', &
      'There could be a problem here :',str,str2)
  ENDIF
!-
  CALL restget_real &
    (fid,vname_q,iim,jjm,llm,itau,def_beha,buff_tmp1)
!-
! 4.0 Transfer the buffer obtained from the restart file
!     into the variable the model expects
!-
  jl=0
  DO jk=1,siz3
    DO jj=1,siz2
      DO ji=1,siz1
        jl=jl+1
        var(ji,jj,jk) = buff_tmp1(jl)
      ENDDO
    ENDDO
  ENDDO
!-------------------------
END SUBROUTINE restget_r3d
!===
SUBROUTINE restget_real &
  (fid,vname_q,iim,jjm,llm,itau,def_beha,var)
!---------------------------------------------------------------------
!- This subroutine is for getting a variable from the restart file.
!- A number of verifications will be made :
!- - Is this the first time we read this variable ?
!- - Are the dimensions correct ?
!- - Is the correct time step present in the file
!- - is a default behaviour possible. If not the model is stoped.
!- Default procedure is to write the content of val_exp on all values.
!-
!- INPUT
!-
!- fid            : Identification of the file
!- vname_q        : Name of the variable to be read
!- iim, jjm ,llm  : Dimensions of the variable that should be read
!- itau           : Time step at whcih we are when we want
!-                  to read the variable
!- def_beha       : If the model can restart without this variable
!-                  then some strange value is given.
!-
!- OUTPUT
!-
!- var            : Variable in which the data is put
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: fid
  CHARACTER(LEN=*) :: vname_q
  INTEGER :: iim,jjm,llm,itau
  LOGICAL :: def_beha
  REAL :: var(:)
!-
  INTEGER :: vid,vnb,ncfid,iret,index,it,ndim,ia
  CHARACTER(LEN=70) str,str2
  CHARACTER(LEN=80) attname
  INTEGER,DIMENSION(4) :: corner,edge
!---------------------------------------------------------------------
  ncfid = netcdf_id(fid,1)
!-
  CALL find_str (varname_in(fid,1:nbvar_in(fid)),vname_q,vnb)
!-
! 1.0 If the variable is not present then ERROR or filled up
!     by default values if allowed
!-
  IF (vnb < 0) THEN
    IF (def_beha) THEN
!-----
      lock_valexp = .TRUE.
      var(:) = val_exp
!----
      str = 'Variable '//TRIM(vname_q) &
          //' is not present in the restart file'
      CALL ipslerr (1,'restget', &
 &      str,'but default values are used to fill in',' ')
!----
      IF (nbvar_in(fid) >= max_var) THEN
        CALL ipslerr (3,'restget', &
         'Too many variables for the restcom module', &
         'Please increase the value of max_var',' ')
      ENDIF
      nbvar_in(fid) = nbvar_in(fid)+1
      vnb = nbvar_in(fid)
      varname_in(fid,vnb) = vname_q
      touched_in(fid,vnb) = .TRUE.
!-----
      CALL restdefv (fid,vname_q,iim,jjm,llm,.TRUE.)
!-----
    ELSE
      str = 'Variable '//TRIM(vname_q) &
          //' is not present in the restart file'
      CALL ipslerr (3,'restget', &
 &      str,'but it is need to restart the model',' ')
    ENDIF
!---
  ELSE
!---
!--  2.0 Check if the variable has not yet been read
!--      and that the time is OK
!---
    vid = varid_in(fid,vnb)
!---
    nbvar_read(fid) = nbvar_read(fid)+1
!---
    IF (touched_in(fid,vnb)) THEN
      str = 'Variable '//TRIM(vname_q) &
          //' has already been read from file'
      CALL ipslerr (3,'restget',str,' ',' ')
    ENDIF
!---
!-- 3.0 get the time step of the restart file
!--     and check if it is correct
!---
    index = -1
    DO it=1,tax_size_in(fid)
      IF (t_index(fid,it) == itau)  index = it
    ENDDO
    IF (index < 0) THEN
      str = 'The time step requested for variable '//TRIM(vname_q)
      CALL ipslerr (3,'restget', &
 &      str,'is not available in the current file',' ')
    ENDIF
!---
!-- 4.0 Read the data. Note that the variables in the restart files
!--     have no time axis is and thus we write -1
!---
    str='Incorrect dimension for '//TRIM(vname_q)
    ndim = 0
    IF (iim > 0) THEN
      ndim = ndim+1
      IF (vardims_in(fid,vnb,ndim) == iim) THEN
        corner(ndim) = 1
        edge(ndim) = iim
      ELSE
        WRITE (str2,'("Incompatibility for iim : ",I6,I6)') &
             iim,vardims_in(fid,vnb,ndim)
        CALL ipslerr (3,'restget',str,str2,' ')
      ENDIF
    ENDIF
!---
    IF (jjm > 0) THEN
      ndim = ndim+1
      IF (vardims_in(fid,vnb,ndim) == jjm) THEN
        corner(ndim) = 1
        edge(ndim) = jjm
      ELSE
        WRITE (str2,'("Incompatibility for jjm : ",I6,I6)') &
             jjm,vardims_in(fid,vnb,ndim)
        CALL ipslerr (3,'restget',str,str2,' ')
      ENDIF
    ENDIF
!---
    IF (llm > 0) THEN
      ndim = ndim+1
      IF (vardims_in(fid,vnb,ndim) == llm) THEN
        corner(ndim) = 1
        edge(ndim) = llm
      ELSE
        WRITE (str2,'("Incompatibility for llm : ",I6,I6)') &
             llm,vardims_in(fid,vnb,ndim)
        CALL ipslerr (3,'restget',str,str2,' ')
      ENDIF
    ENDIF
!---
!-- Time
!---
    ndim = ndim+1
    corner(ndim) = index
!!????? edge(ndim) = index
    edge(ndim) = 1
!---
    iret = NF90_GET_VAR(ncfid,vid,var, &
 &                      start=corner(1:ndim),count=edge(1:ndim))
!---
!-- 5.0 The variable we have just read is created
!--      in the next restart file
!---
    IF (     (netcdf_id(fid,1) /= netcdf_id(fid,2))  &
 &      .AND.(netcdf_id(fid,2) > 0) ) THEN
!-----
      CALL restdefv (fid,vname_q,iim,jjm,llm,.FALSE.)
!-----
      DO ia = 1,varatt_in(fid,vnb)
        iret = NF90_INQ_ATTNAME(ncfid,vid,ia,attname)
        iret = NF90_COPY_ATT(ncfid,vid,attname, &
 &               netcdf_id(fid,2),varid_out(fid,nbvar_out(fid)))
      ENDDO
!-----
      IF (itau_out(fid) >= 0) THEN
        iret = NF90_ENDDEF(netcdf_id(fid,2))
      ENDIF
    ENDIF
!---
  ENDIF
!--------------------------
END SUBROUTINE restget_real
!===
SUBROUTINE restput_opp_r1d &
 & (fid,vname_q,iim,jjm,llm,itau,var,MY_OPERATOR,nbindex,ijndex)
!---------------------------------------------------------------------
!- This subroutine is the interface to restput_real which allows
!- to re-index data onto the original grid of the restart file.
!- The logic we use is still fuzzy in my mind but that is probably
!- only because I have not yet though through everything.
!-
!- In the case iim = nbindex it means that the user attempts
!- to project a vector back onto the original 2D or 3D field.
!- This requires that jjm and llm be equal to 1 or 0,
!- else I would not know what it means.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: fid
  CHARACTER(LEN=*) :: vname_q
  INTEGER :: iim,jjm,llm,itau
  REAL :: var(:)
  CHARACTER(LEN=*) :: MY_OPERATOR
  INTEGER :: nbindex,ijndex(nbindex)
!-
  INTEGER :: req_sz,siz1
  REAL :: scal
  CHARACTER(LEN=7) :: topp
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! 0.0 What size should be the data in the file
!-
  req_sz = 1
  IF ( nbindex == iim .AND. jjm <= 1 .AND. llm <= 1) THEN
    IF (xax_infs(fid,1,1) > 0) req_sz = req_sz*xax_infs(fid,1,1)
    IF (yax_infs(fid,1,1) > 0) req_sz = req_sz*yax_infs(fid,1,1)
    IF (zax_infs(fid,1,1) > 0) req_sz = req_sz*zax_infs(fid,1,1)
  ELSE
    CALL ipslerr (3,'restput_opp_r1d', &
      'Unable to performe an operation on this variable as it has', &
      'a second and third dimension',vname_q)
  ENDIF
!-
! 1.0 Allocate the temporary buffer we need
!     to put the variable in right dimension
!-
  siz1 = SIZE(var)
  CALL rest_alloc (1,siz1,l_dbg,'restput_opp_r1d')
  CALL rest_alloc (2,req_sz,l_dbg,'restput_opp_r1d')
!-
! 2.0 We do the operation needed.
!     It can only be a re-indexing operation.
!     You would not want to change the values in a restart file or ?
!-
  topp = MY_OPERATOR(1:MIN(LEN_TRIM(MY_OPERATOR),7))
!-
  IF (INDEX(indchfun,topp(:LEN_TRIM(topp))) > 0) THEN
    scal = missing_val
    buff_tmp1(1:siz1) = var(:)
    CALL mathop &
 &    (topp,siz1,buff_tmp1,missing_val,nbindex,ijndex, &
 &     scal,req_sz,buff_tmp2)
  ELSE
    CALL ipslerr (3,'restput_opp_r1d', &
 &    'The operation you wish to do on the variable for the ', &
 &    'restart file is not allowed.',topp)
  ENDIF
!-
  CALL restput_real &
 & (fid,vname_q,xax_infs(fid,1,1),yax_infs(fid,1,1), &
 &  zax_infs(fid,1,1),itau,buff_tmp2)
!-----------------------------
END SUBROUTINE restput_opp_r1d
!===
SUBROUTINE restput_opp_r2d &
 & (fid,vname_q,iim,jjm,llm,itau,var,MY_OPERATOR,nbindex,ijndex)
!---------------------------------------------------------------------
!- This subroutine is the interface to restput_real which allows
!- to re-index data onto the original grid of the restart file.
!- The logic we use is still fuzzy in my mind but that is probably
!- only because I have not yet though through everything.
!-
!- In the case iim = nbindex it means that the user attempts
!- to project the first dimension of the matrix back onto a 3D field
!- where jjm will be the third dimension.
!- Here we do not allow for 4D data, thus we will take the first
!- two dimensions in the file and require that llm = 1.
!- These are pretty heavy constraints but I do not know how
!- to make it more general. I need to think about it some more.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: fid
  CHARACTER(LEN=*) :: vname_q
  INTEGER :: iim,jjm,llm,itau
  REAL :: var(:,:)
  CHARACTER(LEN=*) :: MY_OPERATOR
  INTEGER :: nbindex,ijndex(nbindex)
!-
  INTEGER :: jj,req_sz,ist,siz1
  REAL :: scal
  CHARACTER(LEN=7) :: topp
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! 0.0 What size should be the data in the file
!-
  req_sz = 1
  IF ( nbindex == iim .AND. llm <= 1) THEN
    IF (xax_infs(fid,1,1) > 0) req_sz = req_sz*xax_infs(fid,1,1)
    IF (yax_infs(fid,1,1) > 0) req_sz = req_sz*yax_infs(fid,1,1)
  ELSE
    CALL ipslerr (3,'restput_opp_r2d', &
      'Unable to performe an operation on this variable as it has', &
      'a second and third dimension',vname_q)
  ENDIF
!-
  IF (jjm < 1) THEN
    CALL ipslerr (3,'restput_opp_r2d', &
      'Please specify a second dimension which is the', &
      'layer on which the operations are performed',vname_q)
  ENDIF
!-
! 1.0 Allocate the temporary buffer we need
!     to put the variable in right dimension
!-
  siz1 = SIZE(var,1)
  CALL rest_alloc (1,siz1,l_dbg,'restput_opp_r2d')
  CALL rest_alloc (2,req_sz*jjm,l_dbg,'restput_opp_r2d')
!-
! 2.0 We do the operation needed.
!     It can only be a re-indexing operation.
!     You would not want to change the values in a restart file or ?
!-
  topp = MY_OPERATOR(1:MIN(LEN_TRIM(MY_OPERATOR),7))
!-
  IF (INDEX(indchfun,topp(:LEN_TRIM(topp))) > 0) THEN
    scal = missing_val
    DO jj = 1,jjm
      buff_tmp1(1:siz1) = var(:,jj)
      ist = (jj-1)*req_sz+1
      CALL mathop &
 &      (topp,siz1,buff_tmp1,missing_val,nbindex,ijndex, &
 &       scal,req_sz,buff_tmp2(ist:ist+req_sz-1))
    ENDDO
  ELSE
    CALL ipslerr (3,'restput_opp_r2d', &
 &    'The operation you wish to do on the variable for the ', &
 &    'restart file is not allowed.',topp)
  ENDIF
!-
  CALL restput_real &
 & (fid,vname_q,xax_infs(fid,1,1),yax_infs(fid,1,1), &
 &  jjm,itau,buff_tmp2)
!-----------------------------
END SUBROUTINE restput_opp_r2d
!===
SUBROUTINE restput_r1d (fid,vname_q,iim,jjm,llm,itau,var)
!---------------------------------------------------------------------
!- This subroutine serves as an interface to restput_real
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: fid
  CHARACTER(LEN=*) :: vname_q
  INTEGER :: iim,jjm,llm,itau
  REAL :: var(:)
!-
  INTEGER :: ji,jl,req_sz,var_sz,siz1
  CHARACTER(LEN=70) :: str,str2
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! 1.0 Allocate the temporary buffer we need
!     to put the variable in right dimension
!-
  siz1 = SIZE(var)
  var_sz = siz1
  CALL rest_alloc (1,var_sz,l_dbg,'restput_r1d')
!-
! 2.0 Here we could check if the sizes specified agree
!     with the size of the variable provided
!-
  req_sz = 1
  IF (iim > 0) req_sz = req_sz*iim
  IF (jjm > 0) req_sz = req_sz*jjm
  IF (llm > 0) req_sz = req_sz*llm
  IF (req_sz > var_sz) THEN
    WRITE(str, &
 &    '("Size of variable put to the file should be ",I6)') req_sz
    WRITE(str2, &
 &    '("but the provided variable is of size ",I6)') var_sz
    CALL ipslerr (3,'restput_r1d',str,str2,' ')
  ENDIF
  IF (req_sz < var_sz) THEN
    WRITE(str,'("the size of variable put to the file is ",I6)') req_sz
    WRITE(str2,'("but the provided variable is larger ",I6)') var_sz
    CALL ipslerr (2,'restput_r1d', &
      'There could be a problem here :',str,str2)
  ENDIF
!-
! 4.0 Transfer the buffer obtained from the restart file
!     into the variable the model expects
!-
  jl=0
  DO ji=1,siz1
    jl=jl+1
    buff_tmp1(jl) = var(ji)
  ENDDO
!-
  CALL restput_real (fid,vname_q,iim,jjm,llm,itau,buff_tmp1)
!-------------------------
END SUBROUTINE restput_r1d
!===
SUBROUTINE restput_r2d (fid,vname_q,iim,jjm,llm,itau,var)
!---------------------------------------------------------------------
!- This subroutine serves as an interface to restput_real
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: fid
  CHARACTER(LEN=*) :: vname_q
  INTEGER :: iim,jjm,llm,itau
  REAL :: var(:,:)
!-
  INTEGER :: ji,jj,jl,req_sz,var_sz,siz1,siz2
  CHARACTER(LEN=70) :: str,str2
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! 1.0 Allocate the temporary buffer we need
!     to put the variable in right dimension
!-
  siz1 = SIZE(var,1)
  siz2 = SIZE(var,2)
  var_sz = siz1*siz2
  CALL rest_alloc (1,var_sz,l_dbg,'restput_r2d')
!-
! 2.0 Here we could check if the sizes specified agree
!     with the size of the variable provided
!-
  req_sz = 1
  IF (iim > 0) req_sz = req_sz*iim
  IF (jjm > 0) req_sz = req_sz*jjm
  IF (llm > 0) req_sz = req_sz*llm
  IF (req_sz > var_sz) THEN
    WRITE(str, &
&         '("Size of variable put to the file should be ",I6)') req_sz
    WRITE(str2,'("but the provided variable is of size  ",I6)') var_sz
    CALL ipslerr (3,'restput_r2d',str,str2,' ')
  ENDIF
  IF (req_sz < var_sz) THEN
    WRITE(str,'("the size of variable put to the file is ",I6)') req_sz
    WRITE(str2,'("but the provided variable is larger ",I6)')  var_sz
    CALL ipslerr (2,'restput_r2d', &
      'There could be a problem here :',str,str2)
  ENDIF
!-
! 4.0 Transfer the buffer obtained from the restart file
!     into the variable the model expects
!-
  jl=0
  DO jj=1,siz2
    DO ji=1,siz1
      jl=jl+1
      buff_tmp1(jl) = var(ji,jj)
    ENDDO
  ENDDO
!-
  CALL restput_real(fid,vname_q,iim,jjm,llm,itau,buff_tmp1)
!-------------------------
END SUBROUTINE restput_r2d
!===
SUBROUTINE restput_r3d (fid,vname_q,iim,jjm,llm,itau,var)
!---------------------------------------------------------------------
!- This subroutine serves as an interface to restput_real
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: fid
  CHARACTER(LEN=*) :: vname_q
  INTEGER :: iim,jjm,llm,itau
  REAL :: var(:,:,:)
!-
  INTEGER :: ji,jj,jk,jl,req_sz,var_sz,siz1,siz2,siz3
  CHARACTER(LEN=70) :: str,str2
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! 1.0 Allocate the temporary buffer we need
!     to put the variable in right dimension
!-
  siz1 = SIZE(var,1)
  siz2 = SIZE(var,2)
  siz3 = SIZE(var,3)
  var_sz = siz1*siz2*siz3
  CALL rest_alloc (1,var_sz,l_dbg,'restput_r3d')
!-
! 2.0 Here we could check if the sizes specified agree
!     with the size of the variable provided
!-
  req_sz = 1
  IF (iim > 0) req_sz = req_sz*iim
  IF (jjm > 0) req_sz = req_sz*jjm
  IF (llm > 0) req_sz = req_sz*llm
  IF (req_sz > var_sz) THEN
    WRITE(str, &
 &    '("Size of variable put to the file should be ",I6)') req_sz
    WRITE(str2, &
 &    '("but the provided variable is of size  ",I6)')  var_sz
    CALL ipslerr (3,'restput_r3d',str,str2,' ')
  ENDIF
  IF (req_sz < var_sz) THEN
    WRITE(str,'("the size of variable put to the file is ",I6)') req_sz
    WRITE(str2,'("but the provided variable is larger ",I6)')  var_sz
    CALL ipslerr (2,'restput_r3d', &
      'There could be a problem here :',str,str2)
  ENDIF
!-
! 4.0 Transfer the buffer obtained from the restart file
!     into the variable the model expects
!-
  jl=0
  DO jk=1,siz3
    DO jj=1,siz2
      DO ji=1,siz1
        jl=jl+1
        buff_tmp1(jl) = var(ji,jj,jk)
      ENDDO
    ENDDO
  ENDDO
!-
  CALL restput_real(fid,vname_q,iim,jjm,llm,itau,buff_tmp1)
!-------------------------
END SUBROUTINE restput_r3d
!===
SUBROUTINE restput_real (fid,vname_q,iim,jjm,llm,itau,var)
!---------------------------------------------------------------------
!- This subroutine will put a variable into the restart file.
!- But it will do a lot of other things if needed :
!- - Open a file if non is opened for this time-step
!-   and all variables were written.
!- - Add an axis if needed
!- - verify that the variable has the right time step for this file
!- - If it is time for a new file then it is opened
!-   and the old one closed
!- This requires that variables read from the last restart file were all
!- written
!-
!- INPUT
!-
!- fid         : Id of the file in which we will write the variable
!- vname_q     : Name of the variable to be written
!- iim,jjm,llm : Size in 3D of the variable
!- itau        : Time step at which the variable is written
!- var         : Variable
!-
!- OUTPUT
!-
!- NONE
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: vname_q
  INTEGER :: fid,iim,jjm,llm,itau
  REAL :: var(:)
!-
  INTEGER :: iret,vid,ncid,iv,vnb
  INTEGER :: ierr
  REAL :: secsince,one_day,one_year
  INTEGER :: ndims
  INTEGER,DIMENSION(4) :: corner,edge
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
! 0.0 Get some variables
!-
  ncid = netcdf_id(fid,2)
  IF (netcdf_id(fid,2) < 0) THEN
    CALL ipslerr (3,'restput', &
 &    'The output restart file is undefined.',' ',' ')
  ENDIF
  CALL ioget_calendar (one_year,one_day)
!-
! 1.0 Check if the variable is already present
!-
  IF (l_dbg) WRITE(*,*) 'RESTPUT 1.0 : ',TRIM(vname_q)
!-
  CALL find_str (varname_out(fid,1:nbvar_out(fid)),vname_q,vnb)
!-
  IF (l_dbg) THEN
    WRITE(*,*) 'RESTPUT 1.1 : ',varname_out(fid,1:nbvar_out(fid)),vnb
  ENDIF
!-
! 2.0 If variable is not present then declare it
!     and add extra dimensions if needed.
!-
  IF (vnb <= 0) THEN
    CALL restdefv (fid,vname_q,iim,jjm,llm,.TRUE.)
    vnb = nbvar_out(fid)
  ENDIF
  vid = varid_out(fid,vnb)
!-
  IF (l_dbg) WRITE(*,*) 'RESTPUT 2.0 : ',vnb,vid
!-
! 2.1 Is this file already in write mode ?
!     If itau_out is still negative then we have
!     never written to it and we need to go into write mode.
!-
  IF (itau_out(fid) < 0) THEN
    iret = NF90_ENDDEF(ncid)
  ENDIF
!-
! 3.0 Is this itau already on the axis ?
!     If not then check that all variables of previous time is OK.
!-
  IF (l_dbg) WRITE(*,*) 'RESTPUT 3.0 : ',itau,itau_out(fid)
!-
  IF (itau /= itau_out(fid)) THEN
!---
!-- If it is the first time step written on the restart
!-- then we only check the number
!-- Else we see if every variable was written
!---
    IF (tstp_out(fid) == 0) THEN
      IF (nbvar_out(fid) < nbvar_read(fid)) THEN
        WRITE(*,*) "ERROR :",tstp_out(fid), &
                   nbvar_out(fid),nbvar_read(fid)
        CALL ipslerr (1,'restput', &
 &        'There are fewer variables read from the output file', &
 &        'than written onto the input file.', &
 &        'We trust you know what you are doing')
      ENDIF
    ELSE
      ierr = 0
      DO iv=1,nbvar_out(fid)
        IF (.NOT.touched_out(fid,iv)) ierr = ierr+1
      ENDDO
      IF (ierr > 0) THEN
        WRITE(*,*) "ERROR :",nbvar_out(fid)
        CALL ipslerr (1,'restput', &
 &        'There are fewer variables in the output file for this', &
 &        'time step than for the previous one',' ')
      ELSE
        touched_out(fid,:) = .FALSE.
      ENDIF
    ENDIF
!---
    secsince = itau*deltat(fid)
    corner(1) =  tstp_out(fid)+1
    edge(1) = 1
!---
!-- 3.1 Here we add the values to the time axes
!---
    IF (l_dbg) THEN
      WRITE(*,*) 'RESTPUT 3.1 : ',itau,secsince,corner(1),edge(1)
    ENDIF
!---
    iret = NF90_PUT_VAR(ncid,tind_varid_out(fid),itau, &
 &                      start=corner(1:1))
    iret = NF90_PUT_VAR(ncid,tax_varid_out(fid),secsince, &
 &                      start=corner(1:1))
!---
    tstp_out(fid) = tstp_out(fid)+1
    itau_out(fid) = itau
  ENDIF
!-
! 4.0 Variable and time step should be present
!     now so we can dump variable
!-
  ndims = 0
  IF (iim > 0) THEN
    ndims = ndims+1
    corner(ndims) = 1
    edge(ndims) = iim
  ENDIF
  IF (jjm > 0) THEN
    ndims = ndims+1
    corner(ndims) = 1
    edge(ndims) = jjm
  ENDIF
  IF (llm > 0) THEN
    ndims = ndims+1
    corner(ndims) = 1
    edge(ndims) = llm
  ENDIF
  ndims = ndims+1
  corner(ndims) = tstp_out(fid)
  edge(ndims) = 1
!-
  iret = NF90_PUT_VAR(ncid,vid,var, &
 &                    start=corner(1:ndims),count=edge(1:ndims))
!-
  IF (iret /= NF90_NOERR) THEN
    CALL ipslerr (2,'restput_real',NF90_STRERROR(iret), &
 &    'Bug in restput.',&
 &    'Please, verify compatibility between get and put commands.')
  ENDIF
!-
!  5.0 Note that the variables was treated
!-
  touched_out(fid,vnb) = .TRUE.
!---------------------------
END  SUBROUTINE restput_real
!===
SUBROUTINE restdefv (fid,varname,iim,jjm,llm,write_att)
!---------------------------------------------------------------------
! This subroutine adds a variable to the output file.
! The attributes are either taken from.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER ::fid
  CHARACTER(LEN=*) :: varname
  INTEGER :: iim,jjm,llm
  LOGICAL :: write_att
!-
  INTEGER :: dims(4),ic,xloc,ndim,ncfid
  INTEGER :: iret,ax_id
  CHARACTER(LEN=3) :: str
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  ncfid = netcdf_id(fid,2)
  IF (nbvar_out(fid) >= max_var) THEN
    CALL ipslerr (3,'restdefv', &
      'Too many variables for the restcom module', &
      'Please increase the value of max_var',' ')
  ENDIF
  nbvar_out(fid) = nbvar_out(fid)+1
  varname_out(fid,nbvar_out(fid)) = varname
!-
! 0.0 Put the file in define mode if needed
!-
  IF (itau_out(fid) >= 0) THEN
    iret = NF90_REDEF(ncfid)
  ENDIF
!-
! 1.0 Do we have all dimensions and can we go ahead
!-
  IF (l_dbg) THEN
    WRITE(*,*) 'restdefv 1.0 :',TRIM(varname),nbvar_out(fid)
  ENDIF
!-
  ndim = 0
!-
! 1.1 Work on x
!-
  IF (iim > 0) THEN
    ndim = ndim+1
    xloc = 0
    DO ic=1,xax_nb(fid)
      IF (xax_infs(fid,ic,1) == iim) xloc = ic
    ENDDO
!---
    IF (xloc > 0) THEN
      dims(ndim) = xax_infs(fid,xloc,2)
    ELSE
      str='x_'//CHAR(96+xax_nb(fid))
      iret = NF90_DEF_DIM(ncfid,str,iim,ax_id)
      xax_nb(fid) = xax_nb(fid)+1
      xax_infs(fid,xax_nb(fid),1) = iim
      xax_infs(fid,xax_nb(fid),2) = ax_id
      dims(ndim) = ax_id
    ENDIF
  ENDIF
!-
! 1.2 Work on y
!-
  IF (jjm > 0) THEN
    ndim = ndim+1
    xloc = 0
    DO ic=1,yax_nb(fid)
      IF (yax_infs(fid,ic,1) == jjm) xloc = ic
    ENDDO
!---
    IF (xloc > 0) THEN
      dims(ndim) = yax_infs(fid,xloc,2)
    ELSE
      str='y_'//CHAR(96+yax_nb(fid))
      iret = NF90_DEF_DIM(ncfid,str,jjm,ax_id)
      yax_nb(fid) = yax_nb(fid)+1
      yax_infs(fid,yax_nb(fid),1) = jjm
      yax_infs(fid,yax_nb(fid),2) = ax_id
      dims(ndim) = ax_id
    ENDIF
  ENDIF
!-
! 1.3 Work on z
!-
  IF (llm > 0) THEN
    ndim = ndim+1
    xloc = 0
    DO ic=1,zax_nb(fid)
      IF (zax_infs(fid,ic,1) == llm) xloc = ic
    ENDDO
!---
    IF (xloc > 0) THEN
      dims(ndim) = zax_infs(fid,xloc,2)
    ELSE
      str='z_'//CHAR(96+zax_nb(fid))
      iret = NF90_DEF_DIM(ncfid,str,llm,ax_id)
      zax_nb(fid) = zax_nb(fid)+1
      zax_infs(fid,zax_nb(fid),1) = llm
      zax_infs(fid,zax_nb(fid),2) = ax_id
      dims(ndim) = ax_id
    ENDIF
  ENDIF
!-
! 1.4  Time needs to be added
!-
  ndim = ndim+1
  dims(ndim) = tdimid_out(fid)
!-
! 2.0  Declare the variable
!-
  IF (l_dbg) THEN
    WRITE(*,*) 'restdefv 2.0 :',ndim,' :: ',dims(1:ndim),tdimid_out(fid)
  ENDIF
!-
  iret = NF90_DEF_VAR(ncfid,varname,NF90_DOUBLE,dims(1:ndim), &
 &                    varid_out(fid,nbvar_out(fid)))
  IF (iret /= NF90_NOERR) THEN
    CALL ipslerr (3,'restdefv', &
      'Could not define new variable in file', &
      NF90_STRERROR(iret),varname)
  ENDIF
!-
! 3.0 Add the attributes if requested
!-
  IF (write_att) THEN
    IF (rest_units /= 'XXXXX') THEN
      iret =  NF90_PUT_ATT(ncfid,varid_out(fid,nbvar_out(fid)), &
 &                         'units',TRIM(rest_units))
      rest_units = 'XXXXX'
    ENDIF
!---
    IF (rest_lname /= 'XXXXX') THEN
      iret =  NF90_PUT_ATT(ncfid,varid_out(fid,nbvar_out(fid)), &
 &                         'long_name',TRIM(rest_lname))
      rest_lname = 'XXXXX'
    ENDIF
!---
    iret = NF90_PUT_ATT(ncfid,varid_out(fid,nbvar_out(fid)), &
 &                      'missing_value',REAL(missing_val,KIND=4))
!---
    IF (itau_out(fid) >= 0) THEN
      iret = NF90_ENDDEF(ncfid)
    ENDIF
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) &
 &    'restdefv 3.0 : LIST OF VARS ',varname_out(fid,1:nbvar_out(fid))
  ENDIF
!----------------------
END SUBROUTINE restdefv
!===
SUBROUTINE rest_atim (l_msg,c_p)
!---------------------------------------------------------------------
! Called by "c_p", [re]allocate the time axes
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  LOGICAL,INTENT(IN) :: l_msg
  CHARACTER(LEN=*),INTENT(IN) :: c_p
!-
  INTEGER :: i_err,tszij
  INTEGER,ALLOCATABLE :: tmp_index(:,:)
  REAL,ALLOCATABLE :: tmp_julian(:,:)
!---------------------------------------------------------------------
!-
!  Allocate the space we need for the time axes
!-
  IF (.NOT.ALLOCATED(t_index) .AND. .NOT.ALLOCATED(t_julian)) THEN
    IF (l_msg) THEN
      WRITE(*,*) TRIM(c_p)//' : Allocate times axes at :', &
 &               max_file,tax_size_in(nb_fi)
    ENDIF
!---
    ALLOCATE(t_index(max_file,tax_size_in(nb_fi)),STAT=i_err)
    IF (i_err/=0) THEN
      WRITE(*,*) "ERROR IN ALLOCATION of t_index : ",i_err
      CALL ipslerr (3,TRIM(c_p), &
 &      'Problem in allocation of t_index','', &
 &      '(you must increase memory)')
    ENDIF
    t_index (:,:) = 0
!---
    ALLOCATE(t_julian(max_file,tax_size_in(nb_fi)),STAT=i_err)
    IF (i_err/=0) THEN
      WRITE(*,*) "ERROR IN ALLOCATION of t_julian : ",i_err
      CALL ipslerr (3,TRIM(c_p), &
 &      'Problem in allocation of max_file,tax_size_in','', &
 &      '(you must increase memory)')
    ENDIF
    t_julian (:,:) = 0.0
  ELSE IF (    (SIZE(t_index,DIM=2)  < tax_size_in(nb_fi)) &
 &         .OR.(SIZE(t_julian,DIM=2) < tax_size_in(nb_fi)) ) THEN
    IF (l_msg) THEN
      WRITE(*,*) TRIM(c_p)//' : Reallocate times axes at :', &
 &               max_file,tax_size_in(nb_fi)
    ENDIF
!---
    ALLOCATE (tmp_index(max_file,tax_size_in(nb_fi)),STAT=i_err)
    IF (i_err/=0) THEN
      WRITE(*,*) "ERROR IN ALLOCATION of tmp_index : ",i_err
      CALL ipslerr (3,TRIM(c_p), &
 &      'Problem in allocation of tmp_index','', &
 &      '(you must increase memory)')
    ENDIF
    tszij = SIZE(t_index,DIM=2)
    tmp_index(:,1:tszij) = t_index(:,1:tszij)
    DEALLOCATE(t_index)
    ALLOCATE (t_index(max_file,tax_size_in(nb_fi)),STAT=i_err)
    IF (i_err/=0) THEN
      WRITE(*,*) "ERROR IN ALLOCATION of t_index : ",i_err
      CALL ipslerr (3,TRIM(c_p), &
 &     'Problem in reallocation of t_index','', &
 &     '(you must increase memory)')
    ENDIF
    t_index(:,1:tszij) = tmp_index(:,1:tszij)
!---
    ALLOCATE (tmp_julian(max_file,tax_size_in(nb_fi)),STAT=i_err)
    IF (i_err/=0) THEN
      WRITE(*,*) "ERROR IN ALLOCATION of tmp_julian : ",i_err
      CALL ipslerr (3,TRIM(c_p), &
 &     'Problem in allocation of tmp_julian','', &
 &     '(you must increase memory)')
    ENDIF
    tszij = SIZE(t_julian,DIM=2)
    tmp_julian(:,1:tszij) = t_julian(:,1:tszij)
    DEALLOCATE(t_julian)
    ALLOCATE (t_julian(max_file,tax_size_in(nb_fi)),STAT=i_err)
    IF (i_err/=0) THEN
      WRITE(*,*) "ERROR IN ALLOCATION of t_julian : ",i_err
      CALL ipslerr (3,TRIM(c_p), &
 &      'Problem in reallocation of t_julian','', &
 &      '(you must increase memory)')
    ENDIF
    t_julian(:,1:tszij) = tmp_julian(:,1:tszij)
  ENDIF
!-----------------------
END SUBROUTINE rest_atim
!===
SUBROUTINE rest_alloc (i_buff,i_qsz,l_msg,c_p)
!---------------------------------------------------------------------
! Called by "c_p", allocate a temporary buffer
! (buff_tmp[1/2] depending on "i_buff" value) to the size "i_qsz".
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: i_buff,i_qsz
  LOGICAL,INTENT(IN) :: l_msg
  CHARACTER(LEN=*),INTENT(IN) :: c_p
!-
  INTEGER :: i_bsz,i_err
  LOGICAL :: l_alloc1,l_alloc2
  CHARACTER(LEN=9) :: cbn
  CHARACTER(LEN=5) :: c_err
!---------------------------------------------------------------------
  IF      (i_buff == 1) THEN
    IF (ALLOCATED(buff_tmp1)) THEN
      i_bsz = SIZE(buff_tmp1)
    ELSE
      i_bsz = 0
    ENDIF
    l_alloc1 =    (.NOT.ALLOCATED(buff_tmp1)) &
 &            .OR.((ALLOCATED(buff_tmp1)).AND.(i_qsz > i_bsz))
    l_alloc2 = .FALSE.
    cbn = 'buff_tmp1'
  ELSE IF (i_buff == 2) THEN
    IF (ALLOCATED(buff_tmp2)) THEN
      i_bsz = SIZE(buff_tmp2)
    ELSE
      i_bsz = 0
    ENDIF
    l_alloc1 = .FALSE.
    l_alloc2 =    (.NOT.ALLOCATED(buff_tmp2)) &
 &            .OR.((ALLOCATED(buff_tmp2)).AND.(i_qsz > i_bsz))
    cbn = 'buff_tmp2'
  ELSE
    CALL ipslerr (3,'rest_alloc', &
 &    'Called by '//TRIM(c_p),'with a wrong value of i_buff','')
  ENDIF
!-
!-
  IF (l_alloc1.OR.l_alloc2) THEN
    IF (l_msg) THEN
      IF (    (l_alloc1.AND.ALLOCATED(buff_tmp1)) &
 &        .OR.(l_alloc2.AND.ALLOCATED(buff_tmp2)) ) THEN
        WRITE(*,*) TRIM(c_p)//" : re_allocate "//TRIM(cbn)//"=",i_qsz
      ELSE
        WRITE(*,*) TRIM(c_p)//" : allocate "//TRIM(cbn)//"=",i_qsz
      ENDIF
    ENDIF
    IF (l_alloc1) THEN
      IF (ALLOCATED(buff_tmp1)) THEN
        DEALLOCATE(buff_tmp1)
      ENDIF
      ALLOCATE (buff_tmp1(i_qsz),STAT=i_err)
    ELSE
      IF (ALLOCATED(buff_tmp2)) THEN
        DEALLOCATE(buff_tmp2)
      ENDIF
      ALLOCATE (buff_tmp2(i_qsz),STAT=i_err)
    ENDIF
    IF (i_err /= 0) THEN
      WRITE (UNIT=c_err,FMT='(I5)') i_err
      CALL ipslerr (3,TRIM(c_p), &
 &      'Problem in allocation of',TRIM(cbn), &
 &      'Error : '//TRIM(c_err)//' (you must increase memory)')
    ENDIF
  ENDIF
!------------------------
END SUBROUTINE rest_alloc
!===
SUBROUTINE ioconf_setatt (attname,value)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: attname,value
!-
  CHARACTER(LEN=LEN_TRIM(attname)) :: tmp_str
!---------------------------------------------------------------------
  tmp_str = attname
  CALL strlowercase (tmp_str)
!-
  SELECT CASE(tmp_str)
    CASE('units')
      rest_units = value
    CASE('long_name')
      rest_lname = value
    CASE DEFAULT
      CALL ipslerr (2,'ioconf_restatt', &
        'The attribute name provided is unknown',attname,' ')
  END SELECT
!---------------------------
END SUBROUTINE ioconf_setatt
!===
SUBROUTINE ioget_vdim (fid,vname_q,varnbdim_max,varnbdim,vardims)
!---------------------------------------------------------------------
!- This routine allows the user to get the dimensions
!- of a field in the restart file.
!- This is the file which is read.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: fid
  CHARACTER(LEN=*) :: vname_q
  INTEGER,INTENT(IN) :: varnbdim_max
  INTEGER,INTENT(OUT) ::  varnbdim
  INTEGER,DIMENSION(varnbdim_max),INTENT(OUT) :: vardims
!-
  INTEGER :: vnb
!---------------------------------------------------------------------
! Find the index of the variable
  CALL find_str (varname_in(fid,1:nbvar_in(fid)),vname_q,vnb)
!-
  IF (vnb > 0) THEN
    varnbdim = varnbdim_in(fid,vnb)
    IF (varnbdim_max < varnbdim) THEN
      CALL ipslerr (3,'ioget_vdim', &
        'The provided array for the variable dimensions is too small', &
        '','')
    ELSE
      vardims(1:varnbdim) = vardims_in(fid,vnb,1:varnbdim)
    ENDIF
  ELSE
    varnbdim = 0
    CALL ipslerr (2,'ioget_vdim', &
      'Variable '//TRIM(vname_q)//' not found','','')
  ENDIF
!------------------------
END SUBROUTINE ioget_vdim
!===
SUBROUTINE ioget_vname (fid,nbvar,varnames)
!---------------------------------------------------------------------
!- This routine allows the user to extract the list
!- of variables in an opened restart file.
!- This is the file which is read
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: fid
  INTEGER,INTENT(OUT) ::  nbvar
  CHARACTER(LEN=*),INTENT(OUT) :: varnames(:)
!---------------------------------------------------------------------
  nbvar = nbvar_in(fid)
!-
  IF (SIZE(varnames) < nbvar) THEN
    CALL ipslerr (3,'ioget_vname', &
      'The provided array for the variable names is too small','','')
  ELSE
    varnames(1:nbvar) = varname_in(fid,1:nbvar)
  ENDIF
!-------------------------
END SUBROUTINE ioget_vname
!===
SUBROUTINE ioconf_expval (new_exp_val)
!---------------------------------------------------------------------
!- The default value written into the variables which are not
!- in the restart file can only be changed once.
!- This avoids further complications.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  REAL :: new_exp_val
!---------------------------------------------------------------------
  IF (.NOT.lock_valexp) THEN
    lock_valexp = .TRUE.
    val_exp = new_exp_val
  ELSE
    CALL ipslerr (2,'ioconf_expval', &
     'The default value for variable' &
   //'not available in the restart file ', &
     'has already been locked and can not be changed at this point', &
     ' ')
  ENDIF
!---------------------------
END SUBROUTINE ioconf_expval
!===
SUBROUTINE ioget_expval (get_exp_val)
!---------------------------------------------------------------------
!- Once the user has extracted the default value,
!- we lock it so that it can not be changed anymore.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  REAL :: get_exp_val
!---------------------------------------------------------------------
  get_exp_val = val_exp
  lock_valexp = .TRUE.
!--------------------------
END SUBROUTINE ioget_expval
!===
SUBROUTINE restclo (fid)
!---------------------------------------------------------------------
!- This subroutine closes one or any opened restart file.
!-
!- INPUT
!-
!- fid    : File ID in the restcom system (not the netCDF ID)(optional)
!-
!- OUTPUT
!-
!- NONE
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(in),OPTIONAL :: fid
!-
  INTEGER :: iret,ifnc
  CHARACTER(LEN=6) :: n_e
  CHARACTER(LEN=3) :: n_f
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (PRESENT(fid)) THEN
!---
    IF (l_dbg) THEN
      WRITE(*,*) &
        'restclo : Closing specified restart file number :', &
        fid,netcdf_id(fid,1:2)
    ENDIF
!---
    IF (netcdf_id(fid,1) > 0) THEN
      iret = NF90_CLOSE(netcdf_id(fid,1))
      IF (iret /= NF90_NOERR) THEN
        WRITE (n_e,'(I6)') iret
        WRITE (n_f,'(I3)') netcdf_id(fid,1)
        CALL ipslerr (2,'restclo', &
          "Error "//n_e//" in closing file : "//n_f,'',' ')
      ENDIF
      IF (netcdf_id(fid,1) == netcdf_id(fid,2)) THEN
        netcdf_id(fid,2) = -1
      ENDIF
      netcdf_id(fid,1) = -1
    ENDIF
!---
    IF (netcdf_id(fid,2) > 0)  THEN
      iret = NF90_CLOSE(netcdf_id(fid,2))
      IF (iret /= NF90_NOERR) THEN
        WRITE (n_e,'(I6)') iret
        WRITE (n_f,'(I3)') netcdf_id(fid,2)
        CALL ipslerr (2,'restclo', &
          "Error "//n_e//" in closing file : "//n_f,'',' ')
      ENDIF
      netcdf_id(fid,2) = -1
    ENDIF
!---
  ELSE
!---
    IF (l_dbg) WRITE(*,*) 'restclo : Closing all files'
!---
    DO ifnc=1,nb_fi
      IF (netcdf_id(ifnc,1) > 0) THEN
        iret = NF90_CLOSE(netcdf_id(ifnc,1))
        IF (iret /= NF90_NOERR) THEN
          WRITE (n_e,'(I6)') iret
          WRITE (n_f,'(I3)') netcdf_id(ifnc,1)
          CALL ipslerr (2,'restclo', &
            "Error "//n_e//" in closing file : "//n_f,'',' ')
        ENDIF
        IF (netcdf_id(ifnc,1) == netcdf_id(ifnc,2)) THEN
          netcdf_id(ifnc,2) = -1
        ENDIF
        netcdf_id(ifnc,1) = -1
      ENDIF
!-----
      IF (netcdf_id(ifnc,2) > 0) THEN
        iret = NF90_CLOSE(netcdf_id(ifnc,2))
        IF (iret /= NF90_NOERR) THEN
          WRITE (n_e,'(I6)') iret
          WRITE (n_f,'(I3)') netcdf_id(ifnc,2)
          CALL ipslerr (2,'restclo', &
            "Error "//n_e//" in closing file : "//n_f,'',' ')
        END IF
        netcdf_id(ifnc,2) = -1
      ENDIF
    ENDDO
  ENDIF
!---------------------
END SUBROUTINE restclo
!===
!-----------------
END MODULE restcom
