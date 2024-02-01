PROGRAM ncregular
!
!$Id: ncregular.f90 2281 2010-10-15 14:21:13Z smasson $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!---------------------------------------------------------------------
!- This code replaces a 2D surface grid by vectors.
!- Obviously it only works if you have a regular grid.
!-
!- Jan Polcher (polcher@lmd.jussieu.fr)
!- Jacques Bellier (jacques.bellier@cea.fr)
!---------------------------------------------------------------------
  USE netcdf
!-
  IMPLICIT NONE
!-
  INTEGER :: iread, if, in, iv, sz
  INTEGER :: ier, nb_files, iret, ndims, nvars, nb_glat
  INTEGER :: lon_dim_id, lat_dim_id
  INTEGER :: lon_len, lat_len, lon_id, lat_id
  INTEGER :: nav_lon_id, nav_lat_id
  INTEGER :: alloc_stat_lon, alloc_stat_lat
!-
  INTEGER,ALLOCATABLE ::  file_id(:), tax_id(:)
  CHARACTER(LEN=80),ALLOCATABLE :: names(:)
  CHARACTER(LEN=80) :: dim_name
  CHARACTER(LEN=80) :: varname
  CHARACTER(LEN=20) :: xname, yname, lonname, latname
  LOGICAL :: check, regular
!-
  REAL,ALLOCATABLE :: lon(:), lat(:), lon2(:), lat2(:)
  REAL,ALLOCATABLE :: del_lon(:), del_lat(:)
!-
  INTEGER iargc, getarg
  EXTERNAL iargc, getarg
!---------------------------------------------------------------------
  alloc_stat_lon = 0
  alloc_stat_lat = 0
!-
  iread = iargc()
!-
  ALLOCATE (names(iread),stat=ier)
  IF (ier /= 0) THEN
    WRITE (*,*) ' Could not allocate names of size ', iread
    STOP 'nctax'
  ENDIF
!-
  CALL nct_getarg (iread, nb_files, names, check, &
 &                 xname, yname, lonname, latname)
!-
! Allocate space
!-
  ALLOCATE (file_id(nb_files),stat=ier)
  IF (ier /= 0) THEN
    WRITE (*,*) ' Could not allocate file_id of size ', nb_files
    STOP 'nctax'
  ENDIF
!-
  ALLOCATE (tax_id(nb_files),stat=ier)
  IF (ier /= 0) THEN
    WRITE (*,*) ' Could not allocate tax_id of size ', nb_files
    STOP 'nctax'
  ENDIF
!-
  DO if=1,nb_files
!---
    IF (check) THEN
      WRITE(*,*) 'ncregular : ', if,  names(if)
    ENDIF
!---
    iret = NF90_OPEN (names(if),NF90_WRITE,file_id(if))
    iret = NF90_INQUIRE (file_id(if),ndims,nvars,nb_glat,tax_id(if))
!---
!-- Get the IDs of the variables
!---
    lon_len = -9999
    lat_len = -9999
    DO in=1,ndims
!-----
      iret = NF90_INQUIRE_DIMENSION (file_id(if), in, dim_name, sz)
!-----
      IF (     (LEN_TRIM(dim_name) == 1) &
 &        .AND.(INDEX(dim_name,TRIM(xname)) == 1) ) THEN
        lon_dim_id = in
        lon_len = sz
      ENDIF
!-----
      IF (     (LEN_TRIM(dim_name) == 1) &
 &        .AND.(INDEX(dim_name,TRIM(yname)) == 1) ) THEN
        lat_dim_id = in
        lat_len = sz
      ENDIF
!-----
    ENDDO
!---
    IF ( (lon_len == -9999).OR.(lat_len == -9999) ) THEN
      WRITE(*,*) 'ncregular : The specified dimensions were not'
      WRITE(*,*) 'found in file : ',names(if)
      iret = NF90_CLOSE (file_id(if))
      STOP
    ENDIF
!---
    IF (check) THEN
      WRITE(*,*) 'ncregular : lon_dim_id, lon_len',lon_dim_id,lon_len
      WRITE(*,*) 'ncregular : lat_dim_id, lat_len',lat_dim_id,lat_len
    ENDIF
!---
!-- Look for the right variables
!---
    nav_lon_id = -9999
    nav_lat_id = -9999
    DO iv=1,nvars
      iret = NF90_INQUIRE_VARIABLE (file_id(if),iv,name=varname)
      IF (INDEX(varname,TRIM(lonname)) > 0) THEN
        nav_lon_id = iv
      ENDIF
      IF (INDEX(varname,TRIM(latname)) > 0) THEN
        nav_lat_id = iv
      ENDIF
    ENDDO
!---
    IF ( (nav_lon_id == -9999).OR.(nav_lat_id == -9999) ) THEN
      WRITE(*,*) 'ncregular : The specified coordinate fields'
      WRITE(*,*) 'were not found in file : ',names(if)
      iret = NF90_CLOSE (file_id(if))
      STOP
    ENDIF
!---
    IF (check) THEN
      WRITE(*,*) 'ncregular : nav_lon_id :', nav_lon_id
      WRITE(*,*) 'ncregular : nav_lat_id :', nav_lat_id
    ENDIF
!---
!-- Read variables from file and check if regular
!---
!-- Do we have the variable to read the
!---
    IF ( alloc_stat_lon < lon_len) THEN
      IF ( alloc_stat_lon > 0) THEN
        deallocate(lon)
        deallocate(lon2)
        deallocate(del_lon)
      ENDIF
      allocate(lon(lon_len))
      allocate(lon2(lon_len))
      allocate(del_lon(lon_len))
      alloc_stat_lon = lon_len
    ENDIF
!---
    IF ( alloc_stat_lat < lat_len) THEN
      IF ( alloc_stat_lat > 0) THEN
        deallocate(lat)
        deallocate(lat2)
        deallocate(del_lat)
      ENDIF
      allocate(lat(lat_len))
      allocate(lat2(lat_len))
      allocate(del_lat(lat_len))
      alloc_stat_lat = lat_len
    ENDIF
!---
!-- Read data
!---
    iret = NF90_GET_VAR (file_id(if),nav_lon_id,lon, &
 & start=(/1,1/),count=(/lon_len,1/),stride=(/1,1/))
    iret = NF90_GET_VAR (file_id(if),nav_lon_id,lon2, &
 & start=(/1,int(lat_len/2)/),count=(/lon_len,1/),stride=(/1,1/))
    del_lon = lon-lon2
!-
    iret = NF90_GET_VAR (file_id(if),nav_lat_id,lat, &
 & start=(/1,1/),count=(/1,lat_len/),stride=(/lon_len,1/))
    iret = NF90_GET_VAR (file_id(if),nav_lat_id,lat2, &
 & start=(/int(lon_len/2),1/),count=(/1,lat_len/),stride=(/lon_len,1/))
    del_lat = lat-lat2
!-
    regular = (    (MAXVAL(del_lon) < 0.001) &
 &             .OR.(MAXVAL(del_lat) < 0.001) )
!---
!-- Create the new variables
!---
    IF (regular) THEN
      IF (check) THEN
        WRITE(*,*) 'Regular case'
      ENDIF
      iret = NF90_REDEF (file_id(if))
      iret = NF90_RENAME_DIM (file_id(if), lon_dim_id, 'lon')
      iret = NF90_RENAME_DIM (file_id(if), lat_dim_id, 'lat')
      IF (check) THEN
        WRITE(*,*) 'Dimensions renamed'
      ENDIF
      iret = NF90_DEF_VAR (file_id(if), 'lon', NF90_FLOAT, &
 &                         lon_dim_id, lon_id)
      iret = NF90_DEF_VAR (file_id(if), 'lat', NF90_FLOAT, &
 &                         lat_dim_id, lat_id)
      IF (check) THEN
        WRITE(*,*) 'New variables defined'
      ENDIF
!-----
!---- Copy attributes
!-----
      iret = NF90_COPY_ATT (file_id(if),nav_lon_id,'units', &
 &                          file_id(if),lon_id)
      iret = NF90_COPY_ATT (file_id(if),nav_lon_id,'title', &
 &                          file_id(if),lon_id)
      iret = NF90_COPY_ATT (file_id(if),nav_lon_id,'valid_max', &
 &                          file_id(if),lon_id)
      iret = NF90_COPY_ATT (file_id(if),nav_lon_id,'valid_min', &
 &                          file_id(if),lon_id)
!-----
      iret = NF90_COPY_ATT (file_id(if),nav_lat_id,'units', &
 &                          file_id(if),lat_id)
      iret = NF90_COPY_ATT (file_id(if),nav_lat_id,'title', &
 &                          file_id(if),lat_id)
      iret = NF90_COPY_ATT (file_id(if),nav_lat_id,'valid_max', &
 &                          file_id(if),lat_id)
      iret = NF90_COPY_ATT (file_id(if),nav_lat_id,'valid_min', &
 &                          file_id(if),lat_id)
!-----
!---- Go into write mode
!-----
      iret = NF90_ENDDEF (file_id(if))
!-----
!---- Write data
!-----
      iret = NF90_PUT_VAR (file_id(if),lon_id,lon(1:lon_len))
      iret = NF90_PUT_VAR (file_id(if),lat_id,lat(1:lat_len))
!-
      iret = NF90_CLOSE (file_id(if))
    ELSE
      WRITE(*,*) 'ncregular : Your grid is not regular'
      WRITE(*,*) names(if), 'remains unchanged'
      iret = NF90_CLOSE (file_id(if))
    ENDIF
!-
  ENDDO
!--------------------
END PROGRAM ncregular
!-
!===
!-
SUBROUTINE nct_getarg (argx, nb_files, names, check, &
 &                     xname, yname, lonname, latname)
!---------------------------------------------------------------------
!- Read the arguments of nctax.
!---------------------------------------------------------------------
  INTEGER,INTENT(in) :: argx
  INTEGER, INTENT(out) :: nb_files
  CHARACTER(LEN=80),INTENT(out) :: names(argx)
  CHARACTER(LEN=20) :: xname, yname, lonname, latname
!-
  CHARACTER(LEN=80) :: tmp, tmp_arg
  LOGICAL :: check
!---------------------------------------------------------------------
  check = .FALSE.
!-
! Get the number of arguments
!-
  nb_files = 0
!-
  xname = 'x'
  yname = 'y'
  lonname = 'nav_lon'
  latname = 'nav_lat'
!-
! Go through the arguments and analyse them one by one
!-
  IF (check) WRITE(*,*) 'Start going through the arguments'
!-
  IF (argx == 0) THEN
    WRITE(*,*) 'To get usage : nctax -h '
    STOP
  ENDIF
!-
  iread = 1
  DO WHILE (iread <= argx)
    iret = getarg(iread,tmp)
    IF (check) WRITE(*,*) ' iread, tmp :', iread, tmp
    SELECTCASE(tmp)
      CASE('-d')
        WRITE(*,*) 'DEBUG MODE SELECTED'
        check = .TRUE.
        iread = iread+1
      CASE('-h')
        WRITE(*,*) 'Usage : nregular [options] file1 [file2 ...]'
        WRITE(*,*) '    -d : Verbose mode'
        WRITE(*,*) '    -h : This output'
        STOP
      CASE('-dim_lon')
        iread = iread+1
        iret  = getarg(iread,tmp_arg)
        xname = TRIM(tmp_arg)
        iread = iread+1
      CASE('-dim_lat')
        iread = iread+1
        iret  = getarg(iread,tmp_arg)
        yname = TRIM(tmp_arg)
        iread = iread+1
      CASE('-coo_lon')
        iread = iread+1
        iret  = getarg(iread,tmp_arg)
        lonname = TRIM(tmp_arg)
        iread = iread+1
      CASE('-coo_lat')
        iread = iread+1
        iret  = getarg(iread,tmp_arg)
        latname = TRIM(tmp_arg)
        iread = iread+1
      CASE DEFAULT
        IF (check) WRITE(*,*) 'nct_getarg : CASE default'
        IF (INDEX(tmp,'-') /= 1) THEN
          nb_files = nb_files+1
          names(nb_files) = tmp
          iread = iread+1
        ELSE
          WRITE(*,*) "WARNING Unknown option ",tmp
          WRITE(*,*) "For ore information : nctax -h"
        ENDIF
    END SELECT
  ENDDO
!-
  IF (check) THEN
    WRITE(*,*) ' nct_getarg : output >> '
    WRITE(*,*) '>> nb_files : ', nb_files
    WRITE(*,*) '>> names :', (names(ii), ii=1,nb_files)
  ENDIF
!------------------------
END SUBROUTINE nct_getarg
