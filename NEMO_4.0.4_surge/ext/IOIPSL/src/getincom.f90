MODULE getincom
!-
!$Id: getincom.f90 2281 2010-10-15 14:21:13Z smasson $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!---------------------------------------------------------------------
USE errioipsl, ONLY : ipslerr
USE stringop, &
 &   ONLY : nocomma,cmpblank,strlowercase
!-
IMPLICIT NONE
!-
PRIVATE
PUBLIC :: getin_name, getin, getin_dump
!-
!!--------------------------------------------------------------------
!! The "getin_name" routine allows the user to change the name
!! of the definition file in which the data will be read.
!! ("run.def" by default)
!!
!!  SUBROUTINE getin_name (file_name)
!!
!! OPTIONAL INPUT argument
!!
!! (C) file_name :  the name of the file
!!                  in which the data will be read
!!--------------------------------------------------------------------
!-
!-
INTERFACE getin
!!--------------------------------------------------------------------
!! The "getin" routines get a variable.
!! We first check if we find it in the database
!! and if not we get it from the definition file.
!!
!! SUBROUTINE getin (target,ret_val)
!!
!! INPUT
!!
!! (C) target : Name of the variable
!!
!! OUTPUT
!!
!! (I/R/C/L) ret_val : scalar, vector or matrix that will contain
!!                     that will contain the (standard)
!!                     integer/real/character/logical values
!!--------------------------------------------------------------------
  MODULE PROCEDURE getinrs, getinr1d, getinr2d, &
 &                 getinis, getini1d, getini2d, &
 &                 getincs, getinc1d, getinc2d, &
 &                 getinls, getinl1d, getinl2d
END INTERFACE
!-
!!--------------------------------------------------------------------
!! The "getin_dump" routine will dump the content of the database
!! into a file which has the same format as the definition file.
!! The idea is that the user can see which parameters were used
!! and re-use the file for another run.
!!
!!  SUBROUTINE getin_dump (fileprefix)
!!
!! OPTIONAL INPUT argument
!!
!! (C) fileprefix : allows the user to change the name of the file
!!                  in which the data will be archived
!!--------------------------------------------------------------------
!-
  INTEGER,PARAMETER :: max_files=100
  CHARACTER(LEN=100),DIMENSION(max_files),SAVE :: filelist
  INTEGER,SAVE      :: nbfiles
!-
  INTEGER,SAVE :: allread=0
  CHARACTER(LEN=100),SAVE :: def_file = 'run.def'
!-
  INTEGER,PARAMETER :: i_txtslab=1000,l_n=30
  INTEGER,SAVE :: nb_lines,i_txtsize=0
  CHARACTER(LEN=100),SAVE,ALLOCATABLE,DIMENSION(:) :: fichier
  CHARACTER(LEN=l_n),SAVE,ALLOCATABLE,DIMENSION(:) :: targetlist
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: fromfile,compline
!-
  INTEGER,PARAMETER :: n_d_fmt=5,max_msgs=15
  CHARACTER(LEN=6),SAVE :: c_i_fmt = '(I5.5)'
!-
! The data base of parameters
!-
  INTEGER,PARAMETER :: memslabs=200
  INTEGER,PARAMETER :: compress_lim=20
!-
  INTEGER,SAVE :: nb_keys=0
  INTEGER,SAVE :: keymemsize=0
!-
! keystr definition
! name of a key
!-
! keystatus definition
! keystatus = 1 : Value comes from the file defined by 'def_file'
! keystatus = 2 : Default value is used
! keystatus = 3 : Some vector elements were taken from default
!-
! keytype definition
! keytype = 1 : Integer
! keytype = 2 : Real
! keytype = 3 : Character
! keytype = 4 : Logical
!-
  INTEGER,PARAMETER :: k_i=1, k_r=2, k_c=3, k_l=4
!-
! Allow compression for keys (only for integer and real)
! keycompress < 0 : not compressed
! keycompress > 0 : number of repeat of the value
!-
TYPE :: t_key
  CHARACTER(LEN=l_n) :: keystr
  INTEGER :: keystatus, keytype, keycompress, &
 &           keyfromfile, keymemstart, keymemlen
END TYPE t_key
!-
  TYPE(t_key),SAVE,ALLOCATABLE,DIMENSION(:) :: key_tab
!-
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: i_mem
  INTEGER,SAVE :: i_memsize=0, i_mempos=0
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: r_mem
  INTEGER,SAVE :: r_memsize=0, r_mempos=0
  CHARACTER(LEN=100),SAVE,ALLOCATABLE,DIMENSION(:) :: c_mem
  INTEGER,SAVE :: c_memsize=0, c_mempos=0
  LOGICAL,SAVE,ALLOCATABLE,DIMENSION(:) :: l_mem
  INTEGER,SAVE :: l_memsize=0, l_mempos=0
!-
CONTAINS
!-
!=== DEFINITION FILE NAME INTERFACE
!-
SUBROUTINE getin_name (cname)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: cname
!---------------------------------------------------------------------
  IF (allread == 0) THEN
    def_file = ADJUSTL(cname)
  ELSE
    CALL ipslerr (3,'getin_name', &
 &   'The name of the database file (any_name.def)', &
 &   'must be changed *before* any attempt','to read the database.')
  ENDIF
!------------------------
END SUBROUTINE getin_name
!-
!=== INTEGER INTERFACE
!-
SUBROUTINE getinis (target,ret_val)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  INTEGER :: ret_val
!-
  INTEGER,DIMENSION(1) :: tmp_ret_val
  INTEGER :: pos,status=0,fileorig
!---------------------------------------------------------------------
!-
! Do we have this target in our database ?
!-
  CALL get_findkey (1,target,pos)
!-
  tmp_ret_val(1) = ret_val
!-
  IF (pos < 0) THEN
!-- Get the information out of the file
    CALL get_fil (target,status,fileorig,i_val=tmp_ret_val)
!-- Put the data into the database
    CALL get_wdb &
 &   (target,status,fileorig,1,i_val=tmp_ret_val)
  ELSE
!-- Get the value out of the database
    CALL get_rdb (pos,1,target,i_val=tmp_ret_val)
  ENDIF
  ret_val = tmp_ret_val(1)
!---------------------
END SUBROUTINE getinis
!===
SUBROUTINE getini1d (target,ret_val)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  INTEGER,DIMENSION(:) :: ret_val
!-
  INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
  INTEGER,SAVE :: tmp_ret_size = 0
  INTEGER :: pos,size_of_in,status=0,fileorig
!---------------------------------------------------------------------
!-
! Do we have this target in our database ?
!-
  CALL get_findkey (1,target,pos)
!-
  size_of_in = SIZE(ret_val)
  IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
    ALLOCATE (tmp_ret_val(size_of_in))
  ELSE IF (size_of_in > tmp_ret_size) THEN
    DEALLOCATE (tmp_ret_val)
    ALLOCATE (tmp_ret_val(size_of_in))
    tmp_ret_size = size_of_in
  ENDIF
  tmp_ret_val(1:size_of_in) = ret_val(1:size_of_in)
!-
  IF (pos < 0) THEN
!-- Get the information out of the file
    CALL get_fil (target,status,fileorig,i_val=tmp_ret_val)
!-- Put the data into the database
    CALL get_wdb &
 &   (target,status,fileorig,size_of_in,i_val=tmp_ret_val)
  ELSE
!-- Get the value out of the database
    CALL get_rdb (pos,size_of_in,target,i_val=tmp_ret_val)
  ENDIF
  ret_val(1:size_of_in) = tmp_ret_val(1:size_of_in)
!----------------------
END SUBROUTINE getini1d
!===
SUBROUTINE getini2d (target,ret_val)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  INTEGER,DIMENSION(:,:) :: ret_val
!-
  INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
  INTEGER,SAVE :: tmp_ret_size = 0
  INTEGER :: pos,size_of_in,size_1,size_2,status=0,fileorig
  INTEGER :: jl,jj,ji
!---------------------------------------------------------------------
!-
! Do we have this target in our database ?
!-
  CALL get_findkey (1,target,pos)
!-
  size_of_in = SIZE(ret_val)
  size_1 = SIZE(ret_val,1)
  size_2 = SIZE(ret_val,2)
  IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
    ALLOCATE (tmp_ret_val(size_of_in))
  ELSE IF (size_of_in > tmp_ret_size) THEN
    DEALLOCATE (tmp_ret_val)
    ALLOCATE (tmp_ret_val(size_of_in))
    tmp_ret_size = size_of_in
  ENDIF
!-
  jl=0
  DO jj=1,size_2
    DO ji=1,size_1
      jl=jl+1
      tmp_ret_val(jl) = ret_val(ji,jj)
    ENDDO
  ENDDO
!-
  IF (pos < 0) THEN
!-- Get the information out of the file
    CALL get_fil (target,status,fileorig,i_val=tmp_ret_val)
!-- Put the data into the database
    CALL get_wdb &
 &   (target,status,fileorig,size_of_in,i_val=tmp_ret_val)
  ELSE
!-- Get the value out of the database
    CALL get_rdb (pos,size_of_in,target,i_val=tmp_ret_val)
  ENDIF
!-
  jl=0
  DO jj=1,size_2
    DO ji=1,size_1
      jl=jl+1
      ret_val(ji,jj) = tmp_ret_val(jl)
    ENDDO
  ENDDO
!----------------------
END SUBROUTINE getini2d
!-
!=== REAL INTERFACE
!-
SUBROUTINE getinrs (target,ret_val)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  REAL :: ret_val
!-
  REAL,DIMENSION(1) :: tmp_ret_val
  INTEGER :: pos,status=0,fileorig
!---------------------------------------------------------------------
!-
! Do we have this target in our database ?
!-
  CALL get_findkey (1,target,pos)
!-
  tmp_ret_val(1) = ret_val
!-
  IF (pos < 0) THEN
!-- Get the information out of the file
    CALL get_fil (target,status,fileorig,r_val=tmp_ret_val)
!-- Put the data into the database
    CALL get_wdb &
 &   (target,status,fileorig,1,r_val=tmp_ret_val)
  ELSE
!-- Get the value out of the database
    CALL get_rdb (pos,1,target,r_val=tmp_ret_val)
  ENDIF
  ret_val = tmp_ret_val(1)
!---------------------
END SUBROUTINE getinrs
!===
SUBROUTINE getinr1d (target,ret_val)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  REAL,DIMENSION(:) :: ret_val
!-
  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
  INTEGER,SAVE :: tmp_ret_size = 0
  INTEGER :: pos,size_of_in,status=0,fileorig
!---------------------------------------------------------------------
!-
! Do we have this target in our database ?
!-
  CALL get_findkey (1,target,pos)
!-
  size_of_in = SIZE(ret_val)
  IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
    ALLOCATE (tmp_ret_val(size_of_in))
  ELSE IF (size_of_in > tmp_ret_size) THEN
    DEALLOCATE (tmp_ret_val)
    ALLOCATE (tmp_ret_val(size_of_in))
    tmp_ret_size = size_of_in
  ENDIF
  tmp_ret_val(1:size_of_in) = ret_val(1:size_of_in)
!-
  IF (pos < 0) THEN
!-- Get the information out of the file
    CALL get_fil (target,status,fileorig,r_val=tmp_ret_val)
!-- Put the data into the database
    CALL get_wdb &
 &   (target,status,fileorig,size_of_in,r_val=tmp_ret_val)
  ELSE
!-- Get the value out of the database
    CALL get_rdb (pos,size_of_in,target,r_val=tmp_ret_val)
  ENDIF
  ret_val(1:size_of_in) = tmp_ret_val(1:size_of_in)
!----------------------
END SUBROUTINE getinr1d
!===
SUBROUTINE getinr2d (target,ret_val)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  REAL,DIMENSION(:,:) :: ret_val
!-
  REAL,DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
  INTEGER,SAVE :: tmp_ret_size = 0
  INTEGER :: pos,size_of_in,size_1,size_2,status=0,fileorig
  INTEGER :: jl,jj,ji
!---------------------------------------------------------------------
!-
! Do we have this target in our database ?
!-
  CALL get_findkey (1,target,pos)
!-
  size_of_in = SIZE(ret_val)
  size_1 = SIZE(ret_val,1)
  size_2 = SIZE(ret_val,2)
  IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
    ALLOCATE (tmp_ret_val(size_of_in))
  ELSE IF (size_of_in > tmp_ret_size) THEN
    DEALLOCATE (tmp_ret_val)
    ALLOCATE (tmp_ret_val(size_of_in))
    tmp_ret_size = size_of_in
  ENDIF
!-
  jl=0
  DO jj=1,size_2
    DO ji=1,size_1
      jl=jl+1
      tmp_ret_val(jl) = ret_val(ji,jj)
    ENDDO
  ENDDO
!-
  IF (pos < 0) THEN
!-- Get the information out of the file
    CALL get_fil (target,status,fileorig,r_val=tmp_ret_val)
!-- Put the data into the database
    CALL get_wdb &
 &   (target,status,fileorig,size_of_in,r_val=tmp_ret_val)
  ELSE
!-- Get the value out of the database
    CALL get_rdb (pos,size_of_in,target,r_val=tmp_ret_val)
  ENDIF
!-
  jl=0
  DO jj=1,size_2
    DO ji=1,size_1
      jl=jl+1
      ret_val(ji,jj) = tmp_ret_val(jl)
    ENDDO
  ENDDO
!----------------------
END SUBROUTINE getinr2d
!-
!=== CHARACTER INTERFACE
!-
SUBROUTINE getincs (target,ret_val)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  CHARACTER(LEN=*) :: ret_val
!-
  CHARACTER(LEN=100),DIMENSION(1) :: tmp_ret_val
  INTEGER :: pos,status=0,fileorig
!---------------------------------------------------------------------
!-
! Do we have this target in our database ?
!-
  CALL get_findkey (1,target,pos)
!-
  tmp_ret_val(1) = ret_val
!-
  IF (pos < 0) THEN
!-- Get the information out of the file
    CALL get_fil (target,status,fileorig,c_val=tmp_ret_val)
!-- Put the data into the database
    CALL get_wdb &
 &   (target,status,fileorig,1,c_val=tmp_ret_val)
  ELSE
!-- Get the value out of the database
    CALL get_rdb (pos,1,target,c_val=tmp_ret_val)
  ENDIF
  ret_val = tmp_ret_val(1)
!---------------------
END SUBROUTINE getincs
!===
SUBROUTINE getinc1d (target,ret_val)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  CHARACTER(LEN=*),DIMENSION(:) :: ret_val
!-
  CHARACTER(LEN=100),DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
  INTEGER,SAVE :: tmp_ret_size = 0
  INTEGER :: pos,size_of_in,status=0,fileorig
!---------------------------------------------------------------------
!-
! Do we have this target in our database ?
!-
  CALL get_findkey (1,target,pos)
!-
  size_of_in = SIZE(ret_val)
  IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
    ALLOCATE (tmp_ret_val(size_of_in))
  ELSE IF (size_of_in > tmp_ret_size) THEN
    DEALLOCATE (tmp_ret_val)
    ALLOCATE (tmp_ret_val(size_of_in))
    tmp_ret_size = size_of_in
  ENDIF
  tmp_ret_val(1:size_of_in) = ret_val(1:size_of_in)
!-
  IF (pos < 0) THEN
!-- Get the information out of the file
    CALL get_fil (target,status,fileorig,c_val=tmp_ret_val)
!-- Put the data into the database
    CALL get_wdb &
 &   (target,status,fileorig,size_of_in,c_val=tmp_ret_val)
  ELSE
!-- Get the value out of the database
    CALL get_rdb (pos,size_of_in,target,c_val=tmp_ret_val)
  ENDIF
  ret_val(1:size_of_in) = tmp_ret_val(1:size_of_in)
!----------------------
END SUBROUTINE getinc1d
!===
SUBROUTINE getinc2d (target,ret_val)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  CHARACTER(LEN=*),DIMENSION(:,:) :: ret_val
!-
  CHARACTER(LEN=100),DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
  INTEGER,SAVE :: tmp_ret_size = 0
  INTEGER :: pos,size_of_in,size_1,size_2,status=0,fileorig
  INTEGER :: jl,jj,ji
!---------------------------------------------------------------------
!-
! Do we have this target in our database ?
!-
  CALL get_findkey (1,target,pos)
!-
  size_of_in = SIZE(ret_val)
  size_1 = SIZE(ret_val,1)
  size_2 = SIZE(ret_val,2)
  IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
    ALLOCATE (tmp_ret_val(size_of_in))
  ELSE IF (size_of_in > tmp_ret_size) THEN
    DEALLOCATE (tmp_ret_val)
    ALLOCATE (tmp_ret_val(size_of_in))
    tmp_ret_size = size_of_in
  ENDIF
!-
  jl=0
  DO jj=1,size_2
    DO ji=1,size_1
      jl=jl+1
      tmp_ret_val(jl) = ret_val(ji,jj)
    ENDDO
  ENDDO
!-
  IF (pos < 0) THEN
!-- Get the information out of the file
    CALL get_fil (target,status,fileorig,c_val=tmp_ret_val)
!-- Put the data into the database
    CALL get_wdb &
 &   (target,status,fileorig,size_of_in,c_val=tmp_ret_val)
  ELSE
!-- Get the value out of the database
    CALL get_rdb (pos,size_of_in,target,c_val=tmp_ret_val)
  ENDIF
!-
  jl=0
  DO jj=1,size_2
    DO ji=1,size_1
      jl=jl+1
      ret_val(ji,jj) = tmp_ret_val(jl)
    ENDDO
  ENDDO
!----------------------
END SUBROUTINE getinc2d
!-
!=== LOGICAL INTERFACE
!-
SUBROUTINE getinls (target,ret_val)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  LOGICAL :: ret_val
!-
  LOGICAL,DIMENSION(1) :: tmp_ret_val
  INTEGER :: pos,status=0,fileorig
!---------------------------------------------------------------------
!-
! Do we have this target in our database ?
!-
  CALL get_findkey (1,target,pos)
!-
  tmp_ret_val(1) = ret_val
!-
  IF (pos < 0) THEN
!-- Get the information out of the file
    CALL get_fil (target,status,fileorig,l_val=tmp_ret_val)
!-- Put the data into the database
    CALL get_wdb &
 &   (target,status,fileorig,1,l_val=tmp_ret_val)
  ELSE
!-- Get the value out of the database
    CALL get_rdb (pos,1,target,l_val=tmp_ret_val)
  ENDIF
  ret_val = tmp_ret_val(1)
!---------------------
END SUBROUTINE getinls
!===
SUBROUTINE getinl1d (target,ret_val)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  LOGICAL,DIMENSION(:) :: ret_val
!-
  LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
  INTEGER,SAVE :: tmp_ret_size = 0
  INTEGER :: pos,size_of_in,status=0,fileorig
!---------------------------------------------------------------------
!-
! Do we have this target in our database ?
!-
  CALL get_findkey (1,target,pos)
!-
  size_of_in = SIZE(ret_val)
  IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
    ALLOCATE (tmp_ret_val(size_of_in))
  ELSE IF (size_of_in > tmp_ret_size) THEN
    DEALLOCATE (tmp_ret_val)
    ALLOCATE (tmp_ret_val(size_of_in))
    tmp_ret_size = size_of_in
  ENDIF
  tmp_ret_val(1:size_of_in) = ret_val(1:size_of_in)
!-
  IF (pos < 0) THEN
!-- Get the information out of the file
    CALL get_fil (target,status,fileorig,l_val=tmp_ret_val)
!-- Put the data into the database
    CALL get_wdb &
 &   (target,status,fileorig,size_of_in,l_val=tmp_ret_val)
  ELSE
!-- Get the value out of the database
    CALL get_rdb (pos,size_of_in,target,l_val=tmp_ret_val)
  ENDIF
  ret_val(1:size_of_in) = tmp_ret_val(1:size_of_in)
!----------------------
END SUBROUTINE getinl1d
!===
SUBROUTINE getinl2d (target,ret_val)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  LOGICAL,DIMENSION(:,:) :: ret_val
!-
  LOGICAL,DIMENSION(:),ALLOCATABLE,SAVE :: tmp_ret_val
  INTEGER,SAVE :: tmp_ret_size = 0
  INTEGER :: pos,size_of_in,size_1,size_2,status=0,fileorig
  INTEGER :: jl,jj,ji
!---------------------------------------------------------------------
!-
! Do we have this target in our database ?
!-
  CALL get_findkey (1,target,pos)
!-
  size_of_in = SIZE(ret_val)
  size_1 = SIZE(ret_val,1)
  size_2 = SIZE(ret_val,2)
  IF (.NOT.ALLOCATED(tmp_ret_val)) THEN
    ALLOCATE (tmp_ret_val(size_of_in))
  ELSE IF (size_of_in > tmp_ret_size) THEN
    DEALLOCATE (tmp_ret_val)
    ALLOCATE (tmp_ret_val(size_of_in))
    tmp_ret_size = size_of_in
  ENDIF
!-
  jl=0
  DO jj=1,size_2
    DO ji=1,size_1
      jl=jl+1
      tmp_ret_val(jl) = ret_val(ji,jj)
    ENDDO
  ENDDO
!-
  IF (pos < 0) THEN
!-- Get the information out of the file
    CALL get_fil (target,status,fileorig,l_val=tmp_ret_val)
!-- Put the data into the database
    CALL get_wdb &
 &   (target,status,fileorig,size_of_in,l_val=tmp_ret_val)
  ELSE
!-- Get the value out of the database
    CALL get_rdb (pos,size_of_in,target,l_val=tmp_ret_val)
  ENDIF
!-
  jl=0
  DO jj=1,size_2
    DO ji=1,size_1
      jl=jl+1
      ret_val(ji,jj) = tmp_ret_val(jl)
    ENDDO
  ENDDO
!----------------------
END SUBROUTINE getinl2d
!-
!=== Generic file/database INTERFACE
!-
SUBROUTINE get_fil (target,status,fileorig,i_val,r_val,c_val,l_val)
!---------------------------------------------------------------------
!- Subroutine that will extract from the file the values
!- attributed to the keyword target
!-
!- (C) target    : target for which we will look in the file
!- (I) status    : tells us from where we obtained the data
!- (I) fileorig  : index of the file from which the key comes
!- (I) i_val(:)  : INTEGER(nb_to_ret)   values
!- (R) r_val(:)  : REAL(nb_to_ret)      values
!- (L) l_val(:)  : LOGICAL(nb_to_ret)   values
!- (C) c_val(:)  : CHARACTER(nb_to_ret) values
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  INTEGER,INTENT(OUT) :: status,fileorig
  INTEGER,DIMENSION(:),OPTIONAL          :: i_val
  REAL,DIMENSION(:),OPTIONAL             :: r_val
  LOGICAL,DIMENSION(:),OPTIONAL          :: l_val
  CHARACTER(LEN=*),DIMENSION(:),OPTIONAL :: c_val
!-
  INTEGER :: k_typ,nb_to_ret,it,pos,len_str,status_cnt,io_err
  CHARACTER(LEN=n_d_fmt)  :: cnt
  CHARACTER(LEN=80) :: str_READ,str_READ_lower
  CHARACTER(LEN=9)  :: c_vtyp
  LOGICAL,DIMENSION(:),ALLOCATABLE :: found
  LOGICAL :: def_beha,compressed
  CHARACTER(LEN=10) :: c_fmt
  INTEGER :: i_cmpval
  REAL    :: r_cmpval
  INTEGER :: ipos_tr,ipos_fl
!---------------------------------------------------------------------
!-
! Get the type of the argument
  CALL get_qtyp (k_typ,c_vtyp,i_val,r_val,c_val,l_val)
  SELECT CASE (k_typ)
  CASE(k_i)
    nb_to_ret = SIZE(i_val)
  CASE(k_r)
    nb_to_ret = SIZE(r_val)
  CASE(k_c)
    nb_to_ret = SIZE(c_val)
  CASE(k_l)
    nb_to_ret = SIZE(l_val)
  CASE DEFAULT
    CALL ipslerr (3,'get_fil', &
 &   'Internal error','Unknown type of data',' ')
  END SELECT
!-
! Read the file(s)
  CALL getin_read
!-
! Allocate and initialize the memory we need
  ALLOCATE(found(nb_to_ret))
  found(:) = .FALSE.
!-
! See what we find in the files read
  DO it=1,nb_to_ret
!---
!-- First try the target as it is
    CALL get_findkey (2,target,pos)
!---
!-- Another try
!---
    IF (pos < 0) THEN
      WRITE(UNIT=cnt,FMT=c_i_fmt) it
      CALL get_findkey (2,TRIM(target)//'__'//cnt,pos)
    ENDIF
!---
!-- We dont know from which file the target could come.
!-- Thus by default we attribute it to the first file :
    fileorig = 1
!---
    IF (pos > 0) THEN
!-----
      found(it) = .TRUE.
      fileorig = fromfile(pos)
!-----
!---- DECODE
!-----
      str_READ = ADJUSTL(fichier(pos))
      str_READ_lower = str_READ
      CALL strlowercase (str_READ_lower)
!-----
      IF (    (TRIM(str_READ_lower) == 'def')     &
 &        .OR.(TRIM(str_READ_lower) == 'default') ) THEN
        def_beha = .TRUE.
      ELSE
        def_beha = .FALSE.
        len_str = LEN_TRIM(str_READ)
        io_err = 0
        SELECT CASE (k_typ)
        CASE(k_i)
          WRITE (UNIT=c_fmt,FMT='("(I",I3.3,")")') len_str
          READ (UNIT=str_READ(1:len_str), &
 &              FMT=c_fmt,IOSTAT=io_err) i_val(it)
        CASE(k_r)
          READ (UNIT=str_READ(1:len_str), &
 &              FMT=*,IOSTAT=io_err) r_val(it)
        CASE(k_c)
          c_val(it) = str_READ(1:len_str)
        CASE(k_l)
          ipos_tr = -1
          ipos_fl = -1
          ipos_tr = MAX(INDEX(str_READ_lower,'tru'), &
 &                      INDEX(str_READ_lower,'y'))
          ipos_fl = MAX(INDEX(str_READ_lower,'fal'), &
 &                      INDEX(str_READ_lower,'n'))
          IF (ipos_tr > 0) THEN
            l_val(it) = .TRUE.
          ELSE IF (ipos_fl > 0) THEN
            l_val(it) = .FALSE.
          ELSE
            io_err = 100
          ENDIF
        END SELECT
        IF (io_err /= 0) THEN
          CALL ipslerr (3,'get_fil', &
 &         'Target '//TRIM(target), &
 &         'is not of '//TRIM(c_vtyp)//' type',' ')
        ENDIF
      ENDIF
!-----
      IF ( (k_typ == k_i).OR.(k_typ == k_r) ) THEN
!-------
!------ Is this the value of a compressed field ?
        compressed = (compline(pos) > 0)
        IF (compressed) THEN
          IF (compline(pos) /= nb_to_ret) THEN
            CALL ipslerr (2,'get_fil', &
 &           'For key '//TRIM(target)//' we have a compressed field', &
 &           'which does not have the right size.', &
 &           'We will try to fix that.')
          ENDIF
          IF      (k_typ == k_i) THEN
            i_cmpval = i_val(it)
          ELSE IF (k_typ == k_r) THEN
            r_cmpval = r_val(it)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      found(it) = .FALSE.
      def_beha = .FALSE.
      compressed = .FALSE.
    ENDIF
  ENDDO
!-
  IF ( (k_typ == k_i).OR.(k_typ == k_r) ) THEN
!---
!-- If this is a compressed field then we will uncompress it
    IF (compressed) THEN
      DO it=1,nb_to_ret
        IF (.NOT.found(it)) THEN
          IF      (k_typ == k_i) THEN
            i_val(it) = i_cmpval
          ELSE IF (k_typ == k_r) THEN
          ENDIF
          found(it) = .TRUE.
        ENDIF
      ENDDO
    ENDIF
  ENDIF
!-
! Now we set the status for what we found
  IF (def_beha) THEN
    status = 2
    WRITE(*,*) 'USING DEFAULT BEHAVIOUR FOR ',TRIM(target)
  ELSE
    status_cnt = 0
    DO it=1,nb_to_ret
      IF (.NOT.found(it)) THEN
        status_cnt = status_cnt+1
        IF      (status_cnt <= max_msgs) THEN
          WRITE (UNIT=*,FMT='(" USING DEFAULTS : ",A)', &
 &               ADVANCE='NO') TRIM(target)
          IF (nb_to_ret > 1) THEN
            WRITE (UNIT=*,FMT='("__")',ADVANCE='NO')
            WRITE (UNIT=*,FMT=c_i_fmt,ADVANCE='NO') it
          ENDIF
          SELECT CASE (k_typ)
          CASE(k_i)
            WRITE (UNIT=*,FMT=*) "=",i_val(it)
          CASE(k_r)
            WRITE (UNIT=*,FMT=*) "=",r_val(it)
          CASE(k_c)
            WRITE (UNIT=*,FMT=*) "=",c_val(it)
          CASE(k_l)
            WRITE (UNIT=*,FMT=*) "=",l_val(it)
          END SELECT
        ELSE IF (status_cnt == max_msgs+1) THEN
          WRITE (UNIT=*,FMT='(" USING DEFAULTS ... ",A)')
        ENDIF
      ENDIF
    ENDDO
!---
    IF (status_cnt == 0) THEN
      status = 1
    ELSE IF (status_cnt == nb_to_ret) THEN
      status = 2
    ELSE
      status = 3
    ENDIF
  ENDIF
! Deallocate the memory
  DEALLOCATE(found)
!---------------------
END SUBROUTINE get_fil
!===
SUBROUTINE get_rdb (pos,size_of_in,target,i_val,r_val,c_val,l_val)
!---------------------------------------------------------------------
!- Read the required variable in the database
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: pos,size_of_in
  CHARACTER(LEN=*) :: target
  INTEGER,DIMENSION(:),OPTIONAL          :: i_val
  REAL,DIMENSION(:),OPTIONAL             :: r_val
  LOGICAL,DIMENSION(:),OPTIONAL          :: l_val
  CHARACTER(LEN=*),DIMENSION(:),OPTIONAL :: c_val
!-
  INTEGER :: k_typ,k_beg,k_end
  CHARACTER(LEN=9) :: c_vtyp
!---------------------------------------------------------------------
!-
! Get the type of the argument
  CALL get_qtyp (k_typ,c_vtyp,i_val,r_val,c_val,l_val)
  IF (     (k_typ /= k_i).AND.(k_typ /= k_r) &
 &    .AND.(k_typ /= k_c).AND.(k_typ /= k_l) )THEN
    CALL ipslerr (3,'get_rdb', &
 &   'Internal error','Unknown type of data',' ')
  ENDIF
!-
  IF (key_tab(pos)%keytype /= k_typ) THEN
    CALL ipslerr (3,'get_rdb', &
 &   'Wrong data type for keyword '//TRIM(target), &
 &   '(NOT '//TRIM(c_vtyp)//')',' ')
  ENDIF
!-
  IF (key_tab(pos)%keycompress > 0) THEN
    IF (    (key_tab(pos)%keycompress /= size_of_in) &
 &      .OR.(key_tab(pos)%keymemlen /= 1) ) THEN
      CALL ipslerr (3,'get_rdb', &
 &     'Wrong compression length','for keyword '//TRIM(target),' ')
    ELSE
      SELECT CASE (k_typ)
      CASE(k_i)
        i_val(1:size_of_in) = i_mem(key_tab(pos)%keymemstart)
      CASE(k_r)
        r_val(1:size_of_in) = r_mem(key_tab(pos)%keymemstart)
      END SELECT
    ENDIF
  ELSE
    IF (key_tab(pos)%keymemlen /= size_of_in) THEN
      CALL ipslerr (3,'get_rdb', &
 &     'Wrong array length','for keyword '//TRIM(target),' ')
    ELSE
      k_beg = key_tab(pos)%keymemstart
      k_end = k_beg+key_tab(pos)%keymemlen-1
      SELECT CASE (k_typ)
      CASE(k_i)
        i_val(1:size_of_in) = i_mem(k_beg:k_end)
      CASE(k_r)
        r_val(1:size_of_in) = r_mem(k_beg:k_end)
      CASE(k_c)
        c_val(1:size_of_in) = c_mem(k_beg:k_end)
      CASE(k_l)
        l_val(1:size_of_in) = l_mem(k_beg:k_end)
      END SELECT
    ENDIF
  ENDIF
!---------------------
END SUBROUTINE get_rdb
!===
SUBROUTINE get_wdb &
 &  (target,status,fileorig,size_of_in, &
 &   i_val,r_val,c_val,l_val)
!---------------------------------------------------------------------
!- Write data into the data base
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: target
  INTEGER :: status,fileorig,size_of_in
  INTEGER,DIMENSION(:),OPTIONAL          :: i_val
  REAL,DIMENSION(:),OPTIONAL             :: r_val
  LOGICAL,DIMENSION(:),OPTIONAL          :: l_val
  CHARACTER(LEN=*),DIMENSION(:),OPTIONAL :: c_val
!-
  INTEGER :: k_typ
  CHARACTER(LEN=9) :: c_vtyp
  INTEGER :: k_mempos,k_memsize,k_beg,k_end
  LOGICAL :: l_cmp
!---------------------------------------------------------------------
!-
! Get the type of the argument
  CALL get_qtyp (k_typ,c_vtyp,i_val,r_val,c_val,l_val)
  IF (     (k_typ /= k_i).AND.(k_typ /= k_r) &
 &    .AND.(k_typ /= k_c).AND.(k_typ /= k_l) )THEN
    CALL ipslerr (3,'get_wdb', &
 &   'Internal error','Unknown type of data',' ')
  ENDIF
!-
! First check if we have sufficiant space for the new key
  IF (nb_keys+1 > keymemsize) THEN
    CALL getin_allockeys ()
  ENDIF
!-
  SELECT CASE (k_typ)
  CASE(k_i)
    k_mempos = i_mempos; k_memsize = i_memsize;
    l_cmp = (MINVAL(i_val) == MAXVAL(i_val)) &
 &         .AND.(size_of_in > compress_lim)
  CASE(k_r)
    k_mempos = r_mempos; k_memsize = r_memsize;
    l_cmp = (MINVAL(r_val) == MAXVAL(r_val)) &
 &         .AND.(size_of_in > compress_lim)
  CASE(k_c)
    k_mempos = c_mempos; k_memsize = c_memsize;
    l_cmp = .FALSE.
  CASE(k_l)
    k_mempos = l_mempos; k_memsize = l_memsize;
    l_cmp = .FALSE.
  END SELECT
!-
! Fill out the items of the data base
  nb_keys = nb_keys+1
  key_tab(nb_keys)%keystr = target(1:MIN(LEN_TRIM(target),l_n))
  key_tab(nb_keys)%keystatus = status
  key_tab(nb_keys)%keytype = k_typ
  key_tab(nb_keys)%keyfromfile = fileorig
  key_tab(nb_keys)%keymemstart = k_mempos+1
  IF (l_cmp) THEN
    key_tab(nb_keys)%keycompress = size_of_in
    key_tab(nb_keys)%keymemlen = 1
  ELSE
    key_tab(nb_keys)%keycompress = -1
    key_tab(nb_keys)%keymemlen = size_of_in
  ENDIF
!-
! Before writing the actual size lets see if we have the space
  IF (key_tab(nb_keys)%keymemstart+key_tab(nb_keys)%keymemlen &
 &    > k_memsize) THEN
    CALL getin_allocmem (k_typ,key_tab(nb_keys)%keymemlen)
  ENDIF
!-
  k_beg = key_tab(nb_keys)%keymemstart
  k_end = k_beg+key_tab(nb_keys)%keymemlen-1
  SELECT CASE (k_typ)
  CASE(k_i)
    i_mem(k_beg:k_end) = i_val(1:key_tab(nb_keys)%keymemlen)
    i_mempos = k_end
  CASE(k_r)
    r_mem(k_beg:k_end) = r_val(1:key_tab(nb_keys)%keymemlen)
    r_mempos = k_end
  CASE(k_c)
    c_mem(k_beg:k_end) = c_val(1:key_tab(nb_keys)%keymemlen)
    c_mempos = k_end
  CASE(k_l)
    l_mem(k_beg:k_end) = l_val(1:key_tab(nb_keys)%keymemlen)
    l_mempos = k_end
  END SELECT
!---------------------
END SUBROUTINE get_wdb
!-
!===
!-
SUBROUTINE getin_read
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,SAVE :: current
!---------------------------------------------------------------------
  IF (allread == 0) THEN
!-- Allocate a first set of memory.
    CALL getin_alloctxt ()
    CALL getin_allockeys ()
    CALL getin_allocmem (k_i,0)
    CALL getin_allocmem (k_r,0)
    CALL getin_allocmem (k_c,0)
    CALL getin_allocmem (k_l,0)
!-- Start with reading the files
    nbfiles = 1
    filelist(1) = TRIM(def_file)
    current = 1
!--
    DO WHILE (current <= nbfiles)
      CALL getin_readdef (current)
      current = current+1
    ENDDO
    allread = 1
    CALL getin_checkcohe ()
  ENDIF
!------------------------
END SUBROUTINE getin_read
!-
!===
!-
  SUBROUTINE getin_readdef(current)
!---------------------------------------------------------------------
!- This subroutine will read the files and only keep the
!- the relevant information. The information is kept as it
!- found in the file. The data will be analysed later.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: current
!-
  CHARACTER(LEN=100) :: READ_str,NEW_str,last_key,key_str
  CHARACTER(LEN=n_d_fmt) :: cnt
  CHARACTER(LEN=10) :: c_fmt
  INTEGER :: nb_lastkey
!-
  INTEGER :: eof,ptn,len_str,i,it,iund,io_err
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
  eof = 0
  ptn = 1
  nb_lastkey = 0
!-
  IF (check) THEN
    WRITE(*,*) 'getin_readdef : Open file ',TRIM(filelist(current))
  ENDIF
!-
  OPEN (UNIT=22,FILE=filelist(current),STATUS="OLD",IOSTAT=io_err)
  IF (io_err /= 0) THEN
    CALL ipslerr (2,'getin_readdef', &
 &  'Could not open file '//TRIM(filelist(current)),' ',' ')
    RETURN
  ENDIF
!-
  DO WHILE (eof /= 1)
!---
    CALL getin_skipafew (22,READ_str,eof,nb_lastkey)
    len_str = LEN_TRIM(READ_str)
    ptn = INDEX(READ_str,'=')
!---
    IF (ptn > 0) THEN
!---- Get the target
      key_str = TRIM(ADJUSTL(READ_str(1:ptn-1)))
!---- Make sure that a vector keyword has the right length
      iund = INDEX(key_str,'__')
      IF (iund > 0) THEN
        WRITE (UNIT=c_fmt,FMT='("(I",I3.3,")")') &
 &        LEN_TRIM(key_str)-iund-1
        READ(UNIT=key_str(iund+2:LEN_TRIM(key_str)), &
 &           FMT=c_fmt,IOSTAT=io_err) it
        IF ( (io_err == 0).AND.(it > 0) ) THEN
          WRITE(UNIT=cnt,FMT=c_i_fmt) it
          key_str = key_str(1:iund+1)//cnt
        ELSE
          CALL ipslerr (3,'getin_readdef', &
 &         'A very strange key has just been found :', &
 &         TRIM(key_str),' ')
        ENDIF
      ENDIF
!---- Prepare the content
      NEW_str = TRIM(ADJUSTL(READ_str(ptn+1:len_str)))
      CALL nocomma (NEW_str)
      CALL cmpblank (NEW_str)
      NEW_str  = TRIM(ADJUSTL(NEW_str))
      IF (check) THEN
        WRITE(*,*) &
 &        '--> getin_readdef : ',TRIM(key_str),' :: ',TRIM(NEW_str)
      ENDIF
!---- Decypher the content of NEW_str
!-
!---- This has to be a new key word, thus :
      nb_lastkey = 0
!----
      CALL getin_decrypt (current,key_str,NEW_str,last_key,nb_lastkey)
!----
    ELSE IF (len_str > 0) THEN
!---- Prepare the key if we have an old one to which
!---- we will add the line just read
      IF (nb_lastkey > 0) THEN
        iund =  INDEX(last_key,'__')
        IF (iund > 0) THEN
!-------- We only continue a keyword, thus it is easy
          key_str = last_key(1:iund-1)
        ELSE
          IF (nb_lastkey /= 1) THEN
            CALL ipslerr (3,'getin_readdef', &
 &           'We can not have a scalar keyword', &
 &           'and a vector content',' ')
          ENDIF
!-------- The last keyword needs to be transformed into a vector.
          WRITE(UNIT=cnt,FMT=c_i_fmt) 1
          targetlist(nb_lines) = &
 &         last_key(1:MIN(LEN_TRIM(last_key),l_n-n_d_fmt-2))//'__'//cnt
          key_str = last_key(1:LEN_TRIM(last_key))
        ENDIF
      ENDIF
!---- Prepare the content
      NEW_str = TRIM(ADJUSTL(READ_str(1:len_str)))
      CALL getin_decrypt (current,key_str,NEW_str,last_key,nb_lastkey)
    ELSE
!---- If we have an empty line then the keyword finishes
      nb_lastkey = 0
      IF (check) THEN
        WRITE(*,*) 'getin_readdef : Have found an emtpy line '
      ENDIF
    ENDIF
  ENDDO
!-
  CLOSE(UNIT=22)
!-
  IF (check) THEN
    OPEN (UNIT=22,file=TRIM(def_file)//'.test')
    DO i=1,nb_lines
      WRITE(UNIT=22,FMT=*) targetlist(i)," : ",fichier(i)
    ENDDO
    CLOSE(UNIT=22)
  ENDIF
!---------------------------
END SUBROUTINE getin_readdef
!-
!===
!-
SUBROUTINE getin_decrypt(current,key_str,NEW_str,last_key,nb_lastkey)
!---------------------------------------------------------------------
!- This subroutine is going to decypher the line.
!- It essentialy checks how many items are included and
!- it they can be attached to a key.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
! ARGUMENTS
!-
  INTEGER :: current,nb_lastkey
  CHARACTER(LEN=*) :: key_str,NEW_str,last_key
!-
! LOCAL
!-
  INTEGER :: len_str,blk,nbve,starpos
  CHARACTER(LEN=100) :: tmp_str,new_key,mult
  CHARACTER(LEN=n_d_fmt) :: cnt
  CHARACTER(LEN=10) :: c_fmt
!---------------------------------------------------------------------
  len_str = LEN_TRIM(NEW_str)
  blk = INDEX(NEW_str(1:len_str),' ')
  tmp_str = NEW_str(1:len_str)
!-
! If the key is a new file then we take it up. Else
! we save the line and go on.
!-
  IF (INDEX(key_str,'INCLUDEDEF') > 0) THEN
    DO WHILE (blk > 0)
      IF (nbfiles+1 > max_files) THEN
        CALL ipslerr (3,'getin_decrypt', &
 &       'Too many files to include',' ',' ')
      ENDIF
!-----
      nbfiles = nbfiles+1
      filelist(nbfiles) = tmp_str(1:blk)
!-----
      tmp_str = TRIM(ADJUSTL(tmp_str(blk+1:LEN_TRIM(tmp_str))))
      blk = INDEX(tmp_str(1:LEN_TRIM(tmp_str)),' ')
    ENDDO
!---
    IF (nbfiles+1 > max_files) THEN
      CALL ipslerr (3,'getin_decrypt', &
 &     'Too many files to include',' ',' ')
    ENDIF
!---
    nbfiles =  nbfiles+1
    filelist(nbfiles) = TRIM(ADJUSTL(tmp_str))
!---
    last_key = 'INCLUDEDEF'
    nb_lastkey = 1
  ELSE
!-
!-- We are working on a new line of input
!-
    IF (nb_lines+1 > i_txtsize) THEN
      CALL getin_alloctxt ()
    ENDIF
    nb_lines = nb_lines+1
!-
!-- First we solve the issue of conpressed information. Once
!-- this is done all line can be handled in the same way.
!-
    starpos = INDEX(NEW_str(1:len_str),'*')
    IF ( (starpos > 0).AND.(tmp_str(1:1) /= '"') &
 &                    .AND.(tmp_str(1:1) /= "'") ) THEN
!-----
      IF (INDEX(key_str(1:LEN_TRIM(key_str)),'__') > 0) THEN
        CALL ipslerr (3,'getin_decrypt', &
 &       'We can not have a compressed field of values', &
 &       'in a vector notation (TARGET__n).', &
 &       'The key at fault : '//TRIM(key_str))
      ENDIF
!-
!---- Read the multiplied
!-
      mult = TRIM(ADJUSTL(NEW_str(1:starpos-1)))
!---- Construct the new string and its parameters
      NEW_str = TRIM(ADJUSTL(NEW_str(starpos+1:len_str)))
      len_str = LEN_TRIM(NEW_str)
      blk = INDEX(NEW_str(1:len_str),' ')
      IF (blk > 1) THEN
        CALL ipslerr (2,'getin_decrypt', &
 &       'This is a strange behavior','you could report',' ')
      ENDIF
      WRITE (UNIT=c_fmt,FMT='("(I",I5.5,")")') LEN_TRIM(mult)
      READ(UNIT=mult,FMT=c_fmt) compline(nb_lines)
!---
    ELSE
      compline(nb_lines) = -1
    ENDIF
!-
!-- If there is no space wthin the line then the target is a scalar
!-- or the element of a properly written vector.
!-- (ie of the type TARGET__00001)
!-
    IF (    (blk <= 1) &
 &      .OR.(tmp_str(1:1) == '"') &
 &      .OR.(tmp_str(1:1) == "'") ) THEN
!-
      IF (nb_lastkey == 0) THEN
!------ Save info of current keyword as a scalar
!------ if it is not a continuation
        targetlist(nb_lines) = key_str(1:MIN(LEN_TRIM(key_str),l_n))
        last_key = key_str(1:MIN(LEN_TRIM(key_str),l_n))
        nb_lastkey = 1
      ELSE
!------ We are continuing a vector so the keyword needs
!------ to get the underscores
        WRITE(UNIT=cnt,FMT=c_i_fmt) nb_lastkey+1
        targetlist(nb_lines) = &
 &        key_str(1:MIN(LEN_TRIM(key_str),l_n-n_d_fmt-2))//'__'//cnt
        last_key = &
 &        key_str(1:MIN(LEN_TRIM(key_str),l_n-n_d_fmt-2))//'__'//cnt
        nb_lastkey = nb_lastkey+1
      ENDIF
!-----
      fichier(nb_lines) = NEW_str(1:len_str)
      fromfile(nb_lines) = current
    ELSE
!-
!---- If there are blanks whithin the line then we are dealing
!---- with a vector and we need to split it in many entries
!---- with the TARGET__n notation.
!----
!---- Test if the targer is not already a vector target !
!-
      IF (INDEX(TRIM(key_str),'__') > 0) THEN
        CALL ipslerr (3,'getin_decrypt', &
 &       'We have found a mixed vector notation (TARGET__n).', &
 &       'The key at fault : '//TRIM(key_str),' ')
      ENDIF
!-
      nbve = nb_lastkey
      nbve = nbve+1
      WRITE(UNIT=cnt,FMT=c_i_fmt) nbve
!-
      DO WHILE (blk > 0)
!-
!------ Save the content of target__nbve
!-
        fichier(nb_lines) = tmp_str(1:blk)
        new_key = &
 &       key_str(1:MIN(LEN_TRIM(key_str),l_n-n_d_fmt-2))//'__'//cnt
        targetlist(nb_lines) = new_key(1:MIN(LEN_TRIM(new_key),l_n))
        fromfile(nb_lines) = current
!-
        tmp_str = TRIM(ADJUSTL(tmp_str(blk+1:LEN_TRIM(tmp_str))))
        blk = INDEX(TRIM(tmp_str),' ')
!-
        IF (nb_lines+1 > i_txtsize) THEN
          CALL getin_alloctxt ()
        ENDIF
        nb_lines = nb_lines+1
        nbve = nbve+1
        WRITE(UNIT=cnt,FMT=c_i_fmt) nbve
!-
      ENDDO
!-
!---- Save the content of the last target
!-
      fichier(nb_lines) = tmp_str(1:LEN_TRIM(tmp_str))
      new_key = &
 &      key_str(1:MIN(LEN_TRIM(key_str),l_n-n_d_fmt-2))//'__'//cnt
      targetlist(nb_lines) = new_key(1:MIN(LEN_TRIM(new_key),l_n))
      fromfile(nb_lines) = current
!-
      last_key = &
 &      key_str(1:MIN(LEN_TRIM(key_str),l_n-n_d_fmt-2))//'__'//cnt
      nb_lastkey = nbve
!-
    ENDIF
!-
  ENDIF
!---------------------------
END SUBROUTINE getin_decrypt
!-
!===
!-
SUBROUTINE getin_checkcohe ()
!---------------------------------------------------------------------
!- This subroutine checks for redundancies.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: line,n_k,k
!---------------------------------------------------------------------
  DO line=1,nb_lines-1
!-
    n_k = 0
    DO k=line+1,nb_lines
      IF (TRIM(targetlist(line)) == TRIM(targetlist(k))) THEN
        n_k = k
        EXIT
      ENDIF
    ENDDO
!---
!-- IF we have found it we have a problem to solve.
!---
    IF (n_k > 0) THEN
      WRITE(*,*) 'COUNT : ',n_k
      WRITE(*,*) &
 &  'getin_checkcohe : Found a problem on key ',TRIM(targetlist(line))
      WRITE(*,*) &
 &  'getin_checkcohe : The following values were encoutered :'
      WRITE(*,*) &
 &  '                ',TRIM(targetlist(line)),' == ',fichier(line)
      WRITE(*,*) &
 &  '                ',TRIM(targetlist(k)),' == ',fichier(k)
      WRITE(*,*) &
 &  'getin_checkcohe : We will keep only the last value'
      targetlist(line) = ' '
    ENDIF
  ENDDO
!-----------------------------
END SUBROUTINE getin_checkcohe
!-
!===
!-
SUBROUTINE getin_skipafew (unit,out_string,eof,nb_lastkey)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: unit,eof,nb_lastkey
  CHARACTER(LEN=100) :: dummy
  CHARACTER(LEN=100) :: out_string
  CHARACTER(LEN=1) :: first
!---------------------------------------------------------------------
  first="#"
  eof = 0
  out_string = "    "
!-
  DO WHILE (first == "#")
    READ (UNIT=unit,FMT='(A)',ERR=9998,END=7778) dummy
    dummy = TRIM(ADJUSTL(dummy))
    first=dummy(1:1)
    IF (first == "#") THEN
      nb_lastkey = 0
    ENDIF
  ENDDO
  out_string=dummy
!-
  RETURN
!-
9998 CONTINUE
  CALL ipslerr (3,'getin_skipafew','Error while reading file',' ',' ')
!-
7778 CONTINUE
  eof = 1
!----------------------------
END SUBROUTINE getin_skipafew
!-
!===
!-
SUBROUTINE getin_allockeys ()
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  TYPE(t_key),ALLOCATABLE,DIMENSION(:) :: tmp_key_tab
!-
  INTEGER :: ier
  CHARACTER(LEN=20) :: c_tmp
!---------------------------------------------------------------------
  IF (keymemsize == 0) THEN
!---
!-- Nothing exists in memory arrays and it is easy to do.
!---
    WRITE (UNIT=c_tmp,FMT=*) memslabs
    ALLOCATE(key_tab(memslabs),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_allockeys', &
 &     'Can not allocate key_tab', &
 &     'to size '//TRIM(ADJUSTL(c_tmp)),' ')
    ENDIF
    nb_keys = 0
    keymemsize = memslabs
    key_tab(:)%keycompress = -1
!---
  ELSE
!---
!-- There is something already in the memory,
!-- we need to transfer and reallocate.
!---
    WRITE (UNIT=c_tmp,FMT=*) keymemsize
    ALLOCATE(tmp_key_tab(keymemsize),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_allockeys', &
 &     'Can not allocate tmp_key_tab', &
 &     'to size '//TRIM(ADJUSTL(c_tmp)),' ')
    ENDIF
    WRITE (UNIT=c_tmp,FMT=*) keymemsize+memslabs
    tmp_key_tab(1:keymemsize) = key_tab(1:keymemsize)
    DEALLOCATE(key_tab)
    ALLOCATE(key_tab(keymemsize+memslabs),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_allockeys', &
 &     'Can not allocate key_tab', &
 &     'to size '//TRIM(ADJUSTL(c_tmp)),' ')
    ENDIF
    key_tab(:)%keycompress = -1
    key_tab(1:keymemsize) = tmp_key_tab(1:keymemsize)
    DEALLOCATE(tmp_key_tab)
    keymemsize = keymemsize+memslabs
  ENDIF
!-----------------------------
END SUBROUTINE getin_allockeys
!-
!===
!-
SUBROUTINE getin_allocmem (type,len_wanted)
!---------------------------------------------------------------------
!- Allocate the memory of the data base for all 4 types of memory
!- INTEGER / REAL / CHARACTER / LOGICAL
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: type,len_wanted
!-
  INTEGER,ALLOCATABLE :: tmp_int(:)
  REAL,ALLOCATABLE :: tmp_real(:)
  CHARACTER(LEN=100),ALLOCATABLE :: tmp_char(:)
  LOGICAL,ALLOCATABLE :: tmp_logic(:)
  INTEGER :: ier
  CHARACTER(LEN=20) :: c_tmp
!---------------------------------------------------------------------
  SELECT CASE (type)
  CASE(k_i)
    IF (i_memsize == 0) THEN
      ALLOCATE(i_mem(memslabs),stat=ier)
      IF (ier /= 0) THEN
        WRITE (UNIT=c_tmp,FMT=*) memslabs
        CALL ipslerr (3,'getin_allocmem', &
 &       'Unable to allocate db-memory', &
 &       'i_mem to size '//TRIM(ADJUSTL(c_tmp)),' ')
      ENDIF
      i_memsize=memslabs
    ELSE
      ALLOCATE(tmp_int(i_memsize),stat=ier)
      IF (ier /= 0) THEN
        WRITE (UNIT=c_tmp,FMT=*) i_memsize
        CALL ipslerr (3,'getin_allocmem', &
 &       'Unable to allocate tmp_int', &
 &       'to size '//TRIM(ADJUSTL(c_tmp)),' ')
      ENDIF
      tmp_int(1:i_memsize) = i_mem(1:i_memsize)
      DEALLOCATE(i_mem)
      ALLOCATE(i_mem(i_memsize+MAX(memslabs,len_wanted)),stat=ier)
      IF (ier /= 0) THEN
        WRITE (UNIT=c_tmp,FMT=*) i_memsize+MAX(memslabs,len_wanted)
        CALL ipslerr (3,'getin_allocmem', &
 &       'Unable to re-allocate db-memory', &
 &       'i_mem to size '//TRIM(ADJUSTL(c_tmp)),' ')
      ENDIF
      i_mem(1:i_memsize) = tmp_int(1:i_memsize)
      i_memsize = i_memsize+MAX(memslabs,len_wanted)
      DEALLOCATE(tmp_int)
    ENDIF
  CASE(k_r)
    IF (r_memsize == 0) THEN
      ALLOCATE(r_mem(memslabs),stat=ier)
      IF (ier /= 0) THEN
        WRITE (UNIT=c_tmp,FMT=*) memslabs
        CALL ipslerr (3,'getin_allocmem', &
 &       'Unable to allocate db-memory', &
 &       'r_mem to size '//TRIM(ADJUSTL(c_tmp)),' ')
      ENDIF
      r_memsize =  memslabs
    ELSE
      ALLOCATE(tmp_real(r_memsize),stat=ier)
      IF (ier /= 0) THEN
        WRITE (UNIT=c_tmp,FMT=*) r_memsize
        CALL ipslerr (3,'getin_allocmem', &
 &       'Unable to allocate tmp_real', &
 &       'to size '//TRIM(ADJUSTL(c_tmp)),' ')
      ENDIF
      tmp_real(1:r_memsize) = r_mem(1:r_memsize)
      DEALLOCATE(r_mem)
      ALLOCATE(r_mem(r_memsize+MAX(memslabs,len_wanted)),stat=ier)
      IF (ier /= 0) THEN
        WRITE (UNIT=c_tmp,FMT=*) r_memsize+MAX(memslabs,len_wanted)
        CALL ipslerr (3,'getin_allocmem', &
 &       'Unable to re-allocate db-memory', &
 &       'r_mem to size '//TRIM(ADJUSTL(c_tmp)),' ')
      ENDIF
      r_mem(1:r_memsize) = tmp_real(1:r_memsize)
      r_memsize = r_memsize+MAX(memslabs,len_wanted)
      DEALLOCATE(tmp_real)
    ENDIF
  CASE(k_c)
    IF (c_memsize == 0) THEN
      ALLOCATE(c_mem(memslabs),stat=ier)
      IF (ier /= 0) THEN
        WRITE (UNIT=c_tmp,FMT=*) memslabs
        CALL ipslerr (3,'getin_allocmem', &
 &       'Unable to allocate db-memory', &
 &       'c_mem to size '//TRIM(ADJUSTL(c_tmp)),' ')
      ENDIF
      c_memsize = memslabs
    ELSE
      ALLOCATE(tmp_char(c_memsize),stat=ier)
      IF (ier /= 0) THEN
        WRITE (UNIT=c_tmp,FMT=*) c_memsize
        CALL ipslerr (3,'getin_allocmem', &
 &       'Unable to allocate tmp_char', &
 &       'to size '//TRIM(ADJUSTL(c_tmp)),' ')
      ENDIF
      tmp_char(1:c_memsize) = c_mem(1:c_memsize)
      DEALLOCATE(c_mem)
      ALLOCATE(c_mem(c_memsize+MAX(memslabs,len_wanted)),stat=ier)
      IF (ier /= 0) THEN
        WRITE (UNIT=c_tmp,FMT=*) c_memsize+MAX(memslabs,len_wanted)
        CALL ipslerr (3,'getin_allocmem', &
 &       'Unable to re-allocate db-memory', &
 &       'c_mem to size '//TRIM(ADJUSTL(c_tmp)),' ')
      ENDIF
      c_mem(1:c_memsize) = tmp_char(1:c_memsize)
      c_memsize = c_memsize+MAX(memslabs,len_wanted)
      DEALLOCATE(tmp_char)
    ENDIF
  CASE(k_l)
    IF (l_memsize == 0) THEN
      ALLOCATE(l_mem(memslabs),stat=ier)
      IF (ier /= 0) THEN
        WRITE (UNIT=c_tmp,FMT=*) memslabs
        CALL ipslerr (3,'getin_allocmem', &
 &       'Unable to allocate db-memory', &
 &       'l_mem to size '//TRIM(ADJUSTL(c_tmp)),' ')
      ENDIF
      l_memsize = memslabs
    ELSE
      ALLOCATE(tmp_logic(l_memsize),stat=ier)
      IF (ier /= 0) THEN
        WRITE (UNIT=c_tmp,FMT=*) l_memsize
        CALL ipslerr (3,'getin_allocmem', &
 &       'Unable to allocate tmp_logic', &
 &       'to size '//TRIM(ADJUSTL(c_tmp)),' ')
      ENDIF
      tmp_logic(1:l_memsize) = l_mem(1:l_memsize)
      DEALLOCATE(l_mem)
      ALLOCATE(l_mem(l_memsize+MAX(memslabs,len_wanted)),stat=ier)
      IF (ier /= 0) THEN
        WRITE (UNIT=c_tmp,FMT=*) l_memsize+MAX(memslabs,len_wanted)
        CALL ipslerr (3,'getin_allocmem', &
 &       'Unable to re-allocate db-memory', &
 &       'l_mem to size '//TRIM(ADJUSTL(c_tmp)),' ')
      ENDIF
      l_mem(1:l_memsize) = tmp_logic(1:l_memsize)
      l_memsize = l_memsize+MAX(memslabs,len_wanted)
      DEALLOCATE(tmp_logic)
    ENDIF
  CASE DEFAULT
    CALL ipslerr (3,'getin_allocmem','Unknown type of data',' ',' ')
  END SELECT
!----------------------------
END SUBROUTINE getin_allocmem
!-
!===
!-
SUBROUTINE getin_alloctxt ()
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=100),ALLOCATABLE :: tmp_fic(:)
  CHARACTER(LEN=l_n),ALLOCATABLE :: tmp_tgl(:)
  INTEGER,ALLOCATABLE :: tmp_int(:)
!-
  INTEGER :: ier
  CHARACTER(LEN=20) :: c_tmp1,c_tmp2
!---------------------------------------------------------------------
  IF (i_txtsize == 0) THEN
!---
!-- Nothing exists in memory arrays and it is easy to do.
!---
    WRITE (UNIT=c_tmp1,FMT=*) i_txtslab
    ALLOCATE(fichier(i_txtslab),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_alloctxt', &
 &     'Can not allocate fichier', &
 &     'to size '//TRIM(ADJUSTL(c_tmp1)),' ')
    ENDIF
!---
    ALLOCATE(targetlist(i_txtslab),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_alloctxt', &
 &     'Can not allocate targetlist', &
 &     'to size '//TRIM(ADJUSTL(c_tmp1)),' ')
    ENDIF
!---
    ALLOCATE(fromfile(i_txtslab),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_alloctxt', &
 &     'Can not allocate fromfile', &
 &     'to size '//TRIM(ADJUSTL(c_tmp1)),' ')
    ENDIF
!---
    ALLOCATE(compline(i_txtslab),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_alloctxt', &
 &     'Can not allocate compline', &
 &     'to size '//TRIM(ADJUSTL(c_tmp1)),' ')
    ENDIF
!---
    nb_lines = 0
    i_txtsize = i_txtslab
  ELSE
!---
!-- There is something already in the memory,
!-- we need to transfer and reallocate.
!---
    WRITE (UNIT=c_tmp1,FMT=*) i_txtsize
    WRITE (UNIT=c_tmp2,FMT=*) i_txtsize+i_txtslab
    ALLOCATE(tmp_fic(i_txtsize),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_alloctxt', &
 &     'Can not allocate tmp_fic', &
 &     'to size '//TRIM(ADJUSTL(c_tmp1)),' ')
    ENDIF
    tmp_fic(1:i_txtsize) = fichier(1:i_txtsize)
    DEALLOCATE(fichier)
    ALLOCATE(fichier(i_txtsize+i_txtslab),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_alloctxt', &
 &     'Can not allocate fichier', &
 &     'to size '//TRIM(ADJUSTL(c_tmp2)),' ')
    ENDIF
    fichier(1:i_txtsize) = tmp_fic(1:i_txtsize)
    DEALLOCATE(tmp_fic)
!---
    ALLOCATE(tmp_tgl(i_txtsize),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_alloctxt', &
 &     'Can not allocate tmp_tgl', &
 &     'to size '//TRIM(ADJUSTL(c_tmp1)),' ')
    ENDIF
    tmp_tgl(1:i_txtsize) = targetlist(1:i_txtsize)
    DEALLOCATE(targetlist)
    ALLOCATE(targetlist(i_txtsize+i_txtslab),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_alloctxt', &
 &     'Can not allocate targetlist', &
 &     'to size '//TRIM(ADJUSTL(c_tmp2)),' ')
    ENDIF
    targetlist(1:i_txtsize) = tmp_tgl(1:i_txtsize)
    DEALLOCATE(tmp_tgl)
!---
    ALLOCATE(tmp_int(i_txtsize),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_alloctxt', &
 &     'Can not allocate tmp_int', &
 &     'to size '//TRIM(ADJUSTL(c_tmp1)),' ')
    ENDIF
    tmp_int(1:i_txtsize) = fromfile(1:i_txtsize)
    DEALLOCATE(fromfile)
    ALLOCATE(fromfile(i_txtsize+i_txtslab),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_alloctxt', &
 &     'Can not allocate fromfile', &
 &     'to size '//TRIM(ADJUSTL(c_tmp2)),' ')
    ENDIF
    fromfile(1:i_txtsize) = tmp_int(1:i_txtsize)
!---
    tmp_int(1:i_txtsize) = compline(1:i_txtsize)
    DEALLOCATE(compline)
    ALLOCATE(compline(i_txtsize+i_txtslab),stat=ier)
    IF (ier /= 0) THEN
      CALL ipslerr (3,'getin_alloctxt', &
 &     'Can not allocate compline', &
 &     'to size '//TRIM(ADJUSTL(c_tmp2)),' ')
    ENDIF
    compline(1:i_txtsize) = tmp_int(1:i_txtsize)
    DEALLOCATE(tmp_int)
!---
    i_txtsize = i_txtsize+i_txtslab
  ENDIF
!----------------------------
END SUBROUTINE getin_alloctxt
!-
!===
!-
SUBROUTINE getin_dump (fileprefix)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(*),OPTIONAL :: fileprefix
!-
  CHARACTER(LEN=80) :: usedfileprefix
  INTEGER :: ikey,if,iff,iv
  CHARACTER(LEN=20) :: c_tmp
  CHARACTER(LEN=100) :: tmp_str,used_filename
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
  IF (PRESENT(fileprefix)) THEN
    usedfileprefix = fileprefix(1:MIN(LEN_TRIM(fileprefix),80))
  ELSE
    usedfileprefix = "used"
  ENDIF
!-
  DO if=1,nbfiles
!---
    used_filename = TRIM(usedfileprefix)//'_'//TRIM(filelist(if))
    IF (check) THEN
      WRITE(*,*) &
 &      'GETIN_DUMP : opens file : ',TRIM(used_filename),' if = ',if
      WRITE(*,*) 'GETIN_DUMP : NUMBER OF KEYS : ',nb_keys
    ENDIF
    OPEN (UNIT=22,FILE=used_filename)
!---
!-- If this is the first file we need to add the list
!-- of file which belong to it
    IF ( (if == 1).AND.(nbfiles > 1) ) THEN
      WRITE(22,*) '# '
      WRITE(22,*) '# This file is linked to the following files :'
      WRITE(22,*) '# '
      DO iff=2,nbfiles
        WRITE(22,*) 'INCLUDEDEF = ',TRIM(filelist(iff))
      ENDDO
      WRITE(22,*) '# '
    ENDIF
!---
    DO ikey=1,nb_keys
!-----
!---- Is this key from this file ?
      IF (key_tab(ikey)%keyfromfile == if) THEN
!-------
!------ Write some comments
        WRITE(22,*) '#'
        SELECT CASE (key_tab(ikey)%keystatus)
        CASE(1)
          WRITE(22,*) '# Values of ', &
 &          TRIM(key_tab(ikey)%keystr),' comes from ',TRIM(def_file)
        CASE(2)
          WRITE(22,*) '# Values of ', &
 &          TRIM(key_tab(ikey)%keystr),' are all defaults.'
        CASE(3)
          WRITE(22,*) '# Values of ', &
 &          TRIM(key_tab(ikey)%keystr), &
 &          ' are a mix of ',TRIM(def_file),' and defaults.'
        CASE DEFAULT
          WRITE(22,*) '# Dont know from where the value of ', &
 &          TRIM(key_tab(ikey)%keystr),' comes.'
        END SELECT
        WRITE(22,*) '#'
!-------
!------ Write the values
        SELECT CASE (key_tab(ikey)%keytype)
        CASE(k_i)
          IF (key_tab(ikey)%keymemlen == 1) THEN
            IF (key_tab(ikey)%keycompress < 0) THEN
              WRITE(22,*) &
 &              TRIM(key_tab(ikey)%keystr), &
 &              ' = ',i_mem(key_tab(ikey)%keymemstart)
            ELSE
              WRITE(22,*) &
 &              TRIM(key_tab(ikey)%keystr), &
 &              ' = ',key_tab(ikey)%keycompress, &
 &              ' * ',i_mem(key_tab(ikey)%keymemstart)
            ENDIF
          ELSE
            DO iv=0,key_tab(ikey)%keymemlen-1
              WRITE(UNIT=c_tmp,FMT=c_i_fmt) iv+1
              WRITE(22,*) &
 &              TRIM(key_tab(ikey)%keystr), &
 &              '__',TRIM(ADJUSTL(c_tmp)), &
 &              ' = ',i_mem(key_tab(ikey)%keymemstart+iv)
            ENDDO
          ENDIF
        CASE(k_r)
          IF (key_tab(ikey)%keymemlen == 1) THEN
            IF (key_tab(ikey)%keycompress < 0) THEN
              WRITE(22,*) &
 &              TRIM(key_tab(ikey)%keystr), &
 &              ' = ',r_mem(key_tab(ikey)%keymemstart)
            ELSE
              WRITE(22,*) &
 &              TRIM(key_tab(ikey)%keystr), &
 &              ' = ',key_tab(ikey)%keycompress, &
                   & ' * ',r_mem(key_tab(ikey)%keymemstart)
            ENDIF
          ELSE
            DO iv=0,key_tab(ikey)%keymemlen-1
              WRITE(UNIT=c_tmp,FMT=c_i_fmt) iv+1
              WRITE(22,*) &
 &              TRIM(key_tab(ikey)%keystr),'__',TRIM(ADJUSTL(c_tmp)), &
 &              ' = ',r_mem(key_tab(ikey)%keymemstart+iv)
            ENDDO
          ENDIF
        CASE(k_c)
          IF (key_tab(ikey)%keymemlen == 1) THEN
            tmp_str = c_mem(key_tab(ikey)%keymemstart)
            WRITE(22,*) TRIM(key_tab(ikey)%keystr), &
 &              ' = ',TRIM(tmp_str)
          ELSE
            DO iv=0,key_tab(ikey)%keymemlen-1
              WRITE(UNIT=c_tmp,FMT=c_i_fmt) iv+1
              tmp_str = c_mem(key_tab(ikey)%keymemstart+iv)
              WRITE(22,*) &
 &              TRIM(key_tab(ikey)%keystr), &
 &              '__',TRIM(ADJUSTL(c_tmp)), &
 &              ' = ',TRIM(tmp_str)
            ENDDO
          ENDIF
        CASE(k_l)
          IF (key_tab(ikey)%keymemlen == 1) THEN
            IF (l_mem(key_tab(ikey)%keymemstart)) THEN
              WRITE(22,*) TRIM(key_tab(ikey)%keystr),' = TRUE '
            ELSE
              WRITE(22,*) TRIM(key_tab(ikey)%keystr),' = FALSE '
            ENDIF
          ELSE
            DO iv=0,key_tab(ikey)%keymemlen-1
              WRITE(UNIT=c_tmp,FMT=c_i_fmt) iv+1
              IF (l_mem(key_tab(ikey)%keymemstart+iv)) THEN
                WRITE(22,*) TRIM(key_tab(ikey)%keystr),'__', &
 &                          TRIM(ADJUSTL(c_tmp)),' = TRUE '
              ELSE
                WRITE(22,*) TRIM(key_tab(ikey)%keystr),'__', &
 &                          TRIM(ADJUSTL(c_tmp)),' = FALSE '
              ENDIF
            ENDDO
          ENDIF
        CASE DEFAULT
          CALL ipslerr (3,'getin_dump', &
 &         'Unknown type for variable '//TRIM(key_tab(ikey)%keystr), &
 &         ' ',' ')
        END SELECT
      ENDIF
    ENDDO
!-
    CLOSE(UNIT=22)
!-
  ENDDO
!------------------------
END SUBROUTINE getin_dump
!===
SUBROUTINE get_qtyp (k_typ,c_vtyp,i_v,r_v,c_v,l_v)
!---------------------------------------------------------------------
!- Returns the type of the argument (mutually exclusive)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(OUT) :: k_typ
  CHARACTER(LEN=*),INTENT(OUT) :: c_vtyp
  INTEGER,DIMENSION(:),OPTIONAL          :: i_v
  REAL,DIMENSION(:),OPTIONAL             :: r_v
  LOGICAL,DIMENSION(:),OPTIONAL          :: l_v
  CHARACTER(LEN=*),DIMENSION(:),OPTIONAL :: c_v
!---------------------------------------------------------------------
  k_typ = 0
  IF (COUNT((/PRESENT(i_v),PRESENT(r_v),PRESENT(c_v),PRESENT(l_v)/)) &
 &    /= 1) THEN
    CALL ipslerr (3,'get_qtyp', &
 &   'Invalid number of optional arguments','(/= 1)',' ')
  ENDIF
!-
  IF     (PRESENT(i_v)) THEN
    k_typ = k_i
    c_vtyp = 'INTEGER'
  ELSEIF (PRESENT(r_v)) THEN
    k_typ = k_r
    c_vtyp = 'REAL'
  ELSEIF (PRESENT(c_v)) THEN
    k_typ = k_c
    c_vtyp = 'CHARACTER'
  ELSEIF (PRESENT(l_v)) THEN
    k_typ = k_l
    c_vtyp = 'LOGICAL'
  ENDIF
!----------------------
END SUBROUTINE get_qtyp
!===
SUBROUTINE get_findkey (i_tab,c_key,pos)
!---------------------------------------------------------------------
!- This subroutine looks for a key in a table
!---------------------------------------------------------------------
!- INPUT
!-   i_tab  : 1 -> search in key_tab(1:nb_keys)%keystr
!-            2 -> search in targetlist(1:nb_lines)
!-   c_key  : Name of the key we are looking for
!- OUTPUT
!-   pos    : -1 if key not found, else value in the table
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(in) :: i_tab
  CHARACTER(LEN=*),INTENT(in) :: c_key
  INTEGER,INTENT(out) :: pos
!-
  INTEGER :: ikey_max,ikey
  CHARACTER(LEN=l_n) :: c_q_key
!---------------------------------------------------------------------
  pos = -1
  IF     (i_tab == 1) THEN
    ikey_max = nb_keys
  ELSEIF (i_tab == 2) THEN
    ikey_max = nb_lines
  ELSE
    ikey_max = 0
  ENDIF
  IF ( ikey_max > 0 ) THEN
    DO ikey=1,ikey_max
      IF (i_tab == 1) THEN
        c_q_key = key_tab(ikey)%keystr
      ELSE
        c_q_key = targetlist(ikey)
      ENDIF
      IF (TRIM(c_q_key) == TRIM(c_key)) THEN
        pos = ikey
        EXIT
      ENDIF
    ENDDO
  ENDIF
!-------------------------
END SUBROUTINE get_findkey
!===
!------------------
END MODULE getincom
