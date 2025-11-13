PROGRAM flio_rbld
!
!$Id: flio_rbld.f90 3680 2012-11-27 14:42:24Z rblod $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!!--------------------------------------------------------------------
!! PROGRAM flio_rbld
!!
!! PURPOSE :
!!   Recombine the files of MPI version of IOIPSL
!!   along several dimensions.
!!
!! CALLING SEQUENCE :
!!
!!   "flio_rbld" is usually invoked by the script "rebuild"
!!
!!   rebuild -h
!!
!!   rebuild [-v lev] [-f] -o outfile infile[1] ... infile[n]
!!
!! INPUT for "rebuild" :
!!
!!   -h         : help
!!   -v lev     : verbosity level
!!   -f         : force executing mode
!!   -o outfile : name of the recombined file.
!!   infiles    : names of the files that must be recombined.
!!
!! INPUT for "flio_rbld" :
!!
!!  (I) i_v_lev  : verbosity level
!!  (C) c_force  : executing mode (noforce/force)
!!  (I) f_nb     : total number of files
!!  (C) f_nm(:)  : names of the files (input_files output_file)
!!
!!
!! ASSOCIATED MODULES :
!!   IOIPSL(fliocom)
!!
!! RESTRICTIONS :
!!
!!   Cases for character are not coded.
!!
!!   Cases for netCDF variables such as array with more
!!   than 5 dimensions are not coded.
!!
!!   Input files must have the following global attributes :
!!
!!     "DOMAIN_number_total"
!!     "DOMAIN_number"
!!     "DOMAIN_dimensions_ids"
!!     "DOMAIN_size_global"
!!     "DOMAIN_size_local"
!!     "DOMAIN_position_first"
!!     "DOMAIN_position_last"
!!     "DOMAIN_halo_size_start"
!!     "DOMAIN_halo_size_end"
!!     "DOMAIN_type"
!!
!!   NetCDF files must be smaller than 2 Gb.
!!
!!   Character variables should have less than 257 letters
!!
!! EXAMPLE :
!!
!!   rebuild -v -o sst.nc sst_[0-9][0-9][0-9][0-9].nc
!!
!! MODIFICATION HISTORY :
!!   Sebastien Masson   (smasson@jamstec.go.jp)   March 2004
!!   Jacques   Bellier  (Jacques.Bellier@cea.fr)  June  2005
!!--------------------------------------------------------------------
  USE IOIPSL
  USE defprec
!-
  IMPLICIT NONE
!-
! Character length
  INTEGER,PARAMETER :: chlen=256
!-
! DO loops and test related variables
  INTEGER :: i,ia,id,iv,iw,i_i,i_n
  INTEGER :: ik,itmin,itmax,it1,it2,it
  LOGICAL :: l_force,l_uld
!-
! Input arguments related variables
  INTEGER :: i_v_lev
  CHARACTER(LEN=15) :: c_force
  INTEGER :: f_nb,f_nb_in
  CHARACTER(LEN=chlen),DIMENSION(:),ALLOCATABLE :: f_nm
!-
! Domains related variables
  INTEGER :: d_n_t,i_ntd
  INTEGER,DIMENSION(:),ALLOCATABLE :: dom_att,d_d_i,d_s_g
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: d_s_l,d_p_f,d_p_l,d_h_s,d_h_e
  LOGICAL :: l_cgd,l_cof,l_col,l_o_f,l_o_m,l_o_l
  CHARACTER(LEN=chlen) :: c_d_n
!-
! Model files related variables
  LOGICAL :: l_ocf
  INTEGER,DIMENSION(:),ALLOCATABLE :: f_a_id
  INTEGER :: f_id_i1,f_id_i,f_id_o
  INTEGER :: f_d_nb,f_v_nb,f_a_nb,f_d_ul
  INTEGER :: v_a_nb,a_type
  CHARACTER(LEN=chlen),DIMENSION(:),ALLOCATABLE :: &
&  f_d_nm,f_v_nm,f_a_nm,v_a_nm
  CHARACTER(LEN=chlen) :: f_u_nm
  INTEGER,DIMENSION(:),ALLOCATABLE :: v_d_nb,v_d_ul,v_type
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: v_d_i
  INTEGER,DIMENSION(:),ALLOCATABLE :: f_d_i,f_d_l
  INTEGER :: a_l
  INTEGER,DIMENSION(flio_max_var_dims) :: d_i,ib,ie
  INTEGER,DIMENSION(:),ALLOCATABLE :: &
 &  io_i,io_n,ia_sf,io_sf,io_cf,ia_sm,io_sm,io_cm,ia_sl,io_sl,io_cl
  LOGICAL :: l_ex
  CHARACTER(LEN=chlen) :: c_wn1,c_wn2
!-
!?INTEGERS of KIND 1 are not supported on all computers
!?INTEGER(KIND=i_1) :: i1_0d
!?INTEGER(KIND=i_1),DIMENSION(:),ALLOCATABLE :: i1_1d
!?INTEGER(KIND=i_1),DIMENSION(:,:),ALLOCATABLE :: i1_2d
!?INTEGER(KIND=i_1),DIMENSION(:,:,:),ALLOCATABLE :: i1_3d
!?INTEGER(KIND=i_1),DIMENSION(:,:,:,:),ALLOCATABLE :: i1_4d
!?INTEGER(KIND=i_1),DIMENSION(:,:,:,:,:),ALLOCATABLE :: i1_5d
  INTEGER(KIND=i_2) :: i2_0d
  INTEGER(KIND=i_2),DIMENSION(:),ALLOCATABLE :: i2_1d
  INTEGER(KIND=i_2),DIMENSION(:,:),ALLOCATABLE :: i2_2d
  INTEGER(KIND=i_2),DIMENSION(:,:,:),ALLOCATABLE :: i2_3d
  INTEGER(KIND=i_2),DIMENSION(:,:,:,:),ALLOCATABLE :: i2_4d
  INTEGER(KIND=i_2),DIMENSION(:,:,:,:,:),ALLOCATABLE :: i2_5d
  INTEGER(KIND=i_4) :: i4_0d
  INTEGER(KIND=i_4),DIMENSION(:),ALLOCATABLE :: i4_1d
  INTEGER(KIND=i_4),DIMENSION(:,:),ALLOCATABLE :: i4_2d
  INTEGER(KIND=i_4),DIMENSION(:,:,:),ALLOCATABLE :: i4_3d
  INTEGER(KIND=i_4),DIMENSION(:,:,:,:),ALLOCATABLE :: i4_4d
  INTEGER(KIND=i_4),DIMENSION(:,:,:,:,:),ALLOCATABLE :: i4_5d
  REAL(KIND=r_4) :: r4_0d
  REAL(KIND=r_4),DIMENSION(:),ALLOCATABLE :: r4_1d
  REAL(KIND=r_4),DIMENSION(:,:),ALLOCATABLE :: r4_2d
  REAL(KIND=r_4),DIMENSION(:,:,:),ALLOCATABLE :: r4_3d
  REAL(KIND=r_4),DIMENSION(:,:,:,:),ALLOCATABLE :: r4_4d
  REAL(KIND=r_8),DIMENSION(:,:,:,:,:),ALLOCATABLE :: r4_5d
  REAL(KIND=r_8) :: r8_0d
  REAL(KIND=r_8),DIMENSION(:),ALLOCATABLE :: r8_1d
  REAL(KIND=r_8),DIMENSION(:,:),ALLOCATABLE :: r8_2d
  REAL(KIND=r_8),DIMENSION(:,:,:),ALLOCATABLE :: r8_3d
  REAL(KIND=r_8),DIMENSION(:,:,:,:),ALLOCATABLE :: r8_4d
  REAL(KIND=r_8),DIMENSION(:,:,:,:,:),ALLOCATABLE :: r8_5d
!-
! elapsed and cpu time computation variables
  INTEGER :: nb_cc_ini,nb_cc_end,nb_cc_sec,nb_cc_max
  REAL :: t_cpu_ini,t_cpu_end
!---------------------------------------------------------------------
!-
!-------------------
! INPUT arguments
!-------------------
!-
! Retrieve the verbosity level
  READ (UNIT=*,FMT=*) i_v_lev
!-
! Retrieve the executing mode
  READ (UNIT=*,FMT='(A)') c_force
  l_force = (TRIM(c_force)  == 'force')
!-
! Retrieve the number of arguments
  READ (UNIT=*,FMT=*) f_nb
  f_nb_in = f_nb-1
!-
! Retrieve the file names
  ALLOCATE(f_nm(f_nb))
  DO iw=1,f_nb
    READ (UNIT=*,FMT='(A)') f_nm(iw)
  ENDDO
!-
! Allocate and initialize the array of file access identifiers
  ALLOCATE(f_a_id(f_nb_in)); f_a_id(:) = -1;
!-
  IF (i_v_lev >= 1) THEN
    WRITE (UNIT=*,FMT='("")')
    WRITE (UNIT=*,FMT='(" verbosity level : ",I4)') i_v_lev
    WRITE (UNIT=*,FMT='(" executing  mode : ",A)') TRIM(c_force)
    WRITE (UNIT=*,FMT='(" number of args  : ",I4)') f_nb
    WRITE (UNIT=*,FMT='(" Input  files :")')
    DO iw=1,f_nb_in
      WRITE (*,'("   ",A)') TRIM(f_nm(iw))
    ENDDO
    WRITE (UNIT=*,FMT='(" Output file  :")')
    WRITE (*,'("   ",A)') TRIM(f_nm(f_nb))
!-- time initializations
    CALL system_clock &
 &   (count=nb_cc_ini,count_rate=nb_cc_sec,count_max=nb_cc_max)
    CALL cpu_time (t_cpu_ini)
  ENDIF
!-
!---------------------------------------------------
! Retrieve basic informations from the first file
!---------------------------------------------------
!-
! Open the first file
  CALL flrb_of (1,f_id_i)
!-
! Get the attribute "DOMAIN_number_total"
  CALL fliogeta (f_id_i,"?","DOMAIN_number_total",d_n_t)
!-
! Validate the number of input files :
! should be equal to the total number
! of domains used in the simulation
  IF (d_n_t /= f_nb_in) THEN
    IF (l_force) THEN
      iw = 2
    ELSE
      iw = 3
      DEALLOCATE(f_nm,f_a_id)
      CALL flrb_cf (1,.TRUE.)
    ENDIF
    CALL ipslerr (iw,"flio_rbld", &
 &   "The number of input files", &
 &   "is not equal to the number of DOMAINS"," ")
  ENDIF
!-
! Retrieve the basic characteristics of the first input file
  CALL flioinqf &
 & (f_id_i,nb_dim=f_d_nb,nb_var=f_v_nb,nb_gat=f_a_nb,id_uld=f_d_ul)
!-
! Build the list of the names of the
! dimensions/variables/global_attributes and retrieve
! the unlimited_dimension name from the first input file
  ALLOCATE(f_d_nm(f_d_nb),f_v_nm(f_v_nb),f_a_nm(f_a_nb))
  CALL flioinqn (f_id_i,cn_dim=f_d_nm,cn_var=f_v_nm, &
 &                      cn_gat=f_a_nm,cn_uld=f_u_nm)
!-
! Build the list of the dimensions identifiers and lengths
  ALLOCATE(f_d_i(f_d_nb),f_d_l(f_d_nb))
  CALL flioinqf (f_id_i,id_dim=f_d_i,ln_dim=f_d_l)
!-
! Close the file
  CALL flrb_cf (1,.FALSE.)
!-
! Check if the number of needed files is greater than
! the maximum number of simultaneously opened files.
! In that case, open and close model files for each reading,
! otherwise keep the "flio" identifiers of the opened files.
  l_ocf = (f_nb > flio_max_files)
!-
!----------------------------------------------------
! Retrieve domain informations for each input file
!----------------------------------------------------
!-
  DO iw=1,f_nb_in
!---
    CALL flrb_of (iw,f_id_i)
!---
    IF (iw > 1) THEN
      c_wn1 = "DOMAIN_number_total"
      CALL flioinqa (f_id_i,"?",TRIM(c_wn1),l_ex,a_l=a_l)
      IF (l_ex) THEN
        CALL fliogeta (f_id_i,"?",TRIM(c_wn1),i_ntd)
        IF (i_ntd /= d_n_t) THEN
          CALL ipslerr (3,"flio_rbld", &
 &        "File      : "//TRIM(f_nm(iw)), &
 &        "Attribute : "//TRIM(c_wn1), &
 &        "not equal to the one of the first file")
        ENDIF
      ELSE
        CALL ipslerr (3,"flio_rbld", &
 &       "File      : "//TRIM(f_nm(iw)), &
 &       "Attribute : "//TRIM(c_wn1),"not found")
      ENDIF
    ENDIF
!---
    c_wn1 = "DOMAIN_dimensions_ids"
    CALL flioinqa (f_id_i,"?",TRIM(c_wn1),l_ex,a_l=a_l)
    IF (l_ex) THEN
      ALLOCATE(dom_att(a_l))
      CALL fliogeta (f_id_i,"?",TRIM(c_wn1),dom_att)
      IF (iw == 1) THEN
        IF (ANY(dom_att(:) == f_d_ul)) THEN
          CALL ipslerr (3,"flio_rbld", &
 &         "File      : "//TRIM(f_nm(iw)), &
 &         "Attribute : "//TRIM(c_wn1), &
 &         "contains the unlimited dimension")
        ENDIF
        ALLOCATE (d_d_i(a_l))
        d_d_i(:) = dom_att(:)
      ELSEIF (SIZE(dom_att) /= SIZE(d_d_i)) THEN
        CALL ipslerr (3,"flio_rbld", &
 &       "File      : "//TRIM(f_nm(iw)), &
 &       "size of the attribute : "//TRIM(c_wn1), &
 &       "not equal to the one of the first file")
      ELSEIF (ANY(dom_att(:) /= d_d_i(:))) THEN
        CALL ipslerr (3,"flio_rbld", &
 &       "File      : "//TRIM(f_nm(iw)), &
 &       "Attribute : "//TRIM(c_wn1), &
 &       "not equal to the one of the first file")
      ENDIF
      DEALLOCATE(dom_att)
    ELSE
      CALL ipslerr (3,"flio_rbld", &
 &     "File      : "//TRIM(f_nm(iw)), &
 &     "Attribute : "//TRIM(c_wn1),"not found")
    ENDIF
!---
    c_wn1 = "DOMAIN_size_global"
    CALL flioinqa (f_id_i,"?",TRIM(c_wn1),l_ex,a_l=a_l)
    IF (l_ex) THEN
      IF (a_l /= SIZE(d_d_i)) THEN
        CALL ipslerr (3,"flio_rbld", &
 &       "File      : "//TRIM(f_nm(iw)), &
 &       "size of the attribute : "//TRIM(c_wn1), &
 &       "not equal to the size of DOMAIN_dimensions_ids")
      ELSE
        ALLOCATE(dom_att(a_l))
        CALL fliogeta (f_id_i,"?",TRIM(c_wn1),dom_att)
        IF (iw == 1) THEN
          ALLOCATE (d_s_g(a_l))
          d_s_g(:)=dom_att(:)
        ELSEIF (ANY(dom_att(:) /= d_s_g(:))) THEN
          CALL ipslerr (3,"flio_rbld", &
 &         "File      : "//TRIM(f_nm(iw)), &
 &         "Attribute : "//TRIM(c_wn1), &
 &         "not equal to the one of the first file")
        ENDIF
        DEALLOCATE(dom_att)
      ENDIF
    ELSE
      CALL ipslerr (3,"flio_rbld", &
 &     "File      : "//TRIM(f_nm(iw)), &
 &     "Attribute : "//TRIM(c_wn1),"not found")
    ENDIF
!---
    c_wn1 = "DOMAIN_size_local"
    CALL flioinqa (f_id_i,"?",TRIM(c_wn1),l_ex,a_l=a_l)
    IF (l_ex) THEN
      IF (a_l /= SIZE(d_d_i)) THEN
        CALL ipslerr (3,"flio_rbld", &
 &       "File      : "//TRIM(f_nm(iw)), &
 &       "size of the attribute : "//TRIM(c_wn1), &
 &       "not equal to the size of DOMAIN_dimensions_ids")
      ELSE
        ALLOCATE(dom_att(a_l))
        CALL fliogeta (f_id_i,"?",TRIM(c_wn1),dom_att)
        IF (iw == 1) THEN
          ALLOCATE (d_s_l(a_l,f_nb_in))
        ENDIF
        d_s_l(:,iw)=dom_att(:)
        DEALLOCATE(dom_att)
      ENDIF
    ELSE
      CALL ipslerr (3,"flio_rbld", &
 &     "File      : "//TRIM(f_nm(iw)), &
 &     "Attribute : "//TRIM(c_wn1),"not found")
    ENDIF
!---
    c_wn1 = "DOMAIN_position_first"
    CALL flioinqa (f_id_i,"?",TRIM(c_wn1),l_ex,a_l=a_l)
    IF (l_ex) THEN
      IF (a_l /= SIZE(d_d_i)) THEN
        CALL ipslerr (3,"flio_rbld", &
 &       "File      : "//TRIM(f_nm(iw)), &
 &       "size of the attribute : "//TRIM(c_wn1), &
 &       "not equal to the size of DOMAIN_dimensions_ids")
      ELSE
        ALLOCATE(dom_att(a_l))
        CALL fliogeta (f_id_i,"?",TRIM(c_wn1),dom_att)
        IF (iw == 1) THEN
          ALLOCATE (d_p_f(a_l,f_nb_in))
        ENDIF
        d_p_f(:,iw)=dom_att(:)
        DEALLOCATE(dom_att)
      ENDIF
    ELSE
      CALL ipslerr (3,"flio_rbld", &
 &     "File      : "//TRIM(f_nm(iw)), &
 &     "Attribute : "//TRIM(c_wn1),"not found")
    ENDIF
!---
    c_wn1 = "DOMAIN_position_last"
    CALL flioinqa (f_id_i,"?",TRIM(c_wn1),l_ex,a_l=a_l)
    IF (l_ex) THEN
      IF (a_l /= SIZE(d_d_i)) THEN
        CALL ipslerr (3,"flio_rbld", &
 &       "File      : "//TRIM(f_nm(iw)), &
 &       "size of the attribute : "//TRIM(c_wn1), &
 &       "not equal to the size of DOMAIN_dimensions_ids")
      ELSE
        ALLOCATE(dom_att(a_l))
        CALL fliogeta (f_id_i,"?",TRIM(c_wn1),dom_att)
        IF (iw == 1) THEN
          ALLOCATE (d_p_l(a_l,f_nb_in))
        ENDIF
        d_p_l(:,iw)=dom_att(:)
        DEALLOCATE(dom_att)
      ENDIF
    ELSE
      CALL ipslerr (3,"flio_rbld", &
 &     "File      : "//TRIM(f_nm(iw)), &
 &     "Attribute : "//TRIM(c_wn1),"not found")
    ENDIF
!---
    c_wn1 = "DOMAIN_halo_size_start"
    CALL flioinqa (f_id_i,"?",TRIM(c_wn1),l_ex,a_l=a_l)
    IF (l_ex) THEN
      IF (a_l /= SIZE(d_d_i)) THEN
        CALL ipslerr (3,"flio_rbld", &
 &       "File      : "//TRIM(f_nm(iw)), &
 &       "size of the attribute : "//TRIM(c_wn1), &
 &       "not equal to the size of DOMAIN_dimensions_ids")
      ELSE
        ALLOCATE(dom_att(a_l))
        CALL fliogeta (f_id_i,"?",TRIM(c_wn1),dom_att)
        IF (iw == 1) THEN
          ALLOCATE (d_h_s(a_l,f_nb_in))
        ENDIF
        d_h_s(:,iw)=dom_att(:)
        DEALLOCATE(dom_att)
      ENDIF
    ELSE
      CALL ipslerr (3,"flio_rbld", &
 &     "File      : "//TRIM(f_nm(iw)), &
 &     "Attribute : "//TRIM(c_wn1),"not found")
    ENDIF
!---
    c_wn1 = "DOMAIN_halo_size_end"
    CALL flioinqa (f_id_i,"?",TRIM(c_wn1),l_ex,a_l=a_l)
    IF (l_ex) THEN
      IF (a_l /= SIZE(d_d_i)) THEN
        CALL ipslerr (3,"flio_rbld", &
 &       "File      : "//TRIM(f_nm(iw)), &
 &       "size of the attribute : "//TRIM(c_wn1), &
 &       "not equal to the size of DOMAIN_dimensions_ids")
      ELSE
        ALLOCATE(dom_att(a_l))
        CALL fliogeta (f_id_i,"?",TRIM(c_wn1),dom_att)
        IF (iw == 1) THEN
          ALLOCATE (d_h_e(a_l,f_nb_in))
        ENDIF
        d_h_e(:,iw)=dom_att(:)
        DEALLOCATE(dom_att)
      ENDIF
    ELSE
      CALL ipslerr (3,"flio_rbld", &
 &     "File      : "//TRIM(f_nm(iw)), &
 &     "Attribute : "//TRIM(c_wn1),"not found")
    ENDIF
!---
    c_wn1 = "DOMAIN_type"
    c_wn2 = " "
    CALL flioinqa (f_id_i,"?",TRIM(c_wn1),l_ex,a_l=a_l)
    IF (l_ex) THEN
      CALL fliogeta (f_id_i,"?",TRIM(c_wn1),c_wn2)
      CALL strlowercase (c_wn2)
      IF (iw == 1) THEN
        IF (    (TRIM(c_wn2) == "box") &
 &          .OR.(TRIM(c_wn2) == "apple") ) THEN
          c_d_n = c_wn2
        ELSE
          CALL ipslerr (3,"flio_rbld", &
 &         "File      : "//TRIM(f_nm(iw)), &
 &         "Attribute : "//TRIM(c_wn1), &
 &         "type "//TRIM(c_wn2)//" not (yet) supported")
        ENDIF
      ELSEIF (TRIM(c_wn2) /= TRIM(c_d_n)) THEN
        CALL ipslerr (3,"flio_rbld", &
 &       "File      : "//TRIM(f_nm(iw)), &
 &       "Attribute : "//TRIM(c_wn1), &
 &       "not equal to the one of the first file")
      ENDIF
    ELSE
      CALL ipslerr (3,"flio_rbld", &
 &     "File      : "//TRIM(f_nm(iw)), &
 &     "Attribute : "//TRIM(c_wn1),"not found")
    ENDIF
!---
    CALL flrb_cf (iw,l_ocf)
!---
  ENDDO
!-
  IF (i_v_lev >= 2) THEN
    WRITE (UNIT=*,FMT='("")')
    WRITE (*,'(" From the first file : ")')
    WRITE (*,'("   Number of dimensions : ",I2)') f_d_nb
    WRITE (*,'("     Idents  : ",(10(1X,I4),:))') f_d_i(1:f_d_nb)
    WRITE (*,'("     Lengths : ",(10(1X,I4),:))') f_d_l(1:f_d_nb)
    WRITE (*,'("     Names: ")')
    DO i=1,f_d_nb
      WRITE (*,'("       """,A,"""")') TRIM(f_d_nm(i))
    ENDDO
    IF (f_d_ul > 0) THEN
      WRITE (*,'("   Unlimited dimension id : ",I2)') f_d_i(f_d_ul)
    ENDIF
    WRITE (*,'("   Number of variables  : ",I2)') f_v_nb
    WRITE (*,'("     Names: ")')
    DO i=1,f_v_nb
      WRITE (*,'("       """,A,"""")') TRIM(f_v_nm(i))
    ENDDO
    WRITE (*,'("   Number of global attributes : ",I2)') f_a_nb
    WRITE (*,'("     Names: ")')
    DO i=1,f_a_nb
      WRITE (*,'("       """,A,"""")') TRIM(f_a_nm(i))
    ENDDO
  ENDIF
  IF (i_v_lev >= 3) THEN
    WRITE (UNIT=*,FMT='("")')
    WRITE (*,'(" From input files : ")')
    WRITE (*,'("   Total number of DOMAINS : ",I4)') d_n_t
    WRITE (*,'("   DOMAIN_dimensions_ids :",(10(1X,I5),:))') d_d_i(:)
    WRITE (*,'("   DOMAIN_size_global    :",(10(1X,I5),:))') d_s_g(:)
    WRITE (*,'("   DOMAIN_type           : """,(A),"""")') TRIM(c_d_n)
    DO iw=1,f_nb_in
      WRITE (*,'("   File   : ",A)') TRIM(f_nm(iw))
      WRITE (*,'("     d_s_l  :",(10(1X,I5),:))') d_s_l(:,iw)
      WRITE (*,'("     d_p_f  :",(10(1X,I5),:))') d_p_f(:,iw)
      WRITE (*,'("     d_p_l  :",(10(1X,I5),:))') d_p_l(:,iw)
      WRITE (*,'("     d_h_s  :",(10(1X,I5),:))') d_h_s(:,iw)
      IF (TRIM(c_d_n) == "apple") THEN
        IF (COUNT(d_h_s(:,iw) /= 0) > 1) THEN
          CALL ipslerr (3,"flio_rbld", &
 &          "Beginning offset is not yet supported", &
 &          "for more than one dimension"," ")
        ENDIF
      ENDIF
      WRITE (*,'("     d_h_e  :",(10(1X,I5),:))') d_h_e(:,iw)
      IF (TRIM(c_d_n) == "apple") THEN
        IF (COUNT(d_h_e(:,iw) /= 0) > 1) THEN
          CALL ipslerr (3,"flio_rbld", &
 &          "Ending offset is not yet supported", &
 &          "for more than one dimension"," ")
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!-
!---------------------------------------
! Create the dimensionned output file
!---------------------------------------
!-
! Define the dimensions used in the output file
  DO id=1,f_d_nb
    DO i=1,SIZE(d_d_i)
      IF (f_d_i(id) == d_d_i(i)) THEN
        f_d_l(id) = d_s_g(i)
      ENDIF
    ENDDO
  ENDDO
!-
  IF (f_d_ul > 0) THEN
    i = f_d_l(f_d_ul); f_d_l(f_d_ul) = -1;
  ENDIF
!-
! Create the output file
  CALL fliocrfd (TRIM(f_nm(f_nb)),f_d_nm,f_d_l,f_id_o,c_f_n=c_wn1)
!-
  IF (f_d_ul > 0) THEN
    f_d_l(f_d_ul) = i; itmin = 1; itmax = f_d_l(f_d_ul);
  ELSE
    itmin = 1; itmax = 1;
  ENDIF
!-
! open the first input file used to build the output file
!-
  CALL flrb_of (1,f_id_i1)
!-
! define the global attributes in the output file
! copy all global attributes except those beginning by "DOMAIN_"
! eventually actualize the "file_name" attribute
!-
  DO ia=1,f_a_nb
    IF (INDEX(TRIM(f_a_nm(ia)),"DOMAIN_") == 1)  CYCLE
    IF (TRIM(f_a_nm(ia)) == "file_name") THEN
      CALL flioputa (f_id_o,"?",TRIM(f_a_nm(ia)),TRIM(c_wn1))
    ELSE
      CALL fliocpya (f_id_i1,"?",TRIM(f_a_nm(ia)),f_id_o,"?")
    ENDIF
  ENDDO
!-
! define the variables in the output file
!-
  ALLOCATE(v_d_nb(f_v_nb)); v_d_nb(:) = 0;
  ALLOCATE(v_d_ul(f_v_nb)); v_d_ul(:) = 0;
  ALLOCATE(v_type(f_v_nb),v_d_i(flio_max_var_dims,f_v_nb));
  DO iv=1,f_v_nb
!-- get variable informations
    CALL flioinqv &
 &   (f_id_i1,TRIM(f_v_nm(iv)),l_ex,v_t=v_type(iv), &
 &    nb_dims=v_d_nb(iv),id_dims=d_i,nb_atts=v_a_nb)
!-- define the new variable
    IF (v_d_nb(iv) == 0) THEN
      CALL fliodefv &
 &     (f_id_o,TRIM(f_v_nm(iv)),v_t=v_type(iv))
    ELSE
      CALL fliodefv &
 &     (f_id_o,TRIM(f_v_nm(iv)),d_i(1:v_d_nb(iv)),v_t=v_type(iv))
      DO iw=1,v_d_nb(iv)
        IF (f_d_ul > 0) THEN
          IF (d_i(iw) == f_d_ul) THEN
            v_d_ul(iv) = iw
          ENDIF
        ENDIF
      ENDDO
      v_d_i(1:v_d_nb(iv),iv) = d_i(1:v_d_nb(iv))
    ENDIF
!-- copy all variable attributes
    IF (v_a_nb > 0) THEN
      ALLOCATE(v_a_nm(v_a_nb))
      CALL flioinqv (f_id_i1,TRIM(f_v_nm(iv)),l_ex,cn_atts=v_a_nm)
      DO ia=1,v_a_nb
        CALL fliocpya &
 &       (f_id_i1,TRIM(f_v_nm(iv)),TRIM(v_a_nm(ia)), &
 &        f_id_o,TRIM(f_v_nm(iv)))
      ENDDO
      DEALLOCATE(v_a_nm)
    ENDIF
  ENDDO
!-
! update valid_min valid_max attributes values
!-
  CALL flrb_rg
!-
!------------------------
! Fill the output file
!------------------------
!-
  DO ik=1,2
    l_uld = (ik /= 1)
    IF (l_uld) THEN
      it1=itmin; it2=itmax;
    ELSE
      it1=1; it2=1;
    ENDIF
    DO it=it1,it2
      DO iv=1,f_v_nb
        IF (    (.NOT.l_uld.AND.(v_d_ul(iv) > 0)) &
 &          .OR.(l_uld.AND.(v_d_ul(iv) <= 0)) ) THEN
          CYCLE
        ENDIF
        IF (i_v_lev >= 3) THEN
          WRITE (UNIT=*,FMT='("")')
          IF (l_uld) THEN
            WRITE (UNIT=*,FMT=*) "time step     : ",it
          ENDIF
          WRITE (UNIT=*,FMT=*) "variable      : ",TRIM(f_v_nm(iv))
          WRITE (UNIT=*,FMT=*) "var unlim dim : ",v_d_ul(iv)
        ENDIF
!------ do the variable contains dimensions to be recombined ?
        l_cgd = .FALSE.
        i_n = 1
        DO i=1,SIZE(d_d_i)
          l_cgd = ANY(v_d_i(1:v_d_nb(iv),iv) == d_d_i(i))
          l_cgd = l_cgd.AND.ANY(d_s_l(i,1:f_nb_in) /= d_s_g(i))
          IF (l_cgd) THEN
            i_n = f_nb_in
              EXIT
          ENDIF
        ENDDO
        IF (v_d_nb(iv) > 0) THEN
!-------- Allocate io_i,io_n,ia_sm,io_sm,io_cm
          i = v_d_nb(iv)
          ALLOCATE(io_i(i),io_n(i),ia_sm(i),io_sm(i),io_cm(i))
!-------- Default definition of io_i,io_n,io_sm,io_cm
          io_i(:) = 1; io_n(:) = f_d_l(v_d_i(1:v_d_nb(iv),iv));
          ia_sm(:) = 1; io_sm(:) = 1;
          IF (v_d_ul(iv) > 0) THEN
            io_i(v_d_ul(iv))=it
            io_n(v_d_ul(iv))=1
            io_sm(v_d_ul(iv))=it
          ENDIF
          io_cm(:) = io_n(:);
!-------- If needed, allocate offset
          l_o_f = .FALSE.; l_o_m = .TRUE.; l_o_l = .FALSE.;
          IF (TRIM(c_d_n) == "apple") THEN
            ALLOCATE(ia_sf(i),io_sf(i),io_cf(i))
            ALLOCATE(ia_sl(i),io_sl(i),io_cl(i))
            ia_sf(:) = 1; io_sf(:) = 1; io_cf(:) = io_n(:);
            ia_sl(:) = 1; io_sl(:) = 1; io_cl(:) = io_n(:);
            IF (v_d_ul(iv) > 0) THEN
              io_sf(v_d_ul(iv))=it
              io_sl(v_d_ul(iv))=it
            ENDIF
          ENDIF
!-------- Initialize to zero variables data 
          ! approximate dimension
          IF ( it == 1 .AND. l_cgd) THEN
          ! Enter I*J I*J is larger thant total number of single files
            if ( ((f_d_l(1)/(d_s_l(1,1)-3)) * (f_d_l(2)/(d_s_l(2,1)-3) )) .gt. d_n_t ) then 
              CALL ZeroFill (f_id_o, f_v_nm(iv), f_d_l, v_d_nb(iv), v_type(iv), v_d_i(1:v_d_nb(iv),iv))
            endif
          ENDIF
        ENDIF
!------
        DO i_i=1,i_n
          IF (l_cgd) THEN
!---------- the variable contains dimensions to be recombined
!-----------
!---------- open each file containing a small piece of data
            CALL flrb_of (i_i,f_id_i)
!-----------
!---------- do the variable has offset at first/last block ?
            l_cof = .FALSE.; l_col = .FALSE.;
            IF (TRIM(c_d_n) == "apple") THEN
              L_BF: DO id=1,v_d_nb(iv)
                DO i=1,SIZE(d_d_i)
                  IF (v_d_i(id,iv) == d_d_i(i)) THEN
                    l_cof = (d_h_s(i,i_i) /= 0)
                    IF (l_cof)  EXIT L_BF
                  ENDIF
                ENDDO
              ENDDO L_BF
              L_BL: DO id=1,v_d_nb(iv)
                DO i=1,SIZE(d_d_i)
                  IF (v_d_i(id,iv) == d_d_i(i)) THEN
                    l_col = (d_h_e(i,i_i) /= 0)
                    IF (l_col)  EXIT L_BL
                  ENDIF
                ENDDO
              ENDDO L_BL
            ENDIF
!---------- if needed, redefine start and count for dimensions
            l_o_f = .FALSE.; l_o_m = .TRUE.; l_o_l = .FALSE.;
            DO id=1,v_d_nb(iv)
              DO i=1,SIZE(d_d_i)
                IF (v_d_i(id,iv) == d_d_i(i)) THEN
                  io_n(id) = d_p_l(i,i_i)-d_p_f(i,i_i)+1
                  ia_sm(id) = 1
                  io_sm(id) = d_p_f(i,i_i)
                  io_cm(id) = io_n(id)
                  IF     (TRIM(c_d_n) == "box") THEN
                    ia_sm(id) = ia_sm(id)+d_h_s(i,i_i)
                    io_sm(id) = io_sm(id)+d_h_s(i,i_i)
                    io_cm(id) = io_cm(id)-d_h_s(i,i_i)-d_h_e(i,i_i)
                  ELSEIF (TRIM(c_d_n) == "apple") THEN
                    IF (l_cof) THEN
                      IF (d_h_s(i,i_i) /= 0) THEN
                        ia_sf(id) = 1+d_h_s(i,i_i)
                        io_sf(id) = d_p_f(i,i_i)+d_h_s(i,i_i)
                        io_cf(id) = io_n(id)-d_h_s(i,i_i)
                      ELSE
                        io_sf(id) = d_p_f(i,i_i)
                        io_cf(id) = 1
                        ia_sm(id) = ia_sm(id)+1
                        io_sm(id) = io_sm(id)+1
                        io_cm(id) = io_cm(id)-1
                        l_o_f = .TRUE.
                      ENDIF
                    ENDIF
                    IF (l_col) THEN
                      IF (d_h_e(i,i_i) /= 0) THEN
                        ia_sl(id) = 1
                        io_sl(id) = d_p_f(i,i_i)
                        io_cl(id) = io_n(id)-d_h_e(i,i_i)
                      ELSE
                        io_cm(id) = io_cm(id)-1
                        ia_sl(id) = 1+io_n(id)-1
                        io_sl(id) = d_p_f(i,i_i)+io_n(id)-1
                        io_cl(id) = 1
                        l_o_l = .TRUE.
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
            l_o_m = ALL(io_cm > 0)
          ELSE
!---------- the data can be read/write in one piece
            f_id_i = f_id_i1
          ENDIF
!---------
          IF (i_v_lev >= 3) THEN
            WRITE (UNIT=*,FMT=*) &
 &            TRIM(f_nm(i_i))//" - "//TRIM(f_v_nm(iv))
            WRITE (UNIT=*,FMT=*) "io_i  : ",io_i(:)
            WRITE (UNIT=*,FMT=*) "io_n  : ",io_n(:)
            WRITE (UNIT=*,FMT=*) "l_o_f : ",l_o_f
            IF (l_o_f) THEN
              WRITE (UNIT=*,FMT=*) "ia_sf : ",ia_sf(:)
              WRITE (UNIT=*,FMT=*) "io_sf : ",io_sf(:)
              WRITE (UNIT=*,FMT=*) "io_cf : ",io_cf(:)
            ENDIF
            WRITE (UNIT=*,FMT=*) "l_o_m : ",l_o_m
            IF (l_o_m) THEN
              WRITE (UNIT=*,FMT=*) "ia_sm : ",ia_sm(:)
              WRITE (UNIT=*,FMT=*) "io_sm : ",io_sm(:)
              WRITE (UNIT=*,FMT=*) "io_cm : ",io_cm(:)
            ENDIF
            WRITE (UNIT=*,FMT=*) "l_o_l : ",l_o_l
            IF (l_o_l) THEN
              WRITE (UNIT=*,FMT=*) "ia_sl : ",ia_sl(:)
              WRITE (UNIT=*,FMT=*) "io_sl : ",io_sl(:)
              WRITE (UNIT=*,FMT=*) "io_cl : ",io_cl(:)
            ENDIF
          ENDIF
!---------
!-------- Cases according to the type, shape and offsets of the data
!---------
          SELECT CASE (v_type(iv))
!?INTEGERS of KIND 1 are not supported on all computers
!?        CASE (flio_i1) !--- INTEGER 1
!?          SELECT CASE (v_d_nb(iv))
!?          CASE (0) !--- Scalar
!?            CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i1_0d)
!?            CALL flioputv (f_id_o,TRIM(f_v_nm(iv)),i1_0d)
!?          CASE (1) !--- 1d array
!?            ALLOCATE(i1_1d(io_n(1)))
!?            CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i1_1d, &
!? &            start=io_i(:),count=io_n(:))
!?            IF (l_o_f) THEN
!?              ib(1:1) = ia_sf(1:1); ie(1:1) = ib(1:1)+io_cf(1:1)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_1d(ib(1):ie(1)), &
!? &              start=io_sf(:),count=io_cf(:))
!?            ENDIF
!?            IF (l_o_m) THEN
!?              ib(1:1) = ia_sm(1:1); ie(1:1) = ib(1:1)+io_cm(1:1)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_1d(ib(1):ie(1)), &
!? &              start=io_sm(:),count=io_cm(:))
!?            ENDIF
!?            IF (l_o_l) THEN
!?              ib(1:1) = ia_sl(1:1); ie(1:1) = ib(1:1)+io_cl(1:1)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_1d(ib(1):ie(1)), &
!? &              start=io_sl(:),count=io_cl(:))
!?            ENDIF
!?            DEALLOCATE(i1_1d)
!?          CASE (2) !--- 2d array
!?            ALLOCATE(i1_2d(io_n(1),io_n(2)))
!?            CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i1_2d, &
!? &            start=io_i(:),count=io_n(:))
!?            IF (l_o_f) THEN
!?              ib(1:2) = ia_sf(1:2); ie(1:2) = ib(1:2)+io_cf(1:2)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_2d(ib(1):ie(1),ib(2):ie(2)), &
!? &              start=io_sf(:),count=io_cf(:))
!?            ENDIF
!?            IF (l_o_m) THEN
!?              ib(1:2) = ia_sm(1:2); ie(1:2) = ib(1:2)+io_cm(1:2)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_2d(ib(1):ie(1),ib(2):ie(2)), &
!? &              start=io_sm(:),count=io_cm(:))
!?            ENDIF
!?            IF (l_o_l) THEN
!?              ib(1:2) = ia_sl(1:2); ie(1:2) = ib(1:2)+io_cl(1:2)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_2d(ib(1):ie(1),ib(2):ie(2)), &
!? &              start=io_sl(:),count=io_cl(:))
!?            ENDIF
!?            DEALLOCATE(i1_2d)
!?          CASE (3) !--- 3d array
!?            ALLOCATE(i1_3d(io_n(1),io_n(2),io_n(3)))
!?            CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i1_3d, &
!? &            start=io_i(:),count=io_n(:))
!?            IF (l_o_f) THEN
!?              ib(1:3) = ia_sf(1:3); ie(1:3) = ib(1:3)+io_cf(1:3)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
!? &              start=io_sf(:),count=io_cf(:))
!?            ENDIF
!?            IF (l_o_m) THEN
!?              ib(1:3) = ia_sm(1:3); ie(1:3) = ib(1:3)+io_cm(1:3)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
!? &              start=io_sm(:),count=io_cm(:))
!?            ENDIF
!?            IF (l_o_l) THEN
!?              ib(1:3) = ia_sl(1:3); ie(1:3) = ib(1:3)+io_cl(1:3)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
!? &              start=io_sl(:),count=io_cl(:))
!?            ENDIF
!?            DEALLOCATE(i1_3d)
!?          CASE (4) !--- 4d array
!?            ALLOCATE(i1_4d(io_n(1),io_n(2),io_n(3),io_n(4)))
!?            CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i1_4d, &
!? &            start=io_i(:),count=io_n(:))
!?            IF (l_o_f) THEN
!?              ib(1:4) = ia_sf(1:4); ie(1:4) = ib(1:4)+io_cf(1:4)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_4d(ib(1):ie(1),ib(2):ie(2), &
!? &                    ib(3):ie(3),ib(4):ie(4)), &
!? &              start=io_sf(:),count=io_cf(:))
!?            ENDIF
!?            IF (l_o_m) THEN
!?              ib(1:4) = ia_sm(1:4); ie(1:4) = ib(1:4)+io_cm(1:4)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_4d(ib(1):ie(1),ib(2):ie(2), &
!? &                    ib(3):ie(3),ib(4):ie(4)), &
!? &              start=io_sm(:),count=io_cm(:))
!?            ENDIF
!?            IF (l_o_l) THEN
!?              ib(1:4) = ia_sl(1:4); ie(1:4) = ib(1:4)+io_cl(1:4)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_4d(ib(1):ie(1),ib(2):ie(2), &
!? &                    ib(3):ie(3),ib(4):ie(4)), &
!? &              start=io_sl(:),count=io_cl(:))
!?            ENDIF
!?            DEALLOCATE(i1_4d)
!?          CASE (5) !--- 5d array
!?            ALLOCATE(i1_5d(io_n(1),io_n(2),io_n(3),io_n(4),io_n(5)))
!?            CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i1_5d, &
!? &            start=io_i(:),count=io_n(:))
!?            IF (l_o_f) THEN
!?              ib(1:5) = ia_sf(1:5); ie(1:5) = ib(1:5)+io_cf(1:5)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
!? &                    ib(4):ie(4),ib(5):ie(5)), &
!? &              start=io_sf(:),count=io_cf(:))
!?            ENDIF
!?            IF (l_o_m) THEN
!?              ib(1:5) = ia_sm(1:5); ie(1:5) = ib(1:5)+io_cm(1:5)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
!? &                    ib(4):ie(4),ib(5):ie(5)), &
!? &              start=io_sm(:),count=io_cm(:))
!?            ENDIF
!?            IF (l_o_l) THEN
!?              ib(1:5) = ia_sl(1:5); ie(1:5) = ib(1:5)+io_cl(1:5)-1;
!?              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
!? &              i1_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
!? &                    ib(4):ie(4),ib(5):ie(5)), &
!? &              start=io_sl(:),count=io_cl(:))
!?            ENDIF
!?            DEALLOCATE(i1_5d)
!?          END SELECT
!?        CASE (flio_i2) !--- INTEGER 2
          CASE (flio_i1,flio_i2) !--- INTEGER 1/INTEGER 2
            SELECT CASE (v_d_nb(iv))
            CASE (0) !--- Scalar
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i2_0d)
              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)),i2_0d)
            CASE (1) !--- 1d array
              ALLOCATE(i2_1d(io_n(1)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i2_1d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:1) = ia_sf(1:1); ie(1:1) = ib(1:1)+io_cf(1:1)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_1d(ib(1):ie(1)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:1) = ia_sm(1:1); ie(1:1) = ib(1:1)+io_cm(1:1)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_1d(ib(1):ie(1)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:1) = ia_sl(1:1); ie(1:1) = ib(1:1)+io_cl(1:1)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_1d(ib(1):ie(1)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(i2_1d)
            CASE (2) !--- 2d array
              ALLOCATE(i2_2d(io_n(1),io_n(2)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i2_2d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:2) = ia_sf(1:2); ie(1:2) = ib(1:2)+io_cf(1:2)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_2d(ib(1):ie(1),ib(2):ie(2)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:2) = ia_sm(1:2); ie(1:2) = ib(1:2)+io_cm(1:2)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_2d(ib(1):ie(1),ib(2):ie(2)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:2) = ia_sl(1:2); ie(1:2) = ib(1:2)+io_cl(1:2)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_2d(ib(1):ie(1),ib(2):ie(2)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(i2_2d)
            CASE (3) !--- 3d array
              ALLOCATE(i2_3d(io_n(1),io_n(2),io_n(3)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i2_3d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:3) = ia_sf(1:3); ie(1:3) = ib(1:3)+io_cf(1:3)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:3) = ia_sm(1:3); ie(1:3) = ib(1:3)+io_cm(1:3)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:3) = ia_sl(1:3); ie(1:3) = ib(1:3)+io_cl(1:3)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(i2_3d)
            CASE (4) !--- 4d array
              ALLOCATE(i2_4d(io_n(1),io_n(2),io_n(3),io_n(4)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i2_4d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:4) = ia_sf(1:4); ie(1:4) = ib(1:4)+io_cf(1:4)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_4d(ib(1):ie(1),ib(2):ie(2), &
 &                      ib(3):ie(3),ib(4):ie(4)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:4) = ia_sm(1:4); ie(1:4) = ib(1:4)+io_cm(1:4)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_4d(ib(1):ie(1),ib(2):ie(2), &
 &                      ib(3):ie(3),ib(4):ie(4)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:4) = ia_sl(1:4); ie(1:4) = ib(1:4)+io_cl(1:4)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_4d(ib(1):ie(1),ib(2):ie(2), &
 &                      ib(3):ie(3),ib(4):ie(4)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(i2_4d)
            CASE (5) !--- 5d array
              ALLOCATE(i2_5d(io_n(1),io_n(2),io_n(3),io_n(4),io_n(5)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i2_5d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:5) = ia_sf(1:5); ie(1:5) = ib(1:5)+io_cf(1:5)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
 &                      ib(4):ie(4),ib(5):ie(5)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:5) = ia_sm(1:5); ie(1:5) = ib(1:5)+io_cm(1:5)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
 &                      ib(4):ie(4),ib(5):ie(5)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:5) = ia_sl(1:5); ie(1:5) = ib(1:5)+io_cl(1:5)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i2_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
 &                      ib(4):ie(4),ib(5):ie(5)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(i2_5d)
            END SELECT
          CASE (flio_i4) !--- INTEGER 4
            SELECT CASE (v_d_nb(iv))
            CASE (0) !--- Scalar
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i4_0d)
              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)),i4_0d)
            CASE (1) !--- 1d array
              ALLOCATE(i4_1d(io_n(1)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i4_1d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:1) = ia_sf(1:1); ie(1:1) = ib(1:1)+io_cf(1:1)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_1d(ib(1):ie(1)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:1) = ia_sm(1:1); ie(1:1) = ib(1:1)+io_cm(1:1)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_1d(ib(1):ie(1)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:1) = ia_sl(1:1); ie(1:1) = ib(1:1)+io_cl(1:1)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_1d(ib(1):ie(1)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(i4_1d)
            CASE (2) !--- 2d array
              ALLOCATE(i4_2d(io_n(1),io_n(2)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i4_2d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:2) = ia_sf(1:2); ie(1:2) = ib(1:2)+io_cf(1:2)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_2d(ib(1):ie(1),ib(2):ie(2)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:2) = ia_sm(1:2); ie(1:2) = ib(1:2)+io_cm(1:2)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_2d(ib(1):ie(1),ib(2):ie(2)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:2) = ia_sl(1:2); ie(1:2) = ib(1:2)+io_cl(1:2)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_2d(ib(1):ie(1),ib(2):ie(2)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(i4_2d)
            CASE (3) !--- 3d array
              ALLOCATE(i4_3d(io_n(1),io_n(2),io_n(3)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i4_3d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:3) = ia_sf(1:3); ie(1:3) = ib(1:3)+io_cf(1:3)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:3) = ia_sm(1:3); ie(1:3) = ib(1:3)+io_cm(1:3)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:3) = ia_sl(1:3); ie(1:3) = ib(1:3)+io_cl(1:3)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(i4_3d)
            CASE (4) !--- 4d array
              ALLOCATE(i4_4d(io_n(1),io_n(2),io_n(3),io_n(4)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i4_4d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:4) = ia_sf(1:4); ie(1:4) = ib(1:4)+io_cf(1:4)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_4d(ib(1):ie(1),ib(2):ie(2), &
 &                      ib(3):ie(3),ib(4):ie(4)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:4) = ia_sm(1:4); ie(1:4) = ib(1:4)+io_cm(1:4)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_4d(ib(1):ie(1),ib(2):ie(2), &
 &                      ib(3):ie(3),ib(4):ie(4)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:4) = ia_sl(1:4); ie(1:4) = ib(1:4)+io_cl(1:4)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_4d(ib(1):ie(1),ib(2):ie(2), &
 &                      ib(3):ie(3),ib(4):ie(4)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(i4_4d)
            CASE (5) !--- 5d array
              ALLOCATE(i4_5d(io_n(1),io_n(2),io_n(3),io_n(4),io_n(5)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),i4_5d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:5) = ia_sf(1:5); ie(1:5) = ib(1:5)+io_cf(1:5)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
 &                      ib(4):ie(4),ib(5):ie(5)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:5) = ia_sm(1:5); ie(1:5) = ib(1:5)+io_cm(1:5)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
 &                      ib(4):ie(4),ib(5):ie(5)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:5) = ia_sl(1:5); ie(1:5) = ib(1:5)+io_cl(1:5)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                i4_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
 &                      ib(4):ie(4),ib(5):ie(5)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(i4_5d)
            END SELECT
          CASE (flio_r4) !--- REAL 4
            SELECT CASE (v_d_nb(iv))
            CASE (0) !--- Scalar
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),r4_0d)
              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)),r4_0d)
            CASE (1) !--- 1d array
              ALLOCATE(r4_1d(io_n(1)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),r4_1d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:1) = ia_sf(1:1); ie(1:1) = ib(1:1)+io_cf(1:1)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_1d(ib(1):ie(1)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:1) = ia_sm(1:1); ie(1:1) = ib(1:1)+io_cm(1:1)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_1d(ib(1):ie(1)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:1) = ia_sl(1:1); ie(1:1) = ib(1:1)+io_cl(1:1)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_1d(ib(1):ie(1)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(r4_1d)
            CASE (2) !--- 2d array
              ALLOCATE(r4_2d(io_n(1),io_n(2)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),r4_2d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:2) = ia_sf(1:2); ie(1:2) = ib(1:2)+io_cf(1:2)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_2d(ib(1):ie(1),ib(2):ie(2)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:2) = ia_sm(1:2); ie(1:2) = ib(1:2)+io_cm(1:2)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_2d(ib(1):ie(1),ib(2):ie(2)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:2) = ia_sl(1:2); ie(1:2) = ib(1:2)+io_cl(1:2)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_2d(ib(1):ie(1),ib(2):ie(2)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(r4_2d)
            CASE (3) !--- 3d array
              ALLOCATE(r4_3d(io_n(1),io_n(2),io_n(3)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),r4_3d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:3) = ia_sf(1:3); ie(1:3) = ib(1:3)+io_cf(1:3)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:3) = ia_sm(1:3); ie(1:3) = ib(1:3)+io_cm(1:3)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:3) = ia_sl(1:3); ie(1:3) = ib(1:3)+io_cl(1:3)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(r4_3d)
            CASE (4) !--- 4d array
              ALLOCATE(r4_4d(io_n(1),io_n(2),io_n(3),io_n(4)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),r4_4d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:4) = ia_sf(1:4); ie(1:4) = ib(1:4)+io_cf(1:4)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_4d(ib(1):ie(1),ib(2):ie(2), &
 &                      ib(3):ie(3),ib(4):ie(4)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:4) = ia_sm(1:4); ie(1:4) = ib(1:4)+io_cm(1:4)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_4d(ib(1):ie(1),ib(2):ie(2), &
 &                      ib(3):ie(3),ib(4):ie(4)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:4) = ia_sl(1:4); ie(1:4) = ib(1:4)+io_cl(1:4)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_4d(ib(1):ie(1),ib(2):ie(2), &
 &                      ib(3):ie(3),ib(4):ie(4)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(r4_4d)
            CASE (5) !--- 5d array
              ALLOCATE(r4_5d(io_n(1),io_n(2),io_n(3),io_n(4),io_n(5)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),r4_5d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:5) = ia_sf(1:5); ie(1:5) = ib(1:5)+io_cf(1:5)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
 &                      ib(4):ie(4),ib(5):ie(5)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:5) = ia_sm(1:5); ie(1:5) = ib(1:5)+io_cm(1:5)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
 &                      ib(4):ie(4),ib(5):ie(5)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:5) = ia_sl(1:5); ie(1:5) = ib(1:5)+io_cl(1:5)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r4_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
 &                      ib(4):ie(4),ib(5):ie(5)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(r4_5d)
            END SELECT
          CASE (flio_r8) !--- REAL 8
            SELECT CASE (v_d_nb(iv))
            CASE (0) !--- Scalar
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),r8_0d)
              CALL flioputv (f_id_o,TRIM(f_v_nm(iv)),r8_0d)
            CASE (1) !--- 1d array
              ALLOCATE(r8_1d(io_n(1)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),r8_1d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:1) = ia_sf(1:1); ie(1:1) = ib(1:1)+io_cf(1:1)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_1d(ib(1):ie(1)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:1) = ia_sm(1:1); ie(1:1) = ib(1:1)+io_cm(1:1)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_1d(ib(1):ie(1)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:1) = ia_sl(1:1); ie(1:1) = ib(1:1)+io_cl(1:1)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_1d(ib(1):ie(1)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(r8_1d)
            CASE (2) !--- 2d array
              ALLOCATE(r8_2d(io_n(1),io_n(2)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),r8_2d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:2) = ia_sf(1:2); ie(1:2) = ib(1:2)+io_cf(1:2)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_2d(ib(1):ie(1),ib(2):ie(2)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:2) = ia_sm(1:2); ie(1:2) = ib(1:2)+io_cm(1:2)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_2d(ib(1):ie(1),ib(2):ie(2)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:2) = ia_sl(1:2); ie(1:2) = ib(1:2)+io_cl(1:2)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_2d(ib(1):ie(1),ib(2):ie(2)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(r8_2d)
            CASE (3) !--- 3d array
              ALLOCATE(r8_3d(io_n(1),io_n(2),io_n(3)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),r8_3d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:3) = ia_sf(1:3); ie(1:3) = ib(1:3)+io_cf(1:3)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:3) = ia_sm(1:3); ie(1:3) = ib(1:3)+io_cm(1:3)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:3) = ia_sl(1:3); ie(1:3) = ib(1:3)+io_cl(1:3)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_3d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(r8_3d)
            CASE (4) !--- 4d array
              ALLOCATE(r8_4d(io_n(1),io_n(2),io_n(3),io_n(4)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),r8_4d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:4) = ia_sf(1:4); ie(1:4) = ib(1:4)+io_cf(1:4)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_4d(ib(1):ie(1),ib(2):ie(2), &
 &                      ib(3):ie(3),ib(4):ie(4)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:4) = ia_sm(1:4); ie(1:4) = ib(1:4)+io_cm(1:4)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_4d(ib(1):ie(1),ib(2):ie(2), &
 &                      ib(3):ie(3),ib(4):ie(4)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:4) = ia_sl(1:4); ie(1:4) = ib(1:4)+io_cl(1:4)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_4d(ib(1):ie(1),ib(2):ie(2), &
 &                      ib(3):ie(3),ib(4):ie(4)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(r8_4d)
            CASE (5) !--- 5d array
              ALLOCATE(r8_5d(io_n(1),io_n(2),io_n(3),io_n(4),io_n(5)))
              CALL fliogetv (f_id_i,TRIM(f_v_nm(iv)),r8_5d, &
 &              start=io_i(:),count=io_n(:))
              IF (l_o_f) THEN
                ib(1:5) = ia_sf(1:5); ie(1:5) = ib(1:5)+io_cf(1:5)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
 &                      ib(4):ie(4),ib(5):ie(5)), &
 &                start=io_sf(:),count=io_cf(:))
              ENDIF
              IF (l_o_m) THEN
                ib(1:5) = ia_sm(1:5); ie(1:5) = ib(1:5)+io_cm(1:5)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
 &                      ib(4):ie(4),ib(5):ie(5)), &
 &                start=io_sm(:),count=io_cm(:))
              ENDIF
              IF (l_o_l) THEN
                ib(1:5) = ia_sl(1:5); ie(1:5) = ib(1:5)+io_cl(1:5)-1;
                CALL flioputv (f_id_o,TRIM(f_v_nm(iv)), &
 &                r8_5d(ib(1):ie(1),ib(2):ie(2),ib(3):ie(3), &
 &                      ib(4):ie(4),ib(5):ie(5)), &
 &                start=io_sl(:),count=io_cl(:))
              ENDIF
              DEALLOCATE(r8_5d)
            END SELECT
          END SELECT
!-------- eventually close each file containing a small piece of data
          CALL flrb_cf (i_i,l_ocf.AND.l_cgd.AND.(i_i /= 1))
        ENDDO
!------ If needed, deallocate io_* arrays
        IF (v_d_nb(iv) > 0) THEN
          DEALLOCATE(io_i,io_n,ia_sm,io_sm,io_cm)
          IF (TRIM(c_d_n) == "apple") THEN
            DEALLOCATE(ia_sf,io_sf,io_cf)
            DEALLOCATE(ia_sl,io_sl,io_cl)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO
!-
!-------------------
! Ending the work
!-------------------
!-
! Close files
  CALL flrb_cf (0,.TRUE.)
!-
! Deallocate
  DEALLOCATE(f_nm,f_a_id)
  DEALLOCATE(f_d_nm,f_v_nm,f_a_nm)
  DEALLOCATE(f_d_i,f_d_l)
  DEALLOCATE(v_d_nb,v_d_ul,v_type,v_d_i)
  DEALLOCATE(d_d_i,d_s_g)
  DEALLOCATE(d_s_l,d_p_f,d_p_l,d_h_s,d_h_e)
!-
  IF (i_v_lev >= 1) THEN
!-- elapsed and cpu time computation
    CALL cpu_time (t_cpu_end)
    CALL system_clock(count=nb_cc_end)
    WRITE (UNIT=*,FMT='("")')
    WRITE (UNIT=*,fmt='(" elapsed time (s) : ",1PE11.4)') &
 &   REAL(nb_cc_end-nb_cc_ini)/REAL(nb_cc_sec)
    WRITE (UNIT=*,fmt='(" CPU time (s) : ",1PE11.4)') &
 &   t_cpu_end-t_cpu_ini
  ENDIF
!=======
CONTAINS
!=======
SUBROUTINE flrb_of (i_f_n,i_f_i)
!---------------------------------------------------------------------
! Open the file of number "i_f_n" if necessary,
! and returns its identifier in "i_f_i".
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN)  :: i_f_n
  INTEGER,INTENT(OUT) :: i_f_i
!---------------------------------------------------------------------
  IF (f_a_id(i_f_n) < 0) THEN
    CALL flioopfd (TRIM(f_nm(i_f_n)),i_f_i)
    f_a_id(i_f_n) = i_f_i
  ELSE
    i_f_i = f_a_id(i_f_n)
  ENDIF
!---------------------
END SUBROUTINE flrb_of
!===
SUBROUTINE flrb_cf (i_f_n,l_cf)
!---------------------------------------------------------------------
! Close the file of number "i_f_n" if "l_cf" is TRUE.
! Close all files if "i_f_n <= 0".
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: i_f_n
  LOGICAL,INTENT(IN) :: l_cf
!---------------------------------------------------------------------
  IF (i_f_n <= 0) THEN
    CALL flioclo ()
    f_a_id(:) = -1
  ELSE
    IF (l_cf) THEN
      IF (f_a_id(i_f_n) < 0) THEN
        CALL ipslerr (2,"flio_rbld", &
 &       "The file",TRIM(f_nm(i_f_n)),"is already closed")
      ELSE
        CALL flioclo (f_a_id(i_f_n))
        f_a_id(i_f_n) = -1
      ENDIF
    ENDIF
  ENDIF
!---------------------
END SUBROUTINE flrb_cf
!===
SUBROUTINE flrb_rg
!---------------------------------------------------------------------
! Update valid_min valid_max attributes values
!---------------------------------------------------------------------
  INTEGER :: k,j
  LOGICAL :: l_vmin,l_vmax
  INTEGER(KIND=i_4) :: i4_vmin,i4_vmax
  REAL(KIND=r_4) :: r4_vmin,r4_vmax
  REAL(KIND=r_8) :: r8_vmin,r8_vmax
!---------------------------------------------------------------------
  DO k=1,f_v_nb
!-- get attribute informations
    CALL flioinqa &
 &    (f_id_i1,TRIM(f_v_nm(k)),'valid_min',l_vmin,a_t=a_type)
    CALL flioinqa &
 &    (f_id_i1,TRIM(f_v_nm(k)),'valid_max',l_vmax,a_t=a_type)
!---
    IF (l_vmin.OR.l_vmax) THEN
!---- get values of min/max
      SELECT CASE (a_type)
      CASE (flio_i1,flio_i2,flio_i4) !--- INTEGER 1/2/4
        DO j=1,f_nb_in
          CALL flrb_of (j,f_id_i)
          IF (l_vmin) THEN
            CALL fliogeta(f_id_i,TRIM(f_v_nm(k)),"valid_min",i4_0d)
            IF (j == 1) THEN
              i4_vmin = i4_0d
            ELSE
              i4_vmin = MIN(i4_vmin,i4_0d)
            ENDIF
          ENDIF
          IF (l_vmax) THEN
            CALL fliogeta(f_id_i,TRIM(f_v_nm(k)),"valid_max",i4_0d)
            IF (j == 1) THEN
              i4_vmax = i4_0d
            ELSE
              i4_vmax = MAX(i4_vmax,i4_0d)
            ENDIF
          ENDIF
          CALL flrb_cf (j,l_ocf.AND.(f_id_i /= f_id_i1))
        ENDDO
        IF (l_vmin) THEN
          CALL flioputa (f_id_o,TRIM(f_v_nm(k)),"valid_min",i4_vmin)
        ENDIF
        IF (l_vmax) THEN
          CALL flioputa (f_id_o,TRIM(f_v_nm(k)),"valid_max",i4_vmax)
        ENDIF
      CASE (flio_r4) !--- REAL 4
        DO j=1,f_nb_in
          CALL flrb_of (j,f_id_i)
          IF (l_vmin) THEN
            CALL fliogeta(f_id_i,TRIM(f_v_nm(k)),"valid_min",r4_0d)
            IF (j == 1) THEN
              r4_vmin = r4_0d
            ELSE
              r4_vmin = MIN(r4_vmin,r4_0d)
            ENDIF
          ENDIF
          IF (l_vmax) THEN
            CALL fliogeta(f_id_i,TRIM(f_v_nm(k)),"valid_max",r4_0d)
            IF (j == 1) THEN
              r4_vmax = r4_0d
            ELSE
              r4_vmax = MAX(r4_vmax,r4_0d)
            ENDIF
          ENDIF
          CALL flrb_cf (j,l_ocf.AND.(f_id_i /= f_id_i1))
        ENDDO
        IF (l_vmin) THEN
          CALL flioputa (f_id_o,TRIM(f_v_nm(k)),"valid_min",r4_vmin)
        ENDIF
        IF (l_vmax) THEN
          CALL flioputa (f_id_o,TRIM(f_v_nm(k)),"valid_max",r4_vmax)
        ENDIF
      CASE (flio_r8) !--- REAL 8
        DO j=1,f_nb_in
          CALL flrb_of (j,f_id_i)
          IF (l_vmin) THEN
            CALL fliogeta(f_id_i,TRIM(f_v_nm(k)),"valid_min",r8_0d)
            IF (j == 1) THEN
              r8_vmin = r8_0d
            ELSE
              r8_vmin = MIN(r8_vmin,r8_0d)
            ENDIF
          ENDIF
          IF (l_vmax) THEN
            CALL fliogeta(f_id_i,TRIM(f_v_nm(k)),"valid_max",r8_0d)
            IF (j == 1) THEN
              r8_vmax = r8_0d
            ELSE
              r8_vmax = MAX(r8_vmax,r8_0d)
            ENDIF
          ENDIF
          CALL flrb_cf (j,l_ocf.AND.(f_id_i /= f_id_i1))
        ENDDO
        IF (l_vmin) THEN
          CALL flioputa (f_id_o,TRIM(f_v_nm(k)),"valid_min",r8_vmin)
        ENDIF
        IF (l_vmax) THEN
          CALL flioputa (f_id_o,TRIM(f_v_nm(k)),"valid_max",r8_vmax)
        ENDIF
      END SELECT
    ENDIF
  ENDDO
!---------------------
END SUBROUTINE flrb_rg
!===
SUBROUTINE ZeroFill(f_id_o,f_v_nm,f_d_l,v_d_nb,v_type,v_d_i)

  IMPLICIT NONE
! Character length
  INTEGER,PARAMETER :: chlen=256

  INTEGER              :: v_d_nb, v_type ! variable # of dims, variable type, var Unlim dimension
  INTEGER              :: f_id_o         ! Output file ID
  INTEGER,DIMENSION(:) :: f_d_l, v_d_i    ! Global dimensions, variable dimensio ID
  CHARACTER(LEN=chlen) :: f_v_nm         ! Variable name
  INTEGER,DIMENSION(:),ALLOCATABLE :: dims

  INTEGER(KIND=i_2) :: i2_0d
  INTEGER(KIND=i_2), ALLOCATABLE :: i2_1d(:), i2_2d(:,:), i2_3d(:,:,:), i2_4d(:,:,:,:), i2_5d(:,:,:,:,:)
  INTEGER(KIND=i_4) :: i4_0d
  INTEGER(KIND=i_4), ALLOCATABLE :: i4_1d(:), i4_2d(:,:), i4_3d(:,:,:), i4_4d(:,:,:,:), i4_5d(:,:,:,:,:)
  REAL(KIND=r_4) :: r4_0d
  REAL(KIND=r_4), ALLOCATABLE    :: r4_1d(:), r4_2d(:,:), r4_3d(:,:,:), r4_4d(:,:,:,:), r4_5d(:,:,:,:,:)
  REAL(KIND=r_8) :: r8_0d
  REAL(KIND=r_8), ALLOCATABLE    :: r8_1d(:), r8_2d(:,:), r8_3d(:,:,:), r8_4d(:,:,:,:), r8_5d(:,:,:,:,:)
  
  ! write(*,*) ' Into my sub... TOM'
  ! write(*,*) f_id_o, TRIM(f_v_nm), v_d_nb , v_type
  write(*,*) 'Variable: ',TRIM(f_v_nm), ' intiliazed to zero'
  write(*,*)

  ! define variable dimension
  ALLOCATE(dims(v_d_nb)) 
  dims=f_d_l(v_d_i)
  SELECT CASE(v_type) 
    ! INTEGER 1 and 2 
    CASE (flio_i1,flio_i2)
      SELECT CASE (v_d_nb)
       CASE(1)
         ALLOCATE(i2_1d(dims(1)))
         i2_1d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),i2_1d)
         DEALLOCATE(i2_1d)
       CASE(2)
         ALLOCATE(i2_2d(dims(1),dims(2)))
         i2_2d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),i2_2d) 
         DEALLOCATE(i2_2d)
       CASE(3)
         ALLOCATE(i2_3d(dims(1),dims(2),dims(3)))
         i2_3d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),i2_3d)
         DEALLOCATE(i2_3d)
       CASE(4)
         ALLOCATE(i2_4d(dims(1),dims(2),dims(3),dims(4)))
         i2_4d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),i2_4d)
         DEALLOCATE(i2_4d)
       CASE(5)
         ALLOCATE(i2_5d(dims(1),dims(2),dims(3),dims(4),dims(5)))
         i2_5d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),i2_5d)
         DEALLOCATE(i2_5d)
      END SELECT
    ! INTEGER 4
    CASE (flio_i4)
      SELECT CASE (v_d_nb)
       CASE(1)
         ALLOCATE(i4_1d(dims(1)))
         i4_1d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),i4_1d)
         DEALLOCATE(i4_1d)
       CASE(2)
         ALLOCATE(i4_2d(dims(1),dims(2)))
         i4_2d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),i4_2d)
         DEALLOCATE(i4_2d)
       CASE(3)
         ALLOCATE(i4_3d(dims(1),dims(2),dims(3)))
         i4_3d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),i4_3d)
         DEALLOCATE(i4_3d)
       CASE(4)
         ALLOCATE(i4_4d(dims(1),dims(2),dims(3),dims(4)))
         i4_4d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),i4_4d)
         DEALLOCATE(i4_4d)
       CASE(5)
         ALLOCATE(i4_5d(dims(1),dims(2),dims(3),dims(4),dims(5)))
         i4_5d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),i4_5d)
         DEALLOCATE(i4_5d)
      END SELECT
    ! FLOAT 4
    CASE (flio_r4)
      SELECT CASE (v_d_nb)
       CASE(1)
         ALLOCATE(r4_1d(dims(1)))
         r4_1d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),r4_1d) 
         DEALLOCATE(r4_1d)
       CASE(2)
         ALLOCATE(r4_2d(dims(1),dims(2)))
         r4_2d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),r4_2d) 
         DEALLOCATE(r4_2d)
       CASE(3)
         ALLOCATE(r4_3d(dims(1),dims(2),dims(3)))
         r4_3d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),r4_3d) 
         DEALLOCATE(r4_3d)
       CASE(4)
         ALLOCATE(r4_4d(dims(1),dims(2),dims(3),dims(4)))
         r4_4d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),r4_4d)
         DEALLOCATE(r4_4d)
       CASE(5)
         ALLOCATE(r4_5d(dims(1),dims(2),dims(3),dims(4),dims(5)))
         r4_5d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),r4_5d)
         DEALLOCATE(r4_5d)
      END SELECT
    ! FLOAT 8
    CASE (flio_r8)
      SELECT CASE (v_d_nb)
       CASE(1)
         ALLOCATE(r8_1d(dims(1)))
         r8_1d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),r8_1d)
         DEALLOCATE(r8_1d)
       CASE(2)
         ALLOCATE(r8_2d(dims(1),dims(2)))
         r8_2d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),r8_2d) 
         DEALLOCATE(r8_2d)
       CASE(3)
         ALLOCATE(r8_3d(dims(1),dims(2),dims(3)))
         r8_3d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),r8_3d)
         DEALLOCATE(r8_3d)
       CASE(4)
         ALLOCATE(r8_4d(dims(1),dims(2),dims(3),dims(4)))
         r8_4d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),r8_4d)
         DEALLOCATE(r8_4d)
       CASE(5)
         ALLOCATE(r8_5d(dims(1),dims(2),dims(3),dims(4),dims(5)))
         r8_5d=0
         CALL flioputv (f_id_o,TRIM(f_v_nm),r8_5d)
         DEALLOCATE(r8_5d)
      END SELECT
  END SELECT

  DEALLOCATE (dims)

END SUBROUTINE
!===
!--------------------
END PROGRAM flio_rbld
