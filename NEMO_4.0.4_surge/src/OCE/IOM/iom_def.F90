MODULE iom_def
   !!======================================================================
   !!                    ***  MODULE  iom_def ***
   !! IOM variables definitions
   !!======================================================================
   !! History :  9.0  ! 2006 09  (S. Masson) Original code
   !!             -   ! 2007 07  (D. Storkey) Add uldname
   !!            4.0  ! 2017-11 (M. Andrejczuk) Extend IOM interface to write any 3D fields
   !!----------------------------------------------------------------------
   USE par_kind

   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER, PUBLIC ::   jpdom_data          = 1   !: ( 1  :jpiglo, 1  :jpjglo)    !!gm to be suppressed
   INTEGER, PARAMETER, PUBLIC ::   jpdom_global        = 2   !: ( 1  :jpiglo, 1  :jpjglo)
   INTEGER, PARAMETER, PUBLIC ::   jpdom_local         = 3   !: One of the 3 following cases
   INTEGER, PARAMETER, PUBLIC ::   jpdom_local_full    = 4   !: ( 1  :jpi   , 1  :jpi   )
   INTEGER, PARAMETER, PUBLIC ::   jpdom_local_noextra = 5   !: ( 1  :nlci  , 1  :nlcj  )
   INTEGER, PARAMETER, PUBLIC ::   jpdom_local_noovlap = 6   !: (nldi:nlei  ,nldj:nlej  )
   INTEGER, PARAMETER, PUBLIC ::   jpdom_unknown       = 7   !: No dimension checking
   INTEGER, PARAMETER, PUBLIC ::   jpdom_autoglo       = 8   !: 
   INTEGER, PARAMETER, PUBLIC ::   jpdom_autoglo_xy    = 9   !: Automatically set horizontal dimensions only
   INTEGER, PARAMETER, PUBLIC ::   jpdom_autodta       = 10  !: 

   INTEGER, PARAMETER, PUBLIC ::   jp_r8    = 200      !: write REAL(8)
   INTEGER, PARAMETER, PUBLIC ::   jp_r4    = 201      !: write REAL(4)
   INTEGER, PARAMETER, PUBLIC ::   jp_i4    = 202      !: write INTEGER(4)
   INTEGER, PARAMETER, PUBLIC ::   jp_i2    = 203      !: write INTEGER(2)
   INTEGER, PARAMETER, PUBLIC ::   jp_i1    = 204      !: write INTEGER(1)

   INTEGER, PARAMETER, PUBLIC ::   jpmax_files  = 100  !: maximum number of simultaneously opened file
   INTEGER, PARAMETER, PUBLIC ::   jpmax_vars   = 1200 !: maximum number of variables in one file
   INTEGER, PARAMETER, PUBLIC ::   jpmax_dims   =  4   !: maximum number of dimensions for one variable
   INTEGER, PARAMETER, PUBLIC ::   jpmax_digits =  9   !: maximum number of digits for the cpu number in the file name


!$AGRIF_DO_NOT_TREAT
   INTEGER, PUBLIC            ::   iom_open_init = 0   !: used to initialize iom_file(:)%nfid to 0
!XIOS write restart   
   LOGICAL, PUBLIC            ::   lwxios          !: write single file restart using XIOS
   INTEGER, PUBLIC            ::   nxioso          !: type of restart file when writing using XIOS 1 - single, 2 - multiple
!XIOS read restart   
   LOGICAL, PUBLIC            ::   lrxios          !: read single file restart using XIOS
   LOGICAL, PUBLIC            ::   lxios_sini = .FALSE. ! is restart in a single file
   LOGICAL, PUBLIC            ::   lxios_set  = .FALSE. 



   TYPE, PUBLIC ::   file_descriptor
      CHARACTER(LEN=240)                        ::   name     !: name of the file
      INTEGER                                   ::   nfid     !: identifier of the file (0 if closed)
                                                              !: jpioipsl option has been removed)
      INTEGER                                   ::   nvars    !: number of identified varibles in the file
      INTEGER                                   ::   iduld    !: id of the unlimited dimension
      INTEGER                                   ::   lenuld   !: length of the unlimited dimension (number of records in file)
      INTEGER                                   ::   irec     !: writing record position  
      CHARACTER(LEN=32)                         ::   uldname  !: name of the unlimited dimension
      CHARACTER(LEN=32), DIMENSION(jpmax_vars)  ::   cn_var   !: names of the variables
      INTEGER, DIMENSION(jpmax_vars)            ::   nvid     !: id of the variables
      INTEGER, DIMENSION(jpmax_vars)            ::   ndims    !: number of dimensions of the variables
      LOGICAL, DIMENSION(jpmax_vars)            ::   luld     !: variable using the unlimited dimension
      INTEGER, DIMENSION(jpmax_dims,jpmax_vars) ::   dimsz    !: size of variables dimensions 
      REAL(kind=wp), DIMENSION(jpmax_vars)      ::   scf      !: scale_factor of the variables
      REAL(kind=wp), DIMENSION(jpmax_vars)      ::   ofs      !: add_offset of the variables
      INTEGER                                   ::   nlev     ! number of vertical levels
   END TYPE file_descriptor
   TYPE(file_descriptor), DIMENSION(jpmax_files), PUBLIC ::   iom_file !: array containing the info for all opened files
   INTEGER, PARAMETER, PUBLIC                   :: max_rst_fields = 97 !: maximum number of restart variables defined in iom_set_rst_vars
   TYPE, PUBLIC :: RST_FIELD  
    CHARACTER(len=30) :: vname = "NO_NAME" ! names of variables in restart file
    CHARACTER(len=30) :: grid = "NO_GRID"
    LOGICAL           :: active =.FALSE. ! for restart write only: true - write field, false do not write field
   END TYPE RST_FIELD
!$AGRIF_END_DO_NOT_TREAT
   !
   TYPE(RST_FIELD), PUBLIC, SAVE :: rst_wfields(max_rst_fields), rst_rfields(max_rst_fields)
   !
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!======================================================================
END MODULE iom_def
