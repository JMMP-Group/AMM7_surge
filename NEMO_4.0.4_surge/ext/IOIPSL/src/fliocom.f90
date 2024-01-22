MODULE fliocom
!-
!$Id: fliocom.f90 2512 2010-12-23 15:27:09Z smasson $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!---------------------------------------------------------------------
USE netcdf
!-
USE defprec
USE calendar,  ONLY : lock_calendar,ioget_calendar, &
 &                    ioconf_calendar,ju2ymds,ymds2ju
USE errioipsl, ONLY : ipslerr,ipsldbg
USE stringop,  ONLY : strlowercase,str_xfw
!-
IMPLICIT NONE
!-
PRIVATE
!-
PUBLIC :: &
 &  fliocrfd, fliopstc, fliodefv, flioputv, flioputa, &
 &  flioopfd, flioinqf, flioinqn, fliogstc, &
 &  flioinqv, fliogetv, flioinqa, fliogeta, &
 &  fliorenv, fliorena, fliodela, fliocpya, &
 &  flioqstc, fliosync, flioclo,  fliodmpf, &
 &  flio_dom_set,    flio_dom_unset, &
 &  flio_dom_defset, flio_dom_defunset, flio_dom_definq, &
 &  flio_dom_file,   flio_dom_att
!-
!!--------------------------------------------------------------------
!! The following PUBLIC parameters (with "flio_" prefix)
!! are used in the module "fliocom" :
!!
!! flio_max_files     : maximum number of simultaneously opened files
!! flio_max_dims      : maximum number of dimensions for a file
!! flio_max_var_dims  : maximum number of dimensions for a variable
!!
!! FLIO_DOM_NONE    : "named constant" for no_domain identifier
!! FLIO_DOM_DEFAULT : "named constant" for default_domain identifier
!!
!! flio_i  : standard INTEGER external type
!! flio_r  : standard REAL external type
!! flio_c  : CHARACTER external type
!! flio_i1 : INTEGER*1 external type
!! flio_i2 : INTEGER*2 external type
!! flio_i4 : INTEGER*4 external type
!! flio_r4 : REAL*4 external type
!! flio_r8 : REAL*8 external type
!!--------------------------------------------------------------------
  INTEGER,PARAMETER,PUBLIC :: &
 &  flio_max_files=100, flio_max_dims=10, flio_max_var_dims=5
  INTEGER,PARAMETER,PUBLIC :: &
 &  flio_i = -1,        flio_r = -2,        flio_c =nf90_char, &
 &  flio_i1=nf90_int1,  flio_i2=nf90_int2,  flio_i4=nf90_int4, &
 &  flio_r4=nf90_real4, flio_r8=nf90_real8
!-
  INTEGER,PARAMETER,PUBLIC :: FLIO_DOM_NONE    =-1
  INTEGER,PARAMETER,PUBLIC :: FLIO_DOM_DEFAULT = 0
!-
!!--------------------------------------------------------------------
!! The "fliocrfd" routine creates a model file
!! which contains the dimensions needed.
!!
!! SUBROUTINE fliocrfd (f_n,f_d_n,f_d_l,f_i,id_dom,mode,c_f_n)
!!
!! INPUT
!!
!! (C) f_n      : Name of the file to be created
!! (C) f_d_n(:) : Array of (max nb_fd_mx) names of the dimensions
!! (I) f_d_l(:) : Array of (max nb_fd_mx) lengths of the dimensions
!!                For an unlimited dimension, enter a length of -1.
!!                Actually, only one unlimited dimension is supported.
!!
!! OUTPUT
!!
!! (I) f_i  : Model file identifier
!!
!! Optional INPUT arguments
!!
!! (I) id_dom : Identifier of a domain defined by calling
!!              "flio_dom_set". If this argument is present,
!!              and not equal to FLIO_DOM_NONE, it will be
!!              appended to the file name and
!!              the attributes describing the related DOMAIN
!!              will be put in the created file.
!!              This argument can be equal to FLIO_DOM_DEFAULT
!!              (see "flio_dom_defset").
!! (C) mode   : String of (case insensitive) blank-separated words
!!              defining the mode used to create the file.
!!              Supported keywords : REPLACE, 32, 64
!!              If this argument is present with the keyword "REPLACE",
!!              the file will be created in mode "CLOBBER",
!!              else the file will be created in mode "NOCLOBBER".
!!              "32/64" defines the offset mode.
!!              The default offset mode is 64 bits.
!!              Keywords "NETCDF4" and "CLASSIC" are reserved
!!              for future use.
!!
!! Optional OUTPUT arguments
!!
!! (C) c_f_n : Name of the created file.
!!             This name can be different of "f_n",
!!             if a suffix is added to the original name
!!             (".nc" or "DOMAIN_identifier.nc").
!!             The length of "c_f_n" must be sufficient
!!             to receive the created file name.
!!
!!- NOTES
!!
!! The names used to identify the spatio-temporal dimensions
!! (dimension associated to a coordinate variable)
!! are the following :
!!
!!  Axis       Names
!!
!!    x        'x[...]'  'lon[...]'
!!    y        'y[...]'  'lat[...]'
!!    z        'z[...]'  'lev[...]'  'plev[...]'   'depth[...]'
!!    t        't'       'time'      'tstep[...]'  'time_counter[...]'
!!
!! Please, apply these rules so that coordinates are
!! correctly defined.
!!--------------------------------------------------------------------
!-
!!--------------------------------------------------------------------
!! The "fliopstc" routine defines the major coordinates system
!! (spatio-temporal axis) of the model file (created by fliocrfd).
!!
!! SUBROUTINE fliopstc &
!! & (f_i,x_axis,x_axis_2d,y_axis,y_axis_2d,z_axis, &
!! &      t_axis,t_init,t_step,t_calendar)
!!
!! INPUT
!!
!! (I) f_i  : Model file identifier
!!
!! Optional INPUT arguments
!!
!! (R) x_axis(:)      : longitudinal grids
!! (R) x_axis_2d(:,:) : longitudinal grids
!! (R) y_axis(:)      : latitudinal grids
!! (R) y_axis_2d(:,:) : latitudinal grids
!! (R) z_axis(:)      : vertical grid
!! (I) t_axis(:)      : timesteps on the time axis
!! (R) t_init         : date in julian days at the beginning
!! (R) t_step         : timestep in seconds between t_axis steps
!! (C) t_calendar     : calendar
!!
!! [x/y]_axis and [x/y]_axis_2d are mutually exclusive.
!!
!!- NOTES
!!
!! The variables corresponding to the spatio-temporal coordinates
!! are created according to the following characteristics :
!!
!!- Longitude axis     x_axis / x_axis_2d
!!   Variable name     'lon'  / 'nav_lon'
!!   Attributes        Values
!!   'axis'            "X"
!!   'standard_name'   "longitude"
!!   'units'           "degrees_east"
!!   'valid_min'       MINVAL(x_axis/x_axis_2d)
!!   'valid_max'       MAXVAL(x_axis/x_axis_2d)
!!
!!- Latitude axis      y_axis / y_axis_2d
!!   Variable name     'lat'  / 'nav_lat'
!!   Attributes        Values
!!   'axis'            "Y"
!!   'standard_name'   "latitude"
!!   'units'           "degrees_north"
!!   'valid_min'       MINVAL(y_axis/y_axis_2d)
!!   'valid_max'       MAXVAL(y_axis/y_axis_2d)
!!
!!- Vertical axis      z_axis
!!   Variable name     'lev'
!!   Attributes        Values
!!   'axis'            "Z"
!!   'standard_name'   "model_level_number"
!!   'units'           "sigma_level"
!!   'long_name'       "Sigma Levels"
!!   'valid_min'       MINVAL(z_axis)
!!   'valid_max'       MAXVAL(z_axis)
!!
!!- Time axis          t_axis
!!   Variable name     'time'
!!   Attributes        Values
!!   'axis'            "T"
!!   'standard_name'   "time"
!!   'long_name'       "time steps"
!!  ['calendar'        user/default valued]
!!   'units'           calculated
!!
!! If you are not satisfied, it is possible
!! to rename variables ("fliorenv")
!! or overload the values of attributes ("flioputa").
!! Be careful : the new values you use must allow to read variables
!! as coordinates.
!!
!! The dimensions associated to the coordinates variables
!! are searched according to their names (see "fliocrfd")
!!--------------------------------------------------------------------
!-
INTERFACE fliodefv
!!--------------------------------------------------------------------
!! The "fliodefv" routines define a variable in a model file.
!!
!! SUBROUTINE fliodefv &
!! & (f_i,v_n,[v_d],v_t, &
!! &  axis,standard_name,long_name,units, &
!! &  valid_min,valid_max,fillvalue)
!!
!! INPUT
!!
!! (I)  f_i  : Model file identifier
!! (C)  v_n  : Name of variable to be defined
!! (I) [v_d] :
!!             "not present"
!!                --> scalar variable
!!             "array of one or several integers containing
!!              the identifiers of the dimensions of the variable
!!              (in the order specified to "fliocrfd"
!!               or obtained from "flioopfd")"
!!                --> multidimensioned variable
!!
!! Optional INPUT arguments
!!
!! (I) v_t : External type of the variable
!!           "present"     --> see flio_..
!!           "not present" --> type of standard real
!! (C) axis,standard_name,long_name,units : Attributes
!!     (axis should be used only for coordinates)
!! (R) valid_min,valid_max,fillvalue : Attributes
!!--------------------------------------------------------------------
  MODULE PROCEDURE &
 &  fliodv_r0d,fliodv_rnd
END INTERFACE
!-
INTERFACE flioputv
!!--------------------------------------------------------------------
!! The "flioputv" routines put a variable (defined by fliodefv)
!! in a model file.
!!
!! SUBROUTINE flioputv (f_i,v_n,v_v,start,count)
!!
!! INPUT
!!
!! (I) f_i    : model file identifier
!! (C) v_n    : name of the variable to be written
!! (R/I) v_v  : scalar or array (up to flio_max_var_dims dimensions)
!!              containing the (standard) real/integer values
!!
!! Optional INPUT arguments
!!
!! (I) start(:) : array of integers specifying the index
!!                where the first data value will be written
!! (I) count(:) : array of integers specifying the number of
!!                indices that will be written along each dimension
!!                (not present if v_v is a scalar)
!!--------------------------------------------------------------------
!?INTEGERS of KIND 1 are not supported on all computers
  MODULE PROCEDURE &
 & fliopv_i40,fliopv_i41,fliopv_i42,fliopv_i43,fliopv_i44,fliopv_i45, &
 & fliopv_i20,fliopv_i21,fliopv_i22,fliopv_i23,fliopv_i24,fliopv_i25, &
!& fliopv_i10,fliopv_i11,fliopv_i12,fliopv_i13,fliopv_i14,fliopv_i15, &
 & fliopv_r40,fliopv_r41,fliopv_r42,fliopv_r43,fliopv_r44,fliopv_r45, &
 & fliopv_r80,fliopv_r81,fliopv_r82,fliopv_r83,fliopv_r84,fliopv_r85
END INTERFACE
!-
INTERFACE flioputa
!!--------------------------------------------------------------------
!! The "flioputa" routines put a value for an attribute
!! in a model file.
!! If this attribute does not exist, it will be created.
!!
!! SUBROUTINE flioputa (f_i,v_n,a_n,a_v)
!!
!! INPUT
!!
!! (I) f_i  : Model file identifier
!! (C) v_n  : Name of variable to which the attribute is assigned.
!!            If this name is "?", the attribute will be global.
!! (C) a_n  : Name of the attribute to be defined.
!! ( ) a_v  : scalar or array of real (kind 4 or 8) or integer values,
!!            or character string
!!--------------------------------------------------------------------
  MODULE PROCEDURE &
 &  fliopa_r4_0d,fliopa_r4_1d,fliopa_r8_0d,fliopa_r8_1d, &
 &  fliopa_i4_0d,fliopa_i4_1d,fliopa_tx_0d
END INTERFACE
!-
!!--------------------------------------------------------------------
!! The "flioopfd" routine opens an existing model file,
!! and returns the dimensions used in the file and a file identifier.
!! This information can be used to allocate the space needed
!! to extract the data from the file.
!!
!! SUBROUTINE flioopfd (f_n,f_i,mode,nb_dim,nb_var,nb_gat)
!!
!! INPUT
!!
!! (C) f_n     : Name of the file to be opened
!!
!! OUTPUT
!!
!! (I) f_i      : Model file identifier
!!
!! Optional INPUT arguments
!!
!! (C) mode : Access mode to the file.
!!            If this argument is present with the value "WRITE",
!!            the file will be accessed in mode "READ-WRITE",
!!            else the file will be accessed in mode "READ-ONLY".
!!
!! Optional OUTPUT arguments
!!
!! (I) nb_dim : number of dimensions
!! (I) nb_var : number of variables
!! (I) nb_gat : number of global attributes
!!--------------------------------------------------------------------
!-
!!--------------------------------------------------------------------
!! The "flioinqf" routine returns information
!! about an opened model file given its identifier.
!!
!! SUBROUTINE flioinqf &
!! & (f_i,nb_dim,nb_var,nb_gat,id_uld,id_dim,ln_dim)
!!
!! INPUT
!!
!! (I) f_i  : Model file identifier
!!
!! Optional OUTPUT arguments
!!
!! (I) nb_dim    : number of dimensions
!! (I) nb_var    : number of variables
!! (I) nb_gat    : number of global attributes
!! (I) id_uld    : identifier of the unlimited dimension (0 if none)
!! (I) id_dim(:) : identifiers of the dimensions
!! (I) ln_dim(:) : lengths of the dimensions
!!--------------------------------------------------------------------
!-
!!--------------------------------------------------------------------
!! The "flioinqn" routine returns the names
!! of the entities encountered in an opened model file.
!!
!! SUBROUTINE flioinqn &
!! & (f_i,cn_dim,cn_var,cn_gat,cn_uld, &
!! &  id_start,id_count,iv_start,iv_count,ia_start,ia_count)
!!
!! INPUT
!!
!! (I) f_i  : Model file identifier
!!
!! Optional OUTPUT arguments
!!
!! (C) cn_dim(:) : names of dimensions
!! (C) cn_var(:) : names of variables
!! (C) cn_gat(:) : names of global attributes
!! (C) cn_uld    : names of the unlimited dimension
!!
!! Optional INPUT arguments
!!
!! (I) id_start,id_count,iv_start,iv_count,ia_start,ia_count
!!
!!  The prefix ( id       / iv      / ia              ) specifies
!!         the (dimensions/variables/global attributes) entities
!!
!!  The suffix "start" specify the index from which
!!  the first name will be retrieved (1 by default)
!!
!!  The suffix "count" specifies the number of names to be retrieved
!!  (all by default)
!!
!!  If a requested entity is not available, a "?" will be returned.
!!--------------------------------------------------------------------
!-
!!--------------------------------------------------------------------
!! The "fliogstc" routine extracts the major coordinates system
!! (spatio-temporal axis) of the model file (opened by flioopfd).
!!
!! SUBROUTINE fliogstc &
!! & (f_i,x_axis,x_axis_2d,y_axis,y_axis_2d,z_axis, &
!! &      t_axis,t_init,t_step,t_calendar, &
!! &      x_start,x_count,y_start,y_count, &
!! &      z_start,z_count,t_start,t_count)
!!
!! INPUT
!!
!! (I) f_i  : Model file identifier
!!
!! Optional OUTPUT arguments
!!
!! (R) x_axis(:)      : longitudinal grids
!! (R) x_axis_2d(:,:) : longitudinal grids
!! (R) y_axis(:)      : latitudinal grids
!! (R) y_axis_2d(:,:) : latitudinal grids
!! (R) z_axis(:)      : vertical grid
!! (I) t_axis(:)      : timesteps on the time axis
!! (R) t_init         : date in julian days at the beginning
!! (R) t_step         : timestep in seconds between t_axis steps
!! (C) t_calendar     : calendar attribute
!!                      (the value is "not found" if the attribute
!!                       is not present in the model file)
!!
!! [x/y]_axis and [x/y]_axis_2d are mutually exclusive.
!!
!! Optional INPUT arguments
!!
!! (I) x_start,x_count,y_start,y_count,z_start,z_count,t_start,t_count
!!
!!  The prefix (x/y/z/t) specifies the concerned direction.
!!
!!  The suffix "start" specify the index from which
!!  the first data value will be read (1 by default)
!!
!!  The suffix "count" specifies the number of values to be read
!!  (all by default)
!!--------------------------------------------------------------------
!-
!!--------------------------------------------------------------------
!! The "flioinqv" routine returns information about a model
!! variable given its name.
!! This information can be used to allocate the space needed
!! to extract the variable from the file.
!!
!! SUBROUTINE flioinqv &
!! & (f_i,v_n,l_ex,nb_dims,len_dims,id_dims, &
!! &  nb_atts,cn_atts,ia_start,ia_count)
!!
!! INPUT
!!
!! (I) f_i  : Model file identifier
!! (C) v_n  : Name of the variable
!!
!! OUTPUT
!!
!! (L) l_ex  : Existence of the variable
!!
!! Optional OUTPUT arguments
!!
!! (I) v_t          : External type of the variable (see flio_..)
!! (I) nb_dims      : number of dimensions of the variable
!! (I) len_dims(:)  : list of dimension lengths of the variable
!! (I) id_dims(:)   : list of dimension identifiers of the variable
!! (I) nb_atts      : number of attributes of the variable
!! (C) cn_atts(:)   : names of the attributes
!!
!! Optional INPUT arguments
!!
!! (I) ia_start : index of the first attribute whose the name
!!                will be retrieved (1 by default)
!! (I) ia_count : number of names to be retrieved (all by default)
!!
!!  If a requested entity is not available, a "?" will be returned.
!!--------------------------------------------------------------------
!-
INTERFACE fliogetv
!!--------------------------------------------------------------------
!! The "fliogetv" routines get a variable from a model file.
!!
!! SUBROUTINE fliogetv (f_i,v_n,v_v,start,count)
!!
!! INPUT
!!
!! (I) f_i  : Model file identifier
!! (C) v_n  : Name of the variable to be read
!!
!! OUTPUT
!!
!! (R/I) v_v  : scalar or array (up to flio_max_var_dims dimensions)
!!              that will contain the (standard) real/integer values
!!
!! Optional INPUT arguments
!!
!! (I) start(:) : array of integers specifying the index
!!                from which the first data value will be read
!! (I) count(:) : array of integers specifying the number of
!!                indices that will be read along each dimension
!!                (not present if v_v is a scalar)
!!--------------------------------------------------------------------
!?INTEGERS of KIND 1 are not supported on all computers
  MODULE PROCEDURE &
 & fliogv_i40,fliogv_i41,fliogv_i42,fliogv_i43,fliogv_i44,fliogv_i45, &
 & fliogv_i20,fliogv_i21,fliogv_i22,fliogv_i23,fliogv_i24,fliogv_i25, &
!& fliogv_i10,fliogv_i11,fliogv_i12,fliogv_i13,fliogv_i14,fliogv_i15, &
 & fliogv_r40,fliogv_r41,fliogv_r42,fliogv_r43,fliogv_r44,fliogv_r45, &
 & fliogv_r80,fliogv_r81,fliogv_r82,fliogv_r83,fliogv_r84,fliogv_r85
END INTERFACE
!-
!!--------------------------------------------------------------------
!! The "flioinqa" routine returns information about an
!! attribute of a variable given their names, in a model file.
!! Information about a variable includes its existence,
!! and the number of values currently stored in the attribute.
!! For a string-valued attribute, this is the number of
!! characters in the string.
!! This information can be used to allocate the space needed
!! to extract the attribute from the file.
!!
!! SUBROUTINE flioinqa (f_i,v_n,a_n,l_ex,a_t,a_l)
!!
!! INPUT
!!
!! (I) f_i : Model file identifier
!! (C) v_n : Name of variable to which the attribute is assigned.
!!           This name is "?" for a global attribute.
!! (C) a_n : Name of the concerned attribute.
!!
!! OUTPUT
!!
!! (L) l_ex : existence of the variable
!!
!! Optional OUTPUT arguments
!!
!! (I) a_t : external type of the attribute
!! (I) a_l : number of values of the attribute
!!--------------------------------------------------------------------
!-
INTERFACE fliogeta
!!--------------------------------------------------------------------
!! The "fliogeta" routines get a value for an attribute
!! in a model file.
!!
!! SUBROUTINE fliogeta (f_i,v_n,a_n,a_v)
!!
!! INPUT
!!
!! (I) f_i  : Model file identifier
!! (C) v_n  : Name of variable to which the attribute is assigned.
!!            This name is "?" for a global attribute.
!! (C) a_n  : Name of the attribute to be retrieved.
!! ( ) a_v  : scalar or array of real (kind 4 or 8) or integer values,
!!            or character string
!!--------------------------------------------------------------------
  MODULE PROCEDURE &
 &  flioga_r4_0d,flioga_r4_1d,flioga_r8_0d,flioga_r8_1d, &
 &  flioga_i4_0d,flioga_i4_1d,flioga_tx_0d
END INTERFACE
!-
!!--------------------------------------------------------------------
!! The "fliorenv" routine renames a variable, in a model file.
!!
!! SUBROUTINE fliorenv (f_i,v_o_n,v_n_n)
!!
!! INPUT
!!
!! (I) f_i    : Model file identifier
!! (C) v_o_n  : Old name of the variable
!! (C) v_n_n  : New name of the variable
!!--------------------------------------------------------------------
!-
!!--------------------------------------------------------------------
!! The "fliorena" routine renames an attribute
!! of a variable, in a model file.
!!
!! SUBROUTINE fliorena (f_i,v_n,a_o_n,a_n_n)
!!
!! INPUT
!!
!! (I) f_i    : Model file identifier
!! (C) v_n    : Name of variable to which the attribute is assigned.
!!              This name is "?" for a global attribute.
!! (C) a_o_n  : Old name of the concerned attribute.
!! (C) a_n_n  : New name of the concerned attribute.
!!--------------------------------------------------------------------
!-
!!--------------------------------------------------------------------
!! The "fliodela" routine deletes an attribute in a model file.
!!
!! SUBROUTINE fliodela (f_i,v_n,a_n)
!!
!! INPUT
!!
!! (I) f_i  : Model file identifier
!! (C) v_n  : Name of variable to which the attribute is assigned.
!!            This name is "?" for a global attribute.
!! (C) a_n  : Name of the concerned attribute.
!!--------------------------------------------------------------------
!-
!!--------------------------------------------------------------------
!! The "fliocpya" routine copies an attribute
!! from one open model file to another.
!! It can also be used to copy an attribute from
!! one variable to another within the same model file.
!!
!! SUBROUTINE fliocpya (f_i_i,v_n_i,a_n,f_i_o,v_n_o)
!!
!! INPUT
!!
!! (I) f_i_i : Identifier of the input  model file
!! (C) v_n_i : Name of the input variable
!!             This name is "?" for a global attribute.
!! (C) a_n   : Name of the concerned attribute.
!! (I) f_i_o : Identifier of the output model file
!!             It can be the same as the input identifier.
!! (C) v_n_o : Name of the output variable
!!             This name is "?" for a global attribute.
!!--------------------------------------------------------------------
!-
!!--------------------------------------------------------------------
!! The "flioqstc" routine search for a spatio-temporal coordinate
!! in a model file and returns its name.
!!
!! SUBROUTINE flioqstc (f_i,c_type,l_ex,c_name)
!!
!! INPUT
!!
!! (I) f_i     : Model file identifier
!! (C) c_type  : Type of the coordinate ("x"/"y"/"z"/"t")
!!
!! OUTPUT
!!
!! (L) l_ex    : existence of the coordinate
!! (C) c_name  : name of the coordinate
!!
!!- NOTES
!!
!! The following rules are used for searching variables
!! which are spatio-temporal coordinates (x/y/z/t).
!!
!!-- Rule 1 : we look for a variable with one dimension
!!--          and which has the same name as its dimension
!!
!!-- Rule 2 : we look for a correct "axis" attribute
!!
!!  Axis       Axis attribute             Number of dimensions
!!             (case insensitive)
!!
!!    x         X                         1/2
!!    y         Y                         1/2
!!    z         Z                         1
!!    t         T                         1
!!
!!-- Rule 3 : we look for a correct "standard_name" attribute
!!
!!  Axis       Axis attribute          Number of dimensions
!!             (case insensitive)
!!
!!    x         longitude              1/2
!!    y         latitude               1/2
!!    z         model_level_number     1
!!    t         time                   1
!!
!!-- Rule 4 : we look for a specific name
!!
!!  Axis   Names
!!
!!    x    'nav_lon' 'lon'    'longitude'
!!    y    'nav_lat' 'lat'    'latitude'
!!    z    'depth'   'deptht' 'height'      'level'
!!         'lev'     'plev'   'sigma_level' 'layer'
!!    t    'time'    'tstep'  'timesteps'
!!
!!--------------------------------------------------------------------
!-
!!--------------------------------------------------------------------
!! The "fliosync" routine synchronise one or all opened model files,
!! to minimize data loss in case of abnormal termination.
!!
!! SUBROUTINE fliosync (f_i)
!!
!! Optional INPUT arguments
!!
!! (I) f_i  : Model file identifier
!!            If this argument is not present,
!!            all the opened model files are synchronised.
!---------------------------------------------------------------------
!-
!!--------------------------------------------------------------------
!! The "flioclo" routine closes one or all opened model files
!! and frees the space needed to keep information about the files
!!
!! SUBROUTINE flioclo (f_i)
!!
!! Optional INPUT arguments
!!
!! (I) f_i  : Model file identifier
!!            If this argument is not present,
!!            all the opened model files are closed.
!!--------------------------------------------------------------------
!-
!!--------------------------------------------------------------------
!! The "fliodmpf" routine dumps a model file
!! and prints the result on the standard output.
!!
!! SUBROUTINE fliodmpf (f_n)
!!
!! INPUT
!!
!! (C) f_n  : Name of the model file to be dumped
!!--------------------------------------------------------------------
!-
!!--------------------------------------------------------------------
!! This "flio_dom_set" sets up the domain activity of IOIPSL.
!! It stores all the domain information and allows it to be stored
!! in the model file and change the file names.
!!
!! This routine must be called by the user before opening
!! the model file.
!!
!! SUBROUTINE flio_dom_set &
!!  & (dtnb,dnb,did,dsg,dsl,dpf,dpl,dhs,dhe,cdnm,id_dom)
!!
!! INPUT
!!
!! (I) dtnb   : total number of domains
!! (I) dnb    : domain number
!! (I) did(:) : distributed dimensions identifiers
!!              (up to 5 dimensions are supported)
!! (I) dsg(:) : total number of points for each dimension
!! (I) dsl(:) : local number of points for each dimension
!! (I) dpf(:) : position of first local point for each dimension
!! (I) dpl(:) : position of last local point for each dimension
!! (I) dhs(:) : start halo size for each dimension
!! (I) dhe(:) : end halo size for each dimension
!! (C) cdnm   : Model domain definition name.
!!              The names actually supported are :
!!              "BOX", "APPLE", "ORANGE".
!!              These names are case insensitive.
!!
!! OUTPUT argument
!!
!! (I) id_dom : Model domain identifier
!!
!!--------------------------------------------------------------------
!!
!!--------------------------------------------------------------------
!! The "flio_dom_unset" routine unsets one or all set domains
!! and frees the space needed to keep information about the domains
!!
!! This routine should be called by the user to free useless domains.
!!
!! SUBROUTINE flio_dom_unset (id_dom)
!!
!! Optional INPUT arguments
!!
!! (I) id_dom : Model domain identifier
!!      >=1 & <= dom_max_nb : the domain is closed
!!      not present         : all the set model domains are unset
!!--------------------------------------------------------------------
!!
!!--------------------------------------------------------------------
!! The "flio_dom_defset" sets
!! the default domain identifier.
!!
!! SUBROUTINE flio_dom_defset (id_dom)
!!
!! INPUT argument
!!
!! (I) id_dom : Model default domain identifier
!!     ( >=1 & <= dom_max_nb )
!!     This identifier will be able to be taken by calling
!!     "flio_dom_definq" and used to create model files
!!     with the corresponding domain definitions
!!--------------------------------------------------------------------
!!
!!--------------------------------------------------------------------
!! The "flio_dom_defunset" routine unsets
!! the default domain identifier.
!!
!! SUBROUTINE flio_dom_defunset ()
!!
!!--------------------------------------------------------------------
!!
!!--------------------------------------------------------------------
!! The "flio_dom_definq" routine inquires about
!! the default domain identifier.
!! You should call this procedure to safeguard the current
!! default domain identifier if you wish to use locally
!! another default domain, in order to restore it.
!!
!! SUBROUTINE flio_dom_definq (id_dom)
!!
!! OUTPUT argument
!!
!! (I) id_dom : Model default domain identifier
!!     IF no default domain identifier has been set,
!!     the returned value is "FLIO_DOM_NONE".
!!--------------------------------------------------------------------
!-
!---------------------------------------------------------------------
! This is the data we keep concerning each file we open
!---------------------------------------------------------------------
!- For each file
!- (I) nw_id(f_i)   : index to access at this file
!- (I) nw_nd(f_i)   : number of dimensions
!- (I) nw_nv(f_i)   : number of variables
!- (I) nw_na(f_i)   : number of global attributes
!- (I) nw_un(f_i)   : ID of the first unlimited dimension
!- (L) lw_hm(f_i)   : for mode handling (.TRUE. define, .FALSE. data)
!- (I) nw_di(:,f_i) : dimension IDs in the file "f_i"
!- (I) nw_dl(:,f_i) : dimension lengths in the file "f_i"
!- (I) nw_ai(:,f_i) : dimension Ids for the axis in the file "f_i"
!---------------------------------------------------------------------
  INTEGER,PARAMETER :: &
 &  nb_fi_mx=flio_max_files, &
 &  nb_fd_mx=flio_max_dims, &
 &  nb_vd_mx=flio_max_var_dims
  INTEGER,PARAMETER :: nb_ax_mx=4
!-
  INTEGER,PARAMETER :: k_lon=1, k_lat=2, k_lev=3, k_tim=4
!-
  INTEGER,DIMENSION(nb_fi_mx),SAVE :: &
 &  nw_id=-1,nw_nd,nw_nv,nw_na,nw_un
  LOGICAL,DIMENSION(nb_fi_mx),SAVE :: lw_hm
  INTEGER,DIMENSION(nb_fd_mx,nb_fi_mx),SAVE :: nw_di=-1,nw_dl=-1
  INTEGER,DIMENSION(nb_ax_mx,nb_fi_mx),SAVE :: nw_ai=-1
!-
! Maximum number of simultaneously defined domains
  INTEGER,PARAMETER :: dom_max_nb=200
!-
! Maximum number of distributed dimensions for each domain
  INTEGER,PARAMETER :: dom_max_dims=5
!-
! Default domain identifier
  INTEGER,SAVE :: id_def_dom=FLIO_DOM_NONE
!-
! Supported domain definition names
  INTEGER,PARAMETER :: n_dns=3, l_dns=7
  CHARACTER(LEN=l_dns),DIMENSION(n_dns),SAVE :: &
 &  c_dns=(/ "box    ","apple  ","orange "/)
!-
! DOMAINS related variables
  INTEGER,DIMENSION(1:dom_max_nb),SAVE :: &
 &  d_d_n=-1, d_n_t=0, d_n_c=0
  INTEGER,DIMENSION(1:dom_max_dims,1:dom_max_nb),SAVE :: &
 &  d_d_i, d_s_g, d_s_l, d_p_f, d_p_l, d_h_s, d_h_e
  CHARACTER(LEN=l_dns),DIMENSION(1:dom_max_nb),SAVE :: c_d_t
!-
!===
CONTAINS
!===
!-
!---------------------------------------------------------------------
!- Public procedures
!---------------------------------------------------------------------
!-
!===
SUBROUTINE fliocrfd (f_n,f_d_n,f_d_l,f_i,id_dom,mode,c_f_n)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(IN) :: f_n
  CHARACTER(LEN=*),DIMENSION(:),INTENT(IN) :: f_d_n
  INTEGER,DIMENSION(:),INTENT(IN) :: f_d_l
  INTEGER,INTENT(OUT) :: f_i
  INTEGER,OPTIONAL,INTENT(IN) :: id_dom
  CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: mode
  CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: c_f_n
!-
  INTEGER :: i_rc,f_e,idid,ii,m_c,n_u
  CHARACTER(LEN=NF90_MAX_NAME) :: f_nw
  INTEGER,PARAMETER :: l_string=80,l_word=10
  CHARACTER(LEN=l_string) :: c_string
  CHARACTER(LEN=l_word)   :: c_word
  LOGICAL :: l_ok
  INTEGER,PARAMETER :: k_replace=1
  INTEGER,PARAMETER :: k_32=1,k_64=2
!- !? : Code to be activated for NETCDF4
!?  INTEGER,PARAMETER :: k_netcdf4=1,k_classic=1
  INTEGER,PARAMETER :: n_opt=4
  INTEGER,DIMENSION(n_opt) :: i_opt
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->fliocrfd - file name : ",TRIM(f_n)
  ENDIF
!-
! Search for a free local identifier
  f_i = flio_rid()
  IF (f_i < 0) THEN
    CALL ipslerr (3,'fliocrfd', &
 &   'Too many files.','Please increase nb_fi_mx', &
 &   'in module fliocom.f90.')
  ENDIF
!-
! Update the name of the file
  f_nw = f_n
  CALL flio_dom_file (f_nw,id_dom)
!-
! Check the dimensions
  IF (SIZE(f_d_l) /= SIZE(f_d_n)) THEN
    CALL ipslerr (3,'fliocrfd', &
 &   'The number of names is not equal to the number of lengths', &
 &   'for the dimensions of the file',TRIM(f_nw))
  ENDIF
  IF (SIZE(f_d_l) > nb_fd_mx) THEN
    CALL ipslerr (3,'fliocrfd', &
 &   'Too many dimensions','to create the file',TRIM(f_nw))
  ENDIF
!-
! Check the mode
!-
  i_opt(:)=-1
!-
  IF (PRESENT(mode)) THEN
!---
    IF (LEN_TRIM(mode) > l_string) THEN
      CALL ipslerr (3,'fliocrfd', &
 &     '"mode" argument','too long','to be treated')
    ENDIF
    c_string = mode(:)
    CALL strlowercase (c_string)
!---
    DO
      CALL str_xfw  (c_string,c_word,l_ok)
      IF (l_ok) THEN
!- !? : Code to be activated for NETCDF4
        SELECT CASE (TRIM(c_word))
        CASE('replace')
          IF (i_opt(1) > 0) THEN
            CALL ipslerr (3,'fliocrfd', &
 &           'Replace option','already','defined')
          ELSE
            i_opt(1) = k_replace
          ENDIF
!?      CASE('netcdf4')
!?        IF (i_opt(2) > 0) THEN
!?          CALL ipslerr (3,'fliocrfd', &
!? &         'Netcdf4 format','already','defined')
!?        ELSE
!?          i_opt(2) = k_netcdf4
!?        ENDIF
        CASE('32')
          IF (i_opt(3) > 0) THEN
            CALL ipslerr (3,'fliocrfd', &
 &           'Offset format','already','defined')
          ELSE
            i_opt(3) = k_32
          ENDIF
        CASE('64')
          IF (i_opt(3) > 0) THEN
            CALL ipslerr (3,'fliocrfd', &
 &           'Offset format','already','defined')
          ELSE
            i_opt(3) = k_64
          ENDIF
!?      CASE('CLASSIC')
!?        IF (i_opt(4) > 0) THEN
!?          CALL ipslerr (3,'fliocrfd', &
!? &         'Netcdf4 classic format','already','defined')
!?        ELSE
!?          i_opt(4) = k_classic
!?        ENDIF
        CASE DEFAULT
          CALL ipslerr (3,'fliocrfd', &
 &         'Option '//TRIM(c_word),'not','supported')
        END SELECT
      ELSE
        EXIT
      ENDIF
    ENDDO
  ENDIF
!-
  IF (i_opt(1) == k_replace) THEN
    m_c = NF90_CLOBBER
  ELSE
    m_c = NF90_NOCLOBBER
  ENDIF
!-
!- Code to be replaced by the following for NETCDF4
!?  IF (i_opt(2) == k_netcdf4) THEN
!?    m_c = IOR(m_c,NF90_NETCDF4)
!?    IF (i_opt(3) > 0) THEN
!?      CALL ipslerr (3,'fliocrfd', &
!? &     'Netcdf4 format','and offset option','are not compatible')
!?    ELSE IF (i_opt(4) == k_classic) THEN
!?      m_c = IOR(m_c,NF90_CLASSIC_MODEL)
!?    ENDIF
!?   LSE IF (i_opt(4) > 0) THEN
!?    CALL ipslerr (3,'fliocrfd', &
!? &   'Classic option','is reserved','for the Netcdf4 format')
!?  ELSE
    IF (i_opt(3) /= k_32) THEN
      m_c = IOR(m_c,NF90_64BIT_OFFSET)
    ENDIF
!?  ENDIF
!-
! Create file (and enter the definition mode)
  i_rc = NF90_CREATE(f_nw,m_c,f_e)
  lw_hm(f_i) = .TRUE.
  IF (i_rc /= NF90_NOERR) THEN
    CALL ipslerr (3,'fliocrfd', &
 &   'Could not create file :',TRIM(f_nw), &
 &   TRIM(NF90_STRERROR(i_rc))//' (Netcdf)')
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) '  fliocrfd, external model file-id : ',f_e
  ENDIF
!-
! Create dimensions
  n_u = 0
  DO ii=1,SIZE(f_d_l)
    IF (f_d_l(ii) == -1) THEN
      IF (n_u == 0) THEN
        i_rc = NF90_DEF_DIM(f_e,TRIM(f_d_n(ii)),NF90_UNLIMITED,idid)
        n_u = n_u+1
      ELSE
        CALL ipslerr (3,'fliocrfd', &
 &       'Can not handle more than one unlimited dimension', &
 &       'for file :',TRIM(f_nw))
      ENDIF
    ELSE IF (f_d_l(ii) > 0) THEN
      i_rc = NF90_DEF_DIM(f_e,TRIM(f_d_n(ii)),f_d_l(ii),idid)
    ENDIF
    IF ( ((f_d_l(ii) == -1).OR.(f_d_l(ii) > 0)) &
 &      .AND.(i_rc /= NF90_NOERR) ) THEN
      CALL ipslerr (3,'fliocrfd', &
 &     'One dimension can not be defined', &
 &     'for the file :',TRIM(f_nw))
    ENDIF
  ENDDO
!-
! Define "Conventions" global attribute
  i_rc = NF90_PUT_ATT(f_e,NF90_GLOBAL,'Conventions',"CF-1.1")
!-
! Add the DOMAIN attributes if needed
  CALL flio_dom_att (f_e,id_dom)
!-
! Keep the file information
  nw_id(f_i) = f_e
  CALL flio_inf (f_e, &
 &  nb_dims=nw_nd(f_i),id_unlm=nw_un(f_i),nb_atts=nw_na(f_i), &
 &  nn_idm=nw_di(:,f_i),nn_ldm=nw_dl(:,f_i),nn_aid=nw_ai(:,f_i))
!-
! Return the created file name if needed
  IF (PRESENT(c_f_n)) THEN
    IF (LEN(c_f_n) >= LEN_TRIM(f_nw)) THEN
      c_f_n = TRIM(f_nw)
    ELSE
      CALL ipslerr (3,'fliocrfd', &
 &     'the length of "c_f_n" is not sufficient to receive', &
 &     'the name of the created file :',TRIM(f_nw))
    ENDIF
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) '<-fliocrfd'
  ENDIF
!----------------------
END SUBROUTINE fliocrfd
!===
SUBROUTINE fliopstc &
 & (f_i,x_axis,x_axis_2d,y_axis,y_axis_2d,z_axis, &
 &      t_axis,t_init,t_step,t_calendar)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  REAL,DIMENSION(:),OPTIONAL,INTENT(IN)    :: x_axis,y_axis
  REAL,DIMENSION(:,:),OPTIONAL,INTENT(IN)  :: x_axis_2d,y_axis_2d
  REAL,DIMENSION(:),OPTIONAL,INTENT(IN)    :: z_axis
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: t_axis
  CHARACTER(LEN=*),OPTIONAL,INTENT(IN)     :: t_calendar
  REAL,OPTIONAL,INTENT(IN)                 :: t_init,t_step
!-
  INTEGER :: i_rc,f_e
  INTEGER :: lonid,latid,levid,timeid
  INTEGER :: j_yy,j_mo,j_dd,j_hh,j_mn,j_ss
  REAL    :: dt,r_ss,v_min,v_max
  INTEGER :: k,k_1,k_2
  LOGICAL :: l_tmp
  CHARACTER(LEN=20) :: c_tmp1
  CHARACTER(LEN=40) :: c_tmp2
  CHARACTER(LEN=80) :: c_tmp3
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->fliopstc"
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('fliopstc',f_i,f_e)
!-
! Validate the coherence of the arguments
!-
  IF (    (PRESENT(x_axis).AND.PRESENT(x_axis_2d)) &
 &    .OR.(PRESENT(y_axis).AND.PRESENT(y_axis_2d)) ) THEN
    CALL ipslerr (3,'fliopstc', &
 &    'The [x/y]_axis arguments', &
 &    'are not coherent :',&
 &    'can not handle two [x/y]_axis')
  ENDIF
!-
  IF (    PRESENT(x_axis).OR.PRESENT(x_axis_2d) &
 &    .OR.PRESENT(y_axis).OR.PRESENT(y_axis_2d) ) THEN
    k_1=nw_ai(k_lon,f_i); k_2=nw_ai(k_lat,f_i);
  ENDIF
!-
! Define the longitude axis
!-
  IF (PRESENT(x_axis).OR.PRESENT(x_axis_2d)) THEN
!---
    IF (l_dbg) THEN
      WRITE(*,*) '  fliopstc : Define the Longitude axis'
    ENDIF
!---
    IF (PRESENT(x_axis)) THEN
      IF (SIZE(x_axis) /= nw_dl(k_1,f_i)) THEN
        CALL ipslerr (3,'fliopstc', &
 &       'Invalid x_axis dimension :', &
 &       'not equal to the dimension', &
 &       'defined at the creation of the file')
      ENDIF
    ELSE
      IF (    (SIZE(x_axis_2d,DIM=1) /= nw_dl(k_1,f_i)) &
 &        .OR.(SIZE(x_axis_2d,DIM=2) /= nw_dl(k_2,f_i)) ) THEN
        CALL ipslerr (3,'fliopstc', &
 &       'Invalid x_axis_2d dimensions :', &
 &       'not equal to the dimensions', &
 &       'defined at the creation of the file')
      ENDIF
    ENDIF
!---
    CALL flio_hdm (f_i,f_e,.TRUE.)
    IF (PRESENT(x_axis)) THEN
      i_rc = NF90_DEF_VAR(f_e,"lon",NF90_REAL4, &
 &                        nw_di(k_1,f_i),lonid)
      v_min = MINVAL(x_axis)
      v_max = MAXVAL(x_axis)
    ELSE
      i_rc = NF90_DEF_VAR(f_e,"nav_lon",NF90_REAL4, &
 &             nw_di((/k_1,k_2/),f_i),lonid)
      v_min = MINVAL(x_axis_2d)
      v_max = MAXVAL(x_axis_2d)
    ENDIF
    i_rc = NF90_PUT_ATT(f_e,lonid,"axis","X")
    i_rc = NF90_PUT_ATT(f_e,lonid,'standard_name',"longitude")
    i_rc = NF90_PUT_ATT(f_e,lonid,'units',"degrees_east")
    i_rc = NF90_PUT_ATT(f_e,lonid,'valid_min',REAL(v_min,KIND=4))
    i_rc = NF90_PUT_ATT(f_e,lonid,'valid_max',REAL(v_max,KIND=4))
  ENDIF
!-
! Define the Latitude axis
!-
  IF (PRESENT(y_axis).OR.PRESENT(y_axis_2d)) THEN
!---
    IF (l_dbg) THEN
      WRITE(*,*) '  fliopstc : Define the Latitude axis'
    ENDIF
!---
    IF (PRESENT(y_axis)) THEN
      IF (SIZE(y_axis) /= nw_dl(k_2,f_i)) THEN
        CALL ipslerr (3,'fliopstc', &
 &       'Invalid y_axis dimension :', &
 &       'not equal to the dimension', &
 &       'defined at the creation of the file')
      ENDIF
    ELSE
      IF (    (SIZE(y_axis_2d,DIM=1) /= nw_dl(k_1,f_i)) &
 &        .OR.(SIZE(y_axis_2d,DIM=2) /= nw_dl(k_2,f_i)) ) THEN
        CALL ipslerr (3,'fliopstc', &
 &       'Invalid y_axis_2d dimensions :', &
 &       'not equal to the dimensions', &
 &       'defined at the creation of the file')
      ENDIF
    ENDIF
!---
    CALL flio_hdm (f_i,f_e,.TRUE.)
    IF (PRESENT(y_axis)) THEN
      i_rc = NF90_DEF_VAR(f_e,"lat",NF90_REAL4, &
 &                        nw_di(k_2,f_i),latid)
      v_min = MINVAL(y_axis)
      v_max = MAXVAL(y_axis)
    ELSE
      i_rc = NF90_DEF_VAR(f_e,"nav_lat",NF90_REAL4, &
 &             nw_di((/k_1,k_2/),f_i),latid)
      v_min = MINVAL(y_axis_2d)
      v_max = MAXVAL(y_axis_2d)
    ENDIF
    i_rc = NF90_PUT_ATT(f_e,latid,"axis","Y")
    i_rc = NF90_PUT_ATT(f_e,latid,'standard_name',"latitude")
    i_rc = NF90_PUT_ATT(f_e,latid,'units',"degrees_north")
    i_rc = NF90_PUT_ATT(f_e,latid,'valid_min',REAL(v_min,KIND=4))
    i_rc = NF90_PUT_ATT(f_e,latid,'valid_max',REAL(v_max,KIND=4))
  ENDIF
!-
! Define the Vertical axis
!-
  IF (PRESENT(z_axis)) THEN
!---
    IF (l_dbg) THEN
      WRITE(*,*) '  fliopstc : Define the Vertical axis'
    ENDIF
!---
    k_1=nw_ai(k_lev,f_i);
!---
    IF (SIZE(z_axis) /= nw_dl(k_1,f_i)) THEN
      CALL ipslerr (3,'fliopstc', &
 &     'Invalid z_axis dimension :', &
 &     'not equal to the dimension', &
 &     'defined at the creation of the file')
    ENDIF
!---
    v_min = MINVAL(z_axis)
    v_max = MAXVAL(z_axis)
!---
    CALL flio_hdm (f_i,f_e,.TRUE.)
    i_rc = NF90_DEF_VAR(f_e,'lev',NF90_REAL4, &
 &                      nw_di(k_1,f_i),levid)
    i_rc = NF90_PUT_ATT(f_e,levid,"axis","Z")
    i_rc = NF90_PUT_ATT(f_e,levid,'standard_name','model_level_number')
    i_rc = NF90_PUT_ATT(f_e,levid,'units','sigma_level')
    i_rc = NF90_PUT_ATT(f_e,levid,'long_name','Sigma Levels')
    i_rc = NF90_PUT_ATT(f_e,levid,'valid_min',REAL(v_min,KIND=4))
    i_rc = NF90_PUT_ATT(f_e,levid,'valid_max',REAL(v_max,KIND=4))
  ENDIF
!-
! Define the Time axis
!-
  IF (PRESENT(t_axis).AND.PRESENT(t_init).AND.PRESENT(t_step)) THEN
!---
    IF (l_dbg) THEN
      WRITE(*,*) '  fliopstc : Define the Time axis'
    ENDIF
!---
    k_1=nw_ai(k_tim,f_i);
!---
    IF (     (nw_dl(k_1,f_i) /= 0) &
 &      .AND.(SIZE(t_axis) /= nw_dl(k_1,f_i)) ) THEN
      CALL ipslerr (3,'fliopstc', &
 &     'Invalid t_axis dimension :', &
 &     'not equal to the dimension', &
 &     'defined at the creation of the file')
    ENDIF
!-- Retrieve the calendar date
    CALL lock_calendar (old_status=l_tmp)
    IF (PRESENT(t_calendar)) THEN
      CALL ioget_calendar (c_tmp1)
      CALL lock_calendar (new_status=.FALSE.)
      CALL ioconf_calendar (TRIM(t_calendar))
    ENDIF
    CALL ju2ymds (t_init,j_yy,j_mo,j_dd,r_ss)
    IF (PRESENT(t_calendar)) THEN
      CALL lock_calendar (new_status=.FALSE.)
      CALL ioconf_calendar (TRIM(c_tmp1))
    ENDIF
    CALL lock_calendar (new_status=l_tmp)
!--
    k=NINT(r_ss)
    j_hh=k/3600
    k=k-3600*j_hh
    j_mn=k/60
    j_ss=k-60*j_mn
!-- Calculate the step unit
    IF      (ABS(t_step) >= 604800.) THEN
      dt = t_step/604800.
      c_tmp2 = 'weeks'
    ELSE IF (ABS(t_step) >= 86400.) THEN
      dt = t_step/86400.
      c_tmp2 = 'days'
    ELSE IF (ABS(t_step) >=  3600.) THEN
      dt = t_step/3600.
      c_tmp2 = 'hours'
    ELSE IF (ABS(t_step) >=    60.) THEN
      dt = t_step/60.
      c_tmp2 = 'minutes'
    ELSE
      dt = t_step
      c_tmp2 = 'seconds'
    ENDIF
!---
    c_tmp1 = ''
    IF (ABS(dt-NINT(dt)) <= ABS(10.*EPSILON(dt))) THEN
      IF (NINT(dt) /= 1) THEN
        WRITE (UNIT=c_tmp1,FMT='(I15)') NINT(dt)
      ENDIF
    ELSE
      IF (dt < 1.) THEN
       WRITE (UNIT=c_tmp1,FMT='(F8.5)') dt
      ELSE
       WRITE (UNIT=c_tmp1,FMT='(F17.5)') dt
      ENDIF
      DO k=LEN_TRIM(c_tmp1),1,-1
        IF (c_tmp1(k:k) /= '0') THEN
          EXIT
        ELSE
          c_tmp1(k:k) = ' '
        ENDIF
      ENDDO
    ENDIF
    c_tmp2 = TRIM(c_tmp1)//' '//TRIM(c_tmp2)
    WRITE (UNIT=c_tmp3, &
 &   FMT='(A,I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
 &    TRIM(ADJUSTL(c_tmp2))//' since ',j_yy,j_mo,j_dd,j_hh,j_mn,j_ss
!---
    CALL flio_hdm (f_i,f_e,.TRUE.)
    i_rc = NF90_DEF_VAR(f_e,'time',NF90_REAL4, &
 &                      nw_di(k_1,f_i),timeid)
    i_rc = NF90_PUT_ATT(f_e,timeid,"axis",'T')
    i_rc = NF90_PUT_ATT(f_e,timeid,'standard_name','time')
    i_rc = NF90_PUT_ATT(f_e,timeid,'long_name','time steps')
    IF (PRESENT(t_calendar)) THEN
      i_rc = NF90_PUT_ATT(f_e,timeid,'calendar',TRIM(t_calendar))
    ENDIF
    i_rc = NF90_PUT_ATT(f_e,timeid,'units',TRIM(c_tmp3))
  ELSE IF (PRESENT(t_axis).OR.PRESENT(t_init).OR.PRESENT(t_step)) THEN
    CALL ipslerr (3,'fliopstc', &
 &   'For time axis and coordinates', &
 &   'arguments t_axis AND t_init AND t_step', &
 &   'must be PRESENT')
  ENDIF
!-
! Ensuring data mode
!-
  CALL flio_hdm (f_i,f_e,.FALSE.)
!-
! Create the longitude axis
!-
  IF (PRESENT(x_axis).OR.PRESENT(x_axis_2d)) THEN
    IF (l_dbg) THEN
      WRITE(*,*) '  fliopstc : Create the Longitude axis'
    ENDIF
    IF (PRESENT(x_axis)) THEN
      i_rc = NF90_PUT_VAR(f_e,lonid,x_axis(:))
    ELSE
      i_rc = NF90_PUT_VAR(f_e,lonid,x_axis_2d(:,:))
    ENDIF
  ENDIF
!-
! Create the Latitude axis
!-
  IF (PRESENT(y_axis).OR.PRESENT(y_axis_2d)) THEN
    IF (l_dbg) THEN
      WRITE(*,*) '  fliopstc : Create the Latitude axis'
    ENDIF
    IF (PRESENT(y_axis)) THEN
      i_rc = NF90_PUT_VAR(f_e,latid,y_axis(:))
    ELSE
      i_rc = NF90_PUT_VAR(f_e,latid,y_axis_2d(:,:))
    ENDIF
  ENDIF
!-
! Create the Vertical axis
!-
  IF (PRESENT(z_axis)) THEN
    IF (l_dbg) THEN
      WRITE(*,*) '  fliopstc : Create the Vertical axis'
    ENDIF
    i_rc = NF90_PUT_VAR(f_e,levid,z_axis(:))
  ENDIF
!-
! Create the Time axis
!-
  IF (PRESENT(t_axis)) THEN
    IF (l_dbg) THEN
      WRITE(*,*) '  fliopstc : Create the Time axis'
    ENDIF
    i_rc = NF90_PUT_VAR(f_e,timeid,REAL(t_axis(:)))
  ENDIF
!-
! Keep all this information
!-
  CALL flio_inf (f_e,nb_vars=nw_nv(f_i),nb_atts=nw_na(f_i))
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-fliopstc"
  ENDIF
!----------------------
END SUBROUTINE fliopstc
!===
SUBROUTINE fliodv_r0d &
 & (f_i,v_n,v_t, &
 &  axis,standard_name,long_name,units,valid_min,valid_max,fillvalue)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER,OPTIONAL,INTENT(IN) :: v_t
  CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: &
 & axis,standard_name,long_name,units
  REAL,OPTIONAL,INTENT(IN) :: valid_min,valid_max,fillvalue
!---------------------------------------------------------------------
  CALL flio_udv &
 &  (f_i,0,v_n,(/0/),v_t, &
 &   axis,standard_name,long_name,units,valid_min,valid_max,fillvalue)
!------------------------
END SUBROUTINE fliodv_r0d
!===
SUBROUTINE fliodv_rnd &
 & (f_i,v_n,v_d,v_t, &
 &  axis,standard_name,long_name,units,valid_min,valid_max,fillvalue)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER,DIMENSION(:),INTENT(IN) :: v_d
  INTEGER,OPTIONAL,INTENT(IN) :: v_t
  CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: &
 & axis,standard_name,long_name,units
  REAL,OPTIONAL,INTENT(IN) :: valid_min,valid_max,fillvalue
!---------------------------------------------------------------------
  CALL flio_udv &
 &  (f_i,SIZE(v_d),v_n,v_d,v_t, &
 &   axis,standard_name,long_name,units,valid_min,valid_max,fillvalue)
!------------------------
END SUBROUTINE fliodv_rnd
!===
SUBROUTINE flio_udv &
 & (f_i,n_d,v_n,v_d,v_t, &
 &  axis,standard_name,long_name,units,valid_min,valid_max,fillvalue)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i,n_d
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER,DIMENSION(:),INTENT(IN) :: v_d
  INTEGER,OPTIONAL,INTENT(IN) :: v_t
  CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: &
 & axis,standard_name,long_name,units
  REAL,OPTIONAL,INTENT(IN) :: valid_min,valid_max,fillvalue
!-
  INTEGER :: f_e,m_k,i_v,i_rc,ii,idd
  INTEGER,DIMENSION(nb_vd_mx) :: a_i
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->fliodefv ",TRIM(v_n)," ",n_d,"D"
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('fliodefv',f_i,f_e)
!-
  IF (n_d > 0) THEN
    IF (n_d > nb_vd_mx) THEN
      CALL ipslerr (3,'fliodefv', &
 &     'Too many dimensions', &
 &     'required for the variable',TRIM(v_n))
    ENDIF
  ENDIF
!-
  DO ii=1,n_d
    IF ( (v_d(ii) >= 1).AND.(v_d(ii) <= nb_fd_mx) ) THEN
      idd = nw_di(v_d(ii),f_i)
      IF (idd > 0) THEN
        a_i(ii) = idd
      ELSE
        CALL ipslerr (3,'fliodefv', &
 &       'Invalid dimension identifier','(not defined)',' ')
      ENDIF
    ELSE
      CALL ipslerr (3,'fliodefv', &
 &     'Invalid dimension identifier','(not supported)',' ')
    ENDIF
  ENDDO
!-
  i_rc = NF90_INQ_VARID(f_e,v_n,i_v)
  IF (i_rc /= NF90_NOERR) THEN
    CALL flio_hdm (f_i,f_e,.TRUE.)
!---
    IF (PRESENT(v_t)) THEN
      SELECT CASE (v_t)
      CASE(flio_i)
        IF (i_std == i_8) THEN
!-------- I8 not yet supported by NETCDF
!-------- m_k = flio_i8
          m_k = flio_i4
        ELSE
          m_k = flio_i4
        ENDIF
      CASE(flio_r)
        IF (r_std == r_8) THEN
          m_k = flio_r8
        ELSE
          m_k = flio_r4
        ENDIF
      CASE(flio_c,flio_i1,flio_i2,flio_i4,flio_r4,flio_r8)
        m_k = v_t
      CASE DEFAULT
        CALL ipslerr (3,'fliodefv', &
 &        'Variable '//TRIM(v_n),'External type','not supported')
      END SELECT
    ELSE IF (r_std == r_8) THEN
      m_k = flio_r8
    ELSE
      m_k = flio_r4
    ENDIF
!---
    IF (n_d > 0) THEN
      i_rc = NF90_DEF_VAR(f_e,v_n,m_k,a_i(1:n_d),i_v)
    ELSE
      i_rc = NF90_DEF_VAR(f_e,v_n,m_k,i_v)
    ENDIF
    IF (i_rc /= NF90_NOERR) THEN
      CALL ipslerr (3,'fliodefv', &
 &      'Variable '//TRIM(v_n)//' not defined','Error :', &
 &      TRIM(NF90_STRERROR(i_rc)))
    ENDIF
    nw_nv(f_i) = nw_nv(f_i)+1
!---
    IF (PRESENT(axis)) THEN
      i_rc = NF90_PUT_ATT(f_e,i_v,'axis',TRIM(axis))
    ENDIF
    IF (PRESENT(standard_name)) THEN
      i_rc = NF90_PUT_ATT(f_e,i_v,'standard_name',TRIM(standard_name))
    ENDIF
    IF (PRESENT(long_name)) THEN
      i_rc = NF90_PUT_ATT(f_e,i_v,'long_name',TRIM(long_name))
    ENDIF
    IF (PRESENT(units)) THEN
      i_rc = NF90_PUT_ATT(f_e,i_v,'units',TRIM(units))
    ENDIF
    IF (PRESENT(valid_min)) THEN
      SELECT CASE (m_k)
      CASE(flio_i1,flio_i2)
        i_rc = NF90_PUT_ATT(f_e,i_v,'valid_min',NINT(valid_min,KIND=i_2))
      CASE(flio_i4)
        i_rc = NF90_PUT_ATT(f_e,i_v,'valid_min',NINT(valid_min,KIND=i_4))
      CASE(flio_r4)
        i_rc = NF90_PUT_ATT(f_e,i_v,'valid_min',REAL(valid_min,KIND=r_4))
      CASE(flio_r8)
        i_rc = NF90_PUT_ATT(f_e,i_v,'valid_min',REAL(valid_min,KIND=r_8))
      CASE DEFAULT
        CALL ipslerr (2,'fliodefv', &
   &      'Variable '//TRIM(v_n),'attribute valid_min', &
   &      'not supported for this external type')
      END SELECT
    ENDIF
    IF (PRESENT(valid_max)) THEN
      SELECT CASE (m_k)
      CASE(flio_i1,flio_i2)
        i_rc = NF90_PUT_ATT(f_e,i_v,'valid_max',NINT(valid_max,KIND=i_2))
      CASE(flio_i4)
        i_rc = NF90_PUT_ATT(f_e,i_v,'valid_max',NINT(valid_max,KIND=i_4))
      CASE(flio_r4)
        i_rc = NF90_PUT_ATT(f_e,i_v,'valid_max',REAL(valid_max,KIND=r_4))
      CASE(flio_r8)
        i_rc = NF90_PUT_ATT(f_e,i_v,'valid_max',REAL(valid_max,KIND=r_8))
      CASE DEFAULT
        CALL ipslerr (2,'fliodefv', &
   &      'Variable '//TRIM(v_n),'attribute valid_max', &
   &      'not supported for this external type')
      END SELECT
    ENDIF
    IF (PRESENT(fillvalue)) THEN
      SELECT CASE (m_k)
      CASE(flio_i1,flio_i2)
        i_rc = NF90_PUT_ATT(f_e,i_v,'_FillValue',NINT(fillvalue,KIND=i_2))
      CASE(flio_i4)
        i_rc = NF90_PUT_ATT(f_e,i_v,'_FillValue',NINT(fillvalue,KIND=i_4))
      CASE(flio_r4)
        i_rc = NF90_PUT_ATT(f_e,i_v,'_FillValue',REAL(fillvalue,KIND=r_4))
      CASE(flio_r8)
        i_rc = NF90_PUT_ATT(f_e,i_v,'_FillValue',REAL(fillvalue,KIND=r_8))
      CASE DEFAULT
        CALL ipslerr (2,'fliodefv', &
   &      'Variable '//TRIM(v_n),'attribute fillvalue', &
   &      'not supported for this external type')
      END SELECT
    ENDIF
!---
  ELSE
    CALL ipslerr (3,'fliodefv','Variable',TRIM(v_n),'already exist')
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-fliodefv"
  ENDIF
!----------------------
END SUBROUTINE flio_udv
!===
SUBROUTINE fliopv_i40 (f_i,v_n,v_v,start)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,i_40=v_v,start=start)
!------------------------
END SUBROUTINE fliopv_i40
!===
SUBROUTINE fliopv_i41 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),DIMENSION(:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,i_41=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_i41
!===
SUBROUTINE fliopv_i42 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),DIMENSION(:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,i_42=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_i42
!===
SUBROUTINE fliopv_i43 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),DIMENSION(:,:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,i_43=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_i43
!===
SUBROUTINE fliopv_i44 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),DIMENSION(:,:,:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,i_44=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_i44
!===
SUBROUTINE fliopv_i45 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),DIMENSION(:,:,:,:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,i_45=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_i45
!===
SUBROUTINE fliopv_i20 (f_i,v_n,v_v,start)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_2),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,i_20=v_v,start=start)
!------------------------
END SUBROUTINE fliopv_i20
!===
SUBROUTINE fliopv_i21 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_2),DIMENSION(:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,i_21=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_i21
!===
SUBROUTINE fliopv_i22 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_2),DIMENSION(:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,i_22=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_i22
!===
SUBROUTINE fliopv_i23 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_2),DIMENSION(:,:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,i_23=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_i23
!===
SUBROUTINE fliopv_i24 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_2),DIMENSION(:,:,:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,i_24=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_i24
!===
SUBROUTINE fliopv_i25 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_2),DIMENSION(:,:,:,:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,i_25=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_i25
!===
!?INTEGERS of KIND 1 are not supported on all computers
!?SUBROUTINE fliopv_i10 (f_i,v_n,v_v,start)
!?!---------------------------------------------------------------------
!?  IMPLICIT NONE
!?!-
!?  INTEGER,INTENT(IN) :: f_i
!?  CHARACTER(LEN=*),INTENT(IN) :: v_n
!?  INTEGER(KIND=i_1),INTENT(IN) :: v_v
!?  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start
!?!---------------------------------------------------------------------
!?  CALL flio_upv (f_i,v_n,i_10=v_v,start=start)
!?!------------------------
!?END SUBROUTINE fliopv_i10
!?!===
!?SUBROUTINE fliopv_i11 (f_i,v_n,v_v,start,count)
!?!---------------------------------------------------------------------
!?  IMPLICIT NONE
!?!-
!?  INTEGER,INTENT(IN) :: f_i
!?  CHARACTER(LEN=*),INTENT(IN) :: v_n
!?  INTEGER(KIND=i_1),DIMENSION(:),INTENT(IN) :: v_v
!?  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!?!---------------------------------------------------------------------
!?  CALL flio_upv (f_i,v_n,i_11=v_v,start=start,count=count)
!?!------------------------
!?END SUBROUTINE fliopv_i11
!?!===
!?SUBROUTINE fliopv_i12 (f_i,v_n,v_v,start,count)
!?!---------------------------------------------------------------------
!?  IMPLICIT NONE
!?!-
!?  INTEGER,INTENT(IN) :: f_i
!?  CHARACTER(LEN=*),INTENT(IN) :: v_n
!?  INTEGER(KIND=i_1),DIMENSION(:,:),INTENT(IN) :: v_v
!?  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!?!---------------------------------------------------------------------
!?  CALL flio_upv (f_i,v_n,i_12=v_v,start=start,count=count)
!?!------------------------
!?END SUBROUTINE fliopv_i12
!?!===
!?SUBROUTINE fliopv_i13 (f_i,v_n,v_v,start,count)
!?!---------------------------------------------------------------------
!?  IMPLICIT NONE
!?!-
!?  INTEGER,INTENT(IN) :: f_i
!?  CHARACTER(LEN=*),INTENT(IN) :: v_n
!?  INTEGER(KIND=i_1),DIMENSION(:,:,:),INTENT(IN) :: v_v
!?  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!?!---------------------------------------------------------------------
!?  CALL flio_upv (f_i,v_n,i_13=v_v,start=start,count=count)
!?!------------------------
!?END SUBROUTINE fliopv_i13
!?!===
!?SUBROUTINE fliopv_i14 (f_i,v_n,v_v,start,count)
!?!---------------------------------------------------------------------
!?  IMPLICIT NONE
!?!-
!?  INTEGER,INTENT(IN) :: f_i
!?  CHARACTER(LEN=*),INTENT(IN) :: v_n
!?  INTEGER(KIND=i_1),DIMENSION(:,:,:,:),INTENT(IN) :: v_v
!?  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!?!---------------------------------------------------------------------
!?  CALL flio_upv (f_i,v_n,i_14=v_v,start=start,count=count)
!?!------------------------
!?END SUBROUTINE fliopv_i14
!?!===
!?SUBROUTINE fliopv_i15 (f_i,v_n,v_v,start,count)
!?!---------------------------------------------------------------------
!?  IMPLICIT NONE
!?!-
!?  INTEGER,INTENT(IN) :: f_i
!?  CHARACTER(LEN=*),INTENT(IN) :: v_n
!?  INTEGER(KIND=i_1),DIMENSION(:,:,:,:,:),INTENT(IN) :: v_v
!?  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!?!---------------------------------------------------------------------
!?  CALL flio_upv (f_i,v_n,i_15=v_v,start=start,count=count)
!?!------------------------
!?END SUBROUTINE fliopv_i15
!===
SUBROUTINE fliopv_r40 (f_i,v_n,v_v,start)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_4),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,r_40=v_v,start=start)
!------------------------
END SUBROUTINE fliopv_r40
!===
SUBROUTINE fliopv_r41 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_4),DIMENSION(:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,r_41=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_r41
!===
SUBROUTINE fliopv_r42 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_4),DIMENSION(:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,r_42=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_r42
!===
SUBROUTINE fliopv_r43 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_4),DIMENSION(:,:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,r_43=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_r43
!===
SUBROUTINE fliopv_r44 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_4),DIMENSION(:,:,:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,r_44=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_r44
!===
SUBROUTINE fliopv_r45 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_4),DIMENSION(:,:,:,:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,r_45=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_r45
!===
SUBROUTINE fliopv_r80 (f_i,v_n,v_v,start)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_8),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,r_80=v_v,start=start)
!------------------------
END SUBROUTINE fliopv_r80
!===
SUBROUTINE fliopv_r81 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_8),DIMENSION(:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,r_81=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_r81
!===
SUBROUTINE fliopv_r82 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_8),DIMENSION(:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,r_82=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_r82
!===
SUBROUTINE fliopv_r83 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_8),DIMENSION(:,:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,r_83=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_r83
!===
SUBROUTINE fliopv_r84 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_8),DIMENSION(:,:,:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,r_84=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_r84
!===
SUBROUTINE fliopv_r85 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_8),DIMENSION(:,:,:,:,:),INTENT(IN) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_upv (f_i,v_n,r_85=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliopv_r85
!===
SUBROUTINE flio_upv &
 & (f_i,v_n, &
 &  i_40,i_41,i_42,i_43,i_44,i_45, &
 &  i_20,i_21,i_22,i_23,i_24,i_25, &
!? &  i_10,i_11,i_12,i_13,i_14,i_15, &
 &  r_40,r_41,r_42,r_43,r_44,r_45, &
 &  r_80,r_81,r_82,r_83,r_84,r_85, &
 &  start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),INTENT(IN),OPTIONAL :: i_40
  INTEGER(KIND=i_4),DIMENSION(:),INTENT(IN),OPTIONAL :: i_41
  INTEGER(KIND=i_4),DIMENSION(:,:),INTENT(IN),OPTIONAL :: i_42
  INTEGER(KIND=i_4),DIMENSION(:,:,:),INTENT(IN),OPTIONAL :: i_43
  INTEGER(KIND=i_4),DIMENSION(:,:,:,:),INTENT(IN),OPTIONAL :: i_44
  INTEGER(KIND=i_4),DIMENSION(:,:,:,:,:),INTENT(IN),OPTIONAL :: i_45
  INTEGER(KIND=i_2),INTENT(IN),OPTIONAL :: i_20
  INTEGER(KIND=i_2),DIMENSION(:),INTENT(IN),OPTIONAL :: i_21
  INTEGER(KIND=i_2),DIMENSION(:,:),INTENT(IN),OPTIONAL :: i_22
  INTEGER(KIND=i_2),DIMENSION(:,:,:),INTENT(IN),OPTIONAL :: i_23
  INTEGER(KIND=i_2),DIMENSION(:,:,:,:),INTENT(IN),OPTIONAL :: i_24
  INTEGER(KIND=i_2),DIMENSION(:,:,:,:,:),INTENT(IN),OPTIONAL :: i_25
!?INTEGERS of KIND 1 are not supported on all computers
!?INTEGER(KIND=i_1),INTENT(IN),OPTIONAL :: i_10
!?INTEGER(KIND=i_1),DIMENSION(:),INTENT(IN),OPTIONAL :: i_11
!?INTEGER(KIND=i_1),DIMENSION(:,:),INTENT(IN),OPTIONAL :: i_12
!?INTEGER(KIND=i_1),DIMENSION(:,:,:),INTENT(IN),OPTIONAL :: i_13
!?INTEGER(KIND=i_1),DIMENSION(:,:,:,:),INTENT(IN),OPTIONAL :: i_14
!?INTEGER(KIND=i_1),DIMENSION(:,:,:,:,:),INTENT(IN),OPTIONAL :: i_15
  REAL(KIND=r_4),INTENT(IN),OPTIONAL :: r_40
  REAL(KIND=r_4),DIMENSION(:),INTENT(IN),OPTIONAL :: r_41
  REAL(KIND=r_4),DIMENSION(:,:),INTENT(IN),OPTIONAL :: r_42
  REAL(KIND=r_4),DIMENSION(:,:,:),INTENT(IN),OPTIONAL :: r_43
  REAL(KIND=r_4),DIMENSION(:,:,:,:),INTENT(IN),OPTIONAL :: r_44
  REAL(KIND=r_4),DIMENSION(:,:,:,:,:),INTENT(IN),OPTIONAL :: r_45
  REAL(KIND=r_8),INTENT(IN),OPTIONAL :: r_80
  REAL(KIND=r_8),DIMENSION(:),INTENT(IN),OPTIONAL :: r_81
  REAL(KIND=r_8),DIMENSION(:,:),INTENT(IN),OPTIONAL :: r_82
  REAL(KIND=r_8),DIMENSION(:,:,:),INTENT(IN),OPTIONAL :: r_83
  REAL(KIND=r_8),DIMENSION(:,:,:,:),INTENT(IN),OPTIONAL :: r_84
  REAL(KIND=r_8),DIMENSION(:,:,:,:,:),INTENT(IN),OPTIONAL :: r_85
  INTEGER,DIMENSION(:),INTENT(IN),OPTIONAL :: start,count
!-
  INTEGER :: f_e,i_v,i_rc
  CHARACTER(LEN=5) :: cvr_d
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    IF      (PRESENT(i_40)) THEN; cvr_d = "I1 0D";
    ELSE IF (PRESENT(i_41)) THEN; cvr_d = "I1 1D";
    ELSE IF (PRESENT(i_42)) THEN; cvr_d = "I1 2D";
    ELSE IF (PRESENT(i_43)) THEN; cvr_d = "I1 3D";
    ELSE IF (PRESENT(i_44)) THEN; cvr_d = "I1 4D";
    ELSE IF (PRESENT(i_45)) THEN; cvr_d = "I1 5D";
    ELSE IF (PRESENT(i_20)) THEN; cvr_d = "I2 0D";
    ELSE IF (PRESENT(i_21)) THEN; cvr_d = "I2 1D";
    ELSE IF (PRESENT(i_22)) THEN; cvr_d = "I2 2D";
    ELSE IF (PRESENT(i_23)) THEN; cvr_d = "I2 3D";
    ELSE IF (PRESENT(i_24)) THEN; cvr_d = "I2 4D";
    ELSE IF (PRESENT(i_25)) THEN; cvr_d = "I2 5D";
!?  ELSE IF (PRESENT(i_10)) THEN; cvr_d = "I4 0D";
!?  ELSE IF (PRESENT(i_11)) THEN; cvr_d = "I4 1D";
!?  ELSE IF (PRESENT(i_12)) THEN; cvr_d = "I4 2D";
!?  ELSE IF (PRESENT(i_13)) THEN; cvr_d = "I4 3D";
!?  ELSE IF (PRESENT(i_14)) THEN; cvr_d = "I4 4D";
!?  ELSE IF (PRESENT(i_15)) THEN; cvr_d = "I4 5D";
    ELSE IF (PRESENT(r_40)) THEN; cvr_d = "R4 0D";
    ELSE IF (PRESENT(r_41)) THEN; cvr_d = "R4 1D";
    ELSE IF (PRESENT(r_42)) THEN; cvr_d = "R4 2D";
    ELSE IF (PRESENT(r_43)) THEN; cvr_d = "R4 3D";
    ELSE IF (PRESENT(r_44)) THEN; cvr_d = "R4 4D";
    ELSE IF (PRESENT(r_45)) THEN; cvr_d = "R4 5D";
    ELSE IF (PRESENT(r_80)) THEN; cvr_d = "R8 0D";
    ELSE IF (PRESENT(r_81)) THEN; cvr_d = "R8 1D";
    ELSE IF (PRESENT(r_82)) THEN; cvr_d = "R8 2D";
    ELSE IF (PRESENT(r_83)) THEN; cvr_d = "R8 3D";
    ELSE IF (PRESENT(r_84)) THEN; cvr_d = "R8 4D";
    ELSE IF (PRESENT(r_85)) THEN; cvr_d = "R8 5D";
    ENDIF
    WRITE(*,*) "->flioputv ",TRIM(v_n)," ",TRIM(cvr_d)
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('flioputv',f_i,f_e)
!-
! Ensuring data mode
!-
  CALL flio_hdm (f_i,f_e,.FALSE.)
!-
  i_rc = NF90_INQ_VARID(f_e,v_n,i_v)
  IF (i_rc == NF90_NOERR) THEN
    IF      (PRESENT(i_40)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,i_40,start=start)
    ELSE IF (PRESENT(i_41)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,i_41,start=start,count=count)
    ELSE IF (PRESENT(i_42)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,i_42,start=start,count=count)
    ELSE IF (PRESENT(i_43)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,i_43,start=start,count=count)
    ELSE IF (PRESENT(i_44)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,i_44,start=start,count=count)
    ELSE IF (PRESENT(i_45)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,i_45,start=start,count=count)
    ELSE IF (PRESENT(i_20)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,i_20,start=start)
    ELSE IF (PRESENT(i_21)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,i_21,start=start,count=count)
    ELSE IF (PRESENT(i_22)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,i_22,start=start,count=count)
    ELSE IF (PRESENT(i_23)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,i_23,start=start,count=count)
    ELSE IF (PRESENT(i_24)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,i_24,start=start,count=count)
    ELSE IF (PRESENT(i_25)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,i_25,start=start,count=count)
!?  ELSE IF (PRESENT(i_10)) THEN
!?    i_rc = NF90_PUT_VAR(f_e,i_v,i_10,start=start)
!?  ELSE IF (PRESENT(i_11)) THEN
!?    i_rc = NF90_PUT_VAR(f_e,i_v,i_11,start=start,count=count)
!?  ELSE IF (PRESENT(i_12)) THEN
!?    i_rc = NF90_PUT_VAR(f_e,i_v,i_12,start=start,count=count)
!?  ELSE IF (PRESENT(i_13)) THEN
!?    i_rc = NF90_PUT_VAR(f_e,i_v,i_13,start=start,count=count)
!?  ELSE IF (PRESENT(i_14)) THEN
!?    i_rc = NF90_PUT_VAR(f_e,i_v,i_14,start=start,count=count)
!?  ELSE IF (PRESENT(i_15)) THEN
!?    i_rc = NF90_PUT_VAR(f_e,i_v,i_15,start=start,count=count)
    ELSE IF (PRESENT(r_40)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,r_40,start=start)
    ELSE IF (PRESENT(r_41)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,r_41,start=start,count=count)
    ELSE IF (PRESENT(r_42)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,r_42,start=start,count=count)
    ELSE IF (PRESENT(r_43)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,r_43,start=start,count=count)
    ELSE IF (PRESENT(r_44)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,r_44,start=start,count=count)
    ELSE IF (PRESENT(r_45)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,r_45,start=start,count=count)
    ELSE IF (PRESENT(r_80)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,r_80,start=start)
    ELSE IF (PRESENT(r_81)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,r_81,start=start,count=count)
    ELSE IF (PRESENT(r_82)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,r_82,start=start,count=count)
    ELSE IF (PRESENT(r_83)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,r_83,start=start,count=count)
    ELSE IF (PRESENT(r_84)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,r_84,start=start,count=count)
    ELSE IF (PRESENT(r_85)) THEN
      i_rc = NF90_PUT_VAR(f_e,i_v,r_85,start=start,count=count)
    ENDIF
    IF (i_rc /= NF90_NOERR) THEN
      CALL ipslerr (3,'flioputv', &
 &      'Variable '//TRIM(v_n)//' not put','Error :', &
 &      TRIM(NF90_STRERROR(i_rc)))
    ENDIF
  ELSE
    CALL ipslerr (3,'flioputv','Variable',TRIM(v_n),'not defined')
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-flioputv"
  ENDIF
!----------------------
END SUBROUTINE flio_upv
!===
SUBROUTINE fliopa_r4_0d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  REAL(KIND=4),INTENT(IN) :: a_v
!---------------------------------------------------------------------
  CALL flio_upa (f_i,1,v_n,a_n,avr4=(/a_v/))
!--------------------------
END SUBROUTINE fliopa_r4_0d
!===
SUBROUTINE fliopa_r4_1d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  REAL(KIND=4),DIMENSION(:),INTENT(IN) :: a_v
!---------------------------------------------------------------------
  CALL flio_upa (f_i,SIZE(a_v),v_n,a_n,avr4=a_v)
!--------------------------
END SUBROUTINE fliopa_r4_1d
!===
SUBROUTINE fliopa_r8_0d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  REAL(KIND=8),INTENT(IN) :: a_v
!---------------------------------------------------------------------
  CALL flio_upa (f_i,1,v_n,a_n,avr8=(/a_v/))
!--------------------------
END SUBROUTINE fliopa_r8_0d
!===
SUBROUTINE fliopa_r8_1d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  REAL(KIND=8),DIMENSION(:),INTENT(IN) :: a_v
!---------------------------------------------------------------------
  CALL flio_upa (f_i,SIZE(a_v),v_n,a_n,avr8=a_v)
!--------------------------
END SUBROUTINE fliopa_r8_1d
!===
SUBROUTINE fliopa_i4_0d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  INTEGER(KIND=4),INTENT(IN) :: a_v
!---------------------------------------------------------------------
  CALL flio_upa (f_i,1,v_n,a_n,avi4=(/a_v/))
!--------------------------
END SUBROUTINE fliopa_i4_0d
!===
SUBROUTINE fliopa_i4_1d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  INTEGER(KIND=4),DIMENSION(:),INTENT(IN) :: a_v
!---------------------------------------------------------------------
  CALL flio_upa (f_i,SIZE(a_v),v_n,a_n,avi4=a_v)
!--------------------------
END SUBROUTINE fliopa_i4_1d
!===
SUBROUTINE fliopa_tx_0d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  CHARACTER(LEN=*),INTENT(IN) :: a_v
!---------------------------------------------------------------------
  CALL flio_upa (f_i,1,v_n,a_n,avtx=a_v)
!--------------------------
END SUBROUTINE fliopa_tx_0d
!===
SUBROUTINE flio_upa (f_i,l_a,v_n,a_n,avr4,avr8,avi4,avtx)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i,l_a
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  REAL(KIND=4),DIMENSION(:),OPTIONAL,INTENT(IN) :: avr4
  REAL(KIND=8),DIMENSION(:),OPTIONAL,INTENT(IN) :: avr8
  INTEGER(KIND=4),DIMENSION(:),OPTIONAL,INTENT(IN) :: avi4
  CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: avtx
!-
  INTEGER :: f_e,i_v,i_a,i_rc
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->flioputa ",TRIM(v_n)," ",TRIM(a_n)
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('flioputa',f_i,f_e)
!-
  IF (TRIM(v_n) == '?') THEN
    i_v = NF90_GLOBAL
  ELSE
    i_rc = NF90_INQ_VARID(f_e,v_n,i_v)
    IF (i_rc /= NF90_NOERR) THEN
      CALL ipslerr (3,'flioputa', &
       'Variable :',TRIM(v_n),'not found')
    ENDIF
  ENDIF
!-
  i_rc = NF90_INQUIRE_ATTRIBUTE(f_e,i_v,a_n,attnum=i_a)
  IF ( (i_v == NF90_GLOBAL).AND.(i_rc /= NF90_NOERR) ) THEN
    nw_na(f_i) = nw_na(f_i)+1
  ENDIF
  CALL flio_hdm (f_i,f_e,.TRUE.)
  IF      (PRESENT(avr4)) THEN
    i_rc = NF90_PUT_ATT(f_e,i_v,a_n,avr4(1:l_a))
  ELSE IF (PRESENT(avr8)) THEN
    i_rc = NF90_PUT_ATT(f_e,i_v,a_n,avr8(1:l_a))
  ELSE IF (PRESENT(avi4)) THEN
    i_rc = NF90_PUT_ATT(f_e,i_v,a_n,avi4(1:l_a))
  ELSE IF (PRESENT(avtx)) THEN
    i_rc = NF90_PUT_ATT(f_e,i_v,a_n,TRIM(avtx))
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-flioputa"
  ENDIF
!----------------------
END SUBROUTINE flio_upa
!===
SUBROUTINE flioopfd (f_n,f_i,mode,nb_dim,nb_var,nb_gat)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(IN) :: f_n
  INTEGER,INTENT(OUT) :: f_i
  CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: mode
  INTEGER,OPTIONAL,INTENT(OUT) :: nb_dim,nb_var,nb_gat
!-
  INTEGER :: i_rc,f_e,m_c
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) '->flioopfd, file name : ',TRIM(f_n)
  ENDIF
!-
! Search for a free local identifier
!-
  f_i = flio_rid()
  IF (f_i < 0) THEN
    CALL ipslerr (3,'flioopfd', &
      'Too many files.','Please increase nb_fi_mx', &
      'in module fliocom.f90.')
  ENDIF
!-
! Check the mode
!-
  IF (PRESENT(mode)) THEN
    IF (TRIM(mode) == "WRITE") THEN
      m_c = NF90_WRITE
    ELSE
      m_c = NF90_NOWRITE
    ENDIF
  ELSE
    m_c = NF90_NOWRITE
  ENDIF
!-
! Open the file.
!-
  i_rc = NF90_OPEN(TRIM(f_n),m_c,f_e)
  IF (i_rc /= NF90_NOERR) THEN
    CALL ipslerr (3,'flioopfd', &
 &   'Could not open file :',TRIM(f_n), &
 &   TRIM(NF90_STRERROR(i_rc))//' (Netcdf)')
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) '  flioopfd, model file-id : ',f_e
  ENDIF
!-
! Retrieve and keep information about the file
!-
  nw_id(f_i) = f_e
  lw_hm(f_i) = .FALSE.
  CALL flio_inf (f_e, &
 &  nb_dims=nw_nd(f_i),nb_vars=nw_nv(f_i), &
 &  nb_atts=nw_na(f_i),id_unlm=nw_un(f_i), &
 &  nn_idm=nw_di(:,f_i),nn_ldm=nw_dl(:,f_i),nn_aid=nw_ai(:,f_i))
!-
! Return information to the user
!-
  IF (PRESENT(nb_dim)) THEN
    nb_dim = nw_nd(f_i)
  ENDIF
  IF (PRESENT(nb_var)) THEN
    nb_var = nw_nv(f_i)
  ENDIF
  IF (PRESENT(nb_gat)) THEN
    nb_gat = nw_na(f_i)
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,'("   flioopfd - dimensions :",/,(5(1X,I10),:))') &
 &    nw_dl(:,f_i)
    WRITE(*,*) "<-flioopfd"
  ENDIF
!----------------------
END SUBROUTINE flioopfd
!===
SUBROUTINE flioinqf &
 & (f_i,nb_dim,nb_var,nb_gat,id_uld,id_dim,ln_dim)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  INTEGER,OPTIONAL,INTENT(OUT) :: nb_dim,nb_var,nb_gat,id_uld
  INTEGER,OPTIONAL,INTENT(OUT),DIMENSION(:) :: id_dim,ln_dim
!-
  INTEGER :: lll
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->flioinqf"
  ENDIF
!-
  IF ( (f_i < 1).OR.(f_i > nb_fi_mx) ) THEN
    CALL ipslerr (2,'flioinqf', &
 &   'Invalid file identifier',' ',' ')
  ELSE IF (nw_id(f_i) <= 0) THEN
    CALL ipslerr (2,'flioinqf', &
 &   'Unable to inquire about the file :','probably','not opened')
  ELSE
    IF (PRESENT(nb_dim)) THEN
      nb_dim = nw_nd(f_i)
    ENDIF
    IF (PRESENT(nb_var)) THEN
      nb_var = nw_nv(f_i)
    ENDIF
    IF (PRESENT(nb_gat)) THEN
      nb_gat = nw_na(f_i)
    ENDIF
    IF (PRESENT(id_uld)) THEN
      id_uld = nw_un(f_i)
    ENDIF
    IF (PRESENT(id_dim)) THEN
      lll = SIZE(id_dim)
      IF (lll < nw_nd(f_i)) THEN
        CALL ipslerr (2,'flioinqf', &
 &       'Only the first identifiers', &
 &       'of the dimensions','will be returned')
      ENDIF
      lll=MIN(SIZE(id_dim),nw_nd(f_i))
      id_dim(1:lll) = nw_di(1:lll,f_i)
    ENDIF
    IF (PRESENT(ln_dim)) THEN
      lll = SIZE(ln_dim)
      IF (lll < nw_nd(f_i)) THEN
        CALL ipslerr (2,'flioinqf', &
 &       'Only the first lengths', &
 &       'of the dimensions','will be returned')
      ENDIF
      lll=MIN(SIZE(ln_dim),nw_nd(f_i))
      ln_dim(1:lll) = nw_dl(1:lll,f_i)
    ENDIF
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-flioinqf"
  ENDIF
!----------------------
END SUBROUTINE flioinqf
!===
SUBROUTINE flioinqn &
 & (f_i,cn_dim,cn_var,cn_gat,cn_uld, &
 &  id_start,id_count,iv_start,iv_count,ia_start,ia_count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),DIMENSION(:),OPTIONAL,INTENT(OUT) :: &
 & cn_dim,cn_var,cn_gat
  CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: &
 & cn_uld
  INTEGER,OPTIONAL,INTENT(IN) :: &
 & id_start,id_count,iv_start,iv_count,ia_start,ia_count
!-
  INTEGER :: f_e,i_s,i_w,iws,iwc,i_rc
  LOGICAL :: l_ok
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->flioinqn"
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('flioinqn',f_i,f_e)
!-
  IF (PRESENT(cn_dim)) THEN
    l_ok = .TRUE.
    i_s = SIZE(cn_dim)
    DO i_w=1,i_s
      cn_dim(i_w)(:) = '?'
    ENDDO
    IF (PRESENT(id_start)) THEN
      iws = id_start
    ELSE
      iws = 1
    ENDIF
    IF (PRESENT(id_count)) THEN
      iwc = id_count
    ELSE
      iwc = nw_nd(f_i)
    ENDIF
    IF (iws > nw_nd(f_i)) THEN
      l_ok = .FALSE.
      CALL ipslerr (2,'flioinqn', &
 &     'The start index of requested dimensions', &
 &     'is greater than the number of dimensions', &
 &     'in the file')
    ELSE IF (iws < 1) THEN
      l_ok = .FALSE.
      CALL ipslerr (2,'flioinqn', &
 &     'The start index of requested dimensions', &
 &     'is invalid', &
 &     '( < 1 )')
    ENDIF
    IF ((iws+iwc-1) > nw_nd(f_i)) THEN
      CALL ipslerr (2,'flioinqn', &
 &     'The number of requested dimensions', &
 &     'is greater than the number of dimensions', &
 &     'in the file')
    ENDIF
    IF (iwc > i_s) THEN
      CALL ipslerr (2,'flioinqn', &
 &     'The number of dimensions to retrieve', &
 &     'is greater than the size of the array,', &
 &     'only the first dimensions of the file will be returned')
    ELSE IF (iwc < 1) THEN
      l_ok = .FALSE.
      CALL ipslerr (2,'flioinqn', &
 &     'The number of requested dimensions', &
 &     'is invalid', &
 &     '( < 1 )')
    ENDIF
    IF (l_ok) THEN
      DO i_w=1,MIN(iwc,i_s,nw_nd(f_i)-iws+1)
        i_rc = NF90_INQUIRE_DIMENSION(f_e,i_w+iws-1,name=cn_dim(i_w))
      ENDDO
    ENDIF
  ENDIF
!-
  IF (PRESENT(cn_var)) THEN
    l_ok = .TRUE.
    i_s = SIZE(cn_var)
    DO i_w=1,i_s
      cn_var(i_w)(:) = '?'
    ENDDO
    IF (PRESENT(iv_start)) THEN
      iws = iv_start
    ELSE
      iws = 1
    ENDIF
    IF (PRESENT(iv_count)) THEN
      iwc = iv_count
    ELSE
      iwc = nw_nv(f_i)
    ENDIF
    IF (iws > nw_nv(f_i)) THEN
      l_ok = .FALSE.
      CALL ipslerr (2,'flioinqn', &
 &     'The start index of requested variables', &
 &     'is greater than the number of variables', &
 &     'in the file')
    ELSE IF (iws < 1) THEN
      l_ok = .FALSE.
      CALL ipslerr (2,'flioinqn', &
 &     'The start index of requested variables', &
 &     'is invalid', &
 &     '( < 1 )')
    ENDIF
    IF ((iws+iwc-1) > nw_nv(f_i)) THEN
      CALL ipslerr (2,'flioinqn', &
 &     'The number of requested variables', &
 &     'is greater than the number of variables', &
 &     'in the file')
    ENDIF
    IF (iwc > i_s) THEN
      CALL ipslerr (2,'flioinqn', &
 &     'The number of variables to retrieve', &
 &     'is greater than the size of the array,', &
 &     'only the first variables of the file will be returned')
    ELSE IF (iwc < 1) THEN
      l_ok = .FALSE.
      CALL ipslerr (2,'flioinqn', &
 &     'The number of requested variables', &
 &     'is invalid', &
 &     '( < 1 )')
    ENDIF
    IF (l_ok) THEN
      DO i_w=1,MIN(iwc,i_s,nw_nv(f_i)-iws+1)
        i_rc = NF90_INQUIRE_VARIABLE(f_e,i_w+iws-1,name=cn_var(i_w))
      ENDDO
    ENDIF
  ENDIF
!-
  IF (PRESENT(cn_gat)) THEN
    l_ok = .TRUE.
    i_s = SIZE(cn_gat)
    DO i_w=1,i_s
      cn_gat(i_w)(:) = '?'
    ENDDO
    IF (PRESENT(ia_start)) THEN
      iws = ia_start
    ELSE
      iws = 1
    ENDIF
    IF (PRESENT(ia_count)) THEN
      iwc = ia_count
    ELSE
      iwc = nw_na(f_i)
    ENDIF
    IF (iws > nw_na(f_i)) THEN
      l_ok = .FALSE.
      CALL ipslerr (2,'flioinqn', &
 &     'The start index of requested global attributes', &
 &     'is greater than the number of global attributes', &
 &     'in the file')
    ELSE IF (iws < 1) THEN
      l_ok = .FALSE.
      CALL ipslerr (2,'flioinqn', &
 &     'The start index of requested global attributes', &
 &     'is invalid', &
 &     '( < 1 )')
    ENDIF
    IF ((iws+iwc-1) > nw_na(f_i)) THEN
      CALL ipslerr (2,'flioinqn', &
 &     'The number of requested global attributes', &
 &     'is greater than the number of global attributes', &
 &     'in the file')
    ENDIF
    IF (iwc > i_s) THEN
      CALL ipslerr (2,'flioinqn', &
 &     'The number of global attributes to retrieve', &
 &     'is greater than the size of the array,', &
 &     'only the first global attributes of the file will be returned')
    ELSE IF (iwc < 1) THEN
      l_ok = .FALSE.
      CALL ipslerr (2,'flioinqn', &
 &     'The number of requested global attributes', &
 &     'is invalid', &
 &     '( < 1 )')
    ENDIF
    IF (l_ok) THEN
      DO i_w=1,MIN(iwc,i_s,nw_na(f_i)-iws+1)
        i_rc = NF90_INQ_ATTNAME(f_e, &
 &              NF90_GLOBAL,i_w+iws-1,name=cn_gat(i_w))
      ENDDO
    ENDIF
  ENDIF
!-
  IF (PRESENT(cn_uld)) THEN
    cn_uld = '?'
    IF (nw_un(f_i) > 0) THEN
      i_rc = NF90_INQUIRE_DIMENSION(f_e,nw_un(f_i),name=cn_uld)
    ENDIF
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-flioinqn"
  ENDIF
!----------------------
END SUBROUTINE flioinqn
!===
SUBROUTINE fliogstc &
 & (f_i,x_axis,x_axis_2d,y_axis,y_axis_2d,z_axis, &
 &      t_axis,t_init,t_step,t_calendar, &
 &      x_start,x_count,y_start,y_count, &
 &      z_start,z_count,t_start,t_count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  REAL,DIMENSION(:),OPTIONAL,INTENT(OUT)    :: x_axis,y_axis
  REAL,DIMENSION(:,:),OPTIONAL,INTENT(OUT)  :: x_axis_2d,y_axis_2d
  REAL,DIMENSION(:),OPTIONAL,INTENT(OUT)    :: z_axis
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(OUT) :: t_axis
  REAL,OPTIONAL,INTENT(OUT)                 :: t_init,t_step
  CHARACTER(LEN=*),OPTIONAL,INTENT(OUT)     :: t_calendar
  INTEGER,OPTIONAL,INTENT(IN) :: &
 &  x_start,x_count,y_start,y_count,z_start,z_count,t_start,t_count
!-
  INTEGER :: i_rc,f_e,i_v,it_t,nbdim,kv
  INTEGER :: m_x,i_x,l_x,m_y,i_y,l_y,m_z,i_z,l_z,m_t,i_t,l_t
  CHARACTER(LEN=NF90_MAX_NAME) :: name
  CHARACTER(LEN=80) :: units
  CHARACTER(LEN=20) :: c_tmp
  CHARACTER(LEN=1) :: c_1
  REAL    :: r_yy,r_mo,r_dd,r_ss,dtv,dtn
  INTEGER :: j_yy,j_mo,j_dd,j_hh,j_mn,j_ss
  LOGICAL :: l_ok,l_tmp
!-
  REAL,DIMENSION(:),ALLOCATABLE :: v_tmp
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->fliogstc"
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('fliogstc',f_i,f_e)
!-
! Validate the coherence of the arguments
!-
  IF (    (PRESENT(x_axis).AND.PRESENT(x_axis_2d)) &
 &    .OR.(PRESENT(y_axis).AND.PRESENT(y_axis_2d)) ) THEN
    CALL ipslerr (3,'fliogstc', &
 &    'The [x/y]_axis arguments', &
 &    'are not coherent :',&
 &    'can not handle two [x/y]_axis')
  ENDIF
!-
! Retrieve spatio-temporal dimensions
!-
  IF (nw_ai(k_lon,f_i) > 0) THEN
    m_x = nw_dl(nw_ai(k_lon,f_i),f_i);
  ELSE
    m_x = -1;
  ENDIF
  IF (nw_ai(k_lat,f_i) > 0) THEN
    m_y = nw_dl(nw_ai(k_lat,f_i),f_i);
  ELSE
    m_y = -1;
  ENDIF
  IF (nw_ai(k_lev,f_i) > 0) THEN
    m_z = nw_dl(nw_ai(k_lev,f_i),f_i);
  ELSE
    m_z = -1;
  ENDIF
  IF (nw_ai(k_tim,f_i) > 0) THEN
    m_t = nw_dl(nw_ai(k_tim,f_i),f_i);
  ELSE
    m_t = -1;
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,'("   fliogstc - dimensions :",/,(5(1X,I10),:))') &
 &    m_x,m_y,m_z,m_t
  ENDIF
!-
! Initialize the x-y indices
!-
  IF (    PRESENT(x_axis)    &
 &    .OR.PRESENT(x_axis_2d) &
 &    .OR.PRESENT(y_axis_2d) ) THEN
    IF (PRESENT(x_start)) THEN
      i_x = x_start
    ELSE
      i_x = 1
    ENDIF
    IF (PRESENT(x_count)) THEN
      l_x = x_count
    ELSE
      l_x = m_x-i_x+1
    ENDIF
  ENDIF
  IF (    PRESENT(y_axis)    &
 &    .OR.PRESENT(y_axis_2d) &
 &    .OR.PRESENT(x_axis_2d) ) THEN
    IF (PRESENT(y_start)) THEN
      i_y = y_start
    ELSE
      i_y = 1
    ENDIF
    IF (PRESENT(y_count)) THEN
      l_y = y_count
    ELSE
      l_y = m_y-i_y+1
    ENDIF
  ENDIF
  IF (PRESENT(x_axis)) THEN
    IF (m_x <= 0) THEN
      CALL ipslerr (3,'fliogstc', &
 &      'Requested x_axis', &
 &      'but the coordinate is not present','in the file')
    ELSE IF ((i_x+l_x-1) > m_x) THEN
      CALL ipslerr (3,'fliogstc', &
 &      'The requested size for the x_axis', &
 &      'is greater than the size of the coordinate','in the file')
    ENDIF
  ENDIF
  IF (PRESENT(y_axis)) THEN
    IF (m_y <= 0) THEN
      CALL ipslerr (3,'fliogstc', &
 &      'Requested y_axis', &
 &      'but the coordinate is not present','in the file')
    ELSE IF ((i_y+l_y-1) > m_y) THEN
      CALL ipslerr (3,'fliogstc', &
 &      'The requested size for the y_axis', &
 &      'is greater than the size of the coordinate','in the file')
    ENDIF
  ENDIF
  IF (PRESENT(x_axis_2d).OR.PRESENT(y_axis_2d) )THEN
    IF ( (m_x <= 0).OR.(m_y <= 0) ) THEN
      CALL ipslerr (3,'fliogstc', &
 &      'Requested [x/y]_axis_2d', &
 &      'but the coordinates are not iboth present','in the file')
    ELSE IF ( ((i_x+l_x-1) > m_x).OR.((i_y+l_y-1) > m_y) ) THEN
      CALL ipslerr (3,'fliogstc', &
 &      'The requested size for the [x/y]_axis_2d', &
 &      'is greater than the size of the coordinate','in the file')
    ENDIF
  ENDIF
!-
! Ensuring data mode
!-
  CALL flio_hdm (f_i,f_e,.FALSE.)
!-
! Extracting the x coordinate, if needed
!-
  IF (PRESENT(x_axis).OR.PRESENT(x_axis_2d)) THEN
    CALL flio_qax (f_i,'x',i_v,nbdim)
    IF (i_v > 0) THEN
      IF      (nbdim == 1) THEN
        IF (PRESENT(x_axis)) THEN
          i_rc = NF90_GET_VAR(f_e,i_v,x_axis, &
 &                 start=(/i_x/),count=(/l_x/))
        ELSE
          ALLOCATE(v_tmp(l_x))
          i_rc = NF90_GET_VAR(f_e,i_v,v_tmp, &
 &                 start=(/i_x/),count=(/l_x/))
          DO kv=1,l_y
            x_axis_2d(:,kv) = v_tmp(:)
          ENDDO
          DEALLOCATE(v_tmp)
        ENDIF
      ELSE IF (nbdim == 2) THEN
        IF (PRESENT(x_axis)) THEN
          l_ok = .TRUE.
          IF (l_y > 1) THEN
            ALLOCATE(v_tmp(l_y))
            DO kv=i_x,i_x+l_x-1
              i_rc = NF90_GET_VAR(f_e,i_v,v_tmp, &
 &                     start=(/kv,i_y/),count=(/1,l_y/))
              IF (ANY(v_tmp(2:l_y) /= v_tmp(1))) THEN
                l_ok = .FALSE.
                EXIT
              ENDIF
            ENDDO
            DEALLOCATE(v_tmp)
          ENDIF
          IF (l_ok) THEN
            i_rc = NF90_GET_VAR(f_e,i_v,x_axis, &
 &                   start=(/i_x,i_y/),count=(/l_x,1/))
          ELSE
            CALL ipslerr (3,'fliogstc', &
 &            'Requested 1D x_axis', &
 &            'which have 2 not regular dimensions', &
 &            'in the file')
          ENDIF
        ELSE
          i_rc = NF90_GET_VAR(f_e,i_v,x_axis_2d, &
 &                 start=(/i_x,i_y/),count=(/l_x,l_y/))
        ENDIF
      ELSE
        CALL ipslerr (3,'fliogstc', &
 &        'Can not handle x_axis', &
 &        'that have more than 2 dimensions', &
 &        'in the file')
      ENDIF
    ELSE
      CALL ipslerr (3,'fliogstc','No x_axis found','in the file',' ')
    ENDIF
  ENDIF
!-
! Extracting the y coordinate, if needed
!-
  IF (PRESENT(y_axis).OR.PRESENT(y_axis_2d)) THEN
    CALL flio_qax (f_i,'y',i_v,nbdim)
    IF (i_v > 0) THEN
      IF      (nbdim == 1) THEN
        IF (PRESENT(y_axis)) THEN
          i_rc = NF90_GET_VAR(f_e,i_v,y_axis, &
 &                 start=(/i_y/),count=(/l_y/))
        ELSE
          ALLOCATE(v_tmp(l_y))
          i_rc = NF90_GET_VAR(f_e,i_v,v_tmp, &
 &                 start=(/i_y/),count=(/l_y/))
          DO kv=1,l_x
            y_axis_2d(kv,:) = v_tmp(:)
          ENDDO
          DEALLOCATE(v_tmp)
        ENDIF
      ELSE IF (nbdim == 2) THEN
        IF (PRESENT(y_axis)) THEN
          l_ok = .TRUE.
          IF (l_x > 1) THEN
            ALLOCATE(v_tmp(l_x))
            DO kv=i_y,i_y+l_y-1
              i_rc = NF90_GET_VAR(f_e,i_v,v_tmp, &
 &                     start=(/i_x,kv/),count=(/l_x,1/))
              IF (ANY(v_tmp(2:l_x) /= v_tmp(1))) THEN
                l_ok = .FALSE.
                EXIT
              ENDIF
            ENDDO
            DEALLOCATE(v_tmp)
          ENDIF
          IF (l_ok) THEN
            i_rc = NF90_GET_VAR(f_e,i_v,y_axis, &
 &                   start=(/i_x,i_y/),count=(/1,l_y/))
          ELSE
            CALL ipslerr (3,'fliogstc', &
 &            'Requested 1D y_axis', &
 &            'which have 2 not regular dimensions', &
 &            'in the file')
          ENDIF
        ELSE
          i_rc = NF90_GET_VAR(f_e,i_v,y_axis_2d, &
 &                 start=(/i_x,i_y/),count=(/l_x,l_y/))
        ENDIF
      ELSE
        CALL ipslerr (3,'fliogstc', &
 &        'Can not handle y axis', &
 &        'that have more than 2 dimensions', &
 &        'in the file')
      ENDIF
    ELSE
      CALL ipslerr (3,'fliogstc','No y_axis found','in the file',' ')
    ENDIF
  ENDIF
!-
! Extracting the z coordinate, if needed
!-
  IF (PRESENT(z_axis)) THEN
    IF (PRESENT(z_start)) THEN
      i_z = z_start
    ELSE
      i_z = 1
    ENDIF
    IF (PRESENT(z_count)) THEN
      l_z = z_count
    ELSE
      l_z = m_z-i_z+1
    ENDIF
    IF ((i_z+l_z-1) > m_z) THEN
      CALL ipslerr (3,'fliogstc', &
 &      'The requested size for the z axis', &
 &      'is greater than the size of the coordinate',&
 &      'in the file')
    ENDIF
    CALL flio_qax (f_i,'z',i_v,nbdim)
    IF (i_v > 0) THEN
      IF (nbdim == 1) THEN
        i_rc = NF90_GET_VAR(f_e,i_v,z_axis, &
 &               start=(/i_z/),count=(/l_z/))
      ELSE
        CALL ipslerr (3,'fliogstc', &
 &        'Can not handle z_axis', &
 &        'that have more than 1 dimension', &
 &        'in the file')
      ENDIF
    ELSE
      CALL ipslerr (3,'fliogstc','No z_axis found','in the file',' ')
    ENDIF
  ENDIF
!-
! Extracting the t coordinate, if needed
!-
  IF (PRESENT(t_axis).OR.PRESENT(t_init).OR.PRESENT(t_step)) THEN
    CALL flio_qax (f_i,'t',i_v,nbdim)
    IF (i_v < 0) THEN
      CALL ipslerr (3,'fliogstc','No t_axis found','in the file',' ')
    ENDIF
!---
    IF (l_dbg) THEN
      WRITE(*,*) '  fliogstc - get time details'
    ENDIF
!---
!-- Get all the details for the time
!-- Prefered method is '"time_steps" since'
!---
    name=''
    i_rc = NF90_INQUIRE_VARIABLE(f_e,i_v,name=name)
    units=''
    i_rc = NF90_GET_ATT(f_e,i_v,'units',units)
    IF      (INDEX(units,' since ') > 0) THEN
      it_t = 1
    ELSE IF (INDEX(name,'tstep') > 0) THEN
      it_t = 2
    ELSE
      it_t = 0;
    ENDIF
  ENDIF
!-
! Extracting the t coordinate, if needed
!-
  IF (PRESENT(t_axis)) THEN
    IF (PRESENT(t_start)) THEN
      i_t = t_start
    ELSE
      i_t = 1
    ENDIF
    IF (PRESENT(t_count)) THEN
      l_t = t_count
    ELSE
      l_t = m_t-i_t+1
    ENDIF
    IF ((i_t+l_t-1) > m_t) THEN
      CALL ipslerr (3,'fliogstc', &
 &      'The requested size for the t axis', &
 &      'is greater than the size of the coordinate',&
 &      'in the file')
    ENDIF
    ALLOCATE(v_tmp(l_t))
    i_rc = NF90_GET_VAR(f_e,i_v,v_tmp, &
 &           start=(/i_t/),count=(/l_t/))
    t_axis(1:l_t) = NINT(v_tmp(1:l_t))
    DEALLOCATE(v_tmp)
!---
    IF (l_dbg) THEN
      WRITE(*,*) '  fliogstc - first time : ',t_axis(1:1)
    ENDIF
  ENDIF
!-
! Extracting the time at the beginning, if needed
!-
  IF (PRESENT(t_init)) THEN
!-- Find the calendar
    CALL lock_calendar (old_status=l_tmp)
    CALL ioget_calendar (c_tmp)
    units = ''
    i_rc = NF90_GET_ATT(f_e,i_v,'calendar',units)
    IF (i_rc == NF90_NOERR) THEN
      CALL lock_calendar (new_status=.FALSE.)
      CALL ioconf_calendar (TRIM(units))
    ENDIF
    IF (it_t == 1) THEN
      units = ''
      i_rc = NF90_GET_ATT(f_e,i_v,'units',units)
      units = units(INDEX(units,' since ')+7:LEN_TRIM(units))
      READ (units,'(I4.4,5(A,I2.2))') &
 &      j_yy,c_1,j_mo,c_1,j_dd,c_1,j_hh,c_1,j_mn,c_1,j_ss
      r_ss = j_hh*3600.+j_mn*60.+j_ss
      CALL ymds2ju (j_yy,j_mo,j_dd,r_ss,t_init)
    ELSE IF (it_t == 2) THEN
      i_rc = NF90_GET_ATT(f_e,NF90_GLOBAL,'year0',r_yy)
      i_rc = NF90_GET_ATT(f_e,NF90_GLOBAL,'month0',r_mo)
      i_rc = NF90_GET_ATT(f_e,NF90_GLOBAL,'day0',r_dd)
      i_rc = NF90_GET_ATT(f_e,NF90_GLOBAL,'sec0',r_ss)
      j_yy = NINT(r_yy); j_mo = NINT(r_mo); j_dd = NINT(r_dd);
      CALL ymds2ju (j_yy,j_mo,j_dd,r_ss,t_init)
    ELSE
      t_init = 0.
    ENDIF
    CALL lock_calendar (new_status=.FALSE.)
    CALL ioconf_calendar (TRIM(c_tmp))
    CALL lock_calendar (new_status=l_tmp)
    IF (l_dbg) THEN
      WRITE(*,*) '  fliogstc - time_type : '
      WRITE(*,*) it_t
      WRITE(*,*) '  fliogstc - year month day second t_init : '
      WRITE(*,*) j_yy,j_mo,j_dd,r_ss,t_init
    ENDIF
  ENDIF
!-
! Extracting the timestep in seconds, if needed
!-
  IF (PRESENT(t_step)) THEN
    IF      (it_t == 1) THEN
      units = ''
      i_rc = NF90_GET_ATT(f_e,i_v,'units',units)
      units = ADJUSTL(units(1:INDEX(units,' since ')-1))
      dtn = 1.
      IF      (INDEX(units,"week") /= 0) THEN
        kv  = INDEX(units,"week")
        dtv = 604800.
      ELSE IF (INDEX(units,"day")  /= 0) THEN
        kv  = INDEX(units,"day")
        dtv = 86400.
      ELSE IF (INDEX(units,"h")    /= 0) THEN
        kv  = INDEX(units,"h")
        dtv = 3600.
      ELSE IF (INDEX(units,"min")  /= 0) THEN
        kv  = INDEX(units,"min")
        dtv = 60.
      ELSE IF (INDEX(units,"sec")  /= 0) THEN
        kv  = INDEX(units,"sec")
        dtv = 1.
      ELSE IF (INDEX(units,"timesteps") /= 0) THEN
        kv  = INDEX(units,"timesteps")
        i_rc = NF90_GET_ATT(f_e,i_v,'tstep_sec',dtv)
        IF (i_rc /= NF90_NOERR) THEN
          CALL ipslerr (3,'fliogstc','"timesteps" value', &
 &                        'not found','in the file')
        ENDIF
      ELSE
        kv  = 1
        dtv = 1.
      ENDIF
      IF (kv > 1) THEN
        READ (unit=units(1:kv-1),FMT=*) dtn
      ENDIF
      t_step = dtn*dtv
    ELSE IF (it_t == 2) THEN
      i_rc = NF90_GET_ATT(f_e,NF90_GLOBAL,'delta_tstep_sec',t_step)
    ELSE
      t_step = 1.
    ENDIF
  ENDIF
!-
! Extracting the calendar attribute, if needed
!-
  IF (PRESENT(t_calendar)) THEN
    units = ''
    i_rc = NF90_GET_ATT(f_e,i_v,'calendar',units)
    IF (i_rc == NF90_NOERR) THEN
      t_calendar = units
    ELSE
      t_calendar = "not found"
    ENDIF
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-fliogstc"
  ENDIF
!----------------------
END SUBROUTINE fliogstc
!===
SUBROUTINE flioinqv &
 & (f_i,v_n,l_ex,v_t,nb_dims,len_dims,id_dims, &
 &  nb_atts,cn_atts,ia_start,ia_count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  LOGICAL,INTENT(OUT) :: l_ex
  INTEGER,OPTIONAL,INTENT(OUT) :: v_t,nb_dims,nb_atts
  INTEGER,OPTIONAL,INTENT(OUT),DIMENSION(:) :: len_dims,id_dims
  CHARACTER(LEN=*),DIMENSION(:),OPTIONAL,INTENT(OUT) :: cn_atts
  INTEGER,OPTIONAL,INTENT(IN) :: ia_start,ia_count
!-
  INTEGER :: f_e,i_v,n_w,i_s,i_w,iws,iwc,i_rc
  LOGICAL :: l_ok
  INTEGER,DIMENSION(NF90_MAX_VAR_DIMS) :: dim_ids
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->flioinqv ",TRIM(v_n)
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('flioinqv',f_i,f_e)
!-
  i_v = -1
  i_rc = NF90_INQ_VARID(f_e,v_n,i_v)
!-
  l_ex = ( (i_v >= 0).AND.(i_rc == NF90_NOERR) )
!-
  IF (l_ex) THEN
    IF (PRESENT(v_t)) THEN
      i_rc = NF90_INQUIRE_VARIABLE(f_e,i_v,xtype=v_t)
    ENDIF
    n_w = -1
    IF (PRESENT(nb_dims).OR.PRESENT(len_dims)) THEN
      i_rc = NF90_INQUIRE_VARIABLE(f_e,i_v, &
 &             ndims=n_w,dimids=dim_ids)
      IF (PRESENT(nb_dims)) THEN
        nb_dims = n_w
      ENDIF
      IF (PRESENT(len_dims)) THEN
        i_s = SIZE(len_dims)
        len_dims(:) = -1
        IF (i_s < n_w) THEN
          CALL ipslerr (2,'flioinqv', &
 &         'Only the first dimensions of the variable', &
 &         TRIM(v_n),'will be returned')
        ENDIF
        DO i_w=1,MIN(n_w,i_s)
          i_rc = NF90_INQUIRE_DIMENSION(f_e,dim_ids(i_w), &
 &                                      len=len_dims(i_w))
        ENDDO
      ENDIF
      IF (PRESENT(id_dims)) THEN
        i_s = SIZE(id_dims)
        id_dims(:) = -1
        IF (i_s < n_w) THEN
          CALL ipslerr (2,'flioinqv', &
 &         'The number of dimensions to retrieve', &
 &         'is greater than the size of the array,', &
 &         'only the first dimensions of "' &
 &           //TRIM(v_n)//'" will be returned')
        ENDIF
        i_w = MIN(n_w,i_s)
        id_dims(1:i_w) = dim_ids(1:i_w)
      ENDIF
    ENDIF
    IF (PRESENT(nb_atts).OR.PRESENT(cn_atts)) THEN
      i_rc = NF90_INQUIRE_VARIABLE(f_e,i_v,nAtts=n_w)
      IF (PRESENT(nb_atts)) THEN
        nb_atts = n_w
      ENDIF
      IF (PRESENT(cn_atts)) THEN
        l_ok = .TRUE.
        i_s = SIZE(cn_atts)
        DO i_w=1,i_s
          cn_atts(i_w)(:) = '?'
        ENDDO
        IF (PRESENT(ia_start)) THEN
          iws = ia_start
        ELSE
          iws = 1
        ENDIF
        IF (PRESENT(ia_count)) THEN
          iwc = ia_count
        ELSE
          iwc = n_w
        ENDIF
        IF (iws > n_w) THEN
          l_ok = .FALSE.
          CALL ipslerr (2,'flioinqv', &
 &         'The start index of requested attributes', &
 &         'is greater than the number of attributes of', &
 &         '"'//TRIM(v_n)//'"')
        ELSE IF (iws < 1) THEN
          l_ok = .FALSE.
          CALL ipslerr (2,'flioinqv', &
 &         'The start index of requested attributes', &
 &         'is invalid ( < 1 ) for', &
 &         '"'//TRIM(v_n)//'"')
        ENDIF
        IF ((iws+iwc-1) > n_w) THEN
          CALL ipslerr (2,'flioinqv', &
 &         'The number of requested attributes', &
 &         'is greater than the number of attributes of', &
 &         '"'//TRIM(v_n)//'"')
        ENDIF
        IF (iwc > i_s) THEN
          CALL ipslerr (2,'flioinqv', &
 &         'The number of attributes to retrieve', &
 &         'is greater than the size of the array,', &
 &         'only the first attributes of "' &
 &           //TRIM(v_n)//'" will be returned')
        ELSE IF (iwc < 1) THEN
          l_ok = .FALSE.
          CALL ipslerr (2,'flioinqv', &
 &         'The number of requested attributes', &
 &         'is invalid ( < 1 ) for', &
 &         '"'//TRIM(v_n)//'"')
        ENDIF
        IF (l_ok) THEN
          DO i_w=1,MIN(iwc,i_s,n_w-iws+1)
            i_rc = NF90_INQ_ATTNAME(f_e, &
 &                  i_v,i_w+iws-1,name=cn_atts(i_w))
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-flioinqv"
  ENDIF
!----------------------
END SUBROUTINE flioinqv
!===
SUBROUTINE fliogv_i40 (f_i,v_n,v_v,start)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,i_40=v_v,start=start)
!------------------------
END SUBROUTINE fliogv_i40
!===
SUBROUTINE fliogv_i41 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),DIMENSION(:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,i_41=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_i41
!===
SUBROUTINE fliogv_i42 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),DIMENSION(:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,i_42=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_i42
!===
SUBROUTINE fliogv_i43 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),DIMENSION(:,:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,i_43=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_i43
!===
SUBROUTINE fliogv_i44 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),DIMENSION(:,:,:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,i_44=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_i44
!===
SUBROUTINE fliogv_i45 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),DIMENSION(:,:,:,:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,i_45=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_i45
!===
SUBROUTINE fliogv_i20 (f_i,v_n,v_v,start)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_2),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,i_20=v_v,start=start)
!------------------------
END SUBROUTINE fliogv_i20
!===
SUBROUTINE fliogv_i21 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_2),DIMENSION(:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,i_21=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_i21
!===
SUBROUTINE fliogv_i22 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_2),DIMENSION(:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,i_22=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_i22
!===
SUBROUTINE fliogv_i23 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_2),DIMENSION(:,:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,i_23=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_i23
!===
SUBROUTINE fliogv_i24 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_2),DIMENSION(:,:,:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,i_24=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_i24
!===
SUBROUTINE fliogv_i25 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_2),DIMENSION(:,:,:,:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,i_25=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_i25
!===
!?INTEGERS of KIND 1 are not supported on all computers
!?SUBROUTINE fliogv_i10 (f_i,v_n,v_v,start)
!?!---------------------------------------------------------------------
!?  IMPLICIT NONE
!?!-
!?  INTEGER,INTENT(IN) :: f_i
!?  CHARACTER(LEN=*),INTENT(IN) :: v_n
!?  INTEGER(KIND=i_1),INTENT(OUT) :: v_v
!?  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start
!?!---------------------------------------------------------------------
!?  CALL flio_ugv (f_i,v_n,i_10=v_v,start=start)
!?!------------------------
!?END SUBROUTINE fliogv_i10
!?!===
!?SUBROUTINE fliogv_i11 (f_i,v_n,v_v,start,count)
!?!---------------------------------------------------------------------
!?  IMPLICIT NONE
!?!-
!?  INTEGER,INTENT(IN) :: f_i
!?  CHARACTER(LEN=*),INTENT(IN) :: v_n
!?  INTEGER(KIND=i_1),DIMENSION(:),INTENT(OUT) :: v_v
!?  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!?!---------------------------------------------------------------------
!?  CALL flio_ugv (f_i,v_n,i_11=v_v,start=start,count=count)
!?!------------------------
!?END SUBROUTINE fliogv_i11
!?!===
!?SUBROUTINE fliogv_i12 (f_i,v_n,v_v,start,count)
!?!---------------------------------------------------------------------
!?  IMPLICIT NONE
!?!-
!?  INTEGER,INTENT(IN) :: f_i
!?  CHARACTER(LEN=*),INTENT(IN) :: v_n
!?  INTEGER(KIND=i_1),DIMENSION(:,:),INTENT(OUT) :: v_v
!?  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!?!---------------------------------------------------------------------
!?  CALL flio_ugv (f_i,v_n,i_12=v_v,start=start,count=count)
!?!------------------------
!?END SUBROUTINE fliogv_i12
!?!===
!?SUBROUTINE fliogv_i13 (f_i,v_n,v_v,start,count)
!?!---------------------------------------------------------------------
!?  IMPLICIT NONE
!?!-
!?  INTEGER,INTENT(IN) :: f_i
!?  CHARACTER(LEN=*),INTENT(IN) :: v_n
!?  INTEGER(KIND=i_1),DIMENSION(:,:,:),INTENT(OUT) :: v_v
!?  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!?!---------------------------------------------------------------------
!?  CALL flio_ugv (f_i,v_n,i_13=v_v,start=start,count=count)
!?!------------------------
!?END SUBROUTINE fliogv_i13
!?!===
!?SUBROUTINE fliogv_i14 (f_i,v_n,v_v,start,count)
!?!---------------------------------------------------------------------
!?  IMPLICIT NONE
!?!-
!?  INTEGER,INTENT(IN) :: f_i
!?  CHARACTER(LEN=*),INTENT(IN) :: v_n
!?  INTEGER(KIND=i_1),DIMENSION(:,:,:,:),INTENT(OUT) :: v_v
!?  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!?!---------------------------------------------------------------------
!?  CALL flio_ugv (f_i,v_n,i_14=v_v,start=start,count=count)
!?!------------------------
!?END SUBROUTINE fliogv_i14
!?!===
!?SUBROUTINE fliogv_i15 (f_i,v_n,v_v,start,count)
!?!---------------------------------------------------------------------
!?  IMPLICIT NONE
!?!-
!?  INTEGER,INTENT(IN) :: f_i
!?  CHARACTER(LEN=*),INTENT(IN) :: v_n
!?  INTEGER(KIND=i_1),DIMENSION(:,:,:,:,:),INTENT(OUT) :: v_v
!?  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!?!---------------------------------------------------------------------
!?  CALL flio_ugv (f_i,v_n,i_15=v_v,start=start,count=count)
!?!------------------------
!?END SUBROUTINE fliogv_i15
!===
SUBROUTINE fliogv_r40 (f_i,v_n,v_v,start)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_4),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,r_40=v_v,start=start)
!------------------------
END SUBROUTINE fliogv_r40
!===
SUBROUTINE fliogv_r41 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_4),DIMENSION(:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,r_41=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_r41
!===
SUBROUTINE fliogv_r42 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_4),DIMENSION(:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,r_42=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_r42
!===
SUBROUTINE fliogv_r43 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_4),DIMENSION(:,:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,r_43=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_r43
!===
SUBROUTINE fliogv_r44 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_4),DIMENSION(:,:,:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,r_44=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_r44
!===
SUBROUTINE fliogv_r45 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_4),DIMENSION(:,:,:,:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,r_45=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_r45
!===
SUBROUTINE fliogv_r80 (f_i,v_n,v_v,start)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_8),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,r_80=v_v,start=start)
!------------------------
END SUBROUTINE fliogv_r80
!===
SUBROUTINE fliogv_r81 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_8),DIMENSION(:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,r_81=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_r81
!===
SUBROUTINE fliogv_r82 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_8),DIMENSION(:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,r_82=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_r82
!===
SUBROUTINE fliogv_r83 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_8),DIMENSION(:,:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,r_83=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_r83
!===
SUBROUTINE fliogv_r84 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_8),DIMENSION(:,:,:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,r_84=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_r84
!===
SUBROUTINE fliogv_r85 (f_i,v_n,v_v,start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  REAL(KIND=r_8),DIMENSION(:,:,:,:,:),INTENT(OUT) :: v_v
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(IN) :: start,count
!---------------------------------------------------------------------
  CALL flio_ugv (f_i,v_n,r_85=v_v,start=start,count=count)
!------------------------
END SUBROUTINE fliogv_r85
!===
SUBROUTINE flio_ugv &
 & (f_i,v_n, &
 &  i_40,i_41,i_42,i_43,i_44,i_45, &
 &  i_20,i_21,i_22,i_23,i_24,i_25, &
!? &  i_10,i_11,i_12,i_13,i_14,i_15, &
 &  r_40,r_41,r_42,r_43,r_44,r_45, &
 &  r_80,r_81,r_82,r_83,r_84,r_85, &
 &  start,count)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n
  INTEGER(KIND=i_4),INTENT(OUT),OPTIONAL :: i_40
  INTEGER(KIND=i_4),DIMENSION(:),INTENT(OUT),OPTIONAL :: i_41
  INTEGER(KIND=i_4),DIMENSION(:,:),INTENT(OUT),OPTIONAL :: i_42
  INTEGER(KIND=i_4),DIMENSION(:,:,:),INTENT(OUT),OPTIONAL :: i_43
  INTEGER(KIND=i_4),DIMENSION(:,:,:,:),INTENT(OUT),OPTIONAL :: i_44
  INTEGER(KIND=i_4),DIMENSION(:,:,:,:,:),INTENT(OUT),OPTIONAL :: i_45
  INTEGER(KIND=i_2),INTENT(OUT),OPTIONAL :: i_20
  INTEGER(KIND=i_2),DIMENSION(:),INTENT(OUT),OPTIONAL :: i_21
  INTEGER(KIND=i_2),DIMENSION(:,:),INTENT(OUT),OPTIONAL :: i_22
  INTEGER(KIND=i_2),DIMENSION(:,:,:),INTENT(OUT),OPTIONAL :: i_23
  INTEGER(KIND=i_2),DIMENSION(:,:,:,:),INTENT(OUT),OPTIONAL :: i_24
  INTEGER(KIND=i_2),DIMENSION(:,:,:,:,:),INTENT(OUT),OPTIONAL :: i_25
!?INTEGERS of KIND 1 are not supported on all computers
!?INTEGER(KIND=i_1),INTENT(OUT),OPTIONAL :: i_10
!?INTEGER(KIND=i_1),DIMENSION(:),INTENT(OUT),OPTIONAL :: i_11
!?INTEGER(KIND=i_1),DIMENSION(:,:),INTENT(OUT),OPTIONAL :: i_12
!?INTEGER(KIND=i_1),DIMENSION(:,:,:),INTENT(OUT),OPTIONAL :: i_13
!?INTEGER(KIND=i_1),DIMENSION(:,:,:,:),INTENT(OUT),OPTIONAL :: i_14
!?INTEGER(KIND=i_1),DIMENSION(:,:,:,:,:),INTENT(OUT),OPTIONAL :: i_15
  REAL(KIND=r_4),INTENT(OUT),OPTIONAL :: r_40
  REAL(KIND=r_4),DIMENSION(:),INTENT(OUT),OPTIONAL :: r_41
  REAL(KIND=r_4),DIMENSION(:,:),INTENT(OUT),OPTIONAL :: r_42
  REAL(KIND=r_4),DIMENSION(:,:,:),INTENT(OUT),OPTIONAL :: r_43
  REAL(KIND=r_4),DIMENSION(:,:,:,:),INTENT(OUT),OPTIONAL :: r_44
  REAL(KIND=r_4),DIMENSION(:,:,:,:,:),INTENT(OUT),OPTIONAL :: r_45
  REAL(KIND=r_8),INTENT(OUT),OPTIONAL :: r_80
  REAL(KIND=r_8),DIMENSION(:),INTENT(OUT),OPTIONAL :: r_81
  REAL(KIND=r_8),DIMENSION(:,:),INTENT(OUT),OPTIONAL :: r_82
  REAL(KIND=r_8),DIMENSION(:,:,:),INTENT(OUT),OPTIONAL :: r_83
  REAL(KIND=r_8),DIMENSION(:,:,:,:),INTENT(OUT),OPTIONAL :: r_84
  REAL(KIND=r_8),DIMENSION(:,:,:,:,:),INTENT(OUT),OPTIONAL :: r_85
  INTEGER,DIMENSION(:),INTENT(IN),OPTIONAL :: start,count
!-
  INTEGER :: f_e,i_v,i_rc
  CHARACTER(LEN=5) :: cvr_d
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    IF      (PRESENT(i_40)) THEN; cvr_d = "I1 0D";
    ELSE IF (PRESENT(i_41)) THEN; cvr_d = "I1 1D";
    ELSE IF (PRESENT(i_42)) THEN; cvr_d = "I1 2D";
    ELSE IF (PRESENT(i_43)) THEN; cvr_d = "I1 3D";
    ELSE IF (PRESENT(i_44)) THEN; cvr_d = "I1 4D";
    ELSE IF (PRESENT(i_45)) THEN; cvr_d = "I1 5D";
    ELSE IF (PRESENT(i_20)) THEN; cvr_d = "I2 0D";
    ELSE IF (PRESENT(i_21)) THEN; cvr_d = "I2 1D";
    ELSE IF (PRESENT(i_22)) THEN; cvr_d = "I2 2D";
    ELSE IF (PRESENT(i_23)) THEN; cvr_d = "I2 3D";
    ELSE IF (PRESENT(i_24)) THEN; cvr_d = "I2 4D";
    ELSE IF (PRESENT(i_25)) THEN; cvr_d = "I2 5D";
!?  ELSE IF (PRESENT(i_10)) THEN; cvr_d = "I4 0D";
!?  ELSE IF (PRESENT(i_11)) THEN; cvr_d = "I4 1D";
!?  ELSE IF (PRESENT(i_12)) THEN; cvr_d = "I4 2D";
!?  ELSE IF (PRESENT(i_13)) THEN; cvr_d = "I4 3D";
!?  ELSE IF (PRESENT(i_14)) THEN; cvr_d = "I4 4D";
!?  ELSE IF (PRESENT(i_15)) THEN; cvr_d = "I4 5D";
    ELSE IF (PRESENT(r_40)) THEN; cvr_d = "R4 0D";
    ELSE IF (PRESENT(r_41)) THEN; cvr_d = "R4 1D";
    ELSE IF (PRESENT(r_42)) THEN; cvr_d = "R4 2D";
    ELSE IF (PRESENT(r_43)) THEN; cvr_d = "R4 3D";
    ELSE IF (PRESENT(r_44)) THEN; cvr_d = "R4 4D";
    ELSE IF (PRESENT(r_45)) THEN; cvr_d = "R4 5D";
    ELSE IF (PRESENT(r_80)) THEN; cvr_d = "R8 0D";
    ELSE IF (PRESENT(r_81)) THEN; cvr_d = "R8 1D";
    ELSE IF (PRESENT(r_82)) THEN; cvr_d = "R8 2D";
    ELSE IF (PRESENT(r_83)) THEN; cvr_d = "R8 3D";
    ELSE IF (PRESENT(r_84)) THEN; cvr_d = "R8 4D";
    ELSE IF (PRESENT(r_85)) THEN; cvr_d = "R8 5D";
    ENDIF
    WRITE(*,*) "->fliogetv ",TRIM(v_n)," ",TRIM(cvr_d)
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('fliogetv',f_i,f_e)
!-
! Ensuring data mode
!-
  CALL flio_hdm (f_i,f_e,.FALSE.)
!-
  i_rc = NF90_INQ_VARID(f_e,v_n,i_v)
  IF (i_rc == NF90_NOERR) THEN
    IF      (PRESENT(i_40)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,i_40,start=start)
    ELSE IF (PRESENT(i_41)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,i_41,start=start,count=count)
    ELSE IF (PRESENT(i_42)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,i_42,start=start,count=count)
    ELSE IF (PRESENT(i_43)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,i_43,start=start,count=count)
    ELSE IF (PRESENT(i_44)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,i_44,start=start,count=count)
    ELSE IF (PRESENT(i_45)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,i_45,start=start,count=count)
    ELSE IF (PRESENT(i_20)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,i_20,start=start)
    ELSE IF (PRESENT(i_21)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,i_21,start=start,count=count)
    ELSE IF (PRESENT(i_22)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,i_22,start=start,count=count)
    ELSE IF (PRESENT(i_23)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,i_23,start=start,count=count)
    ELSE IF (PRESENT(i_24)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,i_24,start=start,count=count)
    ELSE IF (PRESENT(i_25)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,i_25,start=start,count=count)
!?  ELSE IF (PRESENT(i_10)) THEN
!?    i_rc = NF90_GET_VAR(f_e,i_v,i_10,start=start)
!?  ELSE IF (PRESENT(i_11)) THEN
!?    i_rc = NF90_GET_VAR(f_e,i_v,i_11,start=start,count=count)
!?  ELSE IF (PRESENT(i_12)) THEN
!?    i_rc = NF90_GET_VAR(f_e,i_v,i_12,start=start,count=count)
!?  ELSE IF (PRESENT(i_13)) THEN
!?    i_rc = NF90_GET_VAR(f_e,i_v,i_13,start=start,count=count)
!?  ELSE IF (PRESENT(i_14)) THEN
!?    i_rc = NF90_GET_VAR(f_e,i_v,i_14,start=start,count=count)
!?  ELSE IF (PRESENT(i_15)) THEN
!?    i_rc = NF90_GET_VAR(f_e,i_v,i_15,start=start,count=count)
    ELSE IF (PRESENT(r_40)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,r_40,start=start)
    ELSE IF (PRESENT(r_41)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,r_41,start=start,count=count)
    ELSE IF (PRESENT(r_42)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,r_42,start=start,count=count)
    ELSE IF (PRESENT(r_43)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,r_43,start=start,count=count)
    ELSE IF (PRESENT(r_44)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,r_44,start=start,count=count)
    ELSE IF (PRESENT(r_45)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,r_45,start=start,count=count)
    ELSE IF (PRESENT(r_80)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,r_80,start=start)
    ELSE IF (PRESENT(r_81)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,r_81,start=start,count=count)
    ELSE IF (PRESENT(r_82)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,r_82,start=start,count=count)
    ELSE IF (PRESENT(r_83)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,r_83,start=start,count=count)
    ELSE IF (PRESENT(r_84)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,r_84,start=start,count=count)
    ELSE IF (PRESENT(r_85)) THEN
      i_rc = NF90_GET_VAR(f_e,i_v,r_85,start=start,count=count)
    ENDIF
    IF (i_rc /= NF90_NOERR) THEN
      CALL ipslerr (3,'fliogetv', &
 &      'Variable '//TRIM(v_n)//' not get','Error :', &
 &      TRIM(NF90_STRERROR(i_rc)))
    ENDIF
  ELSE
    CALL ipslerr (3,'fliogetv','Variable',TRIM(v_n),'not found')
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-fliogetv"
  ENDIF
!----------------------
END SUBROUTINE flio_ugv
!===
SUBROUTINE flioinqa (f_i,v_n,a_n,l_ex,a_t,a_l)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  LOGICAL,INTENT(OUT) :: l_ex
  INTEGER,OPTIONAL,INTENT(OUT) :: a_t,a_l
!-
  INTEGER :: i_rc,f_e,i_v,t_ea,l_ea
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->flioinqa ",TRIM(v_n),"-",TRIM(a_n)
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('flioinqa',f_i,f_e)
!-
  IF (TRIM(v_n) == '?') THEN
    i_v = NF90_GLOBAL
  ELSE
    i_rc = NF90_INQ_VARID(f_e,v_n,i_v)
    IF (i_rc /= NF90_NOERR) THEN
      CALL ipslerr (3,'flioinqa', &
       'Variable :',TRIM(v_n),'not found')
    ENDIF
  ENDIF
!-
  i_rc = NF90_INQUIRE_ATTRIBUTE(f_e,i_v,a_n,xtype=t_ea,len=l_ea)
!-
  l_ex = (i_rc == NF90_NOERR)
!-
  IF (l_ex) THEN
    IF (PRESENT(a_t)) THEN
      a_t = t_ea
    ENDIF
    IF (PRESENT(a_l)) THEN
      a_l = l_ea
    ENDIF
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-flioinqa"
  ENDIF
!----------------------
END SUBROUTINE flioinqa
!===
SUBROUTINE flioga_r4_0d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  REAL(KIND=4),INTENT(OUT) :: a_v
!---------------------------------------------------------------------
  CALL flio_uga (f_i,v_n,a_n,avr_4_0=a_v)
!---------------------------
END SUBROUTINE flioga_r4_0d
!===
SUBROUTINE flioga_r4_1d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  REAL(KIND=4),DIMENSION(:),INTENT(OUT) :: a_v
!---------------------------------------------------------------------
  CALL flio_uga (f_i,v_n,a_n,avr_4_1=a_v)
!--------------------------
END SUBROUTINE flioga_r4_1d
!===
SUBROUTINE flioga_r8_0d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  REAL(KIND=8),INTENT(OUT) :: a_v
!---------------------------------------------------------------------
  CALL flio_uga (f_i,v_n,a_n,avr_8_0=a_v)
!---------------------------
END SUBROUTINE flioga_r8_0d
!===
SUBROUTINE flioga_r8_1d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  REAL(KIND=8),DIMENSION(:),INTENT(OUT) :: a_v
!---------------------------------------------------------------------
  CALL flio_uga (f_i,v_n,a_n,avr_8_1=a_v)
!--------------------------
END SUBROUTINE flioga_r8_1d
!===
SUBROUTINE flioga_i4_0d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  INTEGER(KIND=4),INTENT(OUT) :: a_v
!---------------------------------------------------------------------
  CALL flio_uga (f_i,v_n,a_n,avi_4_0=a_v)
!---------------------------
END SUBROUTINE flioga_i4_0d
!===
SUBROUTINE flioga_i4_1d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  INTEGER(KIND=4),DIMENSION(:),INTENT(OUT) :: a_v
!---------------------------------------------------------------------
  CALL flio_uga (f_i,v_n,a_n,avi_4_1=a_v)
!--------------------------
END SUBROUTINE flioga_i4_1d
!===
SUBROUTINE flioga_tx_0d (f_i,v_n,a_n,a_v)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  CHARACTER(LEN=*),INTENT(OUT) :: a_v
!---------------------------------------------------------------------
  CALL flio_uga (f_i,v_n,a_n,avtx=a_v)
!---------------------------
END SUBROUTINE flioga_tx_0d
!===
SUBROUTINE flio_uga &
 & (f_i,v_n,a_n, &
 &  avr_4_0,avr_4_1,avr_8_0,avr_8_1,avi_4_0,avi_4_1,avtx)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
  REAL(KIND=4),OPTIONAL,INTENT(OUT) :: avr_4_0
  REAL(KIND=4),DIMENSION(:),OPTIONAL,INTENT(OUT) :: avr_4_1
  REAL(KIND=8),OPTIONAL,INTENT(OUT) :: avr_8_0
  REAL(KIND=8),DIMENSION(:),OPTIONAL,INTENT(OUT) :: avr_8_1
  INTEGER(KIND=4),OPTIONAL,INTENT(OUT) :: avi_4_0
  INTEGER(KIND=4),DIMENSION(:),OPTIONAL,INTENT(OUT) :: avi_4_1
  CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: avtx
!-
  INTEGER :: f_e,l_ua,i_v,t_ea,l_ea,i_rc
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->fliogeta ",TRIM(v_n)," ",TRIM(a_n)
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('fliogeta',f_i,f_e)
!-
  IF (TRIM(v_n) == '?') THEN
    i_v = NF90_GLOBAL
  ELSE
    i_rc = NF90_INQ_VARID(f_e,v_n,i_v)
    IF (i_rc /= NF90_NOERR) THEN
      CALL ipslerr (3,'fliogeta', &
       'Variable :',TRIM(v_n),'not found')
    ENDIF
  ENDIF
!-
  i_rc = NF90_INQUIRE_ATTRIBUTE(f_e,i_v,a_n,xtype=t_ea,len=l_ea)
  IF (i_rc /= NF90_NOERR) THEN
    CALL ipslerr (3,'fliogeta', &
 &   'Attribute :',TRIM(a_n),'not found')
  ENDIF
!-
  IF ( (.NOT.PRESENT(avtx).AND.(t_ea == NF90_CHAR)) &
 &      .OR.(PRESENT(avtx).AND.(t_ea /= NF90_CHAR)) ) THEN
    CALL ipslerr (3,'fliogeta', &
 &   'The external type of the attribute :',TRIM(a_n), &
 &   'is not compatible with the type of the argument')
  ENDIF
!-
  IF      (PRESENT(avr_4_1)) THEN
    l_ua = SIZE(avr_4_1)
  ELSE IF (PRESENT(avr_8_1)) THEN
    l_ua = SIZE(avr_8_1)
  ELSE IF (PRESENT(avi_4_1)) THEN
    l_ua = SIZE(avi_4_1)
  ELSE IF (PRESENT(avtx)) THEN
    l_ua = LEN(avtx)
  ELSE
    l_ua = 1
  ENDIF
!-
  IF (l_ua < l_ea) THEN
    CALL ipslerr (3,'fliogeta', &
     'Insufficient size of the argument', &
 &   'to receive the values of the attribute :',TRIM(a_n))
  ENDIF
!-
  IF      (PRESENT(avr_4_0)) THEN
    i_rc = NF90_GET_ATT(f_e,i_v,a_n,avr_4_0)
  ELSE IF (PRESENT(avr_4_1)) THEN
    i_rc = NF90_GET_ATT(f_e,i_v,a_n,avr_4_1(1:l_ea))
  ELSE IF (PRESENT(avr_8_0)) THEN
    i_rc = NF90_GET_ATT(f_e,i_v,a_n,avr_8_0)
  ELSE IF (PRESENT(avr_8_1)) THEN
    i_rc = NF90_GET_ATT(f_e,i_v,a_n,avr_8_1(1:l_ea))
  ELSE IF (PRESENT(avi_4_0)) THEN
    i_rc = NF90_GET_ATT(f_e,i_v,a_n,avi_4_0)
  ELSE IF (PRESENT(avi_4_1)) THEN
    i_rc = NF90_GET_ATT(f_e,i_v,a_n,avi_4_1(1:l_ea))
  ELSE IF (PRESENT(avtx)) THEN
    i_rc = NF90_GET_ATT(f_e,i_v,a_n,avtx)
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-fliogeta"
  ENDIF
!----------------------
END SUBROUTINE flio_uga
!===
SUBROUTINE fliorenv (f_i,v_o_n,v_n_n)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_o_n,v_n_n
!-
  INTEGER :: f_e,i_v,i_rc
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) &
 &    "->fliorenv ",TRIM(v_o_n),"->",TRIM(v_n_n)
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('fliorenv',f_i,f_e)
!-
  i_rc = NF90_INQ_VARID(f_e,v_o_n,i_v)
  IF (i_rc /= NF90_NOERR) THEN
    CALL ipslerr (2,'fliorenv', &
     'Variable :',TRIM(v_o_n),'not found')
  ELSE
    CALL flio_hdm (f_i,f_e,.TRUE.)
    i_rc = NF90_RENAME_VAR(f_e,i_v,v_n_n)
    IF (i_rc /= NF90_NOERR) THEN
      CALL ipslerr (2,'fliorenv', &
       'Variable :',TRIM(v_o_n),'can not be renamed')
    ENDIF
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-fliorenv"
  ENDIF
!----------------------
END SUBROUTINE fliorenv
!===
SUBROUTINE fliorena (f_i,v_n,a_o_n,a_n_n)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_o_n,a_n_n
!-
  INTEGER :: f_e,i_v,i_a,i_rc
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) &
 &    "->fliorena ",TRIM(v_n),"-",TRIM(a_o_n),"->",TRIM(a_n_n)
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('fliorena',f_i,f_e)
!-
  IF (TRIM(v_n) == '?') THEN
    i_v = NF90_GLOBAL
  ELSE
    i_rc = NF90_INQ_VARID(f_e,v_n,i_v)
    IF (i_rc /= NF90_NOERR) THEN
      CALL ipslerr (3,'fliorena', &
       'Variable :',TRIM(v_n),'not found')
    ENDIF
  ENDIF
!-
  i_rc = NF90_INQUIRE_ATTRIBUTE(f_e,i_v,a_o_n,attnum=i_a)
  IF (i_rc /= NF90_NOERR) THEN
    CALL ipslerr (2,'fliorena', &
     'Attribute :',TRIM(a_o_n),'not found')
  ELSE
    CALL flio_hdm (f_i,f_e,.TRUE.)
    i_rc = NF90_RENAME_ATT(f_e,i_v,a_o_n,a_n_n)
    IF (i_rc /= NF90_NOERR) THEN
      CALL ipslerr (2,'fliorena', &
       'Attribute :',TRIM(a_o_n),'can not be renamed')
    ENDIF
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-fliorena"
  ENDIF
!----------------------
END SUBROUTINE fliorena
!===
SUBROUTINE fliodela (f_i,v_n,a_n)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: v_n,a_n
!-
  INTEGER :: f_e,i_v,i_a,i_rc
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->fliodela ",TRIM(v_n),"-",TRIM(a_n)
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('fliodela',f_i,f_e)
!-
  IF (TRIM(v_n) == '?') THEN
    i_v = NF90_GLOBAL
  ELSE
    i_rc = NF90_INQ_VARID(f_e,v_n,i_v)
    IF (i_rc /= NF90_NOERR) THEN
      CALL ipslerr (3,'fliodela', &
 &     'Variable :',TRIM(v_n),'not found')
    ENDIF
  ENDIF
!-
  i_rc = NF90_INQUIRE_ATTRIBUTE(f_e,i_v,a_n,attnum=i_a)
  IF (i_rc /= NF90_NOERR) THEN
    CALL ipslerr (2,'fliodela', &
 &   'Attribute :',TRIM(a_n),'not found')
  ELSE
    IF (i_v == NF90_GLOBAL) THEN
      nw_na(f_i) = nw_na(f_i)-1
    ENDIF
    CALL flio_hdm (f_i,f_e,.TRUE.)
    i_rc = NF90_DEL_ATT(f_e,i_v,a_n)
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-fliodela"
  ENDIF
!----------------------
END SUBROUTINE fliodela
!===
SUBROUTINE fliocpya (f_i_i,v_n_i,a_n,f_i_o,v_n_o)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i_i,f_i_o
  CHARACTER(LEN=*),INTENT(IN) :: v_n_i,a_n,v_n_o
!-
  INTEGER :: f_e_i,f_e_o,i_v_i,i_v_o,i_a,i_rc
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->fliocpya - file",f_i_i,"-",TRIM(v_n_i),"-",TRIM(a_n)
    WRITE(*,*) "  copied to file ",f_i_o,"-",TRIM(v_n_o)
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('fliocpya',f_i_i,f_e_i)
  CALL flio_qvid ('fliocpya',f_i_o,f_e_o)
!-
  IF (TRIM(v_n_i) == '?') THEN
    i_v_i = NF90_GLOBAL
  ELSE
    i_rc = NF90_INQ_VARID(f_e_i,v_n_i,i_v_i)
    IF (i_rc /= NF90_NOERR) THEN
      CALL ipslerr (3,'fliocpya', &
 &     'Variable :',TRIM(v_n_i),'not found')
    ENDIF
  ENDIF
!-
  IF (TRIM(v_n_o) == '?') THEN
    i_v_o = NF90_GLOBAL
  ELSE
    i_rc = NF90_INQ_VARID(f_e_o,v_n_o,i_v_o)
    IF (i_rc /= NF90_NOERR) THEN
      CALL ipslerr (3,'fliocpya', &
 &     'Variable :',TRIM(v_n_o),'not found')
    ENDIF
  ENDIF
!-
  i_rc = NF90_INQUIRE_ATTRIBUTE(f_e_i,i_v_i,a_n,attnum=i_a)
  IF (i_rc /= NF90_NOERR) THEN
    CALL ipslerr (3,'fliocpya', &
     'Attribute :',TRIM(a_n),'not found')
  ELSE
    i_rc = NF90_INQUIRE_ATTRIBUTE(f_e_o,i_v_o,a_n,attnum=i_a)
    IF ( (i_v_o == NF90_GLOBAL).AND.(i_rc /= NF90_NOERR) ) THEN
      nw_na(f_i_o) = nw_na(f_i_o)+1
    ENDIF
    CALL flio_hdm (f_i_o,f_e_o,.TRUE.)
    i_rc = NF90_COPY_ATT(f_e_i,i_v_i,a_n,f_e_o,i_v_o)
    IF (i_rc /= NF90_NOERR) THEN
      CALL ipslerr (3,'fliocpya', &
 &      'Attribute '//TRIM(a_n)//' not copied','Error :', &
 &      TRIM(NF90_STRERROR(i_rc)))
    ENDIF
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-fliocpya"
  ENDIF
!----------------------
END SUBROUTINE fliocpya
!===
SUBROUTINE flioqstc (f_i,c_type,l_ex,c_name)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i
  CHARACTER(LEN=*),INTENT(IN) :: c_type
  LOGICAL,INTENT(OUT) :: l_ex
  CHARACTER(LEN=*),INTENT(OUT) :: c_name
!-
  CHARACTER(LEN=1) :: c_ax
  INTEGER :: f_e,idc,ndc,i_rc
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->flioqstc ",TRIM(c_type)
  ENDIF
!-
! Retrieve the external file index
  CALL flio_qvid ('flioqstc',f_i,f_e)
!-
  c_ax = TRIM(c_type)
  IF (    (LEN_TRIM(c_type) == 1) &
 &    .AND.(    (c_ax == 'x').OR.(c_ax == 'y') &
 &          .OR.(c_ax == 'z').OR.(c_ax == 't')) ) THEN
    CALL flio_qax (f_i,c_ax,idc,ndc)
    l_ex = (idc > 0)
    IF (l_ex) THEN
      i_rc = NF90_INQUIRE_VARIABLE(f_e,idc,name=c_name)
    ENDIF
  ELSE
    l_ex = .FALSE.
    CALL ipslerr (2,'flioqstc', &
 &   'The name of the coordinate,',TRIM(c_type),'is not valid')
  ENDIF
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-flioqstc"
  ENDIF
!----------------------
END SUBROUTINE flioqstc
!===
SUBROUTINE fliosync (f_i)
!---------------------------------------------------------------------
  INTEGER,INTENT(in),OPTIONAL :: f_i
!-
  INTEGER :: i_f,f_e,i_rc,i_s,i_e
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->fliosync"
  ENDIF
!-
  IF (PRESENT(f_i)) THEN
    IF ( (f_i >= 1).AND.(f_i <= nb_fi_mx) ) THEN
      i_s = f_i
      i_e = f_i
    ELSE
      i_s = 1
      i_e = 0
      CALL ipslerr (2,'fliosync', &
 &     'Invalid file identifier',' ',' ')
    ENDIF
  ELSE
    i_s = 1
    i_e = nb_fi_mx
  ENDIF
!-
! Ensuring data mode
!-
  CALL flio_hdm (f_i,f_e,.FALSE.)
!-
  DO i_f=i_s,i_e
    f_e = nw_id(i_f)
    IF (f_e > 0) THEN
      IF (l_dbg) THEN
        WRITE(*,*) '  fliosync - synchronising file number ',i_f
      ENDIF
      i_rc = NF90_SYNC(f_e)
    ELSE IF (PRESENT(f_i)) THEN
      CALL ipslerr (2,'fliosync', &
 &     'Unable to synchronise the file :','probably','not opened')
    ENDIF
  ENDDO
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-fliosync"
  ENDIF
!----------------------
END SUBROUTINE fliosync
!===
SUBROUTINE flioclo (f_i)
!---------------------------------------------------------------------
  INTEGER,INTENT(in),OPTIONAL :: f_i
!-
  INTEGER :: i_f,f_e,i_rc,i_s,i_e
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->flioclo"
  ENDIF
!-
  IF (PRESENT(f_i)) THEN
    IF ( (f_i >= 1).AND.(f_i <= nb_fi_mx) ) THEN
      i_s = f_i
      i_e = f_i
    ELSE
      i_s = 1
      i_e = 0
      CALL ipslerr (2,'flioclo', &
 &     'Invalid file identifier',' ',' ')
    ENDIF
  ELSE
    i_s = 1
    i_e = nb_fi_mx
  ENDIF
!-
  DO i_f=i_s,i_e
    f_e = nw_id(i_f)
    IF (f_e > 0) THEN
      IF (l_dbg) THEN
        WRITE(*,*) '  flioclo - closing file number ',i_f
      ENDIF
      i_rc = NF90_CLOSE(f_e)
      nw_id(i_f) = -1
    ELSE IF (PRESENT(f_i)) THEN
      CALL ipslerr (2,'flioclo', &
 &     'Unable to close the file :','probably','not opened')
    ENDIF
  ENDDO
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-flioclo"
  ENDIF
!---------------------
END SUBROUTINE flioclo
!===
SUBROUTINE fliodmpf (f_n)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(IN) :: f_n
!-
  INTEGER :: f_e,n_dims,n_vars,n_atts,i_unlm
  INTEGER :: i_rc,i_n,k_n,t_ea,l_ea
  INTEGER :: tmp_i
  REAL    :: tmp_r
  INTEGER,DIMENSION(:),ALLOCATABLE :: tma_i
  REAL,DIMENSION(:),ALLOCATABLE    :: tma_r
  CHARACTER(LEN=256) :: tmp_c
  INTEGER,DIMENSION(nb_fd_mx) :: n_idim,n_ldim
  INTEGER,DIMENSION(nb_ax_mx) :: n_ai
  CHARACTER(LEN=NF90_MAX_NAME),DIMENSION(nb_fd_mx) :: c_ndim
  INTEGER,DIMENSION(NF90_MAX_VAR_DIMS) :: idimid
  CHARACTER(LEN=NF90_MAX_NAME) :: c_name
!---------------------------------------------------------------------
  i_rc = NF90_OPEN(TRIM(f_n),NF90_NOWRITE,f_e)
  IF (i_rc /= NF90_NOERR) THEN
    CALL ipslerr (3,'fliodmpf', &
 &   'Could not open file :',TRIM(f_n), &
 &   TRIM(NF90_STRERROR(i_rc))//' (Netcdf)')
  ENDIF
!-
  WRITE (*,*) "---"
  WRITE (*,*) "--- File '",TRIM(f_n),"'"
  WRITE (*,*) "---"
!-
  CALL flio_inf &
 &  (f_e,nb_dims=n_dims,nb_vars=n_vars, &
 &       nb_atts=n_atts,id_unlm=i_unlm, &
 &       nn_idm=n_idim,nn_ldm=n_ldim,cc_ndm=c_ndim,nn_aid=n_ai)
!-
  WRITE (*,*) 'External model identifier   : ',f_e
  WRITE (*,*) 'Number of dimensions        : ',n_dims
  WRITE (*,*) 'Number of variables         : ',n_vars
  WRITE (*,*) 'ID unlimited                : ',i_unlm
!-
  WRITE (*,*) "---"
  WRITE (*,*) 'Presumed axis dimensions identifiers :'
  IF (n_ai(k_lon) > 0) THEN
    WRITE (*,*) 'x axis : ',n_ai(k_lon)
  ELSE
    WRITE (*,*) 'x axis : NONE'
  ENDIF
  IF (n_ai(k_lat) > 0) THEN
    WRITE (*,*) 'y axis : ',n_ai(k_lat)
  ELSE
    WRITE (*,*) 'y axis : NONE'
  ENDIF
  IF (n_ai(k_lev) > 0) THEN
    WRITE (*,*) 'z axis : ',n_ai(k_lev)
  ELSE
    WRITE (*,*) 'z axis : NONE'
  ENDIF
  IF (n_ai(k_tim) > 0) THEN
    WRITE (*,*) 't axis : ',n_ai(k_tim)
  ELSE
    WRITE (*,*) 't axis : NONE'
  ENDIF
!-
  WRITE (*,*) "---"
  WRITE (*,*) 'Number of global attributes : ',n_atts
  DO k_n=1,n_atts
    i_rc = NF90_INQ_ATTNAME(f_e,NF90_GLOBAL,k_n,c_name)
    i_rc = NF90_INQUIRE_ATTRIBUTE(f_e,NF90_GLOBAL,c_name, &
 &                                xtype=t_ea,len=l_ea)
    IF      (    (t_ea == NF90_INT4).OR.(t_ea == NF90_INT2) &
             .OR.(t_ea == NF90_INT1) ) THEN
      IF (l_ea > 1) THEN
        ALLOCATE(tma_i(l_ea))
        i_rc = NF90_GET_ATT(f_e,NF90_GLOBAL,c_name,tma_i)
        WRITE (*,'("    ",A," :",/,(5(1X,I10),:))') &
 &        TRIM(c_name),tma_i(1:l_ea)
        DEALLOCATE(tma_i)
      ELSE
        i_rc = NF90_GET_ATT(f_e,NF90_GLOBAL,c_name,tmp_i)
        WRITE(*,*) '   ',TRIM(c_name),' : ',tmp_i
      ENDIF
    ELSE IF ( (t_ea == NF90_REAL4).OR.(t_ea == NF90_REAL8) ) THEN
      IF (l_ea > 1) THEN
        ALLOCATE(tma_r(l_ea))
        i_rc = NF90_GET_ATT(f_e,NF90_GLOBAL,c_name,tma_r)
        WRITE (*,'("    ",A," :",/,(5(1X,1PE11.3),:))') &
 &        TRIM(c_name),tma_r(1:l_ea)
        DEALLOCATE(tma_r)
      ELSE
        i_rc = NF90_GET_ATT(f_e,NF90_GLOBAL,c_name,tmp_r)
        WRITE(*,*) '   ',TRIM(c_name),' : ',tmp_r
      ENDIF
    ELSE
      tmp_c = ''
      i_rc = NF90_GET_ATT(f_e,NF90_GLOBAL,c_name,tmp_c)
      WRITE(*,*) '   ',TRIM(c_name),' : "',TRIM(tmp_c),'"'
    ENDIF
  ENDDO
!-
  DO i_n=1,nb_fd_mx
    IF (n_idim(i_n) > 0) THEN
      WRITE (*,*) "---"
      WRITE (*,*) 'Dimension id   : ',n_idim(i_n)
      WRITE (*,*) 'Dimension name : ',TRIM(c_ndim(i_n))
      WRITE (*,*) 'Dimension size : ',n_ldim(i_n)
    ENDIF
  ENDDO
!-
  DO i_n=1,n_vars
    i_rc = NF90_INQUIRE_VARIABLE(f_e,i_n, &
 &           name=c_name,ndims=n_dims,dimids=idimid,nAtts=n_atts)
    WRITE (*,*) "---"
    WRITE (*,*) "Variable name        : ",TRIM(c_name)
    WRITE (*,*) "Variable identifier  : ",i_n
    WRITE (*,*) "Number of dimensions : ",n_dims
    IF (n_dims > 0) THEN
      WRITE (*,*) "Dimensions ID's      : ",idimid(1:n_dims)
    ENDIF
    WRITE (*,*) "Number of attributes : ",n_atts
    DO k_n=1,n_atts
      i_rc = NF90_INQ_ATTNAME(f_e,i_n,k_n,c_name)
      i_rc = NF90_INQUIRE_ATTRIBUTE(f_e,i_n,c_name, &
 &                                  xtype=t_ea,len=l_ea)
      IF      (    (t_ea == NF90_INT4).OR.(t_ea == NF90_INT2) &
 &             .OR.(t_ea == NF90_INT1) ) THEN
        IF (l_ea > 1) THEN
          ALLOCATE(tma_i(l_ea))
          i_rc = NF90_GET_ATT(f_e,i_n,c_name,tma_i)
          WRITE (*,'("    ",A," :",/,(5(1X,I10),:))') &
 &              TRIM(c_name),tma_i(1:l_ea)
          DEALLOCATE(tma_i)
        ELSE
          i_rc = NF90_GET_ATT(f_e,i_n,c_name,tmp_i)
          WRITE(*,*) '   ',TRIM(c_name),' : ',tmp_i
        ENDIF
      ELSE IF ( (t_ea == NF90_REAL4).OR.(t_ea == NF90_REAL8) ) THEN
        IF (l_ea > 1) THEN
          ALLOCATE(tma_r(l_ea))
          i_rc = NF90_GET_ATT(f_e,i_n,c_name,tma_r)
          WRITE (*,'("    ",A," :",/,(5(1X,1PE11.3),:))') &
 &          TRIM(c_name),tma_r(1:l_ea)
          DEALLOCATE(tma_r)
        ELSE
          i_rc = NF90_GET_ATT(f_e,i_n,c_name,tmp_r)
          WRITE(*,*) '   ',TRIM(c_name),' : ',tmp_r
        ENDIF
      ELSE
        tmp_c = ''
        i_rc = NF90_GET_ATT(f_e,i_n,c_name,tmp_c)
        WRITE(*,*) '   ',TRIM(c_name),' : "',TRIM(tmp_c),'"'
      ENDIF
    ENDDO
  ENDDO
  WRITE (*,*) "---"
!-
  i_rc = NF90_CLOSE(f_e)
!----------------------
END SUBROUTINE fliodmpf
!===
SUBROUTINE flio_dom_set &
 & (dtnb,dnb,did,dsg,dsl,dpf,dpl,dhs,dhe,cdnm,id_dom)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: dtnb,dnb
  INTEGER,DIMENSION(:),INTENT(IN) :: did,dsg,dsl,dpf,dpl,dhs,dhe
  CHARACTER(LEN=*),INTENT(IN) :: cdnm
  INTEGER,INTENT(OUT) :: id_dom
!-
  INTEGER :: k_w,i_w,i_s
  CHARACTER(LEN=l_dns) :: cd_p,cd_w
!---------------------------------------------------------------------
  k_w = flio_dom_rid()
  IF (k_w < 0) THEN
    CALL ipslerr (3,'flio_dom_set', &
 &   'too many domains simultaneously defined', &
 &   'please unset useless domains', &
 &   'by calling flio_dom_unset')
  ENDIF
  id_dom = k_w
!-
  d_n_t(k_w) = dtnb
  d_n_c(k_w) = dnb
!-
  i_s = SIZE(did)
  IF (i_s > dom_max_dims) THEN
    CALL ipslerr (3,'flio_dom_set', &
 &   'too many distributed dimensions', &
 &   'simultaneously defined',' ')
  ENDIF
  d_d_n(k_w) = i_s
  d_d_i(1:i_s,k_w) = did(1:i_s)
!-
  i_w = SIZE(dsg)
  IF (i_w /= i_s) THEN
    CALL ipslerr (3,'flio_dom_set', &
 &   'the size of the DOMAIN_size_global array', &
 &   'is not equal to the size', &
 &   'of the distributed dimensions array')
  ENDIF
  d_s_g(1:i_w,k_w) = dsg(1:i_w)
!-
  i_w = SIZE(dsl)
  IF (i_w /= i_s) THEN
    CALL ipslerr (3,'flio_dom_set', &
 &   'the size of the DOMAIN_size_local array', &
 &   'is not equal to the size', &
 &   'of the distributed dimensions array')
  ENDIF
  d_s_l(1:i_w,k_w) = dsl(1:i_w)
!-
  i_w = SIZE(dpf)
  IF (i_w /= i_s) THEN
    CALL ipslerr (3,'flio_dom_set', &
 &   'the size of the DOMAIN_position_first array', &
 &   'is not equal to the size', &
 &   'of the distributed dimensions array')
  ENDIF
  d_p_f(1:i_w,k_w) = dpf(1:i_w)
!-
  i_w = SIZE(dpl)
  IF (i_w /= i_s) THEN
    CALL ipslerr (3,'flio_dom_set', &
 &   'the size of the DOMAIN_position_last array', &
 &   'is not equal to the size', &
 &   'of the distributed dimensions array')
  ENDIF
  d_p_l(1:i_w,k_w) = dpl(1:i_w)
!-
  i_w = SIZE(dhs)
  IF (i_w /= i_s) THEN
    CALL ipslerr (3,'flio_dom_set', &
 &   'the size of the DOMAIN_halo_size_start array', &
 &   'is not equal to the size', &
 &   'of the distributed dimensions array')
  ENDIF
  d_h_s(1:i_w,k_w) = dhs(1:i_w)
!-
  i_w = SIZE(dhe)
  IF (i_w /= i_s) THEN
    CALL ipslerr (3,'flio_dom_set', &
 &   'the size of the DOMAIN_halo_size_end array', &
 &   'is not equal to the size', &
 &   'of the distributed dimensions array')
  ENDIF
  d_h_e(1:i_w,k_w) = dhe(1:i_w)
!-
  cd_p = "unknown"
  cd_w = cdnm; CALL strlowercase (cd_w)
  DO i_w=1,n_dns
    IF (TRIM(cd_w) == TRIM(c_dns(i_w))) THEN
      cd_p = cd_w; EXIT;
    ENDIF
  ENDDO
  IF (TRIM(cd_p) == "unknown") THEN
    CALL ipslerr (3,'flio_dom_set', &
 &   'DOMAIN_type "'//TRIM(cdnm)//'"', &
 &   'is actually not supported', &
 &   'please use one of the supported names')
  ENDIF
  c_d_t(k_w) = cd_p
!--------------------------
END SUBROUTINE flio_dom_set
!===
SUBROUTINE flio_dom_unset (id_dom)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN),OPTIONAL :: id_dom
!-
  INTEGER :: i_w
!---------------------------------------------------------------------
  IF (PRESENT(id_dom)) THEN
    IF ( (id_dom >= 1).AND.(id_dom <= dom_max_nb) ) THEN
      IF (d_d_n(id_dom) > 0) THEN
        d_d_n(id_dom) = -1
      ELSE
        CALL ipslerr (2,'flio_dom_unset', &
 &       'The domain is not set',' ',' ')
      ENDIF
    ELSE
      CALL ipslerr (2,'flio_dom_unset', &
 &     'Invalid file identifier',' ',' ')
    ENDIF
  ELSE
    DO i_w=1,dom_max_nb
      d_d_n(id_dom) = -1
    ENDDO
  ENDIF
!----------------------------
END SUBROUTINE flio_dom_unset
!===
SUBROUTINE flio_dom_defset (id_dom)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: id_dom
!---------------------------------------------------------------------
  IF ( (id_dom >= 1).AND.(id_dom <= dom_max_nb) ) THEN
    id_def_dom = id_dom
  ELSE
    CALL ipslerr (3,'flio_dom_defset', &
 &   'Invalid domain identifier',' ',' ')
  ENDIF
!-----------------------------
END SUBROUTINE flio_dom_defset
!===
SUBROUTINE flio_dom_defunset ()
!---------------------------------------------------------------------
  IMPLICIT NONE
!---------------------------------------------------------------------
  id_def_dom = FLIO_DOM_NONE
!-------------------------------
END SUBROUTINE flio_dom_defunset
!===
SUBROUTINE flio_dom_definq (id_dom)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(OUT) :: id_dom
!---------------------------------------------------------------------
  id_dom = id_def_dom
!-----------------------------
END SUBROUTINE flio_dom_definq
!===
!-
!---------------------------------------------------------------------
!- Semi-public procedures
!---------------------------------------------------------------------
!-
!===
SUBROUTINE flio_dom_file (f_n,id_dom)
!---------------------------------------------------------------------
!- Update the model file name to include the ".nc" suffix and
!- the DOMAIN number on which this copy of IOIPSL runs, if needed.
!- This routine is called by IOIPSL and not by user anyway.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(INOUT) :: f_n
  INTEGER,OPTIONAL,INTENT(IN) :: id_dom
!-
  INTEGER :: il,iw
  CHARACTER(LEN=4) :: str
!---------------------------------------------------------------------
!-
! Add the ".nc" suffix if needed
  il = LEN_TRIM(f_n)
  IF (f_n(il-2:il) /= '.nc') THEN
    f_n = f_n(1:il)//'.nc'
  ENDIF
!-
! Add the DOMAIN identifier if needed
  IF (PRESENT(id_dom)) THEN
    IF (id_dom == FLIO_DOM_DEFAULT) THEN
      CALL flio_dom_definq (iw)
    ELSE
      iw = id_dom
    ENDIF
    IF (iw /= FLIO_DOM_NONE) THEN
      IF ( (id_dom >= 1).AND.(id_dom <= dom_max_nb) ) THEN
        IF (d_d_n(iw) > 0) THEN
          WRITE(str,'(I4.4)') d_n_c(iw)
          il = INDEX(f_n,'.nc')
          f_n = f_n(1:il-1)//'_'//str//'.nc'
        ELSE
          CALL ipslerr (3,'flio_dom_file', &
 &         'The domain has not been defined', &
 &         'please call flio_dom_set', &
 &         'before calling flio_dom_file')
        ENDIF
      ELSE
        CALL ipslerr (3,'flio_dom_file', &
 &       'Invalid domain identifier',' ',' ')
      ENDIF
    ENDIF
  ENDIF
!---------------------------
END SUBROUTINE flio_dom_file
!===
SUBROUTINE flio_dom_att (f_e,id_dom)
!---------------------------------------------------------------------
!- Add the DOMAIN attributes to the NETCDF file.
!- This routine is called by IOIPSL and not by user anyway.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(in) :: f_e
  INTEGER,OPTIONAL,INTENT(IN) :: id_dom
!-
  INTEGER :: iw,i_rc,i_n
  CHARACTER(LEN=15) :: c_ddim
  INTEGER :: n_idim
  CHARACTER(LEN=NF90_MAX_NAME) :: c_ndim
!---------------------------------------------------------------------
  IF (PRESENT(id_dom)) THEN
    IF (id_dom == FLIO_DOM_DEFAULT) THEN
      CALL flio_dom_definq (iw)
    ELSE
      iw = id_dom
    ENDIF
    IF (iw /= FLIO_DOM_NONE) THEN
      IF ( (id_dom >= 1).AND.(id_dom <= dom_max_nb) ) THEN
        IF (d_d_n(iw) > 0) THEN
          i_rc = NF90_PUT_ATT (f_e,NF90_GLOBAL, &
 &          'DOMAIN_number_total',d_n_t(iw))
          i_rc = NF90_PUT_ATT (f_e,NF90_GLOBAL, &
 &          'DOMAIN_number',d_n_c(iw))
          i_rc = NF90_PUT_ATT (f_e,NF90_GLOBAL, &
 &          'DOMAIN_dimensions_ids',d_d_i(1:d_d_n(iw),iw))
          i_rc = NF90_PUT_ATT (f_e,NF90_GLOBAL, &
 &          'DOMAIN_size_global',d_s_g(1:d_d_n(iw),iw))
          i_rc = NF90_PUT_ATT (f_e,NF90_GLOBAL, &
 &          'DOMAIN_size_local',d_s_l(1:d_d_n(iw),iw))
          i_rc = NF90_PUT_ATT (f_e,NF90_GLOBAL, &
 &          'DOMAIN_position_first',d_p_f(1:d_d_n(iw),iw))
          i_rc = NF90_PUT_ATT (f_e,NF90_GLOBAL, &
 &          'DOMAIN_position_last',d_p_l(1:d_d_n(iw),iw))
          i_rc = NF90_PUT_ATT (f_e,NF90_GLOBAL, &
 &          'DOMAIN_halo_size_start',d_h_s(1:d_d_n(iw),iw))
          i_rc = NF90_PUT_ATT (f_e,NF90_GLOBAL, &
 &          'DOMAIN_halo_size_end',d_h_e(1:d_d_n(iw),iw))
          i_rc = NF90_PUT_ATT (f_e,NF90_GLOBAL, &
 &          'DOMAIN_type',TRIM(c_d_t(iw)))
          i_rc = NF90_INQUIRE (f_e,nDimensions=n_idim)
          DO i_n=1,n_idim
            i_rc = NF90_INQUIRE_DIMENSION (f_e,i_n,name=c_ndim)
            WRITE (UNIT=c_ddim,FMT='("DOMAIN_DIM_N",I3.3)') i_n
            i_rc = NF90_PUT_ATT (f_e,NF90_GLOBAL,c_ddim,TRIM(c_ndim))
          ENDDO
        ELSE
          CALL ipslerr (3,'flio_dom_att', &
 &         'The domain has not been defined', &
 &         'please call flio_dom_set', &
 &         'before calling flio_dom_att')
        ENDIF
      ELSE
        CALL ipslerr (3,'flio_dom_att', &
 &       'Invalid domain identifier',' ',' ')
      ENDIF
    ENDIF
  ENDIF
!--------------------------
END SUBROUTINE flio_dom_att
!===
!-
!---------------------------------------------------------------------
!- Local procedures
!---------------------------------------------------------------------
!-
!===
INTEGER FUNCTION flio_rid()
!---------------------------------------------------------------------
!- returns a free index in nw_id(:)
!---------------------------------------------------------------------
  INTEGER,DIMENSION(1:1) :: nfi
!-
  IF (ANY(nw_id < 0)) THEN
    nfi = MINLOC(nw_id,MASK=nw_id < 0)
    flio_rid = nfi(1)
  ELSE
    flio_rid = -1
  ENDIF
!--------------------
END FUNCTION flio_rid
!===
INTEGER FUNCTION flio_dom_rid()
!---------------------------------------------------------------------
!- returns a free index in d_d_n(:)
!---------------------------------------------------------------------
  INTEGER,DIMENSION(1:1) :: nd
!---------------------------------------------------------------------
  IF (ANY(d_d_n < 0)) THEN
    nd = MINLOC(d_d_n,MASK=d_d_n < 0)
    flio_dom_rid = nd(1)
  ELSE
    flio_dom_rid = -1
  ENDIF
!------------------------
END FUNCTION flio_dom_rid
!===
INTEGER FUNCTION flio_qid(iid)
!---------------------------------------------------------------------
!- returns the external index associated with the internal index "iid"
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: iid
!---------------------------------------------------------------------
  IF ( (iid >= 1).AND.(iid <= nb_fi_mx) ) THEN
    flio_qid = nw_id(iid)
  ELSE
    flio_qid = -1
  ENDIF
!--------------------
END FUNCTION flio_qid
!===
SUBROUTINE flio_qvid (cpg,iid,ixd)
!---------------------------------------------------------------------
!- This subroutine, called by the procedure "cpg",
!- validates and returns the external file index "ixd"
!- associated with the internal file index "iid"
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(IN) :: cpg
  INTEGER,INTENT(IN)  :: iid
  INTEGER,INTENT(OUT) :: ixd
!-
  CHARACTER(LEN=20) :: c_t
!---------------------------------------------------------------------
  ixd = flio_qid(iid)
  IF (ixd < 0) THEN
    WRITE (UNIT=c_t,FMT='(I15)') iid
    CALL ipslerr (3,TRIM(cpg), &
 &    'Invalid internal file index :',TRIM(ADJUSTL(c_t)),' ')
  ENDIF
!-----------------------
END SUBROUTINE flio_qvid
!===
SUBROUTINE flio_hdm (f_i,f_e,lk_hm)
!---------------------------------------------------------------------
!- This subroutine handles the "define/data mode" of NETCDF.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_i,f_e
  LOGICAL,INTENT(IN) :: lk_hm
!-
  INTEGER :: i_rc
!---------------------------------------------------------------------
  i_rc = NF90_NOERR
!-
  IF      ( (.NOT.lw_hm(f_i)).AND.(lk_hm) ) THEN
    i_rc = NF90_REDEF(f_e)
    lw_hm(f_i) = .TRUE.
  ELSE IF ( (lw_hm(f_i)).AND.(.NOT.lk_hm) ) THEN
    i_rc = NF90_ENDDEF(f_e)
    lw_hm(f_i) = .FALSE.
  ENDIF
!-
  IF (i_rc /= NF90_NOERR) THEN
    CALL ipslerr (3,'flio_hdm', &
 &    'Internal error ','in define/data mode :', &
 &    TRIM(NF90_STRERROR(i_rc)))
  ENDIF
!----------------------
END SUBROUTINE flio_hdm
!===
SUBROUTINE flio_inf (f_e, &
 & nb_dims,nb_vars,nb_atts,id_unlm,nn_idm,nn_ldm,nn_aid,cc_ndm)
!---------------------------------------------------------------------
!- This subroutine allows to get some information concerning
!- the model file whose the external identifier is "f_e".
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: f_e
  INTEGER,OPTIONAL,INTENT(OUT) :: nb_dims,nb_vars,nb_atts,id_unlm
  INTEGER,DIMENSION(:),OPTIONAL,INTENT(OUT) :: nn_idm,nn_ldm,nn_aid
  CHARACTER(LEN=*),DIMENSION(:),OPTIONAL,INTENT(OUT) :: cc_ndm
!-
  INTEGER :: nm_dims,nm_vars,nm_atts,nm_unlm,ml
  INTEGER :: i_rc,kv
  CHARACTER(LEN=NF90_MAX_NAME) :: f_d_n
!-
  LOGICAL :: l_dbg
!---------------------------------------------------------------------
  CALL ipsldbg (old_status=l_dbg)
!-
  IF (l_dbg) THEN
    WRITE(*,*) "->flio_inf"
  ENDIF
!-
  i_rc = NF90_INQUIRE(f_e,nDimensions=nm_dims,nVariables=nm_vars, &
 &                    nAttributes=nm_atts,unlimitedDimId=nm_unlm)
!-
  IF (PRESENT(nb_dims))  nb_dims = nm_dims;
  IF (PRESENT(nb_vars))  nb_vars = nm_vars;
  IF (PRESENT(nb_atts))  nb_atts = nm_atts;
  IF (PRESENT(id_unlm))  id_unlm = nm_unlm;
!-
  IF (PRESENT(nn_idm))  nn_idm(:) =  -1;
  IF (PRESENT(nn_ldm))  nn_ldm(:) =   0;
  IF (PRESENT(cc_ndm))  cc_ndm(:) = ' ';
  IF (PRESENT(nn_aid))  nn_aid(:) =  -1;
!-
  DO kv=1,nm_dims
!---
    i_rc = NF90_INQUIRE_DIMENSION(f_e,kv,name=f_d_n,len=ml)
    CALL strlowercase (f_d_n)
    f_d_n = ADJUSTL(f_d_n)
!---
    IF (l_dbg) THEN
      WRITE(*,*) "  flio_inf ",kv,ml," ",TRIM(f_d_n)
    ENDIF
!---
    IF (PRESENT(nn_idm))  nn_idm(kv)=kv;
    IF (PRESENT(nn_ldm))  nn_ldm(kv)=ml;
    IF (PRESENT(cc_ndm))  cc_ndm(kv)=TRIM(f_d_n);
!---
    IF      (    (INDEX(f_d_n,'x') == 1)   &
 &           .OR.(INDEX(f_d_n,'lon') == 1) ) THEN
      IF (PRESENT(nn_aid)) THEN
        IF (nn_aid(k_lon) < 0) THEN
          nn_aid(k_lon)=kv;
        ENDIF
      ENDIF
    ELSE IF (    (INDEX(f_d_n,'y') == 1)   &
 &           .OR.(INDEX(f_d_n,'lat') == 1) ) THEN
      IF (PRESENT(nn_aid)) THEN
        IF (nn_aid(k_lat) < 0) THEN
          nn_aid(k_lat)=kv;
        ENDIF
      ENDIF
    ELSE IF (    (INDEX(f_d_n,'z') == 1)     &
 &           .OR.(INDEX(f_d_n,'lev') == 1)   &
 &           .OR.(INDEX(f_d_n,'plev') == 1)  &
 &           .OR.(INDEX(f_d_n,'depth') == 1) ) THEN
      IF (PRESENT(nn_aid)) THEN
        IF (nn_aid(k_lev) < 0) THEN
          nn_aid(k_lev)=kv;
        ENDIF
      ENDIF
    ELSE IF (    (TRIM(f_d_n) == 't')         &
 &           .OR.(TRIM(f_d_n) == 'time')      &
 &           .OR.(INDEX(f_d_n,'tstep') == 1)  &
 &           .OR.(INDEX(f_d_n,'time_counter') == 1) ) THEN
!---- For the time we certainly need to allow for other names
      IF (PRESENT(nn_aid)) THEN
        IF (nn_aid(k_tim) < 0) THEN
          nn_aid(k_tim)=kv;
        ENDIF
      ENDIF
    ENDIF
!---
  ENDDO
!-
  IF (l_dbg) THEN
    WRITE(*,*) "<-flio_inf"
  ENDIF
!----------------------
END SUBROUTINE flio_inf
!===
SUBROUTINE flio_qax (f_i,axtype,i_v,nbd)
!---------------------------------------------------------------------
!- This subroutine explores the file in order to find
!- an axis (x/y/z/t) according to a number of rules
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER :: f_i,i_v,nbd
  CHARACTER(LEN=*) :: axtype
!-
  INTEGER :: kv,k,n_r,l_d,n_d,i_rc,dimnb
  CHARACTER(LEN=1)  :: c_ax
  CHARACTER(LEN=18) :: c_sn
  CHARACTER(LEN=15),DIMENSION(10) :: c_r
  CHARACTER(LEN=40) :: c_t1,c_t2
!---------------------------------------------------------------------
  i_v = -1; nbd = -1;
!---
!- Keep the name of the axis
!---
  c_ax = TRIM(axtype)
!-
! Validate axis type
!-
  IF (    (LEN_TRIM(axtype) == 1) &
 &    .AND.(    (c_ax == 'x').OR.(c_ax == 'y') &
 &          .OR.(c_ax == 'z').OR.(c_ax == 't')) ) THEN
!---
!-- Define the maximum number of dimensions for the coordinate
!---
      SELECT CASE (c_ax)
      CASE('x')
        l_d = 2
        c_sn = 'longitude'
      CASE('y')
        l_d = 2
        c_sn = 'latitude'
      CASE('z')
        l_d = 1
        c_sn = 'model_level_number'
      CASE('t')
        l_d = 1
        c_sn = 'time'
      END SELECT
!---
!-- Rule 1 : we look for a variable with one dimension
!--          and which has the same name as its dimension (NUG)
!---
    IF (i_v < 0) THEN
      SELECT CASE (c_ax)
      CASE('x')
        k = nw_ai(k_lon,f_i)
      CASE('y')
        k = nw_ai(k_lat,f_i)
      CASE('z')
        k = nw_ai(k_lev,f_i)
      CASE('t')
        k = nw_ai(k_tim,f_i)
      END SELECT
      IF ( (k >= 1).AND.(k <= nb_ax_mx) ) THEN
        dimnb = nw_di(k,f_i)
      ELSE
        dimnb = -1
      ENDIF
!-----
      i_rc = NF90_INQUIRE_DIMENSION(nw_id(f_i),dimnb,name=c_t1)
      IF (i_rc == NF90_NOERR) THEN
        CALL strlowercase (c_t1)
        L_R1: DO kv=1,nw_nv(f_i)
          i_rc = NF90_INQUIRE_VARIABLE &
 &                 (nw_id(f_i),kv,name=c_t2,ndims=n_d)
          IF (n_d == 1) THEN
            CALL strlowercase (c_t2)
            IF (TRIM(c_t1) == TRIM(c_t2)) THEN
              i_v = kv; nbd = n_d;
              EXIT L_R1
            ENDIF
          ENDIF
        ENDDO L_R1
      ENDIF
    ENDIF
!---
!-- Rule 2 : we look for a correct "axis" attribute (CF)
!---
    IF (i_v < 0) THEN
      L_R2: DO kv=1,nw_nv(f_i)
        i_rc = NF90_GET_ATT(nw_id(f_i),kv,'axis',c_t1)
        IF (i_rc == NF90_NOERR) THEN
          CALL strlowercase (c_t1)
          IF (TRIM(c_t1) == c_ax) THEN
            i_rc = NF90_INQUIRE_VARIABLE(nw_id(f_i),kv,ndims=n_d)
            IF (n_d <= l_d) THEN
              i_v = kv; nbd = n_d;
              EXIT L_R2
            ENDIF
          ENDIF
        ENDIF
      ENDDO L_R2
    ENDIF
!---
!-- Rule 3 : we look for a correct "standard_name" attribute (CF)
!---
    IF (i_v < 0) THEN
      L_R3: DO kv=1,nw_nv(f_i)
        i_rc = NF90_GET_ATT(nw_id(f_i),kv,'standard_name',c_t1)
        IF (i_rc == NF90_NOERR) THEN
          CALL strlowercase (c_t1)
          IF (TRIM(c_t1) == TRIM(c_sn)) THEN
            i_rc = NF90_INQUIRE_VARIABLE(nw_id(f_i),kv,ndims=n_d)
            IF (n_d <= l_d) THEN
              i_v = kv; nbd = n_d;
              EXIT L_R3
            ENDIF
          ENDIF
        ENDIF
      ENDDO L_R3
    ENDIF
!---
!-- Rule 4 : we look for a specific name (IOIPSL)
!---
    IF (i_v < 0) THEN
      SELECT CASE (c_ax)
      CASE('x')
        n_r = 3
        c_r(1)='nav_lon'; c_r(2)='lon'; c_r(3)='longitude';
      CASE('y')
        n_r = 3
        c_r(1)='nav_lat'; c_r(2)='lat'; c_r(3)='latitude';
      CASE('z')
        n_r = 8
        c_r(1)='depth'; c_r(2)='deptht'; c_r(3)='height';
        c_r(4)='level'; c_r(5)='lev'; c_r(6)='plev';
        c_r(7)='sigma_level'; c_r(8)='layer';
      CASE('t')
        n_r = 3
        c_r(1)='time'; c_r(2)='tstep'; c_r(3)='timesteps';
      END SELECT
!-----
      L_R4: DO kv=1,nw_nv(f_i)
        i_rc = NF90_INQUIRE_VARIABLE &
 &               (nw_id(f_i),kv,name=c_t1,ndims=n_d)
        IF (i_rc == NF90_NOERR) THEN
          CALL strlowercase (c_t1)
          IF (n_d <= l_d) THEN
            DO k=1,n_r
              IF (TRIM(c_t1) == TRIM(c_r(k))) THEN
                i_v = kv; nbd = n_d;
                EXIT L_R4
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDDO L_R4
    ENDIF
!---
  ENDIF
!----------------------
END SUBROUTINE flio_qax
!-
!===
!-
END MODULE fliocom
