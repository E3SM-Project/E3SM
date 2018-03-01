Module netcdf_nc_interfaces

! Fortran interfaces to netCDF C functions using FORTRAN 2003 C
! Interoperability features. These interfaces are for the base
! netCDF C routines in the libsrc directory

! Written by: Richard Weed, Ph.D.
!             Center for Advanced Vehicular Systems
!             Mississippi State University
!             rweed@cavs.msstate.edu


! License (and other Lawyer Language)

! This software is released under the Apache 2.0 Open Source License. The
! full text of the License can be viewed at :
!
!   http:www.apache.org/licenses/LICENSE-2.0.html
!
! The author grants to the University Corporation for Atmospheric Research
! (UCAR), Boulder, CO, USA the right to revise and extend the software
! without restriction. However, the author retains all copyrights and
! intellectual property rights explicitly stated in or implied by the
! Apache license

! Version 1.:  Sept. 2005 - Initial Cray X1 version
! Version 2.:  May   2006 - Updated to support g95
! Version 3.:  June  2006 - Updated to include netCDF 4 functions
! Version 4.:  April 2009 - Updated to match netCDF 4.0.1 release
! Version 5.:  April 2010 - Updated to match netCDF 4.1.1 release
! Version 6.:  April 2013 - Added nc_inq_path support for fortran 4.4 beta

 USE netcdf_nc_data

 Implicit NONE

!> module procedure interfaces for utility routines

 Interface addCNullChar
  module procedure addCNullChar
 End Interface

 Interface stripCNullChar
  module procedure stripCNullChar
 End Interface

!> Begin explicit interfaces for base nc_ functions. Note that some interfaces
!! expect data to be passed as C_PTR type variables. These data are arrays
!! that could have a NULL pointer passed instead of the array. All strings
!! are passed as C_CHAR interoperable strings. Most data in put routines
!! that map to a void pointer in C are passed as a Type(C_PTR) value.
!! Data from get routine that pass as a void pointer in C are passed as
!! a C_CHAR string.

!! Also note that each interface has an explicit USE ISO_C_BINDING. A better
!! solution is to use the F2003 IMPORT statement (I originally had it this way)
!! However its best to leave the interfaces as is for now because there might
!! be a few compilers out there that support most of the C-interop facility but
!! for some reason haven't implemented IMPORT yet.

!---------------------------------- nc_strerror -------------------------------
Interface
 Function nc_strerror(ncerr) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR

 Integer(KIND=C_INT), VALUE :: ncerr

 Type(C_PTR)                :: nc_strerror

 End Function nc_strerror
End Interface
!---------------------------------- nc_inq_libvers ----------------------------
Interface
 Function nc_inq_libvers() BIND(C)

 USE ISO_C_BINDING, ONLY: C_PTR

 Type(C_PTR) :: nc_inq_libvers

 End Function nc_inq_libvers
End Interface
!---------------------------------- nc_create ---------------------------------
Interface
 Function nc_create(path, cmode, ncidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_CHAR, C_INT

 Character(KIND=C_CHAR), Intent(IN)  :: path(*)
 Integer(KIND=C_INT),    VALUE       :: cmode
 Integer(KIND=C_INT),    Intent(OUT) :: ncidp

 Integer(KIND=C_INT)                 :: nc_create

 End Function nc_create
End Interface
!---------------------------------- nc__create --------------------------------
Interface
 Function nc__create(path, cmode, initialsz, chunksizehintp, ncidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_SIZE_T

 Character(KIND=C_CHAR), Intent(IN)  :: path(*)
 Integer(KIND=C_INT),    VALUE       :: cmode
 Integer(KIND=C_SIZE_T), VALUE       :: initialsz
 Integer(KIND=C_SIZE_T), Intent(IN)  :: chunksizehintp
 Integer(KIND=C_INT),    Intent(OUT) :: ncidp

 Integer(KIND=C_INT)                 :: nc__create

 End Function nc__create
End Interface
!---------------------------------- nc__create_mp -----------------------------
Interface
 Function nc__create_mp(path, cmode, initialsz, basepe, chunksizehintp, ncidp) &
                        BIND(C)

 USE ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_SIZE_T, C_PTR

 Character(KIND=C_CHAR), Intent(IN)  :: path(*)
 Integer(KIND=C_INT),    VALUE       :: cmode
 Integer(KIND=C_SIZE_T), VALUE       :: initialsz
 Integer(KIND=C_SIZE_T), Intent(IN)  :: chunksizehintp
 Type(C_PTR),            VALUE       :: basepe
 Integer(KIND=C_INT),    Intent(OUT) :: ncidp

 Integer(KIND=C_INT)                 :: nc__create_mp

 End Function nc__create_mp
End Interface
!---------------------------------- nc_open -----------------------------------
Interface
 Function nc_open(path, mode, ncidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_CHAR, C_INT

 Character(KIND=C_CHAR), Intent(IN)  :: path(*)
 Integer(KIND=C_INT),    VALUE       :: mode
 Integer(KIND=C_INT),    Intent(OUT) :: ncidp

 Integer(KIND=C_INT)                 :: nc_open

 End Function nc_open
End Interface
!---------------------------------- nc__open ----------------------------------
Interface
 Function nc__open(path, mode, chunksizehintp, ncidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_SIZE_T

 Character(KIND=C_CHAR), Intent(IN)  :: path(*)
 Integer(KIND=C_INT),    VALUE       :: mode
 Integer(KIND=C_SIZE_T), Intent(IN)  :: chunksizehintp
 Integer(KIND=C_INT),    Intent(OUT) :: ncidp

 Integer(KIND=C_INT)                 :: nc__open

 End Function nc__open
End Interface
!---------------------------------- nc__open_mp -------------------------------
Interface
 Function nc__open_mp(path, mode, basepe, chunksizehintp, ncidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_SIZE_T, C_PTR

 Character(KIND=C_CHAR), Intent(IN)  :: path(*)
 Integer(KIND=C_INT),    VALUE       :: mode
 Integer(KIND=C_SIZE_T), Intent(IN)  :: chunksizehintp
 Integer(KIND=C_INT),    Intent(OUT) :: ncidp
 Type(C_PTR),            VALUE       :: basepe

 Integer(KIND=C_INT)                 :: nc__open_mp

 End Function nc__open_mp
End Interface
!---------------------------------- nc_inq_path -----------------------------
Interface
 Function nc_inq_path(ncid, pathlen, path) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: pathlen
 Character(KIND=C_CHAR), Intent(INOUT) :: path(*)

 Integer(KIND=C_INT)                   :: nc_inq_path

 End Function nc_inq_path
End Interface
!---------------------------------- nc_set_fill -------------------------------
Interface
 Function nc_set_fill(ncid, fillmode, old_modep) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid
 Integer(KIND=C_INT), VALUE       :: fillmode
 Integer(KIND=C_INT), Intent(OUT) :: old_modep

 Integer(KIND=C_INT)              :: nc_set_fill

 End Function nc_set_fill
End Interface
!---------------------------------- nc_redef ----------------------------------
Interface
 Function nc_redef(ncid) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: ncid

 Integer(KIND=C_INT)        :: nc_redef

 End Function nc_redef
End Interface
!---------------------------------- nc_enddef ---------------------------------
Interface
 Function nc_enddef(ncid) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: ncid

 Integer(KIND=C_INT)        :: nc_enddef

 End Function nc_enddef
End Interface
!---------------------------------- nc__enddef --------------------------------
Interface
 Function nc__enddef(ncid, h_minfree, v_align, v_minfree, r_align) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T

 Integer(KIND=C_INT),    VALUE :: ncid
 Integer(KIND=C_SIZE_T), VALUE :: h_minfree, v_align, v_minfree, r_align

 Integer(KIND=C_INT)           :: nc__enddef

 End Function nc__enddef
End Interface
!---------------------------------- nc_sync -----------------------------------
Interface
 Function nc_sync(ncid) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: ncid

 Integer(KIND=C_INT)        :: nc_sync

 End Function nc_sync
End Interface
!---------------------------------- nc_abort ----------------------------------
Interface
 Function nc_abort(ncid) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: ncid

 Integer(KIND=C_INT)        :: nc_abort

 End Function nc_abort
End Interface
!---------------------------------- nc_close ----------------------------------
Interface
 Function nc_close(ncid) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: ncid

 Integer(KIND=C_INT)        :: nc_close

 End Function nc_close
End Interface
!---------------------------------- nc_delete --------------------------------
Interface
 Function nc_delete(path) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Character(KIND=C_CHAR), Intent(IN) :: path(*)

 Integer(KIND=C_INT)                :: nc_delete

 End Function nc_delete
End Interface
!---------------------------------- nc_delete_mp -----------------------------
Interface
 Function nc_delete_mp(path, pe) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Character(KIND=C_CHAR), Intent(IN) :: path(*)
 Integer(KIND=C_INT),    VALUE      :: pe

 Integer(KIND=C_INT)                :: nc_delete_mp

 End Function nc_delete_mp
End Interface
!---------------------------------- nc_set_base_pe ----------------------------
Interface
 Function nc_set_base_pe(ncid, pe) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: ncid, pe

 Integer(KIND=C_INT)        :: nc_set_base_pe

 End Function nc_set_base_pe
End Interface
!---------------------------------- nc_inq_base_pe ----------------------------
Interface
 Function nc_inq_base_pe(ncid, pe) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid
 Integer(KIND=C_INT), Intent(OUT) ::  pe

 Integer(KIND=C_INT)              :: nc_inq_base_pe

 End Function nc_inq_base_pe
End Interface
!---------------------------------- nc_inq ------------------------------------
Interface
 Function nc_inq(ncid, ndimsp, nvarsp, ngattsp, unlimdimidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid
 Integer(KIND=C_INT), Intent(OUT) :: ndimsp, nvarsp, ngattsp, unlimdimidp

 Integer(KIND=C_INT)              :: nc_inq

 End Function nc_inq
End Interface
!---------------------------------- nc_inq_ndims ------------------------------
Interface
 Function nc_inq_ndims(ncid, ndimsp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid
 Integer(KIND=C_INT), Intent(OUT) ::  ndimsp

 Integer(KIND=C_INT)              :: nc_inq_ndims

 End Function nc_inq_ndims
End Interface
!---------------------------------- nc_inq_nvars ------------------------------
Interface
 Function nc_inq_nvars(ncid, nvarsp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid
 Integer(KIND=C_INT), Intent(OUT) ::  nvarsp

 Integer(KIND=C_INT)              :: nc_inq_nvars

 End Function nc_inq_nvars
End Interface
!---------------------------------- nc_inq_natts ------------------------------
Interface
 Function nc_inq_natts(ncid, ngattsp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid
 Integer(KIND=C_INT), Intent(OUT) :: ngattsp

 Integer(KIND=C_INT)              :: nc_inq_natts

 End Function nc_inq_natts
End Interface
!---------------------------------- nc_inq_unlimdim ---------------------------
Interface
 Function nc_inq_unlimdim(ncid, unlimdimidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid
 Integer(KIND=C_INT), Intent(OUT) :: unlimdimidp

 Integer(KIND=C_INT)              :: nc_inq_unlimdim

 End Function nc_inq_unlimdim
End Interface
!---------------------------------- nc_inq_format -----------------------------
Interface
 Function nc_inq_format(ncid, formatp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid
 Integer(KIND=C_INT), Intent(OUT) :: formatp

 Integer(KIND=C_INT)              :: nc_inq_format

 End Function nc_inq_format
End Interface
!---------------------------------- nc_def_dim --------------------------------
Interface
 Function nc_def_dim(ncid, name, nlen, idp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Character(KIND=C_CHAR), Intent(IN)    :: name(*)
 Integer(KIND=C_SIZE_T), VALUE         :: nlen
 Integer(KIND=C_INT),    Intent(INOUT) :: idp

 Integer(KIND=C_INT)                   :: nc_def_dim

 End Function nc_def_dim
End Interface
!---------------------------------- nc_inq_dimid ------------------------------
Interface
 Function nc_inq_dimid(ncid, name, idp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Character(KIND=C_CHAR), Intent(IN)    :: name(*)
 Integer(KIND=C_INT),    Intent(INOUT) :: idp

 Integer(KIND=C_INT)                   :: nc_inq_dimid

 End Function nc_inq_dimid
End Interface
!---------------------------------- nc_inq_dim --------------------------------
Interface
 Function nc_inq_dim(ncid, dimid, name, lenp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_INT),    VALUE         :: dimid
 Character(KIND=C_CHAR), Intent(INOUT) :: name(*)
 Integer(KIND=C_SIZE_T), Intent(OUT)   :: lenp

 Integer(KIND=C_INT)                   :: nc_inq_dim

 End Function nc_inq_dim
End Interface
!---------------------------------- nc_inq_dimname ----------------------------
Interface
 Function nc_inq_dimname(ncid, dimid, name) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_INT),    VALUE         :: dimid
 Character(KIND=C_CHAR), Intent(INOUT) :: name(*)

 Integer(KIND=C_INT)                   :: nc_inq_dimname

 End Function nc_inq_dimname
End Interface
!---------------------------------- nc_inq_dimlen -----------------------------
Interface
 Function nc_inq_dimlen(ncid, dimid, lenp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T

 Integer(KIND=C_INT),    VALUE       :: ncid
 Integer(KIND=C_INT),    VALUE       :: dimid
 Integer(KIND=C_SIZE_T), Intent(OUT) :: lenp

 Integer(KIND=C_INT)                 :: nc_inq_dimlen

 End Function nc_inq_dimlen
End Interface
!---------------------------------- nc_rename_dim -----------------------------
Interface
 Function nc_rename_dim(ncid, dimid, name) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid
 Integer(KIND=C_INT),    VALUE      :: dimid
 Character(KIND=C_CHAR), Intent(IN) :: name(*)

 Integer(KIND=C_INT)                :: nc_rename_dim

 End Function nc_rename_dim
End Interface
!---------------------------------- nc_def_var --------------------------------
Interface
 Function nc_def_var(ncid, name, xtype, ndims, dimidsp, varidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Integer(KIND=C_INT),    VALUE       :: xtype
 Integer(KIND=C_INT),    VALUE       :: ndims
 Integer(KIND=C_INT),    Intent(IN)  :: dimidsp(*)
 Integer(KIND=C_INT),    Intent(OUT) :: varidp

 Integer(KIND=C_INT)                 :: nc_def_var

 End Function nc_def_var
End Interface
!---------------------------------- nc_inq_var --------------------------------
Interface
 Function nc_inq_var(ncid, varid, name, xtypep, ndimsp, dimidsp, nattsp) &
                     BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(OUT) :: name(*)
 Integer(KIND=C_INT),    Intent(OUT) :: xtypep
 Integer(KIND=C_INT),    Intent(OUT) :: ndimsp
 Integer(KIND=C_INT),    Intent(OUT) :: dimidsp(*)
 Integer(KIND=C_INT),    Intent(OUT) :: nattsp

 Integer(KIND=C_INT)                 :: nc_inq_var

 End Function nc_inq_var
End Interface
!---------------------------------- nc_inq_varid ------------------------------
Interface
 Function nc_inq_varid(ncid, name, varidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Integer(KIND=C_INT),    Intent(OUT) :: varidp

 Integer(KIND=C_INT)                 :: nc_inq_varid

 End Function nc_inq_varid
End Interface
!---------------------------------- nc_inq_varname ----------------------------
Interface
 Function nc_inq_varname(ncid, varid, name) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(OUT) :: name(*)

 Integer(KIND=C_INT)                 :: nc_inq_varname

 End Function nc_inq_varname
End Interface
!---------------------------------- nc_inq_vartype ----------------------------
Interface
 Function nc_inq_vartype(ncid, varid, xtypep) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Integer(KIND=C_INT), Intent(OUT) :: xtypep

 Integer(KIND=C_INT)              :: nc_inq_vartype

 End Function nc_inq_vartype
End Interface
!---------------------------------- nc_inq_varndims ---------------------------
Interface
 Function nc_inq_varndims(ncid, varid, ndimsp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Integer(KIND=C_INT), Intent(OUT) :: ndimsp

 Integer(KIND=C_INT)              :: nc_inq_varndims

 End Function nc_inq_varndims
End Interface
!---------------------------------- nc_inq_vardimid ---------------------------
Interface
 Function nc_inq_vardimid(ncid, varid, dimidsp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Integer(KIND=C_INT), Intent(OUT) :: dimidsp(*)

 Integer(KIND=C_INT)              :: nc_inq_vardimid

 End Function nc_inq_vardimid
End Interface
!---------------------------------- nc_inq_varnatts ---------------------------
Interface
 Function nc_inq_varnatts(ncid, varid, nattsp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Integer(KIND=C_INT), Intent(OUT) :: nattsp

 Integer(KIND=C_INT)              :: nc_inq_varnatts

 End Function nc_inq_varnatts
End Interface
!---------------------------------- nc_rename_var -----------------------------
Interface
 Function nc_rename_var(ncid, varid, name) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN) :: name(*)

 Integer(KIND=C_INT)                :: nc_rename_var

 End Function nc_rename_var
End Interface
!---------------------------------- nc_put_var_text ---------------------------
Interface
 Function nc_put_var_text(ncid, varid, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN) :: op(*)

 Integer(KIND=C_INT)                :: nc_put_var_text

 End Function nc_put_var_text
End Interface
!---------------------------------- nc_get_var_text ---------------------------
Interface
 Function nc_get_var_text(ncid, varid, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid, varid
 Character(KIND=C_CHAR), Intent(INOUT) :: ip(*)

 Integer(KIND=C_INT)                   :: nc_get_var_text

 End Function nc_get_var_text
End Interface
!---------------------------------- nc_put_var_uchar --------------------------
Interface
 Function nc_put_var_uchar(ncid, varid, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Integer(KIND=CINT1), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_var_uchar

 End Function nc_put_var_uchar
End Interface
!---------------------------------- nc_get_var_uchar --------------------------
Interface
 Function nc_get_var_uchar(ncid, varid, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Integer(KIND=CINT1), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_var_uchar

 End Function nc_get_var_uchar
End Interface
!---------------------------------- nc_put_var_schar --------------------------
Interface
 Function nc_put_var_schar(ncid, varid, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Integer(KIND=CINT1), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_var_schar

 End Function nc_put_var_schar
End Interface
!---------------------------------- nc_get_var_schar --------------------------
Interface
 Function nc_get_var_schar(ncid, varid, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT),  VALUE       :: ncid, varid
 Integer(KIND=CINT1),  Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)               :: nc_get_var_schar

 End Function nc_get_var_schar
End Interface
!---------------------------------- nc_put_var_short --------------------------
Interface
 Function nc_put_var_short(ncid, varid, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT
 USE NETCDF_NC_DATA, ONLY: CINT2

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Integer(KIND=CINT2), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_var_short

 End Function nc_put_var_short
End Interface
!---------------------------------- nc_get_var_short --------------------------
Interface
 Function nc_get_var_short(ncid, varid, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT
 USE NETCDF_NC_DATA, ONLY: CINT2

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Integer(KIND=CINT2), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_var_short

 End Function nc_get_var_short
End Interface
!---------------------------------- nc_put_var_int ----------------------------
Interface
 Function nc_put_var_int(ncid, varid, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT
 USE NETCDF_NC_DATA, ONLY: CINT

 Integer(KIND=C_INT), VALUE     :: ncid, varid
 Integer(KIND=CINT), Intent(IN) :: op(*)

 Integer(KIND=C_INT)            :: nc_put_var_int

 End Function nc_put_var_int
End Interface
!---------------------------------- nc_get_var_int ----------------------------
Interface
 Function nc_get_var_int(ncid, varid, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT
 USE NETCDF_NC_DATA, ONLY: CINT

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Integer(KIND=CINT),  Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_var_int

 End Function nc_get_var_int
End Interface
!---------------------------------- nc_put_var_long ----------------------------
Interface
 Function nc_put_var_long(ncid, varid, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG

 Integer(KIND=C_INT),  VALUE      :: ncid, varid
 Integer(KIND=C_LONG), Intent(IN) :: op(*)

 Integer(KIND=C_INT)              :: nc_put_var_long

 End Function nc_put_var_long
End Interface
!---------------------------------- nc_get_var_long ---------------------------
Interface
 Function nc_get_var_long(ncid, varid, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG

 Integer(KIND=C_INT),  VALUE       :: ncid, varid
 Integer(KIND=C_LONG), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)               :: nc_get_var_long

 End Function nc_get_var_long
End Interface
!---------------------------------- nc_put_var_float --------------------------
Interface
 Function nc_put_var_float(ncid, varid, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Real(KIND=C_FLOAT),  Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_var_float

 End Function nc_put_var_float
End Interface
!---------------------------------- nc_get_var_float --------------------------
Interface
 Function nc_get_var_float(ncid, varid, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT

 Implicit NONE

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Real(KIND=C_FLOAT),  Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_var_float

 End Function nc_get_var_float
End Interface
!---------------------------------- nc_put_var_double -------------------------
Interface
 Function nc_put_var_double(ncid, varid, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Real(KIND=C_DOUBLE), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_var_double

 End Function nc_put_var_double
End Interface
!---------------------------------- nc_get_var_double -------------------------
Interface
 Function nc_get_var_double(ncid, varid, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Real(KIND=C_DOUBLE), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_var_double

 End Function nc_get_var_double
End Interface
!---------------------------------- nc_put_var1_text --------------------------
Interface
 Function nc_put_var1_text(ncid, varid, indexp, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),   VALUE      :: ncid, varid
 Type(C_PTR),           VALUE      :: indexp
 Character(LEN=C_CHAR), Intent(IN) :: op

 Integer(KIND=C_INT)               :: nc_put_var1_text

 End Function nc_put_var1_text
End Interface
!---------------------------------- nc_get_var1_test --------------------------
Interface
 Function nc_get_var1_text(ncid, varid, indexp, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Type(C_PTR),            VALUE       :: indexp
 Character(KIND=C_CHAR), Intent(OUT) :: ip

 Integer(KIND=C_INT)                 :: nc_get_var1_text

 End Function nc_get_var1_text
End Interface
!---------------------------------- nc_put_var1_uchar -------------------------
Interface
 Function nc_put_var1_uchar(ncid, varid, indexp, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT),  VALUE      :: ncid, varid
 Type(C_PTR),          VALUE      :: indexp
 Integer(KIND=CINT1),  Intent(IN) :: op

 Integer(KIND=C_INT)              :: nc_put_var1_uchar

 End Function nc_put_var1_uchar
End Interface
!---------------------------------- nc_get_var1_uchar -------------------------
Interface
 Function nc_get_var1_uchar(ncid, varid, indexp, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT),  VALUE       :: ncid, varid
 Type(C_PTR),          VALUE       :: indexp
 Integer(KIND=CINT1),  Intent(OUT) :: ip

 Integer(KIND=C_INT)               :: nc_get_var1_uchar

 End Function nc_get_var1_uchar
End Interface
!---------------------------------- nc_put_var1_schar -------------------------
Interface
 Function nc_put_var1_schar(ncid, varid, indexp, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT),  VALUE      :: ncid, varid
 Type(C_PTR),          VALUE      :: indexp
 Integer(KIND=CINT1),  Intent(IN) :: op

 Integer(KIND=C_INT)              :: nc_put_var1_schar

 End Function nc_put_var1_schar
End Interface
!---------------------------------- nc_get_var1_schar -------------------------
Interface
 Function nc_get_var1_schar(ncid, varid, indexp, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT),  VALUE       :: ncid, varid
 Type(C_PTR),          VALUE       :: indexp
 Integer(KIND=CINT1),  Intent(OUT) :: ip

 Integer(KIND=C_INT)               :: nc_get_var1_schar

 End Function nc_get_var1_schar
End Interface
!---------------------------------- nc_put_var1_short -------------------------
Interface
 Function nc_put_var1_short(ncid, varid, indexp, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT2

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: indexp
 Integer(KIND=CINT2), Intent(IN) :: op

 Integer(KIND=C_INT)             :: nc_put_var1_short

 End Function nc_put_var1_short
End Interface
!---------------------------------- nc_get_var1_short -------------------------
Interface
 Function nc_get_var1_short(ncid, varid, indexp, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT2

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: indexp
 Integer(KIND=CINT2), Intent(OUT) :: ip

 Integer(KIND=C_INT)              :: nc_get_var1_short

 End Function nc_get_var1_short
End Interface
!---------------------------------- nc_put_var1_int ---------------------------
Interface
 Function nc_put_var1_int(ncid, varid, indexp, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: indexp
 Integer(KIND=C_INT), Intent(IN) :: op

 Integer(KIND=C_INT)             :: nc_put_var1_int

 End Function nc_put_var1_int
End Interface
!---------------------------------- nc_get_var1_int ---------------------------
Interface
 Function nc_get_var1_int(ncid, varid, indexp, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: indexp
 Integer(KIND=C_INT), Intent(OUT) :: ip

 Integer(KIND=C_INT)              :: nc_get_var1_int

 End Function nc_get_var1_int
End Interface
!---------------------------------- nc_put_var1_long --------------------------
Interface
 Function nc_put_var1_long(ncid, varid, indexp, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG, C_PTR

 Integer(KIND=C_INT),  VALUE      :: ncid, varid
 Type(C_PTR),          VALUE      :: indexp
 Integer(KIND=C_LONG), Intent(IN) :: op

 Integer(KIND=C_INT)              :: nc_put_var1_long

 End Function nc_put_var1_long
End Interface
!---------------------------------- nc_get_var1_long --------------------------
Interface
 Function nc_get_var1_long(ncid, varid, indexp, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG, C_PTR

 Integer(KIND=C_INT),  VALUE       :: ncid, varid
 Type(C_PTR),          VALUE       :: indexp
 Integer(KIND=C_LONG), Intent(OUT) :: ip

 Integer(KIND=C_INT)               :: nc_get_var1_long

 End Function nc_get_var1_long
End Interface
!---------------------------------- nc_put_var1_float -------------------------
Interface
 Function nc_put_var1_float(ncid, varid, indexp, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT, C_PTR

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: indexp
 Real(KIND=C_FLOAT),  Intent(IN) :: op

 Integer(KIND=C_INT)             :: nc_put_var1_float

 End Function nc_put_var1_float
End Interface
!---------------------------------- nc_get_var1_float -------------------------
Interface
 Function nc_get_var1_float(ncid, varid, indexp, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT, C_PTR

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: indexp
 Real(KIND=C_FLOAT),  Intent(OUT) :: ip

 Integer(KIND=C_INT)              :: nc_get_var1_float

 End Function nc_get_var1_float
End Interface
!---------------------------------- nc_put_var1_double ------------------------
Interface
 Function nc_put_var1_double(ncid, varid, indexp, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_PTR

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: indexp
 Real(KIND=C_DOUBLE), Intent(IN) :: op

 Integer(KIND=C_INT)             :: nc_put_var1_double

 End Function nc_put_var1_double
End Interface
!---------------------------------- nc_get_var1_double ------------------------
Interface
 Function nc_get_var1_double(ncid, varid, indexp, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_PTR

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: indexp
 Real(KIND=C_DOUBLE), Intent(OUT) :: ip

 Integer(KIND=C_INT)              :: nc_get_var1_double

 End Function nc_get_var1_double
End Interface
!---------------------------------- nc_put_var1 ------------------------------
Interface
 Function nc_put_var1(ncid, varid, indexp, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT), VALUE :: ncid, varid
 Type(C_PTR),         VALUE :: indexp
 Type(C_PTR),         VALUE :: op

 Integer(KIND=C_INT)        :: nc_put_var1

 End Function nc_put_var1
End Interface
!---------------------------------- nc_get_var1 ------------------------------
Interface
 Function nc_get_var1(ncid, varid, indexp, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid, varid
 Type(C_PTR),            VALUE         :: indexp
 Character(KIND=C_CHAR), Intent(INOUT) :: op(*) ! op is actually void * in C

 Integer(KIND=C_INT)                   :: nc_get_var1

 End Function nc_get_var1
End Interface
!---------------------------------- nc_put_vara_text --------------------------
Interface
 Function nc_put_vara_text(ncid, varid, startp, countp, op)  BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Type(C_PTR),            VALUE      :: startp, countp
 Character(KIND=C_CHAR), Intent(IN) :: op(*)

 Integer(KIND=C_INT)                :: nc_put_vara_text

 End Function nc_put_vara_text
End Interface
!---------------------------------- nc_get_vara_text --------------------------
Interface
 Function nc_get_vara_text(ncid, varid, startp, countp, ip)  BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Type(C_PTR),            VALUE       :: startp, countp
 Character(KIND=C_CHAR), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                 :: nc_get_vara_text

 End Function nc_get_vara_text
End Interface
!---------------------------------- nc_put_vara_uchar -------------------------
Interface
 Function nc_put_vara_uchar(ncid, varid, startp, countp, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp
 Integer(KIND=CINT1), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_vara_uchar

 End Function nc_put_vara_uchar
End Interface
!---------------------------------- nc_get_vara_uchar -------------------------
Interface
 Function nc_get_vara_uchar(ncid, varid, startp, countp, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp
 Integer(KIND=CINT1), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_vara_uchar

 End Function nc_get_vara_uchar
End Interface
!---------------------------------- nc_put_vara_schar -------------------------
Interface
 Function nc_put_vara_schar(ncid, varid, startp, countp, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp
 Integer(KIND=CINT1), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_vara_schar

 End Function nc_put_vara_schar
End Interface
!---------------------------------- nc_get_vara_schar -------------------------
Interface
 Function nc_get_vara_schar(ncid, varid, startp, countp, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp
 Integer(KIND=CINT1), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_vara_schar

 End Function nc_get_vara_schar
End Interface
!---------------------------------- nc_put_vara_short -------------------------
Interface
 Function nc_put_vara_short(ncid, varid, startp, countp, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT2

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp
 Integer(KIND=CINT2), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_vara_short

 End Function nc_put_vara_short
End Interface
!---------------------------------- nc_get_vara_short -------------------------
Interface
 Function nc_get_vara_short(ncid, varid, startp, countp, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT2

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp
 Integer(KIND=CINT2), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_vara_short

 End Function nc_get_vara_short
End Interface
!--------------------------------- nc_put_vara_int ----------------------------
Interface
 Function nc_put_vara_int(ncid, varid, startp, countp, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp
 Integer(KIND=CINT),  Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_vara_int

 End Function nc_put_vara_int
End Interface
!--------------------------------- nc_get_vara_int ----------------------------
Interface
 Function nc_get_vara_int(ncid, varid, startp, countp, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp
 Integer(KIND=CINT),  Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_vara_int

 End Function nc_get_vara_int
End Interface
!--------------------------------- nc_put_vara_long ---------------------------
Interface
 Function nc_put_vara_long(ncid, varid, startp, countp, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG, C_PTR

 Integer(KIND=C_INT),  VALUE      :: ncid, varid
 Type(C_PTR),          VALUE      :: startp, countp
 Integer(KIND=C_LONG), Intent(IN) :: op(*)

 Integer(KIND=C_INT)              :: nc_put_vara_long

 End Function nc_put_vara_long
End Interface
!--------------------------------- nc_get_vara_long ---------------------------
Interface
 Function nc_get_vara_long(ncid, varid, startp, countp, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG, C_PTR

 Integer(KIND=C_INT),  VALUE       :: ncid, varid
 Type(C_PTR),          VALUE       :: startp, countp
 Integer(KIND=C_LONG), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)               :: nc_get_vara_long

 End Function nc_get_vara_long
End Interface
!--------------------------------- nc_put_vara_float --------------------------
Interface
 Function nc_put_vara_float(ncid, varid, startp, countp, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT, C_PTR

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp
 Real(KIND=C_FLOAT),  Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_vara_float

 End Function nc_put_vara_float
End Interface
!--------------------------------- nc_get_vara_float --------------------------
Interface
 Function nc_get_vara_float(ncid, varid, startp, countp, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT, C_PTR

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp
 Real(KIND=C_FLOAT),  Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_vara_float

 End Function nc_get_vara_float
End Interface
!--------------------------------- nc_put_vara_double -------------------------
Interface
 Function nc_put_vara_double(ncid, varid, startp, countp, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_PTR

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp
 Real(KIND=C_DOUBLE), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_vara_double

 End Function nc_put_vara_double
End Interface
!--------------------------------- nc_get_vara_double -------------------------
Interface
 Function nc_get_vara_double(ncid, varid, startp, countp, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_PTR

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp
 Real(KIND=C_DOUBLE), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_vara_double

 End Function nc_get_vara_double
End Interface
!---------------------------------- nc_put_vara -------------------------------
Interface
 Function nc_put_vara(ncid, varid, startp, countp, op)  BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR

 Integer(KIND=C_INT), VALUE :: ncid, varid
 Type(C_PTR),         VALUE :: startp, countp
 Type(C_PTR),         VALUE :: op

 Integer(KIND=C_INT)        :: nc_put_vara

 End Function nc_put_vara
End Interface
!---------------------------------- nc_get_vara -------------------------------
Interface
 Function nc_get_vara(ncid, varid, startp, countp, ip)  BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid, varid
 Type(C_PTR),            VALUE         :: startp, countp
 Character(KIND=C_CHAR), Intent(INOUT) :: ip(*)

 Integer(KIND=C_INT)                   :: nc_get_vara

 End Function nc_get_vara
End Interface
!--------------------------------- nc_put_vars_text ---------------------------
Interface
 Function nc_put_vars_text(ncid, varid, startp, countp, stridep, op)  BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Type(C_PTR),            VALUE      :: startp, countp, stridep
 Character(KIND=C_CHAR), Intent(IN) :: op(*)

 Integer(KIND=C_INT)                :: nc_put_vars_text

 End Function nc_put_vars_text
End Interface
!--------------------------------- nc_get_vars_text ---------------------------
Interface
 Function nc_get_vars_text(ncid, varid, startp, countp, stridep, ip)  BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Type(C_PTR),            VALUE       :: startp, countp, stridep
 Character(KIND=C_CHAR), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                 :: nc_get_vars_text

 End Function nc_get_vars_text
End Interface
!--------------------------------- nc_put_vars_uchar --------------------------
Interface
 Function nc_put_vars_uchar(ncid, varid, startp, countp, stridep, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp, stridep
 Integer(KIND=CINT1), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_vars_uchar

 End Function nc_put_vars_uchar
End Interface
!--------------------------------- nc_get_vars_uchar --------------------------
Interface
 Function nc_get_vars_uchar(ncid, varid, startp, countp, stridep, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp, stridep
 Integer(KIND=CINT1), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_vars_uchar

 End Function nc_get_vars_uchar
End Interface
!--------------------------------- nc_put_vars_schar --------------------------
Interface
 Function nc_put_vars_schar(ncid, varid, startp, countp, stridep, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp, stridep
 Integer(KIND=CINT1), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_vars_schar

 End Function nc_put_vars_schar
End Interface
!--------------------------------- nc_get_vars_schar --------------------------
Interface
 Function nc_get_vars_schar(ncid, varid, startp, countp, stridep, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp, stridep
 Integer(KIND=CINT1), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_vars_schar

 End Function nc_get_vars_schar
End Interface
!--------------------------------- nc_put_vars_short --------------------------
Interface
 Function nc_put_vars_short(ncid, varid, startp, countp, stridep, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT2

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp, stridep
 Integer(KIND=CINT2), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_vars_short

 End Function nc_put_vars_short
End Interface
!--------------------------------- nc_get_vars_short --------------------------
Interface
 Function nc_get_vars_short(ncid, varid, startp, countp, stridep, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT,  C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT2

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp, stridep
 Integer(KIND=CINT2), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_vars_short

 End Function nc_get_vars_short
End Interface
!--------------------------------- nc_put_vars_int ----------------------------
Interface
 Function nc_put_vars_int(ncid, varid, startp, countp, stridep, op) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp, stridep
 Integer(KIND=CINT),  Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_vars_int

 End Function nc_put_vars_int
End Interface
!--------------------------------- nc_get_vars_int ----------------------------
Interface
 Function nc_get_vars_int(ncid, varid, startp, countp, stridep, ip) BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp, stridep
 Integer(KIND=CINT),  Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_vars_int

 End Function nc_get_vars_int
End Interface
!--------------------------------- nc_put_vars_long ---------------------------
Interface
 Function nc_put_vars_long(ncid, varid, startp, countp, stridep, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG, C_PTR

 Integer(KIND=C_INT),  VALUE      :: ncid, varid
 Type(C_PTR),          VALUE      :: startp, countp, stridep
 Integer(KIND=C_LONG), Intent(IN) :: op(*)

 Integer(KIND=C_INT)              :: nc_put_vars_long

 End Function nc_put_vars_long
End Interface
!--------------------------------- nc_get_vars_long ---------------------------
Interface
 Function nc_get_vars_long(ncid, varid, startp, countp, stridep, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG, C_PTR

 Integer(KIND=C_INT),  VALUE       :: ncid, varid
 Type(C_PTR),          VALUE       :: startp, countp, stridep
 Integer(KIND=C_LONG), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)               :: nc_get_vars_long

 End Function nc_get_vars_long
End Interface
!--------------------------------- nc_put_vars_float --------------------------
Interface
 Function nc_put_vars_float(ncid, varid, startp, countp, stridep, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT, C_PTR

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp, stridep
 Real(KIND=C_FLOAT),  Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_vars_float

 End Function nc_put_vars_float
End Interface
!--------------------------------- nc_get_vars_float --------------------------
Interface
 Function nc_get_vars_float(ncid, varid, startp, countp, stridep, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT, C_PTR

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp, stridep
 Real(KIND=C_FLOAT),  Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_vars_float

 End Function nc_get_vars_float
End Interface
!--------------------------------- nc_put_vars_double -------------------------
Interface
 Function nc_put_vars_double(ncid, varid, startp, countp, stridep, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_PTR

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp, stridep
 Real(KIND=C_DOUBLE), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_vars_double

 End Function nc_put_vars_double
End Interface
!--------------------------------- nc_get_vars_double -------------------------
Interface
 Function nc_get_vars_double(ncid, varid, startp, countp, stridep, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_PTR

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp, stridep
 Real(KIND=C_DOUBLE), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_vars_double

 End Function nc_get_vars_double
End Interface
!--------------------------------- nc_put_vars --------------------------------
Interface
 Function nc_put_vars(ncid, varid, startp, countp, stridep, op)  BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR

 Integer(KIND=C_INT), VALUE :: ncid, varid
 Type(C_PTR),         VALUE :: startp, countp, stridep
 Type(C_PTR),         VALUE :: op

 Integer(KIND=C_INT)        :: nc_put_vars

 End Function nc_put_vars
End Interface
!--------------------------------- nc_get_vars ---------------------------
Interface
 Function nc_get_vars(ncid, varid, startp, countp, stridep, ip)  BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid, varid
 Type(C_PTR),            VALUE         :: startp, countp, stridep
 Character(KIND=C_CHAR), Intent(INOUT) :: ip(*)

 Integer(KIND=C_INT)                   :: nc_get_vars

 End Function nc_get_vars
End Interface
!--------------------------------- nc_put_varm_text ---------------------------
Interface
! array of characters
 Function nc_put_varm_text(ncid, varid, startp, countp, stridep, imapp,op)  &
                             BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Type(C_PTR),            VALUE      :: startp, countp, stridep, imapp
 Character(KIND=C_CHAR), Intent(IN) :: op(*)

 Integer(KIND=C_INT)                :: nc_put_varm_text

 End Function nc_put_varm_text
End Interface
!--------------------------------- nc_get_varm_text ---------------------------
Interface
 Function nc_get_varm_text(ncid, varid, startp, countp, stridep, imapp,ip)  &
                           BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Type(C_PTR),            VALUE       :: startp, countp, stridep, imapp
 Character(KIND=C_CHAR), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                 :: nc_get_varm_text

 End Function nc_get_varm_text
End Interface
!--------------------------------- nc_put_varm_uchar --------------------------
Interface
 Function nc_put_varm_uchar(ncid, varid, startp, countp, stridep, imapp, op)  &
                            BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp, stridep, imapp
 Integer(KIND=CINT1), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_varm_uchar

 End Function nc_put_varm_uchar
End Interface
!--------------------------------- nc_get_varm_uchar --------------------------
Interface
 Function nc_get_varm_uchar(ncid, varid, startp, countp, stridep, imapp, ip)  &
                            BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp, stridep, imapp
 Integer(KIND=CINT1), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_varm_uchar

 End Function nc_get_varm_uchar
End Interface
!--------------------------------- nc_put_varm_schar --------------------------
Interface
 Function nc_put_varm_schar(ncid, varid, startp, countp, stridep, imapp, op)  &
                            BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp, stridep, imapp
 Integer(KIND=CINT1), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_varm_schar

 End Function nc_put_varm_schar
End Interface
!--------------------------------- nc_get_varm_schar --------------------------
Interface
 Function nc_get_varm_schar(ncid, varid, startp, countp, stridep, imapp, ip)  &
                            BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp, stridep, imapp
 Integer(KIND=CINT1), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_varm_schar

 End Function nc_get_varm_schar
End Interface
!--------------------------------- nc_put_varm_short --------------------------
Interface
 Function nc_put_varm_short(ncid, varid, startp, countp, stridep, imapp, op)  &
                            BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT2

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp, stridep, imapp
 Integer(KIND=CINT2), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_varm_short

 End Function nc_put_varm_short
End Interface
!--------------------------------- nc_get_varm_short --------------------------
Interface
 Function nc_get_varm_short(ncid, varid, startp, countp, stridep, imapp, ip)  &
                            BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT2

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp, stridep, imapp
 Integer(KIND=CINT2), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_varm_short

 End Function nc_get_varm_short
End Interface
!--------------------------------- nc_put_varm_int ----------------------------
Interface
 Function nc_put_varm_int(ncid, varid, startp, countp, stridep, imapp, op)  &
                          BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp, stridep, imapp
 Integer(KIND=CINT),  Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_varm_int

 End Function nc_put_varm_int
End Interface
!--------------------------------- nc_get_varm_int ----------------------------
Interface
 Function nc_get_varm_int(ncid, varid, startp, countp, stridep, imapp, ip)  &
                          BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_PTR
 USE NETCDF_NC_DATA, ONLY: CINT

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp, stridep, imapp
 Integer(KIND=CINT),  Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_varm_int

 End Function nc_get_varm_int
End Interface
!--------------------------------- nc_put_varm_long ---------------------------
Interface
 Function nc_put_varm_long(ncid, varid, startp, countp, stridep, imapp, op)  &
                           BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG, C_PTR

 Integer(KIND=C_INT),  VALUE      :: ncid, varid
 Type(C_PTR),          VALUE      :: startp, countp, stridep, imapp
 Integer(KIND=C_LONG), Intent(IN) :: op(*)

 Integer(KIND=C_INT)              :: nc_put_varm_long

 End Function nc_put_varm_long
End Interface
!--------------------------------- nc_get_varm_long ---------------------------
Interface
 Function nc_get_varm_long(ncid, varid, startp, countp, stridep, imapp, ip)  &
                           BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG, C_PTR

 Integer(KIND=C_INT),  VALUE       :: ncid, varid
 Type(C_PTR),          VALUE       :: startp, countp, stridep, imapp
 Integer(KIND=C_LONG), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)               :: nc_get_varm_long

 End Function nc_get_varm_long
End Interface
!--------------------------------- nc_put_varm_float --------------------------
Interface
 Function nc_put_varm_float(ncid, varid, startp, countp, stridep, imapp, op)  &
                            BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT, C_PTR

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp, stridep, imapp
 Real(KIND=C_FLOAT),  Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_varm_float

 End Function nc_put_varm_float
End Interface
!--------------------------------- nc_get_varm_float --------------------------
Interface
 Function nc_get_varm_float(ncid, varid, startp, countp, stridep, imapp, ip)  &
                            BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT, C_PTR

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp, stridep, imapp
 Real(KIND=C_FLOAT),  Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_varm_float

 End Function nc_get_varm_float
End Interface
!--------------------------------- nc_put_varm_double -------------------------
Interface
 Function nc_put_varm_double(ncid, varid, startp, countp, stridep,imapp, op)  &
                             BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_PTR

 Integer(KIND=C_INT), VALUE      :: ncid, varid
 Type(C_PTR),         VALUE      :: startp, countp, stridep, imapp
 Real(KIND=C_DOUBLE), Intent(IN) :: op(*)

 Integer(KIND=C_INT)             :: nc_put_varm_double

 End Function nc_put_varm_double
End Interface
!--------------------------------- nc_get_varm_double -------------------------
Interface
 Function nc_get_varm_double(ncid, varid, startp, countp, stridep,imapp, ip)  &
                             BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_PTR

 Integer(KIND=C_INT), VALUE       :: ncid, varid
 Type(C_PTR),         VALUE       :: startp, countp, stridep, imapp
 Real(KIND=C_DOUBLE), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)              :: nc_get_varm_double

 End Function nc_get_varm_double
End Interface
!--------------------------------- nc_inq_att --------------------------------
Interface
 Function nc_inq_att(ncid, varid, name, xtypep, lenp)  BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Integer(KIND=C_INT),    Intent(OUT) :: xtypep
 Integer(KIND=C_SIZE_T), Intent(OUT) :: lenp

 Integer(KIND=C_INT)                 :: nc_inq_att

 End Function nc_inq_att
End Interface
!--------------------------------- nc_inq_attid ------------------------------
Interface
 Function nc_inq_attid(ncid, varid, name, attnump)   BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Integer(KIND=C_INT),    Intent(OUT) :: attnump

 Integer(KIND=C_INT)                 :: nc_inq_attid

 End Function nc_inq_attid
End Interface
!--------------------------------- nc_inq_atttype ----------------------------
Interface
 Function nc_inq_atttype(ncid, varid, name, xtypep)   BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Integer(KIND=C_INT),    Intent(OUT) :: xtypep

 Integer(KIND=C_INT)                 :: nc_inq_atttype

 End Function nc_inq_atttype
End Interface
!--------------------------------- nc_inq_attlen -----------------------------
Interface
 Function nc_inq_attlen(ncid, varid, name, lenp)   BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Integer(KIND=C_SIZE_T), Intent(OUT) :: lenp

 Integer(KIND=C_INT)                 :: nc_inq_attlen

 End Function nc_inq_attlen
End Interface
!--------------------------------- nc_inq_attname ----------------------------
Interface
 Function nc_inq_attname(ncid, varid, attnum, name)   BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid, varid, attnum
 Character(KIND=C_CHAR), Intent(INOUT) :: name(*)

 Integer(KIND=C_INT)                   :: nc_inq_attname

 End Function nc_inq_attname
End Interface
!--------------------------------- nc_copy_att -------------------------------
Interface
 Function nc_copy_att(ncid_in, varid_in, name, ncid_out, varid_out )  BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid_in, varid_in, varid_out, &
                                       ncid_out
 Character(KIND=C_CHAR), Intent(IN) :: name(*)

 Integer(KIND=C_INT)                :: nc_copy_att

 End Function nc_copy_att
End Interface
!--------------------------------- nc_rename_att -----------------------------
Interface
 Function nc_rename_att(ncid, varid, name, newname)   BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN) :: name(*), newname(*)

 Integer(KIND=C_INT)                :: nc_rename_att

 End Function nc_rename_att
End Interface
!--------------------------------- nc_del_att --------------------------------
Interface
 Function nc_del_att(ncid, varid, name)  BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN) :: name(*)

 Integer(KIND=C_INT)                :: nc_del_att

 End Function nc_del_att
End Interface
!--------------------------------- nc_put_att_text ---------------------------
Interface
 Function nc_put_att_text(ncid, varid, name, nlen, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Integer(KIND=C_SIZE_T), VALUE      :: nlen
 Character(KIND=C_CHAR), Intent(IN) :: name(*)
 Character(KIND=C_CHAR), Intent(IN) :: op(*)

 Integer(KIND=C_INT)                :: nc_put_att_text

 End Function nc_put_att_text
End Interface
!--------------------------------- nc_get_att_text ---------------------------
Interface
 Function nc_get_att_text(ncid, varid, name, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Character(KIND=C_CHAR), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                 :: nc_get_att_text

 End Function nc_get_att_text
End Interface
!--------------------------------- nc_put_att_uchar --------------------------
Interface
 Function nc_put_att_uchar(ncid, varid, name, xtype, nlen, op)   BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_SIZE_T, C_CHAR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Integer(KIND=C_SIZE_T), VALUE      :: nlen
 Integer(KIND=C_INT),    VALUE      :: xtype
 Character(KIND=C_CHAR), Intent(IN) :: name(*)
 Integer(KIND=CINT1),    Intent(IN) :: op(*)

 Integer(KIND=C_INT)                :: nc_put_att_uchar

 End Function nc_put_att_uchar
End Interface
!--------------------------------- nc_get_att_uchar --------------------------
Interface
 Function nc_get_att_uchar(ncid, varid, name, ip)   BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_CHAR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Integer(KIND=CINT1),    Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                 :: nc_get_att_uchar

 End Function nc_get_att_uchar
End Interface
!--------------------------------- nc_put_att_schar --------------------------
Interface
 Function nc_put_att_schar(ncid, varid, name, xtype, nlen, op)   BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_SIZE_T, C_CHAR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Integer(KIND=C_SIZE_T), VALUE      :: nlen
 Integer(KIND=C_INT),    VALUE      :: xtype
 Character(KIND=C_CHAR), Intent(IN) :: name(*)
 Integer(KIND=CINT1),    Intent(IN) :: op(*)

 Integer(KIND=C_INT)                :: nc_put_att_schar

 End Function nc_put_att_schar
End Interface
!--------------------------------- nc_get_att_schar --------------------------
Interface
 Function nc_get_att_schar(ncid, varid, name, ip)   BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_CHAR
 USE NETCDF_NC_DATA, ONLY: CINT1

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Integer(KIND=CINT1),    Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                 :: nc_get_att_schar

 End Function nc_get_att_schar
End Interface
!--------------------------------- nc_put_att_short --------------------------
Interface
 Function nc_put_att_short(ncid, varid, name, xtype, nlen, op)   BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_SIZE_T, C_CHAR
 USE NETCDF_NC_DATA, ONLY: CINT2

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Integer(KIND=C_SIZE_T), VALUE      :: nlen
 Integer(KIND=C_INT),    VALUE      :: xtype
 Character(KIND=C_CHAR), Intent(IN) :: name(*)
 Integer(KIND=CINT2),    Intent(IN) :: op(*)

 Integer(KIND=C_INT)                :: nc_put_att_short

 End Function nc_put_att_short
End Interface
!--------------------------------- nc_get_att_short --------------------------
Interface
 Function nc_get_att_short(ncid, varid, name, ip)   BIND(C)

 USE ISO_C_BINDING,  ONLY: C_INT, C_CHAR
 USE NETCDF_NC_DATA, ONLY: CINT2

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Integer(KIND=CINT2),    Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                 :: nc_get_att_short

 End Function nc_get_att_short
End Interface
!--------------------------------- nc_put_att_int --------------------------
Interface
 Function nc_put_att_int(ncid, varid, name, xtype, nlen, op)   BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Integer(KIND=C_SIZE_T), VALUE      :: nlen
 Integer(KIND=C_INT),    VALUE      :: xtype
 Character(KIND=C_CHAR), Intent(IN) :: name(*)
 Integer(KIND=C_INT),    Intent(IN) :: op(*)

 Integer(KIND=C_INT)                :: nc_put_att_int

 End Function nc_put_att_int
End Interface
!--------------------------------- nc_get_att_int -----------------------------
Interface
 Function nc_get_att_int(ncid, varid, name, ip)   BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Integer(KIND=C_INT),    Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                 :: nc_get_att_int

 End Function nc_get_att_int
End Interface
!--------------------------------- nc_put_att_long --------------------------
Interface
 Function nc_put_att_long(ncid, varid, name, xtype, nlen, op)   BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_LONG, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Integer(KIND=C_SIZE_T), VALUE      :: nlen
 Integer(KIND=C_INT),    VALUE      :: xtype
 Character(KIND=C_CHAR), Intent(IN) :: name(*)
 Integer(KIND=C_LONG),   Intent(IN) :: op(*)

 Integer(KIND=C_INT)                :: nc_put_att_long

 End Function nc_put_att_long
End Interface
!--------------------------------- nc_get_att_long --------------------------
Interface
 Function nc_get_att_long(ncid, varid, name, ip)   BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Integer(KIND=C_LONG),   Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                 :: nc_get_att_long

 End Function nc_get_att_long
End Interface
!--------------------------------- nc_put_att_float --------------------------
Interface
 Function nc_put_att_float(ncid, varid, name, xtype, nlen, op)   BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_FLOAT, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Integer(KIND=C_SIZE_T), VALUE      :: nlen
 Integer(KIND=C_INT),    VALUE      :: xtype
 Character(KIND=C_CHAR), Intent(IN) :: name(*)
 Real(KIND=C_FLOAT),     Intent(IN) :: op(*)

 Integer(KIND=C_INT) :: nc_put_att_float

 End Function nc_put_att_float
End Interface
!--------------------------------- nc_get_att_float --------------------------
Interface
 Function nc_get_att_float(ncid, varid, name, ip)   BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Real(KIND=C_FLOAT),     Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                 :: nc_get_att_float

 End Function nc_get_att_float
End Interface
!--------------------------------- nc_put_att_double -------------------------
Interface
 Function nc_put_att_double(ncid, varid, name, xtype, nlen, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_DOUBLE, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid, varid
 Integer(KIND=C_SIZE_T), VALUE      :: nlen
 Integer(KIND=C_INT),    VALUE      :: xtype
 Character(KIND=C_CHAR), Intent(IN) :: name(*)
 Real(KIND=C_DOUBLE),    Intent(IN) :: op(*)

 Integer(KIND=C_INT)                :: nc_put_att_double

 End Function nc_put_att_double
End Interface
!------------------------------- nc_get_att_double -------------------------
Interface
 Function nc_get_att_double(ncid, varid, name, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Real(KIND=C_DOUBLE),    Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                 :: nc_get_att_double

 End Function nc_get_att_double
End Interface
!------------------------------- nc_copy_var --------------------------------
Interface
 Function nc_copy_var(ncid_in, varid, ncid_out)  BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: ncid_in, varid, ncid_out

 Integer(KIND=C_INT)        :: nc_copy_var

 End Function nc_copy_var
End Interface
!------------------------------- nc_set_default_format -----------------------
Interface
 Function nc_set_default_format(newform, old_format)  BIND(C)
!
 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: newform
 Integer(KIND=C_INT), Intent(OUT) :: old_format

 Integer(KIND=C_INT)              :: nc_set_default_format

 End Function nc_set_default_format
End Interface
!---------------------------- Start of module procedures ---------------------
CONTAINS

! Utilities to support C interface routines

!----------------------------------- addCNullChar -----------------------------
 Function addCNullChar(string, nlen) Result(cstring)

! Add a C_NULL_CHAR to a string to create a C compatible
! string. Assumes target variable will be of length
! LEN(string)+1. Trailing blanks will be stripped
! from string and length of trimmed string will
! be returned in nlen

! USE ISO_C_BINDING

 Implicit NONE

 Character(LEN=*), Intent(IN)    :: string
 Integer,          Intent(INOUT) :: nlen

 Character(LEN=(LEN(string)+1))  :: cstring

 Integer :: inull


! first check to see if we already have a C NULL char attached
! to string and strip trailing blanks. We will overwrite it if
! we do

 nlen  = LEN_TRIM(string)
 inull = SCAN(string, C_NULL_CHAR)

 If (inull > 1) nlen = inull - 1
 nlen = MAX(1,nlen)  ! Make sure nlen is at least 1

! append null char to trimmed string

 cstring = REPEAT(" ", LEN(cstring)) ! init to blanks
 cstring = string(1:nlen)//C_NULL_CHAR

 End Function addCNullChar
!----------------------------------- stripCNullChar ----------------------------
 Function StripCNullChar(cstring, nlen) Result(string)

! Check cstring for a C NULL char, strip it off and
! return regular string. Limit length of cstring loaded
! into string to nlen

! USE ISO_C_BINDING, ONLY: C_NULL_CHAR

 Implicit NONE

 Character(LEN=*), Intent(IN) :: cstring
 Integer,          Intent(IN) :: nlen

 Character(LEN=nlen)          :: string

 Integer :: ie, inull

 ie    = LEN_TRIM(cstring)
 inull = SCAN(cstring, C_NULL_CHAR)

 If (inull > 1) ie=inull-1
 ie = MAX(1, MIN(ie,nlen)) ! limit ie to 1 or nlen

 string       = REPEAT(" ", nlen)
 string(1:ie) = cstring(1:ie)

 End Function StripCNullChar
!
!----------------------------End of Module netcdf_c_interfaces ----------------
End Module netcdf_nc_interfaces
