Module netcdf4_nc_interfaces

! Fortran interfaces to netCDF4 C functions using FORTRAN 2003 C 
! Interoperability features. These interfaces are for the base
! netCDF C routines in the libsrc4 directory and the results from
! running CPP on the fort-xxx.c routines to get the cfortran.h
! generated interfaces

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

! Version 1.: June.  2007 - Initial version - split from ntecdf_nc_interfaces 
! Version 2.: April, 2009 - Interfaces based on netcdf-4.0.1 source         
! Version 3.: April, 2010 - Interfaces based on netcdf-4.1.1 source         
! Version 4.: Aug,   2013 - Added nc_rename_grp interfaces for netcdf-C 4.3.1

 USE netcdf_nc_interfaces

 Implicit NONE

!--------- Define default C interface parameters from netcdf.h ---------------

! Begin explicit interfaces for nc_ functions. Note that some interfaces
! expect data to be passed as C_PTR type variables. These data are arrays
! that could have a NULL pointer passed instead of the array. Also, note
! the use of generic interfaces to support routines that handle text data
! to allow text to be passed as either a single string or an array of
! single characters

! Also note that each interface has an explicit USE ISO_C_BINDING. A better
! solution is to use the F2003 IMPORT statement (I originally had it this way)
! However its best to leave the interfaces as is for now because there might
! be a few compilers out there that support most of the C-interop facility but
! for some reason haven't implemented IMPORT yet. 

! NETCDF 4 functions supported by FORTRAN interface

!----------------------------- nc_create_par_fortran -------------------------
Interface
 Function nc_create_par_fortran(path, cmode, comm, info, ncidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Character(KIND=C_CHAR), Intent(IN)  :: path(*)
 Integer(KIND=C_INT),    VALUE       :: cmode, comm, info
 Integer(KIND=C_INT),    Intent(OUT) :: ncidp

 Integer(KIND=C_INT)                 :: nc_create_par_fortran

 End Function nc_create_par_fortran
End Interface
!----------------------------- nc_open_par_fortran ---------------------------
Interface
 Function nc_open_par_fortran(path, mode, comm, info, ncidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Character(KIND=C_CHAR), Intent(IN)  :: path(*)
 Integer(KIND=C_INT),    VALUE       :: mode, comm, info
 Integer(KIND=C_INT),    Intent(OUT) :: ncidp

 Integer(KIND=C_INT)                 :: nc_open_par_fortran

 End Function nc_open_par_fortran
End Interface
!------------------------------- nc_var_par_access ----------------------------
Interface
 Function nc_var_par_access(ncid, varid, par_access) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: ncid, varid, par_access

 Integer(KIND=C_INT)        :: nc_var_par_access

 End Function nc_var_par_access
End Interface
!------------------------------- nc_inq_ncid ----------------------------------
Interface
 Function nc_inq_ncid(ncid, name, grp_ncid) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Character(KIND=C_CHAR), Intent(IN)    :: name(*)
 Integer(KIND=C_INT),    Intent(INOUT) :: grp_ncid

 Integer(KIND=C_INT)                   :: nc_inq_ncid

 End Function nc_inq_ncid
End Interface
!------------------------------- nc_inq_grps ----------------------------------
Interface
 Function nc_inq_grps(ncid, numgrps, ncids) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid
 Integer(KIND=C_INT), Intent(INOUT) :: numgrps
 Integer(KIND=C_INT), Intent(INOUT) :: ncids(*)

 Integer(KIND=C_INT)                :: nc_inq_grps

 End Function nc_inq_grps
End Interface
!------------------------------- nc_inq_grpname -------------------------------
Interface
 Function nc_inq_grpname(ncid, name) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Character(KIND=C_CHAR), Intent(INOUT) :: name(*)

 Integer(KIND=C_INT)                   :: nc_inq_grpname

 End Function nc_inq_grpname
End Interface
!------------------------------- nc_inq_grpname_full --------------------------
Interface
 Function nc_inq_grpname_full(ncid, nlen, name) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: nlen
 Character(KIND=C_CHAR), Intent(INOUT) :: name(*)

 Integer(KIND=C_INT)                   :: nc_inq_grpname_full

 End Function nc_inq_grpname_full
End Interface
!------------------------------- nc_inq_grpname_len ---------------------------
Interface
 Function nc_inq_grpname_len(ncid, nlen) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: nlen

 Integer(KIND=C_INT)                   :: nc_inq_grpname_len

 End Function nc_inq_grpname_len
End Interface
!------------------------------- nc_inq_grp_full_ncid -------------------------
Interface
 Function nc_inq_grp_full_ncid(ncid, full_name, grp_ncid) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_INT),    Intent(INOUT) :: grp_ncid
 Character(KIND=C_CHAR), Intent(INOUT) :: full_name(*)

 Integer(KIND=C_INT)                   :: nc_inq_grp_full_ncid

 End Function nc_inq_grp_full_ncid
End Interface
!------------------------------- nc_inq_grp_parent ----------------------------
Interface
 Function nc_inq_grp_parent(ncid, parent_ncid) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid
 Integer(KIND=C_INT), Intent(INOUT) :: parent_ncid

 Integer(KIND=C_INT)                :: nc_inq_grp_parent

 End Function nc_inq_grp_parent
End Interface
!------------------------------- nc_inq_grp_ncid ------------------------------
Interface
 Function nc_inq_grp_ncid(ncid, grp_name, grp_ncid) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Character(KIND=C_CHAR), Intent(IN)    :: grp_name(*)
 Integer(KIND=C_INT),    Intent(INOUT) :: grp_ncid

 Integer(KIND=C_INT)                   :: nc_inq_grp_ncid

 End Function nc_inq_grp_ncid
End Interface
!------------------------------- nc_inq_varids_f ------------------------------
Interface
 Function nc_inq_varids_f(ncid, nvars, varids) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid
 Integer(KIND=C_INT), Intent(INOUT) :: nvars
 Integer(KIND=C_INT), Intent(INOUT) :: varids(*)

 Integer(KIND=C_INT)                :: nc_inq_varids_f

 End Function nc_inq_varids_f
End Interface
!------------------------------- nc_inq_dimids_f ------------------------------
Interface
 Function nc_inq_dimids_f(ncid, ndims, dimids, parent) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid, parent
 Integer(KIND=C_INT), Intent(INOUT) :: ndims
 Integer(KIND=C_INT), Intent(INOUT) :: dimids(*)

 Integer(KIND=C_INT)                :: nc_inq_dimids_f

 End Function nc_inq_dimids_f
End Interface
!------------------------------- nc_inq_typeids -------------------------------
Interface
 Function nc_inq_typeids(ncid, ntypes, typeids) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid
 Integer(KIND=C_INT), Intent(INOUT) :: ntypes
 Integer(KIND=C_INT), Intent(INOUT) :: typeids(*)

 Integer(KIND=C_INT)                :: nc_inq_typeids

 End Function nc_inq_typeids
End Interface
!------------------------------- nc_inq_typeid --------------------------------
Interface
 Function nc_inq_typeid(ncid, name, typeid) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT), VALUE         :: ncid
 Character(KIND=C_CHAR), Intent(IN) :: name(*)
 Integer(KIND=C_INT), Intent(INOUT) :: typeid

 Integer(KIND=C_INT)                :: nc_inq_typeid

 End Function nc_inq_typeid
End Interface
!------------------------------- nc_def_grp -----------------------------------
Interface
 Function nc_def_grp(parent_ncid, name, new_ncid) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: parent_ncid
 Character(KIND=C_CHAR), Intent(IN)    :: name(*)
 Integer(KIND=C_INT),    Intent(INOUT) :: new_ncid

 Integer(KIND=C_INT)                   :: nc_def_grp

 End Function nc_def_grp
End Interface
!------------------------------- nc_rename_grp --------------------------------
Interface
 Function nc_rename_grp(grpid, name) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: grpid
 Character(KIND=C_CHAR), Intent(IN)    :: name(*)

 Integer(KIND=C_INT)                   :: nc_rename_grp

 End Function nc_rename_grp
End Interface
!------------------------------- nc_def_compound ------------------------------
Interface
 Function nc_def_compound(ncid, isize, name, typeidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_SIZE_T), VALUE         :: isize
 Character(KIND=C_CHAR), Intent(IN)    :: name(*)
 Integer(KIND=C_INT),    Intent(INOUT) :: typeidp

 Integer(KIND=C_INT)                   :: nc_def_compound

 End Function nc_def_compound
End Interface
!------------------------------- nc_insert_compound ---------------------------
Interface
 Function nc_insert_compound(ncid, xtype, name, offset, field_typeid) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE      :: ncid
 Integer(KIND=C_INT),    VALUE      :: xtype, field_typeid ! nc_type in C 
 Integer(KIND=C_SIZE_T), VALUE      :: offset
 Character(KIND=C_CHAR), Intent(IN) :: name(*)

 Integer(KIND=C_INT)                :: nc_insert_compound

 End Function nc_insert_compound
End Interface
!------------------------------- nc_insert_array_compound_f -------------------
Interface
 Function nc_insert_array_compound_f(ncid, xtype, name, offset, field_typeid, &
                                     ndims, dim_sizes) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid, ndims
 Integer(KIND=C_INT),    VALUE         :: xtype, field_typeid  ! nc_type in C
 Integer(KIND=C_SIZE_T), VALUE         :: offset
 Character(KIND=C_CHAR), Intent(IN)    :: name(*)
 Integer(KIND=C_INT),    Intent(INOUT) :: dim_sizes(*)

 Integer(KIND=C_INT)                   :: nc_insert_array_compound_f

 End Function nc_insert_array_compound_f
End Interface
!------------------------------- nc_inq_type ----------------------------------
Interface
 Function nc_inq_type(ncid, xtype, name, isize) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_INT),    VALUE         :: xtype ! nc_type in C
 Character(KIND=C_CHAR), Intent(IN)    :: name(*)
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: isize

 Integer(KIND=C_INT)                   :: nc_inq_type

 End Function nc_inq_type
End Interface
!------------------------------- nc_inq_compound -----------------------------
Interface
 Function nc_inq_compound(ncid, xtype, name, isize, nfieldsp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_INT),    VALUE         :: xtype ! nc_type in C
 Character(KIND=C_CHAR), Intent(INOUT) :: name(*)
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: isize, nfieldsp

 Integer(KIND=C_INT)                   :: nc_inq_compound

 End Function nc_inq_compound
End Interface
!------------------------------- nc_inq_compound_name -------------------------
Interface
 Function nc_inq_compound_name(ncid, xtype, name) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_INT),    VALUE         :: xtype ! nc_type in C
 Character(KIND=C_CHAR), Intent(INOUT) :: name(*)

 Integer(KIND=C_INT)                   :: nc_inq_compound_name

 End Function nc_inq_compound_name
End Interface
!------------------------------- nc_inq_compound_size -------------------------
Interface
 Function nc_inq_compound_size(ncid, xtype, isize) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_INT),    VALUE         :: xtype ! nc_type in C
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: isize

 Integer(KIND=C_INT)                   :: nc_inq_compound_size

 End Function nc_inq_compound_size
End Interface
!------------------------------- nc_inq_compound_nfields ----------------------
Interface
 Function nc_inq_compound_nfields(ncid, xtype, nfieldsp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_INT),    VALUE         :: xtype ! nc_type in C
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: nfieldsp

 Integer(KIND=C_INT)                   :: nc_inq_compound_nfields

 End Function nc_inq_compound_nfields
End Interface
!------------------------------- nc_inq_compound_field_f ----------------------
Interface
 Function nc_inq_compound_field_f(ncid, xtype, fieldid, name, offsetp, &
                                  field_typeidp, ndimsp, dim_sizesp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid, fieldid
 Integer(KIND=C_INT),    VALUE         :: xtype  ! nc_type in C
 Integer(KIND=C_INT),    Intent(INOUT) :: field_typeidp  ! nc_type in C
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: offsetp
 Character(KIND=C_CHAR), Intent(INOUT) :: name(*)
 Integer(KIND=C_INT),    Intent(INOUT) :: ndimsp
 Integer(KIND=C_INT),    Intent(INOUT) :: dim_sizesp(*)

 Integer(KIND=C_INT)                   :: nc_inq_compound_field_f

 End Function nc_inq_compound_field_f
End Interface
!------------------------------- nc_inq_compound_fieldoffset ------------------
Interface
 Function nc_inq_compound_fieldoffset(ncid, xtype, fieldid, offsetp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T

 Integer(KIND=C_INT),    VALUE         :: ncid, fieldid
 Integer(KIND=C_INT),    VALUE         :: xtype ! nc_type in C
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: offsetp 

 Integer(KIND=C_INT)                   :: nc_inq_compound_fieldoffset

 End Function nc_inq_compound_fieldoffset
End Interface
!------------------------------- nc_inq_compound_fieldname --------------------
Interface
 Function nc_inq_compound_fieldname(ncid, xtype, fieldid, name) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid, fieldid
 Integer(KIND=C_INT),    VALUE         :: xtype ! nc_type in C
 Character(KIND=C_CHAR), Intent(INOUT) :: name(*)

 Integer(KIND=C_INT)                   :: nc_inq_compound_fieldname

 End Function nc_inq_compound_fieldname
End Interface
!------------------------------- nc_inq_compound_fieldindex -------------------
Interface
 Function nc_inq_compound_fieldindex(ncid, xtype, name, fieldidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_INT),    VALUE         :: xtype ! nc_type in C
 Character(KIND=C_CHAR), Intent(IN)    :: name(*)
 Integer(KIND=C_INT),    Intent(INOUT) :: fieldidp

 Integer(KIND=C_INT)                   :: nc_inq_compound_fieldindex

 End Function nc_inq_compound_fieldindex
End Interface
!------------------------------- nc_inq_compound_fieldtype --------------------
Interface
 Function nc_inq_compound_fieldtype(ncid, xtype, fieldid, field_typeidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid, fieldid
 Integer(KIND=C_INT), VALUE         :: xtype ! nc_type in C
 Integer(KIND=C_INT), Intent(INOUT) :: field_typeidp ! nc_type in C

 Integer(KIND=C_INT)                :: nc_inq_compound_fieldtype

 End Function nc_inq_compound_fieldtype
End Interface
!------------------------------- nc_inq_compound_fieldndims -------------------
Interface
 Function nc_inq_compound_fieldndims(ncid, xtype, fieldid, ndimsp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid, fieldid
 Integer(KIND=C_INT), VALUE         :: xtype ! nc_type in C
 Integer(KIND=C_INT), Intent(INOUT) :: ndimsp

 Integer(KIND=C_INT)                :: nc_inq_compound_fieldndims

 End Function nc_inq_compound_fieldndims
End Interface
!------------------------------- nc_inq_compound_fielddim_sizes ---------------
Interface
 Function nc_inq_compound_fielddim_sizes(ncid, xtype, fieldid, dim_sizes) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid, fieldid
 Integer(KIND=C_INT), VALUE         :: xtype ! nc_type in C
 Integer(KIND=C_INT), Intent(INOUT) :: dim_sizes(*)

 Integer(KIND=C_INT)                :: nc_inq_compound_fielddim_sizes

 End Function nc_inq_compound_fielddim_sizes
End Interface
!------------------------------- nc_def_vlen ----------------------------------
Interface
 Function nc_def_vlen(ncid, name, base_typeid, xtypep) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_INT),    VALUE         :: base_typeid ! nc_type in C
 Integer(KIND=C_INT),    Intent(INOUT) :: xtypep ! nc_type in C 
 Character(KIND=C_CHAR), Intent(IN)    :: name(*)

 Integer(KIND=C_INT)                   :: nc_def_vlen

 End Function nc_def_vlen
End Interface
!------------------------------- nc_inq_vlen ----------------------------------
Interface
 Function nc_inq_vlen(ncid, xtype, name, datum_sizep, base_nc_typep) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_INT),    VALUE         :: xtype ! nc_type in C
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: datum_sizep
 Integer(KIND=C_INT),    Intent(INOUT) :: base_nc_typep ! nc_type in C 
 Character(KIND=C_CHAR), Intent(INOUT) :: name(*)

 Integer(KIND=C_INT)                   :: nc_inq_vlen

 End Function nc_inq_vlen
End Interface
!------------------------------- nc_inq_user_type -----------------------------
Interface
 Function nc_inq_user_type(ncid, xtype, name, isize, base_nc_typep, &
                           nfieldsp, classp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_INT),    VALUE         :: xtype ! nc_type in C
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: isize , nfieldsp
 Integer(KIND=C_INT),    Intent(INOUT) :: base_nc_typep ! nc_type in C 
 Integer(KIND=C_INT),    Intent(INOUT) :: classp
 Character(KIND=C_CHAR), Intent(INOUT) :: name(*)

 Integer(KIND=C_INT)                   :: nc_inq_user_type

 End Function nc_inq_user_type
End Interface
!------------------------------- nc_def_enum ----------------------------------
Interface
 Function nc_def_enum(ncid, base_typeid, name, typeidp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid
 Integer(KIND=C_INT),    VALUE       :: base_typeid ! nc_type in C
 Integer(KIND=C_INT),    Intent(OUT) :: typeidp ! nc_type in C 
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)

 Integer(KIND=C_INT)                 :: nc_def_enum

 End Function nc_def_enum
End Interface
!------------------------------- nc_insert_enum -------------------------------
Interface
 Function nc_insert_enum(ncid, xtype, name, values) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_PTR

 Integer(KIND=C_INT),    VALUE      :: ncid
 Integer(KIND=C_INT),    VALUE      :: xtype ! nc_type in C
 Type(C_PTR),            VALUE      :: values  ! void pointer in C
 Character(KIND=C_CHAR), Intent(IN) :: name(*)

 Integer(KIND=C_INT)                :: nc_insert_enum

 End Function nc_insert_enum
End Interface
!------------------------------- nc_inq_enum ----------------------------------
Interface
 Function nc_inq_enum(ncid, xtype, name, base_nc_typep, base_sizep, &
                      num_membersp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid
 Integer(KIND=C_INT),    VALUE         :: xtype ! nc_type in C
 Integer(KIND=C_INT),    Intent(INOUT) :: base_nc_typep ! nc_type in C 
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: base_sizep, num_membersp 
 Character(KIND=C_CHAR), Intent(INOUT) :: name(*)

 Integer(KIND=C_INT)                   :: nc_inq_enum

 End Function nc_inq_enum
End Interface
!------------------------------- nc_inq_enum_member ---------------------------
Interface
 Function nc_inq_enum_member(ncid, xtype, idx, name, value) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid, idx
 Integer(KIND=C_INT),    VALUE         :: xtype ! nc_type in C
 Character(KIND=C_CHAR), Intent(OUT)   :: value(*)
 Character(KIND=C_CHAR), Intent(INOUT) :: name(*)

 Integer(KIND=C_INT)                   :: nc_inq_enum_member

 End Function nc_inq_enum_member
End Interface
!------------------------------- nc_inq_enum_ident ----------------------------
Interface
 Function nc_inq_enum_ident(ncid, xtype, val, name) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG_LONG, C_CHAR

 Integer(KIND=C_INT),       VALUE         :: ncid
 Integer(KIND=C_INT),       VALUE         :: xtype ! nc_type in C
 Integer(KIND=C_LONG_LONG), VALUE         :: val 
 Character(KIND=C_CHAR),    Intent(INOUT) :: name(*)

 Integer(KIND=C_INT)                      :: nc_inq_enum_ident

 End Function nc_inq_enum_ident
End Interface
!------------------------------- nc_def_opaque --------------------------------
Interface
 Function nc_def_opaque(ncid, isize, name, xtypep) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid
 Integer(KIND=C_SIZE_T), VALUE       :: isize
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Integer(KIND=C_INT),    Intent(OUT) :: xtypep ! nc_type in C

 Integer(KIND=C_INT)                 :: nc_def_opaque

 End Function nc_def_opaque
End Interface
!------------------------------- nc_inq_opaque --------------------------------
Interface
 Function nc_inq_opaque(ncid, xtype, name, sizep) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid
 Integer(KIND=C_INT),    VALUE       :: xtype ! nc_type in C
 Integer(KIND=C_SIZE_T), Intent(OUT) :: sizep
 Character(KIND=C_CHAR), Intent(OUT) :: name(*)

 Integer(KIND=C_INT)                 :: nc_inq_opaque

 End Function nc_inq_opaque
End Interface
!------------------------------- nc_def_var_fill ------------------------------
Interface
 Function nc_def_var_fill(ncid, varid, no_fill, cfill_value_p) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR

 Integer(KIND=C_INT), VALUE :: ncid, varid, no_fill
 Type(C_PTR),         VALUE :: cfill_value_p

 Integer(KIND=C_INT)        :: nc_def_var_fill

 End Function nc_def_var_fill
End Interface
!------------------------------- nc_inq_var_fill ------------------------------
Interface
 Function nc_inq_var_fill(ncid, varid, no_fill, fill_value) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid, varid
 Integer(KIND=C_INT),    Intent(INOUT) :: no_fill
 Character(KIND=C_CHAR), Intent(INOUT) :: fill_value(*)

 Integer(KIND=C_INT)                   :: nc_inq_var_fill

 End Function nc_inq_var_fill
End Interface
!------------------------------- nc_inq_var_szip ------------------------------
Interface
 Function nc_inq_var_szip(ncid, varid, options_mask, pixels_per_block) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid, varid
 Integer(KIND=C_INT), Intent(INOUT) :: options_mask, pixels_per_block 

 Integer(KIND=C_INT)                :: nc_inq_var_szip

 End Function nc_inq_var_szip
End Interface
!------------------------------- nc_def_var_fletcher32 ------------------------
Interface
 Function nc_def_var_fletcher32(ncid, varid, fletcher32) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: ncid, varid, fletcher32

 Integer(KIND=C_INT)        :: nc_def_var_fletcher32

 End Function nc_def_var_fletcher32
End Interface
!------------------------------- nc_inq_var_fletcher32 ------------------------
Interface
 Function nc_inq_var_fletcher32(ncid, varid, fletcher32) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid, varid
 Integer(KIND=C_INT), Intent(INOUT) :: fletcher32

 Integer(KIND=C_INT)                :: nc_inq_var_fletcher32

 End Function nc_inq_var_fletcher32
End Interface
!------------------------------- nc_def_var_deflate ---------------------------
Interface
 Function nc_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: ncid, varid, shuffle, deflate, deflate_level

 Integer(KIND=C_INT)        :: nc_def_var_deflate

 End Function nc_def_var_deflate
End Interface
!------------------------------- nc_inq_var_deflate ---------------------------
Interface
 Function nc_inq_var_deflate(ncid, varid, shuffle, deflate, deflate_level) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid, varid
 Integer(KIND=C_INT), Intent(INOUT) :: shuffle, deflate, deflate_level

 Integer(KIND=C_INT)                :: nc_inq_var_deflate

 End Function nc_inq_var_deflate
End Interface
!------------------------------- nc_def_var_chunking --------------------------
Interface
 Function nc_def_var_chunking(ncid, varid, contiguousp, chunksizesp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T

 Integer(KIND=C_INT),    VALUE         :: ncid, varid, contiguousp
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: chunksizesp

 Integer(KIND=C_INT)                   :: nc_def_var_chunking

 End Function nc_def_var_chunking
End Interface
!------------------------------- nc_inq_var_chunking --------------------------
Interface
 Function nc_inq_var_chunking(ncid, varid, contiguousp, chunksizesp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T

 Integer(KIND=C_INT),    VALUE         :: ncid, varid
 Integer(KIND=C_INT),    Intent(INOUT) :: contiguousp
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: chunksizesp(*)

 Integer(KIND=C_INT)                   :: nc_inq_var_chunking

 End Function nc_inq_var_chunking
End Interface
!------------------------------- nc_def_var_chunking_ints ---------------------
Interface
 Function nc_def_var_chunking_ints(ncid, varid, contiguousp, chunksizesp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR

 Integer(KIND=C_INT), VALUE :: ncid, varid, contiguousp
 Type(C_PTR),         VALUE :: chunksizesp

 Integer(KIND=C_INT)        :: nc_def_var_chunking_ints

 End Function nc_def_var_chunking_ints
End Interface
!------------------------------- nc_inq_var_chunking_ints ---------------------
Interface
 Function nc_inq_var_chunking_ints(ncid, varid, contiguousp, chunksizesp) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid, varid
 Integer(KIND=C_INT), Intent(INOUT) :: contiguousp
 Integer(KIND=C_INT), Intent(INOUT) :: chunksizesp(*)

 Integer(KIND=C_INT)                :: nc_inq_var_chunking_ints

 End Function nc_inq_var_chunking_ints
End Interface
!------------------------------- nc_def_var_endian --------------------------
Interface
 Function nc_def_var_endian(ncid, varid, endiann) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: ncid, varid, endiann

 Integer(KIND=C_INT)        :: nc_def_var_endian

 End Function nc_def_var_endian
End Interface
!------------------------------- nc_inq_var_endian --------------------------
Interface
 Function nc_inq_var_endian(ncid, varid, endiann) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid, varid
 Integer(KIND=C_INT), Intent(INOUT) :: endiann

 Integer(KIND=C_INT)                :: nc_inq_var_endian

 End Function nc_inq_var_endian
End Interface
!------------------------------- nc_put_att -----------------------------------
Interface
 Function nc_put_att(ncid, varid, name, xtype, nlen, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_SIZE_T, C_PTR

 Integer(KIND=C_INT),    VALUE      :: ncid, varid, xtype
 Integer(KIND=C_SIZE_T), VALUE      :: nlen
 Character(KIND=C_CHAR), Intent(IN) :: name(*)
 Type(C_PTR),            VALUE      :: op

 Integer(KIND=C_INT)                :: nc_put_att

 End Function nc_put_att
End Interface
!------------------------------- nc_get_att -----------------------------------
Interface
 Function nc_get_att(ncid, varid, name, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, varid
 Character(KIND=C_CHAR), Intent(IN)  :: name(*)
 Character(KIND=C_CHAR), Intent(OUT) :: op(*)

 Integer(KIND=C_INT)                 :: nc_get_att

 End Function nc_get_att
End Interface
!------------------------------- nc_put_vlen_element --------------------------
Interface
 Function nc_put_vlen_element(ncid, xtype, vlen_element, nlen, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid, xtype 
 Integer(KIND=C_SIZE_T), VALUE         :: nlen 
 Character(KIND=C_CHAR), INTENT(INOUT) :: vlen_element(*)
 Type(C_PTR),            VALUE         :: op

 Integer(KIND=C_INT)                   :: nc_put_vlen_element

 End Function nc_put_vlen_element
End Interface
!------------------------------- nc_get_vlen_element --------------------------
Interface
 Function nc_get_vlen_element(ncid, xtype, vlen_element, nlen, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_SIZE_T, C_PTR

 Integer(KIND=C_INT),    VALUE         :: ncid, xtype 
 Integer(KIND=C_SIZE_T), Intent(INOUT) :: nlen 
 Character(KIND=C_CHAR), INTENT(INOUT) :: vlen_element(*)
 Character(KIND=C_CHAR), Intent(INOUT) :: op(*)

 Integer(KIND=C_INT)                   :: nc_get_vlen_element

 End Function nc_get_vlen_element
End Interface
!------------------------------- nc_free_vlen ---------------------------------
Interface
 Function nc_free_vlen(vl) BIND(C)

 USE ISO_C_BINDING, ONLY: C_PTR, C_INT

 Type(C_PTR), VALUE  :: vl

 Integer(KIND=C_INT) :: nc_free_vlen

 End Function nc_free_vlen
End Interface
!------------------------------- nc_free_vlens -------------------------------
Interface
 Function nc_free_vlens(len, vl) BIND(C)

 USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_SIZE_T
  
 Integer(C_SIZE_T), Intent(IN) :: len
 Type(C_PTR),       VALUE      :: vl

 Integer(KIND=C_INT)           :: nc_free_vlens

 End Function nc_free_vlens
End Interface
!------------------------------- nc_free_string ------------------------------
Interface
 Function nc_free_string(len, vl) BIND(C)

 USE ISO_C_BINDING, ONLY: C_PTR, C_INT, C_SIZE_T
  
 Integer(C_SIZE_T), Intent(IN) :: len
 Type(C_PTR),       VALUE      :: vl

 Integer(KIND=C_INT)           :: nc_free_string

 End Function nc_free_string
End Interface
!------------------------------- nc_put_var1_longlong -------------------------
Interface
 Function nc_put_var1_longlong(ncid, varid, indexp, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG_LONG, C_PTR

 Integer(KIND=C_INT),       VALUE      :: ncid, varid
 Type(C_PTR),               VALUE      :: indexp
 Integer(KIND=C_LONG_LONG), Intent(IN) :: op

 Integer(KIND=C_INT)                   :: nc_put_var1_longlong

 End Function nc_put_var1_longlong
End Interface
!------------------------------- nc_get_var1_longlong -------------------------
Interface
 Function nc_get_var1_longlong(ncid, varid, indexp, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG_LONG, C_PTR

 Integer(KIND=C_INT),       VALUE       :: ncid, varid
 Type(C_PTR),               VALUE       :: indexp
 Integer(KIND=C_LONG_LONG), Intent(OUT) :: ip

 Integer(KIND=C_INT)                    :: nc_get_var1_longlong

 End Function nc_get_var1_longlong
End Interface
!--------------------------------- nc_put_vara_longlong -----------------------
Interface
 Function nc_put_vara_longlong(ncid, varid, startp, countp, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG_LONG, C_PTR

 Integer(KIND=C_INT),       VALUE      :: ncid, varid
 Type(C_PTR),               VALUE      :: startp, countp
 Integer(KIND=C_LONG_LONG), Intent(IN) :: op(*)

 Integer(KIND=C_INT)                   :: nc_put_vara_longlong

 End Function nc_put_vara_longlong
End Interface
!--------------------------------- nc_get_vara_longlong -----------------------
Interface
 Function nc_get_vara_longlong(ncid, varid, startp, countp, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG_LONG, C_PTR

 Integer(KIND=C_INT),       VALUE       :: ncid, varid
 Type(C_PTR),               VALUE       :: startp, countp
 Integer(KIND=C_LONG_LONG), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                    :: nc_get_vara_longlong

 End Function nc_get_vara_longlong
End Interface
!--------------------------------- nc_put_varm_longlong ----------------------
Interface
 Function nc_put_varm_longlong(ncid, varid, startp, countp, stridep, imapp, &
                               op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_LONG_LONG

 Integer(KIND=C_INT),       VALUE      :: ncid, varid
 Type(C_PTR),               VALUE      :: startp, countp, stridep, imapp
 Integer(KIND=C_LONG_LONG), Intent(IN) :: op(*)

 Integer(KIND=C_INT)                   :: nc_put_varm_longlong

 End Function nc_put_varm_longlong
End Interface
!--------------------------------- nc_get_varm_longlong -----------------------
Interface
 Function nc_get_varm_longlong(ncid, varid, startp, countp, stridep, imapp, &
                               ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_LONG_LONG

 Integer(KIND=C_INT),       VALUE       :: ncid, varid
 Type(C_PTR),               VALUE       :: startp, countp, stridep, imapp
 Integer(KIND=C_LONG_LONG), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                    :: nc_get_varm_longlong

 End Function nc_get_varm_longlong
End Interface
!--------------------------------- nc_put_vars_longlong -----------------------
Interface
 Function nc_put_vars_longlong(ncid, varid, startp, countp, stridep, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_LONG_LONG

 Integer(KIND=C_INT),       VALUE      :: ncid, varid
 Type(C_PTR),               VALUE      :: startp, countp, stridep
 Integer(KIND=C_LONG_LONG), Intent(IN) :: op(*)

 Integer(KIND=C_INT)                   :: nc_put_vars_longlong

 End Function nc_put_vars_longlong
End Interface
!--------------------------------- nc_get_vars_longlong -----------------------
Interface
 Function nc_get_vars_longlong(ncid, varid, startp, countp, stridep, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_LONG_LONG

 Integer(KIND=C_INT),       VALUE       :: ncid, varid
 Type(C_PTR),               VALUE       :: startp, countp, stridep
 Integer(KIND=C_LONG_LONG), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                    :: nc_get_vars_longlong

 End Function nc_get_vars_longlong
End Interface
!------------------------------- nc_put_var_longlong --------------------------
Interface
 Function nc_put_var_longlong(ncid, varid, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG_LONG

 Integer(KIND=C_INT),       VALUE      :: ncid, varid
 Integer(KIND=C_LONG_LONG), Intent(IN) :: op(*)

 Integer(KIND=C_INT)                   :: nc_put_var_longlong

 End Function nc_put_var_longlong
End Interface
!------------------------------- nc_get_var_longlong --------------------------
Interface
 Function nc_get_var_longlong(ncid, varid, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_LONG_LONG

 Integer(KIND=C_INT),       VALUE       :: ncid, varid
 Integer(KIND=C_LONG_LONG), Intent(OUT) :: ip(*)

 Integer(KIND=C_INT)                    :: nc_get_var_longlong

 End Function nc_get_var_longlong
End Interface
!------------------------------- nc_set_chunk_cache_ints ----------------------
Interface
 Function nc_set_chunk_cache_ints(size, nelems, preemption) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: size, nelems, preemption 

 Integer(KIND=C_INT)        :: nc_set_chunk_cache_ints

 End Function nc_set_chunk_cache_ints
End Interface
!------------------------------- nc_get_chunk_cache_ints ----------------------
Interface
 Function nc_get_chunk_cache_ints(size, nelems, preemption) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), Intent(INOUT) :: size, nelems, preemption 

 Integer(KIND=C_INT)                :: nc_get_chunk_cache_ints

 End Function nc_get_chunk_cache_ints
End Interface
!------------------------------- nc_set_var_chunk_cache_ints ------------------
Interface
 Function nc_set_var_chunk_cache_ints(ncid, varid, size, nelems, preemption) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: ncid, varid, size, nelems, preemption 

 Integer(KIND=C_INT)        :: nc_set_var_chunk_cache_ints

 End Function nc_set_var_chunk_cache_ints
End Interface
!------------------------------- nc_get_var_chunk_cache_ints ------------------
Interface
 Function nc_get_var_chunk_cache_ints(ncid, varid, size, nelems, preemption) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE         :: ncid, varid
 Integer(KIND=C_INT), Intent(INOUT) :: size, nelems, preemption
 Integer(KIND=C_INT)                :: nc_get_var_chunk_cache_ints

 End Function nc_get_var_chunk_cache_ints
End Interface
!------------------------------- nc_set_chunk_cache ---------------------------
Interface
 Function nc_set_chunk_cache(size, nelems, preemption) BIND(C)

 USE ISO_C_BINDING, ONLY: C_SIZE_T, C_FLOAT, C_INT

 Integer(KIND=C_SIZE_T), VALUE :: size, nelems
 Real(KIND=C_FLOAT),     VALUE :: preemption 

 Integer(KIND=C_INT)           :: nc_set_chunk_cache

 End Function nc_set_chunk_cache
End Interface
!------------------------------- nc_get_chunk_cache ---------------------------
Interface
 Function nc_get_chunk_cache(size, nelems, preemption) BIND(C)

 USE ISO_C_BINDING, ONLY: C_SIZE_T, C_FLOAT, C_INT

 Integer(KIND=C_SIZE_T), Intent(INOUT) :: size, nelems
 Real(KIND=C_FLOAT),     Intent(INOUT) :: preemption 

 Integer(KIND=C_INT)                   :: nc_get_chunk_cache

 End Function nc_get_chunk_cache
End Interface
!---------------------------------- nc_put_var --------------------------------
Interface
 Function nc_put_var(ncid, varid, op) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR

 Integer(KIND=C_INT), VALUE :: ncid, varid
 Type(C_PTR),         VALUE :: op

 Integer(KIND=C_INT)        :: nc_put_var

 End Function nc_put_var
End Interface
!---------------------------------- nc_get_var --------------------------------
Interface
 Function nc_get_var(ncid, varid, ip) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid, varid
 Character(KIND=C_CHAR), Intent(INOUT) :: ip(*)

 Integer(KIND=C_INT)                   :: nc_get_var

 End Function nc_get_var
End Interface

!--------------------------End of Module netcdf4_c_interfaces -----------------
End Module netcdf4_nc_interfaces
