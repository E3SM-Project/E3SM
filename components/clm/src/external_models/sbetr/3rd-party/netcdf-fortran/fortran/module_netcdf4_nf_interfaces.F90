Module netcdf4_nf_interfaces

! Explicit interfaces for Netcdf4 FORTRAN 2003 nf FORTRAN interface routines

! We exclude functions for netCDF4 that pass data of any type to/from void
! pointers in C using a C_CHAR string array. We do provide external statements
! for these routines 

! Written by: Richard Weed, Ph.D.
!             Center for Advanced Vehicular Systems 
!             Mississippi State University
!             rweed@cavs.msstate.edu
!
! The author grants to the University Center for Atmospheric Research
! (UCAR), Boulder, CO, USA the right to revise and extend the software
! without restriction. However, the author retains all copyrights and
! intellectual property rights explicitly stated in or implied by the
! Apache license

! Version 1. April 2009 - Initial module based on netCDF 4.0.1 interfaces
!                         A separate module is required to avoid a
!                         module dependency problem
! Version 2. April 2010 - Updated to Netcdf 4.1.1
! Version 3. Aug.  2013 - Made nf_def_var_fill and nf_inq_var_fill external
!                         to support new nf90_ non-integer fill characters
!                         Added interface for new nf_rename_grp function 
!                         Changed interface type defs to USE netcdf_nf_data

! Most legacy programs don't need to use this module. However, I've created
! it to support FORTRAN programmers who like to provide explicit interfaces
! for all subroutines and functions in their codes. Therefore, this module is
! is primarily for people writing new programs.

!           Explicit interfaces to netCDF 4 specific FORTRAN functions

!-------------------------------- nf_create_par -------------------------------
Interface
 Function nf_create_par (path, cmode, comm, info, ncid) RESULT(status)

 Integer,          Intent(IN)  :: cmode, comm, info
 Character(LEN=*), Intent(IN)  :: path
 Integer,          Intent(OUT) :: ncid
 Integer                       :: status

 End Function nf_create_par
End Interface
!-------------------------------- nf_open_par --------------------------------
Interface
 Function nf_open_par (path, mode, comm, info, ncid) RESULT(status)

 Integer,          Intent(IN)  :: mode, comm, info
 Character(LEN=*), Intent(IN)  :: path
 Integer,          Intent(OUT) :: ncid
 Integer                       :: status

 End Function nf_open_par
End Interface
!-------------------------------- nf_var_par_access -------------------------
Interface
 Function nf_var_par_access( ncid, varid, iaccess) RESULT (status)

 Integer, Intent(IN) :: ncid, varid, iaccess
 Integer             :: status

 End Function nf_var_par_access
End Interface
!-------------------------------- nf_inq_ncid ---------------------------------
Interface
 Function nf_inq_ncid( ncid, name, groupid) RESULT (status)

 Integer,          Intent(IN)  :: ncid
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: groupid
 Integer                       :: status

 End Function nf_inq_ncid
End Interface
!-------------------------------- nf_inq_grps ---------------------------------
Interface
 Function nf_inq_grps( ncid, numgrps, ncids) RESULT (status)

 Integer, Intent(IN)    :: ncid
 Integer, Intent(OUT)   :: numgrps
 Integer, Intent(INOUT) :: ncids(*)
 Integer                :: status

 End Function nf_inq_grps
End Interface
!-------------------------------- nf_inq_grpname ------------------------------
Interface
 Function nf_inq_grpname( ncid, name) RESULT (status)

 Integer,          Intent(IN)  :: ncid
 Character(LEN=*), Intent(OUT) :: name
 Integer                       :: status

 End Function nf_inq_grpname
End Interface
!-------------------------------- nf_inq_grpname_full -------------------------
Interface
 Function nf_inq_grpname_full( ncid, nlen, name) RESULT (status)

 Integer,          Intent(IN)  :: ncid
 Integer,          Intent(OUT) :: nlen 
 Character(LEN=*), Intent(OUT) :: name
 Integer                       :: status

 End Function nf_inq_grpname_full
End Interface
!-------------------------------- nf_inq_grpname_len -------------------------
Interface
 Function nf_inq_grpname_len( ncid, nlen) RESULT (status)

 Integer, Intent(IN)  :: ncid
 Integer, Intent(OUT) :: nlen 
 Integer              :: status

 End Function nf_inq_grpname_len
End Interface
!-------------------------------- nf_inq_grp_parent ---------------------------
Interface
 Function nf_inq_grp_parent( ncid,parent_ncid) RESULT (status)

 Integer, Intent(IN)    :: ncid
 Integer, Intent(INOUT) :: parent_ncid
 Integer                :: status

 End Function nf_inq_grp_parent
End Interface
!-------------------------------- nf_inq_grp_full_ncid ------------------------
Interface
 Function nf_inq_grp_full_ncid( ncid, grp_name, grp_ncid) RESULT (status)

 Integer,          Intent(IN)    :: ncid
 Character(LEN=*), Intent(IN)    :: grp_name
 Integer,          Intent(INOUT) :: grp_ncid
 Integer                         :: status

 End Function nf_inq_grp_full_ncid
End Interface
!-------------------------------- nf_inq_grp_ncid -----------------------------
Interface
 Function nf_inq_grp_ncid( ncid, grp_name, parent_ncid) RESULT (status)

 Integer,          Intent(IN)    :: ncid
 Character(LEN=*), Intent(IN)    :: grp_name
 Integer,          Intent(INOUT) :: parent_ncid
 Integer                         :: status

 End Function nf_inq_grp_ncid
End Interface
!-------------------------------- nf_inq_varids -------------------------------
Interface
 Function nf_inq_varids( ncid, nvars, varids) RESULT (status)

 Integer, Intent(IN)    :: ncid
 Integer, Intent(OUT)   :: nvars 
 Integer, Intent(INOUT) :: varids(*)
 Integer                :: status

 End Function nf_inq_varids
End Interface
!-------------------------------- nf_inq_dimids -------------------------------
Interface
 Function nf_inq_dimids( ncid, ndims, dimids, include_parents) RESULT (status)

 Integer, Intent(IN)    :: ncid, include_parents
 Integer, Intent(OUT)   :: ndims
 Integer, Intent(INOUT) :: dimids(*)
 Integer                :: status

 End Function nf_inq_dimids
End Interface
!-------------------------------- nf_inq_typeids ------------------------------
Interface
 Function nf_inq_typeids( ncid, ntypes, typeids) RESULT (status)

 Integer, Intent(IN)    :: ncid
 Integer, Intent(OUT)   :: ntypes
 Integer, Intent(INOUT) :: typeids(*)
 Integer                :: status

 End Function nf_inq_typeids
End Interface
!-------------------------------- nf_inq_typeid ------------------------------
Interface
 Function nf_inq_typeid(ncid, name, typeid) RESULT (status)

 Integer,          Intent(IN)  :: ncid 
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: typeid
 Integer                       :: status

 End Function nf_inq_typeid
End Interface
!-------------------------------- nf_def_grp ---------------------------------
Interface
 Function nf_def_grp( parent_ncid, name, new_ncid) RESULT (status)

 Integer,          Intent(IN)  :: parent_ncid
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: new_ncid
 Integer                       :: status

 End Function nf_def_grp
End Interface
!-------------------------------- nf_rename_grp -------------------------------
Interface
 Function nf_rename_grp( grpid, name) RESULT (status)

! rename previously defined group

 USE netcdf4_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: grpid
 Character(LEN=*), Intent(IN)  :: name
 Integer                       :: status

 End Function nf_rename_grp
End Interface
!-------------------------------- nf_def_compound -----------------------------
Interface
 Function nf_def_compound( ncid, isize, name, typeid) RESULT (status)

 Integer,          Intent(IN)  :: ncid, isize
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: typeid 
 Integer                       :: status

 End Function nf_def_compound
End Interface
!-------------------------------- nf_insert_compound --------------------------
Interface
 Function nf_insert_compound( ncid, xtype, name, offset, field_typeid) &
                              RESULT (status)

 Integer,          Intent(IN) :: ncid, xtype, field_typeid, offset 
 Character(LEN=*), Intent(IN) :: name
 Integer                      :: status

 End Function nf_insert_compound
End Interface
!-------------------------------- nf_insert_array_compound --------------------
Interface
 Function nf_insert_array_compound( ncid, xtype, name, offset, field_typeid, &
                                    ndims, dim_sizes) RESULT (status)

 Integer,          Intent(IN)    :: ncid, xtype, field_typeid, offset, ndims
 Integer,          Intent(INOUT) :: dim_sizes(*)
 Character(LEN=*), Intent(IN)    :: name
 Integer                         :: status

 End Function nf_insert_array_compound
End Interface
!-------------------------------- nf_inq_type ---------------------------------
Interface
 Function nf_inq_type( ncid, xtype, name, isize) RESULT (status)

 Integer,          Intent(IN)  :: ncid, xtype
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: isize
 Integer                       :: status

 End Function nf_inq_type
End Interface
!-------------------------------- nf_inq_compound -----------------------------
Interface
 Function nf_inq_compound( ncid, xtype, name, isize, nfields) RESULT (status)

 Integer,          Intent(IN)    :: ncid, xtype
 Character(LEN=*), Intent(INOUT) :: name
 Integer,          Intent(INOUT) :: isize, nfields
 Integer                         :: status

 End Function nf_inq_compound
End Interface
!-------------------------------- nf_inq_compound_name ------------------------
Interface
 Function nf_inq_compound_name( ncid, xtype, name) RESULT (status)

 Integer,          Intent(IN)  :: ncid, xtype
 Character(LEN=*), Intent(OUT) :: name
 Integer                       :: status

 End Function nf_inq_compound_name
End Interface
!-------------------------------- nf_inq_compound_size -------------------------
Interface
 Function nf_inq_compound_size( ncid, xtype, isize) RESULT (status)

 Integer, Intent(IN)    :: ncid, xtype
 Integer, Intent(INOUT) :: isize
 Integer                :: status

 End Function nf_inq_compound_size
End Interface
!-------------------------------- nf_inq_compound_nfields ----------------------
Interface
 Function nf_inq_compound_nfields( ncid, xtype, nfields) RESULT (status)

 Integer, Intent(IN)    :: ncid, xtype
 Integer, Intent(INOUT) :: nfields
 Integer                :: status

 End Function nf_inq_compound_nfields
End Interface
!-------------------------------- nf_inq_compound_field -----------------------
Interface
 Function nf_inq_compound_field( ncid, xtype, fieldid, name, offset, &
                                field_typeid, ndims, dim_sizes) RESULT (status)

 Integer,          Intent(IN)  :: ncid, xtype, fieldid
 Character(LEN=*), Intent(OUT) :: name
 Integer,          Intent(OUT) :: offset, field_typeid, ndims
 Integer,          Intent(OUT) :: dim_sizes(*)
 Integer                       :: status

 End Function nf_inq_compound_field
End Interface
!-------------------------------- nf_inq_compound_fieldname -------------------
Interface
 Function nf_inq_compound_fieldname(ncid, xtype, fieldid, name) RESULT(status)

 Integer,          Intent(IN)  :: ncid, xtype, fieldid
 Character(LEN=*), Intent(OUT) :: name
 Integer                       :: status

 End Function nf_inq_compound_fieldname
End Interface
!-------------------------------- nf_inq_compound_fieldindex ------------------
Interface
 Function nf_inq_compound_fieldindex( ncid, xtype, name, fieldid) &
                                      RESULT (status)

 Integer,          Intent(IN)  :: ncid, xtype
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: fieldid
 Integer                       :: status

 End Function nf_inq_compound_fieldindex
End Interface
!-------------------------------- nf_inq_compound_fieldoffset ----------------
Interface
 Function nf_inq_compound_fieldoffset( ncid, xtype, fieldid, offset) &
                                      RESULT (status)

 Integer, Intent(IN)  :: ncid, xtype, fieldid
 Integer, Intent(OUT) :: offset
 Integer              :: status

 End Function nf_inq_compound_fieldoffset
End Interface
!-------------------------------- nf_inq_compound_fieldtype -------------------
Interface
 Function nf_inq_compound_fieldtype( ncid, xtype, fieldid, field_typeid) &
                                     RESULT (status)

 Integer, Intent(IN)  :: ncid, xtype, fieldid
 Integer, Intent(OUT) :: field_typeid
 Integer              :: status

 End Function nf_inq_compound_fieldtype
End Interface
!-------------------------------- nf_inq_compound_fieldndims ------------------
Interface
 Function nf_inq_compound_fieldndims( ncid, xtype, fieldid, ndims) &
                                      RESULT (status)

 Integer, Intent(IN)  :: ncid, xtype, fieldid
 Integer, Intent(OUT) :: ndims
 Integer              :: status

 End Function nf_inq_compound_fieldndims
End Interface
!-------------------------------- nf_inq_compound_fielddim_sizes --------------
Interface
 Function nf_inq_compound_fielddim_sizes( ncid, xtype, fieldid, dim_sizes) &
                                          RESULT (status)

 Integer, Intent(IN)    :: ncid, xtype, fieldid
 Integer, Intent(INOUT) :: dim_sizes(*)
 Integer                :: status

 End Function nf_inq_compound_fielddim_sizes
End Interface
!-------------------------------- nf_def_vlen ---------------------------------
Interface
 Function nf_def_vlen( ncid, name, base_typeid, xtype) RESULT (status)

 Integer,          Intent(IN)  :: ncid, base_typeid 
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: xtype
 Integer                       :: status

 End Function nf_def_vlen
End Interface
!-------------------------------- nf_inq_vlen ---------------------------------
Interface
 Function nf_inq_vlen( ncid, xtype, name, datum_size, base_type) RESULT(status) 

 Integer,          Intent(IN)  :: ncid, xtype
 Character(LEN=*), Intent(OUT) :: name
 Integer,          Intent(OUT) :: datum_size, base_type
 Integer                       :: status

 End Function nf_inq_vlen
End Interface
!-------------------------------- nf_inq_user_type ----------------------------
Interface
 Function nf_inq_user_type( ncid, xtype, name, isize, base_type, nfields, &
                            iclass) RESULT (status)

 Integer,          Intent(IN)    :: ncid, xtype
 Character(LEN=*), Intent(INOUT) :: name
 Integer,          Intent(OUT)   :: isize, nfields, base_type, iclass
 Integer                         :: status

 End Function nf_inq_user_type
End Interface
!-------------------------------- nf_def_enum ---------------------------------
Interface
 Function nf_def_enum( ncid, base_typeid, name, typeid) RESULT (status)

 Integer,          Intent(IN)  :: ncid, base_typeid 
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: typeid
 Integer                       :: status

 End Function nf_def_enum
End Interface
!-------------------------------- nf_insert_enum -------------------------------
! Commented out for now since we pass to/from a void pointer using
! a C_CHAR string which will be non-compatible with different data
! types being passed
!Interface
! Function nf_insert_enum( ncid, xtype, name, value) RESULT (status)

! USE netcdf_nf_data

! Integer,                Intent(IN)         :: ncid, xtype
! Character(LEN=*),       Intent(IN)         :: name
! Character(KIND=C_CHAR), Intent(IN), TARGET :: value(*)
! Integer                                    :: status

! End Function nf_insert_enum
!End Interface
!-------------------------------- nf_inq_enum ---------------------------------
Interface
 Function nf_inq_enum( ncid, xtype, name, base_nf_type, base_size, &
                       num_members) RESULT (status)

! USE netcdf_nf_data

 Integer,          Intent(IN)    :: ncid, xtype
 Character(LEN=*), Intent(INOUT) :: name
 Integer,          Intent(INOUT) :: base_nf_type, base_size, num_members
 Integer                         :: status

 End Function nf_inq_enum
End Interface
!-------------------------------- nf_inq_enum_member --------------------------
! Commented out for now since we pass to/from a void pointer using
! a C_CHAR string which will be non-compatible with different data
! types being passed
!Interface
! Function nf_inq_enum_member( ncid, xtype, idx, name, value) RESULT (status)

! USE netcdf_nf_data

! Integer,                Intent(IN)  :: ncid, xtype, idx
! Character(LEN=*),       Intent(OUT) :: name
! Character(KIND=C_CHAR), Intent(OUT) :: value(*)
! Integer                             :: status

! End Function nf_inq_enum_member
!End Interface
!-------------------------------- nf_inq_enum_ident --------------------------
Interface
 Function nf_inq_enum_ident( ncid, xtype, value, name) RESULT (status)

 Integer,          Intent(IN)    :: ncid, xtype, value 
 Character(LEN=*), Intent(INOUT) :: name
 Integer                         :: status

 End Function nf_inq_enum_ident
End Interface
!-------------------------------- nf_def_opaque -------------------------------
Interface
 Function nf_def_opaque( ncid, isize, name, xtype) RESULT (status)

 Integer,          Intent(IN)  :: ncid, isize
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: xtype
 Integer                       :: status

 End Function nf_def_opaque
End Interface
!-------------------------------- nf_inq_opaque -------------------------------
Interface
 Function nf_inq_opaque( ncid, xtype, name, isize) RESULT (status)

 Integer,          Intent(IN)    :: ncid, xtype
 Character(LEN=*), Intent(INOUT) :: name
 Integer,          Intent(OUT)   :: isize
 Integer                         :: status

 End Function nf_inq_opaque
End Interface
!-------------------------------- nf_def_var_chunking -------------------------
Interface
 Function nf_def_var_chunking( ncid, varid, contiguous, chunksizes) &
                               RESULT(status)

 Integer, Intent(IN)    :: ncid, varid, contiguous
 Integer, Intent(INOUT) :: chunksizes(*)
 Integer                :: status

 End Function nf_def_var_chunking
End Interface
!-------------------------------- nf_inq_var_chunking -------------------------
Interface
 Function nf_inq_var_chunking( ncid, varid, contiguous, chunksizes) &
                               RESULT (status)

 Integer, Intent(IN)    :: ncid, varid
 Integer, Intent(INOUT) :: contiguous
 Integer, Intent(INOUT) :: chunksizes(*)
 Integer                :: status

 End Function nf_inq_var_chunking
End Interface
!-------------------------------- nf_def_var_deflate --------------------------
Interface
 Function nf_def_var_deflate( ncid, varid, shuffle, deflate, deflate_level) &
                               RESULT (status)

 Integer, Intent(IN) :: ncid, varid, shuffle, deflate, deflate_level
 Integer             :: status

 End Function nf_def_var_deflate
End Interface
!-------------------------------- nf_inq_var_deflate -------------------------
Interface
 Function nf_inq_var_deflate( ncid, varid, shuffle, deflate, deflate_level) &
                               RESULT (status)

 Integer, Intent(IN)  :: ncid, varid
 Integer, Intent(OUT) :: shuffle, deflate, deflate_level
 Integer              :: status

 End Function nf_inq_var_deflate
End Interface
!-------------------------------- nf_inq_var_szip -----------------------------
Interface
 Function nf_inq_var_szip(ncid, varid, options_mask, pixels_per_block) RESULT(status)

 Implicit NONE

 Integer, Intent(IN)    :: ncid, varid
 Integer, Intent(INOUT) :: options_mask, pixels_per_block
 Integer                :: status

 End Function nf_inq_var_szip
End Interface
!-------------------------------- nf_def_var_fletcher32 ------------------------
Interface
 Function nf_def_var_fletcher32( ncid, varid, fletcher32) RESULT(status)

 Integer, Intent(IN) :: ncid, varid, fletcher32
 Integer             :: status

 End Function nf_def_var_fletcher32
End Interface
!-------------------------------- nf_inq_var_fletcher32 ------------------------
Interface
 Function nf_inq_var_fletcher32( ncid, varid, fletcher32) RESULT(status)

 Integer, Intent(IN)  :: ncid, varid
 Integer, Intent(OUT) :: fletcher32
 Integer              :: status

 End Function nf_inq_var_fletcher32
End Interface
!-------------------------------- nf_def_var_fill -----------------------------
!Interface
! Function nf_def_var_fill( ncid, varid, no_fill, fill_value) RESULT(status)

! Integer, Intent(IN)    :: ncid, varid, no_fill
! Integer, Intent(IN)    :: fill_value 
! Integer                :: status

! End Function nf_def_var_fill
!End Interface
!-------------------------------- nf_inq_var_fill -----------------------------
!Interface
! Function nf_inq_var_fill( ncid, varid, no_fill, fill_value) RESULT(status)

! Integer, Intent(IN)    :: ncid, varid
! Integer, Intent(OUT)   :: no_fill
! Integer, Intent(INOUT) :: fill_value 
! Integer                :: status

! End Function nf_inq_var_fill
!End Interface
!-------------------------------- nf_def_var_endian ---------------------------
Interface
 Function nf_def_var_endian( ncid, varid, endiann) RESULT(status)

 Integer, Intent(IN) :: ncid, varid, endiann
 Integer             :: status

 End Function nf_def_var_endian
End Interface
!-------------------------------- nf_inq_var_endian ---------------------------
Interface
 Function nf_inq_var_endian( ncid, varid, endiann) RESULT(status)

 Integer, Intent(IN)  :: ncid, varid
 Integer, Intent(OUT) :: endiann
 Integer              :: status

 End Function nf_inq_var_endian
End Interface
!--------------------------------- nf_put_att --------------------------------
! Commented out because we use C_CHAR array to pass data of different
! type to a C void pointer
!Interface
! Function nf_put_att(ncid, varid, name, xtype, nlen, value) RESULT(status)
 
! USE netcdf_nf_data

! Integer,                Intent(IN)         :: ncid, varid, nlen, xtype
! Character(LEN=*),       Intent(IN)         :: name
! Character(KIND=C_CHAR), Intent(IN), TARGET :: value(*)
! Integer                                    :: status

! End Function nf_put_att
!--------------------------------- nf_get_att --------------------------------
! Commented out because we use C_CHAR array to pass data of different
! type to a C void pointer
!Interface
! Function nf_get_att(ncid, varid, name, value) RESULT(status)

! USE netcdf_nf_data

! Implicit NONE

! Integer,                Intent(IN)    :: ncid, varid
! Character(LEN=*),       Intent(IN)    :: name
! Character(KIND=C_CHAR), Intent(INOUT) :: value(*)
! Integer                               :: status

! End Function nf_get_att
!End Interface
!--------------------------------- nf_put_vlen_element ------------------------
! Commented out because we use C_CHAR array to pass data of different
! type to a C void pointer
!Interface
! Function nf_put_vlen_element(ncid, xtype, vlen_element, nlen, value)&
!                               RESULT(status)

! USE netcdf_nf_data

! Integer,                Intent(IN)            :: ncid, xtype, nlen
! Character(LEN=*),       Intent(INOUT), TARGET :: vlen_element
! Character(KIND=C_CHAR), Intent(IN),    TARGET :: value(*)
! Integer                                       :: status

! End Function nf_put_vlen_element
!End Interface
!--------------------------------- nf_get_vlen_element ------------------------
! Commented out because we use C_CHAR array to pass data of different
! type to a C void pointer
!Interface
! Function nf_get_vlen_element(ncid, xtype, vlen_element, nlen, value) &
!                              RESULT(status)

! USE netcdf_nf_data

! Implicit NONE

! Integer,                Intent(IN)            :: ncid, xtype
! Integer,                Intent(INOUT)         :: nlen
! Character(LEN=*),       Intent(INOUT), TARGET :: vlen_element
! Character(KIND=C_CHAR), Intent(INOUT)         :: value(*)
! Integer                                       :: status

! End Function nf_get_vlen_element
!End Interface
!--------------------------------- nf_free_vlenn -------------------------------
! Commented out because we use C_CHAR array to pass data of different
! type to a C void pointer
!Interface
! Function nf_free_vlen(vl) RESULT(status)

! Character(KIND=C_CHAR), Intent(IN), TARGET :: vl(*)
! Integer                                    :: status

! End Function nf_free_vlen
!End Interface
!--------------------------------- nf_put_var ---------------------------------
! Commented out because we use C_CHAR array to pass data of different
! type to a C void pointer
!Interface
! Function nf_put_var(ncid, varid, values) RESULT(status)

! USE netcdf_nf_data

! Integer,                Intent(IN)         :: ncid, varid
! Character(KIND=C_CHAR), Intent(IN), TARGET :: values(*)
! Integer                                    :: status

! End Function nf_put_var
!End Interface
!--------------------------------- nf_get_var ---------------------------------
! Commented out because we use C_CHAR array to pass data of different
! type to a C void pointer
!Interface
! Function nf_get_var(ncid, varid, values) RESULT(status)

! USE netcdf_nf_data

! Integer,                Intent(IN)         :: ncid, varid
! Character(KIND=C_CHAR), Intent(INOUT),     :: values(*)
! Integer                                    :: status

! End Function nf_get_var
!End Interface
!--------------------------------- nf_put_var1_int64 --------------------------
Interface
 Function nf_put_var1_int64(ncid, varid, ndex, ival) RESULT(status)

 USE netcdf_nf_data, ONLY: IK8

 Integer,           Intent(IN) :: ncid, varid
 Integer,           Intent(IN) :: ndex(*)
 Integer(KIND=IK8), Intent(IN) :: ival
 Integer                       :: status

 End Function nf_put_var1_int64
End Interface
!--------------------------------- nf_put_vara_int64 --------------------------
Interface
 Function nf_put_vara_int64(ncid, varid, start, counts, ivals) RESULT(status)

 USE netcdf_nf_data, ONLY: IK8

 Integer,           Intent(IN) :: ncid, varid
 Integer,           Intent(IN) :: start(*), counts(*)
 Integer(KIND=IK8), Intent(IN) :: ivals(*)
 Integer                       :: status

 End Function nf_put_vara_int64
End Interface
!--------------------------------- nf_put_vars_int64 --------------------------
Interface
 Function nf_put_vars_int64(ncid, varid, start, counts, strides, ivals) &
                             RESULT(status)
 
 USE netcdf_nf_data, ONLY: IK8

 Integer,           Intent(IN) :: ncid, varid
 Integer,           Intent(IN) :: start(*), counts(*), strides(*)
 Integer(KIND=IK8), Intent(IN) :: ivals(*)
 Integer                       :: status

 End Function nf_put_vars_int64
End Interface
!--------------------------------- nf_put_varm_int64 -------------------------
Interface
 Function nf_put_varm_int64(ncid, varid, start, counts, strides, maps, &
                            ivals) RESULT(status)

 USE netcdf_nf_data, ONLY: IK8

 Integer,           Intent(IN) :: ncid, varid
 Integer,           Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Integer(KIND=IK8), Intent(IN) :: ivals(*)
 Integer                       :: status

 End Function nf_put_varm_int64
End Interface
!--------------------------------- nf_put_var_int64 --------------------------
Interface
 Function nf_put_var_int64(ncid, varid, ivals) RESULT(status)

 USE netcdf_nf_data, ONLY: IK8

 Integer,           Intent(IN) :: ncid, varid
 Integer(KIND=IK8), Intent(IN) :: ivals(*)
 Integer                       :: status

 End Function nf_put_var_int64
End Interface
!--------------------------------- nf_get_var1_int64 -------------------------
Interface
 Function nf_get_var1_int64(ncid, varid, ndex, ival) RESULT(status)

 USE netcdf_nf_data, ONLY: IK8

 Integer,           Intent(IN)  :: ncid, varid
 Integer,           Intent(IN)  :: ndex(*)
 Integer(KIND=IK8), Intent(OUT) :: ival
 Integer                        :: status

 End Function nf_get_var1_int64
End Interface
!--------------------------------- nf_get_vara_int -------------------------
Interface
 Function nf_get_vara_int64(ncid, varid, start, counts, ivals) RESULT(status)

 USE netcdf_nf_data, ONLY: IK8

 Integer,           Intent(IN)  :: ncid, varid
 Integer,           Intent(IN)  :: start(*), counts(*)
 Integer(KIND=IK8), Intent(OUT) :: ivals(*)
 Integer                        :: status

 End Function nf_get_vara_int64
End Interface
!--------------------------------- nf_get_vars_int64 --------------------------
Interface
 Function nf_get_vars_int64(ncid, varid, start, counts, strides, ivals) &
                             RESULT(status)

 USE netcdf_nf_data, ONLY: IK8

 Integer,           Intent(IN)  :: ncid, varid
 Integer,           Intent(IN)  :: start(*), counts(*), strides(*)
 Integer(KIND=IK8), Intent(OUT) :: ivals(*)
 Integer                        :: status

 End Function nf_get_vars_int64
End Interface
!--------------------------------- nf_get_varm_int64 -------------------------
Interface
 Function nf_get_varm_int64(ncid, varid, start, counts, strides, maps, &
                            ivals) RESULT(status)

 USE netcdf_nf_data, ONLY: IK8

 Integer,           Intent(IN)  :: ncid, varid
 Integer,           Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Integer(KIND=IK8), Intent(OUT) :: ivals(*)
 Integer                        :: status

 End Function nf_get_varm_int64
End Interface
!--------------------------------- nf_get_var_int64 --------------------------
Interface
 Function nf_get_var_int64(ncid, varid, ivals) RESULT(status)

 USE netcdf_nf_data, ONLY: IK8

 Integer,           Intent(IN)  :: ncid, varid
 Integer(KIND=IK8), Intent(OUT) :: ivals(*)
 Integer                        :: status

 End Function nf_get_var_int64
End Interface
!--------------------------------- nf_set_chunk_cache -------------------------
Interface
 Function nf_set_chunk_cache(chunk_size, nelems, preemption) RESULT(status)

 Integer, Intent(IN) :: chunk_size, nelems, preemption
 Integer             :: status

 End Function nf_set_chunk_cache
End Interface
!--------------------------------- nf_get_chunk_cache -------------------------
Interface
 Function nf_get_chunk_cache(chunk_size, nelems, preemption) RESULT(status)

 Integer, Intent(INOUT) :: chunk_size, nelems, preemption
 Integer                :: status

 End Function nf_get_chunk_cache
End Interface
!--------------------------------- nf_set_var_chunk_cache ---------------------
Interface
 Function nf_set_var_chunk_cache(ncid, varid, chunk_size, nelems, preemption) RESULT(status)

! USE netcdf_nf_data

 Implicit NONE

 Integer, Intent(IN) :: ncid, varid, chunk_size, nelems, preemption
 Integer             :: status

 End Function nf_set_var_chunk_cache
End Interface
!--------------------------------- nf_get_var_chunk_cache ---------------------
Interface
 Function nf_get_var_chunk_cache(ncid, varid, chunk_size, nelems, preemption) RESULT(status)

! get chunk cache size. Note this follows the fort-nc4 version which uses
! uses nc_get_var_chunk_cache_ints to avoid size_t issues with fortran. 

! USE netcdf_nf_data

 Implicit NONE

 Integer, Intent(IN)    :: ncid, varid
 Integer, Intent(INOUT) :: chunk_size, nelems, preemption
 Integer                :: status

 End Function nf_get_var_chunk_cache
End Interface

! Declare external values for functions that use C_CHAR strings to pass
! data of different types

 Integer, External :: nf_insert_enum
 Integer, External :: nf_inq_enum_member
 Integer, External :: nf_put_att
 Integer, External :: nf_get_att
 Integer, External :: nf_put_vlen_element
 Integer, External :: nf_get_vlen_element
 Integer, External :: nf_free_vlen
 Integer, External :: nf_free_vlens
 Integer, External :: nf_free_string
 Integer, External :: nf_put_var
 Integer, External :: nf_get_var
 Integer, External :: nf_def_var_fill
 Integer, External :: nf_inq_var_fill

End Module netcdf4_nf_interfaces
