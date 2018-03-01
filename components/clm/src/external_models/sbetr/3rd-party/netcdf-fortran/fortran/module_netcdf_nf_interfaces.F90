Module netcdf_nf_interfaces

! Explicit interfaces for Netcdf FORTRAN 2003 nf FORTRAN interface routines
! Generic interfaces are provided for routines that process text data to
! handle the case where you are using this interface module and are passing
! an array of single characters (ala C) instead of a character string.

! We choose not to provide explicit interfaces for the V2 routines and
! newer functions for netCDF4 that pass data of any type to/from void
! pointers in C using a C_CHAR string array. We do provide external
! statements for the V2 functions as done in netcdf2.inc

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

! Version 1. Sept. 2005 - Initial Cray X1 version
! Version 2. May, 2006  - Updated to support g95
! Version 3. April 2009 - Updated for netCDF 4.0.1
! Version 4. April 2010 - Updated for netCDF 4.1.1
! Version 5. Feb.  2013 - Added nf_inq_path support for fortran 4.4
         
! Most legacy programs don't need to use this module. However, I've created
! it to support FORTRAN programmers who like to provide explicit interfaces
! for all subroutines and functions in their codes. Therefore, this module is
! is primarily for people writting new programs.

 Implicit NONE

#include "nfconfig.inc"

!-------------------- Explicit Interfaces for nf routines ------------------

! Misc functions first
!-------------------------------- nf_inq_libvers -------------------------------
Interface
 Function nf_inq_libvers() RESULT(vermsg)

 Character(LEN=80) :: vermsg

 End Function nf_inq_libvers
End Interface
!-------------------------------- nf_strerror ---------------------------------
Interface
 Function nf_strerror(nerr) RESULT(errmsg)

 Integer, Intent(IN) :: nerr

 Character(LEN=80)   :: errmsg

 End Function nf_strerror
End Interface
!-------------------------------- nf_issyserr ---------------------------------
Interface
 Function nf_issyserr(nerr) RESULT(status)

 Integer, Intent(IN) :: nerr
 Logical             :: status

 End Function nf_issyserr
End Interface

! Control routines
!-------------------------------- nf_create -----------------------------------
Interface
 Function nf_create(path, cmode, ncid) RESULT (status)

 Character(LEN=*), Intent(IN)  :: path
 Integer,          Intent(IN)  :: cmode
 Integer,          Intent(OUT) :: ncid
 Integer                       :: status

 End Function nf_create
End Interface
!-------------------------------- nf__create ----------------------------------
Interface
 Function nf__create(path, cmode, initialsz, chunksizehintp, ncid) &
                        RESULT(status)

 Character(LEN=*), Intent(IN)  :: path
 Integer,          Intent(IN)  :: cmode, initialsz, chunksizehintp
 Integer,          Intent(OUT) :: ncid
 Integer                       :: status

 End Function nf__create
End Interface
!-------------------------------- nf__create_mp -------------------------------
Interface
 Function nf__create_mp(path, cmode, initialsz, basepe, chunksizehintp, ncid) &
                        RESULT(status)

 Character(LEN=*), Intent(IN)  :: path
 Integer,          Intent(IN)  :: cmode, initialsz, chunksizehintp, basepe
 Integer,          Intent(OUT) :: ncid
 Integer                       :: status

 End Function nf__create_mp
End Interface
!-------------------------------- nf_open -------------------------------------
Interface
 Function nf_open(path, mode, ncid) RESULT (status)

 Character(LEN=*), Intent(IN)    :: path
 Integer,          Intent(IN)    :: mode
 Integer,          Intent(INOUT) :: ncid
 Integer                         :: status

 End Function nf_open
End Interface
!-------------------------------- nf__open ------------------------------------
Interface
 Function nf__open(path, mode, chunksizehintp, ncid) RESULT (status)

 Character(LEN=*), Intent(IN)    :: path
 Integer,          Intent(IN)    :: mode, chunksizehintp
 Integer,          Intent(INOUT) :: ncid
 Integer                         :: status

 End Function nf__open
End Interface
!-------------------------------- nf__open_mp ---------------------------------
Interface
 Function nf__open_mp(path, mode, basepe, chunksizehintp, ncid) RESULT (status)

 Character(LEN=*), Intent(IN)    :: path
 Integer,          Intent(IN)    :: mode, chunksizehintp, basepe
 Integer,          Intent(INOUT) :: ncid
 Integer                         :: status

 End Function nf__open_mp
End Interface
!-------------------------------- nf_inq_path ---------------------------------
Interface
 Function nf_inq_path(ncid, pathlen, path) RESULT (status)

 Integer,          Intent(IN)    :: ncid
 Integer,          Intent(INOUT) :: pathlen
 Character(LEN=*), Intent(INOUT) :: path
 Integer                         :: status

 End Function nf_inq_path
End Interface
!-------------------------------- nf_set_fill ---------------------------------
Interface
 Function nf_set_fill(ncid, fillmode, old_mode) RESULT(status)

 Integer, Intent(IN)  :: ncid, fillmode
 Integer, Intent(OUT) :: old_mode
 Integer              :: status

 End Function nf_set_fill
End Interface
!-------------------------------- nf_set_default_format -----------------------
Interface
 Function nf_set_default_format(newform, old_format) RESULT(status)

 Integer, Intent(IN)  :: newform
 Integer, Intent(OUT) :: old_format
 Integer              :: status

 End Function nf_set_default_format
End Interface
!-------------------------------- nf_redef -----------------------------------
Interface
 Function nf_redef(ncid) RESULT(status)

 Integer, Intent(IN) :: ncid
 Integer             :: status

 End Function nf_redef
End Interface
!-------------------------------- nf_enddef -----------------------------------
Interface
 Function nf_enddef(ncid) RESULT(status)

 Integer, Intent(IN) :: ncid
 Integer             :: status

 End Function nf_enddef
End Interface
!-------------------------------- nf__enddef ---------------------------------
Interface
 Function nf__enddef(ncid, h_minfree, v_align, v_minfree, r_align) &
                        RESULT(status)

 Integer, Intent(IN) :: ncid, h_minfree, v_align, v_minfree, r_align
 Integer             :: status

 End Function nf__enddef
End Interface

!-------------------------------- nf_sync -------------------------------------
Interface
 Function nf_sync(ncid) RESULT(status)

 Integer, Intent(IN) :: ncid
 Integer             :: status

 End Function nf_sync
End Interface
!-------------------------------- nf_abort -----------------------------------
Interface
 Function nf_abort(ncid) RESULT(status)

 Integer, Intent(IN) :: ncid
 Integer             :: status

 End Function nf_abort
End Interface
!-------------------------------- nf_close -------------------------------------
Interface
 Function nf_close(ncid) RESULT(status)

 Integer, Intent(IN) :: ncid
 Integer             :: status

 End Function nf_close
End Interface
!-------------------------------- nf_delete -----------------------------------
Interface
 Function nf_delete(path) RESULT(status)

 Character(LEN=*), Intent(IN) :: path
 Integer                      :: status

 End Function nf_delete
End Interface
!-------------------------------- nf_delete_mp ---------------------------------
Interface
 Function nf_delete_mp(path, pe) RESULT(status)

 Character(LEN=*), Intent(IN) :: path
 Integer,          Intent(IN) :: pe
 Integer                      :: status

 End Function nf_delete_mp
End Interface
!-------------------------------- nf_set_base_pe ------------------------------
Interface
 Function nf_set_base_pe(ncid, pe) RESULT(status)

 Integer, Intent(IN) :: ncid, pe
 Integer             :: status

 End Function nf_set_base_pe
End Interface
!-------------------------------- nf_inq_base_pe ------------------------------
Interface
 Function nf_inq_base_pe(ncid, pe) RESULT(status)

 Integer, Intent(IN)  :: ncid
 Integer, Intent(OUT) :: pe
 Integer              :: status

 End Function nf_inq_base_pe
End Interface

! Dimension definition and inquiry functions

!-------------------------------- nf_def_dim ----------------------------------
Interface
 Function nf_def_dim(ncid, name, dlen, dimid) RESULT (status)

 Integer,          Intent(IN)  :: ncid, dlen
 Integer,          Intent(OUT) :: dimid
 Character(LEN=*), Intent(IN)  :: name
 Integer                       :: status

 End Function nf_def_dim
End Interface
!-------------------------------- nf_inq_dim ----------------------------------
Interface
 Function nf_inq_dim(ncid, dimid, name, dlen) RESULT (status)

 Integer,          Intent(IN)   :: ncid, dimid
 Integer,          Intent(OUT)  :: dlen
 Character(LEN=*), Intent(OUT)  :: name
 Integer                        :: status

 End Function nf_inq_dim
End Interface
!-------------------------------- nf_inq_dimid --------------------------------
Interface
 Function nf_inq_dimid(ncid, name, dimid) RESULT (status)

 Integer,          Intent(IN)  :: ncid
 Integer,          Intent(OUT) :: dimid
 Character(LEN=*), Intent(IN)  :: name
 Integer                       :: status

 End Function nf_inq_dimid
End Interface
!-------------------------------- nf_inq_dimlen -------------------------------
Interface
 Function nf_inq_dimlen(ncid, dimid, dlen) RESULT (status)

 Integer, Intent(IN)  :: ncid, dimid
 Integer, Intent(OUT) :: dlen
 Integer              :: status

 End Function nf_inq_dimlen
End Interface
!-------------------------------- nf_inq_dimname ------------------------------
Interface
 Function nf_inq_dimname (ncid, dimid, name) RESULT (status)

 Integer,          Intent(IN)   :: ncid, dimid
 Character(LEN=*), Intent(OUT)  :: name
 Integer                        :: status

 End Function nf_inq_dimname
End Interface
!-------------------------------- nf_rename_dim --------------------------------
Interface
 Function nf_rename_dim(ncid, dimid, name) RESULT (status)

 Integer,          Intent(IN)  :: ncid, dimid
 Character(LEN=*), Intent(IN)  :: name
 Integer                       :: status

 End Function nf_rename_dim
End Interface

! General inquiry functions

!-------------------------------- nf_inq --------------------------------------
Interface
 Function nf_inq(ncid, ndims, nvars, ngatts, unlimdimid) RESULT(status)

 Integer, Intent(IN)  :: ncid
 Integer, Intent(OUT) :: ndims, nvars, ngatts, unlimdimid
 Integer              :: status

 End Function nf_inq
End Interface
!-------------------------------- nf_inq_ndims --------------------------------
Interface
 Function nf_inq_ndims(ncid, ndims) RESULT(status)

 Integer, Intent(IN)  :: ncid
 Integer, Intent(OUT) :: ndims
 Integer              :: status

 End Function nf_inq_ndims
End Interface
!-------------------------------- nf_inq_nvars --------------------------------
Interface
 Function nf_inq_nvars(ncid, nvars) RESULT(status)

 Integer, Intent(IN)  :: ncid
 Integer, Intent(OUT) :: nvars
 Integer              :: status

 End Function nf_inq_nvars
End Interface
!-------------------------------- nf_inq_natts --------------------------------
Interface
 Function nf_inq_natts(ncid, ngatts) RESULT(status)

 Integer, Intent(IN)  :: ncid
 Integer, Intent(OUT) :: ngatts
 Integer              :: status

 End Function nf_inq_natts
End Interface
!-------------------------------- nf_inq_unlimdim -----------------------------
Interface
 Function nf_inq_unlimdim(ncid, unlimdimid) RESULT(status)

 Integer, Intent(IN)  :: ncid
 Integer, Intent(OUT) :: unlimdimid
 Integer              :: status

 End Function nf_inq_unlimdim
End Interface
!-------------------------------- nf_inq_format -------------------------------
Interface
 Function nf_inq_format(ncid, format_type) RESULT(status)

 Integer, Intent(IN)  :: ncid
 Integer, Intent(OUT) :: format_type
 Integer              :: status

 End Function nf_inq_format
End Interface

! General variable functions

!-------------------------------- nf_def_var -----------------------------------
Interface
 Function nf_def_var(ncid, name, xtype, nvdims, vdims, varid) RESULT (status)

 Integer,          Intent(IN)  :: ncid, xtype, nvdims
 Integer,          Intent(IN)  :: vdims(*)
 Integer,          Intent(OUT) :: varid
 Character(LEN=*), Intent(IN)  :: name
 Integer                       :: status

 End Function nf_def_var
End Interface
!-------------------------------- nf_inq_varndims -----------------------------
Interface
 Function nf_inq_varndims(ncid, varid, vndims) RESULT (status)

 Integer, Intent(IN)  :: ncid, varid
 Integer, Intent(OUT) :: vndims
 Integer              :: status

 End Function nf_inq_varndims
End Interface
!-------------------------------- nf_inq_var ----------------------------------
Interface
 Function nf_inq_var(ncid, varid, name, xtype, ndims, dimids, natts) &
                        RESULT (status)

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(OUT) :: name
 Integer,          Intent(OUT) :: dimids(*)
 Integer,          Intent(OUT) :: ndims, xtype, natts
 Integer                       :: status

 End Function nf_inq_var
End Interface
!-------------------------------- nf_inq_vardimid -----------------------------
Interface
 Function nf_inq_vardimid(ncid, varid, dimids) RESULT (status)

 Integer, Intent(IN)  :: ncid, varid
 Integer, Intent(OUT) :: dimids(*)
 Integer              :: status

 End Function nf_inq_vardimid
End Interface
!-------------------------------- nf_inq_varid --------------------------------
Interface
 Function nf_inq_varid(ncid, name, varid) RESULT (status)

 Integer,          Intent(IN)  :: ncid
 Integer,          Intent(OUT) :: varid
 Character(LEN=*), Intent(IN)  :: name
 Integer                       :: status

 End Function nf_inq_varid
End Interface
!-------------------------------- nf_inq_varname ------------------------------
Interface
 Function nf_inq_varname (ncid, varid, name) RESULT (status)

 Integer,          Intent(IN)   :: ncid, varid
 Character(LEN=*), Intent(OUT)  :: name
 Integer                        :: status

 End Function nf_inq_varname
End Interface
!-------------------------------- nf_inq_vartype ------------------------------
Interface
 Function nf_inq_vartype(ncid, varid, xtype) RESULT(status)

 Integer, Intent(IN)  :: ncid, varid
 Integer, Intent(OUT) :: xtype
 Integer              :: status

 End Function nf_inq_vartype
End Interface
!-------------------------------- nf_inq_varnatts -----------------------------
Interface
 Function nf_inq_varnatts(ncid, varid, nvatts) RESULT(status)

 Integer, Intent(IN)  :: ncid, varid
 Integer, Intent(OUT) :: nvatts
 Integer              :: status

 End Function nf_inq_varnatts
End Interface
!-------------------------------- nf_rename_var -------------------------------
Interface
 Function nf_rename_var(ncid, varid, name) RESULT (status)

 Integer,          Intent(IN) :: ncid, varid
 Character(LEN=*), Intent(IN) :: name
 Integer                      :: status

 End Function nf_rename_var
End Interface
!-------------------------------- nf_copy_var ---------------------------------
Interface
 Function nf_copy_var(ncid_in, varid, ncid_out) RESULT(status)

 Integer, Intent(IN) :: ncid_in, varid, ncid_out
 Integer             :: status

 End Function nf_copy_var
End Interface

! General attribute functions

!-------------------------------- nf_inq_att ----------------------------------
Interface
 Function nf_inq_att(ncid, varid, name, xtype, nlen) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(OUT) :: nlen, xtype
 Character(LEN=*), Intent(IN)  :: name
 Integer                       :: status

 End Function nf_inq_att
End Interface
!-------------------------------- nf_inq_atttype ---------------------------
Interface
 Function nf_inq_atttype(ncid, varid, name, xtype) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(OUT) :: xtype
 Character(LEN=*), Intent(IN)  :: name
 Integer                       :: status

 End Function nf_inq_atttype
End Interface
!-------------------------------- nf_inq_attlen -------------------------------
Interface
 Function nf_inq_attlen(ncid, varid, name, nlen) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(OUT) :: nlen
 Character(LEN=*), Intent(IN)  :: name
 Integer                       :: status

 End Function nf_inq_attlen
End Interface
!-------------------------------- nf_inq_attid --------------------------------
Interface
 Function nf_inq_attid(ncid, varid, name, attnum) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(OUT) :: attnum
 Character(LEN=*), Intent(IN)  :: name
 Integer                       :: status

 End Function nf_inq_attid
End Interface
!-------------------------------- nf_inq_attname ------------------------------
Interface
 Function nf_inq_attname(ncid, varid, attnum, name) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid, attnum
 Character(LEN=*), Intent(OUT) :: name
 Integer                       :: status

 End Function nf_inq_attname
End Interface
!-------------------------------- nf_copy_att ---------------------------------
Interface
 Function nf_copy_att(ncid_in, varid_in, name, ncid_out, varid_out) &
                         RESULT(status)

 Integer,          Intent(IN)  :: ncid_in, varid_in, ncid_out, varid_out
 Character(LEN=*), Intent(IN)  :: name
 Integer                       :: status

 End Function nf_copy_att
End Interface
!-------------------------------- nf_rename_att -------------------------------
Interface
 Function nf_rename_att(ncid, varid, name, newname) RESULT(status)

 Integer,          Intent(IN) :: ncid, varid
 Character(LEN=*), Intent(IN) :: name, newname
 Integer                      :: status

 End Function nf_rename_att
End Interface
!-------------------------------- nf_del_att ----------------------------------
Interface
 Function nf_del_att(ncid, varid, name) RESULT(status)

 Integer,          Intent(IN) :: ncid, varid
 Character(LEN=*), Intent(IN) :: name
 Integer                      :: status

 End Function nf_del_att
End Interface

! var1 put and get functions

!--------------------------------- nf_put_var1_text ---------------------------
Interface
 Function nf_put_var1_text(ncid, varid, ndex, chval) RESULT(status)

 Integer,          Intent(IN) :: ncid, varid
 Integer,          Intent(IN) :: ndex(*)
 Character(LEN=1), Intent(IN) :: chval
 Integer                      :: status

 End Function nf_put_var1_text
End Interface
!--------------------------------- nf_put_var1_int1 ------------------------
Interface
 Function nf_put_var1_int1(ncid, varid, ndex, ival) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT1

 Integer,         Intent(IN) :: ncid, varid
 Integer,         Intent(IN) :: ndex(*)
 Integer(NFINT1), Intent(IN) :: ival
 Integer                     :: status

 End Function nf_put_var1_int1
End Interface
!--------------------------------- nf_put_var1_int2 ------------------------
Interface
 Function nf_put_var1_int2(ncid, varid, ndex, ival) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT2

 Integer,         Intent(IN) :: ncid, varid
 Integer,         Intent(IN) :: ndex(*)
 Integer(NFINT2), Intent(IN) :: ival
 Integer                     :: status

 End Function nf_put_var1_int2
End Interface
!--------------------------------- nf_put_var1_int -------------------------
Interface
 Function nf_put_var1_int(ncid, varid, ndex, ival) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT
 Integer,        Intent(IN) :: ncid, varid
 Integer,        Intent(IN) :: ndex(*)
 Integer(NFINT), Intent(IN) :: ival
 Integer                    :: status

 End Function nf_put_var1_int
End Interface
!--------------------------------- nf_put_var1_real ------------------------
Interface
 Function nf_put_var1_real(ncid, varid, ndex, rval) RESULT(status)

 USE netcdf_nf_data, ONLY: NFREAL

 Integer,      Intent(IN) :: ncid, varid
 Integer,      Intent(IN) :: ndex(*)
 Real(NFREAL), Intent(IN) :: rval
 Integer                  :: status

 End Function nf_put_var1_real
End Interface
!--------------------------------- nf_put_var1_double ----------------------
Interface
 Function nf_put_var1_double(ncid, varid, ndex, dval) RESULT(status)

 USE netcdf_nf_data, ONLY: RK8

 Integer,   Intent(IN) :: ncid, varid
 Integer,   Intent(IN) :: ndex(*)
 Real(RK8), Intent(IN) :: dval
 Integer               :: status

 End Function nf_put_var1_double
End Interface
!--------------------------------- nf_get_var1_text ------------------------
Interface
 Function nf_get_var1_text(ncid, varid, ndex, chval) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: ndex(*)
 Character(LEN=1), Intent(OUT) :: chval
 Integer                       :: status

 End Function nf_get_var1_text
End Interface
!--------------------------------- nf_get_var1_int1 ------------------------
Interface
 Function nf_get_var1_int1(ncid, varid, ndex, ival) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT1

 Integer,         Intent(IN)  :: ncid, varid
 Integer,         Intent(IN)  :: ndex(*)
 Integer(NFINT1), Intent(OUT) :: ival
 Integer                      :: status

 End Function nf_get_var1_int1
End Interface
!--------------------------------- nf_get_var1_int2 ------------------------
Interface
 Function nf_get_var1_int2(ncid, varid, ndex, ival) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT2

 Integer,         Intent(IN)  :: ncid, varid
 Integer,         Intent(IN)  :: ndex(*)
 Integer(NFINT2), Intent(OUT) :: ival
 Integer                      :: status

 End Function nf_get_var1_int2
End Interface
!--------------------------------- nf_get_var1_int -------------------------
Interface
 Function nf_get_var1_int(ncid, varid, ndex, ival) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT
 Integer,        Intent(IN)  :: ncid, varid
 Integer,        Intent(IN)  :: ndex(*)
 Integer(NFINT), Intent(OUT) :: ival
 Integer                     :: status

 End Function nf_get_var1_int
End Interface
!--------------------------------- nf_get_var1_real ------------------------
Interface
 Function nf_get_var1_real(ncid, varid, ndex, rval) RESULT(status)

 USE netcdf_nf_data, ONLY: NFREAL

 Integer,      Intent(IN)  :: ncid, varid
 Integer,      Intent(IN)  :: ndex(*)
 Real(NFREAL), Intent(OUT) :: rval
 Integer                   :: status

 End Function nf_get_var1_real
End Interface
!--------------------------------- nf_get_var1_double ----------------------
Interface
 Function nf_get_var1_double(ncid, varid, ndex, rval) RESULT(status)

 USE netcdf_nf_data, ONLY: RK8

 Integer,   Intent(IN)  :: ncid, varid
 Integer,   Intent(IN)  :: ndex(*)
 Real(RK8), Intent(OUT) :: rval
 Integer                :: status

 End Function nf_get_var1_double
End Interface

! var put and get functions

!--------------------------------- nf_put_var_text -------------------------
Interface nf_put_var_text
 Function nf_put_var_text(ncid, varid, text) RESULT(status)

 Integer,          Intent(IN) :: ncid, varid
 Character(LEN=*), Intent(IN) :: text
 Integer                      :: status

 End Function nf_put_var_text
! Array of characters
 Function nf_put_var_text_a(ncid, varid, text) RESULT(status)

 Integer,          Intent(IN) :: ncid, varid
 Character(LEN=1), Intent(IN) :: text(*)
 Integer                      :: status

 End Function nf_put_var_text_a
End Interface
!--------------------------------- nf_put_var_int1 -------------------------
Interface
 Function nf_put_var_int1(ncid, varid, i1vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT1

 Integer,         Intent(IN) :: ncid, varid
 Integer(NFINT1), Intent(IN) :: i1vals(*)
 Integer                     :: status

 End Function nf_put_var_int1
End Interface
!--------------------------------- nf_put_var_int2 -------------------------
Interface
 Function nf_put_var_int2(ncid, varid, i2vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT2

 Integer,         Intent(IN) :: ncid, varid
 Integer(NFINT2), Intent(IN) :: i2vals(*)
 Integer                     :: status

 End Function nf_put_var_int2
End Interface
!--------------------------------- nf_put_var_int --------------------------
Interface
 Function nf_put_var_int(ncid, varid, ivals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT

 Integer,        Intent(IN) :: ncid, varid
 Integer(NFINT), Intent(IN) :: ivals(*)
 Integer             :: status

 End Function nf_put_var_int
End Interface
!--------------------------------- nf_put_var_real -------------------------
Interface
 Function nf_put_var_real(ncid, varid, rvals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFREAL

 Integer,      Intent(IN) :: ncid, varid
 Real(NFREAL), Intent(IN) :: rvals(*)
 Integer                  :: status

 End Function nf_put_var_real
End Interface
!--------------------------------- nf_put_var_double -----------------------
Interface
 Function nf_put_var_double(ncid, varid, dvals) RESULT(status)

 USE netcdf_nf_data, ONLY: RK8

 Integer,   Intent(IN) :: ncid, varid
 Real(RK8), Intent(IN) :: dvals(*)
 Integer               :: status

 End Function nf_put_var_double
End Interface
!--------------------------------- nf_get_var_text ------------------------
Interface nf_get_var_text
 Function nf_get_var_text(ncid, varid, text) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(OUT) :: text
 Integer                       :: status

 End Function nf_get_var_text
! array of characters
 Function nf_get_var_text_a(ncid, varid, text) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=1), Intent(OUT) :: text(*)
 Integer                       :: status

 End Function nf_get_var_text_a
End Interface
!--------------------------------- nf_get_var_int1 -------------------------
Interface
 Function nf_get_var_int1(ncid, varid, i1vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT1

 Integer,         Intent(IN)  :: ncid, varid
 Integer(NFINT1), Intent(OUT) :: i1vals(*)
 Integer                      :: status

 End Function nf_get_var_int1
End Interface
!--------------------------------- nf_get_var_int2 -------------------------
Interface
 Function nf_get_var_int2(ncid, varid, i2vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT2

 Integer,         Intent(IN)  :: ncid, varid
 Integer(NFINT2), Intent(OUT) :: i2vals(*)
 Integer                      :: status

 End Function nf_get_var_int2
End Interface
!--------------------------------- nf_get_var_int --------------------------
Interface
 Function nf_get_var_int(ncid, varid, ivals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT
 
 Integer,        Intent(IN)  :: ncid, varid
 Integer(NFINT), Intent(OUT) :: ivals(*)
 Integer                     :: status

 End Function nf_get_var_int
End Interface
!--------------------------------- nf_get_var_real -------------------------
Interface
 Function nf_get_var_real(ncid, varid, rvals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFREAL

 Integer,      Intent(IN)  :: ncid, varid
 Real(NFREAL), Intent(OUT) :: rvals(*)
 Integer                   :: status

 End Function nf_get_var_real
End Interface
!--------------------------------- nf_get_var_double -----------------------
Interface
 Function nf_get_var_double(ncid, varid, dvals) RESULT(status)

 USE netcdf_nf_data, ONLY: RK8

 Integer,   Intent(IN)  :: ncid, varid
 Real(RK8), Intent(OUT) :: dvals(*)
 Integer                :: status

 End Function nf_get_var_double
End Interface

! vars put and get functions

!--------------------------------- nf_put_vars_text ------------------------
Interface nf_put_vars_text
 Function nf_put_vars_text(ncid, varid, start, counts, strides, text) &
                              RESULT(status)

 Integer,          Intent(IN) :: ncid, varid
 Integer,          Intent(IN) :: start(*), counts(*), strides(*)
 Character(LEN=*), Intent(IN) :: text
 Integer                      :: status

 End Function nf_put_vars_text
! array of characters
 Function nf_put_vars_text_a(ncid, varid, start, counts, strides, text) &
                              RESULT(status)

 Integer,          Intent(IN) :: ncid, varid
 Integer,          Intent(IN) :: start(*), counts(*), strides(*)
 Character(LEN=1), Intent(IN) :: text(*)
 Integer                      :: status

 End Function nf_put_vars_text_a
End Interface
!--------------------------------- nf_put_vars_int1 ------------------------
Interface
 Function nf_put_vars_int1(ncid, varid, start, counts, strides, i1vals) &
                              RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT1

 Integer,         Intent(IN) :: ncid, varid
 Integer,         Intent(IN) :: start(*), counts(*), strides(*)
 Integer(NFINT1), Intent(IN) :: i1vals(*)
 Integer                     :: status

 End Function nf_put_vars_int1
End Interface
!--------------------------------- nf_put_vars_int2 ------------------------
Interface
 Function nf_put_vars_int2(ncid, varid, start, counts, strides, i2vals) &
                              RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT2

 Integer,         Intent(IN) :: ncid, varid
 Integer,         Intent(IN) :: start(*), counts(*), strides(*)
 Integer(NFINT2), Intent(IN) :: i2vals(*)
 Integer                     :: status

 End Function nf_put_vars_int2
End Interface
!--------------------------------- nf_put_vars_int -------------------------
Interface
 Function nf_put_vars_int(ncid, varid, start, counts, strides, ivals) &
                             RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT

 Integer,        Intent(IN) :: ncid, varid
 Integer,        Intent(IN) :: start(*), counts(*), strides(*)
 Integer(NFINT), Intent(IN) :: ivals(*)
 Integer                    :: status

 End Function nf_put_vars_int
End Interface
!--------------------------------- nf_put_vars_real ------------------------
Interface
 Function nf_put_vars_real(ncid, varid, start, counts, strides, rvals) &
                              RESULT(status)

 USE netcdf_nf_data, ONLY: NFREAL

 Integer,      Intent(IN) :: ncid, varid
 Integer,      Intent(IN) :: start(*), counts(*), strides(*)
 Real(NFREAL), Intent(IN) :: rvals(*)
 Integer                  :: status

 End Function nf_put_vars_real
End Interface
!--------------------------------- nf_put_vars_double ----------------------
Interface
 Function nf_put_vars_double(ncid, varid, start, counts, strides, dvals) &
                                RESULT(status)

 USE netcdf_nf_data, ONLY: RK8

 Integer,   Intent(IN) :: ncid, varid
 Integer,   Intent(IN) :: start(*), counts(*), strides(*)
 Real(RK8), Intent(IN) :: dvals(*)
 Integer               :: status

 End Function nf_put_vars_double
End Interface
!--------------------------------- nf_get_vars_text ------------------------
Interface nf_get_vars_text
 Function nf_get_vars_text(ncid, varid, start, counts, strides, text) &
                                RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: start(*), counts(*), strides(*)
 Character(LEN=*), Intent(OUT) :: text
 Integer                       :: status

 End Function nf_get_vars_text
! array of characters
 Function nf_get_vars_text_a(ncid, varid, start, counts, strides, text) &
                                RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: start(*), counts(*), strides(*)
 Character(LEN=1), Intent(OUT) :: text(*)
 Integer                       :: status

 End Function nf_get_vars_text_a
End Interface
!--------------------------------- nf_get_vars_int1 ------------------------
Interface
 Function nf_get_vars_int1(ncid, varid, start, counts, strides, i1vals) &
                              RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT1

 Integer,         Intent(IN)  :: ncid, varid
 Integer,         Intent(IN)  :: start(*), counts(*), strides(*)
 Integer(NFINT1), Intent(OUT) :: i1vals(*)
 Integer                      :: status

 End Function nf_get_vars_int1
End Interface
!--------------------------------- nf_get_vars_int2 ------------------------
Interface
 Function nf_get_vars_int2(ncid, varid, start, counts, strides, i2vals) &
                              RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT2

 Integer,         Intent(IN)  :: ncid, varid
 Integer,         Intent(IN)  :: start(*), counts(*), strides(*)
 Integer(NFINT2), Intent(OUT) :: i2vals(*)
 Integer                      :: status

 End Function nf_get_vars_int2
End Interface
!--------------------------------- nf_get_vars_int -------------------------
Interface
 Function nf_get_vars_int(ncid, varid, start, counts, strides, ivals) &
                             RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT

 Integer,        Intent(IN)  :: ncid, varid
 Integer,        Intent(IN)  :: start(*), counts(*), strides(*)
 Integer(NFINT), Intent(OUT) :: ivals(*)
 Integer                     :: status

 End Function nf_get_vars_int
End Interface
!--------------------------------- nf_get_vars_real ------------------------
Interface
 Function nf_get_vars_real(ncid, varid, start, counts, strides, rvals) &
                              RESULT(status)

 USE netcdf_nf_data, ONLY: NFREAL

 Integer,      Intent(IN)  :: ncid, varid
 Integer,      Intent(IN)  :: start(*), counts(*), strides(*)
 Real(NFREAL), Intent(OUT) :: rvals(*)
 Integer                   :: status

 End Function nf_get_vars_real
End Interface
!--------------------------------- nf_get_vars_double ----------------------
Interface
 Function nf_get_vars_double(ncid, varid, start, counts, strides, dvals) &
                                RESULT(status)

 USE netcdf_nf_data, ONLY: RK8

 Integer,   Intent(IN)  :: ncid, varid
 Integer,   Intent(IN)  :: start(*), counts(*), strides(*)
 Real(RK8), Intent(OUT) :: dvals(*)
 Integer                :: status

 End Function nf_get_vars_double
End Interface

! varm put and get functions

!--------------------------------- nf_put_varm_text ------------------------
Interface nf_put_varm_text
 Function nf_put_varm_text(ncid, varid, start, counts, strides, maps, &
                                text) RESULT(status)



 Integer,          Intent(IN) :: ncid, varid
 Integer,          Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Character(LEN=*), Intent(IN) :: text
 Integer                      :: status

 End Function nf_put_varm_text
! array of characters
 Function nf_put_varm_text_a(ncid, varid, start, counts, strides, maps, &
                                text) RESULT(status)

 Integer,          Intent(IN) :: ncid, varid
 Integer,          Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Character(LEN=1), Intent(IN) :: text(*)
 Integer                      :: status

 End Function nf_put_varm_text_a
End Interface
!--------------------------------- nf_put_varm_int1 ------------------------
Interface
 Function nf_put_varm_int1(ncid, varid, start, counts, strides, maps, &
                              i1vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT1

 Integer,         Intent(IN) :: ncid, varid
 Integer,         Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Integer(NFINT1), Intent(IN) :: i1vals(*)
 Integer                     :: status

 End Function nf_put_varm_int1
End Interface
!--------------------------------- nf_put_varm_int2 ------------------------
Interface
 Function nf_put_varm_int2(ncid, varid, start, counts, strides, maps, &
                              i2vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT2

 Integer,         Intent(IN) :: ncid, varid
 Integer,         Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Integer(NFINT2), Intent(IN) :: i2vals(*)
 Integer                     :: status

 End Function nf_put_varm_int2
End Interface
!--------------------------------- nf_put_varm_int -------------------------
Interface
 Function nf_put_varm_int(ncid, varid, start, counts, strides, maps, &
                             ivals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT

 Integer,        Intent(IN) :: ncid, varid
 Integer,        Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Integer(NFINT), Intent(IN) :: ivals(*)
 Integer                    :: status

 End Function nf_put_varm_int
End Interface

!--------------------------------- nf_put_varm_real ------------------------
Interface
 Function nf_put_varm_real(ncid, varid, start, counts, strides, maps, &
                              rvals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFREAL

 Integer,      Intent(IN) :: ncid, varid
 Integer,      Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Real(NFREAL), Intent(IN) :: rvals(*)
 Integer                  :: status

 End Function nf_put_varm_real
End Interface
!--------------------------------- nf_put_varm_double ----------------------
Interface
 Function nf_put_varm_double(ncid, varid, start, counts, strides, maps, &
                                dvals) RESULT(status)

 USE netcdf_nf_data, ONLY: RK8

 Integer,   Intent(IN) :: ncid, varid
 Integer,   Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Real(RK8), Intent(IN) :: dvals(*)
 Integer               :: status

 End Function nf_put_varm_double
End Interface
!--------------------------------- nf_get_varm_text ------------------------
Interface nf_get_varm_text
 Function nf_get_varm_text(ncid, varid, start, counts, strides, maps, &
                                text) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Character(LEN=*), Intent(OUT) :: text
 Integer                       :: status

 End Function nf_get_varm_text
! array of characters
 Function nf_get_varm_text_a(ncid, varid, start, counts, strides, maps, &
                                text) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Character(LEN=1), Intent(OUT) :: text(*)
 Integer                       :: status

 End Function nf_get_varm_text_a
End Interface
!--------------------------------- nf_get_varm_int1 ------------------------
Interface
 Function nf_get_varm_int1(ncid, varid, start, counts, strides, maps, &
                              i1vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT1

 Integer,         Intent(IN)  :: ncid, varid
 Integer,         Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Integer(NFINT1), Intent(OUT) :: i1vals(*)
 Integer                      :: status

 End Function nf_get_varm_int1
End Interface
!--------------------------------- nf_get_varm_int2 ------------------------
Interface
 Function nf_get_varm_int2(ncid, varid, start, counts, strides, maps, &
                              i2vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT2

 Integer,         Intent(IN)  :: ncid, varid
 Integer,         Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Integer(NFINT2), Intent(OUT) :: i2vals(*)
 Integer                      :: status

 End Function nf_get_varm_int2
End Interface
!--------------------------------- nf_get_varm_int -------------------------
Interface
 Function nf_get_varm_int(ncid, varid, start, counts, strides, maps, &
                             ivals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT

 Integer,        Intent(IN)  :: ncid, varid
 Integer,        Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Integer(NFINT), Intent(OUT) :: ivals(*)
 Integer                     :: status

 End Function nf_get_varm_int
End Interface
!--------------------------------- nf_get_varm_real ------------------------
Interface
 Function nf_get_varm_real(ncid, varid, start, counts, strides, maps, &
                              rvals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFREAL

 Integer,      Intent(IN)  :: ncid, varid
 Integer,      Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Real(NFREAL), Intent(OUT) :: rvals(*)
 Integer                   :: status

 End Function nf_get_varm_real
End Interface
!--------------------------------- nf_get_varm_double ----------------------
Interface
 Function nf_get_varm_double(ncid, varid, start, counts, strides, maps, &
                             dvals) RESULT(status)

 USE netcdf_nf_data, ONLY: RK8

 Integer,   Intent(IN)  :: ncid, varid
 Integer,   Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Real(RK8), Intent(OUT) :: dvals(*)
 Integer                :: status

 End Function nf_get_varm_double
End Interface

! vara put and get routines

!--------------------------------- nf_put_vara_text ------------------------
Interface nf_put_vara_text
 Function nf_put_vara_text(ncid, varid, start, counts, text) RESULT(status)

 Integer,          Intent(IN) :: ncid, varid
 Integer,          Intent(IN) :: start(*), counts(*)
 Character(LEN=*), Intent(IN) :: text
 Integer                      :: status

 End Function nf_put_vara_text
! array of characters
 Function nf_put_vara_text_a(ncid, varid, start, counts, text) RESULT(status)

 Integer,          Intent(IN) :: ncid, varid
 Integer,          Intent(IN) :: start(*), counts(*)
 Character(LEN=1), Intent(IN) :: text(*)
 Integer                      :: status

 End Function nf_put_vara_text_a
End Interface
!--------------------------------- nf_put_vara_int1 ------------------------
Interface
 Function nf_put_vara_int1(ncid, varid, start, counts, i1vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT1

 Integer,         Intent(IN) :: ncid, varid
 Integer,         Intent(IN) :: start(*), counts(*)
 Integer(NFINT1), Intent(IN) :: i1vals(*)
 Integer                     :: status

 End Function nf_put_vara_int1
End Interface
!--------------------------------- nf_put_vara_int2 ------------------------
Interface
 Function nf_put_vara_int2(ncid, varid, start, counts, i2vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT2

 Integer,         Intent(IN) :: ncid, varid
 Integer,         Intent(IN) :: start(*), counts(*)
 Integer(NFINT2), Intent(IN) :: i2vals(*)
 Integer                     :: status

 End Function nf_put_vara_int2
End Interface
!--------------------------------- nf_put_vara_int -------------------------
Interface
 Function nf_put_vara_int(ncid, varid, start, counts, ivals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT

 Integer,        Intent(IN) :: ncid, varid
 Integer,        Intent(IN) :: start(*), counts(*)
 Integer(NFINT), Intent(IN) :: ivals(*)
 Integer                    :: status

 End Function nf_put_vara_int
End Interface
!--------------------------------- nf_put_vara_real ------------------------
Interface
 Function nf_put_vara_real(ncid, varid, start, counts, rvals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFREAL

 Integer,      Intent(IN) :: ncid, varid
 Integer,      Intent(IN) :: start(*), counts(*)
 Real(NFREAL), Intent(IN) :: rvals(*)
 Integer                  :: status

 End Function nf_put_vara_real
End Interface
!--------------------------------- nf_put_vara_double ----------------------
Interface
 Function nf_put_vara_double(ncid, varid, start, counts, dvals) &
                                RESULT(status)

 USE netcdf_nf_data, ONLY: RK8

 Integer,   Intent(IN) :: ncid, varid
 Integer,   Intent(IN) :: start(*), counts(*)
 Real(RK8), Intent(IN) :: dvals(*)
 Integer               :: status

 End Function nf_put_vara_double
End Interface
!--------------------------------- nf_get_vara_text ------------------------
Interface nf_get_vara_text
 Function nf_get_vara_text(ncid, varid, start, counts, text) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: start(*), counts(*)
 Character(LEN=*), Intent(OUT) :: text
 Integer                       :: status

 End Function nf_get_vara_text
! array of characters
 Function nf_get_vara_text_a(ncid, varid, start, counts, text) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: start(*), counts(*)
 Character(LEN=1), Intent(OUT) :: text(*)
 Integer                       :: status

 End Function nf_get_vara_text_a
End Interface
!--------------------------------- nf_get_vara_int1 ------------------------
Interface
 Function nf_get_vara_int1(ncid, varid, start, counts, i1vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT1

 Integer,         Intent(IN)  :: ncid, varid
 Integer,         Intent(IN)  :: start(*), counts(*)
 Integer(NFINT1), Intent(OUT) :: i1vals(*)
 Integer                      :: status

 End Function nf_get_vara_int1
End Interface
!--------------------------------- nf_get_vara_int2 ------------------------
Interface
 Function nf_get_vara_int2(ncid, varid, start, counts, i2vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT2

 Integer,         Intent(IN)  :: ncid, varid
 Integer,         Intent(IN)  :: start(*), counts(*)
 Integer(NFINT2), Intent(OUT) :: i2vals(*)
 Integer                      :: status

 End Function nf_get_vara_int2
End Interface
!--------------------------------- nf_get_vara_int -------------------------
Interface
 Function nf_get_vara_int(ncid, varid, start, counts, ivals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT

 Integer,        Intent(IN)  :: ncid, varid
 Integer,        Intent(IN)  :: start(*), counts(*)
 Integer(NFINT), Intent(OUT) :: ivals(*)
 Integer                     :: status

 End Function nf_get_vara_int
End Interface
!--------------------------------- nf_get_vara_real ------------------------
Interface
 Function nf_get_vara_real(ncid, varid, start, counts, rvals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFREAL

 Integer,      Intent(IN)  :: ncid, varid
 Integer,      Intent(IN)  :: start(*), counts(*)
 Real(NFREAL), Intent(OUT) :: rvals(*)
 Integer                   :: status

 End Function nf_get_vara_real
End Interface
!--------------------------------- nf_get_vara_double ----------------------
Interface
 Function nf_get_vara_double(ncid, varid, start, counts, dvals) &
                                RESULT(status)

 USE netcdf_nf_data, ONLY: RK8

 Integer,   Intent(IN)  :: ncid, varid
 Integer,   Intent(IN)  :: start(*), counts(*)
 Real(RK8), Intent(OUT) :: dvals(*)
 Integer                :: status

 End Function nf_get_vara_double
End Interface
!--------------------------------- nf_put_att_text -------------------------
Interface nf_put_att_text
 Function nf_put_att_text(ncid, varid, name, nlen, text) RESULT(status)

 Integer,          Intent(IN) :: ncid, varid, nlen
 Character(LEN=*), Intent(IN) :: name, text
 Integer                      :: status

 End Function nf_put_att_text
! array of characters
 Function nf_put_att_text_a(ncid, varid, name, nlen, text) RESULT(status)

 Integer,          Intent(IN) :: ncid, varid, nlen
 Character(LEN=*), Intent(IN) :: name
 Character(LEN=1), Intent(IN) :: text(*)
 Integer                      :: status

 End Function nf_put_att_text_a
End Interface
!--------------------------------- nf_put_att_int1 -------------------------
Interface
 Function nf_put_att_int1(ncid, varid, name, xtype, nlen, i1vals) &
                             RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT1

 Integer,          Intent(IN) :: ncid, varid, nlen, xtype
 Character(LEN=*), Intent(IN) :: name
 Integer(NFINT1),  Intent(IN) :: i1vals(*)
 Integer                      :: status

 End Function nf_put_att_int1
End Interface
!--------------------------------- nf_put_att_int2 -------------------------
Interface
 Function nf_put_att_int2(ncid, varid, name, xtype, nlen, i2vals) &
                             RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT2

 Integer,          Intent(IN) :: ncid, varid, nlen, xtype
 Character(LEN=*), Intent(IN) :: name
 Integer(NFINT2),  Intent(IN) :: i2vals(*)
 Integer                      :: status

 End Function nf_put_att_int2
End Interface
!--------------------------------- nf_put_att_int --------------------------
Interface
 Function nf_put_att_int(ncid, varid, name, xtype, nlen, ivals) &
                            RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT

 Integer,          Intent(IN) :: ncid, varid, nlen, xtype
 Character(LEN=*), Intent(IN) :: name
 Integer(NFINT),   Intent(IN) :: ivals(*)
 Integer                      :: status

 End Function nf_put_att_int
End Interface
!--------------------------------- nf_put_att_real -------------------------
Interface
 Function nf_put_att_real(ncid, varid, name, xtype, nlen, rvals) &
                             RESULT(status)

 USE netcdf_nf_data, ONLY: NFREAL

 Integer,          Intent(IN) :: ncid, varid, nlen, xtype
 Character(LEN=*), Intent(IN) :: name
 Real(NFREAL),     Intent(IN) :: rvals(*)
 Integer                      :: status

 End Function nf_put_att_real
End Interface
!--------------------------------- nf_put_att_double -----------------------
Interface
 Function nf_put_att_double(ncid, varid, name, xtype, nlen, dvals) &
                               RESULT(status)

 USE netcdf_nf_data, ONLY: RK8

 Integer,          Intent(IN) :: ncid, varid, nlen, xtype
 Character(LEN=*), Intent(IN) :: name
 Real(RK8),        Intent(IN) :: dvals(*)
 Integer                      :: status

 End Function nf_put_att_double
End Interface
!--------------------------------- nf_get_att_text -------------------------
Interface nf_get_att_text
 Function nf_get_att_text(ncid, varid, name, text) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Character(LEN=*), Intent(OUT) :: text
 Integer                       :: status

 End Function nf_get_att_text
! Array of characters
 Function nf_get_att_text_a(ncid, varid, name, text) RESULT(status)

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Character(LEN=1), Intent(OUT) :: text(*)
 Integer                       :: status

 End Function nf_get_att_text_a
End Interface
!--------------------------------- nf_get_att_int1 -------------------------
Interface
 Function nf_get_att_int1(ncid, varid, name, i1vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT1

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Integer(NFINT1),  Intent(OUT) :: i1vals(*)
 Integer                       :: status

 End Function nf_get_att_int1
End Interface
!--------------------------------- nf_get_att_int2 -------------------------
Interface
 Function nf_get_att_int2(ncid, varid, name, i2vals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT2

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Integer(NFINT2),  Intent(OUT) :: i2vals(*)
 Integer                       :: status

 End Function nf_get_att_int2
End Interface
!--------------------------------- nf_get_att_int --------------------------
Interface
 Function nf_get_att_int(ncid, varid, name, ivals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFINT

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Integer(NFINT),   Intent(OUT) :: ivals(*)
 Integer                       :: status

 End Function nf_get_att_int
End Interface
!--------------------------------- nf_get_att_real -------------------------
Interface
 Function nf_get_att_real(ncid, varid, name, rvals) RESULT(status)

 USE netcdf_nf_data, ONLY: NFREAL

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Real(NFREAL),     Intent(OUT) :: rvals(*)
 Integer                       :: status

 End Function nf_get_att_real
End Interface
!--------------------------------- nf_get_att_double -----------------------
Interface
 Function nf_get_att_double(ncid, varid, name, dvals) RESULT(status)

 USE netcdf_nf_data, ONLY: RK8


 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Real(RK8),        Intent(OUT) :: dvals(*)
 Integer                       :: status

 End Function nf_get_att_double
End Interface

! Externals for functions that use C_CHAR strings to pass data to void
! pointers

 Integer, External :: nf_put_var1
 Integer, External :: nf_get_var1
 Integer, External :: nf_put_vars
 Integer, External :: nf_get_vars
 Integer, External :: nf_put_vara
 Integer, External :: nf_get_vara

#ifndef NO_NETCDF_2
! External definitons for Netcdf2 functions 
 Integer, External :: nccre 
 Integer, External :: ncopn 
 Integer, External :: ncddef
 Integer, External :: ncdid
 Integer, External :: ncvdef
 Integer, External :: ncvid
 Integer, External :: nctlen
 Integer, External :: ncsfil
#endif
!--------------------------- End Module netcdf_nf_interfaces - ----------------

End Module netcdf_nf_interfaces
