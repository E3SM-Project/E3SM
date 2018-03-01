!------ Routines for defining and obtaining info about dataset variables ------

! Replacement for fort-genvar.c

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

! Version 1.: Sept. 2005 - Initial Cray X1 version
! Version 2.: May   2006 - Updated to support g95
! Version 3.: April 2009 - Updated for netCDF 4.0.1
! Version 4.: April 2010 - Updated for netCDF 4.1.1
! Version 5.: May   2014 - Ensure return error status checked from C API calls          
          
!-------------------------------- nf_def_var -------------------------------
 Function nf_def_var(ncid, name, xtype, nvdims, vdims, varid) RESULT (status)

! Define name, datatype, number of dimensions, and dimension ids for a
! dataset variable. Returns variable id

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, xtype, nvdims
 Integer,          Intent(IN)  :: vdims(*)
 Integer,          Intent(OUT) :: varid
 Character(LEN=*), Intent(IN)  :: name

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cnvdims, cvarid, cstatus, cxtype
 Integer(KIND=C_INT)          :: cvdims(NC_MAX_DIMS)
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid   = ncid
 cnvdims = nvdims
 cxtype  = xtype
 
! Check if a C NULL character was appended to name in calling routine

 cname = addCNullChar(name, ie)

! Flip dimids to C order and subtract -1 to yield C ids prior
! to calling nc_def_var. Replaces C utility f2c_dimids

 cvdims = 0
 If (nvdims /= 0) Then 
   cvdims(1:nvdims) = vdims(nvdims:1:-1)-1
 EndIf

 cstatus = nc_def_var(cncid, cname(1:ie+1), cxtype, &
                     cnvdims, cvdims, cvarid)

 If (cstatus == NC_NOERR) Then
    ! Add one to returned C varid to yield FORTRAN id
    varid = cvarid + 1
 EndIf
 status = cstatus

 End Function nf_def_var
!-------------------------------- nf_inq_varndims --------------------------
 Function nf_inq_varndims(ncid, varid, vndims) RESULT (status)

! Get variable dimension for a given varid from NetCDF dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid, varid
 Integer, Intent(OUT) :: vndims

 Integer              :: status

 Integer(KIND=C_INT) :: cncid, cvarid, cvndims, cstatus

 cncid   = ncid
 cvarid  = varid - 1

 cstatus = nc_inq_varndims(cncid, cvarid, cvndims)

 If (cstatus == NC_NOERR) Then
    vndims = cvndims
 EndIf
 status = cstatus

 End Function nf_inq_varndims
!-------------------------------- nf_inq_var ----------------------------
 Function nf_inq_var(ncid, varid, name, xtype, ndims, dimids, natts) &
                     RESULT (status)

! Get variable name, data type, dimension length, dimension ids, and
! number of attributes for a given varid 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(OUT) :: name
 Integer,          Intent(OUT) :: dimids(*)
 Integer,          Intent(OUT) :: ndims, xtype, natts

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cstatus, cndims, cvarid, cnatts
 Integer(KIND=C_INT)          :: cdimids(NC_MAX_DIMS)
 Integer(KIND=C_INT)          :: cxtype
 Character(LEN=NC_MAX_NAME+1) :: tmpname
 Integer                      :: nlen

 cncid  = ncid
 cvarid = varid - 1 ! subtract -1 to yield cvarid

 nlen    = LEN(name)
 tmpname = REPEAT(" ", LEN(tmpname))
 name    = REPEAT(" ", nlen)

 cndims    = 0
 dimids(1) = 0
 xtype     = 0
 natts     = 0
 ndims     = 0

 cstatus = nc_inq_var(cncid, cvarid, tmpname, cxtype, cndims, cdimids, cnatts)

 If (cstatus == NC_NOERR) Then
    xtype = cxtype
    natts = cnatts
    ndims = cndims

    ! Check tmpname for a C null character and strip it and trailing blanks

    name = stripCNullChar(tmpname, nlen)

    ! Reverse order of cdimids and add one to yield FORTRAN id numbers
    ! Replaces c2f_dimids C utility
 
    If (ndims > 0) Then
       dimids(1:ndims) = cdimids(ndims:1:-1)+1
    EndIf
 EndIf

 status = cstatus

 End Function nf_inq_var
!-------------------------------- nf_inq_vardimid -----------------------
 Function nf_inq_vardimid(ncid, varid, dimids) RESULT (status)

! Get variable dimension id for a dimension name

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid, varid
 Integer, Intent(OUT) :: dimids(*)

 Integer              :: status

 Integer(KIND=C_INT) :: cncid, cstatus, cstat2, cndims, cvarid
 Integer(KIND=C_INT) :: cvdimids(NC_MAX_DIMS)
 Integer             :: ndims

 cncid     = ncid
 cvarid    = varid - 1 ! subtract -1 to get C id number
 dimids(1) = 0
 cvdimids  = 0
 cndims    = 0
 ndims     = 0
 
 cstat2  = nc_inq_varndims(cncid, cvarid, cndims)
 cstatus = nc_inq_vardimid(cncid, cvarid, cvdimids)

! Reverse order of cdimids and add 1 to yield FORTRAN id numbers
! Replaces c2f_dimids C utility
 
 If (cstat2 == NC_NOERR .AND. cstatus == NC_NOERR) Then
   ndims = cndims
   If (ndims > 0) Then    
     dimids(1:ndims) = cvdimids(ndims:1:-1)+1
   EndIf
 Else
   ndims = 0
 EndIf
 
 status = cstatus

 End Function nf_inq_vardimid
!-------------------------------- nf_inq_varid --------------------------
 Function nf_inq_varid(ncid, name, varid) RESULT (status)

! Get variable id for a variable name from NetCDF dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: varid

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid = ncid

! Check name for a C NULL character added in calling routine

 cname = addCNullChar(name, ie)

 cstatus = nc_inq_varid(cncid, cname(1:ie+1), cvarid)

 If (cstatus == NC_NOERR) Then
    varid  = cvarid + 1  ! add one to get Fortran id number
 EndIf
 status = cstatus

 End Function nf_inq_varid
!-------------------------------- nf_inq_varname ------------------------
 Function nf_inq_varname (ncid, varid, name) RESULT (status)

! Get variable name for a given varid from NetCDF dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)   :: ncid, varid
 Character(LEN=*), Intent(OUT)  :: name

 Integer                        :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus
 Character(LEN=NC_MAX_NAME+1) :: tmpname 
 Integer                      :: nlen

 cncid   = ncid
 cvarid  = varid - 1 ! subtract one to get C id number

 nlen    = LEN(name)
 tmpname = REPEAT(" ", LEN(tmpname))
 name    = REPEAT(" ", nlen)

! Get tmpname from C interface

 cstatus = nc_inq_varname(cncid, cvarid, tmpname)

 If (cstatus == NC_NOERR) Then
    ! Find first C null character in tmpname if present and set end of string
    name = stripCNullChar(tmpname, nlen)
 EndIf
 status = cstatus

 End Function nf_inq_varname
!-------------------------------- nf_inq_vartype -------------------------
 Function nf_inq_vartype(ncid, varid, xtype) RESULT(status)

! Inquire about numeric type of variable attributes for NetCDF file id ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid, varid
 Integer, Intent(OUT) :: xtype

 Integer              :: status

 Integer(KIND=C_INT) :: cncid, cvarid, cstatus
 Integer(KIND=C_INT) :: cxtype

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C id number

 cstatus = nc_inq_vartype(cncid, cvarid, cxtype)

 If (cstatus == NC_NOERR) Then
    xtype  = cxtype
 EndIf
 status = cstatus

 End Function nf_inq_vartype
!-------------------------------- nf_inq_varnatts -----------------------
 Function nf_inq_varnatts(ncid, varid, nvatts) RESULT(status)

! Inquire about number of variable attributes for NetCDF file id ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid, varid
 Integer, Intent(OUT) :: nvatts

 Integer              :: status

 Integer(KIND=C_INT) :: cncid, cvarid, cnvatts, cstatus

 cncid  = ncid
 cvarid = varid - 1 ! subtract one to get C id number

 cstatus = nc_inq_varnatts(cncid, cvarid, cnvatts)

 If (cstatus == NC_NOERR) Then
    nvatts = cnvatts
 EndIf

 status = cstatus

 End Function nf_inq_varnatts
!-------------------------------- nf_rename_var -------------------------
 Function nf_rename_var(ncid, varid, name) RESULT (status)

! Rename dimension name for a given dimension id

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid 
 Character(LEN=*), Intent(IN) :: name

 Integer                      :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C id number

! Check name for a C NULL character added in calling routine

 cname = addCNullChar(name, ie)

 cstatus = nc_rename_var(cncid, cvarid, cname(1:ie+1))

 status = cstatus

 End Function nf_rename_var
!-------------------------------- nf_copy_var ---------------------------
 Function nf_copy_var(ncid_in, varid, ncid_out) RESULT(status)

! Copies a given variable from dataset ncid_in to dataset ncid_out

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN) :: ncid_in, varid, ncid_out
 Integer             :: status

 Integer(KIND=C_INT) :: cncidin, cvarid, cncidout, cstatus

 cncidin  = ncid_in
 cncidout = ncid_out
 cvarid   = varid - 1 ! Subtract one to get C id number

 cstatus = nc_copy_var(cncidin, cvarid, cncidout)

 status = cstatus

 End Function nf_copy_var
