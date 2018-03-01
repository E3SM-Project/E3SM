!------- Routines for defining, obtaining, etc. global dimension info ------ 

! Replacement for fort-dim.c

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
          
!-------------------------------- nf_def_dim -------------------------------
 Function nf_def_dim(ncid, name, dlen, dimid) RESULT (status)

! Adds new dimensions to the NetCDF dataset given dimension name,
! and length. Returns dimension id

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, dlen
 Integer,          Intent(OUT) :: dimid
 Character(LEN=*), Intent(IN)  :: name

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cdimid, cstatus
 Integer(KIND=C_SIZE_T)       :: cdlen
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid = ncid
 cdlen = dlen

 dimid  = -1
 cdimid = -1

! Check to see if a C null character was appended in FORTRAN

 cname = addCNullChar(name, ie)
 
 cstatus = nc_def_dim(cncid, cname(1:ie+1), cdlen, cdimid)

 If (cstatus == NC_EBADDIM) Then  ! Return dimid=-1
   dimid = -1
 Else                    ! Add 1 to get FORTRAN dimid
   dimid = cdimid+1
 EndIf
 status = cstatus

 End Function nf_def_dim
!-------------------------------- nf_inq_dim -------------------------------
 Function nf_inq_dim(ncid, dimid, name, dlen) RESULT (status)

! Get dimension name and length for a given dimid from NetCDF dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, dimid
 Integer,          Intent(OUT) :: dlen
 Character(LEN=*), Intent(OUT) :: name
 
 Integer                       :: status

 Integer(KIND=C_INT)        :: cncid, cdimid, cstatus
 Integer(KIND=C_SIZE_T)     :: cdlen
 Integer                    :: nlen
 Character(LEN=NC_MAX_NAME) :: tmpname 

 cncid   = ncid
 cdimid  = dimid - 1   ! Subtract 1 to get C dimid
 tmpname = REPEAT(" ", LEN(tmpname))
 name    = REPEAT(" ", LEN(name))
 nlen    = LEN(name)

! Get tmpname and cdlen from C interface

 cstatus = nc_inq_dim(cncid, cdimid, tmpname, cdlen)

 If (cstatus == NC_NOERR) Then
    ! Strip C null char from tmpname if present and set end of string
    name = stripCNullChar(tmpname, nlen)
    dlen   = cdlen
 Endif

 status = cstatus

 End Function nf_inq_dim
!-------------------------------- nf_inq_dimid -----------------------------
 Function nf_inq_dimid(ncid, name, dimid) RESULT (status)

! Get dimension id for a given dimension name from dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid
 Integer,          Intent(OUT) :: dimid
 Character(LEN=*), Intent(IN)  :: name

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cdimid, cstatus
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid = ncid
 dimid  =  0 
 cdimid = -1

! Check to see if a C null character was appended in FORTRAN

 cname = addCNullChar(name, ie)
 
 cstatus = nc_inq_dimid(cncid, cname(1:ie+1), cdimid)

! add one to get FORTRAN dimid if not = -1

 If (cstatus == NC_EBADDIM) Then
   dimid = -1
 Else
   dimid = cdimid + 1
 EndIf
 status = cstatus

 End Function nf_inq_dimid
!-------------------------------- nf_inq_dimlen ----------------------------
 Function nf_inq_dimlen(ncid, dimid, dlen) RESULT (status)

! Get dimension length for a given dimid from NetCDF dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid, dimid
 Integer, Intent(OUT) :: dlen

 Integer              :: status

 Integer(KIND=C_INT)    :: cncid, cdimid, cstatus
 Integer(KIND=C_SIZE_T) :: cdlen

 cncid   = ncid
 cdimid  = dimid - 1 ! Subtract 1 to get C dimid
 dlen    = 0

 cstatus = nc_inq_dimlen(cncid, cdimid, cdlen)

 If (cstatus == NC_NOERR) Then
    dlen   = cdlen
 Endif
 status = cstatus

 End Function nf_inq_dimlen
!-------------------------------- nf_inq_dimname ---------------------------
 Function nf_inq_dimname (ncid, dimid, name) RESULT (status)

! Get dimension name for a given dimid from NetCDF dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, dimid
 Character(LEN=*), Intent(OUT) :: name

 Integer                       :: status

 Integer(KIND=C_INT)        :: cncid, cdimid, cstatus
 Integer                    :: nlen
 Character(LEN=NC_MAX_NAME) :: tmpname 

 cncid   = ncid
 cdimid  = dimid - 1 ! Subtract 1 to get C dimid
 tmpname = REPEAT(" ", LEN(tmpname))
 name    = REPEAT(" ", LEN(name))
 nlen    = LEN(name)

! Get tmpname and cdlen from C interface

 cstatus = nc_inq_dimname(cncid, cdimid, tmpname)

 If (cstatus == NC_NOERR) Then
    ! Strip C null character in tmpname if present and set end of string
    name = stripCNullChar(tmpname, nlen)
 Endif

 status = cstatus

 End Function nf_inq_dimname
!-------------------------------- nf_rename_dim ----------------------------
 Function nf_rename_dim(ncid, dimid, name) RESULT (status)

! Rename dimension name for a given dimension id

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, dimid
 Character(LEN=*), Intent(IN) :: name

 Integer                      :: status

 Integer(KIND=C_INT)          :: cncid, cdimid, cstatus
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cdimid = dimid - 1 ! Subtract 1 to get C dimid

! Check to see if a C null character was appended in FORTRAN

 cname = addCNullChar(name, ie)

 cstatus = nc_rename_dim(cncid, cdimid, cname(1:ie+1))

 status = cstatus

 End Function nf_rename_dim
