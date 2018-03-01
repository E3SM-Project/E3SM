!---------- Routines for defining and obtaining info about attributes --------

! Replacement for fort-genatt.c

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
         
!-------------------------------- nf_inq_att ---------------------------------
 Function nf_inq_att(ncid, varid, name, xtype, nlen) RESULT(status)

! Get attribute data type and length for a given varid and name

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: nlen, xtype

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cstatus, cvarid
 Integer(KIND=C_SIZE_T)       :: cnlen
 Integer(KIND=C_INT)          :: cxtype
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

! Check to see if a C null character was added to name in calling program

 cname = addCNullChar(name, ie)

 cstatus = nc_inq_att(cncid, cvarid, cname(1:ie+1), cxtype, cnlen)

 If (cstatus == NC_NOERR) Then
    xtype = cxtype
    nlen  = cnlen
 EndIf
 status = cstatus

 End Function nf_inq_att
!-------------------------------- nf_inq_atttype ---------------------------
 Function nf_inq_atttype(ncid, varid, name, xtype) RESULT(status)

! Get attribute type for a given varid and name

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: xtype

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cstatus, cvarid
 Integer(KIND=C_INT)          :: cxtype
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

! Check to see if a C null character was added to name in calling program

 cname = addCNullChar(name, ie)

 cstatus = nc_inq_atttype(cncid, cvarid, cname(1:ie+1), cxtype)

 If (cstatus == NC_NOERR) Then
    xtype = cxtype
 EndIf
 status = cstatus

 End Function nf_inq_atttype
!-------------------------------- nf_inq_attlen ----------------------------
 Function nf_inq_attlen(ncid, varid, name, nlen) RESULT(status)

! Get attribute length for a given varid and name

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: nlen

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cstatus, cvarid
 Integer(KIND=C_SIZE_T)       :: cnlen
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

! Check to see if a C null character was added to name in calling program

 cname = addCNullChar(name, ie)

 cstatus = nc_inq_attlen(cncid, cvarid, cname(1:ie+1), cnlen)

 If (cstatus == NC_NOERR) Then
    nlen = cnlen
 EndIf
 status = cstatus

 End Function nf_inq_attlen
!-------------------------------- nf_inq_attid -----------------------------
 Function nf_inq_attid(ncid, varid, name, attnum) RESULT(status)

! Get attribute id for a given varid and name

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Integer,          Intent(OUT) :: attnum

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cstatus, cattnum, cvarid
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

! Check to see if a C null character was added to name in calling program

 cname = addCNullChar(name, ie)

 cstatus = nc_inq_attid(cncid, cvarid, cname(1:ie+1), cattnum)
 
 If (cstatus == NC_NOERR) Then
    attnum = cattnum + 1 ! add 1 to get FORTRAN att id
 EndIf
 status = cstatus

 End Function nf_inq_attid
!-------------------------------- nf_inq_attname ---------------------------
 Function nf_inq_attname(ncid, varid, attnum, name) RESULT(status)

! Get attribute name for a given varid and attribute number

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid, attnum
 Character(LEN=*), Intent(OUT) :: name

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cstatus, cattnum, cvarid
 Character(LEN=(LEN(name)+1)) :: tmpname
 Integer                      :: nlen

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid and att num
 cattnum = attnum - 1
 nlen    = LEN(name)
 name    = REPEAT(" ",nlen)
 tmpname = REPEAT(" ",LEN(tmpname)) ! init to blanks

 cstatus = nc_inq_attname(cncid, cvarid, cattnum, tmpname)

 If (cstatus == NC_NOERR) Then
    ! Strip of any C null characters and load only the part
    ! of tmpname that will fit in name

    name = stripCNullChar(tmpname, nlen)
 EndIf
 status = cstatus

 End Function nf_inq_attname
!-------------------------------- nf_copy_att ------------------------------
 Function nf_copy_att(ncid_in, varid_in, name, ncid_out, varid_out) &
                      RESULT(status)

! Copy attribute name with varid_in from one netcdf file to another
! with new varid_out 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid_in, varid_in, ncid_out, varid_out 
 Character(LEN=*), Intent(IN) :: name

 Integer                      :: status

 Integer(KIND=C_INT)          :: cncidin, cncidout,cvaridin, cvaridout, cstatus
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncidin   = ncid_in
 cvaridin  = varid_in - 1
 cncidout  = ncid_out
 cvaridout = varid_out - 1

! Check to see if a C null character was added to name in calling program

 cname = addCNullChar(name, ie)

 cstatus = nc_copy_att(cncidin, cvaridin, cname(1:ie+1), &
                       cncidout, cvaridout)

 status = cstatus

 End Function nf_copy_att
!-------------------------------- nf_rename_att ----------------------------
 Function nf_rename_att(ncid, varid, name, newname) RESULT(status)

! Rename an attribute to newname givin varid 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid
 Character(LEN=*), Intent(IN) :: name, newname

 Integer                      :: status

 Integer(KIND=C_INT)             :: cncid, cvarid, cstatus
 Character(LEN=(LEN(name)+1))    :: cname 
 Character(LEN=(LEN(newname)+1)) :: cnewname 
 Integer                         :: ie1, ie2, inull

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

! Check to see if a C null character was added to name and newname 
! in calling program

 cname = addCNullChar(name, ie1)

 cnewname = addCNullChar(newname, ie2)

 cstatus = nc_rename_att(cncid, cvarid, cname(1:ie1+1), cnewname(1:ie2+1))

 status = cstatus

 End Function nf_rename_att
!-------------------------------- nf_del_att -------------------------------
 Function nf_del_att(ncid, varid, name) RESULT(status)

! Delete an attribute givne varid and name 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid
 Character(LEN=*), Intent(IN) :: name

 Integer                      :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

! Check to see if a C null character was added to name in calling program

 cname = addCNullChar(name, ie)

 cstatus = nc_del_att(cncid, cvarid, cname(1:ie+1))

 status = cstatus

 End Function nf_del_att
