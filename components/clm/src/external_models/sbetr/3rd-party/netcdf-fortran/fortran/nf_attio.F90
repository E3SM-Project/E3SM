#include "nfconfig.inc"
!---------- Routines to put/get attribute data of various data types ----------

! Replacement for fort-attio.c

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
! The author grants to the University Center for Atmospheric Research
! (UCAR), Boulder, CO, USA the right to revise and extend the software
! without restriction. However, the author retains all copyrights and
! intellectual property rights explicitly stated in or implied by the
! Apache license

! Version 1.: Sept. 2005 - Initial Cray X1 version
! Version 2.: May 2006   - Updated to support g95
! Version 3.: April 2009 - Updated to Netcdf 4.0.1 
! Version 4.: April 2010 - Updated to Netcdf 4.1.1 
! Version 5.: Feb.  2013 - bug fixes for fortran 4.4
          
!--------------------------------- nf_put_att_text ---------------------------
 Function nf_put_att_text(ncid, varid, name, nlen, text) RESULT(status)

! Write variable or global attribute text string to dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid, nlen
 Character(LEN=*), Intent(IN) :: name, text

 Integer                      :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus
 Integer(KIND=C_SIZE_T)       :: cnlen
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid
 cnlen  = nlen

 cname = addCNullChar(name, ie)
 
 cstatus = nc_put_att_text(cncid, cvarid, cname(1:ie+1), cnlen, &
           text)

 status = cstatus

 End Function nf_put_att_text
!--------------------------------- nf_put_att_text_a ------------------------
 Function nf_put_att_text_a(ncid, varid, name, nlen, text) RESULT(status)

! New routine to support passing an array of single characters
! Write variable or global attribute array of characters to dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid, nlen
 Character(LEN=*), Intent(IN) :: name
 Character(LEN=1), Intent(IN) :: text(*)

 Integer                      :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus
 Integer(KIND=C_SIZE_T)       :: cnlen
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid
 cnlen  = nlen

 cname = addCNullChar(name, ie)
 
 cstatus = nc_put_att_text(cncid, cvarid, cname(1:ie+1), cnlen, &
                           text)

 status = cstatus

 End Function nf_put_att_text_a
!--------------------------------- nf_put_att_int1 -------------------------
 Function nf_put_att_int1(ncid, varid, name, xtype, nlen, i1vals) &
                          RESULT(status)

! Write variable or global attribute byte data to dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN) :: ncid, varid, nlen, xtype

 Character(LEN=*),     Intent(IN) :: name
 Integer(KIND=NFINT1), Intent(IN) :: i1vals(*)

 Integer                          :: status

 Integer(KIND=C_INT)             :: cncid, cvarid, cstatus, cxtype
 Integer(KIND=C_SIZE_T)          :: cnlen
 Character(LEN=(LEN(name)+1))    :: cname
 Integer                         :: ie

 If (C_SIGNED_CHAR < 0) Then ! schar not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf
 
 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid
 cnlen  = nlen
 cxtype = xtype

! Check for C null char on name and add one
 
 cname = addCNullChar(name, ie)

#if NF_INT1_IS_C_SIGNED_CHAR 
 cstatus = nc_put_att_schar(cncid, cvarid, cname(1:ie+1), &
                           cxtype, cnlen, i1vals) 
#elif NF_INT1_IS_C_SHORT
 cstatus = nc_put_att_short(cncid, cvarid, cname(1:ie+1), &
                            cxtype, cnlen, i1vals)
#elif NF_INT1_IS_C_INT
 cstatus = nc_put_att_int(cncid, cvarid, cname(1:ie+1), &
                            cxtype, cnlen, i1vals)
#elif NF_INT1_IS_C_LONG
 cstatus = nc_put_att_long(cncid, cvarid, cname(1:ie+1), &
                           cxtype, cnlen, i1vals)
#endif
 status = cstatus

 End Function nf_put_att_int1
!--------------------------------- nf_put_att_int2 -------------------------
 Function nf_put_att_int2(ncid, varid, name, xtype, nlen, i2vals) &
                          RESULT(status)

! Write variable or global attribute 16 bit integer data to dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN) :: ncid, varid, nlen, xtype

 Character(LEN=*),     Intent(IN) :: name
 Integer(KIND=NFINT2), Intent(IN) :: i2vals(*)

 Integer                          :: status

 Integer(KIND=C_INT)             :: cncid, cvarid, cstatus, cxtype
 Integer(KIND=C_SIZE_T)          :: cnlen
 Character(LEN=(LEN(name)+1))    :: cname
 Integer                         :: ie

 If (C_SHORT < 0) Then ! short not supported by processor
   status = NC_EBADTYPE
   Return
 EndIf
 
 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid
 cnlen  = nlen
 cxtype = xtype

 cname = addCNullChar(name, ie)

#if NF_INT2_IS_C_SHORT 
 cstatus = nc_put_att_short(cncid, cvarid, cname(1:ie+1), &
                            cxtype, cnlen, i2vals) 
#elif NF_INT2_IS_C_INT 
 cstatus = nc_put_att_int(cncid, cvarid, cname(1:ie+1), &
                           cxtype, cnlen, i2vals)
#elif NF_INT2_IS_C_LONG 
 cstatus = nc_put_att_long(cncid, cvarid, cname(1:ie+1), &
                           cxtype, cnlen, i2vals)
#endif 
 status = cstatus

 End Function nf_put_att_int2
!--------------------------------- nf_put_att_int --------------------------
 Function nf_put_att_int(ncid, varid, name, xtype, nlen, ivals) &
                         RESULT(status)

! Write variable or global attribute default integer data to dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid, nlen, xtype

 Character(LEN=*), Intent(IN) :: name
 Integer(NFINT),   Intent(IN) :: ivals(*)

 Integer                      :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus, cxtype
 Integer(KIND=C_SIZE_T)       :: cnlen
 Character(LEN=(LEN(name)+1)) :: cname
 Integer :: ie

 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid
 cnlen  = nlen
 cxtype = xtype

! Check for C null char and add one if missing

 cname = addCNullChar(name, ie)

#if NF_INT_IS_C_INT 
 cstatus = nc_put_att_int(cncid, cvarid, cname(1:ie+1), &
                          cxtype, cnlen, ivals) 
#elif NF_INT_IS_C_LONG 
 cstatus = nc_put_att_long(cncid, cvarid, cname(1:ie+1), &
                          cxtype, cnlen, ivals) 
#endif
 status = cstatus

 End Function nf_put_att_int
!--------------------------------- nf_put_att_real -------------------------
 Function nf_put_att_real(ncid, varid, name, xtype, nlen, rvals) &
                          RESULT(status)

! Write variable or global attribute Real(RK4) data to dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid, nlen, xtype

 Character(LEN=*), Intent(IN) :: name
 Real(NFREAL),     Intent(IN) :: rvals(*)

 Integer                      :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus, cxtype
 Integer(KIND=C_SIZE_T)       :: cnlen
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid
 cnlen  = nlen
 cxtype = xtype

! Check for C null char and add one if missing
 
 cname = addCNullChar(name, ie)

#if NF_REAL_IS_C_DOUBLE 
 cstatus = nc_put_att_double(cncid, cvarid, cname(1:ie+1), &
                             cxtype, cnlen, rvals) 
#else
 cstatus = nc_put_att_float(cncid, cvarid, cname(1:ie+1), &
                            cxtype, cnlen, rvals) 
#endif
 status = cstatus

 End Function nf_put_att_real
!--------------------------------- nf_put_att_double -----------------------
 Function nf_put_att_double(ncid, varid, name, xtype, nlen, dvals) &
                               RESULT(status)

! Write variable or global attribute Real(RK8) to dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid, nlen, xtype

 Character(LEN=*), Intent(IN) :: name
 Real(RK8),        Intent(IN) :: dvals(*)

 Integer                      :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus, cxtype
 Integer(KIND=C_SIZE_T)       :: cnlen
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid
 cnlen  = nlen
 cxtype = xtype

! Check for C null char and add one if missing

 cname = addCNullChar(name, ie)
 
 cstatus = nc_put_att_double(cncid, cvarid, cname(1:ie+1), &
                             cxtype, cnlen, dvals) 

 status = cstatus

 End Function nf_put_att_double
!--------------------------------- nf_get_att_text -----------------------
 Function nf_get_att_text(ncid, varid, name, text) RESULT(status)

! Read variable or global attribute character string from dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Character(LEN=*), Intent(OUT) :: text

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid
 text   = REPEAT(" ", LEN(text))

! Check for C null char and add one if missing

 cname = addCNullChar(name, ie)
 
 cstatus = nc_get_att_text(cncid, cvarid, cname(1:ie+1), text)

 status = cstatus

 End Function nf_get_att_text
!--------------------------------- nf_get_att_text_a -----------------------
 Function nf_get_att_text_a(ncid, varid, name, text) RESULT(status)

! New routine to support passing an array of single characters
! Read variable or global attribute array of characters from dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Character(LEN=1), Intent(OUT) :: text(*)

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid

! Check for C null char and add one if missing

 cname = addCNullChar(name, ie)
 
 cstatus = nc_get_att_text(cncid, cvarid, cname(1:ie+1), text)

 status = cstatus

 End Function nf_get_att_text_a
!--------------------------------- nf_get_att_int1 -------------------------
 Function nf_get_att_int1(ncid, varid, name, i1vals) RESULT(status)

! Read variable or global attribute BYTE integer data from dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN)  :: ncid, varid
 Character(LEN=*),     Intent(IN)  :: name
 Integer(KIND=NFINT1), Intent(OUT) :: i1vals(*)

 Integer                          :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 If (C_SIGNED_CHAR < 0) Then ! schar not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf
 
 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid

! Check for C null char and add one if missing

 cname = addCNullChar(name, ie)

#if NF_INT1_IS_C_SIGNED_CHAR 
 cstatus = nc_get_att_schar(cncid, cvarid, cname(1:ie+1), i1vals)
#elif NF_INT1_IS_C_SHORT
 cstatus = nc_get_att_short(cncid, cvarid, cname(1:ie+1), i1vals)
#elif NF_INT1_IS_C_INT
 cstatus = nc_get_att_int(cncid, cvarid, cname(1:ie+1), i1vals)
#elif NF_INT1_IS_C_LONG
 cstatus = nc_get_att_long(cncid, cvarid, cname(1:ie+1), i1vals)
#endif
 status = cstatus

 End Function nf_get_att_int1
!--------------------------------- nf_get_att_int2 --------------------------
 Function nf_get_att_int2(ncid, varid, name, i2vals) RESULT(status)

! Read variable or global attribute 16 bit integer data from dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN)  :: ncid, varid
 Character(LEN=*),     Intent(IN)  :: name
 Integer(KIND=NFINT2), Intent(OUT) :: i2vals(*)

 Integer                           :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 If (C_SHORT < 0) Then ! short not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf
 
 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid

! Check for C null char and add one if missing

 cname = addCNullChar(name, ie)
 
#if NF_INT2_IS_C_SHORT
 cstatus = nc_get_att_short(cncid, cvarid, cname(1:ie+1), i2vals) 
#elif NF_INT2_IS_C_INT
 cstatus = nc_get_att_int(cncid, cvarid, cname(1:ie+1), i2vals) 
#elif NF_INT2_IS_C_LONG
 cstatus = nc_get_att_long(cncid, cvarid, cname(1:ie+1), i2vals) 
#endif
 status = cstatus

 End Function nf_get_att_int2
!--------------------------------- nf_get_att_int ---------------------------
 Function nf_get_att_int(ncid, varid, name, ivals) RESULT(status)

! Read variable or global attribute default Integer data from dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Integer(NFINT),   Intent(OUT) :: ivals(*)

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid

! Check for C null char and add one if missing

 cname = addCNullChar(name, ie)

#if NF_INT_IS_C_INT 
 cstatus = nc_get_att_int(cncid, cvarid, cname(1:ie+1), ivals)
#elif NF_INT_IS_C_LONG
 cstatus = nc_get_att_long(cncid, cvarid, cname(1:ie+1), ivals)
#endif
 status = cstatus

 End Function nf_get_att_int
!--------------------------------- nf_get_att_real -------------------------
 Function nf_get_att_real(ncid, varid, name, rvals) RESULT(status)

! Read variable or global attribute Real(RK4) data from dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Real(NFREAL),     Intent(OUT) :: rvals(*)

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid

! Check for C null char and add one if missing

 cname = addCNullChar(name, ie)

#if NF_REAL_IS_C_DOUBLE 
 cstatus = nc_get_att_double(cncid, cvarid, cname(1:ie+1), rvals) 
#else
 cstatus = nc_get_att_float(cncid, cvarid, cname(1:ie+1), rvals) 
#endif
 status = cstatus

 End Function nf_get_att_real
!--------------------------------- nf_get_att_double -----------------------
 Function nf_get_att_double(ncid, varid, name, dvals) RESULT(status)

! Read variable or global attribute Real(RK8) data from dataset ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(IN)  :: name
 Real(RK8),        Intent(OUT) :: dvals(*)

 Integer                       :: status

 Integer(KIND=C_INT)          :: cncid, cvarid, cstatus
 Character(LEN=(LEN(name)+1)) :: cname
 Integer                      :: ie

 cncid  = ncid
 cvarid = varid -1 ! Subtract 1 to get C varid

! Check for C null char and add one if missing

 cname = addCNullChar(name, ie)
 
 cstatus = nc_get_att_double(cncid, cvarid, cname(1:ie+1), dvals)

 status = cstatus

 End Function nf_get_att_double
