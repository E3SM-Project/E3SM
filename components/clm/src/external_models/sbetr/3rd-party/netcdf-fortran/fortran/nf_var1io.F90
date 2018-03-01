#include "nfconfig.inc"
!----- Routines to put/get single data items of a variety of data types ------

! Replacement for fort-var1io.c

! Written by: Richard Weed, Ph.D
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
!                          Updated to pass ndex as C_PTR variable
! Version 3.: April 2009 - Updated for netCDF 4.0.1
! Version 4.: April 2010 - Updated for netCDF 4.1.1
!                          Added preprocessor test for int and real types
          
!--------------------------------- nf_put_var1_text ------------------------
 Function nf_put_var1_text(ncid, varid, ndex, chval) RESULT(status)

! Write out a single character variable to location vector ndex in dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid
 Integer,          Intent(IN) :: ndex(*)
 Character(LEN=1), Intent(IN) :: chval

 Integer                      :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Integer :: ndims

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cndex  = 0
 
 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! reverse array order and subtract 1 to get C index 
     cndex(1:ndims) = ndex(ndims:1:-1)-1
   EndIf
   cndexptr = C_LOC(cndex)
 EndIf

 cstatus = nc_put_var1_text(cncid, cvarid, cndexptr, chval)

 status = cstatus

 End Function nf_put_var1_text
!--------------------------------- nf_put_var1_int1 ------------------------
 Function nf_put_var1_int1(ncid, varid, ndex, ival) RESULT(status)

! Write out a 8 bit integer variable to location vector ndex in dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN) :: ncid, varid
 Integer,              Intent(IN) :: ndex(*)
 Integer(KIND=NFINT1), Intent(IN) :: ival

 Integer                          :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Integer(KIND=CINT1)            :: cival
 Integer                        :: ndims

 If (C_SIGNED_CHAR < 0) Then ! schar not supported by processor exit
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cival  = ival
 cndex  = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! reverse array order and subtract 1 to get C index 
     cndex(1:ndims) = ndex(ndims:1:-1)-1
   EndIf
   cndexptr = C_LOC(cndex)
 EndIf

#if NF_INT1_IS_C_SIGNED_CHAR
 cstatus = nc_put_var1_schar(cncid, cvarid, cndexptr, cival)
#elif NF_INT1_IS_C_SHORT
 cstatus = nc_put_var1_short(cncid, cvarid, cndexptr, cival)
#elif NF_INT1_IS_C_INT
 cstatus = nc_put_var1_int(cncid, cvarid, cndexptr, cival)
#elif NF_INT1_IS_C_LONG
 cstatus = nc_put_var1_long(cncid, cvarid, cndexptr, cival)
#endif

 status = cstatus

 End Function nf_put_var1_int1
!--------------------------------- nf_put_var1_int2 ------------------------
 Function nf_put_var1_int2(ncid, varid, ndex, ival) RESULT(status)

! Write out a 16 bit integer variable to location vector ndex in dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN) :: ncid, varid
 Integer,              Intent(IN) :: ndex(*)
 Integer(KIND=NFINT2), Intent(IN) :: ival

 Integer                          :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Integer(KIND=CINT2)            :: cival
 Integer                        :: ndims

 If (C_SHORT < 0) Then ! short not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cival  = ival
 cndex  = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! reverse array order and subtract 1 to get C index 
     cndex(1:ndims) = ndex(ndims:1:-1)-1
   EndIf
   cndexptr = C_LOC(cndex)
 EndIf

#if NF_INT2_IS_C_SHORT
 cstatus = nc_put_var1_short(cncid, cvarid, cndexptr, cival)
#elif NF_INT2_IS_C_INT
 cstatus = nc_put_var1_int(cncid, cvarid, cndexptr, cival)
#elif NF_INT2_IS_C_LONG
 cstatus = nc_put_var1_long(cncid, cvarid, cndexptr, cival)
#endif

 status = cstatus

 End Function nf_put_var1_int2
!--------------------------------- nf_put_var1_int -------------------------
 Function nf_put_var1_int(ncid, varid, ndex, ival) RESULT(status)

! Write out a default integer variable to location vector ndex to dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,        Intent(IN) :: ncid, varid
 Integer,        Intent(IN) :: ndex(*)
 Integer(NFINT), Intent(IN) :: ival

 Integer                    :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Integer(KIND=CINT)             :: cival
 Integer                        :: ndims

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cndex  = 0
 cival  = ival

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! reverse array order and subtract 1 to get C index 
     cndex(1:ndims) = ndex(ndims:1:-1)-1
   EndIf
   cndexptr = C_LOC(cndex)
 EndIf

#if NF_INT_IS_C_INT
 cstatus = nc_put_var1_int(cncid, cvarid, cndexptr, cival)
#elif NF_INT_IS_C_LONG
 cstatus = nc_put_var1_long(cncid, cvarid, cndexptr, cival)
#endif

 status = cstatus

 End Function nf_put_var1_int
!--------------------------------- nf_put_var1_real ------------------------
 Function nf_put_var1_real(ncid, varid, ndex, rval) RESULT(status)

! Write out a 32 bit real variable to location vector ndex in dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,      Intent(IN) :: ncid, varid
 Integer,      Intent(IN) :: ndex(*)
 Real(NFREAL), Intent(IN) :: rval

 Integer                  :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Integer                        :: ndims

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cndex  = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! reverse array order and subtract 1 to get C index 
     cndex(1:ndims) = ndex(ndims:1:-1)-1
   EndIf
   cndexptr = C_LOC(cndex)
 EndIf

#if NF_REAL_IS_C_DOUBLE
 cstatus = nc_put_var1_double(cncid, cvarid, cndexptr, rval)
#else
 cstatus = nc_put_var1_float(cncid, cvarid, cndexptr, rval)
#endif

 status = cstatus

 End Function nf_put_var1_real
!--------------------------------- nf_put_var1_double ----------------------
 Function nf_put_var1_double(ncid, varid, ndex, dval) RESULT(status)

! Write out a 64 bit real variable to location vector ndex in dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,   Intent(IN) :: ncid, varid
 Integer,   Intent(IN) :: ndex(*)
 Real(RK8), Intent(IN) :: dval

 Integer               :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Integer                        :: ndims

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cndex  = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims) ! get varid dimension

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then
     cndex(1:ndims) = ndex(ndims:1:-1)-1
   EndIf
   cndexptr = C_LOC(cndex)
 EndIf

 cstatus = nc_put_var1_double(cncid, cvarid, cndexptr, dval)

 status = cstatus

 End Function nf_put_var1_double
!--------------------------------- nf_put_var1 ------------------------
 Function nf_put_var1(ncid, varid, ndex, values) RESULT(status)

! Write out values of any type. We use a C interop character string to
! hold values. Therefore, an explicit interface to nf_put_var1 should
! not be defined in the calling program to avoid rigid TKR conflict
! Just declare it external

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,                Intent(IN)         :: ncid, varid
 Integer,                Intent(IN)         :: ndex(*)
 Character(KIND=C_CHAR), Intent(IN), TARGET :: values(*)
! Type(C_PTR),            VALUE              :: values
 Integer                                    :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Type(C_PTR)                    :: cvaluesptr ! comment for type(C_PTR) values
 Integer                        :: ndims

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cndex  = 0
 
 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! reverse array order and subtract 1 to get C index 
     cndex(1:ndims) = ndex(ndims:1:-1)-1
   EndIf
   cndexptr = C_LOC(cndex)
 EndIf

 cvaluesptr = C_LOC(values)

 cstatus = nc_put_var1(cncid, cvarid, cndexptr, cvaluesptr)
! cstatus = nc_put_var1(cncid, cvarid, cndexptr, values)

 status = cstatus

 End Function nf_put_var1
!--------------------------------- nf_get_var1_text ------------------------
 Function nf_get_var1_text(ncid, varid, ndex, chval) RESULT(status)

! Read in a single character variable from location vector ndex in dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: ndex(*)
 Character(LEN=1), Intent(OUT) :: chval

 Integer                       :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Integer                        :: ndims

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cndex  =  0
 chval  = ' '

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! reverse array order and subtract 1 to get C index 
     cndex(1:ndims) = ndex(ndims:1:-1) -1
   EndIf
   cndexptr = C_LOC(cndex)
  EndIf
 
 cstatus = nc_get_var1_text(cncid, cvarid, cndexptr, chval)

 status = cstatus

 End Function nf_get_var1_text
!--------------------------------- nf_get_var1_int1 ------------------------
 Function nf_get_var1_int1(ncid, varid, ndex, ival) RESULT(status)

! Read in a 8 bit integer variable from location vector ndex in dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN)  :: ncid, varid
 Integer,              Intent(IN)  :: ndex(*)
 Integer(KIND=NFINT1), Intent(OUT) :: ival

 Integer                           :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Integer(KIND=CINT1)            :: cival
 Integer                        :: ndims

 If (C_SIGNED_CHAR < 0) Then ! schar not supported by processor exit
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cndex  = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! reverse array order and subtract 1 to get C index 
     cndex(1:ndims) = ndex(ndims:1:-1)-1
   EndIf
   cndexptr = C_LOC(cndex)
 EndIf

#if NF_INT1_IS_C_SIGNED_CHAR
 cstatus = nc_get_var1_schar(cncid, cvarid, cndexptr, cival)
#elif NF_INT1_IS_C_SHORT
 cstatus = nc_get_var1_short(cncid, cvarid, cndexptr, cival)
#elif NF_INT1_IS_C_INT
 cstatus = nc_get_var1_int(cncid, cvarid, cndexptr, cival)
#elif NF_INT1_IS_C_LONG
 cstatus = nc_get_var1_long(cncid, cvarid, cndexptr, cival)
#endif
 
 ival   = cival
 status = cstatus

 End Function nf_get_var1_int1
!--------------------------------- nf_get_var1_int2 ------------------------
 Function nf_get_var1_int2(ncid, varid, ndex, ival) RESULT(status)

! Read in a 16 bit integer variable from location vector ndex in dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN)  :: ncid, varid
 Integer,              Intent(IN)  :: ndex(*)
 Integer(KIND=NFINT2), Intent(OUT) :: ival

 Integer                           :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Integer(KIND=CINT2)            :: cival
 Integer                        :: ndims

 If (C_SHORT < 0) Then ! short not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cndex  = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! reverse array order and subtract 1 to get C index 
     cndex(1:ndims) = ndex(ndims:1:-1)-1
   EndIf
   cndexptr = C_LOC(cndex)
 EndIf

#if NF_INT2_IS_C_SHORT
 cstatus = nc_get_var1_short(cncid, cvarid, cndexptr, cival)
#elif NF_INT2_IS_C_INT
 cstatus = nc_get_var1_int(cncid, cvarid, cndexptr, cival)
#elif NF_INT2_IS_C_LONG
 cstatus = nc_get_var1_long(cncid, cvarid, cndexptr, cival)
#endif
 
 ival   = cival
 status = cstatus

 End Function nf_get_var1_int2
!--------------------------------- nf_get_var1_int -------------------------
 Function nf_get_var1_int(ncid, varid, ndex, ival) RESULT(status)

! Read in 32 bit integer variable from location vector ndex in dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid, varid
 Integer, Intent(IN)  :: ndex(*)
 Integer, Intent(OUT) :: ival

 Integer              :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Integer(KIND=CINT)             :: cival
 Integer                        :: ndims

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cndex  = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! reverse array order and subtract 1 to get C index 
     cndex(1:ndims) = ndex(ndims:1:-1)-1
   EndIf
   cndexptr = C_LOC(cndex)
 EndIf

#if NF_INT_IS_C_INT
 cstatus = nc_get_var1_int(cncid, cvarid, cndexptr, cival)
#elif NF_INT_IS_C_LONG
 cstatus = nc_get_var1_long(cncid, cvarid, cndexptr, cival)
#endif

 ival   = cival
 status = cstatus

 End Function nf_get_var1_int
!--------------------------------- nf_get_var1_real ------------------------
 Function nf_get_var1_real(ncid, varid, ndex, rval) RESULT(status)

! Read in a 32 bit real variable to location vector ndex in dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,      Intent(IN)  :: ncid, varid
 Integer,      Intent(IN)  :: ndex(*)
 Real(NFREAL), Intent(OUT) :: rval

 Integer                   :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Integer                        :: ndims

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cndex  = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! reverse array order and subtract 1 to get C index 
     cndex(1:ndims) = ndex(ndims:1:-1)-1
   EndIf
   cndexptr = C_LOC(cndex)
 EndIf

#if NF_REAL_IS_C_DOUBLE
 cstatus = nc_get_var1_double(cncid, cvarid, cndexptr, rval)
#else
 cstatus = nc_get_var1_float(cncid, cvarid, cndexptr, rval)
#endif

 status = cstatus

 End Function nf_get_var1_real
!--------------------------------- nf_get_var1_double ----------------------
 Function nf_get_var1_double(ncid, varid, ndex, dval) RESULT(status)

! Read in a 64 bit real variable to location vector ndex in dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,   Intent(IN)  :: ncid, varid
 Integer,   Intent(IN)  :: ndex(*)
 Real(RK8), Intent(OUT) :: dval

 Integer                :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Integer                        :: ndims

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cndex  = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims) ! get varid dimension

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! reverse array order and subtract 1 to get C index 
     cndex(1:ndims) = ndex(ndims:1:-1)-1
   EndIf
   cndexptr = C_LOC(cndex)
 EndIf

 cstatus = nc_get_var1_double(cncid, cvarid, cndexptr, dval)

 status = cstatus

 End Function nf_get_var1_double
!--------------------------------- nf_get_var1 -------------------------------
 Function nf_get_var1(ncid, varid, ndex, values) RESULT(status)

! Read in values of any type. We use a C interop character string to
! hold values. Therefore, an explicit interface to nf_get_var1 should
! not be defined in the calling program to avoid rigid TKR conflict
! Just declare it external

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,                Intent(IN)          :: ncid, varid
 Integer,                Intent(IN)          :: ndex(*)
 Character(KIND=C_CHAR), Intent(OUT), TARGET :: values(*)
! Type(C_PTR),            VALUE               :: values

 Integer                                     :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cndex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cndexptr
 Type(C_PTR)                    :: cvaluesptr
 Integer                        :: ndims

 cncid  = ncid
 cvarid = varid - 1 ! Subtract one to get C varid
 cndex  = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims) ! get varid dimension

 cndexptr = C_NULL_PTR
 ndims    = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! reverse array order and subtract 1 to get C index 
     cndex(1:ndims) = ndex(ndims:1:-1)-1
   EndIf
   cndexptr = C_LOC(cndex)
 EndIf

 cstatus = nc_get_var1(cncid, cvarid, cndexptr, values)

 status = cstatus

 End Function nf_get_var1
