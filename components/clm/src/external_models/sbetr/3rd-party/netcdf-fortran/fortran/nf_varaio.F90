#include "nfconfig.inc"
!--- Array put/get routines for different types for given start and count ----

! Replacement for fort-varaio.c

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
!                          Updated to pass start and counts as C_PTR
!                          variables
! Version 3.: April 2009 - Updated for netCDF 4.0.1
! Version 4.: April 2010 - Updated for netCDF 4.1.1
!                          Added preprocessor test for int and real types
         
!--------------------------------- nf_put_vara_text ----------------------
 Function nf_put_vara_text(ncid, varid, start, counts, text) RESULT(status)

! Write out a character string to dataset for given start and count vectors

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid
 Integer,          Intent(IN) :: start(*), counts(*)
 Character(LEN=*), Intent(IN) :: text

 Integer                      :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

 cstatus = nc_put_vara_text(cncid, cvarid, cstartptr, ccountsptr, text)

 status = cstatus

 End Function nf_put_vara_text
!--------------------------------- nf_put_vara_text_a ----------------------
 Function nf_put_vara_text_a(ncid, varid, start, counts, text) RESULT(status)

! Write out an array of characters to dataset for given start and count vectors

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid
 Integer,          Intent(IN) :: start(*), counts(*)
 Character(LEN=1), Intent(IN) :: text(*)

 Integer                      :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

 cstatus = nc_put_vara_text(cncid, cvarid, cstartptr, ccountsptr, text)

 status = cstatus

 End Function nf_put_vara_text_a
!--------------------------------- nf_put_vara_int1 ------------------------
 Function nf_put_vara_int1(ncid, varid, start, counts, i1vals) RESULT(status)

! Write out 8 bit integer array to dataset for given start and count vectors

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN) :: ncid, varid
 Integer,              Intent(IN) :: start(*), counts(*)
 Integer(KIND=NFINT1), Intent(IN) :: i1vals(*)

 Integer                          :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 If (C_SIGNED_CHAR < 0) Then ! schar not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

#if NF_INT1_IS_C_SIGNED_CHAR
 cstatus = nc_put_vara_schar(cncid, cvarid, cstartptr, ccountsptr, i1vals)
#elif NF_INT1_IS_C_SHORT
 cstatus = nc_put_vara_short(cncid, cvarid, cstartptr, ccountsptr, i1vals)
#elif NF_INT1_IS_C_INT
 cstatus = nc_put_vara_int(cncid, cvarid, cstartptr, ccountsptr, i1vals)
#elif NF_INT1_IS_C_LONG
 cstatus = nc_put_vara_long(cncid, cvarid, cstartptr, ccountsptr, i1vals)
#endif
 
 status = cstatus

 End Function nf_put_vara_int1
!--------------------------------- nf_put_vara_int2 ------------------------
 Function nf_put_vara_int2(ncid, varid, start, counts, i2vals) RESULT(status)

! Write out 16 bit integer array to dataset for given start and count vectors

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN) :: ncid, varid
 Integer,              Intent(IN) :: start(*), counts(*)
 Integer(KIND=NFINT2), Intent(IN) :: i2vals(*)

 Integer                          :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 If (C_SHORT < 0) Then ! short not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

#if NF_INT2_IS_C_SHORT
 cstatus = nc_put_vara_short(cncid, cvarid, cstartptr, ccountsptr, i2vals)
#elif NF_INT2_IS_C_INT
 cstatus = nc_put_vara_int(cncid, cvarid, cstartptr, ccountsptr, i2vals)
#elif NF_INT2_IS_C_LONG
 cstatus = nc_put_vara_long(cncid, cvarid, cstartptr, ccountsptr, i2vals)
#endif

 status = cstatus

 End Function nf_put_vara_int2
!--------------------------------- nf_put_vara_int -------------------------
 Function nf_put_vara_int(ncid, varid, start, counts, ivals) RESULT(status)

! Write out default integer array to dataset for given start and count vectors

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,        Intent(IN) :: ncid, varid
 Integer,        Intent(IN) :: start(*), counts(*)
 Integer(NFINT), Intent(IN) :: ivals(*)

 Integer                    :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

#if NF_INT_IS_C_INT
 cstatus = nc_put_vara_int(cncid, cvarid, cstartptr, ccountsptr, ivals)
#elif NF_INT_IS_C_LONG
 cstatus = nc_put_vara_long(cncid, cvarid, cstartptr, ccountsptr, ivals)
#endif

 status = cstatus

 End Function nf_put_vara_int
!--------------------------------- nf_put_vara_real ------------------------
 Function nf_put_vara_real(ncid, varid, start, counts, rvals) RESULT(status)

! Write out real(RK4) array to dataset for given start and count vectors

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,      Intent(IN) :: ncid, varid
 Integer,      Intent(IN) :: start(*), counts(*)
 Real(NFREAL), Intent(IN) :: rvals(*)

 Integer                  :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

#if NF_REAL_IS_C_DOUBLE
 cstatus = nc_put_vara_double(cncid, cvarid, cstartptr, ccountsptr, rvals)
#else
 cstatus = nc_put_vara_float(cncid, cvarid, cstartptr, ccountsptr, rvals)
#endif

 status = cstatus

 End Function nf_put_vara_real
!--------------------------------- nf_put_vara_double ----------------------
 Function nf_put_vara_double(ncid, varid, start, counts, dvals) &
                                RESULT(status)

! Write out real(RK8) variable to dataset for given start and count vectors

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,   Intent(IN) :: ncid, varid
 Integer,   Intent(IN) :: start(*), counts(*)
 Real(RK8), Intent(IN) :: dvals(*)

 Integer               :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims) ! get varid dimension

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

 cstatus = nc_put_vara_double(cncid, cvarid, cstartptr, ccountsptr, dvals)

 status = cstatus

 End Function nf_put_vara_double
!--------------------------------- nf_put_vara ------------------------------
 Function nf_put_vara(ncid, varid, start, counts, values) RESULT(status)

! Write out an array of any type. We use a C interop character string to
! pass values. Therefore, an explicit interface to nf_put_vara should not
! be used in the calling routine. Just use external. 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,                Intent(IN)         :: ncid, varid
 Integer,                Intent(IN)         :: start(*), counts(*)
 Character(KIND=C_CHAR), Intent(IN), TARGET :: values(*)
! Type(C_PTR),            VALUE              :: values
 Integer                                    :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr, cvaluesptr
! Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

 cvaluesptr = C_LOC(values)

 cstatus = nc_put_vara(cncid, cvarid, cstartptr, ccountsptr, cvaluesptr)
! cstatus = nc_put_vara(cncid, cvarid, cstartptr, ccountsptr, values)

 status = cstatus

 End Function nf_put_vara
!--------------------------------- nf_get_vara_text ----------------------
 Function nf_get_vara_text(ncid, varid, start, counts, text) RESULT(status)

! Read in a character string from dataset for given start and count vectors 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: start(*), counts(*)
 Character(LEN=*), Intent(OUT) :: text

 Integer                       :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0
 text    = REPEAT(" ", LEN(text))

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

 cstatus = nc_get_vara_text(cncid, cvarid, cstartptr, ccountsptr, text)

 status = cstatus

 End Function nf_get_vara_text
!--------------------------------- nf_get_vara_text_a ----------------------
 Function nf_get_vara_text_a(ncid, varid, start, counts, text) RESULT(status)

! Read in an array of characters for given start and count vectors 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: start(*), counts(*)
 Character(LEN=1), Intent(OUT) :: text(*)

 Integer                       :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

 cstatus = nc_get_vara_text(cncid, cvarid, cstartptr, ccountsptr, text)

 status = cstatus

 End Function nf_get_vara_text_a
!--------------------------------- nf_get_vara_int1 ------------------------
 Function nf_get_vara_int1(ncid, varid, start, counts, i1vals) RESULT(status)

! Read in 8 bit integer array from dataset for given start and count vectors 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN)  :: ncid, varid
 Integer,              Intent(IN)  :: start(*), counts(*)
 Integer(KIND=NFINT1), Intent(OUT) :: i1vals(*)

 Integer                           :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 If (C_SIGNED_CHAR < 0) Then ! schar not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

#if NF_INT1_IS_C_SIGNED_CHAR
 cstatus = nc_get_vara_schar(cncid, cvarid, cstartptr, ccountsptr, i1vals)
#elif NF_INT1_IS_C_SHORT
 cstatus = nc_get_vara_short(cncid, cvarid, cstartptr, ccountsptr, i1vals)
#elif NF_INT1_IS_C_INT
 cstatus = nc_get_vara_int(cncid, cvarid, cstartptr, ccountsptr, i1vals)
#elif NF_INT1_IS_C_LONG
 cstatus = nc_get_vara_long(cncid, cvarid, cstartptr, ccountsptr, i1vals)
#endif

 status = cstatus

 End Function nf_get_vara_int1
!--------------------------------- nf_get_vara_int2 ------------------------
 Function nf_get_vara_int2(ncid, varid, start, counts, i2vals) RESULT(status)

! Read in 16 bit integer array from dataset for given start and count vectors

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN)  :: ncid, varid
 Integer,              Intent(IN)  :: start(*), counts(*)
 Integer(KIND=NFINT2), Intent(OUT) :: i2vals(*)

 Integer                           :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 If (C_SHORT < 0) Then ! short not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

#if NF_INT2_IS_C_SHORT
 cstatus = nc_get_vara_short(cncid, cvarid, cstartptr, ccountsptr, i2vals)
#elif NF_INT2_IS_C_INT
 cstatus = nc_get_vara_int(cncid, cvarid, cstartptr, ccountsptr, i2vals)
#elif NF_INT2_IS_C_LONG
 cstatus = nc_get_vara_long(cncid, cvarid, cstartptr, ccountsptr, i2vals)
#endif

 status = cstatus

 End Function nf_get_vara_int2
!--------------------------------- nf_get_vara_int -------------------------
 Function nf_get_vara_int(ncid, varid, start, counts, ivals) RESULT(status)

! Read in default integer array from dataset for given start and count vectors 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,        Intent(IN)  :: ncid, varid
 Integer,        Intent(IN)  :: start(*), counts(*)
 Integer(NFINT), Intent(OUT) :: ivals(*)

 Integer                     :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

#if NF_INT_IS_C_INT
 cstatus = nc_get_vara_int(cncid, cvarid, cstartptr, ccountsptr, ivals)
#elif NF_INT_IS_C_LONG 
 cstatus = nc_get_vara_long(cncid, cvarid, cstartptr, ccountsptr, ivals)
#endif

 status = cstatus

 End Function nf_get_vara_int
!--------------------------------- nf_get_vara_real ------------------------
 Function nf_get_vara_real(ncid, varid, start, counts, rvals) RESULT(status)

! Read in real(RK4) array from dataset for given start and count vectors

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,      Intent(IN)  :: ncid, varid
 Integer,      Intent(IN)  :: start(*), counts(*)
 Real(NFREAL), Intent(OUT) :: rvals(*)

 Integer                   :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

#if NF_REAL_IS_C_DOUBLE
 cstatus = nc_get_vara_double(cncid, cvarid, cstartptr, ccountsptr, rvals)
#else
 cstatus = nc_get_vara_float(cncid, cvarid, cstartptr, ccountsptr, rvals)
#endif

 status = cstatus

 End Function nf_get_vara_real
!--------------------------------- nf_get_vara_double ----------------------
 Function nf_get_vara_double(ncid, varid, start, counts, dvals) &
                             RESULT(status)

! Read in Real(RK8) array from dataset for given start and count vectors

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,   Intent(IN)  :: ncid, varid
 Integer,   Intent(IN)  :: start(*), counts(*)
 Real(RK8), Intent(OUT) :: dvals(*)

 Integer                :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims) ! get varid dimension

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

 cstatus = nc_get_vara_double(cncid, cvarid, cstartptr, ccountsptr, dvals)

 status = cstatus

 End Function nf_get_vara_double
!--------------------------------- nf_get_vara ------------------------------
 Function nf_get_vara(ncid, varid, start, counts, values) RESULT(status)

! Read in an array of any type. We use a C interop character string to
! pass values. Therefore, an explicit interface to nf_put_vara should not
! be used in the calling routine. Just use external. 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,                Intent(IN)            :: ncid, varid
 Integer,                Intent(IN)            :: start(*), counts(*)
 Character(KIND=C_CHAR), Intent(INOUT), TARGET :: values(*)
! Type(C_PTR),            VALUE                 :: values
 Integer                                       :: status

 Integer(KIND=C_INT)            :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1 ! Subtract 1 to get C varid
 cstart  = 0
 ccounts = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! flip array order for C and subtract 1 from start
     cstart(1:ndims)  = start(ndims:1:-1)-1
     ccounts(1:ndims) = counts(ndims:1:-1)
   EndIf
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 EndIf

 cstatus = nc_get_vara(cncid, cvarid, cstartptr, ccountsptr, values)

 status = cstatus

 End Function nf_get_vara
