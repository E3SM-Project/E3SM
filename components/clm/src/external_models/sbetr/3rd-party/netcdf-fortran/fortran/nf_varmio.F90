#include "nfconfig.inc"
!- Array/string put/get routines given start, count, stride, and map vectors - 

! Replacement for fort-varmio.c

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
!                          Updated to pass start, counts, strides, and
!                          maps as C_PTR variables
! Version 3.: April 2009 - Updated for netCDF 4.0.1
! Version 4.: April 2010 - Updated for netCDF 4.1.1
!                          Added preprocessor tests for int and real types

!--------------------------------- nf_put_varm_text ----------------------
 Function nf_put_varm_text(ncid, varid, start, counts, strides, maps, &
                           text) RESULT(status)

! Write out a character string to dataset given start, count, stride and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid
 Integer,          Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Character(LEN=*), Intent(IN) :: text 

 Integer                      :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cmapsptr
 Integer                           :: ndims

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0 

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf

 cstatus = nc_put_varm_text(cncid, cvarid, cstartptr, ccountsptr, &
                            cstridesptr, cmapsptr, text)

 status = cstatus

 End Function nf_put_varm_text
!--------------------------------- nf_put_varm_text_a ----------------------
 Function nf_put_varm_text_a(ncid, varid, start, counts, strides, maps, &
                                text) RESULT(status)

! Write out array of characters to dataset given start, count, stride and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid
 Integer,          Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Character(LEN=1), Intent(IN) :: text(*)

 Integer                      :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr,&
                                      cmapsptr
 Integer                           :: ndims

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf

 cstatus = nc_put_varm_text(cncid, cvarid, cstartptr, ccountsptr, &
                            cstridesptr, cmapsptr, text) 

 status = cstatus

 End Function nf_put_varm_text_a
!--------------------------------- nf_put_varm_int1 ------------------------
 Function nf_put_varm_int1(ncid, varid, start, counts, strides, maps, &
                           i1vals) RESULT(status)

! Write out 8 bit integer array given start, count, stride and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN) :: ncid, varid
 Integer,              Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Integer(KIND=NFINT1), Intent(IN) :: i1vals(*)

 Integer                          :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr,&
                                      cmapsptr
 Integer                           :: ndims

 If (C_SIGNED_CHAR < 0) Then ! schar not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf

#if NF_INT1_IS_C_SIGNED_CHAR
 cstatus = nc_put_varm_schar(cncid, cvarid, cstartptr, ccountsptr, &
                             cstridesptr, cmapsptr, i1vals) 
#elif NF_INT1_IS_C_SHORT
 cstatus = nc_put_varm_short(cncid, cvarid, cstartptr, ccountsptr, &
                             cstridesptr, cmapsptr, i1vals) 
#elif NF_INT1_IS_C_INT
 cstatus = nc_put_varm_int(cncid, cvarid, cstartptr, ccountsptr, &
                           cstridesptr, cmapsptr, i1vals) 
#elif NF_INT1_IS_C_LONG
 cstatus = nc_put_varm_long(cncid, cvarid, cstartptr, ccountsptr, &
                            cstridesptr, cmapsptr, i1vals) 
#endif

 status = cstatus

 End Function nf_put_varm_int1
!--------------------------------- nf_put_varm_int2 ------------------------
 Function nf_put_varm_int2(ncid, varid, start, counts, strides, maps, &
                           i2vals) RESULT(status)

! Write out 16 bit integer array given start, count, stride and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN) :: ncid, varid
 Integer,              Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Integer(KIND=NFINT2), Intent(IN) :: i2vals(*)

 Integer                          :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cmapsptr
 Integer                           :: ndims

 If (C_SHORT < 0) Then ! short not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf

#if NF_INT2_IS_C_SHORT
 cstatus = nc_put_varm_short(cncid, cvarid, cstartptr, ccountsptr, &
                             cstridesptr, cmapsptr, i2vals) 
#elif NF_INT2_IS_C_INT
 cstatus = nc_put_varm_int(cncid, cvarid, cstartptr, ccountsptr, &
                           cstridesptr, cmapsptr, i2vals) 
#elif NF_INT2_IS_C_LONG
 cstatus = nc_put_varm_long(cncid, cvarid, cstartptr, ccountsptr, &
                            cstridesptr, cmapsptr, i2vals)
#endif 

 status = cstatus

 End Function nf_put_varm_int2
!--------------------------------- nf_put_varm_int -------------------------
 Function nf_put_varm_int(ncid, varid, start, counts, strides, maps, &
                             ivals) RESULT(status)

! Write out default integer array given start, count, stride and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,        Intent(IN) :: ncid, varid
 Integer,        Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Integer(NFINT), Intent(IN) :: ivals(*)

 Integer                    :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cmapsptr
 Integer                           :: ndims

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf

#if NF_INT_IS_C_INT
 cstatus = nc_put_varm_int(cncid, cvarid, cstartptr, ccountsptr, &
                           cstridesptr, cmapsptr, ivals) 
#elif NF_INT_IS_C_LONG
 cstatus = nc_put_varm_long(cncid, cvarid, cstartptr, ccountsptr, &
                            cstridesptr, cmapsptr, ivals) 
#endif

 status = cstatus

 End Function nf_put_varm_int
!--------------------------------- nf_put_varm_real ------------------------
 Function nf_put_varm_real(ncid, varid, start, counts, strides, maps, &
                              rvals) RESULT(status)

! Write out 32 bit real array given start, count, stride and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,      Intent(IN) :: ncid, varid
 Integer,      Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Real(NFREAL), Intent(IN) :: rvals(*)

 Integer                  :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cmapsptr
 Integer                           :: ndims

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf

#if NF_REAL_IS_C_DOUBLE
 cstatus = nc_put_varm_double(cncid, cvarid, cstartptr, ccountsptr, &
                              cstridesptr, cmapsptr, rvals) 
#else
 cstatus = nc_put_varm_float(cncid, cvarid, cstartptr, ccountsptr, &
                             cstridesptr, cmapsptr, rvals) 
#endif
 status = cstatus

 End Function nf_put_varm_real
!--------------------------------- nf_put_varm_double ----------------------
 Function nf_put_varm_double(ncid, varid, start, counts, strides, maps, &
                                dvals) RESULT(status)

! Write out 64 bit real array given start, count, stride, and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,   Intent(IN) :: ncid, varid
 Integer,   Intent(IN) :: start(*), counts(*), strides(*), maps(*)
 Real(RK8), Intent(IN) :: dvals(*)

 Integer               :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cmapsptr
 Integer                           :: ndims

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf
  
 cstatus = nc_put_varm_double(cncid, cvarid, cstartptr, ccountsptr, &
                              cstridesptr, cmapsptr, dvals) 

 status = cstatus

 End Function nf_put_varm_double
!--------------------------------- nf_get_varm_text ----------------------
 Function nf_get_varm_text(ncid, varid, start, counts, strides, maps, &
                           text) RESULT(status)

! Read in a character string to dataset given start, count, stride, and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Character(LEN=*), Intent(OUT) :: text 

 Integer                       :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cmapsptr
 Integer                           :: ndims

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0
 text     = REPEAT(" ",LEN(text))

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf

 cstatus = nc_get_varm_text(cncid, cvarid, cstartptr, ccountsptr, &
                            cstridesptr, cmapsptr, text) 

 status = cstatus

 End Function nf_get_varm_text
!--------------------------------- nf_get_varm_text_a ----------------------
 Function nf_get_varm_text_a(ncid, varid, start, counts, strides, maps, &
                             text) RESULT(status)

! Read in array of characters from dataset given start, count, stride and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Character(LEN=1), Intent(OUT) :: text(*) 

 Integer                       :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cmapsptr
 Integer                           :: ndims

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf

 cstatus = nc_get_varm_text(cncid, cvarid, cstartptr, ccountsptr, &
                            cstridesptr, cmapsptr, text) 

 status = cstatus

 End Function nf_get_varm_text_a
!--------------------------------- nf_get_varm_int1 ------------------------
 Function nf_get_varm_int1(ncid, varid, start, counts, strides, maps, &
                           i1vals) RESULT(status)

! Read in 8 bit integer array given start, count, stride and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN)  :: ncid, varid
 Integer,              Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Integer(KIND=NFINT1), Intent(OUT) :: i1vals(*)

 Integer                           :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cmapsptr
 Integer                           :: ndims

 If (C_SIGNED_CHAR < 0) Then ! schar not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf

#if NF_INT1_IS_C_SIGNED_CHAR
 cstatus = nc_get_varm_schar(cncid, cvarid, cstartptr, ccountsptr, &
                             cstridesptr, cmapsptr, i1vals) 
#elif NF_INT1_IS_C_SHORT
 cstatus = nc_get_varm_short(cncid, cvarid, cstartptr, ccountsptr, &
                             cstridesptr, cmapsptr, i1vals) 
#elif NF_INT1_IS_C_INT
 cstatus = nc_get_varm_int(cncid, cvarid, cstartptr, ccountsptr, &
                           cstridesptr, cmapsptr, i1vals) 
#elif NF_INT1_IS_C_LONG
 cstatus = nc_get_varm_long(cncid, cvarid, cstartptr, ccountsptr, &
                            cstridesptr, cmapsptr, i1vals) 
#endif

 status = cstatus

 End Function nf_get_varm_int1
!--------------------------------- nf_get_varm_int2 ------------------------
 Function nf_get_varm_int2(ncid, varid, start, counts, strides, maps, &
                              i2vals) RESULT(status)

! Read in 16 bit integer array given start, count, stride and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN)  :: ncid, varid
 Integer,              Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Integer(KIND=NFINT2), Intent(OUT) :: i2vals(*)

 Integer                           :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cmapsptr
 Integer                           :: ndims

 If (C_SHORT < 0) Then ! short not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf

#if NF_INT2_IS_C_SHORT
 cstatus = nc_get_varm_short(cncid, cvarid, cstartptr, ccountsptr, &
                             cstridesptr, cmapsptr, i2vals) 
#elif NF_INT2_IS_C_INT
 cstatus = nc_get_varm_int(cncid, cvarid, cstartptr, ccountsptr, &
                           cstridesptr, cmapsptr, i2vals) 
#elif NF_INT2_IS_C_LONG
 cstatus = nc_get_varm_long(cncid, cvarid, cstartptr, ccountsptr, &
                            cstridesptr, cmapsptr, i2vals)
#endif 

 status = cstatus

 End Function nf_get_varm_int2
!--------------------------------- nf_get_varm_int -------------------------
 Function nf_get_varm_int(ncid, varid, start, counts, strides, maps, &
                             ivals) RESULT(status)

! Read in default integer array given start, count, stride and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,        Intent(IN)  :: ncid, varid
 Integer,        Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Integer(NFINT), Intent(OUT) :: ivals(*)

 Integer                     :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cmapsptr
 Integer                              :: ndims

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf

#if NF_INT_IS_C_INT
 cstatus = nc_get_varm_int(cncid, cvarid, cstartptr, ccountsptr, &
                           cstridesptr, cmapsptr, ivals) 
#elif NF_INT_IS_C_LONG
 cstatus = nc_get_varm_long(cncid, cvarid, cstartptr, ccountsptr, &
                            cstridesptr, cmapsptr, ivals)
#endif 

 status = cstatus

 End Function nf_get_varm_int
!--------------------------------- nf_get_varm_real ------------------------
 Function nf_get_varm_real(ncid, varid, start, counts, strides, maps, &
                              rvals) RESULT(status)

! Read in 32 bit real array given start, count, stride and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,      Intent(IN)  :: ncid, varid
 Integer,      Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Real(NFREAL), Intent(OUT) :: rvals(*)

 Integer                   :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cmapsptr
 Integer                           :: ndims

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf

#if NF_REAL_IS_C_DOUBLE
 cstatus = nc_get_varm_double(cncid, cvarid, cstartptr, ccountsptr, &
                              cstridesptr, cmapsptr, rvals) 
#else
 cstatus = nc_get_varm_float(cncid, cvarid, cstartptr, ccountsptr, &
                             cstridesptr, cmapsptr, rvals) 
#endif

 status = cstatus

 End Function nf_get_varm_real
!--------------------------------- nf_get_varm_double ----------------------
 Function nf_get_varm_double(ncid, varid, start, counts, strides, maps, &
                             dvals) RESULT(status)

! Read in 64 bit real array given start, count, stride and map

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,   Intent(IN)  :: ncid, varid
 Integer,   Intent(IN)  :: start(*), counts(*), strides(*), maps(*)
 Real(RK8), Intent(OUT) :: dvals(*)

 Integer                :: status

 Integer(KIND=C_INT)               :: cncid, cvarid, cndims, cstat1, cstatus
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cmaps(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cmapsptr
 Integer                           :: ndims

 cncid    = ncid
 cvarid   = varid -1 ! Subtract 1 to get C varid
 cstart   = 0
 ccounts  = 0
 cstrides = 1
 cmaps    = 0

 cstat1 = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cmapsptr    = C_NULL_PTR
 ndims       = cndims

 If (cstat1 == NC_NOERR) Then
   If (ndims > 0) Then ! Flip arrays to C order and subtract 1 from start
     cstart(1:ndims)   = start(ndims:1:-1)-1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1)
     cmaps(1:ndims)    = maps(ndims:1:-1)
   EndIf
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
   cmapsptr    = C_LOC(cmaps)
 EndIf

 cstatus = nc_get_varm_double(cncid, cvarid, cstartptr, ccountsptr, &
                              cstridesptr, cmapsptr, dvals) 

 status = cstatus

 End Function nf_get_varm_double
