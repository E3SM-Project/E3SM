#include "nfconfig.inc"
!------------ Array/string put/get routines for a given varid ----------------

! Replacement for fort-vario.c

! Written by: Richard Weed, Ph.D.
!             Center For Advanced Vehicular Systems 
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
!                          Added preprocessor tests for int and real types

!--------------------------------- nf_put_var_text -----------------------
 Function nf_put_var_text(ncid, varid, text) RESULT(status)

! Write out a character string to dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid
 Character(LEN=*), Intent(IN) :: text

 Integer                      :: status

 Integer(KIND=C_INT) :: cncid, cvarid,  cstatus

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

 cstatus = nc_put_var_text(cncid, cvarid, text)

 status = cstatus

 End Function nf_put_var_text
!--------------------------------- nf_put_var_text_a -----------------------
 Function nf_put_var_text_a(ncid, varid, text) RESULT(status)

! Write out array of characters to dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN) :: ncid, varid
 Character(LEN=1), Intent(IN) :: text(*)

 Integer                      :: status

 Integer(KIND=C_INT) :: cncid, cvarid,  cstatus

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

 cstatus = nc_put_var_text(cncid, cvarid, text)

 status = cstatus

 End Function nf_put_var_text_a
!--------------------------------- nf_put_var_int1 -------------------------
 Function nf_put_var_int1(ncid, varid, i1vals) RESULT(status)

! Write out 8 bit integer array to dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN) :: ncid, varid
 Integer(KIND=NFINT1), Intent(IN) :: i1vals(*)

 Integer                          :: status

 Integer(KIND=C_INT) :: cncid, cvarid,  cstatus

 If (C_SIGNED_CHAR < 0) Then ! schar not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

#if NF_INT1_IS_C_SIGNED_CHAR
 cstatus = nc_put_var_schar(cncid, cvarid, i1vals)
#elif NF_INT1_IS_C_SHORT
 cstatus = nc_put_var_short(cncid, cvarid, i1vals)
#elif NF_INT1_IS_C_INT
 cstatus = nc_put_var_int(cncid, cvarid, i1vals)
#elif NF_INT1_IS_C_LONG
 cstatus = nc_put_var_long(cncid, cvarid, i1vals)
#endif

 status = cstatus

 End Function nf_put_var_int1
!--------------------------------- nf_put_var_int2 -------------------------
 Function nf_put_var_int2(ncid, varid, i2vals) RESULT(status)

! Write out 16 bit integer array to dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN) :: ncid, varid
 Integer(KIND=NFINT2), Intent(IN) :: i2vals(*)

 Integer                          :: status

 Integer(KIND=C_INT) :: cncid, cvarid,  cstatus

 If (C_SHORT < 0) Then ! short not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

#if NF_INT2_IS_C_SHORT
 cstatus = nc_put_var_short(cncid, cvarid, i2vals)
#elif NF_INT2_IS_C_INT
 cstatus = nc_put_var_int(cncid, cvarid, i2vals)
#elif NF_INT2_IS_C_LONG
 cstatus = nc_put_var_long(cncid, cvarid, i2vals)
#endif

 status = cstatus

 End Function nf_put_var_int2
!--------------------------------- nf_put_var_int --------------------------
 Function nf_put_var_int(ncid, varid, ivals) RESULT(status)

! Write out 32 bit integer array to dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,        Intent(IN) :: ncid, varid
 Integer(NFINT), Intent(IN) :: ivals(*)

 Integer                    :: status

 Integer(KIND=C_INT) :: cncid, cvarid,  cstatus

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

#if NF_INT_IS_C_INT
 cstatus = nc_put_var_int(cncid, cvarid, ivals)
#elif NF_INT_IS_C_LONG
 cstatus = nc_put_var_long(cncid, cvarid, ivals)
#endif

 status = cstatus

 End Function nf_put_var_int
!--------------------------------- nf_put_var_real -------------------------
 Function nf_put_var_real(ncid, varid, rvals) RESULT(status)

! Write out 32 bit real array to dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,         Intent(IN) :: ncid, varid
 Real(NFREAL),    Intent(IN) :: rvals(*)
 Integer                     :: status

 Integer(KIND=C_INT) :: cncid, cvarid,  cstatus

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

#if NF_REAL_IS_C_DOUBLE
 cstatus = nc_put_var_double(cncid, cvarid, rvals)
#else
 cstatus = nc_put_var_float(cncid, cvarid, rvals)
#endif

 status = cstatus

 End Function nf_put_var_real
!--------------------------------- nf_put_var_double -----------------------
 Function nf_put_var_double(ncid, varid, dvals) RESULT(status)

! Write out 64 bit real array to dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,   Intent(IN) :: ncid, varid
 Real(RK8), Intent(IN) :: dvals(*)

 Integer               :: status

 Integer(KIND=C_INT) :: cncid, cvarid,  cstatus

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

 cstatus = nc_put_var_double(cncid, cvarid, dvals)

 status = cstatus

 End Function nf_put_var_double
!--------------------------------- nf_get_var_text -----------------------
 Function nf_get_var_text(ncid, varid, text) RESULT(status)

! Read in a character string from dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=*), Intent(OUT) :: text

 Integer                       :: status

 Integer(KIND=C_INT) :: cncid, cvarid,  cstatus

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid
 text   = REPEAT(" ", LEN(text))

 cstatus = nc_get_var_text(cncid, cvarid, text)

 status = cstatus

 End Function nf_get_var_text
!--------------------------------- nf_get_var_text_a -----------------------
 Function nf_get_var_text_a(ncid, varid, text) RESULT(status)

! Read in array of characters from dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Character(LEN=1), Intent(OUT) :: text(*)

 Integer                       :: status

 Integer(KIND=C_INT) :: cncid, cvarid,  cstatus

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

 cstatus = nc_get_var_text(cncid, cvarid, text)

 status = cstatus

 End Function nf_get_var_text_a
!--------------------------------- nf_get_var_int1 -------------------------
 Function nf_get_var_int1(ncid, varid, i1vals) RESULT(status)

! Read in 8 bit integer array from dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN)  :: ncid, varid
 Integer(KIND=NFINT1), Intent(OUT) :: i1vals(*)

 Integer                           :: status

 Integer(KIND=C_INT) :: cncid, cvarid,  cstatus

 If (C_SIGNED_CHAR < 0) Then ! schar not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

#if NF_INT1_IS_C_SIGNED_CHAR
 cstatus = nc_get_var_schar(cncid, cvarid, i1vals)
#elif NF_INT1_IS_C_SHORT
 cstatus = nc_get_var_short(cncid, cvarid, i1vals)
#elif NF_INT1_IS_C_INT
 cstatus = nc_get_var_int(cncid, cvarid, i1vals)
#elif NF_INT1_IS_C_LONG
 cstatus = nc_get_var_long(cncid, cvarid, i1vals)
#endif

 status = cstatus

 End Function nf_get_var_int1
!--------------------------------- nf_get_var_int2 -------------------------
 Function nf_get_var_int2(ncid, varid, i2vals) RESULT(status)

! Read in 16 bit integer array from dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,              Intent(IN)  :: ncid, varid
 Integer(KIND=NFINT2), Intent(OUT) :: i2vals(*)

 Integer                           :: status

 Integer(KIND=C_INT) :: cncid, cvarid,  cstatus

 If (C_SHORT < 0) Then ! short not supported by processor
   status = NC_EBADTYPE
   RETURN
 EndIf

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

#if NF_INT2_IS_C_SHORT
 cstatus = nc_get_var_short(cncid, cvarid, i2vals)
#elif NF_INT2_IS_C_INT
 cstatus = nc_get_var_int(cncid, cvarid, i2vals)
#elif NF_INT2_IS_C_LONG
 cstatus = nc_get_var_long(cncid, cvarid, i2vals)
#endif

 status = cstatus

 End Function nf_get_var_int2
!--------------------------------- nf_get_var_int --------------------------
 Function nf_get_var_int(ncid, varid, ivals) RESULT(status)

! Read in default integer array from dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,        Intent(IN)  :: ncid, varid
 Integer(NFINT), Intent(OUT) :: ivals(*)

 Integer                     :: status

 Integer(KIND=C_INT) :: cncid, cvarid,  cstatus

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

#if NF_INT_IS_C_INT
 cstatus = nc_get_var_int(cncid, cvarid, ivals)
#elif NF_INT_IS_C_LONG
 cstatus = nc_get_var_long(cncid, cvarid, ivals)
#endif

 status = cstatus

 End Function nf_get_var_int
!--------------------------------- nf_get_var_real -------------------------
 Function nf_get_var_real(ncid, varid, rvals) RESULT(status)

! Read in 32 bit real array from dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,      Intent(IN)  :: ncid, varid
 Real(NFREAL), Intent(OUT) :: rvals(*)

 Integer                   :: status

 Integer(KIND=C_INT) :: cncid, cvarid,  cstatus

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

#if NF_REAL_IS_C_DOUBLE
 cstatus = nc_get_var_double(cncid, cvarid, rvals)
#else
 cstatus = nc_get_var_float(cncid, cvarid, rvals)
#endif

 status = cstatus

 End Function nf_get_var_real
!--------------------------------- nf_get_var_double -----------------------
 Function nf_get_var_double(ncid, varid, dvals) RESULT(status)

! Read in 64 bit real array from dataset

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,   Intent(IN)  :: ncid, varid
 Real(RK8), Intent(OUT) :: dvals(*)

 Integer                :: status

 Integer(KIND=C_INT) :: cncid, cvarid, cstatus

 cncid  = ncid
 cvarid = varid - 1 ! Subtract 1 to get C varid

 cstatus = nc_get_var_double(cncid, cvarid, dvals)

 status = cstatus

 End Function nf_get_var_double
