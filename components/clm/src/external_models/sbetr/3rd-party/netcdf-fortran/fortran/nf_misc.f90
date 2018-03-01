!-- Routines for processing error messages, obtaining version numbers, etc. --

! Replacement for fort-misc.c

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
 
!-------------------------------- nf_inq_libvers ---------------------------
 Function nf_inq_libvers() RESULT(vermsg)

! Return string with current version of NetCDF library

 USE netcdf_nc_interfaces

 Implicit NONE

 Character(LEN=80) :: vermsg

 Character(LEN=81), Pointer :: fstrptr
 TYPE(C_PTR)                :: cstrptr
 Integer                    :: inull, ilen

 vermsg = REPEAT(" ", LEN(vermsg)) !initialize vermsg to blanks

! Get C character pointer returned by nc_inq_vers and associate it
! Fortran character pointer (fstrptr). Have to do this when the C
! routine allocates space for the pointer and/or knows where it lives
! not Fortran. This is also how you can pass character data back to
! FORTRAN from C using a C function that returns a character pointer
! instead using a void jacket function and passing the string as a hidden 
! argument. At least this is how cfortran.h appears to do it. 

 
 NULLIFY(fstrptr) ! Nullify fstrptr

 cstrptr = nc_inq_libvers()         ! Get C pointer to version string and

 Call C_F_POINTER(cstrptr, fstrptr) ! associate it with FORTRAN pointer

! Locate first C null character and then set it and remaining characters
! in string to blanks

 ilen  = LEN_TRIM(fstrptr)
 inull = SCAN(fstrptr,C_NULL_CHAR)
 If (inull /= 0) ilen = inull-1
 ilen = MAX(1, MIN(ilen,80)) ! Limit ilen to >=1 and <=80

! Load return value with trimmed fstrptr string

 vermsg(1:ilen) = fstrptr(1:ilen)

 End Function nf_inq_libvers
!-------------------------------- nf_stderror ------------------------------
 Function nf_strerror(ncerr) RESULT(errmsg)

! Returns an error message string given static error code ncerr

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer(KIND=C_INT), Intent(IN) :: ncerr

 Character(LEN=80)               :: errmsg

 Character(LEN=81), Pointer :: fstrptr
 TYPE(C_PTR)                :: cstrptr
 Integer                    :: inull, ilen
 Integer(KIND=C_INT)        :: cncerr

 errmsg = REPEAT(" ", LEN(errmsg)) !initialize errmsg to blanks

! Get C character pointer returned by nc_stderror and associate it
! Fortran character pointer (fstrptr). Have to do this when the C
! routine allocates space for the pointer and/or knows where it lives
! not Fortran. This is also how you can pass character data back to
! FORTRAN from C using a C function that returns a character pointer
! instead using a void jacket function and passing the string as a hidden 
! argument. At least this is how cfortran.h appears to do it. 
 
 NULLIFY(fstrptr) ! Nullify fstrptr

 cncerr  = ncerr

 cstrptr = nc_strerror(cncerr)      ! Return C character pointer and
 Call C_F_POINTER(cstrptr, fstrptr) ! associate C ptr with FORTRAN pointer

! Locate first C null character and then set it and remaining characters
! in string to blanks 

 ilen  = LEN_TRIM(fstrptr)
 inull = SCAN(fstrptr,C_NULL_CHAR)
 If (inull /= 0) ilen = inull-1 
 ilen  = MAX(1, MIN(ilen,80)) ! Limit ilen to >=1 and <=80

! Load return value with trimmed fstrptr string

 errmsg(1:ilen) = fstrptr(1:ilen)

 End Function nf_strerror
!-------------------------------- nf_issyserr ------------------------------
 Function nf_issyserr(nerr) RESULT(status)

! Check to see if nerr is > 0

 Integer, Intent(IN) :: nerr

 Logical             :: status

 status = (nerr > 0) 

 End Function nf_issyserr
