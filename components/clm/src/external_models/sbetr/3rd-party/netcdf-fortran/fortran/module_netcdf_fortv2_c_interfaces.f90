Module netcdf_fortv2_c_interfaces

! Fortran 20003 interfaces to C routines in fort_v2compat.c called by
! the V2 Fortran interfaces. Interface routine names are the same
! as the C routine names.

! Written by : Richard Weed, Ph.D.
!              Center for Advanced Vehicular Systems 
!              Mississipi State University
!              rweed@cavs.msstate.edu


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

! Version 1.: May,   2006 - Initial version 2 interfaces
! Version 2.; April, 2009 - Redone to reflect passing void data types
!                           in C with C_PTR and C_CHAR strings and 
!                           NetCDF 4.0.1
! Version 3.; April, 2010 - Updated to NetCDF 4.1.1

! USE ISO_C_BINDING

! USE NETCDF_NC_DATA,       ONLY: C_PTRDIFF_T
! USE NETCDF_NC_INTERFACES, ONLY: addCNullChar, stripCNullChar
 USE NETCDF_NC_INTERFACES

 Implicit NONE

! The following interfaces are for the netCDF V2 functions. Note that
! the actual C routines return a void pointer for arrays etc. This
! forced me to adopt a commonly used kludge for interfacing old Fortran
! 77 with C, namely, passing the void pointer to an array of C_CHARs.

! Also note that each interface has an explicit USE ISO_C_BINDING. A better
! solution is to use the F2003 IMPORT statement (I originally had it this way)
! However its best to leave the interfaces as is for now because there might
! be a few compilers out there that support most of the C-interop facility but
! for some reason haven't implemented IMPORT yet.

!  Begin fortv2 C interface definitions

!-------------------------------- c_ncpopt ------------------------------------
Interface
 Subroutine c_ncpopt(val) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE :: val

 End Subroutine c_ncpopt
End Interface
!-------------------------------- c_ncgopt ------------------------------------
Interface
 Subroutine c_ncgopt(val) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), Intent(OUT) :: val

 End Subroutine c_ncgopt
End Interface
!-------------------------------- c_nccre -------------------------------------
Interface
 Function c_nccre(pathname, clobmode, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Character(KIND=C_CHAR), Intent(IN)  :: pathname(*)
 Integer(KIND=C_INT),    VALUE       :: clobmode 
 Integer(KIND=C_INT),    Intent(OUT) :: rcode 
 
 Integer(KIND=C_INT)                 :: c_nccre

 End Function c_nccre
End Interface
!-------------------------------- c_ncopn -------------------------------------
Interface
 Function c_ncopn(pathname, rwmode, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Character(KIND=C_CHAR), Intent(IN)  :: pathname(*)
 Integer(KIND=C_INT),    VALUE       :: rwmode 
 Integer(KIND=C_INT),    Intent(OUT) :: rcode
 
 Integer(KIND=C_INT)                 :: c_ncopn

 End Function c_ncopn
End Interface
!-------------------------------- c_ncddef ------------------------------------
Interface
 Function c_ncddef(ncid, dimname, dimlen, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid, dimlen
 Character(KIND=C_CHAR), Intent(IN)  :: dimname(*)
 Integer(KIND=C_INT),    Intent(OUT) :: rcode
 
 Integer(KIND=C_INT)                 :: c_ncddef

 End Function c_ncddef
End Interface
!-------------------------------- c_ncdid -------------------------------------
Interface
 Function c_ncdid(ncid, dimname, rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid
 Character(KIND=C_CHAR), Intent(IN)  :: dimname(*)
 Integer(KIND=C_INT),    Intent(OUT) :: rcode
 
 Integer(KIND=C_INT)                 :: c_ncdid

 End Function c_ncdid
End Interface
!-------------------------------- c_ncvdef ------------------------------------
Interface
 Function c_ncvdef(ncid, varname, datatype, ndims, dimids, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid
 Character(KIND=C_CHAR), Intent(IN)  :: varname(*)
 Integer(KIND=C_INT),    VALUE       :: datatype ! nc_type variable in C
 Integer(KIND=C_INT),    VALUE       :: ndims
 Integer(KIND=C_INT),    Intent(IN)  :: dimids(*)
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 Integer(KIND=C_INT)                 :: c_ncvdef

 End Function c_ncvdef
End Interface
!-------------------------------- c_ncvid -------------------------------------
Interface
 Function c_ncvid(ncid, varname, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid
 Character(KIND=C_CHAR), Intent(IN)  :: varname(*)
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 Integer(KIND=C_INT)                 :: c_ncvid

 End Function c_ncvid
End Interface
!-------------------------------- c_nctlen ------------------------------------
Interface
 Function c_nctlen(datatype, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: datatype ! nc_type var in C
 Integer(KIND=C_INT), Intent(OUT) :: rcode

 Integer(KIND=C_INT)              :: c_nctlen

 End Function c_nctlen
End Interface
!-------------------------------- c_ncclos ------------------------------------
Interface
 Subroutine c_ncclos(ncid, rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid 
 Integer(KIND=C_INT), Intent(OUT) :: rcode

 End Subroutine c_ncclos
End Interface
!-------------------------------- c_ncredf ------------------------------------
Interface
 Subroutine c_ncredf(ncid, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid 
 Integer(KIND=C_INT), Intent(OUT) :: rcode

 End Subroutine c_ncredf
End Interface
!-------------------------------- c_ncendf ------------------------------------
Interface
 Subroutine c_ncendf(ncid, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid 
 Integer(KIND=C_INT), Intent(OUT) :: rcode

 End Subroutine c_ncendf
End Interface
!-------------------------------- c_ncinq -------------------------------------
Interface
 Subroutine c_ncinq(ncid, indims, invars, inatts, irecdim, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid
 Integer(KIND=C_INT), Intent(OUT) :: indims, invars, inatts, irecdim, rcode

 End Subroutine c_ncinq
End Interface
!-------------------------------- c_ncsnc -------------------------------------
Interface
 Subroutine c_ncsnc(ncid, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid 
 Integer(KIND=C_INT), Intent(OUT) :: rcode

 End Subroutine c_ncsnc
End Interface
!-------------------------------- c_ncabor ------------------------------------
Interface
 Subroutine c_ncabor(ncid, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid 
 Integer(KIND=C_INT), Intent(OUT) :: rcode

 End Subroutine c_ncabor
End Interface
!-------------------------------- c_ncdinq -----------------------------------
Interface
 Subroutine c_ncdinq(ncid, dimid, dimname, size, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , dimid
 Character(KIND=C_CHAR), Intent(OUT) :: dimname(*)
 Integer(KIND=C_INT),    Intent(OUT) :: size, rcode

 End Subroutine c_ncdinq
End Interface
!-------------------------------- c_ncdren ------------------------------------
Interface
 Subroutine c_ncdren(ncid, dimid, dimname, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , dimid
 Character(KIND=C_CHAR), Intent(IN)  :: dimname(*)
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncdren
End Interface
!-------------------------------- c_ncviq -------------------------------------
Interface
 Subroutine c_ncvinq(ncid, varid, varname, datatype, indims, dimarray, &
                      inatts, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid , varid
 Character(KIND=C_CHAR), Intent(INOUT) :: varname(*)
 Integer(KIND=C_INT),    Intent(OUT)   :: datatype ! nc_type var in C
 Integer(KIND=C_INT),    Intent(OUT)   :: dimarray(*)
 Integer(KIND=C_INT),    Intent(OUT)   :: indims, inatts, rcode

 End Subroutine c_ncvinq
End Interface
!-------------------------------- c_ncvpt1 ------------------------------------
Interface
 Subroutine c_ncvpt1(ncid, varid, indices, value, rcode) BIND(C) 
 
 USE ISO_C_BINDING, ONLY: C_INT, C_PTR

 Integer(KIND=C_INT), VALUE       :: ncid , varid
 TYPE(C_PTR),         VALUE       :: indices
 Type(C_PTR),         VALUE       :: value
 Integer(KIND=C_INT), Intent(OUT) :: rcode

 End Subroutine c_ncvpt1
End Interface
!-------------------------------- c_ncvp1c ------------------------------------
Interface
 Subroutine c_ncvp1c(ncid, varid, indices, value, rcode) BIND(C) 
 
 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid
 TYPE(C_PTR),            VALUE       :: indices
 Character(KIND=C_CHAR), Intent(IN)  :: value(*) ! void in C
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncvp1c
End Interface
!-------------------------------- c_ncvpt -------------------------------------
Interface
 Subroutine c_ncvpt(ncid, varid, start, count, value, rcode) BIND(C) 
 
 USE ISO_C_BINDING, ONLY: C_INT, C_PTR

 Integer(KIND=C_INT), VALUE       :: ncid , varid
 Type(C_PTR),         VALUE       :: start, count
 Type(C_PTR),         VALUE       :: value
 Integer(KIND=C_INT), Intent(OUT) :: rcode

 End Subroutine c_ncvpt
End Interface
!-------------------------------- c_ncvptc ------------------------------------
Interface
 Subroutine c_ncvptc(ncid, varid, start, count, value, lenstr, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid, lenstr
 Type(C_PTR),            VALUE       :: start, count
 Character(KIND=C_CHAR), Intent(IN)  :: value(*) ! char in C
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncvptc
End Interface
!-------------------------------- c_ncvptg ------------------------------------
Interface
 Subroutine c_ncvptg(ncid, varid, start, count, strides, imap, value, &
                      rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT, C_PTR

 Integer(KIND=C_INT), VALUE       :: ncid , varid
 Type(C_PTR),         VALUE       :: start, count, strides, imap
 Type(C_PTR),         VALUE       :: value
 Integer(KIND=C_INT), Intent(OUT) :: rcode

 End Subroutine c_ncvptg
End Interface
!-------------------------------- c_ncvpgc ------------------------------------
Interface
 Subroutine c_ncvpgc(ncid, varid, start, count, strides, imap, value, &
                        rcode) BIND(C) 
 
 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid
 Type(C_PTR),            VALUE       :: start, count, strides, imap
 Character(KIND=C_CHAR), Intent(IN)  :: value(*) ! char in C
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncvpgc
End Interface
!-------------------------------- c_ncvgt1 ------------------------------------
Interface
 Subroutine c_ncvgt1(ncid, varid, indices, value, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid
 Type(C_PTR),            VALUE       :: indices
 Character(KIND=C_CHAR), Intent(OUT) :: value(*) ! void in C
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncvgt1
End Interface
!-------------------------------- c_ncvg1c ------------------------------------
Interface
 Subroutine c_ncvg1c(ncid, varid, indices, value, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid , varid
 Type(C_PTR),            VALUE         :: indices
 Character(KIND=C_CHAR), Intent(INOUT) :: value(*) ! char in C
 Integer(KIND=C_INT),    Intent(OUT)   :: rcode

 End Subroutine c_ncvg1c
End Interface
!-------------------------------- c_ncvgt -------------------------------------
Interface
 Subroutine c_ncvgt(ncid, varid, start, count, value, rcode) BIND(C)
 
 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid
 Type(C_PTR),            VALUE       :: start, count
 Character(KIND=C_CHAR), Intent(OUT) :: value(*) ! void in C
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncvgt
End Interface
!-------------------------------- c_ncvgtc ------------------------------------
Interface
 Subroutine c_ncvgtc(ncid, varid, start, count, value, lenstr, rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid , varid, lenstr
 Type(C_PTR),            VALUE         :: start, count
 Character(KIND=C_CHAR), Intent(INOUT) :: value(*) ! char in C
 Integer(KIND=C_INT),    Intent(OUT)   :: rcode

 End Subroutine c_ncvgtc
End Interface
!-------------------------------- c_ncvgtg ------------------------------------
Interface
 Subroutine c_ncvgtg(ncid, varid, start, count, strides, imap, value, &
                      rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid
 Type(C_PTR),            VALUE       :: start, count, strides,  imap
 Character(KIND=C_CHAR), Intent(OUT) :: value(*) ! void in C
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncvgtg
End Interface
!-------------------------------- c_ncvggc ------------------------------------
Interface
 Subroutine c_ncvggc(ncid, varid, start, count, strides, imap, value, &
                      rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_PTR, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid
 Type(C_PTR),            VALUE       :: start, count, strides, imap
 Character(KIND=C_CHAR), Intent(OUT) :: value(*) ! char in C
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncvggc
End Interface
!-------------------------------- c_ncvren ------------------------------------
Interface
 Subroutine c_ncvren(ncid, varid, varname, rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid
 Character(KIND=C_CHAR), Intent(IN)  :: varname(*)
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncvren
End Interface
!-------------------------------- c_ncapt -------------------------------------
Interface
 Subroutine c_ncapt(ncid, varid, attname, datatype, attlen, value, &
                     rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR, C_PTR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid
 Character(KIND=C_CHAR), Intent(IN)  :: attname(*)
 Integer(KIND=C_INT),    VALUE       :: datatype ! nc_type var in C
 Integer(KIND=C_SIZE_T), VALUE       :: attlen
 Type(C_PTR),            VALUE       :: value ! void in C
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncapt
End Interface
!-------------------------------- c_ncaptc ------------------------------------
Interface
 Subroutine c_ncaptc(ncid, varid, attname, datatype, attlen, string, &
                      rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid
 Character(KIND=C_CHAR), Intent(IN)  :: attname(*)
 Integer(KIND=C_INT),    VALUE       :: datatype ! nc_type var in C
 Integer(KIND=C_SIZE_T), VALUE       :: attlen
 Character(KIND=C_CHAR), Intent(IN)  :: string(*) ! char in C
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncaptc
End Interface
!-------------------------------- c_ncainq ------------------------------------
Interface
 Subroutine c_ncainq(ncid, varid, attname, datatype, attlen, rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid
 Character(KIND=C_CHAR), Intent(IN)  :: attname(*)
 Integer(KIND=C_INT),    Intent(OUT) :: datatype ! nc_type var in C
 Integer(KIND=C_INT),    Intent(OUT) :: attlen
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncainq
End Interface
!-------------------------------- c_ncagt -------------------------------------
Interface
 Subroutine c_ncagt(ncid, varid, attname, value, rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid
 Character(KIND=C_CHAR), Intent(IN)  :: attname(*)
 Character(KIND=C_CHAR), Intent(OUT) :: value(*) ! void in C
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncagt
End Interface
!-------------------------------- c_ncagtc ------------------------------------
Interface
 Subroutine c_ncagtc(ncid, varid, attname, value, attlen, rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid, attlen
 Character(KIND=C_CHAR), Intent(IN)  :: attname(*)
 Character(KIND=C_CHAR), Intent(OUT) :: value(*) ! char in C
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncagtc
End Interface
!-------------------------------- c_ncacpy ------------------------------------
Interface
 Subroutine c_ncacpy(inncid, invarid, attname, outncid, outvarid, &
                      rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: inncid , invarid, outncid, outvarid 
 Character(KIND=C_CHAR), Intent(IN)  :: attname(*)
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncacpy
End Interface
!-------------------------------- c_ncanam ------------------------------------
Interface
 Subroutine c_ncanam(ncid, varid, attnum, newname, rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE         :: ncid , varid, attnum
 Character(KIND=C_CHAR), Intent(INOUT) :: newname(*)
 Integer(KIND=C_INT),    Intent(OUT)   :: rcode

 End Subroutine c_ncanam
End Interface
!-------------------------------- c_ncaren ------------------------------------
Interface
 Subroutine c_ncaren(ncid, varid, attnam, newname, rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid
 Character(KIND=C_CHAR), Intent(IN)  :: attnam(*), newname(*)
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncaren
End Interface
!-------------------------------- c_ncadel ------------------------------------
Interface
 Subroutine c_ncadel(ncid, varid, attname, rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT, C_CHAR

 Integer(KIND=C_INT),    VALUE       :: ncid , varid
 Character(KIND=C_CHAR), Intent(IN)  :: attname(*)
 Integer(KIND=C_INT),    Intent(OUT) :: rcode

 End Subroutine c_ncadel
End Interface
!-------------------------------- c_ncsfil ------------------------------------
Interface
 Function c_ncsfil(ncid, fillmode, rcode) BIND(C)

 USE ISO_C_BINDING, ONLY: C_INT

 Integer(KIND=C_INT), VALUE       :: ncid , fillmode 
 Integer(KIND=C_INT), Intent(OUT) :: rcode

 Integer(KIND=C_INT)              :: c_ncsfil

 End Function c_ncsfil
End Interface
!---------------------------------v2data_size ---------------------------------
Interface
 Function v2data_size(datatype) BIND(C) 
!
! New function added to fort-v2compat.c
!
 USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T

 Integer(KIND=C_INT),   VALUE :: datatype
 Integer(KIND=C_SIZE_T)       :: v2data_size

 End Function v2data_size  
End Interface

CONTAINS

Subroutine convert_v2_imap(cncid, cvarid, fmap, cmap, inullp)

! Replacement for f2c_v2imap C function. Uses v2data_size to return
! data size defined for C code. A futher test will be made using
! C interop value.s for FORTRAN side
!
!   USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T
!   USE NETCDF_NC_DATA, ONLY: C_PTRDIFF_T
!   USE netcdf_nc_interfaces, ONLY: NC_CHAR, NC_SHORT, NC_INT, NC_FLOAT, &
!                                   NC_BYTE, NC_DOUBLE, NC_MAX_VAR_DIMS, &

! USE netcdf_nc_interfaces, ONLY: nc_inq_vartype, nc_inq_varndims,     &
!                                 nc_inq_vardimid, nc_inq_dimlen,      &
!                                 NC_NOERR , NC_MAX_VAR_DIMS
   
 Implicit NONE

 Integer(KIND=C_INT),       Intent(IN)    :: cncid, cvarid
 Integer(KIND=C_INT),       Intent(IN)    :: fmap(*)
 Integer(KIND=C_PTRDIFF_T), Intent(INOUT) :: cmap(*)
 Integer,                   Intent(OUT)   :: inullp

 Integer(KIND=C_INT)    :: rank, datatype, cstat1, cstat2, cstat3, cstat4
 Integer(KIND=C_SIZE_T) :: total, length, csize
 Integer(KIND=C_INT)    :: dimids(NC_MAX_VAR_DIMS)
 Integer                :: ii, idim

!
 inullp=0

 cstat1 = nc_inq_vartype(cncid, cvarid, datatype)
 cstat2 = nc_inq_varndims(cncid, cvarid, rank)

! Return if nc_inq_vartype or nc_inq_varndims returns an error
! code. Set inullp to trigger use of NULL pointer in calling 
! routine
 
 If (cstat1/=NC_NOERR) Then
   inullp=1
   Return
 EndIf
 If (cstat2/=NC_NOERR) Then
   inullp=1
   Return
 EndIf
 If (rank <= 0) Then
   inullp=1
   Return
 EndIf

 If (fmap(1)==0) Then ! Special Fortran version 2 sematics
   cstat3 = nc_inq_vardimid(cncid, cvarid, dimids)
   If (cstat3 /= NC_NOERR) Then
     inullp=1
     Return
   EndIf
!
   total = 1
   Loop1: Do ii=1, rank
     idim = rank-ii+1
     cmap(idim) = total
     cstat4 = nc_inq_dimlen(cncid, dimids(idim), length)
     If (cstat4 /= NC_NOERR) Then
       inullp=1
       Exit Loop1
     EndIf
     total = total*length
   EndDo Loop1
   If (inullp==1) Return

 Else ! Standard version 2 format - Use KIND parameters to set size

! Get C data type size using v2data_size. Unfortunately, the F03
! standard didn't specify a C_SIZEOF function. This will be
! remedied in the next upgrade to FORTRAN (2008) but for now
! we will rely on a C function to provide the value

   csize = v2data_size(datatype) 
   If (csize <= 0) Then
      inullp=1
      Return 
   EndIf

   cmap(1:rank) = fmap(rank:1:-1) / csize

 EndIf

End Subroutine convert_v2_imap

!-------------------- End module_netcdf_fortv2_c_interfaces -------------------
End Module netcdf_fortv2_c_interfaces
