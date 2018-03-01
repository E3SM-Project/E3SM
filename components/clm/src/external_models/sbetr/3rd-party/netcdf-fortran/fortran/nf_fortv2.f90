!                  Netcdf Version 2 Fortran API 

! These routines replace the cfortran.h defined functions from fort-v2compat.c
! C_CHAR strings and C_PTR types are used as pass-through variables for
! functions that support multiple data types as void pointers in the C routines
!

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

! Version 1: May   2006 - Initial version 2 interfaces
! Version 2: April 2009 - Refactored to pass value data using C_CHAR and C_PTR
!                         strings and pointers and updated for NetCDF 4.0.1
! Version 3: April 2010 - updated for NetCDF 4.1.1

! ------------------------------- ncpopt -------------------------------------- 
 Subroutine ncpopt(ncopts)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer, Intent(IN) :: ncopts
 
 Integer(KIND=C_INT) :: cncopts

 cncopts = ncopts

 Call c_ncpopt(cncopts)

 End Subroutine ncpopt 
! ------------------------------- ncgopt -------------------------------------- 
 Subroutine ncgopt(ncopts)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer, Intent(INOUT) :: ncopts

 Integer(KIND=C_INT)    :: cncopts

 cncopts = 0 

 Call c_ncgopt(cncopts)

 ncopts = cncopts

 End Subroutine ncgopt 
! ------------------------------- nccre --------------------------------------- 
 Function nccre(filename, cmode, rcode) RESULT(ncid)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: filename
 Integer,          Intent(IN)  :: cmode
 Integer,          Intent(OUT) :: rcode

 Integer                       :: ncid

 Character(LEN=(LEN(filename)+1)) :: cfilename
 Integer                          :: ilen
 Integer(KIND=C_INT)              :: ccmode, crcode, cncid 

 ccmode = cmode
 cncid  = 0
 rcode  = 0
 crcode = 0

! check for a null character on end of filename

 cfilename = addCNullChar(filename, ilen)
 
 cncid = c_nccre(cfilename(1:ilen+1), ccmode, crcode )

 rcode = crcode 
 ncid  = cncid

 End Function nccre
! ------------------------------- ncopn --------------------------------------- 
 Function ncopn(filename, rwmode, rcode) RESULT(ncid)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: filename
 Integer,          Intent(IN)  :: rwmode
 Integer,          Intent(OUT) :: rcode

 Integer                       :: ncid

 Character(LEN=(LEN(filename)+1)) :: cfilename
 Integer                          :: ilen
 Integer(KIND=C_INT)              :: crwmode, crcode, cncid 

 crwmode = rwmode
 rcode   = 0
 crcode  = 0
 cncid   = 0

! check for a null character on end of filename

 cfilename = addCNullChar(filename, ilen)
 
 cncid = c_ncopn(cfilename(1:ilen+1), crwmode, crcode )

 rcode = crcode 
 ncid  = cncid

 End Function ncopn
! ------------------------------- ncddef -------------------------------------- 
 Function ncddef(ncid, dimname, dimlen, rcode) RESULT(ndimid)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: dimname 
 Integer,          Intent(IN)  :: ncid, dimlen
 Integer,          Intent(OUT) :: rcode

 Integer                       :: ndimid

 Character(LEN=(LEN(dimname)+1)) :: cdimname
 Integer                         :: ilen
 Integer(KIND=C_INT)             :: cncid, cdimlen, cndimid, crcode

 cncid   = ncid
 cdimlen = dimlen
 cndimid = 0
 rcode   = 0
 ndimid  = 0

! check for a null character on end of dimname

 cdimname = addCNullChar(dimname, ilen)
 
 cndimid = c_ncddef(cncid, cdimname(1:ilen+1), cdimlen, crcode )
 
 rcode  = crcode 
 ndimid = cndimid

 End Function ncddef
! ------------------------------- ncdid --------------------------------------- 
 Function ncdid(ncid, dimname, rcode) RESULT(ndimid)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: dimname 
 Integer,          Intent(IN)  :: ncid
 Integer,          Intent(OUT) :: rcode

 Integer                       :: ndimid

 Character(LEN=(LEN(dimname)+1)) :: cdimname
 Integer                         :: ilen
 Integer(KIND=C_INT)             :: cncid, crcode, cndimid

 cncid   = ncid
 cndimid = 0
 rcode   = 0
 ndimid  = 0

! check for a null character on end of dimname

 cdimname = addCNullChar(dimname, ilen)
 
 cndimid = c_ncdid(cncid, cdimname(1:ilen+1), crcode )

 rcode  = crcode 
 ndimid = cndimid

 End Function ncdid
! ------------------------------- ncvdef -------------------------------------- 
 Function ncvdef(ncid, varname, vartype, nvdims, vdims, rcode) RESULT(nvarid)

 USE netcdf_nc_interfaces, ONLY : NC_MAX_DIMS
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: varname 
 Integer,          Intent(IN)  :: ncid, vartype, nvdims
 Integer,          Intent(IN)  :: vdims(*)
 Integer,          Intent(OUT) :: rcode

 Integer                       :: nvarid

 Character(LEN=(LEN(varname)+1)) :: cvarname
 Integer                         :: ilen
 Integer(KIND=C_INT)             :: cncid, crcode, cnvdims, cvartype, cnvarid
 Integer(KIND=C_INT)             :: cvdims(NC_MAX_DIMS)

 cncid    = ncid
 cnvdims  = nvdims 
 cvartype = vartype
 crcode   = 0
 rcode    = 0
 nvarid   = 0
 cnvarid  = 0
 
! check for a null character on end of varname

 cvarname = addCNullChar(varname, ilen)
 
 ! mimic f2c_dimids

 cvdims = 0
 If (nvdims > 0) Then
   cvdims(1:nvdims) = vdims(nvdims:1:-1) - 1
 EndIf

 cnvarid = c_ncvdef(cncid, cvarname(1:ilen+1), cvartype, &
                    cnvdims, cvdims, crcode )

 rcode  = crcode 
 nvarid = cnvarid

 End Function ncvdef
! ------------------------------- ncvid --------------------------------------- 
 Function ncvid(ncid, varname, rcode) RESULT(nvarid)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: varname 
 Integer,          Intent(IN)  :: ncid
 Integer,          Intent(OUT) :: rcode

 Integer                       :: nvarid

 Character(LEN=(LEN(varname)+1)) :: cvarname
 Integer                         :: ilen
 Integer(KIND=C_INT)             :: cncid, crcode, cnvarid

 cncid   = ncid
 crcode  = 0
 rcode   = 0
 nvarid  = 0
 cnvarid = 0
 
! check for a null character on end of varname

 cvarname = addCNullChar(varname, ilen)
 
 cnvarid = c_ncvid(cncid, cvarname(1:ilen+1), crcode)

 rcode  = crcode 
 nvarid = cnvarid

 End Function ncvid
! ------------------------------- nctlen -------------------------------------- 
 Function nctlen(datatype, rcode) RESULT(nvarlen)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: datatype 
 Integer, Intent(OUT) :: rcode

 Integer              :: nvarlen

 Integer(KIND=C_INT) :: crcode, cnvarlen, cdtype

 cdtype   = datatype
 crcode   = 0
 rcode    = 0
 nvarlen  = 0
 cnvarlen = 0
 
 cnvarlen = c_nctlen(cdtype, crcode)

 rcode   = crcode 
 nvarlen = cnvarlen

 End Function nctlen
! ------------------------------- ncclos ------------------------------------- 
 Subroutine ncclos(ncid, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid 
 Integer, Intent(OUT) :: rcode

 Integer(KIND=C_INT) :: crcode, cncid

 cncid   = ncid
 crcode  = 0
 rcode   = 0
 
 Call c_ncclos(cncid, crcode)

 rcode = crcode 

 End Subroutine ncclos 
! ------------------------------- ncredf ------------------------------------- 
 Subroutine ncredf(ncid, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid 
 Integer, Intent(OUT) :: rcode

 Integer(KIND=C_INT) :: crcode, cncid

 cncid   = ncid
 crcode  = 0
 rcode   = 0
 
 Call c_ncredf(cncid, crcode)

 rcode = crcode 

 End Subroutine ncredf
! ------------------------------- ncendf -------------------------------------- 
 Subroutine ncendf(ncid, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid 
 Integer, Intent(OUT) :: rcode

 Integer(KIND=C_INT) :: crcode, cncid

 cncid   = ncid
 crcode  = 0
 rcode   = 0
 
 Call c_ncendf(cncid, crcode)

 rcode = crcode 

 End Subroutine ncendf
! ------------------------------- ncinq --------------------------------------- 
 Subroutine ncinq(ncid, ndims, nvars, natts, recdim, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid 
 Integer, Intent(OUT) :: ndims, nvars, natts, recdim, rcode

 Integer(KIND=C_INT) :: crcode, cncid, cndims, cnvars, cnatts, crecdim

 cncid   = ncid
 crcode  = 0
 rcode   = 0
 cndims  = 0
 cnvars  = 0
 cnatts  = 0
 ndims   = 0
 nvars   = 0
 natts   = 0
  
 Call c_ncinq(cncid, cndims, cnvars, cnatts, crecdim, crcode)

 ndims = cndims
 nvars = cnvars
 natts = cnatts
 If (crecdim == -1) Then ! no unlimited dimension
   recdim = -1
 Else
   recdim = crecdim + 1  ! shift by plus one for FORTRAN
 EndIf 

 rcode = crcode 

 End Subroutine ncinq
! ------------------------------- ncsnc --------------------------------------- 
 Subroutine ncsnc(ncid, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid 
 Integer, Intent(OUT) :: rcode

 Integer(KIND=C_INT) :: crcode, cncid

 cncid   = ncid
 crcode  = 0
 rcode   = 0
 
 Call c_ncsnc(cncid, crcode)

 rcode = crcode 

 End Subroutine ncsnc
! ------------------------------- ncabor -------------------------------------- 
 Subroutine ncabor(ncid, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid 
 Integer, Intent(OUT) :: rcode

 Integer(KIND=C_INT) :: crcode, cncid

 cncid   = ncid
 crcode  = 0
 rcode   = 0
 
 Call c_ncabor(cncid, crcode)

 rcode = crcode 

 End Subroutine ncabor
! ------------------------------- ncdinq -------------------------------------- 
 Subroutine ncdinq(ncid, dimid, dimname, dimlen, rcode)

 USE netcdf_nc_interfaces, ONLY: NC_MAX_NAME
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, dimid
 Character(LEN=*), Intent(OUT) :: dimname
 Integer,          Intent(OUT) :: dimlen, rcode

 Integer(KIND=C_INT)            :: cncid, crcode, cdimlen, cdimid
 Character(LEN=(NC_MAX_NAME+1)) :: cdimname
 Integer                        :: ilen

 cncid    = ncid
 cdimid   = dimid - 1
 crcode   = 0
 rcode    = 0
 cdimlen  = 0
 cdimname = REPEAT(" ", LEN(cdimname))
 ilen = LEN(dimname)
 
 Call c_ncdinq(cncid, cdimid, cdimname, cdimlen, crcode)

! check for a null character on end of cdimname

 dimname = stripCNullChar(cdimname, ilen)
 
 dimlen          = cdimlen
 rcode           = crcode

 End Subroutine ncdinq
! ------------------------------- ncdren -------------------------------------- 
 Subroutine ncdren(ncid, dimid, dimname, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: dimname
 Integer,          Intent(IN)  :: ncid, dimid
 Integer,          Intent(OUT) :: rcode

 Character(LEN=(LEN(dimname)+1)) :: cdimname
 Integer(KIND=C_INT)             :: cncid, crcode, cdimid
 Integer                         :: ilen

 cncid  = ncid
 cdimid = dimid - 1
 crcode = 0
 rcode  = 0

! check for a null character on end of dimname

 cdimname = addCNullChar(dimname, ilen)
 
 Call c_ncdren(cncid, cdimid, cdimname(1:ilen+1), crcode)

 rcode = crcode 

 End Subroutine ncdren
! ------------------------------- ncvinq -------------------------------------- 
 Subroutine ncvinq(ncid, varid, varname, vartype, nvdims, vdims, &
                   nvatts, rcode)

 USE netcdf_nc_interfaces, ONLY: NC_MAX_DIMS, NC_MAX_NAME
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,          Intent(IN)    :: ncid, varid
 Character(LEN=*), Intent(INOUT) :: varname
 Integer,          Intent(OUT)   :: vartype, nvdims, nvatts, rcode
 Integer,          Intent(INOUT) :: vdims(*)

 Integer(KIND=C_INT)          :: cncid, crcode, cvarid, cvartype, cnvdims, &
                                 cnvatts
 Integer(KIND=C_INT)          :: cvdims(NC_MAX_DIMS)
 Character(LEN=NC_MAX_NAME+1) :: cvarname
 Integer                      :: ilen

 cncid    = ncid
 cvarid   = varid - 1
 crcode   = 0
 rcode    = 0
 cvdims   = 0
 cvdims   = 0
 vartype  = 0
 nvdims   = 0
 nvatts   = 0
 cnvdims  = 0
 cvdims   = 0
 cnvatts  = 0
 cvartype = 0
 cvarname = REPEAT(" ", LEN(cvarname))
 ilen = LEN(varname)

 Call c_ncvinq(cncid, cvarid, cvarname, cvartype, cnvdims, cvdims, cnvatts, &
               crcode)      

 nvdims  = cnvdims
 vartype = cvartype
 nvatts  = cnvatts 
 rcode   = crcode 

! strip C null character from cvarname
   
 varname = stripCNullChar(cvarname, ilen)
 
! convert C dimids to FORTRAN order and rank
! Replaces call to c2f_dimids in C code

 If (nvdims > 0) Then
   vdims(1:nvdims) = cvdims(nvdims:1:-1) + 1
 End If

 End Subroutine ncvinq
! ------------------------------- ncvpt1 -------------------------------------- 
 Subroutine ncvpt1(ncid, varid, mindex, values, rcode) 

 USE netcdf_nc_interfaces, ONLY: NC_MAX_DIMS, NC_NOERR, nc_inq_varndims
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,                Intent(IN)          :: ncid, varid
 Integer,                Intent(IN)          :: mindex(*)
 Character(KIND=C_CHAR), Intent(IN), TARGET  :: values(*)
 Integer,                Intent(OUT)         :: rcode

 Integer(KIND=C_INT)            :: cncid, crcode, cvarid, cstatus, cndims
 Integer(KIND=C_SIZE_T), TARGET :: cmindex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cmindexptr
 Type(C_PTR)                    :: cvaluesptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1
 crcode  = 0
 rcode   = 0
 cmindex = 0
 cndims  = 0
 ndims   = 0

 cstatus = nc_inq_varndims(cncid, cvarid, cndims)

 cmindexptr = C_NULL_PTR
 ndims      = cndims
 
 If (cstatus == NC_NOERR) Then ! mimic f2c_coords in C code
   If (ndims > 0) Then
     cmindex(1:ndims) = mindex(ndims:1:-1) - 1
   Endif
   cmindexptr = C_LOC(cmindex)
 Endif
 
 cvaluesptr = C_LOC(values)

 Call c_ncvpt1(cncid, cvarid, cmindexptr, cvaluesptr, crcode)

 rcode = crcode

 End Subroutine ncvpt1
! ------------------------------- ncvp1c -------------------------------------- 
 Subroutine ncvp1c(ncid, varid, mindex, strings, rcode) 

 USE netcdf_nc_interfaces, ONLY: NC_MAX_DIMS, NC_NOERR, nc_inq_varndims
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: mindex(*)
 Character(LEN=*), Intent(IN)  :: strings
 Integer,          Intent(OUT) :: rcode

 Integer(KIND=C_INT)            :: cncid, crcode, cvarid, cstatus, cndims
 Integer(KIND=C_SIZE_T), TARGET :: cmindex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cmindexptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1
 crcode  = 0
 rcode   = 0
 cmindex = 0
 cndims  = 0
 ndims   = 0

 cstatus = nc_inq_varndims(cncid, cvarid, cndims)

 cmindexptr = C_NULL_PTR
 ndims     = cndims 

 If (cstatus == NC_NOERR) Then ! mimic f2c_coords in C code
   If (ndims > 0) Then
     cmindex(1:ndims) = mindex(ndims:1:-1) - 1
   Endif
   cmindexptr = C_LOC(cmindex)
 Endif

 Call c_ncvp1c(cncid, cvarid, cmindexptr, strings, crcode)

 rcode = crcode

 End Subroutine ncvp1c
! ------------------------------- ncvpt --------------------------------------- 
 Subroutine ncvpt(ncid, varid, start, counts, values, rcode) 

 USE netcdf_nc_interfaces, ONLY: NC_MAX_DIMS, NC_NOERR, nc_inq_varndims
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,                Intent(IN)          :: ncid, varid
 Integer,                Intent(IN)          :: start(*), counts(*)
 Character(KIND=C_CHAR), Intent(IN), TARGET  :: values(*)
 Integer,                Intent(OUT)         :: rcode

 Integer(KIND=C_INT)            :: cncid, crcode, cvarid, cstatus, cndims
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Type(C_PTR)                    :: cvaluesptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1
 crcode  = 0
 rcode   = 0
 cstart  = 0
 ccounts = 0
 cndims  = 0
 ndims   = 0

 cstatus = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims 

 If (cstatus == NC_NOERR) Then ! mimic f2c_coords, etc. in C code
   If (ndims > 0) Then
     cstart(1:ndims)  = start(ndims:1:-1) - 1
     ccounts(1:ndims) = counts(ndims:1:-1)
   Endif
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts) 
 Endif

 cvaluesptr = C_LOC(values)

 Call c_ncvpt(cncid, cvarid, cstartptr, ccountsptr, cvaluesptr, crcode)

 rcode = crcode

 End Subroutine ncvpt
! ------------------------------- ncvptc--------------------------------------- 
 Subroutine ncvptc(ncid, varid, start, counts, strings, lenstr, rcode) 

 USE netcdf_nc_interfaces, ONLY: NC_MAX_DIMS, NC_NOERR, nc_inq_varndims
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,          Intent(IN)    :: ncid, varid, lenstr
 Integer,          Intent(IN)    :: start(*), counts(*)
 Character(LEN=*), Intent(INOUT) :: strings
 Integer,          Intent(OUT)   :: rcode

 Integer(KIND=C_INT)            :: cncid, crcode, cvarid, cstatus, cndims, &
                                   clenstr
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1
 clenstr = lenstr
 crcode  = 0
 rcode   = 0
 cstart  = 0
 ccounts = 0
 cndims  = 0
 ndims   = 0

 cstatus = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims
 
 If (cstatus == NC_NOERR) Then ! mimic f2c_coords, etc. in C code
   If (ndims > 0) Then
     cstart(1:ndims)  = start(ndims:1:-1) - 1
     ccounts(1:ndims) = counts(ndims:1:-1)
   Endif
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 Endif

 Call c_ncvptc(cncid, cvarid, cstartptr, ccountsptr, strings(1:lenstr),&
               clenstr, crcode)

 rcode = crcode

 End Subroutine ncvptc
! ------------------------------- ncvptg -------------------------------------- 
 Subroutine ncvptg(ncid, varid, start, counts, strides, imap, values, &
                   rcode) 

 USE netcdf_nc_interfaces, ONLY: NC_MAX_DIMS, NC_NOERR, nc_inq_varndims
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,                Intent(IN)          :: ncid, varid
 Integer,                Intent(IN)          :: start(*), counts(*), &
                                                strides(*), imap(*)
 Character(KIND=C_CHAR), Intent(IN), TARGET  :: values(*)
 Integer,                Intent(OUT)         :: rcode

 Integer(KIND=C_INT)               :: cncid, crcode, cvarid, cstatus, cndims
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cimap(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cimapptr, &
                                      cstridesptr
 Type(C_PTR)                       :: cvaluesptr
 Integer                           :: ndims, inullp

 cncid   = ncid
 cvarid  = varid - 1
 crcode  = 0
 rcode   = 0
 cstart  = 0
 ccounts = 0
 cndims  = 0
 ndims   = 0
 inullp  = 0

 Call convert_v2_imap(cncid, cvarid, imap, cimap, inullp)
 
 cstatus = nc_inq_varndims(cncid, cvarid, cndims)

 ndims       = cndims 
 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cimapptr    = C_LOC(cimap)
 If (inullp /= 0) cimapptr = C_NULL_PTR
 
 If (cstatus == NC_NOERR) Then ! mimic f2c_coords, etc. in C code
   If (ndims > 0) Then
     cstart(1:ndims)   = start(ndims:1:-1) - 1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1) - 1
   Endif
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
 Endif

 cvaluesptr = C_LOC(values)

 Call c_ncvptg(cncid, cvarid, cstartptr, ccountsptr, cstridesptr, &
               cimapptr, cvaluesptr, crcode)

 rcode = crcode

 End Subroutine ncvptg
! ------------------------------- ncvpgc -------------------------------------- 
 Subroutine ncvpgc(ncid, varid, start, counts, strides, imap, string, rcode) 

 USE netcdf_nc_interfaces, ONLY: NC_MAX_DIMS, NC_NOERR, nc_inq_varndims
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(IN)  :: start(*), counts(*), strides(*), imap(*)
 Character(LEN=*), Intent(IN)  :: string
 Integer,          Intent(OUT) :: rcode

 Integer(KIND=C_INT)               :: cncid, crcode, cvarid, cstatus, cndims
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cimap(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cimapptr
 Integer                           :: ndims, inullp

 cncid   = ncid
 cvarid  = varid - 1
 crcode  = 0
 rcode   = 0
 cstart  = 0
 ccounts = 0
 cndims  = 0
 ndims   = 0
 inullp  = 0

 Call convert_v2_imap(cncid, cvarid, imap, cimap, inullp)

 cstatus = nc_inq_varndims(cncid, cvarid, cndims)

 ndims       = cndims 
 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cimapptr    = C_LOC(cimap)
 If (inullp /= 0) cimapptr = C_NULL_PTR

 If (cstatus == NC_NOERR) Then ! mimic f2c_coords, etc. in C code
   If (ndims > 0) Then
     cstart(1:ndims)   = start(ndims:1:-1) - 1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1) - 1
   Endif
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
 Endif

 Call c_ncvpgc(cncid, cvarid, cstartptr, ccountsptr, cstridesptr, &
               cimapptr, string, crcode)

 rcode = crcode

 End Subroutine ncvpgc
! ------------------------------- ncvgt1 -------------------------------------- 
 Subroutine ncvgt1(ncid, varid, mindex, values, rcode) 

 USE netcdf_nc_interfaces, ONLY: NC_MAX_DIMS, NC_NOERR, nc_inq_varndims
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,                Intent(IN)  :: ncid, varid
 Integer,                Intent(IN)  :: mindex(*)
 Character(KIND=C_CHAR), Intent(OUT) :: values(*)
 Integer,                Intent(OUT) :: rcode

 Integer(KIND=C_INT)            :: cncid, crcode, cvarid, cstatus, cndims
 Integer(KIND=C_SIZE_T), TARGET :: cmindex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cmindexptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1
 crcode  = 0
 rcode   = 0
 cmindex = 0
 cndims  = 0
 ndims   = 0

 cstatus = nc_inq_varndims(cncid, cvarid, cndims)

 cmindexptr = C_NULL_PTR
 ndims      = cndims
 
 If (cstatus == NC_NOERR) Then ! mimic f2c_coords in C code
   If (ndims > 0) Then
     cmindex(1:ndims) = mindex(ndims:1:-1) - 1
   Endif
   cmindexptr = C_LOC(cmindex)
 Endif

 Call c_ncvgt1(cncid, cvarid, cmindexptr, values, crcode)

 rcode = crcode

 End Subroutine ncvgt1
! ------------------------------- ncvg1c -------------------------------------- 
 Subroutine ncvg1c(ncid, varid, mindex, string, rcode) 

 USE netcdf_nc_interfaces, ONLY: NC_MAX_DIMS, NC_NOERR, nc_inq_varndims
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,          Intent(IN)    :: ncid, varid
 Integer,          Intent(IN)    :: mindex(*)
 Character(LEN=*), Intent(INOUT) :: string
 Integer,          Intent(OUT)   :: rcode

 Integer(KIND=C_INT)            :: cncid, crcode, cvarid, cstatus, cndims
 Integer(KIND=C_SIZE_T), TARGET :: cmindex(NC_MAX_DIMS)
 Type(C_PTR)                    :: cmindexptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1
 crcode  = 0
 rcode   = 0
 cmindex = 0
 cndims  = 0
 ndims   = 0

 cstatus = nc_inq_varndims(cncid, cvarid, cndims)

 cmindexptr = C_NULL_PTR
 ndims      = cndims
 
 If (cstatus == NC_NOERR) Then ! mimic f2c_coords in C code
   If (ndims > 0) Then
     cmindex(1:ndims) = mindex(ndims:1:-1) - 1
   Endif
   cmindexptr = C_LOC(cmindex)
 Endif

 Call c_ncvg1c(cncid, cvarid, cmindexptr, string, crcode)

 rcode = crcode

 End Subroutine ncvg1c
! ------------------------------- ncvgt --------------------------------------- 
 Subroutine ncvgt(ncid, varid, start, counts, values, rcode) 

 USE netcdf_nc_interfaces, ONLY: NC_MAX_DIMS, NC_NOERR, nc_inq_varndims
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,                Intent(IN)  :: ncid, varid
 Integer,                Intent(IN)  :: start(*), counts(*)
 Character(KIND=C_CHAR), Intent(OUT) :: values(*)
 Integer,                Intent(OUT) :: rcode

 Integer(KIND=C_INT)            :: cncid, crcode, cvarid, cstatus, cndims
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Integer                        :: ndims

 cncid   = ncid
 cvarid  = varid - 1
 crcode  = 0
 rcode   = 0
 cstart  = 0
 ccounts = 0
 cndims  = 0
 ndims   = 0

 cstatus = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims
 
 If (cstatus == NC_NOERR) Then ! mimic f2c_coords, etc. in C code
   If (ndims > 0) Then
     cstart(1:ndims)  = start(ndims:1:-1) - 1
     ccounts(1:ndims) = counts(ndims:1:-1)
   Endif
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 Endif

 Call c_ncvgt(cncid, cvarid, cstartptr, ccountsptr, values, crcode)

 rcode = crcode

 End Subroutine ncvgt
! ------------------------------- ncvgtc -------------------------------------- 
 Subroutine ncvgtc(ncid, varid, start, counts, string, lenstr, rcode) 

 USE netcdf_nc_interfaces, ONLY: NC_MAX_DIMS, NC_NOERR, nc_inq_varndims
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,          Intent(IN)    :: ncid, varid, lenstr
 Integer,          Intent(IN)    :: start(*), counts(*)
 Character(LEN=*), Intent(INOUT) :: string
 Integer,          Intent(OUT)   :: rcode

 Integer(KIND=C_INT)            :: cncid, crcode, cvarid, cstatus, cndims, &
                                   clenstr
 Integer(KIND=C_SIZE_T), TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Type(C_PTR)                    :: cstartptr, ccountsptr
 Character(LEN=lenstr+1)        :: cstring
 Integer                        :: ndims, slen

 cncid   = ncid
 cvarid  = varid - 1
 clenstr = lenstr
 crcode  = 0
 rcode   = 0
 cstart  = 0
 ccounts = 0
 cndims  = 0
 ndims   = 0
 string  = REPEAT(" ", LEN(string))
 cstring = REPEAT(" ", LEN(cstring))

 cstatus = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr  = C_NULL_PTR
 ccountsptr = C_NULL_PTR
 ndims      = cndims
 
 If (cstatus == NC_NOERR) Then ! mimic f2c_coords, etc. in C code
   If (ndims > 0) Then
     cstart(1:ndims)  = start(ndims:1:-1) - 1
     ccounts(1:ndims) = counts(ndims:1:-1)
   Endif
   cstartptr  = C_LOC(cstart)
   ccountsptr = C_LOC(ccounts)
 Endif

 Call c_ncvgtc(cncid, cvarid, cstartptr, ccountsptr, cstring, clenstr, crcode)

 If (LEN(string) >= lenstr) Then
   string(1:lenstr) = cstring(1:lenstr)
 Else
   slen           = LEN(string)
   string(1:slen) = cstring(1:slen)
 EndIf

 rcode = crcode

 End Subroutine ncvgtc
! ------------------------------- ncvgtg -------------------------------------- 
 Subroutine ncvgtg(ncid, varid, start, counts, strides, imap, values, &
                   rcode) 

 USE netcdf_nc_interfaces, ONLY: NC_MAX_DIMS, NC_NOERR, nc_inq_varndims
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,                Intent(IN)  :: ncid, varid
 Integer,                Intent(IN)  :: start(*), counts(*), strides(*), imap(*)
 Character(KIND=C_CHAR), Intent(OUT) :: values(*)
 Integer,                Intent(OUT) :: rcode

 Integer(KIND=C_INT)               :: cncid, crcode, cvarid, cstatus, cndims
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cimap(NC_MAX_DIMS)
 Type(C_PTR)                       :: cstartptr, ccountsptr, cimapptr, &
                                      cstridesptr
 Integer                           :: ndims, inullp

 cncid   = ncid
 cvarid  = varid - 1
 crcode  = 0
 rcode   = 0
 cstart  = 0
 ccounts = 0
 cndims  = 0
 inullp  = 0

 Call convert_v2_imap(cncid, cvarid, imap, cimap, inullp)
 cstatus = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cimapptr    = C_LOC(cimap)
 ndims       = cndims 
 If (inullp /= 0) cimapptr = C_NULL_PTR

 If (cstatus == NC_NOERR) Then ! mimic f2c_coords, etc. in C code
   If (ndims > 0) Then
     cstart(1:ndims)   = start(ndims:1:-1) - 1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1) - 1
   Endif
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
 Endif

 Call c_ncvgtg(cncid, cvarid, cstartptr, ccountsptr, cstridesptr, &
               cimapptr, values, crcode)

 rcode = crcode

 End Subroutine ncvgtg
! ------------------------------- ncvggc -------------------------------------- 
 Subroutine ncvggc(ncid, varid, start, counts, strides, imap, string, rcode) 

 USE netcdf_nc_interfaces, ONLY: NC_MAX_DIMS, NC_NOERR, nc_inq_varndims
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer,          Intent(IN)    :: ncid, varid
 Integer,          Intent(IN)    :: start(*), counts(*), strides(*), imap(*)
 Character(LEN=*), Intent(INOUT) :: string
 Integer,          Intent(OUT)   :: rcode

 Integer(KIND=C_INT)               :: cncid, crcode, cvarid, cstatus, cndims
 Integer(KIND=C_SIZE_T),    TARGET :: cstart(NC_MAX_DIMS), ccounts(NC_MAX_DIMS)
 Integer(KIND=C_PTRDIFF_T), TARGET :: cstrides(NC_MAX_DIMS), cimap(NC_MAX_DIMS)
 Character(LEN=(LEN(string)+1))    :: cstring
 Type(C_PTR)                       :: cstartptr, ccountsptr, cstridesptr, &
                                      cimapptr
 Integer                           :: ndims, inullp,slen

 cncid   = ncid
 cvarid  = varid - 1
 crcode  = 0
 rcode   = 0
 cstart  = 0
 ccounts = 0
 cndims  = 0
 ndims   = 0
 inullp  = 0
 string  = REPEAT(" ", LEN(string))
 cstring = REPEAT(" ", LEN(cstring)) 

 Call convert_v2_imap(cncid, cvarid, imap, cimap, inullp)

 cstatus = nc_inq_varndims(cncid, cvarid, cndims)

 cstartptr   = C_NULL_PTR
 ccountsptr  = C_NULL_PTR
 cstridesptr = C_NULL_PTR
 cimapptr    = C_LOC(cimap)
 ndims       = cndims 
 If (inullp /= 0) cimapptr = C_NULL_PTR

 If (cstatus == NC_NOERR) Then ! mimic f2c_coords, etc. in C code
   If (ndims > 0) Then
     cstart(1:ndims)   = start(ndims:1:-1) - 1
     ccounts(1:ndims)  = counts(ndims:1:-1)
     cstrides(1:ndims) = strides(ndims:1:-1) - 1
   Endif
   cstartptr   = C_LOC(cstart)
   ccountsptr  = C_LOC(ccounts)
   cstridesptr = C_LOC(cstrides)
 Endif

 Call c_ncvggc(cncid, cvarid, cstartptr, ccountsptr, cstridesptr, &
               cimapptr, cstring, crcode)

 slen           = LEN(string)
 string(1:slen) = cstring(1:slen)

 rcode = crcode

 End Subroutine ncvggc
!-------------------------------- ncvren --------------------------------------
 Subroutine ncvren(ncid, varid, newnam, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: newnam
 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(OUT) :: rcode

 Character(LEN=(LEN(newnam)+1)) :: cnewnam
 Integer(KIND=C_INT)            :: cncid, cvarid, crcode
 Integer                        :: ilen

 cncid  = ncid
 cvarid = varid - 1
 rcode  = 0

! check for a null character on end of newnam

 cnewnam = addCNullChar(newnam, ilen)
 
 Call c_ncvren(cncid, cvarid, cnewnam(1:ilen+1), crcode )

 rcode = crcode 

 End Subroutine ncvren
!-------------------------------- ncapt ---------------------------------------
 Subroutine ncapt(ncid, varid, attnam, attype, attlen, value, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*),       Intent(IN)          :: attnam
 Integer,                Intent(IN)          :: ncid, varid, attype, attlen
 Character(KIND=C_CHAR), Intent(IN), TARGET  :: value(*)
 Integer,                Intent(OUT)         :: rcode

 Integer(KIND=C_INT)            :: cncid, cvarid, cattype, crcode
 Integer(KIND=C_SIZE_T)         :: cattlen
 Type(C_PTR)                    :: cvalueptr
 Character(LEN=(LEN(attnam)+1)) :: cattnam
 Integer                        :: ilen

 cncid   = ncid
 cvarid  = varid - 1
 cattype = attype
 cattlen = attlen
 rcode   = 0

! check for a null character on end of attname

 cattnam = addCNullChar(attnam, ilen)
 
 cvalueptr = C_LOC(value)
 
 Call c_ncapt(cncid, cvarid, cattnam(1:ilen+1), cattype, &
              cattlen, cvalueptr, crcode )

 rcode = crcode 

 End Subroutine ncapt
!-------------------------------- ncaptc ------------------------------------
 Subroutine ncaptc(ncid, varid, attnam, attype, lenstr, string, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: attnam
 Integer,          Intent(IN)  :: ncid, varid, attype, lenstr
 Character(LEN=*), Intent(IN)  :: string 
 Integer,          Intent(OUT) :: rcode

 Integer(KIND=C_INT)            :: cncid, cvarid, cattype, crcode
 Integer(KIND=C_SIZE_T)         :: clenstr
 Character(LEN=(LEN(attnam)+1)) :: cattnam
 Integer                        :: ilen

 cncid   = ncid
 cvarid  = varid - 1
 cattype = attype
 clenstr = lenstr
 rcode   = 0

! check for a null character on end of attname

 cattnam = addCNullChar(attnam, ilen)
 
 Call c_ncaptc(cncid, cvarid, cattnam(1:ilen+1), cattype, &
               clenstr, string, crcode )

 rcode = crcode 

 End Subroutine ncaptc
!-------------------------------- ncainq --------------------------------------
 Subroutine ncainq(ncid, varid, attnam, attype, attlen, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: attnam
 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(OUT) :: attype, attlen, rcode

 Integer(KIND=C_INT)            :: cncid, cvarid, cattype, crcode, cattlen
 Character(LEN=(LEN(attnam)+1)) :: cattnam
 Integer                        :: ilen

 cncid   = ncid
 cvarid  = varid - 1
 cattype = 0
 cattlen = 0
 rcode   = 0

! check for a null character on end of attnam

 cattnam = addCNullChar(attnam, ilen)

 Call c_ncainq(cncid, cvarid, cattnam(1:ilen+1), cattype, &
             cattlen, crcode )

 attype = cattype
 attlen = cattlen
 rcode  = crcode

 End Subroutine ncainq
!-------------------------------- ncagt ---------------------------------------
 Subroutine ncagt(ncid, varid, attnam, values, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*),       Intent(IN)  :: attnam
 Integer,                Intent(IN)  :: ncid, varid
 Character(KIND=C_CHAR), Intent(OUT) :: values(*)
 Integer,                Intent(OUT) :: rcode

 Integer(KIND=C_INT)            :: cncid, cvarid, crcode
 Character(LEN=(LEN(attnam)+1)) :: cattnam
 Integer                        :: ilen

 cncid  = ncid
 cvarid = varid - 1
 rcode  = 0

! check for a null character on end of attnam

 cattnam = addCNullChar(attnam, ilen)
 
 Call c_ncagt(cncid, cvarid, cattnam(1:ilen+1), values, crcode)

 rcode = crcode

 End Subroutine ncagt
!-------------------------------- ncagtc --------------------------------------
 Subroutine ncagtc(ncid, varid, attnam, string, lenstr, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)    :: attnam
 Integer,          Intent(IN)    :: ncid, varid, lenstr
 Character(LEN=*), Intent(INOUT) :: string 
 Integer,          Intent(OUT)   :: rcode

 Integer(KIND=C_INT)            :: cncid, cvarid, crcode
 Character(LEN=(LEN(attnam)+1)) :: cattnam
 Character(LEN=(lenstr+1))      :: cstring
 Integer                        :: ilen

 cncid   = ncid
 cvarid  = varid - 1
 rcode   = 0
 string  = REPEAT(" ", LEN(string))
 cstring = REPEAT(" ", LEN(cstring))

! check for a null character on end of attnam

 cattnam = addCNullChar(attnam, ilen)
 
 Call c_ncagtc(cncid, cvarid, cattnam(1:ilen+1), cstring, lenstr, &
               crcode)

 string(1:lenstr) = cstring(1:lenstr)

 rcode = crcode

 End Subroutine ncagtc
!-------------------------------- ncacpy --------------------------------------
 Subroutine ncacpy(ncid, varid, attnam, outcdf, outvar, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: attnam
 Integer,          Intent(IN)  :: ncid, varid, outcdf, outvar
 Integer,          Intent(OUT) :: rcode

 Integer(KIND=C_INT)            :: cncid, cvarid, coutcdf, coutvar, crcode
 Character(LEN=(LEN(attnam)+1)) :: cattnam
 Integer                        :: ilen

 cncid   = ncid
 cvarid  = varid - 1
 coutcdf = outcdf
 coutvar = outvar-1
 rcode   = 0

! check for a null character on end of attnam

 cattnam = addCNullChar(attnam, ilen)
 
 Call c_ncacpy(cncid, cvarid, cattnam(1:ilen+1), coutcdf, &
               coutvar, crcode)

 rcode = crcode

 End Subroutine ncacpy
!-------------------------------- ncanam --------------------------------------
 Subroutine ncanam(ncid, varid, attnum, attnam, rcode)

 USE netcdf_nc_interfaces, ONLY: NC_MAX_NAME
 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(INOUT) :: attnam
 Integer,          Intent(IN)    :: ncid, varid, attnum
 Integer,          Intent(OUT)   :: rcode

 Integer                      :: ilen
 Integer(KIND=C_INT)          :: cncid, cvarid, cattnum, crcode
 Character(LEN=NC_MAX_NAME+1) :: cattnam

 cncid = ncid
 cvarid = varid - 1
 cattnum = attnum - 1
 rcode = 0
 cattnam = REPEAT(" ", LEN(cattnam))
 ilen = LEN(attnam)

 Call c_ncanam(cncid, cvarid, cattnum, cattnam, crcode)

! check for a null character on end of cattnam

 attnam = stripCNullChar(cattnam, ilen)
 
 rcode = crcode

 End Subroutine ncanam
!-------------------------------- ncaren --------------------------------------
 Subroutine ncaren(ncid, varid, attnam, newnam, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: attnam, newnam
 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(OUT) :: rcode

 Integer(KIND=C_INT)            :: cncid, cvarid, crcode
 Character(LEN=(LEN(attnam)+1)) :: cattnam
 Character(LEN=(LEN(newnam)+1)) :: cnewnam
 Integer                        :: ilen, ilen2

 cncid  = ncid
 cvarid = varid - 1
 rcode  = 0

! check for a null character on end of attnam and newnam

 cattnam = addCNullChar(attnam, ilen)
 
 cnewnam = addCNullChar(newnam, ilen2)
 
 Call c_ncaren(cncid, cvarid, cattnam(1:ilen+1), cnewnam(1:ilen2+1), crcode) 

 rcode = crcode

 End Subroutine ncaren
!-------------------------------- ncadel --------------------------------------
 Subroutine ncadel(ncid, varid, attnam, rcode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: attnam
 Integer,          Intent(IN)  :: ncid, varid
 Integer,          Intent(OUT) :: rcode

 Integer(KIND=C_INT)            :: cncid, cvarid, crcode
 Character(LEN=(LEN(attnam)+1)) :: cattnam
 Integer                        :: ilen

 cncid  = ncid
 cvarid = varid - 1
 rcode  = 0

! check for a null character on end of attnam

 cattnam = addCNullChar(attnam, ilen)
 
 Call c_ncadel(cncid, cvarid, cattnam(1:ilen+1),crcode)

 rcode = crcode

 End Subroutine ncadel
!-------------------------------- ncsfil --------------------------------------
 Function ncsfil(ncid, fillmode, rcode) RESULT(currentmode)

 USE netcdf_fortv2_c_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid, fillmode
 Integer, Intent(OUT) :: rcode

 Integer              :: currentmode

 Integer(KIND=C_INT) :: cncid, cfillmode, crcode, cstatus  

 cncid     = ncid
 cfillmode = fillmode

 cstatus     = c_ncsfil(cncid, cfillmode, crcode)
 rcode       = crcode
 currentmode = cstatus

 End Function ncsfil
