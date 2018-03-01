! ------------ Routines to create/open/close/redefine netcdf files ------------ 

! Replacement for fort-control.c

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

! Version 1.: Sept.  2005 - Initial Cray X1 version
! Version 2.: May,   2006 - Updated to support g95
! Version 3.: April, 2009 - Updated for netcdf 4.0.1
! Version 4.: April, 2010 - Updated for netcdf 4.1.1
! Version 5.: Feb.   2013 - Added nf_inq_path support for fortran 4.4
! Vertion 6.: Nov.   2013 - Added nf_set_log_level support
! Version 7.: May,   2014 - Ensure return error status checked from C API calls
!
!-------------------------------- nf_create --------------------------------
 Function nf_create(path, cmode, ncid) RESULT (status)

! Creates a new NetCDF file given a file name and a creation mode and returns
! the file id and a status flag

 USE netcdf_nc_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: path
 Integer,          Intent(IN)  :: cmode
 Integer,          Intent(OUT) :: ncid
 
 Integer                       :: status

 Integer(KIND=C_INT)          :: ccmode, cncid, cstatus
 Character(LEN=(LEN(path)+1)) :: cpath
 Integer                      :: ie

 ccmode = cmode
 cncid  = 0
 
! Check for C null character on path. We will always add a null
! char so we don't need a second one

 cpath = addCNullChar(path, ie)
 
! Call nc_create to create file

 cstatus = nc_create(cpath(1:ie+1), ccmode, cncid)
 
 If (cstatus == NC_NOERR) Then
    ncid   = cncid 
 EndIf
 status = cstatus

 End Function nf_create
!-------------------------------- nf__create -------------------------------
 Function nf__create(path, cmode, initialsz, chunksizehintp, ncid) &
                        RESULT(status)

! Creates a new NetCDF file and returns the file id and a status flag
! This is an alternate form of nf_create that allows user to input
! two additional tuning parameters

 USE netcdf_nc_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: path
 Integer,          Intent(IN)  :: cmode, initialsz, chunksizehintp
 Integer,          Intent(OUT) :: ncid
 
 Integer                       :: status

 Integer(KIND=C_INT)          :: ccmode, cncid, cstatus
 Integer(KIND=C_SIZE_T)       :: cinit, cchunk
 Character(LEN=(LEN(path)+1)) :: cpath
 Integer                      :: ie

 ccmode = cmode
 cchunk = chunksizehintp
 cinit  = initialsz
 cncid  = 0
 
! Check for C null character on path. We will always add a null
! char so we don't need a second one

 cpath = addCNullChar(path, ie)
 
! Call nc_create to create file

 cstatus = nc__create(cpath(1:ie+1), ccmode, cinit, cchunk, cncid)
 
 If (cstatus == NC_NOERR) Then
    ncid   = cncid 
 EndIf
 status = cstatus

 End Function nf__create
!-------------------------------- nf__create_mp ------------------------------
 Function nf__create_mp(path, cmode, initialsz, basepe, chunksizehintp, ncid) &
                        RESULT(status)

! Creates a new NetCDF file and returns the file id and a status flag
! This is an alternate form of nf__create for shared memory MPP systems 
! two additional tuning parameters

 USE netcdf_nc_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)  :: path
 Integer,          Intent(IN)  :: cmode, initialsz, chunksizehintp, basepe
 Integer,          Intent(OUT) :: ncid
 
 Integer                       :: status

 Integer(KIND=C_INT)          :: ccmode, cncid, cstatus
 Integer(KIND=C_INT), TARGET  :: cbasepe
 Integer(KIND=C_SIZE_T)       :: cinit, cchunk
 Type(C_PTR)                  :: cbasepeptr
 Character(LEN=(LEN(path)+1)) :: cpath
 Integer                      :: ie

 ccmode     = cmode
 cchunk     = chunksizehintp
 cinit      = initialsz
 cncid      = 0
 cbasepe    = basepe
 cbasepeptr = C_LOC(cbasepe)

! Check for C null character on path. We will always add a null
! char so we don't need a second one

 cpath = addCNullChar(path, ie)
 
! Call nc_create_mp to create file for base pe

 cstatus = nc__create_mp(cpath(1:ie+1), ccmode, cinit, cbasepeptr, &
                         cchunk, cncid)
 
 If (cstatus == NC_NOERR) Then
    ncid   = cncid 
 EndIf
 status = cstatus

 End Function nf__create_mp
!-------------------------------- nf_open ----------------------------------
 Function nf_open(path, mode, ncid) RESULT (status)

! Open an existing NetCDF file and return file id and a status flag

 USE netcdf_nc_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)    :: path
 Integer,          Intent(IN)    :: mode
 Integer,          Intent(INOUT) :: ncid
 
 Integer                         :: status

 Integer(KIND=C_INT)          :: cmode, cncid, cstatus
 Character(LEN=(LEN(path)+1)) :: cpath
 Integer                      :: ie

 cmode = mode
 cncid = 0
 
! Check for C null character on path. We will always add a null
! char so we don't need a second one

 cpath = addCNullChar(path, ie) 
 
! Call nc_create to create file

 cstatus = nc_open(cpath(1:ie+1), cmode, cncid)
 
 If (cstatus == NC_NOERR) Then
    ncid   = cncid
 EndIf
 status = cstatus

 End Function nf_open
!-------------------------------- nf__open ---------------------------------
 Function nf__open(path, mode, chunksizehintp, ncid) RESULT (status)

! Open an existing NetCDF file and return file id and a status flag
! Alternate form of nf_open with extra tuning parameter

 USE netcdf_nc_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)    :: path
 Integer,          Intent(IN)    :: mode, chunksizehintp
 Integer,          Intent(INOUT) :: ncid
 
 Integer                         :: status

 Integer(KIND=C_INT)          :: cmode, cncid, cstatus
 Integer(KIND=C_SIZE_T)       :: cchunk
 Character(LEN=(LEN(path)+1)) :: cpath
 Integer                      :: inull, ie

 cmode  = mode
 cchunk = chunksizehintp
 cncid  = 0
 
! Check for C null character in path. A null character is always added
! before we pass path to C we don't need a second one

 cpath = addCNullChar(path,ie)
 
! Call nc_create to create file

 cstatus = nc__open(cpath(1:ie+1), cmode, cchunk, cncid)
 
 If (cstatus == NC_NOERR) Then
    ncid   = cncid
 EndIf
 status = cstatus

 End Function nf__open
!-------------------------------- nf__open_mp --------------------------------
 Function nf__open_mp(path, mode, basepe, chunksizehintp, ncid) RESULT (status)

! Open an existing NetCDF file and return file id and a status flag
! Alternate form of nf__open with parameter to designate basepe on
! shared memory MPP systems. 

 USE netcdf_nc_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN)    :: path
 Integer,          Intent(IN)    :: mode, chunksizehintp, basepe
 Integer,          Intent(INOUT) :: ncid
 
 Integer                         :: status

 Integer(KIND=C_INT)          :: cmode, cncid, cstatus
 Integer(KIND=C_INT), TARGET  :: cbasepe
 Integer(KIND=C_SIZE_T)       :: cchunk
 Type(C_PTR)                  :: cbasepeptr
 Character(LEN=(LEN(path)+1)) :: cpath
 Integer                      :: ie

 cmode      = mode
 cchunk     = chunksizehintp
 cncid      = 0
 cbasepe    = basepe
 cbasepeptr = C_LOC(cbasepe)
 
! Check for C null character in path. A null character is always added
! before we pass path to C we don't need a second one

 cpath = addCNullChar(path, ie) 
 
! Call nc_create to create file

 cstatus = nc__open_mp(cpath(1:ie+1), cmode, cbasepeptr, cchunk, &
                       cncid)
 
 If (cstatus == NC_NOERR) Then
    ncid   = cncid
 EndIf
 status = cstatus

 End Function nf__open_mp
!-------------------------------- nf_inq_path ------------------------------
 Function nf_inq_path(ncid, pathlen, path) RESULT(status)

! Inquire about file pathname and name length

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer,          Intent(IN)    :: ncid
 Integer,          Intent(INOUT) :: pathlen
 Character(LEN=*), Intent(INOUT) :: path

 Integer                         :: status

 Integer(C_INT)             :: cncid, cstatus
 Integer(C_SIZE_T)          :: cpathlen
 Character(LEN=LEN(path)+1) :: tmppath

 cncid   = ncid
 path    = REPEAT(" ", LEN(path))
 tmppath = REPEAT(" ", LEN(tmppath))

 cstatus = nc_inq_path(cncid, cpathlen, tmppath)

 If (cstatus == NC_NOERR) Then
    pathlen = cpathlen
    If (pathlen > LEN(path)) pathlen = LEN(path)
    path = stripCNullchar(tmppath, pathlen)
 EndIf
 status = cstatus

 End Function nf_inq_path
!-------------------------------- nf_set_fill ------------------------------
 Function nf_set_fill(ncid, fillmode, old_mode) RESULT(status)
 
! Sets fill mode for given netcdf file returns old mode if present

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid, fillmode
 Integer, Intent(OUT) :: old_mode

 Integer              :: status

 Integer(KIND=C_INT) :: cncid, cfill, coldmode, cstatus

 cncid    = ncid
 cfill    = fillmode
 coldmode = 0

 cstatus = nc_set_fill(cncid, cfill, coldmode)

 If (cstatus == NC_NOERR) Then
    old_mode = coldmode
 EndIf
 status   = cstatus

 End Function nf_set_fill
!-------------------------------- nf_set_default_format --------------------
 Function nf_set_default_format(newform, old_format) RESULT(status)
 
! Sets new default data format. Used to toggle between 64 bit offset and
! classic mode 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: newform 
 Integer, Intent(OUT) :: old_format

 Integer              :: status

 Integer(KIND=C_INT) :: cnew, cold, cstatus

 cnew = newform

 cstatus = nc_set_default_format(cnew,cold)

 If (cstatus == NC_NOERR) Then
    old_format = cold
 EndIf
 status     = cstatus

 End Function nf_set_default_format
!-------------------------------- nf_redef ---------------------------------
 Function nf_redef(ncid) RESULT(status)
 
! Re-Enter definition mode for NetCDF file id ncid 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN) :: ncid

 Integer             :: status

 Integer(KIND=C_INT) :: cncid, cstatus

 cncid = ncid

 cstatus = nc_redef(cncid)

 status = cstatus

 End Function nf_redef
!-------------------------------- nf_enddef --------------------------------
 Function nf_enddef(ncid) RESULT(status)
 
! Exit definition mode for NetCDF file id ncid

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN) :: ncid

 Integer             :: status

 Integer(KIND=C_INT) :: cncid, cstatus

 cncid = ncid

 cstatus = nc_enddef(cncid)

 status = cstatus

 End Function nf_enddef
!-------------------------------- nf__enddef -------------------------------
 Function nf__enddef(ncid, h_minfree, v_align, v_minfree, r_align) &
                     RESULT(status)
 
! Exit definition mode for NetCDF file id ncid. Alternate version
! with additional tuning parameters

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN) :: ncid, h_minfree, v_align, v_minfree, r_align

 Integer             :: status

 Integer(KIND=C_INT)    :: cncid, cstatus
 Integer(KIND=C_SIZE_T) :: chminfree, cvalign, cvminfree, cralign

 cncid     = ncid
 chminfree = h_minfree
 cvalign   = v_align
 cvminfree = v_minfree
 cralign   = r_align

 cstatus = nc__enddef(cncid, chminfree, cvalign, cvminfree, cralign)

 status = cstatus

 End Function nf__enddef
!-------------------------------- nf_sync ----------------------------------
 Function nf_sync(ncid) RESULT(status)
 
! synch up all open NetCDF files 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN) :: ncid

 Integer             :: status

 Integer(KIND=C_INT) :: cncid, cstatus

 cncid = ncid

 cstatus = nc_sync(cncid)

 status = cstatus

 End Function nf_sync
!-------------------------------- nf_abort ---------------------------------
 Function nf_abort(ncid) RESULT(status)
 
! Abort netCDF file creation and exit 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN) :: ncid

 Integer             :: status

 Integer(KIND=C_INT) :: cncid, cstatus

 cncid = ncid

 cstatus = nc_abort(cncid)

 status = cstatus

 End Function nf_abort
!-------------------------------- nf_close ---------------------------------
 Function nf_close(ncid) RESULT(status)
 
! Close netCDF file id ncid 

 USE netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN) :: ncid

 Integer             :: status

 Integer(KIND=C_INT) :: cncid, cstatus

 cncid   = ncid

 cstatus = nc_close(cncid)

 status  = cstatus

 End Function nf_close
!-------------------------------- nf_delete --------------------------------
 Function nf_delete(path) RESULT(status)
 
! Delete netCDF file id ncid 

 USE netcdf_nc_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN) :: path

 Integer                      :: status

 Integer(KIND=C_INT)          :: cstatus
 Character(LEN=(LEN(path)+1)) :: cpath
 Integer                      :: ie

 cpath = addCNullChar(path,ie)
 
 cstatus = nc_delete(cpath(1:ie+1))

 status = cstatus

 End Function nf_delete
!-------------------------------- nf_delete_mp -------------------------------
 Function nf_delete_mp(path, pe) RESULT(status)
 
! Delete netCDF file id ncid. Alternate form of nf_delete for shared memory
! MPP systems.

 USE netcdf_nc_interfaces

 Implicit NONE

 Character(LEN=*), Intent(IN) :: path
 Integer,          Intent(IN) :: pe

 Integer                      :: status

 Integer(KIND=C_INT)          :: cstatus, cpe
 Character(LEN=(LEN(path)+1)) :: cpath
 Integer                      :: ie

 cpe = pe

 cpath = addCNullChar(path,ie)
 
 cstatus = nc_delete_mp(cpath(1:ie+1), cpe)

 status = cstatus

 End Function nf_delete_mp
!-------------------------------- nf_set_base_pe ------------------------------
 Function nf_set_base_pe(ncid, pe) RESULT(status)

! Sets base pe number on shared memory MPP systems

 Use netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN) :: ncid, pe

 Integer             :: status

 Integer(KIND=C_INT) :: cncid, cpe, cstatus

 cncid = ncid
 cpe   = pe

 cstatus = nc_set_base_pe(cncid, cpe)

 status = cstatus

 End Function nf_set_base_pe
!-------------------------------- nf_inq_base_pe ------------------------------
 Function nf_inq_base_pe(ncid, pe) RESULT(status)

! Gets previously set base pe number on shared memory MPP systems

 Use netcdf_nc_interfaces

 Implicit NONE

 Integer, Intent(IN)  :: ncid
 Integer, Intent(OUT) :: pe

 Integer              :: status

 Integer(KIND=C_INT) :: cncid, cpe, cstatus

 cncid = ncid

 cstatus = nc_inq_base_pe(cncid, cpe)

 If (cstatus == NC_NOERR) Then
    pe     = cpe
 EndIf
 status = cstatus
End Function nf_inq_base_pe
