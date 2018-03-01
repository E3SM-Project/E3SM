Module netcdf_nc_data

! Data types and flags for Fortran2003 interfaces to NetCDF C routines
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

! Version 1.:  Sept. 2005 - Initial Cray X1 version
! Version 2.:  May 2006   - Updated to support g95
! Version 3.:  June 2006  - Updated to include netCDF 4 functions
! Version 4.:  July 2007  - Modified to work with 3.6.2 build system
! Version 5.:  April 2009 - Updated to NetCDF 4.0.1
! Version 6.:  April 2010 - Updated to NetCDF 4.1.1
! Version 7.:  Feb.  2012 - Added support for F2008 Intrinsic kinds
! Version 8.:  Feb.  2013 - Updated for netcdf fortran-4.4
! Version 9.:  April 2014 - Changed C_PTRDIFF_T to UCAR definitions

 USE ISO_C_BINDING  ! All subsequent USE associations of netcdf_nc_data
                    ! will inherit ISO_C_BINDING data

! The following will allow us to use the Fortran 2008 default intrinsic
! kind variables contained in Fortran 2008 ISO_FORTRAN_ENV module when
! compilers support it. Actually most of the major compilers (and even
! the latest gfortran) support these now (Feb. 2012)
 
#ifdef HAVE_F2008
 USE ISO_FORTRAN_ENV, ONLY: REAL32, REAL64, INT8, INT16, INT32, INT64
#endif
 Implicit NONE

#include "nfconfig.inc"

#ifndef HAVE_F2008
 
! Create our own REAL32, REAL64, INT8, INT16, INT32, INT64 if we don't have F2008
! ISO_FORTRAN_ENV module

 Integer, Parameter, PRIVATE :: REAL32 = SELECTED_REAL_KIND(P=6,  R=37)   ! float 
 Integer, Parameter, PRIVATE :: REAL64 = SELECTED_REAL_KIND(P=13, R=307)  ! double
 Integer, Parameter, PRIVATE :: INT8   = SELECTED_INT_KIND( 2)
 Integer, Parameter, PRIVATE :: INT16  = SELECTED_INT_KIND( 4)
 Integer, Parameter, PRIVATE :: INT32  = SELECTED_INT_KIND( 9)            ! int
 Integer, Parameter, PRIVATE :: INT64  = SELECTED_INT_KIND(18)            ! long long
#endif

! Set KIND parameters to shorter names used in f03 interface routines etc.

 Integer, Parameter :: RK4 = REAL32
 Integer, Parameter :: RK8 = REAL64
 Integer, Parameter :: IK1 = INT8
 Integer, Parameter :: IK2 = INT16
 Integer, Parameter :: IK4 = INT32
 Integer, Parameter :: IK8 = INT64
 
! Define processor/compiler dependent parameters for ptrdiff_t, signed char,
! and short types. Note prtdiff_t was not defined in the FORTRAN 2003
! standard as an interoperable type in ISO_C_BINDING but was added as part of
! the recent TS29113 Technical Specification "Futher Interoperability with C" 
! passed in 2012. For now we will make our own using C_INT32_T or C_INT64_T
! but allow users to use the default definition for compilers that support 
! TS29113 (like gfortran 4.8). Default will be C_INTPTR_T 

#ifndef HAVE_TS29113_SUPPORT
#if (SIZEOF_PTRDIFF_T == 4)
 Integer, Parameter :: C_PTRDIFF_T = C_INT32_T
#elif (SIZEOF_PTRDIFF_T == 8)
 Integer, Parameter :: C_PTRDIFF_T = C_INT64_T
#else
 Integer, Parameter :: C_PTRDIFF_T = C_INTPTR_T
#endif
#endif

! Set KIND parameters for 1 and 2 byte integers if the system 
! supports them based on what is set by configure in nfconfig.inc.
! The routines that use these values will issue an EBADTYPE error
! and exit if C_SIGNED_CHAR and C_SHORT are not supported in
! ISO_C_BINDING

! Set the KINDs to default integers otherwise.

! INT1 KINDs

#ifdef NF_INT1_IS_C_SIGNED_CHAR
 Integer, Parameter :: CINT1 = C_SIGNED_CHAR
 Integer, Parameter :: NFINT1 = IK1
#elif NF_INT1_IS_C_SHORT
 Integer, Parameter :: CINT1 = C_SHORT
 Integer, Parameter :: NFINT1 = IK2
#elif NF_INT1_IS_C_INT
 Integer, Parameter :: CINT1 = C_INT
 Integer, Parameter :: NFINT1 = IK4
#elif NF_INT1_IS_C_LONG
 Integer, Parameter :: CINT1 = C_LONG
 Integer, Parameter :: NFINT1 = IK8
#else
 Integer, Parameter :: CINT1 = C_SIGNED_CHAR
 Integer, Parameter :: NFINT1 = IK1
#endif

! INT2 KINDs

#ifdef NF_INT2_IS_C_SHORT
 Integer, Parameter :: CINT2 = C_SHORT
 Integer, Parameter :: NFINT2 = IK2
#elif NF_INT2_IS_C_INT
 Integer, Parameter :: CINT2 = C_INT
 Integer, Parameter :: NFINT2 = IK4
#elif NF_INT2_IS_C_LONG
 Integer, Parameter :: CINT2 = C_LONG
 Integer, Parameter :: NFINT2 = IK8
#else
 Integer, Parameter :: CINT2 = C_SHORT
 Integer, Parameter :: NFINT2 = IK2
#endif

! Set Fortran default integer kind. This
! should take care of the case were default
! integer is a 64 bit int (ala prehistoric
! CRAYS) 

#ifdef NF_INT_IS_C_LONG
 Integer, Parameter :: CINT = C_LONG
 Integer, Parameter :: NFINT = IK8 
#else
 Integer, Parameter :: CINT = C_INT
 Integer, Parameter :: NFINT = IK4
#endif

! Set Fortran default real kind. This should
! take care of the case were the default real
! type is a 64 bit real (ala prehistoric CRAYs) 

#ifdef NF_REAL_IS_C_DOUBLE
  Integer, Parameter :: NFREAL = RK8
#else
  Integer, Parameter :: NFREAL = RK4
#endif
 
!--------- Define default C interface parameters from netcdf.h   ---------------
! This is not a complete impementation of the C header files but 
! defines NC_ values equivalent to the values in the netcdf.inc files
! excludeing the V2 values

!                     NETCDF3 data
!               
! Define enumerator nc_type data as integers

 Integer(KIND=C_INT), Parameter :: NC_NAT    = 0
 Integer(KIND=C_INT), Parameter :: NC_BYTE   = 1
 Integer(KIND=C_INT), Parameter :: NC_CHAR   = 2
 Integer(KIND=C_INT), Parameter :: NC_SHORT  = 3
 Integer(KIND=C_INT), Parameter :: NC_INT    = 4
 Integer(KIND=C_INT), Parameter :: NC_FLOAT  = 5
 Integer(KIND=C_INT), Parameter :: NC_DOUBLE = 6

! Default fill values

 Character(KIND=C_CHAR),      Parameter :: NC_FILL_CHAR   = C_NULL_CHAR 
 Integer(KIND=C_SIGNED_CHAR), Parameter :: NC_FILL_BYTE   = -127_C_SIGNED_CHAR
 Integer(KIND=C_SHORT),       Parameter :: NC_FILL_SHORT  = -32767_C_SHORT
 Integer(KIND=C_INT),         Parameter :: NC_FILL_INT    = -2147483647_C_INT
 Real(KIND=C_FLOAT),          Parameter :: NC_FILL_FLOAT  = 9.9692099683868690E+36
 Real(KIND=C_DOUBLE),         Parameter :: NC_FILL_DOUBLE = 9.9692099683868690D+36

! Mode flags for opening and creating datasets

 Integer(KIND=C_INT), Parameter :: NC_NOWRITE          = 0
 Integer(KIND=C_INT), Parameter :: NC_WRITE            = 1
 Integer(KIND=C_INT), Parameter :: NC_CLOBBER          = 0
 Integer(KIND=C_INT), Parameter :: NC_NOCLOBBER        = 4
 Integer(KIND=C_INT), Parameter :: NC_FILL             = 0
 Integer(KIND=C_INT), Parameter :: NC_NOFILL           = 256
 Integer(KIND=C_INT), Parameter :: NC_LOCK             = 1024
 Integer(KIND=C_INT), Parameter :: NC_SHARE            = 2048
 Integer(KIND=C_INT), Parameter :: NC_STRICT_NC3       = 8
 Integer(KIND=C_INT), Parameter :: NC_64BIT_OFFSET     = 512
 Integer(KIND=C_INT), Parameter :: NC_SIZEHINT_DEFAULT = 0
 Integer(KIND=C_INT), Parameter :: NC_ALIGN_CHUNK      = -1
 Integer(KIND=C_INT), Parameter :: NC_FORMAT_CLASSIC   = 1
 Integer(KIND=C_INT), Parameter :: NC_FORMAT_64BIT     = 2
 Integer(KIND=C_INT), Parameter :: NC_DISKLESS         = 8
 Integer(KIND=C_INT), Parameter :: NC_MMAP             = 16


! Unlimited dimension size argument and global attibute ID

 Integer(KIND=C_INT),  Parameter :: NC_UNLIMITED = 0
 Integer(KIND=C_INT),  Parameter :: NC_GLOBAL    = 0 

! Implementation limits (WARNING!  SHOULD BE THE SAME AS C INTERFACE)

 Integer(KIND=C_INT), Parameter :: NC_MAX_DIMS     = 1024 
 Integer(KIND=C_INT), Parameter :: NC_MAX_ATTRS    = 8192 
 Integer(KIND=C_INT), Parameter :: NC_MAX_VARS     = 8192 
 Integer(KIND=C_INT), Parameter :: NC_MAX_NAME     = 256 
 Integer(KIND=C_INT), Parameter :: NC_MAX_VAR_DIMS = NC_MAX_DIMS

! Error codes

 Integer(KIND=C_INT), Parameter :: NC_NOERR        =  0
 Integer(KIND=C_INT), Parameter :: NC2_ERR         = -1
 Integer(KIND=C_INT), Parameter :: NC_SYSERR       = -31
 Integer(KIND=C_INT), Parameter :: NC_EXDR         = -32
 Integer(KIND=C_INT), Parameter :: NC_EBADID       = -33
 Integer(KIND=C_INT), Parameter :: NC_EBFILE       = -34
 Integer(KIND=C_INT), Parameter :: NC_EEXIST       = -35
 Integer(KIND=C_INT), Parameter :: NC_EINVAL       = -36
 Integer(KIND=C_INT), Parameter :: NC_EPERM        = -37
 Integer(KIND=C_INT), Parameter :: NC_ENOTINDEFINE = -38
 Integer(KIND=C_INT), Parameter :: NC_EINDEFINE    = -39
 Integer(KIND=C_INT), Parameter :: NC_EINVALCOORDS = -40
 Integer(KIND=C_INT), Parameter :: NC_EMAXDIMS     = -41
 Integer(KIND=C_INT), Parameter :: NC_ENAMEINUSE   = -42
 Integer(KIND=C_INT), Parameter :: NC_ENOTATT      = -43
 Integer(KIND=C_INT), Parameter :: NC_EMAXATTS     = -44
 Integer(KIND=C_INT), Parameter :: NC_EBADTYPE     = -45
 Integer(KIND=C_INT), Parameter :: NC_EBADDIM      = -46
 Integer(KIND=C_INT), Parameter :: NC_EUNLIMPOS    = -47
 Integer(KIND=C_INT), Parameter :: NC_EMAXVARS     = -48
 Integer(KIND=C_INT), Parameter :: NC_ENOTVAR      = -49
 Integer(KIND=C_INT), Parameter :: NC_EGLOBAL      = -50
 Integer(KIND=C_INT), Parameter :: NC_ENOTNC       = -51
 Integer(KIND=C_INT), Parameter :: NC_ESTS         = -52
 Integer(KIND=C_INT), Parameter :: NC_EMAXNAME     = -53
 Integer(KIND=C_INT), Parameter :: NC_EUNLIMIT     = -54
 Integer(KIND=C_INT), Parameter :: NC_ENORECVARS   = -55
 Integer(KIND=C_INT), Parameter :: NC_ECHAR        = -56
 Integer(KIND=C_INT), Parameter :: NC_EEDGE        = -57
 Integer(KIND=C_INT), Parameter :: NC_ESTRIDE      = -58
 Integer(KIND=C_INT), Parameter :: NC_EBADNAME     = -59
 Integer(KIND=C_INT), Parameter :: NC_ERANGE       = -60
 Integer(KIND=C_INT), Parameter :: NC_ENOMEM       = -61
 Integer(KIND=C_INT), Parameter :: NC_EVARSIZE     = -62
 Integer(KIND=C_INT), Parameter :: NC_EDIMSIZE     = -63
 Integer(KIND=C_INT), Parameter :: NC_ETRUNC       = -64

! Error handling codes

 Integer(KIND=C_INT), Parameter :: NC_FATAL   = 1
 Integer(KIND=C_INT), Parameter :: NC_VERBOSE = 2

#ifdef USE_NETCDF4
!                          NETCDF4 data
 Integer(KIND=C_INT), Parameter :: NC_FORMAT_NETCDF4         = 3
 Integer(KIND=C_INT), Parameter :: NC_FORMAT_NETCDF4_CLASSIC = 4
 Integer(KIND=C_INT), Parameter :: NC_NETCDF4                = 4096
 Integer(KIND=C_INT), Parameter :: NC_CLASSIC_MODEL          = 256

! extra netcdf4 types
 Integer(KIND=C_INT), Parameter :: NC_LONG     = NC_INT
 Integer(KIND=C_INT), Parameter :: NC_UBYTE    = 7
 Integer(KIND=C_INT), Parameter :: NC_USHORT   = 8 
 Integer(KIND=C_INT), Parameter :: NC_UINT     = 9
 Integer(KIND=C_INT), Parameter :: NC_INT64    = 10 
 Integer(KIND=C_INT), Parameter :: NC_UINT64   = 11 
 Integer(KIND=C_INT), Parameter :: NC_STRING   = 12
 Integer(KIND=C_INT), Parameter :: NC_VLEN     = 13
 Integer(KIND=C_INT), Parameter :: NC_OPAQUE   = 14
 Integer(KIND=C_INT), Parameter :: NC_ENUM     = 15
 Integer(KIND=C_INT), Parameter :: NC_COMPOUND = 16

!extra netcd4 fill values
 Integer(KIND=C_INT),       Parameter :: NC_FILL_UBYTE  = 255
 Integer(KIND=C_INT),       Parameter :: NC_FILL_USHORT = 65535
 Integer(KIND=C_LONG_LONG), Parameter :: NC_FILL_UINT   = 4294967295_C_LONG_LONG
 Integer(KIND=C_LONG_LONG), Parameter :: NC_FILL_INT64  = -9223372036854775806_C_LONG_LONG

! extra netcdf4 variable flags 
 Integer(KIND=C_INT), Parameter :: NC_CHUNK_SEQ      = 0 
 Integer(KIND=C_INT), Parameter :: NC_CHUNK_SUB      = 1 
 Integer(KIND=C_INT), Parameter :: NC_CHUNK_SIZES    = 2 
 Integer(KIND=C_INT), Parameter :: NC_ENDIAN_NATIVE  = 0 
 Integer(KIND=C_INT), Parameter :: NC_ENDIAN_LITTLE  = 1 
 Integer(KIND=C_INT), Parameter :: NC_ENDIAN_BIG     = 2 
 Integer(KIND=C_INT), Parameter :: NC_CHUNKED        = 0
 Integer(KIND=C_INT), Parameter :: NC_NOTCONTIGUOUS  = 0
 Integer(KIND=C_INT), Parameter :: NC_CONTIGUOUS     = 1
 Integer(KIND=C_INT), Parameter :: NC_NOCHECKSUM     = 0
 Integer(KIND=C_INT), Parameter :: NC_FLETCHER32     = 1
 Integer(KIND=C_INT), Parameter :: NC_NOSHUFFLE      = 0
 Integer(KIND=C_INT), Parameter :: NC_SHUFFLE        = 1
 Integer(KIND=C_INT), Parameter :: NC_INDEPENDENT    = 0
 Integer(KIND=C_INT), Parameter :: NC_COLLECTIVE     = 1

! flags for parallel i/o

 Integer(KIND=C_INT), Parameter :: NC_MPIIO          = 8192
 Integer(KIND=C_INT), Parameter :: NC_MPIPOSIX       = 16384
 Integer(KIND=C_INT), Parameter :: NC_PNETCDF        = 32768

 Integer(KIND=C_INT), Parameter :: NC_SZIP_EC_OPTION_MASK = 4
 Integer(KIND=C_INT), Parameter :: NC_SZIP_NN_OPTION_MASK = 32

! extra netcdf4 error flags
 Integer(KIND=C_INT), Parameter :: NC_EHDFERR        = -101
 Integer(KIND=C_INT), Parameter :: NC_ECANTREAD      = -102
 Integer(KIND=C_INT), Parameter :: NC_ECANTWRITE     = -103
 Integer(KIND=C_INT), Parameter :: NC_ECANTCREATE    = -104
 Integer(KIND=C_INT), Parameter :: NC_EFILEMETA      = -105
 Integer(KIND=C_INT), Parameter :: NC_EDIMMETA       = -106
 Integer(KIND=C_INT), Parameter :: NC_EATTMETA       = -107
 Integer(KIND=C_INT), Parameter :: NC_EVARMETA       = -108
 Integer(KIND=C_INT), Parameter :: NC_ENOCOMPOUND    = -109
 Integer(KIND=C_INT), Parameter :: NC_EATTEXISTS     = -110
 Integer(KIND=C_INT), Parameter :: NC_ENOTNC4        = -111
 Integer(KIND=C_INT), Parameter :: NC_ESTRICTNC3     = -112
 Integer(KIND=C_INT), Parameter :: NC_ENOTNC3        = -113
 Integer(KIND=C_INT), Parameter :: NC_ENOPAR         = -114
 Integer(KIND=C_INT), Parameter :: NC_EPARINIT       = -115
 Integer(KIND=C_INT), Parameter :: NC_EBADGRPID      = -116
 Integer(KIND=C_INT), Parameter :: NC_EBADTYPID      = -117
 Integer(KIND=C_INT), Parameter :: NC_ETYPDEFINED    = -118
 Integer(KIND=C_INT), Parameter :: NC_EBADFIELD      = -119
 Integer(KIND=C_INT), Parameter :: NC_EBADCLASS      = -120
 Integer(KIND=C_INT), Parameter :: NC_EMAPTYPE       = -121
 Integer(KIND=C_INT), Parameter :: NC_ELATEFILL      = -122
 Integer(KIND=C_INT), Parameter :: NC_ELATEDEF       = -123
 Integer(KIND=C_INT), Parameter :: NC_EDIMSCALE      = -124
 Integer(KIND=C_INT), Parameter :: NC_ENOGRP         = -125
#endif
!------------------------------------------------------------------------------
End Module netcdf_nc_data
