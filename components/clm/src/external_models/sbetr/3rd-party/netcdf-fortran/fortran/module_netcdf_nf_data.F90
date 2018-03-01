Module netcdf_nf_data

! Module for Netcdf FORTRAN 2003 nf parameters. This includes all the
! error condition parameters, external data types, fill values etc.
! for netCDF2,3,4

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

! Version 1. Sept. 2005 - initial Cray X1 version for netcdf 3.6,0
! Version 2. Apr.  2005 - updated to be consistent with netcdf 3.6.1
! Version 3. Apr.  2009 - updated for netCDF 4.0.1
! Version 4. Apr.  2010 - updated for netCDF 4.1.1
! Version 5. Feb.  2013 - Updated for netCDF 4.4

! This module is provided as a replacement for parts of the netcdf2.inc,
! netcdf3.inc and netcdf4.inc include files. It does not include the 
! external statements for the nf_ functions. The latter are not needed
! if you use the interfaces in module_netcdf_nf_interfaces. If you
! want the externals, just use the include files.

 USE netcdf_nc_data
                                                         
 Implicit NONE

#include "nfconfig.inc"

! Define values found in netcdf2.inc, netcdf3.inc, and netcdf4.inc
! Also some additional values from the NC_DATA interfaces
!                
! Define enumerator nc_type data as integers

 Integer, Parameter ::  NF_NAT    = NC_NAT 
 Integer, Parameter ::  NF_BYTE   = NC_BYTE 
 Integer, Parameter ::  NF_INT1   = NF_BYTE
 Integer, Parameter ::  NF_CHAR   = NC_CHAR 
 Integer, Parameter ::  NF_SHORT  = NC_SHORT 
 Integer, Parameter ::  NF_INT2   = NF_SHORT
 Integer, Parameter ::  NF_INT    = NC_INT 
 Integer, Parameter ::  NF_FLOAT  = NC_FLOAT 
 Integer, Parameter ::  NF_REAL   = NF_FLOAT
 Integer, Parameter ::  NF_DOUBLE = NC_DOUBLE 

! Default fill values

 Integer,          Parameter :: NF_FILL_CHAR  = IACHAR(NC_FILL_CHAR) 
 Integer(IK1),     Parameter :: NF_FILL_BYTE  = NC_FILL_BYTE
 Integer(IK2),     Parameter :: NF_FILL_SHORT = NC_FILL_SHORT
 Integer,          Parameter :: NF_FILL_INT   = NC_FILL_INT

 Real(RK4),        Parameter :: NF_FILL_FLOAT  = NC_FILL_FLOAT
 Real(RK4),        Parameter :: NF_FILL_REAL   = NF_FILL_FLOAT
 Real(RK4),        Parameter :: NF_FILL_REAL4  = NF_FILL_FLOAT 
 Real(RK8),        Parameter :: NF_FILL_DOUBLE = NC_FILL_DOUBLE
 Real(RK8),        Parameter :: NF_FILL_REAL8  = NF_FILL_DOUBLE 

! Mode flags for opening and creating datasets

 Integer, Parameter :: NF_NOWRITE          = NC_NOWRITE 
 Integer, Parameter :: NF_WRITE            = NC_WRITE 
 Integer, Parameter :: NF_CLOBBER          = NC_CLOBBER 
 Integer, Parameter :: NF_NOCLOBBER        = NC_NOCLOBBER 
 Integer, Parameter :: NF_FILL             = NC_FILL 
 Integer, Parameter :: NF_NOFILL           = NC_NOFILL 
 Integer, Parameter :: NF_LOCK             = NC_LOCK 
 Integer, Parameter :: NF_SHARE            = NC_SHARE 
 Integer, Parameter :: NF_STRICT_NC3       = NC_STRICT_NC3 
 Integer, Parameter :: NF_64BIT_OFFSET     = NC_64BIT_OFFSET 
 Integer, Parameter :: NF_SIZEHINT_DEFAULT = NC_SIZEHINT_DEFAULT 
 Integer, Parameter :: NF_ALIGN_CHUNK      = NC_ALIGN_CHUNK
 Integer, Parameter :: NF_FORMAT_CLASSIC   = NC_FORMAT_CLASSIC
 Integer, Parameter :: NF_FORMAT_64BIT     = NC_FORMAT_64BIT 
 Integer, Parameter :: NF_DISKLESS         = NC_DISKLESS 
 Integer, Parameter :: NF_MMAP             = NC_MMAP

! Unlimited dimension size argument and global attibute ID

 Integer,  Parameter :: NF_UNLIMITED = NC_UNLIMITED 
 Integer,  Parameter :: NF_GLOBAL    = NC_GLOBAL 

! Implementation limits (WARNING!  SHOULD BE THE SAME AS C INTERFACE)

 Integer, Parameter :: NF_MAX_DIMS     = NC_MAX_DIMS 
 Integer, Parameter :: NF_MAX_ATTRS    = NC_MAX_ATTRS 
 Integer, Parameter :: NF_MAX_VARS     = NC_MAX_VARS 
 Integer, Parameter :: NF_MAX_NAME     = NC_MAX_NAME 
 Integer, Parameter :: NF_MAX_VAR_DIMS = NC_MAX_DIMS

! Error codes

 Integer, Parameter :: NF_NOERR        = NC_NOERR 
 Integer, Parameter :: NF2_ERR         = NC2_ERR 
 Integer, Parameter :: NF_SYSERR       = NC_SYSERR 
 Integer, Parameter :: NF_EXDR         = NC_EXDR 
 Integer, Parameter :: NF_EBADID       = NC_EBADID 
 Integer, Parameter :: NF_EBFILE       = NC_EBFILE 
 Integer, Parameter :: NF_EEXIST       = NC_EEXIST 
 Integer, Parameter :: NF_EINVAL       = NC_EINVAL 
 Integer, Parameter :: NF_EPERM        = NC_EPERM 
 Integer, Parameter :: NF_ENOTINDEFINE = NC_ENOTINDEFINE
 Integer, Parameter :: NF_EINDEFINE    = NC_EINDEFINE
 Integer, Parameter :: NF_EINVALCOORDS = NC_EINVALCOORDS
 Integer, Parameter :: NF_EMAXDIMS     = NC_EMAXDIMS
 Integer, Parameter :: NF_ENAMEINUSE   = NC_ENAMEINUSE
 Integer, Parameter :: NF_ENOTATT      = NC_ENOTATT
 Integer, Parameter :: NF_EMAXATTS     = NC_EMAXATTS
 Integer, Parameter :: NF_EBADTYPE     = NC_EBADTYPE
 Integer, Parameter :: NF_EBADDIM      = NC_EBADDIM
 Integer, Parameter :: NF_EUNLIMPOS    = NC_EUNLIMPOS
 Integer, Parameter :: NF_EMAXVARS     = NC_EMAXVARS
 Integer, Parameter :: NF_ENOTVAR      = NC_ENOTVAR 
 Integer, Parameter :: NF_EGLOBAL      = NC_EGLOBAL 
 Integer, Parameter :: NF_ENOTNC       = NC_ENOTNC
 Integer, Parameter :: NF_ESTS         = NC_ESTS
 Integer, Parameter :: NF_EMAXNAME     = NC_EMAXNAME
 Integer, Parameter :: NF_EUNLIMIT     = NC_EUNLIMIT
 Integer, Parameter :: NF_ENORECVARS   = NC_ENORECVARS
 
 Integer, Parameter :: NF_ECHAR        = NC_ECHAR 
 Integer, Parameter :: NF_EEDGE        = NC_EEDGE
 Integer, Parameter :: NF_ESTRIDE      = NC_ESTRIDE 
 Integer, Parameter :: NF_EBADNAME     = NC_EBADNAME
 Integer, Parameter :: NF_ERANGE       = NC_ERANGE
 Integer, Parameter :: NF_ENOMEM       = NC_ENOMEM
 Integer, Parameter :: NF_EVARSIZE     = NC_EVARSIZE
 Integer, Parameter :: NF_EDIMSIZE     = NC_EDIMSIZE 
 Integer, Parameter :: NF_ETRUNC       = NC_ETRUNC 

! Error handling codes

 Integer, Parameter :: NF_FATAL   = NC_FATAL 
 Integer, Parameter :: NF_VERBOSE = NC_VERBOSE 

#ifdef USE_NETCDF4

! NETCDF4 parameters 

! data types

 Integer, Parameter :: NF_LONG     = NC_LONG
 Integer, Parameter :: NF_UBYTE    = NC_UBYTE 
 Integer, Parameter :: NF_USHORT   = NC_USHORT 
 Integer, Parameter :: NF_UINT     = NC_UINT 
 Integer, Parameter :: NF_INT64    = NC_INT64 
 Integer, Parameter :: NF_UINT64   = NC_UINT64 
 Integer, Parameter :: NF_STRING   = NC_STRING 
 Integer, Parameter :: NF_VLEN     = NC_VLEN 
 Integer, Parameter :: NF_OPAQUE   = NC_OPAQUE 
 Integer, Parameter :: NF_ENUM     = NC_ENUM 
 Integer, Parameter :: NF_COMPOUND = NC_COMPOUND

! Netcdf4 fill flags - for some reason the F90 values are different

 Integer,      Parameter :: NF_FILL_UBYTE  = NC_FILL_UBYTE 
 Integer,      Parameter :: NF_FILL_UINT1  = NF_FILL_UBYTE 
 Integer,      Parameter :: NF_FILL_USHORT = NC_FILL_USHORT
 Integer,      Parameter :: NF_FILL_UINT2  = NF_FILL_USHORT
 Integer(IK8), Parameter :: NF_FILL_UINT   = NC_FILL_UINT
 Integer(IK8), Parameter :: NF_FILL_INT64  = NC_FILL_INT64 

! new format types
 Integer, Parameter :: NF_FORMAT_NETCDF4         = NC_FORMAT_NETCDF4 
 Integer, Parameter :: NF_FORMAT_NETCDF4_CLASSIC = NC_FORMAT_NETCDF4_CLASSIC
 
! Netcdf4 create mode flags
 Integer, Parameter :: NF_NETCDF4        = NC_NETCDF4 
 Integer, Parameter :: NF_HDF5           = NF_NETCDF4 ! deprecated
 Integer, Parameter :: NF_CLASSIC_MODEL  = NC_CLASSIC_MODEL 
! Netcdf4 variable flags
 Integer, Parameter :: NF_CHUNK_SEQ      = NC_CHUNK_SEQ 
 Integer, Parameter :: NF_CHUNK_SUB      = NC_CHUNK_SUB 
 Integer, Parameter :: NF_CHUNK_SIZES    = NC_CHUNK_SIZES 
 Integer, Parameter :: NF_ENDIAN_NATIVE  = NC_ENDIAN_NATIVE 
 Integer, Parameter :: NF_ENDIAN_LITTLE  = NC_ENDIAN_LITTLE 
 Integer, Parameter :: NF_ENDIAN_BIG     = NC_ENDIAN_BIG 
 Integer, Parameter :: NF_CHUNKED        = NC_CHUNKED 
 Integer, Parameter :: NF_NOTCONTIGUOUS  = NC_NOTCONTIGUOUS 
 Integer, Parameter :: NF_CONTIGUOUS     = NC_CONTIGUOUS 
 Integer, Parameter :: NF_NOCHECKSUM     = NC_NOCHECKSUM 
 Integer, Parameter :: NF_FLETCHER32     = NC_FLETCHER32 
 Integer, Parameter :: NF_NOSHUFFLE      = NC_NOSHUFFLE
 Integer, Parameter :: NF_SHUFFLE        = NC_SHUFFLE 
 Integer, Parameter :: NF_INDEPENDENT    = NC_INDEPENDENT 
 Integer, Parameter :: NF_COLLECTIVE     = NC_COLLECTIVE

! Flags for parallel I/O

 Integer, Parameter :: NF_MPIIO          = NC_MPIIO 
 Integer, Parameter :: NF_MPIPOSIX       = NC_MPIPOSIX 
 Integer, Parameter :: NF_PNETCDF        = NC_PNETCDF

! SZIP flags
 
 Integer, Parameter :: NF_SZIP_EC_OPTION_MASK = NC_SZIP_EC_OPTION_MASK 
 Integer, Parameter :: NF_SZIP_NN_OPTION_MASK = NC_SZIP_NN_OPTION_MASK 

! Netcdf4 error flags

 Integer, Parameter :: NF_EHDFERR        = NC_EHDFERR 
 Integer, Parameter :: NF_ECANTREAD      = NC_ECANTREAD 
 Integer, Parameter :: NF_ECANTWRITE     = NC_ECANTWRITE 
 Integer, Parameter :: NF_ECANTCREATE    = NC_ECANTCREATE 
 Integer, Parameter :: NF_EFILEMETA      = NC_EFILEMETA 
 Integer, Parameter :: NF_EDIMMETA       = NC_EDIMMETA 
 Integer, Parameter :: NF_EATTMETA       = NC_EATTMETA
 Integer, Parameter :: NF_EVARMETA       = NC_EVARMETA 
 Integer, Parameter :: NF_ENOCOMPOUND    = NC_ENOCOMPOUND 
 Integer, Parameter :: NF_EATTEXISTS     = NC_EATTEXISTS 
 Integer, Parameter :: NF_ENOTNC4        = NC_ENOTNC4 
 Integer, Parameter :: NF_ESTRICTNC3     = NC_ESTRICTNC3 
 Integer, Parameter :: NF_ENOTNC3        = NC_ENOTNC3 
 Integer, Parameter :: NF_ENOPAR         = NC_ENOPAR
 Integer, Parameter :: NF_EPARINIT       = NC_EPARINIT 
 Integer, Parameter :: NF_EBADGRPID      = NC_EBADGRPID
 Integer, Parameter :: NF_EBADTYPID      = NC_EBADTYPID
 Integer, Parameter :: NF_ETYPDEFINED    = NC_ETYPDEFINED
 Integer, Parameter :: NF_EBADFIELD      = NC_EBADFIELD
 Integer, Parameter :: NF_EBADCLASS      = NC_EBADCLASS 
 Integer, Parameter :: NF_EMAPTYPE       = NC_EMAPTYPE 
 Integer, Parameter :: NF_ELATEFILL      = NC_ELATEFILL 
 Integer, Parameter :: NF_ELATEDEF       = NC_ELATEDEF 
 Integer, Parameter :: NF_EDIMSCALE      = NC_EDIMSCALE 
 Integer, Parameter :: NF_ENOGRP         = NC_ENOGRP 
#endif

#ifndef NO_NETCDF_2
! V2 interface values

 Integer, Parameter :: NCBYTE     = NF_BYTE 
 Integer, Parameter :: NCCHAR     = NF_CHAR 
 Integer, Parameter :: NCSHORT    = NF_SHORT 
 Integer, Parameter :: NCLONG     = NF_INT 
 Integer, Parameter :: NCFLOAT    = NF_FLOAT 
 Integer, Parameter :: NCDOUBLE   = NF_DOUBLE 

 Integer, Parameter :: NCRDWR     = NF_WRITE 
 Integer, Parameter :: NCCREATE   = 2
 Integer, Parameter :: NCEXCL     = 4
 Integer, Parameter :: NCINDEF    = 8
 Integer, Parameter :: NCNSYNC    = 16
 Integer, Parameter :: NCHSYNC    = 32
 Integer, Parameter :: NCNDIRTY   = 64
 Integer, Parameter :: NCHDIRTY   = 128
 Integer, Parameter :: NCFILL     = NF_FILL 
 Integer, Parameter :: NCNOFILL   = NF_NOFILL 
 Integer, Parameter :: NCLINK     = 32768

 Integer, Parameter :: NCNOWRIT   = NF_NOWRITE 
 Integer, Parameter :: NCWRITE    = NF_WRITE
 Integer, Parameter :: NCCLOB     = NF_CLOBBER
 Integer, Parameter :: NCNOCLOB   = NF_NOCLOBBER

 Integer, Parameter :: NCUNLIM    = NF_UNLIMITED 
 Integer, Parameter :: NCGLOBAL   = NF_GLOBAL 

 Integer, Parameter :: MAXNCOP    = 64
 Integer, Parameter :: MAXNCDIM   = NF_MAX_DIMS 
 Integer, Parameter :: MAXNCATT   = NF_MAX_ATTRS 
 Integer, Parameter :: MAXNCVAR   = NF_MAX_VARS 
 Integer, Parameter :: MAXNCNAM   = NF_MAX_NAME 
 Integer, Parameter :: MAXVDIMS   = MAXNCDIM

 Integer, Parameter :: NCNOERR    = NF_NOERR
 Integer, Parameter :: NCEBADID   = NF_EBADID
 Integer, Parameter :: NCENFILE   = -31
 Integer, Parameter :: NCEEXIST   = NF_EEXIST
 Integer, Parameter :: NCEINVAL   = NF_EINVAL
 Integer, Parameter :: NCEPERM    = NF_EPERM
 Integer, Parameter :: NCENOTIN   = NF_ENOTINDEFINE
 Integer, Parameter :: NCEINDEF   = NF_EINDEFINE
 Integer, Parameter :: NCECOORD   = NF_EINVALCOORDS
 Integer, Parameter :: NCEMAXDS   = NF_EMAXDIMS
 Integer, Parameter :: NCENAME    = NF_ENAMEINUSE
 Integer, Parameter :: NCEMAXAT   = NF_EMAXATTS
 Integer, Parameter :: NCEBADTY   = NF_EBADTYPE
 Integer, Parameter :: NCEBADD    = NF_EBADDIM
 Integer, Parameter :: NCEUNLIM   = NF_EUNLIMPOS
 Integer, Parameter :: NCEMAXVS   = NF_EMAXVARS
 Integer, Parameter :: NCENOTVR   = NF_ENOTVAR
 Integer, Parameter :: NCEGLOB    = NF_EGLOBAL
 Integer, Parameter :: NCNOTNC    = NF_ENOTNC
 Integer, Parameter :: NCESTC     = NF_ESTS
 Integer, Parameter :: NCENTOOL   = NF_EMAXNAME
 Integer, Parameter :: NCFOOBAR   = 32
 Integer, Parameter :: NCSYSERR   = NF_SYSERR 

 Integer, Parameter :: NCFATAL    = NF_FATAL 
 Integer, Parameter :: NCVERBOS   = NF_VERBOSE 

 Integer,      Parameter :: FILCHAR  = NF_FILL_CHAR 
 Integer(IK1), Parameter :: FILBYTE  = NF_FILL_BYTE 
 Integer(IK2), Parameter :: FILSHORT = NF_FILL_SHORT 
 Integer,      Parameter :: FILLONG  = NF_FILL_INT 
 Real(RK4),    Parameter :: FILFLOAT = NF_FILL_FLOAT 
 Real(RK8),    Parameter :: FILDOUB  = NF_FILL_DOUBLE 
#endif

!------------------------------------------------------------------------------
End Module netcdf_nf_data
