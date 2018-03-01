  !
  ! external netcdf data types:
  !
  integer, parameter, public :: &
    nf90_byte   = 1,            &
    nf90_int1   = nf90_byte,    &
    nf90_char   = 2,            &
    nf90_short  = 3,            &
    nf90_int2   = nf90_short,   &
    nf90_int    = 4,            &
    nf90_int4   = nf90_int,     &
    nf90_float  = 5,            &
    nf90_real   = nf90_float,   &
    nf90_real4  = nf90_float,   &
    nf90_double = 6,            &
    nf90_real8  = nf90_double
                        
  !
  ! default fill values:
  !
  character (len = 1),           parameter, public :: &
    nf90_fill_char  = achar(0)
  integer (kind =  OneByteInt),  parameter, public :: &
    nf90_fill_byte  = -127,                           &
    nf90_fill_int1  = nf90_fill_byte
  integer (kind =  TwoByteInt),  parameter, public :: &
    nf90_fill_short = -32767,                         &
    nf90_fill_int2  = nf90_fill_short
  integer (kind = FourByteInt),  parameter, public :: &
    nf90_fill_int   = -2147483647
  real   (kind =  FourByteReal), parameter, public :: &
    nf90_fill_float = 9.9692099683868690e+36,         &
    nf90_fill_real  = nf90_fill_float,                &
    nf90_fill_real4 = nf90_fill_float
  real   (kind = EightByteReal), parameter, public :: &
    nf90_fill_double = 9.9692099683868690e+36,        &
    nf90_fill_real8  = nf90_fill_double

  !
  ! mode flags for opening and creating a netcdf dataset:
  !
  integer, parameter, public :: &
    nf90_nowrite   = 0,         &
    nf90_write     = 1,         &
    nf90_clobber   = 0,         &
    nf90_noclobber = 4,         &
    nf90_fill      = 0,         &
    nf90_nofill    = 256,       &
    nf90_64bit_offset    = 512, &
    nf90_lock      = 1024,      &
    nf90_share     = 2048,      & 
    nf90_diskless  = 8,         &
    nf90_mmap      = 16
  
  integer, parameter, public ::  &
    nf90_sizehint_default = 0,   & 
    nf90_align_chunk      = -1 

  !
  ! size argument for defining an unlimited dimension:
  !
  integer, parameter, public :: nf90_unlimited = 0

  !
  ! global attribute id:
  !
  integer, parameter, public :: nf90_global = 0

  !
  ! implementation limits:
  !
  integer, parameter, public :: &
    nf90_max_dims     = 1024,    &
    nf90_max_attrs    = 8192,   &
    nf90_max_vars     = 8192,   &
    nf90_max_name     = 256,    &
    nf90_max_var_dims = 1024
  
  !
  ! error codes:
  !
  integer, parameter, public :: &
    nf90_noerr        = 0,      & ! No Error
    nf90_ebadid       = -33,    & ! Not a valid ID
    nf90_eexist       = -35,    & ! netcdf file exists && NF90_NOCLOBBER
    nf90_einval       = -36,    & ! Invalid Argument
    nf90_eperm        = -37,    & ! Write to read only
    nf90_enotindefine = -38,    & ! Operation not allowed in data mode
    nf90_eindefine    = -39,    & ! Operation not allowed in define mode
    nf90_einvalcoords = -40,    & ! Index exceeds dimension bound
    nf90_emaxdims     = -41,    & ! nf90_max_dims exceeded
    nf90_enameinuse   = -42,    & ! String match to name in use
    nf90_enotatt      = -43,    & ! Attribute not found
    nf90_emaxatts     = -44,    & ! nf90_max_attrs exceeded
    nf90_ebadtype     = -45,    & ! Not a netcdf data type
    nf90_ebaddim      = -46,    & ! Invalid dimension id or name
    nf90_eunlimpos    = -47,    & ! nf90_unlimited in the wrong index
    nf90_emaxvars     = -48,    & ! nf90_max_vars exceeded
    nf90_enotvar      = -49,    & ! The variable ID is invalid for the specified netCDF dataset.
    nf90_eglobal      = -50,    & ! Action prohibited on nf90_global varid
    nf90_enotnc       = -51,    & ! Not a netcdf file
    nf90_ests         = -52,    & ! In Fortran, string too short
    nf90_emaxname     = -53,    & ! nf90_max_name exceeded
    nf90_eunlimit     = -54,    & ! nf90_unlimited size already in use
    nf90_enorecvars   = -55,    & ! nc_rec op when there are no record vars
    nf90_echar        = -56,    & ! Attempt to convert between text & numbers
    nf90_eedge        = -57,    & ! Start+count exceeds dimension bound
    nf90_estride      = -58,    & ! Illegal stride
    nf90_ebadname     = -59,    & ! Attribute or variable name contains illegal characters
    nf90_erange       = -60,    & ! Math result not representable
    nf90_enomem       = -61,    & ! Memory allocation (malloc) failure
    nf90_evarsize     = -62,    & ! One or more variable sizes violate format constraints 
    nf90_edimsize     = -63,    & ! Invalid dimension size
    nf90_etrunc       = -64,    & ! File likely truncated or possibly corrupted
    nf90_eaxistype    = -65       ! Unknown axis type.
  !
  ! more error codes for DAP
  !
  integer, parameter, public :: &
    nf90_edap         = -66,    & ! Generic DAP error
    nf90_ecurl        = -67,    & ! Generic libcurl error
    nf90_eio          = -68,    & ! Generic IO error
    nf90_enodata      = -69,    & ! Attempt to access variable with no data
    nf90_edapsvc      = -70,    & ! DAP server error
    nf90_edas         = -71,    & ! Malformed or inaccessible DAS
    nf90_edds         = -72,    & ! Malformed or inaccessible DDS
    nf90_edatadds     = -73,    & ! Malformed or inaccessible DATADDS
    nf90_edapurl      = -74,    & ! Malformed DAP URL
    nf90_edapconstraint = -75,  & ! Malformed DAP Constrain
    nf90_etranslation = -76,    & ! Untranslatable construct
    nf904_first_error = -100
  !
  ! error codes for netCDF-4
  !
integer, parameter, public :: &
    nf90_ehdferr      = -101,    & ! Error at HDF5 layer. 
    nf90_ecantread    = -102,    & ! Can't read. 
    nf90_ecantwrite   = -103,    & ! Can't write. 
    nf90_ecantcreate  = -104,    & ! Can't create. 
    nf90_efilemeta    = -105,    & ! Problem with file metadata. 
    nf90_edimmeta     = -106,    & ! Problem with dimension metadata. 
    nf90_eattmeta     = -107,    & ! Problem with attribute metadata. 
    nf90_evarmeta     = -108,    & ! Problem with variable metadata. 
    nf90_enocompound  = -109,    & ! Not a compound type. 
    nf90_eattexists   = -110,    & ! Attribute already exists. 
    nf90_enotnc4      = -111,    & ! Attempting netcdf-4 operation on netcdf-3 file.   
    nf90_estrictnc3   = -112,    & ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
    nf90_enotnc3      = -113,    & ! Attempting netcdf-3 operation on netcdf-4 file.   
    nf90_enopar       = -114,    & ! Parallel operation on file opened for non-parallel access.   
    nf90_eparinit     = -115,    & ! Error initializing for parallel access.   
    nf90_ebadgrpid    = -116,    & ! Bad group ID.   
    nf90_ebadtypid    = -117,    & ! Bad type ID.   
    nf90_etypdefined  = -118,    & ! Type has already been defined and may not be edited. 
    nf90_ebadfield    = -119,    & ! Bad field ID.   
    nf90_ebadclass    = -120,    & ! Bad class.   
    nf90_emaptype     = -121,    & ! Mapped access for atomic types only.   
    nf90_elatefill    = -122,    & ! Attempt to define fill value when data already exists. 
    nf90_elatedef     = -123,    & ! Attempt to define var properties, like deflate, after enddef. 
    nf90_edimscale    = -124,    & ! Probem with HDF5 dimscales.
    nf90_enogrp       = -125,    & ! No group found.
    nf90_estorage     = -126,    & ! Can't specify both contiguous and chunking.
    nf90_ebadchunk    = -127,    & ! Bad chunksize.
    nf90_enotbuilt    = -128,    & ! Attempt to use feature that was not turned on when netCDF was built.
    nf90_ediskless    = -129,    & ! Error in using diskless  access. 
    nf90_ecantextend  = -130,    & ! Attempt to extend dataset during ind. I/O operation. 
    nf90_empi         = -131,    & ! MPI operation failed. 
    nf904_last_error  = -131

  !
  ! error handling modes:
  !
  integer, parameter, public :: &
    nf90_fatal   = 1,           &
    nf90_verbose = 2

  !
  ! format version numbers:
  !
  integer, parameter, public :: &
    nf90_format_classic = 1,    &
    nf90_format_64bit = 2,      &
    nf90_format_netcdf4 = 3,    &
    nf90_format_netcdf4_classic = 4

! extra data types:
integer, parameter, public :: &
     nf90_ubyte = 7, &
     nf90_ushort = 8, &
     nf90_uint = 9, &
     nf90_int64 = 10, &
     nf90_uint64 = 11, &
     nf90_string = 12, &
     nf90_vlen = 13, &
     nf90_opaque = 14, &
     nf90_enum = 15, &
     nf90_compound = 16

                        
! extra default fill values:
integer (kind =  TwoByteInt),  parameter, public :: &
     nf90_fill_ubyte  = 255,                        &
     nf90_fill_uint1  = nf90_fill_ubyte
integer (kind =  FourByteInt),  parameter, public :: &
     nf90_fill_ushort = 65535,                      &
     nf90_fill_uint2  = nf90_fill_ushort
integer (kind = EightByteInt),  parameter, public :: &
     nf90_fill_uint   = 4294967295_EightByteInt

! Extra file create mode flags.
integer, parameter, public :: &
     nf90_netcdf4 = 4096, &
     nf90_hdf5 = 4096, & ! deprecated
     nf90_classic_model = 256

! Flags for parallel access.
integer, parameter, public :: nf90_independent = 0, nf90_collective = 1

! Flags for parallel I/O.
integer, parameter, public :: nf90_mpiio = 8192, nf90_mpiposix = 16384, &
     nf90_pnetcdf = 32768
  
! Extra variable flags.
integer, parameter, public :: &
     nf90_chunk_seq = 0, &
     nf90_chunk_sub = 1, &
     nf90_chunk_sizes = 2, &
     nf90_endian_native = 0, &
     nf90_endian_little = 1, &
     nf90_endian_big = 2, &
     nf90_chunked = 0, &
     nf90_notcontiguous = 0, &
     nf90_contiguous = 1, &
     nf90_nochecksum = 0, &
     nf90_fletcher32 = 1, &
     nf90_noshuffle = 0, &
     nf90_shuffle = 1, &
     nf90_szip_ec_option_mask = 4, &
     nf90_szip_nn_option_mask = 32

! This is the position of NC_NETCDF4 in cmode, counting from the
! right, starting (uncharacteristically for fortran) at 0. It's needed
! for the BTEST function calls.
integer, parameter, private :: NETCDF4_BIT = 12
