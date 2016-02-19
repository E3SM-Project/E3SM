#ifdef USE_PIO2
   integer, parameter :: MPAS_IO_OFFSET_KIND = PIO_OFFSET_KIND
#else
   integer, parameter :: MPAS_IO_OFFSET_KIND = PIO_OFFSET
#endif

   ! File access modes
   integer, parameter :: MPAS_IO_READ  = 1, &
                         MPAS_IO_WRITE = 2

   ! I/O formats
   integer, parameter :: MPAS_IO_NETCDF   = 3, &
                         MPAS_IO_PNETCDF  = 4, &
                         MPAS_IO_NETCDF4  = 5, &
                         MPAS_IO_PNETCDF5 = 6

   ! Field and attribute types
   integer, parameter :: MPAS_IO_REAL     =  7,  &
                         MPAS_IO_DOUBLE   =  8,  &
                         MPAS_IO_INT      =  9,  &
                         MPAS_IO_LOGICAL  = 10,  &
                         MPAS_IO_CHAR     = 11

   ! Field precision
   integer, parameter :: MPAS_IO_SINGLE_PRECISION = 12, &
                         MPAS_IO_DOUBLE_PRECISION = 13, &
#ifdef SINGLE_PRECISION
                         MPAS_IO_NATIVE_PRECISION = MPAS_IO_SINGLE_PRECISION
#else
                         MPAS_IO_NATIVE_PRECISION = MPAS_IO_DOUBLE_PRECISION
#endif

   ! Unlimited / record dimension
   integer, parameter :: MPAS_IO_UNLIMITED_DIM = -123456

   ! Error codes
   integer, parameter :: MPAS_IO_NOERR              =  0, &
                         MPAS_IO_ERR_INVALID_MODE   = -1, &
                         MPAS_IO_ERR_INVALID_FORMAT = -2, &
                         MPAS_IO_ERR_LONG_FILENAME  = -3, &
                         MPAS_IO_ERR_UNINIT_HANDLE  = -4, &
                         MPAS_IO_ERR_PIO            = -5, &
                         MPAS_IO_ERR_DATA_MODE      = -6, &
                         MPAS_IO_ERR_NOWRITE        = -7, &
                         MPAS_IO_ERR_REDEF_DIM      = -8, &
                         MPAS_IO_ERR_REDEF_VAR      = -9, &
                         MPAS_IO_ERR_UNDEFINED_DIM  = -10, &
                         MPAS_IO_ERR_UNDEFINED_VAR  = -11, &
                         MPAS_IO_ERR_REDEF_ATT      = -12, &
                         MPAS_IO_ERR_WRONG_ATT_TYPE = -13, &
                         MPAS_IO_ERR_NO_DECOMP      = -14, &
                         MPAS_IO_ERR_TWO_UNLIMITED_DIMS = -15, &
                         MPAS_IO_ERR_WRONG_MODE         = -16, &
                         MPAS_IO_ERR_NO_UNLIMITED_DIM   = -17, &
                         MPAS_IO_ERR_UNIMPLEMENTED      = -18, &
                         MPAS_IO_ERR_WOULD_CLOBBER      = -19, &
                         MPAS_IO_ERR_NOEXIST_READ       = -20, &
                         MPAS_IO_ERR_MISSING_DIM        = -21

   type MPAS_IO_Handle_type
      logical :: initialized = .false.
      logical :: preexisting_file = .false.
      logical :: data_mode = .false.
      type (file_desc_t) :: pio_file
      character (len=StrKIND) :: filename
      integer :: iomode
      integer :: ioformat
      integer :: pio_unlimited_dimid
      integer :: preexisting_records = 0
      integer (kind=MPAS_IO_OFFSET_KIND) :: frame_number = 1
      type (dimlist_type), pointer :: dimlist_head => null()
      type (dimlist_type), pointer :: dimlist_tail => null()
      type (fieldlist_type), pointer :: fieldlist_head => null()
      type (fieldlist_type), pointer :: fieldlist_tail => null()
      type (attlist_type), pointer :: attlist_head => null()
      type (attlist_type), pointer :: attlist_tail => null()
      type (mpas_io_context_type), pointer :: ioContext => null()
   end type MPAS_IO_Handle_type

   type decomphandle_type
      integer :: field_type
      integer, dimension(:), pointer :: dims
      integer, dimension(:), pointer :: indices
      type (io_desc_t) :: pio_iodesc
   end type decomphandle_type

   type atthandle_type
      character (len=StrKIND) :: attName
      integer :: attType
      integer :: attValueInt
      integer :: precision
      integer, dimension(:), pointer :: attValueIntA => null()
      real (kind=RKIND) :: attValueReal
      real (kind=RKIND), dimension(:), pointer :: attValueRealA => null()
      character (len=StrKIND) :: attValueText
   end type atthandle_type

   type dimhandle_type
      character (len=StrKIND) :: dimname
      logical :: is_unlimited_dim = .false.
      integer :: dimsize
      integer :: dimid
   end type dimhandle_type

   type fieldhandle_type
      character (len=StrKIND) :: fieldname
      integer :: fieldid
      type (Var_desc_t) :: field_desc
      integer :: field_type
      logical :: has_unlimited_dim = .false.
      integer :: ndims
      integer :: precision
      type (dimhandle_type), pointer, dimension(:) :: dims
      type (attlist_type), pointer :: attlist_head => null()
      type (attlist_type), pointer :: attlist_tail => null()
      type (decomphandle_type), pointer :: decomp => null()
   end type fieldhandle_type

   type decomplist_type
      type (decomphandle_type), pointer :: decomphandle
      type (decomplist_type), pointer :: next => null()
   end type decomplist_type

   type attlist_type
      type (atthandle_type), pointer :: atthandle
      type (attlist_type), pointer :: next => null()
   end type attlist_type

   type dimlist_type
      type (dimhandle_type), pointer :: dimhandle
      type (dimlist_type), pointer :: next => null()
   end type dimlist_type

   type fieldlist_type
      type (fieldhandle_type), pointer :: fieldhandle
      type (fieldlist_type), pointer :: next => null()
   end type fieldlist_type

   type mpas_io_context_type
      type (decomplist_type), pointer :: decomp_list => null()
      type (iosystem_desc_t), pointer :: pio_iosystem => null()
      integer :: master_pio_iotype = -999
      type (dm_info), pointer :: dminfo => null()
   end type mpas_io_context_type
