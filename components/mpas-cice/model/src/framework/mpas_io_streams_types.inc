   integer, parameter :: MPAS_STREAM_EXACT_TIME              = 100, &
                         MPAS_STREAM_NEAREST                 = 101, &
                         MPAS_STREAM_LATEST_BEFORE           = 102, &
                         MPAS_STREAM_EARLIEST_AFTER          = 103, &
                         MPAS_STREAM_LATEST_STRICTLY_BEFORE  = 104, &
                         MPAS_STREAM_EARLIEST_STRICTLY_AFTER = 105
   
   integer, parameter :: MPAS_STREAM_NOERR              =  0, &
                         MPAS_STREAM_NOT_INITIALIZED    = -1, &
                         MPAS_STREAM_FIELD_NOT_FOUND    = -2, &
                         MPAS_STREAM_CLOBBER_FILE       = -3, &
                         MPAS_STREAM_CLOBBER_RECORD     = -4, &
                         MPAS_IO_ERR                    = -5

   integer, parameter :: FIELD_0D_INT   =  1, &
                         FIELD_1D_INT   =  2, &
                         FIELD_2D_INT   =  3, &
                         FIELD_3D_INT   =  4, &
                         FIELD_0D_REAL  =  5, &
                         FIELD_1D_REAL  =  6, &
                         FIELD_2D_REAL  =  7, &
                         FIELD_3D_REAL  =  8, &
                         FIELD_4D_REAL  =  9, &
                         FIELD_5D_REAL  =  10, &
                         FIELD_0D_CHAR  =  11, &
                         FIELD_1D_CHAR  =  12

   type field_list_type
      integer :: field_type = -999
      logical :: isDecomposed
      integer :: totalDimSize     ! Total size of outer dimension across all blocks for decomposed fields
      logical, dimension(:), pointer :: isAvailable => null() ! Used for reading var-arrays where one or more 
                                                              !   constitutent arrays may not be present in the input file
      type (field0dInteger), pointer :: int0dField => null()
      type (field1dInteger), pointer :: int1dField => null()
      type (field2dInteger), pointer :: int2dField => null()
      type (field3dInteger), pointer :: int3dField => null()
      type (field0dReal), pointer :: real0dField => null()
      type (field1dReal), pointer :: real1dField => null()
      type (field2dReal), pointer :: real2dField => null()
      type (field3dReal), pointer :: real3dField => null()
      type (field4dReal), pointer :: real4dField => null()
      type (field5dReal), pointer :: real5dField => null()
      type (field0dChar), pointer :: char0dField => null()
      type (field1dChar), pointer :: char1dField => null()
      type (field_list_type), pointer :: next => null()
   end type field_list_type

   type MPAS_Stream_type
      logical :: isInitialized = .false.
      integer :: ioFormat
      integer :: ioDirection
      integer :: defaultPrecision = MPAS_IO_NATIVE_PRECISION
      logical :: clobberRecords = .false.
      ! When blockWrite is true, a stream only writes a single block.
      ! This option is really only useful for writing debugging streams out.
      logical :: blockWrite = .false.
      character(len=StrKIND) :: filename
      type (MPAS_IO_Handle_type) :: fileHandle
      type (att_list_type), pointer :: attList => null()
      type (field_list_type), pointer :: fieldList => null()
   end type MPAS_Stream_type

