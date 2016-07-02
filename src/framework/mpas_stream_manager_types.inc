    integer, public, parameter :: MPAS_STREAM_ERR_FATAL  = 1000, &
                                  MPAS_STREAM_ERR_WARN   = 1001, &
                                  MPAS_STREAM_ERR_SILENT = 1002

    integer, public, parameter :: MPAS_STREAM_MGR_NOERR            =  0, &
                                  MPAS_STREAM_MGR_ERR_CLOBBER_FILE = -1, &
                                  MPAS_STREAM_MGR_ERR_CLOBBER_REC  = -2, &
                                  MPAS_STREAM_MGR_ERROR            = -3

    integer, public, parameter :: MPAS_STREAM_INPUT        = 1, &
                                  MPAS_STREAM_OUTPUT       = 2, &
                                  MPAS_STREAM_INPUT_OUTPUT = 3, &
                                  MPAS_STREAM_NONE         = 4

    integer, public, parameter :: MPAS_STREAM_PROPERTY_ACTIVE        =  5, &
                                  MPAS_STREAM_PROPERTY_IMMUTABLE     =  6, &
                                  MPAS_STREAM_PROPERTY_FILENAME      =  7, &
                                  MPAS_STREAM_PROPERTY_REF_TIME      =  8, &
                                  MPAS_STREAM_PROPERTY_RECORD_INTV   =  9, &
                                  MPAS_STREAM_PROPERTY_PRECISION     = 10, &
                                  MPAS_STREAM_PROPERTY_FILENAME_INTV = 11, &
                                  MPAS_STREAM_PROPERTY_CLOBBER       = 12, &
                                  MPAS_STREAM_PROPERTY_IOTYPE        = 13

    integer, public, parameter :: MPAS_STREAM_CLOBBER_NEVER     = 100, &
                                  MPAS_STREAM_CLOBBER_APPEND    = 101, &
                                  MPAS_STREAM_CLOBBER_TRUNCATE  = 102, &
                                  MPAS_STREAM_CLOBBER_OVERWRITE = 103

    type MPAS_streamManager_type

        integer :: numStreams = 0
        integer :: errorLevel

        type (MPAS_IO_context_type), pointer :: ioContext
        type (MPAS_Clock_type), pointer :: streamClock
        type (MPAS_Pool_type), pointer :: allFields
        type (MPAS_Pool_type), pointer :: allPackages
        type (MPAS_Pool_type), pointer :: allStructs
        type (MPAS_Pool_type), pointer :: defaultAtts

        type (MPAS_stream_list_type), pointer :: currentStream => NULL()
        type (MPAS_stream_list_type), pointer :: streams
        type (MPAS_stream_list_type), pointer :: alarms_in
        type (MPAS_stream_list_type), pointer :: alarms_out

    end type MPAS_streamManager_type


