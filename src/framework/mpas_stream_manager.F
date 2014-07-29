module mpas_stream_manager

#define STREAM_DEBUG_WRITE(M) write(stderrUnit,*) (M)

    use mpas_kind_types
    use mpas_grid_types
    use mpas_timekeeping
    use mpas_io_units
    use mpas_io_streams

    integer, public, parameter :: MPAS_STREAM_ERR_FATAL  = 1000, &
                                  MPAS_STREAM_ERR_WARN   = 1001, &
                                  MPAS_STREAM_ERR_SILENT = 1002

    integer, public, parameter :: MPAS_STREAM_MGR_NOERR = 0

    integer, public, parameter :: MPAS_STREAM_INPUT = 1, &
                                  MPAS_STREAM_OUTPUT = 2, &
                                  MPAS_STREAM_INPUT_OUTPUT = 3, &
                                  MPAS_STREAM_NONE = 4

    integer, public, parameter :: MPAS_STREAM_IMMUTABLE = 5, &
                                  MPAS_STREAM_MUTABLE   = 6

    character(len=*), public, parameter :: MPAS_STREAM_PROPERTY_DIRECTION = 'direction', &
                                           MPAS_STREAM_PROPERTY_IMMUTABLE = 'immutable', &
                                           MPAS_STREAM_PROPERTY_FILENAME = 'filename', &
                                           MPAS_STREAM_PROPERTY_MAX_RECORDS = 'maxRecords', &
                                           MPAS_STREAM_PROPERTY_REF_TIME = 'referenceTime'

    public :: MPAS_streamManager_type, &
              MPAS_stream_mgr_init, &
              MPAS_stream_mgr_finalize, &
              MPAS_stream_mgr_create_stream, &
              MPAS_stream_mgr_destroy_stream, &
              MPAS_stream_mgr_set_property, &
              MPAS_stream_mgr_get_property, &
              MPAS_stream_mgr_add_pkg, &
              MPAS_stream_mgr_remove_pkg, &
              MPAS_stream_mgr_add_field, &
              MPAS_stream_mgr_remove_field, &
              MPAS_stream_mgr_add_alarm, &
              MPAS_stream_mgr_remove_alarm, &
              MPAS_stream_mgr_write, &
              MPAS_stream_mgr_read, &
              MPAS_stream_mgr_expand_filename

    private

    type MPAS_streamManager_type
        integer :: numStreams = 0
        integer :: errorLevel
        type (MPAS_Pool_type), pointer :: streams
        type (MPAS_Pool_type), pointer :: alarms
        type (MPAS_Clock_type), pointer :: streamClock
        type (MPAS_Pool_type), pointer :: allFields
        type (MPAS_Pool_type), pointer :: allPackages
    end type MPAS_streamManager_type

    interface MPAS_stream_mgr_add_att
        module procedure MPAS_stream_mgr_add_att_int
        module procedure MPAS_stream_mgr_add_att_real
        module procedure MPAS_stream_mgr_add_att_char
        module procedure MPAS_stream_mgr_add_att_logical
    end interface


    contains


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_init
    !
    !> \brief Initialize a new MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Instantiates and initializes a streamManager type with a timekeeping
    !>  clock and a pool from which fields may be drawn and added to streams.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_init(manager, clock, allFields, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_init'

        type (MPAS_streamManager_type), pointer :: manager
        type (MPAS_Clock_type), pointer :: clock
        type (MPAS_Pool_type), pointer :: allFields
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_init()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        allocate(manager)
        manager % allFields => allFields
        manager % streamClock => clock
        manager % numStreams = 0
        manager % errorLevel = MPAS_STREAM_ERR_SILENT
        call mpas_pool_create_pool(manager % streams)
        call mpas_pool_create_pool(manager % alarms)


    end subroutine MPAS_stream_mgr_init!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_finalize
    !
    !> \brief Free all memory associated with an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Destroys a streamManager type, freeing all memory that was created as
    !>  part of the manager; the external clock and field pool associated with
    !>  the streamManager are unaffected.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_finalize(manager, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_finalize'

        type (MPAS_streamManager_type), pointer:: manager
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_finalize()') 

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_empty_pool(manager % streams)
        deallocate(manager % streams)
        call mpas_pool_empty_pool(manager % alarms)
        deallocate(manager % alarms)
        deallocate(manager)


    end subroutine MPAS_stream_mgr_finalize!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_create_stream
    !
    !> \brief Instantiate a new stream within an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Creates a new stream within the stream manager. The "direction" 
    !>  argument may be either MPAS_STREAM_INPUT, MPAS_STREAM_OUTPUT, 
    !>  MPAS_STREAM_INPUT_OUTPUT, or MPAS_STREAM_NONE. The "filename" argument 
    !>  is the template of the filenames that are associated with the stream. 
    !>  The "maxRecords" argument specifies the maximum number of records that are 
    !>  stored in each file. Knowing how many records are in each file, and 
    !>  the filename template, a "referenceTime" argument must be provided to 
    !>  specify the first timestamp appearing in any of the files associated with 
    !>  the stream, thereby determining where the "file breaks" will occur between 
    !>  timestamps. If no "referenceTime" is specified, the start time of the 
    !>  clock associated with the stream handler will be used as the reference 
    !>  time.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_create_stream(manager, streamID, direction, filename, maxRecords, referenceTime, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_create_stream'

        type (MPAS_streamManager_type), pointer :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: direction
        character (len=*), intent(in) :: filename
        integer, intent(in) :: maxRecords
        type (MPAS_Time_type), intent(in), optional :: referenceTime
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: attributesStream
        type (MPAS_Pool_type), pointer :: newStream
        type (MPAS_Pool_type), pointer :: attachedPackages
        type (MPAS_Pool_type), pointer :: alarms

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_create_stream()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_create_pool(attributesStream)
        call mpas_pool_create_pool(newStream)
        call mpas_pool_create_pool(attachedPackages)
        call mpas_pool_create_pool(alarms)

        ! Add properties
        call mpas_pool_add_config(newStream, MPAS_STREAM_PROPERTY_DIRECTION, MPAS_STREAM_NONE)
        call mpas_pool_add_config(newStream, MPAS_STREAM_PROPERTY_IMMUTABLE, MPAS_STREAM_MUTABLE)

        call mpas_pool_add_config(newStream, MPAS_STREAM_PROPERTY_FILENAME, filename)
        call mpas_pool_add_config(newStream, MPAS_STREAM_PROPERTY_MAX_RECORDS, maxRecords)

        !DWJ 07/16/2014: Need to deal with adding a time type into the pool
!       if (present(referenceTime)) then
!           call mpas_pool_add_config(newStream, MPAS_STERAM_PROPERTY_REF_TIME, referenceTime)
!       end if
        call mpas_pool_add_subpool(newStream, 'globalAttributes', attributesStream)
        call mpas_pool_add_subpool(newStream, 'attachedPackages', attachedPackages)
        call mpas_pool_add_subpool(newStream, 'alarms', alarms)
        call mpas_pool_add_subpool(manager % streams, streamID, newStream)

        manager % numStreams = manager % numStreams + 1
        write(0,*) 'DEBUG TEST we have this many streams ', manager % numStreams

    end subroutine MPAS_stream_mgr_create_stream!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_destroy_stream
    !
    !> \brief Free all memory associated with a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Destroy the stream, including freeing all memory explicitly associated with the stream.
    !>  This will not deallocate the memory associated with the fields in the stream.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_destroy_stream(manager, streamID, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_destroy_stream'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: stream
        character (len=StrKIND) :: error

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_destroy_stream()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_get_subpool(manager % streams, streamID, stream)

        if (.not. associated(stream)) then
            write(error,*) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not exist in manager.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_remove_subpool(manager % streams, streamID)

        call mpas_pool_empty_pool(stream)
        call mpas_pool_destroy_pool(stream)

        ! DWJ 07/22/2014: Might be leaking memory associated with subpools (alarms, globalAttributes, etc)

    end subroutine MPAS_stream_mgr_destroy_stream!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_set_property
    !
    !> \brief Sets a property of a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Sets the value of a stream property within an MPAS stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_set_property(manager, streamID, propertyName, propertyValue, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_set_property'

        type (MPAS_streamManager_type), pointer :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: propertyName
        integer, intent(in) :: propertyValue
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: stream
        integer, pointer :: curProperty

        character (len=StrKIND) :: error

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_set_property()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_get_subpool(manager % streams, streamID, stream)

        if (.not. associated(stream)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not exist'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_get_config(stream, propertyName, curProperty)

        if (.not.associated(curProperty)) then
            write(error, *) 'ERROR: ', trim(sub), ' property ', trim(propertyName), ' does not exist in stream ', trim(streamID) 
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        curProperty = propertyValue

    end subroutine MPAS_stream_mgr_set_property!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_get_property
    !
    !> \brief Queries the status of a stream property in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Returns the value of a property defined on a particular stream within
    !>  an MPAS stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_get_property(manager, streamID, propertyName, propertyValue, ierr)!{{{

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: propertyName
        integer, intent(out) :: propertyValue
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: stream
        integer, pointer :: curProperty

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_get_property()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_get_subpool(manager % streams, streamID, stream)

        call mpas_pool_get_config(stream, propertyName, curProperty)

        if (.not. associated(curProperty) ) then
            write(stderrUnit, *) 'ERROR: MPAS_stream_mgr_get_property property ', trim(propertyName), ' does not exist in stream ', trim(streamID)
            return
        end if

        propertyValue = curProperty

    end subroutine MPAS_stream_mgr_get_property!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_pkg
    !
    !> \brief Attach a package logical to the specified stream.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Attaches a package logical to a specific stream within an MPAS stream
    !>  manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_pkg(manager, streamID, packageName, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_pkg'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: packageName
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: stream
        type (MPAS_Pool_type), pointer :: streamPackages

        character (len=StrKIND) :: error

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_pkg()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_get_subpool(manager % streams, streamID, stream)

        if (.not. associated(stream)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not exist'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_get_subpool(stream, 'attachedPackages', streamPackages)

        if (.not. associated(streamPackages)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does have an attachedPackages pool'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_add_config(streamPackages, packageName, packageName)

    end subroutine MPAS_stream_mgr_add_pkg!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_remove_pkg
    !
    !> \brief Detaches a package logical from the specified stream.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Removes a package from a stream, so the package no longer controls
    !>  whether or not the stream is active.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_remove_pkg(manager, streamID, packageName, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_remove_pkg'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in), target :: packageName
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: stream
        type (MPAS_Pool_type), pointer :: streamPackages

        character (len=StrKIND) :: error

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_remove_pkg()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_get_subpool(manager % streams, streamID, stream)

        if (.not. associated(stream)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not exist in manager'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_get_subpool(stream, 'attachedPackages', streamPackages)

        if (.not. associated(streamPackages)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not contain an attachedPackages pool.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_remove_config(streamPackages, packageName)

    end subroutine MPAS_stream_mgr_remove_pkg!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_field
    !
    !> \brief Add a field to the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Adds a field from the allFields pool to a stream.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_field(manager, streamID, fieldName, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_field'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: fieldName
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: stream
        type (MPAS_Pool_Member_type), pointer :: newMem

        character (len=StrKIND) :: error

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_field()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_get_subpool(manager % streams, streamID, stream)

        if (.not. associated(stream)) then
           write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not exist.'
           call mpas_stream_mesg(manager % errorLevel, error)
           return
        end if

        allocate(newMem)
        newMem % key = trim(fieldName)
        newMem % keyLen = len_trim(fieldName)
        newMem % contentsType = MPAS_POOL_FIELD

        newMem % data => pool_get_member(manager % allFields, fieldName, MPAS_POOL_FIELD)

        if (.not. pool_add_member(stream, fieldName, newMem)) then
            write(error, *) 'ERROR: ', trim(sub), ' unable to add field ', trim(fieldName), ' to stream ', trim(streamID)
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

    end subroutine MPAS_stream_mgr_add_field!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_remove_field
    !
    !> \brief Remove a field from the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Removes a field from a stream.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_remove_field(manager, streamID, fieldName, ierr)!{{{
    
        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_remove_field'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: fieldName
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: stream

        character (len=StrKIND) :: error

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_remove_field()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_get_subpool(manager % streams, trim(streamID), stream)

        if (.not. associated(stream)) then
            write(error,*) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not exist.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_remove_field(stream, fieldName)

    end subroutine MPAS_stream_mgr_remove_field!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_att_int
    !
    !> \brief Add an integer attribute to the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Add a global integer attribute to the stream within an MPAS stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_att_int(manager, streamID, attName, attVal, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_att_int'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: attName
        integer, intent(in) :: attVal
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: stream, streamAttributes

        character (len=StrKIND) :: error

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_att()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_get_subpool(manager % streams, streamID, stream)

        if (.not. associated(stream)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not exist.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_get_subpool(stream, 'globalAttributes', streamAttributes)

        if (.not. associated(streamAttributes)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not contain a globalAttributes pool.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_add_config(streamAttributes, attName, attVal)


    end subroutine MPAS_stream_mgr_add_att_int!}}}

    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_att_real
    !
    !> \brief Add a real attribute to the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Add a global real attribute to the stream within an MPAS stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_att_real(manager, streamID, attName, attVal, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_att_real'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: attName
        real (kind=RKIND), intent(in) :: attVal
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: stream, streamAttributes

        character (len=StrKIND) :: error

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_att()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_get_subpool(manager % streams, streamID, stream)

        if (.not. associated(stream)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not exist.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_get_subpool(stream, 'globalAttributes', streamAttributes)

        if (.not. associated(streamAttributes)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not contain a globalAttributes pool.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_add_config(streamAttributes, attName, attVal)


    end subroutine MPAS_stream_mgr_add_att_real!}}}



    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_att_char
    !
    !> \brief Add a character attribute to the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Add a global character attribute to the stream within an MPAS stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_att_char(manager, streamID, attName, attVal, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_att_char'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: attName
        character (len=*), intent(in) :: attVal
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: stream, streamAttributes

        character (len=StrKIND) :: error

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_att()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_get_subpool(manager % streams, streamID, stream)

        if (.not. associated(stream)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not exist.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_get_subpool(stream, 'globalAttributes', streamAttributes)

        if (.not. associated(streamAttributes)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not contain a globalAttributes pool.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_add_config(streamAttributes, attName, attVal)


    end subroutine MPAS_stream_mgr_add_att_char!}}}



    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_att_logical
    !
    !> \brief Add a logical attribute to the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Add a global logical attribute to the stream within an MPAS stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_att_logical(manager, streamID, attName, attVal, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_att_logical'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: attName
        logical, intent(in) :: attVal
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: stream, streamAttributes
        character (len=StrKIND) :: attChar

        character (len=StrKIND) :: error

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_att()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_get_subpool(manager % streams, streamID, stream)

        if (.not. associated(stream)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not exist.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_get_subpool(stream, 'globalAttributes', streamAttributes)

        if (.not. associated(streamAttributes)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not contain a globalAttributes pool.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        if (attVal) then
            write(attChar, *) 'T'
        else
            write(attChar, *) 'F'
        end if

        call mpas_pool_add_config(streamAttributes, attName, trim(attChar))


    end subroutine MPAS_stream_mgr_add_att_logical!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_alarm
    !
    !> \brief Add an I/O alarm to a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  This routine will add a stream direction to be associated with an
    !>  alarm. It will not add the alarm to the manager's clock, but it is assumed
    !>  that the alarmID is used in the clock's alarm list.
    !>
    !>  It will create a subpool within the alarms pool that represents the
    !>  alarm (if it doesn't exist already). The pool representing this stream
    !>  will be added to the alarm pool, along with an integer that has the same
    !>  name as the stream whose value will represent the direction the stream
    !>  will be handled when this alarm rings.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_alarm(manager, streamID, alarmID, direction, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_alarm'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: alarmID
        integer, intent(in) :: direction
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: alarm, alarms
        type (MPAS_Pool_type), pointer :: stream

        character (len=StrKIND) :: error

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_alarm()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_get_subpool(manager % streams, streamID, stream)

        if (.not. associated(stream)) then
           write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not exist.'
           call mpas_stream_mesg(manager % errorLevel, error)
           return
        end if

        call mpas_pool_get_subpool(stream, 'alarms', alarms)

        if (.not. associated(alarms)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not have an alarms subpool.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_get_subpool(manager % alarms, alarmID, alarm)

        if (.not. associated(alarm)) then
            call mpas_pool_create_pool(alarm)
            call mpas_pool_add_subpool(manager % alarms, alarmID, alarm)
        end if

        ! Attach stream to alarm
        call mpas_pool_add_subpool(alarm, streamID, stream)
        call mpas_pool_add_config(alarm, streamID, direction)

        ! Attach alarm to stream
        call mpas_pool_add_subpool(alarms, alarmID, alarm)
        call mpas_pool_add_config(alarms, alarmID, direction)

    end subroutine MPAS_stream_mgr_add_alarm!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_remove_alarm
    !
    !> \brief Remove an I/O alarm from a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  This routine will remove the association of a stream to an alarm from
    !>  the stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_remove_alarm(manager, streamID, alarmID, direction, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_remove_alarm'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: alarmID
        integer, intent(in) :: direction
        integer, intent(out), optional :: ierr

        type (MPAS_Pool_type), pointer :: alarm, alarms
        type (MPAS_Pool_type), pointer :: stream

        character (len=StrKIND) :: error

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_remove_alarm()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_get_subpool(manager % alarms, alarmID, alarm)

        if (.not. associated(alarm)) then
            write(error, *) 'ERROR: ', trim(sub), ' alarm ', trim(alarmID), ' does not exist in manager.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_get_subpool(manager % streams, streamID, stream)

        if (.not. associated(stream)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not exist in manager.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        call mpas_pool_get_subpool(stream, 'alarms', alarms)

        if (.not. associated(alarms)) then
            write(error, *) 'ERROR: ', trim(sub), ' stream ', trim(streamID), ' does not have an alarms subpool.'
            call mpas_stream_mesg(manager % errorLevel, error)
            return
        end if

        ! Remove stream from alarm
        call mpas_pool_remove_subpool(alarm, streamID)
        call mpas_pool_remove_config(alarm, streamID)

        ! Remove alarm from stream
        call mpas_pool_remove_subpool(alarms, alarmID)
        call mpas_pool_remove_config(alarms, alarmID)

    end subroutine MPAS_stream_mgr_remove_alarm!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_write
    !
    !> \brief Write streams that are managed by an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  With no optional arguments, writes all streams whose alarms are ringing.
    !>  The "streamID" argument optionally specifies the ID of a particular stream
    !>  to be written; if no other optional arguments are given, the specified
    !>  stream is only written if any of its alarms are ringing.
    !>  The "timeLevel" argument optionally specifies, for fields with multiple
    !>  time levels, the time level from which fields should be written.
    !>  The "mgLevel" argument optionally specifies, for fields that exist for 
    !>  multiple grid levels, the grid level from which fields should be written.
    !>  The "forceWriteNow" argument optionally specifies that all streams -- or 
    !>  the stream specified by the "streamID" argument -- should be written by 
    !>  the call regardless of whether any alarms associated with the stream(s) 
    !>  are ringing.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_write(manager, streamID, timeLevel, mgLevel, forceWriteNow, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in), optional :: streamID
        integer, intent(in), optional :: timeLevel
        integer, intent(in), optional :: mgLevel
        logical, intent(in), optional :: forceWriteNow
        integer, intent(out), optional :: ierr

        integer :: iAlarm, nRinging
        character (len=ShortStrKIND), dimension(MPAS_MAX_ALARMS) :: ringingAlarms

        type (MPAS_Pool_type), pointer :: alarm, alarms
        type (MPAS_Pool_type), pointer :: stream
        type (MPAS_Pool_iterator_type) :: itr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_write()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        ! Handle all streams that are rining
        if (.not. present(streamID)) then
            call mpas_get_clock_ringing_alarms(manager % streamClock, nRinging, ringingAlarms)

            do iAlarm = 1, nRinging
                call mpas_pool_get_subpool(manager % alarms, ringingAlarms(iAlarm), alarm)

                if (associated(alarm)) then
                    
                end if
            end do
        end if

        

    end subroutine MPAS_stream_mgr_write


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_read
    !
    !> \brief Read streams that are managed by an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  With no optional arguments, reads all streams whose alarms are ringing.
    !>  The "streamID" argument optionally specifies the ID of a particular stream
    !>  to be read; if no other optional arguments are given, the specified stream
    !>  is only read if any of its alarms are ringing.
    !>  The "timeLevel" argument optionally specifies, for fields with multiple
    !>  time levels, the time level into which fields should be read.
    !>  The "mgLevel" argument optionally specifies, for fields that exist for 
    !>  multiple grid levels, the grid level into which fields should be read.
    !>  The "when" argument optionally specifies the timestamp from which fields
    !>  are to be read.
    !>  The "whence" argument optionally specifies the method for determining
    !>  the timestamp to read from in case an exact match is not found for the
    !>  read timestamp, which is the current time unless the optional "when"
    !>  argument is given; possible values are MPAS_STREAM_EXACT_TIME, 
    !>  MPAS_STREAM_NEAREST, MPAS_STREAM_LATEST_BEFORE, 
    !>  MPAS_STREAM_LATEST_STRICTLY_BEFORE, MPAS_STREAM_EARLIEST_AFTER, or 
    !>  MPAS_STREAM_EARLIEST_STRICTLY_AFTER.
    !>  The optional output argument "actualWhen" returns the actual time read 
    !>  from a stream in case an exact match for the "when" time is not found, 
    !>  and a nearby time is selected using the "whence" argument.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_read(manager, streamID, timeLevel, mgLevel, rightNow, when, whence, actualWhen, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in), optional :: streamID
        integer, intent(in), optional :: timeLevel
        integer, intent(in), optional :: mgLevel
        logical, intent(in), optional :: rightNow
        character (len=*), intent(in), optional :: when
        integer, intent(in), optional :: whence
        character (len=*), intent(out), optional :: actualWhen
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_read()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR
        if (present(actualWhen)) write(actualWhen,'(a)') '0000-01-01_00:00:00'


    end subroutine MPAS_stream_mgr_read


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mesg
    !
    !> \brief Write an error message (if the level requires it) to 
    !> \author Michael Duda, Doug Jacobsen
    !> \date   07/16/2014
    !> \details Using the input error level, 
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mesg(level, mesg)!{{{
        
        use mpas_dmpar

        integer, intent(in) :: level
        character(len=*), intent(in) :: mesg

        if (level /= MPAS_STREAM_ERR_SILENT) then
            write(stderrUnit, *) trim(mesg)
            if  (level == MPAS_STREAM_ERR_FATAL) then
                call mpas_dmpar_global_abort(mesg)
            end if
        end if

    end subroutine MPAS_stream_mesg!}}}


    !-----------------------------------------------------------------------
    !  routine mpas_stream_mgr_build_stream
    !
    !> \brief This is a utility routine to build a stream type from a pool representing a stream.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   07/23/2014
    !> \details This routine will take as input a pool representing a stream.
    !>  It will then generate a stream type based on this pool, and return that.
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_build_stream(timeStamp, inPool, stream)!{{{
        character (len=*), intent(in) :: timeStamp
        type (MPAS_Pool_type), intent(inout) :: inPool
        type (MPAS_Stream_type), pointer :: stream

        type (MPAS_Pool_iterator_type) :: itr

        character (len=StrKIND), pointer :: filenameTemplate
        integer, pointer :: direction
        integer, pointer :: maxRecords

        call mpas_pool_get_config(inPool, MPAS_STREAM_PROPERTY_FILENAME, filenameTemplate)
        call mpas_pool_get_config(inPool, MPAS_STREAM_PROPERTY_DIRECTION, direction)
        call mpas_pool_get_config(inPool, MPAS_STREAM_PROPERTY_MAX_RECORDS, maxRecords)

        do while ( mpas_pool_get_next_member(inPool, itr) )
            
        end do

    end subroutine MPAS_stream_mgr_build_stream!}}}


   
end module mpas_stream_manager


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! C interface routines for building streams at run-time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine stream_mgr_create_stream_c(mgr, str, i) bind(c)!{{{

    use mpas_c_interfacing, only : mpas_c_to_f_string
    use iso_c_binding, only : c_char, c_int, c_ptr, c_f_pointer
    use mpas_stream_manager, only : MPAS_streamManager_type, MPAS_STREAM_OUTPUT, MPAS_stream_mgr_create_stream

    implicit none

    type (c_ptr) :: mgr
    character(kind=c_char) :: str(*)
    integer(kind=c_int) :: i

    character(len=256) :: fieldname
    integer :: i_fort
    integer :: ierr
    type (MPAS_streamManager_type), pointer :: manager

    call mpas_c_to_f_string(str, fieldname)
    i_fort = i
 write(0,*) 'Doing something with string '//trim(fieldname)
    call c_f_pointer(mgr, manager)
    
    call MPAS_stream_mgr_create_stream(manager, fieldname, MPAS_STREAM_OUTPUT, fieldname, 1, ierr=ierr)

end subroutine stream_mgr_create_stream_c!}}}
