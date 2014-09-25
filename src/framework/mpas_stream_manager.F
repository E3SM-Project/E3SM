module mpas_stream_manager

#define COMMA ,
#define STREAM_DEBUG_WRITE(M) write(stderrUnit,*) M

    use mpas_kind_types
    use mpas_grid_types
    use mpas_timekeeping
    use mpas_io_units
    use mpas_io_streams
    use mpas_stream_list


    integer, public, parameter :: MPAS_STREAM_ERR_FATAL  = 1000, &
                                  MPAS_STREAM_ERR_WARN   = 1001, &
                                  MPAS_STREAM_ERR_SILENT = 1002

    integer, public, parameter :: MPAS_STREAM_MGR_NOERR =  0, &
                                  MPAS_STREAM_MGR_ERROR = -1

    integer, public, parameter :: MPAS_STREAM_INPUT = 1, &
                                  MPAS_STREAM_OUTPUT = 2, &
                                  MPAS_STREAM_INPUT_OUTPUT = 3, &
                                  MPAS_STREAM_NONE = 4

    integer, public, parameter :: MPAS_STREAM_PROPERTY_ACTIVE      =  5, &
                                  MPAS_STREAM_PROPERTY_IMMUTABLE   =  6, &
                                  MPAS_STREAM_PROPERTY_FILENAME    =  7, &
                                  MPAS_STREAM_PROPERTY_MAX_RECORDS =  8, &
                                  MPAS_STREAM_PROPERTY_REF_TIME    =  9, &
                                  MPAS_STREAM_PROPERTY_RECORD_INTV = 10

    public :: MPAS_streamManager_type, &
              MPAS_stream_mgr_init, &
              MPAS_stream_mgr_finalize, &
              MPAS_stream_mgr_create_stream, &
              MPAS_stream_mgr_destroy_stream, &
              MPAS_stream_mgr_get_clock, &
              MPAS_stream_mgr_set_property, &
              MPAS_stream_mgr_get_property, &
              MPAS_stream_mgr_add_pkg, &
              MPAS_stream_mgr_remove_pkg, &
              MPAS_stream_mgr_add_field, &
              MPAS_stream_mgr_remove_field, &
              MPAS_stream_mgr_add_alarm, &
              MPAS_stream_mgr_remove_alarm, &
              MPAS_stream_mgr_reset_alarms, &
              MPAS_stream_mgr_add_att, &
              MPAS_stream_mgr_write, &
              MPAS_stream_mgr_read

    private

    type MPAS_streamManager_type

        integer :: numStreams = 0
        integer :: errorLevel

        type (MPAS_Clock_type), pointer :: streamClock
        type (MPAS_Pool_type), pointer :: allFields
        type (MPAS_Pool_type), pointer :: allPackages
        type (MPAS_Pool_type), pointer :: defaultAtts

        type (MPAS_stream_list_type), pointer :: streams
        type (MPAS_stream_list_type), pointer :: alarms_in
        type (MPAS_stream_list_type), pointer :: alarms_out

    end type MPAS_streamManager_type


    interface MPAS_stream_mgr_set_property
        module procedure MPAS_stream_mgr_set_property_int
        module procedure MPAS_stream_mgr_set_property_char
        module procedure MPAS_stream_mgr_set_property_logical
    end interface


    interface MPAS_stream_mgr_get_property
        module procedure MPAS_stream_mgr_get_property_int
        module procedure MPAS_stream_mgr_get_property_char
        module procedure MPAS_stream_mgr_get_property_logical
    end interface


    interface MPAS_stream_mgr_add_att
        module procedure MPAS_stream_mgr_add_att_int
        module procedure MPAS_stream_mgr_add_att_real
        module procedure MPAS_stream_mgr_add_att_char
        module procedure MPAS_stream_mgr_add_att_logical
    end interface


    !
    ! Used for reindexing connectivity arrays during stream writes by the routines prewrite_reindex() and postwrite_reindex().
    ! Before a stream is written, we set the pointers here to be the heads of linked lists of locally-indexed connectivity fields.
    ! After a stream is written, we reset the arrays for connectivity fields in the stream to these pointers.
    !
    type (field2DInteger), pointer :: cellsOnCell_save
    type (field2DInteger), pointer :: edgesOnCell_save
    type (field2DInteger), pointer :: verticesOnCell_save
    type (field2DInteger), pointer :: cellsOnEdge_save
    type (field2DInteger), pointer :: verticesOnEdge_save
    type (field2DInteger), pointer :: edgesOnEdge_save
    type (field2DInteger), pointer :: cellsOnVertex_save
    type (field2DInteger), pointer :: edgesOnVertex_save


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
    subroutine MPAS_stream_mgr_init(manager, clock, allFields, allPackages, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_init'

        type (MPAS_streamManager_type), pointer :: manager
        type (MPAS_Clock_type), pointer :: clock
        type (MPAS_Pool_type), pointer :: allFields
        type (MPAS_Pool_type), pointer :: allPackages
        integer, intent(out), optional :: ierr

        integer :: err_local


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_init()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        allocate(manager)
        manager % allFields => allFields
        manager % allPackages => allPackages
        manager % streamClock => clock
        manager % numStreams = 0
        manager % errorLevel = MPAS_STREAM_ERR_SILENT

        !
        ! Set up linked list of streams
        !
        call MPAS_stream_list_create(manager % streams, ierr=err_local)
        if (err_local /= MPAS_STREAM_LIST_NOERR) then
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            STREAM_DEBUG_WRITE(' -- ERROR: Problems while creating stream list')
            return 
        end if

        !
        ! Set up linked list of input alarms
        !
        call MPAS_stream_list_create(manager % alarms_in, ierr=err_local)
        if (err_local /= MPAS_STREAM_LIST_NOERR) then
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            STREAM_DEBUG_WRITE(' -- ERROR: Problems while creating input alarm list')
            return 
        end if

        !
        ! Set up linked list of output alarms
        !
        call MPAS_stream_list_create(manager % alarms_out, ierr=err_local)
        if (err_local /= MPAS_STREAM_LIST_NOERR) then
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            STREAM_DEBUG_WRITE(' -- ERROR: Problems while creating output alarm list')
            return 
        end if

        !
        ! Create a pool to hold default global attributes that every stream will have
        !
        call mpas_pool_create_pool(manager % defaultAtts)

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

        integer :: err_local
        type (MPAS_stream_list_type), pointer :: stream_cursor


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_finalize()') 

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Remove all streams
        !
        stream_cursor => manager % streams % head
        do while (associated(stream_cursor))
            STREAM_DEBUG_WRITE(' -- deleting stream '//trim(stream_cursor % name))
            call MPAS_stream_mgr_destroy_stream(manager, stream_cursor % name, ierr=err_local)
            stream_cursor => manager % streams % head
        end do

        !
        ! Free up list of streams
        !
        call MPAS_stream_list_destroy(manager % streams, ierr=err_local)
        if (err_local /= MPAS_STREAM_LIST_NOERR) then
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            STREAM_DEBUG_WRITE(' -- ERROR: Problems while destroying stream list')
        end if

        !
        ! Free up list of input alarms
        !
        call MPAS_stream_list_destroy(manager % alarms_in, ierr=err_local)
        if (err_local /= MPAS_STREAM_LIST_NOERR) then
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            STREAM_DEBUG_WRITE(' -- ERROR: Problems while destroying input alarms list')
        end if

        !
        ! Free up list of output alarms
        !
        call MPAS_stream_list_destroy(manager % alarms_out, ierr=err_local)
        if (err_local /= MPAS_STREAM_LIST_NOERR) then
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            STREAM_DEBUG_WRITE(' -- ERROR: Problems while destroying output alarms list')
        end if

        !
        ! Free up default attribute pool
        !
        call mpas_pool_destroy_pool(manager % defaultAtts)

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
    !>  time. Additionally, the interval between records in the file may be
    !>  specified using the optional "recordInterval" argument; if this argument
    !>  is not supplied, the stream manager will assume that this interval is
    !>  equal to the shortest period of any periodic alarm attached to the stream.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_create_stream(manager, streamID, direction, filename, maxRecords, &
                                             referenceTime, recordInterval, ierr) !{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_create_stream'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: direction
        character (len=*), intent(in) :: filename
        integer, intent(in) :: maxRecords
        type (MPAS_Time_type), intent(in), optional :: referenceTime
        type (MPAS_TimeInterval_type), intent(in), optional :: recordInterval
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: new_stream
        integer :: err_local


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_create_stream() for '//trim(streamID))

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Check that the stream does not already exist
        !
        if (MPAS_stream_list_query(manager % streams, streamID, new_stream, ierr=err_local)) then
            STREAM_DEBUG_WRITE('-- Stream '//trim(streamID)//' already exist in stream manager')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Allocate a stream node to store the new stream
        !
        allocate(new_stream)
        new_stream % name = streamID
        new_stream % direction = direction
        new_stream % valid = .false.
!TODO: ensure that filename does not contain ':' characters, which PNETCDF does not like...
        new_stream % filename_template = filename
        new_stream % maxRecords = maxRecords
        new_stream % nRecords = 0
        allocate(new_stream % referenceTime)
        if (present(referenceTime)) then
            new_stream % referenceTime = referenceTime
        else
            new_stream % referenceTime = mpas_get_clock_time(manager % streamClock, MPAS_START_TIME)
        end if
        if (present(recordInterval)) then
            allocate(new_stream % recordInterval)
            new_stream % recordInterval = recordInterval
        end if
        call MPAS_stream_list_create(new_stream % alarmList_in, ierr=err_local)
        if (err_local /= MPAS_STREAM_LIST_NOERR) then
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            STREAM_DEBUG_WRITE(' -- ERROR: Problems while creating input alarm list')
            deallocate(new_stream)
            return 
        end if
        call MPAS_stream_list_create(new_stream % alarmList_out, ierr=err_local)
        if (err_local /= MPAS_STREAM_LIST_NOERR) then
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            STREAM_DEBUG_WRITE(' -- ERROR: Problems while creating output alarm list')
            deallocate(new_stream)
            return 
        end if
        call mpas_pool_create_pool(new_stream % att_pool)
        call mpas_pool_clone_pool(manager % defaultAtts, new_stream % att_pool)
        call mpas_pool_create_pool(new_stream % field_pool)
        call mpas_pool_create_pool(new_stream % pkg_pool)
        nullify(new_stream % next)


        !
        ! Add stream to list
        !
        call MPAS_stream_list_insert(manager % streams, new_stream, ierr=err_local)
        if (err_local /= MPAS_STREAM_LIST_NOERR) then
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            STREAM_DEBUG_WRITE(' -- ERROR: Problems while adding stream to list')
            return
        end if
        
        manager % numStreams = manager % numStreams + 1

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

        integer :: err_local
        type (MPAS_stream_list_type), pointer :: stream, alarm_cursor, delete_me


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_destroy_stream()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Remove stream from list
        !
        call MPAS_stream_list_remove(manager % streams, streamID, stream, ierr=err_local)
        if (err_local /= MPAS_STREAM_LIST_NOERR) then
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            STREAM_DEBUG_WRITE(' -- ERROR: Problems while removing stream from list')
            return
        end if

        !
        ! Unlink stream from input alarms
        !
        alarm_cursor => stream % alarmList_in % head
        do while (associated(alarm_cursor))
            call MPAS_stream_list_remove(alarm_cursor % xref % streamList, streamID, delete_me, ierr=err_local)
            if (err_local == MPAS_STREAM_LIST_NOERR) then
                deallocate(delete_me)
            else
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                STREAM_DEBUG_WRITE(' -- ERROR: Problems while removing stream from list of input alarm')
                return
            end if
            alarm_cursor => alarm_cursor % next
        end do

        !
        ! Unlink stream from output alarms
        !
        alarm_cursor => stream % alarmList_out % head
        do while (associated(alarm_cursor))
            call MPAS_stream_list_remove(alarm_cursor % xref % streamList, streamID, delete_me, ierr=err_local)
            if (err_local == MPAS_STREAM_LIST_NOERR) then
                deallocate(delete_me)
            else
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                STREAM_DEBUG_WRITE(' -- ERROR: Problems while removing stream from list of output alarm')
                return
            end if
            alarm_cursor => alarm_cursor % next
        end do

        !
        ! Free up stream storage -- reverse of whatever was done when allocating the stream
        !
        call MPAS_stream_list_destroy(stream % alarmList_in, ierr=err_local)
        call MPAS_stream_list_destroy(stream % alarmList_out, ierr=err_local)
        call mpas_pool_destroy_pool(stream % att_pool)
        call mpas_pool_destroy_pool(stream % field_pool)
        call mpas_pool_destroy_pool(stream % pkg_pool)
        if (associated(stream % referenceTime)) then
            deallocate(stream % referenceTime)
        end if
        if (associated(stream % recordInterval)) then
            deallocate(stream % recordInterval)
        end if
        if (stream % valid) then
            call MPAS_closeStream(stream % stream, ierr=err_local)
            if (err_local /= MPAS_STREAM_NOERR) then
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                STREAM_DEBUG_WRITE(' -- ERROR: Problems while closing stream '//trim(stream % name))
            end if
            deallocate(stream % stream)
        end if
        deallocate(stream)
        
        manager % numStreams = manager % numStreams - 1

    end subroutine MPAS_stream_mgr_destroy_stream!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_get_clock
    !
    !> \brief Retrieves the clock used by the stream manager.
    !> \author Michael Duda
    !> \date   22 August 2014
    !> \details
    !>  Returns a pointer to the clock associated with the stream manager, 
    !>  in which any stream alarms should be defined before being added to 
    !>  the stream manager via the MPAS_stream_mgr_add_alarm() routine.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_get_clock(manager, clock, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(in) :: manager
        type (MPAS_Clock_type), pointer :: clock
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_get_clock()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        clock => manager % streamClock

    end subroutine MPAS_stream_mgr_get_clock


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

        type (MPAS_stream_list_type), pointer :: stream
        type (mpas_pool_field_info_type) :: info
        integer, pointer :: test_ptr
        integer :: err_local


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_field()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Check that stream exists
        !
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream, ierr=err_local)) then
            STREAM_DEBUG_WRITE('-- Requested stream '//trim(streamID)//' does not exist in stream manager')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Don't modify an immutable stream
        !
        if (stream % immutable) then
            STREAM_DEBUG_WRITE('-- Requested stream '//trim(streamID)//' is immutable.')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Check that field exists
        !
        info % nDims = -1
        call mpas_pool_get_field_info(manager % allFields, fieldName, info)
        if (info % nDims == -1) then
            STREAM_DEBUG_WRITE('-- Requested field '//trim(fieldName)//' not available')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Check that the field does not already exist in the stream
        !
        nullify(test_ptr)
        call mpas_pool_get_config(stream % field_pool, fieldName, value=test_ptr)
        if (associated(test_ptr)) then
            STREAM_DEBUG_WRITE('-- Requested field '//trim(fieldName)//' already in stream '//trim(streamID))
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Add field to field pool in stream
        !
        call mpas_pool_add_config(stream % field_pool, fieldName, 1)

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

        type (MPAS_stream_list_type), pointer :: stream
        integer, pointer :: test_ptr
        integer :: err_local


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_remove_field()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Check that stream exists
        !
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream, ierr=err_local)) then
            STREAM_DEBUG_WRITE('-- Requested stream '//trim(streamID)//' does not exist in stream manager')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Don't modify an immutable stream
        !
        if (stream % immutable) then
            STREAM_DEBUG_WRITE('-- Requested stream '//trim(streamID)//' is immutable.')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Check that field exists in stream's field pool
        !
        nullify(test_ptr)
        call mpas_pool_get_config(stream % field_pool, fieldName, value=test_ptr)
        if (.not. associated(test_ptr)) then
            STREAM_DEBUG_WRITE('-- Requested field '//trim(fieldName)//' not in stream '//trim(streamID))
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Remove field from stream's field pool
        !
        call mpas_pool_remove_config(stream % field_pool, fieldName)

    end subroutine MPAS_stream_mgr_remove_field!}}}


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

        type (MPAS_stream_list_type), pointer :: stream, new_alarm, new_xref
        integer :: err_local


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_alarm()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Check that stream exists
        !
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream, ierr=err_local)) then
            STREAM_DEBUG_WRITE('-- Requested stream '//trim(streamID)//' does not exist in stream manager')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Check that the specified direction makes sense for the stream
        !
        if (stream % direction == MPAS_STREAM_OUTPUT .and. direction == MPAS_STREAM_INPUT .or. &
            stream % direction == MPAS_STREAM_OUTPUT .and. direction == MPAS_STREAM_INPUT_OUTPUT .or. &
            stream % direction == MPAS_STREAM_INPUT .and. direction == MPAS_STREAM_OUTPUT .or. &
            stream % direction == MPAS_STREAM_INPUT .and. direction == MPAS_STREAM_INPUT_OUTPUT .or. &
            stream % direction == MPAS_STREAM_NONE) then

            STREAM_DEBUG_WRITE('-- Attempting to add an alarm '//trim(alarmID)//' to invalid direction for stream '//trim(streamID))
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Check that alarm exists on clock
        !
        if (.not. mpas_is_alarm_defined(manager % streamClock, alarmID, err_local)) then
            STREAM_DEBUG_WRITE('-- Attempting to add an alarm '//trim(alarmID)//' that does not exist on clock')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Check that the alarm does not already exist for the stream in the specified direction
        !
        if (direction == MPAS_STREAM_INPUT .or. direction == MPAS_STREAM_INPUT_OUTPUT) then
            if (MPAS_stream_list_query(stream % alarmList_in, alarmID, new_alarm, ierr=err_local)) then
                STREAM_DEBUG_WRITE('-- Requested input alarm '//trim(alarmID)//' already on stream '//trim(streamID))
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return
            end if
        end if
        if (direction == MPAS_STREAM_OUTPUT .or. direction == MPAS_STREAM_INPUT_OUTPUT) then
            if (MPAS_stream_list_query(stream % alarmList_out, alarmID, new_alarm, ierr=err_local)) then
                STREAM_DEBUG_WRITE('-- Requested output alarm '//trim(alarmID)//' already on stream '//trim(streamID))
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return
            end if
        end if


        !
        ! Add alarm to alarm to the alarms_in and/or alarms_out list
        ! Add alarm to the alarmList_in and/or alarmList_out list for the field
        !
        if (direction == MPAS_STREAM_INPUT .or. direction == MPAS_STREAM_INPUT_OUTPUT) then

            ! If alarm is not already defined, we need to create a new alarm node
            if (.not. MPAS_stream_list_query(manager % alarms_in, alarmID, new_alarm, ierr=err_local)) then
                allocate(new_alarm)
                new_alarm % name = alarmID
                call MPAS_stream_list_create(new_alarm % streamList, ierr=err_local)
                if (err_local /= MPAS_STREAM_LIST_NOERR) then
                    if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                    STREAM_DEBUG_WRITE(' -- ERROR: Problems while creating stream list for alarm')
                    return 
                end if
                nullify(new_alarm % next)

                call MPAS_stream_list_insert(manager % alarms_in, new_alarm, ierr=err_local)
                if (err_local /= MPAS_STREAM_LIST_NOERR) then
                    if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                    STREAM_DEBUG_WRITE(' -- ERROR: Problems while adding input alarm to list')
                    return
                end if
            end if

            ! Add specified stream to alarm node stream list
            allocate(new_xref)
            new_xref % name = streamID
            new_xref % xref => stream
            call MPAS_stream_list_insert(new_alarm % streamList, new_xref, ierr=err_local)
            if (err_local /= MPAS_STREAM_LIST_NOERR) then
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                STREAM_DEBUG_WRITE(' -- ERROR: Problems while adding stream to alarm stream list')
                return
            end if

            ! Add alarm to stream alarm list
            allocate(new_xref)
            new_xref % name = alarmID
            new_xref % xref => new_alarm
            call MPAS_stream_list_insert(stream % alarmList_in, new_xref, ierr=err_local)
            if (err_local /= MPAS_STREAM_LIST_NOERR) then
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                STREAM_DEBUG_WRITE(' -- ERROR: Problems while adding alarm to stream input alarm list')
                return
            end if
        end if

        if (direction == MPAS_STREAM_OUTPUT .or. direction == MPAS_STREAM_INPUT_OUTPUT) then

            ! If alarm is not already defined, we need to create a new alarm node
            if (.not. MPAS_stream_list_query(manager % alarms_out, alarmID, new_alarm, ierr=err_local)) then
                allocate(new_alarm)
                new_alarm % name = alarmID
                call MPAS_stream_list_create(new_alarm % streamList, ierr=err_local)
                if (err_local /= MPAS_STREAM_LIST_NOERR) then
                    if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                    STREAM_DEBUG_WRITE(' -- ERROR: Problems while creating stream list for alarm')
                    return 
                end if
                nullify(new_alarm % next)

                call MPAS_stream_list_insert(manager % alarms_out, new_alarm, ierr=err_local)
                if (err_local /= MPAS_STREAM_LIST_NOERR) then
                    if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                    STREAM_DEBUG_WRITE(' -- ERROR: Problems while adding output alarm to list')
                    return
                end if
            end if

            ! Add specified stream to alarm node stream list
            allocate(new_xref)
            new_xref % name = streamID
            new_xref % xref => stream
            call MPAS_stream_list_insert(new_alarm % streamList, new_xref, ierr=err_local)
            if (err_local /= MPAS_STREAM_LIST_NOERR) then
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                STREAM_DEBUG_WRITE(' -- ERROR: Problems while adding stream to alarm stream list')
                return
            end if

            ! Add alarm to stream alarm list
            allocate(new_xref)
            new_xref % name = alarmID
            new_xref % xref => new_alarm
            call MPAS_stream_list_insert(stream % alarmList_out, new_xref, ierr=err_local)
            if (err_local /= MPAS_STREAM_LIST_NOERR) then
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                STREAM_DEBUG_WRITE(' -- ERROR: Problems while adding alarm to stream output alarm list')
                return
            end if
        end if

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

        type (MPAS_stream_list_type), pointer :: stream
        type (MPAS_stream_list_type), pointer :: alarmNode
        type (MPAS_stream_list_type), pointer :: streamNode
        integer :: err_local


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_remove_alarm()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Check that stream exists
        !
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream, ierr=err_local)) then
            STREAM_DEBUG_WRITE('-- Requested stream '//trim(streamID)//' does not exist in stream manager')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Unlink alarm from alarmList_in or alarmList_out for stream
        !
        nullify(alarmNode)
        if (direction == MPAS_STREAM_INPUT) then
            call MPAS_stream_list_remove(stream % alarmList_in, alarmID, alarmNode, ierr=ierr)
        else if (direction == MPAS_STREAM_OUTPUT) then
            call MPAS_stream_list_remove(stream % alarmList_out, alarmID, alarmNode, ierr=ierr)
        else
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            STREAM_DEBUG_WRITE('-- Requested to remove alarm from invalid direction from stream '//trim(streamID))
            return
        end if

        !
        ! Remove stream from alarm's streamList in alarms_in or alarms_out
        !
        if (associated(alarmNode)) then
            call MPAS_stream_list_remove(alarmNode % xref % streamList, streamID, streamNode, ierr=ierr)
        else
            if (direction == MPAS_STREAM_INPUT) then
                STREAM_DEBUG_WRITE('-- Input alarm '//trim(alarmID)//' does not exist on stream '//trim(streamID))
            else
                STREAM_DEBUG_WRITE('-- Output alarm '//trim(alarmID)//' does not exist on stream '//trim(streamID))
            end if
            return
        end if 
        if (.not. associated(streamNode)) then
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            STREAM_DEBUG_WRITE('-- Alarm '//trim(alarmID)//' does not have stream '//trim(streamID)//' on its stream list.')
            return
        end if

        !
        ! If the alarm has no associated streams, should we remove it from alarms_in or alarms_out?
        !
        if (MPAS_stream_list_length(alarmNode % xref % streamList) == 0) then
            if (direction == MPAS_STREAM_INPUT) then
                STREAM_DEBUG_WRITE('-- Input alarm '//trim(alarmID)//' has no associated streams and will be deleted.')
                call MPAS_stream_list_remove(manager % alarms_in, alarmID, alarmNode, ierr=ierr)
            else
                STREAM_DEBUG_WRITE('-- Output alarm '//trim(alarmID)//' has no associated streams and will be deleted.')
                call MPAS_stream_list_remove(manager % alarms_out, alarmID, alarmNode, ierr=ierr)
            end if
        end if

    end subroutine MPAS_stream_mgr_remove_alarm!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_reset_alarms
    !
    !> \brief Reset I/O alarms in a stream manager
    !> \author Michael Duda
    !> \date   2 September 2014
    !> \details
    !>  Resets all alarms used by the stream manager. If the optional argument
    !>  'streamID' is provided, only alarms associated with that stream will be
    !>  reset. If the optional 'direction' argument is provided, only alarms 
    !>  associated with that direction will be reset.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_reset_alarms(manager, streamID, direction, ierr)!{{{

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in), optional :: streamID
        integer, intent(in), optional :: direction
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream
        type (MPAS_stream_list_type), pointer :: alarm_cursor
        integer :: local_direction
        integer :: local_ierr

        
        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_reset_alarms()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR


        !
        ! Check for optional direction argument; default direction is both input and output.
        !
        if (present(direction)) then
            local_direction = direction
        else
            local_direction = MPAS_STREAM_INPUT_OUTPUT
        end if


        !
        ! Check for optional streamID argument; default is to handle all alarms in the manager.
        !
        nullify(stream)
        if (present(streamID)) then
            if (.not. MPAS_stream_list_query(manager % streams, streamID, stream, ierr=local_ierr)) then
                STREAM_DEBUG_WRITE('-- Stream '//trim(streamID)//' does not exist in stream manager.')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return
            end if
        end if

        
        if (local_direction == MPAS_STREAM_INPUT .or. local_direction == MPAS_STREAM_INPUT_OUTPUT) then
            if (associated(stream)) then
                alarm_cursor => stream % alarmList_in % head
            else
                alarm_cursor => manager % alarms_in % head
            end if
            do while (associated(alarm_cursor))
                if (mpas_is_alarm_ringing(manager % streamClock, alarm_cursor % name, ierr=local_ierr)) then
                    call mpas_reset_clock_alarm(manager % streamClock, alarm_cursor % name, ierr=local_ierr) 
                end if
                alarm_cursor => alarm_cursor % next
            end do
        end if


        if (local_direction == MPAS_STREAM_OUTPUT .or. local_direction == MPAS_STREAM_INPUT_OUTPUT) then
            if (associated(stream)) then
                alarm_cursor => stream % alarmList_out % head
            else
                alarm_cursor => manager % alarms_out % head
            end if
            do while (associated(alarm_cursor))
                if (mpas_is_alarm_ringing(manager % streamClock, alarm_cursor % name, ierr=local_ierr)) then
                    call mpas_reset_clock_alarm(manager % streamClock, alarm_cursor % name, ierr=local_ierr) 
                end if
                alarm_cursor => alarm_cursor % next
            end do
        end if

    end subroutine MPAS_stream_mgr_reset_alarms!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_set_property_int
    !
    !> \brief Sets a property of a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Sets the value of a stream property within an MPAS stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_set_property_int(manager, streamID, propertyName, propertyValue, ierr) !{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_set_property_int'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: propertyName
        integer, intent(in) :: propertyValue
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_set_property()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
            STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_set_property().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if

        !
        ! Set property
        !
        select case (propertyName)

            case (MPAS_STREAM_PROPERTY_MAX_RECORDS)
                stream_cursor % maxRecords = propertyValue

            case default
                STREAM_DEBUG_WRITE('-- Error: MPAS_stream_mgr_set_property(): No such property ' COMMA propertyName)
                STREAM_DEBUG_WRITE('          or specified property is not of type integer.')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
        end select

    end subroutine MPAS_stream_mgr_set_property_int !}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_set_property_char
    !
    !> \brief Sets a property of a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Sets the value of a stream property within an MPAS stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_set_property_char(manager, streamID, propertyName, propertyValue, ierr) !{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_set_property_char'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: propertyName
        character (len=*), intent(in) :: propertyValue
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_set_property()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
            STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_set_property().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if

        !
        ! Set property
        !
        select case (propertyName)

            case (MPAS_STREAM_PROPERTY_FILENAME)
!TODO: ensure that filename does not contain ':' characters, which PNETCDF does not like...
                stream_cursor % filename_template = propertyValue

            case (MPAS_STREAM_PROPERTY_REF_TIME)
                call mpas_set_time(stream_cursor % referenceTime, dateTimeString=propertyValue)

            case (MPAS_STREAM_PROPERTY_RECORD_INTV)
 
                ! The interval between records may not have been allocated if the optional recordInterval
                !    argument was not provided when the stream was created
                if (.not. associated(stream_cursor % recordInterval)) then
                    allocate(stream_cursor % recordInterval)
                end if
                call mpas_set_timeInterval(stream_cursor % recordInterval, timeString=propertyValue)

            case default
                STREAM_DEBUG_WRITE('-- Error: MPAS_stream_mgr_set_property(): No such property ' COMMA propertyName)
                STREAM_DEBUG_WRITE('          or specified property is not of type character.')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
        end select

    end subroutine MPAS_stream_mgr_set_property_char !}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_set_property_logical
    !
    !> \brief Sets a property of a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Sets the value of a stream property within an MPAS stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_set_property_logical(manager, streamID, propertyName, propertyValue, ierr) !{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_set_property_logical'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: propertyName
        logical, intent(in) :: propertyValue
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_set_property()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
            STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_set_property().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if

        !
        ! Set property
        !
        select case (propertyName)

            case (MPAS_STREAM_PROPERTY_ACTIVE)
                stream_cursor % active_stream = propertyValue

            case (MPAS_STREAM_PROPERTY_IMMUTABLE)
                stream_cursor % immutable = propertyValue

            case default
                STREAM_DEBUG_WRITE('-- Error: MPAS_stream_mgr_set_property(): No such property ' COMMA propertyName)
                STREAM_DEBUG_WRITE('          or specified property is not of type logical.')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
        end select

    end subroutine MPAS_stream_mgr_set_property_logical !}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_get_property_int
    !
    !> \brief Sets a property of a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Sets the value of a stream property within an MPAS stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_get_property_int(manager, streamID, propertyName, propertyValue, ierr) !{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_set_property_int'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: propertyName
        integer, intent(out) :: propertyValue
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_get_property()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
            STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_get_property().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if

        !
        ! Set property
        !
        select case (propertyName)

            case (MPAS_STREAM_PROPERTY_MAX_RECORDS)
                propertyValue = stream_cursor % maxRecords

            case default
                STREAM_DEBUG_WRITE('-- Error: MPAS_stream_mgr_get_property(): No such property ' COMMA propertyName)
                STREAM_DEBUG_WRITE('          or specified property is not of type integer.')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
        end select

    end subroutine MPAS_stream_mgr_get_property_int !}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_get_property_char
    !
    !> \brief Sets a property of a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Sets the value of a stream property within an MPAS stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_get_property_char(manager, streamID, propertyName, propertyValue, ierr) !{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_get_property_char'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: propertyName
        character (len=*), intent(out) :: propertyValue
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        type (MPAS_stream_list_type), pointer :: alarm_cursor
        type (MPAS_timeInterval_type) :: temp_interval, read_interval
        integer :: nAlarms
        integer :: err_local

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_get_property()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
            STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_get_property().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if

        !
        ! Set property
        !
        select case (propertyName)

            case (MPAS_STREAM_PROPERTY_FILENAME)
                propertyValue = stream_cursor % filename_template

            case (MPAS_STREAM_PROPERTY_REF_TIME)
                call mpas_get_time(stream_cursor % referenceTime, dateTimeString=propertyValue)

            case (MPAS_STREAM_PROPERTY_RECORD_INTV)

                ! The interval between records may not have been allocated if the optional recordInterval
                !    argument was not provided when the stream was created. If there is no explicit recordInterval,
                !    assume that the interval is the shortest interval between alarms on the stream; since
                !    recordInterval is only used for reading, use the input alarm list in this check.
                if (.not. associated(stream_cursor % recordInterval)) then
                    alarm_cursor => stream_cursor % alarmList_in % head
                    nAlarms = 0
                    do while (associated(alarm_cursor))
                        temp_interval = mpas_alarm_interval(manager % streamClock, alarm_cursor % name, err_local)
                        if (err_local == 0) then
                            if (nAlarms == 0) then
                                read_interval = temp_interval
                            else if (temp_interval < read_interval) then
                                read_interval = temp_interval
                            end if
                            nAlarms = nAlarms + 1
                        end if
                        alarm_cursor => alarm_cursor % next
                    end do
                    if (nAlarms > 0) then
                        call mpas_get_timeInterval(read_interval, timeString=propertyValue)
                    else
                        propertyValue = '00:00:00'
                    end if
                else
                    call mpas_get_timeInterval(stream_cursor % recordInterval, timeString=propertyValue)
                end if

            case default
                STREAM_DEBUG_WRITE('-- Error: MPAS_stream_mgr_get_property(): No such property ' COMMA propertyName)
                STREAM_DEBUG_WRITE('          or specified property is not of type character.')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
        end select

    end subroutine MPAS_stream_mgr_get_property_char !}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_get_property_logical
    !
    !> \brief Sets a property of a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Sets the value of a stream property within an MPAS stream manager.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_get_property_logical(manager, streamID, propertyName, propertyValue, ierr) !{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_get_property'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: propertyName
        logical, intent(out) :: propertyValue
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_get_property_logical()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
            STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_get_property().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if

        !
        ! Set property
        !
        select case (propertyName)

            case (MPAS_STREAM_PROPERTY_ACTIVE)
                propertyValue = stream_cursor % active_stream

            case (MPAS_STREAM_PROPERTY_IMMUTABLE)
                propertyValue = stream_cursor % immutable

            case default
                STREAM_DEBUG_WRITE('-- Error: MPAS_stream_mgr_get_property(): No such property ' COMMA propertyName)
                STREAM_DEBUG_WRITE('          or specified property is not of type logical.')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
        end select

    end subroutine MPAS_stream_mgr_get_property_logical !}}}


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

        logical, pointer :: package
        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_pkg()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Query pointer to package in the manager-wide package pool
        !
        nullify(package)
        call mpas_pool_get_package(manager % allPackages, packageName, package)
        if (.not. associated(package)) then
            STREAM_DEBUG_WRITE('-- Error: package '//trim(packageName)//' not found in call to MPAS_stream_mgr_add_pkg().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Find requested stream
        !
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
            STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_add_pkg().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if

        !
        ! Add package to the packages pool for the stream
        !
        call mpas_pool_add_package(stream_cursor % pkg_pool, packageName, package)    

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

        logical, pointer :: package
        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_remove_pkg()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
            STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_remove_pkg().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if

        !
        ! Remove package from the packages pool for the stream
        !
        call mpas_pool_remove_package(stream_cursor % pkg_pool, packageName)

    end subroutine MPAS_stream_mgr_remove_pkg!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_att_int
    !
    !> \brief Add an integer attribute to the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Add a global integer attribute to the stream within an MPAS stream manager.
    !>  If the optional streamID argument is not supplied, the attribute will be
    !>  applied to every stream created after the call to add the attribute.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_att_int(manager, attName, attVal, streamID, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_att_int'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: attName
        integer, intent(in) :: attVal
        character (len=*), intent(in), optional :: streamID
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_att()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        if (present(streamID)) then
            if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
                STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_add_att().')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return 
            end if
    
            call mpas_pool_add_config(stream_cursor % att_pool, attName, attVal)
        else
            call mpas_pool_add_config(manager % defaultAtts, attName, attVal)
        end if

    end subroutine MPAS_stream_mgr_add_att_int!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_att_real
    !
    !> \brief Add a real attribute to the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Add a global real attribute to the stream within an MPAS stream manager.
    !>  If the optional streamID argument is not supplied, the attribute will be
    !>  applied to every stream created after the call to add the attribute.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_att_real(manager, attName, attVal, streamID, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_att_real'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: attName
        real (kind=RKIND), intent(in) :: attVal
        character (len=*), intent(in), optional :: streamID
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_att()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        if (present(streamID)) then
            if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
                STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_add_att().')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return 
            end if

            call mpas_pool_add_config(stream_cursor % att_pool, attName, attVal)
        else
            call mpas_pool_add_config(manager % defaultAtts, attName, attVal)
        end if

    end subroutine MPAS_stream_mgr_add_att_real!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_att_char
    !
    !> \brief Add a character attribute to the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Add a global character attribute to the stream within an MPAS stream manager.
    !>  If the optional streamID argument is not supplied, the attribute will be
    !>  applied to every stream created after the call to add the attribute.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_att_char(manager, attName, attVal, streamID, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_att_char'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: attName
        character (len=*), intent(in) :: attVal
        character (len=*), intent(in), optional :: streamID
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_att()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        if (present(streamID)) then
            if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
                STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_add_att().')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return 
            end if

            call mpas_pool_add_config(stream_cursor % att_pool, attName, attVal)
        else
            call mpas_pool_add_config(manager % defaultAtts, attName, attVal)
        end if

    end subroutine MPAS_stream_mgr_add_att_char!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_att_logical
    !
    !> \brief Add a logical attribute to the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Add a global logical attribute to the stream within an MPAS stream manager.
    !>  If the optional streamID argument is not supplied, the attribute will be
    !>  applied to every stream created after the call to add the attribute.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_att_logical(manager, attName, attVal, streamID, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_att_logical'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: attName
        logical, intent(in) :: attVal
        character (len=*), intent(in), optional :: streamID
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_att()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        if (present(streamID)) then
            if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
                STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_add_att().')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return 
            end if

            call mpas_pool_add_config(stream_cursor % att_pool, attName, attVal)
        else
            call mpas_pool_add_config(manager % defaultAtts, attName, attVal)
        end if

    end subroutine MPAS_stream_mgr_add_att_logical!}}}


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

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: local_timeLevel
        integer :: local_mgLevel
        logical :: local_forceWrite
        integer :: local_ierr
        integer :: temp_ierr


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_write()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Use optional arguments or set defaults
        !
        if (present(timeLevel)) then
            local_timeLevel = timeLevel
        else
            local_timeLevel = 1
        end if

        if (present(mgLevel)) then
            local_mgLevel = mgLevel
        else
            local_mgLevel = 1
        end if

        if (present(forceWriteNow)) then
            local_forceWrite = forceWriteNow
        else
            local_forceWrite = .false.
        end if


        !
        ! If a stream is specified, we process just that stream; otherwise,
        !    process all streams
        !
        if (present(streamID)) then
            nullify(stream_cursor)
            if (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=ierr)) then
                STREAM_DEBUG_WRITE('-- Handling write of stream '//trim(stream_cursor % name))

                ! Verify that the stream is an output stream
                if (stream_cursor % direction /= MPAS_STREAM_OUTPUT .and. &
                    stream_cursor % direction /= MPAS_STREAM_INPUT_OUTPUT) then
                    STREAM_DEBUG_WRITE('-- Stream '//streamID//' is not an output stream.')
                    if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                    return 
                end if

                call write_stream(manager, stream_cursor, local_timeLevel, local_mgLevel, local_forceWrite, local_ierr)
            else
                STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_write().')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return 
            end if
        else
            nullify(stream_cursor)
            stream_cursor => manager % streams % head
            do while (associated(stream_cursor))
                STREAM_DEBUG_WRITE('-- Handling write of stream '//trim(stream_cursor % name))

                ! Verify that the stream is an output stream
                if (stream_cursor % direction == MPAS_STREAM_OUTPUT .or. &
                    stream_cursor % direction == MPAS_STREAM_INPUT_OUTPUT) then

                    call write_stream(manager, stream_cursor, local_timeLevel, local_mgLevel, local_forceWrite, temp_ierr)
                    if (temp_ierr /= MPAS_STREAM_MGR_NOERR) then
                        local_ierr = MPAS_STREAM_MGR_ERROR
                    end if

                end if
                stream_cursor => stream_cursor % next
            end do
        end if

        if (present(ierr)) ierr = local_ierr

    end subroutine MPAS_stream_mgr_write


    !-----------------------------------------------------------------------
    !  routine write_stream
    !
    !> \brief Handle the writing of a stream pointed to by the stream list node
    !> \author Michael Duda
    !> \date   2 September 2014
    !> \details
    !>  Private subroutine to handle the details of actually writing a stream.
    !
    !-----------------------------------------------------------------------
    subroutine write_stream(manager, stream, timeLevel, mgLevel, forceWritenow, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        type (MPAS_stream_list_type), intent(inout) :: stream
        integer, intent(in) :: timeLevel
        integer, intent(in) :: mgLevel
        logical, intent(in) :: forceWriteNow
        integer, intent(out) :: ierr

        type (MPAS_stream_list_type), pointer :: alarm_cursor
        type (MPAS_Time_type) :: now_time
        type (MPAS_TimeInterval_type) :: temp_interval
        type (MPAS_TimeInterval_type) :: write_interval
        character (len=StrKIND) :: now_string
        character (len=StrKIND) :: temp_filename
        logical :: ringing_alarm
        integer :: nAlarms
        integer :: local_ierr


        ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Check whether this stream is active
        !
        if (.not. stream % active_stream) then
            STREAM_DEBUG_WRITE('-- Stream '//trim(stream % name)//' is not currently active and will not be written.')
            return
        end if

        !
        ! Check whether all packages for this stream are inactive
        ! Note: if the stream has no packages, it is assumed to be active
        !
        if (.not. stream_active_pkg_check(stream)) then
            STREAM_DEBUG_WRITE('-- Stream '//trim(stream % name)//' has only inactive packages and will not be written.')
            return
        end if

        !
        ! Check whether any of the output alarms for the stream are ringing
        !
        ringing_alarm = .false.
        alarm_cursor => stream % alarmList_out % head
        do while (associated(alarm_cursor))
            if (mpas_is_alarm_ringing(manager % streamClock, alarm_cursor % name, ierr=local_ierr)) then
                ringing_alarm = .true.
                exit
            end if
            alarm_cursor => alarm_cursor % next
        end do

        if ((.not. ringing_alarm) .and. (.not. forceWriteNow)) then
            return
        end if


        !
        ! Figure out the write interval for the stream, defined as the minimum
        !    period of any output alarm associated with the stream
        !
        alarm_cursor => stream % alarmList_out % head
        nAlarms = 0
        do while (associated(alarm_cursor))
            temp_interval = mpas_alarm_interval(manager % streamClock, alarm_cursor % name, local_ierr)
            if (local_ierr == 0) then
                if (nAlarms == 0) then
                    write_interval = temp_interval
                else if (temp_interval < write_interval) then
                    write_interval = temp_interval
                end if
                nAlarms = nAlarms + 1
            end if
            alarm_cursor => alarm_cursor % next
        end do


        !
        ! If the stream is not valid, assume that we have not yet written this
        ! stream, in which case we create the stream from scratch
        !
        if (.not. stream % valid) then
            now_time = mpas_get_clock_time(manager % streamClock, MPAS_NOW, ierr=local_ierr)
            call mpas_get_time(now_time, dateTimeString=now_string, ierr=local_ierr)
            if (stream % maxRecords == 0 .or. nAlarms == 0) then
                call mpas_expand_string(now_string, stream % filename_template, stream % filename)
                stream % nRecords = 1
            else
                call build_filename(stream % referenceTime, now_time, MPAS_STREAM_EXACT_TIME, stream % maxRecords, write_interval, &
                                    stream % filename_template, stream % filename, stream % nRecords, local_ierr)
            end if
            if (local_ierr /= 0) then
                if (local_ierr == 1) then
                    STREAM_DEBUG_WRITE('-- Error: Could not find a file offset that exactly matches')
                    STREAM_DEBUG_WRITE('--        combination of reference time, current time, and interval.')
                end if
                ierr = MPAS_STREAM_MGR_ERROR
                return
            end if

            !
            ! Build stream from pools of fields and attributes
            !
            allocate(stream % stream)
            call MPAS_createStream(stream % stream, stream % filename, MPAS_IO_PNETCDF, MPAS_IO_WRITE, stream % maxRecords, &
                                   ierr=local_ierr)
            if (local_ierr /= MPAS_STREAM_NOERR) then
                ierr = MPAS_STREAM_MGR_ERROR
                return
            end if

            call build_stream(stream, MPAS_STREAM_OUTPUT, manager % allFields, timeLevel, mgLevel, local_ierr)
            if (local_ierr /= MPAS_STREAM_NOERR) then
                ierr = MPAS_STREAM_MGR_ERROR
                return
            end if
            stream % timeLevel = timeLevel

            stream % valid = .true.
        else
            if (stream % maxRecords /= 0 .and. nAlarms /= 0) then
                now_time = mpas_get_clock_time(manager % streamClock, MPAS_NOW, ierr=local_ierr)
                call build_filename(stream % referenceTime, now_time, MPAS_STREAM_EXACT_TIME, stream % maxRecords, write_interval, &
                                    stream % filename_template, temp_filename, stream % nRecords, local_ierr)
                if (local_ierr /= 0) then
                    if (local_ierr == 1) then
                        STREAM_DEBUG_WRITE('-- Error: Could not find a file offset that exactly matches')
                        STREAM_DEBUG_WRITE('--        combination of reference time, current time, and interval.')
                    end if
                    ierr = MPAS_STREAM_MGR_ERROR
                    return
                end if

                if (temp_filename /= stream % filename) then

                    stream % filename = temp_filename

                    !
                    ! Close existing stream
                    !
                    call MPAS_closeStream(stream % stream, ierr=local_ierr)
                    if (local_ierr /= MPAS_STREAM_NOERR) then
                        ierr = MPAS_STREAM_MGR_ERROR
                        return
                    end if

                    !
                    ! Build new stream from pools of fields and attributes
                    !
                    call MPAS_createStream(stream % stream, stream % filename, MPAS_IO_PNETCDF, MPAS_IO_WRITE, &
                                           stream % maxRecords, ierr=local_ierr)
                    if (local_ierr /= MPAS_STREAM_NOERR) then
                        ierr = MPAS_STREAM_MGR_ERROR
                        return
                    end if

                    call build_stream(stream, MPAS_STREAM_OUTPUT, manager % allFields, timeLevel, mgLevel, local_ierr)
                    if (local_ierr /= MPAS_STREAM_NOERR) then
                        ierr = MPAS_STREAM_MGR_ERROR
                        return
                    end if
                    stream % timeLevel = timeLevel

                end if
            else
                stream % nRecords = stream % nRecords + 1
            end if
        end if

        if (timeLevel /= stream % timeLevel) then

            call update_stream(stream, manager % allFields, timeLevel, mgLevel, local_ierr)
            if (local_ierr /= MPAS_STREAM_NOERR) then
                ierr = MPAS_STREAM_MGR_ERROR
                return
            end if
            stream % timeLevel = timeLevel
        end if

        !
        ! For any connectivity arrays in the stream, temporarily convert local indices to global indices
        !
        call prewrite_reindex(manager % allFields, stream % field_pool)

        ! 
        ! Write the stream
        ! 
        call MPAS_writeStream(stream % stream, stream % nRecords, ierr=local_ierr)

        !
        ! Regardless of the error code from MPAS_writeStream, we need to reset global indices in the stream back to local indices
        !
        call postwrite_reindex(manager % allFields, stream % field_pool)

        if (local_ierr /= MPAS_STREAM_NOERR) then
            ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

    end subroutine write_stream


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

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: local_timeLevel
        integer :: local_mgLevel
        logical :: local_rightNow
        character (len=StrKIND) :: local_when
        integer :: local_whence
        integer :: local_ierr
        integer :: temp_ierr
        type (MPAS_Time_type) :: now_time 


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_read()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR
        if (present(actualWhen)) write(actualWhen,'(a)') '0000-01-01_00:00:00'

        !
        ! Use optional arguments or set defaults
        !
        if (present(timeLevel)) then
            local_timeLevel = timeLevel
        else
            local_timeLevel = 1
        end if

        if (present(mgLevel)) then
            local_mgLevel = mgLevel
        else
            local_mgLevel = 1
        end if

        if (present(rightNow)) then
            local_rightNow = rightNow
        else
            local_rightNow = .false.
        end if

        if (present(when)) then
            local_when = when
        else
            now_time = mpas_get_clock_time(manager % streamClock, MPAS_NOW, ierr=local_ierr)
            call mpas_get_time(now_time, dateTimeString=local_when)
        end if

        if (present(whence)) then
            local_whence = whence
        else
            local_whence = MPAS_STREAM_EXACT_TIME
        end if


        !
        ! If a stream is specified, we process just that stream; otherwise,
        !    process all streams
        !
        if (present(streamID)) then
            nullify(stream_cursor)
            if (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=ierr)) then
                STREAM_DEBUG_WRITE('-- Handling read of stream '//trim(stream_cursor % name))

                ! Verify that the stream is an input stream
                if (stream_cursor % direction /= MPAS_STREAM_INPUT .and. stream_cursor % direction /= MPAS_STREAM_INPUT_OUTPUT) then
                    STREAM_DEBUG_WRITE('-- Stream '//streamID//' is not an input stream.')
                    if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                    return 
                end if

                call read_stream(manager, stream_cursor, local_timeLevel, local_mgLevel, local_rightNow, local_when, local_whence, &
                                 actualWhen, local_ierr)
            else
                STREAM_DEBUG_WRITE('-- Stream '//streamID//' does not exist in call to MPAS_stream_mgr_read().')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return 
            end if
        else
            nullify(stream_cursor)
            stream_cursor => manager % streams % head
            do while (associated(stream_cursor))
                STREAM_DEBUG_WRITE('-- Handling read of stream '//trim(stream_cursor % name))

                ! Verify that the stream is an input stream
                if (stream_cursor % direction == MPAS_STREAM_INPUT .or. &
                    stream_cursor % direction /= MPAS_STREAM_INPUT_OUTPUT) then

                    !
                    ! What should be the meaning of actualWhen if we read multiple streams in this call?
                    !
                    call read_stream(manager, stream_cursor, local_timeLevel, local_mgLevel, local_rightNow, &
                                     local_when, local_whence, actualWhen, temp_ierr)
                    if (temp_ierr /= MPAS_STREAM_MGR_NOERR) then
                        local_ierr = MPAS_STREAM_MGR_ERROR
                    end if
                end if

                stream_cursor => stream_cursor % next
            end do
        end if

        if (present(ierr)) ierr = local_ierr

    end subroutine MPAS_stream_mgr_read


    !-----------------------------------------------------------------------
    !  routine read_stream
    !
    !> \brief Handle the reading of a stream pointed to by the stream list node
    !> \author Michael Duda
    !> \date   4 September 2014
    !> \details
    !>  Private subroutine to handle the details of actually reading a stream.
    !
    !-----------------------------------------------------------------------
    subroutine read_stream(manager, stream, timeLevel, mgLevel, forceReadNow, when, whence, actualWhen, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        type (MPAS_stream_list_type), intent(inout) :: stream
        integer, intent(in) :: timeLevel
        integer, intent(in) :: mgLevel
        logical, intent(in) :: forceReadNow
        character (len=*), intent(in) :: when
        integer, intent(in) :: whence
        character (len=*), intent(out), optional :: actualWhen
        integer, intent(out) :: ierr

        type (MPAS_stream_list_type), pointer :: alarm_cursor
        type (MPAS_Time_type) :: now_time
        type (MPAS_TimeInterval_type) :: temp_interval
        type (MPAS_TimeInterval_type) :: read_interval
        character (len=StrKIND) :: temp_filename
        logical :: ringing_alarm
        integer :: nAlarms
        integer :: local_ierr


        ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Check whether this stream is active
        !
        if (.not. stream % active_stream) then
            STREAM_DEBUG_WRITE('-- Stream '//trim(stream % name)//' is not currently active and will not be read.')
            return
        end if

        !
        ! Check whether all packages for this stream are inactive
        ! Note: if the stream has no packages, it is assumed to be active
        !
        if (.not. stream_active_pkg_check(stream)) then
            STREAM_DEBUG_WRITE('-- Stream '//trim(stream % name)//' has only inactive packages and will not be read.')
            return
        end if

        !
        ! Check whether any of the input alarms for the stream are ringing
        !
        ringing_alarm = .false.
        alarm_cursor => stream % alarmList_in % head
        do while (associated(alarm_cursor))
            if (mpas_is_alarm_ringing(manager % streamClock, alarm_cursor % name, ierr=local_ierr)) then
                ringing_alarm = .true.
                exit
            end if
            alarm_cursor => alarm_cursor % next
        end do

        if ((.not. ringing_alarm) .and. (.not. forceReadNow)) then
            return
        end if


        !
        ! Figure out the read interval for the stream, defined as the minimum
        !    period of any input alarm associated with the stream
        !
        if (.not. associated(stream % recordInterval)) then
            alarm_cursor => stream % alarmList_in % head
            nAlarms = 0
            do while (associated(alarm_cursor))
                temp_interval = mpas_alarm_interval(manager % streamClock, alarm_cursor % name, local_ierr)
                if (local_ierr == 0) then
                    if (nAlarms == 0) then
                        read_interval = temp_interval
                    else if (temp_interval < read_interval) then
                        read_interval = temp_interval
                    end if
                    nAlarms = nAlarms + 1
                end if
                alarm_cursor => alarm_cursor % next
            end do
        else
            read_interval = stream % recordInterval

            ! We should handle this stream below as if it were periodic and therefore had a periodic alarm
            nAlarms = 1       
        end if


        !
        ! If the stream is not valid, assume that we have not yet written this
        ! stream, in which case we create the stream from scratch
        !
        if (.not. stream % valid) then
            if (stream % maxRecords == 0 .or. nAlarms == 0) then
                call mpas_expand_string(when, stream % filename_template, stream % filename)
                stream % nRecords = 1
            else
                call mpas_set_time(now_time, dateTimeString=when, ierr=local_ierr)
                call build_filename(stream % referenceTime, now_time, whence, stream % maxRecords, read_interval, &
                                    stream % filename_template, stream % filename, stream % nRecords, local_ierr)
                if (local_ierr /= 0) then
                    if (local_ierr == 1) then
                        STREAM_DEBUG_WRITE('-- Error: Could not find a file offset that exactly matches')
                        STREAM_DEBUG_WRITE('--        combination of reference time, current time, and interval.')
                    end if
                    ierr = MPAS_STREAM_MGR_ERROR
                    return
                end if
            end if

            !
            ! Build stream from pools of fields and attributes
            !
            allocate(stream % stream)
            call MPAS_createStream(stream % stream, stream % filename, MPAS_IO_PNETCDF, MPAS_IO_READ, stream % maxRecords, &
                                   ierr=local_ierr)
            if (local_ierr /= MPAS_STREAM_NOERR) then
                ierr = MPAS_STREAM_MGR_ERROR
                return
            end if

            call build_stream(stream, MPAS_STREAM_INPUT, manager % allFields, timeLevel, mgLevel, local_ierr)
            if (local_ierr /= MPAS_STREAM_NOERR) then
                ierr = MPAS_STREAM_MGR_ERROR
                return
            end if
            stream % timeLevel = timeLevel

            stream % valid = .true.
        else
            if (stream % maxRecords /= 0 .and. nAlarms /= 0) then
                call mpas_set_time(now_time, dateTimeString=when, ierr=local_ierr)
                call build_filename(stream % referenceTime, now_time, whence, stream % maxRecords, read_interval, &
                                    stream % filename_template, temp_filename, stream % nRecords, local_ierr)
                if (local_ierr /= 0) then
                    if (local_ierr == 1) then
                        STREAM_DEBUG_WRITE('-- Error: Could not find a file offset that exactly matches')
                        STREAM_DEBUG_WRITE('--        combination of reference time, current time, and interval.')
                    end if
                    ierr = MPAS_STREAM_MGR_ERROR
                    return
                end if

                if (temp_filename /= stream % filename) then

                    stream % filename = temp_filename

                    !
                    ! Close existing stream
                    !
                    call MPAS_closeStream(stream % stream, ierr=local_ierr)
                    if (local_ierr /= MPAS_STREAM_NOERR) then
                        ierr = MPAS_STREAM_MGR_ERROR
                        return
                    end if

                    !
                    ! Build new stream from pools of fields and attributes
                    !
                    call MPAS_createStream(stream % stream, stream % filename, MPAS_IO_PNETCDF, MPAS_IO_READ, &
                                           stream % maxRecords, ierr=local_ierr)
                    if (local_ierr /= MPAS_STREAM_NOERR) then
                        ierr = MPAS_STREAM_MGR_ERROR
                        return
                    end if

                    call build_stream(stream, MPAS_STREAM_INPUT, manager % allFields, timeLevel, mgLevel, local_ierr)
                    if (local_ierr /= MPAS_STREAM_NOERR) then
                        ierr = MPAS_STREAM_MGR_ERROR
                        return
                    end if
                    stream % timeLevel = timeLevel

                end if
            else

                !
                ! TODO: Here, perhaps we should call build_filename just to ensure we read from the proper record,
                !       rather than assuming the records are at the same interval as calls to read the stream
                !

                stream % nRecords = stream % nRecords + 1
            end if
        end if

        if (timeLevel /= stream % timeLevel) then

            call update_stream(stream, manager % allFields, timeLevel, mgLevel, local_ierr)
            if (local_ierr /= MPAS_STREAM_NOERR) then
                ierr = MPAS_STREAM_MGR_ERROR
                return
            end if
            stream % timeLevel = timeLevel
        end if


        ! 
        ! Read the stream
        ! 
        call MPAS_readStream(stream % stream, stream % nRecords, ierr=local_ierr)
        if (local_ierr /= MPAS_STREAM_NOERR) then
            ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        if (present(actualWhen)) then
            call MPAS_streamTime(stream % stream, stream % nRecords, actualWhen, ierr=local_ierr)
            if (local_ierr /= MPAS_STREAM_NOERR) then
!
! TODO: Add debug prints for all error conditions
!
                ierr = MPAS_STREAM_MGR_ERROR
                return
            end if
        end if

        !
        ! Exchange halos for all decomposed fields in this stream
        !
        call exch_all_halos(manager % allFields, stream % field_pool, stream % timeLevel, local_ierr)

        !
        ! For any connectivity arrays in this stream, convert global indices to local indices
        !
        call postread_reindex(manager % allFields, stream % field_pool)

    end subroutine read_stream


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

        implicit none

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
    !  routine build_filename
    !
    !> \brief Construct the filename that contains a specific time in a stream
    !> \author Michael Duda
    !> \date   21 August 2014
    !> \details 
    !>  Given a filename template and the information necessary to determine the time
    !>  in the stream that matches a time available in any of the files associated with
    !>  the stream, returns a specific filename and a record offset within that file
    !>  matching that time.
    !>  
    !>  Implementation note: When two times in a stream are equidistant from the requested
    !>  time, the MPAS_STREAM_NEAREST option for 'whence' will choose the later time when
    !>  the requested time is before the reference time, and it will choose the earlier time
    !>  when the requested time is after the reference time. In other words, we round the
    !>  requested time towards the reference time.
    !>  
    !>  Return error codes:
    !>  0 no error
    !>  1 requested exact match for time but found remainder when calculating offset
    !-----------------------------------------------------------------------
    subroutine build_filename(ref_time, when, whence, maxRecords, period, filename_template, filename, offset, ierr)

        implicit none

        type (MPAS_Time_type), intent(in) :: ref_time
        type (MPAS_Time_type), intent(in) :: when
        integer, intent(in) :: whence
        integer, intent(in) :: maxRecords
        type (MPAS_TimeInterval_type), intent(in) :: period
        character(len=*), intent(in) :: filename_template
        character(len=*), intent(out) :: filename
        integer, intent(out) :: offset
        integer, intent(out) :: ierr

        character(len=StrKIND) :: when_string
        type (MPAS_Time_type) :: filetime
        type (MPAS_TimeInterval_type) :: zero
        type (MPAS_TimeInterval_type) :: intv, rem
        integer :: nrecs, nfiles
        logical :: in_future

!
! TODO: Use more descriptive error codes so that caller can tell why we couldn't build a filename
!

        ierr = 0

        if (when >= ref_time) then
            intv = when - ref_time
            in_future = .true.
        else 
            intv = ref_time - when
            in_future = .false.
        end if

        call interval_division(intv, period, nrecs, rem) 

        call mpas_set_timeInterval(zero, S=0)

        if (whence == MPAS_STREAM_EXACT_TIME .and. rem /= zero) then
            ierr = 1
            return

        else if (rem /= zero) then
            if ((whence == MPAS_STREAM_LATEST_STRICTLY_BEFORE .or. whence == MPAS_STREAM_LATEST_BEFORE) .and. .not. in_future) then
                nrecs = nrecs + 1
            else if ((whence == MPAS_STREAM_EARLIEST_STRICTLY_AFTER .or. whence == MPAS_STREAM_EARLIEST_AFTER) .and. in_future) then
                nrecs = nrecs + 1
            else if (whence == MPAS_STREAM_NEAREST) then
                if (rem > (period / 2)) then
                    nrecs = nrecs + 1
                end if
            end if

        else if (rem == zero) then
            if (in_future) then
                if (whence == MPAS_STREAM_EARLIEST_STRICTLY_AFTER) then
                    nrecs = nrecs + 1
                else if (whence == MPAS_STREAM_LATEST_STRICTLY_BEFORE) then
                    nrecs = nrecs - 1
                end if
            else
                if (whence == MPAS_STREAM_EARLIEST_STRICTLY_AFTER) then
                    nrecs = nrecs - 1
                else if (whence == MPAS_STREAM_LATEST_STRICTLY_BEFORE) then
                    nrecs = nrecs + 1
                end if
            end if
        end if 

        nfiles = nrecs / maxRecords
        offset = nrecs - (nfiles * maxRecords)

        if (in_future) then
            filetime = ref_time + period * (maxRecords * nfiles)
        else
            if (offset /= 0) then
                nfiles = nfiles + 1
                offset = maxRecords - offset
            end if
            filetime = ref_time - period * (maxRecords * nfiles)
        end if

        ! Switch from a 0-based to a 1-based offset
        offset = offset + 1

        call mpas_get_time(filetime, dateTimeString=when_string)

        call mpas_expand_string(when_string, filename_template, filename)

    end subroutine build_filename


    !-----------------------------------------------------------------------
    !  routine build_stream
    !
    !> \brief This is a utility routine to build a stream type from a pool representing a stream.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   07/23/2014
    !> \details 
    !>  This routine will take as input a pool representing a stream.
    !>  It will then generate a stream type based on this pool, and return that.
    !-----------------------------------------------------------------------
    subroutine build_stream(stream, direction, allFields, timeLevelIn, mgLevelIn, ierr) !{{{

        implicit none

        type (MPAS_stream_list_type), intent(inout) :: stream
        integer, intent(in) :: direction
        type (MPAS_Pool_type), intent(in) :: allFields
        integer, intent(in) :: timeLevelIn
        integer, intent(in) :: mgLevelIn
        integer, intent(out) :: ierr

        type (MPAS_Pool_iterator_type) :: itr
        type (mpas_pool_field_info_type) :: info
        integer :: timeLevel

        type (field5DReal), pointer :: real5d
        type (field4DReal), pointer :: real4d
        type (field3DReal), pointer :: real3d
        type (field2DReal), pointer :: real2d
        type (field1DReal), pointer :: real1d
        type (field0DReal), pointer :: real0d

        type (field3DInteger), pointer :: int3d
        type (field2DInteger), pointer :: int2d
        type (field1DInteger), pointer :: int1d
        type (field0DInteger), pointer :: int0d

        type (field1DChar), pointer :: char1d
        type (field0DChar), pointer :: char0d

        integer, pointer :: intAtt
        logical, pointer :: logAtt
        character (len=StrKIND), pointer :: charAtt
        real (kind=RKIND), pointer :: realAtt

        integer :: local_ierr

        integer, parameter :: idLength = 10
        character (len=idLength) :: file_id


        if (direction == MPAS_STREAM_OUTPUT) then

            !
            ! Write attributes to stream
            !
            call mpas_pool_begin_iteration(stream % att_pool)
            do while (mpas_pool_get_next_member(stream % att_pool, itr))
                if ( itr % memberType == MPAS_POOL_CONFIG) then
                    if ( itr % dataType == MPAS_POOL_REAL ) then
                        call mpas_pool_get_config(stream % att_pool, itr % memberName, realAtt)
                        call mpas_writeStreamAtt(stream % stream, itr % memberName, realAtt, local_ierr)
                    else if ( itr % dataType == MPAS_POOL_INTEGER ) then
                        call mpas_pool_get_config(stream % att_pool, itr % memberName, intAtt)
                        call mpas_writeStreamAtt(stream % stream, itr % memberName, intAtt, local_ierr)
                    else if ( itr % dataType == MPAS_POOL_CHARACTER ) then
                        call mpas_pool_get_config(stream % att_pool, itr % memberName, charAtt)
                        call mpas_writeStreamAtt(stream % stream, itr % memberName, charAtt, local_ierr)
                    else if ( itr % dataType == MPAS_POOL_LOGICAL ) then
                        call mpas_pool_get_config(stream % att_pool, itr % memberName, logAtt)
                        if (logAtt) then
                            call mpas_writeStreamAtt(stream % stream, itr % memberName, 'YES', local_ierr)
                        else
                            call mpas_writeStreamAtt(stream % stream, itr % memberName, 'NO', local_ierr)
                        end if
                    end if

                    if (local_ierr /= MPAS_STREAM_NOERR) then
                        ierr = MPAS_STREAM_MGR_ERROR
                        return
                    end if
                end if
            end do

            !
            ! Generate file_id and write to stream
            !
            call gen_random(idLength, file_id)
write(0,*) 'RANDOM ID '//file_id
            call mpas_writeStreamAtt(stream % stream, 'file_id', file_id, local_ierr)
            if (local_ierr /= MPAS_STREAM_NOERR) then
                ierr = MPAS_STREAM_MGR_ERROR
                return
            end if

        end if


        ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_begin_iteration(stream % field_pool)

        do while ( mpas_pool_get_next_member(stream % field_pool, itr) )

            if (itr % memberType == MPAS_POOL_CONFIG) then

                ! To avoid accidentally matching in case statements below...
                info % fieldType = -1

                call mpas_pool_get_field_info(allFields, itr % memberName, info)

                ! Set time level to read
                if (info % nTimeLevels >= timeLevelIn) then
                    timeLevel = timeLevelIn
                else
                    timeLevel = 1
                end if

                select case (info % fieldType)
                    case (MPAS_POOL_REAL)
                        select case (info % nDims)
                            case (0)
                                call mpas_pool_get_field(allFields, itr % memberName, real0d, timeLevel)
                                call MPAS_streamAddField(stream % stream, real0d)
                            case (1)
                                call mpas_pool_get_field(allFields, itr % memberName, real1d, timeLevel)
                                call MPAS_streamAddField(stream % stream, real1d)
                            case (2)
                                call mpas_pool_get_field(allFields, itr % memberName, real2d, timeLevel)
                                call MPAS_streamAddField(stream % stream, real2d)
                            case (3)
                                call mpas_pool_get_field(allFields, itr % memberName, real3d, timeLevel)
                                call MPAS_streamAddField(stream % stream, real3d)
                            case (4)
                                call mpas_pool_get_field(allFields, itr % memberName, real4d, timeLevel)
                                call MPAS_streamAddField(stream % stream, real4d)
                            case (5)
                                call mpas_pool_get_field(allFields, itr % memberName, real5d, timeLevel)
                                call MPAS_streamAddField(stream % stream, real5d)
                        end select
                    case (MPAS_POOL_INTEGER)
                        select case (info % nDims)
                            case (0)
                                call mpas_pool_get_field(allFields, itr % memberName, int0d, timeLevel)
                                call MPAS_streamAddField(stream % stream, int0d)
                            case (1)
                                call mpas_pool_get_field(allFields, itr % memberName, int1d, timeLevel)
                                call MPAS_streamAddField(stream % stream, int1d)
                            case (2)
                                call mpas_pool_get_field(allFields, itr % memberName, int2d, timeLevel)
                                call MPAS_streamAddField(stream % stream, int2d)
                            case (3)
                                call mpas_pool_get_field(allFields, itr % memberName, int3d, timeLevel)
                                call MPAS_streamAddField(stream % stream, int3d)
                        end select
                    case (MPAS_POOL_CHARACTER)
                        select case (info % nDims)
                            case (0)
                                call mpas_pool_get_field(allFields, itr % memberName, char0d, timeLevel)
                                call MPAS_streamAddField(stream % stream, char0d)
                            case (1)
!                                call mpas_pool_get_field(allFields, itr % memberName, char1d, timeLevel)
!                                call MPAS_streamAddField(stream % stream, char1d)
                                 write(stderrUnit,*) 'Error: In build_stream, unsupported type field1DChar.'
                        end select
                end select

            end if
        end do

    end subroutine build_stream !}}}


    !-----------------------------------------------------------------------
    !  routine update_stream
    !
    !> \brief Updates the time level for fields in a stream
    !> \author Michael Duda, Doug Jacobsen
    !> \date   07/23/2014
    !> \details 
    !>  For an existing stream, updates the time levels for all fields in 
    !>  the stream so that subsequent reads/writes of the stream will read
    !>  from / write to the specified time level.
    !-----------------------------------------------------------------------
    subroutine update_stream(stream, allFields, timeLevelIn, mgLevelIn, ierr) !{{{

        implicit none

        type (MPAS_stream_list_type), intent(inout) :: stream
        type (MPAS_Pool_type), intent(in) :: allFields
        integer, intent(in) :: timeLevelIn
        integer, intent(in) :: mgLevelIn
        integer, intent(out) :: ierr

        type (MPAS_Pool_iterator_type) :: itr
        type (mpas_pool_field_info_type) :: info
        integer :: timeLevel

        type (field5DReal), pointer :: real5d
        type (field4DReal), pointer :: real4d
        type (field3DReal), pointer :: real3d
        type (field2DReal), pointer :: real2d
        type (field1DReal), pointer :: real1d
        type (field0DReal), pointer :: real0d

        type (field3DInteger), pointer :: int3d
        type (field2DInteger), pointer :: int2d
        type (field1DInteger), pointer :: int1d
        type (field0DInteger), pointer :: int0d

        type (field1DChar), pointer :: char1d
        type (field0DChar), pointer :: char0d


        ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_begin_iteration(stream % field_pool)

        do while ( mpas_pool_get_next_member(stream % field_pool, itr) )

            if (itr % memberType == MPAS_POOL_CONFIG) then

                ! To avoid accidentally matching in case statements below...
                info % fieldType = -1

                call mpas_pool_get_field_info(allFields, itr % memberName, info)

                ! Set time level to read
                if (info % nTimeLevels >= timeLevelIn) then
                    timeLevel = timeLevelIn
                else
                    timeLevel = 1
                end if

                select case (info % fieldType)
                    case (MPAS_POOL_REAL)
                        select case (info % nDims)
                            case (0)
                                call mpas_pool_get_field(allFields, itr % memberName, real0d, timeLevel)
                                call MPAS_streamUpdateField(stream % stream, real0d)
                            case (1)
                                call mpas_pool_get_field(allFields, itr % memberName, real1d, timeLevel)
                                call MPAS_streamUpdateField(stream % stream, real1d)
                            case (2)
                                call mpas_pool_get_field(allFields, itr % memberName, real2d, timeLevel)
                                call MPAS_streamUpdateField(stream % stream, real2d)
                            case (3)
                                call mpas_pool_get_field(allFields, itr % memberName, real3d, timeLevel)
                                call MPAS_streamUpdateField(stream % stream, real3d)
                            case (4)
                                call mpas_pool_get_field(allFields, itr % memberName, real4d, timeLevel)
                                call MPAS_streamUpdateField(stream % stream, real4d)
                            case (5)
                                call mpas_pool_get_field(allFields, itr % memberName, real5d, timeLevel)
                                call MPAS_streamUpdateField(stream % stream, real5d)
                        end select
                    case (MPAS_POOL_INTEGER)
                        select case (info % nDims)
                            case (0)
                                call mpas_pool_get_field(allFields, itr % memberName, int0d, timeLevel)
                                call MPAS_streamUpdateField(stream % stream, int0d)
                            case (1)
                                call mpas_pool_get_field(allFields, itr % memberName, int1d, timeLevel)
                                call MPAS_streamUpdateField(stream % stream, int1d)
                            case (2)
                                call mpas_pool_get_field(allFields, itr % memberName, int2d, timeLevel)
                                call MPAS_streamUpdateField(stream % stream, int2d)
                            case (3)
                                call mpas_pool_get_field(allFields, itr % memberName, int3d, timeLevel)
                                call MPAS_streamUpdateField(stream % stream, int3d)
                        end select
                    case (MPAS_POOL_CHARACTER)
                        select case (info % nDims)
                            case (0)
                                call mpas_pool_get_field(allFields, itr % memberName, char0d, timeLevel)
                                call MPAS_streamUpdateField(stream % stream, char0d)
                            case (1)
!                                call mpas_pool_get_field(allFields, itr % memberName, char1d, timeLevel)
!                                call MPAS_streamUpdateField(stream % stream, char1d)
                                 write(stderrUnit,*) 'Error: In update_stream, unsupported type field1DChar.'
                        end select
                end select

            end if
        end do

    end subroutine update_stream !}}}


    !-----------------------------------------------------------------------
    !  routine stream_active_pkg_check
    !
    !> \brief Checks whether a stream has any active packages (or none at all)
    !> \author Michael Duda
    !> \date   23 September 2014
    !> \details 
    !>  This function determines whether a stream has any active packages
    !>  associated with it. If the stream has at least one active package,
    !>  or no packages at all, the function returns true; else, if all packages
    !>  associated with the package are inactive, the function returns false.
    !
    !-----------------------------------------------------------------------
    logical function stream_active_pkg_check(stream)

        implicit none

        type (MPAS_stream_list_type), intent(inout) :: stream

        type (MPAS_Pool_iterator_type) :: itr
        logical, pointer :: pkg_val
        integer :: npkgs


        stream_active_pkg_check = .false.
        npkgs = 0
        call mpas_pool_begin_iteration(stream % pkg_pool)

        do while ( mpas_pool_get_next_member(stream % pkg_pool, itr) )
            if (itr % memberType == MPAS_POOL_PACKAGE) then
                nullify(pkg_val)
                call mpas_pool_get_package(stream % pkg_pool, trim(itr % memberName), pkg_val)
                if (associated(pkg_val)) then
                    npkgs = npkgs + 1
                    stream_active_pkg_check = stream_active_pkg_check .or. pkg_val
                end if
            else
                ! This is unexpected...
                STREAM_DEBUG_WRITE('... found non-package '//trim(itr % memberName)//' in package pool for stream '//trim(stream % name))
            end if
        end do

        if (npkgs == 0) then
            stream_active_pkg_check = .true.
        end if

    end function stream_active_pkg_check


    !-----------------------------------------------------------------------
    !  routine exch_all_halos
    !
    !> \brief Exchange halos of all fields in stream
    !> \author Doug Jacobsen, Michael Duda
    !> \date   09/12/2014
    !> \details
    !>  This routine performs a halo exchange of each decomposed field within a stream.
    !
    !-----------------------------------------------------------------------
    subroutine exch_all_halos(allFields, streamFields, timeLevel, ierr) !{{{

        implicit none

        type (mpas_pool_type), pointer :: allFields
        type (mpas_pool_type), pointer :: streamFields
        integer, intent(in) :: timeLevel
        integer, intent(out) :: ierr

        type (mpas_pool_iterator_type) :: fieldItr
        type (mpas_pool_field_info_type) :: fieldInfo

        type (field1DReal), pointer :: real1DField
        type (field2DReal), pointer :: real2DField
        type (field3DReal), pointer :: real3DField
        type (field4DReal), pointer :: real4DField
        type (field5DReal), pointer :: real5DField
        type (field1DInteger), pointer :: int1DField
        type (field2DInteger), pointer :: int2DField
        type (field3DInteger), pointer :: int3DField


        ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_begin_iteration(streamFields)

        do while ( mpas_pool_get_next_member(streamFields, fieldItr) )

            ! Note: in a stream's field_pool, the names of fields are stored as configs
            if ( fieldItr % memberType == MPAS_POOL_CONFIG ) then

                call mpas_pool_get_field_info(allFields, fieldItr % memberName, fieldInfo)

                if ( fieldInfo % nDims == 1) then
                    if ( fieldInfo % fieldType == MPAS_POOL_REAL ) then
                        if ( timeLevel <= fieldInfo % nTimeLevels ) then
                            call mpas_pool_get_field(allFields, fieldItr % memberName, real1DField, timeLevel)
                        else
                            call mpas_pool_get_field(allFields, fieldItr % memberName, real1DField, 1)
                        end if

                        if ( is_decomposed_dim(real1DField % dimNames(1))) then
                            STREAM_DEBUG_WRITE(' -- Exchange halo for '//trim(fieldItr % memberName))
                            call mpas_dmpar_exch_halo_field(real1DField)
                        end if
                    else if ( fieldInfo % fieldType == MPAS_POOL_INTEGER ) then
                        if ( timeLevel <= fieldInfo % nTimeLevels ) then
                          call mpas_pool_get_field(allFields, fieldItr % memberName, int1DField, timeLevel)
                        else
                          call mpas_pool_get_field(allFields, fieldItr % memberName, int1DField, 1)
                        end if
                        if ( is_decomposed_dim(int1DField % dimNames(1))) then
                            STREAM_DEBUG_WRITE(' -- Exchange halo for '//trim(fieldItr % memberName))
                            call mpas_dmpar_exch_halo_field(int1DField)
                        end if
                    end if

                else if ( fieldInfo % nDims == 2) then
                    if ( fieldInfo % fieldType == MPAS_POOL_REAL ) then
                        if ( timeLevel <= fieldInfo % nTimeLevels ) then
                            call mpas_pool_get_field(allFields, fieldItr % memberName, real2DField, timeLevel)
                        else
                            call mpas_pool_get_field(allFields, fieldItr % memberName, real2DField, 1)
                        end if
                        if ( is_decomposed_dim(real2DField % dimNames(2))) then
                            STREAM_DEBUG_WRITE(' -- Exchange halo for '//trim(fieldItr % memberName))
                            call mpas_dmpar_exch_halo_field(real2DField)
                        end if
                    else if ( fieldInfo % fieldType == MPAS_POOL_INTEGER ) then
                        if ( timeLevel <= fieldInfo % nTimeLevels ) then
                            call mpas_pool_get_field(allFields, fieldItr % memberName, int2DField, timeLevel)
                        else
                            call mpas_pool_get_field(allFields, fieldItr % memberName, int2DField, 1)
                        end if
                        if ( is_decomposed_dim(int2DField % dimNames(2))) then
                            STREAM_DEBUG_WRITE(' -- Exchange halo for '//trim(fieldItr % memberName))
                            call mpas_dmpar_exch_halo_field(int2DField)
                        end if
                    end if

                else if ( fieldInfo % nDims == 3) then
                    if ( fieldInfo % fieldType == MPAS_POOL_REAL ) then
                        if ( timeLevel <= fieldInfo % nTimeLevels ) then
                            call mpas_pool_get_field(allFields, fieldItr % memberName, real3DField, timeLevel)
                        else
                            call mpas_pool_get_field(allFields, fieldItr % memberName, real3DField, 1)
                        end if
                        if ( is_decomposed_dim(real3DField % dimNames(3))) then
                            STREAM_DEBUG_WRITE(' -- Exchange halo for '//trim(fieldItr % memberName))
                            call mpas_dmpar_exch_halo_field(real3DField)
                        end if
                    else if ( fieldInfo % fieldType == MPAS_POOL_INTEGER ) then
                        if ( timeLevel <= fieldInfo % nTimeLevels ) then
                            call mpas_pool_get_field(allFields, fieldItr % memberName, int3DField, timeLevel)
                        else
                            call mpas_pool_get_field(allFields, fieldItr % memberName, int3DField, 1)
                        end if
                        if ( is_decomposed_dim(int3DField % dimNames(3))) then
                            STREAM_DEBUG_WRITE(' -- Exchange halo for '//trim(fieldItr % memberName))
                            call mpas_dmpar_exch_halo_field(int3DField)
                        end if
                    end if

                else if ( fieldInfo % nDims == 4) then
                    if ( fieldInfo % fieldType == MPAS_POOL_REAL ) then
                        if ( timeLevel <= fieldInfo % nTimeLevels ) then
                            call mpas_pool_get_field(allFields, fieldItr % memberName, real4DField, timeLevel)
                        else
                            call mpas_pool_get_field(allFields, fieldItr % memberName, real4DField, 1)
                        end if
                        if ( is_decomposed_dim(real4DField % dimNames(4))) then
                            STREAM_DEBUG_WRITE(' -- Exchange halo for '//trim(fieldItr % memberName))
                            call mpas_dmpar_exch_halo_field(real4DField)
                        end if
                    end if

                else if ( fieldInfo % nDims == 5) then
                    if ( fieldInfo % fieldType == MPAS_POOL_REAL ) then
                        if ( timeLevel <= fieldInfo % nTimeLevels ) then
                            call mpas_pool_get_field(allFields, fieldItr % memberName, real5DField, timeLevel)
                        else
                            call mpas_pool_get_field(allFields, fieldItr % memberName, real5DField, 1)
                        end if
                        if ( is_decomposed_dim(real5DField % dimNames(5))) then
                            STREAM_DEBUG_WRITE(' -- Exchange halo for '//trim(fieldItr % memberName))
                            call mpas_dmpar_exch_halo_field(real5DField)
                        end if
                    end if
                end if

            end if

        end do

    end subroutine exch_all_halos !}}}


    !-----------------------------------------------------------------------
    !  routine is_decomposed_dim
    !
    !> \brief Determines whether a dimension represents a decomposed dimension or not
    !> \author Michael Duda
    !> \date   24 September 2014
    !> \details 
    !>  This function determines whether the name of the input argument is 
    !>  a decompsed dimension or not. Currently in MPAS, the only decomposed
    !>  dimensions are:
    !>      nCells
    !>      nEdges
    !>      nVertices
    !
    !-----------------------------------------------------------------------
    logical function is_decomposed_dim(dimName)

        implicit none

        character(len=*), intent(in) :: dimName

        if (trim(dimName) == 'nCells' .or. &
            trim(dimName) == 'nEdges' .or. &
            trim(dimName) == 'nVertices') then

            is_decomposed_dim = .true.

        else

            is_decomposed_dim = .false.

        end if

    end function is_decomposed_dim


    !-----------------------------------------------------------------------
    !  routine prewrite_reindex
    !
    !> \brief Reindex connectivity fields from local to global index space.
    !> \author Doug Jacobsen, Michael Duda
    !> \date   24 September 2014
    !> \details 
    !>  For any connectivity fields contained in the stream to be written, 
    !>  whose fields include those in the streamFields pool, save the locally
    !>  indexed fields in module variables *_save, and allocate new arrays for
    !>  the fields, which are set to contain global indices.
    !>  This routine should be called immediately before a write of a stream.
    !
    !-----------------------------------------------------------------------
    subroutine prewrite_reindex(allFields, streamFields)

        implicit none

        type (mpas_pool_type), pointer :: allFields
        type (mpas_pool_type), pointer :: streamFields

        type (mpas_pool_iterator_type) :: fieldItr
        type (mpas_pool_field_info_type) :: fieldInfo

        integer, pointer :: nCells, nEdges, nVertices, vertexDegree
        integer, pointer :: maxEdges, maxEdges2, nEdgesSolve, nCellsSolve, nVerticesSolve

        type (field1dInteger), pointer :: nEdgesOnCell, nEdgesOnEdge, indexToCellID, indexToEdgeID, indexToVertexID

        type (field2dInteger), pointer :: cellsOnCell_ptr, edgesOnCell_ptr, verticesOnCell_ptr, &
                                 cellsOnEdge_ptr, verticesOnEdge_ptr, edgesOnEdge_ptr, &
                                 cellsOnVertex_ptr, edgesOnVertex_ptr

        type (field2dInteger), pointer :: cellsOnCell, edgesOnCell, verticesOnCell, &
                                 cellsOnEdge, verticesOnEdge, edgesOnEdge, &
                                 cellsOnVertex, edgesOnVertex

        logical :: handle_cellsOnCell, handle_edgesOnCell, handle_verticesOnCell, handle_cellsOnEdge, handle_verticesOnEdge, &
                   handle_edgesOnEdge, handle_cellsOnVertex, handle_edgesOnVertex

        integer :: i, j


        nullify(cellsOnCell_save)
        nullify(edgesOnCell_save)
        nullify(verticesOnCell_save)
        nullify(cellsOnEdge_save)
        nullify(verticesOnEdge_save)
        nullify(edgesOnEdge_save)
        nullify(cellsOnVertex_save)
        nullify(edgesOnVertex_save)

        nullify(cellsOnCell)
        nullify(edgesOnCell)
        nullify(verticesOnCell)
        nullify(cellsOnEdge)
        nullify(verticesOnEdge)
        nullify(edgesOnEdge)
        nullify(cellsOnVertex)
        nullify(edgesOnVertex)

        !
        ! Determine which connectivity fields exist in this stream
        !
        call mpas_pool_begin_iteration(streamFields)
        do while ( mpas_pool_get_next_member(streamFields, fieldItr) )

            ! Note: in a stream's field_pool, the names of fields are stored as configs
            if ( fieldItr % memberType == MPAS_POOL_CONFIG ) then
                call mpas_pool_get_field_info(allFields, fieldItr % memberName, fieldInfo)

                if (trim(fieldItr % memberName) == 'cellsOnCell') then
                    allocate(cellsOnCell_save)
                    cellsOnCell_ptr => cellsOnCell_save
                    call mpas_pool_get_field(allFields, 'cellsOnCell', cellsOnCell)
                else if (trim(fieldItr % memberName) == 'edgesOnCell') then
                    allocate(edgesOnCell_save)
                    edgesOnCell_ptr => edgesOnCell_save
                    call mpas_pool_get_field(allFields, 'edgesOnCell', edgesOnCell)
                else if (trim(fieldItr % memberName) == 'verticesOnCell') then
                    allocate(verticesOnCell_save)
                    verticesOnCell_ptr => verticesOnCell_save
                    call mpas_pool_get_field(allFields, 'verticesOnCell', verticesOnCell)
                else if (trim(fieldItr % memberName) == 'cellsOnEdge') then
                    allocate(cellsOnEdge_save)
                    cellsOnEdge_ptr => cellsOnEdge_save
                    call mpas_pool_get_field(allFields, 'cellsOnEdge', cellsOnEdge)
                else if (trim(fieldItr % memberName) == 'verticesOnEdge') then
                    allocate(verticesOnEdge_save)
                    verticesOnEdge_ptr => verticesOnEdge_save
                    call mpas_pool_get_field(allFields, 'verticesOnEdge', verticesOnEdge)
                else if (trim(fieldItr % memberName) == 'edgesOnEdge') then
                    allocate(edgesOnEdge_save)
                    edgesOnEdge_ptr => edgesOnEdge_save
                    call mpas_pool_get_field(allFields, 'edgesOnEdge', edgesOnEdge)
                else if (trim(fieldItr % memberName) == 'cellsOnVertex') then
                    allocate(cellsOnVertex_save)
                    cellsOnVertex_ptr => cellsOnVertex_save
                    call mpas_pool_get_field(allFields, 'cellsOnVertex', cellsOnVertex)
                else if (trim(fieldItr % memberName) == 'edgesOnVertex') then
                    allocate(edgesOnVertex_save)
                    edgesOnVertex_ptr => edgesOnVertex_save
                    call mpas_pool_get_field(allFields, 'edgesOnVertex', edgesOnVertex)
                end if
            end if

        end do

        !
        ! Reindex connectivity from local to global index space
        !
        call mpas_pool_get_field(allFields, 'nEdgesOnCell', nEdgesOnCell)
        call mpas_pool_get_field(allFields, 'nEdgesOnEdge', nEdgesOnEdge)
        call mpas_pool_get_field(allFields, 'indexToCellID', indexToCellID)
        call mpas_pool_get_field(allFields, 'indexToEdgeID', indexToEdgeID)
        call mpas_pool_get_field(allFields, 'indexToVertexID', indexToVertexID)

        do while (associated(indexToCellID))

            call mpas_pool_get_dimension(indexToCellID % block % dimensions, 'nCells', nCells)
            call mpas_pool_get_dimension(indexToCellID % block % dimensions, 'nEdges', nEdges)
            call mpas_pool_get_dimension(indexToCellID % block % dimensions, 'nVertices', nVertices)
            call mpas_pool_get_dimension(indexToCellID % block % dimensions, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_dimension(indexToCellID % block % dimensions, 'nEdgesSolve', nEdgesSolve)
            call mpas_pool_get_dimension(indexToCellID % block % dimensions, 'nVerticesSolve', nVerticesSolve)
            call mpas_pool_get_dimension(indexToCellID % block % dimensions, 'maxEdges', maxEdges)
            call mpas_pool_get_dimension(indexToCellID % block % dimensions, 'maxEdges2', maxEdges2)
            call mpas_pool_get_dimension(indexToCellID % block % dimensions, 'vertexDegree', vertexDegree)

            if (associated(cellsOnCell)) then
                nullify(cellsOnCell_ptr % ioinfo)
                cellsOnCell_ptr % array => cellsOnCell % array
                allocate(cellsOnCell % array(maxEdges, nCells+1))

                do i = 1, nCellsSolve
                    do j = 1, nEdgesOnCell % array(i)
                        cellsOnCell % array(j,i) = indexToCellID % array(cellsOnCell_ptr % array(j,i))
                    end do

                    cellsOnCell % array(nEdgesOnCell%array(i)+1:maxEdges,i) = nCells+1
                end do

                cellsOnCell => cellsOnCell % next
                if (associated(cellsOnCell)) then
                    allocate(cellsOnCell_ptr % next)
                    cellsOnCell_ptr => cellsOnCell_ptr % next
                end if
                nullify(cellsOnCell_ptr % next)
            end if

            if (associated(edgesOnCell)) then
                nullify(edgesOnCell_ptr % ioinfo)
                edgesOnCell_ptr % array => edgesOnCell % array
                allocate(edgesOnCell % array(maxEdges, nCells+1))

                do i = 1, nCellsSolve
                    do j = 1, nEdgesOnCell % array(i)
                        edgesOnCell % array(j,i) = indexToEdgeID % array(edgesOnCell_ptr % array(j,i))
                    end do

                    edgesOnCell % array(nEdgesOnCell%array(i)+1:maxEdges,i) = nEdges+1
                end do

                edgesOnCell => edgesOnCell % next
                if (associated(edgesOnCell)) then
                    allocate(edgesOnCell_ptr % next)
                    edgesOnCell_ptr => edgesOnCell_ptr % next
                end if
                nullify(edgesOnCell_ptr % next)
            end if

            if (associated(verticesOnCell)) then
                nullify(verticesOnCell_ptr % ioinfo)
                verticesOnCell_ptr % array => verticesOnCell % array
                allocate(verticesOnCell % array(maxEdges, nCells+1))

                do i = 1, nCellsSolve
                    do j = 1, nEdgesOnCell % array(i)
                        verticesOnCell % array(j,i) = indexToVertexID % array(verticesOnCell_ptr % array(j,i))
                    end do

                    verticesOnCell % array(nEdgesOnCell%array(i)+1:maxEdges,i) = nVertices+1
                end do

                verticesOnCell => verticesOnCell % next
                if (associated(verticesOnCell)) then
                    allocate(verticesOnCell_ptr % next)
                    verticesOnCell_ptr => verticesOnCell_ptr % next
                end if
                nullify(verticesOnCell_ptr % next)
            end if

            if (associated(cellsOnEdge)) then
                nullify(cellsOnEdge_ptr % ioinfo)
                cellsOnEdge_ptr % array => cellsOnEdge % array
                allocate(cellsOnEdge % array(2, nEdges+1))

                do i = 1, nEdgesSolve
                    cellsOnEdge % array(1,i) = indexToCellID % array(cellsOnEdge_ptr % array(1,i))
                    cellsOnEdge % array(2,i) = indexToCellID % array(cellsOnEdge_ptr % array(2,i))
                end do

                cellsOnEdge => cellsOnEdge % next
                if (associated(cellsOnEdge)) then
                    allocate(cellsOnEdge_ptr % next)
                    cellsOnEdge_ptr => cellsOnEdge_ptr % next
                end if
                nullify(cellsOnEdge_ptr % next)
            end if

            if (associated(verticesOnEdge)) then
                nullify(verticesOnEdge_ptr % ioinfo)
                verticesOnEdge_ptr % array => verticesOnEdge % array
                allocate(verticesOnEdge % array(2, nEdges+1))

                do i = 1, nEdgesSolve
                    verticesOnEdge % array(1,i) = indexToVertexID % array(verticesOnEdge_ptr % array(1,i))
                    verticesOnEdge % array(2,i) = indexToVertexID % array(verticesOnEdge_ptr % array(2,i))
                end do

                verticesOnEdge => verticesOnEdge % next
                if (associated(verticesOnEdge)) then
                    allocate(verticesOnEdge_ptr % next)
                    verticesOnEdge_ptr => verticesOnEdge_ptr % next
                end if
                nullify(verticesOnEdge_ptr % next)
            end if

            if (associated(edgesOnEdge)) then
                nullify(edgesOnEdge_ptr % ioinfo)
                edgesOnEdge_ptr % array => edgesOnEdge % array
                allocate(edgesOnEdge % array(maxEdges2, nEdges+1))

                do i = 1, nEdgesSolve
                    do j = 1, nEdgesOnEdge % array(i)
                        edgesOnEdge % array(j,i) = indexToEdgeID % array(edgesOnEdge_ptr % array(j,i))
                    end do

                    edgesOnEdge % array(nEdgesOnEdge%array(i)+1:maxEdges2,i) = nEdges+1
                end do

                edgesOnEdge => edgesOnEdge % next
                if (associated(edgesOnEdge)) then
                    allocate(edgesOnEdge_ptr % next)
                    edgesOnEdge_ptr => edgesOnEdge_ptr % next
                end if
                nullify(edgesOnEdge_ptr % next)
            end if

            if (associated(cellsOnVertex)) then
                nullify(cellsOnVertex_ptr % ioinfo)
                cellsOnVertex_ptr % array => cellsOnVertex % array
                allocate(cellsOnVertex % array(vertexDegree, nVertices+1))

                do i = 1, nVerticesSolve
                    do j = 1, vertexDegree
                        cellsOnVertex % array(j,i) = indexToCellID % array(cellsOnVertex_ptr % array(j,i))
                    end do
                end do

                cellsOnVertex => cellsOnVertex % next
                if (associated(cellsOnVertex)) then
                    allocate(cellsOnVertex_ptr % next)
                    cellsOnVertex_ptr => cellsOnVertex_ptr % next
                end if
                nullify(cellsOnVertex_ptr % next)
            end if

            if (associated(edgesOnVertex)) then
                nullify(edgesOnVertex_ptr % ioinfo)
                edgesOnVertex_ptr % array => edgesOnVertex % array
                allocate(edgesOnVertex % array(vertexDegree, nVertices+1))

                do i = 1, nVerticesSolve
                    do j = 1, vertexDegree
                        edgesOnVertex % array(j,i) = indexToEdgeID % array(edgesOnVertex_ptr % array(j,i))
                    end do
                end do

                edgesOnVertex => edgesOnVertex % next
                if (associated(edgesOnVertex)) then
                    allocate(edgesOnVertex_ptr % next)
                    edgesOnVertex_ptr => edgesOnVertex_ptr % next
                end if
                nullify(edgesOnVertex_ptr % next)
            end if

            nEdgesOnCell => nEdgesOnCell % next
            nEdgesOnEdge => nEdgesOnEdge % next
            indexToCellID => indexToCellID % next
            indexToEdgeID => indexToEdgeID % next
            indexToVertexID => indexToVertexID % next

        end do

    end subroutine prewrite_reindex


    !-----------------------------------------------------------------------
    !  routine postwrite_reindex
    !
    !> \brief Reindex connectivity fields from global to local index space.
    !> \author Doug Jacobsen, Michael Duda
    !> \date   24 September 2014
    !> \details 
    !>  For any connectivity fields contained in the stream to be written, 
    !>  whose fields include those in the streamFields pool, restore the locally
    !>  indexed fields from module variables *_save.
    !>  This routine should be called immediately after a write of a stream.
    !> 
    !>  NB: Even if the write of a stream fails, it is important to stil call
    !>      this routine to reset the connectivity fields to contain local indices.
    !
    !-----------------------------------------------------------------------
    subroutine postwrite_reindex(allFields, streamFields)

        implicit none

        type (mpas_pool_type), pointer :: allFields
        type (mpas_pool_type), pointer :: streamFields

        type (field1dInteger), pointer :: indexToCellID

        type (field2dInteger), pointer :: cellsOnCell_ptr, edgesOnCell_ptr, verticesOnCell_ptr, &
                                 cellsOnEdge_ptr, verticesOnEdge_ptr, edgesOnEdge_ptr, &
                                 cellsOnVertex_ptr, edgesOnVertex_ptr

        type (field2dInteger), pointer :: cellsOnCell, edgesOnCell, verticesOnCell, &
                                 cellsOnEdge, verticesOnEdge, edgesOnEdge, &
                                 cellsOnVertex, edgesOnVertex

        integer :: i, j


        nullify(cellsOnCell)
        nullify(edgesOnCell)
        nullify(verticesOnCell)
        nullify(cellsOnEdge)
        nullify(verticesOnEdge)
        nullify(edgesOnEdge)
        nullify(cellsOnVertex)
        nullify(edgesOnVertex)

        if (associated(cellsOnCell_save)) then
            cellsOnCell_ptr => cellsOnCell_save
            call mpas_pool_get_field(allFields, 'cellsOnCell', cellsOnCell)
        end if 
        if (associated(edgesOnCell_save)) then
            edgesOnCell_ptr => edgesOnCell_save
            call mpas_pool_get_field(allFields, 'edgesOnCell', edgesOnCell)
        end if 
        if (associated(verticesOnCell_save)) then
            verticesOnCell_ptr => verticesOnCell_save
            call mpas_pool_get_field(allFields, 'verticesOnCell', verticesOnCell)
        end if 
        if (associated(cellsOnEdge_save)) then
            cellsOnEdge_ptr => cellsOnEdge_save
            call mpas_pool_get_field(allFields, 'cellsOnEdge', cellsOnEdge)
        end if 
        if (associated(verticesOnEdge_save)) then
            verticesOnEdge_ptr => verticesOnEdge_save
            call mpas_pool_get_field(allFields, 'verticesOnEdge', verticesOnEdge)
        end if 
        if (associated(edgesOnEdge_save)) then
            edgesOnEdge_ptr => edgesOnEdge_save
            call mpas_pool_get_field(allFields, 'edgesOnEdge', edgesOnEdge)
        end if 
        if (associated(cellsOnVertex_save)) then
            cellsOnVertex_ptr => cellsOnVertex_save
            call mpas_pool_get_field(allFields, 'cellsOnVertex', cellsOnVertex)
        end if 
        if (associated(edgesOnVertex_save)) then
            edgesOnVertex_ptr => edgesOnVertex_save
            call mpas_pool_get_field(allFields, 'edgesOnVertex', edgesOnVertex)
        end if 

        !
        ! Reset indices for  connectivity arrays from global to local index space
        !
        call mpas_pool_get_field(allFields, 'indexToCellID', indexToCellID)
        do while (associated(indexToCellID))

            if (associated(cellsOnCell)) then
                deallocate(cellsOnCell % array)
                cellsOnCell % array => cellsOnCell_ptr % array
                nullify(cellsOnCell_ptr % array)
                cellsOnCell_ptr => cellsOnCell_ptr % next
                cellsOnCell => cellsOnCell % next
            end if

            if (associated(edgesOnCell)) then
                deallocate(edgesOnCell % array)
                edgesOnCell % array => edgesOnCell_ptr % array
                nullify(edgesOnCell_ptr % array)
                edgesOnCell_ptr => edgesOnCell_ptr % next
                edgesOnCell => edgesOnCell % next
            end if

            if (associated(verticesOnCell)) then
                deallocate(verticesOnCell % array)
                verticesOnCell % array => verticesOnCell_ptr % array
                nullify(verticesOnCell_ptr % array)
                verticesOnCell_ptr => verticesOnCell_ptr % next
                verticesOnCell => verticesOnCell % next
            end if

            if (associated(cellsOnEdge)) then
                deallocate(cellsOnEdge % array)
                cellsOnEdge % array => cellsOnEdge_ptr % array
                nullify(cellsOnEdge_ptr % array)
                cellsOnEdge_ptr => cellsOnEdge_ptr % next
                cellsOnEdge => cellsOnEdge % next
            end if

            if (associated(verticesOnEdge)) then
                deallocate(verticesOnEdge % array)
                verticesOnEdge % array => verticesOnEdge_ptr % array
                nullify(verticesOnEdge_ptr % array)
                verticesOnEdge_ptr => verticesOnEdge_ptr % next
                verticesOnEdge => verticesOnEdge % next
            end if

            if (associated(edgesOnEdge)) then
                deallocate(edgesOnEdge % array)
                edgesOnEdge % array => edgesOnEdge_ptr % array
                nullify(edgesOnEdge_ptr % array)
                edgesOnEdge_ptr => edgesOnEdge_ptr % next
                edgesOnEdge => edgesOnEdge % next
            end if

            if (associated(cellsOnVertex)) then
                deallocate(cellsOnVertex % array)
                cellsOnVertex % array => cellsOnVertex_ptr % array
                nullify(cellsOnVertex_ptr % array)
                cellsOnVertex_ptr => cellsOnVertex_ptr % next
                cellsOnVertex => cellsOnVertex % next
            end if

            if (associated(edgesOnVertex)) then
                deallocate(edgesOnVertex % array)
                edgesOnVertex % array => edgesOnVertex_ptr % array
                nullify(edgesOnVertex_ptr % array)
                edgesOnVertex_ptr => edgesOnVertex_ptr % next
                edgesOnVertex => edgesOnVertex % next
            end if

            indexToCellID => indexToCellID % next
        end do

        if (associated(cellsOnCell_save))    call mpas_deallocate_field(cellsOnCell_save)
        if (associated(edgesOnCell_save))    call mpas_deallocate_field(edgesOnCell_save)
        if (associated(verticesOnCell_save)) call mpas_deallocate_field(verticesOnCell_save)
        if (associated(cellsOnEdge_save))    call mpas_deallocate_field(cellsOnEdge_save)
        if (associated(verticesOnEdge_save)) call mpas_deallocate_field(verticesOnEdge_save)
        if (associated(edgesOnEdge_save))    call mpas_deallocate_field(edgesOnEdge_save)
        if (associated(cellsOnVertex_save))  call mpas_deallocate_field(cellsOnVertex_save)
        if (associated(edgesOnVertex_save))  call mpas_deallocate_field(edgesOnVertex_save)

        nullify(cellsOnCell_save)
        nullify(edgesOnCell_save)
        nullify(verticesOnCell_save)
        nullify(cellsOnEdge_save)
        nullify(verticesOnEdge_save)
        nullify(edgesOnEdge_save)
        nullify(cellsOnVertex_save)
        nullify(edgesOnVertex_save)

    end subroutine postwrite_reindex


    !-----------------------------------------------------------------------
    !  routine postread_reindex
    !
    !> \brief Reindex connectivity fields from global to local index space.
    !> \author Doug Jacobsen, Michael Duda
    !> \date   24 September 2014
    !> \details 
    !>  For any connectivity fields contained in the stream that was read, 
    !>  whose fields include those in the streamFields pool, convert the
    !>  globally indexed connectivity fields in the stream to local index space.
    !>  This routine should be called immediately after a read of a stream.
    !
    !-----------------------------------------------------------------------
    subroutine postread_reindex(allFields, streamFields)

        implicit none

        type (mpas_pool_type), pointer :: allFields
        type (mpas_pool_type), pointer :: streamFields

        type (mpas_pool_iterator_type) :: fieldItr
        type (mpas_pool_field_info_type) :: fieldInfo

        type (field1DInteger), pointer :: indexToCellID, indexToVertexID, indexToEdgeID, nEdgesOnCell, cursor
        type (field2DInteger), pointer :: int2DField
!TODO: Use a short string kind here?
        character(len=32) :: outDimName, indexSpaceName
        integer, dimension(:,:), pointer :: sortedID
        integer :: innerDim
        integer, pointer :: outerDim, indexSpaceDim
        logical :: skip_field
        integer :: i, j, k


        call mpas_pool_get_field(allFields, 'indexToCellID', indexToCellID)
        call mpas_pool_get_field(allFields, 'indexToEdgeID', indexToEdgeID)
        call mpas_pool_get_field(allFields, 'indexToVertexID', indexToVertexID)

        call mpas_pool_begin_iteration(streamFields)

        do while ( mpas_pool_get_next_member(streamFields, fieldItr) )

            ! Note: in a stream's field_pool, the names of fields are stored as configs
            if ( fieldItr % memberType == MPAS_POOL_CONFIG ) then

                call mpas_pool_get_field_info(allFields, fieldItr % memberName, fieldInfo)

                skip_field = .false.
                if (trim(fieldItr % memberName) == 'cellsOnCell') then

                    STREAM_DEBUG_WRITE('-- Reindexing cellsOnCell')

                    ! Get pointer to the field to be reindexed
                    call mpas_pool_get_field(allFields, 'cellsOnCell', int2DField)

                    ! Set the name of the outer dimension
                    outDimName = 'nCells'

                    ! Set the name of the dimension of the space within which we are indexing dimension
                    indexSpaceName = 'nCells'

                    ! Get pointer to appropriate global index field
                    cursor => indexToCellID

                else if (trim(fieldItr % memberName) == 'edgesOnCell') then

                    STREAM_DEBUG_WRITE('-- Reindexing edgesOnCell')

                    ! Get pointer to the field to be reindexed
                    call mpas_pool_get_field(allFields, 'edgesOnCell', int2DField)

                    ! Set the name of the outer dimension
                    outDimName = 'nCells'

                    ! Set the name of the dimension of the space within which we are indexing dimension
                    indexSpaceName = 'nEdges'

                    ! Get pointer to appropriate global index field
                    cursor => indexToEdgeID

                else if (trim(fieldItr % memberName) == 'verticesOnCell') then

                    STREAM_DEBUG_WRITE('-- Reindexing verticesOnCell')

                    ! Get pointer to the field to be reindexed
                    call mpas_pool_get_field(allFields, 'verticesOnCell', int2DField)

                    ! Set the name of the outer dimension
                    outDimName = 'nCells'

                    ! Set the name of the dimension of the space within which we are indexing dimension
                    indexSpaceName = 'nVertices'

                    ! Get pointer to appropriate global index field
                    cursor => indexToVertexID

                else if (trim(fieldItr % memberName) == 'cellsOnEdge') then

                    STREAM_DEBUG_WRITE('-- Reindexing cellsOnEdge')

                    ! Get pointer to the field to be reindexed
                    call mpas_pool_get_field(allFields, 'cellsOnEdge', int2DField)

                    ! Set the name of the outer dimension
                    outDimName = 'nEdges'

                    ! Set the name of the dimension of the space within which we are indexing dimension
                    indexSpaceName = 'nCells'

                    ! Get pointer to appropriate global index field
                    cursor => indexToCellID

                else if (trim(fieldItr % memberName) == 'verticesOnEdge') then

                    STREAM_DEBUG_WRITE('-- Reindexing verticesOnEdge')

                    ! Get pointer to the field to be reindexed
                    call mpas_pool_get_field(allFields, 'verticesOnEdge', int2DField)

                    ! Set the name of the outer dimension
                    outDimName = 'nEdges'

                    ! Set the name of the dimension of the space within which we are indexing dimension
                    indexSpaceName = 'nVertices'

                    ! Get pointer to appropriate global index field
                    cursor => indexToVertexID

                else if (trim(fieldItr % memberName) == 'edgesOnEdge') then

                    STREAM_DEBUG_WRITE('-- Reindexing edgesOnEdge')

                    ! Get pointer to the field to be reindexed
                    call mpas_pool_get_field(allFields, 'edgesOnEdge', int2DField)

                    ! Set the name of the outer dimension
                    outDimName = 'nEdges'

                    ! Set the name of the dimension of the space within which we are indexing dimension
                    indexSpaceName = 'nEdges'

                    ! Get pointer to appropriate global index field
                    cursor => indexToEdgeID

                else if (trim(fieldItr % memberName) == 'cellsOnVertex') then

                    STREAM_DEBUG_WRITE('-- Reindexing cellsOnVertex')

                    ! Get pointer to the field to be reindexed
                    call mpas_pool_get_field(allFields, 'cellsOnVertex', int2DField)

                    ! Set the name of the outer dimension
                    outDimName = 'nVertices'

                    ! Set the name of the dimension of the space within which we are indexing dimension
                    indexSpaceName = 'nCells'

                    ! Get pointer to appropriate global index field
                    cursor => indexToCellID

                else if (trim(fieldItr % memberName) == 'edgesOnVertex') then

                    STREAM_DEBUG_WRITE('-- Reindexing edgesOnVertex')

                    ! Get pointer to the field to be reindexed
                    call mpas_pool_get_field(allFields, 'edgesOnVertex', int2DField)

                    ! Set the name of the outer dimension
                    outDimName = 'nVertices'

                    ! Set the name of the dimension of the space within which we are indexing dimension
                    indexSpaceName = 'nEdges'

                    ! Get pointer to appropriate global index field
                    cursor => indexToEdgeID

                else
                    skip_field = .true.
                end if

                if (.not. skip_field) then

                    ! Get inner dimension of field to be reindexed (assumed to be block invariant)
                    innerDim = int2DField % dimSizes(1)

                    ! Reindex all blocks for the field
                    do while (associated(int2DField))

                        ! Get outer dimension of field for this block
                        call mpas_pool_get_dimension(cursor % block % dimensions, trim(outDimName), outerDim)
                        call mpas_pool_get_dimension(cursor % block % dimensions, trim(indexSpaceName), indexSpaceDim)

                        ! Set-up reindexing map
                        allocate(sortedID(2,indexSpaceDim))
                        do i = 1, indexSpaceDim
                            sortedID(1,i) = cursor % array(i)
                            sortedID(2,i) = i
                        end do
                        call mpas_quicksort(indexSpaceDim, sortedID)

                        ! Reindex the field
                        do i = 1, outerDim
                            do j = 1, innerDim
                                k = mpas_binary_search(sortedID, 2, 1, indexSpaceDim, int2DField % array(j,i))
                                if (k <= indexSpaceDim) then
                                    int2DField % array(j,i) = sortedID(2,k)
                                else
                                    int2DField % array(j,i) = indexSpaceDim + 1
                                end if
                            end do
                        end do

                        deallocate(sortedID)
                        int2DField => int2DField % next
                        cursor => cursor % next

                    end do

                end if

            end if

        end do

    end subroutine postread_reindex
 

end module mpas_stream_manager


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! C interface routines for building streams at run-time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine stream_mgr_create_stream_c(manager_c, streamID_c, direction_c, filename_c, ref_time_c, rec_intv_c, &
                                      maxRecords_c, immutable_c, ierr_c) bind(c) !{{{

    use mpas_c_interfacing, only : mpas_c_to_f_string
    use iso_c_binding, only : c_char, c_int, c_ptr, c_f_pointer
    use mpas_stream_manager, only : MPAS_streamManager_type, MPAS_STREAM_MGR_NOERR, MPAS_STREAM_PROPERTY_FILENAME, &
                                    MPAS_STREAM_PROPERTY_MAX_RECORDS, MPAS_STREAM_PROPERTY_REF_TIME, &
                                    MPAS_STREAM_PROPERTY_RECORD_INTV, MPAS_stream_mgr_create_stream, MPAS_stream_mgr_set_property
    use mpas_kind_types, only : StrKIND

    implicit none

    type (c_ptr) :: manager_c
    character(kind=c_char) :: streamID_c(*)
    integer(kind=c_int) :: direction_c
    character(kind=c_char) :: filename_c(*)
    character(kind=c_char) :: ref_time_c(*)
    character(kind=c_char) :: rec_intv_c(*)
    integer(kind=c_int) :: maxRecords_c
    integer(kind=c_int) :: immutable_c
    integer(kind=c_int) :: ierr_c

    type (MPAS_streamManager_type), pointer :: manager
    character(len=StrKIND) :: streamID, filename, reference_time, record_interval
    integer :: direction, maxRecords, immutable, ierr


    call c_f_pointer(manager_c, manager)
    call mpas_c_to_f_string(streamID_c, streamID)
    direction = direction_c
    call mpas_c_to_f_string(filename_c, filename)
    call mpas_c_to_f_string(ref_time_c, reference_time)
    call mpas_c_to_f_string(rec_intv_c, record_interval)
    immutable = immutable_c
    maxRecords = maxRecords_c
    
    !
    ! For immutable streams, the stream should have already been defined at this point, and
    !    all we need to do is update the stream's filename template and maximum records per file;
    !    otherwise, we need to create a new stream
    !
    ierr = 0
    if (immutable == 1) then
        call MPAS_stream_mgr_set_property(manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename, ierr=ierr)
        call MPAS_stream_mgr_set_property(manager, streamID, MPAS_STREAM_PROPERTY_MAX_RECORDS, maxRecords, ierr=ierr)
    else
        call MPAS_stream_mgr_create_stream(manager, streamID, direction, filename, maxRecords, ierr=ierr)
    end if

    if (reference_time /= 'none') then
        call MPAS_stream_mgr_set_property(manager, streamID, MPAS_STREAM_PROPERTY_REF_TIME, reference_time, ierr=ierr)
write(0,*) 'For stream '//trim(streamID)//' found reference time '//trim(reference_time)
    end if

    if (record_interval /= 'none') then
        call MPAS_stream_mgr_set_property(manager, streamID, MPAS_STREAM_PROPERTY_RECORD_INTV, record_interval, ierr=ierr)
write(0,*) 'For stream '//trim(streamID)//' found record interval '//trim(record_interval)
    end if

    if (ierr == MPAS_STREAM_MGR_NOERR) then
        ierr_c = 0
    else
        ierr_c = 1
    end if

end subroutine stream_mgr_create_stream_c !}}}


subroutine stream_mgr_add_field_c(manager_c, streamID_c, fieldName_c, ierr_c) bind(c)

    use mpas_c_interfacing, only : mpas_c_to_f_string
    use iso_c_binding, only : c_char, c_int, c_ptr, c_f_pointer
    use mpas_stream_manager, only : MPAS_streamManager_type, MPAS_STREAM_MGR_NOERR, MPAS_stream_mgr_add_field
    use mpas_kind_types, only : StrKIND

    implicit none

    type (c_ptr) :: manager_c
    character(kind=c_char) :: streamID_c(*)
    character(kind=c_char) :: fieldName_c(*)
    integer(kind=c_int) :: ierr_c

    type (MPAS_streamManager_type), pointer :: manager
    character(len=StrKIND) :: streamID, fieldName
    integer :: ierr


    call c_f_pointer(manager_c, manager)
    call mpas_c_to_f_string(streamID_c, streamID)
    call mpas_c_to_f_string(fieldName_c, fieldName)

    call MPAS_stream_mgr_add_field(manager, streamID, fieldName, ierr)

    if (ierr == MPAS_STREAM_MGR_NOERR) then
        ierr_c = 0
    else
        ierr_c = 1
    end if

end subroutine stream_mgr_add_field_c


subroutine stream_mgr_add_alarm_c(manager_c, streamID_c, direction_c, alarmTime_c, alarmInterval_c, ierr_c) bind(c)

    use mpas_c_interfacing, only : mpas_c_to_f_string
    use iso_c_binding, only : c_char, c_int, c_ptr, c_f_pointer
    use mpas_stream_manager, only : MPAS_streamManager_type, MPAS_STREAM_MGR_NOERR, MPAS_STREAM_INPUT, MPAS_STREAM_OUTPUT, &
                                    MPAS_stream_mgr_get_clock, MPAS_stream_mgr_add_alarm
    use mpas_kind_types, only : StrKIND
    use mpas_timekeeping, only : MPAS_Clock_type, MPAS_START_TIME, MPAS_Time_type, MPAS_TimeInterval_type, mpas_add_clock_alarm, &
                                 mpas_set_time, mpas_set_timeInterval, mpas_get_clock_time

    implicit none

    type (c_ptr) :: manager_c
    character(kind=c_char) :: streamID_c(*)
    character(kind=c_char) :: direction_c(*)
    character(kind=c_char) :: alarmTime_c(*)
    character(kind=c_char) :: alarmInterval_c(*)
    integer(kind=c_int) :: ierr_c

    type (MPAS_streamManager_type), pointer :: manager
    type (MPAS_Clock_type), pointer :: clock
    character(len=StrKIND) :: streamID, direction, alarmID, alarmTime, alarmInterval
    type (MPAS_Time_type) :: alarmTime_local
    type (MPAS_TimeInterval_type) :: alarmInterval_local
    integer :: idirection
    integer :: ierr


    ierr = 0

    call c_f_pointer(manager_c, manager)
    call mpas_c_to_f_string(streamID_c, streamID)
    call mpas_c_to_f_string(direction_c, direction)
    call mpas_c_to_f_string(alarmTime_c, alarmTime)
    call mpas_c_to_f_string(alarmInterval_c, alarmInterval)
    write(alarmID, '(a)') trim(streamID)//'_'//trim(direction)

    ! Nothing to do for this stream
    if (trim(alarmInterval) == 'none') then
write(0,*) 'No need to add an alarm to stream '//trim(streamID)//' direction '//trim(direction)
        return
    end if

    if (trim(direction) == 'input') then
        idirection = MPAS_STREAM_INPUT
    else if (trim(direction) == 'output') then
        idirection = MPAS_STREAM_OUTPUT
    end if

write(0,*) 'Adding alarm to stream '//trim(streamID)//' direction '//trim(direction)
write(0,*) '   ID:       '//trim(alarmID)
write(0,*) '   Time:     '//trim(alarmTime)
write(0,*) '   Interval: '//trim(alarmInterval)

    call MPAS_stream_mgr_get_clock(manager, clock)

    if (trim(alarmTime) == 'start') then
        alarmTime_local = mpas_get_clock_time(clock, MPAS_START_TIME, ierr=ierr)
    else
        call mpas_set_time(alarmTime_local, dateTimeString=alarmTime)
    end if

    if (trim(alarmInterval) == 'initial_only') then
        call mpas_add_clock_alarm(clock, alarmID, alarmTime_local, ierr=ierr)
    else
        call mpas_set_timeInterval(alarmInterval_local, timeString=alarmInterval)
        call mpas_add_clock_alarm(clock, alarmID, alarmTime_local, alarmTimeInterval=alarmInterval_local, ierr=ierr)
    end if

    call MPAS_stream_mgr_add_alarm(manager, streamID, alarmID, idirection, ierr=ierr)

    if (ierr == MPAS_STREAM_MGR_NOERR) then
        ierr_c = 0
    else
        ierr_c = 1
    end if

end subroutine stream_mgr_add_alarm_c


subroutine stream_mgr_add_pkg_c(manager_c, streamID_c, package_c, ierr_c) bind(c)

    use mpas_c_interfacing, only : mpas_c_to_f_string
    use iso_c_binding, only : c_char, c_int, c_ptr, c_f_pointer
    use mpas_stream_manager, only : MPAS_streamManager_type, MPAS_STREAM_MGR_NOERR, MPAS_stream_mgr_add_pkg
    use mpas_kind_types, only : StrKIND

    implicit none

    type (c_ptr) :: manager_c
    character(kind=c_char) :: streamID_c(*)
    character(kind=c_char) :: package_c(*)
    integer(kind=c_int) :: ierr_c

    type (MPAS_streamManager_type), pointer :: manager
    character(len=StrKIND) :: streamID, package
    integer :: idirection
    integer :: ierr


    ierr = 0

    call c_f_pointer(manager_c, manager)
    call mpas_c_to_f_string(streamID_c, streamID)
    call mpas_c_to_f_string(package_c, package)
    write(package, '(a)') trim(package)//'Active'

write(0,*) 'Adding package '//trim(package)//' to stream '//trim(streamID)

    call MPAS_stream_mgr_add_pkg(manager, streamID, package, ierr=ierr)

    if (ierr == MPAS_STREAM_MGR_NOERR) then
        ierr_c = 0
    else
        ierr_c = 1
    end if

end subroutine stream_mgr_add_pkg_c
