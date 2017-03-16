module mpas_stream_manager

#define COMMA ,
#define STREAM_DEBUG_WRITE(M) ! write(stderrUnit,*) M
#define STREAM_WARNING_WRITE(M) write(stderrUnit,*) 'WARNING: '//M
#define STREAM_ERROR_WRITE(M) write(stderrUnit,*) 'ERROR: '//M

    use mpas_kind_types
    use mpas_derived_types
    use mpas_field_routines
    use mpas_pool_routines
    use mpas_timekeeping
    use mpas_io_units
    use mpas_io_streams
    use mpas_stream_list
    use mpas_sort
    use mpas_threading
    use mpas_abort, only : mpas_dmpar_global_abort


    public :: MPAS_stream_mgr_init, &
              MPAS_stream_mgr_finalize, &
              MPAS_stream_mgr_create_stream, &
              MPAS_stream_mgr_destroy_stream, &
              MPAS_stream_mgr_get_clock, &
              MPAS_stream_mgr_set_property, &
              MPAS_stream_mgr_get_property, &
              MPAS_stream_mgr_add_pkg, &
              MPAS_stream_mgr_remove_pkg, &
              MPAS_stream_mgr_add_pool, &
              MPAS_stream_mgr_add_field, &
              MPAS_stream_mgr_add_stream_fields, &
              MPAS_stream_mgr_remove_field, &
              MPAS_stream_mgr_add_alarm, &
              MPAS_stream_mgr_remove_alarm, &
              MPAS_stream_mgr_reset_alarms, &
              MPAS_stream_mgr_ringing_alarms, &
              MPAS_stream_mgr_add_att, &
              MPAS_stream_mgr_write, &
              MPAS_stream_mgr_block_write, &
              MPAS_stream_mgr_read, &
              MPAS_stream_mgr_begin_iteration, &
              MPAS_stream_mgr_get_next_stream, &
              MPAS_stream_mgr_get_next_field, &
              MPAS_stream_mgr_stream_exists, &
              MPAS_get_stream_filename, &
              MPAS_build_stream_filename

    private

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
    subroutine MPAS_stream_mgr_init(manager, ioContext, clock, allFields, allPackages, allStructs, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_init'

        type (MPAS_streamManager_type), pointer :: manager
        type (MPAS_IO_context_type), pointer :: ioContext
        type (MPAS_Clock_type), pointer :: clock
        type (MPAS_Pool_type), pointer :: allFields
        type (MPAS_Pool_type), pointer :: allPackages
        type (MPAS_Pool_type), pointer :: allStructs
        integer, intent(out), optional :: ierr

        integer :: err_local, threadNum

        threadNum = mpas_threading_get_thread_num()

        if ( threadNum == 0 ) then
           call seed_random()

           STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_init()')

           if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

           allocate(manager)
           manager % ioContext => ioContext
           manager % allFields => allFields
           manager % allPackages => allPackages
           manager % allStructs => allStructs
           manager % streamClock => clock
           manager % numStreams = 0
           manager % errorLevel = MPAS_STREAM_ERR_SILENT

           !
           ! Set up linked list of streams
           !
           call MPAS_stream_list_create(manager % streams, ierr=err_local)
           if (err_local /= MPAS_STREAM_LIST_NOERR) then
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               STREAM_ERROR_WRITE('Problems while creating stream list')
               return 
           end if

           !
           ! Set up linked list of input alarms
           !
           call MPAS_stream_list_create(manager % alarms_in, ierr=err_local)
           if (err_local /= MPAS_STREAM_LIST_NOERR) then
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               STREAM_ERROR_WRITE('Problems while creating input alarm list')
               return 
           end if

           !
           ! Set up linked list of output alarms
           !
           call MPAS_stream_list_create(manager % alarms_out, ierr=err_local)
           if (err_local /= MPAS_STREAM_LIST_NOERR) then
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               STREAM_ERROR_WRITE('Problems while creating output alarm list')
               return 
           end if

           !
           ! Create a pool to hold default global attributes that every stream will have
           !
           call mpas_pool_create_pool(manager % defaultAtts)
        end if

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

        integer :: err_local, threadNum
        type (MPAS_stream_list_type), pointer :: stream_cursor

        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_finalize()') 

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        if ( threadNum == 0 ) then
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
               STREAM_ERROR_WRITE('Problems while destroying stream list')
           end if

           !
           ! Free up list of input alarms
           !
           call MPAS_stream_list_destroy(manager % alarms_in, ierr=err_local)
           if (err_local /= MPAS_STREAM_LIST_NOERR) then
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               STREAM_ERROR_WRITE('Problems while destroying input alarms list')
           end if

           !
           ! Free up list of output alarms
           !
           call MPAS_stream_list_destroy(manager % alarms_out, ierr=err_local)
           if (err_local /= MPAS_STREAM_LIST_NOERR) then
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               STREAM_ERROR_WRITE('Problems while destroying output alarms list')
           end if

           !
           ! Free up default attribute pool
           !
           call mpas_pool_destroy_pool(manager % defaultAtts)

           deallocate(manager)
        end if

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
    !>  Knowing the interval between files, and 
    !>  the filename template, a "referenceTime" argument must be provided to 
    !>  specify the first timestamp appearing in any of the files associated with 
    !>  the stream, thereby determining where the "file breaks" will occur between 
    !>  timestamps. If no "referenceTime" is specified, the start time of the 
    !>  clock associated with the stream handler will be used as the reference 
    !>  time. Additionally, the interval between records in the file may be
    !>  specified using the optional "recordInterval" argument; if this argument
    !>  is not supplied, the stream manager will assume that this interval is
    !>  equal to the shortest period of any periodic alarm attached to the stream.
    !>  The optional argument 'realPrecision' specifies the precision of
    !>  real-valued fields in the files associated with the stream; this
    !>  argument may take on values MPAS_IO_SINGLE_PRECISION,
    !>  MPAS_IO_DOUBLE_PRECISION, or MPAS_IO_NATIVE_PRECISION; if this argument is
    !>  not supplied, native precision is assumed.
    !>  The optional argument clobberMode determines how the stream manager will
    !>  deal with existing files; possible options include MPAS_STREAM_CLOBBER_NEVER, 
    !>  MPAS_STREAM_CLOBBER_APPEND, MPAS_STREAM_CLOBBER_TRUNCATE, 
    !>  and MPAS_STREAM_CLOBBER_OVERWRITE. The default behavior is to never modify
    !>  existing files (MPAS_STREAM_CLOBBER_NEVER).
    !>  The optional argument ioType specifies the I/O type to use for the
    !>  stream; possible values are defined by constants in the mpas_io module and
    !>  include: MPAS_IO_NETCDF, MPAS_IO_NETCDF4, MPAS_IO_PNETCDF, and
    !>  MPAS_IO_PNETCDF5. If not specified, the io_type will default to
    !>  MPAS_IO_PNETCDF.
    !>  NOTE: This routine does not support regular expressions for StreamID
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_create_stream(manager, streamID, direction, filename, &
                                             filenameInterval, referenceTime, recordInterval, &
                                             realPrecision, clobberMode, ioType, ierr) !{{{

        use mpas_io, only : MPAS_IO_PNETCDF

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_create_stream'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: direction
        character (len=*), intent(in) :: filename
        character (len=*), intent(in), optional :: filenameInterval
        type (MPAS_Time_type), intent(in), optional :: referenceTime
        type (MPAS_TimeInterval_type), intent(in), optional :: recordInterval
        integer, intent(in), optional :: realPrecision
        integer, intent(in), optional :: clobberMode
        integer, intent(in), optional :: ioType
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: new_stream
        integer :: err_local, threadNum


        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_create_stream() for '//trim(streamID))

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        if ( threadNum == 0 ) then
           !
           ! Check that the stream does not already exist
           !
           nullify(new_stream)
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

           ! Filename interval is 'none' by deault, but is set through the set_property routine.
           if (present(filenameInterval)) then
               new_stream % filename_interval = filenameInterval
           else
               new_stream % filename_interval = 'none'
           end if

           new_stream % nRecords = 0
           if (present(clobberMode)) then
               new_stream % clobber_mode = clobberMode
           else
               new_stream % clobber_mode = MPAS_STREAM_CLOBBER_NEVER
           end if
           if (present(ioType)) then
               new_stream % io_type = ioType
           else
               new_stream % io_type = MPAS_IO_PNETCDF
           end if
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
           if (present(realPrecision)) then
               new_stream % precision = realPrecision
           end if
           call MPAS_stream_list_create(new_stream % alarmList_in, ierr=err_local)
           if (err_local /= MPAS_STREAM_LIST_NOERR) then
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               STREAM_ERROR_WRITE('Problems while creating input alarm list')
               deallocate(new_stream)
               return 
           end if
           call MPAS_stream_list_create(new_stream % alarmList_out, ierr=err_local)
           if (err_local /= MPAS_STREAM_LIST_NOERR) then
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               STREAM_ERROR_WRITE('Problems while creating output alarm list')
               deallocate(new_stream)
               return 
           end if
           call mpas_pool_create_pool(new_stream % att_pool)
           call mpas_pool_clone_pool(manager % defaultAtts, new_stream % att_pool)
           call mpas_pool_create_pool(new_stream % field_pool)
           call mpas_pool_create_pool(new_stream % field_pkg_pool)
           call mpas_pool_create_pool(new_stream % pkg_pool)
           nullify(new_stream % next)


           !
           ! Add stream to list
           !
           call MPAS_stream_list_insert(manager % streams, new_stream, ierr=err_local)
           if (err_local /= MPAS_STREAM_LIST_NOERR) then
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               STREAM_ERROR_WRITE('Problems while adding stream to list')
               return
           end if
           
           manager % numStreams = manager % numStreams + 1
        end if

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
    !>  NOTE: This routine does not support regular expressions for StreamID
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_destroy_stream(manager, streamID, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_destroy_stream'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(out), optional :: ierr

        integer :: err_local, threadNum
        type (MPAS_stream_list_type), pointer :: stream, alarm_cursor, delete_me

        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_destroy_stream()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        if ( threadNum == 0 ) then
           !
           ! Remove stream from list
           !
           call MPAS_stream_list_remove(manager % streams, streamID, stream, ierr=err_local)
           if (err_local /= MPAS_STREAM_LIST_NOERR) then
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               STREAM_ERROR_WRITE('Problems while removing stream from list')
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
                   STREAM_ERROR_WRITE('Problems while removing stream from list of input alarm')
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
                   STREAM_ERROR_WRITE('Problems while removing stream from list of output alarm')
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
           call mpas_pool_destroy_pool(stream % field_pkg_pool)
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
                   STREAM_ERROR_WRITE('Problems while closing stream '//trim(stream % name))
               end if
               deallocate(stream % stream)
           end if
           deallocate(stream)
           
           manager % numStreams = manager % numStreams - 1
        end if

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
    subroutine MPAS_stream_mgr_get_clock(manager, clock, ierr) !{{{

        implicit none

        type (MPAS_streamManager_type), intent(in) :: manager
        type (MPAS_Clock_type), pointer :: clock
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_get_clock()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        clock => manager % streamClock

    end subroutine MPAS_stream_mgr_get_clock !}}}

    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_pool
    !
    !> \brief Add a pool of fields to the specified stream in an MPAS stream manager.
    !> \author Doug Jacobsen, Michael Duda
    !> \date   09/15/2014
    !> \details
    !>  Adds a pool from the allStructs pool to the specified stream in an MPAS
    !>   stream manager. Currently, it adds only explicitly named var's and
    !>   var_array's to the stream, but commented code will allow adding all nested
    !>   structs as well. If the optional 'packages' argument is supplied, the
    !>   specified packages will be attached to all var and var_array members
    !>   added to the stream.
    !>   NOTE: This routine does not support regular expressions for StreamID
    !
    !-----------------------------------------------------------------------
    recursive subroutine MPAS_stream_mgr_add_pool(manager, streamID, poolName, packages, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_pool'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: poolName
        character (len=*), intent(in), optional :: packages
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream
        type (mpas_pool_field_info_type) :: info
        integer, pointer :: test_ptr
        integer :: err_local, threadNum

        type (mpas_pool_type), pointer :: fieldPool
        type (mpas_pool_iterator_type) :: poolItr

        type (field0DReal), pointer :: real0DField
        type (field1DReal), pointer :: real1DField
        type (field2DReal), pointer :: real2DField
        type (field3DReal), pointer :: real3DField
        type (field4DReal), pointer :: real4DField
        type (field5DReal), pointer :: real5DField
        type (field0DInteger), pointer :: int0DField
        type (field1DInteger), pointer :: int1DField
        type (field2DInteger), pointer :: int2DField
        type (field3DInteger), pointer :: int3DField
        type (field0DChar), pointer :: char0DField
        type (field1DChar), pointer :: char1DField

        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_pool()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        if ( threadNum == 0 ) then
           !
           ! Check that stream exists
           !
           nullify(stream)
           if (.not. MPAS_stream_list_query(manager % streams, streamID, stream, ierr=err_local)) then
               STREAM_ERROR_WRITE('Requested stream '//trim(streamID)//' does not exist in stream manager')
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               return
           end if

           !
           ! Don't modify an immutable stream
           !
           if (stream % immutable) then
               STREAM_ERROR_WRITE('Requested stream '//trim(streamID)//' is immutable.')
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               return
           end if

           !
           ! Check that the pool exists
           !
           call mpas_pool_get_subpool(manager % allStructs, poolName, fieldPool)
           if (.not. associated(fieldPool) ) then
               STREAM_ERROR_WRITE('Requested pool '//trim(poolName)//' does not exist.')
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               return
           end if

           !
           ! Iterate over pool, adding each field to the stream, and recursively calling this subroutine for each subpool
           !
           call mpas_pool_begin_iteration(fieldPool)
           do while (mpas_pool_get_next_member(fieldPool, poolItr))
               if (poolItr % memberType == MPAS_POOL_SUBPOOL) then
                   STREAM_DEBUG_WRITE('-- Try to add subpool...')
                   ! call mpas_stream_mgr_add_pool(manager, streamId, poolItr % memberName, packages=packages, ierr=ierr)
               else if (poolItr % memberType == MPAS_POOL_FIELD) then
                   if (poolItr % dataType == MPAS_POOL_REAL) then
                       if (poolItr % nDims == 0) then
                           call mpas_pool_get_field(fieldPool, poolItr % memberName, real0DField, 1)
                           call mpas_stream_mgr_add_field(manager, streamID, real0DField % fieldName, packages=packages, ierr=ierr)
                       else if (poolItr % nDims == 1) then
                           call mpas_pool_get_field(fieldPool, poolItr % memberName, real1DField, 1)
                           call mpas_stream_mgr_add_field(manager, streamID, real1DField % fieldName, packages=packages, ierr=ierr)
                       else if (poolItr % nDims == 2) then
                           call mpas_pool_get_field(fieldPool, poolItr % memberName, real2DField, 1)
                           call mpas_stream_mgr_add_field(manager, streamID, real2DField % fieldName, packages=packages, ierr=ierr)
                       else if (poolItr % nDims == 3) then
                           call mpas_pool_get_field(fieldPool, poolItr % memberName, real3DField, 1)
                           call mpas_stream_mgr_add_field(manager, streamID, real3DField % fieldName, packages=packages, ierr=ierr)
                       else if (poolItr % nDims == 4) then
                           call mpas_pool_get_field(fieldPool, poolItr % memberName, real4DField, 1)
                           call mpas_stream_mgr_add_field(manager, streamID, real4DField % fieldName, packages=packages, ierr=ierr)
                       else if (poolItr % nDims == 5) then
                           call mpas_pool_get_field(fieldPool, poolItr % memberName, real5DField, 1)
                           call mpas_stream_mgr_add_field(manager, streamID, real5DField % fieldName, packages=packages, ierr=ierr)
                       end if
                   else if (poolItr % dataType == MPAS_POOL_INTEGER) then
                       if (poolItr % nDims == 0) then
                           call mpas_pool_get_field(fieldPool, poolItr % memberName, int0DField, 1)
                           call mpas_stream_mgr_add_field(manager, streamID, int0DField % fieldName, packages=packages, ierr=ierr)
                       else if (poolItr % nDims == 1) then
                           call mpas_pool_get_field(fieldPool, poolItr % memberName, int1DField, 1)
                           call mpas_stream_mgr_add_field(manager, streamID, int1DField % fieldName, packages=packages, ierr=ierr)
                       else if (poolItr % nDims == 2) then
                           call mpas_pool_get_field(fieldPool, poolItr % memberName, int2DField, 1)
                           call mpas_stream_mgr_add_field(manager, streamID, int2DField % fieldName, packages=packages, ierr=ierr)
                       else if (poolItr % nDims == 3) then
                           call mpas_pool_get_field(fieldPool, poolItr % memberName, int3DField, 1)
                           call mpas_stream_mgr_add_field(manager, streamID, int3DField % fieldName, packages=packages, ierr=ierr)
                       end if
                   else if (poolItr % dataType == MPAS_POOL_CHARACTER) then
                       if (poolItr % nDims == 0) then
                           call mpas_pool_get_field(fieldPool, poolItr % memberName, char0DField, 1)
                           call mpas_stream_mgr_add_field(manager, streamID, char0DField % fieldName, packages=packages, ierr=ierr)
                       else if (poolItr % nDims == 1) then
                           call mpas_pool_get_field(fieldPool, poolItr % memberName, char1DField, 1)
                           call mpas_stream_mgr_add_field(manager, streamID, char1DField % fieldName, packages=packages, ierr=ierr)
                       end if
                   end if
               end if
           end do
        end if

    end subroutine MPAS_stream_mgr_add_pool!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_field
    !
    !> \brief Add a field to the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Adds a field from the allFields pool to a stream. If the optional
    !>  argument 'packages' is present, those packages will be attached to
    !>  the field as well.
    !>  NOTE: This routine does not support regular expressions for StreamID
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_field(manager, streamID, fieldName, packages, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_field'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: fieldName
        character (len=*), intent(in), optional :: packages
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream
        type (mpas_pool_field_info_type) :: info
        character (len=StrKIND) :: field_pkg
        integer, pointer :: test_ptr
        logical :: test_logical
        integer :: err_level
        integer :: err_local, threadNum

        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_field()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        if ( threadNum == 0 ) then
           !
           ! Check that stream exists
           !
           nullify(stream)
           if (.not. MPAS_stream_list_query(manager % streams, streamID, stream, ierr=err_local)) then
               STREAM_ERROR_WRITE('Requested stream '//trim(streamID)//' does not exist in stream manager')
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               return
           end if

           !
           ! Don't modify an immutable stream
           !
           if (stream % immutable) then
               STREAM_ERROR_WRITE('Requested stream '//trim(streamID)//' is immutable.')
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               return
           end if

           !
           ! Check that field exists
           !
           info % nDims = -1
           call mpas_pool_get_field_info(manager % allFields, fieldName, info)
           if (info % nDims == -1) then
               STREAM_ERROR_WRITE('Requested field '//trim(fieldName)//' not available')
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               return
           end if

           !
           ! Check that the field does not already exist in the stream
           !
           nullify(test_ptr)
           err_level = mpas_pool_get_error_level()
           call mpas_pool_set_error_level(MPAS_POOL_SILENT)
           call mpas_pool_get_config(stream % field_pool, fieldName, value=test_ptr)
           call mpas_pool_set_error_level(err_level)
           if (associated(test_ptr)) then
               STREAM_ERROR_WRITE('Requested field '//trim(fieldName)//' already in stream '//trim(streamID))
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               return
           end if

           !
           ! Validate packages
           !
           if (present(packages)) then
               test_logical = parse_package_list(manager % allPackages, trim(packages), err_local)
               if (err_local /= 0) then
                   STREAM_WARNING_WRITE('One or more packages in '''//trim(packages)//''' attached to field '''//trim(fieldName)//''' is undefined')
               end if
           end if

           !
           ! Add field to field pool in stream if the field is activated
           !
           if (info % isActive) then
               call mpas_pool_add_config(stream % field_pool, fieldName, 1)

               if (present(packages)) then
                  STREAM_DEBUG_WRITE('-- Attaching packages '//trim(packages)//' to field '//trim(fieldName)//' in stream '//trim(streamID))
                  write(field_pkg,'(a)') trim(fieldName)//':packages'
                  call mpas_pool_add_config(stream % field_pkg_pool, field_pkg, packages)
               end if
           else
               write(stderrUnit, *) ' * Requested field '//trim(fieldName)//' is deactivated due to packages, or is a scratch variable.'
           end if
        end if

    end subroutine MPAS_stream_mgr_add_field!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_stream_fields
    !
    !> \brief Add all fields from another stream to the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   5 November 2014
    !> \details
    !>  Adds all fields from another specified stream into the new specified stream.
    !>  Both streams need to exist within the same stream manager. If the
    !>  optional 'packages' argument is supplied, the specified packages will be
    !>  attached to all fields added from refStreamID.
    !>  NOTE: This routine does not support regular expressions for StreamID
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_stream_fields(manager, streamID, refStreamID, packages, ierr)!{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_add_stream_fields'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID  !< stream to have fields added to
        character (len=*), intent(in) :: refStreamID  !< stream to supply list of fields to add
        character (len=*), intent(in), optional :: packages
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream, refStream
        type (mpas_pool_field_info_type) :: info
        type (mpas_pool_iterator_type) :: itr
        character (len=StrKIND) :: field_pkg
        integer, pointer :: test_ptr
        logical :: test_logical
        integer :: err_level
        integer :: err_local, threadNum

        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_stream_fields()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Check that reference stream exists
        !
        nullify(refStream)
        if (.not. MPAS_stream_list_query(manager % streams, refStreamID, refStream, ierr=err_local)) then
            STREAM_ERROR_WRITE('Requested reference stream '//trim(refStreamID)//' does not exist in stream manager')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Check that stream exists
        !
        nullify(stream)
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream, ierr=err_local)) then
            STREAM_ERROR_WRITE('Requested stream '//trim(streamID)//' does not exist in stream manager')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Don't modify an immutable stream
        !
        if (stream % immutable) then
            STREAM_ERROR_WRITE('Requested stream '//trim(streamID)//' is immutable.')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Validate packages
        !
        if (present(packages)) then
            test_logical = parse_package_list(manager % allPackages, trim(packages), err_local)
            if (err_local /= 0) then
                STREAM_WARNING_WRITE('One or more packages in '''//trim(packages)//''' attached to fields from stream '''//trim(refStreamID)//''' is undefined')
            end if
        end if

        !
        ! Loop over all fields in refStream and add them one by one to stream
        !
        err_level = mpas_pool_get_error_level()
        call mpas_pool_set_error_level(MPAS_POOL_SILENT)

        if ( threadNum == 0 ) then
           call mpas_pool_begin_iteration(refStream % field_pool)
           do while (mpas_pool_get_next_member(refStream % field_pool, itr))
               if ( itr % memberType == MPAS_POOL_CONFIG ) then
                   if ( itr % dataType == MPAS_POOL_INTEGER ) then

                       !
                       ! Check that field exists
                       !
                       info % nDims = -1
                       call mpas_pool_get_field_info(manager % allFields, itr % memberName, info)
                       if (info % nDims == -1) then
                           STREAM_ERROR_WRITE('Requested field '//trim(itr % memberName)//' not available')
                           if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                           return
                       end if

                       ! Test that the field does not already exist in stream
                       nullify(test_ptr)
                       call mpas_pool_get_config(stream % field_pool, itr % memberName, value=test_ptr)

                       if ( associated(test_ptr) ) then
                           STREAM_ERROR_WRITE('Requested field '//trim(itr % memberName)//' already in stream '//trim(streamID))
                           if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                       end if

                       if ( info % isActive ) then
                           call mpas_pool_add_config(stream % field_pool, itr % memberName, 1)

                           if (present(packages)) then
                              STREAM_DEBUG_WRITE('-- Attaching packages '//trim(packages)//' to field '//trim(itr % memberName)//' in stream '//trim(streamID))
                              write(field_pkg,'(a)') trim(itr % memberName)//':packages'
                              call mpas_pool_add_config(stream % field_pkg_pool, field_pkg, packages)
                           end if
                       else
                           write(stderrUnit, *) ' * Requested field '//trim(itr % memberName)//' is deactivated due to packages, or is a scratch variable.'
                       end if

                   end if
               end if
           end do
        end if
        call mpas_pool_set_error_level(err_level)

    end subroutine MPAS_stream_mgr_add_stream_fields!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_remove_field
    !
    !> \brief Remove a field from the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Removes a field from a stream.
    !>  NOTE: This routine does not support regular expressions for StreamID
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
        integer :: err_local, threadNum

        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_remove_field()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Check that stream exists
        !
        nullify(stream)
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream, ierr=err_local)) then
            STREAM_ERROR_WRITE('Requested stream '//trim(streamID)//' does not exist in stream manager')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Don't modify an immutable stream
        !
        if (stream % immutable) then
            STREAM_ERROR_WRITE('Requested stream '//trim(streamID)//' is immutable.')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Check that field exists in stream's field pool
        !
        nullify(test_ptr)
        call mpas_pool_get_config(stream % field_pool, fieldName, value=test_ptr)
        if (.not. associated(test_ptr)) then
            STREAM_ERROR_WRITE('Requested field '//trim(fieldName)//' not in stream '//trim(streamID))
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Remove field from stream's field pool
        !
        if ( threadNum == 0 ) then
           call mpas_pool_remove_config(stream % field_pool, fieldName)
        end if

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
    !>  NOTE: This routine does not support regular expressions for StreamID
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
        integer :: err_local, threadNum

        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_alarm()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Check that stream exists
        !
        nullify(stream)
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream, ierr=err_local)) then
            STREAM_ERROR_WRITE('Requested stream '//trim(streamID)//' does not exist in stream manager')
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

            STREAM_ERROR_WRITE('Attempting to add an alarm '//trim(alarmID)//' to invalid direction for stream '//trim(streamID))
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Check that alarm exists on clock
        !
        if (.not. mpas_is_alarm_defined(manager % streamClock, alarmID, err_local)) then
            STREAM_ERROR_WRITE('Attempting to add an alarm '//trim(alarmID)//' that does not exist on clock')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Check that the alarm does not already exist for the stream in the specified direction
        !
        if (direction == MPAS_STREAM_INPUT .or. direction == MPAS_STREAM_INPUT_OUTPUT) then
            nullify(new_alarm)
            if (MPAS_stream_list_query(stream % alarmList_in, alarmID, new_alarm, ierr=err_local)) then
                STREAM_ERROR_WRITE('Requested input alarm '//trim(alarmID)//' already on stream '//trim(streamID))
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return
            end if
        end if
        if (direction == MPAS_STREAM_OUTPUT .or. direction == MPAS_STREAM_INPUT_OUTPUT) then
            nullify(new_alarm)
            if (MPAS_stream_list_query(stream % alarmList_out, alarmID, new_alarm, ierr=err_local)) then
                STREAM_ERROR_WRITE('Requested output alarm '//trim(alarmID)//' already on stream '//trim(streamID))
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return
            end if
        end if


        !
        ! Add alarm to alarm to the alarms_in and/or alarms_out list
        ! Add alarm to the alarmList_in and/or alarmList_out list for the field
        !
        if ( threadNum == 0 ) then
           if (direction == MPAS_STREAM_INPUT .or. direction == MPAS_STREAM_INPUT_OUTPUT) then
   
               ! If alarm is not already defined, we need to create a new alarm node
               nullify(new_alarm)
               if (.not. MPAS_stream_list_query(manager % alarms_in, alarmID, new_alarm, ierr=err_local)) then
                   allocate(new_alarm)
                   new_alarm % name = alarmID
                   call MPAS_stream_list_create(new_alarm % streamList, ierr=err_local)
                   if (err_local /= MPAS_STREAM_LIST_NOERR) then
                       if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                       STREAM_ERROR_WRITE('Problems while creating stream list for alarm')
                       return 
                   end if
                   nullify(new_alarm % next)
   
                   call MPAS_stream_list_insert(manager % alarms_in, new_alarm, ierr=err_local)
                   if (err_local /= MPAS_STREAM_LIST_NOERR) then
                       if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                       STREAM_ERROR_WRITE('Problems while adding input alarm to list')
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
                   STREAM_ERROR_WRITE('Problems while adding stream to alarm stream list')
                   return
               end if
   
               ! Add alarm to stream alarm list
               allocate(new_xref)
               new_xref % name = alarmID
               new_xref % xref => new_alarm
               call MPAS_stream_list_insert(stream % alarmList_in, new_xref, ierr=err_local)
               if (err_local /= MPAS_STREAM_LIST_NOERR) then
                   if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                   STREAM_ERROR_WRITE('Problems while adding alarm to stream input alarm list')
                   return
               end if
           end if
   
           if (direction == MPAS_STREAM_OUTPUT .or. direction == MPAS_STREAM_INPUT_OUTPUT) then
   
               ! If alarm is not already defined, we need to create a new alarm node
               nullify(new_alarm)
               if (.not. MPAS_stream_list_query(manager % alarms_out, alarmID, new_alarm, ierr=err_local)) then
                   allocate(new_alarm)
                   new_alarm % name = alarmID
                   call MPAS_stream_list_create(new_alarm % streamList, ierr=err_local)
                   if (err_local /= MPAS_STREAM_LIST_NOERR) then
                       if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                       STREAM_ERROR_WRITE('Problems while creating stream list for alarm')
                       return 
                   end if
                   nullify(new_alarm % next)
   
                   call MPAS_stream_list_insert(manager % alarms_out, new_alarm, ierr=err_local)
                   if (err_local /= MPAS_STREAM_LIST_NOERR) then
                       if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                       STREAM_ERROR_WRITE('Problems while adding output alarm to list')
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
                   STREAM_ERROR_WRITE('Problems while adding stream to alarm stream list')
                   return
               end if
   
               ! Add alarm to stream alarm list
               allocate(new_xref)
               new_xref % name = alarmID
               new_xref % xref => new_alarm
               call MPAS_stream_list_insert(stream % alarmList_out, new_xref, ierr=err_local)
               if (err_local /= MPAS_STREAM_LIST_NOERR) then
                   if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                   STREAM_ERROR_WRITE('Problems while adding alarm to stream output alarm list')
                   return
               end if
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
    !>  NOTE: This routine does not support regular expressions for StreamID
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
        integer :: err_local, threadNum

        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_remove_alarm()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Check that stream exists
        !
        nullify(stream)
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream, ierr=err_local)) then
            STREAM_ERROR_WRITE('Requested stream '//trim(streamID)//' does not exist in stream manager')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        if ( threadNum == 0 ) then
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
               STREAM_ERROR_WRITE('Requested to remove alarm from invalid direction from stream '//trim(streamID))
               return
           end if

           !
           ! Remove stream from alarm's streamList in alarms_in or alarms_out
           !
           if (associated(alarmNode)) then
               call MPAS_stream_list_remove(alarmNode % xref % streamList, streamID, streamNode, ierr=ierr)
           else
               if (direction == MPAS_STREAM_INPUT) then
                   STREAM_ERROR_WRITE('Input alarm '//trim(alarmID)//' does not exist on stream '//trim(streamID))
               else
                   STREAM_ERROR_WRITE('Output alarm '//trim(alarmID)//' does not exist on stream '//trim(streamID))
               end if
               return
           end if 
           if (.not. associated(streamNode)) then
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               STREAM_ERROR_WRITE('Alarm '//trim(alarmID)//' does not have stream '//trim(streamID)//' on its stream list.')
               return
           end if

           !
           ! If the alarm has no associated streams, should we remove it from alarms_in or alarms_out?
           !
           if (MPAS_stream_list_length(alarmNode % xref % streamList) == 0) then
               if (direction == MPAS_STREAM_INPUT) then
                   STREAM_ERROR_WRITE('Input alarm '//trim(alarmID)//' has no associated streams and will be deleted.')
                   call MPAS_stream_list_remove(manager % alarms_in, alarmID, alarmNode, ierr=ierr)
               else
                   STREAM_ERROR_WRITE('Output alarm '//trim(alarmID)//' has no associated streams and will be deleted.')
                   call MPAS_stream_list_remove(manager % alarms_out, alarmID, alarmNode, ierr=ierr)
               end if
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
    !>  'streamID' is provided, only alarms associated with streams that match
    !>  the 'streamID' regular expression will be
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
        integer :: local_ierr, threadNum
        logical :: resetAlarms

        threadNum = mpas_threading_get_thread_num()
        
        if (present(streamID)) then
            STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_reset_alarms() for stream ' // trim(streamID))
        else
            STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_reset_alarms()')
        end if

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR


        !
        ! Check for optional direction argument; default direction is both input and output.
        !
        if (present(direction)) then
            local_direction = direction
        else
            local_direction = MPAS_STREAM_INPUT_OUTPUT
        end if


        if ( threadNum == 0 ) then
            !
            ! Check for optional streamID argument; default is to handle all alarms in the manager.
            !
            nullify(stream)
            if (present(streamID)) then
                resetAlarms = .false.
                do while (MPAS_stream_list_query(manager % streams, streamID, stream, ierr=local_ierr))
                    resetAlarms = .true.
                    if (local_direction == MPAS_STREAM_INPUT .or. local_direction == MPAS_STREAM_INPUT_OUTPUT) then
                        alarm_cursor => stream % alarmList_in % head
                        do while (associated(alarm_cursor))
                            if (mpas_is_alarm_ringing(manager % streamClock, alarm_cursor % name, ierr=local_ierr)) then
                                call mpas_reset_clock_alarm(manager % streamClock, alarm_cursor % name, ierr=local_ierr) 
                            end if
                            alarm_cursor => alarm_cursor % next
                        end do
                    end if

                    if (local_direction == MPAS_STREAM_OUTPUT .or. local_direction == MPAS_STREAM_INPUT_OUTPUT) then
                        alarm_cursor => stream % alarmList_out % head
                        do while (associated(alarm_cursor))
                            if (mpas_is_alarm_ringing(manager % streamClock, alarm_cursor % name, ierr=local_ierr)) then
                                call mpas_reset_clock_alarm(manager % streamClock, alarm_cursor % name, ierr=local_ierr) 
                            end if
                            alarm_cursor => alarm_cursor % next
                        end do
                    end if
                end do

                if ( .not. resetAlarms ) then
                    STREAM_ERROR_WRITE('No stream matching '//trim(streamID)//' exists in stream manager.')
                    if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                    return
                end if
            else

                if (local_direction == MPAS_STREAM_INPUT .or. local_direction == MPAS_STREAM_INPUT_OUTPUT) then
                    alarm_cursor => manager % alarms_in % head
                    do while (associated(alarm_cursor))
                        if (mpas_is_alarm_ringing(manager % streamClock, alarm_cursor % name, ierr=local_ierr)) then
                            call mpas_reset_clock_alarm(manager % streamClock, alarm_cursor % name, ierr=local_ierr) 
                        end if
                        alarm_cursor => alarm_cursor % next
                    end do
                end if


                if (local_direction == MPAS_STREAM_OUTPUT .or. local_direction == MPAS_STREAM_INPUT_OUTPUT) then
                    alarm_cursor => manager % alarms_out % head
                    do while (associated(alarm_cursor))
                        if (mpas_is_alarm_ringing(manager % streamClock, alarm_cursor % name, ierr=local_ierr)) then
                            call mpas_reset_clock_alarm(manager % streamClock, alarm_cursor % name, ierr=local_ierr) 
                        end if
                        alarm_cursor => alarm_cursor % next
                    end do
                end if
            end if
        end if

    end subroutine MPAS_stream_mgr_reset_alarms!}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_ringing_alarms
    !
    !> \brief Test whether any I/O alarms in a stream manager are ringing
    !> \author Michael Duda
    !> \date   30 September 2014
    !> \details
    !>  Tests whether any I/O alarms in a stream manager are ringing; if the optional 
    !>  'streamID' argument is given, only alarms for that stream are tested; if 
    !>  the optional argument 'direction' is given, only alarms for the specified 
    !>  direction are tested. If any of the tested alarms is ringing, the function 
    !>  returns .true.; otherwise, it returns .false..
    !>  Note: This function doesn't support streamID regular expressions
    !
    !-----------------------------------------------------------------------
    logical function MPAS_stream_mgr_ringing_alarms(manager, streamID, direction, ierr) !{{{

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in), optional :: streamID
        integer, intent(in), optional :: direction
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream
        type (MPAS_stream_list_type), pointer :: alarm_cursor
        integer :: local_direction
        integer :: local_ierr, threadNum

        threadNum = mpas_threading_get_thread_num()
        
        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_ringing_alarms()')

        MPAS_stream_mgr_ringing_alarms = .false.

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

        
        if ( threadNum == 0 ) then
           if (local_direction == MPAS_STREAM_INPUT .or. local_direction == MPAS_STREAM_INPUT_OUTPUT) then
               if (associated(stream)) then
                   alarm_cursor => stream % alarmList_in % head
               else
                   alarm_cursor => manager % alarms_in % head
               end if
               do while (associated(alarm_cursor))
                   if (mpas_is_alarm_ringing(manager % streamClock, alarm_cursor % name, ierr=local_ierr)) then
                       MPAS_stream_mgr_ringing_alarms = .true.
                       return
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
                       MPAS_stream_mgr_ringing_alarms = .true.
                       return
                   end if
                   alarm_cursor => alarm_cursor % next
               end do
           end if
        end if

    end function MPAS_stream_mgr_ringing_alarms !}}}


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
    subroutine MPAS_stream_mgr_set_property_int(manager, streamID, propertyName, propertyValue, direction, ierr) !{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_set_property_int'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: propertyName
        integer, intent(in) :: propertyValue
        integer, intent(in), optional :: direction
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local, threadNum
        logical :: setProperties

        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_set_property()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        setProperties = .false.
        nullify(stream_cursor)
        do while (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local))
           setProperties = .true.
           !
           ! Set property
           !
           if ( threadNum == 0 ) then
              select case (propertyName)

                  case (MPAS_STREAM_PROPERTY_PRECISION)
                      stream_cursor % precision = propertyValue

                  case (MPAS_STREAM_PROPERTY_CLOBBER)
                      stream_cursor % clobber_mode = propertyValue

                  case (MPAS_STREAM_PROPERTY_IOTYPE)
                      stream_cursor % io_type = propertyValue

                  case default
                      STREAM_ERROR_WRITE('MPAS_stream_mgr_set_property(): No such property ' COMMA propertyName)
                      STREAM_ERROR_WRITE('    or specified property is not of type integer.')
                      if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
              end select
           end if
        end do

        call mpas_threading_barrier()

        if ( .not. setProperties ) then
            STREAM_ERROR_WRITE('No stream matching '//trim(streamID)//' exists in call to MPAS_stream_mgr_set_property().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if


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
    subroutine MPAS_stream_mgr_set_property_char(manager, streamID, propertyName, propertyValue, direction, ierr) !{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_set_property_char'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: propertyName
        character (len=*), intent(in) :: propertyValue
        integer, intent(in), optional :: direction
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local, threadNum
        logical :: setProperties

        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_set_property()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        setProperties = .false.
        nullify(stream_cursor)
        do while (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local))
            setProperties = .true.
            if ( threadNum == 0 ) then
               !
               ! Set property
               !
               select case (propertyName)

                   case (MPAS_STREAM_PROPERTY_FILENAME)
!TODO: ensure that filename does not contain ':' characters, which PNETCDF does not like...
                       stream_cursor % filename_template = propertyValue

                   case (MPAS_STREAM_PROPERTY_FILENAME_INTV)
                       stream_cursor % filename_interval = propertyValue

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
                       STREAM_ERROR_WRITE('  MPAS_stream_mgr_set_property(): No such property ' COMMA propertyName)
                       STREAM_ERROR_WRITE('    or specified property is not of type character.')
                       if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               end select
            end if
        end do

        call mpas_threading_barrier()

        if ( .not. setProperties ) then
            STREAM_ERROR_WRITE('No stream matching '//trim(streamID)//' exists in call to MPAS_stream_mgr_set_property().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if


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
    subroutine MPAS_stream_mgr_set_property_logical(manager, streamID, propertyName, propertyValue, direction, ierr) !{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_set_property_logical'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: propertyName
        logical, intent(in) :: propertyValue
        integer, intent(in), optional :: direction
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local, threadNum
        logical :: setProperties

        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_set_property()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        setProperties = .false.
        nullify(stream_cursor)
        do while (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local))
            setProperties = .true.

            if ( threadNum == 0 ) then
               !
               ! Set property
               !
               select case (propertyName)

                   case (MPAS_STREAM_PROPERTY_ACTIVE)
                       stream_cursor % active_stream = propertyValue

                   case (MPAS_STREAM_PROPERTY_IMMUTABLE)
                       stream_cursor % immutable = propertyValue

                   case default
                       STREAM_ERROR_WRITE('  MPAS_stream_mgr_set_property(): No such property ' COMMA propertyName)
                       STREAM_ERROR_WRITE('      or specified property is not of type logical.')
                       if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               end select
            end if
        end do

        call mpas_threading_barrier()

        if ( .not. setProperties ) then
            STREAM_ERROR_WRITE('No stream matching '//trim(streamID)//' exists in call to MPAS_stream_mgr_set_property().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if


    end subroutine MPAS_stream_mgr_set_property_logical !}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_get_property_int
    !
    !> \brief Sets a property of a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Retrieves the value of a stream property within an MPAS stream manager.
    !>  NOTE: This routine does not support streamID regular expressions
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_get_property_int(manager, streamID, propertyName, propertyValue, direction, ierr) !{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_set_property_int'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: propertyName
        integer, intent(out) :: propertyValue
        integer, intent(in), optional :: direction
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_get_property()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        nullify(stream_cursor)
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
            STREAM_ERROR_WRITE('Stream '//trim(streamID)//' does not exist in call to MPAS_stream_mgr_get_property().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if

        !
        ! Set property
        !
        select case (propertyName)

            case (MPAS_STREAM_PROPERTY_PRECISION)
                propertyValue = stream_cursor % precision

            case (MPAS_STREAM_PROPERTY_CLOBBER)
                propertyValue = stream_cursor % clobber_mode

            case (MPAS_STREAM_PROPERTY_IOTYPE)
                propertyValue = stream_cursor % io_type

            case default
                STREAM_ERROR_WRITE('MPAS_stream_mgr_get_property(): No such property ' COMMA propertyName)
                STREAM_ERROR_WRITE('    or specified property is not of type integer.')
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
    !>  Retrieves the value of a stream property within an MPAS stream manager.
    !>  NOTE: This routine does not support streamID regular expressions
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_get_property_char(manager, streamID, propertyName, propertyValue, direction, ierr) !{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_get_property_char'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: propertyName
        character (len=*), intent(out) :: propertyValue
        integer, intent(in), optional :: direction
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        type (MPAS_stream_list_type), pointer :: alarm_cursor
        type (MPAS_timeInterval_type) :: temp_interval, interval
        integer :: nAlarms
        integer :: err_local

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_get_property()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        nullify(stream_cursor)
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
            STREAM_ERROR_WRITE('Stream '//trim(streamID)//' does not exist in call to MPAS_stream_mgr_get_property().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if

        !
        ! Set property
        !
        select case (propertyName)

            case (MPAS_STREAM_PROPERTY_FILENAME)
                propertyValue = stream_cursor % filename_template

            case (MPAS_STREAM_PROPERTY_FILENAME_INTV)
                propertyValue = stream_cursor % filename_interval

            case (MPAS_STREAM_PROPERTY_REF_TIME)
                call mpas_get_time(stream_cursor % referenceTime, dateTimeString=propertyValue)

            case (MPAS_STREAM_PROPERTY_RECORD_INTV)

                ! The interval between records may not have been allocated if the optional recordInterval
                !    argument was not provided when the stream was created. If there is no explicit recordInterval,
                !    assume that the interval is the shortest interval between alarms on the stream; since
                !    recordInterval is only used for reading, use the input alarm list in this check.
                if (.not. associated(stream_cursor % recordInterval)) then
            
                    !
                    ! If no direction is specified, return the read interval, since this was the only historic 
                    !    use of the recordInterval for a stream.
                    !
                    if (present(direction)) then
                        if (direction == MPAS_STREAM_OUTPUT) then
                            alarm_cursor => stream_cursor % alarmList_out % head
                        else
                            alarm_cursor => stream_cursor % alarmList_in % head
                        end if
                    else
                        alarm_cursor => stream_cursor % alarmList_in % head
                    end if
                    nAlarms = 0
                    do while (associated(alarm_cursor))
                        temp_interval = mpas_alarm_interval(manager % streamClock, alarm_cursor % name, err_local)
                        if (err_local == 0) then
                            if (nAlarms == 0) then
                                interval = temp_interval
                            else if (temp_interval < interval) then
                                interval = temp_interval
                            end if
                            nAlarms = nAlarms + 1
                        end if
                        alarm_cursor => alarm_cursor % next
                    end do
                    if (nAlarms > 0) then
                        call mpas_get_timeInterval(interval, timeString=propertyValue)
                    else
                        propertyValue = 'none'
                    end if
                else
                    call mpas_get_timeInterval(stream_cursor % recordInterval, timeString=propertyValue)
                end if

            case default
                STREAM_ERROR_WRITE(' MPAS_stream_mgr_get_property(): No such property ' COMMA propertyName)
                STREAM_ERROR_WRITE('   or specified property is not of type character.')
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
    !>  Retrieves the value of a stream property within an MPAS stream manager.
    !>  NOTE: This routine does not support streamID regular expressions
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_get_property_logical(manager, streamID, propertyName, propertyValue, direction, ierr) !{{{

        implicit none

        character (len=*), parameter :: sub = 'MPAS_stream_mgr_get_property'

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: propertyName
        logical, intent(out) :: propertyValue
        integer, intent(in), optional :: direction
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: err_local

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_get_property_logical()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        nullify(stream_cursor)
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local)) then
            STREAM_ERROR_WRITE('Stream '//trim(streamID)//' does not exist in call to MPAS_stream_mgr_get_property().')
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
                STREAM_ERROR_WRITE('MPAS_stream_mgr_get_property(): No such property ' COMMA propertyName)
                STREAM_ERROR_WRITE('    or specified property is not of type logical.')
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
        logical :: addedPackages


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_pkg()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Query pointer to package in the manager-wide package pool
        !
        nullify(package)
        call mpas_pool_get_package(manager % allPackages, packageName, package)
        if (.not. associated(package)) then
            STREAM_ERROR_WRITE('Package '//trim(packageName)//' not found in call to MPAS_stream_mgr_add_pkg().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return
        end if

        !
        ! Find requested stream
        !
        addedPackages = .false.
        nullify(stream_cursor)
        do while (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local))

           addedPackages = .true.

           !
           ! Add package to the packages pool for the stream
           !
           call mpas_pool_add_package(stream_cursor % pkg_pool, packageName, package)    

        end do

        if ( .not. addedPackages ) then
            STREAM_ERROR_WRITE('No stream matching '//trim(streamID)//' exists in call to MPAS_stream_mgr_add_pkg().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if


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
        logical :: removedPackage


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_remove_pkg()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR

        !
        ! Find requested stream
        !
        removedPackage = .false.
        nullify(stream_cursor)
        do while (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local))

            removedPackage = .true.

            !
            ! Remove package from the packages pool for the stream
            !
            call mpas_pool_remove_package(stream_cursor % pkg_pool, packageName)

        end do

        if ( .not. removedPackage ) then
            STREAM_ERROR_WRITE('No stream matching '//trim(streamID)//' exists in call to MPAS_stream_mgr_remove_pkg().')
            if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
            return 
        end if


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
        type (MPAS_pool_type), pointer :: att_pool
        integer, pointer :: queryVal
        integer :: att_type
        integer :: err_level
        integer :: err_local
        logical :: addedAttribute


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_att()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR
        nullify(queryVal)

        if (present(streamID)) then
            addedAttribute = .false.
            nullify(stream_cursor)
            do while (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local))
                addedAttribute = .true.
                att_pool => stream_cursor % att_pool
                err_level = mpas_pool_get_error_level()
                call mpas_pool_set_error_level(MPAS_POOL_SILENT)
                call mpas_pool_get_config(att_pool, attName, queryVal)
                call mpas_pool_set_error_level(err_level)
                if (.not. associated(queryVal)) then
                    !
                    ! Querying the type of the attribute should return MPAS_POOL_FATAL if the attribute really
                    !    does not exist in the pool; otherwise, the attribute exists but was of the wrong type
                    !    in the call above to mpas_pool_get_config()
                    !
                    if (mpas_pool_config_type(att_pool, attName) /= MPAS_POOL_FATAL) then
                        STREAM_ERROR_WRITE('Attribute '//trim(attName)//' in stream '//trim(stream_cursor % name)//' is not of type integer.')
                        if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                        return 
                    end if
                    call mpas_pool_add_config(att_pool, attName, attVal)
                else
                    queryVal = attVal
                end if
            end do

            if ( .not. addedAttribute ) then
                STREAM_ERROR_WRITE('No stream matching '//trim(streamID)//' exists in call to MPAS_stream_mgr_add_att().')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return 
            end if

        else
            att_pool => manager % defaultAtts
            err_level = mpas_pool_get_error_level()
            call mpas_pool_set_error_level(MPAS_POOL_SILENT)
            call mpas_pool_get_config(att_pool, attName, queryVal)
            call mpas_pool_set_error_level(err_level)
            if (.not. associated(queryVal)) then
                !
                ! Querying the type of the attribute should return MPAS_POOL_FATAL if the attribute really
                !    does not exist in the pool; otherwise, the attribute exists but was of the wrong type
                !    in the call above to mpas_pool_get_config()
                !
                if (mpas_pool_config_type(att_pool, attName) /= MPAS_POOL_FATAL) then
                    STREAM_ERROR_WRITE('Attribute '//trim(attName)//' in streamManager is not of type integer.')
                    if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                    return 
                end if
                call mpas_pool_add_config(att_pool, attName, attVal)
            else
                queryVal = attVal
            end if
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
        type (MPAS_pool_type), pointer :: att_pool
        real (kind=RKIND), pointer :: queryVal
        integer :: att_type
        integer :: err_level
        integer :: err_local
        logical :: addedAttribute


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_att()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR
        nullify(queryVal)

        if (present(streamID)) then
            addedAttribute = .false.
            nullify(stream_cursor)
            do while (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local))

                addedAttribute = .true.
                att_pool => stream_cursor % att_pool

                err_level = mpas_pool_get_error_level()
                call mpas_pool_set_error_level(MPAS_POOL_SILENT)
                call mpas_pool_get_config(att_pool, attName, queryVal)
                call mpas_pool_set_error_level(err_level)
                if (.not. associated(queryVal)) then
                    !
                    ! Querying the type of the attribute should return MPAS_POOL_FATAL if the attribute really
                    !    does not exist in the pool; otherwise, the attribute exists but was of the wrong type
                    !    in the call above to mpas_pool_get_config()
                    !
                    if (mpas_pool_config_type(att_pool, attName) /= MPAS_POOL_FATAL) then
                        STREAM_ERROR_WRITE('Attribute '//trim(attName)//' in stream '//trim(stream_cursor % name)//' is not of type real.')
                        if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                        return 
                    end if
                    call mpas_pool_add_config(att_pool, attName, attVal)
                else
                    queryVal = attVal
                end if
            end do

            if ( .not. addedAttribute ) then
                STREAM_ERROR_WRITE('No stream matching '//trim(streamID)//' exists in call to MPAS_stream_mgr_add_att().')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return 
            end if

        else
            att_pool => manager % defaultAtts

            err_level = mpas_pool_get_error_level()
            call mpas_pool_set_error_level(MPAS_POOL_SILENT)
            call mpas_pool_get_config(att_pool, attName, queryVal)
            call mpas_pool_set_error_level(err_level)
            if (.not. associated(queryVal)) then
                !
                ! Querying the type of the attribute should return MPAS_POOL_FATAL if the attribute really
                !    does not exist in the pool; otherwise, the attribute exists but was of the wrong type
                !    in the call above to mpas_pool_get_config()
                !
                if (mpas_pool_config_type(att_pool, attName) /= MPAS_POOL_FATAL) then
                    STREAM_ERROR_WRITE('Attribute '//trim(attName)//' in streamManager is not of type real.')
                    if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                    return 
                end if
                call mpas_pool_add_config(att_pool, attName, attVal)
            else
                queryVal = attVal
            end if

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
        type (MPAS_pool_type), pointer :: att_pool
        character (len=StrKIND), pointer :: queryVal
        integer :: att_type
        integer :: err_level
        integer :: err_local
        logical :: addedAttribute


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_att()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR
        nullify(queryVal)

        if (present(streamID)) then
            addedAttribute = .false.
            nullify(stream_cursor)
            do while (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local))
                addedAttribute = .true.
                att_pool => stream_cursor % att_pool

                err_level = mpas_pool_get_error_level()
                call mpas_pool_set_error_level(MPAS_POOL_SILENT)
                call mpas_pool_get_config(att_pool, attName, queryVal)
                call mpas_pool_set_error_level(err_level)
                if (.not. associated(queryVal)) then
                    !
                    ! Querying the type of the attribute should return MPAS_POOL_FATAL if the attribute really
                    !    does not exist in the pool; otherwise, the attribute exists but was of the wrong type
                    !    in the call above to mpas_pool_get_config()
                    !
                    if (mpas_pool_config_type(att_pool, attName) /= MPAS_POOL_FATAL) then
                        STREAM_ERROR_WRITE('Attribute '//trim(attName)//' in stream '//trim(stream_cursor % name)//' is not of type character.')
                        if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                        return 
                    end if
                    call mpas_pool_add_config(att_pool, attName, attVal)
                else
                    queryVal = attVal
                end if
            end do

            if ( .not. addedAttribute ) then
                STREAM_ERROR_WRITE('No stream matching '//trim(streamID)//' exists in call to MPAS_stream_mgr_add_att().')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return 
            end if

        else
            att_pool => manager % defaultAtts

            err_level = mpas_pool_get_error_level()
            call mpas_pool_set_error_level(MPAS_POOL_SILENT)
            call mpas_pool_get_config(att_pool, attName, queryVal)
            call mpas_pool_set_error_level(err_level)
            if (.not. associated(queryVal)) then
                !
                ! Querying the type of the attribute should return MPAS_POOL_FATAL if the attribute really
                !    does not exist in the pool; otherwise, the attribute exists but was of the wrong type
                !    in the call above to mpas_pool_get_config()
                !
                if (mpas_pool_config_type(att_pool, attName) /= MPAS_POOL_FATAL) then
                    STREAM_ERROR_WRITE('Attribute '//trim(attName)//' in streamManager is not of type character.')
                    if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                    return 
                end if
                call mpas_pool_add_config(att_pool, attName, attVal)
            else
                queryVal = attVal
            end if

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
        type (MPAS_pool_type), pointer :: att_pool
        logical, pointer :: queryVal
        integer :: att_type
        integer :: err_level
        integer :: err_local
        logical :: addedAttribute


        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_att()')

        if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR
        nullify(queryVal)

        if (present(streamID)) then
            addedAttribute = .false.
            nullify(stream_cursor)
            do while (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=err_local))
                addedAttribute = .true.
                att_pool => stream_cursor % att_pool

                err_level = mpas_pool_get_error_level()
                call mpas_pool_set_error_level(MPAS_POOL_SILENT)
                call mpas_pool_get_config(att_pool, attName, queryVal)
                call mpas_pool_set_error_level(err_level)
                if (.not. associated(queryVal)) then
                    !
                    ! Querying the type of the attribute should return MPAS_POOL_FATAL if the attribute really
                    !    does not exist in the pool; otherwise, the attribute exists but was of the wrong type
                    !    in the call above to mpas_pool_get_config()
                    !
                    if (mpas_pool_config_type(att_pool, attName) /= MPAS_POOL_FATAL) then
                        STREAM_ERROR_WRITE('Attribute '//trim(attName)//' in stream '//trim(stream_cursor % name)//' is not of type logical.')
                        if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                        return 
                    end if
                    call mpas_pool_add_config(att_pool, attName, attVal)
                else
                    queryVal = attVal
                end if
            end do

            if ( .not. addedAttribute ) then
                STREAM_ERROR_WRITE('No stream matching '//trim(streamID)//' exists in call to MPAS_stream_mgr_add_att().')
                if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                return 
            end if

        else
            att_pool => manager % defaultAtts

            err_level = mpas_pool_get_error_level()
            call mpas_pool_set_error_level(MPAS_POOL_SILENT)
            call mpas_pool_get_config(att_pool, attName, queryVal)
            call mpas_pool_set_error_level(err_level)
            if (.not. associated(queryVal)) then
                !
                ! Querying the type of the attribute should return MPAS_POOL_FATAL if the attribute really
                !    does not exist in the pool; otherwise, the attribute exists but was of the wrong type
                !    in the call above to mpas_pool_get_config()
                !
                if (mpas_pool_config_type(att_pool, attName) /= MPAS_POOL_FATAL) then
                    STREAM_ERROR_WRITE('Attribute '//trim(attName)//' in streamManger is not of type logical.')
                    if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                    return 
                end if
                call mpas_pool_add_config(att_pool, attName, attVal)
            else
                queryVal = attVal
            end if
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
    !>  The "streamID" argument optionally specifies a regular expression of streams
    !>  to be written; if no other optional arguments are given, the specified
    !>  stream is only written if any of its alarms are ringing.
    !>  The "timeLevel" argument optionally specifies, for fields with multiple
    !>  time levels, the time level from which fields should be written.
    !>  The "mgLevel" argument optionally specifies, for fields that exist for 
    !>  multiple grid levels, the grid level from which fields should be written.
    !>  The "forceWriteNow" argument optionally specifies that all streams -- or 
    !>  the stream specified by the "streamID" argument -- should be written by 
    !>  the call regardless of whether any alarms associated with the stream(s) 
    !>  are ringing. The "writeTime" argument optionally specifies a time stamp
    !>  to be used for expanding a filename template, when it is not passed in,
    !>  the current time of the streamManager's clock is used to expand filename
    !>  templates.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_write(manager, streamID, timeLevel, mgLevel, forceWriteNow, writeTime, ierr) !{{{

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in), optional :: streamID
        integer, intent(in), optional :: timeLevel
        integer, intent(in), optional :: mgLevel
        logical, intent(in), optional :: forceWriteNow
        character (len=*), intent(in), optional :: writeTime
        integer, intent(out), optional :: ierr

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: blockID
        integer :: local_timeLevel
        integer :: local_mgLevel
        logical :: local_forceWrite
        integer :: local_ierr
        type (MPAS_Time_type) :: local_writeTime
        integer :: temp_ierr
        integer :: threadNum
        logical :: wroteStreams

        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_write()')
        local_ierr = MPAS_STREAM_MGR_NOERR

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


        if ( present(writeTime) ) then
            call mpas_set_time(local_writeTime, dateTimeString=writeTime, ierr=ierr)
        else
            local_writeTime = mpas_get_clock_time(manager % streamClock, MPAS_NOW, ierr=local_ierr)
        end if

        ! This routine assumes global writes, so set blockID = -1
        blockID = -1

        if ( threadNum == 0 ) then
           !
           ! If a stream is specified, we process just that stream; otherwise,
           !    process all streams
           !
           if (present(streamID)) then
               wroteStreams = .false.
               nullify(stream_cursor)
               do while (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=ierr))
                   STREAM_DEBUG_WRITE('-- Handling write of stream '//trim(stream_cursor % name))

                   ! Verify that the stream is an output stream
                   if (stream_cursor % direction == MPAS_STREAM_OUTPUT .or. &
                       stream_cursor % direction == MPAS_STREAM_INPUT_OUTPUT) then

                      wroteStreams = .true.
                      stream_cursor % blockWrite = .false.
                      call write_stream(manager, stream_cursor, blockID, local_timeLevel, local_mgLevel, local_forceWrite, local_writeTime, local_ierr)
                   end if
               end do

               if ( .not. wroteStreams ) then
                   STREAM_ERROR_WRITE('No output stream matching '//trim(streamID)//' exists in call to MPAS_stream_mgr_write().')
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

                       stream_cursor % blockWrite = .false.
                       call write_stream(manager, stream_cursor, blockID, local_timeLevel, local_mgLevel, local_forceWrite, local_writeTime, temp_ierr)
                       if (temp_ierr /= MPAS_STREAM_MGR_NOERR) then
                           local_ierr = temp_ierr
                       end if

                   end if
                   stream_cursor => stream_cursor % next
               end do
           end if
        end if

        if (present(ierr)) ierr = local_ierr

    end subroutine MPAS_stream_mgr_write !}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_block_write
    !
    !> \brief Write streams that are managed by an MPAS stream manager, using blockWrite = .true.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   08/25/2016
    !> \details
    !>  This routine performs block writes of any streams requested within an
    !>  MPAS stream manger. Block writes imply each processor should write
    !>  their own file(s). Additionally, block writes disable global reindexing
    !>  before and after writing the stream.
    !>
    !>  With no optional arguments, writes all streams whose alarms are ringing.
    !>  The "streamID" argument optionally specifies a regular expression of streams
    !>  to be written; if no other optional arguments are given, the specified
    !>  stream is only written if any of its alarms are ringing.
    !>  The "timeLevel" argument optionally specifies, for fields with multiple
    !>  time levels, the time level from which fields should be written.
    !>  The "mgLevel" argument optionally specifies, for fields that exist for 
    !>  multiple grid levels, the grid level from which fields should be written.
    !>  The "forceWriteNow" argument optionally specifies that all streams -- or 
    !>  the stream specified by the "streamID" argument -- should be written by 
    !>  the call regardless of whether any alarms associated with the stream(s) 
    !>  are ringing. The "writeTime" argument optionally specifies a time stamp
    !>  to be used for expanding a filename template, when it is not passed in,
    !>  the current time of the streamManager's clock is used to expand filename
    !>  templates.
    !>
    !>  Swapping write types on a stream can cause unexpected behavior. For example:
    !>    - Writing a stream using the normal write routine, followed by the
    !>          block write routine, and again by the normal write routine will cause any
    !>          data in the files prior to the last normal write to be lost.
    !>
    !>  Additionally, this routine is only valid for writing out a single file
    !>  at a time.  If the block write stream is required at multiple time
    !>  levels, the filename template and filename intervals need to be setup to
    !>  give unique files for each record.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_block_write(manager, writeBlock, streamID, timeLevel, mgLevel, forceWriteNow, writeTime, ierr) !{{{

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        type (block_type), intent(in) :: writeBlock
        character (len=*), intent(in), optional :: streamID
        integer, intent(in), optional :: timeLevel
        integer, intent(in), optional :: mgLevel
        logical, intent(in), optional :: forceWriteNow
        character (len=*), intent(in), optional :: writeTime
        integer, intent(out), optional :: ierr

        type (mpas_pool_type), pointer :: backupAllFields, backupAllStructs

        type (MPAS_stream_list_type), pointer :: stream_cursor
        integer :: local_timeLevel
        integer :: local_mgLevel
        logical :: local_forceWrite
        integer :: local_ierr
        type (MPAS_Time_type) :: local_writeTime
        integer :: temp_ierr
        integer :: threadNum
        logical :: wroteStreams

        type (mpas_io_context_type), pointer :: debugContext, backupContext

        threadNum = mpas_threading_get_thread_num()

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_block_write()')
        local_ierr = MPAS_STREAM_MGR_NOERR

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


        if ( present(writeTime) ) then
            call mpas_set_time(local_writeTime, dateTimeString=writeTime, ierr=ierr)
        else
            local_writeTime = mpas_get_clock_time(manager % streamClock, MPAS_NOW, ierr=local_ierr)
        end if


        if ( threadNum == 0 ) then
           ! Setup the debugging ioContext
           allocate(debugContext)
           allocate(debugContext % dminfo)
           call mpas_dmpar_init(debugContext % dminfo, MPI_COMM_SELF)

           backupContext => manager % ioContext
           manager % ioContext => debugContext

           call mpas_io_init(debugContext, 1, 1)

           ! Setup the backup pools, and swap the manager's allFields / allStructs pools
           backupAllFields => manager % allFields
           backupAllStructs => manager % allStructs
           manager % allFields => writeBlock % allFields
           manager % allStructs => writeBlock % allStructs

           !
           ! If a stream is specified, we process just that stream; otherwise,
           !    process all streams
           !
           if (present(streamID)) then
               wroteStreams = .false.
               nullify(stream_cursor)
               do while (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=ierr))
                   STREAM_DEBUG_WRITE('-- Handling local write of stream '//trim(stream_cursor % name))

                   ! Verify that the stream is an output stream
                   if (stream_cursor % direction == MPAS_STREAM_OUTPUT .or. &
                       stream_cursor % direction == MPAS_STREAM_INPUT_OUTPUT) then

                      wroteStreams = .true.
                      stream_cursor % blockWrite = .true.
                      call write_stream(manager, stream_cursor, writeBlock % blockID, local_timeLevel, local_mgLevel, &
                                        local_forceWrite, local_writeTime, local_ierr)
                      stream_cursor % blockWrite = .false.
                      if ( associated(stream_cursor % stream) ) then
                         stream_cursor % valid = .false.
                         call MPAS_closeStream(stream_cursor % stream, ierr=local_ierr)
                         deallocate(stream_cursor % stream)
                      end if
                   end if
               end do

               if ( .not. wroteStreams ) then
                   STREAM_ERROR_WRITE('No output stream matching '//trim(streamID)//' exists in call to MPAS_stream_mgr_block_write().')
                   if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
                   return 
               end if
           else
               nullify(stream_cursor)
               stream_cursor => manager % streams % head
               do while (associated(stream_cursor))
                   STREAM_DEBUG_WRITE('-- Handling local write of stream '//trim(stream_cursor % name))

                   ! Verify that the stream is an output stream
                   if (stream_cursor % direction == MPAS_STREAM_OUTPUT .or. &
                       stream_cursor % direction == MPAS_STREAM_INPUT_OUTPUT) then

                       stream_cursor % blockWrite = .true.
                       call write_stream(manager, stream_cursor, writeBlock % blockID, local_timeLevel, local_mgLevel, &
                                         local_forceWrite, local_writeTime, temp_ierr)
                       stream_cursor % blockWrite = .false.
                       if ( associated(stream_cursor % stream) ) then
                          stream_cursor % valid = .false.
                          call MPAS_closeStream(stream_cursor % stream, ierr=local_ierr)
                          deallocate(stream_cursor % stream)
                       end if
                       if (temp_ierr /= MPAS_STREAM_MGR_NOERR) then
                           local_ierr = temp_ierr
                       end if

                   end if
                   stream_cursor => stream_cursor % next
               end do
           end if
           ! Swap the allFields / allStructs back to the backups.
           manager % allFields => backupAllFields
           manager % allStructs => backupAllStructs

           ! Destroy the IO context, and dminfo that was created
           call MPAS_io_finalize(debugContext, .true.)

           manager % ioContext => backupContext

           call mpas_dmpar_finalize(debugContext % dminfo)

           deallocate(debugContext % dminfo)
           deallocate(debugContext) 
        end if

        call mpas_threading_barrier()

        if (present(ierr)) ierr = local_ierr

    end subroutine MPAS_stream_mgr_block_write !}}}


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
    subroutine write_stream(manager, stream, blockID, timeLevel, mgLevel, forceWritenow, writeTime, ierr) !{{{

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        type (MPAS_stream_list_type), intent(inout) :: stream
        integer, intent(in) :: blockID
        integer, intent(in) :: timeLevel
        integer, intent(in) :: mgLevel
        logical, intent(in) :: forceWriteNow
        type (MPAS_Time_type), intent(in) :: writeTime
        integer, intent(out) :: ierr

        type (MPAS_stream_list_type), pointer :: alarm_cursor
        type (MPAS_Time_type) :: ref_time
        type (MPAS_TimeInterval_type) :: temp_interval
        type (MPAS_TimeInterval_type) :: filename_interval
        character (len=StrKIND) :: now_string, time_string
        character (len=StrKIND) :: temp_filename, actualWhen
        character (len=StrKIND) :: err_string
        logical :: ringing_alarm, recordSeek, swapRecords
        logical :: clobberRecords, clobberFiles, truncateFiles
        integer :: maxRecords, tempRecord
        integer :: local_ierr, threadNum

        threadNum = mpas_threading_get_thread_num()

        ierr = MPAS_STREAM_MGR_NOERR
        if ( threadNum == 0 ) then
           swapRecords = .false.

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
           ! Work out file clobbering options
           !
           if (stream % clobber_mode == MPAS_STREAM_CLOBBER_OVERWRITE) then
               clobberRecords = .true.
           else
               clobberRecords = .false.
           end if

           if (stream % clobber_mode == MPAS_STREAM_CLOBBER_OVERWRITE .or. &
               stream % clobber_mode == MPAS_STREAM_CLOBBER_TRUNCATE .or. &
               stream % clobber_mode == MPAS_STREAM_CLOBBER_APPEND) then
               clobberFiles = .true.
           else
               clobberFiles = .false.
           end if

           if (stream % clobber_mode == MPAS_STREAM_CLOBBER_TRUNCATE) then
               truncateFiles = .true.
           else
               truncateFiles = .false.
           end if

           !
           ! If the stream is not valid, assume that we have not yet written this
           ! stream, in which case we create the stream from scratch
           !
           if (.not. stream % valid) then
               if ( stream % filename_interval /= 'none' ) then
                   call mpas_set_timeInterval(filename_interval, timeString=stream % filename_interval)
                   call mpas_build_stream_filename(stream % referenceTime, writeTime, filename_interval, stream % filename_template, blockID, stream % filename, ierr=local_ierr)
               else
                   call mpas_get_time(stream % referenceTime, dateTimeString=time_string)
                   call mpas_expand_string(time_string, blockID, stream % filename_template, stream % filename)
               end if

               stream % nRecords = 1

               recordSeek = .false.
               ! Based on clobber_mode, determine if it matters if the file exists or not.
               if ( stream % clobber_mode == MPAS_STREAM_CLOBBER_OVERWRITE .or. stream % clobber_mode == MPAS_STREAM_CLOBBER_APPEND ) then
                   STREAM_DEBUG_WRITE(' -- Cobber mode is overwrite or append...')
    
                   ! Check if the file exists
                   inquire(file=trim(stream % filename), exist=recordSeek)
               end if

               !
               ! Build stream from pools of fields and attributes
               !
               allocate(stream % stream)
               call MPAS_createStream(stream % stream, manager % ioContext, stream % filename, stream % io_type, MPAS_IO_WRITE,  &
                                      precision=stream % precision, clobberRecords=clobberRecords, &
                                      clobberFiles=clobberFiles, truncateFiles=truncateFiles, ierr=local_ierr)
               stream % stream % blockWrite = stream % blockWrite
               if (local_ierr /= MPAS_STREAM_NOERR) then
                   if (local_ierr == MPAS_STREAM_CLOBBER_FILE) then
                       !
                       ! We should have only reached this point if clobber_mode =  never_modify
                       !
                       write(err_string,'(a)') 'Writing to stream '''//trim(stream % name)//''' would clobber file '''//&
                                               trim(stream % filename)//''','
                       STREAM_ERROR_WRITE(trim(err_string))
                       write(err_string,'(a)') '    but clobber_mode is set to ''never_modify''.'
                       STREAM_ERROR_WRITE(trim(err_string))
                       ierr = MPAS_STREAM_MGR_ERR_CLOBBER_FILE
                   else
                       ierr = MPAS_STREAM_MGR_ERROR
                   end if
                   return
               end if

               ! File exists on disk, prior to creating stream. Need to seek the record to ensure we're writing to the correct place.
               if ( recordSeek ) then
                   STREAM_DEBUG_WRITE(' -- File exists on disk: ' COMMA trim(stream % filename))
                   call mpas_get_time(writeTime, dateTimeString=now_string)
    
                   ! Look for exact record (in the case of overwriting)
                   ! This also gets the number of records in the file.
                   stream % nRecords = MPAS_seekStream(stream % stream, now_string, MPAS_STREAM_EXACT_TIME, actualWhen, maxRecords, local_ierr)
                   STREAM_DEBUG_WRITE(' -- Seeked record is: ' COMMA stream % nRecords COMMA ' with current records equal to ' COMMA maxRecords COMMA ' and an error of ' COMMA local_ierr)
    
                   if ( stream % nRecords == 0 ) then
                       ! If we didn't find an exact time, set record to point to the end of the file.
                       ! This might result in non-monotonic timestamps in the output file.
                       stream % nRecords = maxRecords + 1
                       STREAM_DEBUG_WRITE(' -- No exact time match found for ' COMMA trim(now_string) COMMA ' appending record instead.')
                       STREAM_DEBUG_WRITE(' -- Setting record to: ' COMMA stream % nRecords)
                   end if
               end if

               call build_stream(stream, MPAS_STREAM_OUTPUT, manager % allFields, manager % allPackages, timeLevel, mgLevel, local_ierr)
               if (local_ierr /= MPAS_STREAM_NOERR) then
                   ierr = MPAS_STREAM_MGR_ERROR
                   return
               end if
               stream % timeLevel = timeLevel

               stream % valid = .true.
           else
               if ( stream % filename_interval /= 'none' ) then
                   call mpas_set_timeInterval(filename_interval, timeString=stream % filename_interval)
                   call mpas_build_stream_filename(stream % referenceTime, writeTime, filename_interval, stream % filename_template, blockID, temp_filename, ierr=local_ierr)
               else
                   call mpas_get_time(stream % referenceTime, dateTimeString=time_string)
                   call mpas_expand_string(time_string, blockID, stream % filename_template, temp_filename)
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

                   recordSeek = .false.
                   ! Based on clobber_mode, determine if it matters if the file exists or not.
                   if ( stream % clobber_mode == MPAS_STREAM_CLOBBER_OVERWRITE .or. stream % clobber_mode == MPAS_STREAM_CLOBBER_APPEND ) then
                       STREAM_DEBUG_WRITE(' -- Cobber mode is overwrite or append...')
           
                       ! Check if the file exists
                       inquire(file=trim(stream % filename), exist=recordSeek)
                   end if

                   stream % nRecords = 1

                   !
                   ! Build new stream from pools of fields and attributes
                   !
                   call MPAS_createStream(stream % stream, manager % ioContext, stream % filename, stream % io_type, MPAS_IO_WRITE, &
                                          precision=stream % precision, clobberRecords=clobberRecords, &
                                          clobberFiles=clobberFiles, truncateFiles=truncateFiles, ierr=local_ierr)
                   stream % stream % blockWrite = stream % blockWrite
                   if (local_ierr /= MPAS_STREAM_NOERR) then
                       if (local_ierr == MPAS_STREAM_CLOBBER_FILE) then
                           !
                           ! We should have only reached this point if clobber_mode =  never_modify
                           !
                           write(err_string,'(a)') 'Writing to stream '''//trim(stream % name)//''' would clobber file '''//&
                                                   trim(stream % filename)//''','
                           STREAM_ERROR_WRITE(trim(err_string))
                           write(err_string,'(a)') '    but clobber_mode is set to ''never_modify''.'
                           STREAM_ERROR_WRITE(trim(err_string))
                           ierr = MPAS_STREAM_MGR_ERR_CLOBBER_FILE
                       else
                           ierr = MPAS_STREAM_MGR_ERROR
                       end if
                       stream % valid = .false.
                       return
                   end if

                   ! File exists on disk, prior to creating stream. Need to seek the record to ensure we're writing to the correct place.
                   if ( recordSeek ) then
                       STREAM_DEBUG_WRITE(' -- File exists on disk: ' COMMA trim(stream % filename))
                       call mpas_get_time(writeTime, dateTimeString=now_string)
           
                       ! Look for exact record (in the case of overwriting)
                       ! This also gets the number of records in the file.
                       stream % nRecords = MPAS_seekStream(stream % stream, now_string, MPAS_STREAM_EXACT_TIME, actualWhen, maxRecords, local_ierr)
                       STREAM_DEBUG_WRITE(' -- Seeked record is: ' COMMA stream % nRecords COMMA ' with current records equal to ' COMMA maxRecords COMMA ' and an error of ' COMMA local_ierr)
           
                       if ( stream % nRecords == 0 ) then
                           ! If we didn't find an exact time, set record to point to the end of the file.
                           ! This might result in non-monotonic timestamps in the output file.
                           stream % nRecords = maxRecords + 1
                           STREAM_DEBUG_WRITE(' -- No exact time match found for ' COMMA trim(now_string) COMMA ' appending record instead.')
                           STREAM_DEBUG_WRITE(' -- Setting record to: ' COMMA stream % nRecords)
                       end if
                   end if

                   call build_stream(stream, MPAS_STREAM_OUTPUT, manager % allFields, manager % allPackages, timeLevel, mgLevel, local_ierr)
                   if (local_ierr /= MPAS_STREAM_NOERR) then
                       ierr = MPAS_STREAM_MGR_ERROR
                       return
                   end if
                   stream % timeLevel = timeLevel
               else
                   stream % nRecords = stream % nRecords + 1
                   if ( stream % clobber_mode == MPAS_STREAM_CLOBBER_OVERWRITE .or. stream % clobber_mode == MPAS_STREAM_CLOBBER_APPEND ) then
                       call mpas_get_time(writeTime, dateTimeString=now_string)
           
                       ! Look for exact record (in the case of overwriting)
                       ! This also gets the number of records in the file.
                       tempRecord = MPAS_seekStream(stream % stream, now_string, MPAS_STREAM_EXACT_TIME, actualWhen, maxRecords, local_ierr)
                       STREAM_DEBUG_WRITE(' -- Seeked record is: ' COMMA tempRecord COMMA ' with current records equal to ' COMMA maxRecords COMMA ' and an error of ' COMMA local_ierr)
           
                       if ( tempRecord /= 0 .and. stream % nRecords < maxRecords ) then
                           ! If we found an exact result
                           ! This might result in non-monotonic timestamps in the output file.
                           swapRecords = .true.
                           maxRecords = stream % nRecords
                           stream % nRecords = tempRecord
                           tempRecord = maxRecords
                           STREAM_DEBUG_WRITE(' -- Exact time match found for ' COMMA trim(now_string) )
                           STREAM_DEBUG_WRITE(' -- Setting record to: ' COMMA stream % nRecords)
                       else if ( tempRecord == 0 .and. stream % nRecords < maxRecords ) then
                           ! If we didn't find an exact time, set record to point to the end of the file.
                           ! This might result in non-monotonic timestamps in the output file.
                           stream % nRecords = maxRecords + 1
                           STREAM_DEBUG_WRITE(' -- No exact time match found for ' COMMA trim(now_string) COMMA ' appending record instead.')
                           STREAM_DEBUG_WRITE(' -- Setting record to: ' COMMA stream % nRecords)
                       end if
                   end if
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
           ! unless this stream is using blockWrite
           !
           if ( .not. stream % blockWrite ) then
              STREAM_DEBUG_WRITE(' -- Prewrite reindex for stream ' COMMA trim(stream % name))
              call prewrite_reindex(manager % allFields, stream % field_pool)
           end if

           ! 
           ! Write the stream
           ! 
           STREAM_DEBUG_WRITE(' -- Writing stream ' COMMA trim(stream % name))
           call MPAS_writeStream(stream % stream, stream % nRecords, ierr=local_ierr)

           !
           ! Regardless of the error code from MPAS_writeStream, we need to reset global indices in the stream back to local indices
           ! unless this stream is using blockWrite
           !
           if ( .not. stream % blockWrite ) then
              STREAM_DEBUG_WRITE(' -- Postwrite reindex for stream ' COMMA trim(stream % name))
              call postwrite_reindex(manager % allFields, stream % field_pool)
           end if

           if (local_ierr /= MPAS_STREAM_NOERR) then
               if (local_ierr == MPAS_STREAM_CLOBBER_RECORD) then
                   !
                   ! We should have only reached this point if clobber_mode =  append
                   !
                   write(err_string,'(a,i4,a)') 'Writing to stream '''//trim(stream % name)//''' would overwrite record ', &
                                                stream % nRecords, ' in file '''//trim(stream % filename)//''','
                   STREAM_ERROR_WRITE(trim(err_string))
                   write(err_string,'(a)') '    but clobber_mode is set to ''append''.'
                   STREAM_ERROR_WRITE(trim(err_string))
                   ierr = MPAS_STREAM_MGR_ERR_CLOBBER_REC
               else
                   ierr = MPAS_STREAM_MGR_ERROR
               end if

               if ( swapRecords ) then
                   stream % nRecords = tempRecord
               end if
               return
           end if

           if ( swapRecords ) then
               stream % nRecords = tempRecord
           end if
        end if

    end subroutine write_stream !}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_read
    !
    !> \brief Read streams that are managed by an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  With no optional arguments, reads all streams whose alarms are ringing.
    !>  The "streamID" argument optionally specifies a regular expression of streams
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
    subroutine MPAS_stream_mgr_read(manager, streamID, timeLevel, mgLevel, rightNow, when, whence, actualWhen, ierr) !{{{

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
        integer :: threadNum
        logical :: readStreams

        threadNum = mpas_threading_get_thread_num()

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


        if ( threadNum == 0 ) then
           !
           ! If a stream is specified, we process just that stream; otherwise,
           !    process all streams
           !
           if (present(streamID)) then
               readStreams = .false.
               nullify(stream_cursor)
               do while (MPAS_stream_list_query(manager % streams, streamID, stream_cursor, ierr=ierr))
                   STREAM_DEBUG_WRITE('-- Handling read of stream '//trim(stream_cursor % name))
   
                   ! Verify that the stream is an input stream
                   if (stream_cursor % direction == MPAS_STREAM_INPUT .or. stream_cursor % direction == MPAS_STREAM_INPUT_OUTPUT) then
                      readStreams = .true.
                      call read_stream(manager, stream_cursor, local_timeLevel, local_mgLevel, local_rightNow, local_when, &
                                       local_whence, actualWhen, local_ierr)
                   end if
               end do

               if ( .not. readStreams ) then
                   STREAM_ERROR_WRITE('No input stream matching '//trim(streamID)//' exists in call to MPAS_stream_mgr_read().')
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
                       stream_cursor % direction == MPAS_STREAM_INPUT_OUTPUT) then
   
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
        end if

        if (present(ierr)) ierr = local_ierr

    end subroutine MPAS_stream_mgr_read !}}}


    !-----------------------------------------------------------------------
    !  routine read_stream
    !
    !> \brief Handle the reading of a stream pointed to by the stream list node
    !> \author Michael Duda, Doug Jacobsen
    !> \date   4 September 2014
    !> \details
    !>  Private subroutine to handle the details of actually reading a stream.
    !
    !-----------------------------------------------------------------------
    subroutine read_stream(manager, stream, timeLevel, mgLevel, forceReadNow, when, whence, actualWhen, ierr) !{{{

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

        integer :: blockID_local
        character (len=StrKIND) :: err_string

        type (MPAS_stream_list_type), pointer :: alarm_cursor
        type (MPAS_Time_type) :: now_time, ref_time, temp_time
        type (MPAS_TimeInterval_type) :: temp_interval
        type (MPAS_TimeInterval_type) :: filename_interval
        character (len=StrKIND) :: temp_filename
        character (len=StrKIND) :: temp_actualWhen
        logical :: ringing_alarm
        integer :: temp_maxRecords
        integer :: local_ierr

        type (MPAS_Time_Type) :: currentTime, filenameTime
        type (MPAS_TimeInterval_Type) :: filenameInterval
        type (MPAS_Time_Type) :: whenTime, firstTime, secondTime
        type (MPAS_TimeInterval_Type) :: firstDiff, secondDiff

        type (MPAS_Stream_type) :: testStream
        character (len=StrKIND) :: test_when
        character (len=StrKIND) :: test_filename
        character (len=StrKIND) :: test_actualWhen
        integer :: test_record, test_maxRecords
        logical :: retestFile, rebuildStream
        integer :: threadNum

        threadNum = mpas_threading_get_thread_num()

        ! For now, reading only supports global reads.
        ! In the future, if reading is supposed to support local reads as well,
        ! blockID can become an input argument and then blockID_local can be
        ! set to blockID
        blockID_local = -1

        ierr = MPAS_STREAM_MGR_NOERR

        if ( threadNum == 0 ) then
           rebuildStream = .false.

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
           ! First we need to build the filename for the current read time.
           !
           if ( stream % filename_interval /= 'none' ) then
               call mpas_set_time(now_time, dateTimeString=when, ierr=local_ierr)
               call mpas_set_timeInterval(filename_interval, timeString=stream % filename_interval)
               
               call mpas_build_stream_filename(stream % referenceTime, now_time, filename_interval, stream % filename_template, blockID_local, temp_filename, ierr=local_ierr)
           else
               call mpas_expand_string(when, blockID_local, stream % filename_template, temp_filename)
           end if

           STREAM_DEBUG_WRITE(' -- Stream filename is: ' COMMA trim(temp_filename) )

           !
           ! If the stream is not valid, assume that we have not yet written this
           ! stream, in which case we create the stream from scratch
           !
           if (.not. stream % valid) then
               stream % filename = temp_filename

               !
               ! Build stream from pools of fields and attributes
               !
               allocate(stream % stream)
               call MPAS_createStream(stream % stream, manager % ioContext, stream % filename, stream % io_type, MPAS_IO_READ, &
                                      precision=stream % precision, ierr=local_ierr)
               if (local_ierr /= MPAS_STREAM_NOERR) then
                   if (local_ierr == MPAS_IO_ERR_NOEXIST_READ) then
                       write(err_string,'(a)') 'Stream '''//trim(stream % name)//''' attempted to read non-existent file '''//trim(stream % filename)//''''
                       STREAM_ERROR_WRITE(trim(err_string))
                       ierr = MPAS_STREAM_MGR_ERROR
                   else
                       ierr = MPAS_STREAM_MGR_ERROR
                   end if
                   return
               end if

               call build_stream(stream, MPAS_STREAM_INPUT, manager % allFields, manager % allPackages, timeLevel, mgLevel, local_ierr)
               if (local_ierr /= MPAS_STREAM_NOERR) then
                   ierr = MPAS_STREAM_MGR_ERROR
                   return
               end if
               stream % timeLevel = timeLevel

               stream % valid = .true.
           else if (temp_filename /= stream % filename) then
              STREAM_DEBUG_WRITE('-- Changing filename from '//trim(stream % filename)//' to '//trim(temp_filename))

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
              call MPAS_createStream(stream % stream, manager % ioContext, stream % filename, stream % io_type, MPAS_IO_READ, precision=stream % precision, ierr=local_ierr)
              if (local_ierr /= MPAS_STREAM_NOERR) then
                  if (local_ierr == MPAS_IO_ERR_NOEXIST_READ) then
                      write(err_string,'(a)') 'Stream '''//trim(stream % name)//''' attempted to read non-existent file '''//trim(stream % filename)//''''
                      STREAM_ERROR_WRITE(trim(err_string))
                      ierr = MPAS_STREAM_MGR_ERROR
                  else
                      ierr = MPAS_STREAM_MGR_ERROR
                  end if
                  return
              end if

              call build_stream(stream, MPAS_STREAM_INPUT, manager % allFields, manager % allPackages, timeLevel, mgLevel, local_ierr)
              if (local_ierr /= MPAS_STREAM_NOERR) then
                  ierr = MPAS_STREAM_MGR_ERROR
                  return
              end if
              stream % timeLevel = timeLevel
           end if

           STREAM_DEBUG_WRITE(' Seeking time of ' COMMA trim(when))

           !
           ! With multiple times per file, we need to get the record number from MPAS_seekStream.
           !
           stream % nRecords = MPAS_seekStream(stream % stream, when, whence, temp_actualWhen, maxRecords=temp_maxRecords, ierr=local_ierr)

           if ( stream % nRecords == 0 .and. temp_maxRecords == 0 ) then
               stream % nRecords = 1
               STREAM_WARNING_WRITE('File ' COMMA trim(stream % filename) COMMA ' does not contain a seekable xtime variable. Forcing a read of the first time record.')
           else if (stream % nRecords /= 0) then
               STREAM_DEBUG_WRITE('   Seeked record is: ' COMMA stream % nRecords COMMA ' out of ' COMMA temp_maxRecords COMMA ' with a time stamp of ' COMMA trim(temp_actualWhen) COMMA ' filename was ' COMMA trim(stream % filename))
           else if (temp_maxRecords /= 0 .and. whence == MPAS_STREAM_EXACT_TIME) then
               STREAM_ERROR_WRITE('File ' COMMA trim(stream % filename) COMMA ' does not contain the time ' COMMA trim(when))
               ierr = MPAS_STREAM_MGR_ERROR
               return
           end if

           retestFile = .false.
           if ( trim(stream % filename_interval) /= 'none' .and. whence /= MPAS_STREAM_EXACT_TIME ) then
              currentTime = mpas_get_clock_time(manager % streamClock, MPAS_NOW, ierr=local_ierr)
              call mpas_set_timeInterval(filenameInterval, timeString=stream % filename_interval, ierr=local_ierr)

              ! Need to handle the case where the time requested was not found.
              !
              ! Things that need to be handled here are when we're at the beginning
              ! or end of a file and we're looking for the next or previous time
              ! record.
              !
              ! Currently this only checks one file each direction (forward or
              ! backward). It will fail finding a file more than one interval away
              ! from when.
              if ( stream % nRecords == 0) then
                 if ( ( whence == MPAS_STREAM_LATEST_BEFORE .and. temp_actualWhen /= when ) .or. whence == MPAS_STREAM_LATEST_STRICTLY_BEFORE ) then
                    ! Subtract filename_interval from when, build new filename, and
                    ! check for a time latest before in that file.
                    filenameTime = currentTime - filenameInterval
                    retestFile = .true.
                    STREAM_DEBUG_WRITE(' Retest latest before...')
                 else if ( ( whence == MPAS_STREAM_EARLIEST_AFTER .and. temp_actualWhen /= when ).or. whence == MPAS_STREAM_EARLIEST_STRICTLY_AFTER ) then
                    ! Add filename_interval from when, build new filename, and
                    ! check for a time latest before in that file.
                    filenameTime = currentTime + filenameInterval
                    retestFile = .true.
                    STREAM_DEBUG_WRITE(' Retest earliest after...')
                 end if
              else
                 ! If time was found, and we were looking for nearest need to make sure nearest isn't in previous or next file.
                 !
                 ! This only needs to be checked if we found the first or last time slice in the file.
                 if ( whence == MPAS_STREAM_NEAREST ) then
                    if ( stream % nRecords == 1 .and. stream % nRecords == temp_maxRecords ) then
                        call mpas_set_time(temp_time, dateTimeString=temp_actualWhen)

                        ! If an exact time was found, read that one, and don't bother re-testing.
                        if ( currentTime == temp_time ) then
                           retestFile = .false.

                        ! If current time is before the time that was read, re-test using the previous file
                        else if ( currentTime < temp_time ) then
                           filenameTime = currentTime - filenameInterval
                           retestFile = .true.
                           STREAM_DEBUG_WRITE('Retest nearest prev file')

                        ! If current time is before the time that was read, re-test using the next file
                        else if ( currentTime > temp_time ) then
                           filenameTime = currentTime + filenameInterval
                           retestFile = .true.
                           STREAM_DEBUG_WRITE('Retest nearest next file')
                        end if
                    else if ( stream % nRecords == 1 ) then
                        ! Subtract filename_interval from when, build new filename, and check for nearest time in that file.
                        ! Compare the two, and keep the one closest to when.
                        filenameTime = currentTime - filenameInterval
                        retestFile = .true.
                        STREAM_DEBUG_WRITE('Retest nearest beginning')
                    else if ( stream % nRecords == temp_maxRecords ) then
                        ! Add filename_interval from when, build new filename, and check for nearest time in that file.
                        ! Compare the two, and keep the one closest to when.
                        filenameTime = currentTime + filenameInterval
                        retestFile = .true.
                        STREAM_DEBUG_WRITE('Retest nearest end')
                    end if
                 end if
              end if
           end if

           if ( retestFile ) then
              STREAM_DEBUG_WRITE(' --- Retesting file... ')
              call mpas_get_time(filenameTime, dateTimeString=test_when)

              call mpas_set_timeInterval(filename_interval, timeString=stream % filename_interval)
              
              call mpas_build_stream_filename(stream % referenceTime, filenameTime, filename_interval, stream % filename_template, blockID_local, test_filename, ierr=local_ierr)

              STREAM_DEBUG_WRITE(' --- Retesting filename is ' COMMA trim(test_filename))

              inquire(file=trim(test_filename), exist=retestFile)

              ! If file exists, the testing stream needs to be built.
              if ( retestFile ) then
                  call mpas_createStream(testStream, manager % ioContext, test_filename, stream % io_type, MPAS_IO_READ, precision=stream % precision, ierr=local_ierr)
              else
                  STREAM_DEBUG_WRITE(' Filename: ' COMMA trim(test_filename) COMMA ' does not exist.')
              end if
           end if

           ! Only continue testing file it if was found.
           if ( retestFile ) then
              test_record = MPAS_seekStream(testStream, when, whence, test_actualWhen, maxRecords=test_maxRecords, ierr=local_ierr)

              STREAM_DEBUG_WRITE(' -- Test record is ' COMMA test_record COMMA ' out of ' COMMA test_maxRecords COMMA ' with a time of ' COMMA trim(test_actualWhen))

              if ( test_record /= 0 ) then
                 if ( whence == MPAS_STREAM_NEAREST ) then
                     call mpas_set_time(whenTime, dateTimeString=when)
                     call mpas_set_time(firstTime, dateTimeString=temp_actualWhen)
                     call mpas_set_time(secondTime, dateTimeString=test_actualWhen)

                     ! Build first diff
                     if ( firstTime > whenTime ) then
                        firstDiff = firstTime - whenTime
                     else
                        firstDiff = whenTime - firstTime
                     end if

                     ! Build second diff
                     if ( secondTime > whenTime ) then
                        secondDiff = secondTime - whenTime
                     else
                        secondDiff = whenTime - secondTime
                     end if

                     ! Compare first and second diff, keeping the closest one to when.
                     ! Only need to rebuild stream if the second* ones are closer.
                     if ( secondDiff == firstDiff ) then

                        ! If times are equidistance, take the later of the two.
                        if ( firstTime > secondTime ) then
                           rebuildStream = .false.
                        else
                           rebuildStream = .true.
                        end if
                     else if ( secondDiff < firstDiff ) then
                        rebuildStream = .true.
                        STREAM_DEBUG_WRITE(' --- New time is closer than old time')
                     else
                        STREAM_DEBUG_WRITE(' --- Old time is closer than test time')
                     end if
                 else if ( stream % nRecords == 0 ) then
                     rebuildStream = .true.
                 end if
              else
                  rebuildStream = .false.
              end if
              call MPAS_closeStream(testStream, ierr=local_ierr)
           end if

           ! Rebuild stream if we need to, because a different file has a closer time.
           if ( rebuildStream ) then
              STREAM_DEBUG_WRITE(' --- rebuilding stream...')
              stream % filename = test_filename

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
              call MPAS_createStream(stream % stream, manager % ioContext, stream % filename, stream % io_type, MPAS_IO_READ, precision=stream % precision, ierr=local_ierr)
              if (local_ierr /= MPAS_STREAM_NOERR) then
                  ierr = MPAS_STREAM_MGR_ERROR
                  return
              end if

              call build_stream(stream, MPAS_STREAM_INPUT, manager % allFields, manager % allPackages, timeLevel, mgLevel, local_ierr)
              if (local_ierr /= MPAS_STREAM_NOERR) then
                  ierr = MPAS_STREAM_MGR_ERROR
                  return
              end if
              stream % timeLevel = timeLevel

              ! Set record number based on test_record from the read we just did.
              stream % nRecords = test_record
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
        end if

    end subroutine read_stream !}}}


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
    !  routine mpas_get_stream_filename
    !
    !> \brief Determine the filename that contains a specific time in a stream
    !> \author Michael Duda, Doug Jacobsen
    !> \date   21 August 2014
    !> \details
    !>  Given a stream manger, a stream name, and an optional time stamp,
    !>  return the filename that would contain that time stamp, or the filename
    !>  that should contain the current time as defined by the clock associated
    !>  with the stream manager.
    !>
    !>  Return error codes:
    !>  0 no error
    !>
    !>  NOTE: This routine does not support streamID regular expressions
    !-----------------------------------------------------------------------
    subroutine mpas_get_stream_filename(manager, streamID, when, blockID, filename, ierr) !{{{

        implicit none

        type (MPAS_streamManager_type), intent(in) :: manager !< Input: Stream manager to get stream from
        character (len=*), intent(in) :: streamID !< Input: Stream name to use for building the filename
        character (len=*), intent(in), optional :: when !< Optional Input: Time file should contain
        integer, intent(in), optional :: blockID !< Input: ID of the block that will be written
        character (len=StrKIND), intent(out) :: filename !< Output: Name of file containing time
        integer, intent(out), optional :: ierr !< Optional Output: Error code

        type (mpas_stream_list_type), pointer :: streamCursor

        integer :: err_local, blockID_local

        character(len=StrKIND) :: when_string
        type (MPAS_TimeInterval_type) :: filename_interval
        type (MPAS_Time_type) :: now_time

        ierr = 0

        if ( present(blockID) ) then
           blockID_local = blockID
        else
           blockID_local = -1
        end if

        nullify(streamCursor)
        if ( mpas_stream_list_query(manager % streams, streamID, streamCursor, err_local) ) then
           if ( present(when) ) then
              call mpas_set_time(now_time, dateTimeString=when, ierr=err_local)
           else
              now_time = mpas_get_clock_time(manager % streamClock, MPAS_NOW, ierr=err_local)
           end if

           !
           ! First we need to build the filename for the current read time.
           !
           if ( streamCursor % filename_interval /= 'none' ) then
               call mpas_set_timeInterval(filename_interval, timeString=streamCursor % filename_interval)
               call mpas_build_stream_filename(streamCursor % referenceTime, now_time, filename_interval, streamCursor % filename_template, blockID_local, filename, ierr=err_local)
           else
               call mpas_get_time(now_time, dateTimeString=when_string, ierr=err_local)
               call mpas_expand_string(when_string, blockID_local, streamCursor % filename_template, filename)
           end if

        else
            if ( present(ierr) ) ierr = MPAS_STREAM_MGR_ERROR
            if ( manager % errorLevel /= MPAS_STREAM_ERR_SILENT ) then
                write(stderrUnit, *) 'ERROR: Stream ' // trim(streamID) // ' does not exist. Filename creation failed...'
            end if
        end if

    end subroutine mpas_get_stream_filename !}}}


    !-----------------------------------------------------------------------
    !  routine mpas_build_stream_filename
    !
    !> \brief Construct the filename that contains a specific time in a stream
    !> \author Michael Duda, Doug Jacobsen
    !> \date   21 August 2014
    !> \details 
    !>  Given a filename template and the information necessary to determine the time
    !>  in the stream that matches a time available in any of the files associated with
    !>  the stream, returns a specific filename that should contain that time.
    !>
    !>  Filenames are assumed to start at the earliest time. In other words, a
    !>  file name is expanded using the earliest time that could possibly be
    !>  stored in the file.
    !>
    !>  This is a low level subroutine to complement the
    !>  mpas_get_stream_Filename routine
    !>  
    !>  Return error codes:
    !>  0 no error
    !-----------------------------------------------------------------------
    subroutine mpas_build_stream_filename(ref_time, when, filename_interval, filename_template, blockID, filename, ierr) !{{{

        implicit none

        type (MPAS_Time_type), intent(in) :: ref_time
        type (MPAS_Time_type), intent(in) :: when
        type (MPAS_TimeInterval_type), intent(in) :: filename_interval
        character(len=*), intent(in) :: filename_template
        integer, intent(in) :: blockID
        character(len=*), intent(out) :: filename
        integer, intent(out) :: ierr

        character(len=StrKIND) :: temp_string
        character(len=StrKIND) :: when_string
        type (MPAS_Time_type) :: filetime
        type (MPAS_TimeInterval_type) :: intv, rem, zeroIntv
        integer :: nrecs, nfiles, irec, direction
        logical :: in_future

        STREAM_DEBUG_WRITE(' ** Building Filename')

        ierr = 0

        ! If current time (when) is further ahead than ref_time
        ! the interval we want is when - ref_time (time between now and the reference time)
        if ( when >= ref_time ) then
           intv = when - ref_time
           direction = 1
        else
           intv = ref_time - when
           direction = -1
        end if

        call mpas_get_time(when, dateTimeString=temp_string)
        STREAM_DEBUG_WRITE(' ** when is: ' COMMA trim(temp_string))

        call mpas_get_time(ref_time, dateTimeString=temp_string)
        STREAM_DEBUG_WRITE(' ** ref_time is: ' COMMA trim(temp_string))

        call mpas_get_timeInterval(intv, timeString=temp_string)
        STREAM_DEBUG_WRITE(' ** intv is: ' COMMA trim(temp_string))

        call mpas_interval_division(ref_time, intv, filename_interval, nrecs, rem) 

!       STREAM_DEBUG_WRITE(' ** Divisions are: ' COMMA nrecs)

        call mpas_set_timeInterval(zeroIntv, s=0)

        if ( rem /= zeroIntv ) then
            ! direction == 1 means when is later than ref_time
            if (direction == 1) then
                filetime = when - rem
            else
                filetime = when + rem
                filetime = filetime - filename_interval
            end if
        else
            filetime = when
        end if
        call mpas_get_time(filetime, dateTimeString=when_string)
        STREAM_DEBUG_WRITE(' ** filetime start is: ' COMMA trim(when_string))

        call mpas_expand_string(when_string, blockID, filename_template, filename)

    end subroutine mpas_build_stream_filename !}}}


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
    subroutine build_stream(stream, direction, allFields, allPackages, timeLevelIn, mgLevelIn, ierr) !{{{

        implicit none

        type (MPAS_stream_list_type), intent(inout) :: stream
        integer, intent(in) :: direction
        type (MPAS_Pool_type), intent(in) :: allFields
        type (MPAS_Pool_type), intent(in) :: allPackages
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

        character (len=StrKIND), pointer :: packages
        logical :: active_field
        integer :: err_level


        if (direction == MPAS_STREAM_OUTPUT) then

            !
            ! Write attributes to stream
            !
            call mpas_pool_begin_iteration(stream % att_pool)
            do while (mpas_pool_get_next_member(stream % att_pool, itr))
                if ( itr % memberType == MPAS_POOL_CONFIG) then
                    if ( itr % dataType == MPAS_POOL_REAL ) then
                        call mpas_pool_get_config(stream % att_pool, itr % memberName, realAtt)
                        call mpas_writeStreamAtt(stream % stream, itr % memberName, realAtt, syncVal=.false., ierr=local_ierr)
                    else if ( itr % dataType == MPAS_POOL_INTEGER ) then
                        call mpas_pool_get_config(stream % att_pool, itr % memberName, intAtt)
                        call mpas_writeStreamAtt(stream % stream, itr % memberName, intAtt, syncVal=.false., ierr=local_ierr)
                    else if ( itr % dataType == MPAS_POOL_CHARACTER ) then
                        call mpas_pool_get_config(stream % att_pool, itr % memberName, charAtt)
                        call mpas_writeStreamAtt(stream % stream, itr % memberName, charAtt, syncVal=.false., ierr=local_ierr)
                    else if ( itr % dataType == MPAS_POOL_LOGICAL ) then
                        call mpas_pool_get_config(stream % att_pool, itr % memberName, logAtt)
                        if (logAtt) then
                            call mpas_writeStreamAtt(stream % stream, itr % memberName, 'YES', syncVal=.false., ierr=local_ierr)
                        else
                            call mpas_writeStreamAtt(stream % stream, itr % memberName, 'NO', syncVal=.false., ierr=local_ierr)
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
            call mpas_writeStreamAtt(stream % stream, 'file_id', file_id, syncVal=.true., ierr=local_ierr)
            if (local_ierr /= MPAS_STREAM_NOERR) then
                ierr = MPAS_STREAM_MGR_ERROR
                return
            end if

        end if


        ierr = MPAS_STREAM_MGR_NOERR

        call mpas_pool_begin_iteration(stream % field_pool)

        FIELD_LOOP: do while ( mpas_pool_get_next_member(stream % field_pool, itr) )

            if (itr % memberType == MPAS_POOL_CONFIG) then

                err_level = mpas_pool_get_error_level()
                call mpas_pool_set_error_level(MPAS_POOL_SILENT)

                nullify(packages)
                call mpas_pool_get_config(stream % field_pkg_pool, trim(itr % memberName)//':packages', packages)
                if (associated(packages)) then
                    active_field = parse_package_list(allPackages, trim(packages))
                else
                    active_field = .true.
                end if
                call mpas_pool_set_error_level(err_level)

                STREAM_DEBUG_WRITE('Is field '''//trim(itr % memberName)//''' active in stream '''//trim(stream % name)//'? ' COMMA active_field)

                if (.not. active_field) cycle FIELD_LOOP

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
                                 call mpas_pool_get_field(allFields, itr % memberName, char1d, timeLevel)
                                 call MPAS_streamAddField(stream % stream, char1d)
                        end select
                end select

            end if
        end do FIELD_LOOP

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
    logical function stream_active_pkg_check(stream) !{{{

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

    end function stream_active_pkg_check !}}}


    !-----------------------------------------------------------------------
    !  routine parse_package_list
    !
    !> \brief Parses a semi-colon-separated list of package names, indicating whether any are active
    !> \author Michael Duda
    !> \date   19 March 2015
    !> \details 
    !>  This function determines whether any of the named strings in 
    !>  the semi-colon-separated list provided in the 'packages' argument are
    !>  active. 
    !>  If any of the packages does not exist in the package pool, the optional
    !>  argument ierr is set to a non-zero value; otherwise, if all packages exist, 
    !>  ierr will be set to zero upon return.
    !
    !-----------------------------------------------------------------------
    logical function parse_package_list(package_pool, packages, ierr) result(active)
 
        implicit none
 
        type (mpas_pool_type), intent(in) :: package_pool
        character (len=*), intent(in) :: packages
        integer, intent(out), optional :: ierr
 
        integer :: i, j, slen
        integer :: err_level
        logical, pointer :: pkg_val
 
 
        if (present(ierr)) ierr = 0
 
        slen = len_trim(packages)
 

        !
        ! No packages
        !
        if (slen == 0) then
            active = .true.
            return
        end if
  
        active = .false.

        err_level = mpas_pool_get_error_level()
        call mpas_pool_set_error_level(MPAS_POOL_SILENT)


        !
        ! Possible semi-colons in 'packages'
        !
        i = 1
        j = index(packages,';')
        do while (j >= i)
            if (j > i) then
                nullify(pkg_val)
                call mpas_pool_get_package(package_pool, packages(i:j-1)//'Active', pkg_val)
                if (associated(pkg_val)) then
                    if (pkg_val) then
                        active = .true.
                        call mpas_pool_set_error_level(err_level)
                        return
                    end if
                else
                    if (present(ierr)) ierr = 1
                end if
            end if
            i = j+1
            j = index(packages(i:slen),';') + i - 1
        end do


        !
        ! No more semi-colons to worry about
        !
        if (i < slen) then
            nullify(pkg_val)
            call mpas_pool_get_package(package_pool, packages(i:slen)//'Active', pkg_val)
            if (associated(pkg_val)) then
                if (pkg_val) then
                    active = .true.
                    call mpas_pool_set_error_level(err_level)
                    return
                end if
            else
                if (present(ierr)) ierr = 1
            end if
        end if

        call mpas_pool_set_error_level(err_level)

    end function parse_package_list


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
                            call mpas_dmpar_field_halo_exch(real1DField % block % domain, fieldItr % memberName, timeLevel=timeLevel)
                        end if
                    else if ( fieldInfo % fieldType == MPAS_POOL_INTEGER ) then
                        if ( timeLevel <= fieldInfo % nTimeLevels ) then
                          call mpas_pool_get_field(allFields, fieldItr % memberName, int1DField, timeLevel)
                        else
                          call mpas_pool_get_field(allFields, fieldItr % memberName, int1DField, 1)
                        end if
                        if ( is_decomposed_dim(int1DField % dimNames(1))) then
                            STREAM_DEBUG_WRITE(' -- Exchange halo for '//trim(fieldItr % memberName))
                            call mpas_dmpar_field_halo_exch(int1DField % block % domain, fieldItr % memberName, timeLevel=timeLevel)
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
                            call mpas_dmpar_field_halo_exch(real2DField % block % domain, fieldItr % memberName, timeLevel=timeLevel)
                        end if
                    else if ( fieldInfo % fieldType == MPAS_POOL_INTEGER ) then
                        if ( timeLevel <= fieldInfo % nTimeLevels ) then
                            call mpas_pool_get_field(allFields, fieldItr % memberName, int2DField, timeLevel)
                        else
                            call mpas_pool_get_field(allFields, fieldItr % memberName, int2DField, 1)
                        end if
                        if ( is_decomposed_dim(int2DField % dimNames(2))) then
                            STREAM_DEBUG_WRITE(' -- Exchange halo for '//trim(fieldItr % memberName))
                            call mpas_dmpar_field_halo_exch(int2DField % block % domain, fieldItr % memberName, timeLevel=timeLevel)
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
                            call mpas_dmpar_field_halo_exch(real3DField % block % domain, fieldItr % memberName, timeLevel=timeLevel)
                        end if
                    else if ( fieldInfo % fieldType == MPAS_POOL_INTEGER ) then
                        if ( timeLevel <= fieldInfo % nTimeLevels ) then
                            call mpas_pool_get_field(allFields, fieldItr % memberName, int3DField, timeLevel)
                        else
                            call mpas_pool_get_field(allFields, fieldItr % memberName, int3DField, 1)
                        end if
                        if ( is_decomposed_dim(int3DField % dimNames(3))) then
                            STREAM_DEBUG_WRITE(' -- Exchange halo for '//trim(fieldItr % memberName))
                            call mpas_dmpar_field_halo_exch(int3DField % block % domain, fieldItr % memberName, timeLevel=timeLevel)
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
                            call mpas_dmpar_field_halo_exch(real4DField % block % domain, fieldItr % memberName, timeLevel=timeLevel)
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
                            call mpas_dmpar_field_halo_exch(real5DField % block % domain, fieldItr % memberName, timeLevel=timeLevel)
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
    logical function is_decomposed_dim(dimName) !{{{

        implicit none

        character(len=*), intent(in) :: dimName

        if (trim(dimName) == 'nCells' .or. &
            trim(dimName) == 'nEdges' .or. &
            trim(dimName) == 'nVertices') then

            is_decomposed_dim = .true.

        else

            is_decomposed_dim = .false.

        end if

    end function is_decomposed_dim !}}}


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
    subroutine prewrite_reindex(allFields, streamFields) !{{{

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

        integer :: i, j, threadNum

        threadNum = mpas_threading_get_thread_num()

        if ( threadNum == 0 ) then
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
        end if

    end subroutine prewrite_reindex !}}}


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
    subroutine postwrite_reindex(allFields, streamFields) !{{{

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

        integer :: i, j, threadNum

        threadNum = mpas_threading_get_thread_num()

        if ( threadNum == 0 ) then
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
        end if

    end subroutine postwrite_reindex !}}}


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
    subroutine postread_reindex(allFields, streamFields) !{{{

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

    end subroutine postread_reindex !}}}


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_begin_iteration
    !
    !> \brief Reset iterator within stream manager
    !> \author Doug Jacobsen, Michael Duda
    !> \date   03/03/2015
    !> \details
    !>  If the optional 'streamID' argument is provided, this routine resets 
    !>  the iterator  within a stream manager so that streams may subsequently 
    !>  be iterated over using the MPAS_stream_mgr_get_next_stream function. 
    !>
    !>  If an optional stream name is provided via the 'streamID' argument, this
    !>  routine will reset the iterator for fields within the specified stream, 
    !>  which may subsequently iterated over using the 
    !>  MPAS_stream_mgr_get_next_field() routine.
    !>
    !>  NOTE: This routine does not support regular expressions for StreamID
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_begin_iteration(manager, streamID, ierr) !{{{

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager     !< Input: Stream manager to begin iteration for
        character (len=*), intent(in), optional :: streamID          !< Input: Name of the stream to iterate over
        integer, intent(out), optional :: ierr                       !< Output: Return error code

        type (MPAS_stream_list_type), pointer :: stream
        integer :: err_local, threadNum

        threadNum = mpas_threading_get_thread_num()


        if (.not. present(streamID)) then

           if ( threadNum == 0 ) then
              nullify(manager % currentStream)
           end if

        else

           !
           ! Check that stream exists
           !
           nullify(stream)
           if (.not. MPAS_stream_list_query(manager % streams, streamID, stream, ierr=err_local)) then
               STREAM_ERROR_WRITE('Requested stream '//trim(streamID)//' does not exist in stream manager')
               if (present(ierr)) ierr = MPAS_STREAM_MGR_ERROR
               return
           end if

           call mpas_pool_begin_iteration(stream % field_pool)

        end if

        !$omp barrier

    end subroutine MPAS_stream_mgr_begin_iteration !}}}


    !-----------------------------------------------------------------------
    !  logical function MPAS_stream_mgr_get_next_stream
    !
    !> \brief Retrieve information about next stream
    !> \author Doug Jacobsen
    !> \date   03/03/2015
    !> \details
    !>  This routine advances the internal iterator to the next stream, and
    !>  returns information about it.
    !
    !-----------------------------------------------------------------------
    logical function MPAS_stream_mgr_get_next_stream(manager, streamID, directionProperty, activeProperty, & !{{{
                                                     immutableProperty, filenameTemplateProperty, &
                                                     referenceTimeProperty, recordIntervalProperty, precisionProperty, &
                                                     filenameIntervalProperty, clobberProperty) result(validStream)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager                   !< Input: Stream manager to iterate over
        character (len=StrKIND), intent(out), optional :: streamID                 !< Output: Name of next stream
        integer, intent(out), optional :: directionProperty                        !< Output: Integer describing the direction of the stream
        logical, intent(out), optional :: activeProperty                           !< Output: Logical describing if the stream is active or not
        logical, intent(out), optional :: immutableProperty                        !< Output: Logical describing if the stream is immutable or not
        character (len=StrKIND), intent(out), optional :: filenameTemplateProperty !< Output: String containing the filename template for the stream
        character (len=StrKIND), intent(out), optional :: referenceTimeProperty    !< Output: String containing the reference time for the stream
        character (len=StrKIND), intent(out), optional :: recordIntervalProperty   !< Output: String containing the record interval for the stream
        integer, intent(out), optional :: precisionProperty                        !< Output: Integer describing the precision of the stream
        character (len=StrKIND), intent(out), optional :: filenameIntervalProperty !< Output: String containing the filename interval for the stream
        integer, intent(out), optional :: clobberProperty                          !< Output: Interger describing the clobber mode of the stream

        integer :: threadNum

        threadNum = mpas_threading_get_thread_num()

        if ( associated(manager % currentStream) ) then
            if (.not. associated(manager % currentStream % next) ) then
                validStream = .false.
                return
            else
                validStream = .true.
                if ( threadNum == 0 ) then
                   manager % currentStream => manager % currentStream % next
                end if
            end if
        else if ( associated(manager % streams % head) ) then
           validStream = .true.
           if ( threadNum == 0 ) then
              manager % currentStream => manager % streams % head
           end if
        else
           validStream = .false.
           return
        end if

        !$omp barrier

        if ( present(streamID) ) then
           streamID = manager % currentStream % name
        end if

        if ( present(directionProperty) ) then
            directionProperty = manager % currentStream % direction
        end if

        if ( present(activeProperty) ) then
            activeProperty = manager % currentStream % active_stream
        end if

        if ( present(immutableProperty) ) then
            immutableProperty = manager % currentStream % immutable
        end if

        if ( present(filenameTemplateProperty) ) then
            filenameTemplateProperty = manager % currentStream % filename_template
        end if

        if ( present(referenceTimeProperty) ) then
            call mpas_get_time(manager % currentStream % referenceTime, dateTimeString=referenceTimeProperty)
        end if

        if ( present(recordIntervalProperty) ) then
            call mpas_get_timeInterval(manager % currentStream % recordInterval, timeString=recordIntervalProperty)
        end if

        if ( present(precisionProperty) ) then
            precisionProperty = manager % currentStream % precision
        end if

        if ( present(filenameIntervalProperty)) then
            filenameIntervalProperty = manager % currentStream % filename_interval
        end if

        if ( present(clobberProperty) ) then
            clobberProperty = manager % currentStream % clobber_mode
        end if

    end function MPAS_stream_mgr_get_next_stream !}}}


    !-----------------------------------------------------------------------
    !  logical function MPAS_stream_mgr_get_next_field
    !
    !> \brief Retrieve the name of the next field in a stream
    !> \author Michael Duda
    !> \date   06 March 2015
    !> \details
    !>  This function queries the name of the next field belonging to the stream
    !>  indicated by the streamID argument. Before the first call to this
    !>  routine for a stream, the routine MPAS_stream_mgr_begin_iteration() must
    !>  have been called to initialize iteration for the specified streamID in
    !>  the stream manager.
    !>
    !>  This function returns .TRUE. if the stream contains another field,
    !>  whether active or not, in which case the output argument fieldName 
    !>  provides the name of this field, and .FALSE. otherwise. If a field name
    !>  is returned, the optional logical argument isActive may be used to
    !>  determine whether the field is currently active in the stream.
    !>  NOTE: This routine does not support regular expressions for StreamID
    !
    !-----------------------------------------------------------------------
    logical function MPAS_stream_mgr_get_next_field(manager, streamID, fieldName, isActive) result(validField) !{{{

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager    !< Input: Stream manager containing stream indicated by streamID
        character (len=*), intent(in) :: streamID                   !< Input: Name of the stream to iterate over
        character (len=StrKIND), intent(out) :: fieldName           !< Output: Name of the next field in streamID
        logical, intent(out), optional :: isActive                  !< Output: Name of the next field in streamID

        type (MPAS_stream_list_type), pointer :: stream
        type (mpas_pool_iterator_type) :: poolItr
        character (len=StrKIND), pointer :: packages
        integer :: err_level
        integer :: err_local


        validField = .false.

        !
        ! Check that stream exists
        !
        nullify(stream)
        if (.not. MPAS_stream_list_query(manager % streams, streamID, stream, ierr=err_local)) then
            STREAM_ERROR_WRITE('Requested stream '//trim(streamID)//' does not exist in stream manager')
            return
        end if

        !
        ! Get the name of the next field, if one exists
        !
        if (mpas_pool_get_next_member(stream % field_pool, poolItr)) then
            fieldName = poolItr % memberName
            validField = .true.
        end if

        !
        ! Optionally, set the isActive variable
        !
        if (validField .and. present(isActive)) then
            err_level = mpas_pool_get_error_level()
            call mpas_pool_set_error_level(MPAS_POOL_SILENT)

            nullify(packages)
            call mpas_pool_get_config(stream % field_pkg_pool, trim(fieldName)//':packages', packages)
            if (associated(packages)) then
                isActive = parse_package_list(manager % allPackages, trim(packages))
            else
                isActive = .true.
            end if
    
            call mpas_pool_set_error_level(err_level)
        end if

    end function MPAS_stream_mgr_get_next_field !}}}


    !-----------------------------------------------------------------------
    !  logical function MPAS_stream_mgr_stream_exists
    !
    !> \brief Determine if a stream exists in a stream manager
    !> \author Doug Jacobsen
    !> \date   07/29/2015
    !> \details
    !>  This function takes a stream manager and the name of a stream, and
    !>  returns a logical describing if the stream exists within the manager or
    !>  not. It gives no information about if the stream will be "handled" at
    !>  any point in time, only if it exists within the manager queried.
    !>  NOTE: This routine does not support regular expressions for StreamID
    !
    !-----------------------------------------------------------------------
    logical function MPAS_stream_mgr_stream_exists(manager, streamID) result(validStream)!{{{

       implicit none

       type (MPAS_streamManager_type), intent(inout) :: manager
       character (len=*), intent(in) :: streamID

       type (mpas_stream_list_type), pointer :: streamPtr

       nullify(streamPtr)
       validStream = mpas_stream_list_query(manager % streams, streamID, streamPtr)

       return

    end function MPAS_stream_mgr_stream_exists!}}}
 

end module mpas_stream_manager


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! C interface routines for building streams at run-time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine stream_mgr_create_stream_c(manager_c, streamID_c, direction_c, filename_c, filename_intv_c, ref_time_c, rec_intv_c, &
                                      immutable_c, precision_c, clobber_c, iotype_c, ierr_c) bind(c) !{{{

    use mpas_c_interfacing, only : mpas_c_to_f_string
    use iso_c_binding, only : c_char, c_int, c_ptr, c_f_pointer
    use mpas_derived_types, only : MPAS_streamManager_type, &
                                        MPAS_STREAM_MGR_NOERR, MPAS_STREAM_PROPERTY_FILENAME, &
                                        MPAS_STREAM_PROPERTY_FILENAME_INTV, MPAS_STREAM_PROPERTY_REF_TIME, &
                                        MPAS_STREAM_PROPERTY_RECORD_INTV, MPAS_STREAM_PROPERTY_PRECISION, &
                                        MPAS_STREAM_PROPERTY_CLOBBER, MPAS_STREAM_CLOBBER_NEVER, MPAS_STREAM_CLOBBER_APPEND, &
                                        MPAS_STREAM_CLOBBER_TRUNCATE, MPAS_STREAM_CLOBBER_OVERWRITE, MPAS_STREAM_PROPERTY_IOTYPE
    use mpas_stream_manager, only : MPAS_stream_mgr_create_stream, MPAS_stream_mgr_set_property
    use mpas_kind_types, only : StrKIND
    use mpas_io, only : MPAS_IO_SINGLE_PRECISION, MPAS_IO_DOUBLE_PRECISION, MPAS_IO_NATIVE_PRECISION, &
                        MPAS_IO_NETCDF, MPAS_IO_PNETCDF, MPAS_IO_NETCDF4, MPAS_IO_PNETCDF5
    use mpas_io_units, only : stderrUnit

    implicit none

    type (c_ptr) :: manager_c
    character(kind=c_char) :: streamID_c(*)
    integer(kind=c_int) :: direction_c
    character(kind=c_char) :: filename_c(*)
    character(kind=c_char) :: filename_intv_c(*)
    character(kind=c_char) :: ref_time_c(*)
    character(kind=c_char) :: rec_intv_c(*)
    integer(kind=c_int) :: immutable_c
    integer(kind=c_int) :: precision_c
    integer(kind=c_int) :: clobber_c
    integer(kind=c_int) :: iotype_c
    integer(kind=c_int) :: ierr_c

    type (MPAS_streamManager_type), pointer :: manager
    character(len=StrKIND) :: streamID, filename, filename_interval, reference_time, record_interval
    integer :: direction, immutable, prec, ierr
    integer :: clobber_mode, iotype

    call c_f_pointer(manager_c, manager)
    call mpas_c_to_f_string(streamID_c, streamID)
    direction = direction_c
    call mpas_c_to_f_string(filename_c, filename)
    call mpas_c_to_f_string(filename_intv_c, filename_interval)
    call mpas_c_to_f_string(ref_time_c, reference_time)
    call mpas_c_to_f_string(rec_intv_c, record_interval)
    immutable = immutable_c

    if (precision_c == 4) then
        prec = MPAS_IO_SINGLE_PRECISION
    else if (precision_c == 8) then
        prec = MPAS_IO_DOUBLE_PRECISION
    else
        prec = MPAS_IO_NATIVE_PRECISION
    end if

    if (clobber_c == 0) then
       clobber_mode = MPAS_STREAM_CLOBBER_NEVER
    else if (clobber_c == 1) then
       clobber_mode = MPAS_STREAM_CLOBBER_APPEND
    else if (clobber_c == 2) then
       clobber_mode = MPAS_STREAM_CLOBBER_TRUNCATE
    else if (clobber_c == 3) then
       clobber_mode = MPAS_STREAM_CLOBBER_OVERWRITE
    else
       clobber_mode = MPAS_STREAM_CLOBBER_NEVER
    end if

    if (iotype_c == 0) then
       iotype = MPAS_IO_PNETCDF
    else if (iotype_c == 1) then
       iotype = MPAS_IO_PNETCDF5
    else if (iotype_c == 2) then
       iotype = MPAS_IO_NETCDF
    else if (iotype_c == 3) then
       iotype = MPAS_IO_NETCDF4
    else
       iotype = MPAS_IO_PNETCDF
    end if

    STREAM_DEBUG_WRITE('Creating stream from c...')
    
    !
    ! For immutable streams, the stream should have already been defined at this point, and
    !    all we need to do is update the stream's filename template;
    !    otherwise, we need to create a new stream
    !
    ierr = 0
    if (immutable == 1) then
        call MPAS_stream_mgr_set_property(manager, streamID, MPAS_STREAM_PROPERTY_FILENAME, filename, ierr=ierr)

        ! If we can't set a property on this immutable stream, most likely the stream wasn't defined in the core's Registry.xml file
        if (ierr /= MPAS_STREAM_MGR_NOERR) then
            write(stderrUnit,*) '********************************************************************************'
            write(stderrUnit,*) ' Error: Stream '''//trim(streamID)//''' was not defined in the Registry.xml file as '
            write(stderrUnit,*) '        an immutable stream. Immutable streams may only be defined in the Registry.xml'
            write(stderrUnit,*) '        file for a core.'
            write(stderrUnit,*) '********************************************************************************'
            ierr_c = 1
            return
        end if
        call MPAS_stream_mgr_set_property(manager, streamID, MPAS_STREAM_PROPERTY_PRECISION, prec, ierr=ierr)
        call MPAS_stream_mgr_set_property(manager, streamID, MPAS_STREAM_PROPERTY_CLOBBER, clobber_mode, ierr=ierr)
        call MPAS_stream_mgr_set_property(manager, streamID, MPAS_STREAM_PROPERTY_IOTYPE, iotype, ierr=ierr)
    else
        call MPAS_stream_mgr_create_stream(manager, streamID, direction, filename, realPrecision=prec, &
                                           clobberMode=clobber_mode, ioType=iotype, ierr=ierr)
    end if

    if (reference_time /= 'initial_time') then
        call MPAS_stream_mgr_set_property(manager, streamID, MPAS_STREAM_PROPERTY_REF_TIME, reference_time, ierr=ierr)
    end if

    if (record_interval /= 'none') then
        call MPAS_stream_mgr_set_property(manager, streamID, MPAS_STREAM_PROPERTY_RECORD_INTV, record_interval, ierr=ierr)
    end if

    if (trim(filename_interval) /= 'none') then
        call MPAS_stream_mgr_set_property(manager, streamID, MPAS_STREAM_PROPERTY_FILENAME_INTV, filename_interval, ierr=ierr)
    end if

    if (ierr == MPAS_STREAM_MGR_NOERR) then
        ierr_c = 0
    else
        ierr_c = 1
    end if

end subroutine stream_mgr_create_stream_c !}}}


subroutine stream_mgr_add_pool_c(manager_c, streamID_c, poolName_c, packages_c, ierr_c) bind(c)!{{{

    use mpas_c_interfacing, only : mpas_c_to_f_string
    use iso_c_binding, only : c_char, c_int, c_ptr, c_f_pointer
    use mpas_derived_types, only : MPAS_streamManager_type, MPAS_STREAM_MGR_NOERR
    use mpas_stream_manager, only : MPAS_stream_mgr_add_pool
    use mpas_kind_types, only : StrKIND

    implicit none

    type (c_ptr) :: manager_c
    character(kind=c_char) :: streamID_c(*)
    character(kind=c_char) :: poolName_c(*)
    character(kind=c_char) :: packages_c(*)
    integer(kind=c_int) :: ierr_c

    type (MPAS_streamManager_type), pointer :: manager
    character(len=StrKIND) :: streamID, poolName, packages
    integer :: ierr


    call c_f_pointer(manager_c, manager)
    call mpas_c_to_f_string(streamID_c, streamID)
    call mpas_c_to_f_string(poolName_c, poolName)
    call mpas_c_to_f_string(packages_c, packages)

    if (len_trim(packages) > 0) then
       call MPAS_stream_mgr_add_pool(manager, streamID, poolName, packages=packages, ierr=ierr)
    else
       call MPAS_stream_mgr_add_pool(manager, streamID, poolName, ierr=ierr)
    end if

    if (ierr == MPAS_STREAM_MGR_NOERR) then
        ierr_c = 0
    else
        ierr_c = 1
    end if

end subroutine stream_mgr_add_pool_c!}}}


subroutine stream_mgr_add_field_c(manager_c, streamID_c, fieldName_c, packages_c, ierr_c) bind(c) !{{{

    use mpas_c_interfacing, only : mpas_c_to_f_string
    use iso_c_binding, only : c_char, c_int, c_ptr, c_f_pointer
    use mpas_derived_types, only : MPAS_streamManager_type, MPAS_STREAM_MGR_NOERR
    use mpas_stream_manager, only : MPAS_stream_mgr_add_field
    use mpas_kind_types, only : StrKIND

    implicit none

    type (c_ptr) :: manager_c
    character(kind=c_char) :: streamID_c(*)
    character(kind=c_char) :: fieldName_c(*)
    character(kind=c_char) :: packages_c(*)
    integer(kind=c_int) :: ierr_c

    type (MPAS_streamManager_type), pointer :: manager
    character(len=StrKIND) :: streamID, fieldName, packages
    integer :: ierr


    call c_f_pointer(manager_c, manager)
    call mpas_c_to_f_string(streamID_c, streamID)
    call mpas_c_to_f_string(fieldName_c, fieldName)
    call mpas_c_to_f_string(packages_c, packages)

    if (len_trim(packages) > 0) then
       call MPAS_stream_mgr_add_field(manager, streamID, fieldName, packages=packages, ierr=ierr)
    else
       call MPAS_stream_mgr_add_field(manager, streamID, fieldName, ierr=ierr)
    end if

    if (ierr == MPAS_STREAM_MGR_NOERR) then
        ierr_c = 0
    else
        ierr_c = 1
    end if

end subroutine stream_mgr_add_field_c !}}}


subroutine stream_mgr_add_immutable_stream_fields_c(manager_c, streamID_c, refStreamID_c, packages_c, ierr_c) bind(c) !{{{

    use mpas_c_interfacing, only : mpas_c_to_f_string
    use iso_c_binding, only : c_char, c_int, c_ptr, c_f_pointer
    use mpas_derived_types, only : MPAS_streamManager_type, MPAS_STREAM_MGR_NOERR, MPAS_STREAM_PROPERTY_IMMUTABLE
    use mpas_stream_manager, only : MPAS_stream_mgr_add_stream_fields, MPAS_stream_mgr_get_property
    use mpas_kind_types, only : StrKIND

    implicit none

    type (c_ptr) :: manager_c
    character(kind=c_char) :: streamID_c(*)    !< stream to add fields to
    character(kind=c_char) :: refStreamID_c(*) !< stream to supply list of fields to add
    character(kind=c_char) :: packages_c(*)
    integer(kind=c_int) :: ierr_c

    type (MPAS_streamManager_type), pointer :: manager
    character(len=StrKIND) :: streamID, refStreamID, packages
    logical :: is_immutable
    integer :: ierr


    call c_f_pointer(manager_c, manager)
    call mpas_c_to_f_string(streamID_c, streamID)
    call mpas_c_to_f_string(refStreamID_c, refStreamID)
    call mpas_c_to_f_string(packages_c, packages)


    call MPAS_stream_mgr_get_property(manager, refStreamID, MPAS_STREAM_PROPERTY_IMMUTABLE, is_immutable, ierr=ierr)
    if (.not. is_immutable) then
       if (ierr == MPAS_STREAM_MGR_NOERR) then
          ierr_c = 0
       else
          ierr_c = 1
       end if
       return ! This stream is not immutable so do not continue.
    endif

    if (len_trim(packages) > 0) then
       call MPAS_stream_mgr_add_stream_fields(manager, streamID, refStreamID, packages=packages, ierr=ierr)
    else
       call MPAS_stream_mgr_add_stream_fields(manager, streamID, refStreamID, ierr=ierr)
    end if

    if (ierr == MPAS_STREAM_MGR_NOERR) then
        ierr_c = 0
    else
        ierr_c = 1
    end if

end subroutine stream_mgr_add_immutable_stream_fields_c !}}}


subroutine stream_mgr_add_alarm_c(manager_c, streamID_c, direction_c, alarmTime_c, alarmInterval_c, ierr_c) bind(c) !{{{

    use mpas_c_interfacing, only : mpas_c_to_f_string
    use iso_c_binding, only : c_char, c_int, c_ptr, c_f_pointer
    use mpas_derived_types, only : MPAS_streamManager_type, MPAS_Clock_type, MPAS_Time_type, MPAS_TimeInterval_type, &
                                   MPAS_STREAM_MGR_NOERR, MPAS_STREAM_INPUT, MPAS_STREAM_OUTPUT, MPAS_START_TIME, &
                                   MPAS_STOP_TIME, MPAS_STREAM_PROPERTY_REF_TIME
    use mpas_stream_manager, only : MPAS_stream_mgr_get_clock, MPAS_stream_mgr_add_alarm, MPAS_stream_mgr_set_property
    use mpas_kind_types, only : StrKIND
    use mpas_timekeeping, only : mpas_add_clock_alarm, mpas_set_time, mpas_get_time, mpas_set_timeInterval, mpas_get_clock_time

    implicit none

    type (c_ptr) :: manager_c
    character(kind=c_char) :: streamID_c(*)
    character(kind=c_char) :: direction_c(*)
    character(kind=c_char) :: alarmTime_c(*)
    character(kind=c_char) :: alarmInterval_c(*)
    integer(kind=c_int) :: ierr_c

    type (MPAS_streamManager_type), pointer :: manager
    type (MPAS_Clock_type), pointer :: clock
    character(len=StrKIND) :: streamID, direction, alarmID, alarmTime, alarmInterval, alarmString
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
        return
    end if

    if (trim(direction) == 'input') then
        idirection = MPAS_STREAM_INPUT
    else if (trim(direction) == 'output') then
        idirection = MPAS_STREAM_OUTPUT
    end if

    call MPAS_stream_mgr_get_clock(manager, clock)

    if (trim(alarmTime) == 'start' .and. trim(alarmInterval) == 'final_only') then
        alarmTime_local = mpas_get_clock_time(clock, MPAS_STOP_TIME, ierr=ierr)
        call mpas_get_time(alarmTime_local, dateTimeString=alarmString)
        call MPAS_stream_mgr_set_property(manager, streamID, MPAS_STREAM_PROPERTY_REF_TIME, alarmString, ierr=ierr)
    else if (trim(alarmTime) == 'start') then
        alarmTime_local = mpas_get_clock_time(clock, MPAS_START_TIME, ierr=ierr)
    else
        call mpas_set_time(alarmTime_local, dateTimeString=alarmTime)
    end if

    if (trim(alarmInterval) == 'initial_only' .or. trim(alarmInterval) == 'final_only') then
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

end subroutine stream_mgr_add_alarm_c !}}}


subroutine stream_mgr_add_pkg_c(manager_c, streamID_c, package_c, ierr_c) bind(c) !{{{

    use mpas_c_interfacing, only : mpas_c_to_f_string
    use iso_c_binding, only : c_char, c_int, c_ptr, c_f_pointer
    use mpas_derived_types, only : MPAS_streamManager_type, MPAS_STREAM_MGR_NOERR
    use mpas_stream_manager, only : MPAS_stream_mgr_add_pkg
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

    call MPAS_stream_mgr_add_pkg(manager, streamID, package, ierr=ierr)

    if (ierr == MPAS_STREAM_MGR_NOERR) then
        ierr_c = 0
    else
        ierr_c = 1
    end if

end subroutine stream_mgr_add_pkg_c !}}}
