module mpas_stream_manager

#define STREAM_DEBUG_WRITE(M) write(stderrUnit,*) (M)

    use mpas_kind_types, only : RKIND
    use mpas_grid_types, only : MPAS_Pool_type
    use mpas_timekeeping, only : MPAS_Clock_type, MPAS_Time_type, MPAS_Alarm_type
    use mpas_io_units, only : stderrUnit

    integer, public, parameter :: MPAS_STREAM_INPUT = 1, &
                                  MPAS_STREAM_OUTPUT = 2, &
                                  MPAS_STREAM_INPUT_OUTPUT = 3, &
                                  MPAS_STREAM_NONE = 4

    integer, public, parameter :: MPAS_STREAM_NOERR = 0

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
              MPAS_stream_mgr_read

    private

    type MPAS_streamManager_type
        integer :: num_streams = 0
    end type MPAS_streamManager_type


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
    subroutine MPAS_stream_mgr_init(manager, clock, allFields, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        type (MPAS_Clock_type), intent(inout) :: clock
        type (MPAS_Pool_type), intent(in) :: allFields
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_init()')

        if (present(ierr)) ierr = MPAS_STREAM_NOERR


    end subroutine MPAS_stream_mgr_init


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
    subroutine MPAS_stream_mgr_finalize(manager, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_finalize()') 

        if (present(ierr)) ierr = MPAS_STREAM_NOERR


    end subroutine MPAS_stream_mgr_finalize


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
    subroutine MPAS_stream_mgr_create_stream(manager, streamID, direction, filename, maxRecords, referenceTime, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: direction
        character (len=*), intent(in) :: filename
        integer, intent(in) :: maxRecords
        type (MPAS_Time_type), intent(in), optional :: referenceTime
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_create_stream()')

        if (present(ierr)) ierr = MPAS_STREAM_NOERR

manager % num_streams = manager % num_streams + 1
write(0,*) 'DEBUG TEST we have this many streams ', manager % num_streams

    end subroutine MPAS_stream_mgr_create_stream


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_destroy_stream
    !
    !> \brief Free all memory associated with a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Blah.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_destroy_stream(manager, streamID, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_destroy_stream()')

        if (present(ierr)) ierr = MPAS_STREAM_NOERR


    end subroutine MPAS_stream_mgr_destroy_stream


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_set_property
    !
    !> \brief Sets a property of a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Blah.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_set_property(manager, streamID, property, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: property
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_set_property()')

        if (present(ierr)) ierr = MPAS_STREAM_NOERR


    end subroutine MPAS_stream_mgr_set_property


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_get_property
    !
    !> \brief Queries the status of a stream property in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Blah.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_get_property(manager, streamID, property, isSet, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: property
        logical, intent(out) :: isSet
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_get_property()')

        if (present(ierr)) ierr = MPAS_STREAM_NOERR

        isSet = .false.

    end subroutine MPAS_stream_mgr_get_property


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_pkg
    !
    !> \brief Attach a package logical to the specified stream.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Blah.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_pkg(manager, streamID, package, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        logical, intent(in), target :: package
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_pkg()')

        if (present(ierr)) ierr = MPAS_STREAM_NOERR


    end subroutine MPAS_stream_mgr_add_pkg


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_remove_pkg
    !
    !> \brief Detaches a package logical from the specified stream.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Blah.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_remove_pkg(manager, streamID, package, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        logical, intent(in), target :: package
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_remove_pkg()')

        if (present(ierr)) ierr = MPAS_STREAM_NOERR


    end subroutine MPAS_stream_mgr_remove_pkg


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_field
    !
    !> \brief Add a field to the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Blah.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_field(manager, streamID, fieldName, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: fieldName
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_field()')

        if (present(ierr)) ierr = MPAS_STREAM_NOERR


    end subroutine MPAS_stream_mgr_add_field


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_remove_field
    !
    !> \brief Remove a field from the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Blah.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_remove_field(manager, streamID, fieldName, ierr)
    
        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        character (len=*), intent(in) :: fieldName
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_remove_field()')

        if (present(ierr)) ierr = MPAS_STREAM_NOERR


    end subroutine MPAS_stream_mgr_remove_field


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_att
    !
    !> \brief Add an attribute to the specified stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Blah.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_att(manager, streamID, att, fieldName, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        real (kind=RKIND), intent(in) :: att
        character (len=*), intent(in) :: fieldName
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_att()')

        if (present(ierr)) ierr = MPAS_STREAM_NOERR


    end subroutine MPAS_stream_mgr_add_att


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_add_alarm
    !
    !> \brief Add an I/O alarm to a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Blah.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_add_alarm(manager, streamID, alarm, alarmID, direction, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        type (MPAS_Alarm_type), target, intent(in) :: alarm
        integer, intent(in) :: alarmID
        integer, intent(in) :: direction
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_add_alarm()')

        if (present(ierr)) ierr = MPAS_STREAM_NOERR


    end subroutine MPAS_stream_mgr_add_alarm


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_mgr_remove_alarm
    !
    !> \brief Remove an I/O alarm from a stream in an MPAS stream manager.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   13 June 2014
    !> \details
    !>  Blah.
    !
    !-----------------------------------------------------------------------
    subroutine MPAS_stream_mgr_remove_alarm(manager, streamID, alarmID, direction, ierr)

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        character (len=*), intent(in) :: streamID
        integer, intent(in) :: alarmID
        integer, intent(in) :: direction
        integer, intent(out), optional :: ierr

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_remove_alarm()')

        if (present(ierr)) ierr = MPAS_STREAM_NOERR


    end subroutine MPAS_stream_mgr_remove_alarm


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

        STREAM_DEBUG_WRITE('-- Called MPAS_stream_mgr_write()')

        if (present(ierr)) ierr = MPAS_STREAM_NOERR


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

        if (present(ierr)) ierr = MPAS_STREAM_NOERR
        if (present(actualWhen)) write(actualWhen,'(a)') '0000-01-01_00:00:00'


    end subroutine MPAS_stream_mgr_read
    
end module mpas_stream_manager


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! C interface routines for building streams at run-time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine stream_mgr_create_stream_c(mgr, str, i) bind(c)

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

end subroutine stream_mgr_create_stream_c
