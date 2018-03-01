! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_atm_diagnostics_utils

    use mpas_derived_types, only : MPAS_streamManager_type 

    private

    public :: mpas_atm_diag_utils_init, &
              mpas_atm_diag_utils_finalize, &
              mpas_field_will_be_written, &
              mpas_stream_inclusion_count


    type (MPAS_streamManager_type), pointer :: streamManager

    contains


    !-----------------------------------------------------------------------
    !  routine MPAS_atm_diag_utils_init
    !
    !> \brief Prepares diagnostics utilities for use by diagnostics modules
    !> \author Michael Duda
    !> \date   12 October 2016
    !> \details
    !>  This routine prepares the diagnostics utilities for use by diagnostics
    !>  modules.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_diag_utils_init(stream_mgr)

        use mpas_derived_types, only : MPAS_streamManager_type

        implicit none

        type (MPAS_streamManager_type), target, intent(inout) :: stream_mgr

        streamManager => stream_mgr

    end subroutine mpas_atm_diag_utils_init


    !-----------------------------------------------------------------------
    !  routine MPAS_atm_diag_utils_finalize
    !
    !> \brief Performs cleanup after diagnostics utilities will no longer be used
    !> \author Michael Duda
    !> \date   12 October 2016
    !> \details
    !>  This routine performs cleanup after diagnostics utilities will no 
    !>  longer be used
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_diag_utils_finalize()

        implicit none

        ! Nothing to do here at present...

    end subroutine mpas_atm_diag_utils_finalize


    !-----------------------------------------------------------------------
    !  routine MPAS_field_will_be_written
    !
    !> \brief Decide if a field will be written in next call to mpas_stream_mgr_write
    !> \author Michael Duda
    !> \date   12 October 2016
    !> \details
    !>  This function queries the stream manager to see whether there are any
    !>  streams that contain the field 'fieldName' and whose output alarms
    !>  are also ringing. If so, the function returns .true..
    !>  The assumption is that, between the call to this function and the next
    !>  call to write all streams with mpas_stream_mgr_write(), the stream
    !>  (or streams) containing the named field will not have their alarms
    !>  externally reset.
    !
    !-----------------------------------------------------------------------
    logical function mpas_field_will_be_written(fieldName)

        use mpas_kind_types, only : StrKIND
        use mpas_derived_types, only : MPAS_STREAM_OUTPUT, MPAS_STREAM_INPUT_OUTPUT
        use mpas_stream_manager, only : mpas_stream_mgr_begin_iteration, mpas_stream_mgr_get_next_stream, &
                                        MPAS_stream_mgr_ringing_alarms, mpas_stream_mgr_get_next_field

        implicit none

        character(len=*), intent(in) :: fieldName

        character (len=StrKIND) :: streamNameItr
        character (len=StrKIND) :: fieldNameItr
        integer :: streamDirection
        logical :: streamActive
        logical :: fieldActive
        integer :: ierr

        mpas_field_will_be_written = .false.

        call mpas_stream_mgr_begin_iteration(streamManager)
        do while (mpas_stream_mgr_get_next_stream(streamManager, streamID = streamNameItr, &
                                                  directionProperty = streamDirection, activeProperty = streamActive))

            if (streamActive .and. ( streamDirection == MPAS_STREAM_OUTPUT .or. streamDirection == MPAS_STREAM_INPUT_OUTPUT )) then

                if (MPAS_stream_mgr_ringing_alarms(streamManager, streamID=streamNameItr, &
                                                   direction=MPAS_STREAM_OUTPUT, ierr=ierr)) then

                    call mpas_stream_mgr_begin_iteration(streamManager, streamID=streamNameItr)
                    do while (mpas_stream_mgr_get_next_field(streamManager, streamNameItr, fieldNameItr, isActive=fieldActive))

                        if (fieldActive .and. (fieldNameItr == fieldName)) then
                            mpas_field_will_be_written = .true.
                            return
                        end if

                    end do
                end if

            end if

        end do

    end function mpas_field_will_be_written


    !-----------------------------------------------------------------------
    !  routine MPAS_stream_inclusion_count
    !
    !> \brief Returns the number of streams containing the specified field
    !> \author Michael Duda
    !> \date   18 October 2016
    !> \details
    !>  This function queries the stream manager to determine how many streams
    !>  contain the specified field. The optional argument 'direction' can be
    !>  used to limit the count to only input streams, only output streams, or
    !>  both input and output streams. By default, the function only counts
    !>  streams that have a non-'none' input/output interval, but even streams
    !>  with input_interval or output_interval equal to 'none' can be considered
    !>  by setting the optional argument 'includeInactive' to .true..
    !
    !-----------------------------------------------------------------------
    integer function mpas_stream_inclusion_count(fieldName, direction, includeInactive)

        use mpas_kind_types, only : StrKIND
        use mpas_derived_types, only : MPAS_STREAM_INPUT, MPAS_STREAM_OUTPUT, MPAS_STREAM_INPUT_OUTPUT, &
                                       MPAS_STREAM_PROPERTY_RECORD_INTV
        use mpas_stream_manager, only : mpas_stream_mgr_begin_iteration, mpas_stream_mgr_get_next_stream, &
                                        MPAS_stream_mgr_ringing_alarms, mpas_stream_mgr_get_next_field, &
                                        mpas_stream_mgr_get_property

        implicit none

        character(len=*), intent(in) :: fieldName
        integer, intent(in), optional :: direction
        logical, intent(in), optional :: includeInactive

        character (len=StrKIND) :: streamNameItr
        character (len=StrKIND) :: fieldNameItr
        character (len=StrKIND) :: recordIntervalIn
        character (len=StrKIND) :: recordIntervalOut
        integer :: streamDirection
        logical :: streamActive
        logical :: fieldActive
        integer :: ierr

        integer :: local_direction
        logical :: local_includeInactive


        if (present(direction)) then
            local_direction = direction
        else
            local_direction = MPAS_STREAM_INPUT_OUTPUT
        end if

        if (present(includeInactive)) then
            local_includeInactive = includeInactive
        else
            local_includeInactive = .false.
        end if


        mpas_stream_inclusion_count = 0

        call mpas_stream_mgr_begin_iteration(streamManager)
        STREAM_LOOP: do while (mpas_stream_mgr_get_next_stream(streamManager, streamID = streamNameItr, &
                                                               directionProperty = streamDirection, activeProperty = streamActive))

            call MPAS_stream_mgr_get_property(streamManager, trim(streamNameItr), MPAS_STREAM_PROPERTY_RECORD_INTV, &
                                              recordIntervalIn, direction=MPAS_STREAM_INPUT, ierr=ierr)

            call MPAS_stream_mgr_get_property(streamManager, trim(streamNameItr), MPAS_STREAM_PROPERTY_RECORD_INTV, &
                                              recordIntervalOut, direction=MPAS_STREAM_OUTPUT, ierr=ierr)

            ! Determine whether this stream is "active" for the purposes of consideration here
            if (.not. local_includeInactive) then 
                if (streamActive) then
                    streamActive = ((local_direction == MPAS_STREAM_INPUT .and. trim(recordIntervalIn) /= 'none') &
                                    .or. (local_direction == MPAS_STREAM_OUTPUT .and. trim(recordIntervalOut) /= 'none') &
                                    .or. (local_direction == MPAS_STREAM_INPUT_OUTPUT .and. ((trim(recordIntervalIn) /= 'none') &
                                                                                             .or. (trim(recordIntervalOut) /= 'none'))))
                end if
            else
                streamActive = .true.
            end if

            if (streamActive .and. ((local_direction == MPAS_STREAM_INPUT .and. streamDirection == MPAS_STREAM_INPUT) &
                                    .or. (local_direction == MPAS_STREAM_OUTPUT .and. streamDirection == MPAS_STREAM_OUTPUT) &
                                    .or. (local_direction == MPAS_STREAM_INPUT_OUTPUT) &
                                    .or. (streamDirection == MPAS_STREAM_INPUT_OUTPUT))) then

                call mpas_stream_mgr_begin_iteration(streamManager, streamID=streamNameItr)
                do while (mpas_stream_mgr_get_next_field(streamManager, streamNameItr, fieldNameItr, isActive=fieldActive))

                    if (fieldActive .and. (fieldNameItr == fieldName)) then
                        mpas_stream_inclusion_count = mpas_stream_inclusion_count + 1
                        cycle STREAM_LOOP
                    end if

                end do

            end if

        end do STREAM_LOOP

    end function mpas_stream_inclusion_count


end module mpas_atm_diagnostics_utils
