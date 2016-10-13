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
              mpas_field_will_be_written


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


end module mpas_atm_diagnostics_utils
