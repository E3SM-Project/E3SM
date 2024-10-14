!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_forcing
!
!> \brief retreives forcing fields
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  Reads in external data and performs time interpolation on that data
!
!-----------------------------------------------------------------------

module mpas_forcing

#define COMMA ,
#define FORCING_DEBUG_WRITE(M) !call mpas_log_write( M )
#define FORCING_WARNING_WRITE(M) call mpas_log_write( M , messageType=MPAS_LOG_WARN)
#define FORCING_ERROR_WRITE(M) call mpas_log_write( M , messageType=MPAS_LOG_ERR)

  use mpas_derived_types
  use mpas_field_routines
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_io_streams
  use mpas_stream_manager
  use mpas_log
  use mpas_io_units

  implicit none

  private
  public :: &
       mpas_forcing_init_group, &
       mpas_forcing_init_field, &
       mpas_forcing_init_field_data, &
       mpas_forcing_get_forcing, &
       mpas_forcing_get_forcing_time, &
       mpas_forcing_write_restart_times, &
       mpas_advance_forcing_clock

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_forcing_init_group
!
!> \brief Add a forcing group to the forcing group list
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details 
!>  Adds a forcing group to the forcing group list. Each forcing group 
!>  contains a forcing clock that can be cycled and a list of associated 
!>  forcing fields that use the forcing group clock. 'forcingGroupHead' 
!>  is a pointer to the forcing group object. 'forcingGroupName' is the 
!>  name of the forcing group to add to the forcing group list. 'domain' 
!>  is the domain type instance that contains the fields that will be read 
!>  into. 'startTimeStr' is the timestamp of the starting time of the 
!>  forcing group clock. 'forcingCycleStart' is the timestamp of starting 
!>  time of the forcing group clock. Specify 'none' for this if no forcing 
!>  is desired. 'forcingCycleDuration' is the 
!>  timestamp of duration of the forcing clock. 'restart' is true if the 
!>  model is restarting. 'forcingCycleStartInclusive' (default: true) is
!>  true if the start time of the forcing cycle is included in the cycle.
!
!-----------------------------------------------------------------------

  subroutine mpas_forcing_init_group(&!{{{
       forcingGroupHead, &
       forcingGroupName, &
       domain, &
       startTimeStr, &
       forcingCycleStart, &
       forcingCycleDuration, &
       restart, &
       forcingCycleStartInclusive)

    type(mpas_forcing_group_type), pointer :: &
         forcingGroupHead ! !< Input: The forcing group linked list

    character(len=*), intent(in) :: &
         forcingGroupName !< Input: the identifying name of the forcing group to be added

    type(domain_type), target :: &
         domain !< Input: the domain to which data will be put

    character(len=*), intent(in) :: &
         startTimeStr, &      !< Input: the forcing start time
         forcingCycleStart, & !< Input: the forcing cycle start time string
         forcingCycleDuration !< Input: the forcing cycle duration time string

    logical, intent(in) :: &
         restart !< Input: whether this is a restarted run

    logical, optional, intent(in) :: &
         forcingCycleStartInclusive !< Input: whether the start of the cycle is inclusive of the cycle

    type(MPAS_Time_type) :: &
         startTime, & ! start time of the forcing clock
         stopTime     ! stop time of forcing clock - !!!! SHOULDNT BE NEEDED

    type(MPAS_TimeInterval_type) :: &
         timeStep, & ! simulation time step
         zeroDuration, & ! a zero duration interval
         timeInterval ! time interval to get to preforcing time

    type(mpas_forcing_group_type), pointer :: &
         forcingGroupNew ! a local instance of the mpas_forcing_group_type object

    FORCING_DEBUG_WRITE('-- Forcing: mpas_forcing_init: startTimeStr: '//trim(startTimeStr))
    FORCING_DEBUG_WRITE('-- Forcing: mpas_forcing_init: forcingCycleStart: '//trim(forcingCycleStart))
    FORCING_DEBUG_WRITE('-- Forcing: mpas_forcing_init: forcingCycleDuration: '//trim(forcingCycleDuration))

    ! loop through the linked list
    if (.not. associated(forcingGroupHead)) then
       allocate(forcingGroupHead)
       nullify(forcingGroupHead % next)
       forcingGroupNew => forcingGroupHead
    else
       forcingGroupNew => forcingGroupHead
       do while (associated(forcingGroupNew % next))
          if (trim(forcingGroupNew % forcingGroupName) == trim(forcingGroupName)) then
             FORCING_ERROR_WRITE('-- Forcing: forcing group name already exists: '//trim(forcingGroupName))
             call mpas_log_write('Forcing: forcing group name already exists', MPAS_LOG_CRIT)
          endif
          forcingGroupNew => forcingGroupNew % next
       enddo
       allocate(forcingGroupNew % next)
       forcingGroupNew => forcingGroupNew % next
       nullify(forcingGroupNew % next)
    endif

    ! set the forcings name
    forcingGroupNew % forcingGroupName = trim(forcingGroupName)

    ! set the forcing group domain
    forcingGroupNew % domain_ptr => domain

    ! add forcing group to restart counter
    call add_forcing_group_to_restart_counter(forcingGroupNew)

    ! create the forcing clock
    call mpas_set_timeInterval(timeStep, dt=0.0_RKIND)
    call mpas_set_time(stopTime, dateTimeString="9999-12-31_23:59:59") ! shouldnt need to do this!
    if (restart) then
       call get_restart_times(forcingGroupNew, timeStep, stopTime)
    else
       call mpas_set_time(startTime, dateTimeString=startTimeStr)
       call mpas_create_clock(forcingGroupNew % forcingClock, startTime=startTime, timeStep=timeStep, stopTime=stopTime)
    endif

    ! determine if we will cycle the forcing clock
    forcingGroupNew % forcingCycleUse = .true.
    if (trim(forcingCycleStart) == "none") then
       forcingGroupNew % forcingCycleUse = .false.
    endif

    ! set forcing clock cycling
    if (forcingGroupNew % forcingCycleUse) then

       ! set the forcing cycle times
       call MPAS_set_time(forcingGroupNew % forcingCycleStart, dateTimeString=forcingCycleStart)
       call MPAS_set_timeInterval(forcingGroupNew % forcingCycleDuration, timeString=forcingCycleDuration)

       ! set the forcing end time
       call MPAS_set_timeInterval(zeroDuration, DD=0)

       if (forcingGroupNew % forcingCycleDuration .gt. zeroDuration) then
          forcingGroupNew % forcingCycleEnd = forcingGroupNew % forcingCycleStart + forcingGroupNew % forcingCycleDuration
       else if (forcingGroupNew % forcingCycleDuration .le. zeroDuration) then
          ! no cycle so set before cycle start
          call MPAS_set_timeInterval(timeInterval, DD=1000)
          forcingGroupNew % forcingCycleEnd = forcingGroupNew % forcingCycleStart - timeInterval
       endif

       ! start time is inclusive or not (default: start is inclusive)
       forcingGroupNew % forcingCycleStartInclusive = .true.
       if (present(forcingCycleStartInclusive)) then
          forcingGroupNew % forcingCycleStartInclusive = forcingCycleStartInclusive
       endif

       ! set the forcing cycle alarm
       call mpas_add_clock_alarm(&
            forcingGroupNew % forcingClock, forcingGroupNew % forcingCycleAlarmID, forcingGroupNew % forcingCycleEnd)

    endif

  end subroutine mpas_forcing_init_group!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_forcing_init_field
!
!> \brief Add a field to the forcing group
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  Add a individual forcing field to the forcing group. All the fields
!>  in a forcing group use the same cycling forcing clock. 'streammanager'
!>  is the streams manager object that will perform the IO. 'forcingGroupHead'
!>  is the forcing group object that the field will be added to and
!>  'forcingGroupName' is the name of that forcing group. 'forcingName' is
!>  an identifying string for the forcing field. 'forcingStreamID' is the
!>  stream identifier associtiated with the field. In registry, the output
!>  field is called 'fieldname' and is in the 'poolname' pool. The interpolation
!>  for this field will use the 'interpolationType' interpolation type
!>  (currently one of 'linear', 'constant' or 'four_point_polynomial').
!>  'forcingReferenceTimeStr' gives a reference time for the forcing data times
!>  of the input data and 'forcingIntervalStr' (optional) gives the fixed.
!>  time interval between input data times. If 'forcingIntervalStr' is not
!>  specified function pointers must be given that calculate a variable
!>  interval. One is for the interval forward from a given forcing data time
!>  ('variable_interval_forward') and the other is for backwards from that
!>  time ('variable_interval_backward').
!
!-----------------------------------------------------------------------

  subroutine mpas_forcing_init_field(&!{{{
       streamManager, &
       forcingGroupHead, &
       forcingGroupName, &
       forcingName, &
       forcingStreamID, &
       poolname, &
       fieldname, &
       interpolationType, &
       forcingReferenceTimeStr, &
       forcingIntervalStr, &
       forcingInitializationType, &
       variable_interval_forward, &
       variable_interval_backward)

    type (MPAS_streamManager_type), intent(inout) :: &
         streamManager ! the stream manager handle

    type(mpas_forcing_group_type), pointer :: &
         forcingGroupHead !< Input: the forcing group linked list

    character(len=*), intent(in) :: &
         forcingGroupName, &     !< Input: the identifying name of the forcing group
         forcingName, &          !< Input: the identifying name of the individual forcing
         forcingStreamID, &      !< Input: the stream ID attached to the forcing field
         poolname, &             !< Input: the pool name with the output field
         fieldname, &            !< Input: the output field name
         interpolationType, &    !< Input: the interpolation type (e.g. 'linear', 'constant')
         forcingReferenceTimeStr !< Input: a reference time for the forcing times

    character(len=*), intent(in), optional :: &
         forcingIntervalStr, &     !< Input: the interval between forcing times
         forcingInitializationType !< Input: the forcing initialization type

    ! forward variable interval function
    interface
       function variable_interval_forward(currentTime) result(variableInterval)
         use mpas_timekeeping
         type(MPAS_Time_type), intent(in) :: currentTime
         type(MPAS_TimeInterval_type) :: variableInterval
       end function variable_interval_forward
    end interface

    ! backward variable interval function
    interface
       function variable_interval_backward(currentTime) result(variableInterval)
         use mpas_timekeeping
         type(MPAS_Time_type), intent(in) :: currentTime
         type(MPAS_TimeInterval_type) :: variableInterval
       end function variable_interval_backward
    end interface

    optional :: variable_interval_forward
    optional :: variable_interval_backward

    type(mpas_forcing_group_type), pointer :: &
         forcingGroup ! forcing group iterating instance

    type(mpas_forcing_stream_type), pointer :: &
         forcingStreamNew ! the new forcing object pointer

    type(MPAS_Time_type) :: &
         forcingAlarmTime ! forcing alarm time

    type(MPAS_TimeInterval_type) :: &
         forcingAlarmInterval ! forcing alarm interval

    ! loop over forcings linked list
    forcingGroup => forcingGroupHead
    do while (associated(forcingGroup))

       ! only get forcings for correct forcings
       if (trim(forcingGroup % forcingGroupName) == trim(forcingGroupName)) then

          FORCING_DEBUG_WRITE('-- Forcing: mpas_forcing_init_forcing '//trim(trim(forcingGroupName))//" $l" COMMA logicArgs=(/associated(forcingGroup % stream)/))

          ! loop through the stream linked list
          if (.not. associated(forcingGroup % stream)) then
             allocate(forcingGroup % stream)
             nullify(forcingGroup % stream % next)
             forcingStreamNew => forcingGroup % stream
          else
             forcingStreamNew => forcingGroup % stream
             do while (associated(forcingStreamNew % next))
                FORCING_DEBUG_WRITE('-- Forcing: mpas_forcing_init_forcing streams: '//trim(forcingStreamID))
                if (trim(forcingStreamNew % forcingStreamID) == trim(forcingStreamID)) then
                   ! forcing stream already exists
                   FORCING_DEBUG_WRITE('-- Forcing: mpas_forcing_init_forcing stream exists: '//trim(forcingStreamID))
                   call add_forcing_field_to_forcing_stream(&
                        forcingStreamNew, &
                        forcingGroup % domain_ptr, &
                        streamManager, &
                        forcingName, &
                        poolname, &
                        fieldname, &
                        interpolationType, &
                        forcingReferenceTimeStr, &
                        forcingIntervalStr, &
                        forcingInitializationType, &
                        variable_interval_forward, &
                        variable_interval_backward)
                   return
                endif
                forcingStreamNew => forcingStreamNew % next
             enddo
             if (trim(forcingStreamNew % forcingStreamID) == trim(forcingStreamID)) then
                ! forcing stream already exists
                FORCING_DEBUG_WRITE('-- Forcing: mpas_forcing_init_forcing stream exists: '//trim(forcingStreamID))
                call add_forcing_field_to_forcing_stream(&
                     forcingStreamNew, &
                     forcingGroup % domain_ptr, &
                     streamManager, &
                     forcingName, &
                     poolname, &
                     fieldname, &
                     interpolationType, &
                     forcingReferenceTimeStr, &
                     forcingIntervalStr, &
                     forcingInitializationType, &
                     variable_interval_forward, &
                     variable_interval_backward)
                return
             endif
             allocate(forcingStreamNew % next)
             forcingStreamNew => forcingStreamNew % next
             nullify(forcingStreamNew % next)
          endif

          FORCING_DEBUG_WRITE('-- Forcing: mpas_forcing_init_forcing create new stream: '//trim(forcingStreamID))
          call create_new_forcing_stream(&
               forcingStreamNew, &
               forcingGroup, &
               streamManager, &
               forcingName, &
               forcingStreamID, &
               interpolationType, &
               forcingReferenceTimeStr, &
               forcingIntervalStr, &
               forcingInitializationType, &
               variable_interval_forward, &
               variable_interval_backward)

          call add_forcing_field_to_forcing_stream(&
               forcingStreamNew, &
               forcingGroup % domain_ptr, &
               streamManager, &
               forcingName, &
               poolname, &
               fieldname, &
               interpolationType, &
               forcingReferenceTimeStr, &
               forcingIntervalStr, &
               forcingInitializationType, &
               variable_interval_forward, &
               variable_interval_backward)

          return
       endif

       forcingGroup => forcingGroup % next
    end do

  end subroutine mpas_forcing_init_field!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  create_new_forcing_stream
!
!> \brief create a new forcing stream
!> \author Adrian K. Turner, LANL
!> \date 24th November 2014
!> \details When adding a forcing field, if its stream is not a member of
!>  the foricng group yet add a new instance to the linked list
!
!-----------------------------------------------------------------------

  subroutine create_new_forcing_stream(&!{{{
       forcingStreamNew, &
       forcingGroup, &
       streamManager, &
       forcingName, &
       forcingStreamID, &
       interpolationType, &
       forcingReferenceTimeStr, &
       forcingIntervalStr, &
       forcingInitializationType, &
       variable_interval_forward, &
       variable_interval_backward)

    type(mpas_forcing_stream_type), pointer :: &
         forcingStreamNew ! the new forcing stream

    type(mpas_forcing_group_type), pointer :: &
         forcingGroup ! the forcing group linekd list

    type (MPAS_streamManager_type), intent(inout) :: &
         streamManager ! the stream manager

    character(len=*), intent(in) :: &
         forcingName, &          !< Input: the name of the forcing field
         forcingStreamID, &      !< Input: the stream ID attached to the forcing field
         interpolationType, &    !< Input: the interpolation type (e.g. 'linear', 'constant')
         forcingReferenceTimeStr !< Input: a reference time for the forcing times

    character(len=*), intent(in), optional :: &
         forcingIntervalStr, &     !< Input: the interval between forcing times
         forcingInitializationType !< Input: the forcing initialization type

    ! forward variable interval function
    interface
       function variable_interval_forward(currentTime) result(variableInterval)
         use mpas_timekeeping
         type(MPAS_Time_type), intent(in) :: currentTime
         type(MPAS_TimeInterval_type) :: variableInterval
       end function variable_interval_forward
    end interface

    ! backward variable interval function
    interface
       function variable_interval_backward(currentTime) result(variableInterval)
         use mpas_timekeeping
         type(MPAS_Time_type), intent(in) :: currentTime
         type(MPAS_TimeInterval_type) :: variableInterval
       end function variable_interval_backward
    end interface

    optional :: variable_interval_forward  ! forward variable interval function
    optional :: variable_interval_backward ! backward variable interval function

    type(MPAS_Time_type) :: &
         forcingAlarmTime

    type(MPAS_TimeInterval_type) :: &
         forcingAlarmInterval

    type(MPAS_clock_type), pointer :: &
         streamClock

    character(len=strKIND) :: &
         streamAlarmID

    integer :: &
         ierr

    FORCING_DEBUG_WRITE('-- Forcing: create_new_forcing_stream streamID: '//trim(forcingStreamID))

    ! set the forcing stream ID
    forcingStreamNew % forcingStreamID = trim(forcingStreamID)

    FORCING_DEBUG_WRITE('-- Forcing: create_new_forcing_stream interpolationType: '//trim(interpolationType))

    ! interpolation type
    forcingStreamNew % interpolationType = trim(interpolationType)
    call interpolation_time_stencil_info(&
         trim(forcingStreamNew % interpolationType), &
         forcingStreamNew % nTimeStencil, &
         forcingStreamNew % nTimeStencilLower, &
         forcingStreamNew % nTimeStencilUpper)
    allocate(forcingStreamNew % forcingTimes(forcingStreamNew % nTimeStencil))

    ! forcing times definition
    FORCING_DEBUG_WRITE('-- Forcing: create_new_forcing_stream forcingIntervalStr: '//trim(forcingIntervalStr))
    FORCING_DEBUG_WRITE('-- Forcing: create_new_forcing_stream forcingReferenceTimeStr: '//trim(forcingReferenceTimeStr))

    ! set the stream forcing reference time
    call mpas_set_time(forcingStreamNew % forcingReferenceTime, dateTimeString=forcingReferenceTimeStr)

    ! set the stream interval
    if (present(forcingIntervalStr) .and. &
         .not. present(variable_interval_forward) .and. &
         .not. present(variable_interval_backward)) then

       ! constant stream interval
       call mpas_set_timeInterval(forcingStreamNew % forcingIntervalConstant, timeString=forcingIntervalStr)

    else if (.not. present(forcingIntervalStr) .and. &
         present(variable_interval_forward) .and. &
         present(variable_interval_backward)) then

       ! variable stream interval
       forcingStreamNew % variable_interval_forward_ptr  => variable_interval_forward
       forcingStreamNew % variable_interval_backward_ptr => variable_interval_backward

    else
       call mpas_log_write('Forcing: create_new_forcing_stream: incorrect forcing interval specification' // &
                            trim(forcingName), MPAS_LOG_CRIT)
    endif

    ! get initial forcing times
    call get_initial_forcing_times(&
         forcingGroup, &
         forcingStreamNew)

    ! set the forcing alarm
    forcingStreamNew % forcingAlarmID = trim(forcingStreamID)

    forcingAlarmTime = forcingStreamNew % forcingTimes(forcingStreamNew % nTimeStencilLower + 1)
    forcingAlarmInterval = forcing_interval(forcingStreamNew, forcingAlarmTime, .true.)

    FORCING_DEBUG_WRITE('-- Forcing: create_new_forcing_stream forcingGroup % forcingGroupName: '//trim(forcingGroup % forcingGroupName))
    call mpas_add_clock_alarm(&
         forcingGroup % forcingClock, forcingStreamNew % forcingAlarmID, forcingAlarmTime, forcingAlarmInterval)

    ! check no alarms defined on stream
    call mpas_stream_mgr_get_clock(streamManager, streamClock, ierr)
    streamAlarmID = trim(forcingStreamID)//'_input'
    if (mpas_is_alarm_defined(streamClock, streamAlarmID, ierr)) then
       call mpas_log_write('Forcing: create_new_forcing_stream: stream has alarm defined'//trim(streamAlarmID), MPAS_LOG_CRIT)
    endif

    ! set the initialization type
    forcingStreamNew % forcingInitializationType = "default"
    if (present(forcingInitializationType)) then
       if (trim(forcingInitializationType) == "default" .or. &
           trim(forcingInitializationType) == "noncycled" .or. &
           trim(forcingInitializationType) == "next") then
          forcingStreamNew % forcingInitializationType = trim(forcingInitializationType)
       else
          call mpas_log_write('Forcing: invalid forcing initialization type: ' // trim(forcingInitializationType) // &
                              'for forcing name: ' // trim(forcingName), MPAS_LOG_CRIT )
       endif
    endif

  end subroutine create_new_forcing_stream!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  add_forcing_field_to_forcing_stream
!
!> \brief add a new forcing field to the forcing stream
!> \author Adrian K. Turner, LANL
!> \date 24th November 2014
!> \details add the details of a new forcing field to the linked list
!>  contained within the forcing stream type
!
!-----------------------------------------------------------------------

  subroutine add_forcing_field_to_forcing_stream(&!{{{
       forcingStream, &
       domain, &
       streamManager, &
       forcingName, &
       poolname, &
       fieldname, &
       interpolationType, &
       forcingReferenceTimeStr, &
       forcingIntervalStr, &
       forcingInitializationType, &
       variable_interval_forward, &
       variable_interval_backward)

    type(mpas_forcing_stream_type), pointer :: &
         forcingStream ! the forcing stream to add the new field to

    type(domain_type), pointer :: &
         domain ! the domain

    type (MPAS_streamManager_type), intent(inout) :: &
         streamManager ! the stream manager

    character(len=*), intent(in) :: &
         forcingName, &          !< Input: the identifying name of the individual forcing
         poolname, &             !< Input: the pool name with the output field
         fieldname, &            !< Input: the output field name
         interpolationType, &    !< Input: the interpolation type (e.g. 'linear', 'constant')
         forcingReferenceTimeStr !< Input: a reference time for the forcing times

    character(len=*), intent(in), optional :: &
         forcingIntervalStr, &     !< Input: the interval between forcing times
         forcingInitializationType !< Input: the forcing initialization type

    ! forward variable interval function
    interface
       function variable_interval_forward(currentTime) result(variableInterval)
         use mpas_timekeeping
         type(MPAS_Time_type), intent(in) :: currentTime
         type(MPAS_TimeInterval_type) :: variableInterval
       end function variable_interval_forward
    end interface

    ! backward variable interval function
    interface
       function variable_interval_backward(currentTime) result(variableInterval)
         use mpas_timekeeping
         type(MPAS_Time_type), intent(in) :: currentTime
         type(MPAS_TimeInterval_type) :: variableInterval
       end function variable_interval_backward
    end interface

    optional :: variable_interval_forward
    optional :: variable_interval_backward

    type(MPAS_TimeInterval_type) :: &
         forcingIntervalConstant

    type(MPAS_Time_type) :: &
         forcingReferenceTime

    type(mpas_forcing_field_type), pointer :: &
         forcingFieldNew

    ! check interpolationType
    if (trim(forcingStream % interpolationType) /= trim(interpolationType)) then
       call mpas_log_write('Forcing: add_forcing_field_to_forcing_stream: incompatible interpolation types in stream: ' // &
           trim(forcingStream % forcingStreamID), MPAS_LOG_CRIT )
    endif

    ! check forcingReferenceTime
    call mpas_set_time(forcingReferenceTime, dateTimeString=forcingReferenceTimeStr)
    if (forcingStream % forcingReferenceTime /= forcingReferenceTime) then
       call mpas_log_write('Forcing: add_forcing_field_to_forcing_stream: incompatible forcing reference times in stream: ' &
                            //trim(forcingStream % forcingStreamID), MPAS_LOG_CRIT)
    endif

    ! check forcing interval optional arguments
    if (present(forcingIntervalStr)) then

       ! constant interval present in forcing specification
       if (present(variable_interval_forward) .or. present(variable_interval_backward)) then
          call mpas_log_write('Forcing: add_forcing_field_to_forcing_stream: overspecified forcing interval: ' &
                               //trim(forcingName), MPAS_LOG_CRIT)
       endif
       if (associated(forcingStream % variable_interval_forward_ptr)) then
          call mpas_log_write('Forcing: add_forcing_field_to_forcing_stream: incompatible forcing interval specification: ' &
                               //trim(forcingName), MPAS_LOG_CRIT)
       endif
       call mpas_set_timeInterval(forcingIntervalConstant, timeString=forcingIntervalStr)
       if (forcingStream % forcingIntervalConstant /= forcingIntervalConstant) then
          call mpas_log_write('Forcing: add_forcing_field_to_forcing_stream: incompatible constant forcing intervals: ' &
                               //trim(forcingName), MPAS_LOG_CRIT)
       endif

    else

       ! variable interval present in forcing specification
       if (.not. present(variable_interval_forward) .or. .not. present(variable_interval_backward)) then
          call mpas_log_write('Forcing: add_forcing_field_to_forcing_stream: underspecified forcing interval: ' &
                               //trim(forcingName), MPAS_LOG_CRIT)
       endif
       if (.not. associated(forcingStream % variable_interval_forward_ptr)) then
          call mpas_log_write('Forcing: add_forcing_field_to_forcing_stream: incompatible forcing interval specification: ' &
                               //trim(forcingName), MPAS_LOG_CRIT)
       endif
       forcingStream % variable_interval_forward_test_ptr  => variable_interval_forward
       forcingStream % variable_interval_backward_test_ptr => variable_interval_backward
       if (.not. associated(forcingStream % variable_interval_forward_test_ptr, &
                            forcingStream % variable_interval_forward_ptr) .or. &
           .not. associated(forcingStream % variable_interval_backward_test_ptr, &
                            forcingStream % variable_interval_backward_ptr)) then
          call mpas_log_write('Forcing: add_forcing_field_to_forcing_stream: incompatible non-constant forcing intervals: ' &
                               //trim(forcingName), MPAS_LOG_CRIT)
       endif

    endif

    ! loop through the stream linked list
    if (.not. associated(forcingStream % field)) then
       allocate(forcingStream % field)
       nullify(forcingStream % field % next)
       forcingFieldNew => forcingStream % field
    else
       forcingFieldNew => forcingStream % field
       do while (associated(forcingFieldNew % next))
          if (trim(forcingFieldNew % forcingName) == trim(forcingName)) then
             ! forcing already exists
             call mpas_log_write('Forcing: add_forcing_field_to_forcing_stream: field ' //trim(forcingStream % forcingStreamID) &
                                  // ' already exists: ' // trim(forcingName), MPAS_LOG_CRIT)
          endif
          forcingFieldNew => forcingFieldNew % next
       enddo
       allocate(forcingFieldNew % next)
       forcingFieldNew => forcingFieldNew % next
       nullify(forcingFieldNew % next)
    endif

    ! check the initialization type
    if (present(forcingInitializationType)) then
       if (trim(forcingInitializationType) /= trim(forcingStream % forcingInitializationType)) then
          FORCING_ERROR_WRITE('-- Forcing: add_forcing_field_to_forcing_stream: inconsistent initialzation type for stream')
          FORCING_ERROR_WRITE('-- Forcing: stream initialization type: '//trim(forcingStream % forcingInitializationType))
          FORCING_ERROR_WRITE('-- Forcing: field initialization type: '//trim(forcingInitializationType))
          FORCING_ERROR_WRITE('-- Forcing: for forcing name: '//trim(forcingName))
          call mpas_log_write('Forcing: inconsistent initialization type', MPAS_LOG_CRIT)
       endif
    endif

    ! set the forcing name
    forcingFieldNew % forcingName = trim(forcingName)

    ! set the field names
    forcingFieldNew % poolname  = trim(poolname)
    forcingFieldNew % fieldname = trim(fieldname)

    ! set up the input fields
    call setup_input_fields(&
         forcingFieldNew, &
         forcingStream, &
         streamManager, &
         domain)

  end subroutine add_forcing_field_to_forcing_stream!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  forcing_interval
!
!> \brief returns the forcing interval given a forcing time
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  Given a forcing time this returns the forcing interval either
!>  forwards or backwards to the next forcing time.
!
!-----------------------------------------------------------------------

  function forcing_interval(forcingStream, currentTime, forward) result(forcingInterval)!{{{

    type(mpas_forcing_stream_type), intent(in) :: &
         forcingStream !< Input: forcing object

    type(MPAS_Time_type), intent(in) :: &
         currentTime !< Input: the forcing time

    logical, optional, intent(in) :: &
         forward !< Input: is the interval forward or backwards from the forcing time

    type(MPAS_TimeInterval_type) :: &
         forcingInterval ! the output forcing interval

    logical :: &
         forward_use

    forward_use = .true.
    if (present(forward)) forward_use = forward

    if (.not. associated(forcingStream % variable_interval_forward_ptr)) then
       ! use constant interval
       forcingInterval = forcingStream % forcingIntervalConstant
    else
       ! use variable interval
       if (forward_use) then
          forcingInterval = forcingStream % variable_interval_forward_ptr(currentTime)
       else
          forcingInterval = forcingStream % variable_interval_backward_ptr(currentTime)
       endif
    end if

  end function forcing_interval!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_initial_forcing_times
!
!> \brief get the initial forcing times
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  Given the forcing reference time and interval determine the initial
!>  forcing times given the current forcing clock time. This searches
!>  for the initial times so reference time is not restricted.
!
!-----------------------------------------------------------------------

  subroutine get_initial_forcing_times(&!{{{
       forcingGroup, &
       forcingStream)

    type(mpas_forcing_group_type), pointer :: &
         forcingGroup ! the forcing group interval

    type(mpas_forcing_stream_type), pointer :: &
         forcingStream ! the individual forcing object

    type(MPAS_Time_type) :: &
         currentTime, & ! the current forcing time from the forcing clock
         forcingTime1, & ! the beginning forcing time of a forcing interval
         forcingTime2    ! the ending forcing time of a forcing interval

    integer :: &
         iTime ! index of forcing times

    ! set the current forcing time
    currentTime = mpas_get_clock_time(forcingGroup % forcingClock, MPAS_NOW)

    ! initialize the forcing data time
    forcingTime1 = forcingStream % forcingReferenceTime
    forcingTime2 = forcingStream % forcingReferenceTime + forcing_interval(forcingStream, forcingStream % forcingReferenceTime)

    if (currentTime .ge. forcingTime1 .and. &
        currentTime .lt. forcingTime2) then
       ! current forcing clock time is in the current interval

       call populate_forcing_times(forcingStream, forcingTime1, forcingStream % nTimeStencilLower)

    else if (currentTime .lt. forcingTime1) then
       ! current forcing clock time earlier than the current interval

       backward_search: do

          forcingTime2 = forcingTime1
          forcingTime1 = forcingTime2 - forcing_interval(forcingStream, forcingTime2, forward=.false.)

          if (currentTime .ge. forcingTime1 .and. &
              currentTime .lt. forcingTime2) then
             call populate_forcing_times(forcingStream, forcingTime1, forcingStream % nTimeStencilLower)
             exit
          endif

       end do backward_search

    else if (currentTime .ge. forcingTime2) then
       ! current forcing clock time later than the current interval

       forward_search: do

          forcingTime1 = forcingTime2
          forcingTime2 = forcingTime1 + forcing_interval(forcingStream, forcingTime1)

          if (currentTime .ge. forcingTime1 .and. &
              currentTime .lt. forcingTime2) then
             call populate_forcing_times(forcingStream, forcingTime1, forcingStream % nTimeStencilLower)
             exit
          endif

       end do forward_search

    endif

  end subroutine get_initial_forcing_times!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  populate_forcing_times
!
!> \brief given one of the forcing times populate the others
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  This subroutine fills the time stencil with forcing times given
!>  one of the forcing times at the forcingSlot slot
!
!-----------------------------------------------------------------------

  subroutine populate_forcing_times(&!{{{
       forcingStream, &
       forcingTimeSlot, &
       forcingSlot)

    type(mpas_forcing_stream_type), pointer :: &
         forcingStream ! forcing object

    type(MPAS_Time_type), intent(in) :: &
         forcingTimeSlot !< Input: the given forcing time

    integer, intent(in) :: &
         forcingSlot !< Input: the slot in the times stencil for forcingTimeSlot

    type(MPAS_Time_type) :: &
         forcingTime ! a forcing time for the other slots

    integer :: &
         iTime ! index of forcing time slots

    character(len=strKIND) :: strout

    ! fill given slot
    forcingStream % forcingTimes(forcingSlot) = forcingTimeSlot

    ! fill out higher
    forcingTime = forcingTimeSlot
    do iTime = forcingSlot+1, forcingStream % nTimeStencil

       forcingTime = forcingTime + forcing_interval(forcingStream, forcingTime)
       forcingStream % forcingTimes(iTime) = forcingTime

    enddo ! iTime

    ! fill out lower
    forcingTime = forcingTimeSlot
    do iTime = forcingSlot-1, 1, -1

       forcingTime = forcingTime - forcing_interval(forcingStream, forcingTime, forward=.false.)
       forcingStream % forcingTimes(iTime) = forcingTime

    enddo ! iTime

  end subroutine populate_forcing_times!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_forcing_init_field_data
!
!> \brief get the initial forcing data
!> \author Adrian K. Turner, LANL
!> \date 24th November 2014
!> \details
!>  Once all the fields have been added to the forcing group object read
!>  in all the data for all initial forcing times for all the fields in
!>  the forcing group. 'forcingGroupHead' is the forcing group object,
!>  'forcingGroupName' is the name of the forcing group to read data into
!>  and 'streamManager' is the stream manager that will perform the IO.
!
!-----------------------------------------------------------------------

  subroutine mpas_forcing_init_field_data(&!{{{
       forcingGroupHead, &
       forcingGroupName, &
       streamManager, &
       restart, &
       interpolateAtInit)

    type(mpas_forcing_group_type), pointer :: &
         forcingGroupHead ! the forcing group linked list

    character(len=*), intent(in) :: &
         forcingGroupName !< Input: the idenifying name of the forcing group

    type (MPAS_streamManager_type), intent(inout) :: &
         streamManager ! the stream manager

    logical, intent(in) :: &
         restart !< Input: whether this is a restarted run

    logical, intent(in) :: &
         interpolateAtInit !< Input: interpolate to forcing variables at time 0.

    type(mpas_forcing_group_type), pointer :: &
         forcingGroup ! forcing group iterator

    type(mpas_forcing_stream_type), pointer :: &
         forcingStream ! forcing stream iterator

    ! loop over forcings linked list
    forcingGroup => forcingGroupHead
    do while (associated(forcingGroup))

       if (trim(forcingGroup % forcingGroupName) == trim(forcingGroupName)) then

          ! loop over streams linked list
          forcingStream => forcingGroup % stream
          do while (associated(forcingStream))

             call get_initial_forcing_data(&
                  forcingGroup, &
                  forcingStream, &
                  streamManager, &
                  restart, &
                  interpolateAtInit)

             forcingStream => forcingStream % next
          end do

       endif

       forcingGroup => forcingGroup % next
    end do

  end subroutine mpas_forcing_init_field_data!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_initial_forcing_data
!
!> \brief read in the initial forcing data
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  Given a set of forcing times read in all the data for each slot in
!>  the time stencil for the fields in one stream
!
!-----------------------------------------------------------------------

  subroutine get_initial_forcing_data(&!{{{
       forcingGroup, &
       forcingStream, &
       streamManager, &
       restart, &
       interpolateAtInit)

    type(mpas_forcing_group_type), pointer :: &
         forcingGroup ! the forcing group object

    type(mpas_forcing_stream_type), pointer :: &
         forcingStream ! the forcing object

    type (MPAS_streamManager_type), intent(inout) :: &
         streamManager ! stream manager

    logical, intent(in) :: &
         restart !< Input: whether this is a restarted run

    logical, intent(in) :: &
         interpolateAtInit !< Input: interpolate to forcing variables at time 0.

    type(MPAS_Time_type) :: &
         forcingTimeCycle, & ! the forcing time after cycling
         currentForcingTime

    character(len=strKIND) :: &
         forcingTimeStr ! the timestamp of the forcing time

    integer :: &
         iTime, & ! index of the forcing time stops
         iTime2, &
         ierr

    FORCING_DEBUG_WRITE('-------------------------------------------')

    ! loop over the forcing time slots
    do iTime = 1, forcingStream % nTimeStencil

       ! get the forcing time of the slot
       forcingTimeCycle = forcingStream % forcingTimes(iTime)

       ! cycle this clock time to ensure its in the correct range
       if (forcingGroup % forcingCycleUse .and. &
            trim(forcingStream % forcingInitializationType) == "default") then

          if (forcingGroup % forcingCycleStartInclusive) then

             ! cycle start time is inclusive to the cycle

             if (forcingTimeCycle .lt. forcingGroup % forcingCycleStart) then
                do
                   forcingTimeCycle = forcingTimeCycle + forcingGroup % forcingCycleDuration
                   if (forcingTimeCycle .ge. forcingGroup % forcingCycleStart .and. &
                        forcingTimeCycle .lt. forcingGroup % forcingCycleEnd) exit
                enddo
             endif

             if (forcingTimeCycle .ge. forcingGroup % forcingCycleEnd) then
                do
                   forcingTimeCycle = forcingTimeCycle - forcingGroup % forcingCycleDuration
                   if (forcingTimeCycle .ge. forcingGroup % forcingCycleStart .and. &
                        forcingTimeCycle .lt. forcingGroup % forcingCycleEnd) exit
                enddo
             endif

          else

             ! cycle start time is not inclusive to the cycle

             if (forcingTimeCycle .le. forcingGroup % forcingCycleStart) then
                do
                   forcingTimeCycle = forcingTimeCycle + forcingGroup % forcingCycleDuration
                   if (forcingTimeCycle .gt. forcingGroup % forcingCycleStart .and. &
                        forcingTimeCycle .le. forcingGroup % forcingCycleEnd) exit
                enddo
             endif

             if (forcingTimeCycle .gt. forcingGroup % forcingCycleEnd) then
                do
                   forcingTimeCycle = forcingTimeCycle - forcingGroup % forcingCycleDuration
                   if (forcingTimeCycle .gt. forcingGroup % forcingCycleStart .and. &
                        forcingTimeCycle .le. forcingGroup % forcingCycleEnd) exit
                enddo
             endif

          endif

       endif

       ! other initialization types
       if (trim(forcingStream % forcingInitializationType) == "next") then

          if (forcingGroup % forcingCycleStartInclusive) then

             if (forcingTimeCycle .lt. forcingGroup % forcingCycleStart) then

                do iTime2 = 1, forcingStream % nTimeStencil
                   if (forcingStream % forcingTimes(iTime2) .ge. forcingGroup % forcingCycleStart) then
                      forcingTimeCycle = forcingStream % forcingTimes(iTime2)
                      exit
                   endif
                enddo ! iTime2

             endif

          else

             if (forcingTimeCycle .le. forcingGroup % forcingCycleStart) then

                do iTime2 = 1, forcingStream % nTimeStencil
                   if (forcingStream % forcingTimes(iTime2) .gt. forcingGroup % forcingCycleStart) then
                      forcingTimeCycle = forcingStream % forcingTimes(iTime2)
                      exit
                   endif
                enddo ! iTime2

             endif

          endif

       endif

       ! load the data into the slot
       call mpas_get_time(forcingTimeCycle, dateTimeString=forcingTimeStr)
       FORCING_DEBUG_WRITE('-- Forcing: get_initial_forcing_data: READ ATTEMPT: '//trim(forcingTimeStr))
       call MPAS_stream_mgr_read(&
            streamManager, &
            forcingStream % forcingStreamID, &
            timeLevel=iTime, &
            when=forcingTimeStr, &
            whence=MPAS_STREAM_EXACT_TIME, &
            rightNow=.true., &
            ierr=ierr)

       if (ierr /= 0) then
          call mpas_log_write('Forcing: get_initial_forcing_data: READ: $i', MPAS_LOG_CRIT, intArgs=(/ierr/))
       endif

       FORCING_DEBUG_WRITE('-- Forcing: get_initial_forcing_data: $i '//trim(forcingStream % forcingStreamID)//" "//trim(forcingTimeStr) COMMA intArgs=(/iTime/))

    enddo ! iTime

    if (interpolateAtInit .and. .not. restart) then

       ! interpolate data
       currentForcingTime = mpas_get_clock_time(forcingGroup % forcingClock, MPAS_NOW)

       call forcing_data_interpolation(&
            forcingGroup % domain_ptr, &
            forcingStream, &
            currentForcingTime)

    endif

  end subroutine get_initial_forcing_data!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_forcing_get_forcing
!
!> \brief Get the forcing data and do interpolation for a forcing group
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  Loop over all the individual forcing fields  in the forcing group
!>  and get the data and perform the time interpolation.
!>  'forcingGroupHead' is the forcing group object, 'forcingGroupName'
!>  is the name of the forcing group to read data into, 'streamManager'
!>  is the stream manager that will perform the IO and 'dt' is the current
!>  time step duration.
!
!-----------------------------------------------------------------------

  subroutine mpas_forcing_get_forcing(&!{{{
       forcingGroupHead, &
       forcingGroupName, &
       streamManager, &
       dt)

    type(mpas_forcing_group_type), pointer :: &
         forcingGroupHead ! the forcing group linked list

    character(len=*), intent(in) :: &
         forcingGroupName !< Input: the idenifying name of the forcing group

    type (MPAS_streamManager_type), intent(inout) :: &
         streamManager ! stream manager

    real(kind=RKIND), intent(in) :: &
         dt !< Input: the current time step

    type(mpas_forcing_group_type), pointer :: &
         forcingGroup ! the forcing group linked list

    type(mpas_forcing_stream_type), pointer :: &
         forcingStream ! forcing stream iterator

    FORCING_DEBUG_WRITE('-- Forcing: mpas_forcing_get_forcings: forcingGroup name: '//trim(forcingGroupName))

    ! loop over forcings linked list
    forcingGroup => forcingGroupHead
    do while (associated(forcingGroup))

       ! only get forcings for correct forcings
       if (trim(forcingGroup % forcingGroupName) == trim(forcingGroupName)) then

          ! advance the forcing time
          call mpas_advance_forcing_clock(forcingGroup, dt)

          ! cycle the forcing clock
          if (forcingGroup % forcingCycleUse) then
             if (forcingGroup % forcingCycleStartInclusive) &
                  call cycle_forcing_clock(forcingGroup)
          endif

          ! loop over individual forcing fields
          forcingStream => forcingGroup % stream
          do while (associated(forcingStream))

             ! get the individual forcing
             call get_forcing(&
                  forcingGroup, &
                  forcingStream, &
                  streamManager)

             forcingStream => forcingStream % next
          enddo

          ! cycle the forcing clock
          if (forcingGroup % forcingCycleUse) then
             if (.not. forcingGroup % forcingCycleStartInclusive) &
                  call cycle_forcing_clock(forcingGroup)
          endif

          exit
       endif

       forcingGroup => forcingGroup % next
    end do

  end subroutine mpas_forcing_get_forcing!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_advance_forcing_clock
!
!> \brief set the forcing clock
!> \author Adrian K. Turner, LANL
!> \date September 29th 2014
!> \details
!>  Before doing the forcing read and interpolation increment the forcing
!>  group clock forward by the timestep
!
!-----------------------------------------------------------------------

  subroutine mpas_advance_forcing_clock(&!{{{
       forcingGroup, &
       dt)

    type(mpas_forcing_group_type), pointer :: &
         forcingGroup ! the forcing group object

    real(kind=RKIND), intent(in) :: &
         dt !< Input: the simulation time step

    type(MPAS_TimeInterval_type) :: &
         timeStep ! time step interval

    ! increment clock with timestep
    call mpas_set_timeInterval(timeStep, dt=dt)
    !call mpas_log_write('time step from mpas_set_timeInterval: $i $i $i', intArgs=(/int(timeStep%ti%basetime%S), int(timeStep%ti%basetime%Sn), int(timeStep%ti%basetime%Sd) /))
    call mpas_advance_clock(forcingGroup % forcingClock, timeStep)

  end subroutine mpas_advance_forcing_clock!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cycle_forcing_clock
!
!> \brief cycle the forcing clock
!> \author Adrian K. Turner, LANL
!> \date September 29th 2014
!> \details
!>  Check to see if it needs cycling and cycle it if it does
!
!-----------------------------------------------------------------------

  subroutine cycle_forcing_clock(&!{{{
       forcingGroup)

    type(mpas_forcing_group_type), pointer :: &
         forcingGroup ! the forcing group object

    type(MPAS_Time_type) :: &
         oldForingClockTime, & ! the initial forcing clock time
         newForcingClockTime ! the forcing clock time after cycling

    integer :: &
         iTime ! forcing times index

    type(mpas_forcing_stream_type), pointer :: &
         forcingStream

    ! check if cycling alarm is ringing
    if (mpas_is_alarm_ringing(forcingGroup % forcingClock, forcingGroup % forcingCycleAlarmID)) then

       FORCING_DEBUG_WRITE('-- Forcing: set_forcing_clock: cycle forcing clock')

       ! if ringing cycle the clock
       oldForingClockTime = mpas_get_clock_time(forcingGroup % forcingClock, MPAS_NOW)
       newForcingClockTime = oldForingClockTime - forcingGroup % forcingCycleDuration
       call mpas_set_clock_time(forcingGroup % forcingClock, newForcingClockTime, MPAS_NOW)

       ! if ringing cycle the forcing times
       forcingStream => forcingGroup % stream
       do while (associated(forcingStream))

          ! cycle all the current forcing times
          do iTime = 1, forcingStream % nTimeStencil
             forcingStream % forcingTimes(iTime) = forcingStream % forcingTimes(iTime) - forcingGroup % forcingCycleDuration
          enddo ! iTime

          forcingStream => forcingStream % next
       enddo

    endif ! forcingCycleAlarmID ringing

  end subroutine cycle_forcing_clock!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_forcing_get_forcing_time
!
!> \brief return the current forcing time
!> \author Adrian K. Turner, LANL
!> \date 24th November 2014
!> \details
!>  Return the current forcing clock time for a forcing group.
!>  'forcingGroupHead' is the forcing group object, 'forcingGroupName'
!>  is the name of the forcing group to get the time for and
!>  'forcingTime' is the output forcing group current clock time.
!
!-----------------------------------------------------------------------

  subroutine mpas_forcing_get_forcing_time(&!{{{
       forcingGroupHead, &
       forcingGroupName, &
       forcingTime)

    type(mpas_forcing_group_type), pointer :: &
         forcingGroupHead ! the forcing group linked list

    character(len=*), intent(in) :: &
         forcingGroupName !< Input: the idenifying name of the forcing group

    type(MPAS_Time_type), intent(out) :: &
         forcingTime ! the current forcing time for the forcing group

    type(mpas_forcing_group_type), pointer :: &
         forcingGroup ! forcing group iterator

    ! loop over forcings linked list
    forcingGroup => forcingGroupHead
    do while (associated(forcingGroup))

       ! only get forcings for correct forcings
       if (trim(forcingGroup % forcingGroupName) == trim(forcingGroupName)) then

          forcingTime = mpas_get_clock_time(forcingGroup % forcingClock, MPAS_NOW)

       endif

       forcingGroup => forcingGroup % next
    end do

  end subroutine mpas_forcing_get_forcing_time!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_forcing
!
!> \brief Perform the forcing on the individual forcing stream
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  Get the data for the new forcing time and perform the time
!>  interpolation for all the fields in a forcing stream
!
!-----------------------------------------------------------------------

  subroutine get_forcing(&!{{{
       forcingGroup, &
       forcingStream, &
       streamManager)

    type(mpas_forcing_group_type), pointer :: &
         forcingGroup ! the forcing group object

    type(mpas_forcing_stream_type), pointer :: &
         forcingStream ! the forcing object

    type (MPAS_streamManager_type), intent(inout) :: &
         streamManager ! the stream manager

    type(MPAS_Time_type) :: &
         currentTime, & ! the current forcing time
         forcingTimeNew ! the next forcing time

    character(len=strKIND) :: &
         forcingTimeCycleStr, & ! the cycled forcing time
         forcingTimeStr         ! time string of forcing

    integer :: &
         iTime, & ! forcing times index
         ierr ! status of stream read

    FORCING_DEBUG_WRITE('-- Forcing: mpas_forcing_get_forcing: $l' COMMA logicArgs=(/mpas_is_alarm_ringing(forcingGroup % forcingClock, forcingStream % forcingAlarmID)/))

    ! get current forcing time
    currentTime = mpas_get_clock_time(forcingGroup % forcingClock, MPAS_NOW)

    ! check alarm to see if ringing for crossed forcing time
    if (mpas_is_alarm_ringing(forcingGroup % forcingClock, forcingStream % forcingAlarmID)) then
       ! need to load new data

       FORCING_DEBUG_WRITE('-- Forcing: mpas_forcing_get_forcing: READ')

       ! reset the forcing alarm
       call mpas_reset_clock_alarm(forcingGroup % forcingClock, forcingStream % forcingAlarmID)

       ! swap times
       do iTime = 1, forcingStream % nTimeStencil-1
          forcingStream % forcingTimes(iTime) = forcingStream % forcingTimes(iTime+1)
       enddo ! iTime

       ! shift data
       call forcing_shift_data(&
            forcingGroup % domain_ptr, &
            forcingStream)

       ! determine new time to load data from
       forcingTimeNew = forcingStream % forcingTimes(forcingStream % nTimeStencil-1) + &
            forcing_interval(forcingStream, forcingStream % forcingTimes(forcingStream % nTimeStencil-1))

       ! add the new time to the final time slot
       forcingStream % forcingTimes(forcingStream % nTimeStencil) = forcingTimeNew

       ! cycle this new time if needed
       if (forcingGroup % forcingCycleUse) then
          if ( (      forcingGroup % forcingCycleStartInclusive .and. forcingTimeNew .ge. forcingGroup % forcingCycleEnd) .or. &
               (.not. forcingGroup % forcingCycleStartInclusive .and. forcingTimeNew .gt. forcingGroup % forcingCycleEnd) ) then
             forcingTimeNew = forcingTimeNew - forcingGroup % forcingCycleDuration
          endif
       endif

       ! load the data into last slot
       call mpas_get_time(forcingTimeNew, dateTimeString=forcingTimeCycleStr)
       call MPAS_stream_mgr_read(&
            streamManager, &
            forcingStream % forcingStreamID, &
            timeLevel=forcingStream % nTimeStencil, &
            when=forcingTimeCycleStr, &
            whence=MPAS_STREAM_EXACT_TIME, &
            rightNow=.true., &
            ierr=ierr)

       if (ierr /= 0) then
          call mpas_log_write('Forcing: mpas_forcing_get_forcing: READ: $i', MPAS_LOG_CRIT, intArgs=(/ierr/))
       endif

       FORCING_DEBUG_WRITE('-- Forcing: mpas_forcing_get_forcing: READ: ' // trim(forcingStream % forcingStreamID) // " " // trim(forcingTimeCycleStr))

    endif ! forcingAlarmID ringing

    ! interpolate data
    call forcing_data_interpolation(&
         forcingGroup % domain_ptr, &
         forcingStream, &
         currentTime)

  end subroutine get_forcing!}}}

!-----------------------------------------------------------------------
! Interpolation
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  interpolation_time_stencil_info
!
!> \brief return time stencil info
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  Given a particular interpolation type return the size of the time
!>  stencil and the number of forcing times earlier and later than the
!>  current forcing time
!
!-----------------------------------------------------------------------

  subroutine interpolation_time_stencil_info(&!{{{
       interpolationType, &
       stencilSize, &
       stencilSizeLower, &
       stencilSizeUpper)

    character(len=*), intent(in) :: &
         interpolationType !< Input: the interpolation type

    integer, intent(out) :: &
         stencilSize, &      !< Output: the time stencil size
         stencilSizeLower, & !< Output: the time stencil size below the current forcing clock time
         stencilSizeUpper    !< Output: the time stencil size above the current forcing clock time

    select case (trim(interpolationType))
    case ("linear")
       stencilSize = 2
       stencilSizeLower = 1
       stencilSizeUpper = 1
    case ("constant")
       stencilSize = 2
       stencilSizeLower = 1
       stencilSizeUpper = 1
    case ("four_point_polynomial")
       stencilSize = 4
       stencilSizeLower = 2
       stencilSizeUpper = 2
    case default
       call mpas_log_write('Forcing: Unknown forcing interpolation type'//trim(interpolationType), MPAS_LOG_CRIT)
    end select

  end subroutine interpolation_time_stencil_info!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_interpolants
!
!> \brief Calculate the interpolants used to perform forcing interpolation
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  Given the current time and forcing times calculate the correct
!>  interpolants given the interpolation type
!
!-----------------------------------------------------------------------

  subroutine get_interpolants(interpolants, forcingStream, currentTime)!{{{

    real(kind=RKIND), dimension(:), intent(out) :: &
         interpolants !< Output: interpolation weights

    type (mpas_forcing_stream_type), intent(in) :: &
         forcingStream !< Input: the forcing object

    type (MPAS_Time_Type), intent(in) :: &
         currentTime !< Input: the current forcing clock time

    select case (trim(forcingStream % interpolationType))
    case ("linear")
       call get_interpolants_linear(interpolants, forcingStream, currentTime)
    case ("constant")
       call get_interpolants_constant(interpolants, forcingStream, currentTime)
    case ("four_point_polynomial")
       call get_interpolants_four_point_polynomial(interpolants, forcingStream, currentTime)
    case default
       call mpas_log_write('Forcing: Unknown forcing interpolation type'//trim(forcingStream % interpolationType), MPAS_LOG_CRIT)
    end select

  end subroutine get_interpolants!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_interpolants_linear
!
!> \brief Calculate the interpolants used to perform forcing interpolation
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  Given the current time and forcing times calculate the correct
!>  interpolants with linear interpolation
!
!-----------------------------------------------------------------------

  subroutine get_interpolants_linear(interpolants, forcingStream, currentTime)!{{{

    real(kind=RKIND), dimension(:), intent(out) :: &
         interpolants !< Output: interpolation weights

    type (mpas_forcing_stream_type), intent(in) :: &
         forcingStream !< Input: the forcing object

    type (MPAS_Time_Type), intent(in) :: &
         currentTime !< Input: the current forcing clock time

    type (MPAS_TimeInterval_Type) :: &
         diff, &  ! time interval between forcing times
         diff1, & ! time interval between current time and earlier forcing time
         diff2    ! time interval between current time and later forcing time

    real(kind=RKIND) :: &
         diffr, &  ! real versions of above
         diffr1, & ! real versions of above
         diffr2    ! real versions of above

    diff = forcingStream % forcingTimes(2) - forcingStream % forcingTimes(1)
    diff1 = currentTime - forcingStream % forcingTimes(1)
    diff2 = forcingStream % forcingTimes(2) - currentTime

    call mpas_get_timeInterval(diff, forcingStream % forcingTimes(1), dt=diffr)
    call mpas_get_timeInterval(diff1, forcingStream % forcingTimes(1), dt=diffr1)
    call mpas_get_timeInterval(diff2, currentTime, dt=diffr2)

    !call mpas_log_write('diffr2 $r, diffr, $r', realArgs=(/ diffr2, diffr /))
    interpolants(1) = diffr2 / diffr
    interpolants(2) = 1.0_RKIND - interpolants(1) !diffr1 / diffr
    

  end subroutine get_interpolants_linear!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_interpolants_constant
!
!> \brief Calculate the interpolants used to perform forcing interpolation
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  Given the current time and forcing times calculate the correct
!>  interpolants with piecewise constant interpolation

!-----------------------------------------------------------------------

  subroutine get_interpolants_constant(interpolants, forcingStream, currentTime)!{{{

    real(kind=RKIND), dimension(:), intent(out) :: &
         interpolants !< Output: interpolation weights

    type (mpas_forcing_stream_type), intent(in) :: &
         forcingStream !< Input: the forcing object

    type (MPAS_Time_Type), intent(in) :: &
         currentTime !< Input: the current forcing clock time

    type (MPAS_TimeInterval_Type) :: &
         diff1, & ! time interval between current time and earlier forcing time
         diff2    ! time interval between current time and later forcing time

    diff1 = currentTime - forcingStream % forcingTimes(1)
    diff2 = forcingStream % forcingTimes(2) - currentTime

    if (diff2 > diff1) then
       interpolants(1) = 1.0_RKIND
       interpolants(2) = 0.0_RKIND
    else
       interpolants(1) = 0.0_RKIND
       interpolants(2) = 1.0_RKIND
    endif

  end subroutine get_interpolants_constant!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_interpolants_four_point_polynomial
!
!> \brief Calculate the interpolants used to perform forcing interpolation
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  Given the forcing times calculate the interpolation weights for
!>  a four point Lagrange polynomial interpolation
!
!-----------------------------------------------------------------------

  subroutine get_interpolants_four_point_polynomial(interpolants, forcingStream, currentTime)!{{{

    real(kind=RKIND), dimension(:), intent(out) :: &
         interpolants !< Output: interpolation weights

    type (mpas_forcing_stream_type), intent(in) :: &
         forcingStream !< Input: the forcing object

    type (MPAS_Time_Type), intent(in) :: &
         currentTime !< Input: the current forcing clock time

    type (MPAS_TimeInterval_Type) :: &
         x0_interval, x1_interval, x2_interval, x3_interval, x4_interval

    real(kind=RKIND) :: &
         x0, x1, x2, x3, x4

    x0_interval = currentTime - currentTime
    x1_interval = forcingStream % forcingTimes(1) - currentTime
    x2_interval = forcingStream % forcingTimes(2) - currentTime
    x3_interval = forcingStream % forcingTimes(3) - currentTime
    x4_interval = forcingStream % forcingTimes(4) - currentTime

    call mpas_get_timeInterval(x0_interval, currentTime, dt=x0)
    call mpas_get_timeInterval(x1_interval, currentTime, dt=x1)
    call mpas_get_timeInterval(x2_interval, currentTime, dt=x2)
    call mpas_get_timeInterval(x3_interval, currentTime, dt=x3)
    call mpas_get_timeInterval(x4_interval, currentTime, dt=x4)

    interpolants(1) = ((x0 - x2) / (x1 - x2)) * ((x0 - x3) / (x1 - x3)) * ((x0 - x4) / (x1 - x4))
    interpolants(2) = ((x0 - x1) / (x2 - x1)) * ((x0 - x3) / (x2 - x3)) * ((x0 - x4) / (x2 - x4))
    interpolants(3) = ((x0 - x1) / (x3 - x1)) * ((x0 - x2) / (x3 - x2)) * ((x0 - x4) / (x3 - x4))
    interpolants(4) = ((x0 - x1) / (x4 - x1)) * ((x0 - x2) / (x4 - x2)) * ((x0 - x3) / (x4 - x3))

  end subroutine get_interpolants_four_point_polynomial!}}}

!-----------------------------------------------------------------------
! arrays
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  setup_input_fields
!
!> \brief creates the input fields
!> \author Adrian K. Turner, LANL
!> \date October 9th, 2014
!> \details
!>  This duplicates the output field with the same number of time levels
!>  as the time stencil. Data is read into here before interpolation.
!
!-----------------------------------------------------------------------

  subroutine setup_input_fields(&!{{{
       forcingField, &
       forcingStream, &
       streamManager, &
       domain)

    type(mpas_forcing_field_type), pointer :: &
         forcingField ! forcing field object

    type(mpas_forcing_stream_type), pointer :: &
         forcingStream ! forcing stream

    type (MPAS_streamManager_type), intent(inout) :: &
         streamManager ! stream manager

    type(domain_type), pointer :: domain

    type (block_type), pointer :: &
         block ! block object

    type (MPAS_pool_type), pointer :: &
         forcingPoolOutput, & ! forcing output pool
         forcingPoolInput, &  ! forcing input pool
         forcingPoolInputPrev

    type (MPAS_pool_field_info_type) :: &
         forcingPoolInfo ! pool info

    ! input and output fields
    type(field0DReal), pointer :: field0DRealInput, field0DRealOutput, field0DRealInputPrev
    type(field1DReal), pointer :: field1DRealInput, field1DRealOutput, field1DRealInputPrev
    type(field2DReal), pointer :: field2DRealInput, field2DRealOutput, field2DRealInputPrev
    type(field3DReal), pointer :: field3DRealInput, field3DRealOutput, field3DRealInputPrev
    type(field4DReal), pointer :: field4DRealInput, field4DRealOutput, field4DRealInputPrev
    type(field5DReal), pointer :: field5DRealInput, field5DRealOutput, field5DRealInputPrev

    type(field0DInteger), pointer :: field0DIntInput, field0DIntOutput, field0DIntInputPrev
    type(field1DInteger), pointer :: field1DIntInput, field1DIntOutput, field1DIntInputPrev
    type(field2DInteger), pointer :: field2DIntInput, field2DIntOutput, field2DIntInputPrev
    type(field3DInteger), pointer :: field3DIntInput, field3DIntOutput, field3DIntInputPrev

    ! field arrays
    type(field0DReal), dimension(:), pointer :: field0DRealArray
    type(field1DReal), dimension(:), pointer :: field1DRealArray
    type(field2DReal), dimension(:), pointer :: field2DRealArray
    type(field3DReal), dimension(:), pointer :: field3DRealArray
    type(field4DReal), dimension(:), pointer :: field4DRealArray
    type(field5DReal), dimension(:), pointer :: field5DRealArray

    type(field0DInteger), dimension(:), pointer :: field0DIntArray
    type(field1DInteger), dimension(:), pointer :: field1DIntArray
    type(field2DInteger), dimension(:), pointer :: field2DIntArray
    type(field3DInteger), dimension(:), pointer :: field3DIntArray

    character(len=strKIND) :: &
         fieldnameInput, & ! input field name
         poolnameInput     ! input pool name

    integer :: &
         iTime, & ! index of forcing times
         ierr ! error flag

    logical :: &
         lInputPoolCreated

    ! get the input field and pool name
    fieldnameInput = trim(forcingField % fieldname)//"_forcing_input"
    poolnameInput  = trim(forcingField % poolname)//"_forcing_input"

    ! remove the ouput field to the stream
    call MPAS_stream_mgr_set_property(&
         streamManager, trim(forcingStream % forcingStreamID), MPAS_STREAM_PROPERTY_IMMUTABLE, .false., ierr=ierr)
    if (ierr /= 0) then
       call mpas_log_write('Forcing: setup_input_fields: MPAS_stream_mgr_set_property: $i', MPAS_LOG_CRIT, intArgs=(/ierr/))
    endif
    call MPAS_stream_mgr_remove_field(&
         streamManager, trim(forcingStream % forcingStreamID), trim(forcingField % fieldname), ierr=ierr)
    if (ierr /= 0) then
       call mpas_log_write('Forcing: setup_input_fields: MPAS_stream_mgr_remove_field: $i', MPAS_LOG_CRIT, intArgs=(/ierr/))
    endif

    ! loop over blocks
    block => domain % blocklist
    do while (associated(block))

       ! get the output pool
       call MPAS_pool_get_subpool(block % structs, trim(forcingField % poolname), forcingPoolOutput)
       forcingPoolInfo % nDims = -1
       call MPAS_pool_get_field_info(forcingPoolOutput, trim(forcingField % fieldname), forcingPoolInfo)
       if (forcingPoolInfo % nDims == -1) then
          call mpas_log_write('Forcing: setup_input_fields: fieldname ' //trim(forcingField % fieldname) // &
               ' does not exist in poolname: ' // trim(forcingField % poolname), MPAS_LOG_CRIT)
       endif

       ! create the new input forcing pool if needed
       lInputPoolCreated = .false.

       nullify(forcingPoolInput)

       if (MPAS_pool_exists(block % structs, trim(poolnameInput), forcingPoolInput)) then
          call MPAS_pool_get_subpool(block % structs, trim(poolnameInput), forcingPoolInput)
       endif

       if (.not. associated(forcingPoolInput)) then
          call MPAS_pool_create_pool(forcingPoolInput)
          lInputPoolCreated = .true.
       endif

       if (forcingPoolInfo % fieldType == MPAS_POOL_REAL) then

          select case (forcingPoolInfo % nDims)
          case (0)

             ! get the output field
             call MPAS_pool_get_field(forcingPoolOutput, trim(forcingField % fieldname), field0DRealOutput, 1)
             ! allocate the time levels for the input field
             allocate(field0DRealArray(forcingStream % nTimeStencil))
             ! loop over the desired time levels
             do iTime = 1, forcingStream % nTimeStencil
                ! get temporary copy of output field
                call MPAS_duplicate_field(field0DRealOutput, field0DRealInput)
                ! set the input time level field to the temporary field
                field0DRealArray(iTime) = field0DRealInput
                ! nullify the next/prev pointers in the input time level field
                nullify(field0DRealArray(iTime) % next)
                nullify(field0DRealArray(iTime) % prev)
                ! if we have a previous block link the fields in the blocks
                if (associated(block % prev)) then
                   ! get the input field from the previous block
                   call mpas_pool_get_subpool(block % prev % structs, poolnameInput, forcingPoolInputPrev)
                   call mpas_pool_get_field(forcingPoolInputPrev, fieldnameInput, field0DRealInputPrev, iTime)
                   ! link the previous and current fields
                   field0DRealInputPrev % next => field0DRealArray(iTime)
                   field0DRealArray(iTime) % prev => field0DRealInputPrev
                endif
             enddo ! iTime
             ! add the input field to the pools
             call MPAS_pool_add_field(forcingPoolInput, fieldnameInput, field0DRealArray)
             call MPAS_pool_add_field(block % allFields, fieldnameInput, field0DRealArray)

          case (1)

             ! get the output field
             call MPAS_pool_get_field(forcingPoolOutput, trim(forcingField % fieldname), field1DRealOutput, 1)
             ! allocate the time levels for the input field
             allocate(field1DRealArray(forcingStream % nTimeStencil))
             ! loop over the desired time levels
             do iTime = 1, forcingStream % nTimeStencil
                ! get temporary copy of output field
                call MPAS_duplicate_field(field1DRealOutput, field1DRealInput)
                ! set the input time level field to the temporary field
                field1DRealArray(iTime) = field1DRealInput
                ! nullify the next/prev pointers in the input time level field
                nullify(field1DRealArray(iTime) % next)
                nullify(field1DRealArray(iTime) % prev)
                ! if we have a previous block link the fields in the blocks
                if (associated(block % prev)) then
                   ! get the input field from the previous block
                   call mpas_pool_get_subpool(block % prev % structs, poolnameInput, forcingPoolInputPrev)
                   call mpas_pool_get_field(forcingPoolInputPrev, fieldnameInput, field1DRealInputPrev, iTime)
                   ! link the previous and current fields
                   field1DRealInputPrev % next => field1DRealArray(iTime)
                   field1DRealArray(iTime) % prev => field1DRealInputPrev
                endif
             enddo ! iTime
             ! add the input field to the pools
             call MPAS_pool_add_field(forcingPoolInput, fieldnameInput, field1DRealArray)
             call MPAS_pool_add_field(block % allFields, fieldnameInput, field1DRealArray)

          case (2)

             ! get the output field
             call MPAS_pool_get_field(forcingPoolOutput, trim(forcingField % fieldname), field2DRealOutput, 1)
             ! allocate the time levels for the input field
             allocate(field2DRealArray(forcingStream % nTimeStencil))
             ! loop over the desired time levels
             do iTime = 1, forcingStream % nTimeStencil
                ! get temporary copy of output field
                call MPAS_duplicate_field(field2DRealOutput, field2DRealInput)
                ! set the input time level field to the temporary field
                field2DRealArray(iTime) = field2DRealInput
                ! nullify the next/prev pointers in the input time level field
                nullify(field2DRealArray(iTime) % next)
                nullify(field2DRealArray(iTime) % prev)
                ! if we have a previous block link the fields in the blocks
                if (associated(block % prev)) then
                   ! get the input field from the previous block
                   call mpas_pool_get_subpool(block % prev % structs, poolnameInput, forcingPoolInputPrev)
                   call mpas_pool_get_field(forcingPoolInputPrev, fieldnameInput, field2DRealInputPrev, iTime)
                   ! link the previous and current fields
                   field2DRealInputPrev % next => field2DRealArray(iTime)
                   field2DRealArray(iTime) % prev => field2DRealInputPrev
                endif
             enddo ! iTime
             ! add the input field to the pools
             call MPAS_pool_add_field(forcingPoolInput, fieldnameInput, field2DRealArray)
             call MPAS_pool_add_field(block % allFields, fieldnameInput, field2DRealArray)

          case (3)

             ! get the output field
             call MPAS_pool_get_field(forcingPoolOutput, trim(forcingField % fieldname), field3DRealOutput, 1)
             ! allocate the time levels for the input field
             allocate(field3DRealArray(forcingStream % nTimeStencil))
             ! loop over the desired time levels
             do iTime = 1, forcingStream % nTimeStencil
                ! get temporary copy of output field
                call MPAS_duplicate_field(field3DRealOutput, field3DRealInput)
                ! set the input time level field to the temporary field
                field3DRealArray(iTime) = field3DRealInput
                ! nullify the next/prev pointers in the input time level field
                nullify(field3DRealArray(iTime) % next)
                nullify(field3DRealArray(iTime) % prev)
                ! if we have a previous block link the fields in the blocks
                if (associated(block % prev)) then
                   ! get the input field from the previous block
                   call mpas_pool_get_subpool(block % prev % structs, poolnameInput, forcingPoolInputPrev)
                   call mpas_pool_get_field(forcingPoolInputPrev, fieldnameInput, field3DRealInputPrev, iTime)
                   ! link the previous and current fields
                   field3DRealInputPrev % next => field3DRealArray(iTime)
                   field3DRealArray(iTime) % prev => field3DRealInputPrev
                endif
             enddo ! iTime
             ! add the input field to the pools
             call MPAS_pool_add_field(forcingPoolInput, fieldnameInput, field3DRealArray)
             call MPAS_pool_add_field(block % allFields, fieldnameInput, field3DRealArray)

          case (4)

             ! get the output field
             call MPAS_pool_get_field(forcingPoolOutput, trim(forcingField % fieldname), field4DRealOutput, 1)
             ! allocate the time levels for the input field
             allocate(field4DRealArray(forcingStream % nTimeStencil))
             ! loop over the desired time levels
             do iTime = 1, forcingStream % nTimeStencil
                ! get temporary copy of output field
                call MPAS_duplicate_field(field4DRealOutput, field4DRealInput)
                ! set the input time level field to the temporary field
                field4DRealArray(iTime) = field4DRealInput
                ! nullify the next/prev pointers in the input time level field
                nullify(field4DRealArray(iTime) % next)
                nullify(field4DRealArray(iTime) % prev)
                ! if we have a previous block link the fields in the blocks
                if (associated(block % prev)) then
                   ! get the input field from the previous block
                   call mpas_pool_get_subpool(block % prev % structs, poolnameInput, forcingPoolInputPrev)
                   call mpas_pool_get_field(forcingPoolInputPrev, fieldnameInput, field4DRealInputPrev, iTime)
                   ! link the previous and current fields
                   field4DRealInputPrev % next => field4DRealArray(iTime)
                   field4DRealArray(iTime) % prev => field4DRealInputPrev
                endif
             enddo ! iTime
             ! add the input field to the pools
             call MPAS_pool_add_field(forcingPoolInput, fieldnameInput, field4DRealArray)
             call MPAS_pool_add_field(block % allFields, fieldnameInput, field4DRealArray)

          case (5)

             ! get the output field
             call MPAS_pool_get_field(forcingPoolOutput, trim(forcingField % fieldname), field5DRealOutput, 1)
             ! allocate the time levels for the input field
             allocate(field5DRealArray(forcingStream % nTimeStencil))
             ! loop over the desired time levels
             do iTime = 1, forcingStream % nTimeStencil
                ! get temporary copy of output field
                call MPAS_duplicate_field(field5DRealOutput, field5DRealInput)
                ! set the input time level field to the temporary field
                field5DRealArray(iTime) = field5DRealInput
                ! nullify the next/prev pointers in the input time level field
                nullify(field5DRealArray(iTime) % next)
                nullify(field5DRealArray(iTime) % prev)
                ! if we have a previous block link the fields in the blocks
                if (associated(block % prev)) then
                   ! get the input field from the previous block
                   call mpas_pool_get_subpool(block % prev % structs, poolnameInput, forcingPoolInputPrev)
                   call mpas_pool_get_field(forcingPoolInputPrev, fieldnameInput, field5DRealInputPrev, iTime)
                   ! link the previous and current fields
                   field5DRealInputPrev % next => field5DRealArray(iTime)
                   field5DRealArray(iTime) % prev => field5DRealInputPrev
                endif
             enddo ! iTime
             ! add the input field to the pools
             call MPAS_pool_add_field(forcingPoolInput, fieldnameInput, field5DRealArray)
             call MPAS_pool_add_field(block % allFields, fieldnameInput, field5DRealArray)

          end select

       else if (forcingPoolInfo % fieldType == MPAS_POOL_INTEGER) then

          select case (forcingPoolInfo % nDims)
          case (0)

             ! get the output field
             call MPAS_pool_get_field(forcingPoolOutput, trim(forcingField % fieldname), field0DIntOutput, 1)
             ! allocate the time levels for the input field
             allocate(field0DIntArray(forcingStream % nTimeStencil))
             ! loop over the desired time levels
             do iTime = 1, forcingStream % nTimeStencil
                ! get temporary copy of output field
                call MPAS_duplicate_field(field0DIntOutput, field0DIntInput)
                ! set the input time level field to the temporary field
                field0DIntArray(iTime) = field0DIntInput
                ! nullify the next/prev pointers in the input time level field
                nullify(field0DIntArray(iTime) % next)
                nullify(field0DIntArray(iTime) % prev)
                ! if we have a previous block link the fields in the blocks
                if (associated(block % prev)) then
                   ! get the input field from the previous block
                   call mpas_pool_get_subpool(block % prev % structs, poolnameInput, forcingPoolInputPrev)
                   call mpas_pool_get_field(forcingPoolInputPrev, fieldnameInput, field0DIntInputPrev, iTime)
                   ! link the previous and current fields
                   field0DIntInputPrev % next => field0DIntArray(iTime)
                   field0DIntArray(iTime) % prev => field0DIntInputPrev
                endif
             enddo ! iTime
             ! add the input field to the pools
             call MPAS_pool_add_field(forcingPoolInput, fieldnameInput, field0DIntArray)
             call MPAS_pool_add_field(block % allFields, fieldnameInput, field0DIntArray)

          case (1)

             ! get the output field
             call MPAS_pool_get_field(forcingPoolOutput, trim(forcingField % fieldname), field1DIntOutput, 1)
             ! allocate the time levels for the input field
             allocate(field1DIntArray(forcingStream % nTimeStencil))
             ! loop over the desired time levels
             do iTime = 1, forcingStream % nTimeStencil
                ! get temporary copy of output field
                call MPAS_duplicate_field(field1DIntOutput, field1DIntInput)
                ! set the input time level field to the temporary field
                field1DIntArray(iTime) = field1DIntInput
                ! nullify the next/prev pointers in the input time level field
                nullify(field1DIntArray(iTime) % next)
                nullify(field1DIntArray(iTime) % prev)
                ! if we have a previous block link the fields in the blocks
                if (associated(block % prev)) then
                   ! get the input field from the previous block
                   call mpas_pool_get_subpool(block % prev % structs, poolnameInput, forcingPoolInputPrev)
                   call mpas_pool_get_field(forcingPoolInputPrev, fieldnameInput, field1DIntInputPrev, iTime)
                   ! link the previous and current fields
                   field1DIntInputPrev % next => field1DIntArray(iTime)
                   field1DIntArray(iTime) % prev => field1DIntInputPrev
                endif
             enddo ! iTime
             ! add the input field to the pools
             call MPAS_pool_add_field(forcingPoolInput, fieldnameInput, field1DIntArray)
             call MPAS_pool_add_field(block % allFields, fieldnameInput, field1DIntArray)

          case (2)

             ! get the output field
             call MPAS_pool_get_field(forcingPoolOutput, trim(forcingField % fieldname), field2DIntOutput, 1)
             ! allocate the time levels for the input field
             allocate(field2DIntArray(forcingStream % nTimeStencil))
             ! loop over the desired time levels
             do iTime = 1, forcingStream % nTimeStencil
                ! get temporary copy of output field
                call MPAS_duplicate_field(field2DIntOutput, field2DIntInput)
                ! set the input time level field to the temporary field
                field2DIntArray(iTime) = field2DIntInput
                ! nullify the next/prev pointers in the input time level field
                nullify(field2DIntArray(iTime) % next)
                nullify(field2DIntArray(iTime) % prev)
                ! if we have a previous block link the fields in the blocks
                if (associated(block % prev)) then
                   ! get the input field from the previous block
                   call mpas_pool_get_subpool(block % prev % structs, poolnameInput, forcingPoolInputPrev)
                   call mpas_pool_get_field(forcingPoolInputPrev, fieldnameInput, field2DIntInputPrev, iTime)
                   ! link the previous and current fields
                   field2DIntInputPrev % next => field2DIntArray(iTime)
                   field2DIntArray(iTime) % prev => field2DIntInputPrev
                endif
             enddo ! iTime
             ! add the input field to the pools
             call MPAS_pool_add_field(forcingPoolInput, fieldnameInput, field2DIntArray)
             call MPAS_pool_add_field(block % allFields, fieldnameInput, field2DIntArray)

          case (3)

             ! get the output field
             call MPAS_pool_get_field(forcingPoolOutput, trim(forcingField % fieldname), field3DIntOutput, 1)
             ! allocate the time levels for the input field
             allocate(field3DIntArray(forcingStream % nTimeStencil))
             ! loop over the desired time levels
             do iTime = 1, forcingStream % nTimeStencil
                ! get temporary copy of output field
                call MPAS_duplicate_field(field3DIntOutput, field3DIntInput)
                ! set the input time level field to the temporary field
                field3DIntArray(iTime) = field3DIntInput
                ! nullify the next/prev pointers in the input time level field
                nullify(field3DIntArray(iTime) % next)
                nullify(field3DIntArray(iTime) % prev)
                ! if we have a previous block link the fields in the blocks
                if (associated(block % prev)) then
                   ! get the input field from the previous block
                   call mpas_pool_get_subpool(block % prev % structs, poolnameInput, forcingPoolInputPrev)
                   call mpas_pool_get_field(forcingPoolInputPrev, fieldnameInput, field3DIntInputPrev, iTime)
                   ! link the previous and current fields
                   field3DIntInputPrev % next => field3DIntArray(iTime)
                   field3DIntArray(iTime) % prev => field3DIntInputPrev
                endif
             enddo ! iTime
             ! add the input field to the pools
             call MPAS_pool_add_field(forcingPoolInput, fieldnameInput, field3DIntArray)
             call MPAS_pool_add_field(block % allFields, fieldnameInput, field3DIntArray)

             !call MPAS_pool_links_pools

          end select

       else

          ! unsupported forcing type
          call mpas_log_write('Forcing: Forcing only supports REAL and INTEGER data types. Field type=$i', MPAS_LOG_CRIT, &
               intArgs=(/forcingPoolInfo % fieldType/))

       endif

       ! add the new input pool to the block
       if (lInputPoolCreated) then
          call MPAS_pool_add_subpool(block % structs, trim(poolnameInput), forcingPoolInput)
       endif

       block => block % next
    enddo ! blocklist

    ! add the input field to the stream
    call MPAS_stream_mgr_add_field(&
         streamManager, trim(forcingStream % forcingStreamID), trim(fieldnameInput), ierr=ierr)
    if (ierr /= 0) then
       call mpas_log_write('Forcing: setup_input_fields: MPAS_stream_mgr_add_field: $i', MPAS_LOG_CRIT, intArgs=(/ierr/))
    endif
    call MPAS_stream_mgr_set_property(&
         streamManager, trim(forcingStream % forcingStreamID), MPAS_STREAM_PROPERTY_IMMUTABLE, .true., ierr=ierr)
    if (ierr /= 0) then
       call mpas_log_write('Forcing: setup_input_fields: MPAS_stream_mgr_set_property: $i', MPAS_LOG_CRIT, intArgs=(/ierr/))
    endif

  end subroutine setup_input_fields!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  forcing_data_interpolation
!
!> \brief Perform forcing interpolation
!> \author Adrian K. Turner, LANL
!> \date September 25th 2014
!> \details
!>  Given two input times and time level data perform interpolation
!>  to produce output data
!
!-----------------------------------------------------------------------

  subroutine forcing_data_interpolation(&!{{{
       domain, &
       forcingStream, &
       currentTime)

    type (domain_type), pointer :: &
         domain ! domain object

    type (mpas_forcing_stream_type), intent(in) :: &
         forcingStream !< Input: forcing object

    type (MPAS_Time_Type), intent(in) :: &
         currentTime !< Input: the current forcing time

    type (MPAS_pool_type), pointer :: &
         forcingPoolInput, & ! the input pool
         forcingPoolOutput   ! the output pool

    type (MPAS_pool_field_info_type) :: &
         forcingPoolInfo ! the pool info object

    real(kind=RKIND), dimension(:), allocatable :: &
         interpolants ! the interpolation weights

    ! data arrays
    real(kind=RKIND),                       pointer :: field0DRealInput, field0DRealOutput
    real(kind=RKIND), dimension(:),         pointer :: field1DRealInput, field1DRealOutput
    real(kind=RKIND), dimension(:,:),       pointer :: field2DRealInput, field2DRealOutput
    real(kind=RKIND), dimension(:,:,:),     pointer :: field3DRealInput, field3DRealOutput
    real(kind=RKIND), dimension(:,:,:,:),   pointer :: field4DRealInput, field4DRealOutput
    real(kind=RKIND), dimension(:,:,:,:,:), pointer :: field5DRealInput, field5DRealOutput

    integer,                   pointer :: field0DIntInput, field0DIntOutput
    integer, dimension(:),     pointer :: field1DIntInput, field1DIntOutput
    integer, dimension(:,:),   pointer :: field2DIntInput, field2DIntOutput
    integer, dimension(:,:,:), pointer :: field3DIntInput, field3DIntOutput

    integer :: &
         iTime ! index of forcing times

    character(len=strKIND) :: &
         poolnameInput, & ! input pool name
         fieldnameInput   ! input pool name

    type(block_type), pointer :: &
         block ! block type

    type(mpas_forcing_field_type), pointer :: &
         forcingField

    ! get the interpolant weights
    allocate(interpolants(forcingStream % nTimeStencil))
    call get_interpolants(interpolants, forcingStream, currentTime)

    ! loop over forcing fields in stream
    forcingField => forcingStream % field
    do while (associated(forcingField))

       ! get input names
       poolnameInput  = trim(forcingField % poolname)//"_forcing_input"
       fieldnameInput = trim(forcingField % fieldname)//"_forcing_input"

       FORCING_DEBUG_WRITE('-- Forcing: forcing_data_interpolation pool: '//trim(poolnameInput)//" field: "//trim(fieldnameInput))

       ! loop over blocks
       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, trim(forcingField % poolname), forcingPoolOutput)
          call MPAS_pool_get_field_info(forcingPoolOutput, trim(forcingField % fieldname), forcingPoolInfo)

          call MPAS_pool_get_subpool(block % structs, trim(poolnameInput), forcingPoolInput)

          if (forcingPoolInfo % fieldType == MPAS_POOL_REAL) then

             select case (forcingPoolInfo % nDims)
             case (0)

                call MPAS_pool_get_array(forcingPoolOutput, trim(forcingField % fieldname), field0DRealOutput, 1)
                call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field0DRealInput, 1)
                field0DRealOutput = field0DRealInput * interpolants(1)
                do iTime = 2, forcingStream % nTimeStencil
                   call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field0DRealInput, iTime)
                   field0DRealOutput = field0DRealOutput + field0DRealInput * interpolants(iTime)
                enddo ! iTime

             case (1)

                call MPAS_pool_get_array(forcingPoolOutput, trim(forcingField % fieldname), field1DRealOutput, 1)
                call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field1DRealInput, 1)
                field1DRealOutput = field1DRealInput * interpolants(1)
                do iTime = 2, forcingStream % nTimeStencil
                   call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field1DRealInput, iTime)
                   field1DRealOutput = field1DRealOutput + field1DRealInput * interpolants(iTime)
                enddo ! iTime

             case (2)

                call MPAS_pool_get_array(forcingPoolOutput, trim(forcingField % fieldname), field2DRealOutput, 1)
                call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field2DRealInput, 1)
                field2DRealOutput = field2DRealInput * interpolants(1)
                do iTime = 2, forcingStream % nTimeStencil
                   call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field2DRealInput, iTime)
                   field2DRealOutput = field2DRealOutput + field2DRealInput * interpolants(iTime)
                enddo ! iTime

             case (3)

                call MPAS_pool_get_array(forcingPoolOutput, trim(forcingField % fieldname), field3DRealOutput, 1)
                call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field3DRealInput, 1)
                field3DRealOutput = field3DRealInput * interpolants(1)
                do iTime = 2, forcingStream % nTimeStencil
                   call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field3DRealInput, iTime)
                   field3DRealOutput = field3DRealOutput + field3DRealInput * interpolants(iTime)
                enddo ! iTime

             case (4)

                call MPAS_pool_get_array(forcingPoolOutput, trim(forcingField % fieldname), field4DRealOutput, 1)
                call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field4DRealInput, 1)
                field4DRealOutput = field4DRealInput * interpolants(1)
                do iTime = 2, forcingStream % nTimeStencil
                   call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field4DRealInput, iTime)
                   field4DRealOutput = field4DRealOutput + field4DRealInput * interpolants(iTime)
                enddo ! iTime

             case (5)

                call MPAS_pool_get_array(forcingPoolOutput, trim(forcingField % fieldname), field5DRealOutput, 1)
                call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field5DRealInput, 1)
                field5DRealOutput = field5DRealInput * interpolants(1)
                do iTime = 2, forcingStream % nTimeStencil
                   call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field5DRealInput, iTime)
                   field5DRealOutput = field5DRealOutput + field5DRealInput * interpolants(iTime)
                enddo ! iTime

             end select

          else if (forcingPoolInfo % fieldType == MPAS_POOL_INTEGER) then

             select case (forcingPoolInfo % nDims)
             case (0)

                call MPAS_pool_get_array(forcingPoolOutput, trim(forcingField % fieldname), field0DIntOutput, 1)
                call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field0DIntInput, 1)
                field0DIntOutput = field0DIntInput * interpolants(1)
                do iTime = 2, forcingStream % nTimeStencil
                   call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field0DIntInput, iTime)
                   field0DIntOutput = field0DIntOutput + field0DIntInput * interpolants(iTime)
                enddo ! iTime

             case (1)

                call MPAS_pool_get_array(forcingPoolOutput, trim(forcingField % fieldname), field1DIntOutput, 1)
                call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field1DIntInput, 1)
                field1DIntOutput = field1DIntInput * interpolants(1)
                do iTime = 2, forcingStream % nTimeStencil
                   call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field1DIntInput, iTime)
                   field1DIntOutput = field1DIntOutput + field1DIntInput * interpolants(iTime)
                enddo ! iTime

             case (2)

                call MPAS_pool_get_array(forcingPoolOutput, trim(forcingField % fieldname), field2DIntOutput, 1)
                call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field2DIntInput, 1)
                field2DIntOutput = field2DIntInput * interpolants(1)
                do iTime = 2, forcingStream % nTimeStencil
                   call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field2DIntInput, iTime)
                   field2DIntOutput = field2DIntOutput + field2DIntInput * interpolants(iTime)
                enddo ! iTime

             case (3)

                call MPAS_pool_get_array(forcingPoolOutput, trim(forcingField % fieldname), field3DIntOutput, 1)
                call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field3DIntInput, 1)
                field3DIntOutput = field3DIntInput * interpolants(1)
                do iTime = 2, forcingStream % nTimeStencil
                   call MPAS_pool_get_array(forcingPoolInput,  trim(fieldnameInput), field3DIntInput, iTime)
                   field3DIntOutput = field3DIntOutput + field3DIntInput * interpolants(iTime)
                enddo ! iTime

             end select

          endif

          block => block % next
       end do

       forcingField => forcingField % next
    end do

    deallocate(interpolants)

  end subroutine forcing_data_interpolation!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  forcing_shift_data
!
!> \brief shift the data in the fields down one time level
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  the input field pointers are shifted so there is room in the last
!>  time level for new data
!
!-----------------------------------------------------------------------

  subroutine forcing_shift_data(&!{{{
       domain, &
       forcingStream)

    type(domain_type), pointer :: &
         domain ! domain type

    type(mpas_forcing_stream_type) :: &
         forcingStream ! forcing object

    type(block_type), pointer :: &
         block ! block type

    type(MPAS_pool_field_info_type) :: &
         forcingPoolInfo ! pool info type

    type(MPAS_pool_type), pointer :: &
         forcingPoolInput ! input pool

    ! input field pointers
    type(field0DReal), pointer :: field0DRealInput1, field0DRealInput2
    type(field1DReal), pointer :: field1DRealInput1, field1DRealInput2
    type(field2DReal), pointer :: field2DRealInput1, field2DRealInput2
    type(field3DReal), pointer :: field3DRealInput1, field3DRealInput2
    type(field4DReal), pointer :: field4DRealInput1, field4DRealInput2
    type(field5DReal), pointer :: field5DRealInput1, field5DRealInput2

    type(field0DInteger), pointer :: field0DIntInput1, field0DIntInput2
    type(field1DInteger), pointer :: field1DIntInput1, field1DIntInput2
    type(field2DInteger), pointer :: field2DIntInput1, field2DIntInput2
    type(field3DInteger), pointer :: field3DIntInput1, field3DIntInput2

    integer :: &
         iTime ! index of forcing times

    character(len=strKIND) :: &
         poolnameInput, & ! input pool name
         fieldnameInput   ! input pool name

    type(MPAS_pool_data_type), pointer :: &
         forcingFieldInput

    type(mpas_forcing_field_type), pointer :: &
         forcingField

    ! loop over forcing fields in stream
    forcingField => forcingStream % field
    do while (associated(forcingField))

       ! get input names
       poolnameInput  = trim(forcingField % poolname)//"_forcing_input"
       fieldnameInput = trim(forcingField % fieldname)//"_forcing_input"

       FORCING_DEBUG_WRITE('-- Forcing: forcing_shift_data pool: '//trim(poolnameInput)//" field: "//trim(fieldnameInput))
       
       ! No need to loop over blocks since mpas_shift_time_levs does this

       block => domain % blocklist

       call MPAS_pool_get_subpool(block % structs, trim(poolnameInput), forcingPoolInput)

       forcingFieldInput => pool_get_member(forcingPoolInput, trim(fieldnameInput), MPAS_POOL_FIELD)

       if (forcingFieldInput % contentsTimeLevs < 2) then
          FORCING_ERROR_WRITE('-- Forcing: forcing_shift_data: too few timelevels Pool: '//trim(poolnameInput)//' Field: '//trim(fieldnameInput))
          call mpas_log_write('Forcing: forcing_shift_data: too few timelevels in Pool: '//trim(poolnameInput)// &
               ' Field: '//trim(fieldnameInput), MPAS_LOG_CRIT)
       endif

       if (forcingFieldInput % contentsType == MPAS_POOL_REAL) then

          select case (forcingFieldInput % contentsDims)
          case (0)
             call mpas_shift_time_levs(forcingFieldInput % r0a)
          case (1)
             call mpas_shift_time_levs(forcingFieldInput % r1a)
          case (2)
             call mpas_shift_time_levs(forcingFieldInput % r2a)
          case (3)
             call mpas_shift_time_levs(forcingFieldInput % r3a)
          case (4)
             call mpas_shift_time_levs(forcingFieldInput % r4a)
          case (5)
             call mpas_shift_time_levs(forcingFieldInput % r5a)
          end select

       else if (forcingFieldInput % contentsType == MPAS_POOL_INTEGER) then

          select case (forcingFieldInput % contentsDims)
          case (0)
             call mpas_shift_time_levs(forcingFieldInput % i0a)
          case (1)
             call mpas_shift_time_levs(forcingFieldInput % i1a)
          case (2)
             call mpas_shift_time_levs(forcingFieldInput % i2a)
          case (3)
             call mpas_shift_time_levs(forcingFieldInput % i3a)
          end select

          endif

       forcingField => forcingField % next
    end do

  end subroutine forcing_shift_data!}}}

!-----------------------------------------------------------------------
! restarts
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_forcing_write_restart_times
!
!> \brief set forcing restart times registry variables
!> \author Adrian K. Turner, LANL
!> \date 9th December 2014
!> \details
!>  loop over the forcing groups in the forcing group object and write
!>  out the forcing clock times to registry variables that are included
!>  in the restart stream. 'forcingGroupHead' is the forcing group
!>  object
!
!-----------------------------------------------------------------------

  subroutine mpas_forcing_write_restart_times(forcingGroupHead)!{{{

    type(mpas_forcing_group_type), pointer :: &
         forcingGroupHead ! forcing group linked list head pointer

    type(mpas_forcing_group_type), pointer :: &
         forcingGroup ! forcing group iterator

    type(MPAS_time_type) :: &
         forcingClockTime

    character(len=strKIND) :: &
         forcingClockTimeStr

    type(block_type), pointer :: &
         block

    type(mpas_pool_type), pointer :: &
         forcingPool

    integer, pointer :: &
         nForcingGroupCounter

    character(len=strKIND), dimension(:), pointer :: &
         forcingGroupNames, &
         forcingGroupRestartTimes

    ! loop over forcing groups
    forcingGroup => forcingGroupHead
    do while (associated(forcingGroup))

       block => forcingGroup % domain_ptr % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "forcing", forcingPool)

          call MPAS_pool_get_array(forcingPool, "nForcingGroupCounter", nforcingGroupCounter)

          nforcingGroupCounter = 0

          block => block % next
       end do

       forcingGroup => forcingGroup % next
    end do

    ! loop over forcing groups
    forcingGroup => forcingGroupHead
    do while (associated(forcingGroup))

       ! get the forcing clock time
       forcingClockTime = MPAS_get_clock_time(forcingGroup % forcingClock, MPAS_NOW)

       call MPAS_get_time(forcingClockTime, dateTimeString=forcingClockTimeStr)

       block => forcingGroup % domain_ptr % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "forcing", forcingPool)

          call MPAS_pool_get_array(forcingPool, "nForcingGroupCounter",     nForcingGroupCounter)
          call MPAS_pool_get_array(forcingPool, "forcingGroupNames",        forcingGroupNames)
          call MPAS_pool_get_array(forcingPool, "forcingGroupRestartTimes", forcingGroupRestartTimes)

          nForcingGroupCounter = nForcingGroupCounter + 1
          forcingGroupNames       (nforcingGroupCounter) = trim(forcingGroup % forcingGroupName)
          forcingGroupRestartTimes(nforcingGroupCounter) = trim(forcingClockTimeStr)

          block => block % next
       end do

       forcingGroup => forcingGroup % next
    end do

  end subroutine mpas_forcing_write_restart_times!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  add_forcing_group_to_restart_counter
!
!> \brief checks number of forcing groups
!> \author Adrian K. Turner, LANL
!> \date 16th May 2017
!> \details
!>  This routine checks that the number of forcing groups per domain
!>  does not exceed the maximum allowed.
!
!-----------------------------------------------------------------------

  subroutine add_forcing_group_to_restart_counter(forcingGroup)

    type(mpas_forcing_group_type), pointer :: &
         forcingGroup ! forcing group

    type(block_type), pointer :: &
         block

    integer, pointer :: &
         nForcingGroupsMax, &
         nforcingGroupCounter

    type(mpas_pool_type), pointer :: &
         forcingPool

    block => forcingGroup % domain_ptr % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nForcingGroupsMax", nForcingGroupsMax)

       call MPAS_pool_get_subpool(block % structs, "forcing", forcingPool)

       call MPAS_pool_get_array(forcingPool, "nForcingGroupCounter", nforcingGroupCounter)

       nforcingGroupCounter = nforcingGroupCounter + 1

       if (nforcingGroupCounter > nForcingGroupsMax) then

          call mpas_log_write(&
               "Init forcing group: "//trim(forcingGroup % forcingGroupName)//&
               ", nForcingGroupsMax too small, nForcingGroupsMax: $i, nforcingGroupCounter: $i", &
               MPAS_LOG_CRIT, &
               intArgs=(/nForcingGroupsMax,nforcingGroupCounter/))

       endif

       block => block % next
    end do

  end subroutine add_forcing_group_to_restart_counter

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_restart_times
!
!> \brief get forcing restart times and create forcing clocks.
!> \author Adrian K. Turner, LANL
!> \date 9th December 2014
!> \details
!>  The forcing clock restart times have been read in from the restart
!>  file. This routine sets the forcing group clock from those times.
!
!-----------------------------------------------------------------------

  subroutine get_restart_times(&!{{{
       forcingGroup, &
       timeStep, &
       stopTime)

    type(mpas_forcing_group_type), pointer :: &
         forcingGroup ! forcing group to restart

    type(MPAS_TimeInterval_type), intent(in) :: &
         timeStep ! simulation time step

    type(MPAS_Time_type), intent(in) :: &
         stopTime     ! stop time of forcing clock - !!!! SHOULDNT BE NEEDED

    type(MPAS_time_type) :: &
         forcingClockTime

    type(mpas_pool_type), pointer :: &
         forcingPool

    integer, pointer :: &
         nForcingGroupsMax

    integer :: &
         iForcingGroup

    character(len=strKIND), dimension(:), pointer :: &
         forcingGroupNames, &
         forcingGroupRestartTimes

    call mpas_pool_get_subpool(forcingGroup % domain_ptr % blocklist % structs, "forcing", forcingPool)

    call MPAS_pool_get_array(forcingPool, "forcingGroupNames",        forcingGroupNames)
    call MPAS_pool_get_array(forcingPool, "forcingGroupRestartTimes", forcingGroupRestartTimes)

    call mpas_pool_get_dimension(forcingGroup % domain_ptr % blocklist % dimensions, "nForcingGroupsMax", nForcingGroupsMax)

    do iForcingGroup = 1, nForcingGroupsMax

       if (trim(forcingGroup % forcingGroupName) == trim(forcingGroupNames(iForcingGroup)(1:ShortStrKind))) then

          ! set the forcing group time
          call MPAS_set_time(forcingClockTime, dateTimeString=trim(forcingGroupRestartTimes(iForcingGroup)(1:ShortStrKind)))
          call mpas_create_clock(forcingGroup % forcingClock, startTime=forcingClockTime, timeStep=timeStep, stopTime=stopTime)

       endif

    enddo ! iForcingGroup

  end subroutine get_restart_times!}}}

!-----------------------------------------------------------------------

end module mpas_forcing
