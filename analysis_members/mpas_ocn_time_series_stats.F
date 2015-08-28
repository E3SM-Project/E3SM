! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! ocn_time_series_stats
!
!> \brief MPAS ocean analysis core member: time_series_stats
!> \author Jon Woodring
!> \date   September 1, 2015
!> \details
!>  Flexible time series averaging, mins, and maxes of fields.
!-----------------------------------------------------------------------
module ocn_time_series_stats
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_dmpar
  use mpas_timekeeping
  use mpas_stream_manager

  use ocn_constants
  use ocn_diagnostics_routines

  implicit none
  private
  save

  ! Public parameters
  !--------------------------------------------------------------------

  ! Public member functions
  !--------------------------------------------------------------------
  public :: ocn_init_time_series_stats, &
         ocn_compute_time_series_stats, &
         ocn_restart_time_series_stats, &
         ocn_finalize_time_series_stats

  ! Private module variables
  !--------------------------------------------------------------------

  ! startup, interval, and restart is done in the outer analysis driver

  ! time buffer type
  ! this keeps track of timers and if and when they need to accumulate
  type time_buffer_type
    ! internal state
    logical :: started_flag, accumulate_flag, reset_flag
    real (kind=RKIND) :: total_accum

    type (MPAS_Time_type) :: start_time
    type (MPAS_TimeInterval_type) :: duration_interval
    type (MPAS_TimeInterval_type) :: repeat_interval
    type (MPAS_TimeInterval_type) :: reset_interval

    ! alarm IDs
    character (len=StrKIND) :: start_alarm_ID
    character (len=StrKIND) :: repeat_alarm_ID
    character (len=StrKIND) :: duration_alarm_ID
    character (len=StrKIND) :: reset_alarm_ID

    ! name of counter
    character (len=StrKIND) :: output_counter
  end type time_buffer_type

  ! time variable type
  ! this keeps track of arrays, array types, and names
  type time_variable_type
    type (mpas_pool_field_info_type) :: info
    character (len=StrKIND) :: input_name
    ! either you have to put a number of buffers per variable
    ! or put the output variables in the buffers (I decided to put it here)
    character (len=StrKIND), dimension(:), allocatable :: output_names
  end type time_variable_type

  ! operation
  integer :: operation
  
  ! stream name
  character (len=StrKIND), pointer :: stream_name

  ! information per variable
  type (time_variable_type), dimension(:), allocatable :: variables

  ! information per buffer
  type (time_buffer_type), dimension(:), allocatable :: buffers

  ! enum of ops and types
  integer, parameter :: AVG_OP = 1
  integer, parameter :: MIN_OP = 2
  integer, parameter :: MAX_OP = 3

  integer, parameter :: START_TIMES = 5
  integer, parameter :: DURATION_INTERVALS = 6
  integer, parameter :: REPEAT_INTERVALS = 7
  integer, parameter :: RESET_INTERVALS = 8

!***********************************************************************
contains



!***********************************************************************
! routine ocn_init_time_series_stats
!
!> \brief Initialize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_init_time_series_stats(domain, err)!{{{
  ! input variables

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: b
  integer :: number_of_variables, number_of_buffers
  character (len=StrKIND) :: instance ! TODO intent(in)
  character (len=StrKIND) :: prefix, op
  character (len=StrKIND), pointer :: stream_name
 
  ! start procedure
  err = 0

  ! string representation
  ! TODO placeholder for some unique ID if this code is replicated
  !    per multiple AMs for multiple streams
  instance = ''
  prefix = 'config_AM_timeSeriesStats' // trim(instance)

  ! get the basic configuration of this stream
  call start_init(domain, prefix, number_of_variables, number_of_buffers, &
    stream_name, op, err)

  ! modify the stream to remove existing vars and add accumulated versions 
  call modify_stream(domain, stream_name, number_of_variables, &
    number_of_buffers, instance, prefix, op, err)

  ! get all of the timing and configuration
  call get_alarms(domain, prefix, number_of_buffers, err)

  ! set all of the alarms based on timers
  call set_alarms(domain % clock, instance, number_of_buffers, err)
end subroutine ocn_init_time_series_stats!}}}



!***********************************************************************
! routine ocn_compute_time_series_stats
!
!> \brief Compute MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_compute_time_series_stats(domain, timeLevel, err)!{{{
  ! input variables
  integer, intent(in) :: timeLevel

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: i, v, b
  real (kind=RKIND), pointer :: counter

  ! start procedure
  err = 0

  ! update the counter
  do b = 1, size(buffers)
    if (buffers(b) % accumulate_flag) then
      if (buffers(b) % reset_flag) then
        buffers(b) % total_accum = 1
      else 
        buffers(b) % total_accum = buffers(b) % total_accum + 1
      end if

      ! update the stream
      call mpas_pool_get_array(domain % blocklist % allFields, &
        buffers(b) % output_counter, counter, 1)
      counter = buffers(b) % total_accum
    end if
  end do

  ! do all of the operations
  do v = 1, size(variables)
    call typed_operate(domain % blocklist, v, operation)
  end do

  ! do all of the time checking and flag setting
  call timer_checking(domain % clock, err)

end subroutine ocn_compute_time_series_stats!}}}



!***********************************************************************
! routine ocn_restart_time_series_stats
!
!> \brief Save restart for MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_restart_time_series_stats(domain, err)!{{{
  ! input variables

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables

  ! start procedure
  err = 0

  ! TODO is there anything needed here?
end subroutine ocn_restart_time_series_stats!}}}



!***********************************************************************
! routine ocn_finalize_time_series_stats
!
!> \brief Finalize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_finalize_time_series_stats(domain, err)!{{{
  ! input variables

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: i, v

  ! start procedure
  err = 0

  ! clean up memory
  if (allocated(buffers)) then
    deallocate(buffers)
  end if
  if (allocated(variables)) then
    do v = 1, size(variables)
      if (allocated(variables(v) % output_names)) &
        then
          deallocate(variables(v) % output_names)
      end if
    end do
    deallocate(variables)
  end if

end subroutine ocn_finalize_time_series_stats!}}}

!
! local subroutines
!

!***********************************************************************
! routine start_init
!
!> \brief Begin the initialization of this analysis member
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!> This will count the number of variables, number of buffers, and
!> also get the stream name and operation strings.
!-----------------------------------------------------------------------
subroutine start_init(domain, prefix, number_of_variables, &
  number_of_buffers, stream_name, op, err)
  ! input variables
  character (len=StrKIND), intent(in) :: prefix

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  character (len=StrKIND), pointer, intent(out) :: stream_name
  character (len=StrKIND), intent(out) :: op
  integer, intent(out) :: number_of_variables, number_of_buffers
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  character (len=StrKIND), pointer :: config_results
  character (len=StrKIND) :: copy, config
  integer :: b
 
  ! start procedure
  err = 0

  ! get the stream name
  config = trim(prefix) // '_stream_name'
  call mpas_pool_get_config(domain % configs, config, stream_name)

  if (stream_name == 'none') then
    call mpas_dmpar_global_abort('Error: stream cannot be "none" ' // &
      'for time series stats.')
  end if

  ! count the number of variables
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    stream_name, err)
  number_of_variables = 0
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    stream_name, copy))
    number_of_variables = number_of_variables + 1
  end do

  ! count the number of buffers
  config = trim(prefix) // '_reference_times' 
  call mpas_pool_get_config(domain % configs, config, config_results)
  copy = config_results
  number_of_buffers = 1
  b = scan(copy, ';')
  do while (b > 0)
    number_of_buffers = number_of_buffers + 1
    copy = copy(b+1:)
    b = scan(copy, ';')
  end do

  ! get our operation
  config = trim(prefix) // '_operation'
  call mpas_pool_get_config(domain % configs, config, config_results)
  if (config_results == 'avg') then
    operation = AVG_OP
    op = 'avg'
  else if (config_results == 'min') then
    operation = MIN_OP
    op = 'min'
  else if (config_results == 'max') then
    operation = MAX_OP
    op = 'max'
  else
    ! error if unknown operation
    call mpas_dmpar_global_abort('Error: unknown operation in time ' // &
      'averaging analysis member configuration.')
  end if

end subroutine start_init



!***********************************************************************
! routine get_alarms
!
!> \brief Read the namelist for timings
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!> This will read the namelist and get the strings and set the clocks
!> for the different timers to be used. The actual alarms are not set.
!-----------------------------------------------------------------------
subroutine get_alarms(domain, prefix, number_of_buffers, err)
  ! input variables
  integer, intent(in) :: number_of_buffers
  character (len=StrKIND) :: prefix

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: b, n
  character (len=StrKIND), pointer :: config_results
  character (len=StrKIND) :: config
  logical :: ok
  type (mpas_timeinterval_type) :: rem, zero

  ! configure start times - we don't have to check ok
  ! because the timer count is based on reference_times tokens
  config = trim(prefix) // '_reference_times' 
  call mpas_pool_get_config(domain % configs, config, config_results)
  call set_times(buffers, number_of_buffers, domain % clock, &
            START_TIMES, config_results, ok, err) 

  ! order matters, don't reorder these following ones!
  ! it matters because times/intervals can be configured to be equal
  ! to other ones

  ! configure reset intervals
  config = trim(prefix) // '_reset_intervals'
  call mpas_pool_get_config(domain % configs, config, config_results)
  call set_times(buffers, number_of_buffers, domain % clock, &
            RESET_INTERVALS, config_results, ok, err) 
  if (.not. ok) then
    call mpas_dmpar_global_abort('Error: number of times in ' // &
      'reset_intervals is not consistent with number of times ' // &
      'in reference_times in time series stats analysis member ' // &
      'configuration.')
  end if

  ! configure repeat intervals
  config = trim(prefix) // '_repeat_intervals'
  call mpas_pool_get_config(domain % configs, config, config_results)
  call set_times(buffers, number_of_buffers, domain % clock, &
            REPEAT_INTERVALS, config_results, ok, err) 
  if (.not. ok) then
    call mpas_dmpar_global_abort('Error: number of times in ' // &
      'repeat_intervals is not consistent with number of times ' // &
      'in reference_times in time series stats analysis member ' // &
      'configuration.')
  end if

  ! configure duration intervals
  config = trim(prefix) // '_duration_intervals'
  call mpas_pool_get_config(domain % configs, config, config_results)
  call set_times(buffers, number_of_buffers, domain % clock, &
            DURATION_INTERVALS, config_results, ok, err) 
  if (.not. ok) then
    call mpas_dmpar_global_abort('Error: number of times in ' // &
      'duration_intervals is not consistent with number of times ' // &
      'in reference_times in time series stats analysis member ' // &
      'configuration.')
  end if

  ! check if some of the time configuration is sensible
  call mpas_set_timeInterval(zero, s=0)

  do b = 1, number_of_buffers
    call mpas_interval_division(buffers(b) % start_time, &
       buffers(b) % repeat_interval, buffers(b) % reset_interval, n, rem) 

    if (n > 1 .or. (n == 1 .and. rem /= zero)) then
      write(stderrUnit,*) 'Warning: repeat_interval > ' // &
        'reset_interval in time averaging analysis member ' // &
        'configuration. Truncating repeat_interval.'
      buffers(b) % repeat_interval = buffers(b) % reset_interval
    end if

    call mpas_interval_division(buffers(b) % start_time, &
       buffers(b) % duration_interval, buffers(b) % repeat_interval, n, rem) 

    if (n > 1 .or. (n == 1 .and. rem /= zero)) then
      write(stderrUnit,*) 'Warning: duration_interval > ' // &
        'repeat_interval in time averaging analysis member ' // &
        'configuration. Truncating duration_interval.'
      buffers(b) % repeat_interval = buffers(b) % reset_interval
    end if
  end do
end subroutine get_alarms



!***********************************************************************
! routine set_alarms
!
!> \brief Set the alarms based on the clocks
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!> Alarms for the different timers are set, such that temporal
!> window alarms are configured.
!-----------------------------------------------------------------------
subroutine set_alarms(clock, instance, number_of_buffers, err)
  ! input variables
  integer, intent(in) :: number_of_buffers
  character (len=StrKIND) :: instance

  ! input/output variables
  type (mpas_clock_type), intent(inout) :: clock

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: b, repeat_n, duration_n, reset_n
  character (len=StrKIND) :: buffer
  type (mpas_time_type) :: current_time, when, &
    duration_time, repeat_time, reset_time
  type (mpas_timeinterval_type) :: elapsed, zero, &
    repeat_rem, duration_rem, reset_rem

  ! get current time
  current_time = mpas_get_clock_time(clock, MPAS_NOW, err)

  ! configure alarms
  do b = 1, number_of_buffers
    write(buffer, '(I0)') b

    ! see if we start in the future or we have already started
    if (current_time >= buffers(b) % start_time) then
       buffers(b) % started_flag = .true.
       ! TODO this needs to be false if do_restart
       buffers(b) % reset_flag = .true.

       ! no start alarm
       buffers(b) % start_alarm_ID = ''
    else
       buffers(b) % started_flag = .false.
       buffers(b) % reset_flag = .false.
 
       ! set the start alarm
       buffers(b) % start_alarm_ID = &
         'tavg_start' // trim(instance) // '_' // buffer
       call mpas_add_clock_alarm(clock, &
         buffers(b) % start_alarm_ID, &
         buffers(b) % start_time, ierr=err)
     end if

    !
    ! determine next alarm times
    !

    ! set next duration time
    when = buffers(b) % start_time + &
      buffers(b) % duration_interval ! duration is offset
    if (current_time > when) then
      elapsed = current_time - when
      call mpas_interval_division(when, elapsed, &
        buffers(b) % repeat_interval, & ! repeat is correct
        duration_n, duration_rem)
      duration_rem = buffers(b) % repeat_interval - duration_rem
      duration_time = current_time + duration_rem ! remainder of repeat
    else
      duration_time = buffers(b) % start_time + buffers(b) % duration_interval
      duration_n = 0
    end if

    ! set next repeat time
    when = buffers(b) % start_time + buffers(b) % repeat_interval
    if (current_time > when) then
      elapsed = current_time - when
      call mpas_interval_division(when, elapsed, &
        buffers(b) % repeat_interval, repeat_n, repeat_rem)
      repeat_rem = buffers(b) % repeat_interval - repeat_rem
      repeat_time = current_time + repeat_rem
    else
      repeat_time = buffers(b) % start_time + buffers(b) % repeat_interval
      repeat_n = 0
    end if

    ! set next reset time
    when = buffers(b) % start_time + buffers(b) % reset_interval
    if (current_time > when) then
      elapsed = current_time - when
      call mpas_interval_division(when, elapsed, &
        buffers(b) % reset_interval, reset_n, reset_rem)
      reset_rem = buffers(b) % reset_interval - reset_rem
      reset_time = current_time + reset_rem
    else
      reset_time = buffers(b) % start_time + buffers(b) % reset_interval
      reset_n = 0
    end if

    ! we're accumulating if we are in a window between duration and repeat
    buffers(b) % accumulate_flag = duration_n == repeat_n

    !
    ! set the reoccurring timers
    !
    buffers(b) % duration_alarm_ID = &
      'tavg_duration' // trim(instance) // '_' // buffer
    call mpas_add_clock_alarm(clock, &
      buffers(b) % duration_alarm_ID, &
      duration_time, & ! duration sets the offset
      buffers(b) % repeat_interval, ierr=err) ! but repeat sets the interval

    buffers(b) % repeat_alarm_ID = &
      'tavg_repeat' // trim(instance) // '_' // buffer
    call mpas_add_clock_alarm(clock, &
      buffers(b) % repeat_alarm_ID, &
      repeat_time, &
      buffers(b) % repeat_interval, ierr=err)

    buffers(b) % reset_alarm_ID = &
      'tavg_reset' // trim(instance) // '_' // buffer
    call mpas_add_clock_alarm(clock, &
      buffers(b) % reset_alarm_ID, &
      reset_time, &
      buffers(b) % reset_interval, ierr=err)
  end do
end subroutine set_alarms



!***********************************************************************
! routine modify_stream
!
!> \brief Remove existing variables and replace them with new ones
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!>  Given a stream name, this will remove the existing variables
!>  in a stream and replace them with similiarly named ones for
!>  their accumulation. It will also add xtime and optionally the mesh.
!-----------------------------------------------------------------------
subroutine modify_stream(domain, stream_name, number_of_variables, &
  number_of_buffers, instance, prefix, op, err)!{{{
  ! input variables
  integer, intent(in) :: number_of_variables, number_of_buffers
  character (len=StrKIND) :: stream_name, instance, prefix, op

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: v, b
  character (len=StrKIND) :: field, buffer, config, var
  logical, pointer :: copy_mesh
  type (Field0DReal), pointer :: src, dst

  ! allocate the variable information
  allocate(variables(number_of_variables))

  ! allocate the state for the buffers
  allocate(buffers(number_of_buffers))

  ! get the old field names 
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    stream_name, err)
  v = 1
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    stream_name, field))
    variables(v) % input_name = field
    v = v + 1
  end do
  
  ! remove the old ones from the stream 
  do v = 1, number_of_variables
    call mpas_stream_mgr_remove_field(domain % streamManager, &
      stream_name, variables(v) % input_name)
  end do

  ! add xtime to the stream
  call mpas_stream_mgr_add_field(domain % streamManager, &
    stream_name, 'xtime', ierr=err)

  ! optionally add mesh to stream
  config = trim(prefix) // '_add_mesh'
  call mpas_pool_get_config(domain % configs, config, copy_mesh)
  if (copy_mesh) then
    call mpas_stream_mgr_begin_iteration(domain % streamManager, &
      'mesh', err)
    do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
      'mesh', field))
      call mpas_stream_mgr_add_field(domain % streamManager, &
        stream_name, field, ierr=err)
    end do
  end if

  ! set up the variables
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    stream_name, err)
  do v = 1, number_of_variables 
    ! allocate space for the names of the outputs
    allocate(variables(v) % output_names(number_of_buffers))
    write(var, '(I0)') v

    ! get the info of the field
    call mpas_pool_get_field_info(domain % blocklist % allFields, &
      variables(v) % input_name, variables(v) % info)

    ! check if we can handle it
    if(.not. &
      ((variables(v) % info % fieldType == MPAS_POOL_REAL) &
      .or. &
      (variables(v) % info % fieldType == MPAS_POOL_INTEGER))) &
      then
      call mpas_dmpar_global_abort('Error: field "' // &
        trim(variables(v) % input_name) // '" listed in the ' // &
        'output stream, for time series stats analysis member ' // &
        'stream, is not real or integer.')
    end if

    ! allocate a number of fields and add field
    do b = 1, number_of_buffers
      ! create the name of the new field
      write(buffer, '(I0)') b
      field = 'time' // trim(instance) // '_' // &
        trim(op) // '_' // trim(buffer) // '_'
      variables(v) % output_names(b) = trim(field) // &
        variables(v) % input_name

      ! create the field and add to pool
      call add_new_field(variables(v) % info, &
        variables(v) % input_name, field, domain % blocklist % allFields)

      ! add the field to the stream
      call mpas_stream_mgr_add_field(domain % streamManager, &
        stream_name, variables(v) % output_names(b), ierr=err)
    end do
  end do ! number_of_variables

  ! add counters to stream
  do b = 1, size(buffers)
    ! create the name of the counter
    write(buffer, '(I0)') b
    field = 'time' // trim(instance) // '_' // trim(buffer) // '_counter'

    ! create counter and add to pool
    call mpas_pool_get_field(domain % blocklist % allFields, &
      'timeSeriesStatsCounter', src, 1)
    call mpas_duplicate_field(src, dst)
    dst % fieldName = field
    call mpas_pool_add_field(domain % blocklist % allFields, &
      field, dst)

    ! add counter to the stream
    buffers(b) % output_counter = field
    call mpas_stream_mgr_add_field(domain % streamManager, &
      stream_name, buffers(b) % output_counter, ierr=err)
  end do

end subroutine modify_stream!}}}



!***********************************************************************
! routine walk_string
!
!> \brief Walk a semicolon delimited string to find substrings
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!>  Walk a string delimited by semicolons and return the first substring
!>  from start index, and modify start to point at the next candidate.
!-----------------------------------------------------------------------
subroutine walk_string(next, substr, ok)!{{{
  ! input variables

  ! input/output variables
  character (len=StrKIND), intent(inout) :: next

  ! output variables
  character (len=StrKIND), intent(out) :: substr
  logical, intent(out) :: ok

  ! local variables
  integer :: i
  character (len=StrKIND) :: copy

  ! make a copy
  copy = trim(next)

  ! if there's anything in it other than whitespace, pass through
  i = verify(copy, ' ')
  ok = i > 0
  if (.not. ok) then
    return
  end if 
  copy = trim(next(i:))

  ! find the first semicolon and split
  i = scan(copy, ';')
  
  ! return that substring and the remainder
  if (i > 0) then
    substr = trim(copy(1:i-1))
    next = trim(copy(i+1:))
  else
    substr = trim(copy)
    next = ''
  end if
    
end subroutine walk_string!}}}



!***********************************************************************
! routine set_times
!
!> \brief Set a list of times
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!>  Walk a list of times delimited by spaces and set the time info
!>  for the buffer structure so that alarms can be set.
!-----------------------------------------------------------------------
subroutine set_times(buffers, number_of_buffers, clock, &
              which, config, ok, err)
  ! input variables
  integer, intent(in) :: number_of_buffers, which
  character (len=StrKIND), pointer, intent(in) :: config

  ! input/output variables
  type (time_buffer_type), dimension(:), intent(inout) :: buffers
  type (MPAS_Clock_type), intent(inout) :: clock

  ! output variables
  logical, intent(out) :: ok
  integer, intent(out) :: err 

  ! local variables
  character (len=StrKIND) :: next, time
  integer :: b 

  ! find the first time in the list
  next = config
  b = 0
  call walk_string(next, time, ok)

  ! while the time string is ok
  do while (ok)
    ! exit if we went over
    b = b + 1
    if (b > number_of_buffers) then
      exit
    end if

    ! set the time
    if (which == START_TIMES) then
      if (time == 'initial_time') then
        buffers(b) % start_time = mpas_get_clock_time(clock, &
          MPAS_START_TIME, err)
      else
        call mpas_set_time(buffers(b) % start_time, &
            dateTimeString=time, ierr=err)
      end if
    else if (which == DURATION_INTERVALS) then
      if (time == 'repeat_interval') then
        buffers(b) % duration_interval = buffers(b) % repeat_interval
      else
        call mpas_set_timeInterval(buffers(b) % duration_interval, &
            timeString=time, ierr=err)
      end if
    else if (which == REPEAT_INTERVALS) then
      if (time == 'reset_interval') then
        buffers(b) % repeat_interval = buffers(b) % reset_interval
      else
        call mpas_set_timeInterval(buffers(b) % repeat_interval, &
            timeString=time, ierr=err)
      end if
    else
      call mpas_set_timeInterval(buffers(b) % reset_interval, &
          timeString=time, ierr=err)
    end if 

    ! get the next time string
    call walk_string(next, time, ok)
  end do

  ! only ok if we parsed out as many as there are number of buffers
  ok = number_of_buffers == b
 end subroutine set_times



!***********************************************************************
! routine add_new_field
!
!> \brief Function to create a new field from an existing field
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!>  This routine conducts all initializations required for
!>  duplicating a field and adding it to the allFields pool.
!-----------------------------------------------------------------------
subroutine add_new_field(info, inname, prefix, pool)!{{{
  ! input variables
  type (mpas_pool_field_info_type), intent(in) :: info
  character (len=StrKIND), intent(in) :: inname, prefix

  ! input/output variables
  type (mpas_pool_type), intent(inout) :: pool 

  ! output variables

  ! local variables

  ! duplicate field and add new field to pool
  if (info % fieldType == MPAS_POOL_REAL) then
    if (info % nDims == 0) then
      call copy_field_0r(inname, pool, prefix)
    else if (info % nDims == 1) then
      call copy_field_1r(inname, pool, prefix)
    else if (info % nDims == 2) then
      call copy_field_2r(inname, pool, prefix)
    else if (info % nDims == 3) then
      call copy_field_3r(inname, pool, prefix)
    else if (info % nDims == 4) then
      call copy_field_4r(inname, pool, prefix)
    else
      call copy_field_5r(inname, pool, prefix)
    end if
  else
    if (info % nDims == 0) then
      call copy_field_0i(inname, pool, prefix)
    else if (info % nDims == 1) then
      call copy_field_1i(inname, pool, prefix)
    else if (info % nDims == 2) then
      call copy_field_2i(inname, pool, prefix)
    else
      call copy_field_3i(inname, pool, prefix)
    end if
  end if

end subroutine add_new_field!}}}



!***********************************************************************
! routine timer_checking 
!
!> \brief Timer functions to determine when to run
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!>  This routine conducts timer checking to determine if it 
!>  needs to run at this particular time. 
!-----------------------------------------------------------------------
subroutine timer_checking(clock, err)!{{{
  ! input variables

  ! input/output variables
  type (mpas_clock_type), intent(inout) :: clock 

  ! output variables
  integer, intent(out) :: err

  ! local variables
  integer :: b

  ! start procedure
  err = 0

  do b = 1, size(buffers)
    ! clear any resets
    if (buffers(b) % reset_flag) then
      if (buffers(b) % accumulate_flag) then
        buffers(b) % reset_flag = .false.
      end if
    end if

    ! see if the started alarm is ringing
    if (trim(buffers(b) % start_alarm_ID) /= '') then
      if (mpas_is_alarm_ringing(clock, &
          buffers(b) % start_alarm_ID, ierr=err)) then
        call mpas_reset_clock_alarm(clock, &
          buffers(b) % start_alarm_ID, ierr=err)
        buffers(b) % started_flag = .true. 
        buffers(b) % reset_flag = .true.
        buffers(b) % accumulate_flag = .true.
      end if
    end if

    ! if we aren't started, cycle to next buffer
    if (.not. buffers(b) % started_flag) then
      cycle
    end if

    ! check various other alarms
    ! see if we need to reset
    if(mpas_is_alarm_ringing(clock, &
      buffers(b) % reset_alarm_ID, ierr=err)) then
      call mpas_reset_clock_alarm(clock, &
        buffers(b) % reset_alarm_ID, ierr=err)
      buffers(b) % reset_flag = .true.
    end if

    ! turn off accumulation
    !
    ! duration needs to be >= 2 * compute_interval 
    ! (a series can only be 2 or more)
    if (mpas_is_alarm_ringing(clock, &
        buffers(b) % duration_alarm_ID, ierr=err)) then
      call mpas_reset_clock_alarm(clock, &
        buffers(b) % duration_alarm_ID, ierr=err)
      buffers(b) % accumulate_flag = .false.
    end if

    ! turn on accumulation 
    ! (this is second, in case the duration and repeat
    ! overlaps on the same timer)
    if (mpas_is_alarm_ringing(clock, &
        buffers(b) % repeat_alarm_ID, ierr=err)) then
      call mpas_reset_clock_alarm(clock, &
        buffers(b) % repeat_alarm_ID, ierr=err)
      buffers(b) % accumulate_flag = .true.
    end if

  end do
end subroutine timer_checking!}}}



!***********************************************************************
! routine typed_operate
!
!> \brief Do the operation, but switch on run-time type
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!>  Since we don't know the type of the array, we need to do some
!>  run-time type switching based on the type of the array.
!-----------------------------------------------------------------------
subroutine typed_operate(block, v, operation)!{{{
  ! input variables
  type (block_type), pointer, intent(in) :: block
  integer, intent(in) :: v, operation

  ! input/output variables

  ! output variables

  ! local variables

  ! switch based on the type, dimensionality, and operation
  if (variables(v) % info % fieldType == MPAS_POOL_REAL) then
    if (variables(v) % info % nDims  == 0) then
      if (operation == AVG_OP) then
        call operate0r_avg(block, variables(v))
      else if (operation == MIN_OP) then
        call operate0r_min(block, variables(v))
      else
        call operate0r_max(block, variables(v))
      end if 
    else if (variables(v) % info % nDims == 1) then
      if (operation == AVG_OP) then
        call operate1r_avg(block, variables(v))
      else if (operation == MIN_OP) then
        call operate1r_min(block, variables(v))
      else
        call operate1r_max(block, variables(v))
      end if 
    else if (variables(v) % info % nDims == 2) then
      if (operation == AVG_OP) then
        call operate2r_avg(block, variables(v))
      else if (operation == MIN_OP) then
        call operate2r_min(block, variables(v))
      else
        call operate2r_max(block, variables(v))
      end if 
    else if (variables(v) % info % nDims == 3) then
      if (operation == AVG_OP) then
        call operate3r_avg(block, variables(v))
      else if (operation == MIN_OP) then
        call operate3r_min(block, variables(v))
      else
        call operate3r_max(block, variables(v))
      end if 
    else if (variables(v) % info % nDims == 4) then
      if (operation == AVG_OP) then
        call operate4r_avg(block, variables(v))
      else if (operation == MIN_OP) then
        call operate4r_min(block, variables(v))
      else
        call operate4r_max(block, variables(v))
      end if 
    else
      if (operation == AVG_OP) then
        call operate5r_avg(block, variables(v))
      else if (operation == MIN_OP) then
        call operate5r_min(block, variables(v))
      else
        call operate5r_max(block, variables(v))
      end if 
    end if
  else 
    if (variables(v) % info % nDims == 0) then
      if (operation == AVG_OP) then
        call operate0i_avg(block, variables(v))
      else if (operation == MIN_OP) then
        call operate0i_min(block, variables(v))
      else
        call operate0i_max(block, variables(v))
      end if 
    else if (variables(v) % info % nDims == 1) then
      if (operation == AVG_OP) then
        call operate1i_avg(block, variables(v))
      else if (operation == MIN_OP) then
        call operate1i_min(block, variables(v))
      else
        call operate1i_max(block, variables(v))
      end if 
    else if (variables(v) % info % nDims == 2) then
      if (operation == AVG_OP) then
        call operate2i_avg(block, variables(v))
      else if (operation == MIN_OP) then
        call operate2i_min(block, variables(v))
      else
        call operate2i_max(block, variables(v))
      end if 
    else
      if (operation == AVG_OP) then
        call operate3i_avg(block, variables(v))
      else if (operation == MIN_OP) then
        call operate3i_min(block, variables(v))
      else
        call operate3i_max(block, variables(v))
      end if 
    end if
  end if
end subroutine typed_operate!}}}



!***********************************************************************
! routine copy_field_X 
!
!> \brief Functions to create a new field from an existing field
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!>  This routine conducts initializations required for
!>  duplicating a field and adding it to the allFields pool based on type.
!-----------------------------------------------------------------------

subroutine copy_field_0r(inname, pool, prefix)!{{{
  character (len=StrKIND), intent(in) :: inname, prefix
  type (mpas_pool_type), intent(inout) :: pool

  type (field0DReal), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = trim(prefix) // dst % fieldName

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(prefix) // &
        dst % constituentNames(i)
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_0r!}}}

subroutine copy_field_1r(inname, pool, prefix)!{{{
  character (len=StrKIND), intent(in) :: inname, prefix
  type (mpas_pool_type), intent(inout) :: pool

  type (field1DReal), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = trim(prefix) // dst % fieldName

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(prefix) // &
        dst % constituentNames(i)
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_1r!}}}

subroutine copy_field_2r(inname, pool, prefix)!{{{
  character (len=StrKIND), intent(in) :: inname, prefix
  type (mpas_pool_type), intent(inout) :: pool

  type (field2DReal), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = trim(prefix) // dst % fieldName

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(prefix) // &
        dst % constituentNames(i)
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_2r!}}}

subroutine copy_field_3r(inname, pool, prefix)!{{{
  character (len=StrKIND), intent(in) :: inname, prefix
  type (mpas_pool_type), intent(inout) :: pool

  type (field3DReal), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = trim(prefix) // dst % fieldName

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(prefix) // &
        dst % constituentNames(i)
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_3r!}}}

subroutine copy_field_4r(inname, pool, prefix)!{{{
  character (len=StrKIND), intent(in) :: inname, prefix
  type (mpas_pool_type), intent(inout) :: pool

  type (field4DReal), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = trim(prefix) // dst % fieldName

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(prefix) // &
        dst % constituentNames(i)
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_4r!}}}

subroutine copy_field_5r(inname, pool, prefix)!{{{
  character (len=StrKIND), intent(in) :: inname, prefix
  type (mpas_pool_type), intent(inout) :: pool

  type (field5DReal), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = trim(prefix) // dst % fieldName

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(prefix) // &
        dst % constituentNames(i)
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_5r!}}}

subroutine copy_field_0i(inname, pool, prefix)!{{{
  character (len=StrKIND), intent(in) :: inname, prefix
  type (mpas_pool_type), intent(inout) :: pool

  type (field0DInteger), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = trim(prefix) // dst % fieldName

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(prefix) // &
        dst % constituentNames(i)
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_0i!}}}

subroutine copy_field_1i(inname, pool, prefix)!{{{
  character (len=StrKIND), intent(in) :: inname, prefix
  type (mpas_pool_type), intent(inout) :: pool

  type (field1DInteger), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = trim(prefix) // dst % fieldName

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(prefix) // &
        dst % constituentNames(i)
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_1i!}}}

subroutine copy_field_2i(inname, pool, prefix)!{{{
  character (len=StrKIND), intent(in) :: inname, prefix
  type (mpas_pool_type), intent(inout) :: pool

  type (field2DInteger), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = trim(prefix) // dst % fieldName

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(prefix) // &
        dst % constituentNames(i)
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_2i!}}}

subroutine copy_field_3i(inname, pool, prefix)!{{{
  character (len=StrKIND), intent(in) :: inname, prefix
  type (mpas_pool_type), intent(inout) :: pool

  type (field3DInteger), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = trim(prefix) // dst % fieldName

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(prefix) // &
        dst % constituentNames(i)
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_3i!}}}


!***********************************************************************
! routine operateX_Y
!
!> \brief Series of subroutines to support operations on run-time types
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details 
!>  These subroutines encapsulate the different opertions that can occur
!>  based on the run-time types. (This would likely be
!>  instantiated generics/templates in other languages.)
!>
!>  Averaging is done by multiplying out and dividing such that
!>  the average state is always in a normalized form -- while
!>  this could (will) cause more error in the long run, it does
!>  mean that other AMs will be able to use this data and it will
!>  always be prenormalized (it also means that we don't have to
!>  have a special case of normalizing the data before writing it
!>  to disk).
!-----------------------------------------------------------------------
subroutine operate0r_avg (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
        out_array = (out_array * &
                (buffers(b) % total_accum - 1) + in_array) &
                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate0r_avg

subroutine operate1r_avg (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
        out_array = (out_array * &
                (buffers(b) % total_accum - 1) + in_array) &
                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate1r_avg

subroutine operate2r_avg (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
        out_array = (out_array * &
                (buffers(b) % total_accum - 1) + in_array) &
                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate2r_avg

subroutine operate3r_avg (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
        out_array = (out_array * &
                (buffers(b) % total_accum - 1) + in_array) &
                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate3r_avg

subroutine operate4r_avg (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:,:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
        out_array = (out_array * &
                (buffers(b) % total_accum - 1) + in_array) &
                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate4r_avg

subroutine operate5r_avg (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:,:,:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
        out_array = (out_array * &
                (buffers(b) % total_accum - 1) + in_array) &
                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate5r_avg

subroutine operate0i_avg (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  integer, pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
        out_array = (out_array * &
                (buffers(b) % total_accum - 1) + in_array) &
                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate0i_avg

subroutine operate1i_avg (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  integer, dimension(:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
        out_array = (out_array * &
                (buffers(b) % total_accum - 1) + in_array) &
                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate1i_avg

subroutine operate2i_avg (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  integer, dimension(:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
        out_array = (out_array * &
                (buffers(b) % total_accum - 1) + in_array) &
                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate2i_avg

subroutine operate3i_avg (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  integer, dimension(:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
        out_array = (out_array * &
                (buffers(b) % total_accum - 1) + in_array) &
                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate3i_avg

subroutine operate0r_min (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate0r_min

subroutine operate1r_min (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate1r_min

subroutine operate2r_min (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate2r_min

subroutine operate3r_min (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate3r_min

subroutine operate4r_min (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:,:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate4r_min

subroutine operate5r_min (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:,:,:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate5r_min

subroutine operate0i_min (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  integer, pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate0i_min

subroutine operate1i_min (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  integer, dimension(:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate1i_min

subroutine operate2i_min (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  integer, dimension(:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate2i_min

subroutine operate3i_min (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  integer, dimension(:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate3i_min

subroutine operate0r_max (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate0r_max

subroutine operate1r_max (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate1r_max

subroutine operate2r_max (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate2r_max

subroutine operate3r_max (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate3r_max

subroutine operate4r_max (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:,:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate4r_max

subroutine operate5r_max (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  real (kind=RKIND), dimension(:,:,:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate5r_max

subroutine operate0i_max (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  integer, pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate0i_max

subroutine operate1i_max (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  integer, dimension(:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate1i_max

subroutine operate2i_max (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  integer, dimension(:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate2i_max

subroutine operate3i_max (start_block, tvar)
  type (block_type), pointer, intent(in) :: start_block
  type (time_variable_type), intent(inout) :: tvar

  integer, dimension(:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block 

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      tvar % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (.not. buffers(b) % accumulate_flag) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          tvar % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag) then
        out_array = in_array
      else 
!        out_array = (out_array * &
!                (buffers(b) % total_accum - 1) + in_array) &
!                / buffers(b) % total_accum ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if 
    end do 

    block => block % next 
  end do 
end subroutine operate3i_max

end module ocn_time_series_stats
! vim: foldmethod=marker
