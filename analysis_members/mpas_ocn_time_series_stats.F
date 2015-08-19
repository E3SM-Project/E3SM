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
!> \date   March 2, 2015
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
    integer :: total_accum

    type (MPAS_Time_type) :: start_time
    type (MPAS_TimeInterval_type) :: duration_interval
    type (MPAS_TimeInterval_type) :: repeat_interval
    type (MPAS_TimeInterval_type) :: reset_interval

    ! alarm IDs
    character (len=StrKIND) :: start_alarm_ID
    character (len=StrKIND) :: repeat_alarm_ID
    character (len=StrKIND) :: duration_alarm_ID
    character (len=StrKIND) :: reset_alarm_ID
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
!> \date    March 2, 2015
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_init_time_series_stats(domain, err)!{{{

  ! input variables
  !-----------------------------------------------------------------

  ! input/output variables
  !-----------------------------------------------------------------
  type (domain_type), intent(inout) :: domain

  ! output variables
  !-----------------------------------------------------------------
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  !-----------------------------------------------------------------
  integer :: v, b
  character (len=StrKIND), pointer :: config_results
  logical, pointer :: copy_mesh
  integer :: number_of_variables, number_of_buffers
  character (len=StrKIND) :: stream_str, prefix_str, &
    config_str, buffer_str, op_str, var_str, field
  logical :: ok

  ! start procedure
  !-----------------------------------------------------------------
  err = 0

  ! TODO do restart

  ! string representation
  ! TODO placeholder for some unique ID if this code is replicated
  !    per multiple AMs for multiple streams
  stream_str = ''
  prefix_str = 'config_AM_timeSeriesStats' // trim(stream_str)

  ! get our operation
  config_str = trim(prefix_str) // '_operation'
  call mpas_pool_get_config(domain % configs, config_str, config_results)
  if (config_results .eq. 'avg') then
    operation = AVG_OP
    op_str = 'avg'
  else if (config_results .eq. 'min') then
    operation = MIN_OP
    op_str = 'min'
  else if (config_results .eq. 'max') then
    operation = MAX_OP
    op_str = 'max'
  else
    ! error if unknown operation
    call mpas_dmpar_global_abort('Error: unknown operation in time ' // &
      'averaging analysis member configuration.')
  end if

  ! count string tokens
  config_str = trim(prefix_str) // '_initial_times' 
  call mpas_pool_get_config(domain % configs, config_str, config_results)
  field = config_results
  number_of_buffers = 1
  b = scan(field, ';')
  do while (b .gt. 0)
    number_of_buffers = number_of_buffers + 1
    field = field(b+1:)
    b = scan(field, ';')
  end do

  ! get the stream name
  config_str = trim(prefix_str) // '_stream_name'
  call mpas_pool_get_config(domain % configs, config_str, stream_name)

  if (stream_name .eq. 'none') then
    call mpas_dmpar_global_abort('Error: stream cannot be "none" ' // &
      'for time series stats.')
  end if

  ! set up all of the timing
  !

  ! allocate the state for the buffers
  allocate(buffers(number_of_buffers))

  ! configure start times
  config_str = trim(prefix_str) // '_initial_times' 
  call mpas_pool_get_config(domain % configs, config_str, config_results)
  call set_times(buffers, number_of_buffers, domain % clock, &
            START_TIMES, config_results, ok, err) 

  ! order matters, don't reorder these following ones!
  ! it matters because times/intervals can be configured to be equal
  ! to other ones

  ! configure reset intervals
  config_str = trim(prefix_str) // '_reset_intervals'
  call mpas_pool_get_config(domain % configs, config_str, config_results)
  call set_times(buffers, number_of_buffers, domain % clock, &
            RESET_INTERVALS, config_results, ok, err) 
  if (.not. ok) then
    call mpas_dmpar_global_abort('Error: number of times in ' // &
      'reset_intervals is not consistent with number of times ' // &
      'in initial_times in time series stats analysis member ' // &
      'configuration.')
  end if

  ! configure repeat intervals
  config_str = trim(prefix_str) // '_repeat_intervals'
  call mpas_pool_get_config(domain % configs, config_str, config_results)
  call set_times(buffers, number_of_buffers, domain % clock, &
            REPEAT_INTERVALS, config_results, ok, err) 
  if (.not. ok) then
    call mpas_dmpar_global_abort('Error: number of times in ' // &
      'repeat_intervals is not consistent with number of times ' // &
      'in initial_times in time series stats analysis member ' // &
      'configuration.')
  end if

  ! configure duration intervals
  config_str = trim(prefix_str) // '_duration_intervals'
  call mpas_pool_get_config(domain % configs, config_str, config_results)
  call set_times(buffers, number_of_buffers, domain % clock, &
            DURATION_INTERVALS, config_results, ok, err) 
  if (.not. ok) then
    call mpas_dmpar_global_abort('Error: number of times in ' // &
      'duration_intervals is not consistent with number of times ' // &
      'in initial_times in time series stats analysis member ' // &
      'configuration.')
  end if

  ! check if the configuration is sensible
  do b = 1, number_of_buffers
    if (buffers(b) % repeat_interval .gt. &
       buffers(b) % reset_interval) then
      write(stderrUnit,*) 'Warning: repeat_interval > ' // &
        'reset_interval in time averaging analysis member ' // &
        'configuration. Truncating repeat_interval.'
      buffers(b) % repeat_interval = buffers(b) % reset_interval
    end if

    if (buffers(b) % duration_interval .gt. &
       buffers(b) % repeat_interval) then
      write(stderrUnit,*) 'Warning: duration_interval > ' // &
        'repeat_interval in time averaging analysis member ' // &
        'configuration. Truncating duration_interval.'
      buffers(b) % repeat_interval = buffers(b) % reset_interval
    end if
  end do

  ! 
  ! OK, if we got this far, then we should be able to allocate memory
  ! and set up the timers and variables that we will analyze 
  !

  ! count the number of variables
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    stream_name, err)
  number_of_variables = 0
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    stream_name, field))
    number_of_variables = number_of_variables + 1
  end do

  ! allocate the variable information
  allocate(variables(number_of_variables))

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
  config_str = trim(prefix_str) // '_add_mesh'
  call mpas_pool_get_config(domain % configs, config_str, copy_mesh)
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
    write(var_str, '(I0)') v

    ! get the info of the field
    call mpas_pool_get_field_info(domain % blocklist % allFields, &
      variables(v) % input_name, variables(v) % info)

    ! check if we can handle it
    if(.not. &
      ((variables(v) % info % fieldType .eq. MPAS_POOL_REAL) &
      .or. &
      (variables(v) % info % fieldType .eq. MPAS_POOL_INTEGER))) &
      then
      call mpas_dmpar_global_abort('Error: field "' // &
        trim(variables(v) % input_name) // '" listed in the ' // &
        'output stream, for time series stats analysis member ' // &
        'stream, is not real or integer.')
    end if

    ! allocate a number of fields and add field
    do b = 1, number_of_buffers
      ! create the name of the new field
      write(buffer_str, '(I0)') b
      field = 'time' // trim(stream_str) // '_' // &
        trim(op_str) // '_' // trim(buffer_str) // '_'
      variables(v) % output_names(b) = trim(field) // &
        variables(v) % input_name

      ! create the field and add to pool
      call add_new_field(variables(v) % info, &
        variables(v) % input_name, field, &
        domain % blocklist % allFields)

      ! add the field to the stream
      call mpas_stream_mgr_add_field(domain % streamManager, &
        stream_name, variables(v) % output_names(b), ierr=err)
    end do

  end do ! number_of_variables

  ! configure alarms
  ! TODO modify the alarms based on do_restart
  do b = 1, number_of_buffers
    write(buffer_str, '(I0)') b
    buffers(b) % start_alarm_ID = &
      'tavg_start' // trim(stream_str) // '_' // buffer_str
    call mpas_add_clock_alarm(domain % clock, &
      buffers(b) % start_alarm_ID, &
      buffers(b) % start_time, ierr=err)

    buffers(b) % repeat_alarm_ID = &
      'tavg_repeat' // trim(stream_str) // '_' // buffer_str
    call mpas_add_clock_alarm(domain % clock, &
      buffers(b) % repeat_alarm_ID, &
      buffers(b) % start_time + &
      buffers(b) % repeat_interval, &
      buffers(b) % repeat_interval, ierr=err)

    buffers(b) % duration_alarm_ID = &
      'tavg_duration' // trim(stream_str) // '_' // buffer_str
    call mpas_add_clock_alarm(domain % clock, &
      buffers(b) % duration_alarm_ID, &
      buffers(b) % start_time + &
      buffers(b) % duration_interval, & 
      buffers(b) % repeat_interval, ierr=err)

    ! reset at start
    buffers(b) % reset_alarm_ID = &
      'tavg_reset' // trim(stream_str) // '_' // buffer_str
    call mpas_add_clock_alarm(domain % clock, &
      buffers(b) % reset_alarm_ID, &
      buffers(b) % start_time, &
      buffers(b) % reset_interval, ierr=err)

  end do

  ! set initial flags
  do b = 1, number_of_buffers
    buffers(b) % started_flag = .false.
    buffers(b) % reset_flag = .false.
    buffers(b) % accumulate_flag = .false.
  end do

  ! set flags based on initial times
  call timer_checking(domain, err)
end subroutine ocn_init_time_series_stats!}}}



!***********************************************************************
! routine ocn_compute_time_series_stats
!
!> \brief Compute MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 2, 2015
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_compute_time_series_stats(domain, timeLevel, err)!{{{
  ! input variables
  !-----------------------------------------------------------------
  integer, intent(in) :: timeLevel

  ! input/output variables
  !-----------------------------------------------------------------
  type (domain_type), intent(inout) :: domain

  ! output variables
  !-----------------------------------------------------------------
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  !-----------------------------------------------------------------
  integer :: i, v, b

  ! start procedure
  !-----------------------------------------------------------------
  err = 0

  ! update the counter
  do b = 1, size(buffers)
    if (buffers(b) % reset_flag) then
      buffers(b) % total_accum = 1
    else if (buffers(b) % accumulate_flag) then
      buffers(b) % total_accum = buffers(b) % total_accum + 1
    end if
  end do

  ! do all of the operations
  do v = 1, size(variables)
    call typed_operate(domain % blocklist, variables(v), operation)
  end do

  ! clear any resets
  do b = 1, size(buffers)
    buffers(b) % reset_flag = .false.
  end do

  ! do all of the time checking and flag setting
  call timer_checking(domain, err)
end subroutine ocn_compute_time_series_stats!}}}



!***********************************************************************
! routine ocn_restart_time_series_stats
!
!> \brief Save restart for MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 2, 2015
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_restart_time_series_stats(domain, err)!{{{

  ! input variables
  !-----------------------------------------------------------------

  ! input/output variables
  !-----------------------------------------------------------------
  type (domain_type), intent(inout) :: domain

  ! output variables
  !-----------------------------------------------------------------
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  !-----------------------------------------------------------------

  ! start procedure
  !-----------------------------------------------------------------
  err = 0

  ! TODO save data to restart and accumulate

end subroutine ocn_restart_time_series_stats!}}}



!***********************************************************************
! routine ocn_finalize_time_series_stats
!
!> \brief Finalize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 2, 2015
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine ocn_finalize_time_series_stats(domain, err)!{{{

  ! input variables
  !-----------------------------------------------------------------

  ! input/output variables
  !-----------------------------------------------------------------
  type (domain_type), intent(inout) :: domain

  ! output variables
  !-----------------------------------------------------------------
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  !-----------------------------------------------------------------
  integer :: i, v

  ! start procedure
  !-----------------------------------------------------------------
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
! routine walk_string
!
!> \brief Walk a semicolon delimited string to find substrings
!> \author  Jon Woodring
!> \date    March 2, 2015
!> \details 
!>  Walk a string delimited by semicolons and return the first substring
!>  from start index, and modify start to point at the next candidate.
!-----------------------------------------------------------------------
subroutine walk_string(next, substr, ok)!{{{
  ! input variables
  !-----------------------------------------------------------------

  ! input/output variables
  !-----------------------------------------------------------------
  character (len=StrKIND), intent(inout) :: next

  ! output variables
  !-----------------------------------------------------------------
  character (len=StrKIND), intent(out) :: substr
  logical, intent(out) :: ok

  ! local variables
  !-----------------------------------------------------------------
  integer :: i
  character (len=StrKIND) :: copy

  ! make a copy
  copy = trim(next)

  ! if there's anything in it other than whitespace, pass through
  i = verify(copy, ' ')
  ok = i .gt. 0
  if (.not. ok) then
    return
  end if 
  copy = trim(next(i:))

  ! find the first semicolon and split
  i = scan(copy, ';')
  
  ! return that substring and the remainder
  if (i .gt. 0) then
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
!> \date    March 2, 2015
!> \details 
!>  Walk a list of times delimited by spaces and set the time info
!>  for the buffer structure so that alarms can be set.
!-----------------------------------------------------------------------
subroutine set_times(buffers, number_of_buffers, clock, &
              which, config_str, ok, err)
  ! input variables
  !-----------------------------------------------------------------
  integer, intent(in) :: number_of_buffers, which
  character (len=StrKIND), pointer, intent(in) :: config_str

  ! input/output variables
  !-----------------------------------------------------------------
  type (time_buffer_type), dimension(:), intent(inout) :: buffers
  type (MPAS_Clock_type), intent(inout) :: clock

  ! output variables
  !-----------------------------------------------------------------
  logical, intent(out) :: ok
  integer, intent(out) :: err 

  ! local variables
  !-----------------------------------------------------------------
  character (len=StrKIND) :: next_str, time_str
  integer :: b 

  ! find the first time in the list
  next_str = config_str
  b = 0
  call walk_string(next_str, time_str, ok)

  ! while the time string is ok
  do while (ok)
    ! exit if we went over
    b = b + 1
    if (b .gt. number_of_buffers) then
      exit
    end if

    ! set the time
    if (which .eq. START_TIMES) then
      if (time_str .eq. 'initial_time') then
        buffers(b) % start_time = mpas_get_clock_time(clock, &
          MPAS_NOW, err)
      else
        call mpas_set_time(buffers(b) % start_time, &
            dateTimeString=time_str, ierr=err)
      end if
    else if (which .eq. DURATION_INTERVALS) then
      if (time_str .eq. 'repeat_interval') then
        buffers(b) % duration_interval = buffers(b) % repeat_interval
      else
        call mpas_set_timeInterval(buffers(b) % duration_interval, &
            timeString=time_str, ierr=err)
      end if
    else if (which .eq. REPEAT_INTERVALS) then
      if (time_str .eq. 'reset_interval') then
        buffers(b) % repeat_interval = buffers(b) % reset_interval
      else
        call mpas_set_timeInterval(buffers(b) % repeat_interval, &
            timeString=time_str, ierr=err)
      end if
    else
      call mpas_set_timeInterval(buffers(b) % reset_interval, &
          timeString=time_str, ierr=err)
    end if 

    ! get the next time string
    call walk_string(next_str, time_str, ok)
  end do

  ! only ok if we parsed out as many as there are number of buffers
  ok = number_of_buffers .eq. b
 end subroutine set_times



!***********************************************************************
! routine add_new_field
!
!> \brief Function to create a new field from an existing field
!> \author  Jon Woodring
!> \date    March 2, 2015
!> \details 
!>  This routine conducts all initializations required for
!>  duplicating a field and adding it to the allFields pool.
!-----------------------------------------------------------------------
subroutine add_new_field(info, inname, prefix, pool)!{{{
  ! input variables
  !-----------------------------------------------------------------
  type (mpas_pool_field_info_type), intent(in) :: info
  character (len=StrKIND), intent(in) :: inname, prefix

  ! input/output variables
  !-----------------------------------------------------------------
  type (mpas_pool_type), intent(inout) :: pool 

  ! output variables
  !-----------------------------------------------------------------

  ! local variables
  !-----------------------------------------------------------------

  ! start procedure
  !-----------------------------------------------------------------

  ! duplicate field and add new field to pool
  if (info % fieldType .eq. MPAS_POOL_REAL) then
    if (info % nDims .eq. 0) then
      call copy_field_0r(inname, pool, prefix)
    else if (info % nDims .eq. 1) then
      call copy_field_1r(inname, pool, prefix)
    else if (info % nDims .eq. 2) then
      call copy_field_2r(inname, pool, prefix)
    else if (info % nDims .eq. 3) then
      call copy_field_3r(inname, pool, prefix)
    else if (info % nDims .eq. 4) then
      call copy_field_4r(inname, pool, prefix)
    else
      call copy_field_5r(inname, pool, prefix)
    end if
  else
    if (info % nDims .eq. 0) then
      call copy_field_0i(inname, pool, prefix)
    else if (info % nDims .eq. 1) then
      call copy_field_1i(inname, pool, prefix)
    else if (info % nDims .eq. 2) then
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
!> \date    March 2, 2015
!> \details 
!>  This routine conducts timer checking to determine if it 
!>  needs to run at this particular time. 
!-----------------------------------------------------------------------
subroutine timer_checking(domain, err)!{{{
  ! input variables
  !-----------------------------------------------------------------

  ! input/output variables
  !-----------------------------------------------------------------
  type (domain_type), intent(inout) :: domain

  ! output variables
  !-----------------------------------------------------------------
  integer, intent(out) :: err

  ! local variables
  !-----------------------------------------------------------------
  integer :: b

  ! start procedure
  !-----------------------------------------------------------------
  err = 0

  do b = 1, size(buffers)
    ! see if the started alarm is ringing
    if (mpas_is_alarm_ringing(domain % clock, &
        buffers(b) % start_alarm_ID, ierr=err)) then
      call mpas_reset_clock_alarm(domain % clock, &
        buffers(b) % start_alarm_ID, ierr=err)
      buffers(b) % started_flag = .true. 
      buffers(b) % accumulate_flag = .true.
    end if

    ! if we aren't started, continue to next buffer
    if (.not. buffers(b) % started_flag) then
      continue
    end if

    ! check various other alarms
    ! see if we need to reset
    if(mpas_is_alarm_ringing(domain % clock, &
      buffers(b) % reset_alarm_ID, ierr=err)) then
      call mpas_reset_clock_alarm(domain % clock, &
        buffers(b) % reset_alarm_ID, ierr=err)
      buffers(b) % reset_flag = .true.
    end if

    ! turn off accumulation
    !
    ! duration needs to be >= 2 * compute_interval 
    ! (a series can only be 2 or more)
    if (mpas_is_alarm_ringing(domain % clock, &
        buffers(b) % duration_alarm_ID, ierr=err)) then
      call mpas_reset_clock_alarm(domain % clock, &
        buffers(b) % duration_alarm_ID, ierr=err)
      buffers(b) % accumulate_flag = .false.
    end if

    ! turn on accumulation 
    ! (this is second, in case the duration and repeat
    ! overlaps on the same timer)
    if (mpas_is_alarm_ringing(domain % clock, &
        buffers(b) % repeat_alarm_ID, ierr=err)) then
      call mpas_reset_clock_alarm(domain % clock, &
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
!> \date    March 2, 2015
!> \details 
!>  Since we don't know the type of the array, we need to do some
!>  run-time type switching based on the type of the array.
!-----------------------------------------------------------------------
subroutine typed_operate(block, tvar, operation)!{{{
  ! input variables
  !-----------------------------------------------------------------
  type (block_type), pointer, intent(in) :: block
  integer, intent(in) :: operation

  ! input/output variables
  !-----------------------------------------------------------------
  type (time_variable_type), intent(inout) :: tvar

  ! output variables
  !-----------------------------------------------------------------

  ! local variables
  !-----------------------------------------------------------------

  ! switch based on the type, dimensionality, and operation
  if (tvar % info % fieldType == MPAS_POOL_REAL) then
    if (tvar % info % nDims  == 0) then
      if (operation .eq. AVG_OP) then
        call operate0r_avg(block, tvar)
      else if (operation .eq. MIN_OP) then
        call operate0r_min(block, tvar)
      else
        call operate0r_max(block, tvar)
      end if 
    else if (tvar % info % nDims == 1) then
      if (operation .eq. AVG_OP) then
        call operate1r_avg(block, tvar)
      else if (operation .eq. MIN_OP) then
        call operate1r_min(block, tvar)
      else
        call operate1r_max(block, tvar)
      end if 
    else if (tvar % info % nDims == 2) then
      if (operation .eq. AVG_OP) then
        call operate2r_avg(block, tvar)
      else if (operation .eq. MIN_OP) then
        call operate2r_min(block, tvar)
      else
        call operate2r_max(block, tvar)
      end if 
    else if (tvar % info % nDims == 3) then
      if (operation .eq. AVG_OP) then
        call operate3r_avg(block, tvar)
      else if (operation .eq. MIN_OP) then
        call operate3r_min(block, tvar)
      else
        call operate3r_max(block, tvar)
      end if 
    else if (tvar % info % nDims == 4) then
      if (operation .eq. AVG_OP) then
        call operate4r_avg(block, tvar)
      else if (operation .eq. MIN_OP) then
        call operate4r_min(block, tvar)
      else
        call operate4r_max(block, tvar)
      end if 
    else
      if (operation .eq. AVG_OP) then
        call operate5r_avg(block, tvar)
      else if (operation .eq. MIN_OP) then
        call operate5r_min(block, tvar)
      else
        call operate5r_max(block, tvar)
      end if 
    end if
  else 
    if (tvar % info % nDims == 0) then
      if (operation .eq. AVG_OP) then
        call operate0i_avg(block, tvar)
      else if (operation .eq. MIN_OP) then
        call operate0i_min(block, tvar)
      else
        call operate0i_max(block, tvar)
      end if 
    else if (tvar % info % nDims == 1) then
      if (operation .eq. AVG_OP) then
        call operate1i_avg(block, tvar)
      else if (operation .eq. MIN_OP) then
        call operate1i_min(block, tvar)
      else
        call operate1i_max(block, tvar)
      end if 
    else if (tvar % info % nDims == 2) then
      if (operation .eq. AVG_OP) then
        call operate2i_avg(block, tvar)
      else if (operation .eq. MIN_OP) then
        call operate2i_min(block, tvar)
      else
        call operate2i_max(block, tvar)
      end if 
    else
      if (operation .eq. AVG_OP) then
        call operate3i_avg(block, tvar)
      else if (operation .eq. MIN_OP) then
        call operate3i_min(block, tvar)
      else
        call operate3i_max(block, tvar)
      end if 
    end if
  end if
end subroutine typed_operate!}}}



!***********************************************************************
! routine copy_field_X 
!
!> \brief Functions to create a new field from an existing field
!> \author  Jon Woodring
!> \date    March 2, 2015
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
!> \date    March 2, 2015
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

        out_array = (out_array * \
                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

        out_array = (out_array * \
                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

        out_array = (out_array * \
                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

        out_array = (out_array * \
                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

        out_array = (out_array * \
                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

        out_array = (out_array * \
                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

        out_array = (out_array * \
                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

        out_array = (out_array * \
                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

        out_array = (out_array * \
                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

        out_array = (out_array * \
                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
      if (buffers(b) % reset_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1)
        out_array = in_array
      else if (buffers(b) % accumulate_flag) then
        call mpas_pool_get_array(block % allFields, &
            tvar % output_names(b), out_array, 1) 

!        out_array = (out_array * \
!                (buffers(b) % total_accum - 1) + in_array) \
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
