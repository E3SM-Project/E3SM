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
  public :: ocn_bootstrap_time_series_stats, &
         ocn_init_time_series_stats, &
         ocn_compute_time_series_stats, &
         ocn_restart_time_series_stats, &
         ocn_finalize_time_series_stats

  ! Private module variables
  !--------------------------------------------------------------------

  type time_series_alarms_type
    type (mpas_time_type) :: start_time
    type (mpas_timeinterval_type) :: duration_interval
    type (mpas_timeinterval_type) :: repeat_interval
    type (mpas_timeinterval_type) :: reset_interval
  end type time_series_alarms_type

  type time_series_variable_type
    ! state per variable, stored in framework
    character (len=StrKIND), pointer :: input_name
    character (len=StrKIND), dimension(:), allocatable :: output_names
  end type time_series_variable_type

  type time_series_buffer_type
    ! state per buffer, stored in framework
    integer, pointer :: started_flag, accumulate_flag, reset_flag

    ! strings for looking up alarms and buffers per buffer
    character (len=StrKIND), pointer :: start_alarm_ID, repeat_alarm_ID, &
      duration_alarm_ID, reset_alarm_ID

    ! counter for accumulation
    integer, pointer :: counter
    character (len=StrKIND), pointer :: xtime_start
    character (len=StrKIND), pointer :: xtime_end
  end type time_series_buffer_type

  type time_series_type
    ! state per instance, stored in framework
    integer, pointer :: operation
    integer, pointer :: number_of_variables
    integer, pointer :: number_of_buffers

    ! allocated on every instance call
    type (time_series_variable_type), dimension(:), allocatable :: variables
    type (time_series_buffer_type), dimension(:), allocatable :: buffers
  end type time_series_type

  ! enum of ops and types
  integer, parameter :: AVG_OP = 1
  integer, parameter :: MIN_OP = 2
  integer, parameter :: MAX_OP = 3

  integer, parameter :: START_TIMES = 11
  integer, parameter :: DURATION_INTERVALS = 12
  integer, parameter :: REPEAT_INTERVALS = 13
  integer, parameter :: RESET_INTERVALS = 14

  character (len=3), parameter :: AVG_TOKEN = 'avg'
  character (len=3), parameter :: MIN_TOKEN = 'min'
  character (len=3), parameter :: MAX_TOKEN = 'max'

  character (len=11), parameter :: TIME_START_PREFIX = 'xtime_start'
  character (len=9), parameter :: TIME_END_PREFIX = 'xtime_end'

  character (len=StrKIND), parameter :: TIME_SERIES_STATS_POOL = &
    'timeSeriesStatsAM'
  character (len=StrKIND), parameter :: ONE_STRING_MEMORY = &
    'timeSeriesStatsOneString'
  character (len=StrKIND), parameter :: ONE_INTEGER_MEMORY = &
    'timeSeriesStatsOneInteger'
  character (len=StrKIND), parameter :: ONE_REAL_MEMORY = &
    'timeSeriesStatsOneReal'

  character (len=StrKIND), parameter :: CONFIG_PREFIX = &
    'config_AM_timeSeriesStats'
  character (len=StrKIND), parameter :: FRAMEWORK_PREFIX = 'time'

  character (len=StrKIND), parameter :: OUTPUT_STREAM_SUFFIX = '_output_stream'
  character (len=StrKIND), parameter :: RESTART_STREAM_SUFFIX = '_restart_stream'
  character (len=StrKIND), parameter :: OPERATION_SUFFIX = '_operation'

  character (len=StrKIND), parameter :: NUMBER_OF_BUFFERS_SUFFIX = &
    '_number_of_buffers'
  character (len=StrKIND), parameter :: NUMBER_OF_VARIABLES_SUFFIX = &
    '_number_of_variables'

  character (len=StrKIND), parameter :: INPUT_NAME_SUFFIX = '_input_name'

  character (len=StrKIND), parameter :: REFERENCE_TIMES_SUFFIX = &
    '_reference_times'
  character (len=StrKIND), parameter :: DURATION_INTERVALS_SUFFIX = &
    '_duration_intervals'
  character (len=StrKIND), parameter :: REPEAT_INTERVALS_SUFFIX = &
    '_repeat_intervals'
  character (len=StrKIND), parameter :: RESET_INTERVALS_SUFFIX = &
    '_reset_intervals'

  character (len=StrKIND), parameter :: STARTED_FLAG_SUFFIX = &
    '_started_flag'
  character (len=StrKIND), parameter :: ACCUMULATE_FLAG_SUFFIX = &
    '_accumulate_flag'
  character (len=StrKIND), parameter :: RESET_FLAG_SUFFIX = &
    '_reset_flag'
  character (len=StrKIND), parameter :: START_ALARM_ID_SUFFIX = &
    '_start_alarm_ID'
  character (len=StrKIND), parameter :: REPEAT_ALARM_ID_SUFFIX = &
    '_repeat_alarm_ID'
  character (len=StrKIND), parameter :: DURATION_ALARM_ID_SUFFIX = &
    '_duration_alarm_ID'
  character (len=StrKIND), parameter :: RESET_ALARM_ID_SUFFIX = &
    '_reset_alarm_ID'
  character (len=StrKIND), parameter :: COUNTER_SUFFIX = &
    '_counter'

  character (len=StrKIND), parameter :: START_ALARM_PREFIX = '_startAlarm_'
  character (len=StrKIND), parameter :: REPEAT_ALARM_PREFIX = '_repeatAlarm_'
  character (len=StrKIND), parameter :: DURATION_ALARM_PREFIX = &
    '_durationAlarm_'
  character (len=StrKIND), parameter :: RESET_ALARM_PREFIX = '_resetAlarm_'

  character (len=StrKIND), parameter :: INITIAL_TIME_TOKEN = 'initial_time'
  character (len=StrKIND), parameter :: REPEAT_INTERVAL_TOKEN = &
    'repeat_interval'
  character (len=StrKIND), parameter :: RESET_INTERVAL_TOKEN = 'reset_interval'

  character (len=StrKIND), parameter :: CURRENT_CORE_NAME = 'MPAS-Ocean'
  character (len=4), parameter :: NONE_TOKEN = 'none'

!***********************************************************************
contains

!***********************************************************************
! routine ocn_bootstrap_time_series_stats
!
!> \brief Bootstrap time_series_stats analysis member
!> \author  Doug Jacobsen
!> \date    10/08/2015
!> \details
!>  This routine performs pre-init configuration of the analysis member.
!>  Specifically, it ensures the streams used for this instance are correctly
!>  configured.
!-----------------------------------------------------------------------
subroutine ocn_bootstrap_time_series_stats(domain, instance, err)!{{{
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: v
  type (time_series_type) :: series
  logical, dimension(:), pointer :: valid_input

  ! start procedure
  err = 0

  ! initial allocation of instance state for this AM from the namelist
  call start_state(domain, instance, series, valid_input, err)

  ! modify the output and restart streams for this AM instance
  ! driver will do a restart read, after this, if necessary to fill values
  call modify_stream(domain, instance, series, valid_input, err)

  ! clean up the instance memory
  do v = 1, series % number_of_variables
    deallocate(series % variables(v) % output_names)
  end do
  deallocate(series % variables)
  deallocate(series % buffers)
end subroutine ocn_bootstrap_time_series_stats!}}}

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
subroutine ocn_init_time_series_stats(domain, instance, err)!{{{
  ! input variables
  character (len=StrKIND), intent(in) :: instance 

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: v, b
  type (time_series_type) :: series
  type (time_series_alarms_type), allocatable, dimension(:) :: alarms
  type (MPAS_Time_type) :: start_intv
  character (len=StrKIND) :: start_xtime

  ! start procedure
  err = 0

  ! coming back from a potential restart read
  ! get all of the state for this instance
  call get_state(domain, instance, series)

  ! get all of the timing configurations from namelist
  allocate(alarms(series % number_of_buffers))
  call get_alarms(domain, instance, series, alarms, err)

  ! set the values of the alarms and current flag states based on timers
  call set_alarms(domain, instance, series, alarms, err)
  deallocate(alarms)

  ! set xtime start if it is still unset, for very first time step
  ! (i.e., no restarts)
  start_intv = mpas_get_clock_time(domain % clock, MPAS_NOW, err)
  call mpas_get_time(start_intv, dateTimeString=start_xtime, ierr=err)

  do b = 1, series % number_of_buffers
    if (trim(series % buffers(b) % xtime_start) == '') then
      series % buffers(b) % xtime_start = start_xtime
    end if
  end do

  ! clean up the instance memory
  do v = 1, series % number_of_variables
    deallocate(series % variables(v) % output_names)
  end do
  deallocate(series % variables)
  deallocate(series % buffers)
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
subroutine ocn_compute_time_series_stats(domain, timeLevel, instance, err)!{{{
  ! input variables
  character (len=StrKIND), intent(in) :: instance
  integer, intent(in) :: timeLevel

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: v, b
  type (time_series_type) :: series
  type (MPAS_TimeInterval_type) :: dt
  type (MPAS_Time_type) :: start_intv, end_intv
  character (len=StrKIND) :: start_xtime, end_xtime
  logical :: unset_xtime

  ! start procedure
  err = 0

  ! get all of the state for this instance to be able to compute
  call get_state(domain, instance, series)

  ! get the strings for the date
  unset_xtime = .true.

  ! update the counter
  do b = 1, series % number_of_buffers
    if (series % buffers(b) % accumulate_flag == 1) then
      if (unset_xtime) then
        end_intv = mpas_get_clock_time(domain % clock, MPAS_NOW, err)
        call mpas_get_time(end_intv, dateTimeString=end_xtime, ierr=err)
        start_intv = end_intv - mpas_get_clock_timestep(domain % clock, err)
        call mpas_get_time(start_intv, dateTimeString=start_xtime, ierr=err)
        unset_xtime = .false.
      end if

      if (series % buffers(b) % reset_flag == 1) then
        series % buffers(b) % xtime_start = start_xtime
        series % buffers(b) % counter = 1
      else
        series % buffers(b) % xtime_end = end_xtime
        series % buffers(b) % counter = series % buffers(b) % counter + 1
      end if
    end if
  end do

  ! do all of the operations
  do v = 1, series % number_of_variables
    call typed_operate(domain % blocklist, &
      series % variables(v), &
      series % buffers, &
      series % operation)
  end do

  ! do all of the time checking and flag setting
  call timer_checking(series, domain % clock, err)

  ! clean up the instance memory
  do v = 1, series % number_of_variables
    deallocate(series % variables(v) % output_names)
  end do
  deallocate(series % variables)
  deallocate(series % buffers)
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
subroutine ocn_restart_time_series_stats(domain, instance, err)!{{{
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables

  ! start procedure
  err = 0

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
subroutine ocn_finalize_time_series_stats(domain, instance, err)!{{{
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables

  ! start procedure
  err = 0

end subroutine ocn_finalize_time_series_stats!}}}

!
! local subroutines
!

!***********************************************************************
! routine debug_state
!
!> \brief Print the internal state of the time series analysis member
!> \author  Jon Woodring
!> \date    May 9, 2016
!> \details
!> Print the state for the time series analysis member. Primarily for
!> internal debugging.
!-----------------------------------------------------------------------
subroutine debug_state(series)
  ! input variables
  type (time_series_type), intent(in) :: series

  ! input/output variables

  ! output variables

  ! local variables
  integer :: v, b

  ! start procedure
  write(*,*) 'operation: ', series % operation
  write(*,*) '# variables: ', series % number_of_variables
  write(*,*) '# buffers: ', series % number_of_buffers

  do b = 1, series % number_of_buffers
    write(*,*) 'buffer #: ', b
    write(*,*) 'started: ', series % buffers(b) % started_flag
    write(*,*) 'accum: ', series % buffers(b) % accumulate_flag
    write(*,*) 'reset: ', series % buffers(b) % reset_flag
    write(*,*) 'start alarm: ', trim(series % buffers(b) % start_alarm_ID)
    write(*,*) 'repeat alarm: ', trim(series % buffers(b) % repeat_alarm_ID)
    write(*,*) 'duration alarm: ', trim(series % buffers(b) % duration_alarm_ID)
    write(*,*) 'reset alarm: ', trim(series % buffers(b) % reset_alarm_ID)
    write(*,*) 'counter: ', series % buffers(b) % counter
    write(*,*) 'xtime start: ', trim(series % buffers(b) % xtime_start)
    write(*,*) 'xtime end: ', trim(series % buffers(b) % xtime_end)
  end do

  do v = 1, series % number_of_variables
    write(*,*) 'variable #: ', v
    write(*,*) 'input: ', trim(series % variables(v) % input_name)
    do b = 1, series % number_of_buffers
      write(*,*) 'buffer #: ', b
      write(*,*) 'output: ', trim(series % variables(v) % output_names(b))
    end do
  end do

end subroutine debug_state

!***********************************************************************
! routine add_new_string
!
!> \brief Allocate a string in the MPAS framework for this AM
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Allocate a new integer in the AM pool and return a pointer to it.
!-----------------------------------------------------------------------
subroutine add_new_string(all_fields, inpool, outpool, field_name, target_ptr)
  ! input variables
  character (len=StrKIND) :: field_name

  ! input/output variables
  type (mpas_pool_type), pointer, intent(inout) :: all_fields, inpool, outpool

  ! output variables
  character (len=StrKIND), pointer, optional :: target_ptr

  ! local variables
  type (field0DChar), pointer :: srcString, dstString

  call mpas_pool_get_field(inpool, ONE_STRING_MEMORY, srcString, 1)
  call mpas_duplicate_field(srcString, dstString)
  dstString % fieldName = field_name
  call mpas_pool_add_field(outpool, dstString % fieldName, dstString)
  call mpas_pool_add_field(all_fields, dstString % fieldName, dstString)
  if (present(target_ptr)) then
    call mpas_pool_get_array(outpool, dstString % fieldName, target_ptr, 1)
  end if
end subroutine add_new_string



!***********************************************************************
! routine add_new_integer
!
!> \brief Allocate an integer in the MPAS framework for this AM
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Allocate a new integer in the AM pool and return a pointer to it.
!-----------------------------------------------------------------------
subroutine add_new_integer(all_fields, inpool, outpool, field_name, target_ptr)
  ! input variables
  character (len=StrKIND) :: field_name

  ! input/output variables
  type (mpas_pool_type), pointer, intent(inout) :: all_fields, inpool, outpool

  ! output variables
  integer, pointer, optional :: target_ptr

  ! local variables
  type (field0DInteger), pointer :: srcInteger, dstInteger

  call mpas_pool_get_field(inpool, ONE_INTEGER_MEMORY, srcInteger, 1)
  call mpas_duplicate_field(srcInteger, dstInteger)
  dstInteger % fieldName = field_name
  call mpas_pool_add_field(outpool, dstInteger % fieldName, dstInteger)
  call mpas_pool_add_field(all_fields, dstInteger % fieldName, dstInteger)
  if (present(target_ptr)) then
    call mpas_pool_get_array(outpool, dstInteger % fieldName, target_ptr, 1)
  end if
end subroutine add_new_integer



!***********************************************************************
! routine add_new_real
!
!> \brief Allocate a real in the MPAS framework for this AM
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Allocate a new real in the AM pool and return a pointer to it.
!-----------------------------------------------------------------------
subroutine add_new_real(all_fields, inpool, outpool, field_name, target_ptr)
  ! input variables
  character (len=StrKIND) :: field_name

  ! input/output variables
  type (mpas_pool_type), pointer, intent(inout) :: all_fields, inpool, outpool

  ! output variables
  real (kind=RKIND), pointer, optional :: target_ptr

  ! local variables
  type (field0DReal), pointer :: srcReal, dstReal

  call mpas_pool_get_field(inpool, ONE_REAL_MEMORY, srcReal, 1)
  call mpas_duplicate_field(srcReal, dstReal)
  dstReal % fieldName = field_name
  call mpas_pool_add_field(outpool, dstReal % fieldName, dstReal)
  call mpas_pool_add_field(all_fields, dstReal % fieldName, dstReal)
  if (present(target_ptr)) then
    call mpas_pool_get_array(outpool, dstReal % fieldName, target_ptr, 1)
  end if
end subroutine add_new_real



!***********************************************************************
! routine check_real_time
!
!> \brief Check to see if an MPAS field can have time series stats applied
!> \author  Jon Woodring
!> \date    December 17, 2015
!> \details
!> Makes sure that time series stats can be applied to a field_name.
!> In particular, checks to see if it is real and has time levels.
!-----------------------------------------------------------------------
logical function check_real_time(all_fields, field_name)
  type (mpas_pool_type), pointer, intent(in) :: all_fields
  character (len=StrKIND), intent(in) :: field_name

  logical :: check
  type (field0DReal), pointer :: r0
  type (field1DReal), pointer :: r1
  type (field2DReal), pointer :: r2
  type (field3DReal), pointer :: r3
  type (field4DReal), pointer :: r4
  type (field5DReal), pointer :: r5
  type(mpas_pool_field_info_type) :: info

  ! get the info of the field
  call mpas_pool_get_field_info(all_fields, field_name, info)
  check = info % fieldType == MPAS_POOL_REAL

  if(.not. check) then
    write(stdErrUnit,*) trim(CURRENT_CORE_NAME) // &
      ' WARNING: field "' // trim(field_name) // '" listed in the ' // &
      'output stream, for time series stats analysis member ' // &
      'stream, is not real. Time series stats will not be applied to ' // &
      'this field.'
  else
    if(info % nDims == 0) then
      call mpas_pool_get_field(all_fields, field_name, r0, 1)
      check = r0 % hasTimeDimension
    else if(info % nDims == 1) then
      call mpas_pool_get_field(all_fields, field_name, r1, 1)
      check = r1 % hasTimeDimension
    else if(info % nDims == 2) then
      call mpas_pool_get_field(all_fields, field_name, r2, 1)
      check = r2 % hasTimeDimension
    else if(info % nDims == 3) then
      call mpas_pool_get_field(all_fields, field_name, r3, 1)
      check = r3 % hasTimeDimension
    else if(info % nDims == 4) then
      call mpas_pool_get_field(all_fields, field_name, r4, 1)
      check = r4 % hasTimeDimension
    else if(info % nDims == 5) then
      call mpas_pool_get_field(all_fields, field_name, r5, 1)
      check = r5 % hasTimeDimension
    end if

    if (.not. check) then
      write(stdErrUnit,*) trim(CURRENT_CORE_NAME) // &
        ' WARNING: field "' // trim(field_name) // '" listed in the ' // &
        'output stream, for time series stats analysis member ' // &
        'stream, does not have a time dimension. Time series stats will ' // &
        'not be applied to this field.'
    end if
  end if

  check_real_time = check

end function check_real_time


!***********************************************************************
! routine get_state
!
!> \brief Get all of the state for this instance.
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details
!> This will allocate and fetch all of the state necessary for this
!> instance that is being run.
!-----------------------------------------------------------------------
subroutine get_state(domain, instance, series)
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  type (time_series_type), intent(out) :: series

  ! local variables
  integer :: v, b
  character (len=StrKIND) :: storage_prefix, var_identifier, &
    buf_identifier, var_prefix, buf_prefix, field_name, op_name
  type (mpas_pool_type), pointer :: amPool

  ! start procedure
  storage_prefix = trim(FRAMEWORK_PREFIX) // trim(instance)
  call mpas_pool_get_subpool(domain % blocklist % structs, &
    TIME_SERIES_STATS_POOL, amPool)

  !
  ! get the base
  !

  ! number_of_variables
  field_name = trim(storage_prefix) // trim(NUMBER_OF_VARIABLES_SUFFIX)
  call mpas_pool_get_array(amPool, field_name, series % number_of_variables, 1)

  ! number_of_buffers
  field_name = trim(storage_prefix) // trim(NUMBER_OF_BUFFERS_SUFFIX)
  call mpas_pool_get_array(amPool, field_name, series % number_of_buffers, 1)

  ! operation
  field_name = trim(storage_prefix) // trim(OPERATION_SUFFIX)
  call mpas_pool_get_array(amPool, field_name, series % operation, 1)

  op_name = operator_naming(series % operation)

  ! create the memory
  allocate(series % variables(series % number_of_variables))
  allocate(series % buffers(series % number_of_buffers))
  do v = 1, series % number_of_variables
    allocate(series % variables(v) % output_names(series % number_of_buffers))
  end do

  !
  ! get the instance values for variables
  !

  do v = 1, series % number_of_variables
    ! identifier
    write(var_identifier, '(I0)') v
    var_prefix = trim(storage_prefix) // '_' // trim(var_identifier)

    ! input_name
    field_name = trim(var_prefix) // trim(INPUT_NAME_SUFFIX)
    call mpas_pool_get_array(amPool, &
      field_name, series % variables(v) % input_name, 1)

    if (series % number_of_buffers > 1) then
      do b = 1, series % number_of_buffers
        write(buf_identifier, '(I0)') b

        ! create output names
        series % variables(v) % output_names(b) = output_naming &
          (storage_prefix, op_name, series % variables(v) % input_name, &
           buf_identifier)
      end do
    else
      series % variables(v) % output_names(1) = output_naming &
        (storage_prefix, op_name, series % variables(v) % input_name)
    end if
  end do

  !
  ! get the instance values for buffers
  !

  do b = 1, series % number_of_buffers
    ! identifier
    write(buf_identifier, '(I0)') b
    buf_prefix = trim(storage_prefix) // '_' // trim(buf_identifier)

    ! started_flag
    field_name = trim(buf_prefix) // trim(STARTED_FLAG_SUFFIX)
    call mpas_pool_get_array(amPool, &
      field_name, series % buffers(b) % started_flag, 1)

    ! accumulate_flag
    field_name = trim(buf_prefix) // trim(ACCUMULATE_FLAG_SUFFIX)
    call mpas_pool_get_array(amPool, &
      field_name, series % buffers(b) % accumulate_flag, 1)

    ! reset_flag
    field_name = trim(buf_prefix) // trim(RESET_FLAG_SUFFIX)
    call mpas_pool_get_array(amPool, &
      field_name, series % buffers(b) % reset_flag, 1)

    ! start_alarm_ID
    field_name = trim(buf_prefix) // trim(START_ALARM_ID_SUFFIX)
    call mpas_pool_get_array(amPool, &
      field_name, series % buffers(b) % start_alarm_ID, 1)

    ! repeat_alarm_ID
    field_name = trim(buf_prefix) // trim(REPEAT_ALARM_ID_SUFFIX)
    call mpas_pool_get_array(amPool, &
      field_name, series % buffers(b) % repeat_alarm_ID, 1)

    ! duration_alarm_ID
    field_name = trim(buf_prefix) // trim(DURATION_ALARM_ID_SUFFIX)
    call mpas_pool_get_array(amPool, &
      field_name, series % buffers(b) % duration_alarm_ID, 1)

    ! reset_alarm_ID
    field_name = trim(buf_prefix) // trim(RESET_ALARM_ID_SUFFIX)
    call mpas_pool_get_array(amPool, &
      field_name, series % buffers(b) % reset_alarm_ID, 1)

    ! counter
    if (series % number_of_buffers > 1) then
      field_name = counter_naming(storage_prefix, buf_identifier)
    else
      field_name = counter_naming(storage_prefix)
    end if
    call mpas_pool_get_array(amPool, &
      field_name, series % buffers(b) % counter, 1)

    ! xtime start
    if (series % number_of_buffers > 1) then
      field_name = trim(TIME_START_PREFIX) // trim(instance) // &
        '_' // buf_identifier
    else
      field_name = trim(TIME_START_PREFIX) // trim(instance)
    end if
    call mpas_pool_get_array(amPool, &
      field_name, series % buffers(b) % xtime_start, 1)

    ! xtime end
    if (series % number_of_buffers > 1) then
      field_name = trim(TIME_END_PREFIX) // trim(instance) // &
        '_' // buf_identifier
    else
      field_name = trim(TIME_END_PREFIX) // trim(instance)
    end if
    call mpas_pool_get_array(amPool, &
      field_name, series % buffers(b) % xtime_end, 1)

  end do

end subroutine get_state

!***********************************************************************
! routine start_state
!
!> \brief Begin the initialization of this analysis member
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details
!> This will count the number of variables, number of buffers, and
!> also get the stream name and operation strings.
!-----------------------------------------------------------------------
subroutine start_state(domain, instance, series, valid_input, err)
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain
  logical, dimension(:), pointer :: valid_input

  ! output variables
  type (time_series_type), intent(out) :: series
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  character (len=StrKIND), pointer :: config_results, output_stream_name
  character (len=StrKIND) :: config, namelist_prefix, storage_prefix, &
    var_identifier, buf_identifier, var_prefix, buf_prefix, field_name
  integer :: b, v
  type (mpas_pool_type), pointer :: amPool

  ! start procedure
  err = 0

  namelist_prefix = trim(CONFIG_PREFIX) // trim(instance)
  storage_prefix = trim(FRAMEWORK_PREFIX) // trim(instance)
  call mpas_pool_get_subpool(domain % blocklist % structs, &
    TIME_SERIES_STATS_POOL, amPool)

  !
  ! allocate some framework memory for instance state
  !

  ! number_of_variables
  field_name = trim(storage_prefix) // trim(NUMBER_OF_VARIABLES_SUFFIX)
  call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
    field_name, series % number_of_variables)

  ! number_of_buffers
  field_name = trim(storage_prefix) // trim(NUMBER_OF_BUFFERS_SUFFIX)
  call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
    field_name, series % number_of_buffers)

  ! operation
  field_name = trim(storage_prefix) // trim(OPERATION_SUFFIX)
  call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
    field_name, series % operation)

  !
  ! assign some instance values
  !

  ! get the stream name
  config = trim(namelist_prefix) // trim(OUTPUT_STREAM_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, output_stream_name)

  if (output_stream_name == NONE_TOKEN) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: stream cannot be "none" ' // &
      'for time series stats.')
  end if

  ! count the number of variables
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)
  b = 0
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    output_stream_name, field_name))
    b = b + 1
  end do

  allocate(valid_input(b))

  ! count the number of variables and mark if valid
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)
  b = 1
  series % number_of_variables = 0
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    output_stream_name, field_name))

    valid_input(b) = check_real_time(domain % blocklist % allFields, &
      field_name)

    if (valid_input(b)) then
      series % number_of_variables = series % number_of_variables + 1
    end if

    b = b + 1
  end do

  if (series % number_of_variables < 1) then
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: there are no fields ' // &
      'in the time series stats output stream "' // trim(output_stream_name) &
      // '" that time series stats can be applied to.')
  end if

  ! count the number of buffers
  config = trim(namelist_prefix) // trim(REFERENCE_TIMES_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, config_results)
  config = config_results
  series % number_of_buffers = 1
  b = scan(config, ';')
  do while (b > 0)
    series % number_of_buffers = series % number_of_buffers + 1
    config = config(b+1:)
    b = scan(config, ';')
  end do

  ! get our operation
  config = trim(namelist_prefix) // trim(OPERATION_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, config_results)

  series % operation = operator_enum(config_results)

  ! create the memory
  allocate(series % variables(series % number_of_variables))
  allocate(series % buffers(series % number_of_buffers))
  do v = 1, series % number_of_variables
    allocate(series % variables(v) % output_names(series % number_of_buffers))
  end do

  !
  ! duplicate memory for storing AM instance state in the framework
  !

  ! create variable space
  b = 1
  do v = 1, series % number_of_variables
    ! identifier
    write(var_identifier, '(I0)') v
    var_prefix = trim(storage_prefix) // '_' // trim(var_identifier)

    ! input_name
    field_name = trim(var_prefix) // trim(INPUT_NAME_SUFFIX)
    call add_new_string(domain % blocklist % allFields, amPool, amPool, &
      field_name, series % variables(v) % input_name)
  end do

  ! create buffer space
  do b = 1, series % number_of_buffers
    ! identifier
    write(buf_identifier, '(I0)') b
    buf_prefix = trim(storage_prefix) // '_' // trim(buf_identifier)

    ! started_flag
    field_name = trim(buf_prefix) // trim(STARTED_FLAG_SUFFIX)
    call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
      field_name, series % buffers(b) % started_flag)

    ! accumulate_flag
    field_name = trim(buf_prefix) // trim(ACCUMULATE_FLAG_SUFFIX)
    call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
      field_name, series % buffers(b) % accumulate_flag)

    ! reset_flag
    field_name = trim(buf_prefix) // trim(RESET_FLAG_SUFFIX)
    call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
      field_name, series % buffers(b) % reset_flag)

    ! start_alarm_ID
    field_name = trim(buf_prefix) // trim(START_ALARM_ID_SUFFIX)
    call add_new_string(domain % blocklist % allFields, amPool, amPool, &
      field_name, series % buffers(b) % start_alarm_ID)

    ! repeat_alarm_ID
    field_name = trim(buf_prefix) // trim(REPEAT_ALARM_ID_SUFFIX)
    call add_new_string(domain % blocklist % allFields, amPool, amPool, &
      field_name, series % buffers(b) % repeat_alarm_ID)

    ! duration_alarm_ID
    field_name = trim(buf_prefix) // trim(DURATION_ALARM_ID_SUFFIX)
    call add_new_string(domain % blocklist % allFields, amPool, amPool, &
      field_name, series % buffers(b) % duration_alarm_ID)

    ! reset_alarm_ID
    field_name = trim(buf_prefix) // trim(RESET_ALARM_ID_SUFFIX)
    call add_new_string(domain % blocklist % allFields, amPool, amPool, &
      field_name, series % buffers(b) % reset_alarm_ID)

    !
    ! counter, xtime_start, and xtime_end is not allocated here,
    ! because it is part of the restart stream and output_stream,
    ! and not just the internal AM state
    !
    ! it is allocated in modify_stream
    !
  end do
end subroutine start_state



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
subroutine modify_stream(domain, instance, series, valid_input, err)!{{{
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain
  type (time_series_type), intent(inout) :: series
  logical, dimension(:), pointer :: valid_input

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: v, b
  logical :: emptyRestartStream
  character (len=StrKIND), pointer :: output_stream_name, restart_stream_name
  character (len=StrKIND) :: fieldName
  character (len=StrKIND) :: field_name, config, op_name
  character (len=StrKIND) :: namelist_prefix, &
    storage_prefix, buf_identifier, buf_prefix
  type (mpas_pool_field_info_type) :: info
  type (mpas_pool_type), pointer :: amPool

  ! start procedure
  err = 0

  namelist_prefix = trim(CONFIG_PREFIX) // trim(instance)
  storage_prefix = trim(FRAMEWORK_PREFIX) // trim(instance)
  call mpas_pool_get_subpool(domain % blocklist % structs, &
    TIME_SERIES_STATS_POOL, amPool)

  ! get the output stream name
  config = trim(namelist_prefix) // trim(OUTPUT_STREAM_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, output_stream_name)

  ! get restart stream name
  config = trim(namelist_prefix) // trim(RESTART_STREAM_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, restart_stream_name)

  op_name = operator_naming(series % operation)

  ! get the old field names, assign to input name, and remove from stream
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)
  b = 1
  v = 1
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    output_stream_name, field_name))

    ! check if we can handle it
    if (valid_input(b)) then
      series % variables(v) % input_name = field_name

      ! remove the old one
      call mpas_stream_mgr_remove_field(domain % streamManager, &
        output_stream_name, series % variables(v) % input_name)

      v = v + 1
    end if

    b = b + 1
  end do

  deallocate(valid_input)


  !
  ! create memory and modify the stream
  !

  ! ensure restart stream is empty
  emptyRestartStream = .true.
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    streamID=restart_stream_name, ierr=err)
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    streamID=restart_stream_name, fieldName=fieldName) .and. emptyRestartStream)
    emptyRestartStream = .false.
  end do

  if (.not. emptyRestartStream) then
    write(stderrUnit, *) trim(CURRENT_CORE_NAME) // ' ERROR: stream named ''' &
      // trim(restart_stream_name) // ''' is not empty, but is used in ' // &
      'an restart instance of the time series stats analysis member. ' // &
      'This restart stream will be built based on the contents of the ''' // &
      trim(output_stream_name) // ''' stream. Please ensure it is empty ' // &
      'in the streams input file.'
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // &
      'ERROR: misconfigured restart stream for time series stats ' // &
      ' analysis member.')
  end if

  ! create and put the counter, xtime_start, xtime_end in the streams
  do b = 1, series % number_of_buffers
    write(buf_identifier, '(I0)') b

    ! allocate counter memory
    if (series % number_of_buffers > 1) then
      field_name = counter_naming(storage_prefix, buf_identifier)
    else
      field_name = counter_naming(storage_prefix)
    end if
    call add_new_integer(domain % blocklist % allFields, amPool, amPool, &
      field_name, series % buffers(b) % counter)

    ! put it in the output and restart
    call mpas_stream_mgr_add_field(domain % streamManager, &
      output_stream_name, field_name, ierr=err)
    call mpas_stream_mgr_add_field(domain % streamManager, &
      restart_stream_name, field_name, ierr=err)

    ! xtime start
    if (series % number_of_buffers > 1) then
      field_name = trim(TIME_START_PREFIX) // trim(instance) // &
        '_' // buf_identifier
    else
      field_name = trim(TIME_START_PREFIX) // trim(instance)
    end if
    call add_new_string(domain % blocklist % allFields, amPool, amPool, &
      field_name, series % buffers(b) % xtime_start)

    ! put it in the output and restart
    call mpas_stream_mgr_add_field(domain % streamManager, &
      output_stream_name, field_name, ierr=err)
    call mpas_stream_mgr_add_field(domain % streamManager, &
      restart_stream_name, field_name, ierr=err)

    ! xtime end
    if (series % number_of_buffers > 1) then
      field_name = trim(TIME_END_PREFIX) // trim(instance) // &
        '_' // buf_identifier
    else
      field_name = trim(TIME_END_PREFIX) // trim(instance)
    end if
    call add_new_string(domain % blocklist % allFields, amPool, amPool, &
      field_name, series % buffers(b) % xtime_end)

    ! put it in the output and restart
    call mpas_stream_mgr_add_field(domain % streamManager, &
      output_stream_name, field_name, ierr=err)
    call mpas_stream_mgr_add_field(domain % streamManager, &
      restart_stream_name, field_name, ierr=err)
  end do

  ! set up the variables
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)
  do v = 1, series % number_of_variables
    ! get the info of the field
    call mpas_pool_get_field_info(domain % blocklist % allFields, &
      series % variables(v) % input_name, info)

    ! allocate a number of fields and add field
    do b = 1, series % number_of_buffers
      write(buf_identifier, '(I0)') b

      if (series % number_of_buffers > 1) then
        field_name = output_naming(storage_prefix, op_name, &
          series % variables(v) % input_name, buf_identifier)
      else
        field_name = output_naming(storage_prefix, op_name, &
          series % variables(v) % input_name)
      end if

      ! create the name of the output var
      series % variables(v) % output_names(b) = field_name

      ! create the field and add to pool
      call add_new_field(info, &
        series % variables(v) % input_name, &
        series % variables(v) % output_names(b), &
        domain % blocklist % allFields, amPool)

      ! add the field to the output stream
      call mpas_stream_mgr_add_field(domain % streamManager, &
        output_stream_name, series % variables(v) % output_names(b), ierr=err)

      ! put it in the restart stream
      call mpas_stream_mgr_add_field(domain % streamManager, &
        restart_stream_name, series % variables(v) % output_names(b), ierr=err)
    end do
  end do ! number_of_variables

end subroutine modify_stream!}}}


!***********************************************************************
! function output_naming
!
!> \brief Given an input name, create a cooresponding output name
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details
!> Code to create consistent output names from input names.
!-----------------------------------------------------------------------
character (len=StrKIND) function output_naming &
(storage_prefix, op_name, input_name, buf_identifier)
  character (len=StrKIND), intent(in) :: storage_prefix, op_name, input_name
  character (len=StrKIND), intent(in), optional :: buf_identifier

  if (present(buf_identifier)) then
    output_naming = trim(storage_prefix) // '_' // trim(op_name) // '_' // &
      trim(input_name) // '_' // trim(buf_identifier)
  else
    output_naming = trim(storage_prefix) // '_' // trim(op_name) // '_' // &
      trim(input_name)
  endif
end function output_naming


!***********************************************************************
! function counter_naming
!
!> \brief Given an buffer number, create a cooresponding counter name
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details
!> Code to create consistent counter names from buffer numbers.
!-----------------------------------------------------------------------
character (len=StrKIND) function counter_naming &
(storage_prefix, buf_identifier)
  character (len=StrKIND), intent(in) :: storage_prefix
  character (len=StrKIND), intent(in), optional :: buf_identifier

  if (present(buf_identifier)) then
    counter_naming = trim(storage_prefix) // trim(COUNTER_SUFFIX) // &
      '_' // trim(buf_identifier)
  else
    counter_naming = trim(storage_prefix) // trim(COUNTER_SUFFIX)
  end if
end function counter_naming


!***********************************************************************
! function operator_naming
!
!> \brief Given an operator number, create a corresponding operator name
!> \author  Jon Woodring
!> \date    May 9, 2016
!> \details
!> Code to create consistent operator names from operator enums.
!-----------------------------------------------------------------------
character (len=StrKIND) function operator_naming (operator_enum)
  integer, intent(in) :: operator_enum

  ! operator
  if (operator_enum == AVG_OP) then
    operator_naming = AVG_TOKEN
  else if (operator_enum == MIN_OP) then
    operator_naming = MIN_TOKEN
  else if (operator_enum == MAX_OP) then
    operator_naming = MAX_TOKEN
  else
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: the impossible happened - ' // &
      'tried to create an operation in the time series stats ' // &
      'analysis member of unknown kind')
  end if
end function operator_naming


!***********************************************************************
! function operator_enum
!
!> \brief Given an operator token, create a corresponding operator enum
!> \author  Jon Woodring
!> \date    May 9, 2016
!> \details
!> Code to create consistent operator enums from operator tokens.
!-----------------------------------------------------------------------
integer function operator_enum (token)
  character (len=StrKIND), intent(in) :: token

  if (trim(token) == AVG_TOKEN) then
    operator_enum = AVG_OP
  else if (trim(token) == MIN_TOKEN) then
    operator_enum = MIN_OP
  else if (trim(token) == MAX_TOKEN) then
    operator_enum = MAX_OP
  else
    call mpas_dmpar_global_abort( &
      trim(CURRENT_CORE_NAME) // ' ERROR: unknown operation "' // &
      trim(token) // '" requested in the time ' // &
      'series stats analysis member configuration')
  end if
end function operator_enum


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
subroutine get_alarms(domain, instance, series, alarms, err)
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain
  type (time_series_type), intent(inout) :: series

  ! output variables
  integer, intent(out) :: err !< Output: error flag
  type (time_series_alarms_type), dimension(:), intent(out) :: alarms

  ! local variables
  character (len=StrKIND), pointer :: config_results
  character (len=StrKIND) :: config, namelist_prefix
  integer :: b, n
  logical :: ok
  type (mpas_timeinterval_type) :: rem, zero

  ! create prefix
  namelist_prefix = trim(CONFIG_PREFIX) // trim(instance)

  ! configure start times - we don't have to check ok
  ! because the timer count is based on reference_times tokens
  config = trim(namelist_prefix) // trim(REFERENCE_TIMES_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, config_results)
  call set_times(series, alarms, domain % clock, START_TIMES, &
    config_results, ok, err)

  ! order matters, don't reorder these following ones!
  ! it matters because times/intervals can be configured to be equal
  ! to other ones

  ! configure reset intervals
  config = trim(namelist_prefix) // trim(RESET_INTERVALS_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, config_results)
  call set_times(series, alarms, domain % clock, RESET_INTERVALS, &
    config_results, ok, err)
  if (.not. ok) then
    call mpas_dmpar_global_abort(trim(CURRENT_CORE_NAME) // &
      ' ERROR: number of listed times in ' // &
      'reset_intervals is not consistent with number of listed times ' // &
      'in reference_times in time series stats analysis member ' // &
      'configuration.')
  end if

  ! configure repeat intervals
  config = trim(namelist_prefix) // trim(REPEAT_INTERVALS_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, config_results)
  call set_times(series, alarms, domain % clock, REPEAT_INTERVALS, &
    config_results, ok, err)
  if (.not. ok) then
    call mpas_dmpar_global_abort(trim(CURRENT_CORE_NAME) // &
      ' ERROR: number of listed times in ' // &
      'repeat_intervals is not consistent with number of listed times ' // &
      'in reference_times in time series stats analysis member ' // &
      'configuration.')
  end if

  ! configure duration intervals
  config = trim(namelist_prefix) // trim(DURATION_INTERVALS_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, config_results)
  call set_times(series, alarms, domain % clock, DURATION_INTERVALS, &
    config_results, ok, err)
  if (.not. ok) then
    call mpas_dmpar_global_abort(trim(CURRENT_CORE_NAME) // &
      ' ERROR: number of listed times in ' // &
      'duration_intervals is not consistent with number of listed times ' // &
      'in reference_times in time series stats analysis member ' // &
      'configuration.')
  end if

  ! check if some of the time configuration is sensible
  call mpas_set_timeInterval(zero, s=0)

  do b = 1, series % number_of_buffers
    call mpas_interval_division(alarms(b) % start_time, &
       alarms(b) % repeat_interval, &
       alarms(b) % reset_interval, n, rem)

    if (n > 1 .or. (n == 1 .and. rem /= zero)) then
      write(stderrUnit,*) trim(CURRENT_CORE_NAME) // &
        'WARNING: repeat_interval > ' // &
        'reset_interval in time series stats analysis member ' // &
        'configuration. Truncating repeat_interval.'
      alarms(b) % repeat_interval = alarms(b) % reset_interval
    end if

    call mpas_interval_division(alarms(b) % start_time, &
       alarms(b) % duration_interval, &
       alarms(b) % repeat_interval, n, rem)

    if (n > 1 .or. (n == 1 .and. rem /= zero)) then
      write(stderrUnit,*) trim(CURRENT_CORE_NAME) // &
        'WARNING: duration_interval > ' // &
        'repeat_interval in time series stats analysis member ' // &
        'configuration. Truncating duration_interval.'
      alarms(b) % repeat_interval = alarms(b) % reset_interval
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
subroutine set_alarms(domain, instance, series, alarms, err)
  ! input variables
  character (len=StrKIND), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain
  type (time_series_type), intent(inout) :: series
  type (time_series_alarms_type), dimension(:), intent(inout) :: alarms

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: b, repeat_n, duration_n, reset_n
  character (len=StrKIND) :: buf_identifier, alarm_prefix
  type (mpas_time_type) :: current_time, when, &
    duration_time, repeat_time, reset_time
  type (mpas_timeinterval_type) :: elapsed, zero, &
    repeat_rem, duration_rem, reset_rem, zero_intv

  ! start procedure
  alarm_prefix = trim(FRAMEWORK_PREFIX) // trim(instance)

  ! get current time
  current_time = mpas_get_clock_time(domain % clock, MPAS_NOW, err)
  call mpas_set_timeInterval(zero_intv, S=0)

  ! configure alarms
  do b = 1, series % number_of_buffers
    write(buf_identifier, '(I0)') b

    ! zero flags
    series % buffers(b) % started_flag = 0
    series % buffers(b) % reset_flag = 0
    series % buffers(b) % accumulate_flag = 0

    ! set start time and flag
    if (current_time >= alarms(b) % start_time) then
       series % buffers(b) % started_flag = 1

       ! no start alarm
       series % buffers(b) % start_alarm_ID = ''
    else
       ! set the start alarm
       series % buffers(b) % start_alarm_ID = trim(alarm_prefix) // &
         trim(START_ALARM_PREFIX) // trim(buf_identifier)
       call mpas_add_clock_alarm(domain % clock, &
         series % buffers(b) % start_alarm_ID, &
         alarms(b) % start_time, ierr=err)
     end if

    ! set next reset time and flag
    when = alarms(b) % start_time + alarms(b) % reset_interval
    if (current_time >= when) then
      elapsed = current_time - when
      call mpas_interval_division(when, elapsed, &
        alarms(b) % reset_interval, reset_n, reset_rem)

      if (reset_rem == zero_intv) then
        ! reset right now
        reset_time = current_time + alarms(b) % reset_interval
        series % buffers(b) % reset_flag = 1
      else
        reset_rem = alarms(b) % reset_interval - reset_rem
        reset_time = current_time + reset_rem
      end if
    else
      reset_time = when
    end if

    ! set next duration time and flag
    when = alarms(b) % start_time + alarms(b) % duration_interval ! is offset
    if (current_time >= when) then
      elapsed = current_time - when
      call mpas_interval_division(when, elapsed, &
        alarms(b) % repeat_interval, & ! repeat is correct
        duration_n, duration_rem)

      if (duration_rem == zero_intv) then
        ! turn off accumulation
        duration_time = current_time + alarms(b) % repeat_interval ! repeat
      else
        duration_rem = alarms(b) % repeat_interval - duration_rem ! repeat
        duration_time = current_time + duration_rem ! remainder of repeat
      end if
    else
      duration_time = when
      duration_n = -1
    end if

    ! set next repeat time and flag
    when = alarms(b) % start_time + alarms(b) % repeat_interval
    if (current_time >= when) then
      elapsed = current_time - when
      call mpas_interval_division(when, elapsed, &
        alarms(b) % repeat_interval, repeat_n, repeat_rem)

      if (repeat_rem == zero_intv) then
        repeat_time = current_time + alarms(b) % repeat_interval
      else
        repeat_rem = alarms(b) % repeat_interval - repeat_rem
        repeat_time = current_time + repeat_rem
      end if
    else
      repeat_time = when
      repeat_n = -1
    end if

    ! accumulate now if in a window (both duration & repeat are untriggered)
    if ((duration_n == repeat_n) .and. &
        (series % buffers(b) % started_flag == 1)) then
      series % buffers(b) % accumulate_flag = 1
    end if

    !
    ! set the reoccurring timers
    !
    series % buffers(b) % duration_alarm_ID = trim(alarm_prefix) // &
      trim(DURATION_ALARM_PREFIX) // trim(buf_identifier)
    call mpas_add_clock_alarm(domain % clock, &
      series % buffers(b) % duration_alarm_ID, &
      duration_time, & ! duration sets the offset
      alarms(b) % repeat_interval, ierr=err) ! but repeat is interval

    series % buffers(b) % repeat_alarm_ID = trim(alarm_prefix) // &
      trim(REPEAT_ALARM_PREFIX) // trim(buf_identifier)
    call mpas_add_clock_alarm(domain % clock, &
      series % buffers(b) % repeat_alarm_ID, &
      repeat_time, &
      alarms(b) % repeat_interval, ierr=err)

    series % buffers(b) % reset_alarm_ID = trim(alarm_prefix) // &
      trim(RESET_ALARM_PREFIX) // trim(buf_identifier)
    call mpas_add_clock_alarm(domain % clock, &
      series % buffers(b) % reset_alarm_ID, &
      reset_time, &
      alarms(b) % reset_interval, ierr=err)
  end do
end subroutine set_alarms



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
subroutine set_times(series, alarms, clock, which, config, ok, err)
  ! input variables
  integer, intent(in) :: which
  character (len=StrKIND), pointer, intent(in) :: config

  ! input/output variables
  type (time_series_type), intent(inout) :: series
  type (MPAS_Clock_type), intent(inout) :: clock
  type (time_series_alarms_type), dimension(:), intent(inout) :: alarms

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
    if (b > series % number_of_buffers) then
      exit
    end if

    ! set the time
    if (which == START_TIMES) then
      if (time == INITIAL_TIME_TOKEN) then
        alarms(b) % start_time = &
          mpas_get_clock_time(clock, MPAS_START_TIME, err)
      else
        call mpas_set_time(alarms(b) % start_time, &
          dateTimeString=time, ierr=err)
      end if
    else if (which == DURATION_INTERVALS) then
      if (time == REPEAT_INTERVAL_TOKEN) then
        alarms(b) % duration_interval = alarms(b) % repeat_interval
      else
        call mpas_set_timeInterval(alarms(b) % duration_interval, &
            timeString=time, ierr=err)
      end if
    else if (which == REPEAT_INTERVALS) then
      if (time == RESET_INTERVAL_TOKEN) then
        alarms(b) % repeat_interval = alarms(b) % reset_interval
      else
        call mpas_set_timeInterval(alarms(b) % repeat_interval, &
            timeString=time, ierr=err)
      end if
    else
      call mpas_set_timeInterval(alarms(b) % reset_interval, &
          timeString=time, ierr=err)
    end if

    ! get the next time string
    call walk_string(next, time, ok)
  end do

  ! only ok if we parsed out as many as there are number of buffers
  ok = series % number_of_buffers == b
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
subroutine add_new_field(info, inname, outname, all_fields, amPool)!{{{
  ! input variables
  type (mpas_pool_field_info_type), intent(in) :: info
  character (len=StrKIND), intent(in) :: inname, outname

  ! input/output variables
  type (mpas_pool_type), intent(inout) :: all_fields, amPool

  ! output variables

  ! local variables

  ! duplicate field and add new field to pool
  if (info % nDims == 0) then
    call copy_field_0r(inname, all_fields, amPool, outname)
  else if (info % nDims == 1) then
    call copy_field_1r(inname, all_fields, amPool, outname)
  else if (info % nDims == 2) then
    call copy_field_2r(inname, all_fields, amPool, outname)
  else if (info % nDims == 3) then
    call copy_field_3r(inname, all_fields, amPool, outname)
  else if (info % nDims == 4) then
    call copy_field_4r(inname, all_fields, amPool, outname)
  else if (info % nDims == 5) then
    call copy_field_5r(inname, all_fields, amPool, outname)
  else
    call mpas_dmpar_global_abort(trim(CURRENT_CORE_NAME) // ' ERROR: ' // &
      'the impossible happened - tried to copy a real field "' // &
      trim(inname) // '" that does not have 0-5 dimensionality')
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
subroutine timer_checking(series, clock, err)!{{{
  ! input variables

  ! input/output variables
  type (time_series_type), intent(inout) :: series
  type (mpas_clock_type), intent(inout) :: clock

  ! output variables
  integer, intent(out) :: err

  ! local variables
  integer :: b

  ! start procedure
  err = 0

  do b = 1, series % number_of_buffers
    ! clear any resets
    if (series % buffers(b) % reset_flag == 1) then
      if (series % buffers(b) % accumulate_flag == 1) then
        series % buffers(b) % reset_flag = 0
      end if
    end if

    ! see if the started alarm is ringing
    if (trim(series % buffers(b) % start_alarm_ID) /= '') then
      if (mpas_is_alarm_ringing(clock, &
          series % buffers(b) % start_alarm_ID, ierr=err)) then
        call mpas_reset_clock_alarm(clock, &
          series % buffers(b) % start_alarm_ID, ierr=err)
        series % buffers(b) % started_flag = 1
        series % buffers(b) % reset_flag = 1
        series % buffers(b) % accumulate_flag = 1

        series % buffers(b) % start_alarm_ID = ''
      end if
    end if

    ! if we aren't started, cycle to next buffer
    if (series % buffers(b) % started_flag == 0) then
      cycle
    end if

    ! check various other alarms
    ! see if we need to reset
    if(mpas_is_alarm_ringing(clock, &
      series % buffers(b) % reset_alarm_ID, ierr=err)) then
      call mpas_reset_clock_alarm(clock, &
        series % buffers(b) % reset_alarm_ID, ierr=err)
      series % buffers(b) % reset_flag = 1
    end if

    ! turn off accumulation
    !
    ! duration needs to be >= 2 * compute_interval
    ! (a series can only be 2 or more)
    if (mpas_is_alarm_ringing(clock, &
        series % buffers(b) % duration_alarm_ID, ierr=err)) then
      call mpas_reset_clock_alarm(clock, &
        series % buffers(b) % duration_alarm_ID, ierr=err)
      series % buffers(b) % accumulate_flag = 0
    end if

    ! turn on accumulation
    ! (this is second, in case the duration and repeat
    ! overlaps on the same timer)
    if (mpas_is_alarm_ringing(clock, &
        series % buffers(b) % repeat_alarm_ID, ierr=err)) then
      call mpas_reset_clock_alarm(clock, &
        series % buffers(b) % repeat_alarm_ID, ierr=err)
      series % buffers(b) % accumulate_flag = 1
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
subroutine typed_operate(block, variable, buffers, operation)!{{{
  ! input variables
  type (block_type), pointer, intent(in) :: block
  integer, intent(in) :: operation
  type (time_series_variable_type), intent(in) :: variable
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers

  ! input/output variables

  ! output variables

  ! local variables
  type (mpas_pool_field_info_type) :: info

  ! get the info
  call mpas_pool_get_field_info(block % allFields, variable % input_name, info)

  ! switch based on the type, dimensionality, and operation
  if (info % nDims  == 0) then
    if (operation == AVG_OP) then
      call operate0r_avg(block, variable, buffers)
    else if (operation == MIN_OP) then
      call operate0r_min(block, variable, buffers)
    else
      call operate0r_max(block, variable, buffers)
    end if
  else if (info % nDims == 1) then
    if (operation == AVG_OP) then
      call operate1r_avg(block, variable, buffers)
    else if (operation == MIN_OP) then
      call operate1r_min(block, variable, buffers)
    else
      call operate1r_max(block, variable, buffers)
    end if
  else if (info % nDims == 2) then
    if (operation == AVG_OP) then
      call operate2r_avg(block, variable, buffers)
    else if (operation == MIN_OP) then
      call operate2r_min(block, variable, buffers)
    else
      call operate2r_max(block, variable, buffers)
    end if
  else if (info % nDims == 3) then
    if (operation == AVG_OP) then
      call operate3r_avg(block, variable, buffers)
    else if (operation == MIN_OP) then
      call operate3r_min(block, variable, buffers)
    else
      call operate3r_max(block, variable, buffers)
    end if
  else if (info % nDims == 4) then
    if (operation == AVG_OP) then
      call operate4r_avg(block, variable, buffers)
    else if (operation == MIN_OP) then
      call operate4r_min(block, variable, buffers)
    else
      call operate4r_max(block, variable, buffers)
    end if
  else if (info % nDims == 5) then
    if (operation == AVG_OP) then
      call operate5r_avg(block, variable, buffers)
    else if (operation == MIN_OP) then
      call operate5r_min(block, variable, buffers)
    else
      call operate5r_max(block, variable, buffers)
    end if
  else
    call mpas_dmpar_global_abort(trim(CURRENT_CORE_NAME) // ' ERROR: ' // &
      'the impossible happened - tried to operate on a real field "' // &
      trim(variable % input_name) // '" that does not have 0-5 dimensionality')
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

subroutine copy_field_0r(inname, all_fields, amPool, outname)!{{{
#include "time_series_inc/copy_field_1.inc"
  type (field0DReal), pointer :: src, dst
#include "time_series_inc/copy_field_2.inc"
end subroutine copy_field_0r!}}}

subroutine copy_field_1r(inname, all_fields, amPool, outname)!{{{
#include "time_series_inc/copy_field_1.inc"
  type (field1DReal), pointer :: src, dst
#include "time_series_inc/copy_field_2.inc"
end subroutine copy_field_1r!}}}

subroutine copy_field_2r(inname, all_fields, amPool, outname)!{{{
#include "time_series_inc/copy_field_1.inc"
  type (field2DReal), pointer :: src, dst
#include "time_series_inc/copy_field_2.inc"
end subroutine copy_field_2r!}}}

subroutine copy_field_3r(inname, all_fields, amPool, outname)!{{{
#include "time_series_inc/copy_field_1.inc"
  type (field3DReal), pointer :: src, dst
#include "time_series_inc/copy_field_2.inc"
end subroutine copy_field_3r!}}}

subroutine copy_field_4r(inname, all_fields, amPool, outname)!{{{
#include "time_series_inc/copy_field_1.inc"
  type (field4DReal), pointer :: src, dst
#include "time_series_inc/copy_field_2.inc"
end subroutine copy_field_4r!}}}

subroutine copy_field_5r(inname, all_fields, amPool, outname)!{{{
#include "time_series_inc/copy_field_1.inc"
  type (field5DReal), pointer :: src, dst
#include "time_series_inc/copy_field_2.inc"
end subroutine copy_field_5r!}}}



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
subroutine operate0r_avg (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = (out_array * &
    (buffers(b) % counter - 1) + in_array) &
    / buffers(b) % counter ;
#include "time_series_inc/operate_3.inc"
end subroutine operate0r_avg

subroutine operate1r_avg (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = (out_array * &
    (buffers(b) % counter - 1) + in_array) &
    / buffers(b) % counter ;
#include "time_series_inc/operate_3.inc"
end subroutine operate1r_avg

subroutine operate2r_avg (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:, :), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = (out_array * &
    (buffers(b) % counter - 1) + in_array) &
    / buffers(b) % counter ;
#include "time_series_inc/operate_3.inc"
end subroutine operate2r_avg

subroutine operate3r_avg (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:, :, :), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = (out_array * &
    (buffers(b) % counter - 1) + in_array) &
    / buffers(b) % counter ;
#include "time_series_inc/operate_3.inc"
end subroutine operate3r_avg

subroutine operate4r_avg (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = (out_array * &
    (buffers(b) % counter - 1) + in_array) &
    / buffers(b) % counter ;
#include "time_series_inc/operate_3.inc"
end subroutine operate4r_avg

subroutine operate5r_avg (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = (out_array * &
    (buffers(b) % counter - 1) + in_array) &
    / buffers(b) % counter ;
#include "time_series_inc/operate_3.inc"
end subroutine operate5r_avg

subroutine operate0r_min (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = min(out_array, in_array) ;
#include "time_series_inc/operate_3.inc"
end subroutine operate0r_min

subroutine operate1r_min (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = min(out_array, in_array) ;
#include "time_series_inc/operate_3.inc"
end subroutine operate1r_min

subroutine operate2r_min (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:, :), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = min(out_array, in_array) ;
#include "time_series_inc/operate_3.inc"
end subroutine operate2r_min

subroutine operate3r_min (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:, :, :), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = min(out_array, in_array) ;
#include "time_series_inc/operate_3.inc"
end subroutine operate3r_min

subroutine operate4r_min (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = min(out_array, in_array) ;
#include "time_series_inc/operate_3.inc"
end subroutine operate4r_min

subroutine operate5r_min (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = min(out_array, in_array) ;
#include "time_series_inc/operate_3.inc"
end subroutine operate5r_min

subroutine operate0r_max (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = max(out_array, in_array) ;
#include "time_series_inc/operate_3.inc"
end subroutine operate0r_max

subroutine operate1r_max (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = max(out_array, in_array) ;
#include "time_series_inc/operate_3.inc"
end subroutine operate1r_max

subroutine operate2r_max (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:, :), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = max(out_array, in_array) ;
#include "time_series_inc/operate_3.inc"
end subroutine operate2r_max

subroutine operate3r_max (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:, :, :), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = max(out_array, in_array) ;
#include "time_series_inc/operate_3.inc"
end subroutine operate3r_max

subroutine operate4r_max (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:, :, :, :), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = max(out_array, in_array) ;
#include "time_series_inc/operate_3.inc"
end subroutine operate4r_max

subroutine operate5r_max (start_block, variable, buffers)
#include "time_series_inc/operate_1.inc"
  real (kind=RKIND), dimension(:, :, :, :, :), pointer :: in_array, out_array
#include "time_series_inc/operate_2.inc"
  out_array = max(out_array, in_array) ;
#include "time_series_inc/operate_3.inc"
end subroutine operate5r_max

end module ocn_time_series_stats
! vim: foldmethod=marker
