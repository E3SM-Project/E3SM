! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! cice_time_series_stats
!
!> \brief MPAS ocean analysis core member: time_series_stats
!> \author Jon Woodring, modified for sea ice by Adrian K. Turner
!> \date   September 1, 2015 (seaice: 10th November 2015)
!> \details
!>  Flexible time series averaging, mins, and maxes of fields.
!-----------------------------------------------------------------------
module cice_time_series_stats
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_dmpar
  use mpas_timekeeping
  use mpas_stream_manager

  !use ocn_constants
  !use ocn_diagnostics_routines

  implicit none
  private
  save

  ! Public parameters
  !--------------------------------------------------------------------

  ! Public member functions
  !--------------------------------------------------------------------
  public :: cice_bootstrap_time_series_stats, &
            cice_init_time_series_stats, &
            cice_precompute_time_series_stats, &
            cice_compute_time_series_stats, &
            cice_restart_time_series_stats, &
            cice_finalize_time_series_stats

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
    real (kind=RKIND), pointer :: counter
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

  integer, parameter :: START_TIMES = 5
  integer, parameter :: DURATION_INTERVALS = 6
  integer, parameter :: REPEAT_INTERVALS = 7
  integer, parameter :: RESET_INTERVALS = 8

  character (len=3), parameter :: AVG_TOKEN = 'avg'
  character (len=3), parameter :: MIN_TOKEN = 'min'
  character (len=3), parameter :: MAX_TOKEN = 'max'

  character (len=4), parameter :: MESH_STREAM = 'mesh'
  character (len=5), parameter :: TIME_STREAM = 'xtime'

  character (len=StrKIND), parameter :: ONE_STRING_MEMORY = &
    'timeSeriesStatsOneString'
  character (len=StrKIND), parameter :: ONE_INTEGER_MEMORY = &
    'timeSeriesStatsOneInteger'
  character (len=StrKIND), parameter :: ONE_REAL_MEMORY = &
    'timeSeriesStatsOneReal'

  character (len=StrKIND), parameter :: CONFIG_PREFIX = &
    'config_AM_timeSeriesStats'
  character (len=StrKIND), parameter :: FRAMEWORK_PREFIX = 'timeSeriesStats'

  character (len=StrKIND), parameter :: OUTPUT_STREAM_SUFFIX = '_output_stream'
  character (len=StrKIND), parameter :: RESTART_STREAM_SUFFIX = '_restart_stream'
  character (len=StrKIND), parameter :: OPERATION_SUFFIX = '_operation'
  character (len=StrKIND), parameter :: ADD_MESH_SUFFIX = '_add_mesh'

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
    '_counter_'

  character (len=StrKIND), parameter :: START_ALARM_PREFIX = '_startAlarm_'
  character (len=StrKIND), parameter :: REPEAT_ALARM_PREFIX = '_repeatAlarm_'
  character (len=StrKIND), parameter :: DURATION_ALARM_PREFIX = &
    '_durationAlarm_'
  character (len=StrKIND), parameter :: RESET_ALARM_PREFIX = '_resetAlarm_'

  character (len=StrKIND), parameter :: INITIAL_TIME_TOKEN = 'initial_time'
  character (len=StrKIND), parameter :: REPEAT_INTERVAL_TOKEN = &
    'repeat_interval'
  character (len=StrKIND), parameter :: RESET_INTERVAL_TOKEN = 'reset_interval'

!***********************************************************************
contains

!***********************************************************************
! routine cice_bootstrap_time_series_stats
!
!> \brief Bootstrap time_series_stats analysis member
!> \author  Doug Jacobsen
!> \date    10/08/2015
!> \details
!>  This routine performs pre-init configuration of the analysis member.
!>  Specifically, it ensures the streams used for this instance are correctly
!>  configured.
!-----------------------------------------------------------------------
subroutine cice_bootstrap_time_series_stats(domain, instance, err)!{{{
  ! input variables
  character(len=*), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: v
  !character (len=StrKIND) :: instance ! TODO intent(in)
  type (time_series_type) :: series
  type (time_series_alarms_type), allocatable, dimension(:) :: alarms

  ! start procedure
  err = 0

  ! TODO placeholder for some unique ID if this code is replicated
  !instance = '' ! TODO to be passed in

  ! initial allocation of instance state for this AM from the namelist
  call start_state(domain, instance, series, err)

  ! modify the output and restart streams for this AM instance
  ! driver will do a restart read, after this, if necessary to fill values
  call modify_stream(domain, instance, series, err)

  ! clean up the instance memory
  do v = 1, series % number_of_variables
    deallocate(series % variables(v) % output_names)
  end do
  deallocate(series % variables)
  deallocate(series % buffers)
end subroutine cice_bootstrap_time_series_stats!}}}

!***********************************************************************
! routine cice_init_time_series_stats
!
!> \brief Initialize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine cice_init_time_series_stats(domain, instance, err)!{{{
  ! input variables
  character(len=*), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: v
  !character (len=StrKIND) :: instance ! TODO intent(in)
  type (time_series_type) :: series
  type (time_series_alarms_type), allocatable, dimension(:) :: alarms

  ! start procedure
  err = 0

  ! TODO placeholder for some unique ID if this code is replicated
  !instance = '' ! TODO to be passed in

  ! coming back from a potential restart read
  ! get all of the state for this instance
  call get_state(domain, instance, series)

  ! get all of the timing configurations from namelist
  allocate(alarms(series % number_of_buffers))
  call get_alarms(domain, instance, series, alarms, err)

  ! set the values of the alarms and current flag states based on timers
  call set_alarms(domain, instance, series, alarms, err)
  deallocate(alarms)

  ! clean up the instance memory
  do v = 1, series % number_of_variables
    deallocate(series % variables(v) % output_names)
  end do
  deallocate(series % variables)
  deallocate(series % buffers)
end subroutine cice_init_time_series_stats!}}}

!***********************************************************************
!
!  routine cice_precompute_time_series_stats
!
!> \brief   Precompute MPAS-CICE analysis member
!> \author  Adrian K. Turner
!> \date    10th November 2015
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-CICE analysis member.
!
!-----------------------------------------------------------------------

   subroutine cice_precompute_time_series_stats(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine cice_precompute_time_series_stats!}}}

!***********************************************************************
! routine cice_compute_time_series_stats
!
!> \brief Compute MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine cice_compute_time_series_stats(domain, instance, timeLevel, err)!{{{
  ! input variables
  character(len=*), intent(in) :: instance
  integer, intent(in) :: timeLevel

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  !character (len=StrKIND) :: instance ! TODO intent(in)
  integer :: v, b
  type (time_series_type) :: series

  ! start procedure
  err = 0

  ! TODO placeholder for some unique ID if this code is replicated
  !instance = '' ! TODO to be passed in

  ! get all of the state for this instance to be able to compute
  call get_state(domain, instance, series)

  ! update the counter
  do b = 1, series % number_of_buffers
    if (series % buffers(b) % accumulate_flag == 1) then
      if (series % buffers(b) % reset_flag == 1) then
        series % buffers(b) % counter = 1
      else
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
end subroutine cice_compute_time_series_stats!}}}



!***********************************************************************
! routine cice_restart_time_series_stats
!
!> \brief Save restart for MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine cice_restart_time_series_stats(domain, instance, err)!{{{
  ! input variables
  character(len=*), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables

  ! start procedure
  err = 0

end subroutine cice_restart_time_series_stats!}}}



!***********************************************************************
! routine cice_finalize_time_series_stats
!
!> \brief Finalize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    September 1, 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
subroutine cice_finalize_time_series_stats(domain, instance, err)!{{{
  ! input variables
  character(len=*), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables

  ! start procedure
  err = 0

end subroutine cice_finalize_time_series_stats!}}}

!
! local subroutines
!

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
  character (len=*), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  type (time_series_type), intent(out) :: series

  ! local variables
  integer :: v, b
  character (len=StrKIND) :: storage_prefix, var_identifier, &
    buf_identifier, var_prefix, buf_prefix, field_name, op_name

  ! start procedure
  storage_prefix = trim(FRAMEWORK_PREFIX) // trim(instance)

  !
  ! get the base
  !

  ! number_of_variables
  field_name = trim(storage_prefix) // trim(NUMBER_OF_VARIABLES_SUFFIX)
  call mpas_pool_get_array(domain % blocklist % allFields, &
    field_name, series % number_of_variables, 1)

  ! number_of_buffers
  field_name = trim(storage_prefix) // trim(NUMBER_OF_BUFFERS_SUFFIX)
  call mpas_pool_get_array(domain % blocklist % allFields, &
    field_name, series % number_of_buffers, 1)

  ! operation
  field_name = trim(storage_prefix) // trim(OPERATION_SUFFIX)
  call mpas_pool_get_array(domain % blocklist % allFields, &
    field_name, series % operation, 1)

  ! operator
  if (series % operation == AVG_OP) then
    op_name = AVG_TOKEN
  else if (series % operation == MIN_OP) then
    op_name = MIN_TOKEN
  else
    op_name = MAX_TOKEN
  end if

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
    call mpas_pool_get_array(domain % blocklist % allFields, &
      field_name, series % variables(v) % input_name, 1)

    do b = 1, series % number_of_buffers
      write(buf_identifier, '(I0)') b

      ! create output names
      series % variables(v) % output_names(b) = output_naming &
        (storage_prefix, op_name, series % variables(v) % input_name, &
         buf_identifier)
    end do
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
    call mpas_pool_get_array(domain % blocklist % allFields, &
      field_name, series % buffers(b) % started_flag, 1)

    ! accumulate_flag
    field_name = trim(buf_prefix) // trim(ACCUMULATE_FLAG_SUFFIX)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      field_name, series % buffers(b) % accumulate_flag, 1)

    ! reset_flag
    field_name = trim(buf_prefix) // trim(RESET_FLAG_SUFFIX)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      field_name, series % buffers(b) % reset_flag, 1)

    ! start_alarm_ID
    field_name = trim(buf_prefix) // trim(START_ALARM_ID_SUFFIX)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      field_name, series % buffers(b) % start_alarm_ID, 1)

    ! repeat_alarm_ID
    field_name = trim(buf_prefix) // trim(REPEAT_ALARM_ID_SUFFIX)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      field_name, series % buffers(b) % repeat_alarm_ID, 1)

    ! duration_alarm_ID
    field_name = trim(buf_prefix) // trim(DURATION_ALARM_ID_SUFFIX)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      field_name, series % buffers(b) % duration_alarm_ID, 1)

    ! reset_alarm_ID
    field_name = trim(buf_prefix) // trim(RESET_ALARM_ID_SUFFIX)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      field_name, series % buffers(b) % reset_alarm_ID, 1)

    ! counter
    field_name = counter_naming(storage_prefix, buf_identifier)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      field_name, series % buffers(b) % counter, 1)
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
subroutine start_state(domain, instance, series, err)
  ! input variables
  character (len=*), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain

  ! output variables
  type (time_series_type), intent(out) :: series
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  character (len=StrKIND), pointer :: config_results, output_stream_name
  character (len=StrKIND) :: config, namelist_prefix, storage_prefix, &
    var_identifier, buf_identifier, var_prefix, buf_prefix
  integer :: b, v
  type (field0DChar), pointer :: srcString, dstString
  type (field0DInteger), pointer :: srcInteger, dstInteger

  ! start procedure
  err = 0

  namelist_prefix = trim(CONFIG_PREFIX) // trim(instance)
  storage_prefix = trim(FRAMEWORK_PREFIX) // trim(instance)

  !
  ! allocate some framework memory for instance state
  !

  ! number_of_variables
  call mpas_pool_get_field(domain % blocklist % allFields, &
    ONE_INTEGER_MEMORY, srcInteger, 1)
  call mpas_duplicate_field(srcInteger, dstInteger)
  dstInteger % fieldName = &
    trim(storage_prefix) // trim(NUMBER_OF_VARIABLES_SUFFIX)
  call mpas_pool_add_field(domain % blocklist % allFields, &
    dstInteger % fieldName, dstInteger)
  call mpas_pool_get_array(domain % blocklist % allFields, &
    dstInteger % fieldName, series % number_of_variables, 1)

  ! number_of_buffers
  call mpas_pool_get_field(domain % blocklist % allFields, &
    ONE_INTEGER_MEMORY, srcInteger, 1)
  call mpas_duplicate_field(srcInteger, dstInteger)
  dstInteger % fieldName = &
    trim(storage_prefix) // trim(NUMBER_OF_BUFFERS_SUFFIX)
  call mpas_pool_add_field(domain % blocklist % allFields, &
    dstInteger % fieldName, dstInteger)
  call mpas_pool_get_array(domain % blocklist % allFields, &
    dstInteger % fieldName, series % number_of_buffers, 1)

  ! operation
  call mpas_pool_get_field(domain % blocklist % allFields, &
    ONE_INTEGER_MEMORY, srcInteger, 1)
  call mpas_duplicate_field(srcInteger, dstInteger)
  dstInteger % fieldName = &
    trim(storage_prefix) // trim(OPERATION_SUFFIX)
  call mpas_pool_add_field(domain % blocklist % allFields, &
    dstInteger % fieldName, dstInteger)
  call mpas_pool_get_array(domain % blocklist % allFields, &
    dstInteger % fieldName, series % operation, 1)

  !
  ! assign some instance values
  !

  ! get the stream name
  config = trim(namelist_prefix) // trim(OUTPUT_STREAM_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, output_stream_name)

  if (output_stream_name == 'none') then
    call mpas_dmpar_global_abort('MPAS-seaice: Error: stream cannot be "none" ' // &
      'for time series stats.')
  end if

  ! count the number of variables
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)
  series % number_of_variables = 0
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    output_stream_name, config))
    series % number_of_variables = series % number_of_variables + 1
  end do

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
  if (config_results == AVG_TOKEN) then
    series % operation = AVG_OP
  else if (config_results == MIN_TOKEN) then
    series % operation = MIN_OP
  else if (config_results == MAX_TOKEN) then
    series % operation = MAX_OP
  else
    ! error if unknown operation
    call mpas_dmpar_global_abort('MPAS-seaice: Error: unknown operation in time ' // &
      'averaging analysis member configuration.')
  end if

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
  do v = 1, series % number_of_variables
    ! identifier
    write(var_identifier, '(I0)') v
    var_prefix = trim(storage_prefix) // '_' // trim(var_identifier)

    ! input_name
    call mpas_pool_get_field(domain % blocklist % allFields, &
      ONE_STRING_MEMORY, srcString, 1)
    call mpas_duplicate_field(srcString, dstString)
    dstString % fieldName = trim(var_prefix) // trim(INPUT_NAME_SUFFIX)
    call mpas_pool_add_field(domain % blocklist % allFields, &
      dstString % fieldName, dstString)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      dstString % fieldName, series % variables(v) % input_name, 1)
  end do

  ! create buffer space
  do b = 1, series % number_of_buffers
    ! identifier
    write(buf_identifier, '(I0)') b
    buf_prefix = trim(storage_prefix) // '_' // trim(buf_identifier)

    ! started_flag
    call mpas_pool_get_field(domain % blocklist % allFields, &
      ONE_INTEGER_MEMORY, srcInteger, 1)
    call mpas_duplicate_field(srcInteger, dstInteger)
    dstInteger % fieldName = trim(buf_prefix) // trim(STARTED_FLAG_SUFFIX)
    call mpas_pool_add_field(domain % blocklist % allFields, &
      dstInteger % fieldName, dstInteger)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      dstInteger % fieldName, series % buffers(b) % started_flag, 1)

    ! accumulate_flag
    call mpas_pool_get_field(domain % blocklist % allFields, &
      ONE_INTEGER_MEMORY, srcInteger, 1)
    call mpas_duplicate_field(srcInteger, dstInteger)
    dstInteger % fieldName = trim(buf_prefix) // trim(ACCUMULATE_FLAG_SUFFIX)
    call mpas_pool_add_field(domain % blocklist % allFields, &
      dstInteger % fieldName, dstInteger)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      dstInteger % fieldName, series % buffers(b) % accumulate_flag, 1)

    ! reset_flag
    call mpas_pool_get_field(domain % blocklist % allFields, &
      ONE_INTEGER_MEMORY, srcInteger, 1)
    call mpas_duplicate_field(srcInteger, dstInteger)
    dstInteger % fieldName = trim(buf_prefix) // trim(RESET_FLAG_SUFFIX)
    call mpas_pool_add_field(domain % blocklist % allFields, &
      dstInteger % fieldName, dstInteger)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      dstInteger % fieldName, series % buffers(b) % reset_flag, 1)

    ! start_alarm_ID
    call mpas_pool_get_field(domain % blocklist % allFields, &
      ONE_STRING_MEMORY, srcString, 1)
    call mpas_duplicate_field(srcString, dstString)
    dstString % fieldName = trim(buf_prefix) // trim(START_ALARM_ID_SUFFIX)
    call mpas_pool_add_field(domain % blocklist % allFields, &
      dstString % fieldName, dstString)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      dstString % fieldName, series % buffers(b) % start_alarm_ID, 1)

    ! repeat_alarm_ID
    call mpas_pool_get_field(domain % blocklist % allFields, &
      ONE_STRING_MEMORY, srcString, 1)
    call mpas_duplicate_field(srcString, dstString)
    dstString % fieldName = trim(buf_prefix) // trim(REPEAT_ALARM_ID_SUFFIX)
    call mpas_pool_add_field(domain % blocklist % allFields, &
      dstString % fieldName, dstString)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      dstString % fieldName, series % buffers(b) % repeat_alarm_ID, 1)

    ! duration_alarm_ID
    call mpas_pool_get_field(domain % blocklist % allFields, &
      ONE_STRING_MEMORY, srcString, 1)
    call mpas_duplicate_field(srcString, dstString)
    dstString % fieldName = trim(buf_prefix) // trim(DURATION_ALARM_ID_SUFFIX)
    call mpas_pool_add_field(domain % blocklist % allFields, &
      dstString % fieldName, dstString)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      dstString % fieldName, series % buffers(b) % duration_alarm_ID, 1)

    ! reset_alarm_ID
    call mpas_pool_get_field(domain % blocklist % allFields, &
      ONE_STRING_MEMORY, srcString, 1)
    call mpas_duplicate_field(srcString, dstString)
    dstString % fieldName = trim(buf_prefix) // trim(RESET_ALARM_ID_SUFFIX)
    call mpas_pool_add_field(domain % blocklist % allFields, &
      dstString % fieldName, dstString)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      dstString % fieldName, series % buffers(b) % reset_alarm_ID, 1)

    !
    ! counter is not allocated here, because it is part of the restart stream,
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
subroutine modify_stream(domain, instance, series, err)!{{{
  ! input variables
  character (len=*), intent(in) :: instance

  ! input/output variables
  type (domain_type), intent(inout) :: domain
  type (time_series_type), intent(inout) :: series

  ! output variables
  integer, intent(out) :: err !< Output: error flag

  ! local variables
  integer :: v, b
  logical :: emptyRestartStream
  character (len=StrKIND), pointer :: output_stream_name, restart_stream_name
  character (len=StrKIND) :: fieldName
  type (field0DReal), pointer :: srcReal, dstReal
  logical, pointer :: copy_mesh
  character (len=StrKIND) :: field_name, config, op_name
  character (len=StrKIND) :: namelist_prefix, &
    storage_prefix, buf_identifier, buf_prefix
  type (mpas_pool_field_info_type) :: info

  ! start procedure
  err = 0

  namelist_prefix = trim(CONFIG_PREFIX) // trim(instance)
  storage_prefix = trim(FRAMEWORK_PREFIX) // trim(instance)

  ! get the output stream name
  config = trim(namelist_prefix) // trim(OUTPUT_STREAM_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, output_stream_name)

  ! get restart stream name
  config = trim(namelist_prefix) // trim(RESTART_STREAM_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, restart_stream_name)

  ! operator
  if (series % operation == AVG_OP) then
    op_name = AVG_TOKEN
  else if (series % operation == MIN_OP) then
    op_name = MIN_TOKEN
  else
    op_name = MAX_TOKEN
  end if

  ! get the old field names
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)
  v = 1
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
    output_stream_name, field_name))
    series % variables(v) % input_name = field_name
    v = v + 1
  end do

  ! remove the old ones from the stream
  do v = 1, series % number_of_variables
    call mpas_stream_mgr_remove_field(domain % streamManager, &
      output_stream_name, series % variables(v) % input_name)
  end do

  !
  ! create memory and modify the stream
  !

  ! add xtime to the output stream
  call mpas_stream_mgr_add_field(domain % streamManager, &
    output_stream_name, TIME_STREAM, ierr=err)

  ! optionally add mesh to output stream
  config = trim(namelist_prefix) // trim(ADD_MESH_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, copy_mesh)
  if (copy_mesh) then
    call mpas_stream_mgr_begin_iteration(domain % streamManager, &
      MESH_STREAM, err)
    do while (mpas_stream_mgr_get_next_field(domain % streamManager, &
      MESH_STREAM, field_name))
      call mpas_stream_mgr_add_field(domain % streamManager, &
        output_stream_name, field_name, ierr=err)
    end do
  end if

  ! ensure restart stream is empty
  emptyRestartStream = .true.
  call mpas_stream_mgr_begin_iteration(domain % streamManager, streamID=restart_stream_name, ierr=err)
  do while (mpas_stream_mgr_get_next_field(domain % streamManager, streamID=restart_stream_name, fieldName=fieldName) .and. emptyRestartStream)
    emptyRestartStream = .false.
  end do

  if (.not. emptyRestartStream) then
    write(stderrUnit, *) 'ERROR: Stream named ''' // trim(restart_stream_name) // ''' is not empty, but is used in '
    write(stderrUnit, *) '       an instance of the time series stats analysis member. This stream will be built'
    write(stderrUnit, *) '       based on the contents of the ''' // trim(output_stream_name) // ''' stream.'
    write(stderrUnit, *) '       Please ensure it is empty in the streams file.'
    call mpas_dmpar_global_abort('MPAS-seaice: ERROR: Misconfigured streams for time series stats analysis member.')
  end if

  ! create and put the counters in the streams
  do b = 1, series % number_of_buffers
    write(buf_identifier, '(I0)') b

    field_name = counter_naming(storage_prefix, buf_identifier)

    ! allocate counter memory
    call mpas_pool_get_field(domain % blocklist % allFields, &
      ONE_REAL_MEMORY, srcReal, 1)
    call mpas_duplicate_field(srcReal, dstReal)
    dstReal % fieldName = field_name
    call mpas_pool_add_field(domain % blocklist % allFields, &
      dstReal % fieldName, dstReal)
    call mpas_pool_get_array(domain % blocklist % allFields, &
      dstReal % fieldName, series % buffers(b) % counter, 1)

    ! put it in the output stream
    call mpas_stream_mgr_add_field(domain % streamManager, &
      output_stream_name, dstReal % fieldName, ierr=err)

    ! put it in the restart stream
    call mpas_stream_mgr_add_field(domain % streamManager, &
      restart_stream_name, dstReal % fieldName, ierr=err)
  end do

  ! set up the variables
  call mpas_stream_mgr_begin_iteration(domain % streamManager, &
    output_stream_name, err)
  do v = 1, series % number_of_variables
    ! get the info of the field
    call mpas_pool_get_field_info(domain % blocklist % allFields, &
      series % variables(v) % input_name, info)

    ! check if we can handle it
    if(.not. ((info % fieldType == MPAS_POOL_REAL) &
       .or. (info % fieldType == MPAS_POOL_INTEGER))) then
      call mpas_dmpar_global_abort('MPAS-seaice: Error: field "' // &
        trim(series % variables(v) % input_name) // '" listed in the ' // &
        'output stream, for time series stats analysis member ' // &
        'stream, is not real or integer.')
    end if

    ! allocate a number of fields and add field
    do b = 1, series % number_of_buffers
      write(buf_identifier, '(I0)') b

      field_name = output_naming(storage_prefix, op_name, &
        series % variables(v) % input_name, buf_identifier)

      ! create the name of the output var
      series % variables(v) % output_names(b) = field_name

      ! create the field and add to pool
      call add_new_field(info, &
        series % variables(v) % input_name, &
        series % variables(v) % output_names(b), &
        domain % blocklist % allFields)

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
  character (len=StrKIND), intent(in) :: storage_prefix, op_name, &
    input_name, buf_identifier

  output_naming = trim(storage_prefix) // '_' // trim(op_name) // '_' // &
    trim(input_name) // '_' // trim(buf_identifier)
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
  character (len=StrKIND), intent(in) :: storage_prefix, buf_identifier

  counter_naming = trim(storage_prefix) // trim(COUNTER_SUFFIX) // &
    trim(buf_identifier)
end function counter_naming



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
  character (len=*), intent(in) :: instance

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
    call mpas_dmpar_global_abort('MPAS-seaice: Error: number of times in ' // &
      'reset_intervals is not consistent with number of times ' // &
      'in reference_times in time series stats analysis member ' // &
      'configuration.')
  end if

  ! configure repeat intervals
  config = trim(namelist_prefix) // trim(REPEAT_INTERVALS_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, config_results)
  call set_times(series, alarms, domain % clock, REPEAT_INTERVALS, &
    config_results, ok, err)
  if (.not. ok) then
    call mpas_dmpar_global_abort('MPAS-seaice: Error: number of times in ' // &
      'repeat_intervals is not consistent with number of times ' // &
      'in reference_times in time series stats analysis member ' // &
      'configuration.')
  end if

  ! configure duration intervals
  config = trim(namelist_prefix) // trim(DURATION_INTERVALS_SUFFIX)
  call mpas_pool_get_config(domain % configs, config, config_results)
  call set_times(series, alarms, domain % clock, DURATION_INTERVALS, &
    config_results, ok, err)
  if (.not. ok) then
    call mpas_dmpar_global_abort('MPAS-seaice: Error: number of times in ' // &
      'duration_intervals is not consistent with number of times ' // &
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
      write(stderrUnit,*) 'Warning: repeat_interval > ' // &
        'reset_interval in time averaging analysis member ' // &
        'configuration. Truncating repeat_interval.'
      alarms(b) % repeat_interval = alarms(b) % reset_interval
    end if

    call mpas_interval_division(alarms(b) % start_time, &
       alarms(b) % duration_interval, &
       alarms(b) % repeat_interval, n, rem)

    if (n > 1 .or. (n == 1 .and. rem /= zero)) then
      write(stderrUnit,*) 'Warning: duration_interval > ' // &
        'repeat_interval in time averaging analysis member ' // &
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
  character (len=*), intent(in) :: instance

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
subroutine add_new_field(info, inname, outname, pool)!{{{
  ! input variables
  type (mpas_pool_field_info_type), intent(in) :: info
  character (len=StrKIND), intent(in) :: inname, outname

  ! input/output variables
  type (mpas_pool_type), intent(inout) :: pool

  ! output variables

  ! local variables

  ! duplicate field and add new field to pool
  if (info % fieldType == MPAS_POOL_REAL) then
    if (info % nDims == 0) then
      call copy_field_0r(inname, pool, outname)
    else if (info % nDims == 1) then
      call copy_field_1r(inname, pool, outname)
    else if (info % nDims == 2) then
      call copy_field_2r(inname, pool, outname)
    else if (info % nDims == 3) then
      call copy_field_3r(inname, pool, outname)
    else if (info % nDims == 4) then
      call copy_field_4r(inname, pool, outname)
    else
      call copy_field_5r(inname, pool, outname)
    end if
  else
    if (info % nDims == 0) then
      call copy_field_0i(inname, pool, outname)
    else if (info % nDims == 1) then
      call copy_field_1i(inname, pool, outname)
    else if (info % nDims == 2) then
      call copy_field_2i(inname, pool, outname)
    else
      call copy_field_3i(inname, pool, outname)
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
  if (info % fieldType == MPAS_POOL_REAL) then
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
    else
      if (operation == AVG_OP) then
        call operate5r_avg(block, variable, buffers)
      else if (operation == MIN_OP) then
        call operate5r_min(block, variable, buffers)
      else
        call operate5r_max(block, variable, buffers)
      end if
    end if
  else
    if (info % nDims == 0) then
      if (operation == AVG_OP) then
        call operate0i_avg(block, variable, buffers)
      else if (operation == MIN_OP) then
        call operate0i_min(block, variable, buffers)
      else
        call operate0i_max(block, variable, buffers)
      end if
    else if (info % nDims == 1) then
      if (operation == AVG_OP) then
        call operate1i_avg(block, variable, buffers)
      else if (operation == MIN_OP) then
        call operate1i_min(block, variable, buffers)
      else
        call operate1i_max(block, variable, buffers)
      end if
    else if (info % nDims == 2) then
      if (operation == AVG_OP) then
        call operate2i_avg(block, variable, buffers)
      else if (operation == MIN_OP) then
        call operate2i_min(block, variable, buffers)
      else
        call operate2i_max(block, variable, buffers)
      end if
    else
      if (operation == AVG_OP) then
        call operate3i_avg(block, variable, buffers)
      else if (operation == MIN_OP) then
        call operate3i_min(block, variable, buffers)
      else
        call operate3i_max(block, variable, buffers)
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

subroutine copy_field_0r(inname, pool, outname)!{{{
  character (len=StrKIND), intent(in) :: inname, outname
  type (mpas_pool_type), intent(inout) :: pool

  type (field0DReal), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = outname

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(outname) // '_' // &
        trim(dst % constituentNames(i))
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_0r!}}}

subroutine copy_field_1r(inname, pool, outname)!{{{
  character (len=StrKIND), intent(in) :: inname, outname
  type (mpas_pool_type), intent(inout) :: pool

  type (field1DReal), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = outname

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(outname) // '_' // &
        trim(dst % constituentNames(i))
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_1r!}}}

subroutine copy_field_2r(inname, pool, outname)!{{{
  character (len=StrKIND), intent(in) :: inname, outname
  type (mpas_pool_type), intent(inout) :: pool

  type (field2DReal), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = outname

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(outname) // '_' // &
        trim(dst % constituentNames(i))
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_2r!}}}

subroutine copy_field_3r(inname, pool, outname)!{{{
  character (len=StrKIND), intent(in) :: inname, outname
  type (mpas_pool_type), intent(inout) :: pool

  type (field3DReal), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = outname

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(outname) // '_' // &
        trim(dst % constituentNames(i))
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_3r!}}}

subroutine copy_field_4r(inname, pool, outname)!{{{
  character (len=StrKIND), intent(in) :: inname, outname
  type (mpas_pool_type), intent(inout) :: pool

  type (field4DReal), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = outname

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(outname) // '_' // &
        trim(dst % constituentNames(i))
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_4r!}}}

subroutine copy_field_5r(inname, pool, outname)!{{{
  character (len=StrKIND), intent(in) :: inname, outname
  type (mpas_pool_type), intent(inout) :: pool

  type (field5DReal), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = outname

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(outname) // '_' // &
        trim(dst % constituentNames(i))
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_5r!}}}

subroutine copy_field_0i(inname, pool, outname)!{{{
  character (len=StrKIND), intent(in) :: inname, outname
  type (mpas_pool_type), intent(inout) :: pool

  type (field0DInteger), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = outname

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(outname) // '_' // &
        trim(dst % constituentNames(i))
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_0i!}}}

subroutine copy_field_1i(inname, pool, outname)!{{{
  character (len=StrKIND), intent(in) :: inname, outname
  type (mpas_pool_type), intent(inout) :: pool

  type (field1DInteger), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = outname

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(outname) // '_' // &
        trim(dst % constituentNames(i))
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_1i!}}}

subroutine copy_field_2i(inname, pool, outname)!{{{
  character (len=StrKIND), intent(in) :: inname, outname
  type (mpas_pool_type), intent(inout) :: pool

  type (field2DInteger), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = outname

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(outname) // '_' // &
        trim(dst % constituentNames(i))
    end do
  end if

  call mpas_pool_add_field(pool, dst % fieldName, dst)
end subroutine copy_field_2i!}}}

subroutine copy_field_3i(inname, pool, outname)!{{{
  character (len=StrKIND), intent(in) :: inname, outname
  type (mpas_pool_type), intent(inout) :: pool

  type (field3DInteger), pointer :: src, dst
  integer :: i

  call mpas_pool_get_field(pool, inname, src, 1)
  call mpas_duplicate_field(src, dst)

  dst % fieldName = outname

  if (dst % isVarArray) then
    do i = 1, size(dst % constituentNames)
      dst % constituentNames(i) = trim(outname) // '_' // &
        trim(dst % constituentNames(i))
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
subroutine operate0r_avg (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
        out_array = (out_array * &
                (buffers(b) % counter - 1) + in_array) &
                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate0r_avg

subroutine operate1r_avg (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
        out_array = (out_array * &
                (buffers(b) % counter - 1) + in_array) &
                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate1r_avg

subroutine operate2r_avg (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
        out_array = (out_array * &
                (buffers(b) % counter - 1) + in_array) &
                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate2r_avg

subroutine operate3r_avg (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
        out_array = (out_array * &
                (buffers(b) % counter - 1) + in_array) &
                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate3r_avg

subroutine operate4r_avg (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:,:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
        out_array = (out_array * &
                (buffers(b) % counter - 1) + in_array) &
                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate4r_avg

subroutine operate5r_avg (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:,:,:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
        out_array = (out_array * &
                (buffers(b) % counter - 1) + in_array) &
                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate5r_avg

subroutine operate0i_avg (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  integer, pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
        out_array = (out_array * &
                (buffers(b) % counter - 1) + in_array) &
                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate0i_avg

subroutine operate1i_avg (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  integer, dimension(:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
        out_array = (out_array * &
                (buffers(b) % counter - 1) + in_array) &
                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate1i_avg

subroutine operate2i_avg (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  integer, dimension(:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
        out_array = (out_array * &
                (buffers(b) % counter - 1) + in_array) &
                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate2i_avg

subroutine operate3i_avg (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  integer, dimension(:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
        out_array = (out_array * &
                (buffers(b) % counter - 1) + in_array) &
                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate3i_avg

subroutine operate0r_min (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate0r_min

subroutine operate1r_min (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate1r_min

subroutine operate2r_min (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate2r_min

subroutine operate3r_min (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate3r_min

subroutine operate4r_min (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:,:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate4r_min

subroutine operate5r_min (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:,:,:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate5r_min

subroutine operate0i_min (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  integer, pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate0i_min

subroutine operate1i_min (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  integer, dimension(:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate1i_min

subroutine operate2i_min (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  integer, dimension(:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate2i_min

subroutine operate3i_min (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  integer, dimension(:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
        out_array = min(out_array, in_array) ;
!        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate3i_min

subroutine operate0r_max (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate0r_max

subroutine operate1r_max (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate1r_max

subroutine operate2r_max (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate2r_max

subroutine operate3r_max (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate3r_max

subroutine operate4r_max (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:,:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate4r_max

subroutine operate5r_max (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  real (kind=RKIND), dimension(:,:,:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate5r_max

subroutine operate0i_max (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  integer, pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate0i_max

subroutine operate1i_max (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  integer, dimension(:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate1i_max

subroutine operate2i_max (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  integer, dimension(:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate2i_max

subroutine operate3i_max (start_block, variable, buffers)
  type (block_type), pointer, intent(in) :: start_block
  type (time_series_buffer_type), dimension(:), intent(in) :: buffers
  type (time_series_variable_type), intent(in) :: variable

  integer, dimension(:,:,:), pointer :: in_array, out_array
  integer :: b
  type (block_type), pointer :: block

  block => start_block
  do while (associated(block))
    call mpas_pool_get_array(block % allFields, &
      variable % input_name, in_array, 1)

    do b = 1, size(buffers)
      if (buffers(b) % accumulate_flag == 0) then
        cycle
      end if

      call mpas_pool_get_array(block % allFields, &
          variable % output_names(b), out_array, 1)

      if (buffers(b) % reset_flag == 1) then
        out_array = in_array
      else
!        out_array = (out_array * &
!                (buffers(b) % counter - 1) + in_array) &
!                / buffers(b) % counter ;
!        out_array = min(out_array, in_array) ;
        out_array = max(out_array, in_array) ;

      end if
    end do

    block => block % next
  end do
end subroutine operate3i_max

end module cice_time_series_stats
! vim: foldmethod=marker
