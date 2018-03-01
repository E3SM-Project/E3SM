  ! abstract interface for variable interface function
  abstract interface
     function variable_interval(currentTime) result(variableInterval)
       import MPAS_Time_type
       import MPAS_TimeInterval_type
       type(MPAS_Time_type), intent(in) :: currentTime
       type(MPAS_TimeInterval_type) :: variableInterval
     end function variable_interval
  end interface

  ! individual forcing field associated with a forcing stream
  type :: mpas_forcing_field_type

     ! variable identificaion
     character(len=strKIND) :: forcingName

     character(len=strKIND) :: poolname  ! pool name of the input field
     character(len=strKIND) :: fieldname ! field name of the output field

     ! linked list next pointer
     type (mpas_forcing_field_type), pointer :: next => null()

  end type mpas_forcing_field_type

  ! forcing stream associated with a forcing group
  type :: mpas_forcing_stream_type

     ! stream identification
     character(len=strKIND) :: forcingStreamID ! the stream ID for this forcing stream

     ! alarm for reading in the next forcing time
     character(len=ShortStrKind) :: forcingAlarmID ! the alarm ID for this forcing stream

     ! interpolation
     character(len=strKIND) :: interpolationType ! the interpolation type (e.g. 'linear', 'constant')
     integer :: nTimeStencil ! the number of forcing data times
     integer :: nTimeStencilLower ! the number of forcing data times less than the current forcing time
     integer :: nTimeStencilUpper ! the number of forcing data times more than the current forcing time
     type(MPAS_time_type), dimension(:), allocatable :: forcingTimes ! the forcing data times
     
     ! forcing times definition
     type(MPAS_TimeInterval_type) :: forcingIntervalConstant ! the forcing interval
     type(MPAS_Time_type) :: forcingReferenceTime ! the forcing reference time

     ! forcing initialization type
     character(len=strKIND) :: forcingInitializationType

     ! optional functions to calculate variable intervals
     procedure (variable_interval), pointer, nopass :: variable_interval_forward_ptr => null ()
     procedure (variable_interval), pointer, nopass :: variable_interval_backward_ptr => null ()

     ! function pointer to new functions to test same as stream ones
     procedure (variable_interval), pointer, nopass :: variable_interval_forward_test_ptr => null ()
     procedure (variable_interval), pointer, nopass :: variable_interval_backward_test_ptr => null ()

     ! linked list of individual forcing fields
     type (mpas_forcing_field_type), pointer :: field => null()

     ! linked list next pointer
     type (mpas_forcing_stream_type), pointer :: next => null()

  end type mpas_forcing_stream_type

  ! collection of forcing steams and fields using the same forcing clock
  type :: mpas_forcing_group_type

     ! the forcings name
     character(len=strKIND) :: forcingGroupName ! the forcing group identifying name

     ! pointer to the forcing group domain
     type(domain_type), pointer :: domain_ptr

     ! forcing clock and cycling times
     type (MPAS_Clock_type) :: forcingClock ! the forcing clock

     logical                       :: forcingCycleUse            ! whether we cycle the forcing clock
     type (MPAS_Time_Type)         :: forcingCycleStart          ! the start year of the forcing cycle
     type (MPAS_Time_Type)         :: forcingCycleEnd            ! the end year of the forcing cycle
     type (MPAS_TimeInterval_Type) :: forcingCycleDuration       ! the duration of the forcing cycle
     logical                       :: forcingCycleStartInclusive ! cycle start time is inclusive to the cycle

      ! the alarm ID for cycling the clock
     character(len=ShortStrKIND) :: forcingCycleAlarmID = "forcingCycleAlarmID"

     ! linked list of individual streams
     type (mpas_forcing_stream_type), pointer :: stream => null()

     ! linked list next pointer
     type (mpas_forcing_group_type), pointer :: next => null()

  end type mpas_forcing_group_type

