   integer, parameter :: MPAS_MAX_ALARMS = 40

   integer, parameter :: MPAS_NOW = 0, &
                         MPAS_START_TIME = 1, &
                         MPAS_STOP_TIME = 2

   integer, parameter :: MPAS_FORWARD = 1, &
                         MPAS_BACKWARD = -1

   integer, parameter :: MPAS_GREGORIAN = 0, &
                         MPAS_GREGORIAN_NOLEAP = 1, &
                         MPAS_360DAY = 2

   type MPAS_Time_type
      type (ESMF_Time) :: t
   end type

   type MPAS_TimeInterval_type
      type (ESMF_TimeInterval) :: ti
   end type

   type MPAS_Alarm_type
      character (len=ShortStrKIND) :: alarmID
      logical :: isRecurring
      logical :: isSet
      type (MPAS_Time_type) :: ringTime
      type (MPAS_Time_type) :: prevRingTime
      type (MPAS_TimeInterval_type) :: ringTimeInterval
      type (MPAS_Alarm_type), pointer :: next => null()
   end type
   
   type MPAS_Clock_type
      integer :: direction
      integer :: nAlarms
      type (ESMF_Clock) :: c
      type (MPAS_Alarm_type), pointer :: alarmListHead => null()
   end type

