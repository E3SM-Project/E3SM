! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_timekeeping

   use mpas_kind_types
   use mpas_io_units
   use ESMF
   use ESMF_BaseMod
   use ESMF_Stubs
   use ESMF_CalendarMod
   use ESMF_ClockMod
   use ESMF_TimeMod
   use ESMF_TimeIntervalMod

   private :: mpas_calibrate_alarms
   private :: mpas_in_ringing_envelope

   integer, parameter :: MPAS_MAX_ALARMS = 20
   integer, parameter :: MPAS_NOW = 0, &
                         MPAS_START_TIME = 1, &
                         MPAS_STOP_TIME = 2
   integer, parameter :: MPAS_FORWARD = 1, &
                         MPAS_BACKWARD = -1
   integer, parameter :: MPAS_GREGORIAN = 0, &
                         MPAS_GREGORIAN_NOLEAP = 1, &
                         MPAS_360DAY = 2

   integer :: TheCalendar 
   integer :: yearWidth

   integer, dimension(12), parameter :: daysInMonth     = (/31,28,31,30,31,30,31,31,30,31,30,31/)
   integer, dimension(12), parameter :: daysInMonthLeap = (/31,29,31,30,31,30,31,31,30,31,30,31/)

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
      type (MPAS_Alarm_type), pointer :: next
   end type
   
   type MPAS_Clock_type
      integer :: direction
      integer :: nAlarms
      type (ESMF_Clock) :: c
      type (MPAS_Alarm_type), pointer :: alarmListHead
   end type

   interface operator (+)
      module procedure add_t_ti
      module procedure add_ti_ti
   end interface

   interface operator (-)
      module procedure sub_t_t
      module procedure sub_t_ti
      module procedure sub_ti_ti
      module procedure neg_ti
   end interface

   interface operator (*)
      module procedure mul_ti_n
   end interface

   interface operator (/)
      module procedure div_ti_n
   end interface

   interface operator (.EQ.)
      module procedure eq_t_t
      module procedure eq_ti_ti
   end interface

   interface operator (.NE.)
      module procedure ne_t_t
      module procedure ne_ti_ti
   end interface

   interface operator (.LT.)
      module procedure lt_t_t
      module procedure lt_ti_ti
   end interface

   interface operator (.GT.)
      module procedure gt_t_t
      module procedure gt_ti_ti
   end interface

   interface operator (.LE.)
      module procedure le_t_t
      module procedure le_ti_ti
   end interface

   interface operator (.GE.)
      module procedure ge_t_t
      module procedure ge_ti_ti
   end interface

   interface abs
      module procedure abs_ti
   end interface


   contains


   subroutine mpas_timekeeping_init(calendar)

      implicit none

      character (len=*), intent(in) :: calendar 

#ifndef MPAS_CESM
      if (trim(calendar) == 'gregorian') then
         TheCalendar = MPAS_GREGORIAN
         call ESMF_Initialize(defaultCalendar=ESMF_CALKIND_GREGORIAN)
      else if (trim(calendar) == 'gregorian_noleap') then
         TheCalendar = MPAS_GREGORIAN_NOLEAP
         call ESMF_Initialize(defaultCalendar=ESMF_CALKIND_NOLEAP)
!     else if (trim(calendar) == '360day') then
!        TheCalendar = MPAS_360DAY
!        call ESMF_Initialize(defaultCalendar=ESMF_CALKIND_360DAY)
      else
         write(stderrUnit,*) 'ERROR: mpas_timekeeping_init: Invalid calendar type'
      end if
#else
      if (trim(calendar) == 'gregorian') then
         TheCalendar = MPAS_GREGORIAN
      else if (trim(calendar) == 'gregorian_noleap') then
         TheCalendar = MPAS_GREGORIAN_NOLEAP
      else
         write(stderrUnit,*) 'ERROR: mpas_timekeeping_init: Invalid calendar type'
      end if
#endif

      yearWidth = 4

   end subroutine mpas_timekeeping_init


   subroutine mpas_timekeeping_finalize()

      implicit none

#ifndef MPAS_CESM      
      call ESMF_Finalize()
#endif

   end subroutine mpas_timekeeping_finalize

   !-----------------------------------------------------------------------
   !  routine mpas_timekeeping_set_year_width
   !
   !> \brief This routine sets the width of the year portion of timestamps.
   !> \author Michael Duda, Doug Jacobsen
   !> \date   07/23/2014
   !> \details This routine sets the width of the year portion of timestamps.
   !>   It can be used to make the year portion of a time stamp or an expanded
   !>   string more than 4 digits, to support years larger than 9999.
   !>
   !-----------------------------------------------------------------------
   subroutine mpas_timekeeping_set_year_width(yearWidthIn)!{{{
      integer, intent(in) :: yearWidthIn

      yearWidth = yearWidthIn

      if (yearWidthIn <= 0) then
          write(stderrUnit,*) 'ERROR: mpas_set_year_width: yearWidth cannot be less than or equal to zero.' 
          ierr = 1
          return
      end if

      yearWidth = yearWidthIn

      call ESMF_setYearWidth(yearWidthIn)

   end subroutine mpas_timekeeping_set_year_width!}}}

   subroutine mpas_create_clock(clock, startTime, timeStep, stopTime, runDuration, ierr)

      implicit none

      type (MPAS_Clock_type), intent(out) :: clock
      type (MPAS_Time_type), intent(in) :: startTime
      type (MPAS_TimeInterval_type), intent(in) :: timeStep
      type (MPAS_Time_type), intent(in), optional :: stopTime
      type (MPAS_TimeInterval_type), intent(in), optional :: runDuration
      integer, intent(out), optional :: ierr

      type (MPAS_Time_type) :: stop_time

      if (present(runDuration)) then
         stop_time = startTime + runDuration
         if (present(stopTime)) then
            if (stopTime /= stop_time) then
               if (present(ierr)) ierr = 1   ! stopTime and runDuration are inconsistent
               write(stderrUnit,*) 'ERROR: MPAS_createClock: stopTime and runDuration are inconsistent'
               return
            end if
         end if
      else if (present(stopTime)) then 
         stop_time = stopTime
      else
         if (present(ierr)) ierr = 1   ! neither stopTime nor runDuration are specified
         write(stderrUnit,*) 'ERROR: MPAS_createClock: neither stopTime nor runDuration are specified'
         return
      end if

      clock % c = ESMF_ClockCreate(TimeStep=timeStep%ti, StartTime=startTime%t, StopTime=stop_time%t, rc=ierr)
      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if
      clock % direction = MPAS_FORWARD
      clock % nAlarms = 0
      nullify(clock % alarmListHead)

   end subroutine mpas_create_clock


   subroutine mpas_destroy_clock(clock, ierr)

      implicit none

      type (MPAS_Clock_type), intent(inout) :: clock
      integer, intent(out), optional :: ierr

      type (MPAS_Alarm_type), pointer :: alarmPtr

      alarmPtr => clock % alarmListHead
      do while (associated(alarmPtr))
         clock % alarmListHead => alarmPtr % next
         deallocate(alarmPtr)
         alarmPtr => clock % alarmListHead
      end do

      call ESMF_ClockDestroy(clock % c, rc=ierr)
      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if

   end subroutine mpas_destroy_clock


   logical function mpas_is_clock_start_time(clock, ierr)

      implicit none

      type (MPAS_Clock_type), intent(in) :: clock
      integer, intent(out), optional :: ierr

      type (ESMF_Time) :: currTime, startTime, stopTime

      call ESMF_ClockGet(clock % c, CurrTime=currTime, rc=ierr)
      call ESMF_ClockGet(clock % c, StartTime=startTime, rc=ierr)
      call ESMF_ClockGet(clock % c, StopTime=stopTime, rc=ierr)
      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if

      if (startTime <= stopTime) then
         mpas_is_clock_start_time = (currTime <= startTime)
      else
         mpas_is_clock_start_time = (currTime >= startTime)
      end if

   end function mpas_is_clock_start_time


   logical function mpas_is_clock_stop_time(clock, ierr)

      implicit none

      type (MPAS_Clock_type), intent(in) :: clock
      integer, intent(out), optional :: ierr

      type (ESMF_Time) :: currTime, startTime, stopTime

      call ESMF_ClockGet(clock % c, CurrTime=currTime, rc=ierr)
      call ESMF_ClockGet(clock % c, StartTime=startTime, rc=ierr)
      call ESMF_ClockGet(clock % c, StopTime=stopTime, rc=ierr)
      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if

      if (startTime <= stopTime) then
         mpas_is_clock_stop_time = (currTime >= stopTime)
      else
         mpas_is_clock_stop_time = (currTime <= stopTime)
      end if

   end function mpas_is_clock_stop_time


   subroutine mpas_set_clock_direction(clock, direction, ierr)

      implicit none

      type (MPAS_Clock_type), intent(inout) :: clock
      integer, intent(in) :: direction
      integer, intent(out), optional :: ierr

      type (MPAS_TimeInterval_type) :: timeStep

      if (direction == MPAS_FORWARD .and. clock % direction == MPAS_FORWARD) return
      if (direction == MPAS_BACKWARD .and. clock % direction == MPAS_BACKWARD) return

      clock % direction = direction
      call ESMF_ClockGet(clock % c, TimeStep=timeStep%ti, rc=ierr)
      timeStep = neg_ti(timeStep)
      call ESMF_ClockSet(clock % c, TimeStep=timeStep%ti, rc=ierr)

      ! specify a valid previousRingTime for each alarm
      call mpas_calibrate_alarms(clock, ierr);

      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if

   end subroutine mpas_set_clock_direction



   integer function mpas_get_clock_direction(clock, ierr)

      implicit none

      type (MPAS_Clock_type), intent(in) :: clock
      integer, intent(out), optional :: ierr

      if (present(ierr)) ierr = 0

      mpas_get_clock_direction = clock % direction

   end function mpas_get_clock_direction


   subroutine mpas_set_clock_timestep(clock, timeStep, ierr)

      implicit none

      type (MPAS_Clock_type), intent(inout) :: clock
      type (MPAS_TimeInterval_type), intent(in) :: timeStep
      integer, intent(out), optional :: ierr

      call ESMF_ClockSet(clock % c, TimeStep=timeStep%ti, rc=ierr)
      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if

   end subroutine mpas_set_clock_timestep


   type (MPAS_TimeInterval_type) function mpas_get_clock_timestep(clock, ierr)

      implicit none

      type (MPAS_Clock_type), intent(in) :: clock
      integer, intent(out), optional :: ierr

      type (MPAS_TimeInterval_type) :: timeStep

      call ESMF_ClockGet(clock % c, TimeStep=timeStep%ti, rc=ierr)
      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if

      mpas_get_clock_timestep = timeStep

   end function mpas_get_clock_timestep


   subroutine mpas_advance_clock(clock, timeStep, ierr)

      implicit none

      type (MPAS_Clock_type), intent(inout) :: clock
      type (MPAS_TimeInterval_type), intent(in), optional :: timeStep
      integer, intent(out), optional :: ierr

      type (ESMF_TimeInterval) :: time_step

      if (present(timeStep)) then
         call ESMF_ClockGet(clock % c, TimeStep=time_step, rc=ierr)
         call ESMF_ClockSet(clock % c, TimeStep=timeStep % ti, rc=ierr)
         call ESMF_ClockAdvance(clock % c, rc=ierr)
         call ESMF_ClockSet(clock % c, TimeStep=time_step, rc=ierr)
      else
         call ESMF_ClockAdvance(clock % c, rc=ierr)
      end if
      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if

   end subroutine mpas_advance_clock


   subroutine mpas_set_clock_time(clock, clock_time, whichTime, ierr)

      implicit none

      type (MPAS_Clock_type), intent(inout) :: clock
      type (MPAS_Time_type), intent(in) :: clock_time
      integer, intent(in) :: whichTime
      integer, intent(out), optional :: ierr

      if (whichTime == MPAS_NOW) then
         call ESMF_ClockSet(clock % c, CurrTime=clock_time%t, rc=ierr)
         call mpas_calibrate_alarms(clock, ierr);
      else if (whichTime == MPAS_START_TIME) then
         call ESMF_ClockSet(clock % c, StartTime=clock_time%t, rc=ierr)
      else if (whichTime == MPAS_STOP_TIME) then
         call ESMF_ClockSet(clock % c, StopTime=clock_time%t, rc=ierr)
      else if (present(ierr)) then
         ierr = 1
      end if
      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if

   end subroutine mpas_set_clock_time


   type (MPAS_Time_type) function mpas_get_clock_time(clock, whichTime, ierr)

      implicit none

      type (MPAS_Clock_type), intent(in) :: clock
      integer, intent(in) :: whichTime
      integer, intent(out), optional :: ierr

      type (MPAS_Time_type) :: clock_time

      if (whichTime == MPAS_NOW) then
         call ESMF_ClockGet(clock % c, CurrTime=clock_time%t, rc=ierr)
      else if (whichTime == MPAS_START_TIME) then
         call ESMF_ClockGet(clock % c, StartTime=clock_time%t, rc=ierr)
      else if (whichTime == MPAS_STOP_TIME) then
         call ESMF_ClockGet(clock % c, StopTime=clock_time%t, rc=ierr)
      else if (present(ierr)) then
         ierr = 1
      end if
      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if

      mpas_get_clock_time = clock_time

   end function mpas_get_clock_time


   subroutine mpas_add_clock_alarm(clock, alarmID, alarmTime, alarmTimeInterval, ierr)
! TODO: possibly add a stop time for recurring alarms

      implicit none

      type (MPAS_Clock_type), intent(inout) :: clock
      character (len=*), intent(in) :: alarmID
      type (MPAS_Time_type), intent(in) :: alarmTime
      type (MPAS_TimeInterval_type), intent(in), optional :: alarmTimeInterval
      integer, intent(out), optional :: ierr

      type (MPAS_Alarm_type), pointer :: alarmPtr

      ! Add a new entry to the linked list of alarms for this clock
      if (.not. associated(clock % alarmListHead)) then
         allocate(clock % alarmListHead)
         nullify(clock % alarmListHead % next)
         alarmPtr => clock % alarmListHead
      else
         alarmPtr => clock % alarmListHead
         do while (associated(alarmPtr % next))
            if (trim(alarmPtr % alarmID) == trim(alarmID)) then
               write(stderrUnit,*) 'OOPS -- we have a duplicate alarmID', trim(alarmID)
               if (present(ierr)) ierr = 1
               return
            end if
            alarmPtr => alarmPtr % next
         end do
            if (trim(alarmPtr % alarmID) == trim(alarmID)) then
               write(stderrUnit,*) 'OOPS -- we have a duplicate alarmID', trim(alarmID)
               if (present(ierr)) ierr = 1
               return
            end if
         allocate(alarmPtr % next)
         alarmPtr => alarmPtr % next
         nullify(alarmPtr % next)
      end if

      alarmPtr % alarmID = trim(alarmID)

      clock % nAlarms = clock % nAlarms + 1

      alarmPtr % isSet = .true.
      alarmPtr % ringTime = alarmTime
      

      if (present(alarmTimeInterval)) then
         alarmPtr % isRecurring = .true.
         alarmPtr % ringTimeInterval = alarmTimeInterval
         if(clock % direction == MPAS_FORWARD) then
            alarmPtr % prevRingTime = alarmTime - alarmTimeInterval
         else
            alarmPtr % prevRingTime = alarmTime + alarmTimeInterval         
         end if
      else
         alarmPtr % isRecurring = .false.
         alarmPtr % prevRingTime = alarmTime
      end if
      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if

   end subroutine mpas_add_clock_alarm


   subroutine mpas_remove_clock_alarm(clock, alarmID, ierr)

      implicit none

      type (MPAS_Clock_type), intent(inout) :: clock
      character (len=*), intent(in) :: alarmID
      integer, intent(out), optional :: ierr

      type (MPAS_Alarm_type), pointer :: alarmPtr
      type (MPAS_Alarm_type), pointer :: alarmParentPtr

      if (present(ierr)) ierr = 0

      alarmPtr => clock % alarmListHead
      alarmParentPtr => alarmPtr
      do while (associated(alarmPtr))
         if (trim(alarmPtr % alarmID) == trim(alarmID)) then
            if (trim(alarmPtr % alarmID) == trim(clock % alarmListHead % alarmID)) then
               clock % alarmListHead => alarmPtr % next
            else
               alarmParentPtr % next => alarmPtr % next
            end if
            deallocate(alarmPtr)
            exit
         end if
         alarmParentPtr => alarmPtr
         alarmPtr => alarmPtr % next
      end do

   end subroutine mpas_remove_clock_alarm


   !-----------------------------------------------------------------------
   !  routine mpas_is_alarm_defined
   !
   !> \brief Check whether an alarm has been defined on a clock
   !> \author Michael Duda
   !> \date   26 August 2014
   !> \details
   !>  For a specified clock and alarm ID, checks whether that alarm ID has
   !>  been defined on the clock and returns the result.
   !
   !-----------------------------------------------------------------------
   logical function mpas_is_alarm_defined(clock, alarmID, ierr)

      implicit none

      type (MPAS_Clock_type), intent(in) :: clock
      character (len=*), intent(in) :: alarmID
      integer, intent(out) :: ierr

      type (MPAS_Alarm_type), pointer :: alarmPtr

      ierr = 0
      mpas_is_alarm_defined = .false.

      alarmPtr => clock % alarmListHead
      do while (associated(alarmPtr))
         if (trim(alarmPtr % alarmID) == trim(alarmID)) then
            mpas_is_alarm_defined = .true.
            return
         end if
         alarmPtr => alarmPtr % next 
      end do

   end function mpas_is_alarm_defined


   !-----------------------------------------------------------------------
   !  routine mpas_alarm_interval
   !
   !> \brief Retrieve the interval for an alarm
   !> \author Michael Duda
   !> \date   4 September 2014
   !> \details
   !>  For a specified clock and alarm ID, returns the time interval
   !>  associated with the alarm.
   !
   !-----------------------------------------------------------------------
   type (MPAS_TimeInterval_type) function mpas_alarm_interval(clock, alarmID, ierr)

      implicit none

      type (MPAS_Clock_type), intent(in) :: clock
      character (len=*), intent(in) :: alarmID
      integer, intent(out) :: ierr

      type (MPAS_Alarm_type), pointer :: alarmPtr

      ierr = 1
      call mpas_set_timeInterval(mpas_alarm_interval, S=0)

      alarmPtr => clock % alarmListHead
      do while (associated(alarmPtr))
         if (trim(alarmPtr % alarmID) == trim(alarmID)) then
            if (alarmPtr % isSet .and. alarmPtr % isRecurring) then
               ierr = 0
               mpas_alarm_interval = alarmPtr % ringTimeInterval
            end if
            return
         end if
         alarmPtr => alarmPtr % next 
      end do

   end function mpas_alarm_interval


   subroutine mpas_print_alarm(clock, alarmID, ierr)

      implicit none

      type (MPAS_Clock_type), intent(in) :: clock
      character (len=*), intent(in) :: alarmID
      integer, intent(out) :: ierr

      type (MPAS_Alarm_type), pointer :: alarmPtr

      character (len=StrKIND) :: printString

      ierr = 0

      alarmPtr => clock % alarmListHead
      do while (associated(alarmPtr))
         if (trim(alarmPtr % alarmID) == trim(alarmID)) then
            write(stderrUnit,*) 'ALARM ', trim(alarmID)

            write(stderrUnit,*) 'isRecurring', alarmPtr % isRecurring
            
            write(stderrUnit,*) 'isSet', alarmPtr % isSet

            call mpas_get_time(alarmPtr % ringTime, dateTimeString=printString, ierr=ierr)
            write(stderrUnit,*) 'ringTime', printString

            call mpas_get_time(alarmPtr % prevRingTime, dateTimeString=printString, ierr=ierr)
            write(stderrUnit,*) 'prevRingTime', printString

            call mpas_get_timeInterval(alarmPtr % ringTimeInterval, timeString=printString, ierr=ierr)
            write(stderrUnit,*) 'ringTimeInterval', printString
            
            exit
         end if
         alarmPtr => alarmPtr % next
      end do

   end subroutine mpas_print_alarm



   logical function mpas_is_alarm_ringing(clock, alarmID, interval, ierr)

      implicit none

      type (MPAS_Clock_type), intent(in) :: clock
      character (len=*), intent(in) :: alarmID
      type (MPAS_TimeInterval_type), intent(in), optional :: interval
      integer, intent(out), optional :: ierr

      type (MPAS_Alarm_type), pointer :: alarmPtr

      if (present(ierr)) ierr = 0

      mpas_is_alarm_ringing = .false.
      
      alarmPtr => clock % alarmListHead
      do while (associated(alarmPtr))
         if (trim(alarmPtr % alarmID) == trim(alarmID)) then
            if (alarmPtr % isSet) then
               if (mpas_in_ringing_envelope(clock, alarmPtr, interval, ierr)) then
                  mpas_is_alarm_ringing = .true.
               end if
            end if
            exit
         end if
         alarmPtr => alarmPtr % next
      end do

   end function mpas_is_alarm_ringing



   subroutine mpas_get_clock_ringing_alarms(clock, nAlarms, alarmList, interval, ierr)

      implicit none

      type (MPAS_Clock_type), intent(in) :: clock
      integer, intent(out) :: nAlarms
      character (len=ShortStrKIND), dimension(MPAS_MAX_ALARMS), intent(out) :: alarmList
      type (MPAS_TimeInterval_type), intent(in), optional :: interval
      integer, intent(out), optional :: ierr

      type (MPAS_Alarm_type), pointer :: alarmPtr

      if (present(ierr)) ierr = 0

      nAlarms = 0

      alarmPtr => clock % alarmListHead
      do while (associated(alarmPtr))
         if (alarmPtr % isSet) then
            if (mpas_in_ringing_envelope(clock, alarmPtr, interval, ierr)) then
               nAlarms = nAlarms + 1
               alarmList(nAlarms) = trim(alarmPtr % alarmID)
            end if
         end if
         alarmPtr => alarmPtr % next
      end do

   end subroutine mpas_get_clock_ringing_alarms


   logical function mpas_in_ringing_envelope(clock, alarmPtr, interval, ierr)

      implicit none
      
      type (MPAS_Clock_type), intent(in) :: clock
      type (MPAS_Alarm_type), pointer :: alarmPtr
      type (MPAS_TimeInterval_type), intent(in), optional :: interval
      integer, intent(out), optional :: ierr
      
      type (MPAS_Time_type) :: alarmNow
      type (MPAS_Time_type) :: alarmThreshold

      alarmNow = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      alarmThreshold = alarmPtr % ringTime 
      
      mpas_in_ringing_envelope = .false.      
               
      if(clock % direction == MPAS_FORWARD) then

         if (present(interval)) then
            alarmNow = alarmNow + interval; 
         end if

         if (alarmPtr % isRecurring) then
            alarmThreshold = alarmPtr % prevRingTime + alarmPtr % ringTimeInterval
         end if

         if (alarmThreshold <= alarmNow) then
            mpas_in_ringing_envelope = .true.
         end if
      else

         if (present(interval)) then
            alarmNow = alarmNow - interval; 
         end if

         if (alarmPtr % isRecurring) then
            alarmThreshold = alarmPtr % prevRingTime - alarmPtr % ringTimeInterval
         end if
            
         if (alarmThreshold >= alarmNow) then
            mpas_in_ringing_envelope = .true.
         end if
      end if

   end function mpas_in_ringing_envelope



   subroutine mpas_reset_clock_alarm(clock, alarmID, interval, ierr)

      implicit none

      type (MPAS_Clock_type), intent(inout) :: clock
      character (len=*), intent(in) :: alarmID
      type (MPAS_TimeInterval_type), intent(in), optional :: interval
      integer, intent(out), optional :: ierr

      type (MPAS_Time_type) :: alarmNow
      type (MPAS_Alarm_type), pointer :: alarmPtr

      if (present(ierr)) ierr = 0

      alarmPtr => clock % alarmListHead
      do while (associated(alarmPtr))
      
         if (trim(alarmPtr % alarmID) == trim(alarmID)) then

            if (mpas_in_ringing_envelope(clock, alarmPtr, interval, ierr)) then

               if (.not. alarmPtr % isRecurring) then
                  alarmPtr % isSet = .false. 
               else
                  alarmNow = mpas_get_clock_time(clock, MPAS_NOW, ierr)

                  if(clock % direction == MPAS_FORWARD) then
                     if (present(interval)) then
                        alarmNow = alarmNow + interval
                     end if

                     do while(alarmPtr % prevRingTime <= alarmNow)
                        alarmPtr % prevRingTime = alarmPtr % prevRingTime + alarmPtr % ringTimeInterval
                     end do
                     alarmPtr % prevRingTime =  alarmPtr % prevRingTime - alarmPtr % ringTimeInterval
                  else
                     if (present(interval)) then
                        alarmNow = alarmNow - interval
                     end if

                     do while(alarmPtr % prevRingTime >= alarmNow)
                        alarmPtr % prevRingTime = alarmPtr % prevRingTime - alarmPtr % ringTimeInterval
                     end do
                     alarmPtr % prevRingTime = alarmPtr % prevRingTime + alarmPtr % ringTimeInterval
                  end if
               end if
            end if
            exit
         end if
         alarmPtr => alarmPtr % next
      end do

   end subroutine mpas_reset_clock_alarm



   ! specify a valid previousRingTime for each alarm
   subroutine mpas_calibrate_alarms(clock, ierr)

      implicit none

      type (MPAS_Clock_type), intent(in) :: clock
      integer, intent(out), optional :: ierr

      type (MPAS_Time_type) :: now
      type (MPAS_Time_type) :: previousRingTime
      type (MPAS_Time_type) :: negativeNeighborRingTime
      type (MPAS_Time_type) :: positiveNeighborRingTime
      type (MPAS_Alarm_type), pointer :: alarmPtr

      now = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      
      alarmPtr => clock % alarmListHead
      do while (associated(alarmPtr))
         
         if (.not. alarmPtr % isRecurring) then
            alarmPtr % isSet = .true.            
         else
         
            previousRingTime = alarmPtr % prevRingTime

            if (previousRingTime <= now) then
            
               do while(previousRingTime <= now)
                  previousRingTime = previousRingTime + alarmPtr % ringTimeInterval
               end do
               positiveNeighborRingTime = previousRingTime
            
               do while(previousRingTime >= now)
                  previousRingTime = previousRingTime - alarmPtr % ringTimeInterval
               end do
               negativeNeighborRingTime = previousRingTime
            
            else

               do while(previousRingTime >= now)
                  previousRingTime = previousRingTime - alarmPtr % ringTimeInterval
               end do
               negativeNeighborRingTime = previousRingTime

               do while(previousRingTime <= now)
                  previousRingTime = previousRingTime + alarmPtr % ringTimeInterval
               end do
               positiveNeighborRingTime = previousRingTime
         
            end if

            if (clock % direction == MPAS_FORWARD) then
               alarmPtr % prevRingTime = negativeNeighborRingTime
            else
               alarmPtr % prevRingTime = positiveNeighborRingTime
            end if

         end if
   
         alarmPtr => alarmPtr % next
         
      end do
   
      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if
   
   end subroutine mpas_calibrate_alarms


   subroutine mpas_set_time(curr_time, YYYY, MM, DD, DoY, H, M, S, S_n, S_d, dateTimeString, ierr)

      implicit none

      type (MPAS_Time_type), intent(out) :: curr_time
      integer, intent(in), optional :: YYYY
      integer, intent(in), optional :: MM
      integer, intent(in), optional :: DD
      integer, intent(in), optional :: DoY
      integer, intent(in), optional :: H
      integer, intent(in), optional :: M
      integer, intent(in), optional :: S
      integer, intent(in), optional :: S_n
      integer, intent(in), optional :: S_d
      character (len=*), intent(in), optional :: dateTimeString
      integer, intent(out), optional :: ierr

      integer, parameter :: integerMaxDigits = 8
      integer :: year, month, day, hour, min, sec
      integer :: numerator, denominator, denominatorPower

      character (len=StrKIND) :: dateTimeString_
      character (len=StrKIND) :: dateSubString
      character (len=StrKIND) :: timeSubString
      character (len=StrKIND) :: secDecSubString
      character(len=StrKIND), pointer, dimension(:) :: subStrings

      if (present(dateTimeString)) then

         dateTimeString_ = dateTimeString
         numerator = 0
         denominator = 1

         call mpas_split_string(dateTimeString_, ".", subStrings)
         if (size(subStrings) == 2) then ! contains second decimals
            dateTimeString_ = subStrings(1)
            secDecSubString = subStrings(2)(:integerMaxDigits)
            deallocate(subStrings)
            denominatorPower = len_trim(secDecSubString)
            if(denominatorPower > 0) then
               read(secDecSubString,*) numerator 
               if(numerator > 0) then
                  denominator = 10**denominatorPower
               end if
            end if
         else if (size(subStrings) /= 1) then
            deallocate(subStrings)
            if (present(ierr)) ierr = 1
            write(stderrUnit,*) 'ERROR: Invalid DateTime string', dateTimeString
            return
         end if

         call mpas_split_string(dateTimeString_, "_", subStrings)

         if(size(subStrings) == 2) then   ! contains a date and time
            dateSubString = subStrings(1)
            timeSubString = subStrings(2)
            deallocate(subStrings)
            
            call mpas_split_string(timeSubString, ":", subStrings)
            
            if (size(subStrings) == 3) then
               read(subStrings(1),*) hour 
               read(subStrings(2),*) min 
               read(subStrings(3),*) sec 
               deallocate(subStrings)
            else
               deallocate(subStrings)
               if (present(ierr)) ierr = 1
               write(stderrUnit,*) 'ERROR: Invalid DateTime string (invalid time substring)', dateTimeString
               return
            end if

         else if(size(subStrings) == 1) then   ! contains only a date- assume all time values are 0 
            dateSubString = subStrings(1)
            deallocate(subStrings)
           
            hour = 0
            min = 0
            sec = 0
         
         else
            deallocate(subStrings)
            if (present(ierr)) ierr = 1
            write(stderrUnit,*) 'ERROR: Invalid DateTime string', dateTimeString
            return
         end if

         call mpas_split_string(dateSubString, "-", subStrings)
            
         if (size(subStrings) == 3) then
            read(subStrings(1),*) year 
            read(subStrings(2),*) month
            read(subStrings(3),*) day
            deallocate(subStrings)
         else
            deallocate(subStrings)
            if (present(ierr)) ierr = 1
            write(stderrUnit,*) 'ERROR: Invalid DateTime string (invalid date substring)', dateTimeString
            return
         end if

         call ESMF_TimeSet(curr_time % t, YY=year, MM=month, DD=day, H=hour, M=min, S=sec, Sn=numerator, Sd=denominator, rc=ierr)

      else
      
         if (present(DoY)) then
            call mpas_get_month_day(YYYY, DoY, month, day)
         
            ! consistency check
            if (present(MM)) then
               if (MM /= month) then
                  if (present(ierr)) ierr = 1
                  write(stderrUnit,*) 'ERROR: MPAS_setTime : DoY and MM are inconsistent - using DoY'
               end if
            end if
            if (present(DD)) then
               if (DD /= day) then
                  if (present(ierr)) ierr = 1
                  write(stderrUnit,*) 'ERROR: MPAS_setTime : DoY and DD are inconsistent - using DoY'
               end if
            end if
         else
            if (present(MM)) then
               month = MM
            else
               if (present(ierr)) ierr = 1
               write(stderrUnit,*) 'ERROR: MPAS_setTime : Neither DoY nor MM are specified'
               return
            end if

            if (present(DD)) then
               day = DD
            else
               if (present(ierr)) ierr = 1
               write(stderrUnit,*) 'ERROR: MPAS_setTime : Neither DoY nor DD are specified'
               return
            end if
         end if

         if (.not. isValidDate(YYYY,month,day)) then
            write(stderrUnit,*) 'ERROR: MPAS_setTime : Invalid date'
            return
         end if

         call ESMF_TimeSet(curr_time % t, YY=YYYY, MM=month, DD=day, H=H, M=M, S=S, Sn=S_n, Sd=S_d, rc=ierr)
      
      end if
      
      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if

   end subroutine mpas_set_time


   subroutine mpas_get_time(curr_time, YYYY, MM, DD, DoY, H, M, S, S_n, S_d, dateTimeString, ierr)

      implicit none

      type (MPAS_Time_type), intent(in) :: curr_time
      integer, intent(out), optional :: YYYY
      integer, intent(out), optional :: MM
      integer, intent(out), optional :: DD
      integer, intent(out), optional :: DoY
      integer, intent(out), optional :: H
      integer, intent(out), optional :: M
      integer, intent(out), optional :: S
      integer, intent(out), optional :: S_n
      integer, intent(out), optional :: S_d
      character (len=StrKIND), intent(out), optional :: dateTimeString
      integer, intent(out), optional :: ierr

      call ESMF_TimeGet(curr_time % t, YY=YYYY, MM=MM, DD=DD, H=H, M=M, S=S, Sn=S_n, Sd=S_d, rc=ierr)
      call ESMF_TimeGet(curr_time % t, dayOfYear=DoY, rc=ierr)
      call ESMF_TimeGet(curr_time % t, timeString=dateTimeString, rc=ierr)
      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if

   end subroutine mpas_get_time


   subroutine mpas_set_timeInterval(interval, YY, MM, DD, H, M, S, S_n, S_d, timeString, dt, ierr)

      implicit none

      type (MPAS_TimeInterval_type), intent(out) :: interval
      integer, intent(in), optional :: YY
      integer, intent(in), optional :: MM
      integer, intent(in), optional :: DD
      integer, intent(in), optional :: H
      integer, intent(in), optional :: M
      integer, intent(in), optional :: S
      integer, intent(in), optional :: S_n
      integer, intent(in), optional :: S_d
      character (len=*), intent(in), optional :: timeString
      real (kind=RKIND), intent(in), optional :: dt
      integer, intent(out), optional :: ierr

      integer, parameter :: integerMaxDigits = 8
!      integer :: days, hours, minutes, seconds
      integer :: numerator, denominator, denominatorPower
      type (MPAS_TimeInterval_type) :: zeroInterval

      integer :: year, month, day, hour, min, sec
      character (len=StrKIND) :: timeString_
      character (len=StrKIND) :: dateSubString
      character (len=StrKIND) :: daySubString
      character (len=StrKIND) :: monthSubString
      character (len=StrKIND) :: yearSubString
      character (len=StrKIND) :: timeSubString
      character (len=StrKIND) :: secDecSubString
      character(len=StrKIND), pointer, dimension(:) :: subStrings

!      if (present(DD)) then
!         days = DD
!      else
!         days = 0
!      end if

!      if (present(H)) then
!         hours = H
!      else
!         hours = 0
!      end if

!      if (present(M)) then
!         minutes = M
!      else
!         minutes = 0
!      end if

!      if (present(S)) then
!         seconds = S
!      else
!         seconds = 0
!      end if


      !
      ! Reduce minute count to something less than one hour
      !
!      do while (minutes > 1440)
!         days = days + 1
!         minutes = minutes - 1440
!      end do
!      do while (minutes > 60)
!         hours = hours + 1
!         minutes = minutes - 60
!      end do
!      do while (minutes < -1440)
!         days = days - 1
!         minutes = minutes + 1440
!      end do
!      do while (minutes < -60)
!         hours = hours - 1
!         minutes = minutes + 60
!      end do

      !
      ! Reduce hour count to something less than one day
      !
!      do while (hours > 24)
!         days = days + 1
!         hours = hours - 24
!      end do
!      do while (hours < -24)
!         days = days - 1
!         hours = hours + 24
!      end do

      !
      ! Any leftover minutes and hours are given to the second count
      !
!      seconds = seconds + hours*3600 + minutes*60

!      call ESMF_TimeIntervalSet(interval % ti, D=days, S=seconds, Sn=S_n, Sd=S_d, rc=ierr)


      if (present(timeString) .or. present(dt)) then


         if(present(dt)) then
            write (timeString_,*) "00:00:", dt         
         else
            timeString_ = timeString
         end if

         numerator = 0
         denominator = 1

         call mpas_split_string(timeString_, ".", subStrings)
         
         if (size(subStrings) == 2) then ! contains second decimals
            timeString_ = subStrings(1)
            secDecSubString = subStrings(2)(:integerMaxDigits)
            deallocate(subStrings)

            denominatorPower = len_trim(secDecSubString)
            if(denominatorPower > 0) then
               read(secDecSubString,*) numerator 
               if(numerator > 0) then
                  denominator = 10**denominatorPower
               end if
            end if
         else if (size(subStrings) /= 1) then
            deallocate(subStrings)
            if (present(ierr)) ierr = 1
            write(stderrUnit,*) 'ERROR: Invalid TimeInterval string ', trim(timeString)
            return
         end if

         call mpas_split_string(timeString_, "_", subStrings)

         if(size(subStrings) == 2) then   ! contains a date and time
            dateSubString = subStrings(1)
            timeSubString = subStrings(2)
            deallocate(subStrings)

            call mpas_split_string(dateSubString, "-", subStrings)

            if(size(subStrings) == 3) then ! Contains year, month, and day
               read(subStrings(1), *), year
               read(subStrings(2), *), month
               read(subStrings(3), *), day
            else if(size(subStrings) == 2) then ! Contains month and day
               year = 0
               read(subStrings(1), *), month
               read(subStrings(2), *), day
            else if(size(subStrings) == 1) then ! Contains day
               year = 0
               month = 0
               read(subStrings(1), *), day
            else ! Error?
               year = 0
               month = 0
               day = 0
               !write(stderrUnit,*) 'ERROR: Invalid TimeInterval string ', trim(timeString)
            end if

            deallocate(subStrings)
         else if(size(subStrings) == 1) then   ! contains only a time- assume year, month, and day are 0
            timeSubString = subStrings(1)
            deallocate(subStrings)
            year = 0
            month = 0
            day = 0
         else
            deallocate(subStrings)
            if (present(ierr)) ierr = 1
            write(stderrUnit,*) 'ERROR: Invalid TimeInterval string ', trim(timeString)
            return
         end if

         call mpas_split_string(timeSubString, ":", subStrings)
            
         if (size(subStrings) == 3) then
            read(subStrings(1),*) hour 
            read(subStrings(2),*) min 
            read(subStrings(3),*) sec 
            deallocate(subStrings)
         else
            deallocate(subStrings)
            if (present(ierr)) ierr = 1
            write(stderrUnit,*) 'ERROR: Invalid TimeInterval string (invalid time substring) ', trim(timeString)
            return
         end if

         call ESMF_TimeIntervalSet(interval % ti, YY=year, MM=month, D=day, H=hour, M=min, S=sec, Sn=numerator, Sd=denominator, rc=ierr)

      else

         call ESMF_TimeIntervalSet(interval % ti, YY=YY, MM=MM, D=DD, H=H, M=M, S=S, Sn=S_n, Sd=S_d, rc=ierr)
      
      end if

!     ! verify that time interval is positive
!     call ESMF_TimeIntervalSet(zeroInterval % ti, D=0, H=0, M=0, S=0, rc=ierr)

!     if (present(ierr)) then
!        if (ierr == ESMF_SUCCESS) ierr = 0
!     end if

!     if (interval <= zeroInterval) then
!        if (present(ierr)) ierr = 1   
!        write(stderrUnit,*) 'ERROR: TimeInterval must be greater than zero: ', trim(timeString) !'ERROR: TimeInterval cannot be negative'
!     end if
      
   end subroutine mpas_set_timeInterval


   subroutine mpas_get_timeInterval(interval, StartTimeIn, DD, H, M, S, S_n, S_d, timeString, dt, ierr)
! TODO: add double-precision seconds

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: interval
      type (MPAS_Time_type), intent(in), optional :: StartTimeIn
      ! For time intervals that require months and/or years, ESMF needs to know the start
      ! time to get a time interval in any format besides the string format.
      integer, intent(out), optional :: DD
      integer, intent(out), optional :: H
      integer, intent(out), optional :: M
      integer, intent(out), optional :: S
      integer, intent(out), optional :: S_n
      integer, intent(out), optional :: S_d
      character (len=StrKIND), intent(out), optional :: timeString
      real (kind=RKIND), intent(out), optional :: dt
      integer, intent(out), optional :: ierr

      integer :: days, seconds, sn, sd



      if (present(StartTimeIn)) then
         call ESMF_TimeIntervalGet(interval % ti, StartTimeIn=StartTimeIn%t, D=days, S=seconds, Sn=sn, Sd=sd, rc=ierr)
      else
          if ( interval % ti % YR /= 0 .or. interval % ti % MM /= 0 ) then
             if (present(ierr)) ierr = 1
             write(stderrUnit,*) 'ERROR: mpas_get_timeInterval cannnot return time interval information for an interval containing months and years without a startTimeIn argument.'
             return
          end if
          call ESMF_TimeIntervalGet(interval % ti, D=days, S=seconds, Sn=sn, Sd=sd, rc=ierr)
      endif

      if (present(dt)) then
         dt = (days * 24 * 60 * 60) + seconds + (sn / sd)
      end if

      if (present(DD)) then
         DD = days
         days = 0
      end if

      if (present(H)) then
         H = (seconds - mod(seconds,3600)) / 3600
         seconds = seconds - H*3600
         H = H + days * 24
         days = 0
      end if

      if (present(M)) then
         M = (seconds - mod(seconds,60)) / 60
         seconds = seconds - M*60
         M = M + days * 1440
         days = 0
      end if

      if (present(S)) then
         S = seconds
      end if

      if (present(S_n)) then
         S_n = sn
      end if

      if (present(S_d)) then
         if (sd == 0) then   ! may only occur if (sn == 0)?
            sd = 1
         end if
         S_d = sd
      end if

      if (present(timeString)) then
         call ESMF_TimeIntervalGet(interval % ti, timeString=timeString, rc=ierr)
      end if

      if (present(ierr)) then
         if (ierr == ESMF_SUCCESS) ierr = 0
      end if

   end subroutine mpas_get_timeInterval


   type (MPAS_Time_type) function add_t_ti(t, ti)

      implicit none

      type (MPAS_Time_type), intent(in) :: t
      type (MPAS_TimeInterval_type), intent(in) :: ti

      add_t_ti % t = t % t + ti % ti

   end function add_t_ti


   type (MPAS_TimeInterval_type) function add_ti_ti(ti1, ti2)

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: ti1, ti2

      add_ti_ti % ti = ti1 % ti + ti2 % ti

   end function add_ti_ti


   type (MPAS_TimeInterval_type) function sub_t_t(t1, t2)

      implicit none

      type (MPAS_Time_type), intent(in) :: t1, t2

      sub_t_t % ti = t1 % t - t2 % t

   end function sub_t_t


   type (MPAS_Time_type) function sub_t_ti(t, ti)

      implicit none

      type (MPAS_Time_type), intent(in) :: t
      type (MPAS_TimeInterval_type), intent(in) :: ti

      sub_t_ti % t = t % t - ti % ti

   end function sub_t_ti


   type (MPAS_TimeInterval_type) function sub_ti_ti(ti1, ti2)

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: ti1, ti2

      sub_ti_ti % ti = ti1 % ti - ti2 % ti

   end function sub_ti_ti


   type (MPAS_TimeInterval_type) function mul_ti_n(ti, n)

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: ti
      integer, intent(in) :: n

      mul_ti_n % ti = ti % ti * n

   end function mul_ti_n


   type (MPAS_TimeInterval_type) function div_ti_n(ti, n)

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: ti
      integer, intent(in) :: n

      div_ti_n % ti = ti % ti / n

   end function div_ti_n

   !-----------------------------------------------------------------------
   !  routine mpas_interval_division
   !
   !> \brief This routine computes the number intervals that fit into another interval.
   !> \author Michael Duda, Doug Jacobsen
   !> \date   10/02/2014
   !> \details This routine is a wrapper to two different methods of computing
   !> the number of intervals that fit into another interval.
   !>
   !-----------------------------------------------------------------------
   subroutine mpas_interval_division(ref_time, num, den, n, rem)

      implicit none

      type (MPAS_Time_type), intent(in) :: ref_time
      type (MPAS_TimeInterval_type), intent(in) :: num
      type (MPAS_TimeInterval_type), intent(in) :: den
      integer, intent(out) :: n
      type (MPAS_TimeInterval_type), intent(out) :: rem

      ! If either time interval contains years or months, need to use slower linear interval computation
      ! Otherwise, we can use faster log interval computation
      if ( num % ti % YR /= 0 .or. num % ti % MM /= 0 .or. den % ti % YR /= 0 .or. den % ti % MM /= 0 ) then
         call mpas_interval_division_linear(ref_time, num, den, n, rem)
      else
         call mpas_interval_division_log(num, den, n, rem)
      end if

   end subroutine mpas_interval_division

   !-----------------------------------------------------------------------
   !  routine mpas_interval_division_log
   !
   !> \brief This routine computes the number intervals that fit into another interval using a log search.
   !> \author Michael Duda, Doug Jacobsen
   !> \date   10/02/2014
   !> \details This routine computes the number of intervals that fit into
   !>   another time interval using a log search. It is preferred over the
   !>   _linear alternative, but only works when the intervals are in terms of days
   !>   or smaller.
   !>
   !-----------------------------------------------------------------------
   subroutine mpas_interval_division_log(num, den, n, rem)

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: num
      type (MPAS_TimeInterval_type), intent(in) :: den
      integer, intent(out) :: n
      type (MPAS_TimeInterval_type), intent(out) :: rem

      type (MPAS_TimeInterval_type) :: temp
      type (MPAS_TimeInterval_type) :: zero
      integer :: nn

      call mpas_set_timeInterval(zero, S=0)

      !
      ! If the numerator is smaller than the denominator, just return the numerator as the remainder
      !
      if (num < den) then
         n = 0
         rem = num
         return
      end if


      !
      ! Avoid division by zero
      !
      if (den == zero) then
         write(stderrUnit,*) 'Error: Attempting to divide by zero.\n'
         n = 0
         rem = zero
         return
      end if


      !
      ! Begin by finding the smallest multiple of the denominator that is at least as large as the numerator and is also a power of two
      !
      temp = den
      nn = 1
      do while (temp <= num)
         temp = temp * 2
         nn = nn * 2
      end do

      !
      ! Dividing by two, we're guaranteed that temp is at most the value of the numerator
      !
      temp = temp / 2
      nn = nn / 2

      !
      ! Work backwards to zero
      !
      n = 0
      rem = num
      do while (nn > 0)
         if (temp <= rem) then
            rem = rem - temp
            n = n + nn
         end if
         nn = nn / 2
         temp = temp / 2
      end do

   end subroutine mpas_interval_division_log


   !-----------------------------------------------------------------------
   !  routine mpas_interval_division_linear
   !
   !> \brief This routine computes the number intervals that fit into another interval using a linear search.
   !> \author Michael Duda, Doug Jacobsen
   !> \date   10/02/2014
   !> \details This routine computes the number of intervals that fit into
   !>   another time interval using a linear search. It is slower than the _log
   !>   alternative, but works when intervals contain months or longer interval
   !>   sections.
   !>
   !-----------------------------------------------------------------------
   subroutine mpas_interval_division_linear(ref_time, num, den, n, rem)

      implicit none

      type (MPAS_Time_type), intent(in) :: ref_time
      type (MPAS_TimeInterval_type), intent(in) :: num
      type (MPAS_TimeInterval_type), intent(in) :: den
      integer, intent(out) :: n
      type (MPAS_TimeInterval_type), intent(out) :: rem

      type (MPAS_Time_type) :: target_time
      type (MPAS_Time_type) :: updated_time

      type (MPAS_TimeInterval_type) :: temp
      type (MPAS_TimeInterval_type) :: zero

      target_time = ref_time + num

      updated_time = ref_time + den

      n = 0

      ! If the denominator is larger than the numerator, return 0 intervals,
      ! and the numerator as the remainder
      if ( target_time < updated_time ) then
         rem = num
         return
      end if

      ! One interval of den already fits into num
      n = n + 1

      do while (target_time > updated_time)
         updated_time = updated_time + den
         n = n + 1
      end do

      ! Here, if updated_time is larger than target time. Need to subtract den once, and compute remainder
      if ( updated_time > target_time ) then
         updated_time = updated_time - den
         n = n - 1
         rem = target_time - updated_time
      else
         call mpas_set_timeInterval(rem, S=0)
      end if

      return
   end subroutine mpas_interval_division_linear


   logical function eq_t_t(t1, t2)

      implicit none

      type (MPAS_Time_type), intent(in) :: t1, t2

      eq_t_t = (t1 % t == t2 % t)

   end function eq_t_t


   logical function ne_t_t(t1, t2)

      implicit none

      type (MPAS_Time_type), intent(in) :: t1, t2

      ne_t_t = (t1 % t /= t2 % t)

   end function ne_t_t


   logical function lt_t_t(t1, t2)

      implicit none

      type (MPAS_Time_type), intent(in) :: t1, t2

      lt_t_t = (t1 % t < t2 % t)

   end function lt_t_t


   logical function gt_t_t(t1, t2)

      implicit none

      type (MPAS_Time_type), intent(in) :: t1, t2

      gt_t_t = (t1 % t > t2 % t)

   end function gt_t_t


   logical function le_t_t(t1, t2)

      implicit none

      type (MPAS_Time_type), intent(in) :: t1, t2

      le_t_t = (t1 % t <= t2 % t)

   end function le_t_t


   logical function ge_t_t(t1, t2)

      implicit none

      type (MPAS_Time_type), intent(in) :: t1, t2

      ge_t_t = (t1 % t >= t2 % t)

   end function ge_t_t


   logical function eq_ti_ti(ti1, ti2)

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: ti1, ti2

      eq_ti_ti = (ti1 % ti == ti2 % ti)

   end function eq_ti_ti


   logical function ne_ti_ti(ti1, ti2)

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: ti1, ti2

      ne_ti_ti = (ti1 % ti /= ti2 % ti)

   end function ne_ti_ti


   logical function lt_ti_ti(ti1, ti2)

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: ti1, ti2

      lt_ti_ti = (ti1 % ti < ti2 % ti)

   end function lt_ti_ti


   logical function gt_ti_ti(ti1, ti2)

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: ti1, ti2

      gt_ti_ti = (ti1 % ti > ti2 % ti)

   end function gt_ti_ti


   logical function le_ti_ti(ti1, ti2)

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: ti1, ti2

      le_ti_ti = (ti1 % ti <= ti2 % ti)

   end function le_ti_ti


   logical function ge_ti_ti(ti1, ti2)

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: ti1, ti2

      ge_ti_ti = (ti1 % ti >= ti2 % ti)

   end function ge_ti_ti


   type (MPAS_TimeInterval_type) function neg_ti(ti)

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: ti

      integer :: rc
      integer :: D, S, Sn, Sd

      call ESMF_TimeIntervalGet(ti % ti, D=D, S=S, Sn=Sn, Sd=Sd, rc=rc)
      D    = -D 
      S    = -S 
      Sn   = -Sn
      call ESMF_TimeIntervalSet(neg_ti % ti, D=D, S=S, Sn=Sn, Sd=Sd, rc=rc)

   end function neg_ti


   type (MPAS_TimeInterval_type) function abs_ti(ti)

      implicit none

      type (MPAS_TimeInterval_type), intent(in) :: ti

      type (MPAS_TimeInterval_type) :: zeroInterval
      integer :: rc
      integer :: D, S, Sn, Sd

      call ESMF_TimeIntervalSet(zeroInterval % ti, D=0, H=0, M=0, S=0, rc=rc)

      if(ti < zeroInterval) then
         call ESMF_TimeIntervalGet(ti % ti, D=D, S=S, Sn=Sn, Sd=Sd, rc=rc)
         D    = -D 
         S    = -S 
         Sn   = -Sn
         call ESMF_TimeIntervalSet(abs_ti % ti, D=D, S=S, Sn=Sn, Sd=Sd, rc=rc)
      else
         abs_ti = ti
      end if

   end function abs_ti


! TODO: Implement this function
!   type (MPAS_TimeInterval_type) function mod(ti1, ti2)
!
!      implicit none
!
!      type (MPAS_TimeInterval_type), intent(in) :: ti1, ti2
!
!      mod % ti = mod(ti1 % ti, ti2 % ti)
!
!   end function mod


   subroutine mpas_split_string(string, delimiter, subStrings)   
      
      implicit none
      
      character(len=*), intent(in) :: string
      character, intent(in) :: delimiter
      character(len=*), pointer, dimension(:) :: subStrings
      
      integer :: i, start, index

      index = 1
      do i = 1, len(string)
         if(string(i:i) == delimiter) then
            index = index + 1
         end if
      end do

      allocate(subStrings(1:index))

      start = 1
      index = 1
      do i = 1, len(string)
         if(string(i:i) == delimiter) then
               subStrings(index) = string(start:i-1) 
               index = index + 1
               start = i + 1
         end if
      end do
      subStrings(index) = string(start:len(string)) 
      
   end subroutine mpas_split_string


    subroutine mpas_get_month_day(YYYY, DoY, month, day)
       
       implicit none

       integer, intent(in) :: YYYY, DoY
       integer, intent(out) :: month, day

       integer, dimension(12) :: dpm
       
       if (isLeapYear(YYYY)) then
          dpm(:) = daysInMonthLeap
       else
          dpm(:) = daysInMonth
       end if

       month = 1
       day = DoY
       do while (day > dpm(month))
          day = day -  dpm(month)
          month = month + 1       
       end do

    end subroutine mpas_get_month_day


   logical function isValidDate(YYYY, MM, DD)
   
      integer, intent(in) :: YYYY, MM, DD
      integer :: daysInMM
      
      isValidDate = .true.

      ! TODO: ???? Gregorian calendar has no year zero, but perhaps 0 = 1 BC ??? 
      !if (YYYY == 0) then
      !   isValidDate = .false.
      !   return
      !end if

      if (MM < 1 .or. MM > 12) then
         isValidDate = .false.
         return
      end if

      if (DD < 1) then
         isValidDate = .false.
         return
      end if

      if(TheCalendar == MPAS_360DAY) then
         daysInMM = 30
      else
         if (TheCalendar == MPAS_GREGORIAN .and. isLeapYear(YYYY)) then
            daysInMM = daysInMonthLeap(MM)
         else
            daysInMM = daysInMonth(MM)        
         end if
      end if
     
      if (DD > daysInMM) then
         isValidDate = .false.
         return
      end if

   end function

    
    logical function isLeapYear(year)

       implicit none

       integer, intent(in) :: year

       isLeapYear = .false.
       
       if (mod(year,4) == 0) then
          if (mod(year,100) == 0) then
             if (mod(year,400) == 0) then
                isLeapYear = .true.
             end if
          else
             isLeapYear = .true.
          end if
       end if

    end function isLeapYear

    !-----------------------------------------------------------------------
    !  routine mpas_expand_string
    !
    !> \brief This is a utility routine that expands a string with a timestamp.
    !> \author Michael Duda, Doug Jacobsen
    !> \date   07/23/2014
    !> \details This routine will take a time stamp, and a string as
    !>   input, and expand the string according to the time stamp provided.
    !>   $Y -> year
    !>   $M -> month
    !>   $D -> day
    !>   $d -> day of year
    !>   $h -> hour
    !>   $m -> minute
    !>   $s -> second
    !>   $g -> multi-grid level
    !-----------------------------------------------------------------------
    subroutine mpas_expand_string(timeStamp, inString, outString)!{{{

        implicit none

        character (len=*), intent(in) :: timeStamp
        character (len=*), intent(in) :: inString
        character (len=StrKIND), intent(out) :: outString

        type (MPAS_Time_type) :: curTime

        integer :: i, curLen
        integer :: year, month, day, hour, minute, second, DoY

        character (len=ShortStrKIND) :: timePart
        character (len=ShortStrKIND) :: yearFormat
        logical :: charExpand

        call mpas_set_time(curTime, dateTimeString=timeStamp)

        call mpas_get_time(curTime, YYYY=year)

        write(yearFormat, '(a,i10,a)') '(i0.',yearWidth,')'

        write(outString,*) ''
        write(timePart,*) ''

        curLen = 0
        charExpand = .false.
        do i = 1, len_trim(inString)
           if (inString(i:i) == '$' ) then
               charExpand = .true.
           else if (inString(i:i) /= '$') then
               if (charExpand) then
                  select case (inString(i:i))
                     case ('Y')
                         call mpas_get_time(curTime, YYYY=year)
                         write(timePart, yearFormat) year
                         outString = trim(outString) // trim(timePart)
                     case ('M')
                         call mpas_get_time(curTime, MM=month)
                         write(timePart, '(i0.2)') month
                         outString = trim(outString) // trim(timePart)
                     case ('D')
                         call mpas_get_time(curTime, DD=day)
                         write(timePart, '(i0.2)') day
                         outString = trim(outString) // trim(timePart)
                     case ('d')
                         call mpas_get_time(curTime, DoY=DoY)
                         write(timePart, '(i0.3)') DoY
                         outString = trim(outString) // trim(timePart)
                     case ('h')
                         call mpas_get_time(curTime, H=hour)
                         write(timePart, '(i0.2)') hour
                         outString = trim(outString) // trim(timePart)
                     case ('m')
                         call mpas_get_time(curTime, M=minute)
                         write(timePart, '(i0.2)') minute
                         outString = trim(outString) // trim(timePart)
                     case ('s')
                         call mpas_get_time(curTime, S=second)
                         write(timePart, '(i0.2)') second
                         outString = trim(outString) // trim(timePart)
!                    case ('G')
                        ! Expands to multi-grid level
                     case default
                        write(stderrUnit, *) 'ERROR: mpas_expand_string option $', inString(i:i), ' is not a valid expansion character.'
                        call mpas_dmpar_global_abort('ERROR: mpas_timekeeping')
                  end select

                  curLen = len_trim(outString)
                  charExpand = .false.
               else
                  outString(curLen+1:curLen+1) = inString(i:i)
                  curLen = curLen+1
               end if
           else
           end if
        end do

    end subroutine mpas_expand_string!}}}
 



end module mpas_timekeeping


subroutine wrf_error_fatal(msg)

   implicit none

   character (len=*) :: msg

   call mpas_dmpar_global_abort('ERROR: mpas_timekeeping: '//trim(msg))

end subroutine wrf_error_fatal
