! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!  ocn_time_averages
!
!> \brief MPAS ocean analysis core member: time_averages
!> \author Jon Woodring
!> \date   March 2, 2015
!> \details
!>  Flexible time averaging, mins, and maxes of fields.
!-----------------------------------------------------------------------
module ocn_time_averages
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
   public :: ocn_init_time_averages, &
             ocn_compute_time_averages, &
             ocn_restart_time_averages, &
             ocn_finalize_time_averages

   ! Private module variables
   !--------------------------------------------------------------------

   ! startup, interval, and restart is done in the outer analysis driver

   ! time buffer type
   ! this keeps track of timers and if and when they need to accumulate
   type time_buffer_type
      ! internal state
      logical :: started_flag, accumulate_flag, reset_flag, delay_reset_flag
      logical :: reset_alarm_armed
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
!  routine walk_string
!
!> \brief   Walk a space delimited string to find substrings
!> \author  Jon Woodring
!> \date    March 2, 2015
!> \details 
!>  Walk a string delimited by spaces and return the first substring
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
 
      ! find the first substring that isn't whitespace
      i = verify(next, ' ')
      ok = i > 0
      ! if we can't find one, stop
      if (.not. ok) then
        return
      end if

      ! make a new string and find the first whitespace
      copy = next(i:)
      i = scan(copy, ' ')
      
      ! return that substring and the remainder
      substr = copy(1:i-1)
      next = copy(i:)
 
   end subroutine walk_string!}}}

!***********************************************************************
!  routine set_times
!
!> \brief   Set a list of times
!> \author  Jon Woodring
!> \date    March 2, 2015
!> \details 
!>  Walk a list of times delimited by spaces and set the time info
!>  for the buffer structure so that alarms can be set.
!-----------------------------------------------------------------------
   subroutine set_times(buffers, number_of_buffers, clock, &
                        which, config_str, inv_str, ok, err)
      ! input variables
      !-----------------------------------------------------------------
      integer, intent(in) :: number_of_buffers, which
      character (len=StrKIND), pointer, intent(in) :: config_str
      character (len=StrKIND), intent(in) :: inv_str

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
      b = 1
      call walk_string(next_str, time_str, ok)
      ! while the time string is ok
      do while (ok)
         ! exit if we went over
         if (b .gt. number_of_buffers) then
           exit
         end if

         ! set the time
         if (which .eq. START_TIMES) then
            if (time_str .eq. 'same_as_simulation') then
               buffers(b) % start_time = mpas_get_clock_time(clock, &
                  MPAS_NOW, err)
            else
               call mpas_set_time(buffers(b) % start_time, &
                    dateTimeString=time_str, ierr=err)
            end if
         else if (which .eq. DURATION_INTERVALS) then
            if (time_str .eq. 'same_as_repeat') then
               buffers(b) % duration_interval = buffers(b) % repeat_interval
            else
               call mpas_set_timeInterval(buffers(b) % duration_interval, &
                    timeString=time_str, ierr=err)
            end if
         else if (which .eq. REPEAT_INTERVALS) then
            if (time_str .eq. 'same_as_reset') then
               buffers(b) % repeat_interval = buffers(b) % reset_interval
            else
               call mpas_set_timeInterval(buffers(b) % repeat_interval, &
                    timeString=time_str, ierr=err)
            end if
         else
            if (time_str .eq. 'same_as_output') then
               buffers(b) % reset_alarm_armed = .false.
            else
               call mpas_set_timeInterval(buffers(b) % reset_interval, &
                    timeString=time_str, ierr=err)
               buffers(b) % reset_alarm_armed = .true.
            end if
         end if 
         ! get the next time string
         b = b + 1
         call walk_string(next_str, time_str, ok)
      end do

      ! only ok if we parsed out as many as there are number of buffers
      ok = number_of_buffers .eq. (b - 1)
    end subroutine set_times



!***********************************************************************
!  routine add_new_field
!
!> \brief   Function to create a new field from an existing field
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
      type (field0DReal), pointer :: r0i, or0
      type (field1DReal), pointer :: r1i, or1
      type (field2DReal), pointer :: r2i, or2
      type (field3DReal), pointer :: r3i, or3
      type (field4DReal), pointer :: r4i, or4
      type (field5DReal), pointer :: r5i, or5
      type (field0DInteger), pointer :: i0i, oi0
      type (field1DInteger), pointer :: i1i, oi1
      type (field2DInteger), pointer :: i2i, oi2
      type (field3DInteger), pointer :: i3i, oi3
      integer :: i

      ! start procedure
      !-----------------------------------------------------------------

! macro
#define COPY_FIELDS(SRC, DST) \
call mpas_pool_get_field(pool, inname, SRC, 1) ;\
call mpas_duplicate_field(SRC, DST) ;\
DST % fieldName = trim(prefix) // DST % fieldName ;\
if (DST % isVarArray) then ;\
   do i = 1, size(DST % constituentNames) ;\
      DST % constituentNames(i) = trim(prefix) // \
         DST % constituentNames(i) ;\
   end do ;\
end if ;\
call mpas_pool_add_field(pool, DST % fieldName, DST)
! end macro

      ! duplicate field and add new field to pool
      if (info % fieldType .eq. MPAS_POOL_REAL) then
         if (info % nDims .eq. 0) then
            COPY_FIELDS(r0i, or0)
         else if (info % nDims .eq. 1) then
            COPY_FIELDS(r1i, or1)
         else if (info % nDims .eq. 2) then
            COPY_FIELDS(r2i, or2)
         else if (info % nDims .eq. 3) then
            COPY_FIELDS(r3i, or3)
         else if (info % nDims .eq. 4) then
            COPY_FIELDS(r4i, or4)
         else
            COPY_FIELDS(r5i, or5)
         end if
      else
         if (info % nDims .eq. 0) then
            COPY_FIELDS(i0i, oi0)
         else if (info % nDims .eq. 1) then
            COPY_FIELDS(i1i, oi1)
         else if (info % nDims .eq. 2) then
            COPY_FIELDS(i2i, oi2)
         else
            COPY_FIELDS(i3i, oi3)
         end if
      end if

   end subroutine add_new_field!}}}

#undef COPY_FIELDS



!***********************************************************************
!  routine ocn_init_time_averages
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 2, 2015
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
   subroutine ocn_init_time_averages(domain, err)!{{{

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
      integer, pointer :: number_of_buffers
      integer :: number_of_variables
      character (len=StrKIND) :: stream_str, prefix_str, &
         config_str, buffer_str, op_str, var_str, field, inv_time
      logical :: ok

      ! start procedure
      !-----------------------------------------------------------------
      err = 0

      ! string representation
      ! TODO placeholder for some unique ID if this code is replicated
      !      per multiple AMs for multiple streams
      stream_str = ''
      prefix_str = 'config_AM_timeAverages' // trim(stream_str)

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
         call mpas_dmpar_global_abort('Error: unknown operation in time averaging analysis member configuration.')
      end if

      ! get the number of individual buffers and set up timers
      config_str = trim(prefix_str) // '_number_of_buffers' 
      call mpas_pool_get_config(domain % configs, config_str, number_of_buffers)

      ! assert number_of_buffers > 0
      if (number_of_buffers .lt. 1) then
         call mpas_dmpar_global_abort('Error: number of buffers < 0 in time averaging analysis member configuration.')
      end if

      ! get the stream name
      config_str = trim(prefix_str) // '_stream_name'
      call mpas_pool_get_config(domain % configs, config_str, stream_name)

      if (stream_name .eq. 'none') then
         call mpas_dmpar_global_abort('Error: stream cannot be "none" for time averages.')
      end if

      ! set up all of the timing
      !
      ! order matters, don't reorder these!
      ! it matters because times/intervals can be configured to be equal

      ! allocate the state for the buffers
      allocate(buffers(number_of_buffers))

      ! get the interval time
      call mpas_stream_mgr_get_property(domain % streamManager, &
           stream_name, MPAS_STREAM_PROPERTY_FILENAME_INTV, &
           inv_time, err)

      ! configure reset intervals
      config_str = trim(prefix_str) // '_reset_intervals'
      call mpas_pool_get_config(domain % configs, config_str, config_results)
      call set_times(buffers, number_of_buffers, domain % clock, &
                     RESET_INTERVALS, config_results, inv_time, ok, err) 
      if (.not. ok) then
         call mpas_dmpar_global_abort('Error: number_of_buffers != number of reset_intervals in time averaging member analysis member configuration.')
      end if
      ! TODO? sanity check that it is <= output_interval

      ! configure repeat intervals
      config_str = trim(prefix_str) // '_repeat_intervals'
      call mpas_pool_get_config(domain % configs, config_str, config_results)
      call set_times(buffers, number_of_buffers, domain % clock, &
                     REPEAT_INTERVALS, config_results, inv_time, ok, err) 
      if (.not. ok) then
         call mpas_dmpar_global_abort('Error: number_of_buffers != number of repeat_intervals in time averaging member analysis member configuration.')
      end if

      ! configure duration intervals
      config_str = trim(prefix_str) // '_duration_intervals'
      call mpas_pool_get_config(domain % configs, config_str, config_results)
      call set_times(buffers, number_of_buffers, domain % clock, &
                     DURATION_INTERVALS, config_results, inv_time, ok, err) 
      if (.not. ok) then
         call mpas_dmpar_global_abort('Error: number_of_buffers != number of duration_intervals in time averaging member analysis member configuration.')
      end if
      ! TODO sanity check to see if >= compute_interval * 2

      ! configure start times
      config_str = trim(prefix_str) // '_initial_times' 
      call mpas_pool_get_config(domain % configs, config_str, config_results)
      call set_times(buffers, number_of_buffers, domain % clock, &
                     START_TIMES, config_results, inv_time, ok, err) 
      if (.not. ok) then
         call mpas_dmpar_global_abort('Error: number_of_buffers != number of initial_times in time averaging member analysis member configuration.')
      end if

      ! check if the configuration is sensible
      do b = 1, number_of_buffers
         if (buffers(b) % repeat_interval .gt. &
             buffers(b) % reset_interval) then
            ! TODO error out
            write(stderrUnit,*) 'Warning: repeat_interval > reset_interval in time averaging analysis member configuration. Truncating repeat_interval.'
            buffers(b) % repeat_interval = buffers(b) % reset_interval
         end if

         if (buffers(b) % duration_interval .gt. &
             buffers(b) % repeat_interval) then
            ! TODO error out
            write(stderrUnit,*) 'Warning: duration_interval > repeat_interval in time averaging analysis member configuration. Truncating duration_interval.'
            buffers(b) % repeat_interval = buffers(b) % reset_interval
         end if
      end do

      ! 
      ! OK, if we got this far, then we should be able to allocate memory
      ! and set up the timers and variables that we will average
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

      !
      ! TODO How to add mesh to stream?
      !
      !! optionally add mesh to stream
      !call mpas_pool_get_config(domain % configs, &
      !   'config_time_averages_copy_mesh', copy_mesh)
      !if (copy_mesh) then
      !   call mpas_stream_mgr_add_stream_fields(manager, stream_name, &
      !      'mesh', err)
      !end if

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
            call mpas_dmpar_global_abort('Error: a field listed in the output stream is not real or integer in time averaging analysis member configuration.')
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
            buffers(b) % start_time, &
            buffers(b) % repeat_interval, ierr=err)

         buffers(b) % duration_alarm_ID = &
            'tavg_duration' // trim(stream_str) // '_' // buffer_str
         call mpas_add_clock_alarm(domain % clock, &
            buffers(b) % duration_alarm_ID, &
            buffers(b) % start_time + &
            buffers(b) % duration_interval, & 
            buffers(b) % repeat_interval, ierr=err)

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
         buffers(b) % delay_reset_flag = .false.
      end do

   end subroutine ocn_init_time_averages!}}}



!***********************************************************************
!  routine timer_checking 
!
!> \brief   Timer functions to determine when to run
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
            ! TODO only reset if not restart
            buffers(b) % reset_flag = .true.
         end if

         ! if we aren't started, continue to next buffer
         if (.not. buffers(b) % started_flag) then
            continue
         end if

         ! check various other alarms
         ! see if we need to reset
         if (buffers(b) % reset_alarm_armed) then
            if(mpas_is_alarm_ringing(domain % clock, &
               buffers(b) % reset_alarm_ID, ierr=err)) then
               call mpas_reset_clock_alarm(domain % clock, &
                  buffers(b) % reset_alarm_ID, ierr=err)
               buffers(b) % reset_flag = .true.
            end if
         else
            if(mpas_stream_mgr_ringing_alarms(domain % streamManager, &
               stream_name)) then
               buffers(b) % reset_flag = .true. 
            end if
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
         !  overlaps on the same timer)
         if (mpas_is_alarm_ringing(domain % clock, &
               buffers(b) % repeat_alarm_ID, ierr=err)) then
            call mpas_reset_clock_alarm(domain % clock, &
               buffers(b) % repeat_alarm_ID, ierr=err)
            buffers(b) % accumulate_flag = .true.
         end if

         ! see if we need to delay resetting to the next time around,
         ! when the reset is on the same time as an output
         ! (we don't want to clear if we are outputting)
         if (buffers(b) % reset_flag) then
            if(mpas_stream_mgr_ringing_alarms(domain % streamManager, &
               stream_name)) then
               buffers(b) % delay_reset_flag = .true. 
            end if
         end if
      end do

   end subroutine timer_checking!}}}

!***********************************************************************
!  macro OPERATE
!
!> \brief   A macro to support operations on different run-time types
!> \author  Jon Woodring
!> \date    March 2, 2015
!> \details 
!>  This macro encapsulates the different opertions that can occur
!>  based on the run-time types. It is written as a macro to cut down
!>  on copy-pasting and duplication errors. (This would likely be
!>  instantiated generics/templates in other languages.)
!-----------------------------------------------------------------------

! first half of the macro
#define FIRST_HALF(SUBNAME) \
subroutine operate ## SUBNAME (start_block, tvar) ;\
type (block_type), pointer, intent(in) :: start_block ;\
\
type (time_variable_type), intent(inout) :: tvar ;\
\
integer :: b ;\
type (block_type), pointer :: block ;

! second half of the macro
#define SECOND_HALF \
block => start_block ;\
do while (associated(block)) ;\
   call mpas_pool_get_array(block % allFields, \
      tvar % input_name, in_array, 1) ;\
\
   do b = 1, size(buffers) ;\
      if (buffers(b) % reset_flag .and. \
          (.not. buffers(b) % delay_reset_flag)) then ;\
         call mpas_pool_get_array(block % allFields, \
              tvar % output_names(b), out_array, 1) ;\
         out_array = in_array ;\
      else if (buffers(b) % accumulate_flag) then ;\
         call mpas_pool_get_array(block % allFields, \
              tvar % output_names(b), out_array, 1) ;

! averaging is done by multiplying out and dividing such that
! the average state is always in a normalized form -- while
! this could (will) cause more error in the long run, it does
! mean that other AMs will be able to use this data and it will
! always be prenormalized (it also means that we don't have to
! have a special case of normalizing the data before writing it
! to disk)
#define AVG_MAC \
         out_array = (out_array * \
                     (buffers(b) % total_accum - 1) + in_array) \
                     / buffers(b) % total_accum ;

#define MIN_MAC \
         out_array = min(out_array, in_array) ;

#define MAX_MAC \
         out_array = max(out_array, in_array) ;

#define END_CAP(SUBNAME) \
      end if ;\
   end do ;\
\
   block => block % next ;\
end do ;\
\
end subroutine operate ## SUBNAME ;

! had to create this as a two part macro because of the comma differences
! and type declaration in 0d vs nd data
! (fpp doesn't seem to like to parse "," correctly as an argument even
!  if you "#define COMMA ,". It was only able to do the commas in the middle
!  of dimension(...) because the , are in a (). Therefore to have two
!  different types of functions, I had to separate them into two macros
!  with different arguments, i.e., I wasn't able to pass one argument
!  to the macro to expand the type definition, because fpp wasn't
!  able to figure out that it was one argument due to ","s. Also,
!  I wasn't able to pass a macro function with an argument for the same
!  reason, as the preprocessor would expand it and get confused by
!  the commas. Quite frequently, I would get an empty argument.)
#define OPERATE_MULTI_AVG(S, A, B) \
FIRST_HALF(S) A, B, pointer :: in_array, out_array ; \
SECOND_HALF AVG_MAC END_CAP(S)
#define OPERATE_SCALAR_AVG(S, A) \
FIRST_HALF(S) A, pointer :: in_array, out_array ; \
SECOND_HALF AVG_MAC END_CAP(S)
#define OPERATE_MULTI_MIN(S, A, B) \
FIRST_HALF(S) A, B, pointer :: in_array, out_array ; \
SECOND_HALF MIN_MAC END_CAP(S)
#define OPERATE_SCALAR_MIN(S, A) \
FIRST_HALF(S) A, pointer :: in_array, out_array ; \
SECOND_HALF MIN_MAC END_CAP(S)
#define OPERATE_MULTI_MAX(S, A, B) \
FIRST_HALF(S) A, B, pointer :: in_array, out_array ; \
SECOND_HALF MAX_MAC END_CAP(S)
#define OPERATE_SCALAR_MAX(S, A) \
FIRST_HALF(S) A, pointer :: in_array, out_array ; \
SECOND_HALF MAX_MAC END_CAP(S)

! here are all the instantiations
OPERATE_SCALAR_AVG(0r_avg, real(kind=RKIND))
OPERATE_MULTI_AVG(1r_avg, real(kind=RKIND), dimension(:))
OPERATE_MULTI_AVG(2r_avg, real(kind=RKIND), dimension(:, :))
OPERATE_MULTI_AVG(3r_avg, real(kind=RKIND), dimension(:, :, :))
OPERATE_MULTI_AVG(4r_avg, real(kind=RKIND), dimension(:, :, :, :))
OPERATE_MULTI_AVG(5r_avg, real(kind=RKIND), dimension(:, :, :, :, :))
OPERATE_SCALAR_AVG(0i_avg, integer)
OPERATE_MULTI_AVG(1i_avg, integer, dimension(:))
OPERATE_MULTI_AVG(2i_avg, integer, dimension(:, :))
OPERATE_MULTI_AVG(3i_avg, integer, dimension(:, :, :))

OPERATE_SCALAR_MIN(0r_min, real(kind=RKIND))
OPERATE_MULTI_MIN(1r_min, real(kind=RKIND), dimension(:))
OPERATE_MULTI_MIN(2r_min, real(kind=RKIND), dimension(:, :))
OPERATE_MULTI_MIN(3r_min, real(kind=RKIND), dimension(:, :, :))
OPERATE_MULTI_MIN(4r_min, real(kind=RKIND), dimension(:, :, :, :))
OPERATE_MULTI_MIN(5r_min, real(kind=RKIND), dimension(:, :, :, :, :))
OPERATE_SCALAR_MIN(0i_min, integer)
OPERATE_MULTI_MIN(1i_min, integer, dimension(:))
OPERATE_MULTI_MIN(2i_min, integer, dimension(:, :))
OPERATE_MULTI_MIN(3i_min, integer, dimension(:, :, :))

OPERATE_SCALAR_MAX(0r_max, real(kind=RKIND))
OPERATE_MULTI_MAX(1r_max, real(kind=RKIND), dimension(:))
OPERATE_MULTI_MAX(2r_max, real(kind=RKIND), dimension(:, :))
OPERATE_MULTI_MAX(3r_max, real(kind=RKIND), dimension(:, :, :))
OPERATE_MULTI_MAX(4r_max, real(kind=RKIND), dimension(:, :, :, :))
OPERATE_MULTI_MAX(5r_max, real(kind=RKIND), dimension(:, :, :, :, :))
OPERATE_SCALAR_MAX(0i_max, integer)
OPERATE_MULTI_MAX(1i_max, integer, dimension(:))
OPERATE_MULTI_MAX(2i_max, integer, dimension(:, :))
OPERATE_MULTI_MAX(3i_max, integer, dimension(:, :, :))

#undef FIRST_HALF
#undef SECOND_HALF
#undef AVG
#undef MIN
#undef MAX
#undef END_CAP
#undef OPERATE_MULTI_AVG
#undef OPERATE_SCALAR_AVG
#undef OPERATE_MULTI_MIN
#undef OPERATE_SCALAR_MIN
#undef OPERATE_MULTI_MAX
#undef OPERATE_SCALAR_MAX



!***********************************************************************
!  routine typed_operate
!
!> \brief   Do the averaging, but switch on run-time type
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
!  routine ocn_compute_time_averages
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 2, 2015
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
   subroutine ocn_compute_time_averages(domain, timeLevel, err)!{{{
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

      ! do all of the time checking and flag setting
      call timer_checking(domain, err)

      ! update number of accumulations, once only
      do b = 1, size(buffers)
         if (buffers(b) % reset_flag .and. &
             (.not. buffers(b) % delay_reset_flag)) then
            buffers(b) % total_accum = 1
         else if (buffers(b) % accumulate_flag) then
            buffers(b) % total_accum = buffers(b) % total_accum + 1
         end if
      end do

      ! do all of the operations
      do v = 1, size(variables)
         call typed_operate(domain % blocklist, variables(v), operation)
      end do

      ! clear resets
      do b = 1, size(buffers)
         if (buffers(b) % delay_reset_flag) then
            buffers(b) % delay_reset_flag = .false.
         else
            buffers(b) % reset_flag = .false.
         end if
      end do

   end subroutine ocn_compute_time_averages!}}}



!***********************************************************************
!  routine ocn_restart_time_averages
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 2, 2015
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
   subroutine ocn_restart_time_averages(domain, err)!{{{

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
 
   end subroutine ocn_restart_time_averages!}}}



!***********************************************************************
!  routine ocn_finalize_time_averages
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Jon Woodring
!> \date    March 2, 2015
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!-----------------------------------------------------------------------
   subroutine ocn_finalize_time_averages(domain, err)!{{{

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

   end subroutine ocn_finalize_time_averages!}}}



end module ocn_time_averages
! vim: foldmethod=marker
