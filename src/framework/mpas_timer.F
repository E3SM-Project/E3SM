! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_timer
!
!> \brief   MPAS Timer module
!> \author  Michael Duda, Doug Jacobsen
!> \date    03/27/13
!> \details 
!> This module provides developers with internal timer routines. These can be
!> use to profile various parts of code within MPAS. Calls to TAU happen in this module as well
!> to provide more detailed profiling.
!
!-----------------------------------------------------------------------
      module mpas_timer

        use mpas_kind_types
        use mpas_grid_types
        use mpas_io_units

        implicit none
        save

#ifdef _PAPI
        include 'f90papi.h'
#endif

        type timer_node
          character (len=StrKIND) :: timer_name
          logical :: running, printable
          integer :: levels, calls
          real (kind=RKIND) :: start_time, end_time, total_time
          real (kind=RKIND) :: max_time, min_time, avg_time
          real (kind=RKIND) :: efficiency
          type (timer_node), pointer :: next
        end type timer_node

        type (timer_node), pointer :: all_timers
        integer :: levels, synced

        type (dm_info), pointer :: domain_info

        public :: mpas_timer_start, &
                  mpas_timer_stop, &
                  mpas_timer_write, &
                  mpas_timer_init

        contains

!***********************************************************************
!
!  routine mpas_timer_start
!
!> \brief   MPAS Timer start routine
!> \author  Doug Jacobsen
!> \date    03/27/13
!> \details 
!> This routine starts a timer. By default, timer_name is searched for in the linked list of timers.
!> If timer_ptr is provided, the search doesn't happen and the pointer to the timer is used (or allocated if not created yet).
!
!-----------------------------------------------------------------------
        subroutine mpas_timer_start(timer_name, clear_timer, timer_ptr)!{{{
#         ifdef _MPI
          include 'mpif.h'
#         endif

          character (len=*), intent (in) :: timer_name !< Input: name of timer, stored as name of timer
          logical, optional, intent(in) :: clear_timer !< Input: flag to clear timer
          type (timer_node), optional, pointer :: timer_ptr !< Output: pointer to store timer in module

          logical :: timer_added, timer_found, string_equal, check_flag
          type (timer_node), pointer :: current, temp

          integer :: clock, hz, usecs

#ifdef MPAS_TAU 
          call tau_start(timer_name)
#endif

          timer_added = .false.
          timer_found = .false.

          if(.not.associated(all_timers)) then
            timer_added = .true.
            allocate(all_timers)
            allocate(all_timers%next)
            levels = 0

            all_timers%timer_name = ''
            current => all_timers%next
            nullify(current%next)
          else 
            current => all_timers%next
            timer_search: do while ((.not.timer_found) .and. associated(current))
              string_equal = (trim(current%timer_name) == trim(timer_name))
              if(string_equal) then
                timer_found = .true.
              else
                current => current%next
              endif
            end do timer_search
          endif

          if(present(timer_ptr)) then
            timer_found = .true.
            if(.not.associated(timer_ptr)) then
              current => all_timers
              find_end_ptr: do while((.not.timer_added) .and. (associated(current%next)))
                current => current%next
              end do find_end_ptr

              allocate(timer_ptr)

              current%next => timer_ptr
              current => timer_ptr
              nullify(timer_ptr%next)
              current%levels = levels
              current%timer_name = timer_name
              current%running = .false.
              current%total_time = 0.0
              current%max_time = 0.0
              current%min_time = 100000000.0
              current%avg_time = 0.0
              current%calls = 0
              current%efficiency = 0.0
            else
              current => timer_ptr
            endif
          endif

          if(.not.timer_found) then
            current => all_timers
            find_end: do while((.not.timer_added) .and. (associated(current%next)))
              current => current%next
            end do find_end

            allocate(current%next)
            current => current%next

            nullify(current%next)
            timer_added = .true.
          endif

          if(timer_added .and. (.not.timer_found)) then
            current%levels = levels
            current%timer_name = timer_name
            current%running = .false.
            current%total_time = 0.0
            current%max_time = 0.0
            current%min_time = 100000000.0
            current%avg_time = 0.0
            current%calls = 0
            current%efficiency = 0.0
          endif

          if((timer_added .or. timer_found) .and. (.not.current%running)) then
            current%running = .true.
            levels = levels + 1

#ifdef _PAPI
            call PAPIF_get_real_usec(usecs, check_flag)
            current%start_time = usecs/1.0e6
#elif _MPI
            current%start_time = MPI_Wtime()
#else
            call system_clock (count=clock)
            call system_clock (count_rate=hz)
            current%start_time = real(clock,kind=RKIND)/real(hz,kind=RKIND)
#endif
          endif

          if(present(clear_timer)) then
            if(clear_timer) then
              current%start_time = 0.0
              current%end_time = 0.0
              current%total_time = 0.0
              current%max_time = 0.0
              current%min_time = 0.0
              current%avg_time = 0.0
              current%calls = 0
              current%running = .false.
            endif
          endif

          if(present(timer_ptr)) then
              timer_ptr => current
          endif
          
        end subroutine mpas_timer_start!}}}
       
!***********************************************************************
!
!  routine mpas_timer_stop
!
!> \brief   MPAS Timer stop routine
!> \author  Doug Jacobsen
!> \date    03/27/13
!> \details 
!> This routine stops a timer. By default, timer_name is searched for in the linked list of timers.
!> If timer_ptr is provided, the search doesn't happen and the pointer to the timer is used.
!
!-----------------------------------------------------------------------
        subroutine mpas_timer_stop(timer_name, timer_ptr)!{{{
#         ifdef _MPI
          include 'mpif.h'
#         endif

          character (len=*), intent(in) :: timer_name !< Input: name of timer to stop
          type (timer_node), pointer, optional :: timer_ptr !< Input: pointer to timer, for stopping

          type (timer_node), pointer :: current
          
          real (kind=RKIND) :: time_temp
          logical :: timer_found, string_equal, check_flag
          integer :: clock, hz, usecs

#ifdef MPAS_TAU 
          call tau_stop(timer_name)
#endif
 
          timer_found = .false.
 
          if(present(timer_ptr)) then
            timer_found = .true.
            current => timer_ptr
          endif
 
          if(.not.associated(all_timers)) then
            print *,' timer_stop :: timer_stop called with no timers initialized'
          else if(.not. timer_found) then
            current => all_timers
            timer_find: do while(.not.timer_found .and. associated(current))
              string_equal = (trim(current%timer_name) == trim(timer_name))

              if(string_equal) then
                timer_found = .true.
              else
                current => current%next
              endif
            end do timer_find
          endif

          if(.not.timer_found) then
            print *,' timer_stop :: timer_stop called with timer_name =', timer_name,' when timer has not been started.'
            stop
          endif

          if(current%running) then
            current%running = .false.
            levels = levels - 1

#ifdef _PAPI
            call PAPIF_get_real_usec(usecs, check_flag)
            current%end_time = usecs/1.0e6
#elif _MPI
            current%end_time = MPI_Wtime()
#else
            call system_clock(count=clock)
            call system_clock(count_rate=hz)
            current%end_time = real(clock,kind=RKIND)/real(hz,kind=RKIND)
#endif
            
            time_temp = current%end_time - current%start_time
            current%total_time = current%total_time + time_temp

            if(time_temp > current%max_time) then
              current%max_time = time_temp
            endif

            if(time_temp < current%min_time) then
              current%min_time = time_temp
            endif

            current%avg_time = current%avg_time + time_temp
            current%calls = current%calls + 1
          endif

        end subroutine mpas_timer_stop!}}}

!***********************************************************************
!
!  recursive routine mpas_timer_write
!
!> \brief   MPAS Timer write routine
!> \author  Doug Jacobsen
!> \date    03/27/13
!> \details 
!> This routine writes all timer output to stdout. It recursively calls 
!> itself until all timers have been written out. Prior to writing timers,
!> this routine calls mpas_timer_sync.
!
!-----------------------------------------------------------------------
        recursive subroutine mpas_timer_write(timer_ptr, total_ptr)!{{{
          type (timer_node), pointer, optional :: timer_ptr !< Input - Optional: Pointer to a specific timer to write out.
          type (timer_node), pointer, optional :: total_ptr !< Input - Optional: Pointer to the total_time timer.
          character (len=StrKIND) :: tname

          logical :: total_found, string_equals
          type (timer_node), pointer :: current, total
          real (kind=RKIND) :: percent
          integer :: i

          total_found = .false.

#ifdef MPAS_SYNC_TIMERS
          if(associated(domain_info) .and. synced == 0) then
            call mpas_timer_sync()
          endif
#endif

          if(present(timer_ptr) .and. (.not.present(total_ptr))) then
            print *,'timer_write :: timer_ptr valid, but total_ptr is not assigned.'
            stop
          else if(present(timer_ptr)) then
            tname = ''
            do i=0,timer_ptr%levels+2
              tname = tname//' '
            end do

            if(timer_ptr%total_time == 0.0d0) then
              timer_ptr%min_time = 0.0d0
              timer_ptr%max_time = 0.0d0
              timer_ptr%avg_time = 0.0d0
              percent = 0.0d0
            else
              timer_ptr%avg_time = timer_ptr%total_time/timer_ptr%calls
              percent = timer_ptr%total_time/total_ptr%total_time
            endif

            write(stdoutUnit,'(i2, 1x, a35, f15.5, i10, 3f15.5, 2f8.2)') timer_ptr%levels, tname(1:timer_ptr%levels)// &
                      timer_ptr%timer_name, timer_ptr%total_time, timer_ptr%calls, timer_ptr%min_time, timer_ptr%max_time, &
                      timer_ptr%avg_time, percent, timer_ptr%efficiency
            return
          endif

          total => all_timers

          find_total: do while((.not.total_found) .and. associated(total))
            string_equals = (trim(total%timer_name) == trim("total time"))
            if(string_equals) then
              total_found = .true.
            else
              total => total%next
            endif
          end do find_total

          if(.not.total_found) then
            print *,' timer_write :: no timer named "total time" found.'
            stop
          end if

          write(stdoutUnit,'(3x, a10, 24x, a15, a10, a13, a15, a15, a12, a12)') 'timer_name', 'total', 'calls', 'min', 'max', &
                    'avg', 'percent', 'efficiency'
          write(stdoutUnit,'(i2, 1x, a35, f15.5, i10, 3f15.5)') total%levels, total%timer_name, total%total_time, total%calls, &
                     total%min_time, total%max_time, total%avg_time

          current => all_timers

          print_timers: do while(associated(current))
            string_equals = (trim(current%timer_name) == trim("total time"))
            string_equals = string_equals .or. (trim(current%timer_name) == trim(" "))

            if(.not.string_equals) then
              call mpas_timer_write(current, total)
              current => current%next
            else
              current => current%next
            endif
          end do print_timers

        end subroutine mpas_timer_write!}}}

!***********************************************************************
!
!  routine mpas_timer_init
!
!> \brief   MPAS Timer init routine
!> \author  Doug Jacobsen
!> \date    03/27/13
!> \details 
!> This routine initializes the mpas_timer setup. It needs to have access to the dminfo object in order to sync timers.
!
!-----------------------------------------------------------------------
        subroutine mpas_timer_init(domain)!{{{
          type (domain_type), intent(in), optional :: domain !< Input - Optional: Domain structure

          if( present(domain) ) then
              domain_info => domain % dminfo
          endif

          synced = 0

        end subroutine mpas_timer_init!}}}

!***********************************************************************
!
!  routine mpas_timer_sync
!
!> \brief   MPAS Timer sync routine
!> \author  Doug Jacobsen
!> \date    03/27/13
!> \details 
!> This routine synchronizes timers across all processors in order to better represent
!> the entire run domain with the timer output.
!
!-----------------------------------------------------------------------
        subroutine mpas_timer_sync()!{{{
          use mpas_dmpar

          type (timer_node), pointer :: current
          real (kind=RKIND) :: all_total_time, all_max_time, all_min_time, all_ave_time

          current => all_timers

          sync_timers: do while(associated(current))
            all_total_time = 0.0
            all_ave_time = 0.0
            all_max_time = 0.0
            all_min_time = 0.0

            call mpas_dmpar_max_real(domain_info, current % total_time, all_total_time)
            call mpas_dmpar_sum_real(domain_info, current % total_time, all_ave_time)
            
            all_ave_time = all_ave_time / domain_info % nprocs

            current % total_time = all_total_time

#ifdef _MPI
            if ( all_total_time > 0.0d0 ) then
               current % efficiency = all_ave_time / all_total_time
            else
               current % efficiency = 1.0
            end if
#else
            current % efficiency = 1.0
#endif

            if ( current % calls > 0 ) then
               current % avg_time = current % total_time / current % calls
            else
               current % avg_time = current % total_time
            end if

            call mpas_dmpar_max_real(domain_info, current % max_time, all_max_time)
            current % max_time = all_max_time

            call mpas_dmpar_min_real(domain_info, current % min_time, all_min_time)
            current % min_time = all_min_time

            current => current % next
          end do sync_timers

          synced = 1

        end subroutine mpas_timer_sync!}}}

      end module mpas_timer

! vim: foldmethod=marker et ts=2
