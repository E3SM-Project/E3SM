! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
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
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This module provides developers with internal timer routines. 
!> Additionally it provides standard interfaces to additional timer libarries,
!> such as tau or gptl.
!> Timers are stored as trees of lists of timers.
!
!-----------------------------------------------------------------------
      module mpas_timer

        use mpas_kind_types
        use mpas_derived_types
        use mpas_dmpar
        use mpas_threading
        use mpas_log

#ifdef MPAS_PERF_MOD_TIMERS
        use perf_mod
#endif

#ifdef MPAS_GPTL_TIMERS
        use gptl
#endif

        implicit none
        save

#ifdef _PAPI
        include 'f90papi.h'
#endif

        type (mpas_timer_root), pointer :: timer_root => null()

        public :: mpas_timer_set_context, &
                  mpas_timer_start, &
                  mpas_timer_stop, &
                  mpas_timer_write, &
                  mpas_timer_init, &
                  mpas_timer_finalize

        contains

!***********************************************************************
!
!  routine mpas_timer_set_context
!
!> \brief   MPAS Timer set context routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine sets the timer_root for the timer infrastructure to point to a specific domain's timer_root.
!> This allows multiple cores to have MPAS timers, and keep their timing information separate.
!
!-----------------------------------------------------------------------
        subroutine mpas_timer_set_context(domain)!{{{

           type (domain_type), intent(in) :: domain

           timer_root => domain % timer_root

        end subroutine mpas_timer_set_context!}}}

!***********************************************************************
!
!  routine mpas_timer_start
!
!> \brief   MPAS Timer start routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine starts a timer named 'timer_name' as a child of whatever the
!> most recently started timer that is still running is.
!
!-----------------------------------------------------------------------
        subroutine mpas_timer_start(timer_name, clear_timer_in)!{{{

          character (len=*), intent (in) :: timer_name !< Input: name of timer, stored as name of timer
          logical, optional, intent(in) :: clear_timer_in !< Input: flag to clear timer

          logical :: setup_timer, timer_found, clear_timer, check_flag

          type (mpas_timer_node), pointer :: current_timer, prev_timer

          character (len=len_trim(timer_name)) :: trimmed_name

          integer :: usecs, nlen, iErr, threadNum, numThreads

          numThreads = mpas_threading_get_max_threads()
          threadNum = mpas_threading_get_thread_num()

          trimmed_name = trim(timer_name)
          nlen = len(trimmed_name)

          clear_timer = .false.

          if ( present(clear_timer_in) ) then
             clear_timer = clear_timer_in
          end if

#ifdef MPAS_TAU_TIMERS
          call tau_start(trimmed_name)
#endif

#ifdef MPAS_PERF_MOD_TIMERS
          call t_startf(trimmed_name)
#endif

#ifdef MPAS_GPTL_TIMERS
          iErr = gptlstart(trimmed_name)
#endif

#ifdef MPAS_NATIVE_TIMERS
             if ( threadNum == 0 ) then
                setup_timer = .false.
   
                ! If the root_timer in timer_root has not been set, create it and use it for the current timer
                if ( .not. associated(timer_root % root_timer) ) then
                   allocate(timer_root % root_timer)
                   current_timer => timer_root % root_timer
                   setup_timer = .true.
                ! If root_timer is set, current_timer should be as well. So look for this timer under current_timer, or 
                ! create it as a child.
                else
                   ! If current_timer is not set, set it to root_timer and search...
                   if ( .not. associated(timer_root % current_timer) ) then
                      timer_root % current_timer => timer_root % root_timer
                   end if
   
                   ! If current_timer doesn't have any children, make a new child that will be this timer
                   if ( .not. associated(timer_root % current_timer % child) ) then
                      allocate(timer_root % current_timer % child)
                      current_timer => timer_root % current_timer % child
                      current_timer % parent => timer_root % current_timer
                      setup_timer = .true.
                   else ! Search through current_timer's children for this timer, or create it as a new child
                      timer_found = .false.
                      current_timer => timer_root % current_timer % child
                      do while (associated(current_timer) .and. .not. timer_found)
                         if ( current_timer % nlen == nlen) then
                            if ( current_timer % timer_name(1:current_timer % nlen) == trimmed_name(1:nlen) ) then
                               timer_found = .true.
                            end if
                         end if
   
                         if ( .not. timer_found ) then
                            prev_timer => current_timer
                            current_timer => current_timer % next
                         end if
                      end do
   
                      if (.not. timer_found) then
                         allocate(prev_timer % next)
                         current_timer => prev_timer % next
                         current_timer % parent => timer_root % current_timer
                         setup_timer = .true.
                      end if
                   end if
                end if
   
                ! Setup timer if needed
                if ( setup_timer ) then
                   current_timer % timer_name = trimmed_name
                   current_timer % nlen = nlen
                   current_timer % printed = .false.
                   allocate(current_timer % running(numThreads))
                   allocate(current_timer % start_time(numThreads))
                   allocate(current_timer % end_time(numThreads))
                   allocate(current_timer % total_time(numThreads))
                   allocate(current_timer % max_time(numThreads))
                   allocate(current_timer % min_time(numThreads))
                   allocate(current_timer % avg_time(numThreads))
                   current_timer % running(:) = .false.
                   current_timer % start_time(:) = 0.0_R8KIND
                   current_timer % end_time(:) = 0.0_R8KIND
                   current_timer % total_time(:) = 0.0_R8KIND
                   current_timer % max_time(:) = -huge(RKIND)
                   current_timer % min_time(:) = huge(RKIND)
                   current_timer % avg_time(:) = 0.0_RKIND
                   current_timer % calls = 0
                end if

                ! Set current timer to be current timer (so timers started after will be children of this timer instance)
                timer_root % current_timer => current_timer
             end if
             call mpas_threading_barrier()

             current_timer => timer_root % current_timer

             if ( clear_timer) then
                if ( threadNum == 0 ) then
                   current_timer % calls = 0
                end if
                current_timer % running(threadNum + 1) = .false.
                current_timer % start_time(threadNum + 1) = 0.0_R8KIND
                current_timer % end_time(threadNum + 1) = 0.0_R8KIND
                current_timer % total_time(threadNum + 1) = 0.0_R8KIND
                current_timer % max_time(threadNum + 1) = -huge(RKIND)
                current_timer % min_time(threadNum + 1) = huge(RKIND)
                current_timer % avg_time(threadNum + 1) = 0.0_RKIND
             end if

             ! Set start time
             ! TODO after making sure the timer structures work as expected...
             if ( threadNum == 0 ) then
                current_timer % calls = current_timer % calls + 1
             end if
             current_timer % running(threadNum + 1) = .true.
#ifdef _PAPI
             call PAPIF_get_real_usec(usecs, check_flag)
             current_timer % start_time(threadNum + 1) = usecs/1.0e6
#else
             call mpas_dmpar_get_time(current_timer % start_time(threadNum + 1))
#endif
#endif

        end subroutine mpas_timer_start!}}}

!***********************************************************************
!
!  routine mpas_timer_stop
!
!> \brief   MPAS Timer stop routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine stops a timer named 'timer_name'. It should be the most
!> recently started timer, so an error is thrown if it is not.
!
!-----------------------------------------------------------------------
        subroutine mpas_timer_stop(timer_name)!{{{

          character (len=*), intent(in) :: timer_name !< Input: name of timer to stop

          character (len=len_trim(timer_name)) :: trimmed_name !< Trimmed timer name

          real (kind=R8KIND) :: temp_time

          integer :: usecs, nlen, iErr, threadNum

          threadNum = mpas_threading_get_thread_num()

          trimmed_name = trim(timer_name)
          nlen = len(trimmed_name)

#ifdef MPAS_TAU_TIMERS
          call tau_stop(trimmed_name)
#endif

#ifdef MPAS_PERF_MOD_TIMERS
          call t_stopf(trimmed_name)
#endif

#ifdef MPAS_GPTL_TIMERS
          iErr = gptlstop(trimmed_name)
#endif

#ifdef MPAS_NATIVE_TIMERS
          ! Timer to stop should be timer_root % current_timer, since you should only be allowed to stop the most 
          ! recently started timer
          if ( .not. associated(timer_root % current_timer) ) then
             call mpas_log_write('Trying to stop a timer when no timer has been started.', MPAS_LOG_CRIT)
          end if

          if ( .not. timer_root % current_timer % nlen == nlen ) then
             if ( .not. timer_root % current_timer % timer_name(1:timer_root % current_timer % nlen) == trimmed_name(1:nlen) ) then
                call mpas_log_write('Trying to stop timer ' // trim(trimmed_name) // ' when ' // &
                           trim(timer_root % current_timer % timer_name) // ' is the most recently started timer.', MPAS_LOG_CRIT)
             end if
          end if

          ! Set stop time
          ! TODO: Compute timer statistics...
          timer_root % current_timer % running(threadNum + 1) = .false.
#ifdef _PAPI
          call PAPIF_get_real_usec(usecs, check_flag)
          current_timer % end_time(threadNum + 1) = usecs/1.0e6
#else
          call mpas_dmpar_get_time(timer_root % current_timer % end_time(threadNum + 1))
#endif
          temp_time = timer_root % current_timer % end_time(threadNum + 1) &
                    - timer_root % current_timer % start_time(threadNum + 1)
          timer_root % current_timer % total_time(threadNum + 1) = timer_root % current_timer % total_time(threadNum + 1) &
                                                                 + temp_time
          timer_root % current_timer % max_time(threadNum + 1) = max(timer_root % current_timer % max_time(threadNum + 1), &
                                                                 real(temp_time, kind=RKIND))
          timer_root % current_timer % min_time(threadNum + 1) = min(timer_root % current_timer % min_time(threadNum + 1), &
                                                                 real(temp_time, kind=RKIND))

          ! Set current_timer to the parent of current_timer
          call mpas_threading_barrier()
          if ( threadNum == 0 ) then
             timer_root % current_timer => timer_root % current_timer % parent
          end if
          call mpas_threading_barrier()
#endif

        end subroutine mpas_timer_stop!}}}

!***********************************************************************
!
!  routine mpas_timer_write_header
!
!> \brief   MPAS Timer write header routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine writes header information about what will be printed for each timer_root.
!
!-----------------------------------------------------------------------
       subroutine mpas_timer_write_header()!{{{

#ifdef MPAS_NATIVE_TIMERS

          character(len=StrKIND) :: msg

          call mpas_log_write('')
          call mpas_log_write('')
          call mpas_log_write(' Timer information: ')
          call mpas_log_write('    Globals are computed across all threads and processors')
          call mpas_log_write('')
          call mpas_log_write(' Columns:')
          call mpas_log_write('    total time: Global max of accumulated time spent in timer')
          call mpas_log_write('    calls: Total number of times this timer was started / stopped.')
          call mpas_log_write('    min: Global min of time spent in a single start / stop')
          call mpas_log_write('    max: Global max of time spent in a single start / stop')
          call mpas_log_write('    avg: Global max of average time spent in a single start / stop')
          call mpas_log_write('    pct_tot: Percent of the timer at level 1')
          call mpas_log_write('    pct_par: Percent of the parent timer (one level up)')
          call mpas_log_write('    par_eff: Parallel efficiency, global average total time / global max total time')
          call mpas_log_write('')
          call mpas_log_write('')

          write(msg,'(3x, a10, 34x, a15, a12, a11, a15, a15, a13, a10, a12)') 'timer_name', 'total', 'calls', 'min', 'max', &
                    'avg', 'pct_tot', 'pct_par', 'par_eff'
          call mpas_log_write(msg)
#endif

        end subroutine mpas_timer_write_header!}}}

!***********************************************************************
!
!  routine mpas_timer_write
!
!> \brief   MPAS Timer write routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine writes all timer output to stdout for the currently set
!> timer_root. It prints timer information in depth first format.
!
!-----------------------------------------------------------------------
       subroutine mpas_timer_write()!{{{
          type (mpas_timer_node), pointer :: current_timer
          character (len=StrKIND) :: indentation
          logical :: next_timer_found
          integer :: iErr, levels, i, numThreads

          real (kind=RKIND) :: percent_inc, percent_exc, efficiency
          real (kind=RKIND) :: global_time, global_max_total, global_ave_total
          character(len=StrKIND) :: msg

#ifdef MPAS_GPTL_TIMERS
#ifdef MPAS_DEBUG
          iErr = gptlpr(timer_root % dminfo % my_proc_id)
#else
          if ( timer_root % dminfo % my_proc_id == 0 ) then
             iErr = gptlpr(timer_root % dminfo % my_proc_id)
          end if
#endif
#endif

#ifdef MPAS_NATIVE_TIMERS
          numThreads = mpas_threading_get_max_threads()
          current_timer => timer_root % root_timer
          levels = 1

          ! Initialize indentation with all spaces.
          indentation = ''
          do i = 1, StrKIND
             indentation(i:i) = ' '
          end do

          do while ( associated(current_timer) )
             if ( current_timer % running(1) ) then
                call mpas_log_write('Timer ' // trim(current_timer % timer_name) // ' is still running.', MPAS_LOG_ERR)
             end if

             if ( .not. current_timer % printed ) then
                ! Compute average timers for each thread
                do i = 1, numThreads
                    current_timer % avg_time(i) = current_timer % total_time(i) / current_timer % calls
                end do

                ! Synchronize timers across threads
                global_ave_total = sum(current_timer % total_time(:))
                global_ave_total = global_ave_total / numThreads
                current_timer % total_time(1) = maxval(current_timer % total_time(:))
                current_timer % max_time(1) = maxval(current_timer % max_time(:))
                current_timer % min_time(1) = minval(current_timer % min_time(:))
                current_timer % avg_time(1) = sum(current_timer % avg_time(:)) / numThreads

                ! Synchronize timers across procs
                call mpas_dmpar_sum_real(timer_root  % dminfo, global_ave_total, global_time)
                global_ave_total = global_time / timer_root % dminfo % nprocs
                call mpas_dmpar_max_real(timer_root % dminfo, real(current_timer % total_time(1), kind=RKIND), global_time)
                global_max_total = global_time
                current_timer % total_time(1) = global_time
                call mpas_dmpar_max_real(timer_root % dminfo, current_timer % max_time(1), global_time)
                current_timer % max_time(1) = global_time
                call mpas_dmpar_min_real(timer_root % dminfo, current_timer % min_time(1), global_time)
                current_timer % min_time(1) = global_time
                call mpas_dmpar_sum_real(timer_root % dminfo, current_timer % avg_time(1), global_time)
                current_timer % avg_time(1) = global_time / timer_root % dminfo % nprocs

                ! Compute percent of run time
                if ( timer_root % root_timer % total_time(1) /= 0.0_RKIND ) then
                   percent_inc = current_timer % total_time(1) / timer_root % root_timer % total_time(1)
                else
                   percent_inc = 0.0_RKIND
                end if

                if ( associated(current_timer % parent) ) then
                   percent_exc = current_timer % total_time(1) / current_timer % parent % total_time(1)
                else
                   percent_exc = 0.0_RKIND
                end if

                ! Compute efficiency
                if ( global_max_total /= 0.0_RKIND ) then
                   efficiency = global_ave_total / global_max_total
                else
                   efficiency = 1.0_RKIND
                end if

                ! Print current_timer
                write(msg,'(i2, 1x, a45, f15.5, i10, 3f15.5, 1x, f8.2, 3x, f8.2, 3x, f8.2)') levels, indentation(1:levels-1) // &
                          current_timer % timer_name, current_timer % total_time(1), current_timer % calls, &
                          current_timer % min_time(1), current_timer % max_time(1), current_timer % avg_time(1), &
                          percent_inc * 100.0_RKIND, percent_exc * 100.0_RKIND, efficiency
                call mpas_log_write(msg)
                current_timer % printed = .true.
             end if


             ! If current_timer has children, move down a level (and increment levels)
             ! If current_timer has no children, print siblings (without incrementing levels)
             ! If current_timer has no children, move back to the parent timer, and decrement levels.
             next_timer_found = .false.
             if ( associated(current_timer % child) ) then
                if ( .not. current_timer % child % printed ) then
                    current_timer => current_timer % child
                    levels = levels + 1
                    next_timer_found = .true.
                end if
             end if

             if ( .not. next_timer_found ) then
                if ( associated(current_timer % next) ) then
                   current_timer => current_timer % next
                else
                   if ( associated(current_timer % parent) ) then
                      current_timer => current_timer % parent
                      levels = levels - 1
                   else
                      nullify(current_timer)
                   end if
                end if
             end if
          end do
#endif

        end subroutine mpas_timer_write!}}}

!***********************************************************************
!
!  routine mpas_timer_init
!
!> \brief   MPAS Timer init routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine initializes the timers within a given domain.
!> Additionally, it sets the context to the timers within the domain, which can
!> be changed later using the mpas_timer_set_context routine.
!
!-----------------------------------------------------------------------
        subroutine mpas_timer_init(domain)!{{{
          type (domain_type), intent(inout) :: domain !< Input/Output: Domain structure

          integer :: iErr

          allocate(domain % timer_root)
          timer_root => domain % timer_root
          timer_root % dminfo => domain % dminfo

#ifdef MPAS_GPTL_TIMERS
          iErr = gptlsetoption(gptloverhead, 0)
          iErr = gptlsetoption(gptlpercent, 0)
          iErr = gptlsetoption(gptlsync_mpi, 1)

          iErr = gptlinitialize()
#endif

        end subroutine mpas_timer_init!}}}

!***********************************************************************
!
!  routine mpas_timer_finalize
!
!> \brief   MPAS Timer finalize routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine destroys all timers for the given domain.
!> Timers are destroyed in depth first format
!
!-----------------------------------------------------------------------
        subroutine mpas_timer_finalize(domain)!{{{
          type (domain_type), intent(inout) :: domain !< Input/Output: Domain structure

#ifdef MPAS_NATIVE_TIMERS
          ! Destroy all children of domain's root_timer
          call mpas_timer_destroy_children(domain % timer_root % root_timer)

          ! Deallocate all arrays of the root_timer, and the timer itself
          deallocate(domain % timer_root % root_timer % start_time)
          deallocate(domain % timer_root % root_timer % end_time)
          deallocate(domain % timer_root % root_timer % total_time)
          deallocate(domain % timer_root % root_timer % max_time)
          deallocate(domain % timer_root % root_timer % min_time)
          deallocate(domain % timer_root % root_timer % avg_time)
          deallocate(domain % timer_root % root_timer % running)
          deallocate(domain % timer_root % root_timer)
#endif

          deallocate(domain % timer_root)

        end subroutine mpas_timer_finalize!}}}

!***********************************************************************
!
!  recursive routine mpas_timer_destroy_chldren
!
!> \brief   MPAS Timer Destroy Children Routine
!> \author  Doug Jacobsen
!> \date    12/22/2015
!> \details
!> This routine destroys all children timers of parent_timer.
!> It does not destroy parent_timer itself.
!
!-----------------------------------------------------------------------
        recursive subroutine mpas_timer_destroy_children(parent_timer)!{{{
           type (mpas_timer_node), pointer :: parent_timer

           type (mpas_timer_node), pointer :: current_timer

           if (.not. associated(parent_timer % child)) then
              return
           end if

           ! First, destroy all children of any sibling in this generation
           current_timer => parent_timer % child
           do while ( associated(current_timer) )
              if ( associated(current_timer % child) ) then
                 call mpas_timer_destroy_children(current_timer)
              end if
              current_timer => current_timer % next
           end do

           ! This generation should not have any children now, so destroy all siblings
           do while ( associated(parent_timer % child) )
              current_timer => parent_timer % child

              if ( associated(current_timer) ) then
                 ! If there is more than one sibling, delete the first one, and set child to be it's sibling
                 if ( associated(current_timer % next) ) then
                    parent_timer % child => current_timer % next
                 else
                    nullify(parent_timer % child)
                 end if

                 ! Deallocate all arrays, and the timer itself
                 deallocate(current_timer % start_time)
                 deallocate(current_timer % end_time)
                 deallocate(current_timer % total_time)
                 deallocate(current_timer % max_time)
                 deallocate(current_timer % min_time)
                 deallocate(current_timer % avg_time)
                 deallocate(current_timer % running)
                 deallocate(current_timer)
              else
                 nullify(parent_timer % child)
              end if
           end do

        end subroutine mpas_timer_destroy_children!}}}

      end module mpas_timer

! vim: foldmethod=marker et ts=2
