!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module ice_timers

!BOP
! !MODULE: timers
!
! !DESCRIPTION:
!  This module contains routine for supporting multiple CPU timers
!  and accumulates time for each individual block and node (task).
!
! !REVISION HISTORY:
!  SVN:$Id: ice_timers.F90 20 2006-09-01 17:09:49Z  $
!
! 2005: Adapted from POP by William Lipscomb
!       Replaced 'stdout' by 'nu_diag'
!
! !USES:

   use ice_kinds_mod
   use ice_constants
   use ice_domain
   use ice_global_reductions
   use ice_exit
   use ice_fileunits, only: nu_diag
   use ice_communicate, only: my_task, master_task

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_ice_timers,     &
             get_ice_timer,       &
             ice_timer_clear,     &
             ice_timer_start,     &
             ice_timer_stop,      &
             ice_timer_print,     &
             ice_timer_print_all, &
             ice_timer_check

!EOP
!BOC

!-----------------------------------------------------------------------
! public timers
!-----------------------------------------------------------------------

   integer (int_kind), public ::      &
      timer_total,            &! total time
      timer_step,             &! time stepping
      timer_dynamics,         &! dynamics
      timer_advect,           &! horizontal advection
      timer_column,           &! column
      timer_thermo,           &! thermodynamics
      timer_sw,               &! radiative transfer
      timer_ridge,            &! ridging
      timer_catconv,          &! category conversions
      timer_couple,           &! coupling
      timer_readwrite,        &! read/write
      timer_diags,            &! diagnostics/history
      timer_hist,             &! diagnostics/history
      timer_bound              ! boundary updates
!      timer_tmp                ! for temporary timings

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      max_timers = 50          ! max number of timers

   type timer_data
      character (char_len) :: &
         name                  ! timer name

      logical (log_kind) ::   &
         in_use,              &! true if timer initialized
         node_started          ! true if any thread has started timer

      integer (int_kind) ::   &
         num_blocks,          &! number of blocks using this timer
         num_nodes,           &! number of nodes  using this timer
         num_starts,          &! number of start requests
         num_stops             ! number of stop requests

      real (dbl_kind) :: &
         node_cycles1,        &! cycle number at start for node timer
         node_cycles2          ! cycle number at stop  for node timer

      real (dbl_kind) ::            &
         node_accum_time       ! accumulated time for node timer

      logical (log_kind), dimension(:), pointer :: &
         block_started         ! true if block timer started

      real (dbl_kind), dimension(:), pointer :: &
         block_cycles1,        &! cycle number at start for block timers
         block_cycles2          ! cycle number at stop  for block timers

      real (dbl_kind), dimension(:), pointer :: &
         block_accum_time       ! accumulated time for block timers

   end type

   type (timer_data), dimension(max_timers) :: &
      all_timers               ! timer data for all timers

   integer (int_kind) ::      & 
      cycles_max               ! max clock cycles allowed by system

   real (dbl_kind) ::               &
      clock_rate               ! clock rate in seconds for each cycle


!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_ice_timers
! !INTERFACE:

 subroutine init_ice_timers

! !DESCRIPTION:
!  This routine initializes machine parameters and timer structures
!  for computing cpu time from F90 intrinsic timer functions.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,                 &! dummy loop counters
      cycles              ! count rate returned by sys_clock

!-----------------------------------------------------------------------
!
!  Call F90 intrinsic system_clock to determine clock rate
!  and maximum cycles for single-processor runs.  If no clock 
!  available, print message.  
!
!-----------------------------------------------------------------------

   cycles = 0

   call system_clock(count_rate=cycles, count_max=cycles_max)

   if (cycles /= 0) then
      clock_rate = c1/real(cycles,kind=dbl_kind)
   else
      clock_rate = c0
      write(nu_diag,'(/,a33,/)') '--- No system clock available ---'
   endif

!-----------------------------------------------------------------------
!
!  initialize timer structures
!
!-----------------------------------------------------------------------

   do n=1,max_timers
      all_timers(n)%name = 'unknown_timer_name'

      all_timers(n)%in_use       = .false.
      all_timers(n)%node_started = .false.

      all_timers(n)%num_blocks   = 0
      all_timers(n)%num_nodes    = 0
      all_timers(n)%num_starts   = 0
      all_timers(n)%num_stops    = 0
      all_timers(n)%node_cycles1 = c0
      all_timers(n)%node_cycles2 = c0

      all_timers(n)%node_accum_time = c0

      nullify(all_timers(n)%block_started)
      nullify(all_timers(n)%block_cycles1)
      nullify(all_timers(n)%block_cycles2)
      nullify(all_timers(n)%block_accum_time)
   end do

   call get_ice_timer(timer_total,    'Total',    nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_step,     'Step',     nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_dynamics, 'Dynamics', nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_advect,   'Advection',nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_column,   'Column',   nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_thermo,   'Thermo',   nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_sw,       'Shortwave',nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_ridge,    'Ridging',  nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_catconv,  'Cat Conv', nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_couple,   'Coupling', nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_readwrite,'ReadWrite',nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_diags,    'Diags    ',nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_hist,     'History  ',nblocks,distrb_info%nprocs)
   call get_ice_timer(timer_bound,    'Bound',    nblocks,distrb_info%nprocs)
!   call get_ice_timer(timer_tmp,      '         ',nblocks,distrb_info%nprocs)

!-----------------------------------------------------------------------
!EOC

   end subroutine init_ice_timers

!***********************************************************************
!BOP
! !IROUTINE: get_ice_timer
! !INTERFACE:

 subroutine get_ice_timer(timer_id, name_choice, num_blocks, num_nodes)

! !DESCRIPTION:
!  This routine initializes a timer with a given name and returns a 
!  timer id.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      name_choice               ! input name for this timer

   integer (int_kind), intent(in) :: &
      num_nodes,               &! number of nodes(tasks) using this timer
      num_blocks                ! number of blocks using this timer
                                ! (can be =1 if timer called outside
                                !  threaded region)

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      timer_id           ! timer number assigned to this timer 

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,                 &! dummy loop index
      srch_error          ! error flag for search

!-----------------------------------------------------------------------
!
!  search for next free timer
!
!-----------------------------------------------------------------------

   srch_error = 1

   srch_loop: do n=1,max_timers
      if (.not. all_timers(n)%in_use) then
         srch_error = 0
         timer_id = n

         all_timers(n)%name       = ' '
         all_timers(n)%name       = name_choice
         all_timers(n)%in_use     = .true.
         all_timers(n)%num_blocks = num_blocks
         all_timers(n)%num_nodes  = num_nodes 

         allocate(all_timers(n)%block_started   (num_blocks), &
                  all_timers(n)%block_cycles1   (num_blocks), &
                  all_timers(n)%block_cycles2   (num_blocks), &
                  all_timers(n)%block_accum_time(num_blocks))

         all_timers(n)%block_started    = .false.
         all_timers(n)%block_cycles1    = c0
         all_timers(n)%block_cycles2    = c0
         all_timers(n)%block_accum_time = c0

         exit srch_loop
      endif
   end do srch_loop

   if (srch_error /= 0) &
      call abort_ice('get_ice_timer: Exceeded maximum number of timers')
                    

!-----------------------------------------------------------------------
!EOC

 end subroutine get_ice_timer

!***********************************************************************
!BOP
! !IROUTINE: ice_timer_clear
! !INTERFACE:

 subroutine ice_timer_clear(timer_id)

! !DESCRIPTION:
!  This routine resets the time for a timer which has already been
!  defined.  NOTE: This routine must be called from outside a threaded
!  region to ensure correct reset of block timers.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      timer_id                ! timer number

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  if the timer has been defined, reset all times to 0
!  otherwise exit with an error
!
!-----------------------------------------------------------------------

   if (all_timers(timer_id)%in_use) then
      all_timers(timer_id)%node_started  = .false.
      all_timers(timer_id)%num_starts    = 0
      all_timers(timer_id)%num_stops     = 0
      all_timers(timer_id)%node_cycles1  = c0
      all_timers(timer_id)%node_cycles2  = c0

      all_timers(timer_id)%node_accum_time = c0

      all_timers(timer_id)%block_started(:)    = .false.
      all_timers(timer_id)%block_cycles1(:)    = c0
      all_timers(timer_id)%block_cycles2(:)    = c0
      all_timers(timer_id)%block_accum_time(:) = c0
   else
      call abort_ice &
                 ('ice_timer_clear: attempt to reset undefined timer')
                    
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_timer_clear

!***********************************************************************
!BOP
! !IROUTINE: ice_timer_start
! !INTERFACE:

 subroutine ice_timer_start(timer_id, block_id)

! !DESCRIPTION:
!  This routine starts a given node timer if it has not already
!  been started by another thread.  If block information is available,
!  the appropriate block timer is also started.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      timer_id                   ! timer number

   integer (int_kind), intent(in), optional :: &
      block_id                 ! optional block id for this block
                               ! this must be the actual local address
                               ! of the block in the distribution
                               ! from which it is called
                               ! (if timer called outside of block
                               ! region, no block info required)

   integer (int_kind) :: &
      cycles                   ! count rate return by sys_clock

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  if timer is defined, start it up
!
!-----------------------------------------------------------------------

   if (all_timers(timer_id)%in_use) then

      !***
      !*** if called from within a block loop, start block timers
      !***

      if (present(block_id)) then

         !*** if block timer already started, stop it first

         if (all_timers(timer_id)%block_started(block_id)) &
            call ice_timer_stop(timer_id, block_id)

         !*** start block timer

         all_timers(timer_id)%block_started(block_id) = .true.
         call system_clock(count=cycles)
         all_timers(timer_id)%block_cycles1(block_id) = real(cycles,kind=dbl_kind)

         !*** start node timer if not already started by
         !*** another thread.  if already started, keep track
         !*** of number of start requests in order to match
         !*** start and stop requests
 
         !$OMP CRITICAL

         if (.not. all_timers(timer_id)%node_started) then
            all_timers(timer_id)%node_started = .true.
            all_timers(timer_id)%num_starts   = 1
            all_timers(timer_id)%num_stops    = 0
            call system_clock(count=cycles)
            all_timers(timer_id)%node_cycles1 = real(cycles,kind=dbl_kind)
         else
            all_timers(timer_id)%num_starts = &
            all_timers(timer_id)%num_starts + 1
         endif

         !$OMP END CRITICAL

      !***
      !*** if called from outside a block loop, start node timer
      !***

      else

         !*** stop timer if already started
         if (all_timers(timer_id)%node_started)  &
                                        call ice_timer_stop(timer_id)

         !*** start node timer

         all_timers(timer_id)%node_started = .true.
         call system_clock(count=cycles)
         all_timers(timer_id)%node_cycles1 = real(cycles,kind=dbl_kind)

      endif
   else
      call abort_ice &
                 ('ice_timer_start: attempt to start undefined timer')
                    
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_timer_start
 
!***********************************************************************
!BOP
! !IROUTINE: ice_timer_stop
! !INTERFACE:

 subroutine ice_timer_stop(timer_id, block_id)

! !DESCRIPTION:
!  This routine stops a given node timer if appropriate.  If block 
!  information is available the appropriate block timer is also stopped.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      timer_id                   ! timer number

   integer (int_kind), intent(in), optional :: &
      block_id                 ! optional block id for this block
                               ! this must be the actual local address
                               ! of the block in the distribution
                               ! from which it is called
                               ! (if timer called outside of block
                               ! region, no block info required)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   real (dbl_kind) :: &
      cycles1, cycles2   ! temps to hold cycle info before correction

   integer (int_kind) :: &
      cycles                   ! count rate returned by sys_clock

!-----------------------------------------------------------------------
!
!  get end cycles
!
!-----------------------------------------------------------------------

   call system_clock(count=cycles)
   cycles2 = real(cycles,kind=dbl_kind)

!-----------------------------------------------------------------------
!
!  if timer is defined, stop it
!
!-----------------------------------------------------------------------

   if (all_timers(timer_id)%in_use) then

      !***
      !*** if called from within a block loop, stop block timer
      !***

      if (present(block_id)) then

         all_timers(timer_id)%block_started(block_id) = .false.

         !*** correct for cycle wraparound and accumulate time

         cycles1 = all_timers(timer_id)%block_cycles1(block_id)
         if (cycles2 >= cycles1) then
            all_timers(timer_id)%block_accum_time(block_id) = &
            all_timers(timer_id)%block_accum_time(block_id) + &
               clock_rate*(cycles2 - cycles1)
         else
            all_timers(timer_id)%block_accum_time(block_id) = &
            all_timers(timer_id)%block_accum_time(block_id) + &
               clock_rate*(cycles_max + cycles2 - cycles1)
         endif

         !*** stop node timer if number of requested stops
         !*** matches the number of starts (to avoid stopping
         !*** a node timer started by multiple threads)
 
         cycles1 = all_timers(timer_id)%node_cycles1

         !$OMP CRITICAL

         all_timers(timer_id)%num_stops = &
         all_timers(timer_id)%num_stops + 1

         if (all_timers(timer_id)%num_starts == &
             all_timers(timer_id)%num_stops) then

            all_timers(timer_id)%node_started = .false.
            if (cycles2 >= cycles1) then
               all_timers(timer_id)%node_accum_time = &
               all_timers(timer_id)%node_accum_time + &
                  clock_rate*(cycles2 - cycles1)
            else
               all_timers(timer_id)%node_accum_time = &
               all_timers(timer_id)%node_accum_time + &
                  clock_rate*(cycles_max + cycles2 - cycles1)
            endif

            all_timers(timer_id)%num_starts   = 0
            all_timers(timer_id)%num_stops    = 0

         endif

         !$OMP END CRITICAL

      !***
      !*** if called from outside a block loop, stop node timer
      !***

      else

         !*** correct for wraparound and accumulate time

         all_timers(timer_id)%node_started = .false.
         cycles1 = all_timers(timer_id)%node_cycles1

         if (cycles2 >= cycles1) then
            all_timers(timer_id)%node_accum_time = &
            all_timers(timer_id)%node_accum_time + &
               clock_rate*(cycles2 - cycles1)
         else
            all_timers(timer_id)%node_accum_time = &
            all_timers(timer_id)%node_accum_time + &
               clock_rate*(cycles_max + cycles2 - cycles1)
         endif

      endif
   else
      call abort_ice &
                 ('ice_timer_stop: attempt to stop undefined timer')
                    
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_timer_stop
 
!***********************************************************************
!BOP
! !IROUTINE: ice_timer_print
! !INTERFACE:

 subroutine ice_timer_print(timer_id,stats)

! !DESCRIPTION:
!  Prints the accumulated time for a given timer and optional
!  statistics for that timer. It is assumed that this routine
!  is called outside of a block loop.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      timer_id                ! timer number

   logical (log_kind), intent(in), optional :: &
      stats                   ! if true, print statistics for node
                              !   and block times for this timer

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,icount           ! dummy loop index and counter

   logical (log_kind) :: &
      lrestart_timer     ! flag to restart timer if timer is running
                         ! when this routine is called

   real (dbl_kind) :: &
      local_time,       &! temp space for holding local timer results
      min_time,         &! minimum accumulated time
      max_time,         &! maximum accumulated time
      mean_time          ! mean    accumulated time

   character (41), parameter :: &
      timer_format = "('Timer ',i3,': ',a9,f11.2,' seconds')"

   character (49), parameter :: &
      stats_fmt1 = "('  Timer stats (node): min = ',f11.2,' seconds')",&
      stats_fmt2 = "('                      max = ',f11.2,' seconds')",&
      stats_fmt3 = "('                      mean= ',f11.2,' seconds')",&
      stats_fmt4 = "('  Timer stats(block): min = ',f11.2,' seconds')"

!-----------------------------------------------------------------------
!
!  if timer has been defined, check to see whether it is currently
!  running.  If it is, stop the timer and print the info.
!
!-----------------------------------------------------------------------

   if (all_timers(timer_id)%in_use) then
      if (all_timers(timer_id)%node_started) then
        call ice_timer_stop(timer_id)
        lrestart_timer = .true.
      else
        lrestart_timer = .false.
      endif

      !*** Find max node time and print that time as default timer
      !*** result

      if (my_task < all_timers(timer_id)%num_nodes) then
         local_time = all_timers(timer_id)%node_accum_time
      else
         local_time = c0
      endif
      max_time = global_maxval(local_time,distrb_info)
      
      if (my_task == master_task) then
        write (nu_diag,timer_format) timer_id, &
              trim(all_timers(timer_id)%name),max_time
      endif

      if (present(stats)) then
      if (stats) then

         !*** compute and print statistics for node timer

         min_time = global_minval(local_time,distrb_info)
         mean_time = global_sum(local_time,distrb_info)/ &
                     real(all_timers(timer_id)%num_nodes,kind=dbl_kind)
         if (my_task == master_task) then
            write (nu_diag,stats_fmt1) min_time
            write (nu_diag,stats_fmt2) max_time
            write (nu_diag,stats_fmt3) mean_time
         endif

         !*** compute and print statistics for block timers
         !*** min block time

         local_time = bignum
         do n=1,all_timers(timer_id)%num_blocks
            local_time = min(local_time, &
                             all_timers(timer_id)%block_accum_time(n))
         end do
         min_time = global_minval(local_time,distrb_info)
         if (min_time == bignum) min_time = c0

         !*** max block time

         local_time = -bignum
         do n=1,all_timers(timer_id)%num_blocks
            local_time = max(local_time, &
                             all_timers(timer_id)%block_accum_time(n))
         end do
         max_time = global_maxval(local_time,distrb_info)
         if (max_time == -bignum) min_time = c0

         !*** mean block time

         local_time = c0
         do n=1,all_timers(timer_id)%num_blocks
            local_time = local_time + &
                         all_timers(timer_id)%block_accum_time(n)
         end do
         icount = global_sum(all_timers(timer_id)%num_blocks, &
                             distrb_info)
         if (icount > 0) mean_time=global_sum(local_time,distrb_info)&
                                   /real(icount,kind=dbl_kind)

         if (my_task == master_task) then
            write (nu_diag,stats_fmt4) min_time
            write (nu_diag,stats_fmt2) max_time
            write (nu_diag,stats_fmt3) mean_time
         endif

      endif
      endif

      if (lrestart_timer) call ice_timer_start(timer_id)
   else
      call abort_ice &
                 ('ice_timer_print: attempt to print undefined timer')
                    
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_timer_print

!***********************************************************************
!BOP
! !IROUTINE: ice_timer_print_all
! !INTERFACE:

 subroutine ice_timer_print_all(stats)

! !DESCRIPTION:
!  Prints the accumulated time for a all timers and optional
!  statistics for that timer. It is assumed that this routine
!  is called outside of a block loop.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (log_kind), intent(in), optional :: &
      stats                   ! if true, print statistics for node
                              !   and block times for this timer

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: n ! dummy loop index

!-----------------------------------------------------------------------
!
!  loop through timers anc call timer_print for each defined timer
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(nu_diag,'(/,a19,/)') 'Timing information:'
   endif

   do n=1,max_timers
      if (all_timers(n)%in_use) then
         if (present(stats)) then
            call ice_timer_print(n,stats)
         else
            call ice_timer_print(n)
         endif
      endif
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_timer_print_all

!***********************************************************************
!BOP
! !IROUTINE: ice_timer_check
! !INTERFACE:

 subroutine ice_timer_check(timer_id,block_id)

! !DESCRIPTION:
!  This routine checks a given timer by stopping and restarting the
!  timer.  This is primarily used to periodically accumulate time in 
!  the timer to prevent timer cycles from wrapping around max_cycles.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      timer_id                   ! timer number

   integer (int_kind), intent(in), optional :: &
      block_id                 ! optional block id for this block
                               ! this must be the actual local address
                               ! of the block in the distribution
                               ! from which it is called
                               ! (if timer called outside of block
                               ! region, no block info required)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  stop and restart the requested timer
!
!-----------------------------------------------------------------------

   if (present(block_id)) then
      call ice_timer_stop (timer_id,block_id)
      call ice_timer_start(timer_id,block_id)
   else
      call ice_timer_stop (timer_id)
      call ice_timer_start(timer_id)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_timer_check

!***********************************************************************

 end module ice_timers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
