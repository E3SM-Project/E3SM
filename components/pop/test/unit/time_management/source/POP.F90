!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 program POP_time_managementTest

!----------------------------------------------------------------------
!
!  this program tests the POP time_management.F90 module over
!  a standard range of options.
!
!----------------------------------------------------------------------

   use kinds_mod
   use io
   use io_tools
   use time_management

   implicit none

!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------
   integer (int_kind) :: &
      my_stop_now,       &! flag id for stop_now flag
      coupled_ts,        &! flag id for coupling timestep
      thats_enough,      &! total number of elapsed days
      iyear_end,         &
      imonth_end,        &
      iday_end

   integer (i4) :: iostat,ierr,nml_error

   logical IsOpen

   namelist /driver_nml/ iyear_end, imonth_end, iday_end

   !*** forcing_coupled 
   character (char_len) ::  &
      coupled_freq_opt, qsw_distrb_opt

   integer (int_kind) ::   &
      coupled_freq_iopt,   &! coupler frequency option
      coupled_freq          ! frequency of coupling


   namelist /coupled_nml/ coupled_freq_opt, coupled_freq, qsw_distrb_opt

!----------------------------------------------------------------------
!
!  Read driver namelist
!
!----------------------------------------------------------------------
   

   if (my_task == master_task) then

     open (nml_in, file=nml_filename, status='old',iostat=nml_error)
     write(6,*) ' nml_error = ', nml_error

      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif

      do while (nml_error > 0)
         read(nml_in, nml=driver_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)

   endif ! master_task

!----------------------------------------------------------------------
!
!  Read coupled namelist
!
!----------------------------------------------------------------------
   

   if (my_task == master_task) then

     open (nml_in, file=nml_filename, status='old',iostat=nml_error)
     write(6,*) ' nml_error = ', nml_error

      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif

      do while (nml_error > 0)
         read(nml_in, nml=coupled_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)

   endif ! master_task

!----------------------------------------------------------------------
!
!  initialize stop_now time flag
!
!----------------------------------------------------------------------
   call access_time_flag('stop_now',my_stop_now)

!----------------------------------------------------------------------
!
!  initialize time management
!
!----------------------------------------------------------------------
   call register_string('init_ts')
   call register_string('info_debug_ge2')
   call init_time1

!----------------------------------------------------------------------
!
!  initialize coupling info (from forcing_coupled)
!
!----------------------------------------------------------------------
   if (my_task == master_task) then
     select case (coupled_freq_opt)

     case ('nday')
       if (coupled_freq == 1) then
         coupled_freq_iopt = freq_opt_nday
       else
         call exit_POP (sigAbort, 'coupled_freq_opt = nday failure')
       endif

     case ('nhour')
       if (coupled_freq <= 24) then
         coupled_freq_iopt = freq_opt_nhour
       else
         call exit_POP (sigAbort, 'coupled_freq_opt = nhour failure')
       endif

     case ('nsecond')
       if (coupled_freq <= seconds_in_day) then
         coupled_freq_iopt = freq_opt_nsecond
       else
         call exit_POP (sigAbort, 'coupled_freq_opt = nsecond failure')
       endif

     case ('nstep')
       if (coupled_freq <= nsteps_per_day) then
         coupled_freq_iopt = freq_opt_nstep
       else
         call exit_POP (sigAbort, 'coupled_freq_opt = nstep failure')
       endif

     case default
         call exit_POP (sigAbort, 'coupled_freq_opt case default failure')
     end select
   end if !master_task

   call init_time_flag('coupled_ts', coupled_ts, owner='pop_driver_routine',  &
                        freq_opt = coupled_freq_iopt,freq = coupled_freq)
   call init_time2

   call document_time_flags

   call ymd2eday (iyear_end, imonth_end, iday_end, thats_enough)

!-----------------------------------------------------------------------
!
!  test the advancement of the model in time
!
!-----------------------------------------------------------------------

   write(6,*) ' '
   write(6,*) ' Begin time-stepping loop'
   write(6,*) ' ========================'
   write(6,*) ' '

   stepping_loop: do while (.not. check_time_flag(my_stop_now) )

      call time_manager(.true., .true.)

      call report

     !if (elapsed_days_this_run > thats_enough)  exit stepping_loop

   enddo stepping_loop


!----------------------------------------------------------------------

 end program POP_time_managementTest

 subroutine report

   use time_management

     if (eoy ) then
!!!!!if (eoy .and. mod(iyear,200) == 0 ) then
        write(6,*) '================== eoy ', iyear, ' ================================'
     endif

     if (eom ) then
        write(6,*) '------------------ eom --------------------------------'
     endif

     if (eod ) then
        write(6,*) '.................. eod ................................'
     endif

     if (eoy .or. eom .or. eod) write(6,*) ' '
      
    !write(6,1100) iyear, imonth, iday, seconds_this_day, eoy, eom, eod

 1100 format (1x, i4.4,'-',i2.2,'-',i2.2, 2x, f8.2, 3L3)

 end subroutine report

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
