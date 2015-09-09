! This program is a test driver for the glc_time_management module.
!
! It allows you to run the time manager for a given set of options, printing the internal
! state of the time manager at each time step.
!
! You can then look at this output to make sure the time manager is working correctly, or
! you can compare it with output from a previous tag for regression testing.

program glc_time_management_test

   use glc_time_management_test_mod
   use glc_time_management, only : time_manager
   use glimmer_paramets, only : stdout

   implicit none

   integer :: n

   call read_time_management_test_namelist

   call init_time_manager
   call report_time_init
   call report_time

   ! cism isn't currently set up to access time flags like 'stop_now', as far as I can
   ! tell (it relies on the coupler telling it when to stop) -- for example, the
   ! access_time_flag subroutine has been removed. So rather than using all the stopping
   ! functionality of the time manager, I am simply requiring that the user specify the
   ! desired number of time steps.
   do n = 1, nsteps_test
      call time_manager
      call report_time
   end do

   write(stdout,'(a,a)') write_prefix, 'SUCCESSFUL TERMINATION OF GLC_TIME_MANAGEMENT_TEST'

end program glc_time_management_test

   
