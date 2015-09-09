! This module contains data and subroutines used by the glc_time_management_test program

module glc_time_management_test_mod
   use glc_constants, only : nml_in
   use glc_kinds_mod
   use glc_communicate, only : my_task, master_task
   use glc_time_management
   use glimmer_paramets, only : stdout
   use writevar_mod

   implicit none
   save

   private

   ! --- parameters ---
   character(len=*), parameter :: test_nml_filename='time_management_test_in'
   integer, parameter :: test_nml_in = 10  ! unit for namelist reads

   ! prefix to write on any line written by the test program
   character(len=*), parameter :: write_prefix = '(time_management_test) '

   ! --- variables read from the time management test namelist ---
   character(len=16) :: runtype_test ! 'initial', 'continue' or 'branch'
   integer :: elapsed_days_test      ! for runtype_test=='continue', the elapsed days at
                                     ! which this run starts
   integer :: nsteps_test            ! number of steps to run for
   real(r8) :: dtt_test              ! if > 0, then use this value for dtt rather than
                                     ! using dt_option & dt_count

   ! --- other module variables ---
   integer :: climate_tstep     ! in the actual code, this is a member of a derived type

   ! --- public variables ---
   public :: write_prefix, nsteps_test

   ! --- public subroutines ---
   public :: read_time_management_test_namelist, &
             init_time_manager, &
             report_time_init, &
             report_time

contains

!***********************************************************************   
   subroutine read_time_management_test_namelist
      ! This subroutine reads the namelist that controls the test program

      namelist /time_management_test_nml/ &
           runtype_test, elapsed_days_test, nsteps_test, dtt_test

      open(test_nml_in, file=test_nml_filename, status='old')
      read(test_nml_in, nml=time_management_test_nml)
      close(test_nml_in)

   end subroutine read_time_management_test_namelist



!***********************************************************************   
   subroutine init_time_manager
      ! Do initialization needed for the time manager.
      ! This includes the relevant code from glc_initMod: glc_initialize

      nml_in = test_nml_in
      runtype = runtype_test

      call init_time1

      ! Allow user to specify their own dtt value, so they're not constrained by the
      ! standard way to specify dtt. If dtt_test <= 0, we use the value of dtt determined
      ! from dt_option / dt_count.
      ! Note that this relies on dtt not having any effect on any other variables in
      ! init_time1
      if (dtt_test > 0) then
         dtt = dtt_test
         dtt_input = dtt_test

         ! steps_per_day and steps_per_year are also set in init_time1 depending on dtt,
         ! but these variables are not important to us for this test driver, so we're not
         ! resetting them here (and note that, for arbitrary dtt, these variables might
         ! be hard to specify)
      end if

      if (my_task==master_task) then
         write(stdout,'(a,a,X,f0.6)') write_prefix, 'dtt =', dtt
      end if

      climate_tstep = nint(dtt/3600._r8)   ! convert from sec to integer hours

      call update_for_restart

      call init_time2

   end subroutine init_time_manager


!***********************************************************************   
   subroutine update_for_restart
      ! If we are testing a restart run, then update time management variables
      ! appropriately, as is done for runtype=='continue' in glc_InitMod: glc_initialize

      integer :: nhour_glint

      nhour_glint = 0     ! number of hours glint has run since start of complete simulation
                          ! must be set to correct value if reading from a restart file

      if (runtype_test == 'continue') then
         nhour_glint = elapsed_days_test * 24
         call ymd2eday (iyear0, imonth0, iday0, elapsed_days0)
         elapsed_days = elapsed_days0 + nhour_glint/24     
         call eday2ymd(elapsed_days, iyear, imonth, iday)
         ihour = 0
         iminute = 0
         isecond = 0
         nsteps_total = nhour_glint / climate_tstep     


         if (my_task == master_task) then
            write(stdout,'(a,a,X,i0)') write_prefix, 'Successfully read restart, nhour_glint =', nhour_glint
            write(stdout,'(a,a,X,i0,X,i0,X,i0,X,i0)') write_prefix, 'Initial eday/y/m/d:', elapsed_days0, iyear0, imonth0, iday0
            write(stdout,'(a,a,X,i0,X,i0,X,i0,X,i0)') write_prefix, 'eday/y/m/d after restart:', elapsed_days, iyear, imonth, iday
            write(stdout,'(a,a,X,i0)') write_prefix, 'nsteps_total =', nsteps_total
            write(stdout,'(a,a)') write_prefix, 'Initialize glint:'
         endif
      endif

      if (my_task==master_task) then
         write(stdout,'(a,a,X,i0)') write_prefix, 'Initialize glint, nhour_glint =', nhour_glint
      endif
   end subroutine update_for_restart


!***********************************************************************   
   subroutine report_time_init
      ! Write out extra reporting information that is done after initialization
      
      write(stdout,'(a,a)') write_prefix, '----- BEGIN REPORT_TIME_INIT -----'
      call writevar(stop_option, "stop_option", write_prefix, stdout)
      call writevar(runid, "runid", write_prefix, stdout)
      call writevar(runtype, "runtype", write_prefix, stdout)
      call writevar(dt_option, "dt_option", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(stop_count, "stop_count", write_prefix, stdout)
      call writevar(stop_iopt, "stop_iopt", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(end_run_at_midnight, "end_run_at_midnight", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(steps_per_year, "steps_per_year", write_prefix, stdout)
      call writevar(steps_per_day, "steps_per_day", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(iyear0, "iyear0", write_prefix, stdout)
      call writevar(imonth0, "imonth0", write_prefix, stdout)
      call writevar(iday0, "iday0", write_prefix, stdout)
      call writevar(ihour0, "ihour0", write_prefix, stdout)
      call writevar(iminute0, "iminute0", write_prefix, stdout)
      call writevar(isecond0, "isecond0", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(iyear_start_run, "iyear_start_run", write_prefix, stdout)
      call writevar(imonth_start_run, "imonth_start_run", write_prefix, stdout)
      call writevar(iday_start_run, "iday_start_run", write_prefix, stdout)
      call writevar(ihour_start_run, "ihour_start_run", write_prefix, stdout)
      call writevar(iminute_start_run, "iminute_start_run", write_prefix, stdout)
      call writevar(isecond_start_run, "isecond_start_run", write_prefix, stdout)
      call writevar(iday_of_year_start_run, "iday_of_year_start_run", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(iyear_end_run, "iyear_end_run", write_prefix, stdout)
      call writevar(imonth_end_run, "imonth_end_run", write_prefix, stdout)
      call writevar(iday_end_run, "iday_end_run", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(elapsed_days_end_run, "elapsed_days_end_run", write_prefix, stdout)
      call writevar(elapsed_days_max, "elapsed_days_max", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(allow_leapyear, "allow_leapyear", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(dtt, "dtt", write_prefix, stdout)
      call writevar(dtt_input, "dtt_input", write_prefix, stdout)
      write(stdout,'(a,a)') write_prefix, '----- END REPORT_TIME_INIT -----'

   end subroutine report_time_init



!***********************************************************************   
   subroutine report_time
      ! Write out current time information
      
      write(stdout,'(a,a)') write_prefix, '----- BEGIN REPORT_TIME -----'
      call writevar(nsteps_total, "nsteps_total", write_prefix, stdout)
      call writevar(nsteps_run, "nsteps_run", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(eod, "eod", write_prefix, stdout)
      call writevar(eom, "eom", write_prefix, stdout)
      call writevar(eoy, "eoy", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(first_step, "first_step", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(midnight, "midnight", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(adjust_nyears, "adjust_nyears", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(new_dtt_value, "new_dtt_value", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(midnight_next, "midnight_next", write_prefix, stdout)
      call writevar(adjust_nyears_next, "adjust_nyears_next", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(iyear, "iyear", write_prefix, stdout)
      call writevar(imonth, "imonth", write_prefix, stdout)
      call writevar(iday, "iday", write_prefix, stdout)
      call writevar(ihour, "ihour", write_prefix, stdout)
      call writevar(iminute, "iminute", write_prefix, stdout)
      call writevar(isecond, "isecond", write_prefix, stdout)
      call writevar(iday_of_year, "iday_of_year", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(imonth_next, "imonth_next", write_prefix, stdout)
      call writevar(iday_next, "iday_next", write_prefix, stdout)
      call writevar(ihour_next, "ihour_next", write_prefix, stdout)
      call writevar(iminute_next, "iminute_next", write_prefix, stdout)
      call writevar(isecond_next, "isecond_next", write_prefix, stdout)
      call writevar(iday_of_year_next, "iday_of_year_next", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(iyear_last, "iyear_last", write_prefix, stdout)
      call writevar(imonth_last, "imonth_last", write_prefix, stdout)
      call writevar(iday_last, "iday_last", write_prefix, stdout)
      call writevar(ihour_last, "ihour_last", write_prefix, stdout)
      call writevar(iday_of_year_last, "iday_of_year_last", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(days_in_year, "days_in_year", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(elapsed_days, "elapsed_days", write_prefix, stdout)
      call writevar(elapsed_days0, "elapsed_days0", write_prefix, stdout)
      call writevar(elapsed_days_jan1, "elapsed_days_jan1", write_prefix, stdout)
      call writevar(elapsed_days_this_run, "elapsed_days_this_run", write_prefix, stdout)
      call writevar(elapsed_days_this_year, "elapsed_days_this_year", write_prefix, stdout)
      call writevar(elapsed_days_init_date, "elapsed_days_init_date", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(elapsed_months, "elapsed_months", write_prefix, stdout)
      call writevar(elapsed_months_this_run, "elapsed_months_this_run", write_prefix, stdout)
      call writevar(elapsed_months_init_date, "elapsed_months_init_date", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(elapsed_years, "elapsed_years", write_prefix, stdout)
      call writevar(elapsed_years_this_run, "elapsed_years_this_run", write_prefix, stdout)
      call writevar(elapsed_years_init_date, "elapsed_years_init_date", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(seconds_this_year, "seconds_this_year", write_prefix, stdout)
      call writevar(seconds_this_day, "seconds_this_day", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(seconds_this_year_next, "seconds_this_year_next", write_prefix, stdout)
      call writevar(seconds_this_day_next, "seconds_this_day_next", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(seconds_in_year, "seconds_in_year", write_prefix, stdout)
      call writevar(hours_in_year, "hours_in_year", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(frac_day, "frac_day", write_prefix, stdout)
      call writevar(tyear, "tyear", write_prefix, stdout)
      call writevar(tmonth, "tmonth", write_prefix, stdout)
      call writevar(tday, "tday", write_prefix, stdout)
      call writevar(thour, "thour", write_prefix, stdout)
      call writevar(tsecond, "tsecond", write_prefix, stdout)
      call writevar(tsecond_old, "tsecond_old", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(newhour, "newhour", write_prefix, stdout)
      call writevar(leapyear, "leapyear", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(tyear00, "tyear00", write_prefix, stdout)
      call writevar(tsecond00, "tsecond00", write_prefix, stdout)
      call writevar(tday00, "tday00", write_prefix, stdout)
      call writevar(thour00, "thour00", write_prefix, stdout)
      write(stdout,'(a)') write_prefix
      call writevar(stepsize, "stepsize", write_prefix, stdout)
      call writevar(stepsize_next, "stepsize_next", write_prefix, stdout)           
      write(stdout,'(a,a)') write_prefix, '----- END REPORT_TIME -----'

   end subroutine report_time

!***********************************************************************   

end module glc_time_management_test_mod
