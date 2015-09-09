!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glc_time_management

!BOP
! !MODULE: glc_time_management

! !DESCRIPTION:
!  This module contains a large number of routines for calendar, time
!  flags and other functions related to model time.
!
! !REVISION HISTORY:
!  SVN:$Id: time_management.F90 923 2006-05-10 22:25:10Z njn01 $
!  Adapted by William Lipscomb from time_management.F90 in POP.
!  Much of the original POP code deleted because not needed here.

!USES:

   use glc_kinds_mod
   use glc_constants
   use glc_communicate, only: my_task, master_task
   use glc_broadcast,   only: broadcast_scalar
   use glc_exit_mod
   use shr_sys_mod

   implicit none
   public
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_time1,        &
             init_time2,        &
             time_manager,              &
             init_time_flag,            &
             set_time_flag,             &
             set_time_flag_last,        &
             check_time_flag,           &
             check_time_flag_freq_opt,  &
             check_time_flag_freq,      &
             time_to_do,                &
             time_to_start,             &
             time_stamp,                &
             int_to_char,               &
             cesm_date_stamp

! !PUBLIC DATA TYPES:

   type time_flag

      character (char_len) :: &
         name                  ! name for flag

      logical (log_kind) ::   &
         value,               &! logical state of flag
         old_value,           &! last state of flag
         default,             &! default state of flag
         has_default           ! T if default defined, F if no default

      integer (int_kind) ::   &
         freq_opt,            &! frequency units for switching flag on
         freq                  ! freq in above units for switching flag

   end type

! !PUBLIC DATA MEMBERS:

!-----------------------------------------------------------------------
!
!  variables for run control
!
!-----------------------------------------------------------------------

   character (char_len) :: &
      stop_option         ,&! specify how to determine stopping time
      runtype             ,&! type of cesm run (initial, continue, branch or hybrid)
      dt_option             ! method to determine tracer timestep size

   character (char_len_long) :: &
      runid                 ! an identifier for the run

   integer (int_kind) :: &
      stop_count        ,&! num of stop_option intervals before stop
                          !   OR date (yyyymmdd) at which model stops
      stop_iopt         ,&! integer value for stop_option
      nsteps_total      ,&! steps (full&half) since beginning of run sequence
      nsteps_run        ,&! steps taken since beginning of this run
      nsteps_per_day      ! integer number of steps per day

   integer (int_kind), private :: &
      stop_now            ,&! time_flag id for stopping
      coupled_ts            ! time_flag id for a coupled timestep

   logical (log_kind)   :: &! this timestep is:
      eod                 ,&!   at the end of the day
      eom                 ,&!   at the end of the month
      eoy                 ,&!   at the end of the year
      first_step          ,&!   first time step
      midnight              !   at midnight
   integer (int_kind)   :: &
      adjust_nyears         !   number of years that we need to increment

   logical (log_kind)   :: &! the next timestep is:
      midnight_next       ,&!   at midnight
      new_dtt_value       ,&! does restart have a new step size
      end_run_at_midnight   ! does model run end at midnight
   integer (int_kind)   :: &
      adjust_nyears_next    !   number of years that we need to increment
 
   real (r8)      ::       &
      steps_per_year      ,&  ! number of timesteps in one year
      steps_per_day       ,&  ! number of timesteps in one day
      dt_tol              ,&  ! used to determine close enough
      dt_tol_year             ! used to determine if seconds_this_year
                              !  is close enough to seconds_in_year

!-----------------------------------------------------------------------
!
!     quantities related to date
!
!-----------------------------------------------------------------------

   character (1) ::        &
      date_separator        ! character to separate year-month-day

   integer (int_kind) ::   &
      iyear               ,&! year    [0,inf)  for present timestep
      imonth              ,&! month   [1,12]          |
      iday                ,&! day     [1,31]          |
      ihour               ,&! hour    [0,23]          |
      iminute             ,&! minute  [0,59]          |
      isecond             ,&! second  [0,59]          |
      iday_of_year          ! day no. [1,365/6]       V

   integer (int_kind) ::   &
      imonth_next         ,&! month   [1,12]    for next timestep
      iday_next           ,&! day     [1,31]          |
      ihour_next          ,&! hour    [0,23]          |
      iminute_next        ,&! minute  [0,59]          |
      isecond_next        ,&! second  [0,59]          |
      iday_of_year_next     ! day no. [1,365/6]       V

   integer (int_kind) ::   &
      iyear_last          ,&! year    [0,inf)   from previous timestep
      imonth_last         ,&! month   [1,12]          |
      iday_last           ,&! day     [1,31]          |
      ihour_last          ,&! hour    [0,23]          |
      iday_of_year_last     ! day no. [1,365/6]       V

   integer (int_kind) ::   &
      iyear0              ,&! initial start date and time
      imonth0             ,&!   for complete run
      iday0               ,&!
      ihour0              ,&!
      iminute0            ,&!
      isecond0              !

   integer (int_kind) ::      &
      iyear_start_run        ,&! initial start date and time
      imonth_start_run       ,&!   for this run              
      iday_start_run         ,&!
      ihour_start_run        ,&!
      iminute_start_run      ,&!
      isecond_start_run      ,&!
      iday_of_year_start_run   !

   integer (int_kind) ::      &
      iyear_end_run          ,&! final date for this run
      imonth_end_run         ,&!  
      iday_end_run             !

   integer (int_kind)   ::     &! number of:
      days_in_year            ,&! days in present year
      elapsed_days            ,&! full days elapsed since    01-01-0000
      elapsed_days0           ,&! full days elapsed between  01-01-0000 
                                !                   and day0 
      elapsed_days_jan1       ,&! full days elapsed prior to 01-01-iyear
      elapsed_days_this_run   ,&! full days elapsed since beginning of
                                !                   this segment of run
      elapsed_days_this_year  ,&! full days elapsed since beginning of yr
      elapsed_days_init_date  ,&! full days elapsed since initial time
      elapsed_days_end_run    ,&! full days elapsed from 01-01-0000 to end
                                !                   of this run
      elapsed_days_max        ,&! maximum number of full days allowed 
      elapsed_months          ,&! full months elapsed since 01-01-0000
      elapsed_months_this_run ,&! full months elapsed since beginning of
                                !                     this segment of run
      elapsed_months_init_date,&! full months elapsed since initial time
      elapsed_years           ,&! full years  elapsed since 01-01-0000
      elapsed_years_this_run  ,&! full years  elapsed since beginning of
                                !                     this segment of run
      elapsed_years_init_date   ! full years  elapsed since initial time
 
   integer (int_kind), parameter :: &
      days_in_leap_year = 366,      & !   days in a leap year
      days_in_norm_year = 365         !   days in a non-leap year

   integer (int_kind), dimension(12) :: &
      days_in_prior_months,  &! cumulative num days in preceeding months
      days_in_month =        &! number of days in each calendar month
        (/31,28,31,  30,31,30, 31,31,30,   31,30,31/)
        !   J  F  M    A  M  J   J  A  S     O  N  D

   real (r8) ::               &
      seconds_this_year      ,&! seconds elapsed since beginning of year
      seconds_this_day       ,&! seconds elapsed this day    
      seconds_this_day_next  ,&! seconds elapsed this day  at next timestep
      seconds_this_year_next ,&! seconds elapsed this year at next timestep
      seconds_in_year        ,&! seconds in one year -- this varies,
                               !         if leap years are allowed
      hours_in_year            ! hours   in one year
 
   real (r8) ::               &
      frac_day               ,&! fraction of the day elapsed today
      tyear                  ,&! decimal elapsed time in years
      tmonth                 ,&! decimal elapsed time in months
      tday                   ,&! decimal elapsed time in days
      thour                  ,&! decimal elapsed time in hours
      tsecond                ,&! decimal elapsed time in seconds
      tsecond_old              ! tsecond from previous timestep

   logical (log_kind) :: &
      newhour           ,&!
      allow_leapyear    ,&! allow leap years?
      leapyear            ! is this a leapyear?

   character (4) ::      &
      cyear               ! character version of year

   character (2) ::      &
      cmonth            ,&! character version of month
      cday              ,&! character version of day
      chour             ,&! character version of hour
      cminute           ,&! character version of minute
      csecond             ! character version of second

   character (3) ::      &
      cmonth3             ! character month in 3-letter form

   character (3), dimension(12), parameter :: &
      month3_all = (/'jan','feb','mar','apr','may','jun', &
                     'jul','aug','sep','oct','nov','dec'/)

   character (2), dimension(12), parameter :: &
      cmonths  = (/'01','02','03','04','05','06',  &
                   '07','08','09','10','11','12'/)
 
   character (2), dimension(31), parameter :: &
      cdays    = (/'01','02','03','04','05','06','07','08','09','10', &
                   '11','12','13','14','15','16','17','18','19','20', &
                   '21','22','23','24','25','26','27','28','29','30', &
                   '31'/)
 
   real (r8), parameter ::            &
      seconds_in_minute =    60.0_r8, &
      seconds_in_hour   =  3600.0_r8, &
      seconds_in_day    = 86400.0_r8, &
      minutes_in_hour   =    60.0_r8

   !*** for forcing calendar

   real (r8), public ::         &
      tyear00                  ,&!
      tsecond00                ,&!
      tday00                   ,&!
      thour00

!-----------------------------------------------------------------------
!
!  parameters for time frequency and start options
!
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &! integer choices for freq option
      freq_opt_never    = 0,        &
      freq_opt_nyear    = 1,        &
      freq_opt_nmonth   = 2,        &
      freq_opt_nday     = 3,        &
      freq_opt_nhour    = 4,        &
      freq_opt_nsecond  = 5,        &
      freq_opt_nstep    = 6

   integer (int_kind), parameter :: &! integer choices for start options
      start_opt_nstep   = 1,        &
      start_opt_nday    = 2,        &
      start_opt_nyear   = 3,        &
      start_opt_date    = 4

   integer (int_kind), parameter :: &
      next_opt_day      = 1,        &
      next_opt_month    = 2,        &
      next_opt_year     = 3,        &
      stop_opt_never    = 0,        &
      stop_opt_sometime = 1 

!-----------------------------------------------------------------------
!
!  user defined time flags
!
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      max_time_flags=99       ! max number of user-defined flags

   type (time_flag), dimension(max_time_flags) :: &
      time_flags              ! array containing user-defined flags

   integer (int_kind) :: &
      num_time_flags = 0

!-----------------------------------------------------------------------
!
!  time-step related constants and variables
!
!-----------------------------------------------------------------------

   real (r8) ::          &
      dt_count          ,&! input count to determine dtt
      dtt               ,&! tracer timestep (sec)
      dtt_input         ,&! tracer timestep (sec) as specified in namelist
                          !   input; may be different from restart value
      stepsize          ,&! size of present timestep (sec)
      stepsize_next       ! size of next timestep (sec)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  a few private variables used only internally
!
!-----------------------------------------------------------------------

   real (r8), private ::          &
      rhour_next,        &! rhour   for next timestep
      rminute_next,      &! rminute for next timestep
      rsecond_next        ! rsecond for next timestep

   logical (kind=log_kind),private   ::   &
      debug_time_management   = .false. 

   ! WJS (12-21-11): xlf has been generating internal compiler errors
   ! for this module, when building on bluefire; adding this unused
   ! variable resolves them
   ! (these errors were generated using IBM XL Fortran for AIX,
   ! V12.1; Version: 12.01.0000.0008)
   logical (log_kind),private   :: dummy_to_make_xlf_happy

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_time1
! !INTERFACE:

 subroutine init_time1

! !DESCRIPTION:
!  Initializes some time manager variables from namelist inputs
!  and sets time step.  Remaining time manager variables are 
!  initialized after restart files are read.
!
! !USES:
   use glc_files, only : nml_filename
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
      nu,                &! i/o unit number
      nm                  ! month index

!-----------------------------------------------------------------------
!
!  namelist input
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      nml_error           ! namelist i/o error flag

   character (char_len) :: &
      message

   namelist /time_manager_nml/                                   &
               runid,           dt_count,       dt_option,       &
               iyear0,          imonth0,        iday0,           &
               ihour0,          iminute0,       isecond0,        &
               stop_option,     stop_count,     date_separator,  &
               allow_leapyear

!-----------------------------------------------------------------------
!
!  Set logical flags to default values.
!
!-----------------------------------------------------------------------

   stop_now    = init_time_flag('stop_now'   ,default=.false.)
   coupled_ts  = init_time_flag('coupled_ts')

   call reset_switches

!-----------------------------------------------------------------------
!
!  set initial values for namelist inputs
!
!-----------------------------------------------------------------------

   runid            = 'unknown_runid'
   allow_leapyear   = .false.
   stop_option      = 'unknown_stop_option'
   stop_count       = -1
   dt_option        = 'steps_per_day'
   dt_count         =  1

   dt_tol     = 1.0e-6
   dt_tol_year= 100.0*dt_tol

   iyear0     = 0
   imonth0    = 1
   iday0      = 1
   ihour0     = 0
   iminute0   = 0
   isecond0   = 0

   date_separator = ' '

!-----------------------------------------------------------------------
!
!  read options from namelist input file
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=time_manager_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)

   if (nml_error /= 0) then
      call exit_glc(sigAbort,'ERROR reading time_manager_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' Time Management:'
      write(stdout,blank_fmt)
      write(stdout,*) ' time_manager_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,time_manager_nml)
      write(stdout,blank_fmt)
   endif

   call broadcast_scalar (runid           , master_task)
   call broadcast_scalar (iyear0          , master_task)
   call broadcast_scalar (imonth0         , master_task)
   call broadcast_scalar (iday0           , master_task)
   call broadcast_scalar (ihour0          , master_task)
   call broadcast_scalar (iminute0        , master_task)
   call broadcast_scalar (isecond0        , master_task)
   call broadcast_scalar (dt_option       , master_task)
   call broadcast_scalar (dt_count        , master_task)
   call broadcast_scalar (stop_option     , master_task)
   call broadcast_scalar (stop_count      , master_task)
   call broadcast_scalar (allow_leapyear  , master_task)
   call broadcast_scalar (date_separator  , master_task)
 

!-----------------------------------------------------------------------
!
!  error checking
!
!-----------------------------------------------------------------------

   ! If the trimmed runid takes up the entire allowed length, there is a good chance that
   ! the intended runid was longer than the allowed length, so abort
   if (len_trim(runid) >= len(runid)) then
      call exit_glc(sigAbort,'runid exceeds max length: '//runid)
   end if

!-----------------------------------------------------------------------
!
!  determine the value for dtt, based upon model input parameters
!
!-----------------------------------------------------------------------

   select case (dt_option)
 
   case('steps_per_year')
      ! WJS (11-22-11): steps_per_year is currently incompatible with allow_leapyear=
      ! .true., since that would require the time step (dtt) to vary each year; this
      ! functionality doesn't exist
      if (allow_leapyear) then
         call exit_glc(sigAbort,'steps_per_year dt option incompatible with allow_leapyear')
      end if

      steps_per_year = dt_count
      steps_per_day  = steps_per_year/days_in_norm_year
      dtt            = seconds_in_day/steps_per_day

   case('steps_per_day')
      steps_per_day  = dt_count
      steps_per_year = steps_per_day*days_in_norm_year
      dtt            = seconds_in_day/steps_per_day

   case('seconds')
      dtt            = dt_count
      steps_per_day  = seconds_in_day/dtt
      steps_per_year = steps_per_day *days_in_norm_year

   case('hours'  )
      dtt            = dt_count*seconds_in_hour
      steps_per_day  = seconds_in_day/dtt
      steps_per_year = steps_per_day*days_in_norm_year

   case default
      call exit_glc(sigAbort,'unknown dt_option')
   end select

  if (verbose .and. my_task==master_task) then
     write(stdout,*) 'dt_option =', trim(dt_option)
     write(stdout,*) 'dt_count =', dt_count
     write(stdout,*) 'seconds_in_day =', seconds_in_day
     write(stdout,*) 'dtt =', dtt
     call shr_sys_flush(stdout)
  endif

   dtt_input = dtt

!-----------------------------------------------------------------------
!
!  check for incompatibility between dtt and allow_leapyear
!
!-----------------------------------------------------------------------

   ! WJS (11-22-11): allow_leapyear = .true. doesn't work with stepsize greater than 1
   ! year, because reduce_seconds doesn't handle this case
   if (allow_leapyear .and. dtt > (seconds_in_day * days_in_norm_year)) then
      call exit_glc(sigAbort,'dt > 1 year incompatible with allow_leapyear')
   end if

!-----------------------------------------------------------------------
!
!  set initial values; some of these may be overwritten by
!  restart input
!
!-----------------------------------------------------------------------

   iyear        = iyear0
   imonth       = imonth0
   iday         = iday0
   ihour        = ihour0
   iminute      = iminute0
   isecond      = isecond0

   nsteps_total = 0

   seconds_this_day = ihour0  *seconds_in_hour   + &
                      iminute0*seconds_in_minute + &
                      isecond0

!-----------------------------------------------------------------------
!
!  define days_in_prior_months; leap-year adjustments are made in
!  subroutine init_timemanager_2
!
!-----------------------------------------------------------------------

   call prior_days (days_in_prior_months, days_in_month)

!-----------------------------------------------------------------------

   call flushm (stdout)
!EOC

 end subroutine init_time1

!***********************************************************************
!BOP
! !IROUTINE: init_time2
! !INTERFACE:

 subroutine init_time2

! !DESCRIPTION:
!  Completes initialization of time manager quantities now that
!  information from restart files is known
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

   integer (int_kind) ::     &
      nm,                    &! month index
      days_in_month_temp,    &! temp for days in month
      days_in_year_end_run,  &! temp for days in last year of run
      ndays_temp              ! temp for number of days

   logical (log_kind) ::     &
      leapyear_test           ! test for leap year

   character (4) ::          &
      cyear_end_run           ! character version of ending year

   character (2) ::          &
      cmonth_end_run,        &! character version of ending month
      cday_end_run            ! character version of ending day

   character (*), parameter :: &
      date_fmt = "(a17, 2x, a4,'-',a2,'-',a2)"
 
!-----------------------------------------------------------------------
!
!  determine the number of days, months, and years elapsed since
!  01-01-0000  and number of days and months elapsed this year
!
!-----------------------------------------------------------------------

   call ymd2eday (iyear , imonth , iday , elapsed_days     )
   call ymd2eday (iyear , 1      , 1    , elapsed_days_jan1)
   call ymd2eday (iyear0, imonth0, iday0, elapsed_days0    )

   elapsed_days_this_year   = elapsed_days - elapsed_days_jan1

   elapsed_years            = iyear
   elapsed_months           = iyear*12 + imonth - 1
   elapsed_days_this_run    = 0
   elapsed_years_this_run   = 0
   elapsed_months_this_run  = 0
   elapsed_days_init_date   = elapsed_days - elapsed_days0
   elapsed_years_init_date  = iyear - iyear0
   elapsed_months_init_date = elapsed_years_init_date*12 - &
                              imonth0 + imonth

   seconds_this_year        = elapsed_days_this_year*seconds_in_day + &
                              seconds_this_day
   iday_of_year             = elapsed_days_this_year + 1

!-----------------------------------------------------------------------
!
!  determine if this is a leap year; set days_in_year and 
!  days_in_prior_months, regardless of value of allow_leapyear
!
!-----------------------------------------------------------------------

   call leap_adjust

!-----------------------------------------------------------------------
!
!  compare the value of dtt selected for this run via namelist input
!  with that from the previous run.  if different, time_manager
!  will execute as if it were not a restart.
!
!-----------------------------------------------------------------------

!lipscomb - TO DO - This may not be needed
   if (dtt /= dtt_input ) then
      new_dtt_value = .true.
      dtt = dtt_input
   endif

!-----------------------------------------------------------------------
!
!  compute seconds_this_year and seconds_the_day for next timestep
!
!-----------------------------------------------------------------------

   stepsize      = dtt
   stepsize_next = dtt
   seconds_this_year_next =  seconds_this_year + stepsize_next
   seconds_this_day_next  =  seconds_this_day  + stepsize_next
   call reduce_seconds (seconds_this_day_next , &
                        seconds_this_year_next, adjust_nyears_next)

!-----------------------------------------------------------------------
!
!  determine tsecond, the number of seconds elapsed since beginning
!  of complete simulation.
!
!-----------------------------------------------------------------------

   tsecond = elapsed_days_init_date*seconds_in_day + seconds_this_day

!-----------------------------------------------------------------------
!
!  check initial iday, imonth, etc values for reasonableness
!
!-----------------------------------------------------------------------

   if ( .not. valid_ymd_hms () ) then
      call exit_glc(sigAbort,'invalid ymd_hms')
   endif

!-----------------------------------------------------------------------
!
!  Compute decimal time in days, months, years, etc
!  NOTE:  newhour is not set initially
!
!-----------------------------------------------------------------------

   call get_tday

   call int_to_char(4,iyear,cyear)
   cday    = cdays     (iday)
   cmonth  = cmonths   (imonth)
   cmonth3 = month3_all(imonth)

   nsteps_run = 0

!-----------------------------------------------------------------------
!
!  set midnight
!
!-----------------------------------------------------------------------

   if (ihour == 0 .and. iminute == 0 .and. isecond == 0) then
      midnight = .true.
   else
      midnight = .false.
   endif

!-----------------------------------------------------------------------
!
!  save iyear, imonth, etc from the beginning of this run
!
!-----------------------------------------------------------------------
      
   iyear_start_run        = iyear
   imonth_start_run       = imonth
   iday_start_run         = iday
   ihour_start_run        = ihour
   iminute_start_run      = iminute
   isecond_start_run      = isecond

   iday_of_year_start_run = iday_of_year

!-----------------------------------------------------------------------
!
!  will this run end exactly at midnight?
!  (this tests only obvious possibilities)
!
!-----------------------------------------------------------------------

   if ( is_near (mod (seconds_in_day, dtt),c0,dt_tol)  .and. &
        is_near (seconds_this_day,         c0,dt_tol) ) then
      end_run_at_midnight = .true.
   else
      end_run_at_midnight = .false.
   endif

!-----------------------------------------------------------------------
!
!  determine iyear, imonth, and iday for the end of this run
!
!-----------------------------------------------------------------------

   stop_iopt = stop_opt_sometime
 
   select case (stop_option)        
 
   case ('never')   !*** coupler or signal catcher stops glc     

      stop_iopt = stop_opt_never     
      iyear_end_run  = 9999
      imonth_end_run = 1        
      iday_end_run   = 1              
      elapsed_days_max = 1e9  
 
   case ('eoy')     !*** stop at end of stop_count years

      if (end_run_at_midnight) then
         iyear_end_run  = iyear + stop_count
         imonth_end_run = 1        
         iday_end_run   = 1              
      else
         if (imonth  == 12 .and. iday == 31) then
            iyear_end_run = iyear + stop_count
         else
            iyear_end_run = iyear + stop_count - 1
         endif 
         imonth_end_run = 12
         iday_end_run   = 31
      endif

   case ('eom')     !*** stop at end of stop_count months
 
      iyear_end_run  = iyear
      imonth_end_run = imonth + stop_count 
 
      call reduce_months (imonth_end_run, iyear_end_run )
 
      if (end_run_at_midnight) then
         iday_end_run = 1
      else
         iday_end_run = days_in_month(imonth_end_run)
      endif
 
   case ('eod')     !*** stop at end of stop_count days

      if (end_run_at_midnight) then
         iyear_end_run  = iyear
         imonth_end_run = imonth 
         iday_end_run   = iday + stop_count
      else
         iyear_end_run  = iyear
         imonth_end_run = imonth 
         iday_end_run   = iday + stop_count - 1
      endif
 
   case ('nyear', 'nyears') !*** stop after stop_count years
                            !*** need not be end of year

      iyear_end_run  = iyear + stop_count
      imonth_end_run = imonth
      iday_end_run   = iday           
      if (allow_leapyear .and. is_leapyear(iyear_end_run)) then
         days_in_year_end_run = days_in_leap_year
      else
         days_in_year_end_run = days_in_norm_year
      endif
      if (is_near(mod(seconds_in_day*days_in_year_end_run, dtt), &  
                  c0, dt_tol) ) then
         end_run_at_midnight = .true.
      else
         end_run_at_midnight = .false.
      endif
 
      case ('nmonth', 'nmonths')  !*** stop after stop_count months
                                  !*** need not be end of month
      iyear_end_run  = iyear   
      imonth_end_run = imonth + stop_count 
      iday_end_run   = iday           

      call reduce_months (imonth_end_run, iyear_end_run )

   case ('nday', 'ndays')    !*** stop after stop_count days
                                !*** identical to 'eod'

      if (end_run_at_midnight) then
         iyear_end_run  = iyear
         imonth_end_run = imonth 
         iday_end_run   = iday + stop_count
      else
         iyear_end_run  = iyear
         imonth_end_run = imonth 
         iday_end_run   = iday + stop_count - 1
      endif
 
   case ('nstep', 'nsteps')   !*** stop after stop_count steps

      ndays_temp     = stop_count/steps_per_day
      iday_end_run   = iday + ndays_temp 
      iyear_end_run  = iyear
      imonth_end_run = imonth 
 
   case ('date')

      call date2ymd (stop_count, iyear_end_run, &
                     imonth_end_run, iday_end_run)
 
   case default
      call exit_glc(sigAbort,'Invalid stop_option: '/&
                                                     &/stop_option)
   end select
 
!-----------------------------------------------------------------------
!
!  if necessary, adjust iyear_end_run, imonth_end_run, iday_end_run
!
!-----------------------------------------------------------------------

   if (is_leapyear(iyear_end_run)) then
      leapyear_test = .true.
   else
      leapyear_test = .false.
   endif

   if (imonth_end_run == 2 .and. stop_option == 'eom' ) then

      if (end_run_at_midnight) then
         imonth_end_run = 3
         iday_end_run   = 1
      else if (leapyear_test) then
         imonth_end_run = 2
         iday_end_run   = 29
      else 
         imonth_end_run = 2
         iday_end_run   = 28
      endif
 
   else if (imonth_end_run == 2 .and. iday_end_run == 29) then
 
      if (.not. leapyear_test) then
         if (end_run_at_midnight) then
            imonth_end_run = 3
            iday_end_run   = 1
         else
            imonth_end_run = 2
            iday_end_run   = 28
         endif
      endif

   else
 
      if (imonth_end_run == 2 .and. leapyear_test) then
         days_in_month_temp = 29
      else
         days_in_month_temp = days_in_month(imonth_end_run)
      endif
 
      do while (iday_end_run > days_in_month_temp)
 
         iday_end_run   = iday_end_run - days_in_month_temp
         imonth_end_run = imonth_end_run + 1
 
         call reduce_months (imonth_end_run, iyear_end_run )
 
         if (allow_leapyear .and. is_leapyear(iyear_end_run)) then
            leapyear_test = .true.
         else
            leapyear_test = .false.
         endif
 
         if (imonth_end_run == 2 .and. is_leapyear(iyear_end_run)) then
            days_in_month_temp = 29
         else
            days_in_month_temp = days_in_month(imonth_end_run)
         endif
 
      enddo

   endif
 
   call ymd2eday (iyear_end_run, imonth_end_run, iday_end_run, &
                  elapsed_days_end_run)
 
   if (stop_iopt /= stop_opt_never)                     &
      elapsed_days_max = elapsed_days_end_run +         &
                         (dtt+dt_tol)/seconds_in_day

   if (elapsed_days_end_run < elapsed_days ) then
      call int_to_char(4, iyear_end_run, cyear_end_run)
      cmonth_end_run = cmonths(imonth_end_run)
      cday_end_run   = cdays  (iday_end_run  )
      if (my_task == master_task) then
         write(stdout,'(a50)') &
            '  Cannot end at a date earlier than starting date.'
         write(stdout,date_fmt) '  Starting date: ', cyear,cmonth,cday
         if (stop_iopt /= stop_opt_never) then
            write(stdout,date_fmt) '  Ending   date: ',           &
                                   cyear_end_run, cmonth_end_run, &
                                   cday_end_run
         else
            write(stdout,'(a17)') '  No ending date.'
            write(stdout,'(a47)') &
               '  Model relies on external signal for stopping.'
         endif
      endif
      call exit_glc(sigAbort,'invalid end date')
   endif
 
!-----------------------------------------------------------------------
!
!  print various time manager options to log (stdout)
!
!-----------------------------------------------------------------------

   call write_time_manager_options

!-----------------------------------------------------------------------
!EOC

 call flushm (stdout)

 end subroutine init_time2

!***********************************************************************
!BOP
! !IROUTINE: time_manager
! !INTERFACE:

 subroutine time_manager

! !DESCRIPTION:
!  This routine updates various time-related variables to their 
!  end-of-step values.  It is called once at the beginning of each 
!  timestep.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  save previous values of tsecond, isec, imin, ihour, etc
!
!-----------------------------------------------------------------------

   tsecond_old = tsecond

   iyear_last        = iyear
   imonth_last       = imonth
   iday_last         = iday
   iday_of_year_last = iday_of_year
   ihour_last        = ihour

!-----------------------------------------------------------------------
!
!  set logical switches to default values
!
!-----------------------------------------------------------------------

   call reset_switches

!-----------------------------------------------------------------------
!
!  increment the timestep counters
!
!-----------------------------------------------------------------------

   nsteps_run   = nsteps_run   + 1
   nsteps_total = nsteps_total + 1

!-----------------------------------------------------------------------
!
!  activate logical switches which depend on nsteps_run, nsteps_total
!
!-----------------------------------------------------------------------

   call set_switches

!-----------------------------------------------------------------------
!
!  determine size of this timestep
!
!-----------------------------------------------------------------------

      stepsize =    dtt

   if (verbose) then   
      write (stdout,*) 'New nsteps_run =', nsteps_run
      write (stdout,*) 'New nsteps_total =', nsteps_total
      write (stdout,*) 'stepsize =', stepsize
      call shr_sys_flush(stdout)
   endif

!-----------------------------------------------------------------------
!
!  compute seconds_this_year and seconds_the_day for this timestep,
!  and adjust seconds counters if necessary
!
!-----------------------------------------------------------------------

   if (nsteps_total == 1 .or.  new_dtt_value) then
      seconds_this_year  = seconds_this_year  + stepsize
      seconds_this_day   = seconds_this_day   + stepsize
      call reduce_seconds (seconds_this_day      , &
                           seconds_this_year     , adjust_nyears)
   else
      seconds_this_year  = seconds_this_year_next
      seconds_this_day   = seconds_this_day_next
      adjust_nyears = adjust_nyears_next
   endif

!-----------------------------------------------------------------------
!
!  if everything is working correctly, then seconds_this_year and
!  seconds_this_day should be their reduced values, regardless of which
!  path was followed in the above conditional; check for this
!
!-----------------------------------------------------------------------

   if (seconds_this_day >= seconds_in_day) then
      call exit_glc(sigAbort,'too large value of seconds_this_day in time_manager')
   end if
   if (seconds_this_year >= seconds_in_year) then
      call exit_glc(sigAbort,'too large value of seconds_this_year in time_manager')
   end if

!-----------------------------------------------------------------------
!
!  determine the size of the next timestep
!
!-----------------------------------------------------------------------

   stepsize_next =    dtt

!-----------------------------------------------------------------------
!
!  compute seconds_this_year and seconds_the_day for next timestep,
!  and adjust seconds counters if necessary
!
!-----------------------------------------------------------------------

   seconds_this_year_next =  seconds_this_year + stepsize_next
   seconds_this_day_next  =  seconds_this_day  + stepsize_next

   call reduce_seconds (seconds_this_day_next , &
                        seconds_this_year_next, adjust_nyears_next)

!-----------------------------------------------------------------------
!
!  compute present year, month, day, hour, minute, and second
!
!-----------------------------------------------------------------------

   call model_date

!-----------------------------------------------------------------------
!
!  compute decimal days, months, and years
!
!-----------------------------------------------------------------------

   call get_tday
   if (ihour /= ihour_last) newhour = .true.

!-----------------------------------------------------------------------
!
!  set all user-defined time flags
!
!-----------------------------------------------------------------------

   call set_time_flag_all

!-----------------------------------------------------------------------
!
!  stop after this timestep?
!
!-----------------------------------------------------------------------

   if (stop_option == 'nstep' .or. stop_option == 'nsteps') then
      if (nsteps_run == stop_count) call set_time_flag(stop_now,.true.)

   else if (stop_iopt /= stop_opt_never .and. eod) then
 
      if (iyear == iyear_end_run .and. imonth == imonth_end_run   &
                                 .and. iday == iday_end_run) then

         call set_time_flag(stop_now,.true.)
 
         if (stop_option == 'eoy' .and. .not. eoy) then
            call set_time_flag(stop_now,.false.)
         endif
         if (stop_option == 'eom' .and. .not. eom) then
            call set_time_flag(stop_now,.false.)
         endif
 
      else if (elapsed_days > elapsed_days_max ) then

         call set_time_flag(stop_now,.true.)

      endif

      if (stop_option == 'eoy' .and.  eoy  .and. &
          elapsed_years_this_run == stop_count)  &
         call set_time_flag(stop_now,.true.)

      if (stop_option == 'eom' .and.  eom  .and. &
          elapsed_months_this_run == stop_count) &
         call set_time_flag(stop_now,.true.)
 
   endif

   if (verbose .and. check_time_flag(stop_now)) then
      write(stdout,*) 'Time manager, stop_now =', check_time_flag(stop_now)
      write(stdout,*) 'stop_option, stop_count =', trim(stop_option), stop_count
   endif

   new_dtt_value = .false.

!-----------------------------------------------------------------------
!     report gkc model time daily
!-----------------------------------------------------------------------

   if (eod) then
    if (my_task == master_task) then
        write(stdout,1000) iyear, cmonth3, iday, seconds_this_day
        call shr_sys_flush(stdout)
    endif
   endif
1000 format (' (time_manager)', ' glc date ', i4.4, '-', a3, '-', &
                                  i2.2,', ', 1pe12.6, ' sec') 

!-----------------------------------------------------------------------
!EOC

 end subroutine time_manager

!***********************************************************************
!BOP
! !IROUTINE: reset_switches
! !INTERFACE:

 subroutine reset_switches

! !DESCRIPTION:
!  Sets most logical switches to default values.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------

   eod                = .false.  ! not end-of-day
   eom                = .false.  ! not end-of-month
   eoy                = .false.  ! not end-of-year

   newhour            = .false.  ! not new hour

   call reset_time_flag_all

!-----------------------------------------------------------------------
!EOC

 end subroutine reset_switches

!***********************************************************************
!BOP
! !IROUTINE: set_switches
! !INTERFACE:

 subroutine set_switches              

! !DESCRIPTION:
!  Determine if logical switches should be set to non-default values
!  for this timestep.  The switches set in this subroutine must depend 
!  ONLY on nsteps\_run, or nsteps\_total.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  if first step, take Euler step
!
!-----------------------------------------------------------------------

   if (first_step) then 
      first_step = .false.
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine set_switches

!***********************************************************************
!BOP
! !IROUTINE: init_time_flag
! !INTERFACE:

 function init_time_flag(flag_name, default, freq_opt, freq)

! !DESCRIPTION:
!  Creates a user-defined time flag with optional default values
!  and frequency.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      flag_name               ! name for this flag

   logical (log_kind), intent(in), optional :: &
      default                 ! default state for this flag

   integer (int_kind), intent(in), optional :: &
      freq_opt,              &! optional freq option for setting flag
      freq                    ! freq in above units  for setting flag

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      init_time_flag          ! flag id which also is integer index 
                              !    into time flag array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,                 &! dummy loop index
      isearch             ! index of flag during search

!-----------------------------------------------------------------------
!
!  search to determine if flag already exists
!
!-----------------------------------------------------------------------

   isearch = 0
   flag_search: do n=1,num_time_flags
      if (trim(time_flags(n)%name) == flag_name) then
         isearch = n
         exit flag_search
      endif
   end do flag_search

!-----------------------------------------------------------------------
!
!  if no flag defined, define new flag and initialize
!
!-----------------------------------------------------------------------

   if (isearch == 0) then  ! no flag exists - define new flag

      num_time_flags = num_time_flags + 1
      isearch = num_time_flags

      time_flags(isearch)%name = flag_name

      time_flags(isearch)%has_default = .false.
      time_flags(isearch)%default     = .false.
      time_flags(isearch)%freq_opt    = freq_opt_never
      time_flags(isearch)%freq        = 0
      time_flags(isearch)%value       = .false.
      time_flags(isearch)%old_value   = .false.
   endif

!-----------------------------------------------------------------------
!
!  set default if requested
!
!  NOTE: If flag previously defined and optional arguments are 
!        present, this will override any previous definition of 
!        optional arguments. user must make sure calls do not 
!        contain optional arguments or else that the last call to 
!        this routine for a specific flag contains desired values.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  set default value of flag
!
!-----------------------------------------------------------------------

   if (present(default)) then
      time_flags(isearch)%has_default = .true.
      time_flags(isearch)%default     = default
      time_flags(isearch)%value       = default
      time_flags(isearch)%old_value   = default
   endif

!-----------------------------------------------------------------------
!
!  define optional frequency for setting flag
!
!-----------------------------------------------------------------------

   if (present(freq_opt)) then
      time_flags(isearch)%freq_opt = freq_opt
      time_flags(isearch)%freq     = freq
   endif

!-----------------------------------------------------------------------
!
!     print time_flag information for debugging purposes
!
!-----------------------------------------------------------------------
 
      if (debug_time_management .and. my_task == master_task) then
        write(stdout,*) ' initialize time_flag(',isearch,')'
        write(stdout,*)  '   name              = '  &
      , trim(time_flags(isearch)%name)
        write(stdout,*) '   has_default        = '  &
      , time_flags(isearch)%has_default
        write(stdout,*) '   default            = '  &
      , time_flags(isearch)%default
        write(stdout,*) '   freq_opt           = '  &
      , time_flags(isearch)%freq_opt
        write(stdout,*) '   freq               = '  &
      , time_flags(isearch)%freq
        write(stdout,*) '   value              = '  &
      , time_flags(isearch)%value
        write(stdout,*) '   old_value          = '  &
      , time_flags(isearch)%old_value
        write(stdout,*) '   '
      endif

!-----------------------------------------------------------------------
!
!  return array index as integer flag id
!
!-----------------------------------------------------------------------

   init_time_flag = isearch

!-----------------------------------------------------------------------
!EOC

 end function init_time_flag

!***********************************************************************
!BOP
! !IROUTINE: set_time_flag
! !INTERFACE:

 subroutine set_time_flag(flag_id, value)

! !DESCRIPTION:
!  Sets the time flag given by flag\_id to the value.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      flag_id                ! index of flag array identifying flag

   logical (log_kind), intent(in) :: &
      value                  ! value requested for flag

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for proper flag id and then set flag
!
!-----------------------------------------------------------------------

   if (flag_id < 1 .or. flag_id > num_time_flags) &
      call exit_glc(sigAbort,'set_time_flag: invalid flag_id')

   time_flags(flag_id)%value = value

!-----------------------------------------------------------------------
!EOC

 end subroutine set_time_flag

!***********************************************************************
!BOP
! !IROUTINE: set_time_flag_last
! !INTERFACE:

 subroutine set_time_flag_last(flag_id, old_value)

! !DESCRIPTION:
!  Sets the old value of time flag given by flag\_id to old\_value.
!
! !REVISION HISTORY:

! !INPUT VARIABLES:

   integer (int_kind), intent(in) :: &
      flag_id                ! index of flag array identifying flag

   logical (log_kind), intent(in) :: &
      old_value              ! old value requested for flag

!-----------------------------------------------------------------------
!
!  check for proper flag id and then set flag
!
!-----------------------------------------------------------------------

   if (flag_id < 1 .or. flag_id > num_time_flags) &
      call exit_glc(sigAbort,'set_time_flag: invalid flag_id')

   time_flags(flag_id)%old_value = old_value

!-----------------------------------------------------------------------
!EOC

 end subroutine set_time_flag_last

!***********************************************************************
!BOP
! !IROUTINE: reset_time_flag
! !INTERFACE:

 subroutine reset_time_flag(flag_id)

! !DESCRIPTION:
!  Sets the time flag given by flag\_id to default value.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      flag_id                ! index of flag array identifying flag

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for proper flag id and then set flag
!
!-----------------------------------------------------------------------

   if (flag_id < 1 .or. flag_id > num_time_flags) &
      call exit_glc(sigAbort,'reset_time_flag: invalid flag_id')

   if (time_flags(flag_id)%has_default) then
      time_flags(flag_id)%value = time_flags(flag_id)%default
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine reset_time_flag

!***********************************************************************
!BOP
! !IROUTINE: reset_time_flag_all
! !INTERFACE:

 subroutine reset_time_flag_all

! !DESCRIPTION:
!  Sets all time flags to default value (if exists).
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

   integer (int_kind) :: n  ! dummy index

!-----------------------------------------------------------------------
!
!  check all flags for default value and set value if default exists
!
!-----------------------------------------------------------------------

   do n=1,num_time_flags
      if (time_flags(n)%has_default) then
         time_flags(n)%value = time_flags(n)%default
      endif
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine reset_time_flag_all

!***********************************************************************
!BOP
! !IROUTINE: check_time_flag
! !INTERFACE:

 function check_time_flag(flag_id)

! !DESCRIPTION:
!  Returns the current value of time flag given by flag\_id.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      flag_id                ! index of flag array identifying flag

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      check_time_flag        ! current value of time flag

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for proper flag id and then return flag value
!
!-----------------------------------------------------------------------

   if (flag_id < 1 .or. flag_id > num_time_flags) then
      call exit_glc(sigAbort,'check_time_flag: invalid flag_id')
   endif

   check_time_flag = time_flags(flag_id)%value

!-----------------------------------------------------------------------
!EOC

 end function check_time_flag

!***********************************************************************
!BOP
! !IROUTINE: check_time_flag_freq_opt
! !INTERFACE:

 function check_time_flag_freq_opt(flag_id)

! !DESCRIPTION:
!  Returns the current frequency of time flag given by flag\_id.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      flag_id                ! index of flag array identifying flag

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      check_time_flag_freq_opt        ! current freqeuncy option of time flag

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for proper flag id and then return flag value
!
!-----------------------------------------------------------------------

   if (flag_id < 1 .or. flag_id > num_time_flags) &
      call exit_glc(sigAbort,'check_time_flag: invalid flag_id')

   check_time_flag_freq_opt = time_flags(flag_id)%freq_opt

!-----------------------------------------------------------------------
!EOC

 end function check_time_flag_freq_opt

!***********************************************************************
!BOP
! !IROUTINE: check_time_flag_freq
! !INTERFACE:

 function check_time_flag_freq(flag_id)

! !DESCRIPTION:
!  Returns the current frequency of time flag given by flag\_id.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      flag_id                ! index of flag array identifying flag

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      check_time_flag_freq        ! current freqeuncy option of time flag

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for proper flag id and then return flag value
!
!-----------------------------------------------------------------------

   if (flag_id < 1 .or. flag_id > num_time_flags) &
      call exit_glc(sigAbort,'check_time_flag_freq: invalid flag_id')

   check_time_flag_freq = time_flags(flag_id)%freq

!-----------------------------------------------------------------------
!EOC

 end function check_time_flag_freq

!***********************************************************************
!BOP
! !IROUTINE: check_time_flag_last
! !INTERFACE:

 function check_time_flag_last(flag_id)

! !DESCRIPTION:
!  Returns the value of time flag given by flag\_id at previous
!  time step
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      flag_id                ! index of flag array identifying flag

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      check_time_flag_last   ! value of time flag at last timestep

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for proper flag id and then return old flag value
!
!-----------------------------------------------------------------------

   if (flag_id < 1 .or. flag_id > num_time_flags) &
      call exit_glc(sigAbort,'check_time_flag_last: invalid flag_id')

   check_time_flag_last = time_flags(flag_id)%old_value

!-----------------------------------------------------------------------
!EOC

 end function check_time_flag_last

!***********************************************************************
!BOP
! !IROUTINE: set_time_flag_all
! !INTERFACE:

 subroutine set_time_flag_all

! !DESCRIPTION:
!  Sets all time flags based on frequency options.
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

   integer (int_kind) :: n  ! dummy index

!-----------------------------------------------------------------------
!
!  if it is time, set time flag and save value from old time
!
!-----------------------------------------------------------------------

   do n=1,num_time_flags

      if (time_flags(n)%freq_opt /= freq_opt_never) then

         time_flags(n)%old_value = time_flags(n)%value
         time_flags(n)%value = time_to_do(time_flags(n)%freq_opt, & 
                                          time_flags(n)%freq)

      endif

   end do

!-----------------------------------------------------------------------

 end subroutine set_time_flag_all

!***********************************************************************
!BOP
! !IROUTINE: time_to_do
! !INTERFACE:

 function time_to_do (in_freq_opt, in_freq)

! !DESCRIPTION:
!  Determines whether it is time to take a particular action based on 
!  input frequency options.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      in_freq_opt,          &! frequency option for this action
      in_freq                ! frequency in above intervals for action

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      time_to_do         ! true if current timestep matches input
                         ! frequency conditions

!EOP
!BOC
!-----------------------------------------------------------------------

   time_to_do = .false.

   select case (in_freq_opt)

   case (freq_opt_nyear)
      if (eoy .and. mod(elapsed_years_init_date,in_freq) == 0) &
         time_to_do = .true.

   case (freq_opt_nmonth)
      if (eom .and. mod(elapsed_months_init_date,in_freq) == 0) &
         time_to_do = .true.

   case (freq_opt_nday)
      if (eod) then
         if (midnight) then
            if (mod(elapsed_days_init_date  ,in_freq) == 0) &
               time_to_do = .true.
         else
            if (mod(elapsed_days_init_date+1,in_freq) == 0) &
               time_to_do = .true.
         endif
      endif

   case (freq_opt_nhour)
      if (newhour .and. mod(ihour,in_freq) == 0) time_to_do = .true.

   case (freq_opt_nsecond)
      if (mod(isecond,in_freq) == 0) time_to_do = .true.

   case (freq_opt_nstep)
      if (mod(nsteps_total,in_freq) == 0) time_to_do = .true.

   case default
   end select

!-----------------------------------------------------------------------
!EOC

 end function time_to_do

!***********************************************************************
!BOP
! !IROUTINE: time_to_start
! !INTERFACE:

 function time_to_start (in_start_opt, in_start)

! !DESCRIPTION:
!  Determines whether it is time to start a particular function based 
!  on input start options.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      in_start_opt,          &! start option for this action
      in_start                ! start after value

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      time_to_start      ! true if current timestep matches input
                         ! start conditions

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      eday_loc           ! temporary value for elapsed days

!-----------------------------------------------------------------------
!
!  check start conditions - do not start if called from initial
!  (nsteps_run == 0) and the condition matches exactly - the 
!  start will instead be triggered during the time step.  This 
!  avoids looking for restarts that do not yet exist.
!
!-----------------------------------------------------------------------

   time_to_start = .false.

   select case (in_start_opt)

   case (start_opt_nstep)
      if (nsteps_total > in_start) then
         time_to_start = .true.
      else if (nsteps_total == in_start .and. nsteps_run /= 0) then
         time_to_start = .true.
      endif

   case (start_opt_nday)
      if (elapsed_days_init_date > in_start) then
         time_to_start = .true.
      else if (elapsed_days_init_date == in_start .and. & 
               nsteps_run /= 0) then
         time_to_start = .true.
      endif

   case (start_opt_nyear)
      if (elapsed_years_init_date > in_start) then
         time_to_start = .true.
      else if (elapsed_years_init_date == in_start .and. & 
               nsteps_run /= 0) then
         time_to_start = .true.
      else if (elapsed_years_init_date == in_start .and. &
               elapsed_days_this_year > 1) then
         time_to_start = .true.
      endif

   case (start_opt_date)
      call date2eday (in_start, eday_loc)
      if (elapsed_days > eday_loc) then
         time_to_start = .true.
      else if (elapsed_days == eday_loc .and. nsteps_run /= 0) then
         time_to_start = .true.
      endif

   case default
      call exit_glc(sigAbort,'unknown start option in time_to_start')
   end select

!-----------------------------------------------------------------------
!EOC

 end function time_to_start

!***********************************************************************
!BOP
! !IROUTINE: model_date
! !INTERFACE:

 subroutine model_date

! !DESCRIPTION:
!  Determines  iyear, imonth, iday, ihour, iminute, isecond, as
!  well as elapsed days, months, and years
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

   real (r8) ::          &
      rhour,             &! number of hours   elapsed today
      rminute,           &! number of minutes beyond the hour
      rsecond,           &! number of seconds beyond the minute
      seconds_today       ! number of seconds elapsed today
 
   integer (int_kind) :: &
      i,                 &! counter
      days_in_prior_years,&! total days in years whose boundaries we have crossed
                           ! in this timestep 
      day_inc,           &! change in the number of days elapsed between timesteps
      month_inc,         &! change in the number of months elapsed between timesteps
      year_inc            ! change in the number of years elapsed between timesteps


   logical (log_kind) ::       &
      next_is_different_day,   &! true if the next timestep is a different day from the current timestep
      next_is_different_month, &! true if the next timestep is a different month from the current timestep
      next_is_month_plus_1      ! true if the next timestep is in the month immediately
                                ! following the month of the current timestep 

!-----------------------------------------------------------------------
!
!  determine iyear
!
!-----------------------------------------------------------------------
   
   days_in_prior_years = 0

   ! We do the following adjustment in a loop to allow the leap year adjustment to be
   ! done properly for each year that we're incrementing (needed for correct
   ! adjustment of days_in_prior_years)
   do i = 1, adjust_nyears
      days_in_prior_years = days_in_prior_years + days_in_year
      iyear = iyear + 1
      
      ! the following sets days_in_year appropriately, among other things
      if (allow_leapyear) call leap_adjust
   end do

!-----------------------------------------------------------------------
!
!  determine iday_of_year, imonth, iday, ihour, iminute, isecond
!                                        rhour, rminute, rsecond
!
!-----------------------------------------------------------------------

   if (nsteps_run == 1) then
      call ymd_hms( seconds_this_year, seconds_this_day, &
                    iday_of_year,                        &
                    imonth, iday   , iday_of_year_last,  &
                    ihour , iminute, isecond,            &
                    rhour , rminute, rsecond,            &
                    midnight       , adjust_nyears)
   else
      iday_of_year  = iday_of_year_next
      imonth        = imonth_next
      iday          = iday_next
      ihour         = ihour_next
      iminute       = iminute_next
      isecond       = isecond_next
      rhour         = rhour_next
      rminute       = rminute_next
      rsecond       = rsecond_next
      midnight      = midnight_next
   endif

!-----------------------------------------------------------------------
!
!  determine iday_of_year, imonth, iday, etc, for next timestep
!
!-----------------------------------------------------------------------
 
   call ymd_hms(seconds_this_year_next,                    &
                seconds_this_day_next,                     &
                iday_of_year_next,                         &
                imonth_next  , iday_next   , iday_of_year, &
                ihour_next   , iminute_next, isecond_next, &
                rhour_next   , rminute_next, rsecond_next, &
                midnight_next, adjust_nyears_next)

!-----------------------------------------------------------------------
!
!  end of day?
!
!-----------------------------------------------------------------------

   next_is_different_day = (adjust_nyears_next > 0 .or. iday_of_year_next /= iday_of_year)
   if (next_is_different_day) then
      if (.not. midnight_next)                     eod = .true.
      if (stepsize_next + dt_tol > seconds_in_day) eod = .true.
   endif

   if (midnight) eod = .true.

!-----------------------------------------------------------------------
!
!  end of month?
!
!-----------------------------------------------------------------------

   next_is_different_month = (adjust_nyears_next > 0 .or. imonth_next /= imonth)
   if (midnight .and. iday == 1) then 
      eom = .true.
   
   else if (next_is_different_month) then
      next_is_month_plus_1 = (adjust_nyears_next == 0 .and. imonth_next == imonth+1) .or. &
                             (adjust_nyears_next == 1 .and. imonth == 12 .and. &
                              imonth_next == 1)
      if (midnight_next .and. iday_next == 1 .and. next_is_month_plus_1) then
         ! the NEXT timestep is considered to be the end of the month
         eom = .false.
      else
         eom = .true.
      end if
   
   else
      eom = .false.
   end if

!-----------------------------------------------------------------------
!
!  elapsed months (integer) 
!
!-----------------------------------------------------------------------

   month_inc = adjust_nyears*12 + imonth - imonth_last
   elapsed_months           = elapsed_months           + month_inc
   elapsed_months_this_run  = elapsed_months_this_run  + month_inc
   elapsed_months_init_date = elapsed_months_init_date + month_inc
 
!-----------------------------------------------------------------------
!
!  end of year?
!
!-----------------------------------------------------------------------

   if (midnight .and. iday_of_year == 1) then
      eoy = .true.
   else if (adjust_nyears_next == 0) then
      eoy = .false.
   else if (adjust_nyears_next == 1) then
      if (midnight_next .and. iday_of_year_next == 1) then
         ! the NEXT timestep is considered to be the end of the year
         eoy = .false.
      else
         eoy = .true.
      end if
   else if (adjust_nyears_next > 1) then 
      eoy = .true.
   else
      call exit_glc(sigAbort,'Unexpected value for adjust_nyears_next in model_date')
   end if

!-----------------------------------------------------------------------
!
!  adjust elapsed years and elapsed days in the year (integer)
!
!-----------------------------------------------------------------------

   year_inc = adjust_nyears
   elapsed_years           = elapsed_years           + year_inc
   elapsed_years_this_run  = elapsed_years_this_run  + year_inc
   elapsed_years_init_date = elapsed_years_init_date + year_inc

   if (year_inc > 0) then
      call ymd2eday (iyear , 1, 1, elapsed_days_jan1)
      elapsed_days_this_year  = elapsed_days - elapsed_days_jan1
   end if

!-----------------------------------------------------------------------
!
!  character values for iyear, imonth, iday
!
!-----------------------------------------------------------------------
      
   if (iyear  /= iyear_last ) call int_to_char(4, iyear, cyear)
 
   if (imonth /= imonth_last) then
      cmonth  = cmonths   (imonth)
      cmonth3 = month3_all(imonth)
   endif

   if (iday   /= iday_last) cday = cdays(iday)                

!-----------------------------------------------------------------------
!
!  elapsed number of days (integer)
!
!-----------------------------------------------------------------------

   day_inc = iday_of_year - iday_of_year_last + days_in_prior_years

   elapsed_days           = elapsed_days           + day_inc
   elapsed_days_this_run  = elapsed_days_this_run  + day_inc
   elapsed_days_this_year = elapsed_days_this_year + day_inc
   elapsed_days_init_date = elapsed_days_init_date + day_inc
 
!-----------------------------------------------------------------------
!
!  has a valid date been selected?
!
!-----------------------------------------------------------------------

   if (.not. valid_ymd_hms()) then
      call exit_glc(sigAbort,'invalid ymd_hms')
   endif

!-----------------------------------------------------------------------
!EOC
 
 end subroutine model_date
 
!***********************************************************************
!BOP
! !IROUTINE: get_tday
! !INTERFACE:

 subroutine get_tday

! !DESCRIPTION:
!  Computes decimal day, month, year, etc.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
! !OUTPUT PARAMETERS:
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  creating floating point values for elapsed time in various units
!
!-----------------------------------------------------------------------

   frac_day= seconds_this_day/seconds_in_day

   tsecond = elapsed_days_init_date*seconds_in_day + seconds_this_day

   tday    = tsecond/seconds_in_day

   tmonth  = real(elapsed_months_init_date,r8) + &
             (real(iday,r8)-c1+frac_day)/days_in_month(imonth)

   tyear   = elapsed_years_init_date + seconds_this_year/seconds_in_year

   thour   = tday*24.0_r8

!-----------------------------------------------------------------------
!
!  define tday00 and thour00 for use in forcing routines.
!  these are the time in days/hours since 01-01-0000.
!
!-----------------------------------------------------------------------

   tyear00   = elapsed_years + seconds_this_year/seconds_in_year
   tday00    = elapsed_days + frac_day
   thour00   = tday00*24.0_r8
   tsecond00 = tday00*seconds_in_day

!-----------------------------------------------------------------------
!EOC

 end subroutine get_tday

!***********************************************************************
!BOP
! !IROUTINE: ymd_hms
! !INTERFACE:

 subroutine ymd_hms(seconds_this_year_loc , seconds_this_day_loc, &
                    iday_of_year_loc,                             &
                    imonth_loc  , iday_loc   , iday_of_year_compare,&
                    ihour_loc   , iminute_loc, isecond_loc,       &
                    rhour_loc   , rminute_loc, rsecond_loc,       &
                    midnight_loc, adjust_nyears_loc)

! !DESCRIPTION:
!  Computes integer values iday\_of\_year, iyear, imonth, iday, ihour, 
!  iminute, isecond.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer  (int_kind), intent(in) :: &
      iday_of_year_compare,           &! day of year to compare to check day change
      adjust_nyears_loc                ! number of years that we need to increment

! !INPUT/OUTPUT PARAMETERS:

   real (r8), intent(inout) :: &
      seconds_this_year_loc   ,&! number of seconds in year
      seconds_this_day_loc      ! number of seconds in day

! !OUTPUT PARAMETERS:

   integer  (int_kind), intent(out) :: &
      imonth_loc,              &! local value of imonth 
      iday_loc,                &! local value of iday
      ihour_loc,               &! local value of ihour
      iminute_loc,             &! local value of iminute
      isecond_loc,             &! local value of isecond
      iday_of_year_loc          ! local value of iday_of_year

   real (r8), intent(out) ::   &
      rhour_loc,               &! real value for hour
      rminute_loc,             &! real value for minute
      rsecond_loc               ! real value for second

   logical (log_kind), intent(out) :: &
      midnight_loc              ! midnight flag
 
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      itest,             &
      nm,                &
      ntest

   real (r8) ::     &
      test_seconds, &
      rtest,        &
      r_ntest

!-----------------------------------------------------------------------
!
!  determine day number   [1,days_in_year]
!
!-----------------------------------------------------------------------
 
   rtest = seconds_this_year_loc/seconds_in_day
   itest =  int (rtest)
   ntest = nint (rtest)
   r_ntest = ntest

   if (is_near(rtest, r_ntest, dt_tol)) then
      iday_of_year_loc = ntest + 1
   else
      iday_of_year_loc = itest + 1
   endif
 
!-----------------------------------------------------------------------
!
!  determine month number [1,12]                
!
!-----------------------------------------------------------------------
 
   imonth_loc = 12

   do nm = 1,11
      if (iday_of_year_loc >  days_in_prior_months(nm)  .and. &
          iday_of_year_loc <= days_in_prior_months(nm+1))     &
         imonth_loc = nm
   enddo

!-----------------------------------------------------------------------
!
!  determine day-of-month number [1,31]
!
!-----------------------------------------------------------------------
 
   iday_loc = iday_of_year_loc - days_in_prior_months(imonth_loc)
 
!-----------------------------------------------------------------------
!
!  determine integer hour, minute, and second
!
!-----------------------------------------------------------------------
 
   call hms (seconds_this_day_loc,                &
             ihour_loc, iminute_loc, isecond_loc, &
             rhour_loc, rminute_loc, rsecond_loc)
 
!-----------------------------------------------------------------------
!
!  midnight?
!
!-----------------------------------------------------------------------
 
   if (ihour_loc == 0 .and. iminute_loc == 0 .and. &
                            isecond_loc == 0) then
      midnight_loc = .true.
   else
      midnight_loc = .false.
   endif

!-----------------------------------------------------------------------
!
!  check for unhandled conditions
!
!-----------------------------------------------------------------------

   ! WJS (11-16-11): These conditions used to be handled, by adjusting iday_loc,
   ! imonth_loc, adjust_year_loc and iday_of_year_loc appropriately. However, I am
   ! concerned that the interactions between these adjustments (in particular, the
   ! adjustment to adjust_year_loc - which is now adjust_nyears_loc) and my change of
   ! adjust_year_loc (logical) to adjust_nyears_loc (integer) weren't being handled
   ! correctly.  Furthermore, I can't see any situation where these conditions will be
   ! triggered. So rather than trying to handle them properly, I am treating them as error
   ! conditions.

   ! Something similar to this condition used to increment iday_loc
   if (iday_of_year_loc == iday_of_year_compare .and. adjust_nyears_loc == 0 .and. &
        midnight_loc) then
      call exit_glc(sigAbort,'Unhandled condition in ymd_hms: midnight condition')
   end if

   ! This condition used to adjust iday_loc, imonth_loc, and adjust_nyears_loc; but I
   ! think it was only necessary because of the possible increase in iday_loc due to the
   ! above, now-unhandled condition
   if (iday_loc > days_in_month(imonth_loc)) then
      call exit_glc(sigAbort,'Unhandled condition in ymd_hms: iday > days_in_month')
   end if

!-----------------------------------------------------------------------
!EOC

 end subroutine ymd_hms

!***********************************************************************
!BOP
! !IROUTINE: hms
! !INTERFACE:

 subroutine hms (seconds_loc,                           &
                 ihour_loc  , iminute_loc, isecond_loc, &
                 rhour_loc  , rminute_loc, rsecond_loc)

! !DESCRIPTION:
!  Determines present hour, minute, and second.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), intent(inout) :: &
      seconds_loc         ! elapsed seconds in current day 
 
! !OUTPUT PARAMETERS:

   integer(log_kind), intent(out) :: &
      ihour_loc,         &! hour in current day
      iminute_loc,       &! minute in current hour
      isecond_loc         ! seconds in current minute

   real (r8), intent(out) :: &
      rhour_loc,         &! real values for the above quantities
      rminute_loc,       &
      rsecond_loc
 
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  compute present hour, minute, and second
!
!-----------------------------------------------------------------------

   rhour_loc   = seconds_loc/seconds_in_hour
   ihour_loc   = rhour_loc

   rminute_loc = (rhour_loc - ihour_loc)*minutes_in_hour
   iminute_loc =  rminute_loc

   rsecond_loc = (rminute_loc - iminute_loc)*seconds_in_minute
   isecond_loc = nint (rsecond_loc)
 
!-----------------------------------------------------------------------
!
!  corrections to second, minute, and/or hour
!
!-----------------------------------------------------------------------

   if (isecond_loc == 60) then
      isecond_loc   =  0
      iminute_loc   =  iminute_loc+ 1
   endif

   if (iminute_loc == 60) then
      iminute_loc   =  0
      ihour_loc     =  ihour_loc + 1
   endif

   if (ihour_loc   == 24) then
      ihour_loc     =  0
   endif

!-----------------------------------------------------------------------
!
!  if h:m:s == 0:00:00, then adjust seconds 
!
!-----------------------------------------------------------------------

   if (ihour_loc == 0  .and. iminute_loc == 0  .and. & 
                             isecond_loc == 0) seconds_loc = c0
 
!-----------------------------------------------------------------------
!EOC

 end subroutine hms

!***********************************************************************
!BOP
! !IROUTINE: reduce_months
! !INTERFACE:

 subroutine reduce_months (imonth_loc, iyear_loc)

! !DESCRIPTION:
!  Reduces imonth such that it never exceeds 12 and
!  increments iyear accordingly.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), intent (inout) :: &
      imonth_loc,          &! current value of imonth
      iyear_loc             ! current value of iyear

!EOP
!BOC
!-----------------------------------------------------------------------
 
   do while (imonth_loc > 12)
      imonth_loc = imonth_loc - 12
      iyear_loc  = iyear_loc  + 1
   enddo
 
!-----------------------------------------------------------------------
!EOC

 end subroutine reduce_months

!***********************************************************************
!BOP
! !IROUTINE: reduce_seconds
! !INTERFACE:

 subroutine reduce_seconds (seconds_this_day_loc, &
                            seconds_this_year_loc, adjust_nyears_loc)

! !DESCRIPTION:
!  Reduce seconds\_this\_day and seconds\_this year, if either
!  exceeds their bounds (eg due to roundoff).
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   real (r8), intent(inout) :: &
      seconds_this_day_loc,    &! current value of seconds_this_day
      seconds_this_year_loc     ! current value of seconds_this_year
 
! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      adjust_nyears_loc         ! number of years that we need to increment

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      ns, ns_end

!-----------------------------------------------------------------------
!
!  if seconds_this_day exceeds the number of seconds in a day, then
!  reset seconds_this_day
!
!-----------------------------------------------------------------------

   if (seconds_this_day_loc >= seconds_in_day) then

      ns_end = nint(seconds_this_day_loc/seconds_in_day)
      do ns = 1, ns_end
         if (seconds_this_day_loc + dt_tol >= seconds_in_day) &
            seconds_this_day_loc = seconds_this_day_loc - &
                                   seconds_in_day
      enddo
 
   endif

!-----------------------------------------------------------------------
!
!  if seconds_this_year exceeds the number of seconds in a year, then
!  reset seconds_this_year
!
!-----------------------------------------------------------------------

   adjust_nyears_loc = 0

   do while (seconds_this_year_loc >= seconds_in_year - stepsize .and.     &
       (seconds_this_year_loc >= seconds_in_year .or.                &
        is_near(seconds_this_year_loc,seconds_in_year,dt_tol_year)))
 
      seconds_this_year_loc = seconds_this_year_loc - seconds_in_year

      if (is_near(seconds_this_year_loc, c0, dt_tol)) then
         seconds_this_year_loc = c0
         seconds_this_day_loc  = c0
      endif
 
      adjust_nyears_loc = adjust_nyears_loc + 1
   end do

!-----------------------------------------------------------------------
!EOC
 
 end subroutine reduce_seconds

!***********************************************************************
!BOP
! !IROUTINE: leap_adjust
! !INTERFACE:

 subroutine leap_adjust

! !DESCRIPTION:
!  Sets leap-year related variables
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!---------------------------------------------------------------------  
!
!  local variables
!
!---------------------------------------------------------------------  

   integer (int_kind) :: nm  ! dummy month index

!-----------------------------------------------------------------------
!
!  is iyear a leap year?
!
!-----------------------------------------------------------------------
 
   leapyear = is_leapyear (iyear)
 
!-----------------------------------------------------------------------
!
!  adjust the number of days in February and in the year
!
!-----------------------------------------------------------------------

   if (leapyear) then
      days_in_month(2)  = 29
      days_in_year      = days_in_leap_year
   else
      days_in_month(2)  = 28
      days_in_year      = days_in_norm_year
   endif

   seconds_in_year     = days_in_year*seconds_in_day
   hours_in_year       = days_in_year*24.0_r8
 
!-----------------------------------------------------------------------
!
!  reset the values of days_in_prior_months(imonth)
!
!-----------------------------------------------------------------------
 
   call prior_days (days_in_prior_months, days_in_month)
 
!-----------------------------------------------------------------------
!EOC

 end subroutine leap_adjust

!***********************************************************************
!BOP
! !IROUTINE: date2ymd
! !INTERFACE:

 subroutine date2ymd (date,year,month,day)

! !DESCRIPTION:
!  Decode the calendar date.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      date            ! Calendar date (integer) in yyyymmdd format

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      year,                           &! Calendar year
      month,                          &! Calendar month
      day                              ! Calendar day

!EOP
!BOC
!-----------------------------------------------------------------------

   if (.not. valid_date(date)) call exit_glc(sigAbort, &
                                             'date2ymd:invalid date')

   year  = int(     date       /10000)
   month = int( mod(date,10000)/  100)
   day   =      mod(date,  100)

!-----------------------------------------------------------------------
!EOC

 end subroutine date2ymd

!***********************************************************************
!BOP
! !IROUTINE: ymd2date
! !INTERFACE:

 subroutine ymd2date (year,month,day,date)

! !DESCRIPTION:
!  Encodes the calendar date.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      year,                          &! Calendar year
      month,                         &! Calendar month
      day                             ! Calendar day

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      date            ! Calendar date (integer) in yyyymmdd format

!EOP
!BOC
!-----------------------------------------------------------------------

   date = 10000*year + 100*month + day

!-----------------------------------------------------------------------
!EOC

 end subroutine ymd2date

!***********************************************************************
!BOP
! !IROUTINE: eday2ymd
! !INTERFACE:

 subroutine eday2ymd (eday,year,month,day)

! !DESCRIPTION:
!  Determines the year, month, and day number from elapsed days
!  since 01-01-0000.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      eday                   ! elapsed day since 01-01-0000

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      year,                           &! calendar year
      month,                          &! calendar month
      day                              ! calendar day

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(0:3) :: &
      days_each_year         ! days in year for 4-year cycle

   integer (int_kind), dimension(12)  :: &
      tdays_in_prior_months, &! temporary days in prior months
      tdays_in_month          ! temporary days in each month

   integer (int_kind) ::    &
      nm,                   &! dummy month index
      test_day,             &!
      nnorm,                &! normal year counters
      nnorm_new,            &! normal year counters
      nleap,                &! leap year counters
      nleap_new,            &! leap year counters
      cycleindex,           &! year index for 4-year cycle
      days,                 &! day counter
      ny,                   &! year index
      max_ny                 ! max possible number of years

   character (char_len) ::  &
      err_string             ! output string if error encountered

!-----------------------------------------------------------------------
!
!  If leap years are not allowed, then compute number of elapsed
!  years and the number of days in the most recent year
!
!-----------------------------------------------------------------------

   if (.not. allow_leapyear) then

!lipscomb - This section of code appears not to work.  Commented out and rewrote.
!     nnorm = eday/days_in_norm_year + 1
!      nleap = 0
!      days  = eday -nnorm*days_in_norm_year
!      tdays_in_prior_months = days_in_prior_months

      year = eday/days_in_norm_year
      days = eday - year*days_in_norm_year
      tdays_in_prior_months = days_in_prior_months

!-----------------------------------------------------------------------
!
!  Compute number of elapsed leap years and "normal" years, and
!  the number of days elapsed in the most recent year
!
!-----------------------------------------------------------------------

   else

      write(stdout,*) 'WARNING: GLINT will not handle leap years correctly!!!'

      !***  First, initialize arrays used to determine date

      days_each_year    = days_in_norm_year
      days_each_year(0) = days_in_leap_year

      tdays_in_month    = days_in_month
      tdays_in_month(2) = 29               ! year 0 value

      call prior_days (tdays_in_prior_months, tdays_in_month)

      days      = 0
      nleap     = 0
      nnorm     = 0
      nleap_new = 0
      nnorm_new = 0

      max_ny    = eday/days_in_norm_year + 1

      !*** Determine the number of elapsed years and the day number of
      !*** the present year [1,days_in_norm_year]

      year_loop: do ny = 0, max_ny

         cycleindex = mod(ny,4)
         year  = nleap + nnorm

         if (cycleindex == 0 .and. &
             (mod(ny,100) /= 0 .or. mod(ny,400) == 0)) then
            nleap_new = nleap + 1
            tdays_in_month(2) = 29
         else
            nnorm_new = nnorm + 1
            tdays_in_month(2) = 28
         endif

         !*** Update Tdays_in_prior_months for the most recent year

         call prior_days (tdays_in_prior_months, tdays_in_month)

         test_day = eday - nnorm*days_in_norm_year - &
                           nleap*days_in_leap_year

         if (test_day <= days_each_year(cycleindex) ) then
            days = test_day
            exit year_loop
         endif

         nnorm = nnorm_new
         nleap = nleap_new

      enddo year_loop

   endif ! .not. allow_leapyear

!-----------------------------------------------------------------------
!
!  Was the number of days this year properly determined?
!
!-----------------------------------------------------------------------

!lipscomb - Modified so that code does not abort when days = 0
 
!!   if (days <= 0 .or. days > days_in_leap_year ) then
   if (days < 0 .or. days > days_in_leap_year ) then
      err_string = char_blank
      write (err_string,'(a,i6)') & 
          'eday2ymd: days undetermined, days = ', days
      call exit_glc(sigAbort,trim(err_string))
   endif

!-----------------------------------------------------------------------
!
!  Determine the day- and month-numbers for this year
!
!-----------------------------------------------------------------------

   month = 0
   day   = 0

   month_loop: do nm = 1,11

      test_day = days - tdays_in_prior_months(nm+1)
      if (test_day < 0) then
         day   = days - tdays_in_prior_months(nm) + 1
         month = nm
         exit month_loop
      endif
   enddo month_loop

   if (month == 0) then
      day   = days - tdays_in_prior_months(12) + 1
      month = 12
      if (day == 32) then
         day   = 1
         month = 1
         year  = year + 1
      endif
   endif

   if (day == 0) day = 1

!-----------------------------------------------------------------------
!EOC

 end subroutine eday2ymd

!***********************************************************************
!BOP
! !IROUTINE: ymd2eday
! !INTERFACE:

 subroutine ymd2eday (year, month, day, eday)

! !DESCRIPTION:
!  Converts calendar date (year, month, day) to elapsed days since
!  01-01-0000.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      year,                          &! calendar year
      month,                         &! calendar month
      day                             ! calendar day

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      eday                   ! elapsed days since 01-01-0000

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind), dimension(0:3) :: &
      days_each_year          ! days in year for 4-year cycle

   integer (int_kind), dimension(12)  :: &
      tdays_in_prior_months, &! temporary days in prior months
      tdays_in_month          ! temporary days in each month

   integer (int_kind) ::     &
      nm,                    &! dummy month index
      num_leapyears           ! leap year counters

!---------------------------------------------------------------------  
!
!  If leap years are not allowed, eday computation is straightforward
!
!---------------------------------------------------------------------  

   if (.not. allow_leapyear) then
      eday = year*days_in_norm_year + &
             days_in_prior_months(month) + day - 1
 
!---------------------------------------------------------------------  
!
!  If leap years are allowed, compute the number of days elapsed
!  in prior months for *this* year and the number of elapsed
!  leap years prior to this year
!
!---------------------------------------------------------------------  

   else

      tdays_in_month    = days_in_month

      num_leapyears  = 1 +  year/4  - year/100 + year/400

      if (is_leapyear(year)) then
         tdays_in_month(2) = 29
         num_leapyears     = num_leapyears - 1
      else
         tdays_in_month(2) = 28
      endif

      call prior_days (tdays_in_prior_months, tdays_in_month)

      !***     Compute elapsed days for this date

      eday =         num_leapyears  *days_in_leap_year + &
             (year - num_leapyears )*days_in_norm_year + &
             tdays_in_prior_months(month) + day - 1

   endif ! .not. allow_leapyear

!---------------------------------------------------------------------  
!EOC

 end subroutine ymd2eday

!***********************************************************************
!BOP
! !IROUTINE: date2eday
! !INTERFACE:

 subroutine date2eday (date,eday)

! !DESCRIPTION:
!  Determine number of elapsed days since 01-01-0000 from the
!  calendar date in (integer) yyyymmdd format.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      date                   ! date in yyyymmdd format

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      eday                   ! elapsed days since 01-01-0000

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      year, month, day   ! year month day indices

!-----------------------------------------------------------------------
!
!  use existing routines for the conversion
!
!-----------------------------------------------------------------------

   if (.not. valid_date(date)) & 
      call exit_glc(sigAbort,'date2eday: invalid date')

   call date2ymd (date, year, month, day)
   call ymd2eday (year, month, day, eday)

!-----------------------------------------------------------------------
!EOC

 end subroutine date2eday

!***********************************************************************
!BOP
! !IROUTINE: eday2date
! !INTERFACE:

 subroutine eday2date (eday,date)

! !DESCRIPTION:
!  Determines calendar date in (integer) yyyymmdd format from the
!  number of elapsed days since 01-01-0000.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      eday                   ! elapsed days since 01-01-0000

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      date                   ! date in yyyymmdd format

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      year, month, day   ! year month day indices

!-----------------------------------------------------------------------
!
!  use existing routines for the conversion
!
!-----------------------------------------------------------------------

   call eday2ymd (eday, year, month, day)
   call ymd2date (year, month, day, date)

!-----------------------------------------------------------------------
!EOC

 end subroutine eday2date

!***********************************************************************
!BOP
! !IROUTINE: prior_days
! !INTERFACE:

 subroutine prior_days (days_in_prior_months_loc,days_in_month_loc)

! !DESCRIPTION:
!  Defines or resets the total number of days in prior months;
!  if leap years are allowed, this routine will be called once per 
!  year.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
 
   integer (int_kind), dimension(12), intent(in) :: &
      days_in_month_loc      ! current num of days in each month

! !OUTPUT PARAMETERS:

   integer (int_kind), dimension(12), intent(out) :: &
      days_in_prior_months_loc  ! number of days in prior months
 
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      nm                 !   local month index
 
!-----------------------------------------------------------------------

   days_in_prior_months_loc(1) = 0

   do nm=2,12
      days_in_prior_months_loc(nm) = &
      days_in_prior_months_loc(nm-1) + days_in_month_loc(nm-1)
   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine prior_days

!***********************************************************************
!BOP
! !IROUTINE: time_stamp
! !INTERFACE:

 subroutine time_stamp (option, order, date_string, time_string, beg_date)

! !DESCRIPTION:
!  Writes character strings containing the date and time stamps
!  mm/dd/yyyy and hh:mm:ss
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      option,   &! string with option for time stamp, 'now', 'last', 'range'
      order      ! string requesting the order of the date stamp ('ymd' or 'mdy')

! !INPUT/OUTPUT PARAMETERS:

   character (*), intent(inout), optional :: &
      date_string,          &! a string to fill with date stamp
      time_string,          &! a string to fill with time stamp
      beg_date               ! date string to use as first date in
                             !   'range' option
 
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (1), parameter :: &
      time_separator=':'

   character (16), parameter ::       &! format strings
      ymd_date_fmt1 = '(i4.4,2(a,i2.2))', &
      ymd_date_fmt2 = '(i4.4,2(i2.2))  ', &
      mdy_date_fmt1 = '(2(i2.2,a),i4.4)', &
      mdy_date_fmt2 = '(2(i2.2),i4.4)  ', &
      time_fmt  = '(i2.2,2(a,i2.2))'

   logical (log_kind) ::  &
      ymd,                &
      mdy
 
   integer (int_kind) :: date_len  ! length of date string

!-----------------------------------------------------------------------
!
!  initialize strings
!
!-----------------------------------------------------------------------

   if (present(date_string)) then
      date_string = ' '
   endif

   if (present(time_string)) then
      time_string = ' '
   endif

   ymd = .false.
   mdy = .false.

!-----------------------------------------------------------------------
!
!  check options
!
!-----------------------------------------------------------------------

   select case (trim(order))
      case ('ymd')
        ymd = .true.
      case ('mdy')
        mdy = .true.
      case default
      call exit_glc(sigAbort,'ERROR selecting order in subroutine time_stamp')
   end select

   select case (trim(option))

!-----------------------------------------------------------------------
!
!  present time
!
!-----------------------------------------------------------------------

   case ('now')

      if (present(date_string)) then
         if (date_separator /= ' ') then
            if (ymd) then
               write (date_string,ymd_date_fmt1) iyear , date_separator, &
                                                 imonth, date_separator, &
                                                 iday
            else if (mdy) then
               write (date_string,mdy_date_fmt1) imonth , date_separator, &
                                                 iday,    date_separator, &
                                                 iyear
            endif
         else
            if (ymd) then
               write (date_string,ymd_date_fmt2) iyear, imonth, iday
            else if (mdy) then
               write (date_string,mdy_date_fmt2) imonth, iday, iyear
            endif
         endif
      endif

      if (present(time_string)) then
         write (time_string,time_fmt) ihour  , time_separator, &
                                      iminute, time_separator, &
                                      isecond
      endif

!-----------------------------------------------------------------------
!
!  last timestep
!
!-----------------------------------------------------------------------

   case ('last')

      if (present(date_string)) then
         if (mdy) &
         call exit_glc(sigAbort,'ERROR time_stamp not supported with option=last & mdy order')

         if (date_separator /= ' ') then
            write (date_string,ymd_date_fmt1) iyear_last ,date_separator, &
                                          imonth_last,date_separator, &
                                          iday_last
         else
            write (date_string,ymd_date_fmt2) iyear, imonth, iday
         endif
      endif

      if (present(time_string)) then
         call exit_glc(sigAbort,'ERROR time_stamp not supported with option=last')
         
!        write (time_string,time_fmt) ihour_last, time_separator, &
!                                     iminute   , time_separator, &
!                                     isecond
      endif

!-----------------------------------------------------------------------
!
!  time range
!
!-----------------------------------------------------------------------

   case ('range')

      if (.not. present(beg_date)) &
         call exit_glc(sigAbort, &
                       'time_stamp: cannot compute range w/o beg date')

      if (present(date_string)) then
         date_string = trim(beg_date)/&
                                      &/'-'
         date_len = len_trim(date_string) + 1
         if (date_separator /= ' ') then
            write (date_string(date_len:),ymd_date_fmt1)              & 
                                          iyear , date_separator, &
                                          imonth, date_separator, &
                                          iday
         else
            write (date_string(date_len:),ymd_date_fmt2) iyear, imonth, iday
         endif
      endif

      if (present(time_string)) then
         call exit_glc(sigAbort,'ERROR time_stamp not supported with option=last')
!        write (time_string,time_fmt) ihour  , time_separator, &
!                                     iminute, time_separator, &
!                                     isecond
      endif

!-----------------------------------------------------------------------

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine time_stamp



!***********************************************************************
!BOP
! !IROUTINE: cesm_date_stamp
! !INTERFACE:
 subroutine cesm_date_stamp (date_string, ymds)

! !DESCRIPTION:
!-----------------------------------------------------------------------
!
!     write a character string containing the date stamp
!        yyyy-mm-dd, yyyy-mm, or yyyy
!
!-----------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

  character (*), intent(in) ::   ymds         ! a string indicating date stamp format

! !OUTPUT PARAMETERS:
  character (*), intent(out) ::   date_string  ! a string to fill with date stamp
 

!EOP
!BOC

  character (4) :: cesm_cyear
  character (2) :: cesm_cmonth
  character (2) :: cesm_cday  
  character (5) :: cesm_csecond
 
  integer (kind=int_kind) ::  &
     iyear_stamp             ,&
     imonth_stamp            ,&
     iday_stamp              ,&
     itotal_second
 
   date_string = char_blank

!---------------------------------------------------------------------
!     set ixxxx_stamp variables to conform to the cesm standard
!---------------------------------------------------------------------
   if (midnight) then
     if (eoy) then
        iyear_stamp  = iyear_last
        imonth_stamp = 12
        iday_stamp   = 31
     elseif (eom) then
        iyear_stamp  = iyear
        imonth_stamp = imonth_last
        iday_stamp   = iday_last
     elseif (eod) then 
        iyear_stamp  = iyear
        imonth_stamp = imonth
        iday_stamp   = iday_last
     endif
   else
     iyear_stamp  = iyear
     imonth_stamp = imonth
     iday_stamp   = iday
   endif
 
   select case (trim(ymds))
     case ('ymds')
!---------------------------------------------------------------------
!    use unmodified ixxx variables if printing ymds information
!---------------------------------------------------------------------
       itotal_second = isecond + 60*iminute + 3600*ihour
       call int_to_char (4,iyear        , cesm_cyear  )
       call int_to_char (2,imonth       , cesm_cmonth )
       call int_to_char (2,iday         , cesm_cday   )
       call int_to_char (5,itotal_second, cesm_csecond)
       write (date_string,1000) cesm_cyear, cesm_cmonth, cesm_cday, &
                                cesm_csecond

     case ('ymd')
        call int_to_char (4,iyear_stamp  , cesm_cyear )
        call int_to_char (2,imonth_stamp , cesm_cmonth)
        call int_to_char (2,iday_stamp   , cesm_cday  )
        write (date_string,1000) cesm_cyear, cesm_cmonth, cesm_cday

     case ('ym')
        call int_to_char (4,iyear_stamp  , cesm_cyear )
        call int_to_char (2,imonth_stamp , cesm_cmonth)
        write (date_string,1000) cesm_cyear, cesm_cmonth

     case ('y')
        call int_to_char (4,iyear_stamp  , cesm_cyear)
        write (date_string,1000) cesm_cyear
        
        case default 
        call exit_glc(sigAbort,'(cesm_date_stamp)')
 
   end select
 

 1000 format (a4,:,'-',a2:,'-',a2,:,'-',a5)
  
!EOC

!-----------------------------------------------------------------------

  end subroutine cesm_date_stamp

!***********************************************************************
!BOP
! !IROUTINE: is_near
! !INTERFACE:

 function is_near (test_value, target, tol)

! !DESCRIPTION:
!  Determines if test\_value is ``near'' the target value.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), intent(in) :: &
      test_value,           &! value to test
      target,               &! value to test against
      tol                    ! tolerance for determining nearness

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: & 
      is_near            ! result (T or F) of nearness test

!-----------------------------------------------------------------------
!
!  just a simple test...
!
!-----------------------------------------------------------------------

   if (abs(test_value - target) <= tol) then
      is_near = .true.
   else
      is_near = .false.
   endif

!-----------------------------------------------------------------------
!EOC

 end function is_near

!***********************************************************************
!BOP
! !IROUTINE: is_leapyear
! !INTERFACE:

 function is_leapyear (iyear_loc)

! !DESCRIPTION:
!  Determines if test\_year is a leapyear.
!
!  Assumptions:
!  \begin{itemize}
!  \item  year = 0 is the first year of the integration
!  \item  standard calendar has 28 days in February, a leap year has 29
!  \end{itemize}
!
!  Algorithm: a year is a leap year if it is:
!  \begin{itemize}
!  \item  divisible by 4,
!  \item  NOT divisible by 100, except if
!  \item  also divisible by 400
!  \end{itemize}
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      iyear_loc              ! input year to test for leapyear
 
! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      is_leapyear  ! logical result with true if test year is leapyear
 
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  check for leap year
!
!-----------------------------------------------------------------------

   is_leapyear = .false.

   if (allow_leapyear .and. mod(iyear_loc,4) == 0  .and.         &
       (mod(iyear_loc,100) /= 0 .or. mod(iyear_loc,400) == 0 ) ) &
      is_leapyear = .true.

!-----------------------------------------------------------------------
!EOC

 end function is_leapyear

!***********************************************************************
!BOP
! !IROUTINE: valid_date
! !INTERFACE:

 function valid_date (date)

! !DESCRIPTION:
!  Determines if a valid year, month & day can be decoded
!  from the calendar date in (integer) yyyymmdd format.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      date                   ! calendar date in yyyymmdd format

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      valid_date         ! logical return value = true if valid date

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      year, month, day   ! year month day indices

!-----------------------------------------------------------------------
!
!  check a variety of possible error conditions
!
!-----------------------------------------------------------------------

   year  = int(     date       /10000)
   month = int( mod(date,10000)/  100)
   day   =      mod(date,  100)

   valid_date = .true.

   if (year  <  0) valid_date = .false.
   if (month <  1) valid_date = .false.
   if (month > 12) valid_date = .false.
   if (day   <  1) valid_date = .false.
   if (day   > days_in_month(month)) valid_date = .false.

!-----------------------------------------------------------------------
!EOC

 end function valid_date

!***********************************************************************
!BOP
! !IROUTINE: valid_ymd_hms()
! !INTERFACE:

 function valid_ymd_hms()

! !DESCRIPTION:
!  Determines if the computed values of iyear,imonth,iday,
!  ihour, iminute, and isecond are within reasonable bounds.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      valid_ymd_hms      ! logical return value = true if current
                         ! year, month, day, hour, minute, second
                         ! are withing valid ranges

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
      
   logical (log_kind) :: &
      valid_year,        &! flags to determine validity of
      valid_month,       &! specific values
      valid_day_b,       &!
      valid_day_e,       &!
      valid_hour,        &!
      valid_minute,      &!
      valid_second,      &!
      valid_eday_run,    &!
      valid_eday_year,   &!
      valid_feb_day       !

   character (char_len) :: err_string ! string for error message

   character (*), parameter :: &
      err_fmt = '(a,i6)'       ! format for error string

!-----------------------------------------------------------------------
!
!  set default condition of true for all flags
!
!-----------------------------------------------------------------------

   valid_ymd_hms     = .true.
   valid_year        = .true.
   valid_month       = .true.
   valid_day_b       = .true.
   valid_day_e       = .true.
   valid_hour        = .true.
   valid_minute      = .true.
   valid_second      = .true.
   valid_eday_run    = .true.
   valid_eday_year   = .true.
   valid_feb_day     = .true.
 
!-----------------------------------------------------------------------
!
!  check a variety of possible error conditions
!
!-----------------------------------------------------------------------

   if (iyear < 0) then
      valid_ymd_hms = .false.
      valid_year    = .false.
   endif

   if (imonth < 0 .or. imonth > 12) then
      valid_ymd_hms = .false.
      valid_month   = .false.
   endif

   if (iday < 1) then
      valid_ymd_hms = .false.
      valid_day_b   = .false.
   endif
 
   if (valid_ymd_hms) then   ! prevents out-of-range reference
      if (iday > days_in_month(imonth)) then
         valid_ymd_hms = .false.
         valid_day_e   = .false.
      endif
   endif
 
   if (ihour < 0 .or. ihour > 24) then
      valid_ymd_hms = .false.
      valid_hour    = .false.
   endif
 
   if (iminute < 0 .or. iminute > 60) then
      valid_ymd_hms = .false.
      valid_minute  = .false.
   endif
 
   if (isecond < 0 .or. isecond > 60) then
      valid_ymd_hms = .false.
      valid_second  = .false.
   endif

   if (elapsed_days_this_year < 0) then
      valid_ymd_hms   = .false.
      valid_eday_year = .false.
   endif
 
   if (elapsed_days_init_date < 0) then
      valid_ymd_hms   = .false.
      valid_eday_run  = .false.
   endif
 
   if (.not. allow_leapyear .and. imonth == 2 .and. iday == 29) then
      valid_ymd_hms = .false.
      valid_feb_day = .false.
   endif
 
!-----------------------------------------------------------------------
!
!  if errors detected, write out message and quit
!
!-----------------------------------------------------------------------

   if (.not. valid_ymd_hms) then

      err_string = char_blank

      if (.not. valid_year) &
         write(err_string,err_fmt) &
              'Invalid date (iyear must be > 0 ): iyear = ', iyear
 
      if (.not. valid_month) &
         write(err_string,err_fmt) &
              'Invalid date ( imonth must be in [1,12] ): imonth = ', &
                                                          imonth
 
      if (.not. valid_day_b) &
         write(err_string,err_fmt) &
              'Invalid date (iday must be greater than 1): iday = ',iday
 
      if (.not. valid_day_e) &
         write(err_string,err_fmt) &
              'Invalid date (iday must be less than days_in_month):'/&
              &/' iday = ',iday 
 
      if (.not. valid_hour) &
         write(err_string,err_fmt) &
              'Invalid date (ihour must be in [0,23] ): ihour = ', ihour
 
      if (.not. valid_minute) &
         write(err_string,err_fmt) &
              'Invalid date (iminute must be in [0,59] ): iminute = ', &
                                                          iminute
 
      if (.not. valid_second) &
         write(err_string,err_fmt) &
              'Invalid date (isecond must be in [0,59] ): isecond = ', &
                                                          isecond
 
      if (.not. valid_eday_run) &
         write(err_string,err_fmt) &
              'Invalid date (elapsed_days_init_date must be > 0 ) ', &
                             elapsed_days_init_date
 
      if (.not. valid_eday_year) &
         write(err_string,err_fmt) &
              'Invalid date (elapsed_days_this_year must be > 0) ', &
                             elapsed_days_this_year
 
      if (.not. valid_feb_day) &
         write(err_string,*) &
              ' Error: initial date contains leap day '/&
              &/' but no leap years are allowed.', iday
 
      call exit_glc(sigAbort,trim(err_string))
 
   endif   ! valid_ymd_hms
 
!-----------------------------------------------------------------------
!EOC

 end function valid_ymd_hms

!***********************************************************************
!BOP
! !IROUTINE: write_time_manager_options
! !INTERFACE:

 subroutine write_time_manager_options

! !DESCRIPTION:
!  Writes all time manager options to stdout.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
! !OUTPUT PARAMETERS:
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      k, ind, nn

   character (1) :: &
      suffix
 
   character (2) ::   &
      cmonth_end_run, &!
      cday_end_run,   &!
      cmonth0,        &!
      cday0            !
 
   character (4) ::   &
      cyear0,         &!
      cyear_end_run    !

   character (*), parameter :: &! output formats
      out_fmt1 = "('       date(month-day-year):',2x,2(a2,'-'),a4)", &
      out_fmt2 = "('                    ',a7,2x,i10)",               &
      out_fmt3 = "('This run will terminate ',/,a)",                 &
      out_fmt4 = "(a, :i7, a,a)",                                    &
      out_fmt5 = "('                   ',a8,2x,f16.3)",              &
      out_fmt6 = "('GLC time step = ',1pe12.6, ' seconds')"

!-----------------------------------------------------------------------
!
!  write only from master task
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then

!-----------------------------------------------------------------------
!
!     write start/current time data
!
!-----------------------------------------------------------------------

      call int_to_char(4, iyear0 , cyear0)
      cmonth0 = cmonths(imonth0)
      cday0   = cdays  (iday0)

      call int_to_char(4, iyear_end_run  , cyear_end_run )
      cmonth_end_run = cmonths(imonth_end_run)
      cday_end_run   = cdays  (iday_end_run)

      write (stdout,blank_fmt)
      write (stdout,ndelim_fmt)
      write (stdout,blank_fmt)
      write (stdout,'(a23)') 'Time Information'
      write (stdout,blank_fmt)
      write (stdout,delim_fmt)

      write (stdout,blank_fmt)
      write (stdout,'(a8,a)') 'Run id: ',trim(runid)

      write (stdout,blank_fmt)
      write (stdout,'(a28)') 'This simulation started from'
      write (stdout,out_fmt1) cmonth0, cday0, cyear0
      write (stdout,out_fmt2) '  hour:', ihour0
      write (stdout,out_fmt2) 'minute:', iminute0
      write (stdout,out_fmt2) 'second:', isecond0

      write (stdout,blank_fmt)
      write (stdout,'(a28)') 'This run        started from'
      write (stdout,out_fmt1) cmonth, cday, cyear
      write (stdout,out_fmt2) '  hour:', ihour
      write (stdout,out_fmt2) 'minute:', iminute
      write (stdout,out_fmt2) 'second:', isecond

      if (nsteps_total /=   0) &
         write (stdout,out_fmt2) '  step:', nsteps_total

      write (stdout,blank_fmt)
 
      if (end_run_at_midnight) then
         write (stdout,out_fmt3)  'at 00:00:00 on'
      else if (dtt > seconds_in_day) then
         write (stdout,out_fmt3)  'at the end of the day on or after'
      else
         write (stdout,out_fmt3)  'at the end of the day on'
      endif
 
      if (stop_count == 1) then
         suffix = ' '
      else
         suffix = 's'
      endif
 
      write (stdout,out_fmt1) cmonth_end_run,cday_end_run,cyear_end_run
 
      select case (stop_option)        
 
      case ('never')
         write (stdout,out_fmt4) 'upon receipt of stop signal' /&
                               &/ ' from external source (eg, cpl)'
      case ('nyear')
         write(stdout,out_fmt4) 'after running for ',stop_count, &
                                ' year', suffix
      case ('nmonth')
         write(stdout,out_fmt4) 'after running for ',stop_count, &
                                ' month', suffix
      case ('nday')
         write(stdout,out_fmt4) 'after running for ',stop_count, &
                                ' day', suffix
      case ('eoy')
         write(stdout,out_fmt4) 'at the end of the year after ', &
                                stop_count, ' year', suffix
      case ('eom')
         write(stdout,out_fmt4) 'at the end of the month after', &
                                stop_count,  ' month', suffix
      case ('eod')
         write (stdout,out_fmt4) 'at the end of the day'
      case ('nstep','nsteps')
         write (stdout,out_fmt4) 'after ', stop_count,' timestep', &
                                 suffix
      case ('date')
         write (stdout,out_fmt4) 'after reaching the specified date'
      case default
      end select
      write (stdout,'(a63)') 'unless a stop signal is received'/&
                          &/' from external source (eg, cpl)'

      write (stdout,blank_fmt)
      write (stdout,'(a28)') 'Starting elapsed time in    '
      write (stdout,out_fmt5) '  years:', tyear
      write (stdout,out_fmt5) ' months:', tmonth
      write (stdout,out_fmt5) '   days:', tday
      write (stdout,out_fmt5) '  hours:', thour
      write (stdout,out_fmt5) 'seconds:', tsecond

!-----------------------------------------------------------------------
!
!     timestep information
!
!-----------------------------------------------------------------------

      write (stdout,blank_fmt)
      if (dt_option == 'auto_dt') then
         write (stdout,'(a45)') &
            'Automatic time step option (auto_dt)  enabled'
      else
         write (stdout,'(a45)') &
            'Automatic time step option (auto_dt) disabled'
      endif

      write (stdout,blank_fmt)
      write (stdout,'(a11,1pe12.6)') 'dt_count = ',dt_count
 
      write (stdout,blank_fmt)
      write (stdout,out_fmt6) dtt

!-----------------------------------------------------------------------
!
!     other options
!
!-----------------------------------------------------------------------

      write (stdout,blank_fmt)
      if (allow_leapyear) then
         write (stdout,'(a22)') 'Leap years     allowed'
      else
         write (stdout,'(a22)') 'Leap years not allowed'
      endif

!-----------------------------------------------------------------------
!
!  end of writes
!
!-----------------------------------------------------------------------

   endif ! (my_task == master_task)

!-----------------------------------------------------------------------
!EOC

 end subroutine write_time_manager_options

!***********************************************************************
!BOP
! !IROUTINE: int_to_char
! !INTERFACE:

 subroutine int_to_char(string_length, int_in, char_out)

! !DESCRIPTION:
!  Converts an integer into a character with a requested length and
!  pads spaces with zeroes.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: & 
      string_length,     &! length of desired output character string
      int_in              ! input integer to be converted

! !OUTPUT PARAMETERS:

   character(string_length), intent(out) :: & 
      char_out            ! character equivalent of input integer

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n,                 &! dummy counter
      ifact,             &! factor of 10 for picking off digits
      iquot, iremaind     ! quotient, remainder for division by ifact

!-----------------------------------------------------------------------
!
!  convert to string by picking off one digit at a time and writing
!  it into a character string
!
!-----------------------------------------------------------------------

   iremaind = int_in

   do n=1,string_length
      ifact    = 10**(string_length - n)   ! power of 10 for leftmost
      iquot    = iremaind/ifact            ! compute leftmost digit
      iremaind = iremaind - iquot*ifact    ! remove digit for next pass

      write(char_out(n:n),'(i1)') iquot    ! write digit to string
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine int_to_char

!***********************************************************************

 end module glc_time_management

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
