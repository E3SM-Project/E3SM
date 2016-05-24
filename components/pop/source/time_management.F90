!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module time_management

!BOP
! !MODULE: time_management

! !DESCRIPTION:
!  This module contains a large number of routines for calendar, time
!  flags and other functions related to model time.
!
! !REVISION HISTORY:
!  SVN:$Id: time_management.F90 26358 2011-01-11 18:10:13Z njn01 $
!
! !USES:

   use kinds_mod
   use constants
   use blocks
   use domain_size
   use domain
   use broadcast
   use grid
   use io
   use io_tools
   use exit_mod
   use registry

   implicit none
   public
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_time1,            &
             init_time2,            &
             time_manager,          &
             init_time_flag,        &
             access_time_flag,      &
             get_time_flag_id,      &
             override_time_flag,    &
             eval_time_flag,        &
             check_time_flag,       &
             check_time_flag_int,   &
             time_to_start,         &
             time_stamp,            &
             int_to_char,           &
             ccsm_date_stamp,       &
             ccsm_char_date_and_time


! !PUBLIC DATA TYPES:

   type time_flag

      character (char_len) ::  &
         name,                 &! name for flag
         owner                  ! name of initializing routine

      logical (log_kind) ::    &
         value,                &! logical state of flag
         old_value,            &! last state of flag
         default,              &! default state of flag
         has_default,          &! T if default defined, F if no default
         has_offset_date,      &! T if flag has reference-date/offset value; F if not
         is_initialized,       &! T if flag has been initialized via init_time_flag; F if not
         is_reserved,          &! T if flag has been "reserved" by non-owner routine
         done                   ! if freq_opt is once, then done is true after one-time event is done
                                !    requires developer to manually set "done"

      integer (int_kind) ::    &
         freq_opt,             &! frequency units for switching flag on
         freq,                 &! freq in above units for switching flag
         offset_year,          &! offset/reference year  (input;optional)
         offset_month,         &! offset/reference month (input;optional)
         offset_day,           &! offset/reference day   (input;optional)
         eyears,               &! elapsed years  between 01-01-0000 and ref date (computed)
         emonths,              &! elapsed months between 01-01-0000 and ref date (computed)
         edays                  ! elapsed days   between 01-01-0000 and ref date (computed)

   end type

! !PUBLIC DATA MEMBERS:

!-----------------------------------------------------------------------
!
!  variables for run control
!
!-----------------------------------------------------------------------

   character (char_len) :: &
      stop_option         ,&! specify how to determine stopping time
      runid               ,&! an identifier for the run
      dt_option             ! method to determine tracer timestep size

   integer (int_kind) :: &
      stop_count        ,&! num of stop_option intervals before stop
                          !   OR date (yyyymmdd) at which model stops
      stop_iopt         ,&! integer value for stop_option
      nsteps_total      ,&! steps (full&half) since beginning of run sequence
      nsteps_run        ,&! steps taken since beginning of this run
      nsteps_per_day    ,&! integer number of steps per day
      len_runid           ! length of character runid

   integer (int_kind) ::      &! variables used with avgfit
      fit_freq              , &! num of intervals/day into which full &
                               !  half timesteps must exactly "fit"
      fullsteps_per_interval, &! num of full timesteps per fitting interval
      halfsteps_per_interval, &! num of half timesteps per fitting interval
      fullsteps_per_day     , &! num of full timesteps per day
      halfsteps_per_day     , &! num of half timesteps per day 
      nsteps_per_interval   , &! number of steps in each 'fit' interval
      nsteps_this_interval     ! number of steps in current 'fit' interval

   integer (int_kind), private :: &
      stop_now            ,&! time_flag id for stopping
      coupled_ts            ! time_flag id for a coupled timestep

   logical (log_kind)   :: &! this timestep is:
      adjust_year         ,&!   step at which year values updated
      eod                 ,&!   at the end of the day
      eom                 ,&!   at the end of the month
      eoy                 ,&!   at the end of the year
      avg_ts              ,&!   an averaging timestep
      back_to_back        ,&!   the second of two avg timesteps in row
      f_euler_ts          ,&!   a forward Euler timestep (first ts)
      first_step          ,&!   first time step
      leapfrogts          ,&!   a leapfrog timestep
      matsuno_ts          ,&!   an Euler-backward timestep
      midnight            ,&!   at midnight
      ice_ts              ,&!   an ice-formation timestep
      sample_qflux_ts       !   time to sample qflux for time avg 

   logical (log_kind)   :: &! the last timestep was:
      eod_last            ,&!   at the end of the day
      eom_last            ,&!   at the end of the month
      eoy_last            ,&!   at the end of the year
      midnight_last       ,&!   at midnight
      avg_ts_last           !   an averaging timestep

   logical (log_kind)   :: &! the next timestep is:
      adjust_year_next    ,&!   step at which year values updated
      eom_next            ,&!   at the end of the month
      midnight_next       ,&!   at midnight
      avg_ts_next         ,&!   an averaging ts?
      back_to_back_next   ,&!   a second avg step in a row
      end_run_at_midnight ,&! does model run end at midnight
      new_dtt_value         ! does restart have a new step size
 
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
      isecond_last        ,&! second  [0,59]          |
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
      days_in_prior_year      ,&! days in prior   year 
      elapsed_days            ,&! full days elapsed since    01-01-0000
      elapsed_days0           ,&! full days elapsed between  01-01-0000 
                                !      and initial start date
      elapsed_days_jan1       ,&! full days elapsed prior to 01-01-iyear
      elapsed_days_this_run   ,&! full days elapsed since beginning of
                                !                   this segment of run
      elapsed_days_this_year  ,&! full days elapsed since beginning of yr
      elapsed_days_init_date  ,&! full days elapsed since initial time
      elapsed_days_end_run    ,&! full days elapsed from 01-01-0000 to end
                                !                   of this run
      elapsed_days_max        ,&! maximum number of full days allowed 
      elapsed_months          ,&! full months elapsed since   01-01-0000
      elapsed_months0         ,&! full months elapsed between 01-01-0000
                                !       and initial start date
      elapsed_months_this_run ,&! full months elapsed since beginning of
                                !                     this segment of run
      elapsed_months_init_date,&! full months elapsed since initial start date
      elapsed_years           ,&! full years  elapsed since   01-01-0000
      elapsed_years0          ,&! full years  elapsed between 01-01-0000 
                                !      and initial start date
      elapsed_years_this_run  ,&! full years  elapsed since beginning of
                                !      this segment of run
      elapsed_years_init_date   ! full years  elapsed since initial start date
 
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
      newday            ,&!
      newhour           ,&!
      newsecond         ,&!
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
      thour00                  ,&!
      thour00_begin_this_year

   real (r8), dimension(12)  :: &
      thour00_midmonth_calendar,&! num hours to middle of calendar month
      thour00_endmonth_calendar,&! num hours to end of calendar month
      thour00_midmonth_equal   ,&! num hours to middle of equal-spaced month
      thour00_endmonth_equal     ! num hours to end of equal-spaced month

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
      freq_opt_nstep    = 6,        &
      freq_opt_once     = 7

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

   logical (log_kind) :: &
      laccel              ! flag for acceleration

   real (r8) ::          &
      dt_count          ,&! input count to determine dtt
      dtt               ,&! tracer timestep (sec)
      dtt_input         ,&! tracer timestep (sec) as specified in namelist
                          !   input; may be different from restart value
      dtu               ,&! momentum timestep (sec)
      dtp               ,&! barotropic timestep (sec)
      c2dtu             ,&!
      c2dtp             ,&!
      c2dtq             ,&!
      dtuxcel           ,&! factor to multiply MOMENTUM timestep
      stepsize          ,&! size of present timestep (sec)
      stepsize_next       ! size of next timestep (sec)

   real (r8), dimension(km) :: &
      dttxcel           ,&! array for depth-dependent acceleration
      dt                ,&! time step at each level
      c2dtt             ,&
      dztxcel           ,&
      dzwxcel

!-----------------------------------------------------------------------
!
!  time-centering and mixing variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      impcor              ! implicit treatment of Coriolis terms

   integer (int_kind), parameter :: &
      tmix_matsuno = 1,  &! use matsuno step for time mixing
      tmix_avg     = 2,  &! use averaging step for time mixing
      tmix_avgbb   = 3,  &! use averaging step for time mixing, with
                          !  back_to_back option to keep time boundaries
      tmix_avgfit  = 4    ! use averaging step for time mixing, 
                          !  selecting the timestep size in such
                          !  a way as to force the end of the day
                          !  (or interval) to coincide with the end of 
                          !  a timestep

   integer (int_kind) :: &
      tmix_iopt,         &! option for which time mixing to use
      time_mix_freq,     &! frequency of mixing
      mix_pass            ! number of passes to perform mixing

   real (r8) :: &
      beta      ! = {alpha,theta} on {leapfrog,Matsuno} steps

   real (r8), parameter ::  &
      alpha = c1/c3,        &! leapfrog grap(ps) time-centering param
      theta = p5,           &! Matsuno grap(ps) time-centering param
      gamma = c1 - c2*alpha  ! for geostrophic balance, otherwise
                             ! coriolis and  surface-pressure gradient
                             ! are not time centered

!-----------------------------------------------------------------------
!
!  variables documenting avgfit progression
!
!-----------------------------------------------------------------------

   logical (kind=log_kind), dimension(:), allocatable :: &
      interval_avg_ts

   real (r8), dimension(:), allocatable :: &
      interval_cum_dayfrac

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  a few private variables used only internally
!
!-----------------------------------------------------------------------

   private :: time_to_do  !access via time-flags only

   character (char_len), private ::  exit_string

   real (r8), private ::          &
      rhour_next,        &! rhour   for next timestep
      rminute_next,      &! rminute for next timestep
      rsecond_next        ! rsecond for next timestep

   integer (int_kind), private, allocatable, dimension(:) ::  &
      hour_at_interval,  &! hour value at end of each interval
      min_at_interval,   &! min  value at end of each interval
      sec_at_interval     ! sec  value at end of each interval

   logical (kind=log_kind),private   ::   &
      debug_time_management   = .false. 

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
      k,                 &! vertical level index
      nu,                &! i/o unit number
      nm                  ! month index

   logical (log_kind)   ::  &
      last_step, error_code ! stepsize error testing

   character (char_len) ::  &
      stepsize_string, message_string ! last step/error condition reporting string
!-----------------------------------------------------------------------
!
!  namelist input
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      nml_error           ! namelist i/o error flag

   character (char_len) :: &
      time_mix_opt,        &! option for time mixing (Matsuno,averaging)
      accel_file,          &! file containing acceleration factors
      message

   namelist /time_manager_nml/                                   &
               runid,          time_mix_opt,    time_mix_freq,   &
               impcor,         laccel,          accel_file,      &
               dtuxcel,        iyear0,          imonth0,         &
               iday0,          ihour0,          iminute0,        &
               isecond0,       dt_option,       dt_count,        &
               stop_option,    stop_count,      date_separator,  &
               allow_leapyear, fit_freq



!-----------------------------------------------------------------------
!
!  Define and initialize time flags
!
!-----------------------------------------------------------------------

   call init_time_flag  ('stop_now' , stop_now,   owner='init_time1', default=.false.)
   call access_time_flag('coupled_ts',coupled_ts)

   call reset_switches

!-----------------------------------------------------------------------
!
!  set initial values  for namelist inputs
!
!-----------------------------------------------------------------------

   runid            = 'unknown_runid'
   impcor           = .true.
   time_mix_opt     = 'avgfit'
   fit_freq         =  1
   time_mix_freq    = 17
   dtuxcel          = c1
   laccel           = .false.
   accel_file       = 'unknown_accel_file'
   allow_leapyear   = .false.
   stop_option      = 'unknown_stop_option'
   stop_count       = -1
   dt_option        = 'auto_dt'
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
      exit_string = 'FATAL ERROR: reading time_manager_nml'
      call document ('init_time1', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
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

   if (my_task == master_task) then
      select case (time_mix_opt)
      case ('matsuno')
         tmix_iopt = tmix_matsuno
      case ('avg')
         tmix_iopt = tmix_avg
      case ('avgbb')
         tmix_iopt = tmix_avgbb
      case ('avgfit')
         tmix_iopt = tmix_avgfit
      case default
         tmix_iopt = -1000
      end select
   endif

   call broadcast_scalar (runid           , master_task)
   call broadcast_scalar (tmix_iopt       , master_task)
   call broadcast_scalar (fit_freq        , master_task)
   call broadcast_scalar (time_mix_freq   , master_task)
   call broadcast_scalar (impcor          , master_task)
   call broadcast_scalar (laccel          , master_task)
   call broadcast_scalar (dtuxcel         , master_task)
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

   if (tmix_iopt == -1000) then
      exit_string = 'FATAL ERROR: unknown option for time mixing'
      call document ('init_time1', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

   if (tmix_iopt == tmix_matsuno) then
      exit_string = 'FATAL ERROR:  matsuno time-mixing option is not supported in CCSM; ' /&
              &/' budget diagnostics are incorrect with matsuno and tavg may be incorrect'
      call document ('init_time1', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif


   len_runid = len_trim(runid)

!-----------------------------------------------------------------------
!
!  determine the value for dtt, based upon model input parameters
!
!-----------------------------------------------------------------------

   select case (dt_option)
 
   case('auto_dt')   
      !*** scale tracer timestep  dt = 1 hr at dx = 2 degrees
      dtt = seconds_in_hour*(180.0_r8/float(nx_global))
      steps_per_day  = seconds_in_day/dtt
      steps_per_year = steps_per_day*days_in_norm_year

   case('steps_per_year')
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
      exit_string = 'FATAL ERROR: unknown dt_option'
      call document ('init_time1', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   end select
 
!-----------------------------------------------------------------------
!
!  modify dtt value, if using avgfit option
!
!-----------------------------------------------------------------------

   if (tmix_iopt  == tmix_avgfit) then
 
     !-----------------------------------------------------------------
     !*** determine the number of full, half, and total number of
     !*** steps in each interval. 
     !*** fit_freq = 1 ==> coupling interval is one day
     !*** fit_freq > 1 ==> coupling every day/fit_freq
     !*** fit_freq < 1 ==> coupling every abs(fit_freq)*day  (not supported)
     !-----------------------------------------------------------------
   
     !*** small values of time_mix_freq are not supported

     if (time_mix_freq <= 3) then
      exit_string = 'FATAL ERROR: time_mix_freq must be > 3'
      call document ('init_time1', exit_string)
      write(exit_string,'(a,1x,i3)') 'FATAL ERROR: unsupported time_mix_freq value: ', time_mix_freq
      call document ('init_time1', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
     endif

     nsteps_per_day         = steps_per_day
     fullsteps_per_interval = nsteps_per_day/fit_freq

     !*** require at least one full step per interval

     if (fullsteps_per_interval < 1) fullsteps_per_interval = 1

     halfsteps_per_interval = (time_mix_freq +fullsteps_per_interval)/  &
                              (time_mix_freq-1)
     nsteps_per_interval    = fullsteps_per_interval +  halfsteps_per_interval

     !*** do not allow last timestep in interval to be a half step

     if (mod(nsteps_per_interval,time_mix_freq) == 0) then
      fullsteps_per_interval = fullsteps_per_interval + 1
      halfsteps_per_interval = (time_mix_freq +fullsteps_per_interval)/ &
                               (time_mix_freq-1)
      nsteps_per_interval    = fullsteps_per_interval +  halfsteps_per_interval
        call document ('init_time1',   &
       'fullsteps_per_interval was adjusted to avoid ending interval on half step')
     endif

     !*** do not allow last timestep in interval to be a half step -- 
     !*** a special case

     if (fullsteps_per_interval == 1 .and. halfsteps_per_interval == 1) then
        fullsteps_per_interval = fullsteps_per_interval + 1
        nsteps_per_interval    = fullsteps_per_interval +  halfsteps_per_interval
        call document ('init_time1',   &
       'fullsteps_per_interval was adjusted to avoid ending interval on half step')
     endif
 
     !*** determine the number of half, full, and total steps in
     !*** each day

     fullsteps_per_day = fit_freq*fullsteps_per_interval
     halfsteps_per_day = fit_freq*halfsteps_per_interval
     nsteps_per_day    = fullsteps_per_day + halfsteps_per_day

     !*** compute modified dtt value
     dtt = seconds_in_day/(fullsteps_per_day + 0.5*halfsteps_per_day)
     steps_per_day  = seconds_in_day/dtt

     !***  allocate arrays used for testing
      allocate (hour_at_interval(fit_freq), min_at_interval(fit_freq),  &
                sec_at_interval(fit_freq))

     !***  test and document tmix_avgfit timestep size and time-stepping logic
     call test_timestep
   else
     nsteps_per_interval = steps_per_day
   endif ! tmix_avgfit


   dtt_input = dtt
   dtp       = dtt
   dtu       = dtt

!-----------------------------------------------------------------------
!
!  multiply tracer timestep by acceleration factor(s)
!  for depth varying acceleration factors, create several arrays
!     to store the variable timestep (dt) and some
!     combination dz/dttxcel arrays for use in convective adjustment.
!
!-----------------------------------------------------------------------

   if (laccel) then

      call get_unit(nu)
      if (my_task == master_task) then

         open(nu, file=accel_file, status = 'old')
         do k = 1,km
            read(nu,*) dttxcel(k)
         end do

         if (dttxcel(1) /= c1) then   !  make sure no accel in top layer
            write(stdout,'(a36)') 'no acceleration allowed in top layer'
            write(stdout,'(a36)') 'resetting acceleration factor to 1.0'
            dttxcel(1) = c1
         endif

         close(nu)
      endif
      call release_unit(nu)

      call broadcast_array(dttxcel, master_task)

   else

      dttxcel = c1

   endif

   do k = 1,km
      dt(k) = dtt*dttxcel(k)
      dztxcel(k) = dz(k)/dttxcel(k)
   enddo
   do k = 1,km-1
      dzwxcel(k) = c1/(dztxcel(k) + dztxcel(k+1))
   enddo
   dzwxcel(km) = c0

   dtu = dtu*dtuxcel
   dtp = dtp*dtuxcel
   if (dtuxcel /= c1 .and. my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,'(a39)') '  MOMENTUM TIMESTEP ACCELERATION ACTIVE'
   endif

!-----------------------------------------------------------------------
!
!  set initial values; some of these may be overwritten by
!  restart input
!
!-----------------------------------------------------------------------

   eom_next          = .false.  
   eom_last          = .false.  
   eod_last          = .false.  
   eoy_last          = .false.  
   midnight_last     = .false.  

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

!-----------------------------------------------------------------------
!
!  Register init_time1
!
!-----------------------------------------------------------------------
   call register_string('init_time1')

   call flushm (stdout)
!EOC

 end subroutine init_time1

!***********************************************************************
!BOP
! !IROUTINE: test_timestep
! !INTERFACE:

 subroutine test_timestep

! !DESCRIPTION:
!  * Tests the timestep computations for the avgfit option. 
!  * Documents the full and half steps taken for one full day.
!  * Identifies the hour, minute, and second at the end of each
!     ocean coupling interval
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
   real (r8)     ::  &
      dummy

   integer (int_kind)   ::  &
      nfit,nn              ! dummy indices


   logical (log_kind)   ::  &
      last_step, error_code ! stepsize error testing

   character (char_len) ::  &
      stepsize_string, message_string ! last step/error condition reporting string

!-----------------------------------------------------------------------
!
!  initialize nstep and seconds counters
!
!-----------------------------------------------------------------------
    nsteps_this_interval = 0
    nsteps_total = 0
    seconds_this_day = c0
    seconds_this_day_next = c0
 
!-----------------------------------------------------------------------
!
!  document avgfit options
!
!-----------------------------------------------------------------------
    nsteps_this_interval = 0
    if (my_task == master_task) then
      write(stdout,blank_fmt) 
      write(stdout,*) 'The following information documents the time-step '
      write(stdout,*) 'distribution for one full model day. '
      write(stdout,blank_fmt) 
      write(stdout,1100) ' time_mix_freq           ', time_mix_freq
      if (fit_freq == 1) then
        write(stdout,1100) ' fullsteps_per_day       ', fullsteps_per_day
        write(stdout,1100) ' halfsteps_per_day       ', halfsteps_per_day
        write(stdout,1100) ' nsteps_per_day          ', nsteps_per_day
      else
        write(stdout,1100) ' intervals per day       ', fit_freq
        write(stdout,1100) ' fullsteps_per_interval  ', fullsteps_per_interval
        write(stdout,1100) ' halfsteps_per_interval  ', halfsteps_per_interval
        write(stdout,1100) ' nsteps_per_interval     ', nsteps_per_interval
        write(stdout,blank_fmt) 
        write(stdout,1100) ' fullsteps_per_day       ', fullsteps_per_day    
        write(stdout,1100) ' halfsteps_per_day       ', halfsteps_per_day    
        write(stdout,1100) ' nsteps_per_day          ', nsteps_per_day
      endif
      write(stdout,1101) ! write column header information
    endif

!-----------------------------------------------------------------------
!
!  cycle through one simulated day, documenting each full/half timestep
!    and testing the last timestep of the interval for correctness
!
!  store values of avg_ts and cum_stepsize for later use by qsw cosz option
!
!-----------------------------------------------------------------------

    if (.not. allocated(interval_avg_ts)) then
      allocate(interval_avg_ts(nsteps_per_interval))
      allocate(interval_cum_dayfrac(0:nsteps_per_interval))
    endif
    interval_cum_dayfrac(0) = c0

    do nfit = 1, fit_freq
    do nn = 1, nsteps_per_interval

      call reset_switches

      nsteps_total = nsteps_total + 1
      nsteps_this_interval = nsteps_this_interval + 1
      if (nsteps_this_interval > nsteps_per_interval) nsteps_this_interval = 1
 
      call set_switches

      if (avg_ts) then
        stepsize = p5*dtt
      else
        stepsize = dtt
      endif

      if (nfit == 1) then
        interval_avg_ts(nn) = avg_ts
        interval_cum_dayfrac(nn) = interval_cum_dayfrac(nn-1) + &
           stepsize / seconds_in_day
      endif

      if (nsteps_total == 1) then
         seconds_this_day = seconds_this_day + stepsize
      else
         seconds_this_day = seconds_this_day_next
      endif

      if (avg_ts_next) then
        stepsize_next = p5*dtt
      else
        stepsize_next = dtt
      endif

      seconds_this_day_next = seconds_this_day + stepsize_next 

      !*** label last step 
      if (nn == nsteps_per_interval ) then
        last_step = .true.
        write(message_string,'(a,i3)')  '<-- LAST STEP in coupling interval ', nfit
        if (nfit == fit_freq) message_string = '<-- LAST STEP in day'
      else
        last_step = .false.
        message_string = ' '
      endif

      !*** label full and half steps
      if (avg_ts) then
        stepsize_string = 'H'
      else
        stepsize_string = 'F'
      endif

      !*** flag any errors
      error_code = .false.
      if (last_step) then
        if (avg_ts) error_code = .true.
        if (seconds_this_day > nfit*86400.0_r8/fit_freq + 0.0000001_r8 .or.  & 
            seconds_this_day < nfit*86400.0_r8/fit_freq - 0.0000001_r8)      &
            error_code = .true.
      endif

      if (error_code) exit_string = trim(message_string)//' -- ERROR!'
       
      if (my_task == master_task) then 
        write(stdout,1102) nsteps_total, trim(stepsize_string), seconds_this_day,   &
                           trim(message_string)
         
        if (last_step .and. nfit == fit_freq) then
          write(stdout,*) ' '
          call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
        endif
      endif ! master_task

     !*** error abort
     if (error_code) then
      exit_string = 'FATAL ERROR: in timestep-size computation'
      call document ('test_timestep', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
     endif
        
    enddo ! nn

!-----------------------------------------------------------------------
!
!  identify and store hour, minute, and second at interval boundaries
!
!-----------------------------------------------------------------------
    call hms (seconds_this_day, hour_at_interval(nfit), min_at_interval(nfit),  &
              sec_at_interval(nfit), dummy,dummy,dummy)
      if(master_task == my_task )  then 
         write(stdout,1103) ' hour, min, sec at end of interval = ',         &
                         hour_at_interval(nfit), min_at_interval(nfit), &
                         sec_at_interval(nfit)
      endif
    enddo ! nfit

      
!-----------------------------------------------------------------------
!
!  reset nstep and seconds counters
!
!-----------------------------------------------------------------------
    nsteps_this_interval = 0
    nsteps_total = 0
    seconds_this_day = c0
    seconds_this_day_next = c0

    
1100 format (1x, a, ' = ', i7)
1101 format (/,5x, 'Step  ', 3x,'Full/', 8x,'Time in',/, &
               5x, 'Number', 3x,'Half ', 8x,'Seconds'/)
1102 format (1x, i6, 8x,a, F25.15:,1x,a)
1103 format (42x, a, 3i4)
 end subroutine test_timestep

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
      ndays_temp,            &! temp for number of days
      nfit

   logical (log_kind) ::     &
      leapyear_test,         &! test for leap year
      error_condition

   character (4) ::          &
      cyear_end_run           ! character version of ending year

   character (2) ::          &
      cmonth_end_run,        &! character version of ending month
      cday_end_run            ! character version of ending day

   character (*), parameter :: &
      date_fmt = "(a17, 2x, a4,'-',a2,'-',a2)"
 
!-----------------------------------------------------------------------
!
!   register init_time2, then determine if both init_time1 
!   and init_ts have been registered
!
!-----------------------------------------------------------------------
   call register_string   ('init_time2')
   call registry_err_check('init_time1', .true., &
                            caller='init_time2' )
   call registry_err_check('init_ts', .true., &
                            caller= 'init_time2' )

!-----------------------------------------------------------------------
!
!  determine the number of days, months, and years elapsed since
!  01-01-0000  and number of days and months elapsed this year
!
!-----------------------------------------------------------------------

   call ymd2eday (iyear , imonth , iday , elapsed_days     )
   call ymd2eday (iyear , 1      , 1    , elapsed_days_jan1)

   elapsed_years0   = iyear0
   elapsed_months0  = elapsed_years0*12 + imonth0 - 1
   call ymd2eday      (iyear0, imonth0, iday0, elapsed_days0)

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
!  determine tsecond, the number of seconds elapsed since beginning
!  of complete simulation.
!
!-----------------------------------------------------------------------

   tsecond = elapsed_days_init_date*seconds_in_day + seconds_this_day

!-----------------------------------------------------------------------
!
!  compare the value of dtt selected for this run via namelist input
!  with that from the previous run.  if different, time_manager
!  will execute as if it were not a restart.
!
!-----------------------------------------------------------------------

   if (dtt /= dtt_input ) then
      new_dtt_value = .true.
      dtt = dtt_input
   endif

!-----------------------------------------------------------------------
!
!  determine if this is a leap year; set days_in_year and 
!  days_in_prior_months, regardless of value of allow_leapyear
!
!-----------------------------------------------------------------------

   call leap_adjust

   if (iyear > 0 .and. is_leapyear(iyear-1) ) then
      days_in_prior_year = days_in_leap_year
   else
      days_in_prior_year = days_in_norm_year
   endif

!-----------------------------------------------------------------------
!
!  check initial iday, imonth, etc values for reasonableness
!
!-----------------------------------------------------------------------

   if ( .not. valid_ymd_hms () ) then
      exit_string = 'FATAL ERROR: invalid ymd_hms'
      call document ('init_time2', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!
!  Compute decimal time in days, months, years, etc
!  NOTE:  newhour and newsecond are not set initially
!
!-----------------------------------------------------------------------

   call get_tday

   call int_to_char(4,iyear,cyear)
   cday    = cdays     (iday)
   cmonth  = cmonths   (imonth)
   cmonth3 = month3_all(imonth)

   nsteps_this_interval = 0
   nsteps_run = 0

!-----------------------------------------------------------------------
!
!  thour00_begin_this_year, thour00_midmonth_{equal,calendar} and 
!       thour00_endmonth_{equal,calendar} are used in forcing routines 
!       with the 'monthly-equal' or  'monthly-calendar'  option for 
!       forcing_data_freq, where 'monthly-equal' designates 12 equally 
!       spaced months of length 365/12 days and 'monthly-calendar' uses 
!       the non-leapyear calendar.
!       
!  thour00_midmonth_{equal,calendar} and 
!  thour00_endmonth_{equal,calendar} 
!       are relative to the beginning of the year, so vary between 
!       0 and 365*24.
!
!-----------------------------------------------------------------------

   thour00_begin_this_year = thour00 - &
                             (seconds_this_year/seconds_in_hour)

   thour00_midmonth_calendar(1) = 24.0_r8*p5*days_in_month(1)
   thour00_endmonth_calendar(1) = 24.0_r8*days_in_month(1)

   thour00_endmonth_equal(1) = hours_in_year/12.0_r8
   thour00_midmonth_equal(1) = p5*thour00_endmonth_equal(1)

   do nm = 2,12

      thour00_endmonth_calendar(nm) = thour00_endmonth_calendar(nm-1) &
                                    + 24.0_r8*days_in_month(nm)
      thour00_midmonth_calendar(nm) = thour00_endmonth_calendar(nm-1) &
                                    + 24.0_r8*p5*days_in_month(nm)

      thour00_endmonth_equal(nm) = thour00_endmonth_equal(1)*nm
      thour00_midmonth_equal(nm) = thour00_midmonth_equal(nm-1) + &
                                   thour00_endmonth_equal(1)

   enddo

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
!  begin error checking -- after restart file has been read
!
!-----------------------------------------------------------------------

   if (tmix_iopt == tmix_avgfit) then 

     !***  does this run start at the beginning of a coupling interval?

      error_condition = .true.

      do nfit = 1, fit_freq
        if (ihour_start_run == hour_at_interval(nfit)  .and.  &
            iminute_start_run == min_at_interval(nfit) .and.  &
            isecond_start_run == sec_at_interval(nfit)) error_condition = .false.
      enddo

      if (error_condition) then
         write(stdout,*) 'ihour_start_run   = ', ihour_start_run
         write(stdout,*) 'iminute_start_run = ', iminute_start_run
         write(stdout,*) 'isecond_start_run = ', isecond_start_run
         call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
         exit_string = 'FATAL ERROR: model run must start at coupling boundary when using avgfit option'
         call document ('init_time2', exit_string)
         call exit_POP (sigAbort, exit_string, out_unit=stdout)
      endif
   endif
 
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

   if (tmix_iopt == tmix_avgfit) end_run_at_midnight = .true.

!-----------------------------------------------------------------------
!
!  determine iyear, imonth, and iday for the end of this run
!
!-----------------------------------------------------------------------

   stop_iopt = stop_opt_sometime
 
   select case (stop_option)        
 
   case ('never')   !*** coupler or signal catcher stops POP     

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
      write(exit_string, '(a,a)') 'FATAL ERROR: Invalid stop_option: ', stop_option
      call document ('init_time2', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
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
      exit_string = 'FATAL ERROR: Invalid end date'
      call document ('init_time2', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
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

 subroutine time_manager (lcoupled, liceform)

! !DESCRIPTION:
!  This routine updates various time-related variables to their 
!  end-of-step values.  It is called once at the beginning of each 
!  timestep.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (log_kind), intent(in) :: &
      lcoupled,       &! flag for when model is coupled
      liceform         ! flag to determine when ice formation is on

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
   isecond_last      = isecond

   eod_last          = eod
   eom_last          = eom
   eoy_last          = eoy
 
   midnight_last     = midnight

   avg_ts_last       = avg_ts
 
!-----------------------------------------------------------------------
!
!  set logical switches and time_flags to their default values
!
!-----------------------------------------------------------------------

   call reset_switches
   call reset_time_flags 

!-----------------------------------------------------------------------
!
!  increment the timestep counters
!
!-----------------------------------------------------------------------

   nsteps_run   = nsteps_run   + 1
   nsteps_total = nsteps_total + 1

   if (tmix_iopt == tmix_avgfit) then
      nsteps_this_interval = nsteps_this_interval + 1
      if (nsteps_this_interval > nsteps_per_interval) &
         nsteps_this_interval = 1
   endif

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

   if (avg_ts .or. back_to_back) then
      stepsize = p5*dtt
   else
      stepsize =    dtt
   endif

!-----------------------------------------------------------------------
!
!  compute seconds_this_year and seconds_the_day for this timestep
!
!-----------------------------------------------------------------------

   if (nsteps_total == 1 .or.  new_dtt_value) then
      seconds_this_year  = seconds_this_year  + stepsize
      seconds_this_day   = seconds_this_day   + stepsize
   else
      seconds_this_year  = seconds_this_year_next
      seconds_this_day   = seconds_this_day_next
   endif

!-----------------------------------------------------------------------
!
!  determine the size of the next timestep
!
!-----------------------------------------------------------------------

   if ( avg_ts_next .or. back_to_back_next ) then
      stepsize_next = p5*dtt
   else
      stepsize_next =    dtt
   endif

!-----------------------------------------------------------------------
!
!  compute seconds_this_year and seconds_the_day for next timestep
!
!-----------------------------------------------------------------------

   seconds_this_year_next =  seconds_this_year + stepsize_next
   seconds_this_day_next  =  seconds_this_day  + stepsize_next

!-----------------------------------------------------------------------
!
!  adjust seconds counters if necessary
!
!-----------------------------------------------------------------------

   if (nsteps_total == 1 .or.  new_dtt_value) then
      call reduce_seconds (seconds_this_day      , &
                           seconds_this_year     , adjust_year)
      call reduce_seconds (seconds_this_day_next , &
                           seconds_this_year_next, adjust_year_next)
   else
      adjust_year = adjust_year_next
      call reduce_seconds (seconds_this_day_next , &
                           seconds_this_year_next, adjust_year_next)
   endif

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
   if (isecond /= isecond_last) newsecond = .true.

!-----------------------------------------------------------------------
!
!  reset thour00_begin_this_year to be the current time if new year.
!
!-----------------------------------------------------------------------

   if (eoy) thour00_begin_this_year = thour00

!-----------------------------------------------------------------------
!
!  now that new date and time are determined, evaluate all user-defined 
!   time flags
!
!-----------------------------------------------------------------------

   call eval_time_flags

!-----------------------------------------------------------------------
!
!  ice formation or sample qflux this timestep?
!  this section must follow call eval_time_flags
!
!-----------------------------------------------------------------------
      
   if (liceform) then
      if (lcoupled) then

         if (check_time_flag(coupled_ts) ) then
            if (tmix_iopt == tmix_matsuno .or. &
                tmix_iopt == tmix_avgfit       ) then
               ice_ts = .true.
            else
               exit_string = 'FATAL ERROR: Cannot use tmix_avg or tmix_avgbb with lcoupled and liceform'
               call document ('time_manager', exit_string)
               call exit_POP (sigAbort, exit_string, out_unit=stdout)
            endif
          
            if (avg_ts) then
               exit_string = 'FATAL ERROR: Cannot have coupled timestep and be an averaging timestep'
               call document ('time_manager', exit_string)
               call exit_POP (sigAbort, exit_string, out_unit=stdout)
            endif
         endif

         if ( tmix_iopt == tmix_avgfit ) then
            if (mod(nsteps_this_interval+1,nsteps_per_interval) == 0) &
               ice_ts = .true.
         endif

      else  ! .not. lcoupled

         if (tmix_iopt == tmix_matsuno) then
            if (mod(nsteps_total,time_mix_freq) == time_mix_freq-1) &
               ice_ts = .true.
         else if (tmix_iopt == tmix_avgfit  .or. &
                  tmix_iopt == tmix_avg          ) then
            ice_ts = .true.
         else
            exit_string = 'FATAL ERROR: tmix_avgbb option is inconsistent with ice_ts'
            call document ('time_manager', exit_string)
            call exit_POP (sigAbort, exit_string, out_unit=stdout)
         endif

      endif ! coupled

      if (ice_ts) sample_qflux_ts = .true.

   endif  ! liceform

!-----------------------------------------------------------------------
!
!  if using Matsuno with a coupled model, take a Matsuno step
!  after a coupling step to guarantee conservation of integrated fluxes
!
!-----------------------------------------------------------------------

   if (tmix_iopt == tmix_matsuno) then
      if (check_time_flag(coupled_ts,old_value=.true.) .and. nsteps_run > 1 ) then
         matsuno_ts = .true.
         leapfrogts = .false.
      endif
   endif

!-----------------------------------------------------------------------
!
!  stop after this timestep?
!
!-----------------------------------------------------------------------

   if (stop_option == 'nstep' .or. stop_option == 'nsteps') then
      if (nsteps_run == stop_count) call override_time_flag(stop_now,value=.true.)

   else if (stop_iopt /= stop_opt_never .and. eod) then
 
      if (iyear == iyear_end_run .and. imonth == imonth_end_run   &
                                 .and. iday == iday_end_run) then

         call override_time_flag(stop_now,value=.true.)
 
         if (stop_option == 'eoy' .and. .not. eoy) then
            call override_time_flag(stop_now,value=.false.)
         endif
         if (stop_option == 'eom' .and. .not. eom) then
            call override_time_flag(stop_now,value=.false.)
         endif
 
      else if (elapsed_days > elapsed_days_max ) then

         call override_time_flag(stop_now,value=.true.)

      endif

      if (stop_option == 'eoy' .and.  eoy  .and. &
          elapsed_years_this_run == stop_count)  &
         call override_time_flag(stop_now,value=.true.)

      if (stop_option == 'eom' .and.  eom  .and. &
          elapsed_months_this_run == stop_count) &
         call override_time_flag(stop_now,value=.true.)
 
   endif

   new_dtt_value = .false.

!-----------------------------------------------------------------------
!     report ocn model time daily
!-----------------------------------------------------------------------

 if (eod .or. debug_time_management .or. registry_match('info_debug_ge2')) then
  if (my_task == master_task) then
    if (iyear <= 9999) then
      write(stdout,1000) iyear, cmonth3, iday, seconds_this_day
    else
      write(stdout,1001) iyear, cmonth3, iday, seconds_this_day
    endif
    call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
  endif
 endif
1000 format (' (time_manager)', ' ocn date ', i4.4, '-', a3, '-', &
                                  i2.2,', ', 1pe12.6, ' sec')
1001 format (' (time_manager)', ' ocn date ', i5.5, '-', a3, '-', &
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

   newday             = .false.  ! not new day
   newhour            = .false.  ! not new hour
   newsecond          = .false.  ! not new second

   leapfrogts         = .true.   ! a leapfrog timestep
   f_euler_ts         = .false.  ! not a forward Euler timestep
   matsuno_ts         = .false.  ! not an Euler-backward timestep
   avg_ts             = .false.  ! not a time-averaging timestep
   avg_ts_next        = .false.  ! not the step before avg step
   back_to_back       = .false.  ! not a second time-averaging step
   back_to_back_next  = .false.  ! not the step before 2nd avg step
   ice_ts             = .false.  ! not an ice timestep
   sample_qflux_ts    = .false.  ! do not sample qflux 


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
      leapfrogts = .false.
      f_euler_ts = .true.
      newday     = .true.
      first_step = .false.
   endif

!-----------------------------------------------------------------------
!
!  set avg_ts flags    (avg or avgbb options)
!
!-----------------------------------------------------------------------

   if (tmix_iopt == tmix_avg .or. tmix_iopt == tmix_avgbb) then

      if (mod(nsteps_total+1,time_mix_freq) == 0) avg_ts_next  = .true.
      if (mod(nsteps_total  ,time_mix_freq) == 0) avg_ts       = .true.

   endif
 
!-----------------------------------------------------------------------
!
!  set back-to-back flags  (avgbb option)   
!
!-----------------------------------------------------------------------

   if (tmix_iopt == tmix_avgbb) then
      if (avg_ts) back_to_back_next  = .true.
      if (nsteps_total > 1 .and. &
          mod(nsteps_total-1,time_mix_freq) == 0) back_to_back = .true.
   endif

!-----------------------------------------------------------------------
!
!  set avg_ts flags (avgfit option)
!
!-----------------------------------------------------------------------
 
   if (tmix_iopt == tmix_avgfit) then

      if      (nsteps_this_interval == 1) then
         avg_ts_next = .true.
      else if (nsteps_this_interval == 2) then
         avg_ts      = .true.
      else if (mod(nsteps_this_interval+1,time_mix_freq) == 0) then
         avg_ts_next = .true.
      else if (mod(nsteps_this_interval  ,time_mix_freq) == 0) then
         avg_ts      = .true.
      endif
 
      !*** no averaging step in the first step of an interval

      if (nsteps_this_interval == nsteps_per_interval) then
         avg_ts_next = .false.   
      endif
 
   endif

!-----------------------------------------------------------------------
!
!  use Euler backward timestep this timestep?
!
!-----------------------------------------------------------------------

   if (tmix_iopt == tmix_matsuno .and. nsteps_total /= 1) then
      if (mod(nsteps_total,time_mix_freq) == 0 .or. &
          time_mix_freq == 1) then
         matsuno_ts = .true.
         leapfrogts = .false.
      endif
   endif

!-----------------------------------------------------------------------
!EOC

  

 end subroutine set_switches

!***********************************************************************
!BOP
! !IROUTINE: init_time_flag
! !INTERFACE:

 subroutine init_time_flag(flag_name, flag_id, owner, default, freq_opt, freq,  &
                           offset_year, offset_month, offset_day)

! !DESCRIPTION:
!  Creates a user-defined time flag with optional default values
!  and frequency.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      flag_name               ! name for this flag

   character (*), intent(in) :: &
      owner                   ! name of routine initializing this flag

   logical (log_kind), intent(in), optional :: &
      default                 ! default state for this flag

   integer (int_kind), intent(in), optional :: &
      freq_opt,                                &! optional freq option for setting flag
      freq,                                    &! freq in above units  for setting flag
      offset_year,                             &! optional reference/offset year
      offset_month,                            &! optional reference/offset month
      offset_day                                ! optional reference/offset day

! !OUTPUT PARAMETERS:
  
   integer (int_kind), intent(out) ::  &
      flag_id

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   logical (log_kind) :: &
      init_time_flag_error

   integer (int_kind) :: &
      isearch

!-----------------------------------------------------------------------
!
!  search time flags; isearch = 0 ==> flag does not yet exist
!
!-----------------------------------------------------------------------

   call search_time_flags(flag_name,isearch)

!-----------------------------------------------------------------------
!
!  only one owner can initialize a flag
!
!-----------------------------------------------------------------------

   if (isearch > 0 ) then
    if ( .not. time_flags(isearch)%is_reserved) then  
      exit_string = 'FATAL ERROR: Cannot initialize an existing time_flag'
      call document ('init_time_flag', exit_string)
      call document ('init_time_flag', 'flag_name = ', trim(flag_name))
      call document ('init_time_flag', 'owner     = ', trim(owner))
      call document ('init_time_flag', 'time_flag = ', isearch)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
    endif
   endif

!-----------------------------------------------------------------------
!
!  get new flag_id, if flag has not already been reserved
!  use existing id if flag has been reserved
!
!-----------------------------------------------------------------------

   if (isearch == 0) then
      call set_num_time_flags (num_time_flags)
      flag_id = num_time_flags
   else
      flag_id = isearch
   endif


!-----------------------------------------------------------------------
!
!  initialize time flag
!
!-----------------------------------------------------------------------

   time_flags(flag_id)%name            = trim(flag_name)
   time_flags(flag_id)%is_initialized  = .true.
   time_flags(flag_id)%owner           = trim(owner)
   time_flags(flag_id)%has_default     = .false.
   time_flags(flag_id)%has_offset_date = .false.
   time_flags(flag_id)%default         = .false.
   time_flags(flag_id)%freq_opt        = freq_opt_never
   time_flags(flag_id)%freq            = 0
   time_flags(flag_id)%value           = .false.
   time_flags(flag_id)%old_value       = .false.
   time_flags(flag_id)%done            = .false.

!-----------------------------------------------------------------------
!
!  set default value of flag, if requested
!
!  NOTE: If flag previously defined and optional arguments are 
!        present, this will override any previous definition of 
!        optional arguments. user must make sure calls do not 
!        contain optional arguments or else that the last call to 
!        this routine for a specific flag contains desired values.
!
!-----------------------------------------------------------------------

   if (present(default)) then
      time_flags(flag_id)%has_default = .true.
      time_flags(flag_id)%default     = default
      time_flags(flag_id)%value       = default
      time_flags(flag_id)%old_value   = default
   endif

!-----------------------------------------------------------------------
!
!  define optional frequency for setting flag
!
!-----------------------------------------------------------------------

   if (present(freq_opt)) then
      time_flags(flag_id)%freq_opt = freq_opt
      time_flags(flag_id)%freq     = freq
   endif

!-----------------------------------------------------------------------
!
!  define optional offset date information
!
!-----------------------------------------------------------------------

   !*** first, do some error checking. Must have all three offset_* values, or none

   init_time_flag_error = .false.

   if (present(offset_year)) then
      if (.not. present(offset_month) .or. .not. present(offset_day  )) init_time_flag_error = .true.
   else 
      if (present(offset_month) .or. present(offset_day)) init_time_flag_error = .true.
   endif

   if (init_time_flag_error) then
      write(exit_string,'(a,a)')  '  time_file name  = ', trim(time_flags(flag_id)%name)
      call document ('init_time_flag', exit_string)
      exit_string = 'FATAL ERROR: missing time_flag ref values'
      call document ('init_time_flag', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

   !*** now compute elapsed reference years, months, days since 01-01-0000
   if (present(offset_year)) then
      time_flags(flag_id)%has_offset_date = .true.
      time_flags(flag_id)%offset_year     = offset_year
      time_flags(flag_id)%offset_month    = offset_month
      time_flags(flag_id)%offset_day      = offset_day
      time_flags(flag_id)%eyears          = offset_year
      time_flags(flag_id)%emonths         = 12*offset_year + offset_month - 1
      call ymd2eday (offset_year,offset_month,offset_day,time_flags(flag_id)%edays)
   endif

 
   if (debug_time_management) call document_time_flag(flag_id)


!-----------------------------------------------------------------------
!EOC

 end subroutine init_time_flag

!***********************************************************************
!BOP
! !IROUTINE: access_time_flag
! !INTERFACE:

 subroutine access_time_flag(flag_name, flag_id)

! !DESCRIPTION:
!  Allows non-owner routines to access time_flag information. If
!  a routine tries to use the time flag prior to its initialization
!  the time flag is "reserved" and can be initialized later.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      flag_name               ! name for this flag

! !OUTPUT PARAMETERS:
  
   integer (int_kind), intent(out) ::  &
      flag_id

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      isearch

!-----------------------------------------------------------------------
!
!  search time flags; isearch = 0 ==> flag does not yet exist
!
!-----------------------------------------------------------------------

   call search_time_flags(flag_name,isearch)

!-----------------------------------------------------------------------
!
!  if time_flag does not yet exist, reserve the time_flag name and assign id
!  if time_flag does exist, then return its id
!
!-----------------------------------------------------------------------

   if (isearch == 0) then  

      !---------------------
      !  get new flag_id
      !---------------------

      call set_num_time_flags (num_time_flags)

      flag_id = num_time_flags

      !-----------------------
      !  reserve time flag
      !-----------------------

      time_flags(flag_id)%name            = trim(flag_name)
      time_flags(flag_id)%is_reserved     = .true.

   else

      flag_id = isearch

   endif

 
!-----------------------------------------------------------------------
!EOC

 end subroutine access_time_flag

!***********************************************************************
!BOP
! !IROUTINE: set_num_time_flags
! !INTERFACE:

 subroutine set_num_time_flags(num_time_flags)

! !DESCRIPTION:
!  Sets the value of num_time_flags 
!
!
! !REVISION HISTORY:
!  same as module


! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind), intent(inout) :: &
      num_time_flags          ! flag id which also is integer index 

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  increment num_time_flags; cannot exceed maximum number of time flags
!
!-----------------------------------------------------------------------

   if (num_time_flags + 1 <= max_time_flags) then
       num_time_flags = num_time_flags + 1
   else
      exit_string = 'FATAL ERROR: num_time_flags exceeds max_time_flags'
      call document ('set_num_time_flags', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine set_num_time_flags

!***********************************************************************
!BOP
! !IROUTINE: get_time_flag_id
! !INTERFACE:

 function get_time_flag_id(flag_name)

! !DESCRIPTION:
!  Returns flag_id for existing flag. Fatal error if flag does not exist.
!
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      flag_name               ! name for this flag

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      get_time_flag_id        ! flag id which also is integer index 
                              !    into time flag array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (char_len) ::  &
      string               ! dummy string

   integer (int_kind) ::  &
      n,                  &! dummy loop index
      isearch             ! index of flag during search

!-----------------------------------------------------------------------
!
!  search for existing time_flag
!
!-----------------------------------------------------------------------

   call search_time_flags(flag_name, isearch)

!-----------------------------------------------------------------------
!
!  if flag does not exist, or if it has not been initialized, fatal error
!
!-----------------------------------------------------------------------

   if (isearch == 0) then
      exit_string = 'FATAL ERROR: Cannot request id of a nonexistent time_flag; model will abort'
      call document ('get_time_flag_id', 'flag_name = ', trim(flag_name))
      call document ('get_time_flag_id', trim(string) )
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   else if (.not. time_flags(isearch)%is_initialized) then
      exit_string = 'FATAL ERROR: Cannot request id of time_flag that has not been initialized; model will abort'
      call document ('get_time_flag_id', 'flag_name = ', trim(flag_name))
      call document ('get_time_flag_id', 'flag_id   = ', isearch)
      call document ('get_time_flag_id', trim(string) )
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!
!  return array index as integer flag id
!
!-----------------------------------------------------------------------

   get_time_flag_id = isearch

!-----------------------------------------------------------------------
!EOC

 end function get_time_flag_id

!***********************************************************************
!BOP
! !IROUTINE: search_time_flags
! !INTERFACE:

 subroutine search_time_flags(flag_name,isearch)

! !DESCRIPTION:
!  Determines a time_flag id; time_flag_id = 0 ==> does not exist
!
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      flag_name               ! name for this flag

! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      isearch                 ! index of flag from search

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  n  ! dummy loop index

!-----------------------------------------------------------------------
!
!  search to determine if flag already exists 
!     isearch = 0 ==> does not exist
!     isearch > 0 == flag id
!
!-----------------------------------------------------------------------

   isearch = 0
   flag_search: do n=1,num_time_flags
      if (trim(time_flags(n)%name) == trim(flag_name)) then
         isearch = n
         exit flag_search
      endif
   end do flag_search


!-----------------------------------------------------------------------
!EOC

 end subroutine search_time_flags

!***********************************************************************
!BOP
! !IROUTINE: document_time_flags
! !INTERFACE:

 subroutine document_time_flags

! !DESCRIPTION:
!  Documents all time flags
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

  integer (int_kind) ::  &
     n

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' Time Flags '
      write(stdout,blank_fmt)
   endif ! master_task
     
   do n=1,num_time_flags
     call document_time_flag(n)
   enddo

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
   endif ! master_task

!-----------------------------------------------------------------------
!EOC

 end subroutine document_time_flags

!***********************************************************************
!BOP
! !IROUTINE: document_time_flag
! !INTERFACE:

 subroutine document_time_flag(flag_id)

! !DESCRIPTION:
!  Documents time flag
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in), optional :: &
      flag_id                ! index of flag array identifying flag
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

  integer (int_kind) ::   n

  character (char_len) ::  string

   if (my_task == master_task) then
      write(stdout,blank_fmt)

      do n = flag_id,flag_id
      if (time_flags(n)%is_initialized) then
 
        select case (time_flags(n)%freq_opt)
           case (freq_opt_never)
             string = '  (never)'
           case (freq_opt_nyear)
             string = '  (nyear)'
           case (freq_opt_nmonth)
             string = '  (nmonth)'
           case (freq_opt_nday)
             string = '  (nday)'
           case (freq_opt_nhour)
             string = '  (nhour)'
           case (freq_opt_nsecond)
             string = '  (nsecond)'
           case (freq_opt_nstep)
             string = '  (nstep)'
           case (freq_opt_once)
             string = '  (once)'
        end select

        write(stdout,*) 'time_flag id = ', n
        write(stdout,*) '   name               = ', trim(time_flags(n)%name)
        write(stdout,*) '   owner              = ', trim(time_flags(n)%owner)
        write(stdout,*) '   has_default        = ', time_flags(n)%has_default
        write(stdout,*) '   default            = ', time_flags(n)%default
        write(stdout,*) '   freq_opt           = ', time_flags(n)%freq_opt, trim(string)
        write(stdout,*) '   freq               = ', time_flags(n)%freq
        write(stdout,*) '   value              = ', time_flags(n)%value
        write(stdout,*) '   old_value          = ', time_flags(n)%old_value
        write(stdout,*) '   done               = ', time_flags(n)%done
        write(stdout,*) '   has_offset_date    = ', time_flags(n)%has_offset_date 
        if (time_flags(n)%has_offset_date) then
          write(stdout,1100) '   Reference date: ', time_flags(n)%offset_year,'-',  &
                                                    time_flags(n)%offset_month,'-', &
                                                    time_flags(n)%offset_day
          write(stdout,*) '   eyears           = ', time_flags(n)%eyears 
          write(stdout,*) '   emonths          = ', time_flags(n)%emonths 
          write(stdout,*) '   edays            = ', time_flags(n)%edays  
        endif
       else
        write(stdout,*) 'time_flag #', n, ' is NOT initialized'
        write(stdout,*) '   name              = ', trim(time_flags(n)%name)
        write(stdout,*) '   is_reserved       = ', time_flags(n)%is_reserved
       endif ! is_initialized
       write(stdout,*) '   '
      enddo ! n

   endif ! master_task

   call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)

1100 format (1x, a, 1x, i4.4, a, i2.2, a, i2.2)

!-----------------------------------------------------------------------
!EOC

 end subroutine document_time_flag

!***********************************************************************
!BOP
! !IROUTINE: override_time_flag
! !INTERFACE:

 subroutine override_time_flag(flag_id, value, old_value, done)

! !DESCRIPTION:
!  Explicitly overrides the current value of the time flag by
!  setting the time flag to the value specified, as follows:
!     if (present(value), time-flag value is set to value
!     if (present(old_value), time-flag old value is set to old_value
!  NOTE: the override is in effect only until the start of the next timestep
 
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      flag_id                         ! index of flag array identifying flag

   logical (log_kind), intent(in), optional :: &
      value,                                   &! value requested for flag
      old_value,                               &
      done

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  check for proper flag id and then set flag
!
!-----------------------------------------------------------------------

   call error_check_time_flag(flag_id,'override_time_flag')

!-----------------------------------------------------------------------
!
!  Set the time-flag as follows:
!     if (present(value),     time_flags(flag_id)%value = value 
!     if (present(old_value), time_flags(flag_id)%old_value = old_value 
!
!-----------------------------------------------------------------------

   if (present(value)) then
     time_flags(flag_id)%value = value
   else if (present(old_value)) then
     time_flags(flag_id)%old_value = old_value
   else if (present(done)) then
     time_flags(flag_id)%done = done
   else
      exit_string = 'FATAL ERROR: must specify value or old_value'
      call document ('override_time_flag', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif
      

!-----------------------------------------------------------------------
!EOC

 end subroutine override_time_flag

!***********************************************************************
!BOP
! !IROUTINE: reset_time_flags
! !INTERFACE:

 subroutine reset_time_flags

! !DESCRIPTION:
!  Reset all time flags to their default values
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

   do n=1,num_time_flags
      call reset_time_flag(n)
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine reset_time_flags


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

   call error_check_time_flag(flag_id, 'reset_time_flag')

   if (time_flags(flag_id)%has_default) then
      time_flags(flag_id)%value = time_flags(flag_id)%default
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine reset_time_flag


!***********************************************************************
!BOP
! !IROUTINE: check_time_flag
! !INTERFACE:

 function check_time_flag(flag_id,old_value,done)

! !DESCRIPTION:
!  Returns the current value of time flag given by flag\_id.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      flag_id                ! index of flag array identifying flag

   logical (log_kind), intent(in), optional :: &
      old_value,                  &! check old value of time flag
      done                         ! check once option -- has it been done?

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      check_time_flag      ! resulting value of time flag (logical)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      check_old_value,   &! check old value, not present value
      check_done          ! check done, not present value

!-----------------------------------------------------------------------
!
!  check for proper flag id and then set flag
!
!-----------------------------------------------------------------------

   call error_check_time_flag(flag_id, 'check_time_flag')

   if (present(old_value)) then
     check_old_value = old_value
   else
     check_old_value = .false.
   endif
      
   if (present(done)) then
     check_done = done
   else
     check_done = .false.
   endif
      
   if (check_old_value) then
      check_time_flag = time_flags(flag_id)%old_value
   else if (check_done) then
      check_time_flag = time_flags(flag_id)%done
   else
      check_time_flag = time_flags(flag_id)%value
   endif

!-----------------------------------------------------------------------
!EOC

 end function check_time_flag

!***********************************************************************
!BOP
! !IROUTINE: check_time_flag_int
! !INTERFACE:

 function check_time_flag_int(flag_id, freq, freq_opt)

! !DESCRIPTION:
!  Returns the current integer value of time flag given by flag\_id.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      flag_id                ! index of flag array identifying flag

   logical (log_kind), intent(in), optional :: &
      freq,                 &! check flag frequency
      freq_opt               ! check flag frequency option


! !OUTPUT PARAMETERS:

   integer (int_kind) :: &
      check_time_flag_int        ! current freqeuncy option of time flag

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   logical (log_kind) :: &
      check_freq,        &! check flag frequency
      check_freq_opt      ! check flag frequency option
 
!-----------------------------------------------------------------------
!
!  check for proper flag id and then return flag value
!
!-----------------------------------------------------------------------

   call error_check_time_flag(flag_id, 'check_time_flag_int')

   if (present(freq)) then
     check_freq = freq
   else
     check_freq = .false.
   endif

   if (present(freq_opt)) then
     check_freq_opt= freq_opt
   else
     check_freq_opt = .false.
   endif

!-----------------------------------------------------------------------
!
!  check one option per call
!
!-----------------------------------------------------------------------

   if (check_freq .and. check_freq_opt) then
      exit_string = 'FATAL ERROR: check one option at a time'
      call document ('check_time_flag_int', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   else if (.not. (check_freq .or. check_freq_opt)) then
      exit_string = 'FATAL ERROR: must check at leaset one option'
      call document ('check_time_flag_int', exit_string)
      call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif
      
   if (check_freq) then
      check_time_flag_int = time_flags(flag_id)%freq
   endif
   
   if (check_freq_opt) then
      check_time_flag_int = time_flags(flag_id)%freq_opt
   endif

!-----------------------------------------------------------------------
!EOC

 end function check_time_flag_int

!***********************************************************************
!BOP
! !IROUTINE: eval_time_flags
! !INTERFACE:

 subroutine eval_time_flags

! !DESCRIPTION:
!  Evaluates all time flags based upon flag frequencies, via the function time_to_do
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
!  if it is time, save value from old time and set new value
!
!-----------------------------------------------------------------------

   do n=1,num_time_flags
     call eval_time_flag(n)
   end do

!-----------------------------------------------------------------------

 end subroutine eval_time_flags

!***********************************************************************
!BOP
! !IROUTINE: eval_time_flag
! !INTERFACE:

 subroutine eval_time_flag(flag_id)

! !DESCRIPTION:
!  Evaluates the requested time flag based upon flag frequency, via the function time_to_do
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
!  check valid time flag
!
!-----------------------------------------------------------------------

  call error_check_time_flag(flag_id, 'eval_time_flag')

!-----------------------------------------------------------------------
!
!  if it is time, save value from old time and set new value
!
!-----------------------------------------------------------------------

  if (time_flags(flag_id)%freq_opt /= freq_opt_never) then

      time_flags(flag_id)%old_value = time_flags(flag_id)%value
      time_flags(flag_id)%value     = time_to_do(flag_id)

  endif

!-----------------------------------------------------------------------
!EOC

 end subroutine eval_time_flag

!***********************************************************************
!BOP
! !IROUTINE: error_check_time_flag
! !INTERFACE:

 subroutine error_check_time_flag(flag_id, calling_routine)

! !DESCRIPTION:
!  Explicitly overrides the current value of the time flag by
!  setting the time flag to the value specified, as follows:
!     if (present(value), time-flag value is set to value
!     if (present(old_value), time-flag old value is set to old_value
!  NOTE: the override is in effect only until the start of the next timestep
 
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      flag_id                         ! index of flag array identifying flag

   character (*) :: calling_routine

!EOP
!BOC

   character (char_len) ::  string

!-----------------------------------------------------------------------
!
!  check for proper flag id 
!
!-----------------------------------------------------------------------


   if (flag_id < 1 .or. flag_id > num_time_flags) then
     write(exit_string,'(a,a)') 'FATAL ERROR: ',trim(calling_routine) // ' ' // trim(time_flags(flag_id)%name)//': invalid flag_id'
     call document ('error_check_time_flag', exit_string)
     call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

!-----------------------------------------------------------------------
!
!  check for initialized time_flag
!
!-----------------------------------------------------------------------

   if (.not. time_flags(flag_id)%is_initialized) then
     write(exit_string,'(a)') 'FATAL ERROR: ' // trim(time_flags(flag_id)%name)//': is not initialized'
     call document ('error_check_time_flag', exit_string)
     call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif
      

!-----------------------------------------------------------------------
!EOC

 end subroutine error_check_time_flag

!***********************************************************************
!BOP
! !IROUTINE: time_to_do
! !INTERFACE:

 function time_to_do (flag_id)

! !DESCRIPTION:
!  Determines whether it is time to take a particular action based on 
!  input frequency options.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      flag_id                         ! flag id

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      time_to_do         ! true if current timestep matches input
                         ! frequency conditions

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      mod_test,          &
      freq_opt,          &
      freq

   logical (log_kind) :: &
      has_offset_date

!-----------------------------------------------------------------------
!
!  check for proper flag id and proper time sequencing, then set time_to_do
!
!-----------------------------------------------------------------------

   call error_check_time_flag(flag_id, 'time_to_do')

   if (.not. registry_match('init_time2')) then
     exit_string = 'FATAL ERROR: init_time2 has not yet been called'
     call document ('time_to_do', exit_string)
     call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

   time_to_do = .false.
 
   freq_opt        = time_flags(flag_id)%freq_opt
   freq            = time_flags(flag_id)%freq
   has_offset_date = time_flags(flag_id)%has_offset_date

   select case (time_flags(flag_id)%freq_opt)

   case (freq_opt_nyear)
   !-----------------------------------------------------------
   ! elapsed whole years between current time and {ref date or start date}
   !-----------------------------------------------------------
      if (has_offset_date) then
        mod_test = elapsed_years - time_flags(flag_id)%eyears 
      else
        mod_test = elapsed_years - elapsed_years0              
      endif
      if (eoy .and. mod(mod_test,freq) == 0) time_to_do = .true.

   case (freq_opt_nmonth)
   !-----------------------------------------------------------
   ! elapsed whole months between current time and {ref date or start date}
   !-----------------------------------------------------------
      if (has_offset_date) then
        mod_test = elapsed_months - time_flags(flag_id)%emonths
      else
        mod_test = elapsed_months-elapsed_months0              
      endif
      if (eom .and. mod(mod_test,freq) == 0) time_to_do = .true.

   case (freq_opt_nday)
   !-----------------------------------------------------------
   ! elapsed whole days between current time and {ref date or start date}
   !-----------------------------------------------------------
      if (has_offset_date) then
        mod_test = elapsed_days - time_flags(flag_id)%edays
      else
        mod_test = elapsed_days - elapsed_days0      
      endif
      if (eod) then
         if (midnight) then
            if (mod(mod_test  ,freq) == 0)  time_to_do = .true.
         else
            if (mod(mod_test+1,freq) == 0)  time_to_do = .true.
         endif
      endif

   case (freq_opt_nhour)
   !-----------------------------------------------------------
   ! convert to hours and test nearest integer
   !-----------------------------------------------------------
      if (has_offset_date) then
        mod_test =  elapsed_days - time_flags(flag_id)%edays*24 + ihour
      else
        mod_test =  elapsed_days*24 + ihour
      endif
      if (newhour .and. mod(mod_test  ,freq) == 0)  time_to_do = .true.

   case (freq_opt_nsecond)
   !-----------------------------------------------------------
   ! convert to seconds and test nearest integer
   !-----------------------------------------------------------
      if (has_offset_date) then
        mod_test =       (elapsed_days - time_flags(flag_id)%edays)*seconds_in_day + seconds_this_day
      else
        mod_test =       elapsed_days*seconds_in_day + seconds_this_day
      endif
      if (mod(mod_test ,freq) == 0)  time_to_do = .true.

   case (freq_opt_nstep)
      if (mod(nsteps_total,freq) == 0) time_to_do = .true.

   case (freq_opt_once)
   !-----------------------------------------------------------
   ! 
   !-----------------------------------------------------------
      if (.not. time_flags(flag_id)%done) time_to_do = .true.

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
     exit_string = 'FATAL ERROR: unknown start option '
     call document ('time_to_start', exit_string)
     call exit_POP (sigAbort, exit_string, out_unit=stdout)
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
      day_inc             ! change in the number of days elapsed
                          !   between timesteps

   logical (log_kind), save ::                 &
      increment_elapsed_months_next = .false., &
      increment_elapsed_years_next  = .false.

!-----------------------------------------------------------------------
!
!  determine iyear
!
!-----------------------------------------------------------------------

   if (adjust_year) then
      iyear = iyear + 1
      days_in_prior_year = days_in_year
      if (allow_leapyear) call leap_adjust
      adjust_year_next = .false.
   endif

!-----------------------------------------------------------------------
!
!  determine iday_of_year, imonth, iday, ihour, iminute, isecond
!                                        rhour, rminute, rsecond
!
!-----------------------------------------------------------------------

   if (nsteps_run == 1) then
      call ymd_hms( seconds_this_year, seconds_this_day, &
                    iday_of_year,                        &
                    imonth, iday   , iday_last,          &
                    ihour , iminute, isecond,            &
                    rhour , rminute, rsecond,            &
                    midnight       , adjust_year)
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
                imonth_next  , iday_next   , iday,         &
                ihour_next   , iminute_next, isecond_next, &
                rhour_next   , rminute_next, rsecond_next, &
                midnight_next, adjust_year_next)

!-----------------------------------------------------------------------
!
!  end of day?
!
!-----------------------------------------------------------------------

   if (iday_of_year_next /= iday_of_year) then
      if (.not. midnight_next)                     eod = .true.
      if (stepsize_next + dt_tol > seconds_in_day) eod = .true.
   endif

   if (midnight) eod = .true.

!-----------------------------------------------------------------------
!
!  newday?   (a timestep can be both eod and newday for dt > 24hrs)
!
!-----------------------------------------------------------------------
 
   if (iday_of_year > iday_of_year_last .and. .not. midnight) &
      newday = .true.
   if (eod_last ) newday = .true.

!-----------------------------------------------------------------------
!
!  end of month?
!
!-----------------------------------------------------------------------

   if (eod) then
 
      if (eom_next) then
         eom      = .true.
         eom_next = .false.
      else
 
         if (imonth_next > imonth_last  .or. &
             imonth_next == 1 .and. imonth_last == 12) then
 
            if (iday <= days_in_month(imonth_last) .and. &
                midnight_next .and. iday_next == 1 ) then
               eom      = .false.
               eom_next = .true.

            elseif (iday <= days_in_month(imonth_last)  .and. &
                                      iday_next >= 1 )  then
               eom      = .true.
               eom_next = .false.

            elseif (midnight .and. midnight_next .and. &
                                   midnight_last) then
               eom      = .false.
               eom_next = .true.
            else
               eom      = .true.
               eom_next = .false.
            endif

         endif
      endif
 
   endif ! eod

!-----------------------------------------------------------------------
!
!  elapsed months (integer) 
!
!-----------------------------------------------------------------------

   if ((eom .and. midnight) .or. increment_elapsed_months_next) then
      elapsed_months           = elapsed_months           + 1
      elapsed_months_this_run  = elapsed_months_this_run  + 1
      elapsed_months_init_date = elapsed_months_init_date + 1
      if (increment_elapsed_months_next) & 
          increment_elapsed_months_next = .false.
   else if (eom) then
      increment_elapsed_months_next = .true.
   else 
      increment_elapsed_months_next = .false.
   endif
 
   if (eom_last) eom = .false.

!-----------------------------------------------------------------------
!
!  end of year?
!
!-----------------------------------------------------------------------

   if (eom .and. imonth_next == 1 .and. imonth_last == 12) eoy = .true.       

!-----------------------------------------------------------------------
!
!  adjust elapsed years and elapsed days in the year (integer)
!
!-----------------------------------------------------------------------

   if ((eoy .and. midnight) .or. increment_elapsed_years_next) then

      elapsed_years           = elapsed_years           + 1
      elapsed_years_this_run  = elapsed_years_this_run  + 1
      elapsed_years_init_date = elapsed_years_init_date + 1

      call ymd2eday (iyear , 1, 1, elapsed_days_jan1)
      elapsed_days_this_year  = elapsed_days - elapsed_days_jan1

      if (increment_elapsed_years_next) & 
          increment_elapsed_years_next = .false.

   else if (eoy) then
      increment_elapsed_years_next = .true.
   else 
      increment_elapsed_years_next = .false.
   endif
 
   if (eoy_last) eoy = .false.

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

   if (iday_of_year >= iday_of_year_last) then
      day_inc = iday_of_year - iday_of_year_last
   else
      day_inc = iday_of_year - iday_of_year_last + days_in_prior_year
   endif

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
     exit_string = 'FATAL ERROR: invalid ymd_hms'
     call document ('model_date', exit_string)
     call exit_POP (sigAbort, exit_string, out_unit=stdout)
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

   tmonth  = real(elapsed_months_init_date,kind=r8) + &
             (real(iday,kind=r8)-c1+frac_day)/days_in_month(imonth)

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
                    imonth_loc  , iday_loc   , iday_compare,      &
                    ihour_loc   , iminute_loc, isecond_loc,       &
                    rhour_loc   , rminute_loc, rsecond_loc,       &
                    midnight_loc, adjust_year_loc)

! !DESCRIPTION:
!  Computes integer values iday\_of\_year, iyear, imonth, iday, ihour, 
!  iminute, isecond.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer  (int_kind), intent(in) :: &
      iday_compare          ! day to compare to check day change

! !INPUT/OUTPUT PARAMETERS:

   logical (log_kind), intent(inout) :: &
      adjust_year_loc           ! year adjustment flag
 
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
!  if midnight, increment iday 
!
!-----------------------------------------------------------------------
 
   if (iday_loc == iday_compare .and. midnight_loc) &
      iday_loc =  iday_loc + 1

!-----------------------------------------------------------------------
!
!  if necessary, adjust month value and year-adjustment indicator
!
!-----------------------------------------------------------------------
 
   if (iday_loc > days_in_month(imonth_loc)) then
      iday_loc = iday_loc - days_in_month(imonth_loc)
      imonth_loc = imonth_loc + 1
      if (imonth_loc == 13 ) then
         imonth_loc = 1
         adjust_year_loc = .true.
      endif
   endif

   iday_of_year_loc = days_in_prior_months(imonth_loc) + iday_loc

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
                            seconds_this_year_loc, adjust_year_loc)

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

   logical (log_kind), intent(out) :: &
      adjust_year_loc         ! flag to signal year adjustment

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
    
   if (seconds_this_year_loc >= seconds_in_year - stepsize .and.     &
       (seconds_this_year_loc >= seconds_in_year .or.                &
        is_near(seconds_this_year_loc,seconds_in_year,dt_tol_year))) &
                                                                then
 
      seconds_this_year_loc = seconds_this_year_loc - seconds_in_year

      if (is_near(seconds_this_year_loc, c0, dt_tol)) then
         seconds_this_year_loc = c0
         seconds_this_day_loc  = c0
      endif
 
      adjust_year_loc  = .true.
   else
      adjust_year_loc  = .false.
   endif

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

   if (.not. valid_date(date)) then
     exit_string = 'FATAL ERROR: invalid date'
     call document ('date2ymd', exit_string)
     call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

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

      nnorm = eday/days_in_norm_year + 1
      nleap = 0
      days  = eday -nnorm*days_in_norm_year
      tdays_in_prior_months = days_in_prior_months

!-----------------------------------------------------------------------
!
!  Compute number of elapsed leap years and "normal" years, and
!  the number of days elapsed in the most recent year
!
!-----------------------------------------------------------------------

   else

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

   if (days <= 0 .or. days > days_in_leap_year ) then
     write (exit_string,'(a,i6)')'FATAL ERROR: days undetermined, days = ', days
     call document ('eday2ymd', exit_string)
     call exit_POP (sigAbort, exit_string, out_unit=stdout)
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

   if (.not. valid_date(date)) then
     write (exit_string,'(a,i6)')'FATAL ERROR: invalid date'
     call document ('date2eday', exit_string)
     call exit_POP (sigAbort, exit_string, out_unit=stdout)
   endif

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
        exit_string = 'FATAL ERROR: selecting order'
        call document ('time_stamp', exit_string)
        call exit_POP (sigAbort, exit_string, out_unit=stdout)
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
         if (mdy) then
           exit_string = 'FATAL ERROR: time_stamp not supported with option=last & mdy order'
           call document ('time_stamp', exit_string)
           call exit_POP (sigAbort, exit_string, out_unit=stdout)
         endif

         if (date_separator /= ' ') then
            write (date_string,ymd_date_fmt1) iyear_last ,date_separator, &
                                          imonth_last,date_separator, &
                                          iday_last
         else
            write (date_string,ymd_date_fmt2) iyear, imonth, iday
         endif
      endif

      if (present(time_string)) then
         exit_string = 'FATAL ERROR: time_stamp not supported with option=last'
         call document ('time_stamp', exit_string)
         call exit_POP (sigAbort, exit_string, out_unit=stdout)
         
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
         call exit_POP(sigAbort, &
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
         exit_string = 'FATAL ERROR: time_stamp not supported with option=last'
         call document ('time_stamp', exit_string)
         call exit_POP (sigAbort, exit_string, out_unit=stdout)
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
! !IROUTINE: ccsm_date_stamp
! !INTERFACE:
 subroutine ccsm_date_stamp (date_string, ymds)

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

  character (4) :: ccsm_cyear
  character (2) :: ccsm_cmonth
  character (2) :: ccsm_cday  
  character (5) :: ccsm_csecond
 
  integer (kind=int_kind) ::  &
     iyear_stamp             ,&
     imonth_stamp            ,&
     iday_stamp              ,&
     itotal_second
 
   date_string = char_blank

!---------------------------------------------------------------------
!     set ixxxx_stamp variables to conform to the ccsm standard
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
       call int_to_char (4,iyear        , ccsm_cyear  )
       call int_to_char (2,imonth       , ccsm_cmonth )
       call int_to_char (2,iday         , ccsm_cday   )
       call int_to_char (5,itotal_second, ccsm_csecond)
       write (date_string,1000) ccsm_cyear, ccsm_cmonth, ccsm_cday, &
                                ccsm_csecond

     case ('ymd')
        call int_to_char (4,iyear_stamp  , ccsm_cyear )
        call int_to_char (2,imonth_stamp , ccsm_cmonth)
        call int_to_char (2,iday_stamp   , ccsm_cday  )
        write (date_string,1000) ccsm_cyear, ccsm_cmonth, ccsm_cday

     case ('ym')
        call int_to_char (4,iyear_stamp  , ccsm_cyear )
        call int_to_char (2,imonth_stamp , ccsm_cmonth)
        write (date_string,1000) ccsm_cyear, ccsm_cmonth

     case ('y')
        call int_to_char (4,iyear_stamp  , ccsm_cyear)
        write (date_string,1000) ccsm_cyear
        
     case default 
        exit_string = 'FATAL ERROR'
        call document ('ccsm_date_stamp', exit_string)
        call exit_POP (sigAbort, exit_string, out_unit=stdout)
 
   end select
 

 1000 format (a4,:,'-',a2:,'-',a2,:,'-',a5)
  
!EOC

!-----------------------------------------------------------------------

  end subroutine ccsm_date_stamp

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

      exit_string = char_blank

      if (.not. valid_year) &
         write(exit_string,err_fmt) &
              'FATAL ERROR: Invalid date (iyear must be > 0 ): iyear = ', iyear
 
      if (.not. valid_month) &
         write(exit_string,err_fmt) &
              'FATAL ERROR: Invalid date ( imonth must be in [1,12] ): imonth = ', &
                                                          imonth
 
      if (.not. valid_day_b) &
         write(exit_string,err_fmt) &
              'FATAL ERROR: Invalid date (iday must be greater than 1): iday = ',iday
 
      if (.not. valid_day_e) &
         write(exit_string,err_fmt) &
              'FATAL ERROR: Invalid date (iday must be less than days_in_month):'/&
              &/' iday = ',iday 
 
      if (.not. valid_hour) &
         write(exit_string,err_fmt) &
              'FATAL ERROR: Invalid date (ihour must be in [0,23] ): ihour = ', ihour
 
      if (.not. valid_minute) &
         write(exit_string,err_fmt) &
              'FATAL ERROR: Invalid date (iminute must be in [0,59] ): iminute = ', &
                                                          iminute
 
      if (.not. valid_second) &
         write(exit_string,err_fmt) &
              'FATAL ERROR: Invalid date (isecond must be in [0,59] ): isecond = ', &
                                                          isecond
 
      if (.not. valid_eday_run) &
         write(exit_string,err_fmt) &
              'FATAL ERROR: Invalid date (elapsed_days_init_date must be > 0 ) ', &
                             elapsed_days_init_date
 
      if (.not. valid_eday_year) &
         write(exit_string,err_fmt) &
              'FATAL ERROR: Invalid date (elapsed_days_this_year must be > 0) ', &
                             elapsed_days_this_year
 
      if (.not. valid_feb_day) then
         write(exit_string,*) &
              'FATAL ERROR: initial date contains leap day '/&
              &/' but no leap years are allowed.', iday
 
         call document ('valid_ymd_hms', exit_string)
         call exit_POP (sigAbort, exit_string, out_unit=stdout)
      endif
 
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
 
   character (3) ::   &
      mix_step

   character (4) ::   &
      cyear0,         &!
      cyear_end_run    !

   character (char_len) :: &
      mix_steps
 
   character (*), parameter :: &! output formats
      out_fmt1 = "('       date(month-day-year):',2x,2(a2,'-'),a4)", &
      out_fmt2 = "('                    ',a7,2x,i10)",               &
      out_fmt3 = "('This run will terminate ',/,a)",                 &
      out_fmt4 = "(a, :i7, a,a)",                                    &
      out_fmt5 = "('                   ',a8,2x,f16.3)",              &
      out_fmt6 = "('Averaging time steps every ',i6,' steps',a)",    &
      out_fmt7 =                                                     &
"('Averaging time steps at the 2nd and ',i5,a2,' step of every day or coupled interval ')",&
      out_fmt8 = "('Surface ',a10,' time step = ',1pe12.6, ' seconds')",&
      out_fmt9 = "('There are ', i6, a6,' steps each day')"

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
 
      if (tmix_iopt == tmix_avgfit) then
         write(stdout,out_fmt9) fullsteps_per_day,' full '
         write(stdout,out_fmt9) halfsteps_per_day,' half '
         write(stdout,out_fmt9) nsteps_per_day,   ' total'
      endif
 
      write (stdout,blank_fmt)
      write (stdout,'(a16,i6)') 'time_mix_freq = ', time_mix_freq
 
      write (stdout,'(a19)') 'Time mixing option:'
      select case (tmix_iopt)
 
      case (tmix_avg)
         write (stdout,'(a23)') '  avg -- time averaging'
 
      case (tmix_avgbb)
         write (stdout,'(a59)') '  avgbb -- time averaging'/&
                               &/' with back-to-back averaging steps'
 
      case (tmix_avgfit)
         write (stdout,'(a26)') '  avgfit -- time averaging'
         write (stdout,'(a71)') '  with timestep chosen to fit'/&
                               &/' exactly into one day or coupling'/&
                               &/' interval'
 
      case (tmix_matsuno)
         write (stdout,'(a25,i6,a6)') &
            'Matsuno time steps every ',time_mix_freq,' steps'
      end select

      write (stdout,blank_fmt)
      select case (tmix_iopt)
      case (tmix_avg)
         write (stdout,out_fmt6) time_mix_freq, ' '
      case (tmix_avgbb)
         write (stdout,out_fmt6) time_mix_freq, &
               ' with back-to-back averaging steps'
      case (tmix_avgfit)
         if (time_mix_freq == 1 .or.           &
             (mod(time_mix_freq,10) == 1 .and. &
              time_mix_freq /= 11)) then
            write (stdout,out_fmt7) time_mix_freq, 'st'
         elseif (time_mix_freq == 2 .or.           &
                 (mod(time_mix_freq,10) == 2 .and. &
                  time_mix_freq /= 12)) then
            write (stdout,out_fmt7) time_mix_freq, 'nd'
         elseif (time_mix_freq == 3 .or.           &
                 (mod(time_mix_freq,10) == 3 .and. &
                  time_mix_freq /= 13)) then
            write (stdout,out_fmt7) time_mix_freq, 'rd'
         else   
            write (stdout,out_fmt7) time_mix_freq, 'th'
         endif
      case (tmix_matsuno)
         write (stdout,'(a25,i6,a6)') &
            'Matsuno time steps every ',time_mix_freq,' steps'
      end select

      write (stdout,blank_fmt)
      write (stdout,out_fmt8) 'tracer    ',dtt
      write (stdout,out_fmt8) 'momentum  ',dtu
      write (stdout,out_fmt8) 'barotropic',dtp

      write (stdout,blank_fmt)
      if (laccel) then
         write (stdout,'(a28)') 'Tracer acceleration  enabled'
         write (stdout,'(a22)') '   k      accel factor'
         write (stdout,'(a22)') '  ---     ------------'
         do k=1,km
            write (stdout,'(2x,i3,7x,f8.3)') k, dttxcel(k)
         end do
      else
         write (stdout,'(a28)') 'Tracer acceleration disabled'
      endif

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

      write (stdout,blank_fmt)
      if (impcor) then
         write (stdout,'(a50)') &
            'Implicit treatment of Coriolis terms (impcor):  ON'
      else
         write (stdout,'(a50)') &
            'Implicit treatment of Coriolis terms (impcor): OFF'
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

!BOP
! !IROUTINE: ccsm_char_date_and_time
! !INTERFACE:

 subroutine ccsm_char_date_and_time

! !DESCRIPTION:
!  This routine converts integer date & time information into character
!  strings.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  set character date and time information
!
!-----------------------------------------------------------------------

   call int_to_char (4,iyear   , cyear  )
   call int_to_char (2,imonth  , cmonth )
   call int_to_char (2,iday    , cday   )
   call int_to_char (2,ihour   , chour  )
   call int_to_char (2,iminute , cminute)
   call int_to_char (2,isecond , csecond)


!-----------------------------------------------------------------------
!EOC

 end subroutine ccsm_char_date_and_time
 
!***********************************************************************

!***********************************************************************

 end module time_management

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
