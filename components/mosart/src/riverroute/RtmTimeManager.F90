module RtmTimeManager

   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_sys_mod , only: shr_sys_abort
   use RtmSpmd     , only: masterproc, iam, mpicom_rof, MPI_INTEGER, MPI_CHARACTER
   use RtmVar      , only: isecspday, iulog, nsrest, nsrContinue
   use RtmIO
   use ESMF


   implicit none
   private

! Public methods

   public ::&
      timemgr_setup,            &! setup startup values
      timemgr_init,             &! time manager initialization
      timemgr_restart,          &! read/write time manager restart info and restart time manager
      advance_timestep,         &! increment timestep number
      get_clock,                &! get the clock from the time-manager
      get_step_size,            &! return step size in seconds
      get_nstep,                &! return timestep number
      get_curr_date,            &! return date components at end of current timestep
      get_prev_date,            &! return date components at beginning of current timestep
      get_start_date,           &! return components of the start date
      get_ref_date,             &! return components of the reference date
      get_curr_time,            &! return components of elapsed time since reference date at end of current timestep
      get_prev_time,            &! return components of elapsed time since reference date at beg of current timestep
      get_calendar,             &! return calendar
      is_first_step,            &! return true on first step of initial run
      is_first_restart_step,    &! return true on first step of restart or branch run
      is_end_curr_day,          &! return true on last timestep in current day
      is_end_curr_month,        &! return true on last timestep in current month
      is_last_step,             &! return true on last timestep
      is_new_year,              &! return true on new year
      is_new_month,             &! return true on new month
      is_new_day,               &! return true on new day
      is_restart                 ! return true if this is a restart run

! Private module data

! Private data for input

   character(len=ESMF_MAXSTR), save :: calendar   = 'NO_LEAP'     ! Calendar to use in date calculations ('NO_LEAP' or 'GREGORIAN')
   integer,  parameter :: uninit_int = -999999999
   real(r8), parameter :: uninit_r8  = -999999999.0

! Input
   integer, save ::&
      dtime          = uninit_int   ! timestep in seconds

! Input from CESM driver
   integer, save ::&
      nelapse       = uninit_int,  &! number of timesteps (or days if negative) to extend a run
      start_ymd     = uninit_int,  &! starting date for run in yearmmdd format
      start_tod     = 0,           &! starting time of day for run in seconds
      stop_ymd      = uninit_int,  &! stopping date for run in yearmmdd format
      stop_tod      = 0,           &! stopping time of day for run in seconds
      ref_ymd       = uninit_int,  &! reference date for time coordinate in yearmmdd format
      ref_tod       = 0             ! reference time of day for time coordinate in seconds
   type(ESMF_Calendar), save  :: &
        tm_cal       ! calendar
   type(ESMF_Clock), save :: &
        tm_clock     ! model clock   
   integer, save ::&                ! Data required to restart time manager:
      rst_nstep     = uninit_int,  &! current step number
      rst_step_days = uninit_int,  &! days component of timestep size
      rst_step_sec  = uninit_int,  &! timestep size seconds
      rst_start_ymd = uninit_int,  &! start date
      rst_start_tod = uninit_int,  &! start time of day
      rst_ref_ymd   = uninit_int,  &! reference date
      rst_ref_tod   = uninit_int,  &! reference time of day
      rst_curr_ymd  = uninit_int,  &! current date
      rst_curr_tod  = uninit_int    ! current time of day
   character(len=ESMF_MAXSTR), save :: &
   rst_calendar    ! Calendar

   logical, save :: tm_first_restart_step = .false.    ! true for first step of a restart or branch run
   integer, save :: cal_type              = uninit_int ! calendar type
   logical, save :: timemgr_set           = .false.    ! true when timemgr initialized
   logical, save :: new_year              = .false.    ! true for new year
   logical, save :: new_month             = .false.    ! true for new month
   logical, save :: new_day               = .false.    ! true for new day

! Private module methods
   private :: timemgr_spmdbcast
   private :: init_calendar
   private :: init_clock
   private :: timemgr_print
   private :: TimeGetymd

contains

!=========================================================================================

subroutine timemgr_setup( calendar_in, start_ymd_in, start_tod_in, ref_ymd_in, &
                          ref_tod_in,  stop_ymd_in,  stop_tod_in,  nelapse_in)

  ! set time manager startup values
  character(len=*), optional, intent(IN) :: calendar_in       ! Calendar type
  integer         , optional, intent(IN) :: nelapse_in        ! Number of step (or days) to advance
  integer         , optional, intent(IN) :: start_ymd_in      ! Start date       (YYYYMMDD)
  integer         , optional, intent(IN) :: start_tod_in      ! Start time of day (sec)
  integer         , optional, intent(IN) :: ref_ymd_in        ! Reference date   (YYYYMMDD)
  integer         , optional, intent(IN) :: ref_tod_in        ! Reference time of day (sec)
  integer         , optional, intent(IN) :: stop_ymd_in       ! Stop date        (YYYYMMDD)
  integer         , optional, intent(IN) :: stop_tod_in       ! Stop time of day (sec)
  character(len=*), parameter :: sub = 'rtm::set_timemgr_init'

  ! timemgr_set is called in timemgr_init and timemgr_restart
  if ( timemgr_set ) then
     call shr_sys_abort( sub//":: timemgr_init or timemgr_restart already called" )
  end if
  if (present(calendar_in) ) calendar  = trim(calendar_in)
  if (present(start_ymd_in)) start_ymd = start_ymd_in
  if (present(start_tod_in)) start_tod = start_tod_in
  if (present(ref_ymd_in)  ) ref_ymd   = ref_ymd_in
  if (present(ref_tod_in)  ) ref_tod   = ref_tod_in
  if (present(stop_ymd_in) ) stop_ymd  = stop_ymd_in
  if (present(stop_tod_in) ) stop_tod  = stop_tod_in
  if (present(nelapse_in)  ) nelapse   = nelapse_in

end subroutine timemgr_setup

!=========================================================================================

subroutine timemgr_init( dtime_in )

  ! Initialize the ESMF time manager from the sync clock
  !
  integer, intent(in) :: dtime_in         ! Time-step (sec)
  !
  integer         :: rc                    ! return code
  integer         :: yr, mon, day, tod     ! Year, month, day, and second as integers
  type(ESMF_Time) :: start_date            ! start date for run
  type(ESMF_Time) :: stop_date             ! stop date for run
  type(ESMF_Time) :: curr_date             ! temporary date used in logic
  type(ESMF_Time) :: ref_date              ! reference date for time coordinate
  type(ESMF_Time) :: current               ! current date (from clock)
  type(ESMF_TimeInterval) :: day_step_size ! day step size
  type(ESMF_TimeInterval) :: step_size     ! timestep size
  logical         :: run_length_specified = .false.
  character(len=*), parameter :: sub = 'rtm::timemgr_init'

  !
  dtime = real(dtime_in)
  call timemgr_spmdbcast( )

  ! Initalize calendar 
  call init_calendar()

  ! Initalize start date.
  if ( start_ymd == uninit_int ) then
     write(iulog,*)sub,': start_ymd must be specified '
     call shr_sys_abort
  end if
  if ( start_tod == uninit_int ) then
     write(iulog,*)sub,': start_tod must be specified '
     call shr_sys_abort
  end if
  start_date = TimeSetymd( start_ymd, start_tod, "start_date" )

  ! Initialize current date
  curr_date = start_date

  ! Initalize stop date.
  stop_date = TimeSetymd( 99991231, stop_tod, "stop_date" )

  call ESMF_TimeIntervalSet( step_size, s=dtime, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting step_size')

  call ESMF_TimeIntervalSet( day_step_size, d=1, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting day_step_size')

  if ( stop_ymd /= uninit_int ) then
     current = TimeSetymd( stop_ymd, stop_tod, "stop_date" )
     if ( current < stop_date ) stop_date = current
     run_length_specified = .true.
  end if
  if ( nelapse /= uninit_int ) then
     if ( nelapse >= 0 ) then
        current = curr_date + step_size*nelapse
     else
        current = curr_date - day_step_size*nelapse
     end if
     if ( current < stop_date ) stop_date = current
     run_length_specified = .true.
  end if
  if ( .not. run_length_specified ) then
     call shr_sys_abort (sub//': Must specify stop_ymd or nelapse')
  end if

  ! Error check 
  if ( stop_date <= start_date ) then
     write(iulog,*)sub, ': stop date must be specified later than start date: '
     call ESMF_TimeGet( start_date, yy=yr, mm=mon, dd=day, s=tod )
     write(iulog,*) ' Start date (yr, mon, day, tod): ', yr, mon, day, tod
     call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
     write(iulog,*) ' Stop date  (yr, mon, day, tod): ', yr, mon, day, tod
     call shr_sys_abort
  end if
  if ( curr_date >= stop_date ) then
     write(iulog,*)sub, ': stop date must be specified later than current date: '
     call ESMF_TimeGet( curr_date, yy=yr, mm=mon, dd=day, s=tod )
     write(iulog,*) ' Current date (yr, mon, day, tod): ', yr, mon, day, tod
     call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
     write(iulog,*) ' Stop date    (yr, mon, day, tod): ', yr, mon, day, tod
     call shr_sys_abort
  end if

  ! Initalize reference date for time coordinate.
  if ( ref_ymd /= uninit_int ) then
     ref_date = TimeSetymd( ref_ymd, ref_tod, "ref_date" )
  else
     ref_date = start_date
  end if
     
  ! Initialize clock
  call init_clock( start_date, ref_date, curr_date, stop_date )

  ! Print configuration summary to log file (stdout).
  if (masterproc) call timemgr_print()

  timemgr_set = .true.

end subroutine timemgr_init

!=========================================================================================

subroutine init_clock( start_date, ref_date, curr_date, stop_date )

  ! Initialize the clock based on the start_date, ref_date, and curr_date
  ! as well as the settings from the namelist specifying the time to stop
  !
  type(ESMF_Time), intent(in) :: start_date  ! start date for run
  type(ESMF_Time), intent(in) :: ref_date    ! reference date for time coordinate
  type(ESMF_Time), intent(in) :: curr_date   ! current date (equal to start_date)
  type(ESMF_Time), intent(in) :: stop_date   ! stop date for run
  !
  character(len=*), parameter :: sub = 'rtm::init_clock'
  type(ESMF_TimeInterval) :: step_size       ! timestep size
  type(ESMF_Time) :: current     ! current date (from clock)
  integer :: yr, mon, day, tod   ! Year, month, day, and second as integers
  integer :: rc                  ! return code
  !
  call ESMF_TimeIntervalSet( step_size, s=dtime, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting step_size')

  ! Initialize the clock

  tm_clock = ESMF_ClockCreate(name="RTM Time-manager clock", timeStep=step_size, startTime=start_date, &
       stopTime=stop_date, refTime=ref_date, rc=rc)
  call chkrc(rc, sub//': error return from ESMF_ClockSetup')

  ! Advance clock to the current time (in case of a restart)

  call ESMF_ClockGet(tm_clock, currTime=current, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_ClockGet')
  do while( curr_date > current )
     call ESMF_ClockAdvance( tm_clock, rc=rc )
     call chkrc(rc, sub//': error return from ESMF_ClockAdvance')
     call ESMF_ClockGet(tm_clock, currTime=current )
     call chkrc(rc, sub//': error return from ESMF_ClockGet')
  end do
end subroutine init_clock

!=========================================================================================

function TimeSetymd( ymd, tod, desc )


  ! Set the time by an integer as YYYYMMDD and integer seconds in the day
  !
  integer, intent(in) :: ymd            ! Year, month, day YYYYMMDD
  integer, intent(in) :: tod            ! Time of day in seconds
  character(len=*), intent(in) :: desc  ! Description of time to set
  !
  type(ESMF_Time) :: TimeSetymd         ! Return value
  !
  character(len=*), parameter :: sub = 'rtm::TimeSetymd'
  integer :: yr, mon, day          ! Year, month, day as integers
  integer :: rc                    ! return code
  !
  if ( (ymd < 0) .or. (tod < 0) .or. (tod > isecspday) )then
     write(iulog,*) sub//': error yymmdd is a negative number or time-of-day out of bounds', &
          ymd, tod
     call shr_sys_abort
  end if
  yr  = ymd / 10000
  mon = (ymd - yr*10000) / 100
  day =  ymd - yr*10000 - mon*100
  call ESMF_TimeSet( TimeSetymd, yy=yr, mm=mon, dd=day, s=tod, &
       calendar=tm_cal, rc=rc)
  call chkrc(rc, sub//': error return from ESMF_TimeSet: setting '//trim(desc))
end function TimeSetymd

!=========================================================================================

integer function TimeGetymd( date, tod )

  ! Get the date and time of day in ymd from ESMF Time.
  !
  type(ESMF_Time), intent(inout) :: date ! Input date to convert to ymd
  integer, intent(out), optional :: tod  ! Time of day in seconds
  !
  character(len=*), parameter :: sub = 'rtm::TimeGetymd'
  integer :: yr, mon, day
  integer :: rc                          ! return code
  !
  call ESMF_TimeGet( date, yy=yr, mm=mon, dd=day, rc=rc)
  call chkrc(rc, sub//': error return from ESMF_TimeGet')
  TimeGetymd = yr*10000 + mon*100 + day
  if ( present( tod ) )then
     call ESMF_TimeGet( date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
     call chkrc(rc, sub//': error return from ESMF_TimeGet')
  end if
  if ( yr < 0 )then
     write(iulog,*) sub//': error year is less than zero', yr
     call shr_sys_abort
  end if
end function TimeGetymd

!=========================================================================================

subroutine timemgr_restart(ncid, flag)

  ! Read/Write information needed on restart to a netcdf file. 
  !
  type(file_desc_t), intent(inout) :: ncid  ! netcdf id
  character(len=*) , intent(in) :: flag     ! 'read' or 'write'
  !
  logical :: run_length_specified = .false.
  integer :: rc                  ! return code
  integer :: yr, mon, day, tod   ! Year, month, day, and second as integers
  logical :: readvar             ! determine if variable is on initial file
  integer :: rst_caltype         ! calendar type
  type(ESMF_Time) :: start_date  ! start date for run
  type(ESMF_Time) :: stop_date   ! stop date for run
  type(ESMF_Time) :: ref_date    ! reference date for run
  type(ESMF_Time) :: curr_date   ! date of data in restart file
  type(ESMF_Time) :: current     ! current date (from clock)
  type(ESMF_TimeInterval) :: day_step_size ! day step size
  type(ESMF_TimeInterval) :: step_size     ! timestep size
  integer, parameter :: noleap = 1
  integer, parameter :: gregorian = 2
  character(len=135) :: varname
  character(len=len(calendar)) :: cal
  character(len=*), parameter :: sub = 'timemgr_restart'
  !
  if (flag == 'write') then
     rst_calendar  = calendar
  else if (flag == 'read') then
     calendar = rst_calendar
  end if
  varname = 'timemgr_rst_type'
  if (flag == 'define') then
     call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
          long_name='calendar type', units='unitless', flag_meanings=(/ "NO_LEAP_C", "GREGORIAN" /), &
          flag_values=(/ noleap, gregorian /), ifill_value=uninit_int )
  else if (flag == 'read' .or. flag == 'write') then
     if (flag== 'write') then
        cal = to_upper(calendar)
        if ( trim(cal) == 'NO_LEAP' ) then
           rst_caltype = noleap
        else if ( trim(cal) == 'GREGORIAN' ) then
           rst_caltype = gregorian
        else
           call shr_sys_abort(sub//'ERROR: unrecognized calendar specified= '//trim(calendar))
        end if
     end if
     call ncd_io(varname=varname, data=rst_caltype, &
          ncid=ncid, flag=flag, readvar=readvar)
     if (flag=='read' .and. .not. readvar) then
        if (is_restart()) then
           call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
        end if
     end if
     if (flag == 'read') then
        if ( rst_caltype == noleap ) then
           calendar = 'NO_LEAP'
        else if ( rst_caltype == gregorian ) then
           calendar = 'GREGORIAN'
        else
           write(iulog,*)sub,': unrecognized calendar type in restart file: ',rst_caltype
           call shr_sys_abort( sub//'ERROR: bad calendar type in restart file')
        end if
     end if
  end if

  if (flag == 'write') then
     call ESMF_ClockGet( tm_clock, startTime=start_date, currTime=curr_date, refTime=ref_date, rc=rc )
     call chkrc(rc, sub//': error return from ESMF_ClockGet')
     rst_step_sec  = dtime
     rst_start_ymd = TimeGetymd( start_date, tod=rst_start_tod )
     rst_ref_ymd   = TimeGetymd( ref_date,   tod=rst_ref_tod   )
     rst_curr_ymd  = TimeGetymd( curr_date,  tod=rst_curr_tod  )
  end if
  
  varname = 'timemgr_rst_step_sec'
  if (flag == 'define') then
     call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
          long_name='seconds component of timestep size', units='sec', nvalid_range=(/0,isecspday/), ifill_value=uninit_int)
  else if (flag == 'read' .or. flag == 'write') then
     call ncd_io(varname=varname, data=rst_step_sec, &
          ncid=ncid, flag=flag, readvar=readvar)
     if (flag=='read' .and. .not. readvar) then
        if (is_restart()) then
           call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
        end if
     end if
     if ( rst_step_sec < 0 .or. rst_step_sec > isecspday ) then
        call shr_sys_abort( sub//'ERROR: '//trim(varname)//' out of range')
     end if
  end if

  varname = 'timemgr_rst_start_ymd'
  if (flag == 'define') then
     call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
          long_name='start date', units='YYYYMMDD', ifill_value=uninit_int)
  else if (flag == 'read' .or. flag == 'write') then
     call ncd_io(varname=varname, data=rst_start_ymd, &
          ncid=ncid, flag=flag, readvar=readvar)
     if (flag=='read' .and. .not. readvar) then
        if (is_restart()) then
           call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
        end if
     end if
  end if

  varname = 'timemgr_rst_start_tod'
  if (flag == 'define') then
     call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
          long_name='start time of day', units='sec', nvalid_range=(/0,isecspday/), ifill_value=uninit_int)
  else if (flag == 'read' .or. flag == 'write') then
     call ncd_io(varname=varname, data=rst_start_tod, &
          ncid=ncid, flag=flag, readvar=readvar)
     if (flag=='read' .and. .not. readvar) then
        if (is_restart()) then
           call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
        end if
     end if
     if ( rst_start_tod < 0 .or. rst_start_tod > isecspday ) then
        call shr_sys_abort( sub//'ERROR: '//trim(varname)//' out of range')
     end if
  end if

  varname = 'timemgr_rst_ref_ymd'
  if (flag == 'define') then
     call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
          long_name='reference date', units='YYYYMMDD', ifill_value=uninit_int)
  else if (flag == 'read' .or. flag == 'write') then
     call ncd_io(varname=varname, data=rst_ref_ymd, &
          ncid=ncid, flag=flag, readvar=readvar)
     if (flag=='read' .and. .not. readvar) then
        if (is_restart()) then
           call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
        end if
     end if
  end if

  varname = 'timemgr_rst_ref_tod'
  if (flag == 'define') then
     call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
          long_name='reference time of day', units='sec', nvalid_range=(/0,isecspday/), ifill_value=uninit_int)
  else if (flag == 'read' .or. flag == 'write') then
     call ncd_io(varname=varname, data=rst_ref_tod, &
          ncid=ncid, flag=flag, readvar=readvar)
     if (flag=='read' .and. .not. readvar) then
        if (is_restart()) then
           call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
        end if
     end if
     if ( rst_start_tod < 0 .or. rst_start_tod > isecspday ) then
        call shr_sys_abort( sub//'ERROR: '//trim(varname)//' out of range')
     end if
  end if

  varname = 'timemgr_rst_curr_ymd'
  if (flag == 'define') then
     call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
          long_name='current date', units='YYYYMMDD', ifill_value=uninit_int)
  else if (flag == 'read' .or. flag == 'write') then
     call ncd_io(varname=varname, data=rst_curr_ymd, &
          ncid=ncid, flag=flag, readvar=readvar)
     if (flag=='read' .and. .not. readvar) then
        if (is_restart()) then
           call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
        end if
     end if
  end if

  varname = 'timemgr_rst_curr_tod'
  if (flag == 'define') then
     call ncd_defvar(ncid=ncid, varname=varname, xtype=ncd_int,  &
          long_name='current time of day', units='sec', nvalid_range=(/0,isecspday/), ifill_value=uninit_int )
  else if (flag == 'read' .or. flag == 'write') then
     call ncd_io(varname=varname, data=rst_curr_tod, &
          ncid=ncid, flag=flag, readvar=readvar)
     if (flag=='read' .and. .not. readvar) then
        if (is_restart()) then
           call shr_sys_abort( sub//'ERROR: '//trim(varname)//' not on file')
        end if
     end if
     if ( rst_curr_tod < 0 .or. rst_curr_tod > isecspday ) then
        call shr_sys_abort( sub//'ERROR: '//trim(varname)//' out of range')
     end if
  end if


  if (flag == 'read') then

     ! Restart the ESMF time manager using the synclock for ending date.
     call timemgr_spmdbcast( )
     
     ! Initialize calendar from restart info
     call init_calendar()
     
     ! Initialize the timestep from restart info
     dtime = rst_step_sec
     
     ! Initialize start date from restart info
     start_date = TimeSetymd( rst_start_ymd, rst_start_tod, "start_date" )
     
     ! Initialize current date from restart info
     curr_date = TimeSetymd( rst_curr_ymd, rst_curr_tod, "curr_date" )
     
     ! Initialize stop date from sync clock or namelist input
     stop_date = TimeSetymd( 99991231, stop_tod, "stop_date" )
     
     call ESMF_TimeIntervalSet( step_size, s=dtime, rc=rc )
     call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting step_size')
     
     call ESMF_TimeIntervalSet( day_step_size, d=1, rc=rc )
     call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet: setting day_step_size')
     
     if    ( stop_ymd /= uninit_int ) then
        current = TimeSetymd( stop_ymd, stop_tod, "stop_date" )
        if ( current < stop_date ) stop_date = current
        run_length_specified = .true.
     else if ( nelapse /= uninit_int ) then
        if ( nelapse >= 0 ) then
           current = curr_date + step_size*nelapse
        else
           current = curr_date - day_step_size*nelapse
        end if
        if ( current < stop_date ) stop_date = current
        run_length_specified = .true.
     end if
     if ( .not. run_length_specified ) then
        call shr_sys_abort (sub//': Must specify stop_ymd or nelapse')
     end if
     
     ! Error check
     if ( stop_date <= start_date ) then
        write(iulog,*)sub, ': stop date must be specified later than start date: '
        call ESMF_TimeGet( start_date, yy=yr, mm=mon, dd=day, s=tod )
        write(iulog,*) ' Start date (yr, mon, day, tod): ', yr, mon, day, tod
        call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
        write(iulog,*) ' Stop date  (yr, mon, day, tod): ', yr, mon, day, tod
        call shr_sys_abort
     end if
     if ( curr_date >= stop_date ) then
        write(iulog,*)sub, ': stop date must be specified later than current date: '
        call ESMF_TimeGet( curr_date, yy=yr, mm=mon, dd=day, s=tod )
        write(iulog,*) ' Current date (yr, mon, day, tod): ', yr, mon, day, tod
        call ESMF_TimeGet( stop_date, yy=yr, mm=mon, dd=day, s=tod )
        write(iulog,*) ' Stop date    (yr, mon, day, tod): ', yr, mon, day, tod
        call shr_sys_abort
     end if
     
     ! Initialize ref date from restart info
     ref_date = TimeSetymd( rst_ref_ymd, rst_ref_tod, "ref_date" )
     
     ! Initialize clock 
     call init_clock( start_date, ref_date, curr_date, stop_date )
     
     ! Set flag that this is the first timestep of the restart run.
     tm_first_restart_step = .true.
     
     ! Print configuration summary to log file (stdout).
     if (masterproc) call timemgr_print()

     timemgr_set = .true.

  end if

end subroutine timemgr_restart

!=========================================================================================

subroutine init_calendar( )

  !---------------------------------------------------------------------------------
  ! Initialize calendar
  !
  ! Local variables
  !
  character(len=*), parameter :: sub = 'rtm::init_calendar'
  type(ESMF_CalKind_Flag) :: cal_type        ! calendar type
  character(len=len(calendar)) :: caltmp
  integer :: rc                              ! return code
  !---------------------------------------------------------------------------------

  caltmp = to_upper(calendar)
  if ( trim(caltmp) == 'NO_LEAP' ) then
     cal_type = ESMF_CALKIND_NOLEAP
  else if ( trim(caltmp) == 'GREGORIAN' ) then
     cal_type = ESMF_CALKIND_GREGORIAN
  else
     write(iulog,*)sub,': unrecognized calendar specified: ',calendar
     call shr_sys_abort
  end if
  tm_cal = ESMF_CalendarCreate( name=caltmp, calkindflag=cal_type, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_CalendarSet')

end subroutine init_calendar

!=========================================================================================

subroutine timemgr_print()

  !---------------------------------------------------------------------------------
  character(len=*), parameter :: sub = 'rtm::timemgr_print'
  integer :: rc
  integer :: yr, mon, day
  integer :: &                   ! Data required to restart time manager:
       nstep     = uninit_int,  &! current step number
       step_sec  = uninit_int,  &! timestep size seconds
       start_yr  = uninit_int,  &! start year
       start_mon = uninit_int,  &! start month
       start_day = uninit_int,  &! start day of month
       start_tod = uninit_int,  &! start time of day
       stop_yr   = uninit_int,  &! stop year
       stop_mon  = uninit_int,  &! stop month
       stop_day  = uninit_int,  &! stop day of month
       stop_tod  = uninit_int,  &! stop time of day
       ref_yr    = uninit_int,  &! reference year
       ref_mon   = uninit_int,  &! reference month
       ref_day   = uninit_int,  &! reference day of month
       ref_tod   = uninit_int,  &! reference time of day
       curr_yr   = uninit_int,  &! current year
       curr_mon  = uninit_int,  &! current month
       curr_day  = uninit_int,  &! current day of month
       curr_tod  = uninit_int    ! current time of day
  integer(ESMF_KIND_I8) :: step_no
  type(ESMF_Time) :: start_date! start date for run
  type(ESMF_Time) :: stop_date ! stop date for run
  type(ESMF_Time) :: curr_date ! date of data in restart file
  type(ESMF_Time) :: ref_date  ! reference date
  type(ESMF_TimeInterval) :: step ! Time-step
  !---------------------------------------------------------------------------------

  call ESMF_ClockGet( tm_clock, startTime=start_date, currTime=curr_date, &
       refTime=ref_date, stopTime=stop_date, timeStep=step, &
       advanceCount=step_no, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_ClockGet')
  nstep = step_no

  write(iulog,*)' ******** RTM Time Manager Configuration ********'

  call ESMF_TimeIntervalGet( step, s=step_sec, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet')

  call ESMF_TimeGet( start_date, yy=start_yr, mm=start_mon, dd=start_day, &
       s=start_tod, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_TimeGet')
  call ESMF_TimeGet( stop_date, yy=stop_yr, mm=stop_mon, dd=stop_day, &
       s=stop_tod, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_TimeGet')
  call ESMF_TimeGet( ref_date, yy=ref_yr, mm=ref_mon, dd=ref_day, s=ref_tod, &
       rc=rc )
  call chkrc(rc, sub//': error return from ESMF_TimeGet')
  call ESMF_TimeGet( curr_date, yy=curr_yr, mm=curr_mon, dd=curr_day, &
       s=curr_tod, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_TimeGet')

  write(iulog,*)'  Calendar type:            ',trim(calendar)
  write(iulog,*)'  Timestep size (seconds):  ', step_sec
  write(iulog,*)'  Start date (yr mon day tod):     ', start_yr, start_mon, &
       start_day, start_tod
  write(iulog,*)'  Stop date (yr mon day tod):      ', stop_yr, stop_mon, &
       stop_day, stop_tod
  write(iulog,*)'  Reference date (yr mon day tod): ', ref_yr, ref_mon, &
       ref_day, ref_tod
  write(iulog,*)'  Current step number:      ', nstep
  write(iulog,*)'  Current date (yr mon day tod):   ', curr_yr, curr_mon, &
       curr_day, curr_tod

  write(iulog,*)' ************************************************'

end subroutine timemgr_print

!=========================================================================================

subroutine advance_timestep(nsub)

  ! Increment the timestep number.

  integer, intent(in) :: nsub  ! subcycling timestep count

  type(ESMF_Time) :: curr_date
  integer :: yy0, mm0, dd0, tod0   ! before advance
  integer :: yy1, mm1, dd1, tod1   ! after advance
  integer :: rc
  character(len=*), parameter :: sub = 'rtm::advance_timestep'

  tm_first_restart_step = .false.
  new_year = .false.
  new_month = .false.
  new_day = .false.

  if (nsub == 1) then
     call get_curr_date(yy0, mm0, dd0, tod0)
     call ESMF_ClockAdvance( tm_clock, rc=rc )
     call chkrc(rc, sub//': error return from ESMF_ClockAdvance')
     call get_curr_date(yy1, mm1, dd1, tod1)

     if (tod0 == 0._r8) then
        new_day = .true.
        if (dd1 == 1) then
           new_month = .true.
           if (mm1 == 1) then
              new_year = .true.
           endif
        endif
     else
        if (tod1 /= 0._r8) then
           if (yy1 /= yy0) new_year = .true.
           if (mm1 /= mm0) new_month = .true.
           if (dd1 /= dd0) new_day = .true.
        endif
     endif
     write(iulog,'(2a,4i6,3l4)') sub," yy,mm,dd = ",yy1,mm1,dd1,tod1,new_year,new_month,new_day
  endif

end subroutine advance_timestep

!=========================================================================================

subroutine get_clock( clock )

  ! Return the ESMF clock

  type(ESMF_Clock), intent(inout) :: clock

  character(len=*), parameter :: sub = 'rtm::get_clock'
  type(ESMF_TimeInterval) :: step_size
  type(ESMF_Time) :: start_date, stop_date, ref_date
  integer :: rc

  call ESMF_ClockGet( tm_clock, timeStep=step_size, startTime=start_date, &
                      stoptime=stop_date, reftime=ref_date, rc=rc )
  call chkrc(rc, sub//': error return from ESMF_ClockGet')
  call ESMF_ClockSet(clock, timeStep=step_size, startTime=start_date, &
                            stoptime=stop_date, reftime=ref_date, rc=rc)
  call chkrc(rc, sub//': error return from ESMF_ClockSet')

end subroutine get_clock

!=========================================================================================

integer function get_step_size()

  ! Return the step size in seconds.
  
  character(len=*), parameter :: sub = 'rtm::get_step_size'
  type(ESMF_TimeInterval) :: step_size       ! timestep size
  integer :: rc
  
  call ESMF_ClockGet(tm_clock, timeStep=step_size, rc=rc)
  call chkrc(rc, sub//': error return from ESMF_ClockGet')

  call ESMF_TimeIntervalGet(step_size, s=get_step_size, rc=rc)
  call chkrc(rc, sub//': error return from ESMF_ClockTimeIntervalGet')
  
end function get_step_size

!=========================================================================================

integer function get_nstep()

  ! Return the timestep number.

   character(len=*), parameter :: sub = 'rtm::get_nstep'
   integer :: rc
   integer(ESMF_KIND_I8) :: step_no

   call ESMF_ClockGet(tm_clock, advanceCount=step_no, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   get_nstep = step_no

end function get_nstep

!=========================================================================================

subroutine get_curr_date(yr, mon, day, tod, offset)

  !-----------------------------------------------------------------------------------------
  ! Return date components valid at end of current timestep with an optional
  ! offset (positive or negative) in seconds.
  
  integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.

   character(len=*), parameter :: sub = 'rtm::get_curr_date'
   integer :: rc
   type(ESMF_Time) :: date
   type(ESMF_TimeInterval) :: off
   !-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, currTime=date, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   if (present(offset)) then
      if (offset > 0) then
         call ESMF_TimeIntervalSet( off, s=offset, rc=rc )
         call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
         date = date + off
      else if (offset < 0) then
         call ESMF_TimeIntervalSet( off, s=-offset, rc=rc )
         call chkrc(rc, sub//': error return from ESMF_TimeIntervalSet')
         date = date - off
      end if
   end if

   call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

end subroutine get_curr_date

!=========================================================================================

subroutine get_prev_date(yr, mon, day, tod)

! Return date components valid at beginning of current timestep.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

! Local variables
   character(len=*), parameter :: sub = 'rtm::get_prev_date'
   integer :: rc
   type(ESMF_Time) :: date
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet(tm_clock, prevTime=date, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

end subroutine get_prev_date

!=========================================================================================

subroutine get_start_date(yr, mon, day, tod)

   ! Return date components valid at beginning of initial run.
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

   character(len=*), parameter :: sub = 'rtm::get_start_date'
   integer :: rc
   type(ESMF_Time) :: date

   call ESMF_ClockGet(tm_clock, startTime=date, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

end subroutine get_start_date

!=========================================================================================

subroutine get_ref_date(yr, mon, day, tod)

! Return date components of the reference date.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

! Local variables
   character(len=*), parameter :: sub = 'rtm::get_ref_date'
   integer :: rc
   type(ESMF_Time) :: date
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet(tm_clock, refTime=date, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   call ESMF_TimeGet(date, yy=yr, mm=mon, dd=day, s=tod, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeGet')

end subroutine get_ref_date

!=========================================================================================

subroutine get_curr_time(days, seconds)

! Return time components valid at end of current timestep.
! Current time is the time interval between the current date and the reference date.

! Arguments
   integer, intent(out) ::&
      days,   &! number of whole days in time interval
      seconds  ! remaining seconds in time interval

! Local variables
   character(len=*), parameter :: sub = 'rtm::get_curr_time'
   integer :: rc
   type(ESMF_Time) :: cdate, rdate
   type(ESMF_TimeInterval) :: diff
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet( tm_clock, currTime=cdate, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   call ESMF_ClockGet( tm_clock, refTime=rdate, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')

   diff = cdate - rdate

   call ESMF_TimeIntervalGet(diff, d=days, s=seconds, rc=rc)
   call chkrc(rc, sub//': error return from ESMF_TimeIntervalGet')

end subroutine get_curr_time

!=========================================================================================

subroutine get_prev_time(days, seconds)

! Return time components valid at beg of current timestep.
! prev time is the time interval between the prev date and the reference date.

! Arguments
   integer, intent(out) ::&
      days,   &! number of whole days in time interval
      seconds  ! remaining seconds in time interval

! Local variables
   character(len=*), parameter :: sub = 'rtm::get_prev_time'
   integer :: rc
   type(ESMF_Time) :: date, ref_date
   type(ESMF_TimeInterval) :: diff
!-----------------------------------------------------------------------------------------

   call ESMF_ClockGet(tm_clock, prevTime=date, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet for prevTime')
   call ESMF_ClockGet(tm_clock, refTime=ref_date, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet for refTime')
   diff = date - ref_date
   call ESMF_TimeIntervalGet( diff, d=days, s=seconds, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_TimeintervalGet')

end subroutine get_prev_time

!=========================================================================================

function get_calendar()

   ! Return calendar

   character(len=ESMF_MAXSTR) :: get_calendar

   get_calendar = calendar

end function get_calendar

!=========================================================================================
 
function is_end_curr_day()

   ! Return true if current timestep is last timestep in current day.
   logical :: is_end_curr_day

   integer ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

   call get_curr_date(yr, mon, day, tod)
   is_end_curr_day = (tod == 0)

end function is_end_curr_day

!=========================================================================================

logical function is_end_curr_month()

  ! Return true if current timestep is last timestep in current month.
   integer ::  yr, mon, day, tod     ! time of day (seconds past 0Z)

   call get_curr_date(yr, mon, day, tod)
   is_end_curr_month = (day == 1  .and.  tod == 0)

end function is_end_curr_month

!=========================================================================================

logical function is_first_step()

  ! Return true on first step of initial run only.
   character(len=*), parameter :: sub = 'rtm::is_first_step'
   integer :: rc
   integer :: nstep
   integer(ESMF_KIND_I8) :: step_no

   call ESMF_ClockGet( tm_clock, advanceCount=step_no, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   nstep = step_no
   is_first_step = (nstep == 0)

end function is_first_step

!=========================================================================================

logical function is_first_restart_step()

   ! Return true on first step of restart run only.
   is_first_restart_step = tm_first_restart_step

end function is_first_restart_step

!=========================================================================================

logical function is_last_step()

  ! Return true on last timestep.
   character(len=*), parameter :: sub = 'rtm::is_last_step'
   type(ESMF_Time) :: stop_date
   type(ESMF_Time) :: curr_date
   type(ESMF_TimeInterval) :: time_step
   integer :: rc

   call ESMF_ClockGet( tm_clock, stopTime=stop_date, &
                       currTime=curr_date, TimeStep=time_step, rc=rc )
   call chkrc(rc, sub//': error return from ESMF_ClockGet')
   if ( curr_date+time_step > stop_date ) then
      is_last_step = .true.
   else
      is_last_step = .false.
   end if

end function is_last_step

!=========================================================================================

subroutine chkrc(rc, mes)
   integer, intent(in)          :: rc   ! return code from time management library
   character(len=*), intent(in) :: mes  ! error message
   if ( rc == ESMF_SUCCESS ) return
   write(iulog,*) mes
   call shr_sys_abort ('CHKRC')
end subroutine chkrc

!=========================================================================================

function to_upper(str)

  ! Convert character string to upper case. Use achar and iachar intrinsics
  ! to ensure use of ascii collating sequence.
  character(len=*), intent(in) :: str ! String to convert to upper case
  character(len=len(str))      :: to_upper

  integer :: i                ! Index
  integer :: aseq             ! ascii collating sequence
  character(len=1) :: ctmp    ! Character temporary
  
  do i = 1, len(str)
     ctmp = str(i:i)
     aseq = iachar(ctmp)
     if ( aseq >= 97  .and.  aseq <= 122 ) ctmp = achar(aseq - 32)
     to_upper(i:i) = ctmp
  end do
  
end function to_upper

!=========================================================================================

logical function is_restart( )
  ! Determine if restart run
  if (nsrest == nsrContinue) then
     is_restart = .true.
  else
     is_restart = .false.
  end if
end function is_restart

!=========================================================================================

logical function is_new_year( )
  ! Determine if new year
  is_new_year = new_year
end function is_new_year

!=========================================================================================

logical function is_new_month( )
  ! Determine if new month
  is_new_month = new_month
end function is_new_month

!=========================================================================================

logical function is_new_day( )
  ! Determine if new day
  is_new_day = new_day
end function is_new_day

!=========================================================================================

subroutine timemgr_spmdbcast( )

  integer :: ier

  call mpi_bcast (dtime, 1, MPI_INTEGER, 0, mpicom_rof, ier)

end subroutine timemgr_spmdbcast

end module RtmTimeManager
