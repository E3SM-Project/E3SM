module seq_timemgr_mod

  ! !DESCRIPTION: A module to create derived types to manage time and clock information

  ! !USES:
  use ESMF, only  : ESMF_Clock, ESMF_Alarm, ESMF_Calendar
  use ESMF, only: operator(<), operator(/=), operator(+), operator(-), operator(*) , operator(>=)
  use ESMF, only: operator(<=), operator(>), operator(==)
  use med_constants_mod, only : CL, IN
  use med_constants_mod, only : seq_timemgr_noleap => med_constants_noleap
  use med_constants_mod, only: seq_timemgr_gregorian=>med_constants_gregorian

  implicit none
  private    ! default private

  ! MEMBER FUNCTIONS:

  ! --- Clock object methods --------------------------------------------------
  public :: seq_timemgr_clockInit              ! Setup the sync clock
  public :: seq_timemgr_EClockGetData          ! Get data from an ESMF clock
  public :: seq_timemgr_EClockDateInSync       ! compare EClock to ymd/tod
  public :: seq_timemgr_EclockPrint            ! Print ESMF clock information
  public :: seq_timemgr_alarmInit              ! initialize an alarm
  public :: seq_timemgr_alarmGet               ! get info about alarm
  public :: seq_timemgr_alarmSetOn             ! Turn an alarm on
  public :: seq_timemgr_alarmSetOff            ! Turn an alarm off
  public :: seq_timemgr_alarmIsOn              ! Is an alarm ringing
  public :: seq_timemgr_ETimeInit              ! Create ESMF_Time object
  public :: seq_timemgr_ETimeGet               ! Query ESMF_Time object

  ! --- For usability, built on interfaces above ---
  public :: seq_timemgr_restartAlarmIsOn       ! Is a restart alarm ringing
  public :: seq_timemgr_stopAlarmIsOn          ! Is a stop alarm ringing
  public :: seq_timemgr_historyAlarmIsOn       ! Is a history alarm ringing
  public :: seq_timemgr_pauseAlarmIsOn         ! Is a pause alarm ringing

  ! --- ESP components need to know about the state of other components
  public :: seq_timemgr_pause_active           ! Pause/resume is enabled
  public :: seq_timemgr_pause_component_index  ! Index of named component
  public :: seq_timemgr_pause_component_active ! .true. is comp should pause

  public :: seq_timemgr_clockPrint ! Print sync clock information

  private:: seq_timemgr_EClockInit
  private:: seq_timemgr_ESMFDebug

  ! PARAMETERS:

  ! History output types
  integer(IN)   ,public            :: seq_timemgr_histavg_type
  integer(IN)   ,public ,parameter :: seq_timemgr_type_other  = -1
  integer(IN)   ,public ,parameter :: seq_timemgr_type_never  =  1
  integer(IN)   ,public ,parameter :: seq_timemgr_type_nhour  =  2
  integer(IN)   ,public ,parameter :: seq_timemgr_type_nday   =  3
  integer(IN)   ,public ,parameter :: seq_timemgr_type_nmonth =  4
  integer(IN)   ,public ,parameter :: seq_timemgr_type_nyear  =  5

  ! Clock and alarm options
  character(len=*), private, parameter :: &
       seq_timemgr_optNONE           = "none"      , &
       seq_timemgr_optNever          = "never"     , &
       seq_timemgr_optNSteps         = "nsteps"    , &
       seq_timemgr_optNStep          = "nstep"     , &
       seq_timemgr_optNSeconds       = "nseconds"  , &
       seq_timemgr_optNSecond        = "nsecond"   , &
       seq_timemgr_optNMinutes       = "nminutes"  , &
       seq_timemgr_optNMinute        = "nminute"   , &
       seq_timemgr_optNHours         = "nhours"    , &
       seq_timemgr_optNHour          = "nhour"     , &
       seq_timemgr_optNDays          = "ndays"     , &
       seq_timemgr_optNDay           = "nday"      , &
       seq_timemgr_optNMonths        = "nmonths"   , &
       seq_timemgr_optNMonth         = "nmonth"    , &
       seq_timemgr_optNYears         = "nyears"    , &
       seq_timemgr_optNYear          = "nyear"     , &
       seq_timemgr_optMonthly        = "monthly"   , &
       seq_timemgr_optYearly         = "yearly"    , &
       seq_timemgr_optDate           = "date"      , &
       seq_timemgr_optIfdays0        = "ifdays0"   , &
       seq_timemgr_optEnd            = "end"       , &
       seq_timemgr_optGLCCouplingPeriod = "glc_coupling_period"

  ! Clock numbers
  integer(IN),private,parameter :: &
       seq_timemgr_nclock_drv  = 1, &
       seq_timemgr_nclock_atm  = 2, &
       seq_timemgr_nclock_lnd  = 3, &
       seq_timemgr_nclock_ocn  = 4, &
       seq_timemgr_nclock_ice  = 5, &
       seq_timemgr_nclock_glc  = 6, &
       seq_timemgr_nclock_wav  = 7, &
       seq_timemgr_nclock_rof  = 8, &
       seq_timemgr_nclock_esp  = 9, &
       max_clocks              = 9

  ! Clock names
  character(len=*), public,parameter :: &
       seq_timemgr_clock_drv  = 'seq_timemgr_clock_drv' , &
       seq_timemgr_clock_atm  = 'seq_timemgr_clock_atm' , &
       seq_timemgr_clock_lnd  = 'seq_timemgr_clock_lnd' , &
       seq_timemgr_clock_ocn  = 'seq_timemgr_clock_ocn' , &
       seq_timemgr_clock_ice  = 'seq_timemgr_clock_ice' , &
       seq_timemgr_clock_glc  = 'seq_timemgr_clock_glc' , &
       seq_timemgr_clock_wav  = 'seq_timemgr_clock_wav' , &
       seq_timemgr_clock_rof  = 'seq_timemgr_clock_rof' , &
       seq_timemgr_clock_esp  = 'seq_timemgr_clock_esp'

  ! Array of clock names
  character(len=8), private,parameter :: seq_timemgr_clocks(max_clocks) = &
       (/'drv     ','atm     ','lnd     ','ocn     ', &
         'ice     ','glc     ','wav     ','rof     ','esp     '/)

  ! Alarm numbers
  integer(IN), private,parameter :: &
       seq_timemgr_nalarm_restart    = 1 , & ! driver and component clock alarm
       seq_timemgr_nalarm_stop       = 2 , & ! driver and component clock alarm
       seq_timemgr_nalarm_datestop   = 3 , & ! driver and component clock alarm
       seq_timemgr_nalarm_history    = 4 , & ! driver and component clock alarm
       seq_timemgr_nalarm_tprof      = 5 , & ! driver and component clock alarm
       seq_timemgr_nalarm_histavg    = 6 , & ! driver and component clock alarm
       seq_timemgr_nalarm_pause      = 7 , &
       seq_timemgr_nalarm_barrier    = 8 , & ! driver and component clock alarm
       max_alarms = seq_timemgr_nalarm_barrier

  ! Alarm names
  character(len=*), public,parameter :: &
       seq_timemgr_alarm_restart    = 'seq_timemgr_alarm_restart ', &
       seq_timemgr_alarm_stop       = 'seq_timemgr_alarm_stop    ', &
       seq_timemgr_alarm_datestop   = 'seq_timemgr_alarm_datestop', &
       seq_timemgr_alarm_history    = 'seq_timemgr_alarm_history ', &
       seq_timemgr_alarm_tprof      = 'seq_timemgr_alarm_tprof   ', &
       seq_timemgr_alarm_histavg    = 'seq_timemgr_alarm_histavg ', &
       seq_timemgr_alarm_pause      = 'seq_timemgr_alarm_pause   ', &
       seq_timemgr_alarm_barrier    = 'seq_timemgr_alarm_barrier '

  ! Active pause - resume components
  logical, private :: pause_active(max_clocks) = .false.

  ! TYPES:

  type EClock_pointer     ! needed for array of pointers
     type(ESMF_Clock),pointer :: EClock => null()
  end type EClock_pointer

  public :: seq_timemgr_type         ! Wrapped clock object
  type seq_timemgr_type
     private
     type(EClock_pointer) :: ECP(max_clocks)               ! ESMF clocks, array of pointers
     type(ESMF_Alarm)     :: EAlarm(max_clocks,max_alarms) ! array of clock alarms
  end type seq_timemgr_type

  ! MODULE DATA

  type(seq_timemgr_type)      :: SyncClock                    ! array of all clocks & alarm
  type(ESMF_Calendar), target :: seq_timemgr_cal              ! calendar
  character(CL)               :: seq_timemgr_calendar         ! calendar string
  integer, parameter          :: SecPerDay = 86400            ! Seconds per day
  integer                     :: seq_timemgr_pause_sig_index  ! Index of pause comp with smallest dt
  logical                     :: seq_timemgr_esp_run_on_pause ! Run ESP component on pause cycle
  logical                     :: seq_timemgr_end_restart      ! write restarts at end of run?
  character(CL)               :: tmpstr
  integer                     :: dbrc
  integer, parameter          :: dbug_flag = 10
  character(len=*), parameter :: sp_str = 'str_undefined'
  character(len=*), parameter :: u_FILE_u =  __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine seq_timemgr_clockInit(driver, logunit, &
       EClock_drv, EClock_atm, EClock_lnd, EClock_ocn, &
       EClock_ice, Eclock_glc, Eclock_rof, EClock_wav, Eclock_esp, rc)

    ! !DESCRIPTION:  Initializes clock
    use med_constants_mod     , only : CS, CL, IN
    use NUOPC                 , only : NUOPC_CompAttributeGet
    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_ClockSet, ESMF_CalendarCreate, ESMF_FAILURE
    use ESMF                  , only : ESMF_Time, ESMF_TimeInterval, ESMF_CalKind_Flag, ESMF_VM, ESMF_VMGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU
    use ESMF                  , only : ESMF_LogFoundError, ESMF_TimeIntervalSet, ESMF_AlarmGet
    use ESMF                  , only : ESMF_CALKIND_NOLEAP, ESMF_CALKIND_GREGORIAN, ESMF_CalKind_Flag
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_mpi_mod           , only : shr_mpi_bcast
    use shr_sys_mod           , only : shr_sys_abort
    use shr_cal_mod           , only : shr_cal_calendarName
    use shr_file_mod          , only : shr_file_getUnit, shr_file_freeUnit
    use netcdf                , only : nf90_open, nf90_nowrite, nf90_noerr, nf90_inq_varid, nf90_get_var, nf90_close

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_GridComp),     intent(inout) :: driver
    integer            ,     intent(in)    :: logunit
    type(ESMF_clock),target, intent(in)    :: EClock_drv   ! drv clock
    type(ESMF_clock),target, intent(in)    :: EClock_atm   ! atm clock
    type(ESMF_clock),target, intent(in)    :: EClock_lnd   ! lnd clock
    type(ESMF_clock),target, intent(in)    :: EClock_ocn   ! ocn clock
    type(ESMF_clock),target, intent(in)    :: EClock_ice   ! ice clock
    type(ESMF_clock),target, intent(in)    :: EClock_glc   ! glc clock
    type(ESMF_clock),target, intent(in)    :: EClock_rof   ! rof clock
    type(ESMF_clock),target, intent(in)    :: EClock_wav   ! wav clock
    type(ESMF_clock),target, intent(in)    :: EClock_esp   ! esp clock
    integer                , intent(out)   :: rc

    !----- local -----
    integer                 :: mpicom       ! MPI communicator
    logical                 :: mastertask
    logical                 :: read_restart
    character(CL)           :: restart_file
    character(CL)           :: restart_pfile
    character(CL)           :: cvalue
    type(ESMF_Time)         :: StartTime            ! Start time
    type(ESMF_Time)         :: RefTime              ! Reference time
    type(ESMF_Time)         :: CurrTime             ! Current time
    type(ESMF_Time)         :: StopTime1            ! Stop time
    type(ESMF_Time)         :: StopTime2            ! Stop time
    type(ESMF_TimeInterval) :: TimeStep             ! Clock time-step
    type(ESMF_CalKind_Flag) :: esmf_caltype         ! local esmf calendar
    integer                 :: n, i                 ! index
    logical                 :: found
    integer                 :: iam, unitn
    integer                 :: min_dt               ! smallest time step
    integer                 :: dtime(max_clocks)    ! time-step to use
    character(CS)           :: calendar             ! Calendar type
    character(CS)           :: stop_option          ! Stop option units
    integer(IN)             :: stop_n               ! Number until stop
    integer(IN)             :: stop_ymd             ! Stop date (YYYYMMDD)
    integer(IN)             :: stop_tod             ! Stop time-of-day
    character(CS)           :: restart_option       ! Restart option units
    integer(IN)             :: restart_n            ! Number until restart interval
    integer(IN)             :: restart_ymd          ! Restart date (YYYYMMDD)
    character(CS)           :: pause_option         ! Pause option units
    integer(IN)             :: pause_n              ! Number between pause intervals
    character(CS)           :: pause_component_list ! Pause - resume components
    character(CS)           :: history_option       ! History option units
    integer(IN)             :: history_n            ! Number until history interval
    integer(IN)             :: history_ymd          ! History date (YYYYMMDD)
    character(CS)           :: histavg_option       ! Histavg option units
    integer(IN)             :: histavg_n            ! Number until histavg interval
    integer(IN)             :: histavg_ymd          ! Histavg date (YYYYMMDD)
    character(CS)           :: barrier_option       ! Barrier option units
    integer(IN)             :: barrier_n            ! Number until barrier interval
    integer(IN)             :: barrier_ymd          ! Barrier date (YYYYMMDD)
    character(CS)           :: tprof_option         ! tprof option units
    integer(IN)             :: tprof_n              ! Number until tprof interval
    integer(IN)             :: tprof_ymd            ! tprof date (YYYYMMDD)
    integer(IN)             :: start_ymd            ! Start date (YYYYMMDD)
    integer(IN)             :: start_tod            ! Start time of day (seconds)
    integer(IN)             :: curr_ymd             ! Current ymd (YYYYMMDD)
    integer(IN)             :: curr_tod             ! Current tod (seconds)
    integer(IN)             :: ref_ymd              ! Reference date (YYYYMMDD)
    integer(IN)             :: ref_tod              ! Reference time of day (seconds)
    integer(IN)             :: atm_cpl_dt           ! Atmosphere coupling interval
    integer(IN)             :: lnd_cpl_dt           ! Land coupling interval
    integer(IN)             :: ice_cpl_dt           ! Sea-Ice coupling interval
    integer(IN)             :: ocn_cpl_dt           ! Ocean coupling interval
    integer(IN)             :: glc_cpl_dt           ! Glc coupling interval
    integer(IN)             :: rof_cpl_dt           ! Runoff coupling interval
    integer(IN)             :: wav_cpl_dt           ! Wav coupling interval
    integer(IN)             :: esp_cpl_dt           ! Esp coupling interval
    character(CS)           :: glc_avg_period       ! Glc avering coupling period
    logical                 :: esp_run_on_pause     ! Run ESP on pause cycle
    logical                 :: end_restart          ! Write restart at end of run
    integer(IN)             :: ierr                 ! Return code
    integer(IN)             :: status, ncid, varid  ! netcdf stuff
    type(ESMF_VM)           :: vm
    character(len=*), parameter :: F0A = "(2A,A)"
    character(len=*), parameter :: F0I = "(2A,I10)"
    character(len=*), parameter :: F0L = "(2A,L3)"
    character(len=*), parameter :: subname = '(seq_timemgr_clockInit) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    call ESMF_GridCompGet(driver, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpicom, localPet=iam, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (iam == 0) then
       mastertask=.true.
    else
       mastertask = .false.
    end if

    call NUOPC_CompAttributeGet(driver, name='read_restart', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart

    !---------------------------------------------------------------------------
    ! Get clock config attributes
    !---------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(driver, name="calendar", value=calendar, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="stop_option", value=stop_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="stop_n", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) stop_n

    call NUOPC_CompAttributeGet(driver, name="stop_ymd", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) stop_ymd

    call NUOPC_CompAttributeGet(driver, name="stop_tod", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) stop_tod

    call NUOPC_CompAttributeGet(driver, name="restart_option", value=restart_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="restart_n", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) restart_n

    call NUOPC_CompAttributeGet(driver, name="restart_ymd", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) restart_ymd

    call NUOPC_CompAttributeGet(driver, name="pause_option", value=pause_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="pause_n", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) pause_n

    ! TODO: currently this is not in namelist_definition_drv.xml
    ! call NUOPC_CompAttributeGet(driver, name="pause_component_list", value=pause_component_list, rc=rc)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    pause_component_list = ' '

    call NUOPC_CompAttributeGet(driver, name="history_option", value=history_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="history_n", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) history_n

    call NUOPC_CompAttributeGet(driver, name="history_ymd", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) history_ymd

    call NUOPC_CompAttributeGet(driver, name="histavg_option", value=histavg_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="histavg_n", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) histavg_n

    call NUOPC_CompAttributeGet(driver, name="histavg_ymd", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) histavg_ymd

    call NUOPC_CompAttributeGet(driver, name="barrier_option", value=barrier_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="barrier_n", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) barrier_n

    call NUOPC_CompAttributeGet(driver, name="barrier_ymd", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) barrier_ymd

    call NUOPC_CompAttributeGet(driver, name="tprof_option", value=tprof_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(driver, name="tprof_n", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) tprof_n

    call NUOPC_CompAttributeGet(driver, name="tprof_ymd", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) tprof_ymd

    call NUOPC_CompAttributeGet(driver, name="start_ymd", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) start_ymd

    call NUOPC_CompAttributeGet(driver, name="start_tod", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) start_tod

    call NUOPC_CompAttributeGet(driver, name="ref_ymd", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ref_ymd

    call NUOPC_CompAttributeGet(driver, name="ref_tod", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ref_tod

    ! These do not appear in namelist_definition_drv.xml and its not clear they should
    ! call NUOPC_CompAttributeGet(driver, name="curr_ymd", value=cvalue, rc=rc)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    ! read(cvalue,*) curr_ymd
    curr_ymd = 0.0

    ! call NUOPC_CompAttributeGet(driver, name="curr_tod", value=cvalue, rc=rc)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    ! read(cvalue,*) curr_tod
    curr_tod = 0.0

    call NUOPC_CompAttributeGet(driver, name="atm_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) atm_cpl_dt

    call NUOPC_CompAttributeGet(driver, name="lnd_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) lnd_cpl_dt

    call NUOPC_CompAttributeGet(driver, name="ice_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ice_cpl_dt

    call NUOPC_CompAttributeGet(driver, name="ocn_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ocn_cpl_dt

    call NUOPC_CompAttributeGet(driver, name="glc_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) glc_cpl_dt

    call NUOPC_CompAttributeGet(driver, name="rof_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) rof_cpl_dt

    call NUOPC_CompAttributeGet(driver, name="wav_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) wav_cpl_dt

    ! TODO: for now - this is not in the namelist_definition_drv.xml file
    ! call NUOPC_CompAttributeGet(driver, name="esp_cpl_dt", value=cvalue, rc=rc)
    ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) call shr_sys_abort()
    ! read(cvalue,*) esp_cpl_dt
    esp_cpl_dt = 0.

    call NUOPC_CompAttributeGet(driver, name="glc_avg_period", value=glc_avg_period, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) glc_avg_period

    call NUOPC_CompAttributeGet(driver, name="esp_run_on_pause", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) esp_run_on_pause

    call NUOPC_CompAttributeGet(driver, name="end_restart", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) end_restart

    if (read_restart) then
       if (iam == 0) then
          call NUOPC_CompAttributeGet(driver, name='driver_restart_file', value=restart_file, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          !--- read rpointer if restart_file is set to sp_str ---
          if (trim(restart_file) == trim(sp_str)) then

             ! Error check on restart_pfile
             call NUOPC_CompAttributeGet(driver, name="driver_restart_pfile", value=restart_pfile, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             if ( len_trim(restart_pfile) == 0 ) then
                rc = ESMF_FAILURE
                call ESMF_LogWrite(trim(subname)//' ERROR driver_restart_pfile must be defined', &
                     ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
                return
             end if

             unitn = shr_file_getUnit()
             call ESMF_LogWrite(trim(subname)//" read rpointer file = "//trim(restart_pfile), &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             open(unitn, file=restart_pfile, form='FORMATTED', status='old',iostat=ierr)
             if (ierr < 0) then
                rc = ESMF_FAILURE
                call ESMF_LogWrite(trim(subname)//' ERROR rpointer file open returns error', &
                     ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
                return
             end if
             read(unitn,'(a)', iostat=ierr) restart_file
             if (ierr < 0) then
                rc = ESMF_FAILURE
                call ESMF_LogWrite(trim(subname)//' ERROR rpointer file read returns error', &
                     ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
                return
             end if
             close(unitn)
             call shr_file_freeUnit( unitn )
             call ESMF_LogWrite(trim(subname)//" read driver restart from file = "//trim(restart_file), &
                  ESMF_LOGMSG_INFO, rc=dbrc)
          endif

          ! tcraig, use netcdf here since it's serial and pio may not have been initialized yet
          status = nf90_open(restart_file, NF90_NOWRITE, ncid)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_open')
          status = nf90_inq_varid(ncid, 'start_ymd', varid)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_inq_varid start_ymd')
          status = nf90_get_var(ncid, varid, start_ymd)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_get_var start_ymd')
          status = nf90_inq_varid(ncid, 'start_tod', varid)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_inq_varid start_tod')
          status = nf90_get_var(ncid, varid, start_tod)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_get_var start_tod')
          status = nf90_inq_varid(ncid, 'ref_ymd', varid)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_inq_varid ref_ymd')
          status = nf90_get_var(ncid, varid, ref_ymd)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_get_var ref_ymd')
          status = nf90_inq_varid(ncid, 'ref_tod', varid)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_inq_varid ref_tod')
          status = nf90_get_var(ncid, varid, ref_tod)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_get_var ref_tod')
          status = nf90_inq_varid(ncid, 'curr_ymd', varid)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_inq_varid curr_ymd')
          status = nf90_get_var(ncid, varid, curr_ymd)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_get_var curr_ymd')
          status = nf90_inq_varid(ncid, 'curr_tod', varid)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_inq_varid curr_tod')
          status = nf90_get_var(ncid, varid, curr_tod)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_get_var curr_tod')
          status = nf90_close(ncid)
          if (status /= nf90_NoErr) call shr_sys_abort(trim(subname)//' ERROR: nf90_close')

          write(tmpstr,*) trim(subname)//" read start_ymd = ",start_ymd
          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
          write(tmpstr,*) trim(subname)//" read start_tod = ",start_tod
          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
          write(tmpstr,*) trim(subname)//" read ref_ymd   = ",ref_ymd
          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
          write(tmpstr,*) trim(subname)//" read ref_tod   = ",ref_tod
          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
          write(tmpstr,*) trim(subname)//" read curr_ymd  = ",curr_ymd
          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
          write(tmpstr,*) trim(subname)//" read curr_tod  = ",curr_tod
          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

       endif

       call shr_mpi_bcast(start_ymd, mpicom)
       call shr_mpi_bcast(start_tod, mpicom)
       call shr_mpi_bcast(  ref_ymd, mpicom)
       call shr_mpi_bcast(  ref_tod, mpicom)
       call shr_mpi_bcast( curr_ymd, mpicom)
       call shr_mpi_bcast( curr_tod, mpicom)

    endif

    !---------------------------------------------------------------------------
    ! Modify input config data as needed
    !---------------------------------------------------------------------------

    if (lnd_cpl_dt == 0) lnd_cpl_dt = atm_cpl_dt ! Copy atm coupling time into lnd
    if (rof_cpl_dt == 0) rof_cpl_dt = atm_cpl_dt ! Copy atm coupling time into rof
    if (ice_cpl_dt == 0) ice_cpl_dt = atm_cpl_dt ! Copy atm coupling time into ice
    if (ocn_cpl_dt == 0) ocn_cpl_dt = atm_cpl_dt ! Copy atm coupling time into ocn
    if (glc_cpl_dt == 0) glc_cpl_dt = atm_cpl_dt ! Copy atm coupling time into glc
    if (wav_cpl_dt == 0) wav_cpl_dt = atm_cpl_dt ! Copy atm coupling time into wav
    if (esp_cpl_dt == 0) esp_cpl_dt = atm_cpl_dt ! Copy atm coupling time into esp

    if ( ref_ymd == 0 ) then
       ref_ymd = start_ymd
       ref_tod = start_tod
    endif
    if ( curr_ymd == 0 ) then
       curr_ymd = start_ymd
       curr_tod = start_tod
    endif
    if ( stop_ymd < 0) then
       stop_ymd = 99990101
       stop_tod = 0
    endif
    if (trim(restart_option) == trim(seq_timemgr_optNone) .or. &
         trim(restart_option) == trim(seq_timemgr_optNever)) then
       if (end_restart) then
          end_restart = .false.
          write(logunit,F0A) trim(subname),' WARNING: overriding end_restart to '// &
               'false based on restart_option '
       endif
    endif
    if (trim(restart_option) == trim(seq_timemgr_optEnd)) then
       restart_option = seq_timemgr_optNone
       write(logunit,F0A) trim(subname),' WARNING: overriding restart_option to '// &
            'none and verifying end_restart flag is true '
       if (.not. end_restart) then
          end_restart = .true.
          write(logunit,F0A) trim(subname),' WARNING: overriding end_restart to '// &
               'true based on restart_option (end) '
       endif
    endif

    !---------------------------------------------------------------------------
    ! Print out the namelist settings
    !---------------------------------------------------------------------------

    if (mastertask) then
       write(logunit,F0A) ' '
       write(logunit,F0A) trim(subname),' Clock Init Settings:'
       write(logunit,F0A) trim(subname),' calendar       = ',trim(calendar)
       write(logunit,F0A) trim(subname),' stop_option    = ',trim(stop_option)
       write(logunit,F0I) trim(subname),' stop_n         = ',stop_n
       write(logunit,F0I) trim(subname),' stop_ymd       = ',stop_ymd
       write(logunit,F0I) trim(subname),' stop_tod       = ',stop_tod
       write(logunit,F0A) trim(subname),' restart_option = ',trim(restart_option)
       write(logunit,F0I) trim(subname),' restart_n      = ',restart_n
       write(logunit,F0I) trim(subname),' restart_ymd    = ',restart_ymd
       write(logunit,F0L) trim(subname),' end_restart    = ',end_restart
       write(logunit,F0A) trim(subname),' pause_option   = ',trim(pause_option)
       write(logunit,F0I) trim(subname),' pause_n        = ',pause_n
       write(logunit,F0A) trim(subname),' pause_component_list = ',trim(pause_component_list)
       write(logunit,F0L) trim(subname),' esp_run_on_pause = ',esp_run_on_pause
       write(logunit,F0A) trim(subname),' history_option = ',trim(history_option)
       write(logunit,F0I) trim(subname),' history_n      = ',history_n
       write(logunit,F0I) trim(subname),' history_ymd    = ',history_ymd
       write(logunit,F0A) trim(subname),' histavg_option = ',trim(histavg_option)
       write(logunit,F0I) trim(subname),' histavg_n      = ',histavg_n
       write(logunit,F0I) trim(subname),' histavg_ymd    = ',histavg_ymd
       write(logunit,F0A) trim(subname),' barrier_option = ',trim(barrier_option)
       write(logunit,F0I) trim(subname),' barrier_n      = ',barrier_n
       write(logunit,F0I) trim(subname),' barrier_ymd    = ',barrier_ymd
       write(logunit,F0A) trim(subname),' tprof_option   = ',trim(tprof_option)
       write(logunit,F0I) trim(subname),' tprof_n        = ',tprof_n
       write(logunit,F0I) trim(subname),' tprof_ymd      = ',tprof_ymd
       write(logunit,F0I) trim(subname),' start_ymd      = ',start_ymd
       write(logunit,F0I) trim(subname),' start_tod      = ',start_tod
       write(logunit,F0I) trim(subname),' ref_ymd        = ',ref_ymd
       write(logunit,F0I) trim(subname),' ref_tod        = ',ref_tod
       write(logunit,F0I) trim(subname),' atm_cpl_dt     = ',atm_cpl_dt
       write(logunit,F0I) trim(subname),' lnd_cpl_dt     = ',lnd_cpl_dt
       write(logunit,F0I) trim(subname),' ice_cpl_dt     = ',ice_cpl_dt
       write(logunit,F0I) trim(subname),' ocn_cpl_dt     = ',ocn_cpl_dt
       write(logunit,F0I) trim(subname),' glc_cpl_dt     = ',glc_cpl_dt
       write(logunit,F0A) trim(subname),' glc_avg_period = ',glc_avg_period
       write(logunit,F0I) trim(subname),' rof_cpl_dt     = ',rof_cpl_dt
       write(logunit,F0I) trim(subname),' wav_cpl_dt     = ',wav_cpl_dt
       write(logunit,F0I) trim(subname),' esp_cpl_dt     = ',esp_cpl_dt
       write(logunit,F0A) ' '

       ! check couling intervals
       if ( atm_cpl_dt <= 0          .or. &
            lnd_cpl_dt /= atm_cpl_dt .or. &
            ice_cpl_dt /= atm_cpl_dt .or. &
            ocn_cpl_dt <= 0          .or. &
            glc_cpl_dt <= 0          .or. &
            rof_cpl_dt <=0           .or. &
            wav_cpl_dt <=0           .or. &
            esp_cpl_dt <=0) then

          write(logunit,*) trim(subname),' ERROR: aliogrwe _cpl_dt = ', &
               atm_cpl_dt, lnd_cpl_dt, ice_cpl_dt, ocn_cpl_dt, glc_cpl_dt, &
               rof_cpl_dt, wav_cpl_dt, esp_cpl_dt

          call shr_sys_abort( subname//': ERROR coupling intervals invalid' )
       end if

       ! check start time date
       if ( (start_ymd < 101) .or. (start_ymd > 99991231)) then
          write(logunit,*) subname,' ERROR: illegal start_ymd',start_ymd
          call shr_sys_abort( subname//': ERROR invalid start_ymd')
       end if

    endif

    ! set module variable seq_timemgr_histavg_type
    if     (trim(histavg_option) == trim(seq_timemgr_optNever) .or. &
            trim(histavg_option) == trim(seq_timemgr_optNone)) then

       seq_timemgr_histavg_type = seq_timemgr_type_never

    elseif (trim(histavg_option) == trim(seq_timemgr_optNHours) .or. &
            trim(histavg_option) == trim(seq_timemgr_optNHour)) then

       seq_timemgr_histavg_type = seq_timemgr_type_nhour

    elseif (trim(histavg_option) == trim(seq_timemgr_optNDays) .or. &
            trim(histavg_option) == trim(seq_timemgr_optNDay)) then

       seq_timemgr_histavg_type = seq_timemgr_type_nday

    elseif (trim(histavg_option) == trim(seq_timemgr_optNMonths) .or. &
            trim(histavg_option) == trim(seq_timemgr_optNMonth)  .or. &
            trim(histavg_option) == trim(seq_timemgr_optMonthly)) then

       seq_timemgr_histavg_type = seq_timemgr_type_nmonth

    elseif (trim(histavg_option) == trim(seq_timemgr_optNYears) .or. &
            trim(histavg_option) == trim(seq_timemgr_optNYear)  .or. &
            trim(histavg_option) == trim(seq_timemgr_optYearly)) then

       seq_timemgr_histavg_type = seq_timemgr_type_nyear

    else

       seq_timemgr_histavg_type = seq_timemgr_type_other

    endif

    ! --- Initialize generic stuff ---
    seq_timemgr_calendar             = shr_cal_calendarName(calendar)
    seq_timemgr_esp_run_on_pause     = esp_run_on_pause
    seq_timemgr_end_restart          = end_restart

    ! --- Figure out which components (if any) are doing pause this run
    rc = 1
    i = 1
    if (trim(pause_component_list) == 'all') then
       pause_active = .true.
    else if (trim(pause_component_list) == 'none') then
       pause_active = .false.
    else
       do
          i = scan(trim(pause_component_list(rc:)), ':') - 1
          if ((i < 0) .and. (len_trim(pause_component_list) >= rc)) then
             i = len_trim(pause_component_list(rc:))
          end if
          if (i > 0) then
             found = .false.
             do n = 1, max_clocks
                if (pause_component_list(rc:rc+i-1) == trim(seq_timemgr_clocks(n))) then
                   pause_active(n) = .true.
                   found = .true.
                   exit
                end if
             end do
             ! Special case for cpl -- synonym for drv
             if ((.not. found) .and. (pause_component_list(rc:rc+i-1) == 'cpl')) then
                pause_active(seq_timemgr_nclock_drv) = .true.
                found = .true.
             end if
             if (.not. found) then
                call shr_sys_abort(subname//': unknown pause component, '//pause_component_list(rc:rc+i-1))
             end if
             rc = rc + i
             if (pause_component_list(rc:rc) == ':') then
                rc = rc + 1
             end if
             if (rc >= len_trim(pause_component_list)) then
                exit
             end if
          else
             exit
          end if
       end do
    end if
    if ( ANY(pause_active) .and.                                              &
         (trim(pause_option) /= seq_timemgr_optNONE)  .and.                   &
         (trim(pause_option) /= seq_timemgr_optNever)) then
       do n = 1, max_clocks
          if (pause_active(n)) then
             write(logunit, '(4a)') subname, ': Pause active for ',           &
                  trim(seq_timemgr_clocks(n)),' component'
          end if
       end do
    end if

    ! --- Create the new calendar if not already set ------
    if ( trim(seq_timemgr_calendar) == trim(seq_timemgr_noleap)) then
       esmf_caltype = ESMF_CALKIND_NOLEAP
    else if ( trim(seq_timemgr_calendar) == trim(seq_timemgr_gregorian)) then
       esmf_caltype = ESMF_CALKIND_GREGORIAN
    else
       write(logunit,*) subname//': unrecognized ESMF calendar specified: '// &
            trim(seq_timemgr_calendar)
       call shr_sys_abort( subname//'ERROR:: bad calendar for ESMF' )
    end if

    seq_timemgr_cal = ESMF_CalendarCreate( name='CMEPS_'//seq_timemgr_calendar, calkindflag=esmf_caltype, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! --- Initialize start, ref, and current date ---

    call seq_timemgr_ETimeInit( StartTime, start_ymd, start_tod, "Start date" )
    call seq_timemgr_ETimeInit( RefTime  , ref_ymd  , ref_tod  , "Reference date" )
    call seq_timemgr_ETimeInit( CurrTime , curr_ymd , curr_tod , "Current date")

    ! --- Figure out what time-stepping interval should be. ---------------

    dtime = 0
    dtime(seq_timemgr_nclock_atm     ) = atm_cpl_dt
    dtime(seq_timemgr_nclock_lnd     ) = lnd_cpl_dt
    dtime(seq_timemgr_nclock_ocn     ) = ocn_cpl_dt
    dtime(seq_timemgr_nclock_ice     ) = ice_cpl_dt
    dtime(seq_timemgr_nclock_glc     ) = glc_cpl_dt
    dtime(seq_timemgr_nclock_rof     ) = rof_cpl_dt
    dtime(seq_timemgr_nclock_wav     ) = wav_cpl_dt
    dtime(seq_timemgr_nclock_esp     ) = esp_cpl_dt

    ! --- this finds the min of dtime excluding the driver value ---
    dtime(seq_timemgr_nclock_drv) = maxval(dtime)
    dtime(seq_timemgr_nclock_drv) = minval(dtime)

    ! --- For figuring pause cycle
    min_dt = maxval(dtime)
    seq_timemgr_pause_sig_index = -1

    do n = 1,max_clocks
       if ( mod(dtime(n),dtime(seq_timemgr_nclock_drv)) /= 0) then
          write(logunit,*) trim(subname),' ERROR: dtime inconsistent = ',dtime
          call shr_sys_abort( subname//' :coupling intervals not compatible' )
       endif
       if (pause_active(n) .and. (dtime(n) < min_dt)) then
          min_dt = dtime(n)
          seq_timemgr_pause_sig_index = n
       end if
    enddo
    if (ANY(pause_active)) then
       if (seq_timemgr_pause_sig_index < 1) then
          write(logunit, *) subname,"ERROR: No pause_sig_index even with active pause"
          call shr_sys_abort(subname//"ERROR: No pause_sig_index even with active pause")
       end if
    else
       ! Don't try to run ESP on non-existent pauses
       seq_timemgr_esp_run_on_pause = .false.
    end if

    ! --- Initialize component and driver clocks and alarms common to components and driver clocks ---
    SyncClock%ECP(seq_timemgr_nclock_drv)%EClock => EClock_drv
    SyncClock%ECP(seq_timemgr_nclock_atm)%EClock => EClock_atm
    SyncClock%ECP(seq_timemgr_nclock_lnd)%EClock => EClock_lnd
    SyncClock%ECP(seq_timemgr_nclock_ocn)%EClock => EClock_ocn
    SyncClock%ECP(seq_timemgr_nclock_ice)%EClock => EClock_ice
    SyncClock%ECP(seq_timemgr_nclock_glc)%EClock => EClock_glc
    SyncClock%ECP(seq_timemgr_nclock_rof)%EClock => EClock_rof
    SyncClock%ECP(seq_timemgr_nclock_wav)%EClock => EClock_wav
    SyncClock%ECP(seq_timemgr_nclock_esp)%EClock => EClock_esp

    do n = 1,max_clocks
       call ESMF_TimeIntervalSet( TimeStep, s=dtime(n), rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call seq_timemgr_EClockInit( TimeStep, StartTime, RefTime, CurrTime, SyncClock%ECP(n)%EClock)

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
            EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_stop),  &
            option  = stop_option,         &
            opt_n   = stop_n,              &
            opt_ymd = stop_ymd,            &
            opt_tod = stop_tod,            &
            RefTime = CurrTime,            &
            alarmname = trim(seq_timemgr_alarm_stop), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
            EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_datestop),  &
            option  = seq_timemgr_optDate, &
            opt_ymd = stop_ymd,            &
            opt_tod = stop_tod,            &
            RefTime = StartTime,           &
            alarmname = trim(seq_timemgr_alarm_datestop), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
            EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_restart),  &
            option  = restart_option,      &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = CurrTime,            &
            alarmname = trim(seq_timemgr_alarm_restart), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
            EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_history),  &
            option  = history_option,      &
            opt_n   = history_n,           &
            opt_ymd = history_ymd,         &
            RefTime = StartTime,           &
            alarmname = trim(seq_timemgr_alarm_history), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
            EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_histavg),  &
            option  = histavg_option,      &
            opt_n   = histavg_n,           &
            opt_ymd = histavg_ymd,         &
            RefTime = StartTime,           &
            alarmname = trim(seq_timemgr_alarm_histavg), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
            EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_barrier),  &
            option  = barrier_option,      &
            opt_n   = barrier_n,           &
            opt_ymd = barrier_ymd,         &
            RefTime = CurrTime,            &
            alarmname = trim(seq_timemgr_alarm_barrier), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
            EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_tprof),  &
            option  = tprof_option,        &
            opt_n   = tprof_n,             &
            opt_ymd = tprof_ymd,           &
            RefTime = StartTime,           &
            alarmname = trim(seq_timemgr_alarm_tprof), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmGet(SyncClock%EAlarm(n,seq_timemgr_nalarm_stop), RingTime=StopTime1, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AlarmGet(SyncClock%EAlarm(n,seq_timemgr_nalarm_datestop), RingTime=StopTime2, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (StopTime2 < StopTime1) then
          call ESMF_ClockSet(SyncClock%ECP(n)%EClock, StopTime=StopTime2, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_ClockSet(SyncClock%ECP(n)%EClock, StopTime=StopTime1, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

       ! Set the pause option if pause/resume is active
       if (pause_active(n)) then
          call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock,                  &
               EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_pause),         &
               option  = pause_option,                                         &
               opt_n   = pause_n,                                              &
               RefTime = CurrTime,                                             &
               alarmname = trim(seq_timemgr_alarm_pause), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock,                  &
               EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_pause),         &
               option  = seq_timemgr_optNever,                                 &
               opt_n   = -1,                                                   &
               RefTime = StartTime,                                            &
               alarmname = trim(seq_timemgr_alarm_pause), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

    enddo

    if (mastertask) then
       call seq_timemgr_clockPrint(SyncClock)
    endif

  end subroutine seq_timemgr_clockInit

  !===============================================================================

  subroutine seq_timemgr_EClockGetData( EClock, &
       curr_yr, curr_mon, curr_day,    &
       curr_ymd, curr_tod, prev_ymd, prev_tod, start_ymd,     &
       start_tod, StepNo, ref_ymd, ref_tod,         &
       stop_ymd, stop_tod, dtime, ECurrTime, alarmcount,      &
       curr_cday, next_cday, curr_time, prev_time, calendar)

    ! !DESCRIPTION: Get various values from the clock.
    use ESMF, only: ESMF_Clock, ESMF_Time, ESMF_TimeInterval
    use ESMF, only: ESMF_ClockGet, ESMF_TimeGet, ESMF_TimeIntervalGet
    use ESMF, only: ESMF_TimeSet, ESMF_TimeIntervalSet
    use med_constants_mod, only : IN, R8, I8
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_ChkErr

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)     , intent(in)            :: EClock     ! Input clock object
    integer(IN) , intent(out), optional :: curr_yr    ! Current year
    integer(IN) , intent(out), optional :: curr_mon   ! Current month
    integer(IN) , intent(out), optional :: curr_day   ! Current day in month
    integer(IN) , intent(out), optional :: curr_ymd   ! Current date YYYYMMDD
    integer(IN) , intent(out), optional :: curr_tod   ! Current time of day (s)
    integer(IN) , intent(out), optional :: prev_ymd   ! Previous date YYYYMMDD
    integer(IN) , intent(out), optional :: prev_tod   ! Previous time of day (s)
    integer(IN) , intent(out), optional :: start_ymd  ! Starting date YYYYMMDD
    integer(IN) , intent(out), optional :: start_tod  ! Starting time-of-day (s)
    integer(IN) , intent(out), optional :: StepNo     ! Number of steps taken
    integer(IN) , intent(out), optional :: ref_ymd    ! Reference date YYYYMMDD
    integer(IN) , intent(out), optional :: ref_tod    ! Reference time-of-day (s)
    integer(IN) , intent(out), optional :: stop_ymd   ! Stop date YYYYMMDD
    integer(IN) , intent(out), optional :: stop_tod   ! Stop time-of-day (s)
    integer(IN) , intent(out), optional :: dtime      ! Time-step (seconds)
    integer(IN) , intent(out), optional :: alarmcount ! Number of Valid Alarms
    type(ESMF_Time)      , intent(out), optional :: ECurrTime  ! Current ESMF time
    real(R8)   , intent(out), optional :: curr_cday  ! current calendar day
    real(R8)   , intent(out), optional :: next_cday  ! current calendar day
    real(R8)   , intent(out), optional :: curr_time  ! time interval between current time and reference date
    real(R8)   , intent(out), optional :: prev_time  ! time interval between previous time and reference date
    character(len=*)     , intent(out), optional :: calendar   ! calendar type

    !----- local -----
    type(ESMF_Time)         :: CurrentTime     ! Current time
    type(ESMF_Time)         :: PreviousTime    ! Previous time
    type(ESMF_Time)         :: StartTime       ! Start time
    type(ESMF_Time)         :: StopTime        ! Stop time
    type(ESMF_Time)         :: RefTime         ! Ref time
    type(ESMF_TimeInterval) :: timeStep        ! Clock, time-step
    type(ESMF_TimeInterval) :: timediff        ! Used to calculate curr_time
    integer(IN)             :: rc              ! Return code
    integer(I8)             :: advSteps        ! Number of time-steps that have advanced
    integer(IN)             :: yy, mm, dd, sec ! Return time values
    integer(IN)             :: ymd             ! Date (YYYYMMDD)
    integer(IN)             :: tod             ! time of day (sec)
    integer(IN)             :: ldtime          ! local dtime
    integer(IN)             :: days            ! number of whole days in time interval
    integer(IN)             :: seconds         ! number of seconds in time interval
    integer(IN)             :: acount          ! number of valid alarms
    real(R8)                :: doy, tmpdoy     ! day of year
    type(ESMF_Time)         :: tmpTime         ! tmp time, needed for next_cday
    type(ESMF_TimeInterval) :: tmpDTime        ! tmp time interval, needed for next_cday
    real(R8), parameter     :: c1 = 1.0_R8
    character(len=*)  , parameter :: subname = '(seq_timemgr_EClockGetData) '
    !-------------------------------------------------------------------------------

    if (present(calendar)) calendar = trim(seq_timemgr_calendar)

    call ESMF_ClockGet( EClock, currTime=CurrentTime, &
         advanceCount=advSteps, prevTime=previousTime, TimeStep=timeStep, &
         startTime=StartTime, stopTime=stopTime, refTime=RefTime, &
         AlarmCount=acount, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( CurrentTime, yy=yy, mm=mm, dd=dd, s=sec, dayofyear_r8=doy, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call seq_timemgr_ETimeGet( CurrentTime, ymd=ymd, tod=tod )
    call ESMF_TimeIntervalGet( timeStep, s=ldtime, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( present(curr_yr)  ) curr_yr  = yy
    if ( present(curr_mon) ) curr_mon = mm
    if ( present(curr_day) ) curr_day = dd
    if ( present(curr_tod) ) curr_tod = tod
    if ( present(curr_ymd) ) curr_ymd = ymd
    if ( present(ECurrTime)) ECurrTime= CurrentTime
    if ( present(StepNo)   ) StepNo   = advSteps
    if ( present(dtime)    ) dtime    = ldtime
    if ( present(curr_cday)) curr_cday = doy
    if ( present(alarmcount)) alarmcount = acount

    if ( present(next_cday)) then
       call ESMF_TimeSet(tmpTime, yy=yy, mm=mm, dd=dd, s=tod, calendar=seq_timemgr_cal, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeIntervalSet( tmpDTime, d=1, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       tmpTime = tmpTime + tmpDTime
       call ESMF_TimeGet(tmpTime, dayOfYear_r8=tmpdoy, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       next_cday = tmpdoy
    endif

    ! ---Current Time (the time interval between the current date and the reference date) ---
    if ( present(curr_time)) then
       timediff = CurrentTime - RefTime
       call ESMF_TimeIntervalGet(timediff, d=days, s=seconds, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       curr_time = days + seconds/real(SecPerDay,R8)
    end if

    ! ---Previous Time (the time interval between the previous date and the reference date) ---
    if ( present(prev_time)) then
       timediff = PreviousTime - RefTime
       call ESMF_TimeIntervalGet(timediff, d=days, s=seconds, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       prev_time = days + seconds/real(SecPerDay,R8)
    end if

    ! --- Previous time --------------------------------------------------------
    if ( present(prev_ymd) .or. present(prev_tod) )then
       call seq_timemgr_ETimeGet( PreviousTime, ymd=ymd, tod=tod )
       if ( present(prev_ymd) ) prev_ymd = ymd
       if ( present(prev_tod) ) prev_tod = tod
    end if

    ! --- If want start date -----------------------------------------------
    if ( present(start_ymd) .or. present(start_tod) )then
       call seq_timemgr_ETimeGet( StartTime, ymd=ymd, tod=tod )
       if ( present(start_ymd) ) start_ymd = ymd
       if ( present(start_tod) ) start_tod = tod
    end if

    ! --- If want stop date -----------------------------------------------
    if ( present(stop_ymd) .or. present(stop_tod) )then
       call seq_timemgr_ETimeGet( stopTime, ymd=ymd, tod=tod )
       if ( present(stop_ymd) ) stop_ymd = ymd
       if ( present(stop_tod) ) stop_tod = tod
    end if

    ! --- If want ref date -----------------------------------------------
    if ( present(ref_ymd) .or. present(ref_tod) )then
       call seq_timemgr_ETimeGet( RefTime, ymd=ymd, tod=tod )
       if ( present(ref_ymd) ) ref_ymd = ymd
       if ( present(ref_tod) ) ref_tod = tod
    end if

  end subroutine seq_timemgr_EClockGetData

  !===============================================================================

  subroutine seq_timemgr_alarmInit( EClock, EAlarm, option, opt_n, opt_ymd, opt_tod, RefTime, alarmname, rc)

    ! !DESCRIPTION: Setup an alarm in a clock
    use shr_sys_mod, only : shr_sys_abort
    use ESMF, only : ESMF_Clock, ESMF_Alarm, ESMF_ClockGet, ESMF_Time, ESMF_TimeGet
    use ESMF, only : ESMF_TimeIntervalSet, ESMF_TimeSet, ESMF_TimeInterval
    use ESMF, only: ESMF_AlarmCreate
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_ChkErr

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)           , intent(INOUT) :: EClock    ! clock
    type(ESMF_Alarm)           , intent(INOUT) :: EAlarm    ! alarm
    character(len=*)           , intent(in)    :: option    ! alarm option
    integer(IN)      ,optional , intent(in)    :: opt_n     ! alarm freq
    integer(IN)      ,optional , intent(in)    :: opt_ymd   ! alarm ymd
    integer(IN)      ,optional , intent(in)    :: opt_tod   ! alarm tod (sec)
    type(ESMF_Time)  ,optional , intent(in)    :: RefTime   ! ref time
    character(len=*) ,optional , intent(in)    :: alarmname ! alarm name
    integer                    , intent(INOUT) :: rc        ! Return code

    !----- local -----
    integer                 :: lymd             ! local ymd
    integer                 :: ltod             ! local tod
    integer                 :: cyy,cmm,cdd,csec ! time info
    integer                 :: nyy,nmm,ndd,nsec ! time info
    character(len=64)       :: lalarmname       ! local alarm name
    logical                 :: update_nextalarm ! update next alarm
    type(ESMF_Time)         :: CurrTime         ! Current Time
    type(ESMF_Time)         :: NextAlarm        ! Next restart alarm time
    type(ESMF_TimeInterval) :: AlarmInterval    ! Alarm interval
    character(len=*), parameter :: subname = '(seq_timemgr_alarmInit): '
    !-------------------------------------------------------------------------------
    ! Notes: The ringtime sent to AlarmCreate MUST be the next alarm time.
    ! If you send an arbitrary but proper ringtime from the past and the ring interval,
    ! the alarm will always go off on the next clock advance and this will cause serious problems.
    ! Even if it makes sense to initialize an alarm with some reference time and the alarm interval,
    ! that reference time has to be advance forward to be >= the current time.  In the logic below
    ! we set an appropriate "NextAlarm" and then we make sure to advance it properly based on the
    ! ring interval.
    !-------------------------------------------------------------------------------

    lalarmname = 'alarm_unknown'
    if (present(alarmname)) then
       lalarmname = trim(alarmname)
    endif

    ltod = 0
    if (present(opt_tod)) then
       ltod = opt_tod
    endif

    lymd = -1
    if (present(opt_ymd)) then
       lymd = opt_ymd
    endif

    call ESMF_ClockGet(EClock, CurrTime=CurrTime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(CurrTime, yy=cyy, mm=cmm, dd=cdd, s=csec, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(CurrTime, yy=nyy, mm=nmm, dd=ndd, s=nsec, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! --- initial guess of next alarm, this will be updated below ---
    if (present(RefTime)) then
       NextAlarm = RefTime
    else
       NextAlarm = CurrTime
    endif

    update_nextalarm  = .true.

    selectcase (trim(option))

    case (seq_timemgr_optNONE)
       !--- tcx seems we need an alarm interval or the alarm create fails,
       !--- problem in esmf_wrf_timemgr?
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=seq_timemgr_cal, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (seq_timemgr_optNever)
       !--- tcx seems we need an alarm interval or the alarm create fails,
       !--- problem in esmf_wrf_timemgr?
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=seq_timemgr_cal, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (seq_timemgr_optDate)
       !--- tcx seems we need an alarm interval or the alarm create fails,
       !--- problem in esmf_wrf_timemgr?
       if (.not. present(opt_ymd)) call shr_sys_abort(subname//trim(option)//' requires opt_ymd')
       if (lymd < 0 .or. ltod < 0) call shr_sys_abort(subname//trim(option)//'opt_ymd, opt_tod invalid')
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call seq_timemgr_ETimeInit(NextAlarm, lymd, ltod, "optDate")
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (seq_timemgr_optIfdays0)
       if (.not. present(opt_ymd)) call shr_sys_abort(subname//trim(option)//' requires opt_ymd')
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=opt_n, s=0, calendar=seq_timemgr_cal, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case (seq_timemgr_optNSteps)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_ClockGet(EClock, TimeStep=AlarmInterval, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNStep)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_ClockGet(EClock, TimeStep=AlarmInterval, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNSeconds)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNSecond)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNMinutes)
       call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNMinute)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNHours)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNHour)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNDays)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNDay)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNMonths)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNMonth)
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optMonthly)
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=1, s=0, calendar=seq_timemgr_cal, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case (seq_timemgr_optNYears)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNYear)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optYearly)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=1, dd=1, s=0, calendar=seq_timemgr_cal, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    case (seq_timemgr_optEnd)
       call shr_sys_abort(subname//'deprecated option '//trim(option))

    case default
       call shr_sys_abort(subname//'unknown option '//trim(option))

    end select

    ! --------------------------------------------------------------------------------
    ! --- AlarmInterval and NextAlarm should be set ---
    ! --------------------------------------------------------------------------------

    ! --- advance Next Alarm so it won't ring on first timestep for
    ! --- most options above. go back one alarminterval just to be careful

    if (update_nextalarm) then
       NextAlarm = NextAlarm - AlarmInterval
       do while (NextAlarm <= CurrTime)
          NextAlarm = NextAlarm + AlarmInterval
       enddo
    endif

    EAlarm = ESMF_AlarmCreate( name=lalarmname, clock=EClock, ringTime=NextAlarm, ringInterval=AlarmInterval, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine seq_timemgr_AlarmInit

  !===============================================================================

  subroutine seq_timemgr_alarmGet( EAlarm, next_ymd, next_tod, prev_ymd, prev_tod,    &
       IntSec, IntMon, IntYrs, name)

    ! !DESCRIPTION: Get informationn from the alarm
    use med_constants_mod, only : IN
    use ESMF, only: ESMF_Alarm, ESMF_Time, ESMF_TimeInterval, ESMF_AlarmGet, ESMF_TimeIntervalGet
    use ESMF, only: ESMF_ALARMLIST_ALL
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_ChkErr

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Alarm)    , intent(INOUT)            :: EAlarm   ! Input Alarm object
    integer(IN), intent(out), optional :: next_ymd ! alarm date yyyymmdd
    integer(IN), intent(out), optional :: next_tod ! alarm tod sec
    integer(IN), intent(out), optional :: prev_ymd ! alarm date yyyymmdd
    integer(IN), intent(out), optional :: prev_tod ! alarm tod sec
    integer(IN), intent(out), optional :: IntSec   ! alarm int sec
    integer(IN), intent(out), optional :: IntMon   ! alarm int mon
    integer(IN), intent(out), optional :: IntYrs   ! alarm int yrs
    character(len=*)    , intent(out), optional :: name     ! alarm name

    !----- local -----
    integer                 :: yy, mm, dd, sec ! Return time values
    integer                 :: ymd             ! Date (YYYYMMDD)
    integer                 :: tod             ! time of day (sec)
    integer                 :: rc              ! error code
    type(ESMF_TimeInterval) :: alarmInterval   ! Alarm interval
    type(ESMF_Time)         :: ringTime        ! Next alarm ring time
    character(len=*), parameter :: subname = '(seq_timemgr_alarmGet) '
    !-------------------------------------------------------------------------------

    if (present(name)) then
       call ESMF_AlarmGet( EAlarm, name=name, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call ESMF_AlarmGet( EAlarm, RingTime=RingTime, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call seq_timemgr_ETimeGet( RingTime, ymd=ymd, tod=tod)
    if ( present(next_ymd) ) next_ymd = ymd
    if ( present(next_tod) ) next_tod = tod

    call ESMF_AlarmGet( EAlarm, PrevRingTime=RingTime, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call seq_timemgr_ETimeGet( RingTime, ymd=ymd, tod=tod)
    if ( present(prev_ymd) ) prev_ymd = ymd
    if ( present(prev_tod) ) prev_tod = tod

    yy = 0
    mm = 0
    dd = 0
    sec = 0
    call ESMF_AlarmGet( EAlarm, RingInterval=AlarmInterval, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet( alarmInterval, yy=yy, mm=mm, d=dd, s=sec, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    sec = sec + dd*(SecPerDay)

    ! --- If want restart next interval information -------------------------
    if ( present(IntSec) ) IntSec = sec
    if ( present(IntMon) ) IntMon = mm
    if ( present(IntYrs) ) IntYrs = yy

  end subroutine seq_timemgr_alarmGet

  !===============================================================================

  subroutine seq_timemgr_AlarmSetOn( EClock, alarmname)

    ! !DESCRIPTION: turn alarm on
    use shr_sys_mod, only : shr_sys_abort
    use ESMF, only : ESMF_Alarm, ESMF_Clock, ESMF_AlarmRingerOn
    use ESMF, only : ESMF_AlarmGet, ESMF_ClockGetAlarmList
    use ESMF, only : ESMF_ALARMLIST_ALL
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_ChkErr
    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock), intent(INOUT) :: EClock      ! clock/alarm
    character(len=*), intent(in), optional :: alarmname  ! alarmname

    !----- local -----
    integer                     :: n
    integer                     :: rc
    logical                     :: found
    logical                     :: set
    character(len=64)           :: name
    type(ESMF_Alarm),pointer    :: EAlarm_list(:)
    integer(IN)        :: AlarmCount ! Number of valid alarms
    character(len=*), parameter :: xalarm = 'unset'
    character(len=*), parameter :: subname = '(seq_timemgr_alarmSetOn) '

    !-------------------------------------------------------------------------------
    ! Notes: The Alarm_list is returned and only a subset of the alarms may
    !   be initialized.  In the esmf_wrf_timemgr, numalarms is not used internally,
    !   and the alarm pointer is valid if it's associated.  If it's not associated
    !   the AlarmGet calls will generally return an error code.  What we really
    !   want is to ignore the unset alarms.  So below, we have to kind of kludge
    !   this up.  We set name=xalarm, a special value, before the AlarmGet call so
    !   if Alarm_list(n) is not associated, the name will remain the value of
    !   xalarm.  Then we check whether it's a valid alarm by first checking
    !   the name vs xalarm.  If name is not xalarm, then it must be a valid alarm
    !   and we either set found to true if we are setting all alarms or we compare
    !   the name returned to the alarm name we're looking for and only set found
    !   to true if the names match.
    !-------------------------------------------------------------------------------

    set = .false.

    call seq_timemgr_EClockGetData(EClock, AlarmCount=AlarmCount)
    allocate(EAlarm_list(AlarmCount))
    call ESMF_ClockGetAlarmList(EClock, alarmListFlag=ESMF_ALARMLIST_ALL, &
         alarmList=EAlarm_list, alarmCount=AlarmCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,AlarmCount
       found = .false.
       if (present(alarmname)) then
          call ESMF_AlarmGet(EAlarm_list(n), name=name, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          if (trim(name) == trim(alarmname)) found = .true.
       else
          found = .true.
       endif
       if (found) then
          set = .true.
          call ESMF_AlarmRingerOn( EAlarm_list(n), rc=rc )
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
    enddo

    if (present(alarmname) .and. .not. set) then
       call shr_sys_abort(subname//' ERROR in alarmname '//trim(alarmname))
    endif
    deallocate(EAlarm_list)

  end subroutine seq_timemgr_AlarmSetOn

  !===============================================================================

  subroutine seq_timemgr_AlarmSetOff( EClock, alarmname, rc)

    ! !DESCRIPTION: turn alarm off
    use med_constants_mod, only : IN
    use shr_sys_mod, only : shr_sys_abort
    use ESMF, only : ESMF_Clock, ESMF_Alarm, ESMF_AlarmRingerOff
    use ESMF, only : ESMF_ClockGetAlarmList, ESMF_AlarmGet
    use ESMF, only : ESMF_ALARMLIST_ALL
    use seq_comm_mct, only: logunit
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_ChkErr
    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(INOUT) :: EClock      ! clock/alarm
    character(len=*), intent(in), optional :: alarmname  ! alarmname
    integer         , intent(INOUT) :: rc

    !----- local -----
    integer                  :: n
    logical                  :: found
    logical                  :: set
    character(len=64)        :: name
    type(ESMF_Alarm),pointer :: EAlarm_list(:)
    integer(IN)     :: AlarmCount ! Number of valid alarms
    character(len=*), parameter :: xalarm = 'unset'
    character(len=*), parameter :: subname = '(seq_timemgr_alarmSetOff) '

    !-------------------------------------------------------------------------------
    ! Notes: The Alarm_list is returned and only a subset of the alarms may
    !   be initialized. We check whether it's a valid alarm by first checking
    !   the name vs xalarm.  If name is not xalarm, then it must be a valid alarm
    !   and we either set found to true if we are setting all alarms or we compare
    !   the name returned to the alarm name we're looking for and only set found
    !   to true if the names match.
    !-------------------------------------------------------------------------------

    set = .false.

    call seq_timemgr_EClockGetData(EClock, AlarmCount=AlarmCount)
    allocate(EAlarm_list(AlarmCount))
    call ESMF_ClockGetAlarmList(EClock, alarmListFlag=ESMF_ALARMLIST_ALL, &
         alarmList=EAlarm_list, alarmCount=AlarmCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,AlarmCount
       found = .false.
       if (present(alarmname)) then
          call ESMF_AlarmGet(EAlarm_list(n), name=name, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          if (trim(name) == trim(alarmname)) found = .true.
       else
          found = .true.
       endif
       if (found) then
          set = .true.
          call ESMF_AlarmRingerOff( EAlarm_list(n), rc=rc )
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
    enddo

    if (present(alarmname) .and. .not. set) then
       write(logunit,*) subname,' ERROR in alarmname ',trim(alarmname)
       call shr_sys_abort()
    endif
    deallocate(EAlarm_list)

  end subroutine seq_timemgr_AlarmSetOff

  !===============================================================================

  logical function seq_timemgr_alarmIsOn( EClock, alarmname, rc)

    ! !DESCRIPTION: check if an alarm is ringing
    use shr_sys_mod, only : shr_sys_abort
    use ESMF, only : ESMF_Clock, ESMF_Time, ESMF_Alarm, ESMF_AlarmIsRinging
    use ESMF, only : ESMF_ClockGetAlarmList, ESMF_AlarmGet, ESMF_ClockGet
    use ESMF, only : ESMF_ALARMLIST_ALL
    use seq_comm_mct, only : logunit
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_ChkErr

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock), intent(in)    :: EClock     ! clock/alarm
    character(len=*), intent(in)    :: alarmname  ! which alarm
    integer         , intent(INOUT) :: rc         ! return code

    !----- local -----
    integer                  :: n
    logical                  :: found
    character(len=64)        :: name
    type(ESMF_Time)          :: ETime1, ETime2
    type(ESMF_Alarm),pointer :: EAlarm_list(:)
    integer(IN)     :: AlarmCount ! Number of valid alarms
    character(len=*), parameter :: xalarm = 'unset'
    character(len=*), parameter :: subname = '(seq_timemgr_alarmIsOn) '

    !-------------------------------------------------------------------------------
    ! Notes:  Because of the esmf_wrf_timemgr implementation with regards to
    ! valid alarms in the alarm_list, we initialize name to xalarm before
    ! querying the alarm name, and if the alarm is not valid, name will not
    ! be updated and we can tell that the alarm is not valid and we should
    ! just ignore it.
    !-------------------------------------------------------------------------------

    seq_timemgr_alarmIsOn = .false.
    found = .false.

    call seq_timemgr_EClockGetData(EClock, AlarmCount=AlarmCount)
    allocate(EAlarm_list(AlarmCount))

    call ESMF_ClockGetAlarmList(EClock, alarmListFlag=ESMF_ALARMLIST_ALL, &
         alarmList=EAlarm_list, alarmCount=AlarmCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1,AlarmCount
       name = trim(xalarm)
       call ESMF_AlarmGet(EAlarm_list(n), name=name, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (trim(name) == trim(alarmname)) then
          found = .true.

          seq_timemgr_alarmIsOn = ESMF_AlarmIsRinging(alarm=EAlarm_list(n),rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          ! --- make sure the datestop will always stop with dates >= stop_date
          if (trim(alarmname) == trim(seq_timemgr_alarm_datestop)) then
             call ESMF_ClockGet(EClock, CurrTime = ETime1, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             call ESMF_AlarmGet(EAlarm_list(n), RingTime = ETime2, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             if (ETime1 >= ETime2) seq_timemgr_alarmIsOn = .true.
          endif

       endif
    enddo

    if (.not.found) then
       write(logunit,*) subname//': ERROR alarm not valid for EClock '//trim(alarmname)
       call shr_sys_abort( subname//'ERROR: alarm invalid '//trim(alarmname) )
    endif
    deallocate(EAlarm_list)

  end function seq_timemgr_alarmIsOn

  !===============================================================================
  logical function seq_timemgr_restartAlarmIsOn( EClock)

    ! !DESCRIPTION: check if restart alarm is ringing
    use ESMF, only : ESMF_Clock
    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock) , intent(in) :: EClock     ! clock/alarm

    !----- local -----
    integer :: rc                    ! return code
    character(len=*), parameter :: subname = '(seq_timemgr_restartAlarmIsOn) '
    !-------------------------------------------------------------------------------

    seq_timemgr_restartAlarmIsOn = seq_timemgr_alarmIsOn(EClock, alarmname=seq_timemgr_alarm_restart, rc=rc)

  end function seq_timemgr_restartAlarmIsOn

  !===============================================================================
  logical function seq_timemgr_stopAlarmIsOn( EClock)

    ! !DESCRIPTION: check if stop alarm is ringing
    use ESMF, only : ESMF_Clock

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock) , intent(in) :: EClock     ! clock/alarm

    !----- local -----
    integer :: rc                    ! return code
    character(len=*), parameter :: subname = '(seq_timemgr_stopAlarmIsOn) '
    !-------------------------------------------------------------------------------

    seq_timemgr_stopAlarmIsOn =  seq_timemgr_alarmIsOn(EClock, alarmname=seq_timemgr_alarm_stop, rc=rc)

  end function seq_timemgr_stopAlarmIsOn

  !===============================================================================
  logical function seq_timemgr_historyAlarmIsOn( EClock)

    ! !DESCRIPTION: check if history alarm is ringing
    use ESMF, only : ESMF_Clock

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock) , intent(in) :: EClock     ! clock/alarm

    !----- local -----
    integer :: rc                    ! return code
    character(len=*), parameter :: subname = '(seq_timemgr_historyAlarmIsOn) '
    !-------------------------------------------------------------------------------

    seq_timemgr_historyAlarmIsOn = seq_timemgr_alarmIsOn(EClock, alarmname=seq_timemgr_alarm_history, rc=rc)

  end function seq_timemgr_historyAlarmIsOn

  !===============================================================================
  logical function seq_timemgr_pauseAlarmIsOn( EClock)

    ! !DESCRIPTION: check if pause alarm is ringing
    use ESMF, only : ESMF_Clock

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock) , intent(in) :: EClock     ! clock/alarm

    !----- local -----
    integer :: rc                    ! return code
    character(len=*), parameter :: subname = '(seq_timemgr_pauseAlarmIsOn) '
    !-------------------------------------------------------------------------------

    seq_timemgr_pauseAlarmIsOn = seq_timemgr_alarmIsOn(EClock, alarmname=seq_timemgr_alarm_pause, rc=rc)

  end function seq_timemgr_pauseAlarmIsOn

  !===============================================================================
  logical function seq_timemgr_pause_active()

    ! !DESCRIPTION: Return .true. if any component is configured for pause/resume

    seq_timemgr_pause_active = ANY(pause_active)

  end function seq_timemgr_pause_active

  !===============================================================================
  integer function seq_timemgr_pause_component_index(component_name)

    ! !DESCRIPTION: Look up a component's internal index for faster processing
    use shr_sys_mod, only : shr_sys_abort

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*), intent(in) :: component_name

    !----- local -----
    integer                     :: ind
    character(len=*), parameter :: subname = '(seq_timemgr_pause_component_index) '
    !-------------------------------------------------------------------------------

    seq_timemgr_pause_component_index = 0
    do ind = 1, max_clocks
       if (trim(component_name) == trim(seq_timemgr_clocks(ind))) then
          seq_timemgr_pause_component_index = ind
          exit
       end if
    end do
    if (seq_timemgr_pause_component_index < 1) then
       if (trim(component_name) == 'cpl') then
          seq_timemgr_pause_component_index = seq_timemgr_nclock_drv
       end if
    end if
    if (seq_timemgr_pause_component_index < 1) then
       call shr_sys_abort(subname//': No index for component '//trim(component_name))
    end if

  end function seq_timemgr_pause_component_index

  !===============================================================================
  logical function seq_timemgr_pause_component_active(component_index)

    ! !DESCRIPTION: Return .true. if component is active in driver pause
    use shr_sys_mod, only : shr_sys_abort

    ! !INPUT/OUTPUT PARAMETERS:
    integer, intent(in) :: component_index

    !----- local -----
    character(len=*), parameter :: subname = '(seq_timemgr_pause_component_active) '
    !-------------------------------------------------------------------------------

    if ((component_index < 1) .or. (component_index > max_clocks)) then
       call shr_sys_abort(subname//': component_index out of range')
    end if
    seq_timemgr_pause_component_active = pause_active(component_index)

  end function seq_timemgr_pause_component_active

  !===============================================================================
  subroutine seq_timemgr_ETimeInit( ETime, ymd, tod, desc )

    use shr_sys_mod           , only : shr_sys_abort
    use ESMF                  , only : ESMF_Time, ESMF_TimeSet
    use shr_cal_mod           , only : shr_cal_date2ymd
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use seq_comm_mct          , only : logunit

    ! !DESCRIPTION: Create the ESMF_Time object corresponding to the given input time, given in
    !               YMD (Year Month Day) and TOD (Time-of-day) format.
    !               Set the time by an integer as YYYYMMDD and integer seconds in the day

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Time) , intent(inout) :: ETime     ! Time
    integer         , intent(in)    :: ymd       ! Year, month, day YYYYMMDD
    integer         , intent(in), optional    :: tod       ! Time of day in seconds
    character(len=*), intent(in), optional    :: desc      ! Description of time to set

    !----- local -----
    character(len=*), parameter :: subname = '(seq_timemgr_ETimeInit) '
    integer :: yr, mon, day          ! Year, month, day as integers
    integer :: ltod                  ! local tod
    character(CL)  :: ldesc ! local desc
    integer :: rc                    ! return code
    !-------------------------------------------------------------------------------

    ltod = 0
    if (present(tod)) then
       ltod = tod
    endif

    ldesc = ''
    if (present(desc)) then
       ldesc = desc
    endif

    if ( (ymd < 0) .or. (ltod < 0) .or. (ltod > SecPerDay) )then
       write(logunit,*) subname//': ERROR yymmdd is a negative number or '// &
            'time-of-day out of bounds', ymd, ltod
       call shr_sys_abort( subname//'ERROR: Bad input' )
    end if

    call shr_cal_date2ymd(ymd,yr,mon,day)

    call ESMF_TimeSet( ETime, yy=yr, mm=mon, dd=day, s=ltod, calendar=seq_timemgr_cal, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine seq_timemgr_ETimeInit

  !===============================================================================
  subroutine seq_timemgr_ETimeGet( ETime, offset, ymd, tod )

    ! !DESCRIPTION: Get the date in YYYYMMDD format from a ESMF time object.
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_ChkErr
    use ESMF, only : ESMF_Time, ESMF_TimeInterval, ESMF_TimeIntervalGet, ESMF_TimeIntervalSet
    use ESMF, only : ESMF_TimeGet
    use shr_cal_mod, only : shr_cal_ymd2date
    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Time),   intent(in)  :: ETime   ! Input ESMF time
    integer, optional, intent(in)  :: offset  ! Offset from input time (sec)
    integer, optional, intent(out) :: ymd     ! date of day
    integer, optional, intent(out) :: tod     ! Time of day

    !----- local -----
    character(len=*), parameter :: subname = '(seq_timemgr_ETimeGet) '
    type(ESMF_Time)         :: ETimeAdd ! ESMF time + offset
    type(ESMF_TimeInterval) :: ETimeOff ! ESMF offset time-interval
    integer                 :: year     ! Year
    integer                 :: month    ! Month
    integer                 :: day      ! Day in month
    integer                 :: sec      ! Day in month
    integer                 :: rc       ! Return code
    !-------------------------------------------------------------------------------

    ETimeAdd = ETime
    if ( present(offset) )then
       if ( offset > 0 )then
          call ESMF_TimeIntervalSet( ETimeOff, s=offset, rc=rc )
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          ETimeAdd = ETime + ETimeOff
       else if ( offset < 0 )then
          call ESMF_TimeIntervalSet( ETimeOff, s=-offset, rc=rc )
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          ETimeAdd = ETime - ETimeOff
       end if
    end if

    call ESMF_TimeGet( ETimeAdd, yy=year, mm=month, dd=day, s=sec, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! shr_cal has restrictions and then "stops", so override that

    if ( present(ymd) ) then
       call shr_cal_ymd2date(year,month,day,ymd)
    endif
    if ( present(tod) ) then
       tod = sec
    endif

  end subroutine seq_timemgr_ETimeGet

  !===============================================================================
  subroutine seq_timemgr_EClockInit( TimeStep, StartTime, RefTime, CurrTime, EClock )

    ! !DESCRIPTION: Setup the ESMF clock
    use med_constants_mod, only : CL
    use ESMF, only: ESMF_Time, ESMF_TimeInterval, ESMF_Clock
    use ESMF, only: ESMF_ClockGet, ESMF_ClockAdvance, ESMF_ClockCreate
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_ChkErr
    use seq_comm_mct, only : loglevel, logunit
    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_TimeInterval), intent(in)  :: TimeStep    ! Time-step of clock
    type(ESMF_Time)        , intent(in)  :: StartTime   ! Start time
    type(ESMF_Time)        , intent(in)  :: RefTime     ! Reference time
    type(ESMF_Time)        , intent(in)  :: CurrTime    ! Current time
    type(ESMF_Clock)       , intent(out) :: EClock      ! Output ESMF clock

    !----- local -----
    integer                     :: rc          ! ESMF return code
    integer                     :: ymd, tod    ! time info
    character(len=CL)  :: description ! Description of this clock
    type(ESMF_Time)             :: clocktime   ! Current time
    character(len=*), parameter :: subname = '(seq_timemgr_EClockInit) '
    !-------------------------------------------------------------------------------

    description = 'ESMF Clock'

    ! ------ Create ESMF Clock with input characteristics -------------------
    ! --- NOTE: StopTime is required in interface but not used, so use  -----
    ! ---       something arbitrary.  Stop handled via alarm            -----

    call seq_timemgr_ETimeInit(clocktime,  99990101, 0, "artificial stop date")

    EClock = ESMF_ClockCreate(name=trim(description), &
         TimeStep=TimeStep, startTime=StartTime, refTime=RefTime, stopTime=clocktime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ------ Advance clock to the current time (in case of a restart) -------
    call ESMF_ClockGet(EClock, currTime=clocktime, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do while( clocktime < CurrTime)
       call ESMF_ClockAdvance( EClock, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockGet( EClock, currTime=clocktime, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    if (clocktime /= CurrTime) then
       if (loglevel > 0) then
          write(logunit,*) trim(subname),' : WARNING clocktime and currtime inconsistent'
          call seq_timemgr_ETimeGet( clocktime, ymd=ymd, tod=tod )
          write(logunit,*) trim(subname),' : clocktime = ',ymd,tod
          call seq_timemgr_ETimeGet( currtime, ymd=ymd, tod=tod )
          write(logunit,*) trim(subname),' : currtime  = ',ymd,tod
       endif
    endif

  end subroutine seq_timemgr_EClockInit

  !===============================================================================
  logical function seq_timemgr_EClockDateInSync( EClock, ymd, tod, prev)

    ! !DESCRIPTION: Check that the given input date/time is in sync with clock time
    use ESMF, only : ESMF_Clock, ESMF_ClockGet, ESMF_Time
    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock), intent(in) :: Eclock   ! Input clock to compare
    integer,          intent(in) :: ymd     ! Date (YYYYMMDD)
    integer,          intent(in) :: tod     ! Time of day (sec)
    logical, optional,intent(in) :: prev    ! If should get previous time

    !----- local -----
    type(ESMF_Time) :: ETime
    integer :: ymd1      ! Date (YYYYMMDD)
    integer :: tod1      ! Time of day
    logical :: previous  ! If need to get previous time for comparison
    integer :: rc        ! error code
    character(len=*), parameter :: subname = "(seq_timemgr_EClockDateInSync) "
    !-------------------------------------------------------------------------------

    previous = .false.
    if ( present(prev) )then
       previous = prev
    end if

    if (previous )then
       call ESMF_ClockGet( EClock, prevTime=ETime, rc=rc)
    else
       call ESMF_ClockGet( EClock, currTime=ETime, rc=rc)
    end if
    call seq_timemgr_ETimeGet( ETime, ymd=ymd1, tod=tod1 )

    ! --- If current dates agree return true -- else false

    if ( (ymd == ymd1) .and. (tod == tod1) )then
       seq_timemgr_EClockDateInSync = .true.
    else
       seq_timemgr_EClockDateInSync = .false.
    end if

  end function seq_timemgr_EClockDateInSync

  !===============================================================================
  subroutine seq_timemgr_clockPrint( SyncClock )

    ! !DESCRIPTION: Print clock information out.
    use med_constants_mod, only : in
    use ESMF, only : ESMF_Alarm, ESMF_ClockGetAlarmList
    use ESMF, only : ESMF_ALARMLIST_ALL
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_ChkErr
    use seq_comm_mct, only : loglevel, logunit
    ! !INPUT/OUTPUT PARAMETERS:
    type(seq_timemgr_type), intent(in) :: SyncClock   ! Input clock to print

    !----- local -----
    integer(IN) :: n
    character(len=*), parameter ::  F06 = "(2A,L3)"
    character(len=*), parameter ::  F07 = "(3A)"
    character(len=*), parameter :: subname = "(seq_timemgr_clockPrint) "
    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------

    if (loglevel <= 0) return

    write(logunit,F07) subname,'calendar      = ', trim(seq_timemgr_calendar)
    write(logunit,F06) subname,'end_restart   = ', seq_timemgr_end_restart
    write(logunit,F07) ''

    do n = 1,max_clocks
       call seq_timemgr_EClockPrint(SyncClock%ECP(n)%EClock, n)
    enddo

  end subroutine seq_timemgr_clockPrint

  !===============================================================================
  subroutine seq_timemgr_EClockPrint( EClock, n )
    use ESMF, only : ESMF_ClockGetAlarmList, ESMF_Clock, ESMF_Alarm
    use ESMF, only : ESMF_ALARMLIST_ALL
    use seq_comm_mct, only : loglevel, logunit
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_ChkErr
    ! !DESCRIPTION: Print clock information out.

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock), intent(in) :: EClock   ! Input clock to print
    integer, intent(in) :: n
    !----- local -----
    integer(IN) :: m
    integer(IN) :: curr_ymd   ! Current date YYYYMMDD
    integer(IN) :: curr_tod   ! Current time of day (s)
    integer(IN) :: StepNo     ! Number of steps taken
    integer(IN) :: start_ymd  ! Starting date YYYYMMDD
    integer(IN) :: start_tod  ! Starting time-of-day (s)
    integer(IN) :: stop_ymd   ! Stop date YYYYMMDD
    integer(IN) :: stop_tod   ! Stop time-of-day (s)
    integer(IN) :: ref_ymd    ! Reference date YYYYMMDD
    integer(IN) :: ref_tod    ! Reference time-of-day (s)
    integer(IN) :: DTime      ! Time-step (seconds)
    integer(IN) :: prev_ymd   ! Prev restart alarm date (YYYYMMDD)
    integer(IN) :: prev_tod   ! Prev restart alarm time-of-day (sec)
    integer(IN) :: next_ymd   ! Next restart alarm date (YYYYMMDD)
    integer(IN) :: next_tod   ! Next restart alarm time-of-day (sec)
    integer(IN) :: IntSec     ! Alarm interval for seconds
    integer(IN) :: IntMon     ! Alarm interval for months
    integer(IN) :: IntYrs     ! Alarm interval for years
    integer(IN) :: AlarmCount ! Number of valid alarms
    character(len=64)    :: alarmname  ! Alarm name
    integer(IN) :: rc         ! error code
    type(ESMF_Alarm), pointer :: EAlarm_list(:)   ! EAlarm list associated with EClock
    character(len=*), parameter :: xalarm = 'unset'
    character(len=*), parameter ::  F06 = "(2A,L3)"
    character(len=*), parameter ::  F07 = "(3A)"
    character(len=*), parameter ::  F08 = "(2A,I8.8,3x,I5.5)"
    character(len=*), parameter ::  F09 = "(2A,2I8,I12)"
    character(len=*), parameter ::  F10 = "(2A,I2,2x,A)"
    character(len=*), parameter :: subname = "(seq_timemgr_EClockPrint) "
    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------

    if (loglevel <= 0) return

       call seq_timemgr_EClockGetData( EClock, curr_ymd=curr_ymd, &
            curr_tod=curr_tod, start_ymd=start_ymd,    &
            start_tod=start_tod, StepNo=StepNo,            &
            ref_ymd=ref_ymd, ref_tod=ref_tod,              &
            stop_ymd=stop_ymd, stop_tod=stop_tod,          &
            dtime = dtime, alarmcount=AlarmCount)
       allocate(EAlarm_list(AlarmCount))
       call ESMF_ClockGetAlarmList(EClock, alarmListFlag=ESMF_ALARMLIST_ALL, &
            alarmList=EAlarm_list, alarmCount=AlarmCount, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       write(logunit,F09) subname,"Clock = "//seq_timemgr_clocks(n),n
       write(logunit,F08) subname,"  Start Time  = ", start_ymd, start_tod
       write(logunit,F08) subname,"  Curr Time   = ", curr_ymd, curr_tod
       write(logunit,F08) subname,"  Ref Time    = ", ref_ymd, ref_tod
       write(logunit,F08) subname,"  Stop Time   = ", stop_ymd, stop_tod
       write(logunit,F09) subname,"  Step number = ", StepNo
       write(logunit,F09) subname,"  Dtime       = ", DTime

       do m = 1,alarmCount
          call seq_timemgr_alarmGet( EAlarm_list(m), &
               next_ymd=next_ymd, next_tod=next_tod, prev_ymd=prev_ymd, prev_tod=prev_tod, &
               IntSec=IntSec, IntMon=IntMon, IntYrs=IntYrs, name=alarmname )
          write(logunit,F10) subname,"  Alarm = ",m,trim(alarmname)
          write(logunit,F08) subname,"    Prev Time   = ", prev_ymd,prev_tod
          write(logunit,F08) subname,"    Next Time   = ", next_ymd,next_tod
          write(logunit,F09) subname,"    Intervl yms = ", IntYrs,IntMon,IntSec
       enddo

       write(logunit,*) ''
       deallocate(EAlarm_list)

  end subroutine seq_timemgr_EClockPrint

  !===============================================================================

  subroutine seq_timemgr_ESMFDebug( EClock, ETime, ETimeInterval, istring )

    ! !DESCRIPTION: Print ESMF stuff for debugging
    use med_constants_mod, only : I8
    use ESMF, only : ESMF_Time, ESMF_TimeInterval, ESMF_TimeGet, ESMF_TimeIntervalGet
    use ESMF, only : ESMF_Clock, ESMF_ClockGet, ESMF_TimeIntervalGet
    use shr_nuopc_methods_mod, only : shr_nuopc_methods_ChkErr
    use seq_comm_mct, only : logunit
    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)        , optional, intent(in)    :: EClock        ! ESMF Clock
    type(ESMF_Time)         , optional, intent(inout) :: ETime         ! ESMF Time
    type(ESMF_TimeInterval) , optional, intent(inout) :: ETimeInterval ! ESMF Time Interval
    character(len=*)        , optional, intent(in)    :: istring

    !----- local -----
    character(len=128)      :: timestring
    integer                 :: yy,mm,dd,s            ! ymds
    type(ESMF_Time)         :: LTime
    type(ESMF_TimeInterval) :: LTimeInterval
    integer(I8)    :: LStep
    integer                 :: rc                    ! return code
    character(len=*), parameter :: subname = '(seq_timemgr_ESMFDebug) '
    !-------------------------------------------------------------------------------
    ! Notes:
    !-------------------------------------------------------------------------------

    if (present(ETime)) then
       write(logunit,*) subname,' ETime ',trim(istring)
       call ESMF_TimeGet(ETime, yy=yy,mm=mm,dd=dd,s=s,timestring=timestring,rc=rc)
       write(logunit,*) subname,rc,'ymds=',yy,mm,dd,s,trim(timestring)
    endif

    if (present(ETimeInterval)) then
       write(logunit,*) subname,' ETimeInterval ',trim(istring)
       call ESMF_TimeIntervalGet(ETimeInterval, yy=yy,mm=mm,d=dd,s=s,timestring=timestring,rc=rc)
       write(logunit,*) subname,rc,'ymds=',yy,mm,dd,s,trim(timestring)
    endif

    if (present(EClock)) then
       write(logunit,*) subname,' EClock ',trim(istring)

       call ESMF_ClockGet( EClock, StartTime=LTime, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(LTime, yy=yy,mm=mm,dd=dd,s=s,timestring=timestring,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       write(logunit,*) subname,rc,'start ymds=',yy,mm,dd,s,trim(timestring)

       call ESMF_ClockGet( EClock, CurrTime=LTime, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(LTime, yy=yy,mm=mm,dd=dd,s=s,timestring=timestring,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       write(logunit,*) subname,rc,'curr ymds=',yy,mm,dd,s,trim(timestring)

       call ESMF_ClockGet( EClock, StopTime=LTime, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(LTime, yy=yy,mm=mm,dd=dd,s=s,timestring=timestring,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       write(logunit,*) subname,rc,'stop ymds=',yy,mm,dd,s,trim(timestring)

       call ESMF_ClockGet( EClock, PrevTime=LTime, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(LTime, yy=yy,mm=mm,dd=dd,s=s,timestring=timestring,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       write(logunit,*) subname,rc,'prev ymds=',yy,mm,dd,s,trim(timestring)

       call ESMF_ClockGet( EClock, RefTime=LTime, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(LTime, yy=yy,mm=mm,dd=dd,s=s,timestring=timestring,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       write(logunit,*) subname,rc,'ref ymds=',yy,mm,dd,s,trim(timestring)

       call ESMF_ClockGet( EClock, TimeStep=LTimeInterval, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeIntervalGet(LTimeInterval, yy=yy,mm=mm,d=dd,s=s,timestring=timestring,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       write(logunit,*) subname,rc,'tint ymds=',yy,mm,dd,s,trim(timestring)

       call ESMF_ClockGet( EClock, AdvanceCount=LStep, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       write(logunit,*) subname,rc,'advcnt =',LStep
    endif

  end subroutine seq_timemgr_ESMFDebug

  !===============================================================================

end module seq_timemgr_mod
