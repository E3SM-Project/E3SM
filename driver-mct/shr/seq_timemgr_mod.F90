!===============================================================================
!
! !MODULE: seq_timemgr_mod --- Time-manager module
!
! !DESCRIPTION:
!
!     A module to create derived types to manage time and clock information 
!     for use with CCSM drivers and models.
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2005-Nov-11 - E. Kluzek - creation as eshr_timemgr_mod
!     2007-Sep-12 - T. Craig - extended
!     2007-Oct-05 - T. Craig - refactored to support concurrent models
!     2007-Nov-15 - T. Craig - refactored for ccsm4 and renamed seq_timemgr_mod
!
! !INTERFACE: ------------------------------------------------------------------

module seq_timemgr_mod

! !USES:
   use ESMF
   use shr_cal_mod
   use SHR_KIND_mod,     only: SHR_KIND_IN, SHR_KIND_R8, SHR_KIND_CS, &
                               SHR_KIND_CL, SHR_KIND_I8
   use seq_comm_mct,     only: logunit, loglevel, seq_comm_iamin, CPLID, &
                               seq_comm_gloroot, seq_comm_iamroot
   use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush

   implicit none

   private    ! default private

! ! PUBLIC TYPES:

   public :: seq_timemgr_type         ! Wrapped clock object

! ! PUBLIC MEMBER FUNCTIONS:

   ! --- Clock object methods --------------------------------------------------
   public :: seq_timemgr_clockInit         ! Setup the sync clock
   public :: seq_timemgr_clockAdvance      ! Advance the sync clock
   public :: seq_timemgr_clockPrint        ! Print sync clock information 

   public :: seq_timemgr_EClockGetData     ! Get data from an ESMF clock

   public :: seq_timemgr_EClockDateInSync  ! compare EClock to ymd/tod
   public :: seq_timemgr_alarmSetOn        ! Turn an alarm on
   public :: seq_timemgr_alarmSetOff       ! Turn an alarm off
   public :: seq_timemgr_alarmIsOn         ! Is an alarm ringing
   public :: seq_timemgr_ETimeInit         ! Create ESMF_Time object
   public :: seq_timemgr_ETimeGet          ! Query ESMF_Time object

   ! --- For usability, built on interfaces above ---
   public :: seq_timemgr_restartAlarmIsOn    ! Is a restart alarm ringing
   public :: seq_timemgr_stopAlarmIsOn       ! Is a stop alarm ringing
   public :: seq_timemgr_historyAlarmIsOn    ! Is a history alarm ringing

   ! --- Data that belongs to the driver (but is here to avoid loops)
   character(SHR_KIND_CS),public :: seq_timemgr_pause_component_list  ! Pause - resume components

 ! ! PUBLIC PARAMETERS:

   integer(SHR_KIND_IN),public :: seq_timemgr_histavg_type 
   integer(SHR_KIND_IN),public,parameter :: seq_timemgr_type_other  = -1
   integer(SHR_KIND_IN),public,parameter :: seq_timemgr_type_never  = 1
   integer(SHR_KIND_IN),public,parameter :: seq_timemgr_type_nhour  = 2
   integer(SHR_KIND_IN),public,parameter :: seq_timemgr_type_nday   = 3
   integer(SHR_KIND_IN),public,parameter :: seq_timemgr_type_nmonth = 4
   integer(SHR_KIND_IN),public,parameter :: seq_timemgr_type_nyear  = 5

   character(SHR_KIND_CL),public,parameter :: seq_timemgr_noleap    = shr_cal_noleap
   character(SHR_KIND_CL),public,parameter :: seq_timemgr_gregorian = shr_cal_gregorian

!  These are public but declared in the private area for clarity

!  clocknames: 
!   character(len=*),public,parameter :: &
!      seq_timemgr_clock_drv 
!      seq_timemgr_clock_atm 
!      seq_timemgr_clock_lnd 
!      seq_timemgr_clock_rof
!      seq_timemgr_clock_ocn 
!      seq_timemgr_clock_ice 
!      seq_timemgr_clock_glc
!      seq_timemgr_clock_wav
!      seq_timemgr_clock_esp

!  alarmnames:
!   character(len=*),public,parameter :: &
!      seq_timemgr_alarm_restart
!      seq_timemgr_alarm_run    
!      seq_timemgr_alarm_stop   
!      seq_timemgr_alarm_datestop
!      seq_timemgr_alarm_history
!      seq_timemgr_alarm_atmrun 
!      seq_timemgr_alarm_lndrun 
!      seq_timemgr_alarm_rofrun 
!      seq_timemgr_alarm_ocnrun 
!      seq_timemgr_alarm_icerun 
!      seq_timemgr_alarm_glcrun 
!      seq_timemgr_alarm_glcrun_avg 
!      seq_timemgr_alarm_wavrun 
!      seq_timemgr_alarm_esprun
!      seq_timemgr_alarm_ocnnext
!      seq_timemgr_alarm_tprof
!      seq_timemgr_alarm_histavg
!      seq_timemgr_alarm_pause
!      seq_timemgr_alarm_barrier

   private:: seq_timemgr_alarmGet
   private:: seq_timemgr_alarmInit
   private:: seq_timemgr_EClockInit
   private:: seq_timemgr_ESMFDebug
   private:: seq_timemgr_ESMFCodeCheck

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
      seq_timemgr_optEnd            = "end"     

   integer(SHR_KIND_IN),private,parameter :: &
      seq_timemgr_nclock_drv  = 1, &
      seq_timemgr_nclock_atm  = 2, &
      seq_timemgr_nclock_lnd  = 3, &
      seq_timemgr_nclock_ocn  = 4, &
      seq_timemgr_nclock_ice  = 5, &
      seq_timemgr_nclock_glc  = 6, &
      seq_timemgr_nclock_wav  = 7, &
      seq_timemgr_nclock_rof  = 8, &
      seq_timemgr_nclock_esp  = 9

   integer(SHR_KIND_IN),private,parameter :: max_clocks = 9
   character(len=*),public,parameter :: &
      seq_timemgr_clock_drv  = 'seq_timemgr_clock_drv' , & 
      seq_timemgr_clock_atm  = 'seq_timemgr_clock_atm' , &
      seq_timemgr_clock_lnd  = 'seq_timemgr_clock_lnd' , &
      seq_timemgr_clock_ocn  = 'seq_timemgr_clock_ocn' , &
      seq_timemgr_clock_ice  = 'seq_timemgr_clock_ice' , &
      seq_timemgr_clock_glc  = 'seq_timemgr_clock_glc' , &
      seq_timemgr_clock_wav  = 'seq_timemgr_clock_wav' , &
      seq_timemgr_clock_rof  = 'seq_timemgr_clock_rof' , &
      seq_timemgr_clock_esp  = 'seq_timemgr_clock_esp' 
   character(len=8),private,parameter :: seq_timemgr_clocks(max_clocks) = &
      (/'drv     ','atm     ','lnd     ','ocn     ', &
        'ice     ','glc     ','wav     ','rof     ','esp     '/)

   ! Alarms on both component clocks and driver clock
   integer(SHR_KIND_IN),private,parameter :: &
      seq_timemgr_nalarm_restart    = 1 , & ! driver and component clock alarm
      seq_timemgr_nalarm_run        = 2 , & ! driver and component clock alarm
      seq_timemgr_nalarm_stop       = 3 , & ! driver and component clock alarm
      seq_timemgr_nalarm_datestop   = 4 , & ! driver and component clock alarm
      seq_timemgr_nalarm_history    = 5 , & ! driver and component clock alarm
      seq_timemgr_nalarm_atmrun     = 6 , & ! driver only clock alarm
      seq_timemgr_nalarm_lndrun     = 7 , & ! driver only clock alarm
      seq_timemgr_nalarm_ocnrun     = 8 , & ! driver only clock alarm
      seq_timemgr_nalarm_icerun     = 9 , & ! driver only clock alarm
      seq_timemgr_nalarm_glcrun     =10 , & ! driver only clock alarm
      seq_timemgr_nalarm_glcrun_avg =11 , & ! driver only clock alarm
      seq_timemgr_nalarm_ocnnext    =12 , & ! driver only clock alarm
      seq_timemgr_nalarm_tprof      =13 , & ! driver and component clock alarm
      seq_timemgr_nalarm_histavg    =14 , & ! driver and component clock alarm
      seq_timemgr_nalarm_rofrun     =15 , & ! driver only clock alarm
      seq_timemgr_nalarm_wavrun     =16 , & ! driver only clock alarm
      seq_timemgr_nalarm_esprun     =17 , & ! driver only clock alarm
      seq_timemgr_nalarm_pause      =18 , &
      seq_timemgr_nalarm_barrier    =19 , & ! driver and component clock alarm
      max_alarms = seq_timemgr_nalarm_barrier

   character(len=*),public,parameter :: &
      seq_timemgr_alarm_restart    = 'seq_timemgr_alarm_restart ', &
      seq_timemgr_alarm_run        = 'seq_timemgr_alarm_run     ', &
      seq_timemgr_alarm_stop       = 'seq_timemgr_alarm_stop    ', &
      seq_timemgr_alarm_datestop   = 'seq_timemgr_alarm_datestop', &
      seq_timemgr_alarm_history    = 'seq_timemgr_alarm_history ', &
      seq_timemgr_alarm_atmrun     = 'seq_timemgr_alarm_atmrun  ', &
      seq_timemgr_alarm_lndrun     = 'seq_timemgr_alarm_lndrun  ', &
      seq_timemgr_alarm_ocnrun     = 'seq_timemgr_alarm_ocnrun  ', &
      seq_timemgr_alarm_icerun     = 'seq_timemgr_alarm_icerun  ', &
      seq_timemgr_alarm_glcrun     = 'seq_timemgr_alarm_glcrun  ', &
      seq_timemgr_alarm_glcrun_avg = 'seq_timemgr_alarm_glcrun_avg' , &
      seq_timemgr_alarm_ocnnext    = 'seq_timemgr_alarm_ocnnext ', &
      seq_timemgr_alarm_tprof      = 'seq_timemgr_alarm_tprof   ', &
      seq_timemgr_alarm_histavg    = 'seq_timemgr_alarm_histavg ', &
      seq_timemgr_alarm_rofrun     = 'seq_timemgr_alarm_rofrun  ', &
      seq_timemgr_alarm_wavrun     = 'seq_timemgr_alarm_wavrun  ', &
      seq_timemgr_alarm_esprun     = 'seq_timemgr_alarm_esprun  ', &
      seq_timemgr_alarm_pause      = 'seq_timemgr_alarm_pause   ', &
      seq_timemgr_alarm_barrier    = 'seq_timemgr_alarm_barrier '

   type EClock_pointer     ! needed for array of pointers
      type(ESMF_Clock),pointer :: EClock => null()
   end type EClock_pointer

   type seq_timemgr_type
      private
      type(EClock_pointer) :: ECP(max_clocks)    ! ESMF clocks, array of pointers
      type(ESMF_Alarm) :: EAlarm(max_clocks,max_alarms) ! array of clock alarms
   end type seq_timemgr_type
 
   ! --- Private local data -------------------------------------------------------

   type(ESMF_Calendar), target, save :: seq_timemgr_cal        ! calendar
   character(SHR_KIND_CL)      ,save :: seq_timemgr_calendar   ! calendar string
   logical :: seq_timemgr_end_restart  ! write restarts at end of run?
   logical, save :: seq_timemgr_setCalendar = .false.          ! if calendar has been set
   integer, parameter :: SecPerDay = 86400         ! Seconds per day

!===============================================================================

contains

!===============================================================================
! !IROUTINE: seq_timemgr_clockInit -- Initializes clocks
!   
! !DESCRIPTION:
!   
!     Initializes clock
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_timemgr_clockInit(SyncClock, nmlfile, restart, restart_file, pioid, mpicom, &
      EClock_drv, EClock_atm, EClock_lnd, EClock_ocn, EClock_ice, Eclock_glc, &
      Eclock_rof, EClock_wav, Eclock_esp)
                                           
! !USES:
  use pio, only : file_desc_T
   use shr_string_mod, only : shr_string_toupper
   use shr_file_mod,   only : shr_file_getunit, shr_file_freeunit
   use shr_mpi_mod,    only : shr_mpi_bcast
   use seq_io_read_mod

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(seq_timemgr_type),  intent(INOUT) :: SyncClock    ! sync clock
   character(len=*),        intent(IN)    :: nmlfile      ! namelist file
   integer,                 intent(IN)    :: mpicom       ! MPI communicator
   logical,                 intent(IN)    :: restart      ! restart logical
   character(len=*),        intent(IN)    :: restart_file
   type(ESMF_clock),target, intent(IN)    :: EClock_drv   ! drv clock
   type(ESMF_clock),target, intent(IN)    :: EClock_atm   ! atm clock
   type(ESMF_clock),target, intent(IN)    :: EClock_lnd   ! lnd clock
   type(ESMF_clock),target, intent(IN)    :: EClock_ocn   ! ocn clock
   type(ESMF_clock),target, intent(IN)    :: EClock_ice   ! ice clock
   type(ESMF_clock),target, intent(IN)    :: EClock_glc   ! glc clock
   type(ESMF_clock),target, intent(IN)    :: EClock_rof   ! rof clock
   type(ESMF_clock),target, intent(IN)    :: EClock_wav   ! wav clock
   type(ESMF_clock),target, intent(IN)    :: EClock_esp   ! esp clock
   type(file_desc_t)                      :: pioid

    !----- local -----
    character(len=*), parameter :: subname = '(seq_timemgr_clockInit) '
    type(ESMF_Time)             :: StartTime          ! Start time
    type(ESMF_Time)             :: RefTime            ! Reference time
    type(ESMF_Time)             :: CurrTime           ! Current time
    type(ESMF_Time)             :: OffsetTime         ! local computed time
    type(ESMF_Time)             :: StopTime1          ! Stop time
    type(ESMF_Time)             :: StopTime2          ! Stop time
    type(ESMF_TimeInterval)     :: TimeStep           ! Clock time-step
    type(ESMF_TimeInterval)     :: AlarmInterval      ! Alarm interval
    type(ESMF_CalKind_Flag)     :: esmf_caltype       ! local esmf calendar
    integer                     :: rc                 ! Return code
    integer                     :: n                  ! index
    integer                     :: dtime(max_clocks)  ! time-step to use
    integer                     :: offset(max_clocks) ! run offset
    integer                     :: unitn              ! i/o unit number
    integer                     :: iam                ! pe rank

    character(SHR_KIND_CS)  :: calendar              ! Calendar type
    character(SHR_KIND_CS)  :: stop_option           ! Stop option units
    integer(SHR_KIND_IN)    :: stop_n                ! Number until stop
    integer(SHR_KIND_IN)    :: stop_ymd              ! Stop date (YYYYMMDD)
    integer(SHR_KIND_IN)    :: stop_tod              ! Stop time-of-day
    character(SHR_KIND_CS)  :: restart_option        ! Restart option units
    integer(SHR_KIND_IN)    :: restart_n             ! Number until restart interval
    integer(SHR_KIND_IN)    :: restart_ymd           ! Restart date (YYYYMMDD)
    character(SHR_KIND_CS)  :: pause_option          ! Pause option units
    integer(SHR_KIND_IN)    :: pause_n               ! Number between pause intervals
    character(SHR_KIND_CS)  :: pause_component_list  ! Pause - resume components
    character(SHR_KIND_CS)  :: history_option        ! History option units
    integer(SHR_KIND_IN)    :: history_n             ! Number until history interval
    integer(SHR_KIND_IN)    :: history_ymd           ! History date (YYYYMMDD)
    character(SHR_KIND_CS)  :: histavg_option        ! Histavg option units
    integer(SHR_KIND_IN)    :: histavg_n             ! Number until histavg interval
    integer(SHR_KIND_IN)    :: histavg_ymd           ! Histavg date (YYYYMMDD)
    character(SHR_KIND_CS)  :: barrier_option        ! Barrier option units
    integer(SHR_KIND_IN)    :: barrier_n             ! Number until barrier interval
    integer(SHR_KIND_IN)    :: barrier_ymd           ! Barrier date (YYYYMMDD)
    character(SHR_KIND_CS)  :: tprof_option          ! tprof option units
    integer(SHR_KIND_IN)    :: tprof_n               ! Number until tprof interval
    integer(SHR_KIND_IN)    :: tprof_ymd             ! tprof date (YYYYMMDD)
    integer(SHR_KIND_IN)    :: start_ymd             ! Start date (YYYYMMDD)
    integer(SHR_KIND_IN)    :: start_tod             ! Start time of day (seconds)
    integer(SHR_KIND_IN)    :: curr_ymd              ! Current ymd (YYYYMMDD)
    integer(SHR_KIND_IN)    :: curr_tod              ! Current tod (seconds)
    integer(SHR_KIND_IN)    :: ref_ymd               ! Reference date (YYYYMMDD)
    integer(SHR_KIND_IN)    :: ref_tod               ! Reference time of day (seconds)
    integer(SHR_KIND_IN)    :: atm_cpl_dt            ! Atmosphere coupling interval
    integer(SHR_KIND_IN)    :: lnd_cpl_dt            ! Land coupling interval
    integer(SHR_KIND_IN)    :: ice_cpl_dt            ! Sea-Ice coupling interval
    integer(SHR_KIND_IN)    :: ocn_cpl_dt            ! Ocean coupling interval
    integer(SHR_KIND_IN)    :: glc_cpl_dt            ! Glc coupling interval
    integer(SHR_KIND_IN)    :: rof_cpl_dt            ! Runoff coupling interval
    integer(SHR_KIND_IN)    :: wav_cpl_dt            ! Wav coupling interval
    integer(SHR_KIND_IN)    :: esp_cpl_dt            ! Esp coupling interval
    integer(SHR_KIND_IN)    :: atm_cpl_offset        ! Atmosphere coupling interval
    integer(SHR_KIND_IN)    :: lnd_cpl_offset        ! Land coupling interval
    integer(SHR_KIND_IN)    :: ice_cpl_offset        ! Sea-Ice coupling interval
    integer(SHR_KIND_IN)    :: ocn_cpl_offset        ! Ocean coupling interval
    integer(SHR_KIND_IN)    :: glc_cpl_offset        ! Glc coupling interval
    integer(SHR_KIND_IN)    :: wav_cpl_offset        ! Wav coupling interval
    integer(SHR_KIND_IN)    :: rof_cpl_offset        ! Runoff coupling interval
    integer(SHR_KIND_IN)    :: esp_cpl_offset        ! Esp coupling interval
    logical                 :: end_restart           ! Write restart at end of run
    integer(SHR_KIND_IN)    :: nlUnit                ! Namelist unit number
    integer(SHR_KIND_IN)    :: ierr                  ! Return code

    character(len=*), parameter ::  F0A = "(2A,A)"
    character(len=*), parameter ::  F0I = "(2A,I10)"
    character(len=*), parameter ::  F0L = "(2A,L3)"

    namelist /seq_timemgr_inparm/  calendar, curr_ymd, curr_tod,  &
         stop_option, stop_n, stop_ymd, stop_tod,                &
         restart_option, restart_n, restart_ymd,                 &
         pause_option,   pause_n,   pause_component_list,        &
         history_option, history_n, history_ymd,                 &
         histavg_option, histavg_n, histavg_ymd,                 &
         barrier_option, barrier_n, barrier_ymd,                 &
         tprof_option, tprof_n, tprof_ymd,                       &
         start_ymd, start_tod, ref_ymd, ref_tod,                 &
         atm_cpl_dt, ocn_cpl_dt, ice_cpl_dt, lnd_cpl_dt,         &
         atm_cpl_offset, lnd_cpl_offset, ocn_cpl_offset,         &
         ice_cpl_offset, glc_cpl_dt, glc_cpl_offset,             &
         wav_cpl_dt, wav_cpl_offset, esp_cpl_dt, esp_cpl_offset, &
         rof_cpl_dt, rof_cpl_offset, end_restart
!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    SyncClock%ECP(seq_timemgr_nclock_drv)%EClock     => EClock_drv
    SyncClock%ECP(seq_timemgr_nclock_atm)%EClock     => EClock_atm
    SyncClock%ECP(seq_timemgr_nclock_lnd)%EClock     => EClock_lnd
    SyncClock%ECP(seq_timemgr_nclock_ocn)%EClock     => EClock_ocn
    SyncClock%ECP(seq_timemgr_nclock_ice)%EClock     => EClock_ice
    SyncClock%ECP(seq_timemgr_nclock_glc)%EClock     => EClock_glc
    SyncClock%ECP(seq_timemgr_nclock_rof)%EClock     => EClock_rof
    SyncClock%ECP(seq_timemgr_nclock_wav)%EClock     => EClock_wav
    SyncClock%ECP(seq_timemgr_nclock_esp)%EClock     => EClock_esp

    call mpi_comm_rank(mpicom,iam,ierr)

    !---------------------------------------------------------------------------
    ! Set syncclock on root pe
    !---------------------------------------------------------------------------

    if (iam == 0) then

       !---------------------------------------------------------------------------
       ! Set namelist defaults
       !---------------------------------------------------------------------------
       calendar         = seq_timemgr_noleap
       stop_option      = ' '
       stop_n           = -1
       stop_ymd         = -1
       stop_tod         = 0
       restart_option   = seq_timemgr_optYearly
       restart_n        = -1
       restart_ymd      = -1
       pause_option     = seq_timemgr_optNever
       pause_n          = -1
       pause_component_list = ' '
       history_option   = seq_timemgr_optNever
       history_n        = -1
       history_ymd      = -1
       histavg_option   = seq_timemgr_optNever
       histavg_n        = -1
       histavg_ymd      = -1
       barrier_option   = seq_timemgr_optNever
       barrier_n        = -1
       barrier_ymd      = -1
       tprof_option     = seq_timemgr_optNever
       tprof_n          = -1
       tprof_ymd        = -1
       start_ymd        = 0
       start_tod        = 0
       ref_ymd          = 0
       ref_tod          = 0
       curr_ymd         = 0
       curr_tod         = 0
       atm_cpl_dt       = 0
       lnd_cpl_dt       = 0
       ice_cpl_dt       = 0
       ocn_cpl_dt       = 0
       glc_cpl_dt       = 0
       rof_cpl_dt       = 0
       wav_cpl_dt       = 0
       esp_cpl_dt       = 0
       atm_cpl_offset   = 0
       lnd_cpl_offset   = 0
       ice_cpl_offset   = 0
       ocn_cpl_offset   = 0
       glc_cpl_offset   = 0
       rof_cpl_offset   = 0
       wav_cpl_offset   = 0
       esp_cpl_offset   = 0
       end_restart      = .true.

       !---------------------------------------------------------------------------
       ! Read in namelist
       !---------------------------------------------------------------------------
       unitn = shr_file_getUnit()
       write(logunit,F0A) trim(subname),' Read in seq_timemgr_inparm namelist from: '//trim(nmlfile)
       open( unitn, file=trim(nmlfile), status='old' )
       ierr = 1
       do while( ierr /= 0 )
          read(unitn,nml=seq_timemgr_inparm,iostat=ierr)
          if (ierr < 0) then
             call shr_sys_abort( subname//':: namelist read returns an'// &
                                 ' end of file or end of record condition' )
          end if
       end do
       close(unitn)
       call shr_file_freeUnit( unitn )
    endif

    !---------------------------------------------------------------------------
    ! Read Restart (seq_io is called on all CPLID pes)
    ! NOTE: slightly messy, seq_io is only valid on CPLID
    !---------------------------------------------------------------------------
    if (restart) then
       if (seq_comm_iamin(CPLID)) then
          call seq_io_read(restart_file,pioid,start_ymd,'seq_timemgr_start_ymd')
          call seq_io_read(restart_file,pioid,start_tod,'seq_timemgr_start_tod')
          call seq_io_read(restart_file,pioid,ref_ymd  ,'seq_timemgr_ref_ymd')
          call seq_io_read(restart_file,pioid,ref_tod  ,'seq_timemgr_ref_tod')
          call seq_io_read(restart_file,pioid,curr_ymd ,'seq_timemgr_curr_ymd')
          call seq_io_read(restart_file,pioid,curr_tod ,'seq_timemgr_curr_tod')
       endif
       !--- Send from CPLID ROOT to GLOBALID ROOT, use bcast as surrogate
       call shr_mpi_bcast(start_ymd,mpicom,pebcast=seq_comm_gloroot(CPLID))
       call shr_mpi_bcast(start_tod,mpicom,pebcast=seq_comm_gloroot(CPLID))
       call shr_mpi_bcast(  ref_ymd,mpicom,pebcast=seq_comm_gloroot(CPLID))
       call shr_mpi_bcast(  ref_tod,mpicom,pebcast=seq_comm_gloroot(CPLID))
       call shr_mpi_bcast( curr_ymd,mpicom,pebcast=seq_comm_gloroot(CPLID))
       call shr_mpi_bcast( curr_tod,mpicom,pebcast=seq_comm_gloroot(CPLID))
    endif

    if (iam == 0) then
       !---------------------------------------------------------------------------
       ! Modify namelist as needed
       !---------------------------------------------------------------------------

       if (lnd_cpl_dt == 0) lnd_cpl_dt = atm_cpl_dt ! Copy atm coupling time into lnd
       if (rof_cpl_dt == 0) rof_cpl_dt = atm_cpl_dt ! Copy atm coupling time into rof
       if (ice_cpl_dt == 0) ice_cpl_dt = atm_cpl_dt ! Copy atm coupling time into ice
       if (ocn_cpl_dt == 0) ocn_cpl_dt = atm_cpl_dt ! Copy atm coupling time into ocn
       if (glc_cpl_dt == 0) glc_cpl_dt = atm_cpl_dt ! Copy atm coupling time into glc
       if (wav_cpl_dt == 0) wav_cpl_dt = atm_cpl_dt ! Copy atm coupling time into wav
       if (esp_cpl_dt == 0) esp_cpl_dt = atm_cpl_dt ! Copy atm coupling time into esp

       if (glc_cpl_avg_dt == 0) then
          ! set default average coupling interval 
          glc_cpl_avg_dt = glc_cpl_dt
       end if

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
       write(logunit,F0I) trim(subname),' glc_cpl_avg_dt = ',glc_cpl_avg_dt
       write(logunit,F0I) trim(subname),' rof_cpl_dt     = ',rof_cpl_dt
       write(logunit,F0I) trim(subname),' wav_cpl_dt     = ',wav_cpl_dt
       write(logunit,F0I) trim(subname),' esp_cpl_dt     = ',esp_cpl_dt
       write(logunit,F0I) trim(subname),' atm_cpl_offset = ',atm_cpl_offset
       write(logunit,F0I) trim(subname),' lnd_cpl_offset = ',lnd_cpl_offset
       write(logunit,F0I) trim(subname),' ice_cpl_offset = ',ice_cpl_offset
       write(logunit,F0I) trim(subname),' ocn_cpl_offset = ',ocn_cpl_offset
       write(logunit,F0I) trim(subname),' glc_cpl_offset = ',glc_cpl_offset
       write(logunit,F0I) trim(subname),' rof_cpl_offset = ',rof_cpl_offset
       write(logunit,F0I) trim(subname),' wav_cpl_offset = ',wav_cpl_offset
       write(logunit,F0I) trim(subname),' esp_cpl_offset = ',esp_cpl_offset
       write(logunit,F0A) ' '

       !---------------------------------------------------------------------------
       ! Check a few things
       !---------------------------------------------------------------------------

       ! --- Coupling intervals ------------------------------------------------
       if ( atm_cpl_dt <= 0 .or. &
            lnd_cpl_dt /= atm_cpl_dt .or. &
            ice_cpl_dt /= atm_cpl_dt .or. &
            ocn_cpl_dt <= 0 .or. glc_cpl_dt <= 0 .or. rof_cpl_dt <=0 .or. &
            wav_cpl_dt <=0 .or. esp_cpl_dt <=0) then
          write(logunit,*) trim(subname),' ERROR: aliogrwe _cpl_dt = ', &
             atm_cpl_dt, lnd_cpl_dt, ice_cpl_dt, ocn_cpl_dt, glc_cpl_dt, &
             rof_cpl_dt, wav_cpl_dt, esp_cpl_dt
          call shr_sys_abort( subname//': ERROR coupling intervals invalid' )
       end if

       ! --- Coupling offsets --------------------------------------------------
       if ( abs(atm_cpl_offset) > atm_cpl_dt .or. &
            abs(lnd_cpl_offset) > lnd_cpl_dt .or. &
            abs(ice_cpl_offset) > ice_cpl_dt .or. &
            abs(glc_cpl_offset) > glc_cpl_dt .or. &
            abs(rof_cpl_offset) > rof_cpl_dt .or. &
            abs(wav_cpl_offset) > wav_cpl_dt .or. &
            abs(esp_cpl_offset) > esp_cpl_dt .or. &
            abs(ocn_cpl_offset) > ocn_cpl_dt) then
          write(logunit,*) trim(subname),' ERROR: aliogrwe _cpl_offset = ', &
             atm_cpl_offset, lnd_cpl_offset, ice_cpl_offset, ocn_cpl_offset, &
             glc_cpl_offset, rof_cpl_offset, wav_cpl_offset, esp_cpl_offset
          call shr_sys_abort( subname//': ERROR coupling offsets invalid' )
       end if

       ! --- Start time date ---------------------------------------------------
       if ( (start_ymd < 101) .or. (start_ymd > 99991231)) then
          write(logunit,*) subname,' ERROR: illegal start_ymd',start_ymd
          call shr_sys_abort( subname//': ERROR invalid start_ymd')
       end if

    endif

    !---------------------------------------------------------------------------
    ! Broadcast namelist data
    !---------------------------------------------------------------------------
    call shr_mpi_bcast( calendar,             mpicom )
    call shr_mpi_bcast( stop_n,               mpicom )
    call shr_mpi_bcast( stop_option,          mpicom )
    call shr_mpi_bcast( stop_ymd,             mpicom )
    call shr_mpi_bcast( stop_tod,             mpicom )
    call shr_mpi_bcast( restart_n,            mpicom )
    call shr_mpi_bcast( restart_option,       mpicom )
    call shr_mpi_bcast( restart_ymd,          mpicom )
    call shr_mpi_bcast( pause_n,              mpicom )
    call shr_mpi_bcast( pause_option,         mpicom )
    call shr_mpi_bcast( pause_component_list, mpicom )
    call shr_mpi_bcast( history_n,            mpicom )
    call shr_mpi_bcast( history_option,       mpicom )
    call shr_mpi_bcast( history_ymd,          mpicom )
    call shr_mpi_bcast( histavg_n,            mpicom )
    call shr_mpi_bcast( histavg_option,       mpicom )
    call shr_mpi_bcast( histavg_ymd,          mpicom )
    call shr_mpi_bcast( tprof_n,              mpicom )
    call shr_mpi_bcast( barrier_n,            mpicom )
    call shr_mpi_bcast( barrier_option,       mpicom )
    call shr_mpi_bcast( barrier_ymd,          mpicom )
    call shr_mpi_bcast( tprof_option,         mpicom )
    call shr_mpi_bcast( tprof_ymd,            mpicom )
    call shr_mpi_bcast( start_ymd,            mpicom )
    call shr_mpi_bcast( start_tod,            mpicom )
    call shr_mpi_bcast( ref_ymd,              mpicom )
    call shr_mpi_bcast( ref_tod,              mpicom )
    call shr_mpi_bcast( curr_ymd,             mpicom )
    call shr_mpi_bcast( curr_tod,             mpicom )
    call shr_mpi_bcast( atm_cpl_dt,           mpicom )
    call shr_mpi_bcast( lnd_cpl_dt,           mpicom )
    call shr_mpi_bcast( ice_cpl_dt,           mpicom )
    call shr_mpi_bcast( ocn_cpl_dt,           mpicom )
    call shr_mpi_bcast( glc_cpl_dt,           mpicom )
    call shr_mpi_bcast( glc_cpl_avg_dt,       mpicom )
    call shr_mpi_bcast( rof_cpl_dt,           mpicom )
    call shr_mpi_bcast( wav_cpl_dt,           mpicom )
    call shr_mpi_bcast( esp_cpl_dt,           mpicom )
    call shr_mpi_bcast( atm_cpl_offset,       mpicom )
    call shr_mpi_bcast( lnd_cpl_offset,       mpicom )
    call shr_mpi_bcast( ice_cpl_offset,       mpicom )
    call shr_mpi_bcast( ocn_cpl_offset,       mpicom )
    call shr_mpi_bcast( glc_cpl_offset,       mpicom )
    call shr_mpi_bcast( rof_cpl_offset,       mpicom )
    call shr_mpi_bcast( wav_cpl_offset,       mpicom )
    call shr_mpi_bcast( esp_cpl_offset,       mpicom )
    call shr_mpi_bcast( end_restart,          mpicom )

    ! --- derive a couple things ---
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
    seq_timemgr_end_restart          = end_restart
    seq_timemgr_pause_component_list = pause_component_list

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

    seq_timemgr_cal = ESMF_CalendarCreate( name='CCSM_'//seq_timemgr_calendar, &
                      calkindflag=esmf_caltype, rc=rc )
    call seq_timemgr_ESMFCodeCheck( rc, subname//': error return from ESMF_CalendarCreate' )

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
    dtime(seq_timemgr_nclock_glc_avg ) = glc_cpl_avg_dt
    dtime(seq_timemgr_nclock_rof     ) = rof_cpl_dt
    dtime(seq_timemgr_nclock_wav     ) = wav_cpl_dt
    dtime(seq_timemgr_nclock_esp     ) = esp_cpl_dt

    ! --- this finds the min of dtime excluding the driver value ---
    dtime(seq_timemgr_nclock_drv) = maxval(dtime)
    dtime(seq_timemgr_nclock_drv) = minval(dtime)

    do n = 1,max_clocks
       if ( mod(dtime(n),dtime(seq_timemgr_nclock_drv)) /= 0) then
          write(logunit,*) trim(subname),' ERROR: dtime inconsistent = ',dtime
          call shr_sys_abort( subname//' :coupling intervals not compatible' )
       endif
    enddo

    ! --- Initialize component and driver clocks and alarms common to components amd drivver clocks ---

    do n = 1,max_clocks
       call ESMF_TimeIntervalSet( TimeStep, s=dtime(n), rc=rc )
       call seq_timemgr_ESMFCodeCheck( rc, subname//': error ESMF_TimeIntervalSet' )

       call seq_timemgr_EClockInit( TimeStep, StartTime, RefTime, CurrTime, SyncClock%ECP(n)%EClock)

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
          EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_run),  &
          option  = seq_timemgr_optNSeconds,         &
          opt_n   = dtime(n), RefTime = CurrTime,    &
          alarmname = trim(seq_timemgr_alarm_run))

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
          EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_stop),  &
          option  = stop_option,         &
          opt_n   = stop_n,              &
          opt_ymd = stop_ymd,            &
          opt_tod = stop_tod,            &
          RefTime = CurrTime,            &
          alarmname = trim(seq_timemgr_alarm_stop))

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
          EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_datestop),  &
          option  = seq_timemgr_optDate, &
          opt_ymd = stop_ymd,            &
          opt_tod = stop_tod,            &
          RefTime = StartTime,            &
          alarmname = trim(seq_timemgr_alarm_datestop))

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
          EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_restart),  &
          option  = restart_option,      &
          opt_n   = restart_n,           &
          opt_ymd = restart_ymd,         &
          RefTime = CurrTime,            &
          alarmname = trim(seq_timemgr_alarm_restart))

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
          EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_history),  &
          option  = history_option,      &
          opt_n   = history_n,           &
          opt_ymd = history_ymd,         &
          RefTime = StartTime,           &
          alarmname = trim(seq_timemgr_alarm_history))

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
          EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_histavg),  &
          option  = histavg_option,      &
          opt_n   = histavg_n,           &
          opt_ymd = histavg_ymd,         &
          RefTime = StartTime,           &
          alarmname = trim(seq_timemgr_alarm_histavg))

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
          EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_barrier),  &
          option  = barrier_option,      &
          opt_n   = barrier_n,           &
          opt_ymd = barrier_ymd,         &
          RefTime = CurrTime,            &
          alarmname = trim(seq_timemgr_alarm_barrier))

       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock, &
          EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_tprof),  &
          option  = tprof_option,        &
          opt_n   = tprof_n,             &
          opt_ymd = tprof_ymd,           &
          RefTime = StartTime,           &
          alarmname = trim(seq_timemgr_alarm_tprof))

       call ESMF_AlarmGet(SyncClock%EAlarm(n,seq_timemgr_nalarm_stop), RingTime=StopTime1, rc=rc )
       call ESMF_AlarmGet(SyncClock%EAlarm(n,seq_timemgr_nalarm_datestop), RingTime=StopTime2, rc=rc )
       if (StopTime2 < StopTime1) then
          call ESMF_ClockSet(SyncClock%ECP(n)%EClock, StopTime=StopTime2)
       else
          call ESMF_ClockSet(SyncClock%ECP(n)%EClock, StopTime=StopTime1)
       endif

       ! Set the pause option if pause/resume is active
       call seq_timemgr_alarmInit(SyncClock%ECP(n)%EClock,                    &
            EAlarm  = SyncClock%EAlarm(n,seq_timemgr_nalarm_pause),           &
            option  = pause_option,                                           &
            opt_n   = pause_n,                                                &
            RefTime = StartTime,                                              &
            alarmname = trim(seq_timemgr_alarm_pause))

    enddo

    ! --------------------------------------------------------------------
    ! Set the timing run alarms, these alarms are synced to the driver
    ! clock and determine when the component clocks are advanced.
    ! We need an offset here of the driver timestep because of the
    ! implementation.  We are advancing the clock first and we want
    ! components to run as soon as possible.  Without the driver offset
    ! the alarms would go off at the last possible timestep, not first.
    ! In addition, we allow the user to set other offsets if desired
    ! via namelist.  tcraig, 10/2007
    ! --------------------------------------------------------------------

    offset(seq_timemgr_nclock_drv)     = 0
    offset(seq_timemgr_nclock_atm)     = atm_cpl_offset
    offset(seq_timemgr_nclock_lnd)     = lnd_cpl_offset
    offset(seq_timemgr_nclock_ocn)     = ocn_cpl_offset
    offset(seq_timemgr_nclock_ice)     = ice_cpl_offset
    offset(seq_timemgr_nclock_glc)     = glc_cpl_offset
    offset(seq_timemgr_nclock_glc_avg) = glc_cpl_offset
    offset(seq_timemgr_nclock_rof)     = rof_cpl_offset
    offset(seq_timemgr_nclock_wav)     = wav_cpl_offset
    offset(seq_timemgr_nclock_esp)     = esp_cpl_offset

    do n = 1,max_clocks
       if (abs(offset(n)) > dtime(n)) then
          write(logunit,*) subname,' ERROR: offset too large',n,dtime(n),offset(n)
          call shr_sys_abort()
       endif

       !--- this is the required driver timestep offset ---
       offset(n) = offset(n) + dtime(seq_timemgr_nclock_drv)

       if (mod(offset(n),dtime(seq_timemgr_nclock_drv)) /= 0) then
          write(logunit,*) subname,' ERROR: offset not multiple',n,dtime(seq_timemgr_nclock_drv),offset(n)
          call shr_sys_abort()
       endif
    enddo

    ! Set component run alarms on driver clock

    call ESMF_TimeIntervalSet( TimeStep, s=offset(seq_timemgr_nclock_atm), rc=rc )
    OffsetTime = CurrTime + TimeStep
    call seq_timemgr_alarmInit(SyncClock%ECP(seq_timemgr_nclock_drv)%EClock, &
       EAlarm  = SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_atmrun),  &
       option  = seq_timemgr_optNSeconds,       &
       opt_n   = dtime(seq_timemgr_nclock_atm), &
       RefTime = OffsetTime,                    &
       alarmname = trim(seq_timemgr_alarm_atmrun))

    call ESMF_TimeIntervalSet( TimeStep, s=offset(seq_timemgr_nclock_lnd), rc=rc )
    OffsetTime = CurrTime + TimeStep
    call seq_timemgr_alarmInit(SyncClock%ECP(seq_timemgr_nclock_drv)%EClock, &
       EAlarm  = SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_lndrun),  &
       option  = seq_timemgr_optNSeconds,       &
       opt_n   = dtime(seq_timemgr_nclock_lnd), &
       RefTime = OffsetTime,                    &
       alarmname = trim(seq_timemgr_alarm_lndrun))

    call ESMF_TimeIntervalSet( TimeStep, s=offset(seq_timemgr_nclock_rof), rc=rc )
    OffsetTime = CurrTime + TimeStep
    call seq_timemgr_alarmInit(SyncClock%ECP(seq_timemgr_nclock_drv)%EClock, &
       EAlarm  = SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_rofrun),  &
       option  = seq_timemgr_optNSeconds,       &
       opt_n   = dtime(seq_timemgr_nclock_rof), &
       RefTime = OffsetTime,                    &
       alarmname = trim(seq_timemgr_alarm_rofrun))

    call ESMF_TimeIntervalSet( TimeStep, s=offset(seq_timemgr_nclock_ice), rc=rc )
    OffsetTime = CurrTime + TimeStep
    call seq_timemgr_alarmInit(SyncClock%ECP(seq_timemgr_nclock_drv)%EClock, &
       EAlarm  = SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_icerun),  &
       option  = seq_timemgr_optNSeconds,       &
       opt_n   = dtime(seq_timemgr_nclock_ice), &
       RefTime = OffsetTime,                    &
       alarmname = trim(seq_timemgr_alarm_icerun))

    call ESMF_TimeIntervalSet( TimeStep, s=offset(seq_timemgr_nclock_wav), rc=rc )
    OffsetTime = CurrTime + TimeStep
    call seq_timemgr_alarmInit(SyncClock%ECP(seq_timemgr_nclock_drv)%EClock, &
       EAlarm  = SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_wavrun),  &
       option  = seq_timemgr_optNSeconds,       &
       opt_n   = dtime(seq_timemgr_nclock_wav), &
       RefTime = OffsetTime,                    &
       alarmname = trim(seq_timemgr_alarm_wavrun))

    call ESMF_TimeIntervalSet( TimeStep, s=offset(seq_timemgr_nclock_esp), rc=rc )
    OffsetTime = CurrTime + TimeStep
    call seq_timemgr_alarmInit(SyncClock%ECP(seq_timemgr_nclock_drv)%EClock, &
       EAlarm  = SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_esprun),  &
       option  = seq_timemgr_optNSeconds,       &
       opt_n   = dtime(seq_timemgr_nclock_esp), &
       RefTime = OffsetTime,                    &
       alarmname = trim(seq_timemgr_alarm_esprun))

    ! --- this is the glcrun alarm (there ^) offset by a -dtime of the driver
    call ESMF_TimeIntervalSet( TimeStep, s=offset(seq_timemgr_nclock_glc), rc=rc )
    OffsetTime = CurrTime + TimeStep
    call ESMF_TimeIntervalSet( TimeStep, s=-offset(seq_timemgr_nclock_drv), rc=rc )
    OffsetTime = OffsetTime + TimeStep
    call seq_timemgr_alarmInit(SyncClock%ECP(seq_timemgr_nclock_drv)%EClock, &
       EAlarm  = SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_glcrun),  &
       option  = seq_timemgr_optNSeconds,       &
       opt_n   = dtime(seq_timemgr_nclock_glc), &
       RefTime = OffsetTime,                    &
       alarmname = trim(seq_timemgr_alarm_glcrun))

    ! --- this is the glcrun_avg alarm (there ^) offset by a -dtime of the driver
    call ESMF_TimeIntervalSet( TimeStep, s=offset(seq_timemgr_nclock_glc_avg), rc=rc )
    OffsetTime = CurrTime + TimeStep
    call ESMF_TimeIntervalSet( TimeStep, s=-offset(seq_timemgr_nclock_drv), rc=rc )
    OffsetTime = OffsetTime + TimeStep
    call seq_timemgr_alarmInit(SyncClock%ECP(seq_timemgr_nclock_drv)%EClock, &
       EAlarm  = SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_glcrun),  &
       option  = seq_timemgr_optNSeconds,       &
       opt_n   = glc_cpl_avg_dt,                &
       RefTime = OffsetTime,                    &
       alarmname = trim(seq_timemgr_alarm_glcrun_avg))

    call ESMF_TimeIntervalSet( TimeStep, s=offset(seq_timemgr_nclock_ocn), rc=rc )
    OffsetTime = CurrTime + TimeStep
    call seq_timemgr_alarmInit(SyncClock%ECP(seq_timemgr_nclock_drv)%EClock, &
       EAlarm  = SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_ocnrun),  &
       option  = seq_timemgr_optNSeconds,       &
       opt_n   = dtime(seq_timemgr_nclock_ocn), &
       RefTime = OffsetTime,                    &
       alarmname = trim(seq_timemgr_alarm_ocnrun))

    ! --- this is the ocnrun alarm (there ^) offset by a -dtime of the driver
    call ESMF_TimeIntervalSet( TimeStep, s=offset(seq_timemgr_nclock_ocn), rc=rc )
    OffsetTime = CurrTime + TimeStep
    call ESMF_TimeIntervalSet( TimeStep, s=-offset(seq_timemgr_nclock_drv), rc=rc )
    OffsetTime = OffsetTime + TimeStep
    call seq_timemgr_alarmInit(SyncClock%ECP(seq_timemgr_nclock_drv)%EClock, &
       EAlarm  = SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_ocnnext),  &
       option  = seq_timemgr_optNSeconds,       &
       opt_n   = dtime(seq_timemgr_nclock_ocn), &
       RefTime = OffsetTime,                    &
       alarmname = trim(seq_timemgr_alarm_ocnnext))

end subroutine seq_timemgr_clockInit

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_EClockGetData -- Get information from the clock
!   
! !DESCRIPTION:
!   
!     Get various values from the clock.
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_timemgr_EClockGetData( EClock, curr_yr, curr_mon, curr_day,    &
                       curr_ymd, curr_tod, prev_ymd, prev_tod, start_ymd,     &
                       start_tod, StepNo, ref_ymd, ref_tod,         &
                       stop_ymd, stop_tod, dtime, ECurrTime, alarmcount,      &
                       curr_cday, next_cday, curr_time, prev_time, calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock),     intent(IN)            :: EClock     ! Input clock object
    integer(SHR_KIND_IN), intent(OUT), optional :: curr_yr    ! Current year
    integer(SHR_KIND_IN), intent(OUT), optional :: curr_mon   ! Current month
    integer(SHR_KIND_IN), intent(OUT), optional :: curr_day   ! Current day in month
    integer(SHR_KIND_IN), intent(OUT), optional :: curr_ymd   ! Current date YYYYMMDD
    integer(SHR_KIND_IN), intent(OUT), optional :: curr_tod   ! Current time of day (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: prev_ymd   ! Previous date YYYYMMDD
    integer(SHR_KIND_IN), intent(OUT), optional :: prev_tod   ! Previous time of day (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: start_ymd  ! Starting date YYYYMMDD
    integer(SHR_KIND_IN), intent(OUT), optional :: start_tod  ! Starting time-of-day (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: StepNo     ! Number of steps taken
    integer(SHR_KIND_IN), intent(OUT), optional :: ref_ymd    ! Reference date YYYYMMDD
    integer(SHR_KIND_IN), intent(OUT), optional :: ref_tod    ! Reference time-of-day (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: stop_ymd   ! Stop date YYYYMMDD
    integer(SHR_KIND_IN), intent(OUT), optional :: stop_tod   ! Stop time-of-day (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: dtime      ! Time-step (seconds)
    integer(SHR_KIND_IN), intent(OUT), optional :: alarmcount ! Number of Valid Alarms
    type(ESMF_Time),      intent(OUT), optional :: ECurrTime  ! Current ESMF time
    real(SHR_KIND_R8)   , intent(OUT), optional :: curr_cday  ! current calendar day
    real(SHR_KIND_R8)   , intent(OUT), optional :: next_cday  ! current calendar day
    real(SHR_KIND_R8)   , intent(OUT), optional :: curr_time  ! time interval between current time 
                                                              ! and reference date
    real(SHR_KIND_R8)   , intent(OUT), optional :: prev_time  ! time interval between previous time
                                                              ! and reference date
    character(len=*)    , intent(OUT), optional :: calendar   ! calendar type

    !----- local -----
    character(len=*), parameter :: subname = '(seq_timemgr_EClockGetData) '
    type(ESMF_Time) :: CurrentTime            ! Current time
    type(ESMF_Time) :: PreviousTime           ! Previous time
    type(ESMF_Time) :: StartTime              ! Start time
    type(ESMF_Time) :: StopTime               ! Stop time
    type(ESMF_Time) :: RefTime                ! Ref time
    type(ESMF_TimeInterval) :: timeStep       ! Clock, time-step
    type(ESMF_TimeInterval) :: timediff       ! Used to calculate curr_time
    integer(SHR_KIND_IN) :: rc                ! Return code
    integer(SHR_KIND_I8) :: advSteps          ! Number of time-steps that have advanced
    integer(SHR_KIND_IN) :: yy, mm, dd, sec   ! Return time values
    integer(SHR_KIND_IN) :: ymd               ! Date (YYYYMMDD)
    integer(SHR_KIND_IN) :: tod               ! time of day (sec)
    integer(SHR_KIND_IN) :: ldtime            ! local dtime
    integer(SHR_KIND_IN) :: intyrs            ! alarm variable
    integer(SHR_KIND_IN) :: intmon            ! alarm variable
    integer(SHR_KIND_IN) :: intsec            ! alarm variable
    integer(SHR_KIND_IN) :: days              ! number of whole days in time interval
    integer(SHR_KIND_IN) :: seconds           ! number of seconds in time interval
    integer(SHR_KIND_IN) :: acount            ! number of valid alarms
    real(SHR_KIND_R8) :: doy, tmpdoy          ! day of year
    real(SHR_KIND_R8),parameter :: c1 = 1.0_SHR_KIND_R8

    type(ESMF_Time) :: tmpTime                ! tmp time, needed for next_cday
    type(ESMF_TimeInterval) :: tmpDTime       ! tmp time interval, needed for next_cday

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------
    if (present(calendar)) calendar = trim(seq_timemgr_calendar)

    call ESMF_ClockGet( EClock, currTime=CurrentTime, &
         advanceCount=advSteps, prevTime=previousTime, TimeStep=timeStep, &
         startTime=StartTime, stopTime=stopTime, refTime=RefTime, &
         AlarmCount=acount, rc=rc )
    call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_ClockGet" )

    call ESMF_TimeGet( CurrentTime, yy=yy, mm=mm, dd=dd, s=sec, dayofyear_r8=doy, rc=rc )
    call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_TimeGet" )
    call seq_timemgr_ETimeGet( CurrentTime, ymd=ymd, tod=tod )
    call ESMF_TimeIntervalGet( timeStep, s=ldtime, rc=rc )
    call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_TimeIntervalGet" )

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
       call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from TimeSet tmpTime")
       call ESMF_TimeIntervalSet( tmpDTime, d=1, rc=rc )
       call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from TimeIntSet tmpDTime")
       tmpTime = tmpTime + tmpDTime
       call ESMF_TimeGet(tmpTime, dayOfYear_r8=tmpdoy, rc=rc)
       call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from TimeGet tmpdoy")
       next_cday = tmpdoy
    endif

    ! ---Current Time (the time interval between the current date and the reference date) ---
    if ( present(curr_time)) then 
       timediff = CurrentTime - RefTime
       call ESMF_TimeIntervalGet(timediff, d=days, s=seconds, rc=rc)
       call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from  TimeIntervalGet timediff")
       curr_time = days + seconds/real(SecPerDay,SHR_KIND_R8)
    end if

    ! ---Previous Time (the time interval between the previous date and the reference date) ---
    if ( present(prev_time)) then 
       timediff = PreviousTime - RefTime
       call ESMF_TimeIntervalGet(timediff, d=days, s=seconds, rc=rc)
       call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from  TimeIntervalGet timediff")
       prev_time = days + seconds/real(SecPerDay,SHR_KIND_R8)
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
!===============================================================================
! !IROUTINE: seq_timemgr_clockAdvance  -- Advance the syncclock
!   
! !DESCRIPTION:
!   
! Advance this clock
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_timemgr_clockAdvance( SyncClock, force_stop, force_stop_ymd, force_stop_tod )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(seq_timemgr_type), intent(INOUT) :: SyncClock    ! Advancing clock
   logical, optional, intent(in) :: force_stop           ! force stop
   integer, optional, intent(in) :: force_stop_ymd       ! force stop ymd
   integer, optional, intent(in) :: force_stop_tod       ! force stop tod

    !----- local -----
    character(len=*), parameter :: subname = '(seq_timemgr_clockAdvance) '
    integer :: n    
    type(ESMF_Time) :: NextAlarm              ! Next restart alarm time
    integer :: rc    ! Return code

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   ! --- set datestop alarm to force_stop alarm ---

   do n = 1,max_clocks
      call seq_timemgr_alarmSetOff(SyncClock%ECP(n)%EClock)
      if (present(force_stop) .and. present(force_stop_ymd) .and. present(force_stop_tod)) then
      if (force_stop) then
         if (n == 1 .and. seq_comm_iamroot(CPLID)) then
            write(logunit,*) subname,'force stop at ',force_stop_ymd, force_stop_tod
         endif
         if (force_stop_ymd < 0 .or. force_stop_tod < 0) then
            call shr_sys_abort(subname//': force_stop_ymd, force_stop_tod invalid')
         endif
         seq_timemgr_end_restart = .true.
         call seq_timemgr_ETimeInit(NextAlarm, force_stop_ymd, force_stop_tod, "optDate")
         CALL ESMF_AlarmSet( SyncClock%EAlarm(n,seq_timemgr_nalarm_datestop),  &
                             name = trim(seq_timemgr_alarm_datestop), &
                             RingTime=NextAlarm,        &
                             rc=rc )
      endif
      endif
   enddo

   ! --- advance driver clock and all driver alarms ---

   call ESMF_ClockAdvance( SyncClock%ECP(seq_timemgr_nclock_drv)%EClock, rc=rc )
   call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from drv ESMF_ClockAdvance")

   ! --- advance other clocks if driver component run alarm is ringing ---

   if (ESMF_AlarmIsRinging(SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_atmrun))) then
      call ESMF_ClockAdvance(SyncClock%ECP(seq_timemgr_nclock_atm)%EClock, rc=rc )
      call seq_timemgr_ESMFCodeCheck(rc, msg=subname//"Error from atm ESMF_ClockAdvance")
   endif

   if (ESMF_AlarmIsRinging(SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_lndrun))) then
      call ESMF_ClockAdvance(SyncClock%ECP(seq_timemgr_nclock_lnd)%EClock, rc=rc )
      call seq_timemgr_ESMFCodeCheck(rc, msg=subname//"Error from lnd ESMF_ClockAdvance")
   endif

   if (ESMF_AlarmIsRinging(SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_rofrun))) then
      call ESMF_ClockAdvance(SyncClock%ECP(seq_timemgr_nclock_rof)%EClock, rc=rc )
      call seq_timemgr_ESMFCodeCheck(rc, msg=subname//"Error from rof ESMF_ClockAdvance")
   endif

   if (ESMF_AlarmIsRinging(SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_ocnrun))) then
      call ESMF_ClockAdvance(SyncClock%ECP(seq_timemgr_nclock_ocn)%EClock, rc=rc )
      call seq_timemgr_ESMFCodeCheck(rc, msg=subname//"Error from ocn ESMF_ClockAdvance")
   endif

   if (ESMF_AlarmIsRinging(SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_icerun))) then
      call ESMF_ClockAdvance(SyncClock%ECP(seq_timemgr_nclock_ice)%EClock, rc=rc )
      call seq_timemgr_ESMFCodeCheck(rc, msg=subname//"Error from ice ESMF_ClockAdvance")
   endif

   if (ESMF_AlarmIsRinging(SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_glcrun))) then
      call ESMF_ClockAdvance(SyncClock%ECP(seq_timemgr_nclock_glc)%EClock, rc=rc )
      call seq_timemgr_ESMFCodeCheck(rc, msg=subname//"Error from glc ESMF_ClockAdvance")
   endif

   if (ESMF_AlarmIsRinging(SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_wavrun))) then
      call ESMF_ClockAdvance(SyncClock%ECP(seq_timemgr_nclock_wav)%EClock, rc=rc )
      call seq_timemgr_ESMFCodeCheck(rc, msg=subname//"Error from wav ESMF_ClockAdvance")
   endif

   if (ESMF_AlarmIsRinging(SyncClock%EAlarm(seq_timemgr_nclock_drv,seq_timemgr_nalarm_esprun))) then
      call ESMF_ClockAdvance(SyncClock%ECP(seq_timemgr_nclock_esp)%EClock, rc=rc )
      call seq_timemgr_ESMFCodeCheck(rc, msg=subname//"Error from esp ESMF_ClockAdvance")
   endif

   if (seq_timemgr_end_restart) then
      do n = 1,max_clocks
         if (seq_timemgr_alarmIsOn(SyncClock%ECP(n)%EClock,seq_timemgr_alarm_stop) .or. &
             seq_timemgr_alarmIsOn(SyncClock%ECP(n)%EClock,seq_timemgr_alarm_datestop)) then
            call seq_timemgr_alarmSetOn(SyncClock%ECP(n)%EClock,seq_timemgr_alarm_restart)
         endif
      enddo
   endif

end subroutine seq_timemgr_clockAdvance

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_alarmInit -- Set an alarm
!   
! !DESCRIPTION:
!   
!     Setup an alarm in a clock
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_timemgr_alarmInit( EClock, EAlarm, option, opt_n, opt_ymd, opt_tod, RefTime, alarmname)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)             , intent(INOUT) :: EClock    ! clock
    type(ESMF_Alarm)             , intent(INOUT) :: EAlarm    ! alarm
    character(len=*)             , intent(IN)    :: option    ! alarm option
    integer(SHR_KIND_IN),optional, intent(IN)    :: opt_n     ! alarm freq
    integer(SHR_KIND_IN),optional, intent(IN)    :: opt_ymd   ! alarm ymd
    integer(SHR_KIND_IN),optional, intent(IN)    :: opt_tod   ! alarm tod (sec)
    type(ESMF_Time)     ,optional, intent(IN)    :: RefTime   ! ref time
    character(len=*)    ,optional, intent(IN)    :: alarmname ! alarm name

    !----- local -----
    character(len=*), parameter :: subname = '(seq_timemgr_alarmInit): '
    integer :: rc                             ! Return code
    integer :: lymd                           ! local ymd
    integer :: ltod                           ! local tod
    integer :: cyy,cmm,cdd,csec               ! time info
    integer :: nyy,nmm,ndd,nsec               ! time info
    character(len=64) :: lalarmname           ! local alarm name
    logical :: update_nextalarm               ! update next alarm
    type(ESMF_Time) :: CurrTime               ! Current Time
    type(ESMF_Time) :: NextAlarm              ! Next restart alarm time
    type(ESMF_Time) :: AltAlarm               ! Alternate alarm time
    type(ESMF_TimeInterval) :: AlarmInterval  ! Alarm interval

!-------------------------------------------------------------------------------
! Notes: This is slightly screwed up because of the way the ESMF alarm
!        initializes.  The ringtime sent to AlarmCreate MUST be the next
!        alarm time.  If you send an arbitrary but proper ringtime from 
!        the past and the ring interval, the alarm will always go off on
!        the next clock advance and this will cause serious problems.
!        So, even if it makes sense to initialize an alarm with some
!        reference time and the alarm interval, that reference time has
!        to be advance forward to be >= the current time.  In the logic
!        below, we set an appropriate "NextAlarm" and then we make sure
!        to advance it properly based on the ring interval.
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
    call ESMF_TimeGet(CurrTime, yy=cyy, mm=cmm, dd=cdd, s=csec, rc=rc )

    ! --- initial guess of next alarm, this will be updated below ---
    if (present(RefTime)) then
       NextAlarm = RefTime
    else
       NextAlarm = CurrTime
    endif
    call ESMF_TimeGet(CurrTime, yy=nyy, mm=nmm, dd=ndd, s=nsec, rc=rc )

    update_nextalarm  = .true.

    selectcase (trim(option))

    case (seq_timemgr_optNONE)
       !--- tcx seems we need an alarm interval or the alarm create fails, 
       !--- problem in esmf_wrf_timemgr?
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=seq_timemgr_cal, rc=rc )
       update_nextalarm  = .false.

    case (seq_timemgr_optNever)
       !--- tcx seems we need an alarm interval or the alarm create fails, 
       !--- problem in esmf_wrf_timemgr?
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=seq_timemgr_cal, rc=rc )
       update_nextalarm  = .false.

    case (seq_timemgr_optDate)
       !--- tcx seems we need an alarm interval or the alarm create fails, 
       !--- problem in esmf_wrf_timemgr?
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (.not. present(opt_ymd)) call shr_sys_abort(subname//trim(option)//' requires opt_ymd')
       if (lymd < 0 .or. ltod < 0) then
          call shr_sys_abort(subname//trim(option)//'opt_ymd, opt_tod invalid')
       endif
       call seq_timemgr_ETimeInit(NextAlarm, lymd, ltod, "optDate")
       update_nextalarm  = .false.

    case (seq_timemgr_optIfdays0)
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (.not. present(opt_ymd)) call shr_sys_abort(subname//trim(option)//' requires opt_ymd')
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=opt_n, s=0, calendar=seq_timemgr_cal, rc=rc )

    case (seq_timemgr_optNSteps)
       call ESMF_ClockGet(EClock, TimeStep=AlarmInterval, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNStep)
       call ESMF_ClockGet(EClock, TimeStep=AlarmInterval, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNSeconds)
       call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNSecond)
       call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNMinutes)
       call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNMinute)
       call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNHours)
       call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNHour)
       call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n

    case (seq_timemgr_optNDays)
       call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n
!       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=cdd, s=0, calendar=seq_timemgr_cal, rc=rc )

    case (seq_timemgr_optNDay)
       call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n
!       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=cdd, s=0, calendar=seq_timemgr_cal, rc=rc )

    case (seq_timemgr_optNMonths)
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n
!       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=1, s=0, calendar=seq_timemgr_cal, rc=rc )

    case (seq_timemgr_optNMonth)
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n
!       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=1, s=0, calendar=seq_timemgr_cal, rc=rc )

    case (seq_timemgr_optMonthly)
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=1, s=0, calendar=seq_timemgr_cal, rc=rc )

    case (seq_timemgr_optNYears)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n
!       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=1, dd=1, s=0, calendar=seq_timemgr_cal, rc=rc )

    case (seq_timemgr_optNYear)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (.not.present(opt_n)) call shr_sys_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_sys_abort(subname//trim(option)//' invalid opt_n')
       AlarmInterval = AlarmInterval * opt_n
!       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=1, dd=1, s=0, calendar=seq_timemgr_cal, rc=rc )

    case (seq_timemgr_optYearly)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=1, dd=1, s=0, calendar=seq_timemgr_cal, rc=rc )

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

    EAlarm = ESMF_AlarmCreate( name=lalarmname, clock=EClock, ringTime=NextAlarm, &
       ringInterval=AlarmInterval, rc=rc)
    call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_AlarmCreate" )

end subroutine seq_timemgr_AlarmInit

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_alarmGet -- Get information from the alarm
!   
! !DESCRIPTION:
!   
!     Get various values from the alarm.
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_timemgr_alarmGet( EAlarm, next_ymd, next_tod, prev_ymd, prev_tod,    &
                                 IntSec, IntMon, IntYrs, name)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Alarm)    , intent(INOUT)            :: EAlarm   ! Input Alarm object
    integer(SHR_KIND_IN), intent(OUT), optional :: next_ymd ! alarm date yyyymmdd
    integer(SHR_KIND_IN), intent(OUT), optional :: next_tod ! alarm tod sec
    integer(SHR_KIND_IN), intent(OUT), optional :: prev_ymd ! alarm date yyyymmdd
    integer(SHR_KIND_IN), intent(OUT), optional :: prev_tod ! alarm tod sec
    integer(SHR_KIND_IN), intent(OUT), optional :: IntSec   ! alarm int sec
    integer(SHR_KIND_IN), intent(OUT), optional :: IntMon   ! alarm int mon
    integer(SHR_KIND_IN), intent(OUT), optional :: IntYrs   ! alarm int yrs
    character(len=*)    , intent(OUT), optional :: name     ! alarm name

    !----- local -----
    character(len=*), parameter :: subname = '(seq_timemgr_alarmGet) '
    integer :: yy, mm, dd, sec                ! Return time values
    integer :: ymd                            ! Date (YYYYMMDD)
    integer :: tod                            ! time of day (sec)
    integer :: rc                             ! error code
    type(ESMF_TimeInterval) :: alarmInterval  ! Alarm interval
    type(ESMF_Time) :: ringTime               ! Next alarm ring time

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    if (present(name)) then
       call ESMF_AlarmGet( EAlarm, name=name, rc=rc)
       call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_AlarmGet name" )
    endif

    call ESMF_AlarmGet( EAlarm, RingTime=RingTime, rc=rc )
    call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_AlarmGet RingTime" )
    call seq_timemgr_ETimeGet( RingTime, ymd=ymd, tod=tod)
    if ( present(next_ymd) ) next_ymd = ymd
    if ( present(next_tod) ) next_tod = tod

    call ESMF_AlarmGet( EAlarm, PrevRingTime=RingTime, rc=rc )
    call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_AlarmGet PrevRingTime")
    call seq_timemgr_ETimeGet( RingTime, ymd=ymd, tod=tod)
    if ( present(prev_ymd) ) prev_ymd = ymd
    if ( present(prev_tod) ) prev_tod = tod

    yy = 0
    mm = 0
    dd = 0
    sec = 0
    call ESMF_AlarmGet( EAlarm, RingInterval=AlarmInterval, rc=rc )
    call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_AlarmGet RingInterval")
    call ESMF_TimeIntervalGet( alarmInterval, yy=yy, mm=mm, d=dd, s=sec, rc=rc )
    call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_TimeIntervalGet" )
    sec = sec + dd*(SecPerDay)

    ! --- If want restart next interval information -------------------------
    if ( present(IntSec) ) IntSec = sec
    if ( present(IntMon) ) IntMon = mm
    if ( present(IntYrs) ) IntYrs = yy

end subroutine seq_timemgr_alarmGet
!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_alarmSetOn -- turn alarm on
!   
! !DESCRIPTION:
!   
!     turn alarm on
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_timemgr_AlarmSetOn( EClock, alarmname)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(INOUT) :: EClock      ! clock/alarm
    character(len=*), intent(IN), optional :: alarmname  ! alarmname

    !----- local -----
    character(len=*), parameter :: subname = '(seq_timemgr_alarmSetOn) '
    character(len=*), parameter :: xalarm = 'unset'
    integer :: n
    integer :: rc
    logical :: found
    logical :: set
    character(len=64) :: name
    type(ESMF_Alarm),pointer :: EAlarm
    type(ESMF_Alarm),pointer :: EAlarm_list(:)
    integer(SHR_KIND_IN) :: AlarmCount ! Number of valid alarms

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
#ifdef USE_ESMF_LIB
    allocate(EAlarm_list(AlarmCount))
    call ESMF_ClockGetAlarmList(EClock, alarmListFlag=ESMF_ALARMLIST_ALL, &
       alarmList=EAlarm_list, alarmCount=AlarmCount, rc=rc)
#else
    call ESMF_ClockGetAlarmList(EClock, EAlarm_list, rc=rc)
#endif
    call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_ClockGetAlarmList" )
    do n = 1,AlarmCount
       found = .false.
       if (present(alarmname)) then
          call ESMF_AlarmGet(EAlarm_list(n), name=name)
          if (trim(name) == trim(alarmname)) found = .true.
       else
          found = .true.
       endif
       if (found) then
          set = .true.
          call ESMF_AlarmRingerOn( EAlarm_list(n), rc=rc )
          call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_AlarmRingerOn" )
       endif
    enddo

    if (present(alarmname) .and. .not. set) then
       write(logunit,*) subname,' ERROR in alarmname ',trim(alarmname)
       call shr_sys_abort()
    endif
#ifdef USE_ESMF_LIB
    deallocate(EAlarm_list)
#endif

end subroutine seq_timemgr_AlarmSetOn

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_alarmSetOff -- turn alarm off
!   
! !DESCRIPTION:
!   
!     turn alarm off
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_timemgr_AlarmSetOff( EClock, alarmname)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(INOUT) :: EClock      ! clock/alarm
    character(len=*), intent(IN), optional :: alarmname  ! alarmname

    !----- local -----
    character(len=*), parameter :: subname = '(seq_timemgr_alarmSetOff) '
    character(len=*), parameter :: xalarm = 'unset'
    integer :: n
    integer :: rc
    logical :: found
    logical :: set
    character(len=64) :: name
    type(ESMF_Alarm),pointer :: EAlarm
    type(ESMF_Alarm),pointer :: EAlarm_list(:)
    integer(SHR_KIND_IN) :: AlarmCount ! Number of valid alarms

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
#ifdef USE_ESMF_LIB
    allocate(EAlarm_list(AlarmCount))
    call ESMF_ClockGetAlarmList(EClock, alarmListFlag=ESMF_ALARMLIST_ALL, &
       alarmList=EAlarm_list, alarmCount=AlarmCount, rc=rc)
#else
    call ESMF_ClockGetAlarmList(EClock, EAlarm_list, rc=rc)
#endif
    call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_ClockGetAlarmList" )
    do n = 1,AlarmCount
       found = .false.
       if (present(alarmname)) then
          call ESMF_AlarmGet(EAlarm_list(n), name=name)
          if (trim(name) == trim(alarmname)) found = .true.
       else
          found = .true.
       endif
       if (found) then
          set = .true.
          call ESMF_AlarmRingerOff( EAlarm_list(n), rc=rc )
          call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_AlarmRingerOff" )
       endif
    enddo

    if (present(alarmname) .and. .not. set) then
       write(logunit,*) subname,' ERROR in alarmname ',trim(alarmname)
       call shr_sys_abort()
    endif
#ifdef USE_ESMF_LIB
    deallocate(EAlarm_list)
#endif

end subroutine seq_timemgr_AlarmSetOff

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_alarmIsOn -- check if an alarm is ringing
!   
! !DESCRIPTION:
!   
!     check if an alarm is ringing
!      
! !INTERFACE: ------------------------------------------------------------------

logical function seq_timemgr_alarmIsOn( EClock, alarmname)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(IN) :: EClock     ! clock/alarm
    character(len=*), intent(IN) :: alarmname  ! which alarm

    !----- local -----
    character(len=*), parameter :: subname = '(seq_timemgr_alarmIsOn) '
    character(len=*), parameter :: xalarm = 'unset'
    integer :: n
    integer :: rc
    logical :: found
    logical :: set
    character(len=64) :: name
    type(ESMF_Time) :: ETime1, ETime2
    type(ESMF_Alarm),pointer :: EAlarm
    type(ESMF_Alarm),pointer :: EAlarm_list(:)
    integer(SHR_KIND_IN) :: AlarmCount ! Number of valid alarms

!-------------------------------------------------------------------------------
! Notes:  Because of the esmf_wrf_timemgr implementation with regards to
!   valid alarms in the alarm_list, we initialize name to xalarm before
!   querying the alarm name, and if the alarm is not valid, name will not
!   be updated and we can tell that the alarm is not valid and we should
!   just ignore it.
!         Use found to verify alarm was valid.  If not, abort
!-------------------------------------------------------------------------------

    seq_timemgr_alarmIsOn = .false.
    found = .false.

    call seq_timemgr_EClockGetData(EClock, AlarmCount=AlarmCount)
#ifdef USE_ESMF_LIB
    allocate(EAlarm_list(AlarmCount))
    call ESMF_ClockGetAlarmList(EClock, alarmListFlag=ESMF_ALARMLIST_ALL, &
       alarmList=EAlarm_list, alarmCount=AlarmCount, rc=rc)
#else
    call ESMF_ClockGetAlarmList(EClock, EAlarm_list, rc=rc)
#endif
    call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_ClockGetAlarmList" )
    do n = 1,AlarmCount
       name = trim(xalarm)
       call ESMF_AlarmGet(EAlarm_list(n), name=name)
       if (trim(name) == trim(alarmname)) then
          found = .true.
          seq_timemgr_alarmIsOn = ESMF_AlarmIsRinging(alarm=EAlarm_list(n),rc=rc)
          call seq_timemgr_ESMFCodeCheck( rc, msg=subname// &
             "Error from ESMF_AlarmIsRinging" )
          ! --- make sure the datestop will always stop with dates >= stop_date
          if (trim(alarmname) == trim(seq_timemgr_alarm_datestop)) then
             call ESMF_ClockGet(EClock, CurrTime = ETime1, rc=rc)
             call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_ClockGet CurrTime" )
             call ESMF_AlarmGet(EAlarm_list(n), RingTime = ETime2, rc=rc)
             call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_AlarmGet RingTime" )
             if (ETime1 >= ETime2) seq_timemgr_alarmIsOn = .true.
          endif
       endif
    enddo

    if (.not.found) then
      write(logunit,*) subname//': ERROR alarm not valid for EClock '//trim(alarmname)
      call shr_sys_abort( subname//'ERROR: alarm invalid '//trim(alarmname) )
    endif
#ifdef USE_ESMF_LIB
    deallocate(EAlarm_list)
#endif

end function seq_timemgr_alarmIsOn

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_restartAlarmIsOn -- check if an alarm is ringing
!   
! !DESCRIPTION:
!   
!     check if an alarm is ringing
!      
! !INTERFACE: ------------------------------------------------------------------

logical function seq_timemgr_restartAlarmIsOn( EClock)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock) , intent(IN) :: EClock     ! clock/alarm

    !----- local -----
    integer :: rc
    character(len=*), parameter :: subname = '(seq_timemgr_restartAlarmIsOn) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    seq_timemgr_restartAlarmIsOn =  &
       seq_timemgr_alarmIsOn(EClock, alarmname=seq_timemgr_alarm_restart)

end function seq_timemgr_restartAlarmIsOn

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_stopAlarmIsOn -- check if an alarm is ringing
!   
! !DESCRIPTION:
!   
!     check if an alarm is ringing
!      
! !INTERFACE: ------------------------------------------------------------------

logical function seq_timemgr_stopAlarmIsOn( EClock)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock) , intent(IN) :: EClock     ! clock/alarm

    !----- local -----
    character(len=*), parameter :: subname = '(seq_timemgr_stopAlarmIsOn) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    seq_timemgr_stopAlarmIsOn =  &
       seq_timemgr_alarmIsOn(EClock, alarmname=seq_timemgr_alarm_stop)

end function seq_timemgr_stopAlarmIsOn

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_historyAlarmIsOn -- check if an alarm is ringing
!   
! !DESCRIPTION:
!   
!     check if an alarm is ringing
!      
! !INTERFACE: ------------------------------------------------------------------

logical function seq_timemgr_historyAlarmIsOn( EClock)

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock) , intent(IN) :: EClock     ! clock/alarm

    !----- local -----
    integer :: rc
    character(len=*), parameter :: subname = '(seq_timemgr_historyAlarmIsOn) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    seq_timemgr_historyAlarmIsOn =  &
       seq_timemgr_alarmIsOn(EClock, alarmname=seq_timemgr_alarm_history)

end function seq_timemgr_historyAlarmIsOn

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_ETimeInit -- Create ESMF_Time object based on YMD values
!   
! !DESCRIPTION:
!   
!     Create the ESMF_Time object corresponding to the given input time, given in
!  YMD (Year Month Day) and TOD (Time-of-day) format.
! Set the time by an integer as YYYYMMDD and integer seconds in the day
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_timemgr_ETimeInit( ETime, ymd, tod, desc )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Time) , intent(inout) :: ETime     ! Time
   integer         , intent(in)    :: ymd       ! Year, month, day YYYYMMDD
   integer         , intent(in), optional    :: tod       ! Time of day in seconds
   character(len=*), intent(in), optional    :: desc      ! Description of time to set

   !----- local -----
   character(len=*), parameter :: subname = '(seq_timemgr_ETimeInit) '
   integer :: yr, mon, day          ! Year, month, day as integers
   integer :: ltod                  ! local tod
   character(SHR_KIND_CL)  :: ldesc ! local desc
   integer :: rc                    ! return code

!-------------------------------------------------------------------------------
! Notes:
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

   call ESMF_TimeSet( ETime, yy=yr, mm=mon, dd=day, s=ltod, &
                      calendar=seq_timemgr_cal, rc=rc )
   call seq_timemgr_ESMFCodeCheck(rc, subname//': error return from '// &
                                 'ESMF_TimeSet: setting '//trim(ldesc))

end subroutine seq_timemgr_ETimeInit

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_ETimeGet -- Get the date in YYYYMMDD from from ESMF Time
!   
! !DESCRIPTION:
!   
!     Get the date in YYYYMMDD format from a ESMF time object.
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_timemgr_ETimeGet( ETime, offset, ymd, tod )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Time),   intent(IN)  :: ETime   ! Input ESMF time
    integer, optional, intent(IN)  :: offset  ! Offset from input time (sec)
    integer, optional, intent(OUT) :: ymd     ! date of day
    integer, optional, intent(OUT) :: tod     ! Time of day

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
! Notes:
!-------------------------------------------------------------------------------

    ETimeAdd = ETime
    if ( present(offset) )then
        if ( offset > 0 )then
           call ESMF_TimeIntervalSet( ETimeOff, s=offset, rc=rc )
           call seq_timemgr_ESMFCodeCheck( rc, msg=subname// &
                                          ": Error from ESMF_TimeIntervalSet" )
           ETimeAdd = ETime + ETimeOff
        else if ( offset < 0 )then
           call ESMF_TimeIntervalSet( ETimeOff, s=-offset, rc=rc )
           call seq_timemgr_ESMFCodeCheck( rc, msg=subname// &
                                          ": Error from ESMF_TimeIntervalSet" )
           ETimeAdd = ETime - ETimeOff
        end if
    end if

    call ESMF_TimeGet( ETimeAdd, yy=year, mm=month, dd=day, s=sec, rc=rc )
    call seq_timemgr_ESMFCodeCheck( rc, msg=subname// &
                                   ": Error from ESMF_TimeGet" )

    ! shr_cal has restrictions and then "stops", so override that

    if ( present(ymd) ) then
          call shr_cal_ymd2date(year,month,day,ymd)
    endif
    if ( present(tod) ) then
          tod = sec
    endif

end subroutine seq_timemgr_ETimeGet

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_EClockInit -- Initialize the ESMF clock in the shared clock
!   
! !DESCRIPTION:
!   
! Private method:
!
! Setup the ESMF clock inside the wrapped CCSM clock
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_timemgr_EClockInit( TimeStep, StartTime, RefTime, CurrTime, EClock )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_TimeInterval), intent(IN)  :: TimeStep    ! Time-step of clock
    type(ESMF_Time)        , intent(IN)  :: StartTime   ! Start time
    type(ESMF_Time)        , intent(IN)  :: RefTime     ! Reference time
    type(ESMF_Time)        , intent(IN)  :: CurrTime    ! Current time
    type(ESMF_Clock)       , intent(OUT) :: EClock      ! Output ESMF clock

    !----- local -----
    character(len=*), parameter :: subname = '(seq_timemgr_EClockInit) '
    integer :: rc                             ! ESMF return code
    character(len=SHR_KIND_CL) :: description ! Description of this clock
    type(ESMF_Time) :: clocktime              ! Current time

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    description = 'CCSM shared Time-manager clock'

    ! ------ Create ESMF Clock with input characteristics -------------------
    ! --- NOTE: StopTime is required in interface but not used, so use  -----
    ! ---       something arbitrary.  Stop handled via alarm            -----

    call seq_timemgr_ETimeInit(clocktime,  99990101, 0, "artificial stop date")

    EClock = ESMF_ClockCreate(name=trim(description), &
       TimeStep=TimeStep, startTime=StartTime, &
       refTime=RefTime, stopTime=clocktime, rc=rc)
    call seq_timemgr_ESMFCodeCheck( rc, msg=subname//': Error from ESMF_ClockCreate')

    ! ------ Advance clock to the current time (in case of a restart) -------
    call ESMF_ClockGet(EClock, currTime=clocktime, rc=rc )
    call seq_timemgr_ESMFCodeCheck(rc, subname//': Error from ESMF_ClockGet')
    do while( clocktime < CurrTime)
       call ESMF_ClockAdvance( EClock, rc=rc )
       call seq_timemgr_ESMFCodeCheck(rc, subname//': Error from ESMF_ClockAdvance')
       call ESMF_ClockGet( EClock, currTime=clocktime )
       call seq_timemgr_ESMFCodeCheck(rc, subname//': Error from ESMF_ClockGet')
    end do

    if (clocktime /= CurrTime) then
         if (loglevel > 0) write(logunit,*) trim(subname), &
            ' : WARNING clocktime and currtime inconsistent'
    endif

end subroutine seq_timemgr_EClockInit

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_EClockDateInSync -- Check that input date in sync with clock
!   
! !DESCRIPTION:
!   
!     Check that the given input date/time is in sync with clock time
!      
! !INTERFACE: ------------------------------------------------------------------

logical function seq_timemgr_EClockDateInSync( EClock, ymd, tod, prev)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock), intent(IN) :: Eclock   ! Input clock to compare
   integer,          intent(IN) :: ymd     ! Date (YYYYMMDD)
   integer,          intent(IN) :: tod     ! Time of day (sec)
   logical, optional,intent(IN) :: prev    ! If should get previous time

   !----- local -----
   character(len=*), parameter :: subname = "(seq_timemgr_EClockDateInSync) "
   type(ESMF_Time) :: ETime
   integer :: ymd1      ! Date (YYYYMMDD)
   integer :: tod1      ! Time of day
   logical :: previous  ! If need to get previous time for comparison
   integer :: rc        ! error code

!-------------------------------------------------------------------------------
! Notes:
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
!===============================================================================
! !IROUTINE: seq_timemgr_clockPrint -- Print clock information out
!   
! !DESCRIPTION:
!   
!      Print clock information out.
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_timemgr_clockPrint( SyncClock )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(seq_timemgr_type), intent(in) :: SyncClock   ! Input clock to print

   character(len=*), parameter :: subname = "(seq_timemgr_clockPrint) "
   integer(SHR_KIND_IN) :: m,n
   integer(SHR_KIND_IN) :: curr_ymd   ! Current date YYYYMMDD
   integer(SHR_KIND_IN) :: curr_tod   ! Current time of day (s)
   integer(SHR_KIND_IN) :: StepNo     ! Number of steps taken
   integer(SHR_KIND_IN) :: start_ymd  ! Starting date YYYYMMDD
   integer(SHR_KIND_IN) :: start_tod  ! Starting time-of-day (s)
   integer(SHR_KIND_IN) :: stop_ymd   ! Stop date YYYYMMDD
   integer(SHR_KIND_IN) :: stop_tod   ! Stop time-of-day (s)
   integer(SHR_KIND_IN) :: ref_ymd    ! Reference date YYYYMMDD
   integer(SHR_KIND_IN) :: ref_tod    ! Reference time-of-day (s)
   integer(SHR_KIND_IN) :: DTime      ! Time-step (seconds)
   integer(SHR_KIND_IN) :: prev_ymd   ! Prev restart alarm date (YYYYMMDD)
   integer(SHR_KIND_IN) :: prev_tod   ! Prev restart alarm time-of-day (sec)
   integer(SHR_KIND_IN) :: next_ymd   ! Next restart alarm date (YYYYMMDD)
   integer(SHR_KIND_IN) :: next_tod   ! Next restart alarm time-of-day (sec)
   integer(SHR_KIND_IN) :: IntSec     ! Alarm interval for seconds
   integer(SHR_KIND_IN) :: IntMon     ! Alarm interval for months
   integer(SHR_KIND_IN) :: IntYrs     ! Alarm interval for years
   integer(SHR_KIND_IN) :: AlarmCount ! Number of valid alarms
   character(len=64)    :: alarmname  ! Alarm name
   character(len=*), parameter :: xalarm = 'unset'
   type(ESMF_Alarm),pointer :: EAlarm_list(:)   ! EAlarm list associated with EClock
   integer(SHR_KIND_IN) :: rc         ! error code

   character(len=*), parameter ::  F06 = "(2A,L3)"
   character(len=*), parameter ::  F07 = "(3A)"
   character(len=*), parameter ::  F08 = "(2A,I8.8,3x,I5.5)"
   character(len=*), parameter ::  F09 = "(2A,2I8,I12)"
   character(len=*), parameter ::  F10 = "(2A,I2,2x,A)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  if (loglevel <= 0) return

  write(logunit,F07) subname,'calendar      = ', trim(seq_timemgr_calendar)
  write(logunit,F06) subname,'end_restart   = ', seq_timemgr_end_restart
  write(logunit,F07) ''

  do n = 1,max_clocks
     call seq_timemgr_EClockGetData( SyncClock%ECP(n)%EClock, curr_ymd=curr_ymd, &
                                curr_tod=curr_tod, start_ymd=start_ymd,    &
                                start_tod=start_tod, StepNo=StepNo,            &
                                ref_ymd=ref_ymd, ref_tod=ref_tod,              &
                                stop_ymd=stop_ymd, stop_tod=stop_tod,          &
                                dtime = dtime, alarmcount=AlarmCount)
#ifdef USE_ESMF_LIB
     allocate(EAlarm_list(AlarmCount))
     call ESMF_ClockGetAlarmList(SyncClock%ECP(n)%EClock, alarmListFlag=ESMF_ALARMLIST_ALL, &
        alarmList=EAlarm_list, alarmCount=AlarmCount, rc=rc)
#else
     call ESMF_ClockGetAlarmList(SyncClock%ECP(n)%EClock, EAlarm_list, rc=rc)
#endif
     call seq_timemgr_ESMFCodeCheck( rc, msg=subname//"Error from ESMF_ClockGetAlarmList" )

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
#ifdef USE_ESMF_LIB
     deallocate(EAlarm_list)
#endif
  enddo

end subroutine seq_timemgr_clockPrint

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_ESMFDebug -- Print ESMF stuff for debugging
!   
! !DESCRIPTION:
!   
! Print ESMF stuff for debugging
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_timemgr_ESMFDebug( EClock, ETime, ETimeInterval, istring )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock), optional, intent(in) :: EClock    ! ESMF Clock
   type(ESMF_Time) , optional, intent(inout) :: ETime     ! ESMF Time
   type(ESMF_TimeInterval), optional, intent(inout) :: ETimeInterval  ! ESMF Time Interval
   character(len=*), optional, intent(in) :: istring

   !----- local -----
   character(len=*), parameter :: subname = '(seq_timemgr_ESMFDebug) '
   character(len=128) :: timestring
   integer :: yy,mm,dd,s            ! ymds
   type(ESMF_Time) :: LTime
   type(ESMF_TimeInterval) :: LTimeInterval
   integer(SHR_KIND_I8) :: LStep
   integer :: rc                    ! return code

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
      call ESMF_ClockGet( EClock, StartTime=LTime )
      call ESMF_TimeGet(LTime, yy=yy,mm=mm,dd=dd,s=s,timestring=timestring,rc=rc)
      write(logunit,*) subname,rc,'start ymds=',yy,mm,dd,s,trim(timestring)
      call ESMF_ClockGet( EClock, CurrTime=LTime )
      call ESMF_TimeGet(LTime, yy=yy,mm=mm,dd=dd,s=s,timestring=timestring,rc=rc)
      write(logunit,*) subname,rc,'curr ymds=',yy,mm,dd,s,trim(timestring)
      call ESMF_ClockGet( EClock, StopTime=LTime )
      call ESMF_TimeGet(LTime, yy=yy,mm=mm,dd=dd,s=s,timestring=timestring,rc=rc)
      write(logunit,*) subname,rc,'stop ymds=',yy,mm,dd,s,trim(timestring)
      call ESMF_ClockGet( EClock, PrevTime=LTime )
      call ESMF_TimeGet(LTime, yy=yy,mm=mm,dd=dd,s=s,timestring=timestring,rc=rc)
      write(logunit,*) subname,rc,'prev ymds=',yy,mm,dd,s,trim(timestring)
      call ESMF_ClockGet( EClock, RefTime=LTime )
      call ESMF_TimeGet(LTime, yy=yy,mm=mm,dd=dd,s=s,timestring=timestring,rc=rc)
      write(logunit,*) subname,rc,'ref ymds=',yy,mm,dd,s,trim(timestring)
      call ESMF_ClockGet( EClock, TimeStep=LTimeInterval )
      call ESMF_TimeIntervalGet(LTimeInterval, yy=yy,mm=mm,d=dd,s=s,timestring=timestring,rc=rc)
      write(logunit,*) subname,rc,'tint ymds=',yy,mm,dd,s,trim(timestring)
      call ESMF_ClockGet( EClock, AdvanceCount=LStep )
      write(logunit,*) subname,rc,'advcnt =',LStep
   endif

end subroutine seq_timemgr_ESMFDebug

!===============================================================================
!===============================================================================
! !IROUTINE: seq_timemgr_ESMFCodeCheck -- Check return-code from ESMF -- abort if not
!   
! !DESCRIPTION:
!   
!     Check ESMF return code and abort if not successful.
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_timemgr_ESMFCodeCheck( rc, msg )

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer, intent(in)          :: rc   ! return code from ESMF
   character(len=*),optional,intent(in) :: msg  ! error message

   character(len=*),parameter :: subname = 'seq_timemgr_ESMFCodeCheck'
!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if ( rc == ESMF_SUCCESS ) return
   if ( present(msg)) then
      write(logunit,*) trim(subname),' error= ',rc,trim(msg)
   else
      write(logunit,*) trim(subname),' error= ',rc
   endif
   call shr_sys_flush(logunit)
   call shr_sys_abort(trim(subname))

end subroutine seq_timemgr_ESMFCodeCheck

!===============================================================================
!===============================================================================

end module seq_timemgr_mod
