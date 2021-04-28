module ocn_comp_mct

  ! !USES:

  use esmf
  use mct_mod
  use perf_mod
  use seq_cdata_mod   , only: seq_cdata, seq_cdata_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_putdata, seq_infodata_getdata
  use seq_comm_mct    , only: seq_comm_inst, seq_comm_name, seq_comm_suffix
  use seq_timemgr_mod , only: seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockGetData
  use shr_kind_mod    , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_strdata_mod , only: shr_strdata_type
  use shr_file_mod    , only: shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel
  use shr_file_mod    , only: shr_file_setlogunit, shr_file_setloglevel, shr_file_setio
  use shr_file_mod    , only: shr_file_freeunit
  use docn_comp_mod   , only: docn_comp_init, docn_comp_run, docn_comp_final
  use docn_shr_mod    , only: docn_shr_read_namelists
  use seq_flds_mod    , only: seq_flds_x2o_fields, seq_flds_o2x_fields

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type(shr_strdata_type)  :: SDOCN
  integer(IN)             :: mpicom              ! mpi communicator
  integer(IN)             :: my_task             ! my task in mpi communicator mpicom
  integer                 :: inst_index          ! number of current instance (ie. 1)
  character(len=16)       :: inst_name           ! fullname of current instance (ie. "lnd_0001")
  character(len=16)       :: inst_suffix         ! char string associated with instance (ie. "_0001" or "")
  integer(IN)             :: logunit             ! logging unit number
  integer(IN)             :: compid              ! mct comp id
  logical                 :: read_restart        ! start from restart
  integer(IN) , parameter :: master_task=0 ! task number of master task
  integer     , parameter :: dbug = 10

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine ocn_init_mct( EClock, cdata, x2o, o2x, NLFilename )

    ! !DESCRIPTION:  initialize docn model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2o, o2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    integer           :: phase                     ! phase of method
    logical           :: ocn_present               ! flag
    logical           :: ocn_prognostic            ! flag
    logical           :: ocnrof_prognostic         ! flag
    integer(IN)       :: shrlogunit                ! original log unit
    integer(IN)       :: shrloglev                 ! original log level
    integer(IN)       :: ierr                      ! error code
    logical           :: scmMode = .false.         ! single column mode
    logical           :: iop_mode = .false.        ! IOP mode
    real(R8)          :: scmLat  = shr_const_SPVAL ! single column lat
    real(R8)          :: scmLon  = shr_const_SPVAL ! single column lon
    character(*), parameter :: F00   = "('(docn_comp_init) ',8a)"
    character(*), parameter :: subName = "(ocn_init_mct) "
    !-------------------------------------------------------------------------------

    ! Set cdata pointers
    call seq_cdata_setptrs(cdata, &
         id=compid, &
         mpicom=mpicom, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    ! Obtain infodata variables
    call seq_infodata_getData(infodata, &
         single_column=scmMode, &
         iop_mode=iop_mode, &
         scmlat=scmlat, scmlon=scmLon, &
         read_restart=read_restart)

    ! Determine instance information
    inst_name   = seq_comm_name(compid)
    inst_index  = seq_comm_inst(compid)
    inst_suffix = seq_comm_suffix(compid)

    ! Determine communicator group
    call mpi_comm_rank(mpicom, my_task, ierr)

    !--- open log file ---
    if (my_task == master_task) then
       logUnit = shr_file_getUnit()
       call shr_file_setIO('ocn_modelio.nml'//trim(inst_suffix),logUnit)
    else
       logUnit = 6
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (logUnit)

    !----------------------------------------------------------------------------
    ! Read docn namelists and set prognostic, present flags in infodata
    !----------------------------------------------------------------------------

    call t_startf('docn_readnml')

    call docn_shr_read_namelists(mpicom, my_task, master_task, &
         inst_index, inst_suffix, inst_name,  &
         logunit, shrlogunit, SDOCN, ocn_present, ocn_prognostic, ocnrof_prognostic)

    call seq_infodata_PutData(infodata, &
         ocn_present=ocn_present, &
         ocn_prognostic=ocn_prognostic, &
         ocnrof_prognostic=ocnrof_prognostic)

    call t_stopf('docn_readnml')

    !----------------------------------------------------------------------------
    ! RETURN if present flag is false
    !----------------------------------------------------------------------------

    if (.not. ocn_present) then
       RETURN
    end if

    ! NOTE: the following will never be called if ocn_present is .false.

    !----------------------------------------------------------------------------
    ! Initialize docn
    !----------------------------------------------------------------------------

    call docn_comp_init(Eclock, x2o, o2x, &
         seq_flds_x2o_fields, seq_flds_o2x_fields, &
         SDOCN, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, inst_name, logunit, read_restart, &
         scmMode, iop_mode, scmlat, scmlon)

    !----------------------------------------------------------------------------
    ! Fill infodata that needs to be returned from docn
    !----------------------------------------------------------------------------

    call seq_infodata_PutData(infodata, &
         ocn_nx=SDOCN%nxg, &
         ocn_ny=SDOCN%nyg )

    !----------------------------------------------------------------------------
    ! diagnostics
    !----------------------------------------------------------------------------

    if (dbug > 1) then
       if (my_task == master_task) then
          call mct_aVect_info(2, o2x, istr="initial diag"//':AV')
       end if
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'docn_comp_init done'
    call shr_sys_flush(logunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine ocn_init_mct

  !===============================================================================

  subroutine ocn_run_mct( EClock, cdata,  x2o, o2x)

    ! !DESCRIPTION: run method for docn model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2o        ! driver -> docn
    type(mct_aVect)             ,intent(inout) :: o2x        ! docn   -> driver

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    integer(IN)                      :: shrlogunit    ! original log unit
    integer(IN)                      :: shrloglev     ! original log level
    character(CL)                    :: case_name     ! case name
    logical                          :: write_restart ! restart alarm is ringing
    integer(IN)                      :: currentYMD    ! model date
    integer(IN)                      :: currentTOD    ! model sec into model date
    character(*), parameter :: subName = "(ocn_run_mct) "
    !-------------------------------------------------------------------------------

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (logUnit)

    call seq_cdata_setptrs(cdata, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    call seq_infodata_GetData(infodata, case_name=case_name)

    write_restart = seq_timemgr_RestartAlarmIsOn(EClock)

    ! For mct - the component clock is advance at the beginning of the time interval
    call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)

    call docn_comp_run(EClock, x2o, o2x, &
       SDOCN, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, read_restart, write_restart, &
       currentYMD, currentTOD, case_name=case_name)

    if (dbug > 1) then
       if (my_task == master_task) then
          call mct_aVect_info(2, o2x, istr="run diag"//':AV')
       end if
    endif

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

  end subroutine ocn_run_mct

  !===============================================================================
  subroutine ocn_final_mct(EClock, cdata, x2o, o2x)

    ! !DESCRIPTION: finalize method for docn model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2o
    type(mct_aVect)             ,intent(inout) :: o2x

    !--- formats ---
    character(*), parameter :: subName = "(ocn_final_mct) "
    !-------------------------------------------------------------------------------

    call docn_comp_final(my_task, master_task, logunit)

  end subroutine ocn_final_mct
  !===============================================================================

end module ocn_comp_mct
