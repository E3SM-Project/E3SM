module lnd_comp_mct

  ! !USES:

  use esmf
  use mct_mod
  use perf_mod
  use seq_cdata_mod   , only: seq_cdata, seq_cdata_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_putdata, seq_infodata_getdata
  use seq_comm_mct    , only: seq_comm_inst, seq_comm_name, seq_comm_suffix
  use shr_kind_mod    , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_strdata_mod , only: shr_strdata_type
  use shr_file_mod    , only: shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel
  use shr_file_mod    , only: shr_file_setlogunit, shr_file_setloglevel, shr_file_setio
  use shr_file_mod    , only: shr_file_freeunit
  use dlnd_comp_mod   , only: dlnd_comp_init, dlnd_comp_run, dlnd_comp_final
  use dlnd_shr_mod    , only: dlnd_shr_read_namelists
  use seq_flds_mod    , only: seq_flds_x2l_fields, seq_flds_l2x_fields

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: lnd_init_mct
  public :: lnd_run_mct
  public :: lnd_final_mct

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type(shr_strdata_type) :: SDLND
  integer(IN)            :: mpicom              ! mpi communicator
  integer(IN)            :: my_task             ! my task in mpi communicator mpicom
  integer                :: inst_index          ! number of current instance (ie. 1)
  character(len=16)      :: inst_name           ! fullname of current instance (ie. "lnd_0001")
  character(len=16)      :: inst_suffix         ! char string associated with instance (ie. "_0001" or "")
  integer(IN)            :: logunit             ! logging unit number
  integer(IN)            :: compid              ! mct comp id

  character(*), parameter :: F00   = "('(dlnd_comp_init) ',8a)"
  integer(IN) , parameter :: master_task=0 ! task number of master task
  character(*), parameter :: subName = "(lnd_init_mct) "

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine lnd_init_mct( EClock, cdata, x2l, l2x, NLFilename )

    ! !DESCRIPTION:  initialize dlnd model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2l, l2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    integer           :: phase                     ! phase of method
    logical           :: lnd_present               ! flag
    logical           :: lnd_prognostic            ! flag
    integer(IN)       :: shrlogunit                ! original log unit
    integer(IN)       :: shrloglev                 ! original log level
    logical           :: read_restart              ! start from restart
    integer(IN)       :: ierr                      ! error code
    logical           :: scmMode = .false.         ! single column mode
    real(R8)          :: scmLat  = shr_const_SPVAL ! single column lat
    real(R8)          :: scmLon  = shr_const_SPVAL ! single column lon
    character(*), parameter :: subName = "(lnd_init_mct) "
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
       call shr_file_setIO('lnd_modelio.nml'//trim(inst_suffix),logUnit)
    else
       logUnit = 6
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    !----------------------------------------------------------------------------
    ! Read dlnd namelists and set prognostic, present flags in infodata
    !----------------------------------------------------------------------------

    call t_startf('dlnd_readnml')

    call dlnd_shr_read_namelists(mpicom, my_task, master_task, &
         inst_index, inst_suffix, inst_name,  &
         logunit, shrlogunit, SDLND, lnd_present, lnd_prognostic)

    call seq_infodata_PutData(infodata, &
         lnd_present=lnd_present, &
         lnd_prognostic=lnd_prognostic)

    call t_stopf('dlnd_readnml')

    !----------------------------------------------------------------------------
    ! RETURN if present flag is false
    !----------------------------------------------------------------------------

    if (.not. lnd_present) then
       RETURN
    end if

    ! NOTE: the following will never be called if lnd_present is .false.

    !----------------------------------------------------------------------------
    ! Initialize dlnd
    !----------------------------------------------------------------------------

    call dlnd_comp_init(Eclock, x2l, l2x, &
         seq_flds_x2l_fields, seq_flds_l2x_fields, &
         SDLND, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, inst_name, logunit, read_restart, &
         scmMode, scmlat, scmlon)

    !----------------------------------------------------------------------------
    ! Fill infodata that needs to be returned from dlnd
    !----------------------------------------------------------------------------

    call seq_infodata_PutData(infodata, &
         lnd_nx=SDLND%nxg, &
         lnd_ny=SDLND%nyg )

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'dlnd_comp_init done'
    call shr_sys_flush(logunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine lnd_init_mct

  !===============================================================================
  subroutine lnd_run_mct( EClock, cdata, x2l, l2x)

    ! !DESCRIPTION: run method for dlnd model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2l
    type(mct_aVect)             ,intent(inout) :: l2x

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    integer(IN)                      :: shrlogunit   ! original log unit
    integer(IN)                      :: shrloglev    ! original log level
    character(CL)                    :: case_name    ! case name
    character(*), parameter :: subName = "(lnd_run_mct) "
    !-------------------------------------------------------------------------------

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    call seq_cdata_setptrs(cdata, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    call seq_infodata_GetData(infodata, case_name=case_name)

    call dlnd_comp_run(EClock, x2l, l2x, &
         SDLND, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, logunit, case_name)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine lnd_run_mct

  !===============================================================================
  subroutine lnd_final_mct(EClock, cdata, x2l, l2x)

    ! !DESCRIPTION: finalize method for dlnd model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2l
    type(mct_aVect)             ,intent(inout) :: l2x

    !--- formats ---
    character(*), parameter :: subName = "(lnd_final_mct) "
    !-------------------------------------------------------------------------------

    call dlnd_comp_final(my_task, master_task, logunit)

  end subroutine lnd_final_mct
  !===============================================================================

end module lnd_comp_mct
