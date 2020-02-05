module rof_comp_mct

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
  use drof_comp_mod   , only: drof_comp_init, drof_comp_run, drof_comp_final
  use drof_shr_mod    , only: drof_shr_read_namelists
  use seq_flds_mod    , only: seq_flds_x2r_fields, seq_flds_r2x_fields

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: rof_init_mct
  public :: rof_run_mct
  public :: rof_final_mct

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type(shr_strdata_type) :: SDROF
  integer(IN)            :: mpicom              ! mpi communicator
  integer(IN)            :: my_task             ! my task in mpi communicator mpicom
  integer                :: inst_index          ! number of current instance (ie. 1)
  character(len=16)      :: inst_name           ! fullname of current instance (ie. "lnd_0001")
  character(len=16)      :: inst_suffix         ! char string associated with instance (ie. "_0001" or "")
  integer(IN)            :: logunit             ! logging unit number
  integer(IN)            :: compid              ! mct comp id

  character(*), parameter :: F00   = "('(drof_comp_init) ',8a)"
  integer(IN) , parameter :: master_task=0 ! task number of master task
  character(*), parameter :: subName = "(rof_init_mct) "

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine rof_init_mct( EClock, cdata, x2r, r2x, NLFilename )

    ! !DESCRIPTION:  initialize drof model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2r, r2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    logical           :: rof_present               ! flag
    logical           :: rof_prognostic            ! flag
    logical           :: rofice_present            ! flag
    logical           :: flood_present             ! flag
    integer(IN)       :: shrlogunit                ! original log unit
    integer(IN)       :: shrloglev                 ! original log level
    logical           :: read_restart              ! start from restart
    integer(IN)       :: ierr                      ! error code
    character(*), parameter :: subName = "(rof_init_mct) "
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
       call shr_file_setIO('rof_modelio.nml'//trim(inst_suffix),logUnit)
    else
       logUnit = 6
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    call t_startf('drof_readnml')

    call drof_shr_read_namelists(mpicom, my_task, master_task, &
         inst_index, inst_suffix, inst_name,  &
         logunit, shrlogunit, SDROF, rof_present, rof_prognostic, rofice_present, flood_present)

    call seq_infodata_PutData(infodata, &
         rof_present=rof_present, &
         rof_prognostic=rof_prognostic, &
         rofice_present=rofice_present, &
         flood_present=flood_present)

    call t_stopf('drof_readnml')

    !----------------------------------------------------------------------------
    ! RETURN if present flag is false
    !----------------------------------------------------------------------------

    if (.not. rof_present) then
       RETURN
    end if

    ! NOTE: the following will never be called if rof_present is .false.

    !----------------------------------------------------------------------------
    ! Initialize drof
    !----------------------------------------------------------------------------

    call drof_comp_init(Eclock, x2r, r2x, &
         seq_flds_x2r_fields, seq_flds_r2x_fields, &
         SDROF, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, inst_name, logunit, read_restart)

    !----------------------------------------------------------------------------
    ! Fill infodata that needs to be returned from drof
    !----------------------------------------------------------------------------

    call seq_infodata_PutData(infodata, &
         rof_nx=SDROF%nxg, &
         rof_ny=SDROF%nyg )

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'drof_comp_init done'
    call shr_sys_flush(logunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine rof_init_mct

  !===============================================================================
  subroutine rof_run_mct( EClock, cdata, x2r, r2x)

    ! !DESCRIPTION: run method for drof model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2r
    type(mct_aVect)             ,intent(inout) :: r2x

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    integer(IN)                      :: shrlogunit   ! original log unit
    integer(IN)                      :: shrloglev    ! original log level
    character(CL)                    :: case_name    ! case name
    character(*), parameter :: subName = "(rof_run_mct) "
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

    call drof_comp_run(EClock, x2r, r2x, &
       SDROF, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, case_name=case_name)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine rof_run_mct

  !===============================================================================
  subroutine rof_final_mct(EClock, cdata, x2r, r2x)

    ! !DESCRIPTION: finalize method for drof model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2r
    type(mct_aVect)             ,intent(inout) :: r2x

    !--- formats ---
    character(*), parameter :: subName = "(rof_final_mct) "
    !-------------------------------------------------------------------------------

    call drof_comp_final(my_task, master_task, logunit)

  end subroutine rof_final_mct
  !===============================================================================

end module rof_comp_mct
