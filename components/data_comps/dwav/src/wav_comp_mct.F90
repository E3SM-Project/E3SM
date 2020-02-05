module wav_comp_mct

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
  use dwav_comp_mod   , only: dwav_comp_init, dwav_comp_run, dwav_comp_final
  use dwav_shr_mod    , only: dwav_shr_read_namelists
  use seq_flds_mod    , only: seq_flds_w2x_fields, seq_flds_x2w_fields

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: wav_init_mct
  public :: wav_run_mct
  public :: wav_final_mct

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type(shr_strdata_type) :: SDWAV
  integer(IN)            :: mpicom              ! mpi communicator
  integer(IN)            :: my_task             ! my task in mpi communicator mpicom
  integer                :: inst_index          ! number of current instance (ie. 1)
  character(len=16)      :: inst_name           ! fullname of current instance (ie. "wav_0001")
  character(len=16)      :: inst_suffix         ! char string associated with instance (ie. "_0001" or "")
  integer(IN)            :: logunit             ! logging unit number
  integer(IN)            :: compid              ! mct comp id

  character(*), parameter :: F00   = "('(dwav_comp_init) ',8a)"
  integer(IN) , parameter :: master_task=0 ! task number of master task
  character(*), parameter :: subName = "(wav_init_mct) "

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine wav_init_mct( EClock, cdata, x2w, w2x, NLFilename )

    ! !DESCRIPTION:  initialize dwav model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2w, w2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    integer           :: phase                     ! phase of method
    logical           :: wav_present               ! flag
    logical           :: wav_prognostic            ! flag
    integer(IN)       :: shrlogunit                ! original log unit
    integer(IN)       :: shrloglev                 ! original log level
    logical           :: read_restart              ! start from restart
    integer(IN)       :: ierr                      ! error code
    logical           :: scmMode = .false.         ! single column mode
    real(R8)          :: scmLat  = shr_const_SPVAL ! single column lat
    real(R8)          :: scmLon  = shr_const_SPVAL ! single column lon
    logical           :: post_assim = .false.      ! Run is post-DA
    character(*), parameter :: subName = "(wav_init_mct) "
    !-------------------------------------------------------------------------------

    ! Set cdata pointers
    call seq_cdata_setptrs(cdata, &
         id=compid,               &
         mpicom=mpicom,           &
         gsMap=gsmap,             &
         dom=ggrid,               &
         infodata=infodata,       &
         post_assimilation=post_assim)

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
       call shr_file_setIO('wav_modelio.nml'//trim(inst_suffix),logUnit)
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
    ! Read dwav namelists and set prognostic, present flags in infodata
    !----------------------------------------------------------------------------

    call t_startf('dwav_readnml')

    call dwav_shr_read_namelists(mpicom, my_task, master_task, &
         inst_index, inst_suffix, inst_name,  &
         logunit, shrlogunit, SDWAV, wav_present, wav_prognostic)

    call seq_infodata_PutData(infodata, &
         wav_present=wav_present, &
         wav_prognostic=wav_prognostic)

    call t_stopf('dwav_readnml')

    !----------------------------------------------------------------------------
    ! RETURN if present flag is false
    !----------------------------------------------------------------------------

    if (.not. wav_present) then
       RETURN
    end if

    ! NOTE: the following will never be called if wav_present is .false.

    ! Diagnostic print statement to test DATA_ASSIMILATION_WAV XML variable
    !   usage (and therefore a proxy for other component types).
    if (my_task == master_task) then
       if (post_assim) then
          write(logunit, *) subName//': Post data assimilation signal'
       else if (read_restart) then
          write(logunit, *) subName//': Restart run'
       else
          write(logunit, *) subName//': Initial run'
       end if
       call shr_sys_flush(logunit)
    end if

    !----------------------------------------------------------------------------
    ! Initialize dwav
    !----------------------------------------------------------------------------

    call dwav_comp_init(Eclock, x2w, w2x, &
         SDWAV, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, inst_name, logunit, read_restart, &
         seq_flds_w2x_fields, seq_flds_x2w_fields)

    !----------------------------------------------------------------------------
    ! Fill infodata that needs to be returned from dwav
    !----------------------------------------------------------------------------

    call seq_infodata_PutData(infodata, &
         wav_nx=SDWAV%nxg, &
         wav_ny=SDWAV%nyg )

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'dwav_comp_init done'
    call shr_sys_flush(logunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine wav_init_mct

  !===============================================================================
  subroutine wav_run_mct( EClock, cdata, x2w, w2x)

    ! !DESCRIPTION: run method for dwav model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2w
    type(mct_aVect)             ,intent(inout) :: w2x

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    integer(IN)                      :: shrlogunit   ! original log unit
    integer(IN)                      :: shrloglev    ! original log level
    logical                          :: read_restart ! start from restart
    character(CL)                    :: case_name    ! case name
    character(*), parameter :: subName = "(wav_run_mct) "
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

    call dwav_comp_run(EClock, x2w, w2x, &
         SDWAV, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, logunit, case_name)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine wav_run_mct

  !===============================================================================
  subroutine wav_final_mct(EClock, cdata, x2w, w2x)

    ! !DESCRIPTION: finalize method for dwav model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2w
    type(mct_aVect)             ,intent(inout) :: w2x

    !--- formats ---
    character(*), parameter :: subName = "(wav_final_mct) "
    !-------------------------------------------------------------------------------

    call dwav_comp_final(my_task, master_task, logunit)

  end subroutine wav_final_mct
  !===============================================================================

end module wav_comp_mct
