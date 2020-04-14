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
  use dead_mct_mod    , only: dead_init_mct, dead_run_mct, dead_final_mct
  use seq_flds_mod    , only: seq_flds_r2x_fields, seq_flds_x2r_fields

  ! !PUBLIC TYPES:
  implicit none
  save
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
  character(CS)          :: myModelName = 'rof' ! user defined model name
  integer(IN)            :: mpicom              ! mpi communicator
  integer(IN)            :: my_task             ! my task in mpi communicator mpicom
  integer                :: inst_index          ! number of current instance (ie. 1)
  character(len=16)      :: inst_name           ! fullname of current instance (ie. "lnd_0001")
  character(len=16)      :: inst_suffix = ""    ! char string associated with instance (ie. "_0001" or "")
  integer(IN)            :: logunit             ! logging unit number
  integer(IN)            :: compid              ! mct comp id
  real(r8) ,  pointer    :: gbuf(:,:)           ! model grid
  integer(IN),parameter  :: master_task=0       ! task number of master task

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine rof_init_mct( EClock, cdata, x2d, d2x, NLFilename )

    ! !DESCRIPTION: initialize dead rof model

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2d, d2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

    !--- local variables ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    integer(IN)                      :: shrlogunit     ! original log unit
    integer(IN)                      :: shrloglev      ! original log level
    integer(IN)                      :: nxg            ! global dim i-direction
    integer(IN)                      :: nyg            ! global dim j-direction
    integer(IN)                      :: phase          ! initialization phase
    integer(IN)                      :: ierr           ! error code
    logical                          :: rof_present    ! if true, component is present
    logical                          :: rof_prognostic ! if true, component is prognostic
    logical                          :: rofice_present
    logical                          :: flood_present
    !-------------------------------------------------------------------------------

    ! Set cdata pointers to derived types (in coupler)
    call seq_cdata_setptrs(cdata, &
         id=compid, &
         mpicom=mpicom, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    ! Obtain infodata variables
    call seq_infodata_getData(infodata, rof_phase=phase)

    ! Determine instance information
    inst_name   = seq_comm_name(compid)
    inst_index  = seq_comm_inst(compid)
    inst_suffix = seq_comm_suffix(compid)

    if (phase == 1) then
       ! Determine communicator group
       call mpi_comm_rank(mpicom, my_task, ierr)

       !--- open log file ---
       if (my_task == master_task) then
          logUnit = shr_file_getUnit()
          call shr_file_setIO('rof_modelio.nml'//trim(inst_suffix),logUnit)
       else
          logUnit = 6
       endif
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    !----------------------------------------------------------------------------
    ! Initialize xrof
    !----------------------------------------------------------------------------

    call dead_init_mct('rof', Eclock, x2d, d2x, &
         seq_flds_x2r_fields, seq_flds_r2x_fields, &
         gsmap, ggrid, gbuf, mpicom, compid, my_task, master_task, &
         inst_index, inst_suffix, inst_name, logunit, nxg, nyg)

    if (nxg == 0 .and. nyg == 0) then
       rof_present = .false.
       rof_prognostic = .false.
       rofice_present = .false.
       flood_present = .false.
    else
       rof_present = .true.
       rof_prognostic = .true.
       rofice_present = .false.
       flood_present = .true.
    end if

    call seq_infodata_PutData( infodata, dead_comps=.true., &
         rof_present=rof_present, &
         rof_prognostic=rof_prognostic, &
         rofice_present=rofice_present, &
         flood_present=flood_present, &
         rof_nx=nxg, rof_ny=nyg)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logunit)

  end subroutine rof_init_mct

  !===============================================================================
  subroutine rof_run_mct(EClock, cdata, x2d, d2x)

    ! !DESCRIPTION: run method for dead rof model

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock) ,intent(inout) :: EClock     ! clock
    type(seq_cdata)  ,intent(inout) :: cdata
    type(mct_aVect)  ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)  ,intent(inout) :: d2x        ! dead   -> driver

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    integer(IN)                      :: shrlogunit     ! original log unit
    integer(IN)                      :: shrloglev      ! original log level
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

    call dead_run_mct('rof', EClock, x2d, d2x, &
       gsmap, ggrid, gbuf, mpicom, compid, my_task, master_task, logunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine rof_run_mct

  !===============================================================================
  subroutine rof_final_mct(EClock, cdata, x2d, d2x)

    implicit none

    ! !DESCRIPTION: finalize method for dead model

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

    !--- formats ---
    character(*), parameter :: subName = "(rof_final_mct) "
    !-------------------------------------------------------------------------------

    call dead_final_mct('rof', my_task, master_task, logunit)

  end subroutine rof_final_mct
  !===============================================================================

end module rof_comp_mct
