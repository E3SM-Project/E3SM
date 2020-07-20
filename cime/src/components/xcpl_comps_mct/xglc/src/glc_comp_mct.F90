module glc_comp_mct

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
  use seq_flds_mod    , only: seq_flds_g2x_fields, seq_flds_x2g_fields

  ! !PUBLIC TYPES:
  implicit none
  save
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: glc_init_mct
  public :: glc_run_mct
  public :: glc_final_mct

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------
  integer(IN)            :: mpicom              ! mpi communicator
  integer(IN)            :: my_task             ! my task in mpi communicator mpicom
  integer                :: inst_index          ! number of current instance (ie. 1)
  character(len=16)      :: inst_name           ! fullname of current instance (ie. "glc_0001")
  character(len=16)      :: inst_suffix = ""    ! char string associated with instance (ie. "_0001" or "")
  integer(IN)            :: logunit             ! logging unit number
  integer(IN)            :: compid              ! mct comp id
  real(r8) ,  pointer    :: gbuf(:,:)           ! model grid
  integer(IN),parameter  :: master_task=0       ! task number of master task

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine glc_init_mct( EClock, cdata, x2d, d2x, NLFilename )

    ! !DESCRIPTION: initialize dead glc model

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
    logical                          :: glc_present    ! if true, component is present
    logical                          :: glc_prognostic ! if true, component is prognostic
    logical                          :: glcocn_present
    logical                          :: glcice_present
    logical                          :: glclnd_present
    !-------------------------------------------------------------------------------

    ! Set cdata pointers to derived types (in coupler)
    call seq_cdata_setptrs(cdata, &
         id=compid, &
         mpicom=mpicom, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    ! Obtain infodata variables
    call seq_infodata_getData(infodata, glc_phase=phase)

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
          call shr_file_setIO('glc_modelio.nml'//trim(inst_suffix),logUnit)
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
    ! Initialize xglc
    !----------------------------------------------------------------------------

    call dead_init_mct('glc', Eclock, x2d, d2x, &
         seq_flds_x2g_fields, seq_flds_g2x_fields, &
         gsmap, ggrid, gbuf, mpicom, compid, my_task, master_task, &
         inst_index, inst_suffix, inst_name, logunit, nxg, nyg)

    if (nxg == 0 .and. nyg == 0) then
       glc_present = .false.
       glc_prognostic = .false.
       glclnd_present = .false.
       glcocn_present = .false.
       glcice_present = .false.
    else
       glc_present = .true.
       glc_prognostic = .true.
       glclnd_present = .true.
       glcocn_present = .false.
       glcice_present = .false.
    end if

    call seq_infodata_PutData( infodata, dead_comps=.true., &
         glc_present=glc_present, &
         glc_prognostic=glc_prognostic, &
         glclnd_present=glclnd_present, &
         glcocn_present=glcocn_present, &
         glcice_present=glcice_present, &
         glc_nx=nxg, glc_ny=nyg)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logunit)

  end subroutine glc_init_mct

  !===============================================================================
  subroutine glc_run_mct(EClock, cdata, x2d, d2x)

    ! !DESCRIPTION: run method for dead glc model

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
    character(*), parameter :: subName = "(glc_run_mct) "
    !-------------------------------------------------------------------------------

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    call seq_cdata_setptrs(cdata, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    call dead_run_mct('glc', EClock, x2d, d2x, &
       gsmap, ggrid, gbuf, mpicom, compid, my_task, master_task, logunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine glc_run_mct

  !===============================================================================
  subroutine glc_final_mct(EClock, cdata, x2d, d2x)

    implicit none

    ! !DESCRIPTION: finalize method for dead model

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

    !--- formats ---
    character(*), parameter :: subName = "(glc_final_mct) "
    !-------------------------------------------------------------------------------

    call dead_final_mct('glc', my_task, master_task, logunit)

  end subroutine glc_final_mct
  !===============================================================================

end module glc_comp_mct
