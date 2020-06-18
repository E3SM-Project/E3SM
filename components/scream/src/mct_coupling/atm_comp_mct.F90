module atm_comp_mct

  ! !USES:
  use scream_scorpio_interface, only: eam_init_pio_1, eam_init_pio_2, eam_history_finalize, eam_history_write
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
  use seq_flds_mod    , only: seq_flds_a2x_fields, seq_flds_x2a_fields
  use seq_timemgr_mod , only: seq_timemgr_EClockGetData
  use iso_c_binding   , only: c_int, c_double

  ! !PUBLIC TYPES:
  implicit none
  save
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: atm_init_mct
  public :: atm_run_mct
  public :: atm_final_mct

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------
  integer                :: mpicom_atm          ! mpi communicator
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
  subroutine atm_init_mct( EClock, cdata, x2d, d2x, NLFilename )
    use iso_c_binding, only: C_CHAR, C_NULL_CHAR
    !
    ! F90<->CXX interfaces
    !
    interface
      subroutine scream_init(f_comm,start_ymd,start_tod,yaml_fname) bind(c)
        use iso_c_binding, only: c_int, C_CHAR
        !
        ! Arguments
        !
        integer (kind=c_int),   intent(in) :: start_tod, start_ymd, f_comm
        character(kind=C_CHAR), target, intent(in) :: yaml_fname(*)
      end subroutine scream_init
    end interface

    ! !DESCRIPTION: initialize dead atm model

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
    logical                          :: atm_present    ! if true, component is present
    logical                          :: atm_prognostic ! if true, component is prognostic
    integer (IN)                     :: start_tod, start_ymd
    character(kind=C_CHAR,len=80)    :: yaml_fname = "data/scream_input.yaml"
    !-------------------------------------------------------------------------------

    ! Set cdata pointers to derived types (in coupler)
    call seq_cdata_setptrs(cdata, &
         id=compid, &
         mpicom=mpicom_atm, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    ! Obtain infodata variables
    call seq_infodata_getData(infodata, atm_phase=phase)
    if (phase > 1) RETURN

    ! Determine instance information
    inst_name   = seq_comm_name(compid)
    inst_index  = seq_comm_inst(compid)
    inst_suffix = seq_comm_suffix(compid)

    if (phase == 1) then
       ! Determine communicator group
       call mpi_comm_rank(mpicom_atm, my_task, ierr)

       !--- open log file ---
       if (my_task == master_task) then
          logUnit = shr_file_getUnit()
          call shr_file_setIO('atm_modelio.nml'//trim(inst_suffix),logUnit)
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
    ! Initialize xatm
    !----------------------------------------------------------------------------

    call dead_init_mct('atm', Eclock, x2d, d2x, &
         seq_flds_x2a_fields, seq_flds_a2x_fields, &
         gsmap, ggrid, gbuf, mpicom_atm, compid, my_task, master_task, &
         inst_index, inst_suffix, inst_name, logunit, nxg, nyg)

    call seq_timemgr_EClockGetData(EClock, start_ymd=start_ymd, start_tod=start_tod)
    call scream_init (mpicom_atm, INT(start_ymd, KIND=c_int), INT(start_tod, KIND=c_int),TRIM(yaml_fname)//C_NULL_CHAR)
    if (nxg == 0 .and. nyg == 0) then
       atm_present = .false.
       atm_prognostic = .false.
    else
       atm_present = .true.
       atm_prognostic = .true.
    end if

    call seq_infodata_PutData( infodata, dead_comps=.true., &
         atm_present=atm_present, &
         atm_prognostic=atm_prognostic, &
         atm_nx=nxg, atm_ny=nyg)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logunit)

    !----------------------------------------------------------------------------
    ! Initialize pio
    !----------------------------------------------------------------------------
    call eam_init_pio_1(mpicom_atm,compid)
    call eam_init_pio_2()

  end subroutine atm_init_mct

  !===============================================================================
  subroutine atm_run_mct(EClock, cdata, x2d, d2x)

    !
    ! F90<->CXX interfaces
    !
    interface
      subroutine scream_run (dt) bind(c)
        use iso_c_binding, only: c_double
        !
        ! Arguments
        !
        real(kind=c_double), intent(in) :: dt
      end subroutine scream_run
    end interface
    ! !DESCRIPTION: run method for dead atm model

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
    real(R8)                         :: nextsw_cday    ! calendar of next atm sw
    character(*), parameter :: subName = "(atm_run_mct) "
    real(kind=c_double)              :: dt_scream
    !-------------------------------------------------------------------------------

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    call seq_cdata_setptrs(cdata, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    call dead_run_mct('atm', EClock, x2d, d2x, &
       gsmap, ggrid, gbuf, mpicom_atm, compid, my_task, master_task, logunit)
    dt_scream = 300.0
    call scream_run( dt_scream )

    ! Set time of next radiadtion computation
    call seq_timemgr_EClockGetData (EClock, next_cday=nextsw_cday)
    call seq_infodata_PutData(infodata, nextsw_cday=nextsw_cday)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)
    !----------------------------------------------------------------------------
    ! Run pio
    !----------------------------------------------------------------------------
    call eam_history_write()

  end subroutine atm_run_mct

  !===============================================================================
  subroutine atm_final_mct(EClock, cdata, x2d, d2x)
    implicit none
    !
    ! F90<->CXX interfaces
    !
    interface
      subroutine scream_finalize () bind(c)
      end subroutine scream_finalize
    end interface


    ! !DESCRIPTION: finalize method for dead model

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

    !--- formats ---
    character(*), parameter :: subName = "(atm_final_mct) "
    !-------------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    ! Finalize pio
    !----------------------------------------------------------------------------
    call eam_history_finalize()
    !----------------------------------------------------------------------------
    ! Finish the rest of ATM model
    !----------------------------------------------------------------------------
    call dead_final_mct('atm', my_task, master_task, logunit)
    call scream_finalize()

  end subroutine atm_final_mct
  !===============================================================================

end module atm_comp_mct
