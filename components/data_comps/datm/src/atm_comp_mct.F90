module atm_comp_mct

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
  use datm_comp_mod   , only: datm_comp_init, datm_comp_run, datm_comp_final
  use datm_shr_mod    , only: datm_shr_read_namelists
  use datm_shr_mod    , only: presaero
  use seq_flds_mod    , only: seq_flds_a2x_fields, seq_flds_x2a_fields

  ! !PUBLIC TYPES:
  implicit none
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
  type(shr_strdata_type) :: SDATM
  integer(IN)            :: mpicom              ! mpi communicator
  integer(IN)            :: my_task             ! my task in mpi communicator mpicom
  integer                :: inst_index          ! number of current instance (ie. 1)
  character(len=16)      :: inst_name           ! fullname of current instance (ie. "lnd_0001")
  character(len=16)      :: inst_suffix = ""    ! char string associated with instance (ie. "_0001" or "")
  integer(IN)            :: logunit             ! logging unit number
  integer(IN)            :: compid              ! mct comp id
  integer(IN),parameter  :: master_task=0       ! task number of master task

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine atm_init_mct( EClock, cdata, x2a, a2x, NLFilename )

    implicit none

    ! !DESCRIPTION: initialize data atm model

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2a, a2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

    !--- local variables ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    integer           :: phase                     ! phase of method
    logical           :: atm_present               ! flag
    logical           :: atm_prognostic            ! flag
    integer(IN)       :: shrlogunit                ! original log unit
    integer(IN)       :: shrloglev                 ! original log level
    logical           :: read_restart              ! start from restart
    integer(IN)       :: ierr                      ! error code
    logical           :: scmMode = .false.         ! single column mode
    real(R8)          :: scmLat  = shr_const_SPVAL ! single column lat
    real(R8)          :: scmLon  = shr_const_SPVAL ! single column lon
    real(R8)          :: orbEccen                  ! orb eccentricity (unit-less)
    real(R8)          :: orbMvelpp                 ! orb moving vernal eq (radians)
    real(R8)          :: orbLambm0                 ! orb mean long of perhelion (radians)
    real(R8)          :: orbObliqr                 ! orb obliquity (radians)
    real(R8)          :: nextsw_cday               ! calendar of next atm sw

    !--- formats ---
    character(*), parameter :: F00   = "('(datm_comp_init) ',8a)"
    integer(IN) , parameter :: master_task=0 ! task number of master task
    character(*), parameter :: subName = "(atm_init_mct) "
    !-------------------------------------------------------------------------------

    ! Set cdata pointers to derived types (in coupler)
    call seq_cdata_setptrs(cdata, &
         id=compid, &
         mpicom=mpicom, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    ! Obtain infodata variables
    call seq_infodata_getData(infodata,&
         atm_phase=phase, &
         single_column=scmMode, &
         scmlat=scmlat, &
         scmlon=scmLon, &
         orb_eccen=orbEccen,&
         orb_mvelpp=orbMvelpp, &
         orb_lambm0=orbLambm0,&
         orb_obliqr=orbObliqr, &
         read_restart=read_restart)

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
    ! Read input namelists and set present and prognostic flags
    !----------------------------------------------------------------------------

    if (phase == 1) then
       call t_startf('datm_readnml')
       call datm_shr_read_namelists(mpicom, my_task, master_task, &
            inst_index, inst_suffix, inst_name, &
            logunit, shrlogunit, SDATM, atm_present, atm_prognostic)

       call seq_infodata_PutData(infodata, &
            atm_present=atm_present, &
            atm_prognostic=atm_prognostic)
       call t_stopf('datm_readnml')
    end if

    !----------------------------------------------------------------------------
    ! RETURN if present flag is false
    !----------------------------------------------------------------------------

    if (phase == 1) then
       if (.not. atm_present) then
          RETURN
       end if
    end if

    ! NOTE: the following will never be called if atm_present is .false.

    !----------------------------------------------------------------------------
    ! Initialize datm
    !----------------------------------------------------------------------------

    call datm_comp_init(Eclock, x2a, a2x, &
         seq_flds_x2a_fields, seq_flds_a2x_fields, &
         SDATM, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, inst_name, logunit, read_restart, &
         scmMode, scmlat, scmlon, &
         orbEccen, orbMvelpp, orbLambm0, orbObliqr, phase, nextsw_cday)

    !----------------------------------------------------------------------------
    ! Fill infodata that needs to be returned from datm
    !----------------------------------------------------------------------------

    if (phase == 1) then
       call seq_infodata_PutData(infodata, &
            atm_nx=SDATM%nxg, &
            atm_ny=SDATM%nyg, &
            atm_aero=presaero, &
            nextsw_cday=nextsw_cday )
    else
       call seq_infodata_PutData(infodata, &
            nextsw_cday=nextsw_cday )
    end if

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'datm_comp_init done'
    call shr_sys_flush(logunit)
    call shr_sys_flush(shrlogunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine atm_init_mct

  !===============================================================================
  subroutine atm_run_mct( EClock, cdata,  x2a, a2x)

    ! !DESCRIPTION: run method for datm model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2a        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: a2x        ! dead   -> driver

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    integer(IN)                      :: shrlogunit     ! original log unit
    integer(IN)                      :: shrloglev      ! original log level
    character(CL)                    :: case_name      ! case name
    real(R8)                         :: orbEccen       ! orb eccentricity (unit-less)
    real(R8)                         :: orbMvelpp      ! orb moving vernal eq (radians)
    real(R8)                         :: orbLambm0      ! orb mean long of perhelion (radians)
    real(R8)                         :: orbObliqr      ! orb obliquity (radians)
    real(R8)                         :: nextsw_cday    ! calendar of next atm sw
    character(*), parameter :: subName = "(atm_run_mct) "
    !-------------------------------------------------------------------------------

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    call seq_cdata_setptrs(cdata, &
         gsMap=gsmap, &
         dom=ggrid, &
         infodata=infodata)

    call seq_infodata_GetData(infodata, &
         case_name=case_name, &
         orb_eccen=orbEccen, &
         orb_mvelpp=orbMvelpp, &
         orb_lambm0=orbLambm0, &
         orb_obliqr=orbObliqr)

    call datm_comp_run( &
         EClock = EClock, &
         x2a = x2a, &
         a2x = a2x, &
         SDATM = SDATM, &
         gsmap = gsmap, &
         ggrid = ggrid, &
         mpicom = mpicom, &
         compid = compid, &
         my_task = my_task, &
         master_task = master_task, &
         inst_suffix = inst_suffix, &
         logunit = logunit, &
         orbEccen = orbEccen, &
         orbMvelpp = orbMvelpp, &
         orbLambm0 = orbLambm0, &
         orbObliqr = orbObliqr, &
         nextsw_cday = nextsw_cday, &
         case_name = case_name)

    call seq_infodata_PutData(infodata, nextsw_cday=nextsw_cday )

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine atm_run_mct

  !===============================================================================
  subroutine atm_final_mct(EClock, cdata, x2a, a2x)

    ! !DESCRIPTION: finalize method for dead atm model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2a
    type(mct_aVect)             ,intent(inout) :: a2x

    !--- formats ---
    character(*), parameter :: subName = "(atm_final_mct) "
    !-------------------------------------------------------------------------------

    call  datm_comp_final(my_task, master_task, logunit)

  end subroutine atm_final_mct
  !===============================================================================

end module atm_comp_mct
