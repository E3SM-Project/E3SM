module ice_comp_mct

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
  use dice_comp_mod   , only: dice_comp_init, dice_comp_run, dice_comp_final
  use dice_shr_mod    , only: dice_shr_read_namelists
  use seq_flds_mod    , only: seq_flds_i2x_fields, seq_flds_x2i_fields, seq_flds_i2o_per_cat
#ifdef HAVE_MOAB
  use seq_comm_mct, only : MPSIID !            iMOAB app id for ice
  use iso_c_binding
  use iMOAB           , only: iMOAB_RegisterApplication
#endif
  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: ice_init_mct
  public :: ice_run_mct
  public :: ice_final_mct

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------
  type(shr_strdata_type) :: SDICE
  integer(IN)            :: mpicom              ! mpi communicator
  integer(IN)            :: my_task             ! my task in mpi communicator mpicom
  integer                :: inst_index          ! number of current instance (ie. 1)
  character(len=16)      :: inst_name           ! fullname of current instance (ie. "lnd_0001")
  character(len=16)      :: inst_suffix         ! char string associated with instance (ie. "_0001" or "")
  integer(IN)            :: logunit             ! logging unit number
  integer(IN)            :: compid              ! mct comp id
  logical                :: read_restart        ! start from restart

  character(*), parameter :: F00   = "('(dice_comp_init) ',8a)"
  integer(IN) , parameter :: master_task=0 ! task number of master task
  character(*), parameter :: subName = "(ice_init_mct) "
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  subroutine ice_init_mct( EClock, cdata, x2i, i2x, NLFilename )
#ifdef HAVE_MOAB
    use shr_stream_mod, only: shr_stream_getDomainInfo, shr_stream_getFile
#endif
    ! !DESCRIPTION: initialize dice model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2i, i2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    logical           :: ice_present               ! flag
    logical           :: ice_prognostic            ! flag
    integer(IN)       :: shrlogunit                ! original log unit
    integer(IN)       :: shrloglev                 ! original log level
    integer(IN)       :: ierr                      ! error code
    logical           :: scmMode = .false.         ! single column mode
    real(R8)          :: scmLat  = shr_const_SPVAL ! single column lat
    real(R8)          :: scmLon  = shr_const_SPVAL ! single column lon
    character(*), parameter :: subName = "(ice_init_mct) "
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
       call shr_file_setIO('ice_modelio.nml'//trim(inst_suffix),logUnit)
    else
       logUnit = 6
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logUnit)

    call t_startf('dice_readnml')

    call dice_shr_read_namelists(mpicom, my_task, master_task, &
         inst_index, inst_suffix, inst_name,  &
         logunit, shrlogunit, SDICE, ice_present, ice_prognostic)

    call seq_infodata_PutData(infodata, &
         ice_present=ice_present, &
         ice_prognostic=ice_prognostic, &
         iceberg_prognostic=.false.)

    call t_stopf('dice_readnml')

    !----------------------------------------------------------------------------
    ! RETURN if present flag is false
    !----------------------------------------------------------------------------

    if (.not. ice_present) then
       RETURN
    end if

    ! NOTE: the following will never be called if ice_present is .false.

    !----------------------------------------------------------------------------
    ! Initialize dice
    !----------------------------------------------------------------------------
#ifdef HAVE_MOAB
    ierr = iMOAB_RegisterApplication(trim("DICE")//C_NULL_CHAR, mpicom, compid, MPSIID)
    if (ierr .ne. 0) then
      write(logunit,*) subname,' error in registering data ice comp'
      call shr_sys_abort(subname//' ERROR in registering data ice comp')
    endif
#endif
    call dice_comp_init(Eclock, x2i, i2x, &
         seq_flds_x2i_fields, seq_flds_i2x_fields, seq_flds_i2o_per_cat, &
         SDICE, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, inst_name, logunit, read_restart, &
         scmMode, scmlat, scmlon)


#ifdef HAVE_MOAB
    if (my_task == master_task) then
       ! send path of ice domain to MOAB coupler.
       write(logunit,*), ' file used for ice domain ', SDICE%domainFile
       call seq_infodata_PutData( infodata, ice_domain=SDICE%domainFile)
    endif
#endif
    !----------------------------------------------------------------------------
    ! Fill infodata that needs to be returned from dice
    !----------------------------------------------------------------------------

    call seq_infodata_PutData(infodata, &
         ice_nx=SDICE%nxg, &
         ice_ny=SDICE%nyg )

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    if (my_task == master_task) write(logunit,F00) 'dice_comp_init done'
    call shr_sys_flush(logunit)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine ice_init_mct

  !===============================================================================
  subroutine ice_run_mct( EClock, cdata,  x2i, i2x)

    ! !DESCRIPTION:  run method for dice model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2i        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: i2x        ! dead   -> driver

    !--- local ---
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap)        , pointer :: gsMap
    type(mct_gGrid)        , pointer :: ggrid
    integer(IN)                      :: shrlogunit   ! original log unit
    integer(IN)                      :: shrloglev    ! original log level
    character(CL)                    :: case_name    ! case name
    character(*), parameter :: subName = "(ice_run_mct) "
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

    call dice_comp_run(EClock, x2i, i2x, &
         seq_flds_i2o_per_cat, &
         SDICE, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, logunit, read_restart, case_name)

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunit)

  end subroutine ice_run_mct

  !===============================================================================
  subroutine ice_final_mct(EClock, cdata, x2d, d2x)

    ! !DESCRIPTION: finalize method for dead ice model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

    !--- formats ---
    character(*), parameter :: subName = "(ice_final_mct) "
    !-------------------------------------------------------------------------------

    call dice_comp_final(my_task, master_task, logunit)

  end subroutine ice_final_mct
  !===============================================================================

end module ice_comp_mct
