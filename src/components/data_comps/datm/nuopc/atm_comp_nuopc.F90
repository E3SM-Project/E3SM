module atm_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for DATM
  !----------------------------------------------------------------------------

  use shr_kind_mod          , only : CXX => shr_kind_CXX
  use med_constants_mod     , only : IN, R8, I8
  use shr_log_mod           , only : shr_log_Unit
  use shr_cal_mod           , only : shr_cal_ymd2date, shr_cal_noleap, shr_cal_gregorian
  use shr_file_mod          , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_file_mod          , only : shr_file_getloglevel, shr_file_setloglevel
  use shr_file_mod          , only : shr_file_setIO, shr_file_getUnit
  use esmFlds               , only : flds_scalar_name
  use esmFlds               , only : flds_scalar_num
  use esmFlds               , only : flds_scalar_index_nx
  use esmFlds               , only : flds_scalar_index_ny
  use esmFlds               , only : flds_scalar_index_nextsw_cday
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Realize
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Concat
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Deactivate
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getnumflds
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getfldinfo
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_Clock_TimePrint
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_SetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_Diagnose
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_Meshinit
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_ArrayToState
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_StateToArray
  use shr_nuopc_time_mod    , only : shr_nuopc_time_alarmInit
  use shr_strdata_mod       , only : shr_strdata_type

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS        => SetServices, &
    model_label_Advance     => label_Advance, &
    model_label_SetRunClock => label_SetRunClock, &
    model_label_Finalize    => label_Finalize

  use datm_shr_mod , only: datm_shr_read_namelists
  use datm_shr_mod , only: iradsw         ! namelist input
  use datm_shr_mod , only: presaero       ! namelist input
  use datm_shr_mod , only: datm_shr_getNextRadCDay
  use datm_comp_mod, only: datm_comp_init, datm_comp_run, datm_comp_final
  use mct_mod

  implicit none
  private ! except

  public  :: SetServices

  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelSetRunClock
  private :: ModelFinalize

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type fld_list_type
    character(len=128) :: stdname
  end type fld_list_type

  integer,parameter          :: fldsMax = 100
  integer                    :: fldsToAtm_num = 0
  integer                    :: fldsFrAtm_num = 0
  type (fld_list_type)       :: fldsToAtm(fldsMax)
  type (fld_list_type)       :: fldsFrAtm(fldsMax)

  character(len=80)          :: myModelName = 'atm'       ! user defined model name
  type(shr_strdata_type)     :: SDATM
  type(mct_gsMap), target    :: gsMap_target
  type(mct_gGrid), target    :: ggrid_target
  type(mct_gsMap), pointer   :: gsMap
  type(mct_gGrid), pointer   :: ggrid
  type(mct_aVect)            :: x2d
  type(mct_aVect)            :: d2x
  integer                    :: compid                    ! mct comp id
  integer                    :: mpicom                    ! mpi communicator
  integer                    :: my_task                   ! my task in mpi communicator mpicom
  integer                    :: inst_index                ! number of current instance (ie. 1)
  character(len=16)          :: inst_name                 ! fullname of current instance (ie. "lnd_0001")
  character(len=16)          :: inst_suffix = ""          ! char string associated with instance (ie. "_0001" or "")
  integer                    :: logunit                   ! logging unit number
  integer    ,parameter      :: master_task=0             ! task number of master task
  integer                    :: localPet
  logical                    :: atm_present               ! flag
  logical                    :: atm_prognostic            ! flag
  logical                    :: unpack_import
  character(len=256)         :: case_name                 ! case name
  character(len=256)         :: tmpstr                    ! tmp string
  integer                    :: dbrc
  integer, parameter         :: dbug = 10
  character(len=*),parameter :: grid_option = "mesh"      ! grid_de, grid_arb, grid_reg, mesh
  character(len=80)          :: calendar                  ! calendar name
  character(len=CXX)         :: flds_a2x = ''
  character(len=CXX)         :: flds_x2a = ''

  !----- formats -----
  character(*),parameter :: modName =  "(atm_comp_nuopc)"
  character(*),parameter :: u_FILE_u = __FILE__

  !===============================================================================
  contains
  !===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter  :: subname=trim(modName)//':(ATMSetServices) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !----------------------------------------------------------------------------
    ! Switch to IPDv01 by filtering all other phaseMap entries
    !----------------------------------------------------------------------------

    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
         acceptStringList=(/"IPDv01p"/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)      :: vm
    integer            :: lmpicom
    character(len=256) :: cvalue
    character(len=80)  :: stdname, shortname
    logical            :: activefld
    integer            :: n,nflds
    integer            :: ierr       ! error code
    integer            :: shrlogunit ! original log unit
    integer            :: shrloglev  ! original log level
    logical            :: isPresent
    character(len=512) :: diro
    character(len=512) :: logfile
    logical            :: flds_co2a  ! use case
    logical            :: flds_co2b  ! use case
    logical            :: flds_co2c  ! use case
    logical            :: flds_wiso  ! use case
    character(len=*),parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localPet=localPet, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call mpi_comm_dup(lmpicom, mpicom, ierr)
    call mpi_comm_rank(mpicom, my_task, ierr)

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="inst_name", value=inst_name, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="inst_index", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) inst_index

    call ESMF_AttributeGet(gcomp, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", value=inst_suffix, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       inst_suffix = ''
    end if

    !----------------------------------------------------------------------------
    ! set logunit and set shr logging to my log file
    !----------------------------------------------------------------------------

    if (my_task == master_task) then
       call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       logunit = shr_file_getUnit()
       open(logunit,file=trim(diro)//"/"//trim(logfile))
    else
       logUnit = 6
    endif

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (logunit)

    !----------------------------------------------------------------------------
    ! Read input namelists and set present and prognostic flags
    !----------------------------------------------------------------------------

    call datm_shr_read_namelists(mpicom, my_task, master_task, &
         inst_index, inst_suffix, inst_name, &
         logunit, shrlogunit, SDATM, atm_present, atm_prognostic)

    ! NOTE: atm_present flag is not needed - since the run sequence
    ! will have no call to this routine for the atm_present flag being
    ! set to false (i.e. null mode) - only the atm_prognostic flag is
    ! needed below

    if (atm_prognostic) then
       unpack_import = .true.
    else
       unpack_import = .false.
    end if

    !--------------------------------
    ! determine necessary toggles for below
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2a
    call ESMF_LogWrite('flds_co2a = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2b', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2b
    call ESMF_LogWrite('flds_co2b = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_co2c
    call ESMF_LogWrite('flds_co2c = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    call NUOPC_CompAttributeGet(gcomp, name='flds_wiso', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_wiso
    call ESMF_LogWrite('flds_wiso = '// trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

    !--------------------------------
    ! advertise import and export fields
    !--------------------------------

    call fld_list_add(fldsFrAtm_num, fldsFrAtm, trim(flds_scalar_name)) 
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_topo'    , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_z'       , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_u'       , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_v'       , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_tbot'    , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_ptem'    , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_shum'    , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_pbot'    , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_dens'    , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_pslv'    , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_rainc' , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_rainl' , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_snowc' , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_snowl' , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_lwdn'  , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swndr' , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swvdr' , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swndf' , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swvdf' , flds_concat=flds_a2x)
    call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swnet' , flds_concat=flds_a2x)  ! only diagnostic
    if (presaero) then
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_bcphidry' , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_bcphodry' , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_bcphiwet' , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_ocphidry' , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_ocphodry' , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_ocphiwet' , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet1'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet2'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet3'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet4'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry1'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry2'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry3'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry4'  , flds_concat=flds_a2x)
    end if
    if (flds_co2a .or. flds_co2b .or. flds_co2c) then
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_co2prog'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_co2diag'  , flds_concat=flds_a2x)
    end if

    call fld_list_add(fldsToAtm_num, fldsToAtm, trim(flds_scalar_name))
    if (atm_prognostic) then
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_anidr'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_avsdf'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_anidf'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_avsdr'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sl_lfrac'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Si_ifrac'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'So_ofrac'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_tref'   , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_qref'   , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_t'      , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'So_t'      , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sl_fv'     , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sl_ram1'   , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sl_snowh'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Si_snowh'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'So_ssq'    , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'So_re'     , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_u10'    , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_taux' , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_tauy' , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_lat'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_sen'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_lwup' , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_evap' , flds_concat=flds_x2a)
    end if

    ! Now advertise these fields
    do n = 1,fldsFrAtm_num
      call NUOPC_Advertise(exportState, standardName=fldsFrAtm(n)%stdname, TransferOfferGeomObject='will provide', rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo
    do n = 1,fldsToAtm_num
      call NUOPC_Advertise(importState, standardName=fldsToAtm(n)%stdname, TransferOfferGeomObject='will provide', rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR)  :: convCIM, purpComp
    type(ESMF_Grid)         :: Egrid
    type(ESMF_TIME)         :: currTime
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_Mesh)         :: Emesh
    type(ESMF_Calendar)     :: esmf_calendar             ! esmf calendar
    type(ESMF_CalKind_Flag) :: esmf_caltype              ! esmf calendar type
    type(ESMF_VM)           :: vm
    integer                 :: nx_global, ny_global
    integer                 :: n
    character(len=256)      :: cvalue
    integer                 :: shrlogunit                ! original log unit
    integer                 :: shrloglev                 ! original log level
    logical                 :: read_restart              ! start from restart
    integer                 :: ierr                      ! error code
    logical                 :: scmMode = .false.         ! single column mode
    real(R8)                :: scmLat  = shr_const_SPVAL ! single column lat
    real(R8)                :: scmLon  = shr_const_SPVAL ! single column lon
    integer                 :: current_ymd               ! model date
    integer                 :: current_year              ! model year
    integer                 :: current_mon               ! model month
    integer                 :: current_day               ! model day
    integer                 :: current_tod               ! model sec into model date
    integer(I8)             :: stepno                    ! step number
    integer                 :: modeldt                   ! integer timestep
    real(R8)                :: nextsw_cday               ! calendar of next atm sw
    logical                 :: connected                 ! is field connected?
    integer                 :: lsize
    integer                 :: iam
    real(r8), pointer       :: lon(:),lat(:)
    integer , pointer       :: gindex(:)
    real(R8)                :: orbEccen                  ! orb eccentricity (unit-less)
    real(R8)                :: orbMvelpp                 ! orb moving vernal eq (radians)
    real(R8)                :: orbLambm0                 ! orb mean long of perhelion (radians)
    real(R8)                :: orbObliqr                 ! orb obliquity (radians)
    character(len=*) , parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    ! TODO: read_restart, scmlat, scmlon, orbeccen, orbmvelpp, orblambm0, orbobliqr needs to be obtained
    ! from the config attributes of the gridded component

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (logUnit)

    !--------------------------------
    ! call datm init routine
    !--------------------------------

    gsmap => gsmap_target
    ggrid => ggrid_target

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlon

    call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlat

    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmMode

    call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart

    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid

    !----------------------------------------------------------------------------
    ! Determine orbital values (these might change dynamically)
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbEccen
    call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbObliqr
    call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbLambm0
    call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbMvelpp

    !----------------------------------------------------------------------------
    ! Determine calendar info
    !----------------------------------------------------------------------------

    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, advanceCount=stepno, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, &
         calkindflag=esmf_caltype, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call ESMF_LogWrite(subname//" ERROR bad ESMF calendar name "//trim(calendar), ESMF_LOGMSG_ERROR, rc=dbrc)
       rc = ESMF_Failure
       return
    end if

    call ESMF_TimeIntervalGet( timeStep, s=modeldt, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call datm_comp_init(clock, x2d, d2x, &
         flds_x2a, flds_a2x, &
         SDATM, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, inst_name, logunit, read_restart, &
         scmMode, scmlat, scmlon, &
         orbEccen, orbMvelpp, orbLambm0, orbObliqr, &
         calendar, modeldt, current_ymd, current_tod, current_mon, atm_prognostic)

    !----------------------------------------------------------------------------
    ! Set nextsw_cday
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart

    if (read_restart) then
       nextsw_cday = datm_shr_getNextRadCDay( current_ymd, current_tod, stepno, modeldt, iradsw, calendar )
    else
       ! For a startup run the nextsw_cday is just the current calendar day
       call ESMF_TimeGet( currTime, dayofyear_r8=nextsw_cday, rc=rc )
    endif

    !--------------------------------
    ! Generate the mesh
    !--------------------------------

    nx_global = SDATM%nxg
    ny_global = SDATM%nyg
    lsize = mct_gsMap_lsize(gsMap, mpicom)
    allocate(lon(lsize))
    allocate(lat(lsize))
    allocate(gindex(lsize))

    call mpi_comm_rank(mpicom, iam, ierr)
    call mct_gGrid_exportRattr(ggrid,'lon',lon,lsize)
    call mct_gGrid_exportRattr(ggrid,'lat',lat,lsize)
    call mct_gsMap_OrderedPoints(gsMap_target, iam, gindex)

    call shr_nuopc_grid_MeshInit(gcomp, nx_global, ny_global, mpicom, gindex, lon, lat, Emesh, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(lon)
    deallocate(lat)
    deallocate(gindex)

    !--------------------------------
    ! realize the actively coupled fields, now that a mesh is established
    ! Note: the following will occur in the realize call
    ! 1) loop over all of the entries in fldsToAtm and fldsFrAtm to creates a field via the following call:
    !    field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=shortname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    ! 2) realizes the field via the following call
    !    call NUOPC_Realize(state, field=field, rc=rc)
    !    where state is either importState or exportState
    !  NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    !  by replacing the advertised fields with the newly created fields of the same name.
    !--------------------------------

    call fld_list_realize( &
         state=ExportState, &
         fldList=fldsFrAtm, &
         numflds=fldsFrAtm_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':datmExport',&
         mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call fld_list_realize( &
         state=importState, &
         fldList=fldsToAtm, &
         numflds=fldsToAtm_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':datmImport',&
         mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Pack export state
    ! Copy from d2x to exportState
    ! Set the coupling scalars
    !--------------------------------

    call shr_nuopc_grid_ArrayToState(d2x%rattr, flds_a2x, exportState, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !TODO: do we still need the nx and ny scalars?
    call shr_nuopc_methods_State_SetScalar(dble(nx_global),flds_scalar_index_nx, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(ny_global),flds_scalar_index_ny, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       if (my_task == master_task) then
          call mct_aVect_info(2, d2x, istr='initial diag'//':AV')
       end if
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"
    call ESMF_AttributeAdd(comp, convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ShortName", "DATM", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", "Climatological Atmosphere Data Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
         "The CIME data models perform the basic function of " // &
         "reading external data, modifying that data, and then " // &
         "sending it to the driver via coupling " // &
         "interfaces. The driver and other models have no " // &
         "fundamental knowledge of whether another component " // &
         "is fully active or just a data model.  In some cases, " // &
         "data models are prognostic and also receive and use " // &
         "some data sent by the driver to the data model.  But " // &
         "in most cases, the data models are not running " // &
         "prognostically and have no need to receive any data " // &
         "from the driver.", &
         convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2010", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Atmosphere", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Name", "TBD", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "EmailAddress", "TBD", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", convention=convCIM, purpose=purpComp, rc=rc)
#endif

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)


  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)        :: clock
    type(ESMF_Time)         :: time
    type(ESMF_State)        :: importState, exportState
    type(ESMF_Alarm)        :: alarm
    integer                 :: shrlogunit    ! original log unit
    integer                 :: shrloglev     ! original log level
    real(r8)                :: nextsw_cday
    logical                 :: write_restart ! restart alarm is ringing
    integer                 :: nextymd       ! model date
    integer                 :: nexttod       ! model sec into model date
    integer                 :: yr            ! year
    integer                 :: mon           ! month
    integer                 :: day           ! day in month
    type(ESMF_Time)         :: currTime
    type(ESMF_Time)         :: nextTime
    type(ESMF_TimeInterval) :: timeStep
    integer(I8)             :: stepno        ! step number
    integer                 :: modeldt       ! model timestep
    real(R8)                :: orbEccen      ! orb eccentricity (unit-less)
    real(R8)                :: orbMvelpp     ! orb moving vernal eq (radians)
    real(R8)                :: orbLambm0     ! orb mean long of perhelion (radians)
    real(R8)                :: orbObliqr     ! orb obliquity (radians)
    character(len=256)      :: cvalue
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !--------------------------------
    ! Reset shr logging to my log file
    !--------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (logunit)

    !--------------------------------
    ! query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 1) then
       if (my_task == master_task) then
          call shr_nuopc_methods_Clock_TimePrint(clock,subname//'clock',rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    endif

    !--------------------------------
    ! Unpack import state
    !--------------------------------

    if (unpack_import) then
       call shr_nuopc_grid_StateToArray(importState, x2d%rattr, flds_x2a, grid_option, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !--------------------------------
    ! Run model
    !--------------------------------

    ! Get orbital parameters (these can be changed by the mediator)
    ! TODO: need to put in capability for these to be modified for variable orbitals
    call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbEccen
    call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbObliqr
    call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbLambm0
    call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbMvelpp

    ! Determine if need to write restarts

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       write_restart = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       write_restart = .false.
    endif

    ! For nuopc - the component clock is advanced at the end of the time interval
    ! For these to match for now - need to advance nuopc one timestep ahead for
    ! shr_strdata time interpolation

    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, advanceCount=stepno, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    nextTime = currTime + timeStep
    call ESMF_TimeGet( nextTime, yy=yr, mm=mon, dd=day, s=nexttod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr, mon, day, nextymd)

    call ESMF_TimeIntervalGet( timeStep, s=modeldt, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Advance the model

    call datm_comp_run(clock, &
         x2a=x2d, &
         a2x=d2x, &
         SDATM=SDATM, &
         gsmap=gsmap, &
         ggrid=ggrid, &
         mpicom=mpicom, &
         compid=compid, &
         my_task=my_task, &
         master_task=master_task, &
         inst_suffix=inst_suffix, &
         logunit=logunit, &
         orbEccen=orbEccen, &
         orbMvelpp=orbMvelpp, &
         orbLambm0=orbLambm0, &
         orbObliqr=orbObliqr, &
         write_restart=write_restart, &
         target_ymd=nextYMD, &
         target_tod=nextTOD, &
         target_mon=mon, &
         calendar=calendar, &
         modeldt=modeldt, &
         case_name=case_name, &
         atm_prognostic=atm_prognostic)

    ! Use nextYMD and nextTOD here since since the component - clock is advance at the END of the time interval
    nextsw_cday = datm_shr_getNextRadCDay( nextYMD, nextTOD, stepno, modeldt, iradsw, calendar )

    !--------------------------------
    ! Pack export state
    !--------------------------------

    call shr_nuopc_grid_ArrayToState(d2x%rattr, flds_a2x, exportState, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       if (my_task == master_task) then
          call mct_aVect_info(2, d2x, istr='run diag'//':AV', pe=localPet)
       end if
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    if (my_task == master_task) then
       call ESMF_ClockPrint(clock, options="currTime", preString="------>Advancing ATM from: ", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockPrint(clock, options="stopTime", preString="--------------------------------> to: ", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option       ! Restart option units
    integer                  :: restart_n            ! Number until restart interval
    integer                  :: restart_ymd          ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    integer                  :: first_time = .true.
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------                                                                                 
    ! set restart alarm
    !--------------------------------                                                                                 

    if (first_time) then
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call shr_nuopc_time_alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       first_time = .false.
    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    implicit none
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call datm_comp_final(my_task, master_task, logunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelFinalize

  !===============================================================================

  subroutine fld_list_add(num, fldlist, stdname, flds_concat)
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname
    character(len=*), optional, intent(inout) :: flds_concat

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(atm_comp_nuopc:fld_list_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
      return
    endif
    fldlist(num)%stdname        = trim(stdname)

    if (present(flds_concat)) then
       if (len_trim(flds_concat) + len_trim(stdname) + 1 >= len(flds_concat)) then
          call ESMF_LogWrite(subname//': ERROR: max len of flds_concat has been exceeded', &
               ESMF_LOGMSG_ERROR, line=__LINE__, file= u_FILE_u, rc=dbrc)
       end if
       if (trim(flds_concat) == '') then
          flds_concat = trim(stdname)
       else
          flds_concat = trim(flds_concat)//':'//trim(stdname)
       end if
    end if

  end subroutine fld_list_add

  !===============================================================================

  subroutine fld_list_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove 
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(atm_comp_nuopc:fld_list_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the field
             field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fld_list_realize

end module atm_comp_nuopc
