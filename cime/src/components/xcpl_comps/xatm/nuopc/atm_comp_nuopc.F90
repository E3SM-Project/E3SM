module atm_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for XATM
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC             , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC             , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC_Model       , only : model_routine_SS        => SetServices
  use NUOPC_Model       , only : model_label_Advance     => label_Advance
  use NUOPC_Model       , only : model_label_SetRunClock => label_SetRunClock
  use NUOPC_Model       , only : model_label_Finalize    => label_Finalize
  use NUOPC_Model       , only : NUOPC_ModelGet
  use shr_sys_mod       , only : shr_sys_abort
  use shr_kind_mod      , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_file_mod      , only : shr_file_getlogunit, shr_file_setlogunit
  use dead_methods_mod  , only : chkerr, state_setscalar,  state_diagnose, alarmInit, memcheck
  use dead_methods_mod  , only : set_component_logging, get_component_instance, log_clock_advance
  use dead_nuopc_mod    , only : dead_grid_lat, dead_grid_lon, dead_grid_index
  use dead_nuopc_mod    , only : dead_init_nuopc, dead_final_nuopc, dead_meshinit
  use dead_nuopc_mod    , only : fld_list_add, fld_list_realize, fldsMax, fld_list_type
  use dead_nuopc_mod    , only : ModelInitPhase, ModelSetRunClock

  implicit none
  private ! except

  public :: SetServices

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(len=CL)      :: flds_scalar_name = ''
  integer                :: flds_scalar_num = 0
  integer                :: flds_scalar_index_nx = 0
  integer                :: flds_scalar_index_ny = 0
  integer                :: flds_scalar_index_nextsw_cday = 0

  integer                :: fldsToAtm_num = 0
  integer                :: fldsFrAtm_num = 0
  type (fld_list_type)   :: fldsToAtm(fldsMax)
  type (fld_list_type)   :: fldsFrAtm(fldsMax)
  integer, parameter     :: gridTofieldMap = 2 ! ungridded dimension is innermost

  real(r8), pointer      :: gbuf(:,:)   ! model info
  real(r8), pointer      :: lat(:)
  real(r8), pointer      :: lon(:)
  integer , allocatable  :: gindex(:)
  integer                :: nxg         ! global dim i-direction
  integer                :: nyg         ! global dim j-direction
  integer                :: my_task     ! my task in mpi communicator mpicom
  integer                :: inst_index  ! number of current instance (ie. 1)
  character(len=5)       :: inst_suffix ! char string associated with instance (ie. "_0001" or "")
  integer                :: logunit     ! logging unit number
  logical                :: mastertask
  integer                :: dbug = 1
  character(*),parameter :: modName =  "(xatm_comp_nuopc)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=ModelInitPhase, phase=0, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, phaseLabelList=(/"IPDv01p1"/), &
         userRoutine=InitializeAdvertise, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, phaseLabelList=(/"IPDv01p3"/), &
         userRoutine=InitializeRealize, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, specRoutine=ModelAdvance, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, specRoutine=ModelSetRunClock, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, specRoutine=ModelFinalize, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)     :: vm
    character(CS)     :: stdname
    integer           :: n
    integer           :: lsize       ! local array size
    integer           :: shrlogunit  ! original log unit
    character(CL)     :: cvalue
    character(len=CL) :: logmsg
    logical           :: isPresent, isSet
    character(len=*),parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=rc)

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localpet=my_task, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    mastertask = (my_task==0)

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! set logunit and set shr logging to my log file
    !----------------------------------------------------------------------------

    call set_component_logging(gcomp, mastertask, logunit, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! Initialize xatm
    !----------------------------------------------------------------------------

    call dead_init_nuopc('atm', inst_suffix, logunit, lsize, gbuf, nxg, nyg)

    allocate(gindex(lsize))
    allocate(lon(lsize))
    allocate(lat(lsize))

    gindex(:) = gbuf(:,dead_grid_index)
    lat(:)    = gbuf(:,dead_grid_lat)
    lon(:)    = gbuf(:,dead_grid_lon)

    !--------------------------------
    ! advertise import and export fields
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldName')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldCount')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNX')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNY')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nextsw_cday
       write(logmsg,*) flds_scalar_index_nextsw_cday
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nextsw_cday = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxNextSwCday')
    endif

    if (nxg /= 0 .and. nyg /= 0) then

       call fld_list_add(fldsFrAtm_num, fldsFrAtm, trim(flds_scalar_name))
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_topo'       )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_z'          )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_u'          )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_v'          )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_tbot'       )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_ptem'       )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_shum'       )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_pbot'       )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_dens'       )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_pslv'       )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_rainc'    )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_rainl'    )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_snowc'    )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_snowl'    )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_lwdn'     )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swndr'    )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swvdr'    )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swndf'    )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swvdf'    )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swnet'    )
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_bcph'  , ungridded_lbound=1, ungridded_ubound=3)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_ocph'  , ungridded_lbound=1, ungridded_ubound=3)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet', ungridded_lbound=1, ungridded_ubound=4)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry', ungridded_lbound=1, ungridded_ubound=4)

       call fld_list_add(fldsToAtm_num, fldsToAtm, trim(flds_scalar_name))
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_anidr'  )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_avsdf'  )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_anidf'  )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_avsdr'  )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sl_lfrac'  )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Si_ifrac'  )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'So_ofrac'  )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_tref'   )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_qref'   )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_t'      )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'So_t'      )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sl_fv'     )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sl_ram1'   )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sl_snowh'  )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Si_snowh'  )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'So_ssq'    )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'So_re'     )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_u10'    )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_taux' )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_tauy' )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_lat'  )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_sen'  )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_lwup' )
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_evap' )

       do n = 1,fldsFrAtm_num
          if(mastertask) write(logunit,*)'Advertising From Xatm ',trim(fldsFrAtm(n)%stdname)
          call NUOPC_Advertise(exportState, standardName=fldsFrAtm(n)%stdname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end do

       do n = 1,fldsToAtm_num
          if(mastertask) write(logunit,*)'Advertising To Xatm',trim(fldsToAtm(n)%stdname)
          call NUOPC_Advertise(importState, standardName=fldsToAtm(n)%stdname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       enddo
    end if

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=rc)

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! input/output arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR) :: convCIM, purpComp
    type(ESMF_Mesh)        :: Emesh
    type(ESMF_Time)        :: nextTime
    real(r8)               :: nextsw_cday
    integer                :: n
    integer                :: shrlogunit                ! original log unit
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize: xatm) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=rc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logUnit)

    !--------------------------------
    ! generate the mesh
    !--------------------------------

    call dead_meshinit(gcomp, nxg, nyg, gindex, lon, lat, Emesh, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! realize the actively coupled fields, now that a mesh is established
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    !--------------------------------

    call fld_list_realize( &
         state=exportState, &
         fldList=fldsFrAtm, &
         numflds=fldsFrAtm_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':datmExport',&
         mesh=Emesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call fld_list_realize( &
         state=importState, &
         fldList=fldsToAtm, &
         numflds=fldsToAtm_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':datmImport',&
         mesh=Emesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Pack export state
    !--------------------------------

    call state_setexport(exportState, rc=rc)

    call State_SetScalar(dble(nxg),flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call State_SetScalar(dble(nyg),flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Set time of next radiation computation

    call ESMF_ClockGetNextTime(clock, nextTime)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(nextTime, dayOfYear_r8=nextsw_cday)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call State_diagnose(exportState,subname//':ES',rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"
    call ESMF_AttributeAdd(comp,  convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ShortName", "XATM", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", "Atmosphere Dead Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
         "The dead models stand in as test model for active components." // &
         "Coupling data is artificially generated ", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2017", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Sea Ice", convention=convCIM, purpose=purpComp, rc=rc)
#endif

    call shr_file_setLogUnit (shrlogunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=rc)

  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)  :: clock
    type(ESMF_State)  :: exportState
    real(r8)          :: nextsw_cday
    integer           :: shrlogunit ! original log unit
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 1) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=rc)
    end if
    call memcheck(subname, 3, mastertask)

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logunit)

    !--------------------------------
    ! Pack export state
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, exportState=exportState, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, &
          flds_scalar_name, flds_scalar_num, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call state_diagnose(exportState,subname//':ES',rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (mastertask) then
          call log_clock_advance(clock, 'XATM', logunit, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       endif
    endif

    call shr_file_setLogUnit (shrlogunit)

    if (dbug > 5) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=rc)
    end if

  end subroutine ModelAdvance

  !===============================================================================

  subroutine state_setexport(exportState, rc)

    ! input/output variables
    type(ESMF_State)  , intent(inout) :: exportState
    integer, intent(out) :: rc

    ! local variables
    integer :: nf, nind
    !--------------------------------------------------

    rc = ESMF_SUCCESS

    ! Start from index 2 in order to Skip the scalar field here  
    do nf = 2,fldsFrAtm_num
       if (fldsFrAtm(nf)%ungridded_ubound == 0) then
          call field_setexport(exportState, trim(fldsFrAtm(nf)%stdname), lon, lat, nf=nf, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          do nind = 1,fldsFrAtm(nf)%ungridded_ubound
             call field_setexport(exportState, trim(fldsFrAtm(nf)%stdname), lon, lat, nf=nf+nind-1, &
                  ungridded_index=nind, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
          end do
       end if
    end do

  end subroutine state_setexport

  !===============================================================================

  subroutine field_setexport(exportState, fldname, lon, lat, nf, ungridded_index, rc)

    use shr_const_mod , only : pi=>shr_const_pi

    ! intput/otuput variables
    type(ESMF_State)  , intent(inout) :: exportState
    character(len=*)  , intent(in)    :: fldname
    real(r8)          , intent(in)    :: lon(:)
    real(r8)          , intent(in)    :: lat(:)
    integer           , intent(in)    :: nf
    integer, optional , intent(in)    :: ungridded_index
    integer           , intent(out)   :: rc

    ! local variables
    integer           :: i, ncomp
    type(ESMF_Field)  :: lfield
    real(r8), pointer :: data1d(:)
    real(r8), pointer :: data2d(:,:)
    !--------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(exportState, itemName=trim(fldname), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ncomp = 1
    if (present(ungridded_index)) then
       call ESMF_FieldGet(lfield, farrayPtr=data2d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (gridToFieldMap == 1) then
          do i = 1,size(data2d, dim=1)
             data2d(i,ungridded_index) = (nf*100) * cos(pi*lat(i)/180.0_R8) * &
                  sin((pi*lon(i)/180.0_R8) - (ncomp-1)*(pi/3.0_R8) ) + (ncomp*10.0_R8)
          end do
       else if (gridToFieldMap == 2) then
          do i = 1,size(data2d, dim=2)
             data2d(ungridded_index,i) = (nf*100) * cos(pi*lat(i)/180.0_R8) * &
                  sin((pi*lon(i)/180.0_R8) - (ncomp-1)*(pi/3.0_R8) ) + (ncomp*10.0_R8)
          end do
       end if
    else
       call ESMF_FieldGet(lfield, farrayPtr=data1d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       do i = 1,size(data1d)
          data1d(i) = (nf*100) * cos(pi*lat(i)/180.0_R8) * &
               sin((pi*lon(i)/180.0_R8) - (ncomp-1)*(pi/3.0_R8) ) + (ncomp*10.0_R8)
       end do
    end if

  end subroutine field_setexport

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=rc)

    call dead_final_nuopc('atm', logunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=rc)

  end subroutine ModelFinalize

end module atm_comp_nuopc
