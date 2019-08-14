module glc_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for XGLC
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC            , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC            , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC_Model      , only : model_routine_SS        => SetServices
  use NUOPC_Model      , only : model_label_Advance     => label_Advance
  use NUOPC_Model      , only : model_label_SetRunClock => label_SetRunClock
  use NUOPC_Model      , only : model_label_Finalize    => label_Finalize
  use NUOPC_Model      , only : NUOPC_ModelGet
  use shr_sys_mod      , only : shr_sys_abort
  use shr_kind_mod     , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_file_mod     , only : shr_file_getlogunit, shr_file_setlogunit
  use dead_methods_mod , only : chkerr, state_setscalar,  state_diagnose, alarmInit, memcheck
  use dead_methods_mod , only : set_component_logging, get_component_instance, log_clock_advance
  use dead_nuopc_mod   , only : dead_grid_lat, dead_grid_lon, dead_grid_index
  use dead_nuopc_mod   , only : dead_init_nuopc, dead_final_nuopc, dead_meshinit
  use dead_nuopc_mod   , only : fld_list_add, fld_list_realize, fldsMax, fld_list_type
  use dead_nuopc_mod   , only : ModelInitPhase, ModelSetRunClock

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

  integer                :: fldsToGlc_num = 0
  integer                :: fldsFrGlc_num = 0
  type (fld_list_type)   :: fldsToGlc(fldsMax)
  type (fld_list_type)   :: fldsFrGlc(fldsMax)
  integer, parameter     :: gridTofieldMap = 2 ! ungridded dimension is innermost

  real(r8), pointer      :: gbuf(:,:)            ! model info
  real(r8), pointer      :: lat(:)
  real(r8), pointer      :: lon(:)
  integer , allocatable  :: gindex(:)
  integer                :: nxg                  ! global dim i-direction
  integer                :: nyg                  ! global dim j-direction
  integer                :: my_task              ! my task in mpi communicator mpicom
  integer                :: inst_index           ! number of current instance (ie. 1)
  character(len=16)      :: inst_suffix = ""     ! char string associated with instance (ie. "_0001" or "")
  integer                :: logunit              ! logging unit number
  integer    ,parameter  :: master_task=0        ! task number of master task
  logical                :: mastertask
  integer                :: dbug = 1
  character(*),parameter :: modName =  "(xglc_comp_nuopc)"
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

    mastertask = my_task == master_task

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! set logunit and set shr logging to my log file
    !----------------------------------------------------------------------------

    call set_component_logging(gcomp, my_task==master_task, logunit, shrlogunit, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! Initialize xglc
    !----------------------------------------------------------------------------

    call dead_init_nuopc('glc', inst_suffix, logunit, lsize, gbuf, nxg, nyg)

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

    if (nxg /= 0 .and. nyg /= 0) then

       call fld_list_add(fldsFrGlc_num, fldsFrGlc, trim(flds_scalar_name))
       call fld_list_add(fldsFrGlc_num, fldsFrGlc, 'Sg_icemask'                )
       call fld_list_add(fldsFrGlc_num, fldsFrGlc, 'Sg_icemask_coupled_fluxes' )
       call fld_list_add(fldsFrGlc_num, fldsFrGlc, 'Sg_ice_covered'            )
       call fld_list_add(fldsFrGlc_num, fldsFrGlc, 'Sg_topo'                   )
       call fld_list_add(fldsFrGlc_num, fldsFrGlc, 'Flgg_hflx'                 )

       call fld_list_add(fldsToGlc_num, fldsToGlc, trim(flds_scalar_name))
       call fld_list_add(fldsToGlc_num, fldsToGlc, 'Sl_tsrf')
       call fld_list_add(fldsToGlc_num, fldsToGlc, 'Sl_topo')
       call fld_list_add(fldsToGlc_num, fldsToGlc, 'Flgg_hflx')

       do n = 1,fldsFrGlc_num
          if (mastertask) write(logunit,*)'Advertising From Xglc ',trim(fldsFrGlc(n)%stdname)
          call NUOPC_Advertise(exportState, standardName=fldsFrglc(n)%stdname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       enddo

       do n = 1,fldsToGlc_num
          if (mastertask) write(logunit,*)'Advertising To Xglc ',trim(fldsToGlc(n)%stdname)
          call NUOPC_Advertise(importState, standardName=fldsToglc(n)%stdname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       enddo
    end if

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=rc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR) :: convCIM, purpComp
    type(ESMF_Mesh)        :: Emesh
    integer                :: shrlogunit                ! original log unit
    integer                :: n
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=rc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logunit)

    !--------------------------------
    ! generate the mesh
    ! grid_option specifies grid or mesh
    !--------------------------------

    call dead_meshinit(gcomp, nxg, nyg, gindex, lon, lat, Emesh, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! realize the actively coupled fields, now that a mesh is established
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    !--------------------------------

    call fld_list_realize( &
         state=ExportState, &
         fldList=fldsFrGlc, &
         numflds=fldsFrGlc_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':dglcExport',&
         mesh=Emesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call fld_list_realize( &
         state=importState, &
         fldList=fldsToGlc, &
         numflds=fldsToGlc_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':dglcImport',&
         mesh=Emesh, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Pack export state
    ! Copy from d2x to exportState
    ! Set the coupling scalars
    !--------------------------------

    call state_setexport(exportState, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call state_setscalar(dble(nxg),flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call state_setscalar(dble(nyg),flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call state_diagnose(exportState,subname//':ES',rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"
    call ESMF_AttributeAdd(comp,  convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ShortName", "XGLC", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", "Land-Ice Dead Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
         "The dead models stand in as test model for active components." // &
         "Coupling data is artificially generated ", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2017", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Land-Ice", convention=convCIM, purpose=purpComp, rc=rc)
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
    integer           :: n
    integer           :: shrlogunit     ! original log unit
    real(r8), pointer :: dataptr(:)
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=rc)
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

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call state_diagnose(exportState,subname//':ES',rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (my_task == master_task) then
          call log_clock_advance(clock, 'XGLC', logunit, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       endif
    endif

    call shr_file_setLogUnit (shrlogunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=rc)

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

    ! Start from index 2 in order to skip the scalar field 
    do nf = 2,fldsFrGlc_num
       if (fldsFrGlc(nf)%ungridded_ubound == 0) then
          call field_setexport(exportState, trim(fldsFrGlc(nf)%stdname), lon, lat, nf=nf, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          do nind = 1,fldsFrGlc(nf)%ungridded_ubound
             call field_setexport(exportState, trim(fldsFrGlc(nf)%stdname), lon, lat, nf=nf+nind-1, &
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

    ncomp = 5
    if (present(ungridded_index)) then
       call ESMF_FieldGet(lfield, farrayPtr=data2d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (gridToFieldMap == 1) then
          do i = 1,size(data2d, dim=1)
             data2d(i,ungridded_index) = (nf*100) &
                  *  cos (pi*lat(i)/180.0_R8) *  cos (pi*lat(i)/180.0_R8) &
                  *  sin (pi*lon(i)/180.0_R8) *  sin (pi*lon(i)/180.0_R8) &
                  + (ncomp*10.0_R8)
          enddo
       else if (gridToFieldMap == 2) then
          do i = 1,size(data2d, dim=2)
             data2d(ungridded_index,i) = (nf*100) &
                  *  cos (pi*lat(i)/180.0_R8) *  cos (pi*lat(i)/180.0_R8) &
                  *  sin (pi*lon(i)/180.0_R8) *  sin (pi*lon(i)/180.0_R8) &
                  + (ncomp*10.0_R8)
          end do
       end if
    else
       call ESMF_FieldGet(lfield, farrayPtr=data1d, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (fldname == 'Sg_icemask' .or.  fldname == 'Sg_icemask_coupled_fluxes' .or. fldname == 'Sg_ice_covered') then
          data1d(:) = 1._r8
       else
          do i = 1,size(data1d)
             data1d(i) = (nf*100) &
                  *  cos (pi*lat(i)/180.0_R8) *  cos (pi*lat(i)/180.0_R8) &
                  *  sin (pi*lon(i)/180.0_R8) *  sin (pi*lon(i)/180.0_R8) &
                  + (ncomp*10.0_R8)
          end do
       end if
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

    call dead_final_nuopc('glc', logunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=rc)

  end subroutine ModelFinalize

  !===============================================================================

end module glc_comp_nuopc
