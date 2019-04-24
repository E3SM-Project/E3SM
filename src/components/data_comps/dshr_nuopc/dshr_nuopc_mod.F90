module dshr_nuopc_mod

  use NUOPC
  use NUOPC_Model
  use ESMF
  use dshr_methods_mod , only : alarmInit, chkerr
  use shr_kind_mod     , only : r8=>shr_kind_r8, cs=>shr_kind_cs, cxx=>shr_kind_cxx
  use shr_string_mod   , only : shr_string_listGetIndex
  use shr_sys_mod      , only : shr_sys_abort

  implicit none
  public

  public :: dshr_fld_add
  public :: dshr_import
  public :: dshr_export
  public :: dshr_realize
  public :: ModelInitPhase   ! TODO: rename to dshr_modelinit
  public :: ModelSetRunClock ! TODO: rename to dshr_setrunclock
  public :: ModelSetMetaData ! TODO rename to dshr_setmetadata

  type fld_list_type
    character(len=128) :: stdname
    integer :: ungridded_lbound = 0
    integer :: ungridded_ubound = 0
  end type fld_list_type

  interface dshr_fld_add ; module procedure &
       dshr_fld_add, &
       dshr_fld_add_model, &
       dshr_fld_add_model_and_data
  end interface dshr_fld_add

  integer                 :: gridTofieldMap = 2 ! ungridded dimension is innermost
  integer     , parameter :: fldsMax = 100
  integer     , parameter :: dbug = 10
  character(*), parameter :: modName =  "(dhsr_nuopc_mod)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dshr_fld_add(med_fld, fldlist_num, fldlist, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    character(len=*)               , intent(in)    :: med_fld
    integer                        , intent(inout) :: fldlist_num
    type(fld_list_type)            , intent(inout) :: fldlist(:)
    integer             , optional , intent(in)    :: ungridded_lbound
    integer             , optional , intent(in)    :: ungridded_ubound

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(dshr_nuopc_mod:dshr_fld_add)'
    ! ----------------------------------------------

    call dshr_fld_list_add(fldlist_num, fldlist, med_fld, ungridded_lbound, ungridded_ubound)

  end subroutine dshr_fld_add

!===============================================================================

  subroutine dshr_fld_add_model(model_fld, model_fld_concat, model_fld_index, &
       fldlist_num, fldlist, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    character(len=*)               , intent(in)    :: model_fld
    character(len=*)               , intent(inout) :: model_fld_concat
    integer             , optional , intent(out)   :: model_fld_index
    integer             , optional , intent(inout) :: fldlist_num
    type(fld_list_type) , optional , intent(inout) :: fldlist(:)
    integer             , optional , intent(in)    :: ungridded_lbound
    integer             , optional , intent(in)    :: ungridded_ubound

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(dshr_nuopc_mod:dshr_fld_add_model)'
    ! ----------------------------------------------

    if (len_trim(model_fld_concat) + len_trim(model_fld) + 1 >= len(model_fld_concat)) then
       call ESMF_LogWrite(subname//': ERROR: max len of model_fld_concat has been exceeded', ESMF_LOGMSG_INFO)
       rc = ESMF_FAILURE
       return
    end if

    if (trim(model_fld_concat) == '') then
       model_fld_concat = trim(model_fld)
    else
       model_fld_concat = trim(model_fld_concat)//':'//trim(model_fld)
    end if

    if (present(model_fld_index)) then
       call shr_string_listGetIndex(trim(model_fld_concat), trim(model_fld),  model_fld_index)
    end if

    !----------------------------------
    ! Update fldlist array if appropriate
    !----------------------------------

    if (present(fldlist_num) .and. present(fldlist)) then
       call dshr_fld_list_add(fldlist_num, fldlist, model_fld, ungridded_lbound, ungridded_ubound)
    end if

  end subroutine dshr_fld_add_model

  !===============================================================================

  subroutine dshr_fld_add_model_and_data( data_fld, data_fld_array, &
       model_fld, model_fld_array, model_fld_concat, model_fld_index, &
       fldlist_num, fldlist, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    character(len=*)               , intent(in)    :: data_fld
    character(len=*)               , pointer       :: data_fld_array(:)
    character(len=*)               , intent(in)    :: model_fld
    character(len=*)               , pointer       :: model_fld_array(:)
    character(len=*)    , optional , intent(inout) :: model_fld_concat
    integer             , optional , intent(out)   :: model_fld_index
    integer             , optional , intent(inout) :: fldlist_num
    type(fld_list_type) , optional , intent(inout) :: fldlist(:)
    integer             , optional , intent(in)    :: ungridded_lbound
    integer             , optional , intent(in)    :: ungridded_ubound

    ! local variables
    integer                     :: rc
    integer                     :: n, oldsize, id
    character(len=CS), pointer  :: new_data_fld_array(:)
    character(len=CS), pointer  :: new_model_fld_array(:)
    character(len=*), parameter :: subname='(dshr_nuopc_mod:dshr_fld_add_model_and_data) '
    ! ----------------------------------------------

    !----------------------------------
    ! Create new data_fld_array and model_fld_array
    ! Model is what the data model sends and receives from the mediator
    ! Data is what the data model obtains from the various streams
    !----------------------------------

    ! 1) determine new index
    if (associated(data_fld_array)) then
       oldsize = size(data_fld_array)
    else
       oldsize = 0
    end if
    id = oldsize + 1

    ! 2) allocate new_data_fld_array and oldavi to one element larger than input
    allocate(new_data_fld_array(id))
    allocate(new_model_fld_array(id))

    ! 3) copy data_fld_array and model_fld_array into first N-1 elements of data_fld_arrays and model_fld_array
    do n = 1,oldsize
       new_data_fld_array(n)  = data_fld_array(n)
       new_model_fld_array(n) = model_fld_array(n)
    end do

    ! 4) deallocate / nullify data_fld_array and model_fld_array
    if (oldsize >  0) then
       deallocate(data_fld_array)
       deallocate(model_fld_array)
       nullify(data_fld_array)
       nullify(model_fld_array)
    end if

    ! 5) point data_fld_array => new_data_fld_array and
    ! model_fld_array => new_model_fld_array and update info for new entry
    data_fld_array  => new_data_fld_array
    model_fld_array => new_model_fld_array
    data_fld_array(id)  = trim(data_fld)
    model_fld_array(id) = trim(model_fld)

    !----------------------------------
    ! Update flds_concat colon delimited string if appropriate
    !----------------------------------

    if (present(model_fld_concat)) then
       if (len_trim(model_fld_concat) + len_trim(model_fld) + 1 >= cxx) then
          call ESMF_LogWrite(subname//': ERROR: max len of model_fld_concat has been exceeded', ESMF_LOGMSG_INFO)
          call shr_sys_abort()
       end if
       if (trim(model_fld_concat) == '') then
          model_fld_concat = trim(model_fld)
       else
          model_fld_concat = trim(model_fld_concat)//':'//trim(model_fld)
       end if

       ! Get model field index if appropriated
       if (present(model_fld_index)) then
          call shr_string_listGetIndex(trim(model_fld_concat), trim(model_fld),  model_fld_index)
       end if
    end if

    !----------------------------------
    ! Update fldlist array if appropriate
    !----------------------------------
    if (present(fldlist_num) .and. present(fldlist)) then
       call dshr_fld_list_add(fldlist_num, fldlist, model_fld, ungridded_lbound, ungridded_ubound)
    end if

  end subroutine dshr_fld_add_model_and_data

  !===============================================================================

  subroutine dshr_fld_list_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname
    integer,          optional, intent(in)    :: ungridded_lbound
    integer,          optional, intent(in)    :: ungridded_ubound

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(dshr_nuopc_mod:fld_list_add)'
    !----------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), ESMF_LOGMSG_INFO)
      rc = ESMF_FAILURE
      return
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine dshr_fld_list_add

  !===============================================================================

  subroutine dshr_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    ! input/output variables
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
    character(len=*),parameter  :: subname='(dshr_nuopc_mod:fld_list_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          else
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                     ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                     ungriddedUbound=(/fldlist(n)%ungridded_ubound/), gridToFieldMap=(/gridToFieldMap/), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
             end if
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO)
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
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

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(dshr_nuopc_mod:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine dshr_realize

  !===============================================================================

  subroutine ModelInitPhase(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine ModelInitPhase

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    ! input/output variables
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
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname='dshr_nuopc_mod:(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart alarm
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelSetMetadata(gcomp, name, rc)

    type(ESMF_GridComp)            :: gcomp
    character(len=*) , intent(in)  :: name
    integer          , intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR)  :: convCIM, purpComp

    rc = ESMF_SUCCESS

    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"
    call ESMF_AttributeAdd(gcomp, convention=convCIM, purpose=purpComp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_AttributeSet(gcomp, "ShortName", trim(name), convention=convCIM, purpose=purpComp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_AttributeSet(gcomp, "LongName", "Climatological SeaIce Data Model", convention=convCIM, purpose=purpComp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_AttributeSet(gcomp, "Description", &
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
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_AttributeSet(gcomp, "ReleaseDate", "2010", convention=convCIM, purpose=purpComp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_AttributeSet(gcomp, "ModelType", "SeaIce", convention=convCIM, purpose=purpComp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_AttributeSet(gcomp, "Name", "TBD", convention=convCIM, purpose=purpComp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_AttributeSet(gcomp, "EmailAddress", "TBD", convention=convCIM, purpose=purpComp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_AttributeSet(gcomp, "ResponsiblePartyRole", "contact", convention=convCIM, purpose=purpComp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine ModelSetMetadata

  !-----------------------------------------------------------------------------

  subroutine dshr_export(array, state, fldname, ungridded_index, rc)

    ! ----------------------------------
    ! copy array data to state fields
    ! ----------------------------------

    ! input/otuput variables
    real(r8)         , intent(inout) :: array(:)
    type(ESMF_State) , intent(inout) :: state
    character(len=*) , intent(in)    :: fldname
    integer, optional, intent(in)    :: ungridded_index
    integer          , intent(out)   :: rc

    ! local variables
    integer           :: lsize, n
    type(ESMF_Field)  :: lfield
    real(R8), pointer :: farray1d(:)
    real(R8), pointer :: farray2d(:,:)
    character(*),parameter :: subName = "(dshr_nuopc_mod: dshr_export)"
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemName=trim(fldname), field=lfield, rc=rc)
    if (.not. ChkErr(rc,__LINE__,u_FILE_u)) then
       call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" copy", ESMF_LOGMSG_INFO)

       lsize = size(array)
       if (present(ungridded_index)) then
          call ESMF_FieldGet(lfield, farrayPtr=farray2d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (gridToFieldMap == 1) then
             do n = 1,lsize
                farray2d(n,ungridded_index) = array(n)
             enddo
          else if (gridToFieldMap == 2) then
             do n = 1,lsize
                farray2d(ungridded_index,n) = array(n)
             end do
          end if
       else
          call ESMF_FieldGet(lfield, farrayPtr=farray1d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n = 1,lsize
             farray1d(n) = array(n)
          enddo
       end if
    end if

  end subroutine dshr_export

  !-----------------------------------------------------------------------------

  subroutine dshr_import(state, fldname, array, ungridded_index, rc)

    ! ----------------------------------
    ! copy state field to array data
    ! ----------------------------------

    ! input/output variables
    type(ESMF_State)  , intent(in)    :: state
    character(len=*)  , intent(in)    :: fldname
    real(r8)          , intent(inout) :: array(:)
    integer, optional , intent(in)    :: ungridded_index
    integer           , intent(out)   :: rc

    ! local variables
    integer           :: lsize, n
    type(ESMF_Field)  :: lfield
    real(R8), pointer :: farray1d(:)
    real(R8), pointer :: farray2d(:,:)
    character(*),parameter :: subName = "(dshr_nuopc_mod: dshr_import)"
    !----------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemName=trim(fldname), field=lfield, rc=rc)
    if (.not. ChkErr(rc,__LINE__,u_FILE_u)) then
       call ESMF_LogWrite(trim(subname)//": fldname = "//trim(fldname)//" copy", ESMF_LOGMSG_INFO)

       lsize = size(array)
       if (present(ungridded_index)) then
          call ESMF_FieldGet(lfield, farrayPtr=farray2d, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (gridToFieldMap == 1) then
             do n = 1,lsize
                array(n) = farray2d(n,ungridded_index)
             enddo
          else if (gridToFieldMap == 2) then
             do n = 1,lsize
                array(n) = farray2d(ungridded_index,n)
             enddo
          end if
       else
          call ESMF_FieldGet(lfield, farrayPtr=farray1d,  rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n = 1,lsize
             array(n) = farray1d(n)
          enddo
       end if
    end if

   end subroutine dshr_import

end module dshr_nuopc_mod
