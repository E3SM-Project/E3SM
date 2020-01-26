module dshr_methods_mod

  use ESMF         , only : operator(<), operator(/=), operator(+)
  use ESMF         , only : operator(-), operator(*) , operator(>=)
  use ESMF         , only : operator(<=), operator(>), operator(==)
  use ESMF         , only : ESMF_LOGERR_PASSTHRU, ESMF_LogFoundError, ESMF_LOGMSG_ERROR, ESMF_MAXSTR
  use ESMF         , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_FAILURE
  use ESMF         , only : ESMF_State, ESMF_StateGet
  use ESMF         , only : ESMF_Field, ESMF_FieldGet
  use ESMF         , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSet
  use ESMF         , only : ESMF_GeomType_Flag, ESMF_FieldStatus_Flag
  use ESMF         , only : ESMF_Mesh, ESMF_MeshGet
  use ESMF         , only : ESMF_GEOMTYPE_MESH, ESMF_GEOMTYPE_GRID, ESMF_FIELDSTATUS_COMPLETE
  use ESMF         , only : ESMF_Clock, ESMF_ClockCreate, ESMF_ClockGet, ESMF_ClockSet 
  use ESMF         , only : ESMF_ClockPrint, ESMF_ClockAdvance 
  use ESMF         , only : ESMF_Alarm, ESMF_AlarmCreate, ESMF_AlarmGet, ESMF_AlarmSet
  use ESMF         , only : ESMF_Calendar, ESMF_CALKIND_NOLEAP, ESMF_CALKIND_GREGORIAN
  use ESMF         , only : ESMF_Time, ESMF_TimeGet, ESMF_TimeSet
  use ESMF         , only : ESMF_TimeInterval, ESMF_TimeIntervalSet, ESMF_TimeIntervalGet
  use ESMF         , only : ESMF_VM, ESMF_VMGet, ESMF_VMBroadcast, ESMF_VMGetCurrent
  use NUOPC        , only : NUOPC_CompAttributeGet
  use NUOPC_Model  , only : NUOPC_ModelGet
  use shr_kind_mod , only : r8 => shr_kind_r8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod  , only : shr_sys_abort
  use shr_file_mod , only : shr_file_setlogunit, shr_file_getLogUnit
  use perf_mod     , only : t_startf, t_stopf 

  implicit none
  private

  public  :: memcheck
  public  :: get_component_instance
  public  :: set_component_logging
  public  :: log_clock_advance
  public  :: state_getscalar
  public  :: state_setscalar
  public  :: state_diagnose
  public  :: state_getfldptr
  public  :: FB_Regrid
  public  :: FB_getFieldN
  public  :: FB_getNameN
  public  :: FB_getFldPtr
  public  :: FB_diagnose
  public  :: FB_fldchk
  public  :: chkerr

  private :: field_getfldptr


  ! used/reused in module
  logical                     :: isPresent
  character(len=1024)         :: msgString
  type(ESMF_FieldStatus_Flag) :: status

  ! Module data
  integer, parameter          :: memdebug_level=1
  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine memcheck(string, level, mastertask)

    ! input/output variables
    character(len=*) , intent(in) :: string
    integer          , intent(in) :: level
    logical          , intent(in) :: mastertask

    ! local vvariables
    integer :: ierr
    integer, external :: GPTLprint_memusage
    !-----------------------------------------------------------------------

    if ((mastertask .and. memdebug_level > level) .or. memdebug_level > level+1) then
       ierr = GPTLprint_memusage(string)
    endif

  end subroutine memcheck

!===============================================================================

  subroutine get_component_instance(gcomp, inst_suffix, inst_index, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    character(len=*) , intent(out) :: inst_suffix
    integer          , intent(out) :: inst_index
    integer          , intent(out) :: rc

    ! local variables
    logical          :: isPresent
    character(len=4) :: cvalue
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", value=inst_suffix, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       cvalue = inst_suffix(2:)
       read(cvalue, *) inst_index
    else
       inst_suffix = ""
       inst_index=1
    endif

  end subroutine get_component_instance

!===============================================================================

  subroutine set_component_logging(gcomp, mastertask, logunit, shrlogunit, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    logical, intent(in)  :: mastertask
    integer, intent(out) :: logunit
    integer, intent(out) :: shrlogunit
    integer, intent(out) :: rc

    ! local variables
    character(len=CL) :: diro
    character(len=CL) :: logfile
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    shrlogunit = 6

    if (mastertask) then
       call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       open(newunit=logunit,file=trim(diro)//"/"//trim(logfile))
    else
       logUnit = 6
    endif

    call shr_file_setLogUnit (logunit)

  end subroutine set_component_logging

!===============================================================================

  subroutine log_clock_advance(clock, component, logunit, rc)

    ! input/output variables
    type(ESMF_Clock)               :: clock
    character(len=*) , intent(in)  :: component
    integer          , intent(in)  :: logunit
    integer          , intent(out) :: rc

    ! local variables
    character(len=CL) :: cvalue, prestring
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    write(prestring, *) "------>Advancing ",trim(component)," from: "
    call ESMF_ClockPrint(clock, options="currTime", unit=cvalue, preString=trim(prestring), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

    call ESMF_ClockPrint(clock, options="stopTime", unit=cvalue, &
         preString="--------------------------------> to: ", rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write(logunit, *) trim(cvalue)

  end subroutine log_clock_advance

!===============================================================================

  subroutine state_getscalar(state, scalar_id, scalar_value, flds_scalar_name, flds_scalar_num, rc)

    ! ----------------------------------------------
    ! Get scalar data from State for a particular name and broadcast it to all other pets
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State), intent(in)     :: state
    integer,          intent(in)     :: scalar_id
    real(r8),         intent(out)    :: scalar_value
    character(len=*), intent(in)     :: flds_scalar_name
    integer,          intent(in)     :: flds_scalar_num
    integer,          intent(inout)  :: rc

    ! local variables
    integer           :: mytask, ierr, len
    type(ESMF_VM)     :: vm
    type(ESMF_Field)  :: field
    real(r8), pointer :: farrayptr(:,:)
    real(r8)          :: tmp(1)
    character(len=*), parameter :: subname='(state_getscalar)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=field, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (mytask == 0) then
      call ESMF_FieldGet(field, farrayPtr = farrayptr, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
        call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO, line=__LINE__, file=u_FILE_u)
        rc = ESMF_FAILURE
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
      endif
      tmp(:) = farrayptr(scalar_id,:)
    endif
    call ESMF_VMBroadCast(vm, tmp, 1, 0, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    scalar_value = tmp(1)

  end subroutine state_getscalar

!================================================================================

  subroutine state_setscalar(scalar_value, scalar_id, State, flds_scalar_name, flds_scalar_num,  rc)

    ! ----------------------------------------------
    ! Set scalar data from State for a particular name
    ! ----------------------------------------------

    ! input/output arguments
    real(r8),         intent(in)     :: scalar_value
    integer,          intent(in)     :: scalar_id
    type(ESMF_State), intent(inout)  :: State
    character(len=*), intent(in)     :: flds_scalar_name
    integer,          intent(in)     :: flds_scalar_num
    integer,          intent(inout)  :: rc

    ! local variables
    integer           :: mytask
    type(ESMF_Field)  :: lfield
    type(ESMF_VM)     :: vm
    real(r8), pointer :: farrayptr(:,:)
    character(len=*), parameter :: subname='(state_setscalar)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=mytask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_StateGet(State, itemName=trim(flds_scalar_name), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (mytask == 0) then
       call ESMF_FieldGet(lfield, farrayPtr = farrayptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       if (scalar_id < 0 .or. scalar_id > flds_scalar_num) then
          call ESMF_LogWrite(trim(subname)//": ERROR in scalar_id", ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       endif
       farrayptr(scalar_id,1) = scalar_value
    endif

  end subroutine state_setscalar

!===============================================================================

  subroutine state_diagnose(State, string, rc)

    ! ----------------------------------------------
    ! Diagnose status of State
    ! ----------------------------------------------

    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: string
    integer         , intent(out) :: rc

    ! local variables
    integer                         :: i,j,n
    type(ESMf_Field)                :: lfield
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    real(r8), pointer               :: dataPtr1d(:)
    real(r8), pointer               :: dataPtr2d(:,:)
    character(len=*),parameter      :: subname='(state_diagnose)'
    ! ----------------------------------------------

    call ESMF_StateGet(state, itemCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    call ESMF_StateGet(state, itemNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1, fieldCount

       call ESMF_StateGet(state, itemName=lfieldnamelist(n), field=lfield, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call field_getfldptr(lfield, fldptr1=dataPtr1d, fldptr2=dataPtr2d, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data
       elseif (lrank == 1) then
          if (size(dataPtr1d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(string)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
          else
             write(msgString,'(A,a)') trim(string)//': '//trim(lfieldnamelist(n))," no data"
          endif
       elseif (lrank == 2) then
          if (size(dataPtr2d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(string)//': '//trim(lfieldnamelist(n)), &
                  minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
          else
             write(msgString,'(A,a)') trim(string)//': '//trim(lfieldnamelist(n))," no data"
          endif
       else
          call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       endif
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    enddo

    deallocate(lfieldnamelist)

  end subroutine state_diagnose

!===============================================================================

  subroutine State_GetFldPtr(State, fldname, fldptr1, fldptr2, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    use ESMF, only : ESMF_State, ESMF_Field, ESMF_StateGet

    ! input/output variables
    type(ESMF_State) ,          intent(in)              :: State
    character(len=*) ,          intent(in)              :: fldname
    real(R8)         , pointer, intent(inout), optional :: fldptr1(:)
    real(R8)         , pointer, intent(inout), optional :: fldptr2(:,:)
    integer          ,          intent(out)             :: rc

    ! local variables
    type(ESMF_Field)           :: lfield
    character(len=*), parameter :: subname='(med_methods_State_GetFldPtr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call Field_GetFldPtr(lfield, fldptr1=fldptr1, fldptr2=fldptr2, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine State_GetFldPtr

!===============================================================================

  subroutine field_getfldptr(field, fldptr1, fldptr2, rank, abort, rc)

    ! ----------------------------------------------
    ! for a field, determine rank and return fldptr1 or fldptr2
    ! abort is true by default and will abort if fldptr is not yet allocated in field
    ! rank returns 0, 1, or 2.  0 means fldptr not allocated and abort=false
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_Field)  , intent(in)              :: field
    real(r8), pointer , intent(inout), optional :: fldptr1(:)
    real(r8), pointer , intent(inout), optional :: fldptr2(:,:)
    integer           , intent(out)  , optional :: rank
    logical           , intent(in)   , optional :: abort
    integer           , intent(out)             :: rc

    ! local variables
    type(ESMF_GeomType_Flag)    :: geomtype
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Mesh)             :: lmesh
    integer                     :: lrank, nnodes, nelements
    logical                     :: labort
    character(len=*), parameter :: subname='(field_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    labort = .true.
    if (present(abort)) then
       labort = abort
    endif
    lrank = -99

    call ESMF_FieldGet(field, status=status, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then

       lrank = 0
       if (labort) then
          call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
          rc = ESMF_FAILURE
          return
       else
          call ESMF_LogWrite(trim(subname)//": WARNING data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       endif

    else

       call ESMF_FieldGet(field, geomtype=geomtype, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (geomtype == ESMF_GEOMTYPE_GRID) then
          call ESMF_FieldGet(field, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       elseif (geomtype == ESMF_GEOMTYPE_MESH) then
          call ESMF_FieldGet(field, rank=lrank, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_FieldGet(field, mesh=lmesh, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          if (nnodes == 0 .and. nelements == 0) lrank = 0
       else  
          call ESMF_LogWrite(trim(subname)//": ERROR geomtype not supported ", &
               ESMF_LOGMSG_INFO, rc=rc)
          rc = ESMF_FAILURE
          return
       endif ! geomtype

       if (lrank == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", &
               ESMF_LOGMSG_INFO)
       elseif (lrank == 1) then
          if (.not.present(fldptr1)) then
             call ESMF_LogWrite(trim(subname)//": ERROR missing rank=1 array ", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
             rc = ESMF_FAILURE
             return
          endif
          call ESMF_FieldGet(field, farrayPtr=fldptr1, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       elseif (lrank == 2) then
          if (.not.present(fldptr2)) then
             call ESMF_LogWrite(trim(subname)//": ERROR missing rank=2 array ", &
                  ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
             rc = ESMF_FAILURE
             return
          endif
          call ESMF_FieldGet(field, farrayPtr=fldptr2, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       else
          call ESMF_LogWrite(trim(subname)//": ERROR in rank ", &
               ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u)
          rc = ESMF_FAILURE
          return
       endif

    endif  ! status

    if (present(rank)) then
       rank = lrank
    endif

  end subroutine field_getfldptr

!===============================================================================

  subroutine FB_diagnose(FB, string, rc)
    ! ----------------------------------------------
    ! Diagnose status of FB
    ! ----------------------------------------------

    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet

    type(ESMF_FieldBundle) , intent(inout)        :: FB
    character(len=*)       , intent(in), optional :: string
    integer                , intent(out)          :: rc

    ! local variables
    integer                         :: i,j,n
    integer                         :: fieldCount, lrank
    character(ESMF_MAXSTR), pointer :: lfieldnamelist(:)
    character(len=CL)               :: lstring
    real(R8), pointer               :: dataPtr1d(:)
    real(R8), pointer               :: dataPtr2d(:,:)
    character(len=*), parameter     :: subname='(FB_diagnose)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    lstring = ''
    if (present(string)) lstring = trim(string) // ' '

    ! Determine number of fields in field bundle and allocate memory for lfieldnamelist
    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldnamelist(fieldCount))

    ! Get the fields in the field bundle
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! For each field in the bundle, get its memory location and print out the field
    do n = 1, fieldCount
       call FB_GetFldPtr(FB, lfieldnamelist(n), fldptr1=dataPtr1d, fldptr2=dataPtr2d, rank=lrank, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (lrank == 0) then
          ! no local data

       elseif (lrank == 1) then
          if (size(dataPtr1d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n))//' ', &
                  minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), " no data"
          endif

       elseif (lrank == 2) then
          if (size(dataPtr2d) > 0) then
             write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n))//' ', &
                  minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
          else
             write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(lfieldnamelist(n)), &
                  " no data"
          endif

       else
          call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
          rc = ESMF_FAILURE
          return
       endif
       call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)
    enddo

    ! Deallocate memory
    deallocate(lfieldnamelist)

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine FB_diagnose

!===============================================================================

  subroutine FB_GetFldPtr(FB, fldname, fldptr1, fldptr2, rank, field, rc)

    ! ----------------------------------------------
    ! Get pointer to a field bundle field
    ! ----------------------------------------------

    use ESMF , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_Field

    type(ESMF_FieldBundle) , intent(in)              :: FB
    character(len=*)       , intent(in)              :: fldname
    real(R8), pointer      , intent(inout), optional :: fldptr1(:)
    real(R8), pointer      , intent(inout), optional :: fldptr2(:,:)
    integer                , intent(out),   optional :: rank
    type(ESMF_Field)       , intent(out),   optional :: field
    integer                , intent(out)              :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer          :: lrank
    character(len=*), parameter :: subname='(FB_GetFldPtr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (.not. FB_FldChk(FB, trim(fldname), rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": ERROR field "//trim(fldname)//" not in FB ", ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
    endif

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), field=lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call Field_GetFldPtr(lfield, fldptr1=fldptr1, fldptr2=fldptr2, rank=lrank, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (present(rank)) rank = lrank
    if (present(field)) field = lfield

  end subroutine FB_GetFldPtr

!===============================================================================

  subroutine FB_Regrid(FBin, FBout, RH, zeroregion, rc)

    ! ----------------------------------------------
    ! Assumes that FBin and FBout contain fields with the same name
    ! ----------------------------------------------

    use ESMF , only : ESMF_FieldBundle, ESMF_FieldBundleGet
    use ESMF , only : ESMF_RouteHandle, ESMF_FieldRegrid, ESMF_Field
    use ESMF , only : ESMF_TERMORDER_SRCSEQ, ESMF_FieldRegridStore, ESMF_SparseMatrixWrite
    use ESMF , only : ESMF_Region_Flag, ESMF_REGION_TOTAL, ESMF_MAXSTR

    ! input/output variables
    type(ESMF_FieldBundle), intent(inout)        :: FBin
    type(ESMF_FieldBundle), intent(inout)        :: FBout
    type(ESMF_RouteHandle), intent(inout)        :: RH
    type(ESMF_Region_Flag), intent(in), optional :: zeroregion
    integer               , intent(out)          :: rc

    ! local
    type(ESMF_Field)           :: field_src
    type(ESMF_Field)           :: field_dst
    integer                    :: fieldcount
    logical                    :: checkflag = .false.
    character(len=8)           :: filename
    type(ESMF_Region_Flag)     :: localzr
    character(CS)              :: fldname
    integer                    :: n
    character(ESMF_MAXSTR), allocatable :: lfieldNameList(:)
    character(len=*),parameter :: subname='(FB_FieldRegrid)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf(subname)

    localzr = ESMF_REGION_TOTAL
    if (present(zeroregion)) then
       localzr = zeroregion
    endif

    call ESMF_FieldBundleGet(FBin, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(lfieldNameList(fieldCount))
    call ESMF_FieldBundleGet(FBin, fieldNameList=lfieldNameList, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    do n = 1,fieldCount
       fldname = trim(lfieldnamelist(n))
       call ESMF_FieldBundleGet(FBin, fieldName=trim(fldname), field=field_src, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleGet(FBout, fieldName=trim(fldname), field=field_dst, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call ESMF_FieldRegrid(field_src, field_dst, routehandle=RH, &
            termorderflag=ESMF_TERMORDER_SRCSEQ, checkflag=.false., zeroregion=localzr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end do

    call t_stopf(subname)

  end subroutine FB_Regrid

!===============================================================================

  subroutine FB_getFieldN(FB, fieldnum, field, rc)

    ! ----------------------------------------------
    ! Get field with number fieldnum in input field bundle FB
    ! ----------------------------------------------

    use ESMF, only : ESMF_Field, ESMF_FieldBundle, ESMF_FieldBundleGet

    ! input/output variables
    type(ESMF_FieldBundle), intent(in)    :: FB
    integer               , intent(in)    :: fieldnum
    type(ESMF_Field)      , intent(inout) :: field
    integer               , intent(out)   :: rc

    ! local variables
    character(len=ESMF_MAXSTR) :: name
    character(len=*),parameter :: subname='(FB_getFieldN)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call FB_getNameN(FB, fieldnum, name, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldBundleGet(FB, fieldName=name, field=field, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

  end subroutine FB_getFieldN

!===============================================================================

  subroutine FB_getNameN(FB, fieldnum, fieldname, rc)

    ! ----------------------------------------------
    ! Get name of field number fieldnum in input field bundle FB
    ! ----------------------------------------------

    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet

    ! input/output variables
    type(ESMF_FieldBundle), intent(in)    :: FB
    integer               , intent(in)    :: fieldnum
    character(len=*)      , intent(out)   :: fieldname
    integer               , intent(out)   :: rc

    ! local variables
    integer                         :: fieldCount
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    character(len=*),parameter      :: subname='(FB_getNameN)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    fieldname = ' '
    call ESMF_FieldBundleGet(FB, fieldCount=fieldCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (fieldnum > fieldCount) then
      call ESMF_LogWrite(trim(subname)//": ERROR fieldnum > fieldCount ", ESMF_LOGMSG_ERROR)
      rc = ESMF_FAILURE
      return
    endif

    allocate(lfieldnamelist(fieldCount))
    call ESMF_FieldBundleGet(FB, fieldNameList=lfieldnamelist, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    fieldname = lfieldnamelist(fieldnum)
    deallocate(lfieldnamelist)

  end subroutine FB_getNameN

!===============================================================================

  logical function FB_FldChk(FB, fldname, rc)

    ! ----------------------------------------------
    ! Determine if field with fldname is in input field bundle
    ! ----------------------------------------------

    use ESMF, only : ESMF_FieldBundle, ESMF_FieldBundleGet
    use ESMF, only : ESMF_FieldBundleIsCreated

    ! input/output variables
    type(ESMF_FieldBundle), intent(in)  :: FB
    character(len=*)      , intent(in)  :: fldname
    integer               , intent(out) :: rc

    ! local variables
    character(len=*), parameter :: subname='(FB_FldChk)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! If field bundle is created determine if fldname is present in field bundle
    FB_FldChk = .false.

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), isPresent=isPresent, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) then
       call ESMF_LogWrite(trim(subname)//" Error checking field: "//trim(fldname), ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    endif

    if (isPresent) then
       FB_FldChk = .true.
    endif

  end function FB_FldChk

!===============================================================================

  subroutine FB_Field_diagnose(FB, fieldname, string, rc)

    ! ----------------------------------------------
    ! Diagnose status of State
    ! ----------------------------------------------

    use ESMF, only : ESMF_FieldBundle

    ! input/output variables
    type(ESMF_FieldBundle), intent(inout)  :: FB
    character(len=*), intent(in)           :: fieldname
    character(len=*), intent(in), optional :: string
    integer         , intent(out)          :: rc

    ! local variables
    integer           :: lrank
    character(len=CS) :: lstring
    real(R8), pointer :: dataPtr1d(:)
    real(R8), pointer :: dataPtr2d(:,:)
    character(len=*),parameter      :: subname='(FB_FieldDiagnose)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    lstring = ''
    if (present(string)) lstring = trim(string)

    call FB_GetFldPtr(FB, fieldname, dataPtr1d, dataPtr2d, lrank, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (lrank == 0) then
       ! no local data
    elseif (lrank == 1) then
       if (size(dataPtr1d) > 0) then
          write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               minval(dataPtr1d), maxval(dataPtr1d), sum(dataPtr1d), size(dataPtr1d)
       else
          write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname)," no data"
       endif
    elseif (lrank == 2) then
       if (size(dataPtr2d) > 0) then
          write(msgString,'(A,3g14.7,i8)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname), &
               minval(dataPtr2d), maxval(dataPtr2d), sum(dataPtr2d), size(dataPtr2d)
       else
          write(msgString,'(A,a)') trim(subname)//' '//trim(lstring)//': '//trim(fieldname)," no data"
       endif
    else
       call ESMF_LogWrite(trim(subname)//": ERROR rank not supported ", ESMF_LOGMSG_ERROR)
       rc = ESMF_FAILURE
       return
    endif
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

  end subroutine FB_Field_diagnose

!===============================================================================

  logical function chkerr(rc, line, file)

    integer, intent(in) :: rc
    integer, intent(in) :: line
    character(len=*), intent(in) :: file

    integer :: lrc

    chkerr = .false.
    lrc = rc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
       chkerr = .true.
    endif
  end function chkerr

end module dshr_methods_mod
