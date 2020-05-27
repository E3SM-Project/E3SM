module dshr_fldlist_mod

  use NUOPC
  use ESMF
  use shr_kind_mod     , only : r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx
  use dshr_methods_mod , only : chkerr

  implicit none
  private

  public :: dshr_fldlist_add
  public :: dshr_fldlist_realize

  type, public :: fldlist_type
    character(len=CS) :: stdname
    integer :: ungridded_lbound = 0
    integer :: ungridded_ubound = 0
    type(fldlist_type), pointer :: next => null()
  end type fldlist_type

  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dshr_fldlist_add(fldlists, fldname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    character(len=*)  , intent(in)    :: fldname
    type(fldlist_type), pointer       :: fldlists
    integer, optional , intent(in)    :: ungridded_lbound
    integer, optional , intent(in)    :: ungridded_ubound

    ! local variables
    type(fldlist_type), pointer :: fldlist_new
    integer :: rc
    ! ----------------------------------------------

    allocate(fldlist_new)
    fldlist_new%next => fldlists
    fldlists => fldlist_new

    fldlist_new%stdname = trim(fldname)
    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist_new%ungridded_lbound = ungridded_lbound
       fldlist_new%ungridded_ubound = ungridded_ubound
    end if
  end subroutine dshr_fldlist_add
  
  !===============================================================================

  subroutine dshr_fldlist_realize(state, fldLists, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(fldlist_type)  , pointer       :: fldLists
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    type(ESMF_Mesh)     , intent(in)    :: mesh
    character(len=*)    , intent(in)    :: tag
    integer             , intent(inout) :: rc

    ! local variables
    integer                     :: n
    type(fldlist_type), pointer :: fldList
    type(ESMF_Field)            :: field
    character(len=CS)           :: stdname
    character(len=*),parameter  :: subname='(dshr_fldList_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    fldList => fldLists ! note that fldlists is the head of the linked list
    do while (associated(fldList))
       stdname = fldList%stdname

       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          else
             ! Create the field
             if (fldList%ungridded_lbound > 0 .and. fldList%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                     ungriddedLbound=(/fldList%ungridded_lbound/), &
                     ungriddedUbound=(/fldList%ungridded_ubound/), gridToFieldMap=(/2/), rc=rc)
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
       fldList => fldList%next

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

  end subroutine dshr_fldlist_realize

end module dshr_fldlist_mod
