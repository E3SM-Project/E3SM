module drof_comp_mod

  use NUOPC                 , only : NUOPC_Advertise
  use ESMF                  , only : ESMF_State, ESMF_SUCCESS, ESMF_STATE, ESMF_LOGMSG_INFO, ESMF_LogWrite, ESMF_Mesh
  use perf_mod              , only : t_startf, t_stopf, t_adj_detailf, t_barrierf
  use shr_kind_mod          , only : r8=>shr_kind_r8, cxx=>shr_kind_cxx, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod           , only : shr_sys_abort
  use shr_strdata_mod       , only : shr_strdata_type, shr_strdata_advance
  use dshr_methods_mod      , only : chkerr, state_getfldptr
  use dshr_dfield_mod       , only : dfield_type, dshr_dfield_add, dshr_dfield_copy
  use dshr_fldlist_mod      , only : fldlist_type, dshr_fldlist_add, dshr_fldlist_realize
  use dshr_nuopc_mod        , only : dshr_get_griddata

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: drof_comp_advertise
  public :: drof_comp_realize
  public :: drof_comp_run

  !--------------------------------------------------------------------------
  ! Module data
  !--------------------------------------------------------------------------

  ! linked lists
  type(fldList_type) , pointer :: fldsImport => null()
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  ! module pointer arrays
  real(r8), pointer :: Forr_rofl(:) => null()
  real(r8), pointer :: Forr_rofi(:) => null()

  character(*)    , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine drof_comp_advertise(importState, exportState, flds_scalar_name, rof_prognostic, rc)

    ! determine export and import fields to advertise to mediator

    ! input/output arguments
    type(ESMF_State)               :: importState
    type(ESMF_State)               :: exportState
    character(len=*) , intent(in)  :: flds_scalar_name
    logical          , intent(in)  :: rof_prognostic
    integer          , intent(out) :: rc

    ! local variables
    integer :: n
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !-------------------
    ! Advrtise export fields
    !-------------------

    call dshr_fldList_add(fldsExport, trim(flds_scalar_name))
    call dshr_fldlist_add(fldsExport, "Forr_rofl")
    call dshr_fldlist_add(fldsExport, "Forr_rofi")

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(drof_comp_advertise): Fr_rof '//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

    ! currently there is no import state to drof

  end subroutine drof_comp_advertise

  !===============================================================================

  subroutine drof_comp_realize(sdat, importState, exportState, flds_scalar_name, flds_scalar_num, mesh, &
       logunit, masterproc, rc)

    !input/output variables
    type(shr_strdata_type) , intent(inout) :: sdat
    type(ESMF_State)       , intent(inout) :: importState
    type(ESMF_State)       , intent(inout) :: exportState
    character(len=*)       , intent(in)    :: flds_scalar_name
    integer                , intent(in)    :: flds_scalar_num
    type(ESMF_Mesh)        , intent(in)    :: mesh
    integer                , intent(in)    :: logunit
    logical                , intent(in)    :: masterproc
    integer                , intent(out)   :: rc

    ! local variables
    character(*), parameter   :: subName = "(drof_comp_realize) "
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! -------------------------------------
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    ! -------------------------------------

    call dshr_fldlist_realize( exportState, fldsExport, flds_scalar_name, flds_scalar_num,  mesh, &
         subname//':drofExport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -------------------------------------
    ! Set pointers to exportState fields
    ! -------------------------------------

    call dshr_dfield_add(dfields, sdat, state_fld='Forr_rofl', strm_fld='rofl', &
         state=exportState, state_ptr=Forr_rofl, logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat, state_fld='Forr_rofi', strm_fld='rofi', &
         state=exportState, state_ptr=Forr_rofi, logunit=logunit, masterproc=masterproc, rc=rc)

  end subroutine drof_comp_realize

  !===============================================================================

  subroutine drof_comp_run(mpicom, my_task, master_task, logunit, target_ymd, target_tod, sdat, rc)

    ! --------------------------
    ! advance drof
    ! --------------------------

    ! input/output variables:
    integer                , intent(in)    :: mpicom           ! mpi communicator
    integer                , intent(in)    :: my_task
    integer                , intent(in)    :: master_task
    integer                , intent(in)    :: logunit
    integer                , intent(in)    :: target_ymd       ! model date
    integer                , intent(in)    :: target_tod       ! model sec into model date
    type(shr_strdata_type) , intent(inout) :: sdat
    integer                , intent(out)   :: rc

    ! local variables
    integer :: n
    character(*), parameter :: subName = "(drof_comp_run) "
    !-------------------------------------------------------------------------------

    call t_startf('DROF_RUN')

    !--------------------
    ! advance drof streams
    !--------------------

    ! time and spatially interpolate to model time and grid
    call t_barrierf('drof_BARRIER',mpicom)
    call t_startf('drof_strdata_advance')
    call shr_strdata_advance(sdat, target_ymd, target_tod, mpicom, 'drof')
    call t_stopf('drof_strdata_advance')

    !--------------------
    ! copy all fields from streams to export state as default
    !--------------------

    ! This automatically will update the fields in the export state
    call t_barrierf('drof_comp_dfield_copy_BARRIER', mpicom)
    call t_startf('drof_dfield_copy')
    call dshr_dfield_copy(dfields,  sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('drof_dfield_copy')

    !-------------------------------------------------
    ! determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('drof_datamode')
    select case (trim(sdat%datamode))
    case('COPYALL')
       ! zero out "special values" of export fields
       do n = 1, size(Forr_rofl)
          if (abs(Forr_rofl(n)) > 1.0e28) Forr_rofl(n) = 0.0_r8
          if (abs(Forr_rofi(n)) > 1.0e28) Forr_rofi(n) = 0.0_r8
       enddo
    end select
    call t_stopf('drof_datamode')

    call t_stopf('DROF_RUN')

  end subroutine drof_comp_run

end module drof_comp_mod
