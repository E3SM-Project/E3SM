module dwav_comp_mod

  use NUOPC                 , only : NUOPC_Advertise
  use ESMF                  , only : ESMF_State, ESMF_SUCCESS, ESMF_STATE, ESMF_Mesh
  use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO
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

  public :: dwav_comp_advertise
  public :: dwav_comp_realize
  public :: dwav_comp_run

  !--------------------------------------------------------------------------
  ! Module data
  !--------------------------------------------------------------------------

  ! linked lists
  type(fldList_type) , pointer :: fldsImport => null()
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  ! module constants
  character(*)    , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dwav_comp_advertise(importState, exportState, flds_scalar_name, wav_prognostic, rc)

    ! determine export and import fields to advertise to mediator

    ! input/output arguments
    type(ESMF_State)               :: importState
    type(ESMF_State)               :: exportState
    character(len=*) , intent(in)  :: flds_scalar_name
    logical          , intent(in)  :: wav_prognostic
    integer          , intent(out) :: rc

    ! local variables
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !-------------------
    ! Advertise export fields
    !-------------------

    call dshr_fldList_add(fldsExport, trim(flds_scalar_name))
    call dshr_fldList_add(fldsExport, 'Sw_lamult' )
    call dshr_fldList_add(fldsExport, 'Sw_ustokes')
    call dshr_fldList_add(fldsExport, 'Sw_vstokes')

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(dwav_comp_advertise): Fr_wav '//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

    ! currently there is no import state to dwav

  end subroutine dwav_comp_advertise

  !===============================================================================

  subroutine dwav_comp_realize(sdat, importState, exportState, flds_scalar_name, flds_scalar_num, mesh, &
       logunit, masterproc, rc)

    ! input/output variables
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
    character(CS), allocatable :: strm_flds(:)
    character(*), parameter    :: subName = "(dwav_comp_realize) "
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! -------------------------------------
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    ! -------------------------------------

    call dshr_fldlist_realize( exportState, fldsExport, flds_scalar_name, flds_scalar_num,  mesh, &
         subname//':dwavExport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create stream-> export state mapping

    call dshr_dfield_add(dfields, sdat, state_fld='Sw_lamult', strm_fld='lamult', state=exportstate, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat, state_fld='Sw_ustokes', strm_fld='ustokes', state=exportstate, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat, state_fld='Sw_vstokes', strm_fld='vstokes', state=exportstate, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine dwav_comp_realize

  !===============================================================================

  subroutine dwav_comp_run(mpicom, my_task, master_task, logunit, target_ymd, target_tod, sdat, rc)

    ! --------------------------
    ! advance dwav
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
    !-------------------------------------------------------------------------------

    call t_startf('DWAV_RUN')

    !--------------------
    ! advance dwav streams
    !--------------------

    ! time and spatially interpolate to model time and grid
    call t_barrierf('dwav_BARRIER',mpicom)
    call t_startf('dwav_strdata_advance')
    call shr_strdata_advance(sdat, target_ymd, target_tod, mpicom, 'dwav')
    call t_stopf('dwav_strdata_advance')

    !--------------------
    ! copy all fields from streams to export state as default
    !--------------------

    ! This automatically will update the fields in the export state
    call t_barrierf('dwav_comp_strdata_copy_BARRIER', mpicom)
    call t_startf('dwav_strdata_copy')
    call dshr_dfield_copy(dfields, sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('dwav_strdata_copy')

    !-------------------------------------------------
    ! determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('dwav_datamode')
    select case (trim(sdat%datamode))
    case('COPYALL')
       ! do nothing
    end select
    call t_stopf('dwav_datamode')

    call t_stopf('DWAV_RUN')

  end subroutine dwav_comp_run

end module dwav_comp_mod
