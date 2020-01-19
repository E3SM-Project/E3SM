module drof_comp_mod

  use NUOPC                 , only : NUOPC_Advertise
  use ESMF                  , only : ESMF_State, ESMF_SUCCESS, ESMF_STATE
  use perf_mod              , only : t_startf, t_stopf, t_adj_detailf, t_barrierf
  use mct_mod               , only : mct_avect_lsize, mct_avect_indexRA
  use shr_kind_mod          , only : r8=>shr_kind_r8, cxx=>shr_kind_cxx, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod           , only : shr_sys_abort
  use shr_strdata_mod       , only : shr_strdata_type  
  use shr_strdata_mod       , only : shr_strdata_advance
  use dshr_methods_mod      , only : chkerr, state_getfldptr
  use dshr_nuopc_mod        , only : fld_list_type, fldsMax, dshr_fld_add 
  use dshr_nuopc_mod        , only : dfield_type, dshr_dfield_add, dshr_streams_copy

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: drof_comp_advertise
  public :: drof_comp_dfields_init
  public :: drof_comp_run

  !--------------------------------------------------------------------------
  ! Module data
  !--------------------------------------------------------------------------

  integer             , public :: fldsToRof_num = 0
  integer             , public :: fldsFrRof_num = 0
  type(fld_list_type) , public :: fldsToRof(fldsMax)
  type(fld_list_type) , public :: fldsFrRof(fldsMax)

  type(dfield_type) :: dfields(fldsMax)
  integer           :: dfields_num

  character(*)    , parameter :: u_FILE_u = &
       __FILE__

  ! module pointer arrays
  real(r8), pointer :: Forr_rofl(:)
  real(r8), pointer :: Forr_rofi(:)

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
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !-------------------
    ! Advrtise export fields
    !-------------------

    fldsFrRof_num=1
    fldsFrRof(1)%stdname = trim(flds_scalar_name)

    call dshr_fld_add("Forr_rofl", fldsFrRof_num, fldsFrRof)
    call dshr_fld_add("Forr_rofi", fldsFrRof_num, fldsFrRof)

    do n = 1,fldsFrRof_num
       call NUOPC_Advertise(exportState, standardName=fldsFrRof(n)%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    ! currently there is no import state to drof

  end subroutine drof_comp_advertise

  !===============================================================================

  subroutine drof_comp_dfields_init(sdat, exportState, rc)

    !input/output variables
    type(shr_strdata_type) , intent(in)     :: sdat
    type(ESMF_State)       , intent(inout ) :: exportState
    integer                , intent(out)    :: rc

    ! local variables
    character(CS), allocatable :: strm_flds(:)
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Note that all pointers are initialized to zero
    call dshr_dfield_add(sdat, exportState, dfields, dfields_num, state_fld='Forr_rofl', strm_fld='rofl', &
         state_ptr=Forr_rofl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(sdat, exportState, dfields, dfields_num, state_fld='Forr_rofi', strm_fld='rofi', &
         state_ptr=Forr_rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine drof_comp_dfields_init

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
    call t_barrierf('drof_comp_streams_copy_BARRIER', mpicom)
    call t_startf('drof_streams_copy')
    call dshr_streams_copy(dfields, dfields_num, sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('drof_streams_copy')

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
