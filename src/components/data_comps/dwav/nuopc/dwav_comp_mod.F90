module dwav_comp_mod

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

  public :: dwav_comp_advertise
  public :: dwav_comp_dfields_init
  public :: dwav_comp_run

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  integer             , public :: fldsToWav_num = 0
  integer             , public :: fldsFrWav_num = 0
  type(fld_list_type) , public :: fldsToWav(fldsMax)
  type(fld_list_type) , public :: fldsFrWav(fldsMax)

  type(dfield_type) :: dfields(fldsMax)
  integer           :: dfields_num

  character(len=*), parameter :: rpfile = 'rpointer.wav'
  character(*)    , parameter :: u_FILE_u = &
       __FILE__

  ! module pointer arrays
  real(r8), pointer :: Sw_lamult(:)
  real(r8), pointer :: Sw_ustokes(:)
  real(r8), pointer :: Sw_vstokes(:)

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
    integer :: n
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !-------------------
    ! Advertise export fields
    !-------------------

    fldsFrWav_num=1
    fldsFrWav(1)%stdname = trim(flds_scalar_name)

    call dshr_fld_add('Sw_lamult' , fldsFrWav_num, fldsFrWav)
    call dshr_fld_add('Sw_ustokes', fldsFrWav_num, fldsFrWav)
    call dshr_fld_add('Sw_vstokes', fldsFrWav_num, fldsFrWav)

    do n = 1,fldsFrWav_num
       call NUOPC_Advertise(exportState, standardName=fldsFrWav(n)%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    ! currently there is no import state to dwav

  end subroutine dwav_comp_advertise

  !===============================================================================

  subroutine dwav_comp_dfields_init(sdat, exportState, rc)

    !input/output variables
    type(shr_strdata_type) , intent(in)     :: sdat
    type(ESMF_State)       , intent(inout ) :: exportState
    integer                , intent(out)    :: rc

    ! local variables
    character(CS), allocatable :: strm_flds(:)
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Note that all pointers are initialized to zero
    call dshr_dfield_add(sdat, exportState, dfields, dfields_num, state_fld='Sw_lamult', strm_fld='lamult', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(sdat, exportState, dfields, dfields_num, state_fld='Sw_ustokes', strm_fld='ustokes', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(sdat, exportState, dfields, dfields_num, state_fld='Sw_vstokes', strm_fld='vstokes', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine dwav_comp_dfields_init

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

    ! local variables
    integer :: n
    character(*), parameter :: subName = "(dwav_comp_run) "
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
    call t_barrierf('dwav_comp_streams_copy_BARRIER', mpicom)
    call t_startf('dwav_streams_copy')
    call dshr_streams_copy(dfields, dfields_num, sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('dwav_streams_copy')

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
