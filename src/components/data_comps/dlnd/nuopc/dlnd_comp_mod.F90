module dlnd_comp_mod

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
  use glc_elevclass_mod     , only : glc_elevclass_as_string, glc_elevclass_init

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: dlnd_comp_advertise
  public :: dlnd_comp_init
  public :: dlnd_comp_run

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  integer             , public :: fldsToLnd_num = 0
  integer             , public :: fldsFrLnd_num = 0
  type(fld_list_type) , public :: fldsToLnd(fldsMax)
  type(fld_list_type) , public :: fldsFrLnd(fldsMax)

  type(dfield_type)           :: dfields(fldsMax)
  integer                     :: dfields_num

  integer                     :: glc_nec
  real(r8), pointer           :: lfrac(:)  
  character(len=*), parameter :: rpfile = 'rpointer.lnd'
  character(*)    , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dlnd_comp_advertise(importState, exportState, flds_scalar_name, glc_nec_in, rc)

    ! determine export and import fields to advertise to mediator

    ! input/output arguments
    type(ESMF_State)               :: importState
    type(ESMF_State)               :: exportState
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(in)  :: glc_nec_in
    integer          , intent(out) :: rc

    ! local variables
    integer :: n
    character(len=CS)  :: data_fld_name
    character(len=CS)  :: model_fld_name
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Set module variable
    glc_nec = glc_nec_in

    call glc_elevclass_init(glc_nec)

    !-------------------
    ! Advertise export fields
    !-------------------

    fldsFrLnd_num=1
    fldsFrLnd(1)%stdname = trim(flds_scalar_name)

    call dshr_fld_add("Sl_lfrin", fldlist_num=fldsFrLnd_num, fldlist=fldsFrLnd)

    ! The following puts all of the elevation class fields as an
    ! undidstributed dimension in the export state field - index1 is bare land - and the total number of
    ! elevation classes not equal to bare land go from index2 -> glc_nec+1
    if (glc_nec > 0) then
       call dshr_fld_add('Sl_tsrf_elev'  , fldsFrLnd_num, fldlist=fldsFrLnd, ungridded_lbound=1, ungridded_ubound=glc_nec+1)
       call dshr_fld_add('Sl_topo_elev'  , fldsFrLnd_num, fldlist=fldsFrLnd, ungridded_lbound=1, ungridded_ubound=glc_nec+1)
       call dshr_fld_add('Flgl_qice_elev', fldsFrLnd_num, fldlist=fldsFrLnd, ungridded_lbound=1, ungridded_ubound=glc_nec+1)
    end if

    do n = 1,fldsFrLnd_num
       call NUOPC_Advertise(exportState, standardName=fldsFrLnd(n)%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    ! TODO: Non snow fields that nead to be added if dlnd is in cplhist mode
    ! "Sl_t        " "Sl_tref     " "Sl_qref     " "Sl_avsdr    "
    ! "Sl_anidr    " "Sl_avsdf    " "Sl_anidf    " "Sl_snowh    "
    ! "Fall_taux   " "Fall_tauy   " "Fall_lat    " "Fall_sen    "
    ! "Fall_lwup   " "Fall_evap   " "Fall_swnet  " "Sl_landfrac "
    ! "Sl_fv       " "Sl_ram1     "
    ! "Fall_flxdst1" "Fall_flxdst2" "Fall_flxdst3" "Fall_flxdst4"

  end subroutine dlnd_comp_advertise

  !===============================================================================

  subroutine dlnd_comp_init(sdat, exportState, rc)

    !input/output variables
    type(shr_strdata_type) , intent(in)     :: sdat
    type(ESMF_State)       , intent(inout ) :: exportState
    integer                , intent(out)    :: rc

    ! local variables
    integer :: n, kf
    integer :: lsize 
    character(len=2) :: nec_str
    character(CS), allocatable :: strm_flds(:)
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Set lfrac in export state - this is not dependent on any streams
    lsize = mct_avect_lsize(SDAT%grid%data)
    allocate(lfrac(lsize))
    call state_getfldptr(exportState, fldname='Sl_lfrin', fldptr1=lfrac, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    kf = mct_aVect_indexRA(SDAT%grid%data,'frac')
    lfrac(:) = SDAT%grid%data%rAttr(kf,:)

    ! Create stream-> export state mapping (note that dfields_num is a module variable)
    dfields_num = 0

    allocate(strm_flds(0:glc_nec))
    do n = 0,glc_nec
       nec_str = glc_elevclass_as_string(n)
       strm_flds(n) = 'tsrf' // trim(nec_str)
    end do
    call dshr_dfield_add(sdat, exportState, dfields, dfields_num, state_fld='Sl_tsrf_elev', strm_flds=strm_flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 0,glc_nec
       nec_str = glc_elevclass_as_string(n)
       strm_flds(n) = 'topo' // trim(nec_str)
    end do
    call dshr_dfield_add(sdat, exportState, dfields, dfields_num, state_fld='Sl_topo_elev', strm_flds=strm_flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 0,glc_nec
       nec_str = glc_elevclass_as_string(n)
       strm_flds(n) = 'qice' // trim(nec_str)
    end do
    call dshr_dfield_add(sdat, exportState, dfields, dfields_num, state_fld='Flgl_qice_elev', strm_flds=strm_flds, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


  end subroutine dlnd_comp_init

  !===============================================================================

  subroutine dlnd_comp_run(mpicom, my_task, master_task, logunit, target_ymd, target_tod, sdat, rc)

    ! --------------------------
    ! advance dlnd
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
    character(*), parameter :: subName = "(dlnd_comp_run) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('DLND_RUN')

    !--------------------
    ! advance dlnd streams  
    !--------------------

    ! time and spatially interpolate to model time and grid
    call t_barrierf('dlnd_BARRIER',mpicom)
    call t_startf('dlnd_strdata_advance')
    call shr_strdata_advance(sdat, target_ymd, target_tod, mpicom, 'dlnd')
    call t_stopf('dlnd_strdata_advance')

    !--------------------
    ! copy all fields from streams to export state as default
    !--------------------

    ! This automatically will update the fields in the export state
    call t_barrierf('dlnd_comp_streams_copy_BARRIER', mpicom)
    call t_startf('dlnd_streams_copy')
    call dshr_streams_copy(dfields, dfields_num, sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('dlnd_streams_copy')

    !-------------------------------------------------
    ! determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('dlnd_datamode')
    select case (trim(sdat%datamode))
    case('COPYALL')
       ! do nothing extra
    end select
    call t_stopf('dlnd_datamode')

    call t_stopf('DLND_RUN')

  end subroutine dlnd_comp_run

end module dlnd_comp_mod
