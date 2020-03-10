module dlnd_comp_mod

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
  use glc_elevclass_mod     , only : glc_elevclass_as_string, glc_elevclass_init

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: dlnd_comp_advertise
  public :: dlnd_comp_realize
  public :: dlnd_comp_run

  !--------------------------------------------------------------------------
  ! Module data
  !--------------------------------------------------------------------------

  ! linked lists
  type(fldList_type) , pointer :: fldsImport => null()
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  ! module pointer arrays
  real(r8), pointer           :: lfrac(:)

  ! module constants
  integer                     :: glc_nec

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
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Set module variable
    glc_nec = glc_nec_in

    call glc_elevclass_init(glc_nec)

    !-------------------
    ! Advertise export fields
    !-------------------

    call dshr_fldList_add(fldsExport, trim(flds_scalar_name))
    call dshr_fldlist_add(fldsExport, "Sl_lfrin")

    ! The following puts all of the elevation class fields as an
    ! undidstributed dimension in the export state field - index1 is bare land - and the total number of
    ! elevation classes not equal to bare land go from index2 -> glc_nec+1
    if (glc_nec > 0) then
       call dshr_fldList_add(fldsExport, 'Sl_tsrf_elev'  , ungridded_lbound=1, ungridded_ubound=glc_nec+1)
       call dshr_fldList_add(fldsExport, 'Sl_topo_elev'  , ungridded_lbound=1, ungridded_ubound=glc_nec+1)
       call dshr_fldList_add(fldsExport, 'Flgl_qice_elev', ungridded_lbound=1, ungridded_ubound=glc_nec+1)
    end if

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(dlnd_comp_advertise): Fr_lnd '//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
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

  subroutine dlnd_comp_realize(sdat, importState, exportState, flds_scalar_name, flds_scalar_num, mesh, &
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
    integer                    :: n, lsize
    character(len=2)           :: nec_str
    character(CS), allocatable :: strm_flds(:)
    character(*), parameter    :: subName = "(drof_comp_realize) "
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! -------------------------------------
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    ! -------------------------------------

    call dshr_fldlist_realize( exportState, fldsExport, flds_scalar_name, flds_scalar_num,  mesh, &
         subname//':drofExport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set lfrac in export state - this is not dependent on any streams

    call state_getfldptr(exportState, fldname='Sl_lfrin', fldptr1=lfrac, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    lsize = size(lfrac)
    call dshr_get_griddata(sdat, 'frac', lfrac)

    ! Create stream-> export state mapping

    allocate(strm_flds(0:glc_nec))
    do n = 0,glc_nec
       nec_str = glc_elevclass_as_string(n)
       strm_flds(n) = 'tsrf' // trim(nec_str)
    end do
    call dshr_dfield_add(dfields, sdat, state_fld='Sl_tsrf_elev', strm_flds=strm_flds, state=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 0,glc_nec
       nec_str = glc_elevclass_as_string(n)
       strm_flds(n) = 'topo' // trim(nec_str)
    end do
    call dshr_dfield_add(dfields, sdat, state_fld='Sl_topo_elev', strm_flds=strm_flds, state=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 0,glc_nec
       nec_str = glc_elevclass_as_string(n)
       strm_flds(n) = 'qice' // trim(nec_str)
    end do
    call dshr_dfield_add(dfields, sdat, state_fld='Flgl_qice_elev', strm_flds=strm_flds, state=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine dlnd_comp_realize

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
    call t_barrierf('dlnd_comp_strdata_copy_BARRIER', mpicom)
    call t_startf('dlnd_strdata_copy')
    call dshr_dfield_copy(dfields, sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('dlnd_strdata_copy')

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
