module atm_comp_mct

  !---------------------------------------------------------------------------
  ! Emulator ATM component — MCT interface module.
  !
  ! Provides the three entry points the MCT driver expects:
  !   atm_init_mct, atm_run_mct, atm_final_mct
  !
  ! All real work happens in C++ (AtmEmulator); this module only wires
  ! the MCT coupling infrastructure.
  !---------------------------------------------------------------------------

  use esmf
  use mct_mod
  use seq_cdata_mod,    only: seq_cdata, seq_cdata_setptrs
  use seq_infodata_mod, only: seq_infodata_type, &
                              seq_infodata_PutData, &
                              seq_infodata_GetData
  use seq_flds_mod,     only: seq_flds_a2x_fields, seq_flds_x2a_fields
  use shr_kind_mod,     only: R8 => SHR_KIND_R8, IN => SHR_KIND_IN, &
                               CS => SHR_KIND_CS, CL => SHR_KIND_CL
  use iso_c_binding

  use eatm_comp_mod,    only: eatm_create_c,        &
                              eatm_get_lsize_c,      &
                              eatm_get_grid_data_c,  &
                              eatm_init_c,           &
                              eatm_run_c,            &
                              eatm_final_c,          &
                              eatm_setup_domain

  implicit none
  private

  public :: atm_init_mct
  public :: atm_run_mct
  public :: atm_final_mct

  ! Module-level state
  integer(IN) :: mpicom       ! MPI communicator
  integer(IN) :: my_task      ! MPI rank
  integer(IN) :: compid       ! component id
  integer(IN) :: nxg          ! global grid nx
  integer(IN) :: nyg          ! global grid ny

  integer(IN), parameter :: master_task = 0

contains

  !===========================================================================
  ! atm_init_mct
  !===========================================================================
  subroutine atm_init_mct(EClock, cdata, x2a, a2x, NLFilename)

    type(ESMF_Clock),           intent(inout) :: EClock
    type(seq_cdata),            intent(inout) :: cdata
    type(mct_aVect),            intent(inout) :: x2a, a2x
    character(len=*), optional, intent(in)    :: NLFilename

    ! locals
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap),         pointer :: gsMap
    type(mct_gGrid),         pointer :: ggrid
    integer(IN)                      :: lsize, ierr
    integer(IN), allocatable         :: gindex(:)
    real(R8),    allocatable         :: lons(:), lats(:)
    real(R8),    allocatable         :: areas(:), masks(:), fracs(:)

    character(*), parameter :: subName = "(atm_init_mct) "

    !--- extract pointers from cdata ---
    call seq_cdata_setptrs(cdata,               &
         id=compid, mpicom=mpicom,              &
         gsMap=gsMap, dom=ggrid,                &
         infodata=infodata)

    call mpi_comm_rank(mpicom, my_task, ierr)

    !--- grid dimensions (hardcoded for skeleton; could read from xml) ---
    nxg = 48
    nyg = 24

    !--- create the C++ AtmEmulator and compute decomposition ---
    call eatm_create_c(int(mpicom, c_int), int(compid, c_int), &
         int(nxg, c_int), int(nyg, c_int))

    !--- query local size from C++ ---
    lsize = int(eatm_get_lsize_c(), IN)

    !--- allocate and retrieve grid data from C++ ---
    allocate(gindex(lsize))
    allocate(lons(lsize), lats(lsize), areas(lsize))
    allocate(masks(lsize), fracs(lsize))

    call eatm_get_grid_data_c(gindex, lons, lats, areas, masks, fracs)

    !--- initialize gsMap (global-to-local index mapping) ---
    call mct_gsMap_init(gsMap, gindex, mpicom, compid, &
         lsize, nxg * nyg)

    !--- initialize MCT domain ---
    call eatm_setup_domain(mpicom, gsMap, ggrid, lsize, &
         gindex, lons, lats, areas, masks, fracs)

    !--- initialize attribute vectors ---
    call mct_aVect_init(a2x, rList=seq_flds_a2x_fields, lsize=lsize)
    call mct_aVect_zero(a2x)

    call mct_aVect_init(x2a, rList=seq_flds_x2a_fields, lsize=lsize)
    call mct_aVect_zero(x2a)

    !--- call C++ emulator initialize ---
    call eatm_init_c()

    !--- tell the coupler we are present and active ---
    call seq_infodata_PutData(infodata, &
         atm_present   = .true.,        &
         atm_prognostic = .true.,       &
         atm_nx = nxg,                  &
         atm_ny = nyg)

    deallocate(gindex, lons, lats, areas, masks, fracs)

  end subroutine atm_init_mct

  !===========================================================================
  ! atm_run_mct
  !===========================================================================
  subroutine atm_run_mct(EClock, cdata, x2a, a2x)

    type(ESMF_Clock), intent(inout) :: EClock
    type(seq_cdata),  intent(inout) :: cdata
    type(mct_aVect),  intent(inout) :: x2a
    type(mct_aVect),  intent(inout) :: a2x

    ! Skeleton: delegate to C++, a2x fields stay zeroed for now
    call eatm_run_c(int(0, c_int))

  end subroutine atm_run_mct

  !===========================================================================
  ! atm_final_mct
  !===========================================================================
  subroutine atm_final_mct(EClock, cdata, x2a, a2x)

    type(ESMF_Clock), intent(inout) :: EClock
    type(seq_cdata),  intent(inout) :: cdata
    type(mct_aVect),  intent(inout) :: x2a
    type(mct_aVect),  intent(inout) :: a2x

    call eatm_final_c()

  end subroutine atm_final_mct

end module atm_comp_mct
