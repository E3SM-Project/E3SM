module eatm_comp_mod

  !---------------------------------------------------------------------------
  ! Emulator ATM component module.
  !
  ! Provides:
  !   - iso_c_binding interfaces to C++ eatm_binding functions
  !   - eatm_setup_domain: MCT domain initialization from C++ grid data
  !---------------------------------------------------------------------------

  use iso_c_binding
  use shr_kind_mod, only: R8 => SHR_KIND_R8, IN => SHR_KIND_IN
  use mct_mod
  use seq_flds_mod, only: seq_flds_dom_coord, seq_flds_dom_other

  implicit none
  private

  ! Public helper
  public :: eatm_setup_domain

  !---------------------------------------------------------------------------
  ! C++ binding interfaces (eatm_binding.cpp)
  !---------------------------------------------------------------------------
  interface

     subroutine eatm_create_c(fcomm, compid, nxg, nyg) &
          bind(C, name="eatm_create_c")
       import :: c_int
       integer(c_int), value, intent(in) :: fcomm
       integer(c_int), value, intent(in) :: compid
       integer(c_int), value, intent(in) :: nxg
       integer(c_int), value, intent(in) :: nyg
     end subroutine eatm_create_c

     function eatm_get_lsize_c() result(lsize) &
          bind(C, name="eatm_get_lsize_c")
       import :: c_int
       integer(c_int) :: lsize
     end function eatm_get_lsize_c

     subroutine eatm_get_grid_data_c(gindex, lons, lats, &
          areas, masks, fracs) &
          bind(C, name="eatm_get_grid_data_c")
       import :: c_int, c_double
       integer(c_int),  intent(out) :: gindex(*)
       real(c_double),  intent(out) :: lons(*)
       real(c_double),  intent(out) :: lats(*)
       real(c_double),  intent(out) :: areas(*)
       real(c_double),  intent(out) :: masks(*)
       real(c_double),  intent(out) :: fracs(*)
     end subroutine eatm_get_grid_data_c

     subroutine eatm_init_c() bind(C, name="eatm_init_c")
     end subroutine eatm_init_c

     subroutine eatm_run_c(dt) bind(C, name="eatm_run_c")
       import :: c_int
       integer(c_int), value, intent(in) :: dt
     end subroutine eatm_run_c

     subroutine eatm_final_c() bind(C, name="eatm_final_c")
     end subroutine eatm_final_c

  end interface

contains

  !---------------------------------------------------------------------------
  ! eatm_setup_domain
  !
  ! Initialize an MCT domain (mct_ggrid) from arrays provided by C++.
  ! Follows the pattern in dead_mct_mod.F90::dead_domain_mct.
  !---------------------------------------------------------------------------
  subroutine eatm_setup_domain(mpicom, gsMap, domain, lsize, &
       gindex, lons, lats, areas, masks, fracs)

    integer(IN),     intent(in)    :: mpicom
    type(mct_gsMap), intent(in)    :: gsMap
    type(mct_ggrid), intent(inout) :: domain
    integer(IN),     intent(in)    :: lsize
    integer(IN),     intent(in)    :: gindex(lsize)
    real(R8),        intent(in)    :: lons(lsize)
    real(R8),        intent(in)    :: lats(lsize)
    real(R8),        intent(in)    :: areas(lsize)
    real(R8),        intent(in)    :: masks(lsize)
    real(R8),        intent(in)    :: fracs(lsize)

    ! locals
    integer(IN)          :: n, my_task, ier
    real(R8), allocatable    :: rdata(:)
    integer(IN), allocatable :: idata(:)

    ! Initialize the mct_gGrid with coordinate and other attributes
    call mct_gGrid_init(GGrid=domain, &
         CoordChars=trim(seq_flds_dom_coord), &
         OtherChars=trim(seq_flds_dom_other), &
         lsize=lsize)
    call mct_aVect_zero(domain%data)

    ! Import global grid indices
    allocate(idata(lsize))
    idata(:) = gindex(:)
    call mct_gGrid_importIAttr(domain, 'GlobGridNum', idata, lsize)
    deallocate(idata)

    ! Import real grid attributes
    allocate(rdata(lsize))

    rdata(:) = lats(:)
    call mct_gGrid_importRAttr(domain, "lat", rdata, lsize)

    rdata(:) = lons(:)
    call mct_gGrid_importRAttr(domain, "lon", rdata, lsize)

    rdata(:) = areas(:)
    call mct_gGrid_importRAttr(domain, "area",  rdata, lsize)
    call mct_gGrid_importRAttr(domain, "aream", rdata, lsize)

    rdata(:) = masks(:)
    call mct_gGrid_importRAttr(domain, "mask", rdata, lsize)

    rdata(:) = fracs(:)
    call mct_gGrid_importRAttr(domain, "frac", rdata, lsize)

    deallocate(rdata)

  end subroutine eatm_setup_domain

end module eatm_comp_mod
