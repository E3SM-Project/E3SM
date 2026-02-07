module eocn_comp_mod

   !---------------------------------------------------------------------------
   ! Emulator OCN component module.
   !
   ! Provides:
   !   - iso_c_binding interfaces to C++ eocn_binding functions
   !   - eocn_setup_domain: MCT domain initialization
   !
   ! Grid-agnostic: decomposition computed here (F90) and passed to C++.
   !---------------------------------------------------------------------------

   use iso_c_binding
   use shr_kind_mod, only: R8 => SHR_KIND_R8, IN => SHR_KIND_IN
   use mct_mod
   use seq_flds_mod, only: seq_flds_dom_coord, seq_flds_dom_other

   implicit none
   private

   ! Public helpers
   public :: eocn_setup_domain
   public :: eocn_compute_grid

   !---------------------------------------------------------------------------
   ! C++ binding interfaces (eocn_binding.cpp)
   !---------------------------------------------------------------------------
   interface

      subroutine eocn_create_c(fcomm, compid) &
         bind(C, name="eocn_create_c")
         import :: c_int
         integer(c_int), value, intent(in) :: fcomm
         integer(c_int), value, intent(in) :: compid
      end subroutine eocn_create_c

      subroutine eocn_set_decomposition_c(lsize, gindex) &
         bind(C, name="eocn_set_decomposition_c")
         import :: c_int
         integer(c_int), value, intent(in) :: lsize
         integer(c_int), intent(in) :: gindex(*)
      end subroutine eocn_set_decomposition_c

      function eocn_get_lsize_c() result(lsize) &
         bind(C, name="eocn_get_lsize_c")
         import :: c_int
         integer(c_int) :: lsize
      end function eocn_get_lsize_c

      subroutine eocn_get_gindex_c(gindex) &
         bind(C, name="eocn_get_gindex_c")
         import :: c_int
         integer(c_int), intent(out) :: gindex(*)
      end subroutine eocn_get_gindex_c

      subroutine eocn_init_c() bind(C, name="eocn_init_c")
      end subroutine eocn_init_c

      subroutine eocn_run_c(dt) bind(C, name="eocn_run_c")
         import :: c_int
         integer(c_int), value, intent(in) :: dt
      end subroutine eocn_run_c

      subroutine eocn_final_c() bind(C, name="eocn_final_c")
      end subroutine eocn_final_c

   end interface

contains

   !---------------------------------------------------------------------------
   ! eocn_compute_grid
   !
   ! Compute a simple contiguous 1-D decomposition and synthetic grid.
   ! Ocean uses same logic as atm for now.
   !---------------------------------------------------------------------------
   subroutine eocn_compute_grid(mpicom, nxg, nyg, lsize, &
      gindex, lons, lats, areas, masks, fracs)

      integer(IN), intent(in)  :: mpicom
      integer(IN), intent(in)  :: nxg, nyg
      integer(IN), intent(out) :: lsize
      integer(IN), allocatable, intent(out) :: gindex(:)
      real(R8), allocatable, intent(out) :: lons(:), lats(:)
      real(R8), allocatable, intent(out) :: areas(:), masks(:), fracs(:)

      ! locals
      integer(IN) :: rank, nprocs, ierr
      integer(IN) :: gsize, base, extra, offset, n, g, ig, jg
      real(R8), parameter :: pi = 3.14159265358979323846_R8
      real(R8), parameter :: deg2rad = pi / 180.0_R8
      real(R8), parameter :: re = 6.37122e6_R8  ! Earth radius (m)
      real(R8) :: dx, ys, yc, yn, dy

      call mpi_comm_rank(mpicom, rank, ierr)
      call mpi_comm_size(mpicom, nprocs, ierr)

      gsize = nxg * nyg

      ! Simple contiguous 1-D decomposition
      base = gsize / nprocs
      extra = mod(gsize, nprocs)
      if (rank < extra) then
         lsize = base + 1
      else
         lsize = base
      endif
      offset = rank * base + min(rank, extra)

      ! Allocate arrays
      allocate(gindex(lsize))
      allocate(lons(lsize), lats(lsize), areas(lsize))
      allocate(masks(lsize), fracs(lsize))

      ! Fill synthetic grid data
      dx = 360.0_R8 / real(nxg, R8) * deg2rad

      do n = 1, lsize
         g = offset + n - 1  ! 0-based global index
         ig = mod(g, nxg)    ! 0-based column
         jg = g / nxg        ! 0-based row

         ys = -90.0_R8 + real(jg, R8) * 180.0_R8 / real(nyg, R8)
         yc = -90.0_R8 + (real(jg, R8) + 0.5_R8) * 180.0_R8 / real(nyg, R8)
         yn = -90.0_R8 + (real(jg, R8) + 1.0_R8) * 180.0_R8 / real(nyg, R8)
         dy = sin(yn * deg2rad) - sin(ys * deg2rad)

         gindex(n) = g + 1   ! 1-based for MCT
         lons(n) = real(ig, R8) * 360.0_R8 / real(nxg, R8)
         lats(n) = yc
         areas(n) = dx * dy * re * re
         masks(n) = 1.0_R8
         fracs(n) = 1.0_R8
      end do

   end subroutine eocn_compute_grid

   !---------------------------------------------------------------------------
   ! eocn_setup_domain
   !
   ! Initialize an MCT domain (mct_ggrid) from arrays.
   !---------------------------------------------------------------------------
   subroutine eocn_setup_domain(mpicom, gsMap, domain, lsize, &
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

   end subroutine eocn_setup_domain

end module eocn_comp_mod
