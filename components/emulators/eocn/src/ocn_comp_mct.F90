module ocn_comp_mct

   !---------------------------------------------------------------------------
   ! Emulator OCN component — MCT interface module.
   !
   ! Provides the three entry points the MCT driver expects:
   !   ocn_init_mct, ocn_run_mct, ocn_final_mct
   !
   ! Grid-agnostic design: grid computed in F90 and passed to C++.
   !---------------------------------------------------------------------------

   use esmf
   use mct_mod
   use seq_cdata_mod,    only: seq_cdata, seq_cdata_setptrs
   use seq_infodata_mod, only: seq_infodata_type, &
      seq_infodata_PutData, &
      seq_infodata_GetData
   use seq_flds_mod,     only: seq_flds_o2x_fields, seq_flds_x2o_fields
   use shr_kind_mod,     only: R8 => SHR_KIND_R8, IN => SHR_KIND_IN, &
      CS => SHR_KIND_CS, CL => SHR_KIND_CL
   use iso_c_binding

   use eocn_comp_mod,    only: eocn_create_c,             &
      eocn_set_decomposition_c,  &
      eocn_init_c,               &
      eocn_run_c,                &
      eocn_final_c,              &
      eocn_compute_grid,         &
      eocn_setup_domain

   implicit none
   private

   public :: ocn_init_mct
   public :: ocn_run_mct
   public :: ocn_final_mct

   ! Module-level state
   integer(IN) :: mpicom       ! MPI communicator
   integer(IN) :: my_task      ! MPI rank
   integer(IN) :: compid       ! component id
   integer(IN) :: nxg          ! global grid nx
   integer(IN) :: nyg          ! global grid ny

   integer(IN), parameter :: master_task = 0

contains

   !===========================================================================
   ! ocn_init_mct
   !===========================================================================
   subroutine ocn_init_mct(EClock, cdata, x2o, o2x, NLFilename)

      type(ESMF_Clock),           intent(inout) :: EClock
      type(seq_cdata),            intent(inout) :: cdata
      type(mct_aVect),            intent(inout) :: x2o, o2x
      character(len=*), optional, intent(in)    :: NLFilename

      ! locals
      type(seq_infodata_type), pointer :: infodata
      type(mct_gsMap),         pointer :: gsMap
      type(mct_gGrid),         pointer :: ggrid
      integer(IN)                      :: lsize, ierr
      integer(IN), allocatable         :: gindex(:)
      real(R8),    allocatable         :: lons(:), lats(:)
      real(R8),    allocatable         :: areas(:), masks(:), fracs(:)

      character(*), parameter :: subName = "(ocn_init_mct) "

      !--- extract pointers from cdata ---
      call seq_cdata_setptrs(cdata,               &
         id=compid, mpicom=mpicom,              &
         gsMap=gsMap, dom=ggrid,                &
         infodata=infodata)

      call mpi_comm_rank(mpicom, my_task, ierr)

      !--- grid dimensions (hardcoded for skeleton; could read from xml) ---
      nxg = 48
      nyg = 24

      !--- create the C++ Ocn emulator ---
      call eocn_create_c(int(mpicom, c_int), int(compid, c_int))

      !--- compute grid in F90 ---
      call eocn_compute_grid(mpicom, nxg, nyg, lsize, &
         gindex, lons, lats, areas, masks, fracs)

      !--- pass decomposition to C++ ---
      call eocn_set_decomposition_c(int(lsize, c_int), gindex)

      !--- initialize gsMap (global-to-local index mapping) ---
      call mct_gsMap_init(gsMap, gindex, mpicom, compid, &
         lsize, nxg * nyg)

      !--- initialize MCT domain ---
      call eocn_setup_domain(mpicom, gsMap, ggrid, lsize, &
         gindex, lons, lats, areas, masks, fracs)

      !--- initialize attribute vectors ---
      call mct_aVect_init(o2x, rList=seq_flds_o2x_fields, lsize=lsize)
      call mct_aVect_zero(o2x)

      call mct_aVect_init(x2o, rList=seq_flds_x2o_fields, lsize=lsize)
      call mct_aVect_zero(x2o)

      !--- call C++ emulator initialize ---
      call eocn_init_c()

      !--- tell the coupler we are present and active ---
      call seq_infodata_PutData(infodata, &
         ocn_present   = .true.,        &
         ocn_prognostic = .true.,       &
         ocn_nx = nxg,                  &
         ocn_ny = nyg)

      deallocate(gindex, lons, lats, areas, masks, fracs)

   end subroutine ocn_init_mct

   !===========================================================================
   ! ocn_run_mct
   !===========================================================================
   subroutine ocn_run_mct(EClock, cdata, x2o, o2x)

      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata),  intent(inout) :: cdata
      type(mct_aVect),  intent(inout) :: x2o
      type(mct_aVect),  intent(inout) :: o2x

      ! Skeleton: delegate to C++, o2x fields stay zeroed for now
      call eocn_run_c(int(0, c_int))

   end subroutine ocn_run_mct

   !===========================================================================
   ! ocn_final_mct
   !===========================================================================
   subroutine ocn_final_mct(EClock, cdata, x2o, o2x)

      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata),  intent(inout) :: cdata
      type(mct_aVect),  intent(inout) :: x2o
      type(mct_aVect),  intent(inout) :: o2x

      call eocn_final_c()

   end subroutine ocn_final_mct

end module ocn_comp_mct
