#include "refining_remapper_rma.hpp"

#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/ekat_pack_utils.hpp>

#include <numeric>

namespace scream
{

RefiningRemapperRMA::
RefiningRemapperRMA (const grid_ptr_type& tgt_grid,
                     const std::string& map_file)
 : HorizInterpRemapperBase(tgt_grid,map_file,InterpType::Refine)
{
  // Nothing to do here
}

RefiningRemapperRMA::
~RefiningRemapperRMA ()
{
  clean_up();
}

void RefiningRemapperRMA::remap_fwd_impl ()
{
  // Start RMA epoch on each field
  for (int i=0; i<m_num_fields; ++i) {
    check_mpi_call(MPI_Win_post(m_mpi_group,0,m_mpi_win[i]),
                   "MPI_Win_post for field: " + m_src_fields[i].name());
    check_mpi_call(MPI_Win_start(m_mpi_group,0,m_mpi_win[i]),
                   "MPI_Win_start for field: " + m_src_fields[i].name());
  }

  // Loop over fields, and grab data
  constexpr HostOrDevice MpiDev = MpiOnDev ? Device : Host;
  const auto& dt = ekat::get_mpi_type<Real>();
  for (int i=0; i<m_num_fields; ++i) {
    const int col_size = m_col_size[i];
    const int col_stride = m_col_stride[i];
    const int col_offset = m_col_offset[i];
    const auto& win = m_mpi_win[i];
    auto ov_data = m_ov_fields[i].get_internal_view_data<Real,MpiDev>();
    for (int icol=0; icol<m_ov_coarse_grid->get_num_local_dofs(); ++icol) {
      const int pid = m_remote_pids[icol];
      const int lid = m_remote_lids[icol];
      check_mpi_call(MPI_Get(ov_data+icol*col_size,col_size,dt,pid,
                             lid*col_stride+col_offset,col_size,dt,win),
                     "MPI_Get for field: " + m_ov_fields[i].name());
    }
  }

  // Close access RMA epoch on each field (exposure is still open)
  for (int i=0; i<m_num_fields; ++i) {
    check_mpi_call(MPI_Win_complete(m_mpi_win[i]),
                   "MPI_Win_complete for field: " + m_ov_fields[i].name());
  }

  // Helpef function, to establish if a field can be handled with packs
  auto can_pack_field = [](const Field& f) {
    const auto& ap = f.get_header().get_alloc_properties();
    return (ap.get_last_extent() % SCREAM_PACK_SIZE) == 0;
  };

  // Loop over each field, perform mat-vec
  constexpr auto COL = ShortFieldTagsNames::COL;
  for (int i=0; i<m_num_fields; ++i) {
    auto& f_tgt = m_tgt_fields[i];

    // Allow to register fields that do not have the COL tag
    // These fields are simply copied from src to tgt.
    if (not f_tgt.get_header().get_identifier().get_layout().has_tag(COL)) {
      f_tgt.deep_copy(m_src_fields[i]);
      continue;
    }

    // Perform the local mat-vec. Recall that in these y=Ax products,
    // x is the overlapped src field, and y is the tgt field.
    const auto& f_ov_src    = m_ov_fields[i];

    // If possible, dispatch kernel with SCREAM_PACK_SIZE
    if (can_pack_field(f_ov_src) and can_pack_field(f_tgt)) {
      local_mat_vec<SCREAM_PACK_SIZE>(f_ov_src,f_tgt);
    } else {
      local_mat_vec<1>(f_ov_src,f_tgt);
    }
  }

  // Close exposure RMA epoch on each field
  for (int i=0; i<m_num_fields; ++i) {
    check_mpi_call(MPI_Win_wait(m_mpi_win[i]),
                   "MPI_Win_post for field: " + m_src_fields[i].name());
  }
}

void RefiningRemapperRMA::setup_mpi_data_structures ()
{
  using namespace ShortFieldTagsNames;

  // Extract some raw mpi info
  const auto mpi_comm  = m_comm.mpi_comm();
  check_mpi_call(MPI_Comm_group(mpi_comm,&m_mpi_group),"MPI_Comm_group");

  // Figure out where data needs to be retrieved from
  const auto ov_src_gids = m_ov_coarse_grid->get_dofs_gids().get_view<const gid_type*,Host>();
  m_src_grid->get_remote_pids_and_lids(ov_src_gids,m_remote_pids,m_remote_lids);

  // TODO: scope out possibility of using sub-groups for start/post calls
  //       (but I'm afraid you can't, b/c start/post may require same groups)

  // Create per-field structures
  m_mpi_win.resize(m_num_fields);
  m_col_size.resize(m_num_fields);
  m_col_stride.resize(m_num_fields);
  m_col_offset.resize(m_num_fields,0);
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f = m_src_fields[i];
    const auto& fh = f.get_header();
    const auto& fap = fh.get_alloc_properties();
    const auto& layout = fh.get_identifier().get_layout();
    m_col_size[i]   = layout.clone().strip_dim(COL).size();

    const int col_stride = m_col_stride[i] = fap.get_num_scalars() / layout.dim(COL);

    // If field has a parent, col_stride and col_offset need to be adjusted
    auto p = fh.get_parent().lock();
    auto win_size = layout.size()*sizeof(Real);
    if (p) {
      EKAT_REQUIRE_MSG (p->get_parent().lock()==nullptr,
          "Error! We do not support remapping of subfields of other subfields.\n");
      const auto& sv_info = fap.get_subview_info();
      m_col_offset[i] = sv_info.slice_idx  * col_stride;
      m_col_stride[i] = sv_info.dim_extent * col_stride;
      win_size *= sv_info.dim_extent;
    }

    auto data = f.get_internal_view_data<Real,Host>();
    check_mpi_call(MPI_Win_create(data,win_size,sizeof(Real),
                                  MPI_INFO_NULL,mpi_comm,&m_mpi_win[i]),
                   "[RefiningRemapperRMA::setup_mpi_data_structures] MPI_Win_create");
#ifndef EKAT_MPI_ERRORS_ARE_FATAL
    check_mpi_call(MPI_Win_set_errhandler(m_mpi_win[i],MPI_ERRORS_RETURN),
                   "[RefiningRemapperRMA::setup_mpi_data_structure] setting MPI_ERRORS_RETURN handler on MPI_Win");
#endif
  }
}

void RefiningRemapperRMA::clean_up ()
{
  // Clear all MPI related structures
  if (m_mpi_group!=MPI_GROUP_NULL) {
    check_mpi_call(MPI_Group_free(&m_mpi_group),"MPI_Group_free");
    m_mpi_group = MPI_GROUP_NULL;
  }
  for (auto& win : m_mpi_win) {
    check_mpi_call(MPI_Win_free(&win),"MPI_Win_free");
  }
  m_mpi_win.clear();
  m_remote_pids.clear();
  m_remote_lids.clear();
  m_col_size.clear();

  HorizInterpRemapperBase::clean_up();
}

} // namespace scream
