#include "coarsening_remapper.hpp"

#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/ekat_pack_utils.hpp>

#include <numeric>

namespace scream
{

CoarseningRemapper::
CoarseningRemapper (const grid_ptr_type& src_grid,
                    // const grid_ptr_type& tgt_grid,
                    const std::string& map_file)
 : AbstractRemapper(src_grid,create_tgt_grid(map_file,src_grid))
 , m_comm (src_grid->get_comm())
{
  using namespace ShortFieldTagsNames;

  // Sanity checks
  EKAT_REQUIRE_MSG (src_grid->type()==GridType::Point,
      "Error! CoarseningRemapper only works on PointGrid grids.\n"
      "  - src grid name: " + src_grid->name() + "\n"
      "  - src_grid_type: " + e2str(src_grid->type()) + "\n");
  EKAT_REQUIRE_MSG (src_grid->is_unique(),
      "Error! CoarseningRemapper requires a unique source grid.\n");

  // This is a coarsening remapper. We only go in one direction
  m_bwd_allowed = false;

  // Create io_grid, containing the indices of the triplets
  // in the map file that this rank has to read
  auto gids_h = get_my_triplets_gids (map_file);
  view_1d<gid_t> gids_d("",gids_h.size());
  Kokkos::deep_copy(gids_d,gids_h);

  auto io_grid = std::make_shared<PointGrid>("",gids_h.size(),0,m_comm);
  io_grid->set_dofs(gids_d);
  const int nlweights = io_grid->get_num_local_dofs();

  // Create CRS matrix views

  // Read in triplets.
  view_1d<gid_t>::HostMirror  row_gids_h("",nlweights);
  view_1d<gid_t>::HostMirror  col_gids_h("",nlweights);
  view_1d<double>::HostMirror S_h ("",nlweights);

  // scream's gids are of type int, while scorpio wants long int as offsets.
  std::vector<scorpio::offset_t> dofs_offsets(nlweights);
  for (int i=0; i<nlweights; ++i) {
    dofs_offsets[i] = gids_h[i];
  }
  const std::string idx_decomp_tag = "int-nnz" + std::to_string(nlweights);
  const std::string val_decomp_tag = "real-nnz" + std::to_string(nlweights);

  scorpio::register_file(map_file,scorpio::FileMode::Read);
  scorpio::get_variable(map_file, "row", "row", {"n_s"}, "int", idx_decomp_tag);
  scorpio::get_variable(map_file, "col", "col", {"n_s"}, "int", idx_decomp_tag);
  scorpio::get_variable(map_file, "S",   "S",   {"n_s"}, "real", val_decomp_tag);
  scorpio::set_dof(map_file,"row",nlweights,dofs_offsets.data());
  scorpio::set_dof(map_file,"col",nlweights,dofs_offsets.data());
  scorpio::set_dof(map_file,"S",nlweights,dofs_offsets.data());
  scorpio::set_decomp(map_file);
  scorpio::grid_read_data_array(map_file,"row",-1,row_gids_h.data(),nlweights);
  scorpio::grid_read_data_array(map_file,"col",-1,col_gids_h.data(),nlweights);
  scorpio::grid_read_data_array(map_file,"S",  -1,S_h.data(),       nlweights);
  scorpio::eam_pio_closefile(map_file);

  // Create an "overlapped" tgt grid, that is, a grid where each rank
  // owns all tgt rows that are affected by at least one of the cols
  // in its src_grid
  std::set<gid_t> ov_tgt_gids;
  for (int i=0; i<nlweights; ++i) {
    ov_tgt_gids.insert(row_gids_h(i)-1);
  }
  const int num_ov_tgt_gids = ov_tgt_gids.size();
  view_1d<int> ov_tgt_gids_d("",num_ov_tgt_gids);
  auto ov_tgt_gids_h = Kokkos::create_mirror_view(ov_tgt_gids_d);
  auto it = ov_tgt_gids.begin();
  for (int i=0; i<num_ov_tgt_gids; ++i, ++it) {
    ov_tgt_gids_h[i] = *it;
  }
  Kokkos::deep_copy(ov_tgt_gids_d,ov_tgt_gids_h);

  auto ov_tgt_grid = std::make_shared<PointGrid>("ov_tgt_grid",num_ov_tgt_gids,0,m_comm);
  ov_tgt_grid->set_dofs(ov_tgt_gids_d);
  m_ov_tgt_grid = ov_tgt_grid;
  const int num_ov_row_gids = m_ov_tgt_grid->get_num_local_dofs();

  // Now we have to create the weights CRS matrix
  m_row_offsets = view_1d<int>("",num_ov_row_gids+1);
  m_col_lids    = view_1d<int>("",nlweights);
  m_weights     = view_1d<Real>("",nlweights);

  // Sort col_gids_h and row_gids_h by row gid. It is easier to sort
  // the array [0,...,n), and use it later to index the 
  std::vector<int> id (nlweights);
  std::iota(id.begin(),id.end(),0);
  auto compare = [&] (const int i, const int j) -> bool {
    return row_gids_h(i) < row_gids_h(j);
  };
  std::sort(id.begin(),id.end(),compare);

  // Create mirror views
  auto row_offsets_h = Kokkos::create_mirror_view(m_row_offsets);
  auto col_lids_h    = Kokkos::create_mirror_view(m_col_lids);
  auto weights_h     = Kokkos::create_mirror_view(m_weights);

  for (int i=0; i<nlweights; ++i) {
    col_lids_h(i) = gid2lid(col_gids_h(id[i])-1,m_src_grid);
    weights_h(i)  = S_h(id[i]);
  }

  Kokkos::deep_copy(m_weights,weights_h);
  Kokkos::deep_copy(m_col_lids,col_lids_h);

  // Compute row offsets
  std::vector<int> row_counts(num_ov_row_gids);
  for (int i=0; i<nlweights; ++i) {
    ++row_counts[gid2lid(row_gids_h(i)-1,m_ov_tgt_grid)];
  }
  std::partial_sum(row_counts.begin(),row_counts.end(),row_offsets_h.data()+1);
  EKAT_REQUIRE_MSG (
      row_offsets_h(num_ov_row_gids)==nlweights,
      "Error! Something went wrong while computing row offsets.\n"
      "  - local nnz       : " + std::to_string(nlweights) + "\n"
      "  - row_offsets(end): " + std::to_string(row_offsets_h(num_ov_row_gids)) + "\n");

  Kokkos::deep_copy(m_row_offsets,row_offsets_h);
}

CoarseningRemapper::
~CoarseningRemapper ()
{
  // We need to free MPI requests
  for (size_t i=0; i<m_send_req.size(); ++i) {
    MPI_Request_free(&m_send_req[i]);
  }
  for (size_t i=0; i<m_recv_req.size(); ++i) {
    MPI_Request_free(&m_recv_req[i]);
  }
}

FieldLayout CoarseningRemapper::
create_src_layout (const FieldLayout& tgt_layout) const
{
  using namespace ShortFieldTagsNames;
  const auto lt = get_layout_type(tgt_layout.tags());
  auto src = FieldLayout::invalid();
  const bool midpoints = tgt_layout.has_tag(LEV);
  const int vec_dim = tgt_layout.is_vector_layout() ? tgt_layout.get_vector_dim() : -1;
  switch (lt) {
    case LayoutType::Scalar2D:
      src = m_src_grid->get_2d_scalar_layout();
      break;
    case LayoutType::Vector2D:
      src = m_src_grid->get_2d_vector_layout(CMP,vec_dim);
      break;
    case LayoutType::Scalar3D:
      src = m_src_grid->get_3d_scalar_layout(midpoints);
      break;
    case LayoutType::Vector3D:
      src = m_src_grid->get_3d_vector_layout(midpoints,CMP,vec_dim);
      break;
    default:
      EKAT_ERROR_MSG ("Layout not supported by CoarseningRemapper: " + e2str(lt) + "\n");
  }
  return src;
}
FieldLayout CoarseningRemapper::
create_tgt_layout (const FieldLayout& src_layout) const
{
  using namespace ShortFieldTagsNames;
  const auto lt = get_layout_type(src_layout.tags());
  auto tgt = FieldLayout::invalid();
  const bool midpoints = src_layout.has_tag(LEV);
  const int vec_dim = src_layout.is_vector_layout() ? src_layout.get_vector_dim() : -1;
  switch (lt) {
    case LayoutType::Scalar2D:
      tgt = m_tgt_grid->get_2d_scalar_layout();
      break;
    case LayoutType::Vector2D:
      tgt = m_tgt_grid->get_2d_vector_layout(CMP,vec_dim);
      break;
    case LayoutType::Scalar3D:
      tgt = m_tgt_grid->get_3d_scalar_layout(midpoints);
      break;
    case LayoutType::Vector3D:
      tgt = m_tgt_grid->get_3d_vector_layout(midpoints,CMP,vec_dim);
      break;
    default:
      EKAT_ERROR_MSG ("Layout not supported by CoarseningRemapper: " + e2str(lt) + "\n");
  }
  return tgt;
}

CoarseningRemapper::grid_ptr_type
CoarseningRemapper::
create_tgt_grid (const std::string& map_file,
                 const grid_ptr_type& src_grid)
{
  EKAT_REQUIRE_MSG (scorpio::is_eam_pio_subsystem_inited(),
      "Error! PIO subsystem was not yet inited.\n");

  scorpio::register_file(map_file,scorpio::FileMode::Read);

  const int ncols_src = scorpio::get_dimlen_c2f(map_file.c_str(),"n_a");
  const int ncols_tgt = scorpio::get_dimlen_c2f(map_file.c_str(),"n_b");

  // Ensure it is indeed a coarsening, and the src grid size matches the src grid
  EKAT_REQUIRE_MSG (ncols_tgt<src_grid->get_num_global_dofs(),
      "Error! CoarseningRemapper requires tgt grid to have less dofs.\n"
      "  - src grid ncols: " + std::to_string(src_grid->get_num_global_dofs()) + "\n"
      "  - map file tgt ncols: " + std::to_string(ncols_tgt) + "\n");
  EKAT_REQUIRE_MSG (ncols_src==src_grid->get_num_global_dofs(),
      "Error! Invalid src grid ncols in CoarseningRemapper map file.\n"
      "  - src grid ncols: " + std::to_string(src_grid->get_num_global_dofs()) + "\n"
      "  - map file src ncols: " + std::to_string(ncols_src) + "\n");

  const int nlevs  = src_grid->get_num_vertical_levels();
  const auto& comm = src_grid->get_comm();

  auto grid = create_point_grid("",ncols_tgt,nlevs,comm);

  scorpio::eam_pio_closefile(map_file);

  return grid;
}

void CoarseningRemapper::
do_register_field (const identifier_type& src, const identifier_type& tgt)
{
  m_src_fields.push_back(field_type(src));
  m_tgt_fields.push_back(field_type(tgt));
}

void CoarseningRemapper::
do_bind_field (const int ifield, const field_type& src, const field_type& tgt)
{
  EKAT_REQUIRE_MSG (
      src.get_header().get_identifier().get_layout().rank()>1 ||
      src.get_header().get_alloc_properties().get_padding()==0,
      "Error! We don't support 2d scalar fields that are padded.\n");
  EKAT_REQUIRE_MSG (
      tgt.get_header().get_identifier().get_layout().rank()>1 ||
      tgt.get_header().get_alloc_properties().get_padding()==0,
      "Error! We don't support 2d scalar fields that are padded.\n");
  m_src_fields[ifield] = src;
  m_tgt_fields[ifield] = tgt;

  // If this was the last field to be bound, we can setup the MPI schedule
  if (this->m_state==RepoState::Closed &&
      (this->m_num_bound_fields+1)==this->m_num_registered_fields) {
    create_ov_tgt_fields ();
    setup_mpi_data_structures ();
  }
}

void CoarseningRemapper::do_registration_ends ()
{
  if (this->m_num_bound_fields==this->m_num_registered_fields) {
    create_ov_tgt_fields ();
    setup_mpi_data_structures ();
  }
}

void CoarseningRemapper::do_remap_fwd () const
{
  // Fire the recv requests right away, so that if some other ranks
  // is done packing before us, we can start receiving their data
  if (not m_recv_req.empty()) {
    int ierr = MPI_Startall(m_recv_req.size(),m_recv_req_ptr);
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something whent wrong while starting persistent recv requests.\n"
        "  - recv rank: " + std::to_string(m_comm.rank()) + "\n");
  }

  // Loop over each field
  constexpr auto can_pack = SCREAM_PACK_SIZE>1;
  for (int i=0; i<m_num_fields; ++i) {
    // First, perform the local mat-vec. Recall that in these y=Ax products,
    // x is the src field, and y is the overlapped tgt field.
    const auto& f_src    = m_src_fields[i];
    const auto& f_ov_tgt = m_ov_tgt_fields[i];

    // Dispatch kernel with the largest possible pack size
    const auto& src_ap = f_src.get_header().get_alloc_properties();
    const auto& ov_tgt_ap = f_ov_tgt.get_header().get_alloc_properties();
    if (can_pack && src_ap.is_compatible<RPack<16>>() &&
                    ov_tgt_ap.is_compatible<RPack<16>>()) {
      local_mat_vec<16>(f_src,f_ov_tgt);
    } else if (can_pack && src_ap.is_compatible<RPack<8>>() &&
                           ov_tgt_ap.is_compatible<RPack<8>>()) {
      local_mat_vec<8>(f_src,f_ov_tgt);
    } else if (can_pack && src_ap.is_compatible<RPack<4>>() &&
                           ov_tgt_ap.is_compatible<RPack<4>>()) {
      local_mat_vec<4>(f_src,f_ov_tgt);
    } else {
      local_mat_vec<1>(f_src,f_ov_tgt);
    }
  }

  // Pack, then fire off the sends
  pack ();
  if (not m_send_req.empty()) {
    int ierr = MPI_Startall(m_send_req.size(),m_send_req_ptr);
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something whent wrong while starting persistent send requests.\n"
        "  - send rank: " + std::to_string(m_comm.rank()) + "\n");
  }

  // Wait for all data to be received, then unpack
  if (not m_recv_req.empty()) {
    int ierr = MPI_Waitall(m_recv_req.size(),m_recv_req_ptr, MPI_STATUSES_IGNORE);
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something whent wrong while waiting on persistent recv requests.\n"
        "  - recv rank: " + std::to_string(m_comm.rank()) + "\n");
  }
  unpack ();

  // Wait for all sends to be completed
  if (not m_send_req.empty()) {
    int ierr = MPI_Waitall(m_send_req.size(),m_send_req_ptr, MPI_STATUSES_IGNORE);
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something whent wrong while waiting on persistent send requests.\n"
        "  - send rank: " + std::to_string(m_comm.rank()) + "\n");
  }

}

template<int PackSize>
void CoarseningRemapper::
local_mat_vec (const Field& x, const Field& y) const
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  using Pack        = ekat::Pack<Real,PackSize>;
  using PackInfo    = ekat::PackInfo<PackSize>;

  const auto& src_layout = x.get_header().get_identifier().get_layout();
  const int rank = src_layout.rank();
  const int nrows = m_ov_tgt_grid->get_num_local_dofs();
  auto row_offsets = m_row_offsets;
  auto col_lids = m_col_lids;
  auto weights = m_weights;
  switch (rank) {
    // Note: in each case, handle 1st contribution to each row separately,
    //       using = instead of +=. This allows to avoid doing an extra
    //       loop to zero out y before the mat-vec.
    case 1:
    {
      auto x_view = x.get_view<const Real*>();
      auto y_view = y.get_view<      Real*>();
      Kokkos::parallel_for(RangePolicy(0,nrows),
                           KOKKOS_LAMBDA(const int& row) {
        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        y_view(row) = weights(beg)*x_view(col_lids(beg));
        for (int icol=beg+1; icol<end; ++icol) {
          y_view(row) += weights(icol)*x_view(col_lids(icol));
        }
      });
      break;
    }
    case 2:
    {
      auto x_view = x.get_view<const Pack**>();
      auto y_view = y.get_view<      Pack**>();
      const int dim1 = PackInfo::num_packs(src_layout.dim(1));
      auto policy = ESU::get_default_team_policy(nrows,dim1);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto row = team.league_rank();

        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team,dim1),
                            [&](const int j){
          y_view(row,j) = weights(beg)*x_view(col_lids(beg),j);
          for (int icol=beg+1; icol<end; ++icol) {
            y_view(row,j) += weights(icol)*x_view(col_lids(icol),j);
          }
        });
      });
      break;
    }
    case 3:
    {
      auto x_view = x.get_view<const Pack***>();
      auto y_view = y.get_view<      Pack***>();
      const int dim1 = src_layout.dim(1);
      const int dim2 = PackInfo::num_packs(src_layout.dim(2));
      auto policy = ESU::get_default_team_policy(nrows,dim1*dim2);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto row = team.league_rank();

        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team,dim1*dim2),
                            [&](const int idx){
          const int j = idx / dim2;
          const int k = idx % dim2;
          y_view(row,j,k) = weights(beg)*x_view(col_lids(beg),j,k);
          for (int icol=beg+1; icol<end; ++icol) {
            y_view(row,j,k) += weights(icol)*x_view(col_lids(icol),j,k);
          }
        });
      });
      break;
    }
  }
}

void CoarseningRemapper::pack () const
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  const int num_send_gids = m_ov_tgt_grid->get_num_local_dofs();
  const auto pid_lid_start = m_send_pid_lids_start;
  const auto lids_pids = m_send_lids_pids;
  const auto buf = m_send_buffer;

  for (int ifield=0; ifield<m_num_fields; ++ifield) {
    const auto& f  = m_ov_tgt_fields[ifield];
    const auto& fl = f.get_header().get_identifier().get_layout();
    const auto lt = get_layout_type(fl.tags());
    const auto f_pid_offsets = ekat::subview(m_send_f_pid_offsets,ifield);

    switch (lt) {
      case LayoutType::Scalar2D:
      {
        auto v = f.get_view<const Real*>();
        Kokkos::parallel_for(RangePolicy(0,num_send_gids),
                             KOKKOS_LAMBDA(const int& i){
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          buf (offset + lidpos) = v(lid);
        });
                             
      } break;
      case LayoutType::Vector2D:
      {
        auto v = f.get_view<const Real**>();
        const int ndims = fl.dim(1);
        auto policy = ESU::get_default_team_policy(num_send_gids,ndims);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int i = team.league_rank();
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team,ndims),
                               [&](const int idim) {
            buf(offset + lidpos*ndims + idim) = v(lid,idim);
          });
        });
      } break;
      case LayoutType::Scalar3D:
      {
        auto v = f.get_view<const Real**>();
        const int nlevs = fl.dims().back();
        auto policy = ESU::get_default_team_policy(num_send_gids,nlevs);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int i = team.league_rank();
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team,nlevs),
                               [&](const int ilev) {
            buf(offset + lidpos*nlevs + ilev) = v(lid,ilev);
          });
        });
      } break;
      case LayoutType::Vector3D:
      {
        auto v = f.get_view<const Real***>();
        const int ndims = fl.dim(1);
        const int nlevs = fl.dims().back();
        auto policy = ESU::get_default_team_policy(num_send_gids,ndims*nlevs);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int i = team.league_rank();
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team,ndims*nlevs),
                               [&](const int idx) {
            const int idim = idx / nlevs;
            const int ilev = idx % nlevs;
            buf(offset + lidpos*ndims*nlevs + idim*nlevs + ilev) = v(lid,idim,ilev);
          });
        });
      } break;

      default:
        EKAT_ERROR_MSG ("Unexpected field rank in CoarseningRemapper::pack.\n"
            "  - MPI rank  : " + std::to_string(m_comm.rank()) + "\n"
            "  - field rank: " + std::to_string(fl.rank()) + "\n");
    }
  }

  // If MPI does not use dev pointers, we need to deep copy from dev to host
  if (not MpiOnDev) {
    Kokkos::deep_copy (m_mpi_send_buffer,m_send_buffer);
  }
}

void CoarseningRemapper::unpack () const
{
  // If MPI does not use dev pointers, we need to deep copy from host to dev
  if (not MpiOnDev) {
    Kokkos::deep_copy (m_recv_buffer,m_mpi_recv_buffer);
  }

  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  const auto pid_lid_start = m_recv_pid_lids_start;
  const auto lids_pids = m_recv_lids_pids;
  const auto buf = m_recv_buffer;
  for (int ifield=0; ifield<m_num_fields; ++ifield) {
    const auto& f  = m_tgt_fields[ifield];
    const auto& fl = f.get_header().get_identifier().get_layout();
    const auto lt = get_layout_type(fl.tags());
    const auto f_pid_offsets = ekat::subview(m_recv_f_pid_offsets,ifield);

    f.deep_copy(0);
    switch (lt) {
      case LayoutType::Scalar2D:
      {
        auto v = f.get_view<Real*>();
        Kokkos::parallel_for(RangePolicy(0,m_total_num_recv_gids),
                             KOKKOS_LAMBDA(const int& i){
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          v(lid) += buf (offset + lidpos);
        });
                             
      } break;
      case LayoutType::Vector2D:
      {
        auto v = f.get_view<Real**>();
        const int ndims = fl.dim(1);
        auto policy = ESU::get_default_team_policy(m_total_num_recv_gids,ndims);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int i = team.league_rank();
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team,ndims),
                               [&](const int idim) {
            v(lid,idim) += buf (offset + lidpos*ndims + idim);
          });
        });
      } break;
      case LayoutType::Scalar3D:
      {
        auto v = f.get_view<Real**>();
        const int nlevs = fl.dims().back();
        auto policy = ESU::get_default_team_policy(m_total_num_recv_gids,nlevs);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int i = team.league_rank();
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team,nlevs),
                               [&](const int ilev) {
            v(lid,ilev) += buf (offset + lidpos*nlevs + ilev);
          });
        });
      } break;
      case LayoutType::Vector3D:
      {
        auto v = f.get_view<Real***>();
        const int ndims = fl.dim(1);
        const int nlevs = fl.dims().back();
        auto policy = ESU::get_default_team_policy(m_total_num_recv_gids,nlevs*ndims);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int i = team.league_rank();
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team,nlevs*ndims),
                               [&](const int idx) {
            const int idim = idx / nlevs;
            const int ilev = idx % nlevs;
            v(lid,idim,ilev) += buf (offset + lidpos*ndims*nlevs + idim*nlevs + ilev);
          });
        });
      } break;

      default:
        EKAT_ERROR_MSG ("Unexpected field rank in CoarseningRemapper::pack.\n"
            "  - MPI rank  : " + std::to_string(m_comm.rank()) + "\n"
            "  - field rank: " + std::to_string(fl.rank()) + "\n");
    }
  }
}


auto CoarseningRemapper::
get_my_triplets_gids (const std::string& map_file) const
  -> view_1d<gid_t>::HostMirror
{
  using namespace ShortFieldTagsNames;

  scorpio::register_file(map_file,scorpio::FileMode::Read);
  // Create a "fake" grid, with as many dofs as the number of triplets in the src file
  const int ngweights = scorpio::get_dimlen_c2f(map_file.c_str(),"n_s");
  const auto io_grid_linear = create_point_grid ("helper",ngweights,1,m_comm);
  const int nlweights = io_grid_linear->get_num_local_dofs();
  const auto& raw_comm = m_comm.mpi_comm();

  gid_t offset = nlweights;
  m_comm.scan(&offset,1,MPI_SUM);
  offset -= nlweights; // scan is inclusive, but we need exclusive

  // Read chunk of triplets col indices
  std::vector<gid_t> cols(nlweights);
  const std::string idx_decomp_tag = "int-nnz" + std::to_string(nlweights);
  scorpio::get_variable(map_file, "col", "col", {"n_s"}, "int", idx_decomp_tag);
  std::vector<scorpio::offset_t> dofs_offsets(nlweights);
  std::iota(dofs_offsets.begin(),dofs_offsets.end(),offset);
  scorpio::set_dof(map_file,"col",nlweights,dofs_offsets.data());
  scorpio::set_decomp(map_file);
  scorpio::grid_read_data_array(map_file,"col",-1,cols.data(),cols.size());
  scorpio::eam_pio_closefile(map_file);

  // Subtract 1 to cols indices
  for (auto& id : cols) {
    --id;
  }

  // Get owners of cols read in
  auto owners_lids = m_src_grid->get_owners_and_lids(cols);

  // Now we need to communicate to each rank which entries of the
  // (global) col array they need to read. Since the io_grid_linear
  // is partitioned linearly, we only need to add an offset to the
  // local id, in order to retrieve the global id
  // gid_t offset = 0;
  // MPI_Exscan (&nlweights,&offset,1,MPI_INT,MPI_SUM,raw_comm);

  std::map<int,std::vector<int>> pid_to_io_grid_gid;
  for (int i=0; i<nlweights; ++i) {
    const auto pid = owners_lids(i,0);
    pid_to_io_grid_gid[pid].push_back(i+offset);
  }

  MPI_Win counter_win, data_win;
  std::vector<int> num_entries(m_comm.size(),0);
  MPI_Win_create (num_entries.data(),m_comm.size()*sizeof(int),sizeof(int),
                  MPI_INFO_NULL,raw_comm,&counter_win);
  MPI_Win_fence(0,counter_win);

  // 1. Communicate number of gids to each pid
  for (const auto& it : pid_to_io_grid_gid) {
    const int n = it.second.size();
    const int pid = it.first;
    MPI_Put (&n,1,MPI_INT,pid,m_comm.rank(),1,MPI_INT,counter_win);
  }
  MPI_Win_fence(0,counter_win);
  MPI_Win_free(&counter_win);

  // 2. Compute offsets of each remote pid into our recv buffer
  int total_num_entries = num_entries[0];
  std::vector<int> pids_offset_in_my_buf (m_comm.size(),0);
  for (int pid=1; pid<m_comm.size(); ++pid) {
    pids_offset_in_my_buf[pid] = pids_offset_in_my_buf[pid-1] + num_entries[pid-1];
    total_num_entries += num_entries[pid];
  }

  // 3. Read from all other pids what MY offset is on their recv buffer
  std::vector<int> my_offset_on_pids (m_comm.size(),-1);
  MPI_Win read_win;
  MPI_Win_create (pids_offset_in_my_buf.data(),sizeof(int)*m_comm.size(),sizeof(int),
                  MPI_INFO_NULL,raw_comm,&read_win);
  MPI_Win_fence(0,read_win);
  for (const auto& it : pid_to_io_grid_gid) {
    const int pid = it.first;
    MPI_Get (&my_offset_on_pids[pid],1,MPI_INT,pid,m_comm.rank(),1,MPI_INT,read_win);
  }
  MPI_Win_fence(0,read_win);
  MPI_Win_free(&read_win);

  // 4. Communicate all gids to their owners
  gid_t* my_triplets_gids_ptr;
  MPI_Win_allocate(total_num_entries*sizeof(gid_t),sizeof(gid_t),
                   MPI_INFO_NULL,raw_comm,&my_triplets_gids_ptr,&data_win);
  MPI_Win_fence(0,data_win);
  for (const auto& it : pid_to_io_grid_gid) {
    const int n = it.second.size();
    const int pid = it.first;
    MPI_Put (it.second.data(),n,MPI_INT,pid,my_offset_on_pids[pid],n,MPI_INT,data_win);
  }
  MPI_Win_fence(0,data_win);

  // 5. Copy from window to a view, before freeing window's memory
  view_1d<gid_t>::HostMirror my_triplets_gids("",total_num_entries);
  std::memcpy(my_triplets_gids.data(),my_triplets_gids_ptr,total_num_entries*sizeof(int));

  // 6. Clean up
  MPI_Win_free(&data_win);

  return my_triplets_gids;
}

void CoarseningRemapper::create_ov_tgt_fields ()
{
  using FL = FieldLayout;
  m_ov_tgt_fields.reserve(m_num_fields);
  const int num_ov_cols = m_ov_tgt_grid->get_num_local_dofs();
  const auto ov_gn = m_ov_tgt_grid->name();
  // for (const auto& f : m_tgt_fields) {
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_src = m_src_fields[i];
    const auto& f_tgt = m_tgt_fields[i];
    const auto& fid = f_tgt.get_header().get_identifier();
    auto tags = fid.get_layout().tags();
    auto dims = fid.get_layout().dims();
    dims[0] = num_ov_cols;
    FieldIdentifier ov_fid (fid.name(),FL(tags,dims),fid.get_units(),ov_gn,fid.data_type());

    // Note: with C++17, emplace_back already returns a ref
    m_ov_tgt_fields.emplace_back(ov_fid);
    auto& ov_f = m_ov_tgt_fields.back();

    // Use same alloc props as src fields, to allow packing in local_mat_vec
    const auto pack_size = f_src.get_header().get_alloc_properties().get_largest_pack_size();
    ov_f.get_header().get_alloc_properties().request_allocation(pack_size);
    ov_f.allocate_view();
  }
}

void CoarseningRemapper::setup_mpi_data_structures ()
{
  using namespace ShortFieldTagsNames;

  const auto mpi_comm  = m_comm.mpi_comm();
  const auto mpi_real  = ekat::get_mpi_type<Real>();
  const auto mpi_gid_t = ekat::get_mpi_type<gid_t>();
  int pos;

  // Pre-compute the amount of data stored in each field on each dof
  std::vector<int> field_col_size (m_num_fields);
  int sum_fields_col_sizes = 0;
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f  = m_ov_tgt_fields[i];  // Doesn't matter if tgt or ov_tgt
    const auto& fl = f.get_header().get_identifier().get_layout();
    field_col_size[i] = fl.size() / fl.dim(0);
    sum_fields_col_sizes += field_col_size[i];
  }

  // --------------------------------------------------------- //
  //                   Setup SEND structures                   //
  // --------------------------------------------------------- //

  // 1. Retrieve pid (and associated lid) of all ov_tgt gids
  //    on the tgt grid
  const auto ov_gids = m_ov_tgt_grid->get_dofs_gids_host();
  auto pids_lids = m_tgt_grid->get_owners_and_lids (ov_gids);

  // 2. Group dofs to send by remote pid
  const int num_ov_gids = ov_gids.size();
  std::map<int,std::vector<int>> send_lids;
  std::map<int,std::vector<gid_t>> send_gids;
  for (int i=0; i<num_ov_gids; ++i) {
    const int pid = pids_lids(i,0);
    send_lids[pid].push_back(i);
    send_gids[pid].push_back(ov_gids(i));
  }
  const int num_send_pids = send_lids.size();
  m_send_lids_pids = view_2d<int>("",num_ov_gids,2);
  m_send_pid_lids_start = view_1d<int>("",m_comm.size());
  auto send_lids_pids_h = Kokkos::create_mirror_view(m_send_lids_pids);
  auto send_pid_lids_start_h = Kokkos::create_mirror_view(m_send_pid_lids_start);
  pos = 0;
  for (int pid=0; pid<m_comm.size(); ++pid) {
    send_pid_lids_start_h(pid) = pos;
    for (auto lid : send_lids[pid]) {
      send_lids_pids_h(pos,0) = lid;
      send_lids_pids_h(pos++,1) = pid;
    }
  }
  Kokkos::deep_copy(m_send_lids_pids,send_lids_pids_h);
  Kokkos::deep_copy(m_send_pid_lids_start,send_pid_lids_start_h);

  // 3. Compute offsets in send buffer for each pid/field pair
  m_send_f_pid_offsets = view_2d<int>("",m_num_fields,m_comm.size());
  m_send_pid_offsets = view_1d<int>("",m_comm.size());
  auto send_f_pid_offsets_h = Kokkos::create_mirror_view(m_send_f_pid_offsets);
  auto send_pid_offsets_h = Kokkos::create_mirror_view(m_send_pid_offsets);
  pos = 0;
  for (int pid=0; pid<m_comm.size(); ++pid) {
    send_pid_offsets_h(pid) = pos;
    for (int i=0; i<m_num_fields; ++i) {
      send_f_pid_offsets_h(i,pid) = pos;
      pos += field_col_size[i]*send_lids[pid].size();
    }
  }
  EKAT_REQUIRE_MSG (pos==num_ov_gids*sum_fields_col_sizes,
      "Error! Something went wrong in CoarseningRemapper::setup_mpi_structures.\n");
  Kokkos::deep_copy (m_send_pid_offsets,send_pid_offsets_h);

  // 4. Allocate send buffers
  m_send_buffer = view_1d<Real>("",sum_fields_col_sizes*num_ov_gids);
  m_mpi_send_buffer = Kokkos::create_mirror_view(decltype(m_mpi_send_buffer)::execution_space(),m_send_buffer);

  // 5. Setup send requests
  m_send_req.reserve(num_send_pids);
  for (const auto& it : send_lids) {
    const int n = it.second.size()*sum_fields_col_sizes;
    if (n==0) {
      continue;
    }

    const int pid = it.first;
    const auto send_ptr = m_mpi_send_buffer.data() + send_pid_offsets_h(pid);

    m_send_req.emplace_back();
    auto& req = m_send_req.back();
    MPI_Send_init (send_ptr, n, mpi_real, pid,
                   0, mpi_comm, &req);
  }
  m_send_req_ptr = m_send_req.data();

  // --------------------------------------------------------- //
  //                   Setup RECV structures                   //
  // --------------------------------------------------------- //

  MPI_Win write_win, read_win;

  // 1. Let each process know how many lids we need to send them
  std::vector<int> num_recv_gids (m_comm.size(),0);
  MPI_Win_create (num_recv_gids.data(),sizeof(int)*m_comm.size(),sizeof(int),
                  MPI_INFO_NULL,mpi_comm,&write_win);
  MPI_Win_fence(0,write_win);
  for (const auto& it : send_gids) {
    const int pid = it.first;
    const int n   = it.second.size();
    if (n>0) {
      MPI_Put (&n,1,MPI_INT,pid,m_comm.rank(),1,MPI_INT,write_win);
    }
  }
  MPI_Win_fence(0,write_win);
  MPI_Win_free(&write_win);
  
  // 2. Allocate space for other pids to communicate the gids.
  //    Then, communicate gids to all other procs.
  m_total_num_recv_gids = 0;
  m_recv_pid_lids_start = view_1d<int>("",m_comm.size()+1);
  auto recv_pid_lids_start_h = Kokkos::create_mirror_view(m_recv_pid_lids_start);
  int num_recv_pids = 0;
  for (int pid=0; pid<m_comm.size(); ++pid) {
    recv_pid_lids_start_h[pid] = m_total_num_recv_gids;
    m_total_num_recv_gids += num_recv_gids[pid];
    if (num_recv_gids[pid]>0) {
      ++num_recv_pids;
    }
  }
  recv_pid_lids_start_h[m_comm.size()] = m_total_num_recv_gids;
  Kokkos::deep_copy(m_recv_pid_lids_start,recv_pid_lids_start_h);
  std::vector<gid_t> recv_gids (m_total_num_recv_gids,-1);

  // 3. Read from all other pids what MY offset is on their recv buffer
  MPI_Win_create (recv_pid_lids_start_h.data(),sizeof(int)*m_comm.size(),sizeof(int),
                  MPI_INFO_NULL,mpi_comm,&read_win);
  std::vector<int> my_offset_on_pids (m_comm.size());
  MPI_Win_fence(0,read_win);
  for (int pid=0; pid<m_comm.size(); ++pid) {
    MPI_Get (&my_offset_on_pids[pid],1,MPI_INT,pid,m_comm.rank(),1,MPI_INT,read_win);
  }
  MPI_Win_fence(0,read_win);
  MPI_Win_free(&read_win);

  // 4. Communicate gids we are sending to all pids
  MPI_Win_create (recv_gids.data(),sizeof(gid_t)*m_total_num_recv_gids,sizeof(gid_t),
                  MPI_INFO_NULL,mpi_comm,&write_win);
  MPI_Win_fence(0,write_win);
  for (const auto& it : send_gids) {
    const int n = it.second.size();
    if (n==0) {
      continue;
    }
    const int pid = it.first;
    MPI_Put (it.second.data(),n,mpi_gid_t,pid,my_offset_on_pids[pid],n,mpi_gid_t,write_win);
  }
  MPI_Win_fence(0,write_win);
  MPI_Win_free(&write_win);

  // 5. Now we know how many gids we recv from each pid, and what they are
  //    Splice gids from all pids, and convert to lid
  m_recv_lids_pids = view_2d<int>("",m_total_num_recv_gids,2);
  auto recv_lids_pids_h = Kokkos::create_mirror_view(m_recv_lids_pids);
  pos = 0;
  for (int pid=0; pid<m_comm.size(); ++pid) {
    const int start = recv_pid_lids_start_h[pid];
    const int end   = recv_pid_lids_start_h[pid+1];
    for (int i=start; i<end; ++i) {
      const auto gid = recv_gids[i];
      const auto lid = gid2lid(gid,m_tgt_grid);
      EKAT_REQUIRE_MSG (lid>=0,
          "Error! Something went wrong while converting gid->lid on tgt grid.\n"
          "  - rank: " + std::to_string(m_comm.rank()) + "\n");
      recv_lids_pids_h(pos,0)   = lid;
      recv_lids_pids_h(pos++,1) = pid;
    }
  }

  // 6. Compute offsets in recv buffer for each pid/field pair
  m_recv_f_pid_offsets = view_2d<int>("",m_num_fields,m_comm.size());
  m_recv_pid_offsets = view_1d<int>("",m_comm.size());
  auto recv_f_pid_offsets_h = Kokkos::create_mirror_view(m_recv_f_pid_offsets);
  auto recv_pid_offsets_h = Kokkos::create_mirror_view(m_recv_pid_offsets);
  pos = 0;
  for (int pid=0; pid<m_comm.size(); ++pid) {
    recv_pid_offsets_h(pid) = pos;
    for (int i=0; i<m_num_fields; ++i) {
      recv_f_pid_offsets_h(i,pid) = pos;
      pos += field_col_size[i]*num_recv_gids[pid];
    }
  }
  EKAT_REQUIRE_MSG (pos==m_total_num_recv_gids*sum_fields_col_sizes,
      "Error! Something went wrong in CoarseningRemapper::setup_mpi_structures.\n");
  Kokkos::deep_copy (m_recv_f_pid_offsets,recv_f_pid_offsets_h);
  Kokkos::deep_copy (m_recv_pid_offsets,recv_pid_offsets_h);

  // 7. Allocate recv buffers
  m_recv_buffer = view_1d<Real>("",sum_fields_col_sizes*m_total_num_recv_gids);
  m_mpi_recv_buffer = Kokkos::create_mirror_view(decltype(m_mpi_recv_buffer)::execution_space(),m_recv_buffer);

  // 8. Setup recv requests
  m_recv_req.reserve(num_recv_pids);
  for (int pid=0; pid<m_comm.size(); ++pid) {
    const int n = num_recv_gids[pid]*sum_fields_col_sizes;
    if (n==0) {
      continue;
    }

    const auto recv_ptr = m_mpi_recv_buffer.data() + recv_pid_offsets_h(pid);

    m_recv_req.emplace_back();
    auto& req = m_recv_req.back();
    MPI_Recv_init (recv_ptr, n, mpi_real, pid,
                   0, mpi_comm, &req);
  }
  m_recv_req_ptr = m_recv_req.data();
}

} // namespace scream
