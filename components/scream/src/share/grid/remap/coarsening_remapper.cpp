#include "coarsening_remapper.hpp"

#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/ekat_pack_utils.hpp>

namespace scream
{

CoarseningRemapper::
CoarseningRemapper (const grid_ptr_type& src_grid,
                    // const grid_ptr_type& tgt_grid,
                    const std::string& map_file)
 : AbstractRemapper(src_grid,create_tgt_grid(map_file,src_grid))
 , m_comm (src_grid->get_comm())
{
  using view_1d_host = AtmosphereInput::view_1d_host;
  using namespace ShortFieldTagsNames;

  scorpio::register_file(map_file,scorpio::FileMode::Read);

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

  auto io_grid = std::make_shared<PointGrid>("",gids_h.size(),2,m_comm);
  io_grid->set_dofs(gids_d);
  const int nlweights = io_grid->get_num_local_dofs();

  // Read my row, col, S entries
  m_weights = view_1d<Real>("",nlweights);
  m_row_col_lids = view_2d<int>("",nlweights,2);
  view_1d<int> row("",nlweights);
  view_1d<int> col("",nlweights);

  // TODO: make AtmosphereInput accept view_1d_host<T>,
  //       so we can read straight into views of ints
  view_1d_host row_h_real("",nlweights), col_h_real("",nlweights);
  auto S_h = Kokkos::create_mirror_view(m_weights);
  ekat::ParameterList reader_params;
  reader_params.set("Filename",map_file);
  std::vector<std::string> fields = {"row","col","S"};
  std::map<std::string,view_1d_host> views = {
    {"row", row_h_real}, {"col", col_h_real}, {"S", S_h}
  };
  std::map<std::string,FieldLayout>  layouts = {
    {"row", FieldLayout({COL},{nlweights})},
    {"col", FieldLayout({COL},{nlweights})},
    {"S", FieldLayout({COL},{nlweights})}
  };

  AtmosphereInput reader(reader_params,io_grid,views,layouts);
  reader.read_variables();

  scorpio::eam_pio_closefile(map_file);

  // Copy triplets to device views
  Kokkos::deep_copy(m_weights,S_h);
  Kokkos::deep_copy(row,row_h_real);
  Kokkos::deep_copy(col,col_h_real);

  auto row_col_lids = m_row_col_lids;
  Kokkos::parallel_for(typename KT::RangePolicy(0,nlweights),
                       KOKKOS_LAMBDA(const int& i){
    row_col_lids(i,0) = row(i);
    row_col_lids(i,1) = col(i);
  });

  // Create an "overlapped" tgt grid, that is, a grid where each rank
  // owns all tgt rows that are affected by at least one of the cols
  // in its src_grid
  auto ov_tgt_grid = std::make_shared<PointGrid>("ov_tgt_grid",nlweights,2,m_comm);
  ov_tgt_grid->set_dofs(row);
  m_ov_tgt_grid = ov_tgt_grid;
}

CoarseningRemapper::
~CoarseningRemapper ()
{
  // We need to free MPI windows and datatypes
  for (auto& w : m_tgt_fields_win) {
    MPI_Win_free(&w);
  }

  for (int i=0; i<m_num_fields; ++i) {
    for (auto& it : m_local_dtypes[i]) {
      MPI_Type_free (&it.second);
    }
    for (auto& it : m_remote_dtypes[i]) {
      MPI_Type_free (&it.second);
    }
  }

  MPI_Group_free (&m_rma_post_group);
  MPI_Group_free (&m_rma_start_group);
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

  const int ncols_src = scorpio::get_dimlen_c2f(map_file.c_str(),"n_b");
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
  // // Fire the recv requests right away, so that if some other ranks
  // // is done packing before us, we can start receiving their data
  // if (not m_recv_req.empty()) {
  //   int ierr = MPI_Startall(m_recv_req.size(),m_recv_req_ptr);
  //   EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
  //       "Error! Something whent wrong while starting persistent recv requests.\n"
  //       "  - recv rank: " + std::to_string(m_comm.rank()) + "\n");
  // }

  // Expose our windows for writing, and tell other pids we could start
  // writing any moment. We do this right away (rather than after
  // completing local_mat_vec) so that
  //  a) other ranks can start writing if they finish local_mat_vec before us
  //  b) even if they finish both local_mat_vec AND rma ops before we finish
  //     local_mat_vec, their MPI_Win_wait call will wait for us.
  // Also, zero out the tgt field *before* the post.
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_tgt = m_ov_tgt_fields[i];

    // Zero out to prepare for accumulation phase
    f_tgt.deep_copy(0);

    MPI_Win_start (m_rma_start_group,0,m_tgt_fields_win[i]);
    MPI_Win_post  (m_rma_post_group, 0,m_tgt_fields_win[i]);
  }

  // Loop over each field
  constexpr auto can_pack = SCREAM_PACK_SIZE>1;
  constexpr HostOrDevice HD = MpiOnDev ? Device : Host;
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
  
    // Now we can do our RMA ops, to transfer data from overlapped
    // target fields to target fields
    const auto& win = m_tgt_fields_win[i];
    const auto ov_tgt_data = f_ov_tgt.get_internal_view_data<const Real,HD>();
    for (const auto pid : m_remote_pids) {
      const auto& loc_dtype = m_local_dtypes[i].at(pid);
      const auto& rem_dtype = m_remote_dtypes[i].at(pid);
      MPI_Accumulate (ov_tgt_data,1,loc_dtype,pid,
                      0,1,rem_dtype,MPI_SUM,win);
    }

    // Let other pids know we're done writing on this window
    MPI_Win_complete(win);
  }

  // Wait for others to be done writing
  // NOTE: there's no point in looping on windows and calling MPI_Win_test
  //       first to not get stuck on a single wait, since we have to wait
  //       on all windows, so we'll have to wait on the slowest anyways.
  for (int i=0; i<m_num_fields; ++i) {
    MPI_Win_wait (m_tgt_fields_win[i]);
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
  const int nnz = m_weights.size();
  auto row_col_lids = m_row_col_lids;
  auto weights = m_weights;
  switch (rank) {
    case 1:
    {
      auto x_view = x.get_view<const Real*>();
      auto y_view = y.get_view<      Real*>();
      Kokkos::parallel_for(RangePolicy(0,nnz),
                           KOKKOS_LAMBDA(const int& i) {
        const auto row = row_col_lids(i,0);
        const auto col = row_col_lids(i,1);
        const auto w = weights(i);
        y_view(row) += w*x_view(col);
      });
      break;
    }
    case 2:
    {
      auto x_view = x.get_view<const Pack**>();
      auto y_view = y.get_view<      Pack**>();
      const int dim1 = PackInfo::num_packs(src_layout.dim(1));
      auto policy = ESU::get_default_team_policy(nnz,dim1);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto i   = team.league_rank();
        const auto row = row_col_lids(i,0);
        const auto col = row_col_lids(i,1);
        const auto w = weights(i);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team,dim1),
                            [&](const int j){
          y_view(row,j) += w*x_view(col,j);
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
      auto policy = ESU::get_default_team_policy(nnz,dim1*dim2);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto i   = team.league_rank();
        const auto row = row_col_lids(i,0);
        const auto col = row_col_lids(i,1);
        const auto w = weights(i);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team,dim1),
                            [&](const int idx){
          const int j = idx / dim2;
          const int k = idx % dim2;
          y_view(row,j,k) += w*x_view(col,j,k);
        });
      });
      break;
    }
  }
}

auto CoarseningRemapper::
get_my_triplets_gids (const std::string& map_file) const
  -> view_1d<gid_t>::HostMirror
{
  using namespace ShortFieldTagsNames;
  using view_1d_host = AtmosphereInput::view_1d_host;

  // Create a "fake" grid, with as many dofs as the number of triplets in the src file
  const int ngweights = scorpio::get_dimlen_c2f(map_file.c_str(),"n_s");
  const auto io_grid_linear = create_point_grid ("helper",ngweights,1,m_comm);
  const int nlweights = io_grid_linear->get_num_local_dofs();
  const auto& raw_comm = m_comm.mpi_comm();

  ekat::ParameterList reader_params;
  reader_params.set("Filename",map_file);
  std::vector<std::string> fields = {"col"};
  std::map<std::string,view_1d_host> views = {
    {"col", view_1d_host("",nlweights)}
  };
  std::map<std::string,FieldLayout>  layouts = {
    {"col", FieldLayout({COL},{nlweights})}
  };

  AtmosphereInput reader(reader_params,io_grid_linear,views,layouts);
  reader.read_variables();

  // TODO: make atm input accept view_1d_host<T>, so we can read directly
  //       into a view of ints
  AbstractGrid::dofs_list_h_type cols ("",nlweights);
  for (int i=0; i<nlweights; ++i) {
    cols(i) = views.at("cols")(i);
  }

  // Get owners of cols read in
  auto owners_lids = m_src_grid->get_owners_and_lids(cols);

  // Now we need to communicate to each rank which entries of the
  // (global) col array they need to read. Since the io_grid_linear
  // is partitioned linearly, we only need to add an offset to the
  // local id, in order to retrieve the global id
  gid_t offset = 0;
  MPI_Exscan (&nlweights,&offset,1,MPI_INT,MPI_SUM,raw_comm);

  std::map<int,std::vector<int>> pid_to_io_grid_gid;
  for (int i=0; i<nlweights; ++i) {
    const auto pid = owners_lids(i,0);
    pid_to_io_grid_gid[pid].push_back(i+offset);
  }

  MPI_Win counter_win, data_win;
  int num_entries = 0;
  MPI_Win_create (&num_entries,sizeof(int),sizeof(int),
                  MPI_INFO_NULL,raw_comm,&counter_win);
  MPI_Win_fence(0,counter_win);

  // First accumulate number of gids into the owner count
  for (const auto& it : pid_to_io_grid_gid) {
    const int n = it.second.size();
    const int pid = it.first;
    MPI_Accumulate (&n,1,MPI_INT,pid,0,1,MPI_INT,MPI_SUM,counter_win);
  }
  MPI_Win_fence(0,counter_win);
  MPI_Win_free(&counter_win);

  // Second, communicate the gids to their owners
  int write_offset = 0;
  MPI_Win_create(&write_offset,sizeof(int),sizeof(int),
                 MPI_INFO_NULL,raw_comm,&counter_win);
  gid_t* my_triplets_gids_ptr;
  MPI_Win_allocate(num_entries*sizeof(gid_t),sizeof(gid_t),
                   MPI_INFO_NULL,raw_comm,&my_triplets_gids_ptr,&data_win);
  MPI_Win_fence(0,data_win);
  MPI_Win_fence(0,counter_win);
  for (const auto& it : pid_to_io_grid_gid) {
    const int n = it.second.size();
    const int pid = it.first;
    int write_start = -1;
    MPI_Get_accumulate (&n,1,MPI_INT,&write_start,1,MPI_INT,pid,
                        0,1,MPI_INT,MPI_SUM,counter_win);

    MPI_Accumulate (&n,1,MPI_INT,pid,write_start,1,MPI_INT,MPI_SUM,data_win);
  }
  MPI_Win_fence(0,data_win);
  MPI_Win_fence(0,counter_win);

  // Copy from window to a view, before freeing window's memory
  view_1d<gid_t>::HostMirror my_triplets_gids("",num_entries);
  std::memcpy(my_triplets_gids.data(),my_triplets_gids_ptr,num_entries*sizeof(int));

  // Clean up
  MPI_Win_free(&counter_win);
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
  const auto raw_comm = m_comm.mpi_comm();

  // 1. Retrieve pid (and associated lid) of all ov_tgt gids
  //    on the tgt grid
  const auto ov_gids = m_ov_tgt_grid->get_dofs_gids_host();
  auto pids_lids = m_tgt_grid->get_owners_and_lids (ov_gids);

  // 2. Group dofs by remote pid (to issue a single RMA call per pid later)
  const int num_ov_dofs = ov_gids.size();
  std::map<int,std::vector<int>> local_offsets, remote_offsets;
  for (int i=0; i<num_ov_dofs; ++i) {
    const int pid = pids_lids(i,0);
    m_remote_pids.push_back(pid);
    const int remote_lid = pids_lids(i,1);
    local_offsets[pid].push_back(i);
    remote_offsets[pid].push_back(remote_lid);
  }

  // 3. Create MPI_Datatype for both local (ov_tgt_field)
  //    and remote (tgt_field)
  m_local_dtypes.resize(m_num_fields);
  m_remote_dtypes.resize(m_num_fields);
  // m_col_dtypes.resize(m_num_fields);
  const auto mpi_real = ekat::get_mpi_type<Real>();
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f = m_tgt_fields[i];
    const auto& ov_f = m_tgt_fields[i];
    const auto& fl  = f.get_header().get_identifier().get_layout();
    const auto& fap = f.get_header().get_alloc_properties();
    const auto& ov_fap = ov_f.get_header().get_alloc_properties();

    int col_size = fl.rank() / fl.dim(0);

    int col_disp = col_size;
    int ov_col_disp = col_size;
    if (fap.get_padding()>0) {
      col_disp = fap.get_alloc_size() / sizeof(Real) / fl.dim(0);
    }
    if (ov_fap.get_padding()>0) {
      ov_col_disp = ov_fap.get_alloc_size() / sizeof(Real) / fl.dim(0);
    }
    // MPI_Type_contiguous(col_size,get_mpi_type<Real>(),&m_col_dtypes[i]);
    // MPI_Type_commit(&m_col_dtypes[i]);

    for (const auto& it : local_offsets) {
      const int pid = it.first;
      auto loc_offsets = it.second;
      auto rem_offsets = remote_offsets.at(pid);

      // Multiply offsets by column displacement
      for (auto& o : loc_offsets) {
        o *= ov_col_disp;
      }
      for (auto& o : rem_offsets) {
        o *= col_disp;
      }
      const std::vector<int> blocks_sizes (loc_offsets.size(),col_size);

      auto& loc_dtype = m_local_dtypes[i][pid];
      auto& rem_dtype = m_remote_dtypes[i][pid];

      // MPI_Type_indexed (loc_offset.size(),ones.data(),loc_offsets.data(),m_col_dtypes[i],&loc_dtype);
      // MPI_Type_indexed (rem_offset.size(),ones.data(),rem_offsets.data(),m_col_dtypes[i],&rem_dtype);
      MPI_Type_indexed (loc_offsets.size(),blocks_sizes.data(),loc_offsets.data(),mpi_real,&loc_dtype);
      MPI_Type_indexed (rem_offsets.size(),blocks_sizes.data(),rem_offsets.data(),mpi_real,&rem_dtype);
      MPI_Type_commit(&loc_dtype);
      MPI_Type_commit(&rem_dtype);
    }
  }

  // 4. Create windows
  m_tgt_fields_win.resize(m_num_fields);
  constexpr HostOrDevice HD = MpiOnDev ? Device : Host;
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f = m_tgt_fields[i];
    const auto& fap = f.get_header().get_alloc_properties();
    auto win_size = fap.get_alloc_size();
    auto data = f.get_internal_view_data<Real,HD>();
    MPI_Win_create (data,win_size,sizeof(Real),
                    MPI_INFO_NULL,raw_comm,&m_tgt_fields_win[i]);
  }

  // 5. Create groups for windows post/start calls
  MPI_Comm_group(m_comm.mpi_comm(),&m_rma_post_group);
  MPI_Group_incl(m_rma_post_group,m_remote_pids.size(),m_remote_pids.data(),&m_rma_start_group);
}

} // namespace scream
