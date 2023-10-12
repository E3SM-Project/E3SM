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
 : AbstractRemapper()
 , m_comm (tgt_grid->get_comm())
{
  using namespace ShortFieldTagsNames;

  // Sanity checks
  EKAT_REQUIRE_MSG (tgt_grid->type()==GridType::Point,
      "Error! RefiningRemapperRMA only works on PointGrid grids.\n"
      "  - tgt grid name: " + tgt_grid->name() + "\n"
      "  - tgt_grid_type: " + e2str(tgt_grid->type()) + "\n");
  EKAT_REQUIRE_MSG (tgt_grid->is_unique(),
      "Error! RefiningRemapperRMA requires a unique target grid.\n");

  // This is a refining remapper. We only go in one direction
  m_bwd_allowed = false;

  // Load (i,j,w) triplets from map file, for all i that are
  // owned on the tgt_grid
  auto my_triplets = get_my_triplets (map_file,tgt_grid);

  // Create an overlapped src map, consisting of all the col gids
  // in the triplets. This is overlapped, since for each gid there
  // may be 2+ ranks owning it.
  std::map<gid_type,int> ov_src_gid2lid;
  for (const auto& t : my_triplets) {
    ov_src_gid2lid.emplace(t.col,ov_src_gid2lid.size());
  }
  int num_ov_src_gids = ov_src_gid2lid.size();
  auto ov_src_grid = std::make_shared<PointGrid>("ov_src_grid",num_ov_src_gids,0,m_comm);
  auto ov_src_gids_h = ov_src_grid->get_dofs_gids().get_view<gid_type*,Host>();
  for (const auto& it : ov_src_gid2lid) {
    ov_src_gids_h[it.second] = it.first;
  }
  ov_src_grid->get_dofs_gids().sync_to_dev();
  m_ov_src_grid = ov_src_grid;

  // Create a unique version of m_ov_src_grid
  auto src_grid_gids = m_ov_src_grid->get_unique_gids();
  const int ngids = src_grid_gids.size();
  const int nlevs  = tgt_grid->get_num_vertical_levels();
  auto src_grid = std::make_shared<PointGrid>("src_grid",ngids,nlevs,m_comm);
  auto src_grid_gids_h = src_grid->get_dofs_gids().get_view<gid_type*,Host>();
  std::memcpy(src_grid_gids_h.data(),src_grid_gids.data(),ngids*sizeof(gid_type));
  src_grid->get_dofs_gids().sync_to_dev();

  // Finally able to set the src and tgt grids
  this->set_grids(src_grid,tgt_grid);

  // 5. Create CRS views and host mirrors
  const int num_my_triplets = my_triplets.size();
  m_row_offsets = view_1d<int>("",tgt_grid->get_num_local_dofs()+1);
  m_col_lids = view_1d<int>("",num_my_triplets);
  m_weights = view_1d<Real>("",num_my_triplets);

  auto row_offsets_h = Kokkos::create_mirror_view(m_row_offsets);
  auto col_lids_h = Kokkos::create_mirror_view(m_col_lids);
  auto weights_h = Kokkos::create_mirror_view(m_weights);

  std::vector<int> num_entries_per_row(tgt_grid->get_num_local_dofs(),0);
  auto gid2lid_row = tgt_grid->get_gid2lid_map();
  for (int i=0; i<num_my_triplets; ++i) {
    const auto& t = my_triplets[i];
    ++num_entries_per_row[gid2lid_row.at(t.row)];
    col_lids_h[i] = ov_src_gid2lid.at(t.col);
    weights_h[i] = t.w;
  }
  row_offsets_h(0) = 0;
  for (int i=0; i<tgt_grid->get_num_local_dofs(); ++i) {
    row_offsets_h(i+1) = row_offsets_h(i) + num_entries_per_row[i];
  }
  Kokkos::deep_copy(m_row_offsets,row_offsets_h);
  Kokkos::deep_copy(m_col_lids,   col_lids_h);
  Kokkos::deep_copy(m_weights,    weights_h);
}

RefiningRemapperRMA::
~RefiningRemapperRMA ()
{
  clean_up();
}

FieldLayout RefiningRemapperRMA::
create_src_layout (const FieldLayout& tgt_layout) const
{
  using namespace ShortFieldTagsNames;
  const auto lt = get_layout_type(tgt_layout.tags());
  auto src = FieldLayout::invalid();
  const bool midpoints = tgt_layout.has_tag(LEV);
  const int vec_dim = tgt_layout.is_vector_layout() ? tgt_layout.dim(CMP) : -1;
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
      EKAT_ERROR_MSG ("Layout not supported by RefiningRemapperRMA: " + e2str(lt) + "\n");
  }
  return src;
}

FieldLayout RefiningRemapperRMA::
create_tgt_layout (const FieldLayout& src_layout) const
{
  using namespace ShortFieldTagsNames;
  const auto lt = get_layout_type(src_layout.tags());
  auto tgt = FieldLayout::invalid();
  const bool midpoints = src_layout.has_tag(LEV);
  const int vec_dim = src_layout.is_vector_layout() ? src_layout.dim(CMP) : -1;
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
      EKAT_ERROR_MSG ("Layout not supported by RefiningRemapperRMA: " + e2str(lt) + "\n");
  }
  return tgt;
}

void RefiningRemapperRMA::
do_register_field (const identifier_type& src, const identifier_type& tgt)
{
  constexpr auto COL = ShortFieldTagsNames::COL;
  EKAT_REQUIRE_MSG (src.get_layout().has_tag(COL),
      "Error! Cannot register a field without COL tag in RefiningRemapperRMA.\n"
      "  - field name: " + src.name() + "\n"
      "  - field layout: " + to_string(src.get_layout()) + "\n");
  m_src_fields.push_back(field_type(src));
  m_tgt_fields.push_back(field_type(tgt));
}

void RefiningRemapperRMA::
do_bind_field (const int ifield, const field_type& src, const field_type& tgt)
{
  EKAT_REQUIRE_MSG (src.data_type()==DataType::RealType,
      "Error! RefiningRemapperRMA only allows fields with RealType data.\n"
      "  - src field name: " + src.name() + "\n"
      "  - src field type: " + e2str(src.data_type()) + "\n");
  EKAT_REQUIRE_MSG (tgt.data_type()==DataType::RealType,
      "Error! RefiningRemapperRMA only allows fields with RealType data.\n"
      "  - tgt field name: " + tgt.name() + "\n"
      "  - tgt field type: " + e2str(tgt.data_type()) + "\n");

  m_src_fields[ifield] = src;
  m_tgt_fields[ifield] = tgt;

  // If this was the last field to be bound, we can setup the MPI schedule
  if (this->m_state==RepoState::Closed &&
      (this->m_num_bound_fields+1)==this->m_num_registered_fields) {
    create_ov_src_fields ();
    setup_mpi_data_structures ();
  }
}

void RefiningRemapperRMA::do_registration_ends ()
{
  if (this->m_num_bound_fields==this->m_num_registered_fields) {
    create_ov_src_fields ();
    setup_mpi_data_structures ();
  }
}

void RefiningRemapperRMA::do_remap_fwd ()
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
    auto ov_data = m_ov_src_fields[i].get_internal_view_data<Real,MpiDev>();
    for (int icol=0; icol<m_ov_src_grid->get_num_local_dofs(); ++icol) {
      const int pid = m_remote_pids[icol];
      const int lid = m_remote_lids[icol];
      check_mpi_call(MPI_Get(ov_data+icol*col_size,col_size,dt,pid,
                             lid*col_stride+col_offset,col_size,dt,win),
                     "MPI_Get for field: " + m_ov_src_fields[i].name());
    }
  }

  // Close access RMA epoch on each field (exposure is still open)
  for (int i=0; i<m_num_fields; ++i) {
    check_mpi_call(MPI_Win_complete(m_mpi_win[i]),
                   "MPI_Win_complete for field: " + m_ov_src_fields[i].name());
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
    const auto& f_ov_src    = m_ov_src_fields[i];

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

template<int PackSize>
void RefiningRemapperRMA::
local_mat_vec (const Field& x, const Field& y) const
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  using Pack        = ekat::Pack<Real,PackSize>;
  using PackInfo    = ekat::PackInfo<PackSize>;

  const auto& src_layout = x.get_header().get_identifier().get_layout();
  const int rank = src_layout.rank();
  const int nrows = m_tgt_grid->get_num_local_dofs();
  auto row_offsets = m_row_offsets;
  auto col_lids = m_col_lids;
  auto weights = m_weights;
  switch (rank) {
    // Note: in each case, handle 1st contribution to each row separately,
    //       using = instead of +=. This allows to avoid doing an extra
    //       loop to zero out y before the mat-vec.
    case 1:
    {
      // Unlike get_view, get_strided_view returns a LayoutStride view,
      // therefore allowing the 1d field to be a subfield of a 2d field
      // along the 2nd dimension.
      auto x_view = x.get_strided_view<const Real*>();
      auto y_view = y.get_strided_view<      Real*>();
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
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1),
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
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2),
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
    default:
    {
      EKAT_ERROR_MSG("Error::refining_remapper::local_mat_vec doesn't support fields of rank 4 or greater");
    }
  }
}

auto RefiningRemapperRMA::
get_my_triplets (const std::string& map_file,
                 const grid_ptr_type& tgt_grid)
 -> std::vector<Triplet>
{
  using namespace ShortFieldTagsNames;
  constexpr int one = 1;

  // 1. Load the map file chunking it evenly across all ranks
  scorpio::register_file(map_file,scorpio::FileMode::Read);

  // 1.1 Create a "helper" grid, with as many dofs as the number
  //     of triplets in the map file, and divided linearly across ranks
  const int ngweights = scorpio::get_dimlen(map_file,"n_s");
  const auto io_grid_linear = create_point_grid ("helper",ngweights,1,m_comm);
  const int nlweights = io_grid_linear->get_num_local_dofs();

  gid_type offset = nlweights;
  m_comm.scan(&offset,1,MPI_SUM);
  offset -= nlweights; // scan is inclusive, but we need exclusive

  // Create a unique decomp tag, which ensures all refining remappers have
  // their own decomposition
  const std::string int_decomp_tag  = "RR::gmtg,int,grid-idx=" + std::to_string(io_grid_linear->get_unique_grid_id());
  const std::string real_decomp_tag = "RR::gmtg,real,grid-idx=" + std::to_string(io_grid_linear->get_unique_grid_id());

  // 1.2 Read a chunk of triplets col indices
  std::vector<gid_type> cols(nlweights);
  std::vector<gid_type> rows(nlweights);
  std::vector<Real>  S(nlweights);

  scorpio::register_variable(map_file, "col", "col", {"n_s"}, "int",  int_decomp_tag);
  scorpio::register_variable(map_file, "row", "row", {"n_s"}, "int",  int_decomp_tag);
  scorpio::register_variable(map_file, "S",   "S",   {"n_s"}, "real", real_decomp_tag);

  std::vector<scorpio::offset_t> dofs_offsets(nlweights);
  std::iota(dofs_offsets.begin(),dofs_offsets.end(),offset);
  scorpio::set_dof(map_file,"col",nlweights,dofs_offsets.data());
  scorpio::set_dof(map_file,"row",nlweights,dofs_offsets.data());
  scorpio::set_dof(map_file,"S"  ,nlweights,dofs_offsets.data());
  scorpio::set_decomp(map_file);

  scorpio::grid_read_data_array(map_file,"col",-1,cols.data(),cols.size());
  scorpio::grid_read_data_array(map_file,"row",-1,rows.data(),rows.size());
  scorpio::grid_read_data_array(map_file,"S"  ,-1,S.data(),S.size());

  scorpio::eam_pio_closefile(map_file);

  // 1.3 Dofs in tgt grid are likely 0-based, while row ids in map file
  // are likely 1-based. To match dofs, we need to offset the row
  // ids we read in.
  int map_file_min_row = std::numeric_limits<int>::max();
  for (int id=0; id<nlweights; id++) {
    map_file_min_row = std::min(rows[id],map_file_min_row);
  }
  int global_map_file_min_row;
  m_comm.all_reduce(&map_file_min_row,&global_map_file_min_row,1,MPI_MIN);

  gid_type row_offset = global_map_file_min_row - tgt_grid->get_global_min_dof_gid();
  for (auto& id : rows) {
    id -= row_offset;
  }

  // 2. Get the owners of the row gids we read in, according to the tgt grid
  std::vector<int> pids, lids;
  tgt_grid->get_remote_pids_and_lids(rows,pids,lids);

  // 3. For each triplet, communicate to the rightful owner that there's one
  //    more triplet for them. In doing that, retrieve the offset at which
  //    the triplet should be written on the remote.
  int num_my_triplets = 0;
  auto win = get_mpi_window (&num_my_triplets,1);
  std::vector<int> write_at(rows.size(),-1);
  check_mpi_call(MPI_Win_fence(0,win),"MPI_Win_fence");
  for (int i=0; i<nlweights; ++i) {
    // Tell pids[i] that we have one triplet for them. Also, get the current
    // value of num_triplets, cause that will tell us where to write when we
    // actually pass the triplet
    check_mpi_call(MPI_Get_accumulate(&one,1,MPI_INT,
                                      &write_at[i],1,MPI_INT,
                                      pids[i],0,1,MPI_INT,MPI_SUM,win),
                   "MPI_Get_accumulate");
  }
  check_mpi_call(MPI_Win_fence(0,win),"MPI_Win_fence");
  check_mpi_call(MPI_Win_free(&win),"MPI_Win_free"); 

  // 4. Ship each triplet to their owner

  // Create data type for a triplet
  auto mpi_gid_t = ekat::get_mpi_type<gid_type>();
  auto mpi_real_t = ekat::get_mpi_type<Real>();
  int lengths[3] = {1,1,1};
  MPI_Aint displacements[3] = {0, offsetof(Triplet,col), offsetof(Triplet,w)};
  MPI_Datatype types[3] = {mpi_gid_t,mpi_gid_t,mpi_real_t};
  MPI_Datatype triplet_mpi_t;
  MPI_Type_create_struct (3,lengths,displacements,types,&triplet_mpi_t);
  MPI_Type_commit(&triplet_mpi_t);

  // Create window, and do RMA stuff
  std::vector<Triplet> my_triplets (num_my_triplets);
  auto triplets_win = get_mpi_window(my_triplets.data(),my_triplets.size());
  check_mpi_call (MPI_Win_fence(0,triplets_win),"MPI_Win_fence");
  for (int i=0; i<nlweights; ++i) {
    Triplet t {rows[i], cols[i], S[i]};
    check_mpi_call (MPI_Put(&t,1,triplet_mpi_t,pids[i],
                            write_at[i],1,triplet_mpi_t,triplets_win),
                    "MPI_Put");
  }
  check_mpi_call (MPI_Win_fence(0,triplets_win),"MPI_Win_fence");
  check_mpi_call (MPI_Win_free(&triplets_win),  "MPI_Win_free");
  MPI_Type_free(&triplet_mpi_t);

  // Sort triplets by the row lids
  auto gid2lid = tgt_grid->get_gid2lid_map();
  auto compare = [&] (const Triplet& lhs, const Triplet& rhs) {
    return gid2lid.at(lhs.row) < gid2lid.at(rhs.row);
  };
  std::sort(my_triplets.begin(),my_triplets.end(),compare);

  return my_triplets;
}

void RefiningRemapperRMA::create_ov_src_fields ()
{
  using FL = FieldLayout;
  m_ov_src_fields.reserve(m_num_fields);
  const int num_ov_cols = m_ov_src_grid->get_num_local_dofs();
  const auto ov_gn = m_ov_src_grid->name();
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_tgt = m_tgt_fields[i];
    const auto& fid = f_tgt.get_header().get_identifier();
    auto tags = fid.get_layout().tags();
    auto dims = fid.get_layout().dims();
    dims[0] = num_ov_cols;
    FieldIdentifier ov_fid (fid.name(),FL(tags,dims),fid.get_units(),ov_gn,DataType::RealType);

    auto& ov_f = m_ov_src_fields.emplace_back(ov_fid);
    // Use same alloc props as tgt fields, to allow packing in local_mat_vec
    const auto pack_size = f_tgt.get_header().get_alloc_properties().get_largest_pack_size();
    ov_f.get_header().get_alloc_properties().request_allocation(pack_size);
    ov_f.allocate_view();
  }
}

void RefiningRemapperRMA::setup_mpi_data_structures ()
{
  using namespace ShortFieldTagsNames;

  // Extract some raw mpi info
  const auto mpi_comm  = m_comm.mpi_comm();
  check_mpi_call(MPI_Comm_group(mpi_comm,&m_mpi_group),"MPI_Comm_group");

  // Figure out where data needs to be retrieved from
  const auto ov_src_gids = m_ov_src_grid->get_dofs_gids().get_view<const gid_type*,Host>();
  m_src_grid->get_remote_pids_and_lids(ov_src_gids,m_remote_pids,m_remote_lids);

  // TODO: scope out possibility of using sub-groups for start/post calls
  //       (but I'm afraid you can't, b/c start/post may require same groups)

  // Create per-field structures
  constexpr auto COL = ShortFieldTagsNames::COL;
  m_mpi_win.resize(m_num_fields);
  m_col_size.resize(m_num_fields);
  m_col_stride.resize(m_num_fields);
  m_col_offset.resize(m_num_fields,0);
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f = m_src_fields[i];
    const auto& fh = f.get_header();
    const auto& layout = fh.get_identifier().get_layout();
    m_col_stride[i] = m_col_size[i] = layout.strip_dim(COL).size();



    // If field has a parent, col_stride and col_offset need to be adjusted
    auto p = fh.get_parent().lock();
    auto win_size = layout.size()*sizeof(Real);
    if (p) {
      EKAT_REQUIRE_MSG (p->get_parent().lock()==nullptr,
          "Error! We do not support remapping of subfields of other subfields.\n");
      const auto& sv_info = fh.get_alloc_properties().get_subview_info();
      m_col_stride[i] = sv_info.dim_extent * m_col_size[i];
      m_col_offset[i] = sv_info.slice_idx  * m_col_size[i];
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

  // Clear all fields
  m_src_fields.clear();
  m_tgt_fields.clear();
  m_ov_src_fields.clear();

  // Reset the state of the base class
  m_state = RepoState::Clean;
  m_num_fields = 0;
  m_num_registered_fields = 0;
  m_fields_are_bound.clear();
  m_num_bound_fields = 0;
}

void RefiningRemapperRMA::
check_mpi_call (int err, const std::string& context) const {
  EKAT_REQUIRE_MSG (err==MPI_SUCCESS,
      "Error! MPI operation encountered an error.\n"
      "  - err code: " + std::to_string(err) + "\n"
      "  - context: " + context + "\n");
}

} // namespace scream
