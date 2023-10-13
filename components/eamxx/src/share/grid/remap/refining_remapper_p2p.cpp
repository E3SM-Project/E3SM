#include "refining_remapper_p2p.hpp"

#include "share/grid/point_grid.hpp"
#include "share/grid/grid_import_export.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/util/scream_utils.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/ekat_pack_utils.hpp>

#include <numeric>

namespace scream
{

RefiningRemapperP2P::
RefiningRemapperP2P (const grid_ptr_type& tgt_grid,
                     const std::string& map_file)
 : AbstractRemapper()
 , m_comm (tgt_grid->get_comm())
{
  using namespace ShortFieldTagsNames;

  // Sanity checks
  EKAT_REQUIRE_MSG (tgt_grid->type()==GridType::Point,
      "Error! RefiningRemapperP2P only works on PointGrid grids.\n"
      "  - tgt grid name: " + tgt_grid->name() + "\n"
      "  - tgt_grid_type: " + e2str(tgt_grid->type()) + "\n");
  EKAT_REQUIRE_MSG (tgt_grid->is_unique(),
      "Error! RefiningRemapperP2P requires a unique target grid.\n");

  // This is a refining remapper. We only go in one direction
  m_bwd_allowed = false;

  // Load (i,j,w) triplets from map file, for all i that are
  // owned on the tgt_grid
  auto my_triplets = get_my_triplets (map_file,m_comm,tgt_grid,OwnedBy::Row);

  // Create an overlapped src map, consisting of all the col gids
  // in the triplets. This is overlapped, since for each gid there
  // may be 2+ ranks owning it.
  std::map<gid_type,int> ov_src_gid2lid;
  int num_my_triplets = my_triplets.size();
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
  m_row_offsets = view_1d<int>("",tgt_grid->get_num_local_dofs()+1);
  m_col_lids = view_1d<int>("",num_my_triplets);
  m_weights = view_1d<Real>("",num_my_triplets);

  auto row_offsets_h = Kokkos::create_mirror_view(m_row_offsets);
  auto col_lids_h    = Kokkos::create_mirror_view(m_col_lids);
  auto weights_h     = Kokkos::create_mirror_view(m_weights);

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

RefiningRemapperP2P::
~RefiningRemapperP2P ()
{
  clean_up();
}

FieldLayout RefiningRemapperP2P::
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
      EKAT_ERROR_MSG ("Layout not supported by RefiningRemapperP2P: " + e2str(lt) + "\n");
  }
  return src;
}

FieldLayout RefiningRemapperP2P::
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
      EKAT_ERROR_MSG ("Layout not supported by RefiningRemapperP2P: " + e2str(lt) + "\n");
  }
  return tgt;
}

void RefiningRemapperP2P::
do_register_field (const identifier_type& src, const identifier_type& tgt)
{
  constexpr auto COL = ShortFieldTagsNames::COL;
  EKAT_REQUIRE_MSG (src.get_layout().has_tag(COL),
      "Error! Cannot register a field without COL tag in RefiningRemapperP2P.\n"
      "  - field name: " + src.name() + "\n"
      "  - field layout: " + to_string(src.get_layout()) + "\n");
  m_src_fields.push_back(field_type(src));
  m_tgt_fields.push_back(field_type(tgt));
}

void RefiningRemapperP2P::
do_bind_field (const int ifield, const field_type& src, const field_type& tgt)
{
  EKAT_REQUIRE_MSG (src.data_type()==DataType::RealType,
      "Error! RefiningRemapperP2P only allows fields with RealType data.\n"
      "  - src field name: " + src.name() + "\n"
      "  - src field type: " + e2str(src.data_type()) + "\n");
  EKAT_REQUIRE_MSG (tgt.data_type()==DataType::RealType,
      "Error! RefiningRemapperP2P only allows fields with RealType data.\n"
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

void RefiningRemapperP2P::do_registration_ends ()
{
  if (this->m_num_bound_fields==this->m_num_registered_fields) {
    create_ov_src_fields ();
    setup_mpi_data_structures ();
  }
}

void RefiningRemapperP2P::do_remap_fwd ()
{
  std::cout << "remap start ...\n";
  // Fire the recv requests right away, so that if some other ranks
  // is done packing before us, we can start receiving their data
  if (not m_recv_req.empty()) {
    check_mpi_call(MPI_Startall(m_recv_req.size(),m_recv_req.data()),
                   "[RefiningRemapperP2P] starting persistent recv requests.\n");
  }

  // Do P2P communications
  std::cout << "pack_and_send start ...\n";
  pack_and_send ();
  std::cout << "recv_and_unpack start ...\n";
  recv_and_unpack ();

  // Perform local-mat vec
  // Helpef function, to establish if a field can be handled with packs
  auto can_pack_field = [](const Field& f) {
    const auto& ap = f.get_header().get_alloc_properties();
    return (ap.get_last_extent() % SCREAM_PACK_SIZE) == 0;
  };

  // Loop over each field, perform mat-vec
  std::cout << "mat-vec start ...\n";
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

  // Wait for all sends to be completed
  if (not m_send_req.empty()) {
    check_mpi_call(MPI_Waitall(m_send_req.size(),m_send_req.data(), MPI_STATUSES_IGNORE),
                   "[RefiningRemapperP2P] waiting on persistent send requests.\n");
  }
  std::cout << "remap completed\n";
}

template<int PackSize>
void RefiningRemapperP2P::
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

// auto RefiningRemapperP2P::
// get_my_triplets (const std::string& map_file,
//                  const grid_ptr_type& tgt_grid)
//  -> std::map<int,std::vector<Triplet>>
// {
//   using namespace ShortFieldTagsNames;

//   // 1. Load the map file chunking it evenly across all ranks
//   scorpio::register_file(map_file,scorpio::FileMode::Read);

//   // 1.1 Create a "helper" grid, with as many dofs as the number
//   //     of triplets in the map file, and divided linearly across ranks
//   const int ngweights = scorpio::get_dimlen(map_file,"n_s");
//   int nlweights = ngweights / m_comm.size();
//   if (m_comm.rank() < (ngweights % m_comm.size())) {
//     nlweights += 1;
//   }

//   gid_type offset = nlweights;
//   m_comm.scan(&offset,1,MPI_SUM);
//   offset -= nlweights; // scan is inclusive, but we need exclusive

//   // Create a unique decomp tag, which ensures all refining remappers have
//   // their own decomposition
//   static int tag_counter = 0;
//   const std::string int_decomp_tag  = "RR::gmtg,int,grid-idx=" + std::to_string(tag_counter++);
//   const std::string real_decomp_tag = "RR::gmtg,real,grid-idx=" + std::to_string(tag_counter++);

//   // 1.2 Read a chunk of triplets col indices
//   std::vector<gid_type> cols(nlweights);
//   std::vector<gid_type> rows(nlweights);
//   std::vector<Real>  S(nlweights);

//   scorpio::register_variable(map_file, "col", "col", {"n_s"}, "int",  int_decomp_tag);
//   scorpio::register_variable(map_file, "row", "row", {"n_s"}, "int",  int_decomp_tag);
//   scorpio::register_variable(map_file, "S",   "S",   {"n_s"}, "real", real_decomp_tag);

//   std::vector<scorpio::offset_t> dofs_offsets(nlweights);
//   std::iota(dofs_offsets.begin(),dofs_offsets.end(),offset);
//   scorpio::set_dof(map_file,"col",nlweights,dofs_offsets.data());
//   scorpio::set_dof(map_file,"row",nlweights,dofs_offsets.data());
//   scorpio::set_dof(map_file,"S"  ,nlweights,dofs_offsets.data());
//   scorpio::set_decomp(map_file);

//   scorpio::grid_read_data_array(map_file,"col",-1,cols.data(),cols.size());
//   scorpio::grid_read_data_array(map_file,"row",-1,rows.data(),rows.size());
//   scorpio::grid_read_data_array(map_file,"S"  ,-1,S.data(),S.size());

//   scorpio::eam_pio_closefile(map_file);

//   // 1.3 Dofs in tgt grid are likely 0-based, while row ids in map file
//   // are likely 1-based. To match dofs, we need to offset the row
//   // ids we read in.
//   int map_file_min_row = std::numeric_limits<int>::max();
//   for (int id=0; id<nlweights; id++) {
//     map_file_min_row = std::min(rows[id],map_file_min_row);
//   }
//   int global_map_file_min_row;
//   m_comm.all_reduce(&map_file_min_row,&global_map_file_min_row,1,MPI_MIN);

//   gid_type row_offset = global_map_file_min_row - tgt_grid->get_global_min_dof_gid();
//   for (auto& id : rows) {
//     id -= row_offset;
//   }

//   // Create a grid based on the row gids I read in (may be duplicated across ranks)
//   std::vector<gid_type> unique_rows;
//   for (auto row : rows) {
//     if (not ekat::contains(unique_rows,row)) {
//       unique_rows.push_back(row);
//     }
//   }
//   auto io_grid = std::make_shared<PointGrid> ("helper",unique_rows.size(),0,m_comm);
//   auto io_grid_gids_h = io_grid->get_dofs_gids().get_view<gid_type*,Host>();
//   int k = 0;
//   for (auto row : rows) {
//     io_grid_gids_h(k++) = row;
//   }
//   io_grid->get_dofs_gids().sync_to_dev();

//   // Create Triplets to export, sorted by row
//   std::map<int,std::vector<Triplet>> io_triplets;
//   auto io_grid_gid2lid = io_grid->get_gid2lid_map();
//   for (int i=0; i<nlweights; ++i) {
//     auto row = rows[i];
//     auto io_lid = io_grid_gid2lid[row];
//     io_triplets[io_lid].emplace_back(rows[i], cols[i], S[i]);
//   }

//   // Create data type for a triplet
//   auto mpi_gid_t = ekat::get_mpi_type<gid_type>();
//   auto mpi_real_t = ekat::get_mpi_type<Real>();
//   int lengths[3] = {1,1,1};
//   MPI_Aint displacements[3] = {0, offsetof(Triplet,col), offsetof(Triplet,w)};
//   MPI_Datatype types[3] = {mpi_gid_t,mpi_gid_t,mpi_real_t};
//   MPI_Datatype mpi_triplet_t;
//   MPI_Type_create_struct (3,lengths,displacements,types,&mpi_triplet_t);
//   MPI_Type_commit(&mpi_triplet_t);

//   // Create import-export
//   GridImportExport imp_exp (tgt_grid,io_grid);
//   std::map<int,std::vector<Triplet>> my_triplets;
//   imp_exp.gather(mpi_triplet_t,io_triplets,my_triplets);

//   return my_triplets;
// }

void RefiningRemapperP2P::create_ov_src_fields ()
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

void RefiningRemapperP2P::setup_mpi_data_structures ()
{
  using namespace ShortFieldTagsNames;

  const int nranks = m_comm.size();

  // Get cumulative col size of each field (to be used to compute offsets)
  m_fields_col_sizes_scan_sum.resize(m_num_fields+1,0);
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f = m_src_fields[i];
    const auto& fl = f.get_header().get_identifier().get_layout();
    const auto& col_size = fl.strip_dim(COL).size();
    m_fields_col_sizes_scan_sum[i+1] = m_fields_col_sizes_scan_sum[i] + col_size;
  }
  auto total_col_size = m_fields_col_sizes_scan_sum.back();

  // ----------- Compute RECV metadata -------------- //

  // Figure out where ov_src cols are received from
  const int ncols_recv = m_ov_src_grid->get_num_local_dofs();

  m_imp_exp = std::make_shared<GridImportExport>(m_src_grid,m_ov_src_grid);

  // We can now compute the offset of each pid in the recv buffer
  m_pids_recv_offsets = view_1d<int>("",nranks+1);
  auto ncols_recv_h = m_imp_exp->num_imports_per_pid_h();
  auto pids_recv_offsets_h = Kokkos::create_mirror_view(m_pids_recv_offsets);
  pids_recv_offsets_h[0] = 0;
  for (int pid=0; pid<nranks; ++pid) {
    pids_recv_offsets_h(pid+1) = pids_recv_offsets_h(pid)
                               + ncols_recv_h(pid)*total_col_size;
  }
  Kokkos::deep_copy(m_pids_recv_offsets,pids_recv_offsets_h);

  // Create the recv buffer(s)
  auto recv_buf_size = ncols_recv*total_col_size;
  m_recv_buffer = decltype(m_recv_buffer)("RefiningRemapperP2P::recv_buf",recv_buf_size);
  m_mpi_recv_buffer = Kokkos::create_mirror_view(decltype(m_mpi_recv_buffer)::execution_space(),m_recv_buffer);

  // ----------- Compute SEND metadata -------------- //

  m_pids_send_offsets = view_1d<int>("",nranks+1);
  auto pids_send_offsets_h = Kokkos::create_mirror_view(m_pids_send_offsets);
  auto ncols_send_h = m_imp_exp->num_exports_per_pid_h();
  pids_send_offsets_h[0] = 0;
  for (int pid=0; pid<nranks; ++pid) {
    pids_send_offsets_h(pid+1) = pids_send_offsets_h(pid)
                               + ncols_send_h(pid)*total_col_size;
  }
  Kokkos::deep_copy(m_pids_send_offsets,pids_send_offsets_h);

  // Create the send buffer(s)
  auto send_buf_size = pids_send_offsets_h(nranks)*total_col_size;
  m_send_buffer = decltype(m_send_buffer)("RefiningRemapperP2P::send_buf",send_buf_size);
  m_mpi_send_buffer = Kokkos::create_mirror_view(decltype(m_mpi_send_buffer)::execution_space(),m_send_buffer);

  // ----------- Create Requests ------------ //

  const auto mpi_comm = m_comm.mpi_comm();
  const auto mpi_real = ekat::get_mpi_type<Real>();
  for (int pid=0; pid<nranks; ++pid) {
    // Send request
    if (ncols_send_h(pid)>0) {
      auto send_ptr = m_mpi_send_buffer.data() + pids_send_offsets_h(pid);
      auto send_count = ncols_send_h(pid)*total_col_size;
      auto& req = m_send_req.emplace_back();
      MPI_Send_init (send_ptr, send_count, mpi_real, pid,
                     0, mpi_comm, &req);
    }
    // Recv request
    if (ncols_recv_h(pid)>0) {
      auto recv_ptr = m_mpi_recv_buffer.data() + pids_recv_offsets_h(pid);
      auto recv_count = ncols_recv_h(pid)*total_col_size;
      auto& req = m_recv_req.emplace_back();
      MPI_Recv_init (recv_ptr, recv_count, mpi_real, pid,
                     0, mpi_comm, &req);
    }
  }
}

void RefiningRemapperP2P::pack_and_send ()
{
  using RangePolicy = typename KT::RangePolicy;
  using TeamMember  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  auto export_pids = m_imp_exp->export_pids();
  auto export_lids = m_imp_exp->export_lids();
  auto ncols_send  = m_imp_exp->num_exports_per_pid();
  auto pids_send_offsets = m_pids_send_offsets;
  auto send_buf = m_send_buffer;
  const int num_exports = export_pids.size();
  const int total_col_size = m_fields_col_sizes_scan_sum.back();
  for (int ifield=0; ifield<m_num_fields; ++ifield) {
    const auto& f = m_src_fields[ifield];
    const auto& fl = f.get_header().get_identifier().get_layout();
    const auto f_col_sizes_scan_sum = m_fields_col_sizes_scan_sum[ifield];
    switch (fl.rank()) {
      case 1:
      {
        const auto v = f.get_strided_view<const Real*>();
        auto pack = KOKKOS_LAMBDA(const int iexp) {
          auto pid = export_pids(iexp);
          auto icol = export_lids(iexp);
          auto pid_offset = pids_send_offsets(pid);
          auto pos_within_pid = iexp - (pid_offset / total_col_size);
          auto offset = pid_offset
                      + ncols_send(pid)*f_col_sizes_scan_sum
                      + pos_within_pid;
          send_buf(offset) = v(icol);
        };
        Kokkos::parallel_for(RangePolicy(0,num_exports),pack);
        break;
      }
      case 2:
      {
        const auto v = f.get_view<const Real**>();
        const int dim1 = fl.dim(1);
        auto policy = ESU::get_default_team_policy(num_exports,dim1);
        auto pack = KOKKOS_LAMBDA(const TeamMember& team) {
          const int iexp = team.league_rank();
          const int icol = export_lids(iexp);
          const int pid  = export_pids(iexp);
          auto pid_offset = pids_send_offsets(pid);
          auto pos_within_pid = iexp - (pid_offset / total_col_size);
          auto offset = pid_offset
                      + ncols_send(pid)*f_col_sizes_scan_sum
                      + pos_within_pid;
          auto col_pack = [&](const int& k) {
            send_buf(offset+k) = v(icol,k);
          };
          auto tvr = Kokkos::TeamVectorRange(team,dim1);
          Kokkos::parallel_for(tvr,col_pack);
        };
        Kokkos::parallel_for(policy,pack);
        break;
      }
      case 3:
      {
        const auto v = f.get_view<const Real***>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dim(2);
        auto policy = ESU::get_default_team_policy(num_exports,dim1*dim2);
        auto pack = KOKKOS_LAMBDA(const TeamMember& team) {
          const int iexp = team.league_rank();
          const int icol = export_pids(iexp);
          const int  pid = export_pids(iexp);
          auto pid_offset = pids_send_offsets(pid);
          auto pos_within_pid = iexp - (pid_offset / total_col_size);
          auto offset = pid_offset
                      + ncols_send(pid)*f_col_sizes_scan_sum
                      + pos_within_pid;
          auto col_pack = [&](const int& idx) {
            const int j = idx / dim2;
            const int k = idx % dim2;
            send_buf(offset+idx) = v(icol,j,k);
          };
          auto tvr = Kokkos::TeamVectorRange(team,dim1*dim2);
          Kokkos::parallel_for(tvr,col_pack);
        };
        Kokkos::parallel_for(policy,pack);
        break;
      }
      default:
        EKAT_ERROR_MSG ("Unexpected field rank in RefiningRemapperP2P::pack.\n"
            "  - MPI rank  : " + std::to_string(m_comm.rank()) + "\n"
            "  - field name: " + f.name() + "\n"
            "  - field rank: " + std::to_string(fl.rank()) + "\n");
    }
  }

  // Wait for all threads to be done packing
  Kokkos::fence();

  // If MPI does not use dev pointers, we need to deep copy from dev to host
  if (not MpiOnDev) {
    Kokkos::deep_copy (m_mpi_send_buffer,m_send_buffer);
  }

  if (not m_send_req.empty()) {
    check_mpi_call(MPI_Startall(m_send_req.size(),m_send_req.data()),
                   "[RefiningRemapperP2P] start persistent send requests.\n");
  }
}

void RefiningRemapperP2P::recv_and_unpack ()
{
  if (not m_recv_req.empty()) {
    check_mpi_call(MPI_Waitall(m_recv_req.size(),m_recv_req.data(), MPI_STATUSES_IGNORE),
                   "[RefiningRemapperP2P] waiting on persistent recv requests.\n");
  }
  // If MPI does not use dev pointers, we need to deep copy from host to dev
  if (not MpiOnDev) {
    Kokkos::deep_copy (m_recv_buffer,m_mpi_recv_buffer);
  }

  using RangePolicy = typename KT::RangePolicy;
  using TeamMember  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  auto import_pids = m_imp_exp->import_pids();
  auto import_lids = m_imp_exp->import_lids();
  auto ncols_recv  = m_imp_exp->num_imports_per_pid();
  auto pids_recv_offsets = m_pids_recv_offsets;
  auto recv_buf = m_recv_buffer;
  const int num_imports = import_pids.size();
  const int total_col_size = m_fields_col_sizes_scan_sum.back();
  for (int ifield=0; ifield<m_num_fields; ++ifield) {
          auto& f  = m_ov_src_fields[ifield];
    const auto& fl = f.get_header().get_identifier().get_layout();
    const auto f_col_sizes_scan_sum = m_fields_col_sizes_scan_sum[ifield];
    switch (fl.rank()) {
      case 1:
      {
        auto v = f.get_view<Real*>();
        auto unpack = KOKKOS_LAMBDA (const int idx) {
          const int pid  = import_pids(idx);
          const int icol = import_lids(idx);
          const auto pid_offset = pids_recv_offsets(pid);
          const auto pos_within_pid = idx - (pid_offset / total_col_size);
          auto offset = pid_offset
                      + ncols_recv(pid)*f_col_sizes_scan_sum
                      + pos_within_pid;
          v(icol) = recv_buf(offset);
        };
        Kokkos::parallel_for(RangePolicy(0,num_imports),unpack);
        break;
      }
      case 2:
      {
        auto v = f.get_view<Real**>();
        const int dim1 = fl.dim(1);
        auto policy = ESU::get_default_team_policy(num_imports,dim1);
        auto unpack = KOKKOS_LAMBDA (const TeamMember& team) {
          const int idx  = team.league_rank();
          const int pid  = import_pids(idx);
          const int icol = import_lids(idx);
          const auto pid_offset = pids_recv_offsets(pid);
          const auto pos_within_pid = idx - (pid_offset / total_col_size);
          auto offset = pid_offset
                      + ncols_recv(pid)*f_col_sizes_scan_sum
                      + pos_within_pid;

          auto col_pack = [&](const int& k) {
            v(icol,k) = recv_buf(offset+k);
          };
          auto tvr = Kokkos::TeamVectorRange(team,dim1);
          Kokkos::parallel_for(tvr,col_pack);
        };
        Kokkos::parallel_for(policy,unpack);
        break;
      }
      case 3:
      {
        auto v = f.get_view<Real***>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dim(2);
        auto policy = ESU::get_default_team_policy(num_imports,dim1*dim2);
        auto unpack = KOKKOS_LAMBDA (const TeamMember& team) {
          const int idx  = team.league_rank();
          const int pid  = import_pids(idx);
          const int icol = import_lids(idx);
          const auto pid_offset = pids_recv_offsets(pid);
          const auto pos_within_pid = idx - (pid_offset / total_col_size);
          auto offset = pid_offset
                            + ncols_recv(pid)*f_col_sizes_scan_sum
                            + pos_within_pid;

          auto col_pack = [&](const int& idx) {
            const int j = idx / dim2;
            const int k = idx % dim2;
            v(icol,j,k) = recv_buf(offset+idx);
          };
          auto tvr = Kokkos::TeamVectorRange(team,dim1);
          Kokkos::parallel_for(tvr,col_pack);
        };
        Kokkos::parallel_for(policy,unpack);
        break;
      }
      default:
        EKAT_ERROR_MSG ("Unexpected field rank in RefiningRemapperP2P::unpack.\n"
            "  - MPI rank  : " + std::to_string(m_comm.rank()) + "\n"
            "  - field name: " + f.name() + "\n"
            "  - field rank: " + std::to_string(fl.rank()) + "\n");
    }
  }
}

void RefiningRemapperP2P::clean_up ()
{
  // Clear all MPI related structures
  m_send_buffer         = view_1d<Real>();
  m_recv_buffer         = view_1d<Real>();
  m_mpi_send_buffer     = mpi_view_1d<Real>();
  m_mpi_recv_buffer     = mpi_view_1d<Real>();
  m_send_req.clear();
  m_recv_req.clear();
  m_imp_exp = nullptr;

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

} // namespace scream
