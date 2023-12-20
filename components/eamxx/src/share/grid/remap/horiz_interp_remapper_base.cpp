#include "horiz_interp_remapper_base.hpp"

#include "share/grid/point_grid.hpp"
#include "share/grid/grid_import_export.hpp"
#include "share/io/scorpio_input.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/ekat_pack_utils.hpp>
#include <numeric>

namespace scream
{

HorizInterpRemapperBase::
HorizInterpRemapperBase (const grid_ptr_type& fine_grid,
                         const std::string& map_file,
                         const InterpType type)
 : m_fine_grid(fine_grid)
 , m_type (type)
 , m_comm (fine_grid->get_comm())
{
  // Sanity checks
  EKAT_REQUIRE_MSG (fine_grid->type()==GridType::Point,
      "Error! Horizontal interpolatory remap only works on PointGrid grids.\n"
      "  - fine grid name: " + fine_grid->name() + "\n"
      "  - fine_grid_type: " + e2str(fine_grid->type()) + "\n");
  EKAT_REQUIRE_MSG (fine_grid->is_unique(),
      "Error! CoarseningRemapper requires a unique source grid.\n");

  // This is a special remapper. We only go in one direction
  m_bwd_allowed = false;

  // Read the map file, loading the triplets this rank needs for the crs matrix
  // in the map file that this rank has to read
  auto my_triplets = get_my_triplets (map_file);

  // Create coarse/ov_coarse grids
  create_coarse_grids (my_triplets);

  // Set src/tgt grid, based on interpolation type
  if (m_type==InterpType::Refine) {
    set_grids (m_coarse_grid,m_fine_grid);
  } else {
    set_grids (m_fine_grid,m_coarse_grid);
  }

  // Create crs matrix
  create_crs_matrix_structures (my_triplets);
}

FieldLayout HorizInterpRemapperBase::
create_src_layout (const FieldLayout& tgt_layout) const
{
  EKAT_REQUIRE_MSG (m_src_grid!=nullptr,
      "Error! Cannot create source layout until the source grid has been set.\n");

  EKAT_REQUIRE_MSG (is_valid_tgt_layout(tgt_layout),
      "[HorizInterpRemapperBase] Error! Input target layout is not valid for this remapper.\n"
      " - input layout: " + to_string(tgt_layout));

  using namespace ShortFieldTagsNames;
  const auto lt = get_layout_type(tgt_layout.tags());
  const bool midpoints = tgt_layout.has_tag(LEV);
  const int vec_dim = tgt_layout.is_vector_layout() ? tgt_layout.dim(CMP) : -1;
  auto src = FieldLayout::invalid();
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

FieldLayout HorizInterpRemapperBase::
create_tgt_layout (const FieldLayout& src_layout) const
{
  EKAT_REQUIRE_MSG (m_tgt_grid!=nullptr,
      "Error! Cannot create target layout until the target grid has been set.\n");

  EKAT_REQUIRE_MSG (is_valid_src_layout(src_layout),
      "[HorizInterpRemapperBase] Error! Input source layout is not valid for this remapper.\n"
      " - input layout: " + to_string(src_layout));

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
      EKAT_ERROR_MSG ("Layout not supported by CoarseningRemapper: " + e2str(lt) + "\n");
  }
  return tgt;
}

void HorizInterpRemapperBase::do_registration_ends ()
{
  if (this->m_num_bound_fields==this->m_num_registered_fields) {
    create_ov_fields ();
    setup_mpi_data_structures ();
  }
}

void HorizInterpRemapperBase::
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

void HorizInterpRemapperBase::
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
    create_ov_fields ();
    setup_mpi_data_structures ();
  }
}

auto HorizInterpRemapperBase::
get_my_triplets (const std::string& map_file) const
 -> std::vector<Triplet>
{
  using gid_type = AbstractGrid::gid_type;
  using namespace ShortFieldTagsNames;

  // 1. Load the map file chunking it evenly across all ranks
  scorpio::register_file(map_file,scorpio::FileMode::Read);

  // 1.1 Create a "helper" grid, with as many dofs as the number
  //     of triplets in the map file, and divided linearly across ranks
  const int ngweights = scorpio::get_dimlen(map_file,"n_s");
  int nlweights = ngweights / m_comm.size();
  if (m_comm.rank() < (ngweights % m_comm.size())) {
    nlweights += 1;
  }

  gid_type offset = nlweights;
  m_comm.scan(&offset,1,MPI_SUM);
  offset -= nlweights; // scan is inclusive, but we need exclusive

  // Create a unique decomp tag, which ensures all refining remappers have
  // their own decomposition
  static int tag_counter = 0;
  const std::string int_decomp_tag  = "RR::gmtg,int,grid-idx=" + std::to_string(tag_counter++);
  const std::string real_decomp_tag = "RR::gmtg,real,grid-idx=" + std::to_string(tag_counter++);

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

  // 1.3 Dofs in grid are likely 0-based, while row/col ids in map file
  // are likely 1-based. To match dofs, we need to offset the row/cols
  // ids we just read in.
  int map_file_min_row = std::numeric_limits<int>::max();
  int map_file_min_col = std::numeric_limits<int>::max();
  for (int id=0; id<nlweights; id++) {
    map_file_min_row = std::min(rows[id],map_file_min_row);
    map_file_min_col = std::min(cols[id],map_file_min_col);
  }
  int global_map_file_min_row, global_map_file_min_col;
  m_comm.all_reduce(&map_file_min_row,&global_map_file_min_row,1,MPI_MIN);
  m_comm.all_reduce(&map_file_min_col,&global_map_file_min_col,1,MPI_MIN);

  gid_type row_offset = global_map_file_min_row;
  gid_type col_offset = global_map_file_min_col;
  if (m_type==InterpType::Refine) {
    row_offset -= m_fine_grid->get_global_min_dof_gid();
  } else {
    col_offset -= m_fine_grid->get_global_min_dof_gid();
  }
  for (auto& id : rows) {
    id -= row_offset;
  }
  for (auto& id : cols) {
    id -= col_offset;
  }

  // Create a grid based on the row gids I read in (may be duplicated across ranks)
  std::vector<gid_type> unique_gids;
  const auto& gids = m_type==InterpType::Refine ? rows : cols;
  for (auto gid : gids) {
    if (not ekat::contains(unique_gids,gid)) {
      unique_gids.push_back(gid);
    }
  }
  auto io_grid = std::make_shared<PointGrid> ("helper",unique_gids.size(),0,m_comm);
  auto io_grid_gids_h = io_grid->get_dofs_gids().get_view<gid_type*,Host>();
  int k = 0;
  for (auto gid : unique_gids) {
    io_grid_gids_h(k++) = gid;
  }
  io_grid->get_dofs_gids().sync_to_dev();

  // Create Triplets to export, sorted by gid
  std::map<int,std::vector<Triplet>> io_triplets;
  auto io_grid_gid2lid = io_grid->get_gid2lid_map();
  for (int i=0; i<nlweights; ++i) {
    auto gid = gids[i];
    auto io_lid = io_grid_gid2lid[gid];
    io_triplets[io_lid].emplace_back(rows[i], cols[i], S[i]);
  }

  // Create data type for a triplet
  auto mpi_gid_t = ekat::get_mpi_type<gid_type>();
  auto mpi_real_t = ekat::get_mpi_type<Real>();
  int lengths[3] = {1,1,1};
  MPI_Aint displacements[3] = {0, offsetof(Triplet,col), offsetof(Triplet,w)};
  MPI_Datatype types[3] = {mpi_gid_t,mpi_gid_t,mpi_real_t};
  MPI_Datatype mpi_triplet_t;
  MPI_Type_create_struct (3,lengths,displacements,types,&mpi_triplet_t);
  MPI_Type_commit(&mpi_triplet_t);

  // Create import-export
  GridImportExport imp_exp (m_fine_grid,io_grid);
  std::map<int,std::vector<Triplet>> my_triplets_map;
  imp_exp.gather(mpi_triplet_t,io_triplets,my_triplets_map);
  MPI_Type_free(&mpi_triplet_t);

  std::vector<Triplet> my_triplets;
  for (auto& it : my_triplets_map) {
    my_triplets.reserve(my_triplets.size()+it.second.size());
    std::move(it.second.begin(),it.second.end(),std::back_inserter(my_triplets));
  }

  return my_triplets;
}

void HorizInterpRemapperBase::
create_coarse_grids (const std::vector<Triplet>& triplets)
{
  const int nlevs = m_fine_grid->get_num_vertical_levels();

  // Gather overlapped coarse grid gids (rows or cols, depending on m_type)
  std::map<gid_type,int> ov_gid2lid;
  bool pickRow = m_type==InterpType::Coarsen;
  for (const auto& t : triplets) {
    ov_gid2lid.emplace(pickRow ? t.row : t.col,ov_gid2lid.size());
  }
  int num_ov_gids = ov_gid2lid.size();

  // Use a temp and then assing, b/c grid_ptr_type is a pointer to const,
  // so you can't modify gids using that pointer
  auto ov_coarse_grid = std::make_shared<PointGrid>("ov_coarse_grid",num_ov_gids,nlevs,m_comm);
  auto ov_coarse_gids_h = ov_coarse_grid->get_dofs_gids().get_view<gid_type*,Host>();
  for (const auto& it : ov_gid2lid) {
    ov_coarse_gids_h[it.second] = it.first;
  }
  ov_coarse_grid->get_dofs_gids().sync_to_dev();
  m_ov_coarse_grid = ov_coarse_grid;

  // Create the unique coarse grid
  auto coarse_gids = m_ov_coarse_grid->get_unique_gids();
  int num_gids = coarse_gids.size();
  m_coarse_grid = std::make_shared<PointGrid>("coarse_grid",num_gids,nlevs,m_comm);
  auto coarse_gids_h = m_coarse_grid->get_dofs_gids().get_view<gid_type*,Host>();
  std::copy(coarse_gids.begin(),coarse_gids.end(),coarse_gids_h.data());
  m_coarse_grid->get_dofs_gids().sync_to_dev();
}

void HorizInterpRemapperBase::
create_crs_matrix_structures (std::vector<Triplet>& triplets)
{
  // Get row/col data depending on interp type
  bool refine = m_type==InterpType::Refine;
  auto row_grid = refine ? m_fine_grid : m_ov_coarse_grid;
  auto col_grid = refine ? m_ov_coarse_grid : m_fine_grid;
  const int num_rows = row_grid->get_num_local_dofs();

  auto col_gid2lid = col_grid->get_gid2lid_map();
  auto row_gid2lid = row_grid->get_gid2lid_map();

  // Sort triplets so that row GIDs appear in the same order as
  // in the row grid. If two row GIDs are the same, use same logic
  // with col
  auto compare = [&] (const Triplet& lhs, const Triplet& rhs) {
    auto lhs_lrow = row_gid2lid.at(lhs.row);
    auto rhs_lrow = row_gid2lid.at(rhs.row);
    auto lhs_lcol = col_gid2lid.at(lhs.col);
    auto rhs_lcol = col_gid2lid.at(rhs.col);
    return lhs_lrow<rhs_lrow or (lhs_lrow==rhs_lrow and lhs_lcol<rhs_lcol);
  };
  std::sort(triplets.begin(),triplets.end(),compare);

  // Alloc views and create mirror views
  const int nnz = triplets.size();
  m_row_offsets = view_1d<int>("",num_rows+1);
  m_col_lids    = view_1d<int>("",nnz);
  m_weights     = view_1d<Real>("",nnz);

  auto row_offsets_h = Kokkos::create_mirror_view(m_row_offsets);
  auto col_lids_h    = Kokkos::create_mirror_view(m_col_lids);
  auto weights_h     = Kokkos::create_mirror_view(m_weights);

  // Fill col ids and weights
  for (int i=0; i<nnz; ++i) {
    col_lids_h(i) = col_gid2lid[triplets[i].col];
    weights_h(i)  = triplets[i].w;
  }
  Kokkos::deep_copy(m_weights,weights_h);
  Kokkos::deep_copy(m_col_lids,col_lids_h);

  // Compute row offsets
  std::vector<int> row_counts(num_rows);
  for (int i=0; i<nnz; ++i) {
    ++row_counts[row_gid2lid[triplets[i].row]];
  }
  std::partial_sum(row_counts.begin(),row_counts.end(),row_offsets_h.data()+1);
  EKAT_REQUIRE_MSG (
      row_offsets_h(num_rows)==nnz,
      "Error! Something went wrong while computing row offsets.\n"
      "  - local nnz       : " + std::to_string(nnz) + "\n"
      "  - row_offsets(end): " + std::to_string(row_offsets_h(num_rows)) + "\n");

  Kokkos::deep_copy(m_row_offsets,row_offsets_h);
}

void HorizInterpRemapperBase::create_ov_fields ()
{
  m_ov_fields.reserve(m_num_fields);
  const auto num_ov_gids = m_ov_coarse_grid->get_num_local_dofs();
  const auto ov_gn = m_ov_coarse_grid->name();
  const auto dt = DataType::RealType;
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f = m_type==InterpType::Refine ? m_tgt_fields[i] : m_src_fields[i];
    const auto& fid = f.get_header().get_identifier();
    auto layout = fid.get_layout();
    layout.set_dimension(0,num_ov_gids);
    FieldIdentifier ov_fid (fid.name(),layout,fid.get_units(),ov_gn,dt);

    auto& ov_f = m_ov_fields.emplace_back(ov_fid);

    // Use same alloc props as fine fields, to allow packing in local_mat_vec
    const auto pack_size = f.get_header().get_alloc_properties().get_largest_pack_size();
    ov_f.get_header().get_alloc_properties().request_allocation(pack_size);
    ov_f.allocate_view();
  }
}

template<int PackSize>
void HorizInterpRemapperBase::
local_mat_vec (const Field& x, const Field& y) const
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  using Pack        = ekat::Pack<Real,PackSize>;
  using PackInfo    = ekat::PackInfo<PackSize>;

  const auto row_grid = m_type==InterpType::Refine ? m_fine_grid : m_ov_coarse_grid;
  const int  nrows    = row_grid->get_num_local_dofs();

  const auto& src_layout = x.get_header().get_identifier().get_layout();
  const int   rank       = src_layout.rank();

  auto row_offsets = m_row_offsets;
  auto col_lids    = m_col_lids;
  auto weights     = m_weights;

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

void HorizInterpRemapperBase::clean_up ()
{
  // Clear all fields
  m_src_fields.clear();
  m_tgt_fields.clear();
  m_ov_fields.clear();

  // Reset the state of the base class
  m_state = RepoState::Clean;
  m_num_fields = 0;
  m_num_registered_fields = 0;
  m_fields_are_bound.clear();
  m_num_bound_fields = 0;
}

// ETI, so derived classes can call this method
template
void HorizInterpRemapperBase::
local_mat_vec<1>(const Field&, const Field&) const;

#if SCREAM_PACK_SIZE>1
template
void HorizInterpRemapperBase::
local_mat_vec<SCREAM_PACK_SIZE>(const Field&, const Field&) const;
#endif

} // namespace scream
