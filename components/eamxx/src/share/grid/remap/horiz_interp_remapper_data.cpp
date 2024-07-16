#include "horiz_interp_remapper_data.hpp"

#include "share/grid/point_grid.hpp"
#include "share/grid/grid_import_export.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include <numeric>

namespace scream {

// --------------- HorizRemapperData ---------------- //

void HorizRemapperData::
build (const std::string& map_file,
       const std::shared_ptr<const AbstractGrid>& fine_grid_in,
       const ekat::Comm& comm_in,
       const InterpType type_in)
{
  comm = comm_in;
  fine_grid = fine_grid_in;
  type = type_in;

  // Gather sparse matrix triplets needed by this rank
  auto my_triplets = get_my_triplets (map_file);

  // Create coarse/ov_coarse grids
  create_coarse_grids (my_triplets);

  // Create crs matrix
  create_crs_matrix_structures (my_triplets);
}

auto HorizRemapperData::
get_my_triplets (const std::string& map_file) const
 -> std::vector<Triplet>
{
  using gid_type = AbstractGrid::gid_type;
  using namespace ShortFieldTagsNames;

  // 1. Load the map file chunking it evenly across all ranks
  scorpio::register_file(map_file,scorpio::FileMode::Read);

  // Inform scorpio that we will provide "int" pointers for row/col indices
  scorpio::change_var_dtype(map_file,"row","int");
  scorpio::change_var_dtype(map_file,"col","int");

  // Decompose n_s dim linearly across ranks
  scorpio::set_dim_decomp (map_file,"n_s");
  int nlweights = scorpio::get_dimlen_local(map_file,"n_s");

  // 1.1 Read a chunk of triplets col indices
  // NOTE: add 1 so that we don't pass nullptr to scorpio read routines (which would trigger
  //       a runtime error). Don't worry though: we never access the last entry of these vectors
  std::vector<gid_type> cols(nlweights+1,-1);
  std::vector<gid_type> rows(nlweights+1,-1);
  std::vector<Real>  S(nlweights+1,0);

  // Figure out if we are reading the right map, that is:
  //  - n_a or n_b matches the fine grid ncols
  //  - the map "direction" (fine->coarse or coarse->fine) matches m_type
  const int n_a = scorpio::get_dimlen(map_file,"n_a");
  const int n_b = scorpio::get_dimlen(map_file,"n_b");
  const int ncols_fine = fine_grid->get_num_global_dofs();
  EKAT_REQUIRE_MSG (n_a==ncols_fine or n_b==ncols_fine,
      "Error! The input map seems incompatible with the remapper fine grid.\n"
      " - map file: " + map_file + "\n"
      " - map file n_a: " + std::to_string(n_a) + "\n"
      " - map file n_b: " + std::to_string(n_b) + "\n"
      " - fine grid ncols: " + std::to_string(ncols_fine) + "\n");
  const bool map_is_coarsening = n_a==ncols_fine;
  EKAT_REQUIRE_MSG (map_is_coarsening==(type==InterpType::Coarsen),
      "Error! The input map seems incompatible with the remapper type.\n"
      " - map file: " + map_file + "\n"
      " - map file n_a: " + std::to_string(n_a) + "\n"
      " - map file n_b: " + std::to_string(n_b) + "\n"
      " - fine grid ncols: " + std::to_string(ncols_fine) + "\n"
      " - remapper type: " + std::string(type==InterpType::Refine ? "refine" : "coarsen") + "\n");

  scorpio::read_var(map_file,"col",cols.data());
  scorpio::read_var(map_file,"row",rows.data());
  scorpio::read_var(map_file,"S"  ,S.data());

  // Previously, we added 1 to their length, to avoid nullptr in scorpio::read.
  // However, we later do range loops on these vectors, so resize them back to nlweights
  cols.resize(nlweights);
  rows.resize(nlweights);
  S.resize(nlweights);

  scorpio::release_file(map_file);

  // 1.2 Dofs in grid are likely 0-based, while row/col ids in map file
  // are likely 1-based. To match dofs, we need to offset the row/cols
  // ids we just read in.
  int map_file_min_row = std::numeric_limits<int>::max();
  int map_file_min_col = std::numeric_limits<int>::max();
  for (int id=0; id<nlweights; id++) {
    map_file_min_row = std::min(rows[id],map_file_min_row);
    map_file_min_col = std::min(cols[id],map_file_min_col);
  }
  int global_map_file_min_row, global_map_file_min_col;
  comm.all_reduce(&map_file_min_row,&global_map_file_min_row,1,MPI_MIN);
  comm.all_reduce(&map_file_min_col,&global_map_file_min_col,1,MPI_MIN);

  gid_type row_offset = global_map_file_min_row;
  gid_type col_offset = global_map_file_min_col;
  if (type==InterpType::Refine) {
    row_offset -= fine_grid->get_global_min_dof_gid();
  } else {
    col_offset -= fine_grid->get_global_min_dof_gid();
  }
  for (auto& id : rows) {
    id -= row_offset;
  }
  for (auto& id : cols) {
    id -= col_offset;
  }

  // Create a grid based on the row gids I read in (may be duplicated across ranks)
  std::vector<gid_type> unique_gids;
  const auto& gids = type==InterpType::Refine ? rows : cols;
  for (auto gid : gids) {
    if (not ekat::contains(unique_gids,gid)) {
      unique_gids.push_back(gid);
    }
  }
  auto io_grid = std::make_shared<PointGrid> ("helper",unique_gids.size(),0,comm);
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
  GridImportExport imp_exp (fine_grid,io_grid);
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

void HorizRemapperData::
create_coarse_grids (const std::vector<Triplet>& triplets)
{
  // Gather overlapped coarse grid gids (rows or cols, depending on type)
  std::map<gid_type,int> ov_gid2lid;
  bool pickRow = type==InterpType::Coarsen;
  for (const auto& t : triplets) {
    ov_gid2lid.emplace(pickRow ? t.row : t.col,ov_gid2lid.size());
  }
  int num_ov_gids = ov_gid2lid.size();

  // Use a temp and then assing, b/c grid_ptr_type is a pointer to const,
  // so you can't modify gids using that pointer
  ov_coarse_grid = std::make_shared<PointGrid>("ov_coarse_grid",num_ov_gids,0,comm);
  auto ov_coarse_gids_h = ov_coarse_grid->get_dofs_gids().get_view<gid_type*,Host>();
  for (const auto& it : ov_gid2lid) {
    ov_coarse_gids_h[it.second] = it.first;
  }
  auto beg = ov_coarse_gids_h.data();
  auto end = beg+ov_coarse_gids_h.size();
  std::sort(beg,end);

  ov_coarse_grid->get_dofs_gids().sync_to_dev();

  // Create the unique coarse grid
  auto coarse_gids = ov_coarse_grid->get_unique_gids();
  int num_gids = coarse_gids.size();
  coarse_grid = std::make_shared<PointGrid>("coarse_grid",num_gids,0,comm);
  auto coarse_gids_h = coarse_grid->get_dofs_gids().get_view<gid_type*,Host>();
  std::copy(coarse_gids.begin(),coarse_gids.end(),coarse_gids_h.data());
  coarse_grid->get_dofs_gids().sync_to_dev();
}

void HorizRemapperData::
create_crs_matrix_structures (std::vector<Triplet>& triplets)
{
  // Get row/col data depending on interp type
  bool refine = type==InterpType::Refine;
  auto row_grid = refine ? fine_grid : ov_coarse_grid;
  auto col_grid = refine ? ov_coarse_grid : fine_grid;
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
  row_offsets = view_1d<int>("",num_rows+1);
  col_lids    = view_1d<int>("",nnz);
  weights     = view_1d<Real>("",nnz);

  auto row_offsets_h = Kokkos::create_mirror_view(row_offsets);
  auto col_lids_h    = Kokkos::create_mirror_view(col_lids);
  auto weights_h     = Kokkos::create_mirror_view(weights);

  // Fill col ids and weights
  for (int i=0; i<nnz; ++i) {
    col_lids_h(i) = col_gid2lid[triplets[i].col];
    weights_h(i)  = triplets[i].w;
  }
  Kokkos::deep_copy(weights,weights_h);
  Kokkos::deep_copy(col_lids,col_lids_h);

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

  Kokkos::deep_copy(row_offsets,row_offsets_h);
}

} // namespace scream
