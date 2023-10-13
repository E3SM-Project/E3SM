#include "horiz_interp_remapper_base.hpp"

#include "share/grid/point_grid.hpp"
#include "share/grid/grid_import_export.hpp"
#include "share/io/scorpio_input.hpp"

#include <numeric>

namespace scream
{

auto HorizInterpRemapperBase::
get_my_triplets (const std::string& map_file,
                 const ekat::Comm&  comm,
                 const std::shared_ptr<const AbstractGrid>& grid,
                 const OwnedBy owned_by) const
 -> std::vector<Triplet>
{
  using gid_type = AbstractGrid::gid_type;
  using namespace ShortFieldTagsNames;

  // 1. Load the map file chunking it evenly across all ranks
  scorpio::register_file(map_file,scorpio::FileMode::Read);

  // 1.1 Create a "helper" grid, with as many dofs as the number
  //     of triplets in the map file, and divided linearly across ranks
  const int ngweights = scorpio::get_dimlen(map_file,"n_s");
  int nlweights = ngweights / comm.size();
  if (comm.rank() < (ngweights % comm.size())) {
    nlweights += 1;
  }

  gid_type offset = nlweights;
  comm.scan(&offset,1,MPI_SUM);
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
  comm.all_reduce(&map_file_min_row,&global_map_file_min_row,1,MPI_MIN);
  comm.all_reduce(&map_file_min_col,&global_map_file_min_col,1,MPI_MIN);

  gid_type row_offset = global_map_file_min_row;
  gid_type col_offset = global_map_file_min_col;
  if (owned_by==OwnedBy::Row) {
    row_offset -= grid->get_global_min_dof_gid();
  } else {
    col_offset -= grid->get_global_min_dof_gid();
  }
  for (auto& id : rows) {
    id -= row_offset;
  }
  for (auto& id : cols) {
    id -= col_offset;
  }

  // Create a grid based on the row gids I read in (may be duplicated across ranks)
  std::vector<gid_type> unique_gids;
  const auto& gids = owned_by==OwnedBy::Row ? rows : cols;
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
  GridImportExport imp_exp (grid,io_grid);
  std::map<int,std::vector<Triplet>> my_triplets_map;
  imp_exp.gather(mpi_triplet_t,io_triplets,my_triplets_map);

  std::vector<Triplet> my_triplets;
  for (auto& it : my_triplets_map) {
    my_triplets.reserve(my_triplets.size()+it.second.size());
    std::move(it.second.begin(),it.second.end(),std::back_inserter(my_triplets));
  }

  return my_triplets;
}

} // namespace scream
