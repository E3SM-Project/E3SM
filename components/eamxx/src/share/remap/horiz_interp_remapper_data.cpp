#include "horiz_interp_remapper_data.hpp"

#include "share/grid/point_grid.hpp"
#include "share/grid/grid_import_export.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <numeric>

namespace scream {

void HorizRemapperData::
build (const std::shared_ptr<const AbstractGrid>& grid,
       const std::string& map_file)
{
  m_input_grid = grid;

  EKAT_REQUIRE_MSG (grid,
      "[HorizRemapperDataRepo::build_data_from_src] Error! Invalid src grid pointer.\n");

  int ncol_a = scorpio::get_dimlen(map_file,"n_a");
  int ncol_b = scorpio::get_dimlen(map_file,"n_b");

  EKAT_REQUIRE_MSG (ncol_a!=ncol_b,
    "[HorizRemapperDataRepo] Error! Source and target grid in the map file MUST have a DIFFERENT number of columns.\n"
    " - map file: " + map_file + "\n"
    " - n_a: " + std::to_string(ncol_a) + "\n"
    " - n_b: " + std::to_string(ncol_b) + "\n"
    "If this is a limiting factor for you, please, contact developers to see if we can relax this assumption.\n");

  // Figure out which direction the remap is going (to or from the input grid)
  int grid_ncol = grid->get_num_global_dofs();
  EKAT_REQUIRE_MSG (grid_ncol==ncol_a or grid_ncol==ncol_b,
    "Error! The number of cols on the input grid does not match either of the map file 'n_a' or 'n_b' dims.\n"
    " - map file: " + map_file + "\n"
    " - grid name: " + grid->name() +"\n"
    " - grid ncol: " + std::to_string(grid_ncol) + "\n"
    " - n_a: " + std::to_string(ncol_a) + "\n"
    " - n_b: " + std::to_string(ncol_b) + "\n");

  m_coarsening = ncol_a>=ncol_b;
  m_built_from_src   = grid_ncol==ncol_a;

  const int nlev = grid->get_num_vertical_levels();
  const auto& comm = grid->get_comm();
  const auto nondim = ekat::units::Units::nondimensional();
  ekat::units::Units deg(nondim,"deg");
  std::string suffix = m_built_from_src ? "_b" : "_a";

  m_generated_grid = create_point_grid(m_built_from_src ? "tgt_grid" : "src_grid",m_built_from_src ? ncol_b : ncol_a,nlev,comm,1);

  // Only read the lat/lon/area vars if they are present. If one is present, we assume they all are
  if (scorpio::has_var(map_file,"yc"+suffix)) {
    m_generated_grid->create_geometry_data("lat", m_generated_grid->get_2d_scalar_layout(),deg);
    m_generated_grid->create_geometry_data("lon", m_generated_grid->get_2d_scalar_layout(),deg);
    m_generated_grid->create_geometry_data("area",m_generated_grid->get_2d_scalar_layout(),nondim);
    m_generated_grid->read_geometry_data(map_file,
                                  {"lat","lon","area"},
                                  {"yc"+suffix,"xc"+suffix,"area"+suffix},
                                  {{"ncol","n"+suffix}});
  }

  // Load sparse matrix triplets, splitting evenly across ranks
  auto triplets = read_mat_triplets(map_file);

  // Gather sparse matrix triplets needed by this rank
  auto my_triplets = get_my_triplets (triplets);

  // Create aux and ov grids
  create_ov_grid (my_triplets);

  // Create crs matrix
  create_crs_matrix_structures (my_triplets);
}

std::vector<Triplet>
HorizRemapperData::
read_mat_triplets (const std::string& map_file)
{
  using gid_type = AbstractGrid::gid_type;
  using namespace ShortFieldTagsNames;

  const auto& comm = m_input_grid->get_comm();

  // Split the triplets evenly across ranks, and read them
  int n_s = scorpio::get_dimlen(map_file,"n_s");
  auto io_grid = create_point_grid("remap_data_io_grid",n_s,0,comm);
  io_grid->reset_field_tag_name(COL,"n_s");

  auto row_h = io_grid->create_geometry_data<gid_type>("row",io_grid->get_2d_scalar_layout()).get_view<const gid_type*,Host>();
  auto col_h = io_grid->create_geometry_data<gid_type>("col",io_grid->get_2d_scalar_layout()).get_view<const gid_type*,Host>();
  auto S_h   = io_grid->create_geometry_data<Real>("S",  io_grid->get_2d_scalar_layout()).get_view<const Real*,Host>();
  io_grid->read_geometry_data(map_file,{"row","col","S"});

  int nlweights = io_grid->get_num_local_dofs();

  std::vector<Triplet> triplets;
  for (int i=0; i<nlweights; ++i) {
    triplets.emplace_back(row_h[i],col_h[i],S_h[i]);
  }

  return triplets;
}

std::vector<Triplet>
HorizRemapperData::
get_my_triplets (const std::vector<Triplet>& triplets)
{
  using gid_type = AbstractGrid::gid_type;

  // Create a grid where the GIDs are the id of rows or cols of the triplets we read
  // We pick row/col based on which side of the remap was the input grid
  std::set<gid_type> unique_gids;
  for (const auto& t : triplets) {
    unique_gids.insert(m_coarsening ? t.col : t.row);
  }
  auto io_grid = std::make_shared<PointGrid> ("helper",unique_gids.size(),0,m_input_grid->get_comm());
  auto io_grid_gids_h = io_grid->get_dofs_gids().get_view<gid_type*,Host>();
  int k = 0;
  for (auto gid : unique_gids) {
    io_grid_gids_h(k++) = gid;
  }
  io_grid->get_dofs_gids().sync_to_dev();

  // Group triplets to export by their gid
  std::map<int,std::vector<Triplet>> io_triplets;
  const auto& io_grid_gid2lid = io_grid->get_gid2lid_map();
  for (const auto& t : triplets) {
    auto io_lid = io_grid_gid2lid.at(m_coarsening ? t.col : t.row);
    io_triplets[io_lid].emplace_back(t.row,t.col,t.w);
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

  // Create import-export and gather the triplets we need for our local mat-vec
  auto fine_grid = m_coarsening ? (m_built_from_src ? m_input_grid : m_generated_grid)
                                : (m_built_from_src ? m_generated_grid : m_input_grid);
  std::map<int,std::vector<Triplet>> my_triplets_map;
  GridImportExport imp_exp (fine_grid,io_grid);
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
create_ov_grid (const std::vector<Triplet>& my_triplets)
{
  using gid_type = AbstractGrid::gid_type;

  // Gather overlapped coarse grid gids (rows or cols, depending on refine vs m_coarsening)
  std::map<gid_type,int> ov_gid2lid;
  for (const auto& t : my_triplets) {
    ov_gid2lid.emplace(m_coarsening ? t.row : t.col, ov_gid2lid.size());
  }
  int num_ov_gids = ov_gid2lid.size();

  m_overlap_grid = std::make_shared<PointGrid>("ov_coarse_grid",num_ov_gids,0,m_input_grid->get_comm());
  auto gids_h = m_overlap_grid->get_dofs_gids().get_view<gid_type*,Host>();
  for (const auto& it : ov_gid2lid) {
    gids_h[it.second] = it.first;
  }
  auto beg = gids_h.data();
  auto end = beg+gids_h.size();
  std::sort(beg,end);
  m_overlap_grid->get_dofs_gids().sync_to_dev();
}

void HorizRemapperData::
create_crs_matrix_structures (std::vector<Triplet>& triplets)
{
  auto fine_grid = m_coarsening ? (m_built_from_src ? m_input_grid : m_generated_grid)
                                : (m_built_from_src ? m_generated_grid : m_input_grid);

  auto src_grid = m_coarsening ? fine_grid : m_overlap_grid;
  auto tgt_grid = m_coarsening ? m_overlap_grid : fine_grid;

  // Get row/col data depending on interp type
  const int num_rows = tgt_grid->get_num_local_dofs();

  const auto& col_gid2lid = src_grid->get_gid2lid_map();
  const auto& row_gid2lid = tgt_grid->get_gid2lid_map();

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
    col_lids_h(i) = col_gid2lid.at(triplets[i].col);
    weights_h(i)  = triplets[i].w;
  }
  Kokkos::deep_copy(m_weights,weights_h);
  Kokkos::deep_copy(m_col_lids,col_lids_h);

  // Compute row offsets
  std::vector<int> row_counts(num_rows);
  for (int i=0; i<nnz; ++i) {
    ++row_counts[row_gid2lid.at(triplets[i].row)];
  }
  std::partial_sum(row_counts.begin(),row_counts.end(),row_offsets_h.data()+1);
  EKAT_REQUIRE_MSG (
      row_offsets_h(num_rows)==nnz,
      "Error! Something went wrong while computing row offsets.\n"
      "  - local nnz       : " + std::to_string(nnz) + "\n"
      "  - row_offsets(end): " + std::to_string(row_offsets_h(num_rows)) + "\n");

  Kokkos::deep_copy(m_row_offsets,row_offsets_h);
}

std::shared_ptr<const HorizRemapperData>
HorizRemapperDataRepo::
get_data (const std::shared_ptr<const AbstractGrid>& grid,
          const std::string& map_file)
{
  auto& data = m_repo[map_file];
  if (auto shared_data = data.lock()) {
    return shared_data;
  }
  
  // Either there was no data for this map file, or the existing weak_ptr was expired.
  // E.g., there WAS a remapper that used this data, but the remapper has since been
  // destroyed. Either way, we can safely (re-)create the data

  auto shared_data = std::make_shared<HorizRemapperData>();
  shared_data->build(grid,map_file);
  data = shared_data;

  return shared_data;
}

} // namespace scream
