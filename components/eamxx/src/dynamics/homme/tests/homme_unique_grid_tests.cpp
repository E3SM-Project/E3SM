#include <catch2/catch.hpp>

#include "dynamics/homme/interface/eamxx_homme_interface.hpp"
#include "share/remap/unique_grid_remapper.hpp"
#include "share/util/gid_directory.hpp"
#include "share/field/field.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/point_grid.hpp"

#include "SimulationParams.hpp"

#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/homme_grids_manager.hpp"

#include <ekat_comm.hpp>
#include <ekat_test_utils.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

extern "C" {
// These are specific C/F calls for these tests (i.e., not part of eamxx_homme_interface.hpp)
void init_test_params_f90 ();
void cleanup_test_f90 ();
}

namespace {

/*
 * Check that the dynamics grid and the physics GLL grid are related as expected.
 *
 * HommeGridsManager sets the physics GLL grid as the unique grid of the dynamics grid. This means that:
 *  - each dof of the dyn grid has a "unique gid", which is the gid of a column of the phys GLL grid;
 *  - each column of the phys GLL grid is the unique column of some dyn dofs (its replicas);
 *  - the replicas of a column are all at the same point of the sphere, which is also the one
 *    of the phys GLL column.
 * The UniqueGridRemapper relies on all of this (and it is tested in share/remap, with hand-made grids).
 */

TEST_CASE("homme_unique_grid", "") {
  using namespace scream;
  using gid_type   = AbstractGrid::gid_type;
  using Directory  = GidDirectory<gid_type>;
  using gid_view_h = typename Directory::gid_view_h;

  constexpr int NQ = HOMMEXX_QSIZE_D;
  constexpr int ne = 2;

  ekat::Comm comm(MPI_COMM_WORLD);

  // Init homme context
  if (!is_parallel_inited_f90()) {
    auto comm_f = MPI_Comm_c2f(MPI_COMM_WORLD);
    init_parallel_f90(comm_f);
  }
  init_test_params_f90 ();

  auto& c = Homme::Context::singleton();
  auto& sp = c.create<Homme::SimulationParams>();
  sp.qsize = NQ;
  set_homme_param("ne",ne);

  // Create the grids
  ekat::ParameterList params;
  params.set<std::string>("physics_grid_type","gll");
  params.set<std::string>("vertical_coordinate_filename","NONE");
  HommeGridsManager gm(comm,params);
  gm.build_grids();

  auto phys_grid = gm.get_grid("physics_gll");
  auto dyn_grid  = gm.get_grid("dynamics");
  REQUIRE (phys_grid->get_num_local_dofs()>0);

  // NOTE: we use CHECK (rather than REQUIRE) for the data checks, so that all ranks keep going in case
  //       of failure on a subset of ranks. Otherwise, the remaining ranks would hang in the collectives.

  // The unique grid of the dyn grid is the physics GLL grid
  CHECK (dyn_grid->get_unique_grid()->name()==phys_grid->name());
  CHECK (dyn_grid->get_num_vertical_levels()==phys_grid->get_num_vertical_levels());

  const int ncols = phys_grid->get_num_local_dofs();
  const int ndofs = dyn_grid->get_num_local_dofs();
  const auto p_gids = phys_grid->get_dofs_gids().get_view<const gid_type*,Host>();
  const auto d_gids = dyn_grid->get_unique_grid_gids().get_view<const gid_type*,Host>();
  const auto l2i    = dyn_grid->get_lid_to_idx_map().get_view<const int**,Host>();

  // --- The two sets of gids coincide (across all ranks) --- //
  {
    // Every dyn unique gid is the gid of a phys column
    Directory phys_dir (comm,gid_view_h(p_gids.data(),ncols));
    CHECK (not phys_dir.has_duplicates());
    const auto res = phys_dir.lookup(gid_view_h(d_gids.data(),ndofs));
    for (int idof=0; idof<ndofs; ++idof) {
      CHECK (res.pids[idof]>=0);
    }

    // Every phys gid is a unique gid of some dyn dofs (that may be owned by different ranks)
    std::vector<gid_type> distinct (d_gids.data(),d_gids.data()+ndofs);
    std::sort(distinct.begin(),distinct.end());
    distinct.erase(std::unique(distinct.begin(),distinct.end()),distinct.end());
    Directory dyn_dir (comm,gid_view_h(distinct.data(),distinct.size()));
    const auto res2 = dyn_dir.lookup(gid_view_h(p_gids.data(),ncols));
    for (int icol=0; icol<ncols; ++icol) {
      CHECK (res2.pids[icol]!=-1);
    }
  }

  // --- The points on the sphere coincide --- //
  {
    constexpr double pi = 3.14159265358979323846;
    // lat/lon in degrees
    auto xyz = [](const double lat_deg, const double lon_deg) {
      const double lat = lat_deg*pi/180, lon = lon_deg*pi/180;
      return std::array<double,3>{std::cos(lat)*std::cos(lon),std::cos(lat)*std::sin(lon),std::sin(lat)};
    };
    auto same_point = [](const std::array<double,3>& a, const std::array<double,3>& b) {
      return std::abs(a[0]-b[0])<1e-4 and std::abs(a[1]-b[1])<1e-4 and std::abs(a[2]-b[2])<1e-4;
    };

    // NOTE: the dyn grid data is in degrees too (despite the grid labeling it as radians)
    const auto d_lat = dyn_grid->get_geometry_data("lat").get_view<const Real***,Host>();
    const auto d_lon = dyn_grid->get_geometry_data("lon").get_view<const Real***,Host>();
    const auto p_lat = phys_grid->get_geometry_data("lat").get_view<const Real*,Host>();
    const auto p_lon = phys_grid->get_geometry_data("lon").get_view<const Real*,Host>();
    auto dyn_xyz = [&](const int idof) {
      return xyz(d_lat(l2i(idof,0),l2i(idof,1),l2i(idof,2)),d_lon(l2i(idof,0),l2i(idof,1),l2i(idof,2)));
    };

    // Each phys column has a replica on the same rank (the unique grid remapper requires it), and it is at the same point
    auto map = compute_unique_grid_dof_map(*dyn_grid,*phys_grid);
    auto map_h = Kokkos::create_mirror_view(map);
    Kokkos::deep_copy(map_h,map);
    for (int icol=0; icol<ncols; ++icol) {
      const int idof = map_h(icol);
      CHECK (d_gids(idof)==p_gids(icol));
      CHECK (same_point(xyz(p_lat(icol),p_lon(icol)),dyn_xyz(idof)));
    }

    // The replicas of a column on this rank are all at the same point
    std::map<gid_type,int> first_dof;
    for (int idof=0; idof<ndofs; ++idof) {
      const auto it = first_dof.emplace(d_gids(idof),idof).first;
      CHECK (same_point(dyn_xyz(idof),dyn_xyz(it->second)));
    }
  }

  Homme::Context::finalize_singleton();
  cleanup_test_f90();
}

} // anonymous namespace
