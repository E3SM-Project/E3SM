#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/remap/do_nothing_remapper.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"

#include <memory>
#include <numeric>

namespace scream {

MeshFreeGridsManager::
MeshFreeGridsManager (const ekat::Comm& comm, const ekat::ParameterList& p)
 : m_params (p)
 , m_comm   (comm)
{
}

MeshFreeGridsManager::remapper_ptr_type
MeshFreeGridsManager::
do_create_remapper (const grid_ptr_type from_grid,
                    const grid_ptr_type to_grid) const
{
  return std::make_shared<DoNothingRemapper>(from_grid,to_grid);
}

void MeshFreeGridsManager::
build_grids ()
{
  auto has_positive_int = [&](const std::string& n) -> bool {
    return m_params.isParameter(n) && (m_params.get<int>(n)>0);
  };
  const bool build_pt = has_positive_int("number_of_global_columns");
  const bool build_se = has_positive_int("number_of_local_elements") &&
                        has_positive_int("number_of_gauss_points");

  const int num_vertical_levels = m_params.get<int>("number_of_vertical_levels");

  if (build_se) {
    // Build a set of completely disconnected spectral elements.
    const int num_local_elems  = m_params.get<int>("number_of_local_elements");
    const int num_gp           = m_params.get<int>("number_of_gauss_points");

    // Set up the degrees of freedom.
    SEGrid::dofs_list_type dofs("", num_local_elems*num_gp*num_gp);
    SEGrid::lid_to_idx_map_type dofs_map("", num_local_elems*num_gp*num_gp, 3);

    auto host_dofs = Kokkos::create_mirror_view(dofs);
    auto host_dofs_map = Kokkos::create_mirror_view(dofs_map);

    // Count unique local dofs. On all elems except the very last one (on rank N),
    // we have num_gp*(num_gp-1) unique dofs;
    int num_local_dofs = num_local_elems*num_gp*num_gp;
    int offset = num_local_dofs*m_comm.rank();

    for (int ie = 0; ie < num_local_elems; ++ie) {
      for (int igp = 0; igp < num_gp; ++igp) {
        for (int jgp = 0; jgp < num_gp; ++jgp) {
          int idof = ie*num_gp*num_gp + igp*num_gp + jgp;
          int gid = offset + idof;
          host_dofs(idof) = gid;
          host_dofs_map(idof, 0) = ie;
          host_dofs_map(idof, 1) = igp;
          host_dofs_map(idof, 2) = jgp;
        }
      }
    }

    // Move the data to the device and set the DOFs.
    Kokkos::deep_copy(dofs, host_dofs);
    Kokkos::deep_copy(dofs_map, host_dofs_map);

    // Create the grid, and set the dofs
    std::shared_ptr<SEGrid> se_grid;
    se_grid = std::make_shared<SEGrid>("SE Grid",num_local_elems,num_gp,num_vertical_levels,m_comm);
    se_grid->setSelfPointer(se_grid);

    se_grid->set_dofs(dofs);
    se_grid->set_lid_to_idx_map(dofs_map);

    add_grid(se_grid);
  }
  if (build_pt) {
    const int num_global_cols  = m_params.get<int>("number_of_global_columns");
    auto pt_grid = create_point_grid("Point Grid",num_global_cols,num_vertical_levels,m_comm);
    add_grid(pt_grid);
    this->alias_grid("Point Grid", "Physics");
  }
}

std::shared_ptr<GridsManager>
create_mesh_free_grids_manager (const ekat::Comm& comm, const int num_local_elems,
                                const int num_gp, const int num_vertical_levels,
                                const int num_global_cols)
{
  ekat::ParameterList gm_params;
  gm_params.set("number_of_global_columns",num_global_cols);
  gm_params.set("number_of_local_elements",num_local_elems);
  gm_params.set("number_of_gauss_points",num_gp);
  gm_params.set("number_of_vertical_levels",num_vertical_levels);
  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  return gm;
}

} // namespace scream
