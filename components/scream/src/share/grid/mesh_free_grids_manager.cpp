#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/remap/do_nothing_remapper.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"

#include <numeric>

namespace scream {

MeshFreeGridsManager::
MeshFreeGridsManager (const ekat::Comm& comm, const ekat::ParameterList& p)
 : m_params (p)
 , m_comm   (comm)
{
  const auto ref_grid = m_params.get<std::string>("Reference Grid","Point Grid");
  EKAT_REQUIRE_MSG (ekat::contains(supported_grids(),ref_grid),
        "Error! MeshFreeGridsManager only supports 'Point Grid' and 'SE Grid' grids.\n"
        "       Requested reference grid: " + ref_grid + "\n");
}

MeshFreeGridsManager::remapper_ptr_type
MeshFreeGridsManager::
do_create_remapper (const grid_ptr_type from_grid,
                    const grid_ptr_type to_grid) const
{
  return std::make_shared<DoNothingRemapper<Real> >(from_grid,to_grid);
}

void MeshFreeGridsManager::
build_grids (const std::set<std::string>& grid_names)
{
  for (const auto& gn : grid_names) {
    EKAT_REQUIRE_MSG (ekat::contains(supported_grids(),gn),
        "Error! MeshFreeGridsManager only supports 'Point Grid' and 'SE Grid' grids.\n"
        "       Requested grid: " + gn + "\n");
  }

  std::string ref_grid = m_params.get<std::string>("Reference Grid");

  const bool build_pt = ekat::contains(grid_names,"Point Grid") || ref_grid=="Point Grid";
  const bool build_se = ekat::contains(grid_names,"SE Grid")    || ref_grid=="SE Grid";

  const auto& gm_params         = m_params.sublist("Mesh Free");
  const int num_vertical_levels = gm_params.get<int>("Number of Vertical Levels");

  if (build_se) {
    // Assume a "strip" of elements:
    //
    //  *---*---*     *---*
    //  |   |   | ... |   *
    //  *---*---*     *---*
    //
    // partitioned evenly across ranks. Use this config to establish
    // a unique numbering for each element, using the following (np=4 here):
    //
    //   1--5--9-13  13-17-21-25
    //   |  |  |  |   |  |  |  |
    //   2--6-10-14  14-18-22-26
    //   |  |  |  |   |  |  |  |
    //   3--7-11-15  15-19-23-27
    //   |  |  |  |   |  |  |  |
    //   4--8-12-16  16-20-24-28

    const int num_local_elems  = gm_params.get<int>("Number of Local Elements");
    const int num_gp           = gm_params.get<int>("Number of Gauss Points");

    // Build the PointGrid first, so we can set it as the 'unique' grid of the SEGrid
    build_pt_from_se (num_local_elems,num_gp,num_vertical_levels);

    // Set up the degrees of freedom.
    SEGrid::dofs_list_type dofs("", num_local_elems*num_gp*num_gp);
    SEGrid::lid_to_idx_map_type dofs_map("", num_local_elems*num_gp*num_gp, 3);

    auto host_dofs = Kokkos::create_mirror_view(dofs);
    auto host_dofs_map = Kokkos::create_mirror_view(dofs_map);

    // Count unique local dofs. On all elems except the very last one (on rank N),
    // we have num_gp*(num_gp-1) unique dofs;
    int num_elem_unique_dofs = num_gp*(num_gp-1);
    int num_local_unique_dofs = num_local_elems*num_elem_unique_dofs;
    int offset = num_local_unique_dofs*m_comm.rank();

    for (int ie = 0; ie < num_local_elems; ++ie) {
      for (int igp = 0; igp < num_gp-1; ++igp) {
        for (int jgp = 0; jgp < num_gp; ++jgp) {
          int idof = ie*num_gp*num_gp + igp*num_gp + jgp;
          int gid = offset + ie*num_elem_unique_dofs + igp*num_gp + jgp;
          host_dofs(idof) = gid;
          host_dofs_map(idof, 0) = ie;
          host_dofs_map(idof, 1) = igp;
          host_dofs_map(idof, 2) = jgp;
        }
      }
      for (int jgp = 0; jgp < num_gp; ++jgp) {
        int idof = ie*num_gp*num_gp + (num_gp-1)*num_gp + jgp;
        int gid = offset + (ie+1)*num_elem_unique_dofs + jgp;
        host_dofs(idof) = gid;
        host_dofs_map(idof, 0) = ie;
        host_dofs_map(idof, 1) = num_gp-1;
        host_dofs_map(idof, 2) = jgp;
      }
    }

    // Move the data to the device and set the DOFs.
    Kokkos::deep_copy(dofs, host_dofs);
    Kokkos::deep_copy(dofs_map, host_dofs_map);

    // Create the grid, and set the dofs
    auto se_grid = std::make_shared<SEGrid>("SE Grid",num_local_elems,num_gp,num_vertical_levels,m_grids["Point Grid"],m_comm);
    se_grid->setSelfPointer(se_grid);

    se_grid->set_dofs(dofs);
    se_grid->set_lid_to_idx_map(dofs_map);

    m_grids["SE Grid"]    = se_grid;
  } else if (build_pt) {
    const int num_global_cols  = gm_params.get<int>("Number of Global Columns");
    auto pt_grid = create_point_grid("Point Grid",num_global_cols,num_vertical_levels,m_comm);
    m_grids["Point Grid"] = pt_grid;
  }
}

void MeshFreeGridsManager::
build_pt_from_se (const int num_local_elems, const int num_gp, const int num_vertical_levels)
{
  // Assume a "strip" of elements:
  //
  //  *---*---*     *---*
  //  |   |   | ... |   *
  //  *---*---*     *---*
  //
  // partitioned evenly across ranks. Use this config to establish
  // a unique numbering for each element, using the following (np=4 here):
  //
  //   1--5--9-13  13-17-21-25
  //   |  |  |  |   |  |  |  |
  //   2--6-10-14  14-18-22-26
  //   |  |  |  |   |  |  |  |
  //   3--7-11-15  15-19-23-27
  //   |  |  |  |   |  |  |  |
  //   4--8-12-16  16-20-24-28
  //
  // Shared dofs are "assigned" to the element on the right, so the PointGrid
  // associated to this SEGrid has np*(np-1)*num_elems dofs, except for the last
  // rank, which has an extra np dofs.

  const int extra_dofs    = m_comm.rank()==(m_comm.size()-1) ? num_gp : 0;
  const int num_cols_base = num_local_elems*num_gp*(num_gp-1);
  const int num_my_cols   = num_cols_base + extra_dofs;
  const int dofs_offset   = m_comm.rank()*num_cols_base;

  PointGrid::dofs_list_type dofs("",num_my_cols);
  auto dofs_h = Kokkos::create_mirror_view(dofs);
  std::iota(dofs_h.data(),dofs_h.data()+num_my_cols,dofs_offset);
  Kokkos::deep_copy(dofs,dofs_h);

  auto pt_grid = std::make_shared<PointGrid>("Point Grid",num_my_cols,num_vertical_levels,m_comm);
  pt_grid->setSelfPointer(pt_grid);
  pt_grid->set_dofs(dofs);

  m_grids["Point Grid"] = pt_grid;
}

std::shared_ptr<GridsManager>
create_mesh_free_grids_manager (const ekat::Comm& comm, const int num_local_elems,
                                const int num_gp, const int num_vertical_levels)
{
  ekat::ParameterList gm_params;
  gm_params.sublist("Mesh Free").set("Number of Local Elements",num_local_elems);
  gm_params.sublist("Mesh Free").set("Number of Gauss Points",num_gp);
  gm_params.sublist("Mesh Free").set("Number of Vertical Levels",num_vertical_levels);
  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_all_grids();
  return gm;
}

} // namespace scream
