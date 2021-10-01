#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/grid/remap/do_nothing_remapper.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"

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

  const auto& gm_params      = m_params.sublist("Mesh Free");
  const int num_vertical_lev = gm_params.get<int>("Number of Vertical Levels");

  if (build_pt) {
    const int num_global_cols  = gm_params.get<int>("Number of Global Columns");
    auto pt_grid = create_point_grid("Point Grid",num_global_cols,num_vertical_lev,m_comm);
    m_grids["Point Grid"] = pt_grid;
  }

  if (build_se) {
    const int num_local_elems  = gm_params.get<int>("Number of Local Elements");
    const int num_gp           = gm_params.get<int>("Number of Gauss Points");
    auto se_grid = std::make_shared<SEGrid>("SE Grid",num_local_elems,num_gp,num_vertical_lev,m_comm);
    m_grids["SE Grid"]    = se_grid;
  }
}

} // namespace scream
