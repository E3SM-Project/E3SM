#ifndef SCREAM_MESH_FREE_GRIDS_MANAGER_HPP
#define SCREAM_MESH_FREE_GRIDS_MANAGER_HPP

#include "share/grid/grids_manager.hpp"

#include "ekat/mpi/ekat_comm.hpp"

namespace scream {

// This class is meant to be used for small unit tests, where we want to
// test some infrastructure of scream, without bothering too much about
// grids-related features. This manager can build a PointGrid and a SEGrid.
// If the SEGrid is build, it corresponds to a "strip" of elements,
// and the PointGrid is built as the 'unique' grid of that SEGrid.

class MeshFreeGridsManager : public GridsManager
{
public:
  using string_pair = std::pair<std::string,std::string>;
  using remap_repo_type = std::map<string_pair,remapper_ptr_type>;

  MeshFreeGridsManager (const ekat::Comm& comm, const ekat::ParameterList& p);

  virtual ~MeshFreeGridsManager () = default;

  std::string name () const { return "Mesh-Free Grids Manager"; }

  std::set<std::string> supported_grids () const {
    std::set<std::string> gnames;

    gnames.insert("Point Grid");
    gnames.insert("SE Grid");

    return gnames;
  }

  void build_grids (const std::set<std::string>& grid_names);

  const grid_repo_type& get_repo () const { return m_grids; }

protected:

  void build_pt_from_se (const int num_local_elems, const int num_gp, const int num_vertical_levels);

  std::string get_reference_grid_name () const {
    return m_params.get<std::string>("Reference Grid");
  }

  remapper_ptr_type
  do_create_remapper (const grid_ptr_type from_grid,
                      const grid_ptr_type to_grid) const;

  grid_repo_type      m_grids;

  remap_repo_type     m_remappers;

  ekat::ParameterList m_params;

  ekat::Comm          m_comm;
};

inline std::shared_ptr<GridsManager>
create_mesh_free_grids_manager (const ekat::Comm& comm, const ekat::ParameterList& p) {
  return std::make_shared<MeshFreeGridsManager>(comm,p);
}

std::shared_ptr<GridsManager>
create_mesh_free_grids_manager (const ekat::Comm& comm, const int num_local_elems,
                                const int num_gp, const int num_vertical_levels);

inline void register_mesh_free_grids_manager () {
  // A physics-only grids manager, in case we run a physics-only test
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("Mesh Free",&create_mesh_free_grids_manager);
}

} // namespace scream

#endif // SCREAM_MESH_FREE_GRIDS_MANAGER_HPP
