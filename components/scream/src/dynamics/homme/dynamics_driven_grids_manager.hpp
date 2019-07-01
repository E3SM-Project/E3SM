#ifndef SCREAM_DYNAMICS_DRIVEN_GRIDS_MANAGER_HPP
#define SCREAM_DYNAMICS_DRIVEN_GRIDS_MANAGER_HPP

#include "share/grid/grids_manager.hpp"

namespace scream
{

// This class is meant to be used for small unit tests, where we want to
// test the Atmosphere Driver (AD) capabilities, without bothering too much
// about grids-related features. This manager lets you set pre-built grids
// in it rather than building them inside the manager.
class DynamicsDrivenGridsManager : public GridsManager
{
public:

  DynamicsDrivenGridsManager (const Comm& comm, const ParameterList& p);

  ~DynamicsDrivenGridsManager () = default;

  std::string name () const { return "Dynamics Driven Grids Manager"; }

protected:

  void build_grid (const std::string& grid_names);

  void build_dynamics_grid ();
  void build_physics_grid  ();

  const repo_type& get_repo () const { return m_grids; }

  repo_type  m_grids;

  AbstractGrid::dofs_map_type m_phys_dofs;
  AbstractGrid::dofs_map_type m_dyn_dofs;
};

inline std::shared_ptr<GridsManager>
create_dynamics_driven_grids_manager (const Comm& comm, const ParameterList& p) {
  return std::make_shared<DynamicsDrivenGridsManager>(comm,p);
}

} // namespace scream

#endif // SCREAM_DYNAMICS_DRIVEN_GRIDS_MANAGER_HPP

