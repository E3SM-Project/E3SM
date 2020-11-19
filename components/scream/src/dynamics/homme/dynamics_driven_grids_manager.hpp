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

  DynamicsDrivenGridsManager (const ekat::Comm& comm, const ekat::ParameterList& p);

  ~DynamicsDrivenGridsManager ();

  std::string name () const { return "Dynamics Driven Grids Manager"; }

  void build_grids (const std::set<std::string>& grid_names,
                    const std::string& reference_grid);

  std::set<std::string> supported_grids () const { return m_valid_grid_names; }
protected:

  remapper_ptr_type
  do_create_remapper (const grid_ptr_type from_grid,
                      const grid_ptr_type to_grid) const;

  void build_dynamics_grid ();
  void build_physics_grid  (const std::string& name);

  const grid_repo_type& get_repo () const { return m_grids; }

  void build_grid_codes ();

  grid_repo_type  m_grids;

  // Admissible grid names
  std::set<std::string> m_valid_grid_names;

  // For each admissible grid name, store an integer code
  std::map<std::string, int> m_grid_codes;
};

inline std::shared_ptr<GridsManager>
create_dynamics_driven_grids_manager (const ekat::Comm& comm, const ekat::ParameterList& p) {
  return std::make_shared<DynamicsDrivenGridsManager>(comm,p);
}

} // namespace scream

#endif // SCREAM_DYNAMICS_DRIVEN_GRIDS_MANAGER_HPP
