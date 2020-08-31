#ifndef SCREAM_PHYSICS_ONLY_GRIDS_MANAGER_HPP
#define SCREAM_PHYSICS_ONLY_GRIDS_MANAGER_HPP

#include "share/grid/grids_manager.hpp"

#include "ekat/mpi/ekat_comm.hpp"

namespace scream
{

// This class is meant to be used for small unit tests, where we want to
// test the Atmosphere Driver (AD) capabilities, without bothering too much
// about grids-related features. This manager lets you set pre-built grids
// in it rather than building them inside the manager.
class PhysicsOnlyGridsManager : public GridsManager
{
public:
  using string_pair = std::pair<std::string,std::string>;
  using remap_repo_type = std::map<string_pair,remapper_ptr_type>;

  PhysicsOnlyGridsManager (const ekat::Comm& comm, const ekat::ParameterList& p);

  virtual ~PhysicsOnlyGridsManager () = default;

  std::string name () const { return "Physics-Only Grids Manager"; }

protected:

  remapper_ptr_type
  do_create_remapper (const grid_ptr_type from_grid,
                      const grid_ptr_type to_grid) const;

  void build_grid (const std::string& grid_name);

  const grid_repo_type& get_repo () const { return m_grids; }

  grid_repo_type  m_grids;

  remap_repo_type m_remappers;

  ekat::ParameterList m_params;

  ekat::Comm          m_comm;
};

inline std::shared_ptr<GridsManager>
create_physics_only_grids_manager (const ekat::Comm& comm, const ekat::ParameterList& p) {
  return std::make_shared<PhysicsOnlyGridsManager>(comm,p);
}

} // namespace scream

#endif // SCREAM_PHYSICS_ONLY_GRIDS_MANAGER_HPP
