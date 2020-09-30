#ifndef SCREAM_GRIDS_MANAGER_HPP
#define SCREAM_GRIDS_MANAGER_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/grid/remap/abstract_remapper.hpp"
#include "share/grid/remap/identity_remapper.hpp"

#include "ekat/util/ekat_factory.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/mpi/ekat_comm.hpp"

#include <map>
#include <set>
#include <memory>

namespace scream
{

class GridsManager
{
public:
  using grid_type         = AbstractGrid;
  using device_type       = grid_type::device_type;
  using grid_ptr_type     = std::shared_ptr<grid_type>;
  using grid_repo_type    = std::map<std::string, grid_ptr_type>;
  using remapper_type     = AbstractRemapper<Real,device_type>;
  using remapper_ptr_type = std::shared_ptr<remapper_type>;

  GridsManager () = default;
  virtual ~GridsManager () = default;

  virtual std::string name () const = 0;

  grid_ptr_type get_grid (const std::string& name) const;

  grid_ptr_type get_reference_grid () const {
    return get_grid("Reference");
  }

  void build_grids (const std::set<std::string>& grid_names,
                    const std::string& reference_grid) {
    for (const auto& name : grid_names) {
      build_grid(name);
    }
    build_grid(reference_grid);
  }

  remapper_ptr_type
  create_remapper (const grid_ptr_type from_grid,
                   const grid_ptr_type to_grid) const {
    EKAT_REQUIRE_MSG( supports_grid(from_grid->name()),
                        "Error! Source grid '" + from_grid->name() + "' is not supported.\n");
    EKAT_REQUIRE_MSG( supports_grid(to_grid->name()),
                        "Error! Target grid '" + to_grid->name() + "' is not supported.\n");

    remapper_ptr_type remapper;

    if (from_grid->name()==to_grid->name()) {
      // We can handle the identity remapper from here
      remapper = std::make_shared<IdentityRemapper<Real,device_type>>(from_grid);
    } else {
      remapper = do_create_remapper(from_grid,to_grid);
    }

    EKAT_REQUIRE_MSG(
      remapper!=nullptr,
      "Error! A remapper from grid '" + from_grid->name() + "' to grid '" + to_grid->name() + "' is not available.\n"
      "       Perhaps you forgot to add it creation to the implementation of the grids manager?\n");

    return remapper;
  }

protected:

  virtual void build_grid (const std::string& grid_name) = 0;

  virtual remapper_ptr_type
  do_create_remapper (const grid_ptr_type from_grid,
                      const grid_ptr_type to_grid) const = 0;

  bool supports_grid (const std::string& grid_name) const {
    const auto& grids = get_repo ();
    return grids.find(grid_name)!=grids.end();
  }

  // This mini-function simply prints the supported grids names as "name1, name2, name3"
  std::string print_supported_grids () const {
    const auto& grids = get_repo ();

    std::string str;
    if (grids.size()==0) {
      return str;
    }
    auto next = ++grids.begin();
    for (auto it=grids.begin(); next!=grids.end(); ++it, ++next) {
      str += it->first + ", ";
    }
    str += next->first;
    return str;
  }

  virtual const grid_repo_type& get_repo () const = 0;
};

inline GridsManager::grid_ptr_type
GridsManager::get_grid(const std::string& name) const
{
  EKAT_REQUIRE_MSG (supports_grid(name),
                      "Error! Grids manager '" + this->name() + "' does not provide grid '" + name + "'.\n"
                      "       Supported grids are: " + print_supported_grids()  + "\n");

  return get_repo().at(name);
}

// A short name for the factory for grid managers
using GridsManagerFactory 
    = ekat::Factory<GridsManager,
                    ekat::CaseInsensitiveString,
                    std::shared_ptr<GridsManager>,
                    const ekat::Comm&,const ekat::ParameterList&>;

} // namespace scream

#endif // SCREAM_GRIDS_MANAGER_HPP
