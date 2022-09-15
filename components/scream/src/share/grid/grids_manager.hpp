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
  using grid_ptr_type     = std::shared_ptr<const grid_type>;
  using grid_repo_type    = std::map<std::string, grid_ptr_type>;
  using remapper_type     = AbstractRemapper;
  using remapper_ptr_type = std::shared_ptr<remapper_type>;

  GridsManager () = default;
  virtual ~GridsManager () = default;

  virtual std::string name () const = 0;

  grid_ptr_type get_grid (const std::string& name) const;

  // Check if the given grid has been built
  bool has_grid (const std::string& grid_name) const;

  virtual void build_grids () = 0;

  remapper_ptr_type
  create_remapper (const grid_ptr_type& from_grid,
                   const grid_ptr_type& to_grid) const;

  remapper_ptr_type
  create_remapper (const std::string& from_grid,
                   const std::string& to_grid) const {
    return create_remapper(get_grid(from_grid),get_grid(to_grid));
  }

  const grid_repo_type& get_repo () const { return m_grids; }

protected:
  using nonconstgrid_ptr_type     = std::shared_ptr<grid_type>;
  using nonconstgrid_repo_type    = std::map<std::string, nonconstgrid_ptr_type>;

  void add_grid (nonconstgrid_ptr_type grid);
  nonconstgrid_ptr_type get_grid_nonconst (const std::string& name) const;

  void alias_grid (const std::string& grid_name, const std::string& grid_alias);

  virtual remapper_ptr_type
  do_create_remapper (const grid_ptr_type from_grid,
                      const grid_ptr_type to_grid) const = 0;

  std::string print_available_grids () const;

private:

  grid_repo_type            m_grids;
  nonconstgrid_repo_type    m_nonconst_grids;
};

// A short name for the factory for grid managers
using GridsManagerFactory 
    = ekat::Factory<GridsManager,
                    ekat::CaseInsensitiveString,
                    std::shared_ptr<GridsManager>,
                    const ekat::Comm&,const ekat::ParameterList&>;

} // namespace scream

#endif // SCREAM_GRIDS_MANAGER_HPP
