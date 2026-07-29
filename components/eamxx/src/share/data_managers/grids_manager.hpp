#ifndef SCREAM_GRIDS_MANAGER_HPP
#define SCREAM_GRIDS_MANAGER_HPP

#include "share/grid/abstract_grid.hpp"

#include <ekat_factory.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_comm.hpp>

#include <map>
#include <set>
#include <memory>

namespace scream
{

class GridsManager
{
public:
  using grid_type              = AbstractGrid;
  using grid_ptr_type          = std::shared_ptr<const grid_type>;
  using grid_repo_type         = std::map<std::string, grid_ptr_type>;
  using nonconstgrid_ptr_type  = std::shared_ptr<grid_type>;
  using nonconstgrid_repo_type = std::map<std::string, nonconstgrid_ptr_type>;

  GridsManager () = default;
  virtual ~GridsManager () = default;

  virtual std::string name () const = 0;

  grid_ptr_type get_grid (const std::string& name) const;
  nonconstgrid_ptr_type get_grid_nonconst (const std::string& name) const;

  // Check if the given grid has been built
  bool has_grid (const std::string& grid_name) const;

  virtual void build_grids () = 0;

  const grid_repo_type& get_repo () const { return m_grids; }

  std::set<std::string> get_grid_names () const;
  std::string print_available_grids () const;

  // Return number of grids in the GridsManager
  int size() const { return m_grids.size(); }

protected:

  void add_nonconst_grid (nonconstgrid_ptr_type grid);
  void add_grid (grid_ptr_type grid);
  void alias_grid (const std::string& grid_name, const std::string& grid_alias);

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
