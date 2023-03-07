#include "share/grid/grids_manager.hpp"

namespace scream
{

auto GridsManager::
get_grid(const std::string& name) const
 -> grid_ptr_type
{
  auto g = get_grid_nonconst(name);
  return g;
}

bool GridsManager::
has_grid (const std::string& grid_name) const
{
  for (const auto& it : m_grids) {
    const auto& g = it.second;
    if (ekat::contains(g->aliases(),grid_name)) {
      return true;
    }
  }
  return false;
}

auto GridsManager::
create_remapper (const grid_ptr_type& from_grid,
                 const grid_ptr_type& to_grid) const
 -> GridsManager::remapper_ptr_type
{
  EKAT_REQUIRE_MSG( has_grid(from_grid->name()),
                    "Error! Source grid '" + from_grid->name() + "' is not supported.\n");
  EKAT_REQUIRE_MSG( has_grid(to_grid->name()),
                    "Error! Target grid '" + to_grid->name() + "' is not supported.\n");

  remapper_ptr_type remapper;

  if (from_grid->name()==to_grid->name()) {
    // We can handle the identity remapper from here
    remapper = std::make_shared<IdentityRemapper>(from_grid);
  } else {
    remapper = do_create_remapper(from_grid,to_grid);
  }

  EKAT_REQUIRE_MSG(
    remapper!=nullptr,
    "Error! A remapper from grid '" + from_grid->name() + "' to grid '" + to_grid->name() + "' is not available.\n"
    "       Perhaps you forgot to add its creation to the implementation of the grids manager?\n");

  return remapper;
}

void GridsManager::
add_grid (nonconstgrid_ptr_type grid)
{
  const auto& name = grid->name();
  EKAT_REQUIRE_MSG (not has_grid(name),
      "Error! Cannot add grid: grid already added.\n"
      "  - grids manager: " + this->name() + "\n"
      "  - grid name    : " + name + "\n");

  m_nonconst_grids[name] = grid;
  m_grids[name] = grid;
}

auto GridsManager::
get_grid_nonconst (const std::string& name) const
 -> nonconstgrid_ptr_type
{
  EKAT_REQUIRE_MSG (has_grid(name),
                      "Error! Grids manager '" + this->name() + "' does not provide grid '" + name + "'.\n"
                      "       Avaialble grids are: " + print_available_grids()  + "\n");

  nonconstgrid_ptr_type g;
  for (const auto& it : m_nonconst_grids) {
    if (it.second->name()==name or
        ekat::contains(it.second->aliases(),name)) {
      g = it.second;
      break;
    }
  }

  EKAT_REQUIRE_MSG (g!=nullptr,
      "Something went wrong while looking up a grid.\n"
      "  - grids manager: " + this->name() + "\n"
      "  - grid name    : " + name   + "\n");

  return g;
  
}

void GridsManager::
alias_grid (const std::string& grid_name,
            const std::string& grid_alias)
{
  EKAT_REQUIRE_MSG (has_grid(grid_name),
      "Error! Cannot alias grid: grid not stored.\n"
      "  - grids manager: " + this->name() + "\n"
      "  - grid name    : " + grid_name + "\n"
      "  - grid alias   : " + grid_alias + "\n");

  if (grid_name==grid_alias or ekat::contains(get_grid(grid_name)->aliases(),grid_alias)) {
    // Unlikely to happen, but should be ok to simply return
    return;
  }

  auto g = get_grid_nonconst(grid_name);
  EKAT_REQUIRE_MSG (not has_grid(grid_alias),
      "Cannot alias grid, since another grid with alias name already exists.\n"
      "  - grids manager: " + this->name() + "\n"
      "  - grid name    : " + grid_name + "\n"
      "  - grid alias   : " + grid_alias + "\n");

  g->add_alias(grid_alias);
}

std::string GridsManager::
print_available_grids () const {
  std::string str;
  if (m_grids.size()==0) {
    return str;
  }
  for (const auto& g : m_grids) {
    if (g.first!=g.second->name()) {
      str += g.first + " (alias of " + g.second->name() + "), ";
    } else {
      str += g.first + ", ";
    }
  }
  str.erase(str.size()-2,2); // Erase trailing ', '
  return str;
}

} // namespace scream
