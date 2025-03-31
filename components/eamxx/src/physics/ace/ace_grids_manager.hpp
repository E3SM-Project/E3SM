#ifndef EAMXX_ACE_GRIDS_MANAGER_HPP
#define EAMXX_ACE_GRIDS_MANAGER_HPP

#include "share/grid/grids_manager.hpp"

namespace scream {

class AceGridsManager : public GridsManager {
 public:
  AceGridsManager(const ekat::Comm &comm, const ekat::ParameterList &p);

  virtual ~AceGridsManager() = default;

  std::string name() const { return "ACE grids_manager"; }

  void build_grids() override;

 protected:
  remapper_ptr_type do_create_remapper(const grid_ptr_type from_grid,
                                       const grid_ptr_type to_grid) const {
    EKAT_ERROR_MSG(
        "Error! AceGridsManager is not capable of creating remappers.\n"
        " - from_grid: " +
        from_grid->name() +
        "\n"
        " - to_grid:   " +
        to_grid->name() + "\n");
    return nullptr;
  }

  ekat::Comm m_comm;
  ekat::ParameterList m_params;
};

inline std::shared_ptr<GridsManager> create_ace_grids_manager(
    const ekat::Comm &comm, const ekat::ParameterList &p) {
  return std::make_shared<AceGridsManager>(comm, p);
}

}  // namespace scream

#endif  // EAMXX_ACE_GRIDS_MANAGER_HPP
