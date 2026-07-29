#ifndef EAMXX_REMAPPER_MANAGER_HPP
#define EAMXX_REMAPPER_MANAGER_HPP

#include "abstract_remapper.hpp"

#include <ekat_factory.hpp>
#include <memory>

namespace scream {

class RemappersFactory : public ekat::Factory<AbstractRemapper,
                                              ekat::CaseInsensitiveString,
                                              std::shared_ptr<AbstractRemapper>,
                                              const ekat::ParameterList&>
{
public:
  using remapper_type     = AbstractRemapper;
  using remapper_ptr_type = std::shared_ptr<remapper_type>;

  RemappersManager (const std::shared_ptr<const GridsManager>& gm);
  virtual ~RemappersManager () = default;

  virtual remapper_ptr_type create_remapper (const grid_ptr_type& from_grid,
                                             const grid_ptr_type& to_grid) const;

  remapper_ptr_type create_remapper (const std::string& from_grid,
                                     const std::string& to_grid) const;
protected:

  std::shared_ptr<const GridsManager> m_grids_manager;
};

// A short name for the factory for grid managers
using RemappersManagerFactory
    = ekat::Factory<RemappersManager,
                    ekat::CaseInsensitiveString,
                    std::shared_ptr<RemappersManager>,
                    std::shared_ptr<const GridsManager>>;

} // namespace scream

#endif // EAMXX_REMAPPER_MANAGER_HPP
