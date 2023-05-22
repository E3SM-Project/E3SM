#ifndef SCREAM_HOMME_GRIDS_MANAGER_HPP
#define SCREAM_HOMME_GRIDS_MANAGER_HPP

#include "share/grid/grids_manager.hpp"

namespace scream
{

// This class is meant to be used for small unit tests, where we want to
// test the Atmosphere Driver (AD) capabilities, without bothering too much
// about grids-related features. This manager lets you set pre-built grids
// in it rather than building them inside the manager.
class HommeGridsManager : public GridsManager
{
public:
  using ci_string = ekat::CaseInsensitiveString;

  HommeGridsManager (const ekat::Comm& comm, const ekat::ParameterList& p);

  ~HommeGridsManager ();

  std::string name () const { return "Homme Grids Manager"; }

  void build_grids ();

#ifndef KOKKOS_ENABLE_CUDA
protected:
#endif

  void build_dynamics_grid ();
  void build_physics_grid  (const ci_string& type,
                            const ci_string& rebalance);

protected:

  remapper_ptr_type
  do_create_remapper (const grid_ptr_type from_grid,
                      const grid_ptr_type to_grid) const;

  void build_pg_codes ();

  // Read vertical coordinates and set them in hommexx's structures
  void initialize_vertical_coordinates (const nonconstgrid_ptr_type& dyn_grid);

  ekat::Comm            m_comm;

  ekat::ParameterList   m_params;

  // For each admissible physics grid type/rebalance, store an integer code
  // We pass these codes to f90, rather than a bunch of strings
  template<typename T>
  using strmap_t = std::map<ekat::CaseInsensitiveString,T>;
  strmap_t<strmap_t<int>>   m_pg_codes;
};

inline std::shared_ptr<GridsManager>
create_homme_grids_manager (const ekat::Comm& comm, const ekat::ParameterList& p) {
  return std::make_shared<HommeGridsManager>(comm,p);
}

} // namespace scream

#endif // SCREAM_HOMME_GRIDS_MANAGER_HPP
