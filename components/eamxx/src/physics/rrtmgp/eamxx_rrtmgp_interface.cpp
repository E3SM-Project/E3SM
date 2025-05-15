#include "eamxx_rrtmgp_interface.hpp"
#include "physics/share/physics_constants.hpp"

namespace scream {

void init_kls ()
{
  // Initialize kokkos
  if(!Kokkos::is_initialized()) { Kokkos::initialize(); }
}

void finalize_kls()
{
  //Kokkos::finalize(); We do the kokkos finalization elsewhere
}

}  // namespace scream
