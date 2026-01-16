#ifndef SCREAM_GW_DRAG_HPP
#define SCREAM_GW_DRAG_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/atm_process/ATMBufferManager.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "physics/gw/gw_functions.hpp"

#include <ekat_parameter_list.hpp>
#include <string>

namespace scream
{

/* Gravity Wave Drag Parameterization Suite

This suite of parameterizations can represent the drag from these sources:
  orographic
  frontogenesis
  deep convection
*/

class GWDrag : public AtmosphereProcess
{

  using KT  = ekat::KokkosTypes<DefaultDevice>;
  using GWF = gw::Functions<Real, DefaultDevice>;
  using PF  = scream::PhysicsFunctions<DefaultDevice>;
  using PC  = scream::physics::Constants<Real>;

  using Scalar   = typename GWF::Scalar;
  using Spack    = typename GWF::Spack;
  using SPackInt = typename GWF::SPackInt;

  public:
    // Constructors
    GWDrag (const ekat::Comm& comm, const ekat::ParameterList& params);

    // The type of subcomponent
    AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

    // The name of the subcomponent
    std::string name () const { return "gw"; }

    // Set the grid
    void create_requests ();

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
  protected:
#endif

    void initialize_impl (const RunType run_type);
    void run_impl        (const double dt);

  protected:

    void finalize_impl   ();

    std::shared_ptr<const AbstractGrid> m_grid;
    int m_ncol;
    int m_nlev;

    // Structures for arguments to ZM
    GWF::gw_runtime_opt gw_opts;

}; // class GWDrag

} // namespace scream

#endif // SCREAM_GW_DRAG_HPP
