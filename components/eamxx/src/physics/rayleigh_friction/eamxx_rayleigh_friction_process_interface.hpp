#ifndef SCREAM_RAYLEIGH_FRICTION_PROCESS_HPP
#define SCREAM_RAYLEIGH_FRICTION_PROCESS_HPP

#include "physics/rayleigh_friction/rayleigh_friction_functions.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include <ekat_parameter_list.hpp>

#include <string>

namespace scream
{

/*
 * The class responsible for applying Rayleigh friction near the model top,
 * ported from EAM (rayleigh_friction.F90). The horizontal winds are damped
 * toward zero and the dissipated kinetic energy is added to the temperature.
 *
 * NOTE: HommeDynamics also applies Rayleigh friction (controlled by the
 * homme::rayleigh_friction_decay_time parameter). Set that to 0 when
 * using this process to avoid applying Rayleigh friction twice.
*/

class RayleighFriction : public AtmosphereProcess
{
  using RFFunctions = rayleigh_friction::Functions<Real, DefaultDevice>;
  using view_1d     = RFFunctions::view_1d<Real>;

public:

  // Constructors
  RayleighFriction (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "rayleigh_friction"; }

  // Set the grid
  void create_requests ();

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif

  void run_impl        (const double dt);

protected:

  void initialize_impl (const RunType run_type);
  void finalize_impl   ();

  // Rayleigh friction parameters
  int  m_rayk0;     // 1-based vertical level at which the profile is centered
  Real m_raykrange; // width of the profile (if 0, set to (rayk0-1)/2)
  Real m_raytau0;   // approximate decay time at model top [days] (if 0, no friction)

  // Decay rate profile [1/s]
  view_1d m_otau;

  // Keep track of field dimensions
  int m_ncols;
  int m_nlevs;

  std::shared_ptr<const AbstractGrid> m_grid;
}; // class RayleighFriction

} // namespace scream

#endif // SCREAM_RAYLEIGH_FRICTION_PROCESS_HPP
