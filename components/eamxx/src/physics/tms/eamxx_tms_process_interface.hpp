#ifndef SCREAM_TMS_PROCESS_HPP
#define SCREAM_TMS_PROCESS_HPP

#include "physics/tms/tms_functions.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"
#include "ekat/ekat_parameter_list.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible for computing/applying the surface drag coefficient
 * and stress associated with subgrid mountains
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
*/

class TurbulentMountainStress : public AtmosphereProcess
{
  using PF           = scream::PhysicsFunctions<DefaultDevice>;
  using TMSFunctions = tms::Functions<Real, DefaultDevice>;
  using Spack        = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using view_2d      = TMSFunctions::view_2d<Spack>;
  using uview_2d     = ekat::Unmanaged<view_2d>;

public:

  // Constructors
  TurbulentMountainStress (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "tms"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Structure for storing local variables initialized using the ATMBufferManager
  struct Buffer {
    static constexpr int num_2d_midpoint_views = 3;
    static constexpr int num_2d_interface_views = 1;

    uview_2d exner, dz, z_mid, z_int;
  };

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif

  void run_impl        (const double dt);

protected:

  void initialize_impl (const RunType run_type);
  void finalize_impl   ();

  // Computes total number of bytes needed for local variables
  size_t requested_buffer_size_in_bytes() const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager &buffer_manager);

  // Struct which contains local variables
  Buffer m_buffer;

  // Keep track of field dimensions and the iteration count
  int m_ncols;
  int m_nlevs;

  std::shared_ptr<const AbstractGrid> m_grid;
}; // class TurbulentMountainStress

} // namespace scream

#endif // SCREAM_TMS_PROCESS_HPP
