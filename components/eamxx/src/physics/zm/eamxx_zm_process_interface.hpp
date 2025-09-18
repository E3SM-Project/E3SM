#ifndef EAMXX_ZM_PROCESS_INTERFACE_HPP
#define EAMXX_ZM_PROCESS_INTERFACE_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/atm_process/ATMBufferManager.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "zm_functions.hpp"

#include <ekat_parameter_list.hpp>

namespace scream
{

// Zhang-McFarlane Deep Convection scheme

class ZMDeepConvection : public AtmosphereProcess
{
  using KT  = ekat::KokkosTypes<DefaultDevice>;
  using ZMF = zm::Functions<Real, DefaultDevice>;
  using PF  = scream::PhysicsFunctions<DefaultDevice>;
  using PC  = scream::physics::Constants<Real>;
  
  using Scalar               = typename ZMF::Scalar;
  using Spack                = typename ZMF::Spack;
  using SPackInt             = typename ZMF::SPackInt;

  public:

    // Constructors
    ZMDeepConvection(const ekat::Comm& comm, const ekat::ParameterList& params);

    // The type of subcomponent
    AtmosphereProcessType type() const override { return AtmosphereProcessType::Physics; }

    // The name of the subcomponent
    std::string name() const override { return "ZM"; }

    // Set the grid
    void set_grids(const std::shared_ptr<const GridsManager> grids_manager) override;

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
  protected:
#endif

    void initialize_impl (const RunType run_type);
    void run_impl        (const double dt);

  protected:

    void finalize_impl   ();

    // Computes total number of bytes needed for local variables
    size_t requested_buffer_size_in_bytes() const;

    // Set local variables using memory provided by the ATMBufferManager
    void init_buffers(const ATMBufferManager &buffer_manager);

    // define ZM process variables
    std::shared_ptr<const AbstractGrid> m_grid;
    int m_pcol;
    int m_ncol;
    int m_nlev;

    // Structures for arguments to ZM
    ZMF::zm_runtime_opt zm_opts;
    ZMF::zm_input_state zm_input;
    ZMF::zm_output_tend zm_output;
    ZMF::zm_output_diag zm_diag;
    
}; // class ZMDeepConvection

} // namespace scream

#endif // EAMXX_ZM_PROCESS_INTERFACE_HPP
