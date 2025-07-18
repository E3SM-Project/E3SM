#ifndef EAMXX_ZM_PROCESS_INTERFACE_HPP
#define EAMXX_ZM_PROCESS_INTERFACE_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"
#include "zm_functions.hpp"

namespace scream
{

// Zhang-McFarlane Deep Convection scheme

class zm_deep_convection : public AtmosphereProcess
{
  using KT  = ekat::KokkosTypes<DefaultDevice>;
  using ZMF = zm::Functions<Real, DefaultDevice>;
  using PF  = scream::PhysicsFunctions<DefaultDevice>;
  using PC  = scream::physics::Constants<Real>;
  
  using Spack                = typename ZMF::Spack;
  using SPackInt             = typename ZMF::SPackInt;

  using view_1d_int          = typename KT::template view_1d<Int>;
  using view_1d              = typename ZMF::view_1d<Real>;
  using view_1d_const        = typename ZMF::view_1d<const Real>;
  using view_2d              = typename ZMF::view_2d<ZMF::Spack>;
  using view_2dl             = typename ZMF::view_2dl<ZMF::Spack>;
  using view_2d_const        = typename ZMF::view_2d<const Spack>;
  using view_3d              = typename ZMF::view_3d<Spack>;
  using view_3d_const        = typename ZMF::view_3d<const Spack>;
  using view_3d_strided      = typename ZMF::view_3d_strided<Spack>;

  using uview_1d  = Unmanaged<view_1d>;
  using uview_2d  = Unmanaged<view_2d>;
  using uview_2dl = Unmanaged<view_2dl>;

  public:

    // Constructors
    zm_deep_convection(const ekat::Comm& comm, const ekat::ParameterList& params);

    // The type of subcomponent
    AtmosphereProcessType type() const override { return AtmosphereProcessType::Physics; }

    // The name of the subcomponent
    std::string name() const override { return "ZM"; }

    // Set the grid
    void set_grids(const std::shared_ptr<const GridsManager> grids_manager) override;

  protected:
    void initialize_impl (const RunType run_type) override;
    void run_impl        (const double dt) override;
    void finalize_impl   () override;

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
    ZMF::zm_buffer_data zm_buff;
    
}; // class zm_deep_convection

} // namespace scream

#endif // EAMXX_ZM_PROCESS_INTERFACE_HPP
