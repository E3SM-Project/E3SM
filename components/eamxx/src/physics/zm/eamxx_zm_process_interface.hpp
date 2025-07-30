#ifndef EAMXX_ZM_PROCESS_INTERFACE_HPP
#define EAMXX_ZM_PROCESS_INTERFACE_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "physics/zm/zm_functions.hpp"

#include <ekat_parameter_list.hpp>

namespace scream
{

// Zhang-McFarlane Deep Convection scheme

class ZMDeepConvection : public AtmosphereProcess
{
  using KT  = ekat::KokkosTypes<DefaultDevice>;
  using ZMF = zm::Functions<Real, DefaultDevice>;
  
  using Spack                = typename ZMF::Spack;
  using SPackInt             = typename ZMF::SPackInt;

  using view_1d_int          = typename KT::template view_1d<Int>;
  using view_1d              = typename ZMF::view_1d<Real>;
  using view_1d_const        = typename ZMF::view_1d<const Real>;
  using view_2d              = typename ZMF::view_2d<ZMF::Spack>;
  using view_2d_const        = typename ZMF::view_2d<const Spack>;
  using view_3d              = typename ZMF::view_3d<Spack>;
  using view_3d_const        = typename ZMF::view_3d<const Spack>;
  using view_3d_strided      = typename ZMF::view_3d_strided<Spack>;

  public:

    // Constructors
    ZMDeepConvection(const ekat::Comm& comm, const ekat::ParameterList& params);

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

    // define ZM process variables
    std::shared_ptr<const AbstractGrid> m_grid;
    int m_ncols;
    int m_nlevs;
    
};

} // namespace scream

#endif // EAMXX_ZM_PROCESS_INTERFACE_HPP
