#ifndef EAMXX_MAM_OPTICS_HPP
#define EAMXX_MAM_OPTICS_HPP

#include <share/atm_process/atmosphere_process.hpp>
#include <share/util/scream_common_physics_functions.hpp>
#include <share/atm_process/ATMBufferManager.hpp>

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/ekat_workspace.hpp>
#include <ekat/logging/ekat_logger.hpp>
#include <mam4xx/mam4.hpp>

#include <string>

#ifndef KOKKOS_ENABLE_CUDA
#define protected public
#define private public
#endif

namespace scream
{

// The process responsible for handling MAM4 aerosol optical properties. The AD
// stores exactly ONE instance of this class in its list of subcomponents.
class MAMOptics final : public scream::AtmosphereProcess {
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // view type that stores optical properties
  using view_3d = typename KT::template view_3d<Real>;

  // a quantity stored in a single vertical column with a single index
  using ColumnView = mam4::ColumnView;

  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;

  // a logger for this process
  using Logger = ekat::logger::Logger<ekat::logger::LogNoFile,
                                      ekat::logger::LogRootRank>;

public:

  // Constructor
  MAMOptics(const ekat::Comm& comm, const ekat::ParameterList& params);

protected:

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // process metadata
  AtmosphereProcessType type() const override;
  std::string name() const override;

  // grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager) override;

  // process behavior
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override;

private:

  Logger logger;

  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  // number of shortwave and longwave radiation bands
  int nswbands_, nlwbands_;

  // MAM4 aerosol particle size description
  mam4::AeroConfig aero_config_;

  // aerosol processes
  //std::unique_ptr<mam4::OpticsProcess> optics_;

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;

  // views that store optical parameters, dimensions are (ncol, nbands, nlev),
  // where nbands is swbands for shortwave quantities and lwbands for longwave
  view_3d aero_g_sw_;   // [-]
  view_3d aero_ssa_sw_; // [-]
  view_3d aero_tau_sw_; // [m] -- assumed SI length units (verify)
  view_3d aero_tau_lw_; // [m] -- assumed SI length units (verify)
}; // MAMOptics

} // namespace scream

#endif // EAMXX_MAM_OPTICS_HPP
