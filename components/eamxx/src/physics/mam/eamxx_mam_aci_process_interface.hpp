#ifndef EAMXX_MAM_ACI_HPP
#define EAMXX_MAM_ACI_HPP

#include <share/atm_process/atmosphere_process.hpp>
#include <share/util/scream_common_physics_functions.hpp>
#include <share/atm_process/ATMBufferManager.hpp>

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/ekat_workspace.hpp>
#include <mam4xx/mam4.hpp>

#include <string>

#ifndef KOKKOS_ENABLE_CUDA
#define protected_except_cuda public
#define private_except_cuda public
#else
#define protected_except_cuda protected
#define private_except_cuda private
#endif

namespace scream
{

// The process responsible for handling MAM4 aerosol optical properties. The AD
// stores exactly ONE instance of this class in its list of subcomponents.
class MAMAci final : public scream::AtmosphereProcess {
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // a quantity stored in a single vertical column with a single index
  using ColumnView = mam4::ColumnView;

  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;

public:

  // Constructor
  MAMAci(const ekat::Comm& comm, const ekat::ParameterList& params);

protected_except_cuda:

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

private_except_cuda:

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
}; // MAMAci

} // namespace scream


#endif // EAMXX_MAM_ACI_HPP
