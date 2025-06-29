#ifndef EAMXX_MAM_INTERPOLATION_MICROPHYSICS_HPP
#define EAMXX_MAM_INTERPOLATION_MICROPHYSICS_HPP

#include <physics/mam/eamxx_mam_generic_process_interface.hpp>
#include <physics/mam/mam_coupling.hpp>
#include <share/util/eamxx_common_physics_functions.hpp>

// For calling MAM4 processes
#include <mam4xx/mam4.hpp>
#include <string>

namespace scream {

class DataInterpolation;


// The process responsible for handling MAM4 aerosol microphysics. The AD
// stores exactly ONE instance of this class in its list of subcomponents.
class MAMInterpolationMicrophysics final : public MAMGenericInterface {
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // views for single- and multi-column data
  using view_1d       = typename KT::template view_1d<Real>;
  using view_2d       = typename KT::template view_2d<Real>;
  using view_3d       = typename KT::template view_3d<Real>;
  using const_view_1d = typename KT::template view_1d<const Real>;
  using const_view_2d = typename KT::template view_2d<const Real>;

  using view_1d_host = typename KT::view_1d<Real>::HostMirror;

  using view_int_2d = typename KT::template view_2d<int>;

  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;

 public:
  // Constructor
  MAMInterpolationMicrophysics(const ekat::Comm &comm, const ekat::ParameterList &params);

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // The name of the subcomponent
  std::string name() const { return "mam_interpolation_aero_microphysics"; }

  // grid
  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

  // management of common atm process memory
  size_t requested_buffer_size_in_bytes() const override;
  void init_buffers(const ATMBufferManager &buffer_manager) override;

  // Initialize variables
  void initialize_impl(const RunType run_type) override;

  // Run the process by one time step
  void run_impl(const double dt) override;

  // Finalize
  void finalize_impl(){/*Do nothing*/};

 private:

  // wet mixing ratios (water species)
  mam_coupling::WetAtmosphere wet_atm_;
  // dry mixing ratios (water species)
  mam_coupling::DryAtmosphere dry_atm_;

  // workspace manager for internal local variables
  mam_coupling::Buffer buffer_;

  std::vector<std::string> m_var_names_oxi;
  std::shared_ptr<DataInterpolation>    m_data_interpolation;

  std::vector<std::string> m_var_names_linoz;
  std::shared_ptr<DataInterpolation>    m_data_interpolation_linoz;

  std::map<std::string, std::vector<std::string>> m_elevated_emis_var_names;
  std::vector<std::shared_ptr<DataInterpolation>> m_data_interpolation_vertical;

};  // MAMInterpolationMicrophysics

}  // namespace scream

#endif  // EAMXX_MAM_MICROPHYSICS_HPP
