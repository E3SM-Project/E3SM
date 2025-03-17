#ifndef EAMXX_MAM_OPTICS_HPP
#define EAMXX_MAM_OPTICS_HPP

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/ekat_workspace.hpp>
#include <mam4xx/mam4.hpp>
#include <physics/mam/eamxx_mam_generic_process_interface.hpp>
#include <physics/mam/mam_aerosol_optics_read_tables.hpp>
#include <physics/mam/mam_coupling.hpp>
#include <share/atm_process/ATMBufferManager.hpp>
#include <share/util/eamxx_common_physics_functions.hpp>
#include <string>

#ifndef KOKKOS_ENABLE_CUDA
#define protected_except_cuda public
#define private_except_cuda public
#else
#define protected_except_cuda protected
#define private_except_cuda private
#endif

namespace scream {

// The process responsible for handling MAM4 aerosol optical properties. The AD
// stores exactly ONE instance of this class in its list of subcomponents.
class MAMOptics final : public MAMGenericInterface {
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // a quantity stored in a single vertical column with a single index
  using ColumnView = mam4::ColumnView;

  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;

 public:
  // Constructor
  MAMOptics(const ekat::Comm &comm, const ekat::ParameterList &params);

  protected_except_cuda :

      // --------------------------------------------------------------------------
      // AtmosphereProcess overrides (see
      // share/atm_process/atmosphere_process.hpp)
      // --------------------------------------------------------------------------
      std::string
      name() const override;

  // grid
  void set_grids(
      const std::shared_ptr<const GridsManager> grids_manager) override;

  // management of common atm process memory
  size_t requested_buffer_size_in_bytes() const override;
  void init_buffers(const ATMBufferManager &buffer_manager) override;

  // process behavior
  void initialize_impl(const RunType run_type) override;
  void run_impl(const double dt) override;
  void finalize_impl() override;

  // number of shortwave and longwave radiation bands
  int nswbands_, nlwbands_;

  // FIXME: move these values to mam_coupling
  mam_coupling::const_view_2d p_int_, p_del_;

  // MAM4 aerosol particle size description
  mam4::AeroConfig aero_config_;

  mam_coupling::view_3d ssa_cmip6_sw_, af_cmip6_sw_, ext_cmip6_sw_;
  // long wave extinction in the units of [1/km]
  mam_coupling::view_3d ext_cmip6_lw_;
  mam4::modal_aer_opt::AerosolOpticsDeviceData aerosol_optics_device_data_;
  // physics grid for column information
  mam_coupling::view_2d work_;
  mam_coupling::view_3d tau_ssa_g_sw_, tau_ssa_sw_, tau_sw_, tau_f_sw_;
  // Mapping from old RRTMG sw bands to new band ordering in RRTMGP
  //  given old index swband (RRTMG) return new index swband RRTMGP
  mam_coupling::view_int_1d get_idx_rrtmgp_from_rrtmg_swbands_;
  // aerosol state variables
  mam_coupling::AerosolState wet_aero_, dry_aero_;
  // wet mixing ratios (water species)
  mam_coupling::WetAtmosphere wet_atm_;
  // dry mixing ratios (water species)
  mam_coupling::DryAtmosphere dry_atm_;
  // workspace manager for internal local variables
  mam_coupling::Buffer buffer_;

  int work_len_;

  int num_2d_scratch_= 8;
};  // MAMOptics

}  // namespace scream

#endif  // EAMXX_MAM_OPTICS_HPP
