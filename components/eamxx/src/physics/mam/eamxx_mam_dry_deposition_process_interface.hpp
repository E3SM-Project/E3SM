#ifndef EAMXX_MAM_DRYDEP_HPP
#define EAMXX_MAM_DRYDEP_HPP

// For declaring dry deposition class derived from atm process class
#include <physics/mam/eamxx_mam_generic_process_interface.hpp>

// For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>

// For AtmosphereInput
#include "share/io/scorpio_input.hpp"

// For component name
#include <string>

namespace scream {

// The process responsible for handling MAM4 dry deposition. The AD
// stores exactly ONE instance of this class in its list of subcomponents.
class MAMDryDep final : public MAMGenericInterface {
 public:
  static constexpr int num_aero_modes = mam_coupling::num_aero_modes();
  static constexpr int aerosol_categories_ =
      mam4::DryDeposition::aerosol_categories;
  static constexpr int n_land_type = mam4::DryDeposition::n_land_type;

  using view_1d       = Field::view_dev_t<Real *>;
  using view_2d       = Field::view_dev_t<Real **>;
  using view_3d       = Field::view_dev_t<Real ***>;
  using view_4d       = Field::view_dev_t<Real ****>;
  using const_view_1d = Field::view_dev_t<const Real *>;
  using const_view_2d = Field::view_dev_t<const Real **>;
  using const_view_3d = Field::view_dev_t<const Real ***>;

 private:
  /* Note on mam4::DryDeposition::aerosol_categories = 4
     used in deposition velocity dimension defined below. These
     correspond to the two attachment states and two moments:
     0 - interstitial aerosol, 0th moment (i.e., number)
     1 - interstitial aerosol, 3rd moment (i.e., volume/mass)
     2 - cloud-borne aerosol,  0th moment (i.e., number)
     3 - cloud-borne aerosol,  3rd moment (i.e., volume/mass)
     see comments in the DryDeposition class in mam4xx.
  */
  // Output deposition velocity of turbulent dry deposition [m/s]
  // Dimensions
  //    [numer of modes, aerosol_categories_, num columns]
  view_3d vlc_trb_;

  // Output deposition velocity of gravitational settling [m/s]
  // Dimensions
  //   [num_modes, aerosol_categories_, num columns, num levels]
  view_4d vlc_grv_;

  // Output deposition velocity, [m/s]
  // fraction landuse weighted sum of vlc_grv and vlc_trb
  // Dimensions
  //   [num_modes, aerosol_categories_, num columns, num levels]
  view_4d vlc_dry_;

  // Output of the the mixing ratio tendencies [kg/kg/s or 1/kg/s]
  // Dimensions
  //   [num columns, num levels, mam4::aero_model::pcnst]
  // Packed the same way qtracers_ is layed out.
  view_3d ptend_q_;

  // Work array to hold the mixing ratios [kg/kg or 1/kg]
  // Dimensions
  //   [num columns, num levels, mam4::aero_model::pcnst]
  // Packs AerosolState::int_aero_nmr
  // and   AerosolState::int_aero_nmr
  // into one array, hence is mixed kg/kg and 1/kg.
  view_3d qtracers_;

  // Work array to hold the air density [kg/m3]
  // Dimensions
  //   [num columns, num levels]
  // Calculated from air pressure at layer midpoint,
  // Constants::r_gas_dry_air and air temperture.
  view_2d rho_;

  // Work array to hold tendency for 1 species [kg/kg/s] or [1/kg/s]
  // Dimensions
  //   [mam4::aero_model::pcnst, num column, num level]
  view_3d dqdt_tmp_;

  // Work array to hold cloud borne aerosols mixing ratios [kg/kg or 1/kg]
  // Dimensions
  //   [mam4::aero_model::pcnst, num column, num level]
  // Filled with Prognostics::n_mode_c and Prognostics::q_aero_c
  view_3d qqcw_;

  // For reading fractional land use file
  std::shared_ptr<AbstractRemapper> horizInterp_;
  std::shared_ptr<AtmosphereInput> dataReader_;
  const_view_2d frac_landuse_;
  view_2d frac_landuse_fm_;
  // aerosol state variables
  mam_coupling::AerosolState wet_aero_, dry_aero_;
  // wet mixing ratios (water species)
  mam_coupling::WetAtmosphere wet_atm_;
  // dry mixing ratios (water species)
  mam_coupling::DryAtmosphere dry_atm_;
  // workspace manager for internal local variables
  mam_coupling::Buffer buffer_;

  int get_len_temporary_views();
  void init_temporary_views();
  int len_temporary_views_{0};

 public:
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // Constructor
  MAMDryDep(const ekat::Comm &comm, const ekat::ParameterList &params);

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // The name of the subcomponent
  std::string name() const override { return "mam_dry_deposition"; }

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
  void finalize_impl() override{/*Do nothing*/};

 private:
  // pre- and postprocessing scratch pads
};  // MAMDryDep

}  // namespace scream

#endif  // EAMXX_MAM_DRYDEP_HPP
