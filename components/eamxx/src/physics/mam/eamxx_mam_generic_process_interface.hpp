#ifndef EAMXX_MAM_GENERIC_PROCESS_HPP
#define EAMXX_MAM_GENERIC_PROCESS_HPP
// For declaring contituent fluxes class derived from atm process
// class
#include <share/atm_process/atmosphere_process.hpp>
#include <share/property_checks/field_within_interval_check.hpp>
// For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>
#include <physics/mam/physical_limits.hpp>
#include <string>

namespace scream {
class MAMGenericInterface : public scream::AtmosphereProcess {
 public:
   using KT = ekat::KokkosTypes<DefaultDevice>;

   // Constructor
  MAMGenericInterface(const ekat::Comm &comm,
                       const ekat::ParameterList &params);

  // void add_invariant_check_for_aerosol();
  void add_aerosol_tracers();
  void add_interval_checks();
  void print_fields_names();
  void populate_wet_and_dry_aero();
  void populate_wet_and_dry_atm();
  void add_tracer_for_wet_and_dry_atm();
  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
  // aerosol state variables
  mam_coupling::AerosolState wet_aero_, dry_aero_;
    // wet mixing ratios (water species)
  mam_coupling::WetAtmosphere wet_atm_;

  // dry mixing ratios (water species)
  mam_coupling::DryAtmosphere dry_atm_;
  // workspace manager for internal local variables
  mam_coupling::Buffer buffer_;
  std::vector<std::string> wet_atm_names_ = {"qv", "qc", "nc", "qi", "ni"};
  std::vector<std::string> dry_atm_names_ = {
        "T_mid",
        "p_mid",
        "p_int",
        "pseudo_density",
        "omega",
        "pbl_height",
        "cldfrac_tot"
    };
  bool check_fields_intervals_{false};
    // Atmosphere processes often have a pre-processing step that constructs
  // required variables from the set of fields stored in the field manager.
  // This functor implements this step, which is called during run_impl.
  struct Preprocess {
    Preprocess() = default;
    // on host: initializes preprocess functor with necessary state data
    void initialize(const int ncol, const int nlev,
                    const mam_coupling::WetAtmosphere &wet_atm,
                    const mam_coupling::AerosolState &wet_aero,
                    const mam_coupling::DryAtmosphere &dry_atm,
                    const mam_coupling::AerosolState &dry_aero) {
      ncol_pre_     = ncol;
      nlev_pre_     = nlev;
      wet_atm_pre_  = wet_atm;
      wet_aero_pre_ = wet_aero;
      dry_atm_pre_  = dry_atm;
      dry_aero_pre_ = dry_aero;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const Kokkos::TeamPolicy<KT::ExeSpace>::member_type &team) const {
      const int i = team.league_rank();  // column index

      compute_dry_mixing_ratios(team, wet_atm_pre_, dry_atm_pre_, i);
      compute_dry_mixing_ratios(team, wet_atm_pre_, wet_aero_pre_,
                                dry_aero_pre_, i);
      team.team_barrier();
      // vertical heights has to be computed after computing dry mixing ratios
      // for atmosphere
      compute_vertical_layer_heights(team, dry_atm_pre_, i);
      compute_updraft_velocities(team, wet_atm_pre_, dry_atm_pre_, i);
      // allows kernels below to use layer heights operator()
      team.team_barrier();
      // set_min_background_mmr(team, dry_aero_pre_,
      //                        i);  // dry_atm_pre_ is the output
    }                             // operator()

    // local variables for preprocess struct
    // number of horizontal columns and vertical levels
    int ncol_pre_, nlev_pre_;

    // local atmospheric and aerosol state data
    mam_coupling::WetAtmosphere wet_atm_pre_;
    mam_coupling::DryAtmosphere dry_atm_pre_;
    mam_coupling::AerosolState wet_aero_pre_, dry_aero_pre_;
  };  // MAMAci::Preprocess

  // Atmosphere processes often have a post-processing step prepares output
  // from this process for the Field Manager. This functor implements this
  // step, which is called during run_impl.
  // Postprocessing functor
  struct Postprocess {
    Postprocess() = default;

    // on host: initializes postprocess functor with necessary state data
    void initialize(const int ncol, const int nlev,
                    const mam_coupling::WetAtmosphere &wet_atm,
                    const mam_coupling::AerosolState &wet_aero,
                    const mam_coupling::DryAtmosphere &dry_atm,
                    const mam_coupling::AerosolState &dry_aero) {
      ncol_post_     = ncol;
      nlev_post_     = nlev;
      wet_atm_post_  = wet_atm;
      wet_aero_post_ = wet_aero;
      dry_atm_post_  = dry_atm;
      dry_aero_post_ = dry_aero;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(
        const Kokkos::TeamPolicy<KT::ExeSpace>::member_type &team) const {
      const int i = team.league_rank();  // column index
      compute_wet_mixing_ratios(team, dry_atm_post_, dry_aero_post_,
                                wet_aero_post_, i);
    }  // operator()

    // number of horizontal columns and vertical levels
    int ncol_post_, nlev_post_;

    // local atmospheric and aerosol state data
    mam_coupling::WetAtmosphere wet_atm_post_;
    mam_coupling::DryAtmosphere dry_atm_post_;
    mam_coupling::AerosolState wet_aero_post_, dry_aero_post_;
  };  // Postprocess

 private:
   // The type of subcomponent
  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------
  AtmosphereProcessType type() const { return AtmosphereProcessType::Physics; }
  std::map<std::string, std::pair<Real, Real>>  limits_aerosol_gas_tracers_;
  void get_aerosol_gas_map();
  const std::pair<Real, Real> get_range(const std::string &field_name);

};  // MAMGenericInterface
}  // namespace scream

#endif  // ifdef EAMXX_MAM_CONSTITUTE_FLUXES_FUNCTIONS_HPP
