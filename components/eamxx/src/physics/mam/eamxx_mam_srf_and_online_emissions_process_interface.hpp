#ifndef EAMXX_MAM_SRF_ONLINE_EMISS_HPP
#define EAMXX_MAM_SRF_ONLINE_EMISS_HPP

#include "share/grid/remap/abstract_remapper.hpp"
#include "share/io/scorpio_input.hpp"

// For MAM4 aerosol configuration
#include <physics/mam/mam_coupling.hpp>
#include <physics/mam/srf_emission.hpp>

// For reading marine organics file
#include <physics/mam/readfiles/marine_organics.hpp>

// For declaring surface and online emission class derived from atm process
// class
#include <share/atm_process/atmosphere_process.hpp>
#include <string>

namespace scream {

// The process responsible for handling MAM4 surface and online emissions. The
// AD stores exactly ONE instance of this class in its list of subcomponents.
class MAMSrfOnlineEmiss final : public scream::AtmosphereProcess {
  using KT            = ekat::KokkosTypes<DefaultDevice>;
  using view_1d       = typename KT::template view_1d<Real>;
  using view_2d       = typename KT::template view_2d<Real>;
  using const_view_1d = typename KT::template view_1d<const Real>;
  using const_view_2d = typename KT::template view_2d<const Real>;

  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  // Wet and dry states of atmosphere
  mam_coupling::WetAtmosphere wet_atm_;
  mam_coupling::DryAtmosphere dry_atm_;

  // buffer for sotring temporary variables
  mam_coupling::Buffer buffer_;

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;

  // Sea surface temoerature [K]
  const_view_1d sst_;

  // Dust fluxes (four values for each col) [kg/m2/s]
  const_view_2d dust_fluxes_;

  // Constituent fluxes of species in [kg/m2/s]
  view_2d constituent_fluxes_;

  // Work array to store fluxes after unit conversions to kg/m2/s
  view_1d fluxes_in_mks_units_;

  // Unified atomic mass unit used for unit conversion (BAD constant)
  static constexpr Real amufac = 1.65979e-23;  // 1.e4* kg / amu

  // For reading soil erodibility file
  std::shared_ptr<AbstractRemapper> serod_horizInterp_;
  std::shared_ptr<AtmosphereInput> serod_dataReader_;
  const_view_1d soil_erodibility_;

 public:
  // For reading surface emissions and marine organics file
  using srfEmissFunc = mam_coupling::srfEmissFunctions<Real, DefaultDevice>;
  using marineOrganicsFunc =
      marine_organics::marineOrganicsFunctions<Real, DefaultDevice>;

  // Constructor
  MAMSrfOnlineEmiss(const ekat::Comm &comm, const ekat::ParameterList &params);

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // The type of subcomponent
  AtmosphereProcessType type() const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name() const { return "mam_srf_online_emissions"; }

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
  // Atmosphere processes often have a pre-processing step that constructs
  // required variables from the set of fields stored in the field manager.
  // This functor implements this step, which is called during run_impl.
  struct Preprocess {
    Preprocess() = default;
    // on host: initializes preprocess functor with necessary state data
    void initialize(const int &ncol, const int &nlev,
                    const mam_coupling::WetAtmosphere &wet_atm,
                    const mam_coupling::DryAtmosphere &dry_atm) {
      ncol_pre_    = ncol;
      nlev_pre_    = nlev;
      wet_atm_pre_ = wet_atm;
      dry_atm_pre_ = dry_atm;
    }
    KOKKOS_INLINE_FUNCTION
    void operator()(
        const Kokkos::TeamPolicy<KT::ExeSpace>::member_type &team) const {
      const int icol = team.league_rank();  // column index

      compute_dry_mixing_ratios(team, wet_atm_pre_, dry_atm_pre_, icol);
      team.team_barrier();
      // vertical heights has to be computed after computing dry mixing ratios
      // for atmosphere
      compute_vertical_layer_heights(team, dry_atm_pre_, icol);
      compute_updraft_velocities(team, wet_atm_pre_, dry_atm_pre_, icol);
    }  // Preprocess operator()

    // local variables for preprocess struct
    // number of horizontal columns and vertical levels
    int ncol_pre_, nlev_pre_;

    // local atmospheric and aerosol state data
    mam_coupling::WetAtmosphere wet_atm_pre_;
    mam_coupling::DryAtmosphere dry_atm_pre_;
  };  // MAMSrfOnlineEmiss::Preprocess
 private:
  // preprocessing scratch pad
  Preprocess preprocess_;

  // Species index (zero-based) in tracer array with "pcnst" dimension
  // FIXME: Remove the hardwired indices and use a function
  // to find them from an array.
  const std::map<std::string, int> spcIndex_in_pcnst_ = {
      {"so2", 12},    {"dms", 13},    {"so4_a1", 15}, {"dst_a1", 19},
      {"ncl_a1", 20}, {"mom_a1", 21}, {"num_a1", 22}, {"so4_a2", 23},
      {"ncl_a2", 25}, {"mom_a2", 26}, {"num_a2", 27}, {"dst_a3", 28},
      {"ncl_a3", 29}, {"num_a3", 35}, {"pom_a4", 36}, {"bc_a4", 37},
      {"mom_a4", 38}, {"num_a4", 39}};

  // A struct carrying all the fields needed to read
  // surface emissions of a species
  struct srf_emiss_ {
    // species name
    std::string species_name;

    // Data file name
    std::string data_file;

    // Sector names in file
    std::vector<std::string> sectors;

    // Data structure for reading interpolation
    std::shared_ptr<AbstractRemapper> horizInterp_;
    std::shared_ptr<AtmosphereInput> dataReader_;
    srfEmissFunc::srfEmissTimeState timeState_;
    srfEmissFunc::srfEmissInput data_start_, data_end_;
    srfEmissFunc::srfEmissOutput data_out_;
  };

  // A vector for carrying emissions for all the species
  std::vector<srf_emiss_> srf_emiss_species_;

  // For reading marine organics file
  std::shared_ptr<AbstractRemapper> morg_horizInterp_;
  std::shared_ptr<AtmosphereInput> morg_dataReader_;
  marineOrganicsFunc::marineOrganicsTimeState morg_timeState_;
  marineOrganicsFunc::marineOrganicsInput morg_data_start_, morg_data_end_;
  marineOrganicsFunc::marineOrganicsOutput morg_data_out_;

  // offset for converting pcnst index to gas_pcnst index
  static constexpr int offset_ =
      mam4::aero_model::pcnst - mam4::gas_chemistry::gas_pcnst;

};  // MAMSrfOnlineEmiss

}  // namespace scream

#endif  // EAMXX_MAM_SRF_ONLINE_EMISS_HPP
