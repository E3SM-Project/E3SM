#ifndef EAMXX_MAM_MICROPHYSICS_HPP
#define EAMXX_MAM_MICROPHYSICS_HPP

#include <physics/mam/eamxx_mam_generic_process_interface.hpp>
#include <physics/mam/mam_coupling.hpp>
#include <share/util/eamxx_common_physics_functions.hpp>

#include "readfiles/tracer_reader_utils.hpp"
// For calling MAM4 processes
#include <mam4xx/mam4.hpp>
#include <string>


namespace scream {

// The process responsible for handling MAM4 aerosol microphysics. The AD
// stores exactly ONE instance of this class in its list of subcomponents.
class MAMMicrophysics final : public MAMGenericInterface {
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
  MAMMicrophysics(const ekat::Comm &comm, const ekat::ParameterList &params);

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // The name of the subcomponent
  std::string name() const { return "mam_aero_microphysics"; }

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
  // Output extra mam4xx diagnostics.
  bool extra_mam4_aero_microphys_diags_ = false;

  // The orbital year, used for zenith angle calculations:
  // If > 0, use constant orbital year for duration of simulation
  // If < 0, use year from timestamp for orbital parameters
  int m_orbital_year;

  // Orbital parameters, used for zenith angle calculations.
  // If >= 0, bypass computation based on orbital year and use fixed parameters
  // If <  0, compute based on orbital year, specified above
  // These variables are required to be double.
  double m_orbital_eccen;  // Eccentricity
  double m_orbital_obliq;  // Obliquity
  double m_orbital_mvelp;  // Vernal Equinox Mean Longitude of Perihelion

  struct Config {
    // stratospheric chemistry parameters
    struct {
      int o3_lbl;   // number of layers with ozone decay from the surface
      Real o3_sfc;  // set from namelist input linoz_sfc
      Real o3_tau;  // set from namelist input linoz_tau
      Real psc_T;   // set from namelist input linoz_psc_T
    } linoz;

    // aqueous chemistry parameters
    mam4::mo_setsox::Config setsox;

    // aero microphysics configuration (see impl/mam4_amicphys.cpp)
    mam4::microphysics::AmicPhysConfig amicphys;
  };
  Config config_;
  // MAM4 aerosol particle size description
  mam4::AeroConfig aero_config_;

  // photolysis rate table (column-independent)
  mam4::mo_photo::PhotoTableData photo_table_;

  // column areas, latitudes, longitudes
  const_view_1d col_latitudes_;

  // surface albedo: shortwave, direct
  const_view_1d d_sfc_alb_dir_vis_;


  view_2d work_photo_table_;
  std::vector<Real> chlorine_values_;
  std::vector<int> chlorine_time_secs_;
  view_3d photo_rates_;
#ifdef USE_OLD_LINOZ_FILE_READ
  // invariants members
  mam_coupling::TracerTimeState trace_time_state_;
  std::shared_ptr<AtmosphereInput> TracerDataReader_;
  std::shared_ptr<AbstractRemapper> TracerHorizInterp_;
  mam_coupling::TracerData tracer_data_;
  std::string oxid_file_name_;
  view_2d cnst_offline_[4];
#endif
  view_3d invariants_;
  // linoz reader
#ifdef USE_OLD_LINOZ_FILE_READ
  mam_coupling::TracerTimeState linoz_time_state_;
  std::shared_ptr<AtmosphereInput> LinozDataReader_;
  std::shared_ptr<AbstractRemapper> LinozHorizInterp_;
  mam_coupling::TracerData linoz_data_;
  std::string linoz_file_name_;
#endif
  // Vertical emission uses 9 files, here I am using std::vector to stote
  // instance of each file.
#ifdef USE_OLD_VERTICAL_FILE_READ
  mam_coupling::TracerTimeState elevated_emiss_time_state_[mam4::gas_chemistry::extcnt];
  std::vector<std::shared_ptr<AtmosphereInput>> ElevatedEmissionsDataReader_;
  std::vector<std::shared_ptr<AbstractRemapper>> ElevatedEmissionsHorizInterp_;
  std::vector<mam_coupling::TracerData> elevated_emis_data_;
  std::map<std::string, std::string> elevated_emis_file_name_;
#endif
  std::vector<std::string> extfrc_lst_;
  std::map<std::string, std::vector<std::string>> elevated_emis_var_names_;

  view_3d extfrc_;
  mam_coupling::ForcingHelper forcings_[mam4::gas_chemistry::extcnt];

  view_1d_host acos_cosine_zenith_host_;
  view_1d acos_cosine_zenith_;

  view_int_2d index_season_lai_;
  // // dq/dt for convection [kg/kg/s]
  view_1d cmfdqr_;
  view_2d work_set_het_;
  // aerosol state variables
  mam_coupling::AerosolState wet_aero_, dry_aero_;
  // wet mixing ratios (water species)
  mam_coupling::WetAtmosphere wet_atm_;
  // dry mixing ratios (water species)
  mam_coupling::DryAtmosphere dry_atm_;
  // workspace manager for internal local variables
  mam_coupling::Buffer buffer_;

  // Work arrays for return values from
  // perform_atmospheric_chemistry_and_microphysics
  view_2d dflx_;
  view_2d dvel_;
  int num_2d_scratch_ = 9;
  int get_len_temporary_views();
  void init_temporary_views();
  int len_temporary_views_{0};

};  // MAMMicrophysics

}  // namespace scream

#endif  // EAMXX_MAM_MICROPHYSICS_HPP
