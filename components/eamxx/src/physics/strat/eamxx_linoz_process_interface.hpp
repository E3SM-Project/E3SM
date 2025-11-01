#ifndef EAMXX_STRAT_LINOZ_HPP
#define EAMXX_STRAT_LINOZ_HPP

#include <mam4xx/mam4.hpp>
#include <string>
#include "share/atm_process/atmosphere_process.hpp"

namespace scream {
class DataInterpolation;
// The process responsible for handling MAM4 aerosol microphysics. The AD
// stores exactly ONE instance of this class in its list of subcomponents.
class STRATLinoz : public AtmosphereProcess {
  using KT = ekat::KokkosTypes<DefaultDevice>;


  // views for single- and multi-column data
  using view_1d       = typename KT::template view_1d<Real>;
  using view_2d       = typename KT::template view_2d<Real>;
  using const_view_1d = typename KT::template view_1d<const Real>;

  using view_1d_host = typename KT::view_1d<Real>::HostMirror;

  // a thread team dispatched to a single vertical column
  using ThreadTeam = mam4::ThreadTeam;

 public:
  // Constructor
  STRATLinoz(const ekat::Comm &comm, const ekat::ParameterList &params);

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // The name of the subcomponent
  std::string name() const { return "strat_linoz"; }

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

  int get_len_temporary_views();
  void init_temporary_views();

  // column areas, latitudes, longitudes
  const_view_1d col_latitudes_;

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

  bool use_prescribed_ozone_;


  // stratospheric chemistry parameters
  mam4::microphysics::LinozConf m_config;

  // names of linoz field
  std::vector<std::string> m_var_names_linoz;

  // data interpolation object for linoz fields
  std::shared_ptr<DataInterpolation>    data_interp_linoz_;
  void set_linoz_reader();

  std::shared_ptr<const AbstractGrid> grid_;
  // number of horizontal columns and vertical levels
  int ncol_, nlev_;

  std::vector<Real> chlorine_values_;
  std::vector<int> chlorine_time_secs_;
  view_2d m_o3_col_dens;
  view_1d_host acos_cosine_zenith_host_;
  view_1d acos_cosine_zenith_;
  //
  view_2d m_vmr;
  view_2d m_o3_col_deltas;

};  // MAMMicrophysics

}  // namespace scream

#endif  // EAMXX_MAM_MICROPHYSICS_HPP
