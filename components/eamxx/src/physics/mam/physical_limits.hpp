// mam4xx: Copyright (c) 2022,
// Battelle Memorial Institute and
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MAM_PHYSICAL_LIMITS_HPP
#define MAM_PHYSICAL_LIMITS_HPP
#include <map>
#include <string>
#include <utility>
#include <share/atm_process/atmosphere_process.hpp>
#include <share/property_checks/field_within_interval_check.hpp>
// For MAM4 aerosol configuration

namespace scream::mam_coupling {

inline const std::pair<Real, Real> &
physical_min_max(const std::string &field_name) {
  static const std::map<std::string, std::pair<Real, Real>> limits = {
      {"T_mid", {100, 500}},
      {"p_mid", {0, 1e10}}, // FIXME
      {"qv", {1e-13, 0.2}}, {"qc", {0, 0.1}},
      {"qi", {0, 0.1}},      {"nc", {0, 0.1e11}},  {"nr", {0, 0.1e10}},
      {"ni", {0, 0.1e10}},   {"nmr", {0, 1e13}},   {"mmr", {-1e-20, 1e-2}},
      {"omega", {-1e10, 1e10}}, // FIXME
      {"p_int", {0, 1e10}}, // FIXME
      {"pseudo_density", {0, 1e10}}, // FIXME
      {"pbl_height", {0, 1e10}}, // FIXME
      {"cldfrac_tot", {0, 1e10}}, // FIXME
      // constituent_fluxes
      {"constituent_fluxes", {0, 1e10}}, // FIXME
      {"phis", {-1e10, 1e10}}, // FIXME
      // aci required
      {"w_variance", {-1e10, 1e10}}, // FIXME
      {"cldfrac_liq", {-1e10, 1e10}}, // FIXME
      {"cldfrac_liq_prev", {-1e10, 1e10}}, // FIXME
      {"eddy_diff_heat", {-1e10, 1e10}}, // FIXME
      {"dgnum", {-1e10, 1e10}}, // FIXME
      // aci compute
      {"ni_activated", {-1e100, 1e100}}, // FIXME
      {"nc_nuceat_tend", {-1e100, 1e100}}, // FIXME
      {"nsource", {-1e10, 1e10}}, // FIXME
      {"ndropmix", {-1e10, 1e10}}, // FIXME
      {"nc_inp_to_aci", {-1e10, 1e10}}, // FIXME
      {"ccn_0p02", {-1e10, 1e10}}, // FIXME
      {"ccn_0p05", {-1e10, 1e10}}, // FIXME
      {"ccn_0p1", {-1e10, 1e10}}, // FIXME
      {"ccn_0p2", {-1e10, 1e10}}, // FIXME
      {"ccn_0p5", {-1e10, 1e10}}, // FIXME
      {"ccn_1p0", {-1e10, 1e10}}, // FIXME
      {"hetfrz_immersion_nucleation_tend", {-1e10, 1e10}}, // FIXME
      {"hetfrz_contact_nucleation_tend", {-1e10, 1e10}}, // FIXME
      {"hetfrz_deposition_nucleation_tend", {-1e10, 1e10}}, // FIXME
      // dry deposition
      {"dgnumwet", {-1e10, 1e10}}, // FIXME
    {"fv", {-1e10, 1e10}}, // FIXME
    {"icefrac", {-1e10, 1e10}}, // FIXME
    {"landfrac", {-1e10, 1e10}}, // FIXME
    {"obklen", {-1e10, 1e10}}, // FIXME
    {"ocnfrac", {-1e10, 1e10}}, // FIXME
    {"ram1", {-1e10, 1e10}}, // FIXME
    {"ustar", {-1e10, 1e10}}, // FIXME
    {"wetdens", {-1e10, 1e10}}, // FIXME
    {"deposition_flux_of_cloud_borne_aerosols", {-1e100, 1e100}}, // FIXME
    {"deposition_flux_of_interstitial_aerosols", {-1e100, 1e100}}, // FIXME
    {"fraction_landuse", {-1e100, 1e100}}, // FIXME
    //microphysics
    {"SW_flux_dn", {-1e10, 1e10}}, // FIXME
    {"horiz_winds", {-1e10, 1e10}}, // FIXME
    {"nevapr", {-1e10, 1e10}}, // FIXME
    {"precip_ice_surf_mass", {-1e10, 1e10}}, // FIXME
    {"precip_liq_surf_mass", {-1e10, 1e10}}, // FIXME
    {"precip_total_tend", {-1e10, 1e10}}, // FIXME
    {"ps", {-1e10, 1e10}}, // FIXME
    {"sfc_alb_dir_vis", {-1e10, 1e10}}, // FIXME
    {"snow_depth_land", {-1e10, 1e10}}, // FIXME
    {"surf_radiative_T", {-1e10, 1e10}}, // FIXME
    // optics
    {"pseudo_density_dry", {-1e10, 1e10}}, // FIXME
    {"aero_g_sw", {-1e10, 1e10}}, // FIXME
    {"aero_ssa_sw", {-1e10, 1e10}}, // FIXME
    {"aero_tau_lw", {-1e10, 1e10}}, // FIXME
    {"aero_tau_sw", {-1e10, 1e10}}, // FIXME
    {"aodvis", {-1e10, 1e10}}, // FIXME
    {"sst", {-1e10, 1e10}}, // FIXME
    {"dstflx", {-1e10, 1e10}}, // FIXME
    // srf_and_online_emissions

    //wetscav
    {"drydep_hydrophilic_bc", {-1e10, 1e10}}, // FIXME
    {"drydep_hydrophilic_oc", {-1e10, 1e10}}, // FIXME
    {"wetdep_dust_bin1", {-1e10, 1e10}}, // FIXME
    {"wetdep_dust_bin2", {-1e10, 1e10}}, // FIXME
    {"wetdep_dust_bin3", {-1e10, 1e10}}, // FIXME
    {"wetdep_dust_bin4", {-1e10, 1e10}}, // FIXME
    {"wetdep_hydrophilic_bc", {-1e10, 1e10}}, // FIXME
    {"wetdep_hydrophilic_oc", {-1e10, 1e10}}, // FIXME
    {"precip_total_tend", {-1e10, 1e10}}, // FIXME
    {"aerdepwetcw", {-1e10, 1e10}}, // FIXME
    {"aerdepwetis", {-1e10, 1e10}}, // FIXME
    {"fracis", {-1e10, 1e10}}, // FIXME
    {"qaerwat", {-1e10, 1e10}} // FIXME
      };

    auto it = limits.find(field_name);
    if (it == limits.end()) {
        EKAT_ERROR_MSG("Error!. "<< field_name << " not found in physical_min_max. ");
    }
  return  it->second;
}
inline Real physical_min(const std::string &field_name) {
  return physical_min_max(field_name).first;
}
inline Real physical_max(const std::string &field_name) {
  return physical_min_max(field_name).second;
}

} // namespace mam4

#endif
