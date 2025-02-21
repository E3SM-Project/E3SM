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
      {"constituent_fluxes", {0, 1e10}}, // FIXME
      {"phis", {0, 1e10}} // FIXME
      };
  return limits.at(field_name);
}
inline Real physical_min(const std::string &field_name) {
  return physical_min_max(field_name).first;
}
inline Real physical_max(const std::string &field_name) {
  return physical_min_max(field_name).second;
}
} // namespace mam4

#endif
