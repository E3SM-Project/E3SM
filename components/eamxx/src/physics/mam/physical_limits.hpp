// mam4xx: Copyright (c) 2022,
// Battelle Memorial Institute and
// National Technology & Engineering Solutions of Sandia, LLC (NTESS)
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MAM_PHYSICAL_LIMITS_HPP
#define MAM_PHYSICAL_LIMITS_HPP
#include <map>
#include <share/atm_process/atmosphere_process.hpp>
#include <share/property_checks/field_within_interval_check.hpp>
#include <string>
#include <utility>
// For MAM4 aerosol configuration

namespace scream::mam_coupling {

inline const std::pair<Real, Real> &physical_min_max(
    const std::string &field_name) {
  static const std::map<std::string, std::pair<Real, Real>> limits = {
      {"nmr", {0, 1e13}}, {"mmr", {-1e-10, 1e-2}}};

  auto it = limits.find(field_name);
  if(it == limits.end()) {
    // NOTE: If we do not find a variable name in physical_min_max,
    // we return a pair (-1, -1) that will bypass the interval check.
    static const std::pair<Real, Real> do_not_check = std::make_pair(-1, -1);
    return do_not_check;
  } else {
    return it->second;
  }
}
inline Real physical_min(const std::string &field_name) {
  return physical_min_max(field_name).first;
}
inline Real physical_max(const std::string &field_name) {
  return physical_min_max(field_name).second;
}

}  // namespace scream::mam_coupling

#endif
