#ifndef EAMXX_AEROCOMCLD_DIAG_UTIL
#define EAMXX_AEROCOMCLD_DIAG_UTIL

#include <map>

namespace scream {

class AeroComCldDiagUtil {
 public:
  std::map<std::string, int> index_map;
  std::map<std::string, std::string> units_map;
  unsigned int size;

  AeroComCldDiagUtil() {
    // Start (post) incrementing size from 0
    size = 0;

    // T_mid
    index_map["T_mid"] = size++;
    units_map["T_mid"] = "K";
    // p_mid
    index_map["p_mid"] = size++;
    units_map["p_mid"] = "Pa";
    // cldfrac_ice
    index_map["cldfrac_ice"] = size++;
    units_map["cldfrac_ice"] = "nondim";
    // cldfrac_liq
    index_map["cldfrac_liq"] = size++;
    units_map["cldfrac_liq"] = "nondim";
    // cdnc
    index_map["cdnc"] = size++;
    units_map["cdnc"] = "#/m3";
    // eff_radius_qc
    index_map["eff_radius_qc"] = size++;
    units_map["eff_radius_qc"] = "micron";
    // eff_radius_qi
    index_map["eff_radius_qi"] = size++;
    units_map["eff_radius_qi"] = "micron";
    // cldfrac_tot
    index_map["cldfrac_tot"] = size++;
    units_map["cldfrac_tot"] = "nondim";
    // nc
    index_map["nc"] = size++;
    units_map["nc"] = "1/kg";
    // ni
    index_map["ni"] = size++;
    units_map["ni"] = "1/kg";
  }
};

constexpr double q_thresh_set() {
  // The total q is checked at the outset of the aerocom_cld routine
  // to determine whether or not there is any cloud. This function sets
  // the arbitrary q threshold to be used.
  return 0.0;
}

constexpr double cldfrac_tot_thresh_set() {
  // The total cldfrac is checked at the outset of the aerocom_cld routine
  // to determine whether or not there is any cloud. This function sets
  // the arbitrary cldfrac threshold to be used, which is slightly more than
  // 0.0 to guard against potential issues near zero.
  return 0.001;
}

}  // namespace scream

#endif  // EAMXX_AEROCOMCLD_DIAG_UTIL
