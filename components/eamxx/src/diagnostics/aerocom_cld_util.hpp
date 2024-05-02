#ifndef EAMXX_AEROCOMCLD_DIAG_UTIL
#define EAMXX_AEROCOMCLD_DIAG_UTIL

#include <map>

namespace scream {

class AeroComCldDiagUtil {
 public:
  std::map<std::string, int> index_map;
  std::map<std::string, std::string> units_map;
  int size;

  AeroComCldDiagUtil() {
    // T_mid
    index_map["T_mid"] = 0;
    units_map["T_mid"] = "K";
    // p_mid
    index_map["p_mid"] = 1;
    units_map["p_mid"] = "Pa";
    // cldfrac_ice
    index_map["cldfrac_ice"] = 2;
    units_map["cldfrac_ice"] = "nondim";
    // cldfrac_liq
    index_map["cldfrac_liq"] = 3;
    units_map["cldfrac_liq"] = "nondim";
    // cdnc
    index_map["cdnc"] = 4;
    units_map["cdnc"] = "#/m3";
    // eff_radius_qc
    index_map["eff_radius_qc"] = 5;
    units_map["eff_radius_qc"] = "micron";
    // eff_radius_qi
    index_map["eff_radius_qi"] = 6;
    units_map["eff_radius_qi"] = "micron";
    // cldfrac_tot
    index_map["cldfrac_tot"] = 7;
    units_map["cldfrac_tot"] = "nondim";
    // nc
    index_map["nc"] = 8;
    units_map["nc"] = "1/kg";
    // ni
    index_map["ni"] = 9;
    units_map["ni"] = "1/kg";

    size = index_map.size();
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
