#ifndef EAMXX_PBL_ENTRAINMENT_BUDGET_DIAG_UTIL
#define EAMXX_PBL_ENTRAINMENT_BUDGET_DIAG_UTIL

#include <map>

namespace scream {
class PBLEntrainmentBudgetDiagUtil {
 public:
  std::map<std::string, int> index_map;
  std::map<std::string, std::string> units_map;
  int size;
  std::string pblinvalg;

  PBLEntrainmentBudgetDiagUtil() {
    // Start (post) incrementing size from 0
    size = 0;

    // pblp
    index_map["p+"] = size++;
    units_map["p+"] = "Pa";
    // tl+
    index_map["tl+"] = size++;
    units_map["tl+"] = "K";
    // tl^
    index_map["tl^"] = size++;
    units_map["tl^"] = "K";
    // tl_ttend
    index_map["tl_ttend"] = size++;
    units_map["tl_ttend"] = "K/s";
    // qt+
    index_map["qt+"] = size++;
    units_map["qt+"] = "kg/kg";
    // qt^
    index_map["qt^"] = size++;
    units_map["qt^"] = "kg/kg";
    // qt_ttend
    index_map["qt_ttend"] = size++;
    units_map["qt_ttend"] = "kg/kg/s";
    // dF
    index_map["dF"] = size++;
    units_map["dF"] = "W/m^2";
    // pblinvalg
    pblinvalg = "temperature-inversion";
  }
};
}  // namespace scream

#endif  // EAMXX_PBL_ENTRAINMENT_BUDGET_DIAG_UTIL
