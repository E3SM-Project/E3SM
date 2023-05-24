#include "share/util/eamxx_fv_phys_rrtmgp_active_gases_workaround.hpp"

namespace scream {

bool fvphyshack;

void fv_phys_rrtmgp_active_gases_init (const ekat::ParameterList& p) {
  const auto& v = p.sublist("atmosphere_processes").sublist("physics")
    .sublist("rrtmgp").get<std::vector<std::string>>("active_gases");
  if (ekat::contains(v, "o3")) {
    TraceGasesWorkaround::singleton().add_active_gas("o3_volume_mix_ratio");
  }
}

void fv_phys_rrtmgp_active_gases_set_restart (const bool restart) {
  TraceGasesWorkaround::singleton().set_restart(restart);
}

}
