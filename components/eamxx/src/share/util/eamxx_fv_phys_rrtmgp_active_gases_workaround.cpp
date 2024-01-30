#include "share/util/eamxx_fv_phys_rrtmgp_active_gases_workaround.hpp"

namespace scream {

bool fvphyshack = false;

void fv_phys_rrtmgp_active_gases_set_restart (const bool restart) {
  TraceGasesWorkaround::singleton().set_restart(restart);
}

}
