#include "ekat/ekat_parameter_list.hpp"

namespace scream {
extern bool fvphyshack;
void fv_phys_rrtmgp_active_gases_init(const ekat::ParameterList& p);
void fv_phys_rrtmgp_active_gases_set_restart(const bool restart);
}
