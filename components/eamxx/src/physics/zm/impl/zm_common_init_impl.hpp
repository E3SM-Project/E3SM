#ifndef ZM_ZM_COMMON_INIT_IMPL_HPP
#define ZM_ZM_COMMON_INIT_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_common_init. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
void Functions<S,D>::zm_common_init()
{
  s_common_init.tau             = 3600.0;
  s_common_init.alfa            = 0.14;
  s_common_init.ke              = 2.5E-6;
  s_common_init.dmpdz           = -0.7e-3;
  s_common_init.tpert_fix       = true;
  s_common_init.tpert_fac       = 2.0;
  s_common_init.tiedke_add      = 0.8;
  s_common_init.c0_lnd          = 0.0020;
  s_common_init.c0_ocn          = 0.0020;
  s_common_init.num_cin         = 1;
  s_common_init.limcnv          = 23; // note - default for E3SMv3 => ne30pg2 w/ L80
  s_common_init.mx_bot_lyr_adj  = 1;
  s_common_init.trig_dcape      = true;
  s_common_init.trig_ull        = true;
  s_common_init.clos_dyn_adj    = true;
  s_common_init.no_deep_pbl     = false;
  // ZM micro parameters
  s_common_init.zm_microp       = true;
  s_common_init.old_snow        = false;
  s_common_init.auto_fac        = 7.0;
  s_common_init.accr_fac        = 1.5;
  s_common_init.micro_dcs       = 150.E-6;
  // MCSP parameters
  s_common_init.mcsp_enabled    = true;
  s_common_init.mcsp_t_coeff    = 0.3;
  s_common_init.mcsp_q_coeff    = 0;
  s_common_init.mcsp_u_coeff    = 0;
  s_common_init.mcsp_v_coeff    = 0;
}

} // namespace zm
} // namespace scream

#endif
