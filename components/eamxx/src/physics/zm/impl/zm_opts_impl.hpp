#ifndef ZM_ZM_OPTS_IMPL_HPP
#define ZM_ZM_OPTS_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_opts_init. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
void Functions<S,D>::zm_opts_init()
{
  s_zm_opts.use_fortran_bridge  = false;
  s_zm_opts.apply_detr_tend     = true;
  s_zm_opts.upper_limit_pref    = 40e2;
  s_zm_opts.plenest             = static_cast<Int>(ZMC::tmax-ZMC::tmin) + 3;
  s_zm_opts.tau                 = 3600;
  s_zm_opts.alfa                = ZMC::alfa;
  s_zm_opts.ke                  = ZMC::ke;
  s_zm_opts.dmpdz               = ZMC::dmpdz;
  s_zm_opts.tpert_fix           = true;
  s_zm_opts.tpert_fac           = 2;
  s_zm_opts.tiedke_add          = ZMC::tiedke_add;
  s_zm_opts.c0_lnd              = ZMC::c0;
  s_zm_opts.c0_ocn              = ZMC::c0;
  s_zm_opts.num_cin             = 1;
  s_zm_opts.limcnv              = 23; // note - default for E3SMv3 => ne30pg2 w/ L80
  s_zm_opts.mx_bot_lyr_adj      = 1;
  s_zm_opts.trig_dcape          = true;
  s_zm_opts.trig_ull            = true;
  s_zm_opts.clos_dyn_adj        = true;
  s_zm_opts.no_deep_pbl         = false;
  // ZM micro parameters
  s_zm_opts.zm_microp           = false;
  s_zm_opts.old_snow            = true;
  // MCSP parameters
  s_zm_opts.mcsp_enabled        = true;
  s_zm_opts.mcsp_t_coeff        = ZMC::MCSP_t_coeff_default;
  s_zm_opts.mcsp_q_coeff        = ZMC::MCSP_q_coeff_default;
  s_zm_opts.mcsp_u_coeff        = ZMC::MCSP_u_coeff_default;
  s_zm_opts.mcsp_v_coeff        = ZMC::MCSP_v_coeff_default;

  //
  // set up table values of saturation vapor pressure
  //

  // Add two to make the table slightly too big, just in case.
  const Int plenest = static_cast<Int>(ZMC::tmax-ZMC::tmin) + 3;

  // Allocate SVP table.
  view_1d<Real> estbl("estbl", plenest);
  Kokkos::parallel_for("zm_calculate_estbl", Kokkos::RangePolicy<typename KT::ExeSpace>(0, plenest), KOKKOS_LAMBDA(const int i) {
    estbl(i) = svp_trans(ZMC::tmin + i);
  });
  s_zm_opts.estbl = estbl;
}

} // namespace zm
} // namespace scream

#endif
