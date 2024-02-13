#ifndef DP_IOP_DEFAULT_OPTS_IMPL_HPP
#define DP_IOP_DEFAULT_OPTS_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace dp {

/*
 * Implementation of dp iop_default_opts. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
void Functions<S,D>::iop_default_opts(Spack& scmlat_out, Spack& scmlon_out, std::string& iopfile_out, bool& single_column_out, bool& scm_iop_srf_prop_out, bool& iop_nudge_tq_out, bool& iop_nudge_uv_out, Spack& iop_nudge_tq_low_out, Spack& iop_nudge_tq_high_out, Spack& iop_nudge_tscale_out, bool& scm_observed_aero_out, bool& iop_dosubsidence_out, bool& scm_multcols_out, bool& dp_crm_out, Spack& iop_perturb_high_out, bool& precip_off_out, bool& scm_zero_non_iop_tracers_out)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace dp
} // namespace scream

#endif
