#ifndef DP_IOP_SETOPTS_IMPL_HPP
#define DP_IOP_SETOPTS_IMPL_HPP

#include "dp_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace dp {

/*
 * Implementation of dp iop_setopts. Clients should NOT
 * #include this file, but include dp_functions.hpp instead.
 */

template<typename S, typename D>
void Functions<S,D>::iop_setopts(const Spack& scmlat_in, const Spack& scmlon_in, const std::string& iopfile_in, const bool& single_column_in, const bool& scm_iop_srf_prop_in, const bool& iop_nudge_tq_in, const bool& iop_nudge_uv_in, const Spack& iop_nudge_tq_low_in, const Spack& iop_nudge_tq_high_in, const Spack& iop_nudge_tscale_in, const bool& scm_observed_aero_in, const bool& iop_dosubsidence_in, const bool& scm_multcols_in, const bool& dp_crm_in, const Spack& iop_perturb_high_in, const bool& precip_off_in, const bool& scm_zero_non_iop_tracers_in)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace dp
} // namespace scream

#endif
