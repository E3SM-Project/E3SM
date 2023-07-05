#include "dp_functions_f90.hpp"

#include "dp_f90.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

#include "share/util/scream_deep_copy.hpp"

#include <random>

using scream::Real;
using scream::Int;

using scream::dp::element_t;
using scream::dp::hvcoord_t;
using scream::dp::hybrid_t;
using scream::dp::timelevel_t;

//
// A C interface to DP fortran calls. The stubs below will link to fortran definitions in dp_iso_c.f90
//

extern "C" {

void advance_iop_forcing_c(Int plev, Int pcnst, Real scm_dt, Real ps_in, Real* u_in, Real* v_in, Real* t_in, Real* q_in, Real* t_phys_frc, Real* u_update, Real* v_update, Real* t_update, Real* q_update);
void advance_iop_nudging_c(Int plev, Real scm_dt, Real ps_in, Real* t_in, Real* q_in, Real* t_update, Real* q_update, Real* relaxt, Real* relaxq);
void advance_iop_subsidence_c(Int plev, Int pcnst, Real scm_dt, Real ps_in, Real* u_in, Real* v_in, Real* t_in, Real* q_in, Real* u_update, Real* v_update, Real* t_update, Real* q_update);
void iop_setinitial_c(Int nelemd, element_t* elem);
void iop_broadcast_c();
void apply_iop_forcing_c(Int nelemd, element_t* elem, hvcoord_t* hvcoord, hybrid_t* hybrid, timelevel_t* tl, Int n, bool t_before_advance, Int nets, Int nete);
void iop_domain_relaxation_c(Int nelemd, Int np, Int nlev, element_t* elem, hvcoord_t hvcoord, hybrid_t hybrid, Int t1, Real* dp, Int nelemd_todo, Int np_todo, Real dt);
void crm_resolved_turb_c(Int nelemd, element_t* elem, hvcoord_t hvcoord, hybrid_t hybrid, Int t1, Int nelemd_todo, Int np_todo);
void iop_default_opts_c(Real* scmlat_out, Real* scmlon_out, char** iopfile_out, bool* single_column_out, bool* scm_iop_srf_prop_out, bool* iop_nudge_tq_out, bool* iop_nudge_uv_out, Real* iop_nudge_tq_low_out, Real* iop_nudge_tq_high_out, Real* iop_nudge_tscale_out, bool* scm_observed_aero_out, bool* iop_dosubsidence_out, bool* scm_multcols_out, bool* dp_crm_out, Real* iop_perturb_high_out, bool* precip_off_out, bool* scm_zero_non_iop_tracers_out);
} // extern "C" : end _c decls

namespace scream {
namespace dp {

//
// Glue functions to call fortran from from C++ with the Data struct
//

void advance_iop_forcing(AdvanceIopForcingData& d)
{
  dp_init(d.plev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  advance_iop_forcing_c(d.plev, d.pcnst, d.scm_dt, d.ps_in, d.u_in, d.v_in, d.t_in, d.q_in, d.t_phys_frc, d.u_update, d.v_update, d.t_update, d.q_update);
  d.transpose<ekat::TransposeDirection::f2c>();
}


void advance_iop_nudging(AdvanceIopNudgingData& d)
{
  dp_init(d.plev, true);
  advance_iop_nudging_c(d.plev, d.scm_dt, d.ps_in, d.t_in, d.q_in, d.t_update, d.q_update, d.relaxt, d.relaxq);
}

void advance_iop_subsidence(AdvanceIopSubsidenceData& d)
{
  dp_init(d.plev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  advance_iop_subsidence_c(d.plev, d.pcnst, d.scm_dt, d.ps_in, d.u_in, d.v_in, d.t_in, d.q_in, d.u_update, d.v_update, d.t_update, d.q_update);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void iop_setinitial(IopSetinitialData& d)
{
  dp_init(d.plev, true);
  iop_setinitial_c(d.nelemd, d.elem);
}

void iop_broadcast(IopBroadcastData& d)
{
  dp_init(d.plev, true);
  iop_broadcast_c();
}

void apply_iop_forcing(ApplyIopForcingData& d)
{
  dp_init(d.plev, true);
  apply_iop_forcing_c(d.nelemd, d.elem, &d.hvcoord, &d.hybrid, &d.tl, d.n, d.t_before_advance, d.nets, d.nete);
}

void iop_domain_relaxation(IopDomainRelaxationData& d)
{
  dp_init(d.nlev, true);
  d.transpose<ekat::TransposeDirection::c2f>();
  iop_domain_relaxation_c(d.nelemd, d.np, d.nlev, d.elem, d.hvcoord, d.hybrid, d.t1, d.dp, d.nelemd_todo, d.np_todo, d.dt);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void crm_resolved_turb(CrmResolvedTurbData& d)
{
  dp_init(d.plev, true);
  crm_resolved_turb_c(d.nelemd, d.elem, d.hvcoord, d.hybrid, d.t1, d.nelemd_todo, d.np_todo);
}

void iop_default_opts(IopDefaultOptsData& d)
{
  dp_init(d.plev, true);
  char cbuff[512];
  char* buffptr = cbuff;
  iop_default_opts_c(&d.scmlat_out, &d.scmlon_out, &buffptr, &d.single_column_out, &d.scm_iop_srf_prop_out, &d.iop_nudge_tq_out, &d.iop_nudge_uv_out, &d.iop_nudge_tq_low_out, &d.iop_nudge_tq_high_out, &d.iop_nudge_tscale_out, &d.scm_observed_aero_out, &d.iop_dosubsidence_out, &d.scm_multcols_out, &d.dp_crm_out, &d.iop_perturb_high_out, &d.precip_off_out, &d.scm_zero_non_iop_tracers_out);
  d.iopfile_out = std::string(buffptr);
}

// end _c impls

//
// _f function definitions. These expect data in C layout
//

void advance_iop_forcing_f(Int plev, Int pcnst, Real scm_dt, Real ps_in, Real* u_in, Real* v_in, Real* t_in, Real* q_in, Real* t_phys_frc, Real* u_update, Real* v_update, Real* t_update, Real* q_update)
{
  // TODO
}
void advance_iop_nudging_f(Int plev, Real scm_dt, Real ps_in, Real* t_in, Real* q_in, Real* t_update, Real* q_update, Real* relaxt, Real* relaxq)
{
  // TODO
}
void advance_iop_subsidence_f(Int plev, Int pcnst, Real scm_dt, Real ps_in, Real* u_in, Real* v_in, Real* t_in, Real* q_in, Real* u_update, Real* v_update, Real* t_update, Real* q_update)
{
  // TODO
}
void iop_setinitial_f(Int nelemd, element_t* elem)
{
  // TODO
}
void iop_broadcast_f()
{
#if 0
  using PF = Functions<Real, DefaultDevice>;

  using Spack   = typename PF::Spack;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    PF::iop_broadcast();
  });
#endif

}
void apply_iop_forcing_f(Int nelemd, element_t* elem, hvcoord_t* hvcoord, hybrid_t hybrid, timelevel_t tl, Int n, bool t_before_advance, Int nets, Int nete)
{
  // TODO
}
void iop_domain_relaxation_f(Int nelemd, Int np, Int nlev, element_t* elem, hvcoord_t hvcoord, hybrid_t hybrid, Int t1, Real* dp, Int nelemd_todo, Int np_todo, Real dt)
{
  // TODO
}
void crm_resolved_turb_f(Int nelemd, element_t* elem, hvcoord_t hvcoord, hybrid_t hybrid, Int t1, Int nelemd_todo, Int np_todo)
{
  // TODO
}
void iop_default_opts_f(Real* scmlat_out, Real* scmlon_out, char** iopfile_out, bool* single_column_out, bool* scm_iop_srf_prop_out, bool* iop_nudge_tq_out, bool* iop_nudge_uv_out, Real* iop_nudge_tq_low_out, Real* iop_nudge_tq_high_out, Real* iop_nudge_tscale_out, bool* scm_observed_aero_out, bool* iop_dosubsidence_out, bool* scm_multcols_out, bool* dp_crm_out, Real* iop_perturb_high_out, bool* precip_off_out, bool* scm_zero_non_iop_tracers_out)
{
#if 0
  using PF = Functions<Real, DefaultDevice>;

  using Spack   = typename PF::Spack;
  using view_1d = typename PF::view_1d<Real>;
  using bview_1d = typename PF::view_1d<bool>;

  view_1d t_d("t_d", 6);
  const auto t_h = Kokkos::create_mirror_view(t_d);

  bview_1d bt_d("bt_d", 10);
  const auto bt_h = Kokkos::create_mirror_view(bt_d);

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack iop_nudge_tq_high_out_(), iop_nudge_tq_low_out_(), iop_nudge_tscale_out_(), iop_perturb_high_out_(), scmlat_out_(), scmlon_out_();
    bool dp_crm_out_(), iop_dosubsidence_out_(), iop_nudge_tq_out_(), iop_nudge_uv_out_(), precip_off_out_(), scm_iop_srf_prop_out_(), scm_multcols_out_(), scm_observed_aero_out_(), scm_zero_non_iop_tracers_out_(), single_column_out_();
    PF::iop_default_opts(scmlat_out_, scmlon_out_, iopfile_out_, single_column_out_, scm_iop_srf_prop_out_, iop_nudge_tq_out_, iop_nudge_uv_out_, iop_nudge_tq_low_out_, iop_nudge_tq_high_out_, iop_nudge_tscale_out_, scm_observed_aero_out_, iop_dosubsidence_out_, scm_multcols_out_, dp_crm_out_, iop_perturb_high_out_, precip_off_out_, scm_zero_non_iop_tracers_out_);
    t_d(0) = iop_nudge_tq_high_out_[0];
    t_d(1) = iop_nudge_tq_low_out_[0];
    t_d(2) = iop_nudge_tscale_out_[0];
    t_d(3) = iop_perturb_high_out_[0];
    t_d(4) = scmlat_out_[0];
    t_d(5) = scmlon_out_[0];
    bt_d(0) = dp_crm_out_;
    bt_d(1) = iop_dosubsidence_out_;
    bt_d(2) = iop_nudge_tq_out_;
    bt_d(3) = iop_nudge_uv_out_;
    bt_d(4) = precip_off_out_;
    bt_d(5) = scm_iop_srf_prop_out_;
    bt_d(6) = scm_multcols_out_;
    bt_d(7) = scm_observed_aero_out_;
    bt_d(8) = scm_zero_non_iop_tracers_out_;
    bt_d(9) = single_column_out_;
  });
  Kokkos::deep_copy(t_h, t_d);
  Kokkos::deep_copy(bt_h, bt_d);
  *iop_nudge_tq_high_out = t_h(0);
  *iop_nudge_tq_low_out = t_h(1);
  *iop_nudge_tscale_out = t_h(2);
  *iop_perturb_high_out = t_h(3);
  *scmlat_out = t_h(4);
  *scmlon_out = t_h(5);
  *dp_crm_out = bt_h(0);
  *iop_dosubsidence_out = bt_h(1);
  *iop_nudge_tq_out = bt_h(2);
  *iop_nudge_uv_out = bt_h(3);
  *precip_off_out = bt_h(4);
  *scm_iop_srf_prop_out = bt_h(5);
  *scm_multcols_out = bt_h(6);
  *scm_observed_aero_out = bt_h(7);
  *scm_zero_non_iop_tracers_out = bt_h(8);
  *single_column_out = bt_h(9);
#endif

}
} // namespace dp
} // namespace scream
