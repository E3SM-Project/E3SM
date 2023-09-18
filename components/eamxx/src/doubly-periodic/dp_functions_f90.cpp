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
void iop_setopts_c(Real scmlat_in, Real scmlon_in, const char** iopfile_in, bool single_column_in, bool scm_iop_srf_prop_in, bool iop_nudge_tq_in, bool iop_nudge_uv_in, Real iop_nudge_tq_low_in, Real iop_nudge_tq_high_in, Real iop_nudge_tscale_in, bool scm_observed_aero_in, bool iop_dosubsidence_in, bool scm_multcols_in, bool dp_crm_in, Real iop_perturb_high_in, bool precip_off_in, bool scm_zero_non_iop_tracers_in);
void setiopupdate_init_c();
void setiopupdate_c();
void readiopdata_c(Int plev, bool iop_update_phase1, Real* hyam, Real* hybm);
void iop_intht_c();
} // extern "C" : end _c decls

namespace scream {
namespace dp {

//
// Glue functions to call fortran from from C++ with the Data struct
//

void advance_iop_forcing(AdvanceIopForcingData& d)
{
  dp_init();
  d.transpose<ekat::TransposeDirection::c2f>();
  advance_iop_forcing_c(d.plev, d.pcnst, d.scm_dt, d.ps_in, d.u_in, d.v_in, d.t_in, d.q_in, d.t_phys_frc, d.u_update, d.v_update, d.t_update, d.q_update);
  d.transpose<ekat::TransposeDirection::f2c>();
}


void advance_iop_nudging(AdvanceIopNudgingData& d)
{
  dp_init();
  advance_iop_nudging_c(d.plev, d.scm_dt, d.ps_in, d.t_in, d.q_in, d.t_update, d.q_update, d.relaxt, d.relaxq);
}

void advance_iop_subsidence(AdvanceIopSubsidenceData& d)
{
  dp_init();
  d.transpose<ekat::TransposeDirection::c2f>();
  advance_iop_subsidence_c(d.plev, d.pcnst, d.scm_dt, d.ps_in, d.u_in, d.v_in, d.t_in, d.q_in, d.u_update, d.v_update, d.t_update, d.q_update);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void iop_setinitial(IopSetinitialData& d)
{
  dp_init();
  //iop_setinitial_c(d.nelemd, d.elem);
}

void iop_broadcast(IopBroadcastData& d)
{
  dp_init();
  iop_broadcast_c();
}

void apply_iop_forcing(ApplyIopForcingData& d)
{
  dp_init();
  apply_iop_forcing_c(d.nelemd, d.elem, &d.hvcoord, &d.hybrid, &d.tl, d.n, d.t_before_advance, d.nets, d.nete);
}

void iop_domain_relaxation(IopDomainRelaxationData& d)
{
  dp_init();
  d.transpose<ekat::TransposeDirection::c2f>();
  iop_domain_relaxation_c(d.nelemd, d.np, d.nlev, d.elem, d.hvcoord, d.hybrid, d.t1, d.dp, d.nelemd_todo, d.np_todo, d.dt);
  d.transpose<ekat::TransposeDirection::f2c>();
}

void crm_resolved_turb(CrmResolvedTurbData& d)
{
  dp_init();
  crm_resolved_turb_c(d.nelemd, d.elem, d.hvcoord, d.hybrid, d.t1, d.nelemd_todo, d.np_todo);
}

void iop_default_opts(IopDefaultOptsData& d)
{
  dp_init();
  char cbuff[512] = "";
  char* buffptr = cbuff;
  iop_default_opts_c(&d.scmlat_out, &d.scmlon_out, &buffptr, &d.single_column_out, &d.scm_iop_srf_prop_out, &d.iop_nudge_tq_out, &d.iop_nudge_uv_out, &d.iop_nudge_tq_low_out, &d.iop_nudge_tq_high_out, &d.iop_nudge_tscale_out, &d.scm_observed_aero_out, &d.iop_dosubsidence_out, &d.scm_multcols_out, &d.dp_crm_out, &d.iop_perturb_high_out, &d.precip_off_out, &d.scm_zero_non_iop_tracers_out);
  d.iopfile_out = std::string(buffptr);
}

void iop_setopts(IopSetoptsData& d)
{
  dp_init();
  const char* cptr = d.iopfile_in.c_str();
  iop_setopts_c(d.scmlat_in, d.scmlon_in, &cptr, d.single_column_in, d.scm_iop_srf_prop_in, d.iop_nudge_tq_in, d.iop_nudge_uv_in, d.iop_nudge_tq_low_in, d.iop_nudge_tq_high_in, d.iop_nudge_tscale_in, d.scm_observed_aero_in, d.iop_dosubsidence_in, d.scm_multcols_in, d.dp_crm_in, d.iop_perturb_high_in, d.precip_off_in, d.scm_zero_non_iop_tracers_in);
}

void setiopupdate_init(SetiopupdateInitData& d)
{
  dp_init();
  setiopupdate_init_c();
}

void setiopupdate(SetiopupdateData& d)
{
  dp_init();
  setiopupdate_c();
}

void readiopdata(ReadiopdataData& d)
{
  dp_init();
  readiopdata_c(d.plev, d.iop_update_phase1, d.hyam, d.hybm);
}

void iop_intht(IopInthtData& d)
{
  dp_init();
  iop_intht_c();
}

// end _c impls

//
// _f function definitions. These expect data in C layout
//

void advance_iop_forcing_f(Int plev, Int pcnst, Real scm_dt, Real ps_in, bool have_u, bool have_v, bool dp_crm, bool use_3dfrc, Real* u_in, Real* v_in, Real* t_in, Real* q_in, Real* t_phys_frc, Real* divt3d, Real* divq3d, Real* divt, Real* divq, Real* wfld, Real* uobs, Real* vobs, Real* hyai, Real* hyam, Real* hybi, Real* hybm, Real* u_update, Real* v_update, Real* t_update, Real* q_update)
{
  using DPF  = Functions<Real, DefaultDevice>;

  using Spack = typename DPF::Spack;
  using view_1d = typename DPF::view_1d<Spack>;
  using view_2d = typename DPF::view_2d<Spack>;
  using KT = typename DPF::KT;
  using ExeSpace = typename KT::ExeSpace;
  using MemberType = typename DPF::MemberType;

  // Some of the workspaces need plev+1 items
  const Int plev_pack = ekat::npack<Spack>(plev);
  const Int plevp_pack = ekat::npack<Spack>(plev+1);

  // Set up views
  std::vector<view_1d> temp_d(AdvanceIopForcingData::NUM_ARRAYS-4);
  std::vector<view_2d> temp_2d_d(4);

  ekat::host_to_device({u_in, v_in, t_in, t_phys_frc, divt3d, divt, wfld, uobs, vobs, hyai, hyam, hybi, hybm, u_update, v_update, t_update},
                       plev, temp_d);

  ekat::host_to_device({ q_in, divq3d, divq, q_update },
                       pcnst, plev, temp_2d_d, true);

  view_1d
    u_in_d       (temp_d[0]),
    v_in_d       (temp_d[1]),
    t_in_d       (temp_d[2]),
    t_phys_frc_d (temp_d[3]),
    divt3d_d     (temp_d[4]),
    divt_d       (temp_d[5]),
    wfld_d       (temp_d[6]),
    uobs_d       (temp_d[7]),
    vobs_d       (temp_d[8]),
    hyai_d       (temp_d[9]),
    hyam_d       (temp_d[10]),
    hybi_d       (temp_d[11]),
    hybm_d       (temp_d[12]),
    u_update_d   (temp_d[13]),
    v_update_d   (temp_d[14]),
    t_update_d   (temp_d[15]);

  view_2d
    q_in_d    (temp_2d_d[0]),
    divq3d_d  (temp_2d_d[1]),
    divq_d    (temp_2d_d[2]),
    q_update_d(temp_2d_d[3]);

  // Call core function from kernel
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, plev_pack);
  ekat::WorkspaceManager<Spack> wsm(plevp_pack, 3, policy);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

      DPF::advance_iop_forcing(
        plev, pcnst, have_u, have_v, dp_crm, use_3dfrc, scm_dt, ps_in,
        u_in_d, v_in_d, t_in_d, q_in_d, t_phys_frc_d, divt3d_d, divq3d_d, divt_d, divq_d, wfld_d, uobs_d, vobs_d, hyai_d, hyam_d, hybi_d, hybm_d,
        team, wsm.get_workspace(team),
        u_update_d, v_update_d, t_update_d, q_update_d);
  });

  // Sync back to host
  std::vector<view_1d> inout_views    = {t_update_d, u_update_d, v_update_d};
  std::vector<view_2d> inout_views_2d = {q_update_d};

  ekat::device_to_host({t_update, u_update, v_update}, plev, inout_views);
  ekat::device_to_host({q_update}, pcnst, plev, inout_views_2d, true);
}

void advance_iop_nudging_f(Int plev, Real scm_dt, Real ps_in, Real* t_in, Real* q_in, Real* tobs, Real* qobs,
                           Real* hyai, Real* hyam, Real* hybi, Real* hybm,
                           Real* t_update, Real* q_update, Real* relaxt, Real* relaxq)
{
  using DPF  = Functions<Real, DefaultDevice>;

  using Spack = typename DPF::Spack;
  using view_1d = typename DPF::view_1d<Spack>;
  using KT = typename DPF::KT;
  using ExeSpace = typename KT::ExeSpace;
  using MemberType = typename DPF::MemberType;

  const Int plev_pack = ekat::npack<Spack>(plev);
  const Int plevp_pack = ekat::npack<Spack>(plev+1);

  // Set up views
  std::vector<view_1d> temp_d(AdvanceIopNudgingData::NUM_ARRAYS);

  ekat::host_to_device({t_in, q_in, tobs, qobs, hyai, hyam, hybi, hybm, t_update, q_update, relaxt, relaxq},
                       plev, temp_d);

  int counter=0;
  view_1d
    t_in_d    (temp_d[counter++]),
    q_in_d    (temp_d[counter++]),
    tobs_d    (temp_d[counter++]),
    qobs_d    (temp_d[counter++]),
    hyai_d    (temp_d[counter++]),
    hyam_d    (temp_d[counter++]),
    hybi_d    (temp_d[counter++]),
    hybm_d    (temp_d[counter++]),
    t_update_d(temp_d[counter++]),
    q_update_d(temp_d[counter++]),
    relaxt_d  (temp_d[counter++]),
    relaxq_d  (temp_d[counter++]);

  // Call core function from kernel
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, plev_pack);
  ekat::WorkspaceManager<Spack> wsm(plevp_pack, 4, policy);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

      DPF::advance_iop_nudging(
        plev, scm_dt, ps_in, t_in_d, q_in_d, tobs_d, qobs_d,
        hyai_d, hyam_d, hybi_d, hybm_d,
        team, wsm.get_workspace(team),
        t_update_d, q_update_d, relaxt_d, relaxq_d);
  });

  // Sync back to host
  std::vector<view_1d> out_views = {t_update_d, q_update_d, relaxt_d, relaxq_d};

  ekat::device_to_host({t_update, q_update, relaxt, relaxq}, plev, out_views);
}

void advance_iop_subsidence_f(Int plev, Int pcnst, Real scm_dt, Real ps_in, Real* u_in, Real* v_in, Real* t_in, Real* q_in, Real* hyai, Real* hyam, Real* hybi, Real* hybm, Real* wfld, Real* u_update, Real* v_update, Real* t_update, Real* q_update)
{
  using DPF  = Functions<Real, DefaultDevice>;

  using Spack = typename DPF::Spack;
  using view_1d = typename DPF::view_1d<Spack>;
  using view_2d = typename DPF::view_2d<Spack>;
  using KT = typename DPF::KT;
  using ExeSpace = typename KT::ExeSpace;
  using MemberType = typename DPF::MemberType;

  // Some of the workspaces need plev+1 items
  const Int plev_pack = ekat::npack<Spack>(plev);
  const Int plevp_pack = ekat::npack<Spack>(plev+1);

  // Set up views
  std::vector<view_1d> temp_d(AdvanceIopSubsidenceData::NUM_ARRAYS-2);
  std::vector<view_2d> temp_2d_d(2);

  ekat::host_to_device({u_in, v_in, t_in, hyai, hyam, hybi, hybm, wfld, u_update, v_update, t_update},
                       plev, temp_d);

  ekat::host_to_device({ q_in, q_update },
                       pcnst, plev, temp_2d_d, true);

  view_1d
    u_in_d       (temp_d[0]),
    v_in_d       (temp_d[1]),
    t_in_d       (temp_d[2]),
    hyai_d       (temp_d[3]),
    hyam_d       (temp_d[4]),
    hybi_d       (temp_d[5]),
    hybm_d       (temp_d[6]),
    wfld_d       (temp_d[7]),
    u_update_d   (temp_d[8]),
    v_update_d   (temp_d[9]),
    t_update_d   (temp_d[10]);

  view_2d
    q_in_d    (temp_2d_d[0]),
    q_update_d(temp_2d_d[1]);

  // Call core function from kernel
  auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, plev_pack);
  ekat::WorkspaceManager<Spack> wsm(plevp_pack, 4, policy);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {

      DPF::advance_iop_subsidence(
        plev, pcnst, scm_dt, ps_in,
        u_in_d, v_in_d, t_in_d, q_in_d, hyai_d, hyam_d, hybi_d, hybm_d, wfld_d,
        team, wsm.get_workspace(team),
        u_update_d, v_update_d, t_update_d, q_update_d);
  });

  // Sync back to host
  std::vector<view_1d> inout_views    = {t_update_d, u_update_d, v_update_d};
  std::vector<view_2d> inout_views_2d = {q_update_d};

  ekat::device_to_host({t_update, u_update, v_update}, plev, inout_views);
  ekat::device_to_host({q_update}, pcnst, plev, inout_views_2d, true);
}

void iop_setinitial_f(Int plev, Int pcnst, Int nelemd, Int np, Int nstep, Real psobs, bool use_replay, bool dynproc, bool have_t, bool have_q, bool have_ps, bool have_u, bool have_v, bool have_numliq, bool have_cldliq, bool have_numice, bool have_cldice, bool scm_zero_non_iop_tracers, bool is_first_restart_step, Real* qmin, Real* uobs, Real* vobs, Real* numliqobs, Real* numiceobs, Real* cldliqobs, Real* cldiceobs, Real* dx_short, tracer_t* tracers, element_t* elem, Real* dyn_dx_size, Real* tobs, Real* qobs)
{
  using DPF  = Functions<Real, DefaultDevice>;

  using Spack   = typename DPF::Spack;
  using Scalarp = ekat::Pack<Real, 1>;
  using view_1d =  typename DPF::view_1d<Spack>;
  using view_1ds = typename DPF::view_1d<Scalarp>;
  using sview_1ds = typename DPF::view_1d<Real>;
  using KT = typename DPF::KT;
  using ExeSpace = typename KT::ExeSpace;
  using MemberType = typename DPF::MemberType;

  // Some of the workspaces need plev+1 items
  const Int plev_pack = ekat::npack<Spack>(plev);

  // Set up views
  std::vector<view_1d>  temp_d(8);
  std::vector<view_1ds> temp2_d(2);

  ekat::host_to_device({uobs, vobs, numliqobs, numiceobs, cldliqobs, cldiceobs, tobs, qobs},
                       plev, temp_d);

  std::vector<Int> scalar_sizes = {nelemd, pcnst};
  ekat::host_to_device({dx_short, qmin}, scalar_sizes, temp2_d);

  view_1d
    uobs_d      (temp_d[0]),
    vobs_d      (temp_d[1]),
    numliqobs_d (temp_d[2]),
    numiceobs_d (temp_d[3]),
    cldliqobs_d (temp_d[4]),
    cldiceobs_d (temp_d[5]),
    tobs_d      (temp_d[6]),
    qobs_d      (temp_d[7]);

  sview_1ds
    dx_short_d(reinterpret_cast<Real*>(temp2_d[0].data()), scalar_sizes[0]),
    qmin_d    (reinterpret_cast<Real*>(temp2_d[1].data()), scalar_sizes[1]);

  // Call core function
  DPF::iop_setinitial(plev, pcnst, nelemd, np, nstep, use_replay, dynproc, have_t, have_ps, have_q, have_u, have_v, have_numliq, have_cldliq, have_numice, have_cldice, scm_zero_non_iop_tracers, is_first_restart_step, qmin_d, uobs_d, vobs_d, numliqobs_d, numiceobs_d, cldliqobs_d, cldiceobs_d, psobs, dx_short_d, *dyn_dx_size, *tracers, *elem, tobs_d, qobs_d);

  // Sync back to host
  std::vector<view_1d> inout_views    = {tobs_d, qobs_d};

  ekat::device_to_host({tobs, qobs}, plev, inout_views);
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
void iop_setopts_f(Real scmlat_in, Real scmlon_in, char** iopfile_in, bool single_column_in, bool scm_iop_srf_prop_in, bool iop_nudge_tq_in, bool iop_nudge_uv_in, Real iop_nudge_tq_low_in, Real iop_nudge_tq_high_in, Real iop_nudge_tscale_in, bool scm_observed_aero_in, bool iop_dosubsidence_in, bool scm_multcols_in, bool dp_crm_in, Real iop_perturb_high_in, bool precip_off_in, bool scm_zero_non_iop_tracers_in)
{
#if 0
  using PF = Functions<Real, DefaultDevice>;

  using Spack   = typename PF::Spack;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    Spack iop_nudge_tq_high_in_(iop_nudge_tq_high_in), iop_nudge_tq_low_in_(iop_nudge_tq_low_in), iop_nudge_tscale_in_(iop_nudge_tscale_in), iop_perturb_high_in_(iop_perturb_high_in), scmlat_in_(scmlat_in), scmlon_in_(scmlon_in);
    PF::iop_setopts(scmlat_in_, scmlon_in_, iopfile_in, single_column_in, scm_iop_srf_prop_in, iop_nudge_tq_in, iop_nudge_uv_in, iop_nudge_tq_low_in_, iop_nudge_tq_high_in_, iop_nudge_tscale_in_, scm_observed_aero_in, iop_dosubsidence_in, scm_multcols_in, dp_crm_in, iop_perturb_high_in_, precip_off_in, scm_zero_non_iop_tracers_in);
  });
#endif

}
void setiopupdate_init_f()
{
#if 0
  using PF = Functions<Real, DefaultDevice>;

  using Spack   = typename PF::Spack;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    PF::setiopupdate_init();
  });
#endif

}
void setiopupdate_f()
{
#if 0
  using PF = Functions<Real, DefaultDevice>;

  using Spack   = typename PF::Spack;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    PF::setiopupdate();
  });
#endif

}
void readiopdata_f(Int plev, bool iop_update_phase1, Real* hyam, Real* hybm)
{
  // TODO
}
void iop_intht_f()
{
#if 0
  using PF = Functions<Real, DefaultDevice>;

  using Spack   = typename PF::Spack;

  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
    PF::iop_intht();
  });
#endif

}
} // namespace dp
} // namespace scream
