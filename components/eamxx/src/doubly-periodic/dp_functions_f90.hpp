#ifndef SCREAM_DP_FUNCTIONS_F90_HPP
#define SCREAM_DP_FUNCTIONS_F90_HPP

#include "share/scream_types.hpp"
#include "physics/share/physics_test_data.hpp"

#include "dp_functions.hpp"
#include "physics_constants.hpp"

#include <vector>
#include <array>
#include <utility>

//
// Bridge functions to call fortran version of dp functions from C++
//

namespace scream {
namespace dp {

struct AdvanceIopForcingData : public PhysicsTestData {
  static constexpr size_t NUM_ARRAYS = 20;

  // Inputs
  Int plev, pcnst;
  Real scm_dt, ps_in;
  bool have_u, have_v, dp_crm, use_3dfrc;
  Real *u_in, *v_in, *t_in, *q_in, *t_phys_frc, *divt3d, *divq3d, *divt,
    *divq, *wfld, *uobs, *vobs, *hyai, *hyam, *hybi, *hybm;

  // Outputs
  Real *u_update, *v_update, *t_update, *q_update;

  AdvanceIopForcingData(Int plev_, Int pcnst_, Real scm_dt_, Real ps_in_, bool have_u_, bool have_v_, bool dp_crm_, bool use_3dfrc_) :
    PhysicsTestData(
      {{ plev_ }, { pcnst_, plev_ }},
      {
        { &u_in, &v_in, &t_in, &t_phys_frc, &divt3d, &divt, &divq, &wfld, &uobs, &vobs, &hyai, &hyam, &hybi, &hybm, &u_update, &v_update, &t_update },
        { &q_in, &divq3d, &divq, &q_update }
      }),
    plev(plev_), pcnst(pcnst_), scm_dt(scm_dt_), ps_in(ps_in_), have_u(have_u_), have_v(have_v_), dp_crm(dp_crm_), use_3dfrc(use_3dfrc_) {}

  PTD_STD_DEF(AdvanceIopForcingData, 8, plev, pcnst, scm_dt, ps_in, have_u, have_v, dp_crm, use_3dfrc);
};


struct AdvanceIopNudgingData : public PhysicsTestData {
  static constexpr size_t NUM_ARRAYS = 12;

  // Inputs
  Int plev;
  Real scm_dt, ps_in;
  Real *t_in, *q_in, *tobs, *qobs, *hyai, *hyam, *hybi, *hybm;

  // Outputs
  Real *t_update, *q_update, *relaxt, *relaxq;

  AdvanceIopNudgingData(Int plev_, Real scm_dt_, Real ps_in_) :
    PhysicsTestData({{ plev_ }}, {{ &t_in, &q_in, &tobs, &qobs, &hyai, &hyam, &hybi, &hybm, &t_update, &q_update, &relaxt, &relaxq }}),
      plev(plev_), scm_dt(scm_dt_), ps_in(ps_in_) {}

  PTD_STD_DEF(AdvanceIopNudgingData, 3, plev, scm_dt, ps_in);
};

struct AdvanceIopSubsidenceData : public PhysicsTestData {
  static constexpr size_t NUM_ARRAYS = 13;

  // Inputs
  Int plev, pcnst;
  Real scm_dt, ps_in;
  Real *u_in, *v_in, *t_in, *q_in, *hyai, *hyam, *hybi, *hybm, *wfld;

  // Outputs
  Real *u_update, *v_update, *t_update, *q_update;

  AdvanceIopSubsidenceData(Int plev_, Int pcnst_, Real scm_dt_, Real ps_in_) :
    PhysicsTestData(
      {{ plev_ }, { plev_, pcnst_ }},
      {
        { &u_in, &v_in, &t_in, &hyai, &hyam, &hybi, &hybm, &wfld, &u_update, &v_update, &t_update },
        { &q_in, &q_update }
      }),
    plev(plev_), pcnst(pcnst_), scm_dt(scm_dt_), ps_in(ps_in_) {}

  PTD_STD_DEF(AdvanceIopSubsidenceData, 4, plev, pcnst, scm_dt, ps_in);
};

struct IopSetinitialData : public PhysicsTestData {
  // Inputs
  Int plev, pcnst, nelemd, np, nstep;

  bool use_replay, dynproc, have_t, have_q, have_ps, have_u, have_v, have_numliq, have_cldliq, have_numice, have_cldice, scm_zero_non_iop_tracers, is_first_restart_step;

  Real psobs;

  Real* qmin, *uobs, *vobs, *numliqobs, *numiceobs, *cldliqobs, *cldiceobs, *dx_short;

  // Inputs/Outputs
  tracer_t tracers;
  element_t elem;

  Real dyn_dx_size;

  Real* tobs, *qobs;

  IopSetinitialData(
    Int plev_, Int pcnst_, Int nelemd_, Int np_, Int nstep_, Real psobs_,
    bool use_replay_, bool dynproc_, bool have_t_, bool have_q_, bool have_ps_, bool have_u_, bool have_v_, bool have_numliq_, bool have_cldliq_, bool have_numice_, bool have_cldice_, bool scm_zero_non_iop_tracers_, bool is_first_restart_step_) :
    PhysicsTestData(
      {{nelemd_}, {pcnst_}, {plev_}},
      {
        {&dx_short},
        {&qmin},
        {&uobs, &vobs, &numliqobs, &numiceobs, &cldliqobs, &cldiceobs, &tobs, &qobs}
      }),
    plev(plev_), pcnst(pcnst_), nelemd(nelemd_), np(np_), nstep(nstep_), psobs(psobs_),
    use_replay(use_replay_), dynproc(dynproc_), have_t(have_t_), have_q(have_q_), have_ps(have_ps_), have_u(have_u_), have_v(have_v_), have_numliq(have_numliq_), have_cldliq(have_cldliq_), have_numice(have_numice_), have_cldice(have_cldice_), scm_zero_non_iop_tracers(scm_zero_non_iop_tracers_), is_first_restart_step(is_first_restart_step_) { }

  PTD_STD_DEF(IopSetinitialData, 19, plev, pcnst, nelemd, np, nstep, psobs, use_replay, dynproc, have_t, have_q, have_ps, have_u, have_v, have_numliq, have_cldliq, have_numice, have_cldice, scm_zero_non_iop_tracers, is_first_restart_step);

  void init()
  {
    tracers.init(nelemd, QSIZE_D);
    elem.init(nelemd, true, true, 2);
  }

  void randomize(std::mt19937_64& engine)
  {
    PhysicsTestData::randomize(engine);
    init();

    tracers.randomize(engine(), 0, 1);
    elem.randomize(engine());

    // Where tobs starts being non-zero is important to the algorithm
    std::uniform_int_distribution<Int> default_int_dist(0, plev);
    Int start_nz = default_int_dist(engine);

    for (Int k = 0; k < start_nz; ++k) {
      tobs[k] = 0;
    }
  }
};

struct IopBroadcastData : public PhysicsTestData {
  // Inputs
  Int plev;

  IopBroadcastData(Int plev_=0) :
    PhysicsTestData({}, {}), plev(plev_) {}

  PTD_STD_DEF(IopBroadcastData, 1, plev);
};

struct ApplyIopForcingData : public PhysicsTestData {
  // Inputs
  Int plev, nelemd, n, nets, nete;
  hybrid_t hybrid;
  timelevel_t tl;
  bool t_before_advance;

  // Inputs/Outputs
  element_t *elem;
  hvcoord_t hvcoord;

  ApplyIopForcingData(Int plev_, Int nelemd_, Int n_, Int nets_, Int nete_, bool t_before_advance_) :
    PhysicsTestData({}, {}), plev(plev_), nelemd(nelemd_), n(n_), nets(nets_), nete(nete_), t_before_advance(t_before_advance_) {}

  PTD_STD_DEF(ApplyIopForcingData, 6, plev, nelemd, n, nets, nete, t_before_advance);
};

struct IopDomainRelaxationData : public PhysicsTestData {
  // Inputs
  Int nelemd, np, nlev, t1, nelemd_todo, np_todo;
  hvcoord_t hvcoord;
  hybrid_t hybrid;
  Real dt;

  // Inputs/Outputs
  element_t *elem;
  Real *dp;

  IopDomainRelaxationData(Int nelemd_, Int np_, Int nlev_, Int t1_, Int nelemd_todo_, Int np_todo_, Real dt_) :
    PhysicsTestData({{ np_, np_, nlev_ }}, {{ &dp }}), nelemd(nelemd_), np(np_), nlev(nlev_), t1(t1_), nelemd_todo(nelemd_todo_), np_todo(np_todo_), dt(dt_) {}

  PTD_STD_DEF(IopDomainRelaxationData, 7, nelemd, np, nlev, t1, nelemd_todo, np_todo, dt);
};

struct CrmResolvedTurbData : public PhysicsTestData {
  // Inputs
  Int plev, nelemd, t1, nelemd_todo, np_todo;
  hvcoord_t hvcoord;
  hybrid_t hybrid;

  // Inputs/Outputs
  element_t *elem;

  CrmResolvedTurbData(Int plev_, Int nelemd_, Int t1_, Int nelemd_todo_, Int np_todo_) :
    PhysicsTestData({}, {}), plev(plev_), nelemd(nelemd_), t1(t1_), nelemd_todo(nelemd_todo_), np_todo(np_todo_) {}

  PTD_STD_DEF(CrmResolvedTurbData, 5, plev, nelemd, t1, nelemd_todo, np_todo);
};

struct IopDefaultOptsData {
  // Inputs
  Int plev;

  // Outputs
  Real scmlat_out, scmlon_out, iop_nudge_tq_low_out, iop_nudge_tq_high_out, iop_nudge_tscale_out, iop_perturb_high_out;
  std::string iopfile_out;
  bool single_column_out, scm_iop_srf_prop_out, iop_nudge_tq_out, iop_nudge_uv_out, scm_observed_aero_out, iop_dosubsidence_out, scm_multcols_out, dp_crm_out, precip_off_out, scm_zero_non_iop_tracers_out;

  void randomize(std::mt19937_64& engine) {}

  IopDefaultOptsData() = default;
};

struct IopSetoptsData {
  // Inputs
  Int plev;
  Real scmlat_in, scmlon_in, iop_nudge_tq_low_in, iop_nudge_tq_high_in, iop_nudge_tscale_in, iop_perturb_high_in;
  std::string iopfile_in;
  bool single_column_in, scm_iop_srf_prop_in, iop_nudge_tq_in, iop_nudge_uv_in, scm_observed_aero_in, iop_dosubsidence_in, scm_multcols_in, dp_crm_in, precip_off_in, scm_zero_non_iop_tracers_in;

  void randomize(std::mt19937_64& engine) {}

  IopSetoptsData() = default;
};

struct SetiopupdateInitData {
  // Inputs
  Int plev;

  void randomize(std::mt19937_64& engine) {}
};

struct SetiopupdateData {
  // Inputs
  Int plev;

  void randomize(std::mt19937_64& engine) {}
};

struct ReadiopdataData : public PhysicsTestData {
  // Inputs
  Int plev;
  bool iop_update_phase1;
  Real *hyam, *hybm;

  ReadiopdataData(Int plev_, bool iop_update_phase1_) :
    PhysicsTestData({{ plev_ }}, {{ &hyam, &hybm }}), plev(plev_), iop_update_phase1(iop_update_phase1_) {}

  PTD_STD_DEF(ReadiopdataData, 2, plev, iop_update_phase1);
};

struct IopInthtData {
    // Inputs
  Int plev;

  void randomize(std::mt19937_64& engine) {}
};

// Glue functions to call fortran from from C++ with the Data struct

void advance_iop_forcing(AdvanceIopForcingData& d);
void advance_iop_nudging(AdvanceIopNudgingData& d);
void advance_iop_subsidence(AdvanceIopSubsidenceData& d);
void iop_setinitial(IopSetinitialData& d);
void iop_broadcast(IopBroadcastData& d);
void apply_iop_forcing(ApplyIopForcingData& d);
void iop_domain_relaxation(IopDomainRelaxationData& d);
void crm_resolved_turb(CrmResolvedTurbData& d);
void iop_default_opts(IopDefaultOptsData& d);
void iop_setopts(IopSetoptsData& d);
void setiopupdate_init(SetiopupdateInitData& d);
void setiopupdate(SetiopupdateData& d);
void readiopdata(ReadiopdataData& d);
void iop_intht(IopInthtData& d);
extern "C" { // _f function decls

void advance_iop_forcing_f(Int plev, Int pcnst, Real scm_dt, Real ps_in, bool have_u, bool have_v, bool dp_crm, bool use_3dfrc, Real* u_in, Real* v_in, Real* t_in, Real* q_in, Real* t_phys_frc, Real* divt3d, Real* divq3d, Real* divt, Real* divq, Real* wfld, Real* uobs, Real* vobs, Real* hyai, Real* hyam, Real* hybi, Real* hybm, Real* u_update, Real* v_update, Real* t_update, Real* q_update);

void advance_iop_nudging_f(Int plev, Real scm_dt, Real ps_in, Real* t_in, Real* q_in, Real* tobs, Real* qobs,
                           Real* hyai, Real* hyam, Real* hybi, Real* hybm,
                           Real* t_update, Real* q_update, Real* relaxt, Real* relaxq);

void advance_iop_subsidence_f(Int plev, Int pcnst, Real scm_dt, Real ps_in, Real* u_in, Real* v_in, Real* t_in, Real* q_in, Real* hyai, Real* hyam, Real* hybi, Real* hybm, Real* wfld, Real* u_update, Real* v_update, Real* t_update, Real* q_update);

void iop_setinitial_f(Int plev, Int pcnst, Int nelemd, Int np, Int nstep, Real psobs, bool use_replay, bool dynproc, bool have_t, bool have_q, bool have_ps, bool have_u, bool have_v, bool have_numliq, bool have_cldliq, bool have_numice, bool have_cldice, bool scm_zero_non_iop_tracers, bool is_first_restart_step, Real* qmin, Real* uobs, Real* vobs, Real* numliqobs, Real* numiceobs, Real* cldliqobs, Real* cldiceobs, Real* dx_short, tracer_t* tracers, element_t* elem, Real* dyn_dx_size, Real* tobs, Real* qobs);

void iop_broadcast_f();
void apply_iop_forcing_f(Int nelemd, element_t* elem, hvcoord_t* hvcoord, hybrid_t hybrid, timelevel_t tl, Int n, bool t_before_advance, Int nets, Int nete);
void iop_domain_relaxation_f(Int nelemd, Int np, Int nlev, element_t* elem, hvcoord_t hvcoord, hybrid_t hybrid, Int t1, Real* dp, Int nelemd_todo, Int np_todo, Real dt);
void crm_resolved_turb_f(Int nelemd, element_t* elem, hvcoord_t hvcoord, hybrid_t hybrid, Int t1, Int nelemd_todo, Int np_todo);
void iop_default_opts_f(Real* scmlat_out, Real* scmlon_out, char** iopfile_out, bool* single_column_out, bool* scm_iop_srf_prop_out, bool* iop_nudge_tq_out, bool* iop_nudge_uv_out, Real* iop_nudge_tq_low_out, Real* iop_nudge_tq_high_out, Real* iop_nudge_tscale_out, bool* scm_observed_aero_out, bool* iop_dosubsidence_out, bool* scm_multcols_out, bool* dp_crm_out, Real* iop_perturb_high_out, bool* precip_off_out, bool* scm_zero_non_iop_tracers_out);
void iop_setopts_f(Real scmlat_in, Real scmlon_in, const char** iopfile_in, bool single_column_in, bool scm_iop_srf_prop_in, bool iop_nudge_tq_in, bool iop_nudge_uv_in, Real iop_nudge_tq_low_in, Real iop_nudge_tq_high_in, Real iop_nudge_tscale_in, bool scm_observed_aero_in, bool iop_dosubsidence_in, bool scm_multcols_in, bool dp_crm_in, Real iop_perturb_high_in, bool precip_off_in, bool scm_zero_non_iop_tracers_in);
void setiopupdate_init_f();
void setiopupdate_f();
void readiopdata_f(Int plev, bool iop_update_phase1, Real* hyam, Real* hybm);
void iop_intht_f();
} // end _f function decls

}  // namespace dp
}  // namespace scream

#endif // SCREAM_DP_FUNCTIONS_F90_HPP
