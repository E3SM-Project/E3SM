#ifndef SCREAM_ZM_FUNCTIONS_F90_HPP
#define SCREAM_ZM_FUNCTIONS_F90_HPP

#include "share/physics/physics_test_data.hpp"
#include "share/core/eamxx_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

// Bridge functions to call fortran version of ZM functions from C++

namespace scream {
namespace zm {

struct IentropyData : public PhysicsTestData {
  // Inputs
  Real s, p, qt, tfg;

  // Outputs
  Real t, qst;

  IentropyData(Real s_, Real p_, Real qt_, Real tfg_, Real t_, Real qst_) :
    PhysicsTestData({}, {}),
    s(s_), p(p_), qt(qt_), tfg(tfg_), t(t_), qst(qst_)
  {}

  PTD_STD_DEF(IentropyData, 6, s, p, qt, tfg, t, qst);
};

struct EntropyData : public PhysicsTestData {
  // Inputs
  Real tk, p, qtot;

  // Outputs
  Real entropy;

  EntropyData(Real tk_, Real p_, Real qtot_, Real entropy_) :
    PhysicsTestData({}, {}),
    tk(tk_), p(p_), qtot(qtot_), entropy(entropy_)
  {}

  PTD_STD_DEF(EntropyData, 4, tk, p, qtot, entropy);
};

struct ZmTransportTracerData : public PhysicsTestData {
  // Inputs
  Int pcols, pver, ncnst, il1g, il2g;
  Int *jt, *mx, *ideep;
  Real dt;
  Real *q, *mu, *md, *du, *eu, *ed, *dp, *fracis, *dpdry;
  bool *doconvtran;

  // Outputs
  Real *dqdt;

  ZmTransportTracerData(Int pcols_, Int pver_, Int ncnst_, Int il1g_, Int il2g_, Real dt_) :
    PhysicsTestData({
      {pcols_, pver_, ncnst_},
      {pcols_, pver_},
      {pcols_},
      {ncnst_}
    },
    {
      {&q, &fracis, &dqdt},
      {&mu, &md, &du, &eu, &ed, &dp, &dpdry}
    },
    {
      {&jt, &mx, &ideep}
    },
    {
      {&doconvtran}
    }),
    pcols(pcols_), pver(pver_), ncnst(ncnst_), il1g(il1g_), il2g(il2g_), dt(dt_)
  {}

  PTD_STD_DEF(ZmTransportTracerData, 6, pcols, pver, ncnst, il1g, il2g, dt);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();
    shift_int_scalar<D>(il1g);
    shift_int_scalar<D>(il2g);
  }

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine);

    // We don't want ideep random, we want it to match current col for now
    for (Int i = 0; i < pcols; ++i) {
      ideep[i] = i;
    }

    // We don't want jt or mx to be random. Those represent cloud tops and bottoms within pver
    for (Int i = 0; i < pcols; ++i) {
      jt[i] = i + 1;
      mx[i] = jt[i] + pver / 2;
      assert(mx[i] < pver);
    }
  }
};

struct ZmTransportMomentumData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, pverp, nwind, il1g, il2g;
  Int *jt, *mx, *ideep;
  Real dt;
  Real *wind_in, *mu, *md, *du, *eu, *ed, *dp;

  // Outputs
  Real *wind_tend, *pguall, *pgdall, *icwu, *icwd, *seten;

  ZmTransportMomentumData(Int pcols_, Int ncol_, Int pver_, Int pverp_, Int nwind_, Int il1g_, Int il2g_, Real dt_) :
    PhysicsTestData({
      {pcols_, pver_, nwind_},
      {pcols_, pver_},
      {pcols_}
    },
    {
      {&wind_in, &wind_tend, &pguall, &pgdall, &icwu, &icwd},
      {&mu, &md, &du, &eu, &ed, &dp, &seten}
    },
    {
      {&jt, &mx, &ideep}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_), pverp(pverp_), nwind(nwind_), il1g(il1g_), il2g(il2g_), dt(dt_)
  {}

  PTD_STD_DEF(ZmTransportMomentumData, 8, pcols, ncol, pver, pverp, nwind, il1g, il2g, dt);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();
    shift_int_scalar<D>(il1g);
    shift_int_scalar<D>(il2g);
  }

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine);

    // We don't want ideep random, we want it to match current col for now
    for (Int i = 0; i < pcols; ++i) {
      ideep[i] = i;
    }

    // We don't want jt or mx to be random. Those represent cloud tops and bottoms within pver
    for (Int i = 0; i < pcols; ++i) {
      jt[i] = i + 1;
      mx[i] = jt[i] + pver / 2;
      assert(mx[i] < pver);
    }
  }
};

struct ComputeDiluteCapeData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, pverp, num_cin, num_msg;
  Int *pblt, *prev_msemax_klev;
  Real *sp_humidity_in, *temperature_in, *zmid, *pmid, *pint, *tpert;
  bool calc_msemax_klev, use_input_tq_mx;

  // Inputs/Outputs
  Int *msemax_klev, *lcl_klev, *eql_klev;
  Real *parcel_qsat, *cape, *q_mx, *t_mx;

  // Outputs
  Real *parcel_temp, *lcl_temperature;

  ComputeDiluteCapeData(Int pcols_, Int ncol_, Int pver_, Int pverp_, Int num_cin_, Int num_msg_, bool calc_msemax_klev_, bool use_input_tq_mx_) :
    PhysicsTestData({
      {pcols_, pver_},
      {pcols_, pverp_},
      {pcols_},
      {pcols_}
    },
    {
      {&sp_humidity_in, &temperature_in, &zmid, &pmid, &parcel_temp, &parcel_qsat},
      {&pint},
      {&tpert, &lcl_temperature, &cape, &q_mx, &t_mx}
    },
    {
      {&pblt, &msemax_klev, &lcl_klev, &eql_klev, &prev_msemax_klev}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_), pverp(pverp_), num_cin(num_cin_), num_msg(num_msg_), calc_msemax_klev(calc_msemax_klev_), use_input_tq_mx(use_input_tq_mx_)
  {}

  PTD_STD_DEF(ComputeDiluteCapeData, 8, pcols, ncol, pver, pverp, num_cin, num_msg, calc_msemax_klev, use_input_tq_mx);

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine, { {pmid, {590, 660}}, {temperature_in, {200, 300}}, {sp_humidity_in, {.004, .018}} });

    // Make sure each column is sorted
    for (Int c = 0; c < pcols; ++c) {
      std::sort(pmid + (c*pver), pmid + ((c+1)*pver));
    }
  }

};

struct FindMseMaxData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, num_msg;
  Int *msemax_top_k;
  Real *temperature, *zmid, *sp_humidity;
  bool pergro_active;

  // Inputs/Outputs
  Int *msemax_klev;
  Real *mse_max_val;

  FindMseMaxData(Int pcols_, Int ncol_, Int pver_, Int num_msg_, bool pergro_active_) :
    PhysicsTestData({
      {pcols_, pver_},
      {pcols_},
      {pcols_}
    },
    {
      {&temperature, &zmid, &sp_humidity},
      {&mse_max_val}
    },
    {
      {&msemax_top_k, &msemax_klev}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_), num_msg(num_msg_), pergro_active(pergro_active_)
  {}

  PTD_STD_DEF(FindMseMaxData, 5, pcols, ncol, pver, num_msg, pergro_active);

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine);

    // We don't want msemax_klev, lcl_klev, or eql_klev to be random
    for (Int i = 0; i < pcols; ++i) {
      msemax_klev[i] = pver / 2 + i;
      msemax_top_k[i]    = i + 1;
    }
  }

};

struct ComputeDiluteParcelData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, num_msg;
  Int *klaunch, *pblt;
  Real *pmid, *temperature, *sp_humidity, *tpert;

  // Inputs/Outputs
  Int *lcl_klev;
  Real *parcel_temp, *parcel_vtemp, *parcel_qsat, *lcl_pmid, *lcl_temperature;

  ComputeDiluteParcelData(Int pcols_, Int ncol_, Int pver_, Int num_msg_) :
    PhysicsTestData({
      {pcols_, pver_},
      {pcols_},
      {pcols_}
    },
    {
      {&pmid, &temperature, &sp_humidity, &parcel_temp, &parcel_vtemp, &parcel_qsat},
      {&tpert, &lcl_pmid, &lcl_temperature}
    },
    {
      {&klaunch, &pblt, &lcl_klev}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_), num_msg(num_msg_)
  {}

  PTD_STD_DEF(ComputeDiluteParcelData, 4, pcols, ncol, pver, num_msg);

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine, { {lcl_pmid, {590, 660}}, {temperature, {200, 300}}, {pmid, {590, 660}}, {sp_humidity, {.004, .018}} });

    // We don't want msemax_klev, lcl_klev, or eql_klev to be random
    for (Int i = 0; i < pcols; ++i) {
      lcl_klev[i]    = pver / 2 - i;
      klaunch[i] = pver / 2 + i;
      pblt[i]    = i + 1;
    }

    // Make sure each column is sorted
    for (Int c = 0; c < pcols; ++c) {
      std::sort(pmid + (c*pver), pmid + ((c+1)*pver));
    }
  }
};

struct ComputeCapeFromParcelData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, pverp, num_cin, num_msg;
  Int *msemax_klev, *lcl_klev;
  Real *temperature, *tv, *sp_humidity, *pint, *lcl_pmid;

  // Inputs/Outputs
  Int *eql_klev;
  Real *parcel_qsat, *parcel_temp, *parcel_vtemp, *cape;

  ComputeCapeFromParcelData(Int pcols_, Int ncol_, Int pver_, Int pverp_, Int num_cin_, Int num_msg_) :
    PhysicsTestData({
      {pcols_, pver_},
      {pcols_, pverp_},
      {pcols_},
      {pcols_}
    },
    {
      {&temperature, &tv, &sp_humidity, &parcel_qsat, &parcel_temp, &parcel_vtemp},
      {&pint},
      {&lcl_pmid, &cape}
    },
    {
      {&msemax_klev, &lcl_klev, &eql_klev}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_), pverp(pverp_), num_cin(num_cin_), num_msg(num_msg_)
  {}

  PTD_STD_DEF(ComputeCapeFromParcelData, 6, pcols, ncol, pver, pverp, num_cin, num_msg);

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine, { {lcl_pmid, {590, 660}}, {temperature, {200, 300}}, {sp_humidity, {.004, .018}} });

    // We don't want msemax_klev, lcl_klev, or eql_klev to be random
    for (Int i = 0; i < pcols; ++i) {
      msemax_klev[i] = pver / 2 + i;
      lcl_klev[i]    = pver / 2 - i;
      eql_klev[i]    = i + 1;
    }
  }
};

struct ZmConvMcspCalculateShearData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver;
  Real *state_pmid, *state_u, *state_v;

  // Outputs
  Real *mcsp_shear;

  ZmConvMcspCalculateShearData(Int pcols_, Int ncol_, Int pver_) :
    PhysicsTestData({
      {pcols_, pver_},
      {pcols_}
    },
    {
      {&state_pmid, &state_u, &state_v},
      {&mcsp_shear}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_)
  {}

  PTD_STD_DEF(ZmConvMcspCalculateShearData, 3, pcols, ncol, pver);

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine, { {state_pmid, {600e2-1000, 600e2+1000}} });

    // Make sure each column is sorted
    for (Int c = 0; c < pcols; ++c) {
      std::sort(state_pmid + (c*pver), state_pmid + ((c+1)*pver));
    }
  }
};

struct ZmConvMcspTendData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, pverp;
  Int *jctop;
  Real ztodt;
  Real *state_pmid, *state_pint, *state_pdel, *state_s, *state_q, *state_u, *state_v, *ptend_zm_s, *ptend_zm_q;

  // Inputs/Outputs
  Real *ptend_s, *ptend_q, *ptend_u, *ptend_v;

  // Outputs
  Real *mcsp_dt_out, *mcsp_dq_out, *mcsp_du_out, *mcsp_dv_out, *mcsp_freq, *mcsp_shear, *zm_depth;

  ZmConvMcspTendData(Int pcols_, Int ncol_, Int pver_, Int pverp_, Real ztodt_) :
    PhysicsTestData({
      {pcols_, pver_},
      {pcols_, pverp_},
      {pcols_},
      {pcols_}
    },
    {
      {&state_pmid, &state_pdel, &state_s, &state_q, &state_u, &state_v, &ptend_zm_s, &ptend_zm_q, &ptend_s, &ptend_q, &ptend_u, &ptend_v, &mcsp_dt_out, &mcsp_dq_out, &mcsp_du_out, &mcsp_dv_out},
      {&state_pint},
      {&mcsp_freq, &mcsp_shear, &zm_depth}
    },
    {
      {&jctop}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_), pverp(pverp_), ztodt(ztodt_)
  {}

  PTD_STD_DEF(ZmConvMcspTendData, 5, pcols, ncol, pver, pverp, ztodt);

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine, { {state_pmid, {600e2-1000, 600e2+1000}}, {state_pint, {1300e2 - 1000, 1300e2 + 1000}}, {state_u, {50, 150}} });

    // Make sure each column is sorted
    for (Int c = 0; c < pcols; ++c) {
      std::sort(state_pmid + (c*pver), state_pmid + ((c+1)*pver));
      std::sort(state_pint + (c*(pverp)), state_pint + ((c+1)*(pverp)));
    }
  }
};

struct ZmConvMainData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, pverp;
  Real time_step;
  Real *t_mid, *q_mid_in, *omega, *p_mid_in, *p_int_in, *p_del_in, *geos, *z_mid_in, *z_int_in, *pbl_hgt, *tpert, *landfrac, *t_star, *q_star;
  bool is_first_step;

  // Outputs
  Int lengath;
  Int *gather_index, *msemax_klev_g, *jctop, *jcbot, *jt;
  Real *prec, *heat, *qtnd, *cape, *dcape, *mcon, *pflx, *zdu, *mflx_up, *entr_up, *detr_up, *mflx_dn, *entr_dn, *p_del, *dsubcld, *ql, *rliq, *rprd, *dlf;

  ZmConvMainData(Int pcols_, Int ncol_, Int pver_, Int pverp_, Real time_step_, bool is_first_step_, Int lengath_) :
    PhysicsTestData({
      {pcols_, pver_},
      {pcols_, pverp_},
      {pcols_},
      {pcols_}
    },
    {
      {&t_mid, &q_mid_in, &omega, &p_mid_in, &p_del_in, &z_mid_in, &t_star, &q_star, &heat, &qtnd, &zdu, &mflx_up, &entr_up, &detr_up, &mflx_dn, &entr_dn, &p_del, &ql, &rprd, &dlf},
      {&p_int_in, &z_int_in, &mcon, &pflx},
      {&geos, &pbl_hgt, &tpert, &landfrac, &prec, &cape, &dcape, &dsubcld, &rliq}
    },
    {
      {&gather_index, &msemax_klev_g, &jctop, &jcbot, &jt}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_), pverp(pverp_), time_step(time_step_), is_first_step(is_first_step_), lengath(lengath_)
  {}

  PTD_STD_DEF(ZmConvMainData, 7, pcols, ncol, pver, pverp, time_step, is_first_step, lengath);
};

struct ZmConvEvapData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, pverp;
  Real time_step;
  Real *p_mid, *p_del, *t_mid, *q_mid, *prdprec, *cldfrc;

  // Inputs/Outputs
  Real *tend_s, *tend_q, *prec;

  // Outputs
  Real *tend_s_snwprd, *tend_s_snwevmlt, *snow, *ntprprd, *ntsnprd, *flxprec, *flxsnow;

  ZmConvEvapData(Int pcols_, Int ncol_, Int pver_, Int pverp_, Real time_step_) :
    PhysicsTestData({
      {pcols_, pver_},
      {pcols_},
      {pcols_, pverp_}
    },
    {
      {&p_mid, &p_del, &t_mid, &q_mid, &prdprec, &cldfrc, &tend_s, &tend_q, &tend_s_snwprd, &tend_s_snwevmlt, &ntprprd, &ntsnprd},
      {&prec, &snow},
      {&flxprec, &flxsnow}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_), pverp(pverp_), time_step(time_step_)
  {}

  PTD_STD_DEF(ZmConvEvapData, 5, pcols, ncol, pver, pverp, time_step);
};

struct ZmCalcFractionalEntrainmentData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, pverp, msg;
  Int *jb, *jt;
  Real *z_mid, *z_int, *dz, *h_env, *h_env_sat;

  // Inputs/Outputs
  Int *j0;
  Real *h_env_min;

  // Outputs
  Real *lambda, *lambda_max;

  ZmCalcFractionalEntrainmentData(Int pcols_, Int ncol_, Int pver_, Int pverp_, Int msg_) :
    PhysicsTestData({
      {pcols_, pver_},
      {pcols_, pverp_},
      {pcols_},
      {pcols_}
    },
    {
      {&z_mid, &dz, &h_env, &h_env_sat, &lambda},
      {&z_int},
      {&h_env_min, &lambda_max}
    },
    {
      {&jb, &jt, &j0}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_), pverp(pverp_), msg(msg_)
  {}

  PTD_STD_DEF(ZmCalcFractionalEntrainmentData, 5, pcols, ncol, pver, pverp, msg);
};

struct ZmDowndraftPropertiesData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, pverp, msg;
  Int *jb, *j0;
  Real *z_int, *dz, *s_mid, *q_mid, *h_env, *lambda, *lambda_max, *qsthat, *hsthat, *gamhat, *rprd, *mflx_up;

  // Inputs/Outputs
  Int *jt, *jd;
  Real *mflx_dn, *entr_dn, *s_dnd, *q_dnd, *h_dnd, *q_dnd_sat, *evp, *totevp;

  ZmDowndraftPropertiesData(Int pcols_, Int ncol_, Int pver_, Int pverp_, Int msg_) :
    PhysicsTestData({
      {pcols_, pverp_},
      {pcols_, pver_},
      {pcols_},
      {pcols_}
    },
    {
      {&z_int},
      {&dz, &s_mid, &q_mid, &h_env, &lambda, &qsthat, &hsthat, &gamhat, &rprd, &mflx_up, &mflx_dn, &entr_dn, &s_dnd, &q_dnd, &h_dnd, &q_dnd_sat, &evp},
      {&lambda_max, &totevp}
    },
    {
      {&jb, &jt, &j0, &jd}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_), pverp(pverp_), msg(msg_)
  {}

  PTD_STD_DEF(ZmDowndraftPropertiesData, 5, pcols, ncol, pver, pverp, msg);
};

struct ZmCloudPropertiesData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, pverp, msg, limcnv;
  Int *jb, *lel;
  Real *p_mid, *z_mid, *z_int, *t_mid, *s_mid, *s_int, *q_mid, *landfrac, *tpert_g;

  // Outputs
  Int *jt, *jlcl, *j0, *jd;
  Real *mflx_up, *entr_up, *detr_up, *mflx_dn, *entr_dn, *mflx_net, *s_upd, *q_upd, *ql, *s_dnd, *q_dnd, *qst, *cu, *evp, *pflx, *rprd;

  ZmCloudPropertiesData(Int pcols_, Int ncol_, Int pver_, Int pverp_, Int msg_, Int limcnv_) :
    PhysicsTestData({
      {pcols_, pver_},
      {pcols_, pverp_},
      {pcols_},
      {pcols_}
    },
    {
      {&p_mid, &z_mid, &t_mid, &s_mid, &s_int, &q_mid, &mflx_up, &entr_up, &detr_up, &mflx_dn, &entr_dn, &mflx_net, &s_upd, &q_upd, &ql, &s_dnd, &q_dnd, &qst, &cu, &evp, &rprd},
      {&z_int, &pflx},
      {&landfrac, &tpert_g}
    },
    {
      {&jb, &lel, &jt, &jlcl, &j0, &jd}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_), pverp(pverp_), msg(msg_), limcnv(limcnv_)
  {}

  PTD_STD_DEF(ZmCloudPropertiesData, 6, pcols, ncol, pver, pverp, msg, limcnv);

  template <ekat::TransposeDirection::Enum D>
  void transition()
  {
    PhysicsTestData::transition<D>();
    shift_int_scalar<D>(limcnv);
  }
};

struct ZmClosureData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, pverp, msg;
  Int *lcl, *lel, *jt, *mx;
  Real cape_threshold_in;
  Real *dsubcld, *z_mid, *z_int, *p_mid, *p_del, *t_mid, *s_mid, *q_mid, *qs, *ql, *s_int, *q_int, *t_pcl_lcl, *t_pcl, *q_pcl_sat, *s_upd, *q_upd, *mflx_net, *detr_up, *mflx_up, *mflx_dn, *q_dnd, *s_dnd, *cape;

  // Outputs
  Real *cld_base_mass_flux;

  ZmClosureData(Int pcols_, Int ncol_, Int pver_, Int pverp_, Int msg_, Real cape_threshold_in_) :
    PhysicsTestData({
      {pcols_},
      {pcols_, pver_},
      {pcols_, pverp_},
      {pcols_}
    },
    {
      {&dsubcld, &t_pcl_lcl, &cape, &cld_base_mass_flux},
      {&z_mid, &p_mid, &p_del, &t_mid, &s_mid, &q_mid, &qs, &ql, &s_int, &q_int, &t_pcl, &q_pcl_sat, &s_upd, &q_upd, &mflx_net, &detr_up, &mflx_up, &mflx_dn, &q_dnd, &s_dnd},
      {&z_int}
    },
    {
      {&lcl, &lel, &jt, &mx}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_), pverp(pverp_), msg(msg_), cape_threshold_in(cape_threshold_in_)
  {}

  PTD_STD_DEF(ZmClosureData, 6, pcols, ncol, pver, pverp, msg, cape_threshold_in);
};

struct ZmCalcOutputTendData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, pverp, msg;
  Int *jt, *mx;
  Real *dsubcld, *p_del, *s_int, *q_int, *s_upd, *q_upd, *mflx_up, *detr_up, *mflx_dn, *s_dnd, *q_dnd, *ql, *evp, *cu;

  // Outputs
  Real *dsdt, *dqdt, *dl;

  ZmCalcOutputTendData(Int pcols_, Int ncol_, Int pver_, Int pverp_, Int msg_) :
    PhysicsTestData({
      {pcols_},
      {pcols_, pver_},
      {pcols_}
    },
    {
      {&dsubcld},
      {&p_del, &s_int, &q_int, &s_upd, &q_upd, &mflx_up, &detr_up, &mflx_dn, &s_dnd, &q_dnd, &ql, &evp, &cu, &dsdt, &dqdt, &dl}
    },
    {
      {&jt, &mx}
    }),
    pcols(pcols_), ncol(ncol_), pver(pver_), pverp(pverp_), msg(msg_)
  {}

  PTD_STD_DEF(ZmCalcOutputTendData, 5, pcols, ncol, pver, pverp, msg);
};

// Glue functions for host test data. We can call either fortran or CXX with this data (_f -> fortran)
void ientropy_f(IentropyData& d);
void ientropy(IentropyData& d);
void entropy_f(EntropyData& d);
void entropy(EntropyData& d);
void zm_transport_tracer_f(ZmTransportTracerData& d);
void zm_transport_tracer(ZmTransportTracerData& d);
void zm_transport_momentum_f(ZmTransportMomentumData& d);
void zm_transport_momentum(ZmTransportMomentumData& d);
void compute_dilute_cape_f(ComputeDiluteCapeData& d);
void compute_dilute_cape(ComputeDiluteCapeData& d);
void find_mse_max_f(FindMseMaxData& d);
void find_mse_max(FindMseMaxData& d);
void compute_dilute_parcel_f(ComputeDiluteParcelData& d);
void compute_dilute_parcel(ComputeDiluteParcelData& d);
void compute_cape_from_parcel_f(ComputeCapeFromParcelData& d);
void compute_cape_from_parcel(ComputeCapeFromParcelData& d);
void zm_conv_mcsp_calculate_shear_f(ZmConvMcspCalculateShearData& d);
void zm_conv_mcsp_calculate_shear(ZmConvMcspCalculateShearData& d);
void zm_conv_mcsp_tend_f(ZmConvMcspTendData& d);
void zm_conv_mcsp_tend(ZmConvMcspTendData& d);
void zm_conv_main_f(ZmConvMainData& d);
void zm_conv_main(ZmConvMainData& d);
void zm_conv_evap_f(ZmConvEvapData& d);
void zm_conv_evap(ZmConvEvapData& d);
void zm_calc_fractional_entrainment_f(ZmCalcFractionalEntrainmentData& d);
void zm_calc_fractional_entrainment(ZmCalcFractionalEntrainmentData& d);
void zm_downdraft_properties_f(ZmDowndraftPropertiesData& d);
void zm_downdraft_properties(ZmDowndraftPropertiesData& d);
void zm_cloud_properties_f(ZmCloudPropertiesData& d);
void zm_cloud_properties(ZmCloudPropertiesData& d);
void zm_closure_f(ZmClosureData& d);
void zm_closure(ZmClosureData& d);
void zm_calc_output_tend_f(ZmCalcOutputTendData& d);
void zm_calc_output_tend(ZmCalcOutputTendData& d);
// End glue function decls

}  // namespace zm
}  // namespace scream

#endif
