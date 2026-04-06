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

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine, {
      {t_mid,    {200.0,   300.0}},    // temperature [K]
      {t_star,   {200.0,   300.0}},
      {q_mid_in, {0.001,   0.02}},     // specific humidity [kg/kg]
      {q_star,   {0.001,   0.02}},
      {omega,    {-0.5,    0.5}},      // vertical pressure velocity [Pa/s]
      {p_mid_in, {10000.0, 100000.0}}, // pressure [Pa]
      {p_int_in, {9000.0,  101000.0}}, // interface pressure [Pa]
      {p_del_in, {500.0,   5000.0}},   // pressure thickness [Pa]
      {geos,     {0.0,     1000.0}},   // surface geopotential [m2/s2]
      {z_mid_in, {0.0,     10000.0}},  // mid-point altitude [m]
      {z_int_in, {0.0,     10500.0}},  // interface altitude [m]
      {pbl_hgt,  {100.0,   2000.0}},   // PBL height [m]
      {tpert,    {0.0,     2.0}},      // temperature perturbation [K]
      {landfrac, {0.0,     1.0}},      // land fraction
    });

    // Sort pressure ascending and z descending per column
    for (Int c = 0; c < pcols; ++c) {
      std::sort(p_mid_in + (c*pver), p_mid_in + ((c+1)*pver));
      std::sort(p_int_in + (c*pverp), p_int_in + ((c+1)*pverp));
      std::sort(z_mid_in + (c*pver), z_mid_in + ((c+1)*pver), std::greater<Real>());
      std::sort(z_int_in + (c*pverp), z_int_in + ((c+1)*pverp), std::greater<Real>());
    }
  }
};

struct ZmConvEvapData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, pverp;
  Real time_step;
  bool old_snow;      // runtime_opt: old snow treatment flag (matches Fortran bridge default: true)
  bool pergro_active; // flag for pergro perturbation in snow fraction
  Real ke;            // runtime_opt: evaporation efficiency (matches Fortran bridge default: 2.5E-6)
  Real *p_mid, *p_del, *t_mid, *q_mid, *prdprec, *cldfrc;

  // Inputs/Outputs
  Real *tend_s, *tend_q, *prec;

  // Outputs
  Real *tend_s_snwprd, *tend_s_snwevmlt, *snow, *ntprprd, *ntsnprd, *flxprec, *flxsnow;

  ZmConvEvapData(Int pcols_, Int ncol_, Int pver_, Int pverp_, Real time_step_,
                 bool old_snow_ = true, bool pergro_active_ = false, Real ke_ = 2.5E-6) :
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
    pcols(pcols_), ncol(ncol_), pver(pver_), pverp(pverp_), time_step(time_step_),
    old_snow(old_snow_), pergro_active(pergro_active_), ke(ke_)
  {}

  PTD_STD_DEF(ZmConvEvapData, 8, pcols, ncol, pver, pverp, time_step, old_snow, pergro_active, ke);

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine, {
      {p_mid,   {10000.0, 100000.0}},  // pressure [Pa]
      {p_del,   {2000.0,  10000.0}},   // pressure thickness [Pa]
      {t_mid,   {200.0,   300.0}},     // temperature [K]
      {q_mid,   {0.001,   0.02}},      // specific humidity [kg/kg]
      {prdprec, {0.0,     1e-4}},      // precipitation production [kg/kg/s]
      {cldfrc,  {0.0,     1.0}},       // cloud fraction
      {tend_s,  {-1.0,    1.0}},       // heating rate [J/kg/s]
      {tend_q,  {-1e-5,   1e-5}},      // water vapor tendency [kg/kg/s]
      {prec,    {0.0,     1e-3}},      // precip rate [kg/m2/s] (before *1000 conversion)
      {flxprec, {0.0,     0.1}},       // precip flux at interfaces [kg/m2/s]
      {flxsnow, {0.0,     0.05}},      // snow flux at interfaces [kg/m2/s]
    });
  }
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

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine, {
      {z_mid,     {0.0,    10000.0}},   // altitude [m]
      {z_int,     {0.0,    11000.0}},   // interface altitude [m]
      {dz,        {100.0,  1000.0}},    // layer thickness [m]
      {h_env,     {3.0e5,  3.5e5}},     // env moist static energy [J/kg]
      {h_env_sat, {3.0e5,  3.5e5}},     // env saturated MSE [J/kg]
      {h_env_min, {3.0e5,  3.4e5}},     // min env MSE [J/kg]
      {lambda,    {0.0,    0.0002}},    // fractional entrainment [1/m]
      {lambda_max,{0.0,    0.0002}},    // max fractional entrainment [1/m]
    });

    // jb (cloud base) > j0 (detrainment onset) > jt (cloud top)
    for (Int i = 0; i < pcols; ++i) {
      jt[i] = pver / 4 + i % (pver / 8);
      j0[i] = pver / 3;
      jb[i] = 3 * pver / 4;
    }

    // Sort z descending and dz positive per column
    for (Int c = 0; c < pcols; ++c) {
      std::sort(z_mid + (c*pver), z_mid + ((c+1)*pver), std::greater<Real>());
      std::sort(z_int + (c*pverp), z_int + ((c+1)*pverp), std::greater<Real>());
    }
  }
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

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine, {
      {z_int,     {0.0,    11000.0}},   // interface altitude [m]
      {dz,        {100.0,  1000.0}},    // layer thickness [m]
      {s_mid,     {200.0,  400.0}},     // dry static energy normalized [K]
      {s_dnd,     {200.0,  400.0}},
      {q_mid,     {0.001,  0.02}},      // specific humidity [kg/kg]
      {q_dnd,     {0.001,  0.02}},
      {q_dnd_sat, {0.001,  0.02}},
      {h_env,     {3.0e5,  3.5e5}},     // env moist static energy [J/kg]
      {h_dnd,     {3.0e5,  3.5e5}},     // downdraft moist static energy [J/kg]
      {lambda,    {0.0,    0.0002}},    // fractional entrainment [1/m]
      {lambda_max,{5e-5,   0.0002}},    // non-zero so downdrafts are active
      {qsthat,    {0.001,  0.02}},      // interface saturation humidity
      {hsthat,    {3.0e5,  3.5e5}},     // interface saturated MSE [J/kg]
      {gamhat,    {0.1,    0.5}},       // interface gamma parameter
      {rprd,      {0.0,    1e-4}},      // rain production rate
      {mflx_up,   {0.0,    1.0}},       // updraft mass flux
      {mflx_dn,   {-1.0,   0.0}},       // downdraft mass flux
      {entr_dn,   {0.0,    0.01}},      // downdraft entrainment rate
      {evp,       {0.0,    1e-4}},      // evaporation rate
      {totevp,    {0.0,    0.1}},       // total evaporation
    });

    // jb (cloud base) >= jd (downdraft init) >= j0 (detrainment onset) >= jt (cloud top)
    for (Int i = 0; i < pcols; ++i) {
      jt[i] = pver / 4 + i % (pver / 8);
      j0[i] = pver / 3;
      jd[i] = jt[i] + 1;
      jb[i] = 3 * pver / 4;
    }

    // Sort z_int descending (high altitude at index 0)
    for (Int c = 0; c < pcols; ++c) {
      std::sort(z_int + (c*pverp), z_int + ((c+1)*pverp), std::greater<Real>());
    }
  }
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

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine, {
      {p_mid,   {600.0,  1000.0}},   // pressure [mb]
      {z_mid,   {0.0,    10000.0}},  // altitude [m]
      {z_int,   {0.0,    11000.0}},  // interface altitude [m]
      {t_mid,   {200.0,  300.0}},    // temperature [K]
      {s_mid,   {200.0,  400.0}},    // dry static energy normalized [K]
      {s_int,   {200.0,  400.0}},
      {q_mid,   {0.001,  0.02}},     // specific humidity [kg/kg]
      {landfrac,{0.0,    1.0}},      // land fraction
      {tpert_g, {0.0,    2.0}},      // PBL temperature perturbation [K]
    });

    // jb (cloud base) > lel (equilibrium/launch level, cloud top)
    for (Int i = 0; i < pcols; ++i) {
      lel[i] = pver / 4 + i % (pver / 8);
      jb[i]  = 3 * pver / 4;
    }

    // Sort pressure ascending and z descending per column
    for (Int c = 0; c < pcols; ++c) {
      std::sort(p_mid + (c*pver), p_mid + ((c+1)*pver));
      std::sort(z_mid + (c*pver), z_mid + ((c+1)*pver), std::greater<Real>());
      std::sort(z_int + (c*pverp), z_int + ((c+1)*pverp), std::greater<Real>());
    }
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

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine, {
      {z_mid,      {0.0,    10000.0}},  // altitude [m]
      {z_int,      {0.0,    11000.0}},  // interface altitude [m]
      {p_mid,      {600.0,  1000.0}},   // pressure [mb]
      {p_del,      {10.0,   100.0}},    // pressure thickness [mb]
      {t_mid,      {200.0,  300.0}},    // temperature [K]
      {s_mid,      {200.0,  400.0}},    // dry static energy normalized [K]
      {s_int,      {200.0,  400.0}},
      {s_upd,      {200.0,  400.0}},
      {s_dnd,      {200.0,  400.0}},
      {q_mid,      {0.001,  0.02}},     // specific humidity [kg/kg]
      {qs,         {0.001,  0.02}},
      {q_int,      {0.001,  0.02}},
      {q_pcl_sat,  {0.001,  0.02}},
      {q_upd,      {0.001,  0.02}},
      {q_dnd,      {0.001,  0.02}},
      {ql,         {0.0,    1e-3}},     // cloud liquid water [kg/kg]
      {t_pcl,      {200.0,  300.0}},    // parcel temperature [K]
      {t_pcl_lcl,  {250.0,  300.0}},    // parcel temperature at LCL [K]
      {mflx_net,   {-1.0,   1.0}},      // net mass flux
      {mflx_up,    {0.0,    1.0}},
      {mflx_dn,    {-1.0,   0.0}},
      {detr_up,    {0.0,    0.01}},
      {cape,       {0.0,    2000.0}},   // CAPE [J/kg]
      {dsubcld,    {50.0,   200.0}},    // sub-cloud thickness [mb]
    });

    // level indices must be valid: jt (cloud top) < j0 < lel < lcl < mx (cloud base)
    for (Int i = 0; i < pcols; ++i) {
      jt[i]  = pver / 4 + i % (pver / 8);
      lel[i] = jt[i] + 1;
      lcl[i] = 3 * pver / 4 - 2;
      mx[i]  = 3 * pver / 4;
    }

    // Sort pressure ascending and z descending per column
    for (Int c = 0; c < pcols; ++c) {
      std::sort(p_mid + (c*pver), p_mid + ((c+1)*pver));
      std::sort(z_mid + (c*pver), z_mid + ((c+1)*pver), std::greater<Real>());
      std::sort(z_int + (c*pverp), z_int + ((c+1)*pverp), std::greater<Real>());
    }
  }
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

  template <typename Engine>
  void randomize(Engine& engine)
  {
    PhysicsTestData::randomize(engine, {
      {p_del,   {10.0,  100.0}},   // pressure thickness [mb]
      {s_int,   {200.0, 400.0}},   // dry static energy normalized [K]
      {s_upd,   {200.0, 400.0}},
      {s_dnd,   {200.0, 400.0}},
      {q_int,   {0.001, 0.02}},    // specific humidity [kg/kg]
      {q_upd,   {0.001, 0.02}},
      {q_dnd,   {0.001, 0.02}},
      {mflx_up, {0.0,   1.0}},     // updraft mass flux
      {mflx_dn, {-1.0,  0.0}},     // downdraft mass flux
      {detr_up, {0.0,   0.01}},    // detrainment rate
      {ql,      {0.0,   1e-3}},    // cloud liquid water
      {evp,     {0.0,   1e-4}},    // evaporation rate
      {cu,      {0.0,   1e-4}},    // condensation rate
      {dsubcld, {50.0,  200.0}},   // sub-cloud pressure thickness [mb]
    });

    // jt (cloud top) and mx (cloud base) must be valid level indices with jt < mx
    for (Int i = 0; i < pcols; ++i) {
      jt[i] = pver / 4 + i % (pver / 4);
      mx[i] = 3 * pver / 4;
    }
  }
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
