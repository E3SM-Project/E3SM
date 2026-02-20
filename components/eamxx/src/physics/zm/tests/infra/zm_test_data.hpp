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
    PhysicsTestData::randomize(engine, { {lcl_pmid, {590, 660}} });

    // We don't want msemax_klev, lcl_klev, or eql_klev to be random
    for (Int i = 0; i < pcols; ++i) {
      lcl_klev[i]    = pver / 2 - i;
      klaunch[i] = pver / 2 + i;
      pblt[i]    = i + 1;
    }
  }

};

struct ComputeCapeFromParcelData : public PhysicsTestData {
  // Inputs
  Int pcols, ncol, pver, pverp, num_cin, num_msg;
  Int *msemax_klev, *lcl_klev;
  Real *temperature, *tv, *zmid, *sp_humidity, *pint, *lcl_pmid;

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
      {&temperature, &tv, &zmid, &sp_humidity, &parcel_qsat, &parcel_temp, &parcel_vtemp},
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
    PhysicsTestData::randomize(engine, { {lcl_pmid, {590, 660}} });

    // We don't want msemax_klev, lcl_klev, or eql_klev to be random
    for (Int i = 0; i < pcols; ++i) {
      msemax_klev[i] = pver / 2 + i;
      lcl_klev[i]    = pver / 2 - i;
      eql_klev[i]    = i + 1;
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
// End glue function decls

}  // namespace zm
}  // namespace scream

#endif
