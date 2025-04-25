#ifndef SCREAM_GW_FUNCTIONS_F90_HPP
#define SCREAM_GW_FUNCTIONS_F90_HPP

#include "physics/gw/gw_functions.hpp"
#include "physics/share/physics_test_data.hpp"
#include "share/eamxx_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

//
// Bridge functions to call fortran version of gw functions from C++
//

namespace scream {
namespace gw {

// The Data struct is special; it is used to do gw initialization, which
// must be called before any gw function.
struct GwInit : public PhysicsTestData {
  // Inputs
  Int pver, pgwv;
  Real dc;
  bool orographic_only, do_molec_diff, tau_0_ubc;
  Int nbot_molec, ktop, kbotbg;
  Real fcrit2, kwv;
  Real *cref, *alpha;

  GwInit(Int pver_, Int pgwv_, Real dc_, bool orographic_only_, bool do_molec_diff_, bool tau_0_ubc_, Int nbot_molec_, Int ktop_, Int kbotbg_, Real fcrit2_, Real kwv_) :
    PhysicsTestData({
      {pgwv_ * 2},
      {pver_ + 1}
    },
    {
      {&cref},
      {&alpha}
    }),
    pver(pver_), pgwv(pgwv_), dc(dc_), orographic_only(orographic_only_), do_molec_diff(do_molec_diff_), tau_0_ubc(tau_0_ubc_), nbot_molec(nbot_molec_), ktop(ktop_), kbotbg(kbotbg_), fcrit2(fcrit2_), kwv(kwv_)
  {
    // Assert valid init data?
    assert(ktop <= pver);
    assert(kbotbg >= 0);
    assert(kbotbg <= ktop);
    assert(pgwv > 0);
    assert(nbot_molec >= 0);
    assert(nbot_molec <= ktop);
  }

  PTD_STD_DEF(GwInit, 11, pver, pgwv, dc, orographic_only, do_molec_diff, tau_0_ubc, nbot_molec, ktop, kbotbg, fcrit2, kwv);
};


struct GwdComputeTendenciesFromStressDivergenceData : public PhysicsTestData {
  // Inputs
  Int ncol, ngwv;
  bool do_taper;
  Real dt, effgw;
  Int *tend_level;
  Real *lat, *dpm, *rdpm, *c, *ubm, *t, *nm, *xv, *yv;
  GwInit init;

  // Inputs/Outputs
  Real *tau;

  // Outputs
  Real *gwut, *utgw, *vtgw;

  GwdComputeTendenciesFromStressDivergenceData(Int ncol_, Int ngwv_, bool do_taper_, Real dt_, Real effgw_, GwInit init_) :
    PhysicsTestData({
      {ncol_},
      {ncol_, init_.pver},
      {ncol_, 2*init_.pgwv},
      {ncol_, 2*init_.pgwv, init_.pver + 1},
      {ncol_, init_.pver, 2*ngwv_},
      {ncol_}
    },
    {
      {&lat, &xv, &yv},
      {&dpm, &rdpm, &ubm, &t, &nm, &utgw, &vtgw},
      {&c},
      {&tau},
      {&gwut}
    },
    {
      {&tend_level}
    }),
    ncol(ncol_), ngwv(ngwv_), do_taper(do_taper_), dt(dt_), effgw(effgw_), init(init_)
  {}

  PTD_STD_DEF_INIT(GwdComputeTendenciesFromStressDivergenceData, 5, ncol, ngwv, do_taper, dt, effgw);
};

// Glue functions to call fortran from from C++ with the Data struct

void gwd_compute_tendencies_from_stress_divergence(GwdComputeTendenciesFromStressDivergenceData& d);
extern "C" { // _f function decls
}

}  // namespace gw
}  // namespace scream

#endif
