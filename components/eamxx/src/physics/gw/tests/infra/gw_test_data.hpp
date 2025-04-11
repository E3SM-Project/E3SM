#ifndef SCREAM_GW_FUNCTIONS_F90_HPP
#define SCREAM_GW_FUNCTIONS_F90_HPP

#include "physics/gw/gw_functions.hpp"
#include "physics/share/physics_test_data.hpp"
#include "share/scream_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

//
// Bridge functions to call fortran version of gw functions from C++
//

namespace scream {
namespace gw {

struct GwdComputeTendenciesFromStressDivergenceData : public PhysicsTestData {
  // Inputs
  Int ncol, pver, pgwv, ngwv;
  bool do_taper;
  Real dt, effgw;
  Int *tend_level;
  Real *lat, *dpm, *rdpm, *c, *ubm, *t, *nm, *xv, *yv;

  // Inputs/Outputs
  Real *tau;

  // Outputs
  Real *gwut, *utgw, *vtgw;

  GwdComputeTendenciesFromStressDivergenceData(Int ncol_, Int pver_, Int pgwv_, Int ngwv_, bool do_taper_, Real dt_, Real effgw_) :
    PhysicsTestData({
      {ncol_},
      {ncol_, pver_},
      {ncol_, 2*pgwv_},
      {ncol_, 2*pgwv_, pver_ + 1},
      {ncol_, pver_, 2*ngwv_},
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
    ncol(ncol_), pver(pver_), pgwv(pgwv_), ngwv(ngwv_), do_taper(do_taper_), dt(dt_), effgw(effgw_)
  {}

  PTD_STD_DEF(GwdComputeTendenciesFromStressDivergenceData, 7, ncol, pver, pgwv, ngwv, do_taper, dt, effgw);
};
// Glue functions to call fortran from from C++ with the Data struct

void gwd_compute_tendencies_from_stress_divergence(GwdComputeTendenciesFromStressDivergenceData& d);
extern "C" { // _f function decls
}

}  // namespace gw
}  // namespace scream

#endif
