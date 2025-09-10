#ifndef SCREAM_ZM_TEST_DATA_FUNCTIONS_F90_HPP
#define SCREAM_ZM_TEST_DATA_FUNCTIONS_F90_HPP

#include "physics/share/physics_test_data.hpp"
#include "share/core/eamxx_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

namespace scream {
namespace zm {

template <typename Engine>
void zm_test_data_generate_profile( Engine& engine, Int pver, Int ncol, Real *zmid, Real *temperature, Real *sp_humidity )
{
  std::uniform_real_distribution<Real> lapse_rate_t(3.0, 9.0);
  std::uniform_real_distribution<Real> lapse_rate_q(0.5e-3, 1.5e-3);

  std::uniform_real_distribution<Real> surface_t(285, 305);
  std::uniform_real_distribution<Real> surface_q(0.5e-3, 1.5e-3);

  std::uniform_real_distribution<Real> perturb_t(-0.1,   0.1);
  std::uniform_real_distribution<Real> perturb_q(-0.1e-3, 0.1e-3);

  for (auto i = decltype(ncol){0}; i < ncol; ++i) {
    for (auto k = decltype(pver){0}; k < pver; ++k) {
      int index = pver*i+k;
      zmid       [index] = k*1000;
      temperature[index] = surface_t(engine) - lapse_rate_t(engine)*(pver-1-k) + perturb_t(engine);
      sp_humidity[index] = surface_q(engine) - lapse_rate_q(engine)*(pver-1-k) + perturb_q(engine);
      // apply limiters
      if (temperature[index]<200) { temperature[index] = 200; }
      if (sp_humidity[index]<  0) { sp_humidity[index] =   0; }
    }
  }
}

}  // namespace zm
}  // namespace scream

#endif
