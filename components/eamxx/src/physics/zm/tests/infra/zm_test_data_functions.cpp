#include "zm_test_data.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_assert.hpp"

#include <random>

using scream::Real;
using scream::Int;

//
// A collection of methods to generate data for ZM testing
//

namespace scream {
namespace zm {

void zm_test_data_generate_profile( Int ncol, Int pver, Real *zmid, Real *temperature, Real *sp_humidity ) {
// void zm_test_data_generate_profile( Int ncol, Int pver, Real* zmid, Real* temperature, Real* sp_humidity ) {
// void zm_test_data_generate_profile( Int ncol, Int pver, Real &zmid, Real &temperature, Real &sp_humidity ) {

  auto engine = Base::get_engine();

  std::uniform_real_distribution<Real> lapse_rate_t(3.0, 9.0);
  std::uniform_real_distribution<Real> lapse_rate_q(0.5e-3, 1.5e-3);

  std::uniform_real_distribution<Real> surface_t(285, 305);
  std::uniform_real_distribution<Real> surface_q(0.5e-3, 1.5e-3);

  std::uniform_real_distribution<Real> perturb_t(-0.1,   0.1);
  std::uniform_real_distribution<Real> perturb_q(-0.1e-3, 0.1e-3);

  for (auto i = decltype(ncol){0}; i < ncol; ++i) {
    for (auto k = decltype(pver){0}; k < pver; ++k) {
      int index = nlev*i+k;
      zmid       [index] = k*1000;
      temperature[index] = surface_t(engine) - lapse_rate_t(engine)*(pver-1-k) + perturb_t(engine);
      sp_humidity[index] = surface_q(engine) - lapse_rate_q(engine)*(pver-1-k) + perturb_q(engine);
      // apply limiters
      if temperature[nlev*i + k]<200 { temperature[nlev*i + k] = 200; }
      if sp_humidity[nlev*i + k]<  0 { sp_humidity[nlev*i + k] =   0; }
    }
  }

}

} // namespace zm
} // namespace scream
