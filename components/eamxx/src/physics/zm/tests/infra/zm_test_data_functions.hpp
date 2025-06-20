#ifndef SCREAM_ZM_TEST_DATA_FUNCTIONS_F90_HPP
#define SCREAM_ZM_TEST_DATA_FUNCTIONS_F90_HPP

#include "physics/share/physics_test_data.hpp"
#include "share/eamxx_types.hpp"

#include <array>
#include <utility>
#include <memory>   // for shared_ptr

namespace scream {
namespace zm {

  void zm_test_data_generate_profile( std::mt19937_64 engine, Int ncol, Int pver, Real* zmid, Real* temperature, Real* sp_humidity );

}  // namespace zm
}  // namespace scream

#endif
