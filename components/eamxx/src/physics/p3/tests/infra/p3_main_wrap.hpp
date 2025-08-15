#ifndef SCREAM_P3_MAIN_WRAP_HPP
#define SCREAM_P3_MAIN_WRAP_HPP

#include "share/eamxx_types.hpp"
#include <memory>
#include <vector>

namespace scream {
namespace p3 {

struct P3Data;

// Returns number of microseconds of p3_main execution
Int p3_main_wrap(const P3Data& d);

int test_p3_init();

int test_p3_ic();

}  // namespace p3
}  // namespace scream

#endif
