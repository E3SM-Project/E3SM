#ifndef SCREAM_P3_MAIN_WRAP_HPP
#define SCREAM_P3_MAIN_WRAP_HPP

#include "share/scream_types.hpp"
#include <memory>
#include <vector>

namespace scream {
namespace p3 {

struct FortranData;

// Returns number of microseconds of p3_main execution
Int p3_main_wrap(const FortranData& d);

int test_p3_init();

int test_p3_ic();

}  // namespace p3
}  // namespace scream

#endif
