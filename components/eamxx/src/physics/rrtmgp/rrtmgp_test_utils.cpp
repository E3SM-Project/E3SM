#include "physics/rrtmgp/rrtmgp_test_utils.hpp"

#include <iostream>
#include <netcdf.h>

namespace rrtmgpTest {

bool file_exists(const char *filename) {
  if (auto file = fopen(filename, "r")) {
    fclose(file);
    return true;
  } else {
  return false;
  }
}

}  // namespace rrtmgp
