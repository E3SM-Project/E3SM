// This is a tiny program that calls p3_init() to generate tables used by
// the p3 unit tests.

#include <cstdio>
#include "physics/p3/p3_f90.hpp"

int main(int argc, char* argv[]) {
  std::printf("Generating p3 rain tables...\n");
  scream::p3::p3_init();
  return 0;
}
