// This is a tiny program that calls p3_init() to generate tables used by p3

#include "physics/p3/p3_data.hpp"

int main(int /* argc */, char** /* argv */) {
  scream::p3::p3_init(/* write_tables = */ true);
  return 0;
}
