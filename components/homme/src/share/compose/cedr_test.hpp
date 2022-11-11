// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_TEST_HPP
#define INCLUDE_CEDR_TEST_HPP

#include "cedr.hpp"
#include "cedr_mpi.hpp"

namespace cedr {
namespace test {
namespace transport1d {

struct Input {
  Int ncells;
  bool verbose;
};

Int run(const mpi::Parallel::Ptr& p, const Input& in);

} // namespace transport1d
} // namespace test
} // namespace cedr

#endif
