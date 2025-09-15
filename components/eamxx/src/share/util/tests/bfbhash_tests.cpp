#include <catch2/catch.hpp>

#include "share/util/eamxx_bfbhash.hpp"

#include <ekat_comm.hpp>
#include <limits>

namespace scream {

template<typename T>
bfbhash::HashType
hash(const scream::bfbhash::HashType h, const T v) {
  auto h2 = h;
  bfbhash::hash(v,h2);
  return h2;
}

TEST_CASE("bfbhash") {
  auto small_d = std::numeric_limits<double>::epsilon();
  auto small_f = std::numeric_limits<float>::epsilon();

  ekat::Comm comm(MPI_COMM_WORLD);

  bfbhash::HashType h1 = 0;

  REQUIRE (h1!=hash(h1,small_d)); // sensitive to small diffs

  auto h2 = hash(h1,small_f);
  if (comm.am_i_root()) {
    h2 = hash(h2,small_f);
  }

  bfbhash::HashType h3 = 0;
  bfbhash::all_reduce_HashType(comm.mpi_comm(),&h2,&h3,1);
  bfbhash::HashType h4 = 0;
  bfbhash::all_reduce_HashType(comm.mpi_comm(),&h1,&h4,1);

  REQUIRE (h3!=h4); // sensitive to small diffs on one rank
}

} // namespace scream
