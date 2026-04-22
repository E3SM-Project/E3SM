#include "share/util/eamxx_bfbhash.hpp"

#include <catch2/catch.hpp>

#include <cmath>
#include <iostream>
#include <limits>

namespace scream {

using namespace bfbhash;

template <typename T>
static void testeq () {
  int cnt = 0;
  for (const T perturb : {std::numeric_limits<T>::epsilon(), T(-1e-3),
                          T(1e-2)*std::numeric_limits<T>::epsilon()}) {
    const T x = M_PI;
    const T y = M_PI*(1 + perturb);
    HashType a = 0, b = 0;
    hash(x, a);
    hash(y, b);
    REQUIRE((a == b) == (x == y));
    if (a == b) ++cnt;
  }
  REQUIRE(cnt == 1);
}

TEST_CASE("bfbhash")
{
  { // Repeated values should not be nullified.
    HashType a = 0;
    hash(M_PI, a);
    hash(M_PI, a);
    REQUIRE(a != 0);
  }

  { // Negated values should not be nullified.
    HashType a = 0;
    hash( M_PI, a);
    hash(-M_PI, a);
    REQUIRE(a != 0);
  }

  { // The hasher is sensitive to diffs that double accum is not.
    double a = M_PI, b = 1e-20, c = a + b;
    HashType x = 0, y = 0;
    hash(a, x);
    hash(a, y); hash(b, y);
    REQUIRE(a == c); // double add doesn't see b
    REQUIRE(x != y); // but the hasher does
  }

  testeq<float>();
  testeq<double>();

  {
    ekat::Comm comm(MPI_COMM_WORLD);
    HashType a = comm.rank();
    HashType b;
    all_reduce_HashType(MPI_COMM_WORLD, &a, &b, 1);
    HashType c = 0;
    for (int i = 0, n = comm.size(); i < n; ++i) hash(HashType(i), c);
    REQUIRE(b == c);
  }
}

} // namespace scream
