#include <catch2/catch.hpp>

#include "share/util/eamxx_array_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

TEST_CASE ("array_utils") {
  using namespace scream;

  auto engine = setup_random_test ();
  using IPDF = std::uniform_int_distribution<int>;
  IPDF pdf(1,10);

  auto total_size = [](const std::vector<int>& v) -> int {
    int s = 1;
    for (int i : v) {
      s *= i;
    }
    return s;
  };

  // Adds one to fastest striding, doing carrying (if possible) based on max dims d
  // Note: cannot use recursion with a pure lambda
  std::function<bool(int*,int,int*)> add_one = [&](int* v, int n, int* d) -> bool{
    // Increase fastest striding index
    ++v[n];

    // If we reached d[n], we need to carry
    if (v[n]>=d[n]) {
      if (n>0) {
        // Try to carry
        v[n] = 0;

        bool b = add_one(v,n-1,d);

        if (not b) {
          // There was no room to carry. Reset v[n]=d[n] and return false
          v[n] = d[n];
          return false;
        }
      } else {
        v[0] = d[0];
        return false;
      }
    }

    return true;
  };

  for (int rank : {1,2,3,4,5,6}) {
    std::vector<int> dims(rank);
    for (int d=0; d<rank; ++d) {
      dims[d] = pdf(engine);
    }

    std::vector<int> ind(rank,0);
    auto s = total_size(dims);
    for (int idx_1d=0; idx_1d<s; ++idx_1d) {
      auto idx_nd = unflatten_idx(dims,idx_1d);    

      REQUIRE (idx_nd==ind);
      add_one(ind.data(),rank-1,dims.data());
    }
  }
}
