#ifndef SCREAM_PARALLEL_FOR_ND_HPP
#define SCREAM_PARALLEL_FOR_ND_HPP

#include <Kokkos_Core.hpp>
#include <ekat_assert.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace scream
{

namespace impl {

// Kokkos has no 0-dim policy, so rank-0 kernels run on a 1-entry RangePolicy, which calls
// the functor with a dummy index. This drops that index, so a rank-0 kernel gets no indices.
template<typename Functor>
struct IgnoreIndex {
  Functor functor;

  KOKKOS_INLINE_FUNCTION
  void operator() (const int) const { functor(); }
};

} // namespace impl

// Run functor(i_0,...,i_{N-1}) for all indices in [0,dims[0])x...x[0,dims[N-1]).
// The rank N must be in [0,6]. For N=0, functor is called once, with no indices.
// NOTE: this does NOT fence.
template<int N, typename ExecSpace, typename Functor>
void parallel_for_nd (const std::vector<int>& dims, const Functor& functor)
{
  static_assert (N>=0 && N<=6, "[parallel_for_nd] Error! Unsupported rank (must be in [0,6]).\n");
  EKAT_REQUIRE_MSG (static_cast<int>(dims.size())==N,
      "[parallel_for_nd] Error! Number of dimensions does not match the rank.\n"
      " - rank: " + std::to_string(N) + "\n"
      " - dims.size(): " + std::to_string(dims.size()) + "\n");

  if constexpr (N==0) {
    Kokkos::RangePolicy<ExecSpace> policy(0,1);
    Kokkos::parallel_for(policy,impl::IgnoreIndex<Functor>{functor});
  } else if constexpr (N==1) {
    Kokkos::RangePolicy<ExecSpace> policy(0,dims[0]);
    Kokkos::parallel_for(policy,functor);
  } else {
    using MDRange = Kokkos::MDRangePolicy<
                      ExecSpace,
                      Kokkos::Rank<N,Kokkos::Iterate::Right,Kokkos::Iterate::Right>
                    >;
    Kokkos::Array<std::int64_t,N> lower, upper;
    for (int i=0; i<N; ++i) {
      lower[i] = 0;
      upper[i] = dims[i];
    }
    Kokkos::parallel_for(MDRange(lower,upper),functor);
  }
}

} // namespace scream

#endif // SCREAM_PARALLEL_FOR_ND_HPP
