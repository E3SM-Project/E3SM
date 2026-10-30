#include <catch2/catch.hpp>

#include "share/util/eamxx_tuple.hpp"
#include "share/util/eamxx_bound_functor.hpp"
#include "share/util/eamxx_parallel_for_nd.hpp"

#include <type_traits>

namespace {

using namespace scream;
using details::Tuple;
using details::make_tuple;
using details::get;
using details::apply;

struct IntDoubleInt { int a; double b; int c; };

TEST_CASE ("tuple") {
  SECTION ("size") {
    // No overhead: the tuple is as big as its entries (so, in particular, no padding for an empty tail)
    static_assert (sizeof(Tuple<double>)==sizeof(double));
    static_assert (sizeof(Tuple<double,double>)==2*sizeof(double));
    static_assert (sizeof(Tuple<double,double,double>)==3*sizeof(double));
    static_assert (sizeof(Tuple<int,double,int>)==sizeof(IntDoubleInt));
    static_assert (std::is_trivially_copyable_v<Tuple<int,double,char>>);
  }

  SECTION ("single_entry") {
    auto t = make_tuple(7);
    static_assert (std::is_same_v<decltype(t),Tuple<int>>);
    REQUIRE (get<0>(t)==7);
  }

  SECTION ("get") {
    auto t = make_tuple(1,2.5,'c');
    static_assert (std::is_same_v<decltype(t),Tuple<int,double,char>>);
    REQUIRE (get<0>(t)==1);
    REQUIRE (get<1>(t)==2.5);
    REQUIRE (get<2>(t)=='c');
  }

  SECTION ("apply") {
    // The tuple entries come first, then the extra arguments
    auto f = KOKKOS_LAMBDA (int a, double b, int i, int j) { return a*1000 + b*100 + i*10 + j; };
    REQUIRE (apply(f,make_tuple(1,2.0),3,4)==1234);

    // Extra arguments are forwarded as they are (e.g., non-const references)
    int result = 0;
    auto assign = KOKKOS_LAMBDA (int a, int& r) { r = a; };
    apply(assign,make_tuple(5),result);
    REQUIRE (result==5);
  }
}

TEST_CASE ("bound_functor") {
  using view_t = Kokkos::View<double**>;
  const int n0 = 3, n1 = 4;

  view_t x("x",n0,n1), y("y",n0,n1);
  Kokkos::deep_copy(x,2.0);
  Kokkos::deep_copy(y,1.0);

  SECTION ("parallel_for") {
    auto w = bind_functor(
      KOKKOS_LAMBDA (const view_t& yv, const view_t& xv, const int i, const int j) {
        yv(i,j) += xv(i,j) + i*10 + j;
      }, y, x);

    parallel_for_nd<2,view_t::execution_space>({n0,n1},w);
    Kokkos::fence();

    auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),y);
    for (int i=0; i<n0; ++i)
      for (int j=0; j<n1; ++j)
        REQUIRE (h(i,j)==3.0 + i*10 + j);
  }

  SECTION ("parallel_reduce") {
    // Kokkos passes the reduction result after the indices: it must reach func as a
    // modifiable reference
    auto w = bind_functor(
      KOKKOS_LAMBDA (const view_t& xv, const int i, const int j, double& sum) {
        sum += xv(i,j);
      }, x);

    double sum = 0;
    using MDRange = Kokkos::MDRangePolicy<view_t::execution_space,Kokkos::Rank<2>>;
    Kokkos::parallel_reduce(MDRange({0,0},{n0,n1}),w,Kokkos::Sum<double>(sum));
    REQUIRE (sum==2.0*n0*n1);
  }
}

// NOTE: the kernel is a functor with a templated operator() (not a generic lambda),
//       since nvcc does not allow generic extended lambdas
template<int N>
struct CountVisitsKernel {
  Kokkos::Array<int,7> ext;
  Kokkos::View<int*> visits;
  Kokkos::View<int*> bad;

  template<typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const Idx&... i) const {
    static_assert (sizeof...(i)==N);
    int k = 0;
    bool ok = true;
    ((ok = ok && static_cast<int>(i)>=0 && static_cast<int>(i)<ext[k++]), ...);
    Kokkos::atomic_inc(&visits(0));
    if (not ok)
      Kokkos::atomic_inc(&bad(0));
  }
};

template<int N>
void check_parallel_for_nd (const std::vector<int>& dims)
{
  using exec_space = Kokkos::DefaultExecutionSpace;

  int size = 1;
  Kokkos::Array<int,7> ext;
  for (int k=0; k<N; ++k) {
    size *= dims[k];
    ext[k] = dims[k];
  }

  // Count the visits (rank 0 has 1 entry, and gets no indices), and the out-of-range indices
  Kokkos::View<int*> visits("visits",1);
  Kokkos::View<int*> bad("bad",1);
  CountVisitsKernel<N> kernel{ext,visits,bad};

  parallel_for_nd<N,exec_space>(dims,kernel);
  Kokkos::fence();

  auto h_visits = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),visits);
  auto h_bad    = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),bad);
  REQUIRE (h_visits(0)==size);
  REQUIRE (h_bad(0)==0);
}

TEST_CASE ("parallel_for_nd") {
  check_parallel_for_nd<0>({});
  check_parallel_for_nd<1>({5});
  check_parallel_for_nd<2>({2,3});
  check_parallel_for_nd<3>({2,3,4});
  check_parallel_for_nd<4>({2,3,4,2});
  check_parallel_for_nd<5>({2,3,2,2,3});
  check_parallel_for_nd<6>({2,2,2,2,2,3});

  // Wrong number of dims
  REQUIRE_THROWS (parallel_for_nd<2,Kokkos::DefaultHostExecutionSpace>({2},[](int,int){}));
}

} // anonymous namespace
