#include <catch2/catch.hpp>

#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include <ekat_pack.hpp>
#include <ekat_subview_utils.hpp>
#include <random>

namespace scream {

TEST_CASE("field_reductions") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using kt = KokkosTypes<DefaultDevice>;

  using P8 = ekat::Pack<Real,8>;

  ekat::Comm comm(MPI_COMM_WORLD);

  std::vector<FieldTag> tags = {COL,LEV};
  std::vector<int> dims = {3,24};

  FieldIdentifier fid ("field_1", {tags,dims}, m/s,"some_grid");
  Field f1(fid);
  f1.get_header().get_alloc_properties().request_allocation(P8::n);
  f1.allocate_view();

  SECTION ("sum") {

    auto v1 = f1.get_strided_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto lsize = fid.get_layout().size();
    auto gsize = lsize*comm.size();
    auto offset = comm.rank()*lsize;
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = offset + idx + 1;
    });
    Kokkos::fence();

    Real lsum = offset*lsize + lsize*(lsize+1) / 2.0;
    Real gsum = gsize*(gsize+1) / 2.0;

    REQUIRE(field_sum(f1).as<Real>()==lsum);
    REQUIRE(field_sum(f1,&comm).as<Real>()==gsum);
  }

  SECTION ("frobenius") {

    auto v1 = f1.get_strided_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto lsize = fid.get_layout().size();
    auto gsize = lsize*comm.size();
    auto offset = comm.rank()*lsize;
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = offset + idx + 1;
    });
    Kokkos::fence();

    // Summing the squares from n=a+1 to n=a+N gives
    // (a+1)^2+(a+2)^2+...
    // which ultimately gives
    // N*a^2 + 2a*(1+2+...+N)
    Real lsum = offset*offset*lsize + 2*offset*lsize*(lsize+1)/2 + lsize*(lsize+1)*(2*lsize+1) / 6.0;
    Real gsum = gsize*(gsize+1)*(2*gsize+1) / 6.0;

    REQUIRE(frobenius_norm(f1).as<Real>()==std::sqrt(lsum));
    REQUIRE(frobenius_norm(f1,&comm).as<Real>()==std::sqrt(gsum));
  }

  SECTION ("inf_norm") {
    using Catch::Matchers::WithinRel;
    std::uniform_real_distribution<Real> rpdf(-1,1);
    auto eps = std::numeric_limits<Real>::epsilon();
    auto tol = 10*eps;

    // Non-negative (0 iff f==0)
    auto seed = get_random_test_seed();
    for (int run=0; run<20; ++run) {
      randomize_uniform(f1,seed++,-1,1);
      REQUIRE (inf_norm(f1,&comm).as<Real>()>=0);
    }
    f1.deep_copy(0);
    REQUIRE (inf_norm(f1,&comm).as<Real>()==0);
    f1.deep_copy(eps);
    REQUIRE (inf_norm(f1,&comm).as<Real>()==eps);

    // Homogeneity
    for (int run=0; run<20; ++run) {
      randomize_uniform(f1,seed++,-1,1);
      auto n1 = inf_norm(f1,&comm).as<Real>();

      std::mt19937_64 engine(seed++);
      auto factor = rpdf(engine);
      comm.broadcast(&factor,1,0);
      f1.scale(factor);
      auto n2 = inf_norm(f1,&comm).as<Real>();
      REQUIRE_THAT (n2,WithinRel(n1*std::abs(factor),tol));

      // Sign invariant
      f1.scale(-1);
      REQUIRE (inf_norm(f1,&comm).as<Real>()==n2);

      // Check norm is the same across ranks
      Real max_across_ranks=-1;
      Real min_across_ranks=-1;
      comm.all_reduce(&n2,&max_across_ranks,1,MPI_MAX);
      comm.all_reduce(&n2,&min_across_ranks,1,MPI_MIN);
      REQUIRE (n2==max_across_ranks);
      REQUIRE (n2==min_across_ranks);
    }

    // Triangular ineq
    auto f2 = f1.clone();
    for (int run=0; run<20; ++run) {
      randomize_uniform(f1,seed++,-1,1);
      randomize_uniform(f2,seed++,-1,1);
      auto n1 = inf_norm(f1,&comm).as<Real>();
      auto n2 = inf_norm(f2,&comm).as<Real>();
      f1.update(f2,1,1);
      auto n3 = inf_norm(f1,&comm).as<Real>();

      REQUIRE(n3<=(n1+n2+tol));
    }
  }

  SECTION ("max") {

    auto v1 = f1.get_strided_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto lsize = fid.get_layout().size();
    auto gsize = lsize*comm.size();
    auto offset = comm.rank()*lsize;
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = offset + idx+1;
    });
    Kokkos::fence();

    Real lmax = offset + lsize;
    Real gmax = gsize;

    REQUIRE(field_max(f1).as<Real>()==lmax);
    REQUIRE(field_max(f1,&comm).as<Real>()==gmax);
  }

  SECTION ("min") {

    auto v1 = f1.get_strided_view<Real**>();
    auto dim0 = fid.get_layout().dim(0);
    auto dim1 = fid.get_layout().dim(1);
    auto lsize = fid.get_layout().size();
    auto offset = comm.rank()*lsize;
    Kokkos::parallel_for(kt::RangePolicy(0,dim0*dim1),
                         KOKKOS_LAMBDA(int idx) {
      int i = idx / dim1;
      int j = idx % dim1;
      v1(i,j) = offset + idx+1;
    });
    Kokkos::fence();

    Real lmin = offset + 1;
    Real gmin = 1;

    REQUIRE(field_min(f1).as<Real>()==lmin);
    REQUIRE(field_min(f1,&comm).as<Real>()==gmin);
  }
}

} // namespace scream
