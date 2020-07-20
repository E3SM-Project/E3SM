#include <catch2/catch.hpp>

#include "ekat/util/scream_lin_interp.hpp"
#include "ekat/util/scream_test_utils.hpp"

#include <random>
#include <vector>
#include <algorithm>

extern "C" {

// This will link to the fortran reference implementation
void linear_interp_c(const scream::Real* x1, const scream::Real* x2, const scream::Real* y1, scream::Real* y2, int km1, int km2, int ncol, scream::Real minthresh);

}

namespace {

using scream::Real;
using vector_2d_t = std::vector<std::vector<Real> >;

const Real* flatten(const vector_2d_t& data)
{
  const size_t dim0 = data.size();
  const size_t dim1 = data[0].size();
  Real* result = new Real[dim0 * dim1];

  for (size_t i = 0; i < dim0; ++i) {
    std::copy(data[i].begin(), data[i].end(), result + dim1 * i);
  }

  return result;
}

template <typename Scalar>
void populate_li_input(int km1, int km2, Scalar* x1_i, Scalar* y1_i, Scalar* x2_i, std::default_random_engine* generator = nullptr)
{
  std::default_random_engine local_generator;
  if (generator == nullptr) {
    generator = &local_generator;
  }
  std::uniform_real_distribution<Real> x_dist(0.0,1.0);
  std::uniform_real_distribution<Real> y_dist(0.0,100.0);

  for (int j = 0; j < km1; ++j) {
    x1_i[j] = x_dist(*generator);
    y1_i[j] = y_dist(*generator);
  }
  for (int j = 0; j < km2; ++j) {
    x2_i[j] = x_dist(*generator);
  }

  // make endpoints same
  if (generator != &local_generator) {
    x1_i[0] = 0.0;
    x2_i[0] = 0.0;
    x1_i[km1-1] = 1.0;
    x2_i[km2-1] = 1.0;
  }

  std::sort(x1_i, x1_i + km1);
  std::sort(x2_i, x2_i + km2);
}

TEST_CASE("lin_interp", "soak") {

  std::default_random_engine generator;
  std::uniform_int_distribution<int> k_dist(10,100);
  const Real minthresh = 0.000001;
  const int ncol = 10;

  // increase iterations for a more-thorough soak
  for (int r = 0; r < 100; ++r) {
    const int km1 = k_dist(generator);
    const int km2 = k_dist(generator);
    vector_2d_t
      x1(ncol, std::vector<Real>(km1)),
      x2(ncol, std::vector<Real>(km2)),
      y1(ncol, std::vector<Real>(km1)),
      y2_base(ncol, std::vector<Real>(km2)),
      y2_cmp(ncol, std::vector<Real>(km2)),
      y2_f90(ncol, std::vector<Real>(km2));

    for (int i = 0; i < ncol; ++i) {
      populate_li_input(km1, km2, x1[i].data(), y1[i].data(), x2[i].data(), &generator);
    }

    using LIV = scream::util::LinInterp<Real>;
    using Pack = scream::pack::BigPack<Real>;
    LIV vect(ncol, km1, km2, minthresh);
    const int km1_pack = scream::pack::npack<Pack>(km1);
    const int km2_pack = scream::pack::npack<Pack>(km2);
    typename LIV::template view_2d<Pack>
      x1kv("x1kv", ncol, km1_pack),
      x2kv("x2kv", ncol, km2_pack),
      y1kv("y1kv", ncol, km1_pack),
      y2kv("y2kv", ncol, km2_pack);

    // Initialize kokkos packed inputs
    {
      auto x1kvm = Kokkos::create_mirror_view(x1kv);
      auto x2kvm = Kokkos::create_mirror_view(x2kv);
      auto y1kvm = Kokkos::create_mirror_view(y1kv);

      for (int i = 0; i < ncol; ++i) {
        for (int j = 0; j < km1_pack; ++j) {
          for (int s = 0; s < Pack::n; ++s) {
            if (j*Pack::n + s < km1) {
              x1kvm(i, j)[s] = x1[i][j*Pack::n + s];
              y1kvm(i, j)[s] = y1[i][j*Pack::n + s];
            }
          }
        }
        for (int j = 0; j < km2_pack; ++j) {
          for (int s = 0; s < Pack::n; ++s) {
            if (j*Pack::n + s < km2) {
              x2kvm(i, j)[s] = x2[i][j*Pack::n + s];
            }
          }
        }
      }

      Kokkos::deep_copy(x1kv, x1kvm);
      Kokkos::deep_copy(x2kv, x2kvm);
      Kokkos::deep_copy(y1kv, y1kvm);
    }

    std::vector<Real> x1f(ncol*km1), y1f(ncol*km1), x2f(ncol*km2), y2f(ncol*km2);

    // Initialize fortran inputs
    {
      const Real* x1flat = flatten(x1);
      const Real* x2flat = flatten(x2);
      const Real* y1flat = flatten(y1);

      scream::util::transpose<scream::util::TransposeDirection::c2f>(x1flat, x1f.data(), ncol, km1);
      scream::util::transpose<scream::util::TransposeDirection::c2f>(y1flat, y1f.data(), ncol, km1);
      scream::util::transpose<scream::util::TransposeDirection::c2f>(x2flat, x2f.data(), ncol, km2);

      delete[] x1flat;
      delete[] x2flat;
      delete[] y1flat;
    }

    // Run fortran and store results in y2_f90
    {
      linear_interp_c(x1f.data(), x2f.data(), y1f.data(), y2f.data(), km1, km2, ncol, minthresh);

      std::vector<Real> y2c(ncol*km2);
      scream::util::transpose<scream::util::TransposeDirection::f2c>(y2f.data(), y2c.data(), ncol, km2);

      for (int i = 0; i < ncol; ++i) {
        for (int j = 0; j < km2; ++j) {
          y2_f90[i][j] = y2c[i*km2 + j];
        }
      }
    }

    // Run LiVect
    {
      Kokkos::parallel_for("lin-interp-ut-vect",
                           vect.m_policy,
                           KOKKOS_LAMBDA(typename LIV::MemberType const& team_member) {
        const int i = team_member.league_rank();
        vect.setup(team_member,
                   scream::util::subview(x1kv, i),
                   scream::util::subview(x2kv, i));
        team_member.team_barrier();
        vect.lin_interp(team_member,
                        scream::util::subview(x1kv, i),
                        scream::util::subview(x2kv, i),
                        scream::util::subview(y1kv, i),
                        scream::util::subview(y2kv, i));
      });
    }

    // Compare results
    {
      auto y2kvm = Kokkos::create_mirror_view(y2kv);
      Kokkos::deep_copy(y2kvm, y2kv);
      for (int i = 0; i < ncol; ++i) {
        for (int j = 0; j < km2; ++j) {
          scream::util::catch2_req_pk_sensitive<scream::util::StrictFP,Pack::n>(y2_f90[i][j], y2kvm(i, j / Pack::n)[j % Pack::n]);
        }
      }
    }
  }
}

} // empty namespace
