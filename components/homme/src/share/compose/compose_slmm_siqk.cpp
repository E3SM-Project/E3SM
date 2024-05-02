#include "compose_slmm_siqk.hpp"

namespace siqk {
namespace sqr {
namespace test {

struct Info {
  Int sum_nits, max_nits, nfails;
};

class TestSphereToRefKernel {
  const Real a_test[9] = {-0.1, -1e-16, 0, 1e-15, 0.1, 0.7, 1, 1-1e-14, 1.1};
  const Int n_a_test = sizeof(a_test)/sizeof(*a_test);

  const Real tol_;
  const bool sphere_;
  mutable ConstVec3s p_;
  mutable ConstIdxs e_;

public:
  typedef Info value_type;

  TestSphereToRefKernel (const ConstVec3s::HostMirror& p_hm,
                         const ConstIdxs::HostMirror& e_hm,
                         const bool sphere,
                         const Real tol = 1e1*std::numeric_limits<Real>::epsilon())
    : tol_(tol), sphere_(sphere)
  {
    { Vec3s p; resize_and_copy(p, p_hm); p_ = p; }
    { Idxs e; resize_and_copy(e, e_hm); e_ = e; }
  }

  Int n () const { return nslices(e_)*square(n_a_test); }
  const Real& tol () const { return tol_; }

  KOKKOS_INLINE_FUNCTION
  void operator() (const Int k, value_type& jinfo) const {
    const Int
      ei = k / square(n_a_test),
      ij = k % square(n_a_test),
      i = ij / n_a_test,
      j = ij % n_a_test;
    const auto cell = slice(e_, ei);
    const Real a_t = a_test[i], b_t = a_test[j];
    Real q[3];
    if (sphere_) sqr::calc_ref_to_sphere(p_, cell, a_t, b_t, q);
    else         sqr::calc_ref_to_plane (p_, cell, a_t, b_t, q);
    Real length_scale = 0;
    for (Int i = 0; i < 4; ++i)
      length_scale = siqk::max(length_scale,
                               SphereGeometry::norm2(slice(p_, cell[i])));
    length_scale = std::sqrt(length_scale);
    Real a, b;
    sqr::Info info;
    if (sphere_)
      sqr::calc_sphere_to_ref(p_, cell, q, a, b, &info, 100, tol_*length_scale);
    else
      sqr::calc_plane_to_ref (p_, cell, q, a, b, &info, 100, tol_*length_scale);
    const Real err = std::sqrt(square(a_t - a) + square(b_t - b));
    // tol is on dx, not (a,b), so adjust slightly.
    if ( ! info.success || err > 1e4*tol_) {
      jinfo.nfails++;
      printf("calc_sphere_to_ref ei %d i %d j %d: nits %d re %1.1e\n",
             ei, i, j, info.n_iterations, err);
    }
    jinfo.sum_nits += info.n_iterations;
    jinfo.max_nits = max(jinfo.max_nits, info.n_iterations);
  }

  KOKKOS_INLINE_FUNCTION
  void init (value_type& info) {
    info.sum_nits = 0;
    info.max_nits = 0;
    info.nfails = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void join (value_type& dst, value_type const& src) const {
    dst.max_nits = max(dst.max_nits, src.max_nits);
    dst.sum_nits += src.sum_nits;
    dst.nfails += src.nfails;
  }
};

Int test_sphere_to_ref (const ConstVec3s::HostMirror& p,
                        const ConstIdxs::HostMirror& e,
                        const bool sphere) {
  TestSphereToRefKernel k(p, e, sphere);
  Info info;
  ko::parallel_reduce(k.n(), k, info);
  return info.nfails;
}

} // namespace test
} // namespace sqr
} // namespace siqk
