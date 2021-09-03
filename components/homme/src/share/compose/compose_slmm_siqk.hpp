#ifndef INCLUDE_COMPOSE_SLMM_SIQK_HPP
#define INCLUDE_COMPOSE_SLMM_SIQK_HPP

#include "compose_homme.hpp"

#include <mpi.h>

#ifndef INCLUDE_SIQK_DEFS_HPP
#define INCLUDE_SIQK_DEFS_HPP

#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <sstream>

#ifdef SIQK_TIME
# include <unistd.h>
# include <sys/time.h>
# include <sys/resource.h>
#endif

namespace siqk {
using homme::Int;
using homme::Real;

namespace ko = Kokkos;

#ifndef pr
#define pr(m) do {                                          \
    int _pid_ = 0;                                          \
    MPI_Comm_rank(MPI_COMM_WORLD, &_pid_);                  \
    std::stringstream _ss_;                                 \
    _ss_.precision(15);                                     \
    _ss_ << "slmm: pid " << _pid_ << " " << m << std::endl; \
    std::cerr << _ss_.str();                                \
  } while (0)
#define prc(m) pr(#m << " | " << (m))
#define puf(m)"(" << #m << " " << (m) << ")"
#define pu(m) << " " << puf(m)
template<typename T>
static void prarr (const std::string& name, const T* const v, const size_t n) {
  std::stringstream ss;
  ss << name << ": ";
  for (size_t i = 0; i < n; ++i) ss << " " << v[i];
  pr(ss.str());
}
#define mprarr(m) siqk::prarr(#m, m.data(), m.size())
#endif

#define SIQK_THROW_IF(condition, message) do {                          \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n" << #condition \
        "\nled to the exception\n" << message << "\n";                  \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)

#define SIQK_STDERR_IF(condition, message) do {                   \
    try { SIQK_THROW_IF(condition, message); }                    \
    catch (const std::logic_error& e) { std::cerr << e.what(); }  \
  } while (0)

KOKKOS_INLINE_FUNCTION static void error (const char* const msg)
{ ko::abort(msg); }

KOKKOS_INLINE_FUNCTION static void message (const char* const msg)
{ printf("%s\n", msg); }

typedef ko::LayoutRight Layout;

// SIQK's array types.
typedef ko::View<Real*[3], Layout> Vec3s;
typedef ko::View<const Real*[3], Layout> ConstVec3s;
typedef ko::View<Real*[6], Layout> Vec6s;
typedef ko::View<const Real*[6], Layout> ConstVec6s;
typedef ko::View<Real*[3], ko::LayoutRight, ko::MemoryTraits<ko::Unmanaged> > RawVec3s;
typedef ko::View<const Real*[3], ko::LayoutRight, ko::MemoryTraits<ko::Unmanaged> > RawConstVec3s;
typedef ko::View<Real*, ko::LayoutRight, ko::MemoryTraits<ko::Unmanaged> > RawArray;
typedef ko::View<const Real*, ko::LayoutRight, ko::MemoryTraits<ko::Unmanaged> > RawConstArray;
typedef ko::View<Int**, Layout> Idxs;
typedef ko::View<const Int**, Layout> ConstIdxs;
typedef ko::View<Int*[8], Layout> Nodes;
typedef ko::View<const Int*[8], Layout> ConstNodes;

// Decorator for a View. UnmanagedView<ViewType> gives the same view as
// ViewType, except the memory is unmanaged.
template <typename ViewT>
using UnmanagedView = ko::View<
  typename ViewT::data_type, typename ViewT::array_layout,
  typename ViewT::device_type, ko::MemoryTraits<ko::Unmanaged> >;

template <typename VT> KOKKOS_FORCEINLINE_FUNCTION
typename VT::value_type*
slice (const VT& v, Int i) {
  assert(i >= 0 && i < v.extent_int(0));
  return v.data() + v.extent_int(1)*i;
}
template <typename VT> KOKKOS_FORCEINLINE_FUNCTION
typename VT::const_value_type*
const_slice (const VT& v, Int i) {
  assert(i >= 0 && i < v.extent_int(0));
  return v.data() + v.extent_int(1)*i;
}

// Number of slices in a 2D array, where each row is a slice.
template <typename A2D> KOKKOS_FORCEINLINE_FUNCTION
Int nslices (const A2D& a) { return static_cast<Int>(a.extent_int(0)); }

// Number of entries in a 2D array's row.
template <typename A2D> KOKKOS_FORCEINLINE_FUNCTION
Int szslice (const A2D& a) { return static_cast<Int>(a.extent_int(1)); }

template <typename V, typename CV>
KOKKOS_INLINE_FUNCTION
static void copy (V dst, CV src, const Int n) {
  for (Int i = 0; i < n; ++i) dst[i] = src[i];
}

template <typename DV, typename SV>
void resize_and_copy (DV& d, const SV& s,
                      typename std::enable_if<DV::rank_dynamic == 1>::type* = 0) {
  ko::resize(d, nslices(s));
  ko::deep_copy(d, s);
}

template <typename DV, typename SV>
void resize_and_copy (DV& d, const SV& s,
                      typename std::enable_if<DV::rank_dynamic == 2>::type* = 0) {
  ko::resize(d, nslices(s), szslice(s));
  ko::deep_copy(d, s);
}

template <typename DV, typename SA>
void hm_resize_and_copy (DV& d, const SA& s, const Int n) {
  ko::resize(d, n);
  auto d_hm = ko::create_mirror_view(d);
  for (Int i = 0; i < n; ++i) d_hm[i] = s[i];
  ko::deep_copy(d, d_hm);
}

// GPU-friendly replacements for std::min/max.
template <typename T> KOKKOS_INLINE_FUNCTION
const T& min (const T& a, const T& b) { return a < b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
const T& max (const T& a, const T& b) { return a > b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
void swap (T& a, T&b) {
  T tmp = a;
  a = b;
  b = tmp;
}
template <typename T> KOKKOS_INLINE_FUNCTION constexpr T square (const T& x) { return x*x; }

template<typename T> KOKKOS_INLINE_FUNCTION
T sign (const T& a) { return a > 0 ? 1 : (a < 0 ? -1 : 0); }

} // namespace siqk

#endif // INCLUDE_SIQK_DEFS_HPP

#ifndef INCLUDE_SIQK_GEOMETRY_HPP
#define INCLUDE_SIQK_GEOMETRY_HPP

//#include "siqk_defs.hpp"
//#include "siqk_quadrature.hpp"

namespace siqk {
// Vectors and points are 2D. Thus, if you're working on planes in 3D, project
// to a 2D space before calling these.
struct PlaneGeometry {
  template <typename V> KOKKOS_INLINE_FUNCTION
  static void scale (const Real& a, V v) {
    v[0] *= a; v[1] *= a;
  }
  template <typename CVA, typename CVB> KOKKOS_INLINE_FUNCTION
  static Real dot (const CVA a, const CVB b) {
    return a[0]*b[0] + a[1]*b[1];
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static Real dot_c_amb (const CV c, const CV a, const CV b) {
    return c[0]*(a[0] - b[0]) + c[1]*(a[1] - b[1]);
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void combine (const CV u, const CV v, const Real& a, V x) {
    const Real& oma = 1 - a;
    x[0] = oma*u[0] + a*v[0];
    x[1] = oma*u[1] + a*v[1];
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void axpy (const Real& a, const CV x, V y) {
    y[0] += a*x[0];
    y[1] += a*x[1];
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static Real norm2 (const CV v) {
    return dot(v, v);
  }
  template <typename V> KOKKOS_INLINE_FUNCTION
  static void normalize (V v) {
    scale(1.0/std::sqrt(norm2(v)), v);
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void edge_normal (const CV e1, const CV e2, V en) {
    en[0] = e1[1] - e2[1];
    en[1] = e2[0] - e1[0];
    normalize(en);
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static bool inside (const CV v, const CV e1, const CV en) {
    return dot_c_amb(en, v, e1) >= 0;
  }
};

// All inputs and outputs are relative to the unit-radius sphere. Vectors and
// points are 3D.
struct SphereGeometry {
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void cross (const CV a, const CV b, V c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
  }
  template <typename CVA, typename CVB> KOKKOS_INLINE_FUNCTION
  static Real dot (const CVA a, const CVB b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static Real norm2 (const CV v) {
    return dot(v, v);
  }
  template <typename V> KOKKOS_INLINE_FUNCTION
  static void scale (const Real& a, V v) {
    v[0] *= a; v[1] *= a; v[2] *= a;
  }
  template <typename V> KOKKOS_INLINE_FUNCTION
  static void normalize (V v) {
    scale(1.0/std::sqrt(norm2(v)), v);
  }
  template <typename CVC, typename CVA, typename CVB> KOKKOS_INLINE_FUNCTION
  static Real dot_c_amb (const CVC c, const CVA a, const CVB b) {
    return c[0]*(a[0] - b[0]) + c[1]*(a[1] - b[1]) + c[2]*(a[2] - b[2]);
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void axpy (const Real& a, const CV x, V y) {
    y[0] += a*x[0];
    y[1] += a*x[1];
    y[2] += a*x[2];
  }
  template <typename CVX, typename CVY, typename VZ> KOKKOS_INLINE_FUNCTION
  static void axpbyz (const Real& a, const CVX x, const Real& b, const CVY y,
                      VZ z) {
    z[0] = a*x[0] + b*y[0];
    z[1] = a*x[1] + b*y[1];
    z[2] = a*x[2] + b*y[2];
  }
  template <typename V, typename CV> KOKKOS_INLINE_FUNCTION
  static void copy (V d, const CV s) {
    d[0] = s[0];
    d[1] = s[1];
    d[2] = s[2];
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void combine (const CV u, const CV v, const Real& a, V x) {
    const Real& oma = 1 - a;
    x[0] = oma*u[0] + a*v[0];
    x[1] = oma*u[1] + a*v[1];
    x[2] = oma*u[2] + a*v[2];
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void edge_normal (const CV a, const CV b, V en) {
    cross(a, b, en);
    normalize(en);
  }
  // Is v inside the line anchored at a with inward-facing normal n?
  template <typename CVV, typename CVA, typename CVN> KOKKOS_INLINE_FUNCTION
  static bool inside (const CVV v, const CVA a, const CVN n) {
    return dot_c_amb(n, v, a) >= 0;
  }
};
} // namespace siqk

#endif // INCLUDE_SIQK_GEOMETRY_HPP

#ifndef INCLUDE_SIQK_SQR_HPP
#define INCLUDE_SIQK_SQR_HPP

//#include "siqk_defs.hpp"
//#include "siqk_intersect.hpp"

namespace siqk {
namespace sqr { // spherical quadrilateral <-> reference square
/* Let p be a 3x4 matrix with p(:,i) the i'th vertex in a spherical quad in CCW
   order. Let (a,b) be coordinates in the reference square [0,1]^2. (Here we
   choose [0,1] instead of [-1,1].) (a,b) = (0,0) corresponds to p(:,1); (1,0)
   is p(:,2); (1,1) is p(:,3); (0,1) is p(:,4).
     The map from reference square to bilinear quad can be written
       T = p*[ 1 -1 1 -1
              -1  1 0  0
              -1  0 0  1
               1  0 0  0]';
       f(a,b) = T(:,1)*a*b + T(:,2)*a + T(:,3)*b + T(:,4);
   The map to the sphere is then completed with
       g(a,b) = norm(f(a,b))
       q = f(a,b) / g(a,b).
   The Jacobian matrix for q is given by
       q_a = f_a/g - (f g_a)/g^2
       g_a = g_f f_a
   and similarly for q_b.
*/

namespace impl {
// Compute T(i,:).
template <typename ConstVec3sT, typename Quad>
KOKKOS_INLINE_FUNCTION
void calc_T_row (const ConstVec3sT& p, const Quad& e, const Int i,
                 Real& t1, Real& t2, Real& t3, Real& t4) {
  t4 = p(e[0],i);
  t3 = -t4 + p(e[3],i);
  t2 = -t4 + p(e[1],i);
  t1 = -t2 + p(e[2],i) - p(e[3],i);
}

// Compute T(:,1)*a*b + T(:,2)*a + T(:,3)*b + T(:,4).
template <typename ConstVec3sT, typename Quad>
KOKKOS_INLINE_FUNCTION
void calc_ref_to_bilinear (const ConstVec3sT& p, const Quad& e,
                           Real a, Real b, Real q[3]) {
  a = 0.5*(a + 1);
  b = 0.5*(b + 1);
  for (Int i = 0; i < 3; ++i) {
    Real t1, t2, t3, t4;
    impl::calc_T_row(p, e, i, t1, t2, t3, t4);
    q[i] = t1*a*b + t2*a + t3*b + t4;
  }
}

// The residual function is r(a,b) = f(a,b)/g(a,b) - q.
template <typename ConstVec3sT, typename Quad>
KOKKOS_INLINE_FUNCTION
void calc_residual (const ConstVec3sT& p, const Quad& e, const Real a,
                    const Real b, const Real q[3], Real r[3],
                    const bool sphere = true) {
  calc_ref_to_bilinear(p, e, a, b, r);
  const Real rnorm = sphere ? std::sqrt(SphereGeometry::norm2(r)) : 1;
  for (Int i = 0; i < 3; ++i)
    r[i] = r[i]/rnorm - q[i];
}

// Compute the Jacobian matrix of the residual function: Jacobian(ref square ->
// sphere).
//   TODO Consider rewriting this in terms of the p=1 basis isoparametric
// interpolation formulation. Better performance? See
// calc_isoparametric_jacobian in slmmir.cpp.
template <typename ConstVec3sT, typename Quad>
KOKKOS_INLINE_FUNCTION
void calc_Jacobian (const ConstVec3sT& p, const Quad& e, Real a, Real b,
                    Real J[6], const bool sphere = true) {
  a = 0.5*(a + 1);
  b = 0.5*(b + 1);  
  Real r[3];
  for (Int i = 0; i < 3; ++i) {
    Real t1, t2, t3, t4;
    calc_T_row(p, e, i, t1, t2, t3, t4);
    r[  i] = t1*a*b + t2*a + t3*b + t4;
    J[  i] = t1*b + t2;
    J[3+i] = t1*a + t3;
  }
  if ( ! sphere) return;
  Real rtJ[2] = {0};
  for (Int j = 0; j < 2; ++j) {
    const Real* const Jj = J + 3*j;
    for (Int i = 0; i < 3; ++i)
      rtJ[j] += r[i]*Jj[i];
  }
  const Real rnorm2 = SphereGeometry::norm2(r), rnorm = std::sqrt(rnorm2);
  for (Int j = 0; j < 2; ++j) {
    Real* const Jj = J + 3*j;
    for (Int i = 0; i < 3; ++i)
      Jj[i] = (Jj[i] - r[i]*rtJ[j]/rnorm2)/rnorm;
  }
}

// Solve J dx = r.
KOKKOS_INLINE_FUNCTION
void solve_Jxr (Real J[6], const Real r[3], Real dx[2]) {
  // QR factorization: J -> J [n1 a; 0 n2].
  const Real n1 = std::sqrt(SphereGeometry::norm2(J));
  SphereGeometry::scale(1/n1, J);
  const Real a = SphereGeometry::dot(J, J+3);
  SphereGeometry::axpy(-a, J, J+3);
  const Real n2 = std::sqrt(SphereGeometry::norm2(J+3));
  SphereGeometry::scale(1/n2, J+3);
  // r -> Q' r.
  Real Qtr[2] = {0};
  for (Int j = 0; j < 2; ++j) {
    const Real* const Jj = J + 3*j;
    for (Int i = 0; i < 3; ++i)
      Qtr[j] += Jj[i]*r[i];
  }
  // dx = R \ (Q' r).
  dx[1] = 2*(Qtr[1] / n2);
  dx[0] = 2*((Qtr[0] - a*dx[1]) / n1);
}
} // namespace impl

struct Info {
  bool success;
  Int n_iterations;
};

template <typename ConstVec3sT, typename Quad>
KOKKOS_INLINE_FUNCTION
void calc_ref_to_sphere (
  // The spherical quad containing the point.
  const ConstVec3sT& p, const Quad& e,
  // (a,b) in [-1,1]
  const Real a, const Real b,
  // The point on the sphere.
  Real q[3])
{
  impl::calc_ref_to_bilinear(p, e, a, b, q);
  SphereGeometry::normalize(q);
}

template <typename ConstVec3sT, typename Quad>
KOKKOS_INLINE_FUNCTION
void calc_sphere_to_ref (
  // The spherical quad containing the point.
  const ConstVec3sT& p, const Quad& e,
  // The point on the sphere.
  const Real q[3],
  // (a,b) in [-1,1]
  Real& a, Real& b,
  // Optional info output.
  Info* const info = nullptr,
  // Max number of iterations before returning with failure.
  const Int max_its = 10,
  // Absolute tolerance for Newton iteration.
  const Real tol = 1e2*std::numeric_limits<Real>::epsilon())
{
  const Real tol2 = square(tol);
  Real rnorm2 = 1;
  a = b = 0;
  Int it = 0;
  for (it = 1; it <= max_its; ++it) { // Newton's method.
    Real r[3], J[6];
    impl::calc_residual(p, e, a, b, q, r);
    rnorm2 = SphereGeometry::norm2(r);
    if (rnorm2 <= tol2) break;
    impl::calc_Jacobian(p, e, a, b, J);
    Real dx[2];
    impl::solve_Jxr(J, r, dx);
    a -= dx[0];
    b -= dx[1];
  }
  if (info) {
    info->success = rnorm2 <= tol2;
    info->n_iterations = it;
  }
}

template <typename ConstVec3sT, typename Quad>
KOKKOS_INLINE_FUNCTION
void calc_ref_to_plane (
  // The quad containing the point.
  const ConstVec3sT& p, const Quad& e,
  // (a,b) in [-1,1]
  const Real a, const Real b,
  // The point on the plane.
  Real q[3])
{
  impl::calc_ref_to_bilinear(p, e, a, b, q);
}

template <typename ConstVec3sT, typename Quad>
KOKKOS_INLINE_FUNCTION
void calc_plane_to_ref (
  // The quad containing the point.
  const ConstVec3sT& p, const Quad& e,
  // The point on the plane.
  const Real q[3],
  // (a,b) in [-1,1]
  Real& a, Real& b,
  // Optional info output.
  Info* const info = nullptr,
  // Max number of iterations before returning with failure.
  const Int max_its = 10,
  // Absolute tolerance for Newton iteration.
  const Real tol = 1e2*std::numeric_limits<Real>::epsilon())
{
  const Real tol2 = square(tol);
  Real rnorm2 = 1;
  a = b = 0;
  if (q[2] == 0) {
    // Init assuming a quad in the x-y plane on an axis-aligned, regular grid.
    for (Int d = 0; d < 2; ++d) {
      // One of the two edges is 0 and the other the one we want if the
      // assumption holds.
      const Real den = siqk::max(std::abs(p(e[1], d) - p(e[0], d)),
                                 std::abs(p(e[2], d) - p(e[1], d)));
      if (den == 0) continue;
      const Real ref = 2*(std::abs(q[d] - p(e[0], d))/den) - 1;
      if (d == 0) a = ref; else b = ref;
    }
  }
  Int it = 0;
  for (it = 1; it <= max_its; ++it) { // Newton's method.
    Real r[3], J[6];
    impl::calc_residual(p, e, a, b, q, r, false);
    rnorm2 = SphereGeometry::norm2(r);
    if (rnorm2 <= tol2) break;
    impl::calc_Jacobian(p, e, a, b, J, false);
    Real dx[2];
    impl::solve_Jxr(J, r, dx);
    a -= dx[0];
    b -= dx[1];
  }
  if (info) {
    info->success = rnorm2 <= tol2;
    info->n_iterations = it;
  }
}

// Ref coords, packed (x,y), CCW, starting from (-1,-1).
KOKKOS_INLINE_FUNCTION
const Real* get_ref_vertices () {
  static const Real c[] = {-1, -1, 1, -1, 1, 1, -1, 1};
  return c;
}

namespace test {
Int test_sphere_to_ref(const ConstVec3s::HostMirror& p,
                       const ConstIdxs::HostMirror& e,
                       const bool sphere = true);
} // namespace test
} // namespace sqr
} // namespace siqk

#endif // INCLUDE_SIQK_SQR_HPP

//> end SIQK

#endif
