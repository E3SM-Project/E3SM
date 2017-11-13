//> begin SIQK
// To make this initial integration less messy, inline SIQK in this translation
// unit. Once we've put together the Compose library, I'll remove this code.

#ifndef INCLUDE_SIQK_DEFS_HPP
#define INCLUDE_SIQK_DEFS_HPP

#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <sstream>

#include <Kokkos_Core.hpp>

#ifdef SIQK_TIME
# include <unistd.h>
# include <sys/time.h>
# include <sys/resource.h>
#endif

// Always want this for GPU.
#define SIQK_NONRECURSIVE

namespace siqk {
namespace ko = Kokkos;

#define pr(m) do {                              \
    std::stringstream _ss_;                     \
    _ss_ << m << std::endl;                     \
    std::cerr << _ss_.str();                    \
  } while (0)
#define prc(m) pr(#m << " | " << (m))
#define puf(m)"(" << #m << " " << (m) << ")"
#define pu(m) << " " << puf(m)
template<typename T>
static void prarr (const std::string& name, const T* const v, const size_t n) {
  std::cerr << name << ": ";
  for (size_t i = 0; i < n; ++i) std::cerr << " " << v[i];
  std::cerr << "\n";
}
#define mprarr(m) siqk::prarr(#m, m.data(), m.size())

#define SIQK_THROW_IF(condition, message) do {                          \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n" << #condition \
        "\nled to the exception\n" << message << "\n";                  \
      throw std::logic_error(_ss_.str());                               \
    }                                                                   \
  } while (0)

#define SIQK_STDERR_IF(condition, message) do { \
  try { SIQK_THROW_IF(condition, message); } \
  catch (const std::logic_error& e) { std::cerr << e.what(); } \
} while (0)

KOKKOS_INLINE_FUNCTION static void error (const char* const msg)
{ ko::abort(msg); }

KOKKOS_INLINE_FUNCTION static void message (const char* const msg)
{ printf("%s\n", msg); }

typedef int Int;
typedef double Real;

#ifdef KOKKOS_HAVE_CUDA
typedef ko::LayoutLeft Layout;
#else
typedef ko::LayoutRight Layout;
#endif

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

// Get the host or device version of the array.
template <typename VT, typename ES> struct InExeSpace {
  typedef VT type;
};
template <typename VT> struct InExeSpace<VT, ko::HostSpace> {
  typedef typename VT::HostMirror type;
};

#ifdef KOKKOS_HAVE_CUDA
// A 1D slice of an array.
template <typename VT> KOKKOS_FORCEINLINE_FUNCTION
ko::View<typename VT::value_type*, ko::LayoutStride, typename VT::device_type,
         ko::MemoryTraits<ko::Unmanaged> >
slice (const VT& v, Int i) { return ko::subview(v, i, ko::ALL()); }
// An explicitly const 1D slice of an array.
template <typename VT> KOKKOS_FORCEINLINE_FUNCTION
ko::View<typename VT::const_value_type*, ko::LayoutStride, typename VT::device_type,
         ko::MemoryTraits<ko::Unmanaged> >
const_slice (const VT& v, Int i) { return ko::subview(v, i, ko::ALL()); }
#else
template <typename VT> KOKKOS_FORCEINLINE_FUNCTION
typename VT::value_type*
slice (const VT& v, Int i) { return v.ptr_on_device() + v.dimension_1()*i; }

template <typename VT> KOKKOS_FORCEINLINE_FUNCTION
typename VT::const_value_type*
const_slice (const VT& v, Int i) { return v.ptr_on_device() + v.dimension_1()*i; }
#endif

// Number of slices in a 2D array, where each row is a slice.
template <typename A2D> KOKKOS_FORCEINLINE_FUNCTION
Int nslices (const A2D& a) { return static_cast<Int>(a.dimension_0()); }

// Number of entries in a 2D array's row.
template <typename A2D> KOKKOS_FORCEINLINE_FUNCTION
Int szslice (const A2D& a) { return static_cast<Int>(a.dimension_1()); }

template <typename V, typename CV>
KOKKOS_INLINE_FUNCTION
static void copy (V dst, CV src, const Int n) {
  for (Int i = 0; i < n; ++i) dst[i] = src[i];
}

template <typename DV, typename SV>
void resize_and_copy (DV& d, const SV& s,
                      typename std::enable_if<DV::rank == 1>::type* = 0) {
  ko::resize(d, nslices(s));
  ko::deep_copy(d, s);
}

template <typename DV, typename SV>
void resize_and_copy (DV& d, const SV& s,
                      typename std::enable_if<DV::rank == 2>::type* = 0) {
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

#ifndef INCLUDE_SIQK_QUADRATURE_HPP
#define INCLUDE_SIQK_QUADRATURE_HPP

//#include "siqk_defs.hpp"

namespace siqk {

/* For the TRISYM entries, see, e.g.,
     Triangular quadrature to use for integration Dunavant, D.A. "High Degree
     Efficient Symmetrical Gaussian Quadrature Rules for the Triangle."
     J. Numer. Meth. Eng., 21, pp 1129-1148.
   and
     Zhang, Linbo, Tao Cui, and Hui Liu. "A set of symmetric quadrature rules on
     triangles and tetrahedra." J. of Computational Mathematics (2009): 89-96.
   For the TRITAYLOR, see
     Day, David M. and Mark A. Taylor, "A new 11 point degree 6 cubature formula
     for the triangle", PAMM 7 (2007)
   and
     Taylor, Mark A., Beth A. Wingate, and Len P. Bos. "A cardinal function
     algorithm for computing multivariate quadrature points." SIAM Journal on
     Numerical Analysis 45.1 (2007): 193-205.
*/

#define SIQK_QUADRATURE_TRISYM_ORDER4_COORD                  \
  {0.108103018168070, 0.445948490915965, 0.445948490915965,  \
   0.445948490915965, 0.108103018168070, 0.445948490915965,  \
   0.445948490915965, 0.445948490915965, 0.108103018168070,  \
   0.816847572980458, 0.091576213509771, 0.091576213509771,  \
   0.091576213509771, 0.816847572980458, 0.091576213509771,  \
   0.091576213509771, 0.091576213509771, 0.816847572980458}
#define SIQK_QUADRATURE_TRISYM_ORDER4_WEIGHT                 \
  {0.223381589678011, 0.223381589678011, 0.223381589678011,  \
   0.109951743655322, 0.109951743655322, 0.109951743655322}

#define SIQK_QUADRATURE_TRISYM_ORDER8_COORD                  \
  {0.333333333333333, 0.333333333333333, 0.333333333333333,  \
   0.081414823414554, 0.459292588292723, 0.459292588292723,  \
   0.459292588292723, 0.081414823414554, 0.459292588292723,  \
   0.459292588292723, 0.459292588292723, 0.081414823414554,  \
   0.658861384496480, 0.170569307751760, 0.170569307751760,  \
   0.170569307751760, 0.658861384496480, 0.170569307751760,  \
   0.170569307751760, 0.170569307751760, 0.658861384496480,  \
   0.898905543365938, 0.050547228317031, 0.050547228317031,  \
   0.050547228317031, 0.898905543365938, 0.050547228317031,  \
   0.050547228317031, 0.050547228317031, 0.898905543365938,  \
   0.008394777409958, 0.263112829634638, 0.728492392955404,  \
   0.008394777409958, 0.728492392955404, 0.263112829634638,  \
   0.263112829634638, 0.008394777409958, 0.728492392955404,  \
   0.263112829634638, 0.728492392955404, 0.008394777409958,  \
   0.728492392955404, 0.263112829634638, 0.008394777409958,  \
   0.728492392955404, 0.008394777409958, 0.263112829634638}
#define SIQK_QUADRATURE_TRISYM_ORDER8_WEIGHT                 \
  {0.144315607677787, 0.095091634267285, 0.095091634267285,  \
   0.095091634267285, 0.103217370534718, 0.103217370534718,  \
   0.103217370534718, 0.032458497623198, 0.032458497623198,  \
   0.032458497623198, 0.027230314174435, 0.027230314174435,  \
   0.027230314174435, 0.027230314174435, 0.027230314174435,  \
   0.027230314174435}

#define SIQK_QUADRATURE_TRITAY_ORDER12_COORD \
  {7.26510255160501828e-02, 9.27348974483949817e-01, 0.00000000000000000e+00, \
   2.11790731803609689e-02, 2.35517332495786824e-02, 9.55269193570060349e-01, \
   1.41841115784669236e-01, 5.40914911362088088e-17, 8.58158884215330708e-01, \
   1.15143666726236216e-02, 9.45475073220970907e-01, 4.30105601064054710e-02, \
   2.77555756156289135e-17, 1.54064601626856063e-01, 8.45935398373143910e-01, \
   3.72684680767588483e-01, -1.88694080537681499e-16, 6.27315319232411683e-01, \
   9.43134911146902510e-01, 2.71109713562557482e-02, 2.97541174968417414e-02, \
   8.44725347421859452e-01, 1.46044961672175677e-01, 9.22969090596487129e-03, \
   8.23277107647898521e-01, 2.11522233831219000e-02, 1.55570668968979586e-01, \
   6.21586880750877868e-01, 1.45665147883470222e-02, 3.63846604460775103e-01, \
   2.21919501597089841e-02, 7.88601719223131714e-01, 1.89206330617159302e-01, \
   2.27722111443204644e-01, 7.49189739790679599e-01, 2.30881487661157569e-02, \
   7.38137544226065284e-02, 7.18714961015890358e-02, 8.54314749475804436e-01, \
   6.43364629415364875e-01, 3.32129083947645065e-01, 2.45062866369900600e-02, \
   2.28091126376529507e-02, 3.61181591189672080e-01, 6.16009296172674969e-01, \
   6.63093778446759319e-01, 2.43458133948799671e-01, 9.34480876044410103e-02, \
   2.51456820638045198e-02, 5.81689214740147453e-01, 3.93165103196048027e-01, \
   4.29837040104380730e-01, 5.44446676271925334e-01, 2.57162836236939363e-02, \
   9.40413011410586863e-02, 8.26003314017559997e-01, 7.99553848413813162e-02, \
   7.94010795132135239e-01, 1.16386499067277244e-01, 8.96027058005875177e-02, \
   7.83496599417470019e-02, 2.03768481077729741e-01, 7.17881858980523258e-01, \
   2.25505520049374242e-01, 6.44132203822605637e-02, 7.10081259568365097e-01, \
   6.43800731623786371e-01, 9.54285858105846096e-02, 2.60770682565629019e-01, \
   5.43837635808460451e-01, 2.44982965093490213e-01, 2.11179399098049336e-01, \
   4.32112641877997194e-01, 7.05667243440369213e-02, 4.97320633777965815e-01, \
   2.55495747579340349e-01, 6.19381257362555782e-01, 1.25122995058103870e-01, \
   1.22162380966293838e-01, 6.27682615680314027e-01, 2.50155003353392136e-01, \
   4.47861373562203791e-01, 4.22605657433460014e-01, 1.29532969004336196e-01, \
   4.09354529674576528e-01, 2.10785259391403995e-01, 3.79860210934019449e-01, \
   1.24718320885524481e-01, 4.08963804491244809e-01, 4.66317874623230710e-01, \
   2.28197277938737758e-01, 2.13777432530059680e-01, 5.58025289531202562e-01, \
   2.88796329020881648e-01, 4.09786577770025306e-01, 3.01417093209092990e-01}

#define SIQK_QUADRATURE_TRITAY_ORDER12_WEIGHT                           \
  {4.888049814660050e-03, 6.675900027367356e-03, 6.845534654343699e-03, \
   7.119751436080721e-03, 7.714492373624846e-03, 9.654708742436301e-03, \
   1.050932673560249e-02, 1.068084365762828e-02, 1.848368581123072e-02, \
   1.854548042160657e-02, 2.062000411968213e-02, 2.168508541701153e-02, \
   2.249074619915818e-02, 2.490407320150775e-02, 2.509917342768508e-02, \
   2.794373431987983e-02, 2.814555860521331e-02, 2.816965445973000e-02, \
   3.052917241207244e-02, 3.057527760403899e-02, 3.957360579297199e-02, \
   4.128188739546268e-02, 4.593784216579169e-02, 4.749957532530720e-02, \
   4.814880503690738e-02, 5.096492487678762e-02, 5.335208304882109e-02, \
   5.414687261316752e-02, 5.943783395113540e-02, 5.998970732710617e-02, \
   6.316454642265663e-02, 7.522206260332436e-02}

#define SIQK_QUADRATURE_TRITAY_ORDER6_COORD                             \
  {4.724686653264358e-02, 5.725498667747682e-02, 8.954981467898796e-01, \
   4.280913872509884e-02, 8.953626400245792e-01, 6.182822125032195e-02, \
   2.921805130458027e-01, 6.844757484565146e-01, 2.334373849768268e-02, \
   8.712234683377076e-01, 6.874625591502949e-02, 6.003027574726293e-02, \
   5.086198608278325e-02, 6.156762055758400e-01, 3.334618083413767e-01, \
   2.128646728100595e-01, 6.279461411977890e-01, 1.591891859921515e-01, \
   2.817957679526839e-01, 6.290913834186361e-02, 6.552950937054525e-01, \
   6.225041026512227e-01, 6.837821192050995e-02, 3.091176854282673e-01, \
   7.604403244598745e-02, 2.875294583743921e-01, 6.364265091796204e-01, \
   5.941924379444020e-01, 3.287835564131346e-01, 7.702400564246337e-02, \
   3.353648085404556e-01, 3.122904050136449e-01, 3.523447864458995e-01}

#define SIQK_QUADRATURE_TRITAY_ORDER6_WEIGHT                            \
  {3.806807185295551e-02, 3.837935530775279e-02, 4.620045674456197e-02, \
   5.346758944419899e-02, 8.375582696574595e-02, 1.016448330255167e-01, \
   1.018615244613670e-01, 1.114218316600018e-01, 1.120094502629461e-01, \
   1.247875714375583e-01, 1.884034888373949e-01}

class TriangleQuadrature {
  const Real trisym_order4_coord_  [ 18] = SIQK_QUADRATURE_TRISYM_ORDER4_COORD;
  const Real trisym_order4_weight_ [  6] = SIQK_QUADRATURE_TRISYM_ORDER4_WEIGHT;
  const Real tritay_order6_coord_  [ 33] = SIQK_QUADRATURE_TRITAY_ORDER6_COORD;
  const Real tritay_order6_weight_ [ 11] = SIQK_QUADRATURE_TRITAY_ORDER6_WEIGHT;
  const Real trisym_order8_coord_  [ 48] = SIQK_QUADRATURE_TRISYM_ORDER8_COORD;
  const Real trisym_order8_weight_ [ 16] = SIQK_QUADRATURE_TRISYM_ORDER8_WEIGHT;
  const Real tritay_order12_coord_ [ 96] = SIQK_QUADRATURE_TRITAY_ORDER12_COORD;
  const Real tritay_order12_weight_[ 32] = SIQK_QUADRATURE_TRITAY_ORDER12_WEIGHT;

public:
  KOKKOS_INLINE_FUNCTION TriangleQuadrature () {}

  KOKKOS_INLINE_FUNCTION
  void get_coef (const int order, RawConstVec3s& coord,
                 RawConstArray& weight) const {
    switch (order) {
    case 4:
      coord = RawConstVec3s(trisym_order4_coord_, 6, 3);
      weight = RawConstArray(trisym_order4_weight_, 6);
      break;
    case 6:
      coord = RawConstVec3s(tritay_order6_coord_, 11, 3);
      weight = RawConstArray(tritay_order6_weight_, 11);
      break;
    case 8:
      coord = RawConstVec3s(trisym_order8_coord_, 16, 3);
      weight = RawConstArray(trisym_order8_weight_, 16);
      break;
    case 12:
      coord = RawConstVec3s(tritay_order12_coord_, 32, 3);
      weight = RawConstArray(tritay_order12_weight_, 32);
      break;
    default:
      ko::abort("TriangleQuadrature::get_coef: order not supported.");
    }
  }
};

} // namespace siqk

#endif // INCLUDE_SIQK_QUADRATURE_HPP

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

  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void edge_normal (const CV e1, const CV e2, V en) {
    en[0] = e1[1] - e2[1];
    en[1] = e2[0] - e1[0];
  }

  template <typename CV> KOKKOS_INLINE_FUNCTION
  static bool inside (const CV v, const CV e1, const CV en) {
    return dot_c_amb(en, v, e1) >= 0;
  }

  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void intersect (const CV v1, const CV v2, const CV e1, const CV en,
                         V intersection) {
    Real a; {
      const Real
        num = dot_c_amb(en, e1, v1),
        den = dot_c_amb(en, v2, v1);
      a = num == 0 || den == 0 ? 0 : num/den;
      a = a < 0 ? 0 : a > 1 ? 1 : a;
    }
    combine(v1, v2, a, intersection);
  }

  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static bool output (const CV v, Int& no, const V vo) {
#ifdef SIQK_DEBUG
    if (no >= nslices(vo)) {
      std::stringstream ss;
      ss << "output: No room in vo; vo.n() is " << nslices(vo) << " but no is "
         << no << "\n";
      message(ss.str().c_str());
    }
#endif
    if (no >= nslices(vo)) return false;
    vo(no,0) = v[0];
    vo(no,1) = v[1];
    ++no;
    return true;
  }

  //todo Handle non-convex case.
  template <typename CV2s>
  KOKKOS_INLINE_FUNCTION
  static Real calc_area (const TriangleQuadrature& , const CV2s& v,
                         const Int n) {
    return calc_area_formula(v, n);
  }

  template <typename CV2s>
  KOKKOS_INLINE_FUNCTION
  static Real calc_area_formula (const CV2s& v, const Int n) {
    Real area = 0;
    for (Int i = 1, ilim = n - 1; i < ilim; ++i)
      area += calc_tri_jacobian(slice(v,0), slice(v,i), slice(v,i+1));
    return 0.5*area;
  }

  template <typename CV, typename CA>
  KOKKOS_INLINE_FUNCTION
  static void bary2coord (const CV v1, const CV v2, const CV v3, const CA alpha,
                          Real u[2]) {
    for (Int k = 0; k < 2; ++k) u[k] = 0;
    axpy(alpha[0], v1, u);
    axpy(alpha[1], v2, u);
    axpy(alpha[2], v3, u);
  }

  template <typename CV>
  KOKKOS_INLINE_FUNCTION
  static Real calc_tri_jacobian (const CV v1, const CV v2, const CV v3) {
    Real r1[2], r2[2];
    r1[0] = v2[0] - v1[0];
    r1[1] = v2[1] - v1[1];
    r2[0] = v3[0] - v1[0];
    r2[1] = v3[1] - v1[1];
    const Real a = r1[0]*r2[1] - r1[1]*r2[0];
    return a;
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
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static Real dot (const CV a, const CV b) {
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
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static Real dot_c_amb (const CV c, const CV a, const CV b) {
    return c[0]*(a[0] - b[0]) + c[1]*(a[1] - b[1]) + c[2]*(a[2] - b[2]);
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void axpy (const Real& a, const CV x, V y) {
    y[0] += a*x[0];
    y[1] += a*x[1];
    y[2] += a*x[2];
  }
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void axpbyz (const Real& a, const CV x, const Real& b, const CV y,
                      V z) {
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
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static bool inside (const CV v, const CV a, const CV n) {
    return dot_c_amb(n, v, a) >= 0;
  }

  /* Let
       en = edge normal
       e1 = edge starting point
       d = en' e1
       v(a) = (1 - a) v1 + a v2.
     Solve n' v = d for a:
       a = (en' (e1 - v1)) / (en' (v2 - v1)).
     Then uvec(v(a)) is the intersection point on the unit sphere. Assume
     intersection exists. (Already filtered by 'inside'.)
  */
  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static void intersect (const CV v1, const CV v2, const CV e1, const CV en,
                         V intersection) {
    Real a; {
      const Real
        num = dot_c_amb(en, e1, v1),
        den = dot_c_amb(en, v2, v1);
      a = num == 0 || den == 0 ? 0 : num/den;
      a = a < 0 ? 0 : a > 1 ? 1 : a;
    }
    combine(v1, v2, a, intersection);
    normalize(intersection);
  }

  template <typename CV, typename V> KOKKOS_INLINE_FUNCTION
  static bool output (const CV v, Int& no, V vo) {
#ifdef SIQK_DEBUG
    if (no >= nslices(vo)) {
      std::stringstream ss;
      ss << "output: No room in vo; vo.n() is " << nslices(vo) << " but no is "
         << no << "\n";
      message(ss.str().c_str());
    }
#endif
    if (no >= nslices(vo)) return false;
    vo(no,0) = v[0];
    vo(no,1) = v[1];
    vo(no,2) = v[2];
    ++no;
    return true;
  }

  //todo Handle non-convex case.
  // This uses a terrible formula, but it's just for testing.
  template <typename CV3s>
  KOKKOS_INLINE_FUNCTION
  static Real calc_area_formula (const CV3s& v, const Int n) {
    Real area = 0;
    for (Int i = 1, ilim = n - 1; i < ilim; ++i) {
      const Real a = calc_arc_length(slice(v,0), slice(v,i));
      const Real b = calc_arc_length(slice(v,i), slice(v,i+1));
      const Real c = calc_arc_length(slice(v,i+1), slice(v,0));
      const Real s = 0.5*(a + b + c);
      const Real d = (std::tan(0.5*s)*std::tan(0.5*(s-a))*
                      std::tan(0.5*(s-b))*std::tan(0.5*(s-c)));
      if (d <= 0) continue;
      area += 4*std::atan(std::sqrt(d));
    }
    return area;
  }
  template <typename CV> KOKKOS_INLINE_FUNCTION
  static Real calc_arc_length (const CV a, const CV b) {
    const Real d = dot(a, b);
    if (d >= 1) return 0;
    return acos(d);
  }

  template <typename CV3s>
  KOKKOS_INLINE_FUNCTION
  static Real calc_area (const TriangleQuadrature& q, const CV3s& v,
                         const Int n) {
    Real area = 0, u[3];
    for (Int i = 1, ilim = n - 1; i < ilim; ++i) {
      Real a = 0;
      RawConstVec3s coord;
      RawConstArray weight;
      q.get_coef(8, coord, weight);
      for (Int k = 0, klim = nslices(coord); k < klim; ++k) {
        const Real jac = calc_tri_jacobian(slice(v,0), slice(v,i), slice(v,i+1),
                                           slice(coord, k), u);
        a += weight[k]*jac;
      }
      area += 0.5*a;
    }
    return area;
  }

  template <typename CV, typename CA>
  KOKKOS_INLINE_FUNCTION
  static Real calc_tri_jacobian (const CV v1, const CV v2, const CV v3,
                                 const CA alpha, Real u[3]) {
    // V(:,i) is vertex i of the spherical triangle on the unit sphere. The
    // coefs
    //     alpha = [a1, a2, 1 - a1 - a2]'
    //           = [1 0; 0 1; -1 -1] [a1, a2]'
    //           = alpha_a a
    // (barycentric coords) give the location
    //     v = V alpha
    // on the planar triangle, and u = uvec(v) is the point on the unit sphere.
    //   For a planar tri in 3D, the jacobian is
    //     v_a = v_alpha alpha_a
    //         = V [1 0; 0 1; -1 -1]
    //     J = norm(cross(v_a(:,1), v_a(:,2))).
    // For a spherical tri with the same vertices,
    //     u = v/(v' v)^{1/2}
    //     u_a = u_alpha alpha_a
    //         = (v'v)^{-1/2} (I - u u') V alpha_a
    //         = (v'v)^{-1/2} (I - u u') v_a
    //     J = norm(cross(u_a(:,1), u_a(:,2))).
    for (Int k = 0; k < 3; ++k) u[k] = 0;
    axpy(alpha[0], v1, u);
    axpy(alpha[1], v2, u);
    axpy(alpha[2], v3, u);
    const auto oovn = 1/std::sqrt(norm2(u));
    scale(oovn, u);
    Real u_a[3][3];
    axpbyz(1, v1, -1, v3, u_a[0]);
    axpbyz(1, v2, -1, v3, u_a[1]);
    for (int i = 0; i < 2; ++i) {
      axpy(-dot(u, u_a[i]), u, u_a[i]);
      scale(oovn, u_a[i]);
    }
    cross(u_a[0], u_a[1], u_a[2]);
    return std::sqrt(norm2(u_a[2]));
  }
};

} // namespace siqk

#endif // INCLUDE_SIQK_GEOMETRY_HPP

#ifndef INCLUDE_SIQK_SEARCH_HPP
#define INCLUDE_SIQK_SEARCH_HPP

//#include "siqk_defs.hpp"
//#include "siqk_geometry.hpp"
#include <vector>

namespace siqk {

// Oct-tree. Might do something else better suited to the sphere later.
template <typename Geo, Int max_depth_ = 10>
class Octree {
public:
  enum { max_depth = max_depth_ };
  typedef Real BoundingBox[6];

  struct Options {
    // Do not go beyond max_depth_ depth, including the root and leaf. With this
    // constraInt, try to go deep enough so that a leaf has no more than
    // max_nelem elements.
    Int max_nelem;
    Options () : max_nelem(8) {}
  };

  // Bounding box for a cluster of points ps (possibly vertices).
  template <typename CV3s, typename BB>
  static void calc_bb (const CV3s& ps, const Int np, BB bb) {
    if (np == 0) return;
    for (Int j = 0; j < 3; ++j)
      bb[j] = bb[j+3] = ps(0,j);
    for (Int i = 1; i < np; ++i) {
      for (Int j = 0; j < 3; ++j) {
        bb[j] = min(bb[j], ps(i,j));
        bb[j+3] = max(bb[j+3], ps(i,j));
      }
    }
    pad_bb(bb);
  }

  template <typename CV3s, typename CIV, typename BB>
  KOKKOS_INLINE_FUNCTION
  static void calc_bb (const CV3s& p, const CIV e, const Int ne, BB ebb) {
    for (Int j = 0; j < 3; ++j)
      ebb[j] = ebb[j+3] = p(e[0], j);
    for (Int i = 1; i < ne; ++i) {
      if (e[i] == -1) break;
      for (Int j = 0; j < 3; ++j) {
        ebb[j] = min(ebb[j], p(e[i], j));
        ebb[j+3] = max(ebb[j+3], p(e[i], j));
      }
    }
    pad_bb(ebb);
  }

  // If a bounding box was constructed from vertices of a spherical polygon,
  // expand it to account for the possible protrusion of the sphere.
  template <typename BB>
  KOKKOS_INLINE_FUNCTION
  static void pad_bb (BB bb) {
    if (std::is_same<Geo, PlaneGeometry>::value) return;
    Real hl = 0.5*std::sqrt(square(bb[3] - bb[0]) + square(bb[4] - bb[1]) +
                            square(bb[5] - bb[2]));
    // Limit the half-length to the circle's radius.
    hl = min(1.0, hl);
    // Max distance from a chord of length 2 hl to the unit circle:
    //     hl = sin theta
    //    pad = 1 - cos theta = 1 - sqrt(1 - sin^2 theta) = 1 - sqrt(1 - hl^2).
    const Real pad = 1 - std::sqrt(1 - square(hl));
    for (Int i = 0; i < 3; ++i) bb[  i] -= pad;
    for (Int i = 0; i < 3; ++i) bb[3+i] += pad;
  }

  template <typename CV3s>
  static void calc_bb (const CV3s& ps, BoundingBox bb) {
    calc_bb(ps, nslices(ps), bb);
  }

  template <typename CV3s, typename CIs, typename V6s>
  static void calc_bb (const CV3s& p, const CIs& e, V6s& ebbs) {
    assert(nslices(ebbs) == nslices(e));
    for (Int k = 0, klim = nslices(e); k < klim; ++k)
      calc_bb(p, slice(e, k), szslice(e), slice(ebbs, k));
  }

  // p is a 3xNp array of points. e is a KxNe array of elements. An entry <0 is
  // ignored. All <0 entries must be at the end of an element's list.
  Octree (const ConstVec3s::HostMirror& p, const ConstIdxs::HostMirror& e,
          const Options& o) {
    init(p, e, o);
  }
  Octree (const ConstVec3s::HostMirror& p, const ConstIdxs::HostMirror& e) {
    Options o;
    init(p, e, o);
  }

  Octree() {}
  void init (const ConstVec3s::HostMirror& p, const ConstIdxs::HostMirror& e) {
    Options o;
    init(p, e, o);
  }

  // Apply f to every element in leaf nodes with which bb overlaps. f must have
  // function
  //     void operator(const Int element).
  template <typename CV, typename Functor>
  KOKKOS_INLINE_FUNCTION
  void apply (const CV bb, Functor& f) const {
    if (nslices(nodes_) == 0) {
      for (Int i = 0; i < offsets_[1]; ++i)
        f(elems_[i]);
      return;
    }
#ifdef SIQK_NONRECURSIVE
    // Non-recursive impl.
    {
      // Stack.
      Real snbb[8*max_depth_];
      Int sni[max_depth_], si[max_depth_];
      Int sp = 0;
      // Args for top-level call.
      copy(snbb, bb_, 8);
      sni[sp] = 0;
      si[sp] = 0;
      while (sp >= 0) {
        // Get stack frame's (nbb, ni, current i) values.
        const Int i = si[sp];
        if (i == 8) {
          --sp;
          continue;
        }
        // Increment stored value of i for next iteration. Current value is
        // stored in 'i' above.
        ++si[sp];
        const Int ni = sni[sp];
        const Real* const nbb = snbb + 8*sp;
        // Can use the next stack frame's bb space for a child bb.
        Real* const child_bb = snbb + 8*(sp+1);
        fill_child_bb(nbb, i, child_bb);
        if ( ! do_bb_overlap(child_bb, bb)) continue;
        Int e = nodes_(ni,i);
        if (e < 0) {
          // Leaf, so apply functor to each element.
          e = std::abs(e + 1);
          for (Int k = offsets_[e]; k < offsets_[e+1]; ++k)
            f(elems_[k]);
        } else if (e > 0) {
          // Recurse.
          ++sp;
          sni[sp] = e;
          si[sp] = 0;
        }
      }
    }
#else
    apply_r(0, bb_, bb, f);
#endif
  }

private:
  /* Each node in the oct-tree contains 8 integers, stored in 'nodes'.

     >0 is an index Into 'nodes', pointing to a child node.

     A <=0 entry in 'nodes' indicates a leaf node. If 0, there are no elements
     in the leaf. If <0, the negative of the entry minus 1 is the index of an
     offset array indexing 'elems'.

     Each segment of 'elems' contains a list of element indices covered by a
     leaf node. Element indices refer to the list of elements the caller
     provides during oct-tree construction.
  */

  // Static data structures holding the completed octree.
  //   nodes(:,i) is a list. The list includes children of node i (>0) and leaf
  // node data (<=0).
  //todo Make these const once ready to do full GPU stuff.
  Nodes nodes_;
  // A leaf node corresponding to -k covers elements
  //     elems[offset[k] : offset[k]-1].
  ko::View<int*> offsets_, elems_;
  // Root node's bounding box.
  BoundingBox bb_;

  // Dynamic data structures for construction phase.
  class IntList {
    Int* const buf_;
    Int i_;
  public:
    IntList (Int* const buf) : buf_(buf), i_(0) {}
    void reset () { i_ = 0; }
    void push (const Int& i) { buf_[i_++] = i; }
    Int* data () { return buf_; }
    Int n () const { return i_; }
    const Int& operator[] (const Int& i) const { return buf_[i]; }
  };

  class DynIntList {
    std::vector<Int> buf_;
  public:
    DynIntList () {}
    void push (const Int& i) { buf_.push_back(i); }
    Int& back () { return buf_.back(); }
    Int& operator[] (const size_t i) {
      if (i >= buf_.size())
        buf_.resize(i+1);
      return buf_[i];
    }
    const Int& operator[] (const size_t i) const { return buf_[i]; }
    Int n () const { return static_cast<Int>(buf_.size()); }
    const Int* data () const { return buf_.data(); }
  };

  // Opposite index slot convention.
  class DynNodes {
    std::vector<Int> buf_;
  public:
    Int n () const { return static_cast<Int>(buf_.size()) >> 3; }
    const Int* data () const { return buf_.data(); }
    Int& operator() (const Int& r, const Int& c) {
      const size_t ec = (c+1) << 3;
      if (ec >= buf_.size())
        buf_.resize(ec);
      return const_cast<Int&>(
        const_cast<const DynNodes*>(this)->operator()(r, c));
    }
    const Int& operator() (const Int& r, const Int& c) const {
      assert(((c << 3) + r) >= 0);
      assert(((c << 3) + r) < (Int) buf_.size());
      return buf_[(c << 3) + r];
    }
  };

  void init (const ConstVec3s::HostMirror& p, const ConstIdxs::HostMirror& e,
             const Options& o) {
    if (nslices(e) == 0) return;
    // Get OT's bounding box.
    calc_bb(p, bb_);
    // Get elements' bounding boxes.
    Vec6s::HostMirror ebbs("ebbs", nslices(e), 6);
    calc_bb(p, e, ebbs);
    // Static element lists for work. Each level has active work space.
    std::vector<Int> buf(max_depth_*nslices(e));
    IntList es(buf.data()), wrk(buf.data() + nslices(e));
    for (Int i = 0, ilim = nslices(e); i < ilim; ++i)
      es.push(i);
    // Dynamic element lists.
    DynIntList offsets, elems;
    offsets[0] = 0;
    // Dynamic node data structure.
    DynNodes nodes;
    // Recurse. We don't care about the return value. If it's 0 and nodes.n() ==
    // 0, we'll detect as much in 'apply'.
    init_r(1, bb_, ebbs, o, es, wrk, offsets, elems, nodes);
    // Build the static data structures.
    if (elems.n() == 0) return;
    init_static_ds(nodes, offsets, elems);
  }

  Int init_r (const Int depth, // Tree's depth at this point, including root.
              const BoundingBox& nbb, // My bounding box.
              const ConstVec6s::HostMirror& ebbs, // All elements' bounding boxes.
              const Options& o, // Options controlling construct of the tree.
              IntList& es, // List of elements in my bounding box.
              IntList& wrk, // Work space to store working element lists.
              DynIntList& offsets, // Offsetss Into elems.
              DynIntList& elems, // Elements belonging to leaf nodes.
              DynNodes& nodes) // Dynamic nodes data structure.
  {
    const Int my_idx = nodes.n(); // My node index.
    // Decide what to do.
    if (es.n() == 0) {
      // I have no elements, so return 0 to indicate I'm a leaf node containing
      // nothing.
      return 0;
    } else if (es.n() <= o.max_nelem || depth == max_depth_) {
      // I'm a leaf node with elements. Store my list of elements and return the
      // storage location.
      const Int os = offsets.back();
      offsets.push(os + es.n());
      for (Int i = 0, n = es.n(); i < n; ++i)
        elems[os + i] = es[i];
      return 1 - offsets.n();
    } else {
      // I'm not a leaf node.
      nodes(0, my_idx) = 0; // Insert myself Into the nodes array.
      for (Int ic = 0; ic < 8; ++ic) {
        BoundingBox child_bb;
        fill_child_bb(nbb, ic, child_bb);
        // Find the elements that are in this child's bb.
        IntList ces(wrk.data());
        for (Int i = 0, n = es.n(); i < n; ++i)
          if (do_bb_overlap(child_bb, slice(ebbs, es[i])))
            ces.push(es[i]);
        // Create some work space.
        IntList cwrk(wrk.data() + ces.n());
        // Recurse.
        const Int child_idx = init_r(depth+1, child_bb, ebbs, o, ces, cwrk,
                                     offsets, elems, nodes);
        nodes(ic, my_idx) = child_idx;
      }
      return my_idx;
    }
  }

  void init_static_ds (const DynNodes nodes, const DynIntList& offsets,
                       const DynIntList& elems) {
    {
      ko::resize(nodes_, nodes.n(), 8);
      auto nodes_hm = ko::create_mirror_view(nodes_);
      for (Int i = 0; i < nodes.n(); ++i)
        for (Int j = 0; j < 8; ++j)
          nodes_hm(i,j) = nodes(j,i);
      ko::deep_copy(nodes_, nodes_hm);
    }
    hm_resize_and_copy(offsets_, offsets, offsets.n());
    hm_resize_and_copy(elems_, elems, elems.n());
  }

  // Using parent bb p, fill child bb c, with child_idx in 0:7.
  template <typename CBB, typename BB>
  KOKKOS_INLINE_FUNCTION
  static void fill_child_bb (const CBB& p, const Int& child_idx, BB& c) {
    const Real m[] = { 0.5*(p[0] + p[3]),
                         0.5*(p[1] + p[4]),
                         0.5*(p[2] + p[5]) };
    switch (child_idx) {
    case 0: c[0] = p[0]; c[1] = p[1]; c[2] = p[2]; c[3] = m[0]; c[4] = m[1]; c[5] = m[2]; break;
    case 1: c[0] = m[0]; c[1] = p[1]; c[2] = p[2]; c[3] = p[3]; c[4] = m[1]; c[5] = m[2]; break;
    case 2: c[0] = m[0]; c[1] = m[1]; c[2] = p[2]; c[3] = p[3]; c[4] = p[4]; c[5] = m[2]; break;
    case 3: c[0] = p[0]; c[1] = m[1]; c[2] = p[2]; c[3] = m[0]; c[4] = p[4]; c[5] = m[2]; break;
    case 4: c[0] = p[0]; c[1] = p[1]; c[2] = m[2]; c[3] = m[0]; c[4] = m[1]; c[5] = p[5]; break;
    case 5: c[0] = m[0]; c[1] = p[1]; c[2] = m[2]; c[3] = p[3]; c[4] = m[1]; c[5] = p[5]; break;
    case 6: c[0] = m[0]; c[1] = m[1]; c[2] = m[2]; c[3] = p[3]; c[4] = p[4]; c[5] = p[5]; break;
    case 7: c[0] = p[0]; c[1] = m[1]; c[2] = m[2]; c[3] = m[0]; c[4] = p[4]; c[5] = p[5]; break;
    default:
      // impossible
      error("fill_child_bb: The impossible has happened.");
    }
  }

  // Do bounding boxes a and b overlap?
  template <typename BB>
  KOKKOS_INLINE_FUNCTION
  static bool do_bb_overlap (const BoundingBox a, const BB b) {
    for (Int i = 0; i < 3; ++i)
      if ( ! do_lines_overlap(a[i], a[i+3], b[i], b[i+3]))
        return false;
    return true;
  }

  KOKKOS_INLINE_FUNCTION
  static bool do_lines_overlap (const Real& a1, const Real& a2,
                                const Real& b1, const Real& b2) {
    return ! (a2 < b1 || a1 > b2);
  }

  template <typename CV, typename Functor> KOKKOS_INLINE_FUNCTION
  void apply_r (const Int ni, const BoundingBox& nbb, const CV bb,
                Functor& f) const {
    for (Int i = 0; i < 8; ++i) {
      BoundingBox child_bb;
      fill_child_bb(nbb, i, child_bb);
      if ( ! do_bb_overlap(child_bb, bb)) continue;
      Int e = nodes_(ni,i);
      if (e > 0)
        apply_r(e, child_bb, bb, f);
      else if (e < 0) {
        e = std::abs(e + 1);
        for (Int k = offsets_[e]; k < offsets_[e+1]; ++k)
          f(elems_[k]);
      }
    }
  }
};

} // namespace siqk

#endif // INCLUDE_SIQK_SEARCH_HPP

#ifndef INCLUDE_SIQK_INTERSECT_HPP
#define INCLUDE_SIQK_INTERSECT_HPP

//#include "siqk_defs.hpp"
//#include "siqk_geometry.hpp"
//#include "siqk_search.hpp"
//#include "siqk_quadrature.hpp"

namespace siqk {

// Sutherland-Hodgmann polygon clipping algorithm. Follow Foley, van Dam,
// Feiner, Hughes Fig 3.49.
namespace sh {
/* A mesh is described by the following arrays:
       p: 3 x #nodes, the array of vertices.
       e: max(#verts) x #elems, the array of element base-0 indices.
       nml: 3 x #edges, the array of edge normals.
       en: max(#verts) x #elems, the array of edge-normal base-0 indices.
     e. e indexes p. e(i,j) == -1 in column j indicates that j:end are not used.
     nml. As a mesh is refined, cancellation error makes an edge normal based
   off of an element's vertices increasingly inaccurate. Roughly, if an edge
   subtends angle phi of the sphere, -log10(phi/(2 pi)) digits are lost in the
   edge normal. Therefore, we compute edge normals offline, since in certain
   meshes, they can be computed by an accurate means. E.g., in a cubed-sphere
   mesh, the whole line of a square face can be used to compute the edge
   normal. Furthermore, there are far fewer unique edge normals than edges.
 */
template <typename ES = ko::DefaultExecutionSpace>
struct Mesh {
  typedef typename InExeSpace<ConstVec3s, ES>::type RealArray;
  typedef typename InExeSpace<ConstIdxs, ES>::type IntArray;

  RealArray p, nml;
  IntArray e, en;

  Mesh () {}

  Mesh (const Mesh<ko::HostSpace>& m) {
    typename InExeSpace<Vec3s, ES>::type tp, tnml;
    typename InExeSpace<Idxs, ES>::type te, ten;
    resize_and_copy(tp, m.p); p = tp;
    resize_and_copy(tnml, m.nml); nml = tnml;
    resize_and_copy(te, m.e); e = te;
    resize_and_copy(ten, m.en); en = ten;
  }
};

// Generally not a user routine.
template <typename geo, typename CV3s, typename V3s, typename CV>
KOKKOS_INLINE_FUNCTION
bool clip_against_edge (
  // Input vertex list.
  const CV3s& vi, const Int ni,
  // Output vertex list.
  V3s& vo, Int& no,
  // One point of the clip edge.
  const CV ce1,
  // Clip edge's inward-facing normal.
  const CV cen)
{
  Real intersection[3];
  no = 0;
  auto s = const_slice(vi, ni-1);
  for (Int j = 0; j < ni; ++j) {
    auto p = const_slice(vi,j);
    if (geo::inside(p, ce1, cen)) {
      if (geo::inside(s, ce1, cen)) {
        if ( ! geo::output(p, no, vo)) return false;
      } else {
        geo::intersect(s, p, ce1, cen, intersection);
        if ( ! geo::output(intersection, no, vo)) return false;
        if ( ! geo::output(p, no, vo)) return false;
      }
    } else if (geo::inside(s, ce1, cen)) {
      geo::intersect(s, p, ce1, cen, intersection);
      if ( ! geo::output(intersection, no, vo)) return false;
    }
    s = p;
  }
  return true;
}

// Efficient user routine that uses the mesh data structure.
//todo An optimization would be to have 2 clip_against_edge routines. One would
// handle the special case of the first vertex list being in (p,e) format.
template <typename geo, typename MeshT, typename CV3s, typename V3s>
KOKKOS_INLINE_FUNCTION
bool clip_against_poly (
  // Clip mesh. m.e(:,cp_e) is the element, and m.en(:,cp_e) is the
  // corresponding list of normal indices.
  const MeshT& m, const Int cp_e,
  // A list of vertices describing the polygon to clip. The vertices must be in
  // a convention-determined order, such as CCW. vi(:,1:ni-1) are valid entries.
  const CV3s& vi, const Int ni,
  // On output, vo(:,0:no-1) are vertices of the clipped polygon. no is 0 if
  // there is no intersection.
  V3s& vo, Int& no,
  // Workspace. Both vo and wrk must be large enough to hold all generated
  // vertices. If they are not, false is returned.
  V3s& wrk)
{
  Int nos[] = { 0, 0 };
  V3s* vs[] = { &vo, &wrk };

  const auto e = slice(m.e, cp_e);
  const auto en = slice(m.en, cp_e);

  auto nv = szslice(m.e); // Number of vertices in clip polygon.
  while (e[nv-1] == -1) --nv;

  no = 0;
  if (nv % 2 == 0) {
    // Make sure the final vertex output list is in the caller's buffer.
    swap(vs[0], vs[1]);
    swap(nos[0], nos[1]);
  }

  if ( ! clip_against_edge<geo>(vi, ni, *vs[0], nos[0], const_slice(m.p, e[0]),
                                const_slice(m.nml, en[0])))
    return false;
  if ( ! nos[0]) return true;

  for (Int ie = 1, ielim = nv - 1; ; ++ie) {
    if ( ! clip_against_edge<geo>(*vs[0], nos[0], *vs[1], nos[1],
                                  const_slice(m.p, e[ie]),
                                  const_slice(m.nml, en[ie])))
      return false;
    if ( ! nos[1]) return true;
    if (ie == ielim) break;
    swap(vs[0], vs[1]);
    swap(nos[0], nos[1]);
  }

  no = nos[1];
  return true;
}

// Not used for real stuff; just a convenient version for testing. In this
// version, clip_poly is a list of clip polygon vertices. This is instead of the
// mesh data structure.
template <typename geo, typename CV3s_CP, typename CV3s_CEN, typename CV3s_VI,
          typename V3s>
KOKKOS_INLINE_FUNCTION
bool clip_against_poly (
  // Clip polygon.
  const CV3s_CP& clip_poly,
  // Clip polygon edges' inward-facing normals.
  const CV3s_CEN& clip_edge_normals,
  const CV3s_VI& vi, const Int ni,
  V3s& vo, Int& no,
  V3s& wrk)
{
  Int nos[] = { 0, 0 };
  V3s* vs[] = { &vo, &wrk };

  no = 0;
  if (nslices(clip_poly) % 2 == 0) {
    // Make sure the final vertex output list is in the caller's buffer.
    swap(vs[0], vs[1]);
    swap(nos[0], nos[1]);
  }

  if ( ! clip_against_edge<geo>(vi, ni, *vs[0], nos[0],
                                const_slice(clip_poly, 0),
                                const_slice(clip_edge_normals, 0)))
    return false;
  if ( ! nos[0]) return true;

  for (Int ie = 1, ielim = nslices(clip_poly) - 1; ; ++ie) {
    if ( ! clip_against_edge<geo>(*vs[0], nos[0], *vs[1], nos[1],
                                  const_slice(clip_poly, ie),
                                  const_slice(clip_edge_normals, ie)))
      return false;
    if ( ! nos[1]) return true;
    if (ie == ielim) break;
    swap(vs[0], vs[1]);
    swap(nos[0], nos[1]);
  }

  no = nos[1];
  return true;
}
} // namespace sh

namespace test {
static constexpr Int max_nvert = 20;
static constexpr Int max_hits = 25; // Covers at least a 2-halo.

// In practice, we want to form high-quality normals using information about the
// mesh.
template <typename geo>
void fill_normals (sh::Mesh<ko::HostSpace>& m) {
  // Count number of edges.
  Int ne = 0;
  for (Int ip = 0; ip < nslices(m.e); ++ip)
    for (Int iv = 0; iv < szslice(m.e); ++iv)
      if (m.e(ip,iv) == -1) break; else ++ne;
  // Fill.
  Idxs::HostMirror en("en", nslices(m.e), szslice(m.e));
  ko::deep_copy(en, -1);
  Vec3s::HostMirror nml("nml", ne, 3);
  Int ie = 0;
  for (Int ip = 0; ip < nslices(m.e); ++ip)
    for (Int iv = 0; iv < szslice(m.e); ++iv)
      if (m.e(ip,iv) == -1)
        break;
      else {
        // Somewhat complicated next node index.
        const Int iv_next = (iv+1 == szslice(m.e) ? 0 :
                             (m.e(ip,iv+1) == -1 ? 0 : iv+1));
        geo::edge_normal(slice(m.p, m.e(ip, iv)), slice(m.p, m.e(ip, iv_next)),
                         slice(nml, ie));
        en(ip,iv) = ie;
        ++ie;
      }
  m.en = en;
  m.nml = nml;
}

//todo The current approach is to do redundant clips so that the hits buffer can
// be small and static. Need to think about this.
template <typename geo>
class AreaOTFunctor {
  const TriangleQuadrature quad_;
  const sh::Mesh<>& cm_;
  const ConstVec3s& p_;
  const ConstIdxs& e_;
  const Int k_; // Index into (p,e).
  //todo More efficient method that also works on GPU.
  Int hits_[max_hits];
  Int nh_;
  Real area_;

public:
  KOKKOS_INLINE_FUNCTION
  AreaOTFunctor (const sh::Mesh<>& cm, const ConstVec3s& p, const ConstIdxs& e,
                 const Int& k)
    : cm_(cm), p_(p), e_(e), k_(k), nh_(0), area_(0)
  {}

  KOKKOS_INLINE_FUNCTION void operator() (const Int mesh_elem_idx) {
    // Check whether we've clipped against this polygon before and there was a
    // non-0 intersection.
    for (Int i = 0; i < nh_; ++i)
      if (hits_[i] == mesh_elem_idx)
        return;
    // We have not, so do the intersection.
    Int no = 0;
    {
      // Area of all overlapping regions.
      // In and out vertex lists.
      Real buf[9*max_nvert];
      RawVec3s
        vi(buf, max_nvert, 3),
        vo(buf + 3*max_nvert, max_nvert, 3),
        wrk(buf + 6*max_nvert, max_nvert, 3);
      Int ni;
      ni = 0;
      for (Int i = 0; i < szslice(e_); ++i) {
        if (e_(k_,i) == -1) break;
        copy(slice(vi, i), slice(p_, e_(k_,i)), 3);
        ++ni;
      }
      sh::clip_against_poly<geo>(cm_, mesh_elem_idx, vi, ni, vo, no, wrk);
      if (no) area_ += geo::calc_area(quad_, vo, no);
    }
    if (no) {
      // Non-0 intersection, so record.
      if (nh_ == max_hits) Kokkos::abort("max_hits is too small.");
      hits_[nh_++] = mesh_elem_idx;
    }
  }

  KOKKOS_INLINE_FUNCTION const Real& area () const { return area_; }
};

template <typename geo, typename OctreeT>
class TestAreaOTKernel {
  const sh::Mesh<> cm_;
  const OctreeT ot_;
  mutable ConstVec3s p_;
  mutable ConstIdxs e_;

public:
  typedef Real value_type;

  TestAreaOTKernel (const sh::Mesh<ko::HostSpace>& cm,
                    const ConstVec3s::HostMirror& p_hm,
                    const ConstIdxs::HostMirror& e_hm, const OctreeT& ot)
    : cm_(cm), ot_(ot)
  {
    { Vec3s p; resize_and_copy(p, p_hm); p_ = p; }
    { Idxs e; resize_and_copy(e, e_hm); e_ = e; }
  }
  
  // Clip the k'th polygon in (p,e) against mesh cm.
  KOKKOS_INLINE_FUNCTION void operator() (const Int k, Real& area) const {
    // Clipped element bounding box.
    Real ebb[6];
    OctreeT::calc_bb(p_, slice(e_, k), szslice(e_), ebb);
    // Get list of possible overlaps.
    AreaOTFunctor<geo> f(cm_, p_, e_, k);
    //todo Team threads.
    ot_.apply(ebb, f);
    area += f.area();
  }

  KOKKOS_INLINE_FUNCTION
  void join (volatile value_type& dst, volatile value_type const& src) const
  { dst += src; }
};

template <typename geo> Real test_area_ot (
  const ConstVec3s::HostMirror& cp, const ConstIdxs::HostMirror& ce,
  const ConstVec3s::HostMirror& p, const ConstIdxs::HostMirror& e)
{
  typedef Octree<geo, 10> OctreeT;

  // Clip mesh and edge normal calculation. (In practice, we'd like to use
  // higher-quality edge normals.)
  sh::Mesh<ko::HostSpace> cm; cm.p = cp; cm.e = ce;
  fill_normals<geo>(cm);

  // Oct-tree over the clip mesh.
  OctreeT ot(cp, ce);

  Real area = 0;
  TestAreaOTKernel<geo, OctreeT> f(cm, p, e, ot);
  ko::parallel_reduce(nslices(e), f, area);
  return area;
}
} // namespace test
} // namespace siqk

#endif // INCLUDE_SIQK_INTERSECT_HPP

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
// In the implementation, (a,b) in [0,1] because convex combinations are used
// throughout; but in the user interface, (a,b) in [-1,1] to agree with the
// definition of the reference square.

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
                           const Real a, const Real b, Real q[3]) {
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
                    const Real b, const Real q[3], Real r[3]) {
  calc_ref_to_bilinear(p, e, a, b, r);
  const Real rnorm = std::sqrt(SphereGeometry::norm2(r));
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
void calc_Jacobian (const ConstVec3sT& p, const Quad& e, const Real a,
                    const Real b, Real J[6]) {
  Real r[3];
  for (Int i = 0; i < 3; ++i) {
    Real t1, t2, t3, t4;
    calc_T_row(p, e, i, t1, t2, t3, t4);
    r[  i] = t1*a*b + t2*a + t3*b + t4;
    J[  i] = t1*b + t2;
    J[3+i] = t1*a + t3;
  }
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
  dx[1] = Qtr[1] / n2;
  dx[0] = (Qtr[0] - a*dx[1]) / n1;
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
  impl::calc_ref_to_bilinear(p, e, 0.5*(a+1), 0.5*(b+1), q);
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
  // Tolerance for Newton iteration.
  const Real tol = 1e2*std::numeric_limits<Real>::epsilon())
{
  const Real tol2 = square(tol);
  Real rnorm2 = 1;
  a = b = 0.5;
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
  a = 2*a - 1;
  b = 2*b - 1;
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
struct Info {
  Int sum_nits, max_nits, nfails;
};

class TestSphereToRefKernel {
  const Real a_test[9] = {-0.1, -1e-16, 0, 1e-15, 0.1, 0.7, 1, 1-1e-14, 1.1};
  const Int n_a_test = sizeof(a_test)/sizeof(*a_test);

  const Real tol_;
  mutable ConstVec3s p_;
  mutable ConstIdxs e_;

public:
  typedef Info value_type;

  TestSphereToRefKernel (const ConstVec3s::HostMirror& p_hm,
                         const ConstIdxs::HostMirror& e_hm,
                         const Real tol = 1e1*std::numeric_limits<Real>::epsilon())
    : tol_(tol)
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
    const Real a_t = 2*a_test[i]-1, b_t = 2*a_test[j]-1;
    Real q[3];
    sqr::calc_ref_to_sphere(p_, slice(e_, ei), a_t, b_t, q);
    Real a, b;
    sqr::Info info;
    sqr::calc_sphere_to_ref(p_, slice(e_, ei), q, a, b, &info, 100, tol_);
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
  void join (volatile value_type& dst, volatile value_type const& src) const {
    dst.max_nits = max(dst.max_nits, src.max_nits);
    dst.sum_nits += src.sum_nits;
    dst.nfails += src.nfails;
  }
};

inline Int test_sphere_to_ref (const ConstVec3s::HostMirror& p,
                               const ConstIdxs::HostMirror& e) {
  TestSphereToRefKernel k(p, e);
  Info info;
  ko::parallel_reduce(k.n(), k, info);
  return info.nfails;
}
} // namespace test
} // namespace sqr
} // namespace siqk

#endif // INCLUDE_SIQK_SQR_HPP

// Unit tests.
#include <limits>

//#include "siqk.hpp"
namespace siqk {

#ifdef IR_MAIN
# define INSTANTIATE_PLANE
#endif

static void make_planar_mesh (Vec3s::HostMirror& p, Idxs::HostMirror& e,
                              const Int n) {
  const Real d = std::sqrt(0.5);
  ko::resize(e, n*n, 4);
  ko::resize(p, (n+1)*(n+1), 3);
  for (Int iy = 0; iy < n+1; ++iy)
    for (Int ix = 0; ix < n+1; ++ix) {
      const auto idx = (n+1)*iy + ix;
      p(idx,0) = 2*(static_cast<Real>(ix)/n - 0.5)*d;
      p(idx,1) = 2*(static_cast<Real>(iy)/n - 0.5)*d;
      p(idx,2) = 0;
    }
  for (Int iy = 0; iy < n; ++iy)
    for (Int ix = 0; ix < n; ++ix) {
      const auto idx = n*iy + ix;
      e(idx,0) = (n+1)*iy + ix;
      e(idx,1) = (n+1)*iy + ix+1;
      e(idx,2) = (n+1)*(iy+1) + ix+1;
      e(idx,3) = (n+1)*(iy+1) + ix;
    }
}

// Row-major R.
inline void form_rotation (const Real axis[3], const Real angle, Real r[9]) {
  const Real nrm = std::sqrt(SphereGeometry::norm2(axis));
  const Real& x = axis[0] / nrm, & y = axis[1] / nrm, & z = axis[2] / nrm,
    & th = angle;
  const Real cth = std::cos(th), sth = std::sin(th), omcth = 1 - cth;
  r[0] = cth + x*x*omcth;
  r[3] = y*x*omcth + z*sth;
  r[6] = z*x*omcth - y*sth;
  r[1] = x*y*omcth - z*sth;
  r[4] = cth + y*y*omcth;
  r[7] = z*y*omcth + x*sth;
  r[2] = x*z*omcth + y*sth;
  r[5] = y*z*omcth - x*sth;
  r[8] = cth + z*z*omcth;
}

template <typename V>
static void rotate (const Real R[9], V p) {
  const Real x = p[0], y = p[1], z = p[2];
  p[0] = R[0]*x + R[1]*y + R[2]*z;
  p[1] = R[3]*x + R[4]*y + R[5]*z;
  p[2] = R[6]*x + R[7]*y + R[8]*z;
}

template <typename V>
static void translate (const Real xlate[3], V p) {
  for (Int i = 0; i < 3; ++i) p[i] += xlate[i];
}

static void transform_planar_mesh (const Real R[9], const Real xlate[3],
                                   Vec3s::HostMirror& p) {
  for (Int i = 0; i < nslices(p); ++i) {
    rotate(R, slice(p, i));
    translate(xlate, slice(p, i));
  }
}

// Remove vertices marked unused and adjust numbering.
static void remove_unused_vertices (Vec3s::HostMirror& p, Idxs::HostMirror& e,
                                    const Real unused) {
  // adjust[i] is the number to subtract from i. Hence if e(ei,0) was originally
  // i, it is adjusted to i - adjust[i].
  std::vector<Int> adjust(nslices(p), 0);
  Int rmcnt = 0;
  for (Int i = 0; i < nslices(p); ++i) {
    if (p(i,0) != unused) continue;
    adjust[i] = 1;
    ++rmcnt;
  }
  // Cumsum.
  for (Int i = 1; i < nslices(p); ++i)
    adjust[i] += adjust[i-1];
  // Adjust e.
  for (Int ei = 0; ei < nslices(e); ++ei)
    for (Int k = 0; k < szslice(e); ++k)
      e(ei,k) -= adjust[e(ei,k)];
  // Remove unused from p.
  Vec3s::HostMirror pc("copy", nslices(p), szslice(p));
  ko::deep_copy(pc, p);
  ko::resize(p, nslices(p) - rmcnt, szslice(p));
  for (Int i = 0, j = 0; i < nslices(pc); ++i) {
    if (pc(i,0) == unused) continue;
    for (Int k = 0; k < szslice(pc); ++k) p(j,k) = pc(i,k);
    ++j;
  }
}

// A very simple cube-sphere mesh with nxn elements per face. At least for now
// I'm not bothering with making the elements well proportioned.
void make_cubesphere_mesh (Vec3s::HostMirror& p, Idxs::HostMirror& e,
                           const Int n) {
  // Transformation of the reference mesh make_planar_mesh to make each of the
  // six faces.
  const Real d = std::sqrt(0.5);
  static Real R[6][9] = {{ 1, 0, 0, 0, 0, 0, 0, 1, 0},  // face 0, -y
                         { 0, 0, 0, 1, 0, 0, 0, 1, 0},  //      1, +x
                         {-1, 0, 0, 0, 0, 0, 0, 1, 0},  //      2, +y
                         { 0, 0, 0,-1, 0, 0, 0, 1, 0},  //      3, -x
                         { 1, 0, 0, 0, 1, 0, 0, 0, 0},  //      4, +z
                         {-1, 0, 0, 0, 1, 0, 0, 0, 0}}; //      5, -z
  static Real xlate[6][3] = {{ 0,-d, 0}, { d, 0, 0}, { 0, d, 0},
                             {-d, 0, 0}, { 0, 0, d}, { 0, 0,-d}};
  // Construct 6 uncoupled faces.
  Vec3s::HostMirror ps[6];
  Vec3s::HostMirror& p_ref = ps[0];
  Idxs::HostMirror es[6];
  Idxs::HostMirror& e_ref = es[0];
  make_planar_mesh(p_ref, e_ref, n);
  ko::resize(e, 6*nslices(e_ref), 4);
  ko::resize(p, 6*nslices(p_ref), 3);
  for (Int i = 1; i < 6; ++i) {
    ko::resize(es[i], nslices(e_ref), 4);
    ko::deep_copy(es[i], e_ref);
    ko::resize(ps[i], nslices(p_ref), 3);
    ko::deep_copy(ps[i], p_ref);
    transform_planar_mesh(R[i], xlate[i], ps[i]);
  }
  transform_planar_mesh(R[0], xlate[0], ps[0]);
  // Pack (p,e), accounting for equivalent vertices. For the moment, keep the p
  // slot for an equivalent vertex to make node numbering simpler, but make the
  // value bogus so we know if there's a problem in the numbering.
  const Real unused = -2;
  ko::deep_copy(p, unused);
  Int p_base = 0, e_base = 0;
  { // -y face
    const Vec3s::HostMirror& fp = ps[0];
    Idxs::HostMirror& fe = es[0];
    for (Int j = 0; j < nslices(fp); ++j)
      for (Int k = 0; k < 3; ++k) p(j,k) = fp(j,k);
    for (Int j = 0; j < nslices(fe); ++j)
      for (Int k = 0; k < 4; ++k) e(j,k) = fe(j,k);
    p_base += nslices(p_ref);
    e_base += nslices(e_ref);
  }
  for (Int fi = 1; fi <= 2; ++fi) { // +x, +y faces
    const Vec3s::HostMirror& fp = ps[fi];
    Idxs::HostMirror& fe = es[fi];
    for (Int j = 0; j < nslices(fp); ++j) {
      if (j % (n+1) == 0) continue; // equiv vertex
      for (Int k = 0; k < 3; ++k) p(p_base+j,k) = fp(j,k);
    }
    for (Int j = 0; j < nslices(fe); ++j) {
      for (Int k = 0; k < 4; ++k) fe(j,k) += p_base;
      // Left 2 vertices of left elem on face fi equiv to right 2 vertices of
      // right elem on face fi-1. Write to the face, then copy to e, so that
      // other faces can use these updated data.
      if (j % n == 0) {
        fe(j,0) = es[fi-1](j+n-1,1);
        fe(j,3) = es[fi-1](j+n-1,2);
      }
      for (Int k = 0; k < 4; ++k) e(e_base+j,k) = fe(j,k);
    }
    p_base += nslices(p_ref);
    e_base += nslices(e_ref);
  }
  { // -x face
    const Vec3s::HostMirror& fp = ps[3];
    Idxs::HostMirror& fe = es[3];
    for (Int j = 0; j < nslices(fp); ++j) {
      if (j % (n+1) == 0 || (j+1) % (n+1) == 0) continue;
      for (Int k = 0; k < 3; ++k) p(p_base+j,k) = fp(j,k);
    }
    for (Int j = 0; j < nslices(fe); ++j) {
      for (Int k = 0; k < 4; ++k) fe(j,k) += p_base;
      if (j % n == 0) {
        fe(j,0) = es[2](j+n-1,1);
        fe(j,3) = es[2](j+n-1,2);
      } else if ((j+1) % n == 0) {
        fe(j,1) = es[0]((j+1)-n,0);
        fe(j,2) = es[0]((j+1)-n,3);
      }
      for (Int k = 0; k < 4; ++k) e(e_base+j,k) = fe(j,k);
    }
    p_base += nslices(p_ref);
    e_base += nslices(e_ref);
  }
  { // +z face
    const Vec3s::HostMirror& fp = ps[4];
    Idxs::HostMirror& fe = es[4];
    for (Int j = n+1; j < nslices(fp) - (n+1); ++j) {
      if (j % (n+1) == 0 || (j+1) % (n+1) == 0) continue;
      for (Int k = 0; k < 3; ++k) p(p_base+j,k) = fp(j,k);
    }
    for (Int j = 0; j < nslices(fe); ++j)
      for (Int k = 0; k < 4; ++k) fe(j,k) += p_base;
    for (Int j = 0; j < n; ++j) { // -y
      fe(j,0) = es[0](n*(n-1)+j,3);
      fe(j,1) = es[0](n*(n-1)+j,2);
    }
    for (Int j = 0; j < n; ++j) { // +y
      fe(n*(n-1)+j,2) = es[2](n*n-1-j,3);
      fe(n*(n-1)+j,3) = es[2](n*n-1-j,2);
    }
    for (Int j = 0, i3 = 0; j < nslices(fe); j += n, ++i3) { // -x
      fe(j,0) = es[3](n*n-1-i3,2);
      fe(j,3) = es[3](n*n-1-i3,3);
    }
    for (Int j = n-1, i1 = 0; j < nslices(fe); j += n, ++i1) { // +x
      fe(j,1) = es[1](n*(n-1)+i1,3);
      fe(j,2) = es[1](n*(n-1)+i1,2);
    }
    for (Int j = 0; j < nslices(fe); ++j)
      for (Int k = 0; k < 4; ++k) e(e_base+j,k) = fe(j,k);
    p_base += nslices(p_ref);
    e_base += nslices(e_ref);
  }
  { // -z face
    const Vec3s::HostMirror& fp = ps[5];
    Idxs::HostMirror& fe = es[5];
    for (Int j = n+1; j < nslices(fp) - (n+1); ++j) {
      if (j % (n+1) == 0 || (j+1) % (n+1) == 0) continue;
      for (Int k = 0; k < 3; ++k) p(p_base+j,k) = fp(j,k);
    }
    for (Int j = 0; j < nslices(fe); ++j)
      for (Int k = 0; k < 4; ++k) fe(j,k) += p_base;
    for (Int j = 0; j < n; ++j) { // -y
      fe(j,0) = es[0](n-1-j,1);
      fe(j,1) = es[0](n-1-j,0);
    }
    for (Int j = 0; j < n; ++j) { // +y
      fe(n*(n-1)+j,2) = es[2](j,1);
      fe(n*(n-1)+j,3) = es[2](j,0);
    }
    for (Int j = 0, i3 = 0; j < nslices(fe); j += n, ++i3) { // -x
      fe(j,0) = es[1](i3,0);
      fe(j,3) = es[1](i3,1);
    }
    for (Int j = n-1, i1 = 0; j < nslices(fe); j += n, ++i1) { // +x
      fe(j,1) = es[3](n-1-i1,1);
      fe(j,2) = es[3](n-1-i1,0);
    }
    for (Int j = 0; j < nslices(fe); ++j)
      for (Int k = 0; k < 4; ++k) e(e_base+j,k) = fe(j,k);
  }
  // Now go back and remove the unused vertices and adjust the numbering.
  remove_unused_vertices(p, e, unused);
  // Project to the unit sphere.
  for (Int i = 0; i < nslices(p); ++i)
    SphereGeometry::normalize(slice(p, i));
}

void calc_elem_ctr (const Vec3s::HostMirror& p, const Idxs::HostMirror& e,
                    const Int ei, Real ctr[3]) {
  for (Int j = 0; j < 3; ++j) ctr[j] = 0;
  Int n = 0;
  for (Int i = 0; i < szslice(e); ++i) {
    if (e(ei,i) < 0) break;
    for (Int j = 0; j < 3; ++j) ctr[j] += p(e(ei,i),j);
    ++n;
  }
  for (Int j = 0; j < 3; ++j) ctr[j] /= n;
}

// Return 0 if all elements' subtri normals point outward relative to the
// sphere.
Int check_elem_normal_against_sphere (const Vec3s::HostMirror& p,
                                      const Idxs::HostMirror& e) {
  Int nerr = 0;
  for (Int ei = 0; ei < nslices(e); ++ei) { // for each element
    Real sphere[3]; // ray through elem ctr
    calc_elem_ctr(p, e, ei, sphere);
    for (Int ti = 0; ti < szslice(e) - 2; ++ti) { // for each tri
      if (e(ei,ti+2) < 0) break;
      Real tri_normal[3]; {
        Real v[2][3];
        for (Int j = 0; j < 2; ++j) {
          SphereGeometry::copy(v[j], slice(p, e(ei,ti+j+1)));
          SphereGeometry::axpy(-1, slice(p, e(ei,0)), v[j]);
        }
        SphereGeometry::cross(v[0], v[1], tri_normal);
      }
      if (SphereGeometry::dot(tri_normal, sphere) <= 0)
        ++nerr;
    }
  }
  return nerr;
}

//> Unit test code.

struct Input {
  Int testno;
  Int n;
  Real angle, xlate, ylate;
  bool geo_sphere;

  Input () {}
  Input(Int argc, char** argv);
  void print(std::ostream& os) const;
};

static void project_onto_sphere (Vec3s::HostMirror& p) {
  for (Int ip = 0; ip < nslices(p); ++ip) {
    p(ip,2) = 1;
    SphereGeometry::normalize(slice(p, ip));
  }
}

static void
perturb_mesh (Vec3s::HostMirror& p, const Real angle, const Real xlate,
              const Real ylate) {
  const Real cr = std::cos(angle), sr = std::sin(angle);
  for (Int ip = 0; ip < nslices(p); ++ip) {
    const Real x = p(ip,0), y = p(ip,1);
    p(ip,0) =  cr*x - sr*y + xlate;
    p(ip,1) = -sr*x + cr*y + ylate;
  }  
}

static void
rotate_mesh (Vec3s::HostMirror& p, const Real axis[3], const Real angle) {
  Real R[9];
  form_rotation(axis, angle, R);
  for (Int i = 0; i < nslices(p); ++i)
    rotate(R, slice(p,i));
}

static void fill_quad (const ConstVec3s::HostMirror& p,
                       Vec3s::HostMirror& poly) {
  const Int n = static_cast<int>(std::sqrt(nslices(p) - 1));
  copy(slice(poly, 0), slice(p, 0), 3);
  copy(slice(poly, 1), slice(p, n), 3);
  copy(slice(poly, 2), slice(p, nslices(p) - 1), 3);
  copy(slice(poly, 3), slice(p, nslices(p) - 1 - n), 3);
}

// Area of the outline of (p,e) clipped against the outline of (cp,ce).
template <typename Geo>
static Real calc_true_area (
  const ConstVec3s::HostMirror& cp, const ConstIdxs::HostMirror& ce,
  const ConstVec3s::HostMirror& p, const ConstIdxs::HostMirror& e)
{
  Vec3s::HostMirror clip_poly("clip_poly", 4, 3), poly("poly", 4, 3),
    nml("nml", 4, 3);
  fill_quad(cp, clip_poly);
  fill_quad(p, poly);
  for (Int i = 0; i < 4; ++i)
    Geo::edge_normal(slice(clip_poly, i), slice(clip_poly, (i+1) % 4),
                     slice(nml, i));
  Vec3s::HostMirror vo("vo", test::max_nvert, 3);
  Int no;
  {
    Vec3s::HostMirror wrk("wrk", test::max_nvert, 3);
    sh::clip_against_poly<Geo>(clip_poly, nml, poly, 4, vo, no, wrk);
  }
  return Geo::calc_area_formula(vo, no);
}

template <typename Geo> void finalize_mesh (Vec3s::HostMirror& p) {}
template <> void finalize_mesh<SphereGeometry> (Vec3s::HostMirror& p) {
  project_onto_sphere(p);
}

template <typename Geo>
static Int
test_area (const Int n, const Real angle, const Real xlate, const Real ylate) {
  Vec3s::HostMirror cp;
  Idxs::HostMirror ce;
  make_planar_mesh(cp, ce, n);

  Vec3s::HostMirror p; resize_and_copy(p, cp);
  Idxs::HostMirror e; resize_and_copy(e, ce);
  perturb_mesh(p, angle, xlate, ylate);

  finalize_mesh<Geo>(cp);
  finalize_mesh<Geo>(p);

  const Real ta = calc_true_area<Geo>(cp, ce, p, e);
  const Real a = test::test_area_ot<Geo>(cp, ce, p, e);

  const Real den = std::max(std::abs(ta), std::abs(a));
  const Real re = den == 0 ? 0 : std::abs(a - ta)/den;

  return re < 1e-8 ? 0 : 1;
}

static Int test_cube (const Input& in) {
  Vec3s::HostMirror cp;
  Idxs::HostMirror ce;
  make_cubesphere_mesh(cp, ce, in.n);
  Vec3s::HostMirror p; resize_and_copy(p, cp);
  Idxs::HostMirror e; resize_and_copy(e, ce);
  Int nerr = 0;
  {
    const Int ne = check_elem_normal_against_sphere(cp, ce);
    if (ne) std::cerr << "FAIL: check_elem_normal_against_sphere\n";
    nerr += ne;
  }
  { // Make a copy, perturb it, and compute the area of the sphere from the
    // overlap mesh.
    Real axis[] = {0.1, -0.3, 0.2};
    rotate_mesh(p, axis, in.angle);
    const Real
      a = test::test_area_ot<SphereGeometry>(cp, ce, p, e),
      ta = 4*M_PI,
      re = std::abs(a - ta)/ta;
    nerr += re < 1e-8 ? 0 : 1;
  }
  // Test ref square <-> spherical quad transformations.
  nerr += sqr::test::test_sphere_to_ref(p, e);
  return nerr;
}

template <typename Geo>
Int run (const Input& in) {
  switch (in.testno) {
  case 0:
    return test_area<Geo>(in.n, in.angle, in.xlate, in.ylate);
  case 1:
    return test_cube(in);
  default:
    return 1;
  }
}

inline bool
eq (const std::string& a, const char* const b1, const char* const b2 = 0) {
  return (a == std::string(b1) || (b2 && a == std::string(b2)) ||
          a == std::string("-") + std::string(b1));
}

Input::Input (Int argc, char** argv)
  : testno(0), n(25), angle(M_PI*1e-1), xlate(1e-1), ylate(1e-1),
    geo_sphere(true)
{
  for (Int i = 1; i < argc; ++i) {
    const std::string& token = argv[i];
    if (eq(token, "--testno")) testno = atoi(argv[++i]);
    else if (eq(token, "-n")) n = atoi(argv[++i]);
    else if (eq(token, "--plane")) geo_sphere = false;
    else if (eq(token, "--xlate")) xlate = atof(argv[++i]);
    else if (eq(token, "--ylate")) ylate = atof(argv[++i]);
    else if (eq(token, "--angle")) angle = atof(argv[++i]);
  }

  print(std::cout);
}

void Input::print (std::ostream& os) const {
  os << "testno " << testno << "\n"
     << "n (-n): " << n << "\n"
     << "planar geometry (--plane): " << ! geo_sphere << "\n"
     << "angle (--angle): " << angle << "\n"
     << "xlate (--xlate): " << xlate << "\n"
     << "ylate (--ylate): " << ylate << "\n";
}

int unittest (const Input& in) {
  Int nerr = 0;
  if (in.geo_sphere)
    nerr += run<SphereGeometry>(in);
  else {
#ifdef INSTANTIATE_PLANE
    nerr += run<PlaneGeometry>(in);
#else
    Kokkos::abort("PlaneGeometry not instantiated.");
#endif
  }
  return nerr;
}

int unittest () {
  static const Int ns[] = {4, 11, 21};
  static const Real values[] = {4.2e-17, 4.2e-11, 4.2e-3, 4.2};
  Input in;
  int nerr = 0;
  for (size_t ni = 0; ni < sizeof(ns)/sizeof(*ns); ++ni) {
    in.n = ns[ni];
    for (size_t ai = 0; ai < sizeof(values)/sizeof(*values); ++ai) {
      in.angle = values[ai];
      in.testno = 1;
      nerr += unittest(in);
      in.testno = 0;
      for (size_t xi = 0; xi < sizeof(values)/sizeof(*values); ++xi) {
        in.xlate = values[xi];
        in.ylate = values[xi];
        nerr += unittest(in);
      }
    }
  }
  return nerr;
}

} // namespace siqk

//> end SIQK

#include <sys/time.h>
#include <mpi.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <list>
#include <vector>
#include <memory>
#include <limits>
#include <algorithm>

#include <Kokkos_Core.hpp>

#define ir_assert(condition) do {                                       \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
#define ir_throw_if(condition, message) do {                            \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\nled to the exception\n" << message << "\n"; \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)

extern "C" void
dpotrf_(const char* uplo, const int* n, double* a, const int* lda, int* info);
extern "C" void
dpotrs_(const char* uplo, const int* n, const int* nrhs, const double* a,
        const int* lda, double* b, const int* ldb, int* info);

namespace ir {
using siqk::Int;
using siqk::Real;
typedef Int Size;

namespace ko = Kokkos;

// A 2D array A can be thought of as having nslices(A) rows and szslice(A)
// columns. A slice can be obtained by
//     auto ak = slice(A, k);
// We use this format for arrays of vertices and adjacency arrays, for
// example. In most or all cases, the intention is to parallelize over slices,
// so a Kokkos operator() will do work on a particular slice.
using siqk::nslices;
using siqk::szslice;
using siqk::slice;

template<typename V> inline Int len (const V& v)
{ return static_cast<Int>(v.dimension_0()); }

template<typename T> inline Int len (const std::vector<T>& v)
{ return static_cast<Int>(v.size()); }

struct Basis {
  Int np, monotone_type;
  Basis () : np(0), monotone_type(-1) {} // invalid values
  Basis (const Int inp, const Int imt) : np(inp), monotone_type(imt) {}
};

// These basis functions follow TempestRemap's monotone_type 0 and 1 basis
// functions.
class GLL {
  const Real oo3 = 1.0/3.0;
  const Real to3 = 2.0/3.0;
  const Real oo6 = 1.0/6.0;
  const Real oo8 = 1.0/8.0;
  const Real sqrt5 = std::sqrt(5.0);
  const Real oosqrt5 = 1.0/sqrt5;
  const Real np2_coord[2] = {-1.0, 1.0};
  const Real np2_wt[2]    = {1.0, 1.0};
  const Real np3_coord[3] = {-1.0, 0.0, 1.0};
  const Real np3_wt[3]    = {oo3, 2.0 - to3, oo3};
  const Real np4_coord[4] = {-1.0, -oosqrt5, oosqrt5, 1.0};
  const Real np4_wt[4]    = {oo6, 1.0 - oo6, 1.0 - oo6, oo6};

public:
  enum { np_max = 4 };

  KOKKOS_INLINE_FUNCTION GLL () {}

  KOKKOS_INLINE_FUNCTION
  void get_coef (const Int np, const Real*& coord, const Real*& wt) const {
    switch (np) {
    case 2:
      coord = np2_coord;
      wt = np2_wt;
      break;
    case 3:
      coord = np3_coord;
      wt = np3_wt;
      break;
    case 4:
      coord = np4_coord;
      wt = np4_wt;
      break;
    default:
      coord = nullptr;
      wt = nullptr;
      siqk::error("GLL::get_coef: order not supported.");
    }
  }
  KOKKOS_INLINE_FUNCTION
  void get_coef (const Basis& b, const Real*& coord, const Real*& wt) const {
    get_coef(b.np, coord, wt);
  }

  // x in [-1, 1].
  KOKKOS_INLINE_FUNCTION
  void eval (const Basis& b, const Real& x, Real* const ge) const {
    if (b.monotone_type < 0 || b.monotone_type > 1)
      siqk::error("GLL::eval: monotone_type not supported.");
    switch (b.np) {
    case 2: {
      ge[0] = 0.5*(1.0 - x);
      ge[1] = 0.5*(1.0 + x);
    } break;
    case 3: {
      const Real x2 = x*x;
      if (b.monotone_type == 0) {
        ge[0] = 0.5*(x2 - x);
        ge[1] = 1.0 - x2;
        ge[2] = 0.5*(x2 + x);
      } else {
        if (x < 0) {
          ge[0] = x2;
          ge[1] = 1 - x2;
          ge[2] = 0;
        } else {
          ge[0] = 0;
          ge[1] = 1 - x2;
          ge[2] = x2;
        }
      }
    } break;
    case 4: {
      const Real x2 = x*x;
      if (b.monotone_type == 0) {
        ge[0] = (1.0 - x)*(5.0*x2 - 1.0)*oo8;
        ge[1] = -sqrt5*oo8*(sqrt5 - 5.0*x)*(x2 - 1.0);
        ge[2] = -sqrt5*oo8*(sqrt5 + 5.0*x)*(x2 - 1.0);
        ge[3] = (1.0 + x)*(5.0*x2 - 1.0)*oo8;
      } else {
        const Real a0 = (1 + sqrt5)/16;
        const Real b0 = (5 + sqrt5)/16;
        const Real c0 = (-5 - 5*sqrt5)/16;
        const Real d0 = (-25 - 5*sqrt5)/16;
        const Real b1 = -0.75*sqrt5;
        const Real d1 = 1.25*sqrt5;
        if (x < -oosqrt5) {
          ge[0] = a0 + x*(b0 + x*(c0 + x*d0));
          ge[1] = 1 - ge[0];
          ge[2] = 0;
          ge[3] = 0;
        } else if (x < oosqrt5) {
          ge[0] = 0;
          ge[1] = 0.5 + x*(b1 + x2*d1);
          ge[2] = 1 - ge[1];
          ge[3] = 0;
        } else {
          ge[0] = 0;
          ge[1] = 0;
          ge[3] = a0 - x*(b0 - x*(c0 - x*d0));
          ge[2] = 1 - ge[3];
        }
      }
    } break;
    default:
      for (int i = 0; i < b.np; ++i) ge[i] = 0;
      siqk::error("GLL::eval: order not supported.");
    }
  }

  KOKKOS_INLINE_FUNCTION
  void eval_derivative (const Basis& b, const Real& x, Real* const ge) const {
    if (b.monotone_type != 0)
      siqk::error("GLL::eval: monotone_type not supported.");
    switch (b.np) {
    case 2: {
      ge[0] = -0.5;
      ge[1] =  0.5;
    } break;
    case 3: {
      const Real x2p = 2*x;
      ge[0] =  0.5*(x2p - 1);
      ge[1] = -x2p;
      ge[2] =  0.5*(x2p + 1);
    } break;
    case 4: {
      ge[0] =  oo8*((10 - 15*x)*x + 1);
      ge[1] = -sqrt5*oo8*((2*sqrt5 - 15*x)*x + 5);
      ge[2] = -sqrt5*oo8*((2*sqrt5 + 15*x)*x - 5);
      ge[3] =  oo8*((10 + 15*x)*x - 1);
    } break;
    default:
      for (int i = 0; i < b.np; ++i) ge[i] = 0;
      siqk::error("GLL::eval: order not supported.");
    }
  }
};

int dpotrf (char uplo, int n, double* a, int lda) {
  int info;
  dpotrf_(&uplo, &n, a, &lda, &info);
  return info;
}

int dpotrs (char uplo, int n, int nrhs, const double* a, int lda, double* bx,
            int ldb) {
  int info;
  dpotrs_(&uplo, &n, &nrhs, const_cast<double*>(a), &lda, bx, &ldb, &info);
  return info;
}

struct Remapper {
  typedef std::shared_ptr<Remapper> Ptr;
  typedef siqk::sh::Mesh<ko::HostSpace> Mesh;

  Remapper (const Int np, const Int nelem)
    : np_(np), np2_(np*np), np4_(np2_*np2_), tq_order_(12)
  {
    mass_mix_.resize(np4_);
    mass_tgt_.resize(np4_);
    local_mesh_.resize(nelem);
  }

  Int np  () const { return np_ ; }
  Int np2 () const { return np2_; }
  Int np4 () const { return np4_; }
  Int tq_order () const { return tq_order_; }

  template <typename Array3D>
  void init_local_mesh_if_needed (const Int ie, const Array3D& corners) {
    ir_assert(ie < static_cast<Int>(local_mesh_.size()));
    if (local_mesh_[ie].p.dimension_0() != 0) return;
    auto& m = local_mesh_[ie];
    const Int
      nd = 3,
      nvert = corners.dimension_1(),
      ncell = corners.dimension_2(),
      N = nvert*ncell;
    m.p = typename Mesh::RealArray("p", N, nd);
    m.e = typename Mesh::IntArray("e", ncell, nvert);
    for (Int ci = 0, k = 0; ci < ncell; ++ci)
      for (Int vi = 0; vi < nvert; ++vi, ++k) {
        for (int j = 0; j < nd; ++j)
          m.p(k,j) = corners(j,vi,ci);
        m.e(ci,vi) = k;
      }
    siqk::test::fill_normals<siqk::SphereGeometry>(m);
  }

  const Mesh& local_mesh (const Int ie) const {
    ir_assert(ie < static_cast<Int>(local_mesh_.size()));
    return local_mesh_[ie];
  };

  std::vector<Real>& rhs_buffer (const Int qsize) {
    rhs_.resize(np2_*qsize);
    return rhs_;
  }

  std::vector<Real>& mass_mix_buffer () { return mass_mix_; }
  std::vector<Real>& mass_tgt_buffer () { return mass_tgt_; }

private:

  const Int np_, np2_, np4_, tq_order_;
  std::vector<Real> mass_tgt_, mass_mix_, rhs_;
  std::vector<Mesh> local_mesh_;
};

static Int test_gll () {
  Int nerr = 0;
  const Real tol = 1e2*std::numeric_limits<Real>::epsilon();
  GLL gll;
  const Real* x, * wt;
  for (Int np = 2; np <= 4; ++np) {
    for (Int monotone_type = 0; monotone_type <= 1; ++monotone_type) {
      const Basis b(np, monotone_type);
      gll.get_coef(b, x, wt);
      Real sum = 0;
      for (Int i = 0; i < b.np; ++i)
        sum += wt[i];
      if (std::abs(2 - sum) > tol) {
        std::cerr << "test_gll " << np << ", " << monotone_type << ": 2 - sum = "
                  << 2 - sum << "\n";
        ++nerr;
      }
      for (Int j = 0; j < b.np; ++j) {
        Real gj[GLL::np_max];
        gll.eval(b, x[j], gj);
        for (Int i = 0; i < b.np; ++i) {
          if (j == i) continue;
          if (std::abs(gj[i]) > tol) {
            std::cerr << "test_gll " << np << ", " << monotone_type << ": gj["
                      << i << "] = " << gj[i] << "\n";
            ++nerr;
          }
        }
      }
    }
  }
  for (Int np = 2; np <= 4; ++np) {
    const Basis b(np, 0);
    Real a[] = {-0.9, -0.7, -0.3, 0.1, 0.2, 0.4, 0.6, 0.8};
    const Real delta = std::sqrt(std::numeric_limits<Real>::epsilon());
    for (size_t ia = 0; ia < sizeof(a)/sizeof(Real); ++ia) {
      Real gj[GLL::np_max], gjp[GLL::np_max], gjm[GLL::np_max];
      gll.eval_derivative(b, a[ia], gj);
      gll.eval(b, a[ia] + delta, gjp);
      gll.eval(b, a[ia] - delta, gjm);
      for (Int i = 0; i < b.np; ++i) {
        const Real fd = (gjp[i] - gjm[i])/(2*delta);
        if (std::abs(fd - gj[i]) >= delta*std::abs(gjp[i]))
          ++nerr;
      }
    }
  }
  return nerr;
}

int unittest () {
  int nerr = 0;
  nerr += test_gll();
  return nerr;
}
} // namespace ir

#ifdef IR_MAIN
int main (int argc, char** argv) {
  int ret = 0;
  Kokkos::initialize(argc, argv);
  try {
    int nerr = siqk::unittest();
    nerr += ir::unittest();
    std::cerr << (nerr ? "FAIL" : "PASS") << "ED\n";
  } catch (const std::exception& e) {
    std::cerr << e.what();
  }
  Kokkos::finalize_all();
  return ret;
}
#endif

namespace homme {
typedef ir::Int Int;
typedef ir::Real Real;

namespace ko = Kokkos;

// Fortran array wrappers.
template <typename T> using FA2 =
  Kokkos::View<T**,    Kokkos::LayoutLeft, Kokkos::HostSpace>;
template <typename T> using FA3 =
  Kokkos::View<T***,  Kokkos::LayoutLeft, Kokkos::HostSpace>;
template <typename T> using FA4 =
  Kokkos::View<T****,  Kokkos::LayoutLeft, Kokkos::HostSpace>;
template <typename T> using FA5 =
  Kokkos::View<T*****, Kokkos::LayoutLeft, Kokkos::HostSpace>;

struct Cartesian3D { Real x, y, z; };

void study (const Int elem_global_id, const Cartesian3D* corners,
            const Int* desc_global_ids, const Cartesian3D* desc_neigh_corners,
            const Int n) {
  int pid;
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  if (pid != 0) return;
  std::stringstream ss;
  const auto out = [&] (const Cartesian3D& c) {
    ss << " (, " << c.x << " " << c.y << " " << c.z << ")";
  };
  ss << "{:gid " << elem_global_id << "\n";
  ss << " :gids (,";
  for (Int i = 0; i < n; ++i) ss << " " << desc_global_ids[i];
  ss << ")\n :me (,";
  for (Int i = 0; i < 4; ++i)
    out(corners[i]);
  ss << ")";
  for (Int j = 0, k = 0; j < n; ++j) {
    ss << "\n :n" << j << " (,";
    for (Int i = 0; i < 4; ++i, ++k)
      out(desc_neigh_corners[k]);
    ss << ")";
  }
  ss << "}\n";
  std::cout << ss.str();
}

static ir::Remapper::Ptr g_remapper;

void ir_init (const Int np, const Int nelem) {
  g_remapper = std::make_shared<ir::Remapper>(np, nelem);
}

// On HSW, this demonstrably is not needed; speedup in ir_project_np4 is a
// result of loop details and local memory. However, I'll keep this in case it's
// of use on KNL.
#define IR_ALIGN(decl) Real decl __attribute__((aligned(64)))

void gll_np4_eval (const Real x, Real y[4]) {
  static constexpr Real oo8 = 1.0/8.0;
  static const Real sqrt5 = std::sqrt(5.0);
  const Real x2 = x*x;
  y[0] = (1.0 - x)*(5.0*x2 - 1.0)*oo8;
  y[1] = -sqrt5*oo8*(sqrt5 - 5.0*x)*(x2 - 1.0);
  y[2] = -sqrt5*oo8*(sqrt5 + 5.0*x)*(x2 - 1.0);
  y[3] = (1.0 + x)*(5.0*x2 - 1.0)*oo8;
}

void ir_project_np4 (
  const Int lev, const Int nets, const Int ie,
  const Int nnc, const Int np, const Int nlev, const Int qsize,
  const Real* dep_points_r, const Real* neigh_corners_r,
  const Real* Qj_src_r, const Real* jac_tgt_r, const Real* rho_tgt_r, const Int tl_np1,
  Real* q_r, Real* q_min_r, Real* q_max_r)
{
  ir_assert(g_remapper);
  using ir::slice;

  static constexpr Int max_num_vertices = 4;
  // Convex quad-quad intersection yields <= 8 intersection points.
  static constexpr Int max_num_intersections = 8;
  static constexpr Int my_np = 4, my_np2 = my_np*my_np,
    my_np4 = my_np2*my_np2;
  static constexpr Int max_qsize = 40;
  const Int ie0 = ie - nets;

  FA2<const Real>
    jac_tgt(jac_tgt_r, np, np);
  FA3<const Real>
    dep_points(reinterpret_cast<const Real*>(dep_points_r), 3, np, np),
    neigh_corners(reinterpret_cast<const Real*>(neigh_corners_r), 3, 4, nnc);
  FA4<const Real>
    Qj_src(Qj_src_r, np, np, qsize+1, nnc),
    rho_tgt(rho_tgt_r, np, np, nlev, tl_np1+1);
  FA4<Real>
    q(q_r, np, np, nlev, qsize);
  FA5<Real>
    q_min(q_min_r, np, np, nlev, qsize, ie0+1),
    q_max(q_max_r, np, np, nlev, qsize, ie0+1);

  IR_ALIGN(M_tgt[my_np4]);
  IR_ALIGN(M_mix[my_np4]);
  IR_ALIGN(rhs_r[my_np2*max_qsize]);
  IR_ALIGN(vi_buf[3*max_num_vertices]);
  IR_ALIGN(vo_buf[3*max_num_intersections]);
  IR_ALIGN(wrk_buf[4*max_num_intersections]);
  IR_ALIGN(sgj[ir::GLL::np_max]);
  IR_ALIGN(sgi[ir::GLL::np_max]);
  IR_ALIGN(tgj[ir::GLL::np_max]);
  IR_ALIGN(tgi[ir::GLL::np_max]);
  IR_ALIGN(svo_coord[2]);
  IR_ALIGN(tvo_coord[2]);

  //todo Move this to ir_init.
  g_remapper->init_local_mesh_if_needed(ie, neigh_corners);

  const Int np2 = np*np, np4 = np2*np2;
  const ir::Basis basis(np, 0);
  const ir::GLL gll;
  siqk::TriangleQuadrature tq;
  siqk::RawConstVec3s tq_bary;
  siqk::RawConstArray tq_w;
  tq.get_coef(g_remapper->tq_order(), tq_bary, tq_w);
  const Int nq = ir::len(tq_w);

  for (Int i = 0, n = qsize*np2; i < n; ++i) rhs_r[i] = 0;
  FA3<Real> rhs(rhs_r, np, np, qsize);

  const auto& m = g_remapper->local_mesh(ie);
  bool first = true;
  for (Int i = 0; i < np4; ++i) M_tgt[i] = 0;
  for (Int sci = 0; sci < nnc; ++sci) { // For each source cell:
    // Intersect.
    const siqk::RawVec3s vi(vi_buf, max_num_vertices, 3);
    siqk::RawVec3s vo(vo_buf, max_num_intersections, 3);
    Int no;
    {
      siqk::RawVec3s wrk(wrk_buf, max_num_intersections, 3);
      static const Int crnr[][2] = {{0, 0}, {np-1, 0}, {np-1, np-1}, {0, np-1}};
      for (int i = 0; i < max_num_vertices; ++i)
        siqk::SphereGeometry::copy(
          slice(vi, i),
          ko::subview(dep_points, ko::ALL(), crnr[i][0], crnr[i][1]));
      siqk::sh::clip_against_poly<siqk::SphereGeometry>(
        g_remapper->local_mesh(ie), sci, vi, max_num_vertices, vo, no, wrk);
    }

    // If the intersection is empty, move on.
    if (no == 0) continue;
    ir_assert(no <= max_num_intersections);

    // Update q_min, q_max with this cell's q source data, since it's in the
    // domain of dependence of the target cell.
    for (Int q = 0; q < qsize; ++q) {
      Real q_min_s, q_max_s;
      q_min_s = q_max_s = Qj_src(0,0,q+1,sci) / Qj_src(0,0,0,sci);
      for (Int j = 0; j < np; ++j)
        for (Int i = 0; i < np; ++i) {
          const Real q_s = Qj_src(i,j,q+1,sci) / Qj_src(i,j,0,sci);
          q_min_s = std::min(q_min_s, q_s);
          q_max_s = std::max(q_max_s, q_s);
        }
      if (first) {
        q_min(0,0,lev,q,ie0) = q_min_s;
        q_max(0,0,lev,q,ie0) = q_max_s;
      } else {
        q_min(0,0,lev,q,ie0) = std::min(q_min(0,0,lev,q,ie0), q_min_s);
        q_max(0,0,lev,q,ie0) = std::max(q_max(0,0,lev,q,ie0), q_max_s);
      }
    }
    first = false;

    // Get overlap element vertices in ref elements.
    Real* const tvo = wrk_buf;
    Real* const svo = tvo + 2*no;
    for (Int k = 0; k < no; ++k) {
      static const Real tol = std::numeric_limits<Real>::epsilon();
      static const int max_nits = 10;
      siqk::sqr::calc_sphere_to_ref(m.p, slice(m.e, sci), slice(vo,k),
                                    svo[2*k], svo[2*k+1],
                                    nullptr, max_nits, tol);
      static const Int vie[] = {0, 1, 2, 3};
      siqk::sqr::calc_sphere_to_ref(vi, vie, slice(vo,k),
                                    tvo[2*k], tvo[2*k+1],
                                    nullptr, max_nits, tol);
    }

    for (Int i = 0; i < my_np4; ++i) M_mix[i] = 0;
    for (Int ktri = 1; ktri < no-1; ++ktri) { // triangles in vo
      const Real s_area = siqk::PlaneGeometry::calc_tri_jacobian(
        svo, svo+2*ktri, svo+2*(ktri+1));
      const Real t_area = siqk::PlaneGeometry::calc_tri_jacobian(
        tvo, tvo+2*ktri, tvo+2*(ktri+1));

      for (Int qp = 0; qp < nq; ++qp) { // quad point
        siqk::PlaneGeometry::bary2coord(tvo, tvo+2*ktri, tvo+2*(ktri+1),
                                        slice(tq_bary, qp), tvo_coord);
        gll_np4_eval(tvo_coord[0], tgi);
        gll_np4_eval(tvo_coord[1], tgj);

        {
          Int os = 0;
          const Real d0 = 0.5 * t_area * tq_w[qp];
          for (Int aj = 0, a_basis_idx = 0; aj < np; ++aj) {
            const Real d1 = d0 * tgj[aj];
            for (Int ai = 0; ai < np; ++ai, ++a_basis_idx) {
              Real d2 = d1 * tgi[ai];
              if (a_basis_idx < 4) {
                const Real d30 = d2 * tgj[0];
                M_tgt[os   ] += d30 * tgi[0]; M_tgt[os+ 1] += d30 * tgi[1];
                M_tgt[os+ 2] += d30 * tgi[2]; M_tgt[os+ 3] += d30 * tgi[3];
              }
              if (a_basis_idx < 8) {
                const Real d31 = d2 * tgj[1];
                M_tgt[os+ 4] += d31 * tgi[0]; M_tgt[os+ 5] += d31 * tgi[1];
                M_tgt[os+ 6] += d31 * tgi[2]; M_tgt[os+ 7] += d31 * tgi[3];
              }
              if (a_basis_idx < 12) {
                const Real d32 = d2 * tgj[2];
                M_tgt[os+ 8] += d32 * tgi[0]; M_tgt[os+ 9] += d32 * tgi[1];
                M_tgt[os+10] += d32 * tgi[2]; M_tgt[os+11] += d32 * tgi[3];
              }
              {
                const Real d33 = d2 * tgj[3];
                M_tgt[os+12] += d33 * tgi[0]; M_tgt[os+13] += d33 * tgi[1];
                M_tgt[os+14] += d33 * tgi[2]; M_tgt[os+15] += d33 * tgi[3];
              }
              os += 16;
            }
          }
        }

        siqk::PlaneGeometry::bary2coord(svo, svo+2*ktri, svo+2*(ktri+1),
                                        slice(tq_bary, qp), svo_coord);
        gll_np4_eval(svo_coord[0], sgi);
        gll_np4_eval(svo_coord[1], sgj);

        {
          Int os = 0;
          const Real d0 = 0.5 * s_area * tq_w[qp];
          for (Int tj = 0; tj < my_np; ++tj) {
            const Real d1 = d0 * tgj[tj];
            for (Int ti = 0; ti < my_np; ++ti) {
              Real d2 = d1 * tgi[ti];
              {
                const Real d30 = d2*sgj[0];
                M_mix[os   ] += d30*sgi[0]; M_mix[os+ 1] += d30*sgi[1];
                M_mix[os+ 2] += d30*sgi[2]; M_mix[os+ 3] += d30*sgi[3];
              } {
                const Real d31 = d2*sgj[1];
                M_mix[os+ 4] += d31*sgi[0]; M_mix[os+ 5] += d31*sgi[1];
                M_mix[os+ 6] += d31*sgi[2]; M_mix[os+ 7] += d31*sgi[3];
              } {
                const Real d32 = d2*sgj[2];
                M_mix[os+ 8] += d32*sgi[0]; M_mix[os+ 9] += d32*sgi[1];
                M_mix[os+10] += d32*sgi[2]; M_mix[os+11] += d32*sgi[3];
              } {
                const Real d33 = d2*sgj[3];
                M_mix[os+12] += d33*sgi[0]; M_mix[os+13] += d33*sgi[1];
                M_mix[os+14] += d33*sgi[2]; M_mix[os+15] += d33*sgi[3];
              }
              os += 16;
            }
          }
        }
      }
    }

    // Accumulate rhs += M_mix Qj_src. This is done separately from the
    // construction of M_mix because q is the slow index.
    for (Int q = 0; q < qsize; ++q) {
      Int os = 0;
      for (Int tj = 0; tj < np; ++tj)
        for (Int ti = 0; ti < np; ++ti) {
          const Real* const v = &Qj_src(0,0,q+1,sci);
          rhs(ti,tj,q) +=
            M_mix[os   ] * v[ 0] + M_mix[os+ 1] * v[ 1] + M_mix[os+ 2] * v[ 2] + M_mix[os+ 3] * v[ 3] +
            M_mix[os+ 4] * v[ 4] + M_mix[os+ 5] * v[ 5] + M_mix[os+ 6] * v[ 6] + M_mix[os+ 7] * v[ 7] +
            M_mix[os+ 8] * v[ 8] + M_mix[os+ 9] * v[ 9] + M_mix[os+10] * v[10] + M_mix[os+11] * v[11] +
            M_mix[os+12] * v[12] + M_mix[os+13] * v[13] + M_mix[os+14] * v[14] + M_mix[os+15] * v[15];
          os += 16;
        }
    }
  }

  // Solve M_tgt Qj_tgt = b.
  int info = ir::dpotrf('L', np2, M_tgt, np2);
  ir_assert(info == 0);
  info = ir::dpotrs('L', np2, qsize, M_tgt, np2, rhs_r, np2);

  // Load q with Qj_tgt / (jac_tgt rho_tgt).
  for (Int tq = 0; tq < qsize; ++tq)
    for (Int j = 0; j < np; ++j)
      for (Int i = 0; i < np; ++i)
        q(i,j,lev,tq) = rhs(i,j,tq) / (rho_tgt(i,j,lev,tl_np1) * jac_tgt(i,j));

  // This method assigns the same q_{min,max} to all GLL points.
  for (Int q = 0; q < qsize; ++q) {
    const Real q_min00 = q_min(0,0,lev,q,ie0);
    for (Int j = 0; j < np; ++j)
      for (Int i = 0; i < np; ++i)
        q_min(i,j,lev,q,ie0) = q_min00;
  }
  for (Int q = 0; q < qsize; ++q) {
    const Real q_max00 = q_max(0,0,lev,q,ie0);
    for (Int j = 0; j < np; ++j)
      for (Int i = 0; i < np; ++i)
        q_max(i,j,lev,q,ie0) = q_max00;
  }
}

#undef IR_ALIGN

void ir_project (
  const Int lev, const Int nets, const Int ie,
  const Int nnc, const Int np, const Int nlev, const Int qsize,
  // Geometry.
  const Cartesian3D* dep_points_r,     // dep_points(1:np, 1:np), unit vectors
  const Cartesian3D* neigh_corners_r,  // neigh_corners(1:4, 1:nnc), unit vectors
  // Fields.
  const Real* Qj_src_r,                //   Qj_src(1:np, 1:np, 2:qsize+1, 1:nnc) and
                                       // rhoj_src(1:np, 1:np, 1        , 1:nnc)
  const Real* jac_tgt_r,               // jac_tgt(1:np, 1:np)
  // rho_tgt(1:np, 1:np, lev, tl_np1+1)
  const Real* rho_tgt_r, const Int tl_np1,
  // Output target tracer mixing ratio and bounds derived from the domain of
  // dependence.
  Real* q_r,                      // q(1:np, 1:np, lev, 1:qsize)
  Real* q_min_r, Real* q_max_r)   // q_{min,max}(1:np, 1:np, lev, 1:qsize, ie-nets+1)
{
  ir_assert(g_remapper);
  if (np == 4) {
    ir_project_np4(lev, nets, ie, nnc, np, nlev, qsize, 
                   reinterpret_cast<const Real*>(dep_points_r),
                   reinterpret_cast<const Real*>(neigh_corners_r),
                   Qj_src_r, jac_tgt_r, rho_tgt_r, tl_np1,
                   q_r, q_min_r, q_max_r);
    return;
  }

  using ir::slice;

  static constexpr Int max_num_vertices = 4;
  // Convex quad-quad intersection yields <= 8 intersection points.
  static constexpr Int max_num_intersections = 8;
  const Int ie0 = ie - nets;

  // In rho_tgt, q, q_{min,max}, the fact that Fortran uses layout left means we
  // don't need to know the rightmost dim.
  FA2<const Real>
    jac_tgt(jac_tgt_r, np, np);
  FA3<const Real>
    dep_points(reinterpret_cast<const Real*>(dep_points_r), 3, np, np),
    neigh_corners(reinterpret_cast<const Real*>(neigh_corners_r), 3, 4, nnc);
  FA4<const Real>
    Qj_src(Qj_src_r, np, np, qsize+1, nnc),
    rho_tgt(rho_tgt_r, np, np, nlev, tl_np1+1);
  FA4<Real>
    q(q_r, np, np, nlev, qsize);
  FA5<Real>
    q_min(q_min_r, np, np, nlev, qsize, ie0+1),
    q_max(q_max_r, np, np, nlev, qsize, ie0+1);

  //todo Move this to ir_init.
  g_remapper->init_local_mesh_if_needed(ie, neigh_corners);

  const Int np2 = np*np, np4 = np2*np2;
  const ir::Basis basis(np, 0);
  const ir::GLL gll;
  siqk::TriangleQuadrature tq;
  siqk::RawConstVec3s tq_bary;
  siqk::RawConstArray tq_w;
  tq.get_coef(g_remapper->tq_order(), tq_bary, tq_w);
  const Int nq = ir::len(tq_w);

  auto& rhs_r = g_remapper->rhs_buffer(qsize);
  for (Int i = 0, n = qsize*np2; i < n; ++i) rhs_r[i] = 0;
  FA3<Real> rhs(rhs_r.data(), np, np, qsize);
  auto& M_mix = g_remapper->mass_mix_buffer();
  auto& M_tgt = g_remapper->mass_tgt_buffer();
  for (Int i = 0; i < np4; ++i) M_tgt[i] = 0;

  const auto& m = g_remapper->local_mesh(ie);
  bool first = true;
  for (Int sci = 0; sci < nnc; ++sci) { // For each source cell:
    // Intersect.
    Real buf[3*max_num_vertices + 3*max_num_intersections + 4*max_num_intersections];
    const siqk::RawVec3s vi(buf, max_num_vertices, 3);
    siqk::RawVec3s vo(buf + 3*max_num_vertices, max_num_intersections, 3);
    Real* const wrkbuf = vo.data() + 3*max_num_intersections;
    Int no;
    {
      siqk::RawVec3s wrk(wrkbuf, max_num_intersections, 3);
      static const Int crnr[][2] = {{0, 0}, {np-1, 0}, {np-1, np-1}, {0, np-1}};
      for (int i = 0; i < max_num_vertices; ++i)
        siqk::SphereGeometry::copy(
          slice(vi, i),
          ko::subview(dep_points, ko::ALL(), crnr[i][0], crnr[i][1]));
      siqk::sh::clip_against_poly<siqk::SphereGeometry>(
        g_remapper->local_mesh(ie), sci, vi, max_num_vertices, vo, no, wrk);
    }

    // If the intersection is empty, move on.
    if (no == 0) continue;
    ir_assert(no <= max_num_intersections);

    // Update q_min, q_max with this cell's q source data, since it's in the
    // domain of dependence of the target cell.
    for (Int q = 0; q < qsize; ++q) {
      Real q_min_s, q_max_s;
      q_min_s = q_max_s = Qj_src(0,0,q+1,sci) / Qj_src(0,0,0,sci);
      for (Int j = 0; j < np; ++j)
        for (Int i = 0; i < np; ++i) {
          const Real q_s = Qj_src(i,j,q+1,sci) / Qj_src(i,j,0,sci);
          q_min_s = std::min(q_min_s, q_s);
          q_max_s = std::max(q_max_s, q_s);
        }
      if (first) {
        q_min(0,0,lev,q,ie0) = q_min_s;
        q_max(0,0,lev,q,ie0) = q_max_s;
      } else {
        q_min(0,0,lev,q,ie0) = std::min(q_min(0,0,lev,q,ie0), q_min_s);
        q_max(0,0,lev,q,ie0) = std::max(q_max(0,0,lev,q,ie0), q_max_s);
      }
    }
    first = false;

    // Get overlap element vertices in ref elements.
    Real* const tvo = wrkbuf;
    Real* const svo = tvo + 2*no;
    for (Int k = 0; k < no; ++k) {
      static const Real tol = std::numeric_limits<Real>::epsilon();
      static const int max_nits = 10;
      siqk::sqr::calc_sphere_to_ref(m.p, slice(m.e, sci), slice(vo,k),
                                    svo[2*k], svo[2*k+1],
                                    nullptr, max_nits, tol);
      static const Int vie[] = {0, 1, 2, 3};
      siqk::sqr::calc_sphere_to_ref(vi, vie, slice(vo,k),
                                    tvo[2*k], tvo[2*k+1],
                                    nullptr, max_nits, tol);
    }

    for (Int i = 0; i < np4; ++i) M_mix[i] = 0;
    for (Int ktri = 1; ktri < no-1; ++ktri) { // triangles in vo
      const Real s_area = siqk::PlaneGeometry::calc_tri_jacobian(
        svo, svo+2*ktri, svo+2*(ktri+1));
      const Real t_area = siqk::PlaneGeometry::calc_tri_jacobian(
        tvo, tvo+2*ktri, tvo+2*(ktri+1));
      Real
        sgj[ir::GLL::np_max], sgi[ir::GLL::np_max],
        tgj[ir::GLL::np_max], tgi[ir::GLL::np_max],
        svo_coord[2], tvo_coord[2];

      for (Int qp = 0; qp < nq; ++qp) { // quad point
        siqk::PlaneGeometry::bary2coord(tvo, tvo+2*ktri, tvo+2*(ktri+1),
                                        slice(tq_bary, qp), tvo_coord);
        gll.eval(basis, tvo_coord[0], tgi);
        gll.eval(basis, tvo_coord[1], tgj);

        Real d0 = 0.5 * t_area * tq_w[qp];
        for (Int aj = 0, a_basis_idx = 0; aj < np; ++aj) {
          const Real d1 = d0 * tgj[aj];
          for (Int ai = 0; ai < np; ++ai, ++a_basis_idx) {
            Real d2 = d1 * tgi[ai];
            for (Int b_basis_idx = a_basis_idx; b_basis_idx < np2; ++b_basis_idx) {
              const Int bj = b_basis_idx / np;
              const Int bi = b_basis_idx % np;
              M_tgt[np2*a_basis_idx + b_basis_idx] += d2 * tgi[bi] *tgj[bj];
            }
          }
        }

        siqk::PlaneGeometry::bary2coord(svo, svo+2*ktri, svo+2*(ktri+1),
                                        slice(tq_bary, qp), svo_coord);
        gll.eval(basis, svo_coord[0], sgi);
        gll.eval(basis, svo_coord[1], sgj);

        d0 = 0.5 * s_area * tq_w[qp];
        Int os = 0;
        for (Int tj = 0; tj < np; ++tj) {
          const Real d1 = d0 * tgj[tj];
          for (Int ti = 0; ti < np; ++ti) {
            Real d2 = d1 * tgi[ti];
            for (Int sj = 0; sj < np; ++sj) {
              const Real d3 = d2 * sgj[sj];
              for (Int si = 0; si < np; ++si)
                M_mix[os++] += d3 * sgi[si];
            }
          }
        }
      }
    }

    // Accumulate rhs += M_mix Qj_src. This is done separately from the
    // construction of M_mix because q is the slow index.
    for (Int q = 0; q < qsize; ++q) {
      Int os = 0;
      for (Int tj = 0; tj < np; ++tj)
        for (Int ti = 0; ti < np; ++ti) {
          auto& r = rhs(ti,tj,q);
          for (Int sj = 0; sj < np; ++sj)
            for (Int si = 0; si < np; ++si)
              r += M_mix[os++] * Qj_src(si,sj,q+1,sci);
        }
    }
  }

  // Solve M_tgt Qj_tgt = b.
  int info = ir::dpotrf('L', np2, M_tgt.data(), np2);
  ir_assert(info == 0);
  info = ir::dpotrs('L', np2, qsize, M_tgt.data(), np2, rhs.data(), np2);

  // Load q with Qj_tgt / (jac_tgt rho_tgt).
  for (Int tq = 0; tq < qsize; ++tq)
    for (Int j = 0; j < np; ++j)
      for (Int i = 0; i < np; ++i)
        q(i,j,lev,tq) = rhs(i,j,tq) / (rho_tgt(i,j,lev,tl_np1) * jac_tgt(i,j));

  // This method assigns the same q_{min,max} to all GLL points.
  for (Int q = 0; q < qsize; ++q)
    for (Int j = 0; j < np; ++j)
      for (Int i = 0; i < np; ++i)
        q_min(i,j,lev,q,ie0) = q_min(0,0,lev,q,ie0);
  for (Int q = 0; q < qsize; ++q)
    for (Int j = 0; j < np; ++j)
      for (Int i = 0; i < np; ++i)
        q_max(i,j,lev,q,ie0) = q_max(0,0,lev,q,ie0);
}
} // namespace homme

extern "C" void ir_init_ (homme::Int* np, homme::Int* nelem) {
  homme::ir_init(*np, *nelem);
}

// Figure shtuff out.
extern "C" void ir_study_ (
  homme::Int* elem_global_id, homme::Cartesian3D* corners,
  homme::Int** desc_global_ids, homme::Cartesian3D** desc_neigh_corners,
  homme::Int* n)
{
  static int gid_max = -1;
  if (*elem_global_id <= gid_max) return;
  gid_max = *elem_global_id;
  homme::study(*elem_global_id, corners, *desc_global_ids, *desc_neigh_corners, *n + 1);
}

extern "C" void ir_project_ (
  homme::Int* lev, homme::Int* ie, homme::Int* nnc, homme::Int* np, homme::Int* nlev,
  homme::Int* qsize, homme::Int* nets, homme::Int* nete,
  homme::Cartesian3D* dep_points, homme::Cartesian3D** neigh_corners,
  homme::Real* neigh_q, homme::Real* J_t, homme::Real* dp3d, homme::Int* tl_np1,
  homme::Real* q, homme::Real* minq, homme::Real* maxq)
{
  homme::ir_project(*lev - 1, *nets - 1, *ie - 1, *nnc, *np, *nlev, *qsize,
                    dep_points, *neigh_corners, neigh_q,
                    J_t, dp3d, *tl_np1 - 1,
                    q, minq, maxq);
}
