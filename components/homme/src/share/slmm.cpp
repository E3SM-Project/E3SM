// Homme doesn't define this, probably b/c Fortran code doesn't use it that
// much, so define it here. Not doing so can slow the code in this file by
// 10-20%.
#ifndef NDEBUG
# define NDEBUG
#endif
// Uncomment this to look for MPI-related memory leaks.
//#define COMPOSE_DEBUG_MPI

#define BUILD_CISL

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
    _ss_ << "slmm: " <<  m << std::endl;        \
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
#endif

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

#define SIQK_QUADRATURE_TRISYM_ORDER14_COORD \
  {0.3333333333333333148296162562473910, 0.3333333333333333148296162562473910, 0.3333333333333333148296162562473910, \
   0.0099797608064584319986778382371995, 0.0099797608064584319986778382371995, 0.9800404783870830804914930922677740, \
   0.0099797608064584319986778382371995, 0.9800404783870830804914930922677740, 0.0099797608064584319986778382371995, \
   0.9800404783870830804914930922677740, 0.0099797608064584319986778382371995, 0.0099797608064584319986778382371995, \
   0.4799778935211884145495275788562139, 0.4799778935211884145495275788562139, 0.0400442129576231709009448422875721, \
   0.4799778935211884145495275788562139, 0.0400442129576231709009448422875721, 0.4799778935211884145495275788562139, \
   0.0400442129576231709009448422875721, 0.4799778935211884145495275788562139, 0.4799778935211884145495275788562139, \
   0.1538119591769669114444951674158801, 0.1538119591769669114444951674158801, 0.6923760816460662326221608964260668, \
   0.1538119591769669114444951674158801, 0.6923760816460662326221608964260668, 0.1538119591769669114444951674158801, \
   0.6923760816460662326221608964260668, 0.1538119591769669114444951674158801, 0.1538119591769669114444951674158801, \
   0.0740234771169878125185448425327195, 0.0740234771169878125185448425327195, 0.8519530457660243749629103149345610, \
   0.0740234771169878125185448425327195, 0.8519530457660243749629103149345610, 0.0740234771169878125185448425327195, \
   0.8519530457660243749629103149345610, 0.0740234771169878125185448425327195, 0.0740234771169878125185448425327195, \
   0.1303546825033299882967696703417460, 0.1303546825033299882967696703417460, 0.7392906349933400234064606593165081, \
   0.1303546825033299882967696703417460, 0.7392906349933400234064606593165081, 0.1303546825033299882967696703417460, \
   0.7392906349933400234064606593165081, 0.1303546825033299882967696703417460, 0.1303546825033299882967696703417460, \
   0.2306172260266531326422523306973744, 0.2306172260266531326422523306973744, 0.5387655479466937347154953386052512, \
   0.2306172260266531326422523306973744, 0.5387655479466937347154953386052512, 0.2306172260266531326422523306973744, \
   0.5387655479466937347154953386052512, 0.2306172260266531326422523306973744, 0.2306172260266531326422523306973744, \
   0.4223320834191477968211358984262915, 0.4223320834191477968211358984262915, 0.1553358331617044063577282031474169, \
   0.4223320834191477968211358984262915, 0.1553358331617044063577282031474169, 0.4223320834191477968211358984262915, \
   0.1553358331617044063577282031474169, 0.4223320834191477968211358984262915, 0.4223320834191477968211358984262915, \
   0.7862373859346609705767150444444269, 0.1906163600319009110428680742188590, 0.0231462540334381183804168813367141, \
   0.7862373859346609705767150444444269, 0.0231462540334381183804168813367141, 0.1906163600319009110428680742188590, \
   0.1906163600319009110428680742188590, 0.7862373859346609705767150444444269, 0.0231462540334381183804168813367141, \
   0.1906163600319009110428680742188590, 0.0231462540334381183804168813367141, 0.7862373859346609705767150444444269, \
   0.0231462540334381183804168813367141, 0.7862373859346609705767150444444269, 0.1906163600319009110428680742188590, \
   0.0231462540334381183804168813367141, 0.1906163600319009110428680742188590, 0.7862373859346609705767150444444269, \
   0.6305521436606074114905595706659369, 0.3623231377435471300962888108188054, 0.0071247185958454584131516185152577, \
   0.6305521436606074114905595706659369, 0.0071247185958454584131516185152577, 0.3623231377435471300962888108188054, \
   0.3623231377435471300962888108188054, 0.6305521436606074114905595706659369, 0.0071247185958454584131516185152577, \
   0.3623231377435471300962888108188054, 0.0071247185958454584131516185152577, 0.6305521436606074114905595706659369, \
   0.0071247185958454584131516185152577, 0.6305521436606074114905595706659369, 0.3623231377435471300962888108188054, \
   0.0071247185958454584131516185152577, 0.3623231377435471300962888108188054, 0.6305521436606074114905595706659369, \
   0.6265773298563063198329814440512564, 0.2907712058836673940653838599246228, 0.0826514642600262861016346960241208, \
   0.6265773298563063198329814440512564, 0.0826514642600262861016346960241208, 0.2907712058836673940653838599246228, \
   0.2907712058836673940653838599246228, 0.6265773298563063198329814440512564, 0.0826514642600262861016346960241208, \
   0.2907712058836673940653838599246228, 0.0826514642600262861016346960241208, 0.6265773298563063198329814440512564, \
   0.0826514642600262861016346960241208, 0.6265773298563063198329814440512564, 0.2907712058836673940653838599246228, \
   0.0826514642600262861016346960241208, 0.2907712058836673940653838599246228, 0.6265773298563063198329814440512564, \
   0.9142099849296254632236014003865421, 0.0711657108777507679819862573822320, 0.0146243041926237687944123422312259, \
   0.9142099849296254632236014003865421, 0.0146243041926237687944123422312259, 0.0711657108777507679819862573822320, \
   0.0711657108777507679819862573822320, 0.9142099849296254632236014003865421, 0.0146243041926237687944123422312259, \
   0.0711657108777507679819862573822320, 0.0146243041926237687944123422312259, 0.9142099849296254632236014003865421, \
   0.0146243041926237687944123422312259, 0.9142099849296254632236014003865421, 0.0711657108777507679819862573822320, \
   0.0146243041926237687944123422312259, 0.0711657108777507679819862573822320, 0.9142099849296254632236014003865421}
#define SIQK_QUADRATURE_TRISYM_ORDER14_WEIGHT \
  {0.0585962852260285965710906452841300,0.0017351512297252675524200649093132,0.0017351512297252675524200649093132, \
   0.0017351512297252675524200649093132,0.0261637825586145227052536910150593,0.0261637825586145227052536910150593, \
   0.0261637825586145227052536910150593,0.0039197292424018289128118119890587,0.0039197292424018289128118119890587, \
   0.0039197292424018289128118119890587,0.0122473597569408669538670864085361,0.0122473597569408669538670864085361, \
   0.0122473597569408669538670864085361,0.0281996285032579604989955157634540,0.0281996285032579604989955157634540, \
   0.0281996285032579604989955157634540,0.0508870871859594883779287499692146,0.0508870871859594883779287499692146, \
   0.0508870871859594883779287499692146,0.0504534399016036000373830461285252,0.0504534399016036000373830461285252, \
   0.0504534399016036000373830461285252,0.0170636442122334523741056244716674,0.0170636442122334523741056244716674, \
   0.0170636442122334523741056244716674,0.0170636442122334523741056244716674,0.0170636442122334523741056244716674, \
   0.0170636442122334523741056244716674,0.0096834664255066003890615178306689,0.0096834664255066003890615178306689, \
   0.0096834664255066003890615178306689,0.0096834664255066003890615178306689,0.0096834664255066003890615178306689, \
   0.0096834664255066003890615178306689,0.0363857559284850029523994408009457,0.0363857559284850029523994408009457, \
   0.0363857559284850029523994408009457,0.0363857559284850029523994408009457,0.0363857559284850029523994408009457, \
   0.0363857559284850029523994408009457,0.0069646633735184126576256424812073,0.0069646633735184126576256424812073, \
   0.0069646633735184126576256424812073,0.0069646633735184126576256424812073,0.0069646633735184126576256424812073, \
   0.0069646633735184126576256424812073}

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
  const Real trisym_order14_coord_ [138] = SIQK_QUADRATURE_TRISYM_ORDER14_COORD;
  const Real trisym_order14_weight_[ 46] = SIQK_QUADRATURE_TRISYM_ORDER14_WEIGHT;

public:
  KOKKOS_INLINE_FUNCTION TriangleQuadrature () {}

  KOKKOS_INLINE_FUNCTION
  void get_coef (const int order, RawConstVec3s& coord,
                 RawConstArray& weight) const {
    switch (order) {
    case 4:
      coord = RawConstVec3s(trisym_order4_coord_, 6);
      weight = RawConstArray(trisym_order4_weight_, 6);
      break;
    case 6:
      coord = RawConstVec3s(tritay_order6_coord_, 11);
      weight = RawConstArray(tritay_order6_weight_, 11);
      break;
    case 8:
      coord = RawConstVec3s(trisym_order8_coord_, 16);
      weight = RawConstArray(trisym_order8_weight_, 16);
      break;
    case 12:
      coord = RawConstVec3s(tritay_order12_coord_, 32);
      weight = RawConstArray(tritay_order12_weight_, 32);
      break;
    case 14:
      coord = RawConstVec3s(trisym_order14_coord_, 46);
      weight = RawConstArray(trisym_order14_weight_, 46);
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

#ifdef BUILD_CISL
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
#endif // BUILD_CISL

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

#ifdef BUILD_CISL

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
#endif // BUILD_CISL
} // namespace sh

namespace test {
static constexpr Int max_nvert = 20;
static constexpr Int max_hits = 25; // Covers at least a 2-halo.

// Inward-oriented normal. In practice, we want to form high-quality normals
// using information about the cubed-sphere mesh. This is a low-quality
// brute-force calculation.
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
  Vec3s::HostMirror nml("nml", ne);
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

#ifdef BUILD_CISL
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
#endif // BUILD_CISL
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
void calc_Jacobian (const ConstVec3sT& p, const Quad& e, Real a, Real b,
                    Real J[6]) {
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
  // Tolerance for Newton iteration.
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
    const Real a_t = a_test[i], b_t = a_test[j];
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

#ifdef SLMM_MAIN
# define INSTANTIATE_PLANE
#endif
#ifdef BUILD_CISL

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
  Vec3s::HostMirror clip_poly("clip_poly", 4), poly("poly", 4),
    nml("nml", 4);
  fill_quad(cp, clip_poly);
  fill_quad(p, poly);
  for (Int i = 0; i < 4; ++i)
    Geo::edge_normal(slice(clip_poly, i), slice(clip_poly, (i+1) % 4),
                     slice(nml, i));
  Vec3s::HostMirror vo("vo", test::max_nvert);
  Int no;
  {
    Vec3s::HostMirror wrk("wrk", test::max_nvert);
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
#endif // BUILD_CISL

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

#ifndef NDEBUG
# define slmm_assert(condition) do {                                    \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
      throw std::logic_error(_ss_.str());                               \
    }                                                                   \
  } while (0)
#else
# define slmm_assert(condition)
#endif
#define slmm_throw_if(condition, message) do {                          \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\nled to the exception\n" << message << "\n"; \
      throw std::logic_error(_ss_.str());                               \
    }                                                                   \
  } while (0)

extern "C" void
dpotrf_(const char* uplo, const int* n, double* a, const int* lda, int* info);
extern "C" void
dpotrs_(const char* uplo, const int* n, const int* nrhs, const double* a,
        const int* lda, double* b, const int* ldb, int* info);

namespace slmm {
using siqk::Int;
using siqk::Real;
typedef Int Size;

namespace ko = Kokkos;

typedef siqk::sh::Mesh<ko::HostSpace> Mesh;

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
{ return static_cast<Int>(v.extent_int(0)); }

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

void copy_vertices (
  const siqk::ConstVec3s::HostMirror& p, const siqk::ConstIdxs::HostMirror& c2n,
  const Int ci, Real* ps)
{
  const auto cell = slice(c2n, ci);
  for (Int i = 0; i < szslice(c2n); ++i) {
    const auto n = slice(p, cell[i]);
    for (Int k = 0; k < 3; ++k) ps[k] = n[k];
    ps += 3;
  }
}

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
        std::cerr << "test_gll " << np << ", " << monotone_type
                  << ": 2 - sum = " << 2 - sum << "\n";
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

std::string to_string (const Mesh& m) {
  std::stringstream ss;
  ss.precision(17);
  ss << "(mesh nnode " << nslices(m.p) << " nelem " << nslices(m.e);
  for (Int ie = 0; ie < nslices(m.e); ++ie) {
    ss << "\n  (elem " << ie;
    ss << " (";
    for (Int iv = 0; iv < szslice(m.e); ++iv) ss << " " << m.e(ie,iv);
    ss << ")";
    for (Int iv = 0; iv < szslice(m.e); ++iv) {
      ss << "\n     (p";
      for (Int d = 0; d < 3; ++d)
        ss << " " << m.p(m.e(ie,iv),d);
      ss << ")";
      ss << "\n     (nml";
      for (Int d = 0; d < 3; ++d)
        ss << " " << m.nml(m.en(ie,iv),d);
      ss << ")";
    }
    ss << "))";
  }
  ss << ")";
  return ss.str();
}

static const Real sqrt5 = std::sqrt(5.0);
static const Real oosqrt5 = 1.0 / sqrt5;

void gll_np4_eval (const Real x, Real y[4]) {
  static constexpr Real oo8 = 1.0/8.0;
  const Real x2 = x*x;
  y[0] = (1.0 - x)*(5.0*x2 - 1.0)*oo8;
  y[1] = -sqrt5*oo8*(sqrt5 - 5.0*x)*(x2 - 1.0);
  y[2] = -sqrt5*oo8*(sqrt5 + 5.0*x)*(x2 - 1.0);
  y[3] = (1.0 + x)*(5.0*x2 - 1.0)*oo8;
}

// Linear interp in each region.
void gll_np4_subgrid_eval (const Real& x, Real y[4]) {
  if (x > 0) {
    gll_np4_subgrid_eval(-x, y);
    std::swap(y[0], y[3]);
    std::swap(y[1], y[2]);    
    return;
  }
  if (x < -oosqrt5) {
    const Real alpha = (x + 1)/(1 - oosqrt5);
    y[0] = 1 - alpha;
    y[1] = alpha;
    y[2] = 0;
    y[3] = 0;
  } else {
    const Real alpha = (x + oosqrt5)/(2*oosqrt5);
    y[0] = 0;
    y[1] = 1 - alpha;
    y[2] = alpha;
    y[3] = 0;
  }
}

// Quadratic interpolant across nodes 1,2,3 -- i.e., excluding node 0 -- of the
// np=4 reference element.
void outer_eval (const Real& x, Real v[4]) {
  static const Real
    xbar = (2*oosqrt5) / (1 + oosqrt5),
    ooxbar = 1 / xbar,
    ybar = 1 / (xbar - 1);
  const Real xn = (x + oosqrt5) / (1 + oosqrt5);
  v[0] = 0;
  v[1] = 1 + ybar*xn*((1 - ooxbar)*xn + ooxbar - xbar);
  v[2] = ybar*ooxbar*xn*(xn - 1);
  v[3] = ybar*xn*(xbar - xn);
}

// In the middle region, use the standard GLL np=4 interpolant; in the two outer
// regions, use an order-reduced interpolant, but one that stabilizes the
// resulting classical SL method stable for the uniform-flow problem on a
// uniform mesh.
void gll_np4_subgrid_exp_eval (const Real& x, Real y[4]) {
  static constexpr Real
    alpha = 0.5527864045000416708,
    v = 0.427*(1 + alpha),
    x2 = 0.4472135954999579277,
    x3 = 1 - x2,
    det = x2*x3*(x2 - x3),
    y2 = alpha,
    y3 = v,
    c1 = (x3*y2 - x2*y3)/det,
    c2 = (-x3*x3*y2 + x2*x2*y3)/det;
  if (x < -oosqrt5 || x > oosqrt5) {
    if (x < -oosqrt5) {
      outer_eval(-x, y);
      std::swap(y[0], y[3]);
      std::swap(y[1], y[2]);
    } else
      outer_eval(x, y);
    Real y4[4];
    gll_np4_eval(x, y4);
    const Real x0 = 1 - std::abs(x);
    const Real a = (c1*x0 + c2)*x0;
    for (int i = 0; i < 4; ++i)
      y[i] = a*y[i] + (1 - a)*y4[i];
  } else
    gll_np4_eval(x, y);
}

// Is v inside, including on, the quad ic?
inline bool is_inside (const siqk::sh::Mesh<siqk::ko::HostSpace>& m,
                       const Real* v, const Real& atol, const Int& ic) {
  using slmm::slice;
  const auto cell = slice(m.e, ic);
  const auto celln = slice(m.en, ic);
  const Int ne = 4;
  assert(szslice(m.e) == ne);
  assert(szslice(m.en) == ne);
  bool inside = true;
  for (Int ie = 0; ie < ne; ++ie)
    if (siqk::SphereGeometry::dot_c_amb(slice(m.nml, celln[ie]),
                                        v, slice(m.p, cell[ie]))
        < -atol) {
      inside = false;
      break;
    }
  return inside;
}

int get_src_cell (const siqk::sh::Mesh<siqk::ko::HostSpace>& m, // Local mesh.
                  const Real* v, // 3D Cartesian point.
                  const Int my_ic = -1) { // Target cell in the local mesh.
  using slmm::len;
  const Int nc = len(m.e);
  Real atol = 0;
  for (Int trial = 0; trial < 3; ++trial) {
    if (trial > 0) {
      if (trial == 1) {
        using slmm::slice;
        // If !inside in the first sweep, pad each cell. Recall we're operating
        // on the unit sphere, so we don't have to worry about a radius in the
        // following.
        //   Get a representative edge length.
        const int ic = my_ic == -1 ? 0 : my_ic;
        const auto cell = slice(m.e, ic);
        Real d[3];
        siqk::SphereGeometry::axpbyz(1, slice(m.p, cell[1]),
                                     -1, slice(m.p, cell[0]),
                                     d);
        const Real L = std::sqrt(siqk::SphereGeometry::norm2(d));
        // We can expect to lose approx. -log10(L) digits due to cancellation in
        // the formation of the edge normal and in dot_c_amb. Multiply by 100
        // for a little extra padding.
        atol = 1e2 * std::numeric_limits<Real>::epsilon() / L;
      } else {
        // Ok, we really didn't do that very well. We're still failing to find
        // the element. Ramp up the atol even more.
        atol = std::max(1e3*atol,
                        std::sqrt(std::numeric_limits<Real>::epsilon()));
      }
    }
    if (my_ic != -1 && is_inside(m, v, atol, my_ic)) return my_ic;
    for (Int ic = 0; ic < nc; ++ic) {
      if (ic == my_ic) continue;
      if (is_inside(m, v, atol, ic)) return ic;
    }
  }
  return -1;
}

namespace nearest_point {
/* Get external segments in preproc step.
   Get approximate nearest point in each segment.
     Project onto plane of segment and normalize.
     If outside of arc, choose nearest endpoint.
     Approx distance as cartesian distance.
   Of the approx dists, take the point associated with smallest.
 */

template <typename ES = ko::DefaultExecutionSpace>
struct MeshNearestPointData {
  typedef ko::View<Int*,ES> Ints;
  // mesh.p(perimp(k),:) is the k'th vertex in the mesh's perimeter. Similarly,
  // mesh.nml(perimnml(k),:) is the k'th edge's normal.
  Ints perimp, perimnml;
};

template <typename V3a, typename V3b>
inline Real calc_dist (const V3a& p0, const V3b& p1) {
  Real len = 0;
  for (Int d = 0; d < 3; ++d) len += siqk::square(p0[d] - p1[d]);
  return std::sqrt(len);    
}

inline Real calc_dist (const Mesh& m, const Int& i0, const Int& i1) {
  return calc_dist(slice(m.p, i0), slice(m.p, i1));
}

void find_external_edges (const Mesh& m, std::vector<Int>& external_edges) {
  const Int ncell = nslices(m.e);
  const Int nedge = 4*ncell;
  // Get the minimum edge length to normalize distance.
  Real min_edge_len = 10;
  for (Int ic = 0; ic < ncell; ++ic) {
    const auto cell = slice(m.e, ic);
    for (Int ie = 0; ie < 4; ++ie)
      min_edge_len = std::min(min_edge_len,
                              calc_dist(m, cell[ie], cell[(ie+1)%4]));
  }
  // If two edges are sufficiently close, they are interpreted as the same and
  // so internal. That leaves just the external edges.
  for (Int ic0 = 0; ic0 < ncell; ++ic0) {
    const auto cell0 = slice(m.e, ic0);
    for (Int ie0 = 0; ie0 < 4; ++ie0) {
      bool fnd_match = false;
      for (Int ic1 = 0; ic1 < ncell; ++ic1) {
        if (ic1 == ic0) continue;
        const auto cell1 = slice(m.e, ic1);
        for (Int ie1 = 0; ie1 < 4; ++ie1) {
          const auto edist = 
            std::min(calc_dist(m, cell0[ie0], cell1[ie1]) +
                     calc_dist(m, cell0[(ie0+1)%4], cell1[(ie1+1)%4]),
                     calc_dist(m, cell0[ie0], cell1[(ie1+1)%4]) +
                     calc_dist(m, cell0[(ie0+1)%4], cell1[ie1]));
          if (edist < 0.01*min_edge_len) {
            fnd_match = true;
            break;
          }
        }
        if (fnd_match) break;
      }
      if ( ! fnd_match) external_edges.push_back(4*ic0 + ie0);
    }
  }
  const Int nee = external_edges.size();
  slmm_throw_if((nslices(m.e) == 9 && nee != 12) ||
                (nslices(m.e) == 8 && nee != 10),
                "The local mesh has " << nslices(m.e)
                << " edges but fill_perim found " << nee
                << " external edges. Mesh:\n"
                << to_string(m));
  // Order edges so the leading point of edge i is the same as the trailing
  // point of edge (i+1) % nee.
  for (Int i = 0; i < nee-2; ++i) {
    const Int ic = external_edges[i] / 4, ie = external_edges[i] % 4;
    const auto icell = slice(m.e, ic);
    Real min_dist = 10;
    Int min_j = -1;
    for (Int j = i+1; j < nee; ++j) {
      const Int jc = external_edges[j] / 4, je = external_edges[j] % 4;
      const auto jcell = slice(m.e, jc);
      // External edges are all oriented CCW, so the leading point of one
      // matches the trailing point of the next.
      const auto d = calc_dist(m, icell[(ie+1)%4], jcell[je]);
      if (d < min_dist) {
        min_dist = d;
        min_j = j;
      }
    }
    if (min_j != i+1)
      std::swap(external_edges[i+1], external_edges[min_j]);
  }
}

void fill_perim (const Mesh& m, MeshNearestPointData<ko::HostSpace>& d) {
  std::vector<Int> external_edges;
  find_external_edges(m, external_edges);
  const Int nee = external_edges.size();
  d.perimp = MeshNearestPointData<ko::HostSpace>::Ints("perimp", nee);
  d.perimnml = MeshNearestPointData<ko::HostSpace>::Ints("perimnml", nee);
  for (Int k = 0; k < nee; ++k) {
    const auto i = external_edges[k];
    const auto ic = i / 4, ie = i % 4;
    d.perimp(k) = m.e(ic, ie);
    d.perimnml(k) = m.en(ic, ie);
  }
}

template <typename V3>
void calc_approx_nearest_point_on_arc (
  const V3& p0, const V3& p1, const V3& nml, const Real* point, Real* nearest)
{
  // Approximation is to project point onto the plane of the arc, then to
  // normalize to the arc. If point is on the arc, then nearest = point, as
  // desired.
  using geo = siqk::SphereGeometry;
  const auto dot = geo::dot(point, nml);
  for (Int d = 0; d < 3; ++d) nearest[d] = point[d] - dot*nml[d];
  geo::normalize(nearest);
  // If this initial value for nearest is outside of the arc, then 'nearest' is
  // set to the nearer of p0 and p1.
  //   Find the nearest point on line in parameterized coord alpha. alpha is in
  // [0,1] if 'nearest' is on the arc segment.
  Real pdiff[3];
  for (Int d = 0; d < 3; ++d) pdiff[d] = p0[d] - p1[d];
  Real p1mn[3];
  for (Int d = 0; d < 3; ++d) p1mn[d] = nearest[d] - p1[d];
  const Real alpha = geo::dot(p1mn, pdiff) / geo::norm2(pdiff);
  if      (alpha <= 0) for (Int d = 0; d < 3; ++d) nearest[d] = p1[d];
  else if (alpha >= 1) for (Int d = 0; d < 3; ++d) nearest[d] = p0[d];
}

void calc (const Mesh& m, const MeshNearestPointData<ko::HostSpace>& d,
           Real* v) {
  using geo = siqk::SphereGeometry;
  const Int nedge = d.perimp.size();
  const auto canpoa = [&] (const Int& ie, Real* vn) {
    calc_approx_nearest_point_on_arc(slice(m.p, d.perimp(ie)),
                                     slice(m.p, d.perimp((ie+1) % nedge)),
                                     slice(m.nml, d.perimnml(ie)),
                                     v, vn);
  };
  Real min_dist2 = 100;
  Int min_ie = -1;
  for (Int ie = 0; ie < nedge; ++ie) {
    Real vn[3];
    canpoa(ie, vn);
    Real d[3];
    for (Int i = 0; i < 3; ++i) d[i] = vn[i] - v[i];
    const Real dist2 = geo::norm2(d);
    if (dist2 < min_dist2) {
      min_ie = ie;
      min_dist2 = dist2;
    }
  }
  {
    Real vn[3];
    canpoa(min_ie, vn);
    for (Int d = 0; d < 3; ++d) v[d] = vn[d];
  }
}

Int test_canpoa () {
  Int nerr = 0;
  using geo = siqk::SphereGeometry;
  Real p0[] = {-0.5,0.5,0}, p1[] = {0.1,0.8,0}, nml[] = {0,0,-1};
  geo::normalize(p0);
  geo::normalize(p1);
  const Real points[][3] = {{-100,1,7}, {0,1,-7}, {11,-1,7}};
  Real correct[3][3];
  geo::copy(correct[0], p0);
  correct[1][0] = 0; correct[1][1] = 1; correct[1][2] = 0;
  geo::copy(correct[2], p1);
  for (Int i = 0; i < 3; ++i) {
    Real nr[3];
    nearest_point::calc_approx_nearest_point_on_arc(p0, p1, nml, points[i], nr);
    Real d[3];
    geo::axpbyz(1, nr, -1, correct[i], d);
    if (std::sqrt(geo::norm2(d)) > 10*std::numeric_limits<Real>::epsilon())
      ++nerr;
  }
  return nerr;
}

Int test_calc (const Mesh& m, const Int& tgt_ic,
               const MeshNearestPointData<ko::HostSpace>& md)
{
  Int nerr = 0;
  // If a point is on the perim, then calc should return it as the nearest
  // point.
  const Int np_perim = len(md.perimp);
  for (Int k = 0; k < np_perim; ++k) {
    const auto p0 = slice(m.p, md.perimp(k));
    const auto p1 = slice(m.p, md.perimp((k+1) % np_perim));
    const auto L = calc_dist(p0, p1);
    Real v[3];
    for (Int d = 0; d < 3; ++d) v[d] = p0[d];
    calc(m, md, v);
    if (calc_dist(p0, v) > std::numeric_limits<Real>::epsilon() / L) ++nerr;
    Real on[3];
    siqk::SphereGeometry::combine(p0, p1, 0.4, on);
    siqk::SphereGeometry::normalize(on);
    for (Int d = 0; d < 3; ++d) v[d] = on[d];
    calc(m, md, v);
    if (calc_dist(on, v) > 1e2*std::numeric_limits<Real>::epsilon() / L) ++nerr;
  }
  return nerr;
}

Int test_fill_perim (const Mesh& m, const Int& tgt_ic,
                     const MeshNearestPointData<ko::HostSpace>& md) {
  using geo = siqk::SphereGeometry;
  Int nerr = 0;
  const auto tgt_cell = slice(m.e, tgt_ic);
  const Int np_perim = len(md.perimp);
  { // Check that no target-cell point is on the perim.
    Int on = 0;
    for (Int ip = 0; ip < 4; ++ip) {
      const auto p = slice(m.p, tgt_cell[ip]);
      for (Int k = 0; k < np_perim; ++k) {
        const auto perim = slice(m.p, md.perimp(k));
        Real diff = 0;
        for (Int i = 0; i < 3; ++i) diff += std::abs(p[i] - perim[i]);
        if (diff == 0) ++on;
      }
    }
    if (on > 0) ++nerr;
  }
  { // Check that adjacent perimeter points agree with the edge's normal.
    for (Int k = 0; k < np_perim; ++k) {
      const auto p0 = slice(m.p, md.perimp(k));
      const auto p1 = slice(m.p, md.perimp((k+1) % np_perim));
      const auto nml = slice(m.nml, md.perimnml(k));
      const auto dot = geo::dot_c_amb(nml, p1, p0);
      const auto L = calc_dist(p0, p1);
      // Normal orthogonal to edge?
      if (std::abs(dot) > 1e2*std::numeric_limits<Real>::epsilon()/L) ++nerr;
      // Oriented inwards?
      Real outward[3];
      geo::cross(p1, p0, outward);
      if (geo::dot(nml, outward) >= 0) ++nerr;
    }
  }
  return nerr;
}
} // namespace nearest_point

typedef nearest_point::MeshNearestPointData<ko::HostSpace> MeshNearestPointData;

void init_nearest_point_data (const Mesh& m, MeshNearestPointData& d) {
  nearest_point::fill_perim(m, d);
}

int get_nearest_point (const Mesh& m, const MeshNearestPointData& d,
                       Real* v, const Int my_ic) {
  nearest_point::calc(m, d, v);
  return get_src_cell(m, v, my_ic);
}

Int unittest (const slmm::Mesh& m, const Int& tgt_elem) {
  Int nerr = 0;
  const Int nc = len(m.e);
  for (Int ic = 0; ic < nc; ++ic) {
    const auto cell = slice(m.e, ic);
    static const Real alphas[] = { 0.01, 0.99, 0, 1 };
    static const int nalphas = sizeof(alphas)/sizeof(*alphas);
    for (Int i = 0; i < nalphas; ++i)
      for (Int j = 0; j < nalphas; ++j) {
        const Real a = alphas[i], oma = 1-a, b = alphas[j], omb = 1-b;
        Real v[3] = {0};
        for (Int d = 0; d < 3; ++d)
          v[d] = (b  *(a*m.p(cell[0], d) + oma*m.p(cell[1], d)) + 
                  omb*(a*m.p(cell[3], d) + oma*m.p(cell[2], d)));
        if (i < 2 && j < 2) {
          if (get_src_cell(m, v, tgt_elem) != ic) ++nerr;
        } else {
          if (get_src_cell(m, v, tgt_elem) == -1) ++nerr;
        }
      }
  }
  nerr += nearest_point::test_canpoa();
  {
    MeshNearestPointData d;
    init_nearest_point_data(m, d);
    nerr += nearest_point::test_fill_perim(m, tgt_elem, d);
    nerr += nearest_point::test_calc(m, tgt_elem, d);
  }
  return nerr;
}

struct Advecter {
  typedef std::shared_ptr<Advecter> Ptr;

  struct Alg {
    enum Enum {
      jct,             // Cell-integrated Jacobian-combined transport.
      qos,             // Cell-integrated quadrature-on-sphere transport.
      csl_gll,         // Intended to mimic original Fortran CSL.
      csl_gll_subgrid, // Stable np=4 subgrid bilinear interp.
      csl_gll_exp,     // Stabilized np=4 subgrid reconstruction.
    };
    static Enum convert (Int alg) {
      switch (alg) {
      case 2: case 29:  return jct;
      case 3: case 39:  return qos;
      case 10: case 18: return csl_gll;
      case 11: case 17: return csl_gll_subgrid;
      case 12: case 19: return csl_gll_exp;
      default: slmm_throw_if(true, "transport_alg " << alg << " is invalid.");
      }
    }
    static bool is_cisl (const Enum& alg) { return alg == jct || alg == qos; }
  };

  Advecter (const Int np, const Int nelem, const Int transport_alg,
            const Int nearest_point_permitted_lev_bdy)
    : alg_(Alg::convert(transport_alg)),
      np_(np), np2_(np*np), np4_(np2_*np2_),
      tq_order_(alg_ == Alg::qos ? 14 : 12),
      nearest_point_permitted_lev_bdy_(nearest_point_permitted_lev_bdy)
  {
    local_mesh_.resize(nelem);
    local_mesh_tgt_elem_.resize(nelem);
    if (Alg::is_cisl(alg_))
      mass_mix_.resize(np4_);
    if (nearest_point_permitted_lev_bdy_ >= 0)
      local_mesh_nearest_point_data_.resize(nelem);
  }

  Int np  () const { return np_ ; }
  Int np2 () const { return np2_; }
  Int np4 () const { return np4_; }
  Int nelem () const { return local_mesh_.size(); }
  Int tq_order () const { return tq_order_; }
  Alg::Enum alg () const { return alg_; }
  bool is_cisl () const { return Alg::is_cisl(alg_); }

  template <typename Array3D>
  void init_local_mesh_if_needed (const Int ie, const Array3D& corners,
                                  const Real* p_inside) {
    slmm_assert(ie < static_cast<Int>(local_mesh_.size()));
    if (local_mesh_[ie].p.extent_int(0) != 0) return;
    auto& m = local_mesh_[ie];
    const Int
      nd = 3,
      nvert = corners.extent_int(1),
      ncell = corners.extent_int(2),
      N = nvert*ncell;
    m.p = typename Mesh::RealArray("p", N);
    m.e = typename Mesh::IntArray("e", ncell, nvert);
    for (Int ci = 0, k = 0; ci < ncell; ++ci)
      for (Int vi = 0; vi < nvert; ++vi, ++k) {
        for (int j = 0; j < nd; ++j)
          m.p(k,j) = corners(j,vi,ci);
        m.e(ci,vi) = k;
      }
    siqk::test::fill_normals<siqk::SphereGeometry>(m);
    local_mesh_tgt_elem_[ie] = slmm::get_src_cell(m, p_inside);
    slmm_assert(local_mesh_tgt_elem_[ie] >= 0 &&
                local_mesh_tgt_elem_[ie] < ncell);
    if (nearest_point_permitted_lev_bdy_ >= 0)
      init_nearest_point_data(local_mesh_[ie],
                              local_mesh_nearest_point_data_[ie]);
  }

  // Check that our ref <-> sphere map agrees with Homme's. p_homme is a GLL
  // point on the sphere. Check that we map it to a GLL ref point.
  void check_ref2sphere (const Int ie, const Real* p_homme) {
    const auto& m = local_mesh(ie);
    const auto& me = local_mesh_tgt_elem(ie);
    Real ref_coord[2];
    siqk::sqr::Info info;
    siqk::sqr::calc_sphere_to_ref(m.p, slice(m.e, me), p_homme,
                                  ref_coord[0], ref_coord[1], &info);
    const slmm::Basis basis(4, 0);
    const slmm::GLL gll;
    const Real* x, * wt;
    gll.get_coef(basis, x, wt);
    int fnd[2] = {0};
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 2; ++j)
        if (std::abs(ref_coord[j] - x[i]) < 1e-12)
          fnd[j] = 1;
    if ( ! fnd[0] || ! fnd[1])
      printf("COMPOSE check_ref2sphere: %1.15e %1.15e %d %d\n",
             ref_coord[0], ref_coord[1], info.success, info.n_iterations);
  }

  void init_M_tgt_if_needed();

  const Mesh& local_mesh (const Int ie) const {
    slmm_assert(ie < static_cast<Int>(local_mesh_.size()));
    return local_mesh_[ie];
  };

  Int local_mesh_tgt_elem (const Int ie) const {
    return local_mesh_tgt_elem_[ie];
  }

  const MeshNearestPointData& nearest_point_data (const Int ie) const {
    slmm_assert(ie < static_cast<Int>(local_mesh_nearest_point_data_.size()));
    return local_mesh_nearest_point_data_[ie];
  }

  std::vector<Real>& rhs_buffer (const Int qsize) {
    rhs_.resize(np2_*qsize);
    return rhs_;
  }

  std::vector<Real>& mass_mix_buffer () { return mass_mix_; }
  const Real* M_tgt (const Int& ie) {
    slmm_assert(ie >= 0 && ie < nelem());
    return alg_ == Alg::jct ?
      mass_tgt_.data() :
      mass_tgt_.data() + ie*np4_;
  }

  bool nearest_point_permitted (const Int& lev) {
    return lev <= nearest_point_permitted_lev_bdy_;
  }

private:
  const Alg::Enum alg_;
  const Int np_, np2_, np4_;
  std::vector<Mesh> local_mesh_;
  std::vector<Int> local_mesh_tgt_elem_;
  // For CISL:
  const Int tq_order_;
  std::vector<Real> mass_tgt_, mass_mix_, rhs_;
  // For recovery from get_src_cell failure:
  Int nearest_point_permitted_lev_bdy_;
  std::vector<MeshNearestPointData> local_mesh_nearest_point_data_;
};

void Advecter::init_M_tgt_if_needed ()  {
  if ( ! Alg::is_cisl(alg_)) return;
  if ( ! mass_tgt_.empty()) return;

  siqk::TriangleQuadrature tq;
  siqk::RawConstVec3s tq_bary;
  siqk::RawConstArray tq_w;
  tq.get_coef(tq_order(), tq_bary, tq_w);
  const Int nq = slmm::len(tq_w);

  const slmm::Basis basis(np_, 0);
  slmm::GLL gll;
  Real tgj[slmm::GLL::np_max], tgi[slmm::GLL::np_max];

  slmm_assert(Alg::is_cisl(alg_));
  if (alg_ == Alg::jct) {
    mass_tgt_.resize(np4_);
    Real* M_tgt = mass_tgt_.data();
    for (Int i = 0; i < np4_; ++i)
      M_tgt[i] = 0;

    const Real* const ps = siqk::sqr::get_ref_vertices();
    for (Int k = 1; k <= 2; ++k) {
      const Real t_area = siqk::PlaneGeometry::calc_tri_jacobian(
        ps, ps+2*k, ps+2*(k+1));

      for (Int qp = 0; qp < nq; ++qp) {
        Real tvo_coord[2];
        siqk::PlaneGeometry::bary2coord(
          ps, ps+2*k, ps+2*(k+1), slice(tq_bary, qp), tvo_coord);
        gll.eval(basis, tvo_coord[0], tgi);
        gll.eval(basis, tvo_coord[1], tgj);

        const Real d0 = 0.5 * t_area * tq_w[qp];
        for (Int aj = 0, a_basis_idx = 0; aj < np_; ++aj) {
          const Real d1 = d0 * tgj[aj];
          for (Int ai = 0; ai < np_; ++ai, ++a_basis_idx) {
            Real d2 = d1 * tgi[ai];
            for (Int b_basis_idx = a_basis_idx; b_basis_idx < np2_; ++b_basis_idx) {
              const Int bj = b_basis_idx / np_;
              const Int bi = b_basis_idx % np_;
              M_tgt[np2_*a_basis_idx + b_basis_idx] += d2 * tgi[bi] *tgj[bj];
            }
          }
        }
      }
    }

    int info = slmm::dpotrf('L', np2_, M_tgt, np2_);
    slmm_assert(info == 0);
  } else {
    const Int ne = nelem();
    slmm_assert(ne > 0);
    const Int sz = np4_ * ne;
    mass_tgt_.resize(sz);
    Real* M_tgt = mass_tgt_.data();

    for (Int ie = 0; ie < ne; ++ie) {
      // Local mesh.
      const auto& mesh = local_mesh_[ie];
      // Target cell in this local mesh.
      const auto ci = local_mesh_tgt_elem_[ie];
      slmm_assert(ci >= 0 && ci < nslices(mesh.e));
      const auto cell = slice(mesh.e, ci);
      // Sphere coordinates of cell corners.
      Real ps[12];
      copy_vertices(mesh.p, mesh.e, ci, ps);

      for (Int i = 0; i < np4_; ++i)
        M_tgt[i] = 0;
      for (Int k = 1; k <= 2; ++k) {
        for (Int qp = 0; qp < nq; ++qp) {
          // Ref -> sphere.
          Real sphere_coord[3];
          const Real jac_reftri2sphere = siqk::SphereGeometry::calc_tri_jacobian(
            ps, ps+3*k, ps+3*(k+1), slice(tq_bary, qp), sphere_coord);
          // Eval basis functions.
          Real a, b;
          siqk::sqr::calc_sphere_to_ref(mesh.p, cell, sphere_coord, a, b);
          gll.eval(basis, a, tgi);
          gll.eval(basis, b, tgj);
          // Integrate.
          const Real d0 = 0.5 * tq_w[qp] * jac_reftri2sphere;
          for (Int aj = 0, a_basis_idx = 0; aj < np_; ++aj) {
            const Real d1 = d0 * tgj[aj];
            for (Int ai = 0; ai < np_; ++ai, ++a_basis_idx) {
              Real d2 = d1 * tgi[ai];
              for (Int b_basis_idx = a_basis_idx; b_basis_idx < np2_; ++b_basis_idx) {
                const Int bj = b_basis_idx / np_;
                const Int bi = b_basis_idx % np_;
                M_tgt[np2_*a_basis_idx + b_basis_idx] += d2 * tgi[bi] *tgj[bj];
              }
            }
          }
        }
      }

      int info = slmm::dpotrf('L', np2_, M_tgt, np2_);
      slmm_assert(info == 0);

      M_tgt += np4_;
    }
  }
}
} // namespace slmm

#ifdef SLMM_MAIN
/*
  mpicxx -O3 -g -DSLMM_MAIN -Wall -pedantic -fopenmp -std=c++11 -I/home/ambradl/lib/kokkos/cpu/include slmm.cpp -L/home/ambradl/lib/kokkos/cpu/lib -lkokkos -ldl -llapack -lblas; if [ $? == 0 ]; then OMP_PROC_BIND=false OMP_NUM_THREADS=2 mpirun -np 14 ./a.out; fi
 */
int main (int argc, char** argv) {
  int ret = 0;
  Kokkos::initialize(argc, argv);
  try {
    int nerr = 0;
#ifdef BUILD_CISL
    nerr += siqk::unittest();
#endif
    nerr += slmm::unittest();
    std::cerr << (nerr ? "FAIL" : "PASS") << "ED\n";
  } catch (const std::exception& e) {
    std::cerr << e.what();
  }
  Kokkos::finalize();
  return ret;
}
#endif

#ifndef SLMM_MAIN
# ifdef HAVE_CONFIG_H
#  include "config.h.c"
# endif
# ifdef CAM
#  define QSIZE_D PCNST
# endif
#else
# define QSIZE_D 64
#endif

namespace homme {
typedef slmm::Int Int;
typedef slmm::Real Real;

namespace ko = Kokkos;

// Fortran array wrappers.
template <typename T> using FA2 =
  Kokkos::View<T**,    Kokkos::LayoutLeft, Kokkos::HostSpace>;
template <typename T> using FA3 =
  Kokkos::View<T***,   Kokkos::LayoutLeft, Kokkos::HostSpace>;
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

static slmm::Advecter::Ptr g_advecter;

void slmm_init (const Int np, const Int nelemd, const Int transport_alg,
                const Int sl_nearest_point_lev) {
  g_advecter = std::make_shared<slmm::Advecter>(np, nelemd, transport_alg,
                                                sl_nearest_point_lev);
}

#ifdef BUILD_CISL
// On HSW, this demonstrably is not needed; speedup in slmm_project_np4 is a
// result of loop details and local memory. However, I'll keep this in case it's
// of use on KNL.
#define IR_ALIGN(decl) Real decl __attribute__((aligned(64)))

void slmm_project_np4 (
  const Int lev, const Int nets, const Int ie,
  const Int nnc, const Int np, const Int nlev, const Int qsize,
  const Real* dep_points_r, const Real* Qj_src_r, const Real* jac_tgt_r,
  const Real* rho_tgt_r, const Int tl_np1, Real* q_r, Real* q_min_r, Real* q_max_r)
{
  slmm_assert(g_advecter);
  using slmm::slice;

  static constexpr Int max_num_vertices = 4;
  // Convex quad-quad intersection yields <= 8 intersection points.
  static constexpr Int max_num_intersections = 8;
  static constexpr Int my_np = 4, my_np2 = my_np*my_np,
    my_np4 = my_np2*my_np2;
  static constexpr Int max_qsize = QSIZE_D;
  const Int ie0 = ie - nets;

  FA2<const Real>
    jac_tgt(jac_tgt_r, np, np);
  FA3<const Real>
    dep_points(reinterpret_cast<const Real*>(dep_points_r), 3, np, np);
  FA4<const Real>
    Qj_src(Qj_src_r, np, np, qsize+1, nnc),
    rho_tgt(rho_tgt_r, np, np, nlev, tl_np1+1);
  FA4<Real>
    q(q_r, np, np, nlev, qsize);
  FA5<Real>
    q_min(q_min_r, np, np, nlev, qsize, ie0+1),
    q_max(q_max_r, np, np, nlev, qsize, ie0+1);

  const Real* const M_tgt = g_advecter->M_tgt(ie);
  IR_ALIGN(M_mix[my_np4]);
  IR_ALIGN(rhs_r[my_np2*max_qsize]);
  IR_ALIGN(vi_buf[3*max_num_vertices]);
  IR_ALIGN(vo_buf[3*max_num_intersections]);
  IR_ALIGN(wrk_buf[4*max_num_intersections]);
  IR_ALIGN(sgj[slmm::GLL::np_max]);
  IR_ALIGN(sgi[slmm::GLL::np_max]);
  IR_ALIGN(tgj[slmm::GLL::np_max]);
  IR_ALIGN(tgi[slmm::GLL::np_max]);
  IR_ALIGN(svo_coord[2]);
  IR_ALIGN(tvo_coord[2]);

  const Int np2 = np*np;
  const slmm::Basis basis(np, 0);
  const slmm::GLL gll;
  siqk::TriangleQuadrature tq;
  siqk::RawConstVec3s tq_bary;
  siqk::RawConstArray tq_w;
  tq.get_coef(g_advecter->tq_order(), tq_bary, tq_w);
  const Int nq = slmm::len(tq_w);

  for (Int i = 0, n = qsize*np2; i < n; ++i) rhs_r[i] = 0;
  FA3<Real> rhs(rhs_r, np, np, qsize);

  const auto& m = g_advecter->local_mesh(ie);
  bool first = true;
  const auto alg = g_advecter->alg();
  for (Int sci = 0; sci < nnc; ++sci) { // For each source cell:
    // Intersect.
    const siqk::RawVec3s vi(vi_buf, max_num_vertices);
    siqk::RawVec3s vo(vo_buf, max_num_intersections);
    Int no;
    {
      siqk::RawVec3s wrk(wrk_buf, max_num_intersections);
      static const Int crnr[][2] = {{0, 0}, {np-1, 0}, {np-1, np-1}, {0, np-1}};
      for (int i = 0; i < max_num_vertices; ++i)
        siqk::SphereGeometry::copy(
          slice(vi, i),
          ko::subview(dep_points, ko::ALL(), crnr[i][0], crnr[i][1]));
      siqk::sh::clip_against_poly<siqk::SphereGeometry>(
        g_advecter->local_mesh(ie), sci, vi, max_num_vertices, vo, no, wrk);
    }

    // If the intersection is empty, move on.
    if (no == 0) continue;
    slmm_assert(no <= max_num_intersections);

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

    for (Int i = 0; i < my_np4; ++i)
      M_mix[i] = 0;
    if (alg == slmm::Advecter::Alg::jct) {
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
      
      for (Int ktri = 1; ktri < no-1; ++ktri) { // triangles in vo
        const Real s_area = siqk::PlaneGeometry::calc_tri_jacobian(
          svo, svo+2*ktri, svo+2*(ktri+1));

        for (Int qp = 0; qp < nq; ++qp) { // quad point
          siqk::PlaneGeometry::bary2coord(tvo, tvo+2*ktri, tvo+2*(ktri+1),
                                          slice(tq_bary, qp), tvo_coord);
          slmm::gll_np4_eval(tvo_coord[0], tgi);
          slmm::gll_np4_eval(tvo_coord[1], tgj);
          siqk::PlaneGeometry::bary2coord(svo, svo+2*ktri, svo+2*(ktri+1),
                                          slice(tq_bary, qp), svo_coord);
          slmm::gll_np4_eval(svo_coord[0], sgi);
          slmm::gll_np4_eval(svo_coord[1], sgj);
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
    } else {
      static const Int tcell[] = {0, 1, 2, 3};
      const auto scell = slice(m.e, sci);
      for (Int ktri = 1; ktri < no-1; ++ktri) {
        for (Int qp = 0; qp < nq; ++qp) {
          Real sphere_coord[3];
          const Real jac_reftri2sphere = siqk::SphereGeometry::calc_tri_jacobian(
            slice(vo, 0), slice(vo, ktri), slice(vo, ktri+1), slice(tq_bary, qp),
            sphere_coord);
          Real ta, tb, sa, sb;
          siqk::sqr::calc_sphere_to_ref(vi, tcell, sphere_coord, ta, tb);
          gll.eval(basis, ta, tgi);
          gll.eval(basis, tb, tgj);
          siqk::sqr::calc_sphere_to_ref(m.p, scell, sphere_coord, sa, sb);
          gll.eval(basis, sa, sgi);
          gll.eval(basis, sb, sgj);
          Real d0 = 0.5 * tq_w[qp] * jac_reftri2sphere;
          for (Int tj = 0, t_basis_idx = 0; tj < my_np; ++tj) {
            const Real d1 = d0 * tgj[tj];
            for (Int ti = 0; ti < my_np; ++ti, ++t_basis_idx) {
              Real d2 = d1 * tgi[ti];
              for (Int sj = 0, s_basis_idx = 0; sj < my_np; ++sj) {
                const Real d3 = d2 * sgj[sj];
                for (Int si = 0; si < my_np; ++si, ++s_basis_idx) {
                  const Real d = d3 * sgi[si];
                  M_mix[my_np2*t_basis_idx + s_basis_idx] += d;
                }
              }
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
  int info = slmm::dpotrs('L', np2, qsize, M_tgt, np2, rhs_r, np2);
  slmm_assert(info == 0);

  // Load q.
  if (alg == slmm::Advecter::Alg::jct) {
    for (Int tq = 0; tq < qsize; ++tq)
      for (Int j = 0; j < np; ++j)
        for (Int i = 0; i < np; ++i)
          q(i,j,lev,tq) = rhs(i,j,tq) / (rho_tgt(i,j,lev,tl_np1) * jac_tgt(i,j));
  } else {
    for (Int tq = 0; tq < qsize; ++tq)
      for (Int j = 0; j < np; ++j)
        for (Int i = 0; i < np; ++i)
          q(i,j,lev,tq) = rhs(i,j,tq) / rho_tgt(i,j,lev,tl_np1);
  }

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

// This version is not thread safe.
void slmm_project (
  const Int lev, const Int nets, const Int ie,
  const Int nnc, const Int np, const Int nlev, const Int qsize,
  // Geometry.
  const Cartesian3D* dep_points_r,     // dep_points(1:np, 1:np)
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
  slmm_assert(g_advecter);

  if (np == 4) {
    // This version is thread safe.
    slmm_project_np4(lev, nets, ie, nnc, np, nlev, qsize, 
                     reinterpret_cast<const Real*>(dep_points_r),
                     Qj_src_r, jac_tgt_r, rho_tgt_r, tl_np1,
                     q_r, q_min_r, q_max_r);
    return;
  }

  using slmm::slice;

  static constexpr Int max_num_vertices = 4;
  // Convex quad-quad intersection yields <= 8 intersection points.
  static constexpr Int max_num_intersections = 8;
  const Int ie0 = ie - nets;

  // In rho_tgt, q, q_{min,max}, the fact that Fortran uses layout left means we
  // don't need to know the rightmost dim.
  FA2<const Real>
    jac_tgt(jac_tgt_r, np, np);
  FA3<const Real>
    dep_points(reinterpret_cast<const Real*>(dep_points_r), 3, np, np);
  FA4<const Real>
    Qj_src(Qj_src_r, np, np, qsize+1, nnc),
    rho_tgt(rho_tgt_r, np, np, nlev, tl_np1+1);
  FA4<Real>
    q(q_r, np, np, nlev, qsize);
  FA5<Real>
    q_min(q_min_r, np, np, nlev, qsize, ie0+1),
    q_max(q_max_r, np, np, nlev, qsize, ie0+1);

  const Int np2 = np*np, np4 = np2*np2;
  const slmm::Basis basis(np, 0);
  const slmm::GLL gll;
  siqk::TriangleQuadrature tq;
  siqk::RawConstVec3s tq_bary;
  siqk::RawConstArray tq_w;
  tq.get_coef(g_advecter->tq_order(), tq_bary, tq_w);
  const Int nq = slmm::len(tq_w);

  auto& rhs_r = g_advecter->rhs_buffer(qsize);
  for (Int i = 0, n = qsize*np2; i < n; ++i) rhs_r[i] = 0;
  FA3<Real> rhs(rhs_r.data(), np, np, qsize);
  const Real* const M_tgt = g_advecter->M_tgt(ie);
  auto& M_mix = g_advecter->mass_mix_buffer();

  const auto& m = g_advecter->local_mesh(ie);
  bool first = true;
  for (Int sci = 0; sci < nnc; ++sci) { // For each source cell:
    // Intersect.
    Real buf[3*max_num_vertices + 3*max_num_intersections + 4*max_num_intersections];
    const siqk::RawVec3s vi(buf, max_num_vertices);
    siqk::RawVec3s vo(buf + 3*max_num_vertices, max_num_intersections);
    Real* const wrkbuf = vo.data() + 3*max_num_intersections;
    Int no;
    {
      siqk::RawVec3s wrk(wrkbuf, max_num_intersections);
      static const Int crnr[][2] = {{0, 0}, {np-1, 0}, {np-1, np-1}, {0, np-1}};
      for (int i = 0; i < max_num_vertices; ++i)
        siqk::SphereGeometry::copy(
          slice(vi, i),
          ko::subview(dep_points, ko::ALL(), crnr[i][0], crnr[i][1]));
      siqk::sh::clip_against_poly<siqk::SphereGeometry>(
        g_advecter->local_mesh(ie), sci, vi, max_num_vertices, vo, no, wrk);
    }

    // If the intersection is empty, move on.
    if (no == 0) continue;
    slmm_assert(no <= max_num_intersections);

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
      Real
        sgj[slmm::GLL::np_max], sgi[slmm::GLL::np_max],
        tgj[slmm::GLL::np_max], tgi[slmm::GLL::np_max],
        svo_coord[2], tvo_coord[2];

      for (Int qp = 0; qp < nq; ++qp) { // quad point
        siqk::PlaneGeometry::bary2coord(tvo, tvo+2*ktri, tvo+2*(ktri+1),
                                        slice(tq_bary, qp), tvo_coord);
        gll.eval(basis, tvo_coord[0], tgi);
        gll.eval(basis, tvo_coord[1], tgj);

        siqk::PlaneGeometry::bary2coord(svo, svo+2*ktri, svo+2*(ktri+1),
                                        slice(tq_bary, qp), svo_coord);
        gll.eval(basis, svo_coord[0], sgi);
        gll.eval(basis, svo_coord[1], sgj);

        {
          const Real d0 = 0.5 * s_area * tq_w[qp];
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
  int info = slmm::dpotrs('L', np2, qsize, M_tgt, np2, rhs.data(), np2);
  slmm_assert(info == 0);

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
#endif // BUILD_CISL

template <int np>
void slmm_csl (
  const Int lev, const Int nets, const Int ie,
  const Int nnc, const Int nlev, const Int qsize,
  // Geometry.
  const Cartesian3D* dep_points_r,     // dep_points(1:3, 1:np, 1:np)
  // Fields.
  const Real* q_src_r,                 // q_src(1:np, 1:np, 1:qsize, 1:nnc)
  const Real* jac_tgt_r,               // jac_tgt(1:np, 1:np) [unused]
  // rho_tgt(1:np, 1:np, lev, tl_np1+1) [unused]
  const Real* rho_tgt_r, const Int tl_np1,
  // Output target tracer mixing ratio and bounds derived from the domain of
  // dependence.
  Real* q_r,                      // q(1:np, 1:np, lev, 1:qsize)
  Real* q_min_r, Real* q_max_r)   // q_{min,max}(1:np, 1:np, lev, 1:qsize, ie-nets+1)
{
  static_assert(np == 4, "SLMM CSL is supported for np 4 only.");
  using slmm::slice;

  const Int np2 = np*np;
  const Int ie0 = ie - nets;
  const auto alg = g_advecter->alg();
  FA3<const Real>
    dep_points(reinterpret_cast<const Real*>(dep_points_r), 3, np2);
  FA4<const Real>
    q_src(q_src_r, np, np, qsize+1, nnc);
  FA4<Real>
    q_tgt(q_r, np, np, nlev, qsize);
  FA5<Real>
    q_min(q_min_r, np, np, nlev, qsize, ie0+1),
    q_max(q_max_r, np, np, nlev, qsize, ie0+1);

  const slmm::Basis basis(np, 0);
  const slmm::GLL gll;

  const auto& m = g_advecter->local_mesh(ie);
  for (Int tvi = 0; tvi < np2; ++tvi) {
    const Int ti = tvi % np;
    const Int tj = tvi / np;

    // Determine which cell the departure point is in.
    const Real* dep_point = dep_points.data() + 3*tvi;
    const Int sci = slmm::get_src_cell(m, dep_point,
                                      g_advecter->local_mesh_tgt_elem(ie));
    slmm_throw_if(sci == -1, "Departure point is outside of halo.");

    // Get reference point.
    Real ref_coord[2];
    siqk::sqr::calc_sphere_to_ref(m.p, slice(m.e, sci), dep_point,
                                  ref_coord[0], ref_coord[1]);

    // Interpolate.
    Real rx[4], ry[4];
    switch (alg) {
    case slmm::Advecter::Alg::csl_gll:
      slmm::gll_np4_eval(ref_coord[0], rx);
      slmm::gll_np4_eval(ref_coord[1], ry);
      break;
    case slmm::Advecter::Alg::csl_gll_subgrid:
      slmm::gll_np4_subgrid_eval(ref_coord[0], rx);
      slmm::gll_np4_subgrid_eval(ref_coord[1], ry);
      break;
    case slmm::Advecter::Alg::csl_gll_exp:
      slmm::gll_np4_subgrid_exp_eval(ref_coord[0], rx);
      slmm::gll_np4_subgrid_exp_eval(ref_coord[1], ry);
      break;
    default:
      assert(0);
    }

    {
      const Real* const qs0 = &q_src(0,0,0,sci);
      for (Int q = 0; q < qsize; ++q) {
        const Real* const qs = qs0 + q*np*np;
        q_tgt(ti,tj,lev,q) =
          (ry[0]*(rx[0]*qs[ 0] + rx[1]*qs[ 1] + rx[2]*qs[ 2] + rx[3]*qs[ 3]) +
           ry[1]*(rx[0]*qs[ 4] + rx[1]*qs[ 5] + rx[2]*qs[ 6] + rx[3]*qs[ 7]) +
           ry[2]*(rx[0]*qs[ 8] + rx[1]*qs[ 9] + rx[2]*qs[10] + rx[3]*qs[11]) +
           ry[3]*(rx[0]*qs[12] + rx[1]*qs[13] + rx[2]*qs[14] + rx[3]*qs[15]));
        //todo Could rm redundant computations of this by storing array
        // of sci and then uniq'ing it.
        Real q_min_s, q_max_s;
        q_min_s = q_max_s = qs[0];
        for (int k = 1; k < 16; ++k) {
          q_min_s = std::min(q_min_s, qs[k]);
          q_max_s = std::max(q_max_s, qs[k]);
        }
        q_min(ti,tj,lev,q,ie0) = q_min_s;
        q_max(ti,tj,lev,q,ie0) = q_max_s;
      }
    }
  }
}

namespace mpi { //todo Share with cedr.
class Parallel {
  MPI_Comm comm_;
public:
  typedef std::shared_ptr<Parallel> Ptr;
  Parallel(MPI_Comm comm) : comm_(comm) {}
  MPI_Comm comm () const { return comm_; }
  Int size() const {
    int sz = 0;
    MPI_Comm_size(comm_, &sz);
    return sz;
  }
  Int rank() const {
    int pid = 0;
    MPI_Comm_rank(comm_, &pid);
    return pid;
  }
  Int root () const { return 0; }
  bool amroot () const { return rank() == root(); }
};

Parallel::Ptr make_parallel (MPI_Comm comm) {
  return std::make_shared<Parallel>(comm);
}

struct Request {
  MPI_Request request;

#ifdef COMPOSE_DEBUG_MPI
  int unfreed;
  Request();
  ~Request();
#endif
};

#ifdef COMPOSE_DEBUG_MPI
Request::Request () : unfreed(0) {}
Request::~Request () {
  if (unfreed) {
    std::stringstream ss;
    ss << "Request is being deleted with unfreed = " << unfreed;
    int fin;
    MPI_Finalized(&fin);
    if (fin) {
      ss << "\n";
      std::cerr << ss.str();
    } else {
      pr(ss.str());
    }
  }
}
#endif

template <typename T> MPI_Datatype get_type();
template <> inline MPI_Datatype get_type<int>() { return MPI_INT; }
template <> inline MPI_Datatype get_type<double>() { return MPI_DOUBLE; }
template <> inline MPI_Datatype get_type<long>() { return MPI_LONG_INT; }

template <typename T>
int isend (const Parallel& p, const T* buf, int count, int dest, int tag,
           Request* ireq) {
  MPI_Datatype dt = get_type<T>();
  MPI_Request ureq;
  MPI_Request* req = ireq ? &ireq->request : &ureq;
  int ret = MPI_Isend(const_cast<T*>(buf), count, dt, dest, tag, p.comm(), req);
  if ( ! ireq) MPI_Request_free(req);
#ifdef COMPOSE_DEBUG_MPI
  else ireq->unfreed++;
#endif
  return ret;
}

template <typename T>
int irecv (const Parallel& p, T* buf, int count, int src, int tag, Request* ireq) {
  MPI_Datatype dt = get_type<T>();
  MPI_Request ureq;
  MPI_Request* req = ireq ? &ireq->request : &ureq;
  int ret = MPI_Irecv(buf, count, dt, src, tag, p.comm(), req);
  if ( ! ireq) MPI_Request_free(req);
#ifdef COMPOSE_DEBUG_MPI
  else ireq->unfreed++;
#endif
  return ret;
}

int waitany (int count, Request* reqs, int* index, MPI_Status* stats = nullptr) {
#ifdef COMPOSE_DEBUG_MPI
  std::vector<MPI_Request> vreqs(count);
  for (int i = 0; i < count; ++i) vreqs[i] = reqs[i].request;
  const auto out = MPI_Waitany(count, vreqs.data(), index,
                               stats ? stats : MPI_STATUS_IGNORE);
  for (int i = 0; i < count; ++i) reqs[i].request = vreqs[i];
  reqs[*index].unfreed--;
  return out;
#else
  return MPI_Waitany(count, reinterpret_cast<MPI_Request*>(reqs), index,
                     stats ? stats : MPI_STATUS_IGNORE);
#endif
}

int waitall (int count, Request* reqs, MPI_Status* stats = nullptr) {
#ifdef COMPOSE_DEBUG_MPI
  std::vector<MPI_Request> vreqs(count);
  for (int i = 0; i < count; ++i) vreqs[i] = reqs[i].request;
  const auto out = MPI_Waitall(count, vreqs.data(),
                               stats ? stats : MPI_STATUS_IGNORE);
  for (int i = 0; i < count; ++i) {
    reqs[i].request = vreqs[i];
    reqs[i].unfreed--;
  }
  return out;
#else
  return MPI_Waitall(count, reinterpret_cast<MPI_Request*>(reqs),
                     stats ? stats : MPI_STATUS_IGNORE);
#endif
}
} // namespace mpi

namespace cslmpi {
// Meta and bulk data for the classical (interpolation) SL method with special
// comm pattern.

#define SLMM_BOUNDS_CHECK
#ifdef SLMM_BOUNDS_CHECK
# define slmm_assert_high(condition) slmm_assert(condition)
#else
# define slmm_assert_high(condition)
#endif

// FixedCapList, ListOfLists, and BufferLayoutArray are simple and somewhat
// problem-specific array data structures for use in CslMpi.
template <typename T>
struct FixedCapList {
  FixedCapList () : n_(0) {}
  FixedCapList (const Int& cap) { slmm_assert_high(cap >= 0); reset_capacity(cap); }
  void reset_capacity (const Int& cap, const bool also_size = false) {
    slmm_assert(cap >= 0);
    d_.resize(cap);
    n_ = also_size ? cap : 0;
  }
  void clear () { n_ = 0; }

  Int n () const { return n_; }
  Int size () const { return n_; }
  Int capacity () const { return d_.size(); }
  const T& operator() (const Int& i) const { slmm_assert_high(i >= 0 && i < n_); return d_[i]; }
  T& operator() (const Int& i) { slmm_assert_high(i >= 0 && i < n_); return d_[i]; }
  void inc () { ++n_; slmm_assert_high(n_ <= static_cast<Int>(d_.size())); }
  void inc (const Int& dn) { n_ += dn; slmm_assert_high(n_ <= static_cast<Int>(d_.size())); }

  const T* data () const { return d_.data(); }
  T* data () { return d_.data(); }  
  const T& back () const { slmm_assert_high(n_ > 0); return d_[n_-1]; }
  T& back () { slmm_assert_high(n_ > 0); return d_[n_-1]; }
  const T* begin () const { return d_.data(); }
  T* begin () { return d_.data(); }
  const T* end () const { return d_.data() + n_; }
  T* end () { return d_.data() + n_; }

 private:
  std::vector<T> d_;
  Int n_;
};

template <typename T>
struct ListOfLists {
  struct List {
    Int n () const { return n_; }

    T& operator() (const Int& i) { slmm_assert_high(i >= 0 && i < n_); return d_[i]; }
    const T& operator() (const Int& i) const { slmm_assert_high(i >= 0 && i < n_); return d_[i]; }

    const T* data () const { return d_; }
    T* data () { return d_; }
    const T* begin () const { return d_; }
    T* begin () { return d_; }
    const T* end () const { return d_ + n_; }
    T* end () { return d_ + n_; }

    void zero () { for (Int i = 0; i < n_; ++i) d_[i] = 0; }

  private:
    friend class ListOfLists<T>;
    List (T* d, const Int& n) : d_(d), n_(n) { slmm_assert_high(n_ >= 0); }
    T* const d_;
    const Int n_;
  };

  ListOfLists () {}
  ListOfLists (const Int nlist, const Int* nlist_per_list) { init(nlist, nlist_per_list); }
  void init (const Int nlist, const Int* nlist_per_list) {
    slmm_assert(nlist >= 0); 
    ptr_.resize(nlist+1);
    ptr_[0] = 0;
    for (Int i = 0; i < nlist; ++i) {
      slmm_assert(nlist_per_list[i] >= 0);
      ptr_[i+1] = ptr_[i] + nlist_per_list[i];
    }
    d_.resize(ptr_.back());
  }

  Int n () const { return static_cast<Int>(ptr_.size()) - 1; }
  List operator() (const Int& i) {
    slmm_assert_high(i >= 0 && i < static_cast<Int>(ptr_.size()) - 1);
    return List(&d_[ptr_[i]], ptr_[i+1] - ptr_[i]);
  }
  const List operator() (const Int& i) const {
    slmm_assert_high(i >= 0 && i < static_cast<Int>(ptr_.size()) - 1);
    return List(const_cast<T*>(&d_[ptr_[i]]), ptr_[i+1] - ptr_[i]);
  }
  T& operator() (const Int& i, const Int& j) {
    slmm_assert_high(i >= 0 && i < static_cast<Int>(ptr_.size()) - 1 &&
                     j >= 0 && j < ptr_[i+1] - ptr_[i]);
    return d_[ptr_[i] + j];
  }
  const T& operator() (const Int& i, const Int& j) const {
    slmm_assert_high(i >= 0 && i < static_cast<Int>(ptr_.size()) - 1 &&
                     j >= 0 && j < ptr_[i+1] - ptr_[i]);
    return d_[ptr_[i] + j];
  }

private:
  friend class BufferLayoutArray;
  std::vector<T> d_;
  std::vector<Int> ptr_;
  T* data () { return d_.data(); }
  const T* data () const { return d_.data(); }
};

struct LayoutTriple {
  Int xptr, qptr, cnt;
  LayoutTriple () : LayoutTriple(0) {}
  LayoutTriple (const Int& val) { xptr = qptr = cnt = 0; }
};

struct BufferLayoutArray {
  struct BufferRankLayoutArray {
    LayoutTriple& operator() (const Int& lidi, const Int& lev) {
      slmm_assert_high(lidi >= 0 && lev >= 0 && lidi*nlev_ + lev < d_.n());
      return d_(lidi*nlev_ + lev);
    }
    const LayoutTriple& operator() (const Int& lidi, const Int& lev) const {
      slmm_assert_high(lidi >= 0 && lev >= 0 && lidi*nlev_ + lev < d_.n());
      return d_(lidi*nlev_ + lev);
    }

  private:
    friend class BufferLayoutArray;
    BufferRankLayoutArray (const ListOfLists<LayoutTriple>::List& d, const Int& nlev)
      : d_(d), nlev_(nlev) {}
    ListOfLists<LayoutTriple>::List d_;
    Int nlev_;
  };

  BufferLayoutArray () : nlev_(0) {}
  BufferLayoutArray (const Int& nrank, const Int* nlid_per_rank, const Int& nlev)
    { init(nrank, nlid_per_rank, nlev); }
  void init (const Int& nrank, const Int* nlid_per_rank, const Int& nlev) {
    slmm_assert(nrank >= 0 && nlev > 0);
    nlev_ = nlev;
    std::vector<Int> ns(nrank);
    for (Int i = 0; i < nrank; ++i) {
      slmm_assert(nlid_per_rank[i] > 0);
      ns[i] = nlid_per_rank[i] * nlev;
    }
    d_.init(nrank, ns.data());
  }

  void zero () {
    const Int ni = d_.n();
#ifdef HORIZ_OPENMP
#   pragma omp for
#endif
    for (Int i = 0; i < ni; ++i) {
      auto&& l = d_(i);
      for (Int j = 0, nj = l.n(); j < nj; ++j)
        l(j) = 0;
    }
  }

  LayoutTriple& operator() (const Int& ri, const Int& lidi, const Int& lev) {
    slmm_assert_high(ri >= 0 && ri < d_.n() &&
                     lidi >= 0 && lev >= 0 &&
                     lidi*nlev_ + lev < d_(ri).n());
    return d_.data()[d_.ptr_[ri] + lidi*nlev_ + lev];
  }
  const LayoutTriple& operator() (const Int& ri, const Int& lidi, const Int& lev) const {
    return const_cast<BufferLayoutArray*>(this)->operator()(ri, lidi, lev);
  }
  BufferRankLayoutArray operator() (const Int& ri) {
    slmm_assert_high(ri >= 0 && ri < d_.n());
    return BufferRankLayoutArray(d_(ri), nlev_);
  }
  const BufferRankLayoutArray operator() (const Int& ri) const {
    slmm_assert_high(ri >= 0 && ri < d_.n());
    return BufferRankLayoutArray(d_(ri), nlev_);
  }

private:
  ListOfLists<LayoutTriple> d_;
  Int nlev_;
};

// Meta and bulk data for the interpolation SL MPI communication pattern.
struct CslMpi {
  typedef std::shared_ptr<CslMpi> Ptr;

  template <typename Datatype>
  using Array = ko::View<Datatype, ko::LayoutRight, ko::HostSpace>;
  typedef Array<Int**> IntArray2D;

  struct GidRank {
    Int
      gid,      // cell global ID
      rank,     // the rank that owns the cell
      rank_idx, // index into list of ranks with whom I communicate, including me
      lid_on_rank,    // the local ID of the cell on the owning rank
      lid_on_rank_idx; // index into list of LIDs on the rank
  };
  // The comm and real data associated with an element patch, the set of
  // elements surrounding an owned cell.
  struct ElemData {
    struct OwnItem {
      short lev;   // level index
      short k;     // linearized GLL index
    };
    struct RemoteItem {
      Int q_extrema_ptr, q_ptr; // pointers into recvbuf
      short lev, k;
    };
    GidRank* me;                     // the owned cell
    FixedCapList<GidRank> nbrs;      // the cell's neighbors (but including me)
    FixedCapList<OwnItem> own;       // points whose q are computed with own rank's data
    FixedCapList<RemoteItem> rmt;    // list computed by a remote rank's data
    IntArray2D src;                  // src(lev,k) = get_src_cell
    Array<Real**[2]> q_extrema;
    const Real* metdet, * qdp, * dp; // the owned cell's data
    Real* q;
  };

  const mpi::Parallel::Ptr p;
  const Int np, np2, nlev, qsize, qsized, nelemd;

  FixedCapList<ElemData> ed; // this rank's owned cells, indexed by LID

  // IDs.
  FixedCapList<Int> ranks, nx_in_rank, mylid_with_comm, mylid_with_comm_tid_ptr;
  ListOfLists <Int> nx_in_lid, lid_on_rank;
  BufferLayoutArray bla;

  // MPI comm data.
  ListOfLists<Real> sendbuf, recvbuf;
  FixedCapList<Int> sendcount;
  FixedCapList<mpi::Request> sendreq, recvreq;

  bool horiz_openmp;
#ifdef HORIZ_OPENMP
  ListOfLists<omp_lock_t> ri_lidi_locks;
#endif

  CslMpi (const mpi::Parallel::Ptr& ip, Int inp, Int inlev, Int iqsize,
          Int iqsized, Int inelemd)
    : p(ip), np(inp), np2(np*np), nlev(inlev), qsize(iqsize), qsized(iqsized),
      nelemd(inelemd)
  {
    slmm_assert(qsized <= QSIZE_D);
  }

  ~CslMpi () {
#ifdef HORIZ_OPENMP
    const Int nrmtrank = static_cast<Int>(ranks.n()) - 1;
    for (Int ri = 0; ri < nrmtrank; ++ri) {
      auto&& locks = ri_lidi_locks(ri);
      for (auto& lock: locks)
        omp_destroy_lock(&lock);  
    }
#endif
  }
};

// Fill in (gid, rank), the list of owning rank per gid.
void collect_gid_rank (CslMpi& cm, const Int* nbr_id_rank, const Int* nirptr) {
  cm.ed.reset_capacity(cm.nelemd, true);
  for (Int i = 0; i < cm.nelemd; ++i) {
    auto& ed = cm.ed(i);
    const Int* nir = nbr_id_rank + nirptr[i];
    const Int nnir = (nirptr[i+1] - nirptr[i]) / 2;
    const Int mygid = nir[0];
    ed.me = nullptr;
    ed.nbrs.reset_capacity(nnir-1, true);
    ed.own.reset_capacity(cm.nlev * cm.np2);
    ed.rmt.reset_capacity(cm.nlev * cm.np2);
    for (Int j = 1; j < nnir; ++j) {
      auto& n = ed.nbrs(j-1);
      n.gid = nir[2*j];
      if (n.gid == mygid) {
        slmm_assert( ! ed.me);
        ed.me = &n;
      }
      n.rank = nir[2*j+1];
      n.lid_on_rank = -1;
      n.rank_idx = -1;
    }
    slmm_assert(ed.me);
  }
}

typedef std::map<Int, std::set<Int> > Rank2Gids;

void get_rank2gids (const CslMpi& cm, Rank2Gids& rank2rmtgids,
                    Rank2Gids& rank2owngids) {
  const Int myrank = cm.p->rank();
  for (Int i = 0; i < cm.nelemd; ++i) {
    const auto& ed = cm.ed(i);
    for (const auto& n: ed.nbrs) {
      if (n.rank == myrank) continue;
      // I need this rmt gid's lid.
      rank2rmtgids[n.rank].insert(n.gid);
      // This rank needs this gid's lid.
      rank2owngids[n.rank].insert(ed.me->gid);
    }
  }
}

// Fill in nbrs.lid_on_rank, the lid on the remote rank corresponding to the gid
// I share but do not own.
void comm_lid_on_rank (CslMpi& cm, const Rank2Gids& rank2rmtgids,
                       const Rank2Gids& rank2owngids,
                       std::map<Int, Int>& gid2rmt_owning_lid) {
  const Int myrank = cm.p->rank();

  std::map<Int, Int> gid2mylid;
  for (Int i = 0; i < cm.nelemd; ++i)
    gid2mylid[cm.ed(i).me->gid] = i;
  
  // Set up to recv remote (gid, lid) lists.
  Int rn = 0;
  for (const auto& e: rank2rmtgids)
    rn += e.second.size();
  const Int nrecv = rank2rmtgids.size();
  std::vector<Int> recv(rn), recvptr(nrecv+1), recvrank(nrecv);
  std::vector<mpi::Request> reqs(nrecv);
  recvptr[0] = 0;
  Int ir = 0;
  rn = 0;
  for (const auto& e: rank2rmtgids) {
    const Int ne = e.second.size();
    const Int rank = e.first;
    slmm_assert(rank != myrank);
    recvrank[ir] = rank;
    mpi::irecv(*cm.p, recv.data() + rn, ne, rank, 42, &reqs[ir++]);
    rn += ne;
    recvptr[ir] = rn;
  }

  // Send my (gid, lid) lists.
  Int sn = 0, i = 0;
  for (const auto& e: rank2owngids)
    sn += e.second.size();
  std::vector<Int> send(sn);
  std::vector<mpi::Request> sendreqs(rank2owngids.size());
  sn = 0;
  for (const auto& e: rank2owngids) {
    // Iteration through a set gives increasing GID, which means the LID list is
    // consistent between communicating ranks.
    Int slot = sn, pgid = -1;
    for (auto gid: e.second) {
      slmm_assert(gid > pgid);
      pgid = gid;
      send[slot++] = gid2mylid[gid];
    }
    const Int ne = e.second.size();
    const Int rank = e.first;
    mpi::isend(*cm.p, send.data() + sn, ne, rank, 42, &sendreqs[i++]);
    sn += ne;
  }

  // Actually recv.
  const Int count = rank2owngids.size();
  for (Int c = 0; c < count; ++c) {
    int idx;
    mpi::waitany(count, reqs.data(), &idx);
    const Int* rmtlids = recv.data() + recvptr[idx];
    const auto& gids = rank2rmtgids.at(recvrank[idx]);
    slmm_assert(recvptr[idx+1] - recvptr[idx] == gids.size());
    Int i = 0;
    for (auto gid: gids)
      gid2rmt_owning_lid[gid] = rmtlids[i++];
  }

  // Fill lid_on_rank and mylid_with_comm.
  std::vector<Int> mylid_with_comm;
  for (Int i = 0; i < cm.nelemd; ++i) {
    auto& ed = cm.ed(i);
    bool has_comm = false;
    for (auto& n: ed.nbrs)
      if (n.rank == myrank) {
        const auto it = gid2mylid.find(n.gid);
        slmm_throw_if(it == gid2mylid.end(),
                      "comm_lid_on_rank: On rank " << myrank << ", gid " << n.gid
                      << " is not in gid2mylid.");
        n.lid_on_rank = it->second;
      } else {
        has_comm = true;
        const auto it = gid2rmt_owning_lid.find(n.gid);
        slmm_throw_if(it == gid2rmt_owning_lid.end(),
                      "comm_lid_on_rank: On rank " << myrank << ", gid " << n.gid
                      << ", which I believe to be owned by rank " << n.rank
                      << ", is not in gid2rmt_owning_lid.");
        n.lid_on_rank = it->second;
      }
    if (has_comm)
      mylid_with_comm.push_back(i);
  }
  cm.mylid_with_comm.reset_capacity(mylid_with_comm.size(), true);
  for (Int i = 0; i < cm.mylid_with_comm.n(); ++i)
    cm.mylid_with_comm(i) = mylid_with_comm[i];

  mpi::waitall(sendreqs.size(), sendreqs.data());
}

// Useful maps between a linear index space 1:K to a set of K unique
// integers. These obviate sorts and use of hash or binary maps during time
// stepping.
void set_idx2_maps (CslMpi& cm, const Rank2Gids& rank2rmtgids,
                    const std::map<Int, Int>& gid2rmt_owning_lid) {
  const Int myrank = cm.p->rank();
  std::map<Int, Int> ranks;
  Int i = 0;
  for (const auto& e: rank2rmtgids)
    if (ranks.find(e.first) == ranks.end())
      ranks[e.first] = i++;
  ranks[myrank] = i;

  cm.ranks.reset_capacity(i+1, true);
  cm.ranks(i) = myrank;
  for (const auto& e: ranks)
    cm.ranks(e.second) = e.first;
  cm.sendcount.reset_capacity(i, true);
  cm.sendreq.reset_capacity(i, true);
  cm.recvreq.reset_capacity(i, true);

  const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
  std::vector<std::map<Int, Int> > lor2idx(nrmtrank);
  std::vector<Int> nlid_on_rank(nrmtrank);
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    const auto& rmtgids = rank2rmtgids.at(cm.ranks(ri));
    nlid_on_rank[ri] = rmtgids.size();
  }
  cm.lid_on_rank.init(nrmtrank, nlid_on_rank.data());
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    const auto& rmtgids = rank2rmtgids.at(cm.ranks(ri));
    auto& lor2idxri = lor2idx[ri];
    Int idx = 0;
    for (const auto& gid: rmtgids) {
      const auto lid = gid2rmt_owning_lid.at(gid);
      lor2idxri[lid] = idx;
      cm.lid_on_rank(ri,idx) = lid;
      ++idx;
    }
    slmm_assert(idx == nlid_on_rank[ri]);
  }

  for (i = 0; i < cm.nelemd; ++i) {
    auto& ed = cm.ed(i);
    ed.me->rank_idx = ranks.at(ed.me->rank);
    ed.me->lid_on_rank_idx = i;
    for (auto& n: ed.nbrs) {
      n.rank_idx = ranks.at(n.rank);
      n.lid_on_rank_idx = n.rank == myrank ?
        i :
        lor2idx[n.rank_idx].at(n.lid_on_rank);
    }
    ed.src = CslMpi::IntArray2D("src", cm.nlev, cm.np2);
    ed.q_extrema = CslMpi::Array<Real**[2]>("q_extrema", cm.qsize, cm.nlev);
  }
}

// In the original MPI pattern that has been in HOMME for years, each owned cell
// has a 1-halo patch of bulk data. The allocations in this routine use
// essentially the same amount of memory, but not more. We could use less if we
// were willing to realloc space at each SL time step.
void alloc_mpi_buffers (CslMpi& cm, const Rank2Gids& rank2rmtgids,
                        const Rank2Gids& rank2owngids) {
  const auto myrank = cm.p->rank();
  // sizeof real, int, single int (b/c of alignment)
  const Int sor = sizeof(Real), soi = sizeof(Int), sosi = sor;
  static_assert(sizeof(Real) >= sizeof(Int),
                "For buffer packing, we require sizeof(Real) >= sizeof(Int)");
  const auto xbufcnt = [&] (const std::set<Int>& rmtgids,
                            const std::set<Int>& owngids) -> Int {
    return (sosi + (2*soi + (2*soi)*cm.nlev)*rmtgids.size() + // meta data
            owngids.size()*cm.nlev*cm.np2*3*sor);             // bulk data
  };
  const auto qbufcnt = [&] (const std::set<Int>& rmtgids,
                            const std::set<Int>& owngids) -> Int {
    return ((rmtgids.size()*2 +      // min/max q
             owngids.size()*cm.np2)* // q
            cm.qsize*cm.nlev*sor);
  };
  const auto bytes2real = [&] (const Int& bytes) {
    return (bytes + sor - 1)/sor;
  };

  slmm_assert(cm.ranks.back() == myrank);
  const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
  std::vector<Int> nlid_per_rank(nrmtrank), sendsz(nrmtrank), recvsz(nrmtrank);
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    const auto& rmtgids = rank2rmtgids.at(cm.ranks(ri));
    const auto& owngids = rank2owngids.at(cm.ranks(ri));
    nlid_per_rank[ri] = rmtgids.size();
    sendsz[ri] = bytes2real(std::max(xbufcnt(rmtgids, owngids),
                                     qbufcnt(owngids, rmtgids)));
    recvsz[ri] = bytes2real(std::max(xbufcnt(owngids, rmtgids),
                                     qbufcnt(rmtgids, owngids)));
  }
  cm.nx_in_rank.reset_capacity(nrmtrank, true);
  cm.nx_in_lid.init(nrmtrank, nlid_per_rank.data());
  cm.bla.init(nrmtrank, nlid_per_rank.data(), cm.nlev);
  cm.sendbuf.init(nrmtrank, sendsz.data());
  cm.recvbuf.init(nrmtrank, recvsz.data());
#ifdef HORIZ_OPENMP
  cm.ri_lidi_locks.init(nrmtrank, nlid_per_rank.data());
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    auto&& locks = cm.ri_lidi_locks(ri);
    for (auto& lock: locks)
      omp_init_lock(&lock);
  }
#endif
}

// At simulation initialization, set up a bunch of stuff to make the work at
// each step as small as possible.
void setup_comm_pattern (CslMpi& cm, const Int* nbr_id_rank, const Int* nirptr) {
  collect_gid_rank(cm, nbr_id_rank, nirptr);
  Rank2Gids rank2rmtgids, rank2owngids;
  get_rank2gids(cm, rank2rmtgids, rank2owngids);
  {
    std::map<Int, Int> gid2rmt_owning_lid;
    comm_lid_on_rank(cm, rank2rmtgids, rank2owngids, gid2rmt_owning_lid);
    set_idx2_maps(cm, rank2rmtgids, gid2rmt_owning_lid);
  }
  alloc_mpi_buffers(cm, rank2rmtgids, rank2owngids);
  if (cm.ed(0).me->rank == 0) std::cout << "COMPOSE> use SL MPI pattern\n";
}

inline int get_tid () {
#ifdef HORIZ_OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

inline int get_num_threads () {
#ifdef HORIZ_OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif
}

// mylid_with_comm(rankidx) is a list of element LIDs that have relations with
// other elements on other ranks. For horizontal threading, need to find the
// subsets that fit within the usual horizontal-threading nets:nete ranges.
void init_mylid_with_comm_threaded (CslMpi& cm, const Int& nets, const Int& nete) {
#ifdef HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
  {
    cm.mylid_with_comm_tid_ptr.reset_capacity(get_num_threads()+1, true);
    cm.horiz_openmp = get_num_threads() > 1;
  }
#ifdef HORIZ_OPENMP
# pragma omp barrier
#endif
  const int tid = get_tid();
  const auto& beg = std::lower_bound(cm.mylid_with_comm.begin(),
                                     cm.mylid_with_comm.end(), nets);
  slmm_assert(cm.p->size() == 1 || beg != cm.mylid_with_comm.end());
  cm.mylid_with_comm_tid_ptr(tid) = beg - cm.mylid_with_comm.begin();
  if (tid == cm.mylid_with_comm_tid_ptr.n() - 2) {
    const auto& end = std::lower_bound(cm.mylid_with_comm.begin(),
                                       cm.mylid_with_comm.end(), nete+1);
    cm.mylid_with_comm_tid_ptr(tid+1) = end - cm.mylid_with_comm.begin();
  }
#ifdef HORIZ_OPENMP
# pragma omp barrier
#endif
}

CslMpi::Ptr init (const mpi::Parallel::Ptr& p,
                  Int np, Int nlev, Int qsize, Int qsized, Int nelemd,
                  const Int* nbr_id_rank, const Int* nirptr) {
  auto cm = std::make_shared<CslMpi>(p, np, nlev, qsize, qsized, nelemd);
  setup_comm_pattern(*cm, nbr_id_rank, nirptr);
  return cm;
}

// Set pointers to HOMME data arrays.
void set_elem_data (CslMpi& cm, const Int ie, const Real* metdet, const Real* qdp,
                    const Real* dp, Real* q, const Int nelem_in_patch) {
  slmm_assert(ie < cm.ed.size());
  slmm_assert(cm.ed(ie).nbrs.size() == nelem_in_patch);
  auto& e = cm.ed(ie);
  e.metdet = metdet;
  e.qdp = qdp;
  e.dp = dp;
  e.q = q;
}

void setup_irecv (CslMpi& cm) {
#ifdef HORIZ_OPENMP
# pragma omp master
#endif
  {
    const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
    for (Int ri = 0; ri < nrmtrank; ++ri) {
      auto&& recvbuf = cm.recvbuf(ri);
      // The count is just the number of slots available, which can be larger
      // than what is actually being received.
      mpi::irecv(*cm.p, recvbuf.data(), recvbuf.n(), cm.ranks(ri), 42,
                 &cm.recvreq(ri));
    }
  }  
}

void isend (CslMpi& cm, const bool want_req = true) {
#ifdef HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
  {
    const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
    for (Int ri = 0; ri < nrmtrank; ++ri)
      mpi::isend(*cm.p, cm.sendbuf(ri).data(), cm.sendcount(ri),
                 cm.ranks(ri), 42, want_req ? &cm.sendreq(ri) : nullptr);
  }
}

void recv_and_wait_on_send (CslMpi& cm) {
#ifdef HORIZ_OPENMP
# pragma omp master
#endif
  {
    mpi::waitall(cm.sendreq.n(), cm.sendreq.data());
    mpi::waitall(cm.recvreq.n(), cm.recvreq.data());
  }
#ifdef HORIZ_OPENMP
# pragma omp barrier
#endif
}

void recv (CslMpi& cm) {
#ifdef HORIZ_OPENMP
# pragma omp master
#endif
  {
    mpi::waitall(cm.recvreq.n(), cm.recvreq.data());
  }
#ifdef HORIZ_OPENMP
# pragma omp barrier
#endif
}

// Find where each departure point is.
void analyze_dep_points (CslMpi& cm, const Int& nets, const Int& nete,
                         const FA4<Real>& dep_points) {
  const auto myrank = cm.p->rank();
  const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
  cm.bla.zero();
#ifdef HORIZ_OPENMP
# pragma omp for
#endif
  for (Int ri = 0; ri < nrmtrank; ++ri)
    cm.nx_in_lid(ri).zero();
  for (Int tci = nets; tci <= nete; ++tci) {
    const auto& mesh = g_advecter->local_mesh(tci);
    const auto tgt_idx = g_advecter->local_mesh_tgt_elem(tci);
    auto& ed = cm.ed(tci);
    ed.own.clear();
    for (Int lev = 0; lev < cm.nlev; ++lev)
      for (Int k = 0; k < cm.np2; ++k) {
        Int sci = slmm::get_src_cell(mesh, &dep_points(0,k,lev,tci), tgt_idx);
        if (sci == -1 && g_advecter->nearest_point_permitted(lev))
          sci = slmm::get_nearest_point(
            mesh, g_advecter->nearest_point_data(tci),
            &dep_points(0,k,lev,tci), tgt_idx);
        if (sci == -1) {
          std::stringstream ss;
          ss.precision(17);
          const auto* v = &dep_points(0,k,lev,tci);
          ss << "Departure point is outside of halo:\n"
             << "  nearest point permitted: "
             << g_advecter->nearest_point_permitted(lev)
             << "\n  elem LID " << tci
             << " elem GID " << ed.me->gid
             << " (lev, k) (" << lev << ", " << k << ")"
             << " v " << v[0] << " " << v[1] << " " << v[2]
             << "\n  tgt_idx " << tgt_idx
             << " local mesh:\n  " << slmm::to_string(mesh) << "\n";
          slmm_throw_if(sci == -1, ss.str());
        }
        ed.src(lev,k) = sci;
        if (ed.nbrs(sci).rank == myrank) {
          ed.own.inc();
          auto& t = ed.own.back();
          t.lev = lev; t.k = k;
        } else {
          const auto ri = ed.nbrs(sci).rank_idx;
          const auto lidi = ed.nbrs(sci).lid_on_rank_idx;
#ifdef HORIZ_OPENMP
          omp_lock_t* lock;
          if (cm.horiz_openmp) {
            lock = &cm.ri_lidi_locks(ri,lidi);
            omp_set_lock(lock);
          }
#endif
          {
            ++cm.nx_in_lid(ri,lidi);
            ++cm.bla(ri,lidi,lev).xptr;
          }
#ifdef HORIZ_OPENMP
          if (cm.horiz_openmp) omp_unset_lock(lock);
#endif
        }
      }
  }
#ifdef HORIZ_OPENMP
# pragma omp barrier
# pragma omp for
#endif
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    auto& nx_in_rank = cm.nx_in_rank(ri);
    nx_in_rank = 0;
    for (Int i = 0, n = cm.lid_on_rank(ri).n(); i < n; ++i)
      nx_in_rank += cm.nx_in_lid(ri,i);
  }
}

static const int nreal_per_2int = (2*sizeof(Int) + sizeof(Real) - 1) / sizeof(Real);

template <typename Buffer>
Int setbuf (Buffer& buf, const Int& os, const Int& i1, const Int& i2) {
  Int* const b = reinterpret_cast<Int*>(&buf(os));
  b[0] = i1;
  b[1] = i2;
  return nreal_per_2int;
}

template <typename Buffer>
Int getbuf (Buffer& buf, const Int& os, Int& i1, Int& i2) {
  const Int* const b = reinterpret_cast<const Int*>(&buf(os));
  i1 = b[0];
  i2 = b[1];
  return nreal_per_2int;
}

/* Pack the departure points (x). We use two passes. We also set up the q
   metadata. Two passes lets us do some efficient tricks that are not available
   with one pass. Departure point and q messages are formatted as follows:
    xs: (#x-in-rank    int
         pad           i
         (lid          i     only packed if #x in lid > 0
          #x-in-lid    i     > 0
          (lev         i     only packed if #x in (lid,lev) > 0
           #x          i     > 0
           x         3 real
            *#x) *#lev) *#lid) *#rank
    qs: (q-extrema    2 qsize r    (min, max) packed together
         q              qsize r
          *#x) *#lev *#lid *#rank
 */
void pack_dep_points_sendbuf_pass1 (CslMpi& cm) {
  const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
#ifdef HORIZ_OPENMP
# pragma omp for
#endif
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    auto&& sendbuf = cm.sendbuf(ri);
    const auto&& lid_on_rank = cm.lid_on_rank(ri);
    Int xos = 0, qos = 0;
    xos += setbuf(sendbuf, xos, cm.nx_in_rank(ri), 0 /* empty space for alignment */);
    auto&& bla = cm.bla(ri);
    for (Int lidi = 0, lidn = cm.lid_on_rank(ri).n(); lidi < lidn; ++lidi) {
      auto nx_in_lid = cm.nx_in_lid(ri,lidi);
      if (nx_in_lid == 0) continue;
      xos += setbuf(sendbuf, xos, lid_on_rank(lidi), nx_in_lid);
      for (Int lev = 0; lev < cm.nlev; ++lev) {
        auto& t = bla(lidi,lev);
        t.qptr = qos;
        slmm_assert_high(t.cnt == 0);
        const Int nx = t.xptr;
        if (nx == 0) {
          t.xptr = -1;
          continue;
        }
        slmm_assert_high(nx > 0);
        const auto dos = setbuf(sendbuf, xos, lev, nx);
        t.xptr = xos + dos;
        xos += dos + 3*nx;
        qos += 2 + nx;
        nx_in_lid -= nx;
      }
      slmm_assert(nx_in_lid == 0);
    }
    cm.sendcount(ri) = xos;
  }
}

void pack_dep_points_sendbuf_pass2 (CslMpi& cm, const FA4<const Real>& dep_points) {
  const auto myrank = cm.p->rank();
  const int tid = get_tid();
  for (Int ptr = cm.mylid_with_comm_tid_ptr(tid),
           end = cm.mylid_with_comm_tid_ptr(tid+1);
       ptr < end; ++ptr) {
    const Int tci = cm.mylid_with_comm(ptr);
    auto& ed = cm.ed(tci);
    ed.rmt.clear();
    for (Int lev = 0; lev < cm.nlev; ++lev) {
      for (Int k = 0; k < cm.np2; ++k) {
        const Int sci = ed.src(lev,k);
        const auto& nbr = ed.nbrs(sci);
        if (nbr.rank == myrank) continue;
        const Int ri = nbr.rank_idx;
        const Int lidi = nbr.lid_on_rank_idx;
        auto&& sb = cm.sendbuf(ri);
#ifdef HORIZ_OPENMP
        omp_lock_t* lock;
        if (cm.horiz_openmp) {
          lock = &cm.ri_lidi_locks(ri,lidi);
          omp_set_lock(lock);
        }
#endif
        Int xptr, qptr, cnt; {
          auto& t = cm.bla(ri,lidi,lev);
          qptr = t.qptr;
          cnt = t.cnt;
          xptr = t.xptr + 3*cnt;
          ++t.cnt;
        }
#ifdef HORIZ_OPENMP
        if (cm.horiz_openmp) omp_unset_lock(lock);
#endif
        slmm_assert_high(xptr > 0);
        for (Int i = 0; i < 3; ++i)
          sb(xptr + i) = dep_points(i,k,lev,tci);
        ed.rmt.inc();
        auto& item = ed.rmt.back();
        item.q_extrema_ptr = cm.qsize * qptr;
        item.q_ptr = item.q_extrema_ptr + cm.qsize*(2 + cnt);
        item.lev = lev;
        item.k = k;
      }
    }
  }
}

template <Int np>
void calc_q_extrema (CslMpi& cm, const Int& nets, const Int& nete) {
  for (Int tci = nets; tci <= nete; ++tci) {
    auto& ed = cm.ed(tci);
    const FA2<const Real> dp(ed.dp, cm.np2, cm.nlev);
    const FA3<const Real> qdp(ed.qdp, cm.np2, cm.nlev, cm.qsize);
    const FA3<Real> q(ed.q, cm.np2, cm.nlev, cm.qsize);
    for (Int iq = 0; iq < cm.qsize; ++iq)
      for (Int lev = 0; lev < cm.nlev; ++lev) {
        const Real* const dp0 = &dp(0,lev);
        const Real* const qdp0 = &qdp(0,lev,iq);
        Real* const q0 = &q(0,lev,iq);
        Real q_min_s, q_max_s;
        q0[0] = qdp0[0] / dp0[0];
        q_min_s = q_max_s = q0[0];
        for (Int k = 1; k < np*np; ++k) {
          q0[k] = qdp0[k] / dp0[k];
          q_min_s = std::min(q_min_s, q0[k]);
          q_max_s = std::max(q_max_s, q0[k]);
        }
        ed.q_extrema(iq,lev,0) = q_min_s;
        ed.q_extrema(iq,lev,1) = q_max_s;
      }
  }  
}

template <Int np>
void calc_q (const CslMpi& cm, const Int& src_lid, const Int& lev,
             const Real* const dep_point, Real* const q_tgt, const bool use_q) {
  const auto& m = g_advecter->local_mesh(src_lid);
  const auto my_idx = g_advecter->local_mesh_tgt_elem(src_lid);

  Real ref_coord[2];
  siqk::sqr::calc_sphere_to_ref(m.p, slmm::slice(m.e, my_idx), dep_point,
                                ref_coord[0], ref_coord[1]);

  // Interpolate.
  Real rx[4], ry[4];
  switch (g_advecter->alg()) {
  case slmm::Advecter::Alg::csl_gll:
    slmm::gll_np4_eval(ref_coord[0], rx);
    slmm::gll_np4_eval(ref_coord[1], ry);
    break;
  case slmm::Advecter::Alg::csl_gll_subgrid:
    slmm::gll_np4_subgrid_eval(ref_coord[0], rx);
    slmm::gll_np4_subgrid_eval(ref_coord[1], ry);
    break;
  case slmm::Advecter::Alg::csl_gll_exp:
    slmm::gll_np4_subgrid_exp_eval(ref_coord[0], rx);
    slmm::gll_np4_subgrid_exp_eval(ref_coord[1], ry);
    break;
  default:
    slmm_assert(0);
  }

  const auto& ed = cm.ed(src_lid);
  const Int levos = np*np*lev;
  const Int np2nlev = np*np*cm.nlev;
  if (use_q) {
    // We can use q from calc_q_extrema.
    const Real* const qs0 = ed.q + levos;
#   pragma ivdep
    for (Int iq = 0; iq < cm.qsize; ++iq) {
      const Real* const qs = qs0 + iq*np2nlev;
      q_tgt[iq] =
        (ry[0]*(rx[0]*qs[ 0] + rx[1]*qs[ 1] + rx[2]*qs[ 2] + rx[3]*qs[ 3]) +
         ry[1]*(rx[0]*qs[ 4] + rx[1]*qs[ 5] + rx[2]*qs[ 6] + rx[3]*qs[ 7]) +
         ry[2]*(rx[0]*qs[ 8] + rx[1]*qs[ 9] + rx[2]*qs[10] + rx[3]*qs[11]) +
         ry[3]*(rx[0]*qs[12] + rx[1]*qs[13] + rx[2]*qs[14] + rx[3]*qs[15]));
    }
  } else {
    // q from calc_q_extrema is being overwritten, so have to use qdp/dp.
    const Real* const dp = ed.dp + levos;
    const Real* const qdp0 = ed.qdp + levos;
#   pragma ivdep
    for (Int iq = 0; iq < cm.qsize; ++iq) {
      const Real* const qdp = qdp0 + iq*np2nlev;
      q_tgt[iq] = (ry[0]*(rx[0]*(qdp[ 0]/dp[ 0]) + rx[1]*(qdp[ 1]/dp[ 1])  +
                          rx[2]*(qdp[ 2]/dp[ 2]) + rx[3]*(qdp[ 3]/dp[ 3])) +
                   ry[1]*(rx[0]*(qdp[ 4]/dp[ 4]) + rx[1]*(qdp[ 5]/dp[ 5])  +
                          rx[2]*(qdp[ 6]/dp[ 6]) + rx[3]*(qdp[ 7]/dp[ 7])) +
                   ry[2]*(rx[0]*(qdp[ 8]/dp[ 8]) + rx[1]*(qdp[ 9]/dp[ 9])  +
                          rx[2]*(qdp[10]/dp[10]) + rx[3]*(qdp[11]/dp[11])) +
                   ry[3]*(rx[0]*(qdp[12]/dp[12]) + rx[1]*(qdp[13]/dp[13])  +
                          rx[2]*(qdp[14]/dp[14]) + rx[3]*(qdp[15]/dp[15])));
    }
  }
}

template <Int np>
void calc_rmt_q (CslMpi& cm) {
  const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
#ifdef HORIZ_OPENMP
# pragma omp for
#endif
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    const auto&& xs = cm.recvbuf(ri);
    auto&& qs = cm.sendbuf(ri);
    Int xos = 0, qos = 0, nx_in_rank, padding;
    xos += getbuf(xs, xos, nx_in_rank, padding);
    if (nx_in_rank == 0) continue;
    // The upper bound is to prevent an inf loop if the msg is corrupted.
    for (Int lidi = 0; lidi < cm.nelemd; ++lidi) {
      Int lid, nx_in_lid;
      xos += getbuf(xs, xos, lid, nx_in_lid);
      const auto& ed = cm.ed(lid);
      for (Int levi = 0; levi < cm.nlev; ++levi) { // same re: inf loop
        Int lev, nx;
        xos += getbuf(xs, xos, lev, nx);
        slmm_assert(nx > 0);
        for (Int iq = 0; iq < cm.qsize; ++iq)
          for (int i = 0; i < 2; ++i)
            qs(qos + 2*iq + i) = ed.q_extrema(iq, lev, i);
        qos += 2*cm.qsize;
        for (Int ix = 0; ix < nx; ++ix) {
          calc_q<np>(cm, lid, lev, &xs(xos), &qs(qos), true);
          xos += 3;
          qos += cm.qsize;
        }
        nx_in_lid -= nx;
        nx_in_rank -= nx;
        if (nx_in_lid == 0) break;
      }
      slmm_assert(nx_in_lid == 0);
      if (nx_in_rank == 0) break;
    }
    slmm_assert(nx_in_rank == 0);
    cm.sendcount(ri) = qos;
  }
}

template <Int np>
void calc_own_q (CslMpi& cm, const Int& nets, const Int& nete,
                 const FA4<const Real>& dep_points,
                 const FA4<Real>& q_min, const FA4<Real>& q_max) {
  for (Int tci = nets; tci <= nete; ++tci) {
    const Int ie0 = tci - nets;
    auto& ed = cm.ed(tci);
    const FA3<Real> q_tgt(ed.q, cm.np2, cm.nlev, cm.qsize);
    for (const auto& e: ed.own) {
      const Int slid = ed.nbrs(ed.src(e.lev, e.k)).lid_on_rank;
      const auto& sed = cm.ed(slid);
      for (Int iq = 0; iq < cm.qsize; ++iq) {
        q_min(e.k, e.lev, iq, ie0) = sed.q_extrema(iq, e.lev, 0);
        q_max(e.k, e.lev, iq, ie0) = sed.q_extrema(iq, e.lev, 1);
      }
      Real qtmp[QSIZE_D];
      calc_q<np>(cm, slid, e.lev, &dep_points(0, e.k, e.lev, tci), qtmp, false);
      for (Int iq = 0; iq < cm.qsize; ++iq)
        q_tgt(e.k, e.lev, iq) = qtmp[iq];
    }
  }
}

void copy_q (CslMpi& cm, const Int& nets,
             const FA4<Real>& q_min, const FA4<Real>& q_max) {
  const auto myrank = cm.p->rank();
  const int tid = get_tid();
  for (Int ptr = cm.mylid_with_comm_tid_ptr(tid),
           end = cm.mylid_with_comm_tid_ptr(tid+1);
       ptr < end; ++ptr) {
    const Int tci = cm.mylid_with_comm(ptr);
    const Int ie0 = tci - nets;
    auto& ed = cm.ed(tci);
    const FA3<Real> q_tgt(ed.q, cm.np2, cm.nlev, cm.qsize);
    for (const auto& e: ed.rmt) {
      slmm_assert(ed.nbrs(ed.src(e.lev, e.k)).rank != myrank);
      const Int ri = ed.nbrs(ed.src(e.lev, e.k)).rank_idx;
      const auto&& recvbuf = cm.recvbuf(ri);
      for (Int iq = 0; iq < cm.qsize; ++iq) {
        q_min(e.k, e.lev, iq, ie0) = recvbuf(e.q_extrema_ptr + 2*iq    );
        q_max(e.k, e.lev, iq, ie0) = recvbuf(e.q_extrema_ptr + 2*iq + 1);
      }
      for (Int iq = 0; iq < cm.qsize; ++iq) {
        slmm_assert(recvbuf(e.q_ptr + iq) != -1);
        q_tgt(e.k, e.lev, iq) = recvbuf(e.q_ptr + iq);
      }
    }
  }
}

/* dep_points is const in principle, but if lev <=
   semi_lagrange_nearest_point_lev, a departure point may be altered if the
   winds take it outside of the comm halo.
 */
template <int np>
void step (
  CslMpi& cm, const Int nets, const Int nete,
  Cartesian3D* dep_points_r,    // dep_points(1:3, 1:np, 1:np)
  Real* q_min_r, Real* q_max_r) // q_{min,max}(1:np, 1:np, lev, 1:qsize, ie-nets+1)
{
  static_assert(np == 4, "SLMM CSL with special MPI is supported for np 4 only.");
  slmm_assert(cm.np == 4);

  const FA4<Real>
    dep_points(reinterpret_cast<Real*>(dep_points_r),
               3, cm.np2, cm.nlev, cm.nelemd);
  const Int nelem = nete - nets + 1;
  const FA4<Real>
    q_min(q_min_r, cm.np2, cm.nlev, cm.qsize, nelem),
    q_max(q_max_r, cm.np2, cm.nlev, cm.qsize, nelem);

  // Partition my elements that communicate with remotes among threads, if I
  // haven't done that yet.
  if (cm.mylid_with_comm_tid_ptr.n() == 0)
    init_mylid_with_comm_threaded(cm, nets, nete);
  // Set up to receive departure point requests from remotes.
  setup_irecv(cm);
  // Determine where my departure points are, and set up requests to remotes as
  // well as to myself to fulfill these.
  analyze_dep_points(cm, nets, nete, dep_points);
  pack_dep_points_sendbuf_pass1(cm);
  pack_dep_points_sendbuf_pass2(cm, dep_points);
  // Send requests.
  isend(cm);
  // While waiting, compute q extrema in each of my elements.
  calc_q_extrema<np>(cm, nets, nete);
  // Wait for the departure point requests. Since this requires a thread
  // barrier, at the same time make sure the send buffer is free for use.
  recv_and_wait_on_send(cm);
  // Compute the requested q for departure points from remotes.
  calc_rmt_q<np>(cm);
  // Send q data.
  isend(cm, false);
  // Set up to receive q for each of my departure point requests sent to
  // remotes. We can't do this until the OpenMP barrier in isend assures that
  // all threads are done with the receive buffer's departure points.
  setup_irecv(cm);
  // While waiting to get my data from remotes, compute q for departure points
  // that have remained in my elements.
  calc_own_q<np>(cm, nets, nete, dep_points, q_min, q_max);
  // Receive remote q data and use this to fill in the rest of my fields.
  recv(cm);
  copy_q(cm, nets, q_min, q_max);
  // Don't need to wait on send buffer again because MPI-level synchronization
  // outside of SL transport assures the send buffer is ready at the next call
  // to step. But do need to dealloc the send requests.
}
} // namespace cslmpi
} // namespace homme

// Valid after slmm_init_local_mesh_ is called.
int slmm_unittest () {
  int nerr = 0, ne;
  {
    ne = 0;
    for (int i = 0; i < homme::g_advecter->nelem(); ++i)
      ne += slmm::unittest(homme::g_advecter->local_mesh(i),
                           homme::g_advecter->local_mesh_tgt_elem(i));
    if (ne)
      fprintf(stderr, "slmm_unittest: slmm::unittest returned %d\n", ne);
    nerr += ne;
  }
  return nerr;
}

#include <cstdlib>

struct Experiment {
  int sl_mpi;
};

template <typename T> T strto(const char* s);
template <> inline int strto (const char* s) { return std::atoi(s); }
template <> inline bool strto (const char* s) { return std::atoi(s); }
template <> inline double strto (const char* s) { return std::atof(s); }
template <> inline std::string strto (const char* s) { return std::string(s); }

Experiment get_options () {
  Experiment e;
  e.sl_mpi = 1;
  char* var_s = std::getenv("slmpi");
  if (var_s) e.sl_mpi = strto<int>(var_s);
  return e;
}

static homme::cslmpi::CslMpi::Ptr g_csl_mpi;

extern "C" {
void slmm_init_impl (
  homme::Int fcomm, homme::Int transport_alg, homme::Int np,
  homme::Int nlev, homme::Int qsize, homme::Int qsized, homme::Int nelemd,
  const homme::Int** nbr_id_rank, const homme::Int** nirptr,
  homme::Int sl_nearest_point_lev)
{
  auto e = get_options();
  homme::slmm_init(np, nelemd, transport_alg, sl_nearest_point_lev - 1);
  if (e.sl_mpi && homme::g_advecter->is_cisl())
    e.sl_mpi = 0;
  if (e.sl_mpi) {
    const auto p = homme::mpi::make_parallel(MPI_Comm_f2c(fcomm));
    g_csl_mpi = homme::cslmpi::init(p, np, nlev, qsize, qsized, nelemd,
                                    *nbr_id_rank, *nirptr);
  }
}

void slmm_get_mpi_pattern (homme::Int* sl_mpi) {
  *sl_mpi = g_csl_mpi ? 1 : 0;
}

void slmm_init_local_mesh_ (
  homme::Int* ie, homme::Cartesian3D** neigh_corners, homme::Int* nnc,
  homme::Cartesian3D* p_inside)
{
  homme::g_advecter->init_local_mesh_if_needed(
    *ie - 1, homme::FA3<const homme::Real>(
      reinterpret_cast<const homme::Real*>(*neigh_corners), 3, 4, *nnc),
    reinterpret_cast<const homme::Real*>(p_inside));
  if (*ie == homme::g_advecter->nelem())
    homme::g_advecter->init_M_tgt_if_needed();
}

void slmm_check_ref2sphere_ (homme::Int* ie, homme::Cartesian3D* p) {
  homme::g_advecter->check_ref2sphere(
    *ie - 1, reinterpret_cast<const homme::Real*>(p));
}

void slmm_advect_ (
  homme::Int* lev, homme::Int* ie, homme::Int* nnc, homme::Int* np, homme::Int* nlev,
  homme::Int* qsize, homme::Int* nets, homme::Int* nete,
  homme::Cartesian3D* dep_points, homme::Real* neigh_q, homme::Real* J_t,
  homme::Real* dp3d, homme::Int* tl_np1, homme::Real* q, homme::Real* minq,
  homme::Real* maxq)
{
  if (homme::g_advecter->is_cisl()) {
#ifdef BUILD_CISL
    homme::slmm_project(*lev - 1, *nets - 1, *ie - 1, *nnc, *np, *nlev, *qsize,
                        dep_points, neigh_q, J_t, dp3d, *tl_np1 - 1, q, minq, maxq);
#else
    throw std::runtime_error("Closed for business to speed up compilation.");
#endif
  } else {
    slmm_assert(np == 4);
    homme::slmm_csl<4>(*lev - 1, *nets - 1, *ie - 1, *nnc, *nlev, *qsize,
                       dep_points, neigh_q, J_t, dp3d, *tl_np1 - 1, q, minq, maxq);
  }
}

void slmm_csl_set_elem_data (
  homme::Int ie, homme::Real* metdet, homme::Real* qdp, homme::Real* dp,
  homme::Real* q, homme::Int nelem_in_patch)
{
  slmm_assert(g_csl_mpi);
  homme::cslmpi::set_elem_data(*g_csl_mpi, ie - 1, metdet, qdp, dp, q,
                               nelem_in_patch);
}

void slmm_csl_ (
  homme::Int* nets, homme::Int* nete, homme::Cartesian3D* dep_points,
  homme::Real* minq, homme::Real* maxq, homme::Int* info)
{
  slmm_assert(g_csl_mpi);
  *info = 0;
  try {
    homme::cslmpi::step<4>(*g_csl_mpi, *nets - 1, *nete - 1,
                           dep_points, minq, maxq);
  } catch (const std::exception& e) {
    std::cerr << e.what();
    *info = -1;
  }
}
} // extern "C"
