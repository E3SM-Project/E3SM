#include <mpi.h>

// Uncomment this to look for MPI-related memory leaks.
//#define COMPOSE_DEBUG_MPI

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

#ifdef _OPENMP
# include <omp.h>
#endif

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

// Both cubed_sphere_map=0 and cubed_sphere_map=2 can use this
// method. (cubed_sphere_map=1 is not impl'ed in Homme.) This method is natural
// for cubed_sphere_map=2, so RRM works automatically. cubed_sphere_map=0 is
// supported in Homme only for regular cubed-sphere meshes, not RRM; in that
// case, edges are great arcs, so again this impl works. In contrast,
// calc_sphere_to_ref has to be specialized on cubed_sphere_map.
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

std::string to_string (const siqk::sh::Mesh<siqk::ko::HostSpace>& m) {
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

namespace nearest_point {
/* Get external segments in preproc step.
   Get approximate nearest point in each segment.
     Project onto plane of segment and normalize.
     If outside of arc, choose nearest endpoint.
     Approx distance as cartesian distance.
   Of the approx dists, take the point associated with smallest.
 */

typedef siqk::sh::Mesh<siqk::ko::HostSpace> Mesh;

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

void init_nearest_point_data (const nearest_point::Mesh& m,
                              MeshNearestPointData& d) {
  nearest_point::fill_perim(m, d);
}

int get_nearest_point (const nearest_point::Mesh& m, const MeshNearestPointData& d,
                       Real* v, const Int my_ic) {
  nearest_point::calc(m, d, v);
  return get_src_cell(m, v, my_ic);
}

Int unittest (const nearest_point::Mesh& m, const Int tgt_elem) {
  Int nerr = 0, ne = 0;
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
          v[d] = (  b*(a*m.p(cell[0], d) + oma*m.p(cell[1], d)) + 
                  omb*(a*m.p(cell[3], d) + oma*m.p(cell[2], d)));
        if (i < 2 && j < 2) {
          if (get_src_cell(m, v, tgt_elem) != ic) ++ne;
        } else {
          if (get_src_cell(m, v, tgt_elem) == -1) ++ne;
        }
      }
  }
  if (ne) pr("slmm::unittest: get_src_cell failed");
  nerr += ne;
  ne = nearest_point::test_canpoa();
  nerr += ne;
  if (ne) pr("slmm::unittest: test_canpoa failed");
  {
    MeshNearestPointData d;
    init_nearest_point_data(m, d);
    ne = nearest_point::test_fill_perim(m, tgt_elem, d);
    if (ne) pr("slmm::unittest: test_fill_perim failed");
    nerr += ne;
    ne = nearest_point::test_calc(m, tgt_elem, d);
    if (ne) pr("slmm::unittest: test_calc failed");
    nerr += ne;
  }
  return nerr;
}

typedef Kokkos::View<Int*, ko::HostSpace> Ints;

//TODO We might switch to one 
// Local mesh patch centered on the target element.
struct LocalMesh : public siqk::sh::Mesh<ko::HostSpace> {
  typedef siqk::sh::Mesh<ko::HostSpace> Super;
  typedef typename Super::IntArray IntArray;
  typedef typename Super::RealArray RealArray;

  // tgt_elem is the index of the target element in this mesh.
  Int tgt_elem;
};

// Wrap call to siqk::sqr::calc_sphere_to_ref. That impl supports only
// cubed_sphere_map=2, and we want to keep it that way. This wrapper supports,
// in addition, cubed_sphere_map=0.
struct SphereToRef {
  void init (const Int cubed_sphere_map, const Int nelem_global,
             const Ints::const_type& lid2facenum) {
    cubed_sphere_map_ = cubed_sphere_map;
    ne_ = static_cast<Int>(std::round(std::sqrt((nelem_global / 6))));
    slmm_throw_if( ! (cubed_sphere_map_ != 0 || 6*ne_*ne_ == nelem_global),
                  "If cubed_sphere_map = 0, then the mesh must be a "
                  "regular cubed-sphere.");
    lid2facenum_ = lid2facenum;
  }

  Real tol () const {
    return 1e3 * ne_ * std::numeric_limits<Real>::epsilon();
  }

  // See siqk::sqr::calc_sphere_to_ref for docs.
  void calc_sphere_to_ref (
    const Int& ie, const LocalMesh& m,
    const Real q[3],
    Real& a, Real& b,
    siqk::sqr::Info* const info = nullptr,
    const Int max_its = 10,
    const Real tol = 1e2*std::numeric_limits<Real>::epsilon()) const
  {
    if (cubed_sphere_map_ == 2)
      siqk::sqr::calc_sphere_to_ref(m.p, slice(m.e, m.tgt_elem), q, a, b,
                                    info, max_its, tol);
    else {
      const Int face = lid2facenum_(ie); //assume: ie corresponds to m.tgt_elem.
      map_sphere_coord_to_face_coord(face-1, q[0], q[1], q[2], a, b);
      a = map_face_coord_to_cell_ref_coord(a);
      b = map_face_coord_to_cell_ref_coord(b);
      if (info) { info->success = true; info->n_iterations = 1; }
    }
  }

  bool check (const Int& ie, const LocalMesh& m) const {
    if (cubed_sphere_map_ != 0) return true;
    const Int face = lid2facenum_(ie); //assume: ie corresponds to m.tgt_elem.
    Real cent[3] = {0};
    const auto cell = slice(m.e, m.tgt_elem);
    for (int i = 0; i < 4; ++i)
      for (int d = 0; d < 3; ++d)
        cent[d] += 0.25*m.p(cell[i], d);
    const Int cf = get_cube_face_idx(cent[0], cent[1], cent[2]) + 1;
    return face == cf;
  }

private:
  Int ne_, cubed_sphere_map_;
  Ints::const_type lid2facenum_;

  // Follow the description given in
  //     coordinate_systems_mod::unit_face_based_cube_to_unit_sphere.
  static Int get_cube_face_idx (const Real& x, const Real& y, const Real& z) {
    const Real ax = std::abs(x), ay = std::abs(y), az = std::abs(z);
    if (ax >= ay) {
      if (ax >= az) return x > 0 ? 0 : 2;
      else return z > 0 ? 5 : 4;
    } else {
      if (ay >= az) return y > 0 ? 1 : 3;
      else return z > 0 ? 5 : 4;
    }
  }

  static void map_sphere_coord_to_face_coord (
    const Int& face_idx, const Real& x, const Real& y, const Real& z,
    Real& fx, Real& fy)
  {
    static constexpr Real theta_max = 0.25*M_PI;
    Real d;
    switch (face_idx) {
    case  0: d = std::abs(x); fx =  y/d; fy =  z/d; break;
    case  1: d = std::abs(y); fx = -x/d; fy =  z/d; break;
    case  2: d = std::abs(x); fx = -y/d; fy =  z/d; break;
    case  3: d = std::abs(y); fx =  x/d; fy =  z/d; break;
    case  4: d = std::abs(z); fx =  y/d; fy =  x/d; break;
    default: d = std::abs(z); fx =  y/d; fy = -x/d;
    }
    fx = std::atan(fx) / theta_max;
    fy = std::atan(fy) / theta_max;
  }

  Real map_face_coord_to_cell_ref_coord (Real a) const {
    a = (0.5*(1 + a))*ne_;
    a = 2*(a - std::floor(a)) - 1;
    return a;
  }
};

// Advecter has purely mesh-local knowledge, with once exception noted below.
struct Advecter {
  typedef std::shared_ptr<Advecter> Ptr;
  typedef std::shared_ptr<const Advecter> ConstPtr;

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
            const Int cubed_sphere_map, const Int nearest_point_permitted_lev_bdy)
    : alg_(Alg::convert(transport_alg)),
      np_(np), np2_(np*np), np4_(np2_*np2_),
      cubed_sphere_map_(cubed_sphere_map),
      tq_order_(alg_ == Alg::qos ? 14 : 12),
      nearest_point_permitted_lev_bdy_(nearest_point_permitted_lev_bdy)
  {
    slmm_throw_if(cubed_sphere_map == 0 && Alg::is_cisl(alg_),
                  "When cubed_sphere_map = 0, SLMM supports only ISL methods.");
    local_mesh_.resize(nelem);
    if (Alg::is_cisl(alg_))
      mass_mix_.resize(np4_);
    if (nearest_point_permitted_lev_bdy_ >= 0)
      local_mesh_nearest_point_data_.resize(nelem);
  }

  void fill_nearest_points_if_needed();

  Int np  () const { return np_ ; }
  Int np2 () const { return np2_; }
  Int np4 () const { return np4_; }
  Int nelem () const { return local_mesh_.size(); }
  Int tq_order () const { return tq_order_; }
  Alg::Enum alg () const { return alg_; }
  bool is_cisl () const { return Alg::is_cisl(alg_); }

  Int cubed_sphere_map () const { return cubed_sphere_map_; }
  const Ints& lid2facenum () const { return lid2facenum_; }

  // nelem_global is used only if cubed_sphere_map = 0, to deduce ne in
  // nelem_global = 6 ne^2. That is b/c cubed_sphere_map = 0 is supported in
  // Homme only for regular meshes (not RRM), and ne is essential to using the
  // efficiency it provides.
  void init_meta_data(const Int nelem_global, const Int* lid2facenum);

  template <typename Array3D>
  void init_local_mesh_if_needed(const Int ie, const Array3D& corners,
                                 const Real* p_inside);

  // Check that our ref <-> sphere map agrees with Homme's. p_homme is a GLL
  // point on the sphere. Check that we map it to a GLL ref point.
  void check_ref2sphere(const Int ie, const Real* p_homme);

  const LocalMesh& local_mesh (const Int ie) const {
    slmm_assert(ie < static_cast<Int>(local_mesh_.size()));
    return local_mesh_[ie];
  }
  LocalMesh& local_mesh (const Int ie) {
    slmm_assert(ie < static_cast<Int>(local_mesh_.size()));
    return local_mesh_[ie];
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

  bool nearest_point_permitted (const Int& lev) const {
    return lev <= nearest_point_permitted_lev_bdy_;
  }

  const SphereToRef& s2r () const { return s2r_; }

private:
  const Alg::Enum alg_;
  const Int np_, np2_, np4_, cubed_sphere_map_;
  std::vector<LocalMesh> local_mesh_;
  // For CISL:
  const Int tq_order_;
  std::vector<Real> mass_tgt_, mass_mix_, rhs_;
  // For recovery from get_src_cell failure:
  Int nearest_point_permitted_lev_bdy_;
  std::vector<MeshNearestPointData> local_mesh_nearest_point_data_;
  // Meta data obtained at initialization that can be used later.
  Ints lid2facenum_;
  SphereToRef s2r_;
};

void Advecter
::init_meta_data (const Int nelem_global, const Int* lid2facenum) {
  const auto nelemd = local_mesh_.size();
  lid2facenum_ = Ints("Advecter::lid2facenum", nelemd);
  std::copy(lid2facenum, lid2facenum + nelemd, lid2facenum_.data());
  s2r_.init(cubed_sphere_map_, nelem_global, lid2facenum_);
}

template <typename Array3D>
void Advecter::init_local_mesh_if_needed (const Int ie, const Array3D& corners,
                                          const Real* p_inside) {
  slmm_assert(ie < static_cast<Int>(local_mesh_.size()));
  if (local_mesh_[ie].p.extent_int(0) != 0) return;
  auto& m = local_mesh_[ie];
  const Int
    nd = 3,
    nvert = corners.extent_int(1),
    ncell = corners.extent_int(2),
    N = nvert*ncell;
  m.p = typename LocalMesh::RealArray("p", N);
  m.e = typename LocalMesh::IntArray("e", ncell, nvert);
  for (Int ci = 0, k = 0; ci < ncell; ++ci)
    for (Int vi = 0; vi < nvert; ++vi, ++k) {
      for (int j = 0; j < nd; ++j)
        m.p(k,j) = corners(j,vi,ci);
      m.e(ci,vi) = k;
    }
  siqk::test::fill_normals<siqk::SphereGeometry>(m);
  m.tgt_elem = slmm::get_src_cell(m, p_inside);
  slmm_assert(m.tgt_elem >= 0 &&
              m.tgt_elem < ncell);
}

void Advecter::fill_nearest_points_if_needed () {
  if (nearest_point_permitted_lev_bdy_ >= 0)
    for (Int ie = 0; ie < static_cast<Int>(local_mesh_.size()); ++ie)
      init_nearest_point_data(local_mesh_[ie],
                              local_mesh_nearest_point_data_[ie]);
}

void Advecter::check_ref2sphere (const Int ie, const Real* p_homme) {
  const auto& m = local_mesh(ie);
  Real ref_coord[2];
  siqk::sqr::Info info;
  const Real tol = s2r_.tol();
  s2r_.calc_sphere_to_ref(ie, m, p_homme, ref_coord[0], ref_coord[1], &info);
  const slmm::Basis basis(4, 0);
  const slmm::GLL gll;
  const Real* x, * wt;
  gll.get_coef(basis, x, wt);
  int fnd[2] = {0};
  Real min[2] = {1,1};
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 2; ++j) {
      const Real d = std::abs(ref_coord[j] - x[i]);
      min[j] = std::min(min[j], d);
      if (d < tol)
        fnd[j] = 1;
    }
  if ( ! fnd[0] || ! fnd[1])
    printf("COMPOSE check_ref2sphere: %1.15e %1.15e (%1.2e %1.2e) %d %d\n",
           ref_coord[0], ref_coord[1], min[0], min[1],
           info.success, info.n_iterations);
  if ( ! s2r_.check(ie, m))
    printf("COMPOSE SphereToRef::check return false: ie = %d\n", ie);
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
    nerr += slmm::unittest();
    std::cerr << (nerr ? "FAIL" : "PASS") << "ED\n";
  } catch (const std::exception& e) {
    std::cerr << e.what();
  }
  Kokkos::finalize();
  return ret;
}
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

static slmm::Advecter::Ptr g_advecter;

void slmm_init (const Int np, const Int nelem, const Int nelemd,
                const Int transport_alg, const Int cubed_sphere_map,
                const Int sl_nearest_point_lev, const Int* lid2facenum) {
  g_advecter = std::make_shared<slmm::Advecter>(
    np, nelemd, transport_alg, cubed_sphere_map, sl_nearest_point_lev);
  g_advecter->init_meta_data(nelem, lid2facenum);
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

int wait (Request* req, MPI_Status* stat = nullptr) {
#ifdef COMPOSE_DEBUG_MPI
  const auto out = MPI_Wait(&req->request, stat ? stat : MPI_STATUS_IGNORE);
  req->unfreed--;
  return out;
#else
  return MPI_Wait(reinterpret_cast<MPI_Request*>(req), stat ? stat : MPI_STATUS_IGNORE);
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
  void init (const Int nlist, const Int* nlist_per_list, T* buf = nullptr) {
    slmm_assert(nlist >= 0); 
    ptr_.resize(nlist+1);
    ptr_[0] = 0;
    for (Int i = 0; i < nlist; ++i) {
      slmm_assert(nlist_per_list[i] >= 0);
      ptr_[i+1] = ptr_[i] + nlist_per_list[i];
    }
    if (buf) {
      d_ = buf;
    } else {
      v_.resize(ptr_.back());
      d_ = v_.data();
    }
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
  std::vector<T> v_;
  std::vector<Int> ptr_;
  T* d_;
  T* data () { return d_; }
  const T* data () const { return d_; }
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
      lid_on_rank,     // the local ID of the cell on the owning rank
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
    Int nin1halo;                    // nbrs[0:n]
    FixedCapList<OwnItem> own;       // points whose q are computed with own rank's data
    FixedCapList<RemoteItem> rmt;    // points computed by a remote rank's data
    IntArray2D src;                  // src(lev,k) = get_src_cell
    Array<Real**[2]> q_extrema;
    const Real* metdet, * qdp, * dp; // the owned cell's data
    Real* q;
  };

  const mpi::Parallel::Ptr p;
  const slmm::Advecter::ConstPtr advecter;
  const Int np, np2, nlev, qsize, qsized, nelemd, halo;

  FixedCapList<ElemData> ed; // this rank's owned cells, indexed by LID

  // IDs.
  FixedCapList<Int> ranks, nx_in_rank, mylid_with_comm, mylid_with_comm_tid_ptr;
  ListOfLists <Int> nx_in_lid, lid_on_rank;
  BufferLayoutArray bla;

  // MPI comm data.
  ListOfLists<Real> sendbuf, recvbuf;
  FixedCapList<Int> sendcount;
  FixedCapList<mpi::Request> sendreq, recvreq;
  FixedCapList<Int> rmt_xs, rmt_qs_extrema, x_bulkdata_offset;
  Int nrmt_xs, nrmt_qs_extrema;

  bool horiz_openmp;
#ifdef HORIZ_OPENMP
  ListOfLists<omp_lock_t> ri_lidi_locks;
#endif

  // temporary work space
  std::vector<Int> nlid_per_rank, sendsz, recvsz;
  Array<Real**> rwork;

  CslMpi (const mpi::Parallel::Ptr& ip, const slmm::Advecter::ConstPtr& advecter,
          Int inp, Int inlev, Int iqsize, Int iqsized, Int inelemd, Int ihalo)
    : p(ip), advecter(advecter),
      np(inp), np2(np*np), nlev(inlev), qsize(iqsize), qsized(iqsized), nelemd(inelemd),
      halo(ihalo)
  {}

  ~CslMpi () {
#ifdef HORIZ_OPENMP
    const Int nrmtrank = static_cast<Int>(ranks.n()) - 1;
    for (Int ri = 0; ri < nrmtrank; ++ri) {
      auto&& locks = ri_lidi_locks(ri);
      for (auto& lock: locks) {
        // The following call is causing a seg fault on at least one Cray KNL
        // system. It doesn't kill the run b/c it occurs after main exits. Not
        // calling this is a memory leak, but it's innocuous in this case. Thus,
        // comment it out:
        //omp_destroy_lock(&lock);
        // It's possible that something in the OpenMP runtime shuts down before
        // this call, and that's the problem. If that turns out to be it, I
        // should make a compose_finalize() call somewhere.
      }
    }
#endif
  }
};

namespace extend_halo {
// Extend halo by one layer. This has two parts: finding neighbor (gid, rank) in
// collect_gid_rank, and extending the Advecter local mesh geometry in
// extend_local_meshes. The two parts have the same comm pattern: in round 1,
// request data for lists of GIDs; in round 2, fulfill these requests.

typedef CslMpi::ElemData ElemData;
typedef Int Gid;
typedef Int Rank;
struct GidRankPair {
  const Gid gid;
  const Rank rank;
  GidRankPair(const Gid& gid_, const Rank& rank_) : gid(gid_), rank(rank_) {}
  bool operator< (const GidRankPair& o) const {
    if (gid < o.gid) return true;
    if (gid > o.gid) return false;
    return rank < o.rank;
  }
};
typedef std::vector<GidRankPair> GidRankPairs;
typedef std::map<Gid, GidRankPairs> Gid2Nbrs;
typedef std::vector<Int> IntBuf;
typedef std::vector<Real> RealBuf;

GidRankPairs all_nbrs_but_me (const ElemData& ed) {
  GidRankPairs gs;
  gs.reserve(ed.nbrs.size() - 1);
  for (const auto& n : ed.nbrs)
    if (&n != ed.me)
      gs.push_back(GidRankPair(n.gid, n.rank));
  return gs;
}

void fill_gid2nbrs (const mpi::Parallel& p, const FixedCapList<ElemData>& eds,
                    Gid2Nbrs& gid2nbrs) {
  static const Int tag = 6;
  const Rank my_rank = p.rank();
  const Int n_owned = eds.size();

  // Fill in the ones we know.
  for (const auto& ed : eds) {
    slmm_assert(ed.me->rank == my_rank);
    gid2nbrs[ed.me->gid] = all_nbrs_but_me(ed);
  }

  std::vector<Rank> ranks;
  Int nrank;
  std::vector<IntBuf> req_sends, req_recvs;
  std::vector<mpi::Request> req_recv_reqs;
  {
    // Find the ranks that know the rest.
    std::map<Rank,Int> rank2rankidx;
    std::map<Gid,Rank> needgid2rank;
    {
      std::set<Rank> unique_ranks;
      for (const auto& item : gid2nbrs)
        for (const auto& n : item.second)
          if (n.rank != my_rank) {
            slmm_assert(gid2nbrs.find(n.gid) == gid2nbrs.end());
            needgid2rank.insert(std::make_pair(n.gid, n.rank));
            unique_ranks.insert(n.rank);
          }
      nrank = unique_ranks.size();
      ranks.insert(ranks.begin(), unique_ranks.begin(), unique_ranks.end());
      Int i = 0;
      for (const auto& rank : ranks)
        rank2rankidx[rank] = i++;
    }

    // Send and receive neighbor queries.
    slmm_assert(ranks.size() == static_cast<size_t>(nrank));
    req_sends.resize(nrank);
    req_recvs.resize(nrank); req_recv_reqs.resize(nrank);
    for (Int i = 0; i < nrank; ++i) {
      auto& r = req_recvs[i];
      r.resize(n_owned+1); // upper bound on #requests, plus 1 for size datum
      mpi::irecv(p, r.data(), r.size(), ranks[i], tag, &req_recv_reqs[i]);
      req_sends[i].push_back(0);
    }
    for (const auto& item : needgid2rank) {
      const auto& gid = item.first;
      const auto& rank = item.second;
      const auto& it = rank2rankidx.find(rank);
      slmm_assert(it != rank2rankidx.end());
      const auto& rank_idx = it->second;
      slmm_assert(gid2nbrs.find(gid) == gid2nbrs.end());
      req_sends[rank_idx].push_back(gid);
    }
    for (Int i = 0; i < nrank; ++i) {
      auto& s = req_sends[i];
      s[0] = s.size() - 1; // #gids in request
      mpi::isend(p, s.data(), s.size(), ranks[i], tag, nullptr);
    }
  }

  // Fullfill queries and receive answers to our queries.
  std::vector<IntBuf> nbr_sends(nrank), nbr_recvs(nrank);
  std::vector<mpi::Request> nbr_send_reqs(nrank), nbr_recv_reqs(nrank);
  for (Int i = 0; i < nrank; ++i) {
    auto& r = nbr_recvs[i];
    // 20 is from dimensions_mod::set_mesh_dimensions; factor of 2 is to get
    // (gid,rank); 1 is for size datum.
    r.resize((20*2 + 1)*(req_sends[i].size() - 1));
    mpi::irecv(p, r.data(), r.size(), ranks[i], tag, &nbr_recv_reqs[i]);
  }
  for (Int k = 0; k < nrank; ++k) {
    Int i;
    mpi::waitany(nrank, req_recv_reqs.data(), &i);
    const auto& r = req_recvs[i];
    auto& s = nbr_sends[i];
    const Int ngid = r[0];
    for (Int j = 0; j < ngid; ++j) {
      const auto gid = r[j+1];
      const auto it = gid2nbrs.find(gid);
      slmm_assert(it != gid2nbrs.end());
      const auto& nbrs = it->second;
      s.push_back(nbrs.size());
      for (const auto& n : nbrs) {
        s.push_back(n.gid);
        s.push_back(n.rank);
      }
    }
    mpi::isend(p, s.data(), s.size(), ranks[i], tag, &nbr_send_reqs[i]);
  }
  for (Int j = 0; j < nrank; ++j) {
    Int i;
    mpi::waitany(nrank, nbr_recv_reqs.data(), &i);
    const auto& s = req_sends[i];
    const Int nrequested = s.size() - 1;
    const auto& r = nbr_recvs[i];
    for (Int si = 0, ri = 0; si < nrequested; ++si) {
      const Gid gid = s[si+1];
      const Int nnbr = r[ri++];
      GidRankPairs nbrs;
      for (Int ni = 0; ni < nnbr; ++ni, ri += 2)
        nbrs.push_back(GidRankPair(r[ri], r[ri+1]));
      slmm_assert(gid2nbrs.find(gid) == gid2nbrs.end());
      gid2nbrs.insert(std::make_pair(gid, nbrs));
    }
  }
  mpi::waitall(nbr_send_reqs.size(), nbr_send_reqs.data());
}

void extend_nbrs (const Gid2Nbrs& gid2nbrs, FixedCapList<ElemData>& eds) {
  for (auto& ed : eds) {
    // Get all <=2-halo neighbors.
    std::set<GidRankPair> new_nbrs;
    for (const auto& n : ed.nbrs) {
      if (&n == ed.me) continue;
      const auto& it = gid2nbrs.find(n.gid);
      slmm_assert(it != gid2nbrs.end());
      const auto& gid_nbrs = it->second;
      for (const auto& gn : gid_nbrs)
        new_nbrs.insert(gn);
    }
    // Remove the already known ones.
    for (const auto& n : ed.nbrs)
      new_nbrs.erase(GidRankPair(n.gid, n.rank));
    // Find me b/c of the pointer reset.
    Int me = -1;
    for (Int i = 0; i < ed.nbrs.size(); ++i)
      if (&ed.nbrs(i) == ed.me) {
        me = i;
        break;
      }
    slmm_assert(me >= 0);
    // Append the, now only new, 2-halo ones.
    Int i = ed.nbrs.size();
    ed.nbrs.reset_capacity(i + new_nbrs.size(), true);
    ed.me = &ed.nbrs(me);
    for (const auto& n : new_nbrs) {
      auto& en = ed.nbrs(i++);
      en.gid = n.gid;
      en.rank = n.rank;
      en.rank_idx = -1;
      en.lid_on_rank = -1;
      en.lid_on_rank_idx = -1;
    }
  }
}

void collect_gid_rank (const mpi::Parallel& p, FixedCapList<ElemData>& eds) {
  Gid2Nbrs gid2nbrs;
  fill_gid2nbrs(p, eds, gid2nbrs);
  extend_nbrs(gid2nbrs, eds);
}

void extend_local_meshes (const mpi::Parallel& p, const FixedCapList<ElemData>& eds,
                          slmm::Advecter& advecter) {
  using slmm::slice;
  using slmm::nslices;
  using slmm::szslice;
  using slmm::len;
  using slmm::LocalMesh;
  using RealArray3 = ko::View<Real***, ko::HostSpace>;

  static const Int tag = 24;
  const Int my_rank = p.rank();
  const Int n_owned = eds.size();

  RealArray3 corners;
  std::map<Gid,Int> owngid2lid, rmtgid2idx;
  {
    // Find relevant ranks.
    std::vector<Int> ranks;
    Int nrank;
    std::vector<IntBuf> req_sends, req_recvs;
    std::vector<mpi::Request> req_recv_reqs;
    {
      std::map<Int,Int> rank2rankidx;
      {
        std::set<Int> uranks;
        for (const auto& ed : eds)
          for (const auto& n : ed.nbrs)
            uranks.insert(n.rank);
        uranks.erase(my_rank);
        ranks.insert(ranks.begin(), uranks.begin(), uranks.end());
        nrank = ranks.size();
        Int i = 0;
        for (const auto& rank : ranks)
          rank2rankidx[rank] = i++;
      }

      // Trade requests.
      req_sends.resize(nrank); req_recvs.resize(nrank);
      req_recv_reqs.resize(nrank);
      for (Int i = 0; i < nrank; ++i) { // Set up recvs.
        auto& r = req_recvs[i];
        r.resize(n_owned+1);
        mpi::irecv(p, r.data(), r.size(), ranks[i], tag, &req_recv_reqs[i]);
      }
      std::set<Gid> unique_rmtgids;
      { // Collect the remote GIDs by rank.
        std::vector<std::set<Gid> > req_gids(nrank);
        for (const auto& ed : eds)
          for (Int i = ed.nin1halo; i < ed.nbrs.size(); ++i) {
            const auto& n = ed.nbrs(i);
            if (n.rank != my_rank) {
              req_gids[rank2rankidx[n.rank]].insert(n.gid);
              unique_rmtgids.insert(n.gid);
            }
          }
        for (Int i = 0; i < nrank; ++i) { // Set up sends.
          auto& s = req_sends[i];
          s.push_back(0);
          s.insert(s.end(), req_gids[i].begin(), req_gids[i].end());
          s[0] = s.size() - 1;
          mpi::isend(p, s.data(), s.size(), ranks[i], tag, nullptr);
        }
        Int i = 0;
        for (const auto& gid : unique_rmtgids)
          rmtgid2idx[gid] = i++;
      }
    }

    for (Int i = 0; i < eds.size(); ++i) {
      const Gid gid = eds(i).me->gid;
      owngid2lid[gid] = i;
    }

    // Fulfill our requests and obtain the data in reply to ours.
    std::vector<RealBuf> geo_recvs(nrank), geo_sends(nrank);
    std::vector<mpi::Request> geo_recv_reqs(nrank), geo_send_reqs(nrank);
    for (Int i = 0; i < nrank; ++i) {
      auto& r = geo_recvs[i];
      r.resize(12*(req_sends[i].size() - 1)); // 4 vertices x 3 dimensions
      mpi::irecv(p, r.data(), r.size(), ranks[i], tag, &geo_recv_reqs[i]);
    }
    for (Int k = 0; k < nrank; ++k) { // Pack cell corner points for requesters.
      Int i;
      mpi::waitany(nrank, req_recv_reqs.data(), &i);
      const auto& r = req_recvs[i];
      auto& s = geo_sends[i];
      const Int ngid = r[0];
      s.resize(12*ngid);
      Int si = 0;
      for (Int j = 0; j < ngid; ++j) {
        const auto gid = r[j+1];
        const auto it = owngid2lid.find(gid);
        slmm_assert(it != owngid2lid.end());
        const auto lid = it->second;
        const auto& local_mesh = advecter.local_mesh(lid);
        const auto cell = slice(local_mesh.e, local_mesh.tgt_elem);
        slmm_assert(szslice(local_mesh.e) == 4 && cell[3] != -1);
        const auto& p = local_mesh.p;
        for (Int v = 0; v < 4; ++v)
          for (Int d = 0; d < 3; ++d)
            s[si++] = p(cell[v],d);
      }
      slmm_assert(si == static_cast<Int>(s.size()));
      mpi::isend(p, s.data(), s.size(), ranks[i], tag, &geo_send_reqs[i]);
    }
    corners = RealArray3("corners", rmtgid2idx.size(), 4, 3);
    for (Int j = 0; j < nrank; ++j) { // Pack cell corner points for me.
      Int i;
      mpi::waitany(nrank, geo_recv_reqs.data(), &i);
      const auto& s = req_sends[i];
      const Int nrequested = s.size() - 1;
      const auto& r = geo_recvs[i];
      for (Int si = 0, ri = 0; si < nrequested; ++si) {
        const Gid gid = s[si+1];
        const auto it = rmtgid2idx.find(gid);
        slmm_assert(it != rmtgid2idx.end());
        const auto idx = it->second;
        for (Int v = 0; v < 4; ++v)
          for (Int d = 0; d < 3; ++d)
            corners(idx,v,d) = r[ri++];
      }
    }
    // Do this here so we can release the memory.
    mpi::waitall(geo_send_reqs.size(), geo_send_reqs.data());
  }

  // Extend the local meshes.
  for (Int lid = 0; lid < eds.size(); ++lid) {
    const auto& ed = eds(lid);
    auto& local_mesh = advecter.local_mesh(lid);
    auto p0 = local_mesh.p;
    auto e0 = local_mesh.e;
    const Int ncell = ed.nbrs.size();
    slmm_assert(szslice(p0) == 3 && szslice(e0) == 4);
    local_mesh.p = LocalMesh::RealArray("p", 4*ncell);
    auto& p = local_mesh.p;
    local_mesh.e = LocalMesh::IntArray("e", ncell, 4);
    auto& e = local_mesh.e;
    // Copy in old data.
    for (Int pi = 0; pi < nslices(p0); ++pi)
      for (Int d = 0; d < 3; ++d)
        p(pi,d) = p0(pi,d);
    for (Int ei = 0; ei < nslices(e0); ++ei)
      for (Int d = 0; d < 4; ++d)
        e(ei,d) = e0(ei,d);
    // Fill in new data.
    for (Int ni = ed.nin1halo; ni < ed.nbrs.size(); ++ni) {
      const auto& n = ed.nbrs(ni);
      const Gid gid = n.gid;
      if (n.rank != my_rank) {
        const auto it = rmtgid2idx.find(gid);
        slmm_assert(it != rmtgid2idx.end());
        const Int idx = it->second;
        for (Int v = 0; v < 4; ++v) {
          const Int slot = 4*ni+v;
          e(ni,v) = slot;
          slmm_assert(slot < nslices(p));
          for (Int d = 0; d < 3; ++d) p(slot,d) = corners(idx,v,d);
        }
      } else {
        const auto it = owngid2lid.find(gid);
        slmm_assert(it != owngid2lid.end());
        const Int lid = it->second;
        const auto& local_mesh_other = advecter.local_mesh(lid);
        const auto cell_other = slice(local_mesh_other.e, local_mesh_other.tgt_elem);
        const auto& p_other = local_mesh_other.p;
        for (Int v = 0; v < 4; ++v) {
          const Int slot = 4*ni+v;
          e(ni,v) = slot;
          slmm_assert(slot < nslices(p));
          for (Int d = 0; d < 3; ++d) p(slot,d) = p_other(cell_other[v],d);
        }
      }
    }
    // Recompute all normals.
    siqk::test::fill_normals<siqk::SphereGeometry>(local_mesh);
  }
}
} // namespace extend_halo

// Fill in (gid, rank), the list of owning rank per gid.
void collect_gid_rank (CslMpi& cm, const Int* nbr_id_rank, const Int* nirptr) {
  cm.ed.reset_capacity(cm.nelemd, true);
  for (Int i = 0; i < cm.nelemd; ++i) {
    auto& ed = cm.ed(i);
    const Int* nir = nbr_id_rank + nirptr[i];
    const Int nnir = (nirptr[i+1] - nirptr[i]) / 2;
    const Int mygid = nir[0];
    ed.me = nullptr;
    ed.nin1halo = nnir-1;
    ed.nbrs.reset_capacity(ed.nin1halo, true);
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
      n.rank_idx = -1;
      n.lid_on_rank = -1;
      n.lid_on_rank_idx = -1;
    }
    slmm_assert(ed.me);
  }
  if (cm.halo == 2) extend_halo::collect_gid_rank(*cm.p, cm.ed);
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
  cm.x_bulkdata_offset.reset_capacity(i, true);

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
// has a 1-halo patch of bulk data. For a 1-halo, allocations in this routine
// use essentially the same amount of memory, but not more. We could use less if
// we were willing to realloc space at each SL time step.
void size_mpi_buffers (CslMpi& cm, const Rank2Gids& rank2rmtgids,
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
  cm.nlid_per_rank.resize(nrmtrank);
  cm.sendsz.resize(nrmtrank);
  cm.recvsz.resize(nrmtrank);
  Int rmt_xs_sz = 0, rmt_qse_sz = 0;
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    const auto& rmtgids = rank2rmtgids.at(cm.ranks(ri));
    const auto& owngids = rank2owngids.at(cm.ranks(ri));
    cm.nlid_per_rank[ri] = rmtgids.size();
    cm.sendsz[ri] = bytes2real(std::max(xbufcnt(rmtgids, owngids),
                                        qbufcnt(owngids, rmtgids)));
    cm.recvsz[ri] = bytes2real(std::max(xbufcnt(owngids, rmtgids),
                                        qbufcnt(rmtgids, owngids)));
    rmt_xs_sz  += 5*cm.np2*cm.nlev*rmtgids.size();
    rmt_qse_sz += 4       *cm.nlev*rmtgids.size();
  }
  cm.rmt_xs.reset_capacity(rmt_xs_sz, true);
  cm.rmt_qs_extrema.reset_capacity(rmt_qse_sz, true);
}

void alloc_mpi_buffers (CslMpi& cm, Real* sendbuf = nullptr, Real* recvbuf = nullptr) {
  const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
  cm.nx_in_rank.reset_capacity(nrmtrank, true);
  cm.nx_in_lid.init(nrmtrank, cm.nlid_per_rank.data());
  cm.bla.init(nrmtrank, cm.nlid_per_rank.data(), cm.nlev);
  cm.sendbuf.init(nrmtrank, cm.sendsz.data(), sendbuf);
  cm.recvbuf.init(nrmtrank, cm.recvsz.data(), recvbuf);
  cm.nlid_per_rank.clear();
  cm.sendsz.clear();
  cm.recvsz.clear();
#ifdef HORIZ_OPENMP
  cm.ri_lidi_locks.init(nrmtrank, cm.nlid_per_rank.data());
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
  size_mpi_buffers(cm, rank2rmtgids, rank2owngids);
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
    const int nthr = get_num_threads();
    cm.rwork = CslMpi::Array<Real**>("rwork", nthr, cm.qsize);
    cm.mylid_with_comm_tid_ptr.reset_capacity(nthr+1, true);
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

CslMpi::Ptr init (const slmm::Advecter::ConstPtr& advecter,
                  const mpi::Parallel::Ptr& p,
                  Int np, Int nlev, Int qsize, Int qsized, Int nelemd,
                  const Int* nbr_id_rank, const Int* nirptr,
                  Int halo) {
  slmm_throw_if(halo < 1 || halo > 2, "halo must be 1 (default) or 2.");
  auto cm = std::make_shared<CslMpi>(p, advecter, np, nlev, qsize, qsized,
                                     nelemd, halo);
  setup_comm_pattern(*cm, nbr_id_rank, nirptr);
  return cm;
}

// For const clarity, take the non-const advecter as an arg, even though cm
// already has a ref to the const'ed one.
void finalize_local_meshes (CslMpi& cm, slmm::Advecter& advecter) {
  if (cm.halo == 2) extend_halo::extend_local_meshes(*cm.p, cm.ed, advecter);
  advecter.fill_nearest_points_if_needed();
}

// Set pointers to HOMME data arrays.
void set_elem_data (CslMpi& cm, const Int ie, const Real* metdet, const Real* qdp,
                    const Real* dp, Real* q, const Int nelem_in_patch) {
  slmm_assert(ie < cm.ed.size());
  slmm_assert(cm.halo > 1 || cm.ed(ie).nbrs.size() == nelem_in_patch);
  auto& e = cm.ed(ie);
  e.metdet = metdet;
  e.qdp = qdp;
  e.dp = dp;
  e.q = q;
}

void setup_irecv (CslMpi& cm, const bool skip_if_empty = false) {
#ifdef HORIZ_OPENMP
# pragma omp master
#endif
  {
    const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
    cm.recvreq.clear();
    for (Int ri = 0; ri < nrmtrank; ++ri) {
      if (skip_if_empty && cm.nx_in_rank(ri) == 0) continue;
      auto&& recvbuf = cm.recvbuf(ri);
      // The count is just the number of slots available, which can be larger
      // than what is actually being received.
      cm.recvreq.inc();
      mpi::irecv(*cm.p, recvbuf.data(), recvbuf.n(), cm.ranks(ri), 42,
                 &cm.recvreq.back());
    }
  }
}

void isend (CslMpi& cm, const bool want_req = true, const bool skip_if_empty = false) {
#ifdef HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
  {
    const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
    for (Int ri = 0; ri < nrmtrank; ++ri) {
      if (skip_if_empty && cm.sendcount(ri) == 0) continue;
      mpi::isend(*cm.p, cm.sendbuf(ri).data(), cm.sendcount(ri),
                 cm.ranks(ri), 42, want_req ? &cm.sendreq(ri) : nullptr);
    }
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

void wait_on_send (CslMpi& cm, const bool skip_if_empty = false) {
#ifdef HORIZ_OPENMP
# pragma omp master
#endif
  {
    for (Int ri = 0; ri < cm.sendreq.n(); ++ri) {
      if (skip_if_empty && cm.sendcount(ri) == 0) continue;
      mpi::wait(&cm.sendreq(ri));
    }
  }
#ifdef HORIZ_OPENMP
# pragma omp barrier
#endif
}

void recv (CslMpi& cm, const bool skip_if_empty = false) {
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
    const auto& mesh = cm.advecter->local_mesh(tci);
    const auto tgt_idx = mesh.tgt_elem;
    auto& ed = cm.ed(tci);
    ed.own.clear();
    for (Int lev = 0; lev < cm.nlev; ++lev)
      for (Int k = 0; k < cm.np2; ++k) {
        Int sci = slmm::get_src_cell(mesh, &dep_points(0,k,lev,tci), tgt_idx);
        if (sci == -1 && cm.advecter->nearest_point_permitted(lev))
          sci = slmm::get_nearest_point(
            mesh, cm.advecter->nearest_point_data(tci),
            &dep_points(0,k,lev,tci), tgt_idx);
        if (sci == -1) {
          std::stringstream ss;
          ss.precision(17);
          const auto* v = &dep_points(0,k,lev,tci);
          ss << "Departure point is outside of halo:\n"
             << "  nearest point permitted: "
             << cm.advecter->nearest_point_permitted(lev)
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
   metadata. Two passes let us do some efficient tricks that are not available
   with one pass. Departure point and q messages are formatted as follows:
    xs: (#x-in-rank         int                                     <-
         x-bulk-data-offset i                                        |
         (lid-on-rank       i     only packed if #x in lid > 0       |
          #x-in-lid         i     > 0                                |- meta data
          (lev              i     only packed if #x in (lid,lev) > 0 |
           #x)              i     > 0                                |
              *#lev) *#lid                                          <-
         x                  3 real                                  <-- bulk data
          *#x-in-rank) *#rank
    qs: (q-extrema    2 qsize r   (min, max) packed together
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
    // metadata offset, x bulk data offset, q bulk data offset
    Int mos = 0, xos = 0, qos = 0, sendcount = 0, cnt;
    cnt = setbuf(sendbuf, mos, 0, 0); // empty space for later
    mos += cnt;
    sendcount += cnt;
    if (cm.nx_in_rank(ri) == 0) {
      setbuf(sendbuf, 0, mos, 0);
      cm.x_bulkdata_offset(ri) = mos;
      cm.sendcount(ri) = sendcount;
      continue;
    }
    auto&& bla = cm.bla(ri);
    for (Int lidi = 0, lidn = lid_on_rank.n(); lidi < lidn; ++lidi) {
      auto nx_in_lid = cm.nx_in_lid(ri,lidi);
      if (nx_in_lid == 0) continue;
      cnt = setbuf(sendbuf, mos, lid_on_rank(lidi), nx_in_lid);
      mos += cnt;
      sendcount += cnt;
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
        const auto dos = setbuf(sendbuf, mos, lev, nx);
        mos += dos;
        sendcount += dos + 3*nx;
        t.xptr = xos;
        xos += 3*nx;
        qos += 2 + nx;
        nx_in_lid -= nx;
      }
      slmm_assert(nx_in_lid == 0);
    }
    setbuf(sendbuf, 0, mos /* offset to x bulk data */, cm.nx_in_rank(ri));
    cm.x_bulkdata_offset(ri) = mos;
    cm.sendcount(ri) = sendcount;
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
          xptr = cm.x_bulkdata_offset(ri) + t.xptr + 3*cnt;
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

static inline Real
calc_q_use_q (const Real rx[4], const Real ry[4], const Real* const qs) {
  return (ry[0]*(rx[0]*qs[ 0] + rx[1]*qs[ 1] + rx[2]*qs[ 2] + rx[3]*qs[ 3]) +
          ry[1]*(rx[0]*qs[ 4] + rx[1]*qs[ 5] + rx[2]*qs[ 6] + rx[3]*qs[ 7]) +
          ry[2]*(rx[0]*qs[ 8] + rx[1]*qs[ 9] + rx[2]*qs[10] + rx[3]*qs[11]) +
          ry[3]*(rx[0]*qs[12] + rx[1]*qs[13] + rx[2]*qs[14] + rx[3]*qs[15]));
}

static inline Real
calc_q_use_qdp (const Real rx[4], const Real ry[4],
                const Real* const dp, const Real* const qdp) {
  return (ry[0]*(rx[0]*(qdp[ 0]/dp[ 0]) + rx[1]*(qdp[ 1]/dp[ 1])  +
                 rx[2]*(qdp[ 2]/dp[ 2]) + rx[3]*(qdp[ 3]/dp[ 3])) +
          ry[1]*(rx[0]*(qdp[ 4]/dp[ 4]) + rx[1]*(qdp[ 5]/dp[ 5])  +
                 rx[2]*(qdp[ 6]/dp[ 6]) + rx[3]*(qdp[ 7]/dp[ 7])) +
          ry[2]*(rx[0]*(qdp[ 8]/dp[ 8]) + rx[1]*(qdp[ 9]/dp[ 9])  +
                 rx[2]*(qdp[10]/dp[10]) + rx[3]*(qdp[11]/dp[11])) +
          ry[3]*(rx[0]*(qdp[12]/dp[12]) + rx[1]*(qdp[13]/dp[13])  +
                 rx[2]*(qdp[14]/dp[14]) + rx[3]*(qdp[15]/dp[15])));
}

template <Int np>
void calc_q (const CslMpi& cm, const Int& src_lid, const Int& lev,
             const Real* const dep_point, Real* const q_tgt, const bool use_q) {
  Real ref_coord[2]; {
    const auto& m = cm.advecter->local_mesh(src_lid);
    cm.advecter->s2r().calc_sphere_to_ref(src_lid, m, dep_point,
                                          ref_coord[0], ref_coord[1]);
  }

  // Interpolate.
  Real rx[4], ry[4];
  switch (cm.advecter->alg()) {
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
  const Int qsize = cm.qsize;
  static const Int blocksize = 8;
  if (use_q) {
    // We can use q from calc_q_extrema.
    const Real* const qs0 = ed.q + levos;
    // It was found that Intel 18 produced code that was not BFB between runs
    // due to this pragma.
    //#pragma ivdep
    for (Int iqo = 0; iqo < qsize; iqo += blocksize) {
      // So, instead, provide the compiler with a clear view of
      // ivdep-ness. Write to tmp here in chunks of blocksize, then
      // move tmp to q_tgt later.
      if (iqo + blocksize <= qsize) {
        Real tmp[blocksize];
        for (Int iqi = 0; iqi < blocksize; ++iqi) {
          const Real* const qs = qs0 + (iqo + iqi)*np2nlev;
          tmp[iqi] = calc_q_use_q(rx, ry, qs);
        }
        for (Int iqi = 0; iqi < blocksize; ++iqi)
          q_tgt[iqo + iqi] = tmp[iqi];
      } else {
        for (Int iq = iqo; iq < qsize; ++iq) {
          const Real* const qs = qs0 + iq*np2nlev;
          q_tgt[iq] = calc_q_use_q(rx, ry, qs);
        }
      }
    }
  } else {
    // q from calc_q_extrema is being overwritten, so have to use qdp/dp.
    const Real* const dp = ed.dp + levos;
    const Real* const qdp0 = ed.qdp + levos;
    // I'm commenting out this pragma, too, to be safe.
    //#pragma ivdep
    for (Int iqo = 0; iqo < qsize; iqo += blocksize) {
      // And I'm using the same technique as above.
      if (iqo + blocksize <= qsize) {
        Real tmp[blocksize];
        for (Int iqi = 0; iqi < blocksize; ++iqi) {
          const Real* const qdp = qdp0 + (iqo + iqi)*np2nlev;
          tmp[iqi] = calc_q_use_qdp(rx, ry, dp, qdp);
        }
        for (Int iqi = 0; iqi < blocksize; ++iqi)
          q_tgt[iqo + iqi] = tmp[iqi];
      } else {
        for (Int iq = iqo; iq < qsize; ++iq) {
          const Real* const qdp = qdp0 + iq*np2nlev;
          q_tgt[iq] = calc_q_use_qdp(rx, ry, dp, qdp);
        }
      }
    }
  }
}

template <Int np>
void calc_rmt_q_pass1 (CslMpi& cm) {
#ifdef HORIZ_OPENMP
# pragma omp master
#endif
  {
    const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
    Int cnt = 0, qcnt = 0;
    for (Int ri = 0; ri < nrmtrank; ++ri) {
      const auto&& xs = cm.recvbuf(ri);
      Int mos = 0, qos = 0, nx_in_rank, xos;
      mos += getbuf(xs, mos, xos, nx_in_rank);
      if (nx_in_rank == 0) {
        cm.sendcount(ri) = 0;
        continue; 
      }
      // The upper bound is to prevent an inf loop if the msg is corrupted.
      for (Int lidi = 0; lidi < cm.nelemd; ++lidi) {
        Int lid, nx_in_lid;
        mos += getbuf(xs, mos, lid, nx_in_lid);
        for (Int levi = 0; levi < cm.nlev; ++levi) { // same re: inf loop
          Int lev, nx;
          mos += getbuf(xs, mos, lev, nx);
          slmm_assert(nx > 0);
          {
            cm.rmt_qs_extrema(4*qcnt + 0) = ri;
            cm.rmt_qs_extrema(4*qcnt + 1) = lid;
            cm.rmt_qs_extrema(4*qcnt + 2) = lev;
            cm.rmt_qs_extrema(4*qcnt + 3) = qos;
            ++qcnt;
            qos += 2;
          }
          for (Int xi = 0; xi < nx; ++xi) {
            cm.rmt_xs(5*cnt + 0) = ri;
            cm.rmt_xs(5*cnt + 1) = lid;
            cm.rmt_xs(5*cnt + 2) = lev;
            cm.rmt_xs(5*cnt + 3) = xos;
            cm.rmt_xs(5*cnt + 4) = qos;
            ++cnt;
            xos += 3;
            ++qos;
          }
          nx_in_lid -= nx;
          nx_in_rank -= nx;
          if (nx_in_lid == 0) break;
        }
        slmm_assert(nx_in_lid == 0);
        if (nx_in_rank == 0) break;
      }
      slmm_assert(nx_in_rank == 0);
      cm.sendcount(ri) = cm.qsize*qos;
    }
    cm.nrmt_xs = cnt;
    cm.nrmt_qs_extrema = qcnt;
  }
#ifdef HORIZ_OPENMP
# pragma omp barrier
#endif  
}

template <Int np>
void calc_rmt_q_pass2 (CslMpi& cm) {
  const Int qsize = cm.qsize;

#ifdef HORIZ_OPENMP
# pragma omp for
#endif
  for (Int it = 0; it < cm.nrmt_qs_extrema; ++it) {
    const Int
      ri = cm.rmt_qs_extrema(4*it), lid = cm.rmt_qs_extrema(4*it + 1),
      lev = cm.rmt_qs_extrema(4*it + 2), qos = qsize*cm.rmt_qs_extrema(4*it + 3);  
    auto&& qs = cm.sendbuf(ri);
    const auto& ed = cm.ed(lid);
    for (Int iq = 0; iq < qsize; ++iq)
      for (int i = 0; i < 2; ++i)
        qs(qos + 2*iq + i) = ed.q_extrema(iq, lev, i);
  }

#ifdef HORIZ_OPENMP
# pragma omp for
#endif
  for (Int it = 0; it < cm.nrmt_xs; ++it) {
    const Int
      ri = cm.rmt_xs(5*it), lid = cm.rmt_xs(5*it + 1), lev = cm.rmt_xs(5*it + 2),
      xos = cm.rmt_xs(5*it + 3), qos = qsize*cm.rmt_xs(5*it + 4);
    const auto&& xs = cm.recvbuf(ri);
    auto&& qs = cm.sendbuf(ri);
    calc_q<np>(cm, lid, lev, &xs(xos), &qs(qos), true);
  }
}

template <Int np>
void calc_rmt_q (CslMpi& cm) {
  calc_rmt_q_pass1<np>(cm);
  calc_rmt_q_pass2<np>(cm);
}

template <Int np>
void calc_own_q (CslMpi& cm, const Int& nets, const Int& nete,
                 const FA4<const Real>& dep_points,
                 const FA4<Real>& q_min, const FA4<Real>& q_max) {
  const int tid = get_tid();
  for (Int tci = 0; tci < cm.nelemd; ++tci) {
    auto& ed = cm.ed(tci);
    const FA3<Real> q_tgt(ed.q, cm.np2, cm.nlev, cm.qsize);
    const Int ned = ed.own.n();
#ifdef HORIZ_OPENMP
    #pragma omp for
#endif
    for (Int idx = 0; idx < ned; ++idx) {
      const auto& e = ed.own(idx);
      const Int slid = ed.nbrs(ed.src(e.lev, e.k)).lid_on_rank;
      const auto& sed = cm.ed(slid);
      for (Int iq = 0; iq < cm.qsize; ++iq) {
        q_min(e.k, e.lev, iq, tci) = sed.q_extrema(iq, e.lev, 0);
        q_max(e.k, e.lev, iq, tci) = sed.q_extrema(iq, e.lev, 1);
      }
      Real* const qtmp = &cm.rwork(tid, 0);
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
    auto& ed = cm.ed(tci);
    const FA3<Real> q_tgt(ed.q, cm.np2, cm.nlev, cm.qsize);
    for (const auto& e: ed.rmt) {
      slmm_assert(ed.nbrs(ed.src(e.lev, e.k)).rank != myrank);
      const Int ri = ed.nbrs(ed.src(e.lev, e.k)).rank_idx;
      const auto&& recvbuf = cm.recvbuf(ri);
      for (Int iq = 0; iq < cm.qsize; ++iq) {
        q_min(e.k, e.lev, iq, tci) = recvbuf(e.q_extrema_ptr + 2*iq    );
        q_max(e.k, e.lev, iq, tci) = recvbuf(e.q_extrema_ptr + 2*iq + 1);
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
    q_min(q_min_r, cm.np2, cm.nlev, cm.qsize, cm.nelemd),
    q_max(q_max_r, cm.np2, cm.nlev, cm.qsize, cm.nelemd);

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
  isend(cm, true /* want_req */, true /* skip_if_empty */);
  // Set up to receive q for each of my departure point requests sent to
  // remotes. We can't do this until the OpenMP barrier in isend assures that
  // all threads are done with the receive buffer's departure points.
  setup_irecv(cm, true /* skip_if_empty */);
  // While waiting to get my data from remotes, compute q for departure points
  // that have remained in my elements.
  calc_own_q<np>(cm, nets, nete, dep_points, q_min, q_max);
  // Receive remote q data and use this to fill in the rest of my fields.
  recv(cm, true /* skip_if_empty */);
  copy_q(cm, nets, q_min, q_max);
  // Wait on send buffer so it's free to be used by others.
  wait_on_send(cm, true /* skip_if_empty */);
}

} // namespace cslmpi
} // namespace homme

// Valid after slmm_init_local_mesh_ is called.
int slmm_unittest () {
  int nerr = 0, ne;
  {
    ne = 0;
    for (int i = 0; i < homme::g_advecter->nelem(); ++i) {
      const auto& m = homme::g_advecter->local_mesh(i);
      ne += slmm::unittest(m, m.tgt_elem);
    }
    if (ne)
      fprintf(stderr, "slmm_unittest: slmm::unittest returned %d\n", ne);
    nerr += ne;
  }
  return nerr;
}

#include <cstdlib>

struct Experiment {
  int sl_mpi, halo;
};

template <typename T> T strto(const char* s);
template <> inline int strto (const char* s) { return std::atoi(s); }
template <> inline bool strto (const char* s) { return std::atoi(s); }
template <> inline double strto (const char* s) { return std::atof(s); }
template <> inline std::string strto (const char* s) { return std::string(s); }

static homme::cslmpi::CslMpi::Ptr g_csl_mpi;

extern "C" {
void slmm_init_impl (
  homme::Int fcomm, homme::Int transport_alg, homme::Int np,
  homme::Int nlev, homme::Int qsize, homme::Int qsized, homme::Int nelem,
  homme::Int nelemd, homme::Int cubed_sphere_map,
  const homme::Int* lid2gid, const homme::Int* lid2facenum,
  const homme::Int* nbr_id_rank, const homme::Int* nirptr,
  homme::Int sl_nearest_point_lev, homme::Int, homme::Int, homme::Int,
  homme::Int)
{
  homme::slmm_init(np, nelem, nelemd, transport_alg, cubed_sphere_map,
                   sl_nearest_point_lev - 1, lid2facenum);
  slmm_throw_if(homme::g_advecter->is_cisl(),
                "CISL code was removed.");
  const auto p = homme::mpi::make_parallel(MPI_Comm_f2c(fcomm));
  g_csl_mpi = homme::cslmpi::init(homme::g_advecter, p, np, nlev, qsize,
                                  qsized, nelemd, nbr_id_rank, nirptr,
                                  2 /* halo */);
}

void slmm_query_bufsz (homme::Int* sendsz, homme::Int* recvsz) {
  slmm_assert(g_csl_mpi);
  homme::Int s = 0, r = 0;
  for (const auto e : g_csl_mpi->sendsz) s += e;
  for (const auto e : g_csl_mpi->recvsz) r += e;
  *sendsz = s;
  *recvsz = r;
}

void slmm_set_bufs (homme::Real* sendbuf, homme::Real* recvbuf,
                    homme::Int, homme::Int) {
  slmm_assert(g_csl_mpi);
  homme::cslmpi::alloc_mpi_buffers(*g_csl_mpi, sendbuf, recvbuf);
}

void slmm_get_mpi_pattern (homme::Int* sl_mpi) {
  *sl_mpi = g_csl_mpi ? 1 : 0;
}

void slmm_init_local_mesh (
  homme::Int ie, homme::Cartesian3D* neigh_corners, homme::Int nnc,
  homme::Cartesian3D* p_inside, homme::Int)
{
  homme::g_advecter->init_local_mesh_if_needed(
    ie - 1, homme::FA3<const homme::Real>(
      reinterpret_cast<const homme::Real*>(neigh_corners), 3, 4, nnc),
    reinterpret_cast<const homme::Real*>(p_inside));
}

void slmm_init_finalize () {
  if (g_csl_mpi)
    homme::cslmpi::finalize_local_meshes(*g_csl_mpi, *homme::g_advecter);
}

void slmm_check_ref2sphere (homme::Int ie, homme::Cartesian3D* p) {
  homme::g_advecter->check_ref2sphere(
    ie - 1, reinterpret_cast<const homme::Real*>(p));
}

void slmm_csl_set_elem_data (
  homme::Int ie, homme::Real* metdet, homme::Real* qdp, homme::Real* dp,
  homme::Real* q, homme::Int nelem_in_patch)
{
  slmm_assert(g_csl_mpi);
  homme::cslmpi::set_elem_data(*g_csl_mpi, ie - 1, metdet, qdp, dp, q,
                               nelem_in_patch);
}

void slmm_csl (
  homme::Int nets, homme::Int nete, homme::Cartesian3D* dep_points,
  homme::Real* minq, homme::Real* maxq, homme::Int* info)
{
  slmm_assert(g_csl_mpi);
  slmm_assert(g_csl_mpi->sendsz.empty()); // alloc_mpi_buffers was called
  *info = 0;
  try {
    homme::cslmpi::step<4>(*g_csl_mpi, nets - 1, nete - 1,
                           dep_points, minq, maxq);
  } catch (const std::exception& e) {
    std::cerr << e.what();
    *info = -1;
  }
}

void slmm_finalize () {
  g_csl_mpi = nullptr;
  homme::g_advecter = nullptr;
}
} // extern "C"
