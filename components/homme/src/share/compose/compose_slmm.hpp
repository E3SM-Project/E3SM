#ifndef INCLUDE_COMPOSE_SLMM_HPP
#define INCLUDE_COMPOSE_SLMM_HPP

#include "compose_kokkos.hpp"
#include "compose_slmm_siqk.hpp"

namespace ko = Kokkos;

#ifndef NDEBUG
# define slmm_assert(condition) do {                                    \
    if ( ! (condition)) {                                               \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": FAIL:\n" << #condition  \
        << "\n";                                                        \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
# define slmm_kernel_assert(condition) do {     \
    if ( ! (condition))                         \
      Kokkos::abort(#condition);                \
  } while (0)
#else
# define slmm_assert(condition)
# define slmm_kernel_assert(condition)
#endif
#define slmm_throw_if(condition, message) do {                          \
    if (condition) {                                                    \
      std::stringstream _ss_;                                           \
      _ss_ << __FILE__ << ":" << __LINE__ << ": The condition:\n"       \
           << #condition "\nled to the exception\n" << message << "\n"; \
        throw std::logic_error(_ss_.str());                             \
    }                                                                   \
  } while (0)
#define slmm_kernel_throw_if(condition, message) do {               \
    if (condition)                                                  \
      Kokkos::abort(#condition " led to the exception\n" message);  \
  } while (0)

#define SLMM_KIF KOKKOS_INLINE_FUNCTION
#define SLMM_KF KOKKOS_FUNCTION

namespace slmm {
using siqk::Int;
using siqk::Real;
typedef Int Size;

#ifdef COMPOSE_TIMERS
struct Timer {
  Timer (const std::string& name_) : name("SLMM_isl_" + name_) { GPTLstart(name.c_str()); }
  ~Timer () { ko::fence(); GPTLstop(name.c_str()); }
private:
  const std::string name;
};
#else
struct Timer {
  Timer (const std::string&) {}
};
#endif

// A 2D array A can be thought of as having nslices(A) rows and szslice(A)
// columns. A slice can be obtained by
//     auto ak = slice(A, k);
// We use this format for arrays of vertices and adjacency arrays, for
// example. In most or all cases, the intention is to parallelize over slices,
// so a Kokkos operator() will do work on a particular slice.
using siqk::nslices;
using siqk::szslice;
using siqk::slice;

template<typename V> SLMM_KIF Int len (const V& v)
{ return static_cast<Int>(v.extent_int(0)); }

template<typename T> inline Int len (const std::vector<T>& v)
{ return static_cast<Int>(v.size()); }

struct Geometry {
  enum Type { sphere = 0, plane = 1 };
};

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

} // namespace slmm

#endif
