#ifndef INCLUDE_SCREAM_PACK
#define INCLUDE_SCREAM_PACK

//TODO
// - bounds checking define

#include "scream_types.hpp"
#include "util/scream_utils.hpp"
#include "scream_macros.hpp"

#include <type_traits>

namespace scream {
namespace pack {

/* API for using "packed" data in scream. Packs are just bundles of N
   scalars wthin a single object. Using packed data makes it much easier
   to get good vectorization with C++.

   Pack is a vectorization pack, and Mask is a conditional mask for Pack::set
   constructed from operators among packs.

   If pack size (PACKN) is 1, then a Pack behaves as a scalar, and a Mask
   behaves roughly as a bool. Mask purposely does not support 'operator bool'
   because it is ambiguous whether operator bool should act as any() or all(),
   so we want the caller to be explicit.
 */

template <int PACKN>
struct Mask {
  // One tends to think a short boolean type would be useful here (e.g., bool or
  // char), but that is bad for vectorization. int or long are best.
  typedef long type;

  // A tag for this struct for type checking.
  enum { masktag = true };
  // Pack and Mask sizes are the same, n.
  enum { n = PACKN };

  KOKKOS_FORCEINLINE_FUNCTION explicit Mask () {}

  // Init all slots of the Mask to 'init'.
  KOKKOS_FORCEINLINE_FUNCTION explicit Mask (const bool& init) {
    //vector_simd // Intel 18 is having an issue with this loop.
    vector_disabled for (int i = 0; i < n; ++i) d[i] = init;
  }

  // Set slot i to val.
  KOKKOS_FORCEINLINE_FUNCTION void set (const int& i, const bool& val) { d[i] = val; }
  // Get slot i.
  KOKKOS_FORCEINLINE_FUNCTION bool operator[] (const int& i) const { return d[i]; }

  // Is any slot true?
  KOKKOS_FORCEINLINE_FUNCTION bool any () const {
    bool b = false;
    vector_simd for (int i = 0; i < n; ++i) if (d[i]) b = true;
    return b;
  }

  // Are all slots true?
  KOKKOS_FORCEINLINE_FUNCTION bool all () const {
    bool b = true;
    vector_simd for (int i = 0; i < n; ++i) if ( ! d[i]) b = false;
    return b;
  }

private:
  type d[n];
};

// Use enable_if and masktag so that we can template on 'Mask' and yet not have
// our operator overloads, in particular, be used for something other than the
// Mask type.
template <typename Mask>
using OnlyMask = typename std::enable_if<Mask::masktag,Mask>::type;
template <typename Mask, typename Return>
using OnlyMaskReturn = typename std::enable_if<Mask::masktag,Return>::type;

// Codify how a user can construct their own loops conditioned on mask slot
// values.
#define scream_masked_loop(mask, s)                         \
  vector_simd for (int s = 0; s < mask.n; ++s) if (mask[s])

#define scream_masked_loop_no_force_vec(mask, s)              \
  vector_ivdep for (int s = 0; s < mask.n; ++s) if (mask[s])

#define scream_masked_loop_no_vec(mask, s)                    \
  vector_novec for (int s = 0; s < mask.n; ++s) if (mask[s])

// Implementation detail for generating binary ops for mask op mask.
#define scream_mask_gen_bin_op_mm(op, impl)                   \
  template <typename Mask> KOKKOS_INLINE_FUNCTION             \
  OnlyMask<Mask> operator op (const Mask& a, const Mask& b) { \
    Mask m(false);                                            \
    vector_simd for (int i = 0; i < Mask::n; ++i)             \
      if (a[i] impl b[i]) m.set(i, true);                     \
    return m;                                                 \
  }

// Implementation detail for generating binary ops for mask op bool.
#define scream_mask_gen_bin_op_mb(op, impl)                   \
  template <typename Mask> KOKKOS_INLINE_FUNCTION             \
  OnlyMask<Mask> operator op (const Mask& a, const bool b) {  \
    Mask m(false);                                            \
    vector_simd for (int i = 0; i < Mask::n; ++i)             \
      if (a[i] impl b) m.set(i, true);                        \
    return m;                                                 \
  }

scream_mask_gen_bin_op_mm(&&, &&)
scream_mask_gen_bin_op_mm(||, ||)
scream_mask_gen_bin_op_mb(&&, &&)
scream_mask_gen_bin_op_mb(||, ||)

// Negate the mask.
template <typename Mask> KOKKOS_INLINE_FUNCTION
OnlyMask<Mask> operator ! (const Mask& m) {
  Mask nm(false);
  vector_simd for (int i = 0; i < Mask::n; ++i) nm.set(i, ! m[i]);
  return nm;
}

// Implementation detail for generating Pack assignment operators. _p means the
// input is a Pack; _s means the input is a scalar.
#define scream_pack_gen_assign_op_p(op)                   \
  KOKKOS_FORCEINLINE_FUNCTION                             \
  Pack& operator op (const Pack& a) {                     \
    vector_simd for (int i = 0; i < n; ++i) d[i] op a[i]; \
    return *this;                                         \
  }
#define scream_pack_gen_assign_op_s(op)                 \
  KOKKOS_FORCEINLINE_FUNCTION                           \
  Pack& operator op (const scalar& a) {                 \
    vector_simd for (int i = 0; i < n; ++i) d[i] op a;  \
    return *this;                                       \
  }
#define scream_pack_gen_assign_op_all(op)       \
  scream_pack_gen_assign_op_p(op)               \
  scream_pack_gen_assign_op_s(op)

// The Pack type. Mask was defined first since it's used in Pack.
template <typename SCALAR, int PACKN>
struct Pack {
  // Pack's tag for type checking.
  enum { packtag = true };
  // Number of slots in the Pack.
  enum { n = PACKN };

  // A definitely non-const version of the input scalar type.
  typedef typename std::remove_const<SCALAR>::type scalar;

  KOKKOS_FORCEINLINE_FUNCTION explicit Pack () {
#ifndef KOKKOS_ENABLE_CUDA
    // Quiet NaNs don't work on Cuda.
    vector_simd for (int i = 0; i < n; ++i)
      d[i] = std::numeric_limits<scalar>::quiet_NaN();
#endif
  }

  // Init all slots to scalar v.
  KOKKOS_FORCEINLINE_FUNCTION explicit Pack (const scalar& v) {
    vector_simd for (int i = 0; i < n; ++i) d[i] = v;
  }

  // Init this Pack from another one.
  template <typename PackIn> KOKKOS_FORCEINLINE_FUNCTION explicit
  Pack (const PackIn& v, typename std::enable_if<PackIn::packtag>::type* = nullptr) {
    static_assert(static_cast<int>(PackIn::n) == static_cast<int>(n),
                  "Pack::n must be the same.");
    vector_simd for (int i = 0; i < n; ++i) d[i] = v[i];
  }

  // Init this Pack from another one, but only where Mask is true; otherwise
  // init to default value.
  template <typename PackIn> KOKKOS_FORCEINLINE_FUNCTION explicit
  Pack (const Mask<Pack::n>& m, const PackIn& p) {
    static_assert(static_cast<int>(PackIn::n) == static_cast<int>(n),
                  "Pack::n must be the same.");
#ifndef KOKKOS_ENABLE_CUDA
    vector_simd for (int i = 0; i < n; ++i)
      d[i] = m[i] ? p[i] : std::numeric_limits<scalar>::quiet_NaN();
#else
    vector_simd for (int i = 0; i < n; ++i)
      d[i] = m[i] ? p[i] : 0;
#endif
  }

  KOKKOS_FORCEINLINE_FUNCTION const scalar& operator[] (const int& i) const { return d[i]; }
  KOKKOS_FORCEINLINE_FUNCTION scalar& operator[] (const int& i) { return d[i]; }

  scream_pack_gen_assign_op_all(=)
  scream_pack_gen_assign_op_all(+=)
  scream_pack_gen_assign_op_all(-=)
  scream_pack_gen_assign_op_all(*=)
  scream_pack_gen_assign_op_all(/=)

  KOKKOS_FORCEINLINE_FUNCTION
  void set (const Mask<n>& mask, const scalar& v) {
    vector_simd for (int i = 0; i < n; ++i) if (mask[i]) d[i] = v;
  }
  template <typename PackIn> KOKKOS_FORCEINLINE_FUNCTION
  void set (const Mask<n>& mask, const PackIn& p,
            typename std::enable_if<PackIn::packtag>::type* = nullptr) {
    static_assert(static_cast<int>(PackIn::n) == static_cast<int>(n),
                  "Pack::n must be the same.");
    vector_simd for (int i = 0; i < n; ++i) if (mask[i]) d[i] = p[i];
  }

private:
  scalar d[n];
};

// Use enable_if and packtag so that we can template on 'Pack' and yet not have
// our operator overloads, in particular, be used for something other than the
// Pack type.
template <typename Pack>
using OnlyPack = typename std::enable_if<Pack::packtag,Pack>::type;
template <typename Pack, typename Return>
using OnlyPackReturn = typename std::enable_if<Pack::packtag,Return>::type;

// Later, we might support type promotion. For now, caller must explicitly
// promote a pack's scalar type in mixed-type arithmetic.

#define scream_pack_gen_bin_op_pp(op)                                   \
  template <typename Pack> KOKKOS_FORCEINLINE_FUNCTION                  \
  OnlyPack<Pack> operator op (const Pack& a, const Pack& b) {           \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a[i] op b[i];  \
    return c;                                                           \
  }
#define scream_pack_gen_bin_op_ps(op)                                   \
  template <typename Pack, typename Scalar> KOKKOS_FORCEINLINE_FUNCTION \
  OnlyPack<Pack> operator op (const Pack& a, const Scalar& b) {         \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a[i] op b;     \
    return c;                                                           \
  }
#define scream_pack_gen_bin_op_sp(op)                                   \
  template <typename Pack, typename Scalar> KOKKOS_FORCEINLINE_FUNCTION \
  OnlyPack<Pack> operator op (const Scalar& a, const Pack& b) {         \
    Pack c;                                                             \
    vector_simd for (int i = 0; i < Pack::n; ++i) c[i] = a op b[i];     \
    return c;                                                           \
  }
#define scream_pack_gen_bin_op_all(op)          \
  scream_pack_gen_bin_op_pp(op)                 \
  scream_pack_gen_bin_op_ps(op)                 \
  scream_pack_gen_bin_op_sp(op)

scream_pack_gen_bin_op_all(+)
scream_pack_gen_bin_op_all(-)
scream_pack_gen_bin_op_all(*)
scream_pack_gen_bin_op_all(/)

#define scream_pack_gen_unary_fn(fn, impl)                            \
  template <typename Pack> KOKKOS_INLINE_FUNCTION                     \
  OnlyPack<Pack> fn (const Pack& p) {                                 \
    Pack s;                                                           \
    vector_simd for (int i = 0; i < Pack::n; ++i) s[i] = impl(p[i]);  \
    return s;                                                         \
  }
#define scream_pack_gen_unary_stdfn(fn) scream_pack_gen_unary_fn(fn, std::fn)
scream_pack_gen_unary_stdfn(abs)
scream_pack_gen_unary_stdfn(exp)
scream_pack_gen_unary_stdfn(log)
scream_pack_gen_unary_stdfn(log10)
scream_pack_gen_unary_stdfn(tgamma)
scream_pack_gen_unary_stdfn(cbrt)


template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<Pack, typename Pack::scalar> min (const Pack& p) {
  typename Pack::scalar v(p[0]);
  vector_disabled for (int i = 0; i < Pack::n; ++i) v = util::min(v, p[i]);
  return v;
}

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<Pack, typename Pack::scalar> max (const Pack& p) {
  typename Pack::scalar v(p[0]);
  vector_simd for (int i = 0; i < Pack::n; ++i) v = util::max(v, p[i]);
  return v;
}

// min(init, min(p(mask)))
template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<Pack, typename Pack::scalar>
min (const Mask<Pack::n>& mask, typename Pack::scalar init, const Pack& p) {
  vector_disabled for (int i = 0; i < Pack::n; ++i)
    if (mask[i]) init = util::min(init, p[i]);
  return init;
}

// max(init, max(p(mask)))
template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<Pack, typename Pack::scalar>
max (const Mask<Pack::n>& mask, typename Pack::scalar init, const Pack& p) {
  vector_simd for (int i = 0; i < Pack::n; ++i)
    if (mask[i]) init = util::max(init, p[i]);
  return init;
}

#define scream_pack_gen_bin_fn_pp(fn, impl)           \
  template <typename Pack> KOKKOS_INLINE_FUNCTION     \
  OnlyPack<Pack> fn (const Pack& a, const Pack& b) {  \
    Pack s;                                           \
    vector_simd for (int i = 0; i < Pack::n; ++i)     \
      s[i] = impl(a[i], b[i]);                        \
    return s;                                         \
  }
#define scream_pack_gen_bin_fn_ps(fn, impl)                         \
  template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION  \
  OnlyPack<Pack> fn (const Pack& a, const Scalar& b) {              \
    Pack s;                                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      s[i] = impl<typename Pack::scalar>(a[i], b);                  \
    return s;                                                       \
  }
#define scream_pack_gen_bin_fn_sp(fn, impl)                         \
  template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION  \
  OnlyPack<Pack> fn (const Scalar& a, const Pack& b) {              \
    Pack s;                                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      s[i] = impl<typename Pack::scalar>(a, b[i]);                  \
    return s;                                                       \
  }
#define scream_pack_gen_bin_fn_all(fn, impl)    \
  scream_pack_gen_bin_fn_pp(fn, impl)           \
  scream_pack_gen_bin_fn_ps(fn, impl)           \
  scream_pack_gen_bin_fn_sp(fn, impl)

scream_pack_gen_bin_fn_all(min, util::min)
scream_pack_gen_bin_fn_all(max, util::max)

// On Intel 17 for KNL, I'm getting a ~1-ulp diff on const Scalar& b. I don't
// understand its source. But, in any case, I'm writing a separate impl here to
// get around that.
//scream_pack_gen_bin_fn_all(pow, std::pow)
template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> pow (const Pack& a, const Scalar/*&*/ b,
                    typename std::enable_if<std::is_floating_point<Scalar>::value>::type* = 0 ) {
  Pack s;
  vector_simd for (int i = 0; i < Pack::n; ++i)
    s[i] = std::pow<typename Pack::scalar>(a[i], b);
  return s;
}

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> pow(const Pack& a, const int b) {
  Pack s;
  vector_simd for (int i = 0; i < Pack::n; ++i)
    s[i] = util::ipow_scalar(a[i], b);
  return s;
}

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> shift_right (const Pack& pm1, const Pack& p) {
  Pack s;
  s[0] = pm1[Pack::n-1];
  vector_simd for (int i = 1; i < Pack::n; ++i) s[i] = p[i-1];
  return s;
}

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> shift_right (const typename Pack::scalar& pm1, const Pack& p) {
  Pack s;
  s[0] = pm1;
  vector_simd for (int i = 1; i < Pack::n; ++i) s[i] = p[i-1];
  return s;
}

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> shift_left (const Pack& pp1, const Pack& p) {
  Pack s;
  s[Pack::n-1] = pp1[0];
  vector_simd for (int i = 0; i < Pack::n-1; ++i) s[i] = p[i+1];
  return s;
}

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> shift_left (const typename Pack::scalar& pp1, const Pack& p) {
  Pack s;
  s[Pack::n-1] = pp1;
  vector_simd for (int i = 0; i < Pack::n-1; ++i) s[i] = p[i+1];
  return s;
}

#define scream_mask_gen_bin_op_pp(op)             \
  template <typename Pack> KOKKOS_INLINE_FUNCTION \
  OnlyPackReturn<Pack, Mask<Pack::n> >            \
  operator op (const Pack& a, const Pack& b) {    \
    Mask<Pack::n> m(false);                       \
    vector_simd for (int i = 0; i < Pack::n; ++i) \
      if (a[i] op b[i]) m.set(i, true);           \
    return m;                                     \
  }
#define scream_mask_gen_bin_op_ps(op)                               \
  template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION  \
  OnlyPackReturn<Pack, Mask<Pack::n> >                              \
  operator op (const Pack& a, const Scalar& b) {                    \
    Mask<Pack::n> m(false);                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      if (a[i] op b) m.set(i, true);                                \
    return m;                                                       \
  }
#define scream_mask_gen_bin_op_sp(op)                               \
  template <typename Pack, typename Scalar> KOKKOS_INLINE_FUNCTION  \
  OnlyPackReturn<Pack, Mask<Pack::n> >                              \
  operator op (const Scalar& a, const Pack& b) {                    \
    Mask<Pack::n> m(false);                                         \
    vector_simd for (int i = 0; i < Pack::n; ++i)                   \
      if (a op b[i]) m.set(i, true);                                \
    return m;                                                       \
  }
#define scream_mask_gen_bin_op_all(op)          \
  scream_mask_gen_bin_op_pp(op)                 \
  scream_mask_gen_bin_op_ps(op)                 \
  scream_mask_gen_bin_op_sp(op)

scream_mask_gen_bin_op_all(==)
scream_mask_gen_bin_op_all(!=)
scream_mask_gen_bin_op_all(>=)
scream_mask_gen_bin_op_all(<=)
scream_mask_gen_bin_op_all(>)
scream_mask_gen_bin_op_all(<)

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPackReturn<Pack,Int> npack(const Int& nscalar) {
  return (nscalar + Pack::n - 1) / Pack::n;
}

template <typename Pack> KOKKOS_INLINE_FUNCTION
OnlyPack<Pack> range (const typename Pack::scalar& start) {
  typedef Pack pack;
  pack p;
  vector_simd for (int i = 0; i < Pack::n; ++i) p[i] = start + i;
  return p;
}

#undef scream_pack_gen_assign_op_p
#undef scream_pack_gen_assign_op_s
#undef scream_pack_gen_assign_op_all
#undef scream_pack_gen_bin_op_pp
#undef scream_pack_gen_bin_op_ps
#undef scream_pack_gen_bin_op_sp
#undef scream_pack_gen_bin_op_all
#undef scream_pack_gen_unary_fn
#undef scream_pack_gen_unary_stdfn
#undef scream_pack_gen_bin_fn_pp
#undef scream_pack_gen_bin_fn_ps
#undef scream_pack_gen_bin_fn_sp
#undef scream_pack_gen_bin_fn_all
#undef scream_mask_gen_bin_op_pp
#undef scream_mask_gen_bin_op_ps
#undef scream_mask_gen_bin_op_sp
#undef scream_mask_gen_bin_op_all
#undef scream_mask_gen_bin_op_mm
#undef scream_mask_gen_bin_op_mb

} // namespace pack

namespace util {

// Specialization of scalar-properties detection helper struct for Pack types

template<typename T, int N>
struct ScalarProperties<pack::Pack<T,N>> {
  // This seems funky. But write down a pow of 2 and a non-pow of 2 in binary (both positive), and you'll see why it works
  static_assert (N>0 && ((N & (N-1))==0), "Error! Packs can only have power of two lengths.\n");

  using scalar_type = typename ScalarProperties<T>::scalar_type;
  static constexpr bool is_pack = true;
};

// Specialization of TypeName struct for Pack type
template<typename T, int N>
struct TypeName<pack::Pack<T,N>> {
  static std::string name () {
    return "Pack<" + TypeName<T>::name() + "," + std::to_string(N) + ">";
  }
};

} // namespace util

} // namespace scream

#endif // INCLUDE_SCREAM_PACK
