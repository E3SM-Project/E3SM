#ifndef SCREAM_COMBINE_OPS_HPP
#define SCREAM_COMBINE_OPS_HPP

#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_scalar_traits.hpp>
#include <ekat_math_utils.hpp>
#include <ekat_pack_where.hpp>
#include <ekat_kernel_assert.hpp>

// For KOKKOS_INLINE_FUNCTION
#include <Kokkos_Core.hpp>
#include <type_traits>

namespace scream {

/*
 * Flags used to specify how to handle a variable update.
 * When computing f(x), there are a few ways to combine
 * the result with the value of the variable were we want
 * to store it:
 *    y = alpha*f(x) + beta*y
 *    y = alpha*f(x) * beta*y
 *    y = beta*y / (alpha*f(x))
 * This enum can be used as template arg in some general functions,
 * so that we can write a single f(x), and then combine:
 *    combine<CM>(f(x),y,alpha,beta);
 * The result has zero overhead compared to any specific version,
 * since the if/switch statements involving the combine mode CM
 * are compiled-away by the compiler
 */

enum class CombineMode {
  Replace,    // out = alpha*in
  Update,     // out = beta*out + alpha*in
  Multiply,   // out = (beta*out)*(alpha*in)
  Divide,     // out = (beta*out)/(alpha*in)
  Max,        // out = max(beta*out,alpha*in)
  Min         // out = min(beta*out,alpha*in)
};

// Small helper functions to combine a new value with an old one.
// The template argument help reducing the number of operations
// performed (the if is resolved at compile time). In the most
// complete form, the function performs (in functional programming notation):
//    result = (op beta*result alpha*newVal) (where op can be +, *, /, max, min)
// This routine should have no overhead compared to a manual
// update (assuming you call it with the proper CM)
// If fill-aware=true, we only perform the combine operation if either
//   a) CM is Replace
//   b) the new value is NOT fill_value
// In other words, with fill_aware=true we IGNORE new values that are equal to fill_value,
// unless we are replacing the content.
// NOTE: the default 'void' for the scalar types is obviously never used, since the
//       type is deduced from the inputs. The reason for the default is to allow to
//       give a default to fill_aware.
template<CombineMode CM, bool fill_aware = false, typename ScalarIn = void, typename ScalarOut = void,
         typename CoeffType = typename ekat::ScalarTraits<ScalarIn>::scalar_type>
KOKKOS_FORCEINLINE_FUNCTION
std::enable_if_t<not ekat::ScalarTraits<ScalarIn>::is_simd and
                 not ekat::ScalarTraits<ScalarOut>::is_simd>
combine (const ScalarIn& newVal, ScalarOut& result,
         const CoeffType alpha, const CoeffType beta)
{
  // Replace is simple, and does not care about fill_aware
  if constexpr (CM==CombineMode::Replace) {
    result = newVal;
    return;
  }

  if constexpr (fill_aware) {
    if (newVal!=constants::fill_value<ScalarIn>) {
      combine<CM,false>(newVal,result,alpha,beta);
    }
  } else {
    // Not a fill-aware combine
    using ekat::impl::max;
    using ekat::impl::min;
    switch (CM) {
      case CombineMode::Update:
        result *= beta;
        result += alpha*newVal;
        break;
      case CombineMode::Multiply:
        result *= (alpha*beta)*newVal;
        break;
      case CombineMode::Divide:
        result /= (alpha/beta) * newVal;
        break;
      case CombineMode::Max:
        result  = max(beta*result,static_cast<const ScalarOut&>(alpha*newVal));
        break;
      case CombineMode::Min:
        result  = min(beta*result,static_cast<const ScalarOut&>(alpha*newVal));
        break;
      default:
        EKAT_KERNEL_ASSERT ("Unsupported/unexpected combine mode.\n");
    }
  }
}

// Special version of combine for pack types
template<CombineMode CM, bool fill_aware = false, typename ScalarIn = void, typename ScalarOut = void,
         typename CoeffType = typename ekat::ScalarTraits<ScalarIn>::scalar_type>
KOKKOS_FORCEINLINE_FUNCTION
std::enable_if_t<ekat::ScalarTraits<ScalarIn>::is_simd or
                 ekat::ScalarTraits<ScalarOut>::is_simd>
combine (const ScalarIn& newVal, ScalarOut& result,
         const CoeffType alpha, const CoeffType beta)
{
  if constexpr (CM==CombineMode::Replace) {
    result = newVal;
    return;
  }
  if constexpr (fill_aware) {
    auto where = ekat::where(newVal!=constants::fill_value<ScalarIn>,result);

    if (where.any()) {
      auto tmp = result;
      combine<CM,false>(newVal,tmp,alpha,beta);
      where = tmp;
    }
  }
}

} // namespace scream

#endif // SCREAM_COMBINE_OPS_HPP
