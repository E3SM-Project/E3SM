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

namespace impl {

// Small helper functions to combine a new value with an old one.
// The template argument help reducing the number of operations
// performed (the if is resolved at compile time). In the most
// complete form, the function performs (in functional programming notation):
//    result = (op beta*result alpha*newVal) (where op can be +, *, /, max, min)
// This routine should have no overhead compared to a manual
// update (assuming you call it with the proper CM)
template<CombineMode CM, typename ScalarIn, typename ScalarOut, typename CoeffType>
KOKKOS_FORCEINLINE_FUNCTION
void combine (const ScalarIn& newVal, ScalarOut& result,
              const CoeffType alpha, const CoeffType beta)
{
  using ekat::impl::max;
  using ekat::impl::min;
  switch (CM) {
    case CombineMode::Replace:
      result = newVal;
      break;
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
      result = max(beta*result,static_cast<const ScalarOut&>(alpha*newVal));
      break;
    case CombineMode::Min:
      result = min(beta*result,static_cast<const ScalarOut&>(alpha*newVal));
      break;
    default:
      EKAT_KERNEL_ASSERT ("Unsupported/unexpected combine mode.\n");
  }
}

} // namespace impl

// This is the function that user will call, which uses the one above internally
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
void combine (const ScalarIn& newVal, ScalarOut& result,
              const CoeffType alpha, const CoeffType beta)
{
  // If not fill-aware, or if CM==Replace, we don't need to check newValue
  if constexpr (not fill_aware or CM==CombineMode::Replace) {
    return impl::combine<CM>(newVal,result,alpha,beta);
  }

  // For the non-simd type case, we can avoid ekat::where, and simply check newVal against fill_value
  if constexpr (ekat::ScalarTraits<ScalarIn>::is_simd) {
    using inner_type = typename ekat::ScalarTraits<ScalarIn>::scalar_type;
    constexpr auto fill_val = constants::fill_value<inner_type>;
    auto where = ekat::where(newVal!=fill_val,result);
    if (where.any()) {
      auto tmp = result;
      impl::combine<CM>(newVal,tmp,alpha,beta);
      where = tmp;
    }
  } else {
    if (newVal!=constants::fill_value<ScalarIn>) {
      impl::combine<CM>(newVal,result,alpha,beta);
    }
  }
}

} // namespace scream

#endif // SCREAM_COMBINE_OPS_HPP
