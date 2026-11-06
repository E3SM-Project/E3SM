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

// Runtime-CM version of the above, which also adds a constant shift gamma (only used
// for CM=Update): result = beta*result + alpha*newVal + gamma.
// The switch is on a value that is uniform across the whole kernel, so the cost is
// negligible for memory-bound kernels, while it avoids instantiating the (expensive)
// kernel launch machinery once per CM.
template<typename ScalarIn, typename ScalarOut, typename CoeffType>
KOKKOS_FORCEINLINE_FUNCTION
void combine (const CombineMode cm, const ScalarIn& newVal, ScalarOut& result,
              const CoeffType alpha, const CoeffType beta, const CoeffType gamma)
{
  switch (cm) {
    case CombineMode::Replace:  combine<CombineMode::Replace> (newVal,result,alpha,beta); break;
    case CombineMode::Update:   combine<CombineMode::Update>  (newVal,result,alpha,beta);
                                result += gamma;                                          break;
    case CombineMode::Multiply: combine<CombineMode::Multiply>(newVal,result,alpha,beta); break;
    case CombineMode::Divide:   combine<CombineMode::Divide>  (newVal,result,alpha,beta); break;
    case CombineMode::Max:      combine<CombineMode::Max>     (newVal,result,alpha,beta); break;
    case CombineMode::Min:      combine<CombineMode::Min>     (newVal,result,alpha,beta); break;
    default:
      EKAT_KERNEL_ASSERT ("Unsupported/unexpected combine mode.\n");
  }
}

// This is similar to the above one, but only does something if newVal is NOT equal to fill_value.
// The only exception is the case CM=Replace, in which case we ALWAYS perform the copy
template<CombineMode CM, typename ScalarIn, typename ScalarOut,
         typename CoeffType = typename ekat::ScalarTraits<ScalarIn>::scalar_type>
KOKKOS_FORCEINLINE_FUNCTION
void combine_fill_aware (const ScalarIn& newVal, ScalarOut& result,
                         const CoeffType alpha, const CoeffType beta)
{
  // If CM==Replace, we don't check for FV. We ALWAYS replace
  if constexpr (CM==CombineMode::Replace) {
    return combine<CM>(newVal,result,alpha,beta);
  }

  // For the non-simd type case, we can avoid ekat::where, and simply check newVal against fill_value
  if constexpr (ekat::ScalarTraits<ScalarIn>::is_simd) {
    using inner_type = typename ekat::ScalarTraits<ScalarIn>::scalar_type;
    constexpr auto fill_val = constants::fill_value<inner_type>;
    auto where = ekat::where(newVal!=fill_val,result);
    if (where.any()) {
      auto tmp = result;
      combine<CM>(newVal,tmp,alpha,beta);
      where = tmp;
    }
  } else {
    if (newVal!=constants::fill_value<ScalarIn>) {
      combine<CM>(newVal,result,alpha,beta);
    }
  }
}

// Runtime-CM version of combine_fill_aware, with the extra shift gamma (CM=Update only).
// NOTE: gamma is applied to result regardless of the value of newVal, unless gamma
//       itself is the fill value (in which case it is skipped).
template<typename ScalarIn, typename ScalarOut,
         typename CoeffType = typename ekat::ScalarTraits<ScalarIn>::scalar_type>
KOKKOS_FORCEINLINE_FUNCTION
void combine_fill_aware (const CombineMode cm, const ScalarIn& newVal, ScalarOut& result,
                         const CoeffType alpha, const CoeffType beta, const CoeffType gamma)
{
  switch (cm) {
    case CombineMode::Replace:  combine_fill_aware<CombineMode::Replace> (newVal,result,alpha,beta); break;
    case CombineMode::Update:   combine_fill_aware<CombineMode::Update>  (newVal,result,alpha,beta);
                                if (gamma!=constants::fill_value<CoeffType>)
                                  result += gamma;                                                 break;
    case CombineMode::Multiply: combine_fill_aware<CombineMode::Multiply>(newVal,result,alpha,beta); break;
    case CombineMode::Divide:   combine_fill_aware<CombineMode::Divide>  (newVal,result,alpha,beta); break;
    case CombineMode::Max:      combine_fill_aware<CombineMode::Max>     (newVal,result,alpha,beta); break;
    case CombineMode::Min:      combine_fill_aware<CombineMode::Min>     (newVal,result,alpha,beta); break;
    default:
      EKAT_KERNEL_ASSERT ("Unsupported/unexpected combine mode.\n");
  }
}

} // namespace scream

#endif // SCREAM_COMBINE_OPS_HPP
