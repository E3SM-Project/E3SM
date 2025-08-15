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

template<CombineMode CM, typename ScalarIn, typename ScalarOut,
         typename CoeffType = typename ekat::ScalarTraits<ScalarIn>::scalar_type>
KOKKOS_FORCEINLINE_FUNCTION
void combine (const ScalarIn& newVal, ScalarOut& result,
              const CoeffType alpha, const CoeffType beta)
{
  using ekat::impl::max;
  using ekat::impl::min;
  switch (CM) {
    case CombineMode::Replace:
      result = alpha*newVal;
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
      result  = max(beta*result,static_cast<const ScalarOut&>(alpha*newVal));
      break;
    case CombineMode::Min:
      result  = min(beta*result,static_cast<const ScalarOut&>(alpha*newVal));
      break;
  }
}

// Special version of combine that ignores newVal if newVal==fill_value.
// Replace is the only combine mode that is allowed to consider fill_val values.
// This is b/c it's the only way we can use this function inside Field method/utils
// in order to set all entries of a Field to fill_val. You can also think of fill_val
// as a special number for which the arithmetic operations are not defined.
// All CM except for Replace involve an arithmetic op between of two numbers,
// so combining with fill_val makes no sense. However, it makes sense to set
// an output variable to fill_val.
template<CombineMode CM, typename ScalarIn, typename ScalarOut,
         typename CoeffType = typename ekat::ScalarTraits<ScalarIn>::scalar_type>
KOKKOS_FORCEINLINE_FUNCTION
std::enable_if_t<ekat::ScalarTraits<ScalarIn>::is_simd or
                 ekat::ScalarTraits<ScalarOut>::is_simd>
fill_aware_combine (const ScalarIn& newVal, ScalarOut& result,
                    const typename ekat::ScalarTraits<ScalarIn>::scalar_type fill_val,
                    const CoeffType alpha, const CoeffType beta)
{
  if constexpr (CM==CombineMode::Replace) {
    return combine<CM>(newVal,result,alpha,beta);
  }

  // The where object will perform the assignment ONLY where the mask is true
  auto where = ekat::where(newVal!=fill_val,result);
  if (where.any()) {
    // TODO: I thought about doing the switch manually, and do stuff like (e.g., for Update)
    //  where *= beta;
    //  where += alpha*newVal
    // but there is no packed version of where.max(rhs), only a scalar version
    // (meaning a version where rhs is a scalar, not a pack).
    // If ekat::where_expression ever implements a packed overload for max/min,
    // we can get rid of the temporary by doing a manual switch.
    auto tmp = result;
    combine<CM>(newVal,tmp,alpha,beta);
    where = tmp;
  }
}

template<CombineMode CM, typename ScalarIn, typename ScalarOut,
         typename CoeffType = typename ekat::ScalarTraits<ScalarIn>::scalar_type>
KOKKOS_FORCEINLINE_FUNCTION
std::enable_if_t<not ekat::ScalarTraits<ScalarIn>::is_simd and
                 not ekat::ScalarTraits<ScalarOut>::is_simd>
fill_aware_combine (const ScalarIn& newVal, ScalarOut& result,
                    const typename ekat::ScalarTraits<ScalarIn>::scalar_type fill_val,
                    const CoeffType alpha, const CoeffType beta)
{
  if constexpr (CM==CombineMode::Replace) {
    return combine<CM>(newVal,result,alpha,beta);
  }

  if (newVal!=fill_val) {
    combine<CM>(newVal,result,alpha,beta);
  }
}

} // namespace scream

#endif // SCREAM_COMBINE_OPS_HPP
