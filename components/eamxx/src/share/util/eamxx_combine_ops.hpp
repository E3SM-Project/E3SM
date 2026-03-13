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

// How to handle FillValue values in combine operations.
// NOTE: CombineMode::Replace will ALWAYS set result=newVal,
//       regardless of fill value handling
enum FillValueHandling : int {
  None = 0,   // treat just like any other number
  IgnoreRhs,  // ignore rhs when equal to FV
  Absorbing   // FV op X = FV, and X op FV = FV, for all X
};

namespace impl {

template<CombineMode CM, typename ScalarIn, typename ScalarOut, typename CoeffType>
KOKKOS_FORCEINLINE_FUNCTION
void combine_fvh_none (const ScalarIn& newVal, ScalarOut& result,
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

template<CombineMode CM, typename ScalarIn, typename ScalarOut, typename CoeffType>
KOKKOS_FORCEINLINE_FUNCTION
void combine_fvh_absorb (const ScalarIn& newVal, ScalarOut& result,
                         const CoeffType alpha, const CoeffType beta)
{
  using ekat::impl::max;
  using ekat::impl::min;

  using inner_type = typename ekat::ScalarTraits<ScalarIn>::scalar_type;
  constexpr auto fill_val = constants::fill_value<inner_type>;

  auto where_fv = ekat::where(newVal==fill_val or result==fill_val,result);
  if (not where_fv.all()) {
    combine_fvh_none<CM>(newVal,result,alpha,beta);
  }
  if (where_fv.any()) {
    where_fv = fill_val;
  }
}

template<CombineMode CM, typename ScalarIn, typename ScalarOut, typename CoeffType>
KOKKOS_FORCEINLINE_FUNCTION
void combine_fvh_ignore (const ScalarIn& newVal, ScalarOut& result,
                         const CoeffType alpha, const CoeffType beta)
{
  using inner_type = typename ekat::ScalarTraits<ScalarIn>::scalar_type;
  constexpr auto fill_val = constants::fill_value<inner_type>;

  auto where_ok = ekat::where(newVal!=fill_val, result);
  if (where_ok.any()) {
    // Compute the combined result in a temporary, then write back only
    // where newVal is not fill_val, leaving fill_val positions unchanged.
    auto tmp = result;
    combine_fvh_none<CM>(newVal, tmp, alpha, beta);
    where_ok = tmp;
  }
}

} // namespace impl

// This is the function that user will call. It dispatches to one of the three
// fill-value handling implementations based on the runtime fvh argument.
// NOTE: CombineMode::Replace always sets result=newVal regardless of fvh.
// NOTE: the default 'void' for the scalar types is obviously never used, since the
//       type is deduced from the inputs. The reason for the default is to allow to
//       give a default to fvh.
template<CombineMode CM, typename ScalarIn = void, typename ScalarOut = void,
         typename CoeffType = typename ekat::ScalarTraits<ScalarIn>::scalar_type>
KOKKOS_FORCEINLINE_FUNCTION
void combine (const ScalarIn& newVal, ScalarOut& result,
              const CoeffType alpha, const CoeffType beta,
              FillValueHandling fvh = FillValueHandling::None)
{
  // Replace always ignores fill value handling: result = newVal unconditionally
  if constexpr (CM == CombineMode::Replace) {
    return impl::combine_fvh_none<CM>(newVal,result,alpha,beta);
  }

  switch (fvh) {
    case FillValueHandling::None:
      impl::combine_fvh_none<CM>(newVal,result,alpha,beta);
      break;
    case FillValueHandling::IgnoreRhs:
      impl::combine_fvh_ignore<CM>(newVal,result,alpha,beta);
      break;
    case FillValueHandling::Absorbing:
      impl::combine_fvh_absorb<CM>(newVal,result,alpha,beta);
      break;
    default:
      EKAT_KERNEL_ASSERT ("Unsupported/unexpected FillValueHandling.\n");
  }
}

} // namespace scream

#endif // SCREAM_COMBINE_OPS_HPP
