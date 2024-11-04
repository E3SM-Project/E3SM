#ifndef SCREAM_COMBINE_OPS_HPP
#define SCREAM_COMBINE_OPS_HPP

#include "share/util/scream_universal_constants.hpp"

// For KOKKOS_INLINE_FUNCTION
#include <Kokkos_Core.hpp>
#include <type_traits>
#include "ekat/ekat_scalar_traits.hpp"

namespace scream {

/*
 * Flags used to specify how to handle a variable update.
 * When computing f(x), there are a few ways to combine
 * the result with the value of the variable were we want
 * to store it:
 *    y = alpha*f(x) + beta*y
 *    y = y*f(x)
 *    y = y/f(x)
 * This enum can be used as template arg in some general functions,
 * so that we can write a single f(x), and then combine:
 *    combine<CM>(f(x),y,alpha,beta);
 * The result has zero overhead compared to any specific version,
 * since the if/switch statements involving the combine mode CM
 * are compiled-away by the compiler
 */

enum class CombineMode {
  ScaleUpdate,  // out = beta*out + alpha*in (most generic case)
  Update,       // out = beta*out + in (special case of ScaleUpdate wiht alpha=1)
  ScaleAdd,     // out = out + alpha*in (special case of ScaleUpdate with beta=1)
  ScaleReplace, // out = alpha*in (special case of ScaleUpdate with beta=0)
  Add,          // out = out + in (special case of ScaleUpdate with alpha=beta=1)
  Rescale,      // out = beta*out
  Replace,      // out = in
  Multiply,     // out = out*in
  Divide        // out = out/in
};

// Functions mostly used for debug purposes. They check whether a combine mode
// requires valid alpha or beta coefficients.
template<CombineMode CM>
KOKKOS_INLINE_FUNCTION
static constexpr bool needsAlpha () {
  return CM==CombineMode::ScaleReplace || CM==CombineMode::ScaleAdd || CM==CombineMode::ScaleUpdate;
}

template<CombineMode CM>
KOKKOS_INLINE_FUNCTION
static constexpr bool needsBeta () {
  return CM==CombineMode::Update || CM==CombineMode::ScaleUpdate || CM==CombineMode::Rescale;
}


// Small helper functions to combine a new value with an old one.
// The template argument help reducing the number of operations
// performed (the if is resolved at compile time). In the most
// complete form, the function performs
//    result = beta*result + alpha*newVal
// This routine should have no overhead compared to a manual
// update (assuming you call it with the proper CM)

template<CombineMode CM, typename ScalarIn, typename ScalarOut,
         typename CoeffType = typename ekat::ScalarTraits<ScalarIn>::scalar_type>
KOKKOS_FORCEINLINE_FUNCTION
void combine (const ScalarIn& newVal, ScalarOut& result,
              const CoeffType alpha = CoeffType(1),
              const CoeffType beta  = CoeffType(0))
{
  switch (CM) {
    case CombineMode::Replace:
      result = newVal;
      break;
    case CombineMode::Rescale:
      result *= beta;
      break;
    case CombineMode::ScaleReplace:
      result = alpha*newVal;
      break;
    case CombineMode::Update:
      result *= beta;
      result += newVal;
      break;
    case CombineMode::ScaleUpdate:
      result *= beta;
      result += alpha*newVal;
      break;
    case CombineMode::ScaleAdd:
      result += alpha*newVal;
      break;
    case CombineMode::Add:
      result += newVal;
      break;
    case CombineMode::Multiply:
      result *= newVal;
      break;
    case CombineMode::Divide:
      result /= newVal;
      break;
  }
}
/* Special version of combine that takes a mask into account */
template<CombineMode CM, typename ScalarIn, typename ScalarOut,
         typename CoeffType = typename ekat::ScalarTraits<ScalarIn>::scalar_type>
KOKKOS_FORCEINLINE_FUNCTION
void combine_and_fill (const ScalarIn& newVal, ScalarOut& result, const ScalarOut fill_val,
              const CoeffType alpha = CoeffType(1),
              const CoeffType beta  = CoeffType(0))
{
  switch (CM) {
    case CombineMode::Replace:
      combine<CM>(newVal,result,alpha,beta);
      break;
    case CombineMode::Rescale:
      if (result != fill_val) {
        combine<CM>(newVal,result,alpha,beta);
      }
      break;
    case CombineMode::ScaleReplace:
      if (newVal == fill_val) {
        result = fill_val;
      } else {
        combine<CM>(newVal,result,alpha,beta);
      }
      break;
    case CombineMode::Update:
    case CombineMode::ScaleUpdate:
    case CombineMode::ScaleAdd:
    case CombineMode::Add:
    case CombineMode::Multiply:
    case CombineMode::Divide:
      if (result == fill_val || newVal == fill_val) {
        result = fill_val;
      } else {
        combine<CM>(newVal,result,alpha,beta);
      }
      break;
  }
}

} // namespace scream

#endif // SCREAM_COMBINE_OPS_HPP
