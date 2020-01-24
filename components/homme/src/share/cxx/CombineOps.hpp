#ifndef HOMMEXX_COMBINE_OPS_HPP
#define HOMMEXX_COMBINE_OPS_HPP

#include "HommexxEnums.hpp"

namespace Homme {

template<CombineMode CM>
KOKKOS_INLINE_FUNCTION
static constexpr bool needsAlpha () {
  return CM==CombineMode::Scale || CM==CombineMode::ScaleAdd || CM==CombineMode::ScaleUpdate;
}

template<CombineMode CM>
KOKKOS_INLINE_FUNCTION
static constexpr bool needsBeta () {
  return CM==CombineMode::Update || CM==CombineMode::ScaleUpdate;
}

// Small helper function to combine a new value with an old one.
// The template argument help reducing the number of operations
// performed (the if is resolved at compile time). In the most
// complete form, the function performs
//    result = beta*result + alpha*newVal
// This routine should have no overhead compared to a manual
// update (assuming you call it with the proper CM)
template<CombineMode CM, typename ScalarIn, typename ScalarOut, typename CoeffType>
KOKKOS_FORCEINLINE_FUNCTION
void combine (const ScalarIn& newVal, ScalarOut& result,
              const CoeffType alpha, const CoeffType beta){
  // Sanity check
  assert ((needsAlpha<CM>() || alpha==CoeffType(1)) &&
          "Error! Input alpha would be discarded by the requested combine mode.\n");
  assert ((needsBeta<CM>() || beta==CoeffType(0)) &&
          "Error! Input beta would be discarded by the requested combine mode.\n");

  switch (CM) {
    case CombineMode::Replace:
      result = newVal;
      break;
    case CombineMode::Scale:
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

} // namespace Homme

#endif // HOMMEXX_COMBINE_OPS_HPP
