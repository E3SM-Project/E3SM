#ifndef P3_FUNCTIONS_EVAPORATE_SUBLIMATE_PRECIP_IMPL.HPP
#define P3_FUNCTIONS_EVAPORATE_SUBLIMATE_PRECIP_IMPL.HPP

#include "p3_functions.hpp"

namespace scream {
  namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::evaporate_sublimate_precip(
