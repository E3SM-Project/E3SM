/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Hash.hpp"

namespace Homme {

void hash (const int tl, const ExecViewManaged<Scalar******>& v, int n5,
           HashType& accum_out) {
  HashType accum;
  Kokkos::parallel_reduce(
    MDRangePolicy<ExecSpace, 6>(
      {0, tl, 0, 0, 0, 0},
      {v.extent_int(0), tl+1, v.extent_int(2), v.extent_int(3), v.extent_int(4), n5}),
    KOKKOS_LAMBDA(int i0, int i1, int i2, int i3, int i4, int i5, HashType& accum) {
      const auto* vcol = &v(i0,i1,i2,i3,i4,0)[0];
      Homme::hash(vcol[i5], accum);
    }, HashReducer<>(accum));
  hash(accum, accum_out);
}

void hash (const int tl, const ExecViewManaged<Scalar*****>& v, int n4,
           HashType& accum_out) {
  HashType accum;
  Kokkos::parallel_reduce(
    MDRangePolicy<ExecSpace, 5>(
      {0, tl, 0, 0, 0},
      {v.extent_int(0), tl+1, v.extent_int(2), v.extent_int(3), n4}),
    KOKKOS_LAMBDA(int i0, int i1, int i2, int i3, int i4, HashType& accum) {
      const auto* vcol = &v(i0,i1,i2,i3,0)[0];
      Homme::hash(vcol[i4], accum);
    }, HashReducer<>(accum));
  hash(accum, accum_out);
}

void hash (const ExecViewManaged<Scalar*****>& v, int n4,
           HashType& accum_out) {
  HashType accum;
  Kokkos::parallel_reduce(
    MDRangePolicy<ExecSpace, 5>(
      {0, 0, 0, 0, 0},
      {v.extent_int(0), v.extent_int(1), v.extent_int(2), v.extent_int(3), n4}),
    KOKKOS_LAMBDA(int i0, int i1, int i2, int i3, int i4, HashType& accum) {
      const auto* vcol = &v(i0,i1,i2,i3,0)[0];
      Homme::hash(vcol[i4], accum);
    }, HashReducer<>(accum));
  hash(accum, accum_out);
}

void hash (const int tl, const ExecViewManaged<Real****>& v,
           HashType& accum_out) {
  HashType accum;
  Kokkos::parallel_reduce(
    MDRangePolicy<ExecSpace, 4>(
      {0, tl, 0, 0},
      {v.extent_int(0), tl+1, v.extent_int(2), v.extent_int(3)}),
    KOKKOS_LAMBDA(int i0, int i1, int i2, int i3, HashType& accum) {
      Homme::hash(v(i0,i1,i2,i3), accum);
    }, HashReducer<>(accum));
  hash(accum, accum_out);
}

} // Homme
