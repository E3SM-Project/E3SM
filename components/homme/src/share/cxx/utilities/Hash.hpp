/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_HASH_HPP
#define HOMMEXX_HASH_HPP

#include <cstdint>

#include "Types.hpp"

/* Utilities to calculate a hash for a given model state. Hash values can be
   compared between runs to find instances of non-BFBness.
 */

namespace Homme {

typedef std::uint64_t HashType;

// Each hash function accumulates v into accum using a hash of its bits.

KOKKOS_INLINE_FUNCTION void hash (const HashType v, HashType& accum) {
  constexpr auto first_bit = 1ULL << 63;
  accum += ~first_bit & v; // no overflow
  accum ^=  first_bit & v; // handle most significant bit  
}

KOKKOS_INLINE_FUNCTION void hash (const double v_, HashType& accum) {
  HashType v;
  std::memcpy(&v, &v_, sizeof(HashType));
  hash(v, accum);
}

// Final n argument for Scalar views is in terms of Real, e.g., NUM_PHYSICAL_LEV
// and not NUM_LEV.
void hash(const int tl,                           // time level
          const ExecViewManaged<Scalar******>& v, // full view
          int n5,                                 // max index of final slot
          HashType& accum);                       // accumulate into input value
void hash(const int tl, const ExecViewManaged<Scalar*****>& v, int n4, HashType& accum);
void hash(const ExecViewManaged<Scalar*****>& v, int n4, HashType& accum);
void hash(const int tl, const ExecViewManaged<Real****>& v, HashType& accum);

// No time level slot.
void hash(const ExecViewManaged<Scalar*****>& v, HashType& accum);

// For Kokkos::parallel_reduce.
template <typename ExecSpace = Kokkos::HostSpace>
struct HashReducer {
  typedef HashReducer reducer;
  typedef HashType value_type;
  typedef Kokkos::View<value_type*, ExecSpace, Kokkos::MemoryUnmanaged> result_view_type;

  KOKKOS_INLINE_FUNCTION HashReducer (value_type& value_) : value(value_) {}
  KOKKOS_INLINE_FUNCTION void join (value_type& dest, const value_type& src) const { hash(src, dest); }
  KOKKOS_INLINE_FUNCTION void init (value_type& val) const { val = 0; }
  KOKKOS_INLINE_FUNCTION value_type& reference () const { return value; }
  KOKKOS_INLINE_FUNCTION bool references_scalar () const { return true; }
  KOKKOS_INLINE_FUNCTION result_view_type view () const { return result_view_type(&value, 1); }

private:
  value_type& value;
};

} // Homme

#endif // HASH_HPP
