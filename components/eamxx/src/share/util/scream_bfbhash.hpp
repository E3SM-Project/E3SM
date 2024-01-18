#ifndef SCREAM_BFBHASH_HPP
#define SCREAM_BFBHASH_HPP

#include <cstdint>

#include <ekat/kokkos/ekat_kokkos_types.hpp>
#include <ekat/mpi/ekat_comm.hpp>

namespace scream {
namespace bfbhash {

typedef std::uint64_t HashType;

KOKKOS_INLINE_FUNCTION void hash (const HashType v, HashType& accum) {
  constexpr auto first_bit = 1ULL << 63;
  accum += ~first_bit & v; // no overflow
  accum ^=  first_bit & v; // handle most significant bit  
}

KOKKOS_INLINE_FUNCTION void hash (const double v_, HashType& accum) {
  static_assert(sizeof(double) == sizeof(HashType),
                "HashType must have size sizeof(double).");
  HashType v;
  std::memcpy(&v, &v_, sizeof(HashType));
  hash(v, accum);
}

KOKKOS_INLINE_FUNCTION void hash (const float v, HashType& accum) {
  hash(double(v), accum);
}

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

int all_reduce_HashType(MPI_Comm comm, const HashType* sendbuf, HashType* rcvbuf,
                        int count);

} // namespace bfbhash
} // namespace scream

#endif
