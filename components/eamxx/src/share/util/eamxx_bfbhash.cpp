#include "share/util/eamxx_bfbhash.hpp"

namespace scream {
namespace bfbhash {

static void reduce_hash (void* invec, void* inoutvec, int* len, MPI_Datatype* /* datatype */) {
  const int n = *len;
  const auto* s = reinterpret_cast<const HashType*>(invec);
  auto* d = reinterpret_cast<HashType*>(inoutvec);
  for (int i = 0; i < n; ++i) hash(s[i], d[i]);
}

int all_reduce_HashType (MPI_Comm comm, const HashType* sendbuf, HashType* rcvbuf,
                         int count) {
  static_assert(sizeof(long long int) == sizeof(HashType),
                "HashType must have size sizeof(long long int).");
  MPI_Op op;
  MPI_Op_create(reduce_hash, true, &op);
  const auto stat = MPI_Allreduce(sendbuf, rcvbuf, count, MPI_LONG_LONG_INT, op, comm);
  MPI_Op_free(&op);
  return stat;
}

} // namespace bfbhash
} // namespace scream
