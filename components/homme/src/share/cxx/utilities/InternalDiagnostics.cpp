/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "InternalDiagnostics.hpp"

#include "Context.hpp"
#include "ElementsState.hpp"
#include "ElementsDerivedState.hpp"
#include "Tracers.hpp"
#include "TimeLevel.hpp"
#include "mpi/Comm.hpp"

namespace Homme {

namespace {
void reduce_hash (void* invec, void* inoutvec, int* len, MPI_Datatype* /* datatype */) {
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
  return MPI_Allreduce(sendbuf, rcvbuf, count, MPI_LONG_LONG_INT, op, comm);
  MPI_Op_free(&op);
}
} // anon

void print_global_state_hash (const std::string& label) {
  const auto& c = Context::singleton();
  const auto& es = c.get<ElementsState>();
  const auto& tr = c.get<Tracers>();
  const auto& tl = c.get<TimeLevel>();
  const auto& comm = c.get<Comm>();
  HashType accum[NUM_TIME_LEVELS + Q_NUM_TIME_LEVELS] = {0};
  hash(es.hash(tl.nm1), accum[0]);
  hash(es.hash(tl.n0 ), accum[1]);
  hash(es.hash(tl.np1), accum[2]);
  hash(tr.hash(tl.n0_qdp ), accum[3]);
  hash(tr.hash(tl.np1_qdp), accum[4]);
  HashType gaccum[NUM_TIME_LEVELS + Q_NUM_TIME_LEVELS];
  all_reduce_HashType(comm.mpi_comm(), accum, gaccum, 5);
  if (comm.root()) {
    for (int i = 0; i < NUM_TIME_LEVELS; ++i)
      fprintf(stderr, "hxxhash> %14d %1d %16lx (E %s)\n",
              tl.nstep, i, gaccum[i], label.c_str());
    for (int i = 0; i < Q_NUM_TIME_LEVELS; ++i)
      fprintf(stderr, "hxxhash> %14d %1d %16lx (T %s)\n",
              tl.nstep, i, gaccum[NUM_TIME_LEVELS+i], label.c_str());
  }
}

} // Homme
