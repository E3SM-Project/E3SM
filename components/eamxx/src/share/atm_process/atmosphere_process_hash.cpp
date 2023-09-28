#include "share/atm_process/atmosphere_process.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/scream_array_utils.hpp"
#include "ekat/ekat_assert.hpp"

#include <cstdint>

namespace scream {

namespace {
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
  const auto stat = MPI_Allreduce(sendbuf, rcvbuf, count, MPI_LONG_LONG_INT, op, comm);
  MPI_Op_free(&op);
  return stat;
}

using ExeSpace = KokkosTypes<DefaultDevice>::ExeSpace;

void hash (const Field::view_dev_t<const Real*>& v,
           const FieldLayout& lo, HashType& accum_out) {
  HashType accum = 0;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<ExeSpace>(0, lo.size()),
    KOKKOS_LAMBDA(const int idx, HashType& accum) {
      hash(v(idx), accum);
    }, HashReducer<>(accum));
  Kokkos::fence();
  hash(accum, accum_out);  
}

void hash (const Field::view_dev_t<const Real**>& v,
           const FieldLayout& lo, HashType& accum_out) {
  HashType accum = 0;
  const auto& dims = lo.extents();
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<ExeSpace>(0, lo.size()),
    KOKKOS_LAMBDA(const int idx, HashType& accum) {
      int i, j;
      unflatten_idx(idx, dims, i, j);
      hash(v(i,j), accum);
    }, HashReducer<>(accum));
  Kokkos::fence();
  hash(accum, accum_out);
}

void hash (const Field::view_dev_t<const Real***>& v,
           const FieldLayout& lo, HashType& accum_out) {
  HashType accum = 0;
  const auto& dims = lo.extents();
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<ExeSpace>(0, lo.size()),
    KOKKOS_LAMBDA(const int idx, HashType& accum) {
      int i, j, k;
      unflatten_idx(idx, dims, i, j, k);
      hash(v(i,j,k), accum);
    }, HashReducer<>(accum));
  Kokkos::fence();
  hash(accum, accum_out);
}

void hash (const Field::view_dev_t<const Real****>& v,
           const FieldLayout& lo, HashType& accum_out) {
  HashType accum = 0;
  const auto& dims = lo.extents();
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<ExeSpace>(0, lo.size()),
    KOKKOS_LAMBDA(const int idx, HashType& accum) {
      int i, j, k, m;
      unflatten_idx(idx, dims, i, j, k, m);
      hash(v(i,j,k,m), accum);
    }, HashReducer<>(accum));
  Kokkos::fence();
  hash(accum, accum_out);
}

void hash (const Field::view_dev_t<const Real*****>& v,
           const FieldLayout& lo, HashType& accum_out) {
  HashType accum = 0;
  const auto& dims = lo.extents();
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<ExeSpace>(0, lo.size()),
    KOKKOS_LAMBDA(const int idx, HashType& accum) {
      int i, j, k, m, n;
      unflatten_idx(idx, dims, i, j, k, m, n);
      hash(v(i,j,k,m,n), accum);
    }, HashReducer<>(accum));
  Kokkos::fence();
  hash(accum, accum_out);
}

void hash (const std::list<FieldGroup>& fgs, HashType& accum_out) {
  
}

void hash (const std::list<Field>& fs, HashType& accum) {
  for (const auto& f : fs) {
    const auto& hd = f.get_header();
    const auto& id = hd.get_identifier();
    if (id.data_type() != DataType::DoubleType) continue;
    const auto& lo = id.get_layout();
    const auto rank = lo.rank();
    switch (rank) {
    case 1: hash(f.get_view<const Real*    >(), lo, accum); break;
    case 2: hash(f.get_view<const Real**   >(), lo, accum); break;
    case 3: hash(f.get_view<const Real***  >(), lo, accum); break;
    case 4: hash(f.get_view<const Real**** >(), lo, accum); break;
    case 5: hash(f.get_view<const Real*****>(), lo, accum); break;
    default: continue;
    }
  }
}
} // namespace anon

void AtmosphereProcess::print_global_state_hash (const std::string& label) const {
  static constexpr int nslot = 3;
  HashType laccum[nslot] = {0};
  hash(m_fields_in, laccum[0]);
  hash(m_groups_in, laccum[0]);
  hash(m_fields_out, laccum[1]);
  hash(m_groups_out, laccum[1]);
  hash(m_internal_fields, laccum[2]);
  HashType gaccum[nslot];
  all_reduce_HashType(m_comm.mpi_comm(), laccum, gaccum, nslot);
  if (m_comm.am_i_root())
    for (int i = 0; i < nslot; ++i)
      fprintf(stderr, "exxhash> %4d-%9.5f %1d %16lx (%s)\n",
              timestamp().get_year(), timestamp().frac_of_year_in_days(),
              i, gaccum[i], label.c_str());
}

void AtmosphereProcess::print_fast_global_state_hash (const std::string& label) const {
  HashType laccum = 0;
  hash(m_fields_in, laccum);
  HashType gaccum;
  all_reduce_HashType(m_comm.mpi_comm(), &laccum, &gaccum, 1);
  if (m_comm.am_i_root())
    fprintf(stderr, "bfbhash> %14d %16lx (%s)\n",
            timestamp().get_num_steps(), gaccum, label.c_str());
}

} // namespace scream
