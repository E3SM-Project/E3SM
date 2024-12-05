#include "share/atm_process/atmosphere_process.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/scream_array_utils.hpp"
#include "share/util/scream_bfbhash.hpp"
#include "ekat/ekat_assert.hpp"

#include <cstdint>

namespace scream {
namespace {

using ExeSpace = KokkosTypes<DefaultDevice>::ExeSpace;
using bfbhash::HashType;

void hash (const Field::view_dev_t<const Real*>& v,
           const FieldLayout& lo, HashType& accum_out) {
  HashType accum = 0;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<ExeSpace>(0, lo.size()),
    KOKKOS_LAMBDA(const int idx, HashType& accum) {
      bfbhash::hash(v(idx), accum);
    }, bfbhash::HashReducer<>(accum));
  Kokkos::fence();
  bfbhash::hash(accum, accum_out);  
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
      bfbhash::hash(v(i,j), accum);
    }, bfbhash::HashReducer<>(accum));
  Kokkos::fence();
  bfbhash::hash(accum, accum_out);
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
      bfbhash::hash(v(i,j,k), accum);
    }, bfbhash::HashReducer<>(accum));
  Kokkos::fence();
  bfbhash::hash(accum, accum_out);
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
      bfbhash::hash(v(i,j,k,m), accum);
    }, bfbhash::HashReducer<>(accum));
  Kokkos::fence();
  bfbhash::hash(accum, accum_out);
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
      bfbhash::hash(v(i,j,k,m,n), accum);
    }, bfbhash::HashReducer<>(accum));
  Kokkos::fence();
  bfbhash::hash(accum, accum_out);
}

void hash (const Field& f, HashType& accum) {
  const auto& hd = f.get_header();
  const auto& id = hd.get_identifier();
  if (id.data_type() != DataType::DoubleType) return;
  const auto& lo = id.get_layout();
  const auto rank = lo.rank();
  switch (rank) {
  case 1: hash(f.get_view<const Real*    >(), lo, accum); break;
  case 2: hash(f.get_view<const Real**   >(), lo, accum); break;
  case 3: hash(f.get_view<const Real***  >(), lo, accum); break;
  case 4: hash(f.get_view<const Real**** >(), lo, accum); break;
  case 5: hash(f.get_view<const Real*****>(), lo, accum); break;
  default: break;
  }  
}

void hash (const std::list<Field>& fs, HashType& accum) {
  for (const auto& f : fs)
    hash(f, accum);
}

void hash (const std::list<FieldGroup>& fgs, HashType& accum) {
  for (const auto& g : fgs)
    for (const auto& e : g.m_fields)
      hash(*e.second, accum);
}

} // namespace anon

// (mem, nmem) describe an arbitrary device array. If non-0, the array will be
// hashed and reported as a fourth line.
void AtmosphereProcess
::print_global_state_hash (const std::string& label, const bool in, const bool out,
                           const bool internal, const Real* mem, const int nmem) const {
  static constexpr int nslot = 4;
  HashType laccum[nslot] = {0};
  hash(m_fields_in, laccum[0]);
  hash(m_groups_in, laccum[0]);
  hash(m_fields_out, laccum[1]);
  hash(m_groups_out, laccum[1]);
  hash(m_internal_fields, laccum[2]);
  const bool hash_array = mem != nullptr;
  if (hash_array) {
    HashType accum = 0;
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExeSpace>(0, nmem),
      KOKKOS_LAMBDA(const int i, HashType& accum) { bfbhash::hash(mem[i], accum); },
      bfbhash::HashReducer<>(accum));
    Kokkos::fence();
    laccum[3] = accum;
  }
  HashType gaccum[nslot];
  const int nr = hash_array ? nslot : nslot-1;
  bfbhash::all_reduce_HashType(m_comm.mpi_comm(), laccum, gaccum, nr);
  const bool show[] = {in, out, internal, hash_array};
  if (m_comm.am_i_root())
    for (int i = 0; i < nslot; ++i)
      if (show[i])
        fprintf(stderr, "exxhash> %4d-%9.5f %1d %16" PRIx64 " (%s)\n",
                timestamp().get_year(), timestamp().frac_of_year_in_days(),
                i, gaccum[i], label.c_str());
}

void AtmosphereProcess::print_fast_global_state_hash (const std::string& label) const {
  HashType laccum = 0;
  hash(m_fields_in, laccum);
  HashType gaccum;
  bfbhash::all_reduce_HashType(m_comm.mpi_comm(), &laccum, &gaccum, 1);
  if (m_comm.am_i_root())
    fprintf(stderr, "bfbhash> %14d %16" PRIx64 " (%s)\n",
            timestamp().get_num_steps(), gaccum, label.c_str());
}

} // namespace scream
