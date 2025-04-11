#include "share/atm_process/atmosphere_process.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_array_utils.hpp"
#include "share/util/eamxx_bfbhash.hpp"
#include "ekat/ekat_assert.hpp"

#include <cstdint>
#include <iomanip>

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
    for (const auto& e : g.m_individual_fields)
      hash(*e.second, accum);
}

} // namespace anon

void AtmosphereProcess
::print_global_state_hash (const std::string& label, const TimeStamp& t,
                           const bool in, const bool out, const bool internal,
                           const Real* mem, const int nmem) const
{
  const bool compute[4] = {in, out, internal, mem!=nullptr};

  std::vector<std::string> hash_names;
  std::vector<HashType> laccum;

  // When calling printf later, how much space does the hash name take (we'll update later)
  int slen = 0;
  // Compute local hashes
  if (m_internal_diagnostics_level==1) {
    // Lump fields together (but keep in/out/internal separated)
    if (compute[0]) {
      laccum.emplace_back();
      hash(m_fields_in, laccum.back());
      hash(m_groups_in, laccum.back());
      hash_names.push_back("inputs");
    }
    if (compute[1]) {
      laccum.emplace_back();
      hash(m_fields_out, laccum.back());
      hash(m_groups_out, laccum.back());
      hash_names.push_back("outputs");
    }
    if (compute[2]) {
      laccum.emplace_back();
      hash(m_internal_fields, laccum.back());
      hash_names.push_back("internals");
    }

    slen = 10;
  } else if (m_internal_diagnostics_level==2) {
    // Hash fields individually. Notice that, if a field is requested individually
    // as well as part of a group, it will be hashed twice (independently)
    auto layout = [](const Field& f) -> std::string {
      const auto& fl = f.get_header().get_identifier().get_layout();
      return " (" + ekat::join(fl.names(),",") + ")";
    };
    if (compute[0]) {
      for (const auto& f : m_fields_in) {
        laccum.emplace_back();
        hash_names.push_back(f.name()+layout(f));
        hash(f,laccum.back());
      }
      for (const auto& g : m_groups_in) {
        for (const auto& [fn,f] : g.m_individual_fields) {
          laccum.emplace_back();
          hash_names.push_back(fn+layout(*f));
          hash(*f,laccum.back());
        }
      }
    }
    if (compute[1]) {
      for (const auto& f : m_fields_out) {
        laccum.emplace_back();
        hash_names.push_back(f.name()+layout(f));
        hash(f,laccum.back());
      }
      for (const auto& g : m_groups_out) {
        for (const auto& [fn,f] : g.m_individual_fields) {
          laccum.emplace_back();
          hash_names.push_back(fn+layout(*f));
          hash(*f,laccum.back());
        }
      }
    }
    if (compute[2]) {
      for (const auto& f : m_internal_fields) {
        laccum.emplace_back();
        hash_names.push_back(f.name()+layout(f));
        hash(f,laccum.back());
      }
    }
  }

  if (compute[3]) {
    laccum.emplace_back();

    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExeSpace>(0, nmem),
      KOKKOS_LAMBDA(const int i, HashType& accum) { bfbhash::hash(mem[i], accum); },
      bfbhash::HashReducer<>(laccum.back()));
    Kokkos::fence();

    hash_names.push_back("mem-buffer");
  }

  // Compute global hashes
  int naccum = laccum.size();
  std::vector<HashType> gaccum(naccum);
  bfbhash::all_reduce_HashType(m_comm.mpi_comm(), laccum.data(), gaccum.data(), naccum);

  // Print
  if (m_comm.am_i_root()) {
    // Find out how long the field names are, to avoid truncating them
    for (const auto& n : hash_names)
      slen = std::max(slen,static_cast<int>(n.size()+1));

    const auto& date = t.get_date();
    const auto& tod  = t.sec_of_day();

    std::stringstream ss;

    ss << "eamxx hash> date="
       << std::setw(4) << std::setfill('0') << date[0] << "-" // Year
       << std::setw(2) << std::setfill('0') << date[1] << "-" // Month
       << std::setw(2) << std::setfill('0') << date[2] << "-" // Day
       << std::setw(5) << std::setfill('0') << tod << " "     // Time of day
       << "(" << label << "), naccum=" << naccum;             // Label and number of accum
    log (ss.str());

    for (int i = 0; i < naccum; ++i) {
      ss.str(""); // Clear content
      ss.clear(); // Clear error flags
      ss << std::setw(slen) << std::setfill(' ') << hash_names[i] << ": "
         << std::hex << std::setfill('0') << std::setw(16) << gaccum[i];
      log (ss.str());
    }
    m_atm_logger->flush();
  }
}

void AtmosphereProcess::
print_fast_global_state_hash (const std::string& label, const TimeStamp& t) const
{
  HashType laccum = 0;
  hash(m_fields_in, laccum);
  HashType gaccum;
  bfbhash::all_reduce_HashType(m_comm.mpi_comm(), &laccum, &gaccum, 1);
  if (m_comm.am_i_root())
    fprintf(stderr, "bfbhash> %14d %16" PRIx64 " (%s)\n",
            t.get_num_steps(), gaccum, label.c_str());
}

} // namespace scream
