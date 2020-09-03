#include "physics_test_data.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"

#include <random>

namespace scream {

PhysicsTestData::PhysicsTestData(Int shcol_, const std::vector<Real**>& ptrs_c, const std::vector<Int**>& idx_c) :
  PhysicsTestData(shcol_, 0, 0, {}, {}, ptrs_c, idx_c)
{}

PhysicsTestData::PhysicsTestData(Int shcol_, Int nlev_, const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_c, const std::vector<Int**>& idx_c) :
  PhysicsTestData(shcol_, nlev_, 0, ptrs, {}, ptrs_c, idx_c) {}

PhysicsTestData::PhysicsTestData(Int shcol_, Int nlev_, Int nlevi_,
                                 const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_i,
                                 const std::vector<Real**>& ptrs_c, const std::vector<Int**>& idx_c) :
  shcol(shcol_),
  nlev(nlev_),
  nlevi(nlevi_),
  m_total(shcol_ * nlev_),
  m_totali(shcol_ * nlevi_),
  m_ptrs(ptrs),
  m_ptrs_i(ptrs_i),
  m_ptrs_c(ptrs_c),
  m_indices_c(idx_c),
  m_data(m_ptrs.size() * m_total + m_ptrs_i.size() * m_totali + m_ptrs_c.size() * shcol, 0),
  m_idx_data(idx_c.size() * shcol, 0)
{
  init_ptrs();
}

PhysicsTestData::PhysicsTestData(const PhysicsTestData &rhs,
                                 const std::vector<Real**>& real_ptrs, const std::vector<Int**>& int_ptrs) :
  shcol(rhs.shcol),
  nlev(rhs.nlev),
  nlevi(rhs.nlevi),
  m_total(rhs.m_total),
  m_totali(rhs.m_totali),
  m_ptrs(rhs.m_ptrs.size()),
  m_ptrs_i(rhs.m_ptrs_i.size()),
  m_ptrs_c(rhs.m_ptrs_c.size()),
  m_indices_c(rhs.m_indices_c.size()),
  m_data(rhs.m_data),
  m_idx_data(rhs.m_idx_data)
{
  const size_t expected_real_ptrs = m_ptrs.size() + m_ptrs_i.size() + m_ptrs_c.size();
  EKAT_REQUIRE_MSG(
    real_ptrs.size() == expected_real_ptrs,
    "PhysicsTestData: not enough Real* members given to copy constructor: " << real_ptrs.size() << " != " << expected_real_ptrs);
  EKAT_REQUIRE_MSG(
    int_ptrs.size() == m_indices_c.size(),
    "PhysicsTestData: not enough Int* members given to copy constructor: " << int_ptrs.size() << " != " << m_indices_c.size());

  size_t real_offset = 0;
  for (auto& item : {&m_ptrs, &m_ptrs_i, &m_ptrs_c}) {
    std::copy(real_ptrs.begin() + real_offset, real_ptrs.begin() + real_offset + item->size(), item->begin());
    real_offset += item->size();
  }

  std::copy(int_ptrs.begin(), int_ptrs.end(), m_indices_c.begin());

  init_ptrs();
}

void PhysicsTestData::init_ptrs()
{
  Int offset       = 0;
  Real *data_begin = m_data.data();

  for (size_t i = 0; i < m_ptrs.size(); ++i) {
    *(m_ptrs[i]) = data_begin + offset;
    offset += m_total;
  }

  for (size_t i = 0; i < m_ptrs_i.size(); ++i) {
    *(m_ptrs_i[i]) = data_begin + offset;
    offset += m_totali;
  }

  for (size_t i = 0; i < m_ptrs_c.size(); ++i) {
    *(m_ptrs_c[i]) = data_begin + offset;
    offset += shcol;
  }

  for (size_t i = 0; i < m_indices_c.size(); ++i) {
    *(m_indices_c[i]) = m_idx_data.data() + shcol*i;
  }
}

void PhysicsTestData::randomize(const std::vector<std::pair<void*, std::pair<Real, Real> > >& ranges)
{
  using range_type = std::remove_const<std::remove_reference<decltype(ranges)>::type >::type;
  using PT = std::pair<decltype(m_ptrs)*,Int>;

  std::default_random_engine generator;

  range_type ranges_copy(ranges);

  for (auto& item : { PT{&m_ptrs, m_total} , PT{&m_ptrs_i, m_totali}, PT{&m_ptrs_c, shcol} }) {
    auto& ptrs = *item.first;
    const Int num_per = item.second;
    for (Real** ptr : ptrs) {
      std::pair<Real, Real> random_range = {0.0, 1.0};
      for (auto itr = ranges_copy.begin(); itr != ranges_copy.end(); ++itr) {
        Real* range_ptr = reinterpret_cast<Real*>(itr->first);
        if (*ptr == range_ptr) {
          random_range = itr->second;
          ranges_copy.erase(itr);
          break;
        }
      }
      std::uniform_real_distribution<Real> data_dist(random_range.first, random_range.second);
      for (Int i = 0; i < num_per; ++i) {
        (*ptr)[i] = data_dist(generator);
      }
    }
  }

  for (Int** ptr : m_indices_c) {
    std::pair<Int, Int> random_range = {0, 1};
    for (auto itr = ranges_copy.begin(); itr != ranges_copy.end(); ++itr) {
      Int* range_ptr = reinterpret_cast<Int*>(itr->first);
      if (*ptr == range_ptr) {
        const Real bottom_range = itr->second.first;
        const Real top_range    = itr->second.second;
        EKAT_REQUIRE_MSG(std::ceil(bottom_range) == bottom_range, "Use of non-round float for integer random range:" << bottom_range);
        EKAT_REQUIRE_MSG(std::ceil(top_range) == top_range, "Use of non-round float for integer random range:" << top_range);
        random_range = std::make_pair(std::lround(bottom_range), std::lround(top_range));
        ranges_copy.erase(itr);
        break;
      }
    }
    std::uniform_int_distribution<Int> data_dist(random_range.first, random_range.second);
    for (Int i = 0; i < shcol; ++i) {
      (*ptr)[i] = data_dist(generator);
    }
  }

  EKAT_REQUIRE_MSG(ranges_copy.empty(), "Some randomize specializations had no matches");
}

} // namespace scream
