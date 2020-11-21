#include "physics_test_data.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"

#include <set>
#include <random>

namespace scream {

PhysicsTestDataGeneric::PhysicsTestDataGeneric(
  const std::vector<std::vector<Int> >& dims,
  const std::vector<std::vector<Real**> > reals,
  const std::vector<std::vector<Int**> > ints,
  const std::vector<std::vector<bool**> > bools) :
  m_dims(dims),
  m_reals(reals),
  m_ints(ints),
  m_totals(dims.size(), 0)
{
  Int total_reals = 0, total_ints = 0;
  EKAT_REQUIRE_MSG(dims.size() == (reals.size() + ints.size()), "Length of member lists did not match length of dimensions");

  for (size_t i = 0; i < dims.size(); ++i) {
    const auto& current_dims = dims[i];
    Int total = 1;
    for (const auto& current_dim : current_dims) {
      total *= current_dim;
    }
    m_totals[i] = total;
    Int* current_total = i < m_reals.size() ? &total_reals : &total_ints;
    const Int num_members = i < m_reals.size() ? m_reals[i].size() : m_ints[i-m_reals.size()].size();
    (*current_total) += (total * num_members);
  }

  m_real_data.resize(total_reals, 0);
  m_int_data.resize(total_ints, 0);

  init_ptrs();
}

std::pair<size_t, size_t> PhysicsTestDataGeneric::get_index(void* member) const
{
  for (size_t i = 0; i < m_reals.size(); ++i) {
    for (size_t j = 0; j < m_reals[i].size(); ++j) {
      if (*(m_reals[i][j]) == reinterpret_cast<Real*>(member)) {
        return std::make_pair(i, j);
      }
    }
  }

  for (size_t i = 0; i < m_ints.size(); ++i) {
    for (size_t j = 0; j < m_ints[i].size(); ++j) {
      if (*(m_ints[i][j]) == reinterpret_cast<Int*>(member)) {
        return std::make_pair(i + m_reals.size(), j);
      }
    }
  }
  EKAT_REQUIRE_MSG(false, "No match for member");
}

Int PhysicsTestDataGeneric::total(void* member) const
{
  return m_totals[get_index(member).first];
}

Int PhysicsTestDataGeneric::dim(void* member, size_t dim_idx) const
{
  const Int idx = get_index(member).first;
  EKAT_REQUIRE_MSG(dim_idx < m_dims[idx].size(), "dim_idx out of bounds");
  return m_dims[idx][dim_idx];
}

void PhysicsTestDataGeneric::init_ptrs()
{
  Int real_offset = 0, int_offset = 0;
  Real* real_data_begin = m_real_data.data();
  Int*  int_data_begin  = m_int_data.data();

  for (size_t i = 0; i < m_dims.size(); ++i) {
    const Int total = m_totals[i];
    if (i < m_reals.size()) {
      for (auto& real_member : m_reals[i]) {
        *real_member = real_data_begin + real_offset;
        real_offset += total;
      }
    }
    else {
      for (auto& int_member : m_ints[i - m_reals.size()]) {
        *int_member = int_data_begin + int_offset;
        int_offset += total;
      }
    }
  }
}

PhysicsTestDataGeneric& PhysicsTestDataGeneric::assignment_impl(const PhysicsTestDataGeneric& rhs)
{
  EKAT_REQUIRE_MSG(m_reals.size() == rhs.m_reals.size() && m_ints.size() == rhs.m_ints.size(),
                   "Assignment between incompatible PhysicsTestData");

  m_dims      = rhs.m_dims;
  m_real_data = rhs.m_real_data;
  m_int_data  = rhs.m_int_data;
  m_totals    = rhs.m_totals;

  init_ptrs();

  return *this;
}

void PhysicsTestDataGeneric::randomize(const std::vector<std::pair<void*, std::pair<Real, Real> > >& ranges)
{
  std::default_random_engine generator;
  std::set<std::pair<size_t, size_t> > matched_indices;

  for (const auto& p : ranges) {
    const auto& range = p.second;
    const auto idx_p = get_index(p.first);
    const size_t idx = idx_p.first;
    const size_t sub_idx = idx_p.second;
    if (idx < m_reals.size()) {
      std::uniform_real_distribution<Real> data_dist(range.first, range.second);
      for (Int i = 0; i < m_totals[idx]; ++i) {
        (*(m_reals[idx][sub_idx]))[i] = data_dist(generator);
      }
    }
    else {
      const Real bottom_range = range.first;
      const Real top_range    = range.second;
      EKAT_REQUIRE_MSG(std::ceil(bottom_range) == bottom_range, "Use of non-round float for integer random range:" << bottom_range);
      EKAT_REQUIRE_MSG(std::ceil(top_range) == top_range, "Use of non-round float for integer random range:" << top_range);
      std::uniform_int_distribution<Int> data_dist(std::lround(bottom_range), std::lround(top_range));
      for (Int i = 0; i < m_totals[idx]; ++i) {
        (*(m_ints[idx - m_reals.size()][sub_idx]))[i] = data_dist(generator);
      }
    }
    bool was_inserted = matched_indices.insert(idx_p).second;
    EKAT_REQUIRE_MSG(was_inserted, "Should not have had duplicate");
  }

  std::uniform_real_distribution<Real> default_real_dist(0.0, 1.0);
  for (size_t i = 0; i < m_reals.size(); ++i) {
    for (size_t j = 0; j < m_reals[i].size(); ++j) {
      auto pair_search = std::make_pair(i, j);
      if (matched_indices.count(pair_search) == 0) {
        for (Int k = 0; k < m_totals[i]; ++k) {
          (*(m_reals[i][j]))[k] = default_real_dist(generator);
        }
      }
    }
  }

  std::uniform_int_distribution<Int> default_int_dist(0, 1);
  for (size_t i = 0; i < m_ints.size(); ++i) {
    for (size_t j = 0; j < m_ints[i].size(); ++j) {
      auto pair_search = std::make_pair(i+m_reals.size(), j);
      if (matched_indices.count(pair_search) == 0) {
        for (Int k = 0; k < m_totals[i+m_reals.size()]; ++k) {
          (*(m_ints[i][j]))[k] = default_int_dist(generator);
        }
      }
    }
  }
}

PhysicsTestData::PhysicsTestData(Int dim1_, const std::vector<Real**>& ptrs_c, const std::vector<Int**>& idx_c) :
  PhysicsTestData(dim1_, 0, 0, {}, {}, ptrs_c, idx_c)
{}

PhysicsTestData::PhysicsTestData(Int dim1_, Int dim2_, const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_c, const std::vector<Int**>& idx_c) :
  PhysicsTestData(dim1_, dim2_, 0, ptrs, {}, ptrs_c, idx_c) {}

PhysicsTestData::PhysicsTestData(Int dim1_, Int dim2_, Int dim3_,
                                 const std::vector<Real**>& ptrs, const std::vector<Real**>& ptrs_i,
                                 const std::vector<Real**>& ptrs_c, const std::vector<Int**>& idx_c,
                                 const std::vector<Real**>& ptrs_3d) :
  PhysicsTestDataGeneric({ {dim1_, dim2_}, {dim1_, dim3_}, {dim1_},  {dim1_, dim2_, dim3_}, {dim1_} },
                         { {ptrs},         {ptrs_i},       {ptrs_c}, {ptrs_3d} },         { {idx_c} } ),
  dim1(dim1_),
  dim2(dim2_),
  dim3(dim3_),
  m_total(dim1 * dim2),
  m_totali(dim1 * dim3),
  m_total3d(dim1 * dim2 * dim3)
{}

PhysicsTestData& PhysicsTestData::assignment_impl(const PhysicsTestData& rhs)
{
  PhysicsTestDataGeneric::assignment_impl(rhs);

  dim1       = rhs.dim1;
  dim2       = rhs.dim2;
  dim3       = rhs.dim3;
  m_total    = rhs.m_total;
  m_totali   = rhs.m_totali;
  m_total3d  = rhs.m_total3d;

  return *this;
}

} // namespace scream
