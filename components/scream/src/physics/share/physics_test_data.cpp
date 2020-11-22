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
  m_reals(dims.begin(), dims.begin() + reals.size(), reals),
  m_ints(dims.begin() + reals.size(), dims.begin() + reals.size() + ints.size(), ints),
  m_bools(dims.begin() + reals.size() + ints.size(), dims.end(), bools),
{}

std::pair<size_t, size_t> PhysicsTestDataGeneric::get_index(void* member) const
{
  std::pair<size_t, size_t> search = m_reals.get_index(reinterpret_cast<Real*>(member));
  if (search.first != -1) {
    return search;
  }
  search = m_ints.get_index(reinterpret_cast<Int*>(member));
  if (search.first != -1) {
    return search;
  }
  search = m_bools.get_index(reinterpret_cast<bool*>(member));
  if (search.first != -1) {
    return search;
  }
  EKAT_REQUIRE_MSG(false, "No match for member");
}

PhysicsTestDataGeneric& PhysicsTestDataGeneric::assignment_impl(const PhysicsTestDataGeneric& rhs)
{
  m_reals.assignment_impl(rhs.m_reals);
  m_ints.assignment_impl(rhs.m_ints);
  m_bools.assignment_impl(rhs.m_bools);

  return *this;
}

void PhysicsTestDataGeneric::randomize(const std::vector<std::pair<void*, std::pair<Real, Real> > >& ranges)
{
  std::default_random_engine generator;

  // generate with default vals
  m_reals.randomize(generator, std::uniform_real_distribution<Real>(0.0, 1.0));
  m_ints.randomize(generator, std::uniform_int_distribution<Int>(0, 1));
  m_bools.randomize(generator, std::uniform_int_distribution<Int>(0, 1));

  // override defauls if user requested something specific
  for (const auto& p : ranges) {
    const auto& range = p.second;
    const Real bottom_range = range.first;
    const Real top_range    = range.second;
    void* member = p.first;

    const auto real_search = get_index(reinterpret_cast<Real*>(member));
    if (real_search.first != -1) {
      m_reals.randomize(generator, std::uniform_real_distribution(range.first, range.second), real_search);
    }
    else {
      const auto int_search = get_index(reinterpret_cast<Int*>(member));
      if (int_search.first != -1) {
        EKAT_REQUIRE_MSG(std::ceil(bottom_range) == bottom_range, "Use of non-round float for integer random range:" << bottom_range);
        EKAT_REQUIRE_MSG(std::ceil(top_range) == top_range, "Use of non-round float for integer random range:" << top_range);
        std::uniform_int_distribution<Int> data_dist(std::lround(bottom_range), std::lround(top_range));

        m_ints.randomize(generator, data_dist, int_search);
      }
      else {
        const auto bool_search = get_index(reinterpret_cast<bool*>(member));
        EKAT_REQUIRE_MSG(bool_search.first != -1, "Failed to find member for randomization");
        EKAT_REQUIRE_MSG(bottom_range == 0.0, "Use of non-zero for bottom of bool random range:" << bottom_range);
        EKAT_REQUIRE_MSG(top_range == 0.0, "Use of non-1 for top of bool random range:" << top_range);
        std::uniform_int_distribution<Int> data_dist(std::lround(bottom_range), std::lround(top_range));
        m_bools.randomize(generator, data_dist, bool_search);
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
