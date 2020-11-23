#include "physics_test_data.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"

#include <set>
#include <random>

namespace scream {

PhysicsTestData::PhysicsTestData(
  const std::vector<std::vector<Int> >& dims,
  const std::vector<std::vector<Real**> >& reals,
  const std::vector<std::vector<Int**> >& ints,
  const std::vector<std::vector<bool**> >& bools) :
  m_reals(dims.begin(), dims.begin() + reals.size(), reals),
  m_ints(dims.begin() + reals.size(), dims.begin() + reals.size() + ints.size(), ints),
  m_bools(dims.begin() + reals.size() + ints.size(), dims.end(), reinterpret_cast<const std::vector<std::vector<char**> >&>(bools))
{
}

PhysicsTestData& PhysicsTestData::assignment_impl(const PhysicsTestData& rhs)
{
  m_reals.assignment_impl(rhs.m_reals);
  m_ints.assignment_impl(rhs.m_ints);
  m_bools.assignment_impl(rhs.m_bools);

  return *this;
}

void PhysicsTestData::randomize(const std::vector<std::pair<void*, std::pair<Real, Real> > >& ranges)
{
  std::default_random_engine generator;
  std::uniform_real_distribution<Real> default_real_dist(0.0, 1.0);
  std::uniform_int_distribution<Int> default_int_dist(0, 1);
  std::uniform_int_distribution<Int> default_bool_dist(0, 1);

  // generate with default vals
  m_reals.randomize(generator, default_real_dist);
  m_ints.randomize(generator, default_int_dist);
  m_bools.randomize(generator, default_bool_dist);

  // override defauls if user requested something specific
  for (const auto& p : ranges) {
    const auto& range = p.second;
    const Real bottom_range = range.first;
    const Real top_range    = range.second;
    void* member = p.first;

    const auto real_search = get_index(reinterpret_cast<Real*>(member));
    if (real_search.first != std::string::npos) {
      std::uniform_real_distribution<Real> real_dist(range.first, range.second);
      m_reals.randomize(generator, real_dist, real_search);
    }
    else {
      const auto int_search = get_index(reinterpret_cast<Int*>(member));
      if (int_search.first != std::string::npos) {
        EKAT_REQUIRE_MSG(std::ceil(bottom_range) == bottom_range, "Use of non-round float for integer random range:" << bottom_range);
        EKAT_REQUIRE_MSG(std::ceil(top_range) == top_range, "Use of non-round float for integer random range:" << top_range);
        std::uniform_int_distribution<Int> data_dist(std::lround(bottom_range), std::lround(top_range));

        m_ints.randomize(generator, data_dist, int_search);
      }
      else {
        const auto bool_search = get_index(reinterpret_cast<bool*>(member));
        EKAT_REQUIRE_MSG(bool_search.first != std::string::npos, "Failed to find member for randomization");
        EKAT_REQUIRE_MSG(bottom_range == 0.0, "Use of non-zero for bottom of bool random range:" << bottom_range);
        EKAT_REQUIRE_MSG(top_range == 0.0, "Use of non-1 for top of bool random range:" << top_range);
        std::uniform_int_distribution<Int> data_dist(std::lround(bottom_range), std::lround(top_range));
        m_bools.randomize(generator, data_dist, bool_search);
      }
    }
  }
}

} // namespace scream
