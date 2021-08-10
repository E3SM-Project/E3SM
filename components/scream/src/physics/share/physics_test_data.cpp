#include "physics_test_data.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_pack_kokkos.hpp"

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

} // namespace scream
