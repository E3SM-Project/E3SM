#include "physics_test_data.hpp"

#include <ekat_assert.hpp>

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

void PhysicsTestData::read(std::ifstream& ifile)
{
  EKAT_REQUIRE_MSG (ifile.good(), "Cannot read from input file. Did you forget to open it?\n");
  m_reals.read(ifile);
  m_ints.read(ifile);
  m_bools.read(ifile);
}

void PhysicsTestData::write(std::ofstream& ofile) const
{
  EKAT_REQUIRE_MSG (ofile.good(), "Cannot write to input file. Did you forget to open it?\n");

  m_reals.write(ofile);
  m_ints.write(ofile);
  m_bools.write(ofile);
}


} // namespace scream
