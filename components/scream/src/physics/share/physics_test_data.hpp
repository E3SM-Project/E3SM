#ifndef SCREAM_PHYSICS_TEST_DATA_HPP
#define SCREAM_PHYSICS_TEST_DATA_HPP

#include "share/scream_types.hpp"

#include "ekat/util/ekat_math_utils.hpp"
#include "ekat/ekat_assert.hpp"

#include <vector>
#include <utility>

//
// PhysicsTestData is meant to offer the client something they can inherit to provide
// convenient handling of arrays of data in the common *Data structs that are used for
// unit-testing and bridging. This class supports up to 4 classes of data
//   * dim1 of Real
//   * (dim1 x dim2) of Real
//   * (dim1 x dim3) of Real
//   * dim1 of Int
//
// Subclasses of PhysicsTestData should look like the following:
/*
struct SHOCGridData : public PhysicsTestData {
  // Inputs
  Real *zt_grid, *zi_grid, *pdel;

  // In/out
  Real *dz_zt, *dz_zi, *rho_zt;

  SHOCGridData(Int shcol_, Int nlev_, Int nlevi_) :
    PhysicsTestData(shcol_, nlev_, nlevi_,
      {&zt_grid, &dz_zt, &pdel, &rho_zt},  // list of (shcol x nlev) members
      {&zi_grid, &dz_zi}) {}               // list of (shcol x nlevi) members

  SHOCGridData(const SHOCGridData &rhs) : PhysicsTestData(rhs,
    {&zt_grid, &dz_zt, &pdel, &rho_zt, &zi_grid, &dz_zi}) {} // list of all PhysicsTestData managed Real members
    // order of these members should match the order they appeared in the constructor

  SHOCGridData &operator=(const SHOCGridData &rhs) { PhysicsTestData::operator=(rhs); return *this; }
};
*/

namespace scream {

// Base class for common test data setup
struct PhysicsTestData
{
  Int shcol, nlev, nlevi;

  PhysicsTestData(Int shcol_,
                  const std::vector<Real**>& ptrs_c,     // [shcol]
                  const std::vector<Int**>& idx_c = {}); // [shcol] (optional)

  PhysicsTestData(Int shcol_, Int nlev_,
                  const std::vector<Real**>& ptrs,        // [shcol x nlev]
                  const std::vector<Real**>& ptrs_c = {}, // [shcol] (optional)
                  const std::vector<Int**>& idx_c = {});  // [shcol] (optional)

  PhysicsTestData(Int shcol_, Int nlev_, Int nlevi_,
                  const std::vector<Real**>& ptrs,        // [schol x nlev]
                  const std::vector<Real**>& ptrs_i,      // [schol x nlevi]
                  const std::vector<Real**>& ptrs_c = {}, // [schol] (optional)
                  const std::vector<Int**>& idx_c = {});  // [schol] (optional)

  PhysicsTestData(const PhysicsTestData &rhs) = delete;

  PhysicsTestData(const PhysicsTestData &rhs,
                  const std::vector<Real**>& real_ptrs, // ALL Real* members, listed in same order as constructor but without breaking them into multiple vectors
                  const std::vector<Int**>& int_ptrs = {}); // ALL Int* members (optional)

  PhysicsTestData &operator=(const PhysicsTestData &rhs) = delete;

  // Since we are also preparing index data, this function is doing more than transposing. It's shifting the
  // format of all data from one language to another
  template <ekat::util::TransposeDirection::Enum D>
  void transpose() {
    std::vector<Real> data(m_data.size());

    // Transpose on the zt grid
    for (size_t i = 0; i < m_ptrs.size(); ++i) {
      ekat::util::transpose<D>(*(m_ptrs[i]), data.data() + (m_total*i) , shcol, nlev);
    }

    // Transpose on the zi grid
    for (size_t i = 0; i < m_ptrs_i.size(); ++i) {
      ekat::util::transpose<D>(*(m_ptrs_i[i]), data.data() + (m_ptrs.size()*m_total) + (m_totali*i), shcol, nlevi);
    }

    // Copy the column only grid
    const Int c_start_offset = m_ptrs.size()*m_total + m_ptrs_i.size()*m_totali;
    std::copy(m_data.begin() + c_start_offset, m_data.end(), data.begin() + c_start_offset);

    m_data = data;

    // Shift the indices. We might not be able to make the assumption that int data represented indices
    for (size_t i = 0; i < m_idx_data.size(); ++i) {
      m_idx_data[i] += (D == ekat::util::TransposeDirection::c2f ? 1 : -1);
      EKAT_ASSERT_MSG(m_idx_data[i] >= 0, "Bad index: " << m_idx_data[i]);
    }
  }

  // Initialize with random values. The default range is 0..1
  // To use non-default ranges, you'll need to provide a pair of pairs, mapping
  // the member to a range.
  // Example, to use a -1 to 1 range for wthl member:
  // d.randomize({ {d.wthl, {-1, 1}} });
  void randomize(const std::vector<std::pair<void*, std::pair<Real, Real> > >& ranges = {});

  Int total() const { return m_total; }
  Int totali() const { return m_totali; }

 private:
  void init_ptrs();

  // Internals
  Int m_total, m_totali;
  std::vector<Real**> m_ptrs, m_ptrs_i, m_ptrs_c;
  std::vector<Int**> m_indices_c;
  std::vector<Real> m_data;
  std::vector<Int> m_idx_data;
};

}  // namespace scream

#endif // SCREAM_PHYSICS_TEST_DATA_HPP
