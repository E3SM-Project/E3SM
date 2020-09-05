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
// Subclasses of PhysicsTestData should look like the following. Note that the copy
// constructor and assignment operators must be defined if you want to be able to copy
// objects of this type. The PTR_DATA_COPY_CTOR and PTD_ASSIGN_OP macros are there to help you do this.
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

  PTD_DATA_COPY_CTOR(SHOCGridData, 3); // 3 => number of arguments constructor expects
  PTD_ASSIGN_OP(SHOCGridData, 0); // 0 => number of scalars
};
*/

// Convenience macros for up to 11 arguments, beyond that, you're on your own :)

#define PTD_ZEROES0
#define PTD_ZEROES1 0
#define PTD_ZEROES2 PTD_ZEROES1, 0
#define PTD_ZEROES3 PTD_ZEROES2, 0
#define PTD_ZEROES4 PTD_ZEROES3, 0
#define PTD_ZEROES5 PTD_ZEROES4, 0
#define PTD_ZEROES6 PTD_ZEROES5, 0
#define PTD_ZEROES7 PTD_ZEROES6, 0
#define PTD_ZEROES8 PTD_ZEROES7, 0
#define PTD_ZEROES9 PTD_ZEROES8, 0
#define PTD_ZEROES10 PTD_ZEROES9, 0

#define PTD_DATA_COPY_CTOR(name, num_args) \
  name(const name& rhs) : name(PTD_ZEROES##num_args) { *this = rhs; }

#define  PTD_ASS0(                                ) ((void) (0))
#define  PTD_ASS1(a                               )                                           a = rhs.a
#define  PTD_ASS2(a, b                            )  PTD_ASS1(a)                            ; b = rhs.b
#define  PTD_ASS3(a, b, c                         )  PTD_ASS2(a, b)                         ; c = rhs.c
#define  PTD_ASS4(a, b, c, d                      )  PTD_ASS3(a, b, c)                      ; d = rhs.d
#define  PTD_ASS5(a, b, c, d, e                   )  PTD_ASS4(a, b, c, d)                   ; e = rhs.e
#define  PTD_ASS6(a, b, c, d, e, f                )  PTD_ASS5(a, b, c, d, e)                ; f = rhs.f
#define  PTD_ASS7(a, b, c, d, e, f, g             )  PTD_ASS6(a, b, c, d, e, f)             ; g = rhs.g
#define  PTD_ASS8(a, b, c, d, e, f, g, h          )  PTD_ASS7(a, b, c, d, e, f, g)          ; h = rhs.h
#define  PTD_ASS9(a, b, c, d, e, f, g, h, i       )  PTD_ASS8(a, b, c, d, e, f, g, h)       ; i = rhs.i
#define PTD_ASS10(a, b, c, d, e, f, g, h, i, j    )  PTD_ASS9(a, b, c, d, e, f, g, h, i)    ; j = rhs.j
#define PTD_ASS11(a, b, c, d, e, f, g, h, i, j, k ) PTD_ASS10(a, b, c, d, e, f, g, h, i, j) ; k = rhs.k

#define PTD_ASSIGN_OP(name, num_scalars, ...)                                  \
  name& operator=(const name& rhs) { PTD_ASS##num_scalars(__VA_ARGS__); assignment_impl(rhs); return *this; }

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

  // Delete this to block subclasses getting the default impls, which would be incorrect
  PhysicsTestData(const PhysicsTestData &rhs) = delete;
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

 protected:

  PhysicsTestData& assignment_impl(const PhysicsTestData& rhs);

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
