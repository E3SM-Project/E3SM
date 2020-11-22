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

  SHOCGridData(Int dim1_, Int dim2_, Int dim3_) :
    PhysicsTestData(dim1_, dim2_, dim3_,
      {&zt_grid, &dz_zt, &pdel, &rho_zt},  // list of (dim1 x dim2) members
      {&zi_grid, &dz_zi}) {}               // list of (dim1 x dim3) members

  PTD_STD_DEF(SHOCGridData, 3, 0); // 3 => number of dimensions (1-3), 0 => number of scalars (0-10)
  // If you have scalars, you'll have to add names after the count
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

#define PTD_ZEROES(a) PTD_ZEROES##a

#define PTD_DATA_COPY_CTOR(name, num_args) \
  name(const name& rhs) : name(PTD_ZEROES(num_args)) { *this = rhs; }

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

#define PTD_DIM_RENAME1(ndim1)                                              Int ndim1() const { return dim1; }
#define PTD_DIM_RENAME2(ndim1, ndim2)        PTD_DIM_RENAME1(ndim1);        Int ndim2() const { return dim2; }
#define PTD_DIM_RENAME3(ndim1, ndim2, ndim3) PTD_DIM_RENAME2(ndim1, ndim2); Int ndim3() const { return dim3; }

#define PTD_DIM_RENAME(num_args, ...)           \
  PTD_DIM_RENAME##num_args(__VA_ARGS__)

#define PTD_ADD_1_0 1
#define PTD_ADD_1_1 2
#define PTD_ADD_1_2 3
#define PTD_ADD_1_3 4
#define PTD_ADD_1_4 5
#define PTD_ADD_1_5 6
#define PTD_ADD_1_6 7
#define PTD_ADD_1_7 8
#define PTD_ADD_1_8 9
#define PTD_ADD_1_9 10

#define PTD_ADD_2_0 2
#define PTD_ADD_2_1 3
#define PTD_ADD_2_2 4
#define PTD_ADD_2_3 5
#define PTD_ADD_2_4 6
#define PTD_ADD_2_5 7
#define PTD_ADD_2_6 8
#define PTD_ADD_2_7 9
#define PTD_ADD_2_8 10
#define PTD_ADD_2_9 11

#define PTD_ADD_3_0 3
#define PTD_ADD_3_1 4
#define PTD_ADD_3_2 5
#define PTD_ADD_3_3 6
#define PTD_ADD_3_4 7
#define PTD_ADD_3_5 8
#define PTD_ADD_3_6 9
#define PTD_ADD_3_7 10
#define PTD_ADD_3_8 11
#define PTD_ADD_3_9 12

#define PTD_STD_DEF(name, dim, num_scalars, ...)        \
  PTD_DATA_COPY_CTOR(name, PTD_ADD_##dim##_##num_scalars);  \
  PTD_ASSIGN_OP(name, num_scalars, __VA_ARGS__)

// For PhysicsTestDataGeneric, the dimensions are also struct scalars
#define PTDG_STD_DEF(name, num_scalars, ...) \
  PTD_DATA_COPY_CTOR(name, num_scalars);     \
  PTD_ASSIGN_OP(name, num_scalars, __VA_ARGS__)

namespace scream {

// Fully Generic Data struct for multi-dimensions reals and ints
class PhysicsTestDataGeneric
{
  template <typename T>
  struct PTDImpl
  {
    template <typename Iterator>
    PTDImpl( Iterator dims_begin, Iterator dims_end,
             const std::vector<std::vector<T**> >& members_list) :
      m_dims_list(dims_begin, dims_end),
      m_members_list(members_list),
      m_totals(m_dims_list.size(), 0)
    {
      EKAT_REQUIRE_MSG(m_dims_list.size() == m_members_list.size(),
                       "Length of member lists did not match length of dimensions");

      // Compute totals
      Int total_total = 0;
      for (size_t i = 0; i < m_dims_list.size(); ++i) {
        const auto& dims    = m_dims_list[i];
        const auto& members = m_members_list[i];

        Int total = 1;
        for (const auto& dim : dims) {
          total *= dim;
        }
        m_totals[i] = total;
        const Int num_members = members.size();
        total_total += (total * num_members);
      }

      m_data.resize(total_total, T());

      init_ptrs();
    }

    std::pair<size_t, size_t> get_index(T* member) const
    {
      for (size_t i = 0; i < m_members_list.size(); ++i) {
        for (size_t j = 0; j < m_members_list[i].size(); ++j) {
          if (*(m_members_list[i][j]) == member) {
            return std::make_pair(i, j);
          }
        }
      }
      return std::make_pair(-1, -1);
    }

    Int total(const size_t& index) const
    {
      EKAT_ASSERT(index < m_totals.size());
      return m_totals[index];
    }

    template <typename Gen, typename Dist>
    void randomize(Gen& generator, Dist& dist)
    {
      for (auto& item : m_data) {
        item = static_cast<T>(dist(generator));
      }
    }

    template <typename Gen, typename Dist>
    void randomize(Gen& generator, Dist& dist, const std::pair<size_t, size_t>& index)
    {
      const Int total = total(index.first);
      T* member = *(m_members_list[index.first][index.second]);
      for (Int i = 0; i < total; ++i) {
        member[i] = static_cast<T>(dist(generator));
      }
    }

    void init_ptrs()
    {
      Int offset = 0;
      for (size_t i = 0; i < m_members_list.size(); ++i) {
        const Int total     = m_totals[i];
        const auto& members = m_members_list[i];

        for (auto& member : members) {
          *member = m_data.data() + offset;
          offset += total;
        }
      }
    }

    void assignment_impl(const PTDImpl<T>& rhs)
    {
      EKAT_REQUIRE_MSG(m_members.size() == rhs.m_members.size(),
                       "Assignment between incompatible PhysicsTestData");

      m_dims_list = rhs.m_dim_list;
      m_data      = rhs.m_data;
      m_totals    = rhs.m_totals;

      init_ptrs();
    }

    template <ekat::TransposeDirection::Enum D>
    void transpose()
    {
      Int offset = 0;
      std::vector<T> new_data(m_data);

      for (size_t i = 0; i < m_dims_list.size(); ++i) {
        const auto& dims      = m_dims_list[i];
        const auto& members   = m_members_list[i];
        const Int num_members = members.size();
        const Int total       = m_totals[i];

        if (dims.size() > 1) { // no need to transpose 1d data
          if (dims.size() == 2) {
            for (auto member : members) {
              ekat::transpose<D>(*member, new_data.data() + offset, dims[0], dims[1]);
              offset += total;
            }
          }
          else if (dims.size() == 3) {
            for (auto member : members) {
              ekat::transpose<D>(*member, new_data.data() + offset, dims[0], dims[1], dims[2]);
              offset += total;
            }
          }
          else {
            EKAT_REQUIRE_MSG(false, "Data dimension > 3 not currently supported");
          }
        }
        else {
          offset += (total * num_members);
        }
      }

      m_data = new_data;
    }

  private:
    std::vector<std::vector<Int> m_dims_list;    // list of dims, one per unique set of dims
    std::vector<std::vector<T**> m_members_list; // list of member pointers, same outer index space as m_dims_list
    std::vector<T>               m_data;         // the member data in a flat vector
    std::vector<Int>             m_totals;       // total sizes of each set of data, same index space as m_dims_list
  };

 public:

  // dims -> the dimensions of real data should come before dimensions of int data
  //         and the dims of int data should come before bool data
  PhysicsTestDataGeneric(
    const std::vector<std::vector<Int> >& dims, // vector of dimensions, each set of dimensions is a vector of Int
    const std::vector<std::vector<Real**> > reals, // vector of pointers to real* members
    const std::vector<std::vector<Int**> > ints = {},  // vector of pointers to int* members
    const std::vector<std::vector<bool**> > bools = {}); // vector of pointers to bool* members

  Int total(const Real* member) const { return m_reals.total(get_index(member).first); }
  Int total(const Int* member) const  { return m_ints.total(get_index(member).first); }
  Int total(const bool* member) const { return m_bools.total(get_index(member).first); }

  // Delete this to block subclasses getting the default impls, which would be incorrect
  PhysicsTestDataGeneric(const PhysicsTestDataGeneric &rhs) = delete;
  PhysicsTestDataGeneric &operator=(const PhysicsTestDataGeneric &rhs) = delete;

  // Initialize with random values. The default range is 0..1
  // To use non-default ranges, you'll need to provide a pair of pairs, mapping
  // the member to a range.
  // Example, to use a -1 to 1 range for wthl member:
  // d.randomize({ {d.wthl, {-1, 1}} });
  void randomize(const std::vector<std::pair<void*, std::pair<Real, Real> > >& ranges = {});

  // Since we are also preparing index data, this function is doing more than transposing. It's shifting the
  // format of all data from one language to another
  template <ekat::TransposeDirection::Enum D>
  void transpose()
  {
    m_reals.transpose<D>();
    m_ints.transpose<D>();
    m_bools.transpose<D>();

    // Shift the indices. We might not be able to make the assumption that int data represented indices
    for (size_t i = 0; i < m_ints.m_data.size(); ++i) {
      m_ints.m_data[i] += (D == ekat::TransposeDirection::c2f ? 1 : -1);
      EKAT_ASSERT_MSG(m_ints.m_data[i] >= 0, "Bad index: " << m_ints.m_data[i]);
    }
  }

 protected:

  PhysicsTestDataGeneric& assignment_impl(const PhysicsTestDataGeneric& rhs);

 private:

  std::pair<size_t, size_t> get_index(const Real* member) const { return m_reals.get_index(member); }
  std::pair<size_t, size_t> get_index(const Int* member)  const { return m_reals.get_index(member); }
  std::pair<size_t, size_t> get_index(const bool* member) const { return m_reals.get_index(member); }

  PTDImpl<Real> m_reals; // manage real data with this member
  PTDImpl<Int>  m_ints;  // manage int data with this member
  PTDImpl<bool> m_bools; // manage bool data with this member
};

// Base class for common test data setup
struct PhysicsTestData : public PhysicsTestDataGeneric
{
  Int dim1, dim2, dim3;

  PhysicsTestData(Int dim1_,
                  const std::vector<Real**>& ptrs_c,     // [dim1]
                  const std::vector<Int**>& idx_c = {}); // [dim1] (optional)

  PhysicsTestData(Int dim1_, Int dim2_,
                  const std::vector<Real**>& ptrs,        // [dim1 x dim2]
                  const std::vector<Real**>& ptrs_c = {}, // [dim1] (optional)
                  const std::vector<Int**>& idx_c = {});  // [dim1] (optional)

  PhysicsTestData(Int dim1_, Int dim2_, Int dim3_,
                  const std::vector<Real**>& ptrs,          // [dim1 x dim2]
                  const std::vector<Real**>& ptrs_i,        // [dim1 x dim3]
                  const std::vector<Real**>& ptrs_c = {},   // [dim1] (optional)
                  const std::vector<Int**>& idx_c = {},     // [dim1] (optional)
                  const std::vector<Real**>& ptrs_3d = {}); // [dim1 x dim2 x dim3] (optional)

  Int total1x2() const { return m_total; }
  Int total1x3() const { return m_totali; }
  Int total3d()  const { return m_total3d; }

 protected:

  PhysicsTestData& assignment_impl(const PhysicsTestData& rhs);

 private:
  // Internals
  Int m_total, m_totali, m_total3d;
  std::vector<Real**> m_ptrs, m_ptrs_i, m_ptrs_c, m_ptrs_3d;
  std::vector<Int**> m_indices_c;
  std::vector<Real> m_data;
  std::vector<Int> m_idx_data;
};

}  // namespace scream

#endif // SCREAM_PHYSICS_TEST_DATA_HPP
