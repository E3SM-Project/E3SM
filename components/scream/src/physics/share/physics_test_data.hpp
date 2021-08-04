#ifndef SCREAM_PHYSICS_TEST_DATA_HPP
#define SCREAM_PHYSICS_TEST_DATA_HPP

#include "share/scream_types.hpp"

#include "ekat/util/ekat_math_utils.hpp"
#include "ekat/ekat_assert.hpp"

#include <random>
#include <vector>
#include <utility>

/*
PhysicsTestData is meant to offer the client something they can inherit to provide
convenient handling of arrays of data in the common *Data structs that are used for
unit-testing and bridging. This class supports storing reals, ints, and bools of
any multidimensionality.

Subclasses of PhysicsTestData should look like the following. Note that the copy
constructor and assignment operators must be defined if you want to be able to copy
objects of this type. The PTD_STD_DEF macro is there to help you do this.

struct SHOCGridData : public PhysicsTestData {
  // Inputs
  Int dim1, dim2, dim3;
  Real *zt_grid, *zi_grid, *pdel;

  // In/out
  Real *dz_zt, *dz_zi, *rho_zt;

  SHOCGridData(Int dim1_, Int dim2_, Int dim3_) :
    PhysicsTestData({ {dim1_, dim2_}, {dim1_, dim3_} },
      {&zt_grid, &dz_zt, &pdel, &rho_zt},  // list of (dim1 x dim2) members
      {&zi_grid, &dz_zi}),                 // list of (dim1 x dim3) members
      dim1(dim1_), dim2(dim2_), dim3(dim3_)// initialize your own scalars
  {}

  PTD_STD_DEF(SHOCGridData, 3, dim1, dim2, dim3); // 3 => number of scalars followed by their names
};
*/

// Convenience macros for up to 11 arguments, beyond that, you're on your own :)

#define PTD_ONES0
#define PTD_ONES1 1
#define PTD_ONES2 PTD_ONES1, 1
#define PTD_ONES3 PTD_ONES2, 1
#define PTD_ONES4 PTD_ONES3, 1
#define PTD_ONES5 PTD_ONES4, 1
#define PTD_ONES6 PTD_ONES5, 1
#define PTD_ONES7 PTD_ONES6, 1
#define PTD_ONES8 PTD_ONES7, 1
#define PTD_ONES9 PTD_ONES8, 1
#define PTD_ONES10 PTD_ONES9, 1
#define PTD_ONES11 PTD_ONES10, 1

#define PTD_ONES(a) PTD_ONES##a

#define PTD_DATA_COPY_CTOR(name, num_args) \
  name(const name& rhs) : name(PTD_ONES(num_args)) { *this = rhs; }

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

#define PTD_STD_DEF(name, num_scalars, ...) \
  PTD_DATA_COPY_CTOR(name, num_scalars);     \
  PTD_ASSIGN_OP(name, num_scalars, __VA_ARGS__)

namespace scream {

// Fully Generic Data struct for multi-dimensions reals and ints
class PhysicsTestData
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

    std::pair<size_t, size_t> get_index(const T* member) const
    {
      for (size_t i = 0; i < m_members_list.size(); ++i) {
        for (size_t j = 0; j < m_members_list[i].size(); ++j) {
          if (*(m_members_list[i][j]) == member) {
            return std::make_pair(i, j);
          }
        }
      }
      return std::make_pair(std::string::npos, std::string::npos);
    }

    Int get_dim(const size_t& dims_index, const size_t& dim_index) const
    {
      EKAT_ASSERT(dims_index < m_dims_list.size());
      EKAT_ASSERT(dim_index  < m_dims_list[dims_index].size());
      return m_dims_list[dims_index][dim_index];
    }

    Int get_total(const size_t& index) const
    {
      EKAT_ASSERT(index < m_totals.size());
      return m_totals[index];
    }

    template <typename Gen, typename Dist>
    void randomize(Gen& generator, Dist& dist)
    {
      for (auto&& item : m_data) {
        item = static_cast<T>(dist(generator));
      }
    }

    template <typename Gen, typename Dist>
    void randomize(Gen& generator, Dist& dist, const std::pair<size_t, size_t>& index)
    {
      const Int total = get_total(index.first);
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

        for (auto member : members) {
          *member = m_data.data() + offset;
          offset += total;
        }
      }
    }

    void assignment_impl(const PTDImpl<T>& rhs)
    {
      EKAT_REQUIRE_MSG(m_members_list.size() == rhs.m_members_list.size(),
                       "Assignment between incompatible PhysicsTestData");

      m_dims_list = rhs.m_dims_list;
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

    std::vector<std::vector<Int> > m_dims_list;    // list of dims, one per unique set of dims
    std::vector<std::vector<T**> > m_members_list; // list of member pointers, same outer index space as m_dims_list
    std::vector<T>                 m_data;         // the member data in a flat vector
    std::vector<Int>               m_totals;       // total sizes of each set of data, same index space as m_dims_list
  };

 public:

  // dims -> the dimensions of real data should come before dimensions of int data
  //         and the dims of int data should come before bool data
  PhysicsTestData(
    const std::vector<std::vector<Int> >& dims, // vector of dimensions, each set of dimensions is a vector of Int
    const std::vector<std::vector<Real**> >& reals, // vector of pointers to real* members
    const std::vector<std::vector<Int**> >& ints = {},  // vector of pointers to int* members
    const std::vector<std::vector<bool**> >& bools = {}); // vector of pointers to bool* members

  Int total(const Real* member) const { return m_reals.get_total(get_index(member).first); }
  Int total(const Int* member) const  { return m_ints.get_total(get_index(member).first); }
  Int total(const bool* member) const { return m_bools.get_total(get_index(member).first); }

  Int dim(const Real* member, const size_t& dim_idx) const { return m_reals.get_dim(get_index(member).first, dim_idx); }
  Int dim(const Int * member, const size_t& dim_idx) const { return m_ints.get_dim(get_index(member).first, dim_idx); }
  Int dim(const bool* member, const size_t& dim_idx) const { return m_bools.get_dim(get_index(member).first, dim_idx); }

  // Delete this to block subclasses getting the default impls, which would be incorrect
  PhysicsTestData(const PhysicsTestData &rhs) = delete;
  PhysicsTestData &operator=(const PhysicsTestData &rhs) = delete;

  // Initialize with random values. The default range is 0..1
  // To use non-default ranges, you'll need to provide a pair of pairs, mapping
  // the member to a range.
  // Example, to use a -1 to 1 range for wthl member:
  // d.randomize({ {d.wthl, {-1, 1}} });
  template <typename Engine>
  void randomize(Engine& engine, const std::vector<std::pair<void*, std::pair<Real, Real> > >& ranges = {})
  {
    std::uniform_real_distribution<Real> default_real_dist(0.0, 1.0);
    std::uniform_int_distribution<Int> default_int_dist(0, 1);
    std::uniform_int_distribution<Int> default_bool_dist(0, 1);

    // generate with default vals
    m_reals.randomize(engine, default_real_dist);
    m_ints.randomize(engine, default_int_dist);
    m_bools.randomize(engine, default_bool_dist);

    // override defaults if user requested something specific
    for (const auto& p : ranges) {
      const auto& range = p.second;
      const Real bottom_range = range.first;
      const Real top_range    = range.second;
      EKAT_REQUIRE_MSG(bottom_range <= top_range, "Expect bottom of range <= top of range");
      void* member = p.first;

      const auto real_search = get_index(reinterpret_cast<Real*>(member));
      if (real_search.first != std::string::npos) {
        std::uniform_real_distribution<Real> real_dist(range.first, range.second);
        m_reals.randomize(engine, real_dist, real_search);
      }
      else {
        const auto int_search = get_index(reinterpret_cast<Int*>(member));
        if (int_search.first != std::string::npos) {
          EKAT_REQUIRE_MSG(std::ceil(bottom_range) == bottom_range, "Use of non-round float for integer random range:" << bottom_range);
          EKAT_REQUIRE_MSG(std::ceil(top_range) == top_range, "Use of non-round float for integer random range:" << top_range);
          std::uniform_int_distribution<Int> data_dist(std::lround(bottom_range), std::lround(top_range));

          m_ints.randomize(engine, data_dist, int_search);
        }
        else {
          const auto bool_search = get_index(reinterpret_cast<bool*>(member));
          EKAT_REQUIRE_MSG(bool_search.first != std::string::npos, "Failed to find member for randomization");
          EKAT_REQUIRE_MSG(bottom_range == 0.0 || bottom_range == 1.0, "Use 0 or 1 for bool ranges, not:" << bottom_range);
          EKAT_REQUIRE_MSG(top_range == 0.0 || top_range == 1.0, "Use 0 or 1 for bool ranges, not:" << top_range);
          std::uniform_int_distribution<Int> data_dist(std::lround(bottom_range), std::lround(top_range));
          m_bools.randomize(engine, data_dist, bool_search);
        }
      }
    }
  }

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

  PhysicsTestData& assignment_impl(const PhysicsTestData& rhs);

 private:

  std::pair<size_t, size_t> get_index(const Real* member) const { return m_reals.get_index(member); }
  std::pair<size_t, size_t> get_index(const Int* member)  const { return m_ints.get_index(member); }
  std::pair<size_t, size_t> get_index(const bool* member) const { return m_bools.get_index(reinterpret_cast<const char*>(member)); }

  PTDImpl<Real> m_reals; // manage real data with this member
  PTDImpl<Int>  m_ints;  // manage int data with this member
  PTDImpl<char> m_bools; // manage bool data with this member, use chars internally to dodge vector<bool> specialization
};

}  // namespace scream

#endif // SCREAM_PHYSICS_TEST_DATA_HPP
