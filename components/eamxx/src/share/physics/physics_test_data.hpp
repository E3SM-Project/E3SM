#ifndef SCREAM_PHYSICS_TEST_DATA_HPP
#define SCREAM_PHYSICS_TEST_DATA_HPP

#include "physics_test_read_write_helpers.hpp"

#include "share/core/eamxx_types.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include <ekat_math_utils.hpp>
#include <ekat_assert.hpp>
#include <ekat_test_utils.hpp>

#include <random>
#include <vector>
#include <utility>
#include <fstream>

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

// Convenience macros for up to 20 arguments, beyond that, you're on your own :)

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
#define PTD_ONES12 PTD_ONES11, 1
#define PTD_ONES13 PTD_ONES12, 1
#define PTD_ONES14 PTD_ONES13, 1
#define PTD_ONES15 PTD_ONES14, 1
#define PTD_ONES16 PTD_ONES15, 1
#define PTD_ONES17 PTD_ONES16, 1
#define PTD_ONES18 PTD_ONES17, 1
#define PTD_ONES19 PTD_ONES18, 1
#define PTD_ONES20 PTD_ONES19, 1

#define PTD_ONES(a) PTD_ONES##a

// We want these to be deep copies. The impl sets all scalars to 1
// and then calls assignment operator. The assignment operator will
// set all the scalars to the correct value. We cannot use default-contstructed
// PTDs, because the list of member pointers for array members must be registered.
#define PTD_DATA_COPY_CTOR(name, num_args) \
  name(const name& rhs) : name(PTD_ONES(num_args)) { *this = rhs; }

#define PTD_DATA_COPY_CTOR_INIT(name, num_args)                         \
  name(const name& rhs) : name(PTD_ONES(num_args), rhs.init) { *this = rhs; }

#define  PTD_ASS0()           ((void) (0))
#define  PTD_ASS1(first)      first = rhs.first; PTD_ASS0()
#define  PTD_ASS2(first, ...) first = rhs.first; PTD_ASS1(__VA_ARGS__)
#define  PTD_ASS3(first, ...) first = rhs.first; PTD_ASS2(__VA_ARGS__)
#define  PTD_ASS4(first, ...) first = rhs.first; PTD_ASS3(__VA_ARGS__)
#define  PTD_ASS5(first, ...) first = rhs.first; PTD_ASS4(__VA_ARGS__)
#define  PTD_ASS6(first, ...) first = rhs.first; PTD_ASS5(__VA_ARGS__)
#define  PTD_ASS7(first, ...) first = rhs.first; PTD_ASS6(__VA_ARGS__)
#define  PTD_ASS8(first, ...) first = rhs.first; PTD_ASS7(__VA_ARGS__)
#define  PTD_ASS9(first, ...) first = rhs.first; PTD_ASS8(__VA_ARGS__)
#define PTD_ASS10(first, ...) first = rhs.first; PTD_ASS9(__VA_ARGS__)
#define PTD_ASS11(first, ...) first = rhs.first; PTD_ASS10(__VA_ARGS__)
#define PTD_ASS12(first, ...) first = rhs.first; PTD_ASS11(__VA_ARGS__)
#define PTD_ASS13(first, ...) first = rhs.first; PTD_ASS12(__VA_ARGS__)
#define PTD_ASS14(first, ...) first = rhs.first; PTD_ASS13(__VA_ARGS__)
#define PTD_ASS15(first, ...) first = rhs.first; PTD_ASS14(__VA_ARGS__)
#define PTD_ASS16(first, ...) first = rhs.first; PTD_ASS15(__VA_ARGS__)
#define PTD_ASS17(first, ...) first = rhs.first; PTD_ASS16(__VA_ARGS__)
#define PTD_ASS18(first, ...) first = rhs.first; PTD_ASS17(__VA_ARGS__)
#define PTD_ASS19(first, ...) first = rhs.first; PTD_ASS18(__VA_ARGS__)
#define PTD_ASS20(first, ...) first = rhs.first; PTD_ASS19(__VA_ARGS__)

// Assign scalars and call assignment_impl, which will handle the array data
#define PTD_ASSIGN_OP(name, num_scalars, ...)                                  \
  name& operator=(const name& rhs) { PTD_ASS##num_scalars(__VA_ARGS__); assignment_impl(rhs); return *this; }

#define PTD_ASSIGN_OP_INIT(name, num_scalars, ...)                      \
  name& operator=(const name& rhs) { PTD_ASS##num_scalars(__VA_ARGS__); assignment_impl(rhs); init = rhs.init; return *this; }

#define PTD_RW_SCALARS(num_scalars, ...)                                \
  void read_scalars (std::ifstream& ifile) {                            \
    EKAT_REQUIRE_MSG (ifile.good(),                                     \
        "Cannot read from input file. Did you forget to open it?\n");   \
    ::scream::impl::read_scalars(ifile,__VA_ARGS__);                    \
  }                                                                     \
  void write_scalars (std::ofstream& ofile) {                           \
    EKAT_REQUIRE_MSG (ofile.good(),                                     \
        "Cannot write to output file. Did you forget to open it?\n");   \
    ::scream::impl::write_scalars(ofile,__VA_ARGS__);                   \
  }

#define PTD_RW_SCALARS_ONLY(num_scalars, ...)                           \
  void read(std::ifstream& ifile) {                                     \
    EKAT_REQUIRE_MSG (ifile.good(),                                     \
        "Cannot read from input file. Did you forget to open it?\n");   \
    ::scream::impl::read_scalars(ifile,__VA_ARGS__);                    \
  }                                                                     \
  void write(std::ofstream& ofile) {                                    \
    EKAT_REQUIRE_MSG (ofile.good(),                                     \
        "Cannot write to output file. Did you forget to open it?\n");   \
    ::scream::impl::write_scalars(ofile,__VA_ARGS__);                   \
  }

#define PTD_RW() \
  void read(std::ifstream& ifile) {           \
    read_scalars(ifile);                      \
    PhysicsTestData::read(ifile);             \
  }                                           \
  void write(std::ofstream& ofile) {          \
    write_scalars(ofile);                     \
    PhysicsTestData::write(ofile);            \
  }

#define PTD_STD_DEF(name, num_scalars, ...) \
  PTD_DATA_COPY_CTOR(name, num_scalars);     \
  PTD_ASSIGN_OP(name, num_scalars, __VA_ARGS__) \
  PTD_RW() \
  PTD_RW_SCALARS(num_scalars, __VA_ARGS__)

#define PTD_STD_DEF_INIT(name, num_scalars, ...) \
  PTD_DATA_COPY_CTOR_INIT(name, num_scalars);     \
  PTD_ASSIGN_OP_INIT(name, num_scalars, __VA_ARGS__) \
  PTD_RW() \
  PTD_RW_SCALARS(num_scalars, __VA_ARGS__)

namespace scream {

template <ekat::TransposeDirection::Enum D>
void shift_int_scalar(int& scalar)
{
  const int shift = (D == ekat::TransposeDirection::c2f ? 1 : -1);
  scalar += shift;

  // Since f90 allows for negative index ranges (-foo:foo), we may
  // have to remove this check.
  EKAT_ASSERT_MSG(scalar >= 0, "Bad index: " << scalar);
}

// Fully Generic Data struct for multi-dimensions reals and ints
class PhysicsTestData
{

  // PTDImpl handles all arrays of data of a certain type. The data within a
  // a single PTDImpl object may belong to several different arrays.
  template <typename T>
  struct PTDImpl
  {

    // Constructor, takes iterators that contain dimension data for members.
    template <typename Iterator>
    PTDImpl( Iterator dims_begin, Iterator dims_end,
             const std::vector<std::vector<T**> >& members_list_by_common_dims) :
      m_dims_list(),
      m_members_list(),
      m_totals(),
      m_data()
    {
      std::vector<std::vector<Int> > dims_slice(dims_begin, dims_end);
      EKAT_REQUIRE_MSG(dims_slice.size() == members_list_by_common_dims.size(),
                       "Length of member lists did not match length of dimensions slice");

      // Compute totals
      Int total_data = 0;
      for (size_t i = 0; i < dims_slice.size(); ++i) {
        const auto& dims    = dims_slice[i];
        const auto& members = members_list_by_common_dims[i];

        Int total = 1;
        for (const auto& dim : dims) {
          total *= dim;
        }

        const Int num_members = members.size();
        total_data += (total * num_members);

        for (auto member : members) {
          m_dims_list.push_back(dims);
          m_members_list.push_back(member);
          m_totals.push_back(total);
        }
      }

      m_data.resize(total_data, T());

      init_ptrs();
    }

    // Get the location of a member within members_list
    size_t get_index(const T* member) const
    {
      for (size_t i = 0; i < m_members_list.size(); ++i) {
        if (*(m_members_list[i]) == member) {
          return i;
        }
      }
      return std::string::npos;
    }

    Int get_dim(const T* member, const size_t& dim_idx) const
    {
      const auto idx = get_index(member);
      EKAT_ASSERT(idx < m_dims_list.size());
      EKAT_ASSERT(dim_idx  < m_dims_list[idx].size());
      return m_dims_list[idx][dim_idx];
    }

    Int get_total(const T* member) const
    {
      const auto idx = get_index(member);
      EKAT_ASSERT(idx < m_totals.size());
      return m_totals[idx];
    }

    template <typename Gen, typename Dist>
    void randomize(Gen& generator, Dist& dist)
    {
      for (auto&& item : m_data) {
        item = static_cast<T>(dist(generator));
      }
    }

    template <typename Gen, typename Dist>
    void randomize(Gen& generator, Dist& dist, size_t idx)
    {
      T** member = m_members_list[idx];
      const Int total = m_totals[idx];
      for (Int i = 0; i < total; ++i) {
        (*member)[i] = static_cast<T>(dist(generator));
      }
    }

    void init_ptrs()
    {
      Int offset = 0;
      for (size_t i = 0; i < m_totals.size(); ++i) {
        const Int total = m_totals[i];
        T**     member  = m_members_list[i];

        T* data = m_data.data() + offset;
        *member = data; // This sets the pointer member of the PTD child struct!

        offset += total;
      }
    }

    void assignment_impl(const PTDImpl<T>& rhs)
    {
      EKAT_REQUIRE_MSG(m_members_list.size() == rhs.m_members_list.size(),
                       "Assignment between incompatible PhysicsTestData?");

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
        const auto& member    = m_members_list[i];
        const Int total       = m_totals[i];

        if (dims.size() > 1) { // no need to transpose 1d data
          if (dims.size() == 2) {
            ekat::transpose<D>(*member, new_data.data() + offset, dims[0], dims[1]);
          }
          else if (dims.size() == 3) {
            ekat::transpose<D>(*member, new_data.data() + offset, dims[0], dims[1], dims[2]);
          }
          else {
            EKAT_REQUIRE_MSG(false, "Data dimension > 3 not currently supported");
          }
        }
        offset += total;
      }

      m_data = new_data;
    }

    template <ekat::TransposeDirection::Enum D>
    void shift_int(const std::vector<int*>& to_skip = {})
    {
      for (size_t i = 0; i < m_members_list.size(); ++i) {
        auto member = m_members_list[i];
        auto it = std::find(to_skip.begin(), to_skip.end(), *member);
        if (it == to_skip.end()) {
          const Int total = m_totals[i];
          for (Int j = 0; j < total; ++j) {
            shift_int_scalar<D>((*member)[j]);
          }
        }
      }
    }

    void read(std::ifstream& ifile)
    {
      impl::read_scalars(ifile, m_data);
    }

    void write(std::ofstream& ofile) const
    {
      impl::write_scalars(ofile, m_data);
    }

    std::vector<std::vector<Int> > m_dims_list;    // list of dims, one vector per array
    std::vector<T**>               m_members_list; // list of member pointers, one item per array
    std::vector<Int>               m_totals;       // total sizes of each set of data, one item per array
    std::vector<T>                 m_data;         // the member data in a flat vector
  };

 public:

  // dims -> the dimensions of real data should come before dimensions of int data
  //         and the dims of int data should come before bool data
  //
  // We use pointers to pointers because we want to set the pointers to the allocated
  // data.
  PhysicsTestData(
    const std::vector<std::vector<Int> >& dims, // vector of dimensions, each set of dimensions is a vector of Int
    const std::vector<std::vector<Real**> >& reals, // vector of pointers to real* members
    const std::vector<std::vector<Int**> >& ints = {},  // vector of pointers to int* members
    const std::vector<std::vector<bool**> >& bools = {}); // vector of pointers to bool* members

  Int total(const Real* member) const { return m_reals.get_total(member); }
  Int total(const Int* member) const  { return m_ints.get_total(member); }
  Int total(const bool* member) const
  { return m_bools.get_total(reinterpret_cast<const char*>(member)); }

  Int dim(const Real* member, const size_t& dim_idx) const { return m_reals.get_dim(member, dim_idx); }
  Int dim(const Int * member, const size_t& dim_idx) const { return m_ints.get_dim(member, dim_idx); }
  Int dim(const bool* member, const size_t& dim_idx) const
  { return m_bools.get_dim(reinterpret_cast<const char*>(member), dim_idx); }

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

      const auto real_search = m_reals.get_index(reinterpret_cast<Real*>(member));
      if (real_search != std::string::npos) {
        std::uniform_real_distribution<Real> real_dist(range.first, range.second);
        m_reals.randomize(engine, real_dist, real_search);
      }
      else {
        const auto int_search = m_ints.get_index(reinterpret_cast<Int*>(member));
        if (int_search != std::string::npos) {
          EKAT_REQUIRE_MSG(std::ceil(bottom_range) == bottom_range, "Use of non-round float for integer random range:" << bottom_range);
          EKAT_REQUIRE_MSG(std::ceil(top_range) == top_range, "Use of non-round float for integer random range:" << top_range);
          std::uniform_int_distribution<Int> data_dist(std::lround(bottom_range), std::lround(top_range));

          m_ints.randomize(engine, data_dist, int_search);
        }
        else {
          const auto bool_search = m_bools.get_index(reinterpret_cast<char*>(member));
          EKAT_REQUIRE_MSG(bool_search != std::string::npos, "Failed to find member for randomization");
          EKAT_REQUIRE_MSG(bottom_range == 0.0 || bottom_range == 1.0, "Use 0 or 1 for bool ranges, not:" << bottom_range);
          EKAT_REQUIRE_MSG(top_range == 0.0 || top_range == 1.0, "Use 0 or 1 for bool ranges, not:" << top_range);
          std::uniform_int_distribution<Int> data_dist(std::lround(bottom_range), std::lround(top_range));
          m_bools.randomize(engine, data_dist, bool_search);
        }
      }
    }
  }

  // Since we are also preparing index data, this function is doing more than transposing. It's shifting the
  // format of all data from one language to another.
  //
  // There is currently no way for this struct to know which integer scalars need
  // to be shifted, so any subclass that has those will need to define their own
  // transition method which will call this one and then adjust their int scalars
  // that represent indices.
  template <ekat::TransposeDirection::Enum D>
  void transition(const std::vector<Int*>& ints_to_skip = {})
  {
    m_reals.transpose<D>();
    m_ints.transpose<D>();
    m_bools.transpose<D>();

    // Shift the indices. We might not be able to make the assumption that int data represented indices.
    // NOTE! This will not shift scalar integers. It is up the children structs to do that
    m_ints.shift_int<D>(ints_to_skip);
  }

  void read(std::ifstream& ifile);

  void write(std::ofstream& ofile) const;

 protected:

  PhysicsTestData& assignment_impl(const PhysicsTestData& rhs);

 private:

  PTDImpl<Real> m_reals; // manage real data with this member
  PTDImpl<Int>  m_ints;  // manage int data with this member
  PTDImpl<char> m_bools; // manage bool data with this member, use chars internally to dodge vector<bool> specialization
};

enum BASELINE_ACTION {
  NONE,
  COMPARE,
  GENERATE
};

/**
 * In $phys_unit_tests_common.hpp, the UnitWrap struct should have an inner struct "Base"
 * that inherits from the struct below. This will ensure common BFB baseline unit tests
 * are set up in a consistent manner.
 */
struct UnitBase
{

  std::string     m_baseline_path;
  std::string     m_test_name;
  BASELINE_ACTION m_baseline_action;
  std::ifstream   m_ifile;
  std::ofstream   m_ofile;

  UnitBase() :
    m_baseline_path(""),
    m_test_name(Catch::getResultCapture().getCurrentTestName()),
    m_baseline_action(NONE)
  {
    auto& ts = ekat::TestSession::get();
    if (ts.flags["c"]) {
      m_baseline_action = COMPARE;
    }
    else if (ts.flags["g"]) {
      m_baseline_action = GENERATE;
    }
    else if (ts.flags["n"]) {
      m_baseline_action = NONE;
    }
    m_baseline_path = ts.params["b"];

    EKAT_REQUIRE_MSG( !(m_baseline_action != NONE && m_baseline_path == ""),
                      "Unit test flags problem: baseline actions were requested but no baseline path was provided");

    std::string baseline_name = m_baseline_path + "/" + m_test_name;

    if (m_baseline_action == COMPARE) {
      m_ifile.open(baseline_name,std::ios::binary);
      EKAT_REQUIRE_MSG(m_ifile.good(), "Missing baselines: " + baseline_name + "\n");
    }
    else if (m_baseline_action == GENERATE) {
      m_ofile.open(baseline_name,std::ios::binary);
      EKAT_REQUIRE_MSG(m_ofile.good(), "Coult not open baseline file for write: " + baseline_name + "\n");
    }
  }

  ~UnitBase() = default;

  std::mt19937_64 get_engine()
  {
    if (m_baseline_action != COMPARE) {
      // We can use any seed
      int seed;
      auto engine = setup_random_test(nullptr, &seed);
      if (m_baseline_action == GENERATE) {
        // Write the seed
        impl::write_scalars(m_ofile,seed);
      }
      return engine;
    }
    else {
      // Read the seed
      int seed;
      impl::read_scalars(m_ifile,seed);
      return setup_random_test(seed, nullptr);
    }
  }
};

}  // namespace scream

#endif // SCREAM_PHYSICS_TEST_DATA_HPP
