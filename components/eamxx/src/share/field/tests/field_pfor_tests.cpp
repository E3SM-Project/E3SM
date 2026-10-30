#include <catch2/catch.hpp>

#include "share/field/field_pfor.hpp"
#include "share/field/field_identifier.hpp"
#include "share/field/field.hpp"

#include <numeric>

namespace {

using namespace scream;
using namespace ekat::units;
using namespace ShortFieldTagsNames;

using details::view_as;
using details::field_pfor;

Field make_field (const std::string& name, const std::vector<FieldTag>& tags,
                  const std::vector<int>& dims, const DataType dt = DataType::DoubleType)
{
  Field f (FieldIdentifier(name,{tags,dims},m,"some_grid",dt));
  f.allocate_view();
  return f;
}

// NOTE: kernels must be functors with a templated operator() (not generic lambdas),
//       since nvcc does not allow generic extended lambdas

// f += g, and g = weighted sum of the indices (to check that we get the indices in the right order)
struct AccumulateAndIndexKernel {
  template<typename FView, typename GView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const FView& fv, const GView& gv, const Idx&... i) const {
    static_assert (sizeof...(i)==FView::rank());
    fv(i...) += gv(i...);

    double w = 1;
    double s = 0;
    ((s += w*i, w *= 10), ...);
    gv(i...) = s;
  }
};

// y = alpha*x + k
struct MixedKernel {
  double alpha;

  template<typename YView, typename XView, typename KView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const YView& yv, const XView& xv, const KView& kv, const Idx&... i) const {
    yv(i...) = alpha*xv(i...) + kv(i...);
  }
};

// A kernel that treats rank-2 fields in a special way
struct RankDependentKernel {
  template<typename FView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const FView& fv, const Idx&... i) const {
    if constexpr (sizeof...(i)==2) {
      const int idx[] = {static_cast<int>(i)...};
      fv(i...) = 100*idx[0] + idx[1];
    }
  }
};

struct NoopKernel {
  template<typename V1, typename V2, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const V1&, const V2&, const Idx&...) const {}
};

TEST_CASE ("field_pfor") {
  SECTION ("rank") {
    // The same, rank-agnostic kernel must work for all ranks. In particular, it must
    // visit each entry exactly once (also for rank 0), and pass the right indices.
    const std::vector<FieldTag> all_tags = {COL,CMP,GP,LEV,ILEV,LEVP};
    for (int rank=0; rank<=6; ++rank) {
      std::vector<FieldTag> tags (all_tags.begin(),all_tags.begin()+rank);
      std::vector<int> dims (rank);
      std::iota(dims.begin(),dims.end(),2);
      auto f = make_field("f",tags,dims);
      auto g = make_field("g",tags,dims);
      f.deep_copy(0.0);
      g.deep_copy(1.0);

      field_pfor(AccumulateAndIndexKernel(), view_as<double>(f), view_as<double>(g));

      f.sync_to_host();
      g.sync_to_host();
      const int size = f.get_header().get_identifier().get_layout().size();
      auto fd = f.get_internal_view_data<const double,Host>();
      auto gd = g.get_internal_view_data<const double,Host>();
      for (int n=0; n<size; ++n) {
        REQUIRE (fd[n]==1.0);

        // Recompute indices from the flat (LayoutRight) index
        int rem = n;
        double w = 1, s = 0;
        std::vector<int> idx (rank);
        for (int k=rank-1; k>=0; --k) {
          idx[k] = rem % dims[k];
          rem /= dims[k];
        }
        for (int k=0; k<rank; ++k) {
          s += w*idx[k];
          w *= 10;
        }
        REQUIRE (gd[n]==s);
      }
    }
  }

  SECTION ("multiple_fields") {
    // Mix strided/non-strided views, different scalar types, and captured scalars
    auto parent = make_field("parent",{COL,CMP},{3,4},DataType::DoubleType);
    auto y = parent.subfield(1,2); // strided
    auto x = make_field("x",{COL},{3},DataType::FloatType);
    auto k = make_field("k",{COL},{3},DataType::IntType);
    parent.deep_copy(0.0);
    x.deep_copy(2.0f);
    k.deep_copy(3);

    field_pfor(MixedKernel{10.0}, view_as<double>(y), view_as<const float>(x), view_as<const int>(k));

    // Only the subfield's slice must have been touched
    parent.sync_to_host();
    auto pv = parent.get_view<const double**,Host>();
    for (int i=0; i<3; ++i) {
      for (int j=0; j<4; ++j) {
        REQUIRE (pv(i,j)==(j==2 ? 23.0 : 0.0));
      }
    }
  }

  SECTION ("rank_dependent_kernel") {
    // A kernel can still branch on the number of indices
    auto f = make_field("f",{COL,LEV},{3,4},DataType::DoubleType);
    field_pfor(RankDependentKernel(), view_as<double>(f));
    f.sync_to_host();
    auto fv = f.get_view<const double**,Host>();
    for (int i=0; i<3; ++i)
      for (int j=0; j<4; ++j)
        REQUIRE (fv(i,j)==100*i+j);
  }

  SECTION ("errors") {
    auto f = make_field("f",{COL},{3},DataType::DoubleType);
    auto g = make_field("g",{COL},{4},DataType::DoubleType);
    auto h = make_field("h",{COL,LEV},{3,4},DataType::DoubleType);
    REQUIRE_THROWS (field_pfor(NoopKernel(),view_as<double>(f),view_as<double>(g)));
    REQUIRE_THROWS (field_pfor(NoopKernel(),view_as<double>(f),view_as<double>(h)));
  }
}

} // anonymous namespace
