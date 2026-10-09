#include <catch2/catch.hpp>

#include "share/field/field_visit.hpp"
#include "share/field/field_identifier.hpp"
#include "share/field/field.hpp"

#include <numeric>
#include <type_traits>

namespace {

using namespace scream;
using namespace ekat::units;
using namespace ShortFieldTagsNames;

using details::view_as;
using details::visit_field_views;

Field make_field (const std::string& name, const std::vector<FieldTag>& tags,
                  const std::vector<int>& dims, const DataType dt)
{
  Field f (FieldIdentifier(name,{tags,dims},m,"some_grid",dt));
  f.allocate_view();
  return f;
}

// Whether the view visited for this (rank-1 or rank-2) field is strided or not
bool visits_strided_view (const Field& f)
{
  bool strided = false;
  visit_field_views(f.get_header().get_identifier().get_layout().rank(),
    [&](const auto& v) {
      using layout_t = typename std::decay_t<decltype(v)>::array_layout;
      strided = std::is_same_v<layout_t,Kokkos::LayoutStride>;
    },
    view_as<double>(f));
  return strided;
}

TEST_CASE ("visit_field_views") {
  SECTION ("rank") {
    // The callback must receive a view of the right rank and extents, for all ranks
    const std::vector<FieldTag> all_tags = {COL,CMP,GP,LEV,ILEV,LEVP};
    for (int rank=0; rank<=6; ++rank) {
      std::vector<FieldTag> tags (all_tags.begin(),all_tags.begin()+rank);
      std::vector<int> dims (rank);
      std::iota(dims.begin(),dims.end(),2);
      auto f = make_field("f",tags,dims,DataType::DoubleType);

      int num_calls = 0;
      visit_field_views(rank,
        [&](const auto& v) {
          ++num_calls;
          REQUIRE (static_cast<int>(v.rank())==rank);
          for (int i=0; i<rank; ++i)
            REQUIRE (v.extent_int(i)==dims[i]);
        },
        view_as<double>(f));
      REQUIRE (num_calls==1);
    }

    // Unsupported ranks are rejected
    auto noop = [](auto) {};
    REQUIRE_THROWS (details::dispatch_on_rank(-1,noop));
    REQUIRE_THROWS (details::dispatch_on_rank(Field::MaxRank+1,noop));
  }

  SECTION ("view_layout") {
    // Subfields that are not contiguous in the parent must get strided views
    auto f = make_field("f",{COL,CMP},{3,4},DataType::DoubleType);
    REQUIRE_FALSE (visits_strided_view(f));
    REQUIRE_FALSE (visits_strided_view(f.subfield(0,1)));
    REQUIRE (visits_strided_view(f.subfield(1,1)));
  }

  SECTION ("multiple_fields") {
    // Views must be passed in the order of the view_as tags, with the requested
    // scalar type, and mixing strided and non-strided views must work.
    auto parent = make_field("parent",{COL,CMP},{3,4},DataType::DoubleType);
    auto y = parent.subfield(1,2); // strided
    auto x = make_field("x",{COL},{3},DataType::FloatType);
    auto k = make_field("k",{COL},{3},DataType::IntType);

    parent.deep_copy(0.0);

    int num_calls = 0;
    visit_field_views(1,
      [&](const auto& yv, const auto& xv, const auto& kv) {
        using y_t = std::decay_t<decltype(yv)>;
        using x_t = std::decay_t<decltype(xv)>;
        using k_t = std::decay_t<decltype(kv)>;
        static_assert (std::is_same_v<typename y_t::value_type,double>);
        static_assert (std::is_same_v<typename x_t::value_type,const float>);
        static_assert (std::is_same_v<typename k_t::value_type,const int>);

        // The callback is instantiated for all layout combos, so layouts must be checked at runtime
        constexpr bool y_strided = std::is_same_v<typename y_t::array_layout,Kokkos::LayoutStride>;
        constexpr bool x_strided = std::is_same_v<typename x_t::array_layout,Kokkos::LayoutStride>;
        constexpr bool k_strided = std::is_same_v<typename k_t::array_layout,Kokkos::LayoutStride>;
        REQUIRE (y_strided);
        REQUIRE_FALSE (x_strided);
        REQUIRE_FALSE (k_strided);

        ++num_calls;
        Kokkos::deep_copy(yv,7.0);
      },
      view_as<double>(y), view_as<const float>(x), view_as<const int>(k));
    REQUIRE (num_calls==1);

    // Writing through the strided view must only affect the subfield's slice
    parent.sync_to_host();
    auto pv = parent.get_view<const double**,Host>();
    for (int i=0; i<3; ++i) {
      for (int j=0; j<4; ++j) {
        REQUIRE (pv(i,j)==(j==2 ? 7.0 : 0.0));
      }
    }
  }
}

} // anonymous namespace
