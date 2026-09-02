#include <catch2/catch.hpp>
#include <numeric>

#include "share/field/field_identifier.hpp"
#include "share/field/field.hpp"
#include "share/field/field_group.hpp"

namespace {

TEST_CASE("field_group") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FID = FieldIdentifier;
  using FL  = FieldLayout;

  constexpr int ncols = 10;
  constexpr int ndims = 4;
  constexpr int nlevs = 8;

  FID fid ("V",FL({COL,CMP,LEV},{ncols,ndims,nlevs}),none,"the_grid");
  Field f (fid);
  f.allocate_view();

  // Create group and set subfields
  FieldGroup g("my_group","my_grid");

  std::vector<std::string> names = {"V_0","V_1","V_2","V_3"};
  g.set_monolithic_field(f,names,1,0);

  // Check const cloning
  auto cg= g.get_const();
  REQUIRE (cg.monolithic_field().is_read_only());
  REQUIRE (cg.individual_fields().size()==g.individual_fields().size());
  REQUIRE (cg.monolithic_field().get_internal_view_data<const Real>()==
            g.monolithic_field().get_internal_view_data<const Real>());
  for (int i=0; i<ndims; ++i) {
    const auto&  f =  g.individual_fields().at("V_"+std::to_string(i));
    const auto& cf = cg.individual_fields().at("V_"+std::to_string(i));
    REQUIRE ( f.get_internal_view_data<const Real>()==
             cf.get_internal_view_data<const Real>());
  }
}

} // anonymous namespace
