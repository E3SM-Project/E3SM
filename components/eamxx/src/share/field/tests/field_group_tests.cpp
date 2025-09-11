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

  FID fid ("V",FL({COL,CMP,LEV},{ncols,ndims,nlevs}),Units::nondimensional(),"the_grid");
  Field f (fid);
  f.allocate_view();

  FieldGroupInfo info("G");
  info.m_monolithic_allocation = true;
  std::vector<Field> f_i;
  for (int i=0; i<ndims; ++i) {
    f_i.push_back(f.get_component(i));
    info.m_fields_names.push_back(f_i[i].name());
    info.m_subview_dim = 1;
    info.m_subview_idx[f_i[i].name()] = i;
  }

  // Create group and set subfields
  FieldGroup g(info);
  g.m_monolithic_field = std::make_shared<Field>(f);
  for (int i=0; i<ndims; ++i) {
    g.m_individual_fields["G_"+std::to_string(i)] = std::make_shared<Field>(f_i[i]);
  }

  // Check const cloning
  auto cg= g.get_const();
  REQUIRE (cg.m_monolithic_field->is_read_only());
  REQUIRE (cg.m_individual_fields.size()==g.m_individual_fields.size());
  REQUIRE (*cg.m_info==*g.m_info);
  REQUIRE (cg.m_monolithic_field->get_internal_view_data<const Real>()==
            g.m_monolithic_field->get_internal_view_data<const Real>());
  for (int i=0; i<ndims; ++i) {
    const auto&  f =  *g.m_individual_fields.at("G_"+std::to_string(i));
    const auto& cf = *cg.m_individual_fields.at("G_"+std::to_string(i));
    REQUIRE ( f.get_internal_view_data<const Real>()==
             cf.get_internal_view_data<const Real>());
  }
}

} // anonymous namespace
