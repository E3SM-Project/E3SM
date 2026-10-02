#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {

TEST_CASE("field_broadcast") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);

  const int nlevs = 7;
  const int ncols = 3;

  auto grid = create_point_grid("physics",ncols*comm.size(),nlevs,comm);
  const int nlocal = grid->get_num_local_dofs();

  auto layout_cl = grid->get_3d_scalar_layout(LEV);
  auto layout_c  = grid->get_2d_scalar_layout();
  auto layout_l  = grid->get_vertical_layout(LEV);

  // The target; its content must not matter
  Field tgt(FieldIdentifier("tgt",layout_cl,m,grid->name()),true);

  int seed = get_random_test_seed(&comm);

  auto& diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params;
  params.set("grid_name",grid->name());
  REQUIRE_THROWS (diag_factory.create("FieldBroadcast",comm,params,grid)); // No 'field_name'/'target_name'
  params.set<std::string>("field_name","x");
  REQUIRE_THROWS (diag_factory.create("FieldBroadcast",comm,params,grid)); // No 'target_name'
  params.set<std::string>("target_name","tgt");

  // Broadcast along COL or LEV
  for (const auto& xl : {layout_c,layout_l}) {
    const bool is_col = xl.tags()[0]==COL;

    Field x(FieldIdentifier("x",xl,kg,grid->name()),true);
    randomize_uniform(x,seed++,1,2);

    auto diag = diag_factory.create("FieldBroadcast",comm,params,grid);
    diag->set_input_field(x);
    diag->set_input_field(tgt);
    diag->initialize();

    auto d = diag->get();
    const auto& dfid = d.get_header().get_identifier();
    REQUIRE (dfid.get_layout()==layout_cl);
    REQUIRE (dfid.get_units()==kg);
    REQUIRE (not d.has_valid_mask());

    diag->compute(util::TimeStamp());
    d.sync_to_host();
    x.sync_to_host();
    auto d_v = d.get_view<const Real**,Host>();
    for (int i=0; i<nlocal; ++i) {
      for (int k=0; k<nlevs; ++k) {
        const Real expected = is_col ? x.get_view<const Real*,Host>()(i)
                                     : x.get_view<const Real*,Host>()(k);
        REQUIRE (d_v(i,k)==expected);
      }
    }
  }

  // Masked input: the mask is broadcasted too, and invalid entries are filled
  {
    Field x(FieldIdentifier("x",layout_l,kg,grid->name()),true);
    randomize_uniform(x,seed++,1,2);
    x.create_valid_mask(Field::MaskInit::Valid);
    x.get_valid_mask().subfield(0,0).deep_copy(0);

    auto diag = diag_factory.create("FieldBroadcast",comm,params,grid);
    diag->set_input_field(x);
    diag->set_input_field(tgt);
    diag->initialize();

    auto d = diag->get();
    REQUIRE (d.has_valid_mask());

    diag->compute(util::TimeStamp());
    d.sync_to_host();
    d.get_valid_mask().sync_to_host();
    x.sync_to_host();
    auto d_v = d.get_view<const Real**,Host>();
    auto m_v = d.get_valid_mask().get_view<const int**,Host>();
    auto x_v = x.get_view<const Real*,Host>();
    for (int i=0; i<nlocal; ++i) {
      for (int k=0; k<nlevs; ++k) {
        if (k==0) {
          REQUIRE (m_v(i,k)==0);
          REQUIRE (d_v(i,k)==constants::fill_value<Real>);
        } else {
          REQUIRE (m_v(i,k)==1);
          REQUIRE (d_v(i,k)==x_v(k));
        }
      }
    }
  }

  // Broadcasting to a layout that is not a superset of the input one is an error
  {
    Field y(FieldIdentifier("x",layout_c,kg,grid->name()),true);
    Field tgt_l(FieldIdentifier("tgt",layout_l,m,grid->name()),true);
    auto diag = diag_factory.create("FieldBroadcast",comm,params,grid);
    diag->set_input_field(y);
    diag->set_input_field(tgt_l);
    REQUIRE_THROWS (diag->initialize());
  }
}

}  // namespace scream
