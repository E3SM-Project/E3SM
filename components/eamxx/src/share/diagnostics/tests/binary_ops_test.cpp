#include <catch2/catch.hpp>

#include <algorithm>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

namespace scream {


TEST_CASE("binary_ops") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager - single column for these tests
  const int nlevs = 201;
  const int ncols = 260;

  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  // Input (randomized) qc, qv
  FieldLayout scalar2d_layout = grid->get_3d_scalar_layout(LEV);
  FieldIdentifier qc_fid("qc", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier qv_fid("qv", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier m_fid("m", scalar2d_layout, kg, grid->name());

  Field qc(qc_fid, true);
  Field qv(qv_fid, true);
  Field m(m_fid, true);

  // Random number generator seed
  int seed = get_random_test_seed();

  // Construct the Diagnostics
  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params_plus, params_times, params_scl_scl;
  REQUIRE_THROWS(diag_factory.create("BinaryOp", comm,
                                     params_plus, grid));  // No 'arg1', 'arg2', or 'binary_op'

  // Set time for qc and randomize its values
  qc.get_header().get_tracking().update_time_stamp(t0);
  qv.get_header().get_tracking().update_time_stamp(t0);
  m.get_header().get_tracking().update_time_stamp(t0);
  randomize_uniform(qc, seed++, 0, 200);
  randomize_uniform(qv, seed++, 0, 200);
  randomize_uniform(m,  seed++, 0, 200);

  // Create random mask fields for qc/qv
  qc.create_valid_mask();
  qv.create_valid_mask();
  randomize_uniform(qc.get_valid_mask(),seed++);
  randomize_uniform(qv.get_valid_mask(),seed++);

  // Create and set up the diagnostic
  params_plus.set("grid_name", grid->name());
  params_plus.set<std::string>("arg1", "qc");
  params_plus.set<std::string>("arg2", "qv");
  params_plus.set<std::string>("binary_op", "starcraft");
  REQUIRE_THROWS(diag_factory.create("BinaryOp", comm, params_plus, grid));  // Unknown binary_op
  params_plus.set<std::string>("binary_op", "plus");
  auto plus_diag = diag_factory.create("BinaryOp", comm, params_plus, grid);

  params_times.set("grid_name", grid->name());
  params_times.set<std::string>("arg1", "m");
  params_times.set<std::string>("arg2", "gravit");
  params_times.set<std::string>("binary_op", "times");
  auto prod_diag = diag_factory.create("BinaryOp", comm, params_times, grid);

  params_scl_scl.set("grid_name", grid->name());
  params_scl_scl.set<std::string>("arg1", "Avogad");
  params_scl_scl.set<std::string>("arg2", "boltz");
  params_scl_scl.set<std::string>("binary_op", "times");
  auto prod_scl_scl = diag_factory.create("BinaryOp", comm, params_scl_scl, grid);

  plus_diag->set_input_field(qc);
  plus_diag->set_input_field(qv);
  plus_diag->initialize();

  prod_diag->set_input_field(m);
  prod_diag->initialize();

  prod_scl_scl->initialize();

  // Run diag
  plus_diag->compute(t0);
  prod_diag->compute(t0);
  prod_scl_scl->compute(t0);

  auto plus_diag_f = plus_diag->get(); plus_diag_f.sync_to_host();
  auto prod_diag_f = prod_diag->get(); prod_diag_f.sync_to_host();
  auto rgas_diag_f = prod_scl_scl->get(); rgas_diag_f.sync_to_host();

  // Check that the output fields have the right values
  const auto g = physics::Constants<Real>::gravit.value;

  // field*constant diag
  m.scale(g);
  REQUIRE (not prod_diag_f.has_valid_mask());
  REQUIRE (views_are_equal(prod_diag_f,m));

  // field+field diag
  REQUIRE (plus_diag_f.has_valid_mask());
  auto tgt_mask = qc.get_valid_mask().clone(CloneFlags::CopyData);
  tgt_mask.scale(qv.get_valid_mask());
  REQUIRE (views_are_equal(tgt_mask,plus_diag_f.get_valid_mask()));

  qv.update(qc,1,1);
  auto tol = std::numeric_limits<Real>::epsilon()*1e2*0;
  REQUIRE (views_are_equal(plus_diag_f,qv,tol));

  // constant*constant diag
  REQUIRE (rgas_diag_f.get_view<Real,Host>()()==physics::Constants<Real>::dictionary().at("Rgas").value);

  // Two fields on different grids are an error
  Field qv_other(FieldIdentifier("qv", scalar2d_layout, kg / kg, "other_grid"), true);
  auto plus_other = diag_factory.create("BinaryOp", comm, params_plus, grid);
  plus_other->set_input_field(qc);
  plus_other->set_input_field(qv_other);
  REQUIRE_THROWS(plus_other->initialize());
}

TEST_CASE("binary_op_constants_do_not_shadow_fields")
{
  // No constant has the name of a field: omega is the vertical velocity
  const auto& dict = physics::Constants<Real>::dictionary();
  REQUIRE (dict.find("omega")==dict.end());
  REQUIRE (dict.find("omega_earth")!=dict.end());

  ekat::Comm comm(MPI_COMM_WORLD);
  auto grid = create_point_grid("physics",16,8,comm);
  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList p;
  p.set("grid_name", grid->name());
  p.set<std::string>("arg1","omega");
  p.set<std::string>("arg2","qv");
  p.set<std::string>("binary_op","times");
  auto diag = diag_factory.create("BinaryOp",comm,p,grid);
  const auto& inputs = diag->get_input_fields_names();
  REQUIRE (std::find(inputs.begin(),inputs.end(),"omega")!=inputs.end());
  REQUIRE (std::find(inputs.begin(),inputs.end(),"qv")!=inputs.end());
}

TEST_CASE("binary_op_literals")
{
  // A numeric literal has no units of its own: in a sum or difference it takes
  // the other operand's, in a product or ratio it is a pure number. A named
  // constant keeps its units in every case.
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});
  auto grid = create_point_grid("physics",16,8,comm);
  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  const auto lay = grid->get_3d_scalar_layout(LEV);
  Field T (FieldIdentifier("T_mid", lay, K,     grid->name()), true);
  Field qv(FieldIdentifier("qv",    lay, kg/kg, grid->name()), true);
  for (auto f : {&T,&qv}) f->get_header().get_tracking().update_time_stamp(t0);
  int seed = get_random_test_seed();
  randomize_uniform(T,  seed++, 200, 300);
  randomize_uniform(qv, seed++, 0, 0.02);

  auto make = [&](const std::string& a1, const std::string& op, const std::string& a2) {
    ekat::ParameterList p;
    p.set("grid_name", grid->name());
    p.set<std::string>("arg1",a1);
    p.set<std::string>("arg2",a2);
    p.set<std::string>("binary_op",op);
    auto d = diag_factory.create("BinaryOp",comm,p,grid);
    for (const auto& n : d->get_input_fields_names()) {
      d->set_input_field(n=="T_mid" ? T : qv);
    }
    d->initialize();
    return d;
  };
  auto units_of = [&](const std::string& a1, const std::string& op, const std::string& a2) {
    return make(a1,op,a2)->get().get_header().get_identifier().get_units();
  };

  // Literals in a sum or difference take the field's units, on either side
  REQUIRE (units_of("T_mid","minus","273.15")==K);
  REQUIRE (units_of("273.15","plus","T_mid")==K);
  REQUIRE (units_of("300","minus","T_mid")==K);
  REQUIRE (units_of("qv","plus","-1e-3")==kg/kg);
  // ...and are pure numbers in a product or ratio
  REQUIRE (units_of("T_mid","times","2")==K);
  REQUIRE (units_of("T_mid","over","100")==K);
  REQUIRE (units_of("2","over","T_mid")==none/K);
  // Two literals have no units to take
  REQUIRE (units_of("2","plus","3")==none);

  // Named constants keep their units: Tmelt is in K
  REQUIRE (units_of("T_mid","minus","Tmelt")==K);
  REQUIRE_THROWS (make("qv","minus","Tmelt"));
  REQUIRE_THROWS (make("Tmelt","plus","qv"));
  // ...and a literal next to a named constant takes the constant's units
  REQUIRE (units_of("Tmelt","minus","40")==K);

  // The value is unchanged by the units bookkeeping
  auto d = make("T_mid","minus","273.15");
  d->compute(t0);
  Field expected = T.clone("expected");
  expected.update(T,1,0,-273.15);
  REQUIRE (views_are_equal(d->get(),expected));
  auto d2 = make("300","minus","T_mid");
  d2->compute(t0);
  expected.update(T,-1,0,300);
  REQUIRE (views_are_equal(d2->get(),expected));
}

TEST_CASE("binary_op_constant_fields")
{
  // An operation on constants alone is a rank-0 field that other operations use as a scalar
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using C = physics::Constants<Real>;

  ekat::Comm comm(MPI_COMM_WORLD);
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});
  auto grid = create_point_grid("physics",16,8,comm);
  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  const auto lay = grid->get_3d_scalar_layout(LEV);
  Field T (FieldIdentifier("T_mid", lay, K, grid->name()), true);
  T.get_header().get_tracking().update_time_stamp(t0);
  randomize_uniform(T, get_random_test_seed(), 200, 300);
  // Rank 0, but not made of constants (think ps.mean('col'))
  Field s (FieldIdentifier("s", FieldLayout({},{}), K, grid->name()), true);
  s.get_header().get_tracking().update_time_stamp(t0);
  s.deep_copy(100000);

  auto make = [&](const std::string& a1, const std::string& op, const std::string& a2,
                  const std::vector<Field>& inputs) {
    ekat::ParameterList p;
    p.set("grid_name", grid->name());
    p.set<std::string>("arg1",a1);
    p.set<std::string>("arg2",a2);
    p.set<std::string>("binary_op",op);
    auto d = diag_factory.create("BinaryOp",comm,p,grid);
    for (const auto& n : d->get_input_fields_names()) {
      for (const auto& f : inputs) {
        if (f.name()==n) d->set_input_field(f);
      }
    }
    d->initialize();
    return d;
  };
  // The Quantity<Real> in the header of a constant (key as in binary_op.cpp), if any
  auto constant_of = [](const Field& f) {
    const auto& h = f.get_header();
    const bool is = h.has_extra_data("constant_value");
    return std::make_pair(is,is ? h.get_extra_data<Quantity<Real>>("constant_value") : Quantity<Real>());
  };

  // Rair/Cpair: a constant, with units (none here)
  auto kappa = make("Rair","over","Cpair",{})->get();
  const Real kappa_v = C::Rair.value / C::Cpair.value;
  REQUIRE (kappa.get_header().get_identifier().get_layout().rank()==0);
  auto [is_const,q] = constant_of(kappa);
  REQUIRE (is_const);
  REQUIRE (q.value==Approx(kappa_v));
  REQUIRE (q.units==none);
  REQUIRE_FALSE (constant_of(s).first);

  // It is a scalar for others, on either side
  Field expected = T.clone("expected");
  auto d1 = make("T_mid","times",kappa.name(),{T,kappa});
  d1->compute(t0);
  expected.deep_copy(T);
  expected.scale(kappa_v);
  REQUIRE (views_are_equal(d1->get(),expected));
  auto d2 = make(kappa.name(),"over","T_mid",{T,kappa});
  d2->compute(t0);
  expected.deep_copy(kappa_v);
  expected.scale_inv(T);
  REQUIRE (views_are_equal(d2->get(),expected));
  REQUIRE (d2->get().get_header().get_identifier().get_units()==none/K);

  // A constant keeps its units, unlike a literal
  REQUIRE_THROWS (make("T_mid","plus",kappa.name(),{T,kappa}));
  // A constant of constants is a constant
  auto kappa2 = make(kappa.name(),"times","2",{kappa})->get();
  REQUIRE (constant_of(kappa2).second.value==Approx(2*kappa_v));

  // Literals alone make a literal: invalid units in the header, none on the field. So T_mid+(2+3) is in K.
  auto five = make("2","plus","3",{})->get();
  REQUIRE (constant_of(five).second.value==5);
  REQUIRE (constant_of(five).second.units==Units::invalid());
  REQUIRE (five.get_header().get_identifier().get_units()==none);
  auto d3 = make("T_mid","plus",five.name(),{T,five});
  d3->compute(t0);
  expected.deep_copy(T);
  expected.update(T,1,0,5);
  REQUIRE (views_are_equal(d3->get(),expected));
  REQUIRE (d3->get().get_header().get_identifier().get_units()==K);
  // With a named constant, it is not a literal
  REQUIRE (constant_of(make("Tmelt","plus","3",{})->get()).second.units==K);

  // Next to a non-constant field with the same layout, a constant stays a field: s/kappa is a true
  // division, not s*(1/kappa), which differs in the last digits for s=100000
  auto d4 = make("s","over",kappa.name(),{s,kappa});
  d4->compute(t0);
  Field expected0 = s.clone("expected0");
  expected0.deep_copy(s);
  expected0.scale_inv(kappa);
  REQUIRE (views_are_equal(d4->get(),expected0));
  REQUIRE_FALSE (constant_of(d4->get()).first);

  // A rank-0 field that is not a constant is not a scalar
  REQUIRE_THROWS_WITH (make("T_mid","plus","s",{T,s}), Catch::Contains("broadcast_like"));
  REQUIRE_THROWS_WITH (make("s","times","T_mid",{T,s}), Catch::Contains("broadcast_like"));
  REQUIRE_THROWS_WITH (make("1","over","0",{}), Catch::Contains("Division by zero"));
}

}  // namespace scream
