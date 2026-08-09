#include "catch2/catch.hpp"

#include "share/diagnostics/register_diagnostics.hpp"

#include "share/io/eamxx_diag_bank.hpp"
#include "share/io/eamxx_io_utils.hpp"
#include "share/grid/point_grid.hpp"

namespace scream {

namespace {

// A diag computing three fields in one pass, which is the case a satellite
// simulator (COSP, say) presents: one expensive call, a whole suite of outputs,
// and no reason for any one of them to be the primary.
class MultiOut : public AbstractDiagnostic
{
public:
  MultiOut (const ekat::Comm& comm, const ekat::ParameterList& params,
            const std::shared_ptr<const AbstractGrid>& grid)
   : AbstractDiagnostic(comm,params,grid)
  {
    m_field_in_names.push_back("T_mid");
  }

  std::string name () const override { return "MultiOut"; }

  static std::vector<std::string> outputs () { return {"sim_a","sim_b","sim_c"}; }

protected:
  void initialize_impl () override {
    using namespace ekat::units;
    const auto layout = m_grid->get_2d_scalar_layout();
    for (const auto& n : outputs()) {
      Field f (FieldIdentifier(n,layout,K,m_grid->name()));
      f.allocate_view();
      add_output(f);
    }
  }

  void compute_impl () override {
    // Each output gets a different value, so a caller mixing them up shows up.
    double v = 1;
    for (const auto& n : outputs()) {
      get(n).deep_copy(v++);
    }
  }
};

// Same, but with a primary output assigned the legacy way, plus extras. This is
// the shape a single-output diag grows into when it starts computing more.
class PrimaryPlusExtra : public AbstractDiagnostic
{
public:
  PrimaryPlusExtra (const ekat::Comm& comm, const ekat::ParameterList& params,
                    const std::shared_ptr<const AbstractGrid>& grid)
   : AbstractDiagnostic(comm,params,grid)
  {
    m_field_in_names.push_back("T_mid");
  }

  std::string name () const override { return "PrimaryPlusExtra"; }

  static std::vector<std::string> outputs () { return {"main_out","side_out"}; }

protected:
  void initialize_impl () override {
    using namespace ekat::units;
    const auto layout = m_grid->get_2d_scalar_layout();

    m_diagnostic_output = Field(FieldIdentifier("main_out",layout,K,m_grid->name()));
    m_diagnostic_output.allocate_view();

    Field side (FieldIdentifier("side_out",layout,K,m_grid->name()));
    side.allocate_view();
    add_output(side);
  }

  void compute_impl () override {
    m_diagnostic_output.deep_copy(1.0);
    get("side_out").deep_copy(2.0);
  }
};

std::shared_ptr<FieldManager>
create_test_fm (const std::shared_ptr<const AbstractGrid>& grid,
                const util::TimeStamp& t0)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  const auto layout = grid->get_3d_scalar_layout(LEV);
  auto fm = std::make_shared<FieldManager>(grid);

  Field f (FieldIdentifier("T_mid",layout,K,grid->name()));
  f.allocate_view();
  f.deep_copy(1.0);
  f.get_header().get_tracking().update_time_stamp(t0);
  fm->add_field(f);

  return fm;
}

} // anonymous namespace

TEST_CASE("multi_output_diags")
{
  ekat::Comm comm(MPI_COMM_WORLD);

  register_diagnostics();
  register_multi_output_diagnostic<MultiOut>("MultiOut",MultiOut::outputs());
  register_multi_output_diagnostic<PrimaryPlusExtra>("PrimaryPlusExtra",
                                                     PrimaryPlusExtra::outputs());

  util::TimeStamp t0 ({2024,1,1},{0,0,0});
  auto grid = create_point_grid("physics",3*comm.size(),10,comm);
  auto fm   = create_test_fm(grid,t0);

  SECTION ("the abstract class hands back every output") {
    auto diag = create_diagnostic("sim_b",grid);
    REQUIRE (std::dynamic_pointer_cast<MultiOut>(diag)!=nullptr);

    diag->set_input_field(fm->get_field("T_mid"));
    diag->initialize();

    REQUIRE (diag->get_output_names()==MultiOut::outputs());
    REQUIRE (diag->get_outputs().size()==3);
    for (const auto& n : MultiOut::outputs()) {
      REQUIRE (diag->has_output(n));
      REQUIRE (diag->get(n).name()==n);
    }
    REQUIRE (not diag->has_output("no_such_output"));
    REQUIRE_THROWS (diag->get("no_such_output"));

    // With no primary assigned, the first output declared plays that role.
    REQUIRE (diag->get().name()=="sim_a");
    REQUIRE (diag->get().is_aliasing(diag->get("sim_a")));
  }

  SECTION ("one compute() timestamps all the outputs") {
    auto diag = create_diagnostic("sim_a",grid);
    diag->set_input_field(fm->get_field("T_mid"));
    diag->initialize();
    diag->compute(t0);

    double expected = 1;
    for (const auto& n : MultiOut::outputs()) {
      auto f = diag->get(n);
      REQUIRE (f.get_header().get_tracking().get_time_stamp()==t0);
      f.sync_to_host();
      auto f_v = f.get_view<const Real*,Host>();
      REQUIRE (f_v(0)==Real(expected));
      ++expected;
    }
  }

  SECTION ("a legacy single-output diag is unchanged") {
    // Nothing declares an output for it, get() is the only accessor it needs,
    // and it answers to its own name only.
    std::map<std::string,diag_ptr_t> repo;
    auto tree = build_diag_tree(lower_to_diag_spec("T_mid.at(lev=3)"),grid,*fm,repo);

    REQUIRE (tree.eval_order.size()==1);
    const auto& diag = tree.eval_order.back();
    REQUIRE (diag->get_output_names().size()==1);
    REQUIRE (diag->has_output(diag->get().name()));
    REQUIRE (tree.top.is_aliasing(diag->get()));
  }

  SECTION ("an output name resolves to the diag computing it") {
    std::map<std::string,diag_ptr_t> repo;
    auto tree = build_diag_tree(lower_to_diag_spec("sim_c"),grid,*fm,repo);

    REQUIRE (tree.eval_order.size()==1);
    REQUIRE (std::dynamic_pointer_cast<MultiOut>(tree.eval_order.back())!=nullptr);

    // The tree hands back the output that was asked for, not the primary.
    REQUIRE (tree.top.name()=="sim_c");
    REQUIRE (tree.top.is_aliasing(tree.eval_order.back()->get("sim_c")));
  }

  SECTION ("sibling outputs share one diag") {
    std::map<std::string,diag_ptr_t> repo;
    auto t1 = build_diag_tree(lower_to_diag_spec("sim_a"),grid,*fm,repo);
    auto t2 = build_diag_tree(lower_to_diag_spec("sim_c"),grid,*fm,repo);

    // One instance: computing the suite twice would be the whole cost of it.
    REQUIRE (t1.eval_order.back()==t2.eval_order.back());
    REQUIRE (t2.eval_order.size()==1);
    REQUIRE (not t2.top.is_aliasing(t1.top));
  }

  SECTION ("outputs can be used inside an expression") {
    std::map<std::string,diag_ptr_t> repo;
    auto tree = build_diag_tree(lower_to_diag_spec("sim_a + sim_b"),grid,*fm,repo);

    // The multi-output diag, then the op consuming two of its outputs.
    REQUIRE (tree.eval_order.size()==2);
    REQUIRE (std::dynamic_pointer_cast<MultiOut>(tree.eval_order.front())!=nullptr);

    for (auto& d : tree.eval_order) {
      d->compute(t0);
    }
    tree.top.sync_to_host();
    auto top_v = tree.top.get_view<const Real*,Host>();
    REQUIRE (top_v(0)==Real(3.0)); // 1 + 2
  }

  SECTION ("a primary plus extras works too") {
    std::map<std::string,diag_ptr_t> repo;
    auto t1 = build_diag_tree(lower_to_diag_spec("main_out"),grid,*fm,repo);
    auto t2 = build_diag_tree(lower_to_diag_spec("side_out"),grid,*fm,repo);

    REQUIRE (t1.eval_order.back()==t2.eval_order.back());

    const auto& diag = t1.eval_order.back();
    REQUIRE (diag->get().name()=="main_out");            // the assigned primary
    REQUIRE (diag->get_output_names()==PrimaryPlusExtra::outputs());
    REQUIRE (t2.top.is_aliasing(diag->get("side_out")));
  }

  SECTION ("the bank builds the suite once and writes each output") {
    DiagBank bank(grid,*fm);
    bank.add("sim_a");
    bank.add("renamed:=sim_c");
    bank.build();

    // Two entries, one diag.
    REQUIRE (bank.entries().size()==2);
    REQUIRE (bank.eval_order().size()==1);
    REQUIRE (bank.entries().at("sim_a").diag==bank.entries().at("renamed").diag);

    REQUIRE (bank.field("sim_a").name()=="sim_a");
    REQUIRE (bank.field("renamed").name()=="renamed");
    REQUIRE (bank.field("renamed").is_aliasing(bank.eval_order().back()->get("sim_c")));
  }

  SECTION ("declaring the same output for two diags is an error") {
    REQUIRE_THROWS (
      register_multi_output_diagnostic<MultiOut>("SomethingElse",{"sim_a"}));
  }
}

} // namespace scream
