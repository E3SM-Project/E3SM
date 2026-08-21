#include "catch2/catch.hpp"

#include "share/diagnostics/register_diagnostics.hpp"

#include "share/io/eamxx_diag_tree.hpp"
#include "share/grid/point_grid.hpp"

namespace scream {

namespace {

// A minimal model FM: a couple of 3d fields, with valid timestamps so that the
// diags can be computed.
std::shared_ptr<FieldManager>
create_test_fm (const std::shared_ptr<const AbstractGrid>& grid,
                const util::TimeStamp& t0)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  const auto layout = grid->get_3d_scalar_layout(LEV);
  auto fm = std::make_shared<FieldManager>(grid);

  for (const auto& name : {"T_mid","qv"}) {
    Field f (FieldIdentifier(name,layout,K,grid->name()));
    f.allocate_view();
    f.deep_copy(1.0);
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
  }
  return fm;
}

DiagTree build (const std::string& expr,
                const std::shared_ptr<const AbstractGrid>& grid,
                const FieldManager& fm,
                std::map<std::string,diag_ptr_t>& repo)
{
  return build_diag_tree(lower_to_diag_spec(expr),grid,fm,repo);
}

} // anonymous namespace

TEST_CASE("build_diag_tree")
{
  ekat::Comm comm(MPI_COMM_WORLD);
  register_diagnostics();

  util::TimeStamp t0 ({2024,1,1},{0,0,0});
  auto grid = create_point_grid("physics",3*comm.size(),10,comm);
  auto fm   = create_test_fm(grid,t0);

  SECTION ("plain field") {
    std::map<std::string,diag_ptr_t> repo;
    auto tree = build("T_mid",grid,*fm,repo);
    REQUIRE (tree.eval_order.empty());   // no diag needed
    REQUIRE (tree.top.name()=="T_mid");
  }

  SECTION ("single diag") {
    std::map<std::string,diag_ptr_t> repo;
    auto tree = build("T_mid.at(lev=3)",grid,*fm,repo);
    REQUIRE (tree.eval_order.size()==1);
    REQUIRE (std::dynamic_pointer_cast<FieldAtLevel>(tree.eval_order[0])!=nullptr);

    // The tree hands back a field named after what was requested, regardless of
    // how the diag chose to name its own output. Note this is the *registered*
    // name (here defaulted to the request), not the canonical one, which would
    // be the parenthesized 'T_mid.at((lev=3))'.
    REQUIRE (tree.top.name()=="T_mid.at(lev=3)");
  }

  SECTION ("nested diags are evaluated inputs-first") {
    std::map<std::string,diag_ptr_t> repo;
    auto tree = build("T_mid.at(lev=3) - qv.at(lev=3)",grid,*fm,repo);

    // Both operands, then the op that consumes them.
    REQUIRE (tree.eval_order.size()==3);
    REQUIRE (std::dynamic_pointer_cast<BinaryOp>(tree.eval_order.back())!=nullptr);
    for (size_t i=0; i<tree.eval_order.size()-1; ++i) {
      REQUIRE (std::dynamic_pointer_cast<FieldAtLevel>(tree.eval_order[i])!=nullptr);
    }
  }

  SECTION ("repeated sub-expressions are built once") {
    std::map<std::string,diag_ptr_t> repo;
    auto tree = build("T_mid.at(lev=3) - T_mid.at(lev=3)",grid,*fm,repo);

    // The shared operand is built once, and appears once in the eval order.
    REQUIRE (tree.eval_order.size()==2);
    REQUIRE (repo.count("T_mid.at((lev=3))")==1);
  }

  SECTION ("diags are shared across trees, but each tree evaluates them") {
    std::map<std::string,diag_ptr_t> repo;
    auto t1 = build("T_mid.at(lev=3)",grid,*fm,repo);
    auto t2 = build("T_mid.at(lev=3).prev()",grid,*fm,repo);

    // t2 reuses t1's diag...
    REQUIRE (t2.eval_order.front()==t1.eval_order.front());
    // ...but still lists it, since t2 must evaluate it too.
    REQUIRE (t2.eval_order.size()==2);
  }

  SECTION ("registered name is what the tree hands back") {
    std::map<std::string,diag_ptr_t> repo;
    auto spec = lower_to_diag_spec("T_mid.at(lev=3)","T3");
    auto tree = build_diag_tree(spec,grid,*fm,repo);
    REQUIRE (tree.top.name()=="T3");

    // Aliasing shares data, so the tree's field IS the diag's field.
    REQUIRE (tree.top.is_aliasing(tree.eval_order.back()->get()));
  }

  SECTION ("the tree can be evaluated") {
    std::map<std::string,diag_ptr_t> repo;
    auto tree = build("T_mid.at(lev=3) - qv.at(lev=3)",grid,*fm,repo);

    // Evaluating in the order the walk produced is enough: no dep is stale.
    for (auto& d : tree.eval_order) {
      d->compute(t0);
    }
    REQUIRE (tree.top.get_header().get_tracking().get_time_stamp().is_valid());
  }

  SECTION ("unknown names are rejected") {
    std::map<std::string,diag_ptr_t> repo;
    REQUIRE_THROWS (build("no_such_field.at(lev=3)",grid,*fm,repo));
    REQUIRE_THROWS (build("no_such_diag",grid,*fm,repo));
  }
}

} // namespace scream
