#include "catch2/catch.hpp"

#include "share/diagnostics/register_diagnostics.hpp"

#include "share/io/eamxx_diag_bank.hpp"
#include "share/grid/point_grid.hpp"

namespace scream {

namespace {

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

} // anonymous namespace

TEST_CASE("diag_bank")
{
  ekat::Comm comm(MPI_COMM_WORLD);
  register_diagnostics();

  util::TimeStamp t0 ({2024,1,1},{0,0,0});
  auto grid = create_point_grid("physics",3*comm.size(),10,comm);
  auto fm   = create_test_fm(grid,t0);

  SECTION ("a request names itself when no name is given") {
    DiagBank bank(grid,*fm);
    bank.add("T_mid.at(lev=3)");
    bank.build();

    REQUIRE (bank.entries().size()==1);
    REQUIRE (bank.field("T_mid.at(lev=3)").name()=="T_mid.at(lev=3)");
    REQUIRE (bank.entries().at("T_mid.at(lev=3)").write);
  }

  SECTION ("a plain model field is a valid request") {
    DiagBank bank(grid,*fm);
    bank.add("T_mid");            // asks for the field itself...
    bank.add("temp:=T_mid");      // ...and this one just renames it
    bank.build();

    REQUIRE (bank.eval_order().empty());          // no diag needed for either
    REQUIRE (bank.entries().at("T_mid").diag==nullptr);
    REQUIRE (bank.field("T_mid").name()=="T_mid");
    REQUIRE (bank.field("temp").name()=="temp");
    REQUIRE (bank.field("temp").is_aliasing(bank.field("T_mid")));
  }

  SECTION (":= gives a request a name of its own") {
    DiagBank bank(grid,*fm);
    bank.add("T3:=T_mid.at(lev=3)");
    bank.build();

    REQUIRE (bank.field("T3").name()=="T3");
    REQUIRE (bank.eval_order().size()==1);
  }

  SECTION ("requests can refer to each other by name") {
    DiagBank bank(grid,*fm);
    bank.add("d1:=T_mid.at(lev=3)",false);  // intermediate: not written
    bank.add("d2:=d1.prev()");

    bank.build();

    REQUIRE_FALSE (bank.entries().at("d1").write);
    REQUIRE       (bank.entries().at("d2").write);

    // d1 is built, and evaluated, before the request that consumes it.
    REQUIRE (bank.eval_order().size()==2);
    REQUIRE (std::dynamic_pointer_cast<FieldAtLevel>(bank.eval_order()[0])!=nullptr);
    REQUIRE (std::dynamic_pointer_cast<FieldPrev>(bank.eval_order()[1])!=nullptr);
  }

  SECTION ("the order requests are added in does not matter") {
    DiagBank bank(grid,*fm);
    bank.add("d2:=d1.prev()");              // consumer added first...
    bank.add("d1:=T_mid.at(lev=3)",false);  // ...producer second
    bank.build();

    REQUIRE (bank.eval_order().size()==2);
    REQUIRE (std::dynamic_pointer_cast<FieldAtLevel>(bank.eval_order()[0])!=nullptr);
  }

  SECTION ("a shared sub-expression is built once and evaluated once") {
    DiagBank bank(grid,*fm);
    bank.add("a:=T_mid.at(lev=3).prev()");
    bank.add("b:=T_mid.at(lev=3) - qv.at(lev=3)");
    bank.build();

    // T_mid.at(lev=3) is common to both requests: 4 diags, not 5.
    REQUIRE (bank.eval_order().size()==4);

    // Both requests are still registered, with distinct fields.
    REQUIRE (bank.field("a").name()=="a");
    REQUIRE (bank.field("b").name()=="b");
  }

  SECTION ("everything can be evaluated in the order given") {
    DiagBank bank(grid,*fm);
    bank.add("d1:=T_mid.at(lev=3)",false);
    bank.add("d2:=d1 - qv.at(lev=3)");
    bank.build();

    for (const auto& d : bank.eval_order()) {
      d->compute(t0);
    }
    REQUIRE (bank.field("d2").get_header().get_tracking().get_time_stamp().is_valid());
  }

  SECTION ("bad requests are rejected") {
    DiagBank bank(grid,*fm);

    bank.add("T3:=T_mid.at(lev=3)");
    REQUIRE_THROWS (bank.add("T3:=qv.at(lev=3)"));  // name used twice
    REQUIRE_THROWS (bank.add("a:=b:=c"));           // malformed
    REQUIRE_THROWS (bank.add(":=x"));               // empty name

    REQUIRE_THROWS (bank.field("nope"));            // never built
  }

  SECTION ("cycles are caught") {
    DiagBank bank(grid,*fm);
    bank.add("a:=b.prev()");
    bank.add("b:=a.prev()");
    REQUIRE_THROWS (bank.build());
  }
}

} // namespace scream
