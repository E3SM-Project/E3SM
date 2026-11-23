#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/grid/point_grid.hpp"

namespace scream {

TEST_CASE("diag_hash")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);
  auto grid = create_point_grid("physics",8,4,comm);
  register_diagnostics();
  auto& factory = DiagnosticFactory::instance();

  const auto lay = grid->get_3d_scalar_layout(LEV);
  Field A  (FieldIdentifier("A",  lay,K,grid->name()), true);
  Field A2 (FieldIdentifier("A",  lay,K,grid->name()), true);  // Same name and layout, other storage
  Field B  (FieldIdentifier("B",  lay,K,grid->name()), true);
  Field uv (FieldIdentifier("uv", grid->get_3d_vector_layout(LEV,2,"cmp"),K,grid->name()), true);

  auto make = [&](const std::string& op, const std::string& list_name = "diag") {
    ekat::ParameterList params(list_name);
    params.set("grid_name", grid->name());
    params.set<std::string>("arg1","p");
    params.set<std::string>("arg2","q");
    params.set<std::string>("binary_op",op);
    return factory.create("BinaryOp",comm,params,grid);
  };
  // The hash of 'p op q', with p and q the fields given (as inputs of that name)
  auto hash_of = [&](const Field& p, const std::string& op, const Field& q) {
    auto d = make(op);
    d->set_input_field(p.alias("p"));
    d->set_input_field(q.alias("q"));
    return d->hash();
  };

  const auto h = hash_of(A,"times",B);
  // Same diag, and the same fields, whatever they are called
  REQUIRE (hash_of(A,"times",B)==h);
  REQUIRE (hash_of(A.alias("something_else"),"times",B)==h);
  // The name of the param list (the request) is not part of it
  auto d = make("times","another_name");
  d->set_input_field(A.alias("p"));
  d->set_input_field(B.alias("q"));
  REQUIRE (d->hash()==h);

  // Not the same: another operation, other fields, the same fields swapped, other storage
  REQUIRE (hash_of(A,"plus",B)!=h);
  REQUIRE (hash_of(B,"times",A)!=h);
  REQUIRE (hash_of(A,"times",A)!=h);
  REQUIRE (hash_of(A2,"times",B)!=h);

  // The components of a field (e.g., U and V of horiz_winds) are different fields...
  const auto u = uv.get_component(0);
  REQUIRE (hash_of(u,"times",B)!=hash_of(uv.get_component(1),"times",B));
  REQUIRE (hash_of(u,"times",B)!=hash_of(uv,"times",B));
  // ...and an alias of one is the same as it
  REQUIRE (hash_of(u.alias("another_name"),"times",B)==hash_of(u,"times",B));

  // All the inputs must be there
  auto incomplete = make("times");
  incomplete->set_input_field(A.alias("p"));
  REQUIRE_THROWS (incomplete->hash());
}

} // namespace scream
