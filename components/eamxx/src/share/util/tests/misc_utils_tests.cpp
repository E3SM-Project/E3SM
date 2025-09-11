#include <catch2/catch.hpp>

#include "share/util/eamxx_family_tracking.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include "share/util/eamxx_utils.hpp"

TEST_CASE("fill_value") {
  using namespace scream::constants;

  // Ensure we have the SAME numerical value for both float and double
  auto fv_d = fill_value<double>;
  auto fv_f = fill_value<float>;

  REQUIRE (fv_d==fv_f);
}

TEST_CASE("contiguous_superset") {
  using namespace scream;

  std::string A = "A";
  std::string B = "B";
  std::string C = "C";
  std::string D = "D";
  std::string E = "E";
  std::string F = "F";
  std::string G = "G";

  using LOLS_type = std::list<std::list<std::string>>;

  // These three lists do not allow a superset from which they can all be
  // contiguously subviewed.
  LOLS_type lol1 = { {A,B}, {B,C}, {A,C} };
  REQUIRE(contiguous_superset(lol1).size()==0);

  // Input inner lists are not sorted
  REQUIRE_THROWS(contiguous_superset(LOLS_type{ {B,A} }));

  // The following should both allow the superset (A,B,C,D,E,F,G)
  // Note: lol3 is simply a shuffled version of lol2
  LOLS_type lol2 = { {A,B,C}, {B,C,D,E}, {C,D}, {C,D,E,F}, {D,E,F,G} };
  LOLS_type lol3 = { {D,E,F,G}, {C,D,E,F}, {A,B,C}, {C,D}, {B,C,D,E} };

  // Flipping a list is still a valid solution, so consider both tgt and its reverse.
  std::list<std::string> tgt = {A,B,C,D,E,F,G};
  std::list<std::string> tgt_rev = tgt;
  tgt_rev.reverse();

  auto superset2 = contiguous_superset(lol2);
  auto superset3 = contiguous_superset(lol3);
  REQUIRE ( (superset2==tgt || superset2==tgt_rev) );
  REQUIRE ( (superset3==tgt || superset3==tgt_rev) );
}

struct Dummy : scream::FamilyTracking<Dummy> {
  Dummy (int i) : value(i) {}
  int value;
};

TEST_CASE("family_tracking") {
  using namespace scream;

  auto child = std::make_shared<Dummy>(2);
  auto parent = std::make_shared<Dummy>(3);

  REQUIRE (not child->get_parent());
  child->create_parent_child_link(parent->shared_from_this());
  REQUIRE (child->get_parent());
  REQUIRE (child->get_parent()==parent);
  REQUIRE (parent->get_children().size()==1);
  REQUIRE (child->get_children().size()==0);
}
