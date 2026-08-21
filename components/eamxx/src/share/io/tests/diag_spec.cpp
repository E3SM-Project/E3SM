#include "catch2/catch.hpp"

#include "share/io/eamxx_diag_spec.hpp"

namespace scream {

namespace {

// Shorthands, to keep the checks below readable.
std::string key (const DiagSpec& s) { return s.factory_key; }
std::string par (const DiagSpec& s, const std::string& p) {
  return s.params.get<std::string>(p);
}
const DiagSpec& kid (const DiagSpec& s, const int i) { return s.children[i]; }

} // anonymous namespace

TEST_CASE("bare names")
{
  // A name we know nothing about stays unresolved: only the FM can tell whether
  // it is a model field or a diag registered under that name.
  const auto f = lower_to_diag_spec("T_mid");
  REQUIRE (f.is_leaf());
  REQUIRE (f.names.canonical=="T_mid");
  REQUIRE (f.children.empty());

  // Names that are a diagnostic all by themselves.
  const auto wp = lower_to_diag_spec("LiqWaterPath");
  REQUIRE (key(wp)=="WaterPath");
  REQUIRE (par(wp,"water_kind")=="Liq");

  const auto np = lower_to_diag_spec("RainNumberPath");
  REQUIRE (key(np)=="NumberPath");
  REQUIRE (par(np,"number_kind")=="Rain");

  const auto pf = lower_to_diag_spec("precip_total_surf_mass_flux");
  REQUIRE (key(pf)=="precip_surf_mass_flux");
  REQUIRE (par(pf,"precip_type")=="total");

  const auto vf = lower_to_diag_spec("ZonalVapFlux");
  REQUIRE (key(vf)=="VaporFlux");
  REQUIRE (par(vf,"wind_component")=="Zonal");

  // PotentialTemperature defaults to the total kind.
  REQUIRE (par(lower_to_diag_spec("PotentialTemperature"),"temperature_kind")=="Tot");
  REQUIRE (par(lower_to_diag_spec("LiqPotentialTemperature"),"temperature_kind")=="Liq");

  const auto vl = lower_to_diag_spec("geopotential_int");
  REQUIRE (key(vl)=="VerticalLayer");
  REQUIRE (par(vl,"diag_name")=="geopotential");
  REQUIRE (par(vl,"vert_location")=="int");

  const auto dz = lower_to_diag_spec("dz");
  REQUIRE (key(dz)=="VerticalLayer");
  REQUIRE (par(dz,"diag_name")=="dz");
  REQUIRE (par(dz,"vert_location")=="mid");

  REQUIRE_THROWS (lower_to_diag_spec("AeroComCldTop")); // disabled for now
}

TEST_CASE("field at")
{
  const auto l = lower_to_diag_spec("T_mid.at(lev=10)");
  REQUIRE (key(l)=="FieldAtLevel");
  REQUIRE (par(l,"vertical_location")=="lev_10");
  REQUIRE (par(l,"field_name")=="T_mid");
  REQUIRE (kid(l,0).names.canonical=="T_mid");

  REQUIRE (par(lower_to_diag_spec("T_mid.at(lev='top')"),"vertical_location")=="model_top");
  REQUIRE (par(lower_to_diag_spec("T_mid.at(lev='bot')"),"vertical_location")=="model_bot");

  const auto p = lower_to_diag_spec("T_mid.at(p=500, units='hPa')");
  REQUIRE (key(p)=="FieldAtPressureLevel");
  REQUIRE (par(p,"pressure_value")=="500");
  REQUIRE (par(p,"pressure_units")=="hPa");

  const auto z = lower_to_diag_spec("T_mid.at(z=20, units='m', ref='sealevel')");
  REQUIRE (key(z)=="FieldAtHeight");
  REQUIRE (par(z,"height_value")=="20");
  REQUIRE (par(z,"height_units")=="m");
  REQUIRE (par(z,"surface_reference")=="sealevel");

  REQUIRE_THROWS (lower_to_diag_spec("T_mid.at(lev='middle')"));       // bad location
  REQUIRE_THROWS (lower_to_diag_spec("T_mid.at(p=500, units='KPa')")); // bad units
  REQUIRE_THROWS (lower_to_diag_spec("T_mid.at(z=1, units='km', ref='surface')"));
  REQUIRE_THROWS (lower_to_diag_spec("T_mid.at(z=1, units='m', ref='ground')"));
  REQUIRE_THROWS (lower_to_diag_spec("T_mid.at(p=500)"));              // missing units
  REQUIRE_THROWS (lower_to_diag_spec("T_mid.at(p=500, unit='hPa')"));  // typo'd keyword
  REQUIRE_THROWS (lower_to_diag_spec("T_mid.at(lev=1, p=500)"));       // ambiguous
  REQUIRE_THROWS (lower_to_diag_spec("T_mid.at(10)"));                 // positional
}

TEST_CASE("contractions")
{
  const auto h = lower_to_diag_spec("T_mid.horiz_avg()");
  REQUIRE (key(h)=="HorizAvg");
  REQUIRE (par(h,"field_name")=="T_mid");

  const auto va = lower_to_diag_spec("T_mid.vert_avg()");
  REQUIRE (key(va)=="VertContract");
  REQUIRE (par(va,"contract_method")=="avg");
  REQUIRE_FALSE (va.params.isParameter("weighting_method"));

  const auto vs = lower_to_diag_spec("T_mid.vert_sum(weight='dp')");
  REQUIRE (par(vs,"contract_method")=="sum");
  REQUIRE (par(vs,"weighting_method")=="dp");

  const auto za = lower_to_diag_spec("T_mid.zonal_avg(bins=36)");
  REQUIRE (key(za)=="ZonalAvg");
  REQUIRE (par(za,"number_of_zonal_bins")=="36");

  const auto d = lower_to_diag_spec("T_mid.derivative(dx='p')");
  REQUIRE (key(d)=="VertDerivative");
  REQUIRE (par(d,"derivative_method")=="p");

  REQUIRE_THROWS (lower_to_diag_spec("T_mid.vert_avg(weight='dx')"));
  REQUIRE_THROWS (lower_to_diag_spec("T_mid.zonal_avg(bins=1.5)"));
  REQUIRE_THROWS (lower_to_diag_spec("T_mid.derivative(dx='x')"));
}

TEST_CASE("binary ops")
{
  const auto b = lower_to_diag_spec("a - b");
  REQUIRE (key(b)=="BinaryOp");
  REQUIRE (par(b,"binary_op")=="minus");
  REQUIRE (par(b,"arg1")=="a");
  REQUIRE (par(b,"arg2")=="b");
  REQUIRE (b.children.size()==2);

  REQUIRE (par(lower_to_diag_spec("a + b"),"binary_op")=="plus");
  REQUIRE (par(lower_to_diag_spec("a * b"),"binary_op")=="times");
  REQUIRE (par(lower_to_diag_spec("a / b"),"binary_op")=="over");

  // Precedence comes from the parser, and parens mean what they say. Neither
  // was true of the regexes, where the left operand was simply greedy.
  const auto def = lower_to_diag_spec("a - b / c");
  REQUIRE (par(def,"binary_op")=="minus");
  REQUIRE (par(kid(def,1),"binary_op")=="over");

  const auto grp = lower_to_diag_spec("(a - b) / c");
  REQUIRE (par(grp,"binary_op")=="over");
  REQUIRE (par(kid(grp,0),"binary_op")=="minus");

  // Comparisons are only meaningful inside .where().
  REQUIRE_THROWS (lower_to_diag_spec("a > b"));
}

TEST_CASE("conditional sampling")
{
  const auto c = lower_to_diag_spec("T_mid.where(qv > 0)");
  REQUIRE (key(c)=="ConditionalSampling");
  REQUIRE (par(c,"field_name")=="T_mid");
  REQUIRE (par(c,"condition_lhs")=="qv");
  REQUIRE (par(c,"condition_cmp")=="gt");
  REQUIRE (par(c,"condition_rhs")=="0");

  REQUIRE (par(lower_to_diag_spec("T_mid.where(qv >= 0)"),"condition_cmp")=="ge");
  REQUIRE (par(lower_to_diag_spec("T_mid.where(qv <= 0)"),"condition_cmp")=="le");
  REQUIRE (par(lower_to_diag_spec("T_mid.where(qv < 0)"),"condition_cmp")=="lt");
  REQUIRE (par(lower_to_diag_spec("T_mid.where(qv == 0)"),"condition_cmp")=="eq");
  REQUIRE (par(lower_to_diag_spec("T_mid.where(qv != 0)"),"condition_cmp")=="ne");

  // 'mask' asks for the mask itself, so it is not an operand...
  const auto m = lower_to_diag_spec("mask.where(bt_prod < 0)");
  REQUIRE (par(m,"field_name")=="mask");
  REQUIRE (m.children.size()==1); // only the lhs

  // ...and 'lev' on the lhs is a reserved name, not a field.
  const auto lv = lower_to_diag_spec("T_mid.where(lev > 3)");
  REQUIRE (par(lv,"condition_lhs")=="lev");
  REQUIRE (lv.children.size()==1); // only the operand

  // The parser allows compound conditions; the diagnostic does not.
  REQUIRE_THROWS (lower_to_diag_spec("T_mid.where(qv > 0 and qc > 0)"));
  REQUIRE_THROWS (lower_to_diag_spec("T_mid.where(qv)"));
}

TEST_CASE("histogram")
{
  const auto h = lower_to_diag_spec("T_mid.histogram([0, 1, 2.5])");
  REQUIRE (key(h)=="Histogram");
  REQUIRE (par(h,"bin_configuration")=="0_1_2.5");
  REQUIRE (par(h,"field_name")=="T_mid");

  REQUIRE_THROWS (lower_to_diag_spec("T_mid.histogram([])"));
  REQUIRE_THROWS (lower_to_diag_spec("T_mid.histogram(0, 1)"));
}

TEST_CASE("tendencies")
{
  const auto p = lower_to_diag_spec("T_mid.prev()");
  REQUIRE (key(p)=="FieldPrev");
  REQUIRE (par(p,"field_name")=="T_mid");

  // 'dt' on the rhs of a division is reserved: X/dt is a rate, not a quotient
  // by a field named dt. The regexes had to fake this with an ordering rule.
  const auto o = lower_to_diag_spec("T_mid/dt");
  REQUIRE (key(o)=="FieldOverDt");
  REQUIRE (par(o,"field_name")=="T_mid");

  // .tend() is shorthand for (X - X.prev())/dt, and lowers to exactly that.
  const auto t = lower_to_diag_spec("T_mid.tend()");
  REQUIRE (key(t)=="FieldOverDt");
  const auto& diff = kid(t,0);
  REQUIRE (key(diff)=="BinaryOp");
  REQUIRE (par(diff,"binary_op")=="minus");
  REQUIRE (par(diff,"arg1")=="T_mid");
  REQUIRE (par(diff,"arg2")=="T_mid.prev()");
  REQUIRE (key(kid(diff,1))=="FieldPrev");

  // Both spellings produce the same canonical names, so they share diags.
  REQUIRE (t.names.canonical==lower_to_diag_spec("(T_mid-T_mid.prev())/dt").names.canonical);
}

TEST_CASE("nesting")
{
  const auto s = lower_to_diag_spec("T_mid.at(p=500, units='hPa').horiz_avg()");
  REQUIRE (key(s)=="HorizAvg");
  REQUIRE (key(kid(s,0))=="FieldAtPressureLevel");

  // The parent refers to its child by the child's canonical name, so that the
  // child's field can be aliased to it when the tree is instantiated.
  REQUIRE (par(s,"field_name")==kid(s,0).names.canonical);
  REQUIRE (kid(s,0).names.canonical=="T_mid.at((p=500), (units='hPa'))");

  // Bare attribute access is the same as a no-arg call.
  REQUIRE (key(lower_to_diag_spec("T_mid.prev"))=="FieldPrev");
}

TEST_CASE("names")
{
  // Internal nodes carry a canonical name only; the root also carries the name
  // it is registered under, which defaults to what the user asked for.
  const auto s = lower_to_diag_spec("a.horiz_avg() * b");
  REQUIRE (s.names.registered=="a.horiz_avg() * b");
  REQUIRE (s.names.canonical=="(a.horiz_avg()*b)");
  REQUIRE (s.names.write);
  REQUIRE (kid(s,0).names.registered=="");
  REQUIRE_FALSE (kid(s,0).names.write);

  // An intermediate declared in the yaml 'aliases' section: it must live in the
  // FM under a nice name, but must not reach the nc file.
  const auto i = lower_to_diag_spec("a.horiz_avg()","bt1",false);
  REQUIRE (i.names.registered=="bt1");
  REQUIRE_FALSE (i.names.write);
}

TEST_CASE("bad input")
{
  REQUIRE_THROWS (lower_to_diag_spec(""));
  REQUIRE_THROWS (lower_to_diag_spec("@"));
  REQUIRE_THROWS (lower_to_diag_spec("T_mid y"));         // trailing input
  REQUIRE_THROWS (lower_to_diag_spec("T_mid.bogus()"));   // unknown method
  REQUIRE_THROWS (lower_to_diag_spec("horiz_avg(T_mid)")); // not method syntax
  REQUIRE_THROWS (lower_to_diag_spec("[1,2]"));           // not a diagnostic
  REQUIRE_THROWS (lower_to_diag_spec("-T_mid"));          // no unary diags
}

} // namespace scream
