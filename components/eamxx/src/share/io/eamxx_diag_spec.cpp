#include "share/io/eamxx_diag_spec.hpp"

#include <edp/ast.hpp>
#include <edp/lexer.hpp>
#include <edp/parser.hpp>
#include <edp/tokens.hpp>

#include <ekat_assert.hpp>

#include <initializer_list>
#include <map>
#include <type_traits>
#include <vector>

namespace scream {

namespace {

using edp::TokenTypes;
using edp::ast::Expression;

// ---------------------------------------------------------------------------
//  AST helpers
// ---------------------------------------------------------------------------

// edp keeps the node variant private, so Expression::visit is the only way in.
// This wraps the "is this node a T?" query that we need everywhere.
template<typename T>
const T* node_as (const Expression& e)
{
  return e.visit([](const auto& n) -> const T* {
    if constexpr (std::is_same_v<std::decay_t<decltype(n)>,T>) {
      return &n;
    } else {
      return nullptr;
    }
  });
}

// The canonical name of a (sub)expression: edp's print is fully parenthesized
// and deterministic, which is exactly what we need out of a wiring/dedup key.
std::string canonical_of (const Expression& e)
{
  return edp::ast::to_string(e);
}

// The value of a leaf token, as the plain string that diags expect in their
// params: 'hPa' -> hPa, 500 -> 500, 0.5 -> 0.5, some_field -> some_field.
std::string leaf_to_string (const Expression& e, const std::string& ctx)
{
  if (auto* s = node_as<edp::ast::StringLiteral>(e)) {
    return s->value;
  }
  if (auto* i = node_as<edp::ast::Identifier>(e)) {
    return i->value;
  }
  if (node_as<edp::ast::IntegerLiteral>(e) or node_as<edp::ast::FloatLiteral>(e)) {
    return edp::ast::to_string(e);
  }
  EKAT_ERROR_MSG (
      "Error! Expected a name, number, or quoted string in a diagnostic expression.\n"
      " - context: " + ctx + "\n"
      " - got    : " + edp::ast::to_string(e) + "\n");
  return "";
}

// ---------------------------------------------------------------------------
//  Call arguments
// ---------------------------------------------------------------------------

// A call's arguments, split into positional ones and keyword ones. edp parses
// 'k=v' as an Assign infix node, so keywords are recovered rather than lexed.
struct CallArgs {
  std::vector<const Expression*>          positional;
  std::map<std::string,const Expression*> kw;
};

CallArgs split_args (const std::vector<edp::ast::ExprPtr>& args)
{
  CallArgs out;
  for (const auto& a : args) {
    const auto* infix = node_as<edp::ast::InfixExpression>(*a);
    const auto* key   = infix && infix->op==TokenTypes::Assign
                      ? node_as<edp::ast::Identifier>(*infix->left) : nullptr;
    if (key) {
      out.kw[key->value] = infix->right.get();
    } else {
      out.positional.push_back(a.get());
    }
  }
  return out;
}

// Reject keywords we do not know about, so that typos ('unit=' for 'units=')
// fail here rather than silently doing the wrong thing.
void check_kw (const CallArgs& args, const std::string& method,
               const std::initializer_list<const char*>& allowed)
{
  for (const auto& [k,v] : args.kw) {
    bool ok = false;
    for (const auto* a : allowed) {
      ok |= (k==a);
    }
    EKAT_REQUIRE_MSG (ok,
        "Error! Unrecognized keyword argument in a diagnostic expression.\n"
        " - method : ." + method + "()\n"
        " - keyword: " + k + "\n");
  }
}

const Expression& kw_required (const CallArgs& args, const std::string& key,
                               const std::string& method)
{
  auto it = args.kw.find(key);
  EKAT_REQUIRE_MSG (it!=args.kw.end(),
      "Error! Missing required keyword argument in a diagnostic expression.\n"
      " - method : ." + method + "()\n"
      " - keyword: " + key + "\n");
  return *it->second;
}

void check_no_positional (const CallArgs& args, const std::string& method)
{
  EKAT_REQUIRE_MSG (args.positional.empty(),
      "Error! Unexpected positional argument in a diagnostic expression.\n"
      " - method: ." + method + "()\n"
      " - hint  : arguments of this method must be given as 'key=value'.\n");
}

// ---------------------------------------------------------------------------
//  Bare diagnostic names
// ---------------------------------------------------------------------------

// Names that mean a diagnostic all by themselves, with a fixed set of params.
// These used to be regexes with capture groups; spelling the combinations out
// costs a few lines and buys exact matching plus greppability.
struct BareDiag {
  const char* name;
  const char* key;
  const char* param1;   // nullptr if the diag takes no params
  const char* value1;
  const char* param2;   // nullptr if the diag takes a single param
  const char* value2;
};

constexpr BareDiag bare_diags[] = {
  {"precip_liq_surf_mass_flux",  "precip_surf_mass_flux", "precip_type",    "liq",     nullptr, nullptr},
  {"precip_ice_surf_mass_flux",  "precip_surf_mass_flux", "precip_type",    "ice",     nullptr, nullptr},
  {"precip_total_surf_mass_flux","precip_surf_mass_flux", "precip_type",    "total",   nullptr, nullptr},

  {"IceWaterPath",               "WaterPath",             "water_kind",     "Ice",     nullptr, nullptr},
  {"LiqWaterPath",               "WaterPath",             "water_kind",     "Liq",     nullptr, nullptr},
  {"RainWaterPath",              "WaterPath",             "water_kind",     "Rain",    nullptr, nullptr},
  {"RimeWaterPath",              "WaterPath",             "water_kind",     "Rime",    nullptr, nullptr},
  {"VapWaterPath",               "WaterPath",             "water_kind",     "Vap",     nullptr, nullptr},

  {"IceNumberPath",              "NumberPath",            "number_kind",    "Ice",     nullptr, nullptr},
  {"LiqNumberPath",              "NumberPath",            "number_kind",    "Liq",     nullptr, nullptr},
  {"RainNumberPath",             "NumberPath",            "number_kind",    "Rain",    nullptr, nullptr},

  {"MeridionalVapFlux",          "VaporFlux",             "wind_component", "Meridional", nullptr, nullptr},
  {"ZonalVapFlux",               "VaporFlux",             "wind_component", "Zonal",   nullptr, nullptr},

  {"PotentialTemperature",       "PotentialTemperature",  "temperature_kind", "Tot",   nullptr, nullptr},
  {"LiqPotentialTemperature",    "PotentialTemperature",  "temperature_kind", "Liq",   nullptr, nullptr},

  {"z_mid",                      "VerticalLayer",         "diag_name", "z",            "vert_location", "mid"},
  {"z_int",                      "VerticalLayer",         "diag_name", "z",            "vert_location", "int"},
  {"geopotential_mid",           "VerticalLayer",         "diag_name", "geopotential", "vert_location", "mid"},
  {"geopotential_int",           "VerticalLayer",         "diag_name", "geopotential", "vert_location", "int"},
  {"height_mid",                 "VerticalLayer",         "diag_name", "height",       "vert_location", "mid"},
  {"height_int",                 "VerticalLayer",         "diag_name", "height",       "vert_location", "int"},
  {"dz",                         "VerticalLayer",         "diag_name", "dz",           "vert_location", "mid"},
};

// ---------------------------------------------------------------------------
//  Lowering
// ---------------------------------------------------------------------------

DiagSpec lower (const Expression& e);

DiagSpec make_spec (const std::string& factory_key, const std::string& canonical)
{
  DiagSpec s;
  s.factory_key    = factory_key;
  s.params         = ekat::ParameterList(canonical);
  s.names.canonical = canonical;
  return s;
}

// A bare name: either one of the fixed diags above, or an unresolved leaf.
DiagSpec lower_bare_name (const std::string& name)
{
  EKAT_REQUIRE_MSG (name!="AeroComCldTop" and name!="AeroComCldBot",
      "Error! AeroComCld diags are disabled for now. Contact developers.\n"
      "      Some recent development made the code produce bad values,\n"
      "      even runtime aborts due to NaNs.\n"
      "      An alternative is to request variables like cdnc_at_cldtop,\n"
      "      which remain unaffected and scientifically valid.\n");

  for (const auto& b : bare_diags) {
    if (name==b.name) {
      auto s = make_spec(b.key,name);
      if (b.param1) s.params.set<std::string>(b.param1,b.value1);
      if (b.param2) s.params.set<std::string>(b.param2,b.value2);
      return s;
    }
  }

  // Unresolved leaf: a model field, or a diag registered under this name.
  return make_spec("",name);
}

// Add 'child' as an operand of 'parent', and store its canonical name in the
// given param, which is how the parent refers to it.
void add_operand (DiagSpec& parent, const std::string& param, DiagSpec child)
{
  parent.params.set<std::string>(param,child.names.canonical);
  parent.children.push_back(std::move(child));
}

DiagSpec lower_at (const Expression& operand, const CallArgs& args,
                   const std::string& canonical)
{
  check_no_positional(args,"at");
  check_kw(args,"at",{"lev","p","z","units","ref"});

  const bool has_lev = args.kw.count("lev")==1;
  const bool has_p   = args.kw.count("p")==1;
  const bool has_z   = args.kw.count("z")==1;
  EKAT_REQUIRE_MSG (int(has_lev)+int(has_p)+int(has_z) == 1,
      "Error! .at() requires exactly one of 'lev', 'p', or 'z'.\n"
      " - expression: " + canonical + "\n");

  if (has_lev) {
    auto s = make_spec("FieldAtLevel",canonical);
    const auto& lev = kw_required(args,"lev","at");
    if (node_as<edp::ast::IntegerLiteral>(lev)) {
      s.params.set<std::string>("vertical_location","lev_" + leaf_to_string(lev,".at(lev=)"));
    } else {
      const auto loc = leaf_to_string(lev,".at(lev=)");
      EKAT_REQUIRE_MSG (loc=="top" or loc=="bot",
          "Error! .at(lev=) takes an integer, or one of 'top'/'bot'.\n"
          " - input: " + loc + "\n");
      s.params.set<std::string>("vertical_location","model_" + loc);
    }
    add_operand(s,"field_name",lower(operand));
    return s;
  }

  if (has_p) {
    auto s = make_spec("FieldAtPressureLevel",canonical);
    const auto units = leaf_to_string(kw_required(args,"units","at"),".at(units=)");
    EKAT_REQUIRE_MSG (units=="hPa" or units=="mb" or units=="Pa",
        "Error! Unsupported pressure units in .at().\n"
        " - input units  : " + units + "\n"
        " - valid choices: hPa, mb, Pa\n");
    s.params.set<std::string>("pressure_value",leaf_to_string(kw_required(args,"p","at"),".at(p=)"));
    s.params.set<std::string>("pressure_units",units);
    add_operand(s,"field_name",lower(operand));
    return s;
  }

  auto s = make_spec("FieldAtHeight",canonical);
  const auto units = leaf_to_string(kw_required(args,"units","at"),".at(units=)");
  EKAT_REQUIRE_MSG (units=="m",
      "Error! Unsupported height units in .at().\n"
      " - input units  : " + units + "\n"
      " - valid choices: m\n");
  const auto ref = leaf_to_string(kw_required(args,"ref","at"),".at(ref=)");
  EKAT_REQUIRE_MSG (ref=="sealevel" or ref=="surface",
      "Error! Unsupported surface reference in .at().\n"
      " - input ref    : " + ref + "\n"
      " - valid choices: sealevel, surface\n");
  s.params.set<std::string>("height_value",leaf_to_string(kw_required(args,"z","at"),".at(z=)"));
  s.params.set<std::string>("height_units",units);
  s.params.set<std::string>("surface_reference",ref);
  add_operand(s,"field_name",lower(operand));
  return s;
}

DiagSpec lower_where (const Expression& operand, const CallArgs& args,
                      const std::string& canonical)
{
  check_kw(args,"where",{});
  EKAT_REQUIRE_MSG (args.positional.size()==1,
      "Error! .where() takes exactly one argument: a comparison.\n"
      " - expression: " + canonical + "\n");

  const auto* cond = node_as<edp::ast::InfixExpression>(*args.positional[0]);
  std::string cmp;
  if (cond) {
    switch (cond->op) {
      case TokenTypes::Equal:        cmp = "eq"; break;
      case TokenTypes::NotEqual:     cmp = "ne"; break;
      case TokenTypes::GreaterThan:  cmp = "gt"; break;
      case TokenTypes::GreaterEqual: cmp = "ge"; break;
      case TokenTypes::LessThan:     cmp = "lt"; break;
      case TokenTypes::LessEq:       cmp = "le"; break;
      default: break;
    }
  }
  EKAT_REQUIRE_MSG (not cmp.empty(),
      "Error! .where() takes a single comparison of the form 'lhs <cmp> rhs'.\n"
      " - expression   : " + canonical + "\n"
      " - valid <cmp>  : ==, !=, >, >=, <, <=\n"
      " - note: compound conditions (and/or) are not supported by the\n"
      "         ConditionalSampling diagnostic.\n");

  auto s = make_spec("ConditionalSampling",canonical);
  const auto lhs = leaf_to_string(*cond->left,".where() lhs");
  const auto rhs = leaf_to_string(*cond->right,".where() rhs");
  s.params.set<std::string>("condition_lhs",lhs);
  s.params.set<std::string>("condition_cmp",cmp);
  s.params.set<std::string>("condition_rhs",rhs);

  // 'mask' is not a field: it asks for the mask itself rather than a sampling
  // of some input, so it must not become an operand.
  const auto* in = node_as<edp::ast::Identifier>(operand);
  if (in and in->value=="mask") {
    s.params.set<std::string>("field_name","mask");
  } else {
    add_operand(s,"field_name",lower(operand));
  }

  // The lhs is a field, unless it is the reserved name 'lev'. The rhs is a
  // field only when it is not a number.
  if (node_as<edp::ast::Identifier>(*cond->left) and lhs!="lev") {
    s.children.push_back(lower(*cond->left));
  }
  if (node_as<edp::ast::Identifier>(*cond->right)) {
    s.children.push_back(lower(*cond->right));
  }
  return s;
}

DiagSpec lower_histogram (const Expression& operand, const CallArgs& args,
                          const std::string& canonical)
{
  check_kw(args,"histogram",{});
  const auto* bins = args.positional.size()==1
                   ? node_as<edp::ast::ArrayExpression>(*args.positional[0]) : nullptr;
  EKAT_REQUIRE_MSG (bins and bins->elements.size()>0,
      "Error! .histogram() takes exactly one argument: a non-empty array of bin edges.\n"
      " - expression: " + canonical + "\n"
      " - example   : X.histogram([0, 1, 2.5])\n");

  // The diag parses the edges back out of an underscore-separated string.
  std::string bin_config;
  for (const auto& b : bins->elements) {
    if (not bin_config.empty()) bin_config += "_";
    bin_config += leaf_to_string(*b,".histogram() bin edge");
  }

  auto s = make_spec("Histogram",canonical);
  s.params.set<std::string>("bin_configuration",bin_config);
  add_operand(s,"field_name",lower(operand));
  return s;
}

// X.tend() is shorthand, and is expanded here rather than being a diag of its
// own. Its canonical name is that of the expansion, so writing either form
// yields the same tree, and hence the same reusable diags.
DiagSpec lower_tend (const Expression& operand, const CallArgs& args,
                     const std::string& canonical)
{
  check_kw(args,"tend",{});
  check_no_positional(args,"tend");

  auto field = lower(operand);
  const auto cf = field.names.canonical;

  auto prev = make_spec("FieldPrev",cf + ".prev()");
  add_operand(prev,"field_name",field);

  auto diff = make_spec("BinaryOp","(" + cf + "-" + prev.names.canonical + ")");
  diff.params.set<std::string>("binary_op","minus");
  diff.params.set<std::string>("arg1",cf);
  diff.params.set<std::string>("arg2",prev.names.canonical);
  diff.children.push_back(lower(operand));
  diff.children.push_back(std::move(prev));

  auto s = make_spec("FieldOverDt","(" + diff.names.canonical + "/dt)");
  add_operand(s,"field_name",std::move(diff));

  // Silence the unused-parameter warning without losing the argument name.
  (void) canonical;
  return s;
}

DiagSpec lower_method (const Expression& operand, const std::string& method,
                       const CallArgs& args, const std::string& canonical)
{
  if (method=="at")        return lower_at(operand,args,canonical);
  if (method=="where")     return lower_where(operand,args,canonical);
  if (method=="histogram") return lower_histogram(operand,args,canonical);
  if (method=="tend")      return lower_tend(operand,args,canonical);

  if (method=="horiz_avg") {
    check_kw(args,method,{});
    check_no_positional(args,method);
    auto s = make_spec("HorizAvg",canonical);
    add_operand(s,"field_name",lower(operand));
    return s;
  }

  if (method=="vert_avg" or method=="vert_sum") {
    check_kw(args,method,{"weight"});
    check_no_positional(args,method);
    auto s = make_spec("VertContract",canonical);
    s.params.set<std::string>("contract_method",method=="vert_avg" ? "avg" : "sum");
    if (args.kw.count("weight")==1) {
      const auto w = leaf_to_string(*args.kw.at("weight"),"." + method + "(weight=)");
      EKAT_REQUIRE_MSG (w=="dp" or w=="dz",
          "Error! Unsupported weighting in ." + method + "().\n"
          " - input weight : " + w + "\n"
          " - valid choices: dp, dz\n");
      s.params.set<std::string>("weighting_method",w);
    }
    add_operand(s,"field_name",lower(operand));
    return s;
  }

  if (method=="zonal_avg") {
    check_kw(args,method,{"bins"});
    check_no_positional(args,method);
    const auto& bins = kw_required(args,"bins",method);
    EKAT_REQUIRE_MSG (node_as<edp::ast::IntegerLiteral>(bins),
        "Error! .zonal_avg(bins=) requires an integer number of bins.\n"
        " - input: " + edp::ast::to_string(bins) + "\n");
    auto s = make_spec("ZonalAvg",canonical);
    s.params.set<std::string>("number_of_zonal_bins",leaf_to_string(bins,".zonal_avg(bins=)"));
    add_operand(s,"field_name",lower(operand));
    return s;
  }

  if (method=="derivative") {
    check_kw(args,method,{"dx"});
    check_no_positional(args,method);
    const auto dx = leaf_to_string(kw_required(args,"dx",method),".derivative(dx=)");
    EKAT_REQUIRE_MSG (dx=="p" or dx=="z",
        "Error! Unsupported differentiation variable in .derivative().\n"
        " - input dx     : " + dx + "\n"
        " - valid choices: p, z\n");
    auto s = make_spec("VertDerivative",canonical);
    s.params.set<std::string>("derivative_method",dx);
    add_operand(s,"field_name",lower(operand));
    return s;
  }

  if (method=="prev" or method=="over_dt") {
    check_kw(args,method,{});
    check_no_positional(args,method);
    auto s = make_spec(method=="prev" ? "FieldPrev" : "FieldOverDt",canonical);
    add_operand(s,"field_name",lower(operand));
    return s;
  }

  EKAT_ERROR_MSG (
      "Error! Unrecognized method in a diagnostic expression.\n"
      " - method    : ." + method + "()\n"
      " - expression: " + canonical + "\n");
  return {};
}

DiagSpec lower_binary (const edp::ast::InfixExpression& infix,
                       const std::string& canonical)
{
  // X/dt means "rate of change", not "divide by a field named dt". With the old
  // regexes this had to be faked by testing _over_dt before the binary ops; here
  // 'dt' is simply a reserved name on the rhs of a division.
  const auto* rhs_id = node_as<edp::ast::Identifier>(*infix.right);
  if (infix.op==TokenTypes::Slash and rhs_id and rhs_id->value=="dt") {
    auto s = make_spec("FieldOverDt",canonical);
    add_operand(s,"field_name",lower(*infix.left));
    return s;
  }

  std::string op;
  switch (infix.op) {
    case TokenTypes::Plus:     op = "plus";  break;
    case TokenTypes::Minus:    op = "minus"; break;
    case TokenTypes::Asterisk: op = "times"; break;
    case TokenTypes::Slash:    op = "over";  break;
    default:
      EKAT_ERROR_MSG (
          "Error! Unsupported operator in a diagnostic expression.\n"
          " - expression   : " + canonical + "\n"
          " - valid choices: +, -, *, /\n"
          " - note: comparisons are only valid inside .where().\n");
  }

  auto s = make_spec("BinaryOp",canonical);
  s.params.set<std::string>("binary_op",op);
  add_operand(s,"arg1",lower(*infix.left));
  add_operand(s,"arg2",lower(*infix.right));
  return s;
}

DiagSpec lower (const Expression& e)
{
  const auto canonical = canonical_of(e);

  if (auto* id = node_as<edp::ast::Identifier>(e)) {
    return lower_bare_name(id->value);
  }

  if (auto* infix = node_as<edp::ast::InfixExpression>(e)) {
    if (infix->op==TokenTypes::Dot) {
      // Attribute access with no call, e.g. 'X.prev'. Same as a call with no
      // arguments, which is the form we document.
      const auto* method = node_as<edp::ast::Identifier>(*infix->right);
      EKAT_REQUIRE_MSG (method,
          "Error! Malformed method access in a diagnostic expression.\n"
          " - expression: " + canonical + "\n");
      return lower_method(*infix->left,method->value,CallArgs{},canonical);
    }
    return lower_binary(*infix,canonical);
  }

  if (auto* call = node_as<edp::ast::FuncExpression>(e)) {
    const auto* callee = node_as<edp::ast::InfixExpression>(*call->function);
    const auto* method = callee && callee->op==TokenTypes::Dot
                       ? node_as<edp::ast::Identifier>(*callee->right) : nullptr;
    EKAT_REQUIRE_MSG (method,
        "Error! Diagnostics are written as methods on a field.\n"
        " - expression: " + canonical + "\n"
        " - example   : X.horiz_avg(), not horiz_avg(X)\n");
    return lower_method(*callee->left,method->value,split_args(call->args),canonical);
  }

  EKAT_ERROR_MSG (
      "Error! Expression is not a valid diagnostic.\n"
      " - expression: " + canonical + "\n");
  return {};
}

} // anonymous namespace

DiagSpec lower_to_diag_spec (const std::string& expr,
                             const std::string& registered,
                             const bool write)
{
  edp::parser::Parser parser {edp::Lexer{expr}};
  const auto ast = parser.parse();
  EKAT_REQUIRE_MSG (ast!=nullptr,
      "Error! Could not parse diagnostic expression.\n"
      " - expression: " + expr + "\n");

  auto spec = lower(*ast);
  spec.names.registered = registered.empty() ? expr : registered;
  spec.names.write = write;
  return spec;
}

} // namespace scream
