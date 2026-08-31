#include "share/io/eamxx_dexpr_diags.hpp"

#include <dexpr/ast.hpp>
#include <dexpr/lexer.hpp>
#include <dexpr/parser.hpp>
#include <dexpr/supported_functions.hpp>
#include <dexpr/tokens.hpp>

#include <ekat_assert.hpp>

#include <map>
#include <string>
#include <type_traits>
#include <vector>

namespace scream {

namespace {

using dexpr::TokenTypes;
namespace ast = dexpr::ast;

// ---------------------------------------------------------------------------
// Small helpers over the AST
// ---------------------------------------------------------------------------

template<typename Node>
const Node* as (const ast::Expression& e) {
  return e.visit([](const auto& node) -> const Node* {
    if constexpr (std::is_same_v<std::decay_t<decltype(node)>,Node>) {
      return &node;
    } else {
      return nullptr;
    }
  });
}

[[noreturn]] void unsupported (const std::string& what, const ast::Expression& e)
{
  EKAT_ERROR_MSG (
      "Error! Unsupported expression: " + what + ".\n"
      " - subexpression: " + ast::to_string(e) + "\n"
      " - for the diagnostics EAMxx can build from an expression, see\n"
      "   components/eamxx/docs/user/diags/expressions.md\n");
  throw std::logic_error("unreachable"); // EKAT_ERROR_MSG already threw
}

// A method call parses as a FuncExpression whose 'function' is a Dot binary
// expression, so the receiver and method name have to be dug back out.
struct Call {
  const ast::Expression* receiver = nullptr;
  std::string name;
  std::vector<const ast::Expression*> positional;
  std::map<std::string,const ast::Expression*> keywords;
};

// Returns false if 'e' is not a method call at all.
bool as_method_call (const ast::Expression& e, Call& call)
{
  const auto* fn = as<ast::FuncExpression>(e);
  if (fn==nullptr) {
    return false;
  }
  const auto* dot = as<ast::BinaryExpression>(*fn->function);
  if (dot==nullptr or dot->op!=TokenTypes::Dot) {
    return false;
  }
  const auto* name = as<ast::Identifier>(*dot->right);
  if (name==nullptr) {
    return false;
  }

  call.receiver = dot->left.get();
  call.name     = name->value;
  for (const auto& arg : fn->args) {
    // A keyword argument is an Assign whose lhs is a name. Anything else is
    // positional, including 'f(1=2)', which validate_calls() rejects on arity.
    const auto* assign = as<ast::BinaryExpression>(*arg);
    const auto* kw = assign!=nullptr && assign->op==TokenTypes::Assign
                   ? as<ast::Identifier>(*assign->left)
                   : nullptr;
    if (kw!=nullptr) {
      call.keywords[kw->value] = assign->right.get();
    } else {
      call.positional.push_back(arg.get());
    }
  }
  return true;
}

// The name a subexpression is known by. Anything but a plain field name comes
// back parenthesized, e.g. "(qc+qv)" -- the name create_diagnostic() resolves
// later, and the resulting field's name. dexpr renders it, so it round-trips.
std::string name_of (const ast::Expression& e) {
  return ast::to_string(e);
}

// dexpr renders literals: a FloatLiteral keeps a double, not its lexeme, and
// dexpr prints the shortest form that reads back.
std::string literal_str (const ast::Expression& e)
{
  if (as<ast::IntegerLiteral>(e) or as<ast::FloatLiteral>(e)) {
    return ast::to_string(e);
  }
  if (const auto* s = as<ast::StringLiteral>(e)) {
    return s->value;
  }
  unsupported("expected a literal",e);
}

int int_arg (const ast::Expression& e, const std::string& ctx)
{
  const auto* i = as<ast::IntegerLiteral>(e);
  if (i==nullptr) {
    // A negative index parses as unary minus; the only place we accept one,
    // since no diagnostic negates a field.
    const auto* u = as<ast::UnaryExpression>(e);
    if (u!=nullptr and u->op==TokenTypes::Minus) {
      const auto* inner = as<ast::IntegerLiteral>(*u->right);
      if (inner!=nullptr) {
        return -inner->value;
      }
    }
    unsupported(ctx + " must be an integer",e);
  }
  return i->value;
}

std::string string_arg (const ast::Expression& e, const std::string& ctx)
{
  const auto* s = as<ast::StringLiteral>(e);
  if (s==nullptr) {
    unsupported(ctx + " must be a quoted string",e);
  }
  return s->value;
}

const ast::Expression* keyword (const Call& call, const std::string& name)
{
  auto it = call.keywords.find(name);
  return it==call.keywords.end() ? nullptr : it->second;
}

// ---------------------------------------------------------------------------
// Translation: one AST node -> one diagnostic
// ---------------------------------------------------------------------------
// Only the root is translated. Operands go down by name, and the customer
// resolving dependencies brings them back here one at a time.

std::string comparison_to_cmp (TokenTypes op)
{
  switch (op) {
    case TokenTypes::GreaterThan:  return "gt";
    case TokenTypes::GreaterEqual: return "ge";
    case TokenTypes::LessThan:     return "lt";
    case TokenTypes::LessEq:       return "le";
    case TokenTypes::Equal:        return "eq";
    case TokenTypes::NotEqual:     return "ne";
    default:                       return "";
  }
}

std::string binary_op_to_diag_op (TokenTypes op)
{
  switch (op) {
    case TokenTypes::Plus:     return "plus";
    case TokenTypes::Minus:    return "minus";
    case TokenTypes::Asterisk: return "times";
    case TokenTypes::Slash:    return "over";
    default:                   return "";
  }
}

void translate_binary (const ast::BinaryExpression& e, const ast::Expression& self,
                       std::string& diag_name, ekat::ParameterList& params)
{
  const auto op = binary_op_to_diag_op(e.op);
  if (op=="") {
    if (comparison_to_cmp(e.op)!="") {
      unsupported("a comparison is only meaningful inside where(..)",self);
    }
    if (e.op==TokenTypes::Dot) {
      unsupported("attribute access with no call",self);
    }
    unsupported("no diagnostic implements the operator '" +
                dexpr::binary_op_to_string(e.op) + "'",self);
  }

  diag_name = "BinaryOp";
  params.set<std::string>("arg1",name_of(*e.left));
  params.set<std::string>("arg2",name_of(*e.right));
  params.set<std::string>("binary_op",op);
}

// ---------------------------------------------------------------------------
// The EAMxx vocabulary
// ---------------------------------------------------------------------------
// Named after xarray where there is an honest analogue. Lives here, not in
// share/dexpr: dexpr owns the grammar, we own the vocabulary.
//
// Sits above translate_call() because the two are a pair -- every entry here
// needs a case there. Nothing enforces that in the type system, so two checks
// do: validate_registry() below, and the create_diag case that builds every
// example.
const dexpr::FunctionRegistry& eamxx_registry ()
{
  static const dexpr::FunctionRegistry reg = [] {
    dexpr::FunctionRegistry r;

    r.add({.name = "isel",
           .desc = "value at a vertical index; -1 is the bottom level",
           .min_positional = 0, .max_positional = 0,
           .keywords = {{"lev",true}},
           .example = "T_mid.isel(lev=10)"});
    r.add({.name = "interp",
           .desc = "value interpolated to a pressure or height coordinate",
           .min_positional = 0, .max_positional = 0,
           .keywords = {{"p_mid",false},{"z_mid",false},{"units",false},{"reference",false}},
           .example = "T_mid.interp(p_mid=500,units='hPa')"});
    r.add({.name = "mean",
           .desc = "average over a dimension ('col' or 'lev')",
           .min_positional = 1, .max_positional = 1,
           .keywords = {{"weights",false}},
           .example = "T_mid.mean('lev',weights='dp')"});
    r.add({.name = "sum",
           .desc = "sum over a dimension ('lev')",
           .min_positional = 1, .max_positional = 1,
           .keywords = {{"weights",false}},
           .example = "T_mid.sum('lev',weights='dz')"});
    r.add({.name = "where",
           .desc = "keep values where a condition holds",
           .min_positional = 1, .max_positional = 1,
           .keywords = {},
           .example = "T_mid.where(qv>0.01)"});
    r.add({.name = "differentiate",
           .desc = "vertical derivative w.r.t. 'p_mid' or 'z_mid'",
           .min_positional = 1, .max_positional = 1,
           .keywords = {},
           .example = "T_mid.differentiate('p_mid')"});
    r.add({.name = "histogram",
           .desc = "counts per bin, given the bin edges",
           .min_positional = 0, .max_positional = 0,
           .keywords = {{"bins",true}},
           .example = "T_mid.histogram(bins=[0,1,2])"});
    r.add({.name = "zonal_mean",
           .desc = "average within latitude bands",
           .min_positional = 0, .max_positional = 0,
           .keywords = {{"bins",true}},
           .example = "T_mid.zonal_mean(bins=20)"});
    r.add({.name = "shift",
           .desc = "value at a previous time step; only time=1 is available",
           .min_positional = 0, .max_positional = 0,
           .keywords = {{"time",true}},
           .example = "T_mid.shift(time=1)"});
    r.add({.name = "over_dt",
           .desc = "value divided by the time step",
           .min_positional = 0, .max_positional = 0,
           .keywords = {},
           .example = "T_mid.over_dt()"});
    r.add({.name = "tend",
           .desc = "tendency over the time step, i.e. (X-X.shift(time=1)).over_dt()",
           .min_positional = 0, .max_positional = 0,
           .keywords = {},
           .example = "T_mid.tend()"});

    // Each example must parse, match the spec that declared it, and call it.
    // A wrong arity or keyword fails here, not when a user first writes it.
    dexpr::validate_registry(r);
    return r;
  }();
  return reg;
}

void translate_call (const Call& call, const ast::Expression& self,
                     std::string& diag_name, ekat::ParameterList& params)
{
  const auto operand = name_of(*call.receiver);
  params.set<std::string>("field_name",operand);

  if (call.name=="isel") {
    // xarray indexing: negative counts from the end. Only -1 resolves without
    // knowing the level count. lev=0 is the model top.
    const auto idx = int_arg(*keyword(call,"lev"),"the 'lev' index");
    std::string location;
    if (idx==-1) {
      location = "model_bot";
    } else {
      EKAT_REQUIRE_MSG (idx>=0,
          "Error! Invalid vertical index in isel().\n"
          " - expression: " + ast::to_string(self) + "\n"
          " - input: " + std::to_string(idx) + "\n"
          " - only -1 (the bottom level) is supported as a negative index,\n"
          "   since the others cannot be resolved without the level count.\n");
      location = "lev_" + std::to_string(idx);
    }
    diag_name = "FieldAtLevel";
    params.set<std::string>("vertical_location",location);
    return;
  }

  if (call.name=="interp") {
    const auto* plev = keyword(call,"p_mid");
    const auto* z    = keyword(call,"z_mid");
    EKAT_REQUIRE_MSG ((plev!=nullptr) != (z!=nullptr),
        "Error! interp() takes exactly one of 'p_mid' or 'z_mid'.\n"
        " - expression: " + ast::to_string(self) + "\n");

    if (plev!=nullptr) {
      const auto* units = keyword(call,"units");
      EKAT_REQUIRE_MSG (keyword(call,"reference")==nullptr,
          "Error! 'reference' applies to interp(z_mid=..), not interp(p_mid=..).\n"
          " - expression: " + ast::to_string(self) + "\n");
      diag_name = "FieldAtPressureLevel";
      params.set<std::string>("pressure_value",literal_str(*plev));
      params.set<std::string>("pressure_units",
          units==nullptr ? std::string("Pa") : string_arg(*units,"'units'"));
    } else {
      const auto* units = keyword(call,"units");
      const auto* ref   = keyword(call,"reference");
      diag_name = "FieldAtHeight";
      params.set<std::string>("height_value",literal_str(*z));
      params.set<std::string>("height_units",
          units==nullptr ? std::string("m") : string_arg(*units,"'units'"));
      params.set<std::string>("surface_reference",
          ref==nullptr ? std::string("sealevel") : string_arg(*ref,"'reference'"));
    }
    return;
  }

  if (call.name=="mean" or call.name=="sum") {
    const auto dim = string_arg(*call.positional[0],"the dimension");
    const auto* weights = keyword(call,"weights");

    if (dim=="col") {
      EKAT_REQUIRE_MSG (call.name=="mean",
          "Error! Only an average is available over 'col'.\n"
          " - expression: " + ast::to_string(self) + "\n");
      EKAT_REQUIRE_MSG (weights==nullptr,
          "Error! Averaging over 'col' is always area weighted; 'weights' does not apply.\n"
          " - expression: " + ast::to_string(self) + "\n");
      diag_name = "HorizAvg";
      return;
    }

    EKAT_REQUIRE_MSG (dim=="lev",
        "Error! Unknown dimension in " + call.name + "().\n"
        " - expression: " + ast::to_string(self) + "\n"
        " - input: '" + dim + "'\n"
        " - valid dimensions: 'col', 'lev'\n");

    diag_name = "VertContract";
    params.set<std::string>("contract_method",call.name=="mean" ? "avg" : "sum");
    if (weights!=nullptr) {
      params.set<std::string>("weighting_method",string_arg(*weights,"'weights'"));
    }
    return;
  }

  if (call.name=="where") {
    const auto& cond = *call.positional[0];
    const auto* cmp = as<ast::BinaryExpression>(cond);
    EKAT_REQUIRE_MSG (cmp!=nullptr and comparison_to_cmp(cmp->op)!="",
        "Error! where() takes a single comparison.\n"
        " - expression: " + ast::to_string(self) + "\n"
        " - condition: " + ast::to_string(cond) + "\n"
        " - valid comparisons: >, >=, <, <=, ==, !=\n"
        " - note: 'and'/'or' are not supported; chain where(..) calls instead.\n");

    diag_name = "ConditionalSampling";
    params.set<std::string>("condition_lhs",name_of(*cmp->left));
    params.set<std::string>("condition_cmp",comparison_to_cmp(cmp->op));
    params.set<std::string>("condition_rhs",name_of(*cmp->right));
    return;
  }

  if (call.name=="differentiate") {
    const auto wrt = string_arg(*call.positional[0],"the coordinate");
    EKAT_REQUIRE_MSG (wrt=="p_mid" or wrt=="z_mid",
        "Error! Unknown coordinate in differentiate().\n"
        " - expression: " + ast::to_string(self) + "\n"
        " - input: '" + wrt + "'\n"
        " - valid coordinates: 'p_mid', 'z_mid'\n");
    diag_name = "VertDerivative";
    params.set<std::string>("derivative_method",wrt=="p_mid" ? "p" : "z");
    return;
  }

  if (call.name=="histogram") {
    const auto* bins = as<ast::ArrayExpression>(*keyword(call,"bins"));
    EKAT_REQUIRE_MSG (bins!=nullptr and bins->elements.size()>=2,
        "Error! histogram() takes 'bins', an array of at least two edges.\n"
        " - expression: " + ast::to_string(self) + "\n");
    // The diag re-splits on '_', so edges must survive: no exponents or signs.
    std::string config;
    for (size_t i=0; i<bins->elements.size(); ++i) {
      const auto edge = literal_str(*bins->elements[i]);
      EKAT_REQUIRE_MSG (edge.find_first_not_of("0123456789.")==std::string::npos,
          "Error! Histogram bin edges must be non-negative decimals.\n"
          " - expression: " + ast::to_string(self) + "\n"
          " - edge: " + edge + "\n");
      config += (i==0 ? "" : "_") + edge;
    }
    diag_name = "Histogram";
    params.set<std::string>("bin_configuration",config);
    return;
  }

  if (call.name=="zonal_mean") {
    const auto bins = int_arg(*keyword(call,"bins"),"'bins'");
    EKAT_REQUIRE_MSG (bins>0,
        "Error! zonal_mean() needs a positive number of bins.\n"
        " - expression: " + ast::to_string(self) + "\n"
        " - bins: " + std::to_string(bins) + "\n");
    diag_name = "ZonalAvg";
    params.set<std::string>("number_of_zonal_bins",std::to_string(bins));
    return;
  }

  if (call.name=="shift") {
    // xarray's sign convention: a positive shift moves data forward, so
    // shift(time=1) is the value from one step BACK. A negative shift would be
    // the value from a step ahead, which a running model cannot know -- worth
    // saying out loud, since isel(lev=-1) makes -1 look like "one back".
    const auto n = int_arg(*keyword(call,"time"),"'time'");
    EKAT_REQUIRE_MSG (n==1,
        "Error! Only shift(time=1) is available.\n"
        " - expression: " + ast::to_string(self) + "\n"
        " - input: " + std::to_string(n) + "\n"
        " - a positive shift looks BACK in time, as in xarray, so time=1 is the\n"
        "   previous step. Negative would look ahead, which we cannot do, and\n"
        "   only one step of history is kept.\n");
    diag_name = "FieldPrev";
    return;
  }

  if (call.name=="over_dt") {
    diag_name = "FieldOverDt";
    return;
  }

  if (call.name=="tend") {
    // Shorthand expanded here, not in the grammar. The operand is the
    // difference expression, which comes back through create_diagnostic().
    // NOTE: the inner parens around 'time=1' are how dexpr renders a keyword
    //       argument, so this matches what name_of() gives for the same
    //       expression written out by hand, and the two share one diag.
    diag_name = "FieldOverDt";
    params.set<std::string>("field_name",
        "(" + operand + "-" + operand + ".shift((time=1)))");
    return;
  }

  // validate_calls() rejects unknown names, so this means the registry and the
  // cases above have drifted apart.
  EKAT_ERROR_MSG (
      "Error! Internal error: '" + call.name + "' is registered but not translated.\n"
      " - expression: " + ast::to_string(self) + "\n");
}

} // anonymous namespace

std::shared_ptr<AbstractDiagnostic>
dexpr_create_diagnostic (const std::string& expr,
                         const std::shared_ptr<const AbstractGrid>& grid)
{
  ast::ExprPtr root;
  try {
    dexpr::parser::Parser parser {dexpr::Lexer{expr}};
    root = parser.parse();
  } catch (const dexpr::parser::ParserError& e) {
    // Matched no legacy pattern and is not a legal identifier, so it can only
    // have been meant as an expression. Say where it went wrong.
    EKAT_ERROR_MSG (
        "Error! '" + expr + "' is neither a registered diagnostic nor a valid expression.\n" +
        std::string(e.what()));
  }

  // A bare identifier is a diag class name or a typo; the caller's to resolve.
  if (as<ast::Identifier>(*root)) {
    return nullptr;
  }

  dexpr::validate_calls(*root,eamxx_registry());

  std::string diag_name;
  ekat::ParameterList params(expr);
  params.set("grid_name",grid->name());

  if (Call call; as_method_call(*root,call)) {
    translate_call(call,*root,diag_name,params);
  } else if (const auto* bin = as<ast::BinaryExpression>(*root)) {
    translate_binary(*bin,*root,diag_name,params);
  } else if (const auto* un = as<ast::UnaryExpression>(*root)) {
    unsupported("no diagnostic implements the unary operator '" +
                dexpr::unary_op_to_string(un->op) + "'",*root);
  } else if (as<ast::FuncExpression>(*root)) {
    unsupported("EAMxx functions are written as methods, X.f(..), not f(X)",*root);
  } else {
    unsupported("an expression must produce a field",*root);
  }

  // Name the field after the request, so customers find what they asked for.
  params.set<std::string>("output_field_name",expr);
  // Mark how it resolved: an expression is not a usable NetCDF variable name,
  // so whoever writes it must require an output name.
  params.set<bool>("from_expression",true);

  return DiagnosticFactory::instance().create(diag_name,grid->get_comm(),params,grid);
}

std::vector<std::string> dexpr_diagnostic_examples ()
{
  std::vector<std::string> out;
  for (const auto& entry : eamxx_registry()) {
    out.push_back(entry.second.example);
  }
  return out;
}

} // namespace scream
