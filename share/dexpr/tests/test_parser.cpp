#include <dexpr/ast.hpp>
#include "catch2/catch_message.hpp"
#include <dexpr/lexer.hpp>
#include <dexpr/parser.hpp>
#include <dexpr/tokens.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>

namespace dexpr {

namespace { // anonymous

std::string parse_to_string(const std::string& input) {
  parser::Parser parser{Lexer{input}};
  const auto expr = parser.parse();
  REQUIRE(expr != nullptr);
  return ast::to_string(*expr);
}

void check_parse(const std::string& input, const std::string& expected) {
  INFO("Input: " << input);
  CHECK(parse_to_string(input) == expected);
}

// Returns the message of the ParserError raised by `input`, or "" if it parsed.
std::string parse_error(const std::string& input) {
  try {
    parse_to_string(input);
  } catch (const parser::ParserError& e) {
    return e.what();
  }
  return "";
}

void check_rejected(const std::string& input, const std::string& substring) {
  INFO("Input: " << input);
  const auto msg = parse_error(input);
  CHECK_FALSE(msg.empty());
  CHECK(msg.find(substring) != std::string::npos);
}

} // anonymous

TEST_CASE("parser: bare literals and identifiers", "[parser]") {
  check_parse("x", "x");
  check_parse("42", "42");
  check_parse("'col'", "'col'");
}


TEST_CASE("parser: field names containing digits", "[parser]") {
  check_parse("bc_a1", "bc_a1");
  check_parse("O3", "O3");
  check_parse("bc_a1 + so4_a2", "(bc_a1+so4_a2)");
  check_parse("O3.mean(dim='lev')", "O3.mean((dim='lev'))");
  check_parse("bc_a1.interp(plev=500, units='hPa')",
              "bc_a1.interp((plev=500), (units='hPa'))");
}

TEST_CASE("parser: simple infix expressions", "[parser]") {
  check_parse("x + 1", "(x+1)");
  check_parse("x * y", "(x*y)");
  check_parse("x = y", "(x=y)");
}

TEST_CASE("parser: multiplication binds tighter than addition", "[parser]") {
  check_parse("1 + 2 * 3", "(1+(2*3))");
  check_parse("1 * 2 + 3", "((1*2)+3)");
}

TEST_CASE("parser: comparison binds tighter than equality", "[parser]") {
  // LessGreater outranks Equal in the precedence table.
  check_parse("a < b == c", "((a<b)==c)");
}

TEST_CASE("parser: boolean operators can be printed", "[parser]") {
  check_parse("x and y", "(x and y)");
  check_parse("x or y", "(x or y)");
}

TEST_CASE("parser: inclusive comparisons can be printed", "[parser]") {
  check_parse("x >= y", "(x>=y)");
  check_parse("x <= y", "(x<=y)");
}

TEST_CASE("parser: division and strict comparison print correctly",
          "[parser]") {
  check_parse("x / y", "(x/y)");
  check_parse("x < y", "(x<y)");
  check_parse("x > y", "(x>y)");
}

TEST_CASE("parser: equality operators print back", "[parser]") {
  check_parse("x == y", "(x==y)");
  check_parse("x != y", "(x!=y)");
  // The printed form has to lex again as the same operator, not as '!'
  CHECK(parse_to_string("(x!=y)") == "(x!=y)");
}

TEST_CASE("parser: arithmetic binds tighter than comparison", "[parser]") {
  check_parse("a + b < c", "((a+b)<c)");
}

TEST_CASE("parser: subtraction is left-associative", "[parser]") {
  check_parse("1 - 2 - 3", "((1-2)-3)");
}

TEST_CASE("parser: parentheses override precedence", "[parser]") {
  check_parse("(1 + 2) * 3", "((1+2)*3)");
  check_parse("(x)", "x");
  check_parse("((1+2))", "(1+2)");
  check_parse("(a+b)*(c-d)", "((a+b)*(c-d))");
  check_parse("(x and y) or z", "((x and y) or z)");
}

TEST_CASE("parser: unary operators bind to a single operand", "[parser]") {
  check_parse("-x", "(-x)");
  check_parse("-42", "(-42)");
  check_parse("not x", "(!x)");
  check_parse("NOT x", "(!x)");
}

TEST_CASE("parser: unary operators stack", "[parser]") {
  check_parse("--x", "(-(-x))");
  check_parse("not not x", "(!(!x))");
}

TEST_CASE("parser: unary binds tighter than infix", "[parser]") {
  check_parse("-x + y", "((-x)+y)");
  check_parse("x + -y", "(x+(-y))");
  check_parse("-x * y", "((-x)*y)");
  check_parse("-x > y", "((-x)>y)");
  check_parse("not x and y", "((!x) and y)");
}

TEST_CASE("parser: unary applies to grouped and called operands", "[parser]") {
  check_parse("-(1+2)", "(-(1+2))");
  check_parse("not (x > 0)", "(!(x>0))");
  // Call outranks Prefix, so the call is the operand of the minus.
  check_parse("-f(x)", "(-f(x))");
  check_parse("f(-x)", "f((-x))");
  check_parse("[-1, -2]", "[(-1), (-2)]");
}

TEST_CASE("parser: call expressions", "[parser]") {
  check_parse("f()", "f()");
  check_parse("f('a', 2, x)", "f('a', 2, x)");
  // The callee is the attribute (a.b), not b applied to a.
  check_parse("a.b(1).c", "a.b(1).c");
}

TEST_CASE("parser: member access and method calls", "[parser]") {
  check_parse("x.foo", "x.foo");
  check_parse("x.sum()", "x.sum()");
  check_parse("x.where(y>0)", "x.where((y>0))");
  // Grouping looser than '.' is still shown.
  check_parse("(a+b).c", "(a+b).c");
  check_parse("f(x).y", "f(x).y");
  check_parse("f(x.y)", "f(x.y)");
}

TEST_CASE("parser: float literals keep their fractional part", "[parser]") {
  check_parse("1.5", "1.5");
  check_parse("1.5 + 2.5", "(1.5+2.5)");
  check_parse("0.1", "0.1");
  check_parse("2.5e-3", "0.0025");
  check_parse("3.14159265", "3.14159265"); // no longer rounded to float
  check_parse("1", "1");
}

TEST_CASE("parser: floats print in a form that lexes back as a float",
          "[parser]") {
  check_parse("1.5e3", "1500.0");
  check_parse("100.0", "100.0");
  check_parse("0.0", "0.0");
  check_parse("1e30", "1e+30");
  check_parse("1.0e-9", "1e-09");
  CHECK(parse_to_string("1500.0") == "1500.0");
  CHECK(parse_to_string("1e+30") == "1e+30");
  CHECK(parse_to_string("1e-09") == "1e-09");
}

// Printing must produce something that parses back to the same value.
TEST_CASE("parser: printing a float is the inverse of parsing it", "[parser]") {
  for (const auto* input : {
           "1.5", "0.1", "0.0025", "273.15", "3.14159265", "1500.0", "0.0",
           "1e+30", "1e-09", "1e+40", "1e-40",
           // Values that only survive at full double precision
           "1.0000000000000002", "2.2250738585072014e-308",
           "1.7976931348623157e+308",
       }) {
    INFO("Input: " << input);
    const auto once = parse_to_string(input);
    const auto twice = parse_to_string(once);
    CHECK(once == twice);
  }
}

TEST_CASE("parser: array literals", "[parser]") {
  check_parse("[]", "[]");
  check_parse("[1]", "[1]");
  check_parse("[1, 2, 3]", "[1, 2, 3]");
  check_parse("[[1,2],[3]]", "[[1, 2], [3]]");
}

TEST_CASE("parser: binary operators are left-associative", "[parser]") {
  check_parse("a.b.c", "a.b.c");
  check_parse("a and b and c", "((a and b) and c)");
}

TEST_CASE("parser: exponentiation is right-associative", "[parser]") {
  check_parse("2 ** 3 ** 2", "(2**(3**2))");
  check_parse("2 ** -x", "(2**(-x))");
}

TEST_CASE("parser: exponentiation outranks unary and arithmetic", "[parser]") {
  check_parse("-x ** 2", "(-(x**2))");
  check_parse("-x ** -y", "(-(x**(-y)))");
  check_parse("x ** 2 * 3", "((x**2)*3)");
  check_parse("x * y ** 2", "(x*(y**2))");
  check_parse("a ** b + c", "((a**b)+c)");
}

TEST_CASE("parser: attribute access outranks unary", "[parser]") {
  check_parse("-x.y", "(-x.y)");
  check_parse("-a.b.c", "(-a.b.c)");
  check_parse("not x.y", "(!x.y)");
}

TEST_CASE("parser: unbalanced delimiters are reported", "[parser]") {
  CHECK_THROWS_AS(parse_to_string("(1 + 2"), parser::ParserError);
  CHECK_THROWS_AS(parse_to_string("f(x"), parser::ParserError);
  CHECK_THROWS_AS(parse_to_string("[1,2"), parser::ParserError);
  // The expected/actual pair is readable rather than run together.
  check_rejected("(1 + 2", "Expected RightParen, got");
}

TEST_CASE("parser: unparseable input throws", "[parser]") {
  CHECK_THROWS_AS(parse_to_string("@"), parser::ParserError);
  check_rejected("@", "Illegal token");
  check_rejected(".x", "Unexpected Prefix Token");
  check_rejected("", "Unexpected Prefix Token");
}

TEST_CASE("parser: input must be consumed in full", "[parser]") {
  // Without this check a typo parses as its leading fragment and the rest is
  // silently dropped -- "x y" would quietly evaluate as "x".
  check_rejected("x y", "Unexpected trailing input");
  check_rejected("1 2", "Unexpected trailing input");
  check_rejected("x, y", "Unexpected trailing input");
  check_rejected("f(x) g(y)", "Unexpected trailing input");
  // Colon has a precedence but no infix handler, so it ends the expression.
  check_rejected("x:y", "Unexpected trailing input");
}

TEST_CASE("parser: illegal tokens are reported, not printed", "[parser]") {
  // These used to write to std::cout and let the parse succeed on the prefix.
  check_rejected("x % y", "Illegal token");
  check_rejected("x; y", "Illegal token");
}

TEST_CASE("parser: out-of-range literals are parser errors", "[parser]") {
  // Previously escaped as a bare std::out_of_range from stoi/stof.
  CHECK_THROWS_AS(parse_to_string("2147483648"), parser::ParserError);
  check_rejected("2147483648", "Integer literal out of range");
  check_rejected("1e400", "Float literal out of range");
  // The largest int still parses.
  check_parse("2147483647", "2147483647");
}

// Literals are held as double, so the range and precision of a threshold
// survive parsing. As float, 1e40 overflowed and 273.15 came back as
// 273.14999389648438.
TEST_CASE("parser: literals keep double range and precision", "[parser]") {
  check_parse("1e40", "1e+40");
  check_parse("1e-40", "1e-40");
  check_parse("273.15", "273.15");
  check_parse("0.1", "0.1");
  // 17 significant digits round-trip; a float would have flattened these two
  // onto the same value.
  check_parse("1.0000000000000002", "1.0000000000000002");
  check_parse("2.2250738585072014e-308", "2.2250738585072014e-308");
}

TEST_CASE("parser: errors say where they happened", "[parser]") {
  //                       1234567
  check_rejected("x + + y", "line 1, column 5");
  check_rejected("(1 + 2", "line 1, column 7");
  check_rejected("x y", "line 1, column 3");
  check_rejected("a @ b", "line 1, column 3");
  // Position is reported for the offending token, not for the start of input
  check_rejected("foo(a, b))", "line 1, column 10");
}

TEST_CASE("some parsed expression", "[parser]") {
  check_parse("x*y.derivative(dx=dy,['col']).where(x>0)",
              "(x*y.derivative((dx=dy), ['col']).where((x>0)))");
}

} // namespace dexpr
